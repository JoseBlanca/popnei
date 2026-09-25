//! The logistic mixed model: a binomial trait with a kinship.
//!
//! [`TheLinearization`] is the piece the fit is made of. A 0/1 trait is
//! turned into a continuous **working trait**, each individual carrying a
//! weight that says how much its 0 or 1 tells us at the fit so far, and a
//! weighted linear mixed model is fitted to that working trait; the
//! working trait and the weights are then made again from the new fit, and
//! so on. One pass of that is a linearization, and it is run with the
//! variance of the random effect of the kinship held where it was. "The
//! logistic mixed model" of `docs/specs/gwas.md` says what the whole fit
//! does, and "How popnei fits it, and why not pyNei's way" why the
//! covariance of the working trait is factored and solved against here
//! rather than inverted.

use popnei_linalg::{TheFirstOperand, TheSecondOperand};

use crate::error::{Error, Result};

use super::distributions::chi2_sf_1df;
use super::dosages::GwasDosages;
use super::grammar_gamma::GrammarGamma;
use super::logistic::{
    LogisticModel, TheSystemOfTheFit, the_chance_of, the_system_of, the_system_that_is_left,
};
use super::result::{Answers, NullModel};
use super::study::{Design, GwasInputShape, GwasModel, TestType};
use super::the_share_that_is_nothing;

/// How many rounds a linearization runs before a linear predictor that is
/// still moving is refused: 200, which is `GLMM_MAX_ITER` of
/// `pynei/gwas.py` and which nobody has measured.
///
/// The panel of `docs/specs/gwas.md`, 200 individuals and two covariates,
/// settles in 5 rounds at the variance of the kinship effect that GMMAT
/// fitted for it and in 26 at a variance of 1e10, which no fit reaches,
/// measured on both backends on 24 September 2026.
pub(super) const ROUNDS_OF_A_LINEARIZATION: usize = 200;

/// How small the largest change in the linear predictor has to be, over
/// its own largest absolute value plus 1, for a linearization to have
/// settled: 1e-6, which is `GLMM_TOL` of `pynei/gwas.py` and which nobody
/// has measured either.
pub(super) const THE_CHANGE_THAT_HAS_SETTLED: f64 = 1e-6;

/// The fit of the logistic mixed model as it stands between two
/// linearizations, with the buffers one linearization is run in.
///
/// The fit starts from the plain logistic null of
/// [`LogisticModel`], fitted with no kinship in it: its coefficients give
/// the first linear predictor and the first fitted chances, and the first
/// working trait and weights come from those. A fit started anywhere else
/// walks a different path through the bracket of the variance of the
/// kinship effect and can stop at another value of it, which "The logistic
/// mixed model" of `docs/specs/gwas.md` says and which no number of this
/// module would show.
///
/// [`TheLinearization::at`] runs one linearization at one value of that
/// variance, `tau` in the formulas. With `w` the weights, `k` the kinship
/// and `d` the design, the covariance of the working trait is
///
/// ```text
/// sigma = tau * k + w⁻¹
/// ```
///
/// with the reciprocals of the weights on its diagonal, and the
/// projection matrix of the working trait is
/// `p = sigma⁻¹ - sigma⁻¹ d (d' sigma⁻¹ d)⁻¹ d' sigma⁻¹`, which takes the
/// covariates out of anything it is applied to and weights it by that
/// covariance. `sigma` is never inverted here: it is factored with a
/// Cholesky, which costs a third of an inverse, and everything the round
/// does with `sigma⁻¹` is a solve against that factorization. What is left
/// when a linearization has settled is what the step on the variance needs
/// and what the score test of every variant is then made against.
pub(crate) struct TheLinearization<'a> {
    /// The trait of each tested individual, 0.0 or 1.0.
    phenotype: &'a [f64],
    /// The design the model is fitted on, `num_individuals` x `num_coefs`,
    /// row after row.
    design: &'a Design<'a>,
    /// The kinship of the tested individuals, `num_individuals` x
    /// `num_individuals`, row after row, in their order.
    kinship: &'a [f64],
    /// How many individuals the study tests.
    num_individuals: usize,
    /// How many columns the design has.
    num_coefs: usize,
    /// The effect of the intercept and of each covariate, one per column
    /// of the design, as a log odds ratio.
    coefs: Vec<f64>,
    /// The linear predictor of each tested individual, the design times
    /// those effects plus the random effect of the kinship, which the
    /// fitted chance of the individual is the logistic curve of.
    linear_predictor: Vec<f64>,
    /// That chance, one per tested individual.
    fitted_chance: Vec<f64>,
    /// `mu (1 - mu)` of each tested individual, which is how much a
    /// binomial trait of that fitted chance varies, and which the diagonal
    /// of the covariance holds the reciprocals of.
    weights: Vec<f64>,
    /// The continuous trait the round fits, one value per tested
    /// individual: the linear predictor plus the trait less the fitted
    /// chance over the weight.
    working_trait: Vec<f64>,
    /// `sigma`, `num_individuals` x `num_individuals`, row after row,
    /// whose lower half becomes its Cholesky factorization.
    covariance: Vec<f64>,
    /// The design laid out one row per column, `num_coefs` x
    /// `num_individuals`, which is the layout the solve against the
    /// covariance takes its right hand sides in: one for each column of
    /// the design.
    design_by_column: Vec<f64>,
    /// `sigma⁻¹ d` the same way round, `num_coefs` x `num_individuals`,
    /// one row for each column of the design.
    of_the_covariance: Vec<f64>,
    /// `d' sigma⁻¹ d`, `num_coefs` x `num_coefs`, whose lower half becomes
    /// its Cholesky factorization.
    dsd: Vec<f64>,
    /// What the design explains of the working trait through the
    /// covariance, `sigma⁻¹ d coefs`, one value per tested individual.
    explained: Vec<f64>,
    /// The working trait through the projection matrix, `p w`, one value
    /// per tested individual, which the step on the variance of the
    /// kinship effect is taken from.
    projected_working: Vec<f64>,
    /// The kinship times it, `k p w`, one value per tested individual.
    of_the_kinship: Vec<f64>,
    /// The linear predictor the round came to, which the change that ends
    /// a linearization is measured against the one before it.
    next_predictor: Vec<f64>,
    /// How many rounds the last linearization ran.
    rounds: usize,
    /// How many rounds every linearization of this fit has run in all,
    /// which is how many times the covariance of the working trait has
    /// been factored.
    ///
    /// It is the count "How popnei fits it, and why not pyNei's way" of
    /// `docs/specs/gwas.md` calls the linearizations of a fit, 22 on the
    /// panel, and it is kept because that is the claim deliverable 3 of
    /// `docs/plans/gwas-logistic.md` checks: the fit factors an
    /// individuals by individuals matrix once per round and inverts one
    /// only at the end. A count taken anywhere but at the call itself
    /// could say so while the code did otherwise.
    factorizations: usize,
}

impl<'a> TheLinearization<'a> {
    /// The fit of `phenotype` over `design` with `kinship` as the
    /// covariance of its random effect, at the plain logistic null it
    /// starts from.
    ///
    /// `phenotype` holds 0.0 or 1.0 for each individual the design has a
    /// row for, and `kinship` is that many rows of that many values, row
    /// after row, already cut to those individuals and in their order.
    /// `null` is the logistic model of the same trait and the same design
    /// fitted with no kinship in it, whose effects and linear predictor
    /// the fit starts at. A null of another design has an effect for a
    /// covariate this one has not, and the first product of a round is
    /// where the buffer that holds them is refused.
    ///
    /// # Errors
    ///
    /// [`Error::GwasInputOfAnotherSize`] when `phenotype` does not hold
    /// one value for each tested individual, when `null` does not hold one
    /// linear predictor for each of them, or when `kinship` does not hold
    /// one row and one column for each of them.
    pub(crate) fn of_the_logistic_null(
        phenotype: &'a [f64],
        design: &'a Design<'a>,
        kinship: &'a [f64],
        null: &LogisticModel,
    ) -> Result<TheLinearization<'a>> {
        let num_individuals = design.num_individuals();
        let num_coefs = design.num_coefs();
        if phenotype.len() != num_individuals || null.linear_predictor().len() != num_individuals {
            return Err(Error::GwasInputOfAnotherSize {
                problem: GwasInputShape::Phenotype {
                    num_values: phenotype.len(),
                    num_individuals,
                },
            });
        }
        // The kinship holds its values in one slice, so a count that does
        // not fit in a `usize` is a count no slice in memory has and the
        // matrix is of another size whichever way the two are compared.
        let of_the_individuals = num_individuals
            .checked_mul(num_individuals)
            .filter(|values| *values == kinship.len());
        let Some(of_the_individuals) = of_the_individuals else {
            return Err(Error::GwasInputOfAnotherSize {
                problem: GwasInputShape::Kinship {
                    num_values: kinship.len(),
                    num_individuals,
                },
            });
        };
        // `d' sigma⁻¹ d` is the columns of the design squared, and the
        // design itself holds more values than that: `Design::of_the_study`
        // refuses a study of no more individuals than its columns plus one,
        // so the rows are the columns plus two at least.
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "the design holds `num_individuals` times `num_coefs` values and \
                      `num_individuals` is `num_coefs` plus 2 at least, so this product is \
                      smaller than the length of the design and fits in a `usize`"
        )]
        let of_the_coefs = num_coefs * num_coefs;
        let mut fitted = TheLinearization {
            phenotype,
            design,
            kinship,
            num_individuals,
            num_coefs,
            coefs: null.coefs().to_vec(),
            linear_predictor: null.linear_predictor().to_vec(),
            fitted_chance: vec![0.0_f64; num_individuals],
            weights: vec![0.0_f64; num_individuals],
            working_trait: vec![0.0_f64; num_individuals],
            covariance: vec![0.0_f64; of_the_individuals],
            design_by_column: vec![0.0_f64; design.values().len()],
            of_the_covariance: vec![0.0_f64; design.values().len()],
            dsd: vec![0.0_f64; of_the_coefs],
            explained: vec![0.0_f64; num_individuals],
            projected_working: vec![0.0_f64; num_individuals],
            of_the_kinship: vec![0.0_f64; num_individuals],
            next_predictor: vec![0.0_f64; num_individuals],
            rounds: 0,
            factorizations: 0,
        };
        for (chance, predicted) in fitted
            .fitted_chance
            .iter_mut()
            .zip(&fitted.linear_predictor)
        {
            *chance = the_chance_of(*predicted);
        }
        // The design is laid out one row per column once, since it is the
        // same matrix at every linearization and the solve against the
        // covariance reads its right hand sides that way round.
        for (column, of_the_column) in fitted
            .design_by_column
            .chunks_exact_mut(num_individuals.max(1))
            .enumerate()
        {
            for (value, row) in of_the_column
                .iter_mut()
                .zip(design.values().chunks_exact(num_coefs.max(1)))
            {
                *value = row.get(column).copied().unwrap_or(f64::NAN);
            }
        }
        Ok(fitted)
    }

    /// One linearization at `genetic_variance`, the variance of the random
    /// effect of the kinship, which is `tau` in the formulas: it runs
    /// rounds of the weighted linear mixed model of the working trait
    /// until the linear predictor settles, and leaves the effects of the
    /// intercept and the covariates, that predictor and the fitted
    /// chances, the Cholesky factorization of the covariance and of
    /// `d' sigma⁻¹ d`, the covariance solved against the design, the
    /// working trait through the projection matrix and the kinship times
    /// it.
    ///
    /// A round is one weighted linear mixed model: the weights and the
    /// working trait are made from the fitted chances the last round left,
    /// the covariance is factored, the effects are solved for, and the new
    /// linear predictor is the design times those effects plus
    /// `genetic_variance` times the kinship times the working trait
    /// through the projection matrix. The linearization stops when the
    /// largest change in that predictor, over its own largest absolute
    /// value plus 1, falls below [`THE_CHANGE_THAT_HAS_SETTLED`].
    ///
    /// A `genetic_variance` of 0 is the boundary where the kinship
    /// explains nothing, and it is an ordinary value here: the covariance
    /// is then the weights alone and the round is one round of the
    /// iteratively reweighted least squares that fitted the logistic null.
    ///
    /// # Errors
    ///
    /// [`Error::GwasFitDidNotSettle`] with the round it reached, which
    /// three things give. The linear predictor is still moving after
    /// [`ROUNDS_OF_A_LINEARIZATION`] rounds. A weight has fallen to 0,
    /// which is the fitted chance of an individual reaching 0 or 1, and
    /// the working trait is then not a finite number; on `panel_called` a
    /// variance of 1e12, which no fit reaches, does it at the round 26. Or
    /// the pivots of `d' sigma⁻¹ d` have collapsed, by the rule of
    /// **Open 5** of `docs/specs/gwas.md`, or the factorization of that
    /// matrix was refused outright, which the weights do as the chances
    /// they come from reach 0 and 1 and which two covariates that carry
    /// nearly the same thing do on their own.
    ///
    /// [`Error::GwasKinshipNotACovariance`] when the covariance cannot be
    /// factored, with the row it stopped at: a weight is at most 0.25, so
    /// the reciprocals put 4 at least on every diagonal entry, and what
    /// takes such a matrix below 0 is a kinship whose own smallest
    /// eigenvalue is below 0 times a variance large enough to reach it.
    /// [`Error::GwasLinalg`] when one of the six products or one of the
    /// three solves could not be done.
    pub(crate) fn at(&mut self, genetic_variance: f64) -> Result<()> {
        self.rounds = 0;
        let mut settled = false;
        while self.rounds < ROUNDS_OF_A_LINEARIZATION && !settled {
            self.rounds = self.rounds.saturating_add(1);
            settled = self.the_round_at(genetic_variance)?;
        }
        match settled {
            true => Ok(()),
            false => Err(Error::GwasFitDidNotSettle {
                model: GwasModel::Glmm,
                rounds: self.rounds,
            }),
        }
    }

    /// The weights and the working trait of the fitted chances the fit
    /// stands at, which every round begins with and which the search over
    /// the variance reads once before the first one, for the variance it
    /// starts at.
    ///
    /// # Errors
    ///
    /// [`Error::GwasFitDidNotSettle`] with the round it reached when a
    /// weight has fallen to 0.
    fn the_working_trait(&mut self) -> Result<()> {
        // A weight of 0, which a fitted chance that has reached 0 or 1
        // gives, leaves the working trait of that individual a value that
        // is not a number and the diagonal of the covariance an infinity,
        // which the linear algebra crate refuses as a matrix that is not
        // finite, naming a matrix the user never saw. It is the fit
        // running away, the same thing the logistic null meets when the
        // design weighted by the weights can no longer be factored, so it
        // is refused here with the round it reached. The reciprocal is
        // read as well as the working trait: a weight small enough for it
        // to overflow leaves a working trait that is finite, since what is
        // divided by the weight there is as small as the weight.
        let mut collapsed = false;
        for (((weight, working), chance), (predicted, measured)) in self
            .weights
            .iter_mut()
            .zip(&mut self.working_trait)
            .zip(&self.fitted_chance)
            .zip(self.linear_predictor.iter().zip(self.phenotype))
        {
            *weight = chance * (1.0 - chance);
            *working = predicted + (measured - chance) / *weight;
            collapsed = collapsed || !working.is_finite() || !(1.0 / *weight).is_finite();
        }
        match collapsed {
            // The search reads the working trait once before any round is
            // run, for the variance it starts at, and the round counter is
            // still the 0 the logistic null left there. So the round named
            // is the one the fit was about to run and never a round it did
            // not run: a refusal that said it did not settle in the 0
            // rounds it was fitted in names a fit that was never made.
            true => Err(Error::GwasFitDidNotSettle {
                model: GwasModel::Glmm,
                rounds: self.rounds.max(1),
            }),
            false => Ok(()),
        }
    }

    /// One round of the linearization at `genetic_variance`, which gives
    /// whether the linear predictor has settled.
    ///
    /// # Errors
    ///
    /// What [`TheLinearization::at`] lists, but for the rounds running
    /// out, which it counts itself.
    fn the_round_at(&mut self, genetic_variance: f64) -> Result<bool> {
        let num_individuals = self.num_individuals;
        let num_coefs = self.num_coefs;
        self.the_working_trait()?;
        for ((row, of_the_kinship), individual) in self
            .covariance
            .chunks_exact_mut(num_individuals.max(1))
            .zip(self.kinship.chunks_exact(num_individuals.max(1)))
            .zip(0..)
        {
            for (value, of_the_pair) in row.iter_mut().zip(of_the_kinship) {
                *value = genetic_variance * of_the_pair;
            }
            // The diagonal carries the reciprocal of that individual's
            // weight, which is what makes the covariance the one of the
            // working trait and not of the random effect alone.
            if let (Some(value), Some(weight)) =
                (row.get_mut(individual), self.weights.get(individual))
            {
                *value += 1.0 / weight;
            }
        }
        self.factorizations = self.factorizations.saturating_add(1);
        the_covariance_that_was_factored(popnei_linalg::cholesky_lower(
            &mut self.covariance,
            num_individuals,
        ))?;
        self.of_the_covariance
            .copy_from_slice(&self.design_by_column);
        popnei_linalg::solve_with_cholesky(
            &self.covariance,
            num_individuals,
            &mut self.of_the_covariance,
            num_coefs,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "solve of the design against the covariance of the working trait",
            source,
        })?;
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &self.design_by_column,
                rows: num_coefs,
            },
            num_individuals,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &self.of_the_covariance,
                cols: num_coefs,
            },
            &mut self.dsd,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design with the design solved against the covariance",
            source,
        })?;
        // The effects of the intercept and the covariates are the solution
        // of `(d' sigma⁻¹ d) coefs = d' sigma⁻¹ w`, so the right hand side
        // goes into the buffer of the effects and the solve overwrites it.
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &self.of_the_covariance,
                rows: num_coefs,
            },
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.working_trait,
                cols: 1,
            },
            &mut self.coefs,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design solved against the covariance with the working \
                        trait",
            source,
        })?;
        // The factorization of `d' sigma⁻¹ d` and the pivots it wrote are
        // one answer: a matrix it refuses and one whose smallest pivot has
        // fallen to nothing against its largest are both the weighted
        // design of the fit having collapsed, which **Open 5** of
        // `docs/specs/gwas.md` measures on the logistic model.
        let of_the_system = match the_system_of(
            popnei_linalg::cholesky_lower(&mut self.dsd, num_coefs),
            "factorization of the design weighted by the covariance of the working trait",
        )? {
            TheSystemOfTheFit::Worked => {
                the_system_that_is_left(&self.dsd, num_coefs, num_individuals)
            }
            TheSystemOfTheFit::RanAway => TheSystemOfTheFit::RanAway,
        };
        match of_the_system {
            TheSystemOfTheFit::Worked => {}
            TheSystemOfTheFit::RanAway => {
                return Err(Error::GwasFitDidNotSettle {
                    model: GwasModel::Glmm,
                    rounds: self.rounds,
                });
            }
        }
        // The factorization above succeeded, so the solve against it
        // cannot find the system singular today. What it answers is read
        // all the same, as both callers in `logistic` read it: a solve
        // whose answer is thrown away leaves a set of effects that nothing
        // solved for and that the round would carry on with.
        match the_system_of(
            popnei_linalg::solve_with_cholesky(&self.dsd, num_coefs, &mut self.coefs, 1),
            "solve of the design weighted by the covariance of the working trait",
        )? {
            TheSystemOfTheFit::Worked => {}
            TheSystemOfTheFit::RanAway => {
                return Err(Error::GwasFitDidNotSettle {
                    model: GwasModel::Glmm,
                    rounds: self.rounds,
                });
            }
        }
        // The working trait through the projection matrix is the working
        // trait solved against the covariance, less what the design
        // explains of it there.
        self.projected_working.copy_from_slice(&self.working_trait);
        popnei_linalg::solve_with_cholesky(
            &self.covariance,
            num_individuals,
            &mut self.projected_working,
            1,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "solve of the working trait against its covariance",
            source,
        })?;
        popnei_linalg::product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &self.of_the_covariance,
                rows: num_individuals,
            },
            num_coefs,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.coefs,
                cols: 1,
            },
            &mut self.explained,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design solved against the covariance with the effects of \
                        the logistic mixed model",
            source,
        })?;
        for (projected, explained) in self.projected_working.iter_mut().zip(&self.explained) {
            *projected -= *explained;
        }
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: self.kinship,
                rows: num_individuals,
            },
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.projected_working,
                cols: 1,
            },
            &mut self.of_the_kinship,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the kinship with the working trait through the projection \
                        matrix",
            source,
        })?;
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: self.design.values(),
                rows: num_individuals,
            },
            num_coefs,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.coefs,
                cols: 1,
            },
            &mut self.next_predictor,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design with the effects of the logistic mixed model",
            source,
        })?;
        for (predicted, of_the_kinship) in self.next_predictor.iter_mut().zip(&self.of_the_kinship)
        {
            *predicted += genetic_variance * of_the_kinship;
        }
        // The change is divided by the largest absolute value of the
        // predictor the round started at and not of the one it came to,
        // which is `_fit_pql_for_tau` of `pynei/gwas.py`, and the plus one
        // keeps a predictor near 0 from asking for more than a change of
        // 1e-6 in absolute terms. It is written as a division and compared
        // with the tolerance, as pyNei writes it, and not as the
        // multiplication that would save it: the two fits are asked to
        // take the same rounds in the same order.
        let of_the_predictor = self
            .linear_predictor
            .iter()
            .fold(0.0_f64, |largest, predicted| largest.max(predicted.abs()));
        let change = self
            .next_predictor
            .iter()
            .zip(&self.linear_predictor)
            .fold(0.0_f64, |largest, (next, predicted)| {
                largest.max((next - predicted).abs())
            });
        self.linear_predictor.copy_from_slice(&self.next_predictor);
        for (chance, predicted) in self.fitted_chance.iter_mut().zip(&self.linear_predictor) {
            *chance = the_chance_of(*predicted);
        }
        Ok(change / (of_the_predictor + 1.0) < THE_CHANGE_THAT_HAS_SETTLED)
    }
}

/// What the factorization of the covariance of the working trait came to:
/// a matrix that is not positive definite is a kinship that is not a
/// covariance, and every other refusal is an error.
///
/// # Errors
///
/// [`Error::GwasKinshipNotACovariance`] when the factorization stopped at
/// a row, with that row, and [`Error::GwasLinalg`] when the call was
/// refused for any other reason.
fn the_covariance_that_was_factored(what_it_gave: popnei_linalg::Result<()>) -> Result<()> {
    match what_it_gave {
        Ok(()) => Ok(()),
        Err(popnei_linalg::Error::Singular { at, .. }) => {
            Err(Error::GwasKinshipNotACovariance { at })
        }
        // The error type of the linear algebra crate is `non_exhaustive`,
        // and a case nobody has written yet lands here, which is where
        // every refusal but a matrix that could not be factored belongs:
        // the call was wrong, a defect of popnei and a `RuntimeError` in
        // Python.
        Err(source) => Err(Error::GwasLinalg {
            operation: "factorization of the covariance of the working trait of the logistic \
                        mixed model",
            source,
        }),
    }
}

/// A refusal of the plain logistic null this fit starts from, under the
/// model the user asked for.
///
/// That null is fitted by [`LogisticModel`], which has no kinship and says
/// so in the one refusal of it that names a model: a fit that did not
/// settle carries which model was being fitted, and its message gives that
/// model's name and that model's remedies. A user who brought a kinship
/// asked for the logistic mixed model, so a null of theirs that walks
/// towards an infinite coefficient is named as that model, with the
/// kinship among the three things to look at. Every other refusal of the
/// null says the same thing for both models and is carried through as it
/// is.
fn the_refusal_of_this_model(of_the_null: Error) -> Error {
    if let Error::GwasFitDidNotSettle {
        model: GwasModel::Glm,
        rounds,
    } = of_the_null
    {
        return Error::GwasFitDidNotSettle {
            model: GwasModel::Glmm,
            rounds,
        };
    }
    of_the_null
}

/// The two values of the variance of the kinship effect the search
/// remembers: one whose derivative asks for a larger variance and one
/// whose derivative asks for a smaller one.
///
/// The answer lies between them, so a Newton step that would leave the
/// interval is replaced by the geometric mean of its two ends. It is what
/// keeps the search from cycling, and `docs/reports/glmm-method/README.md`
/// records that the first fit written for that report, without it, did not
/// converge at all: the two updates, the linearization's and the step on
/// the variance, fight each other, which pyNei's `_GLMMNull` says of its
/// own history too.
///
/// Both ends are `None` until a step has been taken, and the search comes
/// down from a variance that is too large rather than too small, so the
/// end that asks for a smaller variance is usually the first to be known.
struct TheBracket {
    /// The largest variance whose derivative asked for a larger one, and
    /// `None` while no step has asked for one.
    too_small: Option<f64>,
    /// The smallest variance whose derivative asked for a smaller one.
    too_large: Option<f64>,
}

impl TheBracket {
    /// The bracket of a search that has taken no step.
    fn of_a_search_that_has_not_started() -> TheBracket {
        TheBracket {
            too_small: None,
            too_large: None,
        }
    }

    /// `genetic_variance` remembered as the end of the bracket that
    /// `score`, the first derivative of the restricted likelihood there,
    /// says it is: a derivative above 0 asks for a larger variance and one
    /// below it for a smaller one.
    fn the_end_that(&mut self, genetic_variance: f64, score: f64) {
        match score > 0.0 {
            true => self.too_small = Some(genetic_variance),
            false => self.too_large = Some(genetic_variance),
        }
    }

    /// Where the search goes from `genetic_variance` after a Newton
    /// `step`.
    ///
    /// It is that step, unless the bracket is known at both ends and the
    /// step would leave it, and then it is the geometric mean of the two
    /// ends. Before both ends are known there is no interval to stay
    /// inside, and a step that would take the variance to 0 or below, or
    /// to a value that is not a finite number, quarters it instead, which
    /// is a smaller move in the same direction.
    ///
    /// A step that is not finite is what an average information that has
    /// underflowed to 0 gives, and it is the everyday case and not a
    /// remote one: a kinship of all zeros leaves the derivative and the
    /// information both exactly 0 only because every quotient of the
    /// trace's triangular solve rounds to exactly 1, and one unit in the
    /// last place of that solve's right hand side leaves a derivative of
    /// 5e-15 over an information of 0. Without the test for a finite
    /// number here that infinity is the variance the next linearization is
    /// run at, and the covariance is then an infinity times a kinship
    /// entry of 0, which is not a number.
    fn the_variance_after(&self, genetic_variance: f64, step: f64) -> f64 {
        let next = genetic_variance + step;
        match (self.too_small, self.too_large) {
            (Some(too_small), Some(too_large)) => match too_small < next && next < too_large {
                true => next,
                false => (too_small * too_large).sqrt(),
            },
            (None, _) | (_, None) => match next > 0.0 && next.is_finite() {
                true => next,
                false => genetic_variance / 4.0,
            },
        }
    }
}

/// How many steps on the variance of the kinship effect the search takes
/// before a variance that is still moving is refused: 200, which is the
/// same `GLMM_MAX_ITER` of `pynei/gwas.py` that bounds the rounds of a
/// linearization and which nobody has measured.
///
/// The panel of `docs/specs/gwas.md` takes 8 of them.
const STEPS_ON_THE_VARIANCE: usize = 200;

/// How small a step on the variance of the kinship effect has to be for
/// the fit to have settled, as a share of that variance: 1e-6, which is
/// the same `GLMM_TOL` of `pynei/gwas.py` that
/// [`THE_CHANGE_THAT_HAS_SETTLED`] is, and which nobody has measured
/// either.
///
/// A step settles the fit when its absolute value falls below this times
/// the variance plus this, so that a variance near 0 asks for a step below
/// 1e-12 in absolute terms rather than for one that is small against
/// nothing. The same number is the floor below which a variance is taken
/// to be 0, the boundary where the kinship explains nothing.
const THE_STEP_THAT_HAS_SETTLED: f64 = 1e-6;

/// The null model of the logistic mixed model: the variance of the random
/// effect of the kinship, the effects of the intercept and the covariates,
/// the projection matrix every variant is then tested through and the
/// residual of the trait it is tested against.
///
/// The fit is the search over that variance. For a variance held where it
/// is, [`TheLinearization`] runs rounds of a weighted linear mixed model
/// on the working trait until the linear predictor settles; the variance
/// then takes one Newton step from the restricted maximum likelihood, and
/// the whole thing starts again. With `p` the projection matrix, `k` the
/// kinship, `w` the working trait of the linearization that just finished
/// and `pw = p w`, the step is
///
/// ```text
/// score = 0.5 * (pw' k pw - trace(p k))
/// ai    = 0.5 * (k pw)' p (k pw)
/// step  = score / ai
/// ```
///
/// where `ai` is the average information, the average of the observed
/// second derivative and the one expected under the model: the terms that
/// cost the most appear in the two with opposite signs and cancel.
///
/// The fit starts from the plain logistic null of [`LogisticModel`],
/// fitted with no kinship in it, and nowhere else, and the variance starts
/// at half the variance of the first working trait, which is what it would
/// be if the kinship explained all of it, halved. So the search comes down
/// from above, and [`TheBracket`] is what keeps it from cycling.
///
/// Nothing of an individuals by individuals size is inverted while the
/// search runs: `trace(p k)` is the one quantity that seems to need every
/// entry of the covariance's inverse and comes from the identity
/// [`TheStepOnTheVariance::the_trace_of`] gives instead. The inverse is
/// formed once, when the search has settled, because the score test of a
/// variant wants the projection matrix as a matrix.
pub(crate) struct LogisticMixedModel {
    /// The effect of the intercept and of each covariate, one per column
    /// of the design, as a log odds ratio.
    coefs: Vec<f64>,
    /// The variance of the random effect of the kinship the search landed
    /// on.
    genetic_variance: f64,
    /// The projection matrix `p`, `num_individuals` x `num_individuals`,
    /// row after row, which every variant is tested through.
    projection: Vec<f64>,
    /// The residual every variant's numerator is taken against, one value
    /// per tested individual: the trait less the fitted chance.
    ///
    /// "What it gives" of the model in `docs/specs/gwas.md` calls it the
    /// residual `p y` of the score test, which is the name the linear
    /// mixed model's own residual carries, and here it is the trait less
    /// the fitted chance and not the trait through the projection matrix.
    /// What that matrix leaves of the *working* trait is what equals it,
    /// which is the fit's own optimality condition and is the cheapest
    /// evidence there is that the search reached its optimum.
    projected_trait: Vec<f64>,
    /// The largest value of the diagonal of that matrix, which is what a
    /// variant's own squared length is weighted by to say how much of the
    /// variant the projection has left.
    ///
    /// It is the scale of the third place of the meanwhile of **Open 2** of
    /// `docs/specs/gwas.md`, and it is the one the linear mixed model's
    /// score test already takes, so that the two mixed models answer alike
    /// rather than each picking its own. The doc comment of the field of
    /// the same name in `linear_mixed` measures how far under the largest
    /// eigenvalue of the projection, which is what really bounds
    /// `x' p x / x' x`, that diagonal sits: 1.7 times on both panels of
    /// that model, so the threshold is tighter than it was meant to be and
    /// not looser.
    largest_of_the_projection: f64,
    /// How many individuals the study tests.
    num_individuals: usize,
    /// How many rounds the fit ran in all, which is how many times it
    /// factored the covariance of the working trait: 22 on the panel.
    ///
    /// The three counters are in no result and the cargo test of
    /// deliverable 3 of `docs/plans/gwas-logistic.md` is the one thing that
    /// reads them, which is why they are taken at the calls themselves: a
    /// count taken anywhere else could say the fit factors and inverts as
    /// "How popnei fits it, and why not pyNei's way" of
    /// `docs/specs/gwas.md` describes while the code did otherwise, and a
    /// fit that inverted where it should solve would give the same numbers
    /// and nothing else would notice.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the three counters are in no result and the cargo test that counts \
                      the factorizations and the inverses is what reads them"
        )
    )]
    linearizations: usize,
    /// How many steps it took on the variance of the kinship effect: 8 on
    /// the panel.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the three counters are in no result and the cargo test that counts \
                      the factorizations and the inverses is what reads them"
        )
    )]
    steps_on_the_variance: usize,
    /// How many times it formed the inverse of a factorized matrix of that
    /// size: 1, at the end, for the projection matrix.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the three counters are in no result and the cargo test that counts \
                      the factorizations and the inverses is what reads them"
        )
    )]
    inverses: usize,
    /// The GRAMMAR-Gamma approximation, when the study asked for it and
    /// [`LogisticMixedModel::approximate_the_denominator`] estimated its
    /// factor from the first block of the second pass, and `None` when
    /// every variant gets the exact denominator.
    approximation: Option<GrammarGamma>,
    /// The dosages of a block through the projection matrix, `x p`, the
    /// variants that have variance x `num_individuals`. It stays empty
    /// under the approximation, which is the product it does not make.
    projected: Vec<f64>,
    /// The denominator of each variant that has variance, `x' p x`, one per
    /// variant, formed from the row above or approximated.
    den: Vec<f64>,
    /// Each variant times the residual the null left, `x' r`, one per
    /// variant that has variance, which is the `x' p y` of the score test
    /// with the residual this model has in place of `p y`.
    num: Vec<f64>,
    /// The effect of each variant that has variance.
    beta: Vec<f64>,
    /// The standard error of each of those effects.
    se: Vec<f64>,
    /// The p-value of each of those tests.
    p_value: Vec<f64>,
}

impl LogisticMixedModel {
    /// The logistic mixed model of `phenotype` over `design` with
    /// `kinship` as the covariance of its random effect, fitted without
    /// any variant in it.
    ///
    /// `phenotype` holds 0.0 or 1.0 for each individual the design has a
    /// row for, and `kinship` is that many rows of that many values, row
    /// after row, already cut to those individuals and in their order.
    /// The plain logistic null the fit starts from is fitted here, from
    /// the same trait and the same design with no kinship in it, so that a
    /// caller cannot start it anywhere else.
    ///
    /// # Errors
    ///
    /// [`Error::GwasInputOfAnotherSize`] when `phenotype` does not hold
    /// one value for each tested individual or `kinship` does not hold one
    /// row and one column for each of them.
    /// [`Error::GwasFitDidNotSettle`] with the steps it reached when the
    /// variance of the kinship effect is still moving after
    /// [`STEPS_ON_THE_VARIANCE`] of them, and whatever a linearization is
    /// refused for, with the round of that linearization.
    /// [`Error::GwasKinshipNotACovariance`] when the covariance of the
    /// working trait cannot be factored. [`Error::GwasLinalg`] when one of
    /// the products, one of the solves or the one inverse could not be
    /// done.
    pub(crate) fn of_the_study(
        phenotype: &[f64],
        design: &Design<'_>,
        kinship: &[f64],
    ) -> Result<LogisticMixedModel> {
        let num_individuals = design.num_individuals();
        let num_coefs = design.num_coefs();
        let null = LogisticModel::of_the_study(phenotype, design, TestType::Score)
            .map_err(the_refusal_of_this_model)?;
        let mut fitted = TheLinearization::of_the_logistic_null(phenotype, design, kinship, &null)?;
        drop(null);
        let mut step = TheStepOnTheVariance::of(num_individuals, num_coefs);
        // The variance the search starts at is half the variance of the
        // working trait the logistic null leaves, which is what the
        // variance of the kinship effect would be if the kinship explained
        // the whole of it. So the start is too large rather than too
        // small, and the bracket comes down from above.
        fitted.the_working_trait()?;
        let mut genetic_variance = the_variance_of(&fitted.working_trait) / 2.0;
        let mut bracket = TheBracket::of_a_search_that_has_not_started();
        let mut steps_on_the_variance = 0_usize;
        let mut settled = false;
        while steps_on_the_variance < STEPS_ON_THE_VARIANCE && !settled {
            steps_on_the_variance = steps_on_the_variance.saturating_add(1);
            // The variance the next linearization is run at is refused
            // here when it is not a finite number, with the step it
            // reached, and not left to the covariance it would build: an
            // infinity times a kinship entry of 0 is not a number, which
            // the linear algebra crate refuses as a matrix the user never
            // saw, a `RuntimeError` in Python and so a defect of popnei,
            // for a kinship the user brought. What
            // [`TheBracket::the_variance_after`] does with a step that is
            // not finite is what keeps this arm from firing on a kinship
            // whose average information has underflowed; the variance the
            // search starts at is the other way in, since the variance of
            // a working trait of finite values can still overflow.
            if !genetic_variance.is_finite() {
                return Err(Error::GwasFitDidNotSettle {
                    model: GwasModel::Glmm,
                    rounds: steps_on_the_variance,
                });
            }
            fitted.at(genetic_variance)?;
            let (score, information) = step.at(&fitted, genetic_variance)?;
            // A step that is not a finite number is not refused here, and
            // the bracket is what makes that safe: a kinship that says
            // nothing about the trait leaves the derivative and the
            // information both 0 and the step their quotient, and
            // [`TheBracket::the_variance_after`] sends every step that
            // does not land on a finite variance above 0, a NaN and both
            // infinities among them, to a quarter of the variance or to
            // the geometric mean of two ends that are finite. So the
            // variance walks down to the boundary where the kinship
            // explains nothing instead of becoming a number that is not
            // one. Measured on 24 September 2026 on the panel with every
            // genotype called and a kinship of all zeros: 12 steps to a
            // variance of 0, and the same 12 over a kinship of 1e-165
            // times the identity, whose average information underflows.
            let of_the_step = score / information;
            if of_the_step.abs()
                < THE_STEP_THAT_HAS_SETTLED * (genetic_variance + THE_STEP_THAT_HAS_SETTLED)
            {
                settled = true;
            } else {
                bracket.the_end_that(genetic_variance, score);
                let next = bracket.the_variance_after(genetic_variance, of_the_step);
                // A variance below the tolerance is the boundary where the
                // kinship explains nothing, which is 0 and not a small
                // number, and a second one there ends the fit: the search
                // has nowhere below to go.
                match next < THE_STEP_THAT_HAS_SETTLED {
                    true => {
                        settled = genetic_variance < THE_STEP_THAT_HAS_SETTLED;
                        genetic_variance = 0.0;
                    }
                    false => genetic_variance = next,
                }
            }
        }
        if !settled {
            return Err(Error::GwasFitDidNotSettle {
                model: GwasModel::Glmm,
                rounds: steps_on_the_variance,
            });
        }
        // The buffer the trace was taken in is done with, and it is as
        // large as the projection matrix, so it is what the part of the
        // covariance's inverse that the design explains is formed in.
        let of_the_design = step.take_the_buffer_of_the_trace();
        drop(step);
        let mut inverses = 0_usize;
        let projection = the_projection_of(&fitted, of_the_design, &mut inverses)?;
        let projected_trait = fitted
            .phenotype
            .iter()
            .zip(&fitted.fitted_chance)
            .map(|(measured, chance)| measured - chance)
            .collect();
        Ok(LogisticMixedModel {
            coefs: std::mem::take(&mut fitted.coefs),
            genetic_variance,
            largest_of_the_projection: projection
                .chunks_exact(num_individuals.max(1))
                .zip(0..)
                .filter_map(|(row, at)| row.get(at).copied())
                .fold(0.0_f64, f64::max),
            projection,
            projected_trait,
            num_individuals,
            linearizations: fitted.factorizations,
            steps_on_the_variance,
            inverses,
            approximation: None,
            projected: Vec::new(),
            den: Vec::new(),
            num: Vec::new(),
            beta: Vec::new(),
            se: Vec::new(),
            p_value: Vec::new(),
        })
    }

    /// The null model of the result: the effects of the intercept and of
    /// the covariates, as log odds ratios, and the variance of the random
    /// effect of the kinship.
    ///
    /// A logistic model has no free residual variance, since the variance
    /// of a binomial trait is decided by its mean, so `residual_variance`
    /// and the `heritability` that is built from the two are `None`, which
    /// is what "What it gives" of the model in `docs/specs/gwas.md` says.
    #[must_use]
    pub(crate) fn null_model(&self, test: TestType) -> NullModel {
        NullModel {
            model: GwasModel::Glmm,
            test,
            covariate_effects: self.coefs.clone(),
            residual_variance: None,
            genetic_variance: Some(self.genetic_variance),
            heritability: None,
            num_individuals: self.num_individuals,
        }
    }

    /// Estimates the factor of the GRAMMAR-Gamma approximation from
    /// `dosages`, the first block of the second pass, and makes every
    /// variant tested after it take the approximate denominator.
    ///
    /// "The GRAMMAR-Gamma approximation" of `docs/specs/gwas.md` says what
    /// it buys and what it costs: the denominator of a variant stops being
    /// a product with the projection matrix, which grows with the square of
    /// the individuals, and becomes one factor times the squared length of
    /// the variant's centered dosages, which grows with the individuals
    /// alone. The factor stands for a quantity that differs from variant to
    /// variant, so what it costs in accuracy grows with how strongly the
    /// panel is structured.
    ///
    /// # Errors
    ///
    /// What [`GrammarGamma::of_the_first_block`] refuses of that block: no
    /// variant of it that varies among the tested individuals, and a factor
    /// that is not a finite number above 0, which is what a block of
    /// nothing but variants the design explains gives, since the ratio of
    /// such a variant is left out of the mean.
    pub(crate) fn approximate_the_denominator(&mut self, dosages: &GwasDosages) -> Result<()> {
        self.approximation = Some(GrammarGamma::of_the_first_block(
            &self.projection,
            self.num_individuals,
            self.largest_of_the_projection,
            dosages,
        )?);
        Ok(())
    }

    /// The factor of the approximation, and `None` for a model that makes
    /// the exact denominator.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the factor is in no result and the test that compares it with the \
                      one pyNei's `estimate_gamma` gives is what reads it"
        )
    )]
    pub(crate) fn grammar_gamma_factor(&self) -> Option<f64> {
        self.approximation
            .map(|approximation| approximation.factor())
    }

    /// The denominator of every variant of a block that has variance, in
    /// the order of the block, left in `den`.
    ///
    /// Without the approximation it is `x' p x`, the variant through the
    /// projection matrix and then against itself, and the product of the
    /// block with that matrix is the larger half of the cost of a block.
    /// With it, it is the factor times the squared length of the variant's
    /// centered dosages, and no product is made.
    ///
    /// # Errors
    ///
    /// [`Error::GwasVariantsTooLarge`] when the values of the block are
    /// more than a `usize` counts, and [`Error::GwasLinalg`] when the
    /// product could not be done, which is where a block of other
    /// individuals than the null model was fitted over is refused.
    fn the_denominators_of(&mut self, dosages: &GwasDosages, num_vars: usize) -> Result<()> {
        self.den.clear();
        let of_the_variants = dosages.dosages().chunks_exact(self.num_individuals.max(1));
        if let Some(approximation) = self.approximation {
            self.den
                .extend(of_the_variants.map(|of_the_variant| approximation.den_of(of_the_variant)));
            return Ok(());
        }
        let values = num_vars
            .checked_mul(self.num_individuals)
            .ok_or(Error::GwasVariantsTooLarge)?;
        self.projected.resize(values, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: dosages.dosages(),
                rows: num_vars,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.projection,
                cols: self.num_individuals,
            },
            &mut self.projected,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of a block of variants with the projection matrix",
            source,
        })?;
        self.den.extend(
            self.projected
                .chunks_exact(self.num_individuals.max(1))
                .zip(dosages.dosages().chunks_exact(self.num_individuals.max(1)))
                .map(|(row, of_the_variant)| {
                    row.iter()
                        .zip(of_the_variant)
                        .map(|(projected, dosage)| projected * dosage)
                        .sum::<f64>()
                }),
        );
        Ok(())
    }

    /// The score test of every variant of a block that has variance among
    /// the tested individuals, in the order of the block.
    ///
    /// The score test is the only test this model has: a Wald test would
    /// fit one mixed model per variant, and
    /// [`the_model_and_the_test`](super::study::the_model_and_the_test)
    /// refuses a study that asks for one. It asks how steeply the
    /// likelihood would rise if the variant's effect were let off 0,
    /// measured at the null, so it holds the variance of the kinship effect
    /// and the fitted chances where the null left them and needs no fit per
    /// variant. It is what GMMAT's `glmm.score` makes.
    ///
    /// Every variant goes through the projection matrix of the null, which
    /// takes the covariates out of it and weights it by the covariance of
    /// the working trait. With `x` the dosages of a variant, `num` is
    /// `x' r` and `den` is `x' p x`, and then `beta` is `num / den`, `se`
    /// is `1 / sqrt(den)` and `num² / den` is read against a chi square
    /// with one degree of freedom. `r`, the residual the numerator is taken
    /// against, is the trait less the fitted chance and not the trait
    /// through the projection matrix, which is what a logistic mixed model
    /// has in place of the linear one's `p y` and which
    /// [`LogisticMixedModel::of_the_study`] left in `projected_trait`.
    ///
    /// The two products are the whole cost of a block, and the first of
    /// them, the dosages through the projection matrix, is what the
    /// GRAMMAR-Gamma approximation of `docs/specs/gwas.md` stands in for. A
    /// model that [`LogisticMixedModel::approximate_the_denominator`] was
    /// called on makes it no more.
    ///
    /// A variant of which the projection leaves at most the tested
    /// individuals times 2.2e-16 of what there was has no answer, and gets
    /// the three NaNs a variant with no variance gets. What there was is
    /// the variant's own squared length times the largest value of the
    /// diagonal of the projection matrix, which is the scale the linear
    /// mixed model's score test takes, so that the two mixed models answer
    /// a variant there is nothing left to test alike rather than each
    /// picking its own. It is **Open 2** of `docs/specs/gwas.md`.
    ///
    /// That comparison is made against whichever of the two denominators
    /// the study formed, and under the approximation it stops firing, which
    /// "Open 2's threshold under the approximation" of that spec states and
    /// which is a decision of 24 September 2026: the approximate
    /// denominator is a factor above 0 times a sum of squares, so it holds
    /// no cancellation and it is above 0 for every variant that varies,
    /// whatever the projection would have left of that variant. The linear
    /// mixed model's own `test_the_block` has what was measured, and the
    /// two models take the same decision here as they take the same scale.
    ///
    /// # Errors
    ///
    /// [`Error::GwasVariantsTooLarge`] when the values of the block are
    /// more than a `usize` counts, and [`Error::GwasLinalg`] when one of
    /// the two products could not be done, which is where a block of other
    /// individuals than the null model was fitted over is refused.
    pub(crate) fn test_the_block(&mut self, dosages: &GwasDosages) -> Result<Answers<'_>> {
        let num_vars = dosages.num_with_variance();
        self.beta.clear();
        self.se.clear();
        self.p_value.clear();
        if num_vars == 0 {
            // No variant to test, and the products below take a matrix of
            // one column at least. The block still has its rows in the
            // result, with the three NaNs of a variant that has no answer.
            return Ok(self.answers());
        }
        self.the_denominators_of(dosages, num_vars)?;
        self.num.resize(num_vars, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: dosages.dosages(),
                rows: num_vars,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.projected_trait,
                cols: 1,
            },
            &mut self.num,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of a block of variants with the residuals of the null model",
            source,
        })?;
        // The share of what the variant was that the projection has to
        // leave of it for the variant to be worth testing.
        let share_that_is_nothing = the_share_that_is_nothing(self.num_individuals);
        let largest_of_the_projection = self.largest_of_the_projection;
        // The squared length of each variant's dosages is the one the block
        // summed where it wrote the row, on the threads of rayon: reading
        // it here would be a second full read of the block's dosages on
        // this one thread for one value per variant.
        for ((den, num), of_the_dosages) in self
            .den
            .iter()
            .copied()
            .zip(&self.num)
            .zip(dosages.sum_of_squares().iter().copied())
        {
            if den <= share_that_is_nothing * largest_of_the_projection * of_the_dosages {
                self.beta.push(f64::NAN);
                self.se.push(f64::NAN);
                self.p_value.push(f64::NAN);
                continue;
            }
            self.beta.push(num / den);
            self.se.push(1.0 / den.sqrt());
            self.p_value.push(chi2_sf_1df(num * num / den));
        }
        Ok(self.answers())
    }

    /// The three columns of what the block last tested answered.
    #[must_use]
    fn answers(&self) -> Answers<'_> {
        Answers {
            beta: &self.beta,
            se: &self.se,
            p_value: &self.p_value,
        }
    }
}

/// The variance of `values`, the mean of the squares of their distances
/// from their own mean, which is what `numpy.var` gives and what the
/// variance the search starts at is half of.
fn the_variance_of(values: &[f64]) -> f64 {
    let count = values.len() as f64;
    let mean = values.iter().sum::<f64>() / count;
    values
        .iter()
        .map(|value| {
            let from_the_mean = value - mean;
            from_the_mean * from_the_mean
        })
        .sum::<f64>()
        / count
}

/// The buffers one step on the variance of the kinship effect is taken in,
/// made once for the fit and written over at every step.
///
/// What the step needs of the projection matrix is the trace of it times
/// the kinship and its product with one vector, and neither wants that
/// matrix: the linearization leaves the covariance of the working trait
/// factored, and everything here is a solve against that factorization or
/// a product with a handful of columns.
struct TheStepOnTheVariance {
    /// `num_individuals` x `num_individuals`: the right hand sides the
    /// trace's triangular solve is given, one for each individual, which
    /// come back as the solutions.
    of_the_trace: Vec<f64>,
    /// The kinship times the covariance solved against the design,
    /// `k sigma⁻¹ d`, `num_individuals` x `num_coefs`, row after row.
    of_the_kinship_and_the_design: Vec<f64>,
    /// `d' sigma⁻¹ k sigma⁻¹ d`, `num_coefs` x `num_coefs`, which the
    /// solve against the design weighted by the covariance overwrites.
    of_the_design: Vec<f64>,
    /// The kinship times the working trait through the projection, solved
    /// against the covariance, one value per tested individual.
    through_the_covariance: Vec<f64>,
    /// The design against that same vector before it was solved,
    /// `d' sigma⁻¹ k pw`, one value per column of the design, which the
    /// solve against the design weighted by the covariance overwrites.
    of_the_coefs: Vec<f64>,
    /// What the design explains of it, one value per tested individual.
    explained: Vec<f64>,
}

impl TheStepOnTheVariance {
    /// The buffers of a study of `num_individuals` individuals over a
    /// design of `num_coefs` columns.
    fn of(num_individuals: usize, num_coefs: usize) -> TheStepOnTheVariance {
        // Every product below is of counts a slice already in memory
        // holds: the kinship is `num_individuals` times itself and the
        // design is `num_individuals` times `num_coefs`, with
        // `num_individuals` the larger of the two by 2 at least, since
        // `Design::of_the_study` refuses a study of no more individuals
        // than the columns of its design plus one.
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "the kinship holds `num_individuals` times itself values in one \
                      slice and `num_individuals` is `num_coefs` plus 2 at least, so \
                      each of these products is at most that length and fits in a `usize`"
        )]
        let (of_the_individuals, of_the_design, of_the_coefs) = (
            num_individuals * num_individuals,
            num_individuals * num_coefs,
            num_coefs * num_coefs,
        );
        TheStepOnTheVariance {
            of_the_trace: vec![0.0_f64; of_the_individuals],
            of_the_kinship_and_the_design: vec![0.0_f64; of_the_design],
            of_the_design: vec![0.0_f64; of_the_coefs],
            through_the_covariance: vec![0.0_f64; num_individuals],
            of_the_coefs: vec![0.0_f64; num_coefs],
            explained: vec![0.0_f64; num_individuals],
        }
    }

    /// The buffer the trace was taken in, which is as large as the
    /// projection matrix and is what that matrix is then built in.
    fn take_the_buffer_of_the_trace(&mut self) -> Vec<f64> {
        std::mem::take(&mut self.of_the_trace)
    }

    /// The first derivative of the restricted likelihood at
    /// `genetic_variance` and the average information there, from the
    /// linearization that has just settled.
    ///
    /// The Newton step is the first over the second.
    ///
    /// # Errors
    ///
    /// [`Error::GwasLinalg`] when one of the products, one of the solves
    /// or the triangular solve of the trace could not be done.
    fn at(&mut self, fitted: &TheLinearization<'_>, genetic_variance: f64) -> Result<(f64, f64)> {
        let num_individuals = fitted.num_individuals;
        let num_coefs = fitted.num_coefs;
        let of_the_trace = self.the_trace_of(fitted, genetic_variance)?;
        // `p (k pw)` is the covariance solved against that vector less
        // what the design explains of it there, which is the projection
        // matrix applied without forming it.
        self.through_the_covariance
            .copy_from_slice(&fitted.of_the_kinship);
        popnei_linalg::solve_with_cholesky(
            &fitted.covariance,
            num_individuals,
            &mut self.through_the_covariance,
            1,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "solve of the kinship times the projected working trait against the \
                        covariance",
            source,
        })?;
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &fitted.of_the_covariance,
                rows: num_coefs,
            },
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &fitted.of_the_kinship,
                cols: 1,
            },
            &mut self.of_the_coefs,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design solved against the covariance with the kinship \
                        times the projected working trait",
            source,
        })?;
        popnei_linalg::solve_with_cholesky(&fitted.dsd, num_coefs, &mut self.of_the_coefs, 1)
            .map_err(|source| Error::GwasLinalg {
                operation: "solve of the design weighted by the covariance against the kinship \
                            times the projected working trait",
                source,
            })?;
        popnei_linalg::product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &fitted.of_the_covariance,
                rows: num_individuals,
            },
            num_coefs,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.of_the_coefs,
                cols: 1,
            },
            &mut self.explained,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design solved against the covariance with what it \
                        explains of the kinship times the projected working trait",
            source,
        })?;
        let information = 0.5
            * fitted
                .of_the_kinship
                .iter()
                .zip(&self.through_the_covariance)
                .zip(&self.explained)
                .map(|((of_the_kinship, solved), explained)| of_the_kinship * (solved - explained))
                .sum::<f64>();
        let quadratic = fitted
            .projected_working
            .iter()
            .zip(&fitted.of_the_kinship)
            .map(|(projected, of_the_kinship)| projected * of_the_kinship)
            .sum::<f64>();
        Ok((0.5 * (quadratic - of_the_trace), information))
    }

    /// The trace of the projection matrix times the kinship at
    /// `genetic_variance`, without forming either.
    ///
    /// It is
    ///
    /// ```text
    /// trace(p k) = trace(sigma⁻¹ k)
    ///              - trace((d' sigma⁻¹ d)⁻¹ d' sigma⁻¹ k sigma⁻¹ d)
    /// ```
    ///
    /// whose second term is a `num_coefs` by `num_coefs` matrix and costs
    /// nothing. The first looks as though it needs every entry of the
    /// covariance's inverse, and it does not: `tau k = sigma - w⁻¹` with
    /// `tau` the variance and `w⁻¹` the weights on the diagonal, so
    ///
    /// ```text
    /// trace(sigma⁻¹ k) = (n - trace(sigma⁻¹ w⁻¹)) / tau
    /// ```
    ///
    /// and `trace(sigma⁻¹ w⁻¹)` is the sum of the squares of the entries
    /// of `l⁻¹ w^-1/2`, the Cholesky factor of the covariance solved
    /// against with one right hand side for each individual, the lower
    /// half of it being the half a Cholesky fills.
    ///
    /// **That identity divides by `tau`, and `tau` is 0 exactly at the
    /// boundary where the kinship explains nothing.** It is the one
    /// division of this fit whose denominator reaches 0 rather than nearly,
    /// so the case is taken before the division and not after: the
    /// covariance is then the weights alone, `sigma⁻¹` is the weights
    /// themselves, and the trace is their sum against the diagonal of the
    /// kinship.
    ///
    /// # Errors
    ///
    /// [`Error::GwasLinalg`] when the triangular solve, one of the two
    /// products or the solve against the design weighted by the covariance
    /// could not be done.
    fn the_trace_of(
        &mut self,
        fitted: &TheLinearization<'_>,
        genetic_variance: f64,
    ) -> Result<f64> {
        let num_individuals = fitted.num_individuals;
        let num_coefs = fitted.num_coefs;
        let of_the_covariance = match genetic_variance > 0.0 {
            true => {
                for (at, (row, weight)) in self
                    .of_the_trace
                    .chunks_exact_mut(num_individuals.max(1))
                    .zip(&fitted.weights)
                    .enumerate()
                {
                    row.fill(0.0);
                    if let Some(value) = row.get_mut(at) {
                        // The reciprocal of the weight is what went on the
                        // diagonal of the covariance, and this is the
                        // square root of that same number and not the
                        // reciprocal of the square root of the weight,
                        // which rounds elsewhere. At a variance of 0 the
                        // covariance is the weights alone, its factor is
                        // the square roots of their reciprocals, and every
                        // quotient of the solve below is then exactly 1,
                        // so the trace is exactly the tested individuals
                        // and the derivative of a kinship that says
                        // nothing is exactly 0 by construction.
                        *value = (1.0 / weight).sqrt();
                    }
                }
                popnei_linalg::solve_triangular(
                    &fitted.covariance,
                    num_individuals,
                    popnei_linalg::TheHalfThatHoldsTheMatrix::TheLowerHalf,
                    &mut self.of_the_trace,
                    num_individuals,
                )
                .map_err(|source| Error::GwasLinalg {
                    operation: "solve of the weights against the factor of the covariance of \
                                the working trait",
                    source,
                })?;
                let of_the_weights = self
                    .of_the_trace
                    .iter()
                    .map(|value| value * value)
                    .sum::<f64>();
                (num_individuals as f64 - of_the_weights) / genetic_variance
            }
            false => fitted
                .weights
                .iter()
                .zip(fitted.kinship.chunks_exact(num_individuals.max(1)))
                .enumerate()
                .map(|(at, (weight, row))| weight * row.get(at).copied().unwrap_or(f64::NAN))
                .sum::<f64>(),
        };
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: fitted.kinship,
                rows: num_individuals,
            },
            num_individuals,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &fitted.of_the_covariance,
                cols: num_coefs,
            },
            &mut self.of_the_kinship_and_the_design,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the kinship with the design solved against the covariance",
            source,
        })?;
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &fitted.of_the_covariance,
                rows: num_coefs,
            },
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.of_the_kinship_and_the_design,
                cols: num_coefs,
            },
            &mut self.of_the_design,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design solved against the covariance with the kinship \
                        times itself",
            source,
        })?;
        // The solve reads one right hand side per row, and the rows of
        // `d' sigma⁻¹ k sigma⁻¹ d` are its columns, the matrix being
        // symmetric, so its row `j` solved against the design weighted by
        // the covariance holds the entry `j, j` of the product at `j`.
        popnei_linalg::solve_with_cholesky(
            &fitted.dsd,
            num_coefs,
            &mut self.of_the_design,
            num_coefs,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "solve of the design weighted by the covariance against the kinship \
                        between its columns",
            source,
        })?;
        let of_the_design = self
            .of_the_design
            .chunks_exact(num_coefs.max(1))
            .enumerate()
            .map(|(at, row)| row.get(at).copied().unwrap_or(f64::NAN))
            .sum::<f64>();
        Ok(of_the_covariance - of_the_design)
    }
}

/// The inverse of a matrix a Cholesky has factored, counted where it is
/// formed.
///
/// It is the one route this module has to an inverse, so `inverses` is
/// what the fit formed and not a number written down beside it, which is
/// what the doc comment of [`LogisticMixedModel::inverses`] asks of the
/// three counters: a count taken anywhere else could say the fit inverts
/// as "How popnei fits it, and why not pyNei's way" of
/// `docs/specs/gwas.md` describes while the code did otherwise.
///
/// # Errors
///
/// [`Error::GwasLinalg`] when the inverse could not be formed.
fn the_inverse_that_is_counted(
    covariance: &[f64],
    num_individuals: usize,
    into: &mut [f64],
    inverses: &mut usize,
) -> Result<()> {
    *inverses = inverses.saturating_add(1);
    popnei_linalg::invert_with_cholesky(covariance, num_individuals, into).map_err(|source| {
        Error::GwasLinalg {
            operation: "inverse of the covariance of the working trait",
            source,
        }
    })
}

/// The projection matrix of a linearization that has settled,
/// `num_individuals` x `num_individuals`, row after row.
///
/// It is `p = sigma⁻¹ - sigma⁻¹ d (d' sigma⁻¹ d)⁻¹ d' sigma⁻¹`, and it is
/// the one place of the fit where the inverse of an individuals by
/// individuals matrix is formed, because the score test of a variant wants
/// the matrix and not a solve against it. `of_the_design` is the buffer
/// the part the design explains is formed in, which the step on the
/// variance is done with, and `inverses` is the fit's count of the
/// inverses it has formed, which this adds its own to at the call itself.
///
/// The inverse comes back with only its lower half written, the inverse of
/// a symmetric matrix being symmetric, so that half is mirrored into the
/// other before the design is taken out of it.
///
/// # Errors
///
/// [`Error::GwasLinalg`] when the inverse, the solve or the product could
/// not be done.
fn the_projection_of(
    fitted: &TheLinearization<'_>,
    of_the_design: Vec<f64>,
    inverses: &mut usize,
) -> Result<Vec<f64>> {
    let num_individuals = fitted.num_individuals;
    let num_coefs = fitted.num_coefs;
    let mut projection = vec![0.0_f64; fitted.covariance.len()];
    the_inverse_that_is_counted(
        &fitted.covariance,
        num_individuals,
        &mut projection,
        inverses,
    )?;
    let mut rest: &mut [f64] = &mut projection;
    for at in 0..num_individuals {
        let taken = std::mem::take(&mut rest);
        let Some((row, tail)) = taken.split_at_mut_checked(num_individuals) else {
            break;
        };
        for (value, other) in row
            .iter_mut()
            .skip(at.saturating_add(1))
            .zip(tail.chunks_exact(num_individuals.max(1)))
        {
            *value = other.get(at).copied().unwrap_or(f64::NAN);
        }
        rest = tail;
    }
    // The covariance solved against the design is held one row per column
    // of the design, and the solve below takes one right hand side per
    // row, so it is laid out one row per individual first.
    let mut solved = vec![0.0_f64; fitted.design.values().len()];
    for (at, column) in fitted
        .of_the_covariance
        .chunks_exact(num_individuals.max(1))
        .enumerate()
    {
        for (row, value) in solved.chunks_exact_mut(num_coefs.max(1)).zip(column) {
            if let Some(entry) = row.get_mut(at) {
                *entry = *value;
            }
        }
    }
    popnei_linalg::solve_with_cholesky(&fitted.dsd, num_coefs, &mut solved, num_individuals)
        .map_err(|source| Error::GwasLinalg {
            operation: "solve of the design weighted by the covariance against the design \
                        solved against it",
            source,
        })?;
    let mut of_the_design = of_the_design;
    popnei_linalg::product(
        TheFirstOperand::ByTheValuesSummedOver {
            values: &fitted.of_the_covariance,
            rows: num_individuals,
        },
        num_coefs,
        TheSecondOperand::ByTheColumnsOfTheResult {
            values: &solved,
            cols: num_individuals,
        },
        &mut of_the_design,
    )
    .map_err(|source| Error::GwasLinalg {
        operation: "product that takes the design out of the covariance's inverse",
        source,
    })?;
    for (entry, explained) in projection.iter_mut().zip(&of_the_design) {
        *entry -= *explained;
    }
    Ok(projection)
}

/// The linearization of the logistic mixed model on the two panels of
/// `docs/specs/gwas.md`, at the variance of the kinship effect that GMMAT
/// fitted and at the two ends where it has none and where it has more than
/// any fit reaches.
#[cfg(test)]
mod glmm {
    use std::cmp::Ordering;

    use super::{
        LogisticMixedModel, TheBracket, TheLinearization, TheStepOnTheVariance, the_projection_of,
    };
    use crate::block::{BlockReader, Reblock};
    use crate::error::Error;
    use crate::gwas::dosages::{BlockOfThePass, GwasDosages};
    use crate::gwas::linear::lm::{
        THE_HEADER_OF_EIGHT, reader_over, the_study_of, the_trait_and_the_design_of_the_panel,
    };
    use crate::gwas::linear_mixed::lmm::{
        the_individuals_of_the_kinship, the_kinship_of, the_row_of, the_vcf_of_the_panel,
    };
    use crate::gwas::logistic::LogisticModel;
    use crate::gwas::result::Gwas;
    use crate::gwas::study::{Design, GwasInput, GwasModel, TestType, TraitType};
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::Needs;

    /// The variance of the random effect of the kinship that GMMAT 1.5.0's
    /// `glmmkin` fitted for the panel with every genotype called, from
    /// `tests/reference/gwas/gmmat.null_models.tsv`, as
    /// "The logistic mixed model" of `docs/specs/gwas.md` prints it.
    ///
    /// The whole fit lands on its own value of this, which is what the
    /// step on the variance and the tests of the next commit are about.
    /// What these tests hold it at is this number, so that what they read
    /// is the linearization alone.
    const OF_GMMAT_VARIANCE: f64 = 1.508057;

    /// The effects of the intercept, of `cov1` and of `cov2` that the same
    /// fit gives, and how far from them a linearization held at the
    /// variance above may land: 1e-5 absolute, which is the bound "The
    /// logistic mixed model" of `docs/specs/gwas.md` sets on the null
    /// model of this model.
    ///
    /// The three numbers are printed to seven digits, so 5e-7 of the bound
    /// is the printing and the rest is the distance between two fits.
    /// Holding the variance at those seven digits is another 5e-7 of it:
    /// GMMAT's effects are the ones at its own value of the variance,
    /// which the file gives rounded. Measured over the three on 24
    /// September 2026, the worst is 1.97e-6 on both backends, `cov2`,
    /// which is 20 per cent of what is allowed.
    const OF_GMMAT_EFFECTS: ([f64; 3], f64) = ([-1.416464, 0.753476, 1.583210], 1e-5);

    /// How far the variance of the kinship effect the whole fit lands on
    /// may be from GMMAT's [`OF_GMMAT_VARIANCE`]: 1e-5 absolute, which is
    /// the bound "How it is verified" of "The logistic mixed model" of
    /// `docs/specs/gwas.md` sets on the null model of this model.
    ///
    /// It measures the distance between two fits and not popnei's own
    /// arithmetic, so it is the spec's number and is not lowered to two or
    /// three times what popnei reaches:
    /// `docs/reports/glmm-method/README.md` measured pyNei's fit and the
    /// cheaper one alike at 6.3e-6 of GMMAT's variance, and popnei gives
    /// 6.30e-6 on both backends, 63 per cent of what is allowed. What
    /// popnei's own arithmetic is worth here is the distance between the
    /// two backends, 7.4e-16 of the variance, ten orders below the
    /// distance from GMMAT.
    const OF_GMMAT_GENETIC_VARIANCE: f64 = 1e-5;

    /// The variances the trace of the projection matrix times the kinship
    /// is read at, and how far the identity it is taken from may be from
    /// the two matrices multiplied out at each, as a share of the absolute
    /// terms of that sum.
    ///
    /// The trace is a sum of 40000 products of both signs that cancel to
    /// about half of what went into them, so what it is measured against
    /// is those products' absolute values and not its own value, which is
    /// what "How it is verified" of `docs/specs/gwas.md` asks a tolerance
    /// to be measured against.
    ///
    /// One bound does not hold over the range the search walks through.
    /// The identity is `(n - trace(sigma⁻¹ w⁻¹)) / tau` and the difference
    /// of its two terms goes to 0 with the variance while the rounding of
    /// each does not, so the error grows as one over the variance:
    /// measured over both panels on both backends on 24 September 2026,
    /// the worst share is 1.45e-7 at a variance of 1e-6, 7.7e-11 at 1e-4,
    /// 2.7e-12 at 1e-2, 1.1e-13 at 0.1 and 2.1e-14 at GMMAT's 1.508057. A
    /// variance of 0 is not on that curve: the division is not taken
    /// there, and the worst is 6.3e-15. So the 6e-14 this module held the
    /// trace to is the bound at about a variance of 1 and above, and each
    /// variance here carries its own, between 2.8 and 3.6 times the worst
    /// measured at it.
    ///
    /// The fit reads the trace at every step it takes, the small
    /// variances among them: the search comes down from above and the
    /// floor it stops at is 1e-6, where the trace is worth 7 digits and
    /// not 14. What that costs is the derivative of the restricted
    /// likelihood there, and the fit's own answer to a variance below the
    /// floor is the boundary, 0, whatever that derivative says.
    const OF_THE_TRACE: [(f64, f64); 6] = [
        (0.0, 2e-14),
        (1e-6, 5e-7),
        (1e-4, 2.5e-10),
        (1e-2, 9e-12),
        (0.1, 4e-13),
        (OF_GMMAT_VARIANCE, 6e-14),
    ];

    /// How large a column of the design against a column of the projection
    /// matrix may be, as a share of the absolute terms of that sum:
    /// 6e-15.
    ///
    /// Measured over the 600 sums of each panel on 24 September 2026: the
    /// worst is 1.52e-15 on Accelerate, on the panel with every genotype
    /// called, and 1.44e-15 on faer, on the panel with genotypes missing,
    /// so this is 3.9 times where it breaks. The two backends reach their
    /// worst on different panels, which is why each says which.
    const OF_THE_DESIGN_AGAINST_THE_PROJECTION: f64 = 6e-15;

    /// How far the trait less the fitted chance may be from the working
    /// trait through the projection matrix, as a share of the largest
    /// absolute value of that vector: 5e-13.
    ///
    /// The two are equal at the fit's fixed point and not before it, so
    /// what this measures is how near that point the linearization
    /// stopped, which its own rule bounds at 1e-6 of the linear predictor
    /// plus 1 and which the fit comes far nearer than, and not popnei's
    /// arithmetic. It is a share of the residual's own largest value and
    /// not of each value, since a residual near 0 is an individual the fit
    /// explained. Measured over both panels on 24 September 2026: the
    /// worst is 1.172e-13 on Accelerate and 1.175e-13 on faer, both on the
    /// panel with every genotype called, so this is 4.3 times where it
    /// breaks.
    const OF_THE_RESIDUAL_AGAINST_THE_PROJECTED_WORKING_TRAIT: f64 = 5e-13;

    /// How far the effects of a linearization held at a variance of 0 may
    /// be from the plain logistic null's: 5e-15 absolute, where the
    /// largest of the three effects is 1.235.
    ///
    /// It is absolute and not a share of each effect, since the effects of
    /// a fit are near 0 when the covariate does nothing and the bound
    /// would then be reading the rounding of a number that is not there.
    /// Measured on 24 September 2026: the worst is 1.78e-15 on Accelerate
    /// and 2.22e-16 on faer, both at the intercept on both panels, so this
    /// is 2.8 times where it breaks.
    const OF_THE_LOGISTIC_NULL: f64 = 5e-15;

    /// How large the design against the working trait through the
    /// projection matrix may be, as a share of the absolute terms of that
    /// sum: 1.5e-15.
    ///
    /// The entries of `d' p w` are sums of 200 numbers of both signs that
    /// cancel to 0, so what they are measured against is the sum of those
    /// numbers' absolute values and not any value of the result, which is
    /// what "How it is verified" of `docs/specs/gwas.md` asks a tolerance
    /// to be measured against. Measured at GMMAT's variance over both
    /// panels on 24 September 2026: the worst is 3.18e-16 on Accelerate
    /// and 4.25e-16 on faer, both on the panel with genotypes missing, so
    /// this is 3.5 times where it breaks.
    const OF_THE_PROJECTED_WORKING_TRAIT: f64 = 1.5e-15;

    /// How far the covariance applied to the working trait through the
    /// projection may be from that trait less what the design explains of
    /// it, as a share of the largest absolute value of the working trait:
    /// 5e-15.
    ///
    /// Measured at GMMAT's variance over both panels on 24 September 2026:
    /// the worst is 1.80e-15 on Accelerate, on the panel with genotypes
    /// missing, and 1.67e-15 on faer, on the panel with every genotype
    /// called, so this is 2.8 times where it breaks.
    const OF_THE_COVARIANCE_APPLIED_BACK: f64 = 5e-15;

    /// A variance of the kinship effect large enough that the fitted
    /// chance of an individual reaches 0 or 1 and its weight falls to 0,
    /// which no fit of either panel reaches: 1e12, measured on
    /// `panel_called` on 24 September 2026, where 1e10 still settles.
    const A_VARIANCE_THAT_COLLAPSES_THE_WEIGHTS: f64 = 1.0e12;

    /// The round the linearization held at
    /// [`A_VARIANCE_THAT_COLLAPSES_THE_WEIGHTS`] is refused at: 26, which
    /// "The logistic mixed model" of `docs/specs/gwas.md` gives, measured
    /// on `panel_called` on 24 September 2026 on both backends.
    const THE_ROUND_THE_WEIGHTS_COLLAPSE_AT: usize = 26;

    /// A variance large enough that the kinship's own smallest eigenvalue,
    /// -3.4e-15 on `panel_called`, takes the covariance of the working
    /// trait below 0: 1e16.
    ///
    /// A weight is at most 0.25, so the reciprocals put 4 at least on
    /// every diagonal entry of that covariance, and the variance has to
    /// reach 4 over the size of that eigenvalue, 1.2e15, before the
    /// factorization refuses it. That is the arithmetic "A kinship that is
    /// not positive semidefinite" of `docs/specs/gwas.md` gives, and this
    /// is where it lands on the panel whose eigenvalue is the smallest.
    const A_VARIANCE_THAT_IS_NOT_A_COVARIANCE: f64 = 1.0e16;

    /// A variance of the kinship effect at which the pivots of
    /// `d' sigma⁻¹ d` collapse over a kinship of all ones, which no fit
    /// reaches: 1e13.
    ///
    /// A kinship of all ones asks for a random effect that is one number
    /// for every individual, which is what the intercept is, so the
    /// covariance takes that direction out of the design more and more
    /// sharply as the variance grows and the smallest pivot of the design
    /// weighted by it falls to nothing against the largest. Measured on
    /// `panel_called` on 24 September 2026 on both backends: 1e10 still
    /// settles, 1e11 runs out of rounds, 1e12 to 1e15 collapse the pivots
    /// at the first round, and at 1e16 on faer and 1e17 on Accelerate the
    /// covariance can no longer be factored at all.
    const A_VARIANCE_THAT_COLLAPSES_THE_PIVOTS: f64 = 1.0e13;

    /// A variance of the kinship effect at which a linearization of the
    /// eight individuals of the two families over a kinship of all ones
    /// never settles: 1e13.
    ///
    /// It is the same collapse as [`A_VARIANCE_THAT_COLLAPSES_THE_PIVOTS`]
    /// one step earlier: the linear predictor still moves by more than
    /// 1e-6 of its own largest absolute value plus 1 at every round, while
    /// the pivots have not yet fallen far enough for the rule of **Open 5**
    /// of `docs/specs/gwas.md` to stop it. Measured on 24 September 2026 on
    /// both backends: 1e11 settles in 2 rounds on Accelerate and 4 on faer,
    /// 1e12 to 1e14 run the 200 rounds out, and 1e15 collapses the pivots
    /// at the first round.
    const A_VARIANCE_A_LINEARIZATION_NEVER_SETTLES_AT: f64 = 1.0e13;

    /// What is added to the diagonal of a kinship of all ones to make one
    /// the search over the variance never settles on: 1e-14.
    ///
    /// A kinship of all ones on its own is a covariance the two backends
    /// disagree about at the variances the search walks to, faer refusing
    /// its Cholesky where Accelerate accepts it. This much on the diagonal
    /// is enough for both to factor it and far too little for the trait to
    /// tell the random effect from the intercept, so the derivative asks
    /// for a larger variance for ever. Measured on `panel_called` on 24
    /// September 2026: 1e-13, 1e-14 and 1e-15 run the 200 steps out on both
    /// backends, 1e-12 runs the rounds of a linearization out instead, and
    /// at 1e-16 and below faer refuses the covariance.
    const THE_IDENTITY_IN_A_KINSHIP_OF_ALL_ONES: f64 = 1.0e-14;

    /// The trait `binom` and the covariates `cov1` and `cov2` of
    /// `tests/reference/gwas/phenotypes.csv` for the individuals of the
    /// kinship plink2 wrote for `name`, with that kinship.
    fn the_panel(name: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let individuals = the_individuals_of_the_kinship(name);
        assert_eq!(individuals.len(), 200, "the individuals of {name}");
        let kinship = the_kinship_of(name);
        assert_eq!(kinship.len(), 40000, "the entries of the kinship of {name}");
        let (phenotype, design) =
            the_trait_and_the_design_of_the_panel(&individuals, TraitType::Binomial);
        (phenotype, design, kinship)
    }

    /// The design of one of the two panels, checked, and the plain logistic
    /// null the mixed fit starts from.
    fn the_design_and_the_null_of<'a>(
        phenotype: &'a [f64],
        values: &'a [f64],
        kinship: &'a [f64],
        tested: &'a [usize],
    ) -> (Design<'a>, LogisticModel) {
        let study = GwasInput {
            phenotype,
            trait_type: TraitType::Binomial,
            design: values,
            num_coefs: 3,
            kinship: Some(kinship),
            test: Some(TestType::Score),
            use_grammar_gamma_approx: false,
            individuals: tested,
            transform_to_biallelic: false,
        };
        let design = match Design::of_the_study(&study, tested.len()) {
            Ok(design) => design,
            Err(error) => panic!("the design of the panel: {error}"),
        };
        let null = match LogisticModel::of_the_study(phenotype, &design, TestType::Score) {
            Ok(fitted) => fitted,
            Err(error) => panic!("the logistic null of the panel: {error}"),
        };
        (design, null)
    }

    /// At a variance of 0 the kinship explains nothing, the covariance of
    /// the working trait is the weights alone, and the linearization is
    /// the iteratively reweighted least squares that fitted the plain
    /// logistic null: it is already at its own fixed point, so it settles
    /// in one round and gives that null's effects back.
    ///
    /// The effects it is asserted against are the ones
    /// [`LogisticModel::of_the_study`] fitted, which "The logistic model"
    /// of `docs/specs/gwas.md` checks against R's `glm`, and not a second
    /// copy of the same arithmetic: the null is fitted by solving
    /// `d' w d` against the design times the residuals, and this by
    /// factoring an individuals by individuals covariance and solving the
    /// design against it.
    ///
    /// What separates the two is the one further round this takes, whose
    /// step the null had already brought below 1e-8, and the rounding of
    /// two different routes to the same numbers.
    /// [`OF_THE_LOGISTIC_NULL`] is what they are held to.
    #[test]
    fn the_linearization_at_no_kinship_variance_is_the_plain_logistic_fit() {
        let (phenotype, values, kinship) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let (design, null) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
        let of_the_null = null.coefs().to_vec();
        let mut fitted =
            match TheLinearization::of_the_logistic_null(&phenotype, &design, &kinship, &null) {
                Ok(fitted) => fitted,
                Err(error) => panic!("the linearization of the panel: {error}"),
            };
        if let Err(error) = fitted.at(0.0) {
            panic!("the linearization at a variance of 0: {error}");
        }
        assert_eq!(
            fitted.rounds, 1,
            "the rounds the linearization at a variance of 0 took"
        );
        let allowed = OF_THE_LOGISTIC_NULL;
        for (at, (found, expected)) in fitted.coefs.iter().zip(&of_the_null).enumerate() {
            let difference = (found - expected).abs();
            assert!(
                difference <= allowed,
                "the effect {at} at a variance of 0 is {found} and the logistic null gives \
                 {expected}, {difference} away, against the {allowed} allowed"
            );
        }
    }

    /// Held at the variance GMMAT fitted, the linearization gives GMMAT's
    /// three covariate effects, within the 1e-5 absolute that "The
    /// logistic mixed model" of `docs/specs/gwas.md` sets on the null
    /// model of this model.
    ///
    /// It is the linearization alone that this reads: the variance is held
    /// at GMMAT's and not searched for, and the search is the next commit.
    #[test]
    fn the_linearization_at_gmmats_variance_is_gmmats_three_covariate_effects() {
        let (phenotype, values, kinship) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let (design, null) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
        let mut fitted =
            match TheLinearization::of_the_logistic_null(&phenotype, &design, &kinship, &null) {
                Ok(fitted) => fitted,
                Err(error) => panic!("the linearization of the panel: {error}"),
            };
        if let Err(error) = fitted.at(OF_GMMAT_VARIANCE) {
            panic!("the linearization at GMMAT's variance: {error}");
        }
        let (of_gmmat, allowed) = OF_GMMAT_EFFECTS;
        assert_eq!(fitted.coefs.len(), of_gmmat.len(), "the effects fitted");
        for (at, (found, expected)) in fitted.coefs.iter().zip(&of_gmmat).enumerate() {
            let difference = (found - expected).abs();
            assert!(
                difference <= allowed,
                "the effect {at} at GMMAT's variance is {found} and GMMAT gives {expected}, \
                 {difference} away, against the {allowed} allowed"
            );
        }
    }

    /// The projection matrix takes the covariates out of what it is
    /// applied to, so the design against the working trait through it is 0
    /// on both panels, and the covariance applied to that same vector
    /// gives the working trait less what the design explains of it.
    ///
    /// The two are the identities the solves are supposed to satisfy, read
    /// by multiplying back: `d' p w` is 0 and `sigma p w` is
    /// `w - d coefs`. Neither is another way of computing what the code
    /// computed, which is why they can fail.
    ///
    /// Each is bounded by a share of its own scale and not of each value,
    /// and each has its own bound with its own measurement:
    /// [`OF_THE_PROJECTED_WORKING_TRAIT`] and
    /// [`OF_THE_COVARIANCE_APPLIED_BACK`].
    #[test]
    fn the_projection_takes_the_design_out_of_the_working_trait() {
        for name in ["panel_called", "panel"] {
            let (phenotype, values, kinship) = the_panel(name);
            let tested: Vec<usize> = (0..phenotype.len()).collect();
            let (design, null) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
            let mut fitted = match TheLinearization::of_the_logistic_null(
                &phenotype, &design, &kinship, &null,
            ) {
                Ok(fitted) => fitted,
                Err(error) => panic!("the linearization of {name}: {error}"),
            };
            if let Err(error) = fitted.at(OF_GMMAT_VARIANCE) {
                panic!("the linearization of {name} at GMMAT's variance: {error}");
            }
            let num_individuals = phenotype.len();
            let num_coefs = 3;
            for column in 0..num_coefs {
                let mut against = 0.0_f64;
                let mut scale = 0.0_f64;
                for (row, projected) in values
                    .chunks_exact(num_coefs)
                    .zip(&fitted.projected_working)
                {
                    let value = row[column] * projected;
                    against += value;
                    scale += value.abs();
                }
                assert!(
                    against.abs() <= OF_THE_PROJECTED_WORKING_TRAIT * scale,
                    "{name}: the column {column} of the design against the working trait \
                     through the projection is {against}, where the absolute terms of that \
                     sum are {scale} and {OF_THE_PROJECTED_WORKING_TRAIT} of them is allowed"
                );
            }
            let of_the_working = fitted
                .working_trait
                .iter()
                .fold(0.0_f64, |largest, value| largest.max(value.abs()));
            for individual in 0..num_individuals {
                // The covariance of the working trait applied to that
                // trait through the projection, row by row, which the
                // linearization never forms: it holds the factorization of
                // the covariance where this reads the covariance itself.
                let mut applied = 0.0_f64;
                for (other, projected) in fitted.projected_working.iter().enumerate() {
                    let of_the_pair = kinship[individual * num_individuals + other];
                    let mut value = OF_GMMAT_VARIANCE * of_the_pair;
                    if other == individual {
                        value += 1.0 / fitted.weights[individual];
                    }
                    applied += value * projected;
                }
                let explained: f64 = values[individual * num_coefs..]
                    .iter()
                    .take(num_coefs)
                    .zip(&fitted.coefs)
                    .map(|(value, coef)| value * coef)
                    .sum();
                let difference = (applied - (fitted.working_trait[individual] - explained)).abs();
                assert!(
                    difference <= OF_THE_COVARIANCE_APPLIED_BACK * of_the_working,
                    "{name}: the covariance applied to the working trait through the \
                     projection is {applied} for the individual {individual}, where the \
                     working trait less what the design explains of it is \
                     {expected}, {difference} away, against the \
                     {OF_THE_COVARIANCE_APPLIED_BACK} of the working trait's largest value \
                     {of_the_working} allowed",
                    expected = fitted.working_trait[individual] - explained
                );
            }
        }
    }

    /// A variance of the kinship effect large enough to take the fitted
    /// chance of an individual to 0 or 1 leaves that individual with a
    /// weight of 0 and a working trait that is not a number, and the
    /// linearization is refused as a fit that did not settle, with the
    /// round it reached.
    ///
    /// No fit reaches such a variance: on `panel_called` this happens at
    /// 1e12, where the fit lands at 1.508 and where 1e10 still settles,
    /// with a smallest weight of 3.3e-24. What the test guards is the
    /// error the user gets there, which without it is the linear algebra
    /// crate refusing a matrix that is not finite, a `RuntimeError` in
    /// Python and so a defect of popnei, for the data.
    #[test]
    fn a_linearization_whose_weights_have_collapsed_is_refused() {
        let (phenotype, values, kinship) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let (design, null) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
        let mut fitted =
            match TheLinearization::of_the_logistic_null(&phenotype, &design, &kinship, &null) {
                Ok(fitted) => fitted,
                Err(error) => panic!("the linearization of the panel: {error}"),
            };
        match fitted.at(A_VARIANCE_THAT_COLLAPSES_THE_WEIGHTS) {
            Err(Error::GwasFitDidNotSettle { model, rounds }) => {
                assert_eq!(
                    model,
                    crate::gwas::GwasModel::Glmm,
                    "the model the refusal names"
                );
                assert_eq!(
                    rounds, THE_ROUND_THE_WEIGHTS_COLLAPSE_AT,
                    "the round the refusal names"
                );
            }
            other => panic!("the linearization at a variance of 1e12 gave {other:?}"),
        }
    }

    /// A variance large enough for the kinship's own smallest eigenvalue
    /// to take the covariance of the working trait below 0 is refused with
    /// an error that names the kinship, and not as a fit that did not
    /// settle: what the user does about it is look at the matrix they
    /// brought, which missing genotypes can leave with an eigenvalue below
    /// 0.
    #[test]
    fn a_covariance_that_is_not_positive_definite_names_the_kinship() {
        let (phenotype, values, kinship) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let (design, null) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
        let mut fitted =
            match TheLinearization::of_the_logistic_null(&phenotype, &design, &kinship, &null) {
                Ok(fitted) => fitted,
                Err(error) => panic!("the linearization of the panel: {error}"),
            };
        match fitted.at(A_VARIANCE_THAT_IS_NOT_A_COVARIANCE) {
            Err(Error::GwasKinshipNotACovariance { at }) => {
                assert!(at < phenotype.len(), "the row the refusal names, {at}");
            }
            other => panic!("the linearization at a variance of 1e16 gave {other:?}"),
        }
    }

    /// A linearization whose weighted design has collapsed is refused with
    /// the round it reached, by the pivot rule of **Open 5** of
    /// `docs/specs/gwas.md`, and not answered with effects solved for out
    /// of a system that is no longer one.
    ///
    /// The fixture is the panel with every genotype called over a kinship
    /// of all ones, held at
    /// [`A_VARIANCE_THAT_COLLAPSES_THE_PIVOTS`]. What the assertion on the
    /// weights says is that this is the pivot rule and not the other cause
    /// of the same refusal: a weight that has fallen to 0 leaves a
    /// reciprocal that is not finite, and every one of the 200 is finite
    /// here.
    ///
    /// The refusal is at the first round, so the fit is refused before any
    /// number of it reaches the user.
    #[test]
    fn a_linearization_whose_pivots_have_collapsed_is_refused() {
        let (phenotype, values, _) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let kinship = vec![1.0_f64; phenotype.len() * phenotype.len()];
        let (design, null) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
        let mut fitted =
            match TheLinearization::of_the_logistic_null(&phenotype, &design, &kinship, &null) {
                Ok(fitted) => fitted,
                Err(error) => panic!("the linearization over a kinship of all ones: {error}"),
            };
        match fitted.at(A_VARIANCE_THAT_COLLAPSES_THE_PIVOTS) {
            Err(Error::GwasFitDidNotSettle { model, rounds }) => {
                assert_eq!(model, GwasModel::Glmm, "the model the refusal names");
                assert_eq!(rounds, 1, "the round the refusal names");
            }
            other => panic!(
                "the linearization at a variance of {A_VARIANCE_THAT_COLLAPSES_THE_PIVOTS} \
                 gave {other:?}"
            ),
        }
        for (at, weight) in fitted.weights.iter().enumerate() {
            let of_the_covariance = 1.0 / weight;
            assert!(
                of_the_covariance.is_finite(),
                "the weight of the individual {at} is {weight}, whose reciprocal is \
                 {of_the_covariance}, so what the fit was refused for is that weight and \
                 not the pivots"
            );
        }
    }

    /// A linearization whose linear predictor is still moving after the 200
    /// rounds it is given is refused with that count.
    ///
    /// The fixture is the eight individuals of the two families over a
    /// kinship of all ones, held at
    /// [`A_VARIANCE_A_LINEARIZATION_NEVER_SETTLES_AT`], which is one step
    /// below where the pivots collapse there.
    #[test]
    fn a_linearization_that_is_still_moving_after_its_rounds_is_refused() {
        let kinship = [1.0_f64; 64];
        let (design, null) = the_design_and_the_null_of_eight(&kinship);
        let mut fitted = match TheLinearization::of_the_logistic_null(
            &THE_TRAIT_OF_EIGHT,
            &design,
            &kinship,
            &null,
        ) {
            Ok(fitted) => fitted,
            Err(error) => panic!("the linearization of the eight individuals: {error}"),
        };
        match fitted.at(A_VARIANCE_A_LINEARIZATION_NEVER_SETTLES_AT) {
            Err(Error::GwasFitDidNotSettle { model, rounds }) => {
                assert_eq!(model, GwasModel::Glmm, "the model the refusal names");
                // The 200 is the spec's own number and not
                // [`ROUNDS_OF_A_LINEARIZATION`], which a change of that
                // constant would carry this assertion along with.
                assert_eq!(rounds, 200, "the rounds the refusal names");
            }
            other => panic!(
                "the linearization at a variance of \
                 {A_VARIANCE_A_LINEARIZATION_NEVER_SETTLES_AT} gave {other:?}"
            ),
        }
    }

    /// A search whose variance of the kinship effect is still moving after
    /// the 200 steps it is given is refused with that count.
    ///
    /// The fixture is the panel with every genotype called over a kinship
    /// of all ones with [`THE_IDENTITY_IN_A_KINSHIP_OF_ALL_ONES`] on its
    /// diagonal, which is symmetric and positive definite and which the
    /// trait cannot tell from the intercept: the derivative asks for a
    /// larger variance at every step and the fit has no place to stop.
    /// Nothing else here reaches this loop, because a variance that walks
    /// away takes the linearization at it past its own rounds first.
    /// Measured on both backends on 24 September 2026: the 200 steps are
    /// the search's and every linearization in them settles.
    #[test]
    fn a_search_that_is_still_moving_after_its_steps_is_refused() {
        let (phenotype, values, _) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let num_individuals = phenotype.len();
        let mut kinship = vec![1.0_f64; num_individuals * num_individuals];
        for (at, row) in kinship.chunks_exact_mut(num_individuals).enumerate() {
            if let Some(value) = row.get_mut(at) {
                *value += THE_IDENTITY_IN_A_KINSHIP_OF_ALL_ONES;
            }
        }
        let (design, _) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
        match LogisticMixedModel::of_the_study(&phenotype, &design, &kinship) {
            Err(Error::GwasFitDidNotSettle { model, rounds }) => {
                assert_eq!(model, GwasModel::Glmm, "the model the refusal names");
                // The 200 is the spec's own number and not
                // [`STEPS_ON_THE_VARIANCE`], for the reason above.
                assert_eq!(rounds, 200, "the steps the refusal names");
            }
            Err(error) => panic!("the fit over a kinship of all ones gave {error}"),
            Ok(model) => panic!(
                "the fit over a kinship of all ones gave a variance of {}",
                model.genetic_variance
            ),
        }
    }

    /// A Newton step that would take the variance out of the bracket is
    /// replaced by the geometric mean of the bracket's two ends, which is
    /// what keeps the search from cycling.
    ///
    /// `docs/reports/glmm-method/README.md` records that the first fit
    /// written for it, without this, did not converge at all. No fixture
    /// runs that arm: over the whole core suite the step was kept 51 times
    /// and replaced none, measured on 24 September 2026. The ends here are
    /// 1 and 4, so their geometric mean is 2, and a step of 10 from a
    /// variance of 2 would give 12.
    #[test]
    fn a_step_that_would_leave_the_bracket_is_the_geometric_mean_of_its_ends() {
        let mut bracket = TheBracket::of_a_search_that_has_not_started();
        // a derivative above 0 asks for a larger variance, one below it for
        // a smaller one
        bracket.the_end_that(1.0, 0.5);
        bracket.the_end_that(4.0, -0.5);
        for step in [10.0_f64, -10.0] {
            let found = bracket.the_variance_after(2.0, step);
            assert_eq!(
                found.total_cmp(&2.0),
                Ordering::Equal,
                "a step of {step} from a variance of 2, out of the bracket of 1 and 4, gave \
                 {found}"
            );
        }
        let found = bracket.the_variance_after(2.0, 0.5);
        assert_eq!(
            found.total_cmp(&2.5),
            Ordering::Equal,
            "a step of 0.5 from a variance of 2, which stays inside the bracket of 1 and 4, \
             gave {found}"
        );
    }

    /// The whole fit of the panel with every genotype called is GMMAT's:
    /// the variance of the kinship effect and the three covariate effects
    /// within the 1e-5 absolute of "How it is verified" of "The logistic
    /// mixed model" of `docs/specs/gwas.md`, with no residual variance and
    /// no heritability.
    ///
    /// This is the search and not the linearization: nothing here holds
    /// the variance anywhere, and where it lands is what the trace, the
    /// average information and the bracket decide between them. It is the
    /// cargo test of deliverable 1 of `docs/plans/gwas-logistic.md`, whose
    /// own check is the pytest test that the model reaching both packages
    /// brings.
    #[test]
    fn the_fit_of_the_panel_is_gmmats_variance_and_three_covariate_effects() {
        let (phenotype, values, kinship) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let (design, _) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
        let model = match LogisticMixedModel::of_the_study(&phenotype, &design, &kinship) {
            Ok(model) => model,
            Err(error) => panic!("the fit of the panel: {error}"),
        };
        let null = model.null_model(TestType::Score);
        assert_eq!(null.model, GwasModel::Glmm, "the model that was fitted");
        let found = match null.genetic_variance {
            Some(found) => found,
            None => panic!("the fit gave no variance of the kinship effect"),
        };
        let difference = (found - OF_GMMAT_VARIANCE).abs();
        assert!(
            difference <= OF_GMMAT_GENETIC_VARIANCE,
            "the variance of the kinship effect is {found} and GMMAT gives \
             {OF_GMMAT_VARIANCE}, {difference} away, against the \
             {OF_GMMAT_GENETIC_VARIANCE} allowed"
        );
        assert_eq!(
            null.residual_variance, None,
            "the residual variance of a logistic mixed model"
        );
        assert_eq!(
            null.heritability, None,
            "the heritability of a logistic mixed model"
        );
        let (of_gmmat, allowed) = OF_GMMAT_EFFECTS;
        assert_eq!(
            null.covariate_effects.len(),
            of_gmmat.len(),
            "the effects fitted"
        );
        for (at, (found, expected)) in null.covariate_effects.iter().zip(&of_gmmat).enumerate() {
            let difference = (found - expected).abs();
            assert!(
                difference <= allowed,
                "the effect {at} of the fit is {found} and GMMAT gives {expected}, \
                 {difference} away, against the {allowed} allowed"
            );
        }
    }

    /// The fit forms the inverse of a factorized individuals by
    /// individuals matrix once, at the end, and factors one such matrix
    /// once per round: 1 inverse against 22 factorizations on the panel,
    /// over 8 steps on the variance of the kinship effect.
    ///
    /// It is deliverable 3 of `docs/plans/gwas-logistic.md` and the whole
    /// point of fitting this model popnei's way rather than pyNei's, which
    /// `docs/reports/glmm-method/README.md` measured at about twice the
    /// speed: pyNei inverts once per round where this factors, and a
    /// Cholesky factorization costs a third of an inverse. A fit that
    /// inverted where it should solve would give the same numbers and
    /// nothing else here would notice.
    ///
    /// The three counts are taken at the calls themselves. The 22 and the
    /// 8 are what "The logistic mixed model" of `docs/specs/gwas.md` says
    /// the panel takes, and the panel with genotypes missing takes the
    /// same two. They do not pin the rule that ends a linearization:
    /// dropping the plus 1 from the divisor of that rule leaves this
    /// panel's fit bit-identical on both backends, measured on 24
    /// September 2026, because the largest absolute value of its linear
    /// predictor is 4.19 and the plus 1 moves the divisor by a quarter.
    /// [`the_rounds_of_a_linearization_are_measured_against_the_predictor_plus_one`]
    /// is what pins it, on eight individuals whose predictor is near 0.
    #[test]
    fn the_fit_forms_one_inverse_and_factors_the_covariance_once_a_round() {
        let (phenotype, values, kinship) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let (design, _) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
        let model = match LogisticMixedModel::of_the_study(&phenotype, &design, &kinship) {
            Ok(model) => model,
            Err(error) => panic!("the fit of the panel: {error}"),
        };
        assert_eq!(
            model.inverses, 1,
            "the inverses of a factorized individuals by individuals matrix the fit formed"
        );
        assert_eq!(
            model.linearizations, 22,
            "the times the fit factored the covariance of the working trait"
        );
        assert_eq!(
            model.steps_on_the_variance, 8,
            "the steps the fit took on the variance of the kinship effect"
        );
    }

    /// The trace of the projection matrix times the kinship, which the
    /// step on the variance takes from an identity, is the trace of the
    /// two matrices multiplied out.
    ///
    /// It is the one quantity of the step that seems to need every entry
    /// of the covariance's inverse, and the identity is what lets the fit
    /// never form one: the two sides here are that identity, a triangular
    /// solve against the Cholesky factor and a correction of the size of
    /// the design, against the inverse itself. Neither is the other's
    /// arithmetic, which is why this can fail.
    ///
    /// **It is read at a variance of 0 as well as at GMMAT's**, because
    /// the identity divides by that variance and 0 is exactly where the
    /// kinship explains nothing, which is a value the search reaches and
    /// not one it approaches: the case is taken before the division, and a
    /// fit that took it after would give an infinity or a NaN here. This
    /// is the test that pins that ordering, and taking the division at a
    /// variance of 0 leaves the fit over a kinship of all zeros landing on
    /// the same boundary it lands on now.
    ///
    /// It is read at four variances between them as well, because the
    /// error of the identity grows as one over the variance and one bound
    /// does not hold over the range the search walks through:
    /// [`OF_THE_TRACE`] carries the variances and the bound of each, with
    /// what was measured at each.
    ///
    /// Each trace is measured against the absolute terms of the sum it is
    /// and not against its own value, which is a sum of 40000 products of
    /// both signs.
    #[test]
    fn the_trace_from_the_identity_is_the_trace_of_the_two_matrices() {
        for name in ["panel_called", "panel"] {
            let (phenotype, values, kinship) = the_panel(name);
            let tested: Vec<usize> = (0..phenotype.len()).collect();
            let num_individuals = phenotype.len();
            let (design, null) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
            for (variance, allowed) in OF_THE_TRACE {
                let mut fitted = match TheLinearization::of_the_logistic_null(
                    &phenotype, &design, &kinship, &null,
                ) {
                    Ok(fitted) => fitted,
                    Err(error) => panic!("the linearization of {name}: {error}"),
                };
                if let Err(error) = fitted.at(variance) {
                    panic!("the linearization of {name} at a variance of {variance}: {error}");
                }
                let mut step = TheStepOnTheVariance::of(num_individuals, 3);
                let of_the_identity = match step.the_trace_of(&fitted, variance) {
                    Ok(trace) => trace,
                    Err(error) => {
                        panic!("the trace of {name} at a variance of {variance}: {error}")
                    }
                };
                let buffer = step.take_the_buffer_of_the_trace();
                let projection = match the_projection_of(&fitted, buffer, &mut 0) {
                    Ok(projection) => projection,
                    Err(error) => {
                        panic!("the projection of {name} at a variance of {variance}: {error}")
                    }
                };
                let mut of_the_matrices = 0.0_f64;
                let mut scale = 0.0_f64;
                for (row, of_the_kinship) in projection
                    .chunks_exact(num_individuals)
                    .zip(kinship.chunks_exact(num_individuals))
                {
                    for (entry, of_the_pair) in row.iter().zip(of_the_kinship) {
                        of_the_matrices += entry * of_the_pair;
                        scale += (entry * of_the_pair).abs();
                    }
                }
                let difference = (of_the_identity - of_the_matrices).abs();
                assert!(
                    difference <= allowed * scale,
                    "{name} at a variance of {variance}: the trace from the identity is \
                     {of_the_identity} and the two matrices multiplied out give \
                     {of_the_matrices}, {difference} apart, where the absolute terms of that \
                     sum are {scale}, {share} of them, and {allowed} of them is allowed",
                    share = difference / scale
                );
            }
        }
    }

    /// The residual the score test will be taken against, the trait less
    /// the fitted chance, is the working trait through the projection
    /// matrix.
    ///
    /// That equality is the fit's own optimality condition and holds
    /// nowhere but at its fixed point: the covariance applied to the
    /// working trait through the projection is that trait less what the
    /// design explains of it, and the working trait is the linear
    /// predictor plus the trait less the fitted chance over the weight, so
    /// the two statements together leave the reciprocals of the weights
    /// against the difference of the two vectors and nothing else. It is
    /// why "What it gives" of "The logistic mixed model" of
    /// `docs/specs/gwas.md` can say that the residual `p y` of the score
    /// test is simply the trait minus `mu`, and it is the cheapest
    /// evidence there is that the search reached the place it was looking
    /// for, as `y' p y` is for the linear mixed model.
    ///
    /// The working trait through the projection is not kept by the fit, so
    /// it is taken here from a linearization of its own at the variance
    /// the fit landed on, started from the plain logistic null. A
    /// linearization has one fixed point at a variance, so the two paths
    /// meet there, and how nearly is what
    /// [`OF_THE_RESIDUAL_AGAINST_THE_PROJECTED_WORKING_TRAIT`] holds.
    #[test]
    fn the_residual_of_the_fit_is_the_working_trait_through_the_projection() {
        for name in ["panel_called", "panel"] {
            let (phenotype, values, kinship) = the_panel(name);
            let tested: Vec<usize> = (0..phenotype.len()).collect();
            let (design, null) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
            let model = match LogisticMixedModel::of_the_study(&phenotype, &design, &kinship) {
                Ok(model) => model,
                Err(error) => panic!("the fit of {name}: {error}"),
            };
            let mut fitted = match TheLinearization::of_the_logistic_null(
                &phenotype, &design, &kinship, &null,
            ) {
                Ok(fitted) => fitted,
                Err(error) => panic!("the linearization of {name}: {error}"),
            };
            if let Err(error) = fitted.at(model.genetic_variance) {
                panic!("the linearization of {name} at the variance the fit landed on: {error}");
            }
            let of_the_residual = fitted
                .projected_working
                .iter()
                .fold(0.0_f64, |largest, value| largest.max(value.abs()));
            let allowed = OF_THE_RESIDUAL_AGAINST_THE_PROJECTED_WORKING_TRAIT * of_the_residual;
            for (at, (projected, residual)) in fitted
                .projected_working
                .iter()
                .zip(&model.projected_trait)
                .enumerate()
            {
                let difference = (projected - residual).abs();
                assert!(
                    difference <= allowed,
                    "{name}: the working trait through the projection is {projected} for the \
                     individual {at} and the trait less the fitted chance is {residual}, \
                     {difference} apart, against the {allowed} allowed"
                );
            }
        }
    }

    /// The projection matrix the fit leaves takes the design out of
    /// anything it is applied to: every column of the design against every
    /// column of that matrix is 0 on both panels.
    ///
    /// It is the one matrix of the fit that the inverse is formed for, and
    /// it is what every variant of the study is tested through, so what it
    /// is worth is what the score test of the next commit is worth. The
    /// identity is read column by column and not through the trait,
    /// because a matrix that is wrong in one direction alone would survive
    /// any product with one vector.
    ///
    /// [`OF_THE_DESIGN_AGAINST_THE_PROJECTION`] is what each of the 600
    /// sums is held to, as a share of its own absolute terms.
    #[test]
    fn the_projection_of_the_fit_takes_the_design_out_of_every_individual() {
        for name in ["panel_called", "panel"] {
            let (phenotype, values, kinship) = the_panel(name);
            let tested: Vec<usize> = (0..phenotype.len()).collect();
            let num_individuals = phenotype.len();
            let (design, _) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
            let model = match LogisticMixedModel::of_the_study(&phenotype, &design, &kinship) {
                Ok(model) => model,
                Err(error) => panic!("the fit of {name}: {error}"),
            };
            let num_coefs = design.num_coefs();
            for column in 0..num_coefs {
                for of_the_projection in 0..num_individuals {
                    let mut against = 0.0_f64;
                    let mut scale = 0.0_f64;
                    for (row, of_the_individual) in values
                        .chunks_exact(num_coefs)
                        .zip(model.projection.chunks_exact(num_individuals))
                    {
                        let value = row[column] * of_the_individual[of_the_projection];
                        against += value;
                        scale += value.abs();
                    }
                    assert!(
                        against.abs() <= OF_THE_DESIGN_AGAINST_THE_PROJECTION * scale,
                        "{name}: the column {column} of the design against the column \
                         {of_the_projection} of the projection is {against}, where the \
                         absolute terms of that sum are {scale} and \
                         {OF_THE_DESIGN_AGAINST_THE_PROJECTION} of them is allowed"
                    );
                }
            }
        }
    }

    /// A kinship of all zeros explains nothing of the trait, and the fit
    /// walks down to the boundary and stops there with a variance of
    /// exactly 0.
    ///
    /// What it reads is that the search stops there: a second variance at
    /// 0 is what ends the fit, where a fit that only set it to 0 would run
    /// its 200 steps there. Neither the derivative nor the average
    /// information is above 0 here, so the Newton step is not a number,
    /// and what keeps the variance a number is the bracket: a step that
    /// does not land on a finite variance above 0 quarters the variance
    /// instead. That the trace answers a variance of 0 from the weights
    /// and the diagonal of the kinship before its division is not read
    /// here, although the fit walks through it:
    /// [`the_trace_from_the_identity_is_the_trace_of_the_two_matrices`] is
    /// what fails when that division is taken, and this stays green.
    ///
    /// Measured on 24 September 2026 on both backends: 12 steps and 12
    /// rounds, and the effects are the plain logistic null's, which is
    /// what a fit at a variance of 0 is.
    #[test]
    fn a_kinship_that_explains_nothing_lands_at_a_variance_of_zero() {
        let (phenotype, values, _) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let kinship = vec![0.0_f64; phenotype.len() * phenotype.len()];
        let (design, null) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
        let model = match LogisticMixedModel::of_the_study(&phenotype, &design, &kinship) {
            Ok(model) => model,
            Err(error) => panic!("the fit over a kinship of all zeros: {error}"),
        };
        let found = match model.null_model(TestType::Score).genetic_variance {
            Some(found) => found,
            None => panic!("the fit over a kinship of all zeros gave no variance"),
        };
        assert_eq!(
            found.total_cmp(&0.0),
            Ordering::Equal,
            "the variance of the kinship effect over a kinship of all zeros is {found}"
        );
        let allowed = OF_THE_LOGISTIC_NULL;
        for (at, (found, expected)) in model.coefs.iter().zip(null.coefs()).enumerate() {
            let difference = (found - expected).abs();
            assert!(
                difference <= allowed,
                "the effect {at} over a kinship of all zeros is {found} and the logistic null \
                 gives {expected}, {difference} away, against the {allowed} allowed"
            );
        }
    }

    /// A kinship whose magnitude is too small for the variance that fits
    /// it to be a finite number is answered with the boundary, and not
    /// with the linear algebra crate refusing a matrix the user never saw.
    ///
    /// The average information is quadratic in the kinship, so a kinship of
    /// 1e-165 times the identity leaves it below what a `f64` holds while
    /// the derivative is still a number: the Newton step is then an
    /// infinity, and before this was fixed that infinity became the
    /// variance, the covariance became an infinity times a kinship entry of
    /// 0, and the Cholesky refused a matrix that is not finite. That is
    /// [`Error::GwasLinalg`], a `RuntimeError` in Python and so a defect of
    /// popnei, for a kinship that is finite, symmetric and positive
    /// definite.
    ///
    /// Measured on 24 September 2026 on both backends: the default build
    /// raised at 1e-162 and below and answered a variance of 1.27e159 at
    /// 1e-160, and the search now walks down to the boundary in 12 steps
    /// instead. What it lands on is the plain logistic null, which is what
    /// a variance of 0 means, and the effects are held to
    /// [`OF_THE_LOGISTIC_NULL`] against that null's own.
    #[test]
    fn a_kinship_too_small_for_its_variance_to_be_finite_lands_at_the_boundary() {
        let (phenotype, values, _) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let num_individuals = phenotype.len();
        let mut kinship = vec![0.0_f64; num_individuals * num_individuals];
        for (at, row) in kinship.chunks_exact_mut(num_individuals).enumerate() {
            row[at] = 1.0e-165;
        }
        let (design, null) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
        let model = match LogisticMixedModel::of_the_study(&phenotype, &design, &kinship) {
            Ok(model) => model,
            Err(error) => panic!("the fit over a kinship of 1e-165 times the identity: {error}"),
        };
        let found = model.genetic_variance;
        assert_eq!(
            found.total_cmp(&0.0),
            Ordering::Equal,
            "the variance of the kinship effect over a kinship of 1e-165 times the identity \
             is {found}"
        );
        let allowed = OF_THE_LOGISTIC_NULL;
        for (at, (found, expected)) in model.coefs.iter().zip(null.coefs()).enumerate() {
            let difference = (found - expected).abs();
            assert!(
                difference <= allowed,
                "the effect {at} over a kinship of 1e-165 times the identity is {found} and \
                 the logistic null gives {expected}, {difference} away, against the {allowed} \
                 allowed"
            );
        }
    }

    /// A Newton step that is not a finite number quarters the variance,
    /// which is where a step that would take it to 0 or below goes, and an
    /// infinity is as much a step that is not a step as a NaN is.
    ///
    /// The bracket is what the search leans on instead of refusing such a
    /// step itself, and until this was fixed only a NaN and an infinity
    /// below 0 went to the quarter, while an infinity above 0 passed the
    /// test for a variance above 0 and became the variance.
    #[test]
    fn a_step_that_is_not_a_finite_number_quarters_the_variance() {
        let bracket = TheBracket::of_a_search_that_has_not_started();
        for step in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN, -8.0] {
            let found = bracket.the_variance_after(4.0, step);
            assert_eq!(
                found.total_cmp(&1.0),
                Ordering::Equal,
                "a step of {step} from a variance of 4 gave {found}"
            );
        }
        let found = bracket.the_variance_after(4.0, 3.0);
        assert_eq!(
            found.total_cmp(&7.0),
            Ordering::Equal,
            "a step of 3 from a variance of 4 gave {found}"
        );
    }

    /// A weight that is already 0 at the plain logistic null is refused
    /// naming the round the fit was about to run, and never the round 0
    /// that no fit ran.
    ///
    /// The search reads the working trait once before its first
    /// linearization, for the variance it starts at, and the round counter
    /// is the 0 the logistic null left there, so the message said that the
    /// null model did not settle in the 0 rounds it was fitted in. Six
    /// individuals reach it, with a covariate of 1198.97 for one of them
    /// against values between -17 and 0 for the others: the logistic null
    /// fits that individual a chance of 1 and its weight is 0.
    #[test]
    fn a_weight_that_is_zero_at_the_logistic_null_names_the_first_round() {
        let phenotype = [0.0_f64, 1.0, 1.0, 1.0, 0.0, 0.0];
        let covariate = [-16.95_f64, -0.294, 1198.97, -3.188, -6.947, -0.641];
        let values: Vec<f64> = covariate.iter().flat_map(|value| [1.0, *value]).collect();
        let num_individuals = phenotype.len();
        let mut kinship = vec![0.0_f64; num_individuals * num_individuals];
        for (at, row) in kinship.chunks_exact_mut(num_individuals).enumerate() {
            row[at] = 1.0;
        }
        let tested: Vec<usize> = (0..num_individuals).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Binomial,
            design: &values,
            num_coefs: 2,
            kinship: Some(&kinship),
            test: Some(TestType::Score),
            use_grammar_gamma_approx: false,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        let design = match Design::of_the_study(&study, num_individuals) {
            Ok(design) => design,
            Err(error) => panic!("the design of the six individuals: {error}"),
        };
        match LogisticMixedModel::of_the_study(&phenotype, &design, &kinship) {
            Err(Error::GwasFitDidNotSettle { model, rounds }) => {
                assert_eq!(model, GwasModel::Glmm, "the model the refusal names");
                assert_eq!(rounds, 1, "the round the refusal names");
            }
            Err(error) => panic!("the fit of the six individuals gave {error}"),
            Ok(model) => panic!(
                "the fit of the six individuals gave a variance of {}",
                model.genetic_variance
            ),
        }
    }

    /// A study that brought a kinship, and whose plain logistic null runs
    /// away before the first linearization, is refused naming the logistic
    /// mixed model, which is the model the user asked for.
    ///
    /// The fit starts from that null, which is the same trait and the same
    /// design with no kinship in it, and a fit that does not settle names
    /// the model it was fitting: the user was told that a binomial trait
    /// with no kinship is a logistic regression, which is not the study
    /// they made, and was given the two remedies of a fit that has no
    /// kinship instead of the three of this one.
    ///
    /// Six individuals reach it, with a covariate that is the trait
    /// itself: it separates the individuals that have the condition from
    /// the ones that have not, there is no finite effect of it for a fit
    /// to reach, and the plain null walks towards an infinite one. It is
    /// the fixture of `a binomial null model that walks towards an
    /// infinite coefficient` of
    /// `tests/reference/gwas/refusals_of_both_layers.json`, with the
    /// kinship a user brings to mean no relatedness added to it.
    #[test]
    fn a_logistic_null_that_runs_away_names_the_mixed_model() {
        let phenotype = [0.0_f64, 1.0, 0.0, 1.0, 0.0, 1.0];
        let values: Vec<f64> = phenotype.iter().flat_map(|value| [1.0, *value]).collect();
        let num_individuals = phenotype.len();
        let mut kinship = vec![0.0_f64; num_individuals * num_individuals];
        for (at, row) in kinship.chunks_exact_mut(num_individuals).enumerate() {
            row[at] = 1.0;
        }
        let tested: Vec<usize> = (0..num_individuals).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Binomial,
            design: &values,
            num_coefs: 2,
            kinship: Some(&kinship),
            test: Some(TestType::Score),
            use_grammar_gamma_approx: false,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        let design = match Design::of_the_study(&study, num_individuals) {
            Ok(design) => design,
            Err(error) => panic!("the design of the six individuals: {error}"),
        };
        match LogisticMixedModel::of_the_study(&phenotype, &design, &kinship) {
            Err(error @ Error::GwasFitDidNotSettle { .. }) => {
                let message = error.to_string();
                assert!(
                    message.contains("a binomial trait with a kinship is a logistic mixed model"),
                    "the fit of the separated six was refused with {message}"
                );
                assert!(
                    message.contains("it is the kinship to look at"),
                    "the fit of the separated six was refused with {message}"
                );
            }
            Err(error) => panic!("the fit of the separated six gave {error}"),
            Ok(model) => panic!(
                "the fit of the separated six gave a variance of {}",
                model.genetic_variance
            ),
        }
    }

    /// A kinship that is the identity, which is what a user passes to mean
    /// no relatedness, still leaves this model with a variance of its own,
    /// where the linear mixed model's restricted likelihood goes flat and
    /// gives none.
    ///
    /// **Open 3** of `docs/specs/gwas.md` was measured on that model,
    /// where the two variances are those of the kinship and of the
    /// identity and an identity kinship makes them the same thing. This
    /// model has one variance and the other half of its covariance is the
    /// reciprocals of the weights, which differ from individual to
    /// individual, so the identity is not that other half and the two are
    /// told apart. Measured on 24 September 2026 on the panel with every
    /// genotype called: the fit lands at 0.1274 in 6 steps on both
    /// backends. What is asserted is what the spec states, that the fit
    /// gives a variance, and not that number, which no program outside
    /// popnei has produced.
    #[test]
    fn an_identity_kinship_still_leaves_this_model_a_variance_of_its_own() {
        let (phenotype, values, _) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let num_individuals = phenotype.len();
        let mut kinship = vec![0.0_f64; num_individuals * num_individuals];
        for (at, row) in kinship.chunks_exact_mut(num_individuals).enumerate() {
            row[at] = 1.0;
        }
        let (design, _) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
        let model = match LogisticMixedModel::of_the_study(&phenotype, &design, &kinship) {
            Ok(model) => model,
            Err(error) => panic!("the fit over an identity kinship: {error}"),
        };
        let found = model.null_model(TestType::Score).genetic_variance;
        match found {
            Some(found) => assert!(
                found > 0.0 && found.is_finite(),
                "the variance of the kinship effect over an identity kinship is {found}"
            ),
            None => panic!("the fit over an identity kinship gave no variance"),
        }
    }

    /// A kinship of another size than the individuals that are tested is
    /// refused before any round is run, which is what a caller of the core
    /// crate can bring: `calc_gwas` checks the matrix before any model is
    /// fitted.
    #[test]
    fn a_kinship_that_is_not_of_the_tested_individuals_is_refused() {
        let (phenotype, values, kinship) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let (design, null) = the_design_and_the_null_of(&phenotype, &values, &kinship, &tested);
        let shorter = &kinship[..kinship.len() - 200];
        match TheLinearization::of_the_logistic_null(&phenotype, &design, shorter, &null) {
            Err(Error::GwasInputOfAnotherSize { .. }) => {}
            Ok(_) => panic!("a kinship of 39800 values was accepted"),
            Err(error) => panic!("a kinship of 39800 values gave {error}"),
        }
    }

    /// How far `1 / se²` of the score test may be from GMMAT's `VAR`, as a
    /// share of it: 1e-5.
    ///
    /// `VAR` is the variance of the score, which is `x' p x`, the
    /// denominator of the test. It is the bound "How it is verified" of
    /// "The logistic mixed model" of `docs/specs/gwas.md` sets, on the six
    /// literals and on the whole column in Python alike.
    /// `tests/reference/gwas/gmmat.panel_called.glmm.score.tsv` and the file
    /// of the other panel are printed to six significant digits, which
    /// rounds a value by up to 5e-6 of itself, so half of this bound can go
    /// on GMMAT's printing alone and the comparison has twofold headroom at
    /// best.
    ///
    /// Measured over the six variants of both panels on 24 September 2026:
    /// the worst is `var0751` of the panel with every genotype called,
    /// 2.6832488e-6 of `VAR` on Accelerate and 2.6832488e-6 on faer, which
    /// is 27 per cent of what is allowed. The whole of that can be GMMAT's
    /// printing: it writes that `VAR` as 10.563, so half of its last digit
    /// is 4.73e-6 of the value. The bound is the spec's and is not lowered
    /// to two or three times what was measured, as the bounds this module
    /// sets on popnei's own arithmetic are: what it measures is how far
    /// popnei's fit and GMMAT's land apart. What popnei's own arithmetic is
    /// worth here is the distance between the two backends, 3.4e-16 of
    /// `VAR` at that same variant, ten orders below the distance from
    /// GMMAT.
    const OF_GMMAT_SCORE_VARIANCE: f64 = 1e-5;

    /// How far a p-value of the score test may be from GMMAT's, in `log10`:
    /// 1e-4, from the same item of the same spec.
    ///
    /// The p-values of a study span orders of magnitude and what a user
    /// reads is the exponent, so they are compared in `log10`, which is
    /// already a scale. `log10` shrinks a relative difference, so this bound
    /// has more headroom over GMMAT's six printed digits than
    /// [`OF_GMMAT_SCORE_VARIANCE`] has.
    ///
    /// Measured over the six variants of both panels on 24 September 2026:
    /// the worst is `var0052` of the panel with every genotype called,
    /// 5.2446e-6 in `log10` on both backends, which is 5 per cent of what
    /// is allowed. That one is not the printing alone: 5.2446e-6 in `log10`
    /// is 1.21e-5 of the p-value, against the 5e-6 six digits round it by.
    const OF_GMMAT_SCORE_P_VALUE: f64 = 1e-4;

    /// What GMMAT 1.5.0's `glmm.score` gave for six variants of each panel,
    /// from `tests/reference/gwas/gmmat.panel_called.glmm.score.tsv` and
    /// `gmmat.panel.glmm.score.tsv`: the id, then `VAR` and `PVAL` with
    /// every genotype called and the two of the panel with 3 in 100
    /// genotypes missing whole.
    ///
    /// `VAR` is the variance of the score, `x' p x`, which is `1 / se²`.
    /// GMMAT was given both covariates and, for both panels, the kinship of
    /// the panel with every genotype called. A missing genotype takes the
    /// mean dosage of its variant, which GMMAT calls `impute2mean` and which
    /// is popnei's rule too, so the two agree on the second panel. Five of
    /// the six are the causal variants of `causal_vars.csv` and `var0000` is
    /// not causal.
    const OF_GMMAT_SIX: [(&str, f64, f64, f64, f64); 6] = [
        ("var0000", 6.486_64, 0.702_659, 6.430_05, 0.685_719),
        ("var0052", 8.988_34, 0.030_670_3, 8.737_18, 0.029_175_9),
        ("var0629", 6.499_56, 0.089_510_4, 6.193_48, 0.123_19),
        ("var0751", 10.563, 0.014_233_1, 10.153_2, 0.026_766),
        ("var1137", 9.098, 0.093_808, 8.923_84, 0.108_032),
        ("var1188", 9.050_95, 0.026_273_7, 8.706_3, 0.019_680_6),
    ];

    /// The study of the variants of one panel against its binomial trait
    /// `binom`, with both covariates and the kinship plink2 wrote for the
    /// panel with every genotype called.
    ///
    /// The kinship is `panel_called`'s for both panels, because that is what
    /// `tests/reference/gwas/make_reference.py` gave GMMAT: it fits one null
    /// model with it and then scores the variants of each panel against that
    /// fit.
    fn the_study_of_the_panel(name: &str) -> Gwas {
        let path = the_vcf_of_the_panel(name);
        let options = VcfOptions {
            ploidy: 2,
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let individuals = reader.individuals().to_vec();
        assert_eq!(
            individuals,
            the_individuals_of_the_kinship("panel_called"),
            "the individuals of {name} against the ones of the kinship"
        );
        let kinship = the_kinship_of("panel_called");
        let (phenotype, design) =
            the_trait_and_the_design_of_the_panel(&individuals, TraitType::Binomial);
        let tested: Vec<usize> = (0..individuals.len()).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Binomial,
            design: &design,
            num_coefs: 3,
            kinship: Some(&kinship),
            test: None,
            use_grammar_gamma_approx: false,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of {name}: {error}"),
        }
    }

    /// The factor pyNei's `estimate_gamma` gives for each panel under the
    /// logistic mixed model, with the smallest, the largest and the
    /// standard deviation of the 100 ratios it is the mean of.
    ///
    /// Measured on 24 September 2026 with pyNei at the commit
    /// `pyproject.toml` names, numpy 2.5.3 on Accelerate, on the trait
    /// `binom` and the covariates `cov1` and `cov2` of
    /// `tests/reference/gwas/phenotypes.csv` over the kinship plink2 wrote
    /// for the panel with every genotype called.
    ///
    /// The factor is five times smaller than the linear mixed model's on
    /// the same panel, 0.105 against 0.517, because the projection matrix
    /// of this model is weighted by the variance of a binomial trait, which
    /// is at most a quarter. The spread around it is the same shape: the
    /// ratios run from 0.0726 to 0.1220, a standard deviation of 0.00966,
    /// which is 9.2 per cent of the mean and a largest over smallest of
    /// 1.68, against 12.9 per cent and 2.15 for the linear mixed model. So
    /// one factor stands in for a quantity that differs from variant to
    /// variant here as it does there, and a factor that fell outside that
    /// range would be reporting the ratios of some other set of variants.
    const OF_PYNEI_GAMMA_OF_EACH_PANEL: [(&str, f64, f64, f64, f64); 2] = [
        (
            "panel_called",
            0.105_320_734_389_783_95,
            0.072_578_464_929_525_54,
            0.121_993_904_517_353_51,
            0.009_655_522_042_929_226,
        ),
        (
            "panel",
            0.106_058_860_871_757_5,
            0.073_708_044_241_363_11,
            0.123_332_208_807_525_8,
            0.009_666_946_916_712_202,
        ),
    ];

    /// How far the factor of a panel may be from the one pyNei's
    /// `estimate_gamma` gives for it: 1e-8 of it, relative.
    ///
    /// It measures the distance between two fits and the two projection
    /// matrices they build, the factor being read out of pyNei at full
    /// precision above. This model's fit is the one popnei writes
    /// differently from pyNei, a Cholesky and a solve where pyNei inverts,
    /// so what the bound covers here is that route as well as the
    /// arithmetic.
    ///
    /// Measured over the two panels on 24 September 2026: the worst is the
    /// panel with every genotype called, 1.77e-14 of the factor away on
    /// Accelerate and 1.84e-14 on faer, and the panel with 3 in 100
    /// genotypes missing is 1.78e-14 and 1.82e-14. The bound is 5.4 times
    /// the worst of the four.
    ///
    /// It is 100000 times tighter than the 1e-8 the linear mixed model's
    /// own `OF_PYNEI_GAMMA` takes, and the difference is the fit and not
    /// this model's arithmetic: pyNei and popnei walk the search on the
    /// variance of the kinship effect to the same value here, where the
    /// linear mixed model's restricted maximum likelihood has a criterion
    /// that is flat at its minimum and the two land on eigenvalues that
    /// differ in their last bits.
    const OF_PYNEI_GAMMA: f64 = 1e-13;

    /// The factor of each panel under the logistic mixed model is the one
    /// pyNei's `estimate_gamma` gives, estimated from 100 variants of the
    /// first block.
    ///
    /// The first block of either panel holds all 1200 of its variants and
    /// all 1200 vary, so the 100 the factor comes from are the first 100 of
    /// the file and are the ones pyNei took. It is the one check there is
    /// on the factor itself: no program outside the project computes it,
    /// GMMAT making the exact denominator, and the relation to popnei's own
    /// exact answer allows a p-value to be out by a factor of 30.
    #[test]
    fn the_factor_of_each_panel_is_pyneis_over_a_hundred_variants() {
        for (name, gamma, smallest, largest, deviation) in OF_PYNEI_GAMMA_OF_EACH_PANEL {
            let (mut fitted, dosages) = the_null_and_the_first_block_of(name);
            assert_eq!(dosages.num_vars(), 1200, "the first block of {name}");
            assert_eq!(
                dosages.num_with_variance(),
                1200,
                "the variants of that block that vary"
            );
            assert_eq!(fitted.grammar_gamma_factor(), None, "before the estimate");

            if let Err(error) = fitted.approximate_the_denominator(&dosages) {
                panic!("the factor of {name}: {error}");
            }

            let found = fitted
                .grammar_gamma_factor()
                .expect("the factor of a model that approximates");
            let away = (found - gamma).abs() / gamma;
            assert!(
                away <= OF_PYNEI_GAMMA,
                "the factor of {name} is {found} and pyNei gives {gamma}, {away} of it \
                 away, against the {OF_PYNEI_GAMMA} allowed"
            );
            assert!(
                found > smallest && found < largest,
                "the factor of {name} is {found} and the 100 ratios it is the mean of \
                 run from {smallest} to {largest}, a standard deviation of {deviation}"
            );
        }
    }

    /// The logistic mixed model of a panel, fitted as
    /// [`the_study_of_the_panel`] fits it, with the dosages of the first
    /// block of a pass over its variants.
    ///
    /// The pair is what the second pass of the GRAMMAR-Gamma approximation
    /// gives `calc_gwas`: the first block of a pass that asks for the
    /// genotypes alone, with `Reblock` in front of it.
    fn the_null_and_the_first_block_of(name: &str) -> (LogisticMixedModel, GwasDosages) {
        let path = the_vcf_of_the_panel(name);
        let options = VcfOptions {
            ploidy: 2,
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let individuals = reader.individuals().to_vec();
        let kinship = the_kinship_of("panel_called");
        let (phenotype, values) =
            the_trait_and_the_design_of_the_panel(&individuals, TraitType::Binomial);
        let tested: Vec<usize> = (0..individuals.len()).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Binomial,
            design: &values,
            num_coefs: 3,
            kinship: Some(&kinship),
            test: None,
            use_grammar_gamma_approx: false,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        let design = match Design::of_the_study(&study, individuals.len()) {
            Ok(design) => design,
            Err(error) => panic!("the design of {name}: {error}"),
        };
        let fitted = match LogisticMixedModel::of_the_study(&phenotype, &design, &kinship) {
            Ok(fitted) => fitted,
            Err(error) => panic!("the null model of {name}: {error}"),
        };
        reader.set_needs(Needs::GTS);
        let mut blocks = match Reblock::new(&mut reader, None) {
            Ok(blocks) => blocks,
            Err(error) => panic!("the blocks of {name}: {error}"),
        };
        let mut block = match blocks.next_block() {
            Ok(Some(block)) => block,
            Ok(None) => panic!("{name} gave no block"),
            Err(error) => panic!("the first block of {name}: {error}"),
        };
        let mut dosages = GwasDosages::of_a_study();
        if let Err(error) = dosages.read_the_block(
            &mut block,
            &design,
            BlockOfThePass {
                ploidy: 2,
                first_var: 0,
            },
        ) {
            panic!("the dosages of the first block of {name}: {error}");
        }
        (fitted, dosages)
    }

    /// How far the middle of the absolute `log10(p_approx / p_exact)` may
    /// be from 0 over a whole panel, how far the worst variant of it may
    /// be, and how far an effect may be from the exact one as a share of
    /// it: 0.1, 1.5 and 0.5.
    ///
    /// The statistic the first of them bounds is the value at index 600 of
    /// the 1200 absolute log ratios of a panel sorted, which is the middle
    /// one of an even count read at the upper of the two, and it is the
    /// reading the linear mixed model's own `OF_THE_APPROXIMATION` takes
    /// and explains: the spec's number is the median of the **signed** log
    /// ratios, which the pytest and node suites take, and taking the
    /// absolute values first is the stricter of the two.
    ///
    /// They are the three numbers of "What it gives" and "How it is
    /// verified" of the approximation in `docs/specs/gwas.md`, which has
    /// them from `test_grammar_gamma_approx` of pyNei and states them for
    /// the linear mixed model. The spec does not measure the logistic mixed
    /// model, and this test holds it to the same three: what the bound is
    /// about is one factor standing in for a quantity that differs from
    /// variant to variant, which is the same in both models.
    ///
    /// Measured over the 1200 variants of each panel on 24 September 2026,
    /// the same on Accelerate and on faer to six digits: on the panel with
    /// every genotype called the median is 0.00981, the largest 0.524 and
    /// the worst effect 0.393; on the panel with 3 in 100 genotypes missing
    /// they are 0.0107, 0.525 and 0.385. So the worst variant uses 35 per
    /// cent of what its p-value is allowed and 79 per cent of what its
    /// effect is, which is what the linear mixed model's own
    /// `OF_THE_APPROXIMATE_EFFECT` explains: the effect is `num / den` and
    /// the approximation leaves `num` alone, so the share an effect moves
    /// by is exactly the share that variant's own ratio of the two
    /// denominators is from the factor.
    const OF_THE_APPROXIMATION: (f64, f64, f64) = (0.1, 1.5, 0.5);

    /// How far the worst effect of a panel has to move under the
    /// approximation for the approximation to have happened: 0.05 of it.
    ///
    /// The three bounds above are ceilings, and a run that made no
    /// approximation at all would pass all of them, an exact answer being
    /// at no distance from itself. It is not a hypothesis: a reviewer
    /// replaced the approximation with `None` in this model on 24 September
    /// 2026 and the whole core suite stayed green, the pytest and node
    /// suites with it. So the effect is held above a floor as well.
    ///
    /// The floor is 0.05 where the worst effect measured that day is 0.3933
    /// on the panel with every genotype called and 0.385 on the panel with
    /// genotypes missing, a factor of eight of room, and the p-values are
    /// left alone: they are where the approximation moves least, the median
    /// of the absolute log ratios being 0.00981 here, and a floor near that
    /// would go red on a panel the approximation happens to suit.
    const THE_APPROXIMATION_MOVES_THE_EFFECT: f64 = 0.05;

    /// The approximation gives every variant of a panel an answer near the
    /// exact one, and the result says that it was used.
    ///
    /// There is no program outside the project to check this against: GMMAT
    /// makes the exact denominator. So what is compared is the same study
    /// of the same panel with and without the approximation, which is
    /// pyNei's own `test_grammar_gamma_approx` one model along. The bound
    /// is loose, and the factor itself is where an error that does not grow
    /// with the panel would show, which the test against pyNei's
    /// `estimate_gamma` above is for.
    #[test]
    fn the_approximate_answers_of_both_panels_are_near_the_exact_ones() {
        let (of_the_median, of_the_largest, of_the_effects) = OF_THE_APPROXIMATION;
        for name in ["panel_called", "panel"] {
            let exact = the_study_of_the_panel(name);
            let approximated = the_approximated_study_of_the_panel(name);

            assert!(
                approximated.used_grammar_gamma_approx,
                "the study of {name} says it approximated"
            );
            assert!(
                !exact.used_grammar_gamma_approx,
                "the study of {name} that did not approximate says so"
            );
            assert_eq!(approximated.num_vars, 1200, "the variants of {name}");
            assert_eq!(
                approximated.null_model, exact.null_model,
                "the null model of {name} is the same fit either way, since the \
                 approximation is of the test of a variant and not of the fit"
            );
            let mut of_the_p_values: Vec<f64> = Vec::new();
            let mut worst_effect = 0.0_f64;
            for (var, ((p_approx, p_exact), (beta_approx, beta_exact))) in approximated
                .p_value
                .iter()
                .zip(&exact.p_value)
                .zip(approximated.beta.iter().zip(&exact.beta))
                .enumerate()
            {
                assert!(
                    p_approx.is_finite() && p_exact.is_finite(),
                    "the variant {var} of {name} has a p-value under both, {p_approx} \
                     approximated and {p_exact} exact"
                );
                of_the_p_values.push((p_approx / p_exact).log10().abs());
                worst_effect =
                    worst_effect.max((beta_approx - beta_exact).abs() / beta_exact.abs());
            }
            of_the_p_values.sort_by(f64::total_cmp);
            let median = of_the_p_values.get(600).copied().unwrap_or(f64::NAN);
            let largest = of_the_p_values.last().copied().unwrap_or(f64::NAN);
            assert!(
                median <= of_the_median,
                "the middle variant of {name} moves its p-value by {median} in the log, \
                 against the {of_the_median} allowed"
            );
            assert!(
                largest <= of_the_largest,
                "the worst variant of {name} moves its p-value by {largest} in the log, \
                 against the {of_the_largest} allowed"
            );
            assert!(
                worst_effect <= of_the_effects,
                "the worst effect of {name} moves by {worst_effect} of itself, against \
                 the {of_the_effects} allowed"
            );
            assert!(
                worst_effect >= THE_APPROXIMATION_MOVES_THE_EFFECT,
                "the worst effect of {name} moves by {worst_effect} of itself, and an \
                 approximation that happened moves it by \
                 {THE_APPROXIMATION_MOVES_THE_EFFECT} at least"
            );
        }
    }

    /// The study of a panel with the GRAMMAR-Gamma approximation, which is
    /// the study [`the_study_of_the_panel`] makes with a second pass over
    /// the same file beside it.
    fn the_approximated_study_of_the_panel(name: &str) -> Gwas {
        let path = the_vcf_of_the_panel(name);
        let options = VcfOptions {
            ploidy: 2,
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let mut gamma_pass = match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let individuals = reader.individuals().to_vec();
        let kinship = the_kinship_of("panel_called");
        let (phenotype, design) =
            the_trait_and_the_design_of_the_panel(&individuals, TraitType::Binomial);
        let tested: Vec<usize> = (0..individuals.len()).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Binomial,
            design: &design,
            num_coefs: 3,
            kinship: Some(&kinship),
            test: None,
            use_grammar_gamma_approx: true,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        match crate::gwas::calc_gwas(&mut reader, Some(&mut gamma_pass), &study) {
            Ok(result) => result,
            Err(error) => panic!("the approximated study of {name}: {error}"),
        }
    }

    /// The six variants of both panels are GMMAT's under the score test,
    /// which is deliverable 2 of work package 2 of
    /// `docs/plans/gwas-logistic.md`.
    ///
    /// The score test asks how steeply the likelihood rises at an effect of
    /// 0, so it holds the variance of the kinship effect and the fitted
    /// chances at the null and needs no fit with the variant in, which is
    /// the only test this model has. What is compared is `1 / se²` against
    /// the variance of the score GMMAT reports, which is the denominator the
    /// test is built on, and the p-value in `log10`.
    ///
    /// The test the study made is asserted with them, and no test was asked
    /// for: the default of a binomial trait with a kinship is the score
    /// test, since a Wald one would fit one mixed model per variant.
    ///
    /// The second panel is what says that a missing genotype takes the mean
    /// dosage of its variant as GMMAT's `impute2mean` does: 3 in 100 of its
    /// genotypes are missing whole, and every one of the six numbers moves
    /// between the panels.
    #[test]
    fn the_six_variants_of_both_panels_are_gmmats_variance_of_the_score_and_p_value() {
        for (name, of_the_panel) in [("panel_called", 0_usize), ("panel", 1)] {
            let result = the_study_of_the_panel(name);
            assert_eq!(result.num_vars, 1200, "the variants of {name}");
            assert_eq!(
                result.null_model.model,
                GwasModel::Glmm,
                "the model the result of {name} says it fitted"
            );
            assert_eq!(
                result.null_model.test,
                TestType::Score,
                "the test the result of {name} says it made"
            );
            for (id, called_variance, called_p, missing_variance, missing_p) in OF_GMMAT_SIX {
                let (variance, p_value) = match of_the_panel {
                    0 => (called_variance, called_p),
                    _ => (missing_variance, missing_p),
                };
                let var = the_row_of(&result, id);
                let se = result.se[var];
                let found = 1.0 / (se * se);
                let difference = (found - variance).abs();
                assert!(
                    difference <= OF_GMMAT_SCORE_VARIANCE * variance,
                    "1 / se² of {id} of {name} is {found} and GMMAT gives {variance}, \
                     {difference} away, which is {share} of it against the \
                     {OF_GMMAT_SCORE_VARIANCE} allowed",
                    share = difference / variance
                );
                let found = result.p_value[var];
                let difference = (found / p_value).log10().abs();
                assert!(
                    difference <= OF_GMMAT_SCORE_P_VALUE,
                    "the p-value of {id} of {name} is {found} and GMMAT gives {p_value}, \
                     {difference} apart in log10 against the {OF_GMMAT_SCORE_P_VALUE} \
                     allowed"
                );
            }
        }
    }

    /// The kinship of the fixture of **Open 2** below: eight individuals in
    /// two families of four, who are related within a family and not
    /// between them.
    ///
    /// It is the matrix the linear mixed model's own fixture of that rule
    /// takes, so the two score tests are read over the same relatedness.
    const THE_KINSHIP_OF_TWO_FAMILIES: [f64; 64] = [
        1.0, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, //
        0.2, 1.0, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, //
        0.2, 0.2, 1.0, 0.2, 0.0, 0.0, 0.0, 0.0, //
        0.2, 0.2, 0.2, 1.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 0.0, 0.0, 0.0, 1.0, 0.2, 0.2, 0.2, //
        0.0, 0.0, 0.0, 0.0, 0.2, 1.0, 0.2, 0.2, //
        0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 1.0, 0.2, //
        0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 1.0,
    ];

    /// The design of that fixture: the intercept and one covariate, which is
    /// the dosages of the first variant of the fixture, `1 0 2 0 1 2 0 1`,
    /// in units a tenth of theirs, row after row.
    ///
    /// It is the design the plain logistic model's fixture of the same rule
    /// takes, which is what a user gets by putting a genotype in as a
    /// covariate: the covariate carries the same information as the variant
    /// whatever it is multiplied by, so the projection matrix takes the
    /// whole of that variant out and what is left of `x' p x` is rounding.
    const THE_DESIGN_OF_THE_FIRST_VARIANT: [f64; 16] = [
        1.0, 0.1, //
        1.0, 0.0, //
        1.0, 0.2, //
        1.0, 0.0, //
        1.0, 0.1, //
        1.0, 0.2, //
        1.0, 0.0, //
        1.0, 0.1,
    ];

    /// The binomial trait of those eight individuals: five have the
    /// condition and three have not, and neither the intercept alone nor the
    /// covariate separates the two groups, so the fit settles.
    const THE_TRAIT_OF_EIGHT: [f64; 8] = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0];

    /// The same eight individuals with four of them having the condition
    /// and the other four not, which leaves the intercept of the plain
    /// logistic null at 0.36 and the largest absolute value of the linear
    /// predictor of a linearization at 0.62, against the 4.19 of the panel.
    ///
    /// It is the fixture of
    /// [`the_rounds_of_a_linearization_are_measured_against_the_predictor_plus_one`],
    /// and what it is for is that the plus 1 of the stopping rule is a
    /// quarter of the divisor on the panel and more than half of it here.
    const THE_BALANCED_TRAIT_OF_EIGHT: [f64; 8] = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0];

    /// All eight of them are tested.
    const THE_INDIVIDUALS_OF_EIGHT: [usize; 8] = [0, 1, 2, 3, 4, 5, 6, 7];

    /// The design of the eight individuals of the two families, checked,
    /// and the plain logistic null of `phenotype` over it, for a fixture
    /// that holds `kinship`.
    fn the_design_and_the_null_of_eight(kinship: &[f64]) -> (Design<'_>, LogisticModel) {
        let study = GwasInput {
            phenotype: &THE_TRAIT_OF_EIGHT,
            trait_type: TraitType::Binomial,
            design: &THE_DESIGN_OF_THE_FIRST_VARIANT,
            num_coefs: 2,
            kinship: Some(kinship),
            test: Some(TestType::Score),
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS_OF_EIGHT,
            transform_to_biallelic: false,
        };
        let design = match Design::of_the_study(&study, THE_INDIVIDUALS_OF_EIGHT.len()) {
            Ok(design) => design,
            Err(error) => panic!("the design of the eight individuals: {error}"),
        };
        let null = match LogisticModel::of_the_study(&THE_TRAIT_OF_EIGHT, &design, TestType::Score)
        {
            Ok(fitted) => fitted,
            Err(error) => panic!("the logistic null of the eight individuals: {error}"),
        };
        (design, null)
    }

    /// The rounds of a linearization are counted against the largest
    /// absolute value of the linear predictor **plus 1**, so a fit whose
    /// predictor is near 0 is asked for a change below 1e-6 in absolute
    /// terms and not for one that is small against nothing.
    ///
    /// The divisor is the predictor the round started at, which is
    /// `_fit_pql_for_tau` of `pynei/gwas.py`, and the plus 1 is pyNei's
    /// too. Nothing pinned it: the panel's predictor reaches 4.19, where
    /// the plus 1 moves the divisor by a quarter and the whole fit is
    /// bit-identical without it on both backends.
    ///
    /// The fixture is the eight individuals of the two families with four
    /// of them having the condition and four not, whose largest predictor
    /// is 0.62 at a variance of 0.5, so the plus 1 moves the divisor by a
    /// factor of 2.6. Measured on 24 September 2026 on both backends: with
    /// the plus 1 the linearization at that variance settles in 3 rounds
    /// and the whole fit takes 35 linearizations over 15 steps; without it,
    /// 4 rounds, 38 linearizations and 16 steps.
    ///
    /// The counts are what this reads and not the numbers the fit lands
    /// on: the variance of the kinship effect is 1.2245216262670948
    /// under the rule and 1.2245196763662805 without the plus 1, 1.9e-6
    /// apart and 1.6e-6 of its own value. No reference program was run on
    /// this fixture.
    #[test]
    fn the_rounds_of_a_linearization_are_measured_against_the_predictor_plus_one() {
        let study = GwasInput {
            phenotype: &THE_BALANCED_TRAIT_OF_EIGHT,
            trait_type: TraitType::Binomial,
            design: &THE_DESIGN_OF_THE_FIRST_VARIANT,
            num_coefs: 2,
            kinship: Some(&THE_KINSHIP_OF_TWO_FAMILIES),
            test: Some(TestType::Score),
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS_OF_EIGHT,
            transform_to_biallelic: false,
        };
        let design = match Design::of_the_study(&study, THE_INDIVIDUALS_OF_EIGHT.len()) {
            Ok(design) => design,
            Err(error) => panic!("the design of the balanced eight: {error}"),
        };
        let null = match LogisticModel::of_the_study(
            &THE_BALANCED_TRAIT_OF_EIGHT,
            &design,
            TestType::Score,
        ) {
            Ok(fitted) => fitted,
            Err(error) => panic!("the logistic null of the balanced eight: {error}"),
        };
        let mut fitted = match TheLinearization::of_the_logistic_null(
            &THE_BALANCED_TRAIT_OF_EIGHT,
            &design,
            &THE_KINSHIP_OF_TWO_FAMILIES,
            &null,
        ) {
            Ok(fitted) => fitted,
            Err(error) => panic!("the linearization of the balanced eight: {error}"),
        };
        if let Err(error) = fitted.at(0.5) {
            panic!("the linearization of the balanced eight at a variance of 0.5: {error}");
        }
        let largest = fitted
            .linear_predictor
            .iter()
            .fold(0.0_f64, |largest, value| largest.max(value.abs()));
        assert!(
            (0.5..1.0).contains(&largest),
            "the largest absolute value of the linear predictor is {largest}, where the plus \
             1 is what this reads"
        );
        assert_eq!(
            fitted.rounds, 3,
            "the rounds the linearization at a variance of 0.5 took"
        );
        let model = match LogisticMixedModel::of_the_study(
            &THE_BALANCED_TRAIT_OF_EIGHT,
            &design,
            &THE_KINSHIP_OF_TWO_FAMILIES,
        ) {
            Ok(model) => model,
            Err(error) => panic!("the fit of the balanced eight: {error}"),
        };
        assert_eq!(
            model.linearizations, 35,
            "the times the fit of the balanced eight factored the covariance"
        );
        assert_eq!(
            model.steps_on_the_variance, 15,
            "the steps the fit of the balanced eight took on the variance"
        );
    }

    /// A variant that the projection leaves nothing of has no answer under
    /// this model's score test either, which is the meanwhile of **Open 2**
    /// of `docs/specs/gwas.md` reaching the last of its four places.
    ///
    /// The covariate is the first variant's dosages in units a tenth of
    /// theirs, and the projection matrix takes the design out of whatever it
    /// is applied to, so `x' p x` is 0 in exact arithmetic and what is left
    /// is rounding, of whichever sign it fell on. `beta` would be a number
    /// divided by noise and the row would read as a variant that was tested
    /// and showed nothing.
    ///
    /// The threshold is the tested individuals times 2.2e-16 of the
    /// variant's own squared length times the largest value of the diagonal
    /// of the projection matrix, which is the scale `docs/plans/gwas-linear.md`
    /// chose for the linear mixed model's score test and which this model
    /// takes so that the two answer alike. Measured on this fixture on 24
    /// September 2026, which is the first time either mixed model's score
    /// test has reached that rule at all: the collinear variant keeps
    /// 2.082e-17 of `x' p x` on Accelerate and 2.533e-16 on faer, against a
    /// threshold of 2.780e-15, which is the eight tested individuals times
    /// 2.2e-16 of the variant's squared length of 11 times the 0.14228 the
    /// largest diagonal entry of the projection matrix is. The ordinary
    /// variant beside it keeps 0.82216, against a threshold of 4.044e-15.
    /// So the threshold sits 11 times above the largest rounding the two
    /// backends left and fourteen orders of magnitude below a variant that
    /// has something to test. What is left of the collinear variant is
    /// rounding and nothing else, so it moves with any change to the
    /// arithmetic of the fit and these two numbers are of the fit as it
    /// stands; what the fixture rests on is the threshold and the 0.82216,
    /// which do not.
    ///
    /// What the study gave with the threshold taken out, measured the same
    /// day: a `beta` of -14384 with an `se` of 2.684e8 and a p-value of
    /// 0.999957 on Accelerate, and a `beta` of -842.4 with an `se` of
    /// 6.511e7 and a p-value of 0.999990 on faer. Neither is a NaN, so
    /// nothing marks the row: a user reads a variant that was tested and
    /// showed nothing, and the two builds differ by a factor of 17 in the
    /// effect they report for it.
    ///
    /// What this fixture pins is that a threshold is there and not the
    /// scale it is built from. The collinear variant keeps 2.776e-17 of
    /// `x' p x` against a threshold of 2.780e-15 and the ordinary one
    /// keeps 0.82216, so the largest diagonal entry of the projection
    /// matrix, 0.14228, can be dropped from the product and both variants
    /// are still answered the way they are: the threshold would be
    /// 1.954e-14, still 77 times above the rounding and thirteen orders
    /// below the ordinary variant. What says that this model takes the
    /// same scale as the linear mixed model's score test, so that the two
    /// answer alike, is the doc comment of
    /// [`LogisticMixedModel::largest_of_the_projection`] and no fixture
    /// here; a fixture that read it would need a variant whose
    /// `x' p x` fell between the two thresholds, which is a band of a
    /// factor of 7.
    ///
    /// The ordinary variant beside it is asserted to have an answer and not
    /// to any number: no reference program was run on this fixture, and the
    /// six literals of both panels are what say the numbers are right.
    #[test]
    fn a_variant_the_projection_leaves_nothing_of_has_no_answer() {
        let mut vcf = String::from(THE_HEADER_OF_EIGHT);
        for (var, genotypes) in [
            // the dosages of the covariate, allele by allele
            ["0/1", "0/0", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1"],
            ["0/0", "0/0", "0/0", "0/0", "1/1", "1/1", "1/1", "1/1"],
        ]
        .iter()
        .enumerate()
        {
            let pos = var.saturating_add(1).saturating_mul(1000);
            vcf.push_str(&format!("1\t{pos}\tv{var}\tA\tT\t.\t.\t.\tGT"));
            for genotype in genotypes {
                vcf.push('\t');
                vcf.push_str(genotype);
            }
            vcf.push('\n');
        }
        let study = GwasInput {
            phenotype: &THE_TRAIT_OF_EIGHT,
            trait_type: TraitType::Binomial,
            design: &THE_DESIGN_OF_THE_FIRST_VARIANT,
            num_coefs: 2,
            kinship: Some(&THE_KINSHIP_OF_TWO_FAMILIES),
            test: None,
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS_OF_EIGHT,
            transform_to_biallelic: false,
        };
        let mut reader = reader_over(vcf.as_bytes());
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of a variant that is a covariate: {error}"),
        };

        assert_eq!(result.num_vars, 2, "the variants of the fixture");
        // The variant is in the result with its frequency, as every variant
        // that has no answer is: seven of its sixteen alleles are the one
        // that is not the major one.
        let found = result.allele_freq[0];
        assert!(
            (found - 0.4375).abs() <= 1e-12,
            "the frequency of v0 is {found} and seven of its sixteen alleles are the one \
             that is not the major one"
        );
        assert!(
            result.beta[0].is_nan() && result.se[0].is_nan() && result.p_value[0].is_nan(),
            "the variant the projection leaves nothing of was answered with a beta of \
             {beta}, an se of {se} and a p-value of {p_value}",
            beta = result.beta[0],
            se = result.se[0],
            p_value = result.p_value[0]
        );
        assert!(
            result.beta[1].is_finite()
                && result.se[1] > 0.0
                && (0.0..=1.0).contains(&result.p_value[1]),
            "the ordinary variant beside it was answered with a beta of {beta}, an se of \
             {se} and a p-value of {p_value}",
            beta = result.beta[1],
            se = result.se[1],
            p_value = result.p_value[1]
        );
    }
}
