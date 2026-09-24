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

use super::logistic::{
    LogisticModel, TheSystemOfTheFit, the_chance_of, the_system_of, the_system_that_is_left,
};
use super::result::NullModel;
use super::study::{Design, GwasInputShape, GwasModel, TestType};

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
            true => Err(Error::GwasFitDidNotSettle {
                model: GwasModel::Glmm,
                rounds: self.rounds,
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
        the_system_of(
            popnei_linalg::solve_with_cholesky(&self.dsd, num_coefs, &mut self.coefs, 1),
            "solve of the design weighted by the covariance of the working trait",
        )?;
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
    /// inside, and a step that would take the variance to 0 or below
    /// quarters it instead, which is a smaller move in the same direction.
    fn the_variance_after(&self, genetic_variance: f64, step: f64) -> f64 {
        let next = genetic_variance + step;
        match (self.too_small, self.too_large) {
            (Some(too_small), Some(too_large)) => match too_small < next && next < too_large {
                true => next,
                false => (too_small * too_large).sqrt(),
            },
            (None, _) | (_, None) => match next > 0.0 {
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
    /// How many individuals the study tests.
    num_individuals: usize,
    /// How many rounds the fit ran in all, which is how many times it
    /// factored the covariance of the working trait: 22 on the panel.
    linearizations: usize,
    /// How many steps it took on the variance of the kinship effect: 8 on
    /// the panel.
    steps_on_the_variance: usize,
    /// How many times it formed the inverse of a factorized matrix of that
    /// size: 1, at the end, for the projection matrix.
    inverses: usize,
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
        let null = LogisticModel::of_the_study(phenotype, design, TestType::Score)?;
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
            fitted.at(genetic_variance)?;
            let (score, information) = step.at(&fitted, genetic_variance)?;
            // A step that is not a finite number is not refused here, and
            // the bracket is what makes that safe: a kinship that says
            // nothing about the trait leaves the derivative and the
            // information both 0 and the step their quotient, and
            // [`TheBracket::the_variance_after`] sends every step that is
            // not above 0, a NaN among them, to a quarter of the variance
            // or to the geometric mean of two ends that are finite. So the
            // variance walks down to the boundary where the kinship
            // explains nothing instead of becoming a number that is not
            // one. Measured on 24 September 2026 on the panel with every
            // genotype called and a kinship of all zeros: 12 steps to a
            // variance of 0.
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
        let projection = the_projection_of(&fitted, of_the_design)?;
        let projected_trait = fitted
            .phenotype
            .iter()
            .zip(&fitted.fitted_chance)
            .map(|(measured, chance)| measured - chance)
            .collect();
        Ok(LogisticMixedModel {
            coefs: std::mem::take(&mut fitted.coefs),
            genetic_variance,
            projection,
            projected_trait,
            num_individuals,
            linearizations: fitted.factorizations,
            steps_on_the_variance,
            inverses: 1,
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

    /// The projection matrix of the fit, `num_individuals` x
    /// `num_individuals`, row after row, which every variant is tested
    /// through.
    #[must_use]
    pub(crate) fn projection(&self) -> &[f64] {
        &self.projection
    }

    /// The residual every variant's numerator is taken against, one value
    /// per tested individual: the trait less the fitted chance.
    #[must_use]
    pub(crate) fn projected_trait(&self) -> &[f64] {
        &self.projected_trait
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
                        *value = 1.0 / weight.sqrt();
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

/// The projection matrix of a linearization that has settled,
/// `num_individuals` x `num_individuals`, row after row.
///
/// It is `p = sigma⁻¹ - sigma⁻¹ d (d' sigma⁻¹ d)⁻¹ d' sigma⁻¹`, and it is
/// the one place of the fit where the inverse of an individuals by
/// individuals matrix is formed, because the score test of a variant wants
/// the matrix and not a solve against it. `of_the_design` is the buffer
/// the part the design explains is formed in, which the step on the
/// variance is done with.
///
/// The inverse comes back with only its lower half written, the inverse of
/// a symmetric matrix being symmetric, so that half is mirrored into the
/// other before the design is taken out of it.
///
/// # Errors
///
/// [`Error::GwasLinalg`] when the inverse, the solve or the product could
/// not be done.
fn the_projection_of(fitted: &TheLinearization<'_>, of_the_design: Vec<f64>) -> Result<Vec<f64>> {
    let num_individuals = fitted.num_individuals;
    let num_coefs = fitted.num_coefs;
    let mut projection = vec![0.0_f64; fitted.covariance.len()];
    popnei_linalg::invert_with_cholesky(&fitted.covariance, num_individuals, &mut projection)
        .map_err(|source| Error::GwasLinalg {
            operation: "inverse of the covariance of the working trait",
            source,
        })?;
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

    use super::{LogisticMixedModel, TheLinearization, TheStepOnTheVariance, the_projection_of};
    use crate::error::Error;
    use crate::gwas::linear::lm::the_trait_and_the_design_of_the_panel;
    use crate::gwas::linear_mixed::lmm::{the_individuals_of_the_kinship, the_kinship_of};
    use crate::gwas::logistic::LogisticModel;
    use crate::gwas::study::{Design, GwasInput, GwasModel, TestType, TraitType};

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
    /// two backends, 2.7e-15 of the variance, nine orders below the
    /// distance from GMMAT.
    const OF_GMMAT_GENETIC_VARIANCE: f64 = 1e-5;

    /// How far the trace of the projection matrix times the kinship taken
    /// from the identity may be from the two matrices multiplied out, as a
    /// share of the absolute terms of that sum: 6e-14.
    ///
    /// The trace is a sum of 40000 products of both signs that cancel to
    /// about half of what went into them, so what it is measured against
    /// is those products' absolute values and not its own value, which is
    /// what "How it is verified" of `docs/specs/gwas.md` asks a tolerance
    /// to be measured against. Measured over both panels at a variance of
    /// 0 and at GMMAT's on 24 September 2026: the worst is 2.08e-14 on
    /// Accelerate and 2.05e-14 on faer, both on the panel with every
    /// genotype called at GMMAT's variance, so this is 2.9 times where it
    /// breaks.
    const OF_THE_TRACE: f64 = 6e-14;

    /// How large a column of the design against a column of the projection
    /// matrix may be, as a share of the absolute terms of that sum:
    /// 6e-15.
    ///
    /// Measured over the 600 sums of each panel on 24 September 2026: the
    /// worst is 1.94e-15 on Accelerate and 8.97e-16 on faer, on the panel
    /// with every genotype called, so this is 3.1 times where it breaks.
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
    /// worst is 1.170e-13 on Accelerate and 1.177e-13 on faer, both on the
    /// panel with every genotype called, so this is 4.2 times where it
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

    /// The study of one of the two panels: its design, checked, and the
    /// plain logistic null the mixed fit starts from.
    fn the_study_of<'a>(
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
        let (design, null) = the_study_of(&phenotype, &values, &kinship, &tested);
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
        let (design, null) = the_study_of(&phenotype, &values, &kinship, &tested);
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
            let (design, null) = the_study_of(&phenotype, &values, &kinship, &tested);
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
        let (design, null) = the_study_of(&phenotype, &values, &kinship, &tested);
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
                assert!(rounds > 0, "the rounds the refusal names, {rounds}");
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
        let (design, null) = the_study_of(&phenotype, &values, &kinship, &tested);
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
        let (design, _) = the_study_of(&phenotype, &values, &kinship, &tested);
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
    /// same two. They also pin the rule that ends a linearization: the
    /// divisor of its stopping rule is the largest absolute value of the
    /// linear predictor plus 1, and dropping that plus 1 leaves the
    /// variance of the kinship effect where it was to all 18 digits and
    /// moves the three effects by 3.1e-13, which no comparison with GMMAT
    /// can see, while the rounds go from 22 to 23, which this sees.
    /// Measured on both backends on 24 September 2026.
    #[test]
    fn the_fit_forms_one_inverse_and_factors_the_covariance_once_a_round() {
        let (phenotype, values, kinship) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let (design, _) = the_study_of(&phenotype, &values, &kinship, &tested);
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
    /// fit that took it after would give an infinity or a NaN here.
    ///
    /// Each trace is measured against the absolute terms of the sum it is
    /// and not against its own value, which is a sum of 40000 products of
    /// both signs; [`OF_THE_TRACE`] is what they are held to.
    #[test]
    fn the_trace_from_the_identity_is_the_trace_of_the_two_matrices() {
        for name in ["panel_called", "panel"] {
            let (phenotype, values, kinship) = the_panel(name);
            let tested: Vec<usize> = (0..phenotype.len()).collect();
            let num_individuals = phenotype.len();
            let (design, null) = the_study_of(&phenotype, &values, &kinship, &tested);
            for variance in [0.0, OF_GMMAT_VARIANCE] {
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
                let projection = match the_projection_of(&fitted, buffer) {
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
                    difference <= OF_THE_TRACE * scale,
                    "{name} at a variance of {variance}: the trace from the identity is \
                     {of_the_identity} and the two matrices multiplied out give \
                     {of_the_matrices}, {difference} apart, where the absolute terms of that \
                     sum are {scale} and {OF_THE_TRACE} of them is allowed"
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
            let (design, null) = the_study_of(&phenotype, &values, &kinship, &tested);
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
                .zip(model.projected_trait())
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
            let (design, _) = the_study_of(&phenotype, &values, &kinship, &tested);
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
                        .zip(model.projection().chunks_exact(num_individuals))
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
    /// It is the one fixture that takes the whole search through the
    /// boundary, and it takes it through twice: the trace divides by the
    /// variance, so a variance of 0 has to be answered from the weights
    /// and the diagonal of the kinship before the division; and a second
    /// variance at 0 is what ends the fit, where a fit that only set it to
    /// 0 would run its 200 steps there. Neither the derivative nor the
    /// average information is above 0 here, so the Newton step is not a
    /// number, and what keeps the variance a number is the bracket:
    /// a step that is not above 0 quarters the variance instead.
    ///
    /// Measured on 24 September 2026 on both backends: 12 steps and 12
    /// rounds, and the effects are the plain logistic null's, which is
    /// what a fit at a variance of 0 is.
    #[test]
    fn a_kinship_that_explains_nothing_lands_at_a_variance_of_zero() {
        let (phenotype, values, _) = the_panel("panel_called");
        let tested: Vec<usize> = (0..phenotype.len()).collect();
        let kinship = vec![0.0_f64; phenotype.len() * phenotype.len()];
        let (design, null) = the_study_of(&phenotype, &values, &kinship, &tested);
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
        let (design, _) = the_study_of(&phenotype, &values, &kinship, &tested);
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
        let (design, null) = the_study_of(&phenotype, &values, &kinship, &tested);
        let shorter = &kinship[..kinship.len() - 200];
        match TheLinearization::of_the_logistic_null(&phenotype, &design, shorter, &null) {
            Err(Error::GwasInputOfAnotherSize { .. }) => {}
            Ok(_) => panic!("a kinship of 39800 values was accepted"),
            Err(error) => panic!("a kinship of 39800 values gave {error}"),
        }
    }
}
