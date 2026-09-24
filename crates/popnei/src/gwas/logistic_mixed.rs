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
use super::study::{Design, GwasInputShape, GwasModel};

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
        if collapsed {
            return Err(Error::GwasFitDidNotSettle {
                model: GwasModel::Glmm,
                rounds: self.rounds,
            });
        }
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

/// The linearization of the logistic mixed model on the two panels of
/// `docs/specs/gwas.md`, at the variance of the kinship effect that GMMAT
/// fitted and at the two ends where it has none and where it has more than
/// any fit reaches.
#[cfg(test)]
mod glmm {
    use super::TheLinearization;
    use crate::error::Error;
    use crate::gwas::linear::lm::the_trait_and_the_design_of_the_panel;
    use crate::gwas::linear_mixed::lmm::{the_individuals_of_the_kinship, the_kinship_of};
    use crate::gwas::logistic::LogisticModel;
    use crate::gwas::study::{Design, GwasInput, TestType, TraitType};

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
    /// which the file gives rounded. Measured over the three on 25
    /// September 2026, the worst is 1.97e-6 on both backends, `cov2`,
    /// which is 20 per cent of what is allowed.
    const OF_GMMAT_EFFECTS: ([f64; 3], f64) = ([-1.416464, 0.753476, 1.583210], 1e-5);

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
