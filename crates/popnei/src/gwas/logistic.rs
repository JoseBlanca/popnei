//! The logistic model: a binomial trait and no kinship.
//!
//! [`LogisticModel`] is the null model fitted by iteratively reweighted
//! least squares and the two tests of every variant against it: the score
//! test, which needs nothing but that null, and the Wald test, which fits
//! one logistic regression per variant in the buffers of
//! [`TheFitOfOneVariant`]. "The logistic model" of `docs/specs/gwas.md`
//! says what each of them computes and how they are verified.

use popnei_linalg::{TheFirstOperand, TheSecondOperand};

use crate::error::{Error, Result};

use super::distributions::chi2_sf_1df;
use super::dosages::GwasDosages;
use super::result::{Answers, NullModel};
use super::study::{Design, GwasInputShape, GwasModel, TestType};
use super::the_share_that_is_nothing;

/// How many rounds of iteratively reweighted least squares the null model
/// is fitted in before a fit that is still moving is refused: 50, which is
/// `GLM_MAX_ITER` of `pynei/gwas.py` and which nobody has measured.
///
/// The panel of `docs/specs/gwas.md`, 200 individuals and two covariates,
/// settles in 5 rounds. What reaches 50 is a covariate that separates the
/// individuals that have the condition from the ones that have not, whose
/// effect has no finite value to settle at.
const ROUNDS_OF_THE_FIT: usize = 50;

/// How small the largest change in a coefficient has to be for the fit to
/// have settled: 1e-8, which is `GLM_TOL` of `pynei/gwas.py` and which
/// nobody has measured either. R's `glm` stops at 1e-8 of the deviance
/// instead, and the tolerances of "How it is verified" of "The logistic
/// model" of `docs/specs/gwas.md` are the distance between the two fits.
const THE_STEP_THAT_HAS_SETTLED: f64 = 1e-8;

/// How large the effect of a variant may be before its Wald fit has run
/// away: 30, past which the variant gets the three NaNs of one that has no
/// answer.
///
/// It is the threshold of `_wald_test` of `pynei/gwas.py`, inherited from
/// it and never measured. It is read on the effect of the variant alone
/// and not on the effects of the covariates, which "The logistic model" of
/// `docs/specs/gwas.md` says. The effect is a log odds ratio, so 30 is one
/// more copy of the allele multiplying the odds of the condition by 1e13,
/// and what walks there is a variant that separates the individuals that
/// have the condition from the ones that have not, whose effect has no
/// finite value to settle at.
const THE_EFFECT_THAT_HAS_RUN_AWAY: f64 = 30.0;

/// What the mean of the trait is moved by before the log of its odds
/// becomes the first intercept, 1e-6, from `_fit_logistic` of
/// `pynei/gwas.py`.
///
/// A trait of one value is refused before any model is fitted, so the mean
/// lies between 0 and 1 and the log is finite without it; it is kept
/// because "The logistic mixed model" of `docs/specs/gwas.md` starts its
/// own fit from this one's coefficients and a fit started elsewhere walks
/// another path.
const OF_THE_FIRST_INTERCEPT: f64 = 1e-6;

/// The logistic model of a study fitted without any variant in it, with
/// the buffers one block of variants is tested in.
///
/// The chance that an individual has the condition is a logistic curve in
/// the columns of the design, and the fit is iteratively reweighted least
/// squares: with `mu` the fitted chance of each individual, the weight of
/// that individual is `mu (1 - mu)`, and each round solves the design
/// weighted by those against the difference between the trait and `mu`.
/// That is "The logistic model" of `docs/specs/gwas.md`, and it is what
/// R's `glm` fits.
///
/// What the fit leaves for the variants is the weights, the difference
/// between the trait and `mu`, the design weighted by the weights and the
/// Cholesky factorization of `d' w d`, the design weighted by the weights
/// and taken against itself. The score test of a variant needs nothing
/// else, so no matrix is factored and nothing is inverted per variant. The
/// buffers of a block are kept from one block to the next, so this module
/// asks the machine for them once a pass and for none of them per variant.
/// What a pass still allocates per variant is under it, in the linear
/// algebra crate, and only on one of the two backends: on faer, which is
/// what a browser runs, `cholesky_lower` and `solve_with_cholesky` take a
/// scratch buffer on every call, and on BLAS and LAPACK neither does.
/// Counted on 24 September 2026 over one pass of 2000 variants with the
/// Wald test, the whole pass makes 87463 allocations on faer against 9447
/// on BLAS, so the 78016 that separate them are about 39 per variant and
/// all of them are in `crates/popnei-linalg/src/faer.rs`.
pub(crate) struct LogisticModel {
    /// The effect of the intercept and of each covariate, one per column
    /// of the design, as a log odds ratio.
    coefs: Vec<f64>,
    /// The design times those coefficients, one value per tested
    /// individual, which the chance that the individual has the condition
    /// is the logistic curve of.
    linear_predictor: Vec<f64>,
    /// `mu (1 - mu)` of each tested individual, which is how much a
    /// binomial trait of that fitted chance varies.
    weights: Vec<f64>,
    /// The trait less the fitted chance, one value per tested individual,
    /// which every variant is tested against.
    residuals: Vec<f64>,
    /// Each row of the design times the weight of its individual,
    /// `num_individuals` x `num_coefs`, row after row.
    weighted_design: Vec<f64>,
    /// The Cholesky factorization of `d' w d`, `num_coefs` x `num_coefs`,
    /// row after row, of which the lower half is read.
    factored: Vec<f64>,
    /// How many individuals the study tests.
    num_individuals: usize,
    /// How many columns the design has.
    num_coefs: usize,
    /// The dosages of a block against the weighted design, `x' w d`, the
    /// variants that have variance x `num_coefs`, which the solve against
    /// the factorization then overwrites with `(d' w d)⁻¹ d' w x`, the
    /// effect each column of the design has on the variant.
    of_the_design: Vec<f64>,
    /// What the design makes of each variant at those effects, the
    /// variants that have variance x `num_individuals`, which the
    /// denominator of the score test is formed from.
    explained: Vec<f64>,
    /// Each variant's dosages times the trait's residuals, one per variant
    /// that has variance.
    num: Vec<f64>,
    /// The effect of each variant that has variance.
    beta: Vec<f64>,
    /// The standard error of each of those effects.
    se: Vec<f64>,
    /// The p-value of each of those tests.
    p_value: Vec<f64>,
    /// The buffers the Wald test fits one variant in, made with the model
    /// when the study asked for that test and `None` when it asked for
    /// the score test, which reads none of them. They are the trait and
    /// two designs of the individuals by the coefficients with the
    /// variant, 0.96 MB at the 10000 individuals and 5 covariates of
    /// `docs/objectives.md`, and which of the two the model holds is what
    /// [`LogisticModel::test_the_block`] reads to know which test to make.
    of_a_variant: Option<TheFitOfOneVariant>,
}

impl LogisticModel {
    /// The logistic model of `phenotype` over `design`, fitted without any
    /// variant in it.
    ///
    /// `phenotype` holds 0.0 or 1.0 for each individual the design has a
    /// row for, which is what the study tests. Every coefficient starts at
    /// 0 but the intercept, which starts at the log of the odds of the
    /// mean of the trait, and the fit stops when the largest change in a
    /// coefficient falls below 1e-8.
    ///
    /// `test` is the test the study will make of every variant, which
    /// `the_model_and_the_test` of `study` chose before the fit. The Wald
    /// test fits one logistic regression per variant and the buffers it
    /// does that in are made here for it; a study that makes the score
    /// test is given none of them.
    ///
    /// # Errors
    ///
    /// [`Error::GwasInputOfAnotherSize`] when `phenotype` does not hold
    /// one value for each tested individual.
    /// [`Error::GwasFitDidNotSettle`] when the coefficients are still
    /// moving after 50 rounds, which is what a covariate that separates
    /// the individuals that have the condition from the ones that have not
    /// gives; when a round gives a step that is not finite; and when
    /// `d' w d` can no longer be factored, which the weights and the
    /// design itself both do, as
    /// [`LogisticModel::at_the_coefficients`] says. [`Error::GwasLinalg`]
    /// when one of the two products, the factorization for a reason that
    /// is not a singular matrix or the solve against it could not be done.
    pub(crate) fn of_the_study(
        phenotype: &[f64],
        design: &Design<'_>,
        test: TestType,
    ) -> Result<LogisticModel> {
        let num_individuals = design.num_individuals();
        let num_coefs = design.num_coefs();
        if phenotype.len() != num_individuals {
            return Err(Error::GwasInputOfAnotherSize {
                problem: GwasInputShape::Phenotype {
                    num_values: phenotype.len(),
                    num_individuals,
                },
            });
        }
        // `d' w d` is the columns of the design squared, and the design
        // itself holds more values than that: `Design::of_the_study`
        // refuses a study of no more individuals than its columns plus
        // one, so the rows are the columns plus two at least and
        // `num_coefs * num_coefs` is below the length of a slice that is
        // already in memory.
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "the design holds `num_individuals` times `num_coefs` values and \
                      `num_individuals` is `num_coefs` plus 2 at least, so this product is \
                      smaller than the length of the design and fits in a `usize`"
        )]
        let of_the_coefficients = num_coefs * num_coefs;
        let mut fitted = LogisticModel {
            coefs: vec![0.0_f64; num_coefs],
            linear_predictor: vec![0.0_f64; num_individuals],
            weights: vec![0.0_f64; num_individuals],
            residuals: vec![0.0_f64; num_individuals],
            // The design holds one row of its columns for each tested
            // individual, which `Design::of_the_study` has checked, so
            // this is `num_individuals` times `num_coefs` and the product
            // of the two is known to fit.
            weighted_design: vec![0.0_f64; design.values().len()],
            factored: vec![0.0_f64; of_the_coefficients],
            num_individuals,
            num_coefs,
            of_a_variant: match test {
                // The score test reads the weights, the residuals and the
                // factorization this fit leaves and nothing else, so a
                // study that makes it never builds these.
                TestType::Score => None,
                TestType::Wald => Some(TheFitOfOneVariant::of_the_study(phenotype, design)),
            },
            of_the_design: Vec::new(),
            explained: Vec::new(),
            num: Vec::new(),
            beta: Vec::new(),
            se: Vec::new(),
            p_value: Vec::new(),
        };
        // The mean of a binomial trait is the share of the tested
        // individuals that have the condition, and the first intercept is
        // the log of its odds. A trait of one value is refused before any
        // model is fitted, so the mean lies between 0 and 1.
        let mean = phenotype.iter().sum::<f64>() / num_individuals as f64;
        if let Some(intercept) = fitted.coefs.first_mut() {
            *intercept =
                ((mean + OF_THE_FIRST_INTERCEPT) / (1.0 - mean + OF_THE_FIRST_INTERCEPT)).ln();
        }
        let mut step = vec![0.0_f64; num_coefs];
        let mut settled = false;
        let mut rounds = 0_usize;
        while rounds < ROUNDS_OF_THE_FIT && !settled {
            rounds = rounds.saturating_add(1);
            fitted.at_the_coefficients(phenotype, design, rounds)?;
            fitted.the_step_of_the_fit(design, &mut step)?;
            // A step that is not a finite number ends the fit here, where
            // the Wald fit of a variant has the same guard. Added on 24
            // September 2026 as insurance and not for a case that has been
            // seen: no input of this repository reaches it, because the
            // factorization of the round refuses the system first. Without
            // it the step is added, the product of the next round refuses
            // an operand that is not finite, and the user is given
            // `GwasLinalg`, which is a `RuntimeError` in Python and so a
            // defect of popnei, for data of their own.
            if !step.iter().all(|step| step.is_finite()) {
                return Err(Error::GwasFitDidNotSettle {
                    model: GwasModel::Glm,
                    rounds,
                });
            }
            for (coef, step) in fitted.coefs.iter_mut().zip(&step) {
                *coef += *step;
            }
            settled = step
                .iter()
                .all(|step| step.abs() < THE_STEP_THAT_HAS_SETTLED);
        }
        if !settled {
            return Err(Error::GwasFitDidNotSettle {
                model: GwasModel::Glm,
                rounds,
            });
        }
        // The weights, the residuals and the factorization the variants
        // are tested against are the ones of the coefficients the fit
        // stopped at, and the last round left the ones of the round before
        // its own step.
        fitted.at_the_coefficients(phenotype, design, rounds)?;
        Ok(fitted)
    }

    /// The fitted chance of each tested individual at the coefficients the
    /// fit holds, and what is built from it: the weights, the trait less
    /// that chance, the design weighted by the weights and the Cholesky
    /// factorization of the weighted design taken against the design.
    ///
    /// `rounds` is how many rounds the fit has run when this is called,
    /// which the error of a fit that has run away names.
    ///
    /// # Errors
    ///
    /// [`Error::GwasFitDidNotSettle`] when `d' w d` can no longer be
    /// factored, which two different things do. The weights do it when the
    /// fit runs away: a weight falls to 0 when the chance the fit gives an
    /// individual reaches 0 or 1, and a covariate that separates the
    /// individuals that have the condition from the ones that have not is
    /// what takes them there. The design does it alone in the band between
    /// the two tolerances: `Design::of_the_study` refuses columns that are
    /// not independent at the tolerance numpy's rank uses, a Cholesky
    /// factorization's is tighter, and two covariates can be independent
    /// enough for the first and not for the second. "The logistic model"
    /// of `docs/specs/gwas.md` measures that band, and the message of the
    /// error names both causes because they have different remedies.
    /// [`Error::GwasLinalg`] when the product of the
    /// design with the coefficients, the product that gives `d' w d`, or
    /// the factorization for any other reason, could not be done.
    fn at_the_coefficients(
        &mut self,
        phenotype: &[f64],
        design: &Design<'_>,
        rounds: usize,
    ) -> Result<()> {
        let num_individuals = self.num_individuals;
        let num_coefs = self.num_coefs;
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: design.values(),
                rows: num_individuals,
            },
            num_coefs,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.coefs,
                cols: 1,
            },
            &mut self.linear_predictor,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design with the effects of the logistic null model",
            source,
        })?;
        for ((weight, residual), (predicted, measured)) in self
            .weights
            .iter_mut()
            .zip(&mut self.residuals)
            .zip(self.linear_predictor.iter().zip(phenotype))
        {
            let chance = the_chance_of(*predicted);
            *weight = chance * (1.0 - chance);
            *residual = measured - chance;
        }
        for ((weighted, weight), row) in self
            .weighted_design
            .chunks_exact_mut(num_coefs)
            .zip(&self.weights)
            .zip(design.values().chunks_exact(num_coefs))
        {
            for (value, of_the_design) in weighted.iter_mut().zip(row) {
                *value = weight * of_the_design;
            }
        }
        popnei_linalg::product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: design.values(),
                rows: num_coefs,
            },
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.weighted_design,
                cols: num_coefs,
            },
            &mut self.factored,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design with the design weighted by the logistic null model",
            source,
        })?;
        popnei_linalg::cholesky_lower(&mut self.factored, num_coefs).map_err(|source| {
            match matches!(source, popnei_linalg::Error::Singular { .. }) {
                true => Error::GwasFitDidNotSettle {
                    model: GwasModel::Glm,
                    rounds,
                },
                false => Error::GwasLinalg {
                    operation: "factorization of the weighted design of the logistic null model",
                    source,
                },
            }
        })
    }

    /// What the round adds to each coefficient: the solution of `d' w d`
    /// against `d' resid`, which is one weighted least squares fit of the
    /// difference between the trait and the fitted chance.
    ///
    /// `step` holds one value for each column of the design and is
    /// overwritten.
    ///
    /// # Errors
    ///
    /// [`Error::GwasLinalg`] when the product of the design with the
    /// residuals or the solve against the factorization could not be done.
    fn the_step_of_the_fit(&self, design: &Design<'_>, step: &mut [f64]) -> Result<()> {
        popnei_linalg::product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: design.values(),
                rows: self.num_coefs,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.residuals,
                cols: 1,
            },
            step,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design with the residuals of the logistic null model",
            source,
        })?;
        popnei_linalg::solve_with_cholesky(&self.factored, self.num_coefs, step, 1).map_err(
            |source| Error::GwasLinalg {
                operation: "solve of the weighted design of the logistic null model",
                source,
            },
        )
    }

    /// The effects the fit settled at, one for each column of the design,
    /// as log odds ratios.
    ///
    /// "The logistic mixed model" of `docs/specs/gwas.md` starts its own
    /// fit from these and from
    /// [`LogisticModel::linear_predictor`], and from nowhere else: a fit
    /// started elsewhere walks a different path through the bracket of the
    /// variance of the kinship effect and can stop at another value of it.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the logistic mixed model is what reads this, and the study refuses \
                      that model until its step on the variance of the kinship effect and \
                      its score test are written"
        )
    )]
    pub(super) fn coefs(&self) -> &[f64] {
        &self.coefs
    }

    /// The design times those effects, one value per tested individual,
    /// which the fitted chance of the individual is the logistic curve of.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the logistic mixed model is what reads this, and the study refuses \
                      that model until its step on the variance of the kinship effect and \
                      its score test are written"
        )
    )]
    pub(super) fn linear_predictor(&self) -> &[f64] {
        &self.linear_predictor
    }

    /// The null model of the result: the effects of the intercept and of
    /// the covariates, as log odds ratios.
    ///
    /// A logistic model has no residual variance, since the variance of a
    /// binomial trait is decided by its mean, and no kinship, so it has
    /// neither a genetic variance nor a heritability.
    #[must_use]
    pub(crate) fn null_model(&self, test: TestType) -> NullModel {
        NullModel {
            model: GwasModel::Glm,
            test,
            covariate_effects: self.coefs.clone(),
            residual_variance: None,
            genetic_variance: None,
            heritability: None,
            num_individuals: self.num_individuals,
        }
    }

    /// The test of every variant of a block that has variance among the
    /// tested individuals, in the order of the block.
    ///
    /// Both tests are made here, and a study makes the one it asked for:
    /// the Wald test, which is what it gets when it asks for none, fits
    /// one logistic regression per variant, and the score test reads the
    /// null model the study was fitted with and fits nothing. Which of the
    /// two it is was decided before the fit, and what says it here is
    /// whether the model holds the buffers the Wald test needs.
    ///
    /// `design` is the design the null model was fitted on, which the
    /// score test takes each variant through to form its denominator.
    ///
    /// # Errors
    ///
    /// Whatever [`LogisticModel::wald_test_the_block`] or
    /// [`LogisticModel::score_test_the_block`] fails with.
    pub(crate) fn test_the_block(
        &mut self,
        dosages: &GwasDosages,
        design: &Design<'_>,
    ) -> Result<Answers<'_>> {
        match self.of_a_variant.is_some() {
            true => self.wald_test_the_block(dosages),
            false => self.score_test_the_block(dosages, design),
        }
    }

    /// The Wald test of every variant of a block that has variance among
    /// the tested individuals, in the order of the block.
    ///
    /// The Wald test asks how far the effect the model fits for the
    /// variant is from 0 in the units of its own uncertainty, so it needs
    /// the fit: one logistic regression per variant, with the variant in
    /// the design beside the covariates, started at the coefficients of
    /// the null and an effect of 0 for the variant.
    /// [`TheFitOfOneVariant`] is what makes each of them, in buffers that
    /// were made with the model of a study that asked for this test, so
    /// this module allocates nothing for a variant; what the linear
    /// algebra crate allocates under it is on the doc comment of
    /// [`LogisticModel`], and on faer it is about 39 allocations per
    /// variant. A variant whose fit runs away gets the three NaNs of one
    /// that has no answer and the pass goes on.
    ///
    /// # Errors
    ///
    /// Whatever [`TheFitOfOneVariant::fit_and_test`] fails with, which is
    /// [`Error::GwasLinalg`] when a product, the factorization or a solve
    /// of one variant's fit could not be done.
    fn wald_test_the_block(&mut self, dosages: &GwasDosages) -> Result<Answers<'_>> {
        self.beta.clear();
        self.se.clear();
        self.p_value.clear();
        // The buffers are there, [`LogisticModel::test_the_block`] having
        // read that they are to come here at all.
        if let Some(of_a_variant) = self.of_a_variant.as_mut() {
            for of_the_variant in dosages.dosages().chunks_exact(self.num_individuals) {
                of_a_variant.with_the_dosages_of(of_the_variant);
                let (beta, se, p_value) = match of_a_variant.fit_and_test(&self.coefs)? {
                    TheAnswerOfAVariant::Answered { beta, se, p_value } => (beta, se, p_value),
                    TheAnswerOfAVariant::RanAway => (f64::NAN, f64::NAN, f64::NAN),
                };
                self.beta.push(beta);
                self.se.push(se);
                self.p_value.push(p_value);
            }
        }
        Ok(self.answers())
    }

    /// The score test of every variant of a block that has variance among
    /// the tested individuals, in the order of the block.
    ///
    /// The score test asks how steeply the likelihood rises at an effect
    /// of 0, so it needs nothing but the null model. With `x` the dosages
    /// of a variant, `w` the weights, `resid` the trait less the fitted
    /// chance and `d` the design, `num` is `x' resid` and `den` is
    /// `x' w x` less `(x' w d) (d' w d)⁻¹ (d' w x)`, which is the variant
    /// with the covariates taken out of it in the metric the weights give.
    /// The effect is `num / den`, a log odds ratio, it is uncertain by
    /// `1 / sqrt(den)`, and `num² / den` is read against a chi square with
    /// one degree of freedom. The covariates take the place of the
    /// projection matrix of the mixed models, and nothing is inverted per
    /// variant: the whole block is solved against the one factorization
    /// the null model left.
    ///
    /// `den` is formed and not subtracted, which is the second place
    /// popnei departs from pyNei's formula, the linear model's residual
    /// sum of squares being the first. With `b` the effects
    /// `(d' w d)⁻¹ d' w x` the design has on the variant, which the solve
    /// gives, the same quantity is the sum over the individuals of
    /// `w (x - d b)²`, a sum of terms that are 0 or above where the
    /// subtraction is two nearly equal numbers taken from each other.
    /// Measured on 24 September 2026 over six decades on 200 individuals
    /// with a covariate that is a variant's dosages and noise: at 1.3
    /// times the threshold below the subtracted form is out by 4.6e-3 of
    /// itself, and at the threshold by 29 per cent, so the guard would be
    /// reading a quantity whose error is larger than what it tests. It
    /// costs one product of the block per block, the effects against the
    /// design, into a buffer of the variants by the individuals, which is
    /// what the linear model already keeps for its residuals.
    ///
    /// A variant of which the design leaves at most the tested individuals
    /// times 2.2e-16 of `x' w x` has no answer, and gets the three NaNs a
    /// variant with no variance gets. For a variant that is a combination
    /// of the columns of the design what is left is the rounding of that
    /// sum and not a quantity, and `beta` would be a number divided by
    /// noise. It is **Open 2** of `docs/specs/gwas.md`.
    ///
    /// `design` is the design the null model was fitted on, the same one,
    /// which is what a variant is taken through.
    ///
    /// # Errors
    ///
    /// [`Error::GwasVariantsTooLarge`] when the values of the block are
    /// more than a `usize` counts, and [`Error::GwasLinalg`] when one of
    /// the three products or the solve could not be done, which is where a
    /// block of other individuals than the null model was fitted over is
    /// refused.
    fn score_test_the_block(
        &mut self,
        dosages: &GwasDosages,
        design: &Design<'_>,
    ) -> Result<Answers<'_>> {
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
        let of_the_design = num_vars
            .checked_mul(self.num_coefs)
            .ok_or(Error::GwasVariantsTooLarge)?;
        self.of_the_design.resize(of_the_design, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: dosages.dosages(),
                rows: num_vars,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.weighted_design,
                cols: self.num_coefs,
            },
            &mut self.of_the_design,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of a block of variants with the weighted design",
            source,
        })?;
        self.num.resize(num_vars, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: dosages.dosages(),
                rows: num_vars,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.residuals,
                cols: 1,
            },
            &mut self.num,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of a block of variants with the residuals of the null model",
            source,
        })?;
        // The solve overwrites what it is given, and what it was given,
        // `x' w d`, is not read again: what the rows hold afterwards is
        // the effects the columns of the design have on each variant.
        popnei_linalg::solve_with_cholesky(
            &self.factored,
            self.num_coefs,
            &mut self.of_the_design,
            num_vars,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "solve of a block of variants against the weighted design",
            source,
        })?;
        let values = num_vars
            .checked_mul(self.num_individuals)
            .ok_or(Error::GwasVariantsTooLarge)?;
        self.explained.resize(values, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &self.of_the_design,
                rows: num_vars,
            },
            self.num_coefs,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: design.values(),
                cols: self.num_individuals,
            },
            &mut self.explained,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the effects of the design on a block of variants with the \
                        design",
            source,
        })?;
        // The share of what the variant weighed that the design has to
        // leave of it for the variant to be worth testing.
        let share_that_is_nothing = the_share_that_is_nothing(self.num_individuals);
        for ((explained, num), of_the_variant) in self
            .explained
            .chunks_exact(self.num_individuals)
            .zip(&self.num)
            .zip(dosages.dosages().chunks_exact(self.num_individuals))
        {
            let weighted_length = of_the_variant
                .iter()
                .zip(&self.weights)
                .map(|(dosage, weight)| weight * dosage * dosage)
                .sum::<f64>();
            // What the design leaves of the variant, weighed by the
            // weights: a sum of terms that are 0 or above, formed from
            // what the design makes of the variant and not taken from
            // `x' w x` by subtracting the two nearly equal numbers the
            // doc comment measures.
            let den = explained
                .iter()
                .zip(of_the_variant.iter().zip(&self.weights))
                .map(|(explained, (dosage, weight))| {
                    let left = dosage - explained;
                    weight * left * left
                })
                .sum::<f64>();
            if den <= share_that_is_nothing * weighted_length {
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

/// What the Wald fit of one variant came to.
enum TheAnswerOfAVariant {
    /// The effect of the variant as a log odds ratio, how uncertain it is
    /// and the p-value of the test that it is 0.
    Answered {
        /// The effect the fit settled at.
        beta: f64,
        /// Its standard error.
        se: f64,
        /// The chance of a statistic as large as this one when the effect
        /// is 0.
        p_value: f64,
    },
    /// The fit ran away, and the variant gets the three NaNs of one that
    /// has no answer.
    RanAway,
}

/// What became of the square system of one round of a fit weighted by the
/// fitted chances, factored and solved against: `d' w d` of the design
/// with the variant in it for the Wald test of the logistic model, and
/// `d' sigma⁻¹ d` for the linearization of the logistic mixed model.
pub(super) enum TheSystemOfTheFit {
    /// It was factored, or solved, and the fit goes on.
    Worked,
    /// The linear algebra crate refused it as singular, or its pivots have
    /// collapsed, which is the fit running away: the columns of the design
    /// are independent and the weights are what took the matrix there, and
    /// they fall to 0 as the chances the fit gives reach 0 and 1.
    RanAway,
}

/// The chance that an individual has the condition at the linear predictor
/// `predicted`, which is the logistic curve of it.
///
/// A predictor far below 0 overflows the exponential and gives a chance of
/// 0, and one far above it gives 1; both are chances, and the weight of an
/// individual at either is 0.
pub(super) fn the_chance_of(predicted: f64) -> f64 {
    1.0 / (1.0 + (-predicted).exp())
}

/// Whether the system a round factored still has a direction left to move
/// the fit in, read off the pivots of the factorization `factored` of `n` x
/// `n`: the smallest of them against the largest.
///
/// A Cholesky factorization accepts a system whose smallest pivot is above
/// 0 by any margin, and a logistic fit whose weighted design has collapsed
/// can walk far past the point where its system says anything: the steps it
/// solves shrink because the system is nearly singular, not because the
/// coefficients have settled, so the fit declares itself settled while they
/// are still moving and answers with an effect it knows nothing about
/// beside a standard error of 1e7. **Open 5** of `docs/specs/gwas.md` has
/// the eight individuals that gave a p-value of 0.9999996 on one backend
/// and three NaNs on the other, and this is the rule it chose: the same
/// share of a scale that **Open 2** refuses a variant at,
/// `num_individuals` times 2.2e-16, with the largest pivot for the scale.
///
/// The pivots are the squares of the diagonal of what the factorization
/// wrote, since the matrix it factored is `l l'`.
pub(super) fn the_system_that_is_left(
    factored: &[f64],
    n: usize,
    num_individuals: usize,
) -> TheSystemOfTheFit {
    let mut smallest = f64::INFINITY;
    let mut largest = 0.0_f64;
    for (row, values) in factored.chunks_exact(n.max(1)).enumerate().take(n) {
        let diagonal = values.get(row).copied().unwrap_or(f64::NAN);
        let pivot = diagonal * diagonal;
        if pivot < smallest {
            smallest = pivot;
        }
        if pivot > largest {
            largest = pivot;
        }
    }
    match smallest <= the_share_that_is_nothing(num_individuals) * largest {
        true => TheSystemOfTheFit::RanAway,
        false => TheSystemOfTheFit::Worked,
    }
}

/// The buffers one variant is fitted and tested in by the Wald test, made
/// once for the study and written over for every variant of every block.
///
/// The fit is the same iteratively reweighted least squares as the null
/// model's, over the design with the dosages of the variant added to it as
/// its last column, started at the coefficients of the null and an effect
/// of 0 for the variant. What differs is what a fit that does not settle
/// means: the null model is the study, so a null fit that runs away is
/// [`Error::GwasFitDidNotSettle`] and there is nothing to report; one
/// variant of a million running away is an ordinary thing to find in a
/// study, so it gets the three NaNs that "The variants that have no
/// answer" of `docs/specs/gwas.md` gives every variant with no answer, and
/// the pass goes on to the next variant.
struct TheFitOfOneVariant {
    /// The trait of each tested individual, 0.0 or 1.0. The null fit is
    /// given it and does not keep it, and a Wald fit needs it at every
    /// block, so this is the model's own copy of it.
    phenotype: Vec<f64>,
    /// The design of the study with the dosages of the variant as its last
    /// column, `num_individuals` x `num_coefs`, row after row. The columns
    /// of the design are written once for the study and only the last one
    /// changes from one variant to the next.
    design: Vec<f64>,
    /// Each of its rows times the weight of that individual, of the same
    /// size.
    weighted_design: Vec<f64>,
    /// The effect of the intercept, of each covariate and of the variant,
    /// which the rounds move. The last of them is what the variant is
    /// answered with.
    coefs: Vec<f64>,
    /// What the round adds to each of those.
    step: Vec<f64>,
    /// The design with the variant times those coefficients, one value per
    /// tested individual, which the chance that the individual has the
    /// condition is the logistic curve of.
    linear_predictor: Vec<f64>,
    /// `mu (1 - mu)` of each tested individual at that chance.
    weights: Vec<f64>,
    /// The trait less that chance, one value per tested individual.
    residuals: Vec<f64>,
    /// `d' w d` of the design with the variant in it, `num_coefs` x
    /// `num_coefs` row after row, and its Cholesky factorization once a
    /// round has factored it.
    system: Vec<f64>,
    /// The last column of the identity, `num_coefs` values, which is
    /// solved against that factorization to read the last diagonal entry
    /// of the inverse: the variance of the effect of the variant.
    of_the_effect: Vec<f64>,
    /// How many individuals the study tests.
    num_individuals: usize,
    /// How many coefficients the fit has: one for each column of the
    /// design and one for the variant, so 2 at least.
    num_coefs: usize,
}

impl TheFitOfOneVariant {
    /// The buffers of the study, with the columns of the design already in
    /// the ones that hold it.
    fn of_the_study(phenotype: &[f64], design: &Design<'_>) -> TheFitOfOneVariant {
        let num_individuals = design.num_individuals();
        // `Design::of_the_study` refuses a study of no more individuals
        // than its columns plus one, so the individuals are the columns
        // plus two at least and one more column than the design has is
        // fewer columns than there are individuals. The design holds its
        // individuals times its columns values in memory, each of 8 bytes,
        // so that product is below an eighth of what a `usize` counts and
        // neither it plus the individuals, which is the size of the design
        // with the variant, nor the columns squared, which is smaller than
        // that, can overflow one.
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "the columns with the variant are fewer than the individuals, and the \
                      design of that many columns is the design of the study, already in \
                      memory, plus one value for each individual"
        )]
        let num_coefs = design.num_coefs() + 1;
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "the same bound: these two products are at most the values of the \
                      design of the study plus one for each individual"
        )]
        let of_the_design = num_individuals * num_coefs;
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "the columns with the variant are fewer than the individuals, so this \
                      is smaller than the product above"
        )]
        let of_the_system = num_coefs * num_coefs;
        let mut fitted = TheFitOfOneVariant {
            phenotype: phenotype.to_vec(),
            design: vec![0.0_f64; of_the_design],
            weighted_design: vec![0.0_f64; of_the_design],
            coefs: vec![0.0_f64; num_coefs],
            step: vec![0.0_f64; num_coefs],
            linear_predictor: vec![0.0_f64; num_individuals],
            weights: vec![0.0_f64; num_individuals],
            residuals: vec![0.0_f64; num_individuals],
            system: vec![0.0_f64; of_the_system],
            of_the_effect: vec![0.0_f64; num_coefs],
            num_individuals,
            num_coefs,
        };
        for (with_the_variant, of_the_study) in fitted
            .design
            .chunks_exact_mut(num_coefs)
            .zip(design.values().chunks_exact(design.num_coefs()))
        {
            for (value, of_the_design) in with_the_variant.iter_mut().zip(of_the_study) {
                *value = *of_the_design;
            }
        }
        fitted
    }

    /// Puts the dosages of one variant in the last column of the design,
    /// where the fit reads them, leaving the columns of the study's own
    /// design as they are.
    ///
    /// `dosages` holds one value for each tested individual. A shorter
    /// slice leaves the individuals after it with the dosages of the
    /// variant before, which is why the caller reads the block in chunks
    /// of the individuals the model was fitted over.
    fn with_the_dosages_of(&mut self, dosages: &[f64]) {
        for (row, dosage) in self.design.chunks_exact_mut(self.num_coefs).zip(dosages) {
            if let Some(of_the_variant) = row.last_mut() {
                *of_the_variant = *dosage;
            }
        }
    }

    /// The Wald test of the variant whose dosages the design holds: its
    /// effect, how uncertain that effect is and its p-value, or the fit
    /// running away.
    ///
    /// `of_the_null` holds the effects the null model settled at, one for
    /// each column of the study's design, which the fit starts from with
    /// an effect of 0 for the variant. Each round solves the system of the
    /// coefficients plus the variant, and three things mark a fit that has
    /// run away, which "The logistic model" of `docs/specs/gwas.md` names:
    /// a step that is not finite, the system the factorization or the
    /// solve refuses as singular among them; the effect of the variant
    /// passing 30 in absolute value; and a fit still moving after the 50
    /// rounds it is given. A variant that separates the individuals that
    /// have the condition from the ones that have not is what these catch,
    /// and its effect is infinite rather than large.
    ///
    /// The first of the three is read as a value that is **not finite**
    /// and never as an infinity, which is what makes the two linear
    /// algebra backends mark the same variants: `docs/specs/linalg.md`
    /// lets through a diagonal entry whose reciprocal overflows, and there
    /// LAPACK gives an infinity where faer gives a NaN, each reporting
    /// success. The standard error is read the same way, since it comes
    /// out of a solve against that same factorization.
    ///
    /// The standard error is read off one more weighting of the design at
    /// the coefficients the fit stopped at, because the last round left
    /// the system of the coefficients before its own step. `_wald_test` of
    /// `pynei/gwas.py` fits a whole block at once and builds the system of
    /// every variant at every round, so a variant that settled while
    /// others were still moving has the system of its own last
    /// coefficients there too; measured on 24 September 2026 over the 1199
    /// variants of the panel that both libraries answer, the two standard
    /// errors are within 3.369e-15 of each other as a share of themselves,
    /// the worst being the 3.368e-15 of `var1151`.
    ///
    /// # Errors
    ///
    /// [`Error::GwasLinalg`] when one of the products, the factorization
    /// or a solve could not be done for a reason that is not the matrix
    /// being singular, which is the fit running away and not an error.
    fn fit_and_test(&mut self, of_the_null: &[f64]) -> Result<TheAnswerOfAVariant> {
        for (coef, of_the_null) in self.coefs.iter_mut().zip(of_the_null) {
            *coef = *of_the_null;
        }
        if let Some(of_the_variant) = self.coefs.last_mut() {
            *of_the_variant = 0.0;
        }
        let mut settled = false;
        let mut rounds = 0_usize;
        while rounds < ROUNDS_OF_THE_FIT && !settled {
            rounds = rounds.saturating_add(1);
            match self.at_the_coefficients()? {
                TheSystemOfTheFit::Worked => {}
                TheSystemOfTheFit::RanAway => return Ok(TheAnswerOfAVariant::RanAway),
            }
            match self.the_step_of_the_round()? {
                TheSystemOfTheFit::Worked => {}
                TheSystemOfTheFit::RanAway => return Ok(TheAnswerOfAVariant::RanAway),
            }
            if !self.step.iter().all(|step| step.is_finite()) {
                return Ok(TheAnswerOfAVariant::RanAway);
            }
            for (coef, step) in self.coefs.iter_mut().zip(&self.step) {
                *coef += *step;
            }
            if self.the_effect().abs() > THE_EFFECT_THAT_HAS_RUN_AWAY {
                return Ok(TheAnswerOfAVariant::RanAway);
            }
            settled = self
                .step
                .iter()
                .all(|step| step.abs() < THE_STEP_THAT_HAS_SETTLED);
        }
        if !settled {
            return Ok(TheAnswerOfAVariant::RanAway);
        }
        // The variance of the effect is read off the weights of the
        // coefficients the fit stopped at, and the last round left the
        // ones of the round before its own step.
        match self.at_the_coefficients()? {
            TheSystemOfTheFit::Worked => {}
            TheSystemOfTheFit::RanAway => return Ok(TheAnswerOfAVariant::RanAway),
        }
        let Some(se) = self.the_error_of_the_effect()? else {
            return Ok(TheAnswerOfAVariant::RanAway);
        };
        let beta = self.the_effect();
        let statistic = beta / se;
        Ok(TheAnswerOfAVariant::Answered {
            beta,
            se,
            p_value: chi2_sf_1df(statistic * statistic),
        })
    }

    /// The effect the fit has given the variant, which is the last of its
    /// coefficients. The fit has one for each column of the design and one
    /// for the variant, so there is always a last one.
    fn the_effect(&self) -> f64 {
        self.coefs.last().copied().unwrap_or(f64::NAN)
    }

    /// The fitted chance of each tested individual at the coefficients the
    /// fit holds, and what is built from it: the weights, the trait less
    /// that chance, the design weighted by the weights and the Cholesky
    /// factorization of `d' w d`, whose pivots are then read by
    /// [`the_system_that_is_left`].
    ///
    /// # Errors
    ///
    /// [`Error::GwasLinalg`] when the product of the design with the
    /// coefficients, the product that gives `d' w d`, or the factorization
    /// for a reason that is not a singular matrix, could not be done.
    fn at_the_coefficients(&mut self) -> Result<TheSystemOfTheFit> {
        let num_individuals = self.num_individuals;
        let num_coefs = self.num_coefs;
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &self.design,
                rows: num_individuals,
            },
            num_coefs,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.coefs,
                cols: 1,
            },
            &mut self.linear_predictor,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design of a variant with the effects of its logistic fit",
            source,
        })?;
        for ((weight, residual), (predicted, measured)) in self
            .weights
            .iter_mut()
            .zip(&mut self.residuals)
            .zip(self.linear_predictor.iter().zip(&self.phenotype))
        {
            let chance = the_chance_of(*predicted);
            *weight = chance * (1.0 - chance);
            *residual = measured - chance;
        }
        for ((weighted, weight), row) in self
            .weighted_design
            .chunks_exact_mut(num_coefs)
            .zip(&self.weights)
            .zip(self.design.chunks_exact(num_coefs))
        {
            for (value, of_the_design) in weighted.iter_mut().zip(row) {
                *value = weight * of_the_design;
            }
        }
        popnei_linalg::product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &self.design,
                rows: num_coefs,
            },
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.weighted_design,
                cols: num_coefs,
            },
            &mut self.system,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design of a variant with that design weighted by its \
                        logistic fit",
            source,
        })?;
        match the_system_of(
            popnei_linalg::cholesky_lower(&mut self.system, num_coefs),
            "factorization of the weighted design of the logistic fit of a variant",
        )? {
            TheSystemOfTheFit::Worked => Ok(the_system_that_is_left(
                &self.system,
                num_coefs,
                self.num_individuals,
            )),
            TheSystemOfTheFit::RanAway => Ok(TheSystemOfTheFit::RanAway),
        }
    }

    /// What the round adds to each coefficient: the solution of `d' w d`
    /// against `d' resid`, which the fit then adds to the coefficients it
    /// holds. It is left in [`TheFitOfOneVariant::step`].
    ///
    /// # Errors
    ///
    /// [`Error::GwasLinalg`] when the product of the design with the
    /// residuals, or the solve against the factorization for a reason that
    /// is not a singular matrix, could not be done.
    fn the_step_of_the_round(&mut self) -> Result<TheSystemOfTheFit> {
        popnei_linalg::product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &self.design,
                rows: self.num_coefs,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.residuals,
                cols: 1,
            },
            &mut self.step,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design of a variant with the residuals of its logistic fit",
            source,
        })?;
        the_system_of(
            popnei_linalg::solve_with_cholesky(&self.system, self.num_coefs, &mut self.step, 1),
            "solve of the weighted design of the logistic fit of a variant",
        )
    }

    /// How uncertain the effect of the variant is: the square root of the
    /// last diagonal entry of the inverse of `d' w d`, and `None` when the
    /// fit has run away.
    ///
    /// The entry is read by solving the factorization against the last
    /// column of the identity, which gives that column of the inverse and
    /// costs the square of the coefficients where inverting the whole
    /// matrix costs their cube. What comes back is refused when it is not
    /// a finite number above 0: the one way a factorization that was
    /// accepted gives anything else is a diagonal entry whose reciprocal
    /// overflows, where the two backends of `docs/specs/linalg.md` differ
    /// in which of an infinity and a NaN they answer with, and both fail
    /// this.
    ///
    /// # Errors
    ///
    /// [`Error::GwasLinalg`] when the solve could not be done for a reason
    /// that is not a singular matrix.
    fn the_error_of_the_effect(&mut self) -> Result<Option<f64>> {
        self.of_the_effect.fill(0.0);
        if let Some(of_the_variant) = self.of_the_effect.last_mut() {
            *of_the_variant = 1.0;
        }
        match the_system_of(
            popnei_linalg::solve_with_cholesky(
                &self.system,
                self.num_coefs,
                &mut self.of_the_effect,
                1,
            ),
            "solve of the last column of the identity against the logistic fit of a variant",
        )? {
            TheSystemOfTheFit::Worked => {}
            TheSystemOfTheFit::RanAway => return Ok(None),
        }
        let variance = self.of_the_effect.last().copied().unwrap_or(f64::NAN);
        match variance.is_finite() && variance > 0.0 {
            true => Ok(Some(variance.sqrt())),
            false => Ok(None),
        }
    }
}

/// What a factorization or a solve of a fit weighted by the fitted chances
/// came to: the matrix the weights made is refused as singular when the
/// fit has run away, and every other refusal is an error.
///
/// `operation` says what was being done, which the error names.
///
/// # Errors
///
/// [`Error::GwasLinalg`] when the linear algebra crate refused the call
/// for a reason that is not a singular matrix.
pub(super) fn the_system_of(
    what_it_gave: popnei_linalg::Result<()>,
    operation: &'static str,
) -> Result<TheSystemOfTheFit> {
    match what_it_gave {
        Ok(()) => Ok(TheSystemOfTheFit::Worked),
        Err(popnei_linalg::Error::Singular { .. }) => Ok(TheSystemOfTheFit::RanAway),
        Err(source) => Err(Error::GwasLinalg { operation, source }),
    }
}

/// The logistic model against R's score test on the panel with every
/// genotype called, and the three cases of its own: a variant the design
/// leaves nothing of, a fit that does not settle, and the Wald test that
/// is not written yet.
#[cfg(test)]
mod glm {
    use super::TheFitOfOneVariant;
    use crate::block::BlockReader;
    use crate::error::Error;
    use crate::gwas::linear::lm::{
        THE_HEADER_OF_EIGHT, reader_over, the_panel_path, the_reference_path, the_study_of,
        the_trait_and_the_design_of_the_panel,
    };
    use crate::gwas::result::Gwas;
    use crate::gwas::study::{GwasInput, GwasModel, TestType, TraitType};
    use crate::io::vcf::{VcfOptions, VcfReader};

    /// How far the score statistic of one of the six variants may be from
    /// R's: 1e-3 absolute.
    ///
    /// It is the bound of "How it is verified" of "The logistic model" of
    /// `docs/specs/gwas.md`, and it measures the distance between two
    /// fits: R's `glm` stops at 1e-8 of the deviance and popnei's
    /// iteratively reweighted least squares at 1e-8 of the largest change
    /// in a coefficient, so the two stop at different coefficients and
    /// their score statistics differ by more than either arithmetic.
    /// `tests/reference/gwas/r.panel_called.glm.score.tsv` is written at
    /// full precision, so none of it is spent on the printing.
    ///
    /// It is absolute where the statistics of the six run from 1.48 to
    /// 12.48, which "How it is verified" says is the safe way round for
    /// this check to be fragile: on a panel whose statistics were far
    /// larger it would fail a right answer rather than pass a wrong one.
    ///
    /// Measured over the six on 24 September 2026, the same to ten digits
    /// on Accelerate and on faer: the worst is 6.024e-4, `var0629`, which
    /// is 60 per cent of what is allowed, and the next is the 4.521e-4 of
    /// `var0052`. The bound is not lowered to two or three times that, as
    /// a bound on popnei's own arithmetic would be: what it measures is
    /// where R's fit stopped, and the spec fixes it. What popnei's own
    /// arithmetic is worth here is the distance between the two backends,
    /// which is largest at `var0052`, 5.3e-14 of a statistic of 12.48, ten
    /// orders of magnitude below the distance from R.
    const OF_R: f64 = 1e-3;

    /// How far the p-value of one of the six may be from R's: 1e-3 of
    /// itself in `log10`, which is the spec's bound and which comes from
    /// the same distance between the two fits.
    ///
    /// A p-value is compared in `log10` and not as a share of itself
    /// because it runs over orders of magnitude, which "How it is
    /// verified" of "What every model shares" asks for. Measured over the
    /// six on 24 September 2026 on both backends: the worst is 1.430e-4,
    /// `var0629` again, which is 14 per cent of what is allowed.
    const OF_R_IN_LOG10: f64 = 1e-3;

    /// What R 4.6.1's `anova(glm, test = "Rao")` wrote for six variants of
    /// the panel with every genotype called, in
    /// `tests/reference/gwas/r.panel_called.glm.score.tsv`, with `cov1`
    /// and `cov2` as covariates: the id, the score statistic and its
    /// p-value, which are the two quantities "How it is verified" of "The
    /// logistic model" of `docs/specs/gwas.md` compares. Five of them are
    /// the causal variants of `causal_vars.csv` and `var0000` is not
    /// causal.
    ///
    /// R reports the statistic and popnei reports `beta` and `se`, so what
    /// the test compares against the first number is `(beta / se)²`. The
    /// two libraries count the dosages of a variant from different
    /// alleles, R from the alternative one and popnei from the one that is
    /// not the major one among the tested individuals, which turns the
    /// sign of `beta` over for some variants and leaves the statistic and
    /// the p-value as they are.
    const OF_R_SIX: [(&str, f64, f64); 6] = [
        ("var0000", 4.938245, 0.026_268_700),
        ("var0052", 12.484427, 0.000_410_359),
        ("var0629", 9.165576, 0.002_466_100),
        ("var0751", 1.480401, 0.223_711_736),
        ("var1137", 2.911424, 0.087_954_199),
        ("var1188", 10.382961, 0.001_271_835),
    ];

    /// The score test of the six variants of the panel is R's: the
    /// statistic within 1e-3 absolute and the p-value within 1e-3 in
    /// `log10`, which is what the spec asks of the whole columns in
    /// Python.
    ///
    /// The null model comes back with the effects of the intercept and of
    /// the two covariates and with none of the three variances, which is
    /// what a logistic model without a kinship has.
    #[test]
    fn the_six_variants_of_the_panel_are_rs_score_statistic_and_p_value() {
        let result = the_logistic_study_of_the_panel(Some(TestType::Score));
        assert_eq!(result.num_vars, 1200, "the variants of the panel");
        assert_eq!(result.null_model.model, GwasModel::Glm);
        assert_eq!(result.null_model.test, TestType::Score);
        assert_eq!(result.null_model.num_individuals, 200);
        assert_eq!(
            result.null_model.covariate_effects.len(),
            3,
            "the intercept and the two covariates"
        );
        assert_eq!(
            result.null_model.residual_variance, None,
            "a binomial trait has no residual variance of its own"
        );
        assert_eq!(result.null_model.genetic_variance, None, "no kinship");
        assert_eq!(result.null_model.heritability, None, "no kinship");
        let ids = result.ids.as_deref().expect("the ids of the variants");
        for (id, statistic, p_value) in OF_R_SIX {
            let var = match ids.iter().position(|held| held == id) {
                Some(var) => var,
                None => panic!("{id} is not a variant of the panel"),
            };
            let of_popnei = result.beta[var] / result.se[var];
            let of_popnei = of_popnei * of_popnei;
            let difference = (of_popnei - statistic).abs();
            assert!(
                difference <= OF_R,
                "the score statistic of {id} is {of_popnei} and R gives {statistic}, \
                 {difference} away against the {OF_R} allowed"
            );
            let in_log10 = (result.p_value[var] / p_value).log10().abs();
            assert!(
                in_log10 <= OF_R_IN_LOG10,
                "the p-value of {id} is {found} and R gives {p_value}, {in_log10} away in \
                 log10 against the {OF_R_IN_LOG10} allowed",
                found = result.p_value[var]
            );
        }
    }

    /// The design of the fixture of **Open 2** below: the intercept and
    /// one covariate, which is the dosages of the first variant of that
    /// fixture, `1 0 2 0 1 2 0 1`, in units a tenth of theirs, row after
    /// row.
    ///
    /// The covariate carries the same information as the variant whatever
    /// it is multiplied by, and the fit gives it an effect ten times as
    /// large to say the same thing; what the tenth changes is the
    /// arithmetic: the denominator of the first variant's score test lands
    /// on 1.891e-31 on Accelerate and 9.565e-31 on faer at a tenth,
    /// against 1.199e-32 and 1.953e-31 at a covariate equal to the
    /// dosages, all four measured on 24 September 2026 and all four
    /// positive and far under the threshold of **Open 2**, which is what
    /// this fixture is for.
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

    /// The binomial trait of those eight individuals. Five have the
    /// condition and three have not, and neither the intercept alone nor
    /// the covariate separates the two groups, so the fit settles.
    const THE_TRAIT_OF_EIGHT: [f64; 8] = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0];

    /// All eight of them are tested.
    const THE_INDIVIDUALS_OF_EIGHT: [usize; 8] = [0, 1, 2, 3, 4, 5, 6, 7];

    /// What numpy 2.5.3 answers for the second variant of that fixture,
    /// the dosages `0 0 0 0 2 2 2 2`, fitted and tested by the formulas of
    /// "The logistic model" of `docs/specs/gwas.md` with an explicit
    /// inverse of `d' w d`: the effect, the standard error and the
    /// p-value, which is scipy 1.16.2's `chi2.sf` of the statistic
    /// 1.0665823342310243. The threshold that refuses the first variant
    /// leaves this one alone.
    const OF_THE_ORDINARY_VARIANT: (f64, f64, f64) = (
        0.855_370_693_019_237_7,
        0.828_241_853_958_789_1,
        0.301_718_693_520_239_3,
    );

    /// How far the effect and the standard error of that variant may be
    /// from numpy's: 1.5e-15 times the standard error numpy gives it,
    /// 0.8282.
    ///
    /// Both are measured against that standard error and neither against
    /// itself, which is what "How it is verified" of "What every model
    /// shares" of `docs/specs/gwas.md` asks: an effect is what the study
    /// measures and its standard error is the scale it measures it on, and
    /// an effect that cancelled to near 0, which most of a study's do, is
    /// no guide to its own error. The effect of this fixture is 0.8554 and
    /// its standard error 0.8282, so the two bounds happen to be the same
    /// number here; the form is the rule and the fixture is not what
    /// decides it.
    ///
    /// Lowered until it failed on both backends on 24 September 2026: it
    /// breaks at 5e-16, where the effect is 6.70e-16 of the standard error
    /// away on faer, and this is three times that. The other two
    /// measurements are the effect on Accelerate and the standard error on
    /// either, both 2.68e-16, so the whole of what is measured here is a
    /// bit or two of a value of 0.83.
    const OF_NUMPY: f64 = 1.5e-15;

    /// How far the p-value of that variant may be from numpy's, in
    /// `log10`: 1e-14.
    ///
    /// A p-value runs over orders of magnitude, so it is compared in
    /// `log10`, which is already a scale, and not as a share of itself,
    /// which is the same rule of the same spec item.
    ///
    /// Lowered until it failed on both backends on 24 September 2026: it
    /// breaks at 4e-15, where it is 4.44e-15 away on Accelerate and
    /// 4.63e-15 on faer, and this is 2.5 times that. It is an order of
    /// magnitude above what the effect and the standard error spend
    /// because the last step of the test is the distribution and not the
    /// arithmetic: popnei reads the chi square off `erfc` of the `libm`
    /// crate and scipy 1.16.2 computes the same function another way, and
    /// "The two distributions" of `docs/specs/gwas.md` measures the two
    /// 1e-12 of each other apart over its own sample.
    const OF_NUMPYS_P_VALUE: f64 = 1e-14;

    /// A variant that the design leaves nothing of has no answer, which is
    /// the meanwhile of **Open 2** of `docs/specs/gwas.md`.
    ///
    /// The covariate is the first variant's dosages in units a tenth of
    /// theirs, which is what a user gets by putting a genotype in as a
    /// covariate. The design then explains the whole of that variant and
    /// what is left of it is rounding. Measured on 24 September 2026, the
    /// first reproduction of this rule at either score test, with the
    /// denominator formed as the weighted squared length of what the
    /// design leaves of the variant: it comes to 1.891e-31 on Accelerate
    /// and to 9.565e-31 on faer, where `x' w x`, the weighted squared
    /// length the variant had before the covariates were taken out, is
    /// 2.3592 and the threshold that refuses it, the eight tested
    /// individuals times 2.2e-16 of that, is 4.19e-15. Both are positive
    /// and both are sixteen orders of magnitude below the threshold, so
    /// this fixture is what exercises the comparison on both backends; the
    /// subtracted denominator it had before gave 4.44e-16 on Accelerate
    /// and exactly 0 on faer, where a comparison with 0 alone would have
    /// passed.
    ///
    /// What the variant is answered with when the threshold is taken out,
    /// measured the same day on both backends, is why it is there. The
    /// numerator is 0 to the bit, the design being what the fit made its
    /// residuals at right angles to, so the row is a `beta` of 0 and a
    /// p-value of 1 beside an `se` of 2.299e15 on Accelerate and
    /// 1.022e15 on faer, which a user reads as a variant that was tested
    /// and showed nothing and which no filter on a missing effect takes
    /// out.
    ///
    /// The second variant of the fixture is ordinary and is answered: it
    /// is what says that the threshold refuses the first variant and not
    /// every variant of the block.
    #[test]
    fn a_variant_the_design_leaves_nothing_of_has_no_answer() {
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
            kinship: None,
            test: Some(TestType::Score),
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS_OF_EIGHT,
            transform_to_biallelic: false,
        };
        let mut reader = reader_over(vcf.as_bytes());
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of a variant that is a covariate: {error}"),
        };

        assert_eq!(result.num_vars, 2);
        // The variant is in the result with its frequency, as every
        // variant that has no answer is: four of its sixteen alleles are
        // the one that is not the major one.
        assert!(
            (result.allele_freq[0] - 0.4375).abs() <= 1e-12,
            "the frequency of v0 is {found}",
            found = result.allele_freq[0]
        );
        assert!(
            result.beta[0].is_nan() && result.se[0].is_nan() && result.p_value[0].is_nan(),
            "the variant the design leaves nothing of was answered with a beta of {beta}, \
             an se of {se} and a p-value of {p_value}",
            beta = result.beta[0],
            se = result.se[0],
            p_value = result.p_value[0]
        );
        let (beta, se, p_value) = OF_THE_ORDINARY_VARIANT;
        // The effect and the standard error are measured against the
        // standard error, the scale of what the study estimates, and never
        // against the effect, which is what "How it is verified" of "What
        // every model shares" of `docs/specs/gwas.md` asks of every
        // comparison of the two.
        for (found, expected, what) in [
            (result.beta[1], beta, "the effect of v1"),
            (result.se[1], se, "the standard error of v1"),
        ] {
            let difference = (found - expected).abs();
            assert!(
                difference <= OF_NUMPY * se,
                "{what} is {found} and numpy gives {expected}, {difference} away, which is \
                 {share} of the {se} it is uncertain by against the {OF_NUMPY} allowed",
                share = difference / se
            );
        }
        // A p-value runs over orders of magnitude, so it is measured in
        // `log10`, which is already a scale.
        let in_log10 = (result.p_value[1] / p_value).log10().abs();
        assert!(
            in_log10 <= OF_NUMPYS_P_VALUE,
            "the p-value of v1 is {found} and numpy gives {p_value}, {in_log10} away in \
             log10 against the {OF_NUMPYS_P_VALUE} allowed",
            found = result.p_value[1]
        );
    }

    /// A covariate that separates the individuals that have the condition
    /// from the ones that have not is refused, naming the model and the
    /// rounds the fit ran.
    ///
    /// The covariate is 0 to 7 and the four individuals above 3 are the
    /// ones with the condition, so the effect that fits the trait is an
    /// infinite one and the fit walks towards it. numpy 2.5.3 on the same
    /// fixture, solving each round with an LU factorization as pyNei does,
    /// runs the 50 rounds and is still moving, at -299.6 and 84.0. popnei
    /// solves with a Cholesky, which refuses the weighted design once the
    /// fitted chances have reached 0 and 1, and that comes first: measured
    /// on 24 September 2026, at the round 45 on Accelerate and at the round
    /// 43 on faer. The test asserts the rounds are between 1 and the 50 the
    /// fit is given and not that they are 45, since which round the weights
    /// underflow at is the last bits of an exponential and a platform may
    /// put it elsewhere; what the fixture guards is that such a study is
    /// refused and not answered.
    ///
    /// The study is refused before any variant is read, so the one variant
    /// of the VCF is there only to make it a study.
    #[test]
    fn a_null_fit_that_has_not_settled_is_refused() {
        let mut vcf = String::from(THE_HEADER_OF_EIGHT);
        vcf.push_str("1\t1000\tv0\tA\tT\t.\t.\t.\tGT");
        for genotype in ["0/1", "0/0", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1"] {
            vcf.push('\t');
            vcf.push_str(genotype);
        }
        vcf.push('\n');
        let separating = [0.0_f64, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let design: Vec<f64> = (0..8_usize).flat_map(|row| [1.0, row as f64]).collect();
        let study = GwasInput {
            phenotype: &separating,
            trait_type: TraitType::Binomial,
            design: &design,
            num_coefs: 2,
            kinship: None,
            test: Some(TestType::Score),
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS_OF_EIGHT,
            transform_to_biallelic: false,
        };
        let mut reader = reader_over(vcf.as_bytes());
        match the_study_of(&mut reader, &study) {
            Err(Error::GwasFitDidNotSettle { model, rounds }) => {
                assert_eq!(model, GwasModel::Glm, "the model that was being fitted");
                assert!(
                    (1..=50).contains(&rounds),
                    "the fit ran {rounds} rounds of the 50 it is given"
                );
                let said = Error::GwasFitDidNotSettle { model, rounds }.to_string();
                assert!(
                    said.contains("separates"),
                    "the fit was refused with {said}"
                );
            }
            Err(error) => panic!("a covariate that separates the two groups: {error}"),
            Ok(result) => panic!("a study of {} variants was run", result.num_vars),
        }
    }

    /// A Wald fit whose system has collapsed has no answer, although its
    /// steps have fallen below the tolerance and none of the other marks
    /// has fired.
    ///
    /// It is the fixture of **Open 5** of `docs/specs/gwas.md`: eight
    /// individuals, the trait `0 0 0 1 0 1 1 1`, one covariate 0 to 7 and
    /// a variant of dosages `2 1 0 1 2 0 1 0`. No finite effect fits that
    /// variant, and pyNei, which solves each round with an LU
    /// factorization, runs its 50 rounds and is still moving at an
    /// intercept of -38.64, a covariate of 18.97 and an effect of -18.97,
    /// so it marks the fit by the round count and gives the three NaNs.
    /// popnei reached none of the five marks before this rule: the effect
    /// is -18.97, which is under the 30, and it is the intercept that
    /// passes it, which "The logistic model" says no mark reads; and the
    /// nearly singular system makes the steps shrink while the
    /// coefficients are still walking, so the fit declared itself settled
    /// at the round before the count ran out. Measured on 24 September
    /// 2026 through the Python package, it answered `beta`
    /// -18.235836883901765, `se` 38745320.69540999 and `p_value`
    /// 0.9999996244683889 on Accelerate and three NaNs on faer, which is
    /// the two builds disagreeing with the browser in the right.
    ///
    /// What catches it is the pivot of the factorization that has fallen
    /// to the share of the largest that **Open 2** calls nothing. What
    /// says that the rule refuses this fit and not every fit is the panel,
    /// whose 1199 answered variants stay answered: measured on both
    /// backends on 24 September 2026, the smallest pivot of a fit either
    /// reference panel answers is 2.600e-2 of the largest on the panel
    /// with every genotype called and 4.730e-5 on the panel with 3
    /// genotypes missing in 100, against a threshold there of 4.44e-14,
    /// which is nine orders of magnitude of headroom. This fixture has one
    /// variant and no second one, because the covariate it needs nearly
    /// separates the eight individuals on its own and every variant put
    /// beside it runs away too, by this mark or by the marks that were
    /// there before it.
    #[test]
    fn a_wald_fit_whose_system_has_collapsed_has_no_answer() {
        let mut vcf = String::from(THE_HEADER_OF_EIGHT);
        vcf.push_str("1\t1000\tv0\tA\tT\t.\t.\t.\tGT");
        // the dosages 2 1 0 1 2 0 1 0, allele by allele
        for genotype in ["1/1", "0/1", "0/0", "0/1", "1/1", "0/0", "0/1", "0/0"] {
            vcf.push('\t');
            vcf.push_str(genotype);
        }
        vcf.push('\n');
        let phenotype = [0.0_f64, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0];
        let design: Vec<f64> = (0..8_usize).flat_map(|row| [1.0, row as f64]).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Binomial,
            design: &design,
            num_coefs: 2,
            kinship: None,
            test: Some(TestType::Wald),
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS_OF_EIGHT,
            transform_to_biallelic: false,
        };
        let mut reader = reader_over(vcf.as_bytes());
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of a fit that stops before it has settled: {error}"),
        };
        assert_eq!(result.num_vars, 1);
        assert!(
            result.beta[0].is_nan() && result.se[0].is_nan() && result.p_value[0].is_nan(),
            "the variant whose fit collapsed was answered with a beta of {beta}, an se of \
             {se} and a p-value of {p_value}",
            beta = result.beta[0],
            se = result.se[0],
            p_value = result.p_value[0]
        );
    }

    /// The fixture of a fit that settles at an effect past 30: eight
    /// individuals whose trait is `1 0 0 0 0 1 1 1`, whose covariate is
    /// `1 4 2 4 4 2 3 0` and whose variant has the dosages
    /// `1 2 0 1 1 1 1 0`, allele by allele.
    ///
    /// The trait, the covariate and the genotypes were found by running
    /// 200000 fixtures of eight individuals drawn at random through this
    /// model on 24 September 2026 with the mark taken out, and keeping the
    /// ones the fit then answered with an effect past 30. Two of the
    /// 200000 do it.
    const THE_TRAIT_OF_A_FIT_THAT_PASSES_THIRTY: [f64; 8] =
        [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0];

    /// The design of that fixture: the intercept and that covariate, row
    /// after row.
    const THE_DESIGN_OF_A_FIT_THAT_PASSES_THIRTY: [f64; 16] = [
        1.0, 1.0, //
        1.0, 4.0, //
        1.0, 2.0, //
        1.0, 4.0, //
        1.0, 4.0, //
        1.0, 2.0, //
        1.0, 3.0, //
        1.0, 0.0,
    ];

    /// A variant whose fit settles at an effect past 30 has no answer,
    /// which is the only thing in the suite that holds the 30 from above.
    ///
    /// The mark is `|effect| > 30` and every fixture of the panel and of
    /// this module leaves it room on one side only: the largest effect
    /// either panel answers is 2.50, so a mark of 1e9, or of infinity,
    /// would mark the same variants there, the ones that separate the two
    /// groups being caught by the singular factorization a few rounds
    /// later instead. This fixture is on the other side of it. Measured on
    /// 24 September 2026 with the mark taken out: the fit settles at an
    /// effect of 36.4485 on Accelerate and 36.4557 on faer, with a
    /// standard error of 2.0234e7 and a p-value of 0.9999986 on both,
    /// which is an effect a study knows nothing about beside an error of
    /// twenty million, and is what the 30 is there to keep out of a
    /// result. With the mark in, the variant gets the three NaNs at the
    /// round its effect passes 30.
    ///
    /// It is a fit that settles and not one that walks away, so the
    /// factorization never refuses it and the pivot rule of **Open 5** of
    /// `docs/specs/gwas.md` does not catch it either: the 30 is the only
    /// mark of the five that fires here.
    #[test]
    fn a_variant_whose_fit_settles_past_thirty_has_no_answer() {
        let mut vcf = String::from(THE_HEADER_OF_EIGHT);
        vcf.push_str("1\t1000\tv0\tA\tT\t.\t.\t.\tGT");
        for genotype in ["0/1", "1/1", "0/0", "0/1", "0/1", "0/1", "0/1", "0/0"] {
            vcf.push('\t');
            vcf.push_str(genotype);
        }
        vcf.push('\n');
        let study = GwasInput {
            phenotype: &THE_TRAIT_OF_A_FIT_THAT_PASSES_THIRTY,
            trait_type: TraitType::Binomial,
            design: &THE_DESIGN_OF_A_FIT_THAT_PASSES_THIRTY,
            num_coefs: 2,
            kinship: None,
            test: Some(TestType::Wald),
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS_OF_EIGHT,
            transform_to_biallelic: false,
        };
        let mut reader = reader_over(vcf.as_bytes());
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of a fit that settles past 30: {error}"),
        };
        assert_eq!(result.num_vars, 1);
        assert!(
            result.beta[0].is_nan() && result.se[0].is_nan() && result.p_value[0].is_nan(),
            "the variant whose effect passed 30 was answered with a beta of {beta}, an se \
             of {se} and a p-value of {p_value}",
            beta = result.beta[0],
            se = result.se[0],
            p_value = result.p_value[0]
        );
    }

    /// The variance of an effect that is not a finite number above 0 has
    /// no answer, which the standard error is refused at.
    ///
    /// It is the one mark of the five that no data of this repository
    /// reaches, and the pivot rule of **Open 5** of `docs/specs/gwas.md`
    /// puts it further out of reach, since a factorization whose smallest
    /// pivot has fallen that far is already a runaway. So the fit is built
    /// here with the factorization written into it by hand, two
    /// coefficients and a diagonal that a Cholesky would have accepted,
    /// and the standard error is read off it:
    ///
    /// - a diagonal of 1 and 1e-200, where the variance is 1e400, which
    ///   both backends answer as an infinity while reporting success;
    ///   `docs/specs/linalg.md` has the entry they part company on, where
    ///   one gives an infinity and the other a NaN, and the mark is read
    ///   as a value that is not finite so that both are caught;
    /// - a diagonal of 1 and 1e200, where the variance underflows to 0 and
    ///   the standard error would be 0, which would make the statistic
    ///   infinite and the p-value 0, a variant reported as the surest
    ///   finding of the study.
    ///
    /// Both give `None`, which the Wald test answers with the three NaNs
    /// of a variant that has none.
    #[test]
    fn a_variance_that_is_not_a_finite_number_above_zero_has_no_answer() {
        for (diagonal, what) in [
            (1e-200_f64, "a variance that overflows"),
            (1e200_f64, "a variance that underflows to 0"),
        ] {
            let mut fitted = TheFitOfOneVariant {
                phenotype: vec![0.0, 1.0, 0.0, 1.0],
                design: vec![1.0, 0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 3.0],
                weighted_design: vec![0.0; 8],
                coefs: vec![0.0, 0.0],
                step: vec![0.0, 0.0],
                linear_predictor: vec![0.0; 4],
                weights: vec![0.0; 4],
                residuals: vec![0.0; 4],
                // The lower half of the factorization of a system of two
                // coefficients, row after row, which every solve of this
                // crate reads and no factorization of this fixture wrote.
                system: vec![1.0, 0.0, 0.0, diagonal],
                of_the_effect: vec![0.0, 0.0],
                num_individuals: 4,
                num_coefs: 2,
            };
            match fitted.the_error_of_the_effect() {
                Ok(None) => {}
                Ok(Some(se)) => panic!("{what} was answered with a standard error of {se}"),
                Err(error) => panic!("{what}: {error}"),
            }
        }
    }

    /// How far the effect of one of the six variants may be from plink2's,
    /// as a share of the standard error plink2 printed for that variant:
    /// 1e-5.
    ///
    /// It is the bound of "How it is verified" of "The logistic model" of
    /// `docs/specs/gwas.md`, read as that spec item's own rule has every
    /// effect read: a share of the `se` of the variant and never of the
    /// effect itself, because a study is mostly null and an effect that
    /// cancelled to near 0 is no guide to its own error.
    ///
    /// It is not lowered until it fails, as a bound on popnei's own
    /// arithmetic would be. What it measures is where plink2's fit
    /// stopped: plink2 stops its logistic fit earlier than popnei does, so
    /// the two settle at different coefficients, and that distance is
    /// orders of magnitude above both the printing of six significant
    /// digits, which rounds a value by 5e-6 of itself, and the distance
    /// between the two linear algebra backends. The worst of the six is in
    /// the doc comment of the test.
    const OF_PLINK2: f64 = 1e-5;

    /// How far the standard error of one of the six may be from plink2's,
    /// as a share of the standard error itself: 1e-4, which is the spec's
    /// bound and the same rule, the scale of what is estimated.
    const OF_PLINK2S_ERROR: f64 = 1e-4;

    /// How far the p-value of one of the six may be from plink2's, as a
    /// share of it: 5e-3, which is the spec's bound.
    ///
    /// It is 500 times the 1e-5 the linear model is held to against the
    /// same program, and "How it is verified" of "The logistic model" says
    /// why: the two programs stop their logistic fits at different places,
    /// and a run that came within 1e-5 of plink2 on these p-values would
    /// have stopped early too.
    const OF_PLINK2S_P_VALUE: f64 = 5e-3;

    /// What plink2 v2.0.0-a.7.7 wrote for six variants of the panel with
    /// every genotype called, in
    /// `tests/reference/gwas/plink2.panel_called.glm.logistic.hybrid.tsv`,
    /// with `cov1` and `cov2` as covariates: the id, the effect, its
    /// standard error `LOG(OR)_SE` and the p-value `P`, which are the
    /// three columns "How it is verified" of "The logistic model" of
    /// `docs/specs/gwas.md` compares. Five of them are the causal variants
    /// of `causal_vars.csv` and `var0000` is not causal.
    ///
    /// plink2 reports the odds ratio and popnei its logarithm, so the
    /// effect here is the logarithm of the `OR` column and carries more
    /// digits than plink2 printed. The two libraries count the dosages of
    /// a variant from the same allele at all six, plink2 from its `A1` and
    /// popnei from the one that is not the major one among the tested
    /// individuals, so the signs agree; at a variant whose `A1_FREQ` is
    /// above 0.5 they would not, and the effect would be the same number
    /// the other way round.
    const OF_PLINK2_SIX: [(&str, f64, f64, f64); 6] = [
        ("var0000", -0.572_578_694_541_525_8, 0.261_917, 0.028_808_1),
        (
            "var0052",
            -0.852_823_096_430_041_7,
            0.248_979,
            0.000_614_166,
        ),
        ("var0629", -0.949_570_924_908_968, 0.323_553, 0.003_337_39),
        ("var0751", -0.265_916_666_783_602_6, 0.219_207, 0.225_098),
        ("var1137", -0.427_448_481_503_581_95, 0.252_402, 0.090_356_1),
        ("var1188", -0.830_184_139_078_324, 0.264_071, 0.001_667_71),
    ];

    /// The study of the panel with every genotype called against its
    /// binomial trait and the two covariates, with the test the caller
    /// names, which is what both tests of the panel run.
    fn the_logistic_study_of_the_panel(test: Option<TestType>) -> Gwas {
        let path = the_panel_path();
        let options = VcfOptions {
            ploidy: 2,
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let individuals = reader.individuals().to_vec();
        let (phenotype, design) =
            the_trait_and_the_design_of_the_panel(&individuals, TraitType::Binomial);
        let tested: Vec<usize> = (0..individuals.len()).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Binomial,
            design: &design,
            num_coefs: 3,
            kinship: None,
            test,
            use_grammar_gamma_approx: false,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the logistic study of the panel: {error}"),
        }
    }

    /// The Wald test of the six variants of the panel is plink2's: the
    /// effect within 1e-5 of the standard error of that variant, the
    /// standard error within 1e-4 of itself and the p-value within 5e-3 of
    /// itself, which is what the spec asks of the whole columns in Python.
    ///
    /// The study asks for no test, so what it gets is the Wald test, which
    /// is the default of a binomial trait with no kinship.
    ///
    /// Measured over the six on 24 September 2026, on Accelerate and on
    /// faer alike to the three digits given: the worst effect is 2.69e-6
    /// of the standard error of its variant, `var0052`, which is 27 per
    /// cent of what is allowed; the worst standard error is 5.30e-5 of
    /// itself, `var1137`, 53 per cent of its bound; and the worst p-value
    /// is 1.88e-4 of itself, `var1137` again, 4 per cent of its bound. The
    /// two backends differ from each other by 4.4e-16 of a standard error
    /// at the worst of the eighteen numbers, ten orders of magnitude below
    /// what separates either of them from plink2.
    #[test]
    fn the_six_variants_of_the_panel_are_plink2s_wald_effect_error_and_p_value() {
        let result = the_logistic_study_of_the_panel(None);
        assert_eq!(result.num_vars, 1200, "the variants of the panel");
        assert_eq!(result.null_model.model, GwasModel::Glm);
        assert_eq!(result.null_model.test, TestType::Wald);
        assert_eq!(result.null_model.num_individuals, 200);
        assert_eq!(
            result.null_model.covariate_effects.len(),
            3,
            "the intercept and the two covariates"
        );
        assert_eq!(
            result.null_model.residual_variance, None,
            "a binomial trait has no residual variance of its own"
        );
        let ids = result.ids.as_deref().expect("the ids of the variants");
        for (id, beta, se, p_value) in OF_PLINK2_SIX {
            let var = match ids.iter().position(|held| held == id) {
                Some(var) => var,
                None => panic!("{id} is not a variant of the panel"),
            };
            // The effect and the standard error are both measured against
            // the standard error, the scale of what the study estimates,
            // which is what "How it is verified" of "What every model
            // shares" of `docs/specs/gwas.md` asks of every comparison of
            // the two.
            for (found, expected, allowed, what) in [
                (result.beta[var], beta, OF_PLINK2, "the effect"),
                (result.se[var], se, OF_PLINK2S_ERROR, "the standard error"),
            ] {
                let difference = (found - expected).abs();
                assert!(
                    difference <= allowed * se,
                    "{what} of {id} is {found} and plink2 gives {expected}, {difference} \
                     away, which is {share} of the {se} it is uncertain by against the \
                     {allowed} allowed",
                    share = difference / se
                );
            }
            let share = (result.p_value[var] / p_value - 1.0).abs();
            assert!(
                share <= OF_PLINK2S_P_VALUE,
                "the p-value of {id} is {found} and plink2 gives {p_value}, {share} of it \
                 away against the {OF_PLINK2S_P_VALUE} allowed",
                found = result.p_value[var]
            );
        }
    }

    /// The variants of the panel that popnei gives no answer to are exactly
    /// the ones plink2 fell back to a Firth penalized regression for, which
    /// is `var0006` and none of the other 1199.
    ///
    /// A variant that separates the individuals that have the condition
    /// from the ones that have not has no finite effect. popnei gives it
    /// the three NaNs of a variant that has no answer, which is what pyNei
    /// does and what **Open 1** of `docs/specs/gwas.md` leaves with the
    /// owner; plink2 adds a term to the likelihood that pulls the estimate
    /// back from infinity, reports a finite answer and says so in its
    /// `FIRTH?` column. The two programs therefore agree on which variants
    /// have no ordinary answer, and this test is what says that popnei
    /// marks those and no others: a variant wrongly marked loses its
    /// p-value and nothing else in the suite would show it.
    ///
    /// The set is read out of the reference file and not written here as a
    /// literal, so that a reference made again on another machine is
    /// compared against what it says.
    #[test]
    fn the_variants_the_wald_fit_runs_away_on_are_the_ones_plink2_gave_to_firth() {
        let result = the_logistic_study_of_the_panel(Some(TestType::Wald));
        assert_eq!(result.num_vars, 1200, "the variants of the panel");
        let ids = result.ids.as_deref().expect("the ids of the variants");
        let of_popnei: Vec<&str> = ids
            .iter()
            .zip(&result.p_value)
            .filter(|(_, p_value)| p_value.is_nan())
            .map(|(id, _)| id.as_str())
            .collect();
        let of_plink2 = the_variants_plink2_gave_to_firth();
        assert_eq!(
            of_plink2,
            ["var0006"],
            "the variants plink2 marked `FIRTH?` `Y` in its reference file"
        );
        assert_eq!(
            of_popnei, of_plink2,
            "the variants popnei has no answer for against the ones plink2 gave to Firth"
        );
        // The three columns of such a variant are NaN together, and the
        // frequency of its alleles is there as it is for every variant
        // that has no answer: 3 of its 400 alleles are the one that is not
        // the major one, which is `A1_FREQ` 0.015 of the reference file.
        let var = match ids.iter().position(|held| held == "var0006") {
            Some(var) => var,
            None => panic!("var0006 is not a variant of the panel"),
        };
        assert!(
            result.beta[var].is_nan() && result.se[var].is_nan(),
            "var0006 was answered with a beta of {beta} and an se of {se}",
            beta = result.beta[var],
            se = result.se[var]
        );
        assert!(
            (result.allele_freq[var] - 0.015).abs() <= 1e-12,
            "the frequency of var0006 is {found}",
            found = result.allele_freq[var]
        );
    }

    /// The ids of the variants plink2 fell back to a Firth penalized
    /// regression for, which its `FIRTH?` column holds `Y` for, read from
    /// `tests/reference/gwas/plink2.panel_called.glm.logistic.hybrid.tsv`.
    fn the_variants_plink2_gave_to_firth() -> Vec<String> {
        let path = the_reference_path("plink2.panel_called.glm.logistic.hybrid.tsv");
        let text = match std::fs::read_to_string(&path) {
            Ok(text) => text,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let mut lines = text.lines();
        let header: Vec<&str> = lines
            .next()
            .expect("the header of plink2's logistic file")
            .split('\t')
            .collect();
        let column_of = |name: &str| match header.iter().position(|held| *held == name) {
            Some(at) => at,
            None => panic!("plink2's logistic file has no column {name}"),
        };
        let of_the_id = column_of("ID");
        let of_firth = column_of("FIRTH?");
        let mut gave_to_firth = Vec::new();
        for line in lines.filter(|line| !line.is_empty()) {
            let fields: Vec<&str> = line.split('\t').collect();
            let field = |at: usize| match fields.get(at) {
                Some(field) => *field,
                None => panic!("the line `{line}` of plink2's logistic file has no column {at}"),
            };
            if field(of_firth) == "Y" {
                gave_to_firth.push(field(of_the_id).to_owned());
            }
        }
        gave_to_firth
    }

    /// The variant that separates the eight individuals of the Wald
    /// fixture, allele by allele: the five that have the condition carry no
    /// copy of the allele that is not the major one and the three that have
    /// not carry two of it, so no finite effect fits them.
    const THE_SEPARATING_VARIANT: [&str; 8] =
        ["0/0", "1/1", "0/0", "1/1", "0/0", "1/1", "1/1", "1/1"];

    /// What numpy 2.5.3 answers for the ordinary variant of that fixture,
    /// the dosages `0 0 0 0 2 2 2 2`, fitted and tested by the formulas of
    /// "The logistic model" of `docs/specs/gwas.md`: one logistic
    /// regression with the variant in the design, started at the
    /// coefficients of the null fit and an effect of 0 for the variant,
    /// with an explicit inverse of `d' w d` for the standard error. The
    /// p-value is `erfc` of the square root of half the statistic, which
    /// is the chi square with one degree of freedom of "The two
    /// distributions" of the spec.
    ///
    /// The fit settles in 6 rounds, where the separating variant beside it
    /// runs away: its effect passes 30 at the round 29, at -30.48.
    const OF_NUMPYS_WALD_FIT: (f64, f64, f64) = (
        1.000_147_184_039_909_7,
        1.035_791_531_684_167_7,
        0.334_250_715_562_421_6,
    );

    /// How far the effect and the standard error of that variant may be
    /// from numpy's: 1e-12 times the standard error numpy gives it,
    /// 1.0358, which is the scale of what is estimated and not the value
    /// itself, as the same rule of the same spec item asks.
    ///
    /// This one is not lowered until it fails, and the reason is that it
    /// cannot be: measured on 24 September 2026, popnei's effect and
    /// standard error are numpy's to the bit on both backends, so every
    /// bound above 0 passes and none says anything about the arithmetic.
    /// What the number has to leave room for is another machine: a fit of
    /// six rounds reads `exp` six times, `docs/rust_core.md` and the
    /// `coding` skill say that `exp` is not rounded the same on macOS,
    /// Linux and wasm, and a round that starts a bit away from where this
    /// one did settles a bit away from where this one settled. The fit
    /// stops when its step falls below 1e-8 and one more Newton step
    /// squares that, so what a different `exp` can move is about 1e-16 of
    /// the coefficients, and this is four orders of magnitude above it.
    const OF_NUMPYS_WALD_FIT_BOUND: f64 = 1e-12;

    /// How far its p-value may be from numpy's, in `log10`: 1e-12, a
    /// p-value running over orders of magnitude, and the same bound for
    /// the same reason.
    ///
    /// Measured on 24 September 2026 on both backends: 9.64e-17, which is
    /// the last bits of `erfc` read two ways, popnei's from the `libm`
    /// crate and numpy's from `math.erfc` of CPython 3.14.
    const OF_NUMPYS_WALD_P_VALUE: f64 = 1e-12;

    /// A variant whose Wald fit runs away has no answer, and the ordinary
    /// variant beside it in the same block is tested.
    ///
    /// The first variant of the fixture separates the eight individuals
    /// that have the condition from the ones that have not, so its effect
    /// has no finite value to settle at and the fit walks towards one. The
    /// mark that catches it is the effect passing 30: measured with numpy
    /// 2.5.3 on the same fixture on 24 September 2026, the fit reaches
    /// -30.48 at its round 29, where the two marks that would follow, a
    /// step that is not finite and a fit still moving after 50 rounds, are
    /// still ahead of it.
    ///
    /// The third variant is the dosages the covariate of the fixture is a
    /// tenth of, so the design with it in it has two columns that carry
    /// the same information and `d' w d` is no longer a matrix a Cholesky
    /// factorization accepts. That is the first of the three marks, the
    /// system the factorization refuses as singular, and it is reached at
    /// the round 4 on Accelerate and at the round 1 on faer, the two
    /// backends stopping at different rounds and giving the same three
    /// NaNs. numpy 2.5.3 on the same variant, solving each round with an
    /// LU factorization as pyNei does, is answered by a matrix a Cholesky
    /// refuses and reaches the same three NaNs by the third mark instead,
    /// still moving after its 50 rounds at an effect of -9.44. It is the
    /// variant the score test of the fixture above has no answer for
    /// either, by the threshold of **Open 2**.
    ///
    /// The second variant is what says that a runaway is not the block:
    /// it is fitted and tested, and its three numbers are numpy's. The
    /// study is run twice, once asking for the Wald test and once asking
    /// for no test, since the Wald test is what a binomial trait with no
    /// kinship takes when the user asks for none, and both give the same
    /// three rows.
    #[test]
    fn a_variant_whose_wald_fit_runs_away_has_no_answer_and_its_neighbour_is_tested() {
        let mut vcf = String::from(THE_HEADER_OF_EIGHT);
        for (var, genotypes) in [
            THE_SEPARATING_VARIANT,
            ["0/0", "0/0", "0/0", "0/0", "1/1", "1/1", "1/1", "1/1"],
            // the dosages of the covariate, allele by allele
            ["0/1", "0/0", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1"],
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
        for asked in [Some(TestType::Wald), None] {
            let study = GwasInput {
                phenotype: &THE_TRAIT_OF_EIGHT,
                trait_type: TraitType::Binomial,
                design: &THE_DESIGN_OF_THE_FIRST_VARIANT,
                num_coefs: 2,
                kinship: None,
                test: asked,
                use_grammar_gamma_approx: false,
                individuals: &THE_INDIVIDUALS_OF_EIGHT,
                transform_to_biallelic: false,
            };
            let mut reader = reader_over(vcf.as_bytes());
            let result = match the_study_of(&mut reader, &study) {
                Ok(result) => result,
                Err(error) => panic!("the Wald test asked for as {asked:?}: {error}"),
            };
            assert_eq!(result.num_vars, 3);
            assert_eq!(result.null_model.test, TestType::Wald);
            assert!(
                result.beta[0].is_nan() && result.se[0].is_nan() && result.p_value[0].is_nan(),
                "the variant that separates the two groups was answered with a beta of \
                 {beta}, an se of {se} and a p-value of {p_value}",
                beta = result.beta[0],
                se = result.se[0],
                p_value = result.p_value[0]
            );
            // It is in the result with its frequency, as every variant
            // that has no answer is: 6 of its 16 alleles are the one that
            // is not the major one.
            assert!(
                (result.allele_freq[0] - 0.375).abs() <= 1e-12,
                "the frequency of v0 is {found}",
                found = result.allele_freq[0]
            );
            let (beta, se, p_value) = OF_NUMPYS_WALD_FIT;
            for (found, expected, what) in [
                (result.beta[1], beta, "the effect of v1"),
                (result.se[1], se, "the standard error of v1"),
            ] {
                let difference = (found - expected).abs();
                assert!(
                    difference <= OF_NUMPYS_WALD_FIT_BOUND * se,
                    "{what} is {found} and numpy gives {expected}, {difference} away, which \
                     is {share} of the {se} it is uncertain by against the \
                     {OF_NUMPYS_WALD_FIT_BOUND} allowed",
                    share = difference / se
                );
            }
            let in_log10 = (result.p_value[1] / p_value).log10().abs();
            assert!(
                in_log10 <= OF_NUMPYS_WALD_P_VALUE,
                "the p-value of v1 is {found} and numpy gives {p_value}, {in_log10} away in \
                 log10 against the {OF_NUMPYS_WALD_P_VALUE} allowed",
                found = result.p_value[1]
            );
            assert!(
                result.beta[2].is_nan() && result.se[2].is_nan() && result.p_value[2].is_nan(),
                "the variant that is the covariate was answered with a beta of {beta}, an \
                 se of {se} and a p-value of {p_value}",
                beta = result.beta[2],
                se = result.se[2],
                p_value = result.p_value[2]
            );
        }
    }

    /// The header of a VCF of six individuals, which the fixture of a fit
    /// that runs all its rounds needs and which the eight of the others
    /// cannot give: the fits of eight individuals this module has all end
    /// at the factorization before the rounds run out.
    const THE_HEADER_OF_SIX: &str = "##fileformat=VCFv4.2\n\
        ##contig=<ID=1>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti0\ti1\ti2\ti3\ti4\ti5\n";

    /// A null fit still moving after the 50 rounds it is given is refused
    /// by the round count, which is the other of the two ways such a fit
    /// ends and the one pyNei meets.
    ///
    /// The covariate is `0 1 2 3 4 7` and the one individual with the
    /// condition is the one at 7, so the covariate separates the two
    /// groups and no finite effect fits it. What makes this fixture end at
    /// the rounds and not at the factorization is the gap: the individual
    /// with the condition is 3 above the next, so the fit reaches a
    /// predictor that separates them with a small slope and the weights
    /// are still above 0 when the last round runs out. Measured on 24
    /// September 2026 on both backends, it is refused at 50 rounds, where
    /// the eight individuals of
    /// `a_null_fit_that_has_not_settled_is_refused`, whose covariate is 0
    /// to 7, are refused at 45 on Accelerate and 43 on faer by the
    /// factorization. The rounds are asserted and not bounded, because the
    /// count is what says which of the two ways this fixture took: over
    /// the 254 traits of eight individuals with a covariate of 0 to 7, and
    /// over that covariate scaled by 39 factors from an eighth to 5, no
    /// fixture reaches the rounds at all.
    #[test]
    fn a_null_fit_still_moving_after_its_rounds_is_refused_by_their_count() {
        let mut vcf = String::from(THE_HEADER_OF_SIX);
        vcf.push_str("1\t1000\tv0\tA\tT\t.\t.\t.\tGT");
        for genotype in ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1"] {
            vcf.push('\t');
            vcf.push_str(genotype);
        }
        vcf.push('\n');
        let separating = [0.0_f64, 0.0, 0.0, 0.0, 0.0, 1.0];
        let design: [f64; 12] = [
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 2.0, //
            1.0, 3.0, //
            1.0, 4.0, //
            1.0, 7.0,
        ];
        let individuals: [usize; 6] = [0, 1, 2, 3, 4, 5];
        let study = GwasInput {
            phenotype: &separating,
            trait_type: TraitType::Binomial,
            design: &design,
            num_coefs: 2,
            kinship: None,
            test: Some(TestType::Score),
            use_grammar_gamma_approx: false,
            individuals: &individuals,
            transform_to_biallelic: false,
        };
        let mut reader = reader_over(vcf.as_bytes());
        match the_study_of(&mut reader, &study) {
            Err(Error::GwasFitDidNotSettle { model, rounds }) => {
                assert_eq!(model, GwasModel::Glm, "the model that was being fitted");
                assert_eq!(
                    rounds, 50,
                    "the fit is refused by the count of its rounds and not by its \
                     factorization"
                );
            }
            Err(error) => panic!("a covariate that separates the two groups: {error}"),
            Ok(result) => panic!("a study of {} variants was run", result.num_vars),
        }
    }

    /// How many variants a study reads in one block, which is
    /// [`crate::block::MAX_NUM_VARS_PER_BLOCK`]: the pass puts a `Reblock`
    /// before it, so a source has to pass it for the pass to read a second
    /// block.
    const VARS_OF_ONE_BLOCK: usize = 10_000;

    /// A logistic study over more variants than one block holds gives the
    /// same answers in its second block as in its first, with both of its
    /// tests.
    ///
    /// The variants are three patterns over and over, and each is one of
    /// the three answers this model gives: the dosages of the covariate,
    /// which neither test has an answer for, by the threshold of
    /// **Open 2** in the score test and by the singular factorization in
    /// the Wald test; the dosages `0 0 0 0 2 2 2 2`, which both tests
    /// answer and which numpy 2.5.3 gives the numbers of the two fixtures
    /// above for; and a variant every individual is heterozygous at, which
    /// has no variance and so no answer wherever it is. Every variant is
    /// held to the literal of its pattern and not to what the first block
    /// answered, so a study that answered the same wrong thing in every
    /// block would fail this as well.
    ///
    /// What it covers that the fixtures of two and three variants do not
    /// is the buffers of the model and of the dosages being reused from
    /// one block to the next, the answers of the second block being added
    /// after the first's and not over them, and the Wald test's one fit
    /// per variant starting from the null's coefficients again at every
    /// variant of every block. It is the twin of
    /// `a_study_of_more_variants_than_one_block_answers_the_same_in_every_block`
    /// of `linear` and of the one of `linear_mixed`.
    #[test]
    fn a_logistic_study_of_more_variants_than_one_block_answers_the_same_in_every_block() {
        let patterns = [
            // the dosages of the covariate, allele by allele
            ["0/1", "0/0", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1"],
            ["0/0", "0/0", "0/0", "0/0", "1/1", "1/1", "1/1", "1/1"],
            ["0/1", "0/1", "0/1", "0/1", "0/1", "0/1", "0/1", "0/1"],
        ];
        let num_vars = VARS_OF_ONE_BLOCK.saturating_add(101);
        let mut vcf = String::from(THE_HEADER_OF_EIGHT);
        for var in 0..num_vars {
            let pos = var.saturating_add(1).saturating_mul(10);
            // The variants of the second block are on the second
            // chromosome, and the first block fills a whole block.
            let chrom = match var < VARS_OF_ONE_BLOCK {
                true => 1,
                false => 2,
            };
            vcf.push_str(&format!("{chrom}\t{pos}\tv{var}\tA\tT\t.\t.\t.\tGT"));
            for genotype in patterns[var % patterns.len()] {
                vcf.push('\t');
                vcf.push_str(genotype);
            }
            vcf.push('\n');
        }
        for (test, answered) in [
            (TestType::Score, OF_THE_ORDINARY_VARIANT),
            (TestType::Wald, OF_NUMPYS_WALD_FIT),
        ] {
            let study = GwasInput {
                phenotype: &THE_TRAIT_OF_EIGHT,
                trait_type: TraitType::Binomial,
                design: &THE_DESIGN_OF_THE_FIRST_VARIANT,
                num_coefs: 2,
                kinship: None,
                test: Some(test),
                use_grammar_gamma_approx: false,
                individuals: &THE_INDIVIDUALS_OF_EIGHT,
                transform_to_biallelic: false,
            };
            let mut reader = reader_over(vcf.as_bytes());
            let result = match the_study_of(&mut reader, &study) {
                Ok(result) => result,
                Err(error) => panic!("the {test:?} test of {num_vars} variants: {error}"),
            };
            assert_eq!(result.num_vars, num_vars, "the variants of the study");
            assert!(
                num_vars > VARS_OF_ONE_BLOCK,
                "the study has to read more than one block"
            );
            let ids = result.ids.as_deref().expect("the ids of the variants");
            assert_eq!(ids.len(), num_vars, "one id for each variant");
            for var in [0, VARS_OF_ONE_BLOCK, num_vars.saturating_sub(1)] {
                assert_eq!(ids[var], format!("v{var}"), "the id of the variant {var}");
            }
            // The bound of the fixture of two variants for the score test
            // and of the one of three for the Wald test, which are where
            // these three literals come from.
            let allowed = match test {
                TestType::Score => OF_NUMPY,
                TestType::Wald => OF_NUMPYS_WALD_FIT_BOUND,
            };
            let (beta, se, p_value) = answered;
            for var in 0..num_vars {
                match var % patterns.len() {
                    1 => {
                        for (found, expected, what) in [
                            (result.beta[var], beta, "the effect"),
                            (result.se[var], se, "the standard error"),
                        ] {
                            let difference = (found - expected).abs();
                            assert!(
                                difference <= allowed * se,
                                "{what} of the variant {var} of the {test:?} test is \
                                 {found} and numpy gives {expected}, {difference} away, \
                                 which is {share} of the {se} it is uncertain by against \
                                 the {allowed} allowed",
                                share = difference / se
                            );
                        }
                        let in_log10 = (result.p_value[var] / p_value).log10().abs();
                        assert!(
                            in_log10 <= OF_NUMPYS_P_VALUE,
                            "the p-value of the variant {var} of the {test:?} test is \
                             {found} and numpy gives {p_value}, {in_log10} away in log10 \
                             against the {OF_NUMPYS_P_VALUE} allowed",
                            found = result.p_value[var]
                        );
                    }
                    _ => assert!(
                        result.beta[var].is_nan()
                            && result.se[var].is_nan()
                            && result.p_value[var].is_nan(),
                        "the variant {var} of the {test:?} test, which has no answer, was \
                         answered with a beta of {found}, an se of {se} and a p-value of \
                         {p_value}",
                        found = result.beta[var],
                        se = result.se[var],
                        p_value = result.p_value[var]
                    ),
                }
            }
        }
    }

    /// A covariate whose own effect is far past 30 marks no variant, which
    /// is what says the mark is read on the effect of the variant alone.
    ///
    /// `_wald_test` of `pynei/gwas.py` reads the last coefficient of the
    /// fit, which is the variant's, and "The logistic model" of
    /// `docs/specs/gwas.md` says that a covariate whose effect is larger
    /// than 30 marks nothing. Nothing held that: no fixture of this module
    /// and no covariate of either panel has an effect anywhere near it, so
    /// a fit that read the largest coefficient instead of the last would
    /// have passed every test.
    ///
    /// The first covariate of the panel is given here in units a ten
    /// thousandth of its own, which leaves the fit and every answer of the
    /// study where they were and multiplies that covariate's effect by
    /// 1e4. Measured on 24 September 2026 on both backends: its effect in
    /// the null model goes from 0.595018 to 5950.18, and the variants with
    /// no answer are `var0006` and no other, as they are unscaled, with
    /// every effect of the 1199 within 9.437e-16 of the unscaled run.
    #[test]
    fn a_covariate_whose_effect_passes_thirty_marks_no_variant() {
        let path = the_panel_path();
        let options = VcfOptions {
            ploidy: 2,
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let individuals = reader.individuals().to_vec();
        let (phenotype, design) =
            the_trait_and_the_design_of_the_panel(&individuals, TraitType::Binomial);
        // The intercept, the first covariate in units a ten thousandth of
        // its own, and the second as it is.
        let scaled: Vec<f64> = design
            .as_chunks::<3>()
            .0
            .iter()
            .flat_map(|row| [row[0], row[1] * 1e-4, row[2]])
            .collect();
        let tested: Vec<usize> = (0..individuals.len()).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Binomial,
            design: &scaled,
            num_coefs: 3,
            kinship: None,
            test: Some(TestType::Wald),
            use_grammar_gamma_approx: false,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => {
                panic!("the logistic study of the panel with a scaled covariate: {error}")
            }
        };
        let effect_of_the_covariate = result.null_model.covariate_effects[1];
        assert!(
            (effect_of_the_covariate - 5950.18).abs() <= 1e-2,
            "the effect of the scaled covariate in the null model is \
             {effect_of_the_covariate} and it has to be the 5950.18 that puts it far past \
             the 30 of a variant that has run away"
        );
        let ids = result.ids.as_deref().expect("the ids of the variants");
        let with_no_answer: Vec<&str> = ids
            .iter()
            .zip(&result.p_value)
            .filter(|(_, p_value)| p_value.is_nan())
            .map(|(id, _)| id.as_str())
            .collect();
        assert_eq!(
            with_no_answer,
            ["var0006"],
            "the variants with no answer when a covariate's effect is 5950.18"
        );
        // The effects of the study are where they were, which is what
        // says that the scaling moved the covariate's coefficient and
        // nothing else. The bound is the one of the six variants against
        // plink2, read as a share of the standard error of the variant.
        let unscaled = the_logistic_study_of_the_panel(Some(TestType::Wald));
        for (var, (found, expected)) in result.beta.iter().zip(&unscaled.beta).enumerate() {
            if expected.is_nan() {
                continue;
            }
            let difference = (found - expected).abs();
            assert!(
                difference <= OF_PLINK2 * unscaled.se[var],
                "the effect of {id} is {found} where the study with the covariate in its \
                 own units gives {expected}, {difference} away against the {OF_PLINK2} of \
                 its {se} allowed",
                id = ids[var],
                se = unscaled.se[var]
            );
        }
    }
}
