//! The logistic model: a binomial trait and no kinship.
//!
//! [`LogisticModel`] is the null model fitted by iteratively reweighted
//! least squares and the score test of every variant against it. "The
//! logistic model" of `docs/specs/gwas.md` says what it fits and how it is
//! verified. Its Wald test, which fits one logistic regression per
//! variant, is being written and a study that asks for it is refused.

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
/// buffers of a block are kept from one block to the next, so a pass over
/// a million variants asks the machine for them once and allocates nothing
/// for a variant.
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
    /// variants that have variance x `num_coefs`.
    of_the_design: Vec<f64>,
    /// Those rows solved against the factorization, which is
    /// `(d' w d)⁻¹ d' w x` of each variant, of the same size.
    solved: Vec<f64>,
    /// Each variant's dosages times the trait's residuals, one per variant
    /// that has variance.
    num: Vec<f64>,
    /// The effect of each variant that has variance.
    beta: Vec<f64>,
    /// The standard error of each of those effects.
    se: Vec<f64>,
    /// The p-value of each of those tests.
    p_value: Vec<f64>,
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
    /// # Errors
    ///
    /// [`Error::GwasInputOfAnotherSize`] when `phenotype` does not hold
    /// one value for each tested individual.
    /// [`Error::GwasFitDidNotSettle`] when the coefficients are still
    /// moving after 50 rounds, which is what a covariate that separates
    /// the individuals that have the condition from the ones that have not
    /// gives. [`Error::GwasLinalg`] when one of the two products, the
    /// factorization of the weighted design or the solve against it could
    /// not be done, which is where a design the weights have taken to a
    /// matrix that is no longer positive definite is refused.
    pub(crate) fn of_the_study(phenotype: &[f64], design: &Design<'_>) -> Result<LogisticModel> {
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
            reason = "the design holds `num_individuals` times `num_coefs` values and                       `num_individuals` is `num_coefs` plus 2 at least, so this product                       is smaller than the length of the design and fits in a `usize`"
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
            of_the_design: Vec::new(),
            solved: Vec::new(),
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
            for (coef, step) in fitted.coefs.iter_mut().zip(&step) {
                *coef += *step;
            }
            // A step that is not a number is not a step below the
            // tolerance: the comparison is false for it, the round is not
            // the last, and the fit ends with the error below.
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
    /// factored, which is the fit running away: the columns of the design
    /// are independent, `Design::of_the_study` having refused them
    /// otherwise, so the only thing that makes `d' w d` singular is the
    /// weights, and a weight falls to 0 when the chance the fit gives an
    /// individual reaches 0 or 1. A covariate that separates the
    /// individuals that have the condition from the ones that have not is
    /// what takes it there. [`Error::GwasLinalg`] when the product of the
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
            // The logistic curve of the linear predictor, which is the
            // chance that the individual has the condition. A predictor
            // far below 0 overflows the exponential and gives a chance of
            // 0, and one far above it gives 1; both are chances, and the
            // weight of either is 0.
            let chance = 1.0 / (1.0 + (-predicted).exp());
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
    /// The score test is the one this makes. The Wald test of this model
    /// fits one logistic regression per variant and is being written; the
    /// pass refuses a study that asks for it before the null model is
    /// fitted, and this is what says so if one ever reaches here.
    ///
    /// # Errors
    ///
    /// [`Error::GwasModelNotBuilt`] when the Wald test is asked for, and
    /// whatever [`LogisticModel::score_test_the_block`] fails with.
    pub(crate) fn test_the_block(
        &mut self,
        dosages: &GwasDosages,
        test: TestType,
    ) -> Result<Answers<'_>> {
        match test {
            TestType::Score => self.score_test_the_block(dosages),
            TestType::Wald => Err(Error::GwasModelNotBuilt {
                model: GwasModel::Glm,
            }),
        }
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
    /// A variant of which the design leaves at most the tested individuals
    /// times 2.2e-16 of `x' w x` has no answer, and gets the three NaNs a
    /// variant with no variance gets. `den` is 0 or above in exact
    /// arithmetic, and for a variant that is a combination of the columns
    /// of the design it is the rounding of a cancellation, which can fall
    /// below 0: `beta` would be a number divided by noise, `se` the square
    /// root of a negative number and the statistic negative. It is
    /// **Open 2** of `docs/specs/gwas.md`.
    ///
    /// # Errors
    ///
    /// [`Error::GwasVariantsTooLarge`] when the values of the block are
    /// more than a `usize` counts, and [`Error::GwasLinalg`] when one of
    /// the two products or the solve could not be done, which is where a
    /// block of other individuals than the null model was fitted over is
    /// refused.
    fn score_test_the_block(&mut self, dosages: &GwasDosages) -> Result<Answers<'_>> {
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
        // The solve overwrites what it is given, and both what it was
        // given and what it gave are read below, so the rows go into a
        // buffer of their own, which is kept from one block to the next
        // like every other buffer here.
        self.solved.clear();
        self.solved.extend_from_slice(&self.of_the_design);
        popnei_linalg::solve_with_cholesky(
            &self.factored,
            self.num_coefs,
            &mut self.solved,
            num_vars,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "solve of a block of variants against the weighted design",
            source,
        })?;
        // The share of what the variant weighed that the design has to
        // leave of it for the variant to be worth testing.
        let share_that_is_nothing = the_share_that_is_nothing(self.num_individuals);
        for (((of_the_design, solved), num), of_the_variant) in self
            .of_the_design
            .chunks_exact(self.num_coefs)
            .zip(self.solved.chunks_exact(self.num_coefs))
            .zip(&self.num)
            .zip(dosages.dosages().chunks_exact(self.num_individuals))
        {
            let weighted_length = of_the_variant
                .iter()
                .zip(&self.weights)
                .map(|(dosage, weight)| weight * dosage * dosage)
                .sum::<f64>();
            let of_the_covariates = of_the_design
                .iter()
                .zip(solved)
                .map(|(row, solved)| row * solved)
                .sum::<f64>();
            let den = weighted_length - of_the_covariates;
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

/// The logistic model against R's score test on the panel with every
/// genotype called, and the three cases of its own: a variant the design
/// leaves nothing of, a fit that does not settle, and the Wald test that
/// is not written yet.
#[cfg(test)]
mod glm {
    use crate::block::BlockReader;
    use crate::error::Error;
    use crate::gwas::linear::lm::{
        THE_HEADER_OF_EIGHT, reader_over, the_panel_path, the_study_of,
        the_trait_and_the_design_of_the_panel,
    };
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
            test: Some(TestType::Score),
            use_grammar_gamma_approx: false,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the score test of the panel: {error}"),
        };
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
    /// arithmetic. At a covariate equal to the dosages the denominator of
    /// the first variant's score test cancels to exactly 0 on both
    /// backends, and a variant with nothing left is then caught by any
    /// comparison with 0; at a tenth of them the cancellation leaves
    /// 4.4e-16, a positive number, which is the case the threshold of
    /// **Open 2** is for.
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
    /// first reproduction of this rule at either score test: the
    /// denominator of the variant's score test comes to 4.44e-16 on
    /// Accelerate and to exactly 0 on faer, where `x' w x`, the weighted
    /// squared length the variant had before the covariates were taken
    /// out, is 2.3592 and the threshold that refuses it, the eight tested
    /// individuals times 2.2e-16 of that, is 4.19e-15.
    ///
    /// What the variant is answered with when the threshold is taken out,
    /// measured the same day on both backends, is why it is there. The
    /// numerator is 0 to the bit, the design being what the fit made its
    /// residuals at right angles to, so on Accelerate the row is a `beta`
    /// of 0, an `se` of 4.75e7 and a p-value of 1, which a user reads as a
    /// variant that was tested and showed nothing and which no filter on a
    /// missing effect takes out; and on faer, where the denominator is 0,
    /// it is a `beta` of NaN, an `se` of infinity and a p-value of NaN, a
    /// fourth kind of row that `docs/specs/gwas.md` does not describe. The
    /// same fixture with an explicit inverse of `d' w d` in numpy 2.5.3
    /// lands on -4.44e-16 and gives a third row, a `beta` of -0.0 beside
    /// an `se` of NaN.
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
    /// on 24 September 2026, at the round 45 on Accelerate and on faer
    /// alike. The test asserts the rounds are between 1 and the 50 the fit
    /// is given and not that they are 45, since which round the weights
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

    /// The Wald test of a logistic model, which fits one logistic
    /// regression per variant and is being written, is refused, and it is
    /// the test a study of a binomial trait with no kinship gets when it
    /// asks for none.
    #[test]
    fn the_wald_test_of_a_logistic_model_is_refused() {
        let mut vcf = String::from(THE_HEADER_OF_EIGHT);
        vcf.push_str("1\t1000\tv0\tA\tT\t.\t.\t.\tGT");
        for genotype in ["0/0", "0/0", "0/0", "0/0", "1/1", "1/1", "1/1", "1/1"] {
            vcf.push('\t');
            vcf.push_str(genotype);
        }
        vcf.push('\n');
        for asked in [None, Some(TestType::Wald)] {
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
            match the_study_of(&mut reader, &study) {
                Err(Error::GwasModelNotBuilt { model }) => {
                    assert_eq!(model, GwasModel::Glm, "the model of the study");
                }
                Err(error) => panic!("the Wald test asked for as {asked:?}: {error}"),
                Ok(_) => panic!("the Wald test asked for as {asked:?} was run"),
            }
        }
    }
}
