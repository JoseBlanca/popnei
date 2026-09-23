//! The association study of a dataset: which of its variants go with a
//! trait of the individuals, with the effect of each variant on the trait,
//! the uncertainty of that effect, and the chance of an effect at least
//! that large in a dataset where the variant has none, which is its
//! p-value. `docs/specs/gwas.md` says what a study gives and how each of
//! its models is verified.
//!
//! [`calc_gwas`] is the study itself: it fits the null model once, over
//! the individuals that are tested and the design they were given, and
//! then reads the variants in one pass, giving one row for each of them.
//! Of the four models only the linear one is written, a continuous trait
//! with no kinship, and the other three are refused until they are.
//!
//! What is here besides is the two functions that turn the statistic of a
//! test into a p-value, which every model of the module ends in, and the
//! design every model is fitted on. Each of the two is a survival function, the
//! chance that a distribution is beyond a point, which `sf` in their names
//! stands for: `chi2_sf_1df` is the chi square with one degree of freedom
//! that every score test and the Wald test of a logistic model need, and
//! `t_sf_two_sided` the Student t of the Wald test of a linear model and of
//! a linear mixed model, which is written from the regularized incomplete
//! beta function below it.
//!
//! [`Design::of_the_study`] is what a study is refused for before any model
//! is fitted or any variant is read: an individual to test that the source
//! has not, one that is there twice, individuals that are not in the order
//! the source has them, a study of no more individuals than the columns of
//! its design plus one, a phenotype that is not a number or that does not
//! fit the trait, and a design whose columns are not independent. The
//! [`the_model_and_the_test`] chooses which of the four models a study
//! fits and which of the two tests it makes of every variant, and refuses
//! the two pairs of trait, kinship and test that have no test between
//! them. [`GwasDosages`] turns one block into a dosage per tested
//! individual per variant, computed over those individuals and not over
//! the whole panel, with the frequency of each variant among them and
//! whether it has any variance there. [`Gwas`] is what a study gives back,
//! one row per variant, which [`Gwas::add_the_block`] fills block by block
//! with the answers of the variants that have variance and the three NaNs
//! of the ones that have none. [`LinearModel`] is the one model that gives
//! those answers today: the thin QR of the design, the residuals of the
//! trait and the t test of every variant against them. The linear mixed
//! model and the two logistic ones are being written.

use std::fmt;
use std::num::NonZeroUsize;

use popnei_linalg::{TheFirstOperand, TheHalfThatHoldsTheMatrix, TheSecondOperand, ThinQr};

use crate::block::{Block, BlockReader, Reblock};
use crate::error::{Error, Result};
use crate::variant::{
    AlleleCounts, ChromTable, MAX_PLOIDY_OF_THE_VARIANTS, MISSING_CODE, Needs, count_alleles,
    the_codes_of_the_genotypes, the_major_allele, the_row_of,
};

/// The chance that a chi square with one degree of freedom is above `x`.
///
/// This is the p-value of a score test and of the Wald test of a logistic
/// model: the statistic of either is an effect divided by its standard
/// error and then squared, which is a chi square with one degree of freedom
/// when the variant has no effect on the trait. `x` is such a statistic.
///
/// An `x` of 0 or below gives 1.0, as scipy's `chi2.sf` does: nothing
/// exceeds such a statistic, so the chance of exceeding it is 1. That is an
/// answer and not a refusal, and it is written out because the square root
/// of a negative number is NaN and the line below would otherwise give one.
/// A NaN `x` is a different thing, a statistic that was never computed, and
/// it gives NaN. Which caller may arrive with a statistic below 0 is
/// **Open 2** of `docs/specs/gwas.md`.
///
/// The value is `erfc(sqrt(x / 2))`, the complementary error function of
/// the `libm` crate, which gives how much of a standard normal distribution
/// lies further from 0 than a point. It agrees with scipy 1.18.1's
/// `chi2.sf(x, 1)` within 1e-12 of itself over the sample the test of this
/// module uses, and at `x = 100` it gives 1.5239706e-23, so a variant with
/// a strong effect gets the right tail.
#[must_use]
pub fn chi2_sf_1df(x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    libm::erfc((x / 2.0).sqrt())
}

/// The chance that a Student t with `df` degrees of freedom is further
/// from 0 than `t`, both tails added.
///
/// This is the p-value of the Wald test of a linear model and of a linear
/// mixed model: the statistic of either is an effect divided by its
/// standard error, which follows a Student t when the variant has no
/// effect on the trait, and `df` is how many individuals are left once the
/// columns of the design and the variant have been fitted. `t` of either
/// sign gives the same value, and a NaN gives NaN.
///
/// A `df` of 0 or below gives NaN, and so does a NaN `df`. A model with no
/// individual left over has no answer to give, which is not what a
/// statistic of 0 means, and 0 or fewer degrees of freedom would otherwise
/// take `x` to 0 or below and come back as 0.0, the smallest p-value there
/// is. Nothing calls this with such a `df` today, and what keeps it at 1 or
/// above is [`Error::GwasTooFewIndividuals`], which refuses a study of no
/// more individuals than the columns of its design plus one: the degrees
/// of freedom of the Wald test of a linear model are the individuals less
/// those columns less one. So this is a guard and not the repair of a live
/// wrong number.
///
/// The value is the regularized incomplete beta function
/// `I_x(df / 2, 1 / 2)` at `x = df / (df + t * t)`, which is written
/// below. `1 - x`, which the function needs as well, is `t * t / (df + t * t)`
/// here and not `1.0 - x`: for a `t` below about `sqrt(df * eps)`, 5e-07 at
/// 197 degrees of freedom, `x` rounds to 1.0 and `1.0 - x` has lost every
/// digit of a difference that is still the whole of the answer. It agrees
/// with scipy 1.18.1's `2 * t.sf(abs(t), df)` within 1e-10 of itself over
/// the fifteen values of `t` and the five degrees of freedom the test of
/// this module uses, from the 1.6e-96 of `t = 40` at 197 degrees of freedom
/// to the 1 less 8e-08 of `t = 1e-07`.
#[must_use]
pub fn t_sf_two_sided(t: f64, df: f64) -> f64 {
    if df <= 0.0 || df.is_nan() {
        return f64::NAN;
    }
    let t_squared = t * t;
    let x = df / (df + t_squared);
    let one_minus_x = t_squared / (df + t_squared);
    regularized_incomplete_beta(df / 2.0, 0.5, x, one_minus_x)
}

/// The regularized incomplete beta function `I_x(a, b)`, the share of the
/// beta distribution with the shapes `a` and `b` that lies below `x`.
///
/// It is 0 for an `x` at or below 0 and 1 for a `one_minus_x` at or below
/// 0. In between it is the continued fraction of `beta_continued_fraction`,
/// multiplied by the front factor and divided by `a`, while `x` is below
/// `(a + 1) / (a + b + 2)`. The front factor is the exponential of
/// `lgamma(a + b) - lgamma(a) - lgamma(b) + a * ln(x) + b * ln(1 - x)`. At
/// or above `(a + 1) / (a + b + 2)` the fraction converges slowly, and the
/// value is taken from the symmetry `I_x(a, b) = 1 - I_{1-x}(b, a)`
/// instead, which is where `1 - x` does most of its work.
///
/// The caller passes `1 - x` as `one_minus_x` rather than leaving the
/// function to take `1.0 - x`, because a caller that can compute it without
/// subtracting has digits that the subtraction cannot recover: `x` carries
/// 16 digits of itself, and when it is near 1 the difference from 1 keeps
/// only what is left of them. `t_sf_two_sided` is such a caller. A caller
/// with nothing better passes `1.0 - x`, and then the spec's `ln_1p(-x)`
/// and this `ln(one_minus_x)` differ by under 1e-17 in the exponent, which
/// is the size the spec gives for that difference.
///
/// Neither `a` nor `b` may be 0 or negative, and no caller of the module
/// passes one: `a` is half the degrees of freedom of a test and `b` is
/// 1 / 2.
fn regularized_incomplete_beta(a: f64, b: f64, x: f64, one_minus_x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if one_minus_x <= 0.0 {
        return 1.0;
    }
    let log_front =
        libm::lgamma(a + b) - libm::lgamma(a) - libm::lgamma(b) + a * x.ln() + b * one_minus_x.ln();
    let front = log_front.exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * beta_continued_fraction(a, b, x) / a
    } else {
        1.0 - front * beta_continued_fraction(b, a, one_minus_x) / b
    }
}

/// The smallest a denominator of the continued fraction may be: one that
/// has come out at 0 is replaced by this before it divides, which is what
/// lets Lentz's method carry on past a term that vanishes instead of
/// giving an infinity.
///
/// No test covers the five lines that read it, and none can, as "The two
/// distributions" of `docs/specs/gwas.md` sets out: with `b` at 1 / 2,
/// which is what every caller here passes, the first denominator is
/// bounded below by `2 / (a + b + 2)` in both branches, so reaching 1e-300
/// needs `a + b` above about 2e300, a panel of 4e300 individuals. The spec
/// measured the smallest `|c|` or `|d|` over 6009003 calls as
/// 4.027585806198886e-6, which is that bound at a million degrees of
/// freedom. It stays because Numerical Recipes and pyNei have it and
/// because a caller with some other `b` would need it. A reader who finds
/// these lines uncovered should not go looking for the argument that
/// reaches them: there is none.
const CONTINUED_FRACTION_TINY: f64 = 1e-300;

/// How near 1 a round's factor has to come for the fraction to stop. Both
/// this and the 500 rounds below are Numerical Recipes', through pyNei. It
/// caps the work and does not get the digits: with it at 0, so that every
/// call runs its 500 rounds, the spec's sweep moved by 2.3e-13 relative at
/// worst and nothing became not finite.
const CONTINUED_FRACTION_EPS: f64 = 1e-15;

/// How many rounds the fraction may take before it gives what it has.
/// Running out is not an error, and no case of the module's tests comes
/// near it: the 115 calls of the two tests below all stopped, the slowest
/// after 47 rounds, measured on 23 September 2026.
const CONTINUED_FRACTION_MAX_ROUNDS: u32 = 500;

/// The continued fraction of the incomplete beta with the shapes `a` and
/// `b` at `x`, the one of Numerical Recipes evaluated by Lentz's method,
/// which builds the fraction from its front rather than from its far end
/// and so stops as soon as a round no longer changes the value.
fn beta_continued_fraction(a: f64, b: f64, x: f64) -> f64 {
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < CONTINUED_FRACTION_TINY {
        d = CONTINUED_FRACTION_TINY;
    }
    d = 1.0 / d;
    let mut h = d;
    for round in 1..=CONTINUED_FRACTION_MAX_ROUNDS {
        let m = f64::from(round);
        let m2 = 2.0 * m;
        let even = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + even * d;
        if d.abs() < CONTINUED_FRACTION_TINY {
            d = CONTINUED_FRACTION_TINY;
        }
        c = 1.0 + even / c;
        if c.abs() < CONTINUED_FRACTION_TINY {
            c = CONTINUED_FRACTION_TINY;
        }
        d = 1.0 / d;
        h *= d * c;
        let odd = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + odd * d;
        if d.abs() < CONTINUED_FRACTION_TINY {
            d = CONTINUED_FRACTION_TINY;
        }
        c = 1.0 + odd / c;
        if c.abs() < CONTINUED_FRACTION_TINY {
            c = CONTINUED_FRACTION_TINY;
        }
        d = 1.0 / d;
        let delta = d * c;
        h *= delta;
        if (delta - 1.0).abs() < CONTINUED_FRACTION_EPS {
            break;
        }
    }
    h
}

/// What a user measured on each individual, which with the kinship decides
/// which of the four models of `docs/specs/gwas.md` a study fits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraitType {
    /// A measurement, one number per individual, such as the height of a
    /// plant. Without a kinship it is fitted by a linear model and with one
    /// by a linear mixed model.
    Continuous,
    /// 0 or 1: an individual that has a condition and one that has not.
    /// Without a kinship it is fitted by a logistic regression and with one
    /// by a logistic mixed model.
    Binomial,
}

impl TraitType {
    /// The name a user writes for each of the two traits, in the order of
    /// the variants.
    ///
    /// They are here and not in a binding crate so that both of them read
    /// one list: a trait renamed in the core is renamed in Python and in
    /// TypeScript with it, and a third trait would reach the two packages
    /// together. It is what `PopDistMeasure::NAMES` of
    /// [`crate::pop_dists`] does for the measures a user asks for by name.
    pub const NAMES: [&'static str; 2] = ["continuous", "binomial"];

    /// The name a user writes for this trait.
    #[must_use]
    pub fn name(self) -> &'static str {
        let of_the_two = match self {
            TraitType::Continuous => 0,
            TraitType::Binomial => 1,
        };
        // The two names are there, one for each variant of the enum.
        TraitType::NAMES.get(of_the_two).copied().unwrap_or("")
    }

    /// The trait a user named.
    ///
    /// # Errors
    ///
    /// [`Error::GwasTraitOfAnUnknownName`] when the name is of neither
    /// trait, with both of them in the message.
    pub fn of_name(name: &str) -> Result<TraitType> {
        match TraitType::NAMES.iter().position(|known| *known == name) {
            Some(0) => Ok(TraitType::Continuous),
            Some(1) => Ok(TraitType::Binomial),
            Some(_) | None => Err(Error::GwasTraitOfAnUnknownName {
                name: name.to_owned(),
            }),
        }
    }
}

/// Which of the two tests a study makes of every variant.
///
/// Both ask whether the effect of the variant on the trait is 0, and under
/// that they have the same distribution in large samples; they differ in
/// what they cost. Which one each model has and which is its default is in
/// "What every model shares" of `docs/specs/gwas.md`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TestType {
    /// The model is fitted again with the variant in it, and the variant's
    /// effect is measured in its own standard errors away from 0. It costs
    /// a fit per variant.
    Wald,
    /// How steeply the fit would improve if the variant's effect were let
    /// off 0, measured at the null model and against how uncertain that
    /// slope is. It costs no fit per variant.
    Score,
}

impl TestType {
    /// The name a user writes for each of the two tests, in the order of
    /// the variants, which is also the `test` they read in the result.
    ///
    /// They are here for the reason [`TraitType::NAMES`] is.
    pub const NAMES: [&'static str; 2] = ["wald", "score"];

    /// The name a user writes for this test.
    #[must_use]
    pub fn name(self) -> &'static str {
        let of_the_two = match self {
            TestType::Wald => 0,
            TestType::Score => 1,
        };
        // The two names are there, one for each variant of the enum.
        TestType::NAMES.get(of_the_two).copied().unwrap_or("")
    }

    /// The test a user named.
    ///
    /// # Errors
    ///
    /// [`Error::GwasTestOfAnUnknownName`] when the name is of neither
    /// test, with both of them in the message. Which tests a model has is
    /// another question, and [`the_model_and_the_test`] is what answers
    /// it: a name that is of a test popnei makes and that this model has
    /// not is refused there and not here.
    pub fn of_name(name: &str) -> Result<TestType> {
        match TestType::NAMES.iter().position(|known| *known == name) {
            Some(0) => Ok(TestType::Wald),
            Some(1) => Ok(TestType::Score),
            Some(_) | None => Err(Error::GwasTestOfAnUnknownName {
                name: name.to_owned(),
            }),
        }
    }
}

/// What a study is given: the trait of the individuals it tests, the design
/// its models are fitted on, and what the user asked for.
///
/// The individuals tested are those that have a phenotype, and `phenotype`,
/// `design` and `individuals` hold them in the order the source has them,
/// one value, one row and one position each. [`Design::of_the_study`]
/// refuses any other order, since the three are read together and an
/// individual's phenotype would otherwise be measured against another
/// individual's genotypes.
#[derive(Debug, Clone, Copy)]
pub struct GwasInput<'a> {
    /// One value per tested individual: the measurement, or 0.0 or 1.0.
    pub phenotype: &'a [f64],
    /// What was measured.
    pub trait_type: TraitType,
    /// `num_individuals` x `num_coefs`, row after row, the intercept first.
    pub design: &'a [f64],
    /// How many columns the design has: the intercept and one for each
    /// covariate.
    pub num_coefs: usize,
    /// `num_individuals` x `num_individuals`, row after row, already cut to
    /// the individuals that are tested and in their order, or `None` for a
    /// model with no random effect.
    pub kinship: Option<&'a [f64]>,
    /// `None` takes the default for the trait and the kinship.
    pub test: Option<TestType>,
    /// Whether the GRAMMAR-Gamma approximation is used, which a mixed model
    /// can take to spend one product per variant instead of a fit.
    pub use_grammar_gamma_approx: bool,
    /// The positions of the tested individuals among those the reader
    /// gives, in the order `phenotype` and `design` have them.
    pub individuals: &'a [usize],
    /// Whether a variant with more than two alleles among its called
    /// genotypes is read with every allele that is not the major one
    /// counting the same. Without it such a variant is an error.
    pub transform_to_biallelic: bool,
}

/// Which of the three buffers a study was given does not hold the study it
/// was given.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GwasInputShape {
    /// The phenotype does not hold one value for each tested individual.
    Phenotype {
        /// How many values the phenotype holds.
        num_values: usize,
        /// How many individuals are tested.
        num_individuals: usize,
    },
    /// The design does not hold one row of `num_coefs` values for each
    /// tested individual.
    Design {
        /// How many values the design holds.
        num_values: usize,
        /// How many individuals are tested.
        num_individuals: usize,
        /// How many columns the design was said to have.
        num_coefs: usize,
    },
    /// The design has no column, not even the one of ones that fits the
    /// intercept.
    NoCoef,
}

impl fmt::Display for GwasInputShape {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::Phenotype {
                num_values,
                num_individuals,
            } => write!(
                formatter,
                "the phenotype holds {num_values} values and {num_individuals} individuals are tested, and it holds the trait of each of them"
            ),
            Self::Design {
                num_values,
                num_individuals,
                num_coefs,
            } => write!(
                formatter,
                "the design holds {num_values} values and it is {num_individuals} individuals x {num_coefs} columns, one row of {num_coefs} values for each of them"
            ),
            Self::NoCoef => write!(
                formatter,
                "the design has no column, and the intercept gives it a column of ones at least"
            ),
        }
    }
}

/// The matrix every model of a study is fitted on, checked: one row per
/// tested individual and one column per number the model fits, the
/// intercept first and then one for each covariate. It carries the
/// individuals those rows belong to, checked with it.
///
/// [`Design::of_the_study`] is the only way to have one, so a model that
/// takes a `Design` is fitted on individuals the source has, each of them
/// once and in the source's order, on a phenotype that holds a number for
/// each of them and fits the trait, and on columns that are independent
/// and hold numbers.
///
/// The individuals are in it, and not passed beside it, because the
/// dosages of a block are read over them: they are the rows of this matrix
/// in the same order, and [`GwasDosages::read_the_block`] takes a `Design`
/// so that no pass can read a block over positions that nothing checked.
/// `Block::retain_individuals` keeps whatever order it is asked for, so
/// positions that do not rise would put one individual's trait against
/// another individual's genotypes, which is what
/// [`Error::GwasIndividualsOutOfOrder`] exists to prevent.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Design<'a> {
    values: &'a [f64],
    num_individuals: usize,
    num_coefs: usize,
    individuals: &'a [usize],
    num_individuals_of_the_source: usize,
    multiallelic: MultiallelicVariants,
}

impl<'a> Design<'a> {
    /// The design of a study, with everything "Which individuals are
    /// tested, and the design" of `docs/specs/gwas.md` refuses about it,
    /// about the individuals it tests and about their phenotype.
    ///
    /// `num_individuals_of_the_source` is how many individuals the reader
    /// the study reads has, which every position of `input.individuals` is
    /// one of. It is known before the first block is read.
    ///
    /// # Errors
    ///
    /// [`Error::GwasIndividualNotInTheDataset`] when a tested individual is
    /// not one the source has, [`Error::GwasIndividualTestedTwice`] when
    /// one of them is there twice, and [`Error::GwasIndividualsOutOfOrder`]
    /// when they are not in the order the source has them, which is also
    /// what a repeat with another individual between its two halves gives.
    /// [`Error::GwasTooFewIndividuals`] when they are no more than the
    /// columns of the design plus one, which would leave nothing to measure
    /// a variant's uncertainty from.
    /// [`Error::GwasPhenotypeNotFinite`] when a phenotype is not a finite
    /// number, [`Error::GwasPhenotypeNotBinomial`] when a binomial trait
    /// holds a value that is neither 0 nor 1, and
    /// [`Error::GwasPhenotypeOfOneValue`] when every tested individual of
    /// such a trait has the same one.
    /// [`Error::GwasDesignValueNotFinite`] when a value of the design is
    /// not a finite number, naming the individual, the column and the
    /// value. [`Error::GwasCovariatesCollinear`] when the columns of the
    /// design are not independent, and [`Error::GwasLinalg`] when the rank
    /// that finds that out could not be taken.
    /// [`Error::GwasInputOfAnotherSize`] when the phenotype or the design
    /// does not hold one value or one row for each tested individual, or
    /// the design has no column.
    pub(crate) fn of_the_study(
        input: &GwasInput<'a>,
        num_individuals_of_the_source: usize,
    ) -> Result<Design<'a>> {
        let num_individuals = input.individuals.len();
        let num_coefs = input.num_coefs;
        if num_coefs == 0 {
            return Err(Error::GwasInputOfAnotherSize {
                problem: GwasInputShape::NoCoef,
            });
        }
        if input.phenotype.len() != num_individuals {
            return Err(Error::GwasInputOfAnotherSize {
                problem: GwasInputShape::Phenotype {
                    num_values: input.phenotype.len(),
                    num_individuals,
                },
            });
        }
        // A design of more values than a `usize` counts is refused here
        // with the rest: no buffer holds that many, and the product is
        // taken with the method that says so rather than left to wrap.
        if num_individuals
            .checked_mul(num_coefs)
            .is_none_or(|values| input.design.len() != values)
        {
            return Err(Error::GwasInputOfAnotherSize {
                problem: GwasInputShape::Design {
                    num_values: input.design.len(),
                    num_individuals,
                    num_coefs,
                },
            });
        }
        // The design fits one number per column and the variant one more,
        // so the individuals are the columns plus two for one of them to be
        // left over, which is what the uncertainty of the variant's effect
        // is measured from. A `num_coefs` whose plus two does not fit in a
        // `usize` is more columns than any design has and is refused too.
        if num_coefs
            .checked_add(2)
            .is_none_or(|fewest| num_individuals < fewest)
        {
            return Err(Error::GwasTooFewIndividuals {
                num_individuals,
                num_coefs,
            });
        }
        refuse_individuals_that_are_not_the_source_in_order(
            input.individuals,
            num_individuals_of_the_source,
        )?;
        refuse_a_phenotype_that_is_not_the_trait(input.phenotype, input.trait_type)?;
        refuse_a_design_value_that_is_not_finite(input.design, num_coefs)?;
        // The columns of the design have to be independent, and the rank is
        // how many of them are, at the tolerance of numpy's `matrix_rank`,
        // so that a design popnei refuses is a design pyNei refuses.
        let rank =
            popnei_linalg::rank(input.design, num_individuals, num_coefs).map_err(|source| {
                Error::GwasLinalg {
                    operation: "rank of the design",
                    source,
                }
            })?;
        if rank < num_coefs {
            return Err(Error::GwasCovariatesCollinear { num_coefs, rank });
        }
        Ok(Design {
            values: input.design,
            num_individuals,
            num_coefs,
            individuals: input.individuals,
            num_individuals_of_the_source,
            multiallelic: MultiallelicVariants::of_the_study(input),
        })
    }

    /// The design itself, `num_individuals` x `num_coefs`, row after row.
    #[must_use]
    pub(crate) fn values(&self) -> &'a [f64] {
        self.values
    }

    /// How many individuals are tested, which is the rows of the design.
    #[must_use]
    pub(crate) fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// How many columns the design has: the intercept and one for each
    /// covariate.
    #[must_use]
    pub(crate) fn num_coefs(&self) -> usize {
        self.num_coefs
    }

    /// The positions of the tested individuals among those the source has,
    /// which rise and each of which is one the source has.
    #[must_use]
    pub(crate) fn individuals(&self) -> &'a [usize] {
        self.individuals
    }

    /// How many individuals the source has, which every block of it holds
    /// the genotypes of.
    #[must_use]
    pub(crate) fn num_individuals_of_the_source(&self) -> usize {
        self.num_individuals_of_the_source
    }
}

/// The values of the design of a study: a finite number in every row of
/// every column.
///
/// The Python and the TypeScript layers refuse a covariate that is missing
/// or is not a number, so what reaches this is a covariate that came out
/// of a user's own arithmetic as an infinity, and a caller of the core
/// crate. Left in, it would reach the rank, which refuses what it is given
/// and not what it produced, and the user would be told of a defect of
/// popnei where they gave a wrong covariate.
///
/// `design` is `num_individuals` x `num_coefs`, row after row, which the
/// caller has checked, so the row of a value is which individual it
/// belongs to and the rest of the division is which column.
fn refuse_a_design_value_that_is_not_finite(design: &[f64], num_coefs: usize) -> Result<()> {
    for (individual, row) in design.chunks(num_coefs.max(1)).enumerate() {
        for (coef, value) in row.iter().copied().enumerate() {
            if !value.is_finite() {
                return Err(Error::GwasDesignValueNotFinite {
                    individual,
                    coef,
                    value,
                });
            }
        }
    }
    Ok(())
}

/// The individuals a study tests: each of them one the source has, each of
/// them once, and in the order the source has them.
///
/// The three are one walk over the positions, since a position that is not
/// above the one before it is either that one again or one the source has
/// earlier. A repeat with another individual between its two halves comes
/// back as the second of the two and not as the first.
fn refuse_individuals_that_are_not_the_source_in_order(
    individuals: &[usize],
    num_individuals_of_the_source: usize,
) -> Result<()> {
    let mut the_one_before: Option<usize> = None;
    for individual in individuals.iter().copied() {
        if individual >= num_individuals_of_the_source {
            return Err(Error::GwasIndividualNotInTheDataset {
                individual,
                num_individuals: num_individuals_of_the_source,
            });
        }
        if let Some(after) = the_one_before {
            if individual == after {
                return Err(Error::GwasIndividualTestedTwice { individual });
            }
            if individual < after {
                return Err(Error::GwasIndividualsOutOfOrder { individual, after });
            }
        }
        the_one_before = Some(individual);
    }
    Ok(())
}

/// The phenotype of the tested individuals: a finite number for each of
/// them, and for a binomial trait 0 or 1 with somebody in each of the two
/// groups.
///
/// A continuous trait of one value is not refused here. "Which individuals
/// are tested, and the design" of `docs/specs/gwas.md` asks for that
/// refusal of a binomial trait alone, where one of the two groups it
/// compares would be empty.
#[expect(
    clippy::float_cmp,
    reason = "a binomial phenotype is the 0.0 and the 1.0 themselves and not a \
              measurement near either, so what is wanted here is the exact \
              comparison and not one within a tolerance; a value of -0.0 is 0.0 \
              by it, which is the answer for an individual without the condition"
)]
fn refuse_a_phenotype_that_is_not_the_trait(
    phenotype: &[f64],
    trait_type: TraitType,
) -> Result<()> {
    for (position, value) in phenotype.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(Error::GwasPhenotypeNotFinite { position, value });
        }
    }
    match trait_type {
        TraitType::Continuous => Ok(()),
        TraitType::Binomial => {
            for (position, value) in phenotype.iter().copied().enumerate() {
                if value != 0.0 && value != 1.0 {
                    return Err(Error::GwasPhenotypeNotBinomial { position, value });
                }
            }
            // Every value is 0 or 1 by here, so counting those that are 1
            // says whether both groups have somebody in them, and the first
            // value is the one they all have when one of the two is empty.
            let with_the_condition = phenotype.iter().filter(|value| **value == 1.0).count();
            match phenotype.first().copied() {
                Some(value) if with_the_condition == 0 || with_the_condition == phenotype.len() => {
                    Err(Error::GwasPhenotypeOfOneValue { value })
                }
                Some(_) | None => Ok(()),
            }
        }
    }
}

/// Which of the four models of `docs/specs/gwas.md` a study fits, which
/// the trait and the kinship together decide.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GwasModel {
    /// A linear model: a continuous trait and no kinship. The trait is a
    /// straight line in the covariates and the variant.
    Lm,
    /// A linear mixed model: a continuous trait with the kinship as the
    /// covariance of a random effect, which is what a panel with families
    /// in it needs.
    Lmm,
    /// A logistic regression: a binomial trait and no kinship.
    Glm,
    /// A logistic mixed model: a binomial trait with the kinship as the
    /// covariance of a random effect.
    Glmm,
}

impl GwasModel {
    /// The name of each of the four models, in the order of the variants,
    /// which is the `model` of the null model a user reads.
    ///
    /// They are here for the reason [`TraitType::NAMES`] is. No user
    /// writes one: the trait and the kinship are what choose the model.
    pub const NAMES: [&'static str; 4] = ["lm", "lmm", "glm", "glmm"];

    /// The name of this model.
    #[must_use]
    pub fn name(self) -> &'static str {
        let of_the_four = match self {
            GwasModel::Lm => 0,
            GwasModel::Lmm => 1,
            GwasModel::Glm => 2,
            GwasModel::Glmm => 3,
        };
        // The four names are there, one for each variant of the enum.
        GwasModel::NAMES.get(of_the_four).copied().unwrap_or("")
    }

    /// What a study of this model is of, as the message of
    /// [`Error::GwasModelNotBuilt`] names it: the trait, the kinship and
    /// the name the literature gives the model of the two.
    #[must_use]
    pub(crate) fn what_it_is_of(self) -> &'static str {
        match self {
            GwasModel::Lm => "a continuous trait with no kinship is a linear model",
            GwasModel::Lmm => "a continuous trait with a kinship is a linear mixed model",
            GwasModel::Glm => "a binomial trait with no kinship is a logistic regression",
            GwasModel::Glmm => "a binomial trait with a kinship is a logistic mixed model",
        }
    }
}

/// Which model a study fits and which test it makes of every variant: the
/// trait and the kinship give the model, and the test is the one that was
/// asked for or the default of that model.
///
/// The default is the Wald test wherever a fit per variant is cheap, a
/// continuous trait or a binomial one without a kinship, and the score
/// test for a binomial trait with a kinship, where a Wald test would mean
/// one mixed model fit for every variant. The linear mixed model and the
/// logistic regression take either test, and the other two have one each.
///
/// # Errors
///
/// [`Error::GwasScoreTestOfALinearModel`] when the score test is asked of
/// a continuous trait with no kinship, whose only test is the t test of
/// the linear model, and [`Error::GwasWaldTestOfALogisticMixedModel`] when
/// the Wald test is asked of a binomial trait with a kinship.
pub(crate) fn the_model_and_the_test(input: &GwasInput<'_>) -> Result<(GwasModel, TestType)> {
    let model = match (input.trait_type, input.kinship.is_some()) {
        (TraitType::Continuous, false) => GwasModel::Lm,
        (TraitType::Continuous, true) => GwasModel::Lmm,
        (TraitType::Binomial, false) => GwasModel::Glm,
        (TraitType::Binomial, true) => GwasModel::Glmm,
    };
    let test = match (model, input.test) {
        (GwasModel::Lm, Some(TestType::Score)) => {
            return Err(Error::GwasScoreTestOfALinearModel);
        }
        (GwasModel::Glmm, Some(TestType::Wald)) => {
            return Err(Error::GwasWaldTestOfALogisticMixedModel);
        }
        (GwasModel::Lm, None | Some(TestType::Wald)) | (GwasModel::Lmm | GwasModel::Glm, None) => {
            TestType::Wald
        }
        (GwasModel::Lmm | GwasModel::Glm, Some(test)) => test,
        (GwasModel::Glmm, None | Some(TestType::Score)) => TestType::Score,
    };
    Ok((model, test))
}

/// Whether a variant with more than two alleles among the called genotypes
/// of the tested individuals is read or refused, which
/// [`GwasInput::transform_to_biallelic`] chooses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MultiallelicVariants {
    /// Such a variant is an error naming its position among the variants
    /// the reader has given.
    Refused,
    /// Every allele that is not the major one counts the same, so the
    /// dosage of a genotype is how many of its alleles are not the major
    /// one, whichever they are.
    Collapsed,
}

impl MultiallelicVariants {
    /// What the study asked for.
    fn of_the_study(input: &GwasInput<'_>) -> MultiallelicVariants {
        match input.transform_to_biallelic {
            true => MultiallelicVariants::Collapsed,
            false => MultiallelicVariants::Refused,
        }
    }
}

/// What reading one variant gives beside the dosages themselves.
#[derive(Debug, Clone, Copy)]
struct RowDosages {
    /// The frequency of the alleles that are not the major one among the
    /// tested individuals, which is the mean dosage over the ploidy.
    allele_freq: f64,
    /// Whether the called genotypes of the tested individuals hold two
    /// dosages at least, which is what a variant needs to be tested.
    has_variance: bool,
}

/// The buffers one thread keeps while it reads the dosages of the rows of
/// a block, so that nothing is allocated for a variant.
struct DosageScratch {
    /// How often each allele was called among the tested individuals,
    /// which gives the major allele of the variant and how many different
    /// alleles it has.
    allele_counts: AlleleCounts,
    /// The code of the genotype of each tested individual: its dosage, 0
    /// to the ploidy, or [`MISSING_CODE`] for a genotype with an allele
    /// missing.
    codes: Vec<u8>,
}

impl DosageScratch {
    /// The buffers of one thread, for the rows of `num_individuals` tested
    /// individuals.
    fn of(num_individuals: usize) -> DosageScratch {
        DosageScratch {
            allele_counts: [0; 128],
            codes: vec![0; num_individuals],
        }
    }
}

/// The dosages of one variant over the tested individuals written into
/// `row`, with the frequency of the alleles that are not the major one and
/// whether the variant has any variance among them.
///
/// `gts` is the genotypes of one variant, `ploidy` alleles for each tested
/// individual and nobody else, and `row` holds one value for each of them.
/// `position` is which variant of those the reader gave this one is, which
/// the error of a variant with more than two alleles names.
///
/// The dosage of a genotype is how many of its alleles are not the major
/// one, and a genotype with any allele missing takes the mean dosage of
/// its variant, so that once the variant is centered it pulls the
/// individual in no direction. Both the major allele and that mean are of
/// the tested individuals alone, as "What it gives" of
/// `docs/specs/gwas.md` asks, and so is the frequency. The dosages are not
/// divided by anything: `beta` of the result is the effect of one more
/// copy of a non major allele, in the units of the trait, and dividing
/// them by the deviation of the variant would give it in deviations
/// instead. That is where this differs from the standardized row of
/// [`crate::variant`], which every variant of a kinship or of a principal
/// component analysis goes through, and it is why the mean comes back
/// here: the study reports it, as `allele_freq`.
///
/// A variant with no called genotype among the tested individuals has a
/// mean of nothing, which is 0, as `_calc_dosages` of pyNei sets it: its
/// dosages are all 0, its frequency is 0 and it has no variance.
///
/// # Errors
///
/// [`Error::VariantPloidyTooLarge`] when `ploidy` is above
/// [`MAX_PLOIDY_OF_THE_VARIANTS`], which the dosages could not be written
/// one to a byte at. [`Error::GtsNotWholeGenotypes`] when `ploidy` is 0 or
/// `gts` does not hold one genotype of it for each value of `row`.
/// [`Error::VariantWithMoreThanTwoAlleles`] when the tested individuals
/// have more than two different alleles among their called genotypes and
/// the study did not ask for those variants to be read. And whatever the
/// counts of the alleles of one variant refuse, which is
/// [`Error::AlleleBelowTheMissingOne`] and
/// [`Error::MoreAllelesThanACountHolds`].
fn the_dosages_of_a_row(
    gts: &[i8],
    ploidy: usize,
    position: usize,
    multiallelic: MultiallelicVariants,
    scratch: &mut DosageScratch,
    row: &mut [f64],
) -> Result<RowDosages> {
    if ploidy > MAX_PLOIDY_OF_THE_VARIANTS {
        return Err(Error::VariantPloidyTooLarge { ploidy });
    }
    let num_individuals = row.len();
    // One genotype of the ploidy for each tested individual. The ploidy of
    // 0 that `NonZeroUsize` refuses is among these: it would be a genotype
    // of no allele for every one of them.
    let of_a_genotype = match NonZeroUsize::new(ploidy) {
        Some(of_a_genotype) if gts.len() == num_individuals.saturating_mul(ploidy) => of_a_genotype,
        _ => {
            return Err(Error::GtsNotWholeGenotypes {
                num_alleles: gts.len(),
                ploidy,
            });
        }
    };
    let DosageScratch {
        allele_counts,
        codes,
    } = scratch;
    // The alleles are counted for the major allele and for how many
    // different alleles the variant has. How many of them were called is
    // not what the mean is taken over: a genotype with one allele called
    // and one missing has a called allele and no dosage.
    count_alleles(gts, allele_counts)?;
    let num_alleles = allele_counts.iter().filter(|count| **count > 0).count();
    match multiallelic {
        MultiallelicVariants::Refused if num_alleles > 2 => {
            return Err(Error::VariantWithMoreThanTwoAlleles {
                position,
                num_alleles,
            });
        }
        MultiallelicVariants::Refused | MultiallelicVariants::Collapsed => {}
    }
    // The buffer of the codes belongs to the thread and is as long as the
    // rows it has read so far, which are all of the tested individuals:
    // this asks the machine for nothing after the first row.
    codes.resize(num_individuals, 0);
    the_codes_of_the_genotypes(gts, of_a_genotype, the_major_allele(allele_counts), codes);
    let called = the_called_dosages(codes);
    let mean = match called.genotypes {
        0 => 0.0,
        genotypes => called.dosages as f64 / genotypes as f64,
    };
    for (value, code) in row.iter_mut().zip(codes.iter().copied()) {
        *value = match code {
            MISSING_CODE => mean,
            dosage => f64::from(dosage),
        };
    }
    Ok(RowDosages {
        allele_freq: mean / ploidy as f64,
        has_variance: called.highest > called.lowest,
    })
}

/// What the called genotypes of one variant hold: how many they are, the
/// sum of their dosages, and the lowest and the highest of those dosages.
#[derive(Debug, Clone, Copy)]
struct CalledDosages {
    /// How many genotypes were called.
    genotypes: u64,
    /// The sum of their dosages.
    dosages: u64,
    /// The lowest dosage among them, and [`u8::MAX`] when none was called.
    lowest: u8,
    /// The highest dosage among them, and 0 when none was called.
    highest: u8,
}

/// What the called genotypes of one variant hold, read from the code of
/// each of the tested individuals.
///
/// A genotype with an allele missing has no dosage: it is counted in none
/// of the four, so the mean is of the called genotypes and the two
/// extremes are theirs. A variant with no called genotype ends with a
/// lowest of [`u8::MAX`] and a highest of 0, which is why a variant has
/// variance at the strict comparison of the two and not at an inequality.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "a dosage is at most the ploidy, 254, and the genotypes of a variant are at \
              most the alleles of it, which the counts of those alleles checked to be a \
              number a u32 holds, so the sum is below 2^40"
)]
fn the_called_dosages(codes: &[u8]) -> CalledDosages {
    let mut called = CalledDosages {
        genotypes: 0,
        dosages: 0,
        lowest: u8::MAX,
        highest: 0,
    };
    for code in codes.iter().copied() {
        let missing = code == MISSING_CODE;
        called.genotypes += u64::from(!missing);
        called.dosages += u64::from(if missing { 0 } else { code });
        called.lowest = called.lowest.min(if missing { u8::MAX } else { code });
        called.highest = called.highest.max(if missing { 0 } else { code });
    }
    called
}

/// The dosages of the rows of a block over the tested individuals written
/// into `dosages`, and what each row gave, in the order of the block.
///
/// The rows are read on the threads of rayon, as section 3 of
/// `docs/architecture.md` asks: no row reads another and each one writes
/// its own values, so neither the dosages nor the frequencies depend on
/// how many threads there are. The threads are those of the pool the
/// caller is running in, and rayon's global pool only when the caller is
/// in none. Each thread keeps the buffers of one row and allocates nothing
/// per variant.
///
/// `gts` holds the rows of the block, `alleles_per_var` alleles each, over
/// the tested individuals alone; `dosages` holds one value for each of
/// those individuals for each of those rows. `first_var` is which variant
/// of those the reader has given the first row of the block is.
///
/// The error is the one of the first row of the block that has one,
/// wherever the threads found it: each row gives its own result and they
/// are read in the order of the block, so a user who reports a file gets
/// the same message every time.
///
/// # Errors
///
/// What reading the dosages of one row refuses, and
/// [`Error::GwasVariantsTooLarge`] when the position of a variant is
/// beyond what a `usize` counts.
#[cfg(not(target_family = "wasm"))]
fn the_dosages_of_the_rows(
    gts: &[i8],
    alleles_per_var: NonZeroUsize,
    num_individuals: NonZeroUsize,
    ploidy: usize,
    multiallelic: MultiallelicVariants,
    first_var: usize,
    dosages: &mut [f64],
) -> Result<Vec<RowDosages>> {
    use rayon::iter::{IndexedParallelIterator, ParallelIterator};
    use rayon::slice::{ParallelSlice, ParallelSliceMut};

    let rows: Vec<Result<RowDosages>> = gts
        .par_chunks_exact(alleles_per_var.get())
        .zip(dosages.par_chunks_exact_mut(num_individuals.get()))
        .enumerate()
        .map_init(
            || DosageScratch::of(num_individuals.get()),
            |scratch, (var, (gts, row))| {
                let position = first_var
                    .checked_add(var)
                    .ok_or(Error::GwasVariantsTooLarge)?;
                the_dosages_of_a_row(gts, ploidy, position, multiallelic, scratch, row)
            },
        )
        .collect();
    rows.into_iter().collect()
}

/// The same rows, read one after another, which is what WebAssembly does:
/// it has no threads.
///
/// # Errors
///
/// The same as the rows read on threads.
#[cfg(target_family = "wasm")]
fn the_dosages_of_the_rows(
    gts: &[i8],
    alleles_per_var: NonZeroUsize,
    num_individuals: NonZeroUsize,
    ploidy: usize,
    multiallelic: MultiallelicVariants,
    first_var: usize,
    dosages: &mut [f64],
) -> Result<Vec<RowDosages>> {
    the_dosages_of_the_rows_one_by_one(
        gts,
        alleles_per_var,
        num_individuals,
        ploidy,
        multiallelic,
        first_var,
        dosages,
    )
}

/// The rows read one after another into the buffers of one row: what
/// WebAssembly does, and what the test that compares the two ways of
/// reading a block calls.
///
/// # Errors
///
/// The same as the rows read on threads.
#[cfg(any(target_family = "wasm", test))]
fn the_dosages_of_the_rows_one_by_one(
    gts: &[i8],
    alleles_per_var: NonZeroUsize,
    num_individuals: NonZeroUsize,
    ploidy: usize,
    multiallelic: MultiallelicVariants,
    first_var: usize,
    dosages: &mut [f64],
) -> Result<Vec<RowDosages>> {
    let mut scratch = DosageScratch::of(num_individuals.get());
    let mut rows = Vec::new();
    for (var, (gts, row)) in gts
        .chunks_exact(alleles_per_var.get())
        .zip(dosages.chunks_exact_mut(num_individuals.get()))
        .enumerate()
    {
        let position = first_var
            .checked_add(var)
            .ok_or(Error::GwasVariantsTooLarge)?;
        rows.push(the_dosages_of_a_row(
            gts,
            ploidy,
            position,
            multiallelic,
            &mut scratch,
            row,
        )?);
    }
    Ok(rows)
}

/// What the pass over the blocks knows about the block it is handing over:
/// how many alleles the genotype of one individual holds, which the reader
/// gives and which every block of a dataset has, and which variant of
/// those the reader has given the first row of this block is.
///
/// The two are a struct and not two arguments because they are both a
/// count of something and a caller that swapped them would compile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BlockOfThePass {
    /// How many alleles the genotype of one individual holds, which the
    /// reader says its source has. A block of another ploidy is refused:
    /// the frequency of a variant is its mean dosage over this number, so
    /// a block read at another one gives frequencies of nothing.
    pub ploidy: usize,
    /// Which variant of those the reader has given the first row of the
    /// block is, counted from 0. It is what the error of a variant with
    /// more than two alleles names, so it is the reader's count and not
    /// the block's.
    pub first_var: usize,
}

/// The dosages of the variants of one block over the individuals a study
/// tests, with the frequency of the alleles that are not the major one of
/// each variant and whether it has any variance, both over those
/// individuals alone.
///
/// A dosage is how many alleles of a genotype are not the major allele of
/// its variant, so it is a whole number from 0 to the ploidy, and a
/// genotype with any allele missing takes the mean dosage of its variant
/// instead, which is 0 for a variant with no called genotype at all. The
/// major allele and that mean are of the tested individuals, as everything
/// else here is. The dosages are not divided by anything: `beta` of the
/// result is the effect of one more copy of a non major allele in the
/// units of the trait, which the deviation of the variant would turn into
/// deviations.
///
/// The buffers are made as long as a block needs and are kept from one
/// block to the next, so a pass over a million variants asks the machine
/// for them once and allocates nothing for a variant.
#[derive(Debug, Clone)]
pub(crate) struct GwasDosages {
    /// How many variants the block held, which is how many rows of the
    /// result it gives.
    num_vars: usize,
    /// How many individuals the study tests, which is the length of one
    /// row of the dosages.
    num_individuals: usize,
    /// How many variants of the block have variance among them.
    num_with_variance: usize,
    /// The dosages, with the rows of the variants that have variance at
    /// the start, in the order of the block.
    dosages: Vec<f64>,
    /// The frequency of each variant of the block, in its order.
    allele_freq: Vec<f64>,
    /// Whether each variant of the block has variance, in its order.
    has_variance: Vec<bool>,
}

impl GwasDosages {
    /// The buffers of a study that has read no block yet.
    #[must_use]
    pub(crate) fn of_a_study() -> GwasDosages {
        GwasDosages {
            num_vars: 0,
            num_individuals: 0,
            num_with_variance: 0,
            dosages: Vec::new(),
            allele_freq: Vec::new(),
            has_variance: Vec::new(),
        }
    }

    /// The dosages of the variants of `block` over the individuals the
    /// study tests, which replace whatever block was read before.
    ///
    /// The individuals that are not tested leave the block before anything
    /// is counted, so the major allele of a variant, the mean its missing
    /// genotypes take, its frequency and whether it has any variance are
    /// all of the tested individuals and of nobody else. That is what
    /// "What it gives" of `docs/specs/gwas.md` asks for, and it matters as
    /// soon as a phenotype leaves one individual out: the block is the
    /// whole panel and the study is of those that have a trait. The block
    /// is cut down in place, so the caller is left with the genotypes of
    /// the tested individuals and a block that is read twice is a block
    /// whose individuals are already gone.
    ///
    /// `design` carries those individuals, checked, and whether a variant
    /// with more than two alleles is read. `of_the_pass` is the shape the
    /// reader says its source has and where the block sits among the
    /// variants it has given.
    ///
    /// Every block is checked against the reader before a row of it is
    /// read, which is what `alleles_per_var_of` of [`crate::stats`] does
    /// for the two passes of that module, with the same errors: a pass
    /// reads the rows of every block as rows of one run over the variants,
    /// so a block of other individuals or of another ploidy is read one
    /// individual at the place of another, and a block of no variant is a
    /// reader that has stopped without saying so.
    ///
    /// # Errors
    ///
    /// Each of these but the last is a defect of the reader that gave the
    /// block or of the caller, and not of the dataset.
    /// [`Error::FieldsNotInTheBlock`] when the block holds no genotypes,
    /// [`Error::ReaderGaveABlockOfNoVariants`] when it holds no variant,
    /// [`Error::BlocksDoNotFitTogether`] when its individuals or its
    /// ploidy are not the reader's, [`Error::BlockWithNoGenotypeOfAVariant`]
    /// when it holds the genotypes of no individual, and
    /// [`Error::BlockArrayOfAnotherSize`] when its arrays are not of the
    /// size it states, which [`Block::check`] finds and which the rows
    /// that came out are counted against again.
    /// [`Error::BlockTooLarge`] when its individuals times its ploidy are
    /// more than a `usize` counts. What
    /// [`Block::retain_individuals`] refuses of the tested individuals,
    /// which is an individual the block has not and one that is there
    /// twice. And what reading the dosages of one row refuses, with
    /// [`Error::GwasVariantsTooLarge`] when the position of a variant is
    /// beyond what a `usize` counts.
    pub(crate) fn read_the_block(
        &mut self,
        block: &mut Block,
        design: &Design<'_>,
        of_the_pass: BlockOfThePass,
    ) -> Result<()> {
        let missing = Needs::GTS.difference(block.fields());
        if !missing.is_empty() {
            return Err(Error::FieldsNotInTheBlock { fields: missing });
        }
        // The rows of a block are cut out of the sizes it states, so those
        // sizes are checked before anything is read.
        block.check()?;
        if block.num_vars == 0 {
            return Err(Error::ReaderGaveABlockOfNoVariants);
        }
        // The individuals and the ploidy are compared with the reader's
        // before the block is cut down to the tested individuals, since
        // that is what changes the first of the two. A block of another
        // ploidy is the one that no other check catches: its rows would be
        // cut at one width and its genotypes read at another, and when the
        // two disagree enough the rows come out as none at all.
        if block.num_individuals != design.num_individuals_of_the_source()
            || block.ploidy != of_the_pass.ploidy
        {
            return Err(Error::BlocksDoNotFitTogether {
                num_individuals: design.num_individuals_of_the_source(),
                ploidy: of_the_pass.ploidy,
                found_num_individuals: block.num_individuals,
                found_ploidy: block.ploidy,
            });
        }
        if block.alleles_per_var()? == 0 {
            return Err(Error::BlockWithNoGenotypeOfAVariant {
                num_individuals: block.num_individuals,
                ploidy: block.ploidy,
            });
        }
        if !design
            .individuals()
            .iter()
            .copied()
            .eq(0..block.num_individuals)
        {
            block.retain_individuals(design.individuals())?;
        }
        let (Some(num_individuals), Some(alleles_per_var)) = (
            NonZeroUsize::new(block.num_individuals),
            NonZeroUsize::new(block.alleles_per_var()?),
        ) else {
            // The block held the genotypes of one individual at least
            // above, and the tested individuals are one at least, since
            // `retain_individuals` refuses none: neither of these is 0.
            return Err(Error::BlockWithNoGenotypeOfAVariant {
                num_individuals: block.num_individuals,
                ploidy: block.ploidy,
            });
        };
        // The block is of the reader's ploidy and its genotypes are its
        // variants times the alleles of one of them, so this division is
        // exact and the buffer holds one value for each individual of each
        // row.
        let num_values = block
            .num_vars
            .checked_mul(num_individuals.get())
            .ok_or(Error::GwasVariantsTooLarge)?;
        self.dosages.resize(num_values, 0.0);
        let rows = the_dosages_of_the_rows(
            &block.gts,
            alleles_per_var,
            num_individuals,
            of_the_pass.ploidy,
            design.multiallelic,
            of_the_pass.first_var,
            &mut self.dosages,
        )?;
        // The rows are cut out of the genotypes at the width of one
        // variant, so a block whose genotypes are not its variants times
        // that width gives fewer rows than it says it holds, and the
        // variants that are left over would go out of the result with no
        // error. Everything above says that cannot happen here; this is
        // what says so of the rows that actually came out.
        if rows.len() != block.num_vars {
            return Err(Error::BlockArrayOfAnotherSize {
                array: "gts",
                found: block.gts.len(),
                expected: block.num_vars.saturating_mul(alleles_per_var.get()),
            });
        }
        self.num_vars = rows.len();
        self.num_individuals = num_individuals.get();
        self.allele_freq.clear();
        self.allele_freq
            .extend(rows.iter().map(|row| row.allele_freq));
        self.has_variance.clear();
        self.has_variance
            .extend(rows.iter().map(|row| row.has_variance));
        self.num_with_variance = self
            .has_variance
            .iter()
            .filter(|has_variance| **has_variance)
            .count();
        // The rows of the variants that have variance are moved to the
        // start of the buffer, so that the block is one matrix of the
        // variants a model can test and the ones that have no answer are
        // not in it. A block with no variant to leave out moves nothing.
        for (to, (var, _)) in self
            .has_variance
            .iter()
            .enumerate()
            .filter(|(_, has_variance)| **has_variance)
            .enumerate()
        {
            if to != var {
                let from = the_row_of(var, self.num_individuals);
                let start = the_row_of(to, self.num_individuals).start;
                self.dosages.copy_within(from, start);
            }
        }
        Ok(())
    }

    /// How many variants the block held, which is how many rows of the
    /// result it gives: the ones that have no answer are among them. It is
    /// 0 before a block has been read.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "what the block held, which the models that report it will read; \
                      the study itself calls `of_a_study`, `read_the_block`, \
                      `num_with_variance` and `dosages`, and until then the tests of \
                      this module are what read this"
        )
    )]
    pub(crate) fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// How many individuals the study tests, which is the length of one
    /// row of [`GwasDosages::dosages`]. It is 0 before a block has been
    /// read, since it is the block that says which of its individuals were
    /// kept.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "what the block held, which the models that report it will read; \
                      the study itself calls `of_a_study`, `read_the_block`, \
                      `num_with_variance` and `dosages`, and until then the tests of \
                      this module are what read this"
        )
    )]
    pub(crate) fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// How many variants of the block have variance among the tested
    /// individuals, which is how many rows [`GwasDosages::dosages`] holds.
    #[must_use]
    pub(crate) fn num_with_variance(&self) -> usize {
        self.num_with_variance
    }

    /// The dosages of the variants that have variance, one row of
    /// [`GwasDosages::num_individuals`] values for each of them, in the
    /// order of the block. It is what a model tests as one matrix, and it
    /// is empty before a block has been read.
    ///
    /// A dosage is a whole number from 0 to the ploidy, how many alleles
    /// of the genotype are not the major allele of its variant among the
    /// tested individuals, and a genotype with any allele missing holds
    /// the mean dosage of its variant instead, which is what centering the
    /// variant would make 0.
    #[must_use]
    pub(crate) fn dosages(&self) -> &[f64] {
        // The buffer holds one row for every variant of the block, and the
        // rows of the variants that have variance were moved to its start,
        // so it holds this many values at least. A buffer that did not
        // would be a defect of this module, and what it gives then is no
        // value at all and not a longer slice, which a model would read as
        // more variants than the block holds.
        let values = self.num_with_variance.saturating_mul(self.num_individuals);
        self.dosages.get(..values).unwrap_or_default()
    }

    /// The frequency of the alleles that are not the major one, over the
    /// tested individuals: one for each variant of the block, in its
    /// order, the variants that have no answer among them. It is the mean
    /// dosage of the variant over the ploidy, and 0 for a variant with no
    /// called genotype among those individuals.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "what the block held, which the models that report it will read; \
                      the study itself calls `of_a_study`, `read_the_block`, \
                      `num_with_variance` and `dosages`, and until then the tests of \
                      this module are what read this"
        )
    )]
    pub(crate) fn allele_freq(&self) -> &[f64] {
        &self.allele_freq
    }

    /// Whether each variant of the block has variance among the tested
    /// individuals, in the order of the block. A variant that has none has
    /// no answer, as "The variants that have no answer" of
    /// `docs/specs/gwas.md` says.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "what the block held, which the models that report it will read; \
                      the study itself calls `of_a_study`, `read_the_block`, \
                      `num_with_variance` and `dosages`, and until then the tests of \
                      this module are what read this"
        )
    )]
    pub(crate) fn has_variance(&self) -> &[bool] {
        &self.has_variance
    }
}

/// Whether a study used the GRAMMAR-Gamma approximation, which a mixed
/// model can take to spend one product per variant instead of a fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "only a mixed model makes the approximation, and the linear model is \
                  the one that is written, so nothing outside the tests builds `Used`"
    )
)]
pub(crate) enum GrammarGammaApprox {
    /// The study made the approximation.
    Used,
    /// The study did not, which is what every model without a kinship
    /// does, since there is nothing to approximate.
    NotUsed,
}

/// Whether a study makes the GRAMMAR-Gamma approximation when the user
/// asks for nothing: it does not. It is the default of `calc_gwas` of
/// `pynei/gwas.py`, which popnei keeps, and the reason is in "What it
/// gives" of the approximation in `docs/specs/gwas.md`: it costs accuracy
/// where a panel is strongly structured, so a user asks for it.
pub const DEFAULT_USE_GRAMMAR_GAMMA_APPROX: bool = false;

/// The model a study fitted without any variant in it, which every variant
/// is then tested against.
#[derive(Debug, Clone, PartialEq)]
pub struct NullModel {
    /// Which of the four models it is.
    pub model: GwasModel,
    /// Which test is made of every variant.
    pub test: TestType,
    /// One per column of the design: the intercept first and then one for
    /// each covariate.
    pub covariate_effects: Vec<f64>,
    /// What the model left unexplained, and `None` for a binomial trait,
    /// whose variance is decided by its mean.
    pub residual_variance: Option<f64>,
    /// The variance of the random effect of the kinship, and `None`
    /// without a kinship.
    pub genetic_variance: Option<f64>,
    /// The genetic variance over the sum of the two, only for the linear
    /// mixed model.
    pub heritability: Option<f64>,
    /// How many individuals the study tests.
    pub num_individuals: usize,
}

/// What a model answered for the variants of one block that have variance,
/// in the order of the block: one value each for the effect of one more
/// copy of a non major allele, how uncertain that effect is, and the
/// p-value of the test that the effect is 0.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Answers<'a> {
    /// The effect of one more copy of a non major allele.
    pub beta: &'a [f64],
    /// The standard error of that effect.
    pub se: &'a [f64],
    /// The p-value of the test that the effect is 0.
    pub p_value: &'a [f64],
}

/// What a study gives back: one row for each variant it was given, and the
/// null model every one of them was tested against.
#[derive(Debug, Clone)]
pub struct Gwas {
    /// How many variants the study was given, which is the rows of every
    /// column below.
    pub num_vars: usize,
    /// The model fitted without any variant in it.
    pub null_model: NullModel,
    /// The frequency of the alleles that are not the major one, over the
    /// tested individuals.
    pub allele_freq: Vec<f64>,
    /// The effect of one more copy of a non major allele, in the units of
    /// the trait for a continuous one and as a log odds ratio for a
    /// binomial one.
    pub beta: Vec<f64>,
    /// The standard error of that effect.
    pub se: Vec<f64>,
    /// The p-value of the test that the effect is 0.
    pub p_value: Vec<f64>,
    /// Whether the GRAMMAR-Gamma approximation was used.
    pub used_grammar_gamma_approx: bool,
    /// The interned chromosome of each variant, read through
    /// [`Gwas::chrom_table`], and `None` when the source had no such
    /// column.
    pub chroms: Option<Vec<u32>>,
    /// The names of the chromosomes the numbers of `chroms` stand for.
    pub chrom_table: ChromTable,
    /// The position of each variant, and `None` when the source had no
    /// such column.
    pub poss: Option<Vec<u64>>,
    /// The id of each variant, and `None` when the source had no such
    /// column.
    pub ids: Option<Vec<String>>,
}

impl Gwas {
    /// The result of a study that has fitted its null model and read no
    /// block yet.
    ///
    /// The null model is fitted before the pass over the blocks, from the
    /// trait, the design and the kinship alone, as "How it runs" of
    /// `docs/specs/gwas.md` says, so the result is opened with it and the
    /// blocks add their rows to it. `chrom_table` is the one of the reader
    /// the study reads, which the numbers of `chroms` are read through.
    #[must_use]
    pub(crate) fn of_the_null_model(
        null_model: NullModel,
        grammar_gamma_approx: GrammarGammaApprox,
        chrom_table: ChromTable,
    ) -> Gwas {
        Gwas {
            num_vars: 0,
            null_model,
            allele_freq: Vec::new(),
            beta: Vec::new(),
            se: Vec::new(),
            p_value: Vec::new(),
            used_grammar_gamma_approx: matches!(grammar_gamma_approx, GrammarGammaApprox::Used),
            chroms: None,
            chrom_table,
            poss: None,
            ids: None,
        }
    }

    /// It adds one row for each variant of a block: the frequency of every
    /// one of them, and the answer of the ones that have variance.
    ///
    /// A variant whose dosages are all the same among the tested
    /// individuals has no variance and cannot be tested. Its row is still
    /// here, with its frequency, and `beta`, `se` and `p_value` are NaN,
    /// which is "The variants that have no answer" of
    /// `docs/specs/gwas.md`: a variant with one allele lands there, and so
    /// does one where every tested individual is heterozygous and one with
    /// no called genotype at all. `answers` holds one value for each
    /// variant that does have variance, in the order of the block, which
    /// is the order [`GwasDosages::dosages`] holds their rows in.
    ///
    /// The three columns `chroms`, `poss` and `ids` are not filled here:
    /// they are the block's own columns and belong to the pass that reads
    /// the blocks.
    ///
    /// # Errors
    ///
    /// [`Error::GwasAnswersOfAnotherSize`] when a column of `answers` does
    /// not hold one value for each variant of the block that has variance,
    /// and [`Error::GwasVariantsTooLarge`] when the variants of the study
    /// are more than a `usize` counts.
    pub(crate) fn add_the_block(
        &mut self,
        dosages: &GwasDosages,
        answers: Answers<'_>,
    ) -> Result<()> {
        for (column, num_values) in [
            ("beta", answers.beta.len()),
            ("se", answers.se.len()),
            ("p_value", answers.p_value.len()),
        ] {
            if num_values != dosages.num_with_variance {
                return Err(Error::GwasAnswersOfAnotherSize {
                    column,
                    num_values,
                    num_with_variance: dosages.num_with_variance,
                });
            }
        }
        self.num_vars = self
            .num_vars
            .checked_add(dosages.num_vars)
            .ok_or(Error::GwasVariantsTooLarge)?;
        let mut answered = answers
            .beta
            .iter()
            .zip(answers.se)
            .zip(answers.p_value)
            .map(|((beta, se), p_value)| (*beta, *se, *p_value));
        for (allele_freq, has_variance) in dosages.allele_freq.iter().zip(&dosages.has_variance) {
            self.allele_freq.push(*allele_freq);
            let answer = match *has_variance {
                false => Some((f64::NAN, f64::NAN, f64::NAN)),
                true => answered.next(),
            };
            // The answers hold one value for each variant of the block that
            // has variance, which was counted above, so there is one here
            // for every variant that has any.
            let Some((beta, se, p_value)) = answer else {
                return Err(Error::GwasAnswersOfAnotherSize {
                    column: "beta",
                    num_values: answers.beta.len(),
                    num_with_variance: dosages.num_with_variance,
                });
            };
            self.beta.push(beta);
            self.se.push(se);
            self.p_value.push(p_value);
        }
        Ok(())
    }
}

/// The linear model of a study fitted without any variant in it, with the
/// buffers one block of variants is tested in.
///
/// The trait is a straight line in the columns of the design, and the fit
/// is the thin QR of that design, `d = q r`: the coefficients are the `c`
/// of `r c = q' y`, the residuals are `y - q q' y`, what the model left
/// unexplained, and the residual sum of squares is their squared length.
/// That is "The linear model" of `docs/specs/gwas.md`, and it is what
/// plink2's `--glm` fits.
///
/// The `q` and the residuals are kept because every variant is tested
/// against them: the covariates are taken out of the dosages of a block in
/// one product, and the effect of a variant is then the plain slope of the
/// trait's residuals on the variant's. The buffers of a block are kept
/// from one block to the next, so a pass over a million variants asks the
/// machine for them once and allocates nothing for a variant.
pub(crate) struct LinearModel {
    /// The `q` of the thin QR of the design, `num_individuals` x
    /// `num_coefs`, row after row: its columns are of length 1 and at
    /// right angles to each other, and they span what the design spans.
    q: Vec<f64>,
    /// The effect of the intercept and of each covariate, one per column
    /// of the design.
    coefs: Vec<f64>,
    /// What the design left of the trait, one value per tested individual,
    /// which every variant is tested against.
    residuals: Vec<f64>,
    /// The squared length of those residuals, which is what the null model
    /// left unexplained.
    rss: f64,
    /// How many individuals the study tests.
    num_individuals: usize,
    /// How many columns the design has.
    num_coefs: usize,
    /// The individuals less the columns of the design, which the residual
    /// sum of squares is divided by for the variance of the null model. It
    /// is 2 at least, since [`Design::of_the_study`] refuses a study of no
    /// more individuals than its columns plus one.
    degrees_of_freedom_of_the_null: usize,
    /// The same less one more for the variant, which is the degrees of
    /// freedom of the t test of every variant. It is 1 at least.
    degrees_of_freedom: usize,
    /// The dosages of a block times the `q`, the variants that have
    /// variance x `num_coefs`.
    of_the_design: Vec<f64>,
    /// The dosages of a block with the design taken out of them, the
    /// variants that have variance x `num_individuals`.
    residualized: Vec<f64>,
    /// Each variant's residuals times the trait's, one per variant that
    /// has variance.
    num: Vec<f64>,
    /// The effect of each variant that has variance.
    beta: Vec<f64>,
    /// The standard error of each of those effects.
    se: Vec<f64>,
    /// The p-value of each of those tests.
    p_value: Vec<f64>,
}

impl LinearModel {
    /// The linear model of `phenotype` over `design`, fitted without any
    /// variant in it.
    ///
    /// `phenotype` holds one value for each individual the design has a
    /// row for, which is what the study tests.
    ///
    /// # Errors
    ///
    /// [`Error::GwasInputOfAnotherSize`] when `phenotype` does not hold
    /// one value for each tested individual, and [`Error::GwasLinalg`]
    /// when the thin QR of the design, the solve against its `r` or one of
    /// the two products could not be done.
    pub(crate) fn of_the_study(phenotype: &[f64], design: &Design<'_>) -> Result<LinearModel> {
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
        let ThinQr { q, r } = popnei_linalg::thin_qr(design.values(), num_individuals, num_coefs)
            .map_err(|source| Error::GwasLinalg {
            operation: "thin QR of the design",
            source,
        })?;
        // `q' y`, which both the coefficients and the residuals are built
        // from: the residuals are `y - q (q' y)` and the coefficients the
        // solution of `r c = q' y`, so it is computed once and the solve
        // that overwrites it comes after the residuals are taken.
        let mut coefs = vec![0.0_f64; num_coefs];
        popnei_linalg::product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &q,
                rows: num_coefs,
            },
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: phenotype,
                cols: 1,
            },
            &mut coefs,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design's q with the trait",
            source,
        })?;
        let mut residuals = vec![0.0_f64; num_individuals];
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &q,
                rows: num_individuals,
            },
            num_coefs,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &coefs,
                cols: 1,
            },
            &mut residuals,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design's q with the fitted trait",
            source,
        })?;
        for (residual, measured) in residuals.iter_mut().zip(phenotype) {
            *residual = measured - *residual;
        }
        let rss = residuals
            .iter()
            .map(|residual| residual * residual)
            .sum::<f64>();
        popnei_linalg::solve_triangular(
            &r,
            num_coefs,
            TheHalfThatHoldsTheMatrix::TheUpperHalf,
            &mut coefs,
            1,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "solve of the design's r against the trait",
            source,
        })?;
        // The two are written with the plain operator and not with
        // `saturating_sub`: saturating is not the meaning wanted here, and
        // a caller that reached this with fewer individuals would get 0
        // degrees of freedom and an infinite standard error where it
        // should get the error that `Design::of_the_study` raises.
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "`Design::of_the_study` refuses a study of no more individuals than \
                      the columns of its design plus one, and a `Design` is the only way \
                      to reach this, so `num_individuals` is `num_coefs` plus 2 at least"
        )]
        let degrees_of_freedom_of_the_null = num_individuals - num_coefs;
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "the line above is 2 at least, by the same refusal"
        )]
        let degrees_of_freedom = degrees_of_freedom_of_the_null - 1;
        Ok(LinearModel {
            q,
            coefs,
            residuals,
            rss,
            num_individuals,
            num_coefs,
            degrees_of_freedom_of_the_null,
            degrees_of_freedom,
            of_the_design: Vec::new(),
            residualized: Vec::new(),
            num: Vec::new(),
            beta: Vec::new(),
            se: Vec::new(),
            p_value: Vec::new(),
        })
    }

    /// What the fit left unexplained: the squared length of the residuals
    /// of the trait, which each variant's own sum of squares is taken from.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the result of a study carries the residual variance, which is this \
                      over the degrees of freedom of the null, and not this itself; the \
                      worked example of `docs/specs/gwas.md` gives both and the test of \
                      it is what reads this"
        )
    )]
    pub(crate) fn residual_sum_of_squares(&self) -> f64 {
        self.rss
    }

    /// The null model of the result: the effects of the intercept and of
    /// the covariates, and what the fit left unexplained over the
    /// individuals less the columns of the design.
    ///
    /// A linear model has no kinship, so it has neither a genetic variance
    /// nor a heritability.
    #[must_use]
    pub(crate) fn null_model(&self, test: TestType) -> NullModel {
        NullModel {
            model: GwasModel::Lm,
            test,
            covariate_effects: self.coefs.clone(),
            residual_variance: Some(self.rss / self.degrees_of_freedom_of_the_null as f64),
            genetic_variance: None,
            heritability: None,
            num_individuals: self.num_individuals,
        }
    }

    /// The t test of every variant of a block that has variance among the
    /// tested individuals, in the order of the block.
    ///
    /// The covariates are taken out of the dosages of the whole block in
    /// one product, `x - (x q) q'`, and after that the effect of a variant
    /// is the plain slope of the trait's residuals on the variant's: `beta`
    /// is `num / xx`, with `num` the variant's residuals times the trait's
    /// and `xx` the squared length of the variant's. What the variant
    /// leaves unexplained is the squared length of the trait's residuals
    /// less `beta` times the variant's, so each variant gets its own
    /// estimate of the residual variance, which is what makes this a t
    /// test and not a normal one, and `se` is the square root of that over
    /// the degrees of freedom and over `xx`. The p-value is
    /// [`t_sf_two_sided`] of `beta / se`.
    ///
    /// A variant of which the design leaves at most the tested individuals
    /// times 2.2e-16 of its squared length has no answer, and gets the
    /// three NaNs a variant with no variance gets: it is a combination of
    /// the columns of the design, so what is left of it is rounding.
    ///
    /// # Errors
    ///
    /// [`Error::GwasVariantsTooLarge`] when the values of the block are
    /// more than a `usize` counts, and [`Error::GwasLinalg`] when one of
    /// the three products could not be done, which is where a block of
    /// other individuals than the null model was fitted over is refused.
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
        let values = num_vars
            .checked_mul(self.num_individuals)
            .ok_or(Error::GwasVariantsTooLarge)?;
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
                values: &self.q,
                cols: self.num_coefs,
            },
            &mut self.of_the_design,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of a block of variants with the design's q",
            source,
        })?;
        self.residualized.resize(values, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &self.of_the_design,
                rows: num_vars,
            },
            self.num_coefs,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &self.q,
                cols: self.num_individuals,
            },
            &mut self.residualized,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of a block of variants with the design",
            source,
        })?;
        // What the product wrote is the part of each variant the design
        // explains, and what is tested is what is left of the variant.
        for (residualized, dosage) in self.residualized.iter_mut().zip(dosages.dosages()) {
            *residualized = dosage - *residualized;
        }
        self.num.resize(num_vars, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &self.residualized,
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
            operation: "product of a block of variants with the trait's residuals",
            source,
        })?;
        let degrees_of_freedom = self.degrees_of_freedom as f64;
        // The share of its own squared length that a variant has to keep
        // once the design is taken out of it to be worth testing, which is
        // the threshold of **Open 2** of `docs/specs/gwas.md`.
        let share_that_is_nothing = self.num_individuals as f64 * f64::EPSILON;
        for ((row, num), dosages) in self
            .residualized
            .chunks_exact(self.num_individuals)
            .zip(&self.num)
            .zip(dosages.dosages().chunks_exact(self.num_individuals))
        {
            let xx = row.iter().map(|value| value * value).sum::<f64>();
            let of_the_dosages = dosages.iter().map(|value| value * value).sum::<f64>();
            // A variant that is a combination of the columns of the design
            // has nothing left once they are taken out, and what `xx` holds
            // is the rounding of that cancellation: `beta` would be a
            // number divided by noise, large and of whichever sign the
            // rounding chose, and the two backends do not choose the same
            // one. Such a variant has no answer, as one with no variance
            // has.
            if xx <= share_that_is_nothing * of_the_dosages {
                self.beta.push(f64::NAN);
                self.se.push(f64::NAN);
                self.p_value.push(f64::NAN);
                continue;
            }
            let beta = num / xx;
            // What the variant leaves unexplained, formed from its own
            // residuals and not taken from the null model's sum of squares
            // by subtracting `beta * num`: those two quantities agree to
            // their last bits once a variant explains most of what the
            // null left, and the subtraction then gives the rounding of a
            // cancelled sum, which is 0 or negative as often as not. "The
            // linear model" of `docs/specs/gwas.md` measures what that
            // gave: an `se` of 0 at one variant and NaN at another, and
            // the two backends disagreeing about which.
            let rss = row
                .iter()
                .zip(&self.residuals)
                .map(|(value, residual)| {
                    let left = residual - beta * value;
                    left * left
                })
                .sum::<f64>();
            let se = (rss / degrees_of_freedom / xx).sqrt();
            self.beta.push(beta);
            self.se.push(se);
            self.p_value
                .push(t_sf_two_sided(beta / se, degrees_of_freedom));
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

/// The association study of the variants of a reader against a trait of
/// the individuals: one row for each variant, with the frequency of the
/// alleles that are not the major one, the effect of one more copy of one
/// of them on the trait, how uncertain that effect is and the p-value of
/// the test that it is 0.
///
/// `input` says which individuals are tested, with their trait and the
/// design the model is fitted on, and "The Rust interface" of
/// `docs/specs/gwas.md` lays out what each of its fields holds. The trait
/// and the kinship choose the model, and of the four only the linear
/// model, a continuous trait with no kinship, is written: the other three
/// are refused with [`Error::GwasModelNotBuilt`] until they are.
///
/// The null model is fitted before the first block is read, from the trait
/// and the design alone, and then one pass over the blocks tests every
/// variant against what it left unexplained. The pass borrows the reader
/// and does not take it, so that whoever built the chain of filters reads
/// their counts when it returns; it asks for the genotypes and for the
/// chromosome, the position and the id of each variant, and puts
/// [`Reblock`] before it, since a filter leaves blocks of uneven size and
/// the test of a block is matrix work.
///
/// `gamma_pass` is the second pass over the same variants that the
/// GRAMMAR-Gamma approximation reads its first block of. Nothing reads it
/// yet: the approximation stands in for the denominator of a mixed model's
/// test, so a study without a kinship is refused for asking for it and a
/// study with one is refused for its model, and no study that reaches the
/// pass has an approximation to make.
///
/// # Errors
///
/// [`Error::GwasModelNotBuilt`] when the study needs one of the three
/// models that are not written, and
/// [`Error::GwasGrammarGammaWithoutAKinship`] when the approximation was
/// asked for by a study with no kinship. What
/// [`the_model_and_the_test`] refuses of the trait, the kinship and the
/// test, and what [`Design::of_the_study`] refuses of the individuals, the
/// phenotype and the design. [`Error::GwasLinalg`] when the fit or the
/// test of a block could not be done. [`Error::PassGaveNoVariant`] when
/// the reader gives no variant, with what each filter of the pass was
/// given and kept, which is pyNei's refusal of a study with nothing to
/// test. [`Error::FieldsNotInTheBlock`] when a block holds no genotypes,
/// or when it has not a column that the first block of the pass had. What
/// [`GwasDosages::read_the_block`] refuses of a block and of its variants,
/// [`Error::GwasVariantsTooLarge`] when the variants of the study are more
/// than a `usize` counts, and whatever the reader fails with.
#[expect(
    unused_variables,
    reason = "`gamma_pass` keeps the name `docs/specs/gwas.md` gives it, since rustdoc \
              prints the names of the arguments; nothing reads it until a mixed model \
              makes the approximation, and the doc comment above says why"
)]
pub fn calc_gwas<R1: BlockReader, R2: BlockReader>(
    reader: &mut R1,
    gamma_pass: Option<&mut R2>,
    input: &GwasInput<'_>,
) -> Result<Gwas> {
    let (model, test) = the_model_and_the_test(input)?;
    if input.use_grammar_gamma_approx && input.kinship.is_none() {
        return Err(Error::GwasGrammarGammaWithoutAKinship);
    }
    match model {
        GwasModel::Lm => {}
        GwasModel::Lmm | GwasModel::Glm | GwasModel::Glmm => {
            return Err(Error::GwasModelNotBuilt { model });
        }
    }
    let ploidy = reader.ploidy();
    let design = Design::of_the_study(input, reader.individuals().len())?;
    let mut fitted = LinearModel::of_the_study(input.phenotype, &design)?;
    let mut result = Gwas::of_the_null_model(
        fitted.null_model(test),
        GrammarGammaApprox::NotUsed,
        ChromTable::new(),
    );
    // The genotypes and the three columns of the result are what this
    // reads, so a reader over a file leaves the other columns of a variant
    // unparsed.
    reader.set_needs(Needs::GTS | Needs::CHROM_POS | Needs::ID);
    let mut blocks = Reblock::new(reader, None)?;
    let mut dosages = GwasDosages::of_a_study();
    let mut first_var = 0_usize;
    while let Some(mut block) = blocks.next_block()? {
        dosages.read_the_block(&mut block, &design, BlockOfThePass { ploidy, first_var })?;
        let answers = fitted.test_the_block(&dosages)?;
        result.add_the_block(&dosages, answers)?;
        the_columns_of_the_block(&mut result, &block, first_var)?;
        first_var = first_var
            .checked_add(block.num_vars)
            .ok_or(Error::GwasVariantsTooLarge)?;
    }
    // The names of the chromosomes are taken when the pass is over and not
    // before it: a reader over a file interns the name of a variant as it
    // reads the variant, so the table of a reader that has read nothing is
    // empty, and the numbers of `chroms` would stand for no name.
    result.chrom_table = blocks.chroms().clone();
    if result.num_vars == 0 {
        let filters = blocks.filtering_stats();
        return Err(Error::PassGaveNoVariant {
            // The filter nearest the source was given what the source
            // gave; with no filter the pass gave what the source gave,
            // which is nothing.
            num_vars_of_the_source: filters.last().map_or(0, |(_, stats)| stats.vars_processed),
            filters,
        });
    }
    Ok(result)
}

/// Adds the chromosome, the position and the id of the variants of a block
/// to the ones the blocks before it gave.
///
/// Which of the three the result has is the first block's: a source that
/// carries none of them gives a result with none, and `docs/specs/gwas.md`
/// says that a VCF and a vars file both carry the chromosome and the
/// position, so in practice only the id is ever missing. A later block
/// that has not one of them is refused, since its variants would take the
/// values of the variants after them.
///
/// A later block that holds a column the first had not is not read, and
/// nothing of the study can give one: the fields are asked for once before
/// the pass, and [`Reblock`] gives a block whose columns differ from the
/// one before it on its own rather than joining the two.
///
/// # Errors
///
/// [`Error::FieldsNotInTheBlock`] when a block has not a column that the
/// first block of the pass had.
fn the_columns_of_the_block(result: &mut Gwas, block: &Block, first_var: usize) -> Result<()> {
    if first_var == 0 {
        result.chroms = block.chrom.as_ref().map(|_| Vec::new());
        result.poss = block.pos.as_ref().map(|_| Vec::new());
        result.ids = block.id.as_ref().map(|_| Vec::new());
    }
    let mut missing = Needs::empty();
    match (result.chroms.as_mut(), block.chrom.as_ref()) {
        (Some(of_the_pass), Some(of_the_block)) => of_the_pass.extend_from_slice(of_the_block),
        (Some(_), None) => missing |= Needs::CHROM_POS,
        (None, _) => {}
    }
    match (result.poss.as_mut(), block.pos.as_ref()) {
        (Some(of_the_pass), Some(of_the_block)) => of_the_pass.extend_from_slice(of_the_block),
        (Some(_), None) => missing |= Needs::CHROM_POS,
        (None, _) => {}
    }
    match (result.ids.as_mut(), block.id.as_ref()) {
        (Some(of_the_pass), Some(of_the_block)) => of_the_pass.extend_from_slice(of_the_block),
        (Some(_), None) => missing |= Needs::ID,
        (None, _) => {}
    }
    if !missing.is_empty() {
        return Err(Error::FieldsNotInTheBlock { fields: missing });
    }
    Ok(())
}

/// The two distributions against scipy 1.18.1, whose numbers are the
/// literals here, as "How it is verified" of "The two distributions" of
/// `docs/specs/gwas.md` asks.
#[cfg(test)]
mod distributions {
    use super::{chi2_sf_1df, regularized_incomplete_beta, t_sf_two_sided};

    /// `x` and `scipy.stats.chi2.sf(x, 1)` of scipy 1.18.1, each printed
    /// by Python's `repr`, which gives the fewest digits that read the same
    /// `f64` back.
    ///
    /// The thirteen values of `x` that come first are a sample of a chi
    /// square with one degree of freedom: the 1000 draws of
    /// `numpy.random.default_rng(0).chisquare(1, 1000)` of numpy 2.5.3,
    /// sorted, at the ranks 0, 90, 180 and so on to 990, and the largest of
    /// the 1000. They run from 9.8e-08, where the p-value is near 1, to
    /// 10.8. The
    /// three that follow are the 30, 50 and 100 the spec asks for, which
    /// reach a p-value of 1.5e-23 and are the tail a variant with a strong
    /// effect lands in.
    const SCIPY_CHI2_SF_1DF: [(f64, f64); 16] = [
        (9.83015230925879e-08, 0.999749838669647),
        (0.012365365489339817, 0.911458017301763),
        (0.05460850890382662, 0.8152298165740756),
        (0.11322041823323797, 0.7365071071670479),
        (0.20872737725861593, 0.6477667063269484),
        (0.3141949722019434, 0.5751173190337298),
        (0.4625430778164763, 0.49643785458596523),
        (0.6940024026584967, 0.404806342457091),
        (1.054858588171547, 0.3043907579517296),
        (1.605831694008574, 0.20507872868262017),
        (2.5243520191320643, 0.11210082147194526),
        (6.274619438315357, 0.012247942227758789),
        (10.821607681554072, 0.0010032233803823583),
        (30.0, 4.320463057827496e-08),
        (50.0, 1.537459794428033e-12),
        (100.0, 1.5239706048320995e-23),
    ];

    /// The tolerance is the 1e-12 relative of the spec, and not the bits:
    /// `erfc` is one of the functions that are not rounded the same on
    /// every platform, so popnei promises the digits and not the last one.
    /// The largest difference over the sixteen values, on this Mac on 23
    /// September 2026, was 1.8e-14 of scipy's value, at `x = 1.6058`, so
    /// the check has 57 times the headroom it needs and a platform that
    /// rounds `erfc` elsewhere still passes it.
    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "an x at or below 0 returns the literal 1.0, which no \
                  arithmetic has touched, so the two are equal to the bit"
    )]
    fn chi2_sf_1df_matches_scipy_over_a_chi_square_sample_and_in_the_tail() {
        for (x, scipy) in SCIPY_CHI2_SF_1DF {
            let popnei = chi2_sf_1df(x);
            let relative = ((popnei - scipy) / scipy).abs();
            assert!(
                relative < 1e-12,
                "chi2_sf_1df({x}) gave {popnei}, scipy gives {scipy}, \
                 which differ by {relative} of scipy's value"
            );
        }
        for x in [-1.0, -1e-9, 0.0] {
            assert_eq!(
                chi2_sf_1df(x),
                1.0,
                "nothing exceeds a statistic of {x}, so the chance of \
                 exceeding it is 1, which is what scipy's chi2.sf gives"
            );
        }
        assert!(
            chi2_sf_1df(f64::NAN).is_nan(),
            "a statistic that is not a number has no p-value"
        );
    }

    /// The four pairs `(a, b)` the spec asks the incomplete beta to be
    /// checked at, in the order the rows of `SCIPY_BETAINC` have them.
    /// `(98.5, 0.5)` is what a Student t with 197 degrees of freedom uses,
    /// next to the panel's 196: 200 individuals less the three columns of
    /// its design less one for the variant.
    const BETA_PAIRS: [(f64, f64); 4] = [(0.5, 0.5), (10.0, 0.5), (98.5, 0.5), (2.5, 7.0)];

    /// The `x` the incomplete beta is taken at: the 1000 draws of
    /// `numpy.random.default_rng(0).uniform(0, 1, 1000)` of numpy 2.5.3,
    /// sorted, at the ranks 0, 111, 222 and so on to 999, so that the last
    /// of them is the largest of the 1000.
    ///
    /// The function has two branches, the continued fraction below
    /// `(a + 1) / (a + b + 2)` and the symmetry `I_x(a, b) =
    /// 1 - I_{1-x}(b, a)` at or above it, and every pair of `BETA_PAIRS`
    /// takes both over these ten values. For `(98.5, 0.5)` the turn is at
    /// 0.985 and the largest of the 1000, 0.9995, is the only one past it.
    const BETA_X: [f64; 10] = [
        0.00019000160734350402,
        0.12681710226124776,
        0.24209440214479672,
        0.36730140225246877,
        0.4752246242920236,
        0.5776878925178592,
        0.6848100430995354,
        0.7892475482526763,
        0.8889806580193464,
        0.9995013522570269,
    ];

    /// `scipy.special.betainc(a, b, x)` of scipy 1.18.1, one row per pair of
    /// `BETA_PAIRS`, each row over `BETA_X`. The one value that is 0.0 and
    /// the one that is 1.0 are what the front factor of the function
    /// underflows to: at `a = 98.5` and `x = 0.00019` it is `exp(-844)`.
    const SCIPY_BETAINC: [[f64; 10]; 4] = [
        // a = 0.5, b = 0.5
        [
            0.008775513005121197,
            0.2317969726600946,
            0.3274906007616421,
            0.41449691362201163,
            0.4842210445079407,
            0.5496588399779172,
            0.6205121607128349,
            0.6963598067259873,
            0.783746069762884,
            0.9857828301556638,
        ],
        // a = 10.0, b = 0.5
        [
            1.080460509271319e-38,
            2.0155738406315239e-10,
            1.3801539354259724e-07,
            9.655661123183734e-06,
            0.0001375787217554091,
            0.0010620947644470852,
            0.0065538478356862995,
            0.03161946390607385,
            0.1297005104596951,
            0.9214264864085051,
        ],
        // a = 98.5, b = 0.5
        [
            0.0,
            2.794490395778511e-90,
            1.36820610510675e-62,
            1.0163053160450192e-44,
            1.1661642215623547e-33,
            2.9187000208624543e-25,
            6.366584355165142e-18,
            9.119063601149856e-12,
            1.5182519264715282e-06,
            0.7542320646178294,
        ],
        // a = 2.5, b = 7.0
        [
            2.4788133598912145e-08,
            0.16372356028276183,
            0.4889434870152894,
            0.7825611783294392,
            0.922059784639868,
            0.978602728785751,
            0.9966059404816893,
            0.9997577325040361,
            0.9999968246187708,
            1.0,
        ],
    ];

    /// The tolerance is the 1e-12 absolute of the spec and, beside it, 1e-12
    /// of scipy's value wherever that value is above 0. The absolute one
    /// alone checks nothing at 8 of the 40 literals, which are at or below
    /// 1e-12 themselves, seven of them in the row of `(98.5, 0.5)`: an
    /// implementation that gave 0 in the tail passed it. The largest
    /// absolute difference at each pair, on this Mac on 23 September 2026,
    /// was 1.3e-15 at `(0.5, 0.5)`, 2.7e-15 at `(10, 0.5)`, 8.5e-15 at
    /// `(98.5, 0.5)` and 3.3e-16 at `(2.5, 7)`, and the largest relative
    /// one was 1.4e-15, 2.1e-14, 6.6e-14 and 1.2e-15 at the same four, so
    /// the relative bound has 15 times the room it needs where it has
    /// least, at `(98.5, 0.5)`, and a platform that rounds `exp` and
    /// `lgamma` elsewhere still passes both.
    ///
    /// The `x` at or below 0 and at or above 1, which a Student t never
    /// reaches because its `x` is `df / (df + t * t)`, are asserted here
    /// too, at the 0 and the 1 the spec gives them, so that the section's
    /// count of tests is the three it names.
    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "an x outside (0, 1) returns the literal 0.0 or 1.0, which \
                  no arithmetic has touched, so the two are equal to the bit"
    )]
    fn incomplete_beta_matches_scipy_at_the_four_pairs_the_two_models_use() {
        for ((a, b), scipy_row) in BETA_PAIRS.into_iter().zip(SCIPY_BETAINC) {
            for (x, scipy) in BETA_X.into_iter().zip(scipy_row) {
                let popnei = regularized_incomplete_beta(a, b, x, 1.0 - x);
                let difference = (popnei - scipy).abs();
                assert!(
                    difference < 1e-12,
                    "the incomplete beta of a = {a}, b = {b} at x = {x} gave \
                     {popnei}, scipy gives {scipy}, which differ by {difference}"
                );
                if scipy > 0.0 {
                    let relative = difference / scipy;
                    assert!(
                        relative < 1e-12,
                        "the incomplete beta of a = {a}, b = {b} at x = {x} \
                         gave {popnei}, scipy gives {scipy}, which differ by \
                         {relative} of scipy's value"
                    );
                }
            }
            for x in [-1.0, 0.0] {
                assert_eq!(regularized_incomplete_beta(a, b, x, 1.0 - x), 0.0);
            }
            for x in [1.0, 2.0] {
                assert_eq!(regularized_incomplete_beta(a, b, x, 1.0 - x), 1.0);
            }
        }
    }

    /// The degrees of freedom of `SCIPY_T_SF_TWO_SIDED`, one per row, the
    /// five "How it is verified" of the spec names. 5 and 17 are the small
    /// ones pyNei checks, where the t is furthest from a normal. 197 is
    /// what the pair `(98.5, 0.5)` of `BETA_PAIRS` comes from, beside the
    /// panel's own 198 with the intercept alone and 196 with its two
    /// covariates, which is 200 individuals less the columns of the design
    /// less one for the variant. 997 and 9997 are the panels of 1000 and
    /// of 10000 individuals of `docs/objectives.md` counted the same way,
    /// and they are here because the error of the function grows with the
    /// degrees of freedom and the spec claims its 1e-10 up to 9997 and no
    /// further.
    const T_DEGREES_OF_FREEDOM: [f64; 5] = [5.0, 17.0, 197.0, 997.0, 9997.0];

    /// The `t` the two sided tail is taken at: the 1000 draws of
    /// `numpy.random.default_rng(0).standard_normal(1000) * 3` of numpy
    /// 2.5.3, sorted, at the ranks 0, 111, 222 and so on to 999, and then
    /// the 1e-07 and 1e-05 of the paragraph below and the 10, 20 and 40 the
    /// spec asks for. Five of the fifteen are negative, which the function
    /// takes through `t * t`, and the three largest are the tail a variant
    /// with a strong effect lands in: at 197 degrees of freedom `t = 40`
    /// has a p-value of 1.6e-96.
    ///
    /// A `t` of 1e-07 and one of 1e-05 are there because `x` of the
    /// incomplete beta is `df / (df + t * t)`, which rounds to 1.0 for a
    /// `t` below about `sqrt(df * eps)`, 5e-07 at 197 degrees of freedom
    /// and 3e-06 at 9997. Their p-value is 1 less something near 8e-08 and
    /// 8e-06, so an implementation that reads `1 - x` back from a rounded
    /// `x` loses the whole of that difference and gives 1.0.
    const T_VALUES: [f64; 15] = [
        -11.698265190163017,
        -3.6455770776863603,
        -2.2826132625192868,
        -1.4177630302754523,
        -0.6135674651989825,
        0.2141653057738903,
        1.0758242029481555,
        2.156821697445163,
        3.4945919433330843,
        9.19811021714669,
        1e-07,
        1e-05,
        10.0,
        20.0,
        40.0,
    ];

    /// `2 * scipy.stats.t.sf(abs(t), df)` of scipy 1.18.1, one row per
    /// degrees of freedom of `T_DEGREES_OF_FREEDOM`, each row over
    /// `T_VALUES`. The last of the row of 9997 is 0.0 because the value
    /// is below the smallest number an `f64` holds.
    const SCIPY_T_SF_TWO_SIDED: [[f64; 15]; 5] = [
        // df = 5
        [
            8.021986064011076e-05,
            0.014816691179735638,
            0.07130468304194079,
            0.21545726847462635,
            0.5663371298234476,
            0.8388781542853476,
            0.3311547120148151,
            0.08352172539200577,
            0.017384749708466502,
            0.00025487505551227214,
            0.999999924078662,
            0.9999924078662037,
            0.00017094757574296357,
            5.7755163732241715e-06,
            1.841196217177295e-07,
        ],
        // df = 17
        [
            1.484092308678517e-09,
            0.0020008264729398156,
            0.03559999374852078,
            0.17432936212171152,
            0.54762394337886,
            0.8329659145993707,
            0.29703772549565155,
            0.045628366419844114,
            0.0027764075279970296,
            5.1955514732118756e-08,
            0.999999921375654,
            0.9999921375654,
            1.5482821011627456e-08,
            2.9927008195203556e-13,
            2.9320117368673818e-18,
        ],
        // df = 197
        [
            2.398338008739183e-24,
            0.0003413422609230366,
            0.023521078785496523,
            0.15783973504362406,
            0.5402089878596277,
            0.8306396600846631,
            0.2833216667123333,
            0.032230891954663506,
            0.00058638199104543,
            5.316767147362756e-17,
            0.9999999203127337,
            0.9999920312733652,
            2.65918277307826e-19,
            2.5792670821731053e-49,
            1.6256956659911944e-96,
        ],
        // df = 997
        [
            1.0285108008842e-29,
            0.0002805744116466835,
            0.022663015344545432,
            0.15657241277885833,
            0.5396411434338867,
            0.8304619422914968,
            0.2822661515237383,
            0.03125803706343817,
            0.0004956644692266718,
            2.0753391089503426e-19,
            0.9999999202315485,
            0.9999920231548545,
            1.678371723944798e-22,
            4.386345719949788e-75,
            1.7484717633928704e-209,
        ],
        // df = 9997
        [
            2.083192798083594e-31,
            0.00026814639056538145,
            0.022474062079316716,
            0.15629112844735527,
            0.5395151482437246,
            0.8304225182279459,
            0.28203184818277116,
            0.031043323571632118,
            0.00047684358853721274,
            4.3710342345251525e-20,
            0.9999999202135392,
            0.9999920213539206,
            1.9634289994060445e-23,
            2.7678195351970412e-87,
            0.0,
        ],
    ];

    /// The tolerance is the 1e-10 relative of the spec, which holds down to
    /// the 1.6e-96 of `t = 40` at 197 degrees of freedom: a p-value is read
    /// at every size, so this one is relative where the incomplete beta's
    /// is absolute. The largest difference at each degrees of freedom, on
    /// this Mac on 23 September 2026, was 1.8e-15 at 5, 8.3e-15 at 17,
    /// 1.9e-13 at 197, 1.6e-13 at 997 and 2.2e-11 at 9997, all of scipy's
    /// value, so the bound has 4.5 times the room it needs at 9997, where
    /// it has least. The error grows with the degrees of freedom and the
    /// spec claims the 1e-10 no further than 9997.
    ///
    /// A value scipy gives as 0.0 is below the smallest number an `f64`
    /// holds, and there a ratio says nothing; the assertion there is that
    /// popnei underflowed too.
    #[test]
    fn t_sf_two_sided_matches_scipy_from_5_to_9997_degrees_of_freedom() {
        for (df, scipy_row) in T_DEGREES_OF_FREEDOM.into_iter().zip(SCIPY_T_SF_TWO_SIDED) {
            for (t, scipy) in T_VALUES.into_iter().zip(scipy_row) {
                let popnei = t_sf_two_sided(t, df);
                if scipy <= 0.0 {
                    assert!(
                        popnei < f64::MIN_POSITIVE,
                        "t_sf_two_sided({t}, {df}) gave {popnei}, where scipy \
                         underflowed to 0"
                    );
                    continue;
                }
                let relative = ((popnei - scipy) / scipy).abs();
                assert!(
                    relative < 1e-10,
                    "t_sf_two_sided({t}, {df}) gave {popnei}, scipy gives \
                     {scipy}, which differ by {relative} of scipy's value"
                );
            }
        }
        assert!(
            t_sf_two_sided(f64::NAN, 197.0).is_nan(),
            "a statistic that is not a number has no p-value"
        );
        for df in [0.0, -4.0, f64::NAN] {
            assert!(
                t_sf_two_sided(1.0, df).is_nan(),
                "{df} degrees of freedom has no p-value, and before this was \
                 asserted `t_sf_two_sided(1.0, 0.0)` gave 0.0"
            );
        }
    }
}

/// The individuals a study tests, the design and what each of them is
/// refused for, as "Which individuals are tested, and the design" of
/// `docs/specs/gwas.md` states them.
#[cfg(test)]
mod design {
    use std::ptr;

    use super::{Design, GwasInput, GwasInputShape, TraitType};
    use crate::error::Error;

    /// A study of a continuous trait over the individuals at `individuals`
    /// of the source, with no kinship and the default test.
    fn a_study<'a>(
        phenotype: &'a [f64],
        design: &'a [f64],
        num_coefs: usize,
        individuals: &'a [usize],
    ) -> GwasInput<'a> {
        GwasInput {
            phenotype,
            trait_type: TraitType::Continuous,
            design,
            num_coefs,
            kinship: None,
            test: None,
            use_grammar_gamma_approx: false,
            individuals,
            transform_to_biallelic: false,
        }
    }

    /// The four tested individuals of the fixtures below: the phenotype of
    /// each, and the design of the intercept and one covariate, four rows
    /// of two values. Every phenotype and every covariate differs from the
    /// others, so a row that went to another individual can be seen.
    const PHENOTYPE_OF_FOUR: [f64; 4] = [1.5, -0.5, 2.0, 0.25];
    const DESIGN_OF_FOUR: [f64; 8] = [1.0, 0.5, 1.0, -1.5, 1.0, 2.5, 1.0, -0.25];

    /// A design of eight rows and two columns whose columns are at right
    /// angles to each other: the ones of the intercept, and a covariate of
    /// `size` and `-size` in turn. Two columns at right angles have their
    /// own lengths for singular values, `sqrt(8)` and `size * sqrt(8)`, so
    /// the smallest of the two is `size` of the largest and the rank sees
    /// exactly the `size` asked for.
    fn a_design_whose_smallest_singular_value_is(size: f64) -> Vec<f64> {
        [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
            .into_iter()
            .flat_map(|sign| [1.0, sign * size])
            .collect()
    }

    /// The individuals of a study are tested in the order the source has
    /// them, and no other order is taken.
    ///
    /// Their phenotype, their rows of the design and, when the pass over
    /// the variants reads them, their dosages are three lists read row by
    /// row: an order that is not the source's measures one individual's
    /// trait against another individual's genotypes, and every number of
    /// the study is then of nobody. The panel of
    /// `tests/reference/gwas/phenotypes.csv` has its phenotype in the
    /// source's order, so a study that sorted the individuals by any other
    /// key, or that kept them in a set and lost their order, gives the same
    /// answer on it as a study that is right.
    ///
    /// The fixture is four of the six individuals of a source, at the
    /// positions 0, 2, 3 and 5, and the two orders that are refused are the
    /// same four with 2 and 3 exchanged and the four reversed. Sorting
    /// either of the two, or reading it as a set, gives the order that is
    /// kept, so an implementation that did that fails here. The phenotype
    /// and the covariate of each individual differ from the others, and the
    /// design that comes back is asserted to be the buffer that was given,
    /// so no row was moved or dropped.
    #[test]
    fn the_tested_individuals_are_in_the_order_the_source_has_them() {
        let in_the_source_order = [0, 2, 3, 5];
        let study = a_study(&PHENOTYPE_OF_FOUR, &DESIGN_OF_FOUR, 2, &in_the_source_order);
        let design = Design::of_the_study(&study, 6)
            .expect("four of the six individuals of the source, in its order");
        assert_eq!(design.num_individuals(), 4);
        assert_eq!(design.num_coefs(), 2);
        assert!(
            ptr::eq(design.values(), DESIGN_OF_FOUR.as_slice()),
            "the design is the rows that were given, in the order they were given"
        );
        for (out_of_order, expected) in [([0, 3, 2, 5], (2, 3)), ([5, 3, 2, 0], (3, 5))] {
            let study = a_study(&PHENOTYPE_OF_FOUR, &DESIGN_OF_FOUR, 2, &out_of_order);
            match Design::of_the_study(&study, 6) {
                Err(Error::GwasIndividualsOutOfOrder { individual, after }) => {
                    assert_eq!(
                        (individual, after),
                        expected,
                        "the individuals {out_of_order:?} are not in the source's order"
                    );
                }
                other => panic!(
                    "the individuals {out_of_order:?} are the four that are tested in \
                     another order than the source has them, and that gave {other:?}"
                ),
            }
        }
    }

    /// An individual that the source does not have is refused, naming the
    /// position that was asked for and how many individuals there are. The
    /// fixture asks for the position 6 of a source of 6, whose positions
    /// are 0 to 5.
    #[test]
    fn a_tested_individual_the_source_does_not_have_is_refused() {
        let past_the_source = [0, 2, 3, 6];
        let study = a_study(&PHENOTYPE_OF_FOUR, &DESIGN_OF_FOUR, 2, &past_the_source);
        match Design::of_the_study(&study, 6) {
            Err(Error::GwasIndividualNotInTheDataset {
                individual,
                num_individuals,
            }) => {
                assert_eq!((individual, num_individuals), (6, 6));
            }
            other => panic!(
                "the source has the individuals 0 to 5 and the individual 6 was asked \
                 to be tested, and that gave {other:?}"
            ),
        }
    }

    /// An individual that is twice among the ones to test is refused,
    /// naming it. It would weigh twice in the null model and its phenotype
    /// would be read at two rows.
    #[test]
    fn a_tested_individual_given_twice_is_refused() {
        let twice = [0, 2, 2, 5];
        let study = a_study(&PHENOTYPE_OF_FOUR, &DESIGN_OF_FOUR, 2, &twice);
        match Design::of_the_study(&study, 6) {
            Err(Error::GwasIndividualTestedTwice { individual }) => {
                assert_eq!(individual, 2);
            }
            other => {
                panic!("the individual 2 is twice among the ones to test, and that gave {other:?}")
            }
        }
    }

    /// A study of no more individuals than the columns of its design plus
    /// one is refused: the design fits one number per column, the variant
    /// one more, and what is left over is what the uncertainty of the
    /// variant's effect is measured from.
    ///
    /// The fixture is a design of two columns, the intercept and one
    /// covariate, over three individuals, which is the columns plus one,
    /// and the same design over four, which is the fewest that is kept. The
    /// two are either side of the refusal, so a study that counted one
    /// individual more or less fails here.
    #[test]
    fn a_study_with_no_more_individuals_than_the_columns_plus_one_is_refused() {
        let three = [0, 1, 2];
        let phenotype_of_three = [1.5, -0.5, 2.0];
        let design_of_three = [1.0, 0.5, 1.0, -1.5, 1.0, 2.5];
        let study = a_study(&phenotype_of_three, &design_of_three, 2, &three);
        match Design::of_the_study(&study, 6) {
            Err(Error::GwasTooFewIndividuals {
                num_individuals,
                num_coefs,
            }) => {
                assert_eq!((num_individuals, num_coefs), (3, 2));
            }
            other => panic!(
                "three individuals and a design of two columns leave nothing to measure \
                 a variant's uncertainty from, and that gave {other:?}"
            ),
        }
        let four = [0, 1, 2, 3];
        let study = a_study(&PHENOTYPE_OF_FOUR, &DESIGN_OF_FOUR, 2, &four);
        assert!(
            Design::of_the_study(&study, 6).is_ok(),
            "four individuals and a design of two columns leave one individual over, \
             which is the fewest a study is made of"
        );
    }

    /// A phenotype of a binomial trait that is not 0 or 1 is refused,
    /// naming where it is among the tested individuals and what it is. The
    /// same values as a continuous trait are a measurement and are kept, so
    /// the refusal reads the trait it was given.
    #[test]
    fn a_binomial_phenotype_that_is_not_0_or_1_is_refused() {
        let four = [0, 1, 2, 3];
        let phenotype = [0.0, 1.0, 2.0, 1.0];
        let mut study = a_study(&phenotype, &DESIGN_OF_FOUR, 2, &four);
        study.trait_type = TraitType::Binomial;
        match Design::of_the_study(&study, 4) {
            Err(Error::GwasPhenotypeNotBinomial { position, value }) => {
                assert_eq!(position, 2);
                assert_eq!(value.to_bits(), 2.0_f64.to_bits());
            }
            other => panic!(
                "the third of the four individuals has the phenotype 2 of a binomial \
                 trait, which is 0 or 1, and that gave {other:?}"
            ),
        }
        study.trait_type = TraitType::Continuous;
        assert!(
            Design::of_the_study(&study, 4).is_ok(),
            "the same four values of a continuous trait are four measurements"
        );
    }

    /// A binomial trait where every tested individual has the same value is
    /// refused, naming the value: one of the two groups it compares is
    /// empty. Both 0 and 1 are the same refusal, and the same phenotype of
    /// a continuous trait is kept, which is what the spec says for that
    /// trait and what a variant with no variance is tested against.
    #[test]
    fn a_binomial_phenotype_of_one_value_is_refused() {
        let four = [0, 1, 2, 3];
        for same in [0.0, 1.0] {
            let phenotype = [same; 4];
            let mut study = a_study(&phenotype, &DESIGN_OF_FOUR, 2, &four);
            study.trait_type = TraitType::Binomial;
            match Design::of_the_study(&study, 4) {
                Err(Error::GwasPhenotypeOfOneValue { value }) => {
                    assert_eq!(value.to_bits(), same.to_bits());
                }
                other => panic!(
                    "every tested individual has the phenotype {same} of a binomial \
                     trait, and that gave {other:?}"
                ),
            }
            study.trait_type = TraitType::Continuous;
            assert!(
                Design::of_the_study(&study, 4).is_ok(),
                "a continuous trait of one value is refused nowhere in the spec"
            );
        }
        let phenotype = [0.0, 1.0, 1.0, 0.0];
        let mut study = a_study(&phenotype, &DESIGN_OF_FOUR, 2, &four);
        study.trait_type = TraitType::Binomial;
        assert!(
            Design::of_the_study(&study, 4).is_ok(),
            "two individuals have the condition and two have not"
        );
    }

    /// A phenotype that is not a finite number is refused, naming where it
    /// is and what it is. The individuals that are tested are those that
    /// have a phenotype, so a NaN is an individual that should not have
    /// been tested at all, and an infinity would carry through the null
    /// model into the effect of every variant.
    #[test]
    fn a_phenotype_that_is_not_finite_is_refused() {
        let four = [0, 1, 2, 3];
        for (position, not_finite) in [(1, f64::NAN), (3, f64::INFINITY)] {
            let mut phenotype = PHENOTYPE_OF_FOUR;
            phenotype[position] = not_finite;
            let study = a_study(&phenotype, &DESIGN_OF_FOUR, 2, &four);
            match Design::of_the_study(&study, 4) {
                Err(Error::GwasPhenotypeNotFinite {
                    position: where_it_is,
                    value,
                }) => {
                    assert_eq!(where_it_is, position);
                    assert_eq!(value.to_bits(), not_finite.to_bits());
                }
                other => panic!(
                    "the phenotype of the individual {position} is {not_finite}, and \
                     that gave {other:?}"
                ),
            }
        }
    }

    /// A value of the design that is not a finite number is refused,
    /// naming the individual whose row it is in, the column it is in and
    /// the value.
    ///
    /// The Python and the TypeScript layers refuse a covariate that is
    /// missing or is not a number, so what reaches this is a covariate
    /// that came out of a user's own arithmetic as an infinity, and a
    /// caller of the core crate. What the refusal is for is where the
    /// value would go otherwise: the rank refuses what it is given, and
    /// the user would be told that an operation of the linear algebra
    /// could not be done, which is how popnei says it has a defect. That
    /// is asserted here too, by giving the design of the four individuals
    /// a NaN and an infinity in turn and seeing this refusal and not
    /// [`Error::GwasLinalg`].
    ///
    /// The fixture puts the value in the column 1, the covariate, of the
    /// individual 2, and then in the column 0, the intercept, of the
    /// individual 0, so that the two numbers the message carries are told
    /// apart from each other.
    #[test]
    fn a_design_value_that_is_not_a_number_is_refused() {
        let four = [0, 1, 2, 3];
        for (individual, coef, not_finite) in [(2, 1, f64::NAN), (0, 0, f64::INFINITY)] {
            let mut design = DESIGN_OF_FOUR;
            let value = design
                .get_mut(individual * 2 + coef)
                .expect("the value of the design");
            *value = not_finite;
            let study = a_study(&PHENOTYPE_OF_FOUR, &design, 2, &four);
            match Design::of_the_study(&study, 4) {
                Err(Error::GwasDesignValueNotFinite {
                    individual: whose,
                    coef: which,
                    value,
                }) => {
                    assert_eq!((whose, which), (individual, coef));
                    assert_eq!(value.to_bits(), not_finite.to_bits());
                }
                other => panic!(
                    "the column {coef} of the individual {individual} of the design is \
                     {not_finite}, and that gave {other:?}"
                ),
            }
        }
    }

    /// A design with a covariate that is twice another is refused: its
    /// columns are not independent, so its effects are not one set of
    /// numbers but many.
    ///
    /// The fixture is six individuals and three columns, the intercept, a
    /// covariate and that covariate doubled, whose rank is 2 of the 3
    /// columns. The covariate is not constant and has both signs, so the
    /// refusal is the doubling and not the covariate itself.
    #[test]
    fn a_design_whose_covariate_is_twice_another_is_refused() {
        let six = [0, 1, 2, 3, 4, 5];
        let phenotype = [1.5, -0.5, 2.0, 0.25, 3.5, -1.25];
        let design: Vec<f64> = [0.5, 1.5, -2.0, 3.0, 0.25, -1.0]
            .into_iter()
            .flat_map(|covariate| [1.0, covariate, covariate * 2.0])
            .collect();
        let study = a_study(&phenotype, &design, 3, &six);
        match Design::of_the_study(&study, 6) {
            Err(Error::GwasCovariatesCollinear { num_coefs, rank }) => {
                assert_eq!((num_coefs, rank), (3, 2));
            }
            other => panic!(
                "the third column of the design is the second doubled, so two of its \
                 three columns are independent, and that gave {other:?}"
            ),
        }
    }

    /// A design whose smallest singular value is 1e-11 of its largest is
    /// kept, and one whose smallest is 1e-16 of its largest is refused.
    ///
    /// The rank counts the singular values strictly above the largest of
    /// them times the larger dimension times the distance from 1 to the
    /// next `f64`, which `docs/specs/linalg.md` takes from numpy's
    /// `matrix_rank` so that a design popnei refuses is a design pyNei
    /// refuses. For the eight rows and two columns here that tolerance is
    /// 8 * 2.220446049250313e-16 = 1.7763568394002505e-15 of the largest
    /// singular value, so the design that is kept has its smallest 5629
    /// times above it and the one that is refused has its smallest 17.8
    /// times below. The two are 1e5 apart and the tolerance is between
    /// them, which is what says that the rank was taken at numpy's
    /// tolerance and not at a rounder one. Searched on this Mac on 23
    /// September 2026, the fraction popnei turns at is between
    /// 1.7763567728318708e-15 and 1.7763568947297103e-15, which holds that
    /// tolerance to seven digits; the boundary itself is not asserted,
    /// since a singular value the decomposition computes within a bit of
    /// the tolerance falls either side of it by rounding.
    ///
    /// numpy 2.5.3 was run on the same two designs on 23 September 2026 and
    /// gave the same two ranks, 2 and 1, and singular values of
    /// 2.8284271247461903 and `size` times that, which is what two columns
    /// at right angles have.
    ///
    /// 1e-11 is the fraction `docs/specs/linalg.md` measured `matrix_rank`
    /// to give full rank at on a design of 10000 rows and 4 columns, where
    /// the tolerance is 10000 times that distance; on the eight rows here
    /// the same fraction is further from the tolerance, since the tolerance
    /// follows the larger dimension.
    #[test]
    fn a_design_whose_smallest_singular_value_is_1e_11_of_its_largest_is_kept() {
        let eight = [0, 1, 2, 3, 4, 5, 6, 7];
        let phenotype = [1.5, -0.5, 2.0, 0.25, 3.5, -1.25, 0.75, 2.25];
        let kept = a_design_whose_smallest_singular_value_is(1e-11);
        let study = a_study(&phenotype, &kept, 2, &eight);
        assert!(
            Design::of_the_study(&study, 8).is_ok(),
            "the smallest singular value of this design is 1e-11 of its largest, 5629 \
             times the tolerance of the rank"
        );
        let refused = a_design_whose_smallest_singular_value_is(1e-16);
        let study = a_study(&phenotype, &refused, 2, &eight);
        match Design::of_the_study(&study, 8) {
            Err(Error::GwasCovariatesCollinear { num_coefs, rank }) => {
                assert_eq!((num_coefs, rank), (2, 1));
            }
            other => panic!(
                "the smallest singular value of this design is 1e-16 of its largest, \
                 17.8 times below the tolerance of the rank, and that gave {other:?}"
            ),
        }
    }

    /// The phenotype holds one value for each tested individual, the design
    /// one row of its columns for each, and the design has the column of
    /// the intercept at least. None of the three can be reached from Python
    /// or from TypeScript, which build the three from the same individuals,
    /// and a caller of the core crate builds them itself.
    #[test]
    fn a_study_whose_buffers_are_not_of_its_individuals_is_refused() {
        let four = [0, 1, 2, 3];
        let phenotype_of_three = [1.5, -0.5, 2.0];
        let study = a_study(&phenotype_of_three, &DESIGN_OF_FOUR, 2, &four);
        assert_eq!(
            the_shape_refused(Design::of_the_study(&study, 4)),
            GwasInputShape::Phenotype {
                num_values: 3,
                num_individuals: 4,
            }
        );
        let design_of_three = [1.0, 0.5, 1.0, -1.5, 1.0, 2.5];
        let study = a_study(&PHENOTYPE_OF_FOUR, &design_of_three, 2, &four);
        assert_eq!(
            the_shape_refused(Design::of_the_study(&study, 4)),
            GwasInputShape::Design {
                num_values: 6,
                num_individuals: 4,
                num_coefs: 2,
            }
        );
        let study = a_study(&PHENOTYPE_OF_FOUR, &[], 0, &four);
        assert_eq!(
            the_shape_refused(Design::of_the_study(&study, 4)),
            GwasInputShape::NoCoef
        );
    }

    /// Which of the three buffers of a study was of another size, from what
    /// building its design gave.
    fn the_shape_refused(built: crate::error::Result<Design<'_>>) -> GwasInputShape {
        match built {
            Err(Error::GwasInputOfAnotherSize { problem }) => problem,
            other => panic!("the buffers of the study do not hold it, and that gave {other:?}"),
        }
    }
}

/// The panels the tests of the dosages and of the result are read over,
/// and what every one of them is asserted against.
///
/// Each one has individuals that are not tested, which is what a phenotype
/// that leaves somebody out gives and the case that tells a study over the
/// tested individuals apart from one over the whole panel.
#[cfg(test)]
mod fixtures {
    use super::{BlockOfThePass, Design, GwasInput, TraitType};
    use crate::block::Block;
    use crate::variant::MISSING_ALLELE;

    /// An allele that was not called, which is what a `.` of a VCF is.
    pub(super) const MISSING: i8 = MISSING_ALLELE;

    /// How far a value of these tests may be from the literal beside it:
    /// 1e-15 absolute. The two kinds of value it is used on are a
    /// frequency, which lies between 0 and 1, and a dosage, which lies
    /// between 0 and the ploidy, so the bound is 1e-15 of the largest
    /// either can hold and not a share of a value that may be near 0,
    /// which is what "How it is verified" of `docs/specs/gwas.md` asks a
    /// tolerance to be measured against.
    ///
    /// Every one of these numbers is a sum of whole dosages divided by a
    /// count and then by the ploidy, two divisions of numbers below 2^53,
    /// which IEEE 754 rounds the same way everywhere, so the distance is
    /// expected to be 0 and was measured at 0 over every value of every
    /// test here on 23 September 2026, with the bound itself set to 0, on
    /// Accelerate and on faer. What the bound is for is a literal that a
    /// person reads and writes as sixteen digits.
    pub(super) const OF_A_FREQUENCY: f64 = 1e-15;

    /// A study of a continuous trait over the individuals at `individuals`
    /// of the source, with no kinship, the default test and the design of
    /// an intercept and one covariate.
    pub(super) fn a_study<'a>(
        phenotype: &'a [f64],
        design: &'a [f64],
        individuals: &'a [usize],
    ) -> GwasInput<'a> {
        GwasInput {
            phenotype,
            trait_type: TraitType::Continuous,
            design,
            num_coefs: 2,
            kinship: None,
            test: None,
            use_grammar_gamma_approx: false,
            individuals,
            transform_to_biallelic: false,
        }
    }

    /// The same study of the individuals at `individuals`, with a
    /// phenotype and a design of one value and one row for each of them.
    /// The values are not read by anything these tests call, which reads
    /// the genotypes and the individuals; they are here because a study
    /// carries them and because [`super::Design::of_the_study`] is asked
    /// for the design of the first fixture below.
    pub(super) fn the_phenotype_and_the_design_of(individuals: &[usize]) -> (Vec<f64>, Vec<f64>) {
        let phenotype: Vec<f64> = individuals
            .iter()
            .map(|individual| 0.25 + *individual as f64)
            .collect();
        let design: Vec<f64> = individuals
            .iter()
            .flat_map(|individual| [1.0, 0.5 - *individual as f64])
            .collect();
        (phenotype, design)
    }

    /// The design of a study whose source has `num_individuals_of_the_source`
    /// individuals, checked, which is what the dosages of a block are read
    /// over: it carries the tested individuals and the choice about a
    /// variant of more than two alleles.
    pub(super) fn the_design_of<'a>(
        study: &GwasInput<'a>,
        num_individuals_of_the_source: usize,
    ) -> Design<'a> {
        Design::of_the_study(study, num_individuals_of_the_source)
            .expect("the design of the tested individuals")
    }

    /// Where the one block of a fixture sits in the pass that gives it:
    /// the ploidy the reader says its source has, and the first variant.
    pub(super) fn the_first_block_of(ploidy: usize) -> BlockOfThePass {
        BlockOfThePass {
            ploidy,
            first_var: 0,
        }
    }

    /// A block of `rows.len()` variants over `num_individuals`
    /// individuals of `ploidy`, with the genotypes of each variant given
    /// as one row.
    pub(super) fn a_block(num_individuals: usize, ploidy: usize, rows: &[&[i8]]) -> Block {
        Block {
            num_vars: rows.len(),
            num_individuals,
            ploidy,
            gts: rows.concat(),
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        }
    }

    /// Every value of `found` within [`OF_A_FREQUENCY`] of the literal
    /// beside it, and as many values as there are literals.
    pub(super) fn assert_the_values(found: &[f64], expected: &[f64], what: &str) {
        assert_eq!(found.len(), expected.len(), "how many values {what} holds");
        for (position, (found, expected)) in found.iter().zip(expected).enumerate() {
            assert!(
                (found - expected).abs() <= OF_A_FREQUENCY,
                "{what}: the value {position} is {found} and the literal is {expected}"
            );
        }
    }

    /// The four variants of a panel of eight diploid individuals, of which
    /// the four at 0, 2, 4 and 6 have a phenotype and are tested.
    ///
    /// Each variant is one of the four things the individuals that are
    /// left out can do to a study that read them:
    ///
    /// - `v0`, where every tested individual is heterozygous: it has no
    ///   variance among them and cannot be tested, and the four that are
    ///   left out hold `0/0` and `1/1`, so over the panel it has variance
    ///   and would be tested.
    /// - `v1`: its frequency is 0.5 over the tested individuals and 0.25
    ///   over the panel, because the four that are left out are all `0/0`.
    /// - `v2`, where a tested individual has no genotype: it takes the
    ///   mean dosage of the other three, and the four left out are `1/1`,
    ///   which makes 1 the major allele of the panel and not of the tested
    ///   individuals, so every dosage of the panel is counted from the
    ///   other allele.
    /// - `v3`, where no tested individual has a genotype: its dosages are
    ///   the mean of nothing, which is 0, and the four left out do have
    ///   genotypes and would give it a frequency of 0.5.
    pub(super) const OF_EIGHT: [&[i8]; 4] = [
        // v0: 0/1 0/0 0/1 0/0 0/1 1/1 0/1 1/1
        &[0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1],
        // v1: 0/0 0/0 0/1 0/0 1/1 0/0 0/1 0/0
        &[0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0],
        // v2: 0/0 1/1 ./. 1/1 1/1 1/1 0/0 1/1
        &[0, 0, 1, 1, MISSING, MISSING, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1],
        // v3: ./. 0/0 ./. 0/1 ./. 1/1 ./. 0/1
        &[
            MISSING, MISSING, 0, 0, MISSING, MISSING, 0, 1, MISSING, MISSING, 1, 1, MISSING,
            MISSING, 0, 1,
        ],
    ];

    /// Which of the eight individuals of [`OF_EIGHT`] have a phenotype and
    /// are tested: half of them.
    pub(super) const TESTED_OF_EIGHT: [usize; 4] = [0, 2, 4, 6];

    /// All eight, which is the study the same panel gives when nobody is
    /// left out, and what the tested four are asserted to differ from.
    pub(super) const THE_PANEL_OF_EIGHT: [usize; 8] = [0, 1, 2, 3, 4, 5, 6, 7];

    /// The frequency of each variant of [`OF_EIGHT`] over the four tested
    /// individuals, worked out by hand from the fixture.
    ///
    /// The dosage of a genotype is how many of its alleles are not the
    /// major one of its variant among these four, and the frequency is the
    /// mean dosage of the called genotypes over the ploidy, 2.
    ///
    /// - `v0`: their genotypes are `0/1 0/1 0/1 0/1`, so the alleles 0 and
    ///   1 are called four times each and the major one is 0, the lower
    ///   numbered of two that are equally frequent; the dosages are
    ///   1 1 1 1, their mean is 1 and the frequency is 0.5.
    /// - `v1`: `0/0 0/1 1/1 0/1`, four 0s and four 1s, the major allele 0,
    ///   the dosages 0 1 2 1, their mean 4 / 4 = 1 and the frequency 0.5.
    /// - `v2`: `0/0 ./. 1/1 0/0`, four 0s and two 1s, the major allele 0,
    ///   the dosages 0, none, 2 and 0, their mean 2 / 3 and the frequency
    ///   1 / 3 = 0.3333333333333333.
    /// - `v3`: `./. ./. ./. ./.`, no allele called, so the mean of nothing
    ///   is 0 and the frequency is 0.
    pub(super) const FREQUENCIES_OF_THE_TESTED: [f64; 4] = [0.5, 0.5, 0.3333333333333333, 0.0];

    /// The frequency of each variant of [`OF_EIGHT`] over all eight
    /// individuals, worked out by hand the same way, which is what a study
    /// that read the four individuals without a phenotype would give.
    ///
    /// - `v0`: eight 0s and eight 1s, the major allele 0, the dosages
    ///   1 0 1 0 1 2 1 2, their mean 8 / 8 = 1 and the frequency 0.5, the
    ///   one of the four that does not move; what moves is that the
    ///   variant now has variance.
    /// - `v1`: twelve 0s and four 1s, the major allele 0, the dosages
    ///   0 0 1 0 2 0 1 0, their mean 4 / 8 and the frequency 0.25.
    /// - `v2`: four 0s and ten 1s, so the major allele is 1 here and 0
    ///   among the tested four; the dosages are 2 0, none, 0 0 0 2 0,
    ///   their mean 4 / 7 and the frequency 2 / 7 = 0.2857142857142857.
    /// - `v3`: four 0s and four 1s, the major allele 0, the dosages none,
    ///   0, none, 1, none, 2, none, 1, their mean 4 / 4 = 1 and the
    ///   frequency 0.5.
    pub(super) const FREQUENCIES_OF_THE_PANEL: [f64; 4] = [0.5, 0.25, 0.2857142857142857, 0.5];
}

/// The dosages of a block over the individuals a study tests, which are
/// theirs and not the whole panel's, as "What it gives" of
/// `docs/specs/gwas.md` asks.
#[cfg(test)]
mod dosages {
    use super::fixtures::{
        FREQUENCIES_OF_THE_PANEL, FREQUENCIES_OF_THE_TESTED, MISSING, OF_EIGHT, TESTED_OF_EIGHT,
        THE_PANEL_OF_EIGHT, a_block, a_study, assert_the_values, the_design_of, the_first_block_of,
        the_phenotype_and_the_design_of,
    };
    use super::{BlockOfThePass, GwasDosages};
    use crate::error::Error;

    /// The dosages, the major allele, the mean a genotype with an allele
    /// missing takes and the frequency are all of the individuals that
    /// have a phenotype, and the four that have none change every one of
    /// them.
    ///
    /// The fixture is the panel of eight diploid individuals of
    /// `fixtures::OF_EIGHT`, of which the four at 0, 2, 4 and 6 are
    /// tested: the doc comments of `FREQUENCIES_OF_THE_TESTED` and of
    /// `FREQUENCIES_OF_THE_PANEL` work every frequency here out by hand
    /// from the genotypes, and this test asserts both, the study of the
    /// four and the study of all eight, so that neither is read off the
    /// other. Three of the four frequencies move when the individuals
    /// without a phenotype are read, by 0.25, by 0.048 and by 0.5, which
    /// are 2.5e14, 4.8e13 and 5e14 times the 1e-15 the values are held to.
    ///
    /// The fourth, `v0`, keeps its frequency and loses its answer instead:
    /// every tested individual is heterozygous, so it has no variance
    /// among them and cannot be tested, while over the panel it has
    /// variance and would be. A study that read the panel would give it
    /// three numbers where it has none.
    #[test]
    fn the_dosages_and_the_frequencies_are_of_the_tested_individuals_and_not_of_the_panel() {
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        // The four individuals that have a phenotype are four of the eight
        // the source has, in its order, so this is a study that runs, and
        // the design that comes out is what carries them to the block.
        let design = the_design_of(&study, 8);
        let mut block = a_block(8, 2, &OF_EIGHT);
        let mut dosages = GwasDosages::of_a_study();

        dosages
            .read_the_block(&mut block, &design, the_first_block_of(2))
            .expect("the dosages of the block over the four tested individuals");

        assert_eq!(dosages.num_vars(), 4);
        assert_eq!(dosages.num_individuals(), 4, "the tested individuals");
        assert_the_values(
            dosages.allele_freq(),
            &FREQUENCIES_OF_THE_TESTED,
            "the frequencies over the tested individuals",
        );
        assert_eq!(
            dosages.has_variance(),
            [false, true, true, false],
            "v0 is heterozygous in all four and v3 has no genotype among them"
        );
        assert_eq!(dosages.num_with_variance(), 2);
        // The rows of the two variants that have variance, in the order of
        // the block and at the start of the buffer: `v1` is 0/0 0/1 1/1
        // 0/1, whose dosages are 0 1 2 1, and `v2` is 0/0 ./. 1/1 0/0,
        // whose dosages are 0, the mean 2 / 3 of the other three, 2 and 0.
        assert_the_values(
            dosages.dosages(),
            &[0.0, 1.0, 2.0, 1.0, 0.0, 0.6666666666666666, 2.0, 0.0],
            "the dosages of the variants that have variance",
        );

        let (phenotype, design) = the_phenotype_and_the_design_of(&THE_PANEL_OF_EIGHT);
        let study = a_study(&phenotype, &design, &THE_PANEL_OF_EIGHT);
        let design = the_design_of(&study, 8);
        let mut block = a_block(8, 2, &OF_EIGHT);

        dosages
            .read_the_block(&mut block, &design, the_first_block_of(2))
            .expect("the dosages of the block over all eight individuals");

        assert_eq!(dosages.num_individuals(), 8, "the whole panel");
        assert_the_values(
            dosages.allele_freq(),
            &FREQUENCIES_OF_THE_PANEL,
            "the frequencies over the whole panel",
        );
        assert_eq!(
            dosages.has_variance(),
            [true, true, true, true],
            "every variant of the panel has variance, `v0` and `v3` among them"
        );
        // `v2` is where the two studies count the dosages from different
        // alleles: over the panel the major allele is 1, so the first
        // individual, 0/0, has the dosage 2, where over the tested four it
        // has 0.
        assert_the_values(
            dosages.dosages().get(16..24).unwrap_or_default(),
            &[2.0, 0.0, 0.5714285714285714, 0.0, 0.0, 0.0, 2.0, 0.0],
            "the dosages of `v2` over the whole panel",
        );
    }

    /// The three variants of "The worked example" of
    /// `docs/specs/gwas.md`, whose dosages and frequencies the spec gives,
    /// read over its six individuals while two more are in the block and
    /// have no phenotype.
    ///
    /// The spec's dosages are 0 1 2 0 1 2 for `v0`, 0 1 2 0.8 1 0 for
    /// `v1`, where the 0.8 is the genotype of `i3` that was not called
    /// taking the mean `(0 + 1 + 2 + 1 + 0) / 5` of its variant, and
    /// 1 1 1 1 1 1 for `v2`, which has no variance and is not kept; the
    /// frequencies are 0.5, 0.4 and 0.5.
    ///
    /// The two individuals that are left out are `1/1` at all three
    /// variants, which is not what the six hold: it makes 1 the major
    /// allele of every one of the three, so every dosage is counted from
    /// the other allele, and it gives `v2` two dosages and an answer it
    /// does not have among the six. Over all eight the frequencies are
    /// 0.375, 3 / 7 and 0.375, so a study that read them misses all three
    /// of the spec's numbers by 0.125, 0.029 and 0.125.
    #[test]
    fn the_dosages_of_the_worked_example_are_of_its_six_individuals_and_not_of_the_two_beside_them()
    {
        let of_the_worked_example: [&[i8]; 3] = [
            // v0: 0/0 0/1 1/1 0/0 0/1 1/1 | 1/1 1/1
            &[0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1],
            // v1: 0/0 0/1 1/1 ./. 0/1 0/0 | 1/1 1/1
            &[0, 0, 0, 1, 1, 1, MISSING, MISSING, 0, 1, 0, 0, 1, 1, 1, 1],
            // v2: 0/1 0/1 0/1 0/1 0/1 0/1 | 1/1 1/1
            &[0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1],
        ];
        let six = [0, 1, 2, 3, 4, 5];
        let (phenotype, design) = the_phenotype_and_the_design_of(&six);
        let study = a_study(&phenotype, &design, &six);
        let design = the_design_of(&study, 8);
        let mut block = a_block(8, 2, &of_the_worked_example);
        let mut dosages = GwasDosages::of_a_study();

        dosages
            .read_the_block(&mut block, &design, the_first_block_of(2))
            .expect("the dosages of the worked example");

        assert_the_values(
            dosages.allele_freq(),
            &[0.5, 0.4, 0.5],
            "the frequencies of the worked example",
        );
        assert_eq!(dosages.has_variance(), [true, true, false]);
        assert_the_values(
            dosages.dosages(),
            &[0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.8, 1.0, 0.0],
            "the dosages of the worked example",
        );

        let eight = [0, 1, 2, 3, 4, 5, 6, 7];
        let (phenotype, design) = the_phenotype_and_the_design_of(&eight);
        let study = a_study(&phenotype, &design, &eight);
        let design = the_design_of(&study, 8);
        let mut block = a_block(8, 2, &of_the_worked_example);

        dosages
            .read_the_block(&mut block, &design, the_first_block_of(2))
            .expect("the dosages of the worked example and the two beside it");

        assert_the_values(
            dosages.allele_freq(),
            &[0.375, 0.42857142857142855, 0.375],
            "the frequencies of the eight individuals of the block",
        );
        assert_eq!(
            dosages.has_variance(),
            [true, true, true],
            "`v2` has variance once the two individuals that are 1/1 are read"
        );
    }

    /// How many alleles a variant has is counted among the tested
    /// individuals too, so a third allele that only an individual without
    /// a phenotype carries is not a variant of more than two alleles for
    /// the study, and one a tested individual carries is refused with its
    /// position among the variants the reader has given.
    ///
    /// The fixture is six triploid individuals, of which the four at 1, 2,
    /// 3 and 5 are tested, and two variants, the second of which has the
    /// allele 2 in the individual 0, who has no phenotype. The position
    /// asserted is 101 and not 1, because the block is given as the
    /// hundred and first variant of the reader.
    ///
    /// With `transform_to_biallelic` the whole panel is read and every
    /// allele that is not the major one counts the same: the alleles of
    /// the second variant over the six are fourteen 0s, three 1s and one
    /// 2, so the major one is 0 and the dosages are 1 0 0 1 0 2, whose
    /// mean is 4 / 6 and whose frequency over the ploidy of 3 is
    /// 2 / 9 = 0.2222222222222222.
    #[test]
    fn a_third_allele_is_counted_among_the_tested_individuals_and_refused_there() {
        let of_six_triploids: [&[i8]; 2] = [
            // w0: 0/0/0 0/0/1 1/1/1 0/1/1 0/0/0 1/1/1
            &[0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1],
            // w1: 0/0/2 0/0/0 0/0/0 0/0/1 0/0/0 0/1/1
            &[0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1],
        ];
        let four = [1, 2, 3, 5];
        let (phenotype, design) = the_phenotype_and_the_design_of(&four);
        let study = a_study(&phenotype, &design, &four);
        let of_four = the_design_of(&study, 6);
        let mut block = a_block(6, 3, &of_six_triploids);
        let mut dosages = GwasDosages::of_a_study();
        let of_the_pass = BlockOfThePass {
            ploidy: 3,
            first_var: 100,
        };

        dosages
            .read_the_block(&mut block, &of_four, of_the_pass)
            .expect("the allele 2 is the individual 0's, who has no phenotype");

        // `w0` over the four tested individuals is 0/0/1 1/1/1 0/1/1
        // 1/1/1, three 0s and nine 1s, so the major allele is 1 and the
        // dosages are 2 0 1 0, whose mean is 3 / 4 and whose frequency
        // over the ploidy of 3 is 0.25. `w1` is 0/0/0 0/0/0 0/0/1 0/1/1,
        // nine 0s and three 1s, the major allele 0 and the dosages
        // 0 0 1 2, the same mean and the same frequency.
        assert_the_values(
            dosages.allele_freq(),
            &[0.25, 0.25],
            "the frequencies over the four tested individuals",
        );
        assert_the_values(
            dosages.dosages(),
            &[2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0],
            "the dosages over the four tested individuals",
        );

        let six = [0, 1, 2, 3, 4, 5];
        let (phenotype, values) = the_phenotype_and_the_design_of(&six);
        let study = a_study(&phenotype, &values, &six);
        let of_six = the_design_of(&study, 6);
        let mut block = a_block(6, 3, &of_six_triploids);

        match dosages.read_the_block(&mut block, &of_six, of_the_pass) {
            Err(Error::VariantWithMoreThanTwoAlleles {
                position,
                num_alleles,
            }) => {
                assert_eq!((position, num_alleles), (101, 3));
            }
            other => panic!(
                "the individual 0 carries the allele 2 of the second variant of the \
                 block, which is the variant 101 of the reader, and that gave {other:?}"
            ),
        }

        let mut collapsed = a_study(&phenotype, &values, &six);
        collapsed.transform_to_biallelic = true;
        let collapsed = the_design_of(&collapsed, 6);
        let mut block = a_block(6, 3, &of_six_triploids);

        dosages
            .read_the_block(&mut block, &collapsed, of_the_pass)
            .expect("every allele that is not the major one counts the same");

        assert_the_values(
            dosages.allele_freq(),
            &[0.5, 0.2222222222222222],
            "the frequencies over the six individuals with the third allele collapsed",
        );
    }

    /// The rows of a block are read on the threads of rayon, and neither
    /// the dosages nor the frequencies depend on how many there are, nor
    /// on whether the rows were read one after another, which is what
    /// WebAssembly does.
    ///
    /// The block is 300 variants of 40 diploid individuals, of which the
    /// 14 whose position is a multiple of 3 are tested, so a pool of four
    /// shares the rows out in several chunks. The alleles of the genotype
    /// of the individual `k` at the variant `v` are read out of a pattern
    /// of 5 at `(7k + 3v) % 5` and one of 7 at `(11k + 5v) % 7`, both of
    /// them the missing allele, 0 and 1, and every fiftieth variant is
    /// `0/0` in everybody instead.
    ///
    /// What that holds, counted over the 14 tested individuals on 23
    /// September 2026: they take 6 to 8 different genotypes at each of the
    /// variants that are not the fiftieth, of the 9 the two patterns can
    /// give, and never fewer, since a tested individual is a multiple of 3
    /// and the two patterns then run over all of 5 and all of 7; 30.8 per
    /// cent of their genotypes have an allele missing and no variant is
    /// missing in all of them; no row is the row of its neighbour, and the
    /// 300 rows are 36 different ones, so a row read at the place of
    /// another is seen; and the 6 variants that are `0/0` in everybody
    /// have no variance among them, so the two runs have to agree about
    /// which rows are left out of the matrix and where the others moved
    /// to. The fixture before this one gave the tested individuals 3
    /// genotypes and 5 of them the same one at every variant.
    ///
    /// The pools are built here and are not rayon's global one, which has
    /// one thread per core of the machine. rayon is a dependency of the
    /// targets that are not wasm, so this test is compiled for those
    /// alone.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_rows_are_the_same_on_one_thread_on_several_and_read_one_after_another() {
        use std::num::NonZeroUsize;

        use super::{
            MultiallelicVariants, the_dosages_of_the_rows, the_dosages_of_the_rows_one_by_one,
        };

        /// Whether two lists of values hold the same numbers, to the bit,
        /// which is what two ways of reading the same genotypes give: each
        /// row is read on its own and nothing is summed across the rows,
        /// so neither the threads nor their number can move a digit.
        fn the_same_values(one: &[f64], other: &[f64]) -> bool {
            one.len() == other.len()
                && one
                    .iter()
                    .zip(other)
                    .all(|(one, other)| one.to_bits() == other.to_bits())
        }

        // The missing allele, 0 and 1 in a pattern of 5 and one of 7,
        // which no two of the tested individuals read at the same place.
        const OF_THE_FIRST: [i8; 5] = [0, 1, MISSING, 1, 0];
        const OF_THE_SECOND: [i8; 7] = [0, 1, 0, 1, MISSING, 1, 0];
        let genotype = |individual: usize, variant: usize| match variant % 50 {
            0 => [0, 0],
            _ => [
                OF_THE_FIRST[(individual * 7 + variant * 3) % 5],
                OF_THE_SECOND[(individual * 11 + variant * 5) % 7],
            ],
        };
        let of_forty = || {
            let rows: Vec<Vec<i8>> = (0..300)
                .map(|variant| {
                    (0..40)
                        .flat_map(|individual| genotype(individual, variant))
                        .collect()
                })
                .collect();
            let rows: Vec<&[i8]> = rows.iter().map(Vec::as_slice).collect();
            a_block(40, 2, &rows)
        };
        let tested: Vec<usize> = (0..40).filter(|individual| individual % 3 == 0).collect();
        let (phenotype, design) = the_phenotype_and_the_design_of(&tested);
        let study = a_study(&phenotype, &design, &tested);
        let design = the_design_of(&study, 40);
        let read_with = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            let mut block = of_forty();
            let mut dosages = GwasDosages::of_a_study();
            pool.install(|| dosages.read_the_block(&mut block, &design, the_first_block_of(2)))
                .expect("the dosages of the block");
            dosages
        };

        let on_one = read_with(1);
        let on_four = read_with(4);

        assert_eq!(on_one.num_vars(), 300);
        assert_eq!(on_one.num_individuals(), 14);
        assert_eq!(
            on_one.num_with_variance(),
            294,
            "the 6 variants that are 0/0 in everybody have no variance"
        );
        assert!(
            the_same_values(on_one.dosages(), on_four.dosages()),
            "the dosages on one thread and on four"
        );
        assert!(
            the_same_values(on_one.allele_freq(), on_four.allele_freq()),
            "the frequencies on one thread and on four"
        );
        assert_eq!(on_one.has_variance(), on_four.has_variance());

        // The rows read one after another, which is the pass WebAssembly
        // takes, against the rows read on the threads, over the genotypes
        // of the tested individuals alone: the block is compacted to them
        // here as the two passes above had it compacted for them.
        let mut block = of_forty();
        block
            .retain_individuals(&tested)
            .expect("the tested individuals");
        let alleles_per_var =
            NonZeroUsize::new(block.alleles_per_var().expect("the alleles of a variant"))
                .expect("the block holds genotypes");
        let num_individuals = NonZeroUsize::new(block.num_individuals).expect("the individuals");
        let mut on_threads = vec![0.0; block.gts.len() / 2];
        let rows = the_dosages_of_the_rows(
            &block.gts,
            alleles_per_var,
            num_individuals,
            2,
            MultiallelicVariants::Refused,
            0,
            &mut on_threads,
        )
        .expect("the rows on the threads");
        let mut one_by_one = vec![0.0; block.gts.len() / 2];
        let rows_one_by_one = the_dosages_of_the_rows_one_by_one(
            &block.gts,
            alleles_per_var,
            num_individuals,
            2,
            MultiallelicVariants::Refused,
            0,
            &mut one_by_one,
        )
        .expect("the rows one after another");

        assert!(
            the_same_values(&on_threads, &one_by_one),
            "the dosages of every row, on the threads and one after another"
        );
        assert_eq!(rows.len(), 300);
        for (var, (row, one_by_one)) in rows.iter().zip(&rows_one_by_one).enumerate() {
            assert_eq!(
                (row.allele_freq.to_bits(), row.has_variance),
                (one_by_one.allele_freq.to_bits(), one_by_one.has_variance),
                "the variant {var}"
            );
        }
    }

    /// The dosages of a study and the standardized row of
    /// [`crate::variant`] read a variant by the same three rules, and this
    /// is what fails the day one of them is changed and the other is not.
    ///
    /// The two are written separately because a study needs what the other
    /// one does not give: the mean of the variant, which it reports as
    /// `allele_freq`, and dosages that are not divided by the deviation of
    /// the variant, since `beta` is the effect of one more copy of an
    /// allele in the units of the trait. Everything before that division
    /// is the same work, and this asserts that it gives the same answers:
    /// which allele is the major one, and so what every dosage is counted
    /// from; the mean a genotype with an allele missing takes; whether the
    /// variant has variance at all; and the refusal of a variant with more
    /// than two alleles, at the same position and with the same count.
    ///
    /// The fixture is the panel of eight over its four tested individuals,
    /// which holds a variant with no variance among them, two that have
    /// it, one of which has a genotype that was not called, and one with
    /// nothing called at all. Each row of it is read both ways, and the
    /// standardized row is asserted to be this one centered at the mean
    /// this module gives and divided by the deviation worked out from it,
    /// which is `docs/specs/pca.md`'s divisor: the root of the mean square
    /// deviation over all the individuals, the ones whose genotype was not
    /// called counting as no deviation, since they hold the mean.
    ///
    /// The bound is 1e-12 of the largest value of the row, which is the
    /// scale of what is being compared, and not of each value, which is 0
    /// for an individual at the mean. The two paths add the same squares
    /// in a different order, this one over the individuals and the other
    /// over the dosages with a count on each, so the last bits may differ;
    /// measured over this fixture on 23 September 2026 the worst distance
    /// was 0 on Accelerate and on faer, with the bound set to 0.
    #[test]
    fn the_dosages_and_the_standardized_row_of_a_variant_agree() {
        use super::{DosageScratch, MultiallelicVariants, the_dosages_of_a_row};
        use crate::variant::{DosageOptions, DosageScale, RowScratch, the_standardized_row};

        let options = DosageOptions {
            transform_to_biallelic: false,
            scale: DosageScale::OfTheDosages,
        };
        let mut block = a_block(8, 2, &OF_EIGHT);
        block
            .retain_individuals(&TESTED_OF_EIGHT)
            .expect("the four tested individuals");
        let mut mine = DosageScratch::of(4);
        let mut theirs = RowScratch::of(4);
        let mut agreed: Vec<bool> = Vec::new();

        for (var, gts) in block.gts.as_chunks::<8>().0.iter().enumerate() {
            let mut row = [0.0; 4];
            let read = the_dosages_of_a_row(
                gts,
                2,
                var,
                MultiallelicVariants::Refused,
                &mut mine,
                &mut row,
            )
            .expect("the dosages of the row");
            let mut standardized = [0.0; 4];
            let used = the_standardized_row(gts, 2, var, &options, &mut theirs, &mut standardized)
                .expect("the standardized row");

            assert_eq!(
                read.has_variance, used,
                "the variant {var} has variance one way and not the other"
            );
            agreed.push(used);
            if !used {
                continue;
            }
            // The mean of the called dosages, which is the frequency times
            // the ploidy: the two multiply and divide by 2, which is exact.
            let mean = read.allele_freq * 2.0;
            let squares: f64 = row
                .iter()
                .map(|dosage| (dosage - mean) * (dosage - mean))
                .sum();
            let divisor = (squares / 4.0).sqrt();
            let largest = standardized
                .iter()
                .fold(0.0_f64, |largest, value| largest.max(value.abs()));
            for (position, (theirs, mine)) in standardized.iter().zip(&row).enumerate() {
                let from_the_dosage = (mine - mean) / divisor;
                assert!(
                    (from_the_dosage - theirs).abs() <= 1e-12 * largest,
                    "the variant {var} of the individual {position}: the dosage {mine} \
                     centered and divided is {from_the_dosage} and the standardized row \
                     holds {theirs}"
                );
            }
        }

        assert_eq!(
            agreed,
            [false, true, true, false],
            "the variants of the panel that have variance among the four tested"
        );

        // A variant of three alleles: one genotype of the individual 0 is
        // 0/0/2 and the others hold 0 and 1, and both refuse it with its
        // position among the variants the reader gave and the count.
        let block = a_block(
            6,
            3,
            &[&[0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1]],
        );
        let mut mine = DosageScratch::of(6);
        let mut theirs = RowScratch::of(6);
        let mut row = [0.0; 6];
        let mut standardized = [0.0; 6];

        let refused = the_dosages_of_a_row(
            &block.gts,
            3,
            7,
            MultiallelicVariants::Refused,
            &mut mine,
            &mut row,
        );
        let refused_there =
            the_standardized_row(&block.gts, 3, 7, &options, &mut theirs, &mut standardized);

        match (refused, refused_there) {
            (
                Err(Error::VariantWithMoreThanTwoAlleles {
                    position,
                    num_alleles,
                }),
                Err(Error::VariantWithMoreThanTwoAlleles {
                    position: there,
                    num_alleles: alleles_there,
                }),
            ) => {
                assert_eq!((position, num_alleles), (7, 3));
                assert_eq!((there, alleles_there), (7, 3));
            }
            (mine, theirs) => panic!(
                "the variant has three alleles, and the dosages gave {mine:?} and the \
                 standardized row gave {theirs:?}"
            ),
        }
    }

    /// A block that holds no genotype is refused, and so is a block of no
    /// variant, which is a reader that has stopped without saying so.
    #[test]
    fn a_block_that_holds_no_genotype_and_a_block_of_no_variant_are_refused() {
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        let design = the_design_of(&study, 8);
        let mut dosages = GwasDosages::of_a_study();
        let mut block = a_block(8, 2, &OF_EIGHT);
        block.gts = Vec::new();

        match dosages.read_the_block(&mut block, &design, the_first_block_of(2)) {
            Err(Error::FieldsNotInTheBlock { fields }) => {
                assert_eq!(fields, crate::variant::Needs::GTS);
            }
            other => panic!("the block holds no genotype, and that gave {other:?}"),
        }

        let mut block = a_block(8, 2, &[]);

        match dosages.read_the_block(&mut block, &design, the_first_block_of(2)) {
            Err(Error::ReaderGaveABlockOfNoVariants) => {}
            other => panic!("the block holds no variant, and that gave {other:?}"),
        }
    }

    /// A block whose arrays are not of the size it states is refused when
    /// every individual is tested, which is when nothing else looks at it:
    /// the rows are cut out of the genotypes by that size, so a block that
    /// is one allele short is read as three variants where it says four,
    /// and the fourth would leave the result with no error.
    ///
    /// A study of some of the individuals reaches the same refusal through
    /// `Block::retain_individuals`, which checks the block itself before
    /// it moves an allele. A study of all of them cuts nothing down, so
    /// this is the case that says the check is made here as well.
    #[test]
    fn a_block_whose_arrays_are_not_of_its_size_is_refused_with_every_individual_tested() {
        let (phenotype, design) = the_phenotype_and_the_design_of(&THE_PANEL_OF_EIGHT);
        let study = a_study(&phenotype, &design, &THE_PANEL_OF_EIGHT);
        let design = the_design_of(&study, 8);
        let mut dosages = GwasDosages::of_a_study();
        let mut block = a_block(8, 2, &OF_EIGHT);
        block.gts.pop();

        match dosages.read_the_block(&mut block, &design, the_first_block_of(2)) {
            Err(Error::BlockArrayOfAnotherSize {
                array,
                found,
                expected,
            }) => {
                assert_eq!((array, found, expected), ("gts", 63, 64));
            }
            other => panic!("the block is one allele short, and that gave {other:?}"),
        }
    }

    /// A block of another ploidy or of other individuals than the reader
    /// says its source has is refused, naming both shapes.
    ///
    /// The ploidy is the one no other check catches. The rows of a block
    /// are cut at its own individuals times its own ploidy, and the buffer
    /// they are read into is sized at the pass's ploidy, so when the two
    /// disagree the rows and the buffer do not line up. Both of these were
    /// run against the commit before this one on 23 September 2026: the
    /// haploid variant of five individuals below, read at the ploidy 2,
    /// gave `Ok` with every variant of the block gone and nothing to show
    /// it, because the buffer came out shorter than one row and the two
    /// zipped to nothing, so no row ran at all; the diploid block of eight
    /// read at the ploidy 6 got as far as cutting a row and stopped there,
    /// with an error about genotypes that are not whole, which names the
    /// genotypes of the block and not the ploidy the pass was reading at.
    /// Which of the two a block gets turns on how far the sizes are apart
    /// and on how many individuals are tested, and neither is an answer.
    /// The frequency of a variant is its mean dosage over the pass's
    /// ploidy, so even a block that came out whole would be read into
    /// frequencies of nothing.
    #[test]
    fn a_block_of_another_ploidy_or_of_other_individuals_than_the_readers_is_refused() {
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        let of_eight = the_design_of(&study, 8);
        let mut dosages = GwasDosages::of_a_study();
        let mut block = a_block(8, 2, &OF_EIGHT);

        match dosages.read_the_block(&mut block, &of_eight, the_first_block_of(6)) {
            Err(Error::BlocksDoNotFitTogether {
                num_individuals,
                ploidy,
                found_num_individuals,
                found_ploidy,
            }) => {
                assert_eq!((num_individuals, ploidy), (8, 6));
                assert_eq!((found_num_individuals, found_ploidy), (8, 2));
            }
            other => panic!("the block is diploid and the pass reads 6, and that gave {other:?}"),
        }

        // One haploid variant of five individuals, four of which are
        // tested, read by a pass that says its source is diploid.
        let four = [0, 1, 3, 4];
        let (phenotype, design) = the_phenotype_and_the_design_of(&four);
        let study = a_study(&phenotype, &design, &four);
        let of_five = the_design_of(&study, 5);
        let mut block = a_block(5, 1, &[&[0, 1, 1, 1, 0]]);

        match dosages.read_the_block(&mut block, &of_five, the_first_block_of(2)) {
            Err(Error::BlocksDoNotFitTogether {
                ploidy,
                found_ploidy,
                ..
            }) => {
                assert_eq!((ploidy, found_ploidy), (2, 1));
            }
            other => panic!("the block is haploid and the pass reads 2, and that gave {other:?}"),
        }

        // The same four individuals of a source the reader says has nine,
        // which is a block of others: the positions of the tested
        // individuals are positions among the reader's.
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        let of_nine = the_design_of(&study, 9);
        let mut block = a_block(8, 2, &OF_EIGHT);

        match dosages.read_the_block(&mut block, &of_nine, the_first_block_of(2)) {
            Err(Error::BlocksDoNotFitTogether {
                num_individuals,
                found_num_individuals,
                ..
            }) => {
                assert_eq!((num_individuals, found_num_individuals), (9, 8));
            }
            other => panic!(
                "the reader says its source has nine individuals and the block has \
                 eight, and that gave {other:?}"
            ),
        }
    }
}

/// Which model a study fits and which test it makes of every variant, as
/// "What it gives" and "Which individuals are tested, and the design" of
/// `docs/specs/gwas.md` state them.
#[cfg(test)]
mod choice {
    use super::fixtures::{TESTED_OF_EIGHT, a_study, the_phenotype_and_the_design_of};
    use super::{GwasModel, TestType, TraitType, the_model_and_the_test};
    use crate::error::Error;

    /// The kinship of the four tested individuals of the fixtures, which
    /// is what turns a model into a mixed one. Nothing here reads its
    /// values: the choice of the model is made on whether there is one.
    const A_KINSHIP: [f64; 16] = [
        1.0, 0.1, 0.0, 0.0, 0.1, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.2, 0.0, 0.0, 0.2, 1.0,
    ];

    /// What each study of these tests is asked with: the four individuals
    /// of the panel of eight that have a phenotype, so that no fixture
    /// here covers a whole panel either.
    fn the_model_and_the_test_of(
        trait_type: TraitType,
        kinship: Option<&[f64]>,
        test: Option<TestType>,
    ) -> crate::error::Result<(GwasModel, TestType)> {
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let mut study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        study.trait_type = trait_type;
        study.kinship = kinship;
        study.test = test;
        the_model_and_the_test(&study)
    }

    /// A continuous trait with no kinship is a linear model, and its one
    /// test is the Wald test, the t test of the effect it fitted.
    #[test]
    fn a_continuous_trait_with_no_kinship_is_a_linear_model_tested_by_its_t_test() {
        assert_eq!(
            the_model_and_the_test_of(TraitType::Continuous, None, None)
                .expect("the model of a continuous trait with no kinship"),
            (GwasModel::Lm, TestType::Wald)
        );
        assert_eq!(
            the_model_and_the_test_of(TraitType::Continuous, None, Some(TestType::Wald))
                .expect("the test it has, asked for by name"),
            (GwasModel::Lm, TestType::Wald)
        );
    }

    /// A continuous trait with a kinship is a linear mixed model, whose
    /// default is the Wald test and which takes the score test too.
    #[test]
    fn a_continuous_trait_with_a_kinship_is_a_linear_mixed_model_whose_default_is_the_wald_test() {
        assert_eq!(
            the_model_and_the_test_of(TraitType::Continuous, Some(&A_KINSHIP), None)
                .expect("the model of a continuous trait with a kinship"),
            (GwasModel::Lmm, TestType::Wald)
        );
        assert_eq!(
            the_model_and_the_test_of(
                TraitType::Continuous,
                Some(&A_KINSHIP),
                Some(TestType::Score)
            )
            .expect("the score test of a linear mixed model"),
            (GwasModel::Lmm, TestType::Score)
        );
    }

    /// A binomial trait with no kinship is a logistic regression, whose
    /// default is the Wald test and which takes the score test too.
    #[test]
    fn a_binomial_trait_with_no_kinship_is_a_logistic_regression_whose_default_is_the_wald_test() {
        assert_eq!(
            the_model_and_the_test_of(TraitType::Binomial, None, None)
                .expect("the model of a binomial trait with no kinship"),
            (GwasModel::Glm, TestType::Wald)
        );
        assert_eq!(
            the_model_and_the_test_of(TraitType::Binomial, None, Some(TestType::Score))
                .expect("the score test of a logistic regression"),
            (GwasModel::Glm, TestType::Score)
        );
    }

    /// A binomial trait with a kinship is a logistic mixed model, whose
    /// one test is the score test: it costs no fit per variant.
    #[test]
    fn a_binomial_trait_with_a_kinship_is_a_logistic_mixed_model_whose_default_is_the_score_test() {
        assert_eq!(
            the_model_and_the_test_of(TraitType::Binomial, Some(&A_KINSHIP), None)
                .expect("the model of a binomial trait with a kinship"),
            (GwasModel::Glmm, TestType::Score)
        );
        assert_eq!(
            the_model_and_the_test_of(TraitType::Binomial, Some(&A_KINSHIP), Some(TestType::Score))
                .expect("the test it has, asked for by name"),
            (GwasModel::Glmm, TestType::Score)
        );
    }

    /// The score test of a continuous trait with no kinship is refused:
    /// the only test of a linear model is the t test of the effect it
    /// fitted. The same trait with a kinship takes the score test, so the
    /// refusal is the pair and not the test.
    #[test]
    fn the_score_test_of_a_continuous_trait_with_no_kinship_is_refused() {
        match the_model_and_the_test_of(TraitType::Continuous, None, Some(TestType::Score)) {
            Err(Error::GwasScoreTestOfALinearModel) => {}
            other => panic!(
                "a continuous trait with no kinship is a linear model, whose only test is \
                 its t test, and that gave {other:?}"
            ),
        }
    }

    /// The Wald test of a binomial trait with a kinship is refused: it
    /// would fit one logistic mixed model for every variant. The same
    /// trait without a kinship takes the Wald test, so the refusal is the
    /// pair and not the test.
    #[test]
    fn the_wald_test_of_a_binomial_trait_with_a_kinship_is_refused() {
        match the_model_and_the_test_of(TraitType::Binomial, Some(&A_KINSHIP), Some(TestType::Wald))
        {
            Err(Error::GwasWaldTestOfALogisticMixedModel) => {}
            other => panic!(
                "a binomial trait with a kinship is a logistic mixed model, and a Wald \
                 test of it would fit one model for every variant, and that gave {other:?}"
            ),
        }
    }
}

/// The rows a study gives back, and the variants that have no answer, as
/// "The variants that have no answer" of `docs/specs/gwas.md` states them.
#[cfg(test)]
mod result {
    use super::fixtures::{
        FREQUENCIES_OF_THE_TESTED, OF_EIGHT, TESTED_OF_EIGHT, a_block, a_study, assert_the_values,
        the_design_of, the_first_block_of, the_phenotype_and_the_design_of,
    };
    use super::{
        Answers, GrammarGammaApprox, Gwas, GwasDosages, NullModel, the_model_and_the_test,
    };
    use crate::error::Error;
    use crate::variant::ChromTable;

    /// The dosages of the panel of eight over the four individuals that
    /// have a phenotype, of which the second and the third variants have
    /// variance and the first and the fourth have none.
    fn the_dosages_of_the_tested() -> GwasDosages {
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        let design = the_design_of(&study, 8);
        let mut block = a_block(8, 2, &OF_EIGHT);
        let mut dosages = GwasDosages::of_a_study();
        dosages
            .read_the_block(&mut block, &design, the_first_block_of(2))
            .expect("the dosages of the block");
        dosages
    }

    /// A result opened on the null model of that study, with the model and
    /// the test the study chose.
    fn a_result() -> Gwas {
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        let (model, test) = the_model_and_the_test(&study).expect("the model and the test");
        Gwas::of_the_null_model(
            NullModel {
                model,
                test,
                covariate_effects: vec![3.6666666666666683, 1.0],
                residual_variance: Some(3.333333333333334),
                genetic_variance: None,
                heritability: None,
                num_individuals: 4,
            },
            GrammarGammaApprox::NotUsed,
            ChromTable::new(),
        )
    }

    /// The result carries whether the study made the GRAMMAR-Gamma
    /// approximation, which is the one thing about a study that is not in
    /// its null model, and a study that was asked for nothing makes it
    /// not.
    ///
    /// It is a `bool` in the result, which "The Rust interface" of
    /// `docs/specs/gwas.md` fixes, and the two named values on the way in,
    /// so that a call cannot say which it means with a bare `false`.
    #[test]
    fn the_result_says_whether_the_approximation_was_used() {
        let null_model = a_result().null_model;
        for (asked, used) in [
            (GrammarGammaApprox::Used, true),
            (GrammarGammaApprox::NotUsed, false),
        ] {
            let result = Gwas::of_the_null_model(null_model.clone(), asked, ChromTable::new());
            assert_eq!(result.used_grammar_gamma_approx, used);
            assert_eq!(result.num_vars, 0, "no block has been read");
            assert!(result.allele_freq.is_empty());
        }
        let asked_for_nothing = match super::DEFAULT_USE_GRAMMAR_GAMMA_APPROX {
            true => GrammarGammaApprox::Used,
            false => GrammarGammaApprox::NotUsed,
        };
        let result = Gwas::of_the_null_model(null_model, asked_for_nothing, ChromTable::new());
        assert!(
            !result.used_grammar_gamma_approx,
            "a study that is asked for nothing does not approximate"
        );
    }

    /// A variant with no variance among the tested individuals keeps its
    /// frequency and has NaN for its effect, its uncertainty and its
    /// p-value, and the variants that do have variance take the answers of
    /// the model in the order of the block.
    ///
    /// The fixture is the panel of eight of `fixtures::OF_EIGHT` over its
    /// four tested individuals: its first variant is heterozygous in all
    /// four and its fourth has no genotype in any of them, so the answers
    /// the model gives are two and the rows are four. The second block
    /// added is the same one again, which is what says that a block adds
    /// its rows to the ones that are there and does not replace them.
    #[test]
    fn a_variant_with_no_variance_keeps_its_frequency_and_has_three_nans() {
        let dosages = the_dosages_of_the_tested();
        let mut result = a_result();
        let beta = [1.5, -0.25];
        let se = [0.600925212577332, 0.4];
        let p_value = [0.088004892382756, 0.5];

        result
            .add_the_block(
                &dosages,
                Answers {
                    beta: &beta,
                    se: &se,
                    p_value: &p_value,
                },
            )
            .expect("the answers of the two variants that have variance");

        assert_eq!(result.num_vars, 4);
        assert_the_values(
            &result.allele_freq,
            &FREQUENCIES_OF_THE_TESTED,
            "the frequency of every variant of the block",
        );
        for (column, values) in [
            ("beta", &result.beta),
            ("se", &result.se),
            ("p_value", &result.p_value),
        ] {
            assert!(
                values.first().is_some_and(|value| value.is_nan()),
                "the first variant is heterozygous in all four and has no {column}"
            );
            assert!(
                values.get(3).is_some_and(|value| value.is_nan()),
                "the fourth variant has no genotype among the four and no {column}"
            );
        }
        assert_the_values(
            result.beta.get(1..3).unwrap_or_default(),
            &beta,
            "the effects of the two variants that have variance",
        );
        assert_the_values(
            result.se.get(1..3).unwrap_or_default(),
            &se,
            "their standard errors",
        );
        assert_the_values(
            result.p_value.get(1..3).unwrap_or_default(),
            &p_value,
            "their p-values",
        );

        result
            .add_the_block(
                &dosages,
                Answers {
                    beta: &beta,
                    se: &se,
                    p_value: &p_value,
                },
            )
            .expect("the same block again");

        assert_eq!(result.num_vars, 8, "the rows of the two blocks");
        assert_eq!(result.allele_freq.len(), 8);
        assert_eq!(result.p_value.len(), 8);
    }

    /// A model that answered for another number of variants than the block
    /// has with variance is refused, naming the column and the two counts.
    /// It is a defect of popnei and not of a user, and it would otherwise
    /// put the answer of one variant in the row of another from the first
    /// missing value on.
    #[test]
    fn a_model_that_answered_for_another_number_of_variants_is_refused() {
        let dosages = the_dosages_of_the_tested();
        let mut result = a_result();

        match result.add_the_block(
            &dosages,
            Answers {
                beta: &[1.5],
                se: &[0.6, 0.4],
                p_value: &[0.08, 0.5],
            },
        ) {
            Err(Error::GwasAnswersOfAnotherSize {
                column,
                num_values,
                num_with_variance,
            }) => {
                assert_eq!((column, num_values, num_with_variance), ("beta", 1, 2));
            }
            other => panic!("the model answered one effect of two, and that gave {other:?}"),
        }

        match result.add_the_block(
            &dosages,
            Answers {
                beta: &[1.5, -0.25],
                se: &[0.6, 0.4, 0.3],
                p_value: &[0.08, 0.5],
            },
        ) {
            Err(Error::GwasAnswersOfAnotherSize {
                column,
                num_values,
                num_with_variance,
            }) => {
                assert_eq!((column, num_values, num_with_variance), ("se", 3, 2));
            }
            other => {
                panic!("the model answered three standard errors of two, and that gave {other:?}")
            }
        }

        assert_eq!(result.num_vars, 0, "no row was added by either");
        assert!(result.allele_freq.is_empty());
    }
}

/// The linear model against the worked example of `docs/specs/gwas.md`,
/// whose numbers are pyNei's, and against plink2 on the panel with every
/// genotype called.
#[cfg(test)]
mod lm {
    use std::collections::HashMap;
    use std::io::Cursor;
    use std::path::{Path, PathBuf};

    use super::{
        BlockReader, Design, Gwas, GwasInput, GwasModel, LinearModel, TestType, TraitType,
        calc_gwas,
    };
    use crate::error::{Error, Result};
    use crate::io::vcf::{VcfOptions, VcfReader};

    /// How far a number of the worked example may be from pyNei's, as a
    /// share of that number: 3e-15.
    ///
    /// "The worked example" of `docs/specs/gwas.md` asks for 1e-12, and
    /// "How it is verified" of "What every model shares" asks that such a
    /// number be lowered until it fails and then set two or three times
    /// above where it broke. It was, on 23 September 2026 on both
    /// backends: it breaks at 1e-15, where the `beta` of `v1` on faer is
    /// 1.24e-15 away, and the worst difference of the thirteen numbers it
    /// holds, the four of the null model and the nine of the three rows,
    /// is that 1.24e-15 on faer and 7.4e-16 on Accelerate, the `se` of
    /// `v0`. So this is 2.4 times the worst measured, and the spec's 1e-12
    /// had about 800 times it.
    ///
    /// No number of the example is near 0: the smallest is the `beta` of
    /// `v1`, 0.3125, against an `se` of 1.3, so a share of the number and
    /// a share of the scale of what is estimated are the same bound here,
    /// and the rule of the spec that a `beta` is measured against its `se`
    /// changes nothing. The example runs on whole dosages and a trait of
    /// whole numbers, and it reads no reference file, so nothing of it is
    /// rounded away before the comparison.
    const OF_THE_WORKED_EXAMPLE: f64 = 3e-15;

    /// How far a `beta` or an `se` of the panel may be from plink2's, as a
    /// share of the `se` plink2 printed for that variant.
    ///
    /// It is the 1e-5 of "How it is verified" of "The linear model" of
    /// `docs/specs/gwas.md`, which is a share of the `se` and not of the
    /// value, because a study is mostly null and a `beta` that cancelled
    /// to near 0 is no guide to its own error. plink2 prints six
    /// significant digits, so the number it is compared against is itself
    /// rounded by up to 5e-7 absolute, which of the 1.2e-6 this bound
    /// comes to at the smallest `se` of the six is 41 per cent. The
    /// arithmetic gets what is left, and what it took is in the doc
    /// comment of the test.
    ///
    /// This one is not lowered, unlike [`OF_THE_WORKED_EXAMPLE`]: it is a
    /// bound on what plink2 printed and not on what popnei computed, and
    /// tightening it would measure the rounding of six digits.
    const OF_PLINK2: f64 = 1e-5;

    /// How far a p-value of the panel may be from plink2's, as a share of
    /// it. plink2 prints it to six significant digits like the rest, so
    /// the printing alone can take half of this.
    const OF_PLINK2_P_VALUE: f64 = 1e-5;

    /// How far an `allele_freq` of the panel may be from plink2's
    /// `A1_FREQ`: 1e-6 absolute, since it is a frequency and lies between
    /// 0 and 1.
    const OF_PLINK2_FREQUENCY: f64 = 1e-6;

    /// The path of one of the files of `tests/reference/gwas/`, where the
    /// trait, the covariates and plink2's answers are.
    fn the_reference_path(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/gwas")
            .join(name)
    }

    /// The panel with every genotype called, 200 individuals and 1200
    /// biallelic diploid variants, which is the panel plink2 was run on.
    fn the_panel_path() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/kinship/panel_called.vcf.gz")
    }

    /// The study of the variants of `reader`, with no second pass for the
    /// GRAMMAR-Gamma approximation, which a linear model does not make.
    fn the_study_of<R: BlockReader>(reader: &mut R, input: &GwasInput<'_>) -> Result<Gwas> {
        calc_gwas(reader, None::<&mut R>, input)
    }

    /// The worked example of `docs/specs/gwas.md` as a VCF: six diploid
    /// individuals, `i0` to `i5`, and three variants, of which the third
    /// has every individual heterozygous and so no variance. The `./.` of
    /// `i3` at `v1` is the genotype that takes the mean dosage of its
    /// variant, 0.8.
    fn the_worked_example_vcf() -> Vec<u8> {
        let mut vcf = String::from(THE_HEADER);
        for (var, genotypes) in [
            ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1"],
            ["0/0", "0/1", "1/1", "./.", "0/1", "0/0"],
            ["0/1", "0/1", "0/1", "0/1", "0/1", "0/1"],
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
        vcf.into_bytes()
    }

    /// The header of that VCF, six diploid individuals and no variant.
    const THE_HEADER: &str = "##fileformat=VCFv4.2\n\
        ##contig=<ID=1>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti0\ti1\ti2\ti3\ti4\ti5\n";

    /// The trait of the worked example, one value for each of the six
    /// individuals.
    const THE_TRAIT: [f64; 6] = [2.0, 3.0, 5.0, 4.0, 4.0, 7.0];

    /// Its design: the intercept and the one covariate, row after row.
    const THE_DESIGN: [f64; 12] = [
        1.0, 0.0, //
        1.0, 1.0, //
        1.0, 0.0, //
        1.0, 1.0, //
        1.0, 0.0, //
        1.0, 1.0,
    ];

    /// The six individuals of the worked example, which are every
    /// individual of its VCF.
    const THE_INDIVIDUALS: [usize; 6] = [0, 1, 2, 3, 4, 5];

    /// A reader over the bytes of a VCF of diploid individuals.
    fn reader_over(vcf: &[u8]) -> VcfReader<Cursor<Vec<u8>>> {
        match VcfReader::new(Cursor::new(vcf.to_vec()), VcfOptions::default()) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        }
    }

    /// The study of the worked example: a continuous trait, one covariate,
    /// no kinship and the default test.
    fn the_worked_example_study() -> GwasInput<'static> {
        GwasInput {
            phenotype: &THE_TRAIT,
            trait_type: TraitType::Continuous,
            design: &THE_DESIGN,
            num_coefs: 2,
            kinship: None,
            test: None,
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS,
            transform_to_biallelic: false,
        }
    }

    /// `found` within `tolerance` of `expected` as a share of `expected`.
    fn assert_within(found: f64, expected: f64, tolerance: f64, what: &str) {
        let difference = (found - expected).abs();
        assert!(
            difference <= tolerance * expected.abs(),
            "{what} is {found} and the literal is {expected}, {difference} away, which is \
             {share} of it against the {tolerance} allowed",
            share = difference / expected.abs()
        );
    }

    /// `found` within `tolerance` times `scale`, which is the `se` of the
    /// variant for a `beta` and for an `se`.
    fn assert_within_the_scale(found: f64, expected: f64, scale: f64, tolerance: f64, what: &str) {
        let difference = (found - expected).abs();
        assert!(
            difference <= tolerance * scale,
            "{what} is {found} and plink2 gives {expected}, {difference} away, which is \
             {share} of the {scale} it is uncertain by against the {tolerance} allowed",
            share = difference / scale
        );
    }

    /// The null model of the worked example is the one pyNei gives at
    /// commit ef0ca6e: an intercept of 3.6666666666666683 and a covariate
    /// effect of 1.0, a residual sum of squares of 13.333333333333336 over
    /// its 4 degrees of freedom, and so a residual variance of
    /// 3.333333333333334.
    ///
    /// It is fitted here and not through `calc_gwas`, because the result
    /// of a study carries the residual variance and not the sum of squares
    /// the spec's table gives beside it.
    #[test]
    fn the_null_of_the_worked_example_is_pyneis_coefficients_and_sum_of_squares() {
        let study = the_worked_example_study();
        let design = match Design::of_the_study(&study, THE_INDIVIDUALS.len()) {
            Ok(design) => design,
            Err(error) => panic!("the design of the worked example: {error}"),
        };
        let fitted = match LinearModel::of_the_study(&THE_TRAIT, &design) {
            Ok(fitted) => fitted,
            Err(error) => panic!("the null model of the worked example: {error}"),
        };
        let null = fitted.null_model(TestType::Wald);
        assert_eq!(null.model, GwasModel::Lm);
        assert_eq!(null.test, TestType::Wald);
        assert_eq!(null.num_individuals, 6);
        assert_eq!(null.covariate_effects.len(), 2);
        assert_within(
            null.covariate_effects[0],
            3.666_666_666_666_668_3,
            OF_THE_WORKED_EXAMPLE,
            "the intercept",
        );
        assert_within(
            null.covariate_effects[1],
            1.0,
            OF_THE_WORKED_EXAMPLE,
            "the effect of the covariate",
        );
        assert_within(
            fitted.residual_sum_of_squares(),
            13.333_333_333_333_336,
            OF_THE_WORKED_EXAMPLE,
            "the residual sum of squares",
        );
        let residual_variance = null
            .residual_variance
            .expect("the residual variance of a continuous trait");
        assert_within(
            residual_variance,
            3.333_333_333_333_334,
            OF_THE_WORKED_EXAMPLE,
            "the residual variance",
        );
        assert_eq!(null.genetic_variance, None, "a model with no kinship");
        assert_eq!(null.heritability, None, "a model with no kinship");
    }

    /// The three rows of the worked example are pyNei's: `v0` and `v1` get
    /// the frequency, the effect, the standard error and the p-value of
    /// the spec's table, and `v2`, where every individual is
    /// heterozygous, keeps its frequency of 0.5 and has three NaNs.
    ///
    /// This is the one test of the module that needs no reference program
    /// and no reference file. The three columns of the variants come out
    /// of it as well, since the example is read from a VCF that carries
    /// them.
    #[test]
    fn the_worked_example_gives_pyneis_three_rows_and_no_answer_for_the_third() {
        let vcf = the_worked_example_vcf();
        let mut reader = reader_over(&vcf);
        let study = the_worked_example_study();
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of the worked example: {error}"),
        };
        assert_eq!(result.num_vars, 3);
        assert_eq!(
            result.ids.as_deref(),
            Some(["v0".to_owned(), "v1".to_owned(), "v2".to_owned()].as_slice()),
            "the ids of the variants"
        );
        assert_eq!(
            result.poss.as_deref(),
            Some([1000_u64, 2000, 3000].as_slice()),
            "the positions of the variants"
        );
        assert_eq!(
            result.chroms.as_deref(),
            Some([0_u32, 0, 0].as_slice()),
            "the chromosomes of the variants"
        );
        assert_eq!(result.chrom_table.name(0), Some("1"));
        assert!(!result.used_grammar_gamma_approx);
        assert_eq!(result.null_model.model, GwasModel::Lm);
        for (var, (allele_freq, beta, se, p_value)) in [
            (0.5, 1.5, 0.600_925_212_577_332, 0.088_004_892_382_756),
            (0.4, 0.3125, 1.305_204_592_306_424, 0.826_200_867_452_417),
        ]
        .into_iter()
        .enumerate()
        {
            assert_within(
                result.allele_freq[var],
                allele_freq,
                OF_THE_WORKED_EXAMPLE,
                &format!("the frequency of v{var}"),
            );
            assert_within(
                result.beta[var],
                beta,
                OF_THE_WORKED_EXAMPLE,
                &format!("the effect of v{var}"),
            );
            assert_within(
                result.se[var],
                se,
                OF_THE_WORKED_EXAMPLE,
                &format!("the standard error of v{var}"),
            );
            assert_within(
                result.p_value[var],
                p_value,
                OF_THE_WORKED_EXAMPLE,
                &format!("the p-value of v{var}"),
            );
        }
        assert_within(
            result.allele_freq[2],
            0.5,
            OF_THE_WORKED_EXAMPLE,
            "the frequency of v2",
        );
        assert!(result.beta[2].is_nan(), "the effect of v2");
        assert!(result.se[2].is_nan(), "the standard error of v2");
        assert!(result.p_value[2].is_nan(), "the p-value of v2");
    }

    /// How many variants a study of six individuals reads in one block,
    /// which is [`crate::block::MAX_NUM_VARS_PER_BLOCK`]: the pass puts a
    /// `Reblock` before it, so the blocks it sees are of that size
    /// whatever size the reader gives, and a source has to pass it for the
    /// pass to read a second block.
    const VARS_OF_ONE_BLOCK: usize = 10_000;

    /// A study over more variants than one block holds gives the same
    /// answers in its second block as in its first, and one row for every
    /// variant in the order the reader gave them.
    ///
    /// The variants are the three of the worked example over and over, so
    /// every variant of the study has the answer the spec's table gives
    /// for its pattern, whichever block it fell in: the one with no
    /// variance keeps its frequency and has three NaNs wherever it is.
    /// What this covers that the worked example does not is the buffers of
    /// the model and of the dosages being reused, the rows of the second
    /// block being added after the first's and not over them, and the
    /// three columns of the variants growing across the blocks.
    ///
    /// The second block is on another chromosome than the first, so that
    /// the column of chromosomes is asserted where it is built: the
    /// numbers a block holds are read through the table of the reader, and
    /// a study that kept the first block's name for every variant after it
    /// would answer every question this test asks about the ids and the
    /// positions and still be wrong about where a variant is.
    #[test]
    fn a_study_of_more_variants_than_one_block_answers_the_same_in_every_block() {
        let patterns = [
            ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1"],
            ["0/0", "0/1", "1/1", "./.", "0/1", "0/0"],
            ["0/1", "0/1", "0/1", "0/1", "0/1", "0/1"],
        ];
        let num_vars = VARS_OF_ONE_BLOCK.saturating_add(100);
        let mut vcf = String::from(THE_HEADER);
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
        let mut reader = reader_over(vcf.as_bytes());
        let study = the_worked_example_study();
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of {num_vars} variants: {error}"),
        };
        assert_eq!(result.num_vars, num_vars, "the variants of the study");
        assert!(
            num_vars > VARS_OF_ONE_BLOCK,
            "the study has to read more than one block"
        );
        let ids = result.ids.as_deref().expect("the ids of the variants");
        let poss = result.poss.as_deref().expect("the positions");
        let chroms = result.chroms.as_deref().expect("the chromosomes");
        assert_eq!(ids.len(), num_vars, "one id for each variant");
        assert_eq!(poss.len(), num_vars, "one position for each variant");
        assert_eq!(chroms.len(), num_vars, "one chromosome for each variant");
        for var in [0, VARS_OF_ONE_BLOCK, num_vars.saturating_sub(1)] {
            assert_eq!(ids[var], format!("v{var}"), "the id of the variant {var}");
            assert_eq!(
                poss[var],
                var.saturating_add(1).saturating_mul(10) as u64,
                "the position of the variant {var}"
            );
            let name = result.chrom_table.name(chroms[var]);
            let expected = match var < VARS_OF_ONE_BLOCK {
                true => "1",
                false => "2",
            };
            assert_eq!(
                name,
                Some(expected),
                "the chromosome of the variant {var}, which the table of the pass names"
            );
        }
        for (var, (allele_freq, beta, se, p_value)) in (0..num_vars)
            .map(|var| {
                (
                    var,
                    match var % patterns.len() {
                        0 => (0.5, 1.5, 0.600_925_212_577_332, 0.088_004_892_382_756),
                        1 => (0.4, 0.3125, 1.305_204_592_306_424, 0.826_200_867_452_417),
                        _ => (0.5, f64::NAN, f64::NAN, f64::NAN),
                    },
                )
            })
            .collect::<Vec<_>>()
        {
            assert_within(
                result.allele_freq[var],
                allele_freq,
                OF_THE_WORKED_EXAMPLE,
                &format!("the frequency of the variant {var}"),
            );
            if beta.is_nan() {
                assert!(result.beta[var].is_nan(), "the effect of the variant {var}");
                assert!(result.se[var].is_nan(), "the error of the variant {var}");
                assert!(
                    result.p_value[var].is_nan(),
                    "the p-value of the variant {var}"
                );
                continue;
            }
            assert_within(
                result.beta[var],
                beta,
                OF_THE_WORKED_EXAMPLE,
                &format!("the effect of the variant {var}"),
            );
            assert_within(
                result.se[var],
                se,
                OF_THE_WORKED_EXAMPLE,
                &format!("the error of the variant {var}"),
            );
            assert_within(
                result.p_value[var],
                p_value,
                OF_THE_WORKED_EXAMPLE,
                &format!("the p-value of the variant {var}"),
            );
        }
    }

    /// The trait and the two covariates of the panel, read from
    /// `tests/reference/gwas/phenotypes.csv`, in the order `individuals`
    /// has the individuals, which is the order the VCF has them.
    ///
    /// The file holds one line per individual with its name, the
    /// continuous trait, the binomial one, the two covariates and the
    /// subpopulation. What is taken here is the continuous trait and the
    /// two covariates, which is what plink2 was given.
    fn the_trait_and_the_design_of_the_panel(individuals: &[String]) -> (Vec<f64>, Vec<f64>) {
        let path = the_reference_path("phenotypes.csv");
        let text = match std::fs::read_to_string(&path) {
            Ok(text) => text,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let mut lines = text.lines();
        let header: Vec<&str> = lines
            .next()
            .expect("the header of phenotypes.csv")
            .split(',')
            .collect();
        let column_of = |name: &str| match header.iter().position(|held| *held == name) {
            Some(at) => at,
            None => panic!("phenotypes.csv has no column {name}"),
        };
        let of_the_name = column_of("IID");
        let of_the_trait = column_of("cont");
        let of_cov1 = column_of("cov1");
        let of_cov2 = column_of("cov2");
        let mut rows: HashMap<String, [f64; 3]> = HashMap::new();
        for line in lines.filter(|line| !line.is_empty()) {
            let fields: Vec<&str> = line.split(',').collect();
            let number = |at: usize| match fields.get(at).map(|field| field.parse::<f64>()) {
                Some(Ok(number)) => number,
                _ => panic!("the column {at} of the line `{line}` of phenotypes.csv"),
            };
            let name = match fields.get(of_the_name) {
                Some(name) => (*name).to_owned(),
                None => panic!("the line `{line}` of phenotypes.csv has no name"),
            };
            rows.insert(
                name,
                [number(of_the_trait), number(of_cov1), number(of_cov2)],
            );
        }
        let mut phenotype = Vec::with_capacity(individuals.len());
        let mut design = Vec::with_capacity(individuals.len().saturating_mul(3));
        for individual in individuals {
            match rows.get(individual) {
                Some(row) => {
                    phenotype.push(row[0]);
                    design.extend_from_slice(&[1.0, row[1], row[2]]);
                }
                None => panic!("{individual} is not in phenotypes.csv"),
            }
        }
        (phenotype, design)
    }

    /// What plink2 v2.0.0-a.7.7 wrote for six variants of the panel with
    /// every genotype called, in
    /// `tests/reference/gwas/plink2.panel_called.glm.linear.tsv`, with
    /// `cov1` and `cov2` as covariates: the id and then `A1_FREQ`, `BETA`,
    /// `SE` and `P`, the four columns "How it is verified" of "The linear
    /// model" of `docs/specs/gwas.md` compares. Five of them are the
    /// causal variants of `causal_vars.csv` and `var0000` is not causal.
    const OF_PLINK2_SIX: [(&str, f64, f64, f64, f64); 6] = [
        ("var0000", 0.21, -0.424136, 0.139354, 0.00265846),
        ("var0052", 0.3225, -0.697724, 0.122348, 4.28981e-08),
        ("var0629", 0.175, -0.813852, 0.161809, 1.10646e-06),
        ("var0751", 0.41, -0.0963977, 0.130636, 0.461451),
        ("var1137", 0.2625, -0.171137, 0.145987, 0.242511),
        ("var1188", 0.4675, -0.655393, 0.138912, 4.51958e-06),
    ];

    /// The six variants of the panel are plink2's: `beta` and `se` within
    /// 1e-5 of the `se` plink2 printed for the variant, `p_value` within
    /// 1e-5 of itself and `allele_freq` within 1e-6 absolute, which is
    /// what the spec asks of the whole columns in Python.
    ///
    /// Measured over the six on 23 September 2026, the same on Accelerate
    /// and on faer to two digits: the worst `beta` or `se` is 3.7e-6 of
    /// the `se` of its variant, 37 per cent of the 1e-5 allowed, and it is
    /// the `se` of `var0052`; the worst p-value is 2.8e-6 of itself, 28
    /// per cent of its 1e-5; and every `allele_freq` is plink2's
    /// `A1_FREQ` exactly. plink2 prints six significant digits, so its own
    /// rounding of a `beta` or an `se` below 1 is up to 5e-7 absolute,
    /// which at the smallest `se` of the six, the 0.122348 of `var0052`,
    /// is 4.1e-6 of that `se`: the whole of what was measured fits inside
    /// what the printing alone can account for, and the arithmetic has
    /// spent nothing that can be seen here.
    ///
    /// Over all 1200 variants the same run gives 1.94e-5 of the `se`,
    /// which this bound would not hold, and the reason is the printing and
    /// not the arithmetic: `var0482` has a `beta` of 1.0389, and six
    /// significant digits of a value above 1 are rounded by up to 5e-6
    /// absolute, which is 2.6e-5 of its `se` of 0.19. The two variants of
    /// the panel that pass 1e-5 are the two whose `beta` passes 1. The
    /// check over the whole column is the Python one of work package 3 of
    /// `docs/plans/gwas-linear.md`, and it is where that has to be dealt
    /// with; none of the six variants here has a `beta` above 1.
    #[test]
    fn the_six_variants_of_the_panel_are_plink2s_effect_error_and_p_value() {
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
        let (phenotype, design) = the_trait_and_the_design_of_the_panel(&individuals);
        let tested: Vec<usize> = (0..individuals.len()).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Continuous,
            design: &design,
            num_coefs: 3,
            kinship: None,
            test: None,
            use_grammar_gamma_approx: false,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of the panel: {error}"),
        };
        assert_eq!(result.num_vars, 1200, "the variants of the panel");
        assert_eq!(result.null_model.num_individuals, 200);
        let ids = result.ids.as_deref().expect("the ids of the variants");
        for (id, allele_freq, beta, se, p_value) in OF_PLINK2_SIX {
            let var = match ids.iter().position(|held| held == id) {
                Some(var) => var,
                None => panic!("{id} is not a variant of the panel"),
            };
            let found = result.allele_freq[var];
            assert!(
                (found - allele_freq).abs() <= OF_PLINK2_FREQUENCY,
                "the frequency of {id} is {found} and plink2 gives {allele_freq}"
            );
            assert_within_the_scale(
                result.beta[var],
                beta,
                se,
                OF_PLINK2,
                &format!("the effect of {id}"),
            );
            assert_within_the_scale(
                result.se[var],
                se,
                se,
                OF_PLINK2,
                &format!("the standard error of {id}"),
            );
            assert_within(
                result.p_value[var],
                p_value,
                OF_PLINK2_P_VALUE,
                &format!("the p-value of {id}"),
            );
        }
    }

    /// The three models that are not written are refused, each naming the
    /// study that asked for it.
    #[test]
    fn a_model_that_is_not_written_yet_is_refused() {
        let vcf = the_worked_example_vcf();
        let kinship = [0.0_f64; 36];
        for (trait_type, with_a_kinship) in [
            (TraitType::Continuous, true),
            (TraitType::Binomial, false),
            (TraitType::Binomial, true),
        ] {
            let phenotype: Vec<f64> = match trait_type {
                TraitType::Continuous => THE_TRAIT.to_vec(),
                TraitType::Binomial => vec![0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
            };
            let study = GwasInput {
                phenotype: &phenotype,
                trait_type,
                design: &THE_DESIGN,
                num_coefs: 2,
                kinship: match with_a_kinship {
                    true => Some(kinship.as_slice()),
                    false => None,
                },
                test: None,
                use_grammar_gamma_approx: false,
                individuals: &THE_INDIVIDUALS,
                transform_to_biallelic: false,
            };
            let wanted = match (trait_type, with_a_kinship) {
                (TraitType::Continuous, _) => GwasModel::Lmm,
                (TraitType::Binomial, false) => GwasModel::Glm,
                (TraitType::Binomial, true) => GwasModel::Glmm,
            };
            let mut reader = reader_over(&vcf);
            match the_study_of(&mut reader, &study) {
                Err(Error::GwasModelNotBuilt { model }) => {
                    assert_eq!(
                        model, wanted,
                        "{trait_type:?} with a kinship of {with_a_kinship}"
                    );
                    let said = Error::GwasModelNotBuilt { model }.to_string();
                    assert!(
                        said.contains("being written"),
                        "the study was refused with {said}"
                    );
                }
                Err(error) => panic!("{trait_type:?}, kinship {with_a_kinship}: {error}"),
                Ok(_) => panic!("{trait_type:?} with a kinship of {with_a_kinship} was run"),
            }
        }
    }

    /// The GRAMMAR-Gamma approximation stands in for the denominator of a
    /// mixed model's test, so a study with no kinship is refused for
    /// asking for it, which is pyNei's refusal of the same pair.
    #[test]
    fn the_grammar_gamma_approximation_without_a_kinship_is_refused() {
        let vcf = the_worked_example_vcf();
        let mut reader = reader_over(&vcf);
        let study = GwasInput {
            use_grammar_gamma_approx: true,
            ..the_worked_example_study()
        };
        match the_study_of(&mut reader, &study) {
            Err(Error::GwasGrammarGammaWithoutAKinship) => {}
            Err(error) => panic!("the study was refused with {error}"),
            Ok(_) => panic!("the approximation was made without a kinship"),
        }
    }

    /// The trait of the six individuals for a variant that explains all
    /// but `delta` of what the null model left, as literals.
    ///
    /// The trait is `1 + 2 * cov + 3 * dosage + delta * u` over the
    /// covariate `0 1 0 1 0 1` and the dosages `0 1 2 0 1 2` of `v0`, with
    /// `u` the direction `1 0 -1 -1 0 1`, which is at right angles to the
    /// column of ones, to the covariate and to the dosages: `u` sums to 0,
    /// its values at the three individuals of `cov = 1` sum to 0, and its
    /// product with the dosages is `-2 + 2`. So the effect of the variant
    /// is 3 whatever `delta` is, and what the variant leaves unexplained
    /// is `delta² * 4`, `4` being the squared length of `u`.
    ///
    /// The values are written out and not computed, so that the test reads
    /// the same numbers numpy was given: `repr` of a float prints the
    /// shortest text that reads back as the same number, and this is what
    /// numpy 2.5.3 printed for `1 + 2 * cov + 3 * x + delta * u`.
    const THE_TRAIT_A_VARIANT_EXPLAINS: [(f64, [f64; 6], f64); 4] = [
        (
            1e-5,
            [1.00001, 6.0, 6.99999, 2.99999, 4.0, 9.00001],
            5.773_502_691_645_637e-6,
        ),
        (
            1e-7,
            [1.0000001, 6.0, 6.9999999, 2.9999999, 4.0, 9.0000001],
            5.773_502_679_242_529e-8,
        ),
        (
            1e-8,
            [1.00000001, 6.0, 6.99999999, 2.99999999, 4.0, 9.00000001],
            5.773_502_897_178_33e-9,
        ),
        (0.0, [1.0, 6.0, 7.0, 3.0, 4.0, 9.0], 0.0),
    ];

    /// How far the `se` of such a variant may be from the one numpy gives
    /// the same trait: 1e-5 of it.
    ///
    /// It is not tighter because numpy's own number is not. The trait it
    /// was given is the one above, whose values are rounded to the nearest
    /// `f64` and so carry a perturbation of about 1e-16 that is not along
    /// `u`; at a `delta` of 1e-8 that perturbation is 1e-8 of the signal,
    /// so numpy's `se` and the closed form `delta / sqrt(3)` part in their
    /// eighth digit, 5.7735028971e-9 against 5.7735026918e-9. What is
    /// asserted is that popnei answers the rounded trait it was given as
    /// numpy answers it.
    const OF_NUMPYS_STANDARD_ERROR: f64 = 1e-5;

    /// A variant that explains almost all of what the null model left has
    /// a standard error, and the two backends agree on it.
    ///
    /// What a variant leaves unexplained is formed from its residuals and
    /// not by subtracting `beta * num` from the null's sum of squares,
    /// which "The linear model" of `docs/specs/gwas.md` asks for and
    /// measures: the two quantities agree to their last bits once the
    /// variant explains most of the residual, so the subtraction leaves
    /// the rounding of a cancelled sum, of either sign. At a `delta` of
    /// 1e-7 it gave exactly 0, and so an `se` of 0, which is not a
    /// standard error; at 1e-8 it gave -7.1e-15, and so an `se` of NaN,
    /// which is not one either; and which of the two a variant got
    /// differed between Accelerate and faer, so the native build and the
    /// build a browser runs disagreed about whether a variant could be
    /// tested at all.
    ///
    /// The effect is 3 at every `delta` by the construction of the trait,
    /// and the `se` is numpy 2.5.3's, from its own least squares fit of
    /// the trait on the intercept, the covariate and the dosages, with the
    /// standard error taken from the inverse of `d' d`. At a `delta` of 0
    /// the variant explains the trait exactly and there is no number to
    /// assert: what is checked is that the `se` is a finite number at the
    /// size of the rounding and not a NaN, and numpy gives 1.2e-15 there.
    #[test]
    fn a_variant_that_explains_almost_everything_has_a_standard_error() {
        let vcf = the_worked_example_vcf();
        for (delta, phenotype, se_of_numpy) in THE_TRAIT_A_VARIANT_EXPLAINS {
            let study = GwasInput {
                phenotype: &phenotype,
                ..the_worked_example_study()
            };
            let mut reader = reader_over(&vcf);
            let result = match the_study_of(&mut reader, &study) {
                Ok(result) => result,
                Err(error) => panic!("the study at a delta of {delta}: {error}"),
            };
            assert_within(
                result.beta[0],
                3.0,
                1e-12,
                &format!("the effect at a delta of {delta}"),
            );
            let se = result.se[0];
            assert!(
                se.is_finite() && se >= 0.0,
                "the standard error at a delta of {delta} is {se}, and what a variant \
                 leaves unexplained is a squared length, which is 0 or above"
            );
            let p_value = result.p_value[0];
            assert!(
                p_value.is_finite(),
                "the p-value at a delta of {delta} is {p_value}"
            );
            if delta == 0.0 {
                assert!(
                    se < 1e-13,
                    "the variant explains the trait exactly at a delta of 0, and its \
                     standard error is {se}, which is more than the rounding of the sum"
                );
                continue;
            }
            assert_within(
                se,
                se_of_numpy,
                OF_NUMPYS_STANDARD_ERROR,
                &format!("the standard error at a delta of {delta}"),
            );
        }
    }

    /// The header of a VCF of eight diploid individuals and no variant.
    const THE_HEADER_OF_EIGHT: &str = "##fileformat=VCFv4.2\n\
        ##contig=<ID=1>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti0\ti1\ti2\ti3\ti4\ti5\ti6\ti7\n";

    /// The design of the fixture of **Open 2** of `docs/specs/gwas.md`:
    /// eight individuals and a covariate that marks two subpopulations of
    /// four, beside the intercept.
    const THE_DESIGN_OF_TWO_SUBPOPULATIONS: [f64; 16] = [
        1.0, 0.0, //
        1.0, 0.0, //
        1.0, 0.0, //
        1.0, 0.0, //
        1.0, 1.0, //
        1.0, 1.0, //
        1.0, 1.0, //
        1.0, 1.0,
    ];

    /// The trait of those eight individuals.
    const THE_TRAIT_OF_TWO_SUBPOPULATIONS: [f64; 8] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];

    /// All eight of them are tested.
    const THE_INDIVIDUALS_OF_TWO_SUBPOPULATIONS: [usize; 8] = [0, 1, 2, 3, 4, 5, 6, 7];

    /// What numpy 2.5.3 answers for the ordinary variant of that fixture,
    /// the dosages `1 0 2 0 1 2 0 1`, by the same least squares fit as
    /// above: the threshold that refuses the variant beside it leaves this
    /// one alone.
    const OF_THE_ORDINARY_VARIANT: (f64, f64) = (-0.315_789_473_684_210_9, 0.633_330_903_431_213_6);

    /// A variant that the design leaves nothing of has no answer, which is
    /// the meanwhile of **Open 2** of `docs/specs/gwas.md`.
    ///
    /// The variant is twice the covariate, so the design explains all of
    /// it and what is left is rounding: `beta` is a number divided by
    /// noise, and what the study gave before this was 5.9e13 on Accelerate
    /// and -3.0e13 on faer, opposite signs with a p-value that reads as a
    /// variant that was tested and showed nothing. plink2 answers `NA` for
    /// such a variant with `ERRCODE CORR_TOO_HIGH`.
    ///
    /// The threshold is the spec's: the variant is refused when what the
    /// design leaves of its squared length is at most the tested
    /// individuals times 2.2e-16 of what there was. Measured with numpy
    /// 2.5.3 on this fixture, the collinear variant leaves 6.47e-32 of its
    /// squared length against a threshold of 1.78e-15, and the ordinary
    /// variant beside it leaves 0.432, thirteen orders of magnitude apart.
    #[test]
    fn a_variant_the_design_leaves_nothing_of_has_no_answer() {
        let mut vcf = String::from(THE_HEADER_OF_EIGHT);
        for (var, genotypes) in [
            ["0/0", "0/0", "0/0", "0/0", "1/1", "1/1", "1/1", "1/1"],
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
        let study = GwasInput {
            phenotype: &THE_TRAIT_OF_TWO_SUBPOPULATIONS,
            trait_type: TraitType::Continuous,
            design: &THE_DESIGN_OF_TWO_SUBPOPULATIONS,
            num_coefs: 2,
            kinship: None,
            test: None,
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS_OF_TWO_SUBPOPULATIONS,
            transform_to_biallelic: false,
        };
        let mut reader = reader_over(vcf.as_bytes());
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of two subpopulations: {error}"),
        };

        assert_eq!(result.num_vars, 2);
        // The variant is in the result with its frequency, as every
        // variant that has no answer is.
        assert_within(result.allele_freq[0], 0.5, 1e-12, "the frequency of v0");
        assert!(
            result.beta[0].is_nan() && result.se[0].is_nan() && result.p_value[0].is_nan(),
            "the variant the design leaves nothing of was answered with a beta of {beta}, \
             an se of {se} and a p-value of {p_value}",
            beta = result.beta[0],
            se = result.se[0],
            p_value = result.p_value[0]
        );
        let (beta, se) = OF_THE_ORDINARY_VARIANT;
        assert_within(result.beta[1], beta, 1e-12, "the effect of v1");
        assert_within(result.se[1], se, 1e-12, "the standard error of v1");
    }

    /// A study whose reader gives no variant is refused, which is what
    /// pyNei does with "There are no variants to test": the null model is
    /// fitted, and there is nothing to test against it.
    #[test]
    fn a_study_of_no_variant_is_refused() {
        let mut reader = reader_over(THE_HEADER.as_bytes());
        let study = the_worked_example_study();
        match the_study_of(&mut reader, &study) {
            Err(Error::PassGaveNoVariant { .. }) => {}
            Err(error) => panic!("the study was refused with {error}"),
            Ok(result) => panic!("a study of {} variants was run", result.num_vars),
        }
    }
}
