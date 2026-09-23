//! The association study of a dataset: which of its variants go with a
//! trait of the individuals, with the effect of each variant on the trait,
//! the uncertainty of that effect, and the chance of an effect at least
//! that large in a dataset where the variant has none, which is its
//! p-value. `docs/specs/gwas.md` says what a study gives and how each of
//! its models is verified.
//!
//! What is here is the two functions that turn the statistic of a test
//! into a p-value, which every model of the module ends in, and the design
//! every model is fitted on. Each of the two is a survival function, the
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
//! models that give the statistics to the two distributions, and the pass
//! over the variants that feeds them, are being written.

use std::fmt;

use crate::error::{Error, Result};

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
/// is. Nothing calls this with such a `df` today; what will keep it at 1 or
/// above is the refusal of a design with no more rows than columns plus
/// one, which is not written yet, so this is a guard and not the repair of
/// a live wrong number.
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
/// intercept first and then one for each covariate.
///
/// [`Design::of_the_study`] is the only way to have one, so a model that
/// takes a `Design` is fitted on individuals the source has, each of them
/// once and in the source's order, on a phenotype that holds a number for
/// each of them and fits the trait, and on columns that are independent.
#[derive(Debug, Clone, Copy)]
pub struct Design<'a> {
    values: &'a [f64],
    num_individuals: usize,
    num_coefs: usize,
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
    /// [`Error::GwasCovariatesCollinear`] when the columns of the design
    /// are not independent, and [`Error::GwasLinalg`] when the rank that
    /// finds that out could not be taken.
    /// [`Error::GwasInputOfAnotherSize`] when the phenotype or the design
    /// does not hold one value or one row for each tested individual, or
    /// the design has no column.
    pub fn of_the_study(
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
        })
    }

    /// The design itself, `num_individuals` x `num_coefs`, row after row.
    #[must_use]
    pub fn values(&self) -> &'a [f64] {
        self.values
    }

    /// How many individuals are tested, which is the rows of the design.
    #[must_use]
    pub fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// How many columns the design has: the intercept and one for each
    /// covariate.
    #[must_use]
    pub fn num_coefs(&self) -> usize {
        self.num_coefs
    }
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
