//! What a user asks a study for, what a study is refused for before any
//! model is fitted or any variant is read, and which model it fits.
//!
//! [`GwasInput`] is what a caller brings: the individuals to test, their
//! trait, the design the model is fitted on, the kinship and the test.
//! [`TraitType`] and [`TestType`] are what a user names of those, and
//! [`GwasModel`] is which of the four models the study fits, which
//! [`the_model_and_the_test`] chooses from the trait and the kinship
//! together. [`Design`] is the design once it has been checked, and
//! [`Design::of_the_study`] is the only way to have one.

use std::fmt;

use crate::error::{Error, Result};

use super::dosages::MultiallelicVariants;

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
    /// The kinship does not hold one row and one column for each tested
    /// individual.
    Kinship {
        /// How many values the kinship holds.
        num_values: usize,
        /// How many individuals are tested.
        num_individuals: usize,
    },
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
            Self::Kinship {
                num_values,
                num_individuals,
            } => write!(
                formatter,
                "the kinship holds {num_values} values and {num_individuals} individuals are tested, and it is {num_individuals} x {num_individuals}, row after row, cut to them and in their order"
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
/// in the same order, and
/// [`GwasDosages::read_the_block`](super::dosages::GwasDosages::read_the_block)
/// takes a `Design` so that no pass can read a block over positions that
/// nothing checked.
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
    pub(super) multiallelic: MultiallelicVariants,
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
    /// such a trait has the same one, which leaves one of the two groups
    /// it compares empty. [`Error::GwasContinuousPhenotypeOfOneValue`] when
    /// every tested individual of a continuous trait has the same one,
    /// which is a measurement that does not differ between them.
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

/// The kinship a study was given: one row and one column for each tested
/// individual, and a finite number in every cell.
///
/// Neither is a thing a fit would notice, which is why they are checked
/// here and before any model is fitted, as "The Rust interface" of
/// `docs/specs/gwas.md` asks. A matrix of another length is read as
/// another shape and gives numbers; a value that is not finite spreads
/// through the eigendecomposition into every eigenvalue and every
/// eigenvector, and the study comes back with a NaN for every variant and
/// the linear algebra crate's refusal of a matrix at whichever routine met
/// it first, which names no cell.
///
/// `kinship` is `num_individuals` x `num_individuals`, row after row, so
/// the row of a value is one tested individual and the rest of the
/// division is the other.
pub(super) fn refuse_a_kinship_that_is_not_of_the_individuals(
    kinship: &[f64],
    num_individuals: usize,
) -> Result<()> {
    if num_individuals
        .checked_mul(num_individuals)
        .is_none_or(|values| kinship.len() != values)
    {
        return Err(Error::GwasInputOfAnotherSize {
            problem: GwasInputShape::Kinship {
                num_values: kinship.len(),
                num_individuals,
            },
        });
    }
    for (individual, row) in kinship.chunks(num_individuals.max(1)).enumerate() {
        for (other, value) in row.iter().copied().enumerate() {
            if !value.is_finite() {
                return Err(Error::GwasKinshipValueNotFinite {
                    individual,
                    other,
                    value,
                });
            }
        }
    }
    refuse_a_kinship_that_is_not_symmetric(kinship, num_individuals)
}

/// How far a kinship may be from its own transpose before a study refuses
/// it, as a share of its largest absolute entry: 1e-9.
///
/// It is the tolerance the `Kinship` of both packages refuses one at, so a
/// matrix a user could build is refused in the same place whichever layer
/// they came through. The two reference panels' matrices are symmetric to
/// the bit, and what this leaves room for is a matrix a program wrote one
/// triangle of and rounded, which is why it is not exact equality.
const LARGEST_ASYMMETRY_OF_A_KINSHIP: f64 = 1e-9;

/// The kinship a study was given, checked against its own transpose.
///
/// The eigendecomposition reads the lower triangle alone, so a matrix that
/// holds two different numbers for one pair is read as that half mirrored,
/// and a study over it answers with numbers that are of another matrix than
/// the user's. The `Kinship` of both packages checks this when it is built;
/// what reaches here is a frame written into afterwards, which neither
/// package sees, and a caller of the core crate.
///
/// The rows are walked once and the upper triangle is compared with the
/// lower, so nothing of the size of the matrix is allocated: the kinship of
/// 10000 individuals is 800 MB and a check written as a difference with the
/// transpose asks for that much again.
///
/// # Errors
///
/// [`Error::GwasKinshipNotSymmetric`] at the first pair whose two cells are
/// further apart than [`LARGEST_ASYMMETRY_OF_A_KINSHIP`] of the largest
/// absolute entry of the matrix.
fn refuse_a_kinship_that_is_not_symmetric(kinship: &[f64], num_individuals: usize) -> Result<()> {
    let largest = kinship
        .iter()
        .fold(0.0_f64, |largest, value| largest.max(value.abs()));
    let allowed = largest * LARGEST_ASYMMETRY_OF_A_KINSHIP;
    for (individual, row) in kinship.chunks(num_individuals.max(1)).enumerate() {
        for (other, value) in row.iter().copied().enumerate().skip(individual) {
            let and_back = other
                .checked_mul(num_individuals)
                .and_then(|at| at.checked_add(individual))
                .and_then(|at| kinship.get(at).copied());
            let Some(and_back) = and_back else {
                continue;
            };
            if (value - and_back).abs() > allowed {
                return Err(Error::GwasKinshipNotSymmetric {
                    individual,
                    other,
                    value,
                    and_back,
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
/// them, not the same number in all of them, and for a binomial trait 0 or
/// 1 with somebody in each of the two groups.
///
/// A trait of one value is refused for both kinds, which "Which
/// individuals are tested, and the design" of `docs/specs/gwas.md` asks
/// for and which pyNei does for a binomial trait alone. The two get the
/// two errors, each with its own reason: a group of the binomial
/// comparison that nobody is in, and a measurement that does not differ
/// between the individuals.
#[expect(
    clippy::float_cmp,
    reason = "a binomial phenotype is the 0.0 and the 1.0 themselves and not a \
              measurement near either, so what is wanted here is the exact \
              comparison and not one within a tolerance; a value of -0.0 is 0.0 \
              by it, which is the answer for an individual without the condition. \
              A continuous trait is the same number in every individual when \
              every one of its finite values has the bits of the first, two \
              measurements that differ at all differing in a bit, so that \
              comparison is the exact one too"
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
        // Every value is finite by here, so a trait that is the same in
        // every individual is one whose values all equal the first.
        TraitType::Continuous => match phenotype.split_first() {
            Some((first, rest)) if rest.iter().all(|value| *value == *first) => {
                Err(Error::GwasContinuousPhenotypeOfOneValue { value: *first })
            }
            Some(_) | None => Ok(()),
        },
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
    /// empty. Both 0 and 1 are the same refusal, and a phenotype with
    /// somebody in each group is kept. A continuous trait of one value is
    /// refused too, by the test below, and with the other error of the two:
    /// the reason it gives is the one of a measurement that does not
    /// differ, and not an empty group of a condition nobody was asked
    /// about.
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
        }
        let phenotype = [0.0, 1.0, 1.0, 0.0];
        let mut study = a_study(&phenotype, &DESIGN_OF_FOUR, 2, &four);
        study.trait_type = TraitType::Binomial;
        assert!(
            Design::of_the_study(&study, 4).is_ok(),
            "two individuals have the condition and two have not"
        );
    }

    /// A continuous trait where every tested individual has the same value
    /// is refused, naming the value. There is nothing for such a trait to
    /// be associated with, and what a fit gives instead is not an answer:
    /// with no kinship the residual sum of squares is 0 and every `se` is
    /// 0, and with a kinship the genetic variance is fitted at 0 and the
    /// inverse it feeds returns infinities. pyNei refuses a trait of one
    /// value only for a binomial trait, and "Which individuals are tested,
    /// and the design" of `docs/specs/gwas.md` records that difference.
    ///
    /// The fixture is four individuals of the phenotype 2.5, which is
    /// neither 0 nor 1, so the refusal is not the binomial one reached by
    /// another route; and then the same four with one of them at 2.75,
    /// which is kept. The two are either side of the refusal, so a check
    /// that asked for more than one value to differ fails here.
    #[test]
    fn a_continuous_phenotype_of_one_value_is_refused() {
        let four = [0, 1, 2, 3];
        let phenotype = [2.5; 4];
        let study = a_study(&phenotype, &DESIGN_OF_FOUR, 2, &four);
        match Design::of_the_study(&study, 4) {
            Err(Error::GwasContinuousPhenotypeOfOneValue { value }) => {
                assert_eq!(value.to_bits(), 2.5_f64.to_bits());
            }
            other => panic!(
                "every tested individual has the phenotype 2.5 of a continuous trait, \
                 and that gave {other:?}"
            ),
        }
        let one_differs = [2.5, 2.5, 2.75, 2.5];
        let study = a_study(&one_differs, &DESIGN_OF_FOUR, 2, &four);
        assert!(
            Design::of_the_study(&study, 4).is_ok(),
            "one individual of the four has another measurement, which is a trait \
             that differs"
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

/// Which model a study fits and which test it makes of every variant, as
/// "What it gives" and "Which individuals are tested, and the design" of
/// `docs/specs/gwas.md` state them.
#[cfg(test)]
mod choice {
    use super::{GwasModel, TestType, TraitType, the_model_and_the_test};
    use crate::error::Error;
    use crate::gwas::fixtures::{TESTED_OF_EIGHT, a_study, the_phenotype_and_the_design_of};

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
