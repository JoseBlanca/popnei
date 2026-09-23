//! The principal component analysis of a table of numbers.
//!
//! A principal component analysis places the rows of a table, the
//! individuals, on a few axes that hold as much of the variation between
//! them as that many axes can. Each column of the table, a trait, is
//! centered and divided by its standard deviation, which gives Z, and the
//! components are the directions in the space of the traits along which
//! the rows of Z vary most, the first having the largest variance that any
//! direction has, the second the largest among the directions at a right
//! angle to the first, and so on. [`Pca`] holds where each individual falls
//! along each component, how much of the variance each component holds and
//! the weight of each trait in each component.
//!
//! [`pca`] takes the table the user brings and [`pca_of_variants`] the
//! variants a reader gives, each of which becomes one number per
//! individual, its dosage. The variants are read in blocks and nothing of
//! the size variants x individuals is ever held: the individuals x
//! individuals matrix of the products, which is what the components come
//! from, is added up block by block.
//!
//! The products of matrices and the eigendecomposition are those of the
//! crate `popnei-linalg`, which runs them on the BLAS and LAPACK of the
//! system natively, and on faer in WebAssembly and natively when the cargo
//! feature `blas` is off.

use std::cmp::Ordering;
use std::fmt;
use std::num::NonZeroUsize;

use popnei_linalg::{Eigen, TheSecondOperand, add_self_product_lower, eigh_lower, product};

use crate::block::{Block, BlockReader, Reblock};
use crate::error::{Error, Result};
use crate::variant::{AlleleCounts, MISSING_ALLELE, Needs, count_alleles, the_major_allele};

/// Whether the table is centered, which `do_pca` of pyNei does by default
/// and so does popnei.
pub const DEFAULT_CENTER_DATA: bool = true;

/// Whether the table is standardized, which `do_pca` of pyNei does by
/// default and so does popnei. A table that is standardized is centered as
/// well, so turning the centering off turns this off too.
pub const DEFAULT_STANDARDIZE_DATA: bool = true;

/// How close two projections are when the rule that fixes the sign of a
/// component takes them for one absolute value: 64 times the difference
/// between 1 and the next number an `f64` holds, which is 1.4e-14 of the
/// larger of the two.
///
/// Two individuals whose projections are the same number with opposite
/// signs in exact arithmetic do not come out equal bit for bit, and which
/// of the two is the larger is not the same in the two libraries popnei
/// takes its eigendecomposition from: without this tolerance the sign of
/// the whole component would be decided by the last bit and would differ
/// between Python natively, Python under pyodide and TypeScript. "The
/// sign" of `docs/specs/pca.md` has the case it was measured on and why
/// the tolerance is 64 units in the last place.
pub const TIED_PROJECTIONS: f64 = 64.0 * f64::EPSILON;

/// The two steps on the columns of a table before its components are
/// taken.
///
/// Centering takes the mean of each trait from it, which is what makes the
/// components the directions of the variation and not directions that
/// point at the mean of the data. Standardizing then divides each trait by
/// its standard deviation, which puts traits measured in different units
/// on one scale; without it the traits with the largest numbers dominate.
/// Standardizing without centering is an error, since the standard
/// deviation it divides by is the one the trait has once it is centered.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PcaOptions {
    /// Whether the mean of each trait is taken from it.
    pub center: bool,
    /// Whether each trait is divided by its standard deviation, the one
    /// with the number of rows in it and not the number of rows less one,
    /// which is pyNei's.
    pub standardize: bool,
}

/// What the mean or the standard deviation of a trait came out as when it
/// is not a number the analysis can use.
///
/// Each of the three is a trait whose values are too large or too small
/// for the arithmetic of an `f64`, and the user scales that trait or takes
/// it out of the table. They are found before anything is computed from
/// the mean or the deviation, so no analysis is done on the numbers they
/// would give.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraitScale {
    /// The values of the trait sum above the largest `f64`, 1.8e308, so
    /// its mean is an infinity and every centered value of it would be a
    /// NaN.
    MeanNotFinite,
    /// The squares of the deviations of the trait sum above the largest
    /// `f64`, which values of 1e154 give, so its standard deviation is an
    /// infinity and the standardized trait would be a column of zeros,
    /// which is what a trait with no variance gives.
    DeviationNotFinite,
    /// The squares of the deviations of the trait all fall below the
    /// smallest `f64` above 0, 5e-324, which values of 1e-200 give, so its
    /// standard deviation is 0 although its values are not all equal, and
    /// dividing by it would give infinities.
    DeviationOfZero,
}

impl fmt::Display for TraitScale {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let said = match *self {
            Self::MeanNotFinite => {
                "its values sum above the largest f64, so its mean is not finite"
            }
            Self::DeviationNotFinite => {
                "the squares of its deviations sum above the largest f64, so its standard deviation is not finite"
            }
            Self::DeviationOfZero => {
                "the squares of its deviations are all below the smallest f64 above 0, so its standard deviation is 0 although its values are not all equal"
            }
        };
        formatter.write_str(said)
    }
}

/// What a principal component analysis gives.
///
/// Only the components that have variance are here, and `num_comps` is how
/// many of them there are: centering takes one dimension out of the data,
/// so a table of 8 rows and 30 traits has 7 components and not 8. In each
/// component the projection of the largest absolute value is positive,
/// which is the rule that makes the result the same whichever library did
/// the eigendecomposition, and the weights of that component have the sign
/// that rule gave it.
#[derive(Debug, Clone)]
pub struct Pca {
    /// How many columns the data had: the variants the first pass gave,
    /// used or not, which is `num_vars` of the pass stats, or the traits.
    pub num_cols: usize,
    /// How many rows the data had, which is the individuals.
    pub num_rows: usize,
    /// How many components have variance, which is how many are given.
    pub num_comps: usize,
    /// num_rows x num_comps, row after row.
    pub projections: Vec<f64>,
    /// One per component, over the variance of every component of the data.
    pub explained_variance_percent: Vec<f64>,
    /// The positions of the columns that were used: the variants with
    /// variance among those the reader gave, or every trait of a table.
    pub used_cols: Vec<usize>,
    /// How many components the weights are given for, which for a table is
    /// `num_comps`.
    pub num_prin_comps: usize,
    /// The weight of each column that was used in each of the first
    /// `num_prin_comps` components: `num_prin_comps` x `used_cols.len()`,
    /// row after row. Its columns are the columns that were used and not
    /// the `num_cols` the data had, so a variant that was left out is in
    /// neither.
    pub princomps: Vec<f64>,
}

/// The principal components of the table `data`, which is `num_rows` rows
/// of `num_cols` values each, row after row, the rows being the
/// individuals and the columns the traits.
///
/// Every trait is a column of the result's weights, the one with no
/// variance included, which gets a weight of 0 when the table is not
/// standardized. No value may be missing: a table comes whole.
///
/// # Errors
///
/// [`Error::PcaValueNotFinite`] when a value of the table is an infinity
/// or a NaN. [`Error::PcaStandardizeWithoutCentering`] when the options
/// ask for the second and not the first. [`Error::PcaTableTooSmall`] when
/// the table has fewer than 2 rows or no traits.
/// [`Error::PcaTraitsWithNoVariance`] when the table is standardized and a
/// trait has no variance, and [`Error::PcaNoTraitWithVariance`] when no
/// trait of it has any. [`Error::PcaTraitOutOfRange`] when the mean or the
/// standard deviation of a trait is not a number the analysis can use.
/// [`Error::PcaTableOfAnotherSize`] when `data` does not hold exactly
/// `num_rows` times `num_cols` values. [`Error::PcaLinalg`] when the
/// product or the eigendecomposition could not be done.
pub fn pca(data: &[f64], num_rows: usize, num_cols: usize, options: &PcaOptions) -> Result<Pca> {
    if options.standardize && !options.center {
        return Err(Error::PcaStandardizeWithoutCentering);
    }
    if num_rows < 2 || num_cols == 0 {
        return Err(Error::PcaTableTooSmall { num_rows, num_cols });
    }
    // The buffer holds the table and nothing more: a longer one would be
    // analysed on its first values, which is not the table its caller
    // meant.
    let table = match num_rows.checked_mul(num_cols) {
        Some(num_values) if num_values == data.len() => data,
        _ => {
            return Err(Error::PcaTableOfAnotherSize {
                num_values: data.len(),
                num_rows,
                num_cols,
            });
        }
    };
    refuse_a_value_that_is_not_finite(table, num_cols)?;
    let (means, deviations) =
        the_center_and_the_scale_of_each_trait(table, num_rows, num_cols, options)?;
    // The matrix that is decomposed is the product of the smaller of the
    // two sides of the table with itself, Z Z' when there are fewer rows
    // than traits and Z' Z otherwise, and `add_self_product_lower` gives
    // the product of a matrix with itself over its columns. So the copy
    // that is centered and standardized is written with the traits as its
    // rows in the first case and as its columns in the second.
    let layout = if num_rows <= num_cols {
        Layout::TraitsAsRows
    } else {
        Layout::TraitsAsColumns
    };
    let standardized =
        the_standardized_table(table, num_rows, num_cols, &means, &deviations, layout);
    let (num_summed, side) = match layout {
        Layout::TraitsAsRows => (num_cols, num_rows),
        Layout::TraitsAsColumns => (num_rows, num_cols),
    };
    let mut gram = vec![0.0; num_values_of(side, side)];
    add_self_product_lower(&standardized, num_summed, side, &mut gram).map_err(|source| {
        Error::PcaLinalg {
            operation: "product of the table with itself",
            source,
        }
    })?;
    let eigen = eigh_lower(gram, side).map_err(|source| Error::PcaLinalg {
        operation: "eigendecomposition",
        source,
    })?;
    let num_comps = the_components_with_variance(&eigen.values, num_rows, num_cols);
    if num_comps == 0 {
        return Err(Error::PcaNoTraitWithVariance);
    }
    let explained_variance_percent = the_percentages_of(&eigen.values, num_comps);
    let (mut projections, mut princomps) = match layout {
        Layout::TraitsAsRows => the_components_of_the_product_of_the_rows(
            &standardized,
            &eigen,
            num_rows,
            num_cols,
            num_comps,
        )?,
        Layout::TraitsAsColumns => the_components_of_the_product_of_the_traits(
            &standardized,
            &eigen,
            num_rows,
            num_cols,
            num_comps,
        )?,
    };
    fix_the_signs(&mut projections, &mut princomps, num_comps, num_cols);
    Ok(Pca {
        num_cols,
        num_rows,
        num_comps,
        projections,
        explained_variance_percent,
        used_cols: (0..num_cols).collect(),
        num_prin_comps: num_comps,
        princomps,
    })
}

/// Whether every allele that is not the major one counts the same in a
/// variant that has more than two, which neither `do_pca_from_variants` of
/// pyNei nor popnei does unless it is asked for: without it such a variant
/// is an error.
pub const DEFAULT_TRANSFORM_TO_BIALLELIC: bool = false;

/// How many components the weights of the variants are given for. pyNei
/// has no such argument and gives the weight of every variant in every
/// component, which is 0.8 GB for 100000 variants of 1000 individuals and
/// 80 GB for a million variants of 10000; the owner decided on 21
/// September 2026 that popnei gives the first ten.
pub const DEFAULT_NUM_PRIN_COMPS: usize = 10;

/// The largest ploidy the principal components of the variants are taken
/// at, which is one less than the largest the VCF reader takes.
///
/// The first pass writes the genotype of each individual as one byte, its
/// dosage or the missing genotype, which is the loop the compiler
/// vectorizes. A ploidy of 255 has 256 dosages, and those with the missing
/// genotype are one value more than a byte holds.
pub const MAX_PLOIDY_OF_THE_VARIANTS: usize = 254;

/// The most individuals the principal components of the variants are taken
/// on: the largest number whose square is at most the 2147483647 values
/// that the routines of BLAS and LAPACK count a matrix in, which
/// `crates/popnei-linalg` gives as
/// [`THE_MOST_VALUES_OF_A_MATRIX`](popnei_linalg::THE_MOST_VALUES_OF_A_MATRIX).
///
/// The individuals x individuals matrix of that many holds 2147395600
/// values, 17 GB, which no browser tab gives and few machines do.
pub const MAX_INDIVIDUALS_OF_THE_VARIANTS: usize =
    popnei_linalg::THE_MOST_VALUES_OF_A_MATRIX.isqrt();

/// Which size of a dataset is beyond what the principal components of its
/// variants are taken on.
///
/// Each of the four is a number the analysis counts in, and a dataset
/// above it would be read into a number that wrapped. None of them is a
/// dataset of this world: the largest ploidy of an organism is a dozen,
/// and the objectives of popnei reach 10000 individuals and a million
/// variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantsTooLarge {
    /// The ploidy of the dataset, which is above
    /// [`MAX_PLOIDY_OF_THE_VARIANTS`].
    Ploidy(usize),
    /// The individuals of the dataset, which are more than
    /// [`MAX_INDIVIDUALS_OF_THE_VARIANTS`].
    Individuals(usize),
    /// The reader gave more variants than a `usize` counts, which is
    /// 4294967295 in WebAssembly, where a `usize` is 32 bits.
    Variants,
    /// The weights asked for are more values than a `usize` counts: that
    /// many components of that many variants that were used. Fewer
    /// variants reach it than the count of the variants does, since it
    /// grows with the components the weights are asked for.
    Weights {
        /// How many components the weights are given for.
        num_prin_comps: usize,
        /// How many variants were used.
        num_used: usize,
    },
}

impl fmt::Display for VariantsTooLarge {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::Ploidy(ploidy) => write!(
                formatter,
                "its genotypes hold {ploidy} alleles each, and the first pass writes each genotype as one byte, its dosage or the missing genotype, which takes a ploidy of {MAX_PLOIDY_OF_THE_VARIANTS} at most"
            ),
            Self::Individuals(num_individuals) => write!(
                formatter,
                "it has {num_individuals} individuals, and the individuals x individuals matrix of the analysis would hold more values than the 2147483647 the linear algebra counts in, which is {MAX_INDIVIDUALS_OF_THE_VARIANTS} individuals"
            ),
            Self::Variants => write!(
                formatter,
                "the variants given are more than {largest}, which is what this machine counts the columns of the analysis in",
                largest = usize::MAX
            ),
            Self::Weights {
                num_prin_comps,
                num_used,
            } => write!(
                formatter,
                "the weights of {num_prin_comps} components of the {num_used} variants that were used are more values than the {largest} this machine counts them in; ask for fewer components",
                largest = usize::MAX
            ),
        }
    }
}

/// How the variants of the second pass differ from those of the first,
/// which is what a source that changed between the two passes gives.
///
/// The weights of the second pass belong to the variants of the first, the
/// ones the eigenvectors were worked out from, so a pass over other
/// variants has nothing to give.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantsOfTheSecondPass {
    /// The second pass is over another dataset: that many individuals of
    /// that ploidy, where the first pass had the ones the analysis was
    /// done on.
    Dataset {
        /// The individuals of the second pass.
        num_individuals: usize,
        /// The ploidy of the second pass.
        ploidy: usize,
    },
    /// The second pass gave another number of variants.
    Count {
        /// How many variants the second pass gave.
        found: usize,
        /// How many the first pass gave.
        expected: usize,
    },
    /// The variant at that position, among the variants each pass gave,
    /// has variance in the second pass and had none in the first.
    UsedNow(usize),
    /// The variant at that position had variance in the first pass and has
    /// none in the second.
    NotUsedNow(usize),
}

impl fmt::Display for VariantsOfTheSecondPass {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::Dataset {
                num_individuals,
                ploidy,
            } => write!(
                formatter,
                "it reads {num_individuals} individuals of the ploidy {ploidy}, which are not the individuals of the first pass"
            ),
            Self::Count { found, expected } => write!(
                formatter,
                "it gave {found} variants and the first pass gave {expected}"
            ),
            Self::UsedNow(position) => write!(
                formatter,
                "the variant at the position {position} has variance now and had none in the first pass"
            ),
            Self::NotUsedNow(position) => write!(
                formatter,
                "the variant at the position {position} had variance in the first pass and has none now"
            ),
        }
    }
}

/// What the principal components of the variants are taken with.
///
/// The dosage of a genotype, how many of its alleles are not the major
/// one, has a meaning for a variant of two alleles; `transform_to_biallelic`
/// is what says that a variant of more is read with every allele that is
/// not the major one counting the same, and without it such a variant is
/// an error. `num_prin_comps` is how many components the weights of the
/// variants are given for, and with 0 there are no weights and no second
/// pass over the variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VariantPcaOptions {
    /// Whether every allele that is not the major one counts the same.
    pub transform_to_biallelic: bool,
    /// How many components the weights are given for.
    pub num_prin_comps: usize,
}

/// The principal components of the variants of a dataset, which the
/// readers give in blocks.
///
/// Each variant becomes one number per individual, its dosage: how many
/// alleles of the genotype are not the major allele of the variant, which
/// is the most frequent among its called alleles and the lowest numbered
/// of two that are equally frequent. A genotype with an allele missing
/// takes the mean of the dosages of its variant, so that after centering
/// it pulls its individual nowhere, and the standard deviation the variant
/// is divided by has all the individuals in it and not the called ones. A
/// variant whose called genotypes all have one dosage, and one with no
/// called genotype, have no variance and are left out; the ones that were
/// used are [`Pca::used_cols`], their positions among the variants the
/// reader gave.
///
/// The weights of the first `num_prin_comps` components, the weight of
/// each variant that was used in each of them, come from a second pass
/// over the same variants: a weight needs the eigenvectors, which are
/// known when the first pass ends. More components than there are gives
/// those there are, and [`Pca::num_prin_comps`] says how many that was.
///
/// The two readers are over the same variants, and opening one reads no
/// variant. `second_pass` is `None` when `num_prin_comps` is 0 and is not
/// read then. This pass borrows the readers and does not take them, so
/// that whoever built the chain of filters reads their counts when it
/// returns. The function asks each reader for the genotypes alone and puts
/// [`Reblock`] before it, since a filter leaves blocks of uneven size and
/// the product of a block is matrix work.
///
/// # Errors
///
/// [`Error::PcaNoVariants`] when the reader gives no variant and
/// [`Error::PcaNoVariantWithVariance`] when no variant of it has variance.
/// [`Error::PcaVariantWithMoreThanTwoAlleles`] when a variant has more
/// than two different alleles among its called genotypes and
/// `transform_to_biallelic` is false. [`Error::PcaSecondPassMissing`] when
/// `num_prin_comps` is above 0 and no reader was given for the second
/// pass, and [`Error::PcaSecondPassDiffers`] when that reader gives other
/// variants than the first, which is what a source that changed between
/// the two passes gives.
/// [`Error::PcaVariantsTooLarge`] when a size of the dataset is beyond
/// what the analysis counts in, which [`VariantsTooLarge`] lists.
/// [`Error::FieldsNotInTheBlock`] when a block holds no genotypes, and
/// whatever the reader fails with. [`Error::PcaLinalg`] when the product
/// of a block or the eigendecomposition could not be done.
pub fn pca_of_variants<R1: BlockReader, R2: BlockReader>(
    first_pass: &mut R1,
    second_pass: Option<&mut R2>,
    options: &VariantPcaOptions,
) -> Result<Pca> {
    let ploidy = first_pass.ploidy();
    if ploidy > MAX_PLOIDY_OF_THE_VARIANTS {
        return Err(Error::PcaVariantsTooLarge {
            problem: VariantsTooLarge::Ploidy(ploidy),
        });
    }
    let num_individuals = first_pass.individuals().len();
    if num_individuals == 0 {
        // The standardizing reads the rows of a block in chunks of one
        // genotype for each individual, which would be chunks of no
        // allele, and there is nobody to place on the components. No
        // reader of popnei gives such a source, and this function is
        // public.
        return Err(Error::PcaNoIndividual);
    }
    if num_individuals > MAX_INDIVIDUALS_OF_THE_VARIANTS {
        return Err(Error::PcaVariantsTooLarge {
            problem: VariantsTooLarge::Individuals(num_individuals),
        });
    }
    // The weight of a variant needs the eigenvectors, which are known when
    // the first pass ends, so the weights come from a second pass over the
    // same variants. Whoever asked for them and gave no reader for that
    // pass hears it before the variants go by and not after.
    if options.num_prin_comps > 0 && second_pass.is_none() {
        return Err(Error::PcaSecondPassMissing {
            num_prin_comps: options.num_prin_comps,
        });
    }
    first_pass.set_needs(Needs::GTS);
    // A filter leaves blocks of uneven size and the product of a block is
    // matrix work, so the blocks are put back to one size first.
    let mut blocks = Reblock::new(first_pass, None)?;
    let FirstPass {
        gram,
        used_cols,
        num_cols,
    } = the_first_pass(&mut blocks, options, num_individuals, ploidy)?;
    if num_cols == 0 {
        return Err(Error::PcaNoVariants);
    }
    if used_cols.is_empty() {
        return Err(Error::PcaNoVariantWithVariance);
    }
    let eigen = eigh_lower(gram, num_individuals).map_err(|source| Error::PcaLinalg {
        operation: "eigendecomposition",
        source,
    })?;
    // The columns of Z are the variants that were used, so they are the
    // side the threshold of a component with no variance is taken with.
    let num_comps = the_components_with_variance(&eigen.values, num_individuals, used_cols.len());
    if num_comps == 0 {
        return Err(Error::PcaNoVariantWithVariance);
    }
    let explained_variance_percent = the_percentages_of(&eigen.values, num_comps);
    let mut projections = the_projections_of(&eigen, num_individuals, num_comps);
    // The weights of more components than there are are the weights of
    // those there are.
    let num_prin_comps = options.num_prin_comps.min(num_comps);
    let mut princomps = match (num_prin_comps, second_pass) {
        // With 0 there is no second pass and no weight, and the variants
        // that were used are the columns of `princomps` all the same.
        (0, _) => Vec::new(),
        // A `num_prin_comps` above 0 with no second reader was refused
        // before the first variant was read.
        (num_prin_comps, None) => return Err(Error::PcaSecondPassMissing { num_prin_comps }),
        (num_prin_comps, Some(second_pass)) => {
            if second_pass.individuals().len() != num_individuals || second_pass.ploidy() != ploidy
            {
                return Err(Error::PcaSecondPassDiffers {
                    problem: VariantsOfTheSecondPass::Dataset {
                        num_individuals: second_pass.individuals().len(),
                        ploidy: second_pass.ploidy(),
                    },
                });
            }
            second_pass.set_needs(Needs::GTS);
            let mut blocks = Reblock::new(second_pass, None)?;
            let after = AfterTheFirstPass {
                scaled_vectors: &the_scaled_vectors_of(&eigen, num_individuals, num_prin_comps),
                num_prin_comps,
                used_cols: &used_cols,
                num_cols,
            };
            the_weights_of_a_second_pass(&mut blocks, options, num_individuals, ploidy, &after)?
        }
    };
    fix_the_signs(&mut projections, &mut princomps, num_comps, used_cols.len());
    Ok(Pca {
        num_cols,
        num_rows: num_individuals,
        num_comps,
        projections,
        explained_variance_percent,
        used_cols,
        num_prin_comps,
        princomps,
    })
}

/// What the first pass leaves the second: what each block of it is
/// multiplied by, and which variants the first pass used.
struct AfterTheFirstPass<'a> {
    /// The individuals x `num_prin_comps` matrix of the eigenvectors of
    /// the first components divided by sqrt of their eigenvalue, row after
    /// row, which a block of standardized rows times gives the weight of
    /// each of its variants in each of those components.
    scaled_vectors: &'a [f64],
    /// How many components the weights are given for, which is at most how
    /// many components have variance.
    num_prin_comps: usize,
    /// The positions of the variants the first pass used, in order.
    used_cols: &'a [usize],
    /// How many variants the first pass gave, used or not.
    num_cols: usize,
}

/// The weights of the first `num_prin_comps` components, from a second
/// pass over the variants: `num_prin_comps` x the variants the first pass
/// used, row after row.
///
/// Each block is standardized as the first pass standardized it and
/// multiplied by the eigenvectors divided by sqrt(λ), which gives the
/// weight of each variant of the block in each component. The product
/// gives them variant after variant and the result holds them component
/// after component, so the weights of each block are written into the
/// columns that belong to its variants and no copy of the whole matrix is
/// made.
///
/// # Errors
///
/// [`Error::PcaSecondPassDiffers`] when the variants are not those of the
/// first pass, what the standardizing of a block refuses, and
/// [`Error::PcaLinalg`] when a product could not be done.
fn the_weights_of_a_second_pass<R: BlockReader>(
    reader: &mut R,
    options: &VariantPcaOptions,
    num_individuals: usize,
    ploidy: usize,
    after: &AfterTheFirstPass<'_>,
) -> Result<Vec<f64>> {
    let num_used = after.used_cols.len();
    // The weights are the one matrix here whose two sides are not both at
    // most the individuals, so this is the product that is taken with a
    // check: 80 MB for 10 components of a million variants, and more
    // values than a `usize` counts where a `usize` is 32 bits.
    let num_values =
        after
            .num_prin_comps
            .checked_mul(num_used)
            .ok_or(Error::PcaVariantsTooLarge {
                problem: VariantsTooLarge::Weights {
                    num_prin_comps: after.num_prin_comps,
                    num_used,
                },
            })?;
    let mut princomps = vec![0.0; num_values];
    // The buffer a block is standardized into and the one the product
    // writes the weights of its variants into, both kept from one block to
    // the next.
    let mut standardized: Vec<f64> = Vec::new();
    let mut of_the_block: Vec<f64> = Vec::new();
    let mut num_cols = 0_usize;
    // How many variants that were used the passes before this block gave,
    // which is the column of `princomps` its first one goes into.
    let mut used_before = 0_usize;
    while let Some(block) = reader.next_block()? {
        let used = the_standardized_block(
            &block,
            num_individuals,
            ploidy,
            options,
            num_cols,
            &mut standardized,
        )?;
        let kept = the_variants_of_the_first_pass(&used, num_cols, used_before, after)?;
        // The variants of a block that were used are at most all of the
        // ones the first pass used, so this product is at most the values
        // of the weights, which were counted above.
        let needed = num_values_of(kept, after.num_prin_comps);
        if of_the_block.len() < needed {
            of_the_block.resize(needed, 0.0);
        }
        // The weight of each variant of the block in each component: the
        // standardized rows of the block, which are the columns of Z of
        // those variants, times the eigenvectors divided by sqrt(λ). A
        // block whose rows all had no variance has no row here and writes
        // nothing.
        product(
            &standardized,
            kept,
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: after.scaled_vectors,
                cols: after.num_prin_comps,
            },
            &mut of_the_block,
        )
        .map_err(|source| Error::PcaLinalg {
            operation: "product that gives the weights of a block of variants",
            source,
        })?;
        // The product gives the weights variant after variant and the
        // result holds them component after component, so each variant of
        // the block writes its weights into its own column of every
        // component. No copy of the whole matrix is made.
        for (variant, of_the_variant) in of_the_block
            .chunks_exact(after.num_prin_comps)
            .take(kept)
            .enumerate()
        {
            let column = used_before
                .checked_add(variant)
                .ok_or_else(the_variants_are_too_many)?;
            for (of_the_component, weight) in
                princomps.chunks_exact_mut(num_used).zip(of_the_variant)
            {
                // The column is below the variants the first pass used,
                // which this pass has been counting against them, so a
                // column that is not there is a defect of popnei and not
                // a weight to drop in silence.
                let Some(target) = of_the_component.get_mut(column) else {
                    return Err(Error::PcaWeightOutOfPlace { column, num_used });
                };
                *target = *weight;
            }
        }
        used_before = used_before
            .checked_add(kept)
            .ok_or_else(the_variants_are_too_many)?;
        num_cols = num_cols
            .checked_add(block.num_vars)
            .ok_or_else(the_variants_are_too_many)?;
    }
    if num_cols != after.num_cols {
        return Err(Error::PcaSecondPassDiffers {
            problem: VariantsOfTheSecondPass::Count {
                found: num_cols,
                expected: after.num_cols,
            },
        });
    }
    Ok(princomps)
}

/// How many variants of the block have variance, once each of them has
/// been found to be one the first pass used and each of the others one it
/// left out.
///
/// `used` says which variants of the block have variance now,
/// `first_position` is the position of its first variant among the
/// variants this pass has given, and `used_before` how many variants that
/// were used the blocks before it gave. The variants of both passes come
/// in the order of the source, so the one to compare each with is the next
/// of the first pass that has not been matched yet.
///
/// # Errors
///
/// [`Error::PcaSecondPassDiffers`] at the first variant that has variance
/// now and had none, or that had variance and has none, and
/// [`Error::PcaVariantsTooLarge`] when the position of a variant is beyond
/// what a `usize` counts.
fn the_variants_of_the_first_pass(
    used: &[bool],
    first_position: usize,
    used_before: usize,
    after: &AfterTheFirstPass<'_>,
) -> Result<usize> {
    let mut kept = 0_usize;
    for (var, was_used) in used.iter().enumerate() {
        let position = first_position
            .checked_add(var)
            .ok_or_else(the_variants_are_too_many)?;
        let of_the_first_pass = used_before
            .checked_add(kept)
            .and_then(|matched| after.used_cols.get(matched));
        if *was_used {
            if of_the_first_pass != Some(&position) {
                return Err(Error::PcaSecondPassDiffers {
                    problem: VariantsOfTheSecondPass::UsedNow(position),
                });
            }
            kept = kept.checked_add(1).ok_or_else(the_variants_are_too_many)?;
        } else if of_the_first_pass == Some(&position) {
            return Err(Error::PcaSecondPassDiffers {
                problem: VariantsOfTheSecondPass::NotUsedNow(position),
            });
        }
    }
    Ok(kept)
}

/// The eigenvectors of the first `num_comps` components divided by sqrt of
/// their eigenvalue, `num_rows` x `num_comps`, row after row.
///
/// A matrix of rows x `num_rows`, the standardized values of some columns
/// of the data, times this gives the weight of each of those columns in
/// each of the components: the weights of the component j are Z' u_j over
/// sqrt(λ_j).
fn the_scaled_vectors_of(eigen: &Eigen, num_rows: usize, num_comps: usize) -> Vec<f64> {
    let mut scaled = vec![0.0; num_values_of(num_rows, num_comps)];
    // A matrix of no component has no value to write, and this keeps
    // `step_by` below off a step of 0, which panics.
    if num_comps == 0 {
        return scaled;
    }
    for (component, (vector, value)) in eigen
        .vectors
        .chunks_exact(num_rows)
        .zip(&eigen.values)
        .take(num_comps)
        .enumerate()
    {
        let size = value.sqrt();
        for (target, coordinate) in scaled
            .iter_mut()
            .skip(component)
            .step_by(num_comps)
            .zip(vector)
        {
            *target = coordinate / size;
        }
    }
    scaled
}

/// What the first pass over the variants leaves, which is all that is kept
/// from one block to the next besides the buffer a block is standardized
/// into.
struct FirstPass {
    /// The lower half of G = Z Z', the individuals x individuals matrix
    /// whose entry i, j is the sum over the variants that were used of the
    /// standardized value of the individual i times that of the individual
    /// j.
    gram: Vec<f64>,
    /// The positions, among the variants the reader gave, of the ones that
    /// were used, in order.
    used_cols: Vec<usize>,
    /// How many variants the reader gave, used or not.
    num_cols: usize,
}

/// The first pass over the variants: each block is standardized into a
/// buffer of the variants that were used x the individuals, and the
/// product of that buffer with itself is added to G.
///
/// The rows of a block are read on the threads of rayon natively and one
/// after another in WebAssembly, and the product is called from outside
/// rayon, as section 3 of `docs/architecture.md` asks. What each block
/// adds to G does not depend on how many threads there are: each row gives
/// its own values and the product is the linear algebra's.
///
/// # Errors
///
/// [`Error::FieldsNotInTheBlock`] when a block holds no genotypes, what
/// the standardizing of a row refuses, [`Error::PcaVariantsTooLarge`] when
/// the variants given are more than a `usize` counts,
/// [`Error::PcaLinalg`] when a product could not be done, and whatever the
/// reader fails with.
fn the_first_pass<R: BlockReader>(
    reader: &mut R,
    options: &VariantPcaOptions,
    num_individuals: usize,
    ploidy: usize,
) -> Result<FirstPass> {
    let mut gram = vec![0.0; num_values_of(num_individuals, num_individuals)];
    let mut used_cols = Vec::new();
    let mut num_cols = 0_usize;
    // The buffer of one block, which is kept from one block to the next so
    // that a pass over a million variants allocates it once. The rows that
    // are not used are left as they were and nothing reads them.
    let mut standardized: Vec<f64> = Vec::new();
    while let Some(block) = reader.next_block()? {
        let used = the_standardized_block(
            &block,
            num_individuals,
            ploidy,
            options,
            num_cols,
            &mut standardized,
        )?;
        for (var, _) in used.iter().enumerate().filter(|(_, was_used)| **was_used) {
            used_cols.push(
                num_cols
                    .checked_add(var)
                    .ok_or_else(the_variants_are_too_many)?,
            );
        }
        let kept = used.iter().filter(|was_used| **was_used).count();
        add_self_product_lower(&standardized, kept, num_individuals, &mut gram).map_err(
            |source| Error::PcaLinalg {
                operation: "product of a block of variants with itself",
                source,
            },
        )?;
        num_cols = num_cols
            .checked_add(block.num_vars)
            .ok_or_else(the_variants_are_too_many)?;
    }
    Ok(FirstPass {
        gram,
        used_cols,
        num_cols,
    })
}

/// The rows of one block standardized into `standardized`, with the rows
/// of the variants that have variance at its start, in the order of the
/// block; it gives whether each variant of the block was used.
///
/// The rows that were left out are not in those first rows, so the product
/// of a block is over its variants that have variance alone. The buffer is
/// the caller's and is kept from one block to the next: it is made as long
/// as the block needs and the rows that are left out keep whatever they
/// held, which nothing reads.
///
/// `first_position` is the position, among the variants the reader has
/// given, of the first variant of the block, which the error of a variant
/// with more than two alleles names.
///
/// # Errors
///
/// [`Error::FieldsNotInTheBlock`] when the block holds no genotypes, what
/// the standardizing of a row refuses, and
/// [`Error::PcaVariantsTooLarge`] when the position of a variant is beyond
/// what a `usize` counts.
fn the_standardized_block(
    block: &Block,
    num_individuals: usize,
    ploidy: usize,
    options: &VariantPcaOptions,
    first_position: usize,
    standardized: &mut Vec<f64>,
) -> Result<Vec<bool>> {
    let missing = Needs::GTS.difference(block.fields());
    if !missing.is_empty() {
        return Err(Error::FieldsNotInTheBlock { fields: missing });
    }
    let alleles_per_var = block.alleles_per_var()?;
    // The block holds its genotypes, so its rows hold one genotype of the
    // ploidy for each individual: `reblock` checked that the genotypes are
    // the variants of the block times those alleles, so this division is
    // exact, and it is `None` only for a ploidy of 0, which such a block
    // does not have.
    let Some(num_values) = block.gts.len().checked_div(ploidy) else {
        return Err(Error::GtsNotWholeGenotypes {
            num_alleles: block.gts.len(),
            ploidy,
        });
    };
    standardized.resize(num_values, 0.0);
    let used = the_standardized_rows(
        &block.gts,
        alleles_per_var,
        num_individuals,
        ploidy,
        options,
        first_position,
        standardized,
    )?;
    // The rows that were used are moved to the start of the buffer. A
    // block with no row to leave out moves nothing.
    for (to, (var, _)) in used
        .iter()
        .enumerate()
        .filter(|(_, was_used)| **was_used)
        .enumerate()
    {
        if to != var {
            let from = the_row_of(var, num_individuals);
            let start = the_row_of(to, num_individuals).start;
            standardized.copy_within(from, start);
        }
    }
    Ok(used)
}

/// The error of a pass that gave more variants than a `usize` counts,
/// which is 4294967295 in WebAssembly.
fn the_variants_are_too_many() -> Error {
    Error::PcaVariantsTooLarge {
        problem: VariantsTooLarge::Variants,
    }
}

/// Where the row `var` of a buffer of rows of `num_individuals` values
/// begins and ends.
///
/// The buffer holds the variants of the block times the individuals
/// values, which the machine gave, and `var` is below the variants of the
/// block, so neither the product nor the sum carries over.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the buffer holds the variants of the block times the individuals values, which the machine gave, and `var` is below the variants of the block"
)]
fn the_row_of(var: usize, num_individuals: usize) -> std::ops::Range<usize> {
    let start = var * num_individuals;
    start..start + num_individuals
}

/// The rows of a block standardized into `standardized`, and whether each
/// variant was used, in the order of the block.
///
/// The rows are read on the threads of rayon, as section 3 of
/// `docs/architecture.md` asks: no row reads another and each one writes
/// its own values, so neither the values nor the variants that are left
/// out depend on how many threads there are. The threads are those of the
/// pool the caller is running in, and rayon's global pool only when the
/// caller is in none. Each thread keeps the buffers of one row and
/// allocates nothing per variant.
///
/// `gts` holds the rows of the block, `alleles_per_var` alleles each, and
/// `alleles_per_var` is 1 or more; `standardized` holds one value for each
/// individual of each of those rows, and a row that is not used is left as
/// it was.
///
/// The error is the one of the first row of the block that has one,
/// wherever the threads found it: each row gives its own result and they
/// are read in the order of the block, so a user who reports a file gets
/// the same message every time.
///
/// # Errors
///
/// What the standardizing of one row refuses, and
/// [`Error::PcaVariantsTooLarge`] when the position of a variant is beyond
/// what a `usize` counts.
#[cfg(not(target_family = "wasm"))]
fn the_standardized_rows(
    gts: &[i8],
    alleles_per_var: usize,
    num_individuals: usize,
    ploidy: usize,
    options: &VariantPcaOptions,
    first_position: usize,
    standardized: &mut [f64],
) -> Result<Vec<bool>> {
    use rayon::iter::{IndexedParallelIterator, ParallelIterator};
    use rayon::slice::{ParallelSlice, ParallelSliceMut};

    let rows: Vec<Result<bool>> = gts
        .par_chunks_exact(alleles_per_var)
        .zip(standardized.par_chunks_exact_mut(num_individuals))
        .enumerate()
        .map_init(
            || RowScratch::of(num_individuals),
            |scratch, (var, (gts, row))| {
                let position = first_position
                    .checked_add(var)
                    .ok_or_else(the_variants_are_too_many)?;
                the_standardized_row(gts, ploidy, position, options, scratch, row)
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
fn the_standardized_rows(
    gts: &[i8],
    alleles_per_var: usize,
    num_individuals: usize,
    ploidy: usize,
    options: &VariantPcaOptions,
    first_position: usize,
    standardized: &mut [f64],
) -> Result<Vec<bool>> {
    the_standardized_rows_one_by_one(
        gts,
        alleles_per_var,
        num_individuals,
        ploidy,
        options,
        first_position,
        standardized,
    )
}

/// The rows read one after another into the buffers of one row: what
/// WebAssembly does, and what the test that compares the two ways of
/// reading a block calls. `alleles_per_var` is 1 or more, as it is for the
/// rows read on threads.
///
/// # Errors
///
/// What the standardizing of one row refuses, and
/// [`Error::PcaVariantsTooLarge`] when the position of a variant is beyond
/// what a `usize` counts.
#[cfg(any(target_family = "wasm", test))]
fn the_standardized_rows_one_by_one(
    gts: &[i8],
    alleles_per_var: usize,
    num_individuals: usize,
    ploidy: usize,
    options: &VariantPcaOptions,
    first_position: usize,
    standardized: &mut [f64],
) -> Result<Vec<bool>> {
    let mut scratch = RowScratch::of(num_individuals);
    let mut used = Vec::new();
    for (var, (gts, row)) in gts
        .chunks_exact(alleles_per_var)
        .zip(standardized.chunks_exact_mut(num_individuals))
        .enumerate()
    {
        let position = first_position
            .checked_add(var)
            .ok_or_else(the_variants_are_too_many)?;
        used.push(the_standardized_row(
            gts,
            ploidy,
            position,
            options,
            &mut scratch,
            row,
        )?);
    }
    Ok(used)
}

/// Where the traits of the table are in the copy of it that is centered
/// and standardized.
///
/// `add_self_product_lower` gives the product of a matrix with itself over
/// its columns, so the side of the table that is to be the side of that
/// product has to be the columns of the copy. The copy is written in the
/// layout that puts the smaller side there, and nothing is transposed
/// afterwards.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Layout {
    /// The traits are the rows of the copy, which is the transpose of the
    /// table. The product is then the rows x rows Z Z', which is what a
    /// table with more traits than rows takes.
    TraitsAsRows,
    /// The traits are the columns of the copy, which is the table as it
    /// came. The product is the traits x traits Z' Z.
    TraitsAsColumns,
}

/// The values of a matrix of `rows` x `cols`.
///
/// Every matrix of a table has one of its two sides the smaller side of
/// the table, or a count of components, which is at most that side, and
/// the other side at most the other side of the table. So every product
/// there is at most the rows of the table times its traits, which [`pca`]
/// took with `checked_mul` before it built any of them. Every matrix of
/// the variants has both of its sides at most the individuals, which
/// [`pca_of_variants`] refuses above [`MAX_INDIVIDUALS_OF_THE_VARIANTS`],
/// whose square is 2147395600, but for the weights, whose values
/// [`the_weights_of_a_second_pass`] takes with `checked_mul` before it
/// builds either of the two buffers that hold them.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "for a table one side is its smaller side or a count of components, at most that side, and the other is at most the other side, so the product is at most the values of the table; for the variants both sides are at most the 46340 individuals the analysis takes, but for the weights, whose values the second pass took with checked_mul and which are more than the buffer of one block of them holds"
)]
fn num_values_of(rows: usize, cols: usize) -> usize {
    rows * cols
}

/// What each of the first `num_comps` components holds of the variance of
/// all the components, as a percentage, from the eigenvalues of the matrix
/// that was decomposed, which come from the largest.
///
/// The share of the total is taken before the 100, so that an eigenvalue
/// above 1.8e306, which a table of values of 1e153 gives, does not become
/// an infinity on the way to a number between 0 and 100.
fn the_percentages_of(values: &[f64], num_comps: usize) -> Vec<f64> {
    let variance_of_every_component: f64 = values.iter().sum();
    values
        .iter()
        .take(num_comps)
        .map(|value| 100.0 * (value / variance_of_every_component))
        .collect()
}

/// The projections when the matrix that was decomposed is the product of
/// the rows, G = Z Z': the eigenvector u of the eigenvalue λ gives the
/// projections u sqrt(λ). `num_rows` x `num_comps`, row after row.
///
/// `eigen` holds the eigenvectors of G as its rows, `num_rows` values
/// each, and their eigenvalues from the largest.
fn the_projections_of(eigen: &Eigen, num_rows: usize, num_comps: usize) -> Vec<f64> {
    let mut projections = vec![0.0; num_values_of(num_rows, num_comps)];
    // A matrix of no component has no value to write, and this keeps
    // `step_by` below off a step of 0, which panics.
    if num_comps == 0 {
        return projections;
    }
    for (component, (vector, value)) in eigen
        .vectors
        .chunks_exact(num_rows)
        .zip(&eigen.values)
        .take(num_comps)
        .enumerate()
    {
        let size = value.sqrt();
        for (projection, coordinate) in projections
            .iter_mut()
            .skip(component)
            .step_by(num_comps)
            .zip(vector)
        {
            *projection = coordinate * size;
        }
    }
    projections
}

/// Refuses the first value of the table that is an infinity or a NaN, with
/// the row and the trait it is at.
///
/// # Errors
///
/// [`Error::PcaValueNotFinite`] with the place of that value.
fn refuse_a_value_that_is_not_finite(table: &[f64], num_cols: usize) -> Result<()> {
    for (row, values) in table.chunks_exact(num_cols).enumerate() {
        for (col, value) in values.iter().enumerate() {
            if !value.is_finite() {
                return Err(Error::PcaValueNotFinite {
                    row,
                    col,
                    value: *value,
                });
            }
        }
    }
    Ok(())
}

/// The mean of each trait and the standard deviation each one is divided
/// by, in the order of the traits.
///
/// A trait that is not centered has a mean of 0 and one that is not
/// standardized a standard deviation of 1, so that one subtraction and one
/// division give the value of every option. The divisor of the standard
/// deviation is the number of rows and not the number of rows less one,
/// which is pyNei's `data.std(axis=0)`.
///
/// # Errors
///
/// [`Error::PcaTraitsWithNoVariance`] when the table is standardized and
/// the values of a trait are all equal. [`Error::PcaTraitOutOfRange`] when
/// the mean or the standard deviation of a trait is not a number the
/// analysis can use, which [`TraitScale`] lists.
fn the_center_and_the_scale_of_each_trait(
    table: &[f64],
    num_rows: usize,
    num_cols: usize,
    options: &PcaOptions,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let rows = num_rows as f64;
    let mut means = vec![0.0; num_cols];
    if options.center {
        for values in table.chunks_exact(num_cols) {
            for (total, value) in means.iter_mut().zip(values) {
                *total += value;
            }
        }
        for mean in &mut means {
            *mean /= rows;
        }
        for (position, mean) in means.iter().enumerate() {
            if !mean.is_finite() {
                return Err(Error::PcaTraitOutOfRange {
                    position,
                    problem: TraitScale::MeanNotFinite,
                });
            }
        }
    }
    let mut deviations = vec![1.0; num_cols];
    if options.standardize {
        let positions = the_traits_with_no_variance(table, num_cols);
        if !positions.is_empty() {
            return Err(Error::PcaTraitsWithNoVariance {
                positions,
                num_cols,
            });
        }
        let mut squares = vec![0.0; num_cols];
        for values in table.chunks_exact(num_cols) {
            for ((total, value), mean) in squares.iter_mut().zip(values).zip(&means) {
                let deviation = value - mean;
                *total += deviation * deviation;
            }
        }
        for (deviation, total) in deviations.iter_mut().zip(&squares) {
            *deviation = (total / rows).sqrt();
        }
        for (position, deviation) in deviations.iter().enumerate() {
            if !deviation.is_finite() {
                return Err(Error::PcaTraitOutOfRange {
                    position,
                    problem: TraitScale::DeviationNotFinite,
                });
            }
            // The traits whose values are all equal were refused above, so
            // a deviation of 0 here is one whose squares were all too
            // small for an `f64` to hold.
            if *deviation == 0.0 {
                return Err(Error::PcaTraitOutOfRange {
                    position,
                    problem: TraitScale::DeviationOfZero,
                });
            }
        }
    }
    Ok((means, deviations))
}

/// The positions of the traits whose values are all equal, in order.
///
/// The values are compared as they are, which is what
/// `docs/specs/pca.md` asks for: pyNei tests `std == 0` on floats, which a
/// trait of 0.1 repeated can pass with a standard deviation of 1e-17, and
/// the division then gives numbers with no meaning.
#[expect(
    clippy::float_cmp,
    reason = "a trait has no variance when its values are equal as they are, and nothing is computed from them first"
)]
fn the_traits_with_no_variance(table: &[f64], num_cols: usize) -> Vec<usize> {
    let mut rows = table.chunks_exact(num_cols);
    // A table of no row has no trait with values that differ, and the
    // caller has refused one of fewer than 2 rows already.
    let Some(first) = rows.next() else {
        return Vec::new();
    };
    let mut equal_to_the_first = vec![true; num_cols];
    for values in rows {
        for ((same, value), reference) in equal_to_the_first.iter_mut().zip(values).zip(first) {
            if *value != *reference {
                *same = false;
            }
        }
    }
    equal_to_the_first
        .iter()
        .enumerate()
        .filter(|(_, same)| **same)
        .map(|(position, _)| position)
        .collect()
}

/// The table centered and divided by the deviations, in the layout the
/// product of the smaller side needs.
fn the_standardized_table(
    table: &[f64],
    num_rows: usize,
    num_cols: usize,
    means: &[f64],
    deviations: &[f64],
    layout: Layout,
) -> Vec<f64> {
    let mut standardized = vec![0.0; num_values_of(num_rows, num_cols)];
    match layout {
        Layout::TraitsAsRows => {
            for (position, ((values, mean), deviation)) in standardized
                .chunks_exact_mut(num_rows)
                .zip(means)
                .zip(deviations)
                .enumerate()
            {
                for (target, value) in values
                    .iter_mut()
                    .zip(table.iter().skip(position).step_by(num_cols))
                {
                    *target = (value - mean) / deviation;
                }
            }
        }
        Layout::TraitsAsColumns => {
            for (values, row) in standardized
                .chunks_exact_mut(num_cols)
                .zip(table.chunks_exact(num_cols))
            {
                for (((target, value), mean), deviation) in
                    values.iter_mut().zip(row).zip(means).zip(deviations)
                {
                    *target = (value - mean) / deviation;
                }
            }
        }
    }
    standardized
}

/// How many of the eigenvalues, which come from the largest, belong to a
/// component that has variance.
///
/// The threshold is the largest eigenvalue times the larger side of the
/// table times 2.220446049250313e-16, the difference between 1 and the
/// next number an `f64` holds, which is the tolerance numpy's
/// `matrix_rank` has for singular values, used here on eigenvalues.
/// `docs/specs/pca.md` has what was measured with it: the eigenvalue of a
/// component with no variance came out between -2e-16 and 3e-16 times the
/// largest on four tables, and 1.3e-14 times it at 1000 x 20000, where the
/// threshold is 4.4e-12 times it.
///
/// The side and the epsilon are multiplied first, so that a largest
/// eigenvalue near the largest `f64`, which a table of values of 2.5e153
/// gives, does not become an infinity on the way to a threshold that is a
/// small part of it.
fn the_components_with_variance(values: &[f64], num_rows: usize, num_cols: usize) -> usize {
    let Some(largest) = values.first() else {
        return 0;
    };
    let threshold = largest * (num_rows.max(num_cols) as f64 * f64::EPSILON);
    values
        .iter()
        .take_while(|value| **value > threshold)
        .count()
}

/// The projections and the weights when the matrix that was decomposed is
/// the product of the traits, Z' Z, whose eigenvectors are the weights
/// themselves and whose projections are Z times them.
///
/// `standardized` is the rows x traits matrix, and the eigenvectors are
/// over the traits.
///
/// # Errors
///
/// [`Error::PcaLinalg`] when the product that gives the projections could
/// not be done.
fn the_components_of_the_product_of_the_traits(
    standardized: &[f64],
    eigen: &Eigen,
    num_rows: usize,
    num_cols: usize,
    num_comps: usize,
) -> Result<(Vec<f64>, Vec<f64>)> {
    // [`pca`] refuses a table with no component before it gets here, and
    // this keeps `step_by` below off a step of 0, which panics.
    if num_comps == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let princomps: Vec<f64> = eigen
        .vectors
        .chunks_exact(num_cols)
        .take(num_comps)
        .flatten()
        .copied()
        .collect();
    let mut weights_by_trait = vec![0.0; num_values_of(num_cols, num_comps)];
    for (component, weights) in princomps.chunks_exact(num_cols).enumerate() {
        for (target, weight) in weights_by_trait
            .iter_mut()
            .skip(component)
            .step_by(num_comps)
            .zip(weights)
        {
            *target = *weight;
        }
    }
    let mut projections = vec![0.0; num_values_of(num_rows, num_comps)];
    product(
        standardized,
        num_rows,
        num_cols,
        TheSecondOperand::ByTheValuesSummedOver {
            values: &weights_by_trait,
            cols: num_comps,
        },
        &mut projections,
    )
    .map_err(|source| Error::PcaLinalg {
        operation: "product that gives the projections",
        source,
    })?;
    Ok((projections, princomps))
}

/// The projections and the weights when the matrix that was decomposed is
/// the product of the rows, Z Z', whose eigenvector u of the eigenvalue λ
/// gives the projections u sqrt(λ) and the weights Z' u / sqrt(λ).
///
/// `standardized` is the traits x rows matrix, which is Z', and the
/// eigenvectors are over the rows.
///
/// # Errors
///
/// [`Error::PcaLinalg`] when the product that gives the weights could not
/// be done.
fn the_components_of_the_product_of_the_rows(
    standardized: &[f64],
    eigen: &Eigen,
    num_rows: usize,
    num_cols: usize,
    num_comps: usize,
) -> Result<(Vec<f64>, Vec<f64>)> {
    // [`pca`] refuses a table with no component before it gets here, and
    // this keeps `step_by` below off a step of 0, which panics.
    if num_comps == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let projections = the_projections_of(eigen, num_rows, num_comps);
    // The same matrix the second pass over the variants multiplies each of
    // its blocks by: the weights of a table and the weights of the
    // variants are the one formula, Z' u over sqrt(λ).
    let vectors_by_row = the_scaled_vectors_of(eigen, num_rows, num_comps);
    let mut weights_by_trait = vec![0.0; num_values_of(num_cols, num_comps)];
    product(
        standardized,
        num_cols,
        num_rows,
        TheSecondOperand::ByTheValuesSummedOver {
            values: &vectors_by_row,
            cols: num_comps,
        },
        &mut weights_by_trait,
    )
    .map_err(|source| Error::PcaLinalg {
        operation: "product that gives the weights",
        source,
    })?;
    let mut princomps = vec![0.0; num_values_of(num_comps, num_cols)];
    for (component, weights) in princomps.chunks_exact_mut(num_cols).enumerate() {
        for (target, weight) in weights
            .iter_mut()
            .zip(weights_by_trait.iter().skip(component).step_by(num_comps))
        {
            *target = *weight;
        }
    }
    Ok((projections, princomps))
}

/// Gives each component the sign of the rule of `docs/specs/pca.md`: the
/// projection of the largest absolute value is positive, and when two
/// individuals have the same absolute value it is the first of them that
/// is made positive. The weights of a component that is turned round are
/// turned round with it.
///
/// `projections` is the individuals x `num_comps` matrix and `princomps`
/// the `num_comps` x `num_cols` one, both row after row.
fn fix_the_signs(
    projections: &mut [f64],
    princomps: &mut [f64],
    num_comps: usize,
    num_cols: usize,
) {
    for component in 0..num_comps {
        // The projection of the largest absolute value, and the first of
        // them when two are of one size, which a value that has to be
        // above the one kept by more than the tolerance keeps.
        let largest = projections
            .iter()
            .skip(component)
            .step_by(num_comps)
            .copied()
            .fold(0.0_f64, |largest: f64, value| {
                if of_one_size(value, largest) == Ordering::Greater {
                    value
                } else {
                    largest
                }
            });
        if largest < 0.0 {
            for projection in projections.iter_mut().skip(component).step_by(num_comps) {
                *projection = -*projection;
            }
            // The weights are given for the first components alone, so
            // there are none for a component after them and there is
            // nothing to turn round: it is not a weight that went missing.
            if let Some(weights) = princomps.chunks_exact_mut(num_cols).nth(component) {
                for weight in weights {
                    *weight = -*weight;
                }
            }
        }
    }
}

/// How the absolute value of `value` compares with that of `largest` for
/// the sign rule: `Greater` when it is above it by more than
/// [`TIED_PROJECTIONS`], and `Equal` when the two are within it, which is
/// what the rule calls the same absolute value.
///
/// `Less` and `Equal` both leave the rule with the projection it had,
/// which is the first of the two, so the rule does not tell them apart.
fn of_one_size(value: f64, largest: f64) -> Ordering {
    let tolerance = largest.abs() * TIED_PROJECTIONS;
    if value.abs() > largest.abs() + tolerance {
        Ordering::Greater
    } else if value.abs() + tolerance < largest.abs() {
        Ordering::Less
    } else {
        Ordering::Equal
    }
}

/// The buffers one thread keeps while it standardizes the rows of a block,
/// so that nothing is allocated for a variant.
struct RowScratch {
    /// The code of the genotype of each individual: its dosage, 0 to the
    /// ploidy, or [`MISSING_CODE`] for a genotype with an allele missing.
    codes: Vec<u8>,
    /// How often each allele was called in the variant, which gives the
    /// major allele and how many different alleles the variant has.
    allele_counts: AlleleCounts,
    /// How many called genotypes have each dosage, 0 to the ploidy.
    dosage_counts: [u32; 255],
    /// The standardized value of each code, which the second pass over
    /// the row looks up: one entry for each value of a byte, so that the
    /// lookup has no bound to check.
    values: [f64; 256],
}

/// The code of a genotype with an allele missing, which is not a dosage:
/// [`MAX_PLOIDY_OF_THE_VARIANTS`] is what keeps the dosages below it.
const MISSING_CODE: u8 = u8::MAX;

/// How many genotypes are counted with one byte of counters at a time. A
/// counter of a byte holds 255, so a run of 255 genotypes is the longest
/// one whose codes a byte counts without wrapping.
const GENOTYPES_PER_RUN: usize = 255;

/// A run longer than what a byte counts would wrap the counters of
/// [`the_counts_of_the_codes`] with nothing to show it, so the length of a
/// run is checked against the largest number a `u8` holds when this
/// compiles.
const _: () = assert!(GENOTYPES_PER_RUN <= 255, "a byte counts to 255");

impl RowScratch {
    /// The buffers of one thread, for the rows of `num_individuals`
    /// individuals.
    fn of(num_individuals: usize) -> RowScratch {
        RowScratch {
            codes: vec![0; num_individuals],
            allele_counts: [0; 128],
            dosage_counts: [0; 255],
            // A genotype with an allele missing takes the mean of the
            // dosages of its variant, which is 0 once the variant is
            // centered, and this entry is never written again.
            values: [0.0; 256],
        }
    }
}

/// One row of a block standardized into `row`, and whether the variant was
/// used: a variant whose called genotypes all have one dosage, and one
/// with no called genotype, have no variance and are left out, and `row`
/// is then left as it was.
///
/// `gts` is the genotypes of one variant, `ploidy` alleles for each
/// individual, and `row` holds one value for each individual. `position`
/// is which variant of those the reader gave this one is, which the error
/// of a variant with more than two alleles names.
///
/// # Errors
///
/// [`Error::PcaVariantWithMoreThanTwoAlleles`] when the variant has more
/// than two different alleles among its called genotypes and
/// `transform_to_biallelic` is false, and whatever the counts of the
/// alleles of one variant refuse.
fn the_standardized_row(
    gts: &[i8],
    ploidy: usize,
    position: usize,
    options: &VariantPcaOptions,
    scratch: &mut RowScratch,
    row: &mut [f64],
) -> Result<bool> {
    if ploidy > MAX_PLOIDY_OF_THE_VARIANTS {
        return Err(Error::PcaVariantsTooLarge {
            problem: VariantsTooLarge::Ploidy(ploidy),
        });
    }
    let num_individuals = row.len();
    // One genotype of the ploidy for each individual, which is what a row
    // of the genotypes of a block holds. The ploidy of 0 that
    // `NonZeroUsize` refuses is among these: it would be a genotype of no
    // allele for every individual.
    let of_a_genotype = match NonZeroUsize::new(ploidy) {
        Some(of_a_genotype) if gts.len() == num_individuals.saturating_mul(ploidy) => of_a_genotype,
        _ => {
            return Err(Error::GtsNotWholeGenotypes {
                num_alleles: gts.len(),
                ploidy,
            });
        }
    };
    let RowScratch {
        codes,
        allele_counts,
        dosage_counts,
        values,
    } = scratch;
    let called_alleles = count_alleles(gts, allele_counts)?;
    if called_alleles == 0 {
        // A variant with no called genotype has no dosage at all, so it
        // has no variance and is left out, as Open 2 of
        // `docs/specs/pca.md` has it meanwhile; pyNei makes it an error
        // and leaves a variant with one allele out in silence.
        return Ok(false);
    }
    let num_alleles = allele_counts.iter().filter(|count| **count > 0).count();
    if num_alleles > 2 && !options.transform_to_biallelic {
        return Err(Error::PcaVariantWithMoreThanTwoAlleles {
            position,
            num_alleles,
        });
    }
    // The buffer of the codes belongs to the thread and is as long as the
    // rows it has read so far, which are all of the individuals of the
    // source: this asks the machine for nothing after the first row.
    codes.resize(num_individuals, 0);
    the_codes_of_the_genotypes(gts, of_a_genotype, the_major_allele(allele_counts), codes);
    // The dosages of a genotype are 0 to the ploidy, and the counts have
    // one entry for each of them.
    let num_dosages = ploidy.saturating_add(1);
    the_counts_of_the_codes(codes, num_dosages, dosage_counts);
    let Some((mean, deviation)) =
        the_center_and_the_scale_of_the_dosages(dosage_counts, num_dosages, num_individuals)
    else {
        return Ok(false);
    };
    // What each code is worth, which the pass over the row below looks up.
    // The entry of a genotype with an allele missing is the 0 the buffer
    // was built with and is never written: a dosage is below 255, which
    // the ploidy of 254 at most is what keeps.
    for (value, dosage) in values.iter_mut().take(num_dosages).zip(0_u32..) {
        *value = (f64::from(dosage) - mean) / deviation;
    }
    for (target, code) in row.iter_mut().zip(codes.iter()) {
        // The values have an entry for each of the 256 values of a byte,
        // so every code has one and the 0.0 is never taken; it is also
        // what lets the compiler drop the bound of the lookup.
        *target = values.get(usize::from(*code)).copied().unwrap_or(0.0);
    }
    Ok(true)
}

/// The code of the genotype of each individual: its dosage, how many of
/// its alleles are not the major one, or [`MISSING_CODE`] when one allele
/// of it at least was not called.
///
/// `gts` holds one genotype of `of_a_genotype` alleles for each individual
/// and `codes` one byte for each. It is the first of the two passes over a
/// row that `docs/specs/pca.md` writes for the compiler to turn into
/// vector instructions.
///
/// The ploidies 1 to 4 each get the loop with the length of a genotype
/// written into it, because with that length a constant the compiler reads
/// the alleles of several genotypes at once, and with it a number the
/// dataset carries it reads one allele at a time. Every other ploidy takes
/// the loop that reads the length from the dataset. The arms give the same
/// codes: the dosage is a count of alleles and the missing flag is a
/// boolean, and [`the_codes_of_any_ploidy`] is the same body with the same
/// length.
fn the_codes_of_the_genotypes(
    gts: &[i8],
    of_a_genotype: NonZeroUsize,
    major: i8,
    codes: &mut [u8],
) {
    match of_a_genotype.get() {
        1 => the_codes_of_a_ploidy_of::<1>(gts, major, codes),
        2 => the_codes_of_a_ploidy_of::<2>(gts, major, codes),
        3 => the_codes_of_a_ploidy_of::<3>(gts, major, codes),
        4 => the_codes_of_a_ploidy_of::<4>(gts, major, codes),
        of_a_genotype => the_codes_of_any_ploidy(gts, of_a_genotype, major, codes),
    }
}

/// The codes of the genotypes of a row whose genotypes hold `OF_A_GENOTYPE`
/// alleles, which is the body of [`the_codes_of_the_genotypes`] with the
/// length of a genotype known when the code is compiled.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the dosage counts the alleles of one genotype, which are the ploidy, and the caller refused a ploidy above 254"
)]
fn the_codes_of_a_ploidy_of<const OF_A_GENOTYPE: usize>(gts: &[i8], major: i8, codes: &mut [u8]) {
    let (genotypes, _) = gts.as_chunks::<OF_A_GENOTYPE>();
    for (code, genotype) in codes.iter_mut().zip(genotypes) {
        let mut dosage = 0_u8;
        let mut missing = 0_u8;
        for allele in genotype {
            dosage += u8::from(*allele != major);
            missing |= u8::from(*allele == MISSING_ALLELE);
        }
        *code = if missing == 0 { dosage } else { MISSING_CODE };
    }
}

/// The codes of the genotypes of a row whose genotypes hold
/// `of_a_genotype` alleles, a length the dataset carries.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the dosage counts the alleles of one genotype, which are the ploidy, and the caller refused a ploidy above 254"
)]
fn the_codes_of_any_ploidy(gts: &[i8], of_a_genotype: usize, major: i8, codes: &mut [u8]) {
    for (code, genotype) in codes.iter_mut().zip(gts.chunks_exact(of_a_genotype)) {
        let mut dosage = 0_u8;
        let mut missing = 0_u8;
        for allele in genotype {
            dosage += u8::from(*allele != major);
            missing |= u8::from(*allele == MISSING_ALLELE);
        }
        *code = if missing == 0 { dosage } else { MISSING_CODE };
    }
}

/// How many genotypes have each dosage, 0 to the ploidy, written into the
/// first `num_dosages` entries of `counts`. A genotype with an allele
/// missing has no dosage and is counted in none of them.
///
/// The codes are counted in runs of [`GENOTYPES_PER_RUN`] genotypes with
/// counters of one byte, one for each dosage: a counter of a byte counts a
/// run whole without wrapping, and each pass over a run compares the code
/// with one dosage and adds, which is what the compiler turns into vector
/// instructions.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "a run holds 255 codes at most, so a counter of one byte counts it without wrapping, and each total counts the genotypes of the variant, which the counts of its alleles checked to be a number a u32 holds"
)]
fn the_counts_of_the_codes(codes: &[u8], num_dosages: usize, counts: &mut [u32; 255]) {
    for total in counts.iter_mut().take(num_dosages) {
        *total = 0;
    }
    for run in codes.chunks(GENOTYPES_PER_RUN) {
        for (total, dosage) in counts.iter_mut().take(num_dosages).zip(0_u8..) {
            let mut count = 0_u8;
            for code in run {
                count += u8::from(*code == dosage);
            }
            *total += u32::from(count);
        }
    }
}

/// The mean of the called dosages of a variant and the standard deviation
/// it is divided by, or `None` when its called genotypes all have one
/// dosage and it has no variance.
///
/// The divisor of the deviation is all the individuals and not the called
/// genotypes, because the genotypes with an allele missing are in the
/// column with a deviation of 0: a variant with much missing data has a
/// smaller deviation for the same frequencies, and its called genotypes
/// weigh more. The divisor is the individuals and not the individuals less
/// one, which is pyNei's `data.std(axis=0)`.
///
/// A variant that has variance gets a deviation that is finite and above
/// 0, so the values of the buffer a block is standardized into are finite,
/// which is what the product of that buffer with itself asks for: the
/// dosages are whole numbers from 0 to 254, so their sum is exact in an
/// `f64` and the squares of their deviations are at most 64516 times the
/// individuals; and two dosages that differ by 1, the least a variant with
/// variance has, give squares of 1 over the called genotypes at least.
///
/// `counts` holds how many genotypes have each dosage, which
/// [`the_counts_of_the_codes`] wrote, and `num_individuals` is 1 or more.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the called genotypes are at most the individuals of the variant, which the counts of its alleles checked to be a number a u32 holds, and each dosage is 254 at most, so the sum of the dosages is below 2^53"
)]
fn the_center_and_the_scale_of_the_dosages(
    counts: &[u32; 255],
    num_dosages: usize,
    num_individuals: usize,
) -> Option<(f64, f64)> {
    let mut called = 0_u64;
    let mut total = 0_u64;
    let mut dosages_seen = 0_usize;
    for (count, dosage) in counts.iter().take(num_dosages).zip(0_u64..) {
        if *count == 0 {
            continue;
        }
        dosages_seen += 1;
        called += u64::from(*count);
        total += dosage * u64::from(*count);
    }
    // One dosage among the called genotypes is a variant with no variance:
    // one with one allele, and one where every individual is heterozygous,
    // whose major allele frequency is 0.5 and which no filter by frequency
    // would catch. It is decided on the counts, with no float in it, where
    // `_remove_vars_with_no_variance` of pyNei tests `std > 0`.
    if dosages_seen < 2 {
        return None;
    }
    let mean = total as f64 / called as f64;
    let mut squares = 0.0_f64;
    for (count, dosage) in counts.iter().take(num_dosages).zip(0_u64..) {
        let deviation = dosage as f64 - mean;
        squares += f64::from(*count) * deviation * deviation;
    }
    Some((mean, (squares / num_individuals as f64).sqrt()))
}

/// The positions of the first ten traits of a list, as text, with how many
/// more there are: `the position 3`, `the positions 0, 3, 7`, or `the
/// positions 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, and 2 more`.
///
/// The word comes with the numbers because a message that ends in one
/// number, `the traits at 1`, reads as a count of traits to whoever gets
/// it in TypeScript, where no name replaces it.
/// [`Error::PcaTraitsWithNoVariance`] carries the positions and the Python
/// layer puts the name of each trait in their place, as pyNei's message
/// has it, so this is what a reader of the message in Rust or in
/// TypeScript gets.
pub(crate) fn the_positions_listed(positions: &[usize]) -> String {
    let word = if positions.len() == 1 {
        "the position"
    } else {
        "the positions"
    };
    let shown = positions
        .iter()
        .take(10)
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    match positions.len().saturating_sub(10) {
        0 => format!("{word} {shown}"),
        more => format!("{word} {shown}, and {more} more"),
    }
}

/// What the benchmark `standardize_row` calls to time each pass over a row
/// of a block on its own.
///
/// The four passes that standardizing a row is made of are private to this
/// module, and a benchmark is a crate of its own, so nothing outside can
/// call them; the compiler inlines all four into one closure, so a
/// sampling profile of the analysis sees them as one frame. This module is
/// behind the cargo feature `bench-internals`, which is off by default and
/// which nothing of popnei's own builds turn on, and it holds one wrapper
/// for each of them, with the arguments the private function takes, so
/// that what the benchmark times is the code the library runs.
///
/// [`the_standardized_values`] is the one thing here that is not a
/// wrapper. The tail of [`the_standardized_row`] is two loops written
/// inline in that function and not a function of its own, so timing it
/// needs those two loops copied here; a change to them there is a change
/// to them here.
#[cfg(feature = "bench-internals")]
#[doc(hidden)]
pub mod bench_internals {
    use std::num::NonZeroUsize;

    use super::{
        RowScratch, VariantPcaOptions, the_center_and_the_scale_of_the_dosages,
        the_codes_of_the_genotypes as codes_of_the_genotypes,
        the_counts_of_the_codes as counts_of_the_codes, the_standardized_row as standardized_row,
    };
    use crate::error::Result;

    /// The buffers one thread keeps while it standardizes the rows of a
    /// block, so that nothing is allocated for a variant. It is
    /// `RowScratch`, which is private, and it is opaque here: the
    /// benchmark builds one and hands it back.
    pub struct Scratch(RowScratch);

    impl Scratch {
        /// The buffers for the rows of `num_individuals` individuals.
        #[must_use]
        pub fn of(num_individuals: usize) -> Scratch {
            Scratch(RowScratch::of(num_individuals))
        }
    }

    /// One row of a block standardized into `row`, and whether the variant
    /// was used, which is `the_standardized_row` of this module and all
    /// four passes over the row together.
    ///
    /// # Errors
    ///
    /// What `the_standardized_row` gives: a ploidy the analysis refuses,
    /// genotypes that are not whole, a variant of more than two alleles.
    pub fn the_standardized_row(
        gts: &[i8],
        ploidy: usize,
        position: usize,
        options: &VariantPcaOptions,
        scratch: &mut Scratch,
        row: &mut [f64],
    ) -> Result<bool> {
        standardized_row(gts, ploidy, position, options, &mut scratch.0, row)
    }

    /// The code of the genotype of each individual, its dosage or the code
    /// of a genotype with an allele missing, which is
    /// `the_codes_of_the_genotypes` of this module.
    pub fn the_codes_of_the_genotypes(
        gts: &[i8],
        of_a_genotype: NonZeroUsize,
        major: i8,
        codes: &mut [u8],
    ) {
        codes_of_the_genotypes(gts, of_a_genotype, major, codes);
    }

    /// How many genotypes have each dosage, which is
    /// `the_counts_of_the_codes` of this module.
    pub fn the_counts_of_the_codes(codes: &[u8], num_dosages: usize, counts: &mut [u32; 255]) {
        counts_of_the_codes(codes, num_dosages, counts);
    }

    /// The tail of [`the_standardized_row`]: the value of each dosage
    /// written into `values`, and then the value of each code of the row
    /// looked up there and written into `row`. It gives whether the
    /// variant had variance, and leaves `row` as it was when it had none.
    ///
    /// These are the two loops that function ends in, copied, because they
    /// are written inline there and not called.
    pub fn the_standardized_values(
        dosage_counts: &[u32; 255],
        num_dosages: usize,
        codes: &[u8],
        values: &mut [f64; 256],
        row: &mut [f64],
    ) -> bool {
        let Some((mean, deviation)) =
            the_center_and_the_scale_of_the_dosages(dosage_counts, num_dosages, row.len())
        else {
            return false;
        };
        for (value, dosage) in values.iter_mut().take(num_dosages).zip(0_u32..) {
            *value = (f64::from(dosage) - mean) / deviation;
        }
        for (target, code) in row.iter_mut().zip(codes.iter()) {
            *target = values.get(usize::from(*code)).copied().unwrap_or(0.0);
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::num::NonZeroUsize;
    use std::path::{Path, PathBuf};

    use popnei_linalg::eigh_lower;

    use super::{
        AfterTheFirstPass, FirstPass, MAX_INDIVIDUALS_OF_THE_VARIANTS, MAX_PLOIDY_OF_THE_VARIANTS,
        MISSING_CODE, Pca, PcaOptions, RowScratch, TraitScale, VariantPcaOptions,
        VariantsOfTheSecondPass, VariantsTooLarge, fix_the_signs, pca, pca_of_variants,
        the_codes_of_any_ploidy, the_codes_of_the_genotypes, the_components_with_variance,
        the_first_pass, the_scaled_vectors_of, the_standardized_row, the_weights_of_a_second_pass,
    };
    // The two ways of reading the rows of a block are one function in
    // WebAssembly, which has no threads, so the test that compares them,
    // and what only it uses, are of the targets that have them.
    #[cfg(not(target_family = "wasm"))]
    use super::{the_row_of, the_standardized_rows, the_standardized_rows_one_by_one};
    use crate::block::{Block, BlockReader, Reblock};
    use crate::error::{Error, Result};
    use crate::filters::FilteringStats;
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{ChromTable, MISSING_ALLELE, Needs};

    /// The tolerance of "How it is verified" of `docs/specs/pca.md`: every
    /// literal here and in the reference files is written with 12
    /// significant digits, and the largest of them is a projection of iris
    /// of 3.7956, so a literal is within 1e-11 of the number it stands
    /// for.
    const TOLERANCE: f64 = 1e-9;

    /// Centered and standardized, which is what `do_pca` does by default.
    const STANDARDIZED: PcaOptions = PcaOptions {
        center: true,
        standardize: true,
    };

    /// Centered and not standardized.
    const CENTERED: PcaOptions = PcaOptions {
        center: true,
        standardize: false,
    };

    /// Neither centered nor standardized.
    const AS_IT_IS: PcaOptions = PcaOptions {
        center: false,
        standardize: false,
    };

    /// The 3 rows x 5 traits of "How it is verified" of "The PCA of a
    /// table" of `docs/specs/pca.md`, which has fewer rows than traits, so
    /// that the product is the 3 x 3 one of the rows.
    const FIVE_TRAITS: [f64; 15] = [
        1.0, 2.0, 3.0, 4.0, 5.0, //
        2.0, 4.0, 1.0, 3.0, 2.0, //
        5.0, 1.0, 4.0, 2.0, 6.0,
    ];

    /// The table of `test_pca_refuses_traits_with_no_variance` of pyNei,
    /// whose traits `a`, `fixed` and `b` are the columns, the second of
    /// them having no variance.
    const ONE_TRAIT_FIXED: [f64; 9] = [
        1.0, 5.0, 3.0, //
        2.0, 5.0, 1.0, //
        3.0, 5.0, 2.0,
    ];

    /// The 5 rows x 3 traits of "How it is verified", which has more rows
    /// than traits, as iris has, and loses a component, which iris does
    /// not: its second trait has no variance.
    const FIVE_ROWS: [f64; 15] = [
        1.0, 5.0, 3.0, //
        2.0, 5.0, 1.0, //
        3.0, 5.0, 2.0, //
        4.0, 5.0, 9.0, //
        7.0, 5.0, 2.0,
    ];

    /// The 3 rows x 3 traits of "How it is verified" whose values are
    /// multiplied by 1e153 or by 2.5e153, which puts the eigenvalues of
    /// its product near the largest `f64`.
    const NEAR_THE_LARGEST: [f64; 9] = [
        1.0, 2.0, 3.0, //
        2.0, 4.0, 1.0, //
        5.0, 1.0, 4.0,
    ];

    /// The rows of one of the files of `tests/reference/pca/`, without the
    /// first field of each line, which names the row, and without the
    /// first line when the file has a header. The values of R's files
    /// carry leading spaces.
    fn the_reference(name: &str, with_a_header: bool) -> Vec<Vec<f64>> {
        let path = the_reference_path(name);
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("{path}: {error}", path = path.display()));
        text.lines()
            .skip(usize::from(with_a_header))
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                line.split('\t')
                    .skip(1)
                    .map(|field| {
                        field
                            .trim()
                            .parse()
                            .unwrap_or_else(|error| panic!("{name}: `{field}`: {error}"))
                    })
                    .collect()
            })
            .collect()
    }

    /// The same rows, one after another, which is how [`Pca`] holds a
    /// matrix.
    fn row_after_row(rows: Vec<Vec<f64>>) -> Vec<f64> {
        rows.into_iter().flatten().collect()
    }

    /// Each value against the one of the reference, within the tolerance.
    fn assert_close(got: &[f64], expected: &[f64], tolerance: f64, what: &str) {
        assert_eq!(got.len(), expected.len(), "{what}: the count of values");
        for (position, (value, reference)) in got.iter().zip(expected).enumerate() {
            assert!(
                (value - reference).abs() <= tolerance,
                "{what}: the value at {position} is {value} and the reference is {reference}"
            );
        }
    }

    /// What R's `prcomp` gives for iris, from the three files of that
    /// name: the projections, the percentages and the weights, each row
    /// after row.
    fn the_iris_reference(name: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        (
            row_after_row(the_reference(&format!("{name}.r.projections.tsv"), true)),
            row_after_row(the_reference(&format!("{name}.r.percent.tsv"), false)),
            row_after_row(the_reference(&format!("{name}.r.princomps.tsv"), true)),
        )
    }

    /// Every component of the result against a reference held row after
    /// row, with the shape of the result.
    fn assert_the_result_is(
        result: &Pca,
        num_rows: usize,
        num_cols: usize,
        num_comps: usize,
        reference: (&[f64], &[f64], &[f64]),
        what: &str,
    ) {
        assert_eq!(result.num_rows, num_rows, "{what}: the rows");
        assert_eq!(result.num_cols, num_cols, "{what}: the traits");
        assert_eq!(result.num_comps, num_comps, "{what}: the components");
        assert_eq!(
            result.num_prin_comps, num_comps,
            "{what}: the components the weights are given for"
        );
        assert_eq!(
            result.used_cols,
            (0..num_cols).collect::<Vec<usize>>(),
            "{what}: every trait of a table is used"
        );
        assert_close(
            &result.projections,
            reference.0,
            TOLERANCE,
            &format!("{what}: the projections"),
        );
        assert_close(
            &result.explained_variance_percent,
            reference.1,
            TOLERANCE,
            &format!("{what}: the percentages"),
        );
        assert_close(
            &result.princomps,
            reference.2,
            TOLERANCE,
            &format!("{what}: the weights"),
        );
    }

    /// The table of iris, 150 rows x 4 traits, which
    /// `tests/reference/pca/make_reference.py` writes from pyNei's
    /// `test/datasets.py`.
    fn the_iris_table() -> Vec<f64> {
        row_after_row(the_reference("iris.tsv", true))
    }

    /// Iris standardized is the first of the two runs of "How it is
    /// verified" of "The PCA of a table", and it has 150 rows and 4
    /// traits, so its product is the 4 x 4 one of the traits. R's
    /// projections in the file are multiplied by sqrt(n / (n - 1)), so
    /// they are pyNei's and popnei's, the ones of a standard deviation
    /// with n in it.
    #[test]
    fn iris_standardized_gives_the_numbers_of_r() {
        let result = pca(&the_iris_table(), 150, 4, &STANDARDIZED).expect("the analysis of iris");
        let (projections, percent, princomps) = the_iris_reference("iris");
        assert_the_result_is(
            &result,
            150,
            4,
            4,
            (&projections, &percent, &princomps),
            "iris standardized",
        );
    }

    /// The second run of the same part, which R computes with
    /// `scale. = FALSE` and whose projections are R's as they are.
    #[test]
    fn iris_not_standardized_gives_the_numbers_of_r() {
        let result = pca(&the_iris_table(), 150, 4, &CENTERED).expect("the analysis of iris");
        let (projections, percent, princomps) = the_iris_reference("iris_not_standardized");
        assert_the_result_is(
            &result,
            150,
            4,
            4,
            (&projections, &percent, &princomps),
            "iris not standardized",
        );
    }

    /// The table of 3 rows x 5 traits of "How it is verified", which has
    /// fewer rows than traits, so that the matrix that is decomposed is
    /// the 3 x 3 product of the rows and not the 5 x 5 product of the
    /// traits. Centering takes one of the three dimensions of the rows
    /// out, so it has 2 components and not 3.
    #[test]
    fn a_table_of_fewer_rows_than_traits_standardized_gives_the_numbers_of_numpy() {
        let result = pca(&FIVE_TRAITS, 3, 5, &STANDARDIZED).expect("the analysis of the table");
        let projections = [
            -0.347154191646,
            1.62410767053,
            -2.14982434306,
            -0.994055030603,
            2.49697853471,
            -0.630052639924,
        ];
        let percent = [73.1810836106, 26.8189163894];
        let princomps = [
            0.420101972196,
            -0.496432942021,
            0.496432942021,
            -0.317325807736,
            0.47950738561,
            -0.513968704337,
            -0.270671721206,
            0.270671721206,
            0.686274627813,
            0.344001373342,
        ];
        assert_the_result_is(
            &result,
            3,
            5,
            2,
            (&projections, &percent, &princomps),
            "the table of five traits, standardized",
        );
    }

    /// The same table with neither step, which keeps all 3 components,
    /// since nothing takes a dimension out. It is the run that no program
    /// outside the project gives a number for.
    #[test]
    fn a_table_that_is_neither_centered_nor_standardized_keeps_every_component() {
        let result = pca(&FIVE_TRAITS, 3, 5, &AS_IT_IS).expect("the analysis of the table");
        let projections = [
            7.07789159873,
            1.12190772301,
            1.90912901021,
            4.87213205354,
            2.87255560109,
            -1.41801042717,
            8.66041659868,
            -2.53292797372,
            -0.762535387595,
        ];
        let percent = [87.0392022767, 9.31343669027, 3.64736103304];
        let princomps = [
            0.403960199411,
            0.28423522248,
            0.408147561391,
            0.404797068091,
            0.652365999566,
            -0.364035502371,
            0.703323259808,
            -0.244470602227,
            0.504800545607,
            -0.241298733997,
            -0.759913160568,
            -0.419484427674,
            0.201897964367,
            0.297806276491,
            0.342218405411,
        ];
        assert_the_result_is(
            &result,
            3,
            5,
            3,
            (&projections, &percent, &princomps),
            "the table of five traits, as it is",
        );
    }

    /// pyNei's table of three traits, one of which is fixed, without
    /// standardizing, which pyNei gives 3 components and popnei 2: the
    /// third has no variance. The fixed trait is still a column of the
    /// weights, with a weight of 0.
    ///
    /// The second component is given up to its sign, as "How it is
    /// verified" says: its two largest projections are the same number
    /// with opposite signs, equal bit for bit on Accelerate's LAPACK and
    /// one bit apart on faer, so which of the two the sign rule finds is
    /// decided by the rounding of the eigendecomposition and is not the
    /// same on every backend.
    #[test]
    #[expect(
        clippy::approx_constant,
        reason = "the numbers of the spec, written with the 12 significant digits of the reference, and four of them are the first digits of sqrt(2) and of its inverse"
    )]
    fn a_trait_with_no_variance_gets_a_weight_of_0_and_no_component_of_its_own() {
        let result =
            pca(&ONE_TRAIT_FIXED, 3, 3, &CENTERED).expect("the analysis of the fixed trait");
        assert_eq!(result.num_comps, 2);
        assert_eq!(result.num_cols, 3);
        assert_eq!(result.used_cols, vec![0, 1, 2]);
        let first_of_the_projections = [
            result.projections.first().copied().expect("row 0"),
            result.projections.get(2).copied().expect("row 1"),
            result.projections.get(4).copied().expect("row 2"),
        ];
        assert_close(
            &first_of_the_projections,
            &[1.41421356237, -0.707106781187, -0.707106781187],
            TOLERANCE,
            "the projections of the first component",
        );
        assert_close(
            &result.explained_variance_percent,
            &[75.0, 25.0],
            TOLERANCE,
            "the percentages",
        );
        assert_close(
            result.princomps.get(..3).expect("the first component"),
            &[-0.707106781187, 0.0, 0.707106781187],
            TOLERANCE,
            "the weights of the first component",
        );
        // The two largest projections of the second component are the same
        // number with opposite signs, equal bit for bit on Accelerate's
        // LAPACK, one bit apart on faer and two apart on numpy, all of
        // which are inside the tolerance of the sign rule: the first of
        // the two is the positive one on every backend, so the component
        // has one sign and not one sign for each.
        let second: Vec<f64> = result
            .projections
            .iter()
            .skip(1)
            .step_by(2)
            .copied()
            .collect();
        assert_close(
            &second,
            &[0.0, 0.707106781187, -0.707106781187],
            TOLERANCE,
            "the projections of the second component",
        );
        assert_close(
            result.princomps.get(3..6).expect("the second component"),
            &[-0.707106781187, 0.0, -0.707106781187],
            TOLERANCE,
            "the weights of the second component",
        );
    }

    /// The sign rule of "What both analyses compute": in every component
    /// the projection of the largest absolute value is positive. Whichever
    /// sign the eigendecomposition gave, this is what comes out, so that
    /// Python natively, Python under pyodide and TypeScript give the same
    /// numbers.
    #[test]
    fn the_projection_of_the_largest_absolute_value_is_positive_in_every_component() {
        for (data, num_rows, num_cols, options, what) in [
            (&the_iris_table()[..], 150, 4, STANDARDIZED, "iris"),
            (&FIVE_TRAITS[..], 3, 5, STANDARDIZED, "five traits"),
            (&FIVE_TRAITS[..], 3, 5, AS_IT_IS, "five traits as they are"),
            (&ONE_TRAIT_FIXED[..], 3, 3, CENTERED, "one trait fixed"),
        ] {
            let result = pca(data, num_rows, num_cols, &options).expect("the analysis");
            assert!(result.num_comps > 0, "{what}: no component");
            for component in 0..result.num_comps {
                // The first projection of the largest size, where two
                // within the tolerance of the rule are of one size.
                let largest = result
                    .projections
                    .iter()
                    .skip(component)
                    .step_by(result.num_comps)
                    .copied()
                    .fold(0.0_f64, |largest: f64, value| {
                        if super::of_one_size(value, largest) == std::cmp::Ordering::Greater {
                            value
                        } else {
                            largest
                        }
                    });
                assert!(
                    largest > 0.0,
                    "{what}: the first projection of the largest size of component {component} is {largest}"
                );
            }
        }
    }

    /// The other half of the rule: when two individuals have the same
    /// absolute value it is the first of them that is made positive. The
    /// two projections here are exactly opposite, which the
    /// eigendecomposition of a table seldom gives, so the rule is checked
    /// on the step that applies it.
    #[test]
    fn the_sign_rule_breaks_a_tie_towards_the_first_individual() {
        let mut projections = vec![1.0, -1.0];
        let mut princomps = vec![0.5, -0.25];
        fix_the_signs(&mut projections, &mut princomps, 1, 2);
        assert_close(&projections, &[1.0, -1.0], 0.0, "the first is positive");
        assert_close(&princomps, &[0.5, -0.25], 0.0, "the weights are untouched");

        let mut projections = vec![-1.0, 1.0];
        let mut princomps = vec![0.5, -0.25];
        fix_the_signs(&mut projections, &mut princomps, 1, 2);
        assert_close(&projections, &[1.0, -1.0], 0.0, "the component is turned");
        assert_close(&princomps, &[-0.5, 0.25], 0.0, "the weights are turned");
    }

    /// A value that is not finite is refused with the place where it is.
    /// pyNei refuses a NaN and lets an infinity reach numpy's SVD, which
    /// raises `LinAlgError`.
    #[test]
    fn a_value_that_is_not_finite_is_refused_with_its_place() {
        let mut table = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        *table.get_mut(5).expect("the last value") = f64::INFINITY;
        match pca(&table, 2, 3, &CENTERED) {
            Err(Error::PcaValueNotFinite { row, col, value }) => {
                assert_eq!(row, 1);
                assert_eq!(col, 2);
                assert!(value.is_infinite(), "{value}");
            }
            other => panic!("an infinity was taken: {other:?}"),
        }

        *table.get_mut(5).expect("the last value") = 6.0;
        *table.get_mut(1).expect("the second value") = f64::NAN;
        match pca(&table, 2, 3, &CENTERED) {
            Err(Error::PcaValueNotFinite { row, col, value }) => {
                assert_eq!(row, 0);
                assert_eq!(col, 1);
                assert!(value.is_nan(), "{value}");
            }
            other => panic!("a NaN was taken: {other:?}"),
        }
    }

    /// Standardizing divides by the standard deviation of the centered
    /// trait, so it needs the centering. pyNei refuses the pair too.
    #[test]
    fn standardizing_without_centering_is_refused() {
        let options = PcaOptions {
            center: false,
            standardize: true,
        };
        let result = pca(&FIVE_TRAITS, 3, 5, &options);
        assert!(
            matches!(result, Err(Error::PcaStandardizeWithoutCentering)),
            "{result:?}"
        );
    }

    /// One row has no variation to find directions in, and no trait has
    /// nothing to project.
    #[test]
    fn a_table_of_fewer_than_two_rows_or_of_no_traits_is_refused() {
        for (num_rows, num_cols) in [(1_usize, 3_usize), (0, 3), (3, 0)] {
            let result = pca(&ONE_TRAIT_FIXED, num_rows, num_cols, &CENTERED);
            match result {
                Err(Error::PcaTableTooSmall {
                    num_rows: rows,
                    num_cols: cols,
                }) => {
                    assert_eq!(rows, num_rows);
                    assert_eq!(cols, num_cols);
                }
                other => panic!("a table of {num_rows} x {num_cols} was taken: {other:?}"),
            }
        }
    }

    /// The traits that have no variance are refused when the table is
    /// standardized, with the position of each one, which the Python layer
    /// turns into its name. The message names the first ten and says how
    /// many more there are, as pyNei's does.
    #[test]
    fn the_traits_with_no_variance_are_refused_with_their_positions() {
        match pca(&ONE_TRAIT_FIXED, 3, 3, &STANDARDIZED) {
            Err(Error::PcaTraitsWithNoVariance {
                positions,
                num_cols,
            }) => {
                assert_eq!(positions, vec![1]);
                assert_eq!(num_cols, 3);
                let message = Error::PcaTraitsWithNoVariance {
                    positions,
                    num_cols,
                }
                .to_string();
                assert!(message.contains("1 of the 3 traits"), "{message}");
                assert!(message.contains("no variance"), "{message}");
                // The position comes with the word, so that a message
                // that ends in one number does not read as a count.
                assert!(message.contains("at the position 1,"), "{message}");
            }
            other => panic!("a fixed trait was standardized: {other:?}"),
        }

        // Twelve traits of two rows, all of them fixed, which is what a
        // message that names the first ten needs.
        let table = vec![7.0; 24];
        match pca(&table, 2, 12, &STANDARDIZED) {
            Err(Error::PcaTraitsWithNoVariance {
                positions,
                num_cols,
            }) => {
                assert_eq!(positions, (0..12).collect::<Vec<usize>>());
                assert_eq!(num_cols, 12);
                let message = Error::PcaTraitsWithNoVariance {
                    positions,
                    num_cols,
                }
                .to_string();
                assert!(message.contains("12 of the 12 traits"), "{message}");
                assert!(message.contains("at the positions 0, 1,"), "{message}");
                assert!(message.contains("and 2 more"), "{message}");
                assert!(!message.contains("10, 11"), "{message}");
            }
            other => panic!("twelve fixed traits were standardized: {other:?}"),
        }
    }

    /// A buffer that does not hold the values of the table it was said to
    /// be. No argument of the Python or the TypeScript function reaches
    /// this: each binding crate takes the two numbers from the array it
    /// was given.
    #[test]
    fn a_buffer_that_does_not_hold_the_table_is_refused() {
        match pca(&ONE_TRAIT_FIXED, 4, 3, &CENTERED) {
            Err(Error::PcaTableOfAnotherSize {
                num_values,
                num_rows,
                num_cols,
            }) => {
                assert_eq!(num_values, 9);
                assert_eq!(num_rows, 4);
                assert_eq!(num_cols, 3);
            }
            other => panic!("a buffer of 9 values was read as 4 x 3: {other:?}"),
        }

        // One value too many is refused as well: the analysis of the first
        // nine would be of a table the caller did not mean.
        let mut longer = ONE_TRAIT_FIXED.to_vec();
        longer.push(0.0);
        match pca(&longer, 3, 3, &CENTERED) {
            Err(Error::PcaTableOfAnotherSize {
                num_values,
                num_rows,
                num_cols,
            }) => {
                assert_eq!(num_values, 10);
                assert_eq!(num_rows, 3);
                assert_eq!(num_cols, 3);
            }
            other => panic!("a buffer of 10 values was read as 3 x 3: {other:?}"),
        }
    }

    /// The table of 5 rows x 3 traits of "How it is verified", which has
    /// more rows than traits, so that the matrix that is decomposed is the
    /// 3 x 3 product of the traits, and which loses a component, its
    /// second trait having no variance. Iris, the other table with more
    /// rows than traits, keeps all four of its components, so nothing but
    /// this pins that the weights of that side are cut to the components
    /// that have variance.
    #[test]
    fn a_table_of_more_rows_than_traits_drops_the_component_with_no_variance() {
        let result = pca(&FIVE_ROWS, 5, 3, &CENTERED).expect("the analysis of the table");
        let projections = [
            -0.765373820993,
            -2.30958933885,
            -2.58720936972,
            -1.01308818829,
            -1.44494156963,
            -0.179287089168,
            5.62553292806,
            -0.270886092928,
            -0.828008167714,
            3.77285070924,
        ];
        let percent = [66.8261599339, 33.1738400661];
        let princomps = [
            0.15423335048,
            0.0,
            0.988034449602,
            0.988034449602,
            0.0,
            -0.15423335048,
        ];
        assert_the_result_is(
            &result,
            5,
            3,
            2,
            (&projections, &percent, &princomps),
            "the table of five rows",
        );
    }

    /// A table whose eigenvalues are near the largest `f64`. The
    /// percentages are what the table gives unscaled, and the count of the
    /// components does not change with the scale: both are computed in the
    /// order that does not overflow, the share of the total before the 100
    /// and the tolerance of the threshold before the largest eigenvalue.
    /// At 1e153 the largest eigenvalue is 1.4e307, so 100 times it is an
    /// infinity; at 2.5e153 it is 8.9e307, so it times the larger side of
    /// the table is an infinity as well.
    #[test]
    fn a_table_whose_eigenvalues_are_near_the_largest_float_gives_its_percentages() {
        for scale in [1e153, 2.5e153] {
            let table: Vec<f64> = NEAR_THE_LARGEST.iter().map(|value| value * scale).collect();
            let result = pca(&table, 3, 3, &CENTERED).expect("the analysis of the large table");
            assert_eq!(result.num_comps, 2, "the components at the scale {scale}");
            assert_close(
                &result.explained_variance_percent,
                &[78.8675134595, 21.1324865405],
                TOLERANCE,
                &format!("the percentages at the scale {scale}"),
            );
        }
    }

    /// A trait whose values sum above the largest `f64` has a mean that is
    /// not finite, and every centered value of it would be a NaN. It is
    /// found whenever the table is centered, standardized or not.
    #[test]
    fn a_trait_whose_mean_is_not_finite_is_refused() {
        let table = [f64::MAX, 1.0, f64::MAX, 2.0];
        match pca(&table, 2, 2, &CENTERED) {
            Err(Error::PcaTraitOutOfRange { position, problem }) => {
                assert_eq!(position, 0);
                assert_eq!(problem, TraitScale::MeanNotFinite);
                let message = Error::PcaTraitOutOfRange { position, problem }.to_string();
                assert!(message.contains("the trait at the position 0"), "{message}");
                assert!(message.contains("mean"), "{message}");
            }
            other => panic!("a trait whose mean is an infinity was centered: {other:?}"),
        }
    }

    /// A trait whose squared deviations sum above the largest `f64` has a
    /// standard deviation that is not finite, and dividing by it would
    /// make the trait a column of zeros, which looks like a trait with no
    /// variance and would leave the analysis with a weight of 0 and no
    /// word.
    #[test]
    fn a_trait_whose_standard_deviation_is_not_finite_is_refused() {
        let table = [1e154, 1.0, -1e154, 2.0];
        match pca(&table, 2, 2, &STANDARDIZED) {
            Err(Error::PcaTraitOutOfRange { position, problem }) => {
                assert_eq!(position, 0);
                assert_eq!(problem, TraitScale::DeviationNotFinite);
            }
            other => panic!("a trait whose deviation is an infinity was standardized: {other:?}"),
        }
    }

    /// A trait whose squared deviations all fall below the smallest `f64`
    /// above 0 has a standard deviation of 0 although its values differ,
    /// and dividing by it would give infinities, which the linear algebra
    /// would then refuse as a defect of popnei.
    #[test]
    fn a_trait_whose_standard_deviation_falls_to_zero_is_refused() {
        let table = [1e-200, 1.0, 2e-200, 2.0];
        match pca(&table, 2, 2, &STANDARDIZED) {
            Err(Error::PcaTraitOutOfRange { position, problem }) => {
                assert_eq!(position, 0);
                assert_eq!(problem, TraitScale::DeviationOfZero);
                let message = Error::PcaTraitOutOfRange { position, problem }.to_string();
                assert!(message.contains("not all equal"), "{message}");
            }
            other => panic!("a trait whose deviation is 0 was standardized: {other:?}"),
        }
    }

    /// A table with no direction to give: every trait the same value in
    /// every row, which centering turns into zeros, and a table of zeros
    /// that is not centered. pyNei gives 0 for every projection and a
    /// percentage of NaN for every component.
    #[test]
    fn a_table_in_which_no_trait_has_variance_is_refused() {
        let fixed = [1.0, 5.0, 1.0, 5.0, 1.0, 5.0];
        let result = pca(&fixed, 3, 2, &CENTERED);
        assert!(
            matches!(result, Err(Error::PcaNoTraitWithVariance)),
            "{result:?}"
        );
        if let Err(error) = result {
            let message = error.to_string();
            assert!(message.contains("no trait has variance"), "{message}");
        }

        let zeros = [0.0; 6];
        let result = pca(&zeros, 3, 2, &AS_IT_IS);
        assert!(
            matches!(result, Err(Error::PcaNoTraitWithVariance)),
            "{result:?}"
        );
    }

    /// An error of the linalg crate becomes one of popnei, with the
    /// operation that was being done. The values here are finite and their
    /// products are not, so the crate refuses the matrix it is asked to
    /// decompose.
    #[test]
    fn an_error_of_the_linear_algebra_is_wrapped_with_the_operation() {
        let table = [1e200, 2e200, 3e200, 4e200];
        match pca(&table, 2, 2, &AS_IT_IS) {
            Err(Error::PcaLinalg { operation, source }) => {
                let message = Error::PcaLinalg { operation, source }.to_string();
                assert!(message.contains(operation), "{message}");
                assert!(message.contains("finite"), "{message}");
            }
            other => panic!("a product that is not finite was taken: {other:?}"),
        }
    }

    /// The path of one of the files of `tests/reference/pca/`.
    fn the_reference_path(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/pca")
            .join(name)
    }

    /// The bytes of one of the VCFs of `tests/reference/pca/`, which the
    /// tests of the variants build their reader over.
    fn the_reference_vcf(name: &str) -> Vec<u8> {
        let path = the_reference_path(name);
        std::fs::read(&path)
            .unwrap_or_else(|error| panic!("{path}: {error}", path = path.display()))
    }

    /// A reader over the bytes of a VCF, in blocks of `num_vars_per_block`
    /// variants or of the size popnei chooses.
    fn reader_over(vcf: &[u8], num_vars_per_block: Option<usize>) -> VcfReader<Cursor<Vec<u8>>> {
        let options = VcfOptions {
            num_vars_per_block,
            ..VcfOptions::default()
        };
        match VcfReader::new(Cursor::new(vcf.to_vec()), options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        }
    }

    /// The VCF of the genotypes given, one line for each variant, with the
    /// individuals named `i0`, `i1` and so on.
    ///
    /// Every line lists the alleles `A` and `C,G,T`, whichever of them its
    /// genotypes hold: the alleles of a variant are the ones its genotypes
    /// hold and not the ones the VCF lists, which `docs/specs/pca.md` says
    /// and these tests lean on.
    fn vcf_of(num_individuals: usize, rows: &[Vec<String>]) -> Vec<u8> {
        let mut vcf = String::from(
            "##fileformat=VCFv4.2\n##contig=<ID=1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
        );
        for individual in 0..num_individuals {
            vcf.push_str(&format!("\ti{individual}"));
        }
        vcf.push('\n');
        for (var, genotypes) in rows.iter().enumerate() {
            let pos = var.checked_add(1).expect("the position of the variant");
            vcf.push_str(&format!("1\t{pos}\tvar{var}\tA\tC,G,T\t.\t.\t.\tGT"));
            for genotype in genotypes {
                vcf.push('\t');
                vcf.push_str(genotype);
            }
            vcf.push('\n');
        }
        vcf.into_bytes()
    }

    /// The same genotype for every individual of a variant.
    fn every_individual(genotype: &str, num_individuals: usize) -> Vec<String> {
        vec![genotype.to_string(); num_individuals]
    }

    /// The type of the reader of a second pass, which no test here gives:
    /// they all ask for no weights, and the type of a reader that is not
    /// there still has to be named.
    type NoSecondPass = VcfReader<Cursor<Vec<u8>>>;

    /// The principal components of the variants of one reader, with no
    /// second pass and so with no weights.
    fn the_pca_of_the_variants(
        reader: &mut impl BlockReader,
        options: &VariantPcaOptions,
    ) -> Result<Pca> {
        pca_of_variants::<_, NoSecondPass>(reader, None, options)
    }

    /// The options of the first pass: no weights, and a variant of more
    /// than two alleles is an error.
    const NO_WEIGHTS: VariantPcaOptions = VariantPcaOptions {
        transform_to_biallelic: false,
        num_prin_comps: 0,
    };

    /// The same with every allele that is not the major one counting the
    /// same, which is what a variant of more than two alleles takes.
    const BIALLELIC: VariantPcaOptions = VariantPcaOptions {
        transform_to_biallelic: true,
        num_prin_comps: 0,
    };

    /// The worked example of "How it is verified" of "The PCA of the
    /// variants" of `docs/specs/pca.md`, `worked.vcf`: 5 individuals and 5
    /// variants, of which the one with one allele and the one where every
    /// individual is heterozygous have no variance and are left out. The
    /// numbers are R's, from `worked.r.*.tsv`, with the projections
    /// multiplied by sqrt(n / (n - 1)).
    #[test]
    fn the_worked_example_of_the_variants_gives_the_numbers_of_r() {
        let mut reader = reader_over(&the_reference_vcf("worked.vcf"), None);
        let result = the_pca_of_the_variants(&mut reader, &NO_WEIGHTS).expect("the analysis");
        assert_eq!(result.num_rows, 5, "the individuals");
        assert_eq!(result.num_cols, 5, "the variants the pass gave");
        assert_eq!(result.num_comps, 3, "the components with variance");
        assert_eq!(
            result.used_cols,
            vec![0, 1, 4],
            "the variants that were used"
        );
        let projections = [
            -0.762287762540,
            0.994550847306,
            -0.320852228237,
            -0.863187556071,
            -0.875030288320,
            -0.007193635382,
            3.025289193432,
            -0.097439191351,
            -0.021646302880,
            -0.536626318751,
            0.852948920685,
            0.356885801880,
            -0.863187556071,
            -0.875030288320,
            -0.007193635382,
        ];
        assert_close(
            &result.projections,
            &projections,
            TOLERANCE,
            "the projections of the worked example",
        );
        assert_close(
            &result.explained_variance_percent,
            &[76.7440710447, 21.7166910409, 1.53923791438],
            TOLERANCE,
            "the percentages of the worked example",
        );
        // With `num_prin_comps` 0 there is no second pass and no weight,
        // and the variants that were used are there all the same.
        assert_eq!(result.num_prin_comps, 0, "the components of the weights");
        assert!(result.princomps.is_empty(), "{:?}", result.princomps);
    }

    /// `worked3.vcf` is the worked example with a sixth variant of three
    /// alleles, which is the error when every allele that is not the major
    /// one does not count the same. The message names the variant by its
    /// position among those given and says which argument to pass.
    #[test]
    fn a_variant_of_three_alleles_is_the_error_that_names_its_position() {
        let mut reader = reader_over(&the_reference_vcf("worked3.vcf"), None);
        match the_pca_of_the_variants(&mut reader, &NO_WEIGHTS) {
            Err(Error::PcaVariantWithMoreThanTwoAlleles {
                position,
                num_alleles,
            }) => {
                assert_eq!(position, 5);
                assert_eq!(num_alleles, 3);
                let message = Error::PcaVariantWithMoreThanTwoAlleles {
                    position,
                    num_alleles,
                }
                .to_string();
                assert!(message.contains("at the position 5"), "{message}");
                assert!(message.contains("transform_to_biallelic"), "{message}");
            }
            other => panic!("a variant of three alleles was read: {other:?}"),
        }
    }

    /// The same file with every allele that is not the major one counting
    /// the same, which gives 4 variants with variance and 3 components:
    /// pyNei and R give a fourth, whose percentage is 1e-31 and whose
    /// weights are noise. The numbers are the first three components of
    /// `worked3.r.projections.tsv` and `worked3.r.percent.tsv`.
    #[test]
    fn a_variant_of_three_alleles_counts_the_alleles_that_are_not_the_major_one_the_same() {
        let mut reader = reader_over(&the_reference_vcf("worked3.vcf"), None);
        let result = the_pca_of_the_variants(&mut reader, &BIALLELIC).expect("the analysis");
        assert_eq!(result.num_rows, 5, "the individuals");
        assert_eq!(result.num_cols, 6, "the variants the pass gave");
        assert_eq!(result.num_comps, 3, "the components with variance");
        assert_eq!(
            result.used_cols,
            vec![0, 1, 4, 5],
            "the variants that were used"
        );
        let projections: Vec<f64> = the_reference("worked3.r.projections.tsv", true)
            .iter()
            .flat_map(|row| row.iter().take(3).copied())
            .collect();
        assert_close(
            &result.projections,
            &projections,
            TOLERANCE,
            "the projections of worked3",
        );
        let percent: Vec<f64> = the_reference("worked3.r.percent.tsv", false)
            .into_iter()
            .take(3)
            .flatten()
            .collect();
        assert_close(
            &result.explained_variance_percent,
            &percent,
            TOLERANCE,
            "the percentages of worked3",
        );
    }

    /// The panel of "How it is verified", `sim_missing.vcf`: 200
    /// individuals, 1200 variants of two alleles, three subpopulations and
    /// 7128 genotypes missing whole. Centering takes one dimension out, so
    /// it has 199 components and not 200. The numbers are R's, from the
    /// table of the spec.
    #[test]
    fn the_panel_gives_the_numbers_of_r() {
        let mut reader = reader_over(&the_reference_vcf("sim_missing.vcf"), None);
        let result = the_pca_of_the_variants(&mut reader, &NO_WEIGHTS).expect("the analysis");
        assert_eq!(result.num_rows, 200, "the individuals");
        assert_eq!(result.num_cols, 1200, "the variants the pass gave");
        assert_eq!(
            result.used_cols.len(),
            1200,
            "every variant of the panel has variance"
        );
        assert_eq!(result.num_comps, 199, "the components with variance");
        // The first ten components of every individual, which is what
        // `sim_missing.r.projections.tsv` holds, and the two literals of
        // the table of the spec are the first three of the first two rows.
        let of_the_first_ten: Vec<f64> = result
            .projections
            .chunks_exact(result.num_comps)
            .flat_map(|of_an_individual| of_an_individual.iter().take(10).copied())
            .collect();
        assert_close(
            &of_the_first_ten,
            &row_after_row(the_reference("sim_missing.r.projections.tsv", true)),
            TOLERANCE,
            "the projections of the panel",
        );
        let of_s000: Vec<f64> = result.projections.iter().take(3).copied().collect();
        assert_close(
            &of_s000,
            &[1.57303591801, 12.9003037259, -5.07968498438],
            TOLERANCE,
            "the projections of s000",
        );
        let of_s001: Vec<f64> = result
            .projections
            .iter()
            .skip(result.num_comps)
            .take(3)
            .copied()
            .collect();
        assert_close(
            &of_s001,
            &[2.95102060115, 13.3478743206, -6.89340955267],
            TOLERANCE,
            "the projections of s001",
        );
        let percent: Vec<f64> = result
            .explained_variance_percent
            .iter()
            .take(3)
            .copied()
            .collect();
        assert_close(
            &percent,
            &[7.60577910944, 5.55551802152, 1.56605373725],
            TOLERANCE,
            "the percentages of the panel",
        );
    }

    /// The variants with no variance are not among the ones that were
    /// used, which
    /// `test_pca_from_variants_drops_the_monomorphic_variants` of pyNei
    /// asserts: the ones with one allele, the one where every individual
    /// is heterozygous, whose major allele frequency is 0.5 and which no
    /// filter by frequency would catch, and the one with no called
    /// genotype, which pyNei makes an error and popnei leaves out.
    #[test]
    fn the_variants_with_no_variance_are_not_among_the_ones_that_were_used() {
        let mut rows: Vec<Vec<String>> = (0..20_usize)
            .map(|var| {
                (0..8_usize)
                    .map(|individual| {
                        let first = individual.wrapping_add(var) % 2;
                        let second = usize::from(
                            individual.wrapping_mul(2).wrapping_add(var.wrapping_mul(5)) % 3 == 0,
                        );
                        format!("{first}/{second}")
                    })
                    .collect()
            })
            .collect();
        // Two variants with the reference allele alone, one with the
        // alternative one alone, one where every individual is
        // heterozygous and one with no called genotype.
        rows[3] = every_individual("0/0", 8);
        rows[17] = every_individual("0/0", 8);
        rows[11] = every_individual("1/1", 8);
        rows[5] = every_individual("0/1", 8);
        rows[9] = every_individual("./.", 8);
        let vcf = vcf_of(8, &rows);
        let mut reader = reader_over(&vcf, None);
        let result = the_pca_of_the_variants(&mut reader, &NO_WEIGHTS).expect("the analysis");
        assert_eq!(result.num_cols, 20, "the variants the pass gave");
        assert_eq!(
            result.used_cols,
            vec![0, 1, 2, 4, 6, 7, 8, 10, 12, 13, 14, 15, 16, 18, 19],
            "the variants that were used"
        );
        assert_eq!(result.num_rows, 8, "the individuals");
    }

    /// Every variant fixed is the error of a dataset with no direction to
    /// give, which `test_pca_from_variants_with_every_variant_fixed` of
    /// pyNei asserts. One individual gives the same error, since every
    /// variant of one individual has one dosage.
    #[test]
    fn a_dataset_in_which_no_variant_has_variance_is_refused() {
        let vcf = vcf_of(6, &vec![every_individual("0/0", 6); 10]);
        let mut reader = reader_over(&vcf, None);
        let result = the_pca_of_the_variants(&mut reader, &NO_WEIGHTS);
        match result {
            Err(Error::PcaNoVariantWithVariance) => {
                let message = Error::PcaNoVariantWithVariance.to_string();
                assert!(message.contains("nothing to do a PCA with"), "{message}");
            }
            other => panic!("a dataset of fixed variants was analysed: {other:?}"),
        }

        let of_one_individual = vec![
            vec!["0/1".to_string()],
            vec!["0/0".to_string()],
            vec!["1/1".to_string()],
        ];
        let vcf = vcf_of(1, &of_one_individual);
        let mut reader = reader_over(&vcf, None);
        let result = the_pca_of_the_variants(&mut reader, &NO_WEIGHTS);
        assert!(
            matches!(result, Err(Error::PcaNoVariantWithVariance)),
            "{result:?}"
        );
    }

    /// A reader that gives no variant is the other dataset there is
    /// nothing to analyse in: the steps of the variants let none through,
    /// or the source has none.
    #[test]
    fn a_reader_that_gives_no_variant_is_refused() {
        let vcf = vcf_of(3, &[]);
        let mut reader = reader_over(&vcf, None);
        let result = the_pca_of_the_variants(&mut reader, &NO_WEIGHTS);
        match result {
            Err(Error::PcaNoVariants) => {
                let message = Error::PcaNoVariants.to_string();
                assert!(
                    message.contains("no variants to do a PCA with"),
                    "{message}"
                );
            }
            other => panic!("a reader of no variant was analysed: {other:?}"),
        }
    }

    /// An individual with half of its genotypes missing falls on the side
    /// of its own population, which `test_pca_vars_with_missing_gts` of
    /// pyNei asserts: a missing genotype takes the mean of the dosages of
    /// its variant, and counting it as homozygous for the allele that is
    /// not the major one would pull the individual towards the other
    /// population. Two populations of 12 and 8 individuals are fixed for
    /// different alleles, and the first individual has no genotype at half
    /// the variants.
    #[test]
    fn an_individual_with_half_its_genotypes_missing_stays_with_its_population() {
        let of_the_first_pop = 12_usize;
        let num_individuals = 20_usize;
        let rows: Vec<Vec<String>> = (0..40_usize)
            .map(|var| {
                (0..num_individuals)
                    .map(|individual| {
                        if individual == 0 && var < 20 {
                            "./.".to_string()
                        } else if individual < of_the_first_pop {
                            "0/0".to_string()
                        } else {
                            "1/1".to_string()
                        }
                    })
                    .collect()
            })
            .collect();
        let vcf = vcf_of(num_individuals, &rows);
        let mut reader = reader_over(&vcf, None);
        let result = the_pca_of_the_variants(&mut reader, &NO_WEIGHTS).expect("the analysis");
        let first: Vec<f64> = result
            .projections
            .iter()
            .step_by(result.num_comps)
            .copied()
            .collect();
        assert_eq!(
            first.len(),
            num_individuals,
            "one projection of the first component per individual"
        );
        let with_the_missing = first.first().copied().expect("the first individual");
        let of_its_pop: f64 = first
            .get(1..of_the_first_pop)
            .expect("the rest of the first population")
            .iter()
            .sum();
        let of_the_other: f64 = first
            .get(of_the_first_pop..)
            .expect("the second population")
            .iter()
            .sum();
        assert!(
            with_the_missing.signum() == of_its_pop.signum(),
            "the individual with the missing genotypes is at {with_the_missing} and its population at {of_its_pop}"
        );
        assert!(
            with_the_missing.signum() != of_the_other.signum(),
            "the individual with the missing genotypes is at {with_the_missing} and the other population at {of_the_other}"
        );
    }

    /// The dosages of `test_mat012` and of
    /// `test_mat012_keeps_the_missing_gts` of pyNei, as the standardized
    /// values they become: the dosage of a genotype is how many of its
    /// alleles are not the major allele of the variant, a genotype with an
    /// allele missing has none and takes the mean, which is 0 once the
    /// variant is centered, and a variant whose called genotypes all have
    /// one dosage is left out.
    ///
    /// The values are each dosage less the mean of the called dosages,
    /// divided by the standard deviation that has all the individuals in
    /// it. They were worked out with python 3.13 on 22 September 2026 from
    /// the dosages pyNei asserts.
    #[test]
    #[expect(
        clippy::approx_constant,
        reason = "the standardized value of the dosage 0 of that variant, written with the 12 significant digits of the others, is the first digits of the inverse of sqrt(2)"
    )]
    fn the_dosages_are_pyneis_and_a_genotype_with_an_allele_missing_takes_the_mean() {
        // `test_mat012`: four individuals, and the second variant holds
        // the alleles 2 and 3, which pyNei refuses because it counts the
        // alleles of a whole chunk and popnei takes because it counts
        // those of each variant.
        for (what, gts, expected) in [
            (
                "the dosages 1 2 0 0",
                vec![0_i8, 1, 1, 1, 0, 0, 0, 0],
                vec![
                    0.301511344578,
                    1.50755672289,
                    -0.904534033733,
                    -0.904534033733,
                ],
            ),
            (
                "the dosages 1 0 0 2 of the alleles 2 and 3",
                vec![2_i8, 3, 2, 2, 2, 2, 3, 3],
                vec![
                    0.301511344578,
                    -0.904534033733,
                    -0.904534033733,
                    1.50755672289,
                ],
            ),
        ] {
            let mut scratch = RowScratch::of(4);
            let mut row = vec![0.0; 4];
            let used = the_standardized_row(&gts, 2, 0, &NO_WEIGHTS, &mut scratch, &mut row)
                .expect("the standardizing of the row");
            assert!(used, "{what}: the variant has variance");
            assert_close(&row, &expected, TOLERANCE, what);
        }

        // `test_mat012_keeps_the_missing_gts`: the half called genotype
        // `0/.` is missing, as the genotype with both of its alleles
        // missing is.
        let gts = vec![0_i8, 0, 0, 0, 0, 0, -1, -1, 1, 1, 0, -1];
        let mut scratch = RowScratch::of(6);
        let mut row = vec![0.0; 6];
        let used = the_standardized_row(&gts, 2, 0, &NO_WEIGHTS, &mut scratch, &mut row)
            .expect("the standardizing of the row");
        assert!(used, "the variant has variance");
        assert_close(
            &row,
            &[
                -0.707106781187,
                -0.707106781187,
                -0.707106781187,
                0.0,
                2.12132034356,
                0.0,
            ],
            TOLERANCE,
            "the dosages 0 0 0 missing 2 missing",
        );

        // A variant with one allele, whose called genotypes all have the
        // dosage 0, is left out and its row is not written.
        let gts = vec![0_i8; 8];
        let mut scratch = RowScratch::of(4);
        let mut row = vec![7.0; 4];
        let used = the_standardized_row(&gts, 2, 0, &NO_WEIGHTS, &mut scratch, &mut row)
            .expect("the standardizing of the row");
        assert!(!used, "a variant with one allele has no variance");
        assert_close(
            &row,
            &[7.0; 4],
            0.0,
            "the row of a variant that is left out",
        );
    }

    /// Three individuals whose projections are the same absolute value
    /// decide the sign of their component, and the first of them is the
    /// one that comes out positive, whatever the last bits of the
    /// eigendecomposition are.
    ///
    /// Six individuals of three variants, where the individuals 2 and 5
    /// have the same genotype at every variant and the individual 0 has
    /// the opposite dosage at every one, the mean of each variant being 1:
    /// the three projections of the first component are 1.98902864125 in
    /// exact arithmetic, the first of them negative. Accelerate's LAPACK
    /// gives that of the individual 5 one unit in the last place above the
    /// other two, so without the tolerance of [`TIED_PROJECTIONS`] the
    /// largest is positive already and the individual 0 comes out at
    /// -1.989, while numpy's eigendecomposition makes that of the
    /// individual 0 the largest and turns the component round. The numbers
    /// are numpy 2.5.3's, by the route of "What both analyses compute" of
    /// `docs/specs/pca.md`, on 22 September 2026.
    #[test]
    fn the_first_of_the_projections_of_the_largest_size_is_the_positive_one() {
        let rows: Vec<Vec<String>> = [
            ["0/0", "0/0", "1/1", "0/0", "1/1", "1/1"],
            ["0/0", "0/0", "1/1", "0/1", "0/1", "1/1"],
            ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1"],
        ]
        .iter()
        .map(|row| row.iter().map(ToString::to_string).collect())
        .collect();
        let vcf = vcf_of(6, &rows);
        let mut reader = reader_over(&vcf, None);
        let result = the_pca_of_the_variants(&mut reader, &NO_WEIGHTS).expect("the analysis");
        assert_eq!(result.num_comps, 3, "the components with variance");
        assert_close(
            &result.projections,
            &[
                1.98902864125,
                0.0,
                0.209201014086,
                1.28843624895,
                0.866025403784,
                -0.29988669924,
                -1.98902864125,
                0.0,
                -0.209201014086,
                1.28843624895,
                -0.866025403784,
                -0.29988669924,
                -0.58784385666,
                0.0,
                0.808974412566,
                -1.98902864125,
                0.0,
                -0.209201014086,
            ],
            TOLERANCE,
            "the projections of the tied individuals",
        );
        assert_close(
            &result.explained_variance_percent,
            &[86.3022285676, 8.33333333333, 5.36443809907],
            TOLERANCE,
            "the percentages",
        );
    }

    /// The codes of a genotype of 1, 2, 3, 4 and 5 alleles, with the major
    /// allele 0. The ploidies 1 to 4 take the arms of
    /// `the_codes_of_the_genotypes` that have the length of a genotype
    /// written into them and the ploidy of 5 takes the loop that reads it
    /// from the dataset, so this pins every arm of that match. A dosage is
    /// how many alleles of the genotype are not the major one, and a
    /// genotype with one allele missing has no dosage.
    #[test]
    fn the_codes_of_a_genotype_are_its_dosage_at_every_ploidy() {
        let of_one = |ploidy: usize| match NonZeroUsize::new(ploidy) {
            Some(of_a_genotype) => of_a_genotype,
            None => panic!("a ploidy of 0"),
        };
        let cases: [(usize, Vec<i8>, Vec<u8>); 5] = [
            (1, vec![0, 1, -1, 2], vec![0, 1, MISSING_CODE, 1]),
            (
                2,
                vec![0, 0, 0, 1, 1, 1, -1, 0, 0, -1, 1, -1],
                vec![0, 1, 2, MISSING_CODE, MISSING_CODE, MISSING_CODE],
            ),
            (
                3,
                vec![0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, -1],
                vec![0, 2, 3, MISSING_CODE],
            ),
            (
                4,
                vec![0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 2, 1, 1, 1, 1, -1, 0, 0, 0],
                vec![0, 2, 3, 4, MISSING_CODE],
            ),
            (
                5,
                vec![0, 0, 1, 1, 1, 0, 0, 0, 0, -1, 2, 2, 2, 2, 2],
                vec![3, MISSING_CODE, 5],
            ),
        ];
        for (ploidy, gts, expected) in cases {
            let mut codes = vec![0_u8; expected.len()];
            the_codes_of_the_genotypes(&gts, of_one(ploidy), 0, &mut codes);
            assert_eq!(codes, expected, "the codes of a ploidy of {ploidy}");
        }
    }

    /// The arms of `the_codes_of_the_genotypes` that have the length of a
    /// genotype written into them give the same codes as the loop that
    /// reads that length from the dataset, over a row of 131 individuals,
    /// which is not a whole number of the 16 genotypes a vector load of
    /// the first arm reads, so the tail of the loop is walked too. The
    /// alleles repeat with a period of 11, which no ploidy here divides,
    /// so a missing allele falls in every position of a genotype, and the
    /// run of four major alleles the period opens with gives a genotype of
    /// the dosage 0 at each of the four ploidies, which the two assertions
    /// at the end check.
    #[test]
    fn the_codes_of_a_fixed_ploidy_are_the_codes_of_the_loop_that_reads_the_ploidy() {
        const NUM_INDIVIDUALS: usize = 131;
        for ploidy in 1..=4_usize {
            let of_a_genotype = match NonZeroUsize::new(ploidy) {
                Some(of_a_genotype) => of_a_genotype,
                None => panic!("a ploidy of 0"),
            };
            let gts: Vec<i8> = (0..NUM_INDIVIDUALS * ploidy)
                .map(|allele| match allele % 11 {
                    0..=3 | 7 | 8 | 10 => 1,
                    4 | 9 => 0,
                    5 => 2,
                    _ => MISSING_ALLELE,
                })
                .collect();
            let mut of_the_match = vec![0_u8; NUM_INDIVIDUALS];
            let mut of_the_loop = vec![0_u8; NUM_INDIVIDUALS];
            the_codes_of_the_genotypes(&gts, of_a_genotype, 1, &mut of_the_match);
            the_codes_of_any_ploidy(&gts, ploidy, 1, &mut of_the_loop);
            assert_eq!(
                of_the_match, of_the_loop,
                "the codes of a ploidy of {ploidy}"
            );
            assert!(
                of_the_match.contains(&MISSING_CODE),
                "a genotype with an allele missing at a ploidy of {ploidy}"
            );
            assert!(
                of_the_match.contains(&0),
                "a genotype of the major allele only at a ploidy of {ploidy}"
            );
        }
    }

    /// The dosages of a tetraploid variant, which is the only fixture of
    /// a ploidy other than 2 in this module: a dosage is how many alleles
    /// of the genotype are not the major one at any ploidy, so code that
    /// read two alleles of each genotype, or that counted three dosages,
    /// would pass every diploid test here.
    ///
    /// Four individuals of the genotypes of the variant 0 of the
    /// tetraploid table of "How it is verified" of `docs/specs/pca.md`:
    /// the allele 1 was called 9 times and the allele 0 seven, so the
    /// major allele is 1 and the dosages are 4, 2, 1 and 0. The mean is
    /// 1.75 over the four called genotypes and the deviation sqrt(2.1875).
    #[test]
    fn the_dosages_of_a_tetraploid_variant_count_every_allele_of_the_genotype() {
        let gts = [0_i8, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1];
        let mut scratch = RowScratch::of(4);
        let mut row = vec![0.0; 4];
        let used = the_standardized_row(&gts, 4, 0, &NO_WEIGHTS, &mut scratch, &mut row)
            .expect("the standardizing of the row");
        assert!(used, "the variant has variance");
        assert_close(
            &row,
            &[
                1.52127765851,
                0.169030850946,
                -0.507092552837,
                -1.18321595662,
            ],
            TOLERANCE,
            "the standardized dosages 4 2 1 0",
        );
    }

    /// The whole analysis of the tetraploid dataset of "How it is
    /// verified", 4 individuals and 4 variants, whose third variant has
    /// one dosage among its called genotypes and is left out and whose
    /// fourth has a genotype with every allele missing. The numbers are
    /// numpy 2.5.3's, by the route of "What both analyses compute".
    #[test]
    fn a_tetraploid_dataset_gives_the_numbers_of_numpy() {
        let rows: Vec<Vec<String>> = [
            ["0/0/0/0", "0/0/1/1", "0/1/1/1", "1/1/1/1"],
            ["0/0/0/0", "0/0/0/0", "0/0/0/1", "0/0/1/1"],
            ["0/0/1/1", "0/0/1/1", "0/0/1/1", "0/0/1/1"],
            ["0/0/0/1", "./././.", "0/1/1/1", "1/1/1/1"],
        ]
        .iter()
        .map(|row| row.iter().map(ToString::to_string).collect())
        .collect();
        let vcf = vcf_of(4, &rows);
        let tetraploid = VcfOptions {
            ploidy: 4,
            ..VcfOptions::default()
        };
        let build = |options| match VcfReader::new(Cursor::new(vcf.clone()), options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        };
        let mut first_pass = build(tetraploid);
        let mut second_pass = build(tetraploid);
        let result = the_pca_of_two_passes(&mut first_pass, &mut second_pass, &with_weights(3))
            .expect("the analysis");
        assert_eq!(result.num_rows, 4, "the individuals");
        assert_eq!(result.num_cols, 4, "the variants the pass gave");
        assert_eq!(result.num_comps, 3, "the components with variance");
        assert_eq!(
            result.used_cols,
            vec![0, 1, 3],
            "the variants that were used"
        );
        assert_close(
            &result.projections,
            &[
                2.30350257852,
                -0.45498666453,
                -0.0168202039114,
                0.603260051101,
                0.692756790873,
                -0.0540239409669,
                -0.647570806465,
                0.0576947878489,
                0.143573693133,
                -2.25919182316,
                -0.295464914192,
                -0.0727295482544,
            ],
            TOLERANCE,
            "the projections of the tetraploid dataset",
        );
        assert_close(
            &result.explained_variance_percent,
            &[93.2778538477, 6.47960866887, 0.242537483385],
            TOLERANCE,
            "the percentages of the tetraploid dataset",
        );
        assert_close(
            &result.princomps,
            &[
                0.590326503368,
                -0.556614390531,
                0.584546866962,
                -0.327593824194,
                -0.827089115324,
                -0.45673392874,
                -0.737697028442,
                -0.0781281995539,
                0.670596062218,
            ],
            TOLERANCE,
            "the weights of the tetraploid dataset",
        );
    }

    /// A source of no individual has nobody to place on the components,
    /// and the standardizing would read the rows of its blocks in chunks
    /// of no allele. No reader of popnei gives one, so the test builds it.
    #[test]
    fn a_source_of_no_individual_is_refused() {
        let mut reader = OneBlock::of_no_individual();
        let result = the_pca_of_the_variants(&mut reader, &NO_WEIGHTS);
        assert!(matches!(result, Err(Error::PcaNoIndividual)), "{result:?}");
        let message = Error::PcaNoIndividual.to_string();
        assert!(message.contains("no individual"), "{message}");
    }

    /// A reader over one block that a test builds, for the sources no
    /// reader of popnei gives.
    struct OneBlock {
        individuals: Vec<String>,
        chroms: ChromTable,
        block: Option<Block>,
        needs: Needs,
    }

    impl OneBlock {
        /// A source of no individual, whose one block holds one variant of
        /// nobody.
        fn of_no_individual() -> OneBlock {
            OneBlock {
                individuals: Vec::new(),
                chroms: ChromTable::new(),
                block: Some(Block {
                    num_vars: 1,
                    num_individuals: 0,
                    ploidy: 2,
                    gts: Vec::new(),
                    chrom: None,
                    pos: None,
                    id: None,
                    alleles: None,
                    qual: None,
                }),
                needs: Needs::ALL,
            }
        }
    }

    impl BlockReader for OneBlock {
        fn next_block(&mut self) -> Result<Option<Block>> {
            Ok(self.block.take())
        }

        fn individuals(&self) -> &[String] {
            &self.individuals
        }

        fn ploidy(&self) -> usize {
            2
        }

        fn chroms(&self) -> &ChromTable {
            &self.chroms
        }

        fn set_needs(&mut self, needs: Needs) {
            self.needs = needs;
        }

        fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
            Vec::new()
        }
    }

    /// A run of genotypes that all have one dosage fills the counter of a
    /// byte: a run holds 255 genotypes, which is the longest run a byte
    /// counts without wrapping, and 255 genotypes of one dosage is what
    /// counts to 255.
    ///
    /// 260 individuals, the first 255 of the genotype `1/1` and the last
    /// five of `0/0`: the allele 1 was called 510 times and the allele 0
    /// ten, so the major allele is 1 and the dosages are 0 for the 255 and
    /// 2 for the five. The mean is 0.0384615384615 and the deviation
    /// 0.274670324175, which are python 3.13's of 22 September 2026; a
    /// counter that wrapped at the 256th genotype of a dosage would give
    /// another mean.
    #[test]
    fn a_run_of_one_dosage_fills_the_counter_of_a_byte() {
        let mut gts: Vec<i8> = Vec::new();
        for _ in 0..255 {
            gts.extend_from_slice(&[1, 1]);
        }
        for _ in 0..5 {
            gts.extend_from_slice(&[0, 0]);
        }
        let num_individuals = 260;
        let mut scratch = RowScratch::of(num_individuals);
        let mut row = vec![0.0; num_individuals];
        let used = the_standardized_row(&gts, 2, 0, &NO_WEIGHTS, &mut scratch, &mut row)
            .expect("the standardizing of the row");
        assert!(used, "the variant has variance");
        let of_each_dosage = [
            row.first().copied().expect("the first individual"),
            row.get(255).copied().expect("the first of the dosage 2"),
        ];
        assert_close(
            &of_each_dosage,
            &[-0.140028008403, 7.14142842854],
            TOLERANCE,
            "the standardized dosages of a full run",
        );
    }

    /// The dosages of a variant of more than 255 individuals are counted
    /// whole: the counts are taken in runs of 255 genotypes with counters
    /// of one byte, which is the largest run a byte counts without
    /// wrapping, and what each run counted is added to the counts of the
    /// variant.
    ///
    /// 300 individuals, of which 120 have the dosage 0, 90 the dosage 1
    /// and 90 the dosage 2, which the allele 0 is the major one of: 330
    /// of its alleles were called against 270 of the allele 1. The mean
    /// of the dosages is 0.9 and the deviation sqrt(0.69). A count that
    /// dropped a run, or that wrapped, would give another mean.
    #[test]
    fn the_dosages_of_more_than_255_individuals_are_counted_whole() {
        let mut gts: Vec<i8> = Vec::new();
        for (genotype, individuals) in [([0_i8, 0], 120), ([0, 1], 90), ([1, 1], 90)] {
            for _ in 0..individuals {
                gts.extend_from_slice(&genotype);
            }
        }
        let num_individuals = 300;
        let mut scratch = RowScratch::of(num_individuals);
        let mut row = vec![0.0; num_individuals];
        let used = the_standardized_row(&gts, 2, 0, &NO_WEIGHTS, &mut scratch, &mut row)
            .expect("the standardizing of the row");
        assert!(used, "the variant has variance");
        let of_each_dosage = [
            row.first().copied().expect("the first individual"),
            row.get(120).copied().expect("the first of the dosage 1"),
            row.get(210).copied().expect("the first of the dosage 2"),
        ];
        assert_close(
            &of_each_dosage,
            &[-1.08347267777, 0.120385853086, 1.32424438394],
            TOLERANCE,
            "the standardized dosages 0, 1 and 2",
        );
    }

    /// The genotypes of the five variants of the worked example as a
    /// block holds them, variant after variant and inside a variant
    /// individual after individual: a variant of two alleles, one with a
    /// genotype whose alleles are both missing, one with one allele, one
    /// where every individual is heterozygous, and one of the alleles 0
    /// and 2 with a half called genotype.
    ///
    /// Only the test of the two ways of reading the rows of a block uses
    /// it, and that test is of the targets that have threads.
    #[cfg(not(target_family = "wasm"))]
    const THE_WORKED_GENOTYPES: [i8; 50] = [
        0, 0, 0, 1, 1, 1, 0, 0, 0, 1, //
        1, 1, 1, 1, 0, 1, -1, -1, 1, 1, //
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //
        0, 1, 0, 1, 0, 1, 0, 1, 0, 1, //
        0, 2, 0, 0, 2, 2, 0, -1, 0, 0,
    ];

    /// The rows of a block read on the threads of rayon and read one after
    /// another, which is what WebAssembly does, give the same values and
    /// leave the same variants out: no row reads another and each one
    /// writes its own values.
    ///
    /// It is a test of the targets that have threads. In WebAssembly the
    /// two are one function, and a test there would compare it with
    /// itself.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_rows_read_on_threads_are_the_rows_read_one_after_another() {
        // Twenty copies of the five variants, so that the threads have
        // more than one row each.
        let gts = THE_WORKED_GENOTYPES.repeat(20);
        let num_individuals = 5;
        let mut on_threads = vec![0.0; gts.len() / 2];
        let used = the_standardized_rows(
            &gts,
            10,
            num_individuals,
            2,
            &NO_WEIGHTS,
            0,
            &mut on_threads,
        )
        .expect("the rows read on the threads of rayon");
        let mut one_by_one = vec![0.0; gts.len() / 2];
        let used_one_by_one = the_standardized_rows_one_by_one(
            &gts,
            10,
            num_individuals,
            2,
            &NO_WEIGHTS,
            0,
            &mut one_by_one,
        )
        .expect("the rows read one after another");
        assert_eq!(used, used_one_by_one, "the variants that were used");
        assert_eq!(
            used.iter().filter(|was_used| **was_used).count(),
            60,
            "three of the five variants of each copy have variance"
        );
        // The rows that were used hold the same values to the bit: the two
        // ways read each row the same.
        for (var, was_used) in used.iter().enumerate() {
            if !was_used {
                continue;
            }
            let row = the_row_of(var, num_individuals);
            assert_close(
                on_threads.get(row.clone()).expect("the row on the threads"),
                one_by_one.get(row).expect("the row read one by one"),
                0.0,
                &format!("the row of the variant {var}"),
            );
        }
    }

    /// The first pass over the variants of a VCF, in blocks of that many
    /// variants, for the five individuals of the worked example.
    fn the_first_pass_over(vcf: &[u8], num_vars_per_block: usize) -> FirstPass {
        let mut reader = reader_over(vcf, None);
        reader.set_needs(Needs::GTS);
        let mut blocks = match Reblock::new(&mut reader, Some(num_vars_per_block)) {
            Ok(blocks) => blocks,
            Err(error) => panic!("the blocks were not put back to one size: {error}"),
        };
        match the_first_pass(&mut blocks, &NO_WEIGHTS, 5, 2) {
            Ok(pass) => pass,
            Err(error) => panic!("the first pass: {error}"),
        }
    }

    /// The individuals x individuals matrix is added up block by block, so
    /// blocks of another size add the same numbers up in another order and
    /// the last bits of its entries could differ.
    ///
    /// The sizes are tried here and not at `pca_of_variants`, which puts
    /// `reblock` before the reader with the size popnei chooses: that is
    /// 10000 variants for five individuals, so the worked example is one
    /// block of it however the reader cut its own. Blocks of 3 variants
    /// are the ones that leave a variant with no variance before a
    /// variant with one, which is what the buffer of a block has to close
    /// up before its product.
    #[test]
    fn the_matrix_of_the_products_does_not_change_with_the_size_of_the_blocks() {
        let vcf = the_reference_vcf("worked.vcf");
        let whole = the_first_pass_over(&vcf, 5);
        assert_eq!(whole.num_cols, 5, "the variants the pass gave");
        assert_eq!(
            whole.used_cols,
            vec![0, 1, 4],
            "the variants that were used"
        );
        for num_vars_per_block in [1_usize, 2, 3] {
            let pass = the_first_pass_over(&vcf, num_vars_per_block);
            let what = format!("blocks of {num_vars_per_block} variants");
            assert_eq!(pass.num_cols, whole.num_cols, "{what}: the variants given");
            assert_eq!(pass.used_cols, whole.used_cols, "{what}: the variants used");
            assert_close(
                &pass.gram,
                &whole.gram,
                1e-10,
                &format!("{what}: the matrix of the products"),
            );
        }
    }

    /// The result does not change with the size of the blocks the reader
    /// gives, which is what `test_mat012_with_threads_is_the_same_as_without`
    /// of pyNei asserts of the threads. `reblock` puts them back to the
    /// size the analysis works in, and the test above is the one that
    /// changes that size.
    #[test]
    fn the_result_does_not_change_with_the_size_of_the_blocks_of_the_reader() {
        let vcf = the_reference_vcf("worked.vcf");
        let mut of_one_block = reader_over(&vcf, Some(5));
        let whole = the_pca_of_the_variants(&mut of_one_block, &NO_WEIGHTS).expect("the analysis");
        for num_vars_per_block in [1_usize, 2] {
            let mut reader = reader_over(&vcf, Some(num_vars_per_block));
            let result = the_pca_of_the_variants(&mut reader, &NO_WEIGHTS).expect("the analysis");
            let what = format!("blocks of {num_vars_per_block} variants");
            assert_eq!(
                result.num_cols, whole.num_cols,
                "{what}: the variants given"
            );
            assert_eq!(
                result.used_cols, whole.used_cols,
                "{what}: the variants used"
            );
            assert_eq!(result.num_comps, whole.num_comps, "{what}: the components");
            assert_close(
                &result.projections,
                &whole.projections,
                1e-10,
                &format!("{what}: the projections"),
            );
            assert_close(
                &result.explained_variance_percent,
                &whole.explained_variance_percent,
                1e-10,
                &format!("{what}: the percentages"),
            );
        }
    }

    /// The weights need a second pass over the variants, since the weight
    /// of a variant is worked out from the eigenvectors, which are known
    /// when the first pass ends. A `num_prin_comps` above 0 with no reader
    /// for that pass is the error, and it comes before a variant is read.
    #[test]
    fn the_weights_without_a_second_pass_are_refused() {
        let mut reader = reader_over(&the_reference_vcf("worked.vcf"), None);
        let options = VariantPcaOptions {
            transform_to_biallelic: false,
            num_prin_comps: 3,
        };
        match the_pca_of_the_variants(&mut reader, &options) {
            Err(Error::PcaSecondPassMissing { num_prin_comps }) => {
                assert_eq!(num_prin_comps, 3);
                let message = Error::PcaSecondPassMissing { num_prin_comps }.to_string();
                assert!(message.contains("second pass"), "{message}");
                assert!(message.contains("3 components"), "{message}");
            }
            other => panic!("the weights were given with no second pass: {other:?}"),
        }
    }

    /// The principal components of the variants of two readers over the
    /// same variants, which is how the weights are got.
    fn the_pca_of_two_passes(
        first_pass: &mut impl BlockReader,
        second_pass: &mut impl BlockReader,
        options: &VariantPcaOptions,
    ) -> Result<Pca> {
        pca_of_variants(first_pass, Some(second_pass), options)
    }

    /// The same options with the weights of that many components.
    fn with_weights(num_prin_comps: usize) -> VariantPcaOptions {
        VariantPcaOptions {
            transform_to_biallelic: false,
            num_prin_comps,
        }
    }

    /// The genotypes of the five variants of the worked example, as the
    /// lines of a VCF write them.
    fn the_worked_rows() -> Vec<Vec<String>> {
        [
            ["0/0", "0/1", "1/1", "0/0", "0/1"],
            ["1/1", "1/1", "0/1", "./.", "1/1"],
            ["0/0", "0/0", "0/0", "0/0", "0/0"],
            ["0/1", "0/1", "0/1", "0/1", "0/1"],
            ["0/2", "0/0", "2/2", "0/.", "0/0"],
        ]
        .iter()
        .map(|row| row.iter().map(ToString::to_string).collect())
        .collect()
    }

    /// The weights of the worked example, the three components of its
    /// three variants that were used, which R gives in
    /// `worked.r.princomps.tsv` as one row for each component: that is the
    /// layout of `princomps`. With more components than there are, the
    /// weights are those of the components there are.
    #[test]
    fn the_worked_example_gives_the_weights_of_r() {
        let vcf = the_reference_vcf("worked.vcf");
        let princomps = row_after_row(the_reference("worked.r.princomps.tsv", true));
        for num_prin_comps in [3_usize, 10] {
            let mut first_pass = reader_over(&vcf, None);
            let mut second_pass = reader_over(&vcf, None);
            let result = the_pca_of_two_passes(
                &mut first_pass,
                &mut second_pass,
                &with_weights(num_prin_comps),
            )
            .expect("the analysis");
            let what = format!("the weights of {num_prin_comps} components");
            assert_eq!(result.num_comps, 3, "{what}: the components");
            assert_eq!(
                result.num_prin_comps, 3,
                "{what}: the components the weights are given for"
            );
            assert_eq!(result.used_cols, vec![0, 1, 4], "{what}: the variants used");
            assert_close(&result.princomps, &princomps, TOLERANCE, &what);
        }
    }

    /// The weights of `worked3.vcf` with every allele that is not the
    /// major one counting the same: the three components popnei gives, of
    /// the four `worked3.r.princomps.tsv` holds.
    #[test]
    fn the_weights_of_a_variant_of_three_alleles_are_the_ones_of_r() {
        let vcf = the_reference_vcf("worked3.vcf");
        let mut first_pass = reader_over(&vcf, None);
        let mut second_pass = reader_over(&vcf, None);
        let options = VariantPcaOptions {
            transform_to_biallelic: true,
            num_prin_comps: 3,
        };
        let result = the_pca_of_two_passes(&mut first_pass, &mut second_pass, &options)
            .expect("the analysis");
        assert_eq!(result.num_prin_comps, 3, "the components of the weights");
        let princomps: Vec<f64> = row_after_row(the_reference("worked3.r.princomps.tsv", true))
            .into_iter()
            .take(12)
            .collect();
        assert_close(
            &result.princomps,
            &princomps,
            TOLERANCE,
            "the weights of worked3",
        );
    }

    /// The panel with the weights of its first 10 components, which is the
    /// default of `do_pca_from_variants`: 10 components of 1200 variants,
    /// where the two sides of `princomps` differ in count, so a transpose
    /// that swapped them would not even have the shape. The numbers are
    /// R's, from `sim_missing.r.princomps.tsv`, and the five literals of
    /// the table of "How it is verified" are in it.
    #[test]
    fn the_panel_gives_the_weights_of_r() {
        let vcf = the_reference_vcf("sim_missing.vcf");
        let mut first_pass = reader_over(&vcf, None);
        let mut second_pass = reader_over(&vcf, None);
        let result = the_pca_of_two_passes(&mut first_pass, &mut second_pass, &with_weights(10))
            .expect("the analysis");
        assert_eq!(result.num_comps, 199, "the components with variance");
        assert_eq!(result.num_prin_comps, 10, "the components of the weights");
        assert_eq!(
            result.princomps.len(),
            12000,
            "ten components of the 1200 variants"
        );
        let princomps = row_after_row(the_reference("sim_missing.r.princomps.tsv", true));
        assert_close(
            &result.princomps,
            &princomps,
            TOLERANCE,
            "the weights of the panel",
        );
        // The four literals of the table of the spec: the weights of the
        // variants 0 and 1 in the first two components.
        let of_the_first_two = [
            result.princomps.first().copied().expect("variant 0, PC000"),
            result.princomps.get(1).copied().expect("variant 1, PC000"),
            result
                .princomps
                .get(1200)
                .copied()
                .expect("variant 0, PC001"),
            result
                .princomps
                .get(1201)
                .copied()
                .expect("variant 1, PC001"),
        ];
        assert_close(
            &of_the_first_two,
            &[
                0.0665938024571,
                0.0131069975393,
                -0.0296301101543,
                0.0518944546112,
            ],
            TOLERANCE,
            "the weights of the variants 0 and 1",
        );
        // The projections are what the second pass must not change.
        let of_s000: Vec<f64> = result.projections.iter().take(3).copied().collect();
        assert_close(
            &of_s000,
            &[1.57303591801, 12.9003037259, -5.07968498438],
            TOLERANCE,
            "the projections of s000",
        );
    }

    /// A second pass over other variants than the first is refused: its
    /// weights would belong to variants the eigenvectors did not come
    /// from. It is what a source that changed between the two passes
    /// gives, and the message says what differed.
    #[test]
    fn a_second_pass_over_other_variants_is_refused() {
        let rows = the_worked_rows();
        let vcf = vcf_of(5, &rows);

        // A variant more, which has variance: the first pass used three
        // variants and the second has a fourth.
        let mut with_one_more = rows.clone();
        with_one_more.push(every_individual("0/1", 5));
        with_one_more[5] = ["1/1", "0/0", "0/1", "0/0", "1/1"]
            .iter()
            .map(ToString::to_string)
            .collect();
        // A variant fewer, which leaves the count short.
        let with_one_fewer = rows.get(..4).expect("the first four variants").to_vec();
        // A variant that had variance and has none.
        let mut without_the_first = rows.clone();
        without_the_first[0] = every_individual("0/0", 5);

        for (what, changed, expected) in [
            (
                "a variant more",
                vcf_of(5, &with_one_more),
                VariantsOfTheSecondPass::UsedNow(5),
            ),
            (
                "a variant fewer",
                vcf_of(5, &with_one_fewer),
                VariantsOfTheSecondPass::Count {
                    found: 4,
                    expected: 5,
                },
            ),
            (
                "a variant that lost its variance",
                vcf_of(5, &without_the_first),
                VariantsOfTheSecondPass::NotUsedNow(0),
            ),
            (
                "another dataset",
                vcf_of(4, &vcf_rows_of_four_individuals()),
                VariantsOfTheSecondPass::Dataset {
                    num_individuals: 4,
                    ploidy: 2,
                },
            ),
        ] {
            let mut first_pass = reader_over(&vcf, None);
            let mut second_pass = reader_over(&changed, None);
            match the_pca_of_two_passes(&mut first_pass, &mut second_pass, &with_weights(2)) {
                Err(Error::PcaSecondPassDiffers { problem }) => {
                    assert_eq!(problem, expected, "{what}");
                    let message = Error::PcaSecondPassDiffers { problem }.to_string();
                    assert!(message.contains("second pass"), "{what}: {message}");
                    assert!(message.contains("source changed"), "{what}: {message}");
                }
                other => panic!("{what} was taken: {other:?}"),
            }
        }
    }

    /// The five variants of the worked example for four individuals, which
    /// is a second pass over another dataset.
    fn vcf_rows_of_four_individuals() -> Vec<Vec<String>> {
        the_worked_rows()
            .into_iter()
            .map(|row| row.get(..4).expect("four individuals").to_vec())
            .collect()
    }

    /// The weights of the worked example from a second pass over blocks of
    /// that many variants.
    fn the_weights_over_blocks_of(vcf: &[u8], num_vars_per_block: usize) -> Vec<f64> {
        let first = the_first_pass_over(vcf, 5);
        let eigen = match eigh_lower(first.gram.clone(), 5) {
            Ok(eigen) => eigen,
            Err(error) => panic!("the eigendecomposition: {error}"),
        };
        let num_comps = the_components_with_variance(&eigen.values, 5, first.used_cols.len());
        let scaled_vectors = the_scaled_vectors_of(&eigen, 5, num_comps);
        let after = AfterTheFirstPass {
            scaled_vectors: &scaled_vectors,
            num_prin_comps: num_comps,
            used_cols: &first.used_cols,
            num_cols: first.num_cols,
        };
        let mut reader = reader_over(vcf, None);
        reader.set_needs(Needs::GTS);
        let mut blocks = match Reblock::new(&mut reader, Some(num_vars_per_block)) {
            Ok(blocks) => blocks,
            Err(error) => panic!("the blocks were not put back to one size: {error}"),
        };
        match the_weights_of_a_second_pass(&mut blocks, &NO_WEIGHTS, 5, 2, &after) {
            Ok(weights) => weights,
            Err(error) => panic!("the second pass: {error}"),
        }
    }

    /// The weights do not change with the size of the blocks of the second
    /// pass: each block is multiplied on its own and its weights are
    /// written into the columns of its variants, so blocks of another size
    /// cut the same product in another place.
    ///
    /// The sizes are tried where the pass takes the blocks it is given, as
    /// they are for the matrix of the products: `pca_of_variants` puts
    /// `reblock` before each reader with the size popnei chooses, which is
    /// 10000 variants for five individuals. Blocks of 3 variants are the
    /// ones whose second block holds one used variant, which goes into the
    /// third column of the weights.
    #[test]
    fn the_weights_do_not_change_with_the_size_of_the_blocks() {
        let vcf = the_reference_vcf("worked.vcf");
        let whole = the_weights_over_blocks_of(&vcf, 5);
        assert_eq!(whole.len(), 9, "three components of three variants");
        for num_vars_per_block in [1_usize, 2, 3] {
            let weights = the_weights_over_blocks_of(&vcf, num_vars_per_block);
            assert_close(
                &weights,
                &whole,
                1e-10,
                &format!("the weights over blocks of {num_vars_per_block} variants"),
            );
        }
    }

    /// The result does not change with the size of the blocks the second
    /// reader gives, which `reblock` puts back to the size the pass works
    /// in, as it does for the first.
    #[test]
    fn the_weights_do_not_change_with_the_size_of_the_blocks_of_the_reader() {
        let vcf = the_reference_vcf("worked.vcf");
        let princomps = row_after_row(the_reference("worked.r.princomps.tsv", true));
        for num_vars_per_block in [1_usize, 2, 5] {
            let mut first_pass = reader_over(&vcf, Some(num_vars_per_block));
            let mut second_pass = reader_over(&vcf, Some(num_vars_per_block));
            let result = the_pca_of_two_passes(&mut first_pass, &mut second_pass, &with_weights(3))
                .expect("the analysis");
            assert_close(
                &result.princomps,
                &princomps,
                1e-10,
                &format!("the weights over blocks of {num_vars_per_block} variants"),
            );
        }
    }

    /// A dataset whose ploidy or whose individuals are beyond what the
    /// analysis counts in is refused before a variant is read, and the
    /// message says which of the two it is. The third size, more variants
    /// than a `usize` counts, no test reaches.
    #[test]
    fn a_dataset_beyond_what_the_analysis_counts_in_is_refused() {
        let vcf = vcf_of(2, &[]);
        let options = VcfOptions {
            ploidy: 255,
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::new(Cursor::new(vcf), options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        };
        match the_pca_of_the_variants(&mut reader, &NO_WEIGHTS) {
            Err(Error::PcaVariantsTooLarge { problem }) => {
                assert_eq!(problem, VariantsTooLarge::Ploidy(255));
                let message = Error::PcaVariantsTooLarge { problem }.to_string();
                assert!(message.contains("255 alleles each"), "{message}");
                assert!(
                    message.contains(&MAX_PLOIDY_OF_THE_VARIANTS.to_string()),
                    "{message}"
                );
            }
            other => panic!("a ploidy of 255 was analysed: {other:?}"),
        }

        let too_many = MAX_INDIVIDUALS_OF_THE_VARIANTS
            .checked_add(1)
            .expect("one individual more than the analysis takes");
        let vcf = vcf_of(too_many, &[]);
        let mut reader = reader_over(&vcf, None);
        match the_pca_of_the_variants(&mut reader, &NO_WEIGHTS) {
            Err(Error::PcaVariantsTooLarge { problem }) => {
                assert_eq!(problem, VariantsTooLarge::Individuals(too_many));
                let message = Error::PcaVariantsTooLarge { problem }.to_string();
                assert!(message.contains("2147483647"), "{message}");
            }
            other => panic!("46341 individuals were analysed: {other:?}"),
        }
    }
}
