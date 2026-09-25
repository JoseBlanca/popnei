//! The kinship of every pair of the individuals of a dataset: how much more
//! of their genome the two share than two individuals drawn at random from
//! the same panel would.
//!
//! It is the matrix of VanRaden 2008, which plink2's `--make-rel` computes
//! and which GCTA is built on. An entry off the diagonal is twice the
//! coancestry of its pair, about 0.5 for full sibs or for a parent and a
//! child and near 0 for two individuals with no recent ancestor in common;
//! an entry on the diagonal is 1 plus the inbreeding of that individual.
//! Entries below 0 are ordinary and mean a pair less alike than the average
//! pair of the panel, because the whole matrix is measured against that
//! average.
//!
//! Each variant becomes one number per individual, its dosage, how many
//! alleles of the genotype are not the major allele of the variant, and the
//! dosages are centered and divided by the standard deviation the allele
//! frequency of the variant gives it under Hardy Weinberg, which is the
//! **standardized dosage** `z` of `docs/specs/kinship.md`. The entry of the
//! pair `i`, `j` is the sum over the variants of `z[v, i] * z[v, j]`
//! divided by how many of those variants have a called genotype in both of
//! them, the **per pair denominator**. `docs/specs/pca.md` divides by
//! another number, the standard deviation of the dosages themselves, and
//! the two agree only when the genotypes are in Hardy Weinberg proportions;
//! the pass over a row that both calculations make, of
//! [`crate::variant`], takes the divisor from its caller.
//!
//! [`principal_components`] places each individual along the directions in
//! which the panel varies most, from the eigenvectors of that matrix, which
//! is what a user gives an association study as covariates when they
//! account for the structure of the panel without a mixed model.

use std::fmt;
use std::mem;

use popnei_linalg::{add_self_product_lower, eigh_lower};

use crate::block::{Block, BlockReader, Reblock, with_one_block_ahead};
use crate::error::{Error, Result};
use crate::pca::{fix_the_sign_of, the_components_with_variance, the_projections_of};
use crate::phases::{Phase, timed};
use crate::variant::{
    DosageOptions, DosageScale, MISSING_ALLELE, Needs, RowPositions, the_standardized_block,
};

/// The most individuals a kinship is taken of, which is the most any
/// calculation of popnei builds a matrix of the individuals by the
/// individuals for, of [`crate::variant`]. This module checks it at its own
/// entry, before the first block is read.
pub use crate::variant::MAX_INDIVIDUALS_OF_THE_VARIANTS;

/// The kinship of every pair of a set of individuals, which the pass gives
/// away so that the Python binding hands the matrix to numpy without
/// copying it; the wasm binding copies it into a `Float64Array`.
#[derive(Debug, Clone)]
pub struct Kinship {
    /// How many individuals the matrix has on each of its two sides.
    pub num_individuals: usize,
    /// How many variants had variance among these individuals and were
    /// used. A variant whose called genotypes all have one dosage, and one
    /// with no called genotype, are in neither the sum of a pair nor its
    /// denominator.
    pub num_vars: u64,
    /// How many variants the reader gave, used or not, which is the
    /// `num_vars` of the pass stats.
    ///
    /// It is counted after `reblock`, which is where this pass sees the
    /// variants; `reblock` gives every variant it is given, so it is the
    /// count a reader between the pass and its source would make.
    pub num_vars_given: u64,
    /// `num_individuals` x `num_individuals`, row after row, symmetric.
    pub matrix: Vec<f64>,
}

/// Which size of a dataset is beyond what a kinship is taken on.
///
/// Neither is a dataset of this world: the objectives of popnei reach 10000
/// individuals and a million variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KinshipTooLarge {
    /// The individuals of the dataset, which are more than
    /// [`MAX_INDIVIDUALS_OF_THE_VARIANTS`].
    Individuals(usize),
    /// The reader gave more variants than a `usize` counts, which is
    /// 4294967295 in WebAssembly, where a `usize` is 32 bits.
    ///
    /// No test reaches it and none can reach it natively: a reader would
    /// have to give 18446744073709551616 variants, one for every value of
    /// a `u64` and one more. It is here for the browser, where a dataset of
    /// 4295 million variants is a file of that many lines and not a size of
    /// this world either, and for the rule that a count popnei cannot hold
    /// is an error and never a number that wrapped.
    Variants,
}

impl fmt::Display for KinshipTooLarge {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::Individuals(num_individuals) => write!(
                formatter,
                "it has {num_individuals} individuals, and the individuals x individuals matrix of the kinship would hold more values than the 2147483647 the linear algebra counts in, which is {MAX_INDIVIDUALS_OF_THE_VARIANTS} individuals"
            ),
            Self::Variants => write!(
                formatter,
                "the variants given are more than {largest}, which is what this machine counts the variants of the kinship in",
                largest = usize::MAX
            ),
        }
    }
}

/// The kinship of the individuals of a reader, over the variants it gives.
///
/// `individuals` are the positions, among the individuals the reader has,
/// of the ones the matrix is of, in the order it has them, and `None` is
/// all of them in the reader's order. Every frequency, mean and denominator
/// is of those individuals: the kinship of some of them is not the rows and
/// columns of the kinship of the whole panel. `transform_to_biallelic` says
/// that a variant of more than two alleles among its called genotypes is
/// read with every allele that is not the major one counting the same,
/// which the principal components of the variants also take, and without it
/// such a variant is an error.
///
/// One pass over the blocks. Each variant that has variance is standardized
/// and the block is multiplied by itself into an individuals x individuals
/// accumulator; a block with a genotype missing is counted into a second
/// accumulator of the same shape, how many variants each pair had called in
/// both of its individuals. A variant with no variance is in neither, and
/// [`Kinship::num_vars`] counts those that were used.
///
/// The pass borrows the reader and does not take it, so that whoever built
/// the chain of filters reads their counts when it returns. It asks for the
/// genotypes alone and puts [`Reblock`] before it, since a filter leaves
/// blocks of uneven size and the product of a block is matrix work.
///
/// # Errors
///
/// [`Error::PassGaveNoVariant`] when the reader gives no variant, with what
/// each filter of the pass was given and kept, and
/// [`Error::KinshipNoVariantWithVariance`] when no variant of it has
/// variance among these individuals.
/// [`Error::KinshipPairWithNoVariantCalled`] when two of the individuals
/// have no variant called in both, whose entry would be divided by 0.
/// [`Error::VariantWithMoreThanTwoAlleles`] when a variant has more than
/// two different alleles among its called genotypes and
/// `transform_to_biallelic` is false.
/// [`Error::KinshipNoIndividual`] when there is no individual to take the
/// kinship of, and [`Error::KinshipVariantsTooLarge`] when the individuals
/// or the variants are more than the calculation counts in.
/// [`Error::VariantPloidyTooLarge`] when the ploidy is above
/// [`MAX_PLOIDY_OF_THE_VARIANTS`](crate::variant::MAX_PLOIDY_OF_THE_VARIANTS),
/// which the pass over a row raises. [`Error::FieldsNotInTheBlock`] when a
/// block holds no genotypes, what [`Block::retain_individuals`] refuses of
/// `individuals`, [`Error::KinshipLinalg`] when a product could not be
/// done, and whatever the reader fails with.
pub fn calc_kinship<R: BlockReader>(
    reader: &mut R,
    individuals: Option<&[usize]>,
    transform_to_biallelic: bool,
) -> Result<Kinship> {
    let ploidy = reader.ploidy();
    let num_individuals = match individuals {
        Some(individuals) => individuals.len(),
        None => reader.individuals().len(),
    };
    if num_individuals == 0 {
        // The pass over a block reads the rows in chunks of one genotype
        // for each individual, which would be chunks of no allele, and
        // there is no pair to give. No reader of popnei gives such a
        // source, and this function is public.
        return Err(Error::KinshipNoIndividual);
    }
    if num_individuals > MAX_INDIVIDUALS_OF_THE_VARIANTS {
        return Err(Error::KinshipVariantsTooLarge {
            problem: KinshipTooLarge::Individuals(num_individuals),
        });
    }
    // The genotypes are all this reads, so a reader over a file leaves the
    // columns of a variant unparsed.
    reader.set_needs(Needs::GTS);
    // A filter leaves blocks of uneven size and the product of a block is
    // matrix work, so the blocks are put back to one size first.
    let mut blocks = Reblock::new(reader, None)?;
    let options = DosageOptions {
        transform_to_biallelic,
        // The divisor of the kinship, sqrt(ploidy * p * (1 - p)), which is
        // what makes an entry twice a coancestry and what plink2 and GCTA
        // take. The principal components of the variants give the same
        // pass the standard deviation of the dosages instead.
        scale: DosageScale::OfHardyWeinberg,
    };
    let ThePass {
        mut gram,
        denominators,
        num_vars,
        num_vars_given,
    } = the_pass_over_the_blocks(&mut blocks, individuals, num_individuals, ploidy, &options)?;
    if num_vars_given == 0 {
        let filters = blocks.filtering_stats();
        return Err(Error::PassGaveNoVariant {
            // The filter nearest the source was given what the source
            // gave; with no filter the pass gave what the source gave,
            // which is nothing.
            num_vars_of_the_source: filters.last().map_or(0, |(_, stats)| stats.vars_processed),
            filters,
        });
    }
    if num_vars == 0 {
        return Err(Error::KinshipNoVariantWithVariance);
    }
    the_entries_of(&mut gram, num_individuals, &denominators, num_vars)?;
    Ok(Kinship {
        num_individuals,
        num_vars,
        num_vars_given,
        matrix: gram,
    })
}

/// Where each individual falls along the directions in which a panel
/// varies most, taken from the kinship of its individuals.
///
/// A user gives these to an association study as covariates, which is how
/// the structure of a panel is accounted for without a mixed model.
#[derive(Debug, Clone)]
pub struct KinshipPcs {
    /// How many components were given: the `num_pcs` that were asked for,
    /// or the components the kinship has above the tolerance when it has
    /// fewer of them.
    pub num_comps: usize,
    /// The individuals x [`Self::num_comps`] matrix of where each
    /// individual falls along each component, row after row, with the
    /// individuals in the order the kinship has them.
    pub projections: Vec<f64>,
}

/// The principal components of a kinship, `num_pcs` of them at most.
///
/// With `lambda_j` the eigenvalues of the matrix from the largest and
/// `u_j` its eigenvectors, the component `j` is `u_j * sqrt(lambda_j)`:
/// where each individual falls along the direction in which the panel
/// varies the `j`th most.
///
/// A component whose eigenvalue is not above the largest eigenvalue times
/// the individuals times the difference between 1 and the next number an
/// `f64` holds is not given, so a kinship with fewer components than were
/// asked for gives the ones it has and [`KinshipPcs::num_comps`] says how
/// many. The per pair denominators of a dataset with genotypes missing put
/// eigenvalues below 0, which are the length of nothing: pyNei gives
/// `num_pcs` components whatever the eigenvalue, taking the square root of
/// its absolute value. The sign of each component is fixed by the rule of
/// `docs/specs/pca.md`, the projection of the largest absolute value
/// positive, so that the two backends of the eigendecomposition and the
/// three builds of popnei give one answer.
///
/// These are close to the principal components of the variants the kinship
/// was calculated from and they are not the same: a kinship divides each
/// variant by `sqrt(ploidy * p * (1 - p))`, with `p` its allele frequency,
/// and [`crate::pca::pca_of_variants`] by the standard deviation of its
/// dosages, and the two agree only when the genotypes are in Hardy
/// Weinberg proportions.
///
/// The matrix is copied here, since the eigendecomposition writes the
/// eigenvectors over the matrix it is given and the caller keeps its
/// kinship: 800 MB at 10000 individuals.
/// [`principal_components_of`] takes a matrix by value and copies nothing,
/// which is what a caller that does not keep a [`Kinship`] calls.
///
/// # Errors
///
/// [`Error::KinshipNoIndividual`] when the kinship has no individual,
/// which leaves nobody to place along anything,
/// [`Error::KinshipValueNotFinite`] when a value of the matrix is an
/// infinity or a NaN, and [`Error::KinshipLinalg`] when the
/// eigendecomposition could not be done.
pub fn principal_components(kinship: &Kinship, num_pcs: usize) -> Result<KinshipPcs> {
    principal_components_of(kinship.matrix.clone(), kinship.num_individuals, num_pcs)
}

/// The same components, of a matrix this takes over.
///
/// `matrix` is `num_individuals` x `num_individuals`, row after row, and
/// its lower half is what is read: the eigendecomposition writes the
/// eigenvectors over it, so nothing is copied here and the memory of one
/// matrix is what the components cost besides the decomposition's own. A
/// caller that keeps its kinship calls [`principal_components`], which
/// copies.
///
/// Every value is checked to be finite before the decomposition, the upper
/// half among them, although only the lower half is read: a matrix a user
/// built and then wrote a NaN into is a wrong argument, and what the
/// linear algebra would say of it names a matrix `g` and a row of it.
///
/// # Errors
///
/// [`Error::KinshipNoIndividual`] when `num_individuals` is 0, which
/// leaves nobody to place along anything,
/// [`Error::KinshipValueNotFinite`] when a value of the matrix is an
/// infinity or a NaN, with where it is, and [`Error::KinshipLinalg`] when
/// the matrix does not hold one value for each pair of the individuals or
/// the eigendecomposition could not be done.
pub fn principal_components_of(
    matrix: Vec<f64>,
    num_individuals: usize,
    num_pcs: usize,
) -> Result<KinshipPcs> {
    if num_individuals == 0 {
        return Err(Error::KinshipNoIndividual);
    }
    the_values_are_finite(&matrix, num_individuals)?;
    // The eigendecomposition reads the lower half of the matrix and writes
    // the eigenvectors over it, so the matrix this was given is what it
    // works in.
    let eigen = eigh_lower(matrix, num_individuals).map_err(|source| Error::KinshipLinalg {
        operation: "eigendecomposition",
        source,
    })?;
    // The matrix is the individuals by the individuals, so the tolerance
    // of `docs/specs/pca.md`, the largest eigenvalue times the larger side
    // of the matrix times the difference between 1 and the next number an
    // `f64` holds, takes the individuals on both sides.
    let num_comps =
        the_components_with_variance(&eigen.values, num_individuals, num_individuals).min(num_pcs);
    let mut projections = the_projections_of(&eigen, num_individuals, num_comps);
    for component in 0..num_comps {
        // Whether the component was turned round is for a caller that
        // holds the weight of each trait in it, which the principal
        // components of a table do and a kinship does not.
        fix_the_sign_of(&mut projections, component, num_comps);
    }
    Ok(KinshipPcs {
        num_comps,
        projections,
    })
}

/// That every value of the matrix of a kinship is finite, with where the
/// first one that is not is.
///
/// The whole matrix is read and not the lower half alone: a user who wrote
/// a value into a frame after it was checked wrote it somewhere, and an
/// infinity or a NaN above the diagonal says the matrix is wrong as surely
/// as one below it.
///
/// # Errors
///
/// [`Error::KinshipValueNotFinite`] with the row and the column of the
/// value, counted from 0 among the individuals of the kinship.
fn the_values_are_finite(matrix: &[f64], num_individuals: usize) -> Result<()> {
    for (row, of_the_row) in matrix.chunks(num_individuals.max(1)).enumerate() {
        for (col, value) in of_the_row.iter().enumerate() {
            if !value.is_finite() {
                return Err(Error::KinshipValueNotFinite {
                    row,
                    col,
                    value: *value,
                });
            }
        }
    }
    Ok(())
}

/// How many variants each pair of individuals had called in both of them,
/// which the sum of the pair is divided by.
enum Denominators {
    /// No genotype of a variant that was used is missing, so every pair had
    /// the variants that were used, and the denominator is that one count.
    OfEveryPair,
    /// A genotype was missing: the lower half of the individuals x
    /// individuals matrix of how many of the variants that were used are
    /// called in both individuals of each pair, row after row.
    OfThePair(Vec<f64>),
}

/// What the pass over the blocks leaves, which is all that is kept from one
/// block to the next besides the buffers a block is read into.
struct ThePass {
    /// The lower half of the sum over the variants that were used of the
    /// standardized dosage of one individual times that of another,
    /// individuals x individuals, row after row.
    gram: Vec<f64>,
    /// How many variants each pair had called in both.
    denominators: Denominators,
    /// How many variants had variance and were used.
    num_vars: u64,
    /// How many variants the reader gave, used or not.
    num_vars_given: u64,
}

/// The one pass over the blocks: each block is standardized into a buffer
/// of the variants that were used x the individuals, and the product of
/// that buffer with itself is added to the sums of the pairs.
///
/// The genotypes that were called are counted into the denominators only
/// from the first block that has one missing: until then every pair has
/// every variant that was used, and that count is a number and not a
/// matrix. The matrix is then made with the variants of the blocks that
/// went before in **every** entry, since those blocks were called
/// everywhere, and not at zeros, which would lose them.
///
/// `individuals` are the positions of the individuals the kinship is of
/// among those of the blocks, and `None` is all of them in the order of the
/// reader; `num_individuals` is how many of them that is, 1 or more.
///
/// # Errors
///
/// What the standardizing of a block refuses, what
/// [`Block::retain_individuals`] refuses of `individuals`,
/// [`Error::KinshipVariantsTooLarge`] when the variants given are more than
/// a `usize` counts, [`Error::KinshipLinalg`] when a product could not be
/// done, and whatever the reader fails with.
fn the_pass_over_the_blocks<R: BlockReader>(
    blocks: &mut R,
    individuals: Option<&[usize]>,
    num_individuals: usize,
    ploidy: usize,
    options: &DosageOptions,
) -> Result<ThePass> {
    let mut gram = vec![0.0; the_values_of(num_individuals)];
    let mut denominators = Denominators::OfEveryPair;
    let mut num_vars = 0_u64;
    let mut num_vars_given = 0_u64;
    // The buffers of one block, kept from one block to the next so that a
    // pass over a million variants asks for them once. The rows of
    // `standardized` that were not used are left as they were and nothing
    // reads them.
    let mut standardized: Vec<f64> = Vec::new();
    let mut buffers = TheBuffersOfTheDenominators::default();
    // The blocks are read on a thread of its own, one block ahead, so that
    // the read of the next block and the product of the one in hand overlap.
    // Over 100000 variants of 1000 individuals of a vars file,
    // `docs/reports/perf-read-ahead-2026-09-25.md` measured the kinship at
    // 0.198 s without this and 0.184 s with it on 18 cores, and 0.400 s
    // against 0.302 s on one thread: what the pass waits on the handle is
    // 0.008 s of the 0.027 s the read costs at 18 cores, against 0.179 s of
    // work on the same blocks. In wasm, where there is no thread, the blocks
    // come one after another as they did. The chain of readers is lent and
    // not given away, which is what lets `calc_kinship` read its counts when
    // this returns.
    with_one_block_ahead(blocks, |blocks| {
        while let Some(mut block) = timed(Phase::NextBlock, || blocks.next_block())? {
            timed(Phase::Work, || -> Result<()> {
                if let Some(individuals) = individuals {
                    // The frequencies, the means and the denominators are of the
                    // individuals the kinship is of, so the others leave the block
                    // before anything is counted.
                    block.retain_individuals(individuals)?;
                }
                let used = the_standardized_block(
                    &block,
                    num_individuals,
                    ploidy,
                    options,
                    RowPositions {
                        // Where the first variant of this block is among those the
                        // reader has given, which the error of a variant with more
                        // than two alleles names. A pass of more variants than a
                        // `usize` counts is the error of a pass too large, and in
                        // WebAssembly, where a `usize` is 32 bits, it is reachable.
                        first: usize::try_from(num_vars_given)
                            .map_err(|_| the_variants_are_too_many())?,
                        too_many: the_variants_are_too_many,
                    },
                    &mut standardized,
                )?;
                let kept = used.iter().filter(|was_used| **was_used).count();
                add_self_product_lower(&standardized, kept, num_individuals, &mut gram).map_err(
                    |source| Error::KinshipLinalg {
                        operation: "product of a block of variants with itself",
                        source,
                    },
                )?;
                if kept > 0 {
                    the_denominators_of_the_block(
                        &block,
                        &used,
                        kept,
                        num_individuals,
                        num_vars,
                        &mut buffers,
                        &mut denominators,
                    )?;
                }
                num_vars = num_vars
                    .checked_add(the_count_of(kept))
                    .ok_or_else(the_variants_are_too_many)?;
                num_vars_given = num_vars_given
                    .checked_add(the_count_of(block.num_vars))
                    .ok_or_else(the_variants_are_too_many)?;
                Ok(())
            })?;
        }
        Ok(())
    })?;
    Ok(ThePass {
        gram,
        denominators,
        num_vars,
        num_vars_given,
    })
}

/// The buffers the denominators of a block are worked out in, kept from
/// one block to the next so that a pass over a million variants asks the
/// machine for them once.
#[derive(Default)]
struct TheBuffersOfTheDenominators {
    /// The genotypes that were called, one value for each individual of
    /// each variant that was used, which the route of the product
    /// multiplies by itself. It stays empty while no block takes that
    /// route, and it is the largest buffer of the kinship when one does:
    /// 40 MB at 5000 variants of 1000 individuals.
    called: Vec<f64>,
    /// The individuals one variant has a missing genotype in, in growing
    /// order, which the route of the counts fills again for each variant.
    of_the_variant: Vec<usize>,
    /// How many of the variants of the block that were used each
    /// individual has a missing genotype in, which the route of the
    /// counts fills again for each block.
    of_each_individual: Vec<u64>,
}

/// The variants of one block that each pair had called in both, added to
/// the denominators of the blocks before it.
///
/// `used` says which variants of the block have variance, `kept` how many
/// of them that is, 1 or more, and `num_vars_before` how many variants the
/// blocks before this one used. `buffers` holds what the two routes below
/// work in, which is kept from one block to the next.
///
/// A block with no missing genotype adds its variants to every pair
/// alike, so it is a number added to each entry and not a product, which
/// is what `_KinshipCalc.calc_for_chunk` of `pynei/gwas.py` skips too: the
/// product of the genotypes that were called is of the same shape as the
/// product of the standardized dosages and nearly doubles the work of the
/// linear algebra.
///
/// # Errors
///
/// [`Error::KinshipLinalg`] when the product of the genotypes that were
/// called could not be done, and what the sizes of the block refuse.
fn the_denominators_of_the_block(
    block: &Block,
    used: &[bool],
    kept: usize,
    num_individuals: usize,
    num_vars_before: u64,
    buffers: &mut TheBuffersOfTheDenominators,
    denominators: &mut Denominators,
) -> Result<()> {
    let alleles_per_var = block.alleles_per_var()?;
    let sum_of_the_squares = the_sum_of_the_squares_of_the_missing(block, used, alleles_per_var);
    if sum_of_the_squares == 0 {
        // Every pair had every variant of the block, so this block is a
        // number and not a product.
        if let Denominators::OfThePair(of_the_pairs) = denominators {
            let of_the_block = kept as f64;
            // The upper half is left as it was and nothing reads it.
            for count in of_the_pairs.iter_mut() {
                *count += of_the_block;
            }
        }
        return Ok(());
    }
    // The matrix is taken out of the denominators and put back below, so
    // that one call of each route serves both the block that makes it and
    // the blocks that add to it. The blocks before the first one with a
    // genotype missing had no genotype missing, so each of their variants
    // is in the denominator of every pair: the matrix it makes starts at
    // that count in every entry and not at 0, which would lose them.
    let mut of_the_pairs = match *denominators {
        Denominators::OfThePair(ref mut of_the_pairs) => mem::take(of_the_pairs),
        Denominators::OfEveryPair => {
            vec![num_vars_before as f64; the_values_of(num_individuals)]
        }
    };
    // The two routes add the same whole numbers to the same entries, so
    // they leave the matrix the same bit for bit; the test
    // `the_two_routes_to_the_denominators_agree_entry_for_entry` holds
    // them to that.
    let added =
        if the_counts_are_cheaper_than_the_product(sum_of_the_squares, kept, num_individuals) {
            the_counted_denominators_of(
                block,
                used,
                alleles_per_var,
                kept,
                num_individuals,
                buffers,
                &mut of_the_pairs,
            );
            Ok(())
        } else {
            the_called_genotypes_of(
                block,
                used,
                alleles_per_var,
                kept,
                num_individuals,
                &mut buffers.called,
            );
            add_the_called_genotypes(&buffers.called, kept, num_individuals, &mut of_the_pairs)
        };
    *denominators = Denominators::OfThePair(of_the_pairs);
    added
}

/// The sum over the variants of a block that were used of the square of
/// how many individuals each of them has a missing genotype in, which is 0
/// when every one of those genotypes was called.
///
/// It is what the choice between the two routes is made on, because it is
/// the work the route of the counts does: that route touches the entry of
/// every pair of individuals a variant is missing in both of, and those
/// pairs are that square, halved.
///
/// The sum saturates where a `u64` ends instead of raising. A block whose
/// missing genotypes make that many pairs is one the product is cheaper
/// for by a wide margin, and a saturated sum is what sends it there.
fn the_sum_of_the_squares_of_the_missing(
    block: &Block,
    used: &[bool],
    alleles_per_var: usize,
) -> u64 {
    // The rows of the block hold one genotype of its ploidy for each of
    // its individuals, and a ploidy of 0 does not reach here: the pass
    // that standardized the same rows refuses it.
    let ploidy = block.ploidy.max(1);
    let mut total = 0_u64;
    for genotypes in the_rows_used(block, used, alleles_per_var) {
        // A row with nothing missing is left at the whole row read at
        // once, which stops at the first allele that is missing and which
        // the compiler reads several bytes at a time. The genotypes of a
        // row are only cut apart where there is something to count, so a
        // block with every genotype called costs what it did when this
        // was the `any` of one row that had a missing allele.
        if !genotypes.contains(&MISSING_ALLELE) {
            continue;
        }
        let missing = the_count_of(
            genotypes
                .chunks_exact(ploidy)
                .filter(|genotype| genotype.contains(&MISSING_ALLELE))
                .count(),
        );
        total = total.saturating_add(missing.saturating_mul(missing));
    }
    total
}

/// How many of the increments of the route of the counts one entry of the
/// product's result pays for, which is where the two routes cross.
///
/// The product does the individuals squared times the variants that were
/// used multiply-adds, which took 5.1 ps each on the Accelerate of this
/// machine at every shape tried; the route of the counts does the sum of
/// the squares of the missing, halved, read-modify-writes scattered over
/// the matrix of the denominators, and one of those took 0.64 ns at 1000
/// individuals, 0.98 ns at 3000 and 1.33 ns at 6000, since the matrix they
/// walk is 8, 72 and 288 MB and falls out of the caches. So the two routes
/// crossed at a sum of the squares of 1 part in 139 of the variants used
/// times the individuals squared at 1000 individuals, 1 in 289 at 3000 and
/// 1 in 524 at 6000. Measured on 23 September 2026 on an Apple M5 Pro with
/// a trial that ran both routes over one block of 5 million genotypes at
/// twelve rates of missing genotypes, on one thread, and checked at each
/// of them that the two gave the same matrix.
///
/// 512 is taken from the largest of those three, so that the route of the
/// counts is not the slower one at up to 6000 individuals. What that costs
/// is the blocks between 1 part in 512 and 1 in 139 at 1000 individuals,
/// from about 4 to about 7 in 100 genotypes missing, which go to the
/// product where the counts would have been up to twice as fast.
const THE_INCREMENTS_AN_ENTRY_OF_THE_PRODUCT_PAYS_FOR: u64 = 512;

/// Whether counting the denominators of a block is cheaper than the
/// product of its genotypes that were called with itself.
///
/// `sum_of_the_squares` is what [`the_sum_of_the_squares_of_the_missing`]
/// gave for the block and `kept` how many of its variants were used. A
/// product that would not fit in a `u64` saturates, which sends the block
/// to the counts, and a block of that many values does not exist: the
/// individuals are at most 46340 and the variants of a block at most
/// 10000, whose product with the individuals again is 2.1e13.
fn the_counts_are_cheaper_than_the_product(
    sum_of_the_squares: u64,
    kept: usize,
    num_individuals: usize,
) -> bool {
    let individuals = the_count_of(num_individuals);
    let of_the_product = the_count_of(kept)
        .saturating_mul(individuals)
        .saturating_mul(individuals);
    sum_of_the_squares.saturating_mul(THE_INCREMENTS_AN_ENTRY_OF_THE_PRODUCT_PAYS_FOR)
        <= of_the_product
}

/// The variants of a block that each pair had called in both, counted into
/// the lower half of the denominators with no product.
///
/// With `kept` the variants of the block that were used and `M` the ones
/// an individual has a missing genotype in, the pair `i`, `j` gains
///
/// ```text
/// kept - |M_i| - |M_j| + |M_i and M_j|
/// ```
///
/// which is the variants of the block that neither of the two is missing.
/// Every term is a count of whole things and every one of them is far
/// below 2^53, so each entry grows by exactly the whole number the product
/// of the genotypes that were called would have added to it.
///
/// It is the route for a block with few missing genotypes. The third term
/// is what costs: it is one increment for each pair of individuals a
/// variant is missing in both of, which is quadratic in how many
/// individuals a variant is missing in, where the product does the
/// individuals squared for every variant whatever is missing.
fn the_counted_denominators_of(
    block: &Block,
    used: &[bool],
    alleles_per_var: usize,
    kept: usize,
    num_individuals: usize,
    buffers: &mut TheBuffersOfTheDenominators,
    of_the_pairs: &mut [f64],
) {
    // The rows of the block hold one genotype of its ploidy for each of
    // its individuals, and a ploidy of 0 does not reach here: the pass
    // that standardized the same rows refuses it.
    let ploidy = block.ploidy.max(1);
    buffers.of_each_individual.clear();
    buffers.of_each_individual.resize(num_individuals, 0);
    for genotypes in the_rows_used(block, used, alleles_per_var) {
        // A row with nothing missing adds nothing to any pair here: its
        // variant is in the count of every one of them, which the last
        // loop of this pass adds. The whole row is read at once for that,
        // as in the pass that chose this route.
        if !genotypes.contains(&MISSING_ALLELE) {
            continue;
        }
        buffers.of_the_variant.clear();
        for (at, genotype) in genotypes.chunks_exact(ploidy).enumerate() {
            if genotype.contains(&MISSING_ALLELE) {
                buffers.of_the_variant.push(at);
            }
        }
        for one in &buffers.of_the_variant {
            // The block holds one genotype for each individual of the
            // matrix, which the pass that standardized the same rows
            // checked, so every individual a variant is missing in has its
            // count here; a count is at most the variants of the block,
            // which is far from where a `u64` saturates.
            if let Some(count) = buffers.of_each_individual.get_mut(*one) {
                *count = count.saturating_add(1);
            }
        }
        // The individuals a variant is missing in come in growing order,
        // so the pair of the one being read with any of those before it,
        // itself included, is in the lower half of the matrix. The pair of
        // an individual with itself is wanted: the entry of `i` with `i`
        // is kept - |M_i|, which the three terms give.
        for (upto, one) in buffers.of_the_variant.iter().enumerate() {
            // The individual is below the individuals of the matrix, so
            // the row of its pairs starts inside it.
            let Some(row) = one
                .checked_mul(num_individuals)
                .and_then(|start| of_the_pairs.get_mut(start..))
            else {
                continue;
            };
            for other in buffers.of_the_variant.iter().take(upto.saturating_add(1)) {
                if let Some(value) = row.get_mut(*other) {
                    *value += 1.0;
                }
            }
        }
    }
    for ((one, row), of_one) in of_the_pairs
        .chunks_exact_mut(num_individuals.max(1))
        .enumerate()
        .zip(&buffers.of_each_individual)
    {
        // The counts are whole numbers at or below the variants of the
        // block, so each of them is exact in an `f64` and the entry grows
        // by a whole number. The sum of the three terms is at or above 0,
        // though this one term alone can be below it.
        let of_the_row = kept as f64 - *of_one as f64;
        for (value, of_other) in row
            .iter_mut()
            .zip(&buffers.of_each_individual)
            .take(one.saturating_add(1))
        {
            *value += of_the_row - *of_other as f64;
        }
    }
}

/// The product of the genotypes that were called with themselves, added to
/// the lower half of the denominators: the entry of the pair `i`, `j` grows
/// by how many variants of the block are called in both.
///
/// # Errors
///
/// [`Error::KinshipLinalg`] when the product could not be done.
fn add_the_called_genotypes(
    called: &[f64],
    kept: usize,
    num_individuals: usize,
    of_the_pairs: &mut [f64],
) -> Result<()> {
    add_self_product_lower(called, kept, num_individuals, of_the_pairs).map_err(|source| {
        Error::KinshipLinalg {
            operation: "product of the genotypes of a block that were called with themselves",
            source,
        }
    })
}

/// The genotypes of the variants of the block that were used, 1 for a
/// genotype every allele of which was called and 0 for one with an allele
/// missing, into `called`: one value for each individual of each of those
/// variants, in the order of the block.
///
/// A genotype with **any** allele missing has no dosage of its own and
/// takes the mean of its variant, so once the variant is centered it pulls
/// its pairs nowhere; this is what keeps it out of their denominators too.
/// `alleles_per_var` is the alleles of one variant of the block, its
/// individuals times its ploidy.
///
/// The rows are read one after another and not on the threads of rayon:
/// they are a pass over the bytes of the block with no arithmetic in it,
/// beside the pass that standardizes the same rows on the threads.
fn the_called_genotypes_of(
    block: &Block,
    used: &[bool],
    alleles_per_var: usize,
    kept: usize,
    num_individuals: usize,
    called: &mut Vec<f64>,
) {
    // The buffer belongs to the pass and is as long as the largest block it
    // has read, which asks the machine for nothing after the first block.
    called.resize(the_values_of_the_rows(kept, num_individuals), 0.0);
    // The rows of the block hold one genotype of its ploidy for each of
    // its individuals, and a ploidy of 0 does not reach here: the pass that
    // standardized the same rows refuses it.
    let ploidy = block.ploidy.max(1);
    for (row, genotypes) in called
        .chunks_exact_mut(num_individuals.max(1))
        .zip(the_rows_used(block, used, alleles_per_var))
    {
        for (value, genotype) in row.iter_mut().zip(genotypes.chunks_exact(ploidy)) {
            *value = if genotype.iter().all(|allele| *allele != MISSING_ALLELE) {
                1.0
            } else {
                0.0
            };
        }
    }
}

/// The genotypes of each variant of the block that was used, in the order
/// of the block.
///
/// `alleles_per_var` is the alleles of one variant, the individuals of the
/// block times its ploidy; a block of none has no row here, which is what a
/// block that holds no genotype gives.
fn the_rows_used<'a>(
    block: &'a Block,
    used: &'a [bool],
    alleles_per_var: usize,
) -> impl Iterator<Item = &'a [i8]> {
    block
        .gts
        .chunks(alleles_per_var.max(1))
        .zip(used)
        .filter(|(_, was_used)| **was_used)
        .map(|(row, _)| row)
}

/// The entries of the matrix: the sum of each pair divided by how many
/// variants that pair had called in both, and the lower half mirrored into
/// the upper, so that the matrix a user reads is whole.
///
/// `num_vars` is how many variants were used, 1 or more, which is the
/// denominator of every pair when no genotype was missing.
///
/// # Errors
///
/// [`Error::KinshipPairWithNoVariantCalled`] when two individuals have no
/// variant called in both, whose entry would be divided by 0.
fn the_entries_of(
    gram: &mut [f64],
    num_individuals: usize,
    denominators: &Denominators,
    num_vars: u64,
) -> Result<()> {
    match *denominators {
        Denominators::OfEveryPair => {
            let of_every_pair = num_vars as f64;
            for (row, upto) in the_rows_with_their_lower_half(gram, num_individuals) {
                for value in row.iter_mut().take(upto) {
                    *value /= of_every_pair;
                }
            }
        }
        Denominators::OfThePair(ref of_the_pairs) => {
            for ((row, upto), counts) in the_rows_with_their_lower_half(gram, num_individuals)
                .zip(of_the_pairs.chunks_exact(num_individuals))
            {
                // The row is the individual `upto - 1`, since the lower
                // half of a row holds the entries up to the diagonal.
                let other = upto.saturating_sub(1);
                for (one, (value, count)) in row.iter_mut().zip(counts).take(upto).enumerate() {
                    // The counts are sums of the 1 and the 0 of the
                    // genotypes that were called, whole numbers added
                    // exactly, so a pair with one variant is at 1.0.
                    if *count < 1.0 {
                        return Err(Error::KinshipPairWithNoVariantCalled {
                            one,
                            other,
                            num_vars_of_one: the_variants_called_in(
                                of_the_pairs,
                                one,
                                num_individuals,
                            ),
                            num_vars_of_other: the_variants_called_in(
                                of_the_pairs,
                                other,
                                num_individuals,
                            ),
                        });
                    }
                    *value /= *count;
                }
            }
        }
    }
    the_lower_half_mirrored(gram, num_individuals);
    Ok(())
}

/// Each row of the matrix with how many of its values are in the lower
/// half, which are the ones up to the diagonal: 1 for the first row, 2 for
/// the second.
fn the_rows_with_their_lower_half(
    matrix: &mut [f64],
    num_individuals: usize,
) -> impl Iterator<Item = (&mut [f64], usize)> {
    matrix
        .chunks_exact_mut(num_individuals.max(1))
        .zip(1_usize..)
}

/// The lower half of the matrix written into its upper half, so that the
/// entry of the pair `i`, `j` is at both of its places.
///
/// The product of the linear algebra writes the lower half alone and leaves
/// the upper as it was, which is the 0 the matrix was made with.
fn the_lower_half_mirrored(matrix: &mut [f64], num_individuals: usize) {
    // The rows that have been read, which is one for each entry of the
    // lower half of the row being read: the entry `i`, `j` of that row
    // belongs to the row `j` at the column `i`.
    let mut before: Vec<&mut [f64]> = Vec::with_capacity(num_individuals);
    for (at, row) in matrix
        .chunks_exact_mut(num_individuals.max(1))
        .enumerate()
        .take(num_individuals)
    {
        for (earlier, value) in before.iter_mut().zip(row.iter()) {
            // The rows read so far are fewer than the individuals, and
            // every row holds one value for each of them, so each of them
            // has a value at this column.
            if let Some(target) = earlier.get_mut(at) {
                *target = *value;
            }
        }
        before.push(row);
    }
}

/// How many of the variants that were used are called in one individual,
/// which is the entry of that individual with itself, for the error of a
/// pair with no variant called in both.
fn the_variants_called_in(of_the_pairs: &[f64], individual: usize, num_individuals: usize) -> u64 {
    let at = individual
        .checked_mul(num_individuals)
        .and_then(|row| row.checked_add(individual));
    let count = at
        .and_then(|at| of_the_pairs.get(at))
        .copied()
        .unwrap_or(0.0);
    the_whole_number_of(count)
}

/// The largest whole number an `f64` holds one by one, 2^53: a count above
/// it would have been added in steps of more than 1.
const LARGEST_WHOLE_NUMBER: f64 = 9007199254740992.0;

/// A count of variants that the product of the genotypes that were called
/// gave, as the whole number it is.
///
/// The product adds 1 for each variant a pair had called and 0 for each one
/// it did not, so the count is a whole number at or above 0. Anything else
/// is 0, which is what an entry with no room in the matrix reads as.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "the count is checked to be finite, at or above 0 and at most 2^53 before it is read as a whole number, and it is a sum of the 1 and the 0 of the genotypes that were called"
)]
fn the_whole_number_of(count: f64) -> u64 {
    if count.is_finite() && (0.0..=LARGEST_WHOLE_NUMBER).contains(&count) {
        count as u64
    } else {
        0
    }
}

/// How many values the individuals x individuals matrix of a kinship holds.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the individuals are at most MAX_INDIVIDUALS_OF_THE_VARIANTS, 46340, which calc_kinship checks at its entry, and their square is 2147395600, which a usize holds where it is 32 bits too"
)]
fn the_values_of(num_individuals: usize) -> usize {
    num_individuals * num_individuals
}

/// How many values the rows of the variants of one block that were used
/// hold, one for each individual of each of them.
///
/// The rows of a block are the rows the pass over it standardizes, which
/// were asked of the machine already, so the buffer of the genotypes that
/// were called is of a size that fits beside it.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the variants that were used times the individuals are at most the values of the block the pass over it standardized into a buffer of the same shape, which this machine gave"
)]
fn the_values_of_the_rows(num_vars: usize, num_individuals: usize) -> usize {
    num_vars * num_individuals
}

/// A count of the variants of one block as the count over the whole dataset
/// it is added to.
///
/// A `usize` is 64 bits natively and 32 in WebAssembly, and both of them fit
/// in a `u64`. A platform where one did not would give the largest `u64`
/// here, which the sum of the pass refuses as more variants than it counts
/// instead of wrapping.
fn the_count_of(num_vars: usize) -> u64 {
    u64::try_from(num_vars).unwrap_or(u64::MAX)
}

/// The error of a pass that gave more variants than a `usize` counts, which
/// is 4294967295 in WebAssembly.
fn the_variants_are_too_many() -> Error {
    Error::KinshipVariantsTooLarge {
        problem: KinshipTooLarge::Variants,
    }
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::arithmetic_side_effects,
        reason = "small literals in tests: the sizes and the patterns of the fixtures"
    )]

    use std::io::Cursor;
    use std::path::{Path, PathBuf};

    use super::{
        Denominators, Kinship, KinshipTooLarge, MAX_INDIVIDUALS_OF_THE_VARIANTS,
        TheBuffersOfTheDenominators, ThePass, add_the_called_genotypes, calc_kinship,
        the_called_genotypes_of, the_counted_denominators_of,
        the_counts_are_cheaper_than_the_product, the_denominators_of_the_block,
        the_pass_over_the_blocks, the_sum_of_the_squares_of_the_missing, the_values_of,
    };
    use crate::block::{Block, BlockReader};
    use crate::error::{Error, Result};
    use crate::filters::FilteringStats;
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{
        ChromTable, DosageOptions, DosageScale, MISSING_ALLELE as MISSING, Needs,
    };

    /// What the kinship asks of the pass over a row: the divisor its allele
    /// frequency gives a variant under Hardy Weinberg, and a variant of more
    /// than two alleles refused.
    const THE_DOSAGES: DosageOptions = DosageOptions {
        transform_to_biallelic: false,
        scale: DosageScale::OfHardyWeinberg,
    };

    /// A reader over blocks a test built, of `num_individuals` individuals
    /// named `ind0` and on.
    ///
    /// It is how the pass is driven over more than one block: `calc_kinship`
    /// puts `Reblock` before it, and `Reblock` joins anything under 10000
    /// variants into one block, so no VCF small enough for a test reaches
    /// the second block of a pass.
    struct GivenBlocks {
        individuals: Vec<String>,
        ploidy: usize,
        chroms: ChromTable,
        /// The blocks it has not given yet, the next one last.
        left: Vec<Block>,
    }

    impl GivenBlocks {
        fn of(blocks: Vec<Block>, num_individuals: usize, ploidy: usize) -> GivenBlocks {
            let mut left = blocks;
            left.reverse();
            GivenBlocks {
                individuals: (0..num_individuals).map(|at| format!("ind{at}")).collect(),
                ploidy,
                chroms: ChromTable::new(),
                left,
            }
        }
    }

    impl BlockReader for GivenBlocks {
        fn next_block(&mut self) -> Result<Option<Block>> {
            Ok(self.left.pop())
        }

        fn individuals(&self) -> &[String] {
            &self.individuals
        }

        fn ploidy(&self) -> usize {
            self.ploidy
        }

        fn chroms(&self) -> &ChromTable {
            &self.chroms
        }

        fn set_needs(&mut self, _needs: Needs) {}

        fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
            Vec::new()
        }
    }

    /// The pass over the blocks of a reader of the tests, over every
    /// individual of it.
    fn the_pass_over(blocks: &mut GivenBlocks, num_individuals: usize, ploidy: usize) -> ThePass {
        match the_pass_over_the_blocks(blocks, None, num_individuals, ploidy, &THE_DOSAGES) {
            Ok(pass) => pass,
            Err(error) => panic!("the pass over the blocks: {error}"),
        }
    }

    /// The denominators the pass left, or a panic when it left a count.
    fn the_denominators_of(pass: ThePass) -> Vec<f64> {
        match pass.denominators {
            Denominators::OfThePair(of_the_pairs) => of_the_pairs,
            Denominators::OfEveryPair => {
                panic!("the block with a genotype missing left the denominators a count")
            }
        }
    }

    /// A block of the genotypes given, `num_individuals` individuals of the
    /// ploidy `ploidy`, variant after variant, with no column but the
    /// genotypes.
    fn block_of(num_individuals: usize, ploidy: usize, gts: &[i8]) -> Block {
        let of_a_variant = num_individuals
            .checked_mul(ploidy)
            .expect("the alleles of one variant");
        Block {
            num_vars: gts
                .len()
                .checked_div(of_a_variant)
                .expect("the variants of the block"),
            num_individuals,
            ploidy,
            gts: gts.to_vec(),
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        }
    }

    /// What the table of literals of `docs/specs/kinship.md` is held to, one
    /// unit of the last digit plink2 prints for an entry near 1: its text
    /// holds six significant digits, so it rounds such an entry by up to
    /// 5e-6. It is the bound of those eleven numbers and of nothing else;
    /// the whole of each matrix is compared with the `f64` plink2 holds,
    /// where there is room to find an error of the arithmetic.
    const OF_PLINK2: f64 = 1e-5;

    /// What each of the 40000 entries of a panel is held to against the
    /// `f64` of plink2, as a share of the largest absolute entry of that
    /// matrix, which is 1.23 on both panels.
    ///
    /// It is a share of the matrix and not of the entry because an entry is
    /// a sum of products that cancel: it can be as near 0 as the data makes
    /// it while the rounding of its sum stays where it was, so a bound
    /// relative to the entry asks the smallest entries for an accuracy that
    /// no arithmetic has. Measured over the 40000 entries of each panel on
    /// 24 September 2026, the largest difference as a share of the largest
    /// entry is 3.6e-16 and 4.5e-16 with the linear algebra on Accelerate,
    /// and 3.3e-15 and 2.3e-15 on faer, which `--no-default-features` and
    /// both wasm targets build. This bound is thirty times the worst of the
    /// four, and it allows 1.2e-13 at the largest entry where a bound of
    /// 1e-12 relative to the entry, which failed on faer, allowed 1.2e-12.
    const OF_THE_BITS_OF_PLINK2: f64 = 1e-13;

    /// The worked example is whole numbers, which pyNei gives within
    /// 4.4e-16, so nothing of it is near this.
    const OF_THE_WORKED_EXAMPLE: f64 = 1e-12;

    /// The path of one of the files of `tests/reference/kinship/`.
    pub(super) fn the_reference_path(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/kinship")
            .join(name)
    }

    /// The path of one of the files of `tests/reference/dists/`, where the
    /// panel with 3 in 100 of its genotypes missing is.
    pub(super) fn the_dists_path(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/dists")
            .join(name)
    }

    /// The VCF of the genotypes given, one line for each variant, with the
    /// individuals named `i0`, `i1` and so on.
    ///
    /// Every line lists the alleles `A` and `C,G,T`, whichever of them its
    /// genotypes hold: the alleles of a variant are the ones its genotypes
    /// hold and not the ones the VCF lists.
    pub(super) fn vcf_of(num_individuals: usize, rows: &[Vec<String>]) -> Vec<u8> {
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

    /// The genotypes of one variant, written as the VCF writes them.
    pub(super) fn variant(genotypes: &[&str]) -> Vec<String> {
        genotypes
            .iter()
            .map(|genotype| (*genotype).to_string())
            .collect()
    }

    /// A reader over the bytes of a VCF, in blocks of `num_vars_per_block`
    /// variants or of the size popnei chooses.
    pub(super) fn reader_over(
        vcf: &[u8],
        num_vars_per_block: Option<usize>,
    ) -> VcfReader<Cursor<Vec<u8>>> {
        let options = VcfOptions {
            num_vars_per_block,
            ..VcfOptions::default()
        };
        match VcfReader::new(Cursor::new(vcf.to_vec()), options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        }
    }

    /// A reader over the bytes of a VCF of that ploidy, which the reader is
    /// told and does not read from the file.
    fn reader_of_the_ploidy(vcf: &[u8], ploidy: usize) -> VcfReader<Cursor<Vec<u8>>> {
        let options = VcfOptions {
            ploidy,
            ..VcfOptions::default()
        };
        match VcfReader::new(Cursor::new(vcf.to_vec()), options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        }
    }

    /// The kinship of every individual of a VCF held in memory, with a
    /// variant of more than two alleles refused.
    pub(super) fn the_kinship_of(vcf: &[u8], num_vars_per_block: Option<usize>) -> Kinship {
        let mut reader = reader_over(vcf, num_vars_per_block);
        match calc_kinship(&mut reader, None, false) {
            Ok(kinship) => kinship,
            Err(error) => panic!("the kinship was not taken: {error}"),
        }
    }

    /// What the kinship of every individual of a VCF failed with.
    pub(super) fn the_kinship_refused(vcf: &[u8]) -> Error {
        let mut reader = reader_over(vcf, None);
        match calc_kinship(&mut reader, None, false) {
            Ok(kinship) => panic!("the kinship was taken over {} variants", kinship.num_vars),
            Err(error) => error,
        }
    }

    /// The entry of the pair of individuals at those two positions.
    pub(super) fn entry_of(kinship: &Kinship, one: usize, other: usize) -> f64 {
        let at = one
            .checked_mul(kinship.num_individuals)
            .and_then(|row| row.checked_add(other))
            .expect("the place of the pair in the matrix");
        kinship.matrix[at]
    }

    /// The kinship of one of the two reference panels, with the names of
    /// its individuals in the order of the file.
    ///
    /// Both are 200 individuals, `s000` to `s199`, and 1200 biallelic
    /// diploid variants: `panel_called.vcf.gz` with every genotype called
    /// and `panel.vcf.gz` of `docs/specs/dists.md` with 3 in 100 of them
    /// missing whole.
    fn the_kinship_of_the_panel(path: &Path) -> (Vec<String>, Kinship) {
        let options = VcfOptions {
            ploidy: 2,
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::from_path(path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let individuals = reader.individuals().to_vec();
        match calc_kinship(&mut reader, None, false) {
            Ok(kinship) => (individuals, kinship),
            Err(error) => panic!("{path}: {error}", path = path.display()),
        }
    }

    /// The panel every genotype of which is called.
    pub(super) fn the_panel_called() -> (Vec<String>, Kinship) {
        the_kinship_of_the_panel(&the_reference_path("panel_called.vcf.gz"))
    }

    /// The same panel with 3 in 100 of its genotypes missing whole.
    pub(super) fn the_panel_with_genotypes_missing() -> (Vec<String>, Kinship) {
        the_kinship_of_the_panel(&the_dists_path("panel.vcf.gz"))
    }

    /// The kinship of that panel with the reader giving `num_vars_per_block`
    /// variants at a time, where `the_panel_with_genotypes_missing` lets the
    /// reader choose the size.
    fn the_kinship_of_the_panel_in_blocks_of(num_vars_per_block: usize) -> Kinship {
        let options = VcfOptions {
            ploidy: 2,
            num_vars_per_block: Some(num_vars_per_block),
            ..VcfOptions::default()
        };
        let path = the_dists_path("panel.vcf.gz");
        let mut reader = match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        match calc_kinship(&mut reader, None, false) {
            Ok(kinship) => kinship,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        }
    }

    /// Where the individual of that name is in the file.
    pub(super) fn individual_at(individuals: &[String], name: &str) -> usize {
        match individuals.iter().position(|held| held == name) {
            Some(at) => at,
            None => panic!("{name} is not an individual of the file"),
        }
    }

    /// The 40000 entries of one of the panels as plink2 holds them, the
    /// little endian `f64` of `--make-rel square bin`, row after row.
    ///
    /// The text of `--make-rel square` beside it holds six significant
    /// digits, which rounds an entry near 1 by up to 5e-6: it is what the
    /// table of literals of `docs/specs/kinship.md` is read from, and these
    /// are what the whole of the matrix is compared with.
    fn the_bits_of_plink2(name: &str) -> Vec<f64> {
        use std::io::Read;

        let path = the_reference_path(&format!("{name}.plink2.rel.bin.gz"));
        let file = match std::fs::File::open(&path) {
            Ok(file) => file,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let mut bytes = Vec::new();
        if let Err(error) = flate2::read::GzDecoder::new(file).read_to_end(&mut bytes) {
            panic!("{path}: {error}", path = path.display());
        }
        bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|eight| f64::from_le_bytes(*eight))
            .collect()
    }

    /// Every entry of the matrix is plink2's within
    /// [`OF_THE_BITS_OF_PLINK2`] of the largest absolute entry of it, and
    /// the individuals are `s000` to `s199` in the order plink2 wrote them
    /// in, which is the order of the VCF and what `<name>.plink2.rel.id`
    /// says.
    fn assert_the_matrix_is_plink2s(panel: &(Vec<String>, Kinship), name: &str) {
        let (individuals, kinship) = panel;
        for (at, individual) in individuals.iter().enumerate() {
            assert_eq!(individual, &format!("s{at:03}"), "the individual at {at}");
        }
        let of_plink2 = the_bits_of_plink2(name);
        assert_eq!(
            of_plink2.len(),
            kinship.matrix.len(),
            "the entries of {name}"
        );
        // The largest absolute entry of the matrix, which every entry of it
        // is held to a share of: the diagonal of a panel is near 1, and the
        // entries near 0 are differences of sums that cancel and carry the
        // rounding of those sums and not of themselves.
        let largest = of_plink2
            .iter()
            .fold(0.0_f64, |so_far, value| so_far.max(value.abs()));
        let allowed = OF_THE_BITS_OF_PLINK2 * largest;
        for (at, (entry, expected)) in kinship.matrix.iter().zip(&of_plink2).enumerate() {
            let apart = (entry - expected).abs();
            assert!(
                apart <= allowed,
                "the entry {at} of {name} is {entry} and plink2 has {expected}, {apart} apart, where {allowed} is allowed of the largest entry {largest}"
            );
        }
    }

    /// The entry of the two named individuals is the number plink2 wrote
    /// for them, within one unit of the last digit it prints.
    fn assert_the_entry_is(panel: &(Vec<String>, Kinship), one: &str, other: &str, of_plink2: f64) {
        let (individuals, kinship) = panel;
        let entry = entry_of(
            kinship,
            individual_at(individuals, one),
            individual_at(individuals, other),
        );
        assert!(
            (entry - of_plink2).abs() < OF_PLINK2,
            "the entry of {one} and {other} is {entry} and plink2 gives {of_plink2}"
        );
    }

    /// The worked example of "How it is verified" of `docs/specs/kinship.md`:
    /// 4 diploid individuals and 4 variants, of which the one where every
    /// individual is heterozygous and the one with a single allele have no
    /// variance and are left out. The genotype `./.` of `i2` at `v1` takes
    /// the mean dosage of its variant and is in the denominator of no pair.
    pub(super) fn the_worked_example() -> Vec<Vec<String>> {
        vec![
            variant(&["0/0", "0/1", "1/1", "0/1"]),
            variant(&["0/0", "0/1", "./.", "1/1"]),
            variant(&["0/1", "0/1", "0/1", "0/1"]),
            variant(&["0/0", "0/0", "0/0", "0/0"]),
        ]
    }

    /// The matrix of the worked example, row after row: 2 variants, and the
    /// denominator of every pair with `i2` is 1 where the others have 2.
    const OF_THE_WORKED_EXAMPLE_MATRIX: [f64; 16] = [
        2.0, 0.0, -2.0, -1.0, //
        0.0, 0.0, 0.0, 0.0, //
        -2.0, 0.0, 2.0, 0.0, //
        -1.0, 0.0, 0.0, 1.0,
    ];

    /// Every entry of the matrix is the one the table of the spec gives.
    fn assert_it_is_the_worked_example(kinship: &Kinship) {
        assert_eq!(kinship.num_individuals, 4);
        assert_eq!(kinship.num_vars, 2, "the variants that were used");
        assert_eq!(kinship.num_vars_given, 4, "the variants the reader gave");
        for (at, (entry, expected)) in kinship
            .matrix
            .iter()
            .zip(OF_THE_WORKED_EXAMPLE_MATRIX)
            .enumerate()
        {
            assert!(
                (entry - expected).abs() < OF_THE_WORKED_EXAMPLE,
                "the entry {at} of the matrix is {entry} and the spec gives {expected}"
            );
        }
    }

    #[test]
    fn the_worked_example_gives_the_whole_numbers_of_the_spec() {
        let kinship = the_kinship_of(&vcf_of(4, &the_worked_example()), None);

        assert_it_is_the_worked_example(&kinship);
    }

    /// The blocks the reader gives are joined and cut to one size before
    /// the pass, so the size it was asked for changes no entry. A reader of
    /// one variant to a block gives the pass one block of the four here,
    /// since the size of a block of four individuals is 10000 variants.
    #[test]
    fn the_worked_example_read_one_variant_to_a_block_gives_the_same_matrix() {
        let kinship = the_kinship_of(&vcf_of(4, &the_worked_example()), Some(1));

        assert_it_is_the_worked_example(&kinship);
    }

    /// The blocks of a pass are joined and cut to one size before it, and
    /// that size is 10000 variants for a dataset of few individuals, so a
    /// pass over a dataset this small is one block and the denominators of
    /// one block never meet those of another in it. The three ways a block
    /// adds to them are read here instead, and each of the three blocks
    /// holds a variant that was dropped for having no variance, since a
    /// denominator that grew by every variant of a block and not by the
    /// ones that were used moves every entry of the matrix with nothing to
    /// show it:
    ///
    /// - a block with nothing missing before the first one that has a
    ///   genotype missing, which is a count and no matrix;
    /// - the first block with a genotype missing, which makes the matrix
    ///   with the variants of the blocks before it in every entry. Its
    ///   variant that was dropped is the first of the two and is called in
    ///   everyone, so a pass that read the rows of the block as they lie
    ///   and not the ones that were used would count that variant in the
    ///   pairs of `i2` and give them 3;
    /// - a block with nothing missing after it, which adds the variants it
    ///   used to every pair alike.
    #[test]
    fn the_denominators_of_the_blocks_carry_the_blocks_before_them() {
        // Three individuals of two alleles each, and two variants in each
        // block. Every genotype of the first block is called and both of
        // its variants are used; the second block drops its first variant
        // and its second is not called in `i2`; the third block drops its
        // second variant and is called in everyone.
        let called = block_of(3, 2, &[0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1]);
        let missing = block_of(3, 2, &[0, 0, 0, 0, 0, 0, 0, 0, 1, 1, MISSING, MISSING]);
        let called_again = block_of(3, 2, &[0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]);
        let mut of_the_blocks = Denominators::OfEveryPair;
        let mut buffers = TheBuffersOfTheDenominators::default();

        the_denominators_of_the_block(
            &called,
            &[true, true],
            2,
            3,
            0,
            &mut buffers,
            &mut of_the_blocks,
        )
        .expect("the denominators of the first block");
        assert!(
            matches!(of_the_blocks, Denominators::OfEveryPair),
            "a block with nothing missing before the first one that has a genotype missing is a count"
        );

        the_denominators_of_the_block(
            &missing,
            &[false, true],
            1,
            3,
            2,
            &mut buffers,
            &mut of_the_blocks,
        )
        .expect("the denominators of the second block");
        the_denominators_of_the_block(
            &called_again,
            &[true, false],
            1,
            3,
            3,
            &mut buffers,
            &mut of_the_blocks,
        )
        .expect("the denominators of the third block");

        // Four variants were used. All four are called in both of `i0` and
        // `i1`, and the three that are called in `i2` are called in both of
        // it and any other: the entry of a pair with `i2` is 3 and of the
        // others 4.
        let of_the_pairs = match of_the_blocks {
            Denominators::OfThePair(of_the_pairs) => of_the_pairs,
            Denominators::OfEveryPair => {
                panic!("the block with a genotype missing left the denominators a count")
            }
        };
        for (at, expected) in [(0, 4.0), (3, 4.0), (4, 4.0), (6, 3.0), (7, 3.0), (8, 3.0)] {
            let count = of_the_pairs[at];
            assert!(
                (count - expected).abs() < OF_THE_WORKED_EXAMPLE,
                "the entry {at} of the denominators is {count} and the variants called in both are {expected}"
            );
        }
    }

    /// What the blocks before the first one with a genotype missing put in
    /// every entry of the denominators is the variants they **used** and not
    /// the variants they gave.
    ///
    /// The pass is driven over two blocks here, since `calc_kinship` puts
    /// `Reblock` before it and no dataset a test builds reaches a second
    /// block through that. The first block gives two variants of three
    /// individuals and uses one, its first having one allele; the second
    /// gives one variant that `ind2` is not called at. Every pair had the
    /// one variant of the first block, and the pair without `ind2` had the
    /// variant of the second too.
    #[test]
    fn the_denominators_carry_the_variants_the_blocks_used_and_not_the_ones_they_gave() {
        let dropped_and_used = block_of(3, 2, &[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1]);
        let missing = block_of(3, 2, &[0, 0, 1, 1, MISSING, MISSING]);
        let mut blocks = GivenBlocks::of(vec![dropped_and_used, missing], 3, 2);

        let pass = the_pass_over(&mut blocks, 3, 2);

        assert_eq!(pass.num_vars, 2, "the variants that were used");
        assert_eq!(pass.num_vars_given, 3, "the variants the reader gave");
        let of_the_pairs = the_denominators_of(pass);
        for (at, expected) in [(0, 2.0), (3, 2.0), (4, 2.0), (6, 1.0), (7, 1.0), (8, 1.0)] {
            let count = of_the_pairs[at];
            assert!(
                (count - expected).abs() < OF_THE_WORKED_EXAMPLE,
                "the entry {at} of the denominators is {count} and the variants called in both are {expected}"
            );
        }
    }

    /// A block of `num_individuals` individuals of the ploidy 2 and
    /// `num_vars` variants in which the genotype of the individual `ind`
    /// of the variant `var` is missing when `(var * 5 + ind * 7) % period`
    /// is 0, and is `0|0` or `0|1` otherwise.
    ///
    /// A `period` of 32 leaves 4 of 128 individuals missing in each
    /// variant, and one of 3 leaves 42 or 43 of them.
    fn a_block_with_a_genotype_missing_every(
        num_individuals: usize,
        num_vars: usize,
        period: usize,
    ) -> Block {
        let mut gts = Vec::with_capacity(num_individuals * num_vars * 2);
        for var in 0..num_vars {
            for ind in 0..num_individuals {
                if (var * 5 + ind * 7) % period == 0 {
                    gts.push(MISSING);
                    gts.push(MISSING);
                } else {
                    gts.push(0);
                    gts.push(i8::from((var + ind) % 3 == 0));
                }
            }
        }
        block_of(num_individuals, 2, &gts)
    }

    /// The two routes to the denominators of a block leave the same
    /// matrix, entry for entry and not within a tolerance.
    ///
    /// Both add whole numbers far below 2^53 to the same `f64`: the
    /// product of the genotypes that were called sums the 1 and the 0 of
    /// each variant, and the counts add the variants of the block less
    /// the ones either individual is missing plus the ones both are.
    /// Whole numbers that small add exactly in an `f64` in any order, so a
    /// difference of one bit here is a defect and not rounding.
    ///
    /// Two blocks are read, one on each side of where the two routes
    /// cross, and each is checked to take the route it should: the first
    /// has 4 of its 128 individuals missing in each variant and goes to
    /// the counts, the second has 42 or 43 and goes to the product. Both
    /// drop two of their ten variants, so a route that counted the
    /// variants it was given and not the ones that were used would show
    /// here.
    #[test]
    fn the_two_routes_to_the_denominators_agree_entry_for_entry() {
        let num_individuals = 128;
        let used = [true, true, false, true, true, true, true, false, true, true];
        let kept = used.iter().filter(|was_used| **was_used).count();
        // What the blocks before this one used, which both matrices start
        // at in every entry, as the pass builds them.
        let before = 7.0;

        for (period, the_counts_are_cheaper) in [(32_usize, true), (3, false)] {
            let block = a_block_with_a_genotype_missing_every(num_individuals, used.len(), period);
            let alleles_per_var = block.alleles_per_var().expect("the alleles of one variant");
            let squares = the_sum_of_the_squares_of_the_missing(&block, &used, alleles_per_var);
            assert!(
                squares > 0,
                "no genotype of the block of one missing every {period} is missing"
            );
            assert_eq!(
                the_counts_are_cheaper_than_the_product(squares, kept, num_individuals),
                the_counts_are_cheaper,
                "which route the block of a genotype missing every {period} takes, \
                 whose sum of the squares of the missing is {squares}"
            );

            let mut of_the_counts = vec![before; the_values_of(num_individuals)];
            let mut of_the_product = vec![before; the_values_of(num_individuals)];
            let mut buffers = TheBuffersOfTheDenominators::default();
            the_counted_denominators_of(
                &block,
                &used,
                alleles_per_var,
                kept,
                num_individuals,
                &mut buffers,
                &mut of_the_counts,
            );
            let mut called = Vec::new();
            the_called_genotypes_of(
                &block,
                &used,
                alleles_per_var,
                kept,
                num_individuals,
                &mut called,
            );
            add_the_called_genotypes(&called, kept, num_individuals, &mut of_the_product)
                .expect("the product of the genotypes that were called with themselves");

            // How many pairs lost a variant to a missing genotype: a
            // fixture where none did would have the two routes agree on
            // the count of every pair alike and show nothing.
            let mut pairs_that_lost_one = 0_usize;
            for one in 0..num_individuals {
                for other in 0..=one {
                    let at = one * num_individuals + other;
                    assert_eq!(
                        of_the_counts[at].to_bits(),
                        of_the_product[at].to_bits(),
                        "the entry of the pair {one}, {other} of the block of a genotype \
                         missing every {period}: the counts give {} and the product {}",
                        of_the_counts[at],
                        of_the_product[at]
                    );
                    if of_the_counts[at] < before + kept as f64 {
                        pairs_that_lost_one += 1;
                    }
                }
            }
            assert!(
                pairs_that_lost_one > 0,
                "every pair of the block of a genotype missing every {period} had every \
                 variant, so the two routes agreeing says nothing"
            );
        }
    }

    /// The error of a variant with more than two alleles names its place
    /// among the variants the reader gave, and not among the ones that were
    /// used: a user looks for it in their file, where the variants that were
    /// dropped are too. The first block here uses one of its two variants
    /// and the variant of three alleles is the first of the second block, so
    /// the two places are 2 and 1.
    #[test]
    fn a_variant_of_more_than_two_alleles_is_named_by_its_place_among_those_given() {
        let dropped_and_used = block_of(3, 2, &[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1]);
        let three_alleles = block_of(3, 2, &[0, 0, 1, 2, 1, 1]);
        let mut blocks = GivenBlocks::of(vec![dropped_and_used, three_alleles], 3, 2);

        let error = match the_pass_over_the_blocks(&mut blocks, None, 3, 2, &THE_DOSAGES) {
            Ok(pass) => panic!("the pass used {} variants", pass.num_vars),
            Err(error) => error,
        };

        assert!(
            matches!(
                error,
                Error::VariantWithMoreThanTwoAlleles {
                    position: 2,
                    num_alleles: 3
                }
            ),
            "{error}"
        );
    }

    /// The frequencies, the means and the denominators are of the
    /// individuals the kinship was asked for: `v0` has a mean dosage of 0.5
    /// over `i0` and `i3` where it has 1 over the four, so the two entries
    /// are not the ones the matrix of the four has for the same pair.
    #[test]
    fn the_kinship_of_two_individuals_takes_their_own_frequencies() {
        let vcf = vcf_of(4, &the_worked_example());
        let mut reader = reader_over(&vcf, None);

        let kinship = match calc_kinship(&mut reader, Some(&[0, 3]), false) {
            Ok(kinship) => kinship,
            Err(error) => panic!("the kinship was not taken: {error}"),
        };

        assert_eq!(kinship.num_individuals, 2);
        assert_eq!(kinship.num_vars, 2);
        let of_pynei = 4.0 / 3.0;
        for (at, (entry, expected)) in kinship
            .matrix
            .iter()
            .zip([of_pynei, -of_pynei, -of_pynei, of_pynei])
            .enumerate()
        {
            assert!(
                (entry - expected).abs() < OF_THE_WORKED_EXAMPLE,
                "the entry {at} of the matrix is {entry} and pyNei gives {expected}"
            );
        }
    }

    /// The third worked example of "How it is verified" of
    /// `docs/specs/kinship.md`, the one of a ploidy that is not 2: 3
    /// individuals and 2 tetraploid variants, where the ploidy is in the
    /// divisor of each variant twice and the genotype of `i2` at `v1` has
    /// two of its four alleles missing, so it is missing whole and is in the
    /// denominator of no pair. The numbers are pyNei's, read from an array
    /// of genotypes because its VCF parser refuses a ploidy above 2.
    #[test]
    fn the_worked_example_of_four_alleles_to_a_genotype_gives_the_numbers_of_pynei() {
        let vcf = vcf_of(
            3,
            &[
                variant(&["0/0/0/0", "0/0/1/1", "1/1/1/1"]),
                variant(&["0/0/0/0", "0/0/0/1", "0/0/./."]),
            ],
        );
        let mut reader = reader_of_the_ploidy(&vcf, 4);

        let kinship = match calc_kinship(&mut reader, None, false) {
            Ok(kinship) => kinship,
            Err(error) => panic!("the kinship was not taken: {error}"),
        };

        assert_eq!(kinship.num_vars, 2);
        let of_pynei = [
            16.0 / 7.0,
            -2.0 / 7.0,
            -4.0,
            -2.0 / 7.0,
            2.0 / 7.0,
            0.0,
            -4.0,
            0.0,
            4.0,
        ];
        for (at, (entry, expected)) in kinship.matrix.iter().zip(of_pynei).enumerate() {
            assert!(
                (entry - expected).abs() < OF_THE_WORKED_EXAMPLE,
                "the entry {at} of the matrix is {entry} and pyNei gives {expected}"
            );
        }
    }

    #[test]
    fn the_variants_used_of_both_panels_are_the_1200_of_the_files() {
        let (_, called) = the_panel_called();
        let (_, missing) = the_panel_with_genotypes_missing();

        assert_eq!(called.num_vars, 1200);
        assert_eq!(missing.num_vars, 1200);
        assert_eq!(called.num_individuals, 200);
        assert_eq!(missing.num_individuals, 200);
    }

    /// The whole of the matrix, all 40000 entries, against the `f64` plink2
    /// wrote, which is the check of the arithmetic; the eleven tests below
    /// are of the eleven literals of the spec and of the six digits its
    /// table holds.
    #[test]
    fn the_whole_matrix_with_every_genotype_called_is_plink2s() {
        assert_the_matrix_is_plink2s(&the_panel_called(), "panel_called");
    }

    #[test]
    fn the_whole_matrix_with_genotypes_missing_is_plink2s() {
        assert_the_matrix_is_plink2s(&the_panel_with_genotypes_missing(), "panel");
    }

    #[test]
    fn the_diagonal_of_s000_with_every_genotype_called_is_plink2s() {
        assert_the_entry_is(&the_panel_called(), "s000", "s000", 1.09309);
    }

    #[test]
    fn the_diagonal_of_s001_with_every_genotype_called_is_plink2s() {
        assert_the_entry_is(&the_panel_called(), "s001", "s001", 1.22825);
    }

    #[test]
    fn the_full_sibs_s000_and_s001_with_every_genotype_called_are_plink2s() {
        assert_the_entry_is(&the_panel_called(), "s000", "s001", 0.648081);
    }

    #[test]
    fn the_full_sibs_s000_and_s002_with_every_genotype_called_are_plink2s() {
        assert_the_entry_is(&the_panel_called(), "s000", "s002", 0.615611);
    }

    #[test]
    fn the_unrelated_s000_and_s004_with_every_genotype_called_are_plink2s() {
        assert_the_entry_is(&the_panel_called(), "s000", "s004", -0.0945533);
    }

    #[test]
    fn the_unrelated_s000_and_s199_with_every_genotype_called_are_plink2s() {
        assert_the_entry_is(&the_panel_called(), "s000", "s199", -0.0760273);
    }

    #[test]
    fn the_full_sibs_s100_and_s101_with_every_genotype_called_are_plink2s() {
        assert_the_entry_is(&the_panel_called(), "s100", "s101", 0.604995);
    }

    #[test]
    fn the_diagonal_of_s000_with_genotypes_missing_is_plink2s() {
        assert_the_entry_is(&the_panel_with_genotypes_missing(), "s000", "s000", 1.09626);
    }

    #[test]
    fn the_full_sibs_s000_and_s001_with_genotypes_missing_are_plink2s() {
        assert_the_entry_is(
            &the_panel_with_genotypes_missing(),
            "s000",
            "s001",
            0.650379,
        );
    }

    #[test]
    fn the_unrelated_s000_and_s004_with_genotypes_missing_are_plink2s() {
        assert_the_entry_is(
            &the_panel_with_genotypes_missing(),
            "s000",
            "s004",
            -0.103505,
        );
    }

    #[test]
    fn the_full_sibs_s100_and_s101_with_genotypes_missing_are_plink2s() {
        assert_the_entry_is(
            &the_panel_with_genotypes_missing(),
            "s100",
            "s101",
            0.604119,
        );
    }

    /// Every entry above the diagonal holds the entry below it.
    ///
    /// The two halves are not worked out twice: the linear algebra writes
    /// the lower half of the products and the pass copies it into the upper
    /// one, so what this reads is that the copy reached every entry. A
    /// matrix given away with the upper half as it was allocated would hold
    /// 40000 zeros there and fail here, and a copy that missed a row or a
    /// column would fail at its first entry.
    #[test]
    fn the_matrix_of_the_panel_is_symmetric() {
        let (_, kinship) = the_panel_with_genotypes_missing();
        let num_individuals = kinship.num_individuals;

        for one in 0..num_individuals {
            for other in 0..num_individuals {
                let entry = entry_of(&kinship, one, other);
                let mirrored = entry_of(&kinship, other, one);
                assert!(
                    (entry - mirrored).abs() < OF_THE_WORKED_EXAMPLE,
                    "the entry of {one} and {other} is {entry} and of {other} and {one} {mirrored}"
                );
            }
        }
    }

    /// The size of block the reader gave changes no bit of the matrix,
    /// because `reblock` joins and cuts the blocks to one size before the
    /// pass: the panel is 1200 variants and the size for 200 individuals is
    /// 10000, so the pass reads one block whether the reader gave 37
    /// variants at a time or the size it chooses.
    ///
    /// The two are compared bit for bit and not within a tolerance. The sum
    /// over the blocks is in floating point, so a pass that saw the blocks
    /// of the reader as they came would give other last bits here, and that
    /// is the whole of what this reads: 0 of the 40000 entries differ.
    #[test]
    fn the_panel_read_in_blocks_of_37_variants_gives_the_matrix_bit_for_bit() {
        let of_37 = the_kinship_of_the_panel_in_blocks_of(37);
        let of_the_default = the_panel_with_genotypes_missing().1;

        assert_eq!(of_37.num_vars, of_the_default.num_vars);
        assert_eq!(of_37.num_vars_given, of_the_default.num_vars_given);
        assert_eq!(of_37.matrix.len(), of_the_default.matrix.len());
        for (at, (entry, of_the_default)) in
            of_37.matrix.iter().zip(&of_the_default.matrix).enumerate()
        {
            assert_eq!(
                entry.to_bits(),
                of_the_default.to_bits(),
                "the entry {at} is {entry} in blocks of 37 and {of_the_default} in the size popnei chooses"
            );
        }
    }

    /// A variant with three alleles among its called genotypes is an error,
    /// and `transform_to_biallelic` reads it with every allele that is not
    /// the major one counting the same, as the principal components of the
    /// variants do.
    #[test]
    fn a_variant_of_more_than_two_alleles_is_refused_unless_it_is_read_as_biallelic() {
        let vcf = vcf_of(
            4,
            &[
                variant(&["0/0", "0/1", "1/1", "0/1"]),
                variant(&["0/0", "0/1", "1/2", "1/1"]),
            ],
        );
        let mut reader = reader_over(&vcf, None);

        let error = match calc_kinship(&mut reader, None, false) {
            Ok(kinship) => panic!("the kinship was taken over {} variants", kinship.num_vars),
            Err(error) => error,
        };

        let message = error.to_string();
        assert!(
            matches!(
                error,
                Error::VariantWithMoreThanTwoAlleles {
                    position: 1,
                    num_alleles: 3
                }
            ),
            "{message}"
        );
        assert!(message.contains("transform_to_biallelic"), "{message}");

        let mut reader = reader_over(&vcf, None);
        let kinship = match calc_kinship(&mut reader, None, true) {
            Ok(kinship) => kinship,
            Err(error) => panic!("the kinship was not taken: {error}"),
        };
        assert_eq!(kinship.num_vars, 2);
    }

    /// A pass that gave no variant is the error every consumer of a pass
    /// raises, with the counts of the filters in it.
    #[test]
    fn a_pass_that_gave_no_variant_is_the_error_every_consumer_raises() {
        let error = the_kinship_refused(&vcf_of(4, &[]));

        assert!(
            matches!(
                &error,
                Error::PassGaveNoVariant {
                    num_vars_of_the_source: 0,
                    filters,
                } if filters.is_empty()
            ),
            "{error}"
        );
    }

    /// A source with no individual has no pair to give. No reader of popnei
    /// has one, so what reaches it here is a kinship asked for none of the
    /// individuals of the reader.
    #[test]
    fn a_source_with_no_individual_is_refused() {
        let mut reader = reader_over(&vcf_of(4, &the_worked_example()), None);

        let error = match calc_kinship(&mut reader, Some(&[]), false) {
            Ok(kinship) => panic!("the kinship was taken over {} variants", kinship.num_vars),
            Err(error) => error,
        };

        assert!(matches!(error, Error::KinshipNoIndividual), "{error}");
    }

    /// More individuals than the individuals x individuals matrix holds
    /// values for are refused at the entry, before a block is read: the
    /// matrix of 46341 would hold more than the 2147483647 values the
    /// linear algebra counts in. The reader here has four individuals and
    /// is never asked for a block.
    #[test]
    fn more_individuals_than_a_matrix_of_them_counts_in_are_refused_at_the_entry() {
        let too_many = MAX_INDIVIDUALS_OF_THE_VARIANTS
            .checked_add(1)
            .expect("one individual more than the largest matrix holds");
        let mut reader = reader_over(&vcf_of(4, &the_worked_example()), None);

        let error = match calc_kinship(&mut reader, Some(&vec![0; too_many]), false) {
            Ok(kinship) => panic!("the kinship was taken over {} variants", kinship.num_vars),
            Err(error) => error,
        };

        let message = error.to_string();
        assert!(
            matches!(
                error,
                Error::KinshipVariantsTooLarge {
                    problem: KinshipTooLarge::Individuals(num_individuals),
                } if num_individuals == too_many
            ),
            "{message}"
        );
        assert!(message.contains("46341 individuals"), "{message}");
    }

    /// A ploidy above the 254 a dosage is written at is refused by the pass
    /// over a row, which this one gets the error from and writes none of its
    /// own: the message names no calculation, since the two that walk that
    /// pass both raise it.
    #[test]
    fn a_ploidy_above_what_a_dosage_is_written_at_is_refused_by_the_row() {
        let of_255_alleles = ["0"; 255].join("/");
        let vcf = vcf_of(2, &[variant(&[&of_255_alleles, &of_255_alleles])]);
        let mut reader = reader_of_the_ploidy(&vcf, 255);

        let error = match calc_kinship(&mut reader, None, false) {
            Ok(kinship) => panic!("the kinship was taken over {} variants", kinship.num_vars),
            Err(error) => error,
        };

        let message = error.to_string();
        assert!(
            matches!(error, Error::VariantPloidyTooLarge { ploidy: 255 }),
            "{message}"
        );
        assert!(!message.contains("kinship"), "{message}");
        assert!(!message.contains("principal component"), "{message}");
    }

    /// An individual that is not in the dataset is refused by the block,
    /// with its position and how many individuals there are.
    #[test]
    fn an_individual_that_is_not_in_the_dataset_is_refused() {
        let mut reader = reader_over(&vcf_of(4, &the_worked_example()), None);

        let error = match calc_kinship(&mut reader, Some(&[0, 4]), false) {
            Ok(kinship) => panic!("the kinship was taken over {} variants", kinship.num_vars),
            Err(error) => error,
        };

        assert!(
            matches!(
                error,
                Error::IndividualToKeepNotInTheBlock {
                    individual: 4,
                    num_individuals: 4
                }
            ),
            "{error}"
        );
    }
}

/// The cases of "Missing genotypes, variants with no variance, and what
/// pyNei asserts" of `docs/specs/kinship.md`, one test each.
#[cfg(test)]
mod cases {
    use super::calc_kinship;
    use super::tests::{
        entry_of, reader_over, the_kinship_of, the_kinship_refused, variant, vcf_of,
    };
    use crate::error::Error;

    /// A genotype with one allele missing is missing whole: it takes the
    /// mean dosage of its variant and is in the denominator of no pair, so
    /// the worked example with `0/.` in the place of its `./.` gives the
    /// same matrix.
    #[test]
    fn a_half_called_genotype_is_missing_as_a_whole_one_is() {
        let half_called = the_kinship_of(
            &vcf_of(
                4,
                &[
                    variant(&["0/0", "0/1", "1/1", "0/1"]),
                    variant(&["0/0", "0/1", "0/.", "1/1"]),
                    variant(&["0/1", "0/1", "0/1", "0/1"]),
                    variant(&["0/0", "0/0", "0/0", "0/0"]),
                ],
            ),
            None,
        );

        assert_eq!(half_called.num_vars, 2);
        // The denominator of every pair with `i2` is 1 and of the others 2,
        // which is the matrix of the worked example.
        assert!((entry_of(&half_called, 0, 2) + 2.0).abs() < 1e-12);
        assert!((entry_of(&half_called, 0, 3) + 1.0).abs() < 1e-12);
        assert!((entry_of(&half_called, 2, 3)).abs() < 1e-12);
    }

    /// A variant where every individual is heterozygous has one dosage
    /// among its called genotypes, so it has no variance and is left out,
    /// although its major allele frequency is 0.5 and no filter by
    /// frequency catches it. Keeping it would add 0 to every entry and 1 to
    /// every denominator.
    #[test]
    fn a_variant_where_every_individual_is_heterozygous_is_left_out() {
        let kinship = the_kinship_of(
            &vcf_of(
                4,
                &[
                    variant(&["0/0", "0/1", "1/1", "0/1"]),
                    variant(&["0/1", "0/1", "0/1", "0/1"]),
                ],
            ),
            None,
        );

        assert_eq!(kinship.num_vars, 1);
        // The one variant that was used gives `i0` a standardized dosage of
        // -sqrt(2) and `i2` one of sqrt(2), over a denominator of 1.
        assert!((entry_of(&kinship, 0, 2) + 2.0).abs() < 1e-12);
    }

    /// Two individuals with no variant called in both are refused, with
    /// their positions and how many variants each of them has called. pyNei
    /// divides by 0 and leaves the NaN in the matrix.
    ///
    /// `i0` is called at two of the three variants and `i2` at one, so the
    /// two counts of the message are different numbers and the one that
    /// belongs to each of the two is read.
    #[test]
    fn a_pair_with_no_variant_called_in_both_is_refused_naming_the_two() {
        let error = the_kinship_refused(&vcf_of(
            3,
            &[
                variant(&["0/0", "0/1", "./."]),
                variant(&["./.", "0/1", "1/1"]),
                variant(&["0/0", "1/1", "./."]),
            ],
        ));

        let message = error.to_string();
        assert!(
            matches!(
                error,
                Error::KinshipPairWithNoVariantCalled {
                    one: 0,
                    other: 2,
                    num_vars_of_one: 2,
                    num_vars_of_other: 1,
                }
            ),
            "{message}"
        );
        assert!(
            message.contains("2 variants are called in the first and 1 in the second"),
            "{message}"
        );
    }

    /// An individual with no called genotype at all, which is a sequencing
    /// that failed, is named on its own: every pair it is in has no variant
    /// called in both, and what a user has to do is leave that one
    /// individual out, not one of a pair. It is the entry of that individual
    /// with itself that has no variant, and for `i0` it is the first entry
    /// of the matrix that is read.
    #[test]
    fn an_individual_with_no_called_genotype_is_named_on_its_own() {
        let error = the_kinship_refused(&vcf_of(
            3,
            &[
                variant(&["./.", "0/1", "1/1"]),
                variant(&["./.", "0/0", "1/1"]),
            ],
        ));

        let message = error.to_string();
        assert!(
            matches!(
                error,
                Error::KinshipPairWithNoVariantCalled {
                    one: 0,
                    other: 0,
                    num_vars_of_one: 0,
                    num_vars_of_other: 0,
                }
            ),
            "{message}"
        );
        assert!(
            message.contains("has no called genotype among the variants that were used"),
            "{message}"
        );
        assert!(message.contains("leave it out"), "{message}");
        assert!(!message.contains("one of the two"), "{message}");
    }

    /// A dataset in which no variant varies is refused: every variant has
    /// one dosage among its called genotypes, and a kinship measures a pair
    /// against the average pair of the panel, which such a panel has none
    /// of. pyNei raises "No variant varies among the samples, there is no
    /// kinship".
    #[test]
    fn a_dataset_where_no_variant_varies_is_refused() {
        let error = the_kinship_refused(&vcf_of(
            4,
            &[
                variant(&["0/1", "0/1", "0/1", "0/1"]),
                variant(&["0/0", "0/0", "0/0", "0/0"]),
                variant(&["./.", "./.", "./.", "./."]),
            ],
        ));

        assert!(
            matches!(error, Error::KinshipNoVariantWithVariance),
            "{error}"
        );
    }

    /// The four individuals of the worked example with `i1` alone asked
    /// for: one individual and the variants that vary among them all, of
    /// which none varies in one individual, so there is no kinship.
    #[test]
    fn one_individual_alone_has_no_variant_with_variance() {
        let vcf = vcf_of(
            4,
            &[
                variant(&["0/0", "0/1", "1/1", "0/1"]),
                variant(&["0/0", "0/1", "./.", "1/1"]),
            ],
        );
        let mut reader = reader_over(&vcf, None);

        let error = match calc_kinship(&mut reader, Some(&[1]), false) {
            Ok(kinship) => panic!("the kinship was taken over {} variants", kinship.num_vars),
            Err(error) => error,
        };

        assert!(
            matches!(error, Error::KinshipNoVariantWithVariance),
            "{error}"
        );
    }
}

/// The principal components of a kinship, "The principal components of the
/// kinship" of `docs/specs/kinship.md`.
#[cfg(test)]
mod components {
    use super::tests::{
        individual_at, the_kinship_of, the_panel_called, the_panel_with_genotypes_missing,
        the_worked_example, variant, vcf_of,
    };
    use super::{
        Kinship, KinshipPcs, fix_the_sign_of, principal_components, the_components_with_variance,
    };
    use crate::error::Error;

    /// What the eigenvalues and the projections of numpy 2.5.3 are held to,
    /// relative for an eigenvalue and absolute for a projection, which is
    /// what "How it is verified" of `docs/specs/kinship.md` asks. The
    /// projections of the worked example are between 0.2 and 1.5, so the two
    /// bounds are of the same size there.
    const OF_NUMPY: f64 = 1e-9;

    /// The three largest eigenvalues of `panel_called`, from numpy 2.5.3 on
    /// 23 September 2026, written to 15 digits.
    ///
    /// They had 9, and popnei is 4e-16 of itself from what numpy gives
    /// while the third of them rounded to 9 digits is 8.5e-10 away, 85% of
    /// the 1e-9 they are held to: a change that is right and moves an
    /// eigenvalue by 1.5e-10 would have reddened this and the two suites
    /// that assert the same three numbers.
    const THE_EIGENVALUES_OF_THE_PANEL: [f64; 3] = [
        17.269_141_155_457_5,
        12.447_315_235_850_9,
        3.358_712_577_141_36,
    ];

    /// The projections of the two components of the worked example of "How
    /// it is verified" of `docs/specs/kinship.md`, from numpy 2.5.3 on 23
    /// September 2026, with the sign of the rule: 4 individuals x 2
    /// components, row after row. The two of `i1` are 0, since its
    /// standardized dosage is 0 at both variants that were used.
    const THE_COMPONENTS_OF_THE_WORKED_EXAMPLE: [f64; 8] = [
        1.45989777643,
        -0.212872996577, //
        0.0,
        0.0, //
        -1.34910400096,
        -0.550866821327, //
        -0.461372751067,
        0.937211435408,
    ];

    /// The components of a kinship, or a panic with what they failed with.
    fn the_components_of(kinship: &Kinship, num_pcs: usize) -> KinshipPcs {
        match principal_components(kinship, num_pcs) {
            Ok(pcs) => pcs,
            Err(error) => panic!("the components were not taken: {error}"),
        }
    }

    /// Where one individual falls along one component.
    fn the_projection_of(pcs: &KinshipPcs, individual: usize, component: usize) -> f64 {
        let at = individual
            .checked_mul(pcs.num_comps)
            .and_then(|row| row.checked_add(component))
            .expect("the place of the projection");
        pcs.projections[at]
    }

    /// The sum of the squares of the projections of one component, which is
    /// its eigenvalue: a component is `u_j sqrt(lambda_j)` and the
    /// eigenvector `u_j` has length 1. No function of popnei gives an
    /// eigenvalue, so this is how the eigenvalues are read.
    fn the_sum_of_the_squares_of(
        pcs: &KinshipPcs,
        component: usize,
        num_individuals: usize,
    ) -> f64 {
        (0..num_individuals)
            .map(|individual| the_projection_of(pcs, individual, component).powi(2))
            .sum()
    }

    /// The projection of a component that decides its sign: the one of the
    /// largest absolute value, and the first of them when two are exactly of
    /// one size.
    ///
    /// The rule of `docs/specs/pca.md` takes two projections within 64 units
    /// in the last place of each other for one absolute value, which this
    /// does not. No dataset of this module reaches that tolerance: neither
    /// panel has such a pair, the two largest absolute values of a component
    /// being 2.5e-3 of each other at the closest on `panel_called` and
    /// 4.2e-4 on `panel`, measured with numpy 2.5.3 over the 199 components
    /// of each; and the two projections of the kinship of two individuals
    /// below are one number with opposite signs, which both backends give
    /// with identical bits, so their absolute values are compared and found
    /// equal without it.
    /// `the_tolerance_of_the_sign_rule_keeps_the_first_of_two_that_are_of_one_size`
    /// is where the tolerance decides.
    fn the_projection_that_fixes_the_sign_of(
        pcs: &KinshipPcs,
        component: usize,
        num_individuals: usize,
    ) -> f64 {
        let mut largest = 0.0_f64;
        for individual in 0..num_individuals {
            let projection = the_projection_of(pcs, individual, component);
            if projection.abs() > largest.abs() {
                largest = projection;
            }
        }
        largest
    }

    /// Every component of a kinship has the projection that decides its sign
    /// above 0, which is the rule of `docs/specs/pca.md`.
    fn assert_every_component_obeys_the_sign_rule(kinship: &Kinship, pcs: &KinshipPcs, of: &str) {
        assert!(pcs.num_comps > 0, "the components of {of}");
        for component in 0..pcs.num_comps {
            let largest =
                the_projection_that_fixes_the_sign_of(pcs, component, kinship.num_individuals);
            assert!(
                largest > 0.0,
                "the component {component} of {of} has {largest} as the projection of its largest absolute value"
            );
        }
    }

    /// The kinship of the worked example of `docs/specs/kinship.md`: 4
    /// individuals, of which `i1` is 0 against everyone, and 2 variants with
    /// variance of the 4 the reader gives.
    fn the_kinship_of_the_worked_example() -> Kinship {
        the_kinship_of(&vcf_of(4, &the_worked_example()), None)
    }

    /// Two individuals and one variant, `0/0` and `1/1`: the standardized
    /// dosages are -sqrt(2) and sqrt(2) over a denominator of 1, so the
    /// matrix is 2 on the diagonal and -2 off it, and its one component has
    /// two projections of one absolute value.
    fn the_kinship_of_two_individuals() -> Kinship {
        the_kinship_of(&vcf_of(2, &[variant(&["0/0", "1/1"])]), None)
    }

    /// The sum of the squares of the projections of a component is its
    /// eigenvalue, and the three largest of `panel_called` are
    /// 17.2691411554575, 12.4473152358509 and 3.35871257714136 from numpy
    /// 2.5.3.
    #[test]
    fn the_first_three_components_of_the_panel_hold_the_eigenvalues_numpy_gives() {
        let (_, kinship) = the_panel_called();

        let pcs = the_components_of(&kinship, 3);

        assert_eq!(pcs.num_comps, 3, "the components asked for");
        assert_eq!(pcs.projections.len(), 600, "the individuals x components");
        for (component, eigenvalue) in THE_EIGENVALUES_OF_THE_PANEL.into_iter().enumerate() {
            let sum = the_sum_of_the_squares_of(&pcs, component, kinship.num_individuals);
            assert!(
                (sum - eigenvalue).abs() <= OF_NUMPY * eigenvalue,
                "the squares of the component {component} add up to {sum} and numpy gives {eigenvalue}"
            );
        }
    }

    /// The panel with every genotype called, over every component it has.
    #[test]
    fn every_component_of_the_panel_called_obeys_the_sign_rule() {
        let (_, kinship) = the_panel_called();

        let pcs = the_components_of(&kinship, kinship.num_individuals);

        assert_every_component_obeys_the_sign_rule(&kinship, &pcs, "panel_called");
    }

    /// The panel with 3 in 100 of its genotypes missing, whose per pair
    /// denominators put an eigenvalue below 0.
    #[test]
    fn every_component_of_the_panel_with_genotypes_missing_obeys_the_sign_rule() {
        let (_, kinship) = the_panel_with_genotypes_missing();

        let pcs = the_components_of(&kinship, kinship.num_individuals);

        assert_every_component_obeys_the_sign_rule(&kinship, &pcs, "panel");
    }

    /// The rule of `docs/specs/pca.md` gives the first of two projections of
    /// one absolute value the positive sign. The two individuals of one
    /// variant have projections of sqrt(2) and -sqrt(2), the same number
    /// with opposite signs, so which of the two the eigendecomposition
    /// leaves the larger by a bit is what would decide the sign of the
    /// component without the tolerance of the rule.
    #[test]
    fn the_first_of_two_projections_of_one_absolute_value_is_the_positive_one() {
        let kinship = the_kinship_of_two_individuals();

        let pcs = the_components_of(&kinship, 2);

        assert_eq!(pcs.num_comps, 1, "the second eigenvalue is 0");
        let first = the_projection_of(&pcs, 0, 0);
        let second = the_projection_of(&pcs, 1, 0);
        assert!(
            (first - std::f64::consts::SQRT_2).abs() < OF_NUMPY,
            "the first individual is at {first} and sqrt(2) is asked for"
        );
        assert!(
            (second + std::f64::consts::SQRT_2).abs() < OF_NUMPY,
            "the second individual is at {second} and -sqrt(2) is asked for"
        );
    }

    /// The worked example has 4 individuals and eigenvalues 4.16424794,
    /// 1.22713444, -3.2e-16 and -0.39138238 from numpy 2.5.3, so two of them
    /// are above the tolerance and 6 components asked for give 2. pyNei
    /// gives the 6 that were asked for and raises out of pandas above the
    /// individuals.
    #[test]
    fn a_kinship_asked_for_more_components_than_it_has_gives_the_ones_above_the_tolerance() {
        let kinship = the_kinship_of_the_worked_example();

        let pcs = the_components_of(&kinship, 6);

        assert_eq!(pcs.num_comps, 2, "the components above the tolerance");
        assert_eq!(pcs.projections.len(), 8, "4 individuals x 2 components");
    }

    /// The threshold a component has to be above grows with the side of the
    /// matrix: it is the largest eigenvalue times the side times the
    /// epsilon of an `f64`, so a kinship of 200 individuals cuts at
    /// 4.44e-14 of the largest eigenvalue and one of 5 cuts at 1.11e-15 of
    /// it. An eigenvalue of 1.5e-15 of the largest is a component of the
    /// second and of no component of the first.
    ///
    /// The threshold is `crate::pca`'s, which the components of a kinship
    /// take unchanged, and the side of it had no test: taking it out leaves
    /// every test of the principal components and of the kinship green,
    /// where a threshold of 0 reddens ten of them. It is read here and not
    /// beside the threshold because the 47 tests of `pca.rs` are the
    /// evidence that this plan changed no number of the principal
    /// components, and they keep their names and their count.
    #[test]
    fn the_threshold_of_a_component_grows_with_the_side_of_the_matrix() {
        // From the largest, as an eigendecomposition gives them: one
        // eigenvalue of 1.5e-15 of the largest and one below 0.
        let values = [1.0, 1.5e-15, -2e-16];

        assert_eq!(
            the_components_with_variance(&values, 200, 200),
            1,
            "of 200 individuals, whose threshold is 4.44e-14 of the largest eigenvalue"
        );
        assert_eq!(
            the_components_with_variance(&values, 5, 5),
            2,
            "of 5 individuals, whose threshold is 1.11e-15 of it"
        );
    }

    /// Each projection of the two components of the worked example is the
    /// one numpy 2.5.3 gives, with the sign of the rule, which is what says
    /// that a component is the eigenvector times the square root of its
    /// eigenvalue and that the matrix is individuals x components.
    #[test]
    fn the_components_of_the_worked_example_are_the_ones_numpy_gives() {
        let kinship = the_kinship_of_the_worked_example();

        let pcs = the_components_of(&kinship, 2);

        assert_eq!(pcs.num_comps, 2);
        for (at, (projection, expected)) in pcs
            .projections
            .iter()
            .zip(THE_COMPONENTS_OF_THE_WORKED_EXAMPLE)
            .enumerate()
        {
            assert!(
                (projection - expected).abs() < OF_NUMPY,
                "the projection {at} is {projection} and numpy gives {expected}"
            );
        }
    }

    /// A kinship gives the components that were asked for when it has more
    /// of them: the worked example has 2 above the tolerance, and 1 asked
    /// for is 1, the first of them.
    #[test]
    fn a_kinship_asked_for_fewer_components_than_it_has_gives_that_many() {
        let kinship = the_kinship_of_the_worked_example();

        let pcs = the_components_of(&kinship, 1);

        assert_eq!(pcs.num_comps, 1, "the components asked for");
        assert_eq!(pcs.projections.len(), 4, "4 individuals x 1 component");
        for (at, (projection, expected)) in pcs
            .projections
            .iter()
            .zip(THE_COMPONENTS_OF_THE_WORKED_EXAMPLE.into_iter().step_by(2))
            .enumerate()
        {
            assert!(
                (projection - expected).abs() < OF_NUMPY,
                "the projection {at} of the first component is {projection} and numpy gives {expected}"
            );
        }
    }

    /// A kinship measures each pair against the average pair of the panel,
    /// which takes one direction out of it exactly when no genotype is
    /// missing, since each variant is centered; with the per pair
    /// denominators it is measured and not proved, which
    /// `docs/specs/kinship.md` says. The last eigenvalue is 0 or below it
    /// on both panels: -3.44e-15 on `panel_called` and -0.0321 on `panel`
    /// from numpy 2.5.3, against a tolerance of 7.67e-13 on both.
    /// pyNei gives 200 components on either, the last of them the square
    /// root of the absolute value of that eigenvalue.
    #[test]
    fn both_panels_asked_for_a_component_for_each_individual_give_one_fewer() {
        let (_, called) = the_panel_called();
        let (_, missing) = the_panel_with_genotypes_missing();

        let of_the_called = the_components_of(&called, 200);
        let of_the_missing = the_components_of(&missing, 200);

        assert_eq!(of_the_called.num_comps, 199, "panel_called");
        assert_eq!(of_the_missing.num_comps, 199, "panel");
    }

    /// The projection of `s000` and of `s199` on the first component of
    /// `panel_called`, from numpy 2.5.3 on 24 September 2026 with the sign
    /// rule applied: 0.0506331222853770 and -0.2797284741269570.
    ///
    /// They say which individual each row of the matrix of projections
    /// belongs to, which nothing else of this crate reads: the two tests of
    /// a panel above it are a sum over the rows and the sign of the largest
    /// of them, and both are the same numbers when two individuals are
    /// swapped for each other. Swapping the rows 0 and 1 of the
    /// eigenvectors leaves every other test of the crate passing.
    #[test]
    fn the_first_component_of_the_panel_places_two_named_individuals_where_numpy_does() {
        let (individuals, kinship) = the_panel_called();

        let pcs = the_components_of(&kinship, 10);

        assert_eq!(pcs.num_comps, 10);
        for (name, of_numpy) in [
            ("s000", 0.050_633_122_285_377),
            ("s199", -0.279_728_474_126_957),
        ] {
            let at = individual_at(&individuals, name);
            // The projections are the individuals x the components, row
            // after row, so the first component of an individual is where
            // its row starts.
            let projection = pcs.projections[at * pcs.num_comps];
            assert!(
                (projection - of_numpy).abs() < OF_NUMPY,
                "the first component of {name} is {projection} and numpy gives {of_numpy}"
            );
        }
    }

    /// The rule that fixes the sign of a component calls two projections
    /// one absolute value when they are within 64 units in the last place
    /// of each other, and keeps the first of the two: here the second is 32
    /// units above the first in absolute value, so without that tolerance
    /// the second would decide and the component would be left as it is.
    ///
    /// No kinship of this module reaches the tolerance, so nothing else
    /// here reads it: setting it to 0 leaves every test of the kinship
    /// green. It is the rule of `docs/specs/pca.md`, which the components
    /// of a kinship take unchanged, and it is read here because the 47
    /// tests of `pca.rs` are the evidence that this plan changed no number
    /// of the principal components and they keep their names and their
    /// count.
    #[test]
    fn the_tolerance_of_the_sign_rule_keeps_the_first_of_two_that_are_of_one_size() {
        // One component of two individuals: -1 and one number 32 units in
        // the last place above 1.
        let mut projections = [-1.0, 1.0 + 32.0 * f64::EPSILON];

        let turned = fix_the_sign_of(&mut projections, 0, 1);

        assert!(
            turned,
            "the component is turned round by the first of the two"
        );
        assert!(
            projections[0] > 0.0,
            "the projections are {projections:?} and the first of them decides the sign"
        );
    }

    /// The components come in the order of their eigenvalues, from the
    /// largest, which is what makes `PC0` the direction the panel varies
    /// most along.
    ///
    /// The eigenvalue of a component is the sum of the squares of its
    /// projections, since a component is `u_j * sqrt(lambda_j)` and `u_j`
    /// has length 1. All 199 of a panel are read: two components past the
    /// tenth can be swapped for each other with every literal of every
    /// suite still asserting what it did, and a user who asks for 60
    /// components gets two of them in the wrong order.
    #[test]
    fn the_components_of_a_panel_come_in_the_order_of_their_eigenvalues() {
        let (_, called) = the_panel_called();
        let (_, missing) = the_panel_with_genotypes_missing();

        for (name, kinship) in [("panel_called", called), ("panel", missing)] {
            let pcs = the_components_of(&kinship, 200);
            let mut of_the_component: Vec<f64> = vec![0.0; pcs.num_comps];
            for of_the_individual in pcs.projections.chunks_exact(pcs.num_comps) {
                for (sum, projection) in of_the_component.iter_mut().zip(of_the_individual) {
                    *sum += projection * projection;
                }
            }
            for (at, pair) in of_the_component.windows(2).enumerate() {
                let (this, next) = (pair[0], pair[1]);
                assert!(
                    this >= next,
                    "the eigenvalue of the component {at} of {name} is {this} and the one after it {next}"
                );
            }
        }
    }

    /// Asking for no component is not an error and gives none, as asking a
    /// principal component analysis of the variants for none is not.
    #[test]
    fn a_kinship_asked_for_no_component_gives_none() {
        let kinship = the_kinship_of_the_worked_example();

        let pcs = the_components_of(&kinship, 0);

        assert_eq!(pcs.num_comps, 0);
        assert!(pcs.projections.is_empty(), "{:?}", pcs.projections);
    }

    /// A kinship of no individual leaves nobody to place along a component.
    /// A user reaches it from Python with an empty frame, which the checks
    /// of a kinship built by hand let past, since a matrix of no row is
    /// square and names nobody twice.
    #[test]
    fn a_kinship_of_no_individual_is_refused() {
        let kinship = Kinship {
            num_individuals: 0,
            num_vars: 0,
            num_vars_given: 0,
            matrix: Vec::new(),
        };

        let error = match principal_components(&kinship, 3) {
            Ok(pcs) => panic!("{} components were given", pcs.num_comps),
            Err(error) => error,
        };

        assert!(matches!(error, Error::KinshipNoIndividual), "{error}");
    }

    /// A value of the matrix that is not finite is a wrong matrix, and the
    /// message names the row and the column it is at. The frame of a
    /// `Kinship` is checked when the object is built and a user can write
    /// into it afterwards, so this is where such a value arrives; without
    /// the check the linear algebra refuses it and names a matrix `g`, an
    /// internal name of `crates/popnei-linalg`, and Python calls a wrong
    /// matrix a defect of popnei.
    #[test]
    fn a_value_of_the_matrix_that_is_not_finite_is_refused_with_where_it_is() {
        // The matrix is 4 x 4: the entry 4 is the row 1 and the column 0,
        // below the diagonal and read by the eigendecomposition, and the
        // entry 1 is the row 0 and the column 1, above it and read by
        // nothing.
        for (at, row, col) in [(4, 1, 0), (1, 0, 1)] {
            let mut kinship = the_kinship_of_the_worked_example();
            kinship.matrix[at] = f64::NAN;

            let error = match principal_components(&kinship, 2) {
                Ok(pcs) => panic!("{} components were given", pcs.num_comps),
                Err(error) => error,
            };

            let message = error.to_string();
            assert!(
                matches!(
                    error,
                    Error::KinshipValueNotFinite {
                        row: found_row,
                        col: found_col,
                        value,
                    } if found_row == row && found_col == col && value.is_nan()
                ),
                "{message}"
            );
            assert!(
                message.contains("not finite") || message.contains("finite"),
                "{message}"
            );
        }
    }
}
