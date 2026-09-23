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

use std::fmt;

use popnei_linalg::add_self_product_lower;

use crate::block::{Block, BlockReader, Reblock};
use crate::error::{Error, Result};
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
    // The two buffers of one block, kept from one block to the next so that
    // a pass over a million variants asks for them once. The rows that were
    // not used are left as they were and nothing reads them.
    let mut standardized: Vec<f64> = Vec::new();
    let mut called: Vec<f64> = Vec::new();
    while let Some(mut block) = blocks.next_block()? {
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
                first: usize::try_from(num_vars_given).map_err(|_| the_variants_are_too_many())?,
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
                &mut called,
                &mut denominators,
            )?;
        }
        num_vars = num_vars
            .checked_add(the_count_of(kept))
            .ok_or_else(the_variants_are_too_many)?;
        num_vars_given = num_vars_given
            .checked_add(the_count_of(block.num_vars))
            .ok_or_else(the_variants_are_too_many)?;
    }
    Ok(ThePass {
        gram,
        denominators,
        num_vars,
        num_vars_given,
    })
}

/// The variants of one block that each pair had called in both, added to
/// the denominators of the blocks before it.
///
/// `used` says which variants of the block have variance, `kept` how many
/// of them that is, 1 or more, and `num_vars_before` how many variants the
/// blocks before this one used. `called` is the buffer the genotypes that
/// were called are written into, one value for each individual of each
/// variant that was used, which is kept from one block to the next.
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
    called: &mut Vec<f64>,
    denominators: &mut Denominators,
) -> Result<()> {
    let alleles_per_var = block.alleles_per_var()?;
    let any_missing =
        the_rows_used(block, used, alleles_per_var).any(|row| row.contains(&MISSING_ALLELE));
    if !any_missing {
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
    the_called_genotypes_of(block, used, alleles_per_var, kept, num_individuals, called);
    if let Denominators::OfThePair(of_the_pairs) = denominators {
        return add_the_called_genotypes(called, kept, num_individuals, of_the_pairs);
    }
    // The blocks before this one had no genotype missing, so each of their
    // variants is in the denominator of every pair: the matrix starts at
    // that count in every entry and not at 0, which would lose them.
    let mut of_the_pairs = vec![num_vars_before as f64; the_values_of(num_individuals)];
    add_the_called_genotypes(called, kept, num_individuals, &mut of_the_pairs)?;
    *denominators = Denominators::OfThePair(of_the_pairs);
    Ok(())
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
    use std::io::Cursor;
    use std::path::{Path, PathBuf};

    use super::{
        Denominators, Kinship, ThePass, calc_kinship, the_denominators_of_the_block,
        the_pass_over_the_blocks,
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

    /// The entries of the two reference panels are plink2's within one unit
    /// of the last digit it prints for an entry of that size: it writes six
    /// significant digits, so an entry near 1 is rounded by up to 5e-6.
    const OF_PLINK2: f64 = 1e-5;

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
    fn the_panel_called() -> (Vec<String>, Kinship) {
        the_kinship_of_the_panel(&the_reference_path("panel_called.vcf.gz"))
    }

    /// The same panel with 3 in 100 of its genotypes missing whole.
    fn the_panel_with_genotypes_missing() -> (Vec<String>, Kinship) {
        the_kinship_of_the_panel(&the_dists_path("panel.vcf.gz"))
    }

    /// Where the individual of that name is in the file.
    fn individual_at(individuals: &[String], name: &str) -> usize {
        match individuals.iter().position(|held| held == name) {
            Some(at) => at,
            None => panic!("{name} is not an individual of the file"),
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
    fn the_worked_example() -> Vec<Vec<String>> {
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
        let mut buffer = Vec::new();

        the_denominators_of_the_block(
            &called,
            &[true, true],
            2,
            3,
            0,
            &mut buffer,
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
            &mut buffer,
            &mut of_the_blocks,
        )
        .expect("the denominators of the second block");
        the_denominators_of_the_block(
            &called_again,
            &[true, false],
            1,
            3,
            3,
            &mut buffer,
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

    #[test]
    fn the_variants_used_of_both_panels_are_the_1200_of_the_files() {
        let (_, called) = the_panel_called();
        let (_, missing) = the_panel_with_genotypes_missing();

        assert_eq!(called.num_vars, 1200);
        assert_eq!(missing.num_vars, 1200);
        assert_eq!(called.num_individuals, 200);
        assert_eq!(missing.num_individuals, 200);
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

    /// The matrix is the same whichever of them is the major allele, so the
    /// two halves of the matrix hold the same number for a pair.
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

    /// The blocks are joined and cut to one size before the pass, so the
    /// size the reader gave does not change the variants; the sum over the
    /// blocks is in floating point and the boundaries decide its last bits,
    /// which is why this is compared within the tolerance of plink2 and not
    /// bit by bit.
    #[test]
    fn the_panel_read_in_blocks_of_37_variants_gives_the_same_entries() {
        let options = VcfOptions {
            ploidy: 2,
            num_vars_per_block: Some(37),
            ..VcfOptions::default()
        };
        let path = the_dists_path("panel.vcf.gz");
        let mut reader = match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let individuals = reader.individuals().to_vec();
        let kinship = match calc_kinship(&mut reader, None, false) {
            Ok(kinship) => kinship,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };

        let panel = (individuals, kinship);
        assert_eq!(panel.1.num_vars, 1200);
        assert_the_entry_is(&panel, "s000", "s000", 1.09626);
        assert_the_entry_is(&panel, "s000", "s001", 0.650379);
        assert_the_entry_is(&panel, "s100", "s101", 0.604119);
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
    #[test]
    fn a_pair_with_no_variant_called_in_both_is_refused_naming_the_two() {
        let error = the_kinship_refused(&vcf_of(
            3,
            &[
                variant(&["0/0", "0/1", "./."]),
                variant(&["./.", "0/1", "1/1"]),
            ],
        ));

        let message = error.to_string();
        assert!(
            matches!(
                error,
                Error::KinshipPairWithNoVariantCalled {
                    one: 0,
                    other: 2,
                    num_vars_of_one: 1,
                    num_vars_of_other: 1,
                }
            ),
            "{message}"
        );
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
