//! The distances between individuals, and how they are counted over the
//! genotypes of a block.
//!
//! The Kosman distance of two individuals at one variant, d, is how many of
//! their alleles do not pair with an equal allele of the other over the
//! ploidy, and the distance of the pair is the mean of d over the variants
//! at which both genotypes are called. Both come from two integers, the
//! ploidy times the sum of d and how many variants both were called at,
//! which every block adds to.
//!
//! [`KosmanBits`] is the genotypes of one block as sets of bits, one set
//! for each individual and each question that a count of a pair asks, and
//! it gives those two integers for a pair of individuals over that block.
//! `docs/specs/dists.md` has the design, the formula and the numbers the
//! tests assert, and the row `dists` of section 9 of
//! `docs/architecture.md` is where the module sits.

use std::num::NonZeroUsize;

use crate::block::{Block, BlockSize};
use crate::error::{Error, Result};
use crate::variant::{MISSING_ALLELE, Needs};

/// How many variants one word of a set of bits holds, the bits of a `u64`.
const VARS_PER_WORD: usize = 64;

/// The genotypes of one block as sets of bits, one set for each individual
/// and each question that a count of a pair asks of a variant.
///
/// The questions are: is the genotype of the individual called, which is
/// its `called` set; and, for each allele a that the block holds and each
/// count m from 1 to the ploidy k, does its called genotype hold m copies
/// of a or more, which is its `holds` set of a and m. So an individual has
/// 1 + k * A sets, with A the alleles of the block, one more than the
/// largest allele of its genotypes, and each set has one bit for each
/// variant of the block, set where the answer is yes and 0 where the
/// genotype is missing.
///
/// For a pair of individuals, the variants both were called at are the
/// bits set in both `called` sets, and the bits set in both of their
/// `holds` sets, added over the alleles and the counts, are the alleles of
/// one that pair with an equal allele of the other, added over those
/// variants: k - d at each of them.
///
/// The `called` sets of every individual lie in one array and their
/// `holds` sets in another, the sets of an individual side by side in
/// each, so the two sums of a pair are one pass over two contiguous
/// slices and one pass over two more.
#[derive(Debug)]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "nothing in the crate counts the pairs of a block yet: the calculation over \
                  a reader, which builds these sets for every block it takes, is not written"
    )
)]
pub(crate) struct KosmanBits {
    /// How many individuals the block holds, which is how many the sets
    /// are of.
    num_individuals: usize,
    /// How many alleles the genotype of one individual holds, the k of the
    /// formula, by which the sum of d is multiplied.
    ploidy: u32,
    /// How many `u64` one set holds: the variants of the block spread over
    /// words of 64, the last of them partial. The bits of that last word
    /// beyond the variants of the block are set for nobody.
    words_per_set: NonZeroUsize,
    /// How many `u64` the `holds` sets of one individual hold, k * A sets
    /// of `words_per_set` words, which is where the sets of the next
    /// individual start.
    holds_per_individual: NonZeroUsize,
    /// The `called` set of every individual, one after another.
    called: Vec<u64>,
    /// The `holds` sets of every individual, the sets of an individual one
    /// after another and those of the individuals one after another: the
    /// set of the allele a and the count m is the one at a * k + m - 1.
    holds: Vec<u64>,
}

#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "nothing in the crate counts the pairs of a block yet: the calculation over \
                  a reader, which builds these sets for every block it takes, is not written"
    )
)]
impl KosmanBits {
    /// The sets of bits of the genotypes of `block`.
    ///
    /// # Errors
    ///
    /// When the ploidy of the block is 0, when it holds no variant, when
    /// it holds no genotypes, when its genotypes are not as many as its
    /// size says, when a genotype holds an allele below
    /// [`MISSING_ALLELE`], when the sums of a pair over the block would go
    /// above what a `u32` holds, and when the machine does not give the
    /// memory of the sets.
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "the word of a variant and its bit are a division and a remainder by 64, \
                  which is not 0, and a shift by a remainder of 64, which is below the bits \
                  of a u64; the place of an allele is below the ploidy, which is below the \
                  alleles of the block, a number a usize holds; and the set of an allele and \
                  a count, `allele * ploidy + copies - 1`, is 0 at least, since the copies \
                  are 1 at least, and below `A * k`, since the allele is below A and the \
                  copies are at most the ploidy, so its word is below the words of the \
                  `holds` sets of one individual, which were checked above to be a number a \
                  usize holds"
    )]
    pub(crate) fn of_block(block: &Block) -> Result<KosmanBits> {
        if block.ploidy == 0 {
            return Err(Error::GtsNotWholeGenotypes {
                num_alleles: block.gts.len(),
                ploidy: block.ploidy,
            });
        }
        // Every reader of popnei gives blocks of one variant at least, and
        // a block of none would leave the sets with no word to write a bit
        // in.
        if block.num_vars == 0 {
            return Err(Error::ReaderGaveABlockOfNoVariants);
        }
        if !block.fields().contains(Needs::GTS) {
            return Err(Error::FieldsNotInTheBlock { fields: Needs::GTS });
        }
        // The genotypes are the variants times the individuals times the
        // ploidy, so that the rows below are the rows of the block.
        block.check()?;

        // n is at most the variants of the block, and the ploidy times the
        // sum of d at most the ploidy times n, so the two sums of a pair
        // over one block fit in a `u32` when the variants times the ploidy
        // do.
        let sums_too_large = || Error::KosmanSumsTooLarge {
            num_vars: u64::try_from(block.num_vars).unwrap_or(u64::MAX),
            ploidy: block.ploidy,
        };
        let num_vars = u32::try_from(block.num_vars).map_err(|_| sums_too_large())?;
        let ploidy = u32::try_from(block.ploidy).map_err(|_| sums_too_large())?;
        num_vars.checked_mul(ploidy).ok_or_else(sums_too_large)?;

        let too_large_for_the_memory = || Error::BlockTooLarge {
            num_vars_per_block: block.num_vars,
            num_individuals: block.num_individuals,
            ploidy: block.ploidy,
            size: BlockSize::AskedFor,
        };
        let alleles_per_var = block
            .num_individuals
            .checked_mul(block.ploidy)
            .ok_or_else(too_large_for_the_memory)?;
        // A block of no individual holds no genotype, so the fields above
        // gave the error already: what this keeps is the rows of the loop
        // below from being rows of no allele.
        if alleles_per_var == 0 {
            return Err(Error::FieldsNotInTheBlock { fields: Needs::GTS });
        }

        let (smallest_allele, largest_allele) = smallest_and_largest_allele(&block.gts);
        if smallest_allele < MISSING_ALLELE {
            return Err(Error::AlleleBelowTheMissingOne {
                allele: smallest_allele,
            });
        }
        // The alleles of the block are one more than its largest allele.
        // A block whose genotypes are all missing holds none, and it gets
        // the sets of one allele, which no genotype sets a bit in: what
        // that keeps is the `holds` sets of an individual from being no
        // word at all.
        let num_alleles = usize::try_from(largest_allele)
            .map_or(0, |allele| allele.saturating_add(1))
            .max(1);

        // A block of one variant at least has one word in a set at least.
        let words_per_set = NonZeroUsize::new(block.num_vars.div_ceil(VARS_PER_WORD))
            .ok_or(Error::ReaderGaveABlockOfNoVariants)?;
        let holds_per_individual = num_alleles
            .checked_mul(block.ploidy)
            .and_then(|sets| sets.checked_mul(words_per_set.get()))
            .and_then(NonZeroUsize::new)
            .ok_or_else(too_large_for_the_memory)?;
        let mut called = words_of(
            block.num_individuals,
            words_per_set,
            &too_large_for_the_memory,
        )?;
        let mut holds = words_of(
            block.num_individuals,
            holds_per_individual,
            &too_large_for_the_memory,
        )?;

        for (var, row) in block.gts.chunks_exact(alleles_per_var).enumerate() {
            let word_of_the_var = var / VARS_PER_WORD;
            let bit_of_the_var = 1_u64 << (var % VARS_PER_WORD);
            let individuals = called
                .chunks_exact_mut(words_per_set.get())
                .zip(holds.chunks_exact_mut(holds_per_individual.get()));
            for ((called, holds), genotype) in individuals.zip(row.chunks_exact(block.ploidy)) {
                if genotype.contains(&MISSING_ALLELE) {
                    continue;
                }
                // The word of the variant is below the words of a set, and
                // a `called` set is the words of one set, so `get_mut`
                // gives that word.
                if let Some(word) = called.get_mut(word_of_the_var) {
                    *word |= bit_of_the_var;
                }
                for (place, &allele) in genotype.iter().enumerate() {
                    // The copies of this allele up to this place of the
                    // genotype: the m-th time the allele is met is the set
                    // of the genotypes that hold m copies of it or more.
                    let copies = genotype
                        .iter()
                        .take(place + 1)
                        .filter(|&&other| other == allele)
                        .count();
                    // The allele is 0 or more: a genotype with an allele
                    // that was not called was left above, and an allele
                    // below the missing one was refused before the loop.
                    let Ok(allele) = usize::try_from(allele) else {
                        continue;
                    };
                    // The set of the allele and the count is below k * A,
                    // since the allele is below A and the copies are the
                    // ploidy at most, so its word is below the words of
                    // the `holds` sets of one individual.
                    let set = allele * block.ploidy + copies - 1;
                    if let Some(word) = holds.get_mut(set * words_per_set.get() + word_of_the_var) {
                        *word |= bit_of_the_var;
                    }
                }
            }
        }

        Ok(KosmanBits {
            num_individuals: block.num_individuals,
            ploidy,
            words_per_set,
            holds_per_individual,
            called,
            holds,
        })
    }

    /// How many individuals the sets are of.
    pub(crate) fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// The ploidy times the sum of d over the block, and n, the variants
    /// of the block at which both genotypes were called, for the pair of
    /// the individuals `first` and `second`.
    ///
    /// The two numbers are the same for (i, j) and for (j, i). It gives
    /// `None` when the two are one individual and when either of them is
    /// not an individual of the block.
    pub(crate) fn sums_of_the_pair(&self, first: usize, second: usize) -> Option<(u32, u32)> {
        if first == second {
            return None;
        }
        let first = self.individuals_from(first).next()?;
        let second = self.individuals_from(second).next()?;
        Some(sums_of_two(first, second, self.ploidy))
    }

    /// The two sums of every pair of individuals of the block, in the
    /// order of the distance vector: (0, 1), (0, 2), ..., (1, 2), ...
    pub(crate) fn sums_of_the_pairs(&self) -> impl Iterator<Item = (u32, u32)> + '_ {
        self.individuals_from(0)
            .enumerate()
            .flat_map(move |(first, one)| {
                self.individuals_from(first.saturating_add(1))
                    .map(move |other| sums_of_two(one, other, self.ploidy))
            })
    }

    /// The `called` set and the `holds` sets of each individual from the
    /// one at `first` on, in the order of the individuals of the block.
    ///
    /// The individuals before it are left out of each of the two arrays
    /// and not walked over, so the pairs of one individual are reached in
    /// the same time whichever individual it is.
    fn individuals_from(&self, first: usize) -> impl Iterator<Item = (&[u64], &[u64])> {
        self.called
            .chunks_exact(self.words_per_set.get())
            .skip(first)
            .zip(
                self.holds
                    .chunks_exact(self.holds_per_individual.get())
                    .skip(first),
            )
    }
}

/// `num_individuals` times `words_per_individual` words at 0, or the error
/// of a block whose sets the machine does not give the memory for.
///
/// The memory is asked for with `try_reserve_exact`, which gives it back as
/// an error where `vec![0; n]` would end the process.
fn words_of(
    num_individuals: usize,
    words_per_individual: NonZeroUsize,
    too_large: &impl Fn() -> Error,
) -> Result<Vec<u64>> {
    let num_words = num_individuals
        .checked_mul(words_per_individual.get())
        .ok_or_else(too_large)?;
    let mut words: Vec<u64> = Vec::new();
    words
        .try_reserve_exact(num_words)
        .map_err(|_| too_large())?;
    words.resize(num_words, 0);
    Ok(words)
}

/// The smallest and the largest allele of the genotypes, both
/// [`MISSING_ALLELE`] when they are all missing or there are none.
///
/// The largest says how many alleles the block holds, and so how many sets
/// each individual gets; the smallest is looked for in the same pass,
/// because an allele below the missing one, which no reader of popnei
/// gives, has no set of its own and would be counted as an allele that was
/// called.
fn smallest_and_largest_allele(gts: &[i8]) -> (i8, i8) {
    let mut smallest = MISSING_ALLELE;
    let mut largest = MISSING_ALLELE;
    for &allele in gts {
        smallest = smallest.min(allele);
        largest = largest.max(allele);
    }
    (smallest, largest)
}

/// The ploidy times the sum of d, and n, of the two individuals whose sets
/// of bits these are.
///
/// Each of the two is the `called` set and the `holds` sets of one
/// individual of one block, so the `called` sets are of one length and the
/// `holds` sets of another.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the sets are of one block, whose variants times the ploidy `of_block` checked \
              to be a number a u32 holds: the variants called in both are at most the \
              variants of the block, and the alleles that pair at most the ploidy times \
              that, since at a variant called in both the alleles of one that pair with an \
              equal allele of the other are at most the ploidy"
)]
fn sums_of_two(first: (&[u64], &[u64]), second: (&[u64], &[u64]), ploidy: u32) -> (u32, u32) {
    let called_in_both: u32 = first
        .0
        .iter()
        .zip(second.0)
        .map(|(one, other)| (one & other).count_ones())
        .sum();
    let alleles_that_pair: u32 = first
        .1
        .iter()
        .zip(second.1)
        .map(|(one, other)| (one & other).count_ones())
        .sum();
    (ploidy * called_in_both - alleles_that_pair, called_in_both)
}

#[cfg(test)]
mod tests {
    use super::KosmanBits;
    use crate::block::Block;
    use crate::error::Error;
    use crate::variant::MISSING_ALLELE;

    /// A block of the variants given, of `num_individuals` individuals of
    /// the ploidy `ploidy`, with the genotypes and no column: the fields a
    /// calculation of distances asks its reader for. Each row is one
    /// variant, the alleles of one individual after those of the
    /// individual before it.
    fn block_of(variants: &[&[i8]], num_individuals: usize, ploidy: usize) -> Block {
        let mut gts = Vec::new();
        for row in variants {
            gts.extend_from_slice(row);
        }
        Block {
            num_vars: variants.len(),
            num_individuals,
            ploidy,
            gts,
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        }
    }

    /// The four variants of three diploid individuals of the worked
    /// example of "How it is verified" of `docs/specs/dists.md`: 0/0 0/1
    /// 1/1, 0/1 0/1 1/2, 0/0 0/. 2/2 and ./. 1/1 1/1. The third variant
    /// has the half called genotype, which is missing, and the second and
    /// the third hold three alleles where the first and the last hold two.
    fn the_diploid_worked_example() -> Block {
        const M: i8 = MISSING_ALLELE;
        block_of(
            &[
                &[0, 0, 0, 1, 1, 1],
                &[0, 1, 0, 1, 1, 2],
                &[0, 0, 0, M, 2, 2],
                &[M, M, 1, 1, 1, 1],
            ],
            3,
            2,
        )
    }

    /// The three variants of three tetraploid individuals of the same
    /// part of the spec: 0/0/0/1 0/1/1/1 1/1/1/1, 0/0/1/1 0/1/0/1 0/0/2/2
    /// and 0/0/0/0 0/0/./0 1/2/2/2.
    fn the_tetraploid_worked_example() -> Block {
        const M: i8 = MISSING_ALLELE;
        block_of(
            &[
                &[0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1],
                &[0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 2, 2],
                &[0, 0, 0, 0, 0, 0, M, 0, 1, 2, 2, 2],
            ],
            3,
            4,
        )
    }

    /// The four variants of three haploid individuals of the same part of
    /// the spec: 0 0 1, 0 1 2, . 1 1 and 0 0 0.
    fn the_haploid_worked_example() -> Block {
        const M: i8 = MISSING_ALLELE;
        block_of(&[&[0, 0, 1], &[0, 1, 2], &[M, 1, 1], &[0, 0, 0]], 3, 1)
    }

    /// The numbers of the table of the diploid worked example of "How it
    /// is verified" of `docs/specs/dists.md`: twice the sum of d and n of
    /// each of the three pairs. pyNei and `gd.kosman` of R give the same
    /// distances, 0.25, 0.833333 and 0.333333.
    #[test]
    fn the_diploid_worked_example_gives_the_sums_of_the_spec() {
        let bits = KosmanBits::of_block(&the_diploid_worked_example()).unwrap();

        assert_eq!(bits.sums_of_the_pair(0, 1), Some((1, 2)));
        assert_eq!(bits.sums_of_the_pair(0, 2), Some((5, 3)));
        assert_eq!(bits.sums_of_the_pair(1, 2), Some((2, 3)));
    }

    /// The numbers of the table of the tetraploid worked example, which
    /// `gd.kosman` of R gives: four times the sum of d and n. The ploidy
    /// is read here and not taken for 2, and a genotype of four alleles
    /// holds one allele twice and more.
    #[test]
    fn the_tetraploid_worked_example_gives_the_sums_of_the_spec() {
        let bits = KosmanBits::of_block(&the_tetraploid_worked_example()).unwrap();

        assert_eq!(bits.sums_of_the_pair(0, 1), Some((2, 2)));
        assert_eq!(bits.sums_of_the_pair(0, 2), Some((9, 3)));
        assert_eq!(bits.sums_of_the_pair(1, 2), Some((3, 2)));
    }

    /// The numbers of the table of the haploid worked example, which
    /// `gd.kosman` of R gives: the ploidy is 1, so the first number is the
    /// sum of d itself, 0 for two equal alleles and 1 for two different
    /// ones.
    #[test]
    fn the_haploid_worked_example_gives_the_sums_of_the_spec() {
        let bits = KosmanBits::of_block(&the_haploid_worked_example()).unwrap();

        assert_eq!(bits.sums_of_the_pair(0, 1), Some((1, 3)));
        assert_eq!(bits.sums_of_the_pair(0, 2), Some((2, 3)));
        assert_eq!(bits.sums_of_the_pair(1, 2), Some((2, 4)));
    }

    /// The pairs come in the order of the distance vector, (0, 1), (0, 2),
    /// (1, 2), which is the order the three numbers of the spec's table
    /// are in.
    #[test]
    fn the_pairs_of_a_block_come_in_the_order_of_the_distance_vector() {
        let bits = KosmanBits::of_block(&the_diploid_worked_example()).unwrap();

        assert_eq!(
            bits.sums_of_the_pairs().collect::<Vec<_>>(),
            [(1, 2), (5, 3), (2, 3)]
        );
    }

    /// The distance of a pair does not depend on which of the two
    /// individuals is named first.
    #[test]
    fn the_sums_of_a_pair_are_the_same_in_either_order() {
        let bits = KosmanBits::of_block(&the_diploid_worked_example()).unwrap();

        assert_eq!(bits.sums_of_the_pair(2, 0), bits.sums_of_the_pair(0, 2));
        assert_eq!(bits.sums_of_the_pair(2, 1), bits.sums_of_the_pair(1, 2));
    }

    /// A genotype with one allele that was not called is missing, so its
    /// variant counts for no pair the individual is in, and the called
    /// allele is compared with nothing. The third variant of the diploid
    /// worked example is one, and here it is the only variant.
    #[test]
    fn a_half_called_genotype_is_missing_and_its_variant_counts_for_no_pair() {
        const M: i8 = MISSING_ALLELE;
        let block = block_of(&[&[0, M, 0, 0]], 2, 2);

        let bits = KosmanBits::of_block(&block).unwrap();

        assert_eq!(bits.sums_of_the_pair(0, 1), Some((0, 0)));
    }

    /// Every pair of a block whose genotypes are all missing has no
    /// variant called in both, so n is 0 and the sum of d is 0: the pair
    /// gets its distance, if it gets one, from the other blocks.
    #[test]
    fn a_block_whose_genotypes_are_all_missing_gives_no_called_variant_to_any_pair() {
        const M: i8 = MISSING_ALLELE;
        let block = block_of(&[&[M; 6], &[M; 6]], 3, 2);

        let bits = KosmanBits::of_block(&block).unwrap();

        assert_eq!(
            bits.sums_of_the_pairs().collect::<Vec<_>>(),
            [(0, 0), (0, 0), (0, 0)]
        );
    }

    /// A block of one individual has no pair: the sums of every pair of
    /// its individuals are none, and one individual with itself is not a
    /// pair.
    #[test]
    fn a_block_of_one_individual_has_no_pair() {
        let block = block_of(&[&[0, 1], &[0, 0]], 1, 2);

        let bits = KosmanBits::of_block(&block).unwrap();

        assert_eq!(bits.num_individuals(), 1);
        assert_eq!(bits.sums_of_the_pairs().count(), 0);
        assert_eq!(bits.sums_of_the_pair(0, 0), None);
        assert_eq!(bits.sums_of_the_pair(0, 1), None);
    }

    /// A set of more than one word: 100 variants are 64 bits in the first
    /// word and 36 in the second, and the 28 bits of the second word that
    /// are no variant are set for nobody and counted for nobody. The first
    /// individual is 0/0 at every variant, the second 0/1 at every variant,
    /// and the third 1/1 at the even ones and missing at the odd ones, so
    /// the pairs have n of 100, 50 and 50 and d of 0.5, 1 and 0.5 at each
    /// of those variants.
    #[test]
    fn a_block_of_more_than_64_variants_counts_the_variants_of_every_word() {
        const M: i8 = MISSING_ALLELE;
        let rows: Vec<Vec<i8>> = (0..100)
            .map(|var| {
                let third = if var % 2 == 0 { [1, 1] } else { [M, M] };
                vec![0, 0, 0, 1, third[0], third[1]]
            })
            .collect();
        let variants: Vec<&[i8]> = rows.iter().map(Vec::as_slice).collect();

        let bits = KosmanBits::of_block(&block_of(&variants, 3, 2)).unwrap();

        assert_eq!(bits.sums_of_the_pair(0, 1), Some((100, 100)));
        assert_eq!(bits.sums_of_the_pair(0, 2), Some((100, 50)));
        assert_eq!(bits.sums_of_the_pair(1, 2), Some((50, 50)));
    }

    /// Every allele of a multiallelic variant counts as itself, so a
    /// variant of the alleles 0 and 5 beside two variants of the alleles 0
    /// and 1 has the sets of the allele 5 too: 5/5 and 0/5 share one copy
    /// of the 5 and have d of 0.5, as 0/1 and 0/0 do.
    #[test]
    fn a_multiallelic_variant_beside_biallelic_ones_counts_each_of_its_alleles() {
        let block = block_of(&[&[0, 1, 0, 0], &[5, 5, 0, 5], &[0, 1, 1, 0]], 2, 2);

        let bits = KosmanBits::of_block(&block).unwrap();

        assert_eq!(bits.sums_of_the_pair(0, 1), Some((2, 3)));
    }

    /// The genotypes of a block are its variants times its individuals
    /// times its ploidy, and a block whose genotypes are fewer would be
    /// read one individual at the place of another. Only a reader with a
    /// defect builds one.
    #[test]
    fn a_block_whose_genotypes_are_not_as_many_as_its_size_says_is_an_error() {
        let mut block = the_diploid_worked_example();
        block.gts.pop();

        let error = KosmanBits::of_block(&block).unwrap_err();

        assert!(
            matches!(
                error,
                Error::BlockArrayOfAnotherSize {
                    array: "gts",
                    found: 23,
                    expected: 24
                }
            ),
            "{error}"
        );
    }

    /// A block with no genotypes is the error of a field that is not
    /// there: the calculation asks its reader for the genotypes and reads
    /// nothing else.
    #[test]
    fn a_block_with_no_genotypes_is_the_error_of_a_field_that_is_not_there() {
        let mut block = the_diploid_worked_example();
        block.gts = Vec::new();

        let error = KosmanBits::of_block(&block).unwrap_err();

        assert!(
            matches!(error, Error::FieldsNotInTheBlock { fields } if fields == crate::variant::Needs::GTS),
            "{error}"
        );
    }
}
