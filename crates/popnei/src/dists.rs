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
//! [`calc_kosman_sums`] is the calculation a user reaches: it reads a
//! reader of blocks to its end and gives the [`KosmanSums`], the two
//! integers of every pair over all of its variants, which the distance of
//! each pair is worked out from. [`KosmanBits`] is the genotypes of one
//! block as sets of bits, one set for each individual and each question
//! that a count of a pair asks, and it gives those two integers for a pair
//! over that block.
//!
//! `docs/specs/dists.md` has the design, the formula and the numbers the
//! tests assert, and the row `dists` of section 9 of
//! `docs/architecture.md` is where the module sits.

use std::num::NonZeroUsize;

use crate::block::{Block, BlockReader, BlockSize};
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
    pub(crate) fn of_block(block: &Block) -> Result<KosmanBits> {
        KosmanBits::of_block_built(block, HowTheSetsAreBuilt::OnTheThreads)
    }

    /// The same sets, with the bits of one individual written after those
    /// of the individual before it, which is what the test that compares
    /// the two ways of building them calls.
    ///
    /// # Errors
    ///
    /// The same as [`KosmanBits::of_block`].
    #[cfg(test)]
    pub(crate) fn of_block_one_by_one(block: &Block) -> Result<KosmanBits> {
        KosmanBits::of_block_built(block, HowTheSetsAreBuilt::OneAfterAnother)
    }

    /// The sets of bits of the genotypes of `block`, written the way `how`
    /// says.
    ///
    /// # Errors
    ///
    /// The same as [`KosmanBits::of_block`].
    fn of_block_built(block: &Block, how: HowTheSetsAreBuilt) -> Result<KosmanBits> {
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

        let alleles = alleles_of(&block.gts);
        if alleles.smallest < MISSING_ALLELE {
            return Err(Error::AlleleBelowTheMissingOne {
                allele: alleles.smallest,
            });
        }
        let num_alleles = alleles.num_alleles;

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

        let shape = ShapeOfTheSets {
            ploidy: block.ploidy,
            words_per_set,
            holds_per_individual,
            place: &alleles.place,
        };
        match how {
            HowTheSetsAreBuilt::OnTheThreads => {
                write_the_sets(&block.gts, alleles_per_var, shape, &mut called, &mut holds);
            }
            HowTheSetsAreBuilt::OneAfterAnother => {
                write_the_sets_one_by_one(
                    &block.gts,
                    alleles_per_var,
                    shape,
                    &mut called,
                    &mut holds,
                );
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
}

/// Whether the sets of a block are built on the threads of rayon or one
/// individual after another.
///
/// In wasm, which has no threads, both write the same bits one individual
/// after another.
#[derive(Clone, Copy)]
enum HowTheSetsAreBuilt {
    /// The individuals shared out over the threads of the pool the caller
    /// is running in.
    OnTheThreads,
    /// One individual after another, on the thread that asked.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "outside the tests nothing asks for this: `of_block` asks for the \
                      threads, and in wasm, which has none, that is already the same pass \
                      one individual after another"
        )
    )]
    OneAfterAnother,
}

/// How many individuals the sets of one work item are built for, when the
/// sets of a block are built on the threads.
///
/// An item reads, out of every row of the block, the alleles of its own
/// individuals and no others, this many times the ploidy bytes of each
/// row, and writes the sets of those individuals and no others. The
/// smaller it is, the more often the items read the same 128 byte line of
/// the genotypes, 64 diploid individuals being exactly one such line; the
/// larger it is, the fewer items there are for the threads to share, 1000
/// individuals giving 125 items at 8 and 16 at 64. What it also sets is
/// how much of the sets an item writes into at a time: at 64 the words of
/// the item are 202 KB of the 3.16 MB of a block of 1000 individuals.
///
/// It was swept over 8, 14, 32 and 64 on 100000 variants of 1000 diploid
/// individuals, biallelic, 3 in 100 genotypes missing, the blocks handed
/// out from memory, on the owner's Apple M5 Pro. The sets phase took
/// 0.062 s, 0.061 s, 0.058 s and 0.054 s on 18 threads and 0.230 s,
/// 0.192 s, 0.183 s and 0.172 s on one, so the largest of the four won on
/// both and nothing in the sweep says where it stops winning. A dataset of
/// fewer than about 1000 individuals gets fewer items than the machine has
/// threads, and 500 individuals of 5000 variants is the second shape
/// `docs/reports/perf-dists-kosman-2026-09-22.md` asks the benchmark for.
const INDIVIDUALS_PER_ITEM: usize = 64;

/// What the writing of a row needs to know besides the row itself and the
/// sets it writes into: the ploidy of the block, the words of one set and
/// of the `holds` sets of one individual, and where each allele value of
/// the block sits among its alleles.
#[derive(Clone, Copy)]
struct ShapeOfTheSets<'a> {
    /// How many alleles the genotype of one individual holds.
    ploidy: usize,
    /// How many `u64` one set holds.
    words_per_set: NonZeroUsize,
    /// How many `u64` the `holds` sets of one individual hold.
    holds_per_individual: NonZeroUsize,
    /// Where each allele value sits among the alleles of the block, the
    /// `place` of [`AllelesOfTheBlock`].
    place: &'a [u8; ALLELE_VALUES],
}

/// The sets of bits of the individuals of a block, written from its
/// genotypes.
///
/// `gts` holds the rows of the block, `alleles_per_var` alleles each, one
/// individual's after the one before it. `called` holds the `called` set of
/// every individual and `holds` their `holds` sets, an individual's side by
/// side in each, and both come in zeroed.
///
/// Natively the individuals are shared out over the threads of rayon,
/// [`INDIVIDUALS_PER_ITEM`] of them to a work item: the bits of an
/// individual lie in its own words of the two arrays, so no two items write
/// the same word and none of them reads another's, and the bits are the
/// same however many threads there are. Inside an item the rows are read in
/// the order of the variants, so that `gts` is still read from its start to
/// its end. The threads are those of the pool the caller is running in, as
/// in `parse_rows` of the VCF reader. In wasm there are no threads and the
/// same bits are written one individual after another.
#[cfg(not(target_family = "wasm"))]
fn write_the_sets(
    gts: &[i8],
    alleles_per_var: usize,
    shape: ShapeOfTheSets<'_>,
    called: &mut [u64],
    holds: &mut [u64],
) {
    use rayon::iter::{IndexedParallelIterator, ParallelIterator};
    use rayon::slice::ParallelSliceMut;

    // The three sizes are of the same run of individuals, so either all
    // three are there or the sets are written one after another: two of
    // them cut into runs of a different length would pair the words of one
    // individual with the words of another. A block large enough to
    // overflow any of them is a block whose sets the machine did not give
    // the memory of, so this is the arithmetic and not a case a caller
    // meets.
    let sizes = shape
        .words_per_set
        .get()
        .checked_mul(INDIVIDUALS_PER_ITEM)
        .zip(
            shape
                .holds_per_individual
                .get()
                .checked_mul(INDIVIDUALS_PER_ITEM),
        )
        .zip(shape.ploidy.checked_mul(INDIVIDUALS_PER_ITEM));
    let Some(((called_per_item, holds_per_item), alleles_per_item)) = sizes else {
        write_the_sets_one_by_one(gts, alleles_per_var, shape, called, holds);
        return;
    };
    called
        .par_chunks_mut(called_per_item)
        .zip(holds.par_chunks_mut(holds_per_item))
        .enumerate()
        .for_each(|(item, (called_of_the_item, holds_of_the_item))| {
            for (var, row) in gts.chunks_exact(alleles_per_var).enumerate() {
                // The rows are all of `alleles_per_var` alleles and every
                // item takes the same run of every row, so the item that
                // has a run of the words has a run of the alleles.
                if let Some(alleles_of_the_item) = row.chunks(alleles_per_item).nth(item) {
                    write_the_row(
                        alleles_of_the_item,
                        var,
                        shape,
                        called_of_the_item,
                        holds_of_the_item,
                    );
                }
            }
        });
}

/// The same bits, written one individual after another, which is what wasm
/// does: it has no threads.
#[cfg(target_family = "wasm")]
fn write_the_sets(
    gts: &[i8],
    alleles_per_var: usize,
    shape: ShapeOfTheSets<'_>,
    called: &mut [u64],
    holds: &mut [u64],
) {
    write_the_sets_one_by_one(gts, alleles_per_var, shape, called, holds);
}

/// The bits of every individual written one after another, over the rows of
/// the block in the order of the variants.
///
/// It is compiled for every target and not for wasm alone, so that the
/// cargo tests, which run natively, can build the sets of the same block
/// with it and with the threads and compare the bits.
fn write_the_sets_one_by_one(
    gts: &[i8],
    alleles_per_var: usize,
    shape: ShapeOfTheSets<'_>,
    called: &mut [u64],
    holds: &mut [u64],
) {
    for (var, row) in gts.chunks_exact(alleles_per_var).enumerate() {
        write_the_row(row, var, shape, called, holds);
    }
}

/// The bits of one row of the genotypes of a block, for the run of
/// individuals whose sets `called` and `holds` are.
///
/// `alleles` is the alleles of those individuals at the variant `var` of
/// the block, the ploidy of each after the one before it, and `called` and
/// `holds` hold their sets in the same order. A genotype that holds the
/// missing allele sets no bit, so its variant counts as called for no pair
/// the individual is in.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the word of a variant and its bit are a division and a remainder by 64, \
              which is not 0, and a shift by a remainder of 64, which is below the bits \
              of a u64; `at + 1`, how far into a genotype an allele is, is at most the \
              ploidy, which `Block::check` made at most the genotypes of the block, their \
              number being the variants times the individuals times the ploidy; and the \
              set of an allele and a count, `place * ploidy + copies - 1`, is 0 at least, \
              since the copies are 1 at least, and below `A * k`, since the place of the \
              allele among the A alleles of the block is below A and the copies are at \
              most the ploidy, so its word is below the words of the `holds` sets of one \
              individual, which `of_block_built` checked to be a number a usize holds"
)]
fn write_the_row(
    alleles: &[i8],
    var: usize,
    shape: ShapeOfTheSets<'_>,
    called: &mut [u64],
    holds: &mut [u64],
) {
    let word_of_the_var = var / VARS_PER_WORD;
    let bit_of_the_var = 1_u64 << (var % VARS_PER_WORD);
    let individuals = called
        .chunks_exact_mut(shape.words_per_set.get())
        .zip(holds.chunks_exact_mut(shape.holds_per_individual.get()));
    for ((called, holds), genotype) in individuals.zip(alleles.chunks_exact(shape.ploidy)) {
        #[expect(
            clippy::manual_contains,
            reason = "`contains` on a slice of `i8` goes to `memchr`, which searches \
                              a word at a time and pays its setup for a genotype of the \
                              ploidy, 2 bytes here; a profile of 100000 variants of 1000 \
                              diploid individuals had it at 9.4 in 100 of the CPU on one \
                              thread, and `any` inlines to one compare per allele"
        )]
        if genotype.iter().any(|&allele| allele == MISSING_ALLELE) {
            continue;
        }
        // The word of the variant is below the words of a set, and
        // a `called` set is the words of one set, so `get_mut`
        // gives that word.
        if let Some(word) = called.get_mut(word_of_the_var) {
            *word |= bit_of_the_var;
        }
        for (at, &allele) in genotype.iter().enumerate() {
            // The copies of this allele up to this place of the
            // genotype: the m-th time the allele is met is the set
            // of the genotypes that hold m copies of it or more.
            //
            // A genotype of two alleles answers that with one
            // compare: the first allele of the genotype has one
            // copy of itself so far, and the second has two when it
            // equals the first and one otherwise. The general count
            // below is a loop of an unknown length, which the
            // compiler lowers to a 32 byte vector loop, an 8 byte
            // one and a scalar tail, for the one or two alleles a
            // diploid genotype gives it.
            let copies = if shape.ploidy == 2 {
                if at == 1 && genotype.first().is_some_and(|&first| first == allele) {
                    2
                } else {
                    1
                }
            } else {
                genotype
                    .iter()
                    .take(at + 1)
                    .filter(|&&other| other == allele)
                    .count()
            };
            // Every allele that reaches here is 0 to `MAX_ALLELE`,
            // which is 127, so it is one of the 128 values of the
            // table of places and the lookup cannot fail: a
            // genotype that holds the missing allele was skipped
            // above, and a block whose smallest allele is below the
            // missing one was refused before the loop, so no allele
            // of this loop is negative and `cast_unsigned` gives
            // its own value. The `unwrap_or` is what a table of 128
            // places costs against one of 256, and it is never
            // taken; were it taken, the allele would be counted as
            // the first allele of the block.
            let place = usize::from(
                shape
                    .place
                    .get(usize::from(allele.cast_unsigned()))
                    .copied()
                    .unwrap_or(0),
            );
            // The set of the allele and the count is below k * A,
            // since the place of the allele among the alleles of
            // the block is below A and the copies are the ploidy at
            // most, so its word is below the words of the `holds`
            // sets of one individual.
            let set = place * shape.ploidy + copies - 1;
            let word_of_the_set = set * shape.words_per_set.get() + word_of_the_var;
            if let Some(word) = holds.get_mut(word_of_the_set) {
                *word |= bit_of_the_var;
            }
        }
    }
}

impl KosmanBits {
    /// The ploidy times the sum of d over the block, and n, the variants
    /// of the block at which both genotypes were called, for the pair of
    /// the individuals `first` and `second`.
    ///
    /// The two numbers are the same for (i, j) and for (j, i). It gives
    /// `None` when the two are one individual and when either of them is
    /// not an individual of the block.
    ///
    /// The calculation over a reader takes the pairs of a block a row of
    /// the upper triangle at a time and never one pair at a time, so this
    /// is compiled for the tests that assert the two counts of one named
    /// pair and for nothing else.
    #[cfg(test)]
    pub(crate) fn sums_of_one_pair(&self, first: usize, second: usize) -> Option<(u32, u32)> {
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
        (0..self.num_individuals).flat_map(move |first| self.sums_of_the_pairs_of(first))
    }

    /// The two sums of the pairs of the individual `first` with every
    /// individual after it, in the order of the distance vector: for the
    /// individual 1 of a block of four, (1, 2) and then (1, 3).
    ///
    /// These are the pairs of one row of the upper triangle of the square
    /// matrix of the distances, and they are the run of pairs that a thread
    /// takes. It gives nothing when `first` is not an individual of the
    /// block.
    pub(crate) fn sums_of_the_pairs_of(&self, first: usize) -> impl Iterator<Item = (u32, u32)> {
        self.individuals_from(first)
            .next()
            .into_iter()
            .flat_map(move |one| {
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

/// How many allele values a genotype of popnei can hold: 0 to
/// [`MAX_ALLELE`](crate::variant::MAX_ALLELE), which is 127.
const ALLELE_VALUES: usize = 128;

/// The alleles the genotypes of one block hold, each with its place among
/// them, and the smallest allele those genotypes hold.
struct AllelesOfTheBlock {
    /// Where each allele value from 0 to
    /// [`MAX_ALLELE`](crate::variant::MAX_ALLELE) sits among the
    /// alleles the block holds, counting from 0 in the order of the values:
    /// in a block whose genotypes hold the alleles 0 and 127, the 0 is at 0
    /// and the 127 at 1. A value no genotype holds is at 0, which no
    /// genotype of the block looks up. A place is below 128, the values of
    /// the table, so each of them is a number a `u8` holds and the table is
    /// 128 bytes.
    place: [u8; ALLELE_VALUES],
    /// How many allele values the genotypes hold, and so how many sets of
    /// each count each individual gets. It is 1 at least, also for a block
    /// whose genotypes are all missing, so that the `holds` sets of an
    /// individual are one word at least.
    num_alleles: usize,
    /// The smallest allele of the genotypes, [`MISSING_ALLELE`] when they
    /// are all missing or there are none.
    smallest: i8,
}

/// The alleles of the genotypes of a block, in one pass over them.
///
/// The sets of an individual cost in proportion to the alleles the block
/// holds and not to the largest of their values, so the values are ranked:
/// a block of 3000 variants of 1000 individuals with one variant of the
/// alleles 0 and 127 would otherwise give every pair 10240 words to walk in
/// the place of 160, and took 7.74 s where the same block with the allele 1
/// in the place of the 127 takes 0.75 s.
///
/// The smallest allele is looked for in the same pass, because an allele
/// below the missing one, which no reader of popnei gives, has no place
/// among the alleles and would be counted as an allele that was called.
fn alleles_of(gts: &[i8]) -> AllelesOfTheBlock {
    let mut held = [false; ALLELE_VALUES];
    let mut smallest = MISSING_ALLELE;
    for &allele in gts {
        smallest = smallest.min(allele);
        // An allele that was called is 0 to `MAX_ALLELE`, so it is one of
        // the values of the table; the missing one and anything below it
        // are not, and the smallest above is what reports them.
        if let Ok(value) = usize::try_from(allele)
            && let Some(held) = held.get_mut(value)
        {
            *held = true;
        }
    }
    let mut place = [0_u8; ALLELE_VALUES];
    let mut num_alleles: usize = 0;
    for (value, held) in held.iter().enumerate() {
        if *held {
            if let Some(place) = place.get_mut(value) {
                // The places are as many as the values of the table at
                // most, which is 128, so this place is 127 at most and the
                // conversion never takes its 0.
                *place = u8::try_from(num_alleles).unwrap_or(0);
            }
            // The places are as many as the values of the table at most,
            // which is 128.
            num_alleles = num_alleles.saturating_add(1);
        }
    }
    AllelesOfTheBlock {
        place,
        num_alleles: num_alleles.max(1),
        smallest,
    }
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

/// The two sums of every pair of individuals over the variants of
/// `reader`, from which the Kosman distance of each pair is worked out.
///
/// It asks `reader` for the genotypes alone and reads it to its end,
/// building the sets of bits of each block and adding the two integers of
/// every pair of that block to what the blocks before it gave. The sums are
/// integers, so the result is the same to the last bit whatever the size of
/// the blocks and however many threads the pairs were counted on, and the
/// reader needs no [`Reblock`](crate::block::Reblock) before it.
///
/// `reader` is borrowed and not taken, so that whoever built the chain of
/// filters of the pass can read their counts from it when this returns, as
/// `docs/specs/filters.md` says; how many variants the calculation was
/// given is [`KosmanSums::num_vars`].
///
/// # Errors
///
/// When `reader` has no variant, when the sums of a pair would go above
/// what a `u32` holds, when the machine does not give the memory of the two
/// counts of every pair, and when the reader fails, whose error is given on
/// as it is. A block of other individuals or of another ploidy than the
/// reader says its source has is the error of a reader with a defect.
pub fn calc_kosman_sums<R: BlockReader + ?Sized>(reader: &mut R) -> Result<KosmanSums> {
    // The genotypes are all this reads, so a reader over a file leaves the
    // columns of a variant unparsed.
    reader.set_needs(Needs::GTS);
    let num_individuals = reader.individuals().len();
    let ploidy = reader.ploidy();
    // The two counts of every pair are asked of the machine once, when the
    // first block is there: at 10000 individuals they are 400 MB, which a
    // reader with no variant would have asked for and given back.
    let Some(block) = reader.next_block()? else {
        return Err(Error::ReaderGaveNoVariants);
    };
    let mut sums = KosmanSums::at_zero(num_individuals, ploidy)?;
    let mut next = Some(block);
    while let Some(block) = next.take() {
        // Every variant of the block counts here, called in a pair or not:
        // `num_vars` is what a user reads as the variants of the pass. A
        // `usize` is 64 bits natively and 32 in wasm, so the conversion
        // holds; the sum stops at the largest `u64`, which is more variants
        // than any source holds, and the sums of a pair are refused long
        // before it.
        sums.num_vars = sums
            .num_vars
            .saturating_add(u64::try_from(block.num_vars).unwrap_or(u64::MAX));
        let bits = KosmanBits::of_block(&block)?;
        add_the_block(&mut sums, &bits)?;
        // The sets of the block and the block itself are given back before
        // the reader is asked for the next one, so that the memory of two
        // blocks and of two sets of bits is never held at once.
        drop(bits);
        drop(block);
        next = reader.next_block()?;
    }
    Ok(sums)
}

/// For every pair of individuals, the ploidy times the sum of d over the
/// variants of a pass and n, how many of those variants both individuals
/// were called at. The Kosman distance of the pair is the first over the
/// ploidy times the second.
///
/// The pairs are in the order of the distance vector of
/// `docs/specs/dists.md`, (0, 1), (0, 2), ..., (0, N-1), (1, 2), ..., the
/// upper triangle of the square matrix of the distances row by row.
#[derive(Debug)]
pub struct KosmanSums {
    /// How many individuals the source has, which is how many the pairs
    /// are of.
    num_individuals: usize,
    /// How many variants the calculation was given, called in a pair or
    /// not.
    num_vars: u64,
    /// How many alleles the genotype of one individual holds, the k the
    /// sum of d is multiplied by.
    ploidy: usize,
    /// The ploidy times the sum of d and n of every pair, in the order of
    /// the distance vector.
    sums: Vec<(u32, u32)>,
}

impl KosmanSums {
    /// The two counts of every pair of `num_individuals` individuals of the
    /// ploidy `ploidy` at 0, before any block has been added.
    ///
    /// # Errors
    ///
    /// When the individuals make more pairs than a `usize` counts, and when
    /// the machine does not give the memory of the pairs, 8 bytes each. The
    /// memory is asked for with `try_reserve_exact`, which gives it back as
    /// an error where `vec![(0, 0); n]` would end the process.
    fn at_zero(num_individuals: usize, ploidy: usize) -> Result<KosmanSums> {
        let num_pairs = num_pairs_of(num_individuals)
            .ok_or(Error::MorePairsThanAreCounted { num_individuals })?;
        let mut sums: Vec<(u32, u32)> = Vec::new();
        sums.try_reserve_exact(num_pairs)
            .map_err(|_| Error::DistancesOfTooManyIndividuals {
                num_individuals,
                num_pairs,
            })?;
        sums.resize(num_pairs, (0, 0));
        Ok(KosmanSums {
            num_individuals,
            num_vars: 0,
            ploidy,
            sums,
        })
    }

    /// How many individuals the pairs are of.
    #[must_use]
    pub fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// The variants the calculation was given, called in a pair or not.
    #[must_use]
    pub fn num_vars(&self) -> u64 {
        self.num_vars
    }

    /// How many alleles the genotype of one individual holds.
    #[must_use]
    pub fn ploidy(&self) -> usize {
        self.ploidy
    }

    /// The ploidy times the sum of d of the pair, and n, the variants at
    /// which both genotypes were called. The same for (i, j) and (j, i).
    /// `None` when i == j or when either is not an individual.
    #[must_use]
    pub fn sums(&self, i: usize, j: usize) -> Option<(u32, u32)> {
        self.sums.get(self.index_of_the_pair(i, j)?).copied()
    }

    /// The distance of the pair, the first of [`sums`](KosmanSums::sums)
    /// over the ploidy times the second. `None` when n is 0 or below
    /// `min_num_vars`, and where `sums` gives `None`.
    ///
    /// A pair with exactly `min_num_vars` variants keeps its distance.
    #[must_use]
    pub fn dist(&self, i: usize, j: usize, min_num_vars: u32) -> Option<f64> {
        let (k_sum, n) = self.sums(i, j)?;
        self.distance_of(k_sum, n, min_num_vars)
    }

    /// The distance of every pair, in the order of the distance vector.
    /// The binding crates write NaN for a `None`.
    pub fn dists(&self, min_num_vars: u32) -> impl Iterator<Item = Option<f64>> + '_ {
        self.sums
            .iter()
            .map(move |&(k_sum, n)| self.distance_of(k_sum, n, min_num_vars))
    }

    /// The distance of a pair whose two counts these are: `None` when n is
    /// 0 or below `min_num_vars`, and the sum of d over n otherwise.
    ///
    /// The three numbers are whole and far below 2^53, so each is exact in
    /// a `f64`, and there is one division. R's `gd.kosman` and pyNei divide
    /// twice: they add d, which is m over the ploidy, over the variants and
    /// then divide that sum by n. The two ways give the same bits when the
    /// ploidy is a power of two, where m over the ploidy is exact, which is
    /// the case of every dataset the two were compared on, of the ploidies
    /// 1, 2 and 4. At the ploidy 3 they differ in the last place: dividing
    /// by the ploidy and then by n, rather than by their product, gives
    /// another last bit in 4897 of 20000 pairs of an n from 100 to 2000 and
    /// a sum of d drawn at random, and popnei's, with its one rounding, is
    /// the nearer of the two to the exact value.
    fn distance_of(&self, k_sum: u32, n: u32, min_num_vars: u32) -> Option<f64> {
        if n == 0 || n < min_num_vars {
            return None;
        }
        Some(f64::from(k_sum) / (self.ploidy as f64 * f64::from(n)))
    }

    /// Where the pair of the individuals `first` and `second` is in the
    /// vector, in either order, and `None` when the two are one individual
    /// or either of them is not an individual.
    ///
    /// The pairs of the individual i start after those of the individuals
    /// before it, which are all the pairs of the source but those of the
    /// individuals from i on: two triangles of pairs, each of them the
    /// pairs of a number of individuals.
    fn index_of_the_pair(&self, first: usize, second: usize) -> Option<usize> {
        let (first, second) = if first < second {
            (first, second)
        } else {
            (second, first)
        };
        if first == second || second >= self.num_individuals {
            return None;
        }
        let from_first_on = num_pairs_of(self.num_individuals.checked_sub(first)?)?;
        let before_first = num_pairs_of(self.num_individuals)?.checked_sub(from_first_on)?;
        before_first.checked_add(second.checked_sub(first)?.checked_sub(1)?)
    }
}

/// How many pairs `num_individuals` individuals make, `n (n - 1) / 2`, and
/// `None` when that is more than a `usize` holds.
fn num_pairs_of(num_individuals: usize) -> Option<usize> {
    let Some(others) = num_individuals.checked_sub(1) else {
        // No individual makes no pair, which the subtraction below cannot
        // give.
        return Some(0);
    };
    // One of two consecutive numbers is even, so halving that one is exact,
    // and it is halved before the product: `n * (n - 1)` goes above what a
    // `usize` holds for numbers of individuals whose pairs it still counts,
    // 65537 of them in wasm, where a `usize` is 32 bits.
    if num_individuals.is_multiple_of(2) {
        (num_individuals / 2).checked_mul(others)
    } else {
        num_individuals.checked_mul(others / 2)
    }
}

/// The two counts of every pair of the block whose sets of bits `bits` are,
/// added into `sums`.
///
/// # Errors
///
/// When the block is not of the individuals or of the ploidy the sums are
/// of, which only a reader with a defect gives, and when a sum would go
/// above what a `u32` holds.
fn add_the_block(sums: &mut KosmanSums, bits: &KosmanBits) -> Result<()> {
    if bits.num_individuals() != sums.num_individuals
        || u32::try_from(sums.ploidy) != Ok(bits.ploidy)
    {
        return Err(Error::BlocksDoNotFitTogether {
            num_individuals: sums.num_individuals,
            ploidy: sums.ploidy,
            found_num_individuals: bits.num_individuals(),
            found_ploidy: usize::try_from(bits.ploidy).unwrap_or(usize::MAX),
        });
    }
    let (num_vars, ploidy) = (sums.num_vars, sums.ploidy);
    add_the_pairs_of_the_block(&mut sums.sums, bits, num_vars, ploidy)
}

/// The two counts of every pair of the block added into `sums`, the pairs
/// of one individual with the individuals after it on one thread of rayon.
///
/// `sums` is cut into one slice for each of those runs of pairs, so no two
/// threads write in the same place and none of them takes a lock. The runs
/// are as uneven as the upper triangle of a square matrix, the first
/// holding one pair for every other individual and the last holding one,
/// and there is one for each individual: rayon hands them out as the
/// threads ask for them. One work item for each pair was the other way, and
/// 5e7 pairs for 10000 individuals is that many items.
///
/// The counts are integers, so the sums do not depend on how the runs were
/// shared out.
///
/// # Errors
///
/// When a sum would go above what a `u32` holds. Which pair is found first
/// depends on the threads, and every pair that overflows gives the same
/// message, which names the variants and the ploidy and not the pair.
#[cfg(not(target_family = "wasm"))]
fn add_the_pairs_of_the_block(
    sums: &mut [(u32, u32)],
    bits: &KosmanBits,
    num_vars: u64,
    ploidy: usize,
) -> Result<()> {
    use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

    let num_pairs = sums.len();
    let mut rows: Vec<(usize, &mut [(u32, u32)])> = Vec::new();
    // One row for each individual, asked of the machine once for the whole
    // block, as the two counts of the pairs and the sets of bits are.
    rows.try_reserve_exact(bits.num_individuals())
        .map_err(|_| Error::DistancesOfTooManyIndividuals {
            num_individuals: bits.num_individuals(),
            num_pairs,
        })?;
    let mut rest = sums;
    for first in 0..bits.num_individuals() {
        // The pairs of the individual `first` are its pairs with each
        // individual after it. `add_the_block` checked that the block is of
        // the individuals the sums are of, so what is left of the vector
        // holds that run and the runs after it; the smaller of the two is
        // taken so that the split cannot panic.
        let of_the_row = rest.len().min(
            bits.num_individuals()
                .saturating_sub(first)
                .saturating_sub(1),
        );
        if of_the_row == 0 {
            break;
        }
        let (row, tail) = std::mem::take(&mut rest).split_at_mut(of_the_row);
        rows.push((first, row));
        rest = tail;
    }
    rows.par_iter_mut().try_for_each(|(first, row)| {
        for (pair, of_the_block) in row.iter_mut().zip(bits.sums_of_the_pairs_of(*first)) {
            add_the_pair(pair, of_the_block, num_vars, ploidy)?;
        }
        Ok(())
    })
}

/// The same sums, with the pairs added one after another, which is what
/// wasm does: it has no threads.
///
/// # Errors
///
/// When a sum would go above what a `u32` holds.
#[cfg(target_family = "wasm")]
fn add_the_pairs_of_the_block(
    sums: &mut [(u32, u32)],
    bits: &KosmanBits,
    num_vars: u64,
    ploidy: usize,
) -> Result<()> {
    add_the_pairs_one_by_one(sums, bits, num_vars, ploidy)
}

/// The pairs of the block added one after another.
///
/// It is compiled for every target and not for wasm alone, so that the
/// cargo tests, which run natively, can add the same blocks with it and
/// with the threads and compare what the two give.
///
/// # Errors
///
/// When a sum would go above what a `u32` holds.
#[cfg_attr(
    all(not(target_family = "wasm"), not(test)),
    expect(
        dead_code,
        reason = "in wasm it is how a block is added, and natively it is what the test that \
                  compares the two ways of adding one calls; outside the tests and outside \
                  wasm nothing calls it"
    )
)]
fn add_the_pairs_one_by_one(
    sums: &mut [(u32, u32)],
    bits: &KosmanBits,
    num_vars: u64,
    ploidy: usize,
) -> Result<()> {
    for (pair, of_the_block) in sums.iter_mut().zip(bits.sums_of_the_pairs()) {
        add_the_pair(pair, of_the_block, num_vars, ploidy)?;
    }
    Ok(())
}

/// The two counts of one pair over one block added into what the blocks
/// before it gave.
///
/// # Errors
///
/// When either sum would go above what a `u32` holds. `num_vars` and
/// `ploidy` are what the message names: the variants the pass has read so
/// far and the alleles of one genotype.
fn add_the_pair(
    pair: &mut (u32, u32),
    of_the_block: (u32, u32),
    num_vars: u64,
    ploidy: usize,
) -> Result<()> {
    let too_large = || Error::KosmanSumsTooLarge { num_vars, ploidy };
    pair.0 = pair.0.checked_add(of_the_block.0).ok_or_else(too_large)?;
    pair.1 = pair.1.checked_add(of_the_block.1).ok_or_else(too_large)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::path::{Path, PathBuf};

    use super::{KosmanBits, KosmanSums, add_the_block, calc_kosman_sums};
    use crate::block::{Block, BlockReader};
    use crate::error::{Error, Result};
    use crate::filters::FilteringStats;
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{ChromTable, MISSING_ALLELE, Needs};

    /// The allele of a genotype that was not called, short enough to read
    /// a table of genotypes with.
    const M: i8 = MISSING_ALLELE;

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

    /// The rows of a table of genotypes as the slices `block_of` takes.
    fn rows_of<const ALLELES: usize>(variants: &[[i8; ALLELES]]) -> Vec<&[i8]> {
        variants.iter().map(|row| row.as_slice()).collect()
    }

    /// Those rows cut into blocks of `num_vars_per_block` variants, the
    /// last of them shorter, which is how a reader gives them.
    fn blocks_of<const ALLELES: usize>(
        variants: &[[i8; ALLELES]],
        num_individuals: usize,
        ploidy: usize,
        num_vars_per_block: usize,
    ) -> Vec<Block> {
        variants
            .chunks(num_vars_per_block)
            .map(|chunk| block_of(&rows_of(chunk), num_individuals, ploidy))
            .collect()
    }

    /// The four variants of three diploid individuals of the worked
    /// example of "How it is verified" of `docs/specs/dists.md`: 0/0 0/1
    /// 1/1, 0/1 0/1 1/2, 0/0 0/. 2/2 and ./. 1/1 1/1. The third variant
    /// has the half called genotype, which is missing, and the second and
    /// the third hold three alleles where the first and the last hold two.
    const THE_DIPLOID_WORKED_EXAMPLE: [[i8; 6]; 4] = [
        [0, 0, 0, 1, 1, 1],
        [0, 1, 0, 1, 1, 2],
        [0, 0, 0, M, 2, 2],
        [M, M, 1, 1, 1, 1],
    ];

    /// The three variants of three tetraploid individuals of the same
    /// part of the spec: 0/0/0/1 0/1/1/1 1/1/1/1, 0/0/1/1 0/1/0/1 0/0/2/2
    /// and 0/0/0/0 0/0/./0 1/2/2/2.
    const THE_TETRAPLOID_WORKED_EXAMPLE: [[i8; 12]; 3] = [
        [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 2, 2],
        [0, 0, 0, 0, 0, 0, M, 0, 1, 2, 2, 2],
    ];

    /// The four variants of three haploid individuals of the same part of
    /// the spec: 0 0 1, 0 1 2, . 1 1 and 0 0 0.
    const THE_HAPLOID_WORKED_EXAMPLE: [[i8; 3]; 4] = [[0, 0, 1], [0, 1, 2], [M, 1, 1], [0, 0, 0]];

    /// The diploid worked example as one block.
    fn the_diploid_worked_example() -> Block {
        block_of(&rows_of(&THE_DIPLOID_WORKED_EXAMPLE), 3, 2)
    }

    /// The tetraploid worked example as one block.
    fn the_tetraploid_worked_example() -> Block {
        block_of(&rows_of(&THE_TETRAPLOID_WORKED_EXAMPLE), 3, 4)
    }

    /// The haploid worked example as one block.
    fn the_haploid_worked_example() -> Block {
        block_of(&rows_of(&THE_HAPLOID_WORKED_EXAMPLE), 3, 1)
    }

    /// The numbers of the table of the diploid worked example of "How it
    /// is verified" of `docs/specs/dists.md`: twice the sum of d and n of
    /// each of the three pairs. pyNei and `gd.kosman` of R give the same
    /// distances, 0.25, 0.833333 and 0.333333.
    #[test]
    fn the_diploid_worked_example_gives_the_sums_of_the_spec() {
        let bits = KosmanBits::of_block(&the_diploid_worked_example()).unwrap();

        assert_eq!(bits.sums_of_one_pair(0, 1), Some((1, 2)));
        assert_eq!(bits.sums_of_one_pair(0, 2), Some((5, 3)));
        assert_eq!(bits.sums_of_one_pair(1, 2), Some((2, 3)));
    }

    /// The numbers of the table of the tetraploid worked example, which
    /// `gd.kosman` of R gives: four times the sum of d and n. The ploidy
    /// is read here and not taken for 2, and a genotype of four alleles
    /// holds one allele twice and more.
    #[test]
    fn the_tetraploid_worked_example_gives_the_sums_of_the_spec() {
        let bits = KosmanBits::of_block(&the_tetraploid_worked_example()).unwrap();

        assert_eq!(bits.sums_of_one_pair(0, 1), Some((2, 2)));
        assert_eq!(bits.sums_of_one_pair(0, 2), Some((9, 3)));
        assert_eq!(bits.sums_of_one_pair(1, 2), Some((3, 2)));
    }

    /// The numbers of the table of the haploid worked example, which
    /// `gd.kosman` of R gives: the ploidy is 1, so the first number is the
    /// sum of d itself, 0 for two equal alleles and 1 for two different
    /// ones.
    #[test]
    fn the_haploid_worked_example_gives_the_sums_of_the_spec() {
        let bits = KosmanBits::of_block(&the_haploid_worked_example()).unwrap();

        assert_eq!(bits.sums_of_one_pair(0, 1), Some((1, 3)));
        assert_eq!(bits.sums_of_one_pair(0, 2), Some((2, 3)));
        assert_eq!(bits.sums_of_one_pair(1, 2), Some((2, 4)));
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

        assert_eq!(bits.sums_of_one_pair(2, 0), bits.sums_of_one_pair(0, 2));
        assert_eq!(bits.sums_of_one_pair(2, 1), bits.sums_of_one_pair(1, 2));
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

        assert_eq!(bits.sums_of_one_pair(0, 1), Some((0, 0)));
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
        assert_eq!(bits.sums_of_one_pair(0, 0), None);
        assert_eq!(bits.sums_of_one_pair(0, 1), None);
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

        assert_eq!(bits.sums_of_one_pair(0, 1), Some((100, 100)));
        assert_eq!(bits.sums_of_one_pair(0, 2), Some((100, 50)));
        assert_eq!(bits.sums_of_one_pair(1, 2), Some((50, 50)));
    }

    /// Every allele of a multiallelic variant counts as itself, so a
    /// variant of the alleles 0 and 5 beside two variants of the alleles 0
    /// and 1 has the sets of the allele 5 too: 5/5 and 0/5 share one copy
    /// of the 5 and have d of 0.5, as 0/1 and 0/0 do.
    #[test]
    fn a_multiallelic_variant_beside_biallelic_ones_counts_each_of_its_alleles() {
        let block = block_of(&[&[0, 1, 0, 0], &[5, 5, 0, 5], &[0, 1, 1, 0]], 2, 2);

        let bits = KosmanBits::of_block(&block).unwrap();

        assert_eq!(bits.sums_of_one_pair(0, 1), Some((2, 3)));
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

    /// The sets of an individual are as many as the alleles the block
    /// holds and not as many as its largest allele value: a block whose
    /// two alleles are 0 and 127 gives every pair the same counts as the
    /// same block with the allele 1 in the place of the 127, and its
    /// `holds` sets are the same 4 words, 2 alleles times the ploidy 2
    /// times one word, and not 256.
    #[test]
    fn a_block_of_two_far_apart_alleles_has_the_sets_of_two_alleles() {
        let far_apart = block_of(&[&[0, 127, 0, 0], &[0, 0, 127, 127]], 2, 2);
        let beside_each_other = block_of(&[&[0, 1, 0, 0], &[0, 0, 1, 1]], 2, 2);

        let far_apart = KosmanBits::of_block(&far_apart).unwrap();
        let beside_each_other = KosmanBits::of_block(&beside_each_other).unwrap();

        assert_eq!(
            far_apart.sums_of_one_pair(0, 1),
            beside_each_other.sums_of_one_pair(0, 1)
        );
        assert_eq!(far_apart.sums_of_one_pair(0, 1), Some((3, 2)));
        assert_eq!(far_apart.holds_per_individual.get(), 4);
        assert_eq!(
            far_apart.holds_per_individual,
            beside_each_other.holds_per_individual
        );
    }

    /// A block whose genotypes hold the alleles 0, 5 and 9 gets the sets of
    /// three alleles, one for each value it holds and none for the values
    /// between them, and every pair keeps the counts of the spec: 5/5 and
    /// 0/5 share one copy of the 5.
    #[test]
    fn the_sets_of_a_block_are_as_many_as_the_allele_values_it_holds() {
        let block = block_of(&[&[0, 9, 0, 0], &[5, 5, 0, 5], &[0, 9, 9, 0]], 2, 2);

        let bits = KosmanBits::of_block(&block).unwrap();

        // Three alleles times the ploidy 2, one word each.
        assert_eq!(bits.holds_per_individual.get(), 6);
        assert_eq!(bits.sums_of_one_pair(0, 1), Some((2, 3)));
    }

    /// A genotype holds one allele at least, so a block of the ploidy 0 is
    /// not one whose genotypes can be read. Only a reader with a defect
    /// gives one.
    #[test]
    fn a_block_of_the_ploidy_0_is_an_error() {
        let block = block_of(&[&[0, 1]], 2, 0);

        let error = KosmanBits::of_block(&block).unwrap_err();

        assert!(
            matches!(
                error,
                Error::GtsNotWholeGenotypes {
                    num_alleles: 2,
                    ploidy: 0
                }
            ),
            "{error}"
        );
    }

    /// How many genotypes the block below holds, its 100 variants of its
    /// 40 individuals.
    const GENOTYPES_OF_THE_BLOCK_OF_40: usize = 4000;

    /// A block of 100 variants of 40 diploid individuals, drawn by a
    /// generator of its own so that it is the same block on every run and
    /// on every machine.
    ///
    /// The 64 variants a word of a set holds do not divide its 100
    /// variants, so the last word of every set is partial, and the
    /// individuals of a work item do not divide its 40 individuals, so the
    /// last item of the threads is short. One genotype in six is missing
    /// whole and one in six is half called, so a third of them set no bit,
    /// and the alleles are 0, 1 and 2, in homozygous and heterozygous
    /// genotypes.
    fn a_block_of_100_variants_of_40_individuals() -> Block {
        let mut gts: Vec<i8> = Vec::new();
        // A linear congruential generator, of which this test needs
        // nothing but that it draws the same genotypes every time.
        let mut drawn: u32 = 1_234_567;
        for _ in 0..GENOTYPES_OF_THE_BLOCK_OF_40 {
            drawn = drawn.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            match (drawn >> 24) % 6 {
                0 => gts.extend_from_slice(&[M, M]),
                1 => gts.extend_from_slice(&[0, M]),
                2 => gts.extend_from_slice(&[0, 0]),
                3 => gts.extend_from_slice(&[0, 1]),
                4 => gts.extend_from_slice(&[1, 1]),
                _ => gts.extend_from_slice(&[1, 2]),
            }
        }
        Block {
            num_vars: 100,
            num_individuals: 40,
            ploidy: 2,
            gts,
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        }
    }

    /// The bits the threads write are the bits of the pass that writes one
    /// individual after another, which is what wasm runs and what the
    /// worked examples of the spec were checked against. Each individual
    /// owns its own words of the two arrays, so the threads write no word
    /// twice, and the block above makes the two places that could go wrong
    /// reachable: a last word of a set that is partial, and a last work
    /// item that holds fewer individuals than the others.
    #[test]
    fn the_sets_built_on_the_threads_are_the_ones_built_one_after_another() {
        let block = a_block_of_100_variants_of_40_individuals();

        let on_the_threads = KosmanBits::of_block(&block).unwrap();
        let one_after_another = KosmanBits::of_block_one_by_one(&block).unwrap();

        assert_eq!(on_the_threads.called, one_after_another.called);
        assert_eq!(on_the_threads.holds, one_after_another.holds);
        // Two arrays of zeroes would be equal whatever either pass did, so
        // the block has to have set bits for the test to be able to fail.
        assert!(on_the_threads.called.iter().any(|word| *word != 0));
        assert!(on_the_threads.holds.iter().any(|word| *word != 0));
    }

    /// Every block a reader of popnei gives holds one variant at least, and
    /// a block of none would leave the sets with no word to write a bit in.
    #[test]
    fn a_block_of_no_variants_is_an_error() {
        let block = block_of(&[], 3, 2);

        let error = KosmanBits::of_block(&block).unwrap_err();

        assert!(
            matches!(error, Error::ReaderGaveABlockOfNoVariants),
            "{error}"
        );
    }

    /// An allele below the missing one, which no reader of popnei gives, has
    /// no place among the alleles of the block. Read without the guard it
    /// would be dropped while its variant counted as called for the pair,
    /// which is a wrong number and no message, so it is an error.
    #[test]
    fn a_genotype_with_an_allele_below_the_missing_one_is_an_error() {
        let block = block_of(&[&[0, 0, -2, 0]], 2, 2);

        let error = KosmanBits::of_block(&block).unwrap_err();

        assert!(
            matches!(error, Error::AlleleBelowTheMissingOne { allele: -2 }),
            "{error}"
        );
    }

    /// A reader of blocks written for these tests: it gives the blocks it
    /// was built with, keeps what it was last asked to fill, and gives an
    /// error instead of a block at the call a test names, so that a test
    /// sees what the calculation does with the error of a reader.
    struct GivenBlocks {
        individuals: Vec<String>,
        ploidy: usize,
        chroms: ChromTable,
        /// The blocks still to give, the last one first.
        left: Vec<Block>,
        /// The call at which it gives an error instead of a block, counted
        /// from 1.
        fails_at: Option<usize>,
        /// How many times it was asked for a block.
        calls: usize,
        /// What it was last asked to fill.
        needs: Needs,
    }

    impl GivenBlocks {
        /// A reader of `num_individuals` individuals of the ploidy
        /// `ploidy`, named `ind0` and on, that gives `blocks` in their
        /// order.
        fn of(blocks: Vec<Block>, num_individuals: usize, ploidy: usize) -> GivenBlocks {
            let mut left = blocks;
            left.reverse();
            GivenBlocks {
                individuals: (0..num_individuals).map(|at| format!("ind{at}")).collect(),
                ploidy,
                chroms: ChromTable::new(),
                left,
                fails_at: None,
                calls: 0,
                needs: Needs::ALL,
            }
        }

        /// The same reader, whose call number `call` is an error.
        fn failing_at(blocks: Vec<Block>, num_individuals: usize, call: usize) -> GivenBlocks {
            GivenBlocks {
                fails_at: Some(call),
                ..GivenBlocks::of(blocks, num_individuals, 2)
            }
        }
    }

    /// What the reader of the tests says when a test asked it to fail.
    const THE_READER_FAILED: &str = "the reader of the tests failed";

    impl BlockReader for GivenBlocks {
        fn next_block(&mut self) -> Result<Option<Block>> {
            self.calls = self.calls.saturating_add(1);
            if self.fails_at == Some(self.calls) {
                return Err(Error::Io(std::io::Error::other(THE_READER_FAILED)));
            }
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

        fn set_needs(&mut self, needs: Needs) {
            self.needs = needs;
        }

        fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
            Vec::new()
        }
    }

    /// The two counts of every pair, in the order of the distance vector.
    fn every_sum(sums: &KosmanSums) -> Vec<(u32, u32)> {
        let mut every = Vec::new();
        for first in 0..sums.num_individuals() {
            for second in first.saturating_add(1)..sums.num_individuals() {
                match sums.sums(first, second) {
                    Some(of_the_pair) => every.push(of_the_pair),
                    None => panic!("the pair {first}, {second} is not there"),
                }
            }
        }
        every
    }

    /// The distances of the pairs, `None` where a pair has none, compared
    /// within 1e-9, the digits of the shortest values of the spec's table.
    fn assert_the_distances_are(found: &[Option<f64>], expected: &[Option<f64>]) {
        assert_eq!(
            found.len(),
            expected.len(),
            "{found:?} against {expected:?}"
        );
        for (at, (found, expected)) in found.iter().zip(expected).enumerate() {
            match (found, expected) {
                (Some(found), Some(expected)) => assert!(
                    (found - expected).abs() < 1e-9,
                    "the pair {at} has {found} and not {expected}"
                ),
                (None, None) => {}
                _ => panic!("the pair {at} has {found:?} and not {expected:?}"),
            }
        }
    }

    /// The four reference files of the Kosman distances and what R gives
    /// for them live at the root of the repository, beside the Python
    /// tests that read the same files, and not inside this crate. The path
    /// is built from the directory of the manifest, so it holds whether
    /// the tests are run with `cargo test --workspace` or with `cargo test
    /// -p popnei`.
    fn reference_dists(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/dists")
            .join(name)
    }

    /// A VCF reader over one of the four reference files, of the ploidy of
    /// its dataset and with blocks of the size a test asks for.
    fn reader_of_the_reference(
        name: &str,
        ploidy: usize,
        num_vars_per_block: Option<usize>,
    ) -> VcfReader<std::io::BufReader<std::fs::File>> {
        let options = VcfOptions {
            ploidy,
            num_vars_per_block,
            ..VcfOptions::default()
        };
        match VcfReader::from_path(&reference_dists(name), options) {
            Ok(reader) => reader,
            Err(error) => panic!("{name}: {error}"),
        }
    }

    /// The two counts of every pair of one of the reference files, with
    /// the names of its individuals in the order of the file.
    fn sums_of_the_reference(
        name: &str,
        ploidy: usize,
        num_vars_per_block: Option<usize>,
    ) -> (Vec<String>, KosmanSums) {
        let mut reader = reader_of_the_reference(name, ploidy, num_vars_per_block);
        let individuals = reader.individuals().to_vec();
        match calc_kosman_sums(&mut reader) {
            Ok(sums) => (individuals, sums),
            Err(error) => panic!("{name}: {error}"),
        }
    }

    /// Where the individual of that name is in the file.
    fn individual_at(individuals: &[String], name: &str) -> usize {
        match individuals.iter().position(|held| held == name) {
            Some(at) => at,
            None => panic!("{name} is not an individual of the file"),
        }
    }

    /// The pair of the two named individuals has the ploidy times the sum
    /// of d, the n and the distance of one line of the table of "How it is
    /// verified" of `docs/specs/dists.md`: the integers exactly and the
    /// distance within 1e-9.
    fn assert_the_pair_is(
        sums: &KosmanSums,
        individuals: &[String],
        of: (&str, &str),
        k_sum: u32,
        n: u32,
        dist: f64,
    ) {
        let (one, other) = (
            individual_at(individuals, of.0),
            individual_at(individuals, of.1),
        );
        assert_eq!(sums.sums(one, other), Some((k_sum, n)), "{of:?}");
        match sums.dist(one, other, 0) {
            Some(found) => assert!(
                (found - dist).abs() < 1e-9,
                "{of:?} has the distance {found} and not {dist}"
            ),
            None => panic!("{of:?} has no distance"),
        }
    }

    /// The five pairs of the panel in the table of "How it is verified" of
    /// `docs/specs/dists.md`, which `gd.kosman` of R gives: 200 diploid
    /// individuals over 1200 biallelic variants, 3 in 100 genotypes
    /// missing, read with the VCF reader.
    #[test]
    fn the_panel_of_the_reference_has_the_five_pairs_of_the_spec() {
        let (individuals, sums) = sums_of_the_reference("panel.vcf.gz", 2, None);

        assert_eq!(sums.num_individuals(), 200);
        assert_eq!(sums.num_vars(), 1200);
        assert_eq!(sums.ploidy(), 2);
        assert_the_pair_is(
            &sums,
            &individuals,
            ("s000", "s001"),
            372,
            1122,
            0.1657754010695187,
        );
        assert_the_pair_is(
            &sums,
            &individuals,
            ("s000", "s002"),
            376,
            1128,
            0.16666666666666666,
        );
        assert_the_pair_is(
            &sums,
            &individuals,
            ("s198", "s199"),
            351,
            1134,
            0.15476190476190477,
        );
        assert_the_pair_is(
            &sums,
            &individuals,
            ("s010", "s033"),
            804,
            1123,
            0.3579697239536955,
        );
        assert_the_pair_is(
            &sums,
            &individuals,
            ("s116", "s119"),
            310,
            1133,
            0.13680494263018536,
        );
    }

    /// The three pairs of the 4 allele dataset in the same table: 40
    /// diploid individuals over 300 variants of up to four alleles, where
    /// a pair of genotypes can share no allele although both are
    /// heterozygous.
    #[test]
    fn the_four_allele_file_of_the_reference_has_the_three_pairs_of_the_spec() {
        let (individuals, sums) = sums_of_the_reference("four_alleles.vcf.gz", 2, None);

        assert_eq!(sums.num_individuals(), 40);
        assert_eq!(sums.num_vars(), 300);
        assert_the_pair_is(&sums, &individuals, ("i00", "i01"), 325, 269, 0.6040892193);
        assert_the_pair_is(&sums, &individuals, ("i00", "i02"), 304, 264, 0.5757575758);
        assert_the_pair_is(&sums, &individuals, ("i00", "i03"), 347, 272, 0.6378676471);
    }

    /// The three pairs of the tetraploid dataset in the same table: the
    /// ploidy is read from the reader and is 4, so the first integer is
    /// four times the sum of d.
    #[test]
    #[expect(
        clippy::excessive_precision,
        reason = "the distances are the literals of the table of \"How it is verified\" of \
                  docs/specs/dists.md, which are what R printed with 17 digits; the digits \
                  beyond what a f64 keeps read back as the same value, and the literals are \
                  the spec's"
    )]
    fn the_tetraploid_file_of_the_reference_has_the_three_pairs_of_the_spec() {
        let (individuals, sums) = sums_of_the_reference("tetraploid.vcf.gz", 4, None);

        assert_eq!(sums.num_individuals(), 12);
        assert_eq!(sums.num_vars(), 200);
        assert_eq!(sums.ploidy(), 4);
        assert_the_pair_is(&sums, &individuals, ("t00", "t01"), 282, 188, 0.375);
        assert_the_pair_is(
            &sums,
            &individuals,
            ("t00", "t02"),
            284,
            183,
            0.38797814207650272,
        );
        assert_the_pair_is(
            &sums,
            &individuals,
            ("t00", "t03"),
            294,
            183,
            0.40163934426229508,
        );
    }

    /// The three pairs of the haploid dataset in the same table: the
    /// ploidy is 1, so the first integer is the sum of d itself.
    #[test]
    #[expect(
        clippy::excessive_precision,
        reason = "the distances are the literals of the table of \"How it is verified\" of \
                  docs/specs/dists.md, which are what R printed with 17 digits; the digits \
                  beyond what a f64 keeps read back as the same value, and the literals are \
                  the spec's"
    )]
    fn the_haploid_file_of_the_reference_has_the_three_pairs_of_the_spec() {
        let (individuals, sums) = sums_of_the_reference("haploid.vcf.gz", 1, None);

        assert_eq!(sums.num_individuals(), 12);
        assert_eq!(sums.num_vars(), 200);
        assert_eq!(sums.ploidy(), 1);
        assert_the_pair_is(
            &sums,
            &individuals,
            ("h00", "h01"),
            112,
            180,
            0.62222222222222223,
        );
        assert_the_pair_is(
            &sums,
            &individuals,
            ("h00", "h02"),
            123,
            184,
            0.66847826086956519,
        );
        assert_the_pair_is(
            &sums,
            &individuals,
            ("h00", "h03"),
            113,
            179,
            0.63128491620111726,
        );
    }

    /// What `gd.kosman` of R gave for one of the reference files: the
    /// distance, n and k times the sum of d of every pair, in the order of
    /// the distance vector, from the header line on.
    fn gdkosman_of(name: &str) -> Vec<(f64, u32, u32)> {
        let path = reference_dists(name);
        let text = match std::fs::read_to_string(&path) {
            Ok(text) => text,
            Err(error) => panic!("{}: {error}", path.display()),
        };
        text.lines()
            .skip(1)
            .map(|line| {
                let mut fields = line.split('\t');
                let mut number = || match fields.next() {
                    Some(field) => field.to_owned(),
                    None => panic!("{name}: the line `{line}` has fewer than three fields"),
                };
                let (dist, n, k_sum) = (number(), number(), number());
                (
                    dist.parse().unwrap_or_else(|_| panic!("{dist}")),
                    n.parse().unwrap_or_else(|_| panic!("{n}")),
                    k_sum.parse().unwrap_or_else(|_| panic!("{k_sum}")),
                )
            })
            .collect()
    }

    /// Every pair of the four reference files, 20812 of them, has the two
    /// integers that `gd.kosman` of R gives it and a distance within 1e-9
    /// of R's, where the tests above assert the fourteen pairs that the
    /// spec's table names.
    #[test]
    fn every_pair_of_the_four_reference_files_has_the_sums_and_the_distance_of_r() {
        let files = [
            ("panel", 2, 19900),
            ("four_alleles", 2, 780),
            ("tetraploid", 4, 66),
            ("haploid", 1, 66),
        ];
        for (name, ploidy, num_pairs) in files {
            let of_r = gdkosman_of(&format!("{name}.gdkosman.tsv"));
            assert_eq!(of_r.len(), num_pairs, "{name}");
            let (_, sums) = sums_of_the_reference(&format!("{name}.vcf.gz"), ploidy, None);

            let of_popnei = every_sum(&sums);
            let dists: Vec<Option<f64>> = sums.dists(0).collect();
            assert_eq!(of_popnei.len(), num_pairs, "{name}");
            for (at, (dist, n, k_sum)) in of_r.into_iter().enumerate() {
                assert_eq!(
                    of_popnei.get(at),
                    Some(&(k_sum, n)),
                    "{name}, the pair {at}"
                );
                match dists.get(at) {
                    Some(&Some(found)) => assert!(
                        (found - dist).abs() < 1e-9,
                        "{name}, the pair {at}: {found} is not {dist}"
                    ),
                    found => panic!("{name}, the pair {at} has {found:?} and not {dist}"),
                }
            }
        }
    }

    /// The diploid worked example of "How it is verified" of
    /// `docs/specs/dists.md` read from a reader, one variant to a block,
    /// so that the two counts of every pair are added over four blocks:
    /// the distances are 0.25, 0.833333 and 0.333333, which pyNei and
    /// `gd.kosman` give too.
    #[test]
    fn the_diploid_worked_example_read_in_blocks_gives_the_distances_of_the_spec() {
        let blocks = blocks_of(&THE_DIPLOID_WORKED_EXAMPLE, 3, 2, 1);
        let mut reader = GivenBlocks::of(blocks, 3, 2);

        let sums = calc_kosman_sums(&mut reader).expect("the sums");

        assert_eq!(sums.num_vars(), 4);
        assert_eq!(every_sum(&sums), [(1, 2), (5, 3), (2, 3)]);
        assert_the_distances_are(
            &sums.dists(0).collect::<Vec<_>>(),
            &[
                Some(0.25),
                Some(0.8333333333333334),
                Some(0.3333333333333333),
            ],
        );
    }

    /// A pair needs `min_num_vars` variants called in both to get a
    /// distance, and a pair with exactly that many keeps it. The first
    /// pair of the diploid worked example has 2 variants and the other two
    /// have 3, so 3 takes the first pair's distance away and 4 takes them
    /// all.
    #[test]
    fn a_pair_of_fewer_variants_than_min_num_vars_has_no_distance() {
        let mut reader = GivenBlocks::of(vec![the_diploid_worked_example()], 3, 2);

        let sums = calc_kosman_sums(&mut reader).expect("the sums");

        assert_the_distances_are(
            &sums.dists(3).collect::<Vec<_>>(),
            &[None, Some(0.8333333333333334), Some(0.3333333333333333)],
        );
        assert_the_distances_are(&sums.dists(4).collect::<Vec<_>>(), &[None, None, None]);
        assert_eq!(sums.dist(0, 1, 2), Some(0.25));
        assert_eq!(sums.dist(0, 1, 3), None);
    }

    /// The tetraploid worked example read from a reader, in blocks of two
    /// variants: the distances are 0.25, 0.75 and 0.375, which
    /// `gd.kosman` gives.
    #[test]
    fn the_tetraploid_worked_example_read_in_blocks_gives_the_distances_of_the_spec() {
        let blocks = blocks_of(&THE_TETRAPLOID_WORKED_EXAMPLE, 3, 4, 2);
        let mut reader = GivenBlocks::of(blocks, 3, 4);

        let sums = calc_kosman_sums(&mut reader).expect("the sums");

        assert_eq!(sums.ploidy(), 4);
        assert_eq!(every_sum(&sums), [(2, 2), (9, 3), (3, 2)]);
        assert_the_distances_are(
            &sums.dists(0).collect::<Vec<_>>(),
            &[Some(0.25), Some(0.75), Some(0.375)],
        );
    }

    /// The haploid worked example read from a reader, in blocks of three
    /// variants: the distances are 0.333333, 0.666667 and 0.5, which
    /// `gd.kosman` gives.
    #[test]
    fn the_haploid_worked_example_read_in_blocks_gives_the_distances_of_the_spec() {
        let blocks = blocks_of(&THE_HAPLOID_WORKED_EXAMPLE, 3, 1, 3);
        let mut reader = GivenBlocks::of(blocks, 3, 1);

        let sums = calc_kosman_sums(&mut reader).expect("the sums");

        assert_eq!(sums.ploidy(), 1);
        assert_eq!(every_sum(&sums), [(1, 3), (2, 3), (2, 4)]);
        assert_the_distances_are(
            &sums.dists(0).collect::<Vec<_>>(),
            &[
                Some(0.3333333333333333),
                Some(0.6666666666666666),
                Some(0.5),
            ],
        );
    }

    /// The two counts are integers that each block adds to, so the size of
    /// the blocks changes nothing: the 4 allele file in blocks of 7, of
    /// 64, of 65 and of 300 variants gives the same integers, the three
    /// sizes around a word of 64 bits among them.
    #[test]
    fn the_size_of_the_blocks_does_not_change_the_sums() {
        let (_, whole) = sums_of_the_reference("four_alleles.vcf.gz", 2, Some(300));
        let of_the_whole_file = every_sum(&whole);
        assert_eq!(of_the_whole_file.len(), 780);

        for num_vars_per_block in [7, 64, 65] {
            let (_, cut) =
                sums_of_the_reference("four_alleles.vcf.gz", 2, Some(num_vars_per_block));

            assert_eq!(cut.num_vars(), 300, "blocks of {num_vars_per_block}");
            assert_eq!(
                every_sum(&cut),
                of_the_whole_file,
                "blocks of {num_vars_per_block}"
            );
        }
    }

    /// The pairs of a block are counted on the threads of the pool the
    /// caller is in, and the counts are integers, so one thread and four
    /// give the same integers.
    ///
    /// The pools are built here and are not rayon's global one, which has
    /// one thread per core of the machine. rayon is a dependency of the
    /// targets that are not wasm, so this test is compiled for those
    /// alone.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_number_of_threads_does_not_change_the_sums() {
        let in_a_pool = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            pool.install(|| {
                let (_, sums) = sums_of_the_reference("four_alleles.vcf.gz", 2, Some(64));
                every_sum(&sums)
            })
        };

        let on_one = in_a_pool(1);
        assert_eq!(on_one.len(), 780);
        assert_eq!(on_one, in_a_pool(4));
    }

    /// The pairs shared out over the threads are the pairs walked one
    /// after another, which is what wasm does: the two ways of adding a
    /// block give the same integers over the 780 pairs of the 4 allele
    /// file.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_pairs_added_on_the_threads_are_the_ones_added_one_after_another() {
        use super::add_the_pairs_one_by_one;

        let mut reader = reader_of_the_reference("four_alleles.vcf.gz", 2, Some(64));
        let mut one_after_another = KosmanSums::at_zero(40, 2).expect("the sums");
        let mut on_the_threads = KosmanSums::at_zero(40, 2).expect("the sums");

        while let Some(block) = reader.next_block().expect("the block") {
            let bits = KosmanBits::of_block(&block).expect("the sets of bits");
            add_the_pairs_one_by_one(&mut one_after_another.sums, &bits, 0, 2).expect("the sums");
            add_the_block(&mut on_the_threads, &bits).expect("the sums");
        }

        assert_eq!(every_sum(&one_after_another).len(), 780);
        assert_eq!(every_sum(&one_after_another), every_sum(&on_the_threads));
    }

    /// The calculation asks its reader for the genotypes and for nothing
    /// else, so a reader over a file leaves the columns of a variant
    /// unparsed.
    #[test]
    fn the_calculation_asks_its_reader_for_the_genotypes_alone() {
        let mut reader = GivenBlocks::of(vec![the_diploid_worked_example()], 3, 2);

        calc_kosman_sums(&mut reader).expect("the sums");

        assert_eq!(reader.needs, Needs::GTS);
    }

    /// Every variant the calculation was given counts in `num_vars`,
    /// whether or not it was called in a pair: a block whose genotypes are
    /// all missing adds its variants to the count and nothing to the pairs.
    #[test]
    fn num_vars_counts_the_variants_that_no_pair_was_called_at() {
        let blocks = vec![
            the_diploid_worked_example(),
            block_of(&[&[M; 6], &[M; 6]], 3, 2),
        ];
        let mut reader = GivenBlocks::of(blocks, 3, 2);

        let sums = calc_kosman_sums(&mut reader).expect("the sums");

        assert_eq!(sums.num_vars(), 6);
        assert_eq!(every_sum(&sums), [(1, 2), (5, 3), (2, 3)]);
    }

    /// One individual with itself is not a pair, and neither is one with
    /// an individual that is not in the source: both have no counts and no
    /// distance.
    #[test]
    fn an_individual_with_itself_and_one_that_is_not_there_have_no_sums() {
        let mut reader = GivenBlocks::of(vec![the_diploid_worked_example()], 3, 2);

        let sums = calc_kosman_sums(&mut reader).expect("the sums");

        assert_eq!(sums.sums(1, 1), None);
        assert_eq!(sums.sums(0, 3), None);
        assert_eq!(sums.sums(3, 0), None);
        assert_eq!(sums.dist(1, 1, 0), None);
        assert_eq!(sums.dist(0, 3, 0), None);
    }

    /// The two counts of a pair are the same whichever of the two
    /// individuals is named first, and the pairs are in the order of the
    /// distance vector.
    #[test]
    fn the_sums_of_a_pair_over_a_pass_are_the_same_in_either_order() {
        let mut reader = GivenBlocks::of(vec![the_diploid_worked_example()], 3, 2);

        let sums = calc_kosman_sums(&mut reader).expect("the sums");

        assert_eq!(sums.sums(2, 0), sums.sums(0, 2));
        assert_eq!(sums.sums(2, 1), sums.sums(1, 2));
        assert_eq!(sums.sums(0, 2), Some((5, 3)));
    }

    /// A reader with no variant is an error: there is nothing to calculate
    /// over, and the memory of the pairs is never asked for.
    #[test]
    fn a_reader_with_no_variant_is_an_error() {
        let mut reader = GivenBlocks::of(Vec::new(), 3, 2);

        let error = calc_kosman_sums(&mut reader).expect_err("the error");

        assert!(matches!(error, Error::ReaderGaveNoVariants), "{error}");
    }

    /// popnei keeps two `u32` for each pair, and a sum that would go above
    /// what a `u32` holds is an error and does not wrap. The sums are set
    /// by hand here: reaching the end of a `u32` takes more than two
    /// thousand million variants of diploids.
    #[test]
    fn a_sum_above_what_a_u32_holds_is_an_error() {
        let mut sums = KosmanSums::at_zero(3, 2).expect("the sums");
        sums.num_vars = 2_147_483_648;
        sums.sums = vec![(u32::MAX, 5), (0, 0), (0, 0)];
        let bits = KosmanBits::of_block(&the_diploid_worked_example()).expect("the sets of bits");

        let error = add_the_block(&mut sums, &bits).expect_err("the error");

        assert!(
            matches!(
                error,
                Error::KosmanSumsTooLarge {
                    num_vars: 2_147_483_648,
                    ploidy: 2
                }
            ),
            "{error}"
        );
    }

    /// The error of the reader is given on as it is: the calculation adds
    /// nothing of its own to it and does not turn it into an error of its
    /// own.
    #[test]
    fn an_error_of_the_reader_is_given_on_as_it_is() {
        let blocks = blocks_of(&THE_DIPLOID_WORKED_EXAMPLE, 3, 2, 1);
        let mut reader = GivenBlocks::failing_at(blocks, 3, 2);

        let error = calc_kosman_sums(&mut reader).expect_err("the error");

        assert!(
            matches!(&error, Error::Io(of_the_reader) if of_the_reader.to_string() == THE_READER_FAILED),
            "{error}"
        );
        assert_eq!(reader.calls, 2);
    }

    /// A reader that gives a block of other individuals than it says its
    /// source has would have its genotypes read one individual at the
    /// place of another, so it is an error. Only a reader with a defect
    /// gives one.
    #[test]
    fn a_block_of_other_individuals_than_the_reader_says_is_an_error() {
        let blocks = vec![
            the_diploid_worked_example(),
            block_of(&[&[0, 0, 1, 1]], 2, 2),
        ];
        let mut reader = GivenBlocks::of(blocks, 3, 2);

        let error = calc_kosman_sums(&mut reader).expect_err("the error");

        assert!(
            matches!(
                error,
                Error::BlocksDoNotFitTogether {
                    num_individuals: 3,
                    ploidy: 2,
                    found_num_individuals: 2,
                    found_ploidy: 2
                }
            ),
            "{error}"
        );
    }

    /// More individuals than the pairs a `usize` counts is an error before
    /// any memory is asked for: popnei gives each pair a place among the
    /// others, and a place is counted in a `usize`.
    #[test]
    fn more_individuals_than_the_pairs_a_usize_counts_is_an_error() {
        let error = KosmanSums::at_zero(usize::MAX, 2).expect_err("the error");

        assert!(
            matches!(
                error,
                Error::MorePairsThanAreCounted {
                    num_individuals: usize::MAX
                }
            ),
            "{error}"
        );
    }

    /// One individual makes no pair, so the distances of a source of one
    /// are an empty vector, and the variants are counted all the same.
    #[test]
    fn a_source_of_one_individual_has_no_pair() {
        let mut reader = GivenBlocks::of(vec![block_of(&[&[0, 1], &[0, 0]], 1, 2)], 1, 2);

        let sums = calc_kosman_sums(&mut reader).expect("the sums");

        assert_eq!(sums.num_individuals(), 1);
        assert_eq!(sums.num_vars(), 2);
        assert_eq!(sums.dists(0).count(), 0);
    }
}
