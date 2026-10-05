//! The filter of `docs/specs/filters.md` that keeps variants at random,
//! each with the probability a user gives, so that a calculation runs on a
//! sample spread over the whole of a large file.
//!
//! Whether a variant is kept is drawn from SplitMix64, the generator of
//! random numbers that `java.util.SplittableRandom` of Java uses, started at
//! the user's seed in every pass and drawn once for each variant in the
//! order the variants reach the filter, so that every pass of one
//! `Variants` keeps the same variants. [`RandomFilter`] is the filter of one
//! pass, which keeps the variants of a block and counts them, and
//! [`RandomlyFilteredReader`] the reader that puts it over a source.

use std::fmt;

use crate::block::{Block, BlockReader, SourceHeader};
use crate::error::{Error, Result};
use crate::filters::{FilteringStats, RegionSelection};
use crate::variant::{ChromTable, Needs};

/// The name under which the counts of the filter that keeps variants at
/// random reach a Python or a TypeScript user, and by which a chain of
/// readers is asked whether it holds one already.
pub(crate) const RANDOM_KIND: &str = "random";

/// 42, the seed of the filter that keeps variants at random when the user
/// gives none, which the owner chose on 5 October 2026: a call without a
/// seed keeps the same variants every time, so a result can be reproduced
/// from the code alone.
pub const DEFAULT_RANDOM_FILTER_SEED: u64 = 42;

/// SplitMix64, of Steele, Lea and Flood (2014): a state of 64 bits that
/// starts as the seed, and one draw of 64 bits at each call.
///
/// Each draw adds 0x9E3779B97F4A7C15 to the state and mixes the sum into the
/// number it gives, all of it wrapping at 2^64, as `splitmix64.c` of
/// Sebastiano Vigna and `nextLong` of `java.util.SplittableRandom` do.
#[derive(Debug, Clone, Copy)]
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// The generator with its state at `seed`.
    fn new(seed: u64) -> SplitMix64 {
        SplitMix64 { state: seed }
    }

    /// The next draw, a whole number from 0 to 2^64 - 1.
    fn next_draw(&mut self) -> u64 {
        // The wrapping at 2^64 is the arithmetic of the generator, and not
        // an overflow: it is how Java and `splitmix64.c` compute it.
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut mixed = self.state;
        mixed = (mixed ^ (mixed >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        mixed ^ (mixed >> 31)
    }

    /// The next number from 0 to 1, below 1: the top 53 bits of the next
    /// draw divided by 2^53, as `nextDouble` of Java gives it.
    #[expect(
        clippy::cast_precision_loss,
        reason = "the top 53 bits of a draw are below 2^53, which an f64 holds exactly"
    )]
    fn next_number(&mut self) -> f64 {
        // 2^-53, exact in an f64: the product is the division by 2^53.
        const TWO_TO_THE_MINUS_53: f64 = 1.0 / 9_007_199_254_740_992.0;
        (self.next_draw() >> 11) as f64 * TWO_TO_THE_MINUS_53
    }
}

/// The filter of one pass that keeps each variant when the number drawn for
/// it is below the keep rate, and counts what it was given and kept.
///
/// Its generator starts at the seed when it is built and draws one number
/// for each variant of every block it is given, in order, so the variants
/// it keeps are a function of the seed and of the place of each variant
/// among the variants it was given, whatever the size of the blocks. Every
/// pass builds its own, so every pass over the same steps keeps the same
/// variants.
#[derive(Debug)]
pub struct RandomFilter {
    keep_rate: f64,
    seed: u64,
    generator: SplitMix64,
    stats: FilteringStats,
}

impl RandomFilter {
    /// The filter that keeps each variant with the probability `keep_rate`,
    /// with its generator at `seed` and both its counts at 0. At a
    /// `keep_rate` of 1 every variant is kept and at 0 none is.
    ///
    /// # Errors
    ///
    /// When `keep_rate` is NaN, below 0 or above 1, with an error that names
    /// the argument and the value.
    pub fn new(keep_rate: f64, seed: u64) -> Result<RandomFilter> {
        // A NaN is in no range, so this one comparison refuses the three
        // keep rates that are not a number from 0 to 1.
        if !(0.0..=1.0).contains(&keep_rate) {
            return Err(Error::RandomFilterKeepRateOutOfRange { keep_rate });
        }
        Ok(RandomFilter {
            keep_rate,
            seed,
            generator: SplitMix64::new(seed),
            stats: FilteringStats::default(),
        })
    }

    /// The probability with which each variant is kept, from 0 to 1.
    #[must_use]
    pub fn keep_rate(&self) -> f64 {
        self.keep_rate
    }

    /// The seed the generator started at.
    #[must_use]
    pub fn seed(&self) -> u64 {
        self.seed
    }

    /// It draws one number for each variant of the block, in order, keeps
    /// the variants whose number is below the keep rate, in their order, and
    /// compacts the block in place; the variants of the block are added to
    /// the counts, and the ones that stayed to the ones kept. It reads no
    /// field of a variant, so a block without the genotypes is filtered as
    /// one with them. A block of no variants is left as it is and draws no
    /// number.
    ///
    /// # Errors
    ///
    /// When the arrays of the block are not of the size the block states,
    /// which [`Block::check`] finds. Then the block is as it was, nothing
    /// is added to the counts and no number is drawn, so the next block
    /// gets the numbers this one would have.
    pub fn filter_block(&mut self, block: &mut Block) -> Result<()> {
        block.check()?;
        let processed = block.num_vars;
        if processed == 0 {
            return Ok(());
        }
        // The numbers are drawn from a copy of the generator, which takes
        // its place once the block was compacted, so that a block refused
        // draws no number.
        let mut generator = self.generator;
        let keep_rate = self.keep_rate;
        let keep: Vec<bool> = (0..processed)
            .map(|_| generator.next_number() < keep_rate)
            .collect();
        let kept = keep.iter().filter(|keep_it| **keep_it).count();
        block.retain_vars(&keep)?;
        self.generator = generator;
        // A `usize` is 64 bits on the targets popnei builds natively for
        // and 32 in wasm, so every one of them is a `u64` and neither
        // conversion takes the value it saturates at.
        self.stats.vars_processed = self
            .stats
            .vars_processed
            .saturating_add(u64::try_from(processed).unwrap_or(u64::MAX));
        self.stats.vars_kept = self
            .stats
            .vars_kept
            .saturating_add(u64::try_from(kept).unwrap_or(u64::MAX));
        Ok(())
    }

    /// How many variants it was given and how many it kept, over every
    /// block it has taken since it was built.
    #[must_use]
    pub fn stats(&self) -> FilteringStats {
        self.stats
    }
}

/// A reader that gives the variants of its source that a [`RandomFilter`]
/// keeps.
///
/// It takes a block of its source at whatever size it comes, keeps in it
/// the variants the filter keeps and gives it on, so the blocks it gives
/// are of uneven size. What it keeps from one block to the next is the
/// generator of its filter and the two counts, which
/// [`BlockReader::filtering_stats`] gives under `random` before those of
/// its source.
///
/// It keeps the contract of a reader of `docs/specs/block.md`: a block left
/// with no variant is not given and the next one is taken; after an error,
/// of its source or of its filter, it gives `None` at every call and does
/// not call its source again; and a source that gives a block of no
/// variants has a defect and is the error of that.
pub struct RandomlyFilteredReader<R: BlockReader> {
    reader: R,
    filter: RandomFilter,
    /// Whether the source has no more blocks or one of the two, the source
    /// or the filter, gave an error. After any of them there is no block.
    finished: bool,
}

impl<R: BlockReader> RandomlyFilteredReader<R> {
    /// The reader that gives the variants of `reader` that `filter` keeps.
    ///
    /// # Errors
    ///
    /// When `reader` holds a filter that keeps variants at random already,
    /// which its [`BlockReader::filtering_stats`] says: the second would
    /// keep a sample of the sample of the first. The error carries the keep
    /// rate and the seed of `filter`, and not those of the filter that is
    /// set, which the counts of a chain do not say.
    pub fn new(reader: R, filter: RandomFilter) -> Result<RandomlyFilteredReader<R>> {
        if reader
            .filtering_stats()
            .iter()
            .any(|(of_the_chain, _)| *of_the_chain == RANDOM_KIND)
        {
            return Err(Error::RandomFilterThatIsSet {
                keep_rate: filter.keep_rate(),
                seed: filter.seed(),
                that_is_set: None,
            });
        }
        Ok(RandomlyFilteredReader {
            reader,
            filter,
            finished: false,
        })
    }
}

impl<R: BlockReader> BlockReader for RandomlyFilteredReader<R> {
    /// The next block of the source with the variants the filter keeps kept
    /// in it, and the blocks that its filter emptied passed over.
    ///
    /// # Errors
    ///
    /// When the source fails; when a block of the source holds no variant,
    /// which no reader of popnei gives; and when a block of the source has
    /// arrays that are not of its size, which [`RandomFilter::filter_block`]
    /// refuses. After any of them there is no block and the source is not
    /// called again.
    fn next_block(&mut self) -> Result<Option<Block>> {
        if self.finished {
            return Ok(None);
        }
        loop {
            let mut block = match self.reader.next_block() {
                Ok(Some(block)) => block,
                Ok(None) => {
                    self.finished = true;
                    return Ok(None);
                }
                Err(error) => {
                    self.finished = true;
                    return Err(error);
                }
            };
            // A source that gives a block of no variants has a defect, and
            // it is not asked again: over a source that always gives one, a
            // reader that asked again would never come back.
            if block.num_vars == 0 {
                self.finished = true;
                return Err(Error::ReaderGaveABlockOfNoVariants);
            }
            if let Err(error) = self.filter.filter_block(&mut block) {
                self.finished = true;
                return Err(error);
            }
            // A block the filter emptied is not given: the next one is
            // taken, and the source says when there are no more.
            if block.num_vars > 0 {
                return Ok(Some(block));
            }
        }
    }

    fn individuals(&self) -> &[String] {
        self.reader.individuals()
    }

    fn ploidy(&self) -> usize {
        self.reader.ploidy()
    }

    /// The table of the source: a reader over another reader has none of
    /// its own.
    fn chroms(&self) -> &ChromTable {
        self.reader.chroms()
    }

    /// What the consumer asked for, passed on as it is: this filter reads
    /// no field of a variant.
    fn set_needs(&mut self, needs: Needs) {
        self.reader.set_needs(needs);
    }

    /// The counts of this filter, under `random`, and after them those of
    /// the filters between the source and its own source.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        let mut stats = vec![(RANDOM_KIND, self.filter.stats())];
        stats.extend(self.reader.filtering_stats());
        stats
    }

    fn header(&self) -> &SourceHeader {
        self.reader.header()
    }

    /// False, and the source is not offered the regions: the variants it
    /// passed over would draw no number, and the variants after them would
    /// get other numbers than the pass without the regions gives them.
    fn skip_outside(&mut self, _selection: RegionSelection) -> bool {
        false
    }

    /// 0, since this filter hands no regions to its source.
    fn num_skipped(&self) -> u64 {
        0
    }
}

impl<R: BlockReader> fmt::Debug for RandomlyFilteredReader<R> {
    /// The filter and where it has got to. The source is left out, so that
    /// a `RandomlyFilteredReader` over a reader that has no `Debug` has one.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("RandomlyFilteredReader")
            .field("filter", &self.filter)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests;
