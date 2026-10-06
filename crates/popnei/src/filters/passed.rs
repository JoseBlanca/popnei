//! The filter of `docs/specs/filters.md` that keeps the variants whose
//! FILTER column, in the VCF they were read from, was `PASS` or a dot, so
//! that how many variants failed their FILTER is in the counts of a pass
//! beside how many each other filter took out.
//!
//! [`PassedReader`] is the filter, a reader over another reader. It reads
//! the `passed` column of each block, which the VCF reader fills with the
//! function that `only_passed` uses, so that the two never disagree, and a
//! vars file of format 1.2 stores.

use std::fmt;

use crate::block::{Block, BlockReader, SourceHeader};
use crate::error::{Error, Result};
use crate::filters::{FilteringStats, RegionSelection};
use crate::variant::{ChromTable, Needs};

/// The name under which the counts of the filter of the variants that
/// passed reach a Python or a TypeScript user, and by which a chain of
/// readers is asked whether it holds one already.
pub(crate) const PASSED_KIND: &str = "passed";

/// A reader that gives the variants of its source whose `passed` is true,
/// those whose FILTER was `PASS` or a dot.
///
/// It takes a block of its source at whatever size it comes, keeps in place
/// the rows whose `passed` is true and gives the block on, so the blocks it
/// gives are of uneven size. Its counts, which
/// [`BlockReader::filtering_stats`] gives under `passed` before those of its
/// source, are the variants of the blocks it took and the ones it kept: over
/// `many.vcf` read with every variant, 500 given and 475 kept.
///
/// A source whose blocks have no `passed` column is refused at its first
/// block, and not taken as if every variant had passed: a vars file written
/// before format 1.2, or from a source that had no such column, may hold
/// variants that failed.
///
/// It keeps the contract of a reader of `docs/specs/block.md`: a block left
/// with no variant is not given and the next one is taken; after an error
/// it gives `None` at every call and does not call its source again; and a
/// source that gives a block of no variants has a defect and is the error
/// of that.
pub struct PassedReader<R: BlockReader> {
    reader: R,
    /// The variants of the blocks it took, and the ones it kept.
    stats: FilteringStats,
    /// Whether the source has no more blocks or an error was given. After
    /// either there is no block.
    finished: bool,
}

impl<R: BlockReader> PassedReader<R> {
    /// The reader that gives the variants of `reader` that passed their
    /// FILTER.
    ///
    /// # Errors
    ///
    /// When `reader` holds a filter of the variants that passed already,
    /// which its [`BlockReader::filtering_stats`] says: the second would
    /// keep the variants the first kept.
    pub fn new(reader: R) -> Result<PassedReader<R>> {
        if reader
            .filtering_stats()
            .iter()
            .any(|(of_the_chain, _)| *of_the_chain == PASSED_KIND)
        {
            return Err(Error::PassedFilterThatIsSet);
        }
        Ok(PassedReader {
            reader,
            stats: FilteringStats::default(),
            finished: false,
        })
    }
}

impl<R: BlockReader> PassedReader<R> {
    /// It keeps in place the rows of `block` whose `passed` is true, and
    /// adds the variants of the block to the counts, and the ones kept to
    /// the ones kept.
    ///
    /// # Errors
    ///
    /// When the block has no `passed` column, and when its arrays are not of
    /// its size. Then the block is as it was and nothing is added.
    fn keep_the_variants_that_passed(&mut self, block: &mut Block) -> Result<()> {
        block.check()?;
        // `retain_vars` compacts the column it is given as well, so it is
        // given a copy, one for the block and not one for each variant.
        let keep = block.passed.clone().ok_or(Error::PassedNotRecorded)?;
        let processed = block.num_vars;
        let kept = keep.iter().filter(|passed| **passed).count();
        block.retain_vars(&keep)?;
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
}

impl<R: BlockReader> BlockReader for PassedReader<R> {
    /// The next block of the source with the variants that passed kept in
    /// it, and the blocks in which none passed passed over.
    ///
    /// # Errors
    ///
    /// When the source fails; when a block of the source has no `passed`
    /// column, [`Error::PassedNotRecorded`]; when a block of the source holds
    /// no variant, which no reader of popnei gives; and when a block of the
    /// source has arrays that are not of its size, which [`Block::check`]
    /// finds. Then nothing is added to the counts, there is no block after
    /// it and the source is not called again.
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
            if let Err(error) = self.keep_the_variants_that_passed(&mut block) {
                self.finished = true;
                return Err(error);
            }
            // A block in which no variant passed is not given: the next one
            // is taken, and the source says when there are no more.
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

    /// The fields of the consumer and whether each variant passed, which
    /// the filter reads for every variant of every block: so the blocks
    /// this reader gives hold `passed` also when the consumer did not ask
    /// for it.
    fn set_needs(&mut self, needs: Needs) {
        self.reader.set_needs(needs.union(Needs::PASSED));
    }

    /// The counts of this filter, under `passed`, and after them those of
    /// the filters between the source and its own source.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        let mut stats = vec![(PASSED_KIND, self.stats)];
        stats.extend(self.reader.filtering_stats());
        stats
    }

    fn header(&self) -> &SourceHeader {
        self.reader.header()
    }

    /// False, and the source is not offered the regions: the variants it
    /// passed over would never reach this filter, and its counts would be
    /// of fewer variants than the source gave.
    fn skip_outside(&mut self, _selection: RegionSelection) -> bool {
        false
    }

    /// 0, since this filter hands no regions to its source.
    fn num_skipped(&self) -> u64 {
        0
    }
}

impl<R: BlockReader> fmt::Debug for PassedReader<R> {
    /// Where it has got to. The source is left out, so that a
    /// `PassedReader` over a reader that has no `Debug` has one.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("PassedReader")
            .field("stats", &self.stats)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests;
