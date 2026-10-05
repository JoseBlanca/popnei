//! The filter of `docs/specs/filters.md` that keeps the first n variants
//! that reach it and then ends the pass, so that a calculation tried on the
//! start of a large file does not read the rest of it.
//!
//! [`FirstNReader`] is the filter, a reader over another reader;
//! [`stopped_early`] says from the counts of a pass whether it ended that
//! pass; and [`refuse_a_step_after_the_first_n`] refuses a step that takes
//! variants out after it, so that n is the number of variants a calculation
//! gets.

use std::fmt;

use crate::block::{Block, BlockReader, SourceHeader};
use crate::error::{Error, Result};
use crate::filters::{FilteringStats, PassStep, RegionSelection};
use crate::variant::{ChromTable, Needs};

/// The name under which the counts of the filter of the first n reach a
/// Python or a TypeScript user, and by which a chain of readers is asked
/// whether it holds one already.
pub(crate) const FIRST_N_KIND: &str = "first_n";

/// A reader that gives the first `num_vars` variants of its source and then
/// no more, without asking its source for another block.
///
/// It takes a block of its source at whatever size it comes. While the
/// variants it has given and the block together are at most `num_vars`, it
/// gives the block whole; otherwise it keeps the first rows of the block up
/// to `num_vars` and gives them. Once it has given `num_vars` variants,
/// whether the last of them ended a block or fell inside one, every later
/// call gives `None` and its source is not called again, so the source
/// reads no further than the block that held the last of them.
///
/// Its counts, which [`BlockReader::filtering_stats`] gives under
/// `first_n` before those of its source, are of the blocks it took: the
/// variants of those blocks as given, and the ones it gave as kept. With
/// blocks of 7 variants and `num_vars` of 10 it is given 14 and keeps 10.
///
/// It keeps the contract of a reader of `docs/specs/block.md`: after an
/// error of its source it gives `None` at every call and does not call its
/// source again, and a source that gives a block of no variants has a
/// defect and is the error of that.
pub struct FirstNReader<R: BlockReader> {
    reader: R,
    /// How many variants it gives at most.
    num_vars: u64,
    /// The variants of the blocks it took, and the ones it gave.
    stats: FilteringStats,
    /// Whether the source has no more blocks or gave an error. After either
    /// there is no block.
    finished: bool,
}

impl<R: BlockReader> FirstNReader<R> {
    /// The reader that gives the first `num_vars` variants of `reader`.
    ///
    /// # Errors
    ///
    /// When `num_vars` is 0, since a pass it ended would give no variant;
    /// and when `reader` holds a filter of the first n already, which its
    /// [`BlockReader::filtering_stats`] says: two of them keep the first of
    /// the smaller n.
    pub fn new(reader: R, num_vars: u64) -> Result<FirstNReader<R>> {
        if num_vars == 0 {
            return Err(Error::FirstNOfNoVariants);
        }
        if reader
            .filtering_stats()
            .iter()
            .any(|(of_the_chain, _)| *of_the_chain == FIRST_N_KIND)
        {
            return Err(Error::FilterOfAKindThatIsSet { kind: FIRST_N_KIND });
        }
        Ok(FirstNReader {
            reader,
            num_vars,
            stats: FilteringStats::default(),
            finished: false,
        })
    }
}

impl<R: BlockReader> BlockReader for FirstNReader<R> {
    /// The next block of the source, whole or with its first rows alone
    /// when the rest would be more than `num_vars` variants, and `None`
    /// without calling the source once `num_vars` variants were given.
    ///
    /// # Errors
    ///
    /// When the source fails; when a block of the source holds no variant,
    /// which no reader of popnei gives; and when the block to be cut is one
    /// whose arrays are not of its size, which [`Block::retain_vars`]
    /// refuses. After any of them there is no block and the source is not
    /// called again.
    fn next_block(&mut self) -> Result<Option<Block>> {
        // Once the n were given the source is not asked again, so a pass
        // ends at the block that held the n-th variant.
        if self.finished || self.stats.vars_kept >= self.num_vars {
            return Ok(None);
        }
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
        // A source that gives a block of no variants has a defect, and it
        // is not asked again: over a source that always gives one, a reader
        // that asked again would never come back.
        if block.num_vars == 0 {
            self.finished = true;
            return Err(Error::ReaderGaveABlockOfNoVariants);
        }
        // A `usize` is 64 bits on the targets popnei builds natively for and
        // 32 in wasm, so it is always a `u64`, and the value it saturates at
        // is never taken.
        let given = u64::try_from(block.num_vars).unwrap_or(u64::MAX);
        // `vars_kept` is below `num_vars`, which the test above says.
        let wanted = self.num_vars.saturating_sub(self.stats.vars_kept);
        if given > wanted {
            // `wanted` is below the variants of the block, which are a
            // `usize`, so it is one too.
            let rows = usize::try_from(wanted).unwrap_or(block.num_vars);
            let keep: Vec<bool> = (0..block.num_vars).map(|row| row < rows).collect();
            if let Err(error) = block.retain_vars(&keep) {
                self.finished = true;
                return Err(error);
            }
        }
        self.stats.vars_processed = self.stats.vars_processed.saturating_add(given);
        self.stats.vars_kept = self.stats.vars_kept.saturating_add(given.min(wanted));
        Ok(Some(block))
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

    /// The counts of this filter, under `first_n`, and after them those of
    /// the filters between the source and its own source.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        let mut stats = vec![(FIRST_N_KIND, self.stats)];
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

impl<R: BlockReader> fmt::Debug for FirstNReader<R> {
    /// How many it gives and where it has got to. The source is left out,
    /// so that a `FirstNReader` over a reader that has no `Debug` has one.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("FirstNReader")
            .field("num_vars", &self.num_vars)
            .field("stats", &self.stats)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

/// Whether a filter of the first n among `steps` gave its n variants and so
/// ended the pass whose counts are `filtering`, as
/// [`BlockReader::filtering_stats`] of the outermost reader of the chain
/// gives them: true when `steps` has a [`PassStep::FirstN`] of n and the
/// counts under `first_n` have kept n.
///
/// It is true also when the source held nothing after the n-th variant,
/// since the filter does not read on to find out, and false when the
/// filter was given fewer than n, which is a pass that read the whole
/// source. The counts of the other filters of a pass that this ended are of
/// the part of the source that was read, which is what a user is told by
/// it. `steps` are those the pass was built from, which a binding crate
/// keeps from the start of the pass, and not those of the `Variants` when
/// the counts are read: a step added while a pass runs is not in it.
#[must_use]
pub fn stopped_early(steps: &[PassStep], filtering: &[(&'static str, FilteringStats)]) -> bool {
    let Some(num_vars) = steps.iter().find_map(|step| match step {
        PassStep::FirstN(num_vars) => Some(*num_vars),
        PassStep::VarFilter(_) | PassStep::KeepIndividuals(_) | PassStep::Regions(_) => None,
    }) else {
        return false;
    };
    filtering
        .iter()
        .any(|(kind, counts)| *kind == FIRST_N_KIND && counts.vars_kept == num_vars)
}

/// The step of the filter of the first `num_vars` variants, which both
/// binding crates build when a user adds the filter, so that a `num_vars`
/// of 0 is refused at that call and not at the next pass.
///
/// # Errors
///
/// When `num_vars` is 0, with a message that names `num_vars`: a pass the
/// filter ended would give no variant, as [`FirstNReader::new`] says.
pub fn first_n_step(num_vars: u64) -> Result<PassStep> {
    if num_vars == 0 {
        return Err(Error::FirstNOfNoVariants);
    }
    Ok(PassStep::FirstN(num_vars))
}

/// The error of a step that takes variants out, `new`, after a filter of
/// the first n among `set`, the steps that are set already.
///
/// With the filter of the first n set, n is the number of variants every
/// calculation gets, and a filter after it would leave fewer; the owner
/// decided on 5 October 2026 that it is refused. Both binding crates and
/// [`chain_of`](crate::filters::chain_of) call it through
/// [`refuse_a_step`](crate::filters::refuse_a_step), after the refusal of a
/// second filter of a kind.
///
/// # Errors
///
/// When `set` holds a [`PassStep::FirstN`] and `new` is a threshold filter,
/// the filter by linkage disequilibrium or a filter by regions of either
/// kind, with the kind of `new`. The filter of individuals takes no variant
/// out and is not refused; a second filter of the first n is the error of
/// [`refuse_a_second_filter_of_a_kind`](crate::filters::refuse_a_second_filter_of_a_kind),
/// which says what two of them do.
pub fn refuse_a_step_after_the_first_n(set: &[PassStep], new: &PassStep) -> Result<()> {
    let takes_variants_out = match new {
        PassStep::VarFilter(_) | PassStep::Regions(_) => true,
        PassStep::KeepIndividuals(_) | PassStep::FirstN(_) => false,
    };
    if takes_variants_out && set.iter().any(|step| matches!(step, PassStep::FirstN(_))) {
        return Err(Error::StepAfterTheFirstN { kind: new.kind() });
    }
    Ok(())
}

#[cfg(test)]
mod tests;
