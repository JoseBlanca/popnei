//! The filters that keep the variants of a dataset that pass a threshold,
//! and the counts of what each one was given and kept.
//!
//! The three threshold filters, over the missing rate of a variant, over
//! its major allele frequency and over its observed heterozygosity, are
//! being written. What is here is [`FilteringStats`], the two counts of one
//! filter of one pass, which every reader gives for the filters between it
//! and its source through
//! [`BlockReader::filtering_stats`](crate::block::BlockReader::filtering_stats).
//!
//! `docs/specs/filters.md` has the design, and the row `filters` of section
//! 9 of `docs/architecture.md` where the module sits.

/// How many variants a filter was given and how many of them it kept, over
/// every block it has taken since it was built.
///
/// A filter belongs to one pass over the source, one reading of it from its
/// start, so these are the counts of that pass alone. `vars_kept` is at
/// most `vars_processed`. Read while the pass runs, they are of the
/// variants the filter has seen, which can be more than the ones the
/// consumer of the pass has got.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FilteringStats {
    /// The variants the filter was given.
    pub vars_processed: u64,
    /// Those of them that passed its threshold.
    pub vars_kept: u64,
}
