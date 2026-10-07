//! The three statistics that add up over the blocks of a pass, the
//! distributions of the statistics of each variant, the rates of each
//! individual and the density of the variants, given by one pass, with the
//! counts of the FILTER column beside them when they are asked for.
//!
//! "The result so far and the three in one pass" of "The Rust interface" of
//! `docs/specs/stats.md` has the design, and "The three statistics of a file
//! in one pass" of `docs/specs/js_sources.md` what the TypeScript package
//! gives with it.

use super::density::TheDensity;
use super::{
    AfterABlock, DistribsAddedUp, IndividualsCounted, PerIndividualStats, PerVarDistribs,
    PerVarDistribsConfig, SoFar, VarDensity, no_variant_in_the_pass,
};
use crate::block::{Block, BlockReader, with_one_block_ahead};
use crate::error::{Error, Result};
use crate::filters::FilteringStats;
use crate::phases::{Phase, timed};
use crate::variant::Needs;

/// Which of the three statistics a pass of [`calc_variants_summary`] gives,
/// each with the options of its own function, and whether it gives the
/// counts of the FILTER column; one of the four at least, or the pass is
/// refused with [`Error::VariantsSummaryOfNoStatistic`].
#[derive(Debug)]
pub struct VariantsSummaryConfig {
    /// The distributions of the statistics of each variant, with the
    /// options of [`calc_per_var_distribs`](super::calc_per_var_distribs),
    /// or `None` for none.
    pub per_var: Option<PerVarDistribsConfig>,
    /// Whether the rates of each individual of
    /// [`calc_per_individual_stats`](super::calc_per_individual_stats) are
    /// given.
    pub per_individual: bool,
    /// The density of the variants, with the options of
    /// [`calc_var_density`](super::calc_var_density), or `None` for none.
    pub density: Option<VarDensityConfig>,
    /// Whether the counts of the FILTER column, how many variants of the
    /// pass passed their FILTER and how many failed, are given.
    pub filter_column: bool,
}

/// The two arguments of [`calc_var_density`](super::calc_var_density)
/// beside its reader, for the density of a [`calc_variants_summary`].
#[derive(Debug)]
pub struct VarDensityConfig {
    /// The width of a window in base pairs, 1 or more.
    pub window_size: u64,
    /// The lengths of the chromosomes, which replace those of the header of
    /// the source when they are given.
    pub chrom_lengths: Option<Vec<(String, u64)>>,
}

/// What one pass of [`calc_variants_summary`] gives back: each of the three
/// statistics that was asked for, the same to the bit as its own function
/// gives over the same reader, the counts of the FILTER column when they
/// were asked for, and `None` for each that was not.
#[derive(Debug)]
pub struct VariantsSummary {
    /// The distributions of the statistics of each variant.
    pub per_var: Option<PerVarDistribs>,
    /// The rates of each individual.
    pub per_individual: Option<PerIndividualStats>,
    /// The density of the variants along the chromosomes.
    pub density: Option<VarDensity>,
    /// The counts of the FILTER column.
    pub filter_column: Option<FilterColumnCounts>,
}

/// Of the variants of a pass of [`calc_variants_summary`], how many passed
/// their FILTER, whose column in the VCF they were read from was `PASS` or
/// a dot, and how many failed; the two add up to the variants of the pass.
///
/// They are of the variants that reach the summary, after every step of
/// the pass, so with the filter of the variants that passed among the steps
/// `failed` is 0. Over `many.vcf` read with every variant they are 475 and
/// 25.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FilterColumnCounts {
    /// The variants whose FILTER was `PASS` or a dot.
    pub passed: u64,
    /// The variants whose FILTER was anything else.
    pub failed: u64,
}

/// The statistics of `config` over the variants `reader` gives, in one pass
/// over the source through the steps the variants carry, where their own
/// functions take one pass each.
///
/// Each statistic is added up from the blocks by the code its own function
/// runs, and only the pass is shared, so each is what its function gives
/// with the same options. The pass asks the reader for the union of what
/// the statistics asked for ask for: the genotypes when the distributions
/// or the rates are, the chromosome and the position when the density is,
/// and whether each variant passed when the counts of the FILTER column
/// are, so a pass of the density alone reads no genotype.
///
/// The counts of the FILTER column are of the variants that reach the
/// summary, after every step of the chain, and take none of them out; with
/// the counts alone the variants of the pass are those that passed and
/// those that failed.
///
/// `reader` is the outermost reader of the chain of the pass, lent and not
/// taken, so that whoever built the chain reads the counts of its filters
/// from it when this returns. The blocks are read one block ahead on a
/// thread of their own natively, as for the distributions and the rates,
/// and one after another in wasm.
///
/// `after_a_block` is called after each block the pass adds, the last one
/// too, with the statistics over the variants read so far, each what the
/// result so far of its own function would be.
///
/// # Errors
///
/// - [`Error::VariantsSummaryOfNoStatistic`] for a `config` that asks for
///   none of the four, before anything is read.
/// - Those each statistic asked for refuses before its pass, in the order
///   of the fields of `config`: a `poly_threshold` that is not a number from
///   0 to 1, a window of 0 base pairs and lengths that
///   [`calc_var_density`](super::calc_var_density) refuses, and
///   [`Error::FilterColumnNotRecorded`] for the counts of the FILTER column
///   over a source whose header says it keeps no record of whether its
///   variants passed.
/// - What any of them refuses in a block, which ends the pass with no
///   statistic given: when one block is refused by more than one, the error
///   is the one of the distributions, then of the rates, then of the
///   density, then of the counts, which refuse a block without the column
///   with [`Error::FieldsNotInTheBlock`] and [`Needs::PASSED`].
/// - [`Error::PassGaveNoVariant`] for a pass that gave no variant, as each
///   of the three gives it.
/// - What the reader fails with, and the error `after_a_block` returns,
///   which ends the pass there.
pub fn calc_variants_summary<R: BlockReader + ?Sized>(
    reader: &mut R,
    config: &VariantsSummaryConfig,
    after_a_block: AfterABlock<'_, VariantsSummary>,
) -> Result<VariantsSummary> {
    if config.per_var.is_none()
        && !config.per_individual
        && config.density.is_none()
        && !config.filter_column
    {
        return Err(Error::VariantsSummaryOfNoStatistic);
    }
    let mut added_up = AddedUp {
        per_var: config
            .per_var
            .as_ref()
            .map(|per_var| DistribsAddedUp::before_the_pass(per_var, &*reader))
            .transpose()?,
        per_individual: config
            .per_individual
            .then(|| IndividualsCounted::before_the_pass(&*reader)),
        density: config
            .density
            .as_ref()
            .map(|density| {
                TheDensity::before_the_pass(
                    &*reader,
                    density.window_size,
                    density.chrom_lengths.as_deref(),
                )
            })
            .transpose()?,
        filter_column: config
            .filter_column
            .then(|| FilterColumnCounted::before_the_pass(&*reader))
            .transpose()?,
    };
    let mut needs = Needs::empty();
    if added_up.per_var.is_some() || added_up.per_individual.is_some() {
        needs |= Needs::GTS;
    }
    if added_up.density.is_some() {
        needs |= Needs::CHROM_POS;
    }
    if added_up.filter_column.is_some() {
        needs |= Needs::PASSED;
    }
    reader.set_needs(needs);
    // The blocks are read on a thread of its own, one block ahead, as
    // `calc_per_var_distribs` reads them and for its reason; in wasm there
    // is no thread and they come one after another. The chain is lent
    // through a reborrow of its own, as there.
    let mut lent = &mut *reader;
    with_one_block_ahead(&mut lent, |blocks| {
        while let Some(block) = timed(Phase::NextBlock, || blocks.next_block())? {
            timed(Phase::Work, || added_up.add_the_block(&block, &*blocks))?;
            after_a_block(&SummarySoFar {
                added_up: &added_up,
                chain: &*blocks,
            })?;
        }
        Ok(())
    })?;
    if added_up.num_vars() == 0 {
        return Err(no_variant_in_the_pass(reader));
    }
    Ok(added_up.into_the_result())
}

/// What a pass of [`calc_variants_summary`] adds up from block to
/// block: what the function of each statistic that was asked for adds up.
struct AddedUp<'config> {
    per_var: Option<DistribsAddedUp<'config>>,
    per_individual: Option<IndividualsCounted>,
    density: Option<TheDensity>,
    filter_column: Option<FilterColumnCounted>,
}

impl AddedUp<'_> {
    /// It adds `block` into the totals of each statistic, with the names of
    /// the chromosomes of `chain`, which has just given it.
    ///
    /// # Errors
    ///
    /// The first error of the four, in the order of the fields.
    fn add_the_block(&mut self, block: &Block, chain: &dyn BlockReader) -> Result<()> {
        if let Some(per_var) = &mut self.per_var {
            per_var.add_the_block(block)?;
        }
        if let Some(per_individual) = &mut self.per_individual {
            per_individual.add_the_block(block)?;
        }
        if let Some(density) = &mut self.density {
            density.count_the_block(block, chain.chroms())?;
        }
        if let Some(filter_column) = &mut self.filter_column {
            filter_column.count_the_block(block)?;
        }
        Ok(())
    }

    /// The variants added so far, which every statistic asked for has
    /// counted the same: with the counts of the FILTER column alone, those
    /// that passed and those that failed.
    fn num_vars(&self) -> u64 {
        match (
            &self.per_var,
            &self.per_individual,
            &self.density,
            &self.filter_column,
        ) {
            (Some(per_var), _, _, _) => per_var.num_vars,
            (None, Some(per_individual), _, _) => per_individual.num_vars,
            (None, None, Some(density), _) => density.num_vars,
            (None, None, None, Some(filter_column)) => filter_column.num_vars(),
            (None, None, None, None) => 0,
        }
    }

    /// Each statistic when the pass ends, which takes what the rates and
    /// the density counted; the distributions are built from their totals
    /// as they are for a result so far.
    fn into_the_result(self) -> VariantsSummary {
        VariantsSummary {
            per_var: self.per_var.as_ref().map(DistribsAddedUp::the_result),
            per_individual: self.per_individual.map(IndividualsCounted::into_the_result),
            density: self.density.map(TheDensity::into_the_result),
            filter_column: self.filter_column.map(|counted| counted.counts),
        }
    }

    /// Each statistic over the variants added so far.
    fn the_result(&self) -> VariantsSummary {
        VariantsSummary {
            per_var: self.per_var.as_ref().map(DistribsAddedUp::the_result),
            per_individual: self
                .per_individual
                .as_ref()
                .map(IndividualsCounted::the_result),
            density: self.density.as_ref().map(TheDensity::the_result),
            filter_column: self.filter_column.as_ref().map(|counted| counted.counts),
        }
    }
}

/// The counts of the FILTER column a pass of [`calc_variants_summary`]
/// adds up, from the `passed` column of each block.
struct FilterColumnCounted {
    counts: FilterColumnCounts,
}

impl FilterColumnCounted {
    /// The counts before the first block of a pass over `reader`, both 0.
    ///
    /// # Errors
    ///
    /// [`Error::FilterColumnNotRecorded`] when the header of `reader` says
    /// that its source keeps no record of whether its variants passed.
    fn before_the_pass<R: BlockReader + ?Sized>(reader: &R) -> Result<FilterColumnCounted> {
        if !reader.header().keeps_passed {
            return Err(Error::FilterColumnNotRecorded);
        }
        Ok(FilterColumnCounted {
            counts: FilterColumnCounts {
                passed: 0,
                failed: 0,
            },
        })
    }

    /// It adds the variants of `block` that passed to those that passed,
    /// and the others to those that failed.
    ///
    /// # Errors
    ///
    /// When the arrays of the block are not of its size, which
    /// [`Block::check`] finds; [`Error::ReaderGaveABlockOfNoVariants`] for a
    /// block of none; and [`Error::FieldsNotInTheBlock`] with
    /// [`Needs::PASSED`] for a block without the column. Then nothing is
    /// added.
    fn count_the_block(&mut self, block: &Block) -> Result<()> {
        // A column of another length than the block would count variants
        // that are not in it, or leave some of it out.
        block.check()?;
        // As for the three: every reader of popnei gives one variant at
        // least in a block, so a block of none is a defect of the reader.
        if block.num_vars == 0 {
            return Err(Error::ReaderGaveABlockOfNoVariants);
        }
        let Some(passed) = &block.passed else {
            return Err(Error::FieldsNotInTheBlock {
                fields: Needs::PASSED,
            });
        };
        let num_passed = passed.iter().filter(|passed| **passed).count();
        let num_failed = passed.len().saturating_sub(num_passed);
        // A `usize` is 64 bits on the targets popnei builds natively for
        // and 32 in wasm, so each count is a `u64`; and a pass of more than
        // 18446744073709551615 variants reads more rows than any source
        // holds, so neither total, nor the two added, reaches that.
        let as_a_count = |count: usize| u64::try_from(count).unwrap_or(u64::MAX);
        self.counts.passed = self.counts.passed.saturating_add(as_a_count(num_passed));
        self.counts.failed = self.counts.failed.saturating_add(as_a_count(num_failed));
        Ok(())
    }

    /// The variants counted so far, those that passed and those that
    /// failed.
    fn num_vars(&self) -> u64 {
        // No pass reaches 18446744073709551615 variants, as above.
        self.counts.passed.saturating_add(self.counts.failed)
    }
}

/// What a pass of [`calc_variants_summary`] has added up after a block,
/// which builds each statistic from its totals when asked.
struct SummarySoFar<'pass> {
    added_up: &'pass AddedUp<'pass>,
    /// The reader the pass takes its blocks from, which answers with the
    /// counts of the filters as they were when it gave the last block.
    chain: &'pass dyn BlockReader,
}

impl SoFar<VariantsSummary> for SummarySoFar<'_> {
    fn num_vars(&self) -> u64 {
        self.added_up.num_vars()
    }

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        self.chain.filtering_stats()
    }

    fn result(&self) -> Result<VariantsSummary> {
        Ok(self.added_up.the_result())
    }
}

#[cfg(test)]
mod tests;
