//! The three statistics that add up over the blocks of a pass, the
//! distributions of the statistics of each variant, the rates of each
//! individual and the density of the variants, given by one pass.
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
/// each with the options of its own function; one at least, or the pass is
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
/// gives over the same reader, and `None` for each that was not.
#[derive(Debug)]
pub struct VariantsSummary {
    /// The distributions of the statistics of each variant.
    pub per_var: Option<PerVarDistribs>,
    /// The rates of each individual.
    pub per_individual: Option<PerIndividualStats>,
    /// The density of the variants along the chromosomes.
    pub density: Option<VarDensity>,
}

/// The statistics of `config` over the variants `reader` gives, in one pass
/// over the source through the steps the variants carry, where their own
/// functions take one pass each.
///
/// Each statistic is added up from the blocks by the code its own function
/// runs, and only the pass is shared, so each is what its function gives
/// with the same options. The pass asks the reader for the union of what
/// the statistics asked for ask for: the genotypes when the distributions
/// or the rates are, and the chromosome and the position when the density
/// is, so a pass of the density alone reads no genotype.
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
///   none of the three, before anything is read.
/// - Those each statistic asked for refuses before its pass: a
///   `poly_threshold` that is not a number from 0 to 1, a window of 0 base
///   pairs and lengths that [`calc_var_density`](super::calc_var_density)
///   refuses.
/// - What any of them refuses in a block, which ends the pass with no
///   statistic given: when one block is refused by more than one, the error
///   is the one of the distributions, then of the rates, then of the
///   density.
/// - [`Error::PassGaveNoVariant`] for a pass that gave no variant, as each
///   of the three gives it.
/// - What the reader fails with, and the error `after_a_block` returns,
///   which ends the pass there.
pub fn calc_variants_summary<R: BlockReader + ?Sized>(
    reader: &mut R,
    config: &VariantsSummaryConfig,
    after_a_block: AfterABlock<'_, VariantsSummary>,
) -> Result<VariantsSummary> {
    if config.per_var.is_none() && !config.per_individual && config.density.is_none() {
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
    };
    let mut needs = Needs::empty();
    if added_up.per_var.is_some() || added_up.per_individual.is_some() {
        needs |= Needs::GTS;
    }
    if added_up.density.is_some() {
        needs |= Needs::CHROM_POS;
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
}

impl AddedUp<'_> {
    /// It adds `block` into the totals of each statistic, with the names of
    /// the chromosomes of `chain`, which has just given it.
    ///
    /// # Errors
    ///
    /// The first error of the three, in the order of the fields.
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
        Ok(())
    }

    /// The variants added so far, which every statistic asked for has
    /// counted the same.
    fn num_vars(&self) -> u64 {
        match (&self.per_var, &self.per_individual, &self.density) {
            (Some(per_var), _, _) => per_var.num_vars,
            (None, Some(per_individual), _) => per_individual.num_vars,
            (None, None, Some(density)) => density.num_vars,
            (None, None, None) => 0,
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
        }
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
