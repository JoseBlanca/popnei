//! The three statistics of a file in one pass on their way between
//! TypeScript and the core: the distributions of the statistics of each
//! variant, the rates of each individual and the density of the variants
//! along the chromosomes, which `calcPerVarDistribs`,
//! `calcPerIndividualStats` and `calcVarDensity` give in three passes.
//!
//! The arguments of each statistic are those of its own call, checked by
//! the code of that call, and each result crosses as the result of its own
//! call does, built by the same functions of `stats.rs` and `density.rs`,
//! so a page draws one from here as it draws it from its own call. What
//! goes out is [`VariantsSummaryOfAPass`]: each of the three or nothing,
//! and the counts of the one pass, which each of the three also carries.
//!
//! "The three statistics of a file in one pass" of
//! `docs/specs/js_sources.md` has the design.

use wasm_bindgen::JsValue;
use wasm_bindgen::prelude::wasm_bindgen;

use popnei::block::BlockReader;
use popnei::stats::{SoFar, VariantsSummary, VariantsSummaryConfig};

use crate::density::{ArgumentsOfTheDensity, VarDensityOfAPass, density_config_of, windows_of};
use crate::errors::JsPopneiError;
use crate::source::{Consumer, OpenSource, PassCounts, TheResultSoFar, the_run_of};
use crate::stats::{
    ArgumentsOfThePass, DistribsBeforeThePass, DistribsOfThePass, PerIndividualStats,
    PerVarDistribs, TheResultSoFarAsked, distribs_of, rates_of, under_its_name,
};
use crate::steps::{Steps, chain_of};

/// Which of the three statistics a pass of [`variants_summary_of`] gives,
/// each with the arguments of its own call, as they crossed from
/// TypeScript.
pub(crate) struct TheStatisticsAsked {
    /// The distributions of the statistics of each variant, and nothing for
    /// none.
    pub(crate) per_var: Option<ArgumentsOfThePass>,
    /// Whether the rates of each individual are given.
    pub(crate) per_individual: bool,
    /// The density of the variants, and nothing for none.
    pub(crate) density: Option<ArgumentsOfTheDensity>,
}

/// The statistics of `asked` over one pass over `source`, through the steps
/// of `steps`.
///
/// The chain of readers of the pass is built here and stays here, lent to
/// the core, so that the counts of its filters are read when the pass is
/// over: the loop over the blocks is the core's.
///
/// The function of `so_far` is given the three over the variants read so
/// far, each built as its own call builds its result so far, with the
/// counts of the pass as they stand after the block.
///
/// # Errors
///
/// When `asked` asks for none of the three, which the package refuses
/// before the call; those of the distributions, of the rates and of the
/// density that their own calls give, any of which ends the pass with none
/// of the three; and the value the function of `so_far` threw.
pub(crate) fn variants_summary_of(
    source: &dyn OpenSource,
    steps: &Steps,
    asked: &TheStatisticsAsked,
    so_far: TheResultSoFarAsked,
) -> Result<VariantsSummaryOfAPass, JsPopneiError> {
    let per_var = asked
        .per_var
        .as_ref()
        .map(|per_var| DistribsBeforeThePass::of(source, per_var))
        .transpose()?;
    let density = asked.density.as_ref().map(density_config_of).transpose()?;
    the_run_of(source, &Consumer::VariantsSummary, |run| {
        let reader = source.reader(run, None)?;
        let mut chain = chain_of(reader, steps.steps())?;
        let (per_var, under) = match per_var {
            Some(per_var) => {
                let DistribsOfThePass {
                    config,
                    pop_names,
                    hist_bin_edges,
                } = per_var.over(chain.individuals())?;
                (Some(config), (pop_names, hist_bin_edges))
            }
            None => (None, (Vec::new(), Vec::new())),
        };
        let names = TheNamesOfThePass {
            pop_names: under.0,
            hist_bin_edges: under.1,
            // The rates come out in the order of the rows of the blocks,
            // which is the order of these names, as for their own call.
            individuals: chain.individuals().to_vec(),
        };
        let config = VariantsSummaryConfig {
            per_var,
            per_individual: asked.per_individual,
            density,
            filter_column: false,
        };
        let mut told = TheResultSoFar::from_now(so_far.told, so_far.every_seconds)?;
        let given = popnei::stats::calc_variants_summary(
            &mut *chain,
            &config,
            &mut |added_up: &dyn SoFar<VariantsSummary>| {
                let Some(told) = told.as_mut() else {
                    return Ok(());
                };
                told.after_a_block(run, || {
                    let counts = PassCounts::of_the_filters(
                        added_up.num_vars(),
                        steps.steps(),
                        &added_up.filtering_stats(),
                    );
                    let summary = added_up.result().map_err(under_its_name)?;
                    Ok(JsValue::from(summary_of(summary, &names, counts)?))
                })
            },
        )
        .map_err(under_its_name);
        let summary = TheResultSoFar::what_the_pass_gives(told, given)?;
        let counts = PassCounts::of(num_vars_of(&summary), steps.steps(), &*chain);
        summary_of(summary, &names, counts)
    })
}

/// What the results of a pass are given under that the core does not
/// carry: the names of the populations and the edges of the bins of the
/// distributions, empty when they were not asked for, and the names of the
/// individuals of the rates.
struct TheNamesOfThePass {
    pop_names: Vec<String>,
    hist_bin_edges: Vec<f64>,
    individuals: Vec<String>,
}

/// How many variants the pass gave, which each statistic of `summary`
/// counted the same, and 0 for a summary of none, which the core refuses.
fn num_vars_of(summary: &VariantsSummary) -> u64 {
    match (&summary.per_var, &summary.per_individual, &summary.density) {
        (Some(per_var), _, _) => per_var.num_vars,
        (None, Some(per_individual), _) => per_individual.num_vars(),
        (None, None, Some(density)) => density.num_vars(),
        (None, None, None) => 0,
    }
}

/// `summary`, the statistics the core gave over the variants of a pass or
/// of its first blocks, on their way to JavaScript, each built as its own
/// call builds it, under `names` and with `counts`, the counts of the pass.
///
/// # Errors
///
/// Those of [`distribs_of`], [`rates_of`] and [`windows_of`].
fn summary_of(
    summary: VariantsSummary,
    names: &TheNamesOfThePass,
    counts: PassCounts,
) -> Result<VariantsSummaryOfAPass, JsPopneiError> {
    let VariantsSummary {
        per_var,
        per_individual,
        density,
        filter_column: _,
    } = summary;
    Ok(VariantsSummaryOfAPass {
        per_var: per_var
            .map(|distribs| {
                distribs_of(
                    distribs,
                    names.pop_names.clone(),
                    names.hist_bin_edges.clone(),
                    counts.clone(),
                )
            })
            .transpose()?,
        per_individual: per_individual
            .map(|stats| rates_of(&stats, names.individuals.clone(), counts.clone()))
            .transpose()?,
        density: density
            .map(|density| windows_of(&density, counts.clone()))
            .transpose()?,
        counts,
    })
}

/// What one pass of the three statistics gives JavaScript: each of the
/// three that was asked for, as its own call gives it, and the counts of
/// the pass.
///
/// The package takes each of the three out once and frees it when it has
/// read it, as it frees the result of its own call, and frees this object
/// after.
#[wasm_bindgen]
pub struct VariantsSummaryOfAPass {
    per_var: Option<PerVarDistribs>,
    per_individual: Option<PerIndividualStats>,
    density: Option<VarDensityOfAPass>,
    counts: PassCounts,
}

#[wasm_bindgen]
impl VariantsSummaryOfAPass {
    /// The distributions of the statistics of each variant, and nothing
    /// when they were not asked for or were already taken.
    pub fn take_per_var(&mut self) -> Option<PerVarDistribs> {
        self.per_var.take()
    }

    /// The rates of each individual, and nothing when they were not asked
    /// for or were already taken.
    pub fn take_per_individual(&mut self) -> Option<PerIndividualStats> {
        self.per_individual.take()
    }

    /// The density of the variants, and nothing when it was not asked for
    /// or was already taken.
    pub fn take_density(&mut self) -> Option<VarDensityOfAPass> {
        self.density.take()
    }

    /// How many variants the pass gave, and what each filter of it was
    /// given and kept, the outermost filter first.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.counts.clone()
    }
}
