//! How much variety each population holds, on its way between TypeScript
//! and the core.
//!
//! The calculation is the core's,
//! [`popnei::diversity::calc_pop_diversity`], and what this module does is
//! what `pop_dists.rs` does for the distances between populations: it turns
//! the arguments a TypeScript user wrote into what the core takes, builds
//! the chain of readers of the pass from the steps of the `Variants`, lends
//! it to the core, and reads the counts of the filters from that chain when
//! the call is over, since no block of the pass reaches this crate.
//!
//! The populations are looked up against the individuals of the pass, which
//! only the pass knows, and `Pops` of the `stats` module does it, so that a
//! `pops` argument is read the same way here as in `calcPerVarDistribs`.
//! They cross flat, the names of the individuals of every one of them in one
//! array, as those of the statistics do: an array of arrays is not one of
//! the types wasm-bindgen carries.
//!
//! What goes out is one array per column of the result, one value per
//! population in the order the user named them, and the TypeScript package
//! builds the result object out of them. The core gives the totals and this
//! crate gives no mean and no ratio: each of those is a total over a count
//! of variants that the result carries beside it, and the package divides,
//! as the Python one does.
//!
//! `docs/specs/diversity.md` has the design.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::block::BlockReader;
use popnei::diversity::{
    DiversityOptions, DiversityStats, PopDiversity, calc_pop_diversity as diversity_of_the_pops,
};
use popnei::stats::Pops;

use crate::errors::JsPopneiError;
use crate::source::{OpenSource, PassCounts};
use crate::stats::{PopsGiven, pops_of_the_arrays};
use crate::steps::{Steps, chain_of};

/// The arguments of one pass, as they crossed from TypeScript.
///
/// The package has checked that each of them is of the type the core takes,
/// since a number of JavaScript reaches a whole number of the core as 32
/// bits with no error; what is left is what the core says of them, an
/// unknown name of a statistic and a draw of fewer than two alleles among
/// it.
pub(crate) struct ArgumentsOfTheDiversity {
    /// The statistics to calculate, under the names of
    /// [`DiversityStats::NAMES`].
    pub(crate) stats: Vec<String>,
    /// The name of each population, in the order the user gave them, and
    /// `None` when they named none, which is one population of every
    /// individual of the pass.
    pub(crate) pop_names: Option<Vec<String>>,
    /// The names of the individuals of every population, the ones of the
    /// first population first.
    pub(crate) pop_individuals: Vec<String>,
    /// How many individuals each population of `pop_names` holds, which
    /// cuts `pop_individuals` into the names of each of them.
    pub(crate) num_individuals_per_pop: Vec<u32>,
    /// How many called alleles every population is brought down to, and
    /// nothing for a pass that takes no draw.
    pub(crate) num_called_alleles: Option<u32>,
    /// How many called genotypes a population needs at a variant for the
    /// variant to count for it.
    pub(crate) min_num_individuals: u32,
}

/// The diversity of every population of `asked` over one pass over `source`,
/// through the steps of `steps`.
///
/// The chain of readers of the pass stays here, lent to the core, so that
/// the counts of its filters are read when the pass is over: the loop over
/// the blocks is the core's.
///
/// The reader is opened at the size of its own blocks and no `Reblock` is
/// put over the chain. The four counts are integer totals over the variants
/// and are the same numbers whatever the blocks are. F_IS is not: it is the
/// ratio of two sums of float64 that the core adds in the order of the
/// variants within each chunk of a block, so where the blocks are cut moves
/// its last bits. On the panel of `tests/reference/stats/` it moves by
/// 6.7e-16, 5.2e-14 of the value, between the whole file in one block and
/// blocks of 7 variants, measured on 24 September 2026, which is inside the
/// 1e-12 of the value that `docs/specs/diversity.md` compares within and is
/// not nothing.
///
/// # Errors
///
/// When a name of `asked.stats` is of none of the five statistics; when the
/// folded spectrum is asked for with no draw and when the draw is of fewer
/// than two alleles; when a population names an individual the pass does not
/// give, names one twice or names none, and when `pops` holds no population;
/// when the source cannot be read, a wrong line of a VCF among the causes;
/// when the draw is of more alleles than the dataset holds gene copies; when
/// the pass gives no variant; when a variant holds more alleles than a
/// count of them holds; when the pass counted the variants of no population or
/// gave the populations spectra of different numbers of bins, which are both a
/// defect of popnei; and when a count is above what a JavaScript array of
/// counts holds.
pub(crate) fn pop_diversity_of(
    source: &dyn OpenSource,
    steps: &Steps,
    asked: &ArgumentsOfTheDiversity,
) -> Result<PopDiversityOfAPass, JsPopneiError> {
    let options = DiversityOptions {
        stats: the_stats(&asked.stats)?,
        num_called_alleles: asked.num_called_alleles,
        min_num_individuals: asked.min_num_individuals,
    };
    let named = the_pops_given(asked)?;
    let reader = source.reader(None)?;
    let mut chain = chain_of(reader, steps.steps())?;
    let pops = match named {
        Some(named) => Pops::from_names(&named, chain.individuals())?,
        None => Pops::all(chain.individuals().len()),
    };
    let pop_names: Vec<String> = (0..pops.len())
        .map(|pop| pops.name(pop).to_owned())
        .collect();
    let of_each_pop: Vec<&[usize]> = (0..pops.len()).map(|pop| pops.individuals(pop)).collect();
    let diversity = diversity_of_the_pops(&mut *chain, &of_each_pop, &options)?;
    let counts = PassCounts::of(diversity.num_vars_of_the_pass(), &chain.filtering_stats());
    let num_pops = diversity.num_pops();
    // Every pass counts the variants of each of its populations, whatever
    // statistics it was asked for, so a result with no such count is a
    // defect and not a statistic nobody asked for.
    let Some(num_vars_with_data) = of_every_pop(num_pops, |pop| diversity.num_vars(pop)) else {
        return Err(JsPopneiError::Broken(format!(
            "the pass counted the variants of no population, and it was over {num_pops} \
             of them"
        )));
    };
    let Some(num_vars_in_draw) = of_every_pop(num_pops, |pop| diversity.num_vars_in_draw(pop))
    else {
        return Err(JsPopneiError::Broken(format!(
            "the pass counted the variants in the draw for no population, and it was \
             over {num_pops} of them"
        )));
    };
    // The package reads a statistic as the total and the standardized value
    // together and refuses a result that holds one of the two, so the
    // second is there exactly when the first is: the two accessors of the core
    // answer with nothing for the same reason, that nobody asked for the
    // statistic.
    let num_alleles_total = total_for_javascript(
        of_every_pop(num_pops, |pop| diversity.num_alleles(pop)),
        "the alleles one population called",
    )?;
    let private_alleles_total = total_for_javascript(
        of_every_pop(num_pops, |pop| diversity.private_alleles(pop)),
        "the private alleles of one population",
    )?;
    let variable_vars_total = total_for_javascript(
        of_every_pop(num_pops, |pop| diversity.num_variable_vars(pop)),
        "the variable variants of one population",
    )?;
    let spectrum = the_spectrum_of_every_pop(num_pops, &diversity)?;
    Ok(PopDiversityOfAPass {
        pop_names: Some(pop_names),
        num_vars_with_data: Some(for_javascript_each(
            num_vars_with_data,
            "the variants of one population",
        )?),
        num_vars_in_draw: Some(for_javascript_each(
            num_vars_in_draw,
            "the variants in the draw for one population",
        )?),
        num_vars_every_pop: for_javascript(
            diversity.num_vars_every_pop(),
            "the variants that counted for every population",
        )?,
        num_vars_every_pop_in_draw: for_javascript(
            diversity.num_vars_every_pop_in_draw(),
            "the variants every population reached the draw at",
        )?,
        num_alleles_in_draw: of_every_pop(num_pops, |pop| diversity.num_alleles_in_draw(pop)),
        num_alleles_total,
        private_alleles_in_draw: of_every_pop(num_pops, |pop| {
            diversity.private_alleles_in_draw(pop)
        }),
        private_alleles_total,
        variable_vars_in_draw: of_every_pop(num_pops, |pop| {
            diversity.variable_vars_ratio_in_draw(pop)
        }),
        variable_vars_total,
        num_sfs_bins: spectrum.num_bins,
        folded_sfs: spectrum.bins_of_every_pop,
        fis: of_every_pop(num_pops, |pop| diversity.fis(pop)),
        counts,
    })
}

/// The folded spectrum of every population as it crosses to JavaScript: the
/// bins of the first population, then those of the second, in one array, with
/// how many bins each of them has.
///
/// A `Vec<Vec<f64>>` is not one of the types wasm-bindgen carries, so the
/// values cross flat, as the names of the individuals of the populations do,
/// and the bin count is what cuts them: the package hands a user one
/// `Float64Array` per population, so the bins of one population lie together.
struct SpectrumOfEveryPop {
    /// The bins of each population one after another, and nothing when nobody
    /// asked for the spectrum.
    bins_of_every_pop: Option<Vec<f64>>,
    /// How many bins each population has, which is
    /// `num_called_alleles / 2 + 1`, and 0 when nobody asked for the
    /// spectrum.
    num_bins: u32,
}

/// The folded spectrum of every population of `diversity`, flat.
///
/// # Errors
///
/// When the populations do not all have the same bins, and when a population
/// has more bins than a count of them holds, which are both a defect of
/// popnei: the bins of a spectrum are the counts of the rarer allele from 0 to
/// `num_called_alleles` over 2, one draw size for the whole call, and the core
/// refuses a draw of more alleles than the dataset holds gene copies.
fn the_spectrum_of_every_pop(
    num_pops: usize,
    diversity: &PopDiversity,
) -> Result<SpectrumOfEveryPop, JsPopneiError> {
    let Some(of_each_pop) = of_every_pop(num_pops, |pop| diversity.folded_sfs(pop)) else {
        return Ok(SpectrumOfEveryPop {
            bins_of_every_pop: None,
            num_bins: 0,
        });
    };
    let num_bins = of_each_pop.first().map_or(0, |bins| bins.len());
    if of_each_pop.iter().any(|bins| bins.len() != num_bins) {
        return Err(JsPopneiError::Broken(format!(
            "the pass gave the folded spectrum of {num_pops} populations and not the \
             same {num_bins} bins for each of them"
        )));
    }
    Ok(SpectrumOfEveryPop {
        bins_of_every_pop: Some(of_each_pop.concat()),
        num_bins: u32::try_from(num_bins).map_err(|_| {
            JsPopneiError::NotInJavaScript(format!(
                "the folded spectrum of one population came to {num_bins} bins, more \
                 than a JavaScript array of counts holds"
            ))
        })?,
    })
}

/// One value of every population of the result, in their order, and `None`
/// when the statistic was not asked for.
///
/// Every accessor of the core gives `None` both for a population it has not
/// and for a statistic nobody asked for, and the populations here are
/// `0..num_pops` of the result itself, so `None` can only be the second.
fn of_every_pop<T>(num_pops: usize, of_the_pop: impl Fn(usize) -> Option<T>) -> Option<Vec<T>> {
    (0..num_pops).map(of_the_pop).collect()
}

/// The names of the statistics that need no draw, which is what a user who
/// names no statistic in `stats` asks for: the alleles a population called,
/// the private ones among them, the variants that vary in it and F_IS.
///
/// The TypeScript package asks for these when a call gives no `stats`. The
/// four are named in the core alone, and not there as well as in the Python
/// package, so a statistic that needs no draw is added in one place.
#[wasm_bindgen]
#[must_use]
pub fn diversity_stats_without_a_draw() -> Vec<String> {
    DiversityStats::WITHOUT_A_DRAW
        .names()
        .into_iter()
        .map(str::to_owned)
        .collect()
}

/// The statistics a user asked for, out of the names the TypeScript package
/// gives them, which are the core's.
///
/// # Errors
///
/// A name that is of no statistic, which a user reaches by writing one in
/// JavaScript: in TypeScript the five are a union of string literals. The
/// core refuses it and names the five, so the sentence a user reads is
/// written once for both languages.
fn the_stats(names: &[String]) -> Result<DiversityStats, JsPopneiError> {
    let mut asked_for = DiversityStats::empty();
    for name in names {
        asked_for |= DiversityStats::of_name(name)?;
    }
    Ok(asked_for)
}

/// The populations a user named, each with the names of its individuals, out
/// of the flat arrays they crossed in, and `None` when they named none,
/// which is one population of every individual of the pass.
///
/// The names are not looked up here: they are resolved against the
/// individuals the pass gives, which are those of the source after a filter
/// of individuals when the `Variants` has one, and only the pass knows them.
///
/// # Errors
///
/// When the arrays do not hold the individuals of every population, which is
/// a defect of the package: it is what cuts them.
fn the_pops_given(asked: &ArgumentsOfTheDiversity) -> Result<Option<PopsGiven>, JsPopneiError> {
    let Some(names) = asked.pop_names.as_ref() else {
        return Ok(None);
    };
    Ok(Some(pops_of_the_arrays(
        names,
        &asked.pop_individuals,
        &asked.num_individuals_per_pop,
    )?))
}

/// One total of every population as the arrays of a result hold them, or
/// `None` when nobody asked for the statistic.
///
/// # Errors
///
/// When a count is above what a JavaScript array of counts holds.
fn total_for_javascript(
    counts: Option<Vec<u64>>,
    what: &str,
) -> Result<Option<Vec<u32>>, JsPopneiError> {
    counts
        .map(|counts| for_javascript_each(counts, what))
        .transpose()
}

/// Every count of `counts` as the arrays of a result hold them, in their
/// order.
///
/// # Errors
///
/// Those of [`for_javascript`].
fn for_javascript_each(counts: Vec<u64>, what: &str) -> Result<Vec<u32>, JsPopneiError> {
    counts
        .into_iter()
        .map(|count| for_javascript(count, what))
        .collect()
}

/// `count` as the array of counts of JavaScript holds it.
///
/// The counts of a result cross as a `Uint32Array`, which is what
/// `docs/specs/diversity.md` gives them, and the core counts in 64 bits.
///
/// # Errors
///
/// When the count is above 4294967295. A count of variants never reaches it:
/// that is more variants than a file in the memory of a tab holds, which
/// addresses 4 GB, and a variant is a row of a file. A count of alleles can,
/// since it sums the alleles of every variant, and then what a user does is
/// filter the variants or read the same dataset from Python.
fn for_javascript(count: u64, what: &str) -> Result<u32, JsPopneiError> {
    u32::try_from(count).map_err(|_| {
        JsPopneiError::NotInJavaScript(format!(
            "{what} came to {count}, more than a JavaScript array of counts holds"
        ))
    })
}

/// What one pass of the diversity gives JavaScript: the names of the
/// populations in their order, the variants that counted for each of them
/// and those of them that reached the draw, the same two counts over the
/// populations together, the alleles called, the private ones and the
/// variable variants each as a total and a standardized value, the folded
/// spectrum of every population with how many bins each of them has, F_IS,
/// and the counts of the pass.
///
/// Every array leaves the memory of wasm the first time it is asked for, and
/// the call after that gives nothing: the package reads each of them once,
/// into the object a user holds, and frees this.
#[wasm_bindgen]
pub struct PopDiversityOfAPass {
    pop_names: Option<Vec<String>>,
    num_vars_with_data: Option<Vec<u32>>,
    num_vars_in_draw: Option<Vec<u32>>,
    num_vars_every_pop: u32,
    num_vars_every_pop_in_draw: u32,
    num_alleles_total: Option<Vec<u32>>,
    num_alleles_in_draw: Option<Vec<f64>>,
    private_alleles_total: Option<Vec<u32>>,
    private_alleles_in_draw: Option<Vec<f64>>,
    variable_vars_total: Option<Vec<u32>>,
    variable_vars_in_draw: Option<Vec<f64>>,
    /// The bins of the first population, then those of the second, in one
    /// array: a `Vec<Vec<f64>>` is not one of the types wasm-bindgen carries.
    folded_sfs: Option<Vec<f64>>,
    /// How many bins each population has, which cuts `folded_sfs` into the
    /// spectrum of each of them.
    num_sfs_bins: u32,
    fis: Option<Vec<f64>>,
    counts: PassCounts,
}

#[wasm_bindgen]
impl PopDiversityOfAPass {
    /// The names of the populations, in the order the user named them, which
    /// is the order of every array below.
    pub fn pop_names(&mut self) -> Option<Vec<String>> {
        self.pop_names.take()
    }

    /// The variants each population called something at and had
    /// `min_num_individuals` called genotypes in.
    pub fn num_vars_with_data(&mut self) -> Option<Vec<u32>> {
        self.num_vars_with_data.take()
    }

    /// Of those, the ones whose called alleles reached the draw.
    pub fn num_vars_in_draw(&mut self) -> Option<Vec<u32>> {
        self.num_vars_in_draw.take()
    }

    /// The variants that counted for every population, which is the divisor
    /// of the mean private alleles.
    #[must_use]
    pub fn num_vars_every_pop(&self) -> u32 {
        self.num_vars_every_pop
    }

    /// Of those, the ones every population reached the draw at.
    #[must_use]
    pub fn num_vars_every_pop_in_draw(&self) -> u32 {
        self.num_vars_every_pop_in_draw
    }

    /// The alleles each population called, summed over its variants, and
    /// nothing when nobody asked for that statistic.
    pub fn num_alleles_total(&mut self) -> Option<Vec<u32>> {
        self.num_alleles_total.take()
    }

    /// What a draw of `num_called_alleles` is expected to show of them, and
    /// NaN where there is no such value.
    pub fn num_alleles_in_draw(&mut self) -> Option<Vec<f64>> {
        self.num_alleles_in_draw.take()
    }

    /// The alleles each population called that no other population of the
    /// call called at the same variant.
    pub fn private_alleles_total(&mut self) -> Option<Vec<u32>> {
        self.private_alleles_total.take()
    }

    /// The same in a draw of `num_called_alleles`.
    pub fn private_alleles_in_draw(&mut self) -> Option<Vec<f64>> {
        self.private_alleles_in_draw.take()
    }

    /// The variants each population called more than one allele at.
    pub fn variable_vars_total(&mut self) -> Option<Vec<u32>> {
        self.variable_vars_total.take()
    }

    /// The chance that a draw of `num_called_alleles` is not all of one
    /// allele, averaged over the variants in the draw.
    pub fn variable_vars_in_draw(&mut self) -> Option<Vec<f64>> {
        self.variable_vars_in_draw.take()
    }

    /// How many variants of each population a draw of `num_called_alleles` is
    /// expected to show each count of the rarer allele at, and nothing when
    /// nobody asked for the spectrum: the bins of the first population, then
    /// those of the second, which [`PopDiversityOfAPass::num_sfs_bins`] cuts
    /// into one spectrum per population.
    pub fn folded_sfs(&mut self) -> Option<Vec<f64>> {
        self.folded_sfs.take()
    }

    /// How many bins the spectrum of each population has, which is
    /// `num_called_alleles / 2 + 1`, and 0 when nobody asked for the spectrum.
    #[must_use]
    pub fn num_sfs_bins(&self) -> u32 {
        self.num_sfs_bins
    }

    /// The inbreeding coefficient of each population, NaN where it has none.
    pub fn fis(&mut self) -> Option<Vec<f64>> {
        self.fis.take()
    }

    /// How many variants the pass gave, and what each filter of it was given
    /// and kept, the outermost filter first.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.counts.clone()
    }
}
