//! The statistics of the variants and of the individuals, per population,
//! on their way between TypeScript and the core.
//!
//! It builds the two passes of the module. One calculates up to six
//! statistics for every variant and every population and gives back, for
//! each of them, the mean over the variants that had a value and a histogram
//! of them; the other gives the share of the variants at which every
//! individual has no genotype and the share of its called genotypes at which
//! it is heterozygous.
//!
//! Each of the two turns the arguments a TypeScript user wrote into what the
//! core takes, which for the first pass are the statistics they asked for,
//! the bins of the histogram and the thresholds of each statistic and for
//! the second are none; builds the chain of readers of the pass from the
//! steps of the `Variants`, as the writer of the vars file does, and takes
//! from that chain the individuals the pass gives, which are those of the
//! source after a filter of individuals when the variants carry one, under
//! their names for the second pass and as the populations of the first; and
//! reads the counts of the filters from that chain when the pass is over,
//! which are the counts of the pass beside the variants the result counted.
//!
//! What goes out is the arrays of [`PerVarDistribs`] and of
//! [`PerIndividualStats`], which the package puts together into the result
//! objects of `docs/specs/stats.md`: the means of one statistic are one
//! number per population, its histogram counts are the bins of one
//! population after the bins of the one before it, the two rates are one
//! number per individual, and a value the core does not have is NaN.
//!
//! The populations cross flat, as the arguments of the steps do: an array of
//! arrays is not one of the types wasm-bindgen carries, so the names of the
//! individuals of every population come as one array with how many of them
//! each population has beside it.
//!
//! `docs/specs/stats.md` has the design.

use js_sys::Function;
use wasm_bindgen::JsValue;
use wasm_bindgen::prelude::wasm_bindgen;

use popnei::block::BlockReader;
use popnei::stats::{
    ClosedSide, ExpHet, HistBins, Maf, ObsHet, PerVarDistribsConfig, PerVarStat, Pops, SoFar,
};

use crate::errors::JsPopneiError;
use crate::source::{Consumer, OpenSource, PassCounts, TheResultSoFar, the_run_of};
use crate::steps::{Steps, chain_of};

/// The name of the argument that says below which major allele frequency a
/// variant is polymorphic in a population, as a TypeScript user writes it.
const POLY_THRESHOLD: &str = "polyThreshold";

/// The name of the argument that says which kind of bins the histogram
/// holds, as a TypeScript user writes it, and as the core names it in the
/// message that refuses a kind, which is the name a Python user writes.
const BIN_TYPE: &str = "binType";
const BIN_TYPE_IN_THE_CORE: &str = "bin_type";

/// The name of the argument that says what the allele frequencies of the
/// two expected heterozygosities are raised to, as a TypeScript user writes
/// it, and the name the core gives that number, which it calls the ploidy
/// when it is the variants' own.
const PLOIDY: &str = "ploidy";
const THE_EXPONENT_IN_THE_CORE: &str = "exponent";

/// The arguments of the distributions of the statistics of each variant, as
/// they crossed from TypeScript, which `calcPerVarDistribs` and
/// `calcVariantsSummary` give.
///
/// The package has checked that each of them is of the type the core takes,
/// since a number of JavaScript reaches a whole number of the core as 32
/// bits with no error; what is left is what the core says of them, an
/// unknown name of a statistic and a threshold that is no frequency among
/// it.
#[wasm_bindgen]
pub struct ArgumentsOfThePass {
    /// The statistics to calculate, under the names above.
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
    /// How many called genotypes a population needs at a variant to have a
    /// value there.
    pub(crate) min_num_individuals: u32,
    /// The two ends of the histogram, the number of its bins, whether they
    /// are of equal width or of equal ratio, and the edge each bin holds.
    pub(crate) hist_range: (f64, f64),
    pub(crate) num_bins: usize,
    pub(crate) bin_type: String,
    pub(crate) closed: String,
    /// The exponent of the two expected heterozygosities, and the ploidy of
    /// the variants when the user asked for no other.
    pub(crate) ploidy: Option<usize>,
    /// Below this major allele frequency a variant is polymorphic.
    pub(crate) poly_threshold: f64,
}

#[wasm_bindgen]
impl ArgumentsOfThePass {
    /// The arguments of `calcPerVarDistribs` of `docs/specs/stats.md`, as the
    /// package checked them and flat: `stats` holds the name of each
    /// statistic to calculate; the populations are their names, the names of
    /// the individuals of every one of them one after another, and how many
    /// individuals each of them holds, and `pop_names` is nothing when the
    /// user named no population, which is one population of every individual
    /// of the pass; `min_num_individuals` is how many called genotypes a
    /// population needs at a variant to have a value there; `hist_start`,
    /// `hist_end`, `num_bins`, `bin_type` and `closed` are the histogram
    /// every statistic is counted in; `ploidy` is the exponent of the two expected
    /// heterozygosities, and nothing for the ploidy of the variants; and
    /// `poly_threshold` is the major allele frequency below which a variant
    /// is polymorphic in a population.
    #[wasm_bindgen(constructor)]
    #[must_use]
    #[expect(
        clippy::too_many_arguments,
        reason = "the arguments of `calcPerVarDistribs` of `docs/specs/stats.md`, each \
                  one as the package checked it, and the populations flat: an array of \
                  arrays is not one of the types wasm-bindgen carries"
    )]
    pub fn new(
        stats: Vec<String>,
        pop_names: Option<Vec<String>>,
        pop_individuals: Vec<String>,
        num_individuals_per_pop: Vec<u32>,
        min_num_individuals: u32,
        hist_start: f64,
        hist_end: f64,
        num_bins: usize,
        bin_type: String,
        closed: String,
        ploidy: Option<usize>,
        poly_threshold: f64,
    ) -> ArgumentsOfThePass {
        ArgumentsOfThePass {
            stats,
            pop_names,
            pop_individuals,
            num_individuals_per_pop,
            min_num_individuals,
            hist_range: (hist_start, hist_end),
            num_bins,
            bin_type,
            closed,
            ploidy,
            poly_threshold,
        }
    }
}

/// What a consumer of the three that add up over the blocks of a pass, or
/// the summary of the three, is asked to give while the pass runs, as it
/// crossed from TypeScript: the function the package gives for `onSoFar`,
/// nothing when the application gave none, and `soFarEvery`, the seconds
/// between two calls of it.
pub(crate) struct TheResultSoFarAsked {
    pub(crate) told: Option<Function>,
    pub(crate) every_seconds: f64,
}

/// The six per variant statistics of one pass over `source`, through the
/// steps of `steps`.
///
/// The chain of readers of the pass is built here and stays here, lent to
/// the core, so that the counts of its filters are read when the pass is
/// over: the loop over the blocks is the core's.
///
/// The function of `so_far` is given the distributions over the variants
/// read so far, built as the final ones are, with the counts of the pass as
/// they stand after the block.
///
/// # Errors
///
/// When a name of `stats` is of no statistic; when the histogram cannot be
/// made of the range, the number of bins and the kind of bins that were
/// given; when the exponent of the expected heterozygosities is 0 or above
/// 255; when a population names an individual the pass does not give, names
/// one twice or names none, and when `pops` holds no population; when the
/// major allele frequency below which a variant is polymorphic is not a
/// number from 0 to 1; when the source cannot be read; when the pass gives
/// no variant; and the value the function of `so_far` threw.
pub(crate) fn per_var_distribs_of(
    source: &dyn OpenSource,
    steps: &Steps,
    asked: &ArgumentsOfThePass,
    so_far: TheResultSoFarAsked,
) -> Result<PerVarDistribs, JsPopneiError> {
    let before_the_pass = DistribsBeforeThePass::of(source, asked)?;
    the_run_of(source, &Consumer::PerVarDistribs, |run| {
        let reader = source.reader(run, None)?;
        let mut chain = chain_of(reader, steps.steps())?;
        let DistribsOfThePass {
            config,
            pop_names,
            hist_bin_edges,
        } = before_the_pass.over(chain.individuals())?;
        let mut told = TheResultSoFar::from_now(so_far.told, so_far.every_seconds)?;
        let given = popnei::stats::calc_per_var_distribs_with(
            &mut *chain,
            &config,
            &mut |added_up: &dyn SoFar<popnei::stats::PerVarDistribs>| {
                let Some(told) = told.as_mut() else {
                    return Ok(());
                };
                told.after_a_block(run, || {
                    let counts = PassCounts::of_the_filters(
                        added_up.num_vars(),
                        steps.steps(),
                        &added_up.filtering_stats(),
                    );
                    let distribs = added_up.result().map_err(under_its_name)?;
                    Ok(JsValue::from(distribs_of(
                        distribs,
                        pop_names.clone(),
                        hist_bin_edges.clone(),
                        counts,
                    )?))
                })
            },
        )
        .map_err(under_its_name);
        let distribs = TheResultSoFar::what_the_pass_gives(told, given)?;
        let counts = PassCounts::of(distribs.num_vars, steps.steps(), &*chain);
        distribs_of(distribs, pop_names, hist_bin_edges, counts)
    })
}

/// What the distributions of a pass are counted with that is known before
/// the pass starts, out of the arguments a user wrote: the statistics, the
/// bins and their edges, the populations as the user named them, the three
/// statistics that take arguments of their own and the polymorphism
/// threshold.
pub(crate) struct DistribsBeforeThePass {
    stats: Vec<PerVarStat>,
    named: Option<PopsGiven>,
    bins: HistBins,
    hist_bin_edges: Vec<f64>,
    obs_het: ObsHet,
    maf: Maf,
    exp_het: ExpHet,
    poly_threshold: f64,
}

/// What the distributions of a pass are counted with, once the chain of the
/// pass has said which individuals it gives, with the names of the
/// populations and the edges of the bins the result is given under.
pub(crate) struct DistribsOfThePass {
    pub(crate) config: PerVarDistribsConfig,
    pub(crate) pop_names: Vec<String>,
    pub(crate) hist_bin_edges: Vec<f64>,
}

impl DistribsBeforeThePass {
    /// The distributions `asked` asks for over the variants of `source`.
    ///
    /// # Errors
    ///
    /// When a name of the statistics is of no statistic; when the histogram
    /// cannot be made of the range, the number of bins and the kind of bins
    /// that were given; when the exponent of the expected heterozygosities is
    /// 0 or above 255; and when the arrays of the populations do not hold
    /// the individuals of every population, which is a defect of the
    /// package.
    pub(crate) fn of(
        source: &dyn OpenSource,
        asked: &ArgumentsOfThePass,
    ) -> Result<DistribsBeforeThePass, JsPopneiError> {
        let stats = the_stats(&asked.stats)?;
        let (start, end) = asked.hist_range;
        let bins = HistBins::of_kind(&asked.bin_type, start, end, asked.num_bins)
            .map_err(under_its_name)?
            .closed_on(ClosedSide::of_name(&asked.closed)?);
        let named = the_pops_given(asked)?;
        // The ploidy of the variants turns the alleles a population called
        // into called genotypes, for the `min_num_individuals` test, and it
        // is also the exponent of the two expected heterozygosities when the
        // user asks for no other, which the core decides and not this crate.
        let of_the_variants = source.ploidy();
        let obs_het = ObsHet::new(asked.min_num_individuals);
        let maf = Maf::new(of_the_variants, asked.min_num_individuals)?;
        let exp_het = ExpHet::of_the_exponent_asked_for(
            asked.ploidy,
            of_the_variants,
            asked.min_num_individuals,
        )
        .map_err(under_its_name)?;
        // Every statistic of the pass counts its values in these bins, so
        // their edges are the result's and are kept here, where the bins
        // themselves go on to the core.
        let hist_bin_edges = bins.edges().to_vec();
        Ok(DistribsBeforeThePass {
            stats,
            named,
            bins,
            hist_bin_edges,
            obs_het,
            maf,
            exp_het,
            poly_threshold: asked.poly_threshold,
        })
    }

    /// The distributions over `individuals`, those the chain of the pass
    /// gives, which the populations the user named are resolved against.
    ///
    /// # Errors
    ///
    /// When a population names an individual of none of `individuals`,
    /// names one twice or names none, and when the user named no population
    /// at all.
    pub(crate) fn over(self, individuals: &[String]) -> Result<DistribsOfThePass, JsPopneiError> {
        let pops = match self.named {
            Some(named) => Pops::from_names(&named, individuals)?,
            None => Pops::all(individuals.len()),
        };
        let pop_names: Vec<String> = (0..pops.len())
            .map(|pop| pops.name(pop).to_owned())
            .collect();
        Ok(DistribsOfThePass {
            config: PerVarDistribsConfig {
                stats: self.stats,
                pops,
                bins: self.bins,
                obs_het: self.obs_het,
                maf: self.maf,
                exp_het: self.exp_het,
                poly_threshold: self.poly_threshold,
            },
            pop_names,
            hist_bin_edges: self.hist_bin_edges,
        })
    }
}

/// `distribs`, the distributions the core gave over the variants of a pass
/// or of its first blocks, on their way to JavaScript, under the names of
/// the populations and the edges of the bins of the pass and with `counts`,
/// the counts of the pass.
///
/// # Errors
///
/// Those of [`distrib_of`] and [`poly_counts_of`], each a histogram or a
/// count that JavaScript is not given as it is.
pub(crate) fn distribs_of(
    distribs: popnei::stats::PerVarDistribs,
    pop_names: Vec<String>,
    hist_bin_edges: Vec<f64>,
    counts: PassCounts,
) -> Result<PerVarDistribs, JsPopneiError> {
    let popnei::stats::PerVarDistribs {
        obs_het,
        maf,
        exp_het,
        unbiased_exp_het,
        poly_vars_ratio,
        missing_rate,
        num_vars: _,
    } = distribs;
    Ok(PerVarDistribs {
        pop_names,
        hist_bin_edges,
        obs_het: distrib_of(obs_het.as_ref(), PerVarStat::ObsHet)?,
        maf: distrib_of(maf.as_ref(), PerVarStat::Maf)?,
        exp_het: distrib_of(exp_het.as_ref(), PerVarStat::ExpHet)?,
        unbiased_exp_het: distrib_of(unbiased_exp_het.as_ref(), PerVarStat::UnbiasedExpHet)?,
        poly_vars_ratio: poly_counts_of(poly_vars_ratio.as_ref())?,
        missing_rate: distrib_of(missing_rate.as_ref(), PerVarStat::MissingRate)?,
        counts,
    })
}

/// The statistics a user asked for, out of the names the TypeScript package
/// gives them, which are the core's.
///
/// # Errors
///
/// A name that is of no statistic, which a user reaches by writing one in
/// JavaScript: in TypeScript the six are a union of string literals.
fn the_stats(names: &[String]) -> Result<Vec<PerVarStat>, JsPopneiError> {
    let mut asked_for = Vec::with_capacity(names.len());
    for name in names {
        asked_for.push(PerVarStat::of_name(name)?);
    }
    Ok(asked_for)
}

/// The populations a user named, each with the names of its individuals in
/// the order they named them, which is what `Pops::from_names` takes.
pub(crate) type PopsGiven = Vec<(String, Vec<String>)>;

/// The populations a user named, each with the names of its individuals, out
/// of the flat arrays they crossed in, and `None` when they named none.
///
/// The names are not looked up here: they are resolved against the
/// individuals the pass gives, which are those of the source after a filter
/// of individuals when the `Variants` has one, and only the pass knows them.
///
/// # Errors
///
/// When the arrays do not hold the individuals of every population, which is
/// a defect of the package: it is what cuts them.
fn the_pops_given(asked: &ArgumentsOfThePass) -> Result<Option<PopsGiven>, JsPopneiError> {
    let Some(names) = asked.pop_names.as_ref() else {
        return Ok(None);
    };
    Ok(Some(pops_of_the_arrays(
        names,
        &asked.pop_individuals,
        &asked.num_individuals_per_pop,
    )?))
}

/// The populations of `names`, each with the names of its individuals, out
/// of the one array `individuals` holds them all in: the first population
/// takes the first `num_individuals_per_pop[0]` of them, and so on.
///
/// An array of arrays is not one of the types wasm-bindgen carries, so every
/// calculation that takes populations gets them flat and cuts them here.
///
/// The names are not looked up: they are resolved against the individuals
/// the pass gives, which are those of the source after a filter of
/// individuals when the `Variants` has one, and only the pass knows them.
///
/// # Errors
///
/// When the arrays do not hold the individuals of every population, which is
/// a defect of the package: it is what cuts them.
pub(crate) fn pops_of_the_arrays(
    names: &[String],
    individuals: &[String],
    num_individuals_per_pop: &[u32],
) -> Result<PopsGiven, JsPopneiError> {
    let num_pops = names.len();
    let num_counts = num_individuals_per_pop.len();
    if num_pops != num_counts {
        return Err(JsPopneiError::Broken(format!(
            "the pass was given {num_pops} populations and how many individuals \
             {num_counts} of them hold"
        )));
    }
    let mut given = Vec::with_capacity(num_pops);
    let mut first = 0_usize;
    for (name, num_individuals) in names.iter().zip(num_individuals_per_pop) {
        let num_individuals = usize::try_from(*num_individuals).map_err(|_| {
            JsPopneiError::Broken(format!(
                "the population `{name}` holds {num_individuals} individuals, more \
                 than this build counts"
            ))
        })?;
        let past_the_last = first.checked_add(num_individuals).ok_or_else(|| {
            JsPopneiError::Broken(format!(
                "the populations up to `{name}` hold more individuals than this build \
                 counts"
            ))
        })?;
        let of_the_pop = individuals.get(first..past_the_last).ok_or_else(|| {
            JsPopneiError::Broken(format!(
                "the population `{name}` holds {num_individuals} individuals and \
                 the pass was not given the names of every one of them"
            ))
        })?;
        given.push((name.clone(), of_the_pop.to_vec()));
        first = past_the_last;
    }
    Ok(given)
}

/// The mean of each population and its histogram counts, or `None` when
/// nobody asked for the statistic.
///
/// # Errors
///
/// When the histogram of a population does not hold one count for each bin
/// of the distribution, which is a defect of popnei: the package reads the
/// counts of a population by their place, `pop * numBins + bin`, so a
/// population with fewer would give its user the counts of the next one.
/// And when a bin of the histogram counted more variants than a JavaScript
/// array of counts holds, which is more rows than a file of a tab has.
fn distrib_of(
    distrib: Option<&popnei::stats::StatsDistrib>,
    statistic: PerVarStat,
) -> Result<Option<Distrib>, JsPopneiError> {
    let Some(distrib) = distrib else {
        return Ok(None);
    };
    let num_pops = distrib.num_pops();
    let num_bins = distrib.bins().num_bins();
    // A population in which no variant had a value has no mean, and NaN is
    // what the package gives its user for one.
    let mean = (0..num_pops)
        .map(|pop| distrib.mean(pop).unwrap_or(f64::NAN))
        .collect();
    let mut hist_counts = Vec::new();
    for pop in 0..num_pops {
        let of_the_pop = distrib.hist_counts(pop);
        if of_the_pop.len() != num_bins {
            return Err(JsPopneiError::Broken(format!(
                "the histogram of the {name} has {num_bins} bins and holds {given} \
                 counts for one of its {num_pops} populations",
                name = statistic.name(),
                given = of_the_pop.len()
            )));
        }
        for count in of_the_pop {
            hist_counts.push(for_javascript(*count, statistic)?);
        }
    }
    Ok(Some(Distrib { mean, hist_counts }))
}

/// The three counts of the polymorphism ratio and its two ratios, one value
/// per population, or `None` when nobody asked for them.
///
/// # Errors
///
/// When a population counted more variants than a JavaScript array of counts
/// holds.
fn poly_counts_of(
    poly: Option<&popnei::stats::PolyVarsStats>,
) -> Result<Option<PolyCounts>, JsPopneiError> {
    let Some(poly) = poly else {
        return Ok(None);
    };
    let pops = 0..poly.num_pops();
    let mut num_poly = Vec::with_capacity(pops.len());
    let mut num_variable = Vec::with_capacity(pops.len());
    let mut num_vars_with_data = Vec::with_capacity(pops.len());
    for pop in pops.clone() {
        let of_the_ratio = PerVarStat::PolyVarsRatio;
        num_poly.push(for_javascript(poly.num_poly(pop), of_the_ratio)?);
        num_variable.push(for_javascript(poly.num_variable(pop), of_the_ratio)?);
        num_vars_with_data.push(for_javascript(poly.num_vars_with_data(pop), of_the_ratio)?);
    }
    // A ratio whose denominator is 0 has no value, and NaN is what the
    // package gives its user for one.
    let poly_ratio = pops
        .clone()
        .map(|pop| poly.poly_ratio(pop).unwrap_or(f64::NAN))
        .collect();
    let poly_ratio_over_variables = pops
        .map(|pop| poly.poly_ratio_over_variables(pop).unwrap_or(f64::NAN))
        .collect();
    Ok(Some(PolyCounts {
        num_poly,
        num_variable,
        num_vars_with_data,
        poly_ratio,
        poly_ratio_over_variables,
    }))
}

/// `count` as the array of counts of JavaScript holds it.
///
/// The counts of the variants cross as a `Uint32Array`, which is what
/// `docs/specs/stats.md` gives the result, and the core counts in 64 bits.
///
/// # Errors
///
/// When the count is above 4294967295, which is more variants than a file in
/// the memory of a tab holds: that memory addresses 4 GB, and a variant is a
/// row of a file.
fn for_javascript(count: u64, statistic: PerVarStat) -> Result<u32, JsPopneiError> {
    u32::try_from(count).map_err(|_| {
        JsPopneiError::NotInJavaScript(format!(
            "the {name} of one population counted {count} variants, more than a \
             JavaScript array of counts holds",
            name = statistic.name()
        ))
    })
}

/// `error`, with what a user wrote under the name they wrote it in.
///
/// The core names the major allele frequency below which a variant is
/// polymorphic by what it is for, the kind of the bins of a histogram
/// `bin_type`, which is what a Python user writes it as, and what the
/// allele frequencies are raised to the exponent of a statistic of one
/// variant. What a TypeScript user has to look at is the call they wrote,
/// `polyThreshold`, `binType` and `ploidy`, and this crate is what knows
/// those names.
///
/// The other number the core refuses as an exponent is the ploidy of the
/// variants, which a user never writes: a reader refuses a ploidy of 0 or
/// above 255 when the file is opened.
pub(crate) fn under_its_name(error: popnei::Error) -> JsPopneiError {
    if let popnei::Error::PolyThresholdOutOfRange { value } = error {
        return JsPopneiError::Threshold {
            name: POLY_THRESHOLD,
            threshold: value,
        };
    }
    if let popnei::Error::StatPloidyOutOfRange {
        kind: THE_EXPONENT_IN_THE_CORE,
        value,
        largest,
    } = error
    {
        return JsPopneiError::Refused(format!(
            "`{PLOIDY}` is {value}, and it is 1 at least and {largest} at most, the \
             largest ploidy a reader of popnei gives: it is what the allele frequencies \
             of the two expected heterozygosities are raised to, and how many copies \
             the unbiased one draws, which is the ploidy of the variants when it is not \
             given"
        ));
    }
    if matches!(error, popnei::Error::HistBinsOfAnUnknownKind { .. }) {
        // The core writes the name of the argument once, at the start of
        // what it says, and the kind the user wrote comes after it.
        return JsPopneiError::Refused(error.to_string().replacen(
            BIN_TYPE_IN_THE_CORE,
            BIN_TYPE,
            1,
        ));
    }
    JsPopneiError::Core(error)
}

/// The distribution of one statistic on its way to JavaScript.
struct Distrib {
    /// The mean of each population, NaN where no variant of that population
    /// had a value.
    mean: Vec<f64>,
    /// The counts of the bins of each population, the bins of the first
    /// population first.
    hist_counts: Vec<u32>,
}

/// The counts of the polymorphism ratio on their way to JavaScript, one
/// value per population.
struct PolyCounts {
    num_poly: Vec<u32>,
    num_variable: Vec<u32>,
    num_vars_with_data: Vec<u32>,
    /// The two ratios, NaN where the denominator of one is 0.
    poly_ratio: Vec<f64>,
    poly_ratio_over_variables: Vec<f64>,
}

/// What one pass of the per variant statistics gives JavaScript.
///
/// Every array is copied out of the memory of wasm as it is read, and the
/// object itself holds that memory until its `free()` is called, which the
/// package does as soon as it has read every array of it.
#[wasm_bindgen]
pub struct PerVarDistribs {
    pop_names: Vec<String>,
    hist_bin_edges: Vec<f64>,
    obs_het: Option<Distrib>,
    maf: Option<Distrib>,
    exp_het: Option<Distrib>,
    unbiased_exp_het: Option<Distrib>,
    poly_vars_ratio: Option<PolyCounts>,
    missing_rate: Option<Distrib>,
    counts: PassCounts,
}

#[wasm_bindgen]
impl PerVarDistribs {
    /// The name of each population, in the order the user gave them, which
    /// is the order of every array of the result.
    #[must_use]
    pub fn pop_names(&self) -> Vec<String> {
        self.pop_names.clone()
    }

    /// The edges of the bins every histogram counts in, one more than there
    /// are bins.
    #[must_use]
    pub fn hist_bin_edges(&self) -> Vec<f64> {
        self.hist_bin_edges.clone()
    }

    /// The mean observed heterozygosity of each population, and nothing when
    /// nobody asked for that statistic.
    #[must_use]
    pub fn obs_het_mean(&self) -> Option<Vec<f64>> {
        self.obs_het.as_ref().map(|distrib| distrib.mean.clone())
    }

    /// Its histogram counts, the bins of one population after those of the
    /// one before it.
    #[must_use]
    pub fn obs_het_hist_counts(&self) -> Option<Vec<u32>> {
        self.obs_het
            .as_ref()
            .map(|distrib| distrib.hist_counts.clone())
    }

    /// The mean major allele frequency of each population.
    #[must_use]
    pub fn maf_mean(&self) -> Option<Vec<f64>> {
        self.maf.as_ref().map(|distrib| distrib.mean.clone())
    }

    /// Its histogram counts.
    #[must_use]
    pub fn maf_hist_counts(&self) -> Option<Vec<u32>> {
        self.maf.as_ref().map(|distrib| distrib.hist_counts.clone())
    }

    /// The mean plain expected heterozygosity of each population.
    #[must_use]
    pub fn exp_het_mean(&self) -> Option<Vec<f64>> {
        self.exp_het.as_ref().map(|distrib| distrib.mean.clone())
    }

    /// Its histogram counts.
    #[must_use]
    pub fn exp_het_hist_counts(&self) -> Option<Vec<u32>> {
        self.exp_het
            .as_ref()
            .map(|distrib| distrib.hist_counts.clone())
    }

    /// The mean unbiased expected heterozygosity of each population.
    #[must_use]
    pub fn unbiased_exp_het_mean(&self) -> Option<Vec<f64>> {
        self.unbiased_exp_het
            .as_ref()
            .map(|distrib| distrib.mean.clone())
    }

    /// Its histogram counts.
    #[must_use]
    pub fn unbiased_exp_het_hist_counts(&self) -> Option<Vec<u32>> {
        self.unbiased_exp_het
            .as_ref()
            .map(|distrib| distrib.hist_counts.clone())
    }

    /// The mean missing rate of each population, and nothing when nobody
    /// asked for that statistic.
    #[must_use]
    pub fn missing_rate_mean(&self) -> Option<Vec<f64>> {
        self.missing_rate
            .as_ref()
            .map(|distrib| distrib.mean.clone())
    }

    /// Its histogram counts.
    #[must_use]
    pub fn missing_rate_hist_counts(&self) -> Option<Vec<u32>> {
        self.missing_rate
            .as_ref()
            .map(|distrib| distrib.hist_counts.clone())
    }

    /// The polymorphic variants of each population, and nothing when nobody
    /// asked for the polymorphism ratio.
    #[must_use]
    pub fn num_poly(&self) -> Option<Vec<u32>> {
        self.poly_vars_ratio
            .as_ref()
            .map(|poly| poly.num_poly.clone())
    }

    /// The variable variants of each population.
    #[must_use]
    pub fn num_variable(&self) -> Option<Vec<u32>> {
        self.poly_vars_ratio
            .as_ref()
            .map(|poly| poly.num_variable.clone())
    }

    /// The variants that have a major allele frequency in each population.
    #[must_use]
    pub fn num_vars_with_data(&self) -> Option<Vec<u32>> {
        self.poly_vars_ratio
            .as_ref()
            .map(|poly| poly.num_vars_with_data.clone())
    }

    /// The polymorphic variants of each population over the ones with data,
    /// NaN where there are none.
    #[must_use]
    pub fn poly_ratio(&self) -> Option<Vec<f64>> {
        self.poly_vars_ratio
            .as_ref()
            .map(|poly| poly.poly_ratio.clone())
    }

    /// The polymorphic variants of each population over its variable ones,
    /// NaN where there are none.
    #[must_use]
    pub fn poly_ratio_over_variables(&self) -> Option<Vec<f64>> {
        self.poly_vars_ratio
            .as_ref()
            .map(|poly| poly.poly_ratio_over_variables.clone())
    }

    /// How many variants the pass gave, and what each filter of it was given
    /// and kept, the outermost filter first.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.counts.clone()
    }
}

/// The missing rate and the heterozygosity rate of every individual over one
/// pass over `source`, through the steps of `steps`.
///
/// The chain of readers of the pass is built here and stays here, lent to
/// the core, so that the counts of its filters are read when the pass is
/// over: the loop over the blocks is the core's.
///
/// The function of `so_far` is given the rates over the variants read so
/// far, built as the final ones are, with the counts of the pass as they
/// stand after the block.
///
/// # Errors
///
/// When the source cannot be read, a wrong line of a VCF among the causes;
/// when the pass gives no variant; when the chain gave the names of a
/// different number of individuals than the pass gave rates, which is a
/// defect of popnei; and the value the function of `so_far` threw.
pub(crate) fn per_individual_stats_of(
    source: &dyn OpenSource,
    steps: &Steps,
    so_far: TheResultSoFarAsked,
) -> Result<PerIndividualStats, JsPopneiError> {
    the_run_of(source, &Consumer::PerIndividualStats, |run| {
        let reader = source.reader(run, None)?;
        let mut chain = chain_of(reader, steps.steps())?;
        // The rates come out in the order of the rows of the blocks, which is
        // the order of these names: a filter of individuals gives them in the
        // order the user named them.
        let individuals = chain.individuals().to_vec();
        let mut told = TheResultSoFar::from_now(so_far.told, so_far.every_seconds)?;
        let given = popnei::stats::calc_per_individual_stats_with(
            &mut *chain,
            &mut |added_up: &dyn SoFar<popnei::stats::PerIndividualStats>| {
                let Some(told) = told.as_mut() else {
                    return Ok(());
                };
                told.after_a_block(run, || {
                    let counts = PassCounts::of_the_filters(
                        added_up.num_vars(),
                        steps.steps(),
                        &added_up.filtering_stats(),
                    );
                    let stats = added_up.result()?;
                    Ok(JsValue::from(rates_of(
                        &stats,
                        individuals.clone(),
                        counts,
                    )?))
                })
            },
        )
        .map_err(JsPopneiError::from);
        let stats = TheResultSoFar::what_the_pass_gives(told, given)?;
        let counts = PassCounts::of(stats.num_vars(), steps.steps(), &*chain);
        rates_of(&stats, individuals, counts)
    })
}

/// `stats`, the rates the core gave over the variants of a pass or of its
/// first blocks, on their way to JavaScript, under `individuals`, the names
/// the chain of the pass gave, and with `counts`, the counts of the pass.
///
/// # Errors
///
/// When the chain gave the names of a different number of individuals than
/// the pass gave rates, which is a defect of popnei.
pub(crate) fn rates_of(
    stats: &popnei::stats::PerIndividualStats,
    individuals: Vec<String>,
    counts: PassCounts,
) -> Result<PerIndividualStats, JsPopneiError> {
    let num_individuals = stats.num_individuals();
    // The package reads the name of an individual and its two rates at the
    // same place of three arrays, and a name and a rate that are not of the
    // same individual are a wrong number that says nothing about itself.
    if individuals.len() != num_individuals {
        return Err(JsPopneiError::Broken(format!(
            "the pass gave the names of {given} individuals and the rates of \
             {num_individuals}",
            given = individuals.len()
        )));
    }
    let missing_gt_rate = (0..num_individuals)
        .map(|individual| stats.missing_rate(individual))
        .collect();
    // An individual with no called genotype has no heterozygosity rate, and
    // NaN is what the package gives its user for a value the core does not
    // have, as it does for the mean of a population in which no variant had
    // one.
    let obs_het_rate = (0..num_individuals)
        .map(|individual| stats.obs_het_rate(individual).unwrap_or(f64::NAN))
        .collect();
    Ok(PerIndividualStats {
        individuals,
        missing_gt_rate,
        obs_het_rate,
        counts,
    })
}

/// What one pass of the per individual statistics gives JavaScript.
///
/// Every array is copied out of the memory of wasm as it is read, and the
/// object itself holds that memory until its `free()` is called, which the
/// package does as soon as it has read every array of it.
#[wasm_bindgen]
pub struct PerIndividualStats {
    individuals: Vec<String>,
    missing_gt_rate: Vec<f64>,
    obs_het_rate: Vec<f64>,
    counts: PassCounts,
}

#[wasm_bindgen]
impl PerIndividualStats {
    /// The name of each individual the pass gave, in its order, which is the
    /// order of the two arrays of rates.
    #[must_use]
    pub fn individuals(&self) -> Vec<String> {
        self.individuals.clone()
    }

    /// The variants at which each individual has no genotype, a half called
    /// one among them, over the variants of the pass.
    #[must_use]
    pub fn missing_gt_rate(&self) -> Vec<f64> {
        self.missing_gt_rate.clone()
    }

    /// The heterozygous genotypes of each individual over its called ones,
    /// and NaN for an individual that called none of them.
    #[must_use]
    pub fn obs_het_rate(&self) -> Vec<f64> {
        self.obs_het_rate.clone()
    }

    /// How many variants the pass gave, and what each filter of it was given
    /// and kept, the outermost filter first.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.counts.clone()
    }
}

/// How many called genotypes a population needs at a variant to have a value
/// there when the caller says nothing, 20, inherited from pyNei.
#[wasm_bindgen]
#[must_use]
pub fn default_min_num_individuals() -> u32 {
    popnei::stats::DEFAULT_MIN_NUM_INDIVIDUALS
}

/// The two ends of the histogram when the caller says nothing, 0 and 1,
/// which is where these statistics live.
#[wasm_bindgen]
#[must_use]
pub fn default_hist_range() -> Vec<f64> {
    let (start, end) = popnei::stats::DEFAULT_HIST_RANGE;
    vec![start, end]
}

/// How many bins the histogram holds when the caller says nothing, 40.
#[wasm_bindgen]
#[must_use]
pub fn default_num_bins() -> usize {
    popnei::stats::DEFAULT_NUM_BINS
}

/// Which kind of bins the histogram holds when the caller says nothing,
/// those of equal width.
#[wasm_bindgen]
#[must_use]
pub fn default_bin_type() -> String {
    popnei::stats::DEFAULT_BIN_TYPE.to_owned()
}

/// Which edge each bin of the histogram holds when the caller says
/// nothing, the left one, as `numpy.histogram` does.
#[wasm_bindgen]
#[must_use]
pub fn default_closed() -> String {
    popnei::stats::DEFAULT_CLOSED.to_owned()
}

/// The major allele frequency below which a variant is polymorphic in a
/// population when the caller says nothing, 0.95, inherited from pyNei.
#[wasm_bindgen]
#[must_use]
pub fn default_poly_threshold() -> f64 {
    popnei::stats::DEFAULT_POLY_THRESHOLD
}
