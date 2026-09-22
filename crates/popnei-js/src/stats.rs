//! The statistics of the variants, per population, on their way between
//! TypeScript and the core.
//!
//! One pass calculates up to five statistics for every variant and every
//! population and gives back, for each of them, the mean over the variants
//! that had a value and a histogram of them. This module builds that pass:
//! it turns the arguments a TypeScript user wrote into what the core takes,
//! the statistics they asked for, the bins of the histogram and the
//! thresholds of each statistic; it builds the chain of readers of the pass
//! from the steps of the `Variants`, as the writer of the vars file does,
//! and the populations against the individuals that chain gives, which are
//! those of the source after a filter of individuals when the variants carry
//! one; and it reads the counts of the filters from that chain when the pass
//! is over, which are the counts of the pass beside the variants the result
//! counted.
//!
//! What goes out is the arrays of [`PerVarDistribs`], which the package puts
//! together into the result object of `docs/specs/stats.md`: the means of
//! one statistic are one number per population, its histogram counts are the
//! bins of one population after the bins of the one before it, and a
//! population the core has no value for is NaN.
//!
//! The populations cross flat, as the arguments of the steps do: an array of
//! arrays is not one of the types wasm-bindgen carries, so the names of the
//! individuals of every population come as one array with how many of them
//! each population has beside it.
//!
//! `docs/specs/stats.md` has the design.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::block::BlockReader;
use popnei::stats::{ExpHet, HistBins, Maf, ObsHet, PerVarDistribsConfig, PerVarStat, Pops};

use crate::errors::JsPopneiError;
use crate::source::{OpenSource, PassCounts};
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

/// The arguments of one pass, as they crossed from TypeScript.
///
/// The package has checked that each of them is of the type the core takes,
/// since a number of JavaScript reaches a whole number of the core as 32
/// bits with no error; what is left is what the core says of them, an
/// unknown name of a statistic and a threshold that is no frequency among
/// it.
pub(crate) struct ArgumentsOfThePass {
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
    /// The two ends of the histogram, the number of its bins and whether
    /// they are of equal width or of equal ratio.
    pub(crate) hist_range: (f64, f64),
    pub(crate) num_bins: usize,
    pub(crate) bin_type: String,
    /// The exponent of the two expected heterozygosities, and the ploidy of
    /// the variants when the user asked for no other.
    pub(crate) ploidy: Option<usize>,
    /// Below this major allele frequency a variant is polymorphic.
    pub(crate) poly_threshold: f64,
}

/// The five per variant statistics of one pass over `source`, through the
/// steps of `steps`.
///
/// The chain of readers of the pass is built here and stays here, lent to
/// the core, so that the counts of its filters are read when the pass is
/// over: the loop over the blocks is the core's.
///
/// # Errors
///
/// When a name of `stats` is of no statistic; when the histogram cannot be
/// made of the range, the number of bins and the kind of bins that were
/// given; when the exponent of the expected heterozygosities is 0 or above
/// 255; when a population names an individual the pass does not give, names
/// one twice or names none, and when `pops` holds no population; when the
/// major allele frequency below which a variant is polymorphic is not a
/// number from 0 to 1; when the source cannot be read; and when the pass
/// gives no variant.
pub(crate) fn per_var_distribs_of(
    source: &dyn OpenSource,
    steps: &Steps,
    asked: &ArgumentsOfThePass,
) -> Result<PerVarDistribs, JsPopneiError> {
    let stats = the_stats(&asked.stats)?;
    let (start, end) = asked.hist_range;
    let bins =
        HistBins::of_kind(&asked.bin_type, start, end, asked.num_bins).map_err(under_its_name)?;
    let named = the_pops_given(asked)?;
    // The ploidy of the variants turns the alleles a population called into
    // called genotypes, for the `min_num_individuals` test, and it is also
    // the exponent of the two expected heterozygosities when the user asks
    // for no other, which the core decides and not this crate.
    let of_the_variants = source.ploidy();
    let obs_het = ObsHet::new(asked.min_num_individuals);
    let maf = Maf::new(of_the_variants, asked.min_num_individuals)?;
    let exp_het =
        ExpHet::of_the_exponent_asked_for(asked.ploidy, of_the_variants, asked.min_num_individuals)
            .map_err(under_its_name)?;
    // Every statistic of the pass counts its values in these bins, so their
    // edges are the result's and are kept here, where the bins themselves go
    // on to the core.
    let hist_bin_edges = bins.edges().to_vec();
    let reader = source.reader(None)?;
    let mut chain = chain_of(reader, steps.steps())?;
    let pops = match named {
        Some(named) => Pops::from_names(&named, chain.individuals())?,
        None => Pops::all(chain.individuals().len()),
    };
    let pop_names = (0..pops.len())
        .map(|pop| pops.name(pop).to_owned())
        .collect();
    let config = PerVarDistribsConfig {
        stats,
        pops,
        bins,
        obs_het,
        maf,
        exp_het,
        poly_threshold: asked.poly_threshold,
    };
    let distribs =
        popnei::stats::calc_per_var_distribs(&mut *chain, &config).map_err(under_its_name)?;
    let counts = PassCounts::of(distribs.num_vars, &chain.filtering_stats());
    let popnei::stats::PerVarDistribs {
        obs_het,
        maf,
        exp_het,
        unbiased_exp_het,
        poly_vars_ratio,
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
        counts,
    })
}

/// The statistics a user asked for, out of the names the TypeScript package
/// gives them, which are the core's.
///
/// # Errors
///
/// A name that is of no statistic, which a user reaches by writing one in
/// JavaScript: in TypeScript the five are a union of string literals.
fn the_stats(names: &[String]) -> Result<Vec<PerVarStat>, JsPopneiError> {
    let mut asked_for = Vec::with_capacity(names.len());
    for name in names {
        asked_for.push(PerVarStat::of_name(name)?);
    }
    Ok(asked_for)
}

/// The populations a user named, each with the names of its individuals in
/// the order they named them, which is what `Pops::from_names` takes.
type PopsGiven = Vec<(String, Vec<String>)>;

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
    let num_pops = names.len();
    let num_counts = asked.num_individuals_per_pop.len();
    if num_pops != num_counts {
        return Err(JsPopneiError::Broken(format!(
            "the pass was given {num_pops} populations and how many individuals \
             {num_counts} of them hold"
        )));
    }
    let mut given = Vec::with_capacity(num_pops);
    let mut first = 0_usize;
    for (name, num_individuals) in names.iter().zip(&asked.num_individuals_per_pop) {
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
        let of_the_pop = asked
            .pop_individuals
            .get(first..past_the_last)
            .ok_or_else(|| {
                JsPopneiError::Broken(format!(
                    "the population `{name}` holds {num_individuals} individuals and \
                     the pass was not given the names of every one of them"
                ))
            })?;
        given.push((name.clone(), of_the_pop.to_vec()));
        first = past_the_last;
    }
    Ok(Some(given))
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
fn under_its_name(error: popnei::Error) -> JsPopneiError {
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

/// The major allele frequency below which a variant is polymorphic in a
/// population when the caller says nothing, 0.95, inherited from pyNei.
#[wasm_bindgen]
#[must_use]
pub fn default_poly_threshold() -> f64 {
    popnei::stats::DEFAULT_POLY_THRESHOLD
}
