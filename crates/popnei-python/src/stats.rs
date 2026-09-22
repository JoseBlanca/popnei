//! The statistics of the variants and of the individuals, per population,
//! on their way between Python and the core.
//!
//! One pass calculates up to five statistics for every variant and every
//! population and gives back, for each of them, the mean over the variants
//! that had a value and a histogram of them. This module builds that pass:
//! it turns the arguments a Python user wrote into what the core takes, the
//! statistics they asked for, the bins of the histogram and the thresholds
//! of each statistic; it builds the chain of readers of the pass from the
//! steps of the `Variants`, as the writer of the vars file does, and the
//! populations against the individuals that chain gives, which are those of
//! the source after a filter of individuals when the variants carry one;
//! and it reads the counts of the filters from that chain when the pass is
//! over, which are the counts of the pass beside the variants the result
//! counted.
//!
//! What goes out is arrays and tuples, and the Python package builds the
//! frozen dataclasses and the pandas frames out of them: one mean and one
//! column of histogram counts for each population, in the order of the
//! populations, with NaN where the core has no value.
//!
//! `docs/specs/stats.md` has the design.

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::{PyOverflowError, PyValueError};
use pyo3::prelude::*;

use popnei::block::BlockReader;
use popnei::stats::{
    ExpHet, HistBins, LINEAR_BINS, LOGARITHMIC_BINS, Maf, ObsHet, PerVarDistribs,
    PerVarDistribsConfig, PerVarStat, PolyVarsStats, Pops, StatsDistrib,
};

use crate::errors::PyPopneiError;
use crate::source::{PassCounts, count_of, source_of, threshold_of};
use crate::steps::{Steps, chain_of};

/// The name a Python user writes each of the five statistics under, which
/// is also the field of the result that holds it.
const OBS_HET: &str = "obs_het";
const MAF: &str = "maf";
const EXP_HET: &str = "exp_het";
const UNBIASED_EXP_HET: &str = "unbiased_exp_het";
const POLY_VARS_RATIO: &str = "poly_vars_ratio";

/// The name of the argument that says how many called genotypes a
/// population needs at a variant to have a value there, as a Python user
/// writes it.
const MIN_NUM_INDIVIDUALS: &str = "min_num_individuals";

/// The distribution of one statistic on its way to Python: the mean of each
/// population, NaN where no variant of that population had a value, and the
/// histogram counts as bins x populations.
type DistribOfAStat<'py> = (Bound<'py, PyArray1<f64>>, Bound<'py, PyArray2<i64>>);

/// The polymorphism ratio on its way to Python: the polymorphic variants of
/// each population, the variable ones, the ones with data, and the two
/// ratios, NaN where the denominator of a ratio is 0.
type PolyCountsOfAPass<'py> = (
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
);

/// What one pass gives Python: the names of the populations in their order,
/// the edges of the bins, the four distributions and the polymorphism
/// ratio, each `None` when nobody asked for it, and the counts of the pass.
type DistribsOfAPass<'py> = (
    Vec<String>,
    Bound<'py, PyArray1<f64>>,
    Option<DistribOfAStat<'py>>,
    Option<DistribOfAStat<'py>>,
    Option<DistribOfAStat<'py>>,
    Option<DistribOfAStat<'py>>,
    Option<PolyCountsOfAPass<'py>>,
    PassCounts,
);

// The five per variant statistics of one pass over `source`, through the
// steps of `steps`. A `///` comment here would become the `__doc__` of
// `popnei._core.calc_per_var_distribs`, and what a Python user reads belongs
// to the package, which is the API.
#[pyfunction]
#[pyo3(signature = (
    source, steps, stats, pops, min_num_individuals, hist_range, num_bins, bin_type,
    ploidy, poly_threshold,
))]
#[expect(
    clippy::too_many_arguments,
    reason = "the arguments of `calc_per_var_distribs` of `docs/specs/stats.md`, each \
              one taken as the object a Python user wrote so that what is refused names \
              the argument; a struct of them would be built in Python, element by element"
)]
pub(crate) fn calc_per_var_distribs<'py>(
    py: Python<'py>,
    source: &Bound<'py, PyAny>,
    steps: &Bound<'py, Steps>,
    stats: Vec<String>,
    pops: Option<Vec<(String, Vec<String>)>>,
    min_num_individuals: &Bound<'py, PyAny>,
    hist_range: (f64, f64),
    num_bins: &Bound<'py, PyAny>,
    bin_type: &str,
    ploidy: Option<&Bound<'py, PyAny>>,
    poly_threshold: &Bound<'py, PyAny>,
) -> Result<DistribsOfAPass<'py>, PyPopneiError> {
    let source = source_of(source)?;
    let stats = the_stats(&stats)?;
    let min_num_individuals = the_min_num_individuals(min_num_individuals)?;
    let bins = the_bins(bin_type, hist_range, count_of("num_bins", num_bins)?)?;
    // The ploidy of the variants turns the alleles a population called into
    // called genotypes, for the `min_num_individuals` test, and it is the
    // exponent of the two expected heterozygosities unless the user asks
    // for another one.
    let of_the_variants = source.ploidy();
    let exponent = match ploidy {
        Some(asked_for) => count_of("ploidy", asked_for)?,
        None => of_the_variants,
    };
    let obs_het = ObsHet::new(min_num_individuals);
    let maf = Maf::new(of_the_variants, min_num_individuals)?;
    let exp_het = ExpHet::new(exponent, of_the_variants, min_num_individuals)?;
    let poly_threshold = threshold_of("poly_threshold", poly_threshold)?;
    // Every statistic of the pass counts its values in these bins, so their
    // edges are the result's and are kept here, where the bins themselves
    // go on to the core.
    let edges = bins.edges().to_vec();
    let steps = steps.get().of_a_pass()?;
    let path = source.path();
    // A Ctrl-C that was pending when this was called is raised here, before
    // the file is opened.
    py.check_signals()?;
    // The whole source is read inside this one call, which is seconds for a
    // file of hundreds of megabytes, so the interpreter is released for all
    // of it: the loop over the blocks is the core's, and it runs the rows of
    // each block on rayon, whose workers would deadlock on an interpreter
    // this thread held.
    let (distribs, pop_names, filtering) = py
        .detach(|| -> Result<_, popnei::Error> {
            let reader = source.reader(None)?;
            // The chain of the pass stays here, lent to the core, so that
            // the counts of its filters can be read when the call is over.
            let mut chain = chain_of(reader, &steps)?;
            let pops = the_pops(pops.as_deref(), chain.individuals())?;
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
                poly_threshold,
            };
            let distribs = popnei::stats::calc_per_var_distribs(&mut *chain, &config)?;
            let filtering = chain
                .filtering_stats()
                .into_iter()
                .map(|(kind, stats)| (kind, stats.vars_processed, stats.vars_kept))
                .collect();
            Ok((distribs, pop_names, filtering))
        })
        .map_err(|error| PyPopneiError::of_the_file(error, path))?;
    // A Ctrl-C that arrived while the pass ran is still pending: the
    // interpreter was released and no bytecode ran to raise it. It is raised
    // here, before numpy is called, because the first array of a process
    // imports the C API of numpy, that import fails with the exception that
    // is pending, and the numpy crate panics when it does.
    py.check_signals()?;
    let PerVarDistribs {
        obs_het,
        maf,
        exp_het,
        unbiased_exp_het,
        poly_vars_ratio,
        num_vars,
    } = distribs;
    Ok((
        pop_names,
        edges.into_pyarray(py),
        distrib_of(py, obs_het.as_ref())?,
        distrib_of(py, maf.as_ref())?,
        distrib_of(py, exp_het.as_ref())?,
        distrib_of(py, unbiased_exp_het.as_ref())?,
        poly_vars_ratio
            .as_ref()
            .map(|poly| poly_counts_of(py, poly))
            .transpose()?,
        (num_vars, filtering),
    ))
}

/// The statistics a user asked for, out of the names the Python package
/// gives them.
///
/// # Errors
///
/// A name that is of no statistic, which a user reaches by calling
/// `popnei._core` themselves: the package takes the members of its
/// `PerVarStat` and nothing else.
fn the_stats(names: &[String]) -> Result<Vec<PerVarStat>, PyPopneiError> {
    names
        .iter()
        .map(|name| match name.as_str() {
            OBS_HET => Ok(PerVarStat::ObsHet),
            MAF => Ok(PerVarStat::Maf),
            EXP_HET => Ok(PerVarStat::ExpHet),
            UNBIASED_EXP_HET => Ok(PerVarStat::UnbiasedExpHet),
            POLY_VARS_RATIO => Ok(PerVarStat::PolyVarsRatio),
            _ => Err(PyValueError::new_err(format!(
                "`{name}` is not one of the statistics of a variant, which are \
                 `{OBS_HET}`, `{MAF}`, `{EXP_HET}`, `{UNBIASED_EXP_HET}` and \
                 `{POLY_VARS_RATIO}`"
            ))
            .into()),
        })
        .collect()
}

/// The populations of the pass: the ones the user named, looked up among
/// `individuals`, those the pass gives, or the one population of every
/// individual when they named none.
///
/// # Errors
///
/// A name that is not an individual of the pass, a name twice in one
/// population, a population that names no individual, and `pops` with no
/// population at all.
fn the_pops(
    pops: Option<&[(String, Vec<String>)]>,
    individuals: &[String],
) -> Result<Pops, popnei::Error> {
    match pops {
        Some(named) => Pops::from_names(named, individuals),
        None => Ok(Pops::all(individuals.len())),
    }
}

/// The bins of the histogram, of equal width or of equal ratio.
///
/// # Errors
///
/// A `bin_type` that is neither of the two names, and what the core refuses
/// of a range and a number of bins: a histogram of no bin, a range that does
/// not run from a number up to a larger one, and a range of bins of equal
/// ratio that starts at 0 or below.
fn the_bins(
    bin_type: &str,
    (start, end): (f64, f64),
    num_bins: usize,
) -> Result<HistBins, PyPopneiError> {
    let bins = match bin_type {
        LINEAR_BINS => HistBins::linear(start, end, num_bins),
        LOGARITHMIC_BINS => HistBins::logarithmic(start, end, num_bins),
        _ => {
            return Err(PyValueError::new_err(format!(
                "`bin_type` is `{bin_type}`, and the bins of a histogram are \
                 `{LINEAR_BINS}`, of equal width, or `{LOGARITHMIC_BINS}`, of equal \
                 ratio; pyNei spells the first one `lineal`, the Spanish word"
            ))
            .into());
        }
    };
    Ok(bins?)
}

/// The `value` a user gave for `min_num_individuals`.
///
/// It is taken as the object it is and converted here, and not by the
/// signature, because the conversion of pyo3 answers before any rule of
/// popnei: its `OverflowError` names neither the argument nor what the
/// number has to be.
///
/// # Errors
///
/// A number that counts no genotype, a negative one or one above what this
/// machine counts, which is a `ValueError`. What is no whole number at all,
/// `3.1` or `"twenty"`, keeps the `TypeError` of pyo3, which says what it
/// was given.
fn the_min_num_individuals(value: &Bound<'_, PyAny>) -> Result<u32, PyPopneiError> {
    match value.extract::<u32>() {
        Ok(count) => Ok(count),
        Err(error) if error.is_instance_of::<PyOverflowError>(value.py()) => {
            Err(PyValueError::new_err(format!(
                "`{MIN_NUM_INDIVIDUALS}` is {value}, and it is how many called genotypes \
                 a population needs at a variant to have a value there: a whole number \
                 of 0 or more"
            ))
            .into())
        }
        Err(error) => Err(error.into()),
    }
}

/// The mean of each population and the histogram counts as bins x
/// populations, or `None` when nobody asked for the statistic.
///
/// # Errors
///
/// When the counts are not the bins of the distribution times its
/// populations, and when a count is above what a count of a result holds,
/// which are both a defect of popnei.
fn distrib_of<'py>(
    py: Python<'py>,
    distrib: Option<&StatsDistrib>,
) -> Result<Option<DistribOfAStat<'py>>, PyPopneiError> {
    let Some(distrib) = distrib else {
        return Ok(None);
    };
    let num_pops = distrib.num_pops();
    let num_bins = distrib.bins().num_bins();
    // A population with no variant that had a value has no mean, and pandas
    // reads a missing value as NaN, which is where the `Option` of the core
    // becomes one.
    let means: Vec<f64> = (0..num_pops)
        .map(|pop| distrib.mean(pop).unwrap_or(f64::NAN))
        .collect();
    // The counts are taken as the core holds them, the bins of one
    // population after the bins of the one before it, and the array is
    // turned around at the end: what the package builds the frame of is one
    // row per bin and one column per population.
    let mut counts = Vec::new();
    for pop in 0..num_pops {
        for count in distrib.hist_counts(pop) {
            counts.push(of_a_result(*count)?);
        }
    }
    let counts = Array2::from_shape_vec((num_pops, num_bins), counts).map_err(|error| {
        PyPopneiError::Broken {
            message: format!(
                "the histogram of a statistic of {num_pops} populations does not hold \
                 {num_bins} counts for each of them: {error}"
            ),
            path: None,
        }
    })?;
    Ok(Some((
        means.into_pyarray(py),
        counts.reversed_axes().into_pyarray(py),
    )))
}

/// The three counts of the polymorphism ratio and the two ratios, one value
/// per population.
///
/// # Errors
///
/// When a count is above what a count of a result holds.
fn poly_counts_of<'py>(
    py: Python<'py>,
    poly: &PolyVarsStats,
) -> Result<PolyCountsOfAPass<'py>, PyPopneiError> {
    let pops = 0..poly.num_pops();
    let num_poly = of_a_result_each(pops.clone().map(|pop| poly.num_poly(pop)))?;
    let num_variable = of_a_result_each(pops.clone().map(|pop| poly.num_variable(pop)))?;
    let with_data = of_a_result_each(pops.clone().map(|pop| poly.num_vars_with_data(pop)))?;
    // A ratio whose denominator is 0 has no value, and pandas reads a
    // missing value as NaN.
    let poly_ratio: Vec<f64> = pops
        .clone()
        .map(|pop| poly.poly_ratio(pop).unwrap_or(f64::NAN))
        .collect();
    let over_variables: Vec<f64> = pops
        .map(|pop| poly.poly_ratio_over_variables(pop).unwrap_or(f64::NAN))
        .collect();
    Ok((
        num_poly.into_pyarray(py),
        num_variable.into_pyarray(py),
        with_data.into_pyarray(py),
        poly_ratio.into_pyarray(py),
        over_variables.into_pyarray(py),
    ))
}

/// `count` as the arrays of a result hold it, which are the signed 64 bit
/// integers of pyNei's series and frames.
///
/// A user subtracts one count from another, the polymorphic variants of a
/// population from its variable ones or the count of one bin from the count
/// of another, and of unsigned counts they read 18446744073709551600 where
/// the answer is -16.
///
/// # Errors
///
/// When the count is above 9223372036854775807, which is more variants than
/// a file holds: a variant is a row of a file.
fn of_a_result(count: u64) -> Result<i64, PyPopneiError> {
    i64::try_from(count).map_err(|_| PyPopneiError::Broken {
        message: format!("one count of the pass is {count}, more than a result holds"),
        path: None,
    })
}

/// The same for every count of `counts`, in their order.
///
/// # Errors
///
/// Those of [`of_a_result`].
fn of_a_result_each(counts: impl Iterator<Item = u64>) -> Result<Vec<i64>, PyPopneiError> {
    counts.map(of_a_result).collect()
}
