//! How much variety each population holds, on its way between Python and
//! the core.
//!
//! The calculation is the core's, [`popnei::diversity::calc_pop_diversity`],
//! and what this module does is what `pop_dists.rs` does for the distances
//! between populations: it turns the arguments a Python user wrote into what
//! the core takes, builds the chain of readers of the pass from the steps of
//! the `Variants`, lends it to the core with the interpreter released, and
//! reads the counts of the filters from that chain when the call is over,
//! since no block of the pass reaches this crate.
//!
//! The populations are looked up against the individuals of the pass, which
//! only the pass knows, and `Pops` of the `stats` module does it, so that a
//! `pops` argument is read the same way here as in
//! `calc_per_var_distribs`.
//!
//! What goes out is arrays and tuples, one value per population in the order
//! the user named them, and the Python package builds the frozen dataclass
//! and its pandas frames out of them. The core gives the totals and this
//! crate gives no mean and no ratio: each of those is a total over a count
//! of variants that the result carries beside it, and the package divides.
//!
//! `docs/specs/diversity.md` has the design.

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::{PyOverflowError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyBool;

use popnei::block::BlockReader;
use popnei::diversity::{
    DiversityOptions, DiversityStats, calc_pop_diversity as diversity_of_the_pops,
};
use popnei::stats::Pops;

use crate::errors::{PyPopneiError, raise_a_ctrl_c_before_numpy_is_called};
use crate::source::{OpenSource, PassCounts, source_of, written_as};
use crate::stats::{of_a_result, the_min_num_individuals};
use crate::steps::{Step, Steps, chain_of};

/// The name of the argument that says how many called alleles every
/// population is brought down to, as a Python user writes it.
const NUM_CALLED_ALLELES: &str = "num_called_alleles";

/// One count of a result: its total for every population and its
/// standardized value, NaN where there is none, in the order of the
/// populations.
type CountOfEveryPop<'py> = (Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<f64>>);

/// What one pass gives Python: the names of the populations in their order;
/// the variants that counted for each of them and those of them that reached
/// the draw; the same two counts over the populations together, which are
/// the divisors of the private alleles; the alleles called, the private ones
/// and the variable variants, each `None` when nobody asked for it; the
/// folded spectrum as bins x populations; F_IS; and the counts of the pass.
type DiversityOfAPass<'py> = (
    Vec<String>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
    (i64, i64),
    Option<CountOfEveryPop<'py>>,
    Option<CountOfEveryPop<'py>>,
    Option<CountOfEveryPop<'py>>,
    Option<Bound<'py, PyArray2<f64>>>,
    Option<Bound<'py, PyArray1<f64>>>,
    PassCounts,
);

/// The same in Rust, before anything of Python is built: the counts are the
/// core's `u64`, and a statistic nobody asked for is `None`.
struct OfThePass {
    /// The names of the populations, in the order the caller named them.
    pop_names: Vec<String>,
    /// The variants that counted for each population, which every pass has:
    /// the counts of the variants are there whatever statistics were asked
    /// for, and `None` here is a defect of popnei.
    num_vars: Option<Vec<u64>>,
    /// Of those, the ones each population reached `num_called_alleles` at:
    /// the variants in the draw for it, which is 0 everywhere for a pass that
    /// was given no draw.
    num_vars_in_draw: Option<Vec<u64>>,
    /// The variants that counted for every population at once.
    num_vars_every_pop: u64,
    /// Of those, the ones every population reached the draw at, which is the
    /// divisor of the standardized private alleles.
    num_vars_every_pop_in_draw: u64,
    /// The alleles each population called.
    num_alleles: Option<CountOfThePass>,
    /// Of those, the ones no other population called at the same variant.
    private_alleles: Option<CountOfThePass>,
    /// The variants each population called more than one allele at.
    num_variable_vars: Option<CountOfThePass>,
    /// The folded spectrum of every population, laid out as the frame that
    /// holds it.
    folded_sfs: Option<SpectrumOfThePass>,
    /// The inbreeding coefficient of each population, NaN where it has none.
    fis: Option<Vec<f64>>,
    /// The variants the pass gave and what each of its filters counted.
    counts: PassCounts,
}

/// One count of a statistic as the core gives it: the total of every
/// population over the called alleles it has, and the value of every
/// population in a draw of `num_called_alleles`.
///
/// The two are one struct because the core gives both or neither: each
/// accessor is `None` for the same reason, that nobody asked for the
/// statistic.
struct CountOfThePass {
    /// The total of each population, in their order.
    total: Vec<u64>,
    /// The standardized value of each population, NaN where it has none.
    in_draw: Vec<f64>,
}

/// The folded spectrum of every population, laid out as the frame that holds
/// it: one row per count of the rarer allele and one column per population.
///
/// The core keeps the bins of one population together and the frame is indexed
/// by the bin, so the values are read across the populations, which is the one
/// pass that copies them out of the result of the core.
struct SpectrumOfThePass {
    /// The value of every population at the first bin, then at the second, and
    /// so on: the bins x populations array of the result, row after row.
    bins_by_pops: Vec<f64>,
    /// How many bins each population has, which is the rows of that array.
    num_bins: usize,
}

impl SpectrumOfThePass {
    /// The spectrum of every population of `of_each_pop`, in their order.
    ///
    /// `num_bins` is the most bins any one population has, so that a
    /// population of fewer leaves the values short of the shape they are given
    /// and is refused where the array is built, rather than being read as if
    /// it had them.
    fn of_each_pop(of_each_pop: &[&[f64]]) -> SpectrumOfThePass {
        let num_bins = of_each_pop
            .iter()
            .map(|bins| bins.len())
            .max()
            .unwrap_or_default();
        let bins_by_pops = (0..num_bins)
            .flat_map(|bin| {
                of_each_pop
                    .iter()
                    .filter_map(move |bins| bins.get(bin).copied())
            })
            .collect();
        SpectrumOfThePass {
            bins_by_pops,
            num_bins,
        }
    }
}

// How much variety each population of `pops` holds, over one pass of
// `source` through the steps of `steps`. A `///` comment here would become
// the `__doc__` of `popnei._core.calc_pop_diversity`, and what a Python user
// reads belongs to the package, which is the API.
#[pyfunction]
#[pyo3(signature = (source, steps, pops, stats, num_called_alleles, min_num_individuals))]
pub(crate) fn calc_pop_diversity<'py>(
    py: Python<'py>,
    source: &Bound<'py, PyAny>,
    steps: &Bound<'py, Steps>,
    pops: Option<Vec<(String, Vec<String>)>>,
    stats: Vec<String>,
    num_called_alleles: &Bound<'py, PyAny>,
    min_num_individuals: &Bound<'py, PyAny>,
) -> Result<DiversityOfAPass<'py>, PyPopneiError> {
    let source = source_of(source)?;
    let options = DiversityOptions {
        stats: the_stats(&stats)?,
        num_called_alleles: the_num_called_alleles(num_called_alleles)?,
        min_num_individuals: the_min_num_individuals(min_num_individuals)?,
    };
    let steps = steps.get().of_a_pass()?;
    let path = source.path();
    // A Ctrl-C that was pending when this was called is raised here, before
    // the file is opened.
    py.check_signals()?;
    // The whole source is read inside this one call, which is seconds for a
    // file of hundreds of megabytes, so the interpreter is released for all
    // of it: the loop over the blocks is the core's, and it runs the rows of
    // each block on rayon, whose workers would deadlock on a caller that
    // holds the interpreter.
    let of_the_pass = py
        .detach(|| over_the_source(source, &steps, pops.as_deref(), &options))
        // Which of the errors of the pass carries the file it was reading
        // before its message is `popnei::Error::names_the_file` of the core
        // crate, the exhaustive match beside the enum, and not a list here:
        // a second list would have to be kept in step with that one, and the
        // draw above every gene copy of the dataset was in the wrong place
        // in both until 25 September 2026.
        .map_err(|error| PyPopneiError::of_the_file(error, path))?;
    raise_a_ctrl_c_before_numpy_is_called(py)?;
    for_python(py, of_the_pass)
}

/// One pass over `source` through `steps`, and every statistic of `options`
/// for every population of `pops`.
///
/// The interpreter is already released here, so nothing of Python is
/// touched: what comes back is Rust.
///
/// # Errors
///
/// When the source cannot be opened or read, when a population names an
/// individual the pass does not give, and when the core refuses the
/// calculation: the folded spectrum with no draw, a draw of fewer than two
/// alleles, a draw of more alleles than the dataset holds gene copies, a pass
/// that gave no variant, and a variant of more alleles than a count of them
/// holds.
fn over_the_source(
    source: &dyn OpenSource,
    steps: &[Step],
    pops: Option<&[(String, Vec<String>)]>,
    options: &DiversityOptions,
) -> popnei::Result<OfThePass> {
    // The source is opened at the size of its own blocks, and no `Reblock` is
    // put over the chain. The four counts are totals of whole numbers and come
    // out the same whatever the blocks are. F_IS comes out within the
    // tolerance of `docs/specs/diversity.md` and not to the bit: its two means
    // are sums of floats that the pass adds by chunk inside each block, so on
    // the panel it moves 5.2e-14 of the value between the whole file and
    // blocks of 7, which the test of the block sizes in the core measures.
    let reader = source.reader(None)?;
    // The chain of the pass stays here, lent to the core, so that the counts
    // of its filters can be read when the call is over: the loop over the
    // blocks is the core's, and no block of it reaches this crate.
    let mut chain = chain_of(reader, steps)?;
    let pops = match pops {
        Some(named) => Pops::from_names(named, chain.individuals())?,
        None => Pops::all(chain.individuals().len()),
    };
    let pop_names: Vec<String> = (0..pops.len())
        .map(|pop| pops.name(pop).to_owned())
        .collect();
    let of_each_pop: Vec<&[usize]> = (0..pops.len()).map(|pop| pops.individuals(pop)).collect();
    let diversity = diversity_of_the_pops(&mut *chain, &of_each_pop, options)?;
    let num_pops = diversity.num_pops();
    let filtering = chain
        .filtering_stats()
        .into_iter()
        .map(|(kind, stats)| (kind, stats.vars_processed, stats.vars_kept))
        .collect();
    Ok(OfThePass {
        pop_names,
        num_vars: of_every_pop(num_pops, |pop| diversity.num_vars(pop)),
        num_vars_in_draw: of_every_pop(num_pops, |pop| diversity.num_vars_in_draw(pop)),
        num_vars_every_pop: diversity.num_vars_every_pop(),
        num_vars_every_pop_in_draw: diversity.num_vars_every_pop_in_draw(),
        num_alleles: count_of_the_pass(
            num_pops,
            |pop| diversity.num_alleles(pop),
            |pop| diversity.num_alleles_in_draw(pop),
        ),
        private_alleles: count_of_the_pass(
            num_pops,
            |pop| diversity.private_alleles(pop),
            |pop| diversity.private_alleles_in_draw(pop),
        ),
        num_variable_vars: count_of_the_pass(
            num_pops,
            |pop| diversity.num_variable_vars(pop),
            |pop| diversity.variable_vars_ratio_in_draw(pop),
        ),
        folded_sfs: of_every_pop(num_pops, |pop| diversity.folded_sfs(pop))
            .as_deref()
            .map(SpectrumOfThePass::of_each_pop),
        fis: of_every_pop(num_pops, |pop| diversity.fis(pop)),
        counts: (diversity.num_vars_of_the_pass(), filtering),
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

/// One count of a statistic for every population, its total and its
/// standardized value, and `None` when the statistic was not asked for.
///
/// The two accessors of the core answer `None` together, both of them for the
/// statistic nobody asked for, so a count with one of the two is not a case
/// this can give.
fn count_of_the_pass(
    num_pops: usize,
    total: impl Fn(usize) -> Option<u64>,
    in_draw: impl Fn(usize) -> Option<f64>,
) -> Option<CountOfThePass> {
    Some(CountOfThePass {
        total: of_every_pop(num_pops, total)?,
        in_draw: of_every_pop(num_pops, in_draw)?,
    })
}

/// What the pass gave, as the arrays and the tuples the Python package
/// builds its result from.
///
/// # Errors
///
/// When the names of the populations and the counts of the result are not as
/// many, and when a count is above what an array of a result holds, which
/// are both a defect of popnei.
fn for_python<'py>(
    py: Python<'py>,
    of_the_pass: OfThePass,
) -> Result<DiversityOfAPass<'py>, PyPopneiError> {
    let OfThePass {
        pop_names,
        num_vars,
        num_vars_in_draw,
        num_vars_every_pop,
        num_vars_every_pop_in_draw,
        num_alleles,
        private_alleles,
        num_variable_vars,
        folded_sfs,
        fis,
        counts,
    } = of_the_pass;
    let num_pops = pop_names.len();
    // Every pass counts the variants of each of its populations, whatever
    // statistics it was asked for, so a result with no such count is a
    // defect and not a statistic nobody asked for.
    let (Some(num_vars), Some(num_vars_in_draw)) = (num_vars, num_vars_in_draw) else {
        return Err(PyPopneiError::Broken {
            message: format!(
                "the pass counted the variants of no population, and it was over \
                 {num_pops} of them"
            ),
            path: None,
        });
    };
    // The package indexes every frame and every series by these names, and a
    // name and a count that are not of the same population are a wrong
    // number that says nothing about itself.
    if num_vars.len() != num_pops {
        return Err(PyPopneiError::Broken {
            message: format!(
                "the pass gave the names of {num_pops} populations and the counts of \
                 {given}",
                given = num_vars.len()
            ),
            path: None,
        });
    }
    Ok((
        pop_names,
        of_a_result_each(num_vars)?.into_pyarray(py),
        of_a_result_each(num_vars_in_draw)?.into_pyarray(py),
        (
            of_a_result(num_vars_every_pop)?,
            of_a_result(num_vars_every_pop_in_draw)?,
        ),
        count_for_python(py, num_alleles)?,
        count_for_python(py, private_alleles)?,
        count_for_python(py, num_variable_vars)?,
        spectrum_for_python(py, folded_sfs, num_pops)?,
        fis.map(|fis| fis.into_pyarray(py)),
        counts,
    ))
}

/// One count of every population and its standardized value, or `None` when
/// nobody asked for the statistic.
///
/// # Errors
///
/// When a count is above what an array of a result holds.
fn count_for_python<'py>(
    py: Python<'py>,
    count: Option<CountOfThePass>,
) -> Result<Option<CountOfEveryPop<'py>>, PyPopneiError> {
    let Some(CountOfThePass { total, in_draw }) = count else {
        return Ok(None);
    };
    Ok(Some((
        of_a_result_each(total)?.into_pyarray(py),
        in_draw.into_pyarray(py),
    )))
}

/// The folded spectrum of every population as the bins x populations array
/// the Python package indexes by the count of the rarer allele and names by
/// the population, or `None` when nobody asked for it.
///
/// The values were laid out bin by bin where they were copied out of the
/// result of the core, so nothing is copied here: the array takes the vector
/// and numpy takes the array.
///
/// # Errors
///
/// When the populations do not all have the same bins, which is a defect of
/// popnei: the bins of a spectrum are the counts of the rarer allele from 0 to
/// `num_called_alleles` over 2, one draw size for the whole call.
fn spectrum_for_python<'py>(
    py: Python<'py>,
    spectrum: Option<SpectrumOfThePass>,
    num_pops: usize,
) -> Result<Option<Bound<'py, PyArray2<f64>>>, PyPopneiError> {
    let Some(SpectrumOfThePass {
        bins_by_pops,
        num_bins,
    }) = spectrum
    else {
        return Ok(None);
    };
    let spectrum = Array2::from_shape_vec((num_bins, num_pops), bins_by_pops).map_err(|_| {
        PyPopneiError::Broken {
            message: format!(
                "the pass gave the folded spectrum of {num_pops} populations and not \
                 the same {num_bins} bins for each of them"
            ),
            path: None,
        }
    })?;
    Ok(Some(spectrum.into_pyarray(py)))
}

/// Every count of `counts` as the arrays of a result hold them, in their
/// order.
///
/// # Errors
///
/// Those of [`of_a_result`].
fn of_a_result_each(counts: Vec<u64>) -> Result<Vec<i64>, PyPopneiError> {
    counts.into_iter().map(of_a_result).collect()
}

// The names of the statistics that need no draw, which is what a user who
// names no statistic in `stats` asks for: the alleles a population called, the
// private ones among them, the variants that vary in it and F_IS. A `///`
// comment here would become the `__doc__` of
// `popnei._core.diversity_stats_without_a_draw`, and what a Python user reads
// belongs to the package, which is the API.
//
// The Python package builds the default of its `stats` from this list. The
// four are named in the core alone, and not there as well as in the TypeScript
// package, so a statistic that needs no draw is added in one place.
#[pyfunction]
pub fn diversity_stats_without_a_draw() -> Vec<&'static str> {
    DiversityStats::WITHOUT_A_DRAW.names()
}

/// The statistics a user asked for, out of the names the Python package
/// gives them, which are the core's.
///
/// # Errors
///
/// A name that is of no statistic, which a user reaches by calling
/// `popnei._core` themselves: the package takes the members of its
/// `PopDiversityStat` and nothing else. The core refuses it and names the
/// five, so the sentence a user reads is written once for both languages.
fn the_stats(names: &[String]) -> Result<DiversityStats, PyPopneiError> {
    let mut asked_for = DiversityStats::empty();
    for name in names {
        asked_for |= DiversityStats::of_name(name)?;
    }
    Ok(asked_for)
}

/// The `value` a user gave for `num_called_alleles`, and `None` for a pass
/// that takes no draw.
///
/// It is taken as the object it is and converted here, and not by the
/// signature, because the conversion of pyo3 answers before any rule of
/// popnei: it takes `True` as the number 1 with no word, and its
/// `OverflowError` and its `TypeError` name neither the argument nor what
/// the number has to be.
///
/// # Errors
///
/// A negative number and one above what a variant of popnei can have called
/// alleles, which is a `ValueError`; the core refuses a draw below 2 and
/// says why. What is no whole number at all, `3.1`, `"twenty"` or a truth
/// value, which is a `TypeError`: `True` and `False` say nothing about how
/// many alleles to draw.
fn the_num_called_alleles(value: &Bound<'_, PyAny>) -> Result<Option<u32>, PyPopneiError> {
    if value.is_none() {
        return Ok(None);
    }
    // A truth value is a whole number in Python, so it converts to 1 or 0
    // and has to be refused before the conversion is asked for.
    if value.is_instance_of::<PyBool>() {
        return Err(PyTypeError::new_err(no_count_of_alleles(value)).into());
    }
    match value.extract::<u32>() {
        Ok(count) => Ok(Some(count)),
        Err(error) if error.is_instance_of::<PyOverflowError>(value.py()) => {
            Err(PyValueError::new_err(no_count_of_alleles(value)).into())
        }
        Err(_) => Err(PyTypeError::new_err(no_count_of_alleles(value)).into()),
    }
}

/// What a `num_called_alleles` that is no count of called alleles is told,
/// which names the argument and what the user wrote.
fn no_count_of_alleles(value: &Bound<'_, PyAny>) -> String {
    format!(
        "`{NUM_CALLED_ALLELES}` is {given}, and it is how many called alleles every \
         population is brought down to, so that populations of different sizes can be \
         compared: a whole number of 2 or more, or `None` for no draw at all",
        given = written_as(value)
    )
}
