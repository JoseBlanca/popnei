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

use std::path::Path;

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
    /// The variants that counted for every population at once.
    num_vars_every_pop: u64,
    /// The alleles each population called.
    num_alleles: Option<Vec<u64>>,
    /// Of those, the ones no other population called at the same variant.
    private_alleles: Option<Vec<u64>>,
    /// The variants each population called more than one allele at.
    num_variable_vars: Option<Vec<u64>>,
    /// The inbreeding coefficient of each population, NaN where it has none.
    fis: Option<Vec<f64>>,
    /// The variants the pass gave and what each of its filters counted.
    counts: PassCounts,
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
        .map_err(|error| with_its_file(error, path))?;
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
/// alleles, a pass that gave no variant, and a variant of more alleles than
/// a count of them holds.
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
        num_vars_every_pop: diversity.num_vars_every_pop(),
        num_alleles: of_every_pop(num_pops, |pop| diversity.num_alleles(pop)),
        private_alleles: of_every_pop(num_pops, |pop| diversity.private_alleles(pop)),
        num_variable_vars: of_every_pop(num_pops, |pop| diversity.num_variable_vars(pop)),
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
        num_vars_every_pop,
        num_alleles,
        private_alleles,
        num_variable_vars,
        fis,
        counts,
    } = of_the_pass;
    let num_pops = pop_names.len();
    // Every pass counts the variants of each of its populations, whatever
    // statistics it was asked for, so a result with no such count is a
    // defect and not a statistic nobody asked for.
    let Some(num_vars) = num_vars else {
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
    // The draw of a common number of called alleles is work package 3 of
    // `docs/plans/diversity.md`: the core counts no variant in a draw yet
    // and has no standardized value and no spectrum, so the columns that
    // hold them are the missing value here, which is what the spec gives
    // them when there is no draw.
    let none_in_a_draw = || vec![f64::NAN; num_pops];
    Ok((
        pop_names,
        of_a_result_each(num_vars)?.into_pyarray(py),
        vec![0_i64; num_pops].into_pyarray(py),
        (of_a_result(num_vars_every_pop)?, 0),
        count_for_python(py, num_alleles, none_in_a_draw())?,
        count_for_python(py, private_alleles, none_in_a_draw())?,
        count_for_python(py, num_variable_vars, none_in_a_draw())?,
        None,
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
    total: Option<Vec<u64>>,
    in_draw: Vec<f64>,
) -> Result<Option<CountOfEveryPop<'py>>, PyPopneiError> {
    let Some(total) = total else {
        return Ok(None);
    };
    Ok(Some((
        of_a_result_each(total)?.into_pyarray(py),
        in_draw.into_pyarray(py),
    )))
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

/// What the pass failed with, with the file it was reading where that file
/// is part of what went wrong.
///
/// Seven of the nine refusals that "The Rust interface" of
/// `docs/specs/diversity.md` lists are of what a user wrote and are wrong
/// whatever file is read, so they name none, which is what "Errors, and no
/// panics" of `.claude/skills/coding/SKILL.md` asks of an argument that is
/// refused. The list below is those seven and nothing else, in the order that
/// item gives them, so that a reader can count the two against each other.
///
/// A name that is of no statistic is among them although it cannot arrive
/// here today: [`the_stats`] reads a user's `stats` before the pass is built,
/// so that refusal travels back through `?` and never through this function. A
/// list of the reachable cases alone would have to be read together with every
/// call site of this module, and the case is one of the seven whichever call
/// site raises it.
///
/// The two refusals of the item that are not here are of the variants the pass
/// read, and so are the errors of the reader: `PassGaveNoVariant`, where which
/// file it was is what tells a user whether the source held none or the steps
/// kept none, and `MoreAllelesThanACountHolds`, which is a variant of the
/// file. `DiversityMoreVarsThanACountHolds`, a block that says it holds more
/// variants than a count of them holds, is of the file for the same reason and
/// is the one case of the module that item does not carry.
fn with_its_file(error: popnei::Error, path: &Path) -> PyPopneiError {
    if matches!(
        error,
        popnei::Error::DiversitySfsWithoutADraw
            | popnei::Error::DiversityDrawTooSmall { .. }
            | popnei::Error::DiversityWithNoStatistic
            | popnei::Error::DiversityStatOfAnUnknownName { .. }
            | popnei::Error::DiversityPopWithNoIndividual { .. }
            | popnei::Error::DiversityIndividualNotInTheDataset { .. }
            | popnei::Error::DiversityIndividualAskedForTwice { .. }
    ) {
        return PyPopneiError::Core(error);
    }
    PyPopneiError::of_the_file(error, path)
}
