//! The distances between populations on their way to Python: one pass over
//! a source, and every measure that was asked for as an array of one value
//! per pair.
//!
//! The calculation is the core's, [`popnei::pop_dists::calc_pop_dist_sums`],
//! and what this module does is what `dists.rs` does for the distances
//! between individuals: it builds the chain of readers of the pass from the
//! steps of the `Variants`, lends it to the core with the interpreter
//! released, and reads the counts of the filters from that chain when the
//! call is over, since no block of the pass reaches this crate. Those counts
//! and the variants the pass gave are the `pass_stats` of the `PopDists` the
//! package builds.
//!
//! The arguments a Python user wrote are turned here into what the core
//! takes: the measures, which the core names, how the variants are cut into
//! the resampling groups, and how many called genotypes a population needs
//! at a variant. The populations themselves are looked up against the
//! individuals of the pass, which only the pass knows.
//!
//! What a measure has no value for, a pair whose populations share no
//! variant among them, is NaN here, which is where the `Option<f64>` of the
//! core becomes the missing value that numpy and pandas hold.

use std::path::Path;

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::{PyOverflowError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyString};

use popnei::block::BlockReader;
use popnei::pop_dists::{
    GroupId, JackknifeGroups, PopDistMeasure, PopDistOptions, PopDistSums, calc_pop_dist_sums,
};
use popnei::stats::Pops;

use crate::errors::{PyPopneiError, raise_a_ctrl_c_before_numpy_is_called};
use crate::source::{OpenSource, PassCounts, read_only, source_of, written_as};
use crate::stats::{of_a_result, the_min_num_individuals};
use crate::steps::{Step, Steps, chain_of};

/// The name of the argument that says how the variants are cut into the
/// resampling groups, as a Python user writes it.
const JACKKNIFE_GROUP: &str = "jackknife_group";

/// What a user writes in that argument for each variant to be a group of
/// its own.
const PER_VARIANT: &str = "variant";

/// One measure on its way to Python: its value for every pair of
/// populations, NaN where it has none, and its jackknife standard error
/// beside it, `None` when the call asked for no resampling groups.
type MeasureOfThePairs<'py> = (Bound<'py, PyArray1<f64>>, Option<Bound<'py, PyArray1<f64>>>);

/// What one pass gives Python: the names of the populations in their order,
/// which is the order of the pairs; one measure for each that was asked
/// for, in the order they were asked for; how many variants counted for
/// each pair; the f_2 of every pair within every resampling group, groups x
/// pairs, `None` where no groups were asked for; the chromosome and the
/// first and the last position of each group; and the counts of the pass.
type PopDistsOfAPass<'py> = (
    Vec<String>,
    Vec<MeasureOfThePairs<'py>>,
    Bound<'py, PyArray1<i64>>,
    Option<Bound<'py, PyArray2<f64>>>,
    Vec<(String, u64, u64)>,
    PassCounts,
);

/// The same in Rust, before anything of Python is built: the standard
/// errors and the f_2 of the groups are `None` where no groups were asked
/// for, the counts are the core's `u64`, and a group whose chromosome the
/// table of the reader has no name for carries `None` in its place.
struct OfThePass {
    /// The names of the populations, in the order the caller named them.
    pop_names: Vec<String>,
    /// The values of each measure that was asked for and its standard
    /// errors, in the order the measures were asked for.
    of_each_measure: Vec<(Vec<f64>, Option<Vec<f64>>)>,
    /// How many variants counted for each pair, in the order of the
    /// distance vector, and `None` where the core has no count for a pair,
    /// which is a defect: the pairs are built from the populations the core
    /// itself counted over.
    num_vars_of_each_pair: Vec<Option<u64>>,
    /// The f_2 of every pair within every group, the pairs of one group
    /// together, with how many groups there are.
    f2_groups: Option<(usize, Vec<f64>)>,
    /// The chromosome and the two positions of each group.
    group_ids: Vec<(Option<String>, u64, u64)>,
    /// The variants the pass gave and what each of its filters counted.
    counts: PassCounts,
}

// Every measure of `measures` for every pair of the populations of `pops`,
// over one pass of `source` through the steps of `steps`. A `///` comment
// here would become the `__doc__` of `popnei._core.calc_pop_dists`, and what
// a Python user reads belongs to the package, which is the API.
#[pyfunction]
#[pyo3(signature = (source, steps, pops, measures, jackknife_group, min_num_individuals))]
pub(crate) fn calc_pop_dists<'py>(
    py: Python<'py>,
    source: &Bound<'py, PyAny>,
    steps: &Bound<'py, Steps>,
    pops: Vec<(String, Vec<String>)>,
    measures: Vec<String>,
    jackknife_group: &Bound<'py, PyAny>,
    min_num_individuals: &Bound<'py, PyAny>,
) -> Result<PopDistsOfAPass<'py>, PyPopneiError> {
    let source = source_of(source)?;
    // `pops` with no population at all is refused here and not by
    // `Pops::from_names`, whose message asks the caller to leave `pops` out
    // for one population of every individual: that is the statistics, where
    // `pops` is optional, and here it is a required argument. What a user
    // has to do is name two populations, which is what the core says of one
    // population as well.
    if pops.is_empty() {
        return Err(popnei::Error::PopDistsOfFewerThanTwoPops { num_pops: 0 }.into());
    }
    let asked_for = the_measures(&measures)?;
    let options = PopDistOptions {
        min_num_individuals: the_min_num_individuals(min_num_individuals)?,
        groups: the_jackknife_group(jackknife_group)?,
    };
    let steps = steps.get().of_a_pass()?;
    let path = source.path();
    // A Ctrl-C that was pending when this was called is raised here, before
    // the file is opened.
    py.check_signals()?;
    // The whole source is read inside this one call, minutes for a dataset
    // of a million variants, so the interpreter is released for all of it:
    // the loop over the blocks is the core's, and it runs the rows of each
    // block on rayon, whose workers would deadlock on a caller that holds
    // the interpreter.
    let of_the_pass = py
        .detach(|| over_the_source(source, &steps, &pops, &asked_for, &options))
        .map_err(|error| with_its_file(error, path))?;
    raise_a_ctrl_c_before_numpy_is_called(py)?;
    for_python(py, of_the_pass, path)
}

/// What the pass failed with, with the file it was reading where that file
/// is part of what went wrong.
///
/// A `pops` that names one population is wrong whatever file is read, so it
/// names none, which is what "Errors, and no panics" of
/// `.claude/skills/coding/SKILL.md` asks of an argument that is refused.
/// Every other error of the pass is of the variants it read: which file
/// they came from is what a user needs in order to see whether it is the
/// file or the steps that left the calculation with nothing, or what a
/// length that gave too few resampling groups has to become.
fn with_its_file(error: popnei::Error, path: &Path) -> PyPopneiError {
    if matches!(error, popnei::Error::PopDistsOfFewerThanTwoPops { .. }) {
        return PyPopneiError::Core(error);
    }
    PyPopneiError::of_the_file(error, path)
}

/// One pass over `source` through `steps`, and every measure of `asked_for`
/// for every pair of the populations of `pops`.
///
/// The interpreter is already released here, so nothing of Python is
/// touched: what comes back is Rust.
///
/// # Errors
///
/// When the source cannot be opened or read, when a population names an
/// individual the pass does not give, and when the core refuses the
/// calculation: fewer than two populations, a pass that gave no variant,
/// fewer resampling groups than a standard error is built from, and the
/// sums that the machine has not the memory for.
fn over_the_source(
    source: &dyn OpenSource,
    steps: &[Step],
    pops: &[(String, Vec<String>)],
    asked_for: &[PopDistMeasure],
    options: &PopDistOptions,
) -> popnei::Result<OfThePass> {
    // The source is opened at the size of its own blocks: every measure is a
    // ratio of sums that are added over the blocks, so the same numbers come
    // out whatever the size, and no `Reblock` is put over the chain.
    let reader = source.reader(None)?;
    // The chain of the pass stays here, lent to the core, so that the counts
    // of its filters can be read when the call is over: the loop over the
    // blocks is the core's, and no block of it reaches this crate.
    let mut chain = chain_of(reader, steps)?;
    let pops = Pops::from_names(pops, chain.individuals())?;
    let pop_names = (0..pops.len())
        .map(|pop| pops.name(pop).to_owned())
        .collect();
    let sums = calc_pop_dist_sums(&mut *chain, &pops, options)?;
    // The pairs in the order of the distance vector, (0, 1), (0, 2), ...,
    // (1, 2), ..., which is the order of every array of the result.
    let pairs: Vec<(usize, usize)> = (0..sums.num_pops())
        .flat_map(|first| {
            (first..sums.num_pops())
                .skip(1)
                .map(move |second| (first, second))
        })
        .collect();
    let groups_were_asked_for = options.groups != JackknifeGroups::None;
    let of_each_measure = asked_for
        .iter()
        .map(|measure| of_the_measure(&sums, &pairs, *measure, groups_were_asked_for))
        .collect();
    // A pair with no count is not a pair that counted no variant, which is
    // a 0 the core gives: it is a pair the core does not have, and the
    // count of a pair of the distance vector is a count of another pair
    // from there on.
    let num_vars_of_each_pair = pairs
        .iter()
        .map(|(first, second)| sums.num_vars_of(*first, *second))
        .collect();
    let f2_groups = groups_were_asked_for.then(|| f2_of_every_group(&sums, &pairs));
    let group_ids = sums
        .groups()
        .iter()
        .map(|group| named_group(chain.chroms(), *group))
        .collect();
    let filtering = chain
        .filtering_stats()
        .into_iter()
        .map(|(kind, stats)| (kind, stats.vars_processed, stats.vars_kept))
        .collect();
    Ok(OfThePass {
        pop_names,
        of_each_measure,
        num_vars_of_each_pair,
        f2_groups,
        group_ids,
        counts: (sums.num_vars(), filtering),
    })
}

/// One measure of every pair, in the order of `pairs`, with its standard
/// errors beside it where resampling groups were asked for.
///
/// A pair the core has no value for, one whose populations counted no
/// variant together among them, is NaN, and so is a standard error the core
/// has none of, which is a pair whose variants all fell in one group.
fn of_the_measure(
    sums: &PopDistSums,
    pairs: &[(usize, usize)],
    measure: PopDistMeasure,
    groups_were_asked_for: bool,
) -> (Vec<f64>, Option<Vec<f64>>) {
    let values = sums
        .measures(measure)
        .map(|value| value.unwrap_or(f64::NAN))
        .collect();
    let standard_errors = groups_were_asked_for.then(|| {
        pairs
            .iter()
            .map(|(first, second)| {
                sums.standard_error(measure, *first, *second)
                    .unwrap_or(f64::NAN)
            })
            .collect()
    });
    (values, standard_errors)
}

/// The f_2 of every pair within every group, the pairs of one group
/// together, with how many groups there are: a table of groups x pairs that
/// f_3 and f_4 are built from later without reading the genotypes again.
fn f2_of_every_group(sums: &PopDistSums, pairs: &[(usize, usize)]) -> (usize, Vec<f64>) {
    let num_groups = sums.groups().len();
    let values = (0..num_groups)
        .flat_map(|group| {
            pairs.iter().map(move |(first, second)| {
                sums.f2_of_group(group, *first, *second).unwrap_or(f64::NAN)
            })
        })
        .collect();
    (num_groups, values)
}

/// One group with the name of its chromosome, and `None` in its place when
/// the table of the reader holds no name for the number the group carries,
/// which is a defect of popnei.
fn named_group(chroms: &popnei::variant::ChromTable, group: GroupId) -> (Option<String>, u64, u64) {
    (
        chroms.name(group.chrom).map(ToOwned::to_owned),
        group.start,
        group.end,
    )
}

/// What the pass gave, as the arrays and the tuples the Python package
/// builds its result from.
///
/// # Errors
///
/// When a group has no name for its chromosome, when a pair of the distance
/// vector has no count of its variants and when a count of the variants of
/// a pair is above what an array of a result holds, which are all three a
/// defect of popnei; and when an array of numpy cannot be made read only.
fn for_python<'py>(
    py: Python<'py>,
    of_the_pass: OfThePass,
    path: &Path,
) -> Result<PopDistsOfAPass<'py>, PyPopneiError> {
    let OfThePass {
        pop_names,
        of_each_measure,
        num_vars_of_each_pair,
        f2_groups,
        group_ids,
        counts,
    } = of_the_pass;
    let num_pairs = num_vars_of_each_pair.len();
    let of_each_measure = of_each_measure
        .into_iter()
        .map(|(values, standard_errors)| {
            Ok((
                read_only(values.into_pyarray(py))?,
                standard_errors
                    .map(|errors| read_only(errors.into_pyarray(py)))
                    .transpose()?,
            ))
        })
        .collect::<Result<Vec<_>, PyPopneiError>>()?;
    let num_vars = num_vars_of_each_pair
        .into_iter()
        .enumerate()
        .map(|(pair, counted)| match counted {
            Some(counted) => of_a_result(counted),
            None => Err(PyPopneiError::broken_of_the_file(
                format!(
                    "the pass counted the variants of no pair at the place {pair} of the \
                     distance vector"
                ),
                path,
            )),
        })
        .collect::<Result<Vec<i64>, PyPopneiError>>()?;
    let f2_groups = f2_groups
        .map(|(num_groups, values)| the_table_of_the_groups(py, num_groups, num_pairs, values))
        .transpose()?;
    let group_ids = group_ids
        .into_iter()
        .map(|(chrom, start, end)| match chrom {
            Some(chrom) => Ok((chrom, start, end)),
            None => Err(PyPopneiError::broken_of_the_file(
                format!(
                    "the resampling group of the positions {start} to {end} has no name \
                     for its chromosome"
                ),
                path,
            )),
        })
        .collect::<Result<Vec<_>, PyPopneiError>>()?;
    Ok((
        pop_names,
        of_each_measure,
        read_only(num_vars.into_pyarray(py))?,
        f2_groups,
        group_ids,
        counts,
    ))
}

/// The f_2 of every pair within every group as the groups x pairs array the
/// result carries.
///
/// # Errors
///
/// When the values are not one for each pair of each group, which is a
/// defect of popnei, and when the array cannot be made read only.
fn the_table_of_the_groups<'py>(
    py: Python<'py>,
    num_groups: usize,
    num_pairs: usize,
    values: Vec<f64>,
) -> Result<Bound<'py, PyArray2<f64>>, PyPopneiError> {
    let table = Array2::from_shape_vec((num_groups, num_pairs), values).map_err(|error| {
        PyPopneiError::Broken {
            message: format!(
                "the f_2 of {num_pairs} pairs within each of {num_groups} resampling \
                 groups does not hold one value for each of them: {error}"
            ),
            path: None,
        }
    })?;
    read_only(table.into_pyarray(py))
}

/// The measures a user asked for, out of the names the Python package gives
/// them, which are the core's.
///
/// # Errors
///
/// A name that is of none of the seven, which a user reaches by calling
/// `popnei._core` themselves: the package takes the members of its
/// `PopDistMeasure` and nothing else.
fn the_measures(names: &[String]) -> Result<Vec<PopDistMeasure>, PyPopneiError> {
    let mut asked_for = Vec::with_capacity(names.len());
    for name in names {
        asked_for.push(PopDistMeasure::of_name(name)?);
    }
    Ok(asked_for)
}

/// How the variants are cut into the resampling groups, out of the `None`,
/// the `"variant"` or the length in base pairs a user wrote.
///
/// It is taken as the object it is and converted here, and not by the
/// signature, because three kinds of value are one argument and because the
/// conversion of pyo3 answers before any rule of popnei: it takes `True` as
/// the length 1 with no word, and its `OverflowError` and its `TypeError`
/// name neither the argument nor what the three values are.
///
/// # Errors
///
/// A length of 0 base pairs, which is no stretch of a chromosome, a
/// negative one and one above what a position of popnei counts, which are a
/// `ValueError`. What is none of the three kinds, a `2.5`, a truth value or
/// a word that is not `"variant"`, which is a `TypeError` but for the word:
/// a string is the kind the argument takes, and one that is not `"variant"`
/// is the wrong value of that kind.
fn the_jackknife_group(value: &Bound<'_, PyAny>) -> Result<JackknifeGroups, PyPopneiError> {
    if value.is_none() {
        return Ok(JackknifeGroups::None);
    }
    if value.is_instance_of::<PyString>() {
        return match value.extract::<String>() {
            Ok(written) if written == PER_VARIANT => Ok(JackknifeGroups::PerVariant),
            Ok(_) | Err(_) => Err(PyValueError::new_err(no_length(value)).into()),
        };
    }
    // A truth value is a whole number in Python, so `True` would be a group
    // of one base pair, which says nothing about the linkage disequilibrium
    // of the populations being compared.
    if value.is_instance_of::<PyBool>() {
        return Err(PyTypeError::new_err(no_length(value)).into());
    }
    match value.extract::<u64>() {
        Ok(0) => Err(PyValueError::new_err(no_length(value)).into()),
        Ok(length) => Ok(JackknifeGroups::OfBasePairs(length)),
        // A negative length and one above 18446744073709551615, which no
        // position of a chromosome reaches.
        Err(error) if error.is_instance_of::<PyOverflowError>(value.py()) => {
            Err(PyValueError::new_err(no_length(value)).into())
        }
        Err(_) => Err(PyTypeError::new_err(no_length(value)).into()),
    }
}

/// What a `jackknife_group` that is none of the three values is told, which
/// names the argument, what the user wrote and the three.
fn no_length(value: &Bound<'_, PyAny>) -> String {
    format!(
        "`{JACKKNIFE_GROUP}` is {given}, and it says how the variants are cut into the \
         resampling groups the standard errors are built from: a length in base pairs \
         of 1 or more, `\"{PER_VARIANT}\"` for a group of each variant, or `None` for no \
         standard error",
        given = written_as(value)
    )
}
