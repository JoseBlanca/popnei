//! The r² of the variants of a pass on its way to Python: one pass that
//! gives the matrix of every pair with the chromosome and the position of
//! each variant beside it, and one that gives, for each population, how the
//! r² of a pair falls off with the distance between its two variants, in
//! bins of distance.
//!
//! The calculations are the core's, [`popnei::ld::calc_r2_matrix`] and
//! [`popnei::ld::calc_ld_and_dist`], and what
//! this module does is what `calc_pairwise_kosman_dists` of `dists.rs` does
//! for the distances: it builds the chain of readers of the pass from the
//! steps of the `Variants`, lends it to the core with the interpreter
//! released, and reads the counts of the filters from that chain when the
//! call is over, since no block of the pass reaches this crate. How many
//! variants the calculation took comes from the result, and the two together
//! are the `pass_stats` of the `R2Matrix` and of the `LdAndDistPerPop` the
//! package builds.
//!
//! A pair that has no r², which "What it gives" of `docs/specs/ld.md`
//! defines, is NaN in the core already, so nothing of that missing value is
//! decided here. The one that is decided here is a bin with no pair, whose
//! mean and standard deviation the core gives as `None` and which a user
//! reads as NaN.

use std::path::Path;

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use popnei::block::BlockReader;
use popnei::ld::{LdAndDistOptions, LdBins, R2Matrix};

use crate::errors::{PyPopneiError, raise_a_ctrl_c_before_numpy_is_called};
use crate::source::{
    ChromColumn, OpenSource, PassCounts, chrom_column, count_of, distance_of, read_only, source_of,
    threshold_of,
};
use crate::stats::{of_a_result, the_pops};
use crate::steps::{Step, Steps, chain_of};

/// What each filter of a pass was given and kept, under the kind of the
/// filter and in the order of the chain, which is the half of the counts of
/// a pass that the chain of readers holds.
type FilteringCounts = Vec<(&'static str, u64, u64)>;

/// What one pass gives: the matrix of the r² of every pair, and what each
/// filter of it counted.
type TheMatrixOfThePass = (R2Matrix, FilteringCounts);

/// The same on its way to Python: the r² as a square numpy array of float64
/// of the variants of the pass, the name of the chromosome of each of them,
/// their positions, and the counts of the pass.
type R2MatrixForPython<'py> = (
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyTuple>,
    Bound<'py, PyArray1<u64>>,
    PassCounts,
);

// The r² of every pair of the variants that the steps of `steps` keep of
// `source`, with `max_num_vars` the variants the calculation takes before it
// refuses. What it gives back is the matrix row after row, the name of the
// chromosome of each variant, their positions and the counts of the pass. A
// `///` comment here would become the `__doc__` of
// `popnei._core.calc_rogers_huff_r2_matrix`, and what a Python user reads
// belongs to the package, which is the API.
#[pyfunction]
#[pyo3(signature = (source, max_num_vars, steps))]
pub(crate) fn calc_rogers_huff_r2_matrix<'py>(
    py: Python<'py>,
    source: &Bound<'_, PyAny>,
    max_num_vars: &Bound<'_, PyAny>,
    steps: &Bound<'_, Steps>,
) -> Result<R2MatrixForPython<'py>, PyPopneiError> {
    let source = source_of(source)?;
    // The cap is taken as the object a user wrote and converted here, so
    // that `2.5` and `True` are refused under the name of the argument: the
    // conversion of the signature would take `True` as the number 1 and say
    // nothing.
    let max_num_vars = count_of("max_num_vars", max_num_vars)?;
    // A cap of 0 variants is refused here and not by the core, which takes
    // it and stops the pass at its first variant: what a user would read
    // then is "the pass gave 1 variants and `max_num_vars` is 0 ... raise
    // `max_num_vars` or filter the variants", the message of a dataset too
    // large for the cap, where what they asked for is a matrix of no
    // variant. `count_of` leaves the floor of every other count of popnei
    // to the core, which has a message of its own for each of them and
    // none for this one, and `calcRogersHuffR2Matrix` of TypeScript refuses
    // a `maxNumVars` of 0 at the call in the same words.
    if max_num_vars == 0 {
        return Err(PyPopneiError::Count {
            name: "max_num_vars",
            smallest: 1,
            value: max_num_vars.to_string(),
        });
    }
    let steps = steps.get().of_a_pass()?;
    // A Ctrl-C that was pending when this was called is raised here, before
    // the file is opened.
    py.check_signals()?;
    let path = source.path().to_path_buf();
    // The whole source is read inside this one call, and the products of the
    // tiles are the work of every core of the machine, so the interpreter is
    // released for all of it: the threads of rayon deadlock on a caller that
    // holds it. A Ctrl-C that arrives meanwhile is raised when the call is
    // over and not between two blocks, as it is in `Blocks::__next__`: the
    // loop over the blocks is the core's, and a pass that is interrupted
    // loses only itself, since it writes no file and the `Variants` is as it
    // was.
    let calculated = py.detach(|| over_the_source(source, &steps, max_num_vars));
    // The file of the source goes into every error of the pass, and
    // `errors.rs` is what leaves it out of the message of the two that are
    // of the cap a user wrote, as it does for every argument that is refused
    // while a file is being read.
    let (matrix, filtering) =
        calculated.map_err(|error| PyPopneiError::of_the_file(error, &path))?;
    // The variants of the pass, which are the rows of the matrix, go to
    // Python as the `u64` every count of popnei is there: a `usize` is 32
    // bits in WebAssembly and 64 natively, and what a user reads does not
    // depend on that.
    let num_vars = u64::try_from(matrix.num_vars()).map_err(|_| {
        PyPopneiError::broken_of_the_file(
            format!(
                "the pass gave {num_vars} variants, which is more than a count holds",
                num_vars = matrix.num_vars()
            ),
            &path,
        )
    })?;
    // The Ctrl-C that arrived while the interpreter was released is raised
    // before numpy is called: the first array of a process imports the C API
    // of numpy, that import fails with the exception that is pending, and
    // the numpy crate panics when it does, which a user cannot catch.
    py.check_signals()?;
    // The matrix is taken out of the result and not read from it: numpy is
    // given the allocation the core filled, 200 MB at the cap of 5000
    // variants, instead of a second one written from it. Everything the
    // result holds is read from what it gave away, since `given_away`
    // consumes it.
    let matrix = matrix.given_away();
    let chroms = ChromColumn::of(&matrix.chroms, &matrix.chrom_table, &path)?;
    let chroms = chrom_column(py, &chroms, &path)?;
    let poss = read_only(matrix.poss.into_pyarray(py))?;
    let r2 = read_only(the_square_of(py, matrix.num_vars, matrix.r2)?)?;
    Ok((r2, chroms, poss, (num_vars, filtering)))
}

/// One pass over `source` through `steps`, and the r² of every pair of the
/// variants it gave.
///
/// The interpreter is already released here, so nothing of Python is
/// touched: what comes back is Rust.
///
/// # Errors
///
/// When the source cannot be opened or read, when the pass gives no variant,
/// and when the core refuses the calculation: a pass of more variants than
/// `max_num_vars`, a cap whose matrix this machine does not count, a matrix
/// this machine has not the memory for, and a source of more individuals or
/// of a higher ploidy than the sums of the r² are worked out in.
fn over_the_source(
    source: &dyn OpenSource,
    steps: &[Step],
    max_num_vars: usize,
) -> popnei::Result<TheMatrixOfThePass> {
    // The source is opened at the size of its own blocks: the six sums of a
    // pair run over the individuals, which no block cuts, so the same matrix
    // comes out whatever the size, and no `Reblock` is put over the chain.
    let reader = source.reader(None)?;
    // The chain of the pass stays here, lent to the core, so that the counts
    // of its filters can be read when the call is over: the loop over the
    // blocks is the core's, and no block of it reaches this crate.
    let mut chain = chain_of(reader, steps)?;
    // A pass that gave no variant is the core's `PassGaveNoVariant`, which
    // the core builds with the counts it reads from the chain it was lent:
    // they say whether the source had none or the steps kept none, and
    // nothing else would carry them out of a pass that could not be
    // finished.
    let matrix = popnei::ld::calc_r2_matrix(&mut chain, max_num_vars)?;
    Ok((matrix, filtering_of(chain.as_ref())))
}

/// The r² of the matrix as a numpy array of its variants x its variants.
///
/// `values` is the `Vec` the core filled and gave away, and numpy takes it
/// over as it takes the distances of `dists.rs`: nothing of the matrix is
/// copied, and nothing asks this machine for memory here, so the matrix is
/// held once and not twice, 200 MB at the 5000 variants of
/// `DEFAULT_MAX_NUM_VARS` instead of 400 MB. A machine that has not the
/// memory of the matrix is refused by the core, which asks for its own with
/// `try_reserve_exact` and fails with `Error::LdNoMemory`.
///
/// # Errors
///
/// [`PyPopneiError::Broken`] when the core gave a matrix whose values are
/// not the square of its variants, which is a defect of popnei: a user
/// reports it instead of looking for what they typed wrong.
fn the_square_of<'py>(
    py: Python<'py>,
    num_vars: usize,
    values: Vec<f64>,
) -> Result<Bound<'py, PyArray2<f64>>, PyPopneiError> {
    let num_values = values.len();
    let square = Array2::from_shape_vec((num_vars, num_vars), values).map_err(|error| {
        PyPopneiError::Broken {
            message: format!(
                "the r² of the pairs of {num_vars} variants holds {num_values} values: {error}"
            ),
            path: None,
        }
    })?;
    Ok(square.into_pyarray(py))
}

/// What each filter of a chain was given and kept, the outermost filter
/// first, which is the order the package turns around for its user.
fn filtering_of(chain: &dyn BlockReader) -> Vec<(&'static str, u64, u64)> {
    chain
        .filtering_stats()
        .into_iter()
        .map(|(kind, stats)| (kind, stats.vars_processed, stats.vars_kept))
        .collect()
}

/// The bins of one population on their way to Python: the smallest and the
/// largest distance of each bin, both included, how many pairs it holds, the
/// mean of their r² and its standard deviation, and how many variants the
/// population kept at its major allele frequency.
///
/// The five arrays hold one value for each bin, in the order of the
/// distances, and the package makes the rows of a pandas frame out of them.
///
/// The three counts are signed 64 bit integers, which
/// [`crate::stats::of_a_result`] says why: a user subtracts the pairs of
/// one bin from the pairs of another, or the smallest distance of a bin
/// from the largest, and of unsigned counts they read 1.8e19 where the
/// answer is negative.
type BinsOfAPop<'py> = (
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
    u64,
);

/// What one pass of the fall-off gives Python: the names of the populations
/// in their order, the bins of each of them in that same order, and the
/// counts of the pass.
type LdAndDistForPython<'py> = (Vec<String>, Vec<BinsOfAPop<'py>>, PassCounts);

// The fall-off of the r² of a pair of variants with the distance between
// them, in bins of distance and for each population, over the variants that
// the steps of `steps` keep of `source`. A `///` comment here would become
// the `__doc__` of `popnei._core.calc_ld_and_dist_per_pop`, and what a
// Python user reads belongs to the package, which is the API.
#[pyfunction]
#[pyo3(signature = (source, steps, pops, min_dist, max_dist, num_bins, max_allowed_maf))]
#[expect(
    clippy::too_many_arguments,
    reason = "the arguments of `calc_ld_and_dist_per_pop` of `docs/specs/ld.md`, the four \
              of the bins taken as the object a Python user wrote so that what is refused \
              names the argument; a struct of them would be built in Python, element by \
              element"
)]
pub(crate) fn calc_ld_and_dist_per_pop<'py>(
    py: Python<'py>,
    source: &Bound<'py, PyAny>,
    steps: &Bound<'py, Steps>,
    pops: Option<Vec<(String, Vec<String>)>>,
    min_dist: &Bound<'py, PyAny>,
    max_dist: &Bound<'py, PyAny>,
    num_bins: &Bound<'py, PyAny>,
    max_allowed_maf: &Bound<'py, PyAny>,
) -> Result<LdAndDistForPython<'py>, PyPopneiError> {
    let source = source_of(source)?;
    // The four numbers of the bins are taken as the objects a user wrote and
    // converted here, and not by the signature: the conversion of pyo3
    // raises an `OverflowError` for a negative distance and takes `True` as
    // the number 1, and neither of those names the argument. A `min_dist`
    // below 0 cannot reach the core, whose `min_dist` is a `u64`, so
    // `distance_of` is what refuses it, under the name a user wrote and with
    // the number they wrote there.
    let options = LdAndDistOptions {
        min_dist: distance_of("min_dist", 0, min_dist)?,
        max_dist: distance_of("max_dist", 0, max_dist)?,
        num_bins: count_of("num_bins", num_bins)?,
        max_allowed_maf: threshold_of("max_allowed_maf", max_allowed_maf)?,
    };
    let steps = steps.get().of_a_pass()?;
    // A Ctrl-C that was pending when this was called is raised here, before
    // the file is opened.
    py.check_signals()?;
    let path = source.path().to_path_buf();
    // The whole source is read inside this one call, and the products of the
    // tiles of the window are the work of every core of the machine, so the
    // interpreter is released for all of it: the threads of rayon deadlock
    // on a caller that holds it. A Ctrl-C that arrives meanwhile is raised
    // when the call is over and not between two blocks, as it is in
    // `Blocks::__next__`: the loop over the blocks is the core's, and a pass
    // that is interrupted loses only itself, since it writes no file and the
    // `Variants` is as it was.
    let (of_the_pass, pop_names, filtering) = py
        .detach(|| -> popnei::Result<_> {
            // The source is opened at the size of its own blocks: the bins
            // are added up in the order of the variants of the pass, which
            // no block cuts, so the same numbers come out whatever the size,
            // and no `Reblock` is put over the chain.
            let reader = source.reader(None)?;
            // The chain of the pass stays here, lent to the core, so that
            // the counts of its filters can be read when the call is over:
            // the loop over the blocks is the core's, and no block of it
            // reaches this crate.
            let mut chain = chain_of(reader, &steps)?;
            // The names a user wrote are looked up among the individuals the
            // pass gives, which are those of the source after a filter of
            // individuals when the `Variants` carries one, and only the pass
            // knows them. It is the lookup the statistics per population
            // make, so a name that is not an individual of the pass is
            // refused in the same words by both.
            let pops = the_pops(pops.as_deref(), chain.individuals())?;
            let names = (0..pops.len())
                .map(|pop| pops.name(pop).to_owned())
                .collect();
            let of_each_pop: Vec<&[usize]> =
                (0..pops.len()).map(|pop| pops.individuals(pop)).collect();
            let of_the_pass = popnei::ld::calc_ld_and_dist(&mut *chain, &of_each_pop, &options)?;
            let filtering = filtering_of(chain.as_ref());
            Ok((of_the_pass, names, filtering))
        })
        .map_err(|error| PyPopneiError::of_the_file(error, &path))?;
    // The variants the pass gave are the core's count and are not worked out
    // again here: the calculation was given them and counted them with the
    // arithmetic that says what happens on overflow, and a count of this
    // crate beside it would be a second answer to one question.
    let num_vars = of_the_pass.num_vars();
    // The Ctrl-C that arrived while the interpreter was released is raised
    // here, before numpy is called.
    raise_a_ctrl_c_before_numpy_is_called(py)?;
    let mut per_pop = Vec::with_capacity(of_the_pass.num_pops());
    for pop in 0..of_the_pass.num_pops() {
        let bins = of_the_pass.bins_of_pop(pop).ok_or_else(|| {
            PyPopneiError::broken_of_the_file(
                format!(
                    "the pass counted {num_pops} populations and has no bins for the \
                     population {pop}",
                    num_pops = of_the_pass.num_pops()
                ),
                &path,
            )
        })?;
        per_pop.push(the_bins_for_python(py, bins, &path)?);
    }
    Ok((pop_names, per_pop, (num_vars, filtering)))
}

/// The bins of one population as the five arrays and the count the package
/// builds its frame from.
///
/// A bin with no pair has a mean and a standard deviation of NaN, which is
/// where the `None` of the core becomes the missing value that numpy and
/// pandas hold.
///
/// # Errors
///
/// [`PyPopneiError::Broken`] when the core gives no bounds and no count of
/// pairs for a bin it says it has, which is a defect of popnei: a user
/// reports it instead of looking for what they typed wrong; and when a
/// distance or a count of pairs is above what a signed 64 bit integer
/// holds.
fn the_bins_for_python<'py>(
    py: Python<'py>,
    bins: &LdBins,
    path: &Path,
) -> Result<BinsOfAPop<'py>, PyPopneiError> {
    let num_bins = bins.num_bins();
    let mut smallest_dist = Vec::with_capacity(num_bins);
    let mut largest_dist = Vec::with_capacity(num_bins);
    let mut num_pairs = Vec::with_capacity(num_bins);
    let mut mean_r2 = Vec::with_capacity(num_bins);
    let mut sd_r2 = Vec::with_capacity(num_bins);
    for bin in 0..num_bins {
        let (smallest, largest) = bins
            .bounds(bin)
            .ok_or_else(|| not_a_bin("the distances", bin, num_bins, path))?;
        let pairs = bins
            .num_pairs(bin)
            .ok_or_else(|| not_a_bin("the pairs", bin, num_bins, path))?;
        smallest_dist.push(of_a_result(smallest)?);
        largest_dist.push(of_a_result(largest)?);
        num_pairs.push(of_a_result(pairs)?);
        mean_r2.push(bins.mean_r2(bin).unwrap_or(f64::NAN));
        sd_r2.push(bins.sd_r2(bin).unwrap_or(f64::NAN));
    }
    Ok((
        smallest_dist.into_pyarray(py),
        largest_dist.into_pyarray(py),
        num_pairs.into_pyarray(py),
        mean_r2.into_pyarray(py),
        sd_r2.into_pyarray(py),
        bins.num_vars(),
    ))
}

/// What a user is told when the core says it has `num_bins` bins and gives
/// nothing for one of them, `what` naming the values that were asked for.
fn not_a_bin(what: &str, bin: usize, num_bins: usize, path: &Path) -> PyPopneiError {
    PyPopneiError::broken_of_the_file(
        format!("{what} of the bin {bin} of {num_bins} are not in the bins of a population"),
        path,
    )
}
