//! The r² of every pair of the variants of a pass on its way to Python: one
//! pass over a source, and the matrix of the pairs with the chromosome and
//! the position of each variant beside it.
//!
//! The calculation is the core's, [`popnei::ld::calc_r2_matrix`], and what
//! this module does is what `calc_pairwise_kosman_dists` of `dists.rs` does
//! for the distances: it builds the chain of readers of the pass from the
//! steps of the `Variants`, lends it to the core with the interpreter
//! released, and reads the counts of the filters from that chain when the
//! call is over, since no block of the pass reaches this crate. How many
//! variants the calculation took comes from the result, and the two together
//! are the `pass_stats` of the `R2Matrix` the package builds.
//!
//! A pair that has no r², which "What it gives" of `docs/specs/ld.md`
//! defines, is NaN in the core already, so nothing of the missing value is
//! decided here.

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use popnei::block::BlockReader;
use popnei::ld::R2Matrix;

use crate::errors::PyPopneiError;
use crate::source::{
    ChromColumn, OpenSource, PassCounts, chrom_column, count_of, read_only, source_of,
};
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

/// What a pass that could not be finished failed with.
///
/// A pass that gave no variant is kept apart from everything else because
/// the message a user reads is built from the counts of the chain, which the
/// core does not have: it says whether the source had no variant or the
/// steps kept none, and what each filter was given and kept.
enum Refusal {
    /// What the core refused, which the caller gives the file of the source.
    Core(popnei::Error),
    /// The pass gave the calculation no variant, with the counts of each
    /// filter of its chain, the outermost first.
    NoVariant(Vec<(&'static str, u64, u64)>),
}

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
    let (matrix, filtering) = match calculated {
        Ok(calculated) => calculated,
        // The file of the source goes into every error of the pass, and
        // `errors.rs` is what leaves it out of the message of the two that
        // are of the cap a user wrote, as it does for every argument that is
        // refused while a file is being read.
        Err(Refusal::Core(error)) => return Err(PyPopneiError::of_the_file(error, &path)),
        Err(Refusal::NoVariant(filtering)) => {
            return Err(PyPopneiError::NoVariant { path, filtering });
        }
    };
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
    let chroms = ChromColumn::of(matrix.chroms(), matrix.chrom_table(), &path)?;
    let chroms = chrom_column(py, &chroms, &path)?;
    let poss = read_only(matrix.poss().to_vec().into_pyarray(py))?;
    let r2 = read_only(the_square_of(py, &matrix)?)?;
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
) -> Result<TheMatrixOfThePass, Refusal> {
    // The source is opened at the size of its own blocks: the six sums of a
    // pair run over the individuals, which no block cuts, so the same matrix
    // comes out whatever the size, and no `Reblock` is put over the chain.
    let reader = source.reader(None).map_err(Refusal::Core)?;
    // The chain of the pass stays here, lent to the core, so that the counts
    // of its filters can be read when the call is over: the loop over the
    // blocks is the core's, and no block of it reaches this crate.
    let mut chain = chain_of(reader, steps).map_err(Refusal::Core)?;
    let matrix = match popnei::ld::calc_r2_matrix(&mut chain, max_num_vars) {
        Ok(matrix) => matrix,
        Err(error) => {
            // A pass that gave no variant is told with the counts of its
            // filters, which say whether the source had none or the steps
            // kept none and which nothing else would carry out of a pass
            // that could not be finished.
            if matches!(error, popnei::Error::ReaderGaveNoVariants) {
                return Err(Refusal::NoVariant(filtering_of(chain.as_ref())));
            }
            return Err(Refusal::Core(error));
        }
    };
    Ok((matrix, filtering_of(chain.as_ref())))
}

/// The r² of the matrix as a numpy array of its variants x its variants.
///
/// The values the core filled are copied once into the array: the matrix
/// gives them out as a slice it owns, so there is no allocation to hand
/// over to numpy as the distances of `dists.rs` hand theirs. A core that
/// gave its `Vec` away instead would save the copy.
///
/// # Errors
///
/// [`PyPopneiError::NoMemory`] when this machine did not give the memory of
/// the copy, and [`PyPopneiError::Broken`] when the core gave a matrix
/// whose values are not the square of its variants, which is a defect of
/// popnei: a user reports it instead of looking for what they typed wrong.
fn the_square_of<'py>(
    py: Python<'py>,
    matrix: &R2Matrix,
) -> Result<Bound<'py, PyArray2<f64>>, PyPopneiError> {
    let num_vars = matrix.num_vars();
    let num_values = matrix.r2().len();
    let square = Array2::from_shape_vec((num_vars, num_vars), the_r2_copied(matrix.r2())?)
        .map_err(|error| PyPopneiError::Broken {
            message: format!(
                "the r² of the pairs of {num_vars} variants holds {num_values} values: {error}"
            ),
            path: None,
        })?;
    Ok(square.into_pyarray(py))
}

/// The values of the matrix in a `Vec` of their own, which numpy is given.
///
/// The memory is asked for with `try_reserve_exact` and the values are
/// written into what it gave. `to_vec` asks for it through the allocator
/// that ends the process when the memory is not there, `handle_alloc_error`
/// aborting with no traceback and nothing a user could catch, where
/// `docs/specs/ld.md` asks for an error and not a process that ends. The
/// core asks for the memory of its own matrices the same way, and
/// `crates/popnei-js/src/ld.rs` guards this same copy, where a failed
/// allocation is a trap that leaves the module unusable.
///
/// The matrix is held twice while the copy is made, 400 MB at the 5000
/// variants of `DEFAULT_MAX_NUM_VARS` and 14.4 GB at 30000, so what reaches
/// the error is a pass whose matrix fits once and not twice: 30000 variants
/// on a machine with 12 GB free. No test reaches it, since one that did
/// would have to take the memory of the machine it runs on. What a user
/// gets there is a `ValueError` that says how many values could not be
/// held and that the matrix is held twice while it crosses, and their
/// interpreter goes on.
///
/// # Errors
///
/// [`PyPopneiError::NoMemory`] when this machine did not give the values a
/// second time.
fn the_r2_copied(values: &[f64]) -> Result<Vec<f64>, PyPopneiError> {
    let mut copied: Vec<f64> = Vec::new();
    copied
        .try_reserve_exact(values.len())
        .map_err(|_| PyPopneiError::NoMemory {
            message: format!(
                "this machine did not give the memory of the matrix of r², \
                 {num_values} values of 8 bytes: the core holds the matrix and popnei \
                 copies it into the array numpy is given, so it is held twice while \
                 the copy is made. Ask for the matrix of fewer variants, with a filter \
                 on the variants or a lower `max_num_vars`",
                num_values = values.len()
            ),
        })?;
    // The capacity above is the length of the slice, so nothing here asks
    // the machine for memory again.
    copied.extend_from_slice(values);
    Ok(copied)
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
