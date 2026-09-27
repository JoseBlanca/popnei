//! The density of the variants along the chromosomes on its way between
//! Python and the core.
//!
//! It turns the `window_size` and the `chrom_lengths` a Python user wrote
//! into the numbers the core takes, builds the chain of readers of the pass
//! from the steps of the `Variants`, lends it to the core, and reads the
//! counts of its filters when the pass is over. What goes out is arrays: the
//! name of each chromosome with how many windows it has, and the start, the
//! end and the count of every window, which the Python package makes into
//! the frame of `VarDensity`. The core gives each start and end, so nothing
//! here works one out.
//!
//! "The density of the variants along the chromosomes" of
//! `docs/specs/stats.md` has the design.

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use popnei::block::BlockReader;

use crate::errors::PyPopneiError;
use crate::source::{PassCounts, distance_of, read_only, source_of};
use crate::steps::{Steps, chain_of};

/// The name of the argument of the width of a window, as a Python user
/// writes it.
const WINDOW_SIZE: &str = "window_size";

/// The name of the argument of the lengths of the chromosomes, as a Python
/// user writes it.
const CHROM_LENGTHS: &str = "chrom_lengths";

/// What one pass of [`calc_var_density`] gives Python: the name of each
/// chromosome in the order of the result and how many windows it has, then
/// the start, the end and the count of every window, and the counts of the
/// pass.
type DensityOfThePass<'py> = (
    Vec<String>,
    Vec<usize>,
    Bound<'py, PyArray1<u64>>,
    Bound<'py, PyArray1<u64>>,
    Bound<'py, PyArray1<u32>>,
    PassCounts,
);

// The number of variants in each window of `window_size` base pairs along
// each chromosome, over one pass of `source` through the steps of `steps`,
// with the lengths of `chrom_lengths` when it is not None, pairs of a name
// and a length. A `///` comment here would become the `__doc__` of
// `popnei._core.calc_var_density`, and what a Python user reads belongs to
// the package, which is the API.
#[pyfunction]
#[pyo3(signature = (source, steps, window_size, chrom_lengths))]
pub(crate) fn calc_var_density<'py>(
    py: Python<'py>,
    source: &Bound<'py, PyAny>,
    steps: &Bound<'py, Steps>,
    window_size: &Bound<'py, PyAny>,
    chrom_lengths: Option<Vec<(String, Bound<'py, PyAny>)>>,
) -> Result<DensityOfThePass<'py>, PyPopneiError> {
    // A width or a length of 0 goes on to the core, which says what is wrong
    // with it; what is refused here is what is no whole number of 0 or more.
    let window_size = distance_of(WINDOW_SIZE, 1, window_size)?;
    let chrom_lengths = chrom_lengths
        .map(|lengths| {
            lengths
                .into_iter()
                .map(|(chrom, length)| Ok((chrom, distance_of(CHROM_LENGTHS, 1, &length)?)))
                .collect::<Result<Vec<(String, u64)>, PyPopneiError>>()
        })
        .transpose()?;
    let source = source_of(source)?;
    let steps = steps.get().of_a_pass()?;
    let path = source.path();
    // A Ctrl-C that was pending when this was called is raised here, before
    // the file is opened.
    py.check_signals()?;
    // The whole source is read inside this one call, so the interpreter is
    // released for all of it, as for every other pass.
    let (density, filtering) = py
        .detach(|| -> Result<_, popnei::Error> {
            let reader = source.reader(None)?;
            // The chain of the pass stays here, lent to the core, so that
            // the counts of its filters can be read when the call is over.
            let mut chain = chain_of(reader, &steps)?;
            let density = popnei::stats::calc_var_density(
                &mut *chain,
                window_size,
                chrom_lengths.as_deref(),
            )?;
            let filtering = chain
                .filtering_stats()
                .into_iter()
                .map(|(kind, stats)| (kind, stats.vars_processed, stats.vars_kept))
                .collect();
            Ok((density, filtering))
        })
        .map_err(|error| PyPopneiError::of_the_file(error, path))?;
    // A Ctrl-C that arrived while the pass ran is raised before numpy is
    // called, for the reason `calc_per_individual_stats` gives.
    py.check_signals()?;
    let names: Vec<String> = density
        .chroms()
        .iter()
        .map(|chrom| chrom.name.clone())
        .collect();
    // A list of whole numbers of Python, which numpy's `repeat` takes as
    // the `intp` it asks for, where an array of `u64` it refuses to cast: a
    // chromosome has at most `MAX_NUM_WINDOWS` windows.
    let windows_per_chrom: Vec<usize> = density
        .chroms()
        .iter()
        .map(|chrom| chrom.counts.len())
        .collect();
    let num_windows = density.num_windows();
    let mut starts = Vec::with_capacity(num_windows);
    let mut ends = Vec::with_capacity(num_windows);
    let mut num_vars = Vec::with_capacity(num_windows);
    for window in density.windows() {
        starts.push(window.start);
        ends.push(window.end);
        num_vars.push(window.num_vars);
    }
    Ok((
        names,
        windows_per_chrom,
        read_only(starts.into_pyarray(py))?,
        read_only(ends.into_pyarray(py))?,
        read_only(num_vars.into_pyarray(py))?,
        (density.num_vars(), filtering),
    ))
}
