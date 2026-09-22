//! The distances between individuals on their way to Python: one pass over
//! a source, and the distance of every pair as an array.
//!
//! The calculation is the core's, [`popnei::dists::calc_kosman_sums`], and
//! what this module does is what `write_vars` of `vars.rs` does for a file:
//! it builds the chain of readers of the pass from the steps of the
//! `Variants`, lends it to the core with the interpreter released, and reads
//! the counts of the filters from that chain when the call is over, since no
//! block of the pass reaches this crate. How many variants the calculation
//! took comes from the result, and the two together are the `pass_stats` of
//! the `Distances` the package builds.
//!
//! The distance of a pair the calculation gives no distance to is NaN here,
//! which is where the `Option<f64>` of the core becomes the missing value
//! that numpy and pandas hold.

use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

use popnei::block::BlockReader;

use crate::errors::PyPopneiError;
use crate::source::{OpenSource, PassCounts, read_only, source_of};
use crate::steps::{Step, Steps, chain_of};

/// What one pass gives: the distance of every pair, NaN where there is
/// none; the names of the individuals of the source; and the counts of the
/// pass.
type KosmanDists = (Vec<f64>, Vec<String>, PassCounts);

/// The same on its way to Python, with the distances as a numpy array of
/// float64 that holds the allocation the core filled.
type KosmanDistsForPython<'py> = (Bound<'py, PyArray1<f64>>, Vec<String>, PassCounts);

/// What a pass that could not be finished failed with.
///
/// A pass that gave no variant is kept apart from everything else because
/// the message a user reads is built from the counts of the chain, which
/// the core does not have: it says whether the source had no variant or the
/// steps kept none, and what each filter was given and kept.
enum Refusal {
    /// What the core refused, which the caller gives the file of the source.
    Core(popnei::Error),
    /// The pass gave the calculation no variant, with the counts of each
    /// filter of its chain, the outermost first.
    NoVariant(Vec<(&'static str, u64, u64)>),
}

// The Kosman distance of every pair of individuals of `source`, over the
// variants that the steps of `steps` keep, with `min_num_vars` the variants
// a pair needs to get one. What it gives back is the distances in the order
// of the distance vector, NaN for a pair that has none, the names of the
// individuals and the counts of the pass. A `///` comment here would become
// the `__doc__` of `popnei._core.calc_pairwise_kosman_dists`, and what a
// Python user reads belongs to the package, which is the API.
#[pyfunction]
#[pyo3(signature = (source, min_num_vars, steps))]
pub(crate) fn calc_pairwise_kosman_dists<'py>(
    py: Python<'py>,
    source: &Bound<'_, PyAny>,
    min_num_vars: u32,
    steps: &Bound<'_, Steps>,
) -> Result<KosmanDistsForPython<'py>, PyPopneiError> {
    let source = source_of(source)?;
    let steps = steps.get().of_a_pass()?;
    // A Ctrl-C that was pending when this was called is raised here, before
    // the file is opened.
    py.check_signals()?;
    let path = source.path().to_path_buf();
    // The whole source is read inside this one call, minutes for a dataset
    // of a million variants, so the interpreter is released for all of it;
    // it is what lets the core spread the pairs of a block over the threads
    // of rayon, which deadlock on a caller that holds the interpreter.
    let calculated = py.detach(|| over_the_source(source, &steps, min_num_vars));
    let (dists, individuals, counts) = match calculated {
        Ok(calculated) => calculated,
        Err(Refusal::Core(error)) => return Err(PyPopneiError::of_the_file(error, &path)),
        Err(Refusal::NoVariant(filtering)) => {
            return Err(PyPopneiError::NoVariant { path, filtering });
        }
    };
    // The vector of 10000 individuals is 400 MB, and `into_pyarray` hands
    // the allocation the core filled to numpy without copying it.
    let dists = read_only(dists.into_pyarray(py))?;
    Ok((dists, individuals, counts))
}

/// One pass over `source` through `steps`, and the distances of every pair
/// of individuals over the variants it gave.
///
/// The interpreter is already released here, so nothing of Python is
/// touched: what comes back is Rust.
///
/// # Errors
///
/// When the source cannot be opened or read, when the pass gives no
/// variant, and when the core refuses the calculation: sums that no 32 bit
/// number holds, and the two counts of every pair that the machine has not
/// the memory for.
fn over_the_source(
    source: &dyn OpenSource,
    steps: &[Step],
    min_num_vars: u32,
) -> Result<KosmanDists, Refusal> {
    // The source is opened at the size of its own blocks: the calculation
    // adds whole numbers, so the same distances come out whatever the size,
    // and no `Reblock` is put over the chain.
    let reader = source.reader(None).map_err(Refusal::Core)?;
    // The chain of the pass stays here, lent to the core, so that the counts
    // of its filters can be read when the call is over: the loop over the
    // blocks is the core's, and no block of it reaches this crate.
    let mut chain = chain_of(reader, steps).map_err(Refusal::Core)?;
    let sums = match popnei::dists::calc_kosman_sums(&mut chain) {
        Ok(sums) => sums,
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
    let dists = sums
        .dists(min_num_vars)
        .map(|dist| dist.unwrap_or(f64::NAN))
        .collect();
    let individuals = chain.individuals().to_vec();
    Ok((
        dists,
        individuals,
        (sums.num_vars(), filtering_of(chain.as_ref())),
    ))
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
