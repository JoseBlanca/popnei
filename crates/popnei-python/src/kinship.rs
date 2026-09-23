//! The kinship of every pair of the individuals of a source on its way to
//! Python: one pass over the variants, and the matrix of the pairs.
//!
//! The calculation is the core's, [`popnei::kinship::calc_kinship`], and
//! what this module does is what `calc_pairwise_kosman_dists` of `dists.rs`
//! does for the distances: it builds the chain of readers of the pass from
//! the steps of the `Variants`, lends it to the core with the interpreter
//! released, and reads the counts of the filters from that chain when the
//! call is over, since no block of the pass reaches this crate.
//!
//! Two things are counted here that the core does not give back. One is how
//! many variants the pass gave, used or not, which is the `num_vars` of the
//! counts of a pass and which [`CountedVars`] counts as the blocks go by:
//! the `num_vars` of the core's result is how many variants had variance
//! and were used, which is the other number and which the `Kinship` of the
//! package carries under that name. The other is the names of the
//! individuals of the matrix, which the package puts on both sides of its
//! frame: they are those the pass gives, in its order, or the ones a user
//! named, in theirs.
//!
//! The names a user names are turned into their places among the
//! individuals of the pass here, with `popnei::filters::resolve_individuals`,
//! which is what the filter of individuals is given as well: the core takes
//! the places and knows nothing of the names.

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

use popnei::block::{Block, BlockReader};
use popnei::filters::{FilteringStats, resolve_individuals};
use popnei::kinship::Kinship;
use popnei::variant::{ChromTable, Needs};

use crate::errors::PyPopneiError;
use crate::source::{OpenSource, PassCounts, source_of};
use crate::steps::{Step, Steps, chain_of};

/// What one pass gives: the kinship the core calculated, the names of the
/// individuals it is of in the order of its rows, and the counts of the
/// pass.
type TheKinshipOfThePass = (Kinship, Vec<String>, PassCounts);

/// The same on its way to Python: the matrix as a numpy array of
/// individuals x individuals that holds the allocation the core filled, the
/// names of those individuals, how many variants were used and the counts
/// of the pass.
type KinshipForPython<'py> = (Bound<'py, PyArray2<f64>>, Vec<String>, u64, PassCounts);

// The kinship of the individuals of `source` over the variants that the
// steps of `steps` keep, of those `individuals` names or of all of them when
// it is `None`, reading a variant of more than two alleles with every allele
// that is not the major one counting the same when `transform_to_biallelic`
// is true. What it gives back is the matrix row after row, the names of the
// individuals it is of, how many variants had variance and were used, and
// the counts of the pass. A `///` comment here would become the `__doc__` of
// `popnei._core.calc_kinship`, and what a Python user reads belongs to the
// package, which is the API.
#[pyfunction]
#[pyo3(signature = (source, individuals, transform_to_biallelic, steps))]
pub(crate) fn calc_kinship<'py>(
    py: Python<'py>,
    source: &Bound<'_, PyAny>,
    individuals: Option<Vec<String>>,
    transform_to_biallelic: bool,
    steps: &Bound<'_, Steps>,
) -> Result<KinshipForPython<'py>, PyPopneiError> {
    let source = source_of(source)?;
    let steps = steps.get().of_a_pass()?;
    // A Ctrl-C that was pending when this was called is raised here, before
    // the file is opened.
    py.check_signals()?;
    let path = source.path().to_path_buf();
    // The whole source is read inside this one call, and the products of the
    // blocks are matrix work, so the interpreter is released for all of it:
    // the threads of rayon deadlock on a caller that holds it. A Ctrl-C that
    // arrives meanwhile is raised when the call is over and not between two
    // blocks, as it is in `Blocks::__next__`: the loop over the blocks is
    // the core's, and a pass that is interrupted loses only itself, since it
    // writes no file and the `Variants` is as it was.
    let calculated = py.detach(|| {
        over_the_source(
            source,
            &steps,
            individuals.as_deref(),
            transform_to_biallelic,
        )
    });
    // The file of the source goes into every error of the pass, and
    // `errors.rs` is what leaves it out of the message of the ones that are
    // of an argument a user wrote.
    let (kinship, individuals, counts) =
        calculated.map_err(|error| PyPopneiError::of_the_file(error, &path))?;
    // The Ctrl-C that arrived while the interpreter was released is raised
    // before numpy is called: the first array of a process imports the C API
    // of numpy, that import fails with the exception that is pending, and
    // the numpy crate panics when it does, which a user cannot catch.
    py.check_signals()?;
    // The variants that were used go to Python as the `u64` the core counts
    // them in, which is what every count of popnei is there: a `usize` is 32
    // bits in WebAssembly and 64 natively, and what a user reads does not
    // depend on that.
    let num_vars = kinship.num_vars;
    // The matrix of 10000 individuals is 800 MB, and `into_pyarray` hands
    // the allocation the core filled to numpy without copying it.
    let matrix = the_square_of(py, kinship.num_individuals, kinship.matrix)?;
    Ok((matrix, individuals, num_vars, counts))
}

/// One pass over `source` through `steps`, and the kinship of the
/// individuals it gives, or of the ones `individuals` names.
///
/// The interpreter is already released here, so nothing of Python is
/// touched: what comes back is Rust.
///
/// # Errors
///
/// When the source cannot be opened or read, when the pass gives no
/// variant, when a name of `individuals` is not an individual of the pass,
/// is there twice or there is none at all, and when the core refuses the
/// calculation: no variant with variance among those individuals, a pair
/// with no variant called in both of them, a variant of more than two
/// alleles while `transform_to_biallelic` is false, a ploidy or a number of
/// individuals the calculation does not count in, and a product of the
/// linear algebra that did not run.
fn over_the_source(
    source: &dyn OpenSource,
    steps: &[Step],
    individuals: Option<&[String]>,
    transform_to_biallelic: bool,
) -> popnei::Result<TheKinshipOfThePass> {
    // The source is opened at the size of its own blocks, since the core
    // puts a `reblock` over whatever it is given.
    let reader = source.reader(None)?;
    // The chain of the pass stays here, lent to the core, so that the counts
    // of its filters can be read when the call is over: the loop over the
    // blocks is the core's, and no block of it reaches this crate.
    let chain = chain_of(reader, steps)?;
    let mut counted = CountedVars::over(chain);
    // The names of the individuals of the matrix, which the package puts on
    // both sides of its frame: those the pass gives, in its order, or the
    // ones the user named, in theirs, which is the order the core has them
    // in.
    let (positions, names) = match individuals {
        Some(named) => (
            Some(resolve_individuals(named, counted.individuals())?),
            named.to_vec(),
        ),
        None => (None, counted.individuals().to_vec()),
    };
    let kinship =
        popnei::kinship::calc_kinship(&mut counted, positions.as_deref(), transform_to_biallelic)?;
    let num_vars = counted.num_vars();
    Ok((
        kinship,
        names,
        (num_vars, filtering_of(counted.of_the_pass())),
    ))
}

/// The chain of readers of a pass with how many variants it has given
/// counted, which is the `num_vars` of the counts of the pass.
///
/// The kinship gives back how many variants had variance and were used, and
/// a user reads that under `num_vars` of the result. How many the steps let
/// through, which is what the counts of a pass hold, is no number of the
/// core's result, so it is counted here, where the chain is built and lent:
/// the blocks of the pass go through this reader on their way to the core.
struct CountedVars {
    /// The chain the pass reads, this reader's source.
    chain: Box<dyn BlockReader>,
    /// How many variants that chain has given so far.
    num_vars: u64,
}

impl CountedVars {
    /// The chain with its variants counted, before the pass starts.
    fn over(chain: Box<dyn BlockReader>) -> Self {
        Self { chain, num_vars: 0 }
    }

    /// How many variants the chain has given, which after the pass is how
    /// many the steps let through.
    fn num_vars(&self) -> u64 {
        self.num_vars
    }

    /// The chain itself, whose filters are read when the pass is over.
    fn of_the_pass(&self) -> &dyn BlockReader {
        self.chain.as_ref()
    }
}

impl BlockReader for CountedVars {
    fn next_block(&mut self) -> popnei::Result<Option<Block>> {
        let block = self.chain.next_block()?;
        if let Some(ref block) = block {
            // A pass would have to give 18446744073709551615 variants to
            // reach the largest count, which at one variant a nanosecond is
            // 585 years of reading: the count is saturated rather than
            // carried back as an error of its own, which the error type of
            // the core has no case for, and no run of popnei arrives there.
            self.num_vars = self
                .num_vars
                .saturating_add(u64::try_from(block.num_vars).unwrap_or(u64::MAX));
        }
        Ok(block)
    }

    fn individuals(&self) -> &[String] {
        self.chain.individuals()
    }

    fn ploidy(&self) -> usize {
        self.chain.ploidy()
    }

    fn chroms(&self) -> &ChromTable {
        self.chain.chroms()
    }

    fn set_needs(&mut self, needs: Needs) {
        self.chain.set_needs(needs);
    }

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        self.chain.filtering_stats()
    }
}

/// The matrix of the kinship as a numpy array of individuals x individuals.
///
/// `values` is the `Vec` the core filled and gave away, and numpy takes it
/// over as it takes the distances of `dists.rs`: nothing is copied, so the
/// 800 MB of 10000 individuals are held once and not twice.
///
/// # Errors
///
/// [`PyPopneiError::Broken`] when the core gave a matrix whose values are
/// not the square of its individuals, which is a defect of popnei: a user
/// reports it instead of looking for what they typed wrong.
fn the_square_of(
    py: Python<'_>,
    num_individuals: usize,
    values: Vec<f64>,
) -> Result<Bound<'_, PyArray2<f64>>, PyPopneiError> {
    let num_values = values.len();
    let square =
        Array2::from_shape_vec((num_individuals, num_individuals), values).map_err(|error| {
            PyPopneiError::Broken {
                message: format!(
                    "the kinship of {num_individuals} individuals holds {num_values} \
                 values: {error}"
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
