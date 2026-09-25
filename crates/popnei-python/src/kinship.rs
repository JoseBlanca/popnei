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
//! What this module adds to the result of the core is the names of the
//! individuals of the matrix, which the package puts on both sides of its
//! frame: they are those the pass gives, in its order, or the ones a user
//! named, in theirs. The two counts are the core's, `num_vars`, the
//! variants that had variance and were used, and `num_vars_given`, the
//! variants the pass gave, which is the `num_vars` of the counts of a pass.
//!
//! The names a user names are turned into their places among the
//! individuals of the pass here, with `popnei::filters::resolve_individuals`,
//! which is what the filter of individuals is given as well: the core takes
//! the places and knows nothing of the names.
//!
//! [`kinship_principal_components`] is the other half of the module: it
//! places each individual along the directions in which the panel varies
//! most, out of a matrix alone. It reads no source, since the matrix a user
//! holds is the whole input, and a user who built that matrix by hand gets
//! its components as one of a pass does.

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2, PyUntypedArrayMethods as _};
use pyo3::prelude::*;

use popnei::block::BlockReader;
use popnei::filters::resolve_individuals;
use popnei::kinship::Kinship;

use crate::errors::{PyPopneiError, raise_a_ctrl_c_before_numpy_is_called};
use crate::source::{OpenSource, PassCounts, count_of_at_least, source_of};
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
    raise_a_ctrl_c_before_numpy_is_called(py)?;
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
    let mut chain = chain_of(reader, steps)?;
    // The names of the individuals of the matrix, which the package puts on
    // both sides of its frame: those the pass gives, in its order, or the
    // ones the user named, in theirs, which is the order the core has them
    // in.
    let (positions, names) = match individuals {
        Some(named) => (
            Some(resolve_individuals(named, chain.individuals())?),
            named.to_vec(),
        ),
        None => (None, chain.individuals().to_vec()),
    };
    let kinship =
        popnei::kinship::calc_kinship(&mut chain, positions.as_deref(), transform_to_biallelic)?;
    // How many variants the pass gave, used or not, which the core counts
    // and which is the `num_vars` of the counts of the pass.
    let num_vars_given = kinship.num_vars_given;
    Ok((
        kinship,
        names,
        (num_vars_given, filtering_of(chain.as_ref())),
    ))
}

// The principal components of the kinship `matrix`, an individuals x
// individuals float64 array that lies row after row, `num_pcs` of them at
// most. What it gives back is where each individual falls along each
// component, individuals x the components that were given, and how many
// those are: a component whose eigenvalue is not above the tolerance of
// `docs/specs/pca.md` is not given, so a kinship with fewer components than
// were asked for gives the ones it has. A `///` comment here would become
// the `__doc__` of `popnei._core.kinship_principal_components`, and what a
// Python user reads belongs to the package, which is the API.
#[pyfunction]
#[pyo3(signature = (matrix, num_pcs))]
pub(crate) fn kinship_principal_components<'py>(
    py: Python<'py>,
    matrix: PyReadonlyArray2<'py, f64>,
    num_pcs: &Bound<'_, PyAny>,
) -> Result<(Bound<'py, PyArray2<f64>>, usize), PyPopneiError> {
    // A `num_pcs` of 0 is no components and is no error, as asking a
    // principal component analysis for none is not, so 0 is the fewest.
    let num_pcs = count_of_at_least("num_pcs", 0, num_pcs)?;
    let (num_rows, num_columns) = matrix.as_array().dim();
    if num_rows != num_columns {
        return Err(PyPopneiError::MatrixNotSquare {
            name: "matrix",
            num_rows,
            num_columns,
        });
    }
    // The layout is asked of the array itself and not of `as_slice`, which
    // takes an array that lies column after column as well: the core would
    // read the upper half of the matrix as its lower half, which for a
    // matrix that is symmetric only within a tolerance is other numbers.
    if !matrix.is_c_contiguous() {
        return Err(PyPopneiError::ArrayNotContiguous { name: "matrix" });
    }
    let values = matrix
        .as_slice()
        .map_err(|_| PyPopneiError::ArrayNotContiguous { name: "matrix" })?;
    // The values are copied while the interpreter is held, since the array
    // they are in belongs to Python and the closure below has to own what
    // it reads. That one copy is what the components work in: they take the
    // matrix over, and the eigendecomposition writes the eigenvectors over
    // it. A `Kinship` built here would carry two counts nobody gave and the
    // matrix would be copied a second time to protect a kinship that is
    // thrown away.
    let matrix = values.to_vec();
    // The eigendecomposition of a matrix of thousands of individuals takes
    // seconds and the interpreter is of no use to it.
    let pcs = py.detach(|| popnei::kinship::principal_components_of(matrix, num_rows, num_pcs))?;
    raise_a_ctrl_c_before_numpy_is_called(py)?;
    let num_comps = pcs.num_comps;
    let projections = the_projections_of(py, pcs.projections, num_rows, num_comps)?;
    Ok((projections, num_comps))
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

/// The projections as a numpy array of individuals x components, which
/// takes the allocation of the core without copying it.
///
/// # Errors
///
/// [`PyPopneiError::Broken`] when the core gave projections that are not
/// its individuals times its components, which is a defect of popnei: a
/// user reports it instead of looking for what they typed wrong.
fn the_projections_of(
    py: Python<'_>,
    values: Vec<f64>,
    num_individuals: usize,
    num_comps: usize,
) -> Result<Bound<'_, PyArray2<f64>>, PyPopneiError> {
    let num_values = values.len();
    let table = Array2::from_shape_vec((num_individuals, num_comps), values).map_err(|error| {
        PyPopneiError::Broken {
            message: format!(
                "the {num_comps} principal components of a kinship of \
                 {num_individuals} individuals hold {num_values} values: {error}"
            ),
            path: None,
        }
    })?;
    Ok(table.into_pyarray(py))
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
