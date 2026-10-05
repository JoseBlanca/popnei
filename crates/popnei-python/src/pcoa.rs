//! The principal coordinates of distances, and Lingoes' correction, on
//! their way between Python and the core crate.
//!
//! Both are given the distance vector of a `Distances` and how many
//! individuals it is of, which the Python package counts from its names.
//! The core takes the vector as a `Vec` of its own, which it drops before
//! the eigendecomposition in [`pcoa`], so the numpy array is copied into
//! one: the `Distances` of the user keeps its array, and the core cannot
//! drop memory that numpy holds.
//!
//! What goes back are arrays and numbers, which the Python package puts the
//! names on: the projections and the percentage of each component, or the
//! corrected vector, with the constant of the correction and the share of
//! the negative eigenvalues.
//!
//! The principal coordinates of the variants, [`pcoa_of_variants`], are one
//! pass over a source, and what this module does for them is what
//! `calc_pairwise_kosman_dists` of `dists.rs` does: it builds the chain of
//! readers of the pass from the steps of the `Variants`, lends it to the
//! core with the interpreter released, and reads the counts of the filters
//! from that chain when the call is over.

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1};
use pyo3::prelude::*;

use popnei::pca::{PcoaOfVariants, VariantPcoaOptions};

use crate::errors::{PyPopneiError, raise_a_ctrl_c_before_numpy_is_called};
use crate::source::{PassCounts, pass_counts_of, source_of};
use crate::steps::{Steps, chain_of};

/// A principal coordinate analysis as it goes to Python: the projections,
/// individuals x components; the percentage of the variance each component
/// holds; the constant of Lingoes' correction; and the share of the
/// negative eigenvalues of the distances before it.
type PcoaTables<'py> = (
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray1<f64>>,
    f64,
    f64,
);

/// A principal coordinate analysis of the variants as it goes to Python:
/// the tables of [`PcoaTables`] and the counts of the pass.
type VariantPcoaTables<'py> = (
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray1<f64>>,
    f64,
    f64,
    PassCounts,
);

/// Lingoes' correction as it goes to Python: the corrected distance vector,
/// the constant added, and the share of the negative eigenvalues of the
/// distances given.
type LingoesCorrectionForPython<'py> = (Bound<'py, PyArray1<f64>>, f64, f64);

// The principal coordinates of `dist_vector`, the distances of the pairs of
// `num_individuals` individuals in the order (0, 1), (0, 2), ..., (1, 2),
// .... A `///` comment here would become the `__doc__` of
// `popnei._core.pcoa`, and the documentation a Python user reads belongs to
// the package, which is the API.
#[pyfunction]
pub(crate) fn pcoa<'py>(
    py: Python<'py>,
    dist_vector: PyReadonlyArray1<'py, f64>,
    num_individuals: usize,
) -> Result<PcoaTables<'py>, PyPopneiError> {
    let dists = dist_vector.as_array();
    // The eigendecomposition of thousands of individuals takes seconds, and
    // the interpreter is of no use to it. The view is of an array the caller
    // holds and no Python code can reach while this runs; `to_vec` reads it
    // whatever its strides.
    let result = py.detach(|| popnei::pca::pcoa(dists.to_vec(), num_individuals))?;
    raise_a_ctrl_c_before_numpy_is_called(py)?;
    let projections = table_of(
        py,
        result.projections,
        result.num_individuals,
        result.num_comps,
    )?;
    Ok((
        projections,
        result.explained_variance_percent.into_pyarray(py),
        result.lingoes_constant,
        result.negative_eigenvalues_percent,
    ))
}

// Lingoes' correction of `dist_vector`, the distances of the pairs of
// `num_individuals` individuals in the order of `pcoa` above. A `///`
// comment here would become the `__doc__` of
// `popnei._core.correct_dists_by_lingoes`, for the reason `pcoa` gives.
#[pyfunction]
pub(crate) fn correct_dists_by_lingoes<'py>(
    py: Python<'py>,
    dist_vector: PyReadonlyArray1<'py, f64>,
    num_individuals: usize,
) -> Result<LingoesCorrectionForPython<'py>, PyPopneiError> {
    let dists = dist_vector.as_array();
    // The same eigendecomposition as `pcoa`, released for the same reason.
    let corrected =
        py.detach(|| popnei::pca::correct_dists_by_lingoes(dists.to_vec(), num_individuals))?;
    raise_a_ctrl_c_before_numpy_is_called(py)?;
    // The vector of 10000 individuals is 400 MB, and `into_pyarray` hands
    // the allocation of the core to numpy without copying it.
    Ok((
        corrected.dist_vector.into_pyarray(py),
        corrected.constant,
        corrected.negative_eigenvalues_percent,
    ))
}

// The principal coordinates of the Kosman distances of the individuals of
// `source`, over the variants that the steps of `steps` keep, with
// `min_num_vars` the variants a pair needs to get a distance and Lingoes'
// correction applied inside when `correct_by_lingoes` is true. A `///`
// comment here would become the `__doc__` of
// `popnei._core.pcoa_of_variants`, for the reason `pcoa` gives.
#[pyfunction]
#[pyo3(signature = (source, min_num_vars, correct_by_lingoes, steps))]
pub(crate) fn pcoa_of_variants<'py>(
    py: Python<'py>,
    source: &Bound<'_, PyAny>,
    min_num_vars: u32,
    correct_by_lingoes: bool,
    steps: &Bound<'_, Steps>,
) -> Result<VariantPcoaTables<'py>, PyPopneiError> {
    let source = source_of(source)?;
    let steps = steps.get().of_a_pass()?;
    // A Ctrl-C that was pending when this was called is raised here, before
    // the file is opened.
    py.check_signals()?;
    let path = source.path();
    let options = VariantPcoaOptions {
        min_num_vars,
        correct_by_lingoes,
    };
    // The pass over the whole source and the eigendecomposition are both
    // inside this one call, so the interpreter is released for all of it,
    // which also lets the core spread the pairs of a block over the threads
    // of rayon. A Ctrl-C that arrives meanwhile is raised when the call is
    // over, for the reason `calc_pairwise_kosman_dists` gives.
    let (result, counts) = py
        .detach(|| -> Result<(PcoaOfVariants, PassCounts), popnei::Error> {
            // The source is opened at the size of its own blocks: the sums
            // of the pass are whole numbers, the same whatever the size.
            // The chain stays here, lent to the core, so that the counts of
            // its filters can be read when the call is over.
            let mut chain = chain_of(source.reader(None)?, &steps)?;
            let result = popnei::pca::pcoa_of_variants(&mut chain, &options)?;
            let counts = pass_counts_of(result.num_vars, chain.as_ref(), &steps);
            Ok((result, counts))
        })
        .map_err(|error| PyPopneiError::of_the_file(error, path))?;
    raise_a_ctrl_c_before_numpy_is_called(py)?;
    let pcoa = result.pcoa;
    let projections = table_of(py, pcoa.projections, pcoa.num_individuals, pcoa.num_comps)?;
    Ok((
        projections,
        pcoa.explained_variance_percent.into_pyarray(py),
        pcoa.lingoes_constant,
        pcoa.negative_eigenvalues_percent,
        counts,
    ))
}

/// The projections of the result as a numpy array of `rows` x `cols`, which
/// takes the allocation of the core without copying it.
///
/// # Errors
///
/// [`PyPopneiError::Broken`] when the core gave projections that are not
/// its two sides multiplied, which is a defect of popnei: a user reports it
/// instead of looking for what they typed wrong.
fn table_of(
    py: Python<'_>,
    values: Vec<f64>,
    rows: usize,
    cols: usize,
) -> Result<Bound<'_, PyArray2<f64>>, PyPopneiError> {
    let num_values = values.len();
    let table =
        Array2::from_shape_vec((rows, cols), values).map_err(|error| PyPopneiError::Broken {
            message: format!(
                "the principal coordinate analysis gave projections of {rows} x {cols} that \
                 hold {num_values} values: {error}"
            ),
            path: None,
        })?;
    Ok(table.into_pyarray(py))
}
