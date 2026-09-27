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

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1};
use pyo3::prelude::*;

use crate::errors::{PyPopneiError, raise_a_ctrl_c_before_numpy_is_called};

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
