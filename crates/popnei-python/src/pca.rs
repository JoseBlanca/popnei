//! The principal component analysis of a table on its way between Python
//! and the core crate.
//!
//! The table arrives as the float64 array of the values of a pandas frame,
//! individuals x traits, which the Python package has made C contiguous
//! already, and it is read without a copy. What goes back are the three
//! tables of the result as arrays, which the Python package puts the names
//! of the frame on: the projections, the percentage of the variance of each
//! component and the weights of each trait in each component.

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray2, PyUntypedArrayMethods as _};
use pyo3::prelude::*;

use popnei::pca::PcaOptions;

use crate::errors::PyPopneiError;

/// The three tables of a principal component analysis as they go to Python:
/// the projections, individuals x components; the percentage of the
/// variance each component holds; and the weights, components x traits.
///
/// How many components there are is the second dimension of the
/// projections, and the names of the individuals, of the traits and of the
/// components are the Python package's: they are in the frame the user
/// gave, which this crate never sees.
type PcaTables<'py> = (
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray2<f64>>,
);

// The principal components of `data`, a C contiguous float64 array of
// individuals x traits, centered and standardized as the two flags say. A
// `///` comment here would become the `__doc__` of `popnei._core.pca`, and
// the documentation a Python user reads belongs to the package, which is
// the API.
#[pyfunction]
pub(crate) fn pca<'py>(
    py: Python<'py>,
    data: PyReadonlyArray2<'py, f64>,
    center_data: bool,
    standardize_data: bool,
) -> Result<PcaTables<'py>, PyPopneiError> {
    let (num_rows, num_cols) = data.as_array().dim();
    // The layout is asked of the array itself and not of `as_slice`, which
    // takes an array that lies column after column as well and gives its
    // values in that order: the core would read the table transposed and
    // give numbers for a table nobody has.
    if !data.is_c_contiguous() {
        return Err(PyPopneiError::ArrayNotContiguous { name: "data" });
    }
    let table = data
        .as_slice()
        .map_err(|_| PyPopneiError::ArrayNotContiguous { name: "data" })?;
    let options = PcaOptions {
        center: center_data,
        standardize: standardize_data,
    };
    // The decomposition of a table of thousands of individuals takes
    // seconds, and the interpreter is of no use to it: what goes into the
    // closure is the slice of the array, which the caller of this function
    // holds and which no Python code can reach while it runs.
    let result = py.detach(|| popnei::pca::pca(table, num_rows, num_cols, &options))?;
    let num_comps = result.num_comps;
    let num_used = result.used_cols.len();
    let projections = table_of(py, result.projections, num_rows, num_comps)?;
    let explained_variance_percent = result.explained_variance_percent.into_pyarray(py);
    let princomps = table_of(py, result.princomps, result.num_prin_comps, num_used)?;
    Ok((projections, explained_variance_percent, princomps))
}

/// The values of a table of the result as a numpy array of `rows` x `cols`,
/// which takes the allocation of the core without copying it.
///
/// # Errors
///
/// [`PyPopneiError::Broken`] when the core gave a table whose values are
/// not its two sides multiplied, which is a defect of popnei: a user
/// reports it instead of looking for what they typed wrong.
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
                "the principal component analysis gave a table of {rows} x {cols} that \
                 holds {num_values} values: {error}"
            ),
            path: None,
        })?;
    Ok(table.into_pyarray(py))
}
