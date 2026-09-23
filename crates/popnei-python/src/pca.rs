//! The two principal component analyses on their way between Python and
//! the core crate.
//!
//! The table of [`pca`] arrives as the float64 array of the values of a
//! pandas frame, individuals x traits, which the Python package has made C
//! contiguous already, and it is read without a copy. [`pca_of_variants`]
//! is given a source of variants and the steps that were put on it, and
//! opens a reader of its own for each of the two passes over them, as
//! `write_vars` does for its one pass, so that the counts of the filters
//! can be read when the core is done.
//!
//! What goes back are the tables of the result as arrays, which the Python
//! package puts the names on: the projections, the percentage of the
//! variance of each component and the weights of each trait, or of each
//! variant that was used, in each component.

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray2, PyUntypedArrayMethods as _};
use pyo3::prelude::*;

use popnei::block::BlockReader;
use popnei::pca::{Pca, PcaOptions, VariantPcaOptions};

use crate::errors::{PyPopneiError, raise_a_ctrl_c_before_numpy_is_called};
use crate::source::{PassCounts, count_of_at_least, source_of};
use crate::steps::{Steps, chain_of};

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
    raise_a_ctrl_c_before_numpy_is_called(py)?;
    let num_comps = result.num_comps;
    let num_used = result.used_cols.len();
    let projections = table_of(py, result.projections, num_rows, num_comps)?;
    let explained_variance_percent = result.explained_variance_percent.into_pyarray(py);
    let princomps = table_of(py, result.princomps, result.num_prin_comps, num_used)?;
    Ok((projections, explained_variance_percent, princomps))
}

/// What each filter of a pass was given and kept, under the kind of the
/// filter and in the order of the chain, which is the half of the counts of
/// a pass that the chain of readers holds.
type FilteringCounts = Vec<(&'static str, u64, u64)>;

/// The tables of a principal component analysis of the variants as they go
/// to Python: the projections, individuals x components; the percentage of
/// the variance each component holds; the weights, `num_prin_comps` x the
/// variants that were used; the position of each of those variants among
/// the variants the pass gave, from 0; and the counts of the pass.
///
/// The names of the individuals are the Python package's, which has them
/// from the `Variants` its user holds, and so are the names of the
/// components.
type VariantPcaTables<'py> = (
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray1<i64>>,
    PassCounts,
);

// The principal components of the variants of `source` after the `steps`,
// with the weights of the first `num_prin_comps` components. A `///`
// comment here would become the `__doc__` of
// `popnei._core.pca_of_variants`, and the documentation a Python user reads
// belongs to the package, which is the API.
#[pyfunction]
pub(crate) fn pca_of_variants<'py>(
    py: Python<'py>,
    source: &Bound<'_, PyAny>,
    transform_to_biallelic: bool,
    num_prin_comps: &Bound<'_, PyAny>,
    steps: &Bound<'_, Steps>,
) -> Result<VariantPcaTables<'py>, PyPopneiError> {
    let source = source_of(source)?;
    // 0 is no weights and no second pass, so this is the one count of the
    // crate whose smallest is 0 and not 1.
    let num_prin_comps = count_of_at_least("num_prin_comps", 0, num_prin_comps)?;
    let steps = steps.get().of_a_pass()?;
    let path = source.path();
    let options = VariantPcaOptions {
        transform_to_biallelic,
        num_prin_comps,
    };
    // The whole source is read inside this one call, twice when the weights
    // are asked for, which is seconds for a VCF of hundreds of megabytes,
    // so the interpreter is released for all of it.
    //
    // A Ctrl-C that arrives meanwhile is raised when the call is over and
    // not between two blocks, as it is in `Blocks::__next__`: the loop over
    // the blocks is the core's, which `docs/specs/pca.md` has this crate
    // call instead of writing that loop again, and the core has nothing to
    // ask a caller between two blocks. It bites harder here than in
    // `write_vars`: this is minutes of work of the processor on a million
    // variants, the products of the blocks and the eigendecomposition, and
    // not the reading of a file that a disc paces, and a user who wants it
    // stopped waits for it to end. Whether the core takes a callback for
    // that is the owner's to decide and is not in the plan of this module.
    let (result, filtering) = py
        .detach(|| -> Result<(Pca, FilteringCounts), popnei::Error> {
            // The chain of each pass stays here, lent to the core, so that
            // the counts of the filters of the first can be read when the
            // call is over: the loop over the blocks is the core's. The
            // source is opened at the size of its own blocks, since the
            // core puts a `reblock` over whatever it is given.
            let mut first_pass = chain_of(source.reader(None)?, &steps)?;
            // The weight of a variant needs the eigenvectors, which are
            // known when the first pass ends, so it comes from a second
            // pass over the same variants. With no weights asked for there
            // is no second pass and no second reader, which means no second
            // reading of the file.
            let mut second_pass = match num_prin_comps {
                0 => None,
                _ => Some(chain_of(source.reader(None)?, &steps)?),
            };
            let result =
                popnei::pca::pca_of_variants(&mut first_pass, second_pass.as_mut(), &options)?;
            // The two passes go through the same steps and count the same,
            // so the counts are those of the first.
            let filtering = first_pass
                .filtering_stats()
                .into_iter()
                .map(|(kind, stats)| (kind, stats.vars_processed, stats.vars_kept))
                .collect();
            Ok((result, filtering))
        })
        .map_err(|error| PyPopneiError::of_the_file(error, path))?;
    raise_a_ctrl_c_before_numpy_is_called(py)?;
    // The variants the pass gave, used or not, which is the `num_vars` of
    // its counts.
    let num_vars: u64 = counted_for_python(result.num_cols, path)?;
    // The positions of the variants that were used go as signed numbers:
    // they become the columns of a frame, where pandas and pyNei hold them
    // as int64, and a user who takes one number of them from another gets
    // -1 for 0 - 1 and not 18446744073709551615.
    let used_vars = result
        .used_cols
        .iter()
        .map(|position| counted_for_python::<i64>(*position, path))
        .collect::<Result<Vec<i64>, PyPopneiError>>()?;
    let projections = table_of(py, result.projections, result.num_rows, result.num_comps)?;
    let explained_variance_percent = result.explained_variance_percent.into_pyarray(py);
    let princomps = table_of(
        py,
        result.princomps,
        result.num_prin_comps,
        result.used_cols.len(),
    )?;
    Ok((
        projections,
        explained_variance_percent,
        princomps,
        used_vars.into_pyarray(py),
        (num_vars, filtering),
    ))
}

/// A count of the variants of a pass as the number Python reads, which is
/// the same width on every platform: a `usize` is 32 bits in WebAssembly
/// and 64 natively, and what a user gets does not depend on that.
///
/// The counts of the pass go as `u64`, which is what every count of popnei
/// is in Python, and the positions of the variants that were used as `i64`,
/// which is what a frame's columns are in pandas.
///
/// # Errors
///
/// [`PyPopneiError::Broken`] when the count is above what the number holds,
/// which for an `i64` is 9.2e18 variants.
fn counted_for_python<T: TryFrom<usize>>(
    count: usize,
    path: &std::path::Path,
) -> Result<T, PyPopneiError> {
    T::try_from(count).map_err(|_| {
        PyPopneiError::broken_of_the_file(
            format!("the pass over the variants counted {count}, which is more than a count holds"),
            path,
        )
    })
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
