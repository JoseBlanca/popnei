//! What a TypeScript user reaches through `doPca`: the principal components
//! of a table of numbers that the page holds, individuals x traits.
//!
//! The table crosses as one `Float64Array`, row after row, which the code
//! wasm-bindgen generates copies into the memory of wasm before any code of
//! popnei runs, and the arrays of the result cross back the same way. The
//! analysis is the `pca` of the core crate and nothing of it is here.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::pca::{Pca, PcaOptions};

use crate::errors::JsPopneiError;

/// The principal components of a table, on their way to TypeScript.
///
/// Each array leaves the memory of wasm the first time it is asked for, and
/// the call after that gives nothing: the package reads each of them once,
/// into the object a user holds, and frees this. A copy left behind would
/// grow the memory of wasm, which never gives memory back, by the weights of
/// every column a second time.
#[wasm_bindgen]
pub struct PcaOfATable {
    num_comps: usize,
    projections: Option<Vec<f64>>,
    explained_variance_percent: Option<Vec<f64>>,
    num_prin_comps: usize,
    princomps: Option<Vec<f64>>,
}

#[wasm_bindgen]
impl PcaOfATable {
    /// How many components have variance, which is how many are given.
    #[must_use]
    pub fn num_comps(&self) -> usize {
        self.num_comps
    }

    /// Where each row falls along each component, the rows of the table x
    /// `num_comps`, row after row.
    pub fn projections(&mut self) -> Option<Vec<f64>> {
        self.projections.take()
    }

    /// The variance of each component as a percentage of the variance of
    /// every component of the data, the ones with no variance counted, one
    /// number per component.
    pub fn explained_variance_percent(&mut self) -> Option<Vec<f64>> {
        self.explained_variance_percent.take()
    }

    /// How many components the weights are given for, which for a table is
    /// `num_comps`.
    #[must_use]
    pub fn num_prin_comps(&self) -> usize {
        self.num_prin_comps
    }

    /// The weight of each trait in each component, `num_prin_comps` x the
    /// traits of the table, row after row.
    pub fn princomps(&mut self) -> Option<Vec<f64>> {
        self.princomps.take()
    }
}

impl From<Pca> for PcaOfATable {
    /// The three arrays a table's analysis gives, without the two sides of
    /// the table, which the caller wrote, and without `used_cols`, which for
    /// a table is every trait.
    fn from(result: Pca) -> PcaOfATable {
        PcaOfATable {
            num_comps: result.num_comps,
            projections: Some(result.projections),
            explained_variance_percent: Some(result.explained_variance_percent),
            num_prin_comps: result.num_prin_comps,
            princomps: Some(result.princomps),
        }
    }
}

/// The principal components of the table `data`, `num_rows` rows of
/// `num_cols` values each, row after row, the rows being the individuals and
/// the columns the traits.
///
/// `center_data` takes the mean of each trait from it and `standardize_data`
/// then divides each trait by its standard deviation; standardizing without
/// centering is an error. No value may be missing.
///
/// The package checks that `data` holds exactly `num_rows` times `num_cols`
/// values before this is called, so the error the core has for a buffer that
/// is not of its table is not reached from TypeScript.
///
/// # Errors
///
/// When a value of the table is an infinity or a NaN, when the table is to
/// be standardized and not centered, when it has fewer than 2 rows or no
/// traits, when it is standardized and a trait has no variance, when no
/// trait of it has any, when the mean or the standard deviation of a trait
/// is not a number the analysis can use, and when the linear algebra of the
/// analysis could not be done.
#[wasm_bindgen]
pub fn pca(
    data: Vec<f64>,
    num_rows: usize,
    num_cols: usize,
    center_data: bool,
    standardize_data: bool,
) -> Result<PcaOfATable, JsPopneiError> {
    let options = PcaOptions {
        center: center_data,
        standardize: standardize_data,
    };
    Ok(popnei::pca::pca(&data, num_rows, num_cols, &options)?.into())
}

/// Whether the mean of each trait is taken from it when the caller says
/// nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_center_data() -> bool {
    popnei::pca::DEFAULT_CENTER_DATA
}

/// Whether each trait is divided by its standard deviation when the caller
/// says nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_standardize_data() -> bool {
    popnei::pca::DEFAULT_STANDARDIZE_DATA
}
