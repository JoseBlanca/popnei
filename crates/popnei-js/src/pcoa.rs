//! What a TypeScript user reaches through `doPcoa` and
//! `correctDistsByLingoes`: the principal coordinates of the distances of a
//! `Distances`, and Lingoes' correction of them.
//!
//! The distance vector crosses as one `Float64Array`, which the code
//! wasm-bindgen generates copies into the memory of wasm, and the core takes
//! that copy over: the principal coordinates drop it once the matrix they
//! decompose is built. The names of the individuals cross beside it, for the
//! messages of the pairs with no distance and of a distance that is negative
//! or infinite, which the core writes with their positions. Both calculations are the core crate's and nothing of either
//! is here.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::pca::{LingoesCorrection, Pcoa};

use crate::errors::JsPopneiError;
use crate::pca::room_for_the_principal_coordinates_of;

/// The principal coordinates of a `Distances`, on their way to TypeScript.
///
/// Each array leaves the memory of wasm the first time it is asked for, and
/// the call after that gives nothing, as the arrays of the principal
/// components do: the package reads each of them once and frees this.
#[wasm_bindgen]
pub struct PcoaOfDistances {
    num_comps: usize,
    projections: Option<Vec<f64>>,
    explained_variance_percent: Option<Vec<f64>>,
    lingoes_constant: f64,
    negative_eigenvalues_percent: f64,
}

#[wasm_bindgen]
impl PcoaOfDistances {
    /// How many eigenvalues of the matrix of the squared distances are
    /// positive, which is how many components are given.
    #[must_use]
    pub fn num_comps(&self) -> usize {
        self.num_comps
    }

    /// Where each individual falls along each component, the individuals x
    /// `num_comps`, row after row.
    pub fn projections(&mut self) -> Option<Vec<f64>> {
        self.projections.take()
    }

    /// The variance of each component as a percentage of the sum of every
    /// eigenvalue, one number per component.
    pub fn explained_variance_percent(&mut self) -> Option<Vec<f64>> {
        self.explained_variance_percent.take()
    }

    /// The constant of Lingoes' correction, which is 0 here: the principal
    /// coordinates of a `Distances` correct nothing.
    #[must_use]
    pub fn lingoes_constant(&self) -> f64 {
        self.lingoes_constant
    }

    /// The share of the negative eigenvalues before a correction, which is 0
    /// here: distances that have any are refused.
    #[must_use]
    pub fn negative_eigenvalues_percent(&self) -> f64 {
        self.negative_eigenvalues_percent
    }
}

impl From<Pcoa> for PcoaOfDistances {
    /// The arrays and the two numbers of the analysis, without the count of
    /// the individuals, which the caller gave.
    fn from(result: Pcoa) -> PcoaOfDistances {
        PcoaOfDistances {
            num_comps: result.num_comps,
            projections: Some(result.projections),
            explained_variance_percent: Some(result.explained_variance_percent),
            lingoes_constant: result.lingoes_constant,
            negative_eigenvalues_percent: result.negative_eigenvalues_percent,
        }
    }
}

/// Lingoes' correction of a `Distances`, on its way to TypeScript.
///
/// The corrected vector leaves the memory of wasm the first time it is asked
/// for, and the package builds a `Distances` over it with the names and the
/// counts of the one given.
#[wasm_bindgen]
pub struct LingoesCorrectionOfDistances {
    dist_vector: Option<Vec<f64>>,
    constant: f64,
    negative_eigenvalues_percent: f64,
}

#[wasm_bindgen]
impl LingoesCorrectionOfDistances {
    /// sqrt(d² + 2c) for every distance d of the vector given, in its order,
    /// or `undefined` when it was read already.
    pub fn dist_vector(&mut self) -> Option<Vec<f64>> {
        self.dist_vector.take()
    }

    /// c, the absolute value of the most negative eigenvalue of the matrix of
    /// the squared distances, and 0 when the distances were Euclidean.
    #[must_use]
    pub fn constant(&self) -> f64 {
        self.constant
    }

    /// 100 times the sum of the absolute values of the negative eigenvalues
    /// over the sum of every eigenvalue, of the distances given.
    #[must_use]
    pub fn negative_eigenvalues_percent(&self) -> f64 {
        self.negative_eigenvalues_percent
    }
}

impl From<LingoesCorrection> for LingoesCorrectionOfDistances {
    fn from(correction: LingoesCorrection) -> LingoesCorrectionOfDistances {
        LingoesCorrectionOfDistances {
            dist_vector: Some(correction.dist_vector),
            constant: correction.constant,
            negative_eigenvalues_percent: correction.negative_eigenvalues_percent,
        }
    }
}

/// The principal coordinates of `dist_vector`, the distance of every pair of
/// the individuals `names` in the order (0, 1), (0, 2), ..., (1, 2), ...,
/// with NaN for a pair that has no distance.
///
/// The package checks that the vector holds one distance for each pair of
/// the names, which the `Distances` it comes from checked when it was built,
/// so the error the core has for a vector of another length is not reached
/// from TypeScript.
///
/// # Errors
///
/// When the analysis of that many individuals does not fit in the memory of
/// a page, which is above 8695 of them; when there are fewer than 2
/// individuals; when a pair has no distance, whose message names the
/// individuals; when a distance is negative or infinite, whose message names
/// the two individuals; when every distance
/// is 0; when the distances are not Euclidean, whose message names
/// `correctDistsByLingoes`; and when the linear algebra of the analysis could
/// not be done.
#[wasm_bindgen]
pub fn pcoa(dist_vector: Vec<f64>, names: Vec<String>) -> Result<PcoaOfDistances, JsPopneiError> {
    room_for_the_principal_coordinates_of(names.len())?;
    popnei::pca::pcoa(dist_vector, names.len())
        .map(PcoaOfDistances::from)
        .map_err(|error| under_the_names_of_the_individuals(error, &names))
}

/// Lingoes' correction of `dist_vector`, the distance of every pair of the
/// individuals `names` in the order (0, 1), (0, 2), ..., (1, 2), ..., with
/// NaN for a pair that has no distance.
///
/// # Errors
///
/// Those of [`pcoa`] but the distances that are not Euclidean, which this
/// corrects, and when the constant of the correction, which is in the units
/// of the squared distances, is not a normal 64 bit float.
#[wasm_bindgen]
pub fn correct_dists_by_lingoes(
    dist_vector: Vec<f64>,
    names: Vec<String>,
) -> Result<LingoesCorrectionOfDistances, JsPopneiError> {
    room_for_the_principal_coordinates_of(names.len())?;
    popnei::pca::correct_dists_by_lingoes(dist_vector, names.len())
        .map(LingoesCorrectionOfDistances::from)
        .map_err(|error| under_the_names_of_the_individuals(error, &names))
}

/// `error` with the pairs with no distance, or the pair whose distance is
/// negative or infinite, under the names of their individuals, and every
/// other error as it is.
///
/// The core names the individuals by their positions in the order of the
/// distances, which is what it has, and what a user takes out of a
/// `Distances`, or of a `Variants`, is a name. A position the names do not
/// reach is left as the core wrote it rather than named wrongly; the core
/// counts the individuals from the length of the names, so none is.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "the two errors that name individuals by their positions have arms of \
              their own, and any other, one the core adds later among them, crosses \
              as the core wrote it"
)]
pub(crate) fn under_the_names_of_the_individuals(
    error: popnei::Error,
    names: &[String],
) -> JsPopneiError {
    match error {
        popnei::Error::PcoaPairsWithNoDistance {
            num_pairs_with_no_distance,
            num_pairs,
            first_of_the_first,
            second_of_the_first,
            most_often,
            most_often_count,
            from,
        } => match (
            names.get(first_of_the_first),
            names.get(second_of_the_first),
            names.get(most_often),
        ) {
            (Some(first), Some(second), Some(most)) => JsPopneiError::PairsWithNoDistance {
                num_pairs_with_no_distance,
                num_pairs,
                first_of_the_first: first.clone(),
                second_of_the_first: second.clone(),
                most_often: most.clone(),
                most_often_count,
                from,
            },
            _ => JsPopneiError::Core(popnei::Error::PcoaPairsWithNoDistance {
                num_pairs_with_no_distance,
                num_pairs,
                first_of_the_first,
                second_of_the_first,
                most_often,
                most_often_count,
                from,
            }),
        },
        popnei::Error::PcoaDistanceOutOfRange {
            first,
            second,
            value,
        } => match (names.get(first), names.get(second)) {
            (Some(first), Some(second)) => JsPopneiError::DistanceOutOfRange {
                first: first.clone(),
                second: second.clone(),
                value,
            },
            _ => JsPopneiError::Core(popnei::Error::PcoaDistanceOutOfRange {
                first,
                second,
                value,
            }),
        },
        // Every other error of the core names no individual, and crosses
        // with the message the core gives it.
        other => JsPopneiError::Core(other),
    }
}
