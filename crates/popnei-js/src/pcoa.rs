//! What a TypeScript user reaches through `doPcoa`,
//! `correctDistsByLingoes` and `doPcoaFromVariants`: the principal
//! coordinates of the distances of a `Distances`, Lingoes' correction of
//! them, and the principal coordinates of the Kosman distances of the
//! individuals of a source, in one pass over its variants.
//!
//! The distance vector crosses as one `Float64Array`, which the code
//! wasm-bindgen generates copies into the memory of wasm, and the core takes
//! that copy over: the principal coordinates drop it once the matrix they
//! decompose is built. The names of the individuals cross beside it, for the
//! messages of the pairs with no distance and of a distance that is negative
//! or infinite, which the core writes with their positions. The variants
//! are not a vector that crosses: the core reads them block by block from
//! the chain of readers this crate opens over the source and the steps of
//! the `Variants`, as it does for `calcPairwiseKosmanDists`, and no vector
//! of distances is made. The three calculations are the core crate's and
//! nothing of them is here.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::pca::{LingoesCorrection, Pcoa, VariantPcoaOptions};

use crate::errors::JsPopneiError;
use crate::pca::room_for_the_principal_coordinates_of;
use crate::source::{Consumer, OpenSource, PassCounts, the_run_of};
use crate::steps::{Steps, chain_of};

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

/// The principal coordinates of the variants of a source, on their way to
/// TypeScript.
///
/// Each array leaves the memory of wasm the first time it is asked for, as
/// those of the principal coordinates of a `Distances` do. The names of the
/// individuals are not here: the `Variants` a user gave holds the names the
/// pass gives, and the package puts them on the result.
#[wasm_bindgen]
pub struct PcoaOfVariants {
    num_comps: usize,
    projections: Option<Vec<f64>>,
    explained_variance_percent: Option<Vec<f64>>,
    lingoes_constant: f64,
    negative_eigenvalues_percent: f64,
    pass_stats: PassCounts,
}

#[wasm_bindgen]
impl PcoaOfVariants {
    /// How many eigenvalues of the matrix of the squared distances, or of
    /// the corrected one, are positive, which is how many components are
    /// given.
    #[must_use]
    pub fn num_comps(&self) -> usize {
        self.num_comps
    }

    /// Where each individual falls along each component, the individuals of
    /// the pass x `num_comps`, row after row.
    pub fn projections(&mut self) -> Option<Vec<f64>> {
        self.projections.take()
    }

    /// The variance of each component as a percentage of the sum of every
    /// eigenvalue of the matrix decomposed, the corrected one when there was
    /// a correction, one number per component.
    pub fn explained_variance_percent(&mut self) -> Option<Vec<f64>> {
        self.explained_variance_percent.take()
    }

    /// c of Lingoes' correction, and 0 when it was not asked for or the
    /// distances were Euclidean.
    #[must_use]
    pub fn lingoes_constant(&self) -> f64 {
        self.lingoes_constant
    }

    /// The share of the negative eigenvalues of the distances before the
    /// correction, and 0 when it was not asked for, since distances with a
    /// negative eigenvalue are then refused.
    #[must_use]
    pub fn negative_eigenvalues_percent(&self) -> f64 {
        self.negative_eigenvalues_percent
    }

    /// How many variants the pass gave, called in a pair or not, and what
    /// each filter of it was given and kept.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.pass_stats.clone()
    }
}

/// The principal coordinates of the Kosman distances of the individuals of
/// `source`, over the variants the steps of `steps` keep, with no distance
/// for a pair called together at fewer than `min_num_vars` variants, and
/// Lingoes' correction inside when `correct_by_lingoes` is true.
///
/// The page is asked first, for the individuals the pass gives, which the
/// steps say: the matrix that is decomposed is of those, and a
/// `filterIndividuals` that keeps few of a large file makes a small one. No
/// run is opened for an analysis the page cannot hold, so the page is told
/// of no pass. The chain of readers stays here, lent to the core, so that
/// the counts of its filters are read when the analysis returns, and the
/// names of the individuals are the chain's own, taken before the core
/// borrows it, which the errors that name individuals by their positions
/// are written with.
///
/// # Errors
///
/// When the analysis of the individuals of the pass does not fit in the
/// memory of a page, which is above 8695 of them; when there are fewer than
/// 2, or more than the linear algebra decomposes the matrix of; the errors
/// of the pass of the Kosman distances, a pass that gave no variant, sums of
/// a pair beyond a `u32` and a source that cannot be read; when a pair has
/// no distance, whose message names the individuals; when every distance is
/// 0; when the distances are not Euclidean and `correct_by_lingoes` is
/// false, whose message names `correctByLingoes`; when the linear algebra
/// could not be done; and when the eigenvalue that the centering always
/// gives is not found, which is a defect of popnei.
pub(crate) fn pcoa_of_the_variants(
    source: &dyn OpenSource,
    min_num_vars: u32,
    correct_by_lingoes: bool,
    steps: Steps,
) -> Result<PcoaOfVariants, JsPopneiError> {
    room_for_the_principal_coordinates_of(steps.individuals().len())?;
    let options = VariantPcoaOptions {
        min_num_vars,
        correct_by_lingoes,
    };
    the_run_of(source, &Consumer::PcoaOfVariants, |run| {
        let reader = source.reader(run, None)?;
        let mut chain = chain_of(reader, steps.steps())?;
        let names = chain.individuals().to_vec();
        let result = popnei::pca::pcoa_of_variants(&mut chain, &options)
            .map_err(|error| under_the_names_of_the_individuals(error, &names))?;
        let pass_stats = PassCounts::of(result.num_vars, &chain.filtering_stats());
        Ok(PcoaOfVariants {
            num_comps: result.pcoa.num_comps,
            projections: Some(result.pcoa.projections),
            explained_variance_percent: Some(result.pcoa.explained_variance_percent),
            lingoes_constant: result.pcoa.lingoes_constant,
            negative_eigenvalues_percent: result.pcoa.negative_eigenvalues_percent,
            pass_stats,
        })
    })
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
