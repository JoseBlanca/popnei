//! The principal coordinates of distances.
//!
//! A principal coordinate analysis, PCoA, places the individuals on
//! components as a principal component analysis does, from the distance of
//! every pair of them instead of a table of their values: the individuals
//! are put in a space where the straight line between each two of them is
//! as long as their distance, and the components are the directions of
//! that space along which they vary most. It is Gower's method (1966,
//! Biometrika 53: 325), which R's `cmdscale` and ape's `pcoa` compute.
//!
//! With d_ij the distance of the individuals i and j, and n the
//! individuals, the matrix that is decomposed is
//!
//! ```text
//! A_ij = -d_ij² / 2
//! B_ij = A_ij - (mean of row i of A) - (mean of column j of A) + (mean of A)
//! ```
//!
//! and with λ_j its eigenvalues from the largest and u_j its eigenvectors
//! of length 1, the projections of component j are u_j sqrt(λ_j), and its
//! share of the variance 100 λ_j over the sum of every eigenvalue of B.
//!
//! A matrix of distances is Euclidean when some space has points whose
//! straight line distances are those distances, and then no eigenvalue of
//! B is negative. [`pcoa`] refuses one that is not, since a negative
//! eigenvalue is a direction whose projections would be the square root of
//! a negative number.

use std::fmt;

use popnei_linalg::{Eigen, TheFirstOperand, TheSecondOperand, eigh_lower, product};

use crate::block::BlockReader;
use crate::dists::calc_kosman_sums;
use crate::error::{Error, Result};
use crate::pca::{fix_the_sign_of, the_percentages_of, write_the_projections};
use crate::variant::MAX_INDIVIDUALS_OF_THE_VARIANTS;

/// How far along the vector of ones, of length 1, a vector of length 1 may
/// be and still be taken for one at a right angle to it: 1e-6, the
/// tolerance the review of the correction proposed on 27 September 2026.
/// An eigenvector that is at a right angle to it in exact arithmetic comes
/// out along it by a rounding of about the individuals times 2.2e-16, 1e-11
/// at 46340 of them, and one along it by 1, so the two are five orders of
/// magnitude from it on either side. Nothing measured has come near it.
const ALONG_THE_VECTOR_OF_ONES: f64 = 1e-6;

/// What a principal coordinate analysis gives.
///
/// Only the components of the positive eigenvalues of B are here, and
/// `num_comps` is how many there are: an eigenvalue is positive when it is
/// above [`the_threshold_of_the_eigenvalues`], so the 0 that the centering
/// of B always has gives no component. In each component the projection of the largest
/// absolute value is positive, as in the principal components.
#[derive(Debug, Clone)]
pub struct Pcoa {
    /// How many individuals were placed.
    pub num_individuals: usize,
    /// How many eigenvalues of B, or of the corrected B, are positive,
    /// which is how many components are given.
    pub num_comps: usize,
    /// num_individuals x num_comps, row after row.
    pub projections: Vec<f64>,
    /// One per component, over the sum of every eigenvalue of the B
    /// decomposed, the corrected one when there was a correction.
    pub explained_variance_percent: Vec<f64>,
    /// c of Lingoes' correction, 0 when there was none, which is always
    /// from [`pcoa`], or the matrix was Euclidean.
    pub lingoes_constant: f64,
    /// 100 times the sum of |λ| of the negative eigenvalues of B before any
    /// correction over the sum of every eigenvalue of it.
    pub negative_eigenvalues_percent: f64,
}

/// Where the distances of a principal coordinate analysis came from, which
/// decides what the messages of its errors tell the user to do.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PcoaInput {
    /// A distance vector the user gave, the `Distances` of `do_pcoa` and
    /// of `correct_dists_by_lingoes`.
    Distances,
    /// The Kosman distances of a pass over the variants, of
    /// `do_pcoa_from_variants`.
    Variants,
}

/// The principal coordinates of the distance vector `dist_vector` of
/// `num_individuals` individuals, in the order (0, 1), (0, 2), ..., (1,
/// 2), ... of the distances.
///
/// It takes the vector by value and drops it once B is built, before the
/// eigendecomposition. It refuses a matrix that is not Euclidean, and its
/// result has a `lingoes_constant` and a `negative_eigenvalues_percent` of
/// 0.
///
/// # Errors
///
/// [`Error::PcoaTooManyIndividuals`] when there are more than
/// [`MAX_INDIVIDUALS_OF_THE_VARIANTS`] individuals.
/// [`Error::PcoaDistVectorOfAnotherSize`] when the vector does not hold
/// n(n - 1)/2 distances. [`Error::PcoaTooFewIndividuals`] when there are
/// fewer than 2 individuals. [`Error::PcoaPairsWithNoDistance`] when a
/// distance is a NaN, [`Error::PcoaDistanceOutOfRange`] when one is
/// negative or infinite, and [`Error::PcoaAllDistancesZero`] when every
/// one is 0. [`Error::PcoaNotEuclidean`] when B has a negative eigenvalue.
/// [`Error::PcoaLinalg`] when the eigendecomposition could not be done,
/// [`Error::PcoaNoMemory`] when the machine does not give the memory of B
/// or of the projections, and [`Error::PcoaComponentAlongTheVectorOfOnes`]
/// when a component is not at a right angle to the vector of ones, which
/// is a defect of popnei.
pub fn pcoa(dist_vector: Vec<f64>, num_individuals: usize) -> Result<Pcoa> {
    refuse_a_vector_that_cannot_be_analysed(&dist_vector, num_individuals)?;
    let dists = || the_distances_of_the_vector(&dist_vector);
    let largest = the_largest_distance(dists, num_individuals, PcoaInput::Distances)?;
    let centered = the_centered_matrix(
        the_matrix_of(num_individuals)?,
        dists(),
        num_individuals,
        largest,
    );
    drop(dist_vector);
    the_analysis_of(
        centered,
        num_individuals,
        largest,
        WhenNotEuclidean::Refuse,
        PcoaInput::Distances,
    )
}

/// Lingoes' correction of a distance vector: the corrected vector, in the
/// order of the one given, the constant c that was added, and the
/// `negative_eigenvalues_percent` of the vector given.
#[derive(Debug, Clone)]
pub struct LingoesCorrection {
    /// sqrt(d² + 2c) for every distance d of the vector given, or that
    /// vector itself when c is 0.
    pub dist_vector: Vec<f64>,
    /// c, the absolute value of the most negative eigenvalue of B, and 0
    /// when none is below minus the threshold of the components.
    pub constant: f64,
    /// 100 times the sum of |λ| of the negative eigenvalues of B over the
    /// sum of every eigenvalue of it, of the vector given.
    pub negative_eigenvalues_percent: f64,
}

/// Lingoes' correction of the distance vector `dist_vector` of
/// `num_individuals` individuals, in the order (0, 1), (0, 2), ..., (1,
/// 2), ... of the distances, which makes a matrix of distances Euclidean.
///
/// It adds 2c to the square of every distance, c being the absolute value
/// of the most negative eigenvalue of B, which is ape's `pcoa(d, correction
/// = "lingoes")`: the B of the corrected distances has the eigenvectors of
/// B, and every eigenvalue but the 0 of the centering c larger, so the most
/// negative becomes 0 and none is below it. A Euclidean matrix gives a
/// constant of 0 and its own vector back.
///
/// The vector is kept through the eigendecomposition, since the corrected
/// one is written from it. [`pcoa`] drops its vector before and holds the
/// projections after, which are more than the vector.
///
/// # Errors
///
/// The errors of [`pcoa`] but [`Error::PcoaNotEuclidean`], which this
/// corrects, in the same order. [`Error::PcoaLingoesConstantOutOfRange`]
/// when c, which is in the units of the squared distances, is not a normal
/// `f64`: above the largest, 1.8e308, or below the smallest normal one,
/// 2.2e-308.
pub fn correct_dists_by_lingoes(
    dist_vector: Vec<f64>,
    num_individuals: usize,
) -> Result<LingoesCorrection> {
    refuse_a_vector_that_cannot_be_analysed(&dist_vector, num_individuals)?;
    let dists = || the_distances_of_the_vector(&dist_vector);
    let largest = the_largest_distance(dists, num_individuals, PcoaInput::Distances)?;
    let centered = the_centered_matrix(
        the_matrix_of(num_individuals)?,
        dists(),
        num_individuals,
        largest,
    );
    let decomposed = the_decomposition_of(centered, num_individuals)?;
    if decomposed.num_negative == 0 {
        return Ok(LingoesCorrection {
            dist_vector,
            constant: 0.0,
            negative_eigenvalues_percent: 0.0,
        });
    }
    // The eigenvalues are of the distances over the largest, so c over its
    // square: each distance is corrected over the largest as well, and
    // multiplied back, which keeps the square of a distance above 1.3e154
    // from being an infinity on the way.
    let constant_of_the_scaled = decomposed
        .eigen
        .values
        .last()
        .map_or(0.0, |most_negative| most_negative.abs());
    // c is in the units of a squared distance, which the corrected
    // distances are not, so distances of about 1.3e154 and above give a c
    // above the largest f64, an infinity, and distances all below about
    // 1e-154 one below the smallest normal f64, 2.2e-308, which keeps only a
    // few of its significant digits, or is 0 and would say that nothing was
    // corrected.
    let constant = constant_of_the_scaled * largest * largest;
    if !constant.is_normal() {
        return Err(Error::PcoaLingoesConstantOutOfRange { largest });
    }
    let mut corrected = dist_vector;
    for dist in &mut corrected {
        let scaled = *dist / largest;
        *dist = largest * (scaled * scaled + 2.0 * constant_of_the_scaled).sqrt();
    }
    Ok(LingoesCorrection {
        dist_vector: corrected,
        constant,
        negative_eigenvalues_percent: decomposed.negative_eigenvalues_percent,
    })
}

/// What the principal coordinates of the variants are asked for.
#[derive(Debug, Clone, Copy)]
pub struct VariantPcoaOptions {
    /// The `min_num_snps` of Python and TypeScript: a pair of individuals
    /// called together at fewer variants has no distance, and 0 gives one
    /// to every pair called together at one variant at least.
    pub min_num_vars: u32,
    /// Whether Lingoes' correction is applied inside the analysis to
    /// distances that are not Euclidean, which are refused without it.
    pub correct_by_lingoes: bool,
}

/// The principal coordinates of the variants, with how many variants the
/// pass gave, for the pass stats.
#[derive(Debug, Clone)]
pub struct PcoaOfVariants {
    /// The analysis.
    pub pcoa: Pcoa,
    /// How many variants the pass gave, called in a pair or not.
    pub num_vars: u64,
}

/// The principal coordinates of the Kosman distances of the individuals of
/// the variants of `reader`, in one pass.
///
/// The pass is `calc_kosman_sums` of the module `dists`, which borrows the
/// reader, so that the caller reads the counts of the filters from it
/// afterwards. B is built from the sums of each pair and no vector of
/// distances is made; the sums are dropped before the eigendecomposition.
/// A pair called together at fewer than `min_num_vars` variants, or at
/// none, has no distance. Distances that are not Euclidean are refused
/// unless `correct_by_lingoes` is asked for.
///
/// # Errors
///
/// [`Error::PcoaTooManyIndividuals`] and [`Error::PcoaTooFewIndividuals`],
/// before the pass, from the individuals the reader says its source has.
/// The errors of `calc_kosman_sums`: a pass that gave no variant, sums of
/// a pair beyond a `u32`, memory the machine does not give, and the
/// reader's own. Then [`Error::PcoaPairsWithNoDistance`] and
/// [`Error::PcoaAllDistancesZero`], on the sums; [`Error::PcoaNoMemory`] and
/// [`Error::PcoaLinalg`]; and [`Error::PcoaNotEuclidean`] when B has a
/// negative eigenvalue and the correction was not asked for.
pub fn pcoa_of_variants<R: BlockReader + ?Sized>(
    reader: &mut R,
    options: &VariantPcoaOptions,
) -> Result<PcoaOfVariants> {
    let num_individuals = reader.individuals().len();
    refuse_too_many_individuals(num_individuals)?;
    refuse_too_few_individuals(num_individuals)?;
    // B is asked of the machine before the pass, so that the sums, asked
    // for after it, lie above it in the memory of wasm and their room is
    // free again at the top when they are given back, and so that a B the
    // machine cannot hold is refused before the variants are read.
    let matrix = the_matrix_of(num_individuals)?;
    let sums = calc_kosman_sums(reader)?;
    let num_vars = sums.num_vars();
    let dists = || sums.dists(options.min_num_vars);
    let largest = the_largest_distance(dists, num_individuals, PcoaInput::Variants)?;
    let centered = the_centered_matrix(matrix, dists(), num_individuals, largest);
    drop(sums);
    let when_not_euclidean = if options.correct_by_lingoes {
        WhenNotEuclidean::Correct
    } else {
        WhenNotEuclidean::Refuse
    };
    let pcoa = the_analysis_of(
        centered,
        num_individuals,
        largest,
        when_not_euclidean,
        PcoaInput::Variants,
    )?;
    Ok(PcoaOfVariants { pcoa, num_vars })
}

/// What the analysis does with distances that are not Euclidean.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WhenNotEuclidean {
    /// Refuses them, which [`pcoa`] always does.
    Refuse,
    /// Corrects them by Lingoes inside the analysis, from the eigenvalues
    /// of B.
    Correct,
}

/// The principal coordinates of the centered matrix `centered` of the
/// distances divided by `largest`, which, when the distances are not
/// Euclidean, are refused or corrected by Lingoes inside the analysis as
/// `when_not_euclidean` says.
///
/// # Errors
///
/// [`Error::PcoaLinalg`], [`Error::PcoaNotEuclidean`] with `from`, the
/// errors of the correction, and [`Error::PcoaNoMemory`].
fn the_analysis_of(
    centered: Vec<f64>,
    num_individuals: usize,
    largest: f64,
    when_not_euclidean: WhenNotEuclidean,
    from: PcoaInput,
) -> Result<Pcoa> {
    let decomposed = the_decomposition_of(centered, num_individuals)?;
    if decomposed.num_negative == 0 {
        return the_components_of(&decomposed, num_individuals, largest);
    }
    if when_not_euclidean == WhenNotEuclidean::Refuse {
        return Err(Error::PcoaNotEuclidean {
            num_negative: decomposed.num_negative,
            num_individuals,
            negative_eigenvalues_percent: decomposed.negative_eigenvalues_percent,
            from,
        });
    }
    let negative_eigenvalues_percent = decomposed.negative_eigenvalues_percent;
    let (corrected, constant_of_the_scaled) =
        corrected_by_lingoes(decomposed.eigen, num_individuals)?;
    // The Kosman distances are at most 1, so c is an f64 whatever the
    // dataset; it is looked at all the same, as `correct_dists_by_lingoes`
    // does, since a c that is not a normal f64 would be a number with few
    // of its digits or one that reads as no correction.
    let lingoes_constant = constant_of_the_scaled * largest * largest;
    if !lingoes_constant.is_normal() {
        return Err(Error::PcoaLingoesConstantOutOfRange { largest });
    }
    let mut pcoa = the_components_of(&corrected, num_individuals, largest)?;
    pcoa.lingoes_constant = lingoes_constant;
    pcoa.negative_eigenvalues_percent = negative_eigenvalues_percent;
    Ok(pcoa)
}

/// The decomposition of the B of the distances corrected by Lingoes, from
/// the decomposition `eigen` of the B of the distances, and c.
///
/// The corrected B is B + cJ, J being the centering matrix, so it has the
/// eigenvectors of B, and every eigenvalue but the 0 of the centering is c
/// larger: the most negative becomes 0 and the others are above it. The
/// eigenvector of the 0 of the centering is the vector of ones, which a
/// decomposition tells apart from the others only when the eigenvalue 0 has
/// no other eigenvector: two clones, individuals at distance 0 from each
/// other and at the same distance from every other individual, give a
/// second one, and the decomposition may give any two directions of the
/// plane of the two. So the eigenvectors whose eigenvalues are 0 within
/// [`the_threshold_of_the_eigenvalues`] are projected orthogonal to the
/// vector of ones and made orthonormal again, which leaves one fewer, and
/// those take the eigenvalue c; the vector of ones takes the place of the
/// one left over, with the eigenvalue 0, at the end.
///
/// The eigenvalues come back from the largest, as `eigen` has them.
///
/// # Errors
///
/// [`Error::PcoaNoEigenvalueOfTheCentering`] when no eigenvalue of B is 0
/// within the threshold, which the centering always gives and which is
/// then a defect of popnei.
fn corrected_by_lingoes(mut eigen: Eigen, num_individuals: usize) -> Result<(Decomposed, f64)> {
    let side = num_individuals;
    let threshold = the_threshold_of_the_eigenvalues(&eigen.values, side);
    let constant = eigen
        .values
        .last()
        .map_or(0.0, |most_negative| most_negative.abs());
    // The eigenvalues come from the largest, so those that are 0 within
    // the threshold are one run of them.
    let band_start = eigen
        .values
        .iter()
        .take_while(|value| **value > threshold)
        .count();
    let band_size = eigen
        .values
        .iter()
        .skip(band_start)
        .take_while(|value| value.abs() <= threshold)
        .count();
    if band_size == 0 {
        return Err(Error::PcoaNoEigenvalueOfTheCentering {
            num_individuals,
            threshold,
        });
    }
    let band_end = band_start.saturating_add(band_size);
    let band_values = band_start
        .checked_mul(side)
        .zip(band_end.checked_mul(side))
        .and_then(|(from, to)| eigen.vectors.get_mut(from..to))
        .ok_or(Error::PcoaNoEigenvalueOfTheCentering {
            num_individuals,
            threshold,
        })?;
    orthogonal_to_the_vector_of_ones(band_values, side, num_individuals)?;
    for (at, value) in eigen.values.iter_mut().enumerate() {
        *value = if at < band_start || at >= band_end {
            *value + constant
        } else if at.saturating_add(1) < band_end {
            constant
        } else {
            0.0
        };
    }
    // The vector of ones, of the eigenvalue 0, goes after the eigenvalues
    // of B that were negative and are now c larger, 0 and above: one row
    // of the vectors and one value from the last of the band to the end.
    let last_of_the_band = band_end.saturating_sub(1);
    if let Some(values) = eigen.values.get_mut(last_of_the_band..) {
        values.rotate_left(1);
    }
    if let Some(vectors) = last_of_the_band
        .checked_mul(side)
        .and_then(|from| eigen.vectors.get_mut(from..))
    {
        vectors.rotate_left(side);
    }
    let corrected_threshold = the_threshold_of_the_eigenvalues(&eigen.values, side);
    let num_positive = eigen
        .values
        .iter()
        .take_while(|value| **value > corrected_threshold)
        .count();
    Ok((
        Decomposed {
            eigen,
            num_positive,
            num_negative: 0,
            negative_eigenvalues_percent: 0.0,
        },
        constant,
    ))
}

/// Makes the rows of `vectors`, each of `side` values, eigenvectors of
/// the corrected B: each is projected orthogonal to the vector of ones, the
/// mean of its values taken from it, and they are made orthonormal again,
/// which leaves one fewer. The first rows are the orthonormal ones and the
/// last is the vector of ones over the square root of `side`, of length 1.
///
/// The k rows V, once their means are out, span a space of k - 1
/// dimensions when the vector of ones was among the directions they
/// spanned. The eigendecomposition of the k x k matrix V V', whose entries
/// are the products of each two rows, gives it: its eigenvalues are about 1
/// for those k - 1 directions and about 0 for the one the vector of ones
/// took out, and the eigenvector w of each of the k - 1 gives the row w' V
/// over the square root of its eigenvalue, of length 1 and at a right
/// angle to the others. Both the product and the decomposition are those
/// of the crate `popnei-linalg`, which the core leaves its linear algebra
/// to.
///
/// # Errors
///
/// [`Error::PcoaBandWithoutTheVectorOfOnes`] when the rows do not hold the
/// vector of ones: the direction the means took the most from is longer
/// than [`ALONG_THE_VECTOR_OF_ONES`] once they are out.
/// [`Error::PcoaLinalg`] when the product or the decomposition could not
/// be done, and [`Error::PcoaNoMemory`] when the machine does not give the
/// memory of the new rows.
fn orthogonal_to_the_vector_of_ones(
    vectors: &mut [f64],
    side: usize,
    num_individuals: usize,
) -> Result<()> {
    let num_rows = vectors.len().checked_div(side).unwrap_or(0);
    for row in vectors.chunks_exact_mut(side) {
        let mean = row.iter().sum::<f64>() / side as f64;
        for value in row.iter_mut() {
            *value -= mean;
        }
    }
    let num_kept = num_rows.saturating_sub(1);
    let mut products = zeros_or_no_memory(
        num_rows.saturating_mul(num_rows),
        num_individuals,
        "eigenvectors of the eigenvalue 0",
    )?;
    product(
        TheFirstOperand::ByTheRowsOfTheResult {
            values: vectors,
            rows: num_rows,
        },
        side,
        TheSecondOperand::ByTheColumnsOfTheResult {
            values: vectors,
            cols: num_rows,
        },
        &mut products,
    )
    .map_err(|source| Error::PcoaLinalg {
        operation: "product of the eigenvectors of the eigenvalue 0 with themselves",
        source,
    })?;
    let directions = eigh_lower(products, num_rows).map_err(|source| Error::PcoaLinalg {
        operation: "eigendecomposition of the products of the eigenvectors of the eigenvalue 0",
        source,
    })?;
    // The direction of the band the means took the most from has the
    // smallest eigenvalue, the square of its length once they are out: 0
    // when the band held the vector of ones, and 1 when it did not, which
    // rounding that lifted the eigenvalue of the vector of ones out of the
    // band gives, and where the vector of ones would be written over a real
    // eigenvector.
    let length = directions
        .values
        .last()
        .map_or(f64::INFINITY, |smallest| smallest.max(0.0).sqrt());
    if length > ALONG_THE_VECTOR_OF_ONES {
        return Err(Error::PcoaBandWithoutTheVectorOfOnes {
            num_individuals,
            band_size: num_rows,
            length,
        });
    }
    let mut rows = zeros_or_no_memory(
        num_kept.saturating_mul(side),
        num_individuals,
        "eigenvectors of the eigenvalue 0",
    )?;
    let kept = directions
        .vectors
        .get(..num_kept.saturating_mul(num_rows))
        .unwrap_or_default();
    product(
        TheFirstOperand::ByTheRowsOfTheResult {
            values: kept,
            rows: num_kept,
        },
        num_rows,
        TheSecondOperand::ByTheValuesSummedOver {
            values: vectors,
            cols: side,
        },
        &mut rows,
    )
    .map_err(|source| Error::PcoaLinalg {
        operation: "product that gives the eigenvectors of the constant",
        source,
    })?;
    for (row, value) in rows.chunks_exact_mut(side).zip(&directions.values) {
        let length = value.sqrt();
        for entry in row.iter_mut() {
            *entry /= length;
        }
    }
    let of_the_ones = 1.0 / (side as f64).sqrt();
    let mut targets = vectors.chunks_exact_mut(side);
    for (target, row) in targets.by_ref().zip(rows.chunks_exact(side)) {
        target.copy_from_slice(row);
    }
    for target in targets {
        target.fill(of_the_ones);
    }
    Ok(())
}

/// The eigendecomposition of B, with how many of its eigenvalues are
/// negative and the part of the sum of all of them that those are.
#[derive(Debug)]
struct Decomposed {
    /// The eigenvalues of B from the largest and its eigenvectors, of the
    /// distances divided by the largest of them.
    eigen: Eigen,
    /// How many eigenvalues are above [`the_threshold_of_the_eigenvalues`],
    /// which are the components.
    num_positive: usize,
    /// How many eigenvalues are below minus that threshold.
    num_negative: usize,
    /// 100 times the sum of the absolute values of those eigenvalues over
    /// the sum of every eigenvalue, and 0 when there are none.
    negative_eigenvalues_percent: f64,
}

/// The eigendecomposition of the centered matrix `centered` of
/// `num_individuals` individuals, and what its negative eigenvalues are.
///
/// # Errors
///
/// [`Error::PcoaLinalg`] when the eigendecomposition could not be done.
fn the_decomposition_of(centered: Vec<f64>, num_individuals: usize) -> Result<Decomposed> {
    let eigen = eigh_lower(centered, num_individuals).map_err(|source| Error::PcoaLinalg {
        operation: "eigendecomposition",
        source,
    })?;
    let threshold = the_threshold_of_the_eigenvalues(&eigen.values, num_individuals);
    let num_positive = eigen
        .values
        .iter()
        .take_while(|value| **value > threshold)
        .count();
    let (num_negative, sum_of_the_negative) = eigen
        .values
        .iter()
        .filter(|value| **value < -threshold)
        .fold((0_usize, 0.0), |(count, sum), value| {
            (count.saturating_add(1), sum + value.abs())
        });
    let negative_eigenvalues_percent = if num_negative == 0 {
        0.0
    } else {
        let sum_of_every_eigenvalue: f64 = eigen.values.iter().sum();
        100.0 * (sum_of_the_negative / sum_of_every_eigenvalue)
    };
    Ok(Decomposed {
        eigen,
        num_positive,
        num_negative,
        negative_eigenvalues_percent,
    })
}

/// The threshold of the eigenvalues of B: one is positive above it,
/// negative below minus it, and 0 in between. It is the individuals times
/// 2.2e-16, the difference between 1 and the next number an `f64` holds,
/// times the sum of the absolute values of every eigenvalue.
///
/// It is not the threshold of the principal components, the largest
/// eigenvalue times the individuals times 2.2e-16, which is narrower than
/// the rounding of B: each cell of B carries the rounding of its centering,
/// and with that threshold Euclidean matrices of random points were refused
/// or given a component of rounding, and distances corrected by
/// [`correct_dists_by_lingoes`] were refused, as "What it gives" of the
/// principal coordinates in `docs/specs/pca.md` says. The sum of the
/// absolute values is the largest at least and the individuals times it at
/// most, so a component below the individuals squared times 2.2e-16 of the
/// largest, 2.2e-8 of it at 10000 individuals, is not given.
pub(crate) fn the_threshold_of_the_eigenvalues(values: &[f64], num_individuals: usize) -> f64 {
    let sum_of_their_sizes: f64 = values.iter().map(|value| value.abs()).sum();
    num_individuals as f64 * f64::EPSILON * sum_of_their_sizes
}

/// The components of the positive eigenvalues of a decomposition of B,
/// whose distances were divided by `largest`: the projections are
/// multiplied back by it, and the percentages are the same for both.
///
/// # Errors
///
/// [`Error::PcoaComponentAlongTheVectorOfOnes`] when a component is not at
/// a right angle to the vector of ones, and [`Error::PcoaNoMemory`] when
/// the machine does not give the memory of the projections.
fn the_components_of(
    decomposed: &Decomposed,
    num_individuals: usize,
    largest: f64,
) -> Result<Pcoa> {
    // A matrix with no negative eigenvalue has a sum of their sizes of at
    // most the individuals times the largest, so the threshold is at most
    // 46340 squared times 2.2e-16 of the largest, 0.48 of it, and the
    // largest is above it: the largest distance is above 0, so the largest
    // eigenvalue is too.
    let eigen = &decomposed.eigen;
    let num_comps = decomposed.num_positive;
    refuse_a_component_along_the_vector_of_ones(eigen, num_individuals, num_comps)?;
    let num_projections = num_individuals
        .checked_mul(num_comps)
        .ok_or(Error::PcoaNoMemory {
            num_individuals,
            what: "projections",
        })?;
    let mut projections = zeros_or_no_memory(num_projections, num_individuals, "projections")?;
    write_the_projections(eigen, num_individuals, num_comps, &mut projections);
    for projection in &mut projections {
        *projection *= largest;
    }
    for component in 0..num_comps {
        fix_the_sign_of(&mut projections, component, num_comps);
    }
    Ok(Pcoa {
        num_individuals,
        num_comps,
        projections,
        explained_variance_percent: the_percentages_of(&eigen.values, num_comps),
        lingoes_constant: 0.0,
        negative_eigenvalues_percent: 0.0,
    })
}

/// Refuses a component, of the first `num_comps` eigenvectors of `eigen`,
/// that is not at a right angle to the vector of ones: the centering of B
/// takes that vector out of every component, and one along it would put
/// every individual near one projection.
///
/// # Errors
///
/// [`Error::PcoaComponentAlongTheVectorOfOnes`] with the first such
/// component, when its eigenvector is further than
/// [`ALONG_THE_VECTOR_OF_ONES`] along it.
fn refuse_a_component_along_the_vector_of_ones(
    eigen: &Eigen,
    num_individuals: usize,
    num_comps: usize,
) -> Result<()> {
    let root = (num_individuals as f64).sqrt();
    for (component, vector) in eigen
        .vectors
        .chunks_exact(num_individuals)
        .take(num_comps)
        .enumerate()
    {
        let along = vector.iter().sum::<f64>() / root;
        if along.abs() > ALONG_THE_VECTOR_OF_ONES {
            return Err(Error::PcoaComponentAlongTheVectorOfOnes { component, along });
        }
    }
    Ok(())
}

/// The buffer of B of `num_individuals` individuals, `num_individuals` x
/// `num_individuals` zeros.
///
/// # Errors
///
/// [`Error::PcoaNoMemory`] when the machine does not give its memory, 8
/// bytes a cell.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the individuals are at most MAX_INDIVIDUALS_OF_THE_VARIANTS, 46340, checked before, so n x n is at most 2147395600, which a usize of 32 bits holds"
)]
fn the_matrix_of(num_individuals: usize) -> Result<Vec<f64>> {
    zeros_or_no_memory(num_individuals * num_individuals, num_individuals, "matrix")
}

/// `num_values` zeros, the matrix B or the projections of `num_individuals`
/// individuals, or the error of a machine that does not give their memory.
///
/// The memory is asked for with `try_reserve_exact`, which gives it back as
/// an error where `vec![0.0; n]` would end the process: B of 46340
/// individuals is 17 GB, and a Python session that asks for it on a
/// machine that has less would be ended with no message.
///
/// # Errors
///
/// [`Error::PcoaNoMemory`] with the individuals and `what`.
fn zeros_or_no_memory(
    num_values: usize,
    num_individuals: usize,
    what: &'static str,
) -> Result<Vec<f64>> {
    let mut values: Vec<f64> = Vec::new();
    values
        .try_reserve_exact(num_values)
        .map_err(|_| Error::PcoaNoMemory {
            num_individuals,
            what,
        })?;
    values.resize(num_values, 0.0);
    Ok(values)
}

/// Refuses a vector, or the individuals it is said to be of, that no
/// principal coordinate analysis can be done on: more individuals than the
/// linear algebra decomposes the matrix of, a vector that is not of their
/// pairs, and fewer than 2 individuals, in that order.
///
/// # Errors
///
/// [`Error::PcoaTooManyIndividuals`], [`Error::PcoaDistVectorOfAnotherSize`]
/// and [`Error::PcoaTooFewIndividuals`].
fn refuse_a_vector_that_cannot_be_analysed(
    dist_vector: &[f64],
    num_individuals: usize,
) -> Result<()> {
    // Before the length, so that the pairs are counted in a number that
    // holds them, and a caller hears of the size it cannot do before the
    // vector it may not be able to build.
    refuse_too_many_individuals(num_individuals)?;
    if dist_vector.len() != the_number_of_pairs(num_individuals) {
        return Err(Error::PcoaDistVectorOfAnotherSize {
            num_dists: dist_vector.len(),
            num_individuals,
        });
    }
    refuse_too_few_individuals(num_individuals)
}

/// Refuses more individuals than the linear algebra decomposes the
/// individuals x individuals matrix of.
///
/// # Errors
///
/// [`Error::PcoaTooManyIndividuals`].
fn refuse_too_many_individuals(num_individuals: usize) -> Result<()> {
    if num_individuals > MAX_INDIVIDUALS_OF_THE_VARIANTS {
        return Err(Error::PcoaTooManyIndividuals { num_individuals });
    }
    Ok(())
}

/// Refuses fewer than 2 individuals, which have no distance to place.
///
/// # Errors
///
/// [`Error::PcoaTooFewIndividuals`].
fn refuse_too_few_individuals(num_individuals: usize) -> Result<()> {
    if num_individuals < 2 {
        return Err(Error::PcoaTooFewIndividuals { num_individuals });
    }
    Ok(())
}

/// The distances of a vector a user gave, in its order, with `None` for the
/// NaN of a pair that has no distance.
fn the_distances_of_the_vector(dist_vector: &[f64]) -> impl Iterator<Item = Option<f64>> + '_ {
    dist_vector
        .iter()
        .map(|dist| if dist.is_nan() { None } else { Some(*dist) })
}

/// The largest distance of the pairs `dists` gives, once they are found to
/// be distances a principal coordinate analysis can be done on.
///
/// `dists` gives, each time it is called, the distance of every pair of
/// `num_individuals` individuals in the order of the distance vector,
/// `None` for a pair that has none: the vector a user gave, or the Kosman
/// distances worked out from the sums of a pass, of which no vector is
/// built. The individuals are 2 at least and at most
/// [`MAX_INDIVIDUALS_OF_THE_VARIANTS`], and the distances are their pairs.
///
/// # Errors
///
/// [`Error::PcoaPairsWithNoDistance`], [`Error::PcoaDistanceOutOfRange`]
/// and [`Error::PcoaAllDistancesZero`], in that order.
fn the_largest_distance<I: Iterator<Item = Option<f64>>>(
    dists: impl Fn() -> I,
    num_individuals: usize,
    from: PcoaInput,
) -> Result<f64> {
    refuse_the_pairs_with_no_distance(dists(), num_individuals, from)?;
    let mut largest: f64 = 0.0;
    for (dist, (first, second)) in dists().zip(the_pairs(num_individuals)) {
        // The pairs with no distance were refused above.
        let Some(dist) = dist else { continue };
        if dist < 0.0 || dist.is_infinite() {
            return Err(Error::PcoaDistanceOutOfRange {
                first,
                second,
                value: dist,
            });
        }
        largest = largest.max(dist);
    }
    if largest <= 0.0 {
        return Err(Error::PcoaAllDistancesZero);
    }
    Ok(largest)
}

/// Refuses the pairs that have no distance, `None` in `dists`, with how
/// many there are, the first of them and the individual in the most of
/// them.
///
/// # Errors
///
/// [`Error::PcoaPairsWithNoDistance`] when there is one such pair or more.
fn refuse_the_pairs_with_no_distance(
    dists: impl Iterator<Item = Option<f64>>,
    num_individuals: usize,
    from: PcoaInput,
) -> Result<()> {
    let mut num_pairs_with_no_distance = 0_usize;
    let mut first_pair = None;
    // Asked for at the first pair with no distance, which most analyses do
    // not have.
    let mut of_each_individual: Vec<usize> = Vec::new();
    for (dist, (first, second)) in dists.zip(the_pairs(num_individuals)) {
        if dist.is_some() {
            continue;
        }
        num_pairs_with_no_distance = num_pairs_with_no_distance.saturating_add(1);
        first_pair.get_or_insert((first, second));
        if of_each_individual.is_empty() {
            of_each_individual = vec![0_usize; num_individuals];
        }
        // Both are below the individuals, which `the_pairs` gives, and a
        // count is at most the pairs of an individual.
        for individual in [first, second] {
            if let Some(count) = of_each_individual.get_mut(individual) {
                *count = count.saturating_add(1);
            }
        }
    }
    if num_pairs_with_no_distance == 0 {
        return Ok(());
    }
    let num_pairs = the_number_of_pairs(num_individuals);
    // The first of the individuals with the largest count: a later one
    // replaces it only when its count is above.
    let (most_often, most_often_count) = of_each_individual.iter().enumerate().fold(
        (0, 0),
        |(most, most_count), (individual, count)| {
            if *count > most_count {
                (individual, *count)
            } else {
                (most, most_count)
            }
        },
    );
    let (first_of_the_first, second_of_the_first) = first_pair.unwrap_or((0, 0));
    Err(Error::PcoaPairsWithNoDistance {
        num_pairs_with_no_distance,
        num_pairs,
        first_of_the_first,
        second_of_the_first,
        most_often,
        most_often_count,
        from,
    })
}

/// How many pairs `num_individuals` individuals make, n(n - 1)/2, which is
/// the length of their distance vector.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the individuals are at most MAX_INDIVIDUALS_OF_THE_VARIANTS, 46340, checked before, so the product is at most 2147395600, which a usize of 32 bits holds"
)]
fn the_number_of_pairs(num_individuals: usize) -> usize {
    num_individuals * num_individuals.saturating_sub(1) / 2
}

/// The pairs of individuals in the order of the distance vector, (0, 1),
/// (0, 2), ..., (1, 2), ....
fn the_pairs(num_individuals: usize) -> impl Iterator<Item = (usize, usize)> {
    (0..num_individuals).flat_map(move |first| {
        (first..num_individuals)
            .skip(1)
            .map(move |second| (first, second))
    })
}

/// B of the distances `dists`, in the order of the distance vector,
/// divided by `largest`, the largest of them, centered twice, as its lower
/// half: `num_individuals` x `num_individuals`, row after row, of
/// which the entries of column j at most i of row i are written and the
/// others are 0.
///
/// The distances are divided by the largest before they are squared, so
/// that a distance whose square an `f64` would hold as an infinity, above
/// 1.3e154, or as 0, below 1.5e-162 when all of them are, still gives B.
/// The eigenvalues of B are then those of the distances as given over the
/// square of the largest, the percentages are the same, and the
/// projections are the largest times those of B.
///
/// `matrix` is the buffer B is written into, [`the_matrix_of`] the
/// individuals, which the caller asks the machine for when it chooses.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the individuals are at most MAX_INDIVIDUALS_OF_THE_VARIANTS, 46340, checked before, so first + 1 and (first + 1) x n + first are below n x n, which a usize of 32 bits holds"
)]
fn the_centered_matrix(
    mut matrix: Vec<f64>,
    dists: impl Iterator<Item = Option<f64>>,
    num_individuals: usize,
    largest: f64,
) -> Vec<f64> {
    let side = num_individuals;
    // The sum of the values of A of each row, in two parts: those at the
    // pairs where the individual is the first, which are the distances of
    // its segment of the vector, and those where it is the second.
    let mut as_the_second = vec![0.0; side];
    // A pair with no distance was refused before, and a NaN in its place
    // would make the linear algebra refuse B, not give a number.
    let mut dists = dists.map(|dist| dist.unwrap_or(f64::NAN));
    let as_the_first: Vec<f64> = (0..side)
        .map(|first| {
            // The pairs (first, second) for every second above first are
            // the next side - 1 - first distances of the vector, and their
            // cells of the lower half are those of column first below the
            // diagonal, one row apart.
            let column = matrix
                .iter_mut()
                .skip((first + 1) * side + first)
                .step_by(side);
            let mut sum = 0.0;
            for ((cell, dist), of_the_second) in column
                .zip(dists.by_ref().take(side - 1 - first))
                .zip(as_the_second.iter_mut().skip(first + 1))
            {
                let scaled = dist / largest;
                let value = -0.5 * scaled * scaled;
                *cell = value;
                sum += value;
                *of_the_second += value;
            }
            sum
        })
        .collect();
    let row_sums: Vec<f64> = as_the_first
        .iter()
        .zip(&as_the_second)
        .map(|(first, second)| first + second)
        .collect();
    center_the_lower_half(&mut matrix, side, &row_sums);
    // The second centering takes out of the means of the rows and of the
    // columns what the rounding of the first left in them, which the
    // eigenvalue 0 of the centering would otherwise carry.
    let row_sums = the_row_sums_of_the_lower_half(&matrix, side);
    center_the_lower_half(&mut matrix, side, &row_sums);
    matrix
}

/// Takes from each cell of the lower half of the symmetric `matrix`, `side`
/// x `side`, the mean of its row and that of its column, and adds the mean
/// of the matrix, from `row_sums`, the sum of each row of the whole matrix.
fn center_the_lower_half(matrix: &mut [f64], side: usize, row_sums: &[f64]) {
    let row_means: Vec<f64> = row_sums.iter().map(|sum| sum / side as f64).collect();
    let mean: f64 = row_means.iter().sum::<f64>() / side as f64;
    for (row, (row_mean, num_in_the_half)) in
        matrix.chunks_exact_mut(side).zip(row_means.iter().zip(1..))
    {
        for (cell, column_mean) in row.iter_mut().zip(&row_means).take(num_in_the_half) {
            *cell = *cell - row_mean - column_mean + mean;
        }
    }
}

/// The sum of each row of the symmetric `matrix`, `side` x `side`, of which
/// the lower half is read: a cell below the diagonal is of its row and, by
/// the symmetry, of the row of its column.
fn the_row_sums_of_the_lower_half(matrix: &[f64], side: usize) -> Vec<f64> {
    let mut row_sums = vec![0.0; side];
    for (row, index) in matrix.chunks_exact(side).zip(0..) {
        let (of_the_columns, of_this_row_on) = row_sums.split_at_mut(index);
        let mut of_this_row = 0.0;
        for (cell, of_the_column) in row.iter().zip(of_the_columns.iter_mut()) {
            of_this_row += cell;
            *of_the_column += cell;
        }
        // The index is below the side, so the row has its diagonal and the
        // sums have this row.
        if let (Some(diagonal), Some(sum)) = (row.get(index), of_this_row_on.first_mut()) {
            *sum += of_this_row + diagonal;
        }
    }
    row_sums
}

/// What a distance given to a principal coordinate analysis has to be, which
/// the message of [`Error::PcoaDistanceOutOfRange`] ends with. It is `pub`
/// so that each binding crate, which writes that message with the names of
/// the two individuals in the place of their positions, takes it from here
/// and does not keep a copy of its own.
pub const WHAT_A_DISTANCE_HAS_TO_BE: &str = "a principal coordinate analysis needs every distance finite and 0 or above; a negative F_ST or f_2 is of two populations the dataset cannot tell apart";

/// What the message of [`Error::PcoaPairsWithNoDistance`] tells the user
/// to do, which depends on where the distances came from, with the names
/// of Python. It is `pub` so that each binding crate, which writes that
/// message with the names of the individuals in the place of their
/// positions, takes it from here, the TypeScript one rewriting the names
/// of the arguments in camelCase, and does not keep a copy of its own.
pub fn the_remedy_of_the_pairs_with_no_distance(from: PcoaInput) -> &'static str {
    match from {
        PcoaInput::Distances => {
            "a principal coordinate analysis places every individual by its distance to every other, so each of those pairs has to be given a distance or one of its two individuals taken out of the distances"
        }
        PcoaInput::Variants => {
            "those pairs were called together at fewer variants than `min_num_snps`, or at none; take that individual out with `filter_individuals`, lower `min_num_snps`, or run the PCA of the variants, which gives every individual a projection"
        }
    }
}

/// "there is 1 individual" or "there are 0 individuals", for the message
/// of [`Error::PcoaTooFewIndividuals`].
pub(crate) fn the_individuals_there_are(num_individuals: usize) -> String {
    if num_individuals == 1 {
        "there is 1 individual".to_owned()
    } else {
        format!("there are {num_individuals} individuals")
    }
}

/// The verb of a count of pairs: "1 of the 10 pairs has", "2 of the 10
/// pairs have".
pub fn have_or_has(count: usize) -> &'static str {
    if count == 1 { "has" } else { "have" }
}

/// The correction the message of [`Error::PcoaNotEuclidean`] names, which
/// depends on where the distances came from.
pub(crate) fn the_correction_of(from: PcoaInput) -> &'static str {
    match from {
        PcoaInput::Distances => "`correct_dists_by_lingoes`",
        PcoaInput::Variants => "`correct_by_lingoes`",
    }
}

/// A percentage for a message, with three significant digits: 2.98, 0.0417.
/// Two decimals alone would write a small negative part as 0.00 percent,
/// which reads as none.
pub(crate) fn the_percent_shown(percent: f64) -> impl fmt::Display {
    PercentShown(percent)
}

/// A percentage written with three significant digits.
struct PercentShown(f64);

impl fmt::Display for PercentShown {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self(percent) = *self;
        // The digits after the point that give three significant ones: 2
        // for 1 to 10, 1 for 10 to 100, 3 for 0.1 to 1 and so on. The
        // logarithm of a positive finite number is finite, and one of 1e-300
        // asks for 302 digits, which `format` writes.
        let decimals = if percent > 0.0 && percent.is_finite() {
            let magnitude = percent.log10().floor();
            if magnitude >= 2.0 {
                0
            } else {
                // The magnitude is below 2 and at least -324, so 2 minus it
                // is a whole number from 1 to 326.
                #[expect(
                    clippy::cast_possible_truncation,
                    clippy::cast_sign_loss,
                    reason = "2 minus the floor of the logarithm of a positive f64 below 100 is a whole number from 1 to 326"
                )]
                let decimals = (2.0 - magnitude) as usize;
                decimals
            }
        } else {
            2
        };
        write!(formatter, "{percent:.decimals$}")
    }
}

#[cfg(test)]
mod tests {
    use std::path::{Path, PathBuf};

    use super::{
        LingoesCorrection, Pcoa, PcoaInput, PcoaOfVariants, VariantPcoaOptions,
        correct_dists_by_lingoes, pcoa, pcoa_of_variants,
    };
    use popnei_linalg::Eigen;

    use crate::block::{Block, BlockReader};
    use crate::error::{Error, Result};
    use crate::filters::FilteringStats;
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{ChromTable, MAX_INDIVIDUALS_OF_THE_VARIANTS, Needs};

    /// The tolerance of "How it is verified" of the principal coordinates
    /// in `docs/specs/pca.md`: R's numbers are written with 15 significant
    /// digits and compared within 1e-9, as in the principal components.
    const TOLERANCE: f64 = 1e-9;

    /// The ten distances of five individuals, `i1` to `i5`, of `test_pcoa`
    /// of pyNei, the worked example of the spec, in the order of the
    /// distance vector. They are not Euclidean: one eigenvalue of B is
    /// negative.
    const SMALL: [f64; 10] = [0.2, 0.3, 0.9, 0.9, 0.1, 0.8, 0.7, 0.7, 0.8, 0.2];

    /// The path of one of the files of `tests/reference/`.
    fn the_reference_path(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference")
            .join(name)
    }

    /// The lines of a reference file that hold something.
    fn the_lines_of(name: &str) -> Vec<String> {
        let path = the_reference_path(name);
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("{path}: {error}", path = path.display()));
        text.lines()
            .filter(|line| !line.trim().is_empty())
            .map(str::to_owned)
            .collect()
    }

    /// A number of R's files, which carry leading spaces.
    fn the_number(field: &str) -> f64 {
        field
            .trim()
            .parse()
            .unwrap_or_else(|error| panic!("`{field}`: {error}"))
    }

    /// The projections of one of R's files, `*.r.projections.tsv`: one row
    /// for each individual, without the header and the name of each row.
    fn the_projections_of_r(name: &str) -> Vec<Vec<f64>> {
        the_lines_of(&format!("pca/{name}"))
            .iter()
            .skip(1)
            .map(|line| line.split('\t').skip(1).map(the_number).collect())
            .collect()
    }

    /// The one number of each line of one of R's files, the percentages of
    /// `*.r.percent.tsv` or the constant and the percent of the negative
    /// eigenvalues of `*.r.constant.tsv`.
    fn the_column_of_r(name: &str) -> Vec<f64> {
        the_lines_of(&format!("pca/{name}"))
            .iter()
            .map(|line| the_number(line))
            .collect()
    }

    /// The Kosman distances of the 40 individuals of `four_alleles.vcf.gz`
    /// as R's `gd.kosman` gave them, in the order of the distance vector.
    fn four_alleles_kosman() -> Vec<f64> {
        the_lines_of("dists/four_alleles.gdkosman.tsv")
            .iter()
            .skip(1)
            .map(|line| {
                let dist = line.split('\t').next().unwrap_or_else(|| panic!("{line}"));
                the_number(dist)
            })
            .collect()
    }

    /// The result, or the test fails with the error.
    fn the_pcoa_of(dist_vector: Vec<f64>, num_individuals: usize) -> Pcoa {
        pcoa(dist_vector, num_individuals).unwrap_or_else(|error| panic!("{error}"))
    }

    /// Asserts the projections and the percentages of a result against
    /// R's, each projection within `tolerance`.
    fn assert_as_r(result: &Pcoa, projections: &[Vec<f64>], percent: &[f64], tolerance: f64) {
        assert_eq!(result.num_individuals, projections.len());
        assert_eq!(result.num_comps, percent.len());
        assert_eq!(
            Some(result.projections.len()),
            result.num_individuals.checked_mul(result.num_comps)
        );
        for (individual, (ours, of_r)) in result
            .projections
            .chunks_exact(result.num_comps)
            .zip(projections)
            .enumerate()
        {
            assert_eq!(ours.len(), of_r.len(), "the components of {individual}");
            for (component, (ours, of_r)) in ours.iter().zip(of_r).enumerate() {
                assert!(
                    (ours - of_r).abs() <= tolerance,
                    "individual {individual}, component {component}: {ours} against R's {of_r}"
                );
            }
        }
        for (component, (ours, of_r)) in result
            .explained_variance_percent
            .iter()
            .zip(percent)
            .enumerate()
        {
            assert!(
                (ours - of_r).abs() <= TOLERANCE,
                "the percentage of component {component}: {ours} against R's {of_r}"
            );
        }
    }

    /// The worked example is refused, with its one negative eigenvalue of
    /// five and the part of them it is, and the message names the
    /// correction a user of `do_pcoa` has.
    #[test]
    fn distances_that_are_not_euclidean_are_refused_with_their_negative_eigenvalues() {
        match pcoa(SMALL.to_vec(), 5) {
            Err(
                error @ Error::PcoaNotEuclidean {
                    num_negative,
                    num_individuals,
                    negative_eigenvalues_percent,
                    from,
                },
            ) => {
                assert_eq!(num_negative, 1);
                assert_eq!(num_individuals, 5);
                assert!(
                    (negative_eigenvalues_percent - 7.88262807403034).abs() <= TOLERANCE,
                    "{negative_eigenvalues_percent}"
                );
                assert_eq!(from, PcoaInput::Distances);
                let message = error.to_string();
                assert!(message.contains("`correct_dists_by_lingoes`"), "{message}");
                assert!(message.contains("1 of the 5 eigenvalues"), "{message}");
                assert!(message.contains("7.88 percent"), "{message}");
                assert!(!error.names_the_file(), "{message}");
            }
            other => panic!("the worked example was not refused: {other:?}"),
        }
    }

    /// The Kosman distances of `four_alleles.vcf.gz` are Euclidean, and
    /// their 39 components are R's: every eigenvalue but the 0 of the
    /// centering is positive.
    #[test]
    fn the_kosman_distances_of_four_alleles_give_the_components_of_r() {
        let result = the_pcoa_of(four_alleles_kosman(), 40);
        assert_eq!(result.num_comps, 39);
        assert_as_r(
            &result,
            &the_projections_of_r("four_alleles.pcoa.r.projections.tsv"),
            &the_column_of_r("four_alleles.pcoa.r.percent.tsv"),
            TOLERANCE,
        );
        assert!(
            (result.explained_variance_percent[0] - 4.33989428824864).abs() <= TOLERANCE,
            "{}",
            result.explained_variance_percent[0]
        );
    }

    /// The percentages of a Euclidean matrix add up to 100, and a result of
    /// `pcoa`, which corrects nothing, has no constant and no negative
    /// part.
    #[test]
    fn a_euclidean_matrix_has_percentages_that_add_up_to_100_and_no_correction() {
        let result = the_pcoa_of(four_alleles_kosman(), 40);
        let total: f64 = result.explained_variance_percent.iter().sum();
        assert!((total - 100.0).abs() <= 1e-12, "{total}");
        assert_eq!(result.lingoes_constant.to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            result.negative_eigenvalues_percent.to_bits(),
            0.0_f64.to_bits()
        );
    }

    /// The analysis gives the same components, times the same factor, for
    /// distances whose squares would be an infinity or would be 0 in an
    /// `f64`: B is built from the distances divided by the largest of them.
    #[test]
    fn distances_whose_squares_an_f64_cannot_hold_give_the_same_components_scaled() {
        let projections = the_projections_of_r("four_alleles.pcoa.r.projections.tsv");
        let percent = the_column_of_r("four_alleles.pcoa.r.percent.tsv");
        for factor in [1e200, 1e-200] {
            let scaled: Vec<f64> = four_alleles_kosman().iter().map(|d| d * factor).collect();
            let result = the_pcoa_of(scaled, 40);
            let expected: Vec<Vec<f64>> = projections
                .iter()
                .map(|row| row.iter().map(|value| value * factor).collect())
                .collect();
            assert_as_r(&result, &expected, &percent, TOLERANCE * factor);
        }
    }

    /// Two pairs with no distance, in which the individual 3 is in both:
    /// the error has the count, the first pair in the order of the
    /// distance vector, and the individual in the most of them.
    #[test]
    fn pairs_with_no_distance_are_refused_with_their_count_and_positions() {
        let mut dist_vector = SMALL.to_vec();
        // (0, 3) and (1, 3), at the positions 2 and 5 of the vector.
        dist_vector[2] = f64::NAN;
        dist_vector[5] = f64::NAN;
        match pcoa(dist_vector, 5) {
            Err(
                error @ Error::PcoaPairsWithNoDistance {
                    num_pairs_with_no_distance,
                    num_pairs,
                    first_of_the_first,
                    second_of_the_first,
                    most_often,
                    most_often_count,
                    from,
                },
            ) => {
                assert_eq!(num_pairs_with_no_distance, 2);
                assert_eq!(num_pairs, 10);
                assert_eq!((first_of_the_first, second_of_the_first), (0, 3));
                assert_eq!((most_often, most_often_count), (3, 2));
                assert_eq!(from, PcoaInput::Distances);
                let message = error.to_string();
                assert!(message.contains("2 of the 10 pairs"), "{message}");
                assert!(message.contains("given a distance"), "{message}");
                assert!(!message.contains("min_num_snps"), "{message}");
            }
            other => panic!("the pairs with no distance were not refused: {other:?}"),
        }
    }

    /// When two individuals are in as many pairs with no distance, the
    /// first of them in the order of the individuals is the one named,
    /// although the other reaches that count first along the vector.
    #[test]
    fn of_two_individuals_in_as_many_pairs_with_no_distance_the_first_is_named() {
        let mut dist_vector = SMALL.to_vec();
        // (1, 4), (2, 3) and (3, 4), at the positions 6, 7 and 9: 3 and 4
        // are in two of them each, 4 is met first along the vector, and 3
        // comes first among the individuals.
        dist_vector[6] = f64::NAN;
        dist_vector[7] = f64::NAN;
        dist_vector[9] = f64::NAN;
        match pcoa(dist_vector, 5) {
            Err(Error::PcoaPairsWithNoDistance {
                num_pairs_with_no_distance,
                first_of_the_first,
                second_of_the_first,
                most_often,
                most_often_count,
                ..
            }) => {
                assert_eq!(num_pairs_with_no_distance, 3);
                assert_eq!((first_of_the_first, second_of_the_first), (1, 4));
                assert_eq!((most_often, most_often_count), (3, 2));
            }
            other => panic!("the pairs with no distance were not refused: {other:?}"),
        }
    }

    /// A pair with no distance in distances that are not Euclidean is
    /// refused as a pair with no distance: it is looked for before the
    /// eigendecomposition.
    #[test]
    fn a_pair_with_no_distance_is_refused_before_the_matrix_is_found_not_euclidean() {
        let mut dist_vector = SMALL.to_vec();
        dist_vector[9] = f64::NAN;
        match pcoa(dist_vector, 5) {
            Err(Error::PcoaPairsWithNoDistance {
                first_of_the_first,
                second_of_the_first,
                ..
            }) => assert_eq!((first_of_the_first, second_of_the_first), (3, 4)),
            other => panic!("the pair with no distance was not refused: {other:?}"),
        }
    }

    /// A negative distance and an infinite one are refused with their pair
    /// and their value.
    #[test]
    fn a_negative_or_an_infinite_distance_is_refused_with_its_pair() {
        for (position, value, pair) in [(4, -0.05, (1, 2)), (8, f64::INFINITY, (2, 4))] {
            let mut dist_vector = SMALL.to_vec();
            dist_vector[position] = value;
            match pcoa(dist_vector, 5) {
                Err(
                    error @ Error::PcoaDistanceOutOfRange {
                        first,
                        second,
                        value: given,
                    },
                ) => {
                    assert_eq!((first, second), pair);
                    assert_eq!(given.to_bits(), value.to_bits());
                    let message = error.to_string();
                    assert!(message.contains("F_ST"), "{message}");
                }
                other => panic!("the distance {value} was not refused: {other:?}"),
            }
        }
    }

    /// Distances that are all 0 are refused, where B would be 0 and so
    /// would the threshold of its eigenvalues.
    #[test]
    fn distances_that_are_all_zero_are_refused() {
        match pcoa(vec![0.0; 6], 4) {
            Err(error @ Error::PcoaAllDistancesZero) => {
                assert!(
                    error.to_string().contains("nothing to do a PCoA with"),
                    "{error}"
                );
            }
            other => panic!("the distances that are all 0 were not refused: {other:?}"),
        }
    }

    /// No individual and one individual have no distance to be placed by.
    #[test]
    fn fewer_than_two_individuals_are_refused() {
        for num_individuals in [0, 1] {
            match pcoa(Vec::new(), num_individuals) {
                Err(Error::PcoaTooFewIndividuals {
                    num_individuals: said,
                }) => assert_eq!(said, num_individuals),
                other => panic!("{num_individuals} individuals were not refused: {other:?}"),
            }
        }
    }

    /// One individual more than the linear algebra decomposes the matrix
    /// of is refused before its vector is looked at, which would be of
    /// 1073720970 distances.
    #[test]
    fn more_individuals_than_the_linear_algebra_takes_are_refused() {
        let num_individuals = MAX_INDIVIDUALS_OF_THE_VARIANTS + 1;
        match pcoa(Vec::new(), num_individuals) {
            Err(Error::PcoaTooManyIndividuals {
                num_individuals: said,
            }) => assert_eq!(said, num_individuals),
            other => panic!("{num_individuals} individuals were not refused: {other:?}"),
        }
    }

    /// A vector that does not hold n(n - 1)/2 distances is refused with
    /// its length and the individuals, whether it is shorter or longer.
    #[test]
    fn a_vector_of_another_length_than_the_pairs_is_refused() {
        for length in [9, 11] {
            match pcoa(vec![0.5; length], 5) {
                Err(Error::PcoaDistVectorOfAnotherSize {
                    num_dists,
                    num_individuals,
                }) => assert_eq!((num_dists, num_individuals), (length, 5)),
                other => panic!("a vector of {length} was not refused: {other:?}"),
            }
        }
    }

    /// The distances of `small_twin` of the spec: `SMALL` with a sixth
    /// individual, `i6`, at distance 0 from `i5` and at the distances of
    /// `i5` from the others, whose pairs come last of each row of the
    /// vector.
    const TWIN: [f64; 15] = [
        0.2, 0.3, 0.9, 0.9, 0.9, //
        0.1, 0.8, 0.7, 0.7, //
        0.7, 0.8, 0.8, //
        0.2, 0.2, //
        0.0,
    ];

    /// The correction, or the test fails with the error.
    fn the_correction_of(dist_vector: Vec<f64>, num_individuals: usize) -> LingoesCorrection {
        correct_dists_by_lingoes(dist_vector, num_individuals)
            .unwrap_or_else(|error| panic!("{error}"))
    }

    /// The correction of the worked example has the constant of R, the
    /// negative part of the distances given and, for the first pair,
    /// sqrt(0.2² + 2c).
    #[test]
    fn the_correction_of_the_worked_example_has_the_constant_of_r() {
        let correction = the_correction_of(SMALL.to_vec(), 5);
        assert!(
            (correction.constant - 0.0640069399611263).abs() <= TOLERANCE,
            "{}",
            correction.constant
        );
        assert!(
            (correction.negative_eigenvalues_percent - 7.88262807403034).abs() <= TOLERANCE,
            "{}",
            correction.negative_eigenvalues_percent
        );
        assert_eq!(correction.dist_vector.len(), 10);
        assert!(
            (correction.dist_vector[0] - 0.409894962060102).abs() <= TOLERANCE,
            "{}",
            correction.dist_vector[0]
        );
        let of_r = the_column_of_r("small.lingoes.r.constant.tsv");
        assert!((correction.constant - of_r[0]).abs() <= TOLERANCE);
        assert!((correction.negative_eigenvalues_percent - of_r[1]).abs() <= TOLERANCE);
    }

    /// The principal coordinates of the corrected worked example are the
    /// table of the spec, 3 components, with no constant and no negative
    /// part, since `pcoa` corrected nothing.
    #[test]
    fn the_corrected_worked_example_gives_the_table_of_the_spec() {
        let correction = the_correction_of(SMALL.to_vec(), 5);
        let result = the_pcoa_of(correction.dist_vector, 5);
        assert_eq!(result.num_comps, 3);
        let table = [
            vec![-0.431869046368213, -0.0415086068851603, 0.224881556185394],
            vec![-0.283479006142767, -0.158680403878776, -0.138801060767122],
            vec![-0.269028184151739, 0.212468396790800, -0.131478395434635],
            vec![0.492920681079785, 0.195146288470815, 0.0612591236372794],
            vec![0.491455555582935, -0.207425674497679, -0.0158612236209160],
        ];
        let percent = [77.1278402980914, 14.3397713765961, 8.53238832531248];
        assert_as_r(&result, &table, &percent, TOLERANCE);
        assert_as_r(
            &result,
            &the_projections_of_r("small.lingoes.r.projections.tsv"),
            &the_column_of_r("small.lingoes.r.percent.tsv"),
            TOLERANCE,
        );
        assert_eq!(result.lingoes_constant.to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            result.negative_eigenvalues_percent.to_bits(),
            0.0_f64.to_bits()
        );
    }

    /// The twin, whose eigenvalue 0 has two eigenvectors, corrected and
    /// then analysed: the constant of R, the two twins sqrt(2c) apart, and
    /// 4 components, the last of which sets the twins at sqrt(2c)/2 either
    /// side of 0, the first of them positive.
    #[test]
    fn the_corrected_twin_gives_the_four_components_of_r() {
        let correction = the_correction_of(TWIN.to_vec(), 6);
        assert!(
            (correction.constant - 0.072805704185076).abs() <= TOLERANCE,
            "{}",
            correction.constant
        );
        let of_r = the_column_of_r("small_twin.lingoes.r.constant.tsv");
        assert!((correction.constant - of_r[0]).abs() <= TOLERANCE);
        assert!(
            (correction.negative_eigenvalues_percent - of_r[1]).abs() <= TOLERANCE,
            "{}",
            correction.negative_eigenvalues_percent
        );
        let twins = correction.dist_vector[14];
        assert!(
            (twins - 2.0 * 0.190795314650381).abs() <= TOLERANCE,
            "{twins}"
        );
        let result = the_pcoa_of(correction.dist_vector, 6);
        assert_eq!(result.num_comps, 4);
        assert_as_r(
            &result,
            &the_projections_of_r("small_twin.lingoes.r.projections.tsv"),
            &[
                74.4550063085174,
                12.9405131851367,
                7.29289082221065,
                5.31158968413514,
            ],
            TOLERANCE,
        );
        assert!((result.projections[4 * 4 + 3] - 0.190795314650381).abs() <= TOLERANCE);
        assert!((result.projections[5 * 4 + 3] + 0.190795314650381).abs() <= TOLERANCE);
    }

    /// A Euclidean matrix is not corrected: the constant and the negative
    /// part are 0, and the vector is the one given, bit for bit.
    #[test]
    fn a_euclidean_matrix_is_given_back_as_it_was() {
        let dist_vector = four_alleles_kosman();
        let correction = the_correction_of(dist_vector.clone(), 40);
        let of_r = the_column_of_r("four_alleles.lingoes.r.constant.tsv");
        assert_eq!(correction.constant.to_bits(), of_r[0].to_bits());
        assert_eq!(
            correction.negative_eigenvalues_percent.to_bits(),
            of_r[1].to_bits()
        );
        let as_bits = |values: &[f64]| {
            values
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        };
        assert_eq!(as_bits(&correction.dist_vector), as_bits(&dist_vector));
    }

    /// The correction refuses a pair with no distance with the words of a
    /// `Distances`, which name no `min_num_snps`.
    #[test]
    fn the_correction_refuses_a_pair_with_no_distance() {
        let mut dist_vector = SMALL.to_vec();
        dist_vector[1] = f64::NAN;
        match correct_dists_by_lingoes(dist_vector, 5) {
            Err(
                error @ Error::PcoaPairsWithNoDistance {
                    num_pairs_with_no_distance: 1,
                    first_of_the_first: 0,
                    second_of_the_first: 2,
                    from: PcoaInput::Distances,
                    ..
                },
            ) => {
                let message = error.to_string();
                assert!(message.contains("given a distance"), "{message}");
                assert!(!error.names_the_file(), "{message}");
            }
            other => panic!("the pair with no distance was not corrected: {other:?}"),
        }
    }

    /// The correction refuses what `pcoa` refuses before the
    /// eigendecomposition: a negative distance, distances that are all 0,
    /// one individual, more individuals than the linear algebra takes and
    /// a vector of another length.
    #[test]
    fn the_correction_refuses_what_the_analysis_refuses() {
        let mut negative = SMALL.to_vec();
        negative[3] = -1.0;
        assert!(matches!(
            correct_dists_by_lingoes(negative, 5),
            Err(Error::PcoaDistanceOutOfRange {
                first: 0,
                second: 4,
                ..
            })
        ));
        assert!(matches!(
            correct_dists_by_lingoes(vec![0.0; 10], 5),
            Err(Error::PcoaAllDistancesZero)
        ));
        assert!(matches!(
            correct_dists_by_lingoes(Vec::new(), 1),
            Err(Error::PcoaTooFewIndividuals { num_individuals: 1 })
        ));
        assert!(matches!(
            correct_dists_by_lingoes(Vec::new(), MAX_INDIVIDUALS_OF_THE_VARIANTS + 1),
            Err(Error::PcoaTooManyIndividuals { .. })
        ));
        assert!(matches!(
            correct_dists_by_lingoes(SMALL.to_vec(), 6),
            Err(Error::PcoaDistVectorOfAnotherSize {
                num_dists: 10,
                num_individuals: 6,
            })
        ));
    }

    /// Distances whose c is not a normal `f64` are refused with the largest
    /// of them: an infinity at 1e200 times the worked example, a subnormal
    /// at 1e-160 and 1e-161 times it, which keeps a few significant digits
    /// of c, and 0 at 1e-200 times it, which would say that nothing was
    /// corrected. At 1e-153 times it c is 5.2e-308, above the smallest
    /// normal `f64`, 2.2e-308, and is given.
    #[test]
    fn a_constant_beyond_an_f64_is_refused() {
        for factor in [1e200, 1e-160, 1e-161, 1e-200] {
            let scaled: Vec<f64> = SMALL.iter().map(|d| d * factor).collect();
            match correct_dists_by_lingoes(scaled, 5) {
                Err(error @ Error::PcoaLingoesConstantOutOfRange { largest }) => {
                    assert_eq!(largest.to_bits(), (0.9 * factor).to_bits());
                    assert!(!error.names_the_file(), "{error}");
                }
                other => panic!("the constant at {factor:e} was given: {other:?}"),
            }
        }
        let scaled: Vec<f64> = SMALL.iter().map(|d| d * 1e-153).collect();
        let constant = the_correction_of(scaled, 5).constant;
        let expected = 0.0640069399611263e-306;
        assert!(
            ((constant - expected) / expected).abs() <= TOLERANCE,
            "{constant:e}"
        );
    }

    /// A generator of numbers uniform in [0, 1), splitmix64, so that the
    /// random matrices of the tests are the same on every machine.
    struct Uniform(u64);

    impl Uniform {
        fn next(&mut self) -> f64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut bits = self.0;
            bits = (bits ^ (bits >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            bits = (bits ^ (bits >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            bits ^= bits >> 31;
            #[expect(
                clippy::cast_precision_loss,
                reason = "53 bits, which an f64 holds exactly"
            )]
            let value = (bits >> 11) as f64 / (1_u64 << 53) as f64;
            value
        }
    }

    /// The straight line distances of `num_individuals` points of
    /// `num_dims` coordinates uniform in [0, 1), which are Euclidean.
    fn the_distances_of_random_points(
        seed: u64,
        num_individuals: usize,
        num_dims: usize,
    ) -> Vec<f64> {
        let mut uniform = Uniform(seed);
        let points: Vec<Vec<f64>> = (0..num_individuals)
            .map(|_| (0..num_dims).map(|_| uniform.next()).collect())
            .collect();
        let mut dist_vector = Vec::new();
        for (first, of_the_first) in points.iter().enumerate() {
            for of_the_second in points.iter().skip(first).skip(1) {
                let square: f64 = of_the_first
                    .iter()
                    .zip(of_the_second)
                    .map(|(one, other)| (one - other) * (one - other))
                    .sum();
                dist_vector.push(square.sqrt());
            }
        }
        dist_vector
    }

    /// Euclidean matrices of 100 random points in 3000 dimensions have
    /// 99 components, the 0 of the centering being the one left out. With
    /// the threshold of the principal components, the largest eigenvalue
    /// times the individuals times 2.2e-16, the rounding of B took seeds 1
    /// and 5 for not Euclidean and gave seeds 0, 2, 4 and 6 a 100th
    /// component, on LAPACK and on faer.
    #[test]
    fn euclidean_matrices_of_random_points_are_taken_with_one_component_fewer_than_the_points() {
        for seed in 0..8 {
            match pcoa(the_distances_of_random_points(seed, 100, 3000), 100) {
                Ok(result) => assert_eq!(result.num_comps, 99, "seed {seed}"),
                Err(error) => panic!("seed {seed}: {error}"),
            }
        }
    }

    /// Distances corrected by Lingoes are taken by `pcoa`, with n - 2
    /// components: the 0 of the centering and the most negative eigenvalue,
    /// which the correction brings to 0, give none. 40 matrices of 100
    /// individuals at distances uniform from 0.5 to 1; with the threshold of
    /// the principal components seeds 13, 18 and 25 were refused on LAPACK
    /// and 13 on faer, with a message that named the correction.
    #[test]
    fn distances_corrected_by_lingoes_are_taken_by_the_analysis() {
        for seed in 0..40 {
            let mut uniform = Uniform(seed);
            let dist_vector: Vec<f64> = (0..4950).map(|_| 0.5 + 0.5 * uniform.next()).collect();
            let correction = the_correction_of(dist_vector, 100);
            assert!(correction.constant > 0.0, "seed {seed}");
            match pcoa(correction.dist_vector, 100) {
                Ok(result) => assert_eq!(result.num_comps, 98, "seed {seed}"),
                Err(error) => panic!("seed {seed}: {error}"),
            }
        }
    }

    /// The messages are written for the count they carry: one individual
    /// is "there is 1 individual", no individual "there are 0
    /// individuals", and one pair with no distance "has" none.
    #[test]
    fn a_message_of_one_is_written_in_the_singular() {
        let message = |error: Error| error.to_string();
        assert!(
            message(Error::PcoaTooFewIndividuals { num_individuals: 1 })
                .starts_with("there is 1 individual, and"),
            "{}",
            message(Error::PcoaTooFewIndividuals { num_individuals: 1 })
        );
        assert!(
            message(Error::PcoaTooFewIndividuals { num_individuals: 0 })
                .starts_with("there are 0 individuals, and")
        );
        let mut dist_vector = SMALL.to_vec();
        dist_vector[4] = f64::NAN;
        match pcoa(dist_vector, 5) {
            Err(error @ Error::PcoaPairsWithNoDistance { .. }) => {
                let text = error.to_string();
                assert!(
                    text.starts_with("1 of the 10 pairs of individuals has no distance"),
                    "{text}"
                );
                assert!(text.contains("is in 1 of them"), "{text}");
            }
            other => panic!("the pair with no distance was not refused: {other:?}"),
        }
    }

    /// The largest distance of a constant that is refused is written with
    /// an exponent, and not with the 200 digits of 9e199 in full.
    #[test]
    fn the_largest_distance_of_a_refused_constant_is_written_with_an_exponent() {
        let message = Error::PcoaLingoesConstantOutOfRange { largest: 9e199 }.to_string();
        assert!(
            message.starts_with("the largest distance is 9e199,"),
            "{message}"
        );
        let message = Error::PcoaLingoesConstantOutOfRange { largest: 9e-161 }.to_string();
        assert!(
            message.starts_with("the largest distance is 9e-161,"),
            "{message}"
        );
    }

    /// Memory that the machine does not give is an error with the
    /// individuals and what it was for, and not the end of the process that
    /// `vec![0.0; n]` would be: more values than an allocation can hold,
    /// which no machine gives, stand for a machine with too little.
    #[test]
    fn memory_the_machine_does_not_give_is_an_error() {
        match super::zeros_or_no_memory(usize::MAX / 4, 46340, "matrix") {
            Err(
                error @ Error::PcoaNoMemory {
                    num_individuals: 46340,
                    what: "matrix",
                },
            ) => assert!(
                error
                    .to_string()
                    .contains("calculate over fewer individuals"),
                "{error}"
            ),
            other => panic!("the memory was given: {other:?}"),
        }
    }

    /// A pair with no distance comes before a negative distance that comes
    /// before it along the vector: the pairs with no distance are looked
    /// for over the whole vector first, so a NaN at the position 0 and -1
    /// at the position 1 is the error of the pairs with no distance, which
    /// a check of each distance in turn would have made the other one.
    #[test]
    fn a_pair_with_no_distance_is_refused_before_a_negative_distance() {
        let mut dist_vector = SMALL.to_vec();
        dist_vector[0] = f64::NAN;
        dist_vector[1] = -1.0;
        for refused in [
            pcoa(dist_vector.clone(), 5).map(|_| ()),
            correct_dists_by_lingoes(dist_vector, 5).map(|_| ()),
        ] {
            match refused {
                Err(Error::PcoaPairsWithNoDistance {
                    first_of_the_first: 0,
                    second_of_the_first: 1,
                    ..
                }) => {}
                other => panic!("the pair with no distance was not refused first: {other:?}"),
            }
        }
    }

    /// The options of the analysis of the variants without the correction,
    /// with every pair that was called together at one variant at least
    /// given a distance.
    const NOT_CORRECTED: VariantPcoaOptions = VariantPcoaOptions {
        min_num_vars: 0,
        correct_by_lingoes: false,
    };

    /// A VCF reader over one of the diploid reference files of
    /// `tests/reference/`, `dists/panel.vcf.gz` for one, in blocks of
    /// `num_vars_per_block` variants.
    fn reader_of(
        name: &str,
        num_vars_per_block: usize,
    ) -> VcfReader<std::io::BufReader<std::fs::File>> {
        let options = VcfOptions {
            ploidy: 2,
            num_vars_per_block: Some(num_vars_per_block),
            ..VcfOptions::default()
        };
        VcfReader::from_path(&the_reference_path(name), options)
            .unwrap_or_else(|error| panic!("{name}: {error}"))
    }

    /// The principal coordinates of the variants of a reference file.
    fn the_pcoa_of_the_variants(
        name: &str,
        options: &VariantPcoaOptions,
    ) -> Result<PcoaOfVariants> {
        pcoa_of_variants(&mut reader_of(name, 100), options)
    }

    /// The panel, 200 individuals and 1200 variants whose Kosman distances
    /// R's `gd.kosman` gives, is not Euclidean: without the correction it is
    /// refused with 44 negative eigenvalues of 200 and 2.98343616373556
    /// percent, and the message names the option that corrects it. It is
    /// of the file that was read.
    #[test]
    fn the_panel_is_refused_without_the_correction() {
        match the_pcoa_of_the_variants("dists/panel.vcf.gz", &NOT_CORRECTED) {
            Err(
                error @ Error::PcoaNotEuclidean {
                    num_negative,
                    num_individuals,
                    negative_eigenvalues_percent,
                    from,
                },
            ) => {
                assert_eq!((num_negative, num_individuals), (44, 200));
                assert!(
                    (negative_eigenvalues_percent - 2.98343616373556).abs() <= TOLERANCE,
                    "{negative_eigenvalues_percent}"
                );
                assert_eq!(from, PcoaInput::Variants);
                let message = error.to_string();
                assert!(message.contains("`correct_by_lingoes`"), "{message}");
                assert!(!message.contains("correct_dists_by_lingoes"), "{message}");
                assert!(error.names_the_file(), "{message}");
            }
            other => panic!("the panel was not refused: {other:?}"),
        }
    }

    /// With a pair needing 1105 variants called in both, 35 of the 19900
    /// pairs of the panel have no distance, the first of them `s001` and
    /// `s082`, and `s082` is in 17 of them; the message tells a user of the
    /// variants what to change.
    #[test]
    fn the_panel_with_1105_variants_a_pair_has_pairs_with_no_distance() {
        let options = VariantPcoaOptions {
            min_num_vars: 1105,
            correct_by_lingoes: false,
        };
        match the_pcoa_of_the_variants("dists/panel.vcf.gz", &options) {
            Err(
                error @ Error::PcoaPairsWithNoDistance {
                    num_pairs_with_no_distance,
                    num_pairs,
                    first_of_the_first,
                    second_of_the_first,
                    most_often,
                    most_often_count,
                    from,
                },
            ) => {
                assert_eq!((num_pairs_with_no_distance, num_pairs), (35, 19900));
                assert_eq!((first_of_the_first, second_of_the_first), (1, 82));
                assert_eq!((most_often, most_often_count), (82, 17));
                assert_eq!(from, PcoaInput::Variants);
                let message = error.to_string();
                assert!(message.contains("`min_num_snps`"), "{message}");
                assert!(message.contains("`filter_individuals`"), "{message}");
                assert!(error.names_the_file(), "{message}");
            }
            other => panic!("the pairs with no distance were not refused: {other:?}"),
        }
    }

    /// The Kosman distances of `four_alleles.vcf.gz`, 40 individuals and
    /// 300 variants, are Euclidean: 39 components, R's, with no constant
    /// and no negative part, whatever the size of the blocks.
    #[test]
    fn the_variants_of_four_alleles_give_the_components_of_r() {
        for num_vars_per_block in [7, 100, 300] {
            let result = pcoa_of_variants(
                &mut reader_of("dists/four_alleles.vcf.gz", num_vars_per_block),
                &NOT_CORRECTED,
            )
            .unwrap_or_else(|error| panic!("{error}"));
            assert_eq!(result.num_vars, 300);
            assert_eq!(result.pcoa.num_comps, 39);
            assert_as_r(
                &result.pcoa,
                &the_projections_of_r("four_alleles.pcoa.r.projections.tsv"),
                &the_column_of_r("four_alleles.pcoa.r.percent.tsv"),
                TOLERANCE,
            );
            assert_eq!(result.pcoa.lingoes_constant.to_bits(), 0.0_f64.to_bits());
            assert_eq!(
                result.pcoa.negative_eigenvalues_percent.to_bits(),
                0.0_f64.to_bits()
            );
        }
    }

    /// A reader of one individual is refused before any of its blocks is
    /// read, and so is one of more individuals than the linear algebra
    /// takes.
    #[test]
    fn a_reader_of_too_few_or_too_many_individuals_is_refused_before_the_pass() {
        for (num_individuals, refused) in [
            (1, "fewer than 2"),
            (
                MAX_INDIVIDUALS_OF_THE_VARIANTS + 1,
                "more than the linear algebra takes",
            ),
        ] {
            let mut reader = NotToBeRead::of(num_individuals);
            let result = pcoa_of_variants(&mut reader, &NOT_CORRECTED);
            match (num_individuals, &result) {
                (1, Err(Error::PcoaTooFewIndividuals { num_individuals: 1 }))
                | (_, Err(Error::PcoaTooManyIndividuals { .. })) => {}
                _ => panic!("{refused}: {result:?}"),
            }
            assert_eq!(reader.num_blocks_asked_for, 0, "{refused}");
        }
    }

    /// A reader that counts the blocks it was asked for and gives none, for
    /// the refusals that come before the pass.
    struct NotToBeRead {
        individuals: Vec<String>,
        chroms: ChromTable,
        num_blocks_asked_for: usize,
    }

    impl NotToBeRead {
        fn of(num_individuals: usize) -> NotToBeRead {
            NotToBeRead {
                individuals: (0..num_individuals).map(|at| format!("i{at}")).collect(),
                chroms: ChromTable::new(),
                num_blocks_asked_for: 0,
            }
        }
    }

    impl BlockReader for NotToBeRead {
        fn next_block(&mut self) -> Result<Option<Block>> {
            self.num_blocks_asked_for = self.num_blocks_asked_for.saturating_add(1);
            Ok(None)
        }

        fn individuals(&self) -> &[String] {
            &self.individuals
        }

        fn ploidy(&self) -> usize {
            2
        }

        fn chroms(&self) -> &ChromTable {
            &self.chroms
        }

        fn set_needs(&mut self, _needs: Needs) {}

        fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
            Vec::new()
        }

        fn header(&self) -> &crate::block::SourceHeader {
            &crate::block::AN_EMPTY_HEADER
        }

        fn skip_outside(&mut self, _selection: crate::filters::RegionSelection) -> bool {
            false
        }

        fn num_skipped(&self) -> u64 {
            0
        }
    }

    /// The options of the analysis of the variants with the correction.
    const CORRECTED: VariantPcoaOptions = VariantPcoaOptions {
        min_num_vars: 0,
        correct_by_lingoes: true,
    };

    /// The panel corrected inside the analysis gives R's 198 components,
    /// the whole of `panel.lingoes.r.*.tsv`, with the constant and the
    /// negative part of the spec.
    #[test]
    fn the_panel_corrected_gives_the_components_of_r() {
        let result = the_pcoa_of_the_variants("dists/panel.vcf.gz", &CORRECTED)
            .unwrap_or_else(|error| panic!("{error}"));
        let pcoa = &result.pcoa;
        assert_eq!(result.num_vars, 1200);
        assert_eq!(pcoa.num_comps, 198);
        assert!(
            (pcoa.lingoes_constant - 0.014182298472042).abs() <= TOLERANCE,
            "{}",
            pcoa.lingoes_constant
        );
        assert!(
            (pcoa.negative_eigenvalues_percent - 2.98343616373556).abs() <= TOLERANCE,
            "{}",
            pcoa.negative_eigenvalues_percent
        );
        let of_r = the_column_of_r("panel.lingoes.r.constant.tsv");
        assert!((pcoa.lingoes_constant - of_r[0]).abs() <= TOLERANCE);
        assert!((pcoa.negative_eigenvalues_percent - of_r[1]).abs() <= TOLERANCE);
        for (individual, expected) in [
            (
                0,
                [0.0131009923566961, 0.103593034190948, -0.0461610570616341],
            ),
            (
                1,
                [0.0188069490905941, 0.102449570787996, -0.0506243303501242],
            ),
            (
                199,
                [-0.0728375733863349, -0.0163081366424616, 0.0108102147666687],
            ),
        ] {
            for (component, expected) in expected.iter().enumerate() {
                let ours = pcoa.projections[individual * 198 + component];
                assert!(
                    (ours - expected).abs() <= TOLERANCE,
                    "individual {individual}, component {component}: {ours}"
                );
            }
        }
        for (ours, expected) in pcoa.explained_variance_percent.iter().zip([
            9.62407114041929,
            6.65101405201371,
            1.73580899943624,
        ]) {
            assert!((ours - expected).abs() <= TOLERANCE, "{ours}");
        }
        assert_as_r(
            pcoa,
            &the_projections_of_r("panel.lingoes.r.projections.tsv"),
            &the_column_of_r("panel.lingoes.r.percent.tsv"),
            TOLERANCE,
        );
    }

    /// The panel with a clone, `s200`, whose genotypes are those of `s000`:
    /// its eigenvalue 0 has two eigenvectors, the vector of ones and the
    /// direction that sets the two clones apart, which the decomposition
    /// may give mixed. The correction gives the second the constant and
    /// leaves the first at 0, which is R's decomposition of the corrected
    /// matrix: 199 components, the whole of `panel_clone.lingoes.r.*.tsv`,
    /// and on `PC155` the two clones at sqrt(2c)/2 either side of 0.
    #[test]
    fn the_panel_with_a_clone_corrected_gives_the_components_of_r() {
        let result = the_pcoa_of_the_variants("pca/panel_clone.vcf.gz", &CORRECTED)
            .unwrap_or_else(|error| panic!("{error}"));
        let pcoa = &result.pcoa;
        assert_eq!(pcoa.num_comps, 199);
        assert!(
            (pcoa.lingoes_constant - 0.0143697121693818).abs() <= TOLERANCE,
            "{}",
            pcoa.lingoes_constant
        );
        assert!(
            (pcoa.negative_eigenvalues_percent - 2.97950453177523).abs() <= TOLERANCE,
            "{}",
            pcoa.negative_eigenvalues_percent
        );
        let first = pcoa.projections[155];
        let clone = pcoa.projections[200 * 199 + 155];
        assert!((first - 0.0847635303930356).abs() <= TOLERANCE, "{first}");
        assert!((clone + 0.0847635303930342).abs() <= TOLERANCE, "{clone}");
        assert_as_r(
            pcoa,
            &the_projections_of_r("panel_clone.lingoes.r.projections.tsv"),
            &the_column_of_r("panel_clone.lingoes.r.percent.tsv"),
            TOLERANCE,
        );
    }

    /// A Euclidean matrix asked to be corrected is not: the constant is 0,
    /// and the result is that of the analysis without the correction.
    #[test]
    fn the_variants_of_four_alleles_corrected_are_not_changed() {
        let result = the_pcoa_of_the_variants("dists/four_alleles.vcf.gz", &CORRECTED)
            .unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(result.pcoa.num_comps, 39);
        assert_eq!(result.pcoa.lingoes_constant.to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            result.pcoa.negative_eigenvalues_percent.to_bits(),
            0.0_f64.to_bits()
        );
        assert_as_r(
            &result.pcoa,
            &the_projections_of_r("four_alleles.pcoa.r.projections.tsv"),
            &the_column_of_r("four_alleles.pcoa.r.percent.tsv"),
            TOLERANCE,
        );
    }

    /// The distances of individuals that are copies of a few distinct
    /// ones: `copies[i]` is the distinct individual that the individual i
    /// is, and the distinct ones are at distances uniform from 0.5 to 1,
    /// drawn from `seed`. Two copies of one are at
    /// distance 0 from each other and at its distances from the others.
    #[expect(
        clippy::arithmetic_side_effects,
        clippy::needless_range_loop,
        reason = "the indices of a handful of individuals of a test, which fill a square matrix of their distances"
    )]
    fn the_distances_of_copies(seed: u64, copies: &[usize]) -> Vec<f64> {
        let num_distinct = copies.iter().max().map_or(0, |largest| largest + 1);
        let mut uniform = Uniform(seed);
        let mut of_the_distinct = vec![vec![0.0; num_distinct]; num_distinct];
        for first in 0..num_distinct {
            for second in first + 1..num_distinct {
                let dist = 0.5 + 0.5 * uniform.next();
                of_the_distinct[first][second] = dist;
                of_the_distinct[second][first] = dist;
            }
        }
        let mut dist_vector = Vec::new();
        for (first, of_the_first) in copies.iter().enumerate() {
            for of_the_second in copies.iter().skip(first + 1) {
                dist_vector.push(of_the_distinct[*of_the_first][*of_the_second]);
            }
        }
        dist_vector
    }

    /// The principal coordinates of a distance vector corrected by Lingoes
    /// inside the analysis, the route of the variants with
    /// `correct_by_lingoes`, from a vector: the band of the eigenvalue 0 is
    /// the one of the sums of a pass, whatever the distances came from.
    fn corrected_inside(dist_vector: &[f64], num_individuals: usize) -> Pcoa {
        let dists = || super::the_distances_of_the_vector(dist_vector);
        let largest = super::the_largest_distance(dists, num_individuals, PcoaInput::Variants)
            .unwrap_or_else(|error| panic!("{error}"));
        let matrix =
            super::the_matrix_of(num_individuals).unwrap_or_else(|error| panic!("{error}"));
        let centered = super::the_centered_matrix(matrix, dists(), num_individuals, largest);
        super::the_analysis_of(
            centered,
            num_individuals,
            largest,
            super::WhenNotEuclidean::Correct,
            PcoaInput::Variants,
        )
        .unwrap_or_else(|error| panic!("{error}"))
    }

    /// Bands of the eigenvalue 0 of 3 and of 6 eigenvectors: three copies of
    /// one individual among 7, two pairs of copies among 8, and 4 copies of
    /// one and 3 of another among 10. The corrected analysis inside gives
    /// what does not depend on the directions chosen inside the band: the
    /// squared distance of every pair rebuilt from the projections is d² +
    /// 2c, and the components, c and the percentages are those of
    /// `pcoa(correct_dists_by_lingoes(v))`, which decomposes the corrected
    /// matrix again.
    #[test]
    fn a_band_of_several_copies_is_corrected_as_the_corrected_distances_are() {
        for (what, copies) in [
            ("three copies of one", vec![0, 0, 0, 1, 2, 3, 4]),
            ("two pairs of copies", vec![0, 0, 1, 1, 2, 3, 4, 5]),
            (
                "4 copies of one and 3 of another",
                vec![0, 0, 0, 0, 1, 1, 1, 2, 3, 4],
            ),
        ] {
            let num_individuals = copies.len();
            // The first seed whose distinct individuals are not Euclidean,
            // which few of them often are.
            let (dist_vector, correction) = (0..100)
                .map(|seed| {
                    let dist_vector = the_distances_of_copies(seed, &copies);
                    let correction = the_correction_of(dist_vector.clone(), num_individuals);
                    (dist_vector, correction)
                })
                .find(|(_, correction)| correction.constant > 0.0)
                .unwrap_or_else(|| panic!("{what}: no seed below 100 is not Euclidean"));
            let inside = corrected_inside(&dist_vector, num_individuals);
            let of_the_corrected = the_pcoa_of(correction.dist_vector.clone(), num_individuals);
            assert_eq!(inside.num_comps, of_the_corrected.num_comps, "{what}");
            assert!(
                (inside.lingoes_constant - correction.constant).abs() <= TOLERANCE,
                "{what}: {} against {}",
                inside.lingoes_constant,
                correction.constant
            );
            for (ours, of_the_route) in inside
                .explained_variance_percent
                .iter()
                .zip(&of_the_corrected.explained_variance_percent)
            {
                assert!(
                    (ours - of_the_route).abs() <= TOLERANCE,
                    "{what}: {ours} against {of_the_route}"
                );
            }
            let num_comps = inside.num_comps;
            let rows: Vec<&[f64]> = inside.projections.chunks_exact(num_comps).collect();
            let mut at = 0;
            for first in 0..num_individuals {
                for second in first + 1..num_individuals {
                    let rebuilt: f64 = rows[first]
                        .iter()
                        .zip(rows[second])
                        .map(|(one, other)| (one - other) * (one - other))
                        .sum();
                    let expected =
                        dist_vector[at] * dist_vector[at] + 2.0 * inside.lingoes_constant;
                    assert!(
                        (rebuilt - expected).abs() <= TOLERANCE,
                        "{what}, the pair ({first}, {second}): {rebuilt} against {expected}"
                    );
                    at += 1;
                }
            }
        }
    }

    /// Four orthonormal vectors of 4 individuals, row after row: the
    /// vector of ones over 2, the direction that sets the individuals 2 and
    /// 3 apart, as two clones' is, and two more at right angles to both.
    fn four_orthonormal_vectors() -> [[f64; 4]; 4] {
        let half_root = std::f64::consts::FRAC_1_SQRT_2;
        [
            [0.5, 0.5, 0.5, 0.5],
            [0.0, 0.0, half_root, -half_root],
            [half_root, -half_root, 0.0, 0.0],
            [0.5, 0.5, -0.5, -0.5],
        ]
    }

    /// A decomposition of a B of 4 individuals whose rounding lifted the
    /// eigenvalue 0 of the vector of ones above the threshold, to 2e-15,
    /// while that of two clones stays at 0: the band of the eigenvalue 0
    /// holds the clones alone, and the correction refuses it as a defect
    /// instead of writing the vector of ones over the clones' eigenvector.
    #[test]
    fn a_band_without_the_vector_of_ones_is_a_defect() {
        let [ones, clones, first, last] = four_orthonormal_vectors();
        let eigen = Eigen {
            values: vec![1.0, 2e-15, 0.0, -0.5],
            vectors: [first, ones, clones, last].concat(),
        };
        match super::corrected_by_lingoes(eigen, 4) {
            Err(error @ Error::PcoaBandWithoutTheVectorOfOnes { .. }) => {
                assert!(error.to_string().contains("popnei has a defect"), "{error}");
                assert!(error.names_the_file(), "{error}");
            }
            other => panic!("the band without the vector of ones was taken: {other:?}"),
        }
    }

    /// A component whose eigenvector is the vector of ones, which would put
    /// every individual at one projection, is a defect and not a component.
    #[test]
    fn a_component_along_the_vector_of_ones_is_a_defect() {
        let [ones, clones, first, last] = four_orthonormal_vectors();
        let decomposed = super::Decomposed {
            eigen: Eigen {
                values: vec![1.0, 0.5, 0.0, 0.0],
                vectors: [first, ones, clones, last].concat(),
            },
            num_positive: 2,
            num_negative: 0,
            negative_eigenvalues_percent: 0.0,
        };
        match super::the_components_of(&decomposed, 4, 1.0) {
            Err(error @ Error::PcoaComponentAlongTheVectorOfOnes { component: 1, .. }) => {
                assert!(error.to_string().contains("popnei has a defect"), "{error}")
            }
            other => panic!("the component along the vector of ones was given: {other:?}"),
        }
    }

    /// The percentage of a message has three significant digits, so a
    /// small one is not written as 0.
    #[test]
    fn a_percentage_of_a_message_has_three_significant_digits() {
        let shown = |percent: f64| super::the_percent_shown(percent).to_string();
        assert_eq!(shown(7.88262807403034), "7.88");
        assert_eq!(shown(2.98343616373556), "2.98");
        assert_eq!(shown(41.66), "41.7");
        assert_eq!(shown(0.04166), "0.0417");
        assert_eq!(shown(100.0), "100");
    }
}
