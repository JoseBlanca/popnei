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

use popnei_linalg::{Eigen, eigh_lower};

use crate::error::{Error, Result};
use crate::pca::{fix_the_sign_of, the_percentages_of, the_projections_of};
use crate::variant::MAX_INDIVIDUALS_OF_THE_VARIANTS;

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
/// [`Error::PcoaLinalg`] when the eigendecomposition could not be done.
pub fn pcoa(dist_vector: Vec<f64>, num_individuals: usize) -> Result<Pcoa> {
    let largest = the_largest_distance(&dist_vector, num_individuals, PcoaInput::Distances)?;
    let centered = the_centered_matrix(&dist_vector, num_individuals, largest);
    drop(dist_vector);
    let decomposed = the_decomposition_of(centered, num_individuals)?;
    if decomposed.num_negative > 0 {
        return Err(Error::PcoaNotEuclidean {
            num_negative: decomposed.num_negative,
            num_individuals,
            negative_eigenvalues_percent: decomposed.negative_eigenvalues_percent,
            from: PcoaInput::Distances,
        });
    }
    Ok(the_components_of(&decomposed, num_individuals, largest))
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
/// when c, which is in the units of the squared distances, is beyond what
/// an `f64` holds.
pub fn correct_dists_by_lingoes(
    dist_vector: Vec<f64>,
    num_individuals: usize,
) -> Result<LingoesCorrection> {
    let largest = the_largest_distance(&dist_vector, num_individuals, PcoaInput::Distances)?;
    let centered = the_centered_matrix(&dist_vector, num_individuals, largest);
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
    // distances are not, so distances above 1.3e154, or all of them below
    // 1e-150, have a c beyond an f64, an infinity or 0, and a 0 would say
    // that nothing was corrected.
    let constant = constant_of_the_scaled * largest * largest;
    if !(constant.is_finite() && constant > 0.0) {
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

/// The eigendecomposition of B, with how many of its eigenvalues are
/// negative and the part of the sum of all of them that those are.
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
fn the_components_of(decomposed: &Decomposed, num_individuals: usize, largest: f64) -> Pcoa {
    // A matrix with no negative eigenvalue has a sum of their sizes of at
    // most the individuals times the largest, so the threshold is at most
    // 46340 squared times 2.2e-16 of the largest, 0.48 of it, and the
    // largest is above it: the largest distance is above 0, so the largest
    // eigenvalue is too.
    let eigen = &decomposed.eigen;
    let num_comps = decomposed.num_positive;
    let mut projections = the_projections_of(eigen, num_individuals, num_comps);
    for projection in &mut projections {
        *projection *= largest;
    }
    for component in 0..num_comps {
        fix_the_sign_of(&mut projections, component, num_comps);
    }
    Pcoa {
        num_individuals,
        num_comps,
        projections,
        explained_variance_percent: the_percentages_of(&eigen.values, num_comps),
        lingoes_constant: 0.0,
        negative_eigenvalues_percent: 0.0,
    }
}

/// The largest distance of the vector, once the vector and its
/// individuals are found to be ones a principal coordinate analysis can be
/// done on.
///
/// # Errors
///
/// [`Error::PcoaTooManyIndividuals`], [`Error::PcoaDistVectorOfAnotherSize`],
/// [`Error::PcoaTooFewIndividuals`], [`Error::PcoaPairsWithNoDistance`],
/// [`Error::PcoaDistanceOutOfRange`] and [`Error::PcoaAllDistancesZero`],
/// in that order, as [`pcoa`] says.
fn the_largest_distance(
    dist_vector: &[f64],
    num_individuals: usize,
    from: PcoaInput,
) -> Result<f64> {
    // Before the length, so that the pairs are counted in a number that
    // holds them, and a caller hears of the size it cannot do before the
    // vector it may not be able to build.
    if num_individuals > MAX_INDIVIDUALS_OF_THE_VARIANTS {
        return Err(Error::PcoaTooManyIndividuals { num_individuals });
    }
    let num_pairs = the_number_of_pairs(num_individuals);
    if dist_vector.len() != num_pairs {
        return Err(Error::PcoaDistVectorOfAnotherSize {
            num_dists: dist_vector.len(),
            num_individuals,
        });
    }
    if num_individuals < 2 {
        return Err(Error::PcoaTooFewIndividuals { num_individuals });
    }
    refuse_the_pairs_with_no_distance(dist_vector, num_individuals, num_pairs, from)?;
    let mut largest: f64 = 0.0;
    for (dist, (first, second)) in dist_vector.iter().zip(the_pairs(num_individuals)) {
        if *dist < 0.0 || dist.is_infinite() {
            return Err(Error::PcoaDistanceOutOfRange {
                first,
                second,
                value: *dist,
            });
        }
        largest = largest.max(*dist);
    }
    if largest <= 0.0 {
        return Err(Error::PcoaAllDistancesZero);
    }
    Ok(largest)
}

/// Refuses the pairs whose distance is a NaN, with how many there are, the
/// first of them and the individual in the most of them.
///
/// # Errors
///
/// [`Error::PcoaPairsWithNoDistance`] when there is one such pair or more.
fn refuse_the_pairs_with_no_distance(
    dist_vector: &[f64],
    num_individuals: usize,
    num_pairs: usize,
    from: PcoaInput,
) -> Result<()> {
    let num_pairs_with_no_distance = dist_vector.iter().filter(|dist| dist.is_nan()).count();
    if num_pairs_with_no_distance == 0 {
        return Ok(());
    }
    let mut first_pair = None;
    let mut of_each_individual = vec![0_usize; num_individuals];
    for (dist, (first, second)) in dist_vector.iter().zip(the_pairs(num_individuals)) {
        if !dist.is_nan() {
            continue;
        }
        first_pair.get_or_insert((first, second));
        // Both are below the individuals, which `the_pairs` gives, and a
        // count is at most the pairs of an individual.
        for individual in [first, second] {
            if let Some(count) = of_each_individual.get_mut(individual) {
                *count = count.saturating_add(1);
            }
        }
    }
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

/// B of the distances divided by `largest`, the largest of them, centered
/// twice, as its lower half: `num_individuals` x `num_individuals`, row after row, of
/// which the entries of column j at most i of row i are written and the
/// others are 0.
///
/// The distances are divided by the largest before they are squared, so
/// that a distance whose square an `f64` would hold as an infinity, above
/// 1.3e154, or as 0, below 1.5e-162 when all of them are, still gives B.
/// The eigenvalues of B are then those of the distances as given over the
/// square of the largest, the percentages are the same, and the
/// projections are the largest times those of B.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the individuals are at most MAX_INDIVIDUALS_OF_THE_VARIANTS, 46340, checked before, so n x n is at most 2147395600, which a usize of 32 bits holds, and first + 1 and (first + 1) x n + first are below it"
)]
fn the_centered_matrix(dist_vector: &[f64], num_individuals: usize, largest: f64) -> Vec<f64> {
    let side = num_individuals;
    let mut matrix = vec![0.0; side * side];
    // The sum of the values of A of each row, in two parts: those at the
    // pairs where the individual is the first, which are the distances of
    // its segment of the vector, and those where it is the second.
    let mut as_the_second = vec![0.0; side];
    let mut dists = dist_vector.iter();
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

/// What the message of [`Error::PcoaPairsWithNoDistance`] tells the user
/// to do, which depends on where the distances came from.
pub(crate) fn the_remedy_of_the_pairs_with_no_distance(from: PcoaInput) -> &'static str {
    match from {
        PcoaInput::Distances => {
            "a principal coordinate analysis places every individual by its distance to every other, so each of those pairs has to be given a distance or one of its two individuals taken out of the distances"
        }
        PcoaInput::Variants => {
            "those pairs were called together at fewer variants than `min_num_snps`, or at none; take that individual out with `filter_individuals`, lower `min_num_snps`, or run the PCA of the variants, which gives every individual a projection"
        }
    }
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

    use super::{LingoesCorrection, Pcoa, PcoaInput, correct_dists_by_lingoes, pcoa};
    use crate::error::Error;
    use crate::variant::MAX_INDIVIDUALS_OF_THE_VARIANTS;

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

    /// Distances whose c an `f64` does not hold, an infinity at 1e200 times
    /// the worked example and 0 at 1e-200 times it, are refused with the
    /// largest of them, where a c of 0 would say that nothing was
    /// corrected.
    #[test]
    fn a_constant_beyond_an_f64_is_refused() {
        for factor in [1e200, 1e-200] {
            let scaled: Vec<f64> = SMALL.iter().map(|d| d * factor).collect();
            match correct_dists_by_lingoes(scaled, 5) {
                Err(error @ Error::PcoaLingoesConstantOutOfRange { largest }) => {
                    assert_eq!(largest.to_bits(), (0.9 * factor).to_bits());
                    assert!(!error.names_the_file(), "{error}");
                }
                other => panic!("the constant at {factor} was given: {other:?}"),
            }
        }
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
