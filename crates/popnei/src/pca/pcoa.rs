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
use crate::pca::{
    fix_the_sign_of, the_components_with_variance, the_percentages_of, the_projections_of,
    the_threshold_of_variance,
};
use crate::variant::MAX_INDIVIDUALS_OF_THE_VARIANTS;

/// What a principal coordinate analysis gives.
///
/// Only the components of the positive eigenvalues of B are here, and
/// `num_comps` is how many there are: an eigenvalue is positive when it is
/// above the largest times the individuals times 2.2e-16, the threshold of
/// the principal components, so the 0 that the centering of B always has
/// gives no component. In each component the projection of the largest
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
    Ok(the_components_of(
        &decomposed.eigen,
        num_individuals,
        largest,
    ))
}

/// The eigendecomposition of B, with how many of its eigenvalues are
/// negative and the part of the sum of all of them that those are.
struct Decomposed {
    /// The eigenvalues of B from the largest and its eigenvectors, of the
    /// distances divided by the largest of them.
    eigen: Eigen,
    /// How many eigenvalues are below minus the threshold of the
    /// components, which is the largest times the individuals times the
    /// epsilon of an `f64`.
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
    let largest = eigen.values.first().copied().unwrap_or(0.0);
    let threshold = the_threshold_of_variance(largest, num_individuals, num_individuals);
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
        num_negative,
        negative_eigenvalues_percent,
    })
}

/// The components of the positive eigenvalues of a decomposition of B,
/// whose distances were divided by `largest`: the projections are
/// multiplied back by it, and the percentages are the same for both.
fn the_components_of(eigen: &Eigen, num_individuals: usize, largest: f64) -> Pcoa {
    // The largest distance is above 0, so once the distances are divided by
    // it the sum of the diagonal of B, the sum of the squared distances over
    // the individuals, is at least 1 over the individuals, and the largest
    // eigenvalue at least that over the individuals again: the threshold, a
    // part of it of 1e-11 at most, leaves one component at least.
    let num_comps = the_components_with_variance(&eigen.values, num_individuals, num_individuals);
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

/// B of the distances divided by `largest`, the largest of them, as its
/// lower half: `num_individuals` x `num_individuals`, row after row, of
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
    let row_means: Vec<f64> = as_the_first
        .iter()
        .zip(&as_the_second)
        .map(|(first, second)| (first + second) / side as f64)
        .collect();
    let mean: f64 = row_means.iter().sum::<f64>() / side as f64;
    for (row, (row_mean, num_in_the_half)) in
        matrix.chunks_exact_mut(side).zip(row_means.iter().zip(1..))
    {
        for (cell, column_mean) in row.iter_mut().zip(&row_means).take(num_in_the_half) {
            *cell = *cell - row_mean - column_mean + mean;
        }
    }
    matrix
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

    use super::{Pcoa, PcoaInput, pcoa};
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
