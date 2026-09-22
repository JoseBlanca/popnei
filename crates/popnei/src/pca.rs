//! The principal component analysis of a table of numbers.
//!
//! A principal component analysis places the rows of a table, the
//! individuals, on a few axes that hold as much of the variation between
//! them as that many axes can. Each column of the table, a trait, is
//! centered and divided by its standard deviation, which gives Z, and the
//! components are the directions in the space of the traits along which
//! the rows of Z vary most, the first having the largest variance that any
//! direction has, the second the largest among the directions at a right
//! angle to the first, and so on. [`Pca`] holds where each individual falls
//! along each component, how much of the variance each component holds and
//! the weight of each trait in each component.
//!
//! [`pca`] takes the table the user brings. The principal components of the
//! variants of a dataset are the other analysis of `docs/specs/pca.md` and
//! are not written yet.
//!
//! The products of matrices and the eigendecomposition are those of the
//! crate `popnei-linalg`, which runs them on the BLAS and LAPACK of the
//! system natively, and on faer in WebAssembly and natively when the cargo
//! feature `blas` is off.

use std::cmp::Ordering;
use std::fmt;

use popnei_linalg::{Eigen, add_self_product_lower, eigh_lower, product};

use crate::error::{Error, Result};

/// Whether the table is centered, which `do_pca` of pyNei does by default
/// and so does popnei.
pub const DEFAULT_CENTER_DATA: bool = true;

/// Whether the table is standardized, which `do_pca` of pyNei does by
/// default and so does popnei. A table that is standardized is centered as
/// well, so turning the centering off turns this off too.
pub const DEFAULT_STANDARDIZE_DATA: bool = true;

/// The two steps on the columns of a table before its components are
/// taken.
///
/// Centering takes the mean of each trait from it, which is what makes the
/// components the directions of the variation and not directions that
/// point at the mean of the data. Standardizing then divides each trait by
/// its standard deviation, which puts traits measured in different units
/// on one scale; without it the traits with the largest numbers dominate.
/// Standardizing without centering is an error, since the standard
/// deviation it divides by is the one the trait has once it is centered.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PcaOptions {
    /// Whether the mean of each trait is taken from it.
    pub center: bool,
    /// Whether each trait is divided by its standard deviation, the one
    /// with the number of rows in it and not the number of rows less one,
    /// which is pyNei's.
    pub standardize: bool,
}

/// What the mean or the standard deviation of a trait came out as when it
/// is not a number the analysis can use.
///
/// Each of the three is a trait whose values are too large or too small
/// for the arithmetic of an `f64`, and the user scales that trait or takes
/// it out of the table. They are found before anything is computed from
/// the mean or the deviation, so no analysis is done on the numbers they
/// would give.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraitScale {
    /// The values of the trait sum above the largest `f64`, 1.8e308, so
    /// its mean is an infinity and every centered value of it would be a
    /// NaN.
    MeanNotFinite,
    /// The squares of the deviations of the trait sum above the largest
    /// `f64`, which values of 1e154 give, so its standard deviation is an
    /// infinity and the standardized trait would be a column of zeros,
    /// which is what a trait with no variance gives.
    DeviationNotFinite,
    /// The squares of the deviations of the trait all fall below the
    /// smallest `f64` above 0, 5e-324, which values of 1e-200 give, so its
    /// standard deviation is 0 although its values are not all equal, and
    /// dividing by it would give infinities.
    DeviationOfZero,
}

impl fmt::Display for TraitScale {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let said = match *self {
            Self::MeanNotFinite => {
                "its values sum above the largest f64, so its mean is not finite"
            }
            Self::DeviationNotFinite => {
                "the squares of its deviations sum above the largest f64, so its standard deviation is not finite"
            }
            Self::DeviationOfZero => {
                "the squares of its deviations are all below the smallest f64 above 0, so its standard deviation is 0 although its values are not all equal"
            }
        };
        formatter.write_str(said)
    }
}

/// What a principal component analysis gives.
///
/// Only the components that have variance are here, and `num_comps` is how
/// many of them there are: centering takes one dimension out of the data,
/// so a table of 8 rows and 30 traits has 7 components and not 8. In each
/// component the projection of the largest absolute value is positive,
/// which is the rule that makes the result the same whichever library did
/// the eigendecomposition, and the weights of that component have the sign
/// that rule gave it.
#[derive(Debug, Clone)]
pub struct Pca {
    /// How many columns the data had: the variants the first pass gave,
    /// used or not, which is `num_vars` of the pass stats, or the traits.
    pub num_cols: usize,
    /// How many rows the data had, which is the individuals.
    pub num_rows: usize,
    /// How many components have variance, which is how many are given.
    pub num_comps: usize,
    /// num_rows x num_comps, row after row.
    pub projections: Vec<f64>,
    /// One per component, over the variance of every component of the data.
    pub explained_variance_percent: Vec<f64>,
    /// The positions of the columns that were used: the variants with
    /// variance among those the reader gave, or every trait of a table.
    pub used_cols: Vec<usize>,
    /// How many components the weights are given for, which for a table is
    /// `num_comps`.
    pub num_prin_comps: usize,
    /// num_prin_comps x used_cols.len(), row after row.
    pub princomps: Vec<f64>,
}

/// The principal components of the table `data`, which is `num_rows` rows
/// of `num_cols` values each, row after row, the rows being the
/// individuals and the columns the traits.
///
/// Every trait is a column of the result's weights, the one with no
/// variance included, which gets a weight of 0 when the table is not
/// standardized. No value may be missing: a table comes whole.
///
/// # Errors
///
/// [`Error::PcaValueNotFinite`] when a value of the table is an infinity
/// or a NaN. [`Error::PcaStandardizeWithoutCentering`] when the options
/// ask for the second and not the first. [`Error::PcaTableTooSmall`] when
/// the table has fewer than 2 rows or no traits.
/// [`Error::PcaTraitsWithNoVariance`] when the table is standardized and a
/// trait has no variance, and [`Error::PcaNoTraitWithVariance`] when no
/// trait of it has any. [`Error::PcaTraitOutOfRange`] when the mean or the
/// standard deviation of a trait is not a number the analysis can use.
/// [`Error::PcaTableOfAnotherSize`] when `data` does not hold exactly
/// `num_rows` times `num_cols` values. [`Error::PcaLinalg`] when the
/// product or the eigendecomposition could not be done.
pub fn pca(data: &[f64], num_rows: usize, num_cols: usize, options: &PcaOptions) -> Result<Pca> {
    if options.standardize && !options.center {
        return Err(Error::PcaStandardizeWithoutCentering);
    }
    if num_rows < 2 || num_cols == 0 {
        return Err(Error::PcaTableTooSmall { num_rows, num_cols });
    }
    // The buffer holds the table and nothing more: a longer one would be
    // analysed on its first values, which is not the table its caller
    // meant.
    let table = match num_rows.checked_mul(num_cols) {
        Some(num_values) if num_values == data.len() => data,
        _ => {
            return Err(Error::PcaTableOfAnotherSize {
                num_values: data.len(),
                num_rows,
                num_cols,
            });
        }
    };
    refuse_a_value_that_is_not_finite(table, num_cols)?;
    let (means, deviations) =
        the_center_and_the_scale_of_each_trait(table, num_rows, num_cols, options)?;
    // The matrix that is decomposed is the product of the smaller of the
    // two sides of the table with itself, Z Z' when there are fewer rows
    // than traits and Z' Z otherwise, and `add_self_product_lower` gives
    // the product of a matrix with itself over its columns. So the copy
    // that is centered and standardized is written with the traits as its
    // rows in the first case and as its columns in the second.
    let layout = if num_rows <= num_cols {
        Layout::TraitsAsRows
    } else {
        Layout::TraitsAsColumns
    };
    let standardized =
        the_standardized_table(table, num_rows, num_cols, &means, &deviations, layout);
    let (num_summed, side) = match layout {
        Layout::TraitsAsRows => (num_cols, num_rows),
        Layout::TraitsAsColumns => (num_rows, num_cols),
    };
    let mut gram = vec![0.0; num_values_of(side, side)];
    add_self_product_lower(&standardized, num_summed, side, &mut gram).map_err(|source| {
        Error::PcaLinalg {
            operation: "product of the table with itself",
            source,
        }
    })?;
    let eigen = eigh_lower(gram, side).map_err(|source| Error::PcaLinalg {
        operation: "eigendecomposition",
        source,
    })?;
    let num_comps = the_components_with_variance(&eigen.values, num_rows, num_cols);
    if num_comps == 0 {
        return Err(Error::PcaNoTraitWithVariance);
    }
    let variance_of_every_component: f64 = eigen.values.iter().sum();
    let explained_variance_percent = eigen
        .values
        .iter()
        .take(num_comps)
        // The share of the total before the 100, so that an eigenvalue
        // above 1.8e306, which a table of values of 1e153 gives, does not
        // become an infinity on the way to a number between 0 and 100.
        .map(|value| 100.0 * (value / variance_of_every_component))
        .collect();
    let (mut projections, mut princomps) = match layout {
        Layout::TraitsAsRows => the_components_of_the_product_of_the_rows(
            &standardized,
            &eigen,
            num_rows,
            num_cols,
            num_comps,
        )?,
        Layout::TraitsAsColumns => the_components_of_the_product_of_the_traits(
            &standardized,
            &eigen,
            num_rows,
            num_cols,
            num_comps,
        )?,
    };
    fix_the_signs(&mut projections, &mut princomps, num_comps, num_cols);
    Ok(Pca {
        num_cols,
        num_rows,
        num_comps,
        projections,
        explained_variance_percent,
        used_cols: (0..num_cols).collect(),
        num_prin_comps: num_comps,
        princomps,
    })
}

/// Where the traits of the table are in the copy of it that is centered
/// and standardized.
///
/// `add_self_product_lower` gives the product of a matrix with itself over
/// its columns, so the side of the table that is to be the side of that
/// product has to be the columns of the copy. The copy is written in the
/// layout that puts the smaller side there, and nothing is transposed
/// afterwards.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Layout {
    /// The traits are the rows of the copy, which is the transpose of the
    /// table. The product is then the rows x rows Z Z', which is what a
    /// table with more traits than rows takes.
    TraitsAsRows,
    /// The traits are the columns of the copy, which is the table as it
    /// came. The product is the traits x traits Z' Z.
    TraitsAsColumns,
}

/// The values of a matrix of `rows` x `cols`.
///
/// Every matrix of this module has one of its two sides the smaller side
/// of the table, or a count of components, which is at most that side, and
/// the other side at most the other side of the table. So every product
/// here is at most the rows of the table times its traits, which [`pca`]
/// took with `checked_mul` before it built any of them.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "one side is the smaller side of the table or a count of components, at most that side, and the other is at most the other side, so the product is at most the values of the table"
)]
fn num_values_of(rows: usize, cols: usize) -> usize {
    rows * cols
}

/// Refuses the first value of the table that is an infinity or a NaN, with
/// the row and the trait it is at.
///
/// # Errors
///
/// [`Error::PcaValueNotFinite`] with the place of that value.
fn refuse_a_value_that_is_not_finite(table: &[f64], num_cols: usize) -> Result<()> {
    for (row, values) in table.chunks_exact(num_cols).enumerate() {
        for (col, value) in values.iter().enumerate() {
            if !value.is_finite() {
                return Err(Error::PcaValueNotFinite {
                    row,
                    col,
                    value: *value,
                });
            }
        }
    }
    Ok(())
}

/// The mean of each trait and the standard deviation each one is divided
/// by, in the order of the traits.
///
/// A trait that is not centered has a mean of 0 and one that is not
/// standardized a standard deviation of 1, so that one subtraction and one
/// division give the value of every option. The divisor of the standard
/// deviation is the number of rows and not the number of rows less one,
/// which is pyNei's `data.std(axis=0)`.
///
/// # Errors
///
/// [`Error::PcaTraitsWithNoVariance`] when the table is standardized and
/// the values of a trait are all equal. [`Error::PcaTraitOutOfRange`] when
/// the mean or the standard deviation of a trait is not a number the
/// analysis can use, which [`TraitScale`] lists.
fn the_center_and_the_scale_of_each_trait(
    table: &[f64],
    num_rows: usize,
    num_cols: usize,
    options: &PcaOptions,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let rows = num_rows as f64;
    let mut means = vec![0.0; num_cols];
    if options.center {
        for values in table.chunks_exact(num_cols) {
            for (total, value) in means.iter_mut().zip(values) {
                *total += value;
            }
        }
        for mean in &mut means {
            *mean /= rows;
        }
        for (position, mean) in means.iter().enumerate() {
            if !mean.is_finite() {
                return Err(Error::PcaTraitOutOfRange {
                    position,
                    problem: TraitScale::MeanNotFinite,
                });
            }
        }
    }
    let mut deviations = vec![1.0; num_cols];
    if options.standardize {
        let positions = the_traits_with_no_variance(table, num_cols);
        if !positions.is_empty() {
            return Err(Error::PcaTraitsWithNoVariance {
                positions,
                num_cols,
            });
        }
        let mut squares = vec![0.0; num_cols];
        for values in table.chunks_exact(num_cols) {
            for ((total, value), mean) in squares.iter_mut().zip(values).zip(&means) {
                let deviation = value - mean;
                *total += deviation * deviation;
            }
        }
        for (deviation, total) in deviations.iter_mut().zip(&squares) {
            *deviation = (total / rows).sqrt();
        }
        for (position, deviation) in deviations.iter().enumerate() {
            if !deviation.is_finite() {
                return Err(Error::PcaTraitOutOfRange {
                    position,
                    problem: TraitScale::DeviationNotFinite,
                });
            }
            // The traits whose values are all equal were refused above, so
            // a deviation of 0 here is one whose squares were all too
            // small for an `f64` to hold.
            if *deviation == 0.0 {
                return Err(Error::PcaTraitOutOfRange {
                    position,
                    problem: TraitScale::DeviationOfZero,
                });
            }
        }
    }
    Ok((means, deviations))
}

/// The positions of the traits whose values are all equal, in order.
///
/// The values are compared as they are, which is what
/// `docs/specs/pca.md` asks for: pyNei tests `std == 0` on floats, which a
/// trait of 0.1 repeated can pass with a standard deviation of 1e-17, and
/// the division then gives numbers with no meaning.
#[expect(
    clippy::float_cmp,
    reason = "a trait has no variance when its values are equal as they are, and nothing is computed from them first"
)]
fn the_traits_with_no_variance(table: &[f64], num_cols: usize) -> Vec<usize> {
    let mut rows = table.chunks_exact(num_cols);
    // A table of no row has no trait with values that differ, and the
    // caller has refused one of fewer than 2 rows already.
    let Some(first) = rows.next() else {
        return Vec::new();
    };
    let mut equal_to_the_first = vec![true; num_cols];
    for values in rows {
        for ((same, value), reference) in equal_to_the_first.iter_mut().zip(values).zip(first) {
            if *value != *reference {
                *same = false;
            }
        }
    }
    equal_to_the_first
        .iter()
        .enumerate()
        .filter(|(_, same)| **same)
        .map(|(position, _)| position)
        .collect()
}

/// The table centered and divided by the deviations, in the layout the
/// product of the smaller side needs.
fn the_standardized_table(
    table: &[f64],
    num_rows: usize,
    num_cols: usize,
    means: &[f64],
    deviations: &[f64],
    layout: Layout,
) -> Vec<f64> {
    let mut standardized = vec![0.0; num_values_of(num_rows, num_cols)];
    match layout {
        Layout::TraitsAsRows => {
            for (position, ((values, mean), deviation)) in standardized
                .chunks_exact_mut(num_rows)
                .zip(means)
                .zip(deviations)
                .enumerate()
            {
                for (target, value) in values
                    .iter_mut()
                    .zip(table.iter().skip(position).step_by(num_cols))
                {
                    *target = (value - mean) / deviation;
                }
            }
        }
        Layout::TraitsAsColumns => {
            for (values, row) in standardized
                .chunks_exact_mut(num_cols)
                .zip(table.chunks_exact(num_cols))
            {
                for (((target, value), mean), deviation) in
                    values.iter_mut().zip(row).zip(means).zip(deviations)
                {
                    *target = (value - mean) / deviation;
                }
            }
        }
    }
    standardized
}

/// How many of the eigenvalues, which come from the largest, belong to a
/// component that has variance.
///
/// The threshold is the largest eigenvalue times the larger side of the
/// table times 2.220446049250313e-16, the difference between 1 and the
/// next number an `f64` holds, which is the tolerance numpy's
/// `matrix_rank` has for singular values, used here on eigenvalues.
/// `docs/specs/pca.md` has what was measured with it: the eigenvalue of a
/// component with no variance came out between -2e-16 and 3e-16 times the
/// largest on four tables, and 1.3e-14 times it at 1000 x 20000, where the
/// threshold is 4.4e-12 times it.
///
/// The side and the epsilon are multiplied first, so that a largest
/// eigenvalue near the largest `f64`, which a table of values of 2.5e153
/// gives, does not become an infinity on the way to a threshold that is a
/// small part of it.
fn the_components_with_variance(values: &[f64], num_rows: usize, num_cols: usize) -> usize {
    let Some(largest) = values.first() else {
        return 0;
    };
    let threshold = largest * (num_rows.max(num_cols) as f64 * f64::EPSILON);
    values
        .iter()
        .take_while(|value| **value > threshold)
        .count()
}

/// The projections and the weights when the matrix that was decomposed is
/// the product of the traits, Z' Z, whose eigenvectors are the weights
/// themselves and whose projections are Z times them.
///
/// `standardized` is the rows x traits matrix, and the eigenvectors are
/// over the traits.
///
/// # Errors
///
/// [`Error::PcaLinalg`] when the product that gives the projections could
/// not be done.
fn the_components_of_the_product_of_the_traits(
    standardized: &[f64],
    eigen: &Eigen,
    num_rows: usize,
    num_cols: usize,
    num_comps: usize,
) -> Result<(Vec<f64>, Vec<f64>)> {
    // [`pca`] refuses a table with no component before it gets here, and
    // this keeps `step_by` below off a step of 0, which panics.
    if num_comps == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let princomps: Vec<f64> = eigen
        .vectors
        .chunks_exact(num_cols)
        .take(num_comps)
        .flatten()
        .copied()
        .collect();
    let mut weights_by_trait = vec![0.0; num_values_of(num_cols, num_comps)];
    for (component, weights) in princomps.chunks_exact(num_cols).enumerate() {
        for (target, weight) in weights_by_trait
            .iter_mut()
            .skip(component)
            .step_by(num_comps)
            .zip(weights)
        {
            *target = *weight;
        }
    }
    let mut projections = vec![0.0; num_values_of(num_rows, num_comps)];
    product(
        standardized,
        num_rows,
        num_cols,
        &weights_by_trait,
        num_comps,
        &mut projections,
    )
    .map_err(|source| Error::PcaLinalg {
        operation: "product that gives the projections",
        source,
    })?;
    Ok((projections, princomps))
}

/// The projections and the weights when the matrix that was decomposed is
/// the product of the rows, Z Z', whose eigenvector u of the eigenvalue λ
/// gives the projections u sqrt(λ) and the weights Z' u / sqrt(λ).
///
/// `standardized` is the traits x rows matrix, which is Z', and the
/// eigenvectors are over the rows.
///
/// # Errors
///
/// [`Error::PcaLinalg`] when the product that gives the weights could not
/// be done.
fn the_components_of_the_product_of_the_rows(
    standardized: &[f64],
    eigen: &Eigen,
    num_rows: usize,
    num_cols: usize,
    num_comps: usize,
) -> Result<(Vec<f64>, Vec<f64>)> {
    // [`pca`] refuses a table with no component before it gets here, and
    // this keeps `step_by` below off a step of 0, which panics.
    if num_comps == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let mut projections = vec![0.0; num_values_of(num_rows, num_comps)];
    let mut vectors_by_row = vec![0.0; num_values_of(num_rows, num_comps)];
    for (component, (vector, value)) in eigen
        .vectors
        .chunks_exact(num_rows)
        .zip(&eigen.values)
        .take(num_comps)
        .enumerate()
    {
        let size = value.sqrt();
        for ((projection, scaled), coordinate) in projections
            .iter_mut()
            .skip(component)
            .step_by(num_comps)
            .zip(vectors_by_row.iter_mut().skip(component).step_by(num_comps))
            .zip(vector)
        {
            *projection = coordinate * size;
            *scaled = coordinate / size;
        }
    }
    let mut weights_by_trait = vec![0.0; num_values_of(num_cols, num_comps)];
    product(
        standardized,
        num_cols,
        num_rows,
        &vectors_by_row,
        num_comps,
        &mut weights_by_trait,
    )
    .map_err(|source| Error::PcaLinalg {
        operation: "product that gives the weights",
        source,
    })?;
    let mut princomps = vec![0.0; num_values_of(num_comps, num_cols)];
    for (component, weights) in princomps.chunks_exact_mut(num_cols).enumerate() {
        for (target, weight) in weights
            .iter_mut()
            .zip(weights_by_trait.iter().skip(component).step_by(num_comps))
        {
            *target = *weight;
        }
    }
    Ok((projections, princomps))
}

/// Gives each component the sign of the rule of `docs/specs/pca.md`: the
/// projection of the largest absolute value is positive, and when two
/// individuals have the same absolute value it is the first of them that
/// is made positive. The weights of a component that is turned round are
/// turned round with it.
///
/// `projections` is the individuals x `num_comps` matrix and `princomps`
/// the `num_comps` x `num_cols` one, both row after row.
fn fix_the_signs(
    projections: &mut [f64],
    princomps: &mut [f64],
    num_comps: usize,
    num_cols: usize,
) {
    for component in 0..num_comps {
        // The projection of the largest absolute value, and the first of
        // them when two are the same, which the strict comparison keeps.
        let largest = projections
            .iter()
            .skip(component)
            .step_by(num_comps)
            .copied()
            .fold(0.0_f64, |largest: f64, value| {
                if value.abs().total_cmp(&largest.abs()) == Ordering::Greater {
                    value
                } else {
                    largest
                }
            });
        if largest < 0.0 {
            for projection in projections.iter_mut().skip(component).step_by(num_comps) {
                *projection = -*projection;
            }
            if let Some(weights) = princomps.chunks_exact_mut(num_cols).nth(component) {
                for weight in weights {
                    *weight = -*weight;
                }
            }
        }
    }
}

/// The positions of the first ten traits of a list, as text, with how many
/// more there are: `the position 3`, `the positions 0, 3, 7`, or `the
/// positions 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, and 2 more`.
///
/// The word comes with the numbers because a message that ends in one
/// number, `the traits at 1`, reads as a count of traits to whoever gets
/// it in TypeScript, where no name replaces it.
/// [`Error::PcaTraitsWithNoVariance`] carries the positions and the Python
/// layer puts the name of each trait in their place, as pyNei's message
/// has it, so this is what a reader of the message in Rust or in
/// TypeScript gets.
pub(crate) fn the_positions_listed(positions: &[usize]) -> String {
    let word = if positions.len() == 1 {
        "the position"
    } else {
        "the positions"
    };
    let shown = positions
        .iter()
        .take(10)
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    match positions.len().saturating_sub(10) {
        0 => format!("{word} {shown}"),
        more => format!("{word} {shown}, and {more} more"),
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::{Pca, PcaOptions, TraitScale, fix_the_signs, pca};
    use crate::error::Error;

    /// The tolerance of "How it is verified" of `docs/specs/pca.md`: every
    /// literal here and in the reference files is written with 12
    /// significant digits, and the largest of them is a projection of iris
    /// of 3.7956, so a literal is within 1e-11 of the number it stands
    /// for.
    const TOLERANCE: f64 = 1e-9;

    /// Centered and standardized, which is what `do_pca` does by default.
    const STANDARDIZED: PcaOptions = PcaOptions {
        center: true,
        standardize: true,
    };

    /// Centered and not standardized.
    const CENTERED: PcaOptions = PcaOptions {
        center: true,
        standardize: false,
    };

    /// Neither centered nor standardized.
    const AS_IT_IS: PcaOptions = PcaOptions {
        center: false,
        standardize: false,
    };

    /// The 3 rows x 5 traits of "How it is verified" of "The PCA of a
    /// table" of `docs/specs/pca.md`, which has fewer rows than traits, so
    /// that the product is the 3 x 3 one of the rows.
    const FIVE_TRAITS: [f64; 15] = [
        1.0, 2.0, 3.0, 4.0, 5.0, //
        2.0, 4.0, 1.0, 3.0, 2.0, //
        5.0, 1.0, 4.0, 2.0, 6.0,
    ];

    /// The table of `test_pca_refuses_traits_with_no_variance` of pyNei,
    /// whose traits `a`, `fixed` and `b` are the columns, the second of
    /// them having no variance.
    const ONE_TRAIT_FIXED: [f64; 9] = [
        1.0, 5.0, 3.0, //
        2.0, 5.0, 1.0, //
        3.0, 5.0, 2.0,
    ];

    /// The 5 rows x 3 traits of "How it is verified", which has more rows
    /// than traits, as iris has, and loses a component, which iris does
    /// not: its second trait has no variance.
    const FIVE_ROWS: [f64; 15] = [
        1.0, 5.0, 3.0, //
        2.0, 5.0, 1.0, //
        3.0, 5.0, 2.0, //
        4.0, 5.0, 9.0, //
        7.0, 5.0, 2.0,
    ];

    /// The 3 rows x 3 traits of "How it is verified" whose values are
    /// multiplied by 1e153 or by 2.5e153, which puts the eigenvalues of
    /// its product near the largest `f64`.
    const NEAR_THE_LARGEST: [f64; 9] = [
        1.0, 2.0, 3.0, //
        2.0, 4.0, 1.0, //
        5.0, 1.0, 4.0,
    ];

    /// The rows of one of the files of `tests/reference/pca/`, without the
    /// first field of each line, which names the row, and without the
    /// first line when the file has a header. The values of R's files
    /// carry leading spaces.
    fn the_reference(name: &str, with_a_header: bool) -> Vec<Vec<f64>> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/pca")
            .join(name);
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("{path}: {error}", path = path.display()));
        text.lines()
            .skip(usize::from(with_a_header))
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                line.split('\t')
                    .skip(1)
                    .map(|field| {
                        field
                            .trim()
                            .parse()
                            .unwrap_or_else(|error| panic!("{name}: `{field}`: {error}"))
                    })
                    .collect()
            })
            .collect()
    }

    /// The same rows, one after another, which is how [`Pca`] holds a
    /// matrix.
    fn row_after_row(rows: Vec<Vec<f64>>) -> Vec<f64> {
        rows.into_iter().flatten().collect()
    }

    /// Each value against the one of the reference, within the tolerance.
    fn assert_close(got: &[f64], expected: &[f64], tolerance: f64, what: &str) {
        assert_eq!(got.len(), expected.len(), "{what}: the count of values");
        for (position, (value, reference)) in got.iter().zip(expected).enumerate() {
            assert!(
                (value - reference).abs() <= tolerance,
                "{what}: the value at {position} is {value} and the reference is {reference}"
            );
        }
    }

    /// What R's `prcomp` gives for iris, from the three files of that
    /// name: the projections, the percentages and the weights, each row
    /// after row.
    fn the_iris_reference(name: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        (
            row_after_row(the_reference(&format!("{name}.r.projections.tsv"), true)),
            row_after_row(the_reference(&format!("{name}.r.percent.tsv"), false)),
            row_after_row(the_reference(&format!("{name}.r.princomps.tsv"), true)),
        )
    }

    /// Every component of the result against a reference held row after
    /// row, with the shape of the result.
    fn assert_the_result_is(
        result: &Pca,
        num_rows: usize,
        num_cols: usize,
        num_comps: usize,
        reference: (&[f64], &[f64], &[f64]),
        what: &str,
    ) {
        assert_eq!(result.num_rows, num_rows, "{what}: the rows");
        assert_eq!(result.num_cols, num_cols, "{what}: the traits");
        assert_eq!(result.num_comps, num_comps, "{what}: the components");
        assert_eq!(
            result.num_prin_comps, num_comps,
            "{what}: the components the weights are given for"
        );
        assert_eq!(
            result.used_cols,
            (0..num_cols).collect::<Vec<usize>>(),
            "{what}: every trait of a table is used"
        );
        assert_close(
            &result.projections,
            reference.0,
            TOLERANCE,
            &format!("{what}: the projections"),
        );
        assert_close(
            &result.explained_variance_percent,
            reference.1,
            TOLERANCE,
            &format!("{what}: the percentages"),
        );
        assert_close(
            &result.princomps,
            reference.2,
            TOLERANCE,
            &format!("{what}: the weights"),
        );
    }

    /// The table of iris, 150 rows x 4 traits, which
    /// `tests/reference/pca/make_reference.py` writes from pyNei's
    /// `test/datasets.py`.
    fn the_iris_table() -> Vec<f64> {
        row_after_row(the_reference("iris.tsv", true))
    }

    /// Iris standardized is the first of the two runs of "How it is
    /// verified" of "The PCA of a table", and it has 150 rows and 4
    /// traits, so its product is the 4 x 4 one of the traits. R's
    /// projections in the file are multiplied by sqrt(n / (n - 1)), so
    /// they are pyNei's and popnei's, the ones of a standard deviation
    /// with n in it.
    #[test]
    fn iris_standardized_gives_the_numbers_of_r() {
        let result = pca(&the_iris_table(), 150, 4, &STANDARDIZED).expect("the analysis of iris");
        let (projections, percent, princomps) = the_iris_reference("iris");
        assert_the_result_is(
            &result,
            150,
            4,
            4,
            (&projections, &percent, &princomps),
            "iris standardized",
        );
    }

    /// The second run of the same part, which R computes with
    /// `scale. = FALSE` and whose projections are R's as they are.
    #[test]
    fn iris_not_standardized_gives_the_numbers_of_r() {
        let result = pca(&the_iris_table(), 150, 4, &CENTERED).expect("the analysis of iris");
        let (projections, percent, princomps) = the_iris_reference("iris_not_standardized");
        assert_the_result_is(
            &result,
            150,
            4,
            4,
            (&projections, &percent, &princomps),
            "iris not standardized",
        );
    }

    /// The table of 3 rows x 5 traits of "How it is verified", which has
    /// fewer rows than traits, so that the matrix that is decomposed is
    /// the 3 x 3 product of the rows and not the 5 x 5 product of the
    /// traits. Centering takes one of the three dimensions of the rows
    /// out, so it has 2 components and not 3.
    #[test]
    fn a_table_of_fewer_rows_than_traits_standardized_gives_the_numbers_of_numpy() {
        let result = pca(&FIVE_TRAITS, 3, 5, &STANDARDIZED).expect("the analysis of the table");
        let projections = [
            -0.347154191646,
            1.62410767053,
            -2.14982434306,
            -0.994055030603,
            2.49697853471,
            -0.630052639924,
        ];
        let percent = [73.1810836106, 26.8189163894];
        let princomps = [
            0.420101972196,
            -0.496432942021,
            0.496432942021,
            -0.317325807736,
            0.47950738561,
            -0.513968704337,
            -0.270671721206,
            0.270671721206,
            0.686274627813,
            0.344001373342,
        ];
        assert_the_result_is(
            &result,
            3,
            5,
            2,
            (&projections, &percent, &princomps),
            "the table of five traits, standardized",
        );
    }

    /// The same table with neither step, which keeps all 3 components,
    /// since nothing takes a dimension out. It is the run that no program
    /// outside the project gives a number for.
    #[test]
    fn a_table_that_is_neither_centered_nor_standardized_keeps_every_component() {
        let result = pca(&FIVE_TRAITS, 3, 5, &AS_IT_IS).expect("the analysis of the table");
        let projections = [
            7.07789159873,
            1.12190772301,
            1.90912901021,
            4.87213205354,
            2.87255560109,
            -1.41801042717,
            8.66041659868,
            -2.53292797372,
            -0.762535387595,
        ];
        let percent = [87.0392022767, 9.31343669027, 3.64736103304];
        let princomps = [
            0.403960199411,
            0.28423522248,
            0.408147561391,
            0.404797068091,
            0.652365999566,
            -0.364035502371,
            0.703323259808,
            -0.244470602227,
            0.504800545607,
            -0.241298733997,
            -0.759913160568,
            -0.419484427674,
            0.201897964367,
            0.297806276491,
            0.342218405411,
        ];
        assert_the_result_is(
            &result,
            3,
            5,
            3,
            (&projections, &percent, &princomps),
            "the table of five traits, as it is",
        );
    }

    /// pyNei's table of three traits, one of which is fixed, without
    /// standardizing, which pyNei gives 3 components and popnei 2: the
    /// third has no variance. The fixed trait is still a column of the
    /// weights, with a weight of 0.
    ///
    /// The second component is given up to its sign, as "How it is
    /// verified" says: its two largest projections are the same number
    /// with opposite signs, equal bit for bit on Accelerate's LAPACK and
    /// one bit apart on faer, so which of the two the sign rule finds is
    /// decided by the rounding of the eigendecomposition and is not the
    /// same on every backend.
    #[test]
    #[expect(
        clippy::approx_constant,
        reason = "the numbers of the spec, written with the 12 significant digits of the reference, and four of them are the first digits of sqrt(2) and of its inverse"
    )]
    fn a_trait_with_no_variance_gets_a_weight_of_0_and_no_component_of_its_own() {
        let result =
            pca(&ONE_TRAIT_FIXED, 3, 3, &CENTERED).expect("the analysis of the fixed trait");
        assert_eq!(result.num_comps, 2);
        assert_eq!(result.num_cols, 3);
        assert_eq!(result.used_cols, vec![0, 1, 2]);
        let first_of_the_projections = [
            result.projections.first().copied().expect("row 0"),
            result.projections.get(2).copied().expect("row 1"),
            result.projections.get(4).copied().expect("row 2"),
        ];
        assert_close(
            &first_of_the_projections,
            &[1.41421356237, -0.707106781187, -0.707106781187],
            TOLERANCE,
            "the projections of the first component",
        );
        assert_close(
            &result.explained_variance_percent,
            &[75.0, 25.0],
            TOLERANCE,
            "the percentages",
        );
        assert_close(
            result.princomps.get(..3).expect("the first component"),
            &[-0.707106781187, 0.0, 0.707106781187],
            TOLERANCE,
            "the weights of the first component",
        );
        // The second component, up to its sign.
        let second: Vec<f64> = result
            .projections
            .iter()
            .skip(1)
            .step_by(2)
            .map(|value| value.abs())
            .collect();
        assert_close(
            &second,
            &[0.0, 0.707106781187, 0.707106781187],
            TOLERANCE,
            "the projections of the second component",
        );
        let weights: Vec<f64> = result
            .princomps
            .get(3..6)
            .expect("the second component")
            .iter()
            .map(|weight| weight.abs())
            .collect();
        assert_close(
            &weights,
            &[0.707106781187, 0.0, 0.707106781187],
            TOLERANCE,
            "the weights of the second component",
        );
    }

    /// The sign rule of "What both analyses compute": in every component
    /// the projection of the largest absolute value is positive. Whichever
    /// sign the eigendecomposition gave, this is what comes out, so that
    /// Python natively, Python under pyodide and TypeScript give the same
    /// numbers.
    #[test]
    fn the_projection_of_the_largest_absolute_value_is_positive_in_every_component() {
        for (data, num_rows, num_cols, options, what) in [
            (&the_iris_table()[..], 150, 4, STANDARDIZED, "iris"),
            (&FIVE_TRAITS[..], 3, 5, STANDARDIZED, "five traits"),
            (&FIVE_TRAITS[..], 3, 5, AS_IT_IS, "five traits as they are"),
            (&ONE_TRAIT_FIXED[..], 3, 3, CENTERED, "one trait fixed"),
        ] {
            let result = pca(data, num_rows, num_cols, &options).expect("the analysis");
            assert!(result.num_comps > 0, "{what}: no component");
            for component in 0..result.num_comps {
                let largest = result
                    .projections
                    .iter()
                    .skip(component)
                    .step_by(result.num_comps)
                    .copied()
                    .fold(0.0_f64, |largest: f64, value| {
                        if value.abs().total_cmp(&largest.abs()) == std::cmp::Ordering::Greater {
                            value
                        } else {
                            largest
                        }
                    });
                assert!(
                    largest > 0.0,
                    "{what}: the largest projection of component {component} is {largest}"
                );
            }
        }
    }

    /// The other half of the rule: when two individuals have the same
    /// absolute value it is the first of them that is made positive. The
    /// two projections here are exactly opposite, which the
    /// eigendecomposition of a table seldom gives, so the rule is checked
    /// on the step that applies it.
    #[test]
    fn the_sign_rule_breaks_a_tie_towards_the_first_individual() {
        let mut projections = vec![1.0, -1.0];
        let mut princomps = vec![0.5, -0.25];
        fix_the_signs(&mut projections, &mut princomps, 1, 2);
        assert_close(&projections, &[1.0, -1.0], 0.0, "the first is positive");
        assert_close(&princomps, &[0.5, -0.25], 0.0, "the weights are untouched");

        let mut projections = vec![-1.0, 1.0];
        let mut princomps = vec![0.5, -0.25];
        fix_the_signs(&mut projections, &mut princomps, 1, 2);
        assert_close(&projections, &[1.0, -1.0], 0.0, "the component is turned");
        assert_close(&princomps, &[-0.5, 0.25], 0.0, "the weights are turned");
    }

    /// A value that is not finite is refused with the place where it is.
    /// pyNei refuses a NaN and lets an infinity reach numpy's SVD, which
    /// raises `LinAlgError`.
    #[test]
    fn a_value_that_is_not_finite_is_refused_with_its_place() {
        let mut table = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        *table.get_mut(5).expect("the last value") = f64::INFINITY;
        match pca(&table, 2, 3, &CENTERED) {
            Err(Error::PcaValueNotFinite { row, col, value }) => {
                assert_eq!(row, 1);
                assert_eq!(col, 2);
                assert!(value.is_infinite(), "{value}");
            }
            other => panic!("an infinity was taken: {other:?}"),
        }

        *table.get_mut(5).expect("the last value") = 6.0;
        *table.get_mut(1).expect("the second value") = f64::NAN;
        match pca(&table, 2, 3, &CENTERED) {
            Err(Error::PcaValueNotFinite { row, col, value }) => {
                assert_eq!(row, 0);
                assert_eq!(col, 1);
                assert!(value.is_nan(), "{value}");
            }
            other => panic!("a NaN was taken: {other:?}"),
        }
    }

    /// Standardizing divides by the standard deviation of the centered
    /// trait, so it needs the centering. pyNei refuses the pair too.
    #[test]
    fn standardizing_without_centering_is_refused() {
        let options = PcaOptions {
            center: false,
            standardize: true,
        };
        let result = pca(&FIVE_TRAITS, 3, 5, &options);
        assert!(
            matches!(result, Err(Error::PcaStandardizeWithoutCentering)),
            "{result:?}"
        );
    }

    /// One row has no variation to find directions in, and no trait has
    /// nothing to project.
    #[test]
    fn a_table_of_fewer_than_two_rows_or_of_no_traits_is_refused() {
        for (num_rows, num_cols) in [(1_usize, 3_usize), (0, 3), (3, 0)] {
            let result = pca(&ONE_TRAIT_FIXED, num_rows, num_cols, &CENTERED);
            match result {
                Err(Error::PcaTableTooSmall {
                    num_rows: rows,
                    num_cols: cols,
                }) => {
                    assert_eq!(rows, num_rows);
                    assert_eq!(cols, num_cols);
                }
                other => panic!("a table of {num_rows} x {num_cols} was taken: {other:?}"),
            }
        }
    }

    /// The traits that have no variance are refused when the table is
    /// standardized, with the position of each one, which the Python layer
    /// turns into its name. The message names the first ten and says how
    /// many more there are, as pyNei's does.
    #[test]
    fn the_traits_with_no_variance_are_refused_with_their_positions() {
        match pca(&ONE_TRAIT_FIXED, 3, 3, &STANDARDIZED) {
            Err(Error::PcaTraitsWithNoVariance {
                positions,
                num_cols,
            }) => {
                assert_eq!(positions, vec![1]);
                assert_eq!(num_cols, 3);
                let message = Error::PcaTraitsWithNoVariance {
                    positions,
                    num_cols,
                }
                .to_string();
                assert!(message.contains("1 of the 3 traits"), "{message}");
                assert!(message.contains("no variance"), "{message}");
                // The position comes with the word, so that a message
                // that ends in one number does not read as a count.
                assert!(message.contains("at the position 1,"), "{message}");
            }
            other => panic!("a fixed trait was standardized: {other:?}"),
        }

        // Twelve traits of two rows, all of them fixed, which is what a
        // message that names the first ten needs.
        let table = vec![7.0; 24];
        match pca(&table, 2, 12, &STANDARDIZED) {
            Err(Error::PcaTraitsWithNoVariance {
                positions,
                num_cols,
            }) => {
                assert_eq!(positions, (0..12).collect::<Vec<usize>>());
                assert_eq!(num_cols, 12);
                let message = Error::PcaTraitsWithNoVariance {
                    positions,
                    num_cols,
                }
                .to_string();
                assert!(message.contains("12 of the 12 traits"), "{message}");
                assert!(message.contains("at the positions 0, 1,"), "{message}");
                assert!(message.contains("and 2 more"), "{message}");
                assert!(!message.contains("10, 11"), "{message}");
            }
            other => panic!("twelve fixed traits were standardized: {other:?}"),
        }
    }

    /// A buffer that does not hold the values of the table it was said to
    /// be. No argument of the Python or the TypeScript function reaches
    /// this: each binding crate takes the two numbers from the array it
    /// was given.
    #[test]
    fn a_buffer_that_does_not_hold_the_table_is_refused() {
        match pca(&ONE_TRAIT_FIXED, 4, 3, &CENTERED) {
            Err(Error::PcaTableOfAnotherSize {
                num_values,
                num_rows,
                num_cols,
            }) => {
                assert_eq!(num_values, 9);
                assert_eq!(num_rows, 4);
                assert_eq!(num_cols, 3);
            }
            other => panic!("a buffer of 9 values was read as 4 x 3: {other:?}"),
        }

        // One value too many is refused as well: the analysis of the first
        // nine would be of a table the caller did not mean.
        let mut longer = ONE_TRAIT_FIXED.to_vec();
        longer.push(0.0);
        match pca(&longer, 3, 3, &CENTERED) {
            Err(Error::PcaTableOfAnotherSize {
                num_values,
                num_rows,
                num_cols,
            }) => {
                assert_eq!(num_values, 10);
                assert_eq!(num_rows, 3);
                assert_eq!(num_cols, 3);
            }
            other => panic!("a buffer of 10 values was read as 3 x 3: {other:?}"),
        }
    }

    /// The table of 5 rows x 3 traits of "How it is verified", which has
    /// more rows than traits, so that the matrix that is decomposed is the
    /// 3 x 3 product of the traits, and which loses a component, its
    /// second trait having no variance. Iris, the other table with more
    /// rows than traits, keeps all four of its components, so nothing but
    /// this pins that the weights of that side are cut to the components
    /// that have variance.
    #[test]
    fn a_table_of_more_rows_than_traits_drops_the_component_with_no_variance() {
        let result = pca(&FIVE_ROWS, 5, 3, &CENTERED).expect("the analysis of the table");
        let projections = [
            -0.765373820993,
            -2.30958933885,
            -2.58720936972,
            -1.01308818829,
            -1.44494156963,
            -0.179287089168,
            5.62553292806,
            -0.270886092928,
            -0.828008167714,
            3.77285070924,
        ];
        let percent = [66.8261599339, 33.1738400661];
        let princomps = [
            0.15423335048,
            0.0,
            0.988034449602,
            0.988034449602,
            0.0,
            -0.15423335048,
        ];
        assert_the_result_is(
            &result,
            5,
            3,
            2,
            (&projections, &percent, &princomps),
            "the table of five rows",
        );
    }

    /// A table whose eigenvalues are near the largest `f64`. The
    /// percentages are what the table gives unscaled, and the count of the
    /// components does not change with the scale: both are computed in the
    /// order that does not overflow, the share of the total before the 100
    /// and the tolerance of the threshold before the largest eigenvalue.
    /// At 1e153 the largest eigenvalue is 1.4e307, so 100 times it is an
    /// infinity; at 2.5e153 it is 8.9e307, so it times the larger side of
    /// the table is an infinity as well.
    #[test]
    fn a_table_whose_eigenvalues_are_near_the_largest_float_gives_its_percentages() {
        for scale in [1e153, 2.5e153] {
            let table: Vec<f64> = NEAR_THE_LARGEST.iter().map(|value| value * scale).collect();
            let result = pca(&table, 3, 3, &CENTERED).expect("the analysis of the large table");
            assert_eq!(result.num_comps, 2, "the components at the scale {scale}");
            assert_close(
                &result.explained_variance_percent,
                &[78.8675134595, 21.1324865405],
                TOLERANCE,
                &format!("the percentages at the scale {scale}"),
            );
        }
    }

    /// A trait whose values sum above the largest `f64` has a mean that is
    /// not finite, and every centered value of it would be a NaN. It is
    /// found whenever the table is centered, standardized or not.
    #[test]
    fn a_trait_whose_mean_is_not_finite_is_refused() {
        let table = [f64::MAX, 1.0, f64::MAX, 2.0];
        match pca(&table, 2, 2, &CENTERED) {
            Err(Error::PcaTraitOutOfRange { position, problem }) => {
                assert_eq!(position, 0);
                assert_eq!(problem, TraitScale::MeanNotFinite);
                let message = Error::PcaTraitOutOfRange { position, problem }.to_string();
                assert!(message.contains("the trait at the position 0"), "{message}");
                assert!(message.contains("mean"), "{message}");
            }
            other => panic!("a trait whose mean is an infinity was centered: {other:?}"),
        }
    }

    /// A trait whose squared deviations sum above the largest `f64` has a
    /// standard deviation that is not finite, and dividing by it would
    /// make the trait a column of zeros, which looks like a trait with no
    /// variance and would leave the analysis with a weight of 0 and no
    /// word.
    #[test]
    fn a_trait_whose_standard_deviation_is_not_finite_is_refused() {
        let table = [1e154, 1.0, -1e154, 2.0];
        match pca(&table, 2, 2, &STANDARDIZED) {
            Err(Error::PcaTraitOutOfRange { position, problem }) => {
                assert_eq!(position, 0);
                assert_eq!(problem, TraitScale::DeviationNotFinite);
            }
            other => panic!("a trait whose deviation is an infinity was standardized: {other:?}"),
        }
    }

    /// A trait whose squared deviations all fall below the smallest `f64`
    /// above 0 has a standard deviation of 0 although its values differ,
    /// and dividing by it would give infinities, which the linear algebra
    /// would then refuse as a defect of popnei.
    #[test]
    fn a_trait_whose_standard_deviation_falls_to_zero_is_refused() {
        let table = [1e-200, 1.0, 2e-200, 2.0];
        match pca(&table, 2, 2, &STANDARDIZED) {
            Err(Error::PcaTraitOutOfRange { position, problem }) => {
                assert_eq!(position, 0);
                assert_eq!(problem, TraitScale::DeviationOfZero);
                let message = Error::PcaTraitOutOfRange { position, problem }.to_string();
                assert!(message.contains("not all equal"), "{message}");
            }
            other => panic!("a trait whose deviation is 0 was standardized: {other:?}"),
        }
    }

    /// A table with no direction to give: every trait the same value in
    /// every row, which centering turns into zeros, and a table of zeros
    /// that is not centered. pyNei gives 0 for every projection and a
    /// percentage of NaN for every component.
    #[test]
    fn a_table_in_which_no_trait_has_variance_is_refused() {
        let fixed = [1.0, 5.0, 1.0, 5.0, 1.0, 5.0];
        let result = pca(&fixed, 3, 2, &CENTERED);
        assert!(
            matches!(result, Err(Error::PcaNoTraitWithVariance)),
            "{result:?}"
        );
        if let Err(error) = result {
            let message = error.to_string();
            assert!(message.contains("no trait has variance"), "{message}");
        }

        let zeros = [0.0; 6];
        let result = pca(&zeros, 3, 2, &AS_IT_IS);
        assert!(
            matches!(result, Err(Error::PcaNoTraitWithVariance)),
            "{result:?}"
        );
    }

    /// An error of the linalg crate becomes one of popnei, with the
    /// operation that was being done. The values here are finite and their
    /// products are not, so the crate refuses the matrix it is asked to
    /// decompose.
    #[test]
    fn an_error_of_the_linear_algebra_is_wrapped_with_the_operation() {
        let table = [1e200, 2e200, 3e200, 4e200];
        match pca(&table, 2, 2, &AS_IT_IS) {
            Err(Error::PcaLinalg { operation, source }) => {
                let message = Error::PcaLinalg { operation, source }.to_string();
                assert!(message.contains(operation), "{message}");
                assert!(message.contains("finite"), "{message}");
            }
            other => panic!("a product that is not finite was taken: {other:?}"),
        }
    }
}
