//! The projection matrix of a mixed model, kept as a factor a variant is
//! solved against rather than as the matrix itself.
//!
//! Both mixed models take every variant through `p = v⁻¹ - v⁻¹ d
//! (d' v⁻¹ d)⁻¹ d' v⁻¹`, with `v` the covariance of the trait under the
//! null and `d` the design, which takes the covariates out of the variant
//! and weights it by that covariance. Written that way it is an
//! individuals by individuals matrix and a variant costs a product with
//! it, which is 59 to 65 per cent of what the three passes of the two
//! mixed models spend, measured in
//! `docs/reports/perf-gwas-2026-09-24.md`.
//!
//! [`TheProjection`] keeps instead the two pieces `p` is built from. With
//! `l` the lower triangular Cholesky factor of `v`, so that `v = l l'`,
//! and `q` the columns of length 1 at right angles that span what `l⁻¹`
//! makes of the design, `p = m' m` with `m = (i - q q') l⁻¹`. So
//! `x' p x` is `‖m x‖²`, the squared length of one vector of as many
//! values as there are individuals, and `m x` costs a triangular solve
//! where the matrix costs a product: half the arithmetic, a triangular
//! matrix having half the entries. That is the finding H3 of that report.
//!
//! What it buys beside the arithmetic is that the quantities a variant's
//! test turns on stop being differences. `x' p x` through the matrix is a
//! sum of terms of either sign that cancels to rounding for a variant the
//! design explains, and the rounding falls on either side of 0 and not on
//! the same side in the two arithmetic backends; `‖m x‖²` is a sum of
//! squares and is 0 or above whatever the rounding does. The linear mixed
//! model's Wald test gains the same on what the variant leaves of the
//! trait: `y' p y - num² / den` is `‖m y - beta m x‖²`, the trait through
//! the factor less the variant through it times the variant's effect,
//! formed from the residuals of the one variant as the
//! plain linear model already forms its own. Both are **Open 2** of
//! `docs/specs/gwas.md`.

use popnei_linalg::{TheFirstOperand, TheHalfThatHoldsTheMatrix, TheSecondOperand};

use crate::error::{Error, Result};

use super::study::Design;

/// The projection matrix of a mixed model, held as the Cholesky factor of
/// the covariance of the trait and the directions the design spans
/// through it.
///
/// [`TheProjection::through`] is what a block of variants goes through,
/// and it leaves `m x` for each of them, whose squared length is
/// `x' p x`. The matrix `p` is never formed: what is kept is one
/// individuals by individuals triangular factor and one individuals by
/// columns of the design matrix, where `p` was individuals by individuals
/// and dense.
pub(crate) struct TheProjection {
    /// `l`, `num_individuals` x `num_individuals`, row after row, the
    /// lower triangular Cholesky factor of the covariance of the trait
    /// under the null. Only its lower half is read, and its upper half
    /// holds whatever the factorization left there.
    factor: Vec<f64>,
    /// `q`, `num_individuals` x `num_coefs`, row after row: columns of
    /// length 1 at right angles to each other spanning what `l⁻¹` makes
    /// of the design.
    of_the_design: Vec<f64>,
    /// Each row of the block in hand against each column of
    /// `of_the_design`, `num_vars` x `num_coefs`, which is what is taken
    /// out of that row along the columns of the design. It is a field so
    /// that a pass over a million variants asks the machine for it once.
    ///
    /// The row of the block it is taken of is the variant after the solve
    /// against the factor and before the design is taken out, which the
    /// functions below call `t`.
    along_the_design: Vec<f64>,
    /// How many individuals the study tests.
    num_individuals: usize,
    /// How many columns the design has.
    num_coefs: usize,
    /// The largest value of the diagonal of `p` itself, which is what a
    /// variant's own squared length is weighted by to say how much of the
    /// variant the projection has left.
    ///
    /// It is the scale of the third row of the table of **Open 2** of
    /// `docs/specs/gwas.md`. What bounds `x' p x` over `x' x` is the
    /// largest eigenvalue of the projection, and this is not that: no
    /// value of the diagonal is below 0, the matrix being 0 or above as a
    /// quadratic form, and the largest of them is at most that
    /// eigenvalue, 1.697 times below it on `panel_called` and 1.725 times
    /// on `panel`, measured on 24 September 2026. So the threshold built
    /// on it is about 1.7 times tighter than it was meant to be and not
    /// looser, and it costs one walk over the individuals where the
    /// eigenvalue costs a decomposition.
    largest_of_the_diagonal: f64,
}

impl TheProjection {
    /// The projection of a null model whose covariance has already been
    /// factored, `factor` holding `l` in its lower half.
    ///
    /// `factor` is `num_individuals` x `num_individuals`, row after row,
    /// as [`popnei_linalg::cholesky_lower`] leaves it, and `design` is the
    /// design of the same individuals. The buffer is taken and kept, every
    /// variant of the pass being solved against it.
    ///
    /// What is worked out here, once for a study, is `q` and the largest
    /// value of the diagonal of `p`. Both come from `l⁻ᵀ`, which is
    /// written into a buffer of the size of `factor` and dropped before
    /// this returns: the entry i, i of `p = m' m` is the squared length of
    /// `(i - q q')` applied to the column i of `l⁻¹`, and `q q'` takes a
    /// squared length away rather than adding one, `q` being at right
    /// angles to itself, so that entry is `‖row i of l⁻ᵀ‖²` less
    /// `‖row i of l⁻ᵀ q‖²`.
    ///
    /// # Errors
    ///
    /// [`Error::GwasLinalg`] when `factor` holds a 0 on its diagonal, so
    /// that the covariance it factors is singular, when `factor` or the
    /// design is of another size than the study, or when one of the two
    /// solves, the factorization of what `l⁻¹` makes of the design or one
    /// of the two products could not be done.
    pub(crate) fn of_the_factored_covariance(
        factor: Vec<f64>,
        design: &Design<'_>,
    ) -> Result<TheProjection> {
        let num_individuals = design.num_individuals();
        let num_coefs = design.num_coefs();
        // `l⁻ᵀ`: the solve takes one right hand side per row, so the row
        // i of the identity comes back as the column i of `l⁻¹`, which is
        // the row i of `l⁻ᵀ`.
        let mut inverse_of_the_factor = vec![0.0_f64; factor.len()];
        for (at, row) in inverse_of_the_factor
            .chunks_exact_mut(num_individuals.max(1))
            .enumerate()
        {
            if let Some(entry) = row.get_mut(at) {
                *entry = 1.0;
            }
        }
        popnei_linalg::solve_triangular(
            &factor,
            num_individuals,
            TheHalfThatHoldsTheMatrix::TheLowerHalf,
            &mut inverse_of_the_factor,
            num_individuals,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "inverse of the factor of the covariance of the trait",
            source,
        })?;
        // What `l⁻¹` makes of the design, one row per individual. The
        // solve wants one right hand side per row, which is one per column
        // of the design, so the design goes in laid out by column and is
        // turned back after it.
        let mut by_column = vec![0.0_f64; design.values().len()];
        the_matrix_the_other_way_round(design.values(), num_individuals, num_coefs, &mut by_column);
        popnei_linalg::solve_triangular(
            &factor,
            num_individuals,
            TheHalfThatHoldsTheMatrix::TheLowerHalf,
            &mut by_column,
            num_coefs,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "solve of the design against the factor of the covariance",
            source,
        })?;
        let mut of_the_factor = vec![0.0_f64; design.values().len()];
        the_matrix_the_other_way_round(&by_column, num_coefs, num_individuals, &mut of_the_factor);
        let of_the_design = popnei_linalg::thin_qr(&of_the_factor, num_individuals, num_coefs)
            .map_err(|source| Error::GwasLinalg {
                operation: "factorization of what the factor's inverse makes of the design",
                source,
            })?
            .q;
        // `l⁻ᵀ q`, one row per individual, whose row i is what the design
        // takes off the squared length of the column i of `l⁻¹`.
        let mut of_the_design_through_the_factor = vec![0.0_f64; design.values().len()];
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &inverse_of_the_factor,
                rows: num_individuals,
            },
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &of_the_design,
                cols: num_coefs,
            },
            &mut of_the_design_through_the_factor,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the factor's inverse with the directions of the design",
            source,
        })?;
        let largest_of_the_diagonal = inverse_of_the_factor
            .chunks_exact(num_individuals.max(1))
            .zip(of_the_design_through_the_factor.chunks_exact(num_coefs.max(1)))
            .map(|(of_the_individual, of_the_design)| {
                let had = of_the_individual
                    .iter()
                    .map(|value| value * value)
                    .sum::<f64>();
                let explained = of_the_design.iter().map(|value| value * value).sum::<f64>();
                had - explained
            })
            .fold(0.0_f64, f64::max);
        Ok(TheProjection {
            factor,
            of_the_design,
            along_the_design: Vec::new(),
            num_individuals,
            num_coefs,
            largest_of_the_diagonal,
        })
    }

    /// A projection built from its pieces, which is what a test that
    /// wants a chosen `p` and not a fitted one has.
    ///
    /// The two matrices a study's own projection is built from are a
    /// triangular factor and a set of directions at right angles, and a
    /// test that wants `p` to be the identity, or 0, gives them straight:
    /// the identity and one direction of all zeros for the first, since
    /// `m` is then the identity, and the identity twice for the second,
    /// since `m` is then `i - i`. Neither is a projection a fit produces,
    /// which is why this is not the constructor a study calls.
    #[cfg(test)]
    pub(crate) fn of_the_factor_and_the_directions(
        factor: Vec<f64>,
        of_the_design: Vec<f64>,
        num_individuals: usize,
        num_coefs: usize,
        largest_of_the_diagonal: f64,
    ) -> TheProjection {
        TheProjection {
            factor,
            of_the_design,
            along_the_design: Vec::new(),
            num_individuals,
            num_coefs,
            largest_of_the_diagonal,
        }
    }

    /// How many individuals the study tests.
    #[must_use]
    pub(crate) fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// The largest value of the diagonal of `p`.
    #[must_use]
    pub(crate) fn largest_of_the_diagonal(&self) -> f64 {
        self.largest_of_the_diagonal
    }

    /// Every row of `dosages` through the projection, `m x`, left in
    /// `into` as `num_vars` rows of one value per tested individual.
    ///
    /// `dosages` holds one row of `num_individuals` values for each of the
    /// `num_vars` variants, row after row, and may hold more rows than
    /// that, of which the first `num_vars` are taken. `into` is
    /// overwritten and grown to `num_vars` times the individuals, and the
    /// caller keeps it from block to block so that a pass allocates once.
    ///
    /// The squared length of a row is `x' p x` and its product with `m y`
    /// is `x' p y`, so the two quantities a test turns on are a dot
    /// product each and neither of them is a difference.
    ///
    /// # Errors
    ///
    /// [`Error::GwasVariantsTooLarge`] when the values of the block are
    /// more than a `usize` counts or `dosages` holds fewer of them than
    /// `num_vars` rows. [`Error::GwasLinalg`] when the solve or the
    /// product could not be done, which is where a block of other
    /// individuals than the null model was fitted over is refused.
    pub(crate) fn through(
        &mut self,
        dosages: &[f64],
        num_vars: usize,
        into: &mut Vec<f64>,
    ) -> Result<()> {
        let values = num_vars
            .checked_mul(self.num_individuals)
            .ok_or(Error::GwasVariantsTooLarge)?;
        let rows = dosages.get(..values).ok_or(Error::GwasVariantsTooLarge)?;
        into.clear();
        into.extend_from_slice(rows);
        if num_vars == 0 {
            return Ok(());
        }
        popnei_linalg::solve_triangular(
            &self.factor,
            self.num_individuals,
            TheHalfThatHoldsTheMatrix::TheLowerHalf,
            into,
            num_vars,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "solve of a block of variants against the factor of the covariance",
            source,
        })?;
        let along = num_vars
            .checked_mul(self.num_coefs)
            .ok_or(Error::GwasVariantsTooLarge)?;
        self.along_the_design.clear();
        self.along_the_design.resize(along, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: into,
                rows: num_vars,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.of_the_design,
                cols: self.num_coefs,
            },
            &mut self.along_the_design,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of a solved block with the directions of the design",
            source,
        })?;
        the_design_taken_out(
            into,
            &self.along_the_design,
            &self.of_the_design,
            self.num_individuals,
            self.num_coefs,
        );
        Ok(())
    }
}

/// `into` gets the `rows` x `cols` matrix `values` laid out `cols` x
/// `rows`, which is how a solve takes one right hand side per column of a
/// matrix and how the factorization that follows takes it back.
fn the_matrix_the_other_way_round(values: &[f64], rows: usize, cols: usize, into: &mut [f64]) {
    for (at, of_the_row) in values.chunks_exact(cols.max(1)).enumerate().take(rows) {
        for (of_the_col, value) in into.chunks_exact_mut(rows.max(1)).zip(of_the_row) {
            if let Some(entry) = of_the_col.get_mut(at) {
                *entry = *value;
            }
        }
    }
}

/// `t - q (q' t)` for every row `t` of `rows`, in place, which is what the
/// solve against the factor leaves to be done.
///
/// `along` holds `q' t` of each row, `num_coefs` values per row, which the
/// product before this made. Every row is read and written alone and
/// nothing is summed across rows, so the values are the same whatever the
/// threads do and however many of them there are, which is what lets this
/// run on the threads of rayon; the sum over the columns of the design
/// runs in the order of the columns.
#[cfg(not(target_family = "wasm"))]
fn the_design_taken_out(
    rows: &mut [f64],
    along: &[f64],
    of_the_design: &[f64],
    num_individuals: usize,
    num_coefs: usize,
) {
    use rayon::iter::{IndexedParallelIterator, ParallelIterator};
    use rayon::slice::{ParallelSlice, ParallelSliceMut};

    rows.par_chunks_exact_mut(num_individuals.max(1))
        .zip(along.par_chunks_exact(num_coefs.max(1)))
        .for_each(|(row, along)| the_design_taken_out_of(row, along, of_the_design, num_coefs));
}

/// The same rows, one after another, which is what WebAssembly does: it
/// has no threads.
#[cfg(target_family = "wasm")]
fn the_design_taken_out(
    rows: &mut [f64],
    along: &[f64],
    of_the_design: &[f64],
    num_individuals: usize,
    num_coefs: usize,
) {
    the_design_taken_out_one_by_one(rows, along, of_the_design, num_individuals, num_coefs);
}

/// The same rows one after another: what WebAssembly does, and what the
/// test that compares the two ways of taking the design out calls.
#[cfg(any(target_family = "wasm", test))]
fn the_design_taken_out_one_by_one(
    rows: &mut [f64],
    along: &[f64],
    of_the_design: &[f64],
    num_individuals: usize,
    num_coefs: usize,
) {
    for (row, along) in rows
        .chunks_exact_mut(num_individuals.max(1))
        .zip(along.chunks_exact(num_coefs.max(1)))
    {
        the_design_taken_out_of(row, along, of_the_design, num_coefs);
    }
}

/// `t - q (q' t)` for one row, in place, with `along` the `q' t` of that
/// row.
fn the_design_taken_out_of(
    row: &mut [f64],
    along: &[f64],
    of_the_design: &[f64],
    num_coefs: usize,
) {
    for (value, of_the_individual) in row
        .iter_mut()
        .zip(of_the_design.chunks_exact(num_coefs.max(1)))
    {
        let explained = of_the_individual
            .iter()
            .zip(along)
            .map(|(of_the_coef, along)| of_the_coef * along)
            .sum::<f64>();
        *value -= explained;
    }
}

/// The factored projection against the matrix it stands for, on a
/// covariance and a design built here.
///
/// The check on real data, over the 1200 variants of both panels of
/// `docs/specs/gwas.md` and the covariance of a fit, is
/// `the_denominator_from_the_factor_is_the_one_from_the_projection_matrix`
/// of `logistic_mixed`, which has the projection matrix the fit used to
/// build to read against. What is here is the rest of the algebra: the
/// whole matrix and not one quadratic form of it, the numerator, what the
/// Wald test's subtraction gives, the largest value of the diagonal, and
/// that the threads leave the same values as one thread.
#[cfg(test)]
#[expect(
    clippy::arithmetic_side_effects,
    clippy::chunks_exact_to_as_chunks,
    reason = "every index and every count here is built from NUM_INDIVIDUALS and the \
              columns of the design, which are 12 and 2 of this fixture, so nothing \
              can overflow and nothing can be out of range, and a run of that many \
              values is a row of a matrix here as it is everywhere else in this module; \
              a test that walked its own matrices the way the code under it does would \
              not be a second opinion"
)]
mod of_a_factored_covariance {
    use popnei_linalg::{TheFirstOperand, TheSecondOperand};

    use super::{TheProjection, the_design_taken_out, the_design_taken_out_one_by_one};
    use crate::gwas::fixtures::{a_study, the_design_of};
    use crate::gwas::study::Design;

    /// How many individuals the covariance below is of: 12.
    const NUM_INDIVIDUALS: usize = 12;

    /// How many columns the design below has: 2, an intercept and one
    /// covariate.
    const NUM_COEFS: usize = 2;

    /// How far a quantity taken through the factor may be from the same
    /// quantity taken through the projection matrix, as a share of itself:
    /// 1e-12.
    ///
    /// It is the bound of the same comparison on real data in
    /// `logistic_mixed`, which measures it at 2.077e-15 on Accelerate and
    /// 2.016e-15 on faer over 2400 variants. The worst here, over the 144
    /// entries of the matrix and the four quantities of the five variants
    /// the design does not explain, is 1.627e-15 on Accelerate and
    /// 9.037e-16 on faer, measured on 25 September 2026, so this is 615
    /// times where it breaks.
    const OF_THE_TWO_ROUTES: f64 = 1e-12;

    /// A covariance of [`NUM_INDIVIDUALS`] individuals: 1 on the diagonal
    /// and a term that falls off with the distance between two of them
    /// elsewhere, which is positive definite and is not a multiple of the
    /// identity, so that the factor is not triangular by accident and the
    /// design is not taken out of nothing.
    fn the_covariance() -> Vec<f64> {
        (0..NUM_INDIVIDUALS)
            .flat_map(|row| {
                (0..NUM_INDIVIDUALS).map(move |col| match row == col {
                    true => 1.0,
                    // The distance between two individuals is at most
                    // `NUM_INDIVIDUALS`, which is 12, so it is a power
                    // every target holds.
                    false => 0.6_f64.powi(i32::try_from(row.abs_diff(col)).unwrap_or(0)),
                })
            })
            .collect()
    }

    /// The trait and the design of those individuals: an intercept and one
    /// covariate that is neither constant nor a multiple of the trait.
    fn the_trait_and_the_design() -> (Vec<f64>, Vec<f64>) {
        let phenotype: Vec<f64> = (0..NUM_INDIVIDUALS)
            .map(|at| 3.0 + (at as f64) * 0.5 - ((at % 4) as f64))
            .collect();
        let values: Vec<f64> = (0..NUM_INDIVIDUALS)
            .flat_map(|at| [1.0, ((at % 5) as f64) - 2.0])
            .collect();
        (phenotype, values)
    }

    /// Six variants of those individuals, dosages of 0, 1 and 2, of which
    /// the last is a combination of the two columns of the design and so
    /// is one the projection leaves nothing of.
    ///
    /// The five before it are distinct from each other and every one of
    /// them has variance, which the multipliers are chosen for: a variant
    /// of one dosage repeated has an `x' p x` of 0 and would compare the
    /// two routes on nothing.
    fn the_variants(values: &[f64]) -> Vec<f64> {
        let mut rows: Vec<f64> = Vec::new();
        for (var, multiplier) in [1_usize, 2, 4, 5, 7].into_iter().enumerate() {
            rows.extend(
                (0..NUM_INDIVIDUALS)
                    .map(|at| ((at.wrapping_mul(multiplier).wrapping_add(var)) % 3) as f64),
            );
        }
        rows.extend(
            values
                .chunks_exact(NUM_COEFS)
                .map(|of_the_individual| of_the_individual[0] + 0.5 * of_the_individual[1]),
        );
        rows
    }

    /// The projection matrix `p = v⁻¹ - v⁻¹ d (d' v⁻¹ d)⁻¹ d' v⁻¹` of that
    /// covariance and that design, multiplied out, which is what the
    /// factored form is read against.
    fn the_projection_matrix(design: &Design<'_>) -> Vec<f64> {
        let num_coefs = design.num_coefs();
        let mut factored = the_covariance();
        popnei_linalg::cholesky_lower(&mut factored, NUM_INDIVIDUALS).expect("the factorization");
        let mut inverse = vec![0.0_f64; NUM_INDIVIDUALS * NUM_INDIVIDUALS];
        popnei_linalg::invert_with_cholesky(&factored, NUM_INDIVIDUALS, &mut inverse)
            .expect("the inverse of the covariance");
        // The inverse comes back in the lower half alone, as
        // `the_projection_of` of `logistic_mixed` also finds it, so the
        // upper half is filled from it before anything reads a row.
        for row in 0..NUM_INDIVIDUALS {
            for col in (row + 1)..NUM_INDIVIDUALS {
                inverse[row * NUM_INDIVIDUALS + col] = inverse[col * NUM_INDIVIDUALS + row];
            }
        }
        let mut of_the_covariance = vec![0.0_f64; design.values().len()];
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &inverse,
                rows: NUM_INDIVIDUALS,
            },
            NUM_INDIVIDUALS,
            TheSecondOperand::ByTheValuesSummedOver {
                values: design.values(),
                cols: num_coefs,
            },
            &mut of_the_covariance,
        )
        .expect("the covariance's inverse with the design");
        let mut dvd = vec![0.0_f64; num_coefs * num_coefs];
        popnei_linalg::product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: design.values(),
                rows: num_coefs,
            },
            NUM_INDIVIDUALS,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &of_the_covariance,
                cols: num_coefs,
            },
            &mut dvd,
        )
        .expect("the design weighted by the covariance");
        popnei_linalg::cholesky_lower(&mut dvd, num_coefs).expect("its factorization");
        let mut solved = of_the_covariance.clone();
        popnei_linalg::solve_with_cholesky(&dvd, num_coefs, &mut solved, NUM_INDIVIDUALS)
            .expect("the solve");
        let mut explained = vec![0.0_f64; inverse.len()];
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &of_the_covariance,
                rows: NUM_INDIVIDUALS,
            },
            num_coefs,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &solved,
                cols: NUM_INDIVIDUALS,
            },
            &mut explained,
        )
        .expect("what the design explains of the inverse");
        let mut projection = inverse;
        for (entry, of_the_design) in projection.iter_mut().zip(&explained) {
            *entry -= *of_the_design;
        }
        projection
    }

    /// One vector through the projection matrix, `p x`.
    fn through_the_matrix(projection: &[f64], of_the_vector: &[f64]) -> Vec<f64> {
        projection
            .chunks_exact(NUM_INDIVIDUALS)
            .map(|row| {
                row.iter()
                    .zip(of_the_vector)
                    .map(|(entry, value)| entry * value)
                    .sum()
            })
            .collect()
    }

    /// The dot product of two vectors of the same length.
    fn the_product_of(of_the_first: &[f64], of_the_second: &[f64]) -> f64 {
        of_the_first
            .iter()
            .zip(of_the_second)
            .map(|(one, other)| one * other)
            .sum()
    }

    /// The factored projection and the matrix it stands for, built here,
    /// with the trait and the six variants taken through both.
    fn the_two_routes() -> (Vec<f64>, TheProjection, Vec<f64>, Vec<f64>, Vec<f64>) {
        let (phenotype, values) = the_trait_and_the_design();
        let individuals: Vec<usize> = (0..NUM_INDIVIDUALS).collect();
        let study = a_study(&phenotype, &values, &individuals);
        let design = the_design_of(&study, NUM_INDIVIDUALS);
        let matrix = the_projection_matrix(&design);
        let mut factored = the_covariance();
        popnei_linalg::cholesky_lower(&mut factored, NUM_INDIVIDUALS).expect("the factorization");
        let projection = TheProjection::of_the_factored_covariance(factored, &design)
            .expect("the factored projection");
        let variants = the_variants(&values);
        (matrix, projection, phenotype, values, variants)
    }

    /// Every entry of `m' m` is the entry of the projection matrix.
    ///
    /// The matrix is read entry by entry and not through one vector,
    /// because a factor that is wrong in one direction alone would survive
    /// any product with one vector. Each column of the identity goes
    /// through the factor, which gives a column of `m`, and the products
    /// of those columns are the entries of `m' m`.
    #[test]
    fn the_factor_squared_is_the_projection_matrix() {
        let (matrix, mut projection, _, _, _) = the_two_routes();
        let identity: Vec<f64> = (0..NUM_INDIVIDUALS)
            .flat_map(|row| (0..NUM_INDIVIDUALS).map(move |col| f64::from(u8::from(row == col))))
            .collect();
        let mut columns = Vec::new();
        projection
            .through(&identity, NUM_INDIVIDUALS, &mut columns)
            .expect("the identity through the factor");
        let scale = matrix
            .iter()
            .map(|entry| entry.abs())
            .fold(0.0_f64, f64::max);
        for (row, of_the_row) in columns.chunks_exact(NUM_INDIVIDUALS).enumerate() {
            for (col, of_the_col) in columns.chunks_exact(NUM_INDIVIDUALS).enumerate() {
                let from_the_factor = the_product_of(of_the_row, of_the_col);
                let from_the_matrix = matrix[row * NUM_INDIVIDUALS + col];
                let away = (from_the_factor - from_the_matrix).abs() / scale;
                assert!(
                    away <= OF_THE_TWO_ROUTES,
                    "the entry {row}, {col} of m' m is {from_the_factor} and of the \
                     projection matrix {from_the_matrix}, {away} of the largest entry \
                     {scale} apart, against the {OF_THE_TWO_ROUTES} allowed"
                );
            }
        }
    }

    /// `x' p x`, `x' p y`, `y' p y` and what the Wald test's subtraction
    /// gives come out the same through the factor and through the matrix,
    /// for six variants of which one is a combination of the design.
    ///
    /// The last of the four is the quantity **Open 2** of
    /// `docs/specs/gwas.md` is about: `y' p y - num² / den` through the
    /// matrix, `‖m y - beta m x‖²` through the factor. The variant the design
    /// explains is left out of that comparison alone, since there the two
    /// are the rounding of a cancellation and its square and agree about
    /// nothing; that it is 0 or above through the factor and of either
    /// sign through the matrix is the point of forming it, and it is what
    /// this asserts for that variant instead.
    #[test]
    fn the_four_quantities_of_a_test_are_the_same_through_both() {
        let (matrix, mut projection, phenotype, values, variants) = the_two_routes();
        let num_coefs = values.len() / NUM_INDIVIDUALS;
        assert_eq!(num_coefs, NUM_COEFS, "the columns of the design");
        let num_vars = variants.len() / NUM_INDIVIDUALS;
        assert_eq!(num_vars, 6, "the variants of the fixture");
        let mut through = Vec::new();
        projection
            .through(&variants, num_vars, &mut through)
            .expect("the variants through the factor");
        let mut of_the_trait = Vec::new();
        projection
            .through(&phenotype, 1, &mut of_the_trait)
            .expect("the trait through the factor");
        let of_the_matrix = through_the_matrix(&matrix, &phenotype);
        let ypy = the_product_of(&phenotype, &of_the_matrix);
        let of_the_factor = the_product_of(&of_the_trait, &of_the_trait);
        assert!(
            (ypy - of_the_factor).abs() <= OF_THE_TWO_ROUTES * ypy,
            "y' p y is {ypy} through the matrix and {of_the_factor} through the factor"
        );
        for (at, (of_the_variant, row)) in variants
            .chunks_exact(NUM_INDIVIDUALS)
            .zip(through.chunks_exact(NUM_INDIVIDUALS))
            .enumerate()
        {
            let through_the_matrix = through_the_matrix(&matrix, of_the_variant);
            let den = the_product_of(of_the_variant, &through_the_matrix);
            let num = the_product_of(of_the_variant, &of_the_matrix);
            let of_the_factor_den = the_product_of(row, row);
            let of_the_factor_num = the_product_of(row, &of_the_trait);
            // The variant the design explains is the last of the six, and
            // there the two routes agree about nothing. Measured on 25
            // September 2026, against a squared length of 15.25: the
            // matrix gives -4.753e-16 on Accelerate and -1.136e-15 on
            // faer, the rounding of a cancellation, of whichever sign the
            // rounding chose and here below 0 on both; the factor gives
            // 3.722e-30 and 2.169e-30, the square of that rounding, which
            // cannot fall below 0. It is what forming the quantity buys,
            // and it is what this asserts for that variant in place of the
            // four comparisons.
            if at == 5 {
                let had: f64 = of_the_variant.iter().map(|value| value * value).sum();
                assert!(
                    of_the_factor_den >= 0.0 && of_the_factor_den <= 1e-24 * had,
                    "the variant the design explains has an x' p x of {of_the_factor_den} \
                     through the factor, where its own squared length is {had}: it is 0 or \
                     above and it is the square of a rounding"
                );
                assert!(
                    den.abs() <= 1e-12 * had,
                    "the variant the design explains has an x' p x of {den} through the \
                     matrix, where its own squared length is {had}"
                );
                let beta = of_the_factor_num / of_the_factor_den;
                let left: f64 = row
                    .iter()
                    .zip(&of_the_trait)
                    .map(|(value, of_the_trait)| {
                        let left = of_the_trait - beta * value;
                        left * left
                    })
                    .sum();
                assert!(
                    left >= 0.0,
                    "what the variant the design explains leaves of the trait is {left} \
                     through the factor, and a sum of squares is 0 or above"
                );
                continue;
            }
            for (what, from_the_matrix, from_the_factor) in [
                ("x' p x", den, of_the_factor_den),
                ("x' p y", num, of_the_factor_num),
            ] {
                let away = (from_the_matrix - from_the_factor).abs() / from_the_matrix.abs();
                assert!(
                    away <= OF_THE_TWO_ROUTES,
                    "{what} of the variant {at} is {from_the_matrix} through the matrix \
                     and {from_the_factor} through the factor, {away} of itself apart"
                );
            }
            let beta = of_the_factor_num / of_the_factor_den;
            let left: f64 = row
                .iter()
                .zip(&of_the_trait)
                .map(|(value, of_the_trait)| {
                    let left = of_the_trait - beta * value;
                    left * left
                })
                .sum();
            assert!(
                left >= 0.0,
                "what the variant {at} leaves of the trait is {left} through the factor, \
                 and a sum of squares is 0 or above"
            );
            let subtracted = ypy - num * num / den;
            let away = (left - subtracted).abs() / subtracted.abs();
            assert!(
                away <= OF_THE_TWO_ROUTES,
                "what the variant {at} leaves of the trait is {subtracted} through the \
                 subtraction and {left} through the factor, {away} of itself apart"
            );
        }
    }

    /// The largest value of the diagonal of `p` that the factored form
    /// works out is the largest value of the diagonal of the matrix.
    ///
    /// It is the scale of the third row of the table of **Open 2** of
    /// `docs/specs/gwas.md`, and the factored form never builds the matrix
    /// it is the diagonal of.
    #[test]
    fn the_largest_of_the_diagonal_is_the_matrix_own() {
        let (matrix, projection, _, _, _) = the_two_routes();
        let of_the_matrix = matrix
            .chunks_exact(NUM_INDIVIDUALS)
            .zip(0..)
            .filter_map(|(row, at): (&[f64], usize)| row.get(at).copied())
            .fold(0.0_f64, f64::max);
        let found = projection.largest_of_the_diagonal();
        let away = (found - of_the_matrix).abs() / of_the_matrix;
        assert!(
            away <= OF_THE_TWO_ROUTES,
            "the largest value of the diagonal is {found} from the factor and \
             {of_the_matrix} from the matrix, {away} of itself apart"
        );
    }

    /// Taking the design out of the rows of a block on the threads of
    /// rayon leaves the same values as taking it out one row after
    /// another.
    ///
    /// Every row is read and written alone, so the two are the same bits
    /// and not the same within a tolerance, which is what a reduction
    /// across rows would cost. The block spans more rows than the machine
    /// has threads, so rayon splits it.
    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "no row reads another, so the threads make no difference to any bit \
                  and a tolerance here would hide a reduction that crossed a row"
    )]
    fn the_design_taken_out_on_the_threads_is_what_one_thread_gives() {
        let num_coefs = 3_usize;
        let num_vars = 257_usize;
        let of_the_design: Vec<f64> = (0..NUM_INDIVIDUALS * num_coefs)
            .map(|at| 0.1 + (at as f64) * 0.017)
            .collect();
        let rows: Vec<f64> = (0..num_vars * NUM_INDIVIDUALS)
            .map(|at| ((at % 7) as f64) - 3.0)
            .collect();
        let along: Vec<f64> = (0..num_vars * num_coefs)
            .map(|at| 0.5 - ((at % 11) as f64) * 0.03)
            .collect();
        let mut on_the_threads = rows.clone();
        the_design_taken_out(
            &mut on_the_threads,
            &along,
            &of_the_design,
            NUM_INDIVIDUALS,
            num_coefs,
        );
        let mut one_by_one = rows;
        the_design_taken_out_one_by_one(
            &mut one_by_one,
            &along,
            &of_the_design,
            NUM_INDIVIDUALS,
            num_coefs,
        );
        for (at, (of_the_threads, of_one)) in on_the_threads.iter().zip(&one_by_one).enumerate() {
            assert_eq!(
                of_the_threads, of_one,
                "the value {at} is {of_the_threads} on the threads and {of_one} on one"
            );
        }
    }
}
