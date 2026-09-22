//! The linear algebra of popnei.
//!
//! A calculation that reads a block of variants as a matrix, the principal
//! component analysis, the kinship, the genome wide association study,
//! needs a few operations of linear algebra, and this crate is the one
//! place that has them. It holds three: the product of a matrix with
//! itself, [`add_self_product_lower`]; the eigendecomposition of a
//! symmetric matrix, [`eigh_lower`]; and the product of two matrices,
//! [`product`]. `docs/specs/linalg.md` says what each one gives.
//!
//! Every matrix crosses this interface as a `&[f64]` held row after row,
//! the layout of the blocks and of everything the core crate holds, with
//! its numbers of rows and of columns beside it. What runs is the BLAS and
//! LAPACK of the system, the routines numpy calls, and they read a matrix
//! column after column: the buffer of an r x c matrix read that way is its
//! transpose, c x r, so each call is made on the transposes, and a caller
//! sees nothing of it.
//!
//! Those calls are the `unsafe` of popnei and they are here and nowhere
//! else. Each routine is an `unsafe fn` over slices whose lengths nothing
//! checks against the dimensions it is given, so every function of this
//! crate checks the dimensions and the values of what it was given before
//! any routine runs, and the core crate keeps its
//! `#![forbid(unsafe_code)]`.

mod blas;

use thiserror::Error as ThisError;

/// Anything that went wrong in the linear algebra.
///
/// The first two are defects of the caller, found before any routine runs,
/// and the core crate turns them into a `RuntimeError` in Python.
#[derive(Debug, ThisError)]
pub enum Error {
    /// A matrix whose buffer does not hold the values of the dimensions it
    /// was given, or a dimension that the operation does not have: a
    /// number of columns of 0, an `n` of 0, a `g` that is not `cols` x
    /// `cols`, or a matrix larger than the routines of BLAS and LAPACK
    /// take.
    #[error("the argument {argument} does not have the dimensions of the call: {expected}")]
    Dimension {
        /// The name of the argument, as "The Rust interface" of
        /// `docs/specs/linalg.md` spells it.
        argument: &'static str,
        /// What was expected of it, with the numbers of the call.
        expected: String,
    },

    /// A matrix that the operation reads holds a value that is not finite,
    /// an infinity or a NaN.
    ///
    /// It is refused here because the backends do not agree on it:
    /// `dsyevd` on a matrix with a NaN gives NaN eigenvalues and says that
    /// nothing went wrong, and faer's `self_adjoint_eigen` gives its error
    /// of no convergence.
    #[error("the matrix {argument} holds a value that is not finite")]
    NotFinite {
        /// The name of the argument that holds it.
        argument: &'static str,
    },

    /// A decomposition that did not converge.
    #[error("{routine} did not converge: it gave the info {info}")]
    NoConvergence {
        /// The routine of LAPACK that stopped, or `faer`.
        routine: &'static str,
        /// The `info` the routine gave, which is not 0.
        info: i32,
    },
}

/// What every operation of this crate gives back.
pub type Result<T> = std::result::Result<T, Error>;

/// The eigenvalues and the eigenvectors of a symmetric matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct Eigen {
    /// The `n` eigenvalues, from the largest.
    pub values: Vec<f64>,
    /// n x n, row after row: row `j` is the eigenvector of `values[j]`,
    /// of length 1. The sign of a vector is whatever the backend gave, and
    /// a caller that needs a fixed sign fixes it.
    pub vectors: Vec<f64>,
}

/// Adds the lower half of `a'a` to `g`: `g += a'a`, where `a` is `rows` x
/// `cols` and `g` is `cols` x `cols`, both row after row.
///
/// The entry `i`, `j` of `a'a` is the sum over the rows of `a` of the
/// value in column `i` times the value in column `j`. Only the lower half
/// of `g` is read and written, the entries of column `j` at most `i` of
/// row `i`; the upper half is left as it was, and a caller that needs the
/// whole matrix mirrors it.
///
/// `rows` may be 0, and then `g` is left as it was: a block of variants
/// whose rows all had no variance gives such an `a` and adds nothing.
///
/// # Errors
///
/// [`Error::Dimension`] when `cols` is 0, when `a` holds fewer than `rows`
/// times `cols` values or `g` fewer than `cols` times `cols`, or when a
/// dimension is larger than the routines of BLAS and LAPACK take.
/// [`Error::NotFinite`] when `a`, or the lower half of `g`, holds a value
/// that is not finite.
pub fn add_self_product_lower(a: &[f64], rows: usize, cols: usize, g: &mut [f64]) -> Result<()> {
    if cols == 0 {
        return Err(Error::Dimension {
            argument: "cols",
            expected: "1 at least, since the product is the cols x cols matrix g".to_owned(),
        });
    }
    let a = the_matrix_of(a, rows, cols, "a")?;
    let g = the_matrix_of_mut(g, cols, cols, "g")?;
    refuse_a_value_that_is_not_finite(a, "a")?;
    refuse_a_value_that_is_not_finite_in_the_lower_half(g, cols, "g")?;
    if rows == 0 {
        return Ok(());
    }
    crate::blas::add_self_product_lower(a, rows, cols, g)
}

/// The product `c = a b`, where `a` is `rows` x `inner`, `b` is `inner` x
/// `cols` and `c`, which is overwritten, is `rows` x `cols`, all row after
/// row.
///
/// `rows` may be 0, and then `c` holds no rows and nothing is written.
/// `inner` and `cols` are 1 at least: a product with no inner dimension or
/// no column is a defect of the caller, as `docs/specs/linalg.md` has a
/// `cols` of 0 among its errors.
///
/// # Errors
///
/// [`Error::Dimension`] when `inner` or `cols` is 0, when a buffer holds
/// fewer values than its rows times its columns, or when a dimension is
/// larger than the routines of BLAS and LAPACK take. [`Error::NotFinite`]
/// when `a` or `b` holds a value that is not finite.
pub fn product(
    a: &[f64],
    rows: usize,
    inner: usize,
    b: &[f64],
    cols: usize,
    c: &mut [f64],
) -> Result<()> {
    if inner == 0 {
        return Err(Error::Dimension {
            argument: "inner",
            expected: "1 at least, since it is the dimension the product sums over".to_owned(),
        });
    }
    if cols == 0 {
        return Err(Error::Dimension {
            argument: "cols",
            expected: "1 at least, since it is the number of columns of the product".to_owned(),
        });
    }
    let a = the_matrix_of(a, rows, inner, "a")?;
    let b = the_matrix_of(b, inner, cols, "b")?;
    let c = the_matrix_of_mut(c, rows, cols, "c")?;
    refuse_a_value_that_is_not_finite(a, "a")?;
    refuse_a_value_that_is_not_finite(b, "b")?;
    if rows == 0 {
        return Ok(());
    }
    crate::blas::product(a, rows, inner, b, cols, c)
}

/// The eigendecomposition of the symmetric `g` of `n` x `n`, given by its
/// lower half.
///
/// The `n` eigenvalues come from the largest, and row `j` of the vectors
/// is the eigenvector of `values[j]`, of length 1, so that `g v = λ v` and
/// the vectors are at right angles to each other. Only the lower half of
/// `g` is read, the entries of column `j` at most `i` of row `i`; what the
/// upper half holds does not reach the result.
///
/// The buffer of `g` comes back as the eigenvectors, with no copy, which
/// is why `g` is taken by value: a caller that needs `g` afterwards copies
/// it first.
///
/// # Errors
///
/// [`Error::Dimension`] when `n` is 0, when `g` holds fewer than `n` times
/// `n` values, or when `n` is larger than the routines of LAPACK take.
/// [`Error::NotFinite`] when the lower half of `g` holds a value that is
/// not finite. [`Error::NoConvergence`] when the routine stopped.
pub fn eigh_lower(mut g: Vec<f64>, n: usize) -> Result<Eigen> {
    if n == 0 {
        return Err(Error::Dimension {
            argument: "n",
            expected: "1 at least, since g is the n x n matrix to decompose".to_owned(),
        });
    }
    let values = how_many(n, n, "g")?;
    if g.len() < values {
        return Err(Error::Dimension {
            argument: "g",
            expected: format!(
                "{values} values, {n} rows times {n} columns, and it holds {held}",
                held = g.len()
            ),
        });
    }
    g.truncate(values);
    refuse_a_value_that_is_not_finite_in_the_lower_half(&g, n, "g")?;
    crate::blas::eigh_lower(g, n)
}

/// The values of a matrix of `rows` x `cols`.
///
/// # Errors
///
/// [`Error::Dimension`] when they are more values than this machine can
/// hold.
fn how_many(rows: usize, cols: usize, argument: &'static str) -> Result<usize> {
    rows.checked_mul(cols).ok_or_else(|| Error::Dimension {
        argument,
        expected: format!(
            "{rows} rows times {cols} columns, which is more values than this machine can hold"
        ),
    })
}

/// The first `rows` times `cols` values of the buffer, which are the
/// matrix held row after row.
///
/// # Errors
///
/// [`Error::Dimension`] when the buffer holds fewer.
fn the_matrix_of<'a>(
    values: &'a [f64],
    rows: usize,
    cols: usize,
    argument: &'static str,
) -> Result<&'a [f64]> {
    let wanted = how_many(rows, cols, argument)?;
    values.get(..wanted).ok_or_else(|| Error::Dimension {
        argument,
        expected: format!(
            "{wanted} values, {rows} rows times {cols} columns, and it holds {held}",
            held = values.len()
        ),
    })
}

/// The same for a matrix the operation writes.
///
/// # Errors
///
/// [`Error::Dimension`] when the buffer holds fewer.
fn the_matrix_of_mut<'a>(
    values: &'a mut [f64],
    rows: usize,
    cols: usize,
    argument: &'static str,
) -> Result<&'a mut [f64]> {
    let wanted = how_many(rows, cols, argument)?;
    let held = values.len();
    values.get_mut(..wanted).ok_or_else(|| Error::Dimension {
        argument,
        expected: format!("{wanted} values, {rows} rows times {cols} columns, and it holds {held}"),
    })
}

/// # Errors
///
/// [`Error::NotFinite`] when one of the values is an infinity or a NaN.
fn refuse_a_value_that_is_not_finite(values: &[f64], argument: &'static str) -> Result<()> {
    if values.iter().all(|value| value.is_finite()) {
        Ok(())
    } else {
        Err(Error::NotFinite { argument })
    }
}

/// The same for the lower half alone of an `n` x `n` matrix, the entries
/// of column `j` at most `i` of row `i`, which is what the routines read
/// of a matrix given by its lower half. `n` is 1 at least.
///
/// # Errors
///
/// [`Error::NotFinite`] when one of those values is an infinity or a NaN.
fn refuse_a_value_that_is_not_finite_in_the_lower_half(
    values: &[f64],
    n: usize,
    argument: &'static str,
) -> Result<()> {
    let it_is_all_finite = values.chunks_exact(n).enumerate().all(|(row, entries)| {
        entries
            .get(..row.saturating_add(1))
            .is_some_and(|half| half.iter().all(|value| value.is_finite()))
    });
    if it_is_all_finite {
        Ok(())
    } else {
        Err(Error::NotFinite { argument })
    }
}

#[cfg(test)]
mod tests {
    use super::{Eigen, Error, add_self_product_lower, eigh_lower, product};

    /// The A of 2 x 3 of "How it is verified" of `docs/specs/linalg.md`,
    /// rows (1, 2, 0) and (0, 1, 3), row after row.
    const A_OF_2_BY_3: [f64; 6] = [1.0, 2.0, 0.0, 0.0, 1.0, 3.0];

    /// The B of 3 x 2 of the same place, rows (1, 0), (2, 1) and (0, 3).
    const B_OF_3_BY_2: [f64; 6] = [1.0, 0.0, 2.0, 1.0, 0.0, 3.0];

    /// The 3 x 3 symmetric matrix of "How it is verified", rows (4, 1, 0),
    /// (1, 3, 0) and (0, 0, 1), with only its lower half given and the
    /// upper half holding a value that is nothing of the matrix, so that a
    /// call that read the upper half instead would give other eigenvalues.
    fn the_matrix_of_3_by_3() -> Vec<f64> {
        vec![
            4.0, 99.0, 99.0, //
            1.0, 3.0, 99.0, //
            0.0, 0.0, 1.0,
        ]
    }

    /// The vector with the sign that makes its entry of largest absolute
    /// value positive, which is the sign the spec fixes before comparing,
    /// since the sign a backend gives is its own.
    fn with_the_sign_of_the_spec(vector: &[f64]) -> Vec<f64> {
        let largest = vector.iter().copied().fold(0.0_f64, |so_far, value| {
            if value.abs() > so_far.abs() {
                value
            } else {
                so_far
            }
        });
        let sign = if largest < 0.0 { -1.0 } else { 1.0 };
        vector.iter().map(|value| value * sign).collect()
    }

    /// The entries of the two differ by more than `tolerance`, or the two
    /// do not hold the same number of entries.
    fn differ(got: &[f64], expected: &[f64], tolerance: f64) -> bool {
        got.len() != expected.len()
            || got
                .iter()
                .zip(expected)
                .any(|(one, other)| (one - other).abs() > tolerance)
    }

    #[test]
    fn the_self_product_writes_the_lower_half_and_leaves_the_upper_as_it_was() {
        // A'A is (1, 2, 0), (2, 5, 3), (0, 3, 9), and every entry is a sum
        // of at most two products of small integers, so it is exact.
        let mut g = vec![0.0, 7.0, 7.0, 0.0, 0.0, 7.0, 0.0, 0.0, 0.0];
        add_self_product_lower(&A_OF_2_BY_3, 2, 3, &mut g).unwrap();
        assert_eq!(
            g,
            vec![
                1.0, 7.0, 7.0, //
                2.0, 5.0, 7.0, //
                0.0, 3.0, 9.0,
            ]
        );
    }

    #[test]
    fn the_self_product_adds_to_what_g_held() {
        let mut g = vec![0.0; 9];
        add_self_product_lower(&A_OF_2_BY_3, 2, 3, &mut g).unwrap();
        add_self_product_lower(&A_OF_2_BY_3, 2, 3, &mut g).unwrap();
        assert_eq!(
            g,
            vec![
                2.0, 0.0, 0.0, //
                4.0, 10.0, 0.0, //
                0.0, 6.0, 18.0,
            ]
        );
    }

    #[test]
    fn the_self_product_of_an_a_of_no_rows_leaves_g_as_it_was() {
        let mut g = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        add_self_product_lower(&[], 0, 3, &mut g).unwrap();
        assert_eq!(g, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    }

    #[test]
    fn the_product_of_two_matrices_that_are_not_square() {
        // The 2 x 2 matrix with rows (5, 2) and (2, 10), exactly.
        let mut c = vec![0.0; 4];
        product(&A_OF_2_BY_3, 2, 3, &B_OF_3_BY_2, 2, &mut c).unwrap();
        assert_eq!(c, vec![5.0, 2.0, 2.0, 10.0]);
    }

    #[test]
    fn the_product_of_an_a_of_no_rows_writes_nothing() {
        let mut c: Vec<f64> = Vec::new();
        product(&[], 0, 3, &B_OF_3_BY_2, 2, &mut c).unwrap();
        assert!(c.is_empty());
    }

    #[test]
    fn the_eigendecomposition_gives_the_values_from_the_largest_and_the_vectors_as_rows() {
        let Eigen { values, vectors } = eigh_lower(the_matrix_of_3_by_3(), 3).unwrap();
        // (7 + √5) / 2, (7 - √5) / 2 and 1, the twelve digits of the spec.
        assert!(
            !differ(&values, &[4.618_033_988_750, 2.381_966_011_250, 1.0], 1e-12),
            "the eigenvalues are {values:?}"
        );
        let expected = [
            [0.850_650_808_352, 0.525_731_112_119, 0.0],
            [-0.525_731_112_119, 0.850_650_808_352, 0.0],
            [0.0, 0.0, 1.0],
        ];
        for (row, expected_vector) in vectors.as_chunks::<3>().0.iter().zip(&expected) {
            let got = with_the_sign_of_the_spec(row);
            assert!(
                !differ(&got, expected_vector, 1e-12),
                "the eigenvector is {got:?} and not {expected_vector:?}"
            );
        }
        assert_eq!(vectors.len(), 9);
    }

    #[test]
    fn the_self_product_refuses_a_g_that_is_not_cols_by_cols() {
        let mut g = vec![0.0; 8];
        let error = add_self_product_lower(&A_OF_2_BY_3, 2, 3, &mut g).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "g", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_self_product_refuses_an_a_shorter_than_its_rows_times_its_columns() {
        let mut g = vec![0.0; 9];
        let error = add_self_product_lower(&A_OF_2_BY_3[..5], 2, 3, &mut g).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_self_product_refuses_a_cols_of_zero() {
        let mut g: Vec<f64> = Vec::new();
        let error = add_self_product_lower(&[], 2, 0, &mut g).unwrap_err();
        assert!(
            matches!(
                error,
                Error::Dimension {
                    argument: "cols",
                    ..
                }
            ),
            "the error is {error}"
        );
    }

    #[test]
    fn the_self_product_refuses_an_a_with_a_value_that_is_not_finite() {
        let mut g = vec![0.0; 9];
        let a = [1.0, 2.0, 0.0, 0.0, f64::NAN, 3.0];
        let error = add_self_product_lower(&a, 2, 3, &mut g).unwrap_err();
        assert!(
            matches!(error, Error::NotFinite { argument: "a" }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_self_product_refuses_a_g_whose_lower_half_holds_a_value_that_is_not_finite() {
        // The value that is not finite is at row 1, column 0, which the
        // product reads; the same value in the upper half is not read and
        // is not an error, which the next lines check.
        let mut g = vec![0.0; 9];
        *g.get_mut(3).unwrap() = f64::INFINITY;
        let error = add_self_product_lower(&A_OF_2_BY_3, 2, 3, &mut g).unwrap_err();
        assert!(
            matches!(error, Error::NotFinite { argument: "g" }),
            "the error is {error}"
        );
        let mut g = vec![0.0; 9];
        *g.get_mut(1).unwrap() = f64::INFINITY;
        add_self_product_lower(&A_OF_2_BY_3, 2, 3, &mut g).unwrap();
    }

    #[test]
    fn the_product_refuses_a_b_that_does_not_have_the_inner_dimension() {
        let mut c = vec![0.0; 4];
        let error = product(&A_OF_2_BY_3, 2, 3, &B_OF_3_BY_2[..4], 2, &mut c).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "b", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_product_refuses_a_c_shorter_than_its_rows_times_its_columns() {
        let mut c = vec![0.0; 3];
        let error = product(&A_OF_2_BY_3, 2, 3, &B_OF_3_BY_2, 2, &mut c).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "c", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_product_refuses_a_cols_of_zero_and_an_inner_of_zero() {
        let mut c: Vec<f64> = Vec::new();
        let error = product(&A_OF_2_BY_3, 2, 3, &B_OF_3_BY_2, 0, &mut c).unwrap_err();
        assert!(
            matches!(
                error,
                Error::Dimension {
                    argument: "cols",
                    ..
                }
            ),
            "the error is {error}"
        );
        let error = product(&[], 2, 0, &[], 2, &mut c).unwrap_err();
        assert!(
            matches!(
                error,
                Error::Dimension {
                    argument: "inner",
                    ..
                }
            ),
            "the error is {error}"
        );
    }

    #[test]
    fn the_product_refuses_a_value_that_is_not_finite() {
        let mut c = vec![0.0; 4];
        let b = [1.0, 0.0, 2.0, f64::NEG_INFINITY, 0.0, 3.0];
        let error = product(&A_OF_2_BY_3, 2, 3, &b, 2, &mut c).unwrap_err();
        assert!(
            matches!(error, Error::NotFinite { argument: "b" }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_eigendecomposition_refuses_an_n_of_zero() {
        let error = eigh_lower(Vec::new(), 0).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "n", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_eigendecomposition_refuses_a_g_shorter_than_n_times_n() {
        let error = eigh_lower(vec![0.0; 8], 3).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "g", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_eigendecomposition_refuses_a_value_that_is_not_finite_in_the_lower_half() {
        let mut g = the_matrix_of_3_by_3();
        *g.get_mut(3).unwrap() = f64::NAN;
        let error = eigh_lower(g, 3).unwrap_err();
        assert!(
            matches!(error, Error::NotFinite { argument: "g" }),
            "the error is {error}"
        );
    }

    #[test]
    fn a_decomposition_that_did_not_converge_names_the_routine_and_the_info() {
        let error = Error::NoConvergence {
            routine: "dsyevd",
            info: 7,
        };
        assert_eq!(
            error.to_string(),
            "dsyevd did not converge: it gave the info 7"
        );
    }
}
