//! The linear algebra of popnei.
//!
//! A calculation that reads a block of variants as a matrix, the principal
//! component analysis, the kinship, the genome wide association study,
//! needs a few operations of linear algebra, and this crate is the one
//! place that has them. It holds fourteen: the product of a matrix with
//! itself, [`add_self_product_lower`]; the eigendecomposition of a
//! symmetric matrix, [`eigh_lower`]; the Cholesky factorization of a
//! symmetric positive definite one, [`cholesky_lower`], the solve of a
//! system with the matrix it factored, [`solve_with_cholesky`], the log of
//! that matrix's determinant, [`log_determinant_with_cholesky`], and its
//! inverse, [`invert_with_cholesky`]; the thin QR factorization of a
//! design, [`thin_qr`], the solve of a system against a triangular
//! matrix, [`solve_triangular`], which is two of the fourteen because
//! [`TheHalfThatHoldsTheMatrix`] says whether the upper half holds it, as
//! the `r` of that factorization does, or the lower one, as a Cholesky
//! factor does; and how many of a matrix's columns are independent,
//! [`rank`]; and the product of two matrices, [`product`], which is
//! another four, because [`TheFirstOperand`] and [`TheSecondOperand`] each
//! say how one matrix's buffer is laid out and the two together choose
//! among `a b`, `a b'`, `a' b` and `a' b'`.
//! `docs/specs/linalg.md` says what each one gives.
//!
//! Every matrix crosses this interface as a `&[f64]` held row after row,
//! the layout of the blocks and of everything the core crate holds, with
//! its numbers of rows and of columns beside it.
//!
//! Two backends run under that interface and give the same numbers within
//! the tolerance of "How it is verified" of the spec. One is the BLAS and
//! LAPACK of the system, the library numpy is built on too, which is what
//! runs when the target is not WebAssembly and the cargo feature `blas` is
//! on, as it is by default. Which routines of it popnei calls and which of
//! those numpy calls is in `blas.rs`. The other is faer, a linear algebra library written in
//! Rust, which runs on both wasm targets, where there is no BLAS, and
//! natively when the crate is built with `--no-default-features`, so that
//! a machine with no BLAS and no Fortran compiler builds popnei. A caller
//! does not know which one ran.
//!
//! The routines of BLAS and LAPACK read a matrix column after column: the
//! buffer of an r x c matrix read that way is its transpose, c x r, so
//! that backend calls each routine on the transposes. faer is told that
//! the buffers are row major. Neither shows in what a caller gets.
//!
//! The calls to BLAS and LAPACK are the `unsafe` of popnei and they are
//! here and nowhere else. Each routine is an `unsafe fn` over slices whose
//! lengths nothing checks against the dimensions it is given, so every
//! function of this crate checks the dimensions and the values of what it
//! was given before any routine runs, whichever backend is to run it, and
//! the core crate keeps its `#![forbid(unsafe_code)]`.

// The backend is BLAS and LAPACK when the target is not WebAssembly and
// the cargo feature `blas` is on, which it is by default, and faer
// otherwise: on both wasm targets, where there is no BLAS, and natively
// when the crate is built with `--no-default-features`. Each file holds
// the calls of one library and nothing else, and the two are the same
// module to the rest of the crate.
#[cfg(all(feature = "blas", not(target_family = "wasm")))]
#[path = "blas.rs"]
mod backend;
#[cfg(not(all(feature = "blas", not(target_family = "wasm"))))]
#[path = "faer.rs"]
mod backend;

use std::fmt;

use thiserror::Error as ThisError;

/// The largest number of values a matrix of this crate holds, and so the
/// largest any of its dimensions may be.
///
/// The routines of BLAS and LAPACK take every dimension as an `i32` and
/// index their matrices in it, and faer is held to the same limit so that
/// the two backends refuse the same calls. A square matrix of this many
/// values is 46340 x 46340, which is 17 GB of `f64`.
pub const THE_MOST_VALUES_OF_A_MATRIX: usize = 2_147_483_647;

/// Anything that went wrong in the linear algebra.
///
/// The first two are defects of the caller, found before any routine runs,
/// and the core crate turns them into a `RuntimeError` in Python.
#[derive(Debug, ThisError)]
pub enum Error {
    /// A matrix whose buffer does not hold the values of the dimensions it
    /// was given, or a dimension that the operation does not have: a
    /// number of columns of 0, an `inner` of 0, an `n` of 0, a `g` that is
    /// not `cols` x `cols`, or a matrix of more than 2147483647 values,
    /// which is what the routines of BLAS and LAPACK count in.
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

    /// A routine that stopped: it did not converge, or it refused an
    /// argument it was given, which is a defect of popnei.
    #[error(fmt = the_message_of_a_routine_that_stopped)]
    NoConvergence {
        /// The routine of LAPACK that stopped, or `faer`.
        routine: &'static str,
        /// The `info` the routine of LAPACK gave: above 0 when it did not
        /// converge, and below 0 when it refused the argument of that
        /// number, counting from 1. It is 0 when the routine is faer,
        /// which says only that it did not converge.
        info: i32,
    },

    /// A workspace that a routine needs and that this machine has not the
    /// memory for.
    ///
    /// Two operations ask for one, and each backend for its own. The
    /// eigendecomposition on BLAS needs 1 + 6n + 2n² floats besides the
    /// matrix, and a workspace of integers, which is 16 MB at n = 1000 and
    /// 1.6 GB at n = 10000; the inverse on faer needs a scratch of n x n,
    /// 8 MB at n = 1000 and 800 MB at n = 10000. So the crate asks for
    /// that memory instead of taking it, and a machine that has not got it
    /// gets this error where it would otherwise see the process end. No
    /// test of popnei reaches either case.
    #[error("this machine has not the memory for {what}, {values} values")]
    Memory {
        /// What could not be allocated.
        what: &'static str,
        /// How many values it holds.
        values: usize,
    },

    /// A matrix that could not be factored at the row the value names,
    /// counting from 0: the Cholesky reached a diagonal entry that is not
    /// above 0 there, or the solve against a triangular matrix reached one
    /// that is 0, in whichever half the caller named.
    ///
    /// The three operations that are given an `l` that a Cholesky made,
    /// the solve, the log of the determinant and the inverse, read its
    /// diagonal for an entry that is not above 0 and give this same error
    /// at the same row, which is the row the [`cholesky_lower`] that would
    /// have produced that `l` stops at.
    ///
    /// It is not a defect of the caller: the call was right and the matrix
    /// was what the data made it, so the module that called decides what
    /// it means there. Both backends give the same row, counted from 0:
    /// `dpotrf` gives the order of the leading corner, counting from 1,
    /// and faer an index from 0.
    /// The message says where it failed and not why: five operations give
    /// this error and only [`cholesky_lower`] factors anything, so a
    /// message naming a factorization, or calling the matrix singular,
    /// would be wrong for the other four, where the matrix was handed in
    /// as a factorization or, for [`solve_triangular`], was never one.
    /// Drawing no conclusion about the caller's data is what lets a
    /// caller whose matrix was built from the user's, and not given by
    /// them, say what it means itself.
    #[error("the matrix {argument} failed at its row {at}, counting from 0")]
    Singular {
        /// The name of the argument, as "The Rust interface" of
        /// `docs/specs/linalg.md` spells it.
        argument: &'static str,
        /// The row the factorization stopped at, counting from 0.
        at: usize,
    },
}

/// The message of [`Error::NoConvergence`], which says three different
/// things by the sign of the `info` the routine gave.
fn the_message_of_a_routine_that_stopped(
    routine: &&'static str,
    info: &i32,
    formatter: &mut fmt::Formatter<'_>,
) -> fmt::Result {
    if *info > 0 {
        write!(
            formatter,
            "{routine} did not converge: it gave the info {info}"
        )
    } else if *info < 0 {
        write!(
            formatter,
            "{routine} refused its argument {argument}, a defect of popnei",
            argument = info.unsigned_abs()
        )
    } else {
        write!(formatter, "{routine} did not converge")
    }
}

/// What every operation of this crate gives back.
pub type Result<T> = std::result::Result<T, Error>;

/// The eigenvalues and the eigenvectors of a symmetric matrix.
#[derive(Debug, Clone)]
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
/// whose rows all had no variance gives such an `a` and adds nothing. `a`
/// may hold more values than `rows` times `cols`, and then its first
/// `rows` times `cols` are the matrix; `g` holds exactly `cols` times
/// `cols`.
///
/// # Errors
///
/// [`Error::Dimension`] when `cols` is 0, when `a` holds fewer than `rows`
/// times `cols` values, when `g` does not hold exactly `cols` times
/// `cols`, or when a matrix would hold more than 2147483647 values, which
/// is what the routines of BLAS and LAPACK count in.
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
    refuse_a_g_that_is_not_square(g.len(), cols)?;
    refuse_a_value_that_is_not_finite(a, "a")?;
    refuse_a_value_that_is_not_finite_in_the_lower_half(g, cols, "g")?;
    if rows == 0 {
        return Ok(());
    }
    backend::add_self_product_lower(a, rows, cols, g)
}

/// The first operand of a [`product`], with the way it holds its values.
///
/// A product sums over `inner` values for each row of the result, and the
/// first operand holds them either as one row for each row of the result,
/// which is the first matrix of a product as it is usually written, or as
/// one row for each of the values summed over, which is what a caller has
/// whose two matrices are both laid out by the thing they describe and
/// which gives `c = a' b`. Which of the two it is belongs here and not in
/// the name of a function, for the reason [`TheSecondOperand`] gives: the
/// two carry the same numbers, so a call that named the wrong one of two
/// functions would pass every check and give a wrong matrix with no error.
#[derive(Debug, Clone, Copy)]
pub enum TheFirstOperand<'a> {
    /// `rows` x `inner`, one row for each row of the result, which is the
    /// first matrix of a product as it is usually written.
    ByTheRowsOfTheResult {
        /// The values of the matrix, row after row.
        values: &'a [f64],
        /// How many rows it has, which is how many the result has.
        rows: usize,
    },
    /// `inner` x `rows`, one row for each of the values the product sums
    /// over, which gives `c = a' b`.
    ByTheValuesSummedOver {
        /// The values of the matrix, row after row.
        values: &'a [f64],
        /// How many columns it has, which is how many rows the result
        /// has.
        rows: usize,
    },
}

impl TheFirstOperand<'_> {
    /// How many rows the result has, which the operand holds as its rows
    /// or as its columns by the way it is laid out.
    fn rows(self) -> usize {
        match self {
            TheFirstOperand::ByTheRowsOfTheResult { rows, .. }
            | TheFirstOperand::ByTheValuesSummedOver { rows, .. } => rows,
        }
    }
}

/// The second operand of a [`product`], with the way it holds its values.
///
/// A product sums over the `inner` values of each row of its first
/// operand, and the second operand holds those values either as one row
/// for each of them, which is the second matrix of a product as it is
/// usually written, or as one row for each column of the result, which is
/// what a caller has whose two matrices are both laid out by the thing
/// they describe. Which of the two it is belongs here and not in the name
/// of a function, because the two take the same numbers and give
/// different matrices: a call that named the wrong one of two functions
/// would pass every check, since a buffer of `inner` x `cols` values
/// holds `cols` x `inner` as well, and would give a wrong matrix with no
/// error.
#[derive(Debug, Clone, Copy)]
pub enum TheSecondOperand<'a> {
    /// `inner` x `cols`, one row for each of the values the product sums
    /// over, which gives `c = a b`.
    ByTheValuesSummedOver {
        /// The values of the matrix, row after row.
        values: &'a [f64],
        /// How many columns it has, which is how many the result has.
        cols: usize,
    },
    /// `cols` x `inner`, one row for each column of the result, which
    /// gives `c = a b'`. Both backends read it that way inside the
    /// routine and neither pays a copy for it.
    ByTheColumnsOfTheResult {
        /// The values of the matrix, row after row.
        values: &'a [f64],
        /// How many rows it has, which is how many columns the result
        /// has.
        cols: usize,
    },
}

impl TheSecondOperand<'_> {
    /// How many columns the result has, which the operand holds as its
    /// columns or as its rows by the way it is laid out.
    fn cols(self) -> usize {
        match self {
            TheSecondOperand::ByTheValuesSummedOver { cols, .. }
            | TheSecondOperand::ByTheColumnsOfTheResult { cols, .. } => cols,
        }
    }
}

/// The product of `a` with `b` into `c`, which is overwritten and is
/// `rows` x `cols`, all row after row.
///
/// [`TheFirstOperand`] says whether `a` is `rows` x `inner` or `inner` x
/// `rows`, and [`TheSecondOperand`] whether `b` is `inner` x `cols` or
/// `cols` x `inner`, and the two together say which of `a b`, `a b'`,
/// `a' b` and `a' b'` is written, as "The product with its first operand
/// turned" of `docs/specs/linalg.md` lays out. The entry i, j of the
/// result is the sum over the `inner` values of the row i, or the column
/// i, of `a` times the column j, or the row j, of `b`.
///
/// Giving one slice for `a` as
/// [`TheFirstOperand::ByTheRowsOfTheResult`] and for `b` as
/// [`TheSecondOperand::ByTheColumnsOfTheResult`], with `rows` equal to
/// `cols`, is the product of a matrix with its own transpose.
///
/// `rows` may be 0, and then `c` holds no rows and nothing is written.
/// `inner` and `cols` are 1 at least: a product with no inner dimension
/// or no column is a defect of the caller, as `docs/specs/linalg.md` has
/// a `cols` of 0 among its errors. A buffer may hold more values than the
/// rows times the columns of the way it is read, and then its first rows
/// times columns are the matrix.
///
/// # Errors
///
/// [`Error::Dimension`] when `inner` or `cols` is 0, when a buffer holds
/// fewer values than the rows times the columns of the way it is read, or
/// when a matrix would hold more than 2147483647 values, which is what
/// the routines of BLAS and LAPACK count in. [`Error::NotFinite`] when
/// `a` or `b` holds a value that is not finite.
///
/// What these do not catch: `rows`, `inner` and `cols` are the least a
/// buffer may hold and not what it holds, so a call that names fewer
/// rows than its matrix has, or that names the wrong case of
/// [`TheFirstOperand`] or of [`TheSecondOperand`], holds enough values
/// either way and is answered with another matrix and no error. The
/// names of the cases are what a caller has instead of a check.
pub fn product(
    a: TheFirstOperand<'_>,
    inner: usize,
    b: TheSecondOperand<'_>,
    c: &mut [f64],
) -> Result<()> {
    if inner == 0 {
        return Err(Error::Dimension {
            argument: "inner",
            expected: "1 at least, since it is the dimension the product sums over".to_owned(),
        });
    }
    let rows = a.rows();
    let cols = b.cols();
    if cols == 0 {
        return Err(Error::Dimension {
            argument: "cols",
            expected: "1 at least, since it is the number of columns of the product".to_owned(),
        });
    }
    // Each buffer is checked against the dimensions of the way it is read,
    // which for the first operand is rows x inner one way round and inner
    // x rows the other. The two hold the same number of values, so what
    // the way it is read changes is what the error says was expected.
    let values_of_a = match a {
        TheFirstOperand::ByTheRowsOfTheResult { values, rows } => {
            the_matrix_of(values, rows, inner, "a")?
        }
        TheFirstOperand::ByTheValuesSummedOver { values, rows } => {
            the_matrix_of(values, inner, rows, "a")?
        }
    };
    let values_of_b = match b {
        TheSecondOperand::ByTheValuesSummedOver { values, cols } => {
            the_matrix_of(values, inner, cols, "b")?
        }
        TheSecondOperand::ByTheColumnsOfTheResult { values, cols } => {
            the_matrix_of(values, cols, inner, "b")?
        }
    };
    let c = the_matrix_of_mut(c, rows, cols, "c")?;
    refuse_a_value_that_is_not_finite(values_of_a, "a")?;
    refuse_a_value_that_is_not_finite(values_of_b, "b")?;
    if rows == 0 {
        // No row of the result to write. Both backends write nothing for a
        // product of no rows anyway, so no test tells this line from its
        // absence; what it buys is that no routine of a backend is ever
        // reached with a dimension of 0 and the pointer of an empty slice.
        return Ok(());
    }
    match (a, b) {
        (
            TheFirstOperand::ByTheRowsOfTheResult { .. },
            TheSecondOperand::ByTheValuesSummedOver { .. },
        ) => backend::product(values_of_a, rows, inner, values_of_b, cols, c),
        (
            TheFirstOperand::ByTheRowsOfTheResult { .. },
            TheSecondOperand::ByTheColumnsOfTheResult { .. },
        ) => {
            backend::product_with_the_second_turned(values_of_a, rows, inner, values_of_b, cols, c)
        }
        (
            TheFirstOperand::ByTheValuesSummedOver { .. },
            TheSecondOperand::ByTheValuesSummedOver { .. },
        ) => backend::product_with_the_first_turned(values_of_a, rows, inner, values_of_b, cols, c),
        (
            TheFirstOperand::ByTheValuesSummedOver { .. },
            TheSecondOperand::ByTheColumnsOfTheResult { .. },
        ) => backend::product_with_both_turned(values_of_a, rows, inner, values_of_b, cols, c),
    }
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
/// The buffer of `g` comes back as the eigenvectors, which is why `g` is
/// taken by value: a caller that needs `g` afterwards copies it first. It
/// holds exactly `n` times `n` values.
///
/// # Errors
///
/// [`Error::Dimension`] when `n` is 0, when `g` does not hold exactly `n`
/// times `n` values, or when `n` times `n` is more than 2147483647, which
/// is what the routines of BLAS and LAPACK count in.
/// [`Error::NotFinite`] when the lower half of `g` holds a value that is
/// not finite.
/// [`Error::NoConvergence`] when the routine stopped. [`Error::Memory`]
/// when the workspace the routine needs could not be allocated.
pub fn eigh_lower(g: Vec<f64>, n: usize) -> Result<Eigen> {
    if n == 0 {
        return Err(Error::Dimension {
            argument: "n",
            expected: "1 at least, since g is the n x n matrix to decompose".to_owned(),
        });
    }
    refuse_a_g_that_is_not_square(g.len(), n)?;
    refuse_a_value_that_is_not_finite_in_the_lower_half(&g, n, "g")?;
    // Both backends give the eigenvalues from the smallest, LAPACK and
    // faer alike, and both give each eigenvector where popnei's row major
    // buffer reads it as a row. So the turning round is done here, once,
    // and neither backend does it.
    let mut eigen = backend::eigh_lower(g, n)?;
    eigen.values.reverse();
    reverse_the_rows(&mut eigen.vectors, n);
    Ok(eigen)
}

/// Swaps the first row of the matrix with the last, the second with the
/// one before the last, and so on, in the buffer it was given. A matrix of
/// an odd number of rows keeps its middle row where it is. `n` is 1 at
/// least, and a buffer that is not a whole number of rows of `n` keeps
/// what is left over at its end.
fn reverse_the_rows(values: &mut [f64], n: usize) {
    let mut rows = values.chunks_exact_mut(n);
    while let (Some(from_the_top), Some(from_the_bottom)) = (rows.next(), rows.next_back()) {
        from_the_top.swap_with_slice(from_the_bottom);
    }
}

/// The Cholesky factorization of the symmetric positive definite `a` of
/// `n` x `n`, given by its lower half: the lower triangular `l` with
/// `l l' = a`, which overwrites that lower half.
///
/// Only the lower half of `a` is read, the entries of column `j` at most
/// `i` of row `i`; what the upper half holds does not reach the result and
/// is left as it was. The lower half is the factorization when this comes
/// back, so a caller that needs `a` afterwards copies it first. `a` may
/// hold more values than `n` times `n`, and then its first `n` times `n`
/// are the matrix.
///
/// This is the test for a matrix that is positive definite, and the one
/// the association study fits its models through: it stops at the first
/// row whose diagonal entry, once the rows above it have been taken out,
/// is not above 0, and gives [`Error::Singular`] with that row. What that
/// row means is for the module that called, since it is what the data was
/// and not a defect of the call.
///
/// # Errors
///
/// [`Error::Dimension`] when `n` is 0, when `a` holds fewer than `n` times
/// `n` values, or when `n` times `n` is more than 2147483647, which is
/// what the routines of BLAS and LAPACK count in. [`Error::NotFinite`]
/// when the lower half of `a` holds a value that is not finite.
/// [`Error::Singular`] when `a` is not positive definite, with the row the
/// factorization stopped at. [`Error::NoConvergence`] when the routine
/// refused an argument it was given, which is a defect of popnei.
pub fn cholesky_lower(a: &mut [f64], n: usize) -> Result<()> {
    if n == 0 {
        return Err(Error::Dimension {
            argument: "n",
            expected: "1 at least, since a is the n x n matrix to factor".to_owned(),
        });
    }
    let a = the_matrix_of_mut(a, n, n, "a")?;
    refuse_a_value_that_is_not_finite_in_the_lower_half(a, n, "a")?;
    backend::cholesky_lower(a, n)
}

/// The `x` of `a x = b`, where `l` of `n` x `n` is the factorization
/// [`cholesky_lower`] gave of the symmetric positive definite `a`.
///
/// `b` is `sides` x `n`, row after row, one row for each right hand side,
/// and it comes back holding the solutions the same way. That layout is
/// the one the caller already has: where `a` is the coefficients square,
/// each right hand side holds one number for each coefficient and there is
/// one of them for each individual, so the matrix of them is the
/// individuals by the coefficients, which is the buffer the caller holds,
/// and it is passed with `sides` the individuals and `n` the coefficients;
/// nothing is copied and nothing is moved. numpy takes the right hand
/// sides as the columns of a matrix instead and is handed the transpose of
/// that buffer. The solutions come back the same way, one row for each
/// column of the matrix the next product reads, which is [`product`] with
/// that operand [`TheSecondOperand::ByTheColumnsOfTheResult`]. `sides` is
/// 1 at least.
///
/// Only the lower half of `l` is read, the entries of column `j` at most
/// `i` of row `i`; what the upper half holds does not reach the result.
/// Either buffer may hold more values than its dimensions ask for, and
/// then its first `n` times `n`, or `sides` times `n`, are the matrix.
///
/// # Errors
///
/// [`Error::Dimension`] when `n` or `sides` is 0, when `l` holds fewer
/// than `n` times `n` values or `b` fewer than `sides` times `n`, or when
/// either of those counts is more than 2147483647, which is what the
/// routines of BLAS and LAPACK count in. [`Error::NotFinite`] when the
/// lower half of `l`, or `b`, holds a value that is not finite.
/// [`Error::Singular`] when the diagonal of `l` holds an entry that is not
/// above 0, with the first such row, which is the same `Singular` a
/// [`cholesky_lower`] that gave this `l` would have given first: both
/// backends divide by that entry and answer with NaN and infinities,
/// where numpy refuses the system. [`Error::NoConvergence`] when the
/// routine refused an argument it was given, which is a defect of popnei.
///
/// What these do not catch: that `l` is the factorization of the `a` the
/// caller means to solve against. A slice that is some other lower
/// triangular matrix, the matrix `a` itself among them, is solved with and
/// gives a wrong answer and no error. Measured on 23 September 2026, the
/// 3 x 3 with rows (4, 2, 0), (2, 10, 6) and (0, 6, 5) passed where its
/// factorization was wanted gave (0.3848, 0.2304, 0.2160) for the right
/// hand side whose solution is (1, 2, 3). What a caller has instead of a
/// check is that `l` comes from [`cholesky_lower`].
pub fn solve_with_cholesky(l: &[f64], n: usize, b: &mut [f64], sides: usize) -> Result<()> {
    if n == 0 {
        return Err(Error::Dimension {
            argument: "n",
            expected: "1 at least, since l is the n x n factorization to solve with".to_owned(),
        });
    }
    if sides == 0 {
        return Err(Error::Dimension {
            argument: "sides",
            expected: "1 at least, since b holds one row for each right hand side".to_owned(),
        });
    }
    let l = the_matrix_of(l, n, n, "l")?;
    refuse_a_value_that_is_not_finite_in_the_lower_half(l, n, "l")?;
    let b = the_matrix_of_mut(b, sides, n, "b")?;
    refuse_a_value_that_is_not_finite(b, "b")?;
    refuse_a_diagonal_entry(l, n, "l", |entry| entry <= 0.0)?;
    backend::solve_with_cholesky(l, n, b, sides)
}

/// The log of the determinant of the `a` whose factorization `l` of `n` x
/// `n` is, which is twice the sum of the logs of the diagonal of `l`.
///
/// There is no sign to give: the determinant of a positive definite matrix
/// is above 0, so a caller that wanted numpy's `slogdet` has its second
/// value here and its first is always 1.
///
/// Neither backend runs: it is arithmetic over the `n` entries of the
/// diagonal, and the diagonal is all this reads. So the diagonal is what
/// it checks, and what it checks it for is what a [`cholesky_lower`] that
/// gave this `l` would have refused first. `l` may hold more values than
/// `n` times `n`, and then its first `n` times `n` are the matrix.
///
/// # Errors
///
/// [`Error::Dimension`] when `n` is 0, when `l` holds fewer than `n` times
/// `n` values, or when `n` times `n` is more than 2147483647, which is
/// what the routines of BLAS and LAPACK count in. [`Error::NotFinite`]
/// when the diagonal holds a value that is not finite.
/// [`Error::Singular`] when it holds one that is not above 0, with the
/// first such row.
///
/// What these do not catch: that `l` is a factorization at all. What comes
/// back is twice the sum of the logs of the diagonal of whatever slice was
/// given, which for one that is not a factorization is the log of the
/// determinant of that slice read as lower triangular and multiplied by
/// its own transpose, and not of anything the caller
/// holds. Measured on 23 September 2026,
/// the 3 x 3 with rows (4, 2, 0), (2, 10, 6) and (0, 6, 5) passed where
/// its factorization was wanted gave 10.596634733096074, where
/// 3.58351893845611 is the log of that matrix's determinant. What a caller
/// has instead of a check is that `l` comes from [`cholesky_lower`].
pub fn log_determinant_with_cholesky(l: &[f64], n: usize) -> Result<f64> {
    if n == 0 {
        return Err(Error::Dimension {
            argument: "n",
            expected: "1 at least, since l is the n x n factorization to read the diagonal of"
                .to_owned(),
        });
    }
    let l = the_matrix_of(l, n, n, "l")?;
    if !the_diagonal_of(l, n).all(f64::is_finite) {
        return Err(Error::NotFinite { argument: "l" });
    }
    refuse_a_diagonal_entry(l, n, "l", |entry| entry <= 0.0)?;
    // The logs are added in the order of the rows, which is the order the
    // spec's number was taken in and the one a total of floats has to be
    // added in to be the same number on every run.
    let total: f64 = the_diagonal_of(l, n).map(f64::ln).sum();
    Ok(2.0 * total)
}

/// The lower half of the inverse of the `a` whose factorization `l` of `n`
/// x `n` is, written into `inverse` of `n` x `n`.
///
/// `l` and `inverse` are two buffers and not one: the factorization is
/// left as it was, so a caller that has more to do with it, another solve
/// or the log of its determinant, still holds it. Only the lower half of
/// `inverse` is written, the entries of column `j` at most `i` of row `i`,
/// and its upper half is left as it was; the inverse of a symmetric matrix
/// is symmetric, so a caller that needs the whole of it mirrors that half.
///
/// Only the lower half of `l` is read, as [`solve_with_cholesky`] reads
/// it. Either buffer may hold more values than `n` times `n`, and then its
/// first `n` times `n` are the matrix.
///
/// This is the operation that asks for the most memory of its own, and
/// both backends ask: faer inverts into a scratch of `n` x `n`, and BLAS
/// builds the inverse from a triangular solve against the identity and a
/// product, which is two buffers of `n` x `n`, since neither routine may
/// write where it reads and the upper half of `inverse` is the caller's.
/// That is 800 MB and 1.6 GB at the 10000 individuals of
/// `docs/objectives.md`, and the crate asks for the memory instead of
/// taking it, so a machine without it gets [`Error::Memory`] and not the
/// end of the process. Why the BLAS backend does not use LAPACK's
/// `dpotri`, which inverts in place and needs no buffer, is in "Calling
/// the crate from two threads" of `docs/specs/linalg.md`: on Accelerate
/// that routine gives a numerically wrong inverse, about once in every
/// seven while another call of Accelerate is running on another thread.
///
/// # Errors
///
/// [`Error::Dimension`] when `n` is 0, when `l` or `inverse` holds fewer
/// than `n` times `n` values, or when `n` times `n` is more than
/// 2147483647, which is what the routines of BLAS and LAPACK count in.
/// [`Error::NotFinite`] when the lower half of `l` holds a value that is
/// not finite. [`Error::Singular`] when its diagonal holds an entry that
/// is not above 0, with the first such row, which is the same `Singular` a
/// [`cholesky_lower`] that gave this `l` would have given first.
/// [`Error::Memory`] when this machine has not the memory for faer's
/// scratch or for either of the two buffers of the BLAS backend.
/// [`Error::NoConvergence`] when the routine refused an argument
/// it was given, which is a defect of popnei.
///
/// What these do not catch: that `l` is the factorization of the matrix
/// the caller means to invert. Another lower triangular slice is inverted
/// as if it were one, and what comes back is the inverse of some other
/// matrix, with no error, as it is for [`solve_with_cholesky`]. What a
/// caller has instead of a check is that `l` comes from
/// [`cholesky_lower`].
pub fn invert_with_cholesky(l: &[f64], n: usize, inverse: &mut [f64]) -> Result<()> {
    if n == 0 {
        return Err(Error::Dimension {
            argument: "n",
            expected: "1 at least, since l is the n x n factorization to invert with".to_owned(),
        });
    }
    let l = the_matrix_of(l, n, n, "l")?;
    // The dimensions of both buffers are checked before anything else,
    // since a buffer that does not hold them is a defect of the call while
    // a matrix that cannot be worked with is what the data was.
    let inverse = the_matrix_of_mut(inverse, n, n, "inverse")?;
    refuse_a_value_that_is_not_finite_in_the_lower_half(l, n, "l")?;
    refuse_a_diagonal_entry(l, n, "l", |entry| entry <= 0.0)?;
    backend::invert_with_cholesky(l, n, inverse)
}

/// The thin QR factorization of a matrix with at least as many rows as
/// columns, which [`thin_qr`] gives: the two matrices whose product is
/// the matrix it was given.
#[derive(Debug, Clone)]
pub struct ThinQr {
    /// `rows` x `cols`, row after row: the columns are of length 1 and at
    /// right angles to each other. The sign of a column is the backend's.
    pub q: Vec<f64>,
    /// `cols` x `cols`, row after row: upper triangular, its lower half
    /// 0, with `a = q r`.
    pub r: Vec<f64>,
}

/// The thin QR factorization of `a` of `rows` x `cols`, row after row,
/// with `rows` at least `cols` and `cols` 1 at least.
///
/// The `q` of `rows` x `cols` has columns of length 1 that are at right
/// angles to each other, the `r` of `cols` x `cols` is upper triangular
/// with its lower half 0, and `a` is `q r`. Fitting a linear model to
/// more individuals than coefficients is this factorization of the design
/// and then a solve against its `r`, which is the "least squares" that
/// the `linalg` row of section 9 of `docs/architecture.md` names.
///
/// The sign of a column of `q`, and of the row of `r` that goes with it,
/// is whatever the backend gave, as the sign of an eigenvector of
/// [`eigh_lower`] is. Turning both round together leaves `q r` the matrix
/// it was, so a caller that fits a model sees nothing of it, and one that
/// reads a column of `q` by itself fixes the sign it wants.
///
/// `a` may hold more values than `rows` times `cols`, and then its first
/// `rows` times `cols` are the matrix. The whole of it is read, both
/// halves, a design being no more triangular than any other matrix.
///
/// # Errors
///
/// [`Error::Dimension`] when `cols` is 0, when `rows` is below `cols`,
/// when `a` holds fewer than `rows` times `cols` values, or when `rows`
/// times `cols` is more than 2147483647, which is what the routines of
/// BLAS and LAPACK count in. [`Error::NotFinite`] when `a` holds a value
/// that is not finite. [`Error::Memory`] when this machine has not the
/// memory for the copy of `a` the BLAS backend writes, which is as large
/// as `a`; faer is given the buffer as it lies and asks for nothing.
/// [`Error::NoConvergence`] when a routine refused an argument it was
/// given, which is a defect of popnei.
pub fn thin_qr(a: &[f64], rows: usize, cols: usize) -> Result<ThinQr> {
    if cols == 0 {
        return Err(Error::Dimension {
            argument: "cols",
            expected: "1 at least, since a is the rows x cols matrix to factor".to_owned(),
        });
    }
    if rows < cols {
        return Err(Error::Dimension {
            argument: "rows",
            expected: format!(
                "{cols} at least, the columns of a, since the thin QR is of a matrix with at least as many rows as columns, and it is {rows}"
            ),
        });
    }
    let a = the_matrix_of(a, rows, cols, "a")?;
    refuse_a_value_that_is_not_finite(a, "a")?;
    let mut q = vec![0.0_f64; a.len()];
    // `cols` times `cols` is at most the `rows` times `cols` that
    // `the_matrix_of` has just counted, so this refuses nothing that
    // reaches it: it is how the count is made without an arithmetic that
    // could overflow in silence.
    let mut r = vec![0.0_f64; the_values_of(cols, cols, "cols")?];
    backend::thin_qr(a, rows, cols, &mut q, &mut r)?;
    Ok(ThinQr { q, r })
}

/// Which half of the buffer holds the matrix that [`solve_triangular`]
/// is given.
///
/// Which of the two it is belongs here and not in the name of a function,
/// for the reason [`TheFirstOperand`] gives for the operands of a product:
/// the same buffer is a triangular matrix read either way, so a call that
/// named the wrong one of two functions would pass every check and come
/// back with the solution of another system and no error. The same
/// factorization against the same right hand side gives (4, 12, 3) read as
/// the lower half and (6.333333333333333, -4.666666666666666, 27) read as
/// the upper, both of them answers a caller could believe.
///
/// The two cases carry no value, since the half is all they say, and they
/// are an enum and not a `bool` because a `bool` at a call site says
/// nothing about which half it means.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TheHalfThatHoldsTheMatrix {
    /// The entries of column `j` at least `i` of row `i`, the half the
    /// `r` of a thin QR fills, whose diagonal belongs to both halves.
    TheUpperHalf,
    /// The entries of column `j` at most `i` of row `i`, the half a
    /// Cholesky factorization fills.
    TheLowerHalf,
}

/// The `x` of `a x = b` for the triangular `a` of `n` x `n`, held in the
/// half of its buffer the caller names, the other half not being read.
///
/// `b` is `sides` x `n`, row after row, one row for each right hand side,
/// and it comes back holding the solutions the same way, which is the
/// layout [`solve_with_cholesky`] takes and whose doc comment says where
/// that layout comes from. `sides` is 1 at least.
///
/// With [`TheHalfThatHoldsTheMatrix::TheUpperHalf`] this is the second
/// half of fitting a linear model to more individuals than coefficients:
/// [`thin_qr`] of the design gives the `q` and the `r`, and the
/// coefficients are the `c` of `r c = q' y` for the trait `y`. With
/// [`TheHalfThatHoldsTheMatrix::TheLowerHalf`] it is the solve against a
/// factor a Cholesky gave, which is lower triangular: the fit of the null
/// model of the logistic mixed model needs that solve with one right hand
/// side for each individual, and "The solve against a triangular matrix"
/// of `docs/specs/linalg.md` says what it is worth there.
///
/// Only the half named is read, and what the other half holds does not
/// reach the result; the diagonal belongs to both of them. Either buffer
/// may hold more values than its dimensions ask for, and then its first
/// `n` times `n`, or `sides` times `n`, are the matrix.
///
/// # Errors
///
/// [`Error::Dimension`] when `n` or `sides` is 0, when `a` holds fewer
/// than `n` times `n` values or `b` fewer than `sides` times `n`, or when
/// either of those counts is more than 2147483647, which is what the
/// routines of BLAS and LAPACK count in. [`Error::NotFinite`] when the
/// half of `a` the caller named, or `b`, holds a value that is not finite.
/// [`Error::Singular`] when the diagonal of `a` holds a 0, with the first
/// such row: the solve divides by every diagonal entry whichever half it
/// read, and the two backends part company on a 0 there, faer dividing by
/// it and answering with a NaN or an infinity, which of the two by the
/// half and the row, where `dtrtrs` gives an `info`, so the crate reads
/// that diagonal above them both. Only a 0 is refused, and the spec's
/// "The errors the seven add" has what a diagonal entry whose reciprocal
/// overflows gives instead.
/// [`Error::NoConvergence`] when the routine refused an argument it was
/// given, which is a defect of popnei.
///
/// What this does not catch: that `a` is triangular at all. A slice whose
/// other half holds something else is solved against as if that half were
/// 0, and what comes back is the solution of another system, with no
/// error, as it is for [`solve_with_cholesky`] and its `l`.
pub fn solve_triangular(
    a: &[f64],
    n: usize,
    half: TheHalfThatHoldsTheMatrix,
    b: &mut [f64],
    sides: usize,
) -> Result<()> {
    if n == 0 {
        return Err(Error::Dimension {
            argument: "n",
            expected: "1 at least, since a is the n x n triangular matrix to solve against"
                .to_owned(),
        });
    }
    if sides == 0 {
        return Err(Error::Dimension {
            argument: "sides",
            expected: "1 at least, since b holds one row for each right hand side".to_owned(),
        });
    }
    let a = the_matrix_of(a, n, n, "a")?;
    // The value that is not finite is read over the half the caller named,
    // since the other half is nothing of the matrix and does not reach the
    // result.
    match half {
        TheHalfThatHoldsTheMatrix::TheUpperHalf => {
            refuse_a_value_that_is_not_finite_in_the_upper_half(a, n, "a")?;
        }
        TheHalfThatHoldsTheMatrix::TheLowerHalf => {
            refuse_a_value_that_is_not_finite_in_the_lower_half(a, n, "a")?;
        }
    }
    let b = the_matrix_of_mut(b, sides, n, "b")?;
    refuse_a_value_that_is_not_finite(b, "b")?;
    // The diagonal belongs to both halves and the solve divides by every
    // entry of it, so it is read for a 0 whichever half was named.
    refuse_a_diagonal_entry(a, n, "a", |entry| entry == 0.0)?;
    backend::solve_triangular(a, n, half, b, sides)
}

/// The rank of `a` of `rows` x `cols`, row after row: how many of its
/// singular values are strictly above the tolerance numpy takes.
///
/// The singular values of a matrix are the factors by which it stretches
/// space along as many directions at right angles to each other as it has
/// columns, from the largest, and one of them is 0 exactly when a column
/// of the matrix is a combination of the others. So the rank is how many
/// of its columns are independent, at most the smaller of the two
/// dimensions, and it is what refuses a design whose covariates repeat
/// each other before any model is fitted.
///
/// The tolerance is the largest singular value times the larger dimension
/// times `f64::EPSILON`, the distance from 1 to the next `f64` above it,
/// which is 2.220446049250313e-16, and the last two are multiplied
/// together first, as numpy does it. That is the tolerance of
/// `matrix_rank` of numpy 2.5.3, read from that function on 23 September
/// 2026, and popnei takes it so that a design popnei refuses is a design
/// pyNei refuses. A matrix whose values are all 0 has a largest singular
/// value of 0 and so a tolerance of 0, and no value is strictly above
/// that, so its rank is 0.
///
/// `a` may hold more values than `rows` times `cols`, and then its first
/// `rows` times `cols` are the matrix. The whole of it is read, both
/// halves, a design being no more triangular than any other matrix.
///
/// # Errors
///
/// [`Error::Dimension`] when `rows` or `cols` is 0, when `a` holds fewer
/// than `rows` times `cols` values, or when that count, or the workspace
/// the routine asks for, is more than 2147483647, which is what the
/// routines of BLAS and LAPACK count in. [`Error::NotFinite`] when `a`
/// holds a value that is not finite. [`Error::Memory`] when this machine
/// has not the memory for the copy of `a` the BLAS backend writes, which
/// is as large as `a`; faer is given the buffer as it lies and asks for
/// nothing. [`Error::NoConvergence`] when the decomposition did not come
/// out, which is faer with its `SvdError::NoConvergence` and `dgesdd`
/// with an `info` above 0, and when the routine refused an argument it
/// was given, which is `dgesdd` with an `info` below 0 and a defect of
/// popnei; the two are the one case, with the routine and the `info`.
pub fn rank(a: &[f64], rows: usize, cols: usize) -> Result<usize> {
    if rows == 0 {
        return Err(Error::Dimension {
            argument: "rows",
            expected: "1 at least, since a is the rows x cols matrix to take the rank of"
                .to_owned(),
        });
    }
    if cols == 0 {
        return Err(Error::Dimension {
            argument: "cols",
            expected: "1 at least, since a is the rows x cols matrix to take the rank of"
                .to_owned(),
        });
    }
    let a = the_matrix_of(a, rows, cols, "a")?;
    refuse_a_value_that_is_not_finite(a, "a")?;
    let values = backend::singular_values(a, rows, cols)?;
    // Both backends give the values from the largest, and the largest is
    // taken here as the largest and not as the first, so that the
    // tolerance does not turn on that order. There are at most 46340 of
    // them, the smaller dimension of a matrix of 2147483647 values, and
    // the count below fits in a `usize` for the same reason.
    let largest = values.iter().copied().fold(0.0_f64, f64::max);
    // The larger dimension is multiplied by the distance from 1 to the
    // next `f64` first and the largest singular value last, which is
    // numpy's order. The other order is the same `f64` for every matrix
    // whose values a study holds, that distance being a power of two, and
    // it overflows to an infinity once the largest singular value passes
    // about 1e308 divided by the dimension, where no value is above the
    // tolerance and the rank comes back 0: measured on 23 September 2026,
    // the 2 x 2 with 1e308 and 1 on its diagonal is rank 1 this way and 0
    // the other. This product cannot overflow, the dimension times that
    // distance being below 1 and exact.
    let tolerance = largest * (rows.max(cols) as f64 * f64::EPSILON);
    Ok(values.iter().filter(|value| **value > tolerance).count())
}

/// How many values a matrix of `rows` x `cols` holds.
///
/// # Errors
///
/// [`Error::Dimension`] when they are more than
/// [`THE_MOST_VALUES_OF_A_MATRIX`], or more than a `usize` of this machine
/// counts.
fn the_values_of(rows: usize, cols: usize, argument: &'static str) -> Result<usize> {
    let too_many = || Error::Dimension {
        argument,
        expected: format!(
            "{rows} rows times {cols} columns, and a matrix holds at most {THE_MOST_VALUES_OF_A_MATRIX} values, which is what the routines of BLAS and LAPACK count in"
        ),
    };
    let wanted = rows.checked_mul(cols).ok_or_else(too_many)?;
    if wanted > THE_MOST_VALUES_OF_A_MATRIX {
        return Err(too_many());
    }
    Ok(wanted)
}

/// Refuses a `g` that does not hold exactly `n` times `n` values, which is
/// what both operations that take a `g` ask for: the buffer is the whole
/// matrix and not one with room to spare.
///
/// # Errors
///
/// [`Error::Dimension`] when it holds another number of values, and what
/// [`the_values_of`] gives.
fn refuse_a_g_that_is_not_square(held: usize, n: usize) -> Result<()> {
    let wanted = the_values_of(n, n, "g")?;
    if held == wanted {
        Ok(())
    } else {
        Err(Error::Dimension {
            argument: "g",
            expected: format!("{wanted} values, {n} rows times {n} columns, and it holds {held}"),
        })
    }
}

/// The first `rows` times `cols` values of the buffer, which are the
/// matrix held row after row.
///
/// # Errors
///
/// [`Error::Dimension`] when the buffer holds fewer, and what
/// [`the_values_of`] gives.
fn the_matrix_of<'a>(
    values: &'a [f64],
    rows: usize,
    cols: usize,
    argument: &'static str,
) -> Result<&'a [f64]> {
    let wanted = the_values_of(rows, cols, argument)?;
    values.get(..wanted).ok_or_else(|| Error::Dimension {
        argument,
        expected: format!(
            "{wanted} values at least, {rows} rows times {cols} columns, and it holds {held}",
            held = values.len()
        ),
    })
}

/// The same for a matrix the operation writes.
///
/// # Errors
///
/// [`Error::Dimension`] when the buffer holds fewer, and what
/// [`the_values_of`] gives.
fn the_matrix_of_mut<'a>(
    values: &'a mut [f64],
    rows: usize,
    cols: usize,
    argument: &'static str,
) -> Result<&'a mut [f64]> {
    let wanted = the_values_of(rows, cols, argument)?;
    let held = values.len();
    values.get_mut(..wanted).ok_or_else(|| Error::Dimension {
        argument,
        expected: format!(
            "{wanted} values at least, {rows} rows times {cols} columns, and it holds {held}"
        ),
    })
}

/// The bits of the exponent of an `f64`, which are all ones in an infinity
/// and in a NaN and in nothing else, so that a value is finite exactly
/// when its bits and this one are not this one.
const THE_EXPONENT_OF_A_VALUE_THAT_IS_NOT_FINITE: u64 = 0x7ff0_0000_0000_0000;

/// How many values one turn of [`every_value_is_finite`] reads, one to
/// each of its counters. Eight `f64` are four pairs of vector registers on
/// this machine.
const THE_VALUES_OF_ONE_TURN: usize = 8;

/// Whether every value of the buffer is finite.
///
/// It is written without a branch on the value so that the compiler
/// vectorizes it: each turn reads [`THE_VALUES_OF_ONE_TURN`] values and
/// puts into its own counter whether that value is an infinity or a NaN,
/// and the counters are read once at the end. The same test with a branch
/// per value compiles to one value per iteration, which is 1.68 ms of the
/// 12.05 ms of a product of a block of 5000 x 1000 with itself, measured
/// on 22 September 2026 and written in "Speed" of `docs/specs/linalg.md`.
fn every_value_is_finite(values: &[f64]) -> bool {
    let (turns, what_is_left_over) = values.as_chunks::<THE_VALUES_OF_ONE_TURN>();
    let mut what_is_not_finite = [0_u64; THE_VALUES_OF_ONE_TURN];
    for turn in turns {
        for (found, value) in what_is_not_finite.iter_mut().zip(turn) {
            *found |= u64::from(
                value.to_bits() & THE_EXPONENT_OF_A_VALUE_THAT_IS_NOT_FINITE
                    == THE_EXPONENT_OF_A_VALUE_THAT_IS_NOT_FINITE,
            );
        }
    }
    what_is_not_finite.iter().all(|found| *found == 0)
        && what_is_left_over.iter().all(|value| value.is_finite())
}

/// Refuses a matrix that holds an infinity or a NaN, which the backends do
/// not treat the same way and which a routine turns into a result that
/// says nothing of where it came from.
///
/// # Errors
///
/// [`Error::NotFinite`] when one of the values is one of those.
fn refuse_a_value_that_is_not_finite(values: &[f64], argument: &'static str) -> Result<()> {
    if every_value_is_finite(values) {
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
            .is_some_and(every_value_is_finite)
    });
    if it_is_all_finite {
        Ok(())
    } else {
        Err(Error::NotFinite { argument })
    }
}

/// The same for the upper half alone of an `n` x `n` matrix, the entries
/// of column `j` at least `i` of row `i`, which is what
/// [`solve_triangular`] reads when the caller names that half. `n` is 1 at
/// least.
///
/// # Errors
///
/// [`Error::NotFinite`] when one of those values is an infinity or a NaN.
fn refuse_a_value_that_is_not_finite_in_the_upper_half(
    values: &[f64],
    n: usize,
    argument: &'static str,
) -> Result<()> {
    let it_is_all_finite = values
        .chunks_exact(n)
        .enumerate()
        .all(|(row, entries)| entries.get(row..).is_some_and(every_value_is_finite));
    if it_is_all_finite {
        Ok(())
    } else {
        Err(Error::NotFinite { argument })
    }
}

/// Refuses a matrix whose diagonal holds an entry that the operation
/// cannot work with, and names the first row that holds one. The matrix is
/// `n` x `n` row after row and holds exactly that many values, and `n` is
/// 1 at least.
///
/// `it_is_refused` is what that operation cannot take: an entry that is
/// not above 0 for the `l` of a Cholesky factorization, which is what
/// [`cholesky_lower`] stops at and what the three operations that read
/// such an `l` give the same [`Error::Singular`] for.
///
/// The diagonal is read here, above the backends, and not left to them,
/// because on an `l` whose diagonal holds a 0 the two do not agree and,
/// where they agree, both are wrong. Measured on 23 September 2026 on such
/// an `l`: `dpotri` gave an `info` of 2 while faer's inverse gave no error
/// at all and wrote infinities and NaN into the buffer; and `dpotrs` gave
/// a solution of NaN, an infinity and an infinity with its sign turned
/// round while faer's solve gave three NaN, each with no error, where
/// numpy refuses the same system. No `l` that [`cholesky_lower`] gave is
/// such a matrix, since that is what it stops at, but a caller holds its
/// buffers apart and can pass one that never was a factorization.
///
/// # Errors
///
/// [`Error::Singular`] with the first row whose diagonal entry is refused.
fn refuse_a_diagonal_entry(
    values: &[f64],
    n: usize,
    argument: &'static str,
    it_is_refused: impl Fn(f64) -> bool,
) -> Result<()> {
    // The entry `i`, `i` of a matrix held row after row is the value `i`
    // times `n` plus `i` of the buffer, so the diagonal is one value in
    // every `n` plus 1 from the first, and the last of them is the last
    // value of the buffer. The addition cannot overflow: `the_matrix_of`
    // has refused every `n` whose square is above 2147483647, which leaves
    // `n` at 46340 at most.
    match the_diagonal_of(values, n).position(it_is_refused) {
        Some(row) => Err(Error::Singular { argument, at: row }),
        None => Ok(()),
    }
}

/// The `n` entries of the diagonal of a matrix of `n` x `n` values held
/// row after row, from the first row, with `n` 1 at least and the buffer
/// cut to exactly those values.
fn the_diagonal_of(values: &[f64], n: usize) -> impl Iterator<Item = f64> {
    values.iter().copied().step_by(n.saturating_add(1))
}

#[cfg(test)]
mod tests {
    use super::{
        Eigen, Error, TheFirstOperand, TheHalfThatHoldsTheMatrix, TheSecondOperand, ThinQr,
        add_self_product_lower, cholesky_lower, eigh_lower, invert_with_cholesky,
        log_determinant_with_cholesky, product, rank, reverse_the_rows, solve_triangular,
        solve_with_cholesky, thin_qr,
    };

    /// The A of 2 x 3 of "How it is verified" of `docs/specs/linalg.md`,
    /// rows (1, 2, 0) and (0, 1, 3), row after row.
    const A_OF_2_BY_3: [f64; 6] = [1.0, 2.0, 0.0, 0.0, 1.0, 3.0];

    /// The first B of the same place, 3 x 2, rows (1, 0), (2, 1) and
    /// (0, 3). A times it is symmetric.
    const B_OF_3_BY_2: [f64; 6] = [1.0, 0.0, 2.0, 1.0, 0.0, 3.0];

    /// The second B, 3 x 2, rows (1, 1), (2, 0) and (0, 3). A times it is
    /// not symmetric, so a backend that wrote the transpose of C would
    /// give another answer.
    const B_THAT_IS_NOT_SYMMETRIC: [f64; 6] = [1.0, 1.0, 2.0, 0.0, 0.0, 3.0];

    /// The third B, 3 x 1, rows (1), (0) and (2). With it the three
    /// dimensions of the product are all different.
    const B_OF_3_BY_1: [f64; 3] = [1.0, 0.0, 2.0];

    /// The B of the first case of `product_with_the_second_turned` of
    /// the same place, 2 x 3, rows (1, 1, 0) and (0, 2, 1). A times its
    /// transpose is not symmetric.
    const B_OF_2_BY_3: [f64; 6] = [1.0, 1.0, 0.0, 0.0, 2.0, 1.0];

    /// The B of its second case, 1 x 3, the row (2, 0, 1). With it the
    /// three dimensions of the product are all different.
    const B_OF_1_BY_3: [f64; 3] = [2.0, 0.0, 1.0];

    /// The A of 2 x 3 written the other way round, 3 x 2 with rows
    /// (1, 0), (2, 1) and (0, 3), which is what a caller holds that has
    /// its first operand laid out by the values the product sums over.
    const A_OF_2_BY_3_THE_OTHER_WAY_ROUND: [f64; 6] = [1.0, 0.0, 2.0, 1.0, 0.0, 3.0];

    /// The B that is not symmetric written the other way round, 2 x 3
    /// with rows (1, 2, 0) and (1, 0, 3), which is what a caller holds
    /// that has its second operand laid out by the columns of the result.
    const B_THAT_IS_NOT_SYMMETRIC_THE_OTHER_WAY_ROUND: [f64; 6] = [1.0, 2.0, 0.0, 1.0, 0.0, 3.0];

    /// The 2 x 2 that the four combinations all write, rows (5, 1) and
    /// (2, 9), which is not symmetric.
    const THE_PRODUCT_THAT_IS_NOT_SYMMETRIC: [f64; 4] = [5.0, 1.0, 2.0, 9.0];

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
    /// do not hold the same number of entries. A value that is not a
    /// number differs from every value, so that a result of NaN fails the
    /// test that compares it instead of passing every comparison.
    fn differ(got: &[f64], expected: &[f64], tolerance: f64) -> bool {
        got.len() != expected.len()
            || got.iter().zip(expected).any(|(one, other)| {
                let difference = (one - other).abs();
                !difference.is_finite() || difference > tolerance
            })
    }

    /// The two differ by more than `tolerance` times the size of the value
    /// that was expected, which is how the spec compares an eigenvalue: an
    /// eigenvalue of 361 and one of 0.79 are asked for the same number of
    /// digits and not for the same absolute error. A value that is not a
    /// number differs, as above.
    fn differ_in_their_digits(got: f64, expected: f64, tolerance: f64) -> bool {
        let difference = (got - expected).abs();
        !difference.is_finite() || difference > tolerance * expected.abs()
    }

    /// The first `how_many` numbers of the generator of "How it is
    /// verified" of `docs/specs/linalg.md`, so that the test and numpy
    /// make the same matrix: a 64 bit state that starts at 7 with its
    /// lowest bit set, and for each number `s ^= s << 13; s ^= s >> 7;
    /// s ^= s << 17`, the shifts dropping the bits that leave the 64, and
    /// the number is `(s >> 11) / 2^53 - 0.5`.
    fn the_numbers_of_the_generator(how_many: usize) -> Vec<f64> {
        let mut state = 7_u64;
        (0..how_many)
            .map(|_| {
                state ^= state.wrapping_shl(13);
                state ^= state.wrapping_shr(7);
                state ^= state.wrapping_shl(17);
                // 2^53, below which a count is an exact f64.
                state.wrapping_shr(11) as f64 / 9007199254740992.0 - 0.5
            })
            .collect()
    }

    /// The 3 x 3 symmetric positive definite matrix of "How the seven are
    /// verified" of `docs/specs/linalg.md`, rows (4, 2, 0), (2, 10, 6) and
    /// (0, 6, 5), with only its lower half given and the upper half
    /// holding a value that is nothing of the matrix, so that a call that
    /// read the upper half instead would give another factorization.
    fn the_matrix_to_factor() -> Vec<f64> {
        vec![
            4.0, 99.0, 99.0, //
            2.0, 10.0, 99.0, //
            0.0, 6.0, 5.0,
        ]
    }

    /// The right hand side of "How the seven are verified", (8, 40, 27),
    /// whose solution against the 3 x 3 above is (1, 2, 3).
    const THE_RIGHT_HAND_SIDE: [f64; 3] = [8.0, 40.0, 27.0];

    /// The second right hand side of the same place, (4, 2, 0), which is
    /// the first column of that 3 x 3, so its solution is (1, 0, 0).
    const THE_SECOND_RIGHT_HAND_SIDE: [f64; 3] = [4.0, 2.0, 0.0];

    /// The order of the matrix below, which is the number of individuals
    /// the spec checks the two backends at.
    const THE_LARGE_CASE: usize = 1000;

    /// The columns of the Z that matrix is built from, which are the
    /// variants of those individuals.
    const THE_VARIANTS_OF_THE_LARGE_CASE: usize = 1200;

    /// The G of "How the seven are verified", the lower half of ZZ' for
    /// the Z of 1000 rows and 1200 columns of the generator, which is the
    /// matrix the eigendecomposition is checked on as well: ZZ' is A'A for
    /// A = Z', of 1200 rows and 1000 columns, and the numbers in the order
    /// the generator gives them are the rows of that A, one after another.
    fn the_matrix_of_1000_by_1000_of_the_generator() -> Vec<f64> {
        let a = the_numbers_of_the_generator(THE_VARIANTS_OF_THE_LARGE_CASE * THE_LARGE_CASE);
        let mut g = vec![0.0_f64; THE_LARGE_CASE * THE_LARGE_CASE];
        add_self_product_lower(&a, THE_VARIANTS_OF_THE_LARGE_CASE, THE_LARGE_CASE, &mut g).unwrap();
        g
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
    fn the_self_product_reads_the_first_rows_of_an_a_that_holds_more() {
        // The same A with a seventh value, which is no row of a 2 x 3.
        let a = [1.0, 2.0, 0.0, 0.0, 1.0, 3.0, 1000.0];
        let mut g = vec![0.0; 9];
        add_self_product_lower(&a, 2, 3, &mut g).unwrap();
        assert_eq!(
            g,
            vec![
                1.0, 0.0, 0.0, //
                2.0, 5.0, 0.0, //
                0.0, 3.0, 9.0,
            ]
        );
    }

    #[test]
    fn the_product_of_two_matrices_that_are_not_square() {
        // C holds other values first, so a product that added to C instead
        // of overwriting it would leave them in the result.
        let mut c = vec![7.0; 4];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &B_OF_3_BY_2,
                cols: 2,
            },
            &mut c,
        )
        .unwrap();
        assert_eq!(c, vec![5.0, 2.0, 2.0, 10.0]);
    }

    #[test]
    fn the_product_of_a_b_whose_result_is_not_symmetric() {
        let mut c = vec![7.0; 4];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &B_THAT_IS_NOT_SYMMETRIC,
                cols: 2,
            },
            &mut c,
        )
        .unwrap();
        assert_eq!(c, vec![5.0, 1.0, 2.0, 9.0]);
    }

    #[test]
    fn the_product_of_matrices_whose_three_dimensions_are_all_different() {
        let mut c = vec![7.0; 2];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &B_OF_3_BY_1,
                cols: 1,
            },
            &mut c,
        )
        .unwrap();
        assert_eq!(c, vec![1.0, 6.0]);
    }

    #[test]
    fn the_product_reads_and_writes_the_first_values_of_buffers_that_hold_more() {
        let b = [1.0, 0.0, 2.0, 1.0, 0.0, 3.0, 1000.0];
        let mut c = vec![7.0; 5];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &b,
                cols: 2,
            },
            &mut c,
        )
        .unwrap();
        assert_eq!(c, vec![5.0, 2.0, 2.0, 10.0, 7.0]);
    }

    #[test]
    fn the_product_of_an_a_of_no_rows_writes_nothing() {
        let mut c = vec![7.0, 7.0];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &[],
                rows: 0,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &B_OF_3_BY_2,
                cols: 2,
            },
            &mut c,
        )
        .unwrap();
        assert_eq!(c, vec![7.0, 7.0]);
    }

    #[test]
    fn the_product_by_the_columns_of_the_result_of_two_matrices_whose_result_is_not_symmetric() {
        // C holds other values first, so a product that added to C instead
        // of overwriting it would leave them in the result. The result is
        // not symmetric, so a backend that wrote the transpose of C would
        // fail this.
        let mut c = vec![7.0; 4];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &B_OF_2_BY_3,
                cols: 2,
            },
            &mut c,
        )
        .unwrap();
        assert_eq!(c, vec![3.0, 4.0, 1.0, 5.0]);
    }

    #[test]
    fn the_product_by_the_columns_of_the_result_of_matrices_whose_three_dimensions_are_all_different()
     {
        let mut c = vec![7.0; 2];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &B_OF_1_BY_3,
                cols: 1,
            },
            &mut c,
        )
        .unwrap();
        assert_eq!(c, vec![2.0, 3.0]);
    }

    #[test]
    fn the_product_by_the_columns_of_the_result_of_a_matrix_with_itself_is_the_product_of_the_two_layouts()
     {
        // A A', which is the first case of `product` above written the
        // other way round: the B of 3 x 2 there is this A. The two
        // functions are asserted to give it alike, which is what catches
        // one of them reading an operand the way the other does.
        let mut c = vec![7.0; 4];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &A_OF_2_BY_3,
                cols: 2,
            },
            &mut c,
        )
        .unwrap();
        assert_eq!(c, vec![5.0, 2.0, 2.0, 10.0]);
        let mut of_the_product = vec![7.0; 4];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &B_OF_3_BY_2,
                cols: 2,
            },
            &mut of_the_product,
        )
        .unwrap();
        assert_eq!(c, of_the_product);
    }

    #[test]
    fn the_product_by_the_columns_of_the_result_of_an_a_of_no_rows_writes_nothing() {
        let mut c = vec![7.0, 7.0];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &[],
                rows: 0,
            },
            3,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &B_OF_2_BY_3,
                cols: 2,
            },
            &mut c,
        )
        .unwrap();
        assert_eq!(c, vec![7.0, 7.0]);
    }

    /// The four combinations of "How it is verified" of "The product with
    /// its first operand turned", each given the matrices its two cases
    /// say it holds, write the same 2 x 2, rows (5, 1) and (2, 9). That
    /// matrix is not symmetric, so a backend that wrote the transpose of
    /// the result fails here; and the four pairs differ in which buffer is
    /// written the other way round, so a backend that read an operand the
    /// way another combination reads it gives another matrix for at least
    /// one of the four.
    ///
    /// It is also what tells the two cases of the first operand apart: a
    /// buffer holds rows times inner values whichever way round it is
    /// read, so no check of a length can, and only the matrix can.
    #[test]
    fn the_same_matrix_comes_out_of_the_product_four_ways() {
        // Each `c` holds other values first, so a product that added to it
        // instead of overwriting it would leave them in the result.
        let mut of_a_b = vec![7.0; 4];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &B_THAT_IS_NOT_SYMMETRIC,
                cols: 2,
            },
            &mut of_a_b,
        )
        .unwrap();
        assert_eq!(of_a_b, THE_PRODUCT_THAT_IS_NOT_SYMMETRIC);

        let mut of_a_turned_and_b = vec![7.0; 4];
        product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &A_OF_2_BY_3_THE_OTHER_WAY_ROUND,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &B_THAT_IS_NOT_SYMMETRIC,
                cols: 2,
            },
            &mut of_a_turned_and_b,
        )
        .unwrap();
        assert_eq!(of_a_turned_and_b, THE_PRODUCT_THAT_IS_NOT_SYMMETRIC);

        let mut of_a_and_b_turned = vec![7.0; 4];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &B_THAT_IS_NOT_SYMMETRIC_THE_OTHER_WAY_ROUND,
                cols: 2,
            },
            &mut of_a_and_b_turned,
        )
        .unwrap();
        assert_eq!(of_a_and_b_turned, THE_PRODUCT_THAT_IS_NOT_SYMMETRIC);

        let mut of_both_turned = vec![7.0; 4];
        product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &A_OF_2_BY_3_THE_OTHER_WAY_ROUND,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &B_THAT_IS_NOT_SYMMETRIC_THE_OTHER_WAY_ROUND,
                cols: 2,
            },
            &mut of_both_turned,
        )
        .unwrap();
        assert_eq!(of_both_turned, THE_PRODUCT_THAT_IS_NOT_SYMMETRIC);
    }

    #[test]
    fn the_product_by_the_columns_of_the_result_refuses_the_dimensions_the_product_refuses() {
        // A `b` that does not hold its cols rows of inner values, a `c`
        // shorter than its rows times its columns, and the two dimensions
        // that are 1 at least.
        let mut c = vec![0.0; 4];
        let error = product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &B_OF_2_BY_3[..4],
                cols: 2,
            },
            &mut c,
        )
        .unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "b", .. }),
            "the error is {error}"
        );
        let mut short = vec![0.0; 3];
        let error = product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &B_OF_2_BY_3,
                cols: 2,
            },
            &mut short,
        )
        .unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "c", .. }),
            "the error is {error}"
        );
        let error = product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &B_OF_2_BY_3,
                cols: 0,
            },
            &mut c,
        )
        .unwrap_err();
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
        let error = product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &[],
                rows: 2,
            },
            0,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &[],
                cols: 2,
            },
            &mut c,
        )
        .unwrap_err();
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
        let mut with_one_that_is_not_finite = B_OF_2_BY_3;
        with_one_that_is_not_finite[3] = f64::NAN;
        let error = product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &with_one_that_is_not_finite,
                cols: 2,
            },
            &mut c,
        )
        .unwrap_err();
        assert!(
            matches!(error, Error::NotFinite { argument: "b" }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_eigendecomposition_gives_the_values_from_the_largest_and_the_vectors_as_rows() {
        let Eigen { values, vectors } = eigh_lower(the_matrix_of_3_by_3(), 3).unwrap();
        // (7 + √5) / 2, (7 - √5) / 2 and 1, all the digits numpy prints.
        assert!(
            !differ(&values, &[4.618033988749895, 2.381966011250105, 1.0], 1e-12),
            "the eigenvalues are {values:?}"
        );
        let expected = [
            [0.8506508083520399, 0.5257311121191335, 0.0],
            [-0.5257311121191335, 0.8506508083520399, 0.0],
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
    fn the_eigendecomposition_of_the_1000_by_1000_matrix_of_the_generator() {
        // The G of "How it is verified" is ZZ' for a Z of 1000 rows and
        // 1200 columns whose z(i, c) is the (c * 1000 + i)-th number of the
        // generator. ZZ' is A'A for A = Z', of 1200 rows and 1000 columns,
        // and the numbers in the order the generator gives them are the
        // rows of that A, one after another. The matrix is of the size the
        // two backends run at, and it has full rank, so no two of its
        // eigenvalues are closer than 0.003 and each eigenvector is the
        // only one of its eigenvalue up to its sign.
        const INDIVIDUALS: usize = 1000;
        const VARIANTS: usize = 1200;
        let a = the_numbers_of_the_generator(VARIANTS * INDIVIDUALS);
        let mut g = vec![0.0_f64; INDIVIDUALS * INDIVIDUALS];
        add_self_product_lower(&a, VARIANTS, INDIVIDUALS, &mut g).unwrap();

        let trace: f64 = g
            .as_chunks::<INDIVIDUALS>()
            .0
            .iter()
            .enumerate()
            .map(|(row, entries)| entries[row])
            .sum();
        assert!(
            !differ_in_their_digits(trace, 99996.3873081677, 1e-12),
            "the trace of the product is {trace}"
        );

        let Eigen { values, vectors } = eigh_lower(g, INDIVIDUALS).unwrap();
        assert_eq!(values.len(), INDIVIDUALS);
        let sum: f64 = values.iter().sum();
        assert!(
            !differ_in_their_digits(sum, 99996.38730816769, 1e-12),
            "the eigenvalues add up to {sum}"
        );
        for (got, expected) in
            values
                .iter()
                .zip([361.9125119011332, 359.66517178659313, 356.4439329949563])
        {
            assert!(
                !differ_in_their_digits(*got, expected, 1e-12),
                "an eigenvalue among the three largest is {got} and not {expected}"
            );
        }
        let smallest = values[INDIVIDUALS - 1];
        assert!(
            !differ_in_their_digits(smallest, 0.7933289215408484, 1e-12),
            "the smallest eigenvalue is {smallest}"
        );

        // An eigenvector is less well determined than its eigenvalue by the
        // gap to its neighbours, so its entries are compared to 1e-9.
        let first = with_the_sign_of_the_spec(&vectors[..INDIVIDUALS]);
        assert!(
            !differ(
                &first[..3],
                &[
                    0.018111301861995926,
                    -0.004576022169100455,
                    -0.00518748985091827
                ],
                1e-9
            ),
            "the first three entries of the eigenvector of the largest eigenvalue are {:?}",
            &first[..3]
        );
    }

    #[test]
    fn the_rows_of_an_odd_number_are_turned_round_about_the_middle_one() {
        let mut values = vec![
            1.0, 2.0, //
            3.0, 4.0, //
            5.0, 6.0,
        ];
        reverse_the_rows(&mut values, 2);
        assert_eq!(values, vec![5.0, 6.0, 3.0, 4.0, 1.0, 2.0]);
    }

    #[test]
    fn the_self_product_refuses_a_g_that_is_not_cols_by_cols() {
        for held in [8_usize, 10] {
            let mut g = vec![0.0; held];
            let error = add_self_product_lower(&A_OF_2_BY_3, 2, 3, &mut g).unwrap_err();
            assert!(
                matches!(error, Error::Dimension { argument: "g", .. }),
                "the error for a g of {held} values is {error}"
            );
        }
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
        // Row 1, column 0, below the diagonal, and row 1, column 1, on it:
        // the product reads both. The same value in the upper half is not
        // read and is not an error, which the last lines check.
        for entry in [3_usize, 4] {
            let mut g = vec![0.0; 9];
            *g.get_mut(entry).unwrap() = f64::INFINITY;
            let error = add_self_product_lower(&A_OF_2_BY_3, 2, 3, &mut g).unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "g" }),
                "the error for the entry {entry} is {error}"
            );
        }
        let mut g = vec![0.0; 9];
        *g.get_mut(1).unwrap() = f64::INFINITY;
        add_self_product_lower(&A_OF_2_BY_3, 2, 3, &mut g).unwrap();
    }

    /// The two matrices of these tests hold 6 and 9 values, and the scans
    /// read 8 at a turn, so a value that is not finite in them is found by
    /// the few values left over at the end and never by the turns. These
    /// two put one in every place of a matrix of 2 x 12 and of the lower
    /// half of a g of 12 x 12, which the turns do read.
    #[test]
    fn the_self_product_refuses_an_a_whose_value_that_is_not_finite_is_at_any_place() {
        let rows = 2;
        let cols = 12;
        for place in 0..rows * cols {
            let mut a = vec![1.0; rows * cols];
            *a.get_mut(place).unwrap() = if place % 2 == 0 {
                f64::NAN
            } else {
                f64::NEG_INFINITY
            };
            let mut g = vec![0.0; cols * cols];
            let error = add_self_product_lower(&a, rows, cols, &mut g).unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "a" }),
                "the error for the place {place} is {error}"
            );
        }
        let a = vec![1.0; rows * cols];
        let mut g = vec![0.0; cols * cols];
        add_self_product_lower(&a, rows, cols, &mut g).unwrap();
    }

    #[test]
    fn the_self_product_refuses_a_g_whose_value_that_is_not_finite_is_at_any_place_of_its_lower_half()
     {
        let n = 12;
        let a = vec![1.0; 2 * n];
        for row in 0..n {
            for column in 0..n {
                let mut g = vec![0.0; n * n];
                *g.get_mut(row * n + column).unwrap() = f64::INFINITY;
                let result = add_self_product_lower(&a, 2, n, &mut g);
                if column <= row {
                    let error = result.unwrap_err();
                    assert!(
                        matches!(error, Error::NotFinite { argument: "g" }),
                        "the error for the row {row} and the column {column} is {error}"
                    );
                } else {
                    result.unwrap_or_else(|error| {
                        panic!("the row {row} and the column {column} gave {error}")
                    });
                }
            }
        }
    }

    #[test]
    fn the_product_refuses_a_b_that_does_not_have_the_inner_dimension() {
        let mut c = vec![0.0; 4];
        let error = product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &B_OF_3_BY_2[..4],
                cols: 2,
            },
            &mut c,
        )
        .unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "b", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_product_refuses_a_c_shorter_than_its_rows_times_its_columns() {
        let mut c = vec![0.0; 3];
        let error = product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &B_OF_3_BY_2,
                cols: 2,
            },
            &mut c,
        )
        .unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "c", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_product_refuses_a_cols_of_zero_and_an_inner_of_zero() {
        let mut c: Vec<f64> = Vec::new();
        let error = product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &B_OF_3_BY_2,
                cols: 0,
            },
            &mut c,
        )
        .unwrap_err();
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
        let error = product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &[],
                rows: 2,
            },
            0,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &[],
                cols: 2,
            },
            &mut c,
        )
        .unwrap_err();
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
        let error = product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &A_OF_2_BY_3,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &b,
                cols: 2,
            },
            &mut c,
        )
        .unwrap_err();
        assert!(
            matches!(error, Error::NotFinite { argument: "b" }),
            "the error is {error}"
        );
    }

    /// The first operand of a product is checked as the second one is:
    /// its length against its rows times the inner dimension, and its
    /// values for one that is not finite. Both are made whichever way the
    /// second operand is held, and both name `a` and not `b`.
    #[test]
    fn the_product_refuses_an_a_that_is_not_its_rows_times_the_inner_dimension() {
        let of_the_values = TheSecondOperand::ByTheValuesSummedOver {
            values: &B_OF_3_BY_2,
            cols: 2,
        };
        let of_the_columns = TheSecondOperand::ByTheColumnsOfTheResult {
            values: &B_OF_2_BY_3,
            cols: 2,
        };
        for b in [of_the_values, of_the_columns] {
            let mut c = vec![0.0; 4];
            let error = product(
                TheFirstOperand::ByTheRowsOfTheResult {
                    values: &A_OF_2_BY_3[..5],
                    rows: 2,
                },
                3,
                b,
                &mut c,
            )
            .unwrap_err();
            assert!(
                matches!(error, Error::Dimension { argument: "a", .. }),
                "the error for {b:?} is {error}"
            );
            // The same call with the whole of `a` is no error.
            product(
                TheFirstOperand::ByTheRowsOfTheResult {
                    values: &A_OF_2_BY_3,
                    rows: 2,
                },
                3,
                b,
                &mut c,
            )
            .unwrap();
        }
    }

    #[test]
    fn the_product_refuses_a_value_that_is_not_finite_in_a() {
        let of_the_values = TheSecondOperand::ByTheValuesSummedOver {
            values: &B_OF_3_BY_2,
            cols: 2,
        };
        let of_the_columns = TheSecondOperand::ByTheColumnsOfTheResult {
            values: &B_OF_2_BY_3,
            cols: 2,
        };
        for b in [of_the_values, of_the_columns] {
            for place in 0..A_OF_2_BY_3.len() {
                let mut a = A_OF_2_BY_3;
                *a.get_mut(place).unwrap() = if place % 2 == 0 {
                    f64::NAN
                } else {
                    f64::NEG_INFINITY
                };
                let mut c = vec![0.0; 4];
                let error = product(
                    TheFirstOperand::ByTheRowsOfTheResult {
                        values: &a,
                        rows: 2,
                    },
                    3,
                    b,
                    &mut c,
                )
                .unwrap_err();
                assert!(
                    matches!(error, Error::NotFinite { argument: "a" }),
                    "the error for the place {place} of {b:?} is {error}"
                );
            }
        }
    }

    /// The two combinations that turn the first operand, each on a result
    /// whose three dimensions differ. The four cases above are all 2 x 2,
    /// so the rows and the columns of the result are the same number and a
    /// backend that gave a routine the one where the other goes writes
    /// every one of them; these two are what catches that. Both take the A
    /// of 2 x 3 written the other way round, with `inner` 3 and one
    /// column, and their two matrices are the ones "How it is verified" of
    /// "The product with its first operand turned" gives.
    #[test]
    fn the_product_that_turns_its_first_operand_writes_a_result_that_is_not_square() {
        let mut of_a_turned_and_b = vec![7.0; 2];
        product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &A_OF_2_BY_3_THE_OTHER_WAY_ROUND,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &B_OF_3_BY_1,
                cols: 1,
            },
            &mut of_a_turned_and_b,
        )
        .unwrap();
        assert_eq!(of_a_turned_and_b, vec![1.0, 6.0]);

        let mut of_both_turned = vec![7.0; 2];
        product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &A_OF_2_BY_3_THE_OTHER_WAY_ROUND,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &B_OF_1_BY_3,
                cols: 1,
            },
            &mut of_both_turned,
        )
        .unwrap();
        assert_eq!(of_both_turned, vec![2.0, 3.0]);
    }

    /// The `a` of no rows of the same two combinations: nothing is
    /// written, as it is for the two that do not turn it.
    #[test]
    fn the_product_that_turns_its_first_operand_of_an_a_of_no_rows_writes_nothing() {
        for b in [
            TheSecondOperand::ByTheValuesSummedOver {
                values: &B_THAT_IS_NOT_SYMMETRIC,
                cols: 2,
            },
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &B_THAT_IS_NOT_SYMMETRIC_THE_OTHER_WAY_ROUND,
                cols: 2,
            },
        ] {
            let mut c = vec![7.0, 7.0];
            product(
                TheFirstOperand::ByTheValuesSummedOver {
                    values: &[],
                    rows: 0,
                },
                3,
                b,
                &mut c,
            )
            .unwrap();
            assert_eq!(c, vec![7.0, 7.0], "the result for {b:?}");
        }
    }

    /// The two combinations whose first operand is held by the values
    /// summed over refuse what the two that hold it by the rows of the
    /// result refuse, and each buffer is checked against the dimensions of
    /// the way it is read: `a` as inner x rows here, which the message
    /// says, and `b` as its own case says.
    #[test]
    fn the_product_that_turns_its_first_operand_refuses_the_dimensions_the_others_refuse() {
        let of_the_values = TheSecondOperand::ByTheValuesSummedOver {
            values: &B_THAT_IS_NOT_SYMMETRIC,
            cols: 2,
        };
        let of_the_values_one_short = TheSecondOperand::ByTheValuesSummedOver {
            values: &B_THAT_IS_NOT_SYMMETRIC[..5],
            cols: 2,
        };
        let of_the_columns = TheSecondOperand::ByTheColumnsOfTheResult {
            values: &B_THAT_IS_NOT_SYMMETRIC_THE_OTHER_WAY_ROUND,
            cols: 2,
        };
        let of_the_columns_one_short = TheSecondOperand::ByTheColumnsOfTheResult {
            values: &B_THAT_IS_NOT_SYMMETRIC_THE_OTHER_WAY_ROUND[..5],
            cols: 2,
        };
        for (b, one_short, the_dimensions_of_b) in [
            (
                of_the_values,
                of_the_values_one_short,
                "3 rows times 2 columns",
            ),
            (
                of_the_columns,
                of_the_columns_one_short,
                "2 rows times 3 columns",
            ),
        ] {
            let a = TheFirstOperand::ByTheValuesSummedOver {
                values: &A_OF_2_BY_3_THE_OTHER_WAY_ROUND,
                rows: 2,
            };
            // The call with nothing wrong gives the matrix, so each case
            // below is refused for the one thing that case changes.
            let mut c = vec![0.0; 4];
            product(a, 3, b, &mut c).unwrap();
            assert_eq!(c, THE_PRODUCT_THAT_IS_NOT_SYMMETRIC);

            let error = product(
                TheFirstOperand::ByTheValuesSummedOver {
                    values: &A_OF_2_BY_3_THE_OTHER_WAY_ROUND[..5],
                    rows: 2,
                },
                3,
                b,
                &mut c,
            )
            .unwrap_err();
            assert!(
                matches!(error, Error::Dimension { argument: "a", .. })
                    && error.to_string().contains("3 rows times 2 columns"),
                "the error for {b:?} is {error}"
            );

            let error = product(a, 3, one_short, &mut c).unwrap_err();
            assert!(
                matches!(error, Error::Dimension { argument: "b", .. })
                    && error.to_string().contains(the_dimensions_of_b),
                "the error for {b:?} is {error}"
            );

            let mut shorter_than_the_result = vec![0.0; 3];
            let error = product(a, 3, b, &mut shorter_than_the_result).unwrap_err();
            assert!(
                matches!(error, Error::Dimension { argument: "c", .. }),
                "the error for {b:?} is {error}"
            );
        }

        // The two dimensions that are 1 at least, which no case of an
        // operand changes.
        let mut c = vec![0.0; 4];
        let error = product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &A_OF_2_BY_3_THE_OTHER_WAY_ROUND,
                rows: 2,
            },
            3,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &B_THAT_IS_NOT_SYMMETRIC_THE_OTHER_WAY_ROUND,
                cols: 0,
            },
            &mut c,
        )
        .unwrap_err();
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
        let error = product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &[],
                rows: 2,
            },
            0,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &[],
                cols: 2,
            },
            &mut c,
        )
        .unwrap_err();
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

    /// The same two combinations refuse a value that is not finite in
    /// either operand, wherever in the buffer it sits.
    #[test]
    fn the_product_that_turns_its_first_operand_refuses_a_value_that_is_not_finite() {
        let mut of_the_values_with_one_that_is_not_finite = B_THAT_IS_NOT_SYMMETRIC;
        of_the_values_with_one_that_is_not_finite[3] = f64::NAN;
        let mut of_the_columns_with_one_that_is_not_finite =
            B_THAT_IS_NOT_SYMMETRIC_THE_OTHER_WAY_ROUND;
        of_the_columns_with_one_that_is_not_finite[3] = f64::NEG_INFINITY;
        let of_the_values = TheSecondOperand::ByTheValuesSummedOver {
            values: &B_THAT_IS_NOT_SYMMETRIC,
            cols: 2,
        };
        let of_the_columns = TheSecondOperand::ByTheColumnsOfTheResult {
            values: &B_THAT_IS_NOT_SYMMETRIC_THE_OTHER_WAY_ROUND,
            cols: 2,
        };
        for (b, with_one_that_is_not_finite) in [
            (
                of_the_values,
                TheSecondOperand::ByTheValuesSummedOver {
                    values: &of_the_values_with_one_that_is_not_finite,
                    cols: 2,
                },
            ),
            (
                of_the_columns,
                TheSecondOperand::ByTheColumnsOfTheResult {
                    values: &of_the_columns_with_one_that_is_not_finite,
                    cols: 2,
                },
            ),
        ] {
            // The call with nothing wrong gives the matrix, so the test
            // reaches the combination it is about and not only the check
            // of the values, which comes before the two are chosen.
            let mut c = vec![0.0; 4];
            product(
                TheFirstOperand::ByTheValuesSummedOver {
                    values: &A_OF_2_BY_3_THE_OTHER_WAY_ROUND,
                    rows: 2,
                },
                3,
                b,
                &mut c,
            )
            .unwrap();
            assert_eq!(c, THE_PRODUCT_THAT_IS_NOT_SYMMETRIC);

            for place in 0..A_OF_2_BY_3_THE_OTHER_WAY_ROUND.len() {
                let mut a = A_OF_2_BY_3_THE_OTHER_WAY_ROUND;
                a[place] = if place % 2 == 0 {
                    f64::NAN
                } else {
                    f64::INFINITY
                };
                let mut c = vec![0.0; 4];
                let error = product(
                    TheFirstOperand::ByTheValuesSummedOver {
                        values: &a,
                        rows: 2,
                    },
                    3,
                    b,
                    &mut c,
                )
                .unwrap_err();
                assert!(
                    matches!(error, Error::NotFinite { argument: "a" }),
                    "the error for the place {place} of {b:?} is {error}"
                );
            }
            let mut c = vec![0.0; 4];
            let error = product(
                TheFirstOperand::ByTheValuesSummedOver {
                    values: &A_OF_2_BY_3_THE_OTHER_WAY_ROUND,
                    rows: 2,
                },
                3,
                with_one_that_is_not_finite,
                &mut c,
            )
            .unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "b" }),
                "the error for {b:?} is {error}"
            );
        }
    }

    #[test]
    fn a_dimension_above_what_the_routines_count_in_is_refused() {
        // 2^31, one more than the largest an i32 holds. The check comes
        // before the one of the length of the buffer, so an empty slice
        // reaches it, and it is made whichever backend would run.
        let mut c: Vec<f64> = Vec::new();
        let error = product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &[],
                rows: 1 << 31,
            },
            1,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &[],
                cols: 1,
            },
            &mut c,
        )
        .unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error is {error}"
        );
        let error = product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &[],
                rows: 1 << 31,
            },
            1,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &[],
                cols: 1,
            },
            &mut c,
        )
        .unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error of the second operand held the other way is {error}"
        );
        let error = product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &[],
                rows: 1 << 31,
            },
            1,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &[],
                cols: 1,
            },
            &mut c,
        )
        .unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error of the first operand held the other way is {error}"
        );
        let error = product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &[],
                rows: 1 << 31,
            },
            1,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &[],
                cols: 1,
            },
            &mut c,
        )
        .unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error of the two operands held the other way is {error}"
        );
        // With no rows the a of the self product holds no values whatever
        // its number of columns is, and what the columns are too many for
        // is the cols x cols g.
        let mut g: Vec<f64> = Vec::new();
        let error = add_self_product_lower(&[], 0, 1 << 31, &mut g).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "g", .. }),
            "the error is {error}"
        );
        let error = eigh_lower(Vec::new(), 1 << 31).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "g", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn a_number_of_values_that_no_usize_counts_is_refused() {
        // 2^40 times 2^40 is 2^80, which no usize of any machine holds, so
        // this is the arm where the two dimensions do not multiply.
        let mut g: Vec<f64> = Vec::new();
        let error = add_self_product_lower(&[], 1 << 40, 1 << 40, &mut g).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
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
    fn the_eigendecomposition_refuses_a_g_that_is_not_n_by_n() {
        for held in [8_usize, 10] {
            let error = eigh_lower(vec![0.0; held], 3).unwrap_err();
            assert!(
                matches!(error, Error::Dimension { argument: "g", .. }),
                "the error for a g of {held} values is {error}"
            );
        }
    }

    #[test]
    fn the_eigendecomposition_refuses_a_value_that_is_not_finite_in_the_lower_half() {
        // Below the diagonal and on it.
        for entry in [3_usize, 4] {
            let mut g = the_matrix_of_3_by_3();
            *g.get_mut(entry).unwrap() = f64::NAN;
            let error = eigh_lower(g, 3).unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "g" }),
                "the error for the entry {entry} is {error}"
            );
        }
    }

    #[test]
    fn a_routine_that_stopped_says_which_of_the_three_things_happened() {
        let did_not_converge = Error::NoConvergence {
            routine: "dsyevd",
            info: 7,
        };
        assert_eq!(
            did_not_converge.to_string(),
            "dsyevd did not converge: it gave the info 7"
        );
        let refused_an_argument = Error::NoConvergence {
            routine: "dsyevd",
            info: -2,
        };
        assert_eq!(
            refused_an_argument.to_string(),
            "dsyevd refused its argument 2, a defect of popnei"
        );
        let faer = Error::NoConvergence {
            routine: "faer",
            info: 0,
        };
        assert_eq!(faer.to_string(), "faer did not converge");
    }

    #[test]
    fn a_workspace_that_could_not_be_allocated_says_what_it_was() {
        // No call of the crate reaches this case on a machine that has the
        // memory, so the message is checked on an error built by hand.
        let error = Error::Memory {
            what: "the workspace of floats of dsyevd",
            values: 200060001,
        };
        assert_eq!(
            error.to_string(),
            "this machine has not the memory for the workspace of floats of dsyevd, 200060001 values"
        );
    }

    #[test]
    fn the_cholesky_of_the_3_by_3_writes_the_lower_half_and_leaves_the_upper_as_it_was() {
        // The factorization of "How the seven are verified" has rows
        // (2, 0, 0), (1, 3, 0) and (0, 2, 1). Every entry is a small
        // integer, the square roots are of 4, 9 and 1 and the entries
        // below the diagonal are sums of at most two products of small
        // integers, so the arithmetic is exact and the assertion is too.
        let mut a = the_matrix_to_factor();
        cholesky_lower(&mut a, 3).unwrap();
        assert_eq!(
            a,
            vec![
                2.0, 99.0, 99.0, //
                1.0, 3.0, 99.0, //
                0.0, 2.0, 1.0,
            ]
        );
    }

    #[test]
    fn the_cholesky_of_a_matrix_that_is_not_positive_definite_is_singular_at_the_row_it_stopped_at()
    {
        // The 3 x 3 of "The errors the seven add", rows (4, 2, 0),
        // (2, 1, 0) and (0, 0, 1), whose leading 2 x 2 has a determinant
        // of 0: the row it stops at is the middle one of the three, so a
        // backend that counted from 1, as LAPACK does, or from the other
        // end gives another number.
        let mut a = vec![
            4.0, 99.0, 99.0, //
            2.0, 1.0, 99.0, //
            0.0, 0.0, 1.0,
        ];
        let error = cholesky_lower(&mut a, 3).unwrap_err();
        assert!(
            matches!(
                error,
                Error::Singular {
                    argument: "a",
                    at: 1
                }
            ),
            "the error is {error}"
        );
        assert_eq!(
            error.to_string(),
            "the matrix a failed at its row 1, counting from 0"
        );
    }

    #[test]
    fn the_cholesky_is_singular_at_the_first_row_and_at_the_last_when_that_is_where_it_stops() {
        // Every other fixture of the crate stops at the middle row of
        // three, so the row a backend gives could be a constant and pass
        // them all. The first: the 3 x 3 whose first diagonal entry is -4,
        // which no square root reaches, stops at the row 0. The last: the
        // 3 x 3 to factor with its last diagonal entry -5, where the rows
        // above it take 0 and 4 out of that -5 and leave -9, stops at the
        // row 2. Both were run on the two backends.
        let the_matrices_and_their_rows = [
            (
                vec![
                    -4.0, 99.0, 99.0, //
                    2.0, 10.0, 99.0, //
                    0.0, 6.0, 5.0,
                ],
                0,
            ),
            (
                vec![
                    4.0, 99.0, 99.0, //
                    2.0, 10.0, 99.0, //
                    0.0, 6.0, -5.0,
                ],
                2,
            ),
        ];
        for (mut a, row) in the_matrices_and_their_rows {
            let error = cholesky_lower(&mut a, 3).unwrap_err();
            assert!(
                matches!(error, Error::Singular { argument: "a", at } if at == row),
                "the error for the matrix that stops at the row {row} is {error}"
            );
        }
    }

    #[test]
    fn the_cholesky_of_the_1000_by_1000_matrix_of_the_generator() {
        let mut g = the_matrix_of_1000_by_1000_of_the_generator();
        let trace: f64 = g
            .as_chunks::<THE_LARGE_CASE>()
            .0
            .iter()
            .enumerate()
            .map(|(row, entries)| entries[row])
            .sum();
        assert!(
            !differ_in_their_digits(trace, 99996.3873081677, 1e-12),
            "the trace of the matrix to factor is {trace}"
        );

        cholesky_lower(&mut g, THE_LARGE_CASE).unwrap();
        let first = g[0];
        let last = g[THE_LARGE_CASE * THE_LARGE_CASE - 1];
        assert!(
            !differ_in_their_digits(first, 10.135944716832457, 1e-12),
            "the first entry of the diagonal of the factorization is {first}"
        );
        assert!(
            !differ_in_their_digits(last, 4.115426436421405, 1e-12),
            "the last entry of the diagonal of the factorization is {last}"
        );
    }

    #[test]
    fn the_cholesky_reads_and_writes_the_first_values_of_an_a_that_holds_more() {
        let mut a = the_matrix_to_factor();
        a.push(7.0);
        cholesky_lower(&mut a, 3).unwrap();
        assert_eq!(
            a,
            vec![
                2.0, 99.0, 99.0, //
                1.0, 3.0, 99.0, //
                0.0, 2.0, 1.0, //
                7.0,
            ]
        );
    }

    #[test]
    fn the_cholesky_refuses_an_n_of_zero() {
        let error = cholesky_lower(&mut [], 0).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "n", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_cholesky_refuses_an_a_shorter_than_n_times_n() {
        let mut a = [0.0; 8];
        let error = cholesky_lower(&mut a, 3).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_cholesky_refuses_a_dimension_above_what_the_routines_count_in() {
        // 2^31, one more than the largest an i32 holds. The check comes
        // before the one of the length of the buffer, so an empty slice
        // reaches it, and it is made whichever backend would run.
        let error = cholesky_lower(&mut [], 1 << 31).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_cholesky_refuses_a_value_that_is_not_finite_in_the_lower_half_alone() {
        // Below the diagonal and on it.
        for entry in [3_usize, 4] {
            let mut a = the_matrix_to_factor();
            a[entry] = f64::NAN;
            let error = cholesky_lower(&mut a, 3).unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "a" }),
                "the error for the entry {entry} is {error}"
            );
        }
        // And above it, where nothing is read: the factorization is the
        // one of the matrix and the value is left where it was, which is
        // what says that the upper half never reaches a routine.
        let mut a = the_matrix_to_factor();
        a[1] = f64::NAN;
        cholesky_lower(&mut a, 3).unwrap();
        let lower_half = vec![a[0], a[3], a[4], a[6], a[7], a[8]];
        assert_eq!(lower_half, vec![2.0, 1.0, 3.0, 0.0, 2.0, 1.0]);
        assert!(a[1].is_nan(), "the entry above the diagonal is {}", a[1]);
    }

    #[test]
    fn the_solve_with_the_cholesky_of_the_3_by_3_gives_the_one_right_hand_side() {
        // The right hand side (8, 40, 27) of "How the seven are verified"
        // and its solution (1, 2, 3). The upper half of the buffer still
        // holds the 99 the factorization left there, so a call that read
        // it instead of the lower half would give something else.
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        let mut b = THE_RIGHT_HAND_SIDE;
        solve_with_cholesky(&l, 3, &mut b, 1).unwrap();
        assert!(
            !differ(&b, &[1.0, 2.0, 3.0], 1e-14),
            "the solution is {b:?}"
        );
    }

    #[test]
    fn the_solve_with_the_cholesky_of_the_3_by_3_gives_the_two_right_hand_sides() {
        // The two right hand sides of "How the seven are verified", one
        // row each: (8, 40, 27) gives (1, 2, 3) and (4, 2, 0), which is
        // the first column of the matrix, gives (1, 0, 0). A backend that
        // read the rows of `b` as its columns gives something else for
        // both, and `sides` is 2 against an `n` of 3, so a call that swapped
        // the two dimensions would not fit the buffer either.
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        let mut b: Vec<f64> = THE_RIGHT_HAND_SIDE
            .into_iter()
            .chain(THE_SECOND_RIGHT_HAND_SIDE)
            .collect();
        solve_with_cholesky(&l, 3, &mut b, 2).unwrap();
        assert!(
            !differ(&b, &[1.0, 2.0, 3.0, 1.0, 0.0, 0.0], 1e-14),
            "the solutions are {b:?}"
        );
    }

    #[test]
    fn the_solve_with_the_cholesky_of_the_1000_by_1000_matrix_of_the_generator() {
        // The x of G x = v for v the vector of 1000 ones, of which "How
        // the seven are verified" gives the first three entries and the
        // sum over the 1000.
        let mut g = the_matrix_of_1000_by_1000_of_the_generator();
        cholesky_lower(&mut g, THE_LARGE_CASE).unwrap();
        let mut b = vec![1.0; THE_LARGE_CASE];
        solve_with_cholesky(&g, THE_LARGE_CASE, &mut b, 1).unwrap();
        let the_first_three = [
            -0.3054837936349659,
            -0.04576083734778211,
            -0.21314634692119025,
        ];
        for (entry, (got, expected)) in b.iter().zip(the_first_three).enumerate() {
            assert!(
                !differ_in_their_digits(*got, expected, 1e-11),
                "the entry {entry} of the solution is {got}"
            );
        }
        let total: f64 = b.iter().sum();
        assert!(
            !differ_in_their_digits(total, 73.9565335781636, 1e-11),
            "the sum of the solution is {total}"
        );
    }

    #[test]
    fn the_solve_with_the_cholesky_reads_the_lower_half_of_l_alone() {
        // A value that is nothing of the factorization above its diagonal
        // reaches neither the check nor the routine: the solution is the
        // one of the matrix.
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        l[1] = f64::NAN;
        let mut b = THE_RIGHT_HAND_SIDE;
        solve_with_cholesky(&l, 3, &mut b, 1).unwrap();
        assert!(
            !differ(&b, &[1.0, 2.0, 3.0], 1e-14),
            "the solution is {b:?}"
        );
    }

    #[test]
    fn the_solve_with_the_cholesky_reads_and_writes_the_first_values_of_buffers_that_hold_more() {
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        l.push(7.0);
        let mut b = vec![8.0, 40.0, 27.0, 7.0];
        solve_with_cholesky(&l, 3, &mut b, 1).unwrap();
        assert!(
            !differ(&b, &[1.0, 2.0, 3.0, 7.0], 1e-14),
            "the buffer of the solution is {b:?}"
        );
    }

    #[test]
    fn the_solve_with_the_cholesky_refuses_an_n_of_zero() {
        let error = solve_with_cholesky(&[], 0, &mut [], 1).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "n", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_solve_with_the_cholesky_refuses_a_sides_of_zero() {
        let error = solve_with_cholesky(&[0.0; 9], 3, &mut [], 0).unwrap_err();
        assert!(
            matches!(
                error,
                Error::Dimension {
                    argument: "sides",
                    ..
                }
            ),
            "the error is {error}"
        );
    }

    #[test]
    fn the_solve_with_the_cholesky_refuses_an_l_shorter_than_n_times_n() {
        let error = solve_with_cholesky(&[0.0; 8], 3, &mut [0.0; 3], 1).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "l", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_solve_with_the_cholesky_refuses_a_b_shorter_than_sides_times_n() {
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        let error = solve_with_cholesky(&l, 3, &mut [0.0; 5], 2).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "b", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_solve_with_the_cholesky_refuses_a_dimension_above_what_the_routines_count_in() {
        // 2^31, one more than the largest an i32 holds. The check comes
        // before the one of the length of the buffer, so empty slices
        // reach it, and it is made whichever backend would run.
        let error = solve_with_cholesky(&[], 1 << 31, &mut [], 1).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "l", .. }),
            "the error of the n is {error}"
        );
        let error = solve_with_cholesky(&[0.0; 9], 3, &mut [], 1 << 31).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "b", .. }),
            "the error of the sides is {error}"
        );
    }

    #[test]
    fn the_solve_with_the_cholesky_refuses_a_value_that_is_not_finite_in_the_lower_half_of_l() {
        // Below the diagonal and on it. Above it is the test that reads
        // the lower half alone.
        for entry in [3_usize, 4] {
            let mut l = the_matrix_to_factor();
            cholesky_lower(&mut l, 3).unwrap();
            l[entry] = f64::INFINITY;
            let mut b = THE_RIGHT_HAND_SIDE;
            let error = solve_with_cholesky(&l, 3, &mut b, 1).unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "l" }),
                "the error for the entry {entry} is {error}"
            );
        }
    }

    #[test]
    fn the_solve_with_the_cholesky_refuses_a_value_that_is_not_finite_in_b() {
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        let mut b = [8.0, f64::NAN, 27.0];
        let error = solve_with_cholesky(&l, 3, &mut b, 1).unwrap_err();
        assert!(
            matches!(error, Error::NotFinite { argument: "b" }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_log_determinant_with_the_cholesky_of_the_3_by_3_is_the_log_of_36() {
        // The diagonal of the factorization is 2, 3 and 1, and twice the
        // sum of their logs is the f64 numpy's slogdet gives for the
        // matrix, 3.58351893845611, which is the log of its determinant of
        // 36. It is not asserted to the bit, although it lands on that f64
        // here: `ln` is not rounded the same on every platform, and the
        // libm of wasm32-unknown-unknown already differs from this machine
        // on the log of 3. The spec's tolerance of 1e-15 relative leaves
        // about eight of the units in the last place of the answer.
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        let logarithm = log_determinant_with_cholesky(&l, 3).unwrap();
        assert!(
            !differ_in_their_digits(logarithm, 3.58351893845611, 1e-15),
            "the log of the determinant is {logarithm}"
        );
    }

    #[test]
    fn the_log_determinant_with_the_cholesky_reads_the_diagonal_alone() {
        // A value that is not finite off the diagonal, above it and below
        // it, changes neither the number nor the check: the operation
        // reads the diagonal and nothing else.
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        l[1] = f64::NAN;
        l[3] = f64::INFINITY;
        let logarithm = log_determinant_with_cholesky(&l, 3).unwrap();
        assert!(
            !differ_in_their_digits(logarithm, 3.58351893845611, 1e-15),
            "the log of the determinant is {logarithm}"
        );
    }

    #[test]
    fn the_log_determinant_with_the_cholesky_of_the_1000_by_1000_matrix_of_the_generator() {
        let mut g = the_matrix_of_1000_by_1000_of_the_generator();
        cholesky_lower(&mut g, THE_LARGE_CASE).unwrap();
        let logarithm = log_determinant_with_cholesky(&g, THE_LARGE_CASE).unwrap();
        assert!(
            !differ_in_their_digits(logarithm, 3963.7986384485084, 1e-13),
            "the log of the determinant is {logarithm}"
        );
    }

    #[test]
    fn the_log_determinant_of_a_diagonal_entry_that_is_not_above_zero_is_singular_at_that_row() {
        // A diagonal entry of 0 and one below 0, each at the middle row of
        // the three, which is the row a `cholesky_lower` that gave such an
        // `l` would have stopped at.
        for entry in [0.0, -3.0] {
            let l = vec![
                2.0, 99.0, 99.0, //
                1.0, entry, 99.0, //
                0.0, 2.0, 1.0,
            ];
            let error = log_determinant_with_cholesky(&l, 3).unwrap_err();
            assert!(
                matches!(
                    error,
                    Error::Singular {
                        argument: "l",
                        at: 1
                    }
                ),
                "the error for the diagonal entry {entry} is {error}"
            );
            assert_eq!(
                error.to_string(),
                "the matrix l failed at its row 1, counting from 0"
            );
        }
    }

    #[test]
    fn the_three_operations_given_the_cholesky_are_singular_at_the_row_of_the_diagonal_entry() {
        // A diagonal entry that is not above 0 at each of the three rows
        // in turn, so that the row the error names is the row that holds
        // it and not a constant. The solve, the log of the determinant and
        // the inverse read the same diagonal and give the same `Singular`,
        // which is what a `cholesky_lower` that had produced such an `l`
        // would have given at that row. Without the check both backends
        // divide by the entry: measured on 23 September 2026, `dpotrs`
        // gave NaN and two infinities and faer's solve three NaN, each
        // with no error, where numpy refuses the system.
        let the_factorizations_and_their_rows = [
            (
                vec![
                    -1.0, 99.0, 99.0, //
                    1.0, 3.0, 99.0, //
                    0.0, 2.0, 1.0,
                ],
                0,
            ),
            (
                vec![
                    2.0, 99.0, 99.0, //
                    1.0, 0.0, 99.0, //
                    0.0, 2.0, 1.0,
                ],
                1,
            ),
            (
                vec![
                    2.0, 99.0, 99.0, //
                    1.0, 3.0, 99.0, //
                    0.0, 2.0, -1.0,
                ],
                2,
            ),
        ];
        for (l, row) in the_factorizations_and_their_rows {
            let mut b = THE_RIGHT_HAND_SIDE;
            let error = solve_with_cholesky(&l, 3, &mut b, 1).unwrap_err();
            assert!(
                matches!(error, Error::Singular { argument: "l", at } if at == row),
                "the error of the solve at the row {row} is {error}"
            );
            // The right hand side is left as it was, the check coming
            // before the routine, so no NaN reaches the caller's buffer.
            assert_eq!(b.to_vec(), THE_RIGHT_HAND_SIDE.to_vec());

            let error = log_determinant_with_cholesky(&l, 3).unwrap_err();
            assert!(
                matches!(error, Error::Singular { argument: "l", at } if at == row),
                "the error of the log of the determinant at the row {row} is {error}"
            );

            let mut inverse = vec![THE_VALUE_THE_INVERSE_HELD; 9];
            let error = invert_with_cholesky(&l, 3, &mut inverse).unwrap_err();
            assert!(
                matches!(error, Error::Singular { argument: "l", at } if at == row),
                "the error of the inverse at the row {row} is {error}"
            );
            assert_eq!(inverse, vec![THE_VALUE_THE_INVERSE_HELD; 9]);
        }
    }

    #[test]
    fn the_log_determinant_refuses_a_value_that_is_not_finite_in_the_diagonal() {
        for entry in [f64::NAN, f64::INFINITY] {
            let l = vec![
                2.0, 99.0, 99.0, //
                1.0, 3.0, 99.0, //
                0.0, 2.0, entry,
            ];
            let error = log_determinant_with_cholesky(&l, 3).unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "l" }),
                "the error for the diagonal entry {entry} is {error}"
            );
        }
    }

    #[test]
    fn the_log_determinant_refuses_an_n_of_zero() {
        let error = log_determinant_with_cholesky(&[], 0).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "n", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_log_determinant_refuses_an_l_shorter_than_n_times_n() {
        let error = log_determinant_with_cholesky(&[0.0; 8], 3).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "l", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_log_determinant_refuses_a_dimension_above_what_the_routines_count_in() {
        // 2^31, one more than the largest an i32 holds. No backend runs
        // here, and the limit holds all the same, so that the four
        // operations of a factorization refuse the same call.
        let error = log_determinant_with_cholesky(&[], 1 << 31).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "l", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_log_determinant_reads_the_first_values_of_an_l_that_holds_more() {
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        l.push(7.0);
        let logarithm = log_determinant_with_cholesky(&l, 3).unwrap();
        assert!(
            !differ_in_their_digits(logarithm, 3.58351893845611, 1e-15),
            "the log of the determinant is {logarithm}"
        );
    }

    #[test]
    fn the_cholesky_of_a_1_by_1_is_the_square_root_of_its_one_entry() {
        // The smallest matrix the four operations take, which is one
        // number: a of 4 factors into an l of 2. The dimensions are all 1
        // here, so a backend that swapped two of them passes, and what the
        // case is for is that the loops over the rows and the halves
        // behave when there is one row and no entry off the diagonal.
        let mut a = vec![4.0];
        cholesky_lower(&mut a, 1).unwrap();
        assert_eq!(a, vec![2.0]);
    }

    #[test]
    fn the_solve_with_the_cholesky_of_a_1_by_1_divides_by_its_entry_twice() {
        // 4 x = 8 gives 2 and 4 x = 4 gives 1, the two right hand sides
        // one row each.
        let mut l = vec![4.0];
        cholesky_lower(&mut l, 1).unwrap();
        let mut b = vec![8.0, 4.0];
        solve_with_cholesky(&l, 1, &mut b, 2).unwrap();
        assert!(!differ(&b, &[2.0, 1.0], 1e-14), "the solutions are {b:?}");
    }

    #[test]
    fn the_log_determinant_with_the_cholesky_of_a_1_by_1_is_the_log_of_its_entry() {
        // The determinant of the 1 x 1 whose entry is 4 is 4, and its log
        // is 1.3862943611198906, which is what Python 3.14's `math.log`
        // gives for 4 on this machine on 23 September 2026.
        let mut l = vec![4.0];
        cholesky_lower(&mut l, 1).unwrap();
        let logarithm = log_determinant_with_cholesky(&l, 1).unwrap();
        assert!(
            !differ_in_their_digits(logarithm, 1.3862943611198906, 1e-15),
            "the log of the determinant is {logarithm}"
        );
    }

    #[test]
    fn inverting_with_the_cholesky_of_a_1_by_1_is_one_over_its_entry() {
        // One over 4 is 0.25, which every f64 holds exactly.
        let mut l = vec![4.0];
        cholesky_lower(&mut l, 1).unwrap();
        let mut inverse = vec![THE_VALUE_THE_INVERSE_HELD];
        invert_with_cholesky(&l, 1, &mut inverse).unwrap();
        assert_eq!(inverse, vec![0.25]);
    }
    /// The lower half of the inverse of the 3 x 3 of "How the seven are
    /// verified", row after row: 7/18, -5/18, 5/9, 1/3, -2/3 and 1. They
    /// are the exact values and not the ones numpy prints, whose `inv`
    /// goes through an LU and gives 0.38888888888888884 where 7/18 is
    /// 0.3888888888888889 and 0.9999999999999998 for the entry that is 1,
    /// up to 2 units in the last place away from what a Cholesky gives.
    const THE_LOWER_HALF_OF_THE_INVERSE_OF_THE_3_BY_3: [f64; 6] = [
        7.0 / 18.0,
        -5.0 / 18.0,
        5.0 / 9.0,
        1.0 / 3.0,
        -2.0 / 3.0,
        1.0,
    ];

    /// The value the buffer of the inverse holds before a call, which is
    /// nothing of the matrix and nothing of the factorization either, so
    /// that a call that read what it was given, or that wrote into the
    /// wrong buffer, gives something else.
    const THE_VALUE_THE_INVERSE_HELD: f64 = -7.0;

    /// The entries of the lower half of an `n` x `n` matrix held row after
    /// row, the entries of column `j` at most `i` of row `i`, in that
    /// order, which is the order the spec writes the six of the 3 x 3 in.
    fn the_lower_half_of(matrix: &[f64], n: usize) -> Vec<f64> {
        matrix
            .chunks_exact(n)
            .enumerate()
            .flat_map(|(row, entries)| entries[..=row].to_vec())
            .collect()
    }

    #[test]
    fn inverting_with_the_cholesky_of_the_3_by_3_writes_the_lower_half_and_leaves_the_upper_as_it_was()
     {
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        let mut inverse = vec![THE_VALUE_THE_INVERSE_HELD; 9];
        invert_with_cholesky(&l, 3, &mut inverse).unwrap();
        for (entry, (got, expected)) in the_lower_half_of(&inverse, 3)
            .into_iter()
            .zip(THE_LOWER_HALF_OF_THE_INVERSE_OF_THE_3_BY_3)
            .enumerate()
        {
            assert!(
                !differ_in_their_digits(got, expected, 1e-15),
                "the entry {entry} of the lower half of the inverse is {got}"
            );
        }
        let upper_half = vec![inverse[1], inverse[2], inverse[5]];
        assert_eq!(upper_half, vec![THE_VALUE_THE_INVERSE_HELD; 3]);
    }

    #[test]
    fn inverting_with_the_cholesky_leaves_the_factorization_as_it_was() {
        // `l` and `inverse` are two buffers and not one, so the caller
        // still holds the factorization after the call: a backend that
        // inverted in place, which `dpotri` does, would leave the inverse
        // in `l` instead.
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        let the_factorization = l.clone();
        let mut inverse = vec![THE_VALUE_THE_INVERSE_HELD; 9];
        invert_with_cholesky(&l, 3, &mut inverse).unwrap();
        assert_eq!(l, the_factorization);
    }

    #[test]
    fn inverting_with_the_cholesky_reads_the_lower_half_of_l_alone() {
        // A value that is nothing of the factorization above its diagonal
        // reaches neither the check nor the routine: the inverse is the
        // one of the matrix.
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        l[1] = f64::NAN;
        let mut inverse = vec![THE_VALUE_THE_INVERSE_HELD; 9];
        invert_with_cholesky(&l, 3, &mut inverse).unwrap();
        assert!(
            !differ(
                &the_lower_half_of(&inverse, 3),
                &THE_LOWER_HALF_OF_THE_INVERSE_OF_THE_3_BY_3,
                1e-15
            ),
            "the lower half of the inverse is {:?}",
            the_lower_half_of(&inverse, 3)
        );
    }

    #[test]
    fn inverting_with_the_cholesky_of_the_1000_by_1000_matrix_of_the_generator() {
        // The first and the last entries of the diagonal of the inverse
        // and its trace, which "How the seven are verified" gives within
        // 1e-11; the comparison here is relative, which is the stricter of
        // the two for the two entries, of about 0.06, and the same for the
        // trace of about 59.8. The two backends were measured at 5.4e-15
        // and 4.9e-15 relative away from these three numbers.
        let mut g = the_matrix_of_1000_by_1000_of_the_generator();
        cholesky_lower(&mut g, THE_LARGE_CASE).unwrap();
        let mut inverse = vec![0.0_f64; THE_LARGE_CASE * THE_LARGE_CASE];
        invert_with_cholesky(&g, THE_LARGE_CASE, &mut inverse).unwrap();
        let first = inverse[0];
        let last = inverse[THE_LARGE_CASE * THE_LARGE_CASE - 1];
        assert!(
            !differ_in_their_digits(first, 0.06230734831937398, 1e-11),
            "the first entry of the diagonal of the inverse is {first}"
        );
        assert!(
            !differ_in_their_digits(last, 0.059043258015696806, 1e-11),
            "the last entry of the diagonal of the inverse is {last}"
        );
        let trace: f64 = inverse
            .as_chunks::<THE_LARGE_CASE>()
            .0
            .iter()
            .enumerate()
            .map(|(row, entries)| entries[row])
            .sum();
        assert!(
            !differ_in_their_digits(trace, 59.78707893196584, 1e-11),
            "the trace of the inverse is {trace}"
        );
    }

    #[test]
    fn inverting_with_the_cholesky_reads_and_writes_the_first_values_of_buffers_that_hold_more() {
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        l.push(7.0);
        let mut inverse = vec![THE_VALUE_THE_INVERSE_HELD; 9];
        inverse.push(5.0);
        invert_with_cholesky(&l, 3, &mut inverse).unwrap();
        assert!(
            !differ(
                &the_lower_half_of(&inverse[..9], 3),
                &THE_LOWER_HALF_OF_THE_INVERSE_OF_THE_3_BY_3,
                1e-15
            ),
            "the lower half of the inverse is {:?}",
            the_lower_half_of(&inverse[..9], 3)
        );
        assert_eq!(inverse[9..], [5.0]);
    }

    #[test]
    fn inverting_with_a_diagonal_entry_of_the_cholesky_that_is_not_above_zero_is_singular_at_that_row()
     {
        // A diagonal entry of 0 and one below 0, each at the middle row of
        // the three, which is the row a `cholesky_lower` that gave such an
        // `l` would have stopped at. The two backends do not agree on this
        // `l`, which is why the crate reads the diagonal above them:
        // measured on 23 September 2026, `dpotri` gave an `info` of 2 for
        // the 0 and faer's inverse gave no error and wrote infinities and
        // NaN into the buffer.
        for entry in [0.0, -3.0] {
            let l = vec![
                2.0, 99.0, 99.0, //
                1.0, entry, 99.0, //
                0.0, 2.0, 1.0,
            ];
            let mut inverse = vec![THE_VALUE_THE_INVERSE_HELD; 9];
            let error = invert_with_cholesky(&l, 3, &mut inverse).unwrap_err();
            assert!(
                matches!(
                    error,
                    Error::Singular {
                        argument: "l",
                        at: 1
                    }
                ),
                "the error for the diagonal entry {entry} is {error}"
            );
            assert_eq!(inverse, vec![THE_VALUE_THE_INVERSE_HELD; 9]);
        }
    }

    #[test]
    fn inverting_with_the_cholesky_refuses_an_n_of_zero() {
        let error = invert_with_cholesky(&[], 0, &mut []).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "n", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn inverting_with_the_cholesky_refuses_an_l_shorter_than_n_times_n() {
        let error = invert_with_cholesky(&[0.0; 8], 3, &mut [0.0; 9]).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "l", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn inverting_with_the_cholesky_refuses_an_inverse_shorter_than_n_times_n() {
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        let error = invert_with_cholesky(&l, 3, &mut [0.0; 8]).unwrap_err();
        assert!(
            matches!(
                error,
                Error::Dimension {
                    argument: "inverse",
                    ..
                }
            ),
            "the error is {error}"
        );
    }

    #[test]
    fn inverting_with_the_cholesky_refuses_a_dimension_above_what_the_routines_count_in() {
        // 2^31, one more than the largest an i32 holds. The check comes
        // before the one of the length of the buffer, so empty slices
        // reach it, and it is made whichever backend would run.
        let error = invert_with_cholesky(&[], 1 << 31, &mut []).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "l", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn inverting_with_the_cholesky_refuses_a_value_that_is_not_finite_in_the_lower_half_of_l() {
        // Below the diagonal and on it. Above it is the test that reads
        // the lower half alone.
        for entry in [3_usize, 4] {
            let mut l = the_matrix_to_factor();
            cholesky_lower(&mut l, 3).unwrap();
            l[entry] = f64::INFINITY;
            let mut inverse = vec![THE_VALUE_THE_INVERSE_HELD; 9];
            let error = invert_with_cholesky(&l, 3, &mut inverse).unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "l" }),
                "the error for the entry {entry} is {error}"
            );
        }
    }

    /// The design of "How the seven are verified" of
    /// `docs/specs/linalg.md`, the 4 x 2 of an intercept and one
    /// covariate, rows (1, 1), (1, 2), (1, 3) and (1, 4), row after row.
    const THE_DESIGN_OF_4_BY_2: [f64; 8] = [
        1.0, 1.0, //
        1.0, 2.0, //
        1.0, 3.0, //
        1.0, 4.0,
    ];

    /// The `r` of that design, rows (2, 5) and (0, 2.23606797749979),
    /// whose second diagonal entry is the square root of 5, with the sign
    /// of each column taken so that the diagonal is positive.
    const THE_R_OF_THE_DESIGN: [f64; 4] = [
        2.0,
        5.0, //
        0.0,
        2.23606797749979,
    ];

    /// The `q` of that design with the same sign taken, row after row: its
    /// first column is (0.5, 0.5, 0.5, 0.5) and its second
    /// (-0.6708203932499368, -0.22360679774997894, 0.223606797749979,
    /// 0.6708203932499369), which is the covariate less its mean divided
    /// by the length of that.
    const THE_Q_OF_THE_DESIGN: [f64; 8] = [
        0.5,
        -0.6708203932499368, //
        0.5,
        -0.22360679774997894, //
        0.5,
        0.223606797749979, //
        0.5,
        0.6708203932499369,
    ];

    /// How far a factorization of that design may be from the two matrices
    /// above: Accelerate gave -2.0 for the first entry of `r` and faer
    /// -1.9999999999999998, which is one unit in the last place.
    const THE_TOLERANCE_OF_THE_THIN_QR: f64 = 1e-14;

    /// The factorization with the sign of each column of `q`, and of the
    /// row of `r` that goes with it, taken so that the diagonal entry of
    /// `r` in that row is positive. `cols` is the columns of both.
    ///
    /// The sign a backend gives is its own, and turning a column of `q`
    /// and the row of `r` round together leaves their product as it was,
    /// so this is what the spec compares after. On 23 September 2026
    /// LAPACK, faer and numpy 2.5.3 all gave the negative diagonal for the
    /// design above, and a backend that chose the other sign meets the
    /// same assertions through this.
    fn with_the_diagonal_of_r_positive(factorization: &ThinQr, cols: usize) -> ThinQr {
        let signs: Vec<f64> = factorization
            .r
            .chunks_exact(cols)
            .enumerate()
            .map(|(row, entries)| match entries.get(row) {
                Some(diagonal) if *diagonal < 0.0 => -1.0,
                _ => 1.0,
            })
            .collect();
        ThinQr {
            q: factorization
                .q
                .chunks_exact(cols)
                .flat_map(|row| row.iter().zip(&signs).map(|(entry, sign)| entry * sign))
                .collect(),
            r: factorization
                .r
                .chunks_exact(cols)
                .zip(&signs)
                .flat_map(|(row, sign)| row.iter().map(move |entry| entry * sign))
                .collect(),
        }
    }

    #[test]
    fn the_thin_qr_of_the_design_of_4_by_2_gives_the_r_and_the_q_of_the_spec() {
        let factorization = thin_qr(&THE_DESIGN_OF_4_BY_2, 4, 2).unwrap();
        let factorization = with_the_diagonal_of_r_positive(&factorization, 2);
        assert!(
            !differ(
                &factorization.r,
                &THE_R_OF_THE_DESIGN,
                THE_TOLERANCE_OF_THE_THIN_QR
            ),
            "the r is {r:?}",
            r = factorization.r
        );
        assert!(
            !differ(
                &factorization.q,
                &THE_Q_OF_THE_DESIGN,
                THE_TOLERANCE_OF_THE_THIN_QR
            ),
            "the q is {q:?}",
            q = factorization.q
        );
    }

    #[test]
    fn the_thin_qr_writes_the_lower_half_of_r_as_zero() {
        // What the routines of LAPACK leave below the diagonal there is
        // the vectors the factorization is built from, and what faer
        // leaves is its own, so the crate writes that half itself. It is
        // read before the sign is fixed, which multiplies by 1 or -1 and
        // would leave a 0 where it found one.
        let factorization = thin_qr(&THE_DESIGN_OF_4_BY_2, 4, 2).unwrap();
        let the_lower_half: Vec<f64> = factorization
            .r
            .as_chunks::<2>()
            .0
            .iter()
            .enumerate()
            .flat_map(|(row, entries)| entries.iter().take(row).copied())
            .collect();
        assert_eq!(
            the_lower_half,
            vec![0.0],
            "the r is {r:?}",
            r = factorization.r
        );
    }

    #[test]
    fn the_thin_qr_of_the_design_multiplies_back_to_the_design() {
        // `a = q r` is what the factorization is for, and it holds
        // whichever sign the backend chose, so this is asserted on what
        // the backend gave. The two matrices have different dimensions, 4
        // x 2 and 2 x 2, so a backend that wrote one where the other goes
        // fails here as well as above.
        let factorization = thin_qr(&THE_DESIGN_OF_4_BY_2, 4, 2).unwrap();
        let mut design = vec![0.0_f64; 8];
        product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &factorization.q,
                rows: 4,
            },
            2,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &factorization.r,
                cols: 2,
            },
            &mut design,
        )
        .unwrap();
        assert!(
            !differ(&design, &THE_DESIGN_OF_4_BY_2, THE_TOLERANCE_OF_THE_THIN_QR),
            "q r is {design:?}"
        );
    }

    #[test]
    fn the_thin_qr_reads_the_first_values_of_an_a_that_holds_more() {
        let mut a = THE_DESIGN_OF_4_BY_2.to_vec();
        a.push(7.0);
        let factorization = thin_qr(&a, 4, 2).unwrap();
        let factorization = with_the_diagonal_of_r_positive(&factorization, 2);
        assert!(
            !differ(
                &factorization.r,
                &THE_R_OF_THE_DESIGN,
                THE_TOLERANCE_OF_THE_THIN_QR
            ),
            "the r is {r:?}",
            r = factorization.r
        );
        assert!(
            !differ(
                &factorization.q,
                &THE_Q_OF_THE_DESIGN,
                THE_TOLERANCE_OF_THE_THIN_QR
            ),
            "the q is {q:?}",
            q = factorization.q
        );
    }

    #[test]
    fn the_thin_qr_refuses_a_cols_of_zero() {
        let error = thin_qr(&[], 4, 0).unwrap_err();
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
    fn the_thin_qr_refuses_an_a_of_fewer_rows_than_columns() {
        // The 2 x 4 written where the 4 x 2 goes, which holds the eight
        // values either way and which no other check would catch.
        let error = thin_qr(&THE_DESIGN_OF_4_BY_2, 2, 4).unwrap_err();
        assert!(
            matches!(
                error,
                Error::Dimension {
                    argument: "rows",
                    ..
                }
            ),
            "the error is {error}"
        );
    }

    #[test]
    fn the_thin_qr_refuses_an_a_shorter_than_rows_times_cols() {
        let a = [0.0; 7];
        let error = thin_qr(&a, 4, 2).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_thin_qr_refuses_a_dimension_above_what_the_routines_count_in() {
        // 2^31 rows of one column, one value more than the largest an i32
        // holds. The check comes before the one of the length of the
        // buffer, so an empty slice reaches it, and it is made whichever
        // backend would run.
        let error = thin_qr(&[], 1 << 31, 1).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_thin_qr_refuses_a_value_that_is_not_finite_anywhere_in_a() {
        // The whole of a design is read, both halves, so the entry of the
        // first row and second column, which a check of the lower half
        // alone would walk past, is refused as the last entry is.
        for entry in [1_usize, 7] {
            let mut a = THE_DESIGN_OF_4_BY_2;
            a[entry] = f64::NAN;
            let error = thin_qr(&a, 4, 2).unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "a" }),
                "the error for the entry {entry} is {error}"
            );
        }
    }

    /// The three right hand sides of "How the seven are verified" of
    /// `docs/specs/linalg.md`, one row each: the `q' y` of the traits
    /// (1, 3, 5, 7), (4, 7, 10, 13) and (1, 2, 3, 4), which are twice the
    /// covariate less 1, three times it plus 1, and the covariate itself.
    const THE_RIGHT_HAND_SIDES_OF_THE_THREE_FITS: [f64; 6] = [
        8.0,
        4.472_135_954_999_58, //
        17.0,
        6.708_203_932_499_37, //
        5.0,
        2.236_067_977_499_79,
    ];

    /// The coefficients of those three fits, (-1, 2), (1, 3) and (0, 1),
    /// one row each, as `b` comes back holding them. Every fit is exact,
    /// so a backend that read `r` the wrong way round gives something
    /// else.
    const THE_COEFFICIENTS_OF_THE_THREE_FITS: [f64; 6] = [
        -1.0, 2.0, //
        1.0, 3.0, //
        0.0, 1.0,
    ];

    /// How far the coefficients may be from those: numpy 2.5.3 gives
    /// 0.9999999999999991 and 3.0000000000000004 for the second fit, 9e-16
    /// relative away from the exact one.
    const THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE: f64 = 1e-14;

    #[test]
    fn the_triangular_solve_of_the_r_of_the_design_gives_the_coefficients_of_the_fit() {
        let mut coefficients = [0.0_f64; 2];
        coefficients.copy_from_slice(&THE_RIGHT_HAND_SIDES_OF_THE_THREE_FITS[..2]);
        solve_triangular(
            &THE_R_OF_THE_DESIGN,
            2,
            TheHalfThatHoldsTheMatrix::TheUpperHalf,
            &mut coefficients,
            1,
        )
        .unwrap();
        assert!(
            !differ(
                &coefficients,
                &THE_COEFFICIENTS_OF_THE_THREE_FITS[..2],
                THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE
            ),
            "the coefficients are {coefficients:?}"
        );
    }

    #[test]
    fn the_triangular_solve_of_the_r_of_the_design_gives_the_coefficients_of_three_fits() {
        // Three right hand sides against an `r` of 2 x 2: the three rows
        // tell `sides` from `n`, which one right hand side of two numbers
        // cannot, and they catch a backend that read the rows of `b` as
        // its columns, which would solve three other systems here.
        let mut coefficients = THE_RIGHT_HAND_SIDES_OF_THE_THREE_FITS;
        solve_triangular(
            &THE_R_OF_THE_DESIGN,
            2,
            TheHalfThatHoldsTheMatrix::TheUpperHalf,
            &mut coefficients,
            3,
        )
        .unwrap();
        assert!(
            !differ(
                &coefficients,
                &THE_COEFFICIENTS_OF_THE_THREE_FITS,
                THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE
            ),
            "the coefficients are {coefficients:?}"
        );
    }

    /// The same `r` with the sign the backends give it, rows (-2, -5) and
    /// (0, -2.23606797749979), and the first right hand side with it,
    /// (-8, -4.47213595499958). The test of the thin QR takes the sign of
    /// each column of `q` so that the diagonal of `r` is positive, and
    /// "How the seven are verified" of `docs/specs/linalg.md` records that
    /// Accelerate gave -2.0 for the first entry of that `r`: what a caller
    /// which solves against an `r` as [`thin_qr`] gives it holds is this,
    /// a diagonal below 0. Turning `r` and `q' y` round together leaves
    /// the coefficients the fit has, so these are the literals of the fit
    /// above with their signs turned and no new number.
    const THE_R_OF_THE_DESIGN_WITH_THE_SIGN_THE_BACKENDS_GIVE: [f64; 4] = [
        -2.0,
        -5.0, //
        0.0,
        -2.23606797749979,
    ];

    /// That right hand side, the `q' y` of the trait (1, 3, 5, 7) against
    /// that `q`, one row.
    const THE_RIGHT_HAND_SIDE_WITH_THE_SIGN_THE_BACKENDS_GIVE: [f64; 2] = [-8.0, -4.47213595499958];

    #[test]
    fn the_triangular_solve_of_an_r_whose_diagonal_is_negative_gives_the_coefficients() {
        // The diagonal is read for a 0 and for nothing else, and this is
        // the case that says so: an `r` as a backend gives it, whose
        // diagonal is negative, and the coefficients (-1, 2) all the same.
        // A check that refused an entry at most 0, which is what the three
        // operations that read the `l` of a Cholesky use, would refuse
        // every least squares fit here with `Singular` and no test of the
        // twelve would have failed.
        let mut coefficients = THE_RIGHT_HAND_SIDE_WITH_THE_SIGN_THE_BACKENDS_GIVE;
        solve_triangular(
            &THE_R_OF_THE_DESIGN_WITH_THE_SIGN_THE_BACKENDS_GIVE,
            2,
            TheHalfThatHoldsTheMatrix::TheUpperHalf,
            &mut coefficients,
            1,
        )
        .unwrap();
        assert!(
            !differ(
                &coefficients,
                &THE_COEFFICIENTS_OF_THE_THREE_FITS[..2],
                THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE
            ),
            "the coefficients are {coefficients:?}"
        );
    }

    #[test]
    fn the_triangular_solve_reads_the_upper_half_of_a_alone() {
        // A NaN below the diagonal is neither refused nor read: the check
        // for a value that is not finite walks the upper half alone, and
        // the answer is the one the upper half gives.
        let mut a = THE_R_OF_THE_DESIGN;
        a[2] = f64::NAN;
        let mut coefficients = [0.0_f64; 2];
        coefficients.copy_from_slice(&THE_RIGHT_HAND_SIDES_OF_THE_THREE_FITS[..2]);
        solve_triangular(
            &a,
            2,
            TheHalfThatHoldsTheMatrix::TheUpperHalf,
            &mut coefficients,
            1,
        )
        .unwrap();
        assert!(
            !differ(
                &coefficients,
                &THE_COEFFICIENTS_OF_THE_THREE_FITS[..2],
                THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE
            ),
            "the coefficients are {coefficients:?}"
        );
    }

    #[test]
    fn the_triangular_solve_reads_and_writes_the_first_values_of_buffers_that_hold_more() {
        let mut a = THE_R_OF_THE_DESIGN.to_vec();
        a.push(7.0);
        let mut b = THE_RIGHT_HAND_SIDES_OF_THE_THREE_FITS[..2].to_vec();
        b.push(9.0);
        solve_triangular(&a, 2, TheHalfThatHoldsTheMatrix::TheUpperHalf, &mut b, 1).unwrap();
        assert!(
            !differ(
                &b[..2],
                &THE_COEFFICIENTS_OF_THE_THREE_FITS[..2],
                THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE
            ),
            "the solution is {b:?}"
        );
        assert_eq!(
            b.get(2),
            Some(&9.0),
            "the value after the right hand sides was written"
        );
    }

    #[test]
    fn the_triangular_solve_of_an_upper_half_with_a_zero_in_its_diagonal_is_singular_there() {
        // The `r` with rows (2, 5) and (0, 0) of "How the seven are
        // verified", which faer would divide by and answer a NaN or an
        // infinity for, by the half and the row, and the same 0 moved to
        // the first row, so that the row the error names is read and is
        // not the last row of the matrix.
        for (row, a) in [(1_usize, [2.0, 5.0, 0.0, 0.0]), (0, [0.0, 5.0, 0.0, 2.0])] {
            let mut b = [0.0_f64; 2];
            b.copy_from_slice(&THE_RIGHT_HAND_SIDES_OF_THE_THREE_FITS[..2]);
            let error = solve_triangular(&a, 2, TheHalfThatHoldsTheMatrix::TheUpperHalf, &mut b, 1)
                .unwrap_err();
            assert!(
                matches!(error, Error::Singular { argument: "a", at } if at == row),
                "the error for the row {row} is {error}"
            );
        }
    }

    #[test]
    fn the_triangular_solve_refuses_an_n_of_zero() {
        let mut b = [0.0_f64; 0];
        let error = solve_triangular(&[], 0, TheHalfThatHoldsTheMatrix::TheUpperHalf, &mut b, 1)
            .unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "n", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_triangular_solve_refuses_a_sides_of_zero() {
        let mut b = [0.0_f64; 0];
        let error = solve_triangular(
            &THE_R_OF_THE_DESIGN,
            2,
            TheHalfThatHoldsTheMatrix::TheUpperHalf,
            &mut b,
            0,
        )
        .unwrap_err();
        assert!(
            matches!(
                error,
                Error::Dimension {
                    argument: "sides",
                    ..
                }
            ),
            "the error is {error}"
        );
    }

    #[test]
    fn the_triangular_solve_refuses_an_a_shorter_than_n_times_n() {
        let a = [0.0_f64; 3];
        let mut b = [0.0_f64; 2];
        let error = solve_triangular(&a, 2, TheHalfThatHoldsTheMatrix::TheUpperHalf, &mut b, 1)
            .unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_triangular_solve_refuses_a_b_shorter_than_sides_times_n() {
        let mut b = [0.0_f64; 5];
        let error = solve_triangular(
            &THE_R_OF_THE_DESIGN,
            2,
            TheHalfThatHoldsTheMatrix::TheUpperHalf,
            &mut b,
            3,
        )
        .unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "b", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_triangular_solve_refuses_a_dimension_above_what_the_routines_count_in() {
        // 46341 rows of 46341 is 2147488281 values, more than the routines
        // count in. The check comes before the one of the length of the
        // buffer, so an empty slice reaches it, and it is made whichever
        // backend would run.
        let mut b = [0.0_f64; 0];
        let error = solve_triangular(
            &[],
            46341,
            TheHalfThatHoldsTheMatrix::TheUpperHalf,
            &mut b,
            1,
        )
        .unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_triangular_solve_refuses_a_value_that_is_not_finite_in_the_upper_half_of_a() {
        // The two entries of the diagonal and the one above it, which are
        // the three places of a 2 x 2 that the solve reads.
        for entry in [0_usize, 1, 3] {
            let mut a = THE_R_OF_THE_DESIGN;
            a[entry] = f64::INFINITY;
            let mut b = [0.0_f64; 2];
            b.copy_from_slice(&THE_RIGHT_HAND_SIDES_OF_THE_THREE_FITS[..2]);
            let error = solve_triangular(&a, 2, TheHalfThatHoldsTheMatrix::TheUpperHalf, &mut b, 1)
                .unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "a" }),
                "the error for the entry {entry} is {error}"
            );
        }
    }

    #[test]
    fn the_triangular_solve_refuses_a_value_that_is_not_finite_anywhere_in_b() {
        // The whole of `b` is read, so a value that is not finite is
        // refused wherever it sits among the three right hand sides.
        for entry in 0..6 {
            let mut b = THE_RIGHT_HAND_SIDES_OF_THE_THREE_FITS;
            b[entry] = f64::NAN;
            let error = solve_triangular(
                &THE_R_OF_THE_DESIGN,
                2,
                TheHalfThatHoldsTheMatrix::TheUpperHalf,
                &mut b,
                3,
            )
            .unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "b" }),
                "the error for the entry {entry} is {error}"
            );
        }
    }

    /// The `l` of the Cholesky of "How the seven are verified" of
    /// `docs/specs/linalg.md`, rows (2, 0, 0), (1, 3, 0) and (0, 2, 1),
    /// row after row, with the upper half of the buffer holding the mirror
    /// of the lower one, values that are nothing of that `l`: a call that
    /// read the upper half instead solves against the rows (2, 1, 0),
    /// (0, 3, 2) and (0, 0, 1), which is another system, and
    /// `the_triangular_solve_reads_the_half_the_caller_named` asserts what
    /// each of the two gives.
    const THE_L_WHOSE_UPPER_HALF_IS_NOTHING_OF_IT: [f64; 9] = [
        2.0, 1.0, 0.0, //
        1.0, 3.0, 2.0, //
        0.0, 2.0, 1.0,
    ];

    /// The two right hand sides of that solve, (8, 40, 27) and (4, 2, 0),
    /// one row each, which are the two the solve with the Cholesky is
    /// checked on.
    const THE_RIGHT_HAND_SIDES_OF_THE_LOWER_HALF: [f64; 6] = [
        8.0, 40.0, 27.0, //
        4.0, 2.0, 0.0,
    ];

    /// What the lower half gives for them, (4, 12, 3) and (2, 0, 0), one
    /// row each. Every entry is a small whole number, and they are
    /// compared within [`THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE`] all the
    /// same and not to the bit, which is what "How the seven are
    /// verified" of the spec asks for, because the two backends do not
    /// give the same bits here. Measured on 23 September 2026, `dtrtrs`
    /// gives exactly 12 and 3 for the second and the third entries and
    /// faer gives 11.999999999999998 and 3.0000000000000036, 1.5e-16 and
    /// 1.2e-15 relative away: faer scales each term of a row by the
    /// reciprocal of the diagonal entry and adds them, where the routine
    /// subtracts first and divides once, and `40 * (1/3) - 4 * (1/3)` is
    /// 11.999999999999998 where `(40 - 4) * (1/3)` is 12 exactly.
    const THE_SOLUTIONS_OF_THE_LOWER_HALF: [f64; 6] = [
        4.0, 12.0, 3.0, //
        2.0, 0.0, 0.0,
    ];

    /// What the upper half of the same buffer gives for the first of those
    /// right hand sides: the answer a caller that named the wrong half
    /// would come back with, from numpy 2.5.3 on 23 September 2026.
    const THE_SOLUTION_THE_UPPER_HALF_OF_IT_GIVES: [f64; 3] =
        [6.333333333333333, -4.666666666666666, 27.0];

    #[test]
    fn the_triangular_solve_against_the_lower_half_gives_the_solution() {
        let mut b = [0.0_f64; 3];
        b.copy_from_slice(&THE_RIGHT_HAND_SIDES_OF_THE_LOWER_HALF[..3]);
        solve_triangular(
            &THE_L_WHOSE_UPPER_HALF_IS_NOTHING_OF_IT,
            3,
            TheHalfThatHoldsTheMatrix::TheLowerHalf,
            &mut b,
            1,
        )
        .unwrap();
        assert!(
            !differ(
                &b,
                &THE_SOLUTIONS_OF_THE_LOWER_HALF[..3],
                THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE
            ),
            "the solution of the lower half is {b:?}"
        );
    }

    #[test]
    fn the_triangular_solve_against_the_lower_half_gives_the_solutions_of_two_right_hand_sides() {
        // Two right hand sides against an `l` of 3 x 3: `sides` is 2 and
        // `n` is 3, so a call that read the one dimension for the other
        // would not come back with these, and the two rows catch a backend
        // that read the rows of `b` as its columns.
        let mut b = THE_RIGHT_HAND_SIDES_OF_THE_LOWER_HALF;
        solve_triangular(
            &THE_L_WHOSE_UPPER_HALF_IS_NOTHING_OF_IT,
            3,
            TheHalfThatHoldsTheMatrix::TheLowerHalf,
            &mut b,
            2,
        )
        .unwrap();
        assert!(
            !differ(
                &b,
                &THE_SOLUTIONS_OF_THE_LOWER_HALF,
                THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE
            ),
            "the solutions of the lower half are {b:?}"
        );
    }

    #[test]
    fn the_triangular_solve_reads_the_half_the_caller_named() {
        // The same buffer and the same right hand side read as the one
        // half and as the other: the lower half gives (4, 12, 3) and the
        // upper (6.333333333333333, -4.666666666666666, 27), both of them
        // answers a caller could believe. Asserting the two is what says
        // the half is read at all, since no length tells them apart.
        let mut against_the_lower_half = [0.0_f64; 3];
        against_the_lower_half.copy_from_slice(&THE_RIGHT_HAND_SIDES_OF_THE_LOWER_HALF[..3]);
        let mut against_the_upper_half = against_the_lower_half;
        solve_triangular(
            &THE_L_WHOSE_UPPER_HALF_IS_NOTHING_OF_IT,
            3,
            TheHalfThatHoldsTheMatrix::TheLowerHalf,
            &mut against_the_lower_half,
            1,
        )
        .unwrap();
        solve_triangular(
            &THE_L_WHOSE_UPPER_HALF_IS_NOTHING_OF_IT,
            3,
            TheHalfThatHoldsTheMatrix::TheUpperHalf,
            &mut against_the_upper_half,
            1,
        )
        .unwrap();
        assert!(
            !differ(
                &against_the_lower_half,
                &THE_SOLUTIONS_OF_THE_LOWER_HALF[..3],
                THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE
            ),
            "the solution of the lower half is {against_the_lower_half:?}"
        );
        assert!(
            !differ(
                &against_the_upper_half,
                &THE_SOLUTION_THE_UPPER_HALF_OF_IT_GIVES,
                THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE
            ),
            "the solution of the upper half is {against_the_upper_half:?}"
        );
    }

    #[test]
    fn the_triangular_solve_of_the_lower_half_gives_numbers_that_are_not_whole() {
        // The third right hand side of "How the seven are verified",
        // (2, 6, 1), whose solution is the one of the three that is not a
        // triple of whole numbers, so it is asserted within the tolerance
        // and not exactly.
        let mut b = [2.0_f64, 6.0, 1.0];
        solve_triangular(
            &THE_L_WHOSE_UPPER_HALF_IS_NOTHING_OF_IT,
            3,
            TheHalfThatHoldsTheMatrix::TheLowerHalf,
            &mut b,
            1,
        )
        .unwrap();
        assert!(
            !differ(
                &b,
                &[1.0, 1.6666666666666665, -2.333333333333333],
                THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE
            ),
            "the solution of the lower half is {b:?}"
        );
    }

    #[test]
    fn the_triangular_solve_reads_the_lower_half_of_a_alone() {
        // A NaN above the diagonal is neither refused nor read when the
        // lower half is the one named: the check for a value that is not
        // finite walks the half the caller named, and the answer is the one
        // the lower half gives.
        let mut a = THE_L_WHOSE_UPPER_HALF_IS_NOTHING_OF_IT;
        a[1] = f64::NAN;
        let mut b = [0.0_f64; 3];
        b.copy_from_slice(&THE_RIGHT_HAND_SIDES_OF_THE_LOWER_HALF[..3]);
        solve_triangular(&a, 3, TheHalfThatHoldsTheMatrix::TheLowerHalf, &mut b, 1).unwrap();
        assert!(
            !differ(
                &b,
                &THE_SOLUTIONS_OF_THE_LOWER_HALF[..3],
                THE_TOLERANCE_OF_THE_TRIANGULAR_SOLVE
            ),
            "the solution of the lower half is {b:?}"
        );
    }

    #[test]
    fn the_triangular_solve_refuses_a_value_that_is_not_finite_in_the_lower_half_of_a() {
        // Every one of the six places of a 3 x 3 that the solve against
        // the lower half reads: the three entries of the diagonal and the
        // three below it, the last row among them.
        for entry in [0_usize, 3, 4, 6, 7, 8] {
            let mut a = THE_L_WHOSE_UPPER_HALF_IS_NOTHING_OF_IT;
            a[entry] = f64::INFINITY;
            let mut b = [0.0_f64; 3];
            b.copy_from_slice(&THE_RIGHT_HAND_SIDES_OF_THE_LOWER_HALF[..3]);
            let error = solve_triangular(&a, 3, TheHalfThatHoldsTheMatrix::TheLowerHalf, &mut b, 1)
                .unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "a" }),
                "the error for the entry {entry} is {error}"
            );
        }
    }

    #[test]
    fn the_triangular_solve_of_a_lower_half_with_a_zero_in_its_diagonal_is_singular_at_that_row() {
        // The `l` with rows (2, 0) and (5, 0) of "How the seven are
        // verified", which faer would divide by and answer a NaN or an
        // infinity for, by the half and the row, and the same 0 moved to
        // the first row, so that the row the error names is read and is
        // not the last row of the matrix. The diagonal is read for either
        // half, since it belongs to both.
        //
        // What guards the crate's own check here is the faer run, `cargo
        // test -p popnei-linalg --no-default-features`: `dtrtrs` gives an
        // `info` for this case of its own accord and the BLAS backend maps
        // it to this same error, so on that backend the test passes with
        // `refuse_a_diagonal_entry` gone. The review of this work package
        // measured both on 23 September 2026.
        for (row, a) in [(1_usize, [2.0, 0.0, 5.0, 0.0]), (0, [0.0, 0.0, 5.0, 2.0])] {
            let mut b = [8.0_f64, 40.0];
            let error = solve_triangular(&a, 2, TheHalfThatHoldsTheMatrix::TheLowerHalf, &mut b, 1)
                .unwrap_err();
            assert!(
                matches!(error, Error::Singular { argument: "a", at } if at == row),
                "the error for the row {row} is {error}"
            );
        }
    }

    /// The 4 x 3 of "How the seven are verified" of
    /// `docs/specs/linalg.md`, rows (1, 1, 2), (1, 2, 3), (1, 3, 4) and
    /// (1, 4, 5), whose third column is the sum of the first two, row
    /// after row. Its singular values are 9.344132686098556,
    /// 0.8289658283575813 and 3.651382431893325e-17, against a tolerance
    /// of 8.299257002607302e-15.
    const THE_DESIGN_OF_4_BY_3_OF_A_COLUMN_THAT_REPEATS: [f64; 12] = [
        1.0, 1.0, 2.0, //
        1.0, 2.0, 3.0, //
        1.0, 3.0, 4.0, //
        1.0, 4.0, 5.0,
    ];

    /// The 4 x 2 of the same place whose covariate is the constant 5, row
    /// after row. Its singular values are 10.19803902718557 and 0.
    const THE_DESIGN_OF_A_CONSTANT_COVARIATE: [f64; 8] = [
        1.0, 5.0, //
        1.0, 5.0, //
        1.0, 5.0, //
        1.0, 5.0,
    ];

    #[test]
    fn the_rank_of_the_design_of_an_intercept_and_one_covariate_is_2() {
        assert_eq!(rank(&THE_DESIGN_OF_4_BY_2, 4, 2).unwrap(), 2);
    }

    #[test]
    fn the_rank_of_a_design_whose_third_column_is_the_sum_of_the_first_two_is_2() {
        // A 4 x 3 and not a square matrix, so the two dimensions are told
        // apart: the same twelve values read as the 3 x 4 have rank 3,
        // from numpy 2.5.3 on 23 September 2026.
        assert_eq!(
            rank(&THE_DESIGN_OF_4_BY_3_OF_A_COLUMN_THAT_REPEATS, 4, 3).unwrap(),
            2
        );
    }

    #[test]
    fn the_rank_of_a_design_whose_covariate_is_constant_is_1() {
        assert_eq!(rank(&THE_DESIGN_OF_A_CONSTANT_COVARIATE, 4, 2).unwrap(), 1);
    }

    #[test]
    fn the_rank_counts_a_singular_value_above_the_tolerance_and_not_one_below_it() {
        // The pair that pins the tolerance itself, which for a 2 x 2 whose
        // largest singular value is 1 is 4.440892098500626e-16: the three
        // designs above give their counts at a wrong threshold too, and
        // these two do not. numpy 2.5.3 gave 2 and 1 for them on 23
        // September 2026.
        for (entry, wanted) in [(5e-16, 2_usize), (4e-16, 1)] {
            let a = [1.0, 0.0, 0.0, entry];
            assert_eq!(
                rank(&a, 2, 2).unwrap(),
                wanted,
                "the rank of the 2 x 2 whose second singular value is {entry}"
            );
        }
    }

    #[test]
    fn the_rank_of_a_matrix_of_zeros_is_0() {
        // The largest singular value is 0, so the tolerance is 0 and no
        // value is strictly above it. numpy 2.5.3 gives 0 for the same
        // matrix.
        let a = [0.0_f64; 6];
        assert_eq!(rank(&a, 3, 2).unwrap(), 0);
    }

    #[test]
    fn the_rank_reads_the_first_values_of_an_a_that_holds_more() {
        let mut a = THE_DESIGN_OF_A_CONSTANT_COVARIATE.to_vec();
        // A ninth value that would make the covariate no longer constant
        // if it were read.
        a.push(7.0);
        assert_eq!(rank(&a, 4, 2).unwrap(), 1);
    }

    #[test]
    fn the_rank_refuses_a_rows_of_zero() {
        let error = rank(&[], 0, 2).unwrap_err();
        assert!(
            matches!(
                error,
                Error::Dimension {
                    argument: "rows",
                    ..
                }
            ),
            "the error is {error}"
        );
    }

    #[test]
    fn the_rank_refuses_a_cols_of_zero() {
        let error = rank(&[], 4, 0).unwrap_err();
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
    fn the_rank_refuses_an_a_shorter_than_rows_times_cols() {
        let a = [0.0_f64; 11];
        let error = rank(&a, 4, 3).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_rank_refuses_a_dimension_above_what_the_routines_count_in() {
        // 2^31 rows of one column, one value more than the largest an i32
        // holds. The check comes before the one of the length of the
        // buffer, so an empty slice reaches it, and it is made whichever
        // backend would run.
        let error = rank(&[], 1 << 31, 1).unwrap_err();
        assert!(
            matches!(error, Error::Dimension { argument: "a", .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_rank_refuses_a_value_that_is_not_finite_anywhere_in_a() {
        // The whole of the matrix is read, so the first value and the last
        // are both refused.
        for entry in [0_usize, 11] {
            let mut a = THE_DESIGN_OF_4_BY_3_OF_A_COLUMN_THAT_REPEATS;
            a[entry] = f64::NAN;
            let error = rank(&a, 4, 3).unwrap_err();
            assert!(
                matches!(error, Error::NotFinite { argument: "a" }),
                "the error for the entry {entry} is {error}"
            );
        }
    }

    /// The second design of "How the seven are verified" of
    /// `docs/specs/linalg.md`, the 4 x 3 with rows (1, 1, 2), (1, 2, 5),
    /// (1, 3, 1) and (1, 4, 9), row after row. It has three columns, so
    /// the half of its `r` below the diagonal is three entries and not the
    /// one the 4 x 2 above has.
    const THE_DESIGN_OF_THREE_COLUMNS: [f64; 12] = [
        1.0, 1.0, 2.0, //
        1.0, 2.0, 5.0, //
        1.0, 3.0, 1.0, //
        1.0, 4.0, 9.0,
    ];

    /// The `r` of that design with the diagonal made positive, rows
    /// (2, 5, 8.5), (0, 2.23606797749979, 3.801315561749642) and
    /// (0, 0, 4.9295030175464944).
    const THE_R_OF_THE_DESIGN_OF_THREE_COLUMNS: [f64; 9] = [
        2.0,
        5.0,
        8.5, //
        0.0,
        2.236_067_977_499_79,
        3.801_315_561_749_642, //
        0.0,
        0.0,
        4.929_503_017_546_494_4,
    ];

    /// The `q` of that design with the same sign taken, row after row: its
    /// third column is (0.06085806194501853, 0.3245763303734317,
    /// -0.8317268465819189, 0.44629245426346875) and its first two are
    /// those of the 4 x 2 above.
    const THE_Q_OF_THE_DESIGN_OF_THREE_COLUMNS: [f64; 12] = [
        0.5,
        -0.670_820_393_249_936_8,
        0.060_858_061_945_018_53, //
        0.5,
        -0.223_606_797_749_978_94,
        0.324_576_330_373_431_7, //
        0.5,
        0.223_606_797_749_979,
        -0.831_726_846_581_918_9, //
        0.5,
        0.670_820_393_249_936_9,
        0.446_292_454_263_468_75,
    ];

    /// The 2 x 2 design of the same place, rows (1, 1) and (1, 2), which
    /// has as many columns as rows: the smallest matrix `thin_qr` takes
    /// and the one a check of `rows` below `cols` that read `rows` at most
    /// `cols` would refuse.
    const THE_DESIGN_OF_2_BY_2: [f64; 4] = [
        1.0, 1.0, //
        1.0, 2.0,
    ];

    /// The `r` of that design with the diagonal made positive, rows
    /// (1.4142135623730951, 2.1213203435596424) and
    /// (0, 0.7071067811865475). numpy 2.5.3 gives the first row with its
    /// sign the other way round, which is the sign the spec writes it
    /// with.
    #[expect(
        clippy::approx_constant,
        reason = "these are the entries numpy 2.5.3 gives for the r of that design, which the spec writes out, and two of them land on the f64 of the square root of 2 and of its reciprocal; a test asserts the number the spec has and not a constant that happens to equal it"
    )]
    const THE_R_OF_THE_DESIGN_OF_2_BY_2: [f64; 4] = [
        1.414_213_562_373_095_1,
        2.121_320_343_559_642_4, //
        0.0,
        0.707_106_781_186_547_5,
    ];

    #[test]
    fn the_thin_qr_of_a_design_of_three_columns_gives_the_r_and_the_q_of_the_spec() {
        let factorization = thin_qr(&THE_DESIGN_OF_THREE_COLUMNS, 4, 3).unwrap();
        let factorization = with_the_diagonal_of_r_positive(&factorization, 3);
        assert!(
            !differ(
                &factorization.r,
                &THE_R_OF_THE_DESIGN_OF_THREE_COLUMNS,
                THE_TOLERANCE_OF_THE_THIN_QR
            ),
            "the r is {r:?}",
            r = factorization.r
        );
        assert!(
            !differ(
                &factorization.q,
                &THE_Q_OF_THE_DESIGN_OF_THREE_COLUMNS,
                THE_TOLERANCE_OF_THE_THIN_QR
            ),
            "the q is {q:?}",
            q = factorization.q
        );
    }

    #[test]
    fn the_thin_qr_of_a_design_of_three_columns_writes_every_entry_below_the_diagonal_as_zero() {
        // Three entries and not the one a design of two columns has, so a
        // backend that zeroed the first of them and left the rest fails
        // here: what `dorgqr` leaves below the diagonal is the vectors of
        // the `q` and what faer leaves is its own, and the doc comment of
        // `ThinQr` says that half is 0.
        let factorization = thin_qr(&THE_DESIGN_OF_THREE_COLUMNS, 4, 3).unwrap();
        let the_lower_half: Vec<f64> = factorization
            .r
            .as_chunks::<3>()
            .0
            .iter()
            .enumerate()
            .flat_map(|(row, entries)| entries.iter().take(row).copied())
            .collect();
        assert_eq!(
            the_lower_half,
            vec![0.0, 0.0, 0.0],
            "the r is {r:?}",
            r = factorization.r
        );
    }

    #[test]
    fn the_thin_qr_takes_a_design_of_as_many_columns_as_rows() {
        // `rows` at least `cols` is what `thin_qr` asks for, so the 2 x 2
        // is taken and not refused: a check that read `rows` above `cols`
        // would give `Dimension` here.
        let factorization = thin_qr(&THE_DESIGN_OF_2_BY_2, 2, 2).unwrap();
        let factorization = with_the_diagonal_of_r_positive(&factorization, 2);
        assert!(
            !differ(
                &factorization.r,
                &THE_R_OF_THE_DESIGN_OF_2_BY_2,
                THE_TOLERANCE_OF_THE_THIN_QR
            ),
            "the r is {r:?}",
            r = factorization.r
        );
    }

    #[test]
    fn the_rank_takes_the_larger_of_the_two_dimensions_into_its_tolerance() {
        // The two 2 x 2 matrices above pin the threshold and not which
        // dimension it comes from, the larger and the smaller being one
        // number there. These are 4 x 2, where the tolerance is
        // 8.881784197001252e-16 and the one the smaller dimension would
        // give is 4.440892098500626e-16: the 6e-16 lies between them, so a
        // rank that took the smaller gives 2 where numpy 2.5.3 gives 1.
        for (entry, wanted) in [(6e-16, 1_usize), (1e-15, 2)] {
            let a = [
                1.0, 0.0, //
                0.0, entry, //
                0.0, 0.0, //
                0.0, 0.0,
            ];
            assert_eq!(
                rank(&a, 4, 2).unwrap(),
                wanted,
                "the rank of the 4 x 2 whose second singular value is {entry}"
            );
        }
    }

    #[test]
    fn the_rank_of_a_matrix_of_more_columns_than_rows_is_3() {
        // Every other matrix here has at least as many rows as columns,
        // which is what `thin_qr` asks for and what `rank` does not. The
        // twelve values of the 4 x 3 above read as a 3 x 4 are rows
        // (1, 1, 2, 1), (2, 3, 1, 3) and (4, 1, 4, 5), and both backends
        // and numpy 2.5.3 give 3 for it.
        assert_eq!(
            rank(&THE_DESIGN_OF_4_BY_3_OF_A_COLUMN_THAT_REPEATS, 3, 4).unwrap(),
            3
        );
    }

    #[test]
    fn the_rank_of_a_matrix_whose_largest_singular_value_is_near_the_largest_f64() {
        // The tolerance is the largest singular value times the product of
        // the larger dimension with the distance from 1 to the next `f64`,
        // and in that order. Multiplying the largest singular value by the
        // dimension first overflows here to an infinity, above which no
        // value lies, and the rank comes back 0. Both backends and numpy
        // 2.5.3 give 1 on 23 September 2026.
        let a = [
            1e308, 0.0, //
            0.0, 1.0,
        ];
        assert_eq!(rank(&a, 2, 2).unwrap(), 1);
    }
}
