//! The linear algebra of popnei.
//!
//! A calculation that reads a block of variants as a matrix, the principal
//! component analysis, the kinship, the genome wide association study,
//! needs a few operations of linear algebra, and this crate is the one
//! place that has them. It holds nine: the product of a matrix with
//! itself, [`add_self_product_lower`]; the eigendecomposition of a
//! symmetric matrix, [`eigh_lower`]; the Cholesky factorization of a
//! symmetric positive definite one, [`cholesky_lower`], the solve of a
//! system with the matrix it factored, [`solve_with_cholesky`], and the
//! log of that matrix's determinant, [`log_determinant_with_cholesky`];
//! and the product of two matrices, [`product`], which is the other four,
//! because [`TheFirstOperand`] and [`TheSecondOperand`] each say how one
//! matrix's buffer is laid out and the two together choose among `a b`,
//! `a b'`, `a' b` and `a' b'`.
//! `docs/specs/linalg.md` says what each one gives.
//!
//! Every matrix crosses this interface as a `&[f64]` held row after row,
//! the layout of the blocks and of everything the core crate holds, with
//! its numbers of rows and of columns beside it.
//!
//! Two backends run under that interface and give the same numbers within
//! the tolerance of "How it is verified" of the spec. One is the BLAS and
//! LAPACK of the system, the routines numpy calls, which is what runs when
//! the target is not WebAssembly and the cargo feature `blas` is on, as it
//! is by default. The other is faer, a linear algebra library written in
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
    /// The eigendecomposition of an n x n matrix needs 1 + 6n + 2n² floats
    /// besides the matrix, which is 16 MB at n = 1000 and 1.6 GB at
    /// n = 10000, so the crate asks for that memory instead of taking it,
    /// and a machine that has not got it gets this error where it would
    /// otherwise see the process end. No test of popnei reaches this case.
    #[error("this machine has not the memory for {what}, {values} values")]
    Memory {
        /// What could not be allocated.
        what: &'static str,
        /// How many values it holds.
        values: usize,
    },

    /// A matrix that could not be factored at the row the value names,
    /// counting from 0: the Cholesky reached a diagonal entry that is not
    /// above 0 there, or the solve against an upper triangular matrix
    /// reached one that is 0.
    ///
    /// It is not a defect of the caller: the call was right and the matrix
    /// was what the data made it, so the module that called decides what
    /// it means there. Both backends give the same row, counted from 0:
    /// `dpotrf` gives the order of the leading corner, counting from 1,
    /// and faer an index from 0.
    #[error(
        "the matrix {argument} is singular: the factorization stopped at its row {at}, counting from 0"
    )]
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
/// the one the caller already has: a buffer of `n` rows and `c` columns
/// held row after row is the same buffer as one of `c` rows and `n`
/// columns held row after row, so a caller whose right hand sides are the
/// columns of a matrix of `n` rows passes the buffer it holds, as `c`
/// right hand sides of `n` numbers each, and nothing is copied and nothing
/// is moved. The solutions come back one row for each column of the matrix
/// the next product reads, which is [`product`] with that operand
/// [`TheSecondOperand::ByTheColumnsOfTheResult`]. `sides` is 1 at least.
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
/// [`Error::NoConvergence`] when the routine refused an argument it was
/// given, which is a defect of popnei.
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
pub fn log_determinant_with_cholesky(l: &[f64], n: usize) -> Result<f64> {
    if n == 0 {
        return Err(Error::Dimension {
            argument: "n",
            expected: "1 at least, since l is the n x n factorization to read the diagonal of"
                .to_owned(),
        });
    }
    let l = the_matrix_of(l, n, n, "l")?;
    // The entry `i`, `i` of a matrix held row after row is the value `i`
    // times `n` plus `i` of the buffer, so the diagonal is one value in
    // every `n` plus 1 from the first, and the last of them is the last
    // value of the buffer. The addition cannot overflow: `the_matrix_of`
    // has refused every `n` whose square is above 2147483647, which leaves
    // `n` at 46340 at most.
    let the_diagonal = || l.iter().copied().step_by(n.saturating_add(1));
    if !the_diagonal().all(f64::is_finite) {
        return Err(Error::NotFinite { argument: "l" });
    }
    if let Some((row, _)) = the_diagonal().enumerate().find(|(_, entry)| *entry <= 0.0) {
        return Err(Error::Singular {
            argument: "l",
            at: row,
        });
    }
    // The logs are added in the order of the rows, which is the order the
    // spec's number was taken in and the one a total of floats has to be
    // added in to be the same number on every run.
    let total: f64 = the_diagonal().map(f64::ln).sum();
    Ok(2.0 * total)
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

#[cfg(test)]
mod tests {
    use super::{
        Eigen, Error, TheFirstOperand, TheSecondOperand, add_self_product_lower, cholesky_lower,
        eigh_lower, log_determinant_with_cholesky, product, reverse_the_rows, solve_with_cholesky,
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
            "the matrix a is singular: the factorization stopped at its row 1, counting from 0"
        );
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
        // sum of their logs lands on the f64 numpy's slogdet gives for the
        // matrix, 3.58351893845611, which is the log of its determinant of
        // 36: the tolerance of the comparison is 0.
        let mut l = the_matrix_to_factor();
        cholesky_lower(&mut l, 3).unwrap();
        let logarithm = log_determinant_with_cholesky(&l, 3).unwrap();
        assert!(
            !differ(&[logarithm], &[3.58351893845611], 0.0),
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
            !differ(&[logarithm], &[3.58351893845611], 0.0),
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
                "the matrix l is singular: the factorization stopped at its row 1, counting from 0"
            );
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
}
