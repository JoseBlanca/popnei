//! The BLAS and LAPACK backend: the routines of the library of the
//! system, `dsyrk`, `dgemm`, which the four products of this module call,
//! `dsyevd`, `dpotrf`, `dpotrs`, `dpotri`, `dgeqrf`, `dorgqr`, `dtrtrs`
//! and `dgesdd`, the ones numpy calls.
//!
//! Every matrix reaches this module row after row, and these routines read
//! a matrix column after column. The buffer of an r x c matrix read that
//! way is its transpose, c x r, so each routine is called on the
//! transposes: the lower half of `g` in popnei's layout is the upper half
//! for the routine, so `uplo` is `U`; `a'a` is `a a'` of the transposed
//! view, so `trans` is `N`; and `c = a b` is `c' = b' a'`, so `dgemm` gets
//! the buffer of `b` as its first operand and that of `a` as its second.
//! `c = a b'` is `c' = b a'` the same way, and there the first operand of
//! the routine is the transpose of what the buffer of `b` is in its view,
//! so that call has `transa` `T` and `transb` `N`. `c = a' b` is
//! `c' = b' a`, which turns the routine's second operand instead, so it
//! has `transa` `N` and `transb` `T`; and `c = a' b'` is `c' = b a`,
//! whose two flags are both `T`. None of the four copies a buffer.
//!
//! The thin QR and the singular values are the two functions here that do
//! copy: `dgeqrf`, `dorgqr` and `dgesdd` are much slower on the wide
//! matrix that the buffer of a design is in their view than on the tall
//! one it is, so each writes the transpose of that buffer into one of its
//! own and calls them on that.
//!
//! The functions here are given slices whose lengths the caller has
//! already cut to the dimensions, and they check nothing else: the checks
//! are in `lib.rs`, where they hold for every backend. Each `unsafe` block
//! says why the slices are long enough for the dimensions it passes.

use crate::{Eigen, Error, Result};

// The two crates below hold no code of their own: each emits the argument
// that links the library holding the routines, `-framework Accelerate` on
// this Mac. Nothing of this crate names them, so without these two lines
// the linker would not see them and every routine would be undefined.
extern crate blas_src as _;
extern crate lapack_src as _;

/// Adds the lower half of `a'a` to `g`, with `a` of exactly `rows` x
/// `cols` values and `g` of exactly `cols` x `cols`, both row after row
/// and `rows` and `cols` 1 at least.
///
/// # Errors
///
/// [`Error::Dimension`] when a dimension is larger than the `i32` the
/// routine takes.
pub(crate) fn add_self_product_lower(
    a: &[f64],
    rows: usize,
    cols: usize,
    g: &mut [f64],
) -> Result<()> {
    let order = the_i32_of(cols, "cols")?;
    let inner = the_i32_of(rows, "rows")?;
    // SAFETY: with `uplo` U, `trans` N, `n` = cols, `k` = rows and `lda` =
    // cols, dsyrk reads `a` as a column major matrix of cols rows and rows
    // columns, which is the cols * rows values of the transpose of `a`,
    // and `a` holds exactly that many; and with `beta` 1 and `ldc` = cols
    // it reads and writes the upper triangle of `g` as a column major
    // matrix of cols x cols, which is cols * cols values, and `g` holds
    // exactly that many. Neither dimension is 0 and both fit in the `i32`
    // the routine takes, which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of BLAS are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::blas::dsyrk(b'U', b'N', order, inner, 1.0, a, order, 1.0, g, order);
    }
    Ok(())
}

/// Writes `a b` into `c`, with `a` of exactly `rows` x `inner` values, `b`
/// of `inner` x `cols` and `c` of `rows` x `cols`, all row after row and
/// every dimension 1 at least.
///
/// # Errors
///
/// [`Error::Dimension`] when a dimension is larger than the `i32` the
/// routine takes.
pub(crate) fn product(
    a: &[f64],
    rows: usize,
    inner: usize,
    b: &[f64],
    cols: usize,
    c: &mut [f64],
) -> Result<()> {
    let m = the_i32_of(cols, "cols")?;
    let n = the_i32_of(rows, "rows")?;
    let k = the_i32_of(inner, "inner")?;
    // SAFETY: `c = a b` in popnei's layout is `c' = b' a'` in the layout
    // the routine reads, so with both `trans` N, `m` = cols, `n` = rows
    // and `k` = inner the routine reads `b` as a column major matrix of
    // cols rows and inner columns with `lda` = cols, which is the inner *
    // cols values of `b`; `a` as one of inner rows and rows columns with
    // `ldb` = inner, which is the rows * inner values of `a`; and, with
    // `beta` 0, writes a column major matrix of cols rows and rows columns
    // with `ldc` = cols, which is the rows * cols values of `c`. Each of
    // the three slices holds exactly the values of its dimensions, no
    // dimension is 0 and all three fit in the `i32` the routine takes,
    // which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of BLAS are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::blas::dgemm(b'N', b'N', m, n, k, 1.0, b, m, a, k, 0.0, c, m);
    }
    Ok(())
}

/// Writes `a b'` into `c`, with `a` of exactly `rows` x `inner` values,
/// `b` of `cols` x `inner` and `c` of `rows` x `cols`, all row after row
/// and every dimension 1 at least.
///
/// # Errors
///
/// [`Error::Dimension`] when a dimension is larger than the `i32` the
/// routine takes.
pub(crate) fn product_with_the_second_turned(
    a: &[f64],
    rows: usize,
    inner: usize,
    b: &[f64],
    cols: usize,
    c: &mut [f64],
) -> Result<()> {
    let m = the_i32_of(cols, "cols")?;
    let n = the_i32_of(rows, "rows")?;
    let k = the_i32_of(inner, "inner")?;
    // SAFETY: `c = a b'` in popnei's layout is `c' = b a'` in the layout
    // the routine reads, where `b` is the transpose of the column major
    // matrix the buffer of `b` is, so `transa` is T and `transb` is N.
    // With `m` = cols, `n` = rows and `k` = inner the routine reads the
    // buffer of `b` as a column major matrix of inner rows and cols
    // columns with `lda` = inner, which `transa` T turns into the cols x
    // inner first operand and which is the cols * inner values of `b`;
    // the buffer of `a` as one of inner rows and rows columns with `ldb`
    // = inner, which is the rows * inner values of `a`; and, with `beta`
    // 0, writes a column major matrix of cols rows and rows columns with
    // `ldc` = cols, which is the rows * cols values of `c`. Each of the
    // three slices holds exactly the values of its dimensions, no
    // dimension is 0 and all three fit in the `i32` the routine takes,
    // which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of BLAS are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::blas::dgemm(b'T', b'N', m, n, k, 1.0, b, k, a, k, 0.0, c, m);
    }
    Ok(())
}

/// Writes `a' b` into `c`, with `a` of exactly `inner` x `rows` values,
/// `b` of `inner` x `cols` and `c` of `rows` x `cols`, all row after row
/// and every dimension 1 at least.
///
/// # Errors
///
/// [`Error::Dimension`] when a dimension is larger than the `i32` the
/// routine takes.
pub(crate) fn product_with_the_first_turned(
    a: &[f64],
    rows: usize,
    inner: usize,
    b: &[f64],
    cols: usize,
    c: &mut [f64],
) -> Result<()> {
    let m = the_i32_of(cols, "cols")?;
    let n = the_i32_of(rows, "rows")?;
    let k = the_i32_of(inner, "inner")?;
    // SAFETY: `c = a' b` in popnei's layout is `c' = b' a` in the layout
    // the routine reads, where `b'` is the column major matrix the buffer
    // of `b` is and `a` is the transpose of the one the buffer of `a` is,
    // so `transa` is N and `transb` is T. With `m` = cols, `n` = rows and
    // `k` = inner the routine reads the buffer of `b` as a column major
    // matrix of cols rows and inner columns with `lda` = cols, which is
    // the inner * cols values of `b`; the buffer of `a` as one of rows
    // rows and inner columns with `ldb` = rows, which `transb` T turns
    // into the inner x rows second operand and which is the inner * rows
    // values of `a`; and, with `beta` 0, writes a column major matrix of
    // cols rows and rows columns with `ldc` = cols, which is the rows *
    // cols values of `c`. Each of the three slices holds exactly the
    // values of its dimensions, no dimension is 0 and all three fit in
    // the `i32` the routine takes, which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of BLAS are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::blas::dgemm(b'N', b'T', m, n, k, 1.0, b, m, a, n, 0.0, c, m);
    }
    Ok(())
}

/// Writes `a' b'` into `c`, with `a` of exactly `inner` x `rows` values,
/// `b` of `cols` x `inner` and `c` of `rows` x `cols`, all row after row
/// and every dimension 1 at least.
///
/// # Errors
///
/// [`Error::Dimension`] when a dimension is larger than the `i32` the
/// routine takes.
pub(crate) fn product_with_both_turned(
    a: &[f64],
    rows: usize,
    inner: usize,
    b: &[f64],
    cols: usize,
    c: &mut [f64],
) -> Result<()> {
    let m = the_i32_of(cols, "cols")?;
    let n = the_i32_of(rows, "rows")?;
    let k = the_i32_of(inner, "inner")?;
    // SAFETY: `c = a' b'` in popnei's layout is `c' = b a` in the layout
    // the routine reads, and there each of the two operands is the
    // transpose of the column major matrix its buffer is, so both `trans`
    // flags are T. With `m` = cols, `n` = rows and `k` = inner the routine
    // reads the buffer of `b` as a column major matrix of inner rows and
    // cols columns with `lda` = inner, which `transa` T turns into the
    // cols x inner first operand and which is the cols * inner values of
    // `b`; the buffer of `a` as one of rows rows and inner columns with
    // `ldb` = rows, which `transb` T turns into the inner x rows second
    // operand and which is the inner * rows values of `a`; and, with
    // `beta` 0, writes a column major matrix of cols rows and rows columns
    // with `ldc` = cols, which is the rows * cols values of `c`. Each of
    // the three slices holds exactly the values of its dimensions, no
    // dimension is 0 and all three fit in the `i32` the routine takes,
    // which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of BLAS are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::blas::dgemm(b'T', b'T', m, n, k, 1.0, b, k, a, n, 0.0, c, m);
    }
    Ok(())
}

/// The eigendecomposition of the symmetric `g`, of exactly `n` x `n`
/// values row after row with its lower half filled and `n` 1 at least.
///
/// The eigenvalues come back from the smallest, which is the order the
/// routine gives them in, and each eigenvector is a row of the buffer:
/// `lib.rs` turns both round together.
///
/// # Errors
///
/// [`Error::Dimension`] when the workspace the routine asks for is larger
/// than the `i32` the routine takes its length as.
/// [`Error::NoConvergence`] when the routine gave an `info` other than 0.
/// [`Error::Memory`] when a workspace could not be allocated.
pub(crate) fn eigh_lower(mut g: Vec<f64>, n: usize) -> Result<Eigen> {
    let order = the_i32_of(n, "n")?;
    let mut values = vec![0.0_f64; n];
    let mut info = 0_i32;

    // The routine says how much it wants to work in when it is called with
    // the two lengths at -1, which is how LAPACK is asked, and it writes
    // the two numbers into the first entry of each workspace.
    let mut floats_asked = [0.0_f64; 1];
    let mut integers_asked = [0_i32; 1];
    // SAFETY: with the two lengths at -1 the routine writes the first
    // entry of `work` and of `iwork` and reads nothing else of them, and
    // each holds one value; it reads and writes nothing of `a`, of `w` or
    // of `info` other than to store the sizes, and `a` holds n * n values,
    // `w` holds n and `info` is one integer. `n` fits in the `i32` the
    // routine takes, which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dsyevd(
            b'V',
            b'U',
            order,
            &mut g,
            order,
            &mut values,
            &mut floats_asked,
            -1,
            &mut integers_asked,
            -1,
            &mut info,
        );
    }
    if info != 0 {
        return Err(Error::NoConvergence {
            routine: "dsyevd",
            info,
        });
    }

    let (floats, integers) = the_workspace_of(
        n,
        floats_asked.first().copied().unwrap_or(0.0),
        integers_asked.first().copied().unwrap_or(0),
    )?;
    let lwork = the_length_of_a_workspace(floats)?;
    let liwork = the_length_of_a_workspace(integers)?;
    // The workspace of an n of 10000 is 1.6 GB, so it is asked for and not
    // taken: `vec!` on a machine that has not the memory ends the process,
    // and the core crate gives an error instead, as its reader of blocks
    // does for the columns it allocates.
    let mut work: Vec<f64> = Vec::new();
    work.try_reserve_exact(floats).map_err(|_| Error::Memory {
        what: "the workspace of floats of dsyevd",
        values: floats,
    })?;
    work.resize(floats, 0.0);
    let mut iwork: Vec<i32> = Vec::new();
    iwork
        .try_reserve_exact(integers)
        .map_err(|_| Error::Memory {
            what: "the workspace of integers of dsyevd",
            values: integers,
        })?;
    iwork.resize(integers, 0);
    // SAFETY: with `jobz` V and `uplo` U the routine reads the upper
    // triangle of `a` as a column major matrix of n x n with `lda` = n,
    // which is the lower half of `g` in popnei's layout, and overwrites
    // the whole of it with the eigenvectors, n * n values, which is what
    // `g` holds; it writes the n eigenvalues into `w`, which holds n; and
    // it works in the first `lwork` values of `work` and the first
    // `liwork` of `iwork`, which hold exactly that many, `lwork` and
    // `liwork` being the larger of what the query above asked for and the
    // minimum the routine documents. All four lengths fit in the `i32` the
    // routine takes, which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dsyevd(
            b'V',
            b'U',
            order,
            &mut g,
            order,
            &mut values,
            &mut work,
            lwork,
            &mut iwork,
            liwork,
            &mut info,
        );
    }
    if info != 0 {
        return Err(Error::NoConvergence {
            routine: "dsyevd",
            info,
        });
    }

    // The routine gives the eigenvalues from the smallest and each
    // eigenvector as a column of the matrix it read column after column,
    // which in popnei's buffer is a row, so the buffer is already the
    // eigenvectors a row each. `lib.rs` turns the values and the rows
    // round together.
    Ok(Eigen { values, vectors: g })
}

/// The Cholesky factorization of the symmetric positive definite `a`, of
/// exactly `n` x `n` values row after row with its lower half filled and
/// `n` 1 at least, which overwrites that lower half with the lower
/// triangular `l`.
///
/// # Errors
///
/// [`Error::Dimension`] when `n` is larger than the `i32` the routine
/// takes. [`Error::Singular`] when `a` is not positive definite, with the
/// row the routine stopped at counted from 0. [`Error::NoConvergence`]
/// when the routine refused an argument, which is a defect of popnei.
pub(crate) fn cholesky_lower(a: &mut [f64], n: usize) -> Result<()> {
    let order = the_i32_of(n, "n")?;
    let mut info = 0_i32;
    // SAFETY: with `uplo` U, `n` = n and `lda` = n the routine reads the
    // upper triangle of `a` as a column major matrix of n x n, which is
    // the lower half of `a` in popnei's layout, and overwrites that same
    // triangle with the factorization; both are inside the n * n values of
    // `a`, which is what it holds. It reads and writes nothing else of
    // `a`, and `info` is one integer. `n` is not 0 and fits in the `i32`
    // the routine takes, which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dpotrf(b'U', order, a, order, &mut info);
    }
    if info == 0 {
        return Ok(());
    }
    // An `info` above 0 is the order of the leading corner that is not
    // positive definite, counting from 1, and the row the crate gives is
    // that corner's last, counting from 0, which is one less. An `info`
    // below 0 is an argument the routine refused, and that is the arm the
    // conversion fails in, since no negative number is a count.
    match usize::try_from(info) {
        Ok(corner) => Err(Error::Singular {
            argument: "a",
            at: corner.saturating_sub(1),
        }),
        Err(_) => Err(Error::NoConvergence {
            routine: "dpotrf",
            info,
        }),
    }
}

/// The `x` of `a x = b` for the factorization `l` of exactly `n` x `n`
/// values row after row with its lower half filled, and `b` of exactly
/// `sides` x `n` values row after row, one row for each right hand side,
/// which comes back holding the solutions the same way. `n` and `sides`
/// are 1 at least.
///
/// The buffer of `b` read column after column is the n x `sides` matrix
/// whose columns are the right hand sides, which is what the routine
/// takes, so the layout of this crate is the layout `dpotrs` wants and
/// nothing is copied here either.
///
/// # Errors
///
/// [`Error::Dimension`] when a dimension is larger than the `i32` the
/// routine takes. [`Error::NoConvergence`] when the routine refused an
/// argument it was given, which is a defect of popnei.
pub(crate) fn solve_with_cholesky(l: &[f64], n: usize, b: &mut [f64], sides: usize) -> Result<()> {
    let order = the_i32_of(n, "n")?;
    let right_hand_sides = the_i32_of(sides, "sides")?;
    let mut info = 0_i32;
    // SAFETY: with `uplo` U, `n` = n and `lda` = n the routine reads the
    // upper triangle of `l` as a column major matrix of n x n, which is
    // the lower half of `l` in popnei's layout and is inside the n * n
    // values `l` holds; and with `nrhs` = sides and `ldb` = n it reads and
    // writes `b` as a column major matrix of n rows and sides columns,
    // which is the sides * n values `b` holds. It writes nothing else, and
    // `info` is one integer. Neither dimension is 0 and both fit in the
    // `i32` the routine takes, which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dpotrs(b'U', order, right_hand_sides, l, order, b, order, &mut info);
    }
    if info == 0 {
        Ok(())
    } else {
        // `dpotrs` gives an `info` other than 0 for one reason, an
        // argument it refused, which is a defect of popnei and never a
        // matrix it could not solve with.
        Err(Error::NoConvergence {
            routine: "dpotrs",
            info,
        })
    }
}

/// The lower half of the inverse of the `a` whose factorization `l` is,
/// with `l` of exactly `n` x `n` values row after row with its lower half
/// filled and `inverse` of exactly `n` x `n`, and `n` 1 at least. The
/// upper half of `inverse` is left as it was.
///
/// `dpotri` inverts a factorization where it lies, so the lower half of
/// the factorization is copied into the buffer the caller gave for the
/// inverse and the routine works there: that copy is what keeps `l` as it
/// was, which the interface of the crate promises, and the routine asks
/// for no workspace besides it. The copy is of that half alone, and not of
/// the whole buffer, because the upper half of `inverse` is the caller's
/// and is left as it was, and because the routine reads nothing else.
///
/// # Errors
///
/// [`Error::Dimension`] when `n` is larger than the `i32` the routine
/// takes. [`Error::Singular`] when the diagonal of `l` holds a 0 at the
/// row the error names, which the routine would divide by; `lib.rs` reads
/// that diagonal before either backend runs, since faer's inverse does not
/// look at it, so no caller of the crate reaches this one.
/// [`Error::NoConvergence`] when the routine refused an argument it was
/// given, which is a defect of popnei.
pub(crate) fn invert_with_cholesky(l: &[f64], n: usize, inverse: &mut [f64]) -> Result<()> {
    let order = the_i32_of(n, "n")?;
    // Both buffers hold exactly n * n values, `lib.rs` having cut them to
    // the dimensions, so each is n rows of n, and the lower half is the
    // entries of column `j` at most `i` of row `i`. What is left of the
    // row after those is the caller's and is not written.
    for (row, (into, from)) in inverse
        .chunks_exact_mut(n)
        .zip(l.chunks_exact(n))
        .enumerate()
    {
        for (into, from) in into.iter_mut().zip(from).take(row.saturating_add(1)) {
            *into = *from;
        }
    }
    let mut info = 0_i32;
    // SAFETY: with `uplo` U, `n` = n and `lda` = n the routine reads and
    // writes the upper triangle of `inverse` as a column major matrix of
    // n x n, which is the lower half of `inverse` in popnei's layout and
    // is inside the n * n values it holds, the copy of `l` just written
    // there. It reads and writes nothing else of that buffer and nothing
    // at all of `l`, and `info` is one integer. `n` is not 0 and fits in
    // the `i32` the routine takes, which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dpotri(b'U', order, inverse, order, &mut info);
    }
    if info == 0 {
        return Ok(());
    }
    // An `info` above 0 is the row of the diagonal entry of the
    // factorization that is 0, counting from 1, which no `l` that
    // `cholesky_lower` gave holds, since it stops at the first entry that
    // is not above 0. An `info` below 0 is an argument the routine
    // refused, and that is the arm the conversion fails in, since no
    // negative number is a count.
    match usize::try_from(info) {
        Ok(row) => Err(Error::Singular {
            argument: "l",
            at: row.saturating_sub(1),
        }),
        Err(_) => Err(Error::NoConvergence {
            routine: "dpotri",
            info,
        }),
    }
}

/// The thin QR of `a`, of exactly `rows` x `cols` values row after row
/// with `rows` at least `cols` and `cols` 1 at least: `q` of exactly
/// `rows` x `cols` values and the upper triangular `r` of exactly `cols`
/// x `cols`, both written row after row and the lower half of `r` set
/// to 0.
///
/// The buffer of `a` read column after column is the `cols` x `rows`
/// matrix, the wide one, and the two routines are much slower on it than
/// on the tall one: a design of 10000 x 5 took 2.98 ms that way and 0.165
/// ms through the copy below, measured on 23 September 2026 and written
/// in "What the seven of the GWAS cost" of `docs/specs/linalg.md`. So
/// this backend writes the transpose of `a` into a buffer of its own,
/// `rows` x `cols` values and 400 KB at that size, which is `a` held
/// column after column, and calls them on that.
///
/// `dgeqrf` leaves the factorization in that buffer: its upper triangle
/// is the `r`, and below the diagonal are the vectors that `dorgqr`
/// builds the `q` from, which it writes over the whole of the buffer. So
/// `r` is read out between the two calls.
///
/// # Errors
///
/// [`Error::Dimension`] when a dimension, or the workspace a routine
/// asks for, is larger than the `i32` the routines take.
/// [`Error::NoConvergence`] when a routine refused an argument it was
/// given, which is a defect of popnei: neither has another reason to give
/// an `info` other than 0.
pub(crate) fn thin_qr(
    a: &[f64],
    rows: usize,
    cols: usize,
    q: &mut [f64],
    r: &mut [f64],
) -> Result<()> {
    let m = the_i32_of(rows, "rows")?;
    let n = the_i32_of(cols, "cols")?;
    let mut column_major = the_column_major_copy_of(a, rows, cols);
    // One coefficient for each column, which is what the factorization
    // keeps beside the vectors it leaves in the matrix.
    let mut coefficients = vec![0.0_f64; cols];
    let mut info = 0_i32;

    // Each routine says how much it wants to work in when it is called
    // with the length of its workspace at -1, which is how LAPACK is
    // asked, and writes that number into the first entry of the workspace
    // it was given. It is what `eigh_lower` above asks `dsyevd`.
    let mut asked = [0.0_f64; 1];
    // SAFETY: with `lwork` at -1 the routine writes the first entry of
    // `work` and reads nothing else of it, and `asked` holds one value;
    // it reads and writes nothing of `a`, of `tau` or of `info` other
    // than to store that size, and `column_major` holds rows * cols
    // values, `coefficients` holds cols and `info` is one integer.
    // Neither dimension is 0 and both fit in the `i32` the routine takes,
    // which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dgeqrf(
            m,
            n,
            &mut column_major,
            m,
            &mut coefficients,
            &mut asked,
            -1,
            &mut info,
        );
    }
    if info != 0 {
        return Err(Error::NoConvergence {
            routine: "dgeqrf",
            info,
        });
    }
    let mut work =
        vec![0.0_f64; the_workspace_of_the_thin_qr(cols, asked.first().copied().unwrap_or(0.0))];
    let lwork = the_length_of_the_workspace_of_the_thin_qr(work.len())?;
    // SAFETY: with `m` = rows, `n` = cols and `lda` = rows the routine
    // reads and overwrites `a` as a column major matrix of rows x cols,
    // which is the rows * cols values `column_major` holds; it writes one
    // coefficient for each of the min(rows, cols) columns into `tau`,
    // which is cols of them since rows is at least cols, and
    // `coefficients` holds cols; and it works in the first `lwork` values
    // of `work`, which holds exactly that many. It writes nothing else,
    // and `info` is one integer. Neither dimension is 0 and all three
    // lengths fit in the `i32` the routine takes, which `the_i32_of` and
    // `the_length_of_the_workspace_of_the_thin_qr` have just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dgeqrf(
            m,
            n,
            &mut column_major,
            m,
            &mut coefficients,
            &mut work,
            lwork,
            &mut info,
        );
    }
    if info != 0 {
        return Err(Error::NoConvergence {
            routine: "dgeqrf",
            info,
        });
    }

    // The leading cols x cols corner of the factorization, before
    // `dorgqr` writes the `q` over it. What is below its diagonal is the
    // vectors of that `q` and no part of `r`, so that half is set to 0,
    // which is what the interface of the crate gives.
    write_the_rows_of(&column_major, rows, cols, r);
    for (row, entries) in r.chunks_exact_mut(cols).enumerate() {
        for entry in entries.iter_mut().take(row) {
            *entry = 0.0;
        }
    }

    let mut asked = [0.0_f64; 1];
    // SAFETY: with `lwork` at -1 the routine writes the first entry of
    // `work` and reads nothing else of it, and `asked` holds one value;
    // it reads and writes nothing of `a`, of `tau` or of `info` other
    // than to store that size, and the three slices hold what the call
    // below says they hold. No dimension is 0 and all three fit in the
    // `i32` the routine takes, which `the_i32_of` has already checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dorgqr(
            m,
            n,
            n,
            &mut column_major,
            m,
            &coefficients,
            &mut asked,
            -1,
            &mut info,
        );
    }
    if info != 0 {
        return Err(Error::NoConvergence {
            routine: "dorgqr",
            info,
        });
    }
    let mut work =
        vec![0.0_f64; the_workspace_of_the_thin_qr(cols, asked.first().copied().unwrap_or(0.0))];
    let lwork = the_length_of_the_workspace_of_the_thin_qr(work.len())?;
    // SAFETY: with `m` = rows, `n` = cols, `k` = cols and `lda` = rows
    // the routine reads the vectors that `dgeqrf` left in `a` as a column
    // major matrix of rows x cols and overwrites it with the first cols
    // columns of the `q`, which is the rows * cols values `column_major`
    // holds; it reads one coefficient for each of the k = cols vectors
    // from `tau`, which holds cols; and it works in the first `lwork`
    // values of `work`, which holds exactly that many. It writes nothing
    // else, and `info` is one integer. No dimension is 0, `k` is at most
    // `n` and `n` at most `m` since rows is at least cols, and all three
    // lengths fit in the `i32` the routine takes, which `the_i32_of` and
    // `the_length_of_the_workspace_of_the_thin_qr` have just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dorgqr(
            m,
            n,
            n,
            &mut column_major,
            m,
            &coefficients,
            &mut work,
            lwork,
            &mut info,
        );
    }
    if info != 0 {
        return Err(Error::NoConvergence {
            routine: "dorgqr",
            info,
        });
    }
    write_the_rows_of(&column_major, rows, cols, q);
    Ok(())
}

/// The `rows` x `cols` matrix that the buffer holds row after row, held
/// column after column instead, which is its transpose written out and
/// what the routines of the thin QR are fast on. The buffer holds exactly
/// `rows` x `cols` values and `rows` and `cols` are 1 at least.
///
/// The column `j` of a matrix held row after row is one value in every
/// `cols` from the value `j`, and those are the values of the column `j`
/// of the copy, one after another.
fn the_column_major_copy_of(a: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    let mut column_major = vec![0.0_f64; a.len()];
    for (column, into) in column_major.chunks_exact_mut(rows).enumerate() {
        for (into, from) in into.iter_mut().zip(a.iter().skip(column).step_by(cols)) {
            *into = *from;
        }
    }
    column_major
}

/// Writes into `into`, row after row and `cols` to a row, as many first
/// rows of the matrix that `column_major` holds column after column with
/// `rows` rows as `into` has room for: the whole of it for a `q` of
/// `rows` x `cols`, and the leading `cols` x `cols` corner for an `r`.
///
/// The row `i` of a matrix held column after column with `rows` rows is
/// one value in every `rows` from the value `i`, which is what each row
/// of `into` takes `cols` of.
fn write_the_rows_of(column_major: &[f64], rows: usize, cols: usize, into: &mut [f64]) {
    for (row, entries) in into.chunks_exact_mut(cols).enumerate() {
        for (entry, from) in entries
            .iter_mut()
            .zip(column_major.iter().skip(row).step_by(rows))
        {
            *entry = *from;
        }
    }
}

/// How many floats `dgeqrf` and `dorgqr` work in for a matrix of `cols`
/// columns: the larger of what the query asked for and `cols`, the
/// minimum both routines document.
///
/// The query writes its number as an `f64`. An infinity, a NaN, a
/// negative number and one above the `i32` the length is passed as are
/// left out, and the minimum stands; a value that is a count is taken by
/// its whole part, the fraction that the routine cannot have meant being
/// dropped. It is what `the_workspace_of` below does with the query of
/// `dsyevd`.
fn the_workspace_of_the_thin_qr(cols: usize, asked: f64) -> usize {
    let the_most_a_length_holds = f64::from(i32::MAX);
    let asked = if asked.is_finite() && asked >= 0.0 && asked <= the_most_a_length_holds {
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "the line above has checked that the value is finite, is not negative and is at most i32::MAX, so it is a count this machine holds; what the cast drops is the fraction, which the routine cannot have meant, and the minimum below stands when what is left is smaller than it"
        )]
        let asked = asked as usize;
        asked
    } else {
        0
    };
    asked.max(cols)
}

/// The length of a workspace of the thin QR as the `i32` the two routines
/// take it as.
///
/// [`the_workspace_of_the_thin_qr`] gives no number above that `i32`, and
/// `cols` is at most 46340 since `lib.rs` refuses a matrix of more than
/// 2147483647 values, so this refuses nothing that reaches it: it is how
/// the conversion is made without an `as` that could truncate in silence.
///
/// # Errors
///
/// [`Error::Dimension`] when the workspace is larger than that.
fn the_length_of_the_workspace_of_the_thin_qr(values: usize) -> Result<i32> {
    i32::try_from(values).map_err(|_| Error::Dimension {
        argument: "cols",
        expected: format!(
            "small enough that the workspace dgeqrf and dorgqr ask for, {values} values here, is at most the {largest} their lengths are passed as",
            largest = i32::MAX
        ),
    })
}

/// How many floats and how many integers `dsyevd` works in for a matrix of
/// `n` x `n` with its eigenvectors: the larger of what the query asked
/// for and the minimum the routine documents, 1 + 6n + 2n² floats and
/// 3 + 5n integers.
///
/// The query writes the number of floats as an `f64`. An infinity, a NaN,
/// a negative number and one above the `i32` the length is passed as are
/// left out, and the minimum stands; a value that is a count is taken by
/// its whole part, the fraction that the routine cannot have meant being
/// dropped. The same for the integers, whose query is an `i32` already.
///
/// # Errors
///
/// [`Error::Dimension`] when the minimum for `n` is more values than this
/// machine can hold.
fn the_workspace_of(n: usize, floats_asked: f64, integers_asked: i32) -> Result<(usize, usize)> {
    let too_large = || Error::Dimension {
        argument: "n",
        expected: format!(
            "small enough for the workspace of dsyevd, 1 + 6n + 2n² floats, to fit in this machine, and n is {n}"
        ),
    };
    let squared = n.checked_mul(n).ok_or_else(too_large)?;
    let floats_at_least = 2_usize
        .checked_mul(squared)
        .and_then(|values| values.checked_add(6_usize.checked_mul(n)?))
        .and_then(|values| values.checked_add(1))
        .ok_or_else(too_large)?;
    let integers_at_least = 5_usize
        .checked_mul(n)
        .and_then(|values| values.checked_add(3))
        .ok_or_else(too_large)?;

    let the_most_a_length_holds = f64::from(i32::MAX);
    let floats = if floats_asked.is_finite()
        && floats_asked >= 0.0
        && floats_asked <= the_most_a_length_holds
    {
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "the line above has checked that the value is finite, is not negative and is at most i32::MAX, so it is a count this machine holds; what the cast drops is the fraction, which the routine cannot have meant, and the minimum below stands when what is left is smaller than it"
        )]
        let asked = floats_asked as usize;
        asked
    } else {
        0
    };
    let integers = usize::try_from(integers_asked).unwrap_or(0);
    Ok((floats.max(floats_at_least), integers.max(integers_at_least)))
}

/// The dimension as the `i32` the routines of BLAS and LAPACK take.
///
/// `lib.rs` has already refused every dimension and every number of values
/// of a matrix above that, for both backends, so this fails for no call it
/// lets through.
///
/// # Errors
///
/// [`Error::Dimension`] when it is larger than that.
fn the_i32_of(value: usize, argument: &'static str) -> Result<i32> {
    i32::try_from(value).map_err(|_| Error::Dimension {
        argument,
        expected: format!(
            "at most {largest}, which is what the routines of BLAS and LAPACK take, and it is {value}",
            largest = i32::MAX
        ),
    })
}

/// The length of a workspace of `dsyevd` as the `i32` the routine takes
/// it as. It is the one number of a call that `lib.rs` cannot check,
/// because the routine says how large it is.
///
/// # Errors
///
/// [`Error::Dimension`] when it is larger than that, which for the 2n²
/// floats of the workspace happens at an `n` of about 32768.
fn the_length_of_a_workspace(values: usize) -> Result<i32> {
    i32::try_from(values).map_err(|_| Error::Dimension {
        argument: "n",
        expected: format!(
            "small enough that the workspace dsyevd asks for, {values} values here, is at most the {largest} its length is passed as",
            largest = i32::MAX
        ),
    })
}

/// The `x` of `r x = b` for the upper triangular `r` of exactly `n` x `n`
/// values row after row with its upper half filled, and `b` of exactly
/// `sides` x `n` values row after row, one row for each right hand side,
/// which comes back holding the solutions the same way. `n` and `sides`
/// are 1 at least.
///
/// The buffer of `r` read column after column is its transpose, which is
/// lower triangular, so the routine is told `uplo` L and `trans` T, which
/// asks it for the solve against the transpose of what it read: the two
/// turns undo each other and what is solved against is popnei's `r`. The
/// buffer of `b` read that way is the n x `sides` matrix whose columns are
/// the right hand sides, which is what the routine takes, as it is for
/// `dpotrs` above, so nothing is copied here either.
///
/// # Errors
///
/// [`Error::Dimension`] when a dimension is larger than the `i32` the
/// routine takes. [`Error::Singular`] when the diagonal of `r` holds a 0
/// at the row the error names, which the routine would divide by;
/// `lib.rs` reads that diagonal before either backend runs, since faer
/// does not look at it, so no caller of the crate reaches this one.
/// [`Error::NoConvergence`] when the routine refused an argument it was
/// given, which is a defect of popnei.
pub(crate) fn solve_upper_triangular(
    r: &[f64],
    n: usize,
    b: &mut [f64],
    sides: usize,
) -> Result<()> {
    let order = the_i32_of(n, "n")?;
    let right_hand_sides = the_i32_of(sides, "sides")?;
    let mut info = 0_i32;
    // SAFETY: with `uplo` L, `trans` T, `diag` N, `n` = n and `lda` = n
    // the routine reads the lower triangle of `r` as a column major
    // matrix of n x n, which is the upper half of `r` in popnei's layout
    // and is inside the n * n values `r` holds; and with `nrhs` = sides
    // and `ldb` = n it reads and writes `b` as a column major matrix of n
    // rows and sides columns, which is the sides * n values `b` holds. It
    // writes nothing else, and `info` is one integer. Neither dimension
    // is 0 and both fit in the `i32` the routine takes, which
    // `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dtrtrs(
            b'L',
            b'T',
            b'N',
            order,
            right_hand_sides,
            r,
            order,
            b,
            order,
            &mut info,
        );
    }
    if info == 0 {
        return Ok(());
    }
    // An `info` above 0 is the row of the diagonal entry that is 0,
    // counting from 1, which `lib.rs` has already refused. An `info` below
    // 0 is an argument the routine refused, and that is the arm the
    // conversion fails in, since no negative number is a count.
    match usize::try_from(info) {
        Ok(row) => Err(Error::Singular {
            argument: "r",
            at: row.saturating_sub(1),
        }),
        Err(_) => Err(Error::NoConvergence {
            routine: "dtrtrs",
            info,
        }),
    }
}

/// The singular values of `a`, of exactly `rows` x `cols` values row after
/// row and both dimensions 1 at least, from the largest: as many as the
/// smaller dimension.
///
/// The buffer of `a` read column after column is the `cols` x `rows`
/// matrix, the wide one, and the routine is much slower on it than on the
/// tall one: a design of 10000 x 5 took 1.54 ms that way and 0.145 ms
/// through the copy below, measured on 23 September 2026 and written in
/// "What the seven of the GWAS cost" of `docs/specs/linalg.md`. So this
/// backend writes the transpose of `a` into a buffer of its own, which is
/// `a` held column after column, and calls the routine on that, as the
/// thin QR above does. The routine overwrites that buffer, and `a` itself
/// is left as it was.
///
/// `jobz` is `N`, which computes the values and neither of the two
/// matrices of vectors, so the two buffers they would go in hold one
/// value each and are never written.
///
/// # Errors
///
/// [`Error::Dimension`] when a dimension, or the workspace the routine
/// asks for, is larger than the `i32` the routine takes.
/// [`Error::NoConvergence`] when the routine gave an `info` other than 0,
/// which is either a decomposition that did not come out or an argument
/// it refused.
pub(crate) fn singular_values(a: &[f64], rows: usize, cols: usize) -> Result<Vec<f64>> {
    let m = the_i32_of(rows, "rows")?;
    let n = the_i32_of(cols, "cols")?;
    let mut column_major = the_column_major_copy_of(a, rows, cols);
    let smallest = rows.min(cols);
    let mut values = vec![0.0_f64; smallest];
    // The routine asks for 8 integers for each of the smaller dimension.
    // The multiplication overflows for no matrix that reaches here:
    // `lib.rs` has refused every one of more than 2147483647 values, which
    // leaves the smaller dimension at 46340.
    let integers_wanted = smallest.checked_mul(8).ok_or_else(|| Error::Dimension {
        argument: "rows",
        expected: format!(
            "small enough for the workspace of integers of dgesdd, 8 for each of the smaller dimension, to fit in this machine, and the dimensions are {rows} and {cols}"
        ),
    })?;
    let mut integers = vec![0_i32; integers_wanted];
    // `jobz` N writes neither matrix of vectors, and their two leading
    // dimensions are passed as 1, the smallest the routine takes.
    let mut no_left_vectors = [0.0_f64; 1];
    let mut no_right_vectors = [0.0_f64; 1];
    let mut info = 0_i32;

    // The routine says how much it wants to work in when it is called
    // with the length of its workspace at -1, which is how LAPACK is
    // asked, and writes that number into the first entry of the workspace
    // it was given. It is what `eigh_lower` above asks `dsyevd` and what
    // the thin QR asks its two routines.
    let mut asked = [0.0_f64; 1];
    // SAFETY: with `lwork` at -1 the routine writes the first entry of
    // `work` and reads nothing else of it, and `asked` holds one value;
    // it reads and writes nothing of `a`, of `s`, of `u`, of `vt` or of
    // `info` other than to store that size, and `column_major` holds rows
    // * cols values, `values` holds min(rows, cols), the two buffers of
    // vectors hold one value each, `integers` holds 8 * min(rows, cols)
    // and `info` is one integer. Neither dimension is 0 and both fit in
    // the `i32` the routine takes, which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dgesdd(
            b'N',
            m,
            n,
            &mut column_major,
            m,
            &mut values,
            &mut no_left_vectors,
            1,
            &mut no_right_vectors,
            1,
            &mut asked,
            -1,
            &mut integers,
            &mut info,
        );
    }
    if info != 0 {
        return Err(Error::NoConvergence {
            routine: "dgesdd",
            info,
        });
    }

    let floats =
        the_workspace_of_the_singular_values(rows, cols, asked.first().copied().unwrap_or(0.0))?;
    let lwork = the_length_of_the_workspace_of_the_singular_values(floats)?;
    let mut work = vec![0.0_f64; floats];
    // SAFETY: with `jobz` N, `m` = rows, `n` = cols and `lda` = rows the
    // routine reads and overwrites `a` as a column major matrix of rows x
    // cols, which is the rows * cols values `column_major` holds; it
    // writes the min(rows, cols) singular values into `s`, which holds
    // that many; it writes neither `u` nor `vt`, whose leading dimensions
    // are 1 and whose buffers hold one value each; and it works in the
    // first `lwork` values of `work`, which holds exactly that many, and
    // in the 8 * min(rows, cols) of `iwork`, which is what `integers`
    // holds. It writes nothing else, and `info` is one integer. Neither
    // dimension is 0 and all three lengths fit in the `i32` the routine
    // takes, which `the_i32_of` and
    // `the_length_of_the_workspace_of_the_singular_values` have just
    // checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dgesdd(
            b'N',
            m,
            n,
            &mut column_major,
            m,
            &mut values,
            &mut no_left_vectors,
            1,
            &mut no_right_vectors,
            1,
            &mut work,
            lwork,
            &mut integers,
            &mut info,
        );
    }
    if info != 0 {
        return Err(Error::NoConvergence {
            routine: "dgesdd",
            info,
        });
    }
    Ok(values)
}

/// How many floats `dgesdd` works in for a matrix of `rows` x `cols` with
/// `jobz` `N`: the larger of what the query asked for and the minimum the
/// routine documents, 3m + max(M, 7m) for m the smaller dimension and M
/// the larger.
///
/// The query writes its number as an `f64`. An infinity, a NaN, a
/// negative number and one above the `i32` the length is passed as are
/// left out, and the minimum stands; a value that is a count is taken by
/// its whole part, the fraction that the routine cannot have meant being
/// dropped. It is what `the_workspace_of` below does with the query of
/// `dsyevd`.
///
/// # Errors
///
/// [`Error::Dimension`] when the minimum for those dimensions is more
/// values than this machine can hold.
fn the_workspace_of_the_singular_values(rows: usize, cols: usize, asked: f64) -> Result<usize> {
    let too_large = || Error::Dimension {
        argument: "rows",
        expected: format!(
            "small enough for the workspace of dgesdd, 3m + max(M, 7m) floats for the smaller dimension m and the larger M, to fit in this machine, and they are {rows} and {cols}"
        ),
    };
    let smaller = rows.min(cols);
    let larger = rows.max(cols);
    let at_least = 3_usize
        .checked_mul(smaller)
        .and_then(|values| values.checked_add(larger.max(7_usize.checked_mul(smaller)?)))
        .ok_or_else(too_large)?;

    let the_most_a_length_holds = f64::from(i32::MAX);
    let asked = if asked.is_finite() && asked >= 0.0 && asked <= the_most_a_length_holds {
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "the line above has checked that the value is finite, is not negative and is at most i32::MAX, so it is a count this machine holds; what the cast drops is the fraction, which the routine cannot have meant, and the minimum below stands when what is left is smaller than it"
        )]
        let asked = asked as usize;
        asked
    } else {
        0
    };
    Ok(asked.max(at_least))
}

/// The length of the workspace of `dgesdd` as the `i32` the routine takes
/// it as.
///
/// [`the_workspace_of_the_singular_values`] gives no number above that
/// `i32` from the query, and the minimum it can give instead is at most
/// 10 times 46340, since `lib.rs` refuses a matrix of more than 2147483647
/// values, so this refuses nothing that reaches it: it is how the
/// conversion is made without an `as` that could truncate in silence.
///
/// # Errors
///
/// [`Error::Dimension`] when the workspace is larger than that.
fn the_length_of_the_workspace_of_the_singular_values(values: usize) -> Result<i32> {
    i32::try_from(values).map_err(|_| Error::Dimension {
        argument: "rows",
        expected: format!(
            "small enough that the workspace dgesdd asks for, {values} values here, is at most the {largest} its length is passed as",
            largest = i32::MAX
        ),
    })
}

#[cfg(test)]
mod tests {
    use super::the_workspace_of;

    /// The minimum `dsyevd` documents for an n of 3: 1 + 6n + 2n² = 37
    /// floats and 3 + 5n = 18 integers.
    const THE_MINIMUM_FOR_3: (usize, usize) = (37, 18);

    #[test]
    fn a_workspace_the_query_did_not_give_a_number_for_is_the_minimum() {
        assert_eq!(
            the_workspace_of(3, f64::NAN, -1).unwrap(),
            THE_MINIMUM_FOR_3
        );
        assert_eq!(
            the_workspace_of(3, f64::INFINITY, 0).unwrap(),
            THE_MINIMUM_FOR_3
        );
    }

    #[test]
    fn a_workspace_the_query_gave_a_negative_number_for_is_the_minimum() {
        assert_eq!(the_workspace_of(3, -5.0, -7).unwrap(), THE_MINIMUM_FOR_3);
    }

    #[test]
    fn a_workspace_the_query_asked_more_for_is_what_it_asked() {
        assert_eq!(the_workspace_of(3, 100.0, 40).unwrap(), (100, 40));
    }

    #[test]
    fn a_workspace_the_query_gave_a_fraction_for_is_its_whole_part() {
        assert_eq!(the_workspace_of(3, 100.9, 40).unwrap(), (100, 40));
    }

    #[test]
    fn a_workspace_the_query_asked_more_than_a_length_holds_for_is_the_minimum() {
        let above_the_largest_length = f64::from(i32::MAX) * 2.0;
        assert_eq!(
            the_workspace_of(3, above_the_largest_length, 40).unwrap(),
            (37, 40)
        );
    }
}
