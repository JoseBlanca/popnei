//! The BLAS and LAPACK backend: the routines of the library of the
//! system, `dsyrk` and `dgemm`, which the four products of this module
//! call, and `dsyevd`, `dpotrf`, `dpotrs`, `dgeqrf`, `dorgqr`, `dtrtrs`
//! and `dgesdd`, which the other operations call. The inverse is the one
//! operation here with no routine of its own: `dpotri`, which LAPACK has
//! for it, is the one call of this backend that Accelerate computes
//! wrongly under concurrency, so the inverse is built from the triangular
//! solve and the product, as [`invert_with_cholesky`] says and as
//! "Calling the crate from two threads" of `docs/specs/linalg.md` has it.
//!
//! Two of those numpy does not call: `dpotrs` and `dtrtrs`, of which
//! `_umath_linalg` of numpy 2.5.3 exports neither, checked with `nm` on
//! 23 September 2026. numpy solves a system and inverts a matrix
//! through an LU of a general square matrix where this crate goes through
//! a Cholesky of a symmetric positive definite one, and it solves
//! `r c = q' y` with that same general solve where this crate has
//! `dtrtrs`. "Why a Cholesky where numpy uses an LU" of
//! `docs/specs/linalg.md` says why.
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

use crate::{Eigen, Error, Result, TheHalfThatHoldsTheMatrix};

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
    // exactly that many. Neither dimension is 0, which `lib.rs` refuses
    // above both backends, and both fit in the `i32` the routine takes,
    // which `the_i32_of` has just checked.
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
    // dimension is 0, which `lib.rs` refuses above both backends, and all
    // three fit in the `i32` the routine takes, which `the_i32_of` has
    // just checked.
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
    // inner first operand and which is the cols * inner values of `b`; the
    // buffer of `a` as one of inner rows and rows columns with `ldb` =
    // inner, which is the rows * inner values of `a`; and, with `beta` 0,
    // writes a column major matrix of cols rows and rows columns with
    // `ldc` = cols, which is the rows * cols values of `c`. Each of the
    // three slices holds exactly the values of its dimensions, no
    // dimension is 0, which `lib.rs` refuses above both backends, and all
    // three fit in the `i32` the routine takes, which `the_i32_of` has
    // just checked.
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
    // values of its dimensions, no dimension is 0, which `lib.rs` refuses
    // above both backends, and all three fit in the `i32` the routine
    // takes, which `the_i32_of` has just checked.
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
    // dimension is 0, which `lib.rs` refuses above both backends, and all
    // three fit in the `i32` the routine takes, which `the_i32_of` has
    // just checked.
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
    // `a`, and `info` is one integer. `n` is not 0, which `lib.rs` refuses
    // above both backends, and it fits in the `i32` the routine takes,
    // which `the_i32_of` has just checked.
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
    // `info` is one integer. Neither dimension is 0, which `lib.rs`
    // refuses above both backends, and both fit in the `i32` the routine
    // takes, which `the_i32_of` has just checked.
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
/// The inverse is built from the triangular solve and the product above
/// and not from LAPACK's `dpotri`, which is the one routine of this
/// backend that Accelerate computes wrongly: measured on this Mac on 24
/// September 2026, 2000 inversions of one 200 x 200 matrix gave between
/// 213 and 546 answers different from the rest beside seven threads doing
/// nothing but `product`, the worst by 2.4e-4 on entries that reach
/// 1.7e-2, and 2, 2, 1 and 0 different with nothing else of popnei
/// running, since Accelerate's own threads inside the routine are enough.
/// The route below gave 0 of 2000 in both cases. "Calling the crate from
/// two threads" of `docs/specs/linalg.md` has the whole of what was
/// measured, on which routines and at which sizes, and what this route
/// costs in time.
///
/// The route is the identity solved against `l`, whose answer is the
/// inverse of `l`, and then that answer times its own transpose: `a` is
/// `l l'`, so the inverse of `a` is `l⁻¹' l⁻¹`. The solve is given the
/// identity as `n` right hand sides and gives back one solution to a row,
/// which is `l⁻¹` held column after column, and that buffer read as this
/// crate reads one is `l⁻¹'`, so the product wanted is that buffer times
/// its transpose, which is [`product_with_the_second_turned`].
///
/// It costs two buffers of `n` x `n` values, 320 KB each for the 200
/// individuals of popnei's panels, that `dpotri` did not need: the solve
/// overwrites the right hand sides it is given, so the identity is a
/// buffer of its own, and the product may not write where it reads, nor
/// into `inverse`, whose upper half is the caller's.
///
/// # Errors
///
/// [`Error::Dimension`] when `n` is larger than the `i32` the routines
/// take. [`Error::Memory`] when this machine has not the memory for
/// either buffer. [`Error::Singular`] when the diagonal of `l` holds a 0
/// at the row the error names, which the solve would divide by; `lib.rs`
/// reads that diagonal before either backend runs, since faer's inverse
/// does not look at it, so no caller of the crate reaches this one.
/// [`Error::NoConvergence`] when a routine refused an argument it was
/// given, which is a defect of popnei.
pub(crate) fn invert_with_cholesky(l: &[f64], n: usize, inverse: &mut [f64]) -> Result<()> {
    // The dimension is checked here and not left to the calls below so
    // that the error names `n`, which is the argument the caller passed,
    // and not the `sides` the solve is given it as.
    the_i32_of(n, "n")?;
    // `l` holds exactly n * n values, `lib.rs` having cut it to the
    // dimensions, so its length is the length of each buffer here.
    let mut turned = the_buffer_of(l.len(), "the identity the inverse is solved against")?;
    for (row, entries) in turned.chunks_exact_mut(n).enumerate() {
        if let Some(entry) = entries.get_mut(row) {
            *entry = 1.0;
        }
    }
    if let Err(error) = solve_triangular(
        l,
        n,
        TheHalfThatHoldsTheMatrix::TheLowerHalf,
        &mut turned,
        n,
    ) {
        // The solve names the triangular matrix `a`, which here is the
        // factorization the caller gave as `l`.
        if let Error::Singular { at, .. } = &error {
            return Err(Error::Singular {
                argument: "l",
                at: *at,
            });
        }
        return Err(error);
    }
    let mut whole = the_buffer_of(l.len(), "the inverse the product of the solve writes")?;
    product_with_the_second_turned(&turned, n, n, &turned, n, &mut whole)?;
    // Both buffers hold exactly n * n values, so each is n rows of n, and
    // the lower half is the entries of column `j` at most `i` of row `i`.
    // What is left of the row after those is the caller's and is not
    // written.
    for (row, (into, from)) in inverse
        .chunks_exact_mut(n)
        .zip(whole.chunks_exact(n))
        .enumerate()
    {
        for (into, from) in into.iter_mut().zip(from).take(row.saturating_add(1)) {
            *into = *from;
        }
    }
    Ok(())
}

/// A buffer of `values` floats, every one 0, whose memory is asked for and
/// not taken, as [`the_column_major_copy_of`] asks for its own.
///
/// # Errors
///
/// [`Error::Memory`] when this machine has not the memory for it.
fn the_buffer_of(values: usize, what: &'static str) -> Result<Vec<f64>> {
    let mut buffer: Vec<f64> = Vec::new();
    buffer
        .try_reserve_exact(values)
        .map_err(|_| Error::Memory { what, values })?;
    buffer.resize(values, 0.0);
    Ok(buffer)
}

/// The thin QR of `a`, of exactly `rows` x `cols` values row after row
/// with `rows` at least `cols` and `cols` 1 at least: `q` of exactly
/// `rows` x `cols` values and the upper triangular `r` of exactly `cols`
/// x `cols`, both written row after row and the lower half of `r` set
/// to 0.
///
/// The buffer of `a` read column after column is the `cols` x `rows`
/// matrix, the wide one, and the two routines are much slower on it than
/// on the tall one: for a design of 10000 x 5 "What the seven of the GWAS
/// cost" of `docs/specs/linalg.md` measured the routines at 2.98 ms that
/// way and 0.165 ms through the copy below, on 23 September 2026. Those
/// are the routines alone and not this function, which adds the checks of
/// `lib.rs` and its allocations and which three runs on the same machine
/// gave 0.1573, 0.1855 and 0.207 ms. So this backend writes the transpose
/// of `a` into a buffer of its own, `rows` x `cols` values and 400 KB at
/// that size, which is `a` held column after column, and calls them on
/// that.
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
/// [`Error::Memory`] when this machine has not the memory for the column
/// major copy. [`Error::NoConvergence`] when a routine refused an
/// argument it was given, which is a defect of popnei: neither has
/// another reason to give an `info` other than 0.
pub(crate) fn thin_qr(
    a: &[f64],
    rows: usize,
    cols: usize,
    q: &mut [f64],
    r: &mut [f64],
) -> Result<()> {
    let m = the_i32_of(rows, "rows")?;
    let n = the_i32_of(cols, "cols")?;
    let mut column_major = the_column_major_copy_of(
        a,
        rows,
        cols,
        "the column major copy of a for dgeqrf and dorgqr",
    )?;
    // One coefficient for each column, which is what the factorization
    // keeps beside the vectors it leaves in the matrix.
    let mut coefficients = vec![0.0_f64; cols];
    let mut info = 0_i32;

    // Each routine says how much it wants to work in when it is called
    // with the length of its workspace at -1, which is how LAPACK is
    // asked, and writes that number into the first entry of the workspace
    // it was given. It is what `eigh_lower` above asks `dsyevd`.
    let mut asked = [0.0_f64; 1];
    // SAFETY: with `lwork` at -1 the routine factors nothing. It writes
    // the size it wants into the first entry of `work`, which is `asked`
    // and holds one value, and it writes `info`, which is one integer; it
    // reads and writes nothing of `a` and nothing of `tau`. Those two are
    // still passed as the buffers of the call below, with `m` = rows, `n`
    // = cols and `lda` = rows: `column_major` holds the rows * cols values
    // of a column major matrix of rows x cols and `coefficients` holds the
    // cols of `tau`. Neither dimension is 0, which `lib.rs` refuses above
    // both backends, and both fit in the `i32` the routine takes, which
    // `the_i32_of` has just checked.
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
    // SAFETY: with `lwork` at -1 the routine builds no `q`. It writes the
    // size it wants into the first entry of `work`, which is `asked` and
    // holds one value, and it writes `info`, which is one integer; it
    // reads and writes nothing of `a` and nothing of `tau`. Those two are
    // still passed as the buffers of the call below, with `m` = rows,
    // `n` = cols, `k` = cols and `lda` = rows: `column_major` holds the
    // rows * cols values of a column major matrix of rows x cols and
    // `coefficients` holds the cols of `tau`. No dimension is 0, `k` is
    // at most `n` and `n` at most `m` since rows is at least cols, and
    // all three fit in the `i32` the routine takes, which `the_i32_of`
    // has already checked.
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
/// what the routines of the thin QR and of the singular values are fast
/// on. The buffer holds exactly `rows` x `cols` values and `rows` and
/// `cols` are 1 at least. `what` names the copy in the error, as the
/// routines that are given it.
///
/// The column `j` of a matrix held row after row is one value in every
/// `cols` from the value `j`, and those are the values of the column `j`
/// of the copy, one after another.
///
/// The copy is as large as the matrix, 400 KB for a design of 10000 x 5
/// and 17 GB for the largest matrix this crate takes, so it is asked for
/// and not taken with `vec!`, which ends the process on a machine that
/// has not the memory, as the workspaces of the eigendecomposition above
/// are. A caller that could hold the matrix can usually hold the copy, so
/// no call the association study makes reaches that error, and
/// [`thin_qr`] and `rank` are public and nothing bounds what else is
/// passed to them.
///
/// # Errors
///
/// [`Error::Memory`] when this machine has not the memory for the copy.
fn the_column_major_copy_of(
    a: &[f64],
    rows: usize,
    cols: usize,
    what: &'static str,
) -> Result<Vec<f64>> {
    let mut column_major: Vec<f64> = Vec::new();
    column_major
        .try_reserve_exact(a.len())
        .map_err(|_| Error::Memory {
            what,
            values: a.len(),
        })?;
    column_major.resize(a.len(), 0.0);
    for (column, into) in column_major.chunks_exact_mut(rows).enumerate() {
        for (into, from) in into.iter_mut().zip(a.iter().skip(column).step_by(cols)) {
            *into = *from;
        }
    }
    Ok(column_major)
}

/// Writes into `into`, row after row and `cols` to a row, as many first
/// rows of the matrix that `column_major` holds column after column with
/// `rows` rows as `into` has room for: the whole of it for a `q` of
/// `rows` x `cols`, and the leading `cols` x `cols` corner for an `r`.
///
/// The row `i` of a matrix held column after column with `rows` rows is
/// one value in every `rows` from the value `i`, and there are `cols` of
/// them for every `i` below `rows`, which is what fills each row of
/// `into`.
///
/// `cols` is at most `rows`, so `into` has no row that the buffer has no
/// `i` for. Above that the `r` of `cols` x `cols` would have some, the
/// `zip` below would stop at the shorter side, and those rows would keep
/// the zeros they were allocated with and nothing would say so.
/// [`thin_qr`] refuses a matrix of fewer rows than columns before either
/// of its calls to this.
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

/// The `x` of `a x = b` for the triangular `a` of exactly `n` x `n`
/// values row after row, held in the half of that buffer `half` names,
/// and `b` of exactly `sides` x `n` values row after row, one row for each
/// right hand side, which comes back holding the solutions the same way.
/// `n` and `sides` are 1 at least.
///
/// The buffer of `a` read column after column is its transpose, so the
/// half popnei was given is the other half in the routine's view: `uplo`
/// is L for popnei's upper half and U for its lower one. `trans` is T,
/// which asks the routine for the solve against the transpose of what it
/// read, so the two turns undo each other and what is solved against is
/// popnei's `a`. The buffer of `b` read that way is the n x `sides` matrix
/// whose columns are the right hand sides, which is what the routine
/// takes, as it is for `dpotrs` above, so nothing is copied here either.
///
/// # Errors
///
/// [`Error::Dimension`] when a dimension is larger than the `i32` the
/// routine takes. [`Error::Singular`] when the diagonal of `a` holds a 0
/// at the row the error names, which the routine would divide by;
/// `lib.rs` reads that diagonal before either backend runs, since faer
/// does not look at it, so no caller of the crate reaches this one.
/// [`Error::NoConvergence`] when the routine refused an argument it was
/// given, which is a defect of popnei.
pub(crate) fn solve_triangular(
    a: &[f64],
    n: usize,
    half: TheHalfThatHoldsTheMatrix,
    b: &mut [f64],
    sides: usize,
) -> Result<()> {
    let order = the_i32_of(n, "n")?;
    let right_hand_sides = the_i32_of(sides, "sides")?;
    // The half of popnei's buffer is the other half of the routine's
    // column major view, as the doc comment above says, so popnei's upper
    // half is the routine's lower triangle and its lower half the
    // routine's upper one.
    let uplo = match half {
        TheHalfThatHoldsTheMatrix::TheUpperHalf => b'L',
        TheHalfThatHoldsTheMatrix::TheLowerHalf => b'U',
    };
    let mut info = 0_i32;
    // SAFETY: with `trans` T, `diag` N, `n` = n and `lda` = n the routine
    // reads the triangle `uplo` names of `a` as a column major matrix of n
    // x n, which is the half of `a` the caller named in popnei's layout
    // and is inside the n * n values `a` holds; and with `nrhs` = sides
    // and `ldb` = n it reads and writes `b` as a column major matrix of n
    // rows and sides columns, which is the sides * n values `b` holds. It
    // writes nothing else, and `info` is one integer. Neither dimension is
    // 0, which `lib.rs` refuses above both backends, and both fit in the
    // `i32` the routine takes, which `the_i32_of` has just checked.
    #[expect(
        unsafe_code,
        reason = "the routines of LAPACK are declared as unsafe functions over slices whose lengths nothing checks against the dimensions, which is why they are called here and nowhere else in popnei"
    )]
    unsafe {
        ::lapack::dtrtrs(
            uplo,
            b'T',
            b'N',
            order,
            right_hand_sides,
            a,
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
            argument: "a",
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
/// tall one: for a design of 10000 x 5 "What the seven of the GWAS cost"
/// of `docs/specs/linalg.md` measured the routine at 1.54 ms that way and
/// 0.145 ms through the copy below, on 23 September 2026. Those are the
/// routine alone and not this function, which adds the checks of `lib.rs`
/// and its allocations and which measured 0.1881 ms on the same machine.
/// So this backend writes the transpose of `a` into a buffer of its own,
/// which is `a` held column after column, and calls the routine on that,
/// as the thin QR above does. The routine overwrites that buffer, and `a`
/// itself is left as it was.
///
/// `jobz` is `N`, which computes the values and neither of the two
/// matrices of vectors, so the two buffers they would go in hold one
/// value each and are never written.
///
/// # Errors
///
/// [`Error::Dimension`] when a dimension, or the workspace the routine
/// asks for, is larger than the `i32` the routine takes.
/// [`Error::Memory`] when this machine has not the memory for the column
/// major copy. [`Error::NoConvergence`] when the routine gave an `info`
/// other than 0, which is either a decomposition that did not come out or
/// an argument it refused, the second being a defect of popnei.
pub(crate) fn singular_values(a: &[f64], rows: usize, cols: usize) -> Result<Vec<f64>> {
    let m = the_i32_of(rows, "rows")?;
    let n = the_i32_of(cols, "cols")?;
    let mut column_major =
        the_column_major_copy_of(a, rows, cols, "the column major copy of a for dgesdd")?;
    let smallest = rows.min(cols);
    let mut values = vec![0.0_f64; smallest];
    // The routine asks for 8 integers for each of the smaller dimension.
    // The multiplication overflows for no matrix that reaches here:
    // `lib.rs` has refused every one of more than 2147483647 values, which
    // leaves the smaller dimension at 46340.
    let integers_wanted = smallest.checked_mul(8).ok_or_else(|| Error::Dimension {
        argument: the_smaller_dimension_of(rows, cols),
        expected: format!(
            "small enough for the workspace of integers of dgesdd, 8 for each of the smaller dimension, to fit in this machine, and the dimensions are {rows} rows and {cols} columns"
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
    // SAFETY: with `lwork` at -1 the routine decomposes nothing. It writes
    // the size it wants into the first entry of `work`, which is `asked`
    // and holds one value, and it writes `info`, which is one integer; it
    // reads and writes nothing of `a`, of `s`, of `u`, of `vt` or of
    // `iwork`. Those are still passed as the buffers of the call below,
    // with `jobz` N, `m` = rows, `n` = cols, `lda` = rows and the two
    // leading dimensions of the vectors 1: `column_major` holds the rows *
    // cols values of a column major matrix of rows x cols, `values` holds
    // the min(rows, cols) of `s`, the two buffers of vectors hold one
    // value each and `integers` holds the 8 * min(rows, cols) of `iwork`.
    // Neither dimension is 0, which `lib.rs` refuses above both backends,
    // and both fit in the `i32` the routine takes, which `the_i32_of` has
    // just checked.
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
    let lwork = the_length_of_the_workspace_of_the_singular_values(
        floats,
        the_larger_dimension_of(rows, cols),
    )?;
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
        argument: the_larger_dimension_of(rows, cols),
        expected: format!(
            "small enough for the workspace of dgesdd, 3m + max(M, 7m) floats for the smaller dimension m and the larger M, to fit in this machine, and they are {rows} rows and {cols} columns"
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
/// it as. `argument` is the dimension that set it, which is the larger of
/// the two.
///
/// This one does refuse matrices that reach it. What
/// [`the_workspace_of_the_singular_values`] takes from the query is at
/// most that `i32`, but the minimum it gives instead is 3m + max(M, 7m)
/// for the smaller dimension m and the larger M, and `lib.rs` lets
/// through every matrix of at most 2147483647 values, the one row of
/// 2147483647 columns among them, whose minimum is 2147483650 and is
/// above what the length holds. So the conversion is made here and not
/// with an `as` that would truncate it in silence.
///
/// # Errors
///
/// [`Error::Dimension`] when the workspace is larger than that.
fn the_length_of_the_workspace_of_the_singular_values(
    values: usize,
    argument: &'static str,
) -> Result<i32> {
    i32::try_from(values).map_err(|_| Error::Dimension {
        argument,
        expected: format!(
            "small enough that the workspace dgesdd asks for, {values} values here, is at most the {largest} its length is passed as",
            largest = i32::MAX
        ),
    })
}

/// The name of the larger of the two dimensions of a matrix, which is the
/// one at fault when the workspace of floats of `dgesdd`, 3m + max(M, 7m),
/// does not fit: `rows` when the two are equal.
fn the_larger_dimension_of(rows: usize, cols: usize) -> &'static str {
    if rows >= cols { "rows" } else { "cols" }
}

/// The name of the smaller of the two, which is the one at fault when the
/// workspace of integers of `dgesdd`, 8 for each of it, does not fit:
/// `rows` when the two are equal.
fn the_smaller_dimension_of(rows: usize, cols: usize) -> &'static str {
    if rows <= cols { "rows" } else { "cols" }
}

#[cfg(test)]
mod tests {
    use super::{
        add_self_product_lower, cholesky_lower, invert_with_cholesky, product, solve_with_cholesky,
        the_workspace_of, the_workspace_of_the_singular_values, the_workspace_of_the_thin_qr,
    };
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, Ordering};

    /// The individuals of the covariance the inverse is measured on, which
    /// is the 200 of popnei's panels.
    const THE_INDIVIDUALS: usize = 200;

    /// The rows of the matrix that covariance is the self product of, five
    /// for each individual, which is what leaves it well conditioned.
    const THE_ROWS_OF_THE_COVARIANCE: usize = 1000;

    /// How many inversions the test below makes beside the other threads.
    const THE_INVERSIONS: usize = 400;

    /// How many threads do nothing but products while it inverts.
    const THE_THREADS_IN_BLAS: usize = 7;

    /// The covariance the inverse is measured on: A'A for the A of 1000
    /// rows and 200 columns of the generator of "How the seven are
    /// verified" of `docs/specs/linalg.md`, a 200 x 200 symmetric positive
    /// definite matrix whose condition number is 6.5, so that nothing the
    /// test sees comes from a matrix that is hard to invert. Only its
    /// lower half is filled, which is what a Cholesky and an inverse read.
    fn the_covariance_of_the_generator() -> Vec<f64> {
        let mut state = 7_u64;
        let a: Vec<f64> = (0..THE_ROWS_OF_THE_COVARIANCE * THE_INDIVIDUALS)
            .map(|_| {
                state ^= state.wrapping_shl(13);
                state ^= state.wrapping_shr(7);
                state ^= state.wrapping_shl(17);
                // 2^53, below which a count is an exact f64.
                state.wrapping_shr(11) as f64 / 9007199254740992.0 - 0.5
            })
            .collect();
        let mut covariance = vec![0.0_f64; THE_INDIVIDUALS * THE_INDIVIDUALS];
        add_self_product_lower(
            &a,
            THE_ROWS_OF_THE_COVARIANCE,
            THE_INDIVIDUALS,
            &mut covariance,
        )
        .unwrap();
        covariance
    }

    /// The identity of 200 x 200, row after row, which is the 200 right
    /// hand sides a solve inverts a 200 x 200 matrix with. Its diagonal is
    /// one entry in every 201 of the buffer, a row and a column further on
    /// each time.
    fn the_identity_of_the_individuals() -> Vec<f64> {
        let mut identity = vec![0.0_f64; THE_INDIVIDUALS * THE_INDIVIDUALS];
        for entry in identity.iter_mut().step_by(THE_INDIVIDUALS + 1) {
            *entry = 1.0;
        }
        identity
    }

    /// How far the lower halves of two `n` x `n` matrices are apart, the
    /// largest difference of an entry of column `j` at most `i` in row
    /// `i`.
    fn the_largest_difference_of_the_lower_halves(one: &[f64], other: &[f64], n: usize) -> f64 {
        let mut largest = 0.0_f64;
        for (row, (entries, others)) in one.chunks_exact(n).zip(other.chunks_exact(n)).enumerate() {
            for (entry, other) in entries.iter().zip(others).take(row.saturating_add(1)) {
                largest = largest.max((entry - other).abs());
            }
        }
        largest
    }

    /// What this guards: on Accelerate, LAPACK's `dpotri`, which inverts a
    /// Cholesky factorization in one call, gives a numerically wrong
    /// answer, most often while another call of Accelerate is running on
    /// another thread of the same process, and `invert_with_cholesky` is
    /// written as a triangular solve and a product for that reason. The
    /// test inverts one matrix again and again while other threads do
    /// nothing but products, and every answer has to be the answer the
    /// same route gave before any of those threads was started, bit for
    /// bit.
    ///
    /// Measured on this Mac on 24 September 2026, on this matrix and with
    /// these seven threads: the route through `dpotri` gave between 213
    /// and 546 answers different from the rest in seven runs of 2000
    /// inversions, the worst off by 2.4e-4 on entries that reach 1.7e-2,
    /// and the route the backend has now gave 0 of 2000 in five runs. At a
    /// rate of about 1 in 8 the 400 inversions here would have to miss 400
    /// times running to pass on the old route, and 20 runs of the test
    /// against that route failed 20 times. The rate depends on the size,
    /// which is why the test is at the 200 of popnei's panels and not at a
    /// size that would make it quicker: with everything else the same,
    /// 2000 inversions of a 50 x 50 and of a 100 x 100 gave 0 wrong, of a
    /// 400 x 400, 27, and of an 800 x 800, 1127.
    ///
    /// The reference the answers are compared with is taken before the
    /// other threads are started, so it is what one thread computes. That
    /// the route gives the right matrix at all is what the comparison with
    /// `solve_with_cholesky` checks, which inverts through `dpotrs` and
    /// the identity instead and is a different routine of LAPACK.
    #[test]
    fn the_inverse_is_the_same_while_other_threads_are_in_blas() {
        let mut l = the_covariance_of_the_generator();
        cholesky_lower(&mut l, THE_INDIVIDUALS).unwrap();
        let mut reference = vec![0.0_f64; THE_INDIVIDUALS * THE_INDIVIDUALS];
        invert_with_cholesky(&l, THE_INDIVIDUALS, &mut reference).unwrap();

        // The same inverse through `dpotrs`: the solution of a x = i for
        // the identity is the inverse of a, and the buffer comes back with
        // one solution to a row, which for a symmetric inverse is the
        // inverse row after row.
        let mut through_the_solve = the_identity_of_the_individuals();
        solve_with_cholesky(&l, THE_INDIVIDUALS, &mut through_the_solve, THE_INDIVIDUALS).unwrap();
        let apart = the_largest_difference_of_the_lower_halves(
            &reference,
            &through_the_solve,
            THE_INDIVIDUALS,
        );
        assert!(
            apart < 1e-15,
            "the inverse and the solve against the identity are {apart:e} apart, and the two routes measured 3.1e-17 on 24 September 2026"
        );

        let stop = Arc::new(AtomicBool::new(false));
        let mut threads = Vec::new();
        for _ in 0..THE_THREADS_IN_BLAS {
            let stop = Arc::clone(&stop);
            threads.push(std::thread::spawn(move || {
                let a = the_covariance_of_the_generator();
                let mut c = vec![0.0_f64; THE_INDIVIDUALS * THE_INDIVIDUALS];
                while !stop.load(Ordering::Relaxed) {
                    product(
                        &a,
                        THE_INDIVIDUALS,
                        THE_INDIVIDUALS,
                        &a,
                        THE_INDIVIDUALS,
                        &mut c,
                    )
                    .unwrap();
                }
            }));
        }
        let mut wrong = 0_usize;
        let mut worst = 0.0_f64;
        let mut inverse = vec![0.0_f64; THE_INDIVIDUALS * THE_INDIVIDUALS];
        for _ in 0..THE_INVERSIONS {
            invert_with_cholesky(&l, THE_INDIVIDUALS, &mut inverse).unwrap();
            let apart =
                the_largest_difference_of_the_lower_halves(&reference, &inverse, THE_INDIVIDUALS);
            if apart > 0.0 {
                wrong = wrong.saturating_add(1);
                worst = worst.max(apart);
            }
        }
        stop.store(true, Ordering::Relaxed);
        for thread in threads {
            thread.join().unwrap();
        }
        assert_eq!(
            wrong, 0,
            "{wrong} of {THE_INVERSIONS} inversions differ from the one the same route gave on one thread, the worst by {worst:e}"
        );
    }

    /// The minimum `dgeqrf` and `dorgqr` document for a matrix of 3
    /// columns, which is those 3 floats.
    const THE_MINIMUM_OF_THE_THIN_QR_FOR_3: usize = 3;

    #[test]
    fn a_workspace_of_the_thin_qr_the_query_did_not_give_a_number_for_is_the_minimum() {
        assert_eq!(
            the_workspace_of_the_thin_qr(3, f64::NAN),
            THE_MINIMUM_OF_THE_THIN_QR_FOR_3
        );
        assert_eq!(
            the_workspace_of_the_thin_qr(3, f64::INFINITY),
            THE_MINIMUM_OF_THE_THIN_QR_FOR_3
        );
    }

    #[test]
    fn a_workspace_of_the_thin_qr_the_query_gave_a_negative_number_for_is_the_minimum() {
        assert_eq!(
            the_workspace_of_the_thin_qr(3, -5.0),
            THE_MINIMUM_OF_THE_THIN_QR_FOR_3
        );
    }

    #[test]
    fn a_workspace_of_the_thin_qr_the_query_asked_more_for_is_what_it_asked() {
        assert_eq!(the_workspace_of_the_thin_qr(3, 100.0), 100);
    }

    #[test]
    fn a_workspace_of_the_thin_qr_the_query_gave_a_fraction_for_is_its_whole_part() {
        assert_eq!(the_workspace_of_the_thin_qr(3, 100.9), 100);
    }

    #[test]
    fn a_workspace_of_the_thin_qr_the_query_asked_more_than_a_length_holds_for_is_the_minimum() {
        let above_the_largest_length = f64::from(i32::MAX) * 2.0;
        assert_eq!(
            the_workspace_of_the_thin_qr(3, above_the_largest_length),
            THE_MINIMUM_OF_THE_THIN_QR_FOR_3
        );
    }

    /// The dimensions the workspace of `dgesdd` is asked for at, 4 rows
    /// and 2 columns, and the minimum it documents for them with `jobz`
    /// `N`: 3m + max(M, 7m) for m the smaller dimension and M the larger,
    /// which is 6 + max(4, 14) = 20.
    const THE_DIMENSIONS_OF_THE_SINGULAR_VALUES: (usize, usize) = (4, 2);
    const THE_MINIMUM_OF_THE_SINGULAR_VALUES: usize = 20;

    #[test]
    fn a_workspace_of_the_singular_values_the_query_did_not_give_a_number_for_is_the_minimum() {
        let (rows, cols) = THE_DIMENSIONS_OF_THE_SINGULAR_VALUES;
        assert_eq!(
            the_workspace_of_the_singular_values(rows, cols, f64::NAN).unwrap(),
            THE_MINIMUM_OF_THE_SINGULAR_VALUES
        );
        assert_eq!(
            the_workspace_of_the_singular_values(rows, cols, f64::INFINITY).unwrap(),
            THE_MINIMUM_OF_THE_SINGULAR_VALUES
        );
    }

    #[test]
    fn a_workspace_of_the_singular_values_the_query_gave_a_negative_number_for_is_the_minimum() {
        let (rows, cols) = THE_DIMENSIONS_OF_THE_SINGULAR_VALUES;
        assert_eq!(
            the_workspace_of_the_singular_values(rows, cols, -5.0).unwrap(),
            THE_MINIMUM_OF_THE_SINGULAR_VALUES
        );
    }

    #[test]
    fn a_workspace_of_the_singular_values_the_query_asked_more_for_is_what_it_asked() {
        let (rows, cols) = THE_DIMENSIONS_OF_THE_SINGULAR_VALUES;
        assert_eq!(
            the_workspace_of_the_singular_values(rows, cols, 100.0).unwrap(),
            100
        );
    }

    #[test]
    fn a_workspace_of_the_singular_values_the_query_gave_a_fraction_for_is_its_whole_part() {
        let (rows, cols) = THE_DIMENSIONS_OF_THE_SINGULAR_VALUES;
        assert_eq!(
            the_workspace_of_the_singular_values(rows, cols, 100.9).unwrap(),
            100
        );
    }

    #[test]
    fn a_workspace_of_the_singular_values_the_query_asked_more_than_a_length_holds_for_is_the_minimum()
     {
        let (rows, cols) = THE_DIMENSIONS_OF_THE_SINGULAR_VALUES;
        let above_the_largest_length = f64::from(i32::MAX) * 2.0;
        assert_eq!(
            the_workspace_of_the_singular_values(rows, cols, above_the_largest_length).unwrap(),
            THE_MINIMUM_OF_THE_SINGULAR_VALUES
        );
    }

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
