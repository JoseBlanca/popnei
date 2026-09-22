//! The BLAS and LAPACK backend: the three routines of the library of the
//! system, `dsyrk`, `dgemm` and `dsyevd`, the ones numpy calls.
//!
//! Every matrix reaches this module row after row, and these routines read
//! a matrix column after column. The buffer of an r x c matrix read that
//! way is its transpose, c x r, so each routine is called on the
//! transposes: the lower half of `g` in popnei's layout is the upper half
//! for the routine, so `uplo` is `U`; `a'a` is `a a'` of the transposed
//! view, so `trans` is `N`; and `c = a b` is `c' = b' a'`, so `dgemm` gets
//! the buffer of `b` as its first operand and that of `a` as its second.
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

/// The eigendecomposition of the symmetric `g`, of exactly `n` x `n`
/// values row after row with its lower half filled and `n` 1 at least.
///
/// # Errors
///
/// [`Error::Dimension`] when `n`, or the workspace the routine asks for,
/// is larger than the `i32` the routine takes.
/// [`Error::NoConvergence`] when the routine gave an `info` other than 0.
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
    let lwork = the_i32_of(floats, "the workspace of dsyevd")?;
    let liwork = the_i32_of(integers, "the workspace of dsyevd")?;
    let mut work = vec![0.0_f64; floats];
    let mut iwork = vec![0_i32; integers];
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

    // The routine gives the eigenvalues from the smallest, and each
    // eigenvector as a column of the matrix it read column after column,
    // which in popnei's buffer is a row. So the values are turned round
    // and the rows with them.
    values.reverse();
    reverse_the_rows(&mut g, n);
    Ok(Eigen { values, vectors: g })
}

/// Swaps row 0 of the n x n matrix with row n - 1, row 1 with row n - 2,
/// and so on, in the buffer it was given.
fn reverse_the_rows(values: &mut [f64], n: usize) {
    let Some(the_top_half) = (n / 2).checked_mul(n) else {
        return;
    };
    let Some((front, back)) = values.split_at_mut_checked(the_top_half) else {
        return;
    };
    for (from_the_top, from_the_bottom) in front.chunks_exact_mut(n).zip(back.rchunks_exact_mut(n))
    {
        from_the_top.swap_with_slice(from_the_bottom);
    }
}

/// How many floats and how many integers `dsyevd` works in for a matrix of
/// `n` x `n` with its eigenvectors: the larger of what the query asked
/// for and the minimum the routine documents, 1 + 6n + 2n² floats and
/// 3 + 5n integers.
///
/// The query writes the number of floats as an `f64`, and a value that is
/// not a whole count, an infinity, a NaN, a negative number or one above
/// the `i32` the length is passed as, is left out and the minimum stands.
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
            reason = "the line above has checked that the value is finite, is not negative and is at most i32::MAX, so it is a count this machine holds; the fraction the routine cannot have written is dropped, and the minimum below stands when that leaves too little"
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
