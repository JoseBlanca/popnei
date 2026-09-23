//! The BLAS and LAPACK backend: the routines of the library of the
//! system, `dsyrk`, `dgemm`, which the two products of this module call,
//! and `dsyevd`, the ones numpy calls.
//!
//! Every matrix reaches this module row after row, and these routines read
//! a matrix column after column. The buffer of an r x c matrix read that
//! way is its transpose, c x r, so each routine is called on the
//! transposes: the lower half of `g` in popnei's layout is the upper half
//! for the routine, so `uplo` is `U`; `a'a` is `a a'` of the transposed
//! view, so `trans` is `N`; and `c = a b` is `c' = b' a'`, so `dgemm` gets
//! the buffer of `b` as its first operand and that of `a` as its second.
//! `c = a b'` is `c' = b a'` the same way, and there the first operand is
//! the transpose of what the buffer of `b` is in the routine's view, so
//! that call is the one whose `transa` is `T`.
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
pub(crate) fn product_by_transpose(
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
