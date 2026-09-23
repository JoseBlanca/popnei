//! The faer backend: the operations of this crate on faer, a linear
//! algebra library written in Rust, which runs where there is no BLAS to
//! link and natively when the crate is built with `--no-default-features`.
//!
//! faer reads a matrix the way it is told to, so here the buffers are
//! given to it as what they are, row after row, and its functions are
//! called as written: no transpose and no turn of the halves, which is
//! what the BLAS backend beside this file needs. The transposes here are
//! the operands that the caller holds the other way round, the second of
//! `a b'`, the first of `a' b` and both of `a' b'`, each of them a matrix
//! reference over the same values and not a copy.
//!
//! The functions here are given slices whose lengths the caller has
//! already cut to the dimensions, and they check nothing else: the checks
//! are in `lib.rs`, where they hold for whichever backend runs. Nothing
//! here is `unsafe`.

use faer::dyn_stack::{MemBuffer, MemStack, StackReq};
use faer::linalg::cholesky::llt::factor::{
    LltError, LltRegularization, cholesky_in_place, cholesky_in_place_scratch,
};
use faer::linalg::cholesky::llt::inverse::{
    inverse as inverse_with_cholesky_of_faer, inverse_scratch,
};
use faer::linalg::cholesky::llt::solve::{solve_in_place_scratch, solve_in_place_with_conj};
use faer::linalg::evd::EvdError;
use faer::linalg::matmul::matmul;
use faer::linalg::matmul::triangular::{BlockStructure, matmul as triangular_matmul};
use faer::linalg::triangular_solve::solve_upper_triangular_in_place;
use faer::{Accum, Conj, MatMut, MatRef, Par, Side};

use crate::{Eigen, Error, Result};

/// The threads a product runs on: the global pool of rayon, the one every
/// parallel loop of popnei runs on, which `RAYON_NUM_THREADS` sizes; and
/// one thread in WebAssembly, where there are none and faer's `rayon`
/// feature is off. It is the rule faer's own default follows, which is
/// what the eigendecomposition below takes, since that one reads it
/// itself and takes no argument for it.
#[cfg(not(target_family = "wasm"))]
fn the_threads() -> Par {
    Par::rayon(0)
}

/// The same in WebAssembly.
#[cfg(target_family = "wasm")]
fn the_threads() -> Par {
    Par::Seq
}

/// Adds the lower half of `a'a` to `g`, with `a` of exactly `rows` x
/// `cols` values and `g` of exactly `cols` x `cols`, both row after row
/// and `rows` and `cols` 1 at least.
///
/// # Errors
///
/// None: faer refuses nothing that the checks of `lib.rs` let through.
/// The signature is the one of the BLAS backend, which does.
#[expect(
    clippy::unnecessary_wraps,
    reason = "the two backends have the same signature, and the BLAS one fails when a dimension is larger than the i32 its routines take"
)]
pub(crate) fn add_self_product_lower(
    a: &[f64],
    rows: usize,
    cols: usize,
    g: &mut [f64],
) -> Result<()> {
    let a = MatRef::from_row_major_slice(a, rows, cols);
    let g = MatMut::from_row_major_slice_mut(g, cols, cols);
    triangular_matmul(
        g,
        BlockStructure::TriangularLower,
        Accum::Add,
        a.transpose(),
        BlockStructure::Rectangular,
        a,
        BlockStructure::Rectangular,
        1.0,
        the_threads(),
    );
    Ok(())
}

/// Writes `a b` into `c`, with `a` of exactly `rows` x `inner` values, `b`
/// of `inner` x `cols` and `c` of `rows` x `cols`, all row after row and
/// every dimension 1 at least.
///
/// # Errors
///
/// None, as for the product above.
#[expect(
    clippy::unnecessary_wraps,
    reason = "the two backends have the same signature, and the BLAS one fails when a dimension is larger than the i32 its routines take"
)]
pub(crate) fn product(
    a: &[f64],
    rows: usize,
    inner: usize,
    b: &[f64],
    cols: usize,
    c: &mut [f64],
) -> Result<()> {
    let a = MatRef::from_row_major_slice(a, rows, inner);
    let b = MatRef::from_row_major_slice(b, inner, cols);
    let c = MatMut::from_row_major_slice_mut(c, rows, cols);
    matmul(c, Accum::Replace, a, b, 1.0, the_threads());
    Ok(())
}

/// Writes `a b'` into `c`, with `a` of exactly `rows` x `inner` values,
/// `b` of `cols` x `inner` and `c` of `rows` x `cols`, all row after row
/// and every dimension 1 at least.
///
/// # Errors
///
/// None, as for the two products above.
#[expect(
    clippy::unnecessary_wraps,
    reason = "the two backends have the same signature, and the BLAS one fails when a dimension is larger than the i32 its routines take"
)]
pub(crate) fn product_with_the_second_turned(
    a: &[f64],
    rows: usize,
    inner: usize,
    b: &[f64],
    cols: usize,
    c: &mut [f64],
) -> Result<()> {
    let a = MatRef::from_row_major_slice(a, rows, inner);
    let b = MatRef::from_row_major_slice(b, cols, inner);
    let c = MatMut::from_row_major_slice_mut(c, rows, cols);
    // The transpose of a matrix reference is another reference over the
    // same values, read the other way round: faer walks them as they lie
    // and nothing is copied.
    matmul(c, Accum::Replace, a, b.transpose(), 1.0, the_threads());
    Ok(())
}

/// Writes `a' b` into `c`, with `a` of exactly `inner` x `rows` values,
/// `b` of `inner` x `cols` and `c` of `rows` x `cols`, all row after row
/// and every dimension 1 at least.
///
/// # Errors
///
/// None, as for the products above.
#[expect(
    clippy::unnecessary_wraps,
    reason = "the two backends have the same signature, and the BLAS one fails when a dimension is larger than the i32 its routines take"
)]
pub(crate) fn product_with_the_first_turned(
    a: &[f64],
    rows: usize,
    inner: usize,
    b: &[f64],
    cols: usize,
    c: &mut [f64],
) -> Result<()> {
    let a = MatRef::from_row_major_slice(a, inner, rows);
    let b = MatRef::from_row_major_slice(b, inner, cols);
    let c = MatMut::from_row_major_slice_mut(c, rows, cols);
    // The transpose is another reference over the same values, as in the
    // product above: nothing is copied.
    matmul(c, Accum::Replace, a.transpose(), b, 1.0, the_threads());
    Ok(())
}

/// Writes `a' b'` into `c`, with `a` of exactly `inner` x `rows` values,
/// `b` of `cols` x `inner` and `c` of `rows` x `cols`, all row after row
/// and every dimension 1 at least.
///
/// # Errors
///
/// None, as for the products above.
#[expect(
    clippy::unnecessary_wraps,
    reason = "the two backends have the same signature, and the BLAS one fails when a dimension is larger than the i32 its routines take"
)]
pub(crate) fn product_with_both_turned(
    a: &[f64],
    rows: usize,
    inner: usize,
    b: &[f64],
    cols: usize,
    c: &mut [f64],
) -> Result<()> {
    let a = MatRef::from_row_major_slice(a, inner, rows);
    let b = MatRef::from_row_major_slice(b, cols, inner);
    let c = MatMut::from_row_major_slice_mut(c, rows, cols);
    // Both operands are read the other way round, and both transposes are
    // references over the same values: nothing is copied here either.
    matmul(
        c,
        Accum::Replace,
        a.transpose(),
        b.transpose(),
        1.0,
        the_threads(),
    );
    Ok(())
}

/// The Cholesky factorization of the symmetric positive definite `a`, of
/// exactly `n` x `n` values row after row with its lower half filled and
/// `n` 1 at least, which overwrites that lower half with the lower
/// triangular `l`.
///
/// faer is given the regularization it calls its default, which is the one
/// that refuses: told to regularize, it patches a pivot that is not
/// positive and factors on, and [`Error::Singular`] would never be raised.
///
/// # Errors
///
/// [`Error::Singular`] when `a` is not positive definite, with the row
/// faer stopped at, which it counts from 0 as this crate does.
pub(crate) fn cholesky_lower(a: &mut [f64], n: usize) -> Result<()> {
    let a = MatMut::from_row_major_slice_mut(a, n, n);
    // The scratch faer asks for here is n values, one column of the
    // matrix, which is 8 KB at the 1000 individuals of the spec.
    let mut scratch = MemBuffer::new(cholesky_in_place_scratch::<f64>(
        n,
        the_threads(),
        Default::default(),
    ));
    cholesky_in_place(
        a,
        LltRegularization::default(),
        the_threads(),
        MemStack::new(&mut scratch),
        Default::default(),
    )
    .map(|_| ())
    .map_err(|error| match error {
        LltError::NonPositivePivot { index } => Error::Singular {
            argument: "a",
            at: index,
        },
    })
}

/// The `x` of `a x = b` for the factorization `l` of exactly `n` x `n`
/// values row after row with its lower half filled, and `b` of exactly
/// `sides` x `n` values row after row, one row for each right hand side,
/// which comes back holding the solutions the same way. `n` and `sides`
/// are 1 at least.
///
/// faer takes the right hand sides as the columns of a matrix of `n` rows,
/// so the buffer of `b` is given to it as the `sides` x `n` matrix it is
/// and then turned the other way round. That transpose is another
/// reference over the same values, as the ones of the products above are:
/// nothing is copied and nothing is moved.
///
/// # Errors
///
/// None: faer refuses nothing that the checks of `lib.rs` let through, and
/// a factorization it was given is one it does not read for a pivot. The
/// signature is the one of the BLAS backend, which fails when a dimension
/// is larger than the `i32` its routines take.
#[expect(
    clippy::unnecessary_wraps,
    reason = "the two backends have the same signature, and the BLAS one fails when a dimension is larger than the i32 its routines take"
)]
pub(crate) fn solve_with_cholesky(l: &[f64], n: usize, b: &mut [f64], sides: usize) -> Result<()> {
    let l = MatRef::from_row_major_slice(l, n, n);
    let b = MatMut::from_row_major_slice_mut(b, sides, n).transpose_mut();
    // The scratch faer asks for here is nothing at all in faer 0.24.4: the
    // solve walks the two triangles in place and `solve_in_place_scratch`
    // gives an empty request, measured at 0 bytes for the n of 5 and the
    // 10000 right hand sides the association study solves at. So this
    // allocation is `MemBuffer::new` and not the `try_new` of the inverse,
    // whose scratch is n x n and which gives `Error::Memory` for it.
    let mut scratch = MemBuffer::new(solve_in_place_scratch::<f64>(n, sides, the_threads()));
    solve_in_place_with_conj(l, Conj::No, b, the_threads(), MemStack::new(&mut scratch));
    Ok(())
}

/// The lower half of the inverse of the `a` whose factorization `l` is,
/// with `l` of exactly `n` x `n` values row after row with its lower half
/// filled and `inverse` of exactly `n` x `n`, and `n` 1 at least. The
/// upper half of `inverse` is left as it was, faer writing the lower
/// triangle of its result alone.
///
/// faer inverts the triangle into a matrix of its own and multiplies it by
/// itself into `inverse`, so it asks for a scratch of `n` x `n` values, 8
/// MB at the 1000 individuals the spec measured it at and 800 MB at the
/// 10000 of `docs/objectives.md`. That is the one allocation of this
/// backend large enough that a machine can be without it and go on, so it
/// is asked for with `try_new` and not taken with `MemBuffer::new`, which
/// ends the process when it fails, and what a caller gets instead is
/// [`Error::Memory`], as it does for the workspaces the BLAS backend asks
/// for the eigendecomposition.
///
/// # Errors
///
/// [`Error::Memory`] when this machine has not the memory for that
/// scratch. faer refuses nothing else that the checks of `lib.rs` let
/// through, and a factorization it was given is one it does not read for a
/// pivot: a diagonal entry of 0 makes it divide by 0 and write infinities
/// and NaN where `dpotri` of the BLAS backend gives an `info`, which is
/// why `lib.rs` reads that diagonal before either backend runs.
pub(crate) fn invert_with_cholesky(l: &[f64], n: usize, inverse: &mut [f64]) -> Result<()> {
    let l = MatRef::from_row_major_slice(l, n, n);
    let inverse = MatMut::from_row_major_slice_mut(inverse, n, n);
    let request = inverse_scratch::<f64>(n, the_threads());
    let mut scratch = MemBuffer::try_new(request).map_err(|_| Error::Memory {
        what: "the scratch of the inverse of faer",
        values: the_values_of_a_request(request),
    })?;
    inverse_with_cholesky_of_faer(inverse, l, the_threads(), MemStack::new(&mut scratch));
    Ok(())
}

/// How many `f64` a request of faer holds: faer asks for a number of
/// bytes and [`Error::Memory`] carries a number of values, as its doc
/// comment says. A request that is not a whole number of them is as many
/// values as it takes.
fn the_values_of_a_request(request: StackReq) -> usize {
    request.size_bytes().div_ceil(size_of::<f64>())
}

/// The thin QR of `a`, of exactly `rows` x `cols` values row after row
/// with `rows` at least `cols` and `cols` 1 at least: `q` of exactly
/// `rows` x `cols` values and the upper triangular `r` of exactly `cols`
/// x `cols`, both written row after row and the lower half of `r` set
/// to 0.
///
/// faer is given the buffer as it lies and copies nothing of it, where
/// the BLAS backend writes the transpose of `a` into a buffer of its own
/// because its two routines are much slower on the wide matrix that the
/// buffer is in their view. faer builds the two matrices of its own,
/// `compute_thin_Q` the `q` and `thin_R` the `r`, and what it leaves
/// below the diagonal of that `r` is its own, so a 0 is written there.
///
/// Both run on the threads faer's own default takes, the global pool of
/// rayon natively and one thread in WebAssembly, as the
/// eigendecomposition below does: neither takes an argument for them.
///
/// # Errors
///
/// None: faer refuses nothing that the checks of `lib.rs` let through.
/// The signature is the one of the BLAS backend, which fails when a
/// dimension, or a workspace its routines ask for, is larger than the
/// `i32` they take.
#[expect(
    clippy::unnecessary_wraps,
    reason = "the two backends have the same signature, and the BLAS one fails when a dimension is larger than the i32 its routines take"
)]
pub(crate) fn thin_qr(
    a: &[f64],
    rows: usize,
    cols: usize,
    q: &mut [f64],
    r: &mut [f64],
) -> Result<()> {
    let factorization = MatRef::from_row_major_slice(a, rows, cols).qr();
    let orthogonal = factorization.compute_thin_Q();
    let triangular = factorization.thin_R();
    for (row, entries) in q.chunks_exact_mut(cols).enumerate() {
        for (column, entry) in entries.iter_mut().enumerate() {
            *entry = orthogonal[(row, column)];
        }
    }
    for (row, entries) in r.chunks_exact_mut(cols).enumerate() {
        for (column, entry) in entries.iter_mut().enumerate() {
            *entry = if column < row {
                0.0
            } else {
                triangular[(row, column)]
            };
        }
    }
    Ok(())
}

/// The eigendecomposition of the symmetric `g`, of exactly `n` x `n`
/// values row after row with its lower half filled and `n` 1 at least.
///
/// The eigenvalues come back from the smallest, which is the order faer
/// gives them in, and each eigenvector is a row of the buffer: `lib.rs`
/// turns both round together.
///
/// faer builds the eigenvectors in a matrix of its own, one eigenvector
/// per column, so this backend writes them into the buffer of `g` a row at
/// a time, while the BLAS one has `dsyevd` overwrite that buffer in place.
/// What a caller sees is the same.
///
/// # Errors
///
/// [`Error::NoConvergence`] when faer reached its limit of iterations.
pub(crate) fn eigh_lower(g: Vec<f64>, n: usize) -> Result<Eigen> {
    let decomposition = MatRef::from_row_major_slice(&g, n, n)
        .self_adjoint_eigen(Side::Lower)
        .map_err(|error| match error {
            EvdError::NoConvergence => Error::NoConvergence {
                routine: "faer",
                info: 0,
            },
        })?;
    // faer gives the eigenvalues from the smallest, as LAPACK does, and
    // the eigenvector of the j-th as the column j of a matrix of its own,
    // so column j goes into row j of the buffer and the order is left as
    // it is: `lib.rs` turns the values and the rows round together.
    let eigenvalues = decomposition.S();
    let eigenvalues = eigenvalues.column_vector();
    let eigenvectors = decomposition.U();
    let values: Vec<f64> = (0..n).map(|column| eigenvalues[column]).collect();
    let mut vectors = g;
    for (row, column) in vectors.chunks_exact_mut(n).zip(0..n) {
        for (entry, entry_of_the_column) in row.iter_mut().zip(0..n) {
            *entry = eigenvectors[(entry_of_the_column, column)];
        }
    }
    Ok(Eigen { values, vectors })
}

/// The `x` of `r x = b` for the upper triangular `r` of exactly `n` x `n`
/// values row after row with its upper half filled, and `b` of exactly
/// `sides` x `n` values row after row, one row for each right hand side,
/// which comes back holding the solutions the same way. `n` and `sides`
/// are 1 at least.
///
/// faer takes the right hand sides as the columns of a matrix of `n` rows,
/// so the buffer of `b` is given to it as the `sides` x `n` matrix it is
/// and then turned the other way round, as it is for the solve with the
/// Cholesky above: that transpose is another reference over the same
/// values and nothing is copied.
///
/// # Errors
///
/// None: faer refuses nothing that the checks of `lib.rs` let through. It
/// divides by a diagonal entry of 0 as it finds it and answers with an
/// infinity, where `dtrtrs` of the BLAS backend gives an `info`, which is
/// why `lib.rs` reads that diagonal before either backend runs. The
/// signature is the one of that backend, which fails when a dimension is
/// larger than the `i32` its routine takes.
#[expect(
    clippy::unnecessary_wraps,
    reason = "the two backends have the same signature, and the BLAS one fails when a dimension is larger than the i32 its routines take"
)]
pub(crate) fn solve_upper_triangular(
    r: &[f64],
    n: usize,
    b: &mut [f64],
    sides: usize,
) -> Result<()> {
    let r = MatRef::from_row_major_slice(r, n, n);
    let b = MatMut::from_row_major_slice_mut(b, sides, n).transpose_mut();
    solve_upper_triangular_in_place(r, b, the_threads());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{inverse_scratch, solve_in_place_scratch, the_threads, the_values_of_a_request};

    /// The size of the matrix the association study solves at, and how
    /// many right hand sides it gives it: five coefficients and one right
    /// hand side for each of 10000 individuals.
    const THE_SIZE_THE_GWAS_SOLVES_AT: (usize, usize) = (5, 10000);

    #[test]
    fn the_scratch_of_the_solve_asks_for_no_memory_at_the_size_the_gwas_solves_at() {
        // faer 0.24.4 walks the two triangles of the factorization in
        // place and asks for nothing, measured here and not read from its
        // documentation. That is why the buffer for it is taken with
        // `MemBuffer::new`, which ends the process when an allocation
        // fails, and not with the `try_new` that the inverse needs for its
        // scratch of n x n. A faer that asked for memory here would fail
        // this test, and then the solve would need `Error::Memory` too,
        // which is a point for the spec and for the owner.
        let (n, sides) = THE_SIZE_THE_GWAS_SOLVES_AT;
        let request = solve_in_place_scratch::<f64>(n, sides, the_threads());
        assert_eq!(request.size_bytes(), 0);
    }

    /// The n the spec measured faer's scratch for the inverse at, which is
    /// the 1000 individuals of "How the seven are verified".
    const THE_SIZE_THE_SPEC_MEASURED_THE_INVERSE_AT: usize = 1000;

    #[test]
    fn the_scratch_of_inverting_with_the_cholesky_is_the_n_by_n_the_spec_measured() {
        // 8000000 bytes at n = 1000, which is n x n values of the 8 bytes
        // of an `f64` and the number "The inverse of a factorized matrix"
        // of `docs/specs/linalg.md` records. This is the allocation of
        // this backend that is asked for with `try_new`, since at the
        // 10000 individuals of `docs/objectives.md` it is 800 MB, and what
        // `Error::Memory` carries for it is those values and not the
        // bytes.
        let n = THE_SIZE_THE_SPEC_MEASURED_THE_INVERSE_AT;
        let request = inverse_scratch::<f64>(n, the_threads());
        assert_eq!(request.size_bytes(), 8_000_000);
        assert_eq!(the_values_of_a_request(request), n * n);
    }
}
