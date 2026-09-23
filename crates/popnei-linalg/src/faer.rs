//! The faer backend: the six operations on faer, a linear algebra
//! library written in Rust, which runs where there is no BLAS to link and
//! natively when the crate is built with `--no-default-features`.
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

use faer::linalg::evd::EvdError;
use faer::linalg::matmul::matmul;
use faer::linalg::matmul::triangular::{BlockStructure, matmul as triangular_matmul};
use faer::{Accum, MatMut, MatRef, Par, Side};

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
pub(crate) fn product_by_transpose(
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
pub(crate) fn product_of_the_transpose(
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
pub(crate) fn product_of_both_transposes(
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
