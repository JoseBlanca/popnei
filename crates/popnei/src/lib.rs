//! popnei: population genetics over the variants of a dataset, the
//! successor of the Python library pyNei.
//!
//! This is the core crate, where every calculation lives. It is pure Rust:
//! no Python and no JavaScript reach it, and it builds for the two targets
//! of WebAssembly that popnei ships, `wasm32-unknown-unknown` for a web
//! application and `wasm32-unknown-emscripten` for the wheel of pyodide.
//! Python and TypeScript call it through a binding crate each, as section 8
//! of `docs/architecture.md` lays out.
//!
//! A variant is one site of the genome with its alleles and the genotype of
//! every individual at it, and the variants flow in blocks: a run of
//! consecutive variants held as arrays, with the genotypes of all of them
//! in one. The `block` module has the block, the trait that everything
//! giving blocks implements and the reader that puts blocks back to a size;
//! `io` the VCF reader, which parses the lines of a file into the rows of a
//! block; `variant` what the other modules say about one variant: which
//! fields a consumer wants, the table of the chromosome names, the view of
//! one variant of a block, and the row helpers over it, which count its
//! alleles and its genotypes and turn it into one standardized dosage per
//! individual; `filters` the variants that a user keeps
//! by a threshold, with the counts of what each filter was given and kept;
//! and `stats` the populations a statistic is calculated for, each a named
//! set of individuals, with the pass over the variants that gives, for each
//! population, the mean and the histogram of five statistics of a variant:
//! the observed heterozygosity, the major allele frequency, the expected
//! heterozygosity, plain and unbiased, and how many variants vary. `pca`
//! gives the principal components of a table of numbers, individuals by
//! traits, and of the variants of a reader; `dists` holds the Kosman
//! distance of every pair of individuals over the variants of a reader,
//! counted from the genotypes of each block as sets of bits; `ld` reads
//! the genotypes as dosages, how many alleles of a genotype are not the
//! major allele of its variant, which is what r², how much the genotype of
//! one variant says about the genotype of another, is worked out from;
//! `pop_dists` the seven measures of how far apart two populations are,
//! which one pass over the variants gives from the counts of each
//! population at each of them; `kinship` gives how much more of their
//! genome every pair of individuals shares than two drawn at random from
//! the same panel would, which is the matrix a mixed model of an
//! association study takes as the covariance of its random effect; and
//! `gwas` tests every variant against a trait of the individuals, giving
//! the effect of each variant on the trait, the uncertainty of that effect
//! and its p-value, of which the linear model is written so far. The
//! modules that follow them are being written, and `docs/architecture.md`
//! has their order.
//!
//! The linear algebra those modules need, the products of matrices and
//! the eigendecomposition, is not a module here but a crate beside this
//! one, `popnei-linalg`, because the calls it makes to BLAS and LAPACK
//! are `unsafe` and this crate forbids that. `docs/specs/linalg.md` says
//! what it gives and which backend runs where.

#![forbid(unsafe_code)]

pub mod block;
pub mod dists;
pub mod error;
pub mod filters;
pub mod gwas;
pub mod io;
pub mod kinship;
pub mod ld;
pub mod pca;
pub mod pop_dists;
pub mod stats;
pub mod variant;

pub use error::{Error, Result};

/// The version of this crate, `major.minor.patch`, as its manifest gives
/// it.
///
/// The Python and the TypeScript packages publish it as their own version,
/// so that a user who reports a result names the code that gave it.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::version;

    /// The literal is the version of the workspace manifest, and it changes
    /// with it.
    #[test]
    fn version_is_the_one_of_the_manifest() {
        assert_eq!(version(), "0.1.0");
    }
}
