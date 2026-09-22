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
//! fields a consumer wants, the table of the chromosome names and the view
//! of one variant of a block; and `filters` the variants that a user keeps
//! by a threshold, with the counts of what each filter was given and kept.
//! The modules that calculate over blocks are being written: `dists`, the
//! distances between individuals, holds so far the genotypes of a block as
//! the sets of bits that the counts of a pair are read from.
//! `docs/architecture.md` has the order of the ones that follow.

#![forbid(unsafe_code)]

pub mod block;
pub mod dists;
pub mod error;
pub mod filters;
pub mod io;
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
