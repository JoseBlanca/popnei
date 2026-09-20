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
//! every individual at it. The `variant` module has it, with the trait that
//! anything giving variants implements, and `io` the VCF reader, which
//! fills it from a VCF. The `block` module holds a run of consecutive
//! variants as arrays, which is what the calculations with matrices consume
//! and the only way genotypes leave the core; the modules that calculate
//! over them are being written, and `docs/architecture.md` has their order.

#![forbid(unsafe_code)]

pub mod block;
pub mod error;
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
