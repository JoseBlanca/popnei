//! The JavaScript binding crate: it translates between JavaScript and the
//! core crate `popnei` and holds no calculation of its own, as sections 8
//! and 11 of `docs/architecture.md` say.
//!
//! It is compiled for `wasm32-unknown-unknown`, WebAssembly with no
//! operating system under it, and the `wasm-bindgen` command line reads
//! the wasm file it gives and writes the JavaScript that calls into it,
//! with the TypeScript declarations of every exported function. The
//! TypeScript package `js/popnei` is what an application sees: it carries
//! the names of the Python API in camelCase and the result objects.
//!
//! Every function here is marked with `#[wasm_bindgen]`, which generates
//! the entry point that the generated JavaScript calls and the TypeScript
//! declaration of the function. What it generates passes the lints of the
//! workspace as it is, for the wasm target and for the native one, the
//! `unsafe_code` that is denied among them, so no function here silences
//! one: an `#[expect(unsafe_code)]` added here is itself an error, "this
//! lint expectation is unfulfilled".

use wasm_bindgen::prelude::wasm_bindgen;

pub mod dists;
pub mod errors;
pub mod gwas;
pub mod kinship;
pub mod ld;
pub mod pca;
pub mod pop_dists;
pub mod source;
pub mod stats;
pub mod steps;
pub mod vars;
pub mod vcf;

/// The version of the core crate, `major.minor.patch`.
///
/// The TypeScript package publishes it as its own `version`, so that a
/// user who reports a result names the code that gave it.
#[wasm_bindgen]
#[must_use]
pub fn version() -> String {
    popnei::version().to_owned()
}
