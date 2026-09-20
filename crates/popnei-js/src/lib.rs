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
//! the `extern "C"` entry point that WebAssembly calls. That generated
//! code is `unsafe`, so each of them silences the lint that denies unsafe
//! code in the workspace, with the reason on the function itself.

use wasm_bindgen::prelude::wasm_bindgen;

/// The version of the core crate, `major.minor.patch`.
///
/// The TypeScript package publishes it as its own `version`, so that a
/// user who reports a result names the code that gave it.
#[wasm_bindgen]
#[must_use]
pub fn version() -> String {
    popnei::version().to_owned()
}
