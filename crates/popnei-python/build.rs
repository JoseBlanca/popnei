//! The link arguments an extension module needs, which pyo3 leaves to the
//! crate that builds one.
//!
//! On macOS a shared library whose symbols are resolved by whoever loads it
//! is linked with `-undefined dynamic_lookup`, and without it `cargo build`
//! ends in "Undefined symbols for architecture arm64" over every `Py_`
//! symbol. On `wasm32-unknown-emscripten` and on Linux the call adds
//! nothing under the rustc this repository uses. maturin and pyodide-build
//! pass the same arguments themselves; a second copy of a link argument
//! changes nothing.

fn main() {
    pyo3_build_config::add_extension_module_link_args();
}
