//! The Python binding crate: it translates between Python and the core
//! crate `popnei` and holds no calculation of its own, as section 8 of
//! `docs/architecture.md` and `.claude/skills/coding/pyo3.md` say.
//!
//! What it builds is the private module `popnei._core`, which only the
//! Python package `python/popnei` imports. The package is what a user
//! sees: it carries the signatures of pyNei, the defaults and the
//! docstrings.

use pyo3::prelude::*;

/// The version of the core crate, `major.minor.patch`.
///
/// The Python package publishes it as its own `__version__`, so that a
/// user who reports a result names the code that gave it.
#[pyfunction]
fn version() -> &'static str {
    popnei::version()
}

#[pymodule]
mod _core {
    #[pymodule_export]
    use super::version;
}
