//! The errors of the core crate on their way to a Python exception.
//!
//! `impl From<popnei::Error> for PyErr` cannot be written here, because
//! neither type belongs to this crate, so the functions of the crate fail
//! with [`PyPopneiError`], which does belong to it, and pyo3 turns that one
//! into the exception. `?` on a call of the core crate works everywhere and
//! no call site has a `map_err`, which is what
//! `.claude/skills/coding/pyo3.md` asks for.

use pyo3::exceptions::{PyOSError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;

/// What a function of this crate fails with.
pub(crate) enum PyPopneiError {
    /// Something the core crate refused: an argument it takes, or what it
    /// found in a file.
    Core(popnei::Error),
    /// An argument this crate refuses on its own, before the core sees it,
    /// with the message a Python user reads.
    Argument(String),
    /// Something that cannot happen unless this crate has a defect: a lock
    /// a panic left broken, or a chromosome whose number is not in the
    /// table of the reader that gave it.
    Broken(String),
}

impl From<popnei::Error> for PyPopneiError {
    fn from(error: popnei::Error) -> PyPopneiError {
        PyPopneiError::Core(error)
    }
}

impl From<PyPopneiError> for PyErr {
    fn from(error: PyPopneiError) -> PyErr {
        match error {
            PyPopneiError::Core(error) => exception_of(error),
            PyPopneiError::Argument(message) => PyValueError::new_err(message),
            PyPopneiError::Broken(message) => PyRuntimeError::new_err(message),
        }
    }
}

/// The exception of one error of the core crate, which a pyNei user
/// recognises: `ValueError` for an argument that is wrong or a file whose
/// content popnei cannot read, `OSError` for the file system.
///
/// An `OSError` built with the number of the error, the message and the
/// path is the `FileNotFoundError` or the `PermissionError` of that number,
/// which Python itself raises for the same causes, and it carries the path
/// in `filename`.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "popnei::Error is non_exhaustive, so a match on it outside the core crate \
              has to have a wildcard arm; a case that a later module adds is a ValueError, \
              which is what every case that is not of the file system is"
)]
fn exception_of(error: popnei::Error) -> PyErr {
    let message = error.to_string();
    match error {
        popnei::Error::FileNotOpened { path, source } => {
            let path = path.to_string_lossy().into_owned();
            match source.raw_os_error() {
                Some(number) => PyOSError::new_err((number, message, path)),
                None => PyOSError::new_err(message),
            }
        }
        popnei::Error::Io(_) => PyOSError::new_err(message),
        _ => PyValueError::new_err(message),
    }
}
