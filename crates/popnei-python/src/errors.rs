//! The errors of the core crate on their way to a Python exception.
//!
//! `impl From<popnei::Error> for PyErr` cannot be written here, because
//! neither type belongs to this crate, so the functions of the crate fail
//! with [`PyPopneiError`], which does belong to it, and pyo3 turns that one
//! into the exception. `?` on a call of the core crate works everywhere and
//! no call site has a `map_err`, which is what
//! `.claude/skills/coding/pyo3.md` asks for.

use std::path::Path;

use pyo3::exceptions::{PyOSError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;

/// What a function of this crate fails with.
pub(crate) enum PyPopneiError {
    /// Something the core crate refused: an argument it takes, or what it
    /// found in a file.
    Core(popnei::Error),
    /// The bytes of a file could not be read, with the file they were read
    /// from, which the core does not carry there and this crate knows.
    NotRead {
        /// Why the file system refused the read.
        source: std::io::Error,
        /// The file that was being read.
        path: String,
    },
    /// An argument that says how many of something there are, the ploidy or
    /// the variants of a block, and holds a number that counts nothing: a
    /// negative one, or one above what this machine counts, which in wasm
    /// is 4295 million.
    Count {
        /// The name of the argument, as a Python user writes it.
        name: &'static str,
        /// What was given for it.
        value: i64,
    },
    /// Something that cannot happen unless this crate has a defect: a lock
    /// a panic left broken, or a chromosome whose number is not in the
    /// table of the reader that gave it.
    Broken(String),
}

impl PyPopneiError {
    /// The error of the core, with `path` in it when it is the file system
    /// that refused: a directory where a VCF was asked for, a file that was
    /// taken away while it was read, a disc that answers no more.
    #[expect(
        clippy::wildcard_enum_match_arm,
        reason = "popnei::Error is non_exhaustive, so a match on it outside the core \
                  crate has to have a wildcard arm; every case but the one of a read \
                  that failed goes on as it is"
    )]
    pub(crate) fn of_the_file(error: popnei::Error, path: &Path) -> PyPopneiError {
        match error {
            popnei::Error::Io(source) => PyPopneiError::NotRead {
                source,
                path: path.to_string_lossy().into_owned(),
            },
            other => PyPopneiError::Core(other),
        }
    }
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
            PyPopneiError::NotRead { source, path } => os_error(
                source.raw_os_error(),
                format!("the file {path} could not be read: {source}"),
                path,
            ),
            PyPopneiError::Count { name, value } => PyValueError::new_err(format!(
                "`{name}` is {value}, and it says how many of something there are: 0 or \
                 more, and at most {largest}",
                largest = usize::MAX
            )),
            PyPopneiError::Broken(message) => PyRuntimeError::new_err(message),
        }
    }
}

/// The exception of one error of the core crate, which a pyNei user
/// recognises: `ValueError` for an argument that is wrong or a file whose
/// content popnei cannot read, `OSError` for the file system.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "popnei::Error is non_exhaustive, so a match on it outside the core crate \
              has to have a wildcard arm; a case that a later module adds is a ValueError, \
              which is what every case that is not of the file system is"
)]
fn exception_of(error: popnei::Error) -> PyErr {
    let message = error.to_string();
    match error {
        popnei::Error::FileNotOpened { path, source } => os_error(
            source.raw_os_error(),
            message,
            path.to_string_lossy().into_owned(),
        ),
        // The path of a read that failed is put in by `of_the_file`, which
        // the calls that have it use. This one is left for a source that is
        // not a file, which Python has none of yet.
        popnei::Error::Io(_) => PyOSError::new_err(message),
        _ => PyValueError::new_err(message),
    }
}

/// The error of the file system as Python raises it itself: built with the
/// number the system gave, it is the `FileNotFoundError`, the
/// `IsADirectoryError` or the `PermissionError` of that number, and it
/// carries the file in `filename`, where the standard library puts it.
///
/// A cause that no number came with, which is an error of Rust's own, is an
/// `OSError` with the message alone.
fn os_error(number: Option<i32>, message: String, path: String) -> PyErr {
    match number {
        Some(number) => PyOSError::new_err((number, message, path)),
        None => PyOSError::new_err(message),
    }
}
