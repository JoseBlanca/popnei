//! The errors of the core crate on their way to a Python exception.
//!
//! `impl From<popnei::Error> for PyErr` cannot be written here, because
//! neither type belongs to this crate, so every function of the crate fails
//! with [`PyPopneiError`], which does belong to it, and pyo3 turns that one
//! into the exception. `?` on a call of the core crate works everywhere, and
//! a call site maps an error by hand only to add what the core does not
//! have, the path of the file, which is what
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
        /// What was given for it, as Python prints it: an integer of Python
        /// is of any size, so the number that was refused does not always
        /// fit in one of Rust.
        value: String,
    },
    /// Something that cannot happen unless this crate has a defect: a lock
    /// a panic left broken, or a chromosome whose number is not in the
    /// table of the reader that gave it.
    Broken(String),
    /// An exception the interpreter itself raised, on its way back to it as
    /// it is: the `KeyboardInterrupt` of a Ctrl-C that `check_signals`
    /// found between two blocks, and what building a tuple of the names of
    /// the chromosomes of a block raised. This crate reads them and does
    /// not choose them, so it carries them back untouched.
    Python(PyErr),
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

impl From<PyErr> for PyPopneiError {
    /// So that `?` on a call of pyo3, `check_signals` or the building of a
    /// tuple, works in a function that fails with this type, which every
    /// function of this crate does.
    fn from(error: PyErr) -> PyPopneiError {
        PyPopneiError::Python(error)
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
            // No largest number is named here. What the largest is depends
            // on the argument, 255 for a ploidy, and the core says it of
            // each: a bound of this crate beside it would give a user two
            // limits for one argument, and the one they read first would be
            // the one that is not theirs.
            PyPopneiError::Count { name, value } => PyValueError::new_err(format!(
                "`{name}` is {value}, and it says how many of something there are: a \
                 whole number of 1 or more that this machine can count"
            )),
            PyPopneiError::Broken(message) => PyRuntimeError::new_err(message),
            PyPopneiError::Python(error) => error,
        }
    }
}

/// The exception of one error of the core crate, which a pyNei user
/// recognises: `ValueError` for an argument that is wrong or a file whose
/// content popnei cannot read, `OSError` for the file system, and
/// `RuntimeError` for the two that say a reader of the core has a defect.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "popnei::Error is non_exhaustive, so a match on it outside the core crate \
              has to have a wildcard arm; a case that a later module adds is a ValueError, \
              which is what every case that is neither of the file system nor a defect of \
              a reader is"
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
        // The three cases with which `docs/specs/block.md` says that a
        // reader has a defect: blocks of one source that do not hold the
        // same dataset, a block whose arrays are not of its size, and a
        // block of no variants. Nothing a user asks for gives them, so they
        // are the `RuntimeError` of `PyPopneiError::Broken` and not the
        // `ValueError` of the rest, and a user who gets one reports it
        // instead of looking for what they typed wrong.
        popnei::Error::BlocksDoNotFitTogether { .. }
        | popnei::Error::BlockArrayOfAnotherSize { .. }
        | popnei::Error::ReaderGaveABlockOfNoVariants => PyRuntimeError::new_err(message),
        _ => PyValueError::new_err(message),
    }
}

/// The error of the file system as Python raises it itself: built with the
/// number the system gave, it is the `FileNotFoundError`, the
/// `IsADirectoryError` or the `PermissionError` of that number, and it
/// carries the file in `filename`, where the standard library puts it.
///
/// A cause that no number came with, which is an error of Rust's own, a
/// gzip stream that ends in the middle among them, is an `OSError` whose
/// `errno` is `None` and whose `filename` is the file all the same: it is
/// the file a user needs, and which of the two ways the read failed is not
/// theirs to tell apart.
fn os_error(number: Option<i32>, message: String, path: String) -> PyErr {
    PyOSError::new_err((number, message, path))
}
