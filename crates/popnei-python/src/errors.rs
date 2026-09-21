//! The errors of the core crate on their way to a Python exception.
//!
//! `impl From<popnei::Error> for PyErr` cannot be written here, because
//! neither type belongs to this crate, so every function of the crate fails
//! with [`PyPopneiError`], which does belong to it, and pyo3 turns that one
//! into the exception. `?` on a call of the core crate works everywhere, and
//! a call site maps an error by hand only to add what the core does not
//! have, the path of the file, which is what
//! `.claude/skills/coding/pyo3.md` asks for.

use std::path::{Path, PathBuf};

use pyo3::exceptions::{PyOSError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;

/// What a function of this crate fails with.
pub(crate) enum PyPopneiError {
    /// Something the core crate refused: an argument it takes, or what it
    /// found in a file.
    Core(popnei::Error),
    /// Something the core refused while a file was being read, with the
    /// file, which the core does not carry and this crate knows: a reader
    /// is built over bytes, and only the call that opened the path knows
    /// which file they are.
    OfTheFile {
        /// What the core refused.
        error: popnei::Error,
        /// The file that was being read, as the caller gave it and not as
        /// text: the name of a file is bytes under Linux and macOS, and
        /// pyo3 gives a `PathBuf` back to Python as the standard library
        /// does, with the bytes that are not UTF-8 as the surrogates of
        /// `os.fsdecode`, so that `open(error.filename)` opens the file
        /// the caller asked for.
        path: PathBuf,
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
    /// The error of the core with the file it happened in: what a user is
    /// told then names that file, in `OSError.filename` where the exception
    /// is one of the file system, and a user who reads a directory of VCFs
    /// knows which one to look at.
    pub(crate) fn of_the_file(error: popnei::Error, path: &Path) -> PyPopneiError {
        PyPopneiError::OfTheFile {
            error,
            path: path.to_path_buf(),
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
            PyPopneiError::Core(error) => exception_of(error, None),
            PyPopneiError::OfTheFile { error, path } => exception_of(error, Some(path)),
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

/// The exception of one error of the core crate, by the convention the
/// owner gave on 21 September 2026: an `OSError` for a file that cannot be
/// read, that was cut short or that is corrupted; a `RuntimeError` for a
/// defect of popnei; and a `ValueError` for a wrong input of a function,
/// which a file whose content is not what a VCF holds is.
///
/// `path` is the file the error happened in, for the calls that read one.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "popnei::Error is non_exhaustive, so a match on it outside the core crate \
              has to have a wildcard arm; a case that a later module adds is a ValueError, \
              which is what every case that is neither of the file system nor a defect of \
              popnei is"
)]
fn exception_of(error: popnei::Error, path: Option<PathBuf>) -> PyErr {
    let message = error.to_string();
    match error {
        // The two that carry a cause of the file system. The path goes to
        // `filename` and not into the message: Python prints an `OSError`
        // with the file at its end, so a message that named it too would
        // say it twice.
        popnei::Error::FileNotOpened {
            path: of_the_core,
            source,
        } => os_error(
            source.raw_os_error(),
            format!("the file could not be opened: {source}"),
            path.or(Some(of_the_core)),
        ),
        popnei::Error::Io(source) => os_error(
            source.raw_os_error(),
            format!("the file could not be read: {source}"),
            path,
        ),
        // A file that was cut short and one that is corrupted are errors of
        // the file and not of what a user wrote, so they are an `OSError`
        // too, with no number: nothing of the system refused anything, and
        // what is wrong is in the bytes of the file.
        popnei::Error::VcfBgzipEndMissing | popnei::Error::VcfBgzipCorrupted { .. } => {
            os_error(None, message, path)
        }
        // The cases that say popnei has a defect: the three with which
        // `docs/specs/block.md` says that a reader has one, blocks of a
        // source that do not hold the same dataset, a block whose arrays
        // are not of its size and a block of no variants; and the parse of
        // a batch of lines that did not come back, which a panic inside it
        // leaves behind. Nothing a user asks for gives them, so a user who
        // gets one reports it instead of looking for what they typed wrong.
        popnei::Error::BlocksDoNotFitTogether { .. }
        | popnei::Error::BlockArrayOfAnotherSize { .. }
        | popnei::Error::ReaderGaveABlockOfNoVariants
        | popnei::Error::VcfParseNotFinished { .. } => {
            PyRuntimeError::new_err(of_the_file(message, path))
        }
        // The arguments a user writes: how many variants a block holds,
        // and how many alleles a genotype of the file has, which the reader
        // is given when the file is opened because it needs it to read the
        // first genotype. What is wrong with them is wrong whatever file is
        // read, so they name no file although they are refused while one is
        // being opened.
        popnei::Error::BlockOfNoVariants
        | popnei::Error::BlockTooLarge { .. }
        | popnei::Error::VcfPloidyOutOfRange { .. } => PyValueError::new_err(message),
        _ => PyValueError::new_err(of_the_file(message, path)),
    }
}

/// The message with the file it happened in before it, which is what a user
/// who reads a directory of VCFs needs in order to know which one to look
/// at. The core has the line, the column and the value, and not the file: a
/// reader is built over bytes, and the call that opened the path is where
/// the two meet.
fn of_the_file(message: String, path: Option<PathBuf>) -> String {
    match path {
        // A message is text, and the name of a file is bytes under Linux
        // and macOS. A byte that is not UTF-8 is shown here as the
        // replacement character, which is what a reader of the message
        // needs; the name itself travels whole in `filename`, where the
        // exception is an `OSError`, and a caller who has to open the file
        // again has it from the call they made.
        Some(path) => format!("{path}: {message}", path = path.to_string_lossy()),
        None => message,
    }
}

/// The error of the file system as Python raises it itself: built with the
/// number the system gave, it is the `FileNotFoundError`, the
/// `IsADirectoryError` or the `PermissionError` of that number, and it
/// carries the file in `filename`, where the standard library puts it and
/// where it prints it, after the message.
///
/// A cause that no number came with, a gzip stream that ends in the middle
/// or a file that bgzip wrote and that was cut, is an `OSError` whose
/// `errno` is `None` and whose `filename` is the file all the same: it is
/// the file a user needs, and which of the two ways the read failed is not
/// theirs to tell apart.
///
/// The file goes in as an `OsString`, which pyo3 gives to Python as the
/// text the standard library would, the bytes that are not UTF-8 as the
/// surrogates of `os.fsdecode`, so that `error.filename` is the path the
/// caller wrote and `open(error.filename)` opens their file. A `PathBuf`
/// would arrive as a `pathlib.Path`, which no `OSError` of Python carries
/// and which is not what the caller gave when they gave a `str`.
fn os_error(number: Option<i32>, message: String, path: Option<PathBuf>) -> PyErr {
    match path {
        Some(path) => PyOSError::new_err((number, message, path.into_os_string())),
        // A source that is not a file, which Python has none of yet.
        None => PyOSError::new_err(message),
    }
}
