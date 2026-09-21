//! What a Python user reaches through `write_vars`: the variants of a
//! source written into a vars file, the arrow file popnei keeps them in.
//!
//! The core writes into anything that takes bytes and has no path of its
//! own, so the file is made here: this crate opens the path, refuses one
//! that a file is already at, and takes away what was written when the core
//! gives an error, so that a call that failed leaves nothing behind and can
//! be made again once what was wrong is fixed.

use std::fs::File;
use std::io::{BufWriter, ErrorKind};
use std::path::{Path, PathBuf};

use pyo3::prelude::*;

use popnei::io::vcf::VcfReader;

use crate::errors::PyPopneiError;
use crate::vcf::{VcfSource, count_of};

// Every variant of `source` into a vars file at `path`, one batch of
// `num_vars_per_block` variants after another, and `None` for the size
// popnei chooses for the individuals of the source. A `///` comment here
// would become the `__doc__` of `popnei._core.write_vars`, and what a
// Python user reads belongs to the package, which is the API.
#[pyfunction]
#[pyo3(signature = (source, path, num_vars_per_block))]
pub(crate) fn write_vars(
    py: Python<'_>,
    source: &VcfSource,
    path: PathBuf,
    num_vars_per_block: Option<&Bound<'_, PyAny>>,
) -> Result<(), PyPopneiError> {
    let num_vars_per_block = num_vars_per_block
        .map(|asked_for| count_of("num_vars_per_block", asked_for))
        .transpose()?;
    // A Ctrl-C that was pending when this was called is raised here, before
    // a file is made: what a user stopped leaves no file at the path.
    py.check_signals()?;
    let vcf = source.path().to_path_buf();
    let options = source.options();
    let file = file_at(&path)?;
    // The whole source is read inside this one call, which is seconds for a
    // VCF of hundreds of megabytes, so the interpreter is released for all
    // of it. A Ctrl-C that arrives meanwhile is raised when the call is
    // over and not between two blocks, as it is in `Blocks::__next__`: the
    // loop over the blocks is the core's, which `docs/specs/io_vars.md` has
    // this crate call instead of writing that loop again.
    let written = py.detach(|| -> Result<(), popnei::Error> {
        let reader = VcfReader::from_path(&vcf, options)?;
        let sink = popnei::io::vars::write_vars(reader, BufWriter::new(file), num_vars_per_block)?;
        // What the buffer still holds goes to the file here, where its
        // error is read; a buffer that is dropped writes it and loses it.
        let file = sink
            .into_inner()
            .map_err(|failure| not_written(failure.into_error()))?;
        // The bytes reach the disc here and not when the file is closed,
        // where nothing reads what the close said: a file system that only
        // then says that it is full would leave a file that is not whole
        // after a call that returned and said nothing.
        file.sync_all().map_err(not_written)?;
        Ok(())
    });
    if let Err(error) = written {
        let taken = take_away(&path);
        // The file the error names is the one it is about: a wrong line of
        // the VCF names the VCF, and a disc that filled up names the file
        // that was being written, which the core keeps apart from an error
        // of a source that could not be read.
        let refusal = of_the_file_it_is_about(error, &vcf, &path);
        return Err(with_what_was_left(refusal, &path, taken));
    }
    // A Ctrl-C that arrived while the file was being written is still
    // pending: the interpreter was released and no bytecode ran to raise
    // it. A user who stopped the call is told so and is left with no file,
    // as every call that fails leaves none, so that the call they make
    // again finds the path free.
    if let Err(interrupted) = py.check_signals() {
        let taken = take_away(&path);
        return Err(with_what_was_left(interrupted.into(), &path, taken));
    }
    Ok(())
}

/// `error` with the file it is about: the vars file when the write is what
/// failed, and the source that was being read otherwise.
fn of_the_file_it_is_about(error: popnei::Error, vcf: &Path, vars: &Path) -> PyPopneiError {
    let of_the_write = matches!(error, popnei::Error::VarsFileNotWritten { .. });
    PyPopneiError::of_the_file(error, if of_the_write { vars } else { vcf })
}

/// What the file system said while the vars file was being written, as the
/// error of a vars file that could not be written, which is the case the
/// core keeps for it. The number the system gave travels with it, since
/// that number is what the exception of Python is built with.
fn not_written(failure: std::io::Error) -> popnei::Error {
    popnei::Error::VarsFileNotWritten {
        problem: failure.to_string(),
        source: Some(failure),
    }
}

/// The file to write the variants into, made at `path`.
///
/// # Errors
///
/// When a file is already at the path, which is a wrong argument and not an
/// error of the file system: `write_vars` writes no file over another one.
/// And when the file system refuses to make the file, a directory that is
/// not there or one that cannot be written in.
fn file_at(path: &Path) -> Result<File, PyPopneiError> {
    // `create_new` asks and makes in one call: a path that is looked at
    // first and written afterwards is a file of somebody else in between.
    File::create_new(path).map_err(|error| {
        if error.kind() == ErrorKind::AlreadyExists {
            PyPopneiError::PathTaken {
                path: path.to_path_buf(),
            }
        } else {
            // Everything else the file system refused, a directory that is
            // not there or one that cannot be written in, is the error of
            // a file that could not be opened, with the number the system
            // gave, which is what makes Python raise the exception of that
            // cause with the path in `filename`.
            PyPopneiError::Core(popnei::Error::FileNotOpened {
                path: path.to_path_buf(),
                source: error,
            })
        }
    })
}

/// The file that the call was writing, taken away from the path.
///
/// # Errors
///
/// When the file is still there: a directory whose permissions changed
/// while the file was being written is one way, and it is the caller that
/// tells the user, since it is their error that is being reported.
fn take_away(path: &Path) -> Result<(), std::io::Error> {
    if let Err(problem) = std::fs::remove_file(path) {
        // A file that is not there any more is a file that was taken away:
        // something else removed it, and the path is free for the call the
        // user makes again, which is what they are told about.
        if problem.kind() != ErrorKind::NotFound {
            return Err(problem);
        }
    }
    Ok(())
}

/// `error`, and with it the file that is still at `path` when it could not
/// be taken away, which the user's next call at that path would refuse.
fn with_what_was_left(
    error: PyPopneiError,
    path: &Path,
    taken: Result<(), std::io::Error>,
) -> PyPopneiError {
    match taken {
        Ok(()) => error,
        Err(problem) => PyPopneiError::LeftBehind {
            error: Box::new(error),
            path: path.to_path_buf(),
            problem: problem.to_string(),
        },
    }
}
