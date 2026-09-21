//! What a Python user reaches through `open_vars` and `write_vars`: a vars
//! file, the arrow file popnei keeps its variants in, read as a source of
//! variants and written from one.
//!
//! [`VarsSource`] is a vars file that was opened. It reads the schema and
//! the footer of the file when it is built, so a file that is not a vars
//! file fails at `open_vars`, and the names of the individuals and the
//! ploidy come from the `popnei` key of that schema. Every pass over it
//! opens the file again and goes through the `Blocks` of `source.rs`, the
//! one a VCF goes through.
//!
//! The core writes into anything that takes bytes and has no path of its
//! own, so the file `write_vars` makes is made here: this crate opens the
//! path, refuses one that a file is already at, and takes away what was
//! written when the core gives an error, so that a call that failed leaves
//! nothing behind and can be made again once what was wrong is fixed.

use std::fs::File;
use std::io::{BufWriter, ErrorKind};
use std::path::{Path, PathBuf};

use pyo3::prelude::*;

use popnei::block::BlockReader;
use popnei::io::vars::VarsReader;

use crate::errors::PyPopneiError;
use crate::source::{Blocks, OpenSource, blocks_of, count_of, source_of};

// A vars file that was opened: its path, and the individuals and the ploidy
// its schema named. A `///` here would become the `__doc__` of the class,
// and what a Python user reads belongs to the package, which is the API.
#[pyclass(frozen, module = "popnei._core")]
pub(crate) struct VarsSource {
    path: PathBuf,
    individuals: Vec<String>,
    ploidy: usize,
}

#[pymethods]
impl VarsSource {
    // The names of the individuals, in the order of their genotypes in the
    // rows of a block.
    fn individuals(&self) -> Vec<String> {
        self.individuals.clone()
    }

    // How many alleles the genotype of one individual holds.
    fn ploidy(&self) -> usize {
        self.ploidy
    }

    // One pass over the file: it is opened again, and its blocks hold
    // `fields` besides the genotypes, `num_vars_per_block` variants each.
    #[pyo3(signature = (fields, num_vars_per_block))]
    fn blocks(
        &self,
        py: Python<'_>,
        fields: Vec<String>,
        num_vars_per_block: Option<&Bound<'_, PyAny>>,
    ) -> Result<Blocks, PyPopneiError> {
        blocks_of(py, self, fields, num_vars_per_block)
    }
}

impl OpenSource for VarsSource {
    fn path(&self) -> &Path {
        &self.path
    }

    /// The size the caller asks for is not passed on: the reader gives each
    /// batch of the file as a block, at the size the file was written with,
    /// and the `Reblock` that every pass ends with cuts them where the
    /// caller wants them.
    fn reader(
        &self,
        _num_vars_per_block: Option<usize>,
    ) -> Result<Box<dyn BlockReader>, popnei::Error> {
        Ok(Box::new(VarsReader::from_path(&self.path)?))
    }
}

// The vars file at `path`. It reads the schema and the footer, so the
// individuals, the ploidy and the batches are known when it returns and a
// file that is not a vars file fails here. A `///` comment would become the
// `__doc__` of `popnei._core.open_vars`, and what a Python user reads
// belongs to the package, which is the API.
#[pyfunction]
pub(crate) fn open_vars(py: Python<'_>, path: PathBuf) -> Result<VarsSource, PyPopneiError> {
    let (individuals, ploidy) = py
        .detach(|| -> Result<_, popnei::Error> {
            // The schema and the footer are read when the reader is built
            // and no batch is, so a file whose blocks would need more
            // memory than this machine gives is opened all the same.
            let reader = VarsReader::from_path(&path)?;
            let metadata = reader.metadata();
            Ok((metadata.individuals.clone(), metadata.ploidy))
        })
        .map_err(|error| PyPopneiError::of_the_file(error, &path))?;
    Ok(VarsSource {
        path,
        individuals,
        ploidy,
    })
}

// Every variant of `source` into a vars file at `path`, one batch of
// `num_vars_per_block` variants after another, and `None` for the size
// popnei chooses for the individuals of the source. `source` is a VCF that
// `open_vcf` opened or a vars file that `open_vars` did, which is what a
// `Variants` of the package holds. A `///` comment here would become the
// `__doc__` of `popnei._core.write_vars`, and what a Python user reads
// belongs to the package, which is the API.
#[pyfunction]
#[pyo3(signature = (source, path, num_vars_per_block))]
pub(crate) fn write_vars(
    py: Python<'_>,
    source: &Bound<'_, PyAny>,
    path: PathBuf,
    num_vars_per_block: Option<&Bound<'_, PyAny>>,
) -> Result<(), PyPopneiError> {
    let source = source_of(source)?;
    let num_vars_per_block = num_vars_per_block
        .map(|asked_for| count_of("num_vars_per_block", asked_for))
        .transpose()?;
    // A Ctrl-C that was pending when this was called is raised here, before
    // a file is made: what a user stopped leaves no file at the path.
    py.check_signals()?;
    let read = source.path().to_path_buf();
    let file = file_at(&path)?;
    // The whole source is read inside this one call, which is seconds for a
    // VCF of hundreds of megabytes, so the interpreter is released for all
    // of it. A Ctrl-C that arrives meanwhile is raised when the call is
    // over and not between two blocks, as it is in `Blocks::__next__`: the
    // loop over the blocks is the core's, which `docs/specs/io_vars.md` has
    // this crate call instead of writing that loop again.
    let written = py.detach(|| -> Result<(), popnei::Error> {
        // The source is opened at the size of its own blocks: the core puts
        // a `reblock` of `num_vars_per_block` over whatever it is given, so
        // the batches of the file hold that many variants whichever source
        // they came from.
        let reader = source.reader(None)?;
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
        // the VCF that was being read names that VCF, and a disc that
        // filled up names the file that was being written, which the core
        // keeps apart from an error of a source that could not be read.
        let refusal = of_the_file_it_is_about(error, &read, &path);
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

/// `error` with the file it is about: the vars file that was being written
/// when the write is what failed, and the file that was being read
/// otherwise.
fn of_the_file_it_is_about(error: popnei::Error, read: &Path, written: &Path) -> PyPopneiError {
    let of_the_write = matches!(error, popnei::Error::VarsFileNotWritten { .. });
    PyPopneiError::of_the_file(error, if of_the_write { written } else { read })
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

/// The number the system gives for a directory where a file was asked for,
/// `EISDIR`, which is 21 on macOS, on Linux and in emscripten, the systems
/// popnei runs on. `File::create_new` on a directory says instead that
/// something is already at the path, which is true and which would tell a
/// user to take away the directory they meant to write into.
const A_DIRECTORY_IS_THERE: i32 = 21;

/// The file to write the variants into, made at `path`.
///
/// # Errors
///
/// When a file is already at the path, which is a wrong argument and not an
/// error of the file system: `write_vars` writes no file over another one.
/// And when the file system refuses to make the file, a directory at the
/// path, a directory that is not there or one that cannot be written in.
fn file_at(path: &Path) -> Result<File, PyPopneiError> {
    // `create_new` asks and makes in one call: a path that is looked at
    // first and written afterwards is a file of somebody else in between.
    File::create_new(path).map_err(|error| {
        if error.kind() == ErrorKind::AlreadyExists {
            if path.is_dir() {
                // The path is asked about after the call and not before it,
                // so nothing is done on what the answer says: a directory
                // that became a file meanwhile is told as a file that is
                // already there, which it is.
                return PyPopneiError::Core(popnei::Error::FileNotOpened {
                    path: path.to_path_buf(),
                    source: std::io::Error::from_raw_os_error(A_DIRECTORY_IS_THERE),
                });
            }
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
