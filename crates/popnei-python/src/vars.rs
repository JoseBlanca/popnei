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
        sink.into_inner()
            .map_err(|error| popnei::Error::Io(error.into_error()))?;
        Ok(())
    });
    if let Err(error) = written {
        take_away(&path);
        // The file the error names is the VCF: what a user fixes is the
        // line of it that popnei could not read. An error of the file that
        // was being written, a disk with no room left among them, names the
        // VCF too, since the core says that a write failed and not which of
        // the two files it was reading or writing when it did.
        return Err(PyPopneiError::of_the_file(error, &vcf));
    }
    // A Ctrl-C that arrived while the file was being written is still
    // pending: the interpreter was released and no bytecode ran to raise
    // it. A user who stopped the call is told so and is left with no file,
    // as every call that fails leaves none, so that the call they make
    // again finds the path free.
    if let Err(interrupted) = py.check_signals() {
        take_away(&path);
        return Err(interrupted.into());
    }
    Ok(())
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
fn take_away(path: &Path) {
    // What is reported is the error the caller asked about, the line of
    // their VCF that popnei could not read, and a file that could not be
    // removed does not take its place. A path that keeps a file of a call
    // that failed is refused by the next call, which says that a file is
    // already there.
    let _ = std::fs::remove_file(path);
}
