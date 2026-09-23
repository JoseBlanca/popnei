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

use pyo3::create_exception;
use pyo3::exceptions::{PyOSError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;

create_exception!(
    popnei._core,
    TraitsWithNoVariance,
    PyValueError,
    "The traits of a table that is to be standardized and that have no \
     variance. `args[0]` is what the core says, which names them by their \
     position, and `args[1]` their positions among the traits, from 0.\n\n\
     `popnei.do_pca` catches it and raises the `ValueError` its user reads, \
     whose message names those traits as the frame names them. The core has \
     the positions and not the names, and this class is how they reach the \
     layer that has the frame. It derives from `ValueError`, so a user who \
     catches that one catches this one as well."
);

create_exception!(
    popnei._core,
    TraitOutOfRange,
    PyValueError,
    "A trait whose mean or whose standard deviation is not a number a \
     principal component analysis can use, its values being too large or \
     too small for the arithmetic of an f64. `args[0]` is what the core \
     says, `args[1]` the position of the trait among the traits, from 0, \
     and `args[2]` which of the three it is: `mean_not_finite`, \
     `deviation_not_finite` or `deviation_of_zero`.\n\n\
     `popnei.do_pca` catches it and raises the `ValueError` its user reads, \
     with the name the frame gives that trait, as it does for the traits \
     with no variance. It derives from `ValueError` as well."
);

create_exception!(
    popnei._core,
    KinshipPairWithNoVariantCalled,
    PyValueError,
    "Two individuals of a kinship that have no variant called in both of \
     them, so that the sum of their pair would be divided by no variant at \
     all. `args[0]` is what the core says, with the file that was read, \
     which names the two by their position; `args[1]` and `args[2]` are \
     those positions among the individuals of the kinship, from 0, and are \
     one position twice when an individual has no called genotype at all; \
     and `args[3]` and `args[4]` are how many of the variants that were \
     used are called in each of the two.\n\n\
     `popnei.calc_kinship` catches it and raises the `ValueError` its user \
     reads, whose message names the two individuals. The core has their \
     positions and not their names, and this class is how they reach the \
     layer that has the names, as it is for the two traits above. It \
     derives from `ValueError`, so a user who catches that one catches this \
     one as well."
);

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
    /// negative one, one above what this machine counts, which in wasm is
    /// 4295 million, or a 0 where the argument takes 1 at least and the
    /// core has no message of its own for it.
    Count {
        /// The name of the argument, as a Python user writes it.
        name: &'static str,
        /// The fewest of them the argument takes, which is 1 for every one
        /// of them but the components the weights of a principal component
        /// analysis are given for, where 0 is no weights.
        smallest: usize,
        /// What was given for it, as Python prints it: an integer of Python
        /// is of any size, so the number that was refused does not always
        /// fit in one of Rust.
        value: String,
    },
    /// A threshold of a filter that is not a number from 0 to 1, under the
    /// name of the argument a user wrote it in: the core refuses it and
    /// names the filter by its kind, `maf`, and what a user has to look at
    /// is the call they wrote, `filter_by_maf(1.5)`.
    Threshold {
        /// The name of the argument, as a Python user writes it,
        /// `max_allowed_maf`.
        name: &'static str,
        /// What was given for it, which is NaN, below 0 or above 1, as
        /// Python prints it: a whole number of Python is of any size, so a
        /// threshold that was refused does not always fit in one of Rust.
        value: String,
    },
    /// An array that does not lie in memory row after row, which the core
    /// reads as a slice of values and cannot take, under the name of the
    /// argument a user wrote it in. The Python package makes every array C
    /// contiguous before the call, so a user reaches it only through
    /// `popnei._core`, and the message says what makes one.
    ArrayNotContiguous {
        /// The name of the argument, as a Python user writes it.
        name: &'static str,
    },
    /// An array of two dimensions that is not square, given to a call that
    /// takes a matrix of the individuals by the individuals, under the name
    /// of the argument a user wrote it in. The package builds that array
    /// from the frame of a `Kinship`, which is square by the checks of that
    /// class, so a user reaches it only through `popnei._core`.
    MatrixNotSquare {
        /// The name of the argument, as a Python user writes it.
        name: &'static str,
        /// How many rows the array has.
        num_rows: usize,
        /// How many columns it has.
        num_columns: usize,
    },
    /// A path that a file is already at, given to a call that writes one.
    /// This crate refuses it before the core is called and writes nothing,
    /// which is what `docs/specs/io_vars.md` asks of `write_vars`, as in
    /// pyNei.
    PathTaken {
        /// The path the caller gave.
        path: PathBuf,
    },
    /// What went wrong in a call that was writing a file, with the file it
    /// could not take away afterwards, a directory whose permissions
    /// changed under it among the causes. What went wrong is what the user
    /// reads, and the file that was left is a note on it: their next call
    /// finds that path taken and would say only that.
    LeftBehind {
        /// What went wrong, which is the exception the user gets.
        error: Box<PyPopneiError>,
        /// The file that is still at the path.
        path: PathBuf,
        /// Why it could not be taken away, as the system said it.
        problem: String,
    },
    /// Something that cannot happen unless this crate has a defect: a lock
    /// a panic left broken, or a chromosome whose number is not in the
    /// table of the reader that gave it.
    Broken {
        /// What went wrong, for whoever reports it.
        message: String,
        /// The file that was being read, where there is one.
        path: Option<PathBuf>,
    },
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

    /// A defect of this crate that was found while `path` was being read,
    /// which the message names as every error of a file does.
    pub(crate) fn broken_of_the_file(message: String, path: &Path) -> PyPopneiError {
        PyPopneiError::Broken {
            message,
            path: Some(path.to_path_buf()),
        }
    }
}

/// Raises the Ctrl-C that arrived while the interpreter was released, which
/// is still pending: no bytecode ran to raise it.
///
/// It is raised before numpy is called, because the first array of a
/// process imports the C API of numpy, that import fails with the exception
/// that is pending, and the numpy crate panics when it does: a user who
/// asked for a Ctrl-C would get a `PanicException`, which no `except` of
/// theirs catches and which ends the session.
///
/// Every call that releases the interpreter for a whole pass over a source
/// and then builds an array of what it found calls this between the two.
///
/// # Errors
///
/// The `KeyboardInterrupt` of that Ctrl-C, on its way back as it is.
pub(crate) fn raise_a_ctrl_c_before_numpy_is_called(py: Python<'_>) -> Result<(), PyPopneiError> {
    py.check_signals()?;
    Ok(())
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
            PyPopneiError::Count {
                name,
                smallest,
                value,
            } => PyValueError::new_err(format!(
                "`{name}` is {value}, and it says how many of something there are: a \
                 whole number of {smallest} or more that this machine can count"
            )),
            // The threshold of a filter, which is the number a user wrote
            // in the call that adds it: the message names the argument, and
            // the rule it broke is the core's, which refuses the same
            // thresholds when a pass builds its filters.
            PyPopneiError::Threshold { name, value } => PyValueError::new_err(format!(
                "`{name}` is {value}, and a threshold is a number from 0 to 1, both \
                 included: the number of the variant it is compared with is one count of \
                 the variant divided by another"
            )),
            // The values of an array that does not lie row after row are
            // not a slice, and the core takes a slice: what a user does
            // about it is to make an array that is contiguous, which is
            // what the message says.
            PyPopneiError::ArrayNotContiguous { name } => PyValueError::new_err(format!(
                "`{name}` does not lie in memory row after row, and popnei reads the \
                 values of an array as they lie: `numpy.ascontiguousarray({name})` \
                 gives one that does"
            )),
            PyPopneiError::MatrixNotSquare {
                name,
                num_rows,
                num_columns,
            } => PyValueError::new_err(format!(
                "`{name}` is {num_rows} by {num_columns}, and the matrix of a kinship                  is a square one of the individuals by the individuals"
            )),
            // A file that is already at the path is a wrong argument of the
            // call and not an error of the file system, so it is a
            // `ValueError`, whose message starts with the path as that of
            // every error of a file does.
            PyPopneiError::PathTaken { path } => PyValueError::new_err(of_the_file(
                "a file is already there, and popnei writes no file over another one: \
                 write to another path, or take that file away"
                    .to_owned(),
                Some(path),
            )),
            PyPopneiError::LeftBehind {
                error,
                path,
                problem,
            } => left_behind(PyErr::from(*error), &path, &problem),
            PyPopneiError::Broken { message, path } => {
                PyRuntimeError::new_err(of_the_file(message, path))
            }
            PyPopneiError::Python(error) => error,
        }
    }
}

/// `raised` with a note that says that the file the call was writing is
/// still at `path`, because it could not be taken away.
///
/// A note is text that Python keeps in `__notes__` and prints under the
/// message of the exception, which is where what a user has to do about a
/// second thing goes: what went wrong stays the exception they see, with
/// its kind and its message.
fn left_behind(raised: PyErr, path: &Path, problem: &str) -> PyErr {
    let note = format!(
        "the file that was being written is still at {path}, because it could not be \
         taken away: {problem}. The call made again at that path is refused while it \
         is there",
        path = path.to_string_lossy()
    );
    Python::attach(|py| {
        // A note that could not be added does not take the place of what
        // went wrong, which is what the user asked about and what this
        // returns either way.
        let _ = raised.value(py).call_method1("add_note", (note,));
    });
    raised
}

/// The exception of one error of the core crate, by the convention the
/// owner gave on 21 September 2026: an `OSError` for a file that cannot be
/// read, that was cut short or that is corrupted; a `RuntimeError` for a
/// defect of popnei; and a `ValueError` for a wrong input of a function,
/// which a file whose content is not what a VCF holds is. The memory a
/// calculation asked this machine for and was not given is a `ValueError`
/// as well, for want of a fourth exception that says it.
///
/// `path` is the file the error happened in, for the calls that read one.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "popnei::Error is non_exhaustive, so a match on it outside the core crate \
              has to have a wildcard arm; a case that a later module adds is a ValueError \
              with the file it was read from before its message, which is what a wrong \
              input found in a file is, and the cases of a module that are a defect of \
              popnei or that name no file are the ones listed by name above"
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
            format!("the file could not be opened: {}", what_went_wrong(&source)),
            path.or(Some(of_the_core)),
        ),
        popnei::Error::Io(source) => os_error(
            source.raw_os_error(),
            format!("the file could not be read: {}", what_went_wrong(&source)),
            path,
        ),
        // The vars file that a call was writing and that the file system
        // or arrow-rs refused, which is an error of that file and not of
        // the source the call was reading: `path` is the file being
        // written wherever this case travels, and the message of the core
        // says already that it could not be written. The number is the
        // system's when the file system is what refused, so that Python
        // raises the exception of that number, and there is none when
        // arrow-rs refused what it was handed.
        popnei::Error::VarsFileNotWritten { ref source, .. } => {
            let number = source.as_ref().and_then(std::io::Error::raw_os_error);
            os_error(number, without_the_number(message, number), path)
        }
        // A file that was cut short and one that is corrupted are errors of
        // the file and not of what a user wrote, so they are an `OSError`
        // too, with no number: nothing of the system refused anything, and
        // what is wrong is in the bytes of the file. A vars file that was
        // damaged after it was written is one of the two when its bytes no
        // longer decode: it ends before what it says it holds, or a batch
        // of it cannot be decoded. Damage that does decode, into content
        // the format does not allow, is a `ValueError` below, as an allele
        // of `gts` below the missing one is: the exception follows what is
        // wrong with the content, and the reader cannot tell a file a disc
        // changed from one another program wrote badly.
        popnei::Error::VcfBgzipEndMissing
        | popnei::Error::VcfBgzipCorrupted { .. }
        | popnei::Error::VarsFileCutShort { .. }
        | popnei::Error::VarsBatchNotRead { .. } => os_error(None, message, path),
        // The cases that say popnei has a defect: the three with which
        // `docs/specs/block.md` says that a reader has one, blocks of a
        // source that do not hold the same dataset, a block whose arrays
        // are not of its size and a block of no variants; the number of
        // values a filter of popnei gave `Block::retain_vars`, which is one
        // for each variant of the block it was given; and the parse of a
        // batch of lines that did not come back, which a panic inside it
        // leaves behind. Nothing a user asks for gives them, so a user who
        // gets one reports it instead of looking for what they typed wrong.
        // The three of the vars file writer are of the same kind: a block
        // that does not hold the individuals of the file, one whose columns
        // are not those of the first block written, and a chromosome number
        // that the table given with the block has no name for. `write_vars`
        // gives the writer the individuals, the fields and the table of one
        // reader, so a user reaches them only through a reader with a
        // defect. The three of the counts of one variant are of that kind
        // too, which `docs/specs/variant.md` says in "The Rust interface":
        // the counts have no function in Python, so the genotypes they
        // refuse, the ploidy they were given and a variant of more alleles
        // than a count of them holds are a reader's and not a user's.
        // The three of `Block::retain_individuals` are of that kind as
        // well, the indices of the individuals a filter of individuals
        // keeps: they come from `resolve_individuals`, which refuses the
        // name behind an index at or beyond the individuals of the block,
        // behind one that is there twice, and a call that names no
        // individual at all, so a user who gets one of them has read what a
        // reader with a defect built. A block that holds the genotypes of
        // no individual, which a calculation over the variants refuses, is
        // one more of a reader's: the VCF reader refuses a header with no
        // individual and a ploidy of 0, which are the two ways a block
        // comes out like that.
        popnei::Error::IndividualToKeepNotInTheBlock { .. }
        | popnei::Error::IndividualToKeepTwice { .. }
        | popnei::Error::NoIndividualToKeep
        | popnei::Error::GtsNotWholeGenotypes { .. }
        | popnei::Error::MoreAllelesThanACountHolds { .. }
        | popnei::Error::AlleleBelowTheMissingOne { .. }
        | popnei::Error::IndividualBeyondTheVariant { .. }
        | popnei::Error::BlocksDoNotFitTogether { .. }
        | popnei::Error::BlockArrayOfAnotherSize { .. }
        | popnei::Error::ReaderGaveABlockOfNoVariants
        | popnei::Error::BlockWithNoGenotypeOfAVariant { .. }
        | popnei::Error::KeepOfAnotherSize { .. }
        | popnei::Error::VcfParseNotFinished { .. }
        | popnei::Error::VarsBlockDoesNotFit { .. }
        | popnei::Error::VarsBlockColumns { .. }
        | popnei::Error::VarsChromNameMissing { .. }
        // The two of the principal component analysis that no argument of
        // `do_pca` gives, which is what "Errors and the cases pyNei asserts"
        // of `docs/specs/pca.md` says of them: a buffer that does not hold
        // the rows times the traits it was said to hold, which this crate
        // takes from the array itself, and an operation of the linear
        // algebra that did not run, which is left with a table whose
        // products are not finite and a machine with too little memory for
        // the workspace of the eigendecomposition.
        | popnei::Error::PcaTableOfAnotherSize { .. }
        | popnei::Error::PcaLinalg { .. }
        // The three of the principal components of the variants that no
        // argument of `do_pca_from_variants` gives: a second pass over the
        // variants that was not made, which this crate opens a reader for
        // whenever the weights are asked for; a second pass that read other
        // variants than the first, which is what a source that changed
        // between the two gives; and a weight that had no column to go in,
        // which the second pass counts against the variants of the first as
        // it goes, so nothing a user writes reaches it.
        | popnei::Error::PcaSecondPassMissing { .. }
        | popnei::Error::PcaSecondPassDiffers { .. }
        | popnei::Error::PcaWeightOutOfPlace { .. }
        // The one of the kinship that no argument of `calc_kinship` gives:
        // a product of the linear algebra that did not run, which is the
        // standardized dosages of a block with themselves or the genotypes
        // that were called with themselves. Every size of both products is
        // checked before they are asked for, the individuals at the entry of
        // the calculation and the variants of the block as it is read, so
        // what is left is a defect of popnei or a backend that refused the
        // work, and a user reports it.
        | popnei::Error::KinshipLinalg { .. }
        // The four of the r² of two sets of variants that no argument of
        // `calc_rogers_huff_r2_matrix` gives, for the same reason as the
        // two of the principal component analysis above: a range of
        // variants that is not in the dosages, which the tiles of the
        // products and the window of the filter by linkage disequilibrium
        // ask for; two sets of dosages built over different individuals of
        // the block, which one call of this crate builds both of; a buffer
        // for the r² that does not hold one value for each pair, which
        // this crate holds and a user never sees; and a product of the
        // linear algebra that did not run, which is left with a result of
        // more values than the routines of BLAS and LAPACK count in, a
        // size the cap of `calc_r2_matrix` refuses before a user reaches
        // it.
        | popnei::Error::LdRowsNotInTheDosages { .. }
        | popnei::Error::LdDosagesOfOtherIndividuals { .. }
        | popnei::Error::LdR2OfAnotherSize { .. }
        | popnei::Error::LdLinalg { .. }
        // The one of the distances between populations that no argument of
        // `calc_pop_dists` gives: the sums of a resampling group that do
        // not hold one place for each pair of the populations, which every
        // variant of a block is added into. The places are made from the
        // populations the pass counts over, so a user who gets it reports
        // it instead of looking for what they typed wrong.
        | popnei::Error::PopDistSumsOfAnotherSize { .. }
        // The plain filter of a threshold built for the criterion of the
        // filter by linkage disequilibrium, which it does not answer:
        // whether a variant passes that one turns on the variants kept
        // before it and not on the variant alone. No call of a user reaches
        // it. This crate builds the filter that does answer each criterion
        // when a step is added, and `chain_of` builds the readers of a
        // pass, so a user who gets this one has found a defect of a caller
        // of the core crate and reports it instead of looking at what they
        // wrote. It is named here because the arm below would make it the
        // `ValueError` of an argument a user wrote, which there is none of.
        | popnei::Error::VarFilterOfTheLdCriterion
        // The three of the association study that no argument of
        // `calc_gwas` gives: the phenotype, the design and the positions of
        // the tested individuals not holding the same individuals, which
        // this crate is given built from one list of them; a model that
        // answered for another number of variants than the block it was
        // given holds; and an operation of the linear algebra that did not
        // run, which is the rank of the design, the thin QR a model is
        // fitted with, a solve against it or a product of a block with
        // something the null model holds. Every size and every value of the
        // three is checked before the linear algebra is called, so what is
        // left is a defect of popnei or a backend that refused the work.
        // "The Rust interface" of `docs/specs/gwas.md` has them, each as
        // the `RuntimeError` it is here.
        | popnei::Error::GwasInputOfAnotherSize { .. }
        | popnei::Error::GwasAnswersOfAnotherSize { .. }
        | popnei::Error::GwasLinalg { .. } => {
            PyRuntimeError::new_err(of_the_file(message, path))
        }
        // The two errors of a trait that the layer holding the frame names:
        // this crate has the positions of those traits and not their names,
        // and `popnei.do_pca` catches these exceptions and raises the
        // `ValueError` a user reads. Each of them carries what the core
        // says as well, so that a caller of `popnei._core` reads a message
        // and not a list of numbers. Both derive from `ValueError`, so
        // nothing of a user's changes when one reaches them.
        popnei::Error::PcaTraitsWithNoVariance { positions, .. } => {
            TraitsWithNoVariance::new_err((message, positions))
        }
        popnei::Error::PcaTraitOutOfRange { position, problem } => {
            TraitOutOfRange::new_err((message, position, name_of(problem)))
        }
        // The pair of individuals of a kinship that the layer holding their
        // names names: the core has where each of the two is among the
        // individuals of the kinship and not what they are called, and
        // `popnei.calc_kinship` raises the `ValueError` a user reads, whose
        // message names them. It carries what the core says as well, with
        // the file the variants were read from, so that a caller of
        // `popnei._core` reads a message and not four numbers.
        popnei::Error::KinshipPairWithNoVariantCalled {
            one,
            other,
            num_vars_of_one,
            num_vars_of_other,
        } => KinshipPairWithNoVariantCalled::new_err((
            of_the_file(message, path),
            one,
            other,
            num_vars_of_one,
            num_vars_of_other,
        )),
        // The arguments a user writes: how many variants a block holds,
        // and how many alleles a genotype of the file has, which the reader
        // is given when the file is opened because it needs it to read the
        // first genotype. The three of `docs/specs/filters.md` are of the
        // same kind: the threshold of a filter that is not a number from 0
        // to 1, a second filter of a kind the variants are filtered by
        // already, and a window of the filter by linkage disequilibrium
        // that is no base pairs wide, all three of which a user gets at the
        // call that adds the filter. The four of the filter of individuals
        // are of it too: a name that is
        // not an individual of the variants, a name that is there twice, a
        // call that names none, and a second filter of individuals, all of
        // them what a user wrote in the call that adds the step. The four
        // of `pops`, the populations a statistic is calculated for, are the
        // same kind of thing in the argument of the call that calculates
        // it: a name that is not an individual of the variants, a name
        // twice in one population, a population that names no individual,
        // and `pops` with no population at all. The three of the
        // histogram of a statistic are of the same kind, in `hist_kwargs`:
        // a histogram of no bin, a range that does not run from a number up
        // to a larger one, and a range of bins of equal ratio that starts at
        // 0 or below. So is the ploidy or the exponent of a statistic of one
        // variant that is 0 or above 255, which the pass of the statistics
        // names its `ploidy` argument before it comes here. The threshold
        // below which a variant counts as polymorphic in a population is one
        // more: `poly_threshold` is a number from 0 to 1, which is where a
        // major allele frequency lies. What is wrong with them is wrong
        // whatever file is read, so they name no file although some of them
        // are refused while one is being opened.
        popnei::Error::BlockOfNoVariants
        | popnei::Error::BlockTooLarge { .. }
        | popnei::Error::VcfPloidyOutOfRange { .. }
        | popnei::Error::VarFilterThresholdOutOfRange { .. }
        | popnei::Error::VarFilterOfAKindThatIsSet { .. }
        | popnei::Error::LdFilterMaxDistTooSmall { .. }
        | popnei::Error::IndividualNotInTheSource { .. }
        | popnei::Error::IndividualNamedTwice { .. }
        | popnei::Error::NoIndividualNamed
        | popnei::Error::FilterOfIndividualsThatIsSet { .. }
        | popnei::Error::IndividualOfAPopNotInThePass { .. }
        | popnei::Error::IndividualNamedTwiceInAPop { .. }
        | popnei::Error::PopWithNoIndividual { .. }
        | popnei::Error::NoPop
        | popnei::Error::HistWithNoBin
        | popnei::Error::HistRangeNotGoingUp { .. }
        | popnei::Error::HistLogRangeNotAboveZero { .. }
        | popnei::Error::StatPloidyOutOfRange { .. }
        | popnei::Error::PolyThresholdOutOfRange { .. }
        // The four of the table of a principal component analysis that a
        // user writes: a value of it that is not finite, a table to be
        // standardized and not centered, one of fewer than 2 rows or of no
        // traits, and one in which no trait has variance once it is
        // centered, which has no direction to give. The table comes from
        // the user and not from a file, so they name none.
        | popnei::Error::PcaValueNotFinite { .. }
        | popnei::Error::PcaStandardizeWithoutCentering
        | popnei::Error::PcaTableTooSmall { .. }
        | popnei::Error::PcaNoTraitWithVariance
        // The five of the r² of a set of variants that are wrong whatever
        // file is read: an index that is not an individual of the dataset
        // and one given twice, which are the individuals of a population
        // as a user writes them; and the three sizes the calculation
        // cannot be done at, dosages of more values than the linear
        // algebra counts in, a variant of more alleles than the sums come
        // out of exactly, and a matrix this machine has not the memory
        // for. What a user does about each of the last three is calculate
        // over fewer variants or over fewer individuals, whichever file
        // they read.
        | popnei::Error::LdIndividualNotInTheDataset { .. }
        | popnei::Error::LdIndividualAskedForTwice { .. }
        | popnei::Error::LdDosagesTooLarge { .. }
        | popnei::Error::LdTooManyAllelesInAVariant { .. }
        | popnei::Error::LdNoMemory { .. }
        // The `max_num_vars` of the matrix of every pair that is more
        // variants than this machine counts the pairs of, which is the one
        // number a user writes at that call and nothing of any file: it is
        // looked at before the pass, so the same number is refused whatever
        // the source holds.
        | popnei::Error::LdMaxNumVarsTooLarge { .. }
        // The nine of the association study that are of what a user wrote
        // and are wrong whatever file is read: a phenotype or a covariate
        // that is not a finite number, which the package lets through as a
        // value that came out of the user's own arithmetic as an infinity;
        // a binomial trait whose phenotype is not 0 or 1 or is one value
        // for everybody; covariates that are not independent, a copy of one
        // another or a constant; the score test asked of a linear model and
        // the Wald test of a logistic mixed one, which are the two pairs no
        // model has; the GRAMMAR-Gamma approximation asked for by a study
        // with no kinship; and a study whose trait and kinship ask for one
        // of the three models that are not written. "The Rust interface" of
        // `docs/specs/gwas.md` has them, each as the `ValueError` it is
        // here.
        | popnei::Error::GwasPhenotypeNotFinite { .. }
        | popnei::Error::GwasPhenotypeNotBinomial { .. }
        | popnei::Error::GwasPhenotypeOfOneValue { .. }
        | popnei::Error::GwasDesignValueNotFinite { .. }
        | popnei::Error::GwasKinshipValueNotFinite { .. }
        | popnei::Error::GwasCovariatesCollinear { .. }
        | popnei::Error::GwasScoreTestOfALinearModel
        | popnei::Error::GwasWaldTestOfALogisticMixedModel
        | popnei::Error::GwasGrammarGammaWithoutAKinship
        | popnei::Error::GwasModelNotBuilt { .. }
        // The name of a trait and the name of a test that are of neither
        // of the two, which a user writes in `trait` and in `test`.
        | popnei::Error::GwasTraitOfAnUnknownName { .. }
        | popnei::Error::GwasTestOfAnUnknownName { .. } => PyValueError::new_err(message),
        // The six that the dataset a user gave is wrong for: four of the
        // principal components of the variants and two of the pass over a
        // row that those components and the kinship share. Of the
        // components: no variants, which the steps of a `Variants` can
        // leave; no variant with variance, which one individual gives; a
        // source of no individual, which is nobody to place on the axes and
        // which no source of popnei is, since one that names no individual
        // is refused when it is opened; and a dataset of a size the
        // analysis cannot count in, which "Errors and the cases pyNei
        // asserts" of `docs/specs/pca.md` lists. Of the row: a variant of
        // more than two different alleles among its called genotypes with
        // `transform_to_biallelic` false, and a ploidy the dosages could
        // not be written one to a byte at. Each names the file the variants
        // were read from, as every error of a file does.
        popnei::Error::PcaNoVariants
        | popnei::Error::PcaNoVariantWithVariance
        | popnei::Error::VariantWithMoreThanTwoAlleles { .. }
        | popnei::Error::VariantPloidyTooLarge { .. }
        | popnei::Error::PcaNoIndividual
        | popnei::Error::PcaVariantsTooLarge { .. }
        // The three of the kinship that the dataset a user gave is wrong
        // for, which are of that same kind: no variant with variance among
        // the individuals it was asked for, which one individual gives and
        // which pyNei raises for as well; a source with no individual,
        // which is nobody to give
        // a kinship of and which no source of popnei is, since one that
        // names no individual is refused when it is opened; and a dataset
        // of a size the calculation cannot count in, more individuals than
        // the matrix of the linear algebra holds or more variants than this
        // machine counts. `docs/specs/kinship.md` has them in "The Rust
        // interface", each as the `ValueError` it is here, and the fourth,
        // a pair with no variant called in both, is the arm above, which
        // carries the two to the layer that has their names.
        | popnei::Error::KinshipNoVariantWithVariance
        | popnei::Error::KinshipNoIndividual
        | popnei::Error::KinshipVariantsTooLarge { .. }
        // The pass that gave more variants than `max_num_vars`, which is
        // of that same kind, a dataset larger than the calculation takes:
        // the cap a user wrote and the variants the file holds decide it
        // together, so a user who runs over a directory of files needs to
        // know which of them the cap was too low for, and the message says
        // both numbers and the memory the matrix of those variants would
        // have needed. The cap that no matrix could be held under, above,
        // names no file, because that one is wrong before any file is
        // opened.
        | popnei::Error::LdTooManyVars { .. }
        // The nine of the distances between populations, which "The Rust
        // interface" of `docs/specs/dists.md` lists. Four are of what a
        // user wrote and name no file, since what is wrong with them is
        // wrong whatever file is read: a measure under a name that is of
        // none of the seven, resampling groups of 0 base pairs, fewer than
        // two populations, and populations that make more pairs than this
        // machine counts. The other five are of the variants that were
        // read: fewer resampling groups than a standard error is built
        // from, a variant whose position goes back and one of a chromosome
        // that the variants before it had left, the six sums of every pair
        // and group that the machine has not the memory for, and a source
        // whose genotypes hold more alleles than popnei reads. Which of the
        // two a case is, `with_its_file` of `pop_dists.rs` decides: it is
        // the call that knows whether a file was being read.
        | popnei::Error::PopDistMeasureOfAnUnknownName { .. }
        | popnei::Error::JackknifeGroupOfNoBasePairs
        | popnei::Error::PopDistsOfFewerThanTwoPops { .. }
        | popnei::Error::PopDistsOfTooManyPops { .. }
        | popnei::Error::TooFewJackknifeGroups { .. }
        | popnei::Error::JackknifeGroupsVariantGoesBack { .. }
        | popnei::Error::JackknifeGroupsChromComesBack { .. }
        | popnei::Error::PopDistSumsTooLarge { .. }
        | popnei::Error::PopDistsPloidyOutOfRange { .. }
        // The five of the association study that the dataset a user gave is
        // wrong for: fewer tested individuals than the columns of the
        // design plus two, which leaves nothing to measure the uncertainty
        // of a variant from; the three the core makes of the positions of
        // those individuals, one that is not in the source, one that is
        // there twice and an order that is not the source's, which the
        // package cannot reach, since it builds those positions by walking
        // the individuals the pass gives; and a source of more variants
        // than this machine counts them in.
        | popnei::Error::GwasTooFewIndividuals { .. }
        | popnei::Error::GwasIndividualNotInTheDataset { .. }
        | popnei::Error::GwasIndividualTestedTwice { .. }
        | popnei::Error::GwasIndividualsOutOfOrder { .. }
        | popnei::Error::GwasVariantsTooLarge => {
            PyValueError::new_err(of_the_file(message, path))
        }
        // Everything else is a wrong input of a function, which a file
        // whose content is not what the format holds is, and it names the
        // file it was found in: the wrong data lines and headers of the VCF
        // reader, and the thirteen cases of the vars file that "The Rust
        // interface" of `docs/specs/io_vars.md` lists as a `ValueError`,
        // among them a `qual` that is a value and is not finite, an allele
        // of `gts` below the missing one, and a file whose genotypes hold
        // no allele, which `open_vars` gives for a `popnei` key that names
        // no individual. The variant that the filter by linkage
        // disequilibrium was given and whose position does not rise within
        // its chromosome is one of them: what is wrong is the order the
        // variants come in the file, which that filter is the one reader of
        // popnei to refuse. The block with more text
        // or more alleles in one column than a column of a batch takes is
        // one no call from Python reaches: 2147483647 bytes of text or
        // alleles in one block is more memory than a machine gives. A pass
        // that gave no variant is here too: which file was read is what a
        // user needs in order to see whether it is the file that holds
        // none or the steps that kept none of what it holds, and the
        // message says which of the two it was.
        _ => PyValueError::new_err(of_the_file(message, path)),
    }
}

/// The name Python reads one of the three scales of a trait under, as the
/// kind of a filter travels under its name: the layer that has the frame
/// writes the message, and it chooses the words by this.
fn name_of(problem: popnei::pca::TraitScale) -> &'static str {
    match problem {
        popnei::pca::TraitScale::MeanNotFinite => "mean_not_finite",
        popnei::pca::TraitScale::DeviationNotFinite => "deviation_not_finite",
        popnei::pca::TraitScale::DeviationOfZero => "deviation_of_zero",
    }
}

/// What an error of the input says, without the number of the system that
/// Rust writes at the end of it, `No such file or directory (os error 2)`.
/// Python prints an `OSError` with that number before the message,
/// `[Errno 2]`, and a user reads it once.
fn what_went_wrong(source: &std::io::Error) -> String {
    without_the_number(source.to_string(), source.raw_os_error())
}

/// `said` without the ` (os error 2)` that Rust writes at the end of what
/// an error of the file system says, when `number` is that number.
fn without_the_number(said: String, number: Option<i32>) -> String {
    let Some(number) = number else {
        return said;
    };
    match said.strip_suffix(&format!(" (os error {number})")) {
        Some(without_it) => without_it.to_owned(),
        None => said,
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
