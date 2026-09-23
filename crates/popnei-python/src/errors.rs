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
use pyo3::exceptions::{PyMemoryError, PyOSError, PyRuntimeError, PyValueError};
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
    /// Memory this crate asked the machine for and was not given, with what
    /// a user reads of it: what could not be held, how large it is and what
    /// they do about it.
    ///
    /// The one case today is the copy of the matrix of r² that numpy is
    /// given, which `ld.rs` asks for with `try_reserve_exact`. The core
    /// refuses the memory of its own matrices the same way, and both are a
    /// `MemoryError` in Python.
    NoMemory {
        /// What a user reads, which the call site writes because it is the
        /// one that knows what was being held.
        message: String,
    },

    /// A pass that gave a calculation no variant, with the file it read and
    /// what each filter of its chain was given and kept, the outermost
    /// filter first.
    ///
    /// The core says that the reader gave no variant and no more: it is
    /// given a chain of readers and does not know whether the source held no
    /// variant or the steps kept none, and the counts of a pass that could
    /// not be finished reach nobody otherwise, as "A pass that was not
    /// finished" of `docs/specs/filters.md` says. This crate holds the chain,
    /// so it reads them from it and builds the message.
    NoVariant {
        /// The file the variants were read from.
        path: PathBuf,
        /// The kind of each filter of the pass, how many variants it was
        /// given and how many it kept, the outermost filter first.
        filtering: Vec<(&'static str, u64, u64)>,
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
            // A pass that gave no variant is a wrong input of the
            // calculation and not a result of NaN, so it is a `ValueError`,
            // whose message starts with the path as that of every error of a
            // file does.
            PyPopneiError::NoVariant { path, filtering } => {
                PyValueError::new_err(of_the_file(no_variant_message(&filtering), Some(path)))
            }
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
            // Memory the machine did not give is the `MemoryError` of
            // Python, as the memory the core asked for and was not given
            // is: the call is right and popnei is not broken, and the same
            // call on a machine with the memory free gives the result. It
            // names no file, since what could not be held is the size of
            // the calculation and not what any file holds.
            PyPopneiError::NoMemory { message } => PyMemoryError::new_err(message),
            PyPopneiError::Broken { message, path } => {
                PyRuntimeError::new_err(of_the_file(message, path))
            }
            PyPopneiError::Python(error) => error,
        }
    }
}

/// What a user reads of a pass that gave a calculation no variant: which of
/// the two it was, the source having none or the steps keeping none, and
/// what each filter was given and kept.
///
/// The filters come in the order of the chain, the outermost first, and the
/// message names them in the order of the steps, which is the one the user
/// wrote them in. The last of the chain is the innermost, the filter the
/// source feeds, so a pass whose innermost filter was given no variant is a
/// pass over a source that has none.
///
/// The wording is the one `crates/popnei-js` gives a TypeScript user, word
/// for word, so that the two languages say the same of the same pass:
/// `the source has no variant, and a calculation needs 1 variant at least`,
/// or `the steps kept no variant of the N the source gave, and a calculation
/// needs 1 variant at least`, and after either, when the pass has filters,
/// `: the filter `kind` was given n variants and kept m`, joined with `, `.
fn no_variant_message(filtering: &[(&'static str, u64, u64)]) -> String {
    // The last filter of the chain is the innermost, the one the source
    // feeds, so what it was given is what the source gave. A pass with no
    // filter reaches the calculation from the source itself.
    let from_the_source = match filtering.last() {
        Some(&(_, vars_processed, _)) => vars_processed,
        None => 0,
    };
    let what_happened = if from_the_source == 0 {
        "the source has no variant".to_owned()
    } else {
        format!("the steps kept no variant of the {from_the_source} the source gave")
    };
    if filtering.is_empty() {
        return format!("{what_happened}, and a calculation needs 1 variant at least");
    }
    let counts: Vec<String> = filtering
        .iter()
        .rev()
        .map(|&(kind, vars_processed, vars_kept)| {
            format!(
                "the filter `{kind}` was given {vars_processed} variants and kept \
                 {vars_kept}"
            )
        })
        .collect();
    format!(
        "{what_happened}, and a calculation needs 1 variant at least: {counts}",
        counts = counts.join(", ")
    )
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
/// which a file whose content is not what a VCF holds is. A `MemoryError`
/// is the fourth, for the memory a calculation asked this machine for and
/// was not given, which is neither of the three: the call is right, popnei
/// is not broken, and the same call on a machine with the memory free
/// gives the result.
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
        popnei::Error::GtsNotWholeGenotypes { .. }
        | popnei::Error::MoreAllelesThanACountHolds { .. }
        | popnei::Error::AlleleBelowTheMissingOne { .. }
        | popnei::Error::BlocksDoNotFitTogether { .. }
        | popnei::Error::BlockArrayOfAnotherSize { .. }
        | popnei::Error::ReaderGaveABlockOfNoVariants
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
        | popnei::Error::VarFilterOfTheLdCriterion => {
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
        // The arguments a user writes: how many variants a block holds,
        // and how many alleles a genotype of the file has, which the reader
        // is given when the file is opened because it needs it to read the
        // first genotype. The three of `docs/specs/filters.md` are of the
        // same kind: the threshold of a filter that is not a number from 0
        // to 1, a second filter of a kind the variants are filtered by
        // already, and a window of the filter by linkage disequilibrium
        // that is no base pairs wide, all three of which a user gets at the
        // call that adds the filter. What is wrong with them is wrong
        // whatever file is read, so they name no file although some of them
        // are refused while one is being opened.
        popnei::Error::BlockOfNoVariants
        | popnei::Error::BlockTooLarge { .. }
        | popnei::Error::VcfPloidyOutOfRange { .. }
        | popnei::Error::VarFilterThresholdOutOfRange { .. }
        | popnei::Error::VarFilterOfAKindThatIsSet { .. }
        | popnei::Error::LdFilterMaxDistTooSmall { .. }
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
        // The four of the r² of a set of variants that are wrong whatever
        // file is read: an index that is not an individual of the dataset
        // and one given twice, which are the individuals of a population
        // as a user writes them; and the two sizes the calculation cannot
        // be done at, dosages of more values than the linear algebra
        // counts in and a variant of more alleles than the sums come out
        // of exactly. What a user does about those two is calculate over
        // fewer variants or over fewer individuals, whichever file they
        // read.
        | popnei::Error::LdIndividualNotInTheDataset { .. }
        | popnei::Error::LdIndividualAskedForTwice { .. }
        | popnei::Error::LdDosagesTooLarge { .. }
        | popnei::Error::LdTooManyAllelesInAVariant { .. }
        // The two of the `max_num_vars` of the matrix of every pair, which
        // is the one number a user writes at that call: a pass that gave
        // more variants than it, and a cap of more variants than this
        // machine counts the pairs of. Both are found while a file is being
        // read, and neither is about the file: what a user does about the
        // first is raise the cap or filter the variants, whichever file
        // they read.
        | popnei::Error::LdTooManyVars { .. }
        | popnei::Error::LdMaxNumVarsTooLarge { .. } => PyValueError::new_err(message),
        // The matrix of the r², or one of the matrices it is worked out
        // through, that this machine did not give the memory of, which is
        // the `MemoryError` Python has for an allocation that was not
        // given. It is not a `ValueError`: the call says what the user
        // meant, and it runs on a machine that has the memory free, so a
        // caller who catches it takes fewer variants or fewer individuals
        // instead of looking for what they typed wrong. It names no file,
        // since what it refuses is the size of the calculation and not
        // what any file holds.
        popnei::Error::LdNoMemory { .. } => PyMemoryError::new_err(message),
        // The five of the principal components of the variants that the
        // dataset a user gave is wrong for: no variants, which the steps of
        // a `Variants` can leave; no variant with variance, which one
        // individual gives; a variant of more than two different alleles
        // among its called genotypes with `transform_to_biallelic` false; a
        // source of no individual, which is nobody to place on the axes and
        // which no source of popnei is, since one that names no individual
        // is refused when it is opened; and a dataset of a size the
        // analysis cannot count in, which "Errors and the cases pyNei
        // asserts" of `docs/specs/pca.md` lists. Each names the file the
        // variants were read from, as every error of a file does.
        popnei::Error::PcaNoVariants
        | popnei::Error::PcaNoVariantWithVariance
        | popnei::Error::PcaVariantWithMoreThanTwoAlleles { .. }
        | popnei::Error::PcaNoIndividual
        | popnei::Error::PcaVariantsTooLarge { .. } => {
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
        // alleles in one block is more memory than a machine gives.
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
