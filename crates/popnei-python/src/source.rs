//! What every source of variants shares: one pass over it, and the columns
//! of the blocks that pass gives Python.
//!
//! A source is a file with what is needed to read it, a VCF with its
//! options in `vcf.rs` and a vars file in `vars.rs`. Each is a class of
//! Python that the package's `Variants` holds, each reads what says what the
//! file holds when it is built, so a file that is not of its format fails at
//! `open_vcf` or at `open_vars`, and each opens the file again at every
//! pass, which is what lets a user give the same `Variants` to one
//! calculation after another. What they have in common is [`OpenSource`]:
//! the path the errors of a pass name, and the reader of one pass.
//!
//! [`Blocks`] is that pass, whichever source it came from: it owns the chain
//! of readers of the pass, the source with a filter over it for each step of
//! the `Variants` and a `Reblock` at its end, which gives the blocks the
//! size that was asked for. It counts the variants of the blocks it gives
//! and reads the counts of the filters from that chain, which is the
//! `PassStats` of `docs/specs/variant.md`.
//!
//! The columns of a block leave as they are in the core: the genotypes as a
//! numpy array that holds the allocation the core filled, the positions and
//! the qualities the same way, and the chromosomes, the ids and the alleles
//! as tuples. The Python package builds the `Block` dataclass out of them,
//! as `docs/specs/block.md` describes it.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use numpy::ndarray::Array3;
use numpy::{IntoPyArray, PyArray1, PyArray3};
use pyo3::exceptions::{PyOverflowError, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyString, PyTuple};

use popnei::block::{AllelesColumn, Block, BlockReader, Reblock, needs_of_the_fields};
use popnei::variant::{ChromTable, Needs};

use crate::errors::PyPopneiError;
use crate::steps::{Steps, chain_of};
use crate::vars::VarsSource;
use crate::vcf::VcfSource;

/// The counts of one pass on their way to Python: how many variants it has
/// given, and, for each filter of its chain, its kind, how many variants it
/// was given and how many it kept.
///
/// The filters come in the order of the chain, the outermost first, which is
/// the reverse of the order of the steps: the package turns them around, as
/// "How it runs" of the counts of `docs/specs/filters.md` says.
pub(crate) type PassCounts = (u64, Vec<(&'static str, u64, u64)>);

/// The columns of one block on their way to Python: the genotypes, and then
/// the chromosomes, the positions, the ids, the alleles and the qualities,
/// each `None` when it was not asked for.
type BlockColumns<'py> = (
    Bound<'py, PyArray3<i8>>,
    Option<Bound<'py, PyTuple>>,
    Option<Bound<'py, PyArray1<u64>>>,
    Option<Bound<'py, PyTuple>>,
    Option<Bound<'py, PyTuple>>,
    Option<Bound<'py, PyArray1<f32>>>,
);

/// A file of variants that was opened, which every pass reads again.
///
/// It is `Sync` because a source is a frozen class of Python, which any
/// thread can hold, and because the reader of a pass is built with the
/// interpreter released.
pub(crate) trait OpenSource: Sync {
    /// The file the source reads, which the errors of a pass over it name.
    fn path(&self) -> &Path;

    /// The reader of one pass over the source, which opens the file again.
    ///
    /// `num_vars_per_block` is the size the caller will ask the blocks for,
    /// which a source whose reader can give them at that size is built
    /// with, so that the `Reblock` over it has nothing to cut or to join. A
    /// source that gives its blocks as they are in the file leaves it to
    /// that `Reblock`.
    ///
    /// # Errors
    ///
    /// When the file cannot be opened or what says what it holds, the
    /// header of a VCF or the schema of a vars file, cannot be read.
    fn reader(
        &self,
        num_vars_per_block: Option<usize>,
    ) -> Result<Box<dyn BlockReader>, popnei::Error>;
}

/// The source that `object` is, which is one of the two classes a
/// `Variants` of the package holds.
///
/// # Errors
///
/// When it is neither, which the package refuses before it calls this crate
/// and which a user reaches by calling `popnei._core` themselves.
pub(crate) fn source_of<'a>(
    object: &'a Bound<'_, PyAny>,
) -> Result<&'a dyn OpenSource, PyPopneiError> {
    if let Ok(vcf) = object.cast::<VcfSource>() {
        return Ok(vcf.get());
    }
    if let Ok(vars) = object.cast::<VarsSource>() {
        return Ok(vars.get());
    }
    // The type is named and no article is put before it: `a int` and `a
    // NoneType` are what one written here would give.
    Err(PyTypeError::new_err(format!(
        "`source` is of the type `{what}`, and the variants to write come from a source \
         that popnei opened: give what `open_vcf` or `open_vars` gives",
        what = object.get_type().name()?
    ))
    .into())
}

/// One pass over `source`, through the steps of `steps`, whose blocks hold
/// `fields` besides the genotypes, `num_vars_per_block` variants each.
///
/// The steps are taken as they are here, when the pass starts: one added
/// while it runs holds from the next pass, as `docs/specs/filters.md` says.
///
/// # Errors
///
/// When a name of `fields` is of no field and when `num_vars_per_block`
/// counts no variant, which are refused before the file is opened, and when
/// the source cannot be opened or read.
pub(crate) fn blocks_of(
    py: Python<'_>,
    source: &dyn OpenSource,
    fields: Vec<String>,
    num_vars_per_block: Option<&Bound<'_, PyAny>>,
    steps: &Steps,
) -> Result<Blocks, PyPopneiError> {
    let needs = needs_of_the_fields(fields.iter().map(String::as_str))?;
    let num_vars_per_block = num_vars_per_block
        .map(|asked_for| count_of("num_vars_per_block", asked_for))
        .transpose()?;
    let steps = steps.of_a_pass()?;
    let path = source.path();
    let reader = py
        .detach(|| -> Result<Box<dyn BlockReader>, popnei::Error> {
            // The source is asked for the size the user wants, so a reader
            // that can give it has nothing for the `Reblock` over it to cut
            // or to join and every block goes through with no copy. That
            // `Reblock` is there for the sources that give another size,
            // the vars file whose batches were written at one size and a
            // filter among them, and it is what `docs/specs/block.md` puts
            // at the end of every `iter_blocks`.
            let reader = source.reader(num_vars_per_block)?;
            // The fields are asked of the whole chain and not of the source
            // alone: a filter asks its source for what it was asked for and
            // for the genotypes, which it needs itself.
            let mut chain = chain_of(reader, &steps)?;
            chain.set_needs(needs.union(Needs::GTS));
            Ok(Box::new(Reblock::new(chain, num_vars_per_block)?))
        })
        .map_err(|error| PyPopneiError::of_the_file(error, path))?;
    Ok(Blocks {
        pass: Mutex::new(Pass {
            reader,
            finished: false,
            num_vars: 0,
        }),
        path: path.to_path_buf(),
    })
}

/// The reader of one pass, whether the pass is over and how many variants
/// it has given: they are read and written together, under one lock,
/// because a pass that is over gives no block whatever its reader would
/// say, and the count is of the blocks that reader gave.
struct Pass {
    reader: Box<dyn BlockReader>,
    /// Whether the reader has no more blocks or a block was lost with an
    /// error. After either there is no block.
    finished: bool,
    /// The variants of the blocks the pass has given, which is the
    /// `num_vars` a user reads in its counts. A block that was lost with an
    /// error is not among them: it never reached the user.
    num_vars: u64,
}

// One pass over a source of variants, which gives them block by block.
#[pyclass(frozen, module = "popnei._core")]
pub(crate) struct Blocks {
    pass: Mutex<Pass>,
    /// The file the reader reads, for the errors of the file system, which
    /// carry it where Python keeps it, `OSError.filename`.
    path: PathBuf,
}

#[pymethods]
impl Blocks {
    fn __iter__(this: PyRef<'_, Self>) -> PyRef<'_, Self> {
        this
    }

    fn __next__<'py>(&self, py: Python<'py>) -> Result<Option<BlockColumns<'py>>, PyPopneiError> {
        // A block of the default size is a few million genotypes, so a user
        // who asks for the blocks of a big file waits here, and a Ctrl-C
        // between two blocks is how they stop. It costs no block, so the
        // pass is not over after it; everything from the read on loses the
        // block it happened in, and a pass that went on would give the
        // variants that follow as if nothing had happened, which
        // `docs/specs/block.md` asks of every reader that it not do.
        py.check_signals()?;
        self.columns_of_the_next_block(py)
            .inspect_err(|_| self.finish(py))
    }

    // How many variants the pass has given, and what each filter of it was
    // given and kept, the outermost filter first. It is read while the pass
    // runs too, and it then holds what has been read up to there.
    fn pass_stats(&self, py: Python<'_>) -> Result<PassCounts, PyPopneiError> {
        // The interpreter is released while the lock is waited for: the
        // thread that reads a block holds that lock for the whole read,
        // inside its own `detach`, and a thread that waited for it with the
        // interpreter in hand would stop every other thread of the process
        // for as long as that read takes, seconds for a block of a big
        // file.
        py.detach(|| {
            let pass = self.pass.lock().map_err(|_| {
                PyPopneiError::broken_of_the_file(
                    "the counts of this pass cannot be read: a panic left the reader half \
                     way through a block"
                        .to_string(),
                    &self.path,
                )
            })?;
            let filtering = pass
                .reader
                .filtering_stats()
                .into_iter()
                .map(|(kind, stats)| (kind, stats.vars_processed, stats.vars_kept))
                .collect();
            Ok((pass.num_vars, filtering))
        })
    }
}

impl Blocks {
    /// The columns of the next block, or `None` when the source has no more
    /// variants. What ends the pass at an error is [`Blocks::__next__`],
    /// which calls this one.
    fn columns_of_the_next_block<'py>(
        &self,
        py: Python<'py>,
    ) -> Result<Option<BlockColumns<'py>>, PyPopneiError> {
        let read = py.detach(|| self.next_block())?;
        // A Ctrl-C that arrived while the block was read is still pending:
        // the interpreter was released and no bytecode ran to raise it. It
        // is raised here, before numpy is called, because the first array of
        // a process imports the C API of numpy, that import fails with the
        // exception that is pending, and the numpy crate panics when it
        // does: a user who asked for a Ctrl-C would get a `PanicException`,
        // which no `except` of theirs catches and which ends the session.
        py.check_signals()?;
        let Some((block, chroms)) = read else {
            return Ok(None);
        };
        let Block {
            num_vars,
            num_individuals,
            ploidy,
            gts,
            // The numbers of the chromosomes were turned into their names
            // while the reader that holds the table was at hand.
            chrom: _,
            pos,
            id,
            alleles,
            qual,
        } = block;
        // The three dimensions are given to the array as it is built and
        // not by reshaping one of a single dimension, whose array would
        // stay under the block's as a writable view of the same genotypes.
        let gts =
            Array3::from_shape_vec((num_vars, num_individuals, ploidy), gts).map_err(|error| {
                PyPopneiError::broken_of_the_file(
                    format!(
                        "a block of {num_vars} variants of {num_individuals} individuals \
                         of the ploidy {ploidy} does not hold that many genotypes: {error}"
                    ),
                    &self.path,
                )
            })?;
        let gts = read_only(gts.into_pyarray(py))?;
        let chrom = chroms
            .map(|chroms| chrom_column(py, &chroms, &self.path))
            .transpose()?;
        let id = id.map(|ids| id_column(py, &ids)).transpose()?;
        let alleles = alleles
            .map(|column| alleles_column(py, &column))
            .transpose()?;
        let pos = pos.map(|pos| read_only(pos.into_pyarray(py))).transpose()?;
        let qual = qual
            .map(|qual| read_only(qual.into_pyarray(py)))
            .transpose()?;
        self.counted(py, num_vars);
        Ok(Some((gts, chrom, pos, id, alleles, qual)))
    }

    /// The `num_vars` variants of a block that is going out, added to the
    /// count of the pass.
    ///
    /// They are counted here, where the block is the user's, and not where
    /// it was read: a block that was lost with an error, and one the pass
    /// had read when a Ctrl-C arrived, never reached them. The counts of
    /// the filters are another matter, since a filter sits under the
    /// `Reblock` and has counted that block already, which "A pass that was
    /// not finished" of `docs/specs/filters.md` says a user reads.
    ///
    /// A lock that a panic left broken is the end of the pass, which the
    /// read of the next block reports: a count that was not added is not
    /// what the user is told about then.
    ///
    /// The interpreter is released while the lock is waited for, as it is
    /// in [`Blocks::pass_stats`]: another thread may hold it for a whole
    /// block.
    fn counted(&self, py: Python<'_>, num_vars: usize) {
        py.detach(|| self.count(num_vars));
    }

    /// The `num_vars` variants added to the count, with the interpreter
    /// already released.
    fn count(&self, num_vars: usize) {
        if let Ok(mut pass) = self.pass.lock() {
            // A variant is a row of a file, so a pass of the
            // 18446744073709551615 variants this count holds is more rows
            // than any file system takes: the sum cannot reach its end. The
            // conversion cannot fail either: a `usize` is 64 bits natively
            // and 32 in wasm, and both fit in a `u64`.
            pass.num_vars = pass
                .num_vars
                .saturating_add(u64::try_from(num_vars).unwrap_or(u64::MAX));
        }
    }

    /// The pass is over, and every call after this one gives no block.
    ///
    /// A lock that a panic left broken is already the end of the pass: every
    /// read of it is the error of a reader that cannot be read any more.
    ///
    /// The interpreter is released while the lock is waited for, as it is in
    /// [`Blocks::pass_stats`]: another thread may hold it for a whole block.
    fn finish(&self, py: Python<'_>) {
        py.detach(|| {
            if let Ok(mut pass) = self.pass.lock() {
                pass.finished = true;
            }
        });
    }

    /// The next block of the reader, with the chromosomes of its variants,
    /// or `None` when the source has no more variants.
    ///
    /// The names are taken after the block was given, as
    /// `docs/specs/block.md` says: the table of the reader grows while the
    /// file is read, and a block holds numbers of it. They are copied
    /// because the interpreter is released here and the lock is not held
    /// while the Python objects are built: a thread that waited for the
    /// lock with the interpreter in hand would never get it back from this
    /// one.
    ///
    /// The block is checked here, before its genotypes go to numpy as one
    /// array of variants x individuals x ploidy: a block whose arrays are
    /// not of its size would be read one genotype at the place of another,
    /// with nothing to show it.
    fn next_block(&self) -> Result<Option<(Block, Option<ChromColumn>)>, PyPopneiError> {
        let mut pass = self.pass.lock().map_err(|_| {
            PyPopneiError::broken_of_the_file(
                "the blocks of this pass cannot be read any more: a panic left the \
                 reader half way through a block"
                    .to_string(),
                &self.path,
            )
        })?;
        if pass.finished {
            return Ok(None);
        }
        let block = pass.next_block(&self.path);
        if !matches!(block, Ok(Some(_))) {
            pass.finished = true;
        }
        block
    }
}

impl Pass {
    /// The next block of the reader, with the chromosomes of its variants,
    /// read with `path` at hand for the errors of the file system. What ends
    /// the pass is [`Blocks::next_block`], which calls this one.
    fn next_block(
        &mut self,
        path: &Path,
    ) -> Result<Option<(Block, Option<ChromColumn>)>, PyPopneiError> {
        let Some(block) = self
            .reader
            .next_block()
            .map_err(|error| PyPopneiError::of_the_file(error, path))?
        else {
            return Ok(None);
        };
        block
            .check()
            .map_err(|error| PyPopneiError::of_the_file(error, path))?;
        let chroms = match block.chrom.as_deref() {
            Some(numbers) => Some(ChromColumn::of(numbers, self.reader.chroms(), path)?),
            None => None,
        };
        Ok(Some((block, chroms)))
    }
}

/// The `value` that was given for the argument `name`, as a number of
/// things: the one place where an argument that counts something crosses
/// from Python.
///
/// The object is taken as it is and converted here, and not by the
/// signature, because an integer of Python is of any size: converting it in
/// the signature raises the `OverflowError` of pyo3, "Python int too large
/// to convert to C long", which names neither the argument nor what is
/// wrong with it, and it does so before any code of ours runs.
///
/// # Errors
///
/// When the object is a whole number that counts nothing, a negative one or
/// one above what this machine counts, which is the error that names the
/// argument and the value. An object that is not a whole number at all,
/// `2.5` or `"two"`, keeps the `TypeError` of pyo3, which says what it was
/// given.
pub(crate) fn count_of(
    name: &'static str,
    value: &Bound<'_, PyAny>,
) -> Result<usize, PyPopneiError> {
    match value.extract::<usize>() {
        Ok(count) => Ok(count),
        Err(error) if error.is_instance_of::<PyOverflowError>(value.py()) => {
            Err(PyPopneiError::Count {
                name,
                value: value.to_string(),
            })
        }
        Err(error) => Err(error.into()),
    }
}

/// The array, which nothing writes into any more: a block is frozen, and
/// its arrays hold the memory the core filled.
fn read_only<'py, T>(array: Bound<'py, T>) -> Result<Bound<'py, T>, PyPopneiError> {
    array
        .as_any()
        .getattr("flags")?
        .setattr("writeable", false)?;
    Ok(array)
}

/// The chromosomes of the variants of one block: the name of each
/// chromosome the block holds, once, and which of those names each variant
/// has.
struct ChromColumn {
    names: Vec<String>,
    /// One index into `names` for each variant of the block.
    of_each_variant: Vec<usize>,
}

impl ChromColumn {
    /// The chromosomes that `numbers`, the column of a block, names in
    /// `chroms`, the table of the reader that filled it.
    ///
    /// The table of a de novo assembly holds 10^4 scaffolds or more and a
    /// block holds a few of them, so what is copied is the name of every
    /// chromosome of the block and not the table.
    fn of(numbers: &[u32], chroms: &ChromTable, path: &Path) -> Result<ChromColumn, PyPopneiError> {
        let mut names = Vec::new();
        let mut of_each_variant = Vec::with_capacity(numbers.len());
        let mut where_each_number_went: HashMap<u32, usize> = HashMap::new();
        for number in numbers {
            let index = match where_each_number_went.get(number) {
                Some(index) => *index,
                None => {
                    let Some(name) = chroms.name(*number) else {
                        return Err(PyPopneiError::broken_of_the_file(
                            format!(
                                "the chromosome number {number} of a block is not in \
                                 the table of the reader that gave it"
                            ),
                            path,
                        ));
                    };
                    names.push(name.to_owned());
                    let index = names.len().saturating_sub(1);
                    where_each_number_went.insert(*number, index);
                    index
                }
            };
            of_each_variant.push(index);
        }
        Ok(ChromColumn {
            names,
            of_each_variant,
        })
    }
}

/// The name of the chromosome of every variant of a block.
fn chrom_column<'py>(
    py: Python<'py>,
    chroms: &ChromColumn,
    path: &Path,
) -> Result<Bound<'py, PyTuple>, PyPopneiError> {
    // One Python string for each chromosome, which the variants of that
    // chromosome share: a block of 10000 variants of one chromosome holds
    // one name and not 10000.
    let names: Vec<Bound<'py, PyString>> = chroms
        .names
        .iter()
        .map(|name| PyString::new(py, name))
        .collect();
    let of_each_variant = chroms
        .of_each_variant
        .iter()
        .map(|index| {
            names.get(*index).cloned().ok_or_else(|| {
                PyPopneiError::broken_of_the_file(
                    format!("the chromosome {index} of a block has no name beside it"),
                    path,
                )
            })
        })
        .collect::<Result<Vec<_>, PyPopneiError>>()?;
    Ok(PyTuple::new(py, of_each_variant)?)
}

/// The id of every variant of a block, `None` for a variant that has none,
/// which the core gives as an empty id.
fn id_column<'py>(py: Python<'py>, ids: &[String]) -> Result<Bound<'py, PyTuple>, PyPopneiError> {
    Ok(PyTuple::new(
        py,
        ids.iter().map(|id| (!id.is_empty()).then_some(id.as_str())),
    )?)
}

/// The alleles of every variant of a block, the reference one first.
fn alleles_column<'py>(
    py: Python<'py>,
    column: &AllelesColumn,
) -> Result<Bound<'py, PyTuple>, PyPopneiError> {
    let of_each_variant = (0..column.num_vars())
        .map(|var| {
            PyTuple::new(
                py,
                (0..column.num_alleles(var)).map(|allele| column.allele(var, allele)),
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyTuple::new(py, of_each_variant)?)
}
