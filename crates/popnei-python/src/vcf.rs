//! What a Python user reaches through `open_vcf` and `Variants`: a VCF with
//! its options, and the blocks of its variants.
//!
//! Two classes. [`VcfSource`] holds the path of a VCF and the options it is
//! read with, and it reads the header when it is built, so a file that is
//! not a VCF fails at `open_vcf`. [`Blocks`] is one pass over that file: it
//! owns a reader of blocks of the core with a `Reblock` at its end, which
//! gives the blocks the size that was asked for, and every call of
//! `VcfSource::blocks` opens the file again, which is what lets a user give
//! the same `Variants` to one calculation after another.
//!
//! The columns of a block leave as they are in the core: the genotypes as a
//! numpy array that holds the allocation the core filled, the positions and
//! the qualities the same way, and the chromosomes, the ids and the alleles
//! as tuples. The Python package builds the `Block` dataclass out of them,
//! as `docs/specs/block.md` describes it.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Mutex;

use numpy::ndarray::Array3;
use numpy::{IntoPyArray, PyArray1, PyArray3};
use pyo3::prelude::*;
use pyo3::types::{PyString, PyTuple};

use popnei::block::{
    AllelesColumn, Block, BlockReader, CollectedBlocks, Reblock, needs_of_the_fields,
};
use popnei::io::vcf::{VcfOptions, VcfReader};
use popnei::variant::{ChromTable, Needs};

use crate::errors::PyPopneiError;

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

// A VCF that was opened: its path, the options it is read with, and the
// individuals its header named. A `///` here would become the `__doc__` of
// the class, and what a Python user reads belongs to the package, which is
// the API.
#[pyclass(frozen, module = "popnei._core")]
pub(crate) struct VcfSource {
    path: PathBuf,
    options: VcfOptions,
    individuals: Vec<String>,
}

#[pymethods]
impl VcfSource {
    // The names of the individuals, in the order of the columns of the
    // VCF.
    fn individuals(&self) -> Vec<String> {
        self.individuals.clone()
    }

    // How many alleles the genotype of one individual holds.
    fn ploidy(&self) -> usize {
        self.options.ploidy
    }

    // One pass over the file: it is opened again, and its blocks hold
    // `fields` besides the genotypes, `num_vars_per_block` variants each.
    #[pyo3(signature = (fields, num_vars_per_block))]
    fn blocks(
        &self,
        py: Python<'_>,
        fields: Vec<String>,
        num_vars_per_block: Option<i64>,
    ) -> Result<Blocks, PyPopneiError> {
        let needs = needs_of_the_fields(fields.iter().map(String::as_str))?;
        let num_vars_per_block = num_vars_per_block
            .map(|asked_for| count_of("num_vars_per_block", asked_for))
            .transpose()?;
        let path = &self.path;
        let options = self.options;
        let reader = py
            .detach(|| -> Result<Box<dyn BlockReader>, popnei::Error> {
                // The source is asked for the size the user wants, so the
                // `Reblock` over it has nothing to cut or to join and every
                // block goes through with no copy. It is there for the
                // sources that give another size, a filter among them, and
                // it is what `docs/specs/block.md` puts at the end of every
                // `iter_blocks`. The source today is the adapter of the
                // core over the reader of single variants; task 2.4 of
                // `docs/plans/block-readers.md` puts the VCF reader itself
                // there, which is the line that builds `blocks`.
                let reader = VcfReader::from_path(path, options)?;
                let blocks = CollectedBlocks::new(reader, needs, num_vars_per_block)?;
                Ok(Box::new(Reblock::new(blocks, num_vars_per_block)?))
            })
            .map_err(|error| PyPopneiError::of_the_file(error, path))?;
        Ok(Blocks {
            reader: Mutex::new(reader),
            path: path.clone(),
        })
    }
}

// One pass over a VCF, which gives its variants block by block.
#[pyclass(frozen, module = "popnei._core")]
pub(crate) struct Blocks {
    reader: Mutex<Box<dyn BlockReader>>,
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
        // who asks for the blocks of a big VCF waits here, and a Ctrl-C
        // between two blocks is how they stop.
        py.check_signals()?;
        let Some((block, chroms)) = py.detach(|| self.next_block())? else {
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
                PyPopneiError::Broken(format!(
                    "a block of {num_vars} variants of {num_individuals} individuals of \
                     the ploidy {ploidy} does not hold that many genotypes: {error}"
                ))
            })?;
        let gts = read_only(gts.into_pyarray(py))?;
        let chrom = chroms.map(|chroms| chrom_column(py, &chroms)).transpose()?;
        let id = id.map(|ids| id_column(py, &ids)).transpose()?;
        let alleles = alleles
            .map(|column| alleles_column(py, &column))
            .transpose()?;
        let pos = pos.map(|pos| read_only(pos.into_pyarray(py))).transpose()?;
        let qual = qual
            .map(|qual| read_only(qual.into_pyarray(py)))
            .transpose()?;
        Ok(Some((gts, chrom, pos, id, alleles, qual)))
    }
}

impl Blocks {
    /// The next block of the reader, with the chromosomes of its variants,
    /// or `None` when the VCF has no more variants.
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
        let mut reader = self.reader.lock().map_err(|_| {
            PyPopneiError::Broken(
                "the blocks of this pass cannot be read any more: a panic left the \
                 reader half way through a block"
                    .to_string(),
            )
        })?;
        let Some(block) = reader
            .next_block()
            .map_err(|error| PyPopneiError::of_the_file(error, &self.path))?
        else {
            return Ok(None);
        };
        block.check()?;
        let chroms = match block.chrom.as_deref() {
            Some(numbers) => Some(ChromColumn::of(numbers, reader.chroms())?),
            None => None,
        };
        Ok(Some((block, chroms)))
    }
}

// The VCF at `path`, read with `ploidy` alleles in every genotype and, when
// `only_passed` is true, without the variants that failed a filter. It
// reads the header, so the individuals are known when it returns.
#[pyfunction]
pub(crate) fn open_vcf(
    py: Python<'_>,
    path: PathBuf,
    ploidy: i64,
    only_passed: bool,
) -> Result<VcfSource, PyPopneiError> {
    let options = VcfOptions {
        ploidy: count_of("ploidy", ploidy)?,
        only_passed,
    };
    let individuals = py
        .detach(|| -> Result<_, popnei::Error> {
            // The header is read when the reader is built, and the names it
            // gave are asked of the adapter, which is what this crate holds
            // a reader through. Task 2.4 of
            // `docs/plans/block-readers.md` puts the VCF reader itself in
            // the place of the adapter, and it answers this itself. No
            // variant is read here, so the fields and the size of a block
            // are the ones that ask for nothing.
            let reader = VcfReader::from_path(&path, options)?;
            let reader = CollectedBlocks::new(reader, Needs::GTS, None)?;
            Ok(reader.individuals().to_vec())
        })
        .map_err(|error| PyPopneiError::of_the_file(error, &path))?;
    Ok(VcfSource {
        path,
        options,
        individuals,
    })
}

/// The `value` that was given for the argument `name`, as a number of
/// things: the one place where an argument that counts something crosses
/// from Python, so that a number that counts nothing is refused by the name
/// a user wrote and not by the conversion, whose `OverflowError` names
/// neither the argument nor what is wrong with it.
fn count_of(name: &'static str, value: i64) -> Result<usize, PyPopneiError> {
    usize::try_from(value).map_err(|_| PyPopneiError::Count { name, value })
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
    fn of(numbers: &[u32], chroms: &ChromTable) -> Result<ChromColumn, PyPopneiError> {
        let mut names = Vec::new();
        let mut of_each_variant = Vec::with_capacity(numbers.len());
        let mut where_each_number_went: HashMap<u32, usize> = HashMap::new();
        for number in numbers {
            let index = match where_each_number_went.get(number) {
                Some(index) => *index,
                None => {
                    let Some(name) = chroms.name(*number) else {
                        return Err(PyPopneiError::Broken(format!(
                            "the chromosome number {number} of a block is not in the \
                             table of the reader that gave it"
                        )));
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
                PyPopneiError::Broken(format!(
                    "the chromosome {index} of a block has no name beside it"
                ))
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
