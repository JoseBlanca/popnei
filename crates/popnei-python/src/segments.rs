//! The writer of segments: a vars file written from the dicts of arrays that
//! a user's code yields, one genomic segment at a time, "The writer of
//! segments" of `docs/specs/io_vars.md`.
//!
//! [`SegmentReader`] is the source of the pass. It takes each item from the
//! user's iterator when the core asks for a block, checks it, and builds a
//! block with every column and every declared field, so that every block of
//! the pass holds the same ones. The pass is `write_vars` of the core, as
//! for any other source, inside the file that [`write_a_new_file`] makes
//! and takes away when the pass fails.
//!
//! The core cannot carry an exception of Python. When the user's code
//! raises one, or a segment is refused here, the reader keeps that error and
//! ends the pass with [`popnei::Error::SegmentsFailed`]; the call then
//! raises the error it kept in the place of that one, the user's exception
//! as it was raised, with its traceback.

use std::collections::HashSet;
use std::io::BufWriter;
use std::path::PathBuf;

use numpy::{
    Element, PyArrayDescrMethods, PyArrayDyn, PyArrayMethods, PyUntypedArray, PyUntypedArrayMethods,
};
use pyo3::exceptions::{PyKeyError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyFloat, PyInt, PyIterator, PyMapping, PyString};

use popnei::block::{
    AllelesColumn, Block, BlockReader, ExtraField, ExtraFieldDecl, ExtraKind, ExtraType,
    ExtraValues, SourceHeader,
};
use popnei::filters::{FilteringStats, RegionSelection};
use popnei::io::vcf::MAX_PLOIDY;
use popnei::variant::{ChromTable, Needs};

use crate::errors::PyPopneiError;
use crate::source::{PassCounts, count_of};
use crate::vars::{of_the_file_it_is_about, to_the_disc, write_a_new_file};

/// The most alleles one variant of a segment can have, numbered 0 to 127,
/// the most that the byte of an allele of `gts` holds. The owner kept it on
/// 9 November 2026, "A call that fails" of the writer of segments of
/// `docs/specs/io_vars.md`; alleles of two bytes would double the genotypes
/// in memory and in the file.
const MAX_ALLELES_OF_A_VARIANT: usize = 128;

/// The keys of a segment that are not a field, in the order of the table of
/// "Its Python function" of the writer of segments.
const KEYS: [&str; 6] = ["chrom", "pos", "alleles", "gts", "id", "qual"];

// Every variant of the segments that `segments` yields into a vars file made
// at `path`, of `individuals` of the ploidy `ploidy`, with the extra fields
// `var_fields` and `gt_fields`, each a name and the name numpy gives its
// type, in the order of the file; the lengths of the chromosomes, batches
// of `num_vars_per_block` variants, `None` for the size popnei chooses, and
// the check of the order unless `allow_unsorted` is true. What it gives back
// is the counts of the pass, the variants written and no filter. A `///`
// comment here would become the `__doc__` of
// `popnei._core.write_vars_from_segments`, and what a Python user reads
// belongs to the package, which is the API.
#[pyfunction]
#[pyo3(signature = (
    segments, path, individuals, ploidy, var_fields, gt_fields, chrom_lengths,
    num_vars_per_block, allow_unsorted,
))]
#[expect(
    clippy::too_many_arguments,
    reason = "the arguments of `write_vars_from_segments` of the package, which the spec \
              gives, each passed on as the package received it"
)]
pub(crate) fn write_vars_from_segments(
    py: Python<'_>,
    segments: &Bound<'_, PyAny>,
    path: PathBuf,
    individuals: Vec<String>,
    ploidy: &Bound<'_, PyAny>,
    var_fields: Vec<(String, String)>,
    gt_fields: Vec<(String, String)>,
    chrom_lengths: Option<Vec<(String, Bound<'_, PyAny>)>>,
    num_vars_per_block: Option<&Bound<'_, PyAny>>,
    allow_unsorted: bool,
) -> Result<PassCounts, PyPopneiError> {
    // Everything that goes into the schema is checked before the file is
    // made, so that a wrong argument leaves the path as it was.
    individuals_fit(&individuals)?;
    let ploidy = count_of("ploidy", ploidy)?;
    // `count_of` leaves the 0 to the core, which would refuse it once the
    // file is made.
    if ploidy == 0 || ploidy > MAX_PLOIDY {
        return Err(popnei::Error::VcfPloidyOutOfRange {
            ploidy,
            largest: MAX_PLOIDY,
        }
        .into());
    }
    let mut fields = decls_of(&var_fields, ExtraKind::Variant)?;
    fields.extend(decls_of(&gt_fields, ExtraKind::Genotype)?);
    popnei::io::vars::check_extra_fields(&fields)?;
    let chrom_lengths = chrom_lengths
        .unwrap_or_default()
        .iter()
        .map(|(chrom, length)| {
            name_fits("chrom_lengths", "a chromosome", chrom)?;
            Ok((chrom.clone(), length_of(chrom, length)?))
        })
        .collect::<Result<Vec<(String, u64)>, PyPopneiError>>()?;
    let num_vars_per_block = num_vars_per_block
        .map(|asked_for| count_of("num_vars_per_block", asked_for))
        .transpose()?;
    // `count_of` leaves the 0 to the core, whose message would not name
    // the argument.
    if num_vars_per_block == Some(0) {
        return Err(PyPopneiError::Count {
            name: "num_vars_per_block",
            smallest: 1,
            value: "0".to_owned(),
        });
    }
    // An object that cannot be iterated is the `TypeError` of Python,
    // before the file is made.
    let segments = segments.try_iter()?.unbind();
    let numpy = Numpy::new(py)?;
    let reader = SegmentReader {
        segments,
        numpy,
        header: SourceHeader {
            individuals,
            chrom_lengths,
            vcf_meta_lines: None,
            keeps_passed: false,
            extra_fields: fields,
        },
        ploidy,
        chroms: ChromTable::new(),
        num_segments: 0,
        kept: None,
        ended: false,
    };
    write_a_new_file(py, &path, |file| {
        let mut reader = reader;
        let written = popnei::io::vars::write_vars(
            &mut reader,
            BufWriter::new(file),
            num_vars_per_block,
            allow_unsorted,
        );
        let (sink, num_vars) = match written {
            Ok(written) => written,
            // What the reader kept is what ended the pass: the exception of
            // the user's code, a segment refused or a Ctrl-C. Any other
            // error is the core's, and the only file of the call is the one
            // being written.
            Err(error) => {
                return Err(reader
                    .kept
                    .take()
                    .unwrap_or_else(|| of_the_file_it_is_about(error, &path, &path)));
            }
        };
        to_the_disc(sink).map_err(|error| PyPopneiError::of_the_file(error, &path))?;
        // No filter ran, so the counts are the variants written alone.
        Ok((num_vars, Vec::new(), false))
    })
}

/// That `individuals` can be the individuals of a vars file: one at least,
/// each of its own name, which [`name_fits`].
///
/// # Errors
///
/// [`PyPopneiError::Argument`] for no individual, for a name given twice
/// and for one that [`name_fits`] refuses.
fn individuals_fit(individuals: &[String]) -> Result<(), PyPopneiError> {
    if individuals.is_empty() {
        return Err(PyPopneiError::Argument {
            name: "individuals",
            problem: "names no individual, and a vars file holds the genotypes of one at least"
                .to_owned(),
        });
    }
    let mut named = HashSet::with_capacity(individuals.len());
    for name in individuals {
        name_fits("individuals", "an individual", name)?;
        if !named.insert(name.as_str()) {
            return Err(PyPopneiError::Argument {
                name: "individuals",
                problem: format!(
                    "names `{name}` twice, and each individual of a file has its own name"
                ),
            });
        }
    }
    Ok(())
}

/// That `name`, which the argument `argument` gives to `what`, is a name a
/// VCF that `write_vcf` writes from the file can hold: not empty, and with
/// no tab, newline, carriage return or space.
///
/// # Errors
///
/// [`PyPopneiError::Argument`] for a name that is not.
fn name_fits(argument: &'static str, what: &str, name: &str) -> Result<(), PyPopneiError> {
    let problem = if name.is_empty() {
        format!("names {what} whose name is an empty text")
    } else if let Some(problem) = text_problem(name, "") {
        format!("names `{name}`, which {problem}")
    } else {
        return Ok(());
    };
    Err(PyPopneiError::Argument {
        name: argument,
        problem,
    })
}

/// The fields of one kind as the core declares them, from the name of each
/// and the name numpy gives its type, which the package reads from whatever
/// spelling the user wrote.
///
/// # Errors
///
/// [`PyPopneiError::Argument`] for a type that is not one of the five.
fn decls_of(
    fields: &[(String, String)],
    kind: ExtraKind,
) -> Result<Vec<ExtraFieldDecl>, PyPopneiError> {
    fields
        .iter()
        .map(|(name, ty)| {
            let Some(ty) = ExtraType::ALL
                .into_iter()
                .find(|known| known.to_string() == *ty)
            else {
                return Err(PyPopneiError::Argument {
                    name: match kind {
                        ExtraKind::Variant => "var_fields",
                        ExtraKind::Genotype => "gt_fields",
                    },
                    problem: format!(
                        "gives the field `{name}` the type {ty}, and a field is of one of \
                         int32, int64, uint64, float32 and float64"
                    ),
                });
            };
            Ok(ExtraFieldDecl {
                name: name.clone(),
                kind,
                ty,
            })
        })
        .collect()
}

/// The length that `chrom_lengths` gives `chrom`, a whole number above 0.
///
/// # Errors
///
/// [`PyPopneiError::Argument`] for anything else, a truth value among them,
/// which Python counts as 1 or 0.
fn length_of(chrom: &str, length: &Bound<'_, PyAny>) -> Result<u64, PyPopneiError> {
    let read = if length.is_instance_of::<PyBool>() {
        None
    } else {
        length.extract::<u64>().ok().filter(|length| *length > 0)
    };
    read.ok_or_else(|| PyPopneiError::Argument {
        name: "chrom_lengths",
        problem: format!(
            "gives `{chrom}` the length {length}, and a length is a whole number of base \
             pairs above 0 that 64 bits hold"
        ),
    })
}

/// The source of the pass of the writer of segments: the blocks it builds
/// from the items of the user's iterator.
///
/// It holds Python objects and takes the interpreter back for each item it
/// asks for, so it is `Send`, as `BlockReader` asks, and every touch of
/// those objects is inside `Python::attach`.
pub(crate) struct SegmentReader {
    /// The iterator of the user's segments.
    segments: Py<PyIterator>,
    /// What of numpy the arrays of a segment are read with.
    numpy: Numpy,
    /// The individuals, the lengths of the chromosomes and the extra
    /// fields, the variant fields first, which the core writes into the
    /// `popnei` key of the file.
    header: SourceHeader,
    /// The ploidy of every genotype.
    ploidy: usize,
    /// The chromosomes of the segments read so far.
    chroms: ChromTable,
    /// How many items have been taken, those of no variant among them.
    num_segments: u64,
    /// What ended the pass, which the call raises in the place of the error
    /// of the core: the exception of the user's code or of the interpreter,
    /// or a segment refused here.
    kept: Option<PyPopneiError>,
    /// Whether the iterator is over or failed: every call after gives no
    /// block, as the contract of `BlockReader` asks.
    ended: bool,
}

impl BlockReader for SegmentReader {
    fn next_block(&mut self) -> popnei::Result<Option<Block>> {
        if self.ended {
            return Ok(None);
        }
        match Python::attach(|py| self.next_segment(py)) {
            Ok(Some(block)) => Ok(Some(block)),
            Ok(None) => {
                self.ended = true;
                Ok(None)
            }
            Err(error) => {
                self.ended = true;
                self.kept = Some(error);
                Err(popnei::Error::SegmentsFailed)
            }
        }
    }

    fn individuals(&self) -> &[String] {
        &self.header.individuals
    }

    fn ploidy(&self) -> usize {
        self.ploidy
    }

    fn chroms(&self) -> &ChromTable {
        &self.chroms
    }

    /// Every block holds every column and every field whatever is asked
    /// for: the writer, the one consumer, asks for all of them.
    fn set_needs(&mut self, _needs: Needs) {}

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        Vec::new()
    }

    fn header(&self) -> &SourceHeader {
        &self.header
    }

    fn skip_outside(&mut self, _selection: RegionSelection) -> bool {
        false
    }

    fn num_skipped(&self) -> u64 {
        0
    }
}

/// Where a variant is, for the message of a segment that is refused: the
/// segment and the variant of it, each counted from 1, and its chromosome
/// and its position when those are known and are not what is wrong.
struct Place {
    segment: u64,
    var: Option<u64>,
    at: Option<(String, u64)>,
}

impl Place {
    /// The error of the segment at this place, with what is wrong.
    fn refuses(self, problem: String) -> PyPopneiError {
        PyPopneiError::Segment {
            segment: self.segment,
            var: self.var,
            at: self.at,
            problem,
        }
    }
}

/// The functions and the types of numpy that the arrays of a segment are
/// read with, looked up once for the pass.
struct Numpy {
    /// `numpy.asarray`, which makes an array of Python objects of a list.
    asarray: Py<PyAny>,
    /// `numpy.ascontiguousarray`, which every array of a segment goes
    /// through, and which leaves an array that lies row after row as it is.
    ascontiguousarray: Py<PyAny>,
    /// `numpy.ma.MaskedArray`, which a segment cannot give.
    masked_array: Py<PyAny>,
    /// `numpy.bool_`, the truth value of numpy.
    bool_: Py<PyAny>,
    /// `numpy.integer`, which every integer scalar of numpy is.
    integer: Py<PyAny>,
    /// `numpy.floating`, which every float scalar of numpy is.
    floating: Py<PyAny>,
}

impl Numpy {
    /// What of numpy the pass reads with.
    ///
    /// # Errors
    ///
    /// What importing numpy raised.
    fn new(py: Python<'_>) -> Result<Numpy, PyPopneiError> {
        let numpy = py.import("numpy")?;
        let of = |name: &str| -> PyResult<Py<PyAny>> { Ok(numpy.getattr(name)?.unbind()) };
        Ok(Numpy {
            asarray: of("asarray")?,
            ascontiguousarray: of("ascontiguousarray")?,
            masked_array: numpy.getattr("ma")?.getattr("MaskedArray")?.unbind(),
            bool_: of("bool_")?,
            integer: of("integer")?,
            floating: of("floating")?,
        })
    }
}

/// A number of an array of a segment, widened to a type that holds every
/// value of the numpy types a segment can give, or an integer of Python
/// beyond those.
enum Number {
    /// An integer of a numpy array of integers, or of Python.
    Int(i128),
    /// A float of a numpy array of floats, or of Python.
    Float(f64),
    /// An integer of Python that an `i128` does not hold.
    Big {
        /// The `float64` that holds it exactly, when one does.
        exactly: Option<f64>,
        /// The integer as Python writes it.
        shown: String,
    },
}

/// A numpy type a segment can give its numbers in.
trait ToNumber: Element + Copy {
    fn number(self) -> Number;
}

macro_rules! integers_to_numbers {
    ($($ty:ty),*) => {
        $(impl ToNumber for $ty {
            fn number(self) -> Number {
                Number::Int(i128::from(self))
            }
        })*
    };
}

integers_to_numbers!(i8, i16, i32, i64, u8, u16, u32, u64);

impl ToNumber for f32 {
    fn number(self) -> Number {
        Number::Float(f64::from(self))
    }
}

impl ToNumber for f64 {
    fn number(self) -> Number {
        Number::Float(self)
    }
}

/// What a segment holds once its keys were read: the value of each key it
/// has.
struct SegmentValues<'py> {
    chrom: Bound<'py, PyAny>,
    pos: Bound<'py, PyAny>,
    alleles: Bound<'py, PyAny>,
    gts: Bound<'py, PyAny>,
    id: Option<Bound<'py, PyAny>>,
    qual: Option<Bound<'py, PyAny>>,
    fields: Vec<Bound<'py, PyAny>>,
}

impl SegmentReader {
    /// The block of the next segment that has variants, or `None` when the
    /// iterator is over. A Ctrl-C pending is raised before each item is
    /// asked for.
    ///
    /// # Errors
    ///
    /// The exception of the user's iterator, the `KeyboardInterrupt` of a
    /// Ctrl-C, and the refusal of a segment.
    fn next_segment(&mut self, py: Python<'_>) -> Result<Option<Block>, PyPopneiError> {
        loop {
            py.check_signals()?;
            let mut segments = self.segments.bind(py).clone();
            let Some(item) = segments.next() else {
                return Ok(None);
            };
            let item = item?;
            // A segment is something the user's code computed, and no
            // program yields the 18446744073709551615 that would saturate.
            self.num_segments = self.num_segments.saturating_add(1);
            if let Some(block) = self.block_of(py, &item)? {
                return Ok(Some(block));
            }
        }
    }

    /// The place of the segment being read, with no variant.
    fn segment(&self) -> Place {
        Place {
            segment: self.num_segments,
            var: None,
            at: None,
        }
    }

    /// The place of the variant `var` of the segment being read, counted
    /// from 0, with its chromosome and its position when `at` is true.
    fn variant(&self, var: usize, chrom: &[u32], pos: &[u64], at: bool) -> Place {
        let number = u64::try_from(var).unwrap_or(u64::MAX).saturating_add(1);
        let at = if at {
            chrom
                .get(var)
                .and_then(|number| self.chroms.name(*number))
                .zip(pos.get(var))
                .map(|(name, pos)| (name.to_owned(), *pos))
        } else {
            None
        };
        Place {
            segment: self.num_segments,
            var: Some(number),
            at,
        }
    }

    /// The block of the segment `item`, or `None` for a segment of no
    /// variants.
    ///
    /// # Errors
    ///
    /// The refusal of what "A call that fails" of the writer of segments
    /// lists: a `TypeError` for an item that is not a mapping, and a
    /// `ValueError` for everything else.
    fn block_of(
        &mut self,
        py: Python<'_>,
        item: &Bound<'_, PyAny>,
    ) -> Result<Option<Block>, PyPopneiError> {
        let Ok(segment) = item.cast::<PyMapping>() else {
            return Err(PyPopneiError::SegmentNotAMapping {
                segment: self.num_segments,
                found: item.get_type().name()?.to_string(),
            });
        };
        let values = self.values_of(segment)?;
        let chrom = self.chrom_of(&values.chrom)?;
        let pos = self.pos_of(py, &values.pos, &chrom)?;
        let num_vars = pos.len();
        let chrom = match chrom {
            Chroms::One(number) => vec![number; num_vars],
            Chroms::Each(numbers) => {
                if numbers.len() != num_vars {
                    return Err(self.segment().refuses(format!(
                        "`chrom` holds {} chromosomes and `pos` {num_vars} positions, and \
                         every value of a segment has one for each of its variants",
                        numbers.len()
                    )));
                }
                numbers
            }
        };
        if num_vars == 0 {
            return Ok(None);
        }
        let alleles = self.alleles_of(&values.alleles, &chrom, &pos)?;
        let gts = self.gts_of(py, &values.gts, &alleles, &chrom, &pos)?;
        let id = match &values.id {
            Some(id) => self.ids_of(id, &chrom, &pos)?,
            None => vec![String::new(); num_vars],
        };
        let qual = match &values.qual {
            Some(qual) => self.quals_of(py, qual, &chrom, &pos)?,
            None => vec![f32::NAN; num_vars],
        };
        let mut extra = Vec::with_capacity(self.header.extra_fields.len());
        for (decl, value) in self.header.extra_fields.iter().zip(&values.fields) {
            let values = self.field_of(py, decl, value, &chrom, &pos)?;
            extra.push(ExtraField {
                decl: decl.clone(),
                values,
            });
        }
        Ok(Some(Block {
            num_vars,
            num_individuals: self.header.individuals.len(),
            ploidy: self.ploidy,
            gts,
            chrom: Some(chrom),
            pos: Some(pos),
            id: Some(id),
            alleles: Some(alleles),
            qual: Some(qual),
            passed: None,
            vcf_text: None,
            // A source with no extra field gives blocks without them, as
            // every reader of popnei does.
            extra: (!extra.is_empty()).then_some(extra),
        }))
    }

    /// The value of each key of `segment`.
    ///
    /// # Errors
    ///
    /// A key that is neither one of the six nor a declared field, a key that
    /// is not a text, and a required key that is missing.
    fn values_of<'py>(
        &self,
        segment: &Bound<'py, PyMapping>,
    ) -> Result<SegmentValues<'py>, PyPopneiError> {
        for key in segment.keys()? {
            let known = key.cast::<PyString>().ok().and_then(|key| {
                let key = key.to_str().ok()?;
                let known = KEYS.contains(&key)
                    || self
                        .header
                        .extra_fields
                        .iter()
                        .any(|field| field.name == key);
                Some(known)
            });
            if known != Some(true) {
                let key = match key.cast::<PyString>() {
                    Ok(text) => text.to_string(),
                    Err(_) => key.repr()?.to_string(),
                };
                return Err(self.segment().refuses(format!(
                    "has the key `{key}`, which is neither one of chrom, pos, alleles, gts, id \
                     and qual nor a declared field{}",
                    self.declared()
                )));
            }
        }
        let required = |name: &str| -> Result<Bound<'py, PyAny>, PyPopneiError> {
            self.value_at(segment, name)?.ok_or_else(|| {
                self.segment()
                    .refuses(format!("has no `{name}`, which every segment has"))
            })
        };
        Ok(SegmentValues {
            chrom: required("chrom")?,
            pos: required("pos")?,
            alleles: required("alleles")?,
            gts: required("gts")?,
            // `None` is the key left out, which "Its Python function"
            // of the writer of segments asks for.
            id: self.value_at(segment, "id")?.filter(|id| !id.is_none()),
            qual: self
                .value_at(segment, "qual")?
                .filter(|qual| !qual.is_none()),
            fields: self
                .header
                .extra_fields
                .iter()
                .map(|field| required(&field.name))
                .collect::<Result<_, _>>()?,
        })
    }

    /// The declared fields, as the message of a key that is not one names
    /// them.
    fn declared(&self) -> String {
        if self.header.extra_fields.is_empty() {
            return ", and no field was declared".to_owned();
        }
        let names: Vec<String> = self
            .header
            .extra_fields
            .iter()
            .map(|field| format!("`{}`", field.name))
            .collect();
        format!(", which are {}", names.join(", "))
    }

    /// The value of `name` in `segment`, or `None` when it has none.
    ///
    /// # Errors
    ///
    /// What the mapping raised other than a `KeyError`.
    fn value_at<'py>(
        &self,
        segment: &Bound<'py, PyMapping>,
        name: &str,
    ) -> Result<Option<Bound<'py, PyAny>>, PyPopneiError> {
        match segment.get_item(name) {
            Ok(value) => Ok(Some(value)),
            Err(error) if error.is_instance_of::<PyKeyError>(segment.py()) => Ok(None),
            Err(error) => Err(error.into()),
        }
    }

    /// The chromosome of the segment, one text for all its variants or one
    /// for each, as numbers of the table of the reader.
    ///
    /// # Errors
    ///
    /// A value that is neither a text nor a sequence of texts, an empty
    /// text, and a text that holds a tab, a newline, a carriage return or a
    /// space.
    fn chrom_of(&mut self, chrom: &Bound<'_, PyAny>) -> Result<Chroms, PyPopneiError> {
        if let Ok(name) = chrom.cast::<PyString>() {
            let name = name.to_str()?;
            if let Some(problem) = chrom_problem(name) {
                return Err(self.segment().refuses(problem));
            }
            return Ok(Chroms::One(self.chroms.intern(name)));
        }
        let names = self.texts_of(chrom, "chrom")?;
        let mut numbers = Vec::with_capacity(names.len());
        for (var, name) in names.iter().enumerate() {
            let Ok(name) = name.cast::<PyString>() else {
                return Err(self.variant(var, &[], &[], false).refuses(format!(
                    "its chromosome is {name}, and a chromosome is a text"
                )));
            };
            let name = name.to_str()?;
            if let Some(problem) = chrom_problem(name) {
                return Err(self.variant(var, &[], &[], false).refuses(problem));
            }
            numbers.push(self.chroms.intern(name));
        }
        Ok(Chroms::Each(numbers))
    }

    /// The items of `value`, a sequence of something for each variant.
    ///
    /// # Errors
    ///
    /// A text, which Python would read as a sequence of characters, and
    /// anything that is not a sequence with a length.
    fn texts_of<'py>(
        &self,
        value: &Bound<'py, PyAny>,
        key: &str,
    ) -> Result<Vec<Bound<'py, PyAny>>, PyPopneiError> {
        let not_a_sequence = || {
            self.segment().refuses(format!(
                "`{key}` is {value}, and it is a sequence of one value for each variant"
            ))
        };
        if value.is_instance_of::<PyString>() {
            return Err(not_a_sequence());
        }
        let Ok(length) = value.len() else {
            return Err(not_a_sequence());
        };
        let Ok(items) = value.try_iter() else {
            return Err(not_a_sequence());
        };
        let mut values = Vec::with_capacity(length);
        for item in items {
            values.push(item?);
        }
        Ok(values)
    }

    /// The positions of the segment.
    ///
    /// # Errors
    ///
    /// A value that is not a sequence of integers, and a negative one.
    fn pos_of(
        &self,
        py: Python<'_>,
        pos: &Bound<'_, PyAny>,
        chrom: &Chroms,
    ) -> Result<Vec<u64>, PyPopneiError> {
        let array = self.array_of(py, pos, "pos", Numbers::Integers)?;
        let num_vars = self.shape_fits(&array, "pos", &[None])?;
        let mut positions = Vec::with_capacity(num_vars);
        self.each_number(&array, "pos", &mut |var, number| {
            let position = match number {
                Number::Int(value) => u64::try_from(value).ok(),
                Number::Float(_) | Number::Big { .. } => None,
            };
            let Some(position) = position else {
                let chrom = match chrom {
                    Chroms::One(number) => Some(*number),
                    Chroms::Each(numbers) => numbers.get(var).copied(),
                };
                let on = chrom
                    .and_then(|number| self.chroms.name(number))
                    .map(|name| format!(" on {name}"))
                    .unwrap_or_default();
                return Err(self.variant(var, &[], &[], false).refuses(format!(
                    "its position{on} is {}, and a position is a whole number of 0 or more",
                    shown(&number)
                )));
            };
            positions.push(position);
            Ok(())
        })?;
        Ok(positions)
    }

    /// The alleles of the segment, as a column of the block.
    ///
    /// # Errors
    ///
    /// A value that is not one sequence of texts for each variant, a variant
    /// with no allele or with more than [`MAX_ALLELES_OF_A_VARIANT`], an
    /// allele that is empty, is `.`, is given twice in its variant, or holds
    /// a tab, a newline, a carriage return, a space or a comma.
    fn alleles_of(
        &self,
        alleles: &Bound<'_, PyAny>,
        chrom: &[u32],
        pos: &[u64],
    ) -> Result<AllelesColumn, PyPopneiError> {
        let of_each = self.texts_of(alleles, "alleles")?;
        if of_each.len() != pos.len() {
            return Err(self.segment().refuses(format!(
                "`alleles` holds the alleles of {} variants and `pos` {} positions, and \
                 every value of a segment has one for each of its variants",
                of_each.len(),
                pos.len()
            )));
        }
        let mut column = AllelesColumn::with_num_vars(pos.len()).map_err(|_| {
            self.segment()
                .refuses("holds more alleles than this machine gives the memory for".to_owned())
        })?;
        // The texts of one variant, kept from one variant to the next so
        // that the column is filled with no allocation for each variant.
        let mut texts: Vec<Bound<'_, PyString>> = Vec::new();
        for (var, of_the_var) in of_each.iter().enumerate() {
            let here = || self.variant(var, chrom, pos, true);
            let not_texts = || {
                here().refuses(format!(
                    "its alleles are {of_the_var}, and they are a sequence of texts, the \
                     reference allele first"
                ))
            };
            if of_the_var.is_instance_of::<PyString>() {
                return Err(not_texts());
            }
            let Ok(items) = of_the_var.try_iter() else {
                return Err(not_texts());
            };
            texts.clear();
            for item in items {
                let item = item?;
                let Ok(text) = item.cast_into::<PyString>() else {
                    return Err(not_texts());
                };
                if texts.len() == MAX_ALLELES_OF_A_VARIANT {
                    return Err(here().refuses(format!(
                        "it has more than {MAX_ALLELES_OF_A_VARIANT} alleles, the most that \
                         one byte of `gts` numbers, 0 to 127"
                    )));
                }
                let allele = text.to_str()?;
                if allele.is_empty() {
                    return Err(here().refuses(
                        "one of its alleles is an empty text, and an allele is a text of one \
                         character at least"
                            .to_owned(),
                    ));
                }
                if allele == "." {
                    return Err(here().refuses(
                        "one of its alleles is `.`, which a VCF reads as an allele that is \
                         missing"
                            .to_owned(),
                    ));
                }
                if let Some(problem) = text_problem(allele, ",") {
                    return Err(here().refuses(format!("the allele `{allele}` {problem}")));
                }
                // At most 127 alleles before it, so the comparisons of a
                // variant are bounded and allocate nothing.
                if texts
                    .iter()
                    .any(|before| before.to_str().is_ok_and(|before| before == allele))
                {
                    return Err(here().refuses(format!(
                        "the allele `{allele}` is given twice, and each allele of a variant is \
                         another one"
                    )));
                }
                texts.push(text);
            }
            if texts.is_empty() {
                return Err(
                    here().refuses("it has no allele, and a variant has one at least".to_owned())
                );
            }
            column.push(texts.iter().filter_map(|text| text.to_str().ok()));
        }
        Ok(column)
    }

    /// The genotypes of the segment, as `gts` of the block.
    ///
    /// # Errors
    ///
    /// A value that is not an array of integers of the variants by the
    /// individuals by the ploidy, and an allele below -1 or not among the
    /// alleles of its variant.
    fn gts_of(
        &self,
        py: Python<'_>,
        gts: &Bound<'_, PyAny>,
        alleles: &AllelesColumn,
        chrom: &[u32],
        pos: &[u64],
    ) -> Result<Vec<i8>, PyPopneiError> {
        let array = self.array_of(py, gts, "gts", Numbers::Integers)?;
        let num_individuals = self.header.individuals.len();
        self.shape_fits(
            &array,
            "gts",
            &[Some(pos.len()), Some(num_individuals), Some(self.ploidy)],
        )?;
        // The individuals and the ploidy are each 1 at least, and their
        // product was checked to fit when the shape matched an array that
        // exists.
        let per_var = num_individuals.saturating_mul(self.ploidy).max(1);
        let mut genotypes = Vec::with_capacity(per_var.saturating_mul(pos.len()));
        self.each_number(&array, "gts", &mut |index, number| {
            let var = index.checked_div(per_var).unwrap_or(0);
            let num_alleles = alleles.num_alleles(var);
            let allele = match &number {
                Number::Int(value) => *value,
                Number::Float(_) => {
                    return Err(self.variant(var, chrom, pos, true).refuses(format!(
                        "`gts` holds the float {}, and it holds integers",
                        shown(&number)
                    )));
                }
                Number::Big { shown, .. } if shown.starts_with('-') => i128::MIN,
                Number::Big { .. } => i128::MAX,
            };
            let problem = if allele < -1 {
                format!(
                    "`gts` holds the allele {}, and the alleles of a genotype are numbered \
                     from 0, with -1 for one that is missing",
                    shown(&number)
                )
            } else if let Some(allele) = i8::try_from(allele)
                .ok()
                .filter(|allele| usize::try_from(*allele).map_or(true, |at| at < num_alleles))
            {
                genotypes.push(allele);
                return Ok(());
            } else {
                format!(
                    "`gts` holds the allele {}, which is not among the {num_alleles} alleles \
                     of the variant, numbered from 0",
                    shown(&number)
                )
            };
            Err(self.variant(var, chrom, pos, true).refuses(problem))
        })?;
        Ok(genotypes)
    }

    /// The ids of the segment.
    ///
    /// # Errors
    ///
    /// A value that is not one text or `None` for each variant, and an id
    /// that holds a tab, a newline, a carriage return, a space or a
    /// semicolon.
    fn ids_of(
        &self,
        ids: &Bound<'_, PyAny>,
        chrom: &[u32],
        pos: &[u64],
    ) -> Result<Vec<String>, PyPopneiError> {
        let ids = self.texts_of(ids, "id")?;
        if ids.len() != pos.len() {
            return Err(self.segment().refuses(format!(
                "`id` holds {} ids and `pos` {} positions, and every value of a segment has \
                 one for each of its variants",
                ids.len(),
                pos.len()
            )));
        }
        let mut texts = Vec::with_capacity(ids.len());
        for (var, id) in ids.iter().enumerate() {
            if id.is_none() {
                texts.push(String::new());
                continue;
            }
            let Ok(text) = id.cast::<PyString>() else {
                return Err(self.variant(var, chrom, pos, true).refuses(format!(
                    "its id is {id}, and an id is a text, or None for a variant with none"
                )));
            };
            let text = text.to_str()?;
            if let Some(problem) = text_problem(text, ";") {
                return Err(self
                    .variant(var, chrom, pos, true)
                    .refuses(format!("the id `{text}` {problem}")));
            }
            texts.push(text.to_owned());
        }
        Ok(texts)
    }

    /// The qualities of the segment.
    ///
    /// # Errors
    ///
    /// A value that is not a sequence of numbers, a quality that is an
    /// infinity, which the reader of the file would refuse, and one that a
    /// `float32` does not hold by the rule of a field of floats.
    fn quals_of(
        &self,
        py: Python<'_>,
        qual: &Bound<'_, PyAny>,
        chrom: &[u32],
        pos: &[u64],
    ) -> Result<Vec<f32>, PyPopneiError> {
        let array = self.array_of(py, qual, "qual", Numbers::Any)?;
        self.shape_fits(&array, "qual", &[Some(pos.len())])?;
        let mut quals = Vec::with_capacity(pos.len());
        self.each_number(&array, "qual", &mut |var, number| {
            let quality = match to_f32(&number) {
                Ok(quality) if !quality.is_infinite() => quality,
                Ok(_) | Err(_) => {
                    return Err(self.variant(var, chrom, pos, true).refuses(format!(
                        "its quality is {}, and a quality is a finite number that a float32 \
                         holds, an integer exactly, or NaN for a variant with none",
                        shown_as_given(&number)
                    )));
                }
            };
            quals.push(quality);
            Ok(())
        })?;
        Ok(quals)
    }

    /// The values of the field `decl` in the segment, of its type, by the
    /// rule of "Its Python function" of the writer of segments.
    ///
    /// # Errors
    ///
    /// An array of another shape than the variants, or the variants by the
    /// individuals; floats given to a field of integers; an integer that
    /// the type of the field does not hold, or a field of floats does not
    /// hold exactly; and a finite float beyond the largest of a field of
    /// floats.
    fn field_of(
        &self,
        py: Python<'_>,
        decl: &ExtraFieldDecl,
        value: &Bound<'_, PyAny>,
        chrom: &[u32],
        pos: &[u64],
    ) -> Result<ExtraValues, PyPopneiError> {
        let name = decl.name.as_str();
        let holds = match decl.ty {
            ExtraType::I32 | ExtraType::I64 | ExtraType::U64 => Numbers::Integers,
            ExtraType::F32 | ExtraType::F64 => Numbers::Any,
        };
        let array = self.array_of(py, value, name, holds)?;
        let num_individuals = self.header.individuals.len();
        let per_var = match decl.kind {
            ExtraKind::Variant => {
                self.shape_fits(&array, name, &[Some(pos.len())])?;
                1
            }
            ExtraKind::Genotype => {
                self.shape_fits(&array, name, &[Some(pos.len()), Some(num_individuals)])?;
                num_individuals.max(1)
            }
        };
        let total = per_var.saturating_mul(pos.len());
        let out_of_type = |index: usize, number: &Number, problem: NotHeld| {
            let var = index.checked_div(per_var).unwrap_or(0);
            let ty = decl.ty;
            self.variant(var, chrom, pos, true).refuses(match problem {
                NotHeld::OutOfRange => format!(
                    "`{name}` holds {}, which its type, {ty}, does not hold",
                    shown(number)
                ),
                NotHeld::NotExactly => format!(
                    "`{name}` holds the integer {}, which its type, {ty}, does not hold \
                         exactly: the nearest {ty} is another number",
                    shown(number)
                ),
                NotHeld::Overflows => format!(
                    "`{name}` holds {}, beyond the largest {ty}, which would make it an \
                         infinity",
                    shown(number)
                ),
            })
        };
        Ok(match decl.ty {
            ExtraType::I32 => {
                ExtraValues::I32(self.integers_of_the_type(&array, name, total, &out_of_type)?)
            }
            ExtraType::I64 => {
                ExtraValues::I64(self.integers_of_the_type(&array, name, total, &out_of_type)?)
            }
            ExtraType::U64 => {
                ExtraValues::U64(self.integers_of_the_type(&array, name, total, &out_of_type)?)
            }
            ExtraType::F32 => {
                let mut values = Vec::with_capacity(total);
                self.each_number(&array, name, &mut |index, number| {
                    values.push(
                        to_f32(&number).map_err(|problem| out_of_type(index, &number, problem))?,
                    );
                    Ok(())
                })?;
                ExtraValues::F32(values)
            }
            ExtraType::F64 => {
                let mut values = Vec::with_capacity(total);
                self.each_number(&array, name, &mut |index, number| {
                    values.push(
                        to_f64(&number).map_err(|problem| out_of_type(index, &number, problem))?,
                    );
                    Ok(())
                })?;
                ExtraValues::F64(values)
            }
        })
    }

    /// `value` as a numpy array that lies row after row, in the byte order
    /// of this machine, of a kind that `holds` takes; a float of another
    /// size than 4 or 8 bytes becomes a `float64`.
    ///
    /// A value that is not a numpy array, a list among them, becomes an
    /// array of its Python objects, whose numbers are read one by one: numpy
    /// would make `[2**64 - 1, 5]` an array of `float64`, which does not hold
    /// the first, `[1, 2.5]` one that rounds nothing a field of integers
    /// would refuse, and `[1, True]` one of the integers 1 and 1.
    ///
    /// # Errors
    ///
    /// A masked array, a scalar where an array of the variants is expected,
    /// an array of a kind that `holds` does not take, named by the type the
    /// user gave, and what numpy refuses to make an array of, its
    /// `ValueError` and its `TypeError`. Any other exception raised while
    /// the value is read, one of the user's own `__array__` among them,
    /// goes back as it was raised.
    fn array_of<'py>(
        &self,
        py: Python<'py>,
        value: &Bound<'py, PyAny>,
        key: &str,
        holds: Numbers,
    ) -> Result<Bound<'py, PyUntypedArray>, PyPopneiError> {
        let not_an_array = |error: PyErr| -> PyPopneiError {
            if error.is_instance_of::<PyValueError>(py) || error.is_instance_of::<PyTypeError>(py) {
                self.segment().refuses(format!(
                    "`{key}` is not an array of numbers popnei reads: {error}"
                ))
            } else {
                error.into()
            }
        };
        if value.is_instance(self.numpy.masked_array.bind(py))? {
            return Err(self.segment().refuses(format!(
                "`{key}` is a masked array of numpy, whose mask popnei would lose: \
                 `{key}.filled(...)` with the missing value of your choice gives the array \
                 to write"
            )));
        }
        let array = match value.cast::<PyUntypedArray>() {
            Ok(array) => array.clone(),
            Err(_) => {
                let objects = PyDict::new(py);
                objects.set_item("dtype", "object")?;
                self.numpy
                    .asarray
                    .bind(py)
                    .call((value,), Some(&objects))
                    .map_err(not_an_array)?
                    .cast_into::<PyUntypedArray>()
                    .map_err(|error| not_an_array(error.into()))?
            }
        };
        if array.ndim() == 0 {
            return Err(self.segment().refuses(format!(
                "`{key}` is {}, one value, and it holds one value for each variant",
                value.repr()?
            )));
        }
        let dtype = array.dtype();
        let takes: &[u8] = match holds {
            Numbers::Integers => b"iuO",
            Numbers::Any => b"iufO",
        };
        if !array.is_empty() && !takes.contains(&dtype.kind()) {
            return Err(self.segment().refuses(format!(
                "`{key}` is of the type {dtype}, and it holds {}",
                match holds {
                    Numbers::Integers => "integers",
                    Numbers::Any => "numbers",
                }
            )));
        }
        let array = if dtype.is_native_byteorder() == Some(false) {
            let native = dtype.call_method1("newbyteorder", ("=",))?;
            array
                .call_method1("astype", (native,))?
                .cast_into::<PyUntypedArray>()
                .map_err(|error| not_an_array(error.into()))?
        } else {
            array
        };
        let dtype = array.dtype();
        let array = if dtype.kind() == b'f' && !matches!(dtype.itemsize(), 4 | 8) {
            array
                .call_method1("astype", ("float64",))?
                .cast_into::<PyUntypedArray>()
                .map_err(|error| not_an_array(error.into()))?
        } else {
            array
        };
        self.numpy
            .ascontiguousarray
            .bind(py)
            .call1((array,))
            .map_err(not_an_array)?
            .cast_into::<PyUntypedArray>()
            .map_err(|error| not_an_array(error.into()))
    }

    /// The length of the first dimension of `array`, once its shape is
    /// checked against `shape`, where `None` takes any length.
    ///
    /// # Errors
    ///
    /// Another number of dimensions, or another length of one of them.
    fn shape_fits(
        &self,
        array: &Bound<'_, PyUntypedArray>,
        key: &str,
        shape: &[Option<usize>],
    ) -> Result<usize, PyPopneiError> {
        let found = array.shape();
        let fits = found.len() == shape.len()
            && found
                .iter()
                .zip(shape)
                .all(|(found, expected)| expected.is_none_or(|expected| *found == expected));
        if !fits {
            let expected: Vec<String> = shape
                .iter()
                .map(|length| length.map_or_else(|| "any".to_owned(), |length| length.to_string()))
                .collect();
            return Err(self.segment().refuses(format!(
                "`{key}` is of the shape {found:?}, and it has to be of the shape ({}): \
                 the variants of the segment, given by `pos`, the individuals and the ploidy",
                expected.join(", ")
            )));
        }
        Ok(found.first().copied().unwrap_or(0))
    }

    /// The integers of `array`, the values of the field `key`, `total` of
    /// them, each of the type `T`.
    ///
    /// # Errors
    ///
    /// `out_of_type` of the first that `T` does not hold, with its place
    /// among the values, and what [`SegmentReader::each_number`] refuses.
    fn integers_of_the_type<T: TryFrom<i128>>(
        &self,
        array: &Bound<'_, PyUntypedArray>,
        key: &str,
        total: usize,
        out_of_type: &dyn Fn(usize, &Number, NotHeld) -> PyPopneiError,
    ) -> Result<Vec<T>, PyPopneiError> {
        let mut values = Vec::with_capacity(total);
        self.each_number(array, key, &mut |index, number| {
            let value = match &number {
                Number::Int(value) => T::try_from(*value).ok(),
                Number::Float(_) | Number::Big { .. } => None,
            };
            let Some(value) = value else {
                return Err(out_of_type(index, &number, NotHeld::OutOfRange));
            };
            values.push(value);
            Ok(())
        })?;
        Ok(values)
    }

    /// `each` of every number of `array`, the values of the key `key`, with
    /// its place among them, in the order they lie, which is row after row.
    ///
    /// # Errors
    ///
    /// What `each` gives; of an array of Python objects, a value that is
    /// not a number or is a truth value; and an array of a type that is
    /// none of the ten of numpy that a segment can give its numbers in, nor
    /// of Python objects, which [`SegmentReader::array_of`] refused before.
    fn each_number(
        &self,
        array: &Bound<'_, PyUntypedArray>,
        key: &str,
        each: &mut dyn FnMut(usize, Number) -> Result<(), PyPopneiError>,
    ) -> Result<(), PyPopneiError> {
        macro_rules! of_the_types {
            ($($ty:ty),*) => {
                $(if let Some(done) = numbers_of::<$ty>(array, each) {
                    return done;
                })*
            };
        }
        of_the_types!(i8, i16, i32, i64, u8, u16, u32, u64, f32, f64);
        if let Some(done) = self.numbers_of_objects(array, key, each) {
            return done;
        }
        Err(self.segment().refuses(format!(
            "`{key}` is of the type {}, and popnei reads arrays of integers and floats",
            array.dtype()
        )))
    }

    /// [`SegmentReader::each_number`] for an array of Python objects, the
    /// values of a list, or `None` when `array` is not of objects. An
    /// integer, of Python or of numpy, is read as it is, and a float as a
    /// float, which a key of integers refuses with its value.
    fn numbers_of_objects(
        &self,
        array: &Bound<'_, PyUntypedArray>,
        key: &str,
        each: &mut dyn FnMut(usize, Number) -> Result<(), PyPopneiError>,
    ) -> Option<Result<(), PyPopneiError>> {
        let typed = array.cast::<PyArrayDyn<Py<PyAny>>>().ok()?;
        let py = array.py();
        Some((|| {
            let readonly = typed
                .try_readonly()
                .map_err(|error| PyPopneiError::Broken {
                    message: format!("an array of a segment could not be read: {error}"),
                    path: None,
                })?;
            let values = readonly.as_slice().map_err(|error| PyPopneiError::Broken {
                message: format!("an array of a segment could not be read: {error}"),
                path: None,
            })?;
            for (index, value) in values.iter().enumerate() {
                let value = value.bind(py);
                each(index, self.number_of(value, key)?)?;
            }
            Ok(())
        })())
    }

    /// The number that `value`, a Python object of a list, is.
    ///
    /// # Errors
    ///
    /// A truth value, which Python counts as 1 or 0, and anything that is
    /// not an integer or a float of Python or of numpy.
    fn number_of(&self, value: &Bound<'_, PyAny>, key: &str) -> Result<Number, PyPopneiError> {
        let py = value.py();
        let refused = |what: &str| -> Result<Number, PyPopneiError> {
            Err(self
                .segment()
                .refuses(format!("`{key}` holds {}, {what}", value.repr()?)))
        };
        if value.is_instance_of::<PyBool>() || value.is_instance(self.numpy.bool_.bind(py))? {
            return refused("a truth value, and it holds numbers");
        }
        if value.is_instance_of::<PyInt>() || value.is_instance(self.numpy.integer.bind(py))? {
            if let Ok(integer) = value.extract::<i128>() {
                return Ok(Number::Int(integer));
            }
            // Python compares an integer with a float exactly, so the
            // float is kept only when it is that integer.
            let exactly = match value.extract::<f64>() {
                Ok(float) if value.eq(float)? => Some(float),
                Ok(_) | Err(_) => None,
            };
            return Ok(Number::Big {
                exactly,
                shown: value.str()?.to_string(),
            });
        }
        if value.is_instance_of::<PyFloat>() || value.is_instance(self.numpy.floating.bind(py))? {
            return Ok(Number::Float(value.extract::<f64>()?));
        }
        refused("which is not a number")
    }
}

/// Why a field of a type does not hold a number.
#[derive(Clone, Copy)]
enum NotHeld {
    /// An integer beyond the type of a field of integers, or a float given
    /// to one.
    OutOfRange,
    /// An integer whose nearest float of a field of floats is another
    /// number.
    NotExactly,
    /// A finite float beyond the largest of a field of floats.
    Overflows,
}

/// What an array of a segment holds: integers, the positions, the
/// genotypes and a field of integers, or numbers of any kind, the
/// qualities and a field of floats.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Numbers {
    Integers,
    Any,
}

/// The chromosomes of a segment: one for every variant, given as one text,
/// or one for each.
enum Chroms {
    One(u32),
    Each(Vec<u32>),
}

/// What is wrong with `text` as a chromosome, an id or an allele, which a
/// VCF that `write_vcf` writes from the file would split: a tab, a newline,
/// a carriage return or a space, or one of the characters of `also`.
fn text_problem(text: &str, also: &str) -> Option<String> {
    let found = text.chars().find(|character| {
        matches!(character, '\t' | '\n' | '\r' | ' ') || also.contains(*character)
    })?;
    Some(format!(
        "holds {}, which would break the columns or the lists of a VCF that `write_vcf` \
         writes from the file",
        match found {
            '\t' => "a tab".to_owned(),
            '\n' => "a newline".to_owned(),
            '\r' => "a carriage return".to_owned(),
            ' ' => "a space".to_owned(),
            other => format!("`{other}`"),
        }
    ))
}

/// What is wrong with `name` as a chromosome: an empty text, or what
/// [`text_problem`] finds.
fn chrom_problem(name: &str) -> Option<String> {
    if name.is_empty() {
        return Some("the chromosome is an empty text, and a chromosome has a name".to_owned());
    }
    text_problem(name, "").map(|problem| format!("the chromosome `{name}` {problem}"))
}

/// A number as a message shows it.
fn shown(number: &Number) -> String {
    match number {
        Number::Int(value) => value.to_string(),
        Number::Float(value) => format!("{value:?}"),
        Number::Big { shown, .. } => shown.clone(),
    }
}

/// A number as a message shows it, with `the integer` before an integer, so
/// that a reader tells `16777217` given as an integer from a float.
fn shown_as_given(number: &Number) -> String {
    match number {
        Number::Int(_) | Number::Big { .. } => format!("the integer {}", shown(number)),
        Number::Float(_) => shown(number),
    }
}

/// 2^127, the first `float64` above every `i128`.
const TWO_TO_THE_127: f64 = 170_141_183_460_469_231_731_687_303_715_884_105_728.0;

/// Whether `float`, the nearest float of the integer `integer`, is that
/// integer. A float at 2^127 or above is not, since no `i128` is that
/// large, and the cast back to `i128` would saturate there.
#[expect(
    clippy::cast_possible_truncation,
    reason = "the float is whole, the nearest float of an integer, and below 2^127 in \
              absolute value, so it is an i128 with nothing cut"
)]
fn is_exactly(integer: i128, float: f64) -> bool {
    (-TWO_TO_THE_127..TWO_TO_THE_127).contains(&float) && float as i128 == integer
}

/// `number` as a `float64`, by the rule of a field of floats of "Its
/// Python function" of the writer of segments.
///
/// # Errors
///
/// [`NotHeld::NotExactly`] for an integer that no `float64` is.
#[expect(
    clippy::cast_precision_loss,
    reason = "the nearest float64 of the integer, which is then checked to be it"
)]
fn to_f64(number: &Number) -> Result<f64, NotHeld> {
    match number {
        Number::Int(value) => {
            let float = *value as f64;
            is_exactly(*value, float)
                .then_some(float)
                .ok_or(NotHeld::NotExactly)
        }
        Number::Float(value) => Ok(*value),
        Number::Big { exactly, .. } => exactly.ok_or(NotHeld::NotExactly),
    }
}

/// `number` as a `float32`, by the rule of a field of floats of "Its
/// Python function" of the writer of segments: a float rounded to the
/// nearest `float32`, which is what declaring a field `float32` asks for,
/// and an integer only when a `float32` is it.
///
/// # Errors
///
/// [`NotHeld::NotExactly`] for an integer that no `float32` is, and
/// [`NotHeld::Overflows`] for a finite float beyond the largest `float32`.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    reason = "the rounding to the nearest float32 is what a field declared float32 asks for, \
              the rule of \"Its Python function\" of the writer of segments, and an integer \
              or a float that it does not hold is then refused"
)]
fn to_f32(number: &Number) -> Result<f32, NotHeld> {
    match number {
        Number::Int(value) => {
            let float = *value as f32;
            is_exactly(*value, f64::from(float))
                .then_some(float)
                .ok_or(NotHeld::NotExactly)
        }
        Number::Float(value) => {
            let float = *value as f32;
            if value.is_finite() && float.is_infinite() {
                return Err(NotHeld::Overflows);
            }
            Ok(float)
        }
        Number::Big { exactly, .. } => {
            let float = exactly.ok_or(NotHeld::NotExactly)? as f32;
            (exactly == &Some(f64::from(float)))
                .then_some(float)
                .ok_or(NotHeld::NotExactly)
        }
    }
}

/// [`SegmentReader::each_number`] for an array of the type `T`, or `None`
/// when `array` is
/// of another.
fn numbers_of<T: ToNumber>(
    array: &Bound<'_, PyUntypedArray>,
    each: &mut dyn FnMut(usize, Number) -> Result<(), PyPopneiError>,
) -> Option<Result<(), PyPopneiError>> {
    let typed = array.cast::<PyArrayDyn<T>>().ok()?;
    Some((|| {
        let readonly = typed
            .try_readonly()
            .map_err(|error| PyPopneiError::Broken {
                message: format!("an array of a segment could not be read: {error}"),
                path: None,
            })?;
        // `numpy.ascontiguousarray` made it lie row after row, which is
        // asked first because `as_slice` gives the slice of an array that
        // lies column after column too.
        if !typed.is_c_contiguous() {
            return Err(PyPopneiError::Broken {
                message: "an array of a segment does not lie row after row after \
                          `numpy.ascontiguousarray`"
                    .to_owned(),
                path: None,
            });
        }
        let values = readonly.as_slice().map_err(|error| PyPopneiError::Broken {
            message: format!("an array of a segment could not be read: {error}"),
            path: None,
        })?;
        for (index, value) in values.iter().enumerate() {
            each(index, value.number())?;
        }
        Ok(())
    })())
}
