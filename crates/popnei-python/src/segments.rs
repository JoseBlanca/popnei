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
use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyIterator, PyMapping, PyString};

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
        .map(|(chrom, length)| Ok((chrom.clone(), length_of(chrom, length)?)))
        .collect::<Result<Vec<(String, u64)>, PyPopneiError>>()?;
    let num_vars_per_block = num_vars_per_block
        .map(|asked_for| count_of("num_vars_per_block", asked_for))
        .transpose()?;
    // An object that cannot be iterated is the `TypeError` of Python,
    // before the file is made.
    let segments = segments.try_iter()?.unbind();
    let ascontiguousarray = py.import("numpy")?.getattr("ascontiguousarray")?.unbind();
    let reader = SegmentReader {
        segments,
        ascontiguousarray,
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
/// each of its own name.
///
/// # Errors
///
/// [`PyPopneiError::Argument`] for no individual and for a name given twice.
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
    /// `numpy.ascontiguousarray`, which every array of a segment goes
    /// through: it turns a list into an array, and leaves an array that
    /// lies row after row as it is.
    ascontiguousarray: Py<PyAny>,
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

/// A number of a numpy array of one of the integer types or one of the
/// float ones, widened to a type that holds every value of all of them.
#[derive(Clone, Copy)]
enum Number {
    Int(i128),
    Float(f64),
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
                return Err(self.segment().refuses(format!(
                    "has the key {key}, which is neither one of chrom, pos, alleles, gts, id \
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
            id: self.value_at(segment, "id")?,
            qual: self.value_at(segment, "qual")?,
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
    /// A value that is neither a text nor a sequence of texts, and a text
    /// that holds a tab, a newline, a carriage return or a space.
    fn chrom_of(&mut self, chrom: &Bound<'_, PyAny>) -> Result<Chroms, PyPopneiError> {
        if let Ok(name) = chrom.cast::<PyString>() {
            let name = name.to_str()?;
            if let Some(problem) = text_problem(name, "") {
                return Err(self
                    .segment()
                    .refuses(format!("the chromosome `{name}` {problem}")));
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
            if let Some(problem) = text_problem(name, "") {
                return Err(self
                    .variant(var, &[], &[], false)
                    .refuses(format!("the chromosome `{name}` {problem}")));
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
        self.integers_only(&array, "pos")?;
        let mut positions = Vec::with_capacity(num_vars);
        each_number(&array, &mut |var, number| {
            let position = match number {
                Number::Int(value) => u64::try_from(value).ok(),
                Number::Float(_) => None,
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
                    shown(number)
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
    /// with no allele or with more than [`MAX_ALLELES_OF_A_VARIANT`], and an
    /// allele that is empty or holds a tab, a newline, a carriage return, a
    /// space or a comma.
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
                if let Some(problem) = text_problem(allele, ",") {
                    return Err(here().refuses(format!("the allele `{allele}` {problem}")));
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
        self.integers_only(&array, "gts")?;
        // The individuals and the ploidy are each 1 at least, and their
        // product was checked to fit when the shape matched an array that
        // exists.
        let per_var = num_individuals.saturating_mul(self.ploidy).max(1);
        let mut genotypes = Vec::with_capacity(per_var.saturating_mul(pos.len()));
        each_number(&array, &mut |index, number| {
            let var = index.checked_div(per_var).unwrap_or(0);
            let num_alleles = alleles.num_alleles(var);
            let allele = match number {
                Number::Int(value) => value,
                Number::Float(_) => i128::MIN,
            };
            let problem = if allele < -1 {
                format!(
                    "`gts` holds the allele {}, and the alleles of a genotype are numbered \
                     from 0, with -1 for one that is missing",
                    shown(number)
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
                    shown(number)
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
    /// A value that is not a sequence of numbers, and a quality that is an
    /// infinity, or that a `float32` cannot hold, which the reader of the
    /// file would refuse.
    fn quals_of(
        &self,
        py: Python<'_>,
        qual: &Bound<'_, PyAny>,
        chrom: &[u32],
        pos: &[u64],
    ) -> Result<Vec<f32>, PyPopneiError> {
        let array = self.array_of(py, qual, "qual", Numbers::Any)?;
        self.shape_fits(&array, "qual", &[Some(pos.len())])?;
        self.numbers_only(&array, "qual")?;
        let mut quals = Vec::with_capacity(pos.len());
        each_number(&array, &mut |var, number| {
            let quality = to_f32(number);
            if quality.is_infinite() {
                return Err(self.variant(var, chrom, pos, true).refuses(format!(
                    "its quality is {}, and a quality is a finite number that a float32 \
                     holds, or NaN for a variant with none",
                    shown(number)
                )));
            }
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
    /// individuals; floats given to a field of integers; and an integer that
    /// the type of the field does not hold.
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
        match holds {
            Numbers::Integers => self.integers_only(&array, name)?,
            Numbers::Any => self.numbers_only(&array, name)?,
        }
        let out_of_type = |index: usize, number: Number| {
            let var = index.checked_div(per_var).unwrap_or(0);
            self.variant(var, chrom, pos, true).refuses(format!(
                "`{name}` holds {}, which its type, {}, does not hold",
                shown(number),
                decl.ty
            ))
        };
        Ok(match decl.ty {
            ExtraType::I32 => ExtraValues::I32(integers_of_the_type(&array, total, &out_of_type)?),
            ExtraType::I64 => ExtraValues::I64(integers_of_the_type(&array, total, &out_of_type)?),
            ExtraType::U64 => ExtraValues::U64(integers_of_the_type(&array, total, &out_of_type)?),
            ExtraType::F32 => {
                let mut values = Vec::with_capacity(total);
                each_number(&array, &mut |_, number| {
                    values.push(to_f32(number));
                    Ok(())
                })?;
                ExtraValues::F32(values)
            }
            ExtraType::F64 => {
                let mut values = Vec::with_capacity(total);
                each_number(&array, &mut |_, number| {
                    values.push(to_f64(number));
                    Ok(())
                })?;
                ExtraValues::F64(values)
            }
        })
    }

    /// `value` as a numpy array that lies row after row, with the type it
    /// has, or the one numpy gives a list; a float of another size than 4
    /// or 8 bytes becomes a `float64`.
    ///
    /// When it is to hold `Numbers::Integers` and `value` is not an array, a list that
    /// numpy would make an array of floats is made an array of its Python
    /// objects instead, whose numbers are read one by one: numpy makes
    /// `[2**64 - 1, 5]` an array of `float64`, which does not hold the
    /// first, and `[1, 2.5]` one that rounds nothing a field of integers
    /// would refuse.
    ///
    /// # Errors
    ///
    /// What numpy refuses to make an array of, a list of rows of other
    /// lengths among them.
    fn array_of<'py>(
        &self,
        py: Python<'py>,
        value: &Bound<'py, PyAny>,
        key: &str,
        holds: Numbers,
    ) -> Result<Bound<'py, PyUntypedArray>, PyPopneiError> {
        let not_an_array = |problem: String| {
            self.segment().refuses(format!(
                "`{key}` is not an array of numbers popnei reads: {problem}"
            ))
        };
        let ascontiguousarray = self.ascontiguousarray.bind(py);
        let array = ascontiguousarray
            .call1((value,))
            .map_err(|error| not_an_array(error.to_string()))?
            .cast_into::<PyUntypedArray>()
            .map_err(|error| not_an_array(error.to_string()))?;
        let array = if holds == Numbers::Integers
            && !value.is_instance_of::<PyUntypedArray>()
            && !array.is_empty()
            && !matches!(array.dtype().kind(), b'i' | b'u')
        {
            let objects = PyDict::new(py);
            objects.set_item("dtype", "object")?;
            ascontiguousarray
                .call((value,), Some(&objects))
                .map_err(|error| not_an_array(error.to_string()))?
                .cast_into::<PyUntypedArray>()
                .map_err(|error| not_an_array(error.to_string()))?
        } else {
            array
        };
        let dtype = array.dtype();
        if dtype.kind() == b'f' && !matches!(dtype.itemsize(), 4 | 8) {
            let widened = array.call_method1("astype", ("float64",))?;
            return widened
                .cast_into::<PyUntypedArray>()
                .map_err(|error| not_an_array(error.to_string()));
        }
        Ok(array)
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

    /// That `array` holds integers, or Python objects whose numbers are
    /// checked one by one as they are read. An array of no value holds no
    /// wrong one, whatever its type: numpy makes `[]` an array of
    /// `float64`, and a segment of no variants is written as `"pos": []`.
    ///
    /// # Errors
    ///
    /// An array of another type, floats among them even when they are
    /// whole, since a float in a field of integers is a mistake numpy would
    /// round in silence.
    fn integers_only(
        &self,
        array: &Bound<'_, PyUntypedArray>,
        key: &str,
    ) -> Result<(), PyPopneiError> {
        if array.is_empty() || matches!(array.dtype().kind(), b'i' | b'u' | b'O') {
            return Ok(());
        }
        Err(self.segment().refuses(format!(
            "`{key}` is of the type {}, and it holds integers",
            array.dtype()
        )))
    }

    /// That `array` holds integers or floats, or no value.
    ///
    /// # Errors
    ///
    /// An array of another type.
    fn numbers_only(
        &self,
        array: &Bound<'_, PyUntypedArray>,
        key: &str,
    ) -> Result<(), PyPopneiError> {
        if array.is_empty() || matches!(array.dtype().kind(), b'i' | b'u' | b'f') {
            return Ok(());
        }
        Err(self.segment().refuses(format!(
            "`{key}` is of the type {}, and it holds numbers",
            array.dtype()
        )))
    }
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

/// A number as a message shows it.
fn shown(number: Number) -> String {
    match number {
        Number::Int(value) => value.to_string(),
        Number::Float(value) => format!("{value:?}"),
    }
}

/// `number` as a `float64`, rounded to the nearest one when it is an
/// integer above 2^53, which is what declaring a field `float64` asks for.
fn to_f64(number: Number) -> f64 {
    match number {
        Number::Int(value) => value as f64,
        Number::Float(value) => value,
    }
}

/// `number` as a `float32`, rounded to the nearest one, which is what
/// declaring a field `float32` asks for. A value beyond the largest
/// `float32` becomes an infinity.
#[expect(
    clippy::cast_possible_truncation,
    reason = "the rounding to the nearest float32 is what a field declared float32 asks for, \
              the rule of \"Its Python function\" of the writer of segments"
)]
fn to_f32(number: Number) -> f32 {
    match number {
        Number::Int(value) => value as f32,
        Number::Float(value) => value as f32,
    }
}

/// The integers of `array`, `total` of them, each of the type `T`.
///
/// # Errors
///
/// `out_of_type` of the first that `T` does not hold, with its place among
/// the values.
fn integers_of_the_type<T: TryFrom<i128>>(
    array: &Bound<'_, PyUntypedArray>,
    total: usize,
    out_of_type: &dyn Fn(usize, Number) -> PyPopneiError,
) -> Result<Vec<T>, PyPopneiError> {
    let mut values = Vec::with_capacity(total);
    each_number(array, &mut |index, number| {
        let value = match number {
            Number::Int(value) => T::try_from(value).ok(),
            Number::Float(_) => None,
        };
        let Some(value) = value else {
            return Err(out_of_type(index, number));
        };
        values.push(value);
        Ok(())
    })?;
    Ok(values)
}

/// `each` of every number of `array`, with its place among them, in the
/// order they lie, which is row after row.
///
/// # Errors
///
/// What `each` gives, and an array of a type that is none of the ten of
/// numpy that a segment can give its numbers in, nor of Python objects;
/// the callers checked the
/// kind of the array before, so that is a byte order other than the one of
/// this machine.
fn each_number(
    array: &Bound<'_, PyUntypedArray>,
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
    if let Some(done) = numbers_of_objects(array, each) {
        return done;
    }
    Err(PyPopneiError::Argument {
        name: "segments",
        problem: format!(
            "holds an array of the type {}, and popnei reads arrays of integers and floats \
             in the byte order of this machine",
            array.dtype()
        ),
    })
}

/// [`each_number`] for an array of the type `T`, or `None` when `array` is
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

/// [`each_number`] for an array of Python objects, the integers of a list
/// that numpy would have made floats of, or `None` when `array` is not of
/// objects. An integer of Python is read as it is, a float as a float,
/// which a field of integers refuses with its value.
fn numbers_of_objects(
    array: &Bound<'_, PyUntypedArray>,
    each: &mut dyn FnMut(usize, Number) -> Result<(), PyPopneiError>,
) -> Option<Result<(), PyPopneiError>> {
    let typed = array.cast::<PyArrayDyn<Py<PyAny>>>().ok()?;
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
            let value = value.bind(array.py());
            let number = if let Ok(integer) = value.extract::<i128>() {
                Number::Int(integer)
            } else if let Ok(float) = value.extract::<f64>() {
                Number::Float(float)
            } else {
                return Err(PyPopneiError::Argument {
                    name: "segments",
                    problem: format!("holds {value}, where a number was expected"),
                });
            };
            each(index, number)?;
        }
        Ok(())
    })())
}
