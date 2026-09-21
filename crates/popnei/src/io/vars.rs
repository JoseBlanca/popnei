//! The vars file: popnei's own file of variants, and what it says about
//! itself.
//!
//! A vars file is one arrow IPC file, also called feather v2, with one
//! record batch for each block of variants and the columns `chrom`, `pos`,
//! `id`, `alleles`, `qual` and `gts`, so that any program with an arrow
//! library opens it as a table. What the file says about itself travels in
//! two keys whose values are json: `popnei`, in the schema, holds what is
//! known before the first variant, and `popnei_batches`, in the footer,
//! what is known only after the last one, one entry for each batch.
//! [`VarsMetadata`] is the value of the first and [`BatchInfo`] one entry
//! of the second.
//!
//! [`write_vars`] writes the variants of any reader of blocks into such a
//! file, one batch for each block, with its buffers compressed with lz4,
//! and [`VarsWriter`] is what it does it with, for a caller that has the
//! blocks and not a reader.
//!
//! `docs/specs/io_vars.md` has the format, the writer and the reader. What
//! is here is the two keys and the writer; the reader is being written.

use std::collections::HashMap;
use std::io::Write;
use std::sync::Arc;

use arrow_array::builder::{ListBuilder, StringBuilder};
use arrow_array::{
    ArrayRef, FixedSizeListArray, Float32Array, Int8Array, ListArray, RecordBatch, StringArray,
    UInt64Array,
};
use arrow_buffer::{NullBuffer, ScalarBuffer};
use arrow_ipc::CompressionType;
use arrow_ipc::writer::{FileWriter, IpcWriteOptions};
use arrow_schema::{ArrowError, DataType, Field, Schema, SchemaRef};
use serde_json::{Map, Value};

use crate::block::{AllelesColumn, Block, BlockReader, BlockSize, Reblock, size_of_the_blocks};
use crate::error::{Error, Result};
use crate::variant::{ChromTable, Needs};

/// The key of the schema of a vars file whose value says what is known
/// before its first variant.
const POPNEI_KEY: &str = "popnei";

/// The key of the footer of a vars file whose value says, for each batch,
/// how many variants it holds and where they are. It is in the footer and
/// not in the schema because arrow-rs writes the schema before the first
/// batch and the footer when the file is finished.
const POPNEI_BATCHES_KEY: &str = "popnei_batches";

/// The version of the format that popnei writes, `major.minor`.
///
/// A reader refuses a file whose major version, the part before the dot, is
/// not [`FORMAT_VERSION_READ`], and reads a file with any minor version, a
/// later one than its own too, ignoring the keys and the columns it does
/// not know: that is what lets a later version of the format add a column
/// without making the files or the readers that are there useless.
pub const FORMAT_VERSION: &str = "1.0";

/// The major version of the format that popnei reads, the part of
/// [`FORMAT_VERSION`] before the dot.
pub const FORMAT_VERSION_READ: &str = "1";

/// What a vars file says about itself before its first variant, the value
/// of the `popnei` key of its schema.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VarsMetadata {
    /// The whole string, "1.0". Only the part before the dot is checked.
    pub format_version: String,
    /// The names of the individuals, in the order of their genotypes in
    /// every row of `gts`.
    pub individuals: Vec<String>,
    /// How many alleles the genotype of one individual holds.
    pub ploidy: usize,
    /// How many variants a batch of the file holds, the last one aside.
    pub num_vars_per_block: usize,
}

/// What the footer of a vars file says of one of its batches, one entry of
/// the `popnei_batches` key.
///
/// It is a batch and not a block: the blocks that `iter_blocks` gives from a
/// file are cut at the size the caller asks for, which does not have to be
/// that of its batches.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BatchInfo {
    /// How many variants the batch holds, which the footer of an arrow file
    /// does not say: with it the variants of a file are counted without
    /// reading a batch.
    pub num_vars: usize,
    /// One for each chromosome with a variant in the batch, in the order in
    /// which they first appear. Empty in a file with no chrom and pos
    /// columns.
    pub regions: Vec<Region>,
}

/// The variants of one chromosome in one batch of a vars file, as the
/// smallest and the largest of their positions.
///
/// A function that is asked for a region skips, without decompressing them,
/// the batches with no entry that overlaps it, and that is right for a file
/// whose variants are in any order: the two positions are the smallest and
/// the largest of the batch and not those of its first and its last
/// variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    /// The name of the chromosome, as the `chrom` column holds it.
    pub chrom: String,
    /// The smallest and the largest position of the variants of `chrom` in
    /// the batch, both included.
    pub min_pos: u64,
    /// The largest one.
    pub max_pos: u64,
}

/// The value of the `popnei` key of the schema, as the json that goes into
/// the file.
pub(crate) fn metadata_as_json(metadata: &VarsMetadata) -> String {
    let individuals: Vec<String> = metadata
        .individuals
        .iter()
        .map(|name| json_text(name))
        .collect();
    json_object(&[
        ("format_version", json_text(&metadata.format_version)),
        ("individuals", json_array(&individuals)),
        ("ploidy", metadata.ploidy.to_string()),
        (
            "num_vars_per_block",
            metadata.num_vars_per_block.to_string(),
        ),
    ])
}

/// The value of the `popnei_batches` key of the footer, one entry for each
/// batch in the order of the batches.
pub(crate) fn batches_as_json(batches: &[BatchInfo]) -> String {
    let entries: Vec<String> = batches.iter().map(batch_as_json).collect();
    json_array(&entries)
}

/// What the `popnei` key of the schema of a vars file says, from its value.
///
/// # Errors
///
/// The source is not a vars file when that value is not a json object or
/// when one of its four keys is missing or does not hold what that key
/// holds.
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "the reader of the vars file is its caller and is being written; the tests of \
                  this module call it already, so this holds for the build without them, and \
                  the lint itself asks for it to go when the reader lands"
    )
)]
pub(crate) fn metadata_from_json(text: &str) -> Result<VarsMetadata> {
    let value = json_of_the_popnei_key(text)?;
    let Some(object) = value.as_object() else {
        return Err(not_a_vars_file(format!(
            "the value of the `{POPNEI_KEY}` key of its schema is not a json object"
        )));
    };
    Ok(VarsMetadata {
        format_version: format_version_of(object)?,
        individuals: individuals_of(object)?,
        ploidy: count_of(object, "ploidy")?,
        num_vars_per_block: count_of(object, "num_vars_per_block")?,
    })
}

/// What the `popnei_batches` key of the footer of a vars file says of each
/// of its batches, from its value.
///
/// # Errors
///
/// The source is not a vars file when that value is not a json array of
/// entries, or when an entry does not hold the number of variants of its
/// batch or holds a region that is not one.
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "the reader of the vars file is its caller and is being written; the tests of \
                  this module call it already, so this holds for the build without them, and \
                  the lint itself asks for it to go when the reader lands"
    )
)]
pub(crate) fn batches_from_json(text: &str) -> Result<Vec<BatchInfo>> {
    let value: Value = serde_json::from_str(text).map_err(|problem| {
        not_a_vars_file(format!(
            "the value of the `{POPNEI_BATCHES_KEY}` key of its footer is not json: {problem}"
        ))
    })?;
    let Some(entries) = value.as_array() else {
        return Err(not_a_vars_file(format!(
            "the value of the `{POPNEI_BATCHES_KEY}` key of its footer is not a json array of one entry for each batch"
        )));
    };
    entries.iter().map(batch_from_json).collect()
}

/// One entry of the `popnei_batches` key. A batch with no regions, which a
/// file without the `chrom` and the `pos` columns has, is written with
/// `num_vars` alone: there is nothing to say about where its variants are.
fn batch_as_json(batch: &BatchInfo) -> String {
    let num_vars = ("num_vars", batch.num_vars.to_string());
    if batch.regions.is_empty() {
        return json_object(&[num_vars]);
    }
    let regions: Vec<String> = batch.regions.iter().map(region_as_json).collect();
    json_object(&[num_vars, ("regions", json_array(&regions))])
}

/// One region of one entry of the `popnei_batches` key.
fn region_as_json(region: &Region) -> String {
    json_object(&[
        ("chrom", json_text(&region.chrom)),
        ("min_pos", region.min_pos.to_string()),
        ("max_pos", region.max_pos.to_string()),
    ])
}

/// One json object, the keys in the order they are given and their values
/// already json.
///
/// The keys of a `serde_json::Map` come out in the order of the alphabet,
/// and the two keys of a vars file are written in the order in which
/// `docs/specs/io_vars.md` gives them, which is the order a person who
/// opens the file with another program reads them in.
fn json_object(entries: &[(&str, String)]) -> String {
    let written: Vec<String> = entries
        .iter()
        .map(|(key, value)| format!("{key}:{value}", key = json_text(key)))
        .collect();
    format!("{{{joined}}}", joined = written.join(","))
}

/// One json array of values that are already json.
fn json_array(values: &[String]) -> String {
    format!("[{joined}]", joined = values.join(","))
}

/// One text as json, with the quotes and the escapes that json asks for,
/// which the name of an individual or of a chromosome can need.
fn json_text(text: &str) -> String {
    Value::String(text.to_owned()).to_string()
}

/// The value of the `popnei` key, parsed as json.
fn json_of_the_popnei_key(text: &str) -> Result<Value> {
    serde_json::from_str(text).map_err(|problem| {
        not_a_vars_file(format!(
            "the value of the `{POPNEI_KEY}` key of its schema is not json: {problem}"
        ))
    })
}

/// The version of the format the file gives, whatever it is: which versions
/// popnei reads is not read here but where the file is opened.
fn format_version_of(object: &Map<String, Value>) -> Result<String> {
    let Some(value) = object.get("format_version") else {
        return Err(missing_in_the_popnei_key("format_version"));
    };
    match value.as_str() {
        Some(version) => Ok(version.to_owned()),
        None => Err(not_a_vars_file(format!(
            "`format_version` of the `{POPNEI_KEY}` key of its schema is {value} and not a text"
        ))),
    }
}

/// The names of the individuals, in the order the file gives them. That two
/// of them are the same name is not read here but where the file is opened.
fn individuals_of(object: &Map<String, Value>) -> Result<Vec<String>> {
    let Some(value) = object.get("individuals") else {
        return Err(missing_in_the_popnei_key("individuals"));
    };
    let Some(names) = value.as_array() else {
        return Err(not_a_vars_file(format!(
            "`individuals` of the `{POPNEI_KEY}` key of its schema is {value} and not a json array of names"
        )));
    };
    names
        .iter()
        .map(|name| match name.as_str() {
            Some(name) => Ok(name.to_owned()),
            None => Err(not_a_vars_file(format!(
                "the `individuals` of the `{POPNEI_KEY}` key of its schema hold {name}, which is not the name of an individual"
            ))),
        })
        .collect()
}

/// A value of the `popnei` key that says how many of something there are,
/// the ploidy or the variants of a batch: a whole number of 1 or more that
/// this machine can count.
fn count_of(object: &Map<String, Value>, key: &str) -> Result<usize> {
    let Some(value) = object.get(key) else {
        return Err(missing_in_the_popnei_key(key));
    };
    let count = value
        .as_u64()
        .and_then(|number| usize::try_from(number).ok())
        .filter(|count| *count >= 1);
    count.ok_or_else(|| {
        not_a_vars_file(format!(
            "`{key}` of the `{POPNEI_KEY}` key of its schema is {value}, and it says how many of something there are: a whole number of 1 or more that this machine can count"
        ))
    })
}

/// One entry of the `popnei_batches` key, what the footer says of one
/// batch.
fn batch_from_json(value: &Value) -> Result<BatchInfo> {
    let Some(object) = value.as_object() else {
        return Err(not_a_vars_file(format!(
            "an entry of the `{POPNEI_BATCHES_KEY}` key of its footer is {value} and not a json object"
        )));
    };
    let num_vars = num_vars_of(object)?;
    let regions = match object.get("regions") {
        // A batch of a file without the `chrom` and the `pos` columns has
        // nothing to say about where its variants are, and its entry holds
        // the number of them alone.
        None => Vec::new(),
        Some(value) => {
            let Some(regions) = value.as_array() else {
                return Err(not_a_vars_file(format!(
                    "the `regions` of an entry of the `{POPNEI_BATCHES_KEY}` key of its footer are {value} and not a json array"
                )));
            };
            regions
                .iter()
                .map(region_from_json)
                .collect::<Result<_>>()?
        }
    };
    Ok(BatchInfo { num_vars, regions })
}

/// How many variants one entry of the `popnei_batches` key says its batch
/// holds: a whole number of 0 or more, since a batch of no variants is one
/// that another arrow program can write.
fn num_vars_of(object: &Map<String, Value>) -> Result<usize> {
    let Some(value) = object.get("num_vars") else {
        return Err(not_a_vars_file(format!(
            "an entry of the `{POPNEI_BATCHES_KEY}` key of its footer has no `num_vars`"
        )));
    };
    let num_vars = value
        .as_u64()
        .and_then(|number| usize::try_from(number).ok());
    num_vars.ok_or_else(|| {
        not_a_vars_file(format!(
            "`num_vars` of an entry of the `{POPNEI_BATCHES_KEY}` key of its footer is {value}, and it says how many variants a batch holds: a whole number of 0 or more that this machine can count"
        ))
    })
}

/// One region of one entry of the `popnei_batches` key.
fn region_from_json(value: &Value) -> Result<Region> {
    let Some(object) = value.as_object() else {
        return Err(not_a_vars_file(format!(
            "a region of an entry of the `{POPNEI_BATCHES_KEY}` key of its footer is {value} and not a json object"
        )));
    };
    let Some(chrom) = object.get("chrom") else {
        return Err(missing_in_a_region("chrom"));
    };
    let Some(chrom) = chrom.as_str() else {
        return Err(not_a_vars_file(format!(
            "the `chrom` of a region of an entry of the `{POPNEI_BATCHES_KEY}` key of its footer is {chrom}, which is not the name of a chromosome"
        )));
    };
    Ok(Region {
        chrom: chrom.to_owned(),
        min_pos: position_of_a_region(object, "min_pos")?,
        max_pos: position_of_a_region(object, "max_pos")?,
    })
}

/// One of the two positions of a region: a whole number that a position of
/// a variant fits in.
fn position_of_a_region(object: &Map<String, Value>, key: &str) -> Result<u64> {
    let Some(value) = object.get(key) else {
        return Err(missing_in_a_region(key));
    };
    value.as_u64().ok_or_else(|| {
        not_a_vars_file(format!(
            "the `{key}` of a region of an entry of the `{POPNEI_BATCHES_KEY}` key of its footer is {value}, which is not the position of a variant"
        ))
    })
}

/// One of the four values of the `popnei` key that the file does not have.
fn missing_in_the_popnei_key(key: &str) -> Error {
    not_a_vars_file(format!(
        "the value of the `{POPNEI_KEY}` key of its schema has no `{key}`"
    ))
}

/// One of the three values of a region that the file does not have.
fn missing_in_a_region(key: &str) -> Error {
    not_a_vars_file(format!(
        "a region of an entry of the `{POPNEI_BATCHES_KEY}` key of its footer has no `{key}`"
    ))
}

/// The error of a source that is not a vars file, whatever of one it lacks.
fn not_a_vars_file(problem: String) -> Error {
    Error::NotAVarsFile { problem }
}

// The name of each column of a vars file, in the order in which a `Block`
// holds them, which is the order of the columns of the file.
const CHROM_COLUMN: &str = "chrom";
const POS_COLUMN: &str = "pos";
const ID_COLUMN: &str = "id";
const ALLELES_COLUMN: &str = "alleles";
const QUAL_COLUMN: &str = "qual";
const GTS_COLUMN: &str = "gts";

/// The name arrow gives the values of a list, which is what pyarrow, and
/// so pandas and polars, write and read.
const ITEM_FIELD: &str = "item";

/// What the writer has to write on: the sink itself while no block has
/// arrived, since the columns of the file are those of the first block and
/// an arrow file starts with its schema, and the arrow writer once a block
/// has fixed them.
enum Sink<W: Write> {
    /// No block has been written, and nothing of the file either.
    BeforeTheFirstBlock(W),
    /// The arrow writer over the sink. It is boxed because it is far
    /// larger than a sink, and every writer would otherwise be of its
    /// size.
    Started(Box<FileWriter<W>>),
    /// The header of the file could not be written and the sink went with
    /// it, so there is nothing left to write on.
    Gone,
}

/// The writer of a vars file: it takes the blocks of a reader and writes
/// each of them as one record batch of an arrow IPC file.
///
/// [`VarsWriter::new`] takes what the `popnei` key of the schema says, and
/// the columns of the file are those of the first block written, which its
/// [`Block::fields`] gives: a source that gives the genotypes alone makes
/// a file with a `gts` column only. Every batch of an arrow file shares one
/// schema, so a later block of other columns is an error.
///
/// The memory it uses is one block, and the vector of genotypes of each
/// block becomes the buffer of its `gts` column with no copy.
///
/// [`write_vars`] is what a caller with a reader uses; this one is for a
/// caller that has the blocks.
pub struct VarsWriter<W: Write> {
    sink: Sink<W>,
    /// What the `popnei` key of the schema says, which is written before
    /// the first batch.
    metadata: VarsMetadata,
    /// How many alleles one variant holds, the individuals times the
    /// ploidy, which is the width of the `gts` column.
    alleles_per_var: i32,
    /// How the buffers of every batch are compressed, lz4.
    options: IpcWriteOptions,
    /// What the first block written fixed, or `None` while no block has
    /// been written.
    columns: Option<WrittenColumns>,
    /// One entry for each batch written, in their order, which
    /// [`VarsWriter::finish`] writes into the footer.
    batches: Vec<BatchInfo>,
}

/// What the first block written fixed: the columns of the file and the
/// arrow schema that every one of its batches has.
struct WrittenColumns {
    fields: Needs,
    schema: SchemaRef,
}

impl<W: Write> VarsWriter<W> {
    /// The writer of a vars file of those individuals on `sink`.
    ///
    /// `num_vars_per_block` is what the `popnei` key will say, the size the
    /// caller gives the blocks it writes; it is not checked against them,
    /// and the last block of a file is shorter than the others.
    ///
    /// # Errors
    ///
    /// When the genotypes of one variant, the individuals times the ploidy,
    /// are more than the `gts` column of an arrow file holds, 2147483647
    /// alleles, which is a block of more memory than a machine gives. And
    /// when arrow-rs does not take lz4, which no build of popnei reaches,
    /// since the crate takes `arrow-ipc` with its `lz4` feature on.
    pub fn new(
        sink: W,
        individuals: &[String],
        ploidy: usize,
        num_vars_per_block: usize,
    ) -> Result<VarsWriter<W>> {
        let num_individuals = individuals.len();
        let alleles_per_var = num_individuals
            .checked_mul(ploidy)
            .and_then(|alleles| i32::try_from(alleles).ok())
            // A block of one variant of so many individuals is already more
            // memory than a machine gives, 2 GB of genotypes, so the error
            // is that of a block that does not fit, of the size this writer
            // was given, which is a size its caller asked for.
            .ok_or(Error::BlockTooLarge {
                num_vars_per_block,
                num_individuals,
                ploidy,
                size: BlockSize::AskedFor,
            })?;
        let options = IpcWriteOptions::default()
            .try_with_compression(Some(CompressionType::LZ4_FRAME))
            .map_err(not_written)?;
        Ok(VarsWriter {
            sink: Sink::BeforeTheFirstBlock(sink),
            metadata: VarsMetadata {
                format_version: FORMAT_VERSION.to_owned(),
                individuals: individuals.to_vec(),
                ploidy,
                num_vars_per_block,
            },
            alleles_per_var,
            options,
            columns: None,
            batches: Vec::new(),
        })
    }

    /// It writes `block` as one batch, of whatever size the block has.
    ///
    /// `chroms` is the table of the reader the block came from, which holds
    /// the names behind its chromosome numbers: the file holds the name of
    /// the chromosome of every variant as text.
    ///
    /// # Errors
    ///
    /// When the block holds other individuals or another ploidy than the
    /// writer was built for, when it does not pass [`Block::check`], when
    /// its columns are not those of the first block written, and when
    /// `chroms` has no name for one of its chromosome numbers; in each of
    /// them nothing of the block is written. And when the sink fails.
    pub fn write_block(&mut self, block: Block, chroms: &ChromTable) -> Result<()> {
        self.fits(&block)?;
        // The rows of the batch are read out of the arrays by their place
        // in them, so the arrays have to be of the size the block says.
        block.check()?;
        let fields = block.fields();
        let schema = match &self.columns {
            Some(written) if written.fields != fields => {
                return Err(Error::VarsBlockColumns {
                    first: written.fields,
                    found: fields,
                });
            }
            Some(written) => Arc::clone(&written.schema),
            None => Arc::new(self.schema_of(fields)),
        };
        let num_vars = block.num_vars;
        let (arrays, regions) = self.arrays_of(block, chroms, fields)?;
        let batch = RecordBatch::try_new(Arc::clone(&schema), arrays).map_err(not_written)?;
        self.writer(&schema)?.write(&batch).map_err(not_written)?;
        self.columns = Some(WrittenColumns { fields, schema });
        self.batches.push(BatchInfo { num_vars, regions });
        Ok(())
    }

    /// It writes the footer, with the `popnei_batches` key, and gives the
    /// sink back.
    ///
    /// A writer that was given no block writes a file with the two keys, a
    /// `gts` column and no batch, which reads back as no variants.
    ///
    /// # Errors
    ///
    /// When the sink fails.
    pub fn finish(self) -> Result<W> {
        let batches = batches_as_json(&self.batches);
        let mut writer = match self.sink {
            Sink::Started(writer) => *writer,
            // A source with no variants: the columns of the file are the
            // one column every file has.
            Sink::BeforeTheFirstBlock(sink) => {
                let schema = schema_of(Needs::GTS, self.alleles_per_var, &self.metadata);
                FileWriter::try_new_with_options(sink, &schema, self.options)
                    .map_err(not_written)?
            }
            Sink::Gone => return Err(the_sink_is_gone()),
        };
        writer.write_metadata(POPNEI_BATCHES_KEY, batches);
        writer.into_inner().map_err(not_written)
    }

    /// That the block holds the individuals and the ploidy the `popnei` key
    /// of the file names, which every batch of it holds.
    fn fits(&self, block: &Block) -> Result<()> {
        let num_individuals = self.metadata.individuals.len();
        let ploidy = self.metadata.ploidy;
        if block.num_individuals == num_individuals && block.ploidy == ploidy {
            return Ok(());
        }
        Err(Error::VarsBlockDoesNotFit {
            num_individuals,
            ploidy,
            found_num_individuals: block.num_individuals,
            found_ploidy: block.ploidy,
        })
    }

    /// The schema of a file whose blocks hold `fields`, with the `popnei`
    /// key of this writer.
    fn schema_of(&self, fields: Needs) -> Schema {
        schema_of(fields, self.alleles_per_var, &self.metadata)
    }

    /// The columns of one batch, in the order of the columns of the file,
    /// and where the variants of that batch are.
    ///
    /// The block goes in by value: the vector of its genotypes becomes the
    /// buffer of the `gts` column with no copy.
    fn arrays_of(
        &self,
        block: Block,
        chroms: &ChromTable,
        fields: Needs,
    ) -> Result<(Vec<ArrayRef>, Vec<Region>)> {
        let Block {
            num_vars: _,
            num_individuals: _,
            ploidy: _,
            gts,
            chrom,
            pos,
            id,
            alleles,
            qual,
        } = block;
        let mut arrays: Vec<ArrayRef> = Vec::new();
        let mut regions = Vec::new();
        // `Block::fields` reports the chromosome and the position only when
        // the block holds both columns.
        if let (Some(chrom), Some(pos)) = (chrom, pos) {
            let (names, where_they_are) = chrom_column(&chrom, &pos, chroms)?;
            regions = where_they_are;
            arrays.push(Arc::new(names));
            arrays.push(Arc::new(UInt64Array::new(ScalarBuffer::from(pos), None)));
        }
        if let Some(id) = id {
            arrays.push(Arc::new(id_column(&id)));
        }
        if let Some(alleles) = alleles {
            arrays.push(Arc::new(alleles_column(&alleles)));
        }
        if let Some(qual) = qual {
            arrays.push(Arc::new(qual_column(qual)));
        }
        // A block of no variants holds the genotypes although its `gts` is
        // empty, and one built without them does not, so which it is is
        // read from the fields and not from the vector.
        if fields.contains(Needs::GTS) {
            arrays.push(gts_column(gts, self.alleles_per_var)?);
        }
        Ok((arrays, regions))
    }

    /// The arrow writer, made over the sink with the schema of the first
    /// block when the first block is written.
    ///
    /// # Errors
    ///
    /// When the header of the file cannot be written, which leaves the
    /// writer with no sink.
    fn writer(&mut self, schema: &SchemaRef) -> Result<&mut FileWriter<W>> {
        match std::mem::replace(&mut self.sink, Sink::Gone) {
            Sink::BeforeTheFirstBlock(sink) => {
                let writer = FileWriter::try_new_with_options(sink, schema, self.options.clone())
                    .map_err(not_written)?;
                self.sink = Sink::Started(Box::new(writer));
            }
            Sink::Started(writer) => self.sink = Sink::Started(writer),
            Sink::Gone => {}
        }
        match &mut self.sink {
            Sink::Started(writer) => Ok(writer),
            Sink::BeforeTheFirstBlock(_) | Sink::Gone => Err(the_sink_is_gone()),
        }
    }
}

/// Every variant of `reader` into a vars file on `sink`, one batch for each
/// block, and the sink back.
///
/// It asks `reader` for every field, so a file written from a VCF holds its
/// six columns and can stand in for it in any later analysis, and it puts a
/// [`Reblock`] of `num_vars_per_block` variants over it, `None` for
/// [`default_num_vars_per_block`](crate::block::default_num_vars_per_block)
/// for the individuals of `reader`, which is then the number that the
/// `popnei` key says. A source with no variants gives a file with the two
/// keys, a `gts` column and no batch.
///
/// This is what both binding crates call. The Python binding crate opens
/// the file, refuses a path that exists, and removes the file when this
/// returns an error.
///
/// # Errors
///
/// When `num_vars_per_block` is 0 or the blocks of that size do not fit in
/// the memory of the machine, when the reader fails, when a block of it is
/// not one a vars file can hold, and when the sink fails. The bytes that
/// were written before the error are on the sink: a file that a failed call
/// was writing is not one that can be read, and it is the caller that
/// removes it.
pub fn write_vars<R: BlockReader, W: Write>(
    mut reader: R,
    sink: W,
    num_vars_per_block: Option<usize>,
) -> Result<W> {
    // Every field, so that a file written from a VCF holds its six columns
    // whether or not the user will read them.
    reader.set_needs(Needs::ALL);
    let individuals = reader.individuals().to_vec();
    let ploidy = reader.ploidy();
    // The size the `popnei` key says is the one the blocks are cut to,
    // which for a caller that asked for none is the one popnei chose for
    // these individuals.
    let (num_vars_per_block, _) =
        size_of_the_blocks(num_vars_per_block, individuals.len(), ploidy)?;
    let mut writer = VarsWriter::new(sink, &individuals, ploidy, num_vars_per_block)?;
    let mut blocks = Reblock::new(reader, Some(num_vars_per_block))?;
    while let Some(block) = blocks.next_block()? {
        writer.write_block(block, blocks.chroms())?;
    }
    writer.finish()
}

/// The schema of a vars file whose blocks hold `fields`, with the `popnei`
/// key of its metadata: the columns of the table of "What it holds" of
/// `docs/specs/io_vars.md` that those fields fill, in the order it gives
/// them.
fn schema_of(fields: Needs, alleles_per_var: i32, metadata: &VarsMetadata) -> Schema {
    let mut columns = Vec::new();
    if fields.contains(Needs::CHROM_POS) {
        columns.push(Field::new(CHROM_COLUMN, DataType::Utf8, false));
        columns.push(Field::new(POS_COLUMN, DataType::UInt64, false));
    }
    if fields.contains(Needs::ID) {
        columns.push(Field::new(ID_COLUMN, DataType::Utf8, true));
    }
    if fields.contains(Needs::ALLELES) {
        columns.push(Field::new(ALLELES_COLUMN, alleles_type(), false));
    }
    if fields.contains(Needs::QUAL) {
        columns.push(Field::new(QUAL_COLUMN, DataType::Float32, true));
    }
    if fields.contains(Needs::GTS) {
        columns.push(Field::new(GTS_COLUMN, gts_type(alleles_per_var), false));
    }
    let popnei = HashMap::from([(POPNEI_KEY.to_owned(), metadata_as_json(metadata))]);
    Schema::new(columns).with_metadata(popnei)
}

/// The arrow type of the `alleles` column, one list of texts for each
/// variant. The values of the list are named and take nulls as pyarrow
/// writes them, so that a file of popnei and a file of pyarrow have the
/// same column; no allele of a variant is a null.
fn alleles_type() -> DataType {
    DataType::List(Arc::new(Field::new(ITEM_FIELD, DataType::Utf8, true)))
}

/// The arrow type of the `gts` column, the `alleles_per_var` alleles of
/// each variant. A fixed size list keeps no offsets, so the column is one
/// flat buffer of variants x individuals x ploidy signed bytes.
fn gts_type(alleles_per_var: i32) -> DataType {
    DataType::FixedSizeList(
        Arc::new(Field::new(ITEM_FIELD, DataType::Int8, true)),
        alleles_per_var,
    )
}

/// The `chrom` column of one batch, the name of the chromosome of every
/// variant, and where the variants of each of its chromosomes are: the
/// smallest and the largest of their positions, in the order in which the
/// chromosomes first appear.
///
/// The two positions are the smallest and the largest and not those of the
/// first and the last variant, so a caller that skips the batches outside a
/// region skips none that has a variant in it, whether the file is sorted
/// or not.
///
/// # Errors
///
/// When `chroms` has no name for a number of the block.
fn chrom_column(
    chrom: &[u32],
    pos: &[u64],
    chroms: &ChromTable,
) -> Result<(StringArray, Vec<Region>)> {
    let mut names = StringBuilder::new();
    let mut regions: Vec<Region> = Vec::new();
    // Which region the variant before went into. The variants of a block
    // run along one chromosome, so the entry of a name is searched for
    // only where the number changes.
    let mut last: Option<(u32, usize)> = None;
    for (number, position) in chrom.iter().copied().zip(pos.iter().copied()) {
        let Some(name) = chroms.name(number) else {
            return Err(Error::VarsChromNameMissing { number });
        };
        names.append_value(name);
        let at = match last {
            Some((seen, at)) if seen == number => at,
            _ => region_of(&mut regions, name, position),
        };
        if let Some(region) = regions.get_mut(at) {
            region.min_pos = region.min_pos.min(position);
            region.max_pos = region.max_pos.max(position);
        }
        last = Some((number, at));
    }
    Ok((names.finish(), regions))
}

/// Where the region of `chrom` is among the ones found so far, which is one
/// more entry when that chromosome has no variant in the batch yet.
fn region_of(regions: &mut Vec<Region>, chrom: &str, position: u64) -> usize {
    if let Some(at) = regions.iter().position(|region| region.chrom == chrom) {
        return at;
    }
    let at = regions.len();
    regions.push(Region {
        chrom: chrom.to_owned(),
        min_pos: position,
        max_pos: position,
    });
    at
}

/// The `id` column of one batch. A block holds an empty id for a variant
/// that has none and the file holds a null, which is what any other program
/// that opens it takes for a value that is not there.
fn id_column(ids: &[String]) -> StringArray {
    let mut column = StringBuilder::new();
    for id in ids {
        match id.is_empty() {
            true => column.append_null(),
            false => column.append_value(id),
        }
    }
    column.finish()
}

/// The `alleles` column of one batch, the reference allele of each variant
/// first and then its alternative ones.
fn alleles_column(alleles: &AllelesColumn) -> ListArray {
    let mut column = ListBuilder::new(StringBuilder::new());
    for var in 0..alleles.num_vars() {
        for allele in 0..alleles.num_alleles(var) {
            column.values().append_value(alleles.allele(var, allele));
        }
        column.append(true);
    }
    column.finish()
}

/// The `qual` column of one batch. A block holds a NaN for a variant with
/// no quality and the file holds a null, as it does for an id that is not
/// there.
fn qual_column(qual: Vec<f32>) -> Float32Array {
    let there: NullBuffer = qual.iter().map(|quality| !quality.is_nan()).collect();
    Float32Array::new(ScalarBuffer::from(qual), Some(there))
}

/// The genotypes of one block as the `gts` column of one batch: one flat
/// buffer of the alleles of every variant, `alleles_per_var` of them in
/// each.
///
/// The vector of the block becomes that buffer with no copy, which is what
/// keeps the memory of the writer to one block.
///
/// # Errors
///
/// When the alleles are not `alleles_per_var` for each variant, which
/// [`Block::check`] is what says before the block reaches here.
fn gts_column(gts: Vec<i8>, alleles_per_var: i32) -> Result<ArrayRef> {
    let alleles = Int8Array::new(ScalarBuffer::from(gts), None);
    let column = FixedSizeListArray::try_new(
        Arc::new(Field::new(ITEM_FIELD, DataType::Int8, true)),
        alleles_per_var,
        Arc::new(alleles),
        None,
    )
    .map_err(not_written)?;
    Ok(Arc::new(column))
}

/// What arrow-rs said while the vars file was being written, as the error
/// of the crate.
///
/// What fails while a file is written is the output, so the
/// `std::io::Error` of the system is kept as it is, with the number it
/// carries: that number is what a binding crate builds the exception of its
/// language with.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "arrow-rs has twenty cases of its error and popnei tells one of them, the error of \
              the output, from every other"
)]
fn not_written(problem: ArrowError) -> Error {
    match problem {
        ArrowError::IoError(_, failure) => Error::Io(failure),
        other => Error::Io(std::io::Error::other(other)),
    }
}

/// The error of a writer whose sink is gone, which is what is left after
/// the header of the file could not be written.
fn the_sink_is_gone() -> Error {
    Error::Io(std::io::Error::other(
        "the header of the vars file could not be written, and the writer has nothing left to \
         write on",
    ))
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::sync::Arc;

    use arrow_array::cast::AsArray;
    use arrow_array::types::{Float32Type, Int8Type, UInt64Type};
    use arrow_array::{Array, FixedSizeListArray, RecordBatch};
    use arrow_ipc::reader::FileReader;
    use arrow_schema::{DataType, Field};

    use super::{
        BatchInfo, FORMAT_VERSION, FORMAT_VERSION_READ, POPNEI_BATCHES_KEY, POPNEI_KEY, Region,
        VarsMetadata, VarsWriter, batches_as_json, batches_from_json, gts_column, metadata_as_json,
        metadata_from_json, write_vars,
    };
    use crate::block::{AllelesColumn, Block, BlockReader};
    use crate::error::{Error, Result};
    use crate::variant::{ChromTable, MISSING_ALLELE, Needs};

    /// One row of the table of `cases.vcf` of "How it is verified" of
    /// `docs/specs/io_vcf.md`, which the blocks of these tests are built
    /// from: three individuals of the ploidy 2.
    struct Row {
        chrom: &'static str,
        pos: u64,
        /// Empty for a variant with no id, which the file holds as a null.
        id: &'static str,
        alleles: &'static [&'static str],
        /// `None` for a variant whose QUAL is `.`, which is a NaN in the
        /// column of a block and a null in the file.
        qual: Option<f32>,
        gts: &'static [i8],
    }

    /// How many individuals the blocks of these tests have and how many
    /// alleles the genotype of each holds, which are those of `cases.vcf`.
    const CASES_INDIVIDUALS: usize = 3;
    const CASES_PLOIDY: usize = 2;

    /// The four bytes that every lz4 frame starts with, which arrow writes
    /// before each buffer of a batch it compressed with lz4.
    const LZ4_FRAME_MARK: [u8; 4] = [0x04, 0x22, 0x4d, 0x18];

    const MISSING: i8 = MISSING_ALLELE;

    /// The four variants of `cases.vcf`, in the order of the file.
    const CASES: [Row; 4] = [
        Row {
            chrom: "chr1",
            pos: 100,
            id: "rs1",
            alleles: &["A", "T"],
            qual: Some(29.5),
            gts: &[0, 0, 0, 1, 1, 1],
        },
        Row {
            chrom: "chr1",
            pos: 200,
            id: "",
            alleles: &["A", "T"],
            qual: None,
            gts: &[MISSING, MISSING, 0, 1, MISSING, 0],
        },
        Row {
            chrom: "chr1",
            pos: 300,
            id: "",
            alleles: &["A", "G", "T"],
            qual: Some(67.0),
            gts: &[1, 2, 2, 1, 2, 2],
        },
        Row {
            chrom: "chr1",
            pos: 400,
            id: "",
            alleles: &["T"],
            qual: Some(47.0),
            gts: &[0, 0, 0, 0, 0, 0],
        },
    ];

    /// Four variants of two chromosomes in no order, the second cargo test
    /// of "How it is verified" of the writer: what the footer says of a
    /// batch of them is the smallest and the largest position of each
    /// chromosome and not the positions of its first and its last variant.
    const NOT_SORTED: [Row; 4] = [
        Row {
            chrom: "chr1",
            pos: 300,
            id: "",
            alleles: &["A", "T"],
            qual: None,
            gts: &[0, 0, 0, 1, 1, 1],
        },
        Row {
            chrom: "chr2",
            pos: 50,
            id: "",
            alleles: &["A", "T"],
            qual: None,
            gts: &[0, 0, 0, 1, 1, 1],
        },
        Row {
            chrom: "chr1",
            pos: 100,
            id: "",
            alleles: &["A", "T"],
            qual: None,
            gts: &[0, 0, 0, 1, 1, 1],
        },
        Row {
            chrom: "chr2",
            pos: 60,
            id: "",
            alleles: &["A", "T"],
            qual: None,
            gts: &[0, 0, 0, 1, 1, 1],
        },
    ];

    /// A block built by hand from the rows of one of those tables, with
    /// every column, which is what a reader of such a VCF gives. The names
    /// of the chromosomes are interned into `chroms` in the order in which
    /// they first appear, as a reader does.
    fn block_of(rows: &[&Row], chroms: &mut ChromTable) -> Block {
        let mut gts = Vec::new();
        let mut chrom = Vec::new();
        let mut pos = Vec::new();
        let mut id = Vec::new();
        let mut qual = Vec::new();
        let mut alleles = AllelesColumn::with_num_vars(rows.len()).expect("the alleles");
        for row in rows {
            gts.extend_from_slice(row.gts);
            chrom.push(chroms.intern(row.chrom));
            pos.push(row.pos);
            id.push(row.id.to_owned());
            qual.push(row.qual.unwrap_or(f32::NAN));
            let texts: Vec<String> = row.alleles.iter().map(|text| (*text).to_owned()).collect();
            alleles.push(&texts);
        }
        Block {
            num_vars: rows.len(),
            num_individuals: CASES_INDIVIDUALS,
            ploidy: CASES_PLOIDY,
            gts,
            chrom: Some(chrom),
            pos: Some(pos),
            id: Some(id),
            alleles: Some(alleles),
            qual: Some(qual),
        }
    }

    /// A block of the four variants of `cases.vcf`, with every column.
    fn cases_block(chroms: &mut ChromTable) -> Block {
        let rows: Vec<&Row> = CASES.iter().collect();
        block_of(&rows, chroms)
    }

    /// The three individuals of `cases.vcf`.
    fn cases_individuals() -> Vec<String> {
        vec!["ind1".to_owned(), "ind2".to_owned(), "ind3".to_owned()]
    }

    /// A reader of blocks written for these tests: it gives the blocks it
    /// was built with, with the table of chromosome names they were built
    /// from, so that a test writes a vars file from variants it holds as
    /// literals and not from a file.
    struct GivenBlocks {
        individuals: Vec<String>,
        ploidy: usize,
        chroms: ChromTable,
        /// The blocks still to give, the last one first.
        left: Vec<Block>,
        /// What it was last asked to fill.
        needs: Needs,
    }

    impl GivenBlocks {
        /// A reader of the three individuals of `cases.vcf` that gives
        /// `blocks`, in their order, with the names of `chroms`.
        fn of(blocks: Vec<Block>, chroms: ChromTable) -> GivenBlocks {
            let mut left = blocks;
            left.reverse();
            GivenBlocks {
                individuals: cases_individuals(),
                ploidy: CASES_PLOIDY,
                chroms,
                left,
                needs: Needs::ALL,
            }
        }
    }

    impl BlockReader for GivenBlocks {
        fn next_block(&mut self) -> Result<Option<Block>> {
            Ok(self.left.pop())
        }

        fn individuals(&self) -> &[String] {
            &self.individuals
        }

        fn ploidy(&self) -> usize {
            self.ploidy
        }

        fn chroms(&self) -> &ChromTable {
            &self.chroms
        }

        fn set_needs(&mut self, needs: Needs) {
            self.needs = needs;
        }
    }

    /// A vars file opened as any other program with an arrow library opens
    /// it: its schema, its batches and the key of its footer.
    struct FileRead {
        /// The name of each column, its arrow type and whether it takes
        /// nulls, in the order of the columns of the file.
        columns: Vec<(String, DataType, bool)>,
        /// What the `popnei` key of the schema says.
        metadata: VarsMetadata,
        /// What the `popnei_batches` key of the footer says.
        batches: Vec<BatchInfo>,
        rows: Vec<RecordBatch>,
    }

    /// The bytes of a vars file, read back with the `FileReader` of
    /// arrow-rs, which is what stands for another program here until the
    /// reader of popnei is written.
    fn file_read(bytes: Vec<u8>) -> FileRead {
        let read =
            FileReader::try_new(Cursor::new(bytes), None).expect("the bytes are an arrow file");
        let schema = read.schema();
        let columns = schema
            .fields()
            .iter()
            .map(|field| {
                (
                    field.name().clone(),
                    field.data_type().clone(),
                    field.is_nullable(),
                )
            })
            .collect();
        let metadata = metadata_from_json(
            schema
                .metadata()
                .get(POPNEI_KEY)
                .expect("the file has the `popnei` key"),
        )
        .expect("the `popnei` key says what the file holds");
        let batches = batches_from_json(
            read.custom_metadata()
                .get(POPNEI_BATCHES_KEY)
                .expect("the file has the `popnei_batches` key"),
        )
        .expect("the `popnei_batches` key says what each batch holds");
        let rows = read
            .map(|batch| batch.expect("a batch of the file"))
            .collect();
        FileRead {
            columns,
            metadata,
            batches,
            rows,
        }
    }

    /// The arrow type of the `alleles` column, a list of texts.
    fn alleles_type() -> DataType {
        DataType::List(Arc::new(Field::new_list_field(DataType::Utf8, true)))
    }

    /// The arrow type of the `gts` column of a file of that many alleles
    /// for each variant.
    fn gts_type(alleles_per_var: i32) -> DataType {
        DataType::FixedSizeList(
            Arc::new(Field::new_list_field(DataType::Int8, true)),
            alleles_per_var,
        )
    }

    /// The names of the chromosomes of a batch, one for each variant.
    fn chroms_of(batch: &RecordBatch) -> Vec<String> {
        batch
            .column_by_name("chrom")
            .expect("the `chrom` column")
            .as_string::<i32>()
            .iter()
            .map(|name| name.expect("a chromosome").to_owned())
            .collect()
    }

    /// The positions of a batch.
    fn positions_of(batch: &RecordBatch) -> Vec<Option<u64>> {
        batch
            .column_by_name("pos")
            .expect("the `pos` column")
            .as_primitive::<UInt64Type>()
            .iter()
            .collect()
    }

    /// The ids of a batch, `None` where the file holds a null.
    fn ids_of(batch: &RecordBatch) -> Vec<Option<String>> {
        batch
            .column_by_name("id")
            .expect("the `id` column")
            .as_string::<i32>()
            .iter()
            .map(|id| id.map(str::to_owned))
            .collect()
    }

    /// The alleles of each variant of a batch.
    fn alleles_of(batch: &RecordBatch) -> Vec<Vec<String>> {
        batch
            .column_by_name("alleles")
            .expect("the `alleles` column")
            .as_list::<i32>()
            .iter()
            .map(|alleles| {
                alleles
                    .expect("the alleles of a variant")
                    .as_string::<i32>()
                    .iter()
                    .map(|allele| allele.expect("an allele").to_owned())
                    .collect()
            })
            .collect()
    }

    /// The qualities of a batch, `None` where the file holds a null.
    fn quals_of(batch: &RecordBatch) -> Vec<Option<f32>> {
        batch
            .column_by_name("qual")
            .expect("the `qual` column")
            .as_primitive::<Float32Type>()
            .iter()
            .collect()
    }

    /// The genotypes of a batch, variant after variant.
    fn gts_of(batch: &RecordBatch) -> Vec<i8> {
        batch
            .column_by_name("gts")
            .expect("the `gts` column")
            .as_any()
            .downcast_ref::<FixedSizeListArray>()
            .expect("the `gts` column is a fixed size list")
            .values()
            .as_primitive::<Int8Type>()
            .values()
            .to_vec()
    }

    /// How many variants each batch of a file holds.
    fn num_rows_of(batches: &[RecordBatch]) -> Vec<usize> {
        batches.iter().map(RecordBatch::num_rows).collect()
    }

    /// The vars file of the four variants of `cases.vcf`, written from a
    /// reader that gives them in one block, with batches of
    /// `num_vars_per_block` variants.
    fn cases_written_in_batches_of(num_vars_per_block: usize) -> Vec<u8> {
        let mut chroms = ChromTable::new();
        let block = cases_block(&mut chroms);
        let reader = GivenBlocks::of(vec![block], chroms);
        write_vars(reader, Vec::new(), Some(num_vars_per_block)).expect("the file was written")
    }

    /// The four variants of the table of `cases.vcf`, in batches of three:
    /// every column of the file is the one of the table of "What it holds"
    /// of `docs/specs/io_vars.md`, with the type and the nulls it gives.
    #[test]
    fn the_four_variants_of_cases_vcf_are_written_as_two_batches_with_every_column() {
        let file = file_read(cases_written_in_batches_of(3));

        assert_eq!(
            file.columns,
            vec![
                ("chrom".to_owned(), DataType::Utf8, false),
                ("pos".to_owned(), DataType::UInt64, false),
                ("id".to_owned(), DataType::Utf8, true),
                ("alleles".to_owned(), alleles_type(), false),
                ("qual".to_owned(), DataType::Float32, true),
                // 3 individuals of the ploidy 2.
                ("gts".to_owned(), gts_type(6), false),
            ]
        );
        assert_eq!(num_rows_of(&file.rows), [3, 1]);

        let first = &file.rows[0];
        assert_eq!(chroms_of(first), ["chr1", "chr1", "chr1"]);
        assert_eq!(positions_of(first), [Some(100), Some(200), Some(300)]);
        // The variants of the table with no id, which the file holds as a
        // null: the last three of the four.
        assert_eq!(ids_of(first), [Some("rs1".to_owned()), None, None]);
        assert_eq!(
            alleles_of(first),
            [vec!["A", "T"], vec!["A", "T"], vec!["A", "G", "T"]]
        );
        // The variant whose QUAL is a dot is a null too.
        assert_eq!(quals_of(first), [Some(29.5), None, Some(67.0)]);
        assert_eq!(
            gts_of(first),
            [
                0, 0, 0, 1, 1, 1, MISSING, MISSING, 0, 1, MISSING, 0, 1, 2, 2, 1, 2, 2
            ]
        );

        let last = &file.rows[1];
        assert_eq!(chroms_of(last), ["chr1"]);
        assert_eq!(positions_of(last), [Some(400)]);
        assert_eq!(ids_of(last), [None]);
        assert_eq!(alleles_of(last), [vec!["T"]]);
        assert_eq!(quals_of(last), [Some(47.0)]);
        assert_eq!(gts_of(last), [0, 0, 0, 0, 0, 0]);

        // The nulls of the whole file, which is what a program that opens
        // it counts: three ids and one quality of the four variants.
        let nulls = |name: &str| -> usize {
            file.rows
                .iter()
                .map(|batch| batch.column_by_name(name).expect("the column").null_count())
                .sum()
        };
        assert_eq!(nulls("id"), 3);
        assert_eq!(nulls("qual"), 1);
        assert_eq!(nulls("chrom"), 0);
        assert_eq!(nulls("pos"), 0);
        assert_eq!(nulls("alleles"), 0);
        assert_eq!(nulls("gts"), 0);
    }

    /// The two keys of the file say what it holds before its first variant
    /// and where the variants of each batch are, which is what a reader
    /// knows as soon as the file is opened.
    #[test]
    fn the_two_keys_of_the_file_say_what_it_holds_and_where_the_variants_of_each_batch_are() {
        let file = file_read(cases_written_in_batches_of(3));

        assert_eq!(
            file.metadata,
            VarsMetadata {
                format_version: FORMAT_VERSION.to_owned(),
                individuals: cases_individuals(),
                ploidy: 2,
                num_vars_per_block: 3,
            }
        );
        assert_eq!(
            file.batches,
            vec![
                BatchInfo {
                    num_vars: 3,
                    regions: vec![Region {
                        chrom: "chr1".to_owned(),
                        min_pos: 100,
                        max_pos: 300,
                    }],
                },
                BatchInfo {
                    num_vars: 1,
                    regions: vec![Region {
                        chrom: "chr1".to_owned(),
                        min_pos: 400,
                        max_pos: 400,
                    }],
                },
            ]
        );
    }

    /// The variants of a file are in no order in general, and the entry of
    /// a batch holds the smallest and the largest position of each of its
    /// chromosomes, in the order in which they first appear, so that a
    /// caller that skips the batches outside a region skips none that has a
    /// variant in it.
    #[test]
    fn the_regions_of_a_batch_of_variants_that_are_not_sorted_are_their_smallest_and_largest() {
        let mut chroms = ChromTable::new();
        let rows: Vec<&Row> = NOT_SORTED.iter().collect();
        let block = block_of(&rows, &mut chroms);
        let reader = GivenBlocks::of(vec![block], chroms);
        let bytes = write_vars(reader, Vec::new(), Some(4)).expect("the file was written");
        let file = file_read(bytes);

        assert_eq!(num_rows_of(&file.rows), [4]);
        assert_eq!(
            chroms_of(&file.rows[0]),
            ["chr1", "chr2", "chr1", "chr2"],
            "the chromosome of every variant is written as its name"
        );
        assert_eq!(
            file.batches,
            vec![BatchInfo {
                num_vars: 4,
                regions: vec![
                    Region {
                        chrom: "chr1".to_owned(),
                        min_pos: 100,
                        max_pos: 300,
                    },
                    Region {
                        chrom: "chr2".to_owned(),
                        min_pos: 50,
                        max_pos: 60,
                    },
                ],
            }]
        );
    }

    /// A file holds the columns its source could fill: a source whose
    /// blocks carry the genotypes alone gives a file of one column, and its
    /// batches have nothing to say about where their variants are.
    #[test]
    fn a_source_of_the_genotypes_alone_gives_a_file_of_one_column_and_batches_with_no_regions() {
        let mut chroms = ChromTable::new();
        let mut block = cases_block(&mut chroms);
        block.chrom = None;
        block.pos = None;
        block.id = None;
        block.alleles = None;
        block.qual = None;
        assert_eq!(block.fields(), Needs::GTS);
        let reader = GivenBlocks::of(vec![block], chroms);
        let bytes = write_vars(reader, Vec::new(), Some(3)).expect("the file was written");
        let file = file_read(bytes);

        assert_eq!(file.columns, vec![("gts".to_owned(), gts_type(6), false)]);
        assert_eq!(num_rows_of(&file.rows), [3, 1]);
        assert_eq!(
            file.batches,
            vec![
                BatchInfo {
                    num_vars: 3,
                    regions: Vec::new(),
                },
                BatchInfo {
                    num_vars: 1,
                    regions: Vec::new(),
                },
            ]
        );
    }

    /// A source with no variants is written and is not an error, as a VCF
    /// with no variants is read and is not one.
    #[test]
    fn a_source_with_no_variants_gives_a_file_with_both_keys_and_no_batch() {
        let reader = GivenBlocks::of(Vec::new(), ChromTable::new());
        let bytes = write_vars(reader, Vec::new(), Some(3)).expect("the file was written");
        let file = file_read(bytes);

        assert_eq!(file.columns, vec![("gts".to_owned(), gts_type(6), false)]);
        assert!(file.rows.is_empty());
        assert_eq!(file.batches, Vec::new());
        assert_eq!(file.metadata.individuals, cases_individuals());
        assert_eq!(file.metadata.ploidy, 2);
        assert_eq!(file.metadata.num_vars_per_block, 3);
    }

    /// The ploidy of the file is the one of its blocks, and the `gts`
    /// column holds the individuals times the ploidy alleles for each
    /// variant: a writer that wrote the individuals alone would give the
    /// same file for a diploid and for a tetraploid source.
    #[test]
    fn a_tetraploid_source_gives_a_file_of_twelve_alleles_for_each_variant() {
        let mut chroms = ChromTable::new();
        let number = chroms.intern("chr1");
        let block = Block {
            num_vars: 2,
            num_individuals: 3,
            ploidy: 4,
            gts: vec![
                0, 0, 1, 1, 0, 1, 1, 1, MISSING, MISSING, MISSING, MISSING, 0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 1,
            ],
            chrom: Some(vec![number, number]),
            pos: Some(vec![100, 200]),
            id: None,
            alleles: None,
            qual: None,
        };
        let reader = GivenBlocks {
            ploidy: 4,
            ..GivenBlocks::of(vec![block], chroms)
        };
        let bytes = write_vars(reader, Vec::new(), Some(2)).expect("the file was written");
        let file = file_read(bytes);

        assert_eq!(
            file.columns,
            vec![
                ("chrom".to_owned(), DataType::Utf8, false),
                ("pos".to_owned(), DataType::UInt64, false),
                // 3 individuals of the ploidy 4.
                ("gts".to_owned(), gts_type(12), false),
            ]
        );
        assert_eq!(file.metadata.ploidy, 4);
        assert_eq!(num_rows_of(&file.rows), [2]);
        assert_eq!(
            gts_of(&file.rows[0]).get(8..12),
            Some([MISSING; 4].as_slice())
        );
        assert_eq!(file.batches[0].num_vars, 2);
    }

    /// The memory the writer uses is one block: the vector of genotypes of
    /// the block it was given becomes the buffer of the `gts` column of the
    /// batch, at the address the vector had.
    #[test]
    fn the_genotypes_of_a_block_become_the_buffer_of_the_gts_column_with_no_copy() {
        let mut chroms = ChromTable::new();
        let block = cases_block(&mut chroms);
        let address = block.gts.as_ptr().addr();

        let column = gts_column(block.gts, 6).expect("the column of the genotypes");

        let alleles = column
            .as_any()
            .downcast_ref::<FixedSizeListArray>()
            .expect("the column is a fixed size list");
        assert_eq!(alleles.len(), 4);
        let values = alleles.values().as_primitive::<Int8Type>();
        assert_eq!(values.values().inner().as_ptr().addr(), address);
    }

    /// The buffers of every batch are compressed with lz4, which is the
    /// compression popnei writes and the one every build of it reads, and
    /// a caller that asks for no size of block gets the one popnei chooses
    /// for the individuals of the source.
    ///
    /// Arrow compresses each buffer of a batch on its own and writes the
    /// four bytes that mark the start of an lz4 frame before each, so those
    /// four bytes in the file are what says that the batch is compressed
    /// and with which of the two compressions of the format.
    #[test]
    fn the_buffers_of_the_file_are_compressed_with_lz4_and_popnei_chooses_the_size_of_the_blocks() {
        let mut chroms = ChromTable::new();
        let rows: Vec<&Row> = vec![&CASES[0]; 1000];
        let block = block_of(&rows, &mut chroms);
        let reader = GivenBlocks::of(vec![block], chroms);
        let bytes = write_vars(reader, Vec::new(), None).expect("the file was written");

        let frames = bytes
            .windows(LZ4_FRAME_MARK.len())
            .filter(|window| *window == LZ4_FRAME_MARK)
            .count();
        // The batch has fourteen buffers and each is compressed on its
        // own, the six that hold the values of the six columns among them.
        assert!(frames >= 6, "the file holds {frames} lz4 frames");
        // The genotypes of the thousand variants, which are the same
        // variant, are 1000 x 3 x 2 = 6000 bytes with no compression, and
        // the file holds no run of zeros of that length.
        let longest_run_of_zeros = bytes
            .split(|byte| *byte != 0)
            .map(<[u8]>::len)
            .max()
            .unwrap_or(0);
        assert!(
            longest_run_of_zeros < 1000,
            "the file holds a run of {longest_run_of_zeros} zeros"
        );

        let file = file_read(bytes);
        assert_eq!(num_rows_of(&file.rows), [1000]);
        // No size was asked for, so the blocks hold the number of variants
        // popnei chooses for 3 individuals, the largest it chooses.
        assert_eq!(file.metadata.num_vars_per_block, 10_000);
    }

    /// A vars file of the blocks given, written with `write_block` and the
    /// table of the names of their chromosomes, with what the writer said
    /// of the block at which it stopped.
    fn written_block_by_block(
        blocks: Vec<Block>,
        chroms: &ChromTable,
        individuals: &[String],
        ploidy: usize,
    ) -> (Vec<u8>, Option<Error>) {
        let mut writer =
            VarsWriter::new(Vec::new(), individuals, ploidy, 3).expect("the writer was built");
        let mut refused = None;
        for block in blocks {
            if let Err(error) = writer.write_block(block, chroms) {
                refused = Some(error);
                break;
            }
        }
        let bytes = writer.finish().expect("the file was finished");
        (bytes, refused)
    }

    /// The writer is built for the individuals and the ploidy that the
    /// `popnei` key of the file names, and every batch holds them: a block
    /// of others is refused and nothing of it is written.
    #[test]
    fn a_block_of_other_individuals_or_another_ploidy_is_refused_and_nothing_of_it_is_written() {
        let mut chroms = ChromTable::new();
        let mut block = cases_block(&mut chroms);
        block.num_individuals = 4;
        block.ploidy = 2;
        // 4 variants x 4 individuals x 2 alleles, so the block is of its
        // own size and what is wrong with it is the individuals alone.
        block.gts = vec![0; 32];
        let (bytes, refused) =
            written_block_by_block(vec![block], &chroms, &cases_individuals(), 2);
        let Some(Error::VarsBlockDoesNotFit {
            num_individuals,
            ploidy,
            found_num_individuals,
            found_ploidy,
        }) = refused
        else {
            panic!("the block of four individuals was taken: {refused:?}");
        };
        assert_eq!(
            (num_individuals, ploidy, found_num_individuals, found_ploidy),
            (3, 2, 4, 2)
        );
        assert!(file_read(bytes).rows.is_empty());

        let mut chroms = ChromTable::new();
        let mut block = cases_block(&mut chroms);
        block.ploidy = 1;
        // 4 variants x 3 individuals x 1 allele.
        block.gts = vec![0; 12];
        let (bytes, refused) =
            written_block_by_block(vec![block], &chroms, &cases_individuals(), 2);
        let Some(Error::VarsBlockDoesNotFit { found_ploidy, .. }) = refused else {
            panic!("the block of the ploidy 1 was taken: {refused:?}");
        };
        assert_eq!(found_ploidy, 1);
        assert!(file_read(bytes).rows.is_empty());
    }

    /// The genotypes of a batch are one flat buffer, so a block whose
    /// arrays are not of its size would be written with every genotype
    /// after the fault at the place of another: the writer calls
    /// `Block::check` and writes nothing of such a block.
    #[test]
    fn a_block_whose_arrays_are_not_of_its_size_is_refused_and_nothing_of_it_is_written() {
        let mut chroms = ChromTable::new();
        let mut block = cases_block(&mut chroms);
        block.gts.pop();
        let (bytes, refused) =
            written_block_by_block(vec![block], &chroms, &cases_individuals(), 2);

        let Some(Error::BlockArrayOfAnotherSize {
            array,
            found,
            expected,
        }) = refused
        else {
            panic!("the block short of an allele was taken: {refused:?}");
        };
        // 4 variants x 3 individuals x 2 alleles.
        assert_eq!((array, found, expected), ("gts", 23, 24));
        assert!(file_read(bytes).rows.is_empty());
    }

    /// The file holds the name of the chromosome of every variant as text,
    /// so a block whose chromosome number the table given with it has no
    /// name for cannot be written.
    #[test]
    fn a_chromosome_number_with_no_name_is_refused_and_nothing_of_the_block_is_written() {
        let mut chroms = ChromTable::new();
        let mut block = cases_block(&mut chroms);
        // A table of one name, `chr1`, and a block that holds a second
        // number: the reader that gave it has a defect.
        block.chrom = Some(vec![0, 0, 7, 0]);
        let (bytes, refused) =
            written_block_by_block(vec![block], &chroms, &cases_individuals(), 2);

        let Some(Error::VarsChromNameMissing { number }) = refused else {
            panic!("the block of a number with no name was taken: {refused:?}");
        };
        assert_eq!(number, 7);
        assert!(file_read(bytes).rows.is_empty());
    }

    /// Every batch of an arrow file has the columns of one schema, which
    /// the first block written fixes: a later block of other columns is
    /// refused, the message names the field the two differ in, and the file
    /// holds the batches written before it.
    #[test]
    fn a_block_of_other_columns_than_the_first_is_refused_and_names_the_field() {
        let mut chroms = ChromTable::new();
        let first = block_of(&[&CASES[0], &CASES[1]], &mut chroms);
        let mut second = block_of(&[&CASES[2], &CASES[3]], &mut chroms);
        second.qual = None;
        let (bytes, refused) =
            written_block_by_block(vec![first, second], &chroms, &cases_individuals(), 2);

        let Some(error) = refused else {
            panic!("the block without the qualities was taken");
        };
        let Error::VarsBlockColumns { first, found } = &error else {
            panic!("the error is {error}");
        };
        assert_eq!(*first, Needs::ALL);
        assert_eq!(*found, Needs::ALL.difference(Needs::QUAL));
        let message = error.to_string();
        assert!(message.contains("differ in `qual`"), "{message}");

        // The block that was taken is in the file, and the one that was
        // refused is in neither the batches nor the footer.
        let file = file_read(bytes);
        assert_eq!(num_rows_of(&file.rows), [2]);
        assert_eq!(positions_of(&file.rows[0]), [Some(100), Some(200)]);
        assert_eq!(file.batches.len(), 1);
        assert_eq!(file.batches[0].num_vars, 2);
    }

    /// The four values of "What it holds" of `docs/specs/io_vars.md`, of a
    /// file of the three individuals of the `cases.vcf` table of
    /// `docs/specs/io_vcf.md`.
    fn metadata_of_cases() -> VarsMetadata {
        VarsMetadata {
            format_version: FORMAT_VERSION.to_owned(),
            individuals: vec!["ind1".to_owned(), "ind2".to_owned(), "ind3".to_owned()],
            ploidy: 2,
            num_vars_per_block: 3,
        }
    }

    /// The text of the key is what another program that opens the file
    /// reads, so the test holds it as a literal and not only what it parses
    /// back to.
    #[test]
    fn the_popnei_key_is_written_with_its_four_values_and_parsed_back() {
        let metadata = metadata_of_cases();
        let written = metadata_as_json(&metadata);
        assert_eq!(
            written,
            r#"{"format_version":"1.0","individuals":["ind1","ind2","ind3"],"ploidy":2,"num_vars_per_block":3}"#
        );
        let parsed = metadata_from_json(&written).expect("the key written is parsed back");
        assert_eq!(parsed, metadata);
    }

    /// The ploidy and the variants of a batch are two numbers of one key,
    /// and a file of another ploidy than 2 tells whether each is read from
    /// its own name. A name with a quotation mark in it is json and not
    /// text, and a reader of the file would stop at it.
    #[test]
    fn the_popnei_key_of_another_ploidy_and_of_a_name_that_needs_an_escape_is_parsed_back() {
        let metadata = VarsMetadata {
            format_version: "1.7".to_owned(),
            individuals: vec![r#"the "first" one"#.to_owned(), "ind\\2".to_owned()],
            ploidy: 4,
            num_vars_per_block: 100,
        };
        let written = metadata_as_json(&metadata);
        assert_eq!(
            written,
            r#"{"format_version":"1.7","individuals":["the \"first\" one","ind\\2"],"ploidy":4,"num_vars_per_block":100}"#
        );
        let parsed = metadata_from_json(&written).expect("the key written is parsed back");
        assert_eq!(parsed, metadata);
    }

    /// The entry of "What it holds", the batch of 100 variants of two
    /// chromosomes, which is the third batch of the file of `many.vcf` of
    /// "How it is verified" of the writer.
    #[test]
    fn the_batches_key_is_written_as_the_example_of_the_spec_and_parsed_back() {
        let batches = vec![BatchInfo {
            num_vars: 100,
            regions: vec![
                Region {
                    chrom: "chr1".to_owned(),
                    min_pos: 8400,
                    max_pos: 10213,
                },
                Region {
                    chrom: "chr2".to_owned(),
                    min_pos: 10250,
                    max_pos: 12063,
                },
            ],
        }];
        let written = batches_as_json(&batches);
        assert_eq!(
            written,
            r#"[{"num_vars":100,"regions":[{"chrom":"chr1","min_pos":8400,"max_pos":10213},{"chrom":"chr2","min_pos":10250,"max_pos":12063}]}]"#
        );
        let parsed = batches_from_json(&written).expect("the key written is parsed back");
        assert_eq!(parsed, batches);
    }

    /// A file whose source gave the genotypes alone has no `chrom` and no
    /// `pos` column, so there is nothing to say about where the variants of
    /// its batches are. A batch of no variants is one that no writer of
    /// popnei makes and another arrow program can.
    #[test]
    fn a_batch_with_no_regions_is_written_with_the_number_of_its_variants_alone() {
        let batches = vec![
            BatchInfo {
                num_vars: 3,
                regions: Vec::new(),
            },
            BatchInfo {
                num_vars: 0,
                regions: Vec::new(),
            },
        ];
        let written = batches_as_json(&batches);
        assert_eq!(written, r#"[{"num_vars":3},{"num_vars":0}]"#);
        let parsed = batches_from_json(&written).expect("the key written is parsed back");
        assert_eq!(parsed, batches);

        let parsed = batches_from_json("[]").expect("a file with no batch has an empty key");
        assert_eq!(parsed, Vec::new());
    }

    /// What a user is told when the file is another arrow file, or a vars
    /// file whose keys another program wrote: the message names the key
    /// that is not what a vars file holds.
    #[test]
    fn a_key_whose_value_is_not_json_is_not_a_vars_file() {
        let problem = problem_of(metadata_from_json("{not json"));
        assert!(problem.contains("`popnei`"), "{problem}");
        assert!(problem.contains("schema"), "{problem}");
        assert!(problem.contains("not json"), "{problem}");

        let problem = problem_of(batches_from_json("[{\"num_vars\": 3}"));
        assert!(problem.contains("`popnei_batches`"), "{problem}");
        assert!(problem.contains("footer"), "{problem}");
        assert!(problem.contains("not json"), "{problem}");

        let problem = problem_of(metadata_from_json("3"));
        assert!(problem.contains("not a json object"), "{problem}");
        let problem = problem_of(batches_from_json(r#"{"num_vars": 3}"#));
        assert!(problem.contains("not a json array"), "{problem}");
    }

    /// Each of the four values of the `popnei` key is read from its own
    /// name, and the message names the one the file does not have.
    #[test]
    fn a_popnei_key_that_lacks_one_of_its_four_values_names_the_one_that_is_missing() {
        let whole = metadata_as_json(&metadata_of_cases());
        for key in [
            "format_version",
            "individuals",
            "ploidy",
            "num_vars_per_block",
        ] {
            let without = whole.replace(&format!("\"{key}\""), "\"another_key\"");
            let problem = problem_of(metadata_from_json(&without));
            assert!(problem.contains(&format!("has no `{key}`")), "{problem}");
        }
        assert!(metadata_from_json(&whole).is_ok());
    }

    /// A value of the `popnei` key that is not what that key holds is the
    /// same for a reader as one that is not there: the file was not written
    /// by popnei, and the message says which value it is about.
    #[test]
    fn a_popnei_key_whose_values_are_not_what_it_holds_is_not_a_vars_file() {
        let problem = problem_of(metadata_from_json(
            r#"{"format_version":1.0,"individuals":["ind1"],"ploidy":2,"num_vars_per_block":3}"#,
        ));
        assert!(problem.contains("`format_version`"), "{problem}");
        assert!(problem.contains("not a text"), "{problem}");

        let problem = problem_of(metadata_from_json(
            r#"{"format_version":"1.0","individuals":"ind1","ploidy":2,"num_vars_per_block":3}"#,
        ));
        assert!(problem.contains("`individuals`"), "{problem}");
        assert!(problem.contains("not a json array"), "{problem}");

        let problem = problem_of(metadata_from_json(
            r#"{"format_version":"1.0","individuals":[7],"ploidy":2,"num_vars_per_block":3}"#,
        ));
        assert!(
            problem.contains("not the name of an individual"),
            "{problem}"
        );

        // A ploidy of 0 says that a genotype holds no allele, and a batch
        // of 0 variants is a size no reader takes: neither is a number the
        // key can hold, and a negative one is not a count at all.
        for value in ["0", "-2", "\"two\""] {
            let problem = problem_of(metadata_from_json(&format!(
                r#"{{"format_version":"1.0","individuals":["ind1"],"ploidy":{value},"num_vars_per_block":3}}"#
            )));
            assert!(problem.contains("`ploidy`"), "{problem}");
            assert!(problem.contains("1 or more"), "{problem}");
        }
        let problem = problem_of(metadata_from_json(
            r#"{"format_version":"1.0","individuals":["ind1"],"ploidy":2,"num_vars_per_block":0}"#,
        ));
        assert!(problem.contains("`num_vars_per_block`"), "{problem}");
        assert!(problem.contains("1 or more"), "{problem}");
    }

    /// The footer is what says how many variants the file holds and where
    /// they are, and a reader that took an entry it could not read for an
    /// empty one would give a file of fewer variants than it has.
    #[test]
    fn an_entry_of_the_footer_that_is_not_one_is_not_a_vars_file() {
        let problem = problem_of(batches_from_json("[3]"));
        assert!(problem.contains("not a json object"), "{problem}");

        let problem = problem_of(batches_from_json(r#"[{"regions":[]}]"#));
        assert!(problem.contains("has no `num_vars`"), "{problem}");

        let problem = problem_of(batches_from_json(r#"[{"num_vars":-1}]"#));
        assert!(problem.contains("`num_vars`"), "{problem}");
        assert!(problem.contains("0 or more"), "{problem}");

        let problem = problem_of(batches_from_json(r#"[{"num_vars":3,"regions":{}}]"#));
        assert!(problem.contains("`regions`"), "{problem}");
        assert!(problem.contains("not a json array"), "{problem}");

        let problem = problem_of(batches_from_json(
            r#"[{"num_vars":3,"regions":[{"min_pos":1,"max_pos":9}]}]"#,
        ));
        assert!(problem.contains("has no `chrom`"), "{problem}");

        let problem = problem_of(batches_from_json(
            r#"[{"num_vars":3,"regions":[{"chrom":"chr1","min_pos":1}]}]"#,
        ));
        assert!(problem.contains("has no `max_pos`"), "{problem}");

        let problem = problem_of(batches_from_json(
            r#"[{"num_vars":3,"regions":[{"chrom":"chr1","min_pos":"start","max_pos":9}]}]"#,
        ));
        assert!(problem.contains("`min_pos`"), "{problem}");
        assert!(
            problem.contains("not the position of a variant"),
            "{problem}"
        );
    }

    /// A user who gets one of these has a file that another program wrote,
    /// or one that was damaged, and what they do about it is in the
    /// message: which column, which variant, which batch.
    #[test]
    fn the_message_of_a_file_popnei_cannot_read_names_what_is_wrong_with_it() {
        let message = Error::VarsFormatVersion {
            found: "2.0".to_owned(),
            read: FORMAT_VERSION_READ,
        }
        .to_string();
        assert!(message.contains("2.0"), "{message}");
        assert!(message.contains("starts with 1"), "{message}");

        let message = Error::VarsColumnType {
            column: "pos",
            found: "Int32".to_owned(),
            expected: "UInt64".to_owned(),
        }
        .to_string();
        assert!(message.contains("`pos`"), "{message}");
        assert!(message.contains("Int32"), "{message}");
        assert!(message.contains("UInt64"), "{message}");

        let message = Error::VarsGtsWidth {
            found: 7,
            expected: 6,
            num_individuals: 3,
            ploidy: 2,
        }
        .to_string();
        assert!(message.contains("7 alleles"), "{message}");
        assert!(message.contains("3 individuals"), "{message}");
        assert!(message.contains("ploidy 2"), "{message}");
        assert!(message.contains("hold 6"), "{message}");

        let message = Error::VarsNullValue {
            column: "pos",
            var: 42,
        }
        .to_string();
        assert!(message.contains("`pos`"), "{message}");
        assert!(message.contains("variant 42"), "{message}");

        let message = Error::VarsBatchNumVars {
            batch: 2,
            found: 97,
            expected: 100,
        }
        .to_string();
        assert!(message.contains("batch 2"), "{message}");
        assert!(message.contains("97 variants"), "{message}");
        assert!(message.contains("says 100"), "{message}");

        // A user whose file is compressed with zstd has to write it again,
        // and the message says with what.
        let message = Error::VarsZstd.to_string();
        assert!(message.contains("zstd"), "{message}");
        assert!(message.contains("lz4"), "{message}");

        let message = Error::VarsBatchesDoNotMatch {
            found: 1,
            expected: 2,
        }
        .to_string();
        assert!(message.contains("1 entries"), "{message}");
        assert!(message.contains("2 batches"), "{message}");

        let message = Error::VarsIndividualTwice {
            name: "ind2".to_owned(),
        }
        .to_string();
        assert!(message.contains("`ind2`"), "{message}");

        let message = Error::VarsFileCutShort {
            problem: "the footer of the file is 8 bytes and 4 were read".to_owned(),
        }
        .to_string();
        assert!(message.contains("cut short"), "{message}");
        assert!(message.contains("4 were read"), "{message}");

        let message = Error::VarsBatchNotRead {
            batch: 3,
            problem: "lz4 decompression failed".to_owned(),
        }
        .to_string();
        assert!(message.contains("batch 3"), "{message}");
        assert!(message.contains("lz4 decompression failed"), "{message}");
    }

    /// Nobody gets one of these by what they write, so the message is for
    /// whoever reports the defect: it says what the writer was given and
    /// what it was built for.
    #[test]
    fn the_message_of_a_defect_of_the_writer_names_what_does_not_fit() {
        let message = Error::VarsBlockDoesNotFit {
            num_individuals: 3,
            ploidy: 2,
            found_num_individuals: 4,
            found_ploidy: 2,
        }
        .to_string();
        assert!(message.contains("3 individuals"), "{message}");
        assert!(message.contains("4 individuals"), "{message}");

        // The field the two blocks differ in is what says which reader gave
        // the second one and what it stopped giving.
        let message = Error::VarsBlockColumns {
            first: Needs::ALL,
            found: Needs::ALL.difference(Needs::QUAL),
        }
        .to_string();
        assert!(message.contains("differ in `qual`"), "{message}");
        assert!(
            message.contains("`gts`, `chrom and pos`, `id`, `alleles`, `qual`"),
            "{message}"
        );

        let message = Error::VarsChromNameMissing { number: 5 }.to_string();
        assert!(message.contains("number 5"), "{message}");
    }

    /// What a call that refused the source says, for a test that looks at
    /// the message. Anything else is a failure of the test itself.
    fn problem_of<T: std::fmt::Debug>(result: crate::error::Result<T>) -> String {
        match result {
            Err(Error::NotAVarsFile { problem }) => problem,
            other => panic!("the source was not refused as a vars file: {other:?}"),
        }
    }
}
