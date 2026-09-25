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
//! blocks and not a reader. [`VarsReader`] opens such a file, says what it
//! holds and gives each of its batches as a block.
//!
//! `docs/specs/io_vars.md` has the format, the writer and the reader.

#[cfg(not(target_family = "wasm"))]
use std::collections::VecDeque;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, ErrorKind, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Arc;

use arrow_array::builder::{ListBuilder, StringBuilder};
use arrow_array::{
    Array, ArrayRef, FixedSizeListArray, Float32Array, Int8Array, ListArray, RecordBatch,
    StringArray, UInt64Array,
};
use arrow_buffer::{Buffer, NullBuffer, ScalarBuffer};
use arrow_ipc::convert::try_fb_to_schema;
use arrow_ipc::reader::FileDecoder;
use arrow_ipc::writer::{FileWriter, IpcWriteOptions};
use arrow_ipc::{
    Block as ArrowBlock, CompressionType, Footer, MessageHeader, MetadataVersion,
    RecordBatch as BatchMessage, root_as_footer, root_as_message,
};
use arrow_schema::{ArrowError, DataType, Field, Fields, Schema, SchemaRef};
use serde_json::{Map, Value};

use crate::block::{AllelesColumn, Block, BlockReader, BlockSize, Reblock, size_of_the_blocks};
use crate::error::{Error, Result};
use crate::filters::FilteringStats;
use crate::variant::{ChromTable, MISSING_ALLELE, Needs};

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
/// not `FORMAT_VERSION_READ`, and reads a file with any minor version, a
/// later one than its own too, ignoring the keys and the columns it does
/// not know: that is what lets a later version of the format add a column
/// without making the files or the readers that are there useless.
pub const FORMAT_VERSION: &str = "1.0";

/// The major version of the format that popnei reads, the part of
/// [`FORMAT_VERSION`] before the dot.
pub(crate) const FORMAT_VERSION_READ: &str = "1";

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

/// How many bytes of text one column of a batch holds, which is
/// `i32::MAX`: arrow keeps where each text of a column ends in a 32 bit
/// number. A dataset with more text than that in one block is written in
/// more batches, which is what a smaller `num_vars_per_block` gives.
const MAX_COLUMN_BYTES: u64 = 2_147_483_647;

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
    /// When there is no individual or the ploidy is 0, which is a file of
    /// genotypes of no allele, and when `num_vars_per_block` is 0: the
    /// `popnei` key of a vars file says one individual at least, a ploidy
    /// of 1 at least and batches of 1 variant at least, and a reader of
    /// popnei refuses a file whose key says less.
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
        // What the `popnei` key of the file would say, which a reader of
        // popnei refuses when it says no individual, a genotype of no
        // allele or a batch of no variant.
        if num_individuals == 0 || ploidy == 0 {
            return Err(Error::VarsFileOfNoGenotypes {
                num_individuals,
                ploidy,
            });
        }
        if num_vars_per_block == 0 {
            return Err(Error::BlockOfNoVariants);
        }
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
            let (names, where_they_are) = chrom_column(&chrom, &pos, chroms, MAX_COLUMN_BYTES)?;
            regions = where_they_are;
            arrays.push(Arc::new(names));
            arrays.push(Arc::new(UInt64Array::new(ScalarBuffer::from(pos), None)));
        }
        if let Some(id) = id {
            arrays.push(Arc::new(id_column(&id, MAX_COLUMN_BYTES)?));
        }
        if let Some(alleles) = alleles {
            arrays.push(Arc::new(alleles_column(&alleles, MAX_COLUMN_BYTES)?));
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
/// block, and the sink back with how many variants were written.
///
/// That count is what a Python or a TypeScript user reads as the variants
/// of the pass, `num_vars` of the `PassStats` of `docs/specs/variant.md`.
/// The binding crate cannot count them itself, as it does in an
/// `iter_blocks`, because the loop over the blocks is here; the counts of
/// the filters of the pass it reads from the chain of readers it keeps,
/// which it lends here as `&mut reader`.
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
) -> Result<(W, u64)> {
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
    let mut num_vars: u64 = 0;
    while let Some(block) = blocks.next_block()? {
        let of_the_block = u64::try_from(block.num_vars).unwrap_or(u64::MAX);
        writer.write_block(block, blocks.chroms())?;
        // A variant is a row of the file that is being written, so a pass
        // of the 18446744073709551615 variants this count holds is more
        // rows than any file system takes: the sum cannot reach its end.
        // The conversion above cannot fail either: a `usize` is 64 bits
        // natively and 32 in wasm, and both fit in a `u64`.
        num_vars = num_vars.saturating_add(of_the_block);
    }
    Ok((writer.finish()?, num_vars))
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
/// variant. The values of the list carry the name pyarrow gives the values
/// of any list and say that they hold no null, which is what is true of
/// them: no allele of a variant is a null. pyarrow writes that field as one
/// that can hold nulls, and the reader of popnei takes either.
fn alleles_type() -> DataType {
    DataType::List(Arc::new(alleles_field()))
}

/// The field arrow gives the values inside the `alleles` column: the name
/// pyarrow gives the values of any list, and no null, which is what the
/// alleles of a block are. A field that can hold nulls costs a bit for each
/// value in the file.
fn alleles_field() -> Field {
    Field::new(ITEM_FIELD, DataType::Utf8, false)
}

/// The arrow type of the `gts` column, the `alleles_per_var` alleles of
/// each variant. A fixed size list keeps no offsets, so the column is one
/// flat buffer of variants x individuals x ploidy signed bytes.
fn gts_type(alleles_per_var: i32) -> DataType {
    DataType::FixedSizeList(Arc::new(gts_field()), alleles_per_var)
}

/// The field arrow gives the alleles inside the `gts` column: the name
/// pyarrow gives the values of any list, and no null, since an allele that
/// was not called is the -1 of the genotypes and not a value that is not
/// there. A field that can hold nulls costs a bit for each allele in the
/// file, which arrow-rs writes as a mask of ones.
fn gts_field() -> Field {
    Field::new(ITEM_FIELD, DataType::Int8, false)
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
/// When the names hold more than `largest` bytes, and when `chroms` has no
/// name for a number of the block.
fn chrom_column(
    chrom: &[u32],
    pos: &[u64],
    chroms: &ChromTable,
    largest: u64,
) -> Result<(StringArray, Vec<Region>)> {
    // A number the table has no name for counts as no text here; the loop
    // below is what refuses the block for it, with the number in the error.
    text_fits(
        CHROM_COLUMN,
        chrom
            .iter()
            .map(|number| chroms.name(*number).unwrap_or("").len()),
        largest,
    )?;
    let mut names = StringBuilder::new();
    let mut regions: Vec<Region> = Vec::new();
    for (number, position) in chrom.iter().copied().zip(pos.iter().copied()) {
        let Some(name) = chroms.name(number) else {
            return Err(Error::VarsChromNameMissing { number });
        };
        names.append_value(name);
        // The chromosomes of a batch are a handful, so the entry of this
        // one is looked for among them and never by an index that could
        // be one past them.
        match regions.iter_mut().find(|region| region.chrom == name) {
            Some(region) => {
                region.min_pos = region.min_pos.min(position);
                region.max_pos = region.max_pos.max(position);
            }
            None => regions.push(Region {
                chrom: name.to_owned(),
                min_pos: position,
                max_pos: position,
            }),
        }
    }
    Ok((names.finish(), regions))
}

/// The `id` column of one batch. A block holds an empty id for a variant
/// that has none and the file holds a null, which is what any other program
/// that opens it takes for a value that is not there.
///
/// # Errors
///
/// When the ids hold more than `largest` bytes.
fn id_column(ids: &[String], largest: u64) -> Result<StringArray> {
    text_fits(ID_COLUMN, ids.iter().map(String::len), largest)?;
    let mut column = StringBuilder::new();
    for id in ids {
        match id.is_empty() {
            true => column.append_null(),
            false => column.append_value(id),
        }
    }
    Ok(column.finish())
}

/// The `alleles` column of one batch, the reference allele of each variant
/// first and then its alternative ones.
///
/// # Errors
///
/// When the alleles hold more than `largest` bytes, and when they are more
/// than `largest` alleles: arrow keeps where the alleles of each variant
/// end in a 32 bit number as it keeps where each text ends, and a block of
/// empty alleles holds more entries than bytes.
fn alleles_column(alleles: &AllelesColumn, largest: u64) -> Result<ListArray> {
    let texts = (0..alleles.num_vars()).flat_map(|var| {
        (0..alleles.num_alleles(var)).map(move |allele| alleles.allele(var, allele).len())
    });
    text_fits(ALLELES_COLUMN, texts, largest)?;
    let mut num_alleles: u64 = 0;
    for var in 0..alleles.num_vars() {
        num_alleles =
            num_alleles.saturating_add(u64::try_from(alleles.num_alleles(var)).unwrap_or(u64::MAX));
    }
    count_fits(ALLELES_COLUMN, "alleles", num_alleles, largest)?;
    // The builder writes the field of the values of the column, which is
    // the one the schema of the file gives.
    let mut column = ListBuilder::new(StringBuilder::new()).with_field(Arc::new(alleles_field()));
    for var in 0..alleles.num_vars() {
        for allele in 0..alleles.num_alleles(var) {
            column.values().append_value(alleles.allele(var, allele));
        }
        column.append(true);
    }
    Ok(column.finish())
}

/// That the texts of one column of a block fit in one column of a batch.
///
/// `lengths` gives the bytes of each text of the column. The bytes are
/// counted before a builder of arrow-rs is given any of them, because the
/// builder panics at the text that goes past the limit and a panic there
/// would leave a half written file behind.
///
/// # Errors
///
/// When the texts hold more than `largest` bytes: the error names the
/// column, both counts and the way out, a smaller `num_vars_per_block`.
fn text_fits(
    column: &'static str,
    lengths: impl Iterator<Item = usize>,
    largest: u64,
) -> Result<()> {
    let mut found: u64 = 0;
    for length in lengths {
        found = found.saturating_add(u64::try_from(length).unwrap_or(u64::MAX));
    }
    count_fits(column, "bytes of text", found, largest)
}

/// That a count of one column of a block fits in one column of a batch.
///
/// # Errors
///
/// When it is more than `largest`: the error names the column, what was
/// counted, both numbers and the way out, a smaller `num_vars_per_block`.
fn count_fits(column: &'static str, counted: &'static str, found: u64, largest: u64) -> Result<()> {
    if found > largest {
        return Err(Error::VarsTextTooLarge {
            column,
            counted,
            found,
            largest,
        });
    }
    Ok(())
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
        Arc::new(gts_field()),
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
/// The error the file system gave is kept with the number it carries, since
/// that number is what a binding crate builds the exception of its language
/// with; what arrow-rs says of anything else goes in as text.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "arrow-rs has twenty cases of its error and popnei tells one of them, the error of \
              the file system, from every other"
)]
fn not_written(problem: ArrowError) -> Error {
    match problem {
        ArrowError::IoError(_, failure) => Error::VarsFileNotWritten {
            problem: failure.to_string(),
            source: Some(failure),
        },
        other => Error::VarsFileNotWritten {
            problem: other.to_string(),
            source: None,
        },
    }
}

/// The error of a writer whose sink is gone, which is what is left after
/// the header of the file could not be written.
fn the_sink_is_gone() -> Error {
    Error::VarsFileNotWritten {
        problem: "the header of the file could not be written, and the writer has nothing left \
                  to write on"
            .to_owned(),
        source: None,
    }
}

/// The six bytes an arrow IPC file starts with, and the six it ends with.
const ARROW_MAGIC: [u8; 6] = *b"ARROW1";

/// How many bytes the trailer of an arrow IPC file holds: the length of the
/// footer as four bytes, and then [`ARROW_MAGIC`] again.
const TRAILER_BYTES: u64 = 10;

/// How many bytes the message of a batch of an arrow IPC file starts with:
/// the four of the mark of a continuation and the four that say how long the
/// message is. arrow-rs reads them before anything else of a batch, so a
/// batch of fewer bytes than these is refused before it is handed over.
const MESSAGE_START_BYTES: usize = 8;

/// The four bytes that a message of an arrow IPC file may start with, which
/// say that a length follows. A message that has them starts 8 bytes in, and
/// one that has not, 4.
const CONTINUATION_MARK: [u8; 4] = [0xff; 4];

/// How many bytes the length of an uncompressed buffer takes at the start of
/// a compressed one: arrow writes it before the compressed bytes, and -1
/// there says that what follows is not compressed.
const UNCOMPRESSED_LENGTH_BYTES: usize = 8;

/// What -1 in those bytes says: the bytes that follow are not compressed.
const NOT_COMPRESSED: i64 = -1;

/// The most an lz4 frame gives back for each byte it holds. The format
/// encodes a match of 255 bytes in one, so a buffer of `n` bytes holds at
/// most 255 times that many, and popnei refuses a buffer that says more:
/// arrow-rs asks the machine for the memory of what the buffer says before
/// it decompresses, and a damaged length there ends the process. `zstd`
/// gives back far more for a byte, and no build of popnei decompresses it.
const LZ4_BYTES_FOR_A_BYTE: u64 = 255;

/// What is allowed above that for the header and the footer of the frame
/// and for the smallest buffers, where the ratio alone is too tight.
const LZ4_FRAME_SLACK: u64 = 1024;

/// How many values a column of a list holds at most, which is `i32::MAX`:
/// arrow keeps where the values of each row of such a column end in a 32
/// bit number, and the last of those numbers is how many values there are.
const MAX_VALUES_OF_A_LIST: u64 = 2_147_483_647;

/// The number the system gives for a directory where a file was asked for,
/// `EISDIR`, which is 21 on macOS, on Linux and in emscripten, the systems
/// popnei runs on. Opening a directory succeeds on those systems and only
/// the first read of it fails, so [`VarsReader::from_path`] gives this
/// number itself.
const A_DIRECTORY_IS_THERE: i32 = 21;

/// One column of the table of "What it holds" of `docs/specs/io_vars.md`,
/// the six a vars file can have and popnei knows. A column of any other
/// name is ignored, which is what lets a later version of the format add
/// one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VarsColumn {
    /// The name of the chromosome of each variant.
    Chrom,
    /// The position of each variant.
    Pos,
    /// The id of each variant, null for a variant with none.
    Id,
    /// The alleles of each variant, the reference one first.
    Alleles,
    /// The quality of each variant, null for a variant with none.
    Qual,
    /// The genotypes, the individuals times the ploidy for each variant.
    Gts,
}

impl VarsColumn {
    /// The column of that name, or `None` for a name popnei does not know.
    fn of_the_name(name: &str) -> Option<VarsColumn> {
        match name {
            CHROM_COLUMN => Some(VarsColumn::Chrom),
            POS_COLUMN => Some(VarsColumn::Pos),
            ID_COLUMN => Some(VarsColumn::Id),
            ALLELES_COLUMN => Some(VarsColumn::Alleles),
            QUAL_COLUMN => Some(VarsColumn::Qual),
            GTS_COLUMN => Some(VarsColumn::Gts),
            _ => None,
        }
    }

    /// Its name in the file.
    fn name(self) -> &'static str {
        match self {
            VarsColumn::Chrom => CHROM_COLUMN,
            VarsColumn::Pos => POS_COLUMN,
            VarsColumn::Id => ID_COLUMN,
            VarsColumn::Alleles => ALLELES_COLUMN,
            VarsColumn::Qual => QUAL_COLUMN,
            VarsColumn::Gts => GTS_COLUMN,
        }
    }

    /// The arrow type it holds, as the table of the spec writes it, which
    /// is the type the message of a column of another one gives.
    fn arrow_type(self) -> &'static str {
        match self {
            VarsColumn::Chrom | VarsColumn::Id => "Utf8",
            VarsColumn::Pos => "UInt64",
            VarsColumn::Alleles => "List<Utf8>",
            VarsColumn::Qual => "Float32",
            VarsColumn::Gts => "FixedSizeList<Int8>",
        }
    }

    /// Whether the column of the file is of the type it holds.
    ///
    /// What is compared is the type of the values: for the two lists, the
    /// type inside the list, and for `gts` that the width the list gives is
    /// a count. Neither the name of the field arrow gives the values of a
    /// list nor whether that field takes nulls is compared: popnei writes
    /// the `item` that takes nulls that pyarrow writes, and arrow programs
    /// differ in both, so a file whose list holds texts under another name
    /// is read.
    fn holds(self, found: &DataType) -> bool {
        match self {
            VarsColumn::Chrom | VarsColumn::Id => *found == DataType::Utf8,
            VarsColumn::Pos => *found == DataType::UInt64,
            VarsColumn::Alleles => {
                matches!(found, DataType::List(item) if *item.data_type() == DataType::Utf8)
            }
            VarsColumn::Qual => *found == DataType::Float32,
            VarsColumn::Gts => matches!(
                found,
                DataType::FixedSizeList(item, width)
                    if *item.data_type() == DataType::Int8 && *width >= 0
            ),
        }
    }

    /// Where this column's place in the file is kept.
    fn place_in(self, columns: &mut VarsColumns) -> &mut Option<usize> {
        match self {
            VarsColumn::Chrom => &mut columns.chrom,
            VarsColumn::Pos => &mut columns.pos,
            VarsColumn::Id => &mut columns.id,
            VarsColumn::Alleles => &mut columns.alleles,
            VarsColumn::Qual => &mut columns.qual,
            VarsColumn::Gts => &mut columns.gts,
        }
    }
}

/// Where each column popnei knows is in the file, by its place in the
/// schema, and `None` for one the file does not have: a file holds only the
/// columns its source could fill. The places are what the projection that a
/// [`Needs`] becomes is built from, which is why they come from the schema
/// of the file and not from the table of the spec.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct VarsColumns {
    /// Where `chrom` is.
    chrom: Option<usize>,
    /// Where `pos` is.
    pos: Option<usize>,
    /// Where `id` is.
    id: Option<usize>,
    /// Where `alleles` is.
    alleles: Option<usize>,
    /// Where `qual` is.
    qual: Option<usize>,
    /// Where `gts` is.
    gts: Option<usize>,
}

impl VarsColumns {
    /// Where the columns popnei knows are in that schema, each checked
    /// against the type it holds.
    ///
    /// # Errors
    ///
    /// When a column popnei knows is of another arrow type, and when the
    /// `gts` column holds another number of alleles for each variant than
    /// the individuals and the ploidy of the `popnei` key give.
    fn of(schema: &Schema, metadata: &VarsMetadata) -> Result<VarsColumns> {
        let mut columns = VarsColumns {
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
            gts: None,
        };
        for (place, field) in schema.fields().iter().enumerate() {
            let Some(column) = VarsColumn::of_the_name(field.name()) else {
                continue;
            };
            let found = field.data_type();
            if !column.holds(found) {
                return Err(Error::VarsColumnType {
                    column: column.name(),
                    found: type_name(found),
                    expected: column.arrow_type().to_owned(),
                });
            }
            if let DataType::FixedSizeList(_, width) = found
                && column == VarsColumn::Gts
            {
                gts_width_fits(*width, metadata)?;
            }
            // A file with two columns of one name is read with the first of
            // them, as any program that asks a table for a column by its
            // name reads it; the second is checked against its type too.
            column.place_in(&mut columns).get_or_insert(place);
        }
        Ok(columns)
    }
}

/// That the `gts` column holds the alleles of the individuals that the
/// `popnei` key names: the width of that column is what turns its flat
/// buffer into variants, so a file whose two say different things is
/// refused and not read one allele beside another.
///
/// # Errors
///
/// When the width is not the individuals times the ploidy, with both
/// numbers.
fn gts_width_fits(width: i32, metadata: &VarsMetadata) -> Result<()> {
    let num_individuals = metadata.individuals.len();
    let ploidy = metadata.ploidy;
    // The individuals are a vector of this machine and the ploidy comes
    // from the file, so their product can carry past what a `usize` counts;
    // a width, which an arrow file keeps in 32 bits, is never the product
    // that saturated, so such a file is refused here.
    let expected = num_individuals.saturating_mul(ploidy);
    let found = usize::try_from(width).unwrap_or(usize::MAX);
    if found != expected {
        return Err(Error::VarsGtsWidth {
            found,
            expected,
            num_individuals,
            ploidy,
        });
    }
    Ok(())
}

/// The arrow type of a column as the table of "What it holds" of
/// `docs/specs/io_vars.md` writes it, so that the type a file has and the
/// type popnei reads are written the same way in one message.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "arrow has forty types and the two lists are what popnei writes in its own notation; \
              every other one is written as arrow writes it"
)]
fn type_name(found: &DataType) -> String {
    match found {
        DataType::List(item) => format!("List<{inside}>", inside = type_name(item.data_type())),
        DataType::FixedSizeList(item, width) => {
            format!(
                "FixedSizeList<{inside}>[{width}]",
                inside = type_name(item.data_type())
            )
        }
        other => other.to_string(),
    }
}

/// Where one batch of a vars file is, from the footer of the arrow file:
/// the message of the batch and then its buffers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BatchAt {
    /// The byte of the file where the message of the batch starts.
    offset: u64,
    /// How many bytes that message holds.
    metadata_len: u64,
    /// How many bytes the buffers after it hold.
    body_len: u64,
}

/// The reader of a vars file: what the file says about itself, from its
/// schema and its footer, and each of its batches as a block.
///
/// [`VarsReader::new`] reads the schema and the footer, so the names of the
/// individuals, the ploidy, the variants of the file and the regions of
/// every batch are known when it returns, before any batch is read, and a
/// source that is not a vars file fails there and not at the first block.
///
/// The source is anything that can be read and seeked in, natively a file
/// and in a tab the bytes the user picked: the footer of an arrow file is
/// at its end, so the reader seeks there when it is opened and then to each
/// batch in turn.
///
/// It is a [`BlockReader`]: the batches come out as blocks, as they are in
/// the file, with only the columns that
/// [`set_needs`](BlockReader::set_needs) asks for decompressed.
pub struct VarsReader<R: Read + Seek> {
    /// The bytes of the file.
    source: R,
    /// The columns of the file as arrow gives them, with the column popnei
    /// does not know among them, which is what a batch is decoded with.
    schema: SchemaRef,
    /// Which version of the messages of arrow the file was written with,
    /// which the decoder of a batch is built with.
    version: MetadataVersion,
    /// Where each batch of the file is, in the order of the batches.
    ///
    /// A batch is read with the `FileDecoder` of arrow-rs, which takes the
    /// bytes of one batch and the columns to decompress, and not with its
    /// `FileReader`, which takes the columns when it is built and holds the
    /// source: a consumer changes which fields it asks for between two
    /// blocks, and `docs/specs/io_vars.md` asks that the columns of a batch
    /// be chosen when that batch is read.
    blocks: Vec<BatchAt>,
    /// Where each column popnei knows is in the schema.
    columns: VarsColumns,
    /// What the `popnei` key of the schema says.
    metadata: VarsMetadata,
    /// What the `popnei_batches` key of the footer says of each batch.
    batches: Vec<BatchInfo>,
    /// The variants of the whole file, the sum of those of its batches.
    num_vars: usize,
    /// Which fields the consumer asks for, which the projection of the
    /// next batch that is read is built from.
    needs: Needs,
    /// Which batch of the file is the next one whose bytes are read.
    next: usize,
    /// The batches that were decoded and not yet given, in the order of the
    /// file, each with its place among the batches, and the first error of a
    /// decode as the last entry after them.
    ///
    /// The batches of a window are decoded at once, and a consumer takes one
    /// block per call, so the ones that are decoded and not asked for yet
    /// wait here. The error goes in as an entry of its own, at the end,
    /// because the blocks of the batches before it are given first, as they
    /// were when each batch was decoded on its own.
    #[cfg(not(target_family = "wasm"))]
    decoded: VecDeque<(usize, Result<RecordBatch>)>,
    /// The names of the chromosomes of the variants that were given, each
    /// with its number.
    chroms: ChromTable,
    /// How many variants the batches before each batch of the file hold,
    /// one entry for each batch, which is what the variant of the error of a
    /// null is counted from.
    ///
    /// It is the prefix sum of what the `popnei_batches` key of the footer
    /// says each batch holds, worked out when the file is opened, and not a
    /// count of the rows of the batches that were read: a batch is decoded
    /// before the ones in front of it have been. It is the same number,
    /// because a batch whose rows are not as many as its entry of the footer
    /// says is refused before any batch after it is given.
    vars_before: Vec<u64>,
    /// Whether the reader gave its last block or an error. After either,
    /// every call gives no block.
    finished: bool,
}

impl<R: Read + Seek> VarsReader<R> {
    /// The reader over the vars file in `source`, whose schema and footer
    /// it reads.
    ///
    /// # Errors
    ///
    /// When the source is not a vars file: it does not start as an arrow
    /// IPC file, its footer is not one of an arrow file, its schema has no
    /// `popnei` key or that key does not hold its four values, or its
    /// footer has no `popnei_batches` key. When the version of the format
    /// it gives does not start with the one popnei reads, when a column
    /// popnei knows is of another arrow type, when the `gts` column holds
    /// another number of alleles for each variant than the `popnei` key
    /// gives, when that key names one individual twice, and when the
    /// entries of the footer are not as many as the batches of the file.
    /// When the file starts as an arrow file and was cut short. And when
    /// the source cannot be read.
    pub fn new(mut source: R) -> Result<VarsReader<R>> {
        let file_len = source.seek(SeekFrom::End(0))?;
        starts_as_an_arrow_file(&mut source, file_len)?;
        let footer_bytes = footer_of(&mut source, file_len)?;
        let footer = root_as_footer(&footer_bytes).map_err(|problem| {
            not_a_vars_file(format!(
                "it starts as an arrow IPC file and its footer is not one of an arrow file: {problem}"
            ))
        })?;
        let schema = Arc::new(schema_of_the_footer(&footer)?);
        // The version is read before the columns: a file of a version
        // popnei does not know holds whatever that version says, and what a
        // user does about it is get a later popnei and not a file of other
        // columns.
        let metadata = metadata_of_the_schema(&schema)?;
        let columns = VarsColumns::of(&schema, &metadata)?;
        // Every vars file holds the genotypes, so a file without them is
        // not one: a reader of this version cannot tell a file whose column
        // was lost from one that a later version of the format wrote with
        // no such column.
        if columns.gts.is_none() {
            return Err(not_a_vars_file(format!(
                "it has no `{GTS_COLUMN}` column, which every vars file has"
            )));
        }
        let blocks = batches_of_the_footer(&footer, file_len)?;
        let batches = batch_info_of_the_footer(&footer, blocks.len())?;
        let num_vars = num_vars_of_the_file(&batches, &metadata)?;
        let vars_before = vars_before_each_batch(&batches);
        Ok(VarsReader {
            source,
            schema,
            version: footer.version(),
            blocks,
            columns,
            metadata,
            batches,
            num_vars,
            needs: Needs::ALL,
            next: 0,
            #[cfg(not(target_family = "wasm"))]
            decoded: VecDeque::new(),
            chroms: ChromTable::new(),
            vars_before,
            finished: false,
        })
    }

    /// What the `popnei` key of the schema of the file says: the version of
    /// the format, the names of the individuals, the ploidy and how many
    /// variants a batch holds.
    #[must_use]
    pub fn metadata(&self) -> &VarsMetadata {
        &self.metadata
    }

    /// What the `popnei_batches` key of the footer says, one entry for each
    /// batch of the file, in the order of the batches: how many variants it
    /// holds and where they are.
    #[must_use]
    pub fn batches(&self) -> &[BatchInfo] {
        &self.batches
    }

    /// The variants of the whole file, the sum of those of its batches,
    /// which the footer says without a batch being read.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }
}

/// Which columns of the file a consumer that asks for `needs` has
/// decompressed, and what each of them fills, in the order of the columns of
/// the file.
///
/// A field that is asked for and whose column the file lacks is not there,
/// and the block comes without that column, which is the rule of
/// `docs/specs/variant.md` for a source that has no such field. The
/// chromosome and the position are one field and go together, so a file with
/// one of the two columns and not the other gives neither.
fn projection_of(needs: Needs, columns: &VarsColumns) -> Vec<(VarsColumn, usize)> {
    let mut wanted: Vec<(VarsColumn, usize)> = Vec::new();
    if needs.contains(Needs::CHROM_POS)
        && let (Some(chrom), Some(pos)) = (columns.chrom, columns.pos)
    {
        wanted.push((VarsColumn::Chrom, chrom));
        wanted.push((VarsColumn::Pos, pos));
    }
    for (field, place, column) in [
        (Needs::ID, columns.id, VarsColumn::Id),
        (Needs::ALLELES, columns.alleles, VarsColumn::Alleles),
        (Needs::QUAL, columns.qual, VarsColumn::Qual),
        (Needs::GTS, columns.gts, VarsColumn::Gts),
    ] {
        if needs.contains(field)
            && let Some(place) = place
        {
            wanted.push((column, place));
        }
    }
    wanted.sort_by_key(|(_, place)| *place);
    wanted
}

/// Where a batch is, for the messages of the errors it gives.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BatchPlace {
    /// Which batch of the file it is, counted from 1.
    batch: u64,
    /// How many variants the batches before it hold, which the variant of
    /// the error of a null is counted from.
    vars_before: u64,
    /// How many variants its entry of the footer says it holds, which the
    /// error of a batch the machine has no memory for names.
    num_vars: usize,
}

/// The most batches of a vars file that are decoded at once, on the threads
/// of rayon: the window is this many or the threads of the pool the reader
/// runs in, whichever is fewer, so a pool of one thread decodes one batch at
/// a time and a pool of 18 decodes eight.
///
/// The decompression of the buffers is what bounds a pass over a vars file:
/// of the 5.5 ms a batch of 5000 variants of 1000 individuals takes, 4.99 ms
/// is lz4. The genotypes of such a batch are one lz4 frame of three
/// independent blocks of 4 MiB, so decompressing the blocks of one batch on
/// the threads takes that batch from 4.534 to 1.972 ms and no further.
/// Decoding whole batches at once is the only route that scales past it, and
/// this is how many are in flight.
///
/// What it costs is memory. The reader holds the compressed bytes of the
/// batches of a window while they are decoded, and then their decoded arrow
/// buffers until a consumer has asked for each of them: for 1000 individuals
/// and 5000 variants in a batch, 3.9 MB compressed and 11.4 MB decoded for
/// each batch of the window. On `bigcalled.vars`, 100000 variants of 1000
/// individuals in 20 batches, `/usr/bin/time -l` of the `vars_file`
/// benchmark reads 495.6 MB of maximum resident set size before this and
/// 623.0 MB after it at 18 threads, 538.4 MB at one.
///
/// The value is the knee of a sweep of 1, 2, 4, 8 and 18 at 18 threads on
/// the owner's Apple M5 Pro of 18 cores, over that file, the best of 10 runs
/// of the genotypes alone: 107.8 ms at 1, 62.3 at 2, 37.7 at 4, 27.8 at 8
/// and 24.2 at 18. A window of 18 is another 10% of the best time for 2.25
/// times the memory of a window of 8, and its median over four sets of runs,
/// 26.8 to 30.7 ms, is no better than the 28.7 to 29.9 ms of the window of
/// 8. A window smaller than the pool leaves the threads above its size with
/// nothing to do, which is what the 24.2 ms at 18 are: the report of the
/// performance review of the vars reader has the whole sweep.
///
/// Why the pool bounds the window too: a window of eight decodes about 80 MB
/// of arrow buffers before a consumer touches the first of them, and this
/// machine has 128 KB of first level data cache for each performance core
/// and one 16 MB second level shared by six of them, so with one thread
/// every block has left the caches before its genotypes are read, where one
/// batch at a time decoded each block and then had it consumed while it was
/// still warm. Over `bigcalled.vars`, the Kosman distance of every pair of
/// its 1000 individuals in a pool of one thread, the best of three runs of
/// `cargo bench --bench kosman_dists`, went from 0.880 s with one batch at a
/// time to 1.026 s with a window of eight, 17% more, where the pass of the
/// `vars_file` benchmark, whose consumer only sums the genotypes, lost 1.1%.
/// Bounding the window by the pool gives the one thread its 0.880 s back and
/// keeps what the 18 threads won. The two binaries were run one after the other, four times
/// over, because the load of this machine moves a median by a fifth; the
/// three sets the load left alone read, on one thread, 0.885, 0.880 and
/// 0.884 s with one batch at a time against 0.891, 0.883 and 0.882 s with the
/// window bounded by the pool, and on 18 threads 0.218, 0.217 and 0.218 s
/// against 0.134, 0.128 and 0.129 s. The fourth set was thrown away: the same
/// calculation over blocks already in memory, which no reader is in, read
/// 0.785 to 1.401 s within it.
#[cfg(not(target_family = "wasm"))]
const BATCHES_AT_ONCE: usize = 8;

/// How many variants the batches before each batch hold, one entry for each
/// batch, which is the prefix sum of what the footer says each one holds.
///
/// The count of a whole file was checked to fit in a `usize` when the file
/// was opened, so no sum here is above what a `u64` holds and none of them
/// saturates.
fn vars_before_each_batch(batches: &[BatchInfo]) -> Vec<u64> {
    let mut before = Vec::with_capacity(batches.len());
    let mut so_far: u64 = 0;
    for batch in batches {
        before.push(so_far);
        so_far = so_far.saturating_add(u64::try_from(batch.num_vars).unwrap_or(u64::MAX));
    }
    before
}

impl<R: Read + Seek> VarsReader<R> {
    /// Where the batch at `index` is, counted from 0: which batch of the
    /// file it is, how many variants the batches before it hold and how many
    /// its entry of the footer says it holds.
    ///
    /// # Errors
    ///
    /// When the file has no entry of the footer for that batch, which
    /// `new` refuses a file for.
    fn place_of(&self, index: usize) -> Result<BatchPlace> {
        // The two vectors are of one entry for each entry of the footer, so
        // either both are there or the footer has no entry for this batch,
        // which is what the error says.
        let (Some(info), Some(vars_before)) =
            (self.batches.get(index), self.vars_before.get(index))
        else {
            return Err(Error::VarsBatchesDoNotMatch {
                found: self.batches.len(),
                expected: self.blocks.len(),
            });
        };
        Ok(BatchPlace {
            batch: counted_from_one(index),
            vars_before: *vars_before,
            num_vars: info.num_vars,
        })
    }

    /// The bytes of the batch at `at`, read from the source.
    ///
    /// The source is one `Read + Seek`, so this is the part of reading a
    /// batch that no thread but the caller's does.
    ///
    /// # Errors
    ///
    /// When the batch holds more bytes than this machine counts, and when
    /// the source cannot be read.
    fn bytes_of_the_batch(&mut self, at: BatchAt, place: BatchPlace) -> Result<Vec<u8>> {
        // Both lengths were checked to lie inside the file when it was
        // opened, so their sum is one of its bytes.
        let len = at.metadata_len.saturating_add(at.body_len);
        // A batch of more bytes than this machine counts, which under wasm,
        // where a `usize` is 32 bits, is 4 GB: the file is more than this
        // build of popnei reads, and the way out is smaller batches.
        let Ok(len) = usize::try_from(len) else {
            return Err(block_too_large(place.num_vars, &self.metadata));
        };
        bytes_at(&mut self.source, at.offset, len)
    }

    /// The block of the batch that was decoded at `index`, with the checks
    /// that were always made on the thread that asks for it.
    ///
    /// # Errors
    ///
    /// When the batch holds another number of variants than its entry of the
    /// footer says, and when one of its columns holds a null where every
    /// variant has a value.
    fn block_of_a_decoded_batch(
        &mut self,
        index: usize,
        batch: &RecordBatch,
        wanted: &[(VarsColumn, usize)],
    ) -> Result<Option<Block>> {
        let place = self.place_of(index)?;
        let num_vars = batch.num_rows();
        if num_vars != place.num_vars {
            return Err(Error::VarsBatchNumVars {
                batch: place.batch,
                found: num_vars,
                expected: place.num_vars,
            });
        }
        // A batch of no variants is not given as a block, since
        // `docs/specs/block.md` says that a reader never gives one: the
        // caller takes the next batch, as a filter does with a block it
        // emptied.
        if num_vars == 0 {
            return Ok(None);
        }
        Ok(Some(block_of_the_batch(
            batch,
            wanted,
            &self.metadata,
            &mut self.chroms,
            place,
        )?))
    }
}

#[cfg(not(target_family = "wasm"))]
impl<R: Read + Seek> VarsReader<R> {
    /// The next batch of the file as a block, the batches of no variant
    /// passed over, and `None` when there are no more.
    ///
    /// The bytes of a window of batches, as many as [`BATCHES_AT_ONCE`] says
    /// at the most, are read from the source, in the order of the file, and
    /// those batches are decoded on the threads of the pool the caller is
    /// in; the blocks then come out of
    /// that window one call at a time, in the order of the file, and the
    /// checks of each batch and the building of each block are made on the
    /// thread that asks for it. The chromosome table hands its numbers out
    /// in the order the names are first seen, and nine other places of
    /// popnei index their results by a running count of the variants, so a
    /// block out of order would be a wrong result and not a slower one.
    ///
    /// # Errors
    ///
    /// When a batch cannot be read, when it holds another number of
    /// variants than its entry of the footer, and when one of its columns
    /// holds a null where every variant has a value. The error is that of
    /// the first batch of the window in the order of the file that gave one,
    /// and the blocks of the batches before it are given first.
    fn next_batch(&mut self) -> Result<Option<Block>> {
        // The columns to decompress are chosen when the batch is read, so a
        // consumer that asks for other fields is served from here on.
        let wanted = projection_of(self.needs, &self.columns);
        let places: Vec<usize> = wanted.iter().map(|(_, place)| *place).collect();
        loop {
            if self.decoded.is_empty() {
                self.decode_the_next_batches(&places);
            }
            let Some((index, decoded)) = self.decoded.pop_front() else {
                return Ok(None);
            };
            let batch = decoded?;
            if let Some(block) = self.block_of_a_decoded_batch(index, &batch, &wanted)? {
                return Ok(Some(block));
            }
        }
    }

    /// It reads the bytes of the next batches of the file, as many as the
    /// window holds, decodes them on the threads of rayon and puts what each
    /// gave into the queue, in the order of the file.
    ///
    /// The window is [`BATCHES_AT_ONCE`] or the threads of the pool this runs
    /// in, whichever is fewer, for the reason that constant gives.
    /// `rayon::current_num_threads` is the pool of the caller: 1 inside a
    /// pool of one thread, and the threads of the global pool outside any
    /// `install`, which is one for each core of the machine.
    ///
    /// The bytes are read one batch after another, because they come from one
    /// `Read + Seek`. The results are collected in the order of the file and
    /// walked in that order, and the first error stops the queue, so the
    /// error a consumer is given is the one of the first batch of the window
    /// that failed and not of whichever thread failed first: `try_for_each`
    /// or `find_any` would make the message depend on the threads.
    ///
    /// It gives nothing back: an error of a batch goes into the queue behind
    /// the batches before it, whose blocks a consumer is given first, as it
    /// was when the batches were read one at a time.
    fn decode_the_next_batches(&mut self, places: &[usize]) {
        use rayon::iter::{IntoParallelIterator, ParallelIterator};

        let mut fetched: Vec<(usize, BatchPlace, BatchAt, Vec<u8>)> = Vec::new();
        let mut failed: Option<(usize, Error)> = None;
        let window = BATCHES_AT_ONCE.min(rayon::current_num_threads());
        while fetched.len() < window {
            let Some(at) = self.blocks.get(self.next).copied() else {
                break;
            };
            let index = self.next;
            let fetch = self
                .place_of(index)
                .and_then(|place| Ok((place, self.bytes_of_the_batch(at, place)?)));
            let (place, bytes) = match fetch {
                Ok(fetched) => fetched,
                Err(problem) => {
                    failed = Some((index, problem));
                    break;
                }
            };
            // The batches of a file are as many as the machine counts, so
            // this never saturates.
            self.next = self.next.saturating_add(1);
            fetched.push((index, place, at, bytes));
        }
        // Neither of these is of the reader, so the decode touches nothing
        // that is: the schema is shared and the version is a number.
        let schema = &self.schema;
        let version = self.version;
        let decoded: Vec<(usize, Result<RecordBatch>)> = fetched
            .into_par_iter()
            .map(|(index, place, at, bytes)| {
                (index, batch_of(schema, version, places, at, bytes, place))
            })
            .collect();
        for (index, batch) in decoded {
            let failed = batch.is_err();
            self.decoded.push_back((index, batch));
            if failed {
                return;
            }
        }
        if let Some((index, problem)) = failed {
            self.decoded.push_back((index, Err(problem)));
        }
    }

    /// It throws the batches that were decoded and not given away, and the
    /// next batch to read is the first of them again.
    ///
    /// A batch of the queue was decoded with the projection of the fields
    /// that were asked for before, so a consumer that asks for others is
    /// served from the next block on, which is what `set_needs` says. The few
    /// batches of the window are decoded again, which costs one window and is
    /// right.
    fn throw_the_queue_away(&mut self) {
        if let Some((index, _)) = self.decoded.front() {
            self.next = *index;
        }
        self.decoded.clear();
    }
}

#[cfg(target_family = "wasm")]
impl<R: Read + Seek> VarsReader<R> {
    /// The next batch of the file as a block, the batches of no variant
    /// passed over, and `None` when there are no more.
    ///
    /// One batch at a time, which is what wasm does: it has no threads.
    ///
    /// # Errors
    ///
    /// When a batch cannot be read, when it holds another number of
    /// variants than its entry of the footer, and when one of its columns
    /// holds a null where every variant has a value.
    fn next_batch(&mut self) -> Result<Option<Block>> {
        // The columns to decompress are chosen when the batch is read, so a
        // consumer that asks for other fields is served from here on.
        let wanted = projection_of(self.needs, &self.columns);
        let places: Vec<usize> = wanted.iter().map(|(_, place)| *place).collect();
        loop {
            let Some(at) = self.blocks.get(self.next).copied() else {
                return Ok(None);
            };
            let index = self.next;
            let place = self.place_of(index)?;
            let bytes = self.bytes_of_the_batch(at, place)?;
            // The batches of a file are as many as the machine counts, so
            // this never saturates.
            self.next = self.next.saturating_add(1);
            let batch = batch_of(&self.schema, self.version, &places, at, bytes, place)?;
            if let Some(block) = self.block_of_a_decoded_batch(index, &batch, &wanted)? {
                return Ok(Some(block));
            }
        }
    }
}

/// One batch of a vars file, from `bytes`, with the columns of `places`
/// decompressed and the buffers of the rest walked past.
///
/// Nothing of the reader is read here: the schema and the version are what
/// the file was opened with, `at` and `place` are where the batch is, and
/// `bytes` are its bytes. That is what lets the batches of a window be
/// decoded at once on the threads of rayon.
///
/// # Errors
///
/// When the bytes of the batch are not what the file says they are, and when
/// they are compressed with zstd.
fn batch_of(
    schema: &SchemaRef,
    version: MetadataVersion,
    places: &[usize],
    at: BatchAt,
    bytes: Vec<u8>,
    place: BatchPlace,
) -> Result<RecordBatch> {
    let (Ok(metadata_len), Ok(body_len)) =
        (i32::try_from(at.metadata_len), i64::try_from(at.body_len))
    else {
        return Err(batch_of_other_bytes(
            place.batch,
            format!(
                "its footer says the batch is {metadata} and {body} bytes, which an arrow file does not hold",
                metadata = at.metadata_len,
                body = at.body_len
            ),
        ));
    };
    // arrow-rs reads the four bytes that mark a continuation and the
    // four that say how long the message is before anything else, so a
    // batch of fewer bytes than those eight is refused here.
    if bytes.len() < MESSAGE_START_BYTES {
        return Err(batch_of_other_bytes(
            place.batch,
            format!(
                "it is {found} bytes and the message of a batch of an arrow file starts with {MESSAGE_START_BYTES}",
                found = bytes.len()
            ),
        ));
    }
    // What the message of the batch says about its buffers, checked
    // against the bytes that are there before arrow-rs reads any of
    // them by their place.
    message_fits(&bytes, metadata_len, schema, place)?;
    let decoder = FileDecoder::new(Arc::clone(schema), version).with_projection(places.to_vec());
    // The net under the checks above: arrow-rs reads a length of the
    // message and panics where the bytes it points at are not there,
    // and a file that was damaged in a way those checks do not see must
    // not end the session of a user. Nothing of the reader is given to
    // arrow-rs, so what it leaves behind is the batch alone, and the
    // reader is finished after the error either way. Under wasm, where
    // a panic ends the program and unwinds nothing, this catches
    // nothing and the checks above are the whole defence.
    let read = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        decoder.read_record_batch(
            &ArrowBlock::new(0, metadata_len, body_len),
            &Buffer::from(bytes),
        )
    }));
    match read {
        Ok(Ok(Some(batch))) => Ok(batch),
        Ok(Ok(None)) => Err(batch_of_other_bytes(
            place.batch,
            "the message where it starts is not one of a batch".to_owned(),
        )),
        Ok(Err(problem)) => Err(batch_not_read(&problem, place.batch)),
        Err(_) => Err(batch_of_other_bytes(
            place.batch,
            "arrow-rs did not come back from reading it".to_owned(),
        )),
    }
}

impl<R: Read + Seek + Send> BlockReader for VarsReader<R> {
    /// The next batch of the file as a block, which holds one variant at
    /// least, and `None` when there are no more batches and at every call
    /// after that.
    ///
    /// The batches come out as they are in the file, so a file written with
    /// the default size of block is read back in blocks of that size with no
    /// `reblock`. A batch of no variants is not given: the next one is
    /// taken.
    ///
    /// # Errors
    ///
    /// When the source cannot be read, when a batch is not what the file
    /// says it is, when it holds another number of variants than its entry
    /// of the footer, when one of its columns holds a null where every
    /// variant has a value, and when the file is compressed with zstd. The
    /// batch the error happened in is lost, the blocks before it were
    /// given, and every call after it gives `None`.
    fn next_block(&mut self) -> Result<Option<Block>> {
        if self.finished {
            return Ok(None);
        }
        match self.next_batch() {
            Ok(Some(block)) => Ok(Some(block)),
            Ok(None) => {
                self.finished = true;
                Ok(None)
            }
            Err(error) => {
                self.finished = true;
                Err(error)
            }
        }
    }

    fn individuals(&self) -> &[String] {
        &self.metadata.individuals
    }

    fn ploidy(&self) -> usize {
        self.metadata.ploidy
    }

    fn chroms(&self) -> &ChromTable {
        &self.chroms
    }

    /// The fields the blocks from the next one on hold. The projection of a
    /// batch is chosen when that batch is read, so a change holds from the
    /// next block.
    fn set_needs(&mut self, needs: Needs) {
        self.needs = needs;
        // The batches that were decoded and not given were decoded with the
        // projection of the fields that were asked for before this call.
        #[cfg(not(target_family = "wasm"))]
        self.throw_the_queue_away();
    }

    /// None: a source has no filter over it.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        Vec::new()
    }
}

impl VarsReader<BufReader<File>> {
    /// The reader over the vars file at `path`, for the callers that have a
    /// path and not bytes.
    ///
    /// # Errors
    ///
    /// When the file cannot be opened, with the path in the error, and a
    /// directory at the path among them: opening a directory succeeds on
    /// macOS and on Linux and only the first read of it fails, so the
    /// reader asks the file it opened what it is. And everything
    /// [`VarsReader::new`] fails with.
    pub fn from_path(path: &Path) -> Result<Self> {
        let file = File::open(path).map_err(|problem| not_opened(path, problem))?;
        // What was opened is asked and not the path, so what happens at the
        // path between the two calls does not change the answer.
        let what_it_is = file
            .metadata()
            .map_err(|problem| not_opened(path, problem))?;
        if what_it_is.is_dir() {
            return Err(not_opened(
                path,
                std::io::Error::from_raw_os_error(A_DIRECTORY_IS_THERE),
            ));
        }
        VarsReader::new(BufReader::new(file))
    }
}

/// That the message of a batch says of its buffers and of its rows what the
/// bytes that came with it can hold.
///
/// arrow-rs reads the offsets and the lengths of that message by their
/// place, and asks the machine for the memory that the length of a
/// compressed buffer says before it decompresses it, so a file that was
/// damaged there reaches a panic inside it, and a length of
/// 144115188075855871 bytes ends the process, which nothing catches. What
/// is checked: the message is one of a batch; its rows are the variants its
/// entry of the footer gives; every buffer lies inside the body; every
/// compressed buffer says a length that the column it belongs to can hold,
/// which for the genotypes is the rows times the individuals times the
/// ploidy, and, where the schema does not give that number, one that lz4
/// can give from the bytes it holds; and no field node says more values
/// than its column can hold, which for the genotypes is one for each
/// allele, the rows times the individuals times the ploidy, and, for a
/// column popnei does not walk, what the body of the batch can hold.
///
/// # Errors
///
/// The batch could not be read, with the batch and what does not fit, and
/// the batch holds another number of variants than its entry of the footer
/// when its message says so.
fn message_fits(bytes: &[u8], metadata_len: i32, schema: &Schema, place: BatchPlace) -> Result<()> {
    let damaged = |problem: String| batch_of_other_bytes(place.batch, problem);
    // A message that starts with the mark of a continuation has its length
    // after it, and one that has not starts with that length.
    let starts_at = match bytes.get(..CONTINUATION_MARK.len()) {
        Some(start) if start == CONTINUATION_MARK => MESSAGE_START_BYTES,
        Some(_) | None => CONTINUATION_MARK.len(),
    };
    let Some(message) = bytes.get(starts_at..) else {
        return Err(damaged(format!(
            "it is {found} bytes and its message starts at the byte {starts_at}",
            found = bytes.len()
        )));
    };
    let message = root_as_message(message)
        .map_err(|problem| damaged(format!("its message is not one of arrow: {problem}")))?;
    if message.header_type() != MessageHeader::RecordBatch {
        return Err(damaged(
            "the message where it starts is not one of a batch".to_owned(),
        ));
    }
    let Some(batch) = message.header_as_record_batch() else {
        return Err(damaged(
            "its message says it is a batch and holds none".to_owned(),
        ));
    };
    let rows = batch.length();
    let found = usize::try_from(rows)
        .map_err(|_| damaged(format!("its message says it holds {rows} rows")))?;
    if found != place.num_vars {
        return Err(Error::VarsBatchNumVars {
            batch: place.batch,
            found,
            expected: place.num_vars,
        });
    }
    let body = u64::try_from(metadata_len)
        .ok()
        .and_then(|message_len| {
            u64::try_from(bytes.len())
                .ok()
                .map(|whole| whole.saturating_sub(message_len))
        })
        .unwrap_or(0);
    // A conversion that fails leaves no rows and so no bound, which
    // refuses the batch: a check that cannot work out its bound must not
    // let the bytes through.
    let rows = u64::try_from(found).unwrap_or(0);
    let compressed_with_lz4 = batch
        .compression()
        .is_some_and(|how| how.codec() == CompressionType::LZ4_FRAME);
    let holds = what_the_buffers_hold(schema.fields(), rows);
    buffers_fit(
        &batch,
        bytes,
        metadata_len,
        body,
        compressed_with_lz4,
        &holds,
        &damaged,
    )?;
    // A field node says how many values a column holds, which arrow-rs
    // turns into the length of an array. Two things bound it, and the
    // smaller one is what the node is held to.
    //
    // The column itself is the tighter of the two wherever popnei walks the
    // type, and it is what the compressed bytes of the batch cannot give: a
    // batch is decompressed before its values are counted, so a file whose
    // genotypes repeat holds more values than the batch holds bits.
    //
    // The body is the other, and it is the one that holds when the schema
    // gives none, at a column whose type popnei does not walk and at every
    // column after it. Without it those columns have no bound at all, since
    // the rows of a batch are what its entry of the footer says and no byte
    // of the file bounds them.
    let in_the_body = values_the_body_holds(body, compressed_with_lz4);
    let values = what_the_nodes_hold(schema.fields(), rows);
    for (at, node) in batch.nodes().into_iter().flatten().enumerate() {
        let length = node.length();
        let null_count = node.null_count();
        if length < 0 || null_count < 0 || null_count > length {
            return Err(damaged(format!(
                "a column of its message says it holds {length} values, {null_count} of them not there"
            )));
        }
        let at_most = match values.get(at).copied() {
            Some(NodeHolds::Values(most)) => most.min(in_the_body),
            Some(NodeHolds::NotBounded) | None => in_the_body,
        };
        if u64::try_from(length).unwrap_or(u64::MAX) > at_most {
            return Err(damaged(format!(
                "a column of its message says it holds {length} values and that column holds {at_most}"
            )));
        }
    }
    Ok(())
}

/// How many values the body of a batch can hold, whatever the schema says
/// of its columns: a value takes a bit at the very least once the batch is
/// decompressed, and an lz4 frame gives at most [`LZ4_BYTES_FOR_A_BYTE`]
/// for each byte it holds.
///
/// It is loose by a factor of eight or more for every column of a file
/// popnei writes, whose exact bound the schema gives. What it is for is the
/// columns the schema does not reach: a column of a type popnei does not
/// walk, and every column after it, which otherwise nothing bounds, because
/// the rows of a batch are what its entry of the footer says and no byte of
/// the file bounds them.
fn values_the_body_holds(body: u64, compressed_with_lz4: bool) -> u64 {
    let decompressed = if compressed_with_lz4 {
        body.saturating_mul(LZ4_BYTES_FOR_A_BYTE)
            .saturating_add(LZ4_FRAME_SLACK)
    } else {
        body
    };
    // A bit for each value.
    decompressed.saturating_mul(8)
}

/// What one buffer of a batch holds at most once it is decompressed.
///
/// arrow-rs asks the machine for the memory a compressed buffer says it
/// holds before it decompresses it, so what the column of that buffer can
/// hold is the bound that matters: under wasm a `gts` buffer that says
/// 3000000000 bytes is a trap that ends the tab, and one that says
/// 2000000000 leaves the memory of the tab grown for its life, although
/// both give an error in the end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BufferHolds {
    /// That many bytes, which the schema of the file and the rows of the
    /// batch give.
    Bytes(u64),
    /// What the schema does not say: the texts of a column of texts, the
    /// values of a list, whose count its offsets give, and every buffer of
    /// a column whose type popnei does not walk. What bounds those is what
    /// lz4 gives for the bytes the buffer holds, and the 2147483647 bytes
    /// the offsets of a column of texts address.
    WhatLz4Gives,
}

/// What each buffer of a batch of `fields` of `rows` rows holds at most, in
/// the order the IPC format lays them out: a depth first walk of the
/// columns, the buffer of the nulls of each before the rest of its own, and
/// the buffers of the values of a list after those of the list.
///
/// The walk stops at a column whose type popnei does not know the buffers
/// of, since it cannot say where the buffers of the columns after it start:
/// those are not in the list and are bounded by what lz4 gives.
fn what_the_buffers_hold(fields: &Fields, rows: u64) -> Vec<BufferHolds> {
    let mut holds = Vec::new();
    for field in fields {
        if !buffers_of_the_field(field, Some(rows), &mut holds) {
            break;
        }
    }
    holds
}

/// The buffers of one column, and of the columns inside it, after the ones
/// already in `holds`; `false` when popnei does not know the buffers of its
/// type, which leaves `holds` without them.
///
/// `rows` is `None` for the values of a list, which are as many as the
/// offsets of that list say and not as many as the batch has rows: their
/// buffers are counted, so the ones after them keep their place, and none
/// of them is bounded.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "arrow has forty types and a version of arrow-rs adds more; what popnei walks is \
              the handful named here and the ones of a fixed width, and every other one stops \
              the walk"
)]
fn buffers_of_the_field(field: &Field, rows: Option<u64>, holds: &mut Vec<BufferHolds>) -> bool {
    let bytes = |per_row: u64, and: u64| match rows {
        Some(rows) => BufferHolds::Bytes(rows.saturating_mul(per_row).saturating_add(and)),
        None => BufferHolds::WhatLz4Gives,
    };
    // The nulls of a column are a bit for each row, which arrow writes as
    // no buffer at all when there is none.
    let nulls = match rows {
        Some(rows) => BufferHolds::Bytes(bits_in_bytes(rows)),
        None => BufferHolds::WhatLz4Gives,
    };
    match field.data_type() {
        // No buffer at all: every row of such a column is a null.
        DataType::Null => true,
        DataType::Boolean => {
            holds.push(nulls);
            holds.push(nulls);
            true
        }
        DataType::Utf8 | DataType::Binary => {
            holds.push(nulls);
            // One offset of 32 bits for each row and one after the last.
            holds.push(bytes(4, 4));
            holds.push(BufferHolds::WhatLz4Gives);
            true
        }
        DataType::LargeUtf8 | DataType::LargeBinary => {
            holds.push(nulls);
            holds.push(bytes(8, 8));
            holds.push(BufferHolds::WhatLz4Gives);
            true
        }
        DataType::FixedSizeBinary(width) => {
            holds.push(nulls);
            holds.push(bytes(u64::try_from(*width).unwrap_or(0), 0));
            true
        }
        DataType::List(inside) | DataType::Map(inside, _) => {
            holds.push(nulls);
            holds.push(bytes(4, 4));
            buffers_of_the_field(inside, None, holds)
        }
        DataType::LargeList(inside) => {
            holds.push(nulls);
            holds.push(bytes(8, 8));
            buffers_of_the_field(inside, None, holds)
        }
        DataType::FixedSizeList(inside, width) => {
            holds.push(nulls);
            let values = rows.map(|rows| rows.saturating_mul(u64::try_from(*width).unwrap_or(0)));
            buffers_of_the_field(inside, values, holds)
        }
        DataType::Struct(inside) => {
            holds.push(nulls);
            for field in inside {
                if !buffers_of_the_field(field, rows, holds) {
                    return false;
                }
            }
            true
        }
        // Every type of arrow that holds a value of a fixed width: the
        // numbers, the dates, the times and the decimals. `primitive_width`
        // is what says which, and the types it gives no width for are the
        // ones whose buffers popnei does not walk.
        other => match other.primitive_width() {
            Some(width) => {
                holds.push(nulls);
                holds.push(bytes(u64::try_from(width).unwrap_or(0), 0));
                true
            }
            // A dictionary, a union, a view or a type arrow adds later:
            // popnei does not know how many buffers it takes, so the
            // buffers of it and of the columns after it are not bounded.
            None => false,
        },
    }
}

/// How many bytes that many bits take.
fn bits_in_bytes(bits: u64) -> u64 {
    bits.saturating_add(7).saturating_div(8)
}

/// How many values one column of a batch holds at most, which arrow-rs
/// turns into the length of the array it builds for that column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NodeHolds {
    /// That many values, which the schema of the file and the rows of the
    /// batch give.
    Values(u64),
    /// What the schema does not say: the values of a large list, whose
    /// count its 64 bit offsets give and nothing in the schema bounds.
    /// What bounds those is [`values_the_body_holds`].
    NotBounded,
}

/// How many values each column of a batch of `fields` of `rows` rows holds
/// at most, in the order the IPC format lays its field nodes out: one node
/// for each column, and the nodes of the columns inside a column after its
/// own, in a depth first walk.
///
/// A node with no entry in the list is the one of a column popnei does not
/// walk, or of a column after one: popnei cannot say how many nodes such a
/// column takes, so it cannot say which column the nodes after it belong
/// to. Those are bounded by [`values_the_body_holds`] as well.
fn what_the_nodes_hold(fields: &Fields, rows: u64) -> Vec<NodeHolds> {
    let mut holds = Vec::new();
    for field in fields {
        if !nodes_of_the_field(field, NodeHolds::Values(rows), &mut holds) {
            break;
        }
    }
    holds
}

/// The node of one column, which holds `values` values at most, and the
/// nodes of the columns inside it, after the ones already in `holds`;
/// `false` when popnei does not know how many nodes its type takes, which
/// leaves the columns after it out of `holds`.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "arrow has forty types and a version of arrow-rs adds more; what popnei walks is \
              the handful named here and the ones of a fixed width, and every other one stops \
              the walk"
)]
fn nodes_of_the_field(field: &Field, values: NodeHolds, holds: &mut Vec<NodeHolds>) -> bool {
    holds.push(values);
    match field.data_type() {
        DataType::Null
        | DataType::Boolean
        | DataType::Utf8
        | DataType::Binary
        | DataType::LargeUtf8
        | DataType::LargeBinary
        | DataType::FixedSizeBinary(_) => true,
        // The values of a list are as many as its offsets say, which the
        // reader has not decompressed: what bounds them is the largest
        // number those offsets hold.
        DataType::List(inside) | DataType::Map(inside, _) => {
            nodes_of_the_field(inside, NodeHolds::Values(MAX_VALUES_OF_A_LIST), holds)
        }
        // The offsets of a large list are 64 bits, so nothing in the schema
        // says how many values it holds.
        DataType::LargeList(inside) => nodes_of_the_field(inside, NodeHolds::NotBounded, holds),
        // The genotypes are one of these: every row holds the same number
        // of values, the individuals times the ploidy. A width that is not
        // a number of values, which a damaged schema gives, leaves no room
        // for a value and refuses the column.
        DataType::FixedSizeList(inside, width) => {
            let inside_values = match values {
                NodeHolds::Values(values) => {
                    NodeHolds::Values(values.saturating_mul(u64::try_from(*width).unwrap_or(0)))
                }
                NodeHolds::NotBounded => NodeHolds::NotBounded,
            };
            nodes_of_the_field(inside, inside_values, holds)
        }
        // Every column of a struct holds the rows of the struct itself.
        DataType::Struct(inside) => {
            for field in inside {
                if !nodes_of_the_field(field, values, holds) {
                    return false;
                }
            }
            true
        }
        // Every type of arrow that holds a value of a fixed width, as in
        // `buffers_of_the_field`: one node and nothing inside it. A
        // dictionary, a union, a view or a type arrow adds later takes a
        // number of nodes popnei does not know, so the columns after it are
        // not bounded.
        other => other.primitive_width().is_some(),
    }
}

/// That every buffer of the message lies inside the body of the batch, and
/// that a compressed one says a length that lz4 can give from the bytes it
/// holds.
///
/// # Errors
///
/// The batch could not be read, with what does not fit.
fn buffers_fit(
    batch: &BatchMessage<'_>,
    bytes: &[u8],
    metadata_len: i32,
    body: u64,
    compressed_with_lz4: bool,
    holds: &[BufferHolds],
    damaged: &impl Fn(String) -> Error,
) -> Result<()> {
    for (place, buffer) in batch.buffers().into_iter().flatten().enumerate() {
        let (offset, length) = (buffer.offset(), buffer.length());
        let ends_at = u64::try_from(offset)
            .ok()
            .zip(u64::try_from(length).ok())
            .and_then(|(offset, length)| offset.checked_add(length))
            .filter(|end| *end <= body);
        if ends_at.is_none() {
            return Err(damaged(format!(
                "a buffer of its message is {length} bytes at the byte {offset} of a body of {body}"
            )));
        }
        if !compressed_with_lz4 || length == 0 {
            continue;
        }
        // The bytes of a buffer start after the message of the batch, and
        // the first eight of them say how long it is once it is
        // decompressed.
        let says = starts_at_in(bytes, metadata_len, offset)
            .and_then(|at| bytes.get(at..at.saturating_add(UNCOMPRESSED_LENGTH_BYTES)))
            .and_then(|start| <[u8; UNCOMPRESSED_LENGTH_BYTES]>::try_from(start).ok())
            .map(i64::from_le_bytes);
        let Some(says) = says else {
            return Err(damaged(format!(
                "a buffer of its message is {length} bytes at the byte {offset} and the bytes that say how long it is when it is decompressed are not there"
            )));
        };
        if says == NOT_COMPRESSED || says == 0 {
            continue;
        }
        let compressed = u64::try_from(length)
            .unwrap_or(0)
            .saturating_sub(u64::try_from(UNCOMPRESSED_LENGTH_BYTES).unwrap_or(0));
        let loose = compressed
            .saturating_mul(LZ4_BYTES_FOR_A_BYTE)
            .saturating_add(LZ4_FRAME_SLACK)
            .min(MAX_COLUMN_BYTES);
        let at_most = match holds.get(place).copied() {
            Some(BufferHolds::Bytes(bytes)) => bytes,
            Some(BufferHolds::WhatLz4Gives) | None => loose,
        };
        if says < 0 || u64::try_from(says).unwrap_or(u64::MAX) > at_most {
            return Err(damaged(format!(
                "a buffer of its message holds {compressed} bytes compressed with lz4 and says it is {says} bytes decompressed, and its column holds {at_most}"
            )));
        }
    }
    Ok(())
}

/// Where the bytes of a buffer of the body are in the bytes of the batch:
/// after the message, at the offset the buffer gives.
fn starts_at_in(bytes: &[u8], metadata_len: i32, offset: i64) -> Option<usize> {
    let at = usize::try_from(metadata_len)
        .ok()?
        .checked_add(usize::try_from(offset).ok()?)?;
    (at < bytes.len()).then_some(at)
}

/// One batch of a vars file as a block: the columns of `wanted`, which are
/// the ones that were decompressed and are in the batch in that order, and
/// the individuals and the ploidy of the `popnei` key.
///
/// The names of the chromosomes are interned into `chroms` as they are read,
/// so a number means the order of first appearance among the variants that
/// were given, as in the VCF reader, and the table is looked up only when the
/// name differs from that of the variant before.
///
/// # Errors
///
/// When a column holds a null where every variant has a value, the
/// chromosome, the position, the alleles and the genotypes; when a column of
/// the batch is not what the schema of the file says it is, which is a file
/// that was damaged; and when the machine does not give the memory of a
/// column of the block.
fn block_of_the_batch(
    batch: &RecordBatch,
    wanted: &[(VarsColumn, usize)],
    metadata: &VarsMetadata,
    chroms: &mut ChromTable,
    place: BatchPlace,
) -> Result<Block> {
    let num_vars = batch.num_rows();
    let mut block = Block {
        num_vars,
        num_individuals: metadata.individuals.len(),
        ploidy: metadata.ploidy,
        gts: Vec::new(),
        chrom: None,
        pos: None,
        id: None,
        alleles: None,
        qual: None,
    };
    // The texts of the alleles of one variant, written over for the next
    // one: nothing is allocated for each variant of the block.
    let mut texts: Vec<String> = Vec::new();
    for (index, (column, _)) in wanted.iter().enumerate() {
        let Some(array) = batch.columns().get(index) else {
            return Err(column_not_read(
                *column,
                place.batch,
                "is not among the columns that were read",
            ));
        };
        match column {
            VarsColumn::Chrom => {
                block.chrom = Some(chrom_numbers(array, num_vars, metadata, chroms, place)?);
            }
            VarsColumn::Pos => block.pos = Some(positions(array, num_vars, metadata, place)?),
            VarsColumn::Id => block.id = Some(ids(array, num_vars, metadata, place)?),
            VarsColumn::Alleles => {
                block.alleles = Some(alleles_of_the_batch(
                    array, num_vars, &mut texts, metadata, place,
                )?);
            }
            VarsColumn::Qual => block.qual = Some(qualities(array, num_vars, metadata, place)?),
            VarsColumn::Gts => block.gts = genotypes(array, num_vars, metadata, place)?,
        }
    }
    Ok(block)
}

/// The number of the chromosome of each variant of the batch, in `chroms`.
///
/// # Errors
///
/// When a variant has no chromosome, and when the column is not one of
/// texts.
fn chrom_numbers(
    array: &ArrayRef,
    num_vars: usize,
    metadata: &VarsMetadata,
    chroms: &mut ChromTable,
    place: BatchPlace,
) -> Result<Vec<u32>> {
    let column = column_of::<StringArray>(array, VarsColumn::Chrom, place)?;
    let mut numbers = reserved_column(num_vars, metadata)?;
    // The name of the variant before and its number: the variants of one
    // chromosome come one after another in a sorted file, so the table is
    // looked up once for the run.
    let mut before: Option<(&str, u32)> = None;
    for (row, name) in column.iter().enumerate() {
        let Some(name) = name else {
            return Err(null_value(VarsColumn::Chrom, place, row));
        };
        let number = match before {
            Some((seen, number)) if seen == name => number,
            Some(_) | None => {
                let number = chroms.intern(name);
                before = Some((name, number));
                number
            }
        };
        numbers.push(number);
    }
    Ok(numbers)
}

/// The position of each variant of the batch, copied out of the column.
///
/// # Errors
///
/// When a variant has no position, and when the column is not one of
/// unsigned 64 bit numbers.
fn positions(
    array: &ArrayRef,
    num_vars: usize,
    metadata: &VarsMetadata,
    place: BatchPlace,
) -> Result<Vec<u64>> {
    let column = column_of::<UInt64Array>(array, VarsColumn::Pos, place)?;
    if let Some(row) = first_null(column) {
        return Err(null_value(VarsColumn::Pos, place, row));
    }
    let Some(values) = column.values().get(..num_vars) else {
        return Err(column_not_read(
            VarsColumn::Pos,
            place.batch,
            "holds fewer values than the batch has variants",
        ));
    };
    let mut positions = reserved_column(num_vars, metadata)?;
    positions.extend_from_slice(values);
    Ok(positions)
}

/// The id of each variant of the batch, the empty id where the file holds a
/// null.
///
/// # Errors
///
/// When the column is not one of texts, and when the machine does not give
/// its memory.
fn ids(
    array: &ArrayRef,
    num_vars: usize,
    metadata: &VarsMetadata,
    place: BatchPlace,
) -> Result<Vec<String>> {
    let column = column_of::<StringArray>(array, VarsColumn::Id, place)?;
    let mut ids: Vec<String> = reserved_column(num_vars, metadata)?;
    for id in column.iter() {
        ids.push(id.unwrap_or("").to_owned());
    }
    Ok(ids)
}

/// The quality of each variant of the batch, a NaN where the file holds a
/// null, which is what a block holds for a variant with no quality.
///
/// # Errors
///
/// When the column is not one of 32 bit floats, and when the machine does
/// not give its memory.
fn qualities(
    array: &ArrayRef,
    num_vars: usize,
    metadata: &VarsMetadata,
    place: BatchPlace,
) -> Result<Vec<f32>> {
    let column = column_of::<Float32Array>(array, VarsColumn::Qual, place)?;
    let mut qualities = reserved_column(num_vars, metadata)?;
    for (row, quality) in column.iter().enumerate() {
        let Some(quality) = quality else {
            // The variant has no quality, which a block holds as a NaN.
            qualities.push(f32::NAN);
            continue;
        };
        // A value that is not finite: the NaN of the column of a block is
        // what says that a variant has no quality, and an infinite quality
        // is a probability of no variant of 0.
        if !quality.is_finite() {
            return Err(Error::VarsQualityNotFinite {
                found: quality,
                var: place.vars_before.saturating_add(counted_from_one(row)),
            });
        }
        qualities.push(quality);
    }
    Ok(qualities)
}

/// The alleles of each variant of the batch, the reference one first.
///
/// `texts` is the buffer the texts of one variant are read into and written
/// over for the next one.
///
/// # Errors
///
/// When a variant has no alleles or one of its alleles is a null, when the
/// column is not a list of texts, and when the machine does not give the
/// memory of the column of the block.
fn alleles_of_the_batch(
    array: &ArrayRef,
    num_vars: usize,
    texts: &mut Vec<String>,
    metadata: &VarsMetadata,
    place: BatchPlace,
) -> Result<AllelesColumn> {
    let column = column_of::<ListArray>(array, VarsColumn::Alleles, place)?;
    let values = column
        .values()
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| {
            column_not_read(
                VarsColumn::Alleles,
                place.batch,
                "is not a list of the texts of the alleles",
            )
        })?;
    let mut alleles =
        AllelesColumn::with_num_vars(num_vars).map_err(|_| block_too_large(num_vars, metadata))?;
    // One window for each variant: where its alleles start in the texts of
    // the column and where they end.
    for (row, window) in column.offsets().windows(2).enumerate() {
        let (Some(start), Some(end)) = (window.first(), window.get(1)) else {
            break;
        };
        if column.is_null(row) {
            return Err(null_value(VarsColumn::Alleles, place, row));
        }
        let ends = usize::try_from(*start)
            .ok()
            .zip(usize::try_from(*end).ok())
            .filter(|(start, end)| end >= start && *end <= values.len());
        let Some((start, end)) = ends else {
            return Err(column_not_read(
                VarsColumn::Alleles,
                place.batch,
                "says its alleles are at bytes that are not in it",
            ));
        };
        let num_alleles = end.saturating_sub(start);
        for (allele, index) in (start..end).enumerate() {
            if values.is_null(index) {
                return Err(null_value(VarsColumn::Alleles, place, row));
            }
            // The index is one of the column: the window it comes from was
            // checked against the length of the texts above.
            let text = values.value(index);
            match texts.get_mut(allele) {
                Some(held) => {
                    held.clear();
                    held.push_str(text);
                }
                None => texts.push(text.to_owned()),
            }
        }
        alleles.push(texts.get(..num_alleles).unwrap_or(&[]));
    }
    Ok(alleles)
}

/// The genotypes of the batch, the alleles of every variant one after
/// another, copied out of the buffer arrow decompressed into the vector of
/// the block.
///
/// # Errors
///
/// When a variant has no genotypes or one of its alleles is a null, when the
/// column is not a fixed size list of signed bytes or holds another number
/// of alleles for each variant than the `popnei` key of the file gives, when
/// an allele of it is below [`MISSING_ALLELE`], and when the machine does
/// not give the memory of the genotypes of the block.
fn genotypes(
    array: &ArrayRef,
    num_vars: usize,
    metadata: &VarsMetadata,
    place: BatchPlace,
) -> Result<Vec<i8>> {
    let column = column_of::<FixedSizeListArray>(array, VarsColumn::Gts, place)?;
    if let Some(row) = first_null(column) {
        return Err(null_value(VarsColumn::Gts, place, row));
    }
    let num_individuals = metadata.individuals.len();
    let ploidy = metadata.ploidy;
    let expected = num_individuals
        .checked_mul(ploidy)
        .ok_or_else(|| block_too_large(num_vars, metadata))?;
    let found = usize::try_from(column.value_length()).unwrap_or(usize::MAX);
    // The width of the column was checked against the `popnei` key when the
    // file was opened, and it is what turns the flat buffer into variants.
    if found != expected {
        return Err(Error::VarsGtsWidth {
            found,
            expected,
            num_individuals,
            ploidy,
        });
    }
    let values = column
        .values()
        .as_any()
        .downcast_ref::<Int8Array>()
        .ok_or_else(|| {
            column_not_read(
                VarsColumn::Gts,
                place.batch,
                "does not hold the alleles as signed bytes",
            )
        })?;
    if let Some(index) = first_null(values) {
        return Err(null_value(
            VarsColumn::Gts,
            place,
            index.checked_div(found).unwrap_or(0),
        ));
    }
    let wanted = num_vars
        .checked_mul(found)
        .ok_or_else(|| block_too_large(num_vars, metadata))?;
    let Some(alleles) = values.values().get(..wanted) else {
        return Err(column_not_read(
            VarsColumn::Gts,
            place.batch,
            "holds fewer alleles than the variants of the batch",
        ));
    };
    // No reader of popnei gives an allele below the missing one, which the
    // owner decided on 21 September 2026 and `docs/specs/variant.md` says
    // under "An allele that no reader gives": the two counts of one variant
    // refuse such an allele as a defect of the reader, and a pass that
    // counts nothing would put it in the array of a user. The alleles of
    // the file are signed bytes, so a damaged one and a file of another
    // program can say -2.
    if let Some((at, allele)) = first_allele_below_the_missing_one(alleles) {
        return Err(Error::VarsAlleleBelowMissing {
            found: allele,
            var: place
                .vars_before
                .saturating_add(counted_from_one(at.checked_div(found).unwrap_or(0))),
        });
    }
    let mut genotypes = Vec::new();
    genotypes
        .try_reserve_exact(wanted)
        .map_err(|_| block_too_large(num_vars, metadata))?;
    genotypes.extend_from_slice(alleles);
    Ok(genotypes)
}

/// The first allele of `alleles` below [`MISSING_ALLELE`], with its place
/// in the slice, and `None` when every one of them is an allele.
///
/// Whether there is one is the smallest of them in one pass, which the
/// compiler reduces over the lanes of a vector register; a comparison
/// written for each allele on its own does not vectorise, and this reads
/// every allele of every batch. Which one it is and where it is come from
/// one walk of the same slice, which only a batch that is refused pays
/// for, so the allele the error names is the allele at the place it names,
/// whatever either pass is changed into.
fn first_allele_below_the_missing_one(alleles: &[i8]) -> Option<(usize, i8)> {
    if alleles.iter().copied().fold(i8::MAX, i8::min) >= MISSING_ALLELE {
        return None;
    }
    alleles
        .iter()
        .copied()
        .enumerate()
        .find(|(_, allele)| *allele < MISSING_ALLELE)
}

/// One column of the batch as the array it holds.
///
/// # Errors
///
/// When it is not that array, which the schema of the file said it is, so
/// the file was damaged after it was written.
fn column_of<A: Array + 'static>(
    array: &ArrayRef,
    column: VarsColumn,
    place: BatchPlace,
) -> Result<&A> {
    array.as_any().downcast_ref::<A>().ok_or_else(|| {
        column_not_read(
            column,
            place.batch,
            "is not of the arrow type that the schema of the file gives it",
        )
    })
}

/// The first row of the column with no value, or `None` when every row has
/// one.
fn first_null(column: &dyn Array) -> Option<usize> {
    let nulls = column.nulls()?;
    nulls.iter().position(|there| !there)
}

/// The memory of one column of a block of `num_vars` variants, empty.
///
/// The memory is asked for with `try_reserve_exact`, which gives it back as
/// an error: `Vec::with_capacity` ends the process when the machine has not
/// the memory, and a batch of a file that a caller of popnei wrote reaches
/// it.
///
/// # Errors
///
/// When the machine does not give it.
fn reserved_column<T>(num_vars: usize, metadata: &VarsMetadata) -> Result<Vec<T>> {
    let mut column = Vec::new();
    column
        .try_reserve_exact(num_vars)
        .map_err(|_| block_too_large(num_vars, metadata))?;
    Ok(column)
}

/// The error of a block of a batch that the machine does not give the memory
/// for.
///
/// The blocks of a vars file are its batches, built whole whatever size the
/// caller asked its blocks to be, so the size in the message is the one the
/// file fixed and the way out is to write the file again with a smaller
/// `num_vars_per_block`.
fn block_too_large(num_vars: usize, metadata: &VarsMetadata) -> Error {
    Error::BlockTooLarge {
        num_vars_per_block: num_vars,
        num_individuals: metadata.individuals.len(),
        ploidy: metadata.ploidy,
        size: BlockSize::FixedByAFile,
    }
}

/// The error of a null in a column where every variant has a value, with the
/// column and the variant, counted from 1 over the whole file.
fn null_value(column: VarsColumn, place: BatchPlace, row: usize) -> Error {
    Error::VarsNullValue {
        column: column.name(),
        var: place.vars_before.saturating_add(counted_from_one(row)),
    }
}

/// That place of a file, counted from 1, as the messages count.
fn counted_from_one(place: usize) -> u64 {
    u64::try_from(place).unwrap_or(u64::MAX).saturating_add(1)
}

/// The error of a column of a batch that is not what the schema of the file
/// says it is, which is a file that was damaged after it was written.
fn column_not_read(column: VarsColumn, batch: u64, problem: &str) -> Error {
    batch_of_other_bytes(
        batch,
        format!("its `{name}` column {problem}", name = column.name()),
    )
}

/// The error of a batch whose bytes are not what the file says they are.
fn batch_of_other_bytes(batch: u64, problem: String) -> Error {
    Error::VarsBatchNotRead { batch, problem }
}

/// What arrow-rs said about a batch it could not read, as the error of the
/// crate.
///
/// A file whose buffers are compressed with zstd is the case of its own that
/// says so: no build of popnei carries the zstd crate, and arrow-rs asks for
/// it with an error that says the feature is off, which is told from every
/// other by its kind and by the name of the compression in it. The test on
/// `tests/reference/vars/zstd.vars` is what holds this when arrow-rs changes
/// what it says.
fn batch_not_read(problem: &ArrowError, batch: u64) -> Error {
    if let ArrowError::InvalidArgumentError(said) = problem
        && said.contains("zstd")
    {
        return Error::VarsZstd;
    }
    batch_of_other_bytes(batch, problem.to_string())
}

/// The file at `path` could not be opened, with the path, which a binding
/// crate builds the exception of its language with.
fn not_opened(path: &Path, problem: std::io::Error) -> Error {
    Error::FileNotOpened {
        path: path.to_path_buf(),
        source: problem,
    }
}

/// That the source starts with the six bytes every arrow IPC file starts
/// with.
///
/// # Errors
///
/// The source is not a vars file when it does not, which the bytes of a VCF
/// are; a source that does and ends before what it says it holds is a file
/// that was cut short, which is another error and another exception.
fn starts_as_an_arrow_file<R: Read + Seek>(source: &mut R, file_len: u64) -> Result<()> {
    let not_one = || {
        not_a_vars_file(format!(
            "it does not start with the {bytes} bytes `ARROW1` that an arrow IPC file starts with",
            bytes = ARROW_MAGIC.len()
        ))
    };
    if file_len < magic_bytes() {
        return Err(not_one());
    }
    let start = bytes_at(source, 0, ARROW_MAGIC.len())?;
    if start != ARROW_MAGIC {
        return Err(not_one());
    }
    Ok(())
}

/// The bytes of the footer of the arrow file, which its last four bytes
/// before the mark of its end say the length of.
///
/// # Errors
///
/// The file was cut short when it does not end with the mark an arrow file
/// ends with, and when the footer it says it has does not fit in it: both
/// are a file whose bytes ran out after its header was written.
fn footer_of<R: Read + Seek>(source: &mut R, file_len: u64) -> Result<Vec<u8>> {
    let Some(trailer_at) = file_len.checked_sub(TRAILER_BYTES) else {
        return Err(cut_short(format!(
            "the file is {file_len} bytes and the {TRAILER_BYTES} of the trailer that says where \
             its footer is are not in it"
        )));
    };
    let trailer = bytes_at(source, trailer_at, usize_of(TRAILER_BYTES)?)?;
    let (Some(said), Some(magic)) = (trailer.get(..4), trailer.get(4..)) else {
        return Err(cut_short(format!(
            "the trailer of the file is {TRAILER_BYTES} bytes and fewer were read"
        )));
    };
    if magic != ARROW_MAGIC {
        return Err(cut_short(
            "it does not end with the six bytes `ARROW1` that come after the footer of an arrow \
             file, so the bytes after its last batch are missing"
                .to_owned(),
        ));
    }
    let footer_len = <[u8; 4]>::try_from(said).map_or(-1, i32::from_le_bytes);
    let Some(footer_at) = u64::try_from(footer_len)
        .ok()
        .and_then(|len| trailer_at.checked_sub(len))
        .filter(|at| *at >= magic_bytes())
    else {
        return Err(cut_short(format!(
            "the file says its footer is {footer_len} bytes and it holds {trailer_at} before the \
             four bytes of that number"
        )));
    };
    let len = trailer_at.saturating_sub(footer_at);
    bytes_at(source, footer_at, usize_of(len)?)
}

/// The columns of the file, from its footer.
///
/// # Errors
///
/// The source is not a vars file when its footer says nothing of the
/// columns, when they are not ones arrow reads, and when the file was
/// written on a machine whose bytes are in the other order, which is what
/// the columns of every batch of it would be read in.
fn schema_of_the_footer(footer: &Footer<'_>) -> Result<Schema> {
    let Some(columns) = footer.schema() else {
        return Err(not_a_vars_file(
            "its footer says nothing of the columns of the file".to_owned(),
        ));
    };
    if !columns.endianness().equals_to_target_endianness() {
        return Err(not_a_vars_file(
            "it was written on a machine that holds the bytes of a number in the other order"
                .to_owned(),
        ));
    }
    try_fb_to_schema(columns).map_err(|problem| {
        not_a_vars_file(format!(
            "the columns its footer gives are not ones of an arrow file: {problem}"
        ))
    })
}

/// What the `popnei` key of the schema says, with the version and the names
/// of the individuals checked.
///
/// # Errors
///
/// The source is not a vars file when its schema has no `popnei` key and
/// for everything [`metadata_from_json`] refuses. The version is refused
/// when its first part is not the one popnei reads, with the whole version
/// in the message. And two individuals of one name are refused, as they are
/// for the VCF reader.
fn metadata_of_the_schema(schema: &Schema) -> Result<VarsMetadata> {
    let Some(value) = schema.metadata().get(POPNEI_KEY) else {
        return Err(not_a_vars_file(format!(
            "its schema has no `{POPNEI_KEY}` key"
        )));
    };
    let metadata = metadata_from_json(value)?;
    let major = metadata.format_version.split('.').next().unwrap_or("");
    if major != FORMAT_VERSION_READ {
        return Err(Error::VarsFormatVersion {
            found: metadata.format_version.clone(),
            read: FORMAT_VERSION_READ,
        });
    }
    // A file of no individual holds the genotypes of nobody, and every
    // source of popnei has one individual at least.
    if metadata.individuals.is_empty() {
        return Err(Error::VarsFileOfNoGenotypes {
            num_individuals: 0,
            ploidy: metadata.ploidy,
        });
    }
    let mut seen: HashSet<&str> = HashSet::with_capacity(metadata.individuals.len());
    for name in &metadata.individuals {
        if !seen.insert(name.as_str()) {
            return Err(Error::VarsIndividualTwice { name: name.clone() });
        }
    }
    Ok(metadata)
}

/// Where each batch of the file is, from its footer.
///
/// # Errors
///
/// The source is not a vars file when its footer says nothing of its
/// batches or puts one at a byte that is not one of a file, and it was cut
/// short when a batch it names ends past the end of the file.
fn batches_of_the_footer(footer: &Footer<'_>, file_len: u64) -> Result<Vec<BatchAt>> {
    let Some(blocks) = footer.recordBatches() else {
        return Err(not_a_vars_file(
            "its footer says nothing of the batches of the file".to_owned(),
        ));
    };
    blocks
        .iter()
        .map(|block| batch_at(block, file_len))
        .collect()
}

/// Where one batch is, from one entry of the footer.
///
/// # Errors
///
/// As [`batches_of_the_footer`], for one batch.
fn batch_at(block: &ArrowBlock, file_len: u64) -> Result<BatchAt> {
    let (Ok(offset), Ok(metadata_len), Ok(body_len)) = (
        u64::try_from(block.offset()),
        u64::try_from(block.metaDataLength()),
        u64::try_from(block.bodyLength()),
    ) else {
        return Err(not_a_vars_file(
            "its footer puts a batch at a byte of the file that is not one".to_owned(),
        ));
    };
    let ends_at = metadata_len
        .checked_add(body_len)
        .and_then(|bytes| offset.checked_add(bytes));
    if ends_at.is_none_or(|end| end > file_len) {
        return Err(cut_short(format!(
            "its footer says a batch of {metadata_len} and {body_len} bytes starts at the byte \
             {offset}, and the file is {file_len} bytes"
        )));
    }
    Ok(BatchAt {
        offset,
        metadata_len,
        body_len,
    })
}

/// What the `popnei_batches` key of the footer says of each batch, checked
/// against the batches the file has.
///
/// # Errors
///
/// The source is not a vars file when its footer has no `popnei_batches`
/// key and for everything [`batches_from_json`] refuses, and the entries of
/// that key are refused when they are not as many as the batches, since no
/// entry could then be trusted to be that of its batch.
fn batch_info_of_the_footer(footer: &Footer<'_>, num_batches: usize) -> Result<Vec<BatchInfo>> {
    let no_key = || not_a_vars_file(format!("its footer has no `{POPNEI_BATCHES_KEY}` key"));
    let mut value = None;
    for entry in footer.custom_metadata().into_iter().flatten() {
        if entry.key() == Some(POPNEI_BATCHES_KEY) {
            value = entry.value();
        }
    }
    let Some(value) = value else {
        return Err(no_key());
    };
    let batches = batches_from_json(value)?;
    if batches.len() != num_batches {
        return Err(Error::VarsBatchesDoNotMatch {
            found: batches.len(),
            expected: num_batches,
        });
    }
    Ok(batches)
}

/// The variants of the whole file, the sum of those of its batches.
///
/// # Errors
///
/// The file holds more variants than this machine counts, which under wasm,
/// where a `usize` is 32 bits, is 4295 million: the error is that of a block
/// the machine cannot hold, of the batch that carried the count past it.
fn num_vars_of_the_file(batches: &[BatchInfo], metadata: &VarsMetadata) -> Result<usize> {
    let mut num_vars: usize = 0;
    for batch in batches {
        num_vars = num_vars
            .checked_add(batch.num_vars)
            .ok_or_else(|| block_too_large(batch.num_vars, metadata))?;
    }
    Ok(num_vars)
}

/// `len` bytes of the source from `offset`, which the length of the file
/// says are there.
///
/// # Errors
///
/// The file was cut short when the bytes run out although the length of the
/// file said they were there, which is a file that changed while it was
/// being read, and [`Error::Io`] when the source fails or the machine does
/// not give the memory of the bytes.
fn bytes_at<R: Read + Seek>(source: &mut R, offset: u64, len: usize) -> Result<Vec<u8>> {
    source.seek(SeekFrom::Start(offset))?;
    // The length comes from the file, so the room for it is asked for and
    // not taken: `vec![0; len]` ends the process when the machine has not
    // the memory, and the VCF reader asks the same way for the bytes it
    // reads.
    let mut bytes: Vec<u8> = Vec::new();
    bytes.try_reserve_exact(len).map_err(|_| {
        Error::Io(std::io::Error::new(
            ErrorKind::OutOfMemory,
            format!("the {len} bytes of a part of the vars file were not given"),
        ))
    })?;
    bytes.resize(len, 0);
    source.read_exact(&mut bytes).map_err(|problem| {
        if problem.kind() == ErrorKind::UnexpectedEof {
            cut_short(format!(
                "{len} bytes were asked for at the byte {offset} of the file and the bytes ran out"
            ))
        } else {
            Error::Io(problem)
        }
    })?;
    Ok(bytes)
}

/// How many bytes an arrow file starts with, as the machine counts bytes.
fn magic_bytes() -> u64 {
    u64::try_from(ARROW_MAGIC.len()).unwrap_or(u64::MAX)
}

/// A length of the trailer or of the footer of the file, as a length this
/// machine can hold at once.
///
/// # Errors
///
/// When the length is more than a `usize`, which under wasm, where a `usize`
/// is 32 bits, is 4 GB: the file is more than this build of popnei holds. An
/// arrow file keeps the length of its footer in 32 bits, so no footer
/// reaches it; the length of a batch, which the footer keeps in 64, is
/// checked where the batch is read.
fn usize_of(len: u64) -> Result<usize> {
    usize::try_from(len).map_err(|_| {
        Error::Io(std::io::Error::new(
            ErrorKind::OutOfMemory,
            format!(
                "a part of the vars file is {len} bytes, which this machine does not hold at once"
            ),
        ))
    })
}

/// The error of a vars file that starts as an arrow file and ends before
/// what it says it holds.
fn cut_short(problem: String) -> Error {
    Error::VarsFileCutShort { problem }
}

#[cfg(test)]
mod tests {
    use std::cell::RefCell;
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufReader, Cursor, ErrorKind, Write};
    use std::path::{Path, PathBuf};
    use std::rc::Rc;
    use std::sync::{Arc, Mutex};
    use std::thread::ThreadId;

    use arrow_array::builder::{Int8Builder, ListBuilder, StringBuilder};
    use arrow_array::cast::AsArray;
    use arrow_array::types::{Float32Type, Int8Type, Int32Type, UInt64Type};
    use arrow_array::{
        Array, ArrayRef, DictionaryArray, FixedSizeListArray, Float32Array, Float64Array,
        Int8Array, Int32Array, ListArray, RecordBatch, StringArray, UInt64Array,
    };
    use arrow_buffer::{NullBuffer, ScalarBuffer};
    use arrow_ipc::reader::FileReader;
    use arrow_ipc::writer::FileWriter;
    use arrow_ipc::{root_as_footer, root_as_message};
    use arrow_schema::{DataType, Field, Fields, Schema};

    use super::{
        ALLELES_COLUMN, BatchAt, BatchInfo, BatchPlace, CHROM_COLUMN, CONTINUATION_MARK,
        FORMAT_VERSION, FORMAT_VERSION_READ, GTS_COLUMN, ITEM_FIELD, MAX_COLUMN_BYTES,
        MAX_VALUES_OF_A_LIST, MESSAGE_START_BYTES, NOT_COMPRESSED, NodeHolds, POPNEI_BATCHES_KEY,
        POPNEI_KEY, POS_COLUMN, QUAL_COLUMN, Region, UNCOMPRESSED_LENGTH_BYTES, VarsColumn,
        VarsColumns, VarsMetadata, VarsReader, VarsWriter, alleles_column, batches_as_json,
        batches_from_json, batches_of_the_footer, block_of_the_batch, block_too_large,
        chrom_column, counted_from_one, footer_of, gts_field, id_column, metadata_as_json,
        metadata_from_json, num_vars_of_the_file, projection_of, schema_of, what_the_nodes_hold,
        write_vars,
    };
    use crate::block::{AllelesColumn, Block, BlockReader, BlockSize};
    use crate::error::{Error, Result};
    use crate::filters::FilteringStats;
    use crate::io::vcf::{VcfOptions, VcfReader};
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
    ///
    /// They start on the chromosome that comes second in the alphabet, so
    /// that the regions of the batch, which are in the order in which the
    /// chromosomes first appear, are told apart from the same regions in
    /// any other order.
    const NOT_SORTED: [Row; 4] = [
        Row {
            chrom: "chr2",
            pos: 300,
            id: "",
            alleles: &["A", "T"],
            qual: None,
            gts: &[0, 0, 0, 1, 1, 1],
        },
        Row {
            chrom: "chr1",
            pos: 50,
            id: "",
            alleles: &["A", "T"],
            qual: None,
            gts: &[0, 0, 0, 1, 1, 1],
        },
        Row {
            chrom: "chr2",
            pos: 100,
            id: "",
            alleles: &["A", "T"],
            qual: None,
            gts: &[0, 0, 0, 1, 1, 1],
        },
        Row {
            chrom: "chr1",
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
        /// What it was last asked to fill, which starts as nothing: a
        /// reader is asked for its fields before it is read, and the test
        /// that holds this sees whether it was asked at all. It is shared
        /// so that the test reads it after the reader was moved into
        /// `write_vars`.
        asked_for: Arc<Mutex<Needs>>,
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
                asked_for: Arc::new(Mutex::new(Needs::empty())),
            }
        }

        /// What the reader was asked to fill, which the test keeps while
        /// the reader itself is given away.
        fn asked_for(&self) -> Arc<Mutex<Needs>> {
            Arc::clone(&self.asked_for)
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
            if let Ok(mut asked_for) = self.asked_for.lock() {
                *asked_for = needs;
            }
        }

        fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
            Vec::new()
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

    /// The arrow type of the `alleles` column, a list of texts that holds
    /// no null.
    fn alleles_type() -> DataType {
        DataType::List(Arc::new(Field::new_list_field(DataType::Utf8, false)))
    }

    /// The arrow type of the `gts` column of a file of that many alleles
    /// for each variant, whose alleles hold no null.
    fn gts_type(alleles_per_var: i32) -> DataType {
        DataType::FixedSizeList(
            Arc::new(Field::new_list_field(DataType::Int8, false)),
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
        write_vars(reader, Vec::new(), Some(num_vars_per_block))
            .expect("the file was written")
            .0
    }

    /// `write_vars` asks its reader for every field, so that a file written
    /// from a VCF holds its six columns, the ids, the alleles and the
    /// qualities among them, whether or not the user will read them, and
    /// can stand in for the VCF in any later analysis.
    #[test]
    fn write_vars_asks_its_reader_for_every_field() {
        let mut chroms = ChromTable::new();
        let block = cases_block(&mut chroms);
        let reader = GivenBlocks::of(vec![block], chroms);
        let asked_for = reader.asked_for();
        assert_eq!(
            *asked_for.lock().expect("what the reader was asked for"),
            Needs::empty(),
            "the reader was asked for fields before the call"
        );

        write_vars(reader, Vec::new(), Some(3)).expect("the file was written");

        assert_eq!(
            *asked_for.lock().expect("what the reader was asked for"),
            Needs::ALL
        );
    }

    /// `write_vars` says how many variants it wrote, which is what a Python
    /// or a TypeScript user reads as the variants of the pass: the four of
    /// `cases.vcf`, the 500 of `many.vcf` whatever the size of the batches,
    /// and none for a source that has no variants.
    #[test]
    fn write_vars_says_how_many_variants_it_wrote() {
        let mut chroms = ChromTable::new();
        let block = cases_block(&mut chroms);
        let reader = GivenBlocks::of(vec![block], chroms);
        let (_, num_vars) = write_vars(reader, Vec::new(), Some(3)).expect("the file");
        assert_eq!(num_vars, 4);

        let (_, num_vars) =
            write_vars(many_vcf_reader(None), Vec::new(), Some(100)).expect("the file");
        assert_eq!(num_vars, 500);
        let (_, num_vars) = write_vars(many_vcf_reader(None), Vec::new(), None).expect("the file");
        assert_eq!(num_vars, 500);

        let of_no_variants = GivenBlocks::of(Vec::new(), ChromTable::new());
        let (_, num_vars) = write_vars(of_no_variants, Vec::new(), Some(3)).expect("the file");
        assert_eq!(num_vars, 0);
    }

    /// A reader that is borrowed and not taken is what `write_vars` reads
    /// when its caller keeps the chain of the pass, which is how a binding
    /// crate reads the counts of the filters of that pass when the call
    /// returns, as `docs/specs/filters.md` has it: the file is the one the
    /// same reader written gives.
    #[test]
    fn write_vars_writes_the_same_file_from_a_reader_it_borrows() {
        let mut chroms = ChromTable::new();
        let block = cases_block(&mut chroms);
        let mut reader = GivenBlocks::of(vec![block], chroms);

        let (bytes, num_vars) =
            write_vars(&mut reader, Vec::new(), Some(3)).expect("the file was written");

        assert_eq!(num_vars, 4);
        assert_eq!(bytes, cases_written_in_batches_of(3));
        assert!(reader.filtering_stats().is_empty());
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
        let (bytes, _) = write_vars(reader, Vec::new(), Some(4)).expect("the file was written");
        let file = file_read(bytes);

        assert_eq!(num_rows_of(&file.rows), [4]);
        assert_eq!(
            chroms_of(&file.rows[0]),
            ["chr2", "chr1", "chr2", "chr1"],
            "the chromosome of every variant is written as its name"
        );
        // `chr2` is the first region because it is where the first variant
        // of the block is, which is not the order of the alphabet.
        assert_eq!(
            file.batches,
            vec![BatchInfo {
                num_vars: 4,
                regions: vec![
                    Region {
                        chrom: "chr2".to_owned(),
                        min_pos: 100,
                        max_pos: 300,
                    },
                    Region {
                        chrom: "chr1".to_owned(),
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
        let (bytes, _) = write_vars(reader, Vec::new(), Some(3)).expect("the file was written");
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
        let (bytes, _) = write_vars(reader, Vec::new(), Some(3)).expect("the file was written");
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
        let (bytes, _) = write_vars(reader, Vec::new(), Some(2)).expect("the file was written");
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
    ///
    /// It is made at the columns of a batch and not at the column of the
    /// genotypes alone, because what the writer is given is a block: a copy
    /// of the vector on the way from the block to that column is what this
    /// finds.
    #[test]
    fn the_genotypes_of_a_block_become_the_buffer_of_the_gts_column_with_no_copy() {
        let mut chroms = ChromTable::new();
        let block = cases_block(&mut chroms);
        let address = block.gts.as_ptr().addr();
        let fields = block.fields();
        let writer: VarsWriter<Vec<u8>> =
            VarsWriter::new(Vec::new(), &cases_individuals(), 2, 3).expect("the writer");

        let (arrays, _) = writer
            .arrays_of(block, &chroms, fields)
            .expect("the columns of the batch");

        // The genotypes are the last of the six columns of the file.
        let column = arrays.last().expect("the columns of the batch");
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
        let (bytes, _) = write_vars(reader, Vec::new(), None).expect("the file was written");

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

    /// The `popnei` key of a vars file names one individual at least, a
    /// ploidy of 1 at least and batches of 1 variant at least, and
    /// `metadata_from_json` of this module refuses a key that says less. A
    /// writer built with any of the three would write a file that no
    /// reader of popnei opens, and a file of no individual would have no
    /// `gts` column, which every vars file has.
    #[test]
    fn a_writer_of_no_genotype_or_of_batches_of_no_variant_is_refused() {
        let individuals = cases_individuals();
        let writer: Result<VarsWriter<Vec<u8>>> = VarsWriter::new(Vec::new(), &individuals, 2, 3);
        assert!(writer.is_ok());

        let error = match VarsWriter::new(Vec::new(), &[], 2, 3) {
            Ok(_) => panic!("a writer of no individual was built"),
            Err(error) => error,
        };
        let Error::VarsFileOfNoGenotypes {
            num_individuals,
            ploidy,
        } = &error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((*num_individuals, *ploidy), (0, 2));
        let message = error.to_string();
        assert!(message.contains("0 individuals"), "{message}");
        assert!(message.contains("one individual at least"), "{message}");

        let error = match VarsWriter::new(Vec::new(), &individuals, 0, 3) {
            Ok(_) => panic!("a writer of the ploidy 0 was built"),
            Err(error) => error,
        };
        let Error::VarsFileOfNoGenotypes {
            num_individuals,
            ploidy,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((num_individuals, ploidy), (3, 0));

        let error = match VarsWriter::new(Vec::new(), &individuals, 2, 0) {
            Ok(_) => panic!("a writer of batches of no variant was built"),
            Err(error) => error,
        };
        // The case every reader that takes a size gives, which
        // `write_vars` gives for the same number.
        assert!(matches!(error, Error::BlockOfNoVariants), "{error}");
    }

    /// Arrow keeps where each text of a column of a batch ends in a 32 bit
    /// number, and arrow-rs panics at the text that goes past it, which
    /// would leave a half written file behind. A block of more text than
    /// one column holds is refused instead, with the way out in the
    /// message, and nothing of it is written.
    ///
    /// The limit is an argument of the three columns of texts, so that this
    /// test does not have to hold 2 GiB of text.
    #[test]
    fn a_column_of_more_text_than_an_arrow_column_holds_is_refused() {
        assert_eq!(MAX_COLUMN_BYTES, u64::from(i32::MAX.unsigned_abs()));

        let mut chroms = ChromTable::new();
        let block = cases_block(&mut chroms);
        let chrom = block.chrom.as_deref().expect("the chromosomes");
        let pos = block.pos.as_deref().expect("the positions");
        let ids = block.id.as_deref().expect("the ids");
        let alleles = block.alleles.as_ref().expect("the alleles");

        // The four variants of `cases.vcf` are 16 bytes of chromosome
        // names, `chr1` four times; 3 bytes of ids, `rs1` and the three
        // that are empty; and 8 bytes of alleles, `A`, `T`, `A`, `T`, `A`,
        // `G`, `T` and `T`.
        assert!(chrom_column(chrom, pos, &chroms, 16).is_ok());
        assert!(id_column(ids, 3).is_ok());
        assert!(alleles_column(alleles, 8).is_ok());

        let error = match chrom_column(chrom, pos, &chroms, 15) {
            Ok(_) => panic!("the chromosome names were taken"),
            Err(error) => error,
        };
        let Error::VarsTextTooLarge {
            column,
            found,
            largest,
            ..
        } = &error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((*column, *found, *largest), ("chrom", 16, 15));
        // What a user does about it is write the file in smaller batches.
        let message = error.to_string();
        assert!(message.contains("`chrom`"), "{message}");
        assert!(message.contains("num_vars_per_block"), "{message}");

        let error = match id_column(ids, 2) {
            Ok(_) => panic!("the ids were taken"),
            Err(error) => error,
        };
        let Error::VarsTextTooLarge {
            column,
            found,
            largest,
            ..
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((column, found, largest), ("id", 3, 2));

        let error = match alleles_column(alleles, 7) {
            Ok(_) => panic!("the alleles were taken"),
            Err(error) => error,
        };
        let Error::VarsTextTooLarge {
            column,
            counted,
            found,
            largest,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(
            (column, counted, found, largest),
            ("alleles", "bytes of text", 8, 7)
        );
    }

    /// Arrow keeps where the alleles of each variant end in a 32 bit number
    /// too, and the builder of arrow-rs panics at the entry that goes past
    /// it, so the alleles of a block are counted as well as their bytes: a
    /// block whose alleles are empty texts holds no byte and one entry for
    /// each of them.
    #[test]
    fn a_block_of_more_alleles_than_a_list_column_holds_is_refused() {
        let mut empty = AllelesColumn::with_num_vars(2).expect("the alleles");
        empty.push(&["".to_owned(), "".to_owned(), "".to_owned()]);
        empty.push(&["".to_owned(), "".to_owned()]);
        // Five alleles of no byte at all.
        assert!(alleles_column(&empty, 5).is_ok());

        let error = match alleles_column(&empty, 4) {
            Ok(_) => panic!("the five alleles were taken"),
            Err(error) => error,
        };

        let Error::VarsTextTooLarge {
            column,
            counted,
            found,
            largest,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(
            (column, counted, found, largest),
            ("alleles", "alleles", 5, 4)
        );
    }

    /// The number the system gives when a disc fills up, `ENOSPC`, which
    /// is 28 on macOS and on Linux.
    const NO_SPACE_LEFT: i32 = 28;

    /// A sink that keeps what it was written into a buffer the test holds
    /// and fails once it has taken `takes` bytes, which is what a disc that
    /// fills up under the writer does. It fails with the number of the
    /// system, as a file would.
    #[derive(Clone)]
    struct SinkThatFills {
        written: Rc<RefCell<Vec<u8>>>,
        takes: usize,
    }

    impl SinkThatFills {
        fn of(takes: usize) -> SinkThatFills {
            SinkThatFills {
                written: Rc::new(RefCell::new(Vec::new())),
                takes,
            }
        }
    }

    impl Write for SinkThatFills {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            let mut written = self.written.borrow_mut();
            let left = self.takes.saturating_sub(written.len());
            if left == 0 {
                return Err(std::io::Error::from_raw_os_error(NO_SPACE_LEFT));
            }
            let taken = bytes.len().min(left);
            written.extend_from_slice(&bytes[..taken]);
            Ok(taken)
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    /// A sink that fails is a vars file that could not be written and not a
    /// source that could not be read: which of the two files of the call
    /// went wrong is what a user reads from the exception, and the number
    /// the system gave is what the exception is built with.
    #[test]
    fn a_sink_that_fails_is_the_error_of_a_vars_file_that_could_not_be_written() {
        let mut chroms = ChromTable::new();
        let block = cases_block(&mut chroms);
        let reader = GivenBlocks::of(vec![block], chroms);
        let sink = SinkThatFills::of(200);

        let error = match write_vars(reader, sink.clone(), Some(3)) {
            Ok(_) => panic!("the file was written on a sink that fails"),
            Err(error) => error,
        };

        let Error::VarsFileNotWritten { problem, source } = &error else {
            panic!("the error is {error}");
        };
        assert_eq!(
            source.as_ref().and_then(std::io::Error::raw_os_error),
            Some(NO_SPACE_LEFT)
        );
        assert!(problem.contains("os error 28"), "{problem}");
        let message = error.to_string();
        assert!(message.contains("could not be written"), "{message}");

        // What the doc comment of `write_vars` says the caller is left
        // with: the bytes that were written are on the sink, and they are
        // not a file that can be read.
        let written = sink.written.borrow().clone();
        assert_eq!(written.len(), 200);
        assert_eq!(written.get(..6), Some(b"ARROW1".as_slice()));
        assert!(FileReader::try_new(Cursor::new(written), None).is_err());
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

    /// A reader of popnei refuses a file whose version does not start with
    /// the part it reads, so the version popnei writes starts with it.
    /// Nothing but this ties the two constants together.
    #[test]
    fn the_version_that_is_written_starts_with_the_part_that_is_read() {
        assert_eq!(FORMAT_VERSION.split('.').next(), Some(FORMAT_VERSION_READ));
        assert_eq!(FORMAT_VERSION, "1.0");
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

    /// The pieces of a vars file built with arrow-rs and not with popnei's
    /// writer: one field and one array for each column, the text of the
    /// `popnei` key of the schema, the text of the `popnei_batches` key of
    /// the footer, and how many batches the file holds, the columns written
    /// once for each. A test changes one of them and [`FileParts::written`]
    /// gives the bytes of the file that comes out, which is how the files
    /// another arrow program writes, and popnei's writer does not, are made.
    ///
    /// The file is written with no compression, which popnei reads as it
    /// reads the lz4 that it writes.
    struct FileParts {
        columns: Vec<(Field, ArrayRef)>,
        popnei: Option<String>,
        popnei_batches: Option<String>,
        num_batches: usize,
    }

    impl FileParts {
        /// The six columns of the four variants of `cases.vcf` in one
        /// batch, with the two keys popnei's writer writes for them.
        fn of_cases() -> FileParts {
            let mut chroms = ChromTable::new();
            let block = cases_block(&mut chroms);
            let fields = block.fields();
            let metadata = metadata_of_cases();
            // 3 individuals of the ploidy 2.
            let schema = schema_of(fields, 6, &metadata);
            let writer: VarsWriter<Vec<u8>> =
                VarsWriter::new(Vec::new(), &cases_individuals(), 2, 3).expect("the writer");
            let (arrays, regions) = writer
                .arrays_of(block, &chroms, fields)
                .expect("the columns of the batch");
            let columns = schema
                .fields()
                .iter()
                .map(|field| field.as_ref().clone())
                .zip(arrays)
                .collect();
            FileParts {
                columns,
                popnei: Some(metadata_as_json(&metadata)),
                popnei_batches: Some(batches_as_json(&[BatchInfo {
                    num_vars: 4,
                    regions,
                }])),
                num_batches: 1,
            }
        }

        /// The field and the array of the column of that name, for a test
        /// that puts another one there.
        fn column(&mut self, name: &str) -> &mut (Field, ArrayRef) {
            self.columns
                .iter_mut()
                .find(|(field, _)| field.name() == name)
                .unwrap_or_else(|| panic!("the `{name}` column"))
        }

        /// The bytes of the file these pieces make.
        fn written(&self) -> Vec<u8> {
            let fields: Vec<Field> = self
                .columns
                .iter()
                .map(|(field, _)| field.clone())
                .collect();
            let arrays: Vec<ArrayRef> = self
                .columns
                .iter()
                .map(|(_, array)| Arc::clone(array))
                .collect();
            let mut schema = Schema::new(fields);
            if let Some(popnei) = &self.popnei {
                schema =
                    schema.with_metadata(HashMap::from([(POPNEI_KEY.to_owned(), popnei.clone())]));
            }
            let batch = RecordBatch::try_new(Arc::new(schema.clone()), arrays)
                .expect("the columns are a batch");
            let mut writer =
                FileWriter::try_new(Vec::new(), &schema).expect("the file was started");
            for _ in 0..self.num_batches {
                writer.write(&batch).expect("the batch was written");
            }
            if let Some(batches) = &self.popnei_batches {
                writer.write_metadata(POPNEI_BATCHES_KEY, batches.clone());
            }
            writer.into_inner().expect("the file was finished")
        }
    }

    /// The reader over those bytes, which is how a test opens a vars file
    /// it holds in memory.
    fn opened(bytes: Vec<u8>) -> Result<VarsReader<Cursor<Vec<u8>>>> {
        VarsReader::new(Cursor::new(bytes))
    }

    /// The path of one of the reference files, which live at the root of the
    /// repository, beside the Python tests that read the same ones, and not
    /// inside this crate.
    fn reference(kind: &str, name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference")
            .join(kind)
            .join(name)
    }

    /// The reader over `tests/reference/vars/zstd.vars`, the four variants
    /// of `cases.vcf` in one batch compressed with zstd, which
    /// `tests/reference/vars/make_reference.py` writes with pyarrow because
    /// no build of popnei can write one.
    fn zstd_vars() -> VarsReader<Cursor<Vec<u8>>> {
        let path = reference("vars", "zstd.vars");
        let bytes = match std::fs::read(&path) {
            Ok(bytes) => bytes,
            Err(problem) => panic!("{path:?}: {problem}"),
        };
        opened(bytes).expect("zstd.vars is a vars file")
    }

    /// The reader over `many.vcf`, the 500 variants of 50 diploid
    /// individuals of `docs/specs/io_vcf.md`, with the variants that failed
    /// a filter among them.
    fn many_vcf_reader(num_vars_per_block: Option<usize>) -> VcfReader<BufReader<File>> {
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block,
        };
        let path = reference("vcf", "many.vcf");
        match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(problem) => panic!("{path:?}: {problem}"),
        }
    }

    /// The bytes of the vars file of `many.vcf`, written from the VCF
    /// reader with batches of that many variants.
    fn many_vcf_written(num_vars_per_block: Option<usize>) -> Vec<u8> {
        write_vars(many_vcf_reader(None), Vec::new(), num_vars_per_block)
            .expect("many.vcf was written as a vars file")
            .0
    }

    /// Where the variants of those rows are, one entry for each of their
    /// chromosomes, for the footer of a file built in a test.
    fn regions_of(rows: &[Row]) -> Vec<Region> {
        let mut regions: Vec<Region> = Vec::new();
        for row in rows {
            match regions.iter_mut().find(|region| region.chrom == row.chrom) {
                Some(region) => {
                    region.min_pos = region.min_pos.min(row.pos);
                    region.max_pos = region.max_pos.max(row.pos);
                }
                None => regions.push(Region {
                    chrom: row.chrom.to_owned(),
                    min_pos: row.pos,
                    max_pos: row.pos,
                }),
            }
        }
        regions
    }

    /// The bytes of a file of two batches of the four variants of
    /// `cases.vcf`, whose `pos` column is declared as one that holds nulls,
    /// which is what another program writes, with `positions` as the
    /// positions of its second batch and the entries of its footer saying
    /// `says` variants.
    fn cases_in_two_batches(positions: UInt64Array, says: [usize; 2]) -> Vec<u8> {
        let mut parts = FileParts::of_cases();
        *parts.column(POS_COLUMN) = (
            Field::new(POS_COLUMN, DataType::UInt64, true),
            Arc::clone(&parts.column(POS_COLUMN).1),
        );
        let fields: Vec<Field> = parts
            .columns
            .iter()
            .map(|(field, _)| field.clone())
            .collect();
        let arrays: Vec<ArrayRef> = parts
            .columns
            .iter()
            .map(|(_, array)| Arc::clone(array))
            .collect();
        let later: Vec<ArrayRef> = parts
            .columns
            .iter()
            .map(|(field, array)| match field.name() == POS_COLUMN {
                true => Arc::new(positions.clone()) as ArrayRef,
                false => Arc::clone(array),
            })
            .collect();
        let popnei = parts.popnei.clone().expect("the `popnei` key");
        let schema = Arc::new(
            Schema::new(fields).with_metadata(HashMap::from([(POPNEI_KEY.to_owned(), popnei)])),
        );
        let mut writer = FileWriter::try_new(Vec::new(), &schema).expect("the file was started");
        for arrays in [arrays, later] {
            let batch =
                RecordBatch::try_new(Arc::clone(&schema), arrays).expect("the columns are a batch");
            writer.write(&batch).expect("the batch was written");
        }
        let batches: Vec<BatchInfo> = says
            .iter()
            .map(|num_vars| BatchInfo {
                num_vars: *num_vars,
                regions: regions_of(&CASES),
            })
            .collect();
        writer.write_metadata(POPNEI_BATCHES_KEY, batches_as_json(&batches));
        writer.into_inner().expect("the file was finished")
    }

    /// The bytes of the file with the entry of its first batch in the
    /// footer changed to say that the batch is `metadata_len` and
    /// `body_len` bytes.
    ///
    /// The footer of an arrow file holds one entry of 24 bytes for each
    /// batch: where it starts, 8 bytes; how long its message is, 4; four
    /// bytes of padding; and how long its body is, 8. The entry of the file
    /// is found by the three numbers it holds, which the reader gives.
    fn batch_entry_changed(bytes: Vec<u8>, metadata_len: i32, body_len: i64) -> Vec<u8> {
        let at = opened(bytes.clone())
            .expect("the file is a vars file")
            .blocks[0];
        let mut entry = Vec::new();
        entry.extend_from_slice(&i64::try_from(at.offset).expect("the offset").to_le_bytes());
        entry.extend_from_slice(
            &i32::try_from(at.metadata_len)
                .expect("the message")
                .to_le_bytes(),
        );
        entry.extend_from_slice(&[0; 4]);
        entry.extend_from_slice(&i64::try_from(at.body_len).expect("the body").to_le_bytes());
        let found: Vec<usize> = bytes
            .windows(entry.len())
            .enumerate()
            .filter(|(_, window)| *window == entry.as_slice())
            .map(|(place, _)| place)
            .collect();
        assert_eq!(found.len(), 1, "the entry of the batch in the file");
        let mut changed = bytes;
        let message_at = found[0].saturating_add(8);
        let body_at = found[0].saturating_add(16);
        changed[message_at..message_at.saturating_add(4)]
            .copy_from_slice(&metadata_len.to_le_bytes());
        changed[body_at..body_at.saturating_add(8)].copy_from_slice(&body_len.to_le_bytes());
        changed
    }

    /// The bytes of a file whose batches hold that many of the four
    /// variants of `cases.vcf`, in their order, with the entry of the
    /// footer of each batch to match. A size of 0 gives a batch of no
    /// variants, which no writer of popnei makes and another arrow program
    /// can.
    fn cases_in_batches_of(sizes: &[usize]) -> Vec<u8> {
        let parts = FileParts::of_cases();
        let fields: Vec<Field> = parts
            .columns
            .iter()
            .map(|(field, _)| field.clone())
            .collect();
        let arrays: Vec<ArrayRef> = parts
            .columns
            .iter()
            .map(|(_, array)| Arc::clone(array))
            .collect();
        let popnei = parts.popnei.clone().expect("the `popnei` key");
        let schema = Arc::new(
            Schema::new(fields).with_metadata(HashMap::from([(POPNEI_KEY.to_owned(), popnei)])),
        );
        let whole = RecordBatch::try_new(Arc::clone(&schema), arrays).expect("the four variants");
        let mut writer = FileWriter::try_new(Vec::new(), &schema).expect("the file was started");
        let mut batches = Vec::new();
        let mut first = 0;
        for size in sizes {
            writer
                .write(&whole.slice(first, *size))
                .expect("the batch was written");
            let rows = CASES.get(first..first.saturating_add(*size)).unwrap_or(&[]);
            batches.push(BatchInfo {
                num_vars: *size,
                regions: regions_of(rows),
            });
            first = first.saturating_add(*size);
        }
        writer.write_metadata(POPNEI_BATCHES_KEY, batches_as_json(&batches));
        writer.into_inner().expect("the file was finished")
    }

    /// What the reader said about bytes it refused. That it took them is a
    /// failure of the test itself.
    fn refused(bytes: Vec<u8>) -> Error {
        match opened(bytes) {
            Ok(_) => panic!("the bytes were opened as a vars file"),
            Err(error) => error,
        }
    }

    /// Whether the reader takes those bytes, with nothing of the reader
    /// kept, for a test that looks at what it says of the ones it refuses.
    fn opening(bytes: Vec<u8>) -> Result<()> {
        opened(bytes).map(|_| ())
    }

    /// What a vars file says about itself is in its schema and its footer,
    /// so it is known as soon as the file is opened: the names of the
    /// individuals, the ploidy, the variants of the file and where those of
    /// each batch are.
    #[test]
    fn a_vars_file_says_what_it_holds_and_where_its_variants_are_before_any_batch_is_read() {
        let reader = opened(cases_written_in_batches_of(3)).expect("the file is a vars file");

        assert_eq!(
            reader.metadata(),
            &VarsMetadata {
                format_version: FORMAT_VERSION.to_owned(),
                individuals: cases_individuals(),
                ploidy: 2,
                num_vars_per_block: 3,
            }
        );
        assert_eq!(
            reader.batches(),
            [
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
        // The four variants of the two batches, counted with no batch read.
        assert_eq!(reader.num_vars(), 4);
    }

    /// The footer of an arrow file says where each of its batches is, which
    /// is what the reader seeks to when it reads one, and the file holds
    /// the bytes of every one of them.
    #[test]
    fn the_reader_knows_where_each_batch_of_the_file_is() {
        let bytes = cases_written_in_batches_of(3);
        let file_len = u64::try_from(bytes.len()).expect("the length of the file");
        let reader = opened(bytes).expect("the file is a vars file");

        assert_eq!(reader.schema.fields().len(), 6);
        assert_eq!(reader.blocks.len(), 2);
        let first = reader.blocks[0];
        let second = reader.blocks[1];
        // The batches come after the six bytes of the header, one after
        // another, and the last one ends before the footer.
        assert!(first.offset >= 6, "the first batch is at {first:?}");
        assert!(
            second.offset >= first.offset + first.metadata_len + first.body_len,
            "the batches are at {first:?} and {second:?}"
        );
        assert!(
            second.offset + second.metadata_len + second.body_len < file_len,
            "the last batch is at {second:?} of a file of {file_len} bytes"
        );
    }

    /// A source that does not start as an arrow file is not a vars file and
    /// is a wrong input of the call, which the bytes of a VCF are.
    #[test]
    fn bytes_that_are_not_an_arrow_file_are_not_a_vars_file() {
        let vcf = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\n".to_vec();
        let problem = problem_of(opening(vcf));
        assert!(problem.contains("ARROW1"), "{problem}");

        let problem = problem_of(opening(Vec::new()));
        assert!(problem.contains("ARROW1"), "{problem}");
    }

    /// A file that starts as an arrow file and ends before what it says it
    /// holds is refused and not read as a file of fewer variants: the
    /// variants after the cut are not in it, and a download that stopped is
    /// what a user has to know about.
    #[test]
    fn a_vars_file_that_was_cut_short_is_refused_and_not_read_as_a_file_of_fewer_variants() {
        let whole = cases_written_in_batches_of(3);
        let len = whole.len();
        assert!(len > 400, "the file of the test is {len} bytes");
        // The four bytes before the mark of the end of an arrow file say
        // how many bytes its footer holds.
        let said = <[u8; 4]>::try_from(&whole[len - 10..len - 6]).expect("the four bytes");
        let footer_len = usize::try_from(i32::from_le_bytes(said)).expect("the footer");

        for (what, cut) in [
            ("the middle of a batch", 200),
            ("the middle of its footer", len - 10 - footer_len / 2),
            ("its last byte", len - 1),
        ] {
            let error = refused(whole[..cut].to_vec());
            let Error::VarsFileCutShort { problem } = &error else {
                panic!("the file cut at {what} gave {error}");
            };
            let message = error.to_string();
            assert!(message.contains("cut short"), "{what}: {message}");
            assert!(!problem.is_empty(), "{what}: the error says nothing");
        }

        // The six bytes of the header alone, with nothing of the footer.
        let error = refused(b"ARROW1".to_vec());
        assert!(
            matches!(error, Error::VarsFileCutShort { .. }),
            "the header alone gave {error}"
        );
    }

    /// An arrow file of another program has no `popnei` key, and so has a
    /// vars file of pyNei: what is refused is the source and not the key,
    /// since a user with such a file has the wrong file.
    #[test]
    fn an_arrow_file_whose_schema_has_no_popnei_key_is_not_a_vars_file() {
        let mut parts = FileParts::of_cases();
        parts.popnei = None;

        let problem = problem_of(opening(parts.written()));

        assert!(problem.contains(&format!("`{POPNEI_KEY}`")), "{problem}");
        assert!(problem.contains("schema"), "{problem}");
    }

    /// A reader refuses a file whose major version is not the one it reads,
    /// with the version the file gives, so that a user of an old popnei and
    /// a file of a later format is told to get a later popnei.
    #[test]
    fn a_file_of_another_major_version_is_refused_with_the_version_it_gives() {
        let mut parts = FileParts::of_cases();
        parts.popnei = Some(metadata_as_json(&VarsMetadata {
            format_version: "2.0".to_owned(),
            ..metadata_of_cases()
        }));

        let error = refused(parts.written());

        let Error::VarsFormatVersion { found, read } = &error else {
            panic!("the file of the version 2.0 gave {error}");
        };
        assert_eq!((found.as_str(), *read), ("2.0", FORMAT_VERSION_READ));
        let message = error.to_string();
        assert!(message.contains("2.0"), "{message}");
    }

    /// A file of a later minor version is read, the keys and the columns it
    /// holds that popnei does not know ignored: that is what lets a later
    /// version of the format add a column without making the files or the
    /// readers that are there useless.
    #[test]
    fn a_file_of_a_later_minor_version_is_read() {
        let mut parts = FileParts::of_cases();
        parts.popnei = Some(metadata_as_json(&VarsMetadata {
            format_version: "1.7".to_owned(),
            ..metadata_of_cases()
        }));

        let reader = opened(parts.written()).expect("a file of the version 1.7 is read");

        assert_eq!(reader.metadata().format_version, "1.7");
        assert_eq!(reader.num_vars(), 4);
    }

    /// A column popnei knows is read as the type that column holds and as
    /// no other: a file whose numbers are of another type would be read as
    /// other numbers or not at all, so it is refused with the column, the
    /// type it has and the type popnei reads.
    #[test]
    fn a_column_of_another_type_names_the_column_and_both_types() {
        let mut parts = FileParts::of_cases();
        let positions: ArrayRef = Arc::new(Int32Array::from(vec![100, 200, 300, 400]));
        *parts.column(POS_COLUMN) = (Field::new(POS_COLUMN, DataType::Int32, false), positions);
        let error = refused(parts.written());
        let Error::VarsColumnType {
            column,
            found,
            expected,
        } = &error
        else {
            panic!("the file whose `pos` is Int32 gave {error}");
        };
        assert_eq!(
            (*column, found.as_str(), expected.as_str()),
            ("pos", "Int32", "UInt64")
        );

        let mut parts = FileParts::of_cases();
        let qualities: ArrayRef = Arc::new(Float64Array::from(vec![
            Some(29.5),
            None,
            Some(67.0),
            Some(47.0),
        ]));
        *parts.column(QUAL_COLUMN) = (Field::new(QUAL_COLUMN, DataType::Float64, true), qualities);
        let error = refused(parts.written());
        let Error::VarsColumnType {
            found, expected, ..
        } = &error
        else {
            panic!("the file whose `qual` is Float64 gave {error}");
        };
        assert_eq!((found.as_str(), expected.as_str()), ("Float64", "Float32"));

        // What is inside a list is compared too: a file of one number for
        // each allele is not a file of the texts of the alleles.
        let mut parts = FileParts::of_cases();
        let mut alleles = ListBuilder::new(Int8Builder::new());
        for row in 0..4 {
            alleles.values().append_value(row);
            alleles.append(true);
        }
        let inside = Arc::new(Field::new(ITEM_FIELD, DataType::Int8, true));
        *parts.column(ALLELES_COLUMN) = (
            Field::new(ALLELES_COLUMN, DataType::List(inside), false),
            Arc::new(alleles.finish()),
        );
        let error = refused(parts.written());
        let Error::VarsColumnType {
            found, expected, ..
        } = &error
        else {
            panic!("the file whose alleles are numbers gave {error}");
        };
        assert_eq!(
            (found.as_str(), expected.as_str()),
            ("List<Int8>", "List<Utf8>")
        );
    }

    /// The width of the `gts` column is what turns its flat buffer into
    /// variants, so a file whose width and whose `popnei` key say different
    /// things is refused and not read one allele beside another.
    #[test]
    fn a_gts_column_of_another_width_than_the_individuals_and_the_ploidy_is_refused() {
        let mut parts = FileParts::of_cases();
        let inside = Arc::new(Field::new(ITEM_FIELD, DataType::Int8, true));
        // 4 variants of 7 alleles each, where the three individuals of the
        // ploidy 2 of the `popnei` key hold 6.
        let alleles = Int8Array::new(ScalarBuffer::from(vec![0_i8; 28]), None);
        let column = FixedSizeListArray::try_new(Arc::clone(&inside), 7, Arc::new(alleles), None)
            .expect("the genotypes of seven alleles");
        *parts.column(GTS_COLUMN) = (
            Field::new(GTS_COLUMN, DataType::FixedSizeList(inside, 7), false),
            Arc::new(column),
        );

        let error = refused(parts.written());

        let Error::VarsGtsWidth {
            found,
            expected,
            num_individuals,
            ploidy,
        } = error
        else {
            panic!("the file of seven alleles for each variant gave {error}");
        };
        assert_eq!((found, expected, num_individuals, ploidy), (7, 6, 3, 2));
    }

    /// Arrow gives the values inside a list a field of their own, and arrow
    /// programs differ in the name of that field and in whether it takes
    /// nulls. What a reader compares is what the list holds, so a file
    /// whose alleles and whose genotypes are under another name is read.
    #[test]
    fn the_name_of_the_field_inside_a_list_is_not_what_the_reader_compares() {
        let mut parts = FileParts::of_cases();
        let (_, alleles) = parts.column(ALLELES_COLUMN).clone();
        let alleles = alleles
            .as_any()
            .downcast_ref::<ListArray>()
            .expect("the alleles are a list")
            .clone();
        let (_, offsets, texts, nulls) = alleles.into_parts();
        // The name of another program, and the nulls pyarrow allows the
        // values of a list, where popnei writes `item` and no null.
        let inside = Arc::new(Field::new("element", DataType::Utf8, true));
        *parts.column(ALLELES_COLUMN) = (
            Field::new(ALLELES_COLUMN, DataType::List(Arc::clone(&inside)), false),
            Arc::new(ListArray::new(inside, offsets, texts, nulls)),
        );

        let (_, gts) = parts.column(GTS_COLUMN).clone();
        let gts = gts
            .as_any()
            .downcast_ref::<FixedSizeListArray>()
            .expect("the genotypes are a fixed size list")
            .clone();
        let (_, width, alleles, nulls) = gts.into_parts();
        let inside = Arc::new(Field::new("element", DataType::Int8, true));
        let column = FixedSizeListArray::try_new(Arc::clone(&inside), width, alleles, nulls)
            .expect("the genotypes under another name");
        *parts.column(GTS_COLUMN) = (
            Field::new(GTS_COLUMN, DataType::FixedSizeList(inside, width), false),
            Arc::new(column),
        );

        let bytes = parts.written();
        let reader = opened(bytes.clone()).expect("the file is read");

        assert_eq!(reader.num_vars(), 4);
        assert_eq!(reader.columns.alleles, Some(3));
        assert_eq!(reader.columns.gts, Some(5));

        // And its variants are read: the blocks of such a file are those of
        // a file of popnei.
        let (blocks, chroms) = blocks_read(bytes, Needs::ALL);
        let expected: Vec<ReadRow> = CASES.iter().map(row_read).collect();
        assert_eq!(rows_of(&blocks, &chroms), expected);
    }

    /// A column popnei does not know is ignored, which is what lets a later
    /// version of the format add one, and the columns it knows are read at
    /// the place they have in the file and not at the place the table of
    /// the spec gives them.
    #[test]
    fn a_column_popnei_does_not_know_is_ignored_and_the_others_are_read_where_they_are() {
        let mut parts = FileParts::of_cases();
        let depth: ArrayRef = Arc::new(Int32Array::from(vec![Some(10), None, Some(30), Some(40)]));
        parts
            .columns
            .insert(0, (Field::new("depth", DataType::Int32, true), depth));

        let reader = opened(parts.written()).expect("the file with a seventh column is read");

        assert_eq!(reader.schema.fields().len(), 7);
        assert_eq!(
            reader.columns,
            VarsColumns {
                chrom: Some(1),
                pos: Some(2),
                id: Some(3),
                alleles: Some(4),
                qual: Some(5),
                gts: Some(6),
            }
        );
        assert_eq!(reader.num_vars(), 4);
    }

    /// A file holds the columns its source could fill, so a file whose
    /// source gave the genotypes alone has one column and its batches say
    /// nothing about where their variants are.
    #[test]
    fn a_file_of_the_genotypes_alone_is_opened_and_its_batches_have_no_region() {
        let mut chroms = ChromTable::new();
        let mut block = cases_block(&mut chroms);
        block.chrom = None;
        block.pos = None;
        block.id = None;
        block.alleles = None;
        block.qual = None;
        let reader = GivenBlocks::of(vec![block], chroms);
        let (bytes, _) = write_vars(reader, Vec::new(), Some(3)).expect("the file was written");

        let reader = opened(bytes).expect("the file of one column is a vars file");

        assert_eq!(
            reader.columns,
            VarsColumns {
                chrom: None,
                pos: None,
                id: None,
                alleles: None,
                qual: None,
                gts: Some(0),
            }
        );
        assert_eq!(reader.num_vars(), 4);
        assert_eq!(
            reader.batches(),
            [
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

    /// A source with no variants is written and is read back as no
    /// variants, as a VCF with no variants is read and is not an error:
    /// the file has the two keys, one column and no batch.
    #[test]
    fn a_file_with_no_variants_is_opened_and_says_it_holds_none() {
        let reader = GivenBlocks::of(Vec::new(), ChromTable::new());
        let (bytes, _) = write_vars(reader, Vec::new(), Some(3)).expect("the file was written");

        let reader = opened(bytes).expect("the file with no batch is a vars file");

        assert_eq!(reader.num_vars(), 0);
        assert_eq!(reader.batches(), []);
        assert_eq!(reader.blocks, Vec::new());
        assert_eq!(reader.metadata().individuals, cases_individuals());
    }

    /// The key of the footer is what says how many variants the file holds
    /// and where they are, and a file without it is not a vars file: an
    /// arrow file of another program has the columns and none of that.
    #[test]
    fn a_file_whose_footer_has_no_popnei_batches_key_is_not_a_vars_file() {
        let mut parts = FileParts::of_cases();
        parts.popnei_batches = None;

        let problem = problem_of(opening(parts.written()));

        assert!(
            problem.contains(&format!("`{POPNEI_BATCHES_KEY}`")),
            "{problem}"
        );
        assert!(problem.contains("footer"), "{problem}");
    }

    /// There is one entry of the footer for each batch, and a file with
    /// another number of one than of the other is refused with both counts:
    /// no entry of such a file can be trusted to be that of its batch.
    #[test]
    fn a_footer_whose_entries_are_not_as_many_as_the_batches_is_refused_with_both_counts() {
        let mut parts = FileParts::of_cases();
        // The same columns written twice, with the one entry of a file of
        // one batch left in the footer.
        parts.num_batches = 2;

        let error = refused(parts.written());

        let Error::VarsBatchesDoNotMatch { found, expected } = error else {
            panic!("the file of two batches and one entry gave {error}");
        };
        assert_eq!((found, expected), (1, 2));
    }

    /// The names of the individuals are how a user asks for one, so two of
    /// one name are refused, with the name, as they are for the VCF reader:
    /// they would be one individual for the user and two columns of
    /// genotypes in the file.
    #[test]
    fn two_individuals_of_one_name_are_refused_with_the_name() {
        let mut parts = FileParts::of_cases();
        parts.popnei = Some(metadata_as_json(&VarsMetadata {
            individuals: vec!["ind1".to_owned(), "ind2".to_owned(), "ind1".to_owned()],
            ..metadata_of_cases()
        }));

        let error = refused(parts.written());

        let Error::VarsIndividualTwice { name } = &error else {
            panic!("the file of two individuals of one name gave {error}");
        };
        assert_eq!(name, "ind1");
    }

    /// `from_path` opens the file at the path for a caller that has a path,
    /// and what cannot be opened carries the path, a directory among them:
    /// opening a directory succeeds on macOS and on Linux and only the
    /// first read of it fails, so a reader that did not ask would say that
    /// the directory is not a vars file.
    #[test]
    fn from_path_opens_the_file_at_the_path_and_refuses_a_directory_and_a_path_that_is_not_there() {
        let path = std::env::temp_dir().join(format!(
            "popnei-{process}-from-path.vars",
            process = std::process::id()
        ));
        std::fs::write(&path, cases_written_in_batches_of(3)).expect("the file was written");
        let reader = VarsReader::from_path(&path).expect("the file at the path is a vars file");
        assert_eq!(reader.num_vars(), 4);
        assert_eq!(reader.metadata().individuals, cases_individuals());
        drop(reader);
        std::fs::remove_file(&path).expect("the file was taken away");

        let directory = Path::new(env!("CARGO_MANIFEST_DIR"));
        let error = match VarsReader::from_path(directory) {
            Ok(_) => panic!("a directory was opened as a vars file"),
            Err(error) => error,
        };
        let Error::FileNotOpened {
            path: named,
            source,
        } = &error
        else {
            panic!("the directory gave {error}");
        };
        assert_eq!(named, directory);
        // `EISDIR`, which is what makes it the `IsADirectoryError` of
        // Python.
        assert_eq!(source.raw_os_error(), Some(21));

        let missing = directory.join("there-is-no-such-vars-file.vars");
        let error = match VarsReader::from_path(&missing) {
            Ok(_) => panic!("a path with no file at it was opened"),
            Err(error) => error,
        };
        let Error::FileNotOpened {
            path: named,
            source,
        } = &error
        else {
            panic!("the path with no file at it gave {error}");
        };
        assert_eq!(*named, missing);
        assert_eq!(source.kind(), ErrorKind::NotFound);
    }

    /// A `popnei` key that names no individual is refused when the file is
    /// opened, with the error the writer asked for such a file gives: the
    /// genotypes of no individual hold no allele, and the blocks such a file
    /// gives are variants of nobody.
    #[test]
    fn a_popnei_key_that_names_no_individual_is_refused_with_both_numbers() {
        let mut parts = FileParts::of_cases();
        parts.popnei = Some(metadata_as_json(&VarsMetadata {
            individuals: Vec::new(),
            ..metadata_of_cases()
        }));
        // The `gts` column of such a file, four variants of no allele,
        // which is what its width of 0 says.
        let inside = Arc::new(Field::new(ITEM_FIELD, DataType::Int8, true));
        let column = FixedSizeListArray::try_new_with_length(
            Arc::clone(&inside),
            0,
            Arc::new(Int8Array::from(Vec::<i8>::new())),
            None,
            4,
        )
        .expect("the genotypes of no allele");
        *parts.column(GTS_COLUMN) = (
            Field::new(GTS_COLUMN, DataType::FixedSizeList(inside, 0), false),
            Arc::new(column),
        );

        let error = refused(parts.written());

        let Error::VarsFileOfNoGenotypes {
            num_individuals,
            ploidy,
        } = error
        else {
            panic!("the file of no individual gave {error}");
        };
        assert_eq!((num_individuals, ploidy), (0, 2));
    }

    /// Every vars file has a `gts` column, so a file without one is not a
    /// vars file: a reader of this version cannot tell a file whose column
    /// was lost from one that a later version of the format wrote without
    /// it.
    #[test]
    fn a_file_with_no_gts_column_is_not_a_vars_file() {
        let mut parts = FileParts::of_cases();
        parts
            .columns
            .retain(|(field, _)| field.name() != GTS_COLUMN);

        let problem = problem_of(opening(parts.written()));

        assert!(problem.contains(&format!("`{GTS_COLUMN}`")), "{problem}");
    }

    /// One variant as a block gives it, for the comparison with the rows of
    /// the tables of `cases.vcf` and of the variants that are not sorted:
    /// the name of its chromosome and not its number, and `None` for the
    /// quality of a variant that has none, which a block holds as a NaN.
    ///
    /// A field is `None` when the block has no such column, so a block that
    /// lost a column is not one that holds another value there.
    #[derive(Debug, PartialEq)]
    struct ReadRow {
        chrom: Option<String>,
        pos: Option<u64>,
        id: Option<String>,
        alleles: Option<Vec<String>>,
        qual: Option<Option<f32>>,
        gts: Vec<i8>,
    }

    /// The variants of the blocks, in their order, with the names behind the
    /// chromosome numbers taken from `chroms`, the table of the reader that
    /// gave them.
    fn rows_of(blocks: &[Block], chroms: &ChromTable) -> Vec<ReadRow> {
        let mut rows = Vec::new();
        for block in blocks {
            for variant in block.variants() {
                rows.push(ReadRow {
                    chrom: variant.chrom().map(|number| {
                        chroms
                            .name(number)
                            .unwrap_or_else(|| panic!("the name of the chromosome {number}"))
                            .to_owned()
                    }),
                    pos: variant.pos(),
                    id: variant.id().map(str::to_owned),
                    alleles: variant.num_alleles().map(|num_alleles| {
                        // An allele of an empty text, which no source of
                        // popnei gives and a damaged file holds, reads as
                        // none: it is compared as the empty text it is.
                        (0..num_alleles)
                            .map(|allele| variant.allele(allele).unwrap_or("").to_owned())
                            .collect()
                    }),
                    qual: variant
                        .qual()
                        .map(|quality| (!quality.is_nan()).then_some(quality)),
                    gts: variant.gts().to_vec(),
                });
            }
        }
        rows
    }

    /// One row of the tables of the variants as a block would give it, which
    /// is what a round trip through the file compares with.
    fn row_read(row: &Row) -> ReadRow {
        ReadRow {
            chrom: Some(row.chrom.to_owned()),
            pos: Some(row.pos),
            id: Some(row.id.to_owned()),
            alleles: Some(row.alleles.iter().map(|text| (*text).to_owned()).collect()),
            qual: Some(row.qual),
            gts: row.gts.to_vec(),
        }
    }

    /// Every block of a reader, until it has no more or it fails.
    fn blocks_of(reader: &mut impl BlockReader) -> Result<Vec<Block>> {
        let mut blocks = Vec::new();
        while let Some(block) = reader.next_block()? {
            blocks.push(block);
        }
        Ok(blocks)
    }

    /// How many variants each block holds.
    fn num_vars_of(blocks: &[Block]) -> Vec<usize> {
        blocks.iter().map(|block| block.num_vars).collect()
    }

    /// The blocks of the file in those bytes, asked for `needs`, and the
    /// names of the chromosomes of the reader that gave them.
    fn blocks_read(bytes: Vec<u8>, needs: Needs) -> (Vec<Block>, ChromTable) {
        let mut reader = opened(bytes).expect("the bytes are a vars file");
        reader.set_needs(needs);
        let blocks = match blocks_of(&mut reader) {
            Ok(blocks) => blocks,
            Err(error) => panic!("the blocks of the file: {error}"),
        };
        for block in &blocks {
            block.check().expect("a block of the reader");
        }
        (blocks, reader.chroms)
    }

    /// What the reader said about a file whose blocks it refused. That it
    /// gave them is a failure of the test itself.
    fn refused_at_the_block(bytes: Vec<u8>) -> Error {
        let mut reader = opened(bytes).expect("the bytes are a vars file");
        match blocks_of(&mut reader) {
            Ok(blocks) => panic!("the file gave {count} blocks", count = blocks.len()),
            Err(error) => error,
        }
    }

    /// The reader is a source: no filter stands between it and the file,
    /// before a block is read and after the last one.
    #[test]
    fn a_vars_file_reader_gives_no_filtering_stats() {
        let mut reader = opened(cases_written_in_batches_of(3)).expect("the bytes are a vars file");
        assert!(reader.filtering_stats().is_empty());
        let blocks = blocks_of(&mut reader).expect("the blocks of the file");
        assert_eq!(num_vars_of(&blocks), [3, 1]);
        assert!(reader.filtering_stats().is_empty());
    }

    /// The four variants of the table of `cases.vcf`, written in batches of
    /// three and read back: two blocks, of 3 variants and of 1, with every
    /// field of every variant as it went in, the empty id of the last three
    /// among them.
    #[test]
    fn the_four_variants_of_cases_vcf_are_read_back_as_two_blocks_with_every_field() {
        let (blocks, chroms) = blocks_read(cases_written_in_batches_of(3), Needs::ALL);

        assert_eq!(num_vars_of(&blocks), [3, 1]);
        let expected: Vec<ReadRow> = CASES.iter().map(row_read).collect();
        assert_eq!(rows_of(&blocks, &chroms), expected);
        for block in &blocks {
            assert_eq!(block.fields(), Needs::ALL);
            assert_eq!(block.num_individuals, CASES_INDIVIDUALS);
            assert_eq!(block.ploidy, CASES_PLOIDY);
        }
    }

    /// The variants of a file are in no order in general, and the reader
    /// gives them in the order of the file: the numbers of the chromosomes
    /// are the order in which their names first appear among the variants
    /// that were given, as in the VCF reader.
    #[test]
    fn the_chromosome_numbers_of_the_blocks_are_the_order_in_which_the_names_first_appear() {
        let mut chroms = ChromTable::new();
        let rows: Vec<&Row> = NOT_SORTED.iter().collect();
        let block = block_of(&rows, &mut chroms);
        let reader = GivenBlocks::of(vec![block], chroms);
        let (bytes, _) = write_vars(reader, Vec::new(), Some(4)).expect("the file was written");

        let (blocks, chroms) = blocks_read(bytes, Needs::ALL);

        assert_eq!(num_vars_of(&blocks), [4]);
        let expected: Vec<ReadRow> = NOT_SORTED.iter().map(row_read).collect();
        assert_eq!(rows_of(&blocks, &chroms), expected);
        // `chr2` is the chromosome of the first variant of the file, so it
        // is the number 0, although it comes second in the alphabet.
        assert_eq!(chroms.name(0), Some("chr2"));
        assert_eq!(chroms.name(1), Some("chr1"));
        assert_eq!(
            blocks[0].chrom.as_deref(),
            Some([0, 1, 0, 1].as_slice()),
            "the number of a name is looked up once and used for the variants that follow it"
        );
    }

    /// A file holds the columns its source could fill, and a field whose
    /// column the file lacks gives a block without that column: a source
    /// whose blocks carried the genotypes alone gives a file of one column,
    /// and the blocks read back hold the genotypes and nothing else,
    /// although every field was asked for.
    #[test]
    fn a_file_of_the_genotypes_alone_gives_blocks_of_the_genotypes_alone() {
        let mut chroms = ChromTable::new();
        let mut block = cases_block(&mut chroms);
        block.chrom = None;
        block.pos = None;
        block.id = None;
        block.alleles = None;
        block.qual = None;
        let reader = GivenBlocks::of(vec![block], chroms);
        let (bytes, _) = write_vars(reader, Vec::new(), Some(3)).expect("the file was written");

        let (blocks, chroms) = blocks_read(bytes, Needs::ALL);

        assert_eq!(num_vars_of(&blocks), [3, 1]);
        let expected: Vec<ReadRow> = CASES
            .iter()
            .map(|row| ReadRow {
                chrom: None,
                pos: None,
                id: None,
                alleles: None,
                qual: None,
                gts: row.gts.to_vec(),
            })
            .collect();
        assert_eq!(rows_of(&blocks, &chroms), expected);
        for block in &blocks {
            assert_eq!(block.fields(), Needs::GTS);
        }
    }

    /// The reader gives the batches of the file as they are, so a file
    /// written with batches of a hundred variants is read back in blocks of
    /// a hundred with no `reblock` in between.
    #[test]
    fn a_file_written_with_batches_of_a_hundred_gives_blocks_of_a_hundred() {
        let bytes = many_vcf_written(Some(100));

        let (blocks, _) = blocks_read(bytes, Needs::ALL);

        assert_eq!(num_vars_of(&blocks), [100, 100, 100, 100, 100]);
    }

    /// `many.vcf` through the VCF reader, written and read back: the
    /// genotypes, the chromosome names and the positions of the blocks of
    /// the file are those of the blocks of the VCF, variant by variant.
    ///
    /// What the VCF reader gives is checked against bcftools and against
    /// pyNei in `docs/specs/io_vcf.md` and in `docs/specs/block.md`, so this
    /// carries those checks over to the file.
    #[test]
    fn the_blocks_of_many_vcf_are_read_back_from_the_file_they_were_written_to() {
        let mut vcf = many_vcf_reader(None);
        let from_the_vcf = blocks_of(&mut vcf).expect("the blocks of many.vcf");
        let expected = rows_of(&from_the_vcf, vcf.chroms());
        // The 500 variants of `many.vcf`, which is read with the variants
        // that failed a filter among them, in one block of the size popnei
        // chooses for 50 individuals.
        assert_eq!(num_vars_of(&from_the_vcf), [500]);

        let (blocks, chroms) = blocks_read(many_vcf_written(None), Needs::ALL);

        assert_eq!(num_vars_of(&blocks), [500]);
        let read = rows_of(&blocks, &chroms);
        assert_eq!(read.len(), expected.len());
        for (var, (read, expected)) in read.iter().zip(expected.iter()).enumerate() {
            assert_eq!(read, expected, "the variant {var} of many.vcf");
        }
        // The two chromosomes of the file, in the order of the VCF.
        assert_eq!(chroms.name(0), Some("chr1"));
        assert_eq!(chroms.name(1), Some("chr2"));
    }

    /// Genotypes that repeat compress below one bit each, and a file of
    /// them is read back: 200 diploid individuals whose every genotype is
    /// `0/1` give, for 5000 variants, 2000000 alleles in a file that lz4
    /// leaves under 250000 bytes.
    ///
    /// Issue 2 of the repository, of 24 September 2026: popnei wrote such a
    /// file and then refused to read it, saying it was damaged, because it
    /// bounded the values that a column of a batch says it holds by the
    /// bits of that batch as it lies on disk, where it is compressed, and
    /// not by what the column can hold once it is decompressed.
    #[test]
    fn a_file_whose_genotypes_compress_below_a_bit_each_is_read_back() {
        const NUM_VARS: usize = 5000;
        const NUM_INDIVIDUALS: usize = 200;
        // 5000 variants x 200 individuals x 2 alleles.
        const ALLELES: usize = 2_000_000;
        let individuals: Vec<String> = (0..NUM_INDIVIDUALS).map(|at| format!("ind{at}")).collect();
        let mut chroms = ChromTable::new();
        let chrom = chroms.intern("chr1");
        let mut alleles = AllelesColumn::with_num_vars(NUM_VARS).expect("the alleles");
        let two = ["A".to_owned(), "T".to_owned()];
        for _ in 0..NUM_VARS {
            alleles.push(&two);
        }
        let block = Block {
            num_vars: NUM_VARS,
            num_individuals: NUM_INDIVIDUALS,
            ploidy: 2,
            gts: [0i8, 1].into_iter().cycle().take(ALLELES).collect(),
            chrom: Some(vec![chrom; NUM_VARS]),
            pos: Some((1u64..).take(NUM_VARS).collect()),
            id: Some(vec![String::new(); NUM_VARS]),
            alleles: Some(alleles),
            qual: Some(vec![30.0; NUM_VARS]),
        };
        let expected = block.gts.clone();
        let reader = GivenBlocks {
            individuals,
            ploidy: 2,
            chroms,
            left: vec![block],
            asked_for: Arc::new(Mutex::new(Needs::empty())),
        };

        let bytes = write_vars(reader, Vec::new(), Some(NUM_VARS))
            .expect("the blocks were written as a vars file")
            .0;

        // Under a bit for each of the 2000000 alleles, which is what the
        // file has to be for this test to hold the case of the issue.
        assert!(
            bytes.len() < 250_000,
            "the file is {found} bytes, so its genotypes did not compress below a bit each",
            found = bytes.len()
        );
        let (blocks, chroms) = blocks_read(bytes, Needs::ALL);
        assert_eq!(num_vars_of(&blocks), [NUM_VARS]);
        let block = blocks.first().expect("the block of the file");
        assert_eq!(block.gts, expected);
        assert_eq!(chroms.name(0), Some("chr1"));
    }

    /// The reader is asked for the fields the consumer wants and gives a
    /// block with those columns alone, and a change of them holds from the
    /// next block, as `docs/specs/block.md` asks of every reader: the
    /// projection of each batch is chosen when that batch is read.
    #[test]
    fn the_genotypes_alone_asked_for_give_blocks_with_no_other_column_and_a_change_holds_at_once() {
        let mut reader = opened(cases_written_in_batches_of(3)).expect("the file is a vars file");
        reader.set_needs(Needs::GTS);

        let first = reader
            .next_block()
            .expect("the first block")
            .expect("a block");
        assert_eq!(first.fields(), Needs::GTS);
        assert_eq!(first.num_vars, 3);
        assert_eq!(
            first.gts,
            [
                0, 0, 0, 1, 1, 1, MISSING, MISSING, 0, 1, MISSING, 0, 1, 2, 2, 1, 2, 2
            ]
        );
        assert_eq!(first.chrom, None);
        assert_eq!(first.pos, None);
        assert_eq!(first.id, None);
        assert_eq!(first.qual, None);

        reader.set_needs(Needs::CHROM_POS | Needs::GTS);
        let second = reader
            .next_block()
            .expect("the second block")
            .expect("a block");
        assert_eq!(second.fields(), Needs::CHROM_POS | Needs::GTS);
        assert_eq!(second.pos.as_deref(), Some([400].as_slice()));
        assert_eq!(second.id, None);
        assert!(reader.next_block().expect("the end of the file").is_none());
    }

    /// Only the columns a consumer asks for are decompressed: a `Needs`
    /// becomes a list of the places of those columns in the schema, which
    /// arrow-rs walks past the buffers of the rest with, without
    /// decompressing them.
    ///
    /// What shows it is `zstd.vars`, whose buffers no build of popnei can
    /// decompress: a consumer that asks for no field gets its four variants
    /// and no column, and one that asks for the genotypes gets the error of
    /// a file compressed with zstd.
    #[test]
    fn only_the_columns_that_were_asked_for_are_decompressed() {
        let mut reader = zstd_vars();
        reader.set_needs(Needs::empty());

        let block = reader
            .next_block()
            .expect("no buffer of the batch was decompressed")
            .expect("a block");

        assert_eq!(block.num_vars, 4);
        assert_eq!(block.fields(), Needs::empty());
        assert!(block.gts.is_empty());
        assert_eq!(block.chrom, None);
    }

    /// The columns that a `Needs` asks the reader to decompress are those
    /// of the file, in the order the file has them, and a field whose
    /// column the file lacks is not among them.
    #[test]
    fn the_projection_of_a_needs_is_the_columns_of_the_file_it_asks_for() {
        let six = VarsColumns {
            chrom: Some(0),
            pos: Some(1),
            id: Some(2),
            alleles: Some(3),
            qual: Some(4),
            gts: Some(5),
        };
        assert_eq!(
            projection_of(Needs::ALL, &six),
            vec![
                (VarsColumn::Chrom, 0),
                (VarsColumn::Pos, 1),
                (VarsColumn::Id, 2),
                (VarsColumn::Alleles, 3),
                (VarsColumn::Qual, 4),
                (VarsColumn::Gts, 5),
            ]
        );
        assert_eq!(projection_of(Needs::GTS, &six), vec![(VarsColumn::Gts, 5)]);
        assert_eq!(
            projection_of(Needs::CHROM_POS, &six),
            vec![(VarsColumn::Chrom, 0), (VarsColumn::Pos, 1)]
        );
        assert_eq!(projection_of(Needs::empty(), &six), Vec::new());

        // A file of one column: every field is asked for and the only one
        // the file has is given.
        let one = VarsColumns {
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
            gts: Some(0),
        };
        assert_eq!(projection_of(Needs::ALL, &one), vec![(VarsColumn::Gts, 0)]);
        assert_eq!(projection_of(Needs::ID, &one), Vec::new());

        // The chromosome and the position are one field: a file with one of
        // the two columns and not the other gives neither, since a block
        // holds both or neither.
        let half = VarsColumns {
            chrom: Some(0),
            pos: None,
            ..one
        };
        assert_eq!(projection_of(Needs::ALL, &half), vec![(VarsColumn::Gts, 0)]);

        // The columns are taken at the place they have in the file, which a
        // column popnei does not know moves.
        let moved = VarsColumns {
            chrom: Some(3),
            pos: Some(4),
            id: None,
            alleles: None,
            qual: None,
            gts: Some(1),
        };
        assert_eq!(
            projection_of(Needs::ALL, &moved),
            vec![
                (VarsColumn::Gts, 1),
                (VarsColumn::Chrom, 3),
                (VarsColumn::Pos, 4),
            ]
        );
    }

    /// A batch of no variants is not given as a block, since
    /// `docs/specs/block.md` says that a reader never gives one: the reader
    /// takes the next batch, as a filter does with a block it emptied. No
    /// writer of popnei makes such a batch and another arrow program can.
    #[test]
    fn a_batch_of_no_variants_between_two_others_is_passed_over() {
        let (blocks, chroms) = blocks_read(cases_in_batches_of(&[2, 0, 2]), Needs::ALL);

        assert_eq!(num_vars_of(&blocks), [2, 2]);
        let expected: Vec<ReadRow> = CASES.iter().map(row_read).collect();
        assert_eq!(rows_of(&blocks, &chroms), expected);
    }

    /// A file with no variants is read and is not an error, as a VCF with no
    /// variants is: it holds the two keys, one column and no batch, and its
    /// reader gives no block at the first call and at every call after.
    #[test]
    fn a_file_with_no_variants_gives_no_block() {
        let reader = GivenBlocks::of(Vec::new(), ChromTable::new());
        let (bytes, _) = write_vars(reader, Vec::new(), Some(3)).expect("the file was written");
        let mut reader = opened(bytes).expect("the file with no batch is a vars file");

        assert!(reader.next_block().expect("the empty file").is_none());
        assert!(reader.next_block().expect("the empty file again").is_none());
    }

    /// A file written with no compression is read as the lz4 that popnei
    /// writes is: another arrow program writes one, and the compression of a
    /// batch is what the batch itself says.
    #[test]
    fn a_file_written_with_no_compression_is_read() {
        let (blocks, chroms) = blocks_read(FileParts::of_cases().written(), Needs::ALL);

        assert_eq!(num_vars_of(&blocks), [4]);
        let expected: Vec<ReadRow> = CASES.iter().map(row_read).collect();
        assert_eq!(rows_of(&blocks, &chroms), expected);
    }

    /// The two binding crates hold their reader as a `Box<dyn BlockReader>`,
    /// because neither a pyo3 class nor a wasm-bindgen class can be generic,
    /// so the reader of a vars file is one.
    #[test]
    fn a_vars_reader_gives_its_blocks_as_a_boxed_block_reader() {
        let reader = opened(cases_written_in_batches_of(3)).expect("the file is a vars file");
        let mut boxed: Box<dyn BlockReader> = Box::new(reader);

        boxed.set_needs(Needs::GTS);
        let blocks = blocks_of(&mut boxed).expect("the blocks of the file");

        assert_eq!(boxed.individuals(), cases_individuals());
        assert_eq!(boxed.ploidy(), 2);
        assert_eq!(num_vars_of(&blocks), [3, 1]);
    }

    /// A null in a column where every variant has a value is an error that
    /// names the column and the variant, counted from 1 over the whole file
    /// and not inside its batch, so that the message points at the variant a
    /// user counts.
    #[test]
    fn a_null_position_names_the_column_and_the_variant() {
        let mut parts = FileParts::of_cases();
        let positions: ArrayRef = Arc::new(UInt64Array::from(vec![
            Some(100),
            Some(200),
            None,
            Some(400),
        ]));
        *parts.column(POS_COLUMN) = (Field::new(POS_COLUMN, DataType::UInt64, true), positions);

        let error = refused_at_the_block(parts.written());

        let Error::VarsNullValue { column, var } = error else {
            panic!("the file whose third position is a null gave {error}");
        };
        assert_eq!((column, var), ("pos", 3));
        let message = Error::VarsNullValue { column, var }.to_string();
        assert!(message.contains("`pos`"), "{message}");
        assert!(message.contains("variant 3"), "{message}");
    }

    /// A null in the chromosome, in the alleles or in the genotypes is an
    /// error too, and a null inside one of the two lists is a null of that
    /// column: arrow gives the values inside a list a field of their own
    /// that can hold nulls, and no value popnei writes is one.
    #[test]
    fn a_null_chromosome_a_null_allele_and_a_null_genotype_are_errors_with_their_column() {
        let mut parts = FileParts::of_cases();
        let names: ArrayRef = Arc::new(StringArray::from(vec![
            Some("chr1"),
            None,
            Some("chr1"),
            Some("chr1"),
        ]));
        *parts.column(CHROM_COLUMN) = (Field::new(CHROM_COLUMN, DataType::Utf8, true), names);
        let error = refused_at_the_block(parts.written());
        let Error::VarsNullValue { column, var } = error else {
            panic!("the file whose second chromosome is a null gave {error}");
        };
        assert_eq!((column, var), ("chrom", 2));

        // A null among the alleles of the first variant, where the file
        // holds the list and the allele inside it is the null.
        let mut parts = FileParts::of_cases();
        let mut alleles = ListBuilder::new(StringBuilder::new());
        for row in 0..4 {
            match row {
                0 => {
                    alleles.values().append_value("A");
                    alleles.values().append_null();
                }
                _ => {
                    alleles.values().append_value("A");
                    alleles.values().append_value("T");
                }
            }
            alleles.append(true);
        }
        // The field of the values is declared as holding nulls, which is
        // what another program writes and what a null allele needs: popnei
        // writes that field as holding none.
        let inside = Arc::new(Field::new(ITEM_FIELD, DataType::Utf8, true));
        *parts.column(ALLELES_COLUMN) = (
            Field::new(ALLELES_COLUMN, DataType::List(inside), false),
            Arc::new(alleles.finish()),
        );
        let error = refused_at_the_block(parts.written());
        let Error::VarsNullValue { column, var } = error else {
            panic!("the file whose first variant has a null allele gave {error}");
        };
        assert_eq!((column, var), ("alleles", 1));

        // A null among the genotypes of the last variant, which is an
        // allele that was not called written as a null and not as the -1
        // that popnei writes.
        let mut parts = FileParts::of_cases();
        let mut genotypes = Int8Builder::new();
        for allele in 0..24 {
            match allele {
                20 => genotypes.append_null(),
                _ => genotypes.append_value(0),
            }
        }
        let inside = Arc::new(Field::new(ITEM_FIELD, DataType::Int8, true));
        let column =
            FixedSizeListArray::try_new(Arc::clone(&inside), 6, Arc::new(genotypes.finish()), None)
                .expect("the genotypes with a null");
        *parts.column(GTS_COLUMN) = (
            Field::new(GTS_COLUMN, DataType::FixedSizeList(inside, 6), false),
            Arc::new(column),
        );
        let error = refused_at_the_block(parts.written());
        let Error::VarsNullValue { column, var } = error else {
            panic!("the file whose last variant has a null genotype gave {error}");
        };
        assert_eq!((column, var), ("gts", 4));
    }

    /// The blocks of a vars file are its batches, whose size the file was
    /// written with, so a block the machine has no memory for tells the
    /// caller to write the file again with a smaller `num_vars_per_block`
    /// and not to ask for fewer variants in a block, which changes nothing
    /// of what the reader builds. A footer that says the file holds more
    /// variants than this machine counts is the same error.
    #[test]
    fn a_batch_the_machine_has_no_memory_for_says_to_write_the_file_again() {
        let metadata = metadata_of_cases();
        let message = block_too_large(1000, &metadata).to_string();
        assert!(message.contains("1000 variants"), "{message}");
        assert!(message.contains("3 individuals"), "{message}");
        assert!(message.contains("num_vars_per_block"), "{message}");
        assert!(message.contains("write the file again"), "{message}");

        let batches = [
            BatchInfo {
                num_vars: usize::MAX,
                regions: Vec::new(),
            },
            BatchInfo {
                num_vars: 2,
                regions: Vec::new(),
            },
        ];
        let error = match num_vars_of_the_file(&batches, &metadata) {
            Ok(num_vars) => panic!("the file says it holds {num_vars} variants"),
            Err(error) => error,
        };
        let Error::BlockTooLarge {
            num_vars_per_block,
            size,
            ..
        } = error
        else {
            panic!("the footer of more variants than a `usize` counts gave {error}");
        };
        assert_eq!((num_vars_per_block, size), (2, BlockSize::FixedByAFile));
    }

    /// popnei writes `chrom`, `pos`, `alleles` and `gts` as columns with no
    /// nulls, and arrow-rs refuses a batch of such a column that holds one
    /// before popnei sees it: a null there is a batch that could not be
    /// read, with what arrow-rs said, and not the error with the column and
    /// the variant, which is for a file whose schema declares the column
    /// nullable.
    #[test]
    fn a_null_in_a_column_whose_schema_says_it_has_none_is_a_batch_that_was_not_read() {
        let parts = FileParts::of_cases();
        let popnei = parts.popnei.clone().expect("the `popnei` key");
        // The schema of the file is the one popnei writes, where `pos` has
        // no nulls.
        let file = Arc::new(
            Schema::new(
                parts
                    .columns
                    .iter()
                    .map(|(field, _)| field.clone())
                    .collect::<Vec<Field>>(),
            )
            .with_metadata(HashMap::from([(POPNEI_KEY.to_owned(), popnei)])),
        );
        assert!(
            !file.field(1).is_nullable(),
            "the `pos` column of popnei is written with no nulls"
        );
        // The batch is built with `pos` declared nullable, because arrow-rs
        // builds no batch of a column that its schema says has no nulls and
        // that holds one. The file is written with the schema above, so it
        // is the file another program makes when it writes a null into such
        // a column.
        let arrays: Vec<ArrayRef> = parts
            .columns
            .iter()
            .enumerate()
            .map(|(place, (_, array))| match place {
                1 => Arc::new(UInt64Array::from(vec![
                    Some(100),
                    None,
                    Some(300),
                    Some(400),
                ])),
                _ => Arc::clone(array),
            })
            .collect();
        let mut nullable: Vec<Field> = file.fields().iter().map(|f| f.as_ref().clone()).collect();
        nullable[1] = Field::new(POS_COLUMN, DataType::UInt64, true);
        let batch = RecordBatch::try_new(Arc::new(Schema::new(nullable)), arrays)
            .expect("the batch whose `pos` is declared nullable");
        let mut writer = FileWriter::try_new(Vec::new(), &file).expect("the file was started");
        writer.write(&batch).expect("the batch was written");
        writer.write_metadata(
            POPNEI_BATCHES_KEY,
            batches_as_json(&[BatchInfo {
                num_vars: 4,
                regions: regions_of(&CASES),
            }]),
        );
        let bytes = writer.into_inner().expect("the file was finished");

        let error = refused_at_the_block(bytes);

        let Error::VarsBatchNotRead { batch, problem } = &error else {
            panic!("the null in a column with no nulls gave {error}");
        };
        assert_eq!(*batch, 1);
        assert!(problem.contains("pos"), "{problem}");
    }

    /// A quality that is a value and is not a finite number is refused with
    /// the value and the variant: a NaN in the column of a block is what
    /// says that the variant has no quality, so a NaN that is a value in
    /// the file would be read as a variant that has none, and an infinite
    /// quality is a probability of no variant of 0, which is not what phred
    /// scaling says. The VCF reader refuses both.
    #[test]
    fn a_quality_that_is_a_value_and_is_not_finite_is_refused_with_the_variant() {
        for (row, quality) in [(0, f32::NAN), (2, f32::INFINITY), (3, f32::NEG_INFINITY)] {
            let mut parts = FileParts::of_cases();
            let mut qualities = vec![Some(29.5), None, Some(67.0), Some(47.0)];
            qualities[row] = Some(quality);
            *parts.column(QUAL_COLUMN) = (
                Field::new(QUAL_COLUMN, DataType::Float32, true),
                Arc::new(Float32Array::from(qualities)),
            );

            let error = refused_at_the_block(parts.written());

            let Error::VarsQualityNotFinite { found, var } = error else {
                panic!("the file whose quality is {quality} gave {error}");
            };
            assert!(!found.is_finite(), "the error says the quality is {found}");
            assert_eq!(var, counted_from_one(row));
        }
    }

    /// A null `id` is the empty id and a null `qual` is a variant with no
    /// quality, which a block holds as a NaN: both are what any other
    /// program that opens the file writes for a value that is not there, and
    /// neither is an error.
    #[test]
    fn a_null_id_is_the_empty_id_and_a_null_quality_is_not_a_number() {
        let (blocks, _) = blocks_read(cases_written_in_batches_of(4), Needs::ALL);

        let ids = blocks[0].id.as_ref().expect("the ids of the block");
        assert_eq!(ids, &["rs1", "", "", ""]);
        let quals = blocks[0].qual.as_ref().expect("the qualities of the block");
        assert!(quals[1].is_nan(), "the second quality is {}", quals[1]);
        // The same qualities with the NaN of the variant that has none as a
        // `None`, which is how the quality of a variant is compared without
        // comparing two NaNs.
        let there: Vec<Option<f32>> = quals
            .iter()
            .map(|quality| (!quality.is_nan()).then_some(*quality))
            .collect();
        assert_eq!(there, [Some(29.5), None, Some(67.0), Some(47.0)]);
    }

    /// A batch that does not hold the variants its entry of the footer gives
    /// is refused when it is read, with the batch and both counts, so that
    /// the number of variants the file announces, which is read from those
    /// entries with no batch read, is never a wrong one that goes unnoticed.
    #[test]
    fn a_batch_of_another_number_of_variants_than_its_entry_of_the_footer_is_refused() {
        let mut parts = FileParts::of_cases();
        parts.popnei_batches = Some(batches_as_json(&[BatchInfo {
            num_vars: 3,
            regions: vec![Region {
                chrom: "chr1".to_owned(),
                min_pos: 100,
                max_pos: 300,
            }],
        }]));

        let error = refused_at_the_block(parts.written());

        let Error::VarsBatchNumVars {
            batch,
            found,
            expected,
        } = error
        else {
            panic!("the batch of four variants whose entry says three gave {error}");
        };
        assert_eq!((batch, found, expected), (1, 4, 3));
    }

    /// No build of popnei carries the zstd crate, so a file whose buffers
    /// are compressed with zstd is refused: it opens, because the two keys
    /// are not compressed, and the error comes with the first block, because
    /// arrow decompresses a batch when it reads it.
    #[test]
    fn a_file_compressed_with_zstd_opens_and_gives_its_error_at_the_first_block() {
        let mut reader = zstd_vars();

        // The two keys of the file are read although no batch can be.
        assert_eq!(reader.num_vars(), 4);
        assert_eq!(reader.metadata().individuals, cases_individuals());

        let error = match reader.next_block() {
            Ok(block) => panic!("the file compressed with zstd gave {block:?}"),
            Err(error) => error,
        };
        assert!(
            matches!(error, Error::VarsZstd),
            "the file compressed with zstd gave {error}"
        );
        let message = error.to_string();
        assert!(message.contains("zstd"), "{message}");
        assert!(message.contains("lz4"), "{message}");
        assert!(
            reader
                .next_block()
                .expect("the reader gave its error once")
                .is_none()
        );
    }

    /// A batch whose bytes are not what the file says they are is a file
    /// that was damaged after it was written, which the reader refuses with
    /// the batch and what arrow-rs said, instead of giving the variants it
    /// can still read.
    #[test]
    fn a_batch_that_arrow_cannot_decode_is_a_file_that_was_damaged() {
        let mut bytes = cases_written_in_batches_of(4);
        let at = opened(bytes.clone())
            .expect("the file is a vars file")
            .blocks[0];
        // The bytes of the message of the batch, which say where its buffers
        // are and how long they are.
        let start = usize::try_from(at.offset).expect("the offset of the batch");
        for byte in 8..24 {
            if let Some(place) = bytes.get_mut(start.saturating_add(byte)) {
                *place = 0xff;
            }
        }

        let error = refused_at_the_block(bytes);

        let Error::VarsBatchNotRead { batch, problem } = &error else {
            panic!("the file whose batch was damaged gave {error}");
        };
        assert_eq!(*batch, 1);
        assert!(!problem.is_empty(), "the error says nothing");
        let message = error.to_string();
        assert!(message.contains("damaged"), "{message}");
    }

    /// How many bytes one field node of the message of a batch takes: two
    /// 64 bit numbers, how many values its column holds and how many of
    /// those are nulls.
    const FIELD_NODE_BYTES: usize = 16;

    /// How many bytes the first of that pair takes.
    const FIELD_NODE_LENGTH_BYTES: usize = 8;

    /// Where the field nodes of the message of a batch start in the bytes
    /// of the file, and what each one says: a pair of the values its column
    /// holds and the nulls among them.
    ///
    /// The nodes are found by the bytes of the whole vector, and not by
    /// those of one node, because a buffer of the same message is a pair of
    /// 64 bit numbers too and can hold the same pair.
    fn nodes_of_the_message(bytes: &[u8], at: BatchAt) -> (usize, Vec<(i64, i64)>) {
        let start = usize::try_from(at.offset).expect("the offset of the batch");
        let end = start
            .checked_add(usize::try_from(at.metadata_len).expect("the message of the batch"))
            .expect("where the message of the batch ends");
        let message = bytes.get(start..end).expect("the message of the batch");
        let starts_at = match message.get(..CONTINUATION_MARK.len()) {
            Some(mark) if mark == CONTINUATION_MARK => MESSAGE_START_BYTES,
            Some(_) | None => CONTINUATION_MARK.len(),
        };
        let parsed = root_as_message(message.get(starts_at..).expect("the message"))
            .expect("the message of the batch is one of arrow");
        let batch = parsed
            .header_as_record_batch()
            .expect("the message is one of a batch");
        let nodes: Vec<(i64, i64)> = batch
            .nodes()
            .expect("the field nodes of the batch")
            .iter()
            .map(|node| (node.length(), node.null_count()))
            .collect();
        let pattern: Vec<u8> = nodes
            .iter()
            .flat_map(|(length, nulls)| length.to_le_bytes().into_iter().chain(nulls.to_le_bytes()))
            .collect();
        let where_they_are: Vec<usize> = message
            .windows(pattern.len())
            .enumerate()
            .filter(|(_, window)| *window == pattern.as_slice())
            .map(|(at, _)| at)
            .collect();
        assert_eq!(
            where_they_are.len(),
            1,
            "the field nodes of the message are in {count} places of it",
            count = where_they_are.len()
        );
        let nodes_at = start
            .checked_add(*where_they_are.first().expect("the field nodes"))
            .expect("where the field nodes are in the file");
        (nodes_at, nodes)
    }

    /// Those bytes with the field node at `which` saying that its column
    /// holds `says` values, which is how a test damages one column of one
    /// batch and leaves the rest of the file as it was.
    fn a_node_that_says(mut bytes: Vec<u8>, at: BatchAt, which: usize, says: i64) -> Vec<u8> {
        let (nodes_at, nodes) = nodes_of_the_message(&bytes, at);
        assert!(
            which < nodes.len(),
            "the message has {count} field nodes",
            count = nodes.len()
        );
        let length_at = nodes_at
            .checked_add(which.checked_mul(FIELD_NODE_BYTES).expect("the field node"))
            .expect("where the field node is in the file");
        let ends_at = length_at
            .checked_add(FIELD_NODE_LENGTH_BYTES)
            .expect("where the length of the field node ends");
        bytes
            .get_mut(length_at..ends_at)
            .expect("the length of the field node")
            .copy_from_slice(&says.to_le_bytes());
        bytes
    }

    /// The bound on the values a column says it holds is what keeps a
    /// damaged length from reaching arrow-rs, which builds an array of that
    /// length: the alleles of the genotypes of a batch are its rows times
    /// the individuals times the ploidy, and a message that says one more
    /// is refused.
    #[test]
    fn a_column_that_says_more_values_than_it_can_hold_is_a_file_that_was_damaged() {
        assert_eq!(MAX_VALUES_OF_A_LIST, u64::from(i32::MAX.unsigned_abs()));
        let bytes = cases_written_in_batches_of(4);
        let at = opened(bytes.clone())
            .expect("the file is a vars file")
            .blocks[0];
        // The eight field nodes of the batch, of which the last is the
        // alleles of the genotypes: 4 variants of 3 diploid individuals
        // hold 24 alleles.
        let (_, nodes) = nodes_of_the_message(&bytes, at);
        assert_eq!(
            nodes,
            [
                (4, 0),
                (4, 0),
                (4, 3),
                (4, 0),
                (8, 0),
                (4, 1),
                (4, 0),
                (24, 0)
            ]
        );

        let error = refused_at_the_block(a_node_that_says(bytes, at, 7, 25));

        let Error::VarsBatchNotRead { batch, problem } = &error else {
            panic!("the file whose node says 25 alleles gave {error}");
        };
        assert_eq!(*batch, 1);
        assert!(problem.contains("25 values"), "{problem}");
        assert!(problem.contains("holds 24"), "{problem}");
    }

    /// A column that comes after one whose type popnei does not walk keeps
    /// a bound: popnei cannot say how many field nodes such a column takes,
    /// so it cannot say which column the nodes after it belong to and the
    /// schema gives them none. What is left is the body of the batch, which
    /// holds at most 255 bytes for each of its own and a value for each bit
    /// of those.
    ///
    /// Without that bound nothing holds those columns, since the rows of a
    /// batch are what its entry of the footer says and no byte of the file
    /// bounds them: arrow-rs reaches `integer overflow computing expected
    /// number of expected values in FixedListSize`, an `expect` that panics
    /// in a release build too, which `catch_unwind` holds natively and
    /// which ends a browser tab.
    #[test]
    fn a_column_after_one_popnei_does_not_walk_is_bounded_by_the_bytes_of_the_batch() {
        // A column of a type popnei does not walk, before the six it knows,
        // which the reader ignores and which stops the walk of the schema.
        // pyarrow writes one for any categorical of pandas.
        let mut parts = FileParts::of_cases();
        let names: DictionaryArray<Int32Type> =
            vec!["one", "two", "one", "two"].into_iter().collect();
        parts.columns.insert(
            0,
            (
                Field::new("kind", names.data_type().clone(), false),
                Arc::new(names),
            ),
        );
        let bytes = parts.written();
        let at = opened(bytes.clone())
            .expect("the file with a dictionary column is read")
            .blocks[0];
        // The nine field nodes, of which the first is the dictionary that
        // stops the walk and the last the alleles of the genotypes.
        let (_, nodes) = nodes_of_the_message(&bytes, at);
        assert_eq!(nodes.len(), 9);
        assert_eq!(nodes.last(), Some(&(24, 0)));

        let error = refused_at_the_block(a_node_that_says(bytes, at, 8, 3_074_457_345_618_258_603));

        let Error::VarsBatchNotRead { batch, problem } = &error else {
            panic!("the file whose node says 3074457345618258603 alleles gave {error}");
        };
        assert_eq!(*batch, 1);
        assert!(problem.contains("3074457345618258603 values"), "{problem}");
        // The batch was refused by its own bytes and never handed to
        // arrow-rs, which is the whole of the defence in a browser tab.
        assert!(
            !problem.contains("arrow-rs"),
            "the length reached arrow-rs: {problem}"
        );
    }

    /// The walk of the columns gives one bound for each field node the IPC
    /// format lays out, in its order, for the types a file another program
    /// wrote can hold and popnei's own writer never makes: a struct, whose
    /// columns hold the rows of the struct itself; a large list, whose
    /// values its 64 bit offsets count and nothing in the schema bounds; a
    /// map, which is a list of its entries; and a column of a type popnei
    /// does not walk, which ends the list there because popnei cannot say
    /// how many nodes it takes.
    ///
    /// popnei's own six columns reach four of the arms and no test reaches
    /// the others, so a walk that counted the nodes of one of them wrong
    /// would put every bound after it on another column.
    #[test]
    fn the_walk_gives_one_bound_for_each_field_node_of_the_types_popnei_does_not_write() {
        let inside = Fields::from(vec![
            Field::new("one", DataType::UInt64, false),
            Field::new("two", DataType::Utf8, true),
        ]);
        let fields = Fields::from(vec![
            Field::new("plain", DataType::Float32, true),
            Field::new("both", DataType::Struct(inside), false),
            Field::new(
                "many",
                DataType::LargeList(Arc::new(Field::new(ITEM_FIELD, DataType::Int8, false))),
                false,
            ),
            Field::new(
                "pairs",
                DataType::Map(
                    Arc::new(Field::new(
                        "entries",
                        DataType::Struct(Fields::from(vec![
                            Field::new("key", DataType::Utf8, false),
                            Field::new("value", DataType::Int8, true),
                        ])),
                        false,
                    )),
                    false,
                ),
                false,
            ),
            Field::new(
                "six",
                DataType::FixedSizeList(Arc::new(gts_field()), 6),
                false,
            ),
        ]);

        let holds = what_the_nodes_hold(&fields, 100);

        // plain; both and its two columns, which hold the rows of the
        // struct; many and its values, which nothing in the schema bounds;
        // pairs, its entries and the key and the value of an entry, which
        // its 32 bit offsets bound; six and its 600 alleles.
        assert_eq!(
            holds,
            [
                NodeHolds::Values(100),
                NodeHolds::Values(100),
                NodeHolds::Values(100),
                NodeHolds::Values(100),
                NodeHolds::Values(100),
                NodeHolds::NotBounded,
                NodeHolds::Values(100),
                NodeHolds::Values(MAX_VALUES_OF_A_LIST),
                NodeHolds::Values(MAX_VALUES_OF_A_LIST),
                NodeHolds::Values(MAX_VALUES_OF_A_LIST),
                NodeHolds::Values(100),
                NodeHolds::Values(600),
            ]
        );

        // A column of a type popnei does not walk ends the list where it
        // is: it keeps its own bound, and nothing after it has one.
        let with_a_dictionary = Fields::from(vec![
            Field::new("plain", DataType::Float32, true),
            Field::new(
                "kind",
                DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8)),
                false,
            ),
            Field::new("after", DataType::UInt64, false),
        ]);

        assert_eq!(
            what_the_nodes_hold(&with_a_dictionary, 100),
            [NodeHolds::Values(100), NodeHolds::Values(100)]
        );
    }

    /// After an error a reader gives no block at every call, as
    /// `docs/specs/block.md` asks: one that went on would give the variants
    /// that follow a batch it could not read as if nothing had happened.
    #[test]
    fn a_reader_that_failed_gives_no_block_at_every_call_after() {
        let mut parts = FileParts::of_cases();
        let positions: ArrayRef = Arc::new(UInt64Array::from(vec![
            Some(100),
            None,
            Some(300),
            Some(400),
        ]));
        *parts.column(POS_COLUMN) = (Field::new(POS_COLUMN, DataType::UInt64, true), positions);
        parts.num_batches = 2;
        parts.popnei_batches = Some(batches_as_json(&[
            BatchInfo {
                num_vars: 4,
                regions: Vec::new(),
            },
            BatchInfo {
                num_vars: 4,
                regions: Vec::new(),
            },
        ]));
        let mut reader = opened(parts.written()).expect("the file is a vars file");

        let error = match reader.next_block() {
            Ok(block) => panic!("the batch with a null position gave {block:?}"),
            Err(error) => error,
        };
        assert!(
            matches!(error, Error::VarsNullValue { .. }),
            "the batch with a null position gave {error}"
        );

        assert!(
            reader
                .next_block()
                .expect("the reader gave its error once")
                .is_none()
        );
        assert!(
            reader
                .next_block()
                .expect("and at every call after")
                .is_none()
        );
    }

    /// arrow-rs reads the four bytes of the mark of a continuation and the
    /// four of the length of the message by their place, so a batch of
    /// fewer bytes than those eight panics inside it. The reader refuses
    /// such a batch before arrow-rs is given it, with the batch and what is
    /// wrong.
    #[test]
    fn a_batch_of_fewer_bytes_than_the_start_of_a_message_is_refused_before_arrow_reads_it() {
        // The entry of the footer says the batch is 4 bytes of message and
        // no body, where a message starts with 8.
        let bytes = batch_entry_changed(cases_written_in_batches_of(4), 4, 0);

        let error = refused_at_the_block(bytes);

        let Error::VarsBatchNotRead { batch, problem } = &error else {
            panic!("the batch of four bytes gave {error}");
        };
        assert_eq!(*batch, 1);
        assert!(problem.contains('8'), "{problem}");
    }

    /// A variant with no alleles at all, and one with no genotypes, are the
    /// list itself being a null and not a value inside it: a reader that
    /// looked only inside would give a variant with no allele and
    /// genotypes read out of the padding of the column.
    #[test]
    fn a_null_row_of_the_alleles_and_of_the_genotypes_is_an_error_of_that_column() {
        let mut parts = FileParts::of_cases();
        // The values of the list hold no null, as popnei writes them; the
        // null is the row of the column, which is a variant with no
        // alleles at all.
        let mut alleles = ListBuilder::new(StringBuilder::new()).with_field(Arc::new(Field::new(
            ITEM_FIELD,
            DataType::Utf8,
            false,
        )));
        for row in 0..4 {
            match row {
                1 => alleles.append_null(),
                _ => {
                    alleles.values().append_value("A");
                    alleles.values().append_value("T");
                    alleles.append(true);
                }
            }
        }
        *parts.column(ALLELES_COLUMN) = (
            Field::new(ALLELES_COLUMN, alleles_type(), true),
            Arc::new(alleles.finish()),
        );
        let error = refused_at_the_block(parts.written());
        let Error::VarsNullValue { column, var } = error else {
            panic!("the file whose second variant has no alleles gave {error}");
        };
        assert_eq!((column, var), ("alleles", 2));

        let mut parts = FileParts::of_cases();
        let inside = Arc::new(Field::new(ITEM_FIELD, DataType::Int8, false));
        let there: NullBuffer = [true, true, true, false].into_iter().collect();
        let column = FixedSizeListArray::try_new(
            Arc::clone(&inside),
            6,
            Arc::new(Int8Array::from(vec![0_i8; 24])),
            Some(there),
        )
        .expect("the genotypes whose last variant has none");
        *parts.column(GTS_COLUMN) = (
            Field::new(GTS_COLUMN, DataType::FixedSizeList(inside, 6), true),
            Arc::new(column),
        );
        let error = refused_at_the_block(parts.written());
        let Error::VarsNullValue { column, var } = error else {
            panic!("the file whose last variant has no genotypes gave {error}");
        };
        assert_eq!((column, var), ("gts", 4));
    }

    /// The variant of the error of a null is counted from 1 over the whole
    /// file and not inside its batch, and the batch of an error is counted
    /// from 1 over the file: both are read in the second batch of a file of
    /// two, where a count that started again at each batch is seen.
    #[test]
    fn the_variant_and_the_batch_of_an_error_are_counted_over_the_whole_file() {
        // The third variant of the second batch of two batches of four,
        // which is the variant 7 of the file.
        let with_a_null = UInt64Array::from(vec![Some(100), Some(200), None, Some(400)]);
        let error = refused_at_the_block(cases_in_two_batches(with_a_null, [4, 4]));
        let Error::VarsNullValue { column, var } = error else {
            panic!("the null position of the second batch gave {error}");
        };
        assert_eq!((column, var), ("pos", 7));

        // The entry of the second batch says three variants and the batch
        // holds four.
        let positions = UInt64Array::from(vec![100, 200, 300, 400]);
        let error = refused_at_the_block(cases_in_two_batches(positions, [4, 3]));
        let Error::VarsBatchNumVars {
            batch,
            found,
            expected,
        } = error
        else {
            panic!("the second batch of another number of variants gave {error}");
        };
        assert_eq!((batch, found, expected), (2, 4, 3));
    }

    /// A file with two columns of one name is read with the first of them,
    /// as any program that asks a table for a column by its name reads it.
    #[test]
    fn a_file_with_two_columns_of_one_name_is_read_with_the_first() {
        let mut parts = FileParts::of_cases();
        let (field, array) = parts.column(POS_COLUMN).clone();
        // A second `pos` column, of other positions, after the first.
        let other: ArrayRef = Arc::new(UInt64Array::from(vec![900, 901, 902, 903]));
        parts.columns.push((field, other));
        assert_eq!(parts.columns.len(), 7);

        let (blocks, _) = blocks_read(parts.written(), Needs::ALL);

        let positions = blocks[0].pos.as_deref().expect("the positions");
        let first = array
            .as_any()
            .downcast_ref::<UInt64Array>()
            .expect("the first `pos` column")
            .values()
            .to_vec();
        assert_eq!(positions, first);
        assert_eq!(positions, [100, 200, 300, 400]);
    }

    /// The columns of a file are read at the place they have in it and not
    /// at the place the table of the spec gives them, values and all: a
    /// reader that took them in the order of the spec would give the
    /// positions of one column and the ids of another.
    #[test]
    fn a_file_whose_columns_are_in_another_order_is_read_with_the_values_of_each() {
        let mut parts = FileParts::of_cases();
        parts.columns.reverse();
        assert_eq!(parts.columns[0].0.name(), GTS_COLUMN);

        let (blocks, chroms) = blocks_read(parts.written(), Needs::ALL);

        assert_eq!(num_vars_of(&blocks), [4]);
        let expected: Vec<ReadRow> = CASES.iter().map(row_read).collect();
        assert_eq!(rows_of(&blocks, &chroms), expected);
    }

    /// How many panics the net of the reader caught in the sweep over four
    /// values of each byte of the file, on 24 September 2026 with arrow-rs
    /// 60: one about a buffer that is not long enough for the values its
    /// message says, and one about a buffer whose length is not a whole
    /// number of the values of its column. Both are asserts of arrow-rs
    /// that only the walk of the schema its decoder does would see.
    ///
    /// It was 14 until the bound on the values a column says it holds
    /// became the column itself, on 24 September 2026: a column of a batch
    /// holds its rows, which is far below the bits of the batch that the
    /// bound was before, so three more of these files are refused before
    /// arrow-rs reads them.
    const PANICS_CAUGHT: u64 = 11;

    /// The variants of the file in those bytes, with every field, or the
    /// error it gave: what a sweep over a damaged file reads.
    fn rows_or_error(bytes: Vec<u8>) -> Result<Vec<ReadRow>> {
        let mut reader = opened(bytes)?;
        reader.set_needs(Needs::ALL);
        let blocks = blocks_of(&mut reader)?;
        for block in &blocks {
            block.check()?;
        }
        Ok(rows_of(&blocks, reader.chroms()))
    }

    /// What a sweep over the bytes of a file found: how many of the files
    /// it made gave an error, how many were read as the whole file and how
    /// many were read as another file with no error, and where arrow-rs
    /// panicked inside the reader, which the net of the reader turns into
    /// an error.
    #[derive(Debug)]
    struct Sweep {
        errors: u64,
        the_same: u64,
        others: u64,
        /// The place of each panic, `file:line`, with how many times it was
        /// reached.
        panics: Vec<(String, u64)>,
    }

    impl Sweep {
        /// How many panics the net of the reader caught.
        fn caught(&self) -> u64 {
            self.panics.iter().map(|(_, count)| *count).sum()
        }

        /// That every panic was one of arrow-rs, which the reader catches:
        /// one of popnei's own code is a defect and not a net that held.
        fn panics_are_of_arrow(&self) {
            for (place, count) in &self.panics {
                assert!(
                    !place.contains("crates/popnei"),
                    "popnei panicked at {place}, {count} times"
                );
            }
        }
    }

    /// Every byte of `whole` set in turn to each of `values`, read with
    /// every field and counted.
    ///
    /// A file that gives no error and other variants is not a failure here:
    /// a byte of a compressed buffer that decompresses into other genotypes
    /// is what a checksum of the format would catch, which the file has
    /// not, and `docs/specs/io_vars.md` says so. What must not happen is a
    /// panic that comes out of the reader, which fails the test where it
    /// happens, or an abort, which kills the test binary.
    ///
    /// The hook a panic runs is the one of the whole process, and cargo
    /// runs the other tests of this binary in threads of that process
    /// while the sweep runs, so a test that fails elsewhere at that moment
    /// panics into this hook: it was counted as a panic the reader caught,
    /// the sweep failed as well, and its message named the file of the
    /// other test. Two reviewers got `popnei panicked at
    /// crates/popnei/src/filters.rs:1588:9, 1 times` that way on 21
    /// September 2026, by breaking a test of the filters. So each panic is
    /// counted under the thread it happened in, and the sweep keeps the
    /// ones of its own: the reader of a vars file uses no threads, so every
    /// panic of the files this makes is on this thread.
    fn swept(whole: &[u8], values: &[u8], expected: &[ReadRow]) -> Sweep {
        let seen: Arc<Mutex<HashMap<(ThreadId, String), u64>>> =
            Arc::new(Mutex::new(HashMap::new()));
        let hook = std::panic::take_hook();
        let writing = Arc::clone(&seen);
        std::panic::set_hook(Box::new(move |panic| {
            let place = panic
                .location()
                .map_or_else(|| "nowhere".to_owned(), |at| format!("{at}"));
            if let Ok(mut seen) = writing.lock() {
                let counted = seen
                    .entry((std::thread::current().id(), place))
                    .or_insert(0);
                *counted = counted.saturating_add(1);
            }
        }));
        let sweeping = std::thread::current().id();
        let mut errors: u64 = 0;
        let mut the_same: u64 = 0;
        let mut others: u64 = 0;
        for (at, byte) in whole.iter().enumerate() {
            for value in values {
                if value == byte {
                    continue;
                }
                let mut bytes = whole.to_vec();
                bytes[at] = *value;
                match rows_or_error(bytes) {
                    Err(_) => errors = errors.saturating_add(1),
                    Ok(rows) if rows == expected => the_same = the_same.saturating_add(1),
                    Ok(_) => others = others.saturating_add(1),
                }
            }
        }
        std::panic::set_hook(hook);
        let panics = seen
            .lock()
            .map(|seen| {
                let mut places: Vec<(String, u64)> = seen
                    .iter()
                    .filter(|((thread, _), _)| *thread == sweeping)
                    .map(|((_, at), count)| (at.clone(), *count))
                    .collect();
                places.sort();
                places
            })
            .unwrap_or_default();
        Sweep {
            errors,
            the_same,
            others,
            panics,
        }
    }

    /// No change of one byte of a vars file reaches a panic or an abort.
    ///
    /// A reviewer of the work package changed every byte of such a file to
    /// each of the 255 other values on 21 September 2026: of 851190 files,
    /// 70243 panicked inside arrow-rs and 2854 ended the process, which
    /// `Vec::with_capacity` does with the length a damaged compressed
    /// buffer says. The reader checks the message of a batch against the
    /// bytes that came with it before arrow-rs reads any of them, and holds
    /// what that does not see in `catch_unwind`.
    ///
    /// On 24 September 2026 the four values of each byte gave 16296 files
    /// in 0.13 s: 6946 errors, 8988 read as the whole file, 362 read as
    /// another file with no error, 11 panics caught and no abort. The check
    /// of an allele below the missing one moved 24 files from the third
    /// count to the first.
    ///
    /// Every byte is set to four values here, the low bit and the high bit
    /// flipped, 0 and 255; the sweep over all 255 is the ignored test that
    /// follows.
    ///
    /// The panics the net catches are counted and not refused: what is left
    /// of them is an assert of arrow-rs about a buffer that its message
    /// does not fit, which only the walk of the schema that its decoder
    /// does would see.
    #[test]
    fn no_change_of_one_byte_of_a_vars_file_reaches_a_panic() {
        let whole = cases_written_in_batches_of(3);
        let expected = rows_or_error(whole.clone()).expect("the whole file is read");
        assert_eq!(expected.len(), 4);

        let values = [0x00, 0xff];
        let mut sweep = swept(&whole, &values, &expected);
        // The two values that depend on the byte, which `swept` takes as
        // they are, so the sweep is run again for each of them.
        for bit in [0x01_u8, 0x80] {
            for (at, byte) in whole.iter().enumerate() {
                let mut bytes = whole.clone();
                bytes[at] = byte ^ bit;
                match rows_or_error(bytes) {
                    Err(_) => sweep.errors = sweep.errors.saturating_add(1),
                    Ok(rows) if rows == expected => {
                        sweep.the_same = sweep.the_same.saturating_add(1);
                    }
                    Ok(_) => sweep.others = sweep.others.saturating_add(1),
                }
            }
        }

        // A test in which nothing is read, or nothing refused, is one that
        // cannot fail.
        assert!(sweep.errors > 100, "{} errors", sweep.errors);
        assert!(sweep.the_same > 100, "{} files read whole", sweep.the_same);
        sweep.panics_are_of_arrow();
        assert!(
            sweep.caught() <= PANICS_CAUGHT,
            "{} panics were caught, and {PANICS_CAUGHT} were counted on 21 September 2026: {:?}",
            sweep.caught(),
            sweep.panics
        );
        println!("{sweep:?}");
    }

    /// The same sweep with every byte set to each of the 255 other values.
    /// It is run by hand: 1299990 files and 9.2 s in the profile of the
    /// tests on the owner's Apple M5 Pro on 24 September 2026, where it
    /// gave 553104 errors, 726033 files read as the whole one, 20853 read
    /// as another file with no error, which is what a checksum of the
    /// format would catch and nothing else does, 1939 panics of arrow-rs
    /// that the net caught and no abort. The panics were 2783 until the
    /// bound on the values a column says it holds became the column
    /// itself. The check of an allele below the
    /// missing one moved 3049 files, 13 in 100 of the 23902 that were read
    /// as another file before it, from the third count to the first.
    ///
    ///     cargo test -p popnei --lib \
    ///         no_change_of_any_byte_of_a_vars_file -- --ignored
    #[test]
    #[ignore = "1299990 files; the sweep over four values of each byte is the one that runs with the suite"]
    fn no_change_of_any_byte_of_a_vars_file_reaches_a_panic() {
        let whole = cases_written_in_batches_of(3);
        let expected = rows_or_error(whole.clone()).expect("the whole file is read");
        let values: Vec<u8> = (0..=255).collect();

        let sweep = swept(&whole, &values, &expected);

        assert!(sweep.errors > 10000, "{} errors", sweep.errors);
        sweep.panics_are_of_arrow();
        // What the reader gives with no error and other variants: a byte of
        // a compressed buffer that decompresses into other genotypes, which
        // only a checksum of the format would catch.
        assert!(
            sweep.others < sweep.errors,
            "{} files read as another file",
            sweep.others
        );
        println!("{sweep:?}");
    }

    /// Where the length of the last buffer of the batch `batch` of the
    /// file, counted from 0, is declared: the eight bytes arrow writes
    /// before the bytes it compressed, which say how long that buffer is
    /// once it is decompressed. The last buffer of a batch of the six
    /// columns holds the genotypes, which is the large one.
    fn the_length_of_the_genotypes(bytes: &[u8], batch: usize) -> usize {
        let at = opened(bytes.to_vec())
            .expect("the file is a vars file")
            .blocks[batch];
        let offset = usize::try_from(at.offset).expect("the offset of the batch");
        let metadata_len = usize::try_from(at.metadata_len).expect("the message of the batch");
        let starts_at = offset.saturating_add(MESSAGE_START_BYTES);
        let message = root_as_message(&bytes[starts_at..]).expect("the message of the batch");
        let batch = message
            .header_as_record_batch()
            .expect("the message is one of a batch");
        let last = batch
            .buffers()
            .expect("the buffers of the batch")
            .iter()
            .next_back()
            .expect("the last buffer");
        let in_the_body = usize::try_from(last.offset()).expect("the offset of the buffer");
        offset
            .saturating_add(metadata_len)
            .saturating_add(in_the_body)
    }

    /// The byte of the first allele of the batch `batch` of the file,
    /// counted from 0.
    ///
    /// lz4 makes the genotypes of a file of four variants no smaller, so
    /// arrow writes that buffer as it is and says -1 where its length would
    /// be: each allele of the batch is one byte of the file, in the order
    /// the variants are in, and the test that fails here is one whose file
    /// grew until lz4 had something to take out of it.
    fn the_first_allele_of(bytes: &[u8], batch: usize) -> usize {
        let at = the_length_of_the_genotypes(bytes, batch);
        let says = <[u8; UNCOMPRESSED_LENGTH_BYTES]>::try_from(
            &bytes[at..at.saturating_add(UNCOMPRESSED_LENGTH_BYTES)],
        )
        .map(i64::from_le_bytes)
        .expect("what the buffer of the genotypes says");
        assert_eq!(
            says, NOT_COMPRESSED,
            "the genotypes of the batch {batch} are compressed"
        );
        at.saturating_add(UNCOMPRESSED_LENGTH_BYTES)
    }

    /// An allele below the missing one is neither an allele of a variant
    /// nor a missing genotype, and the reader refuses the batch that holds
    /// one instead of handing it over. The owner decided on 21 September
    /// 2026 that such an allele "is never allowed", after a reviewer wrote
    /// a vars file with popnei, changed one byte of its genotypes to 254
    /// and got the genotype `[-2, 0]` out of `iter_blocks` with no error.
    ///
    /// The genotypes of the four variants of `cases.vcf` are `0/0 0/1 1/1`,
    /// `./. 0/1 ./0`, `1/2 2/1 2/2` and `0/0 0/0 0/0`, two batches of two
    /// variants here, so the last allele of the first batch is of the
    /// variant 2 and the first of the second batch is of the variant 3.
    #[test]
    fn an_allele_below_the_missing_one_is_refused_with_the_allele_and_its_variant() {
        let whole = cases_written_in_batches_of(2);

        for (batch, allele_of_the_batch, allele, var) in [
            (0_usize, 11_usize, -2_i8, 2_u64),
            (1, 0, -3, 3),
            (1, 6, i8::MIN, 4),
        ] {
            let at = the_first_allele_of(&whole, batch).saturating_add(allele_of_the_batch);
            let mut bytes = whole.clone();
            bytes[at] = allele.to_le_bytes()[0];

            let error = refused_at_the_block(bytes);

            let Error::VarsAlleleBelowMissing { found, var: of } = error else {
                panic!("the allele {allele} of the variant {var} gave {error}");
            };
            assert_eq!((found, of), (allele, var));
        }

        // The file as it was written is read, and the missing allele, which
        // the second variant holds three of, is not refused.
        let (blocks, _) = blocks_read(whole, Needs::GTS);
        assert_eq!(
            blocks[0].gts,
            vec![
                0,
                0,
                0,
                1,
                1,
                1,
                MISSING_ALLELE,
                MISSING_ALLELE,
                0,
                1,
                MISSING_ALLELE,
                0
            ]
        );
    }

    /// Those bytes with the length that the buffer of the genotypes
    /// declares changed to `says`.
    fn genotypes_that_say(bytes: &[u8], says: i64) -> Vec<u8> {
        let at = the_length_of_the_genotypes(bytes, 0);
        let mut changed = bytes.to_vec();
        changed[at..at.saturating_add(UNCOMPRESSED_LENGTH_BYTES)]
            .copy_from_slice(&says.to_le_bytes());
        changed
    }

    /// A buffer of a batch says how long it is once it is decompressed, and
    /// arrow-rs asks the machine for that memory before it decompresses it:
    /// under wasm a length of 3000000000 is a trap that ends the tab, and
    /// one of 2000000000 leaves the memory of the tab grown for its life.
    /// So the reader holds each buffer to what its column can hold, which
    /// the schema and the rows of the batch give: the genotypes are the
    /// variants times the individuals times the ploidy bytes and no more.
    #[test]
    fn a_buffer_that_says_it_holds_more_than_its_column_is_refused_before_arrow_reads_it() {
        // A file of a thousand variants, whose genotypes are 6000 bytes and
        // are compressed: arrow writes the buffers it could not make
        // smaller, which the four variants of `cases.vcf` give, as they are
        // and says -1 where their length would be.
        let mut chroms = ChromTable::new();
        let rows: Vec<&Row> = vec![&CASES[0]; 1000];
        let block = block_of(&rows, &mut chroms);
        let reader = GivenBlocks::of(vec![block], chroms);
        let (whole, _) = write_vars(reader, Vec::new(), Some(1000)).expect("the file was written");
        let at = the_length_of_the_genotypes(&whole, 0);
        let says = <[u8; UNCOMPRESSED_LENGTH_BYTES]>::try_from(
            &whole[at..at.saturating_add(UNCOMPRESSED_LENGTH_BYTES)],
        )
        .map(i64::from_le_bytes)
        .expect("what the buffer of the genotypes says");
        // The thousand variants of the three individuals of the ploidy 2.
        assert_eq!(says, 6000);

        for says in [6001, 3_000_000_000] {
            let error = refused_at_the_block(genotypes_that_say(&whole, says));
            let Error::VarsBatchNotRead { batch, problem } = &error else {
                panic!("the buffer that says {says} bytes gave {error}");
            };
            assert_eq!(*batch, 1);
            assert!(problem.contains(&says.to_string()), "{problem}");
            // The error is the one of the reader, given before arrow-rs
            // asked the machine for the memory the buffer says, and not the
            // one arrow-rs gives once it has.
            assert!(problem.contains("its column holds"), "{problem}");
        }

        // And the file whose buffer says what it holds is read.
        let (blocks, _) = blocks_read(whole, Needs::ALL);
        assert_eq!(num_vars_of(&blocks), [1000]);
    }

    /// A batch whose arrays are a window into longer ones is read through
    /// that window: arrow gives an array an offset when it is sliced, and
    /// the values of a column are taken through it and not from the start of
    /// its buffer.
    ///
    /// The decoder of arrow-rs gives the arrays of a batch of a file with no
    /// offset, so the file cannot hold one and this is made at the block
    /// that one batch becomes.
    #[test]
    fn a_batch_whose_arrays_are_a_window_into_longer_ones_is_read_through_it() {
        let parts = FileParts::of_cases();
        let fields: Vec<Field> = parts
            .columns
            .iter()
            .map(|(field, _)| field.clone())
            .collect();
        let arrays: Vec<ArrayRef> = parts
            .columns
            .iter()
            .map(|(_, array)| array.slice(1, 2))
            .collect();
        let batch = RecordBatch::try_new(Arc::new(Schema::new(fields)), arrays)
            .expect("the window is a batch");
        let wanted = projection_of(
            Needs::ALL,
            &VarsColumns {
                chrom: Some(0),
                pos: Some(1),
                id: Some(2),
                alleles: Some(3),
                qual: Some(4),
                gts: Some(5),
            },
        );
        let mut chroms = ChromTable::new();

        let block = block_of_the_batch(
            &batch,
            &wanted,
            &metadata_of_cases(),
            &mut chroms,
            BatchPlace {
                batch: 1,
                vars_before: 0,
                num_vars: 2,
            },
        )
        .expect("the window of two variants is a block");

        block.check().expect("the block of the window");
        assert_eq!(block.num_vars, 2);
        let expected: Vec<ReadRow> = CASES[1..3].iter().map(row_read).collect();
        assert_eq!(rows_of(&[block], &chroms), expected);
    }

    /// The chromosome that every batch of the file of many batches holds,
    /// which the first batch is the first to hold.
    const SHARED_CHROM: &str = "chr0";

    /// One batch of the file of many batches: two variants, the first on
    /// `chrom`, a name that no batch before it held, and the second on the
    /// chromosome every batch holds.
    ///
    /// Every field of the batch depends on which batch it is, so a block that
    /// came out of another batch is not equal to it: the batches of a file
    /// whose variants were all the same would be read in any order and
    /// nothing would say so.
    fn two_variants_of_a_batch(batch: u64, chrom: &str, chroms: &mut ChromTable) -> Block {
        let mut alleles = AllelesColumn::with_num_vars(2).expect("the alleles");
        alleles.push(&["A".to_owned(), "T".to_owned()]);
        alleles.push(&["C".to_owned(), "G".to_owned(), "GG".to_owned()]);
        let of_its_own = chroms.intern(chrom);
        let shared = chroms.intern(SHARED_CHROM);
        let first_allele = i8::try_from(batch % 3).expect("the first allele of the batch");
        let quality = 29.5 + f32::from(u16::try_from(batch).expect("the quality of the batch"));
        Block {
            num_vars: 2,
            num_individuals: CASES_INDIVIDUALS,
            ploidy: CASES_PLOIDY,
            gts: vec![
                first_allele,
                1,
                MISSING,
                1,
                1,
                0,
                1,
                0,
                0,
                0,
                MISSING,
                MISSING,
            ],
            chrom: Some(vec![of_its_own, shared]),
            pos: Some(vec![
                100_u64.saturating_add(batch),
                1000_u64.saturating_add(batch),
            ]),
            // The second variant has no id, which the file holds as a null.
            id: Some(vec![format!("rs{batch}"), String::new()]),
            alleles: Some(alleles),
            // The second variant has no quality, which is a NaN in the block
            // and a null in the file.
            qual: Some(vec![quality, f32::NAN]),
        }
    }

    /// A vars file of `num_batches` batches of two variants, where the name
    /// of the chromosome of the first variant of each batch first appears in
    /// that batch: `chr0` in the first, `chr1` in the second, and so on, so
    /// that the table of the names of a pass is in the order of the batches
    /// and a pass that read them in another order numbered them otherwise.
    fn a_file_of_many_batches(num_batches: u64) -> Vec<u8> {
        let mut chroms = ChromTable::new();
        let mut blocks = Vec::new();
        for batch in 0..num_batches {
            let chrom = format!("chr{batch}");
            blocks.push(two_variants_of_a_batch(batch, &chrom, &mut chroms));
        }
        let (bytes, refused) =
            written_block_by_block(blocks, &chroms, &cases_individuals(), CASES_PLOIDY);
        assert!(refused.is_none(), "the file was written: {refused:?}");
        bytes
    }

    /// Everything a pass over a vars file gave, as values a test compares:
    /// every variant of every block in the order they came out, the numbers
    /// behind their chromosomes block by block, and the names of the table in
    /// the order of their numbers.
    ///
    /// The names go in beside the numbers because the two can be wrong
    /// together: a pass that read the batches in another order would hand the
    /// numbers out in that order, and each variant would still find its own
    /// name behind its own number.
    #[derive(Debug, PartialEq)]
    struct WholePass {
        rows: Vec<ReadRow>,
        chrom_numbers: Vec<Vec<u32>>,
        chrom_names: Vec<String>,
    }

    /// One pass over the vars file in `bytes` with every field asked for.
    fn the_whole_pass(bytes: &[u8]) -> WholePass {
        let mut reader = VarsReader::new(Cursor::new(bytes)).expect("the reader");
        reader.set_needs(Needs::ALL);
        let blocks = blocks_of(&mut reader).expect("the blocks of the file");
        let chroms = reader.chroms();
        WholePass {
            rows: rows_of(&blocks, chroms),
            chrom_numbers: blocks
                .iter()
                .map(|block| block.chrom.clone().unwrap_or_default())
                .collect(),
            chrom_names: (0..chroms.len())
                .map(|number| {
                    let number = u32::try_from(number).expect("the number of a chromosome");
                    chroms.name(number).unwrap_or("").to_owned()
                })
                .collect(),
        }
    }

    /// The batches of a vars file are decompressed on the threads of the pool
    /// the caller is in, and the blocks come out in the order of the file
    /// whatever that pool is: a pool of one thread and a pool of four give
    /// the same variants, field by field, in the same order, with the same
    /// chromosome numbers and the same table of names.
    ///
    /// The file holds more batches than the reader decompresses at once, so
    /// several windows of batches are read, and the name of the chromosome of
    /// the first variant of each batch first appears in that batch, so a
    /// batch decoded before the one in front of it would number the names in
    /// another order and this would see it. The numbers of a table are handed
    /// out in the order the names are first seen, and nine places of popnei
    /// index their results by a running count of the variants, so a block out
    /// of order is a wrong result and not a slower one.
    ///
    /// The pools are built here and are not rayon's global one, which has one
    /// thread per core of the machine. rayon is a dependency of the targets
    /// that are not wasm, so this test is compiled for those alone.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_number_of_threads_does_not_change_the_blocks_of_a_vars_file() {
        use super::BATCHES_AT_ONCE;

        let num_batches = u64::try_from(BATCHES_AT_ONCE)
            .expect("how many batches are decompressed at once")
            .saturating_mul(2)
            .saturating_add(1);
        let bytes = a_file_of_many_batches(num_batches);
        let in_a_pool = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            pool.install(|| the_whole_pass(&bytes))
        };

        let on_one = in_a_pool(1);

        // Two variants of each batch, and the names of the chromosomes are
        // the one every batch holds and one of each batch after the first.
        assert_eq!(
            on_one.rows.len(),
            usize::try_from(num_batches).expect("the batches") * 2
        );
        assert_eq!(
            on_one.chrom_names.len(),
            usize::try_from(num_batches).expect("the batches")
        );
        assert_eq!(on_one.chrom_names.first().map(String::as_str), Some("chr0"));
        let of_the_last_batch = format!("chr{last}", last = num_batches.saturating_sub(1));
        assert_eq!(
            on_one.chrom_names.last().map(String::as_str),
            Some(of_the_last_batch.as_str()),
            "the name of the chromosome of the last batch is the last of the table"
        );
        assert_eq!(on_one, in_a_pool(4));
        assert_eq!(on_one, in_a_pool(18));
    }

    /// Where each batch of the file in `bytes` is, from its footer, which is
    /// what a test that damages one batch of a file needs.
    fn batches_at(bytes: &[u8]) -> Vec<BatchAt> {
        let file_len = u64::try_from(bytes.len()).expect("the length of the file");
        let mut source = Cursor::new(bytes);
        let footer_bytes = footer_of(&mut source, file_len).expect("the footer");
        let footer = root_as_footer(&footer_bytes).expect("the footer of an arrow file");
        batches_of_the_footer(&footer, file_len).expect("where the batches are")
    }

    /// The bytes of the file with the eight bytes that the message of the
    /// batch at `batch`, counted from 0, starts with written over with
    /// zeroes, which is a batch arrow-rs cannot read.
    fn a_batch_damaged(bytes: &[u8], batch: usize) -> Vec<u8> {
        let at = batches_at(bytes)[batch];
        let offset = usize::try_from(at.offset).expect("where the batch starts");
        let ends = offset
            .checked_add(MESSAGE_START_BYTES)
            .expect("where the message of the batch ends");
        let mut damaged = bytes.to_vec();
        for byte in &mut damaged[offset..ends] {
            *byte = 0;
        }
        damaged
    }

    /// The batches of a file are decompressed several at once, and the error
    /// the reader gives is still the one of the first batch in the order of
    /// the file: a file with a damaged batch gives the error of that batch,
    /// and a file with two damaged batches gives the error of the earlier of
    /// the two and not of whichever thread failed first.
    ///
    /// The blocks before the damaged batch are given, as they are today, and
    /// every call after the error gives none.
    ///
    /// Away from wasm the window is the smaller of [`BATCHES_AT_ONCE`] and
    /// the threads of the pool, so the pool is built with as many threads as
    /// that constant: a machine of one core would otherwise give a window of
    /// one batch, and then the order of the errors of a window is not what is
    /// being read. wasm has no threads and reads one batch at a time.
    #[test]
    fn the_first_damaged_batch_in_the_order_of_the_file_is_the_error_the_reader_gives() {
        #[cfg(not(target_family = "wasm"))]
        rayon::ThreadPoolBuilder::new()
            .num_threads(super::BATCHES_AT_ONCE)
            .build()
            .expect("the pool")
            .install(the_error_of_the_first_damaged_batch);
        #[cfg(target_family = "wasm")]
        the_error_of_the_first_damaged_batch();
    }

    /// The body of the test above, which away from wasm is run in a pool of
    /// as many threads as the window holds batches.
    fn the_error_of_the_first_damaged_batch() {
        let bytes = a_file_of_many_batches(17);
        let the_sixth = a_batch_damaged(&bytes, 5);
        let the_sixth_and_the_eighth = a_batch_damaged(&the_sixth, 7);

        let mut reader = VarsReader::new(Cursor::new(&the_sixth)).expect("the reader");
        let blocks = blocks_of(&mut reader);
        let of_the_sixth = match blocks {
            Err(Error::VarsBatchNotRead { batch, problem }) => {
                assert_eq!(batch, 6);
                problem
            }
            other => panic!("the sixth batch of the file is damaged: {other:?}"),
        };
        // The reader is finished after an error: it gave the five blocks
        // before the damaged batch and gives nothing from here on.
        assert!(
            reader
                .next_block()
                .expect("no block after the error")
                .is_none(),
            "the reader is finished after the error"
        );

        let mut reader =
            VarsReader::new(Cursor::new(&the_sixth_and_the_eighth)).expect("the reader");
        match blocks_of(&mut reader) {
            Err(Error::VarsBatchNotRead { batch, problem }) => {
                assert_eq!(batch, 6, "the earlier of the two damaged batches");
                assert_eq!(problem, of_the_sixth, "the error of the sixth batch alone");
            }
            other => panic!("the sixth and the eighth batches are damaged: {other:?}"),
        }
    }

    /// The five blocks before the damaged batch are given before the error,
    /// which is what says that a reader of a window of batches does not lose
    /// the blocks it decoded before the fault.
    ///
    /// Away from wasm the pool is built with as many threads as the window
    /// holds batches, for the reason the test above gives.
    #[test]
    fn the_blocks_before_a_damaged_batch_are_given_before_the_error() {
        #[cfg(not(target_family = "wasm"))]
        rayon::ThreadPoolBuilder::new()
            .num_threads(super::BATCHES_AT_ONCE)
            .build()
            .expect("the pool")
            .install(the_blocks_before_a_damaged_batch);
        #[cfg(target_family = "wasm")]
        the_blocks_before_a_damaged_batch();
    }

    /// The body of the test above, which away from wasm is run in a pool of
    /// as many threads as the window holds batches.
    fn the_blocks_before_a_damaged_batch() {
        let bytes = a_file_of_many_batches(17);
        let damaged = a_batch_damaged(&bytes, 5);
        let whole = the_whole_pass(&bytes);

        let mut reader = VarsReader::new(Cursor::new(&damaged)).expect("the reader");
        reader.set_needs(Needs::ALL);
        let mut given = Vec::new();
        let error = loop {
            match reader.next_block() {
                Ok(Some(block)) => given.push(block),
                Ok(None) => panic!("the file ended and its sixth batch is damaged"),
                Err(error) => break error,
            }
        };

        assert!(
            matches!(error, Error::VarsBatchNotRead { batch: 6, .. }),
            "the error of the sixth batch: {error:?}"
        );
        assert_eq!(given.len(), 5);
        assert_eq!(
            rows_of(&given, reader.chroms()),
            whole.rows[..10],
            "the ten variants of the five batches before the damaged one"
        );
    }
}
