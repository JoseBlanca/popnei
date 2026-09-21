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
//! `docs/specs/io_vars.md` has the format, the writer and the reader. What
//! is here is the two keys; the writer and the reader are being written.

use serde_json::{Map, Value};

use crate::error::{Error, Result};

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
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "the writer of the vars file is its caller and is being written; the tests of \
                  this module call it already, so this holds for the build without them, and \
                  the lint itself asks for it to go when the writer lands"
    )
)]
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
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "the writer of the vars file is its caller and is being written; the tests of \
                  this module call it already, so this holds for the build without them, and \
                  the lint itself asks for it to go when the writer lands"
    )
)]
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

#[cfg(test)]
mod tests {
    use super::{
        BatchInfo, FORMAT_VERSION, FORMAT_VERSION_READ, Region, VarsMetadata, batches_as_json,
        batches_from_json, metadata_as_json, metadata_from_json,
    };
    use crate::error::Error;
    use crate::variant::Needs;

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
