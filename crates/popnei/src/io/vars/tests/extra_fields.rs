//! The extra fields of a vars file of format 1.3: their columns and the two
//! entries of the `popnei` key that the writer writes, what the reader gives
//! of them and what it refuses, and the bounds of their buffers. "The extra
//! fields" of `docs/specs/io_vars.md` has what they are and "How it is
//! verified" there the cases.

use std::io::Cursor;
use std::sync::Arc;

use arrow_array::cast::AsArray;
use arrow_array::types::{Float32Type, UInt64Type};
use arrow_array::{Array, ArrayRef, FixedSizeListArray, Float32Array, Int64Array, UInt64Array};
use arrow_buffer::NullBuffer;
use arrow_schema::{DataType, Field, Fields};

use super::{
    CASES, CASES_INDIVIDUALS, FileParts, GivenBlocks, ReadRow, Row, block_of, blocks_of,
    blocks_read, cases_block, cases_individuals, file_read, metadata_of_cases, num_vars_of, opened,
    refused, refused_at_the_block, rows_of, swept_with, the_length_of_the_genotypes,
};
use crate::block::{
    Block, BlockReader, ExtraField, ExtraFieldDecl, ExtraKind, ExtraType, ExtraValues, SourceHeader,
};
use crate::error::{Error, Result};
use crate::filters::{FilteringStats, RegionSelection};
use crate::io::vars::{
    BufferHolds, NodeHolds, UNCOMPRESSED_LENGTH_BYTES, VarsMetadata, VarsReader, VarsWriter,
    a_name_of_a_field, extra_type_of, metadata_as_json, metadata_from_json, what_the_buffers_hold,
    what_the_nodes_hold, write_vars,
};
use crate::variant::{ChromTable, Needs};

/// The end of the region of each of the four variants of `cases.vcf`, a
/// variant field of `uint64`.
const END: [u64; 4] = [150, 250, 350, 450];

/// A value for each genotype of the four variants of `cases.vcf`, three
/// individuals each, a genotype field of `float32`. Every value is another,
/// so a value read at the place of another is seen.
const GP: [f32; 12] = [
    0.9, 0.05, 0.5, 0.25, 0.75, 1.0, 0.125, 0.375, 0.625, 0.0, 0.875, 0.3,
];

/// The variant field of the tests.
fn end_decl() -> ExtraFieldDecl {
    ExtraFieldDecl {
        name: "END".to_owned(),
        kind: ExtraKind::Variant,
        ty: ExtraType::U64,
    }
}

/// The genotype field of the tests.
fn gp_decl() -> ExtraFieldDecl {
    ExtraFieldDecl {
        name: "GP".to_owned(),
        kind: ExtraKind::Genotype,
        ty: ExtraType::F32,
    }
}

/// The two fields of the tests, the genotype field first, which is the
/// order the blocks carry them in and not the order of the columns of the
/// file, which hold the variant fields first.
fn gp_and_end() -> Vec<ExtraFieldDecl> {
    vec![gp_decl(), end_decl()]
}

/// The values of `GP` of the variants `first` to `last`, not included.
fn gp_of(first: usize, last: usize) -> ExtraField {
    ExtraField {
        decl: gp_decl(),
        values: ExtraValues::F32(
            GP[first.saturating_mul(CASES_INDIVIDUALS)..last.saturating_mul(CASES_INDIVIDUALS)]
                .to_vec(),
        ),
    }
}

/// The values of `END` of the variants `first` to `last`, not included.
fn end_of(first: usize, last: usize) -> ExtraField {
    ExtraField {
        decl: end_decl(),
        values: ExtraValues::U64(END[first..last].to_vec()),
    }
}

/// The block of the four variants of `cases.vcf` with `GP` and `END`, in
/// that order.
fn cases_block_with_fields(chroms: &mut ChromTable) -> Block {
    let mut block = cases_block(chroms);
    block.extra = Some(vec![gp_of(0, 4), end_of(0, 4)]);
    block
}

/// A reader of the blocks it was given whose header declares `fields`, as
/// the reader of a vars file with extra fields does.
struct WithFields {
    blocks: GivenBlocks,
    header: SourceHeader,
}

impl WithFields {
    /// The reader of those blocks of the individuals of `cases.vcf`, with
    /// the names of `chroms`, whose header declares `fields`.
    fn of(blocks: Vec<Block>, chroms: ChromTable, fields: Vec<ExtraFieldDecl>) -> WithFields {
        WithFields {
            blocks: GivenBlocks::of(blocks, chroms),
            header: SourceHeader {
                individuals: cases_individuals(),
                chrom_lengths: Vec::new(),
                vcf_meta_lines: None,
                keeps_passed: false,
                extra_fields: fields,
            },
        }
    }
}

impl BlockReader for WithFields {
    fn next_block(&mut self) -> Result<Option<Block>> {
        self.blocks.next_block()
    }

    fn individuals(&self) -> &[String] {
        self.blocks.individuals()
    }

    fn ploidy(&self) -> usize {
        self.blocks.ploidy()
    }

    fn chroms(&self) -> &ChromTable {
        self.blocks.chroms()
    }

    fn set_needs(&mut self, needs: Needs) {
        self.blocks.set_needs(needs);
    }

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

/// The vars file of the four variants of `cases.vcf` with `GP` and `END`,
/// in batches of 3 variants, so of 3 and 1.
fn cases_with_fields_written() -> Vec<u8> {
    let mut chroms = ChromTable::new();
    let block = cases_block_with_fields(&mut chroms);
    let reader = WithFields::of(vec![block], chroms, gp_and_end());
    write_vars(reader, Vec::new(), Some(3), false)
        .expect("the file was written")
        .0
}

/// The extra fields of each block.
fn extra_of(blocks: &[Block]) -> Vec<Option<Vec<ExtraField>>> {
    blocks.iter().map(|block| block.extra.clone()).collect()
}

/// The arrow type of the column of a genotype field of `float32` of that
/// many individuals, whose values hold no null.
fn gp_type(num_individuals: i32) -> DataType {
    DataType::FixedSizeList(
        Arc::new(Field::new_list_field(DataType::Float32, false)),
        num_individuals,
    )
}

/// A block with the fields of the tests goes through the writer and the
/// reader as it was: the file has a column for each field after `passed`,
/// the variant field first, the two entries of the `popnei` key list them
/// with their types, and the blocks read back hold their values, cut at the
/// batches, in the order of the file.
#[test]
fn extra_fields_of_a_block_go_through_the_writer_and_the_reader_as_they_were() {
    let bytes = cases_with_fields_written();

    // As another program opens it.
    let file = file_read(bytes.clone());
    let columns: Vec<(String, DataType, bool)> = file.columns[5..].to_vec();
    assert_eq!(
        columns,
        vec![
            ("gts".to_owned(), super::gts_type(6), false),
            ("END".to_owned(), DataType::UInt64, false),
            ("GP".to_owned(), gp_type(3), false),
        ]
    );
    assert_eq!(file.metadata.format_version, "1.3");
    assert_eq!(file.metadata.extra_fields, vec![end_decl(), gp_decl()]);
    let ends: Vec<Vec<u64>> = file
        .rows
        .iter()
        .map(|batch| {
            batch
                .column_by_name("END")
                .expect("the `END` column")
                .as_primitive::<UInt64Type>()
                .values()
                .to_vec()
        })
        .collect();
    assert_eq!(ends, vec![vec![150, 250, 350], vec![450]]);
    let gps: Vec<Vec<f32>> = file
        .rows
        .iter()
        .map(|batch| {
            batch
                .column_by_name("GP")
                .expect("the `GP` column")
                .as_any()
                .downcast_ref::<FixedSizeListArray>()
                .expect("a fixed size list")
                .values()
                .as_primitive::<Float32Type>()
                .values()
                .to_vec()
        })
        .collect();
    assert_eq!(gps, vec![GP[..9].to_vec(), GP[9..].to_vec()]);

    // As popnei reads it.
    let reader = opened(bytes.clone()).expect("the file is a vars file");
    assert_eq!(reader.metadata().extra_fields, vec![end_decl(), gp_decl()]);
    assert_eq!(reader.header().extra_fields, vec![end_decl(), gp_decl()]);
    let (blocks, chroms) = blocks_read(bytes, Needs::ALL);
    assert_eq!(num_vars_of(&blocks), [3, 1]);
    assert_eq!(
        extra_of(&blocks),
        vec![
            Some(vec![end_of(0, 3), gp_of(0, 3)]),
            Some(vec![end_of(3, 4), gp_of(3, 4)]),
        ]
    );
    for block in &blocks {
        assert_eq!(block.fields(), Needs::ALL.difference(Needs::PASSED));
    }
    let expected: Vec<ReadRow> = CASES.iter().map(super::row_read).collect();
    assert_eq!(rows_of(&blocks, &chroms), expected);
}

/// The two entries of the `popnei` key are written as the example of the
/// spec, the name and the type of each field, and parsed back; a file with
/// no field has both, as `[]`.
#[test]
fn extra_fields_are_written_in_the_popnei_key_and_parsed_back() {
    let metadata = VarsMetadata {
        extra_fields: vec![
            end_decl(),
            gp_decl(),
            ExtraFieldDecl {
                name: "NODE_A".to_owned(),
                kind: ExtraKind::Genotype,
                ty: ExtraType::I32,
            },
        ],
        ..metadata_of_cases()
    };
    let written = metadata_as_json(&metadata);
    assert!(
        written.ends_with(
            r#""chrom_lengths":[],"var_fields":[["END","uint64"]],"gt_fields":[["GP","float32"],["NODE_A","int32"]]}"#
        ),
        "{written}"
    );
    assert_eq!(metadata_from_json(&written).expect("parsed back"), metadata);
}

/// A key of 1.2 or before has neither entry and no extra field.
#[test]
fn extra_fields_of_a_key_of_1_2_are_none() {
    let key = r#"{"format_version":"1.2","individuals":["ind1"],"ploidy":2,"num_vars_per_block":3,"chrom_lengths":[]}"#;
    let metadata = metadata_from_json(key).expect("the key of 1.2");
    assert_eq!(metadata.extra_fields, Vec::new());
}

/// A field of each of the five types is read back as it was written,
/// values at the ends of their ranges among them: a reader that read one
/// type as another, or a value through a narrower one, gives other values.
#[test]
fn extra_fields_of_each_of_the_five_types_are_read_back_as_they_were_written() {
    let decl = |name: &str, kind: ExtraKind, ty: ExtraType| ExtraFieldDecl {
        name: name.to_owned(),
        kind,
        ty,
    };
    let fields = vec![
        ExtraField {
            decl: decl("NODE_A", ExtraKind::Genotype, ExtraType::I32),
            values: ExtraValues::I32(vec![i32::MIN, -1, 0, 1, 2, i32::MAX]),
        },
        ExtraField {
            decl: decl("LEN", ExtraKind::Variant, ExtraType::I64),
            values: ExtraValues::I64(vec![-5_000_000_000, i64::MAX]),
        },
        ExtraField {
            decl: decl("_end.2", ExtraKind::Variant, ExtraType::U64),
            values: ExtraValues::U64(vec![u64::MAX, 5_000_000_000]),
        },
        ExtraField {
            decl: decl("GP", ExtraKind::Genotype, ExtraType::F32),
            values: ExtraValues::F32(vec![f32::MIN_POSITIVE, -0.5, 1e30, 0.0, 3.25, -7.0]),
        },
        ExtraField {
            decl: decl("DIV", ExtraKind::Genotype, ExtraType::F64),
            values: ExtraValues::F64(vec![1e-300, -2.5, 0.1, 1e300, f64::MAX, -0.0]),
        },
    ];
    let declared: Vec<ExtraFieldDecl> = fields.iter().map(|field| field.decl.clone()).collect();
    let mut chroms = ChromTable::new();
    let rows: Vec<&Row> = CASES[..2].iter().collect();
    let mut block = block_of(&rows, &mut chroms);
    block.extra = Some(fields.clone());
    let reader = WithFields::of(vec![block], chroms, declared);
    let (bytes, _) = write_vars(reader, Vec::new(), None, false).expect("the file was written");

    let (blocks, _) = blocks_read(bytes, Needs::ALL);

    // The variant fields first, each kind in the order it was declared.
    let in_the_file: Vec<ExtraField> = [1, 2, 0, 3, 4]
        .iter()
        .map(|at| fields[*at].clone())
        .collect();
    assert_eq!(extra_of(&blocks), vec![Some(in_the_file)]);
}

/// The floats of a field go through the writer and the reader bit for bit:
/// a NaN, a NaN with a payload, -0.0 and the two infinities, which `==`
/// cannot tell from 0.0 or from another NaN, are compared with `to_bits`.
#[test]
fn extra_field_floats_go_through_the_writer_and_the_reader_bit_for_bit() {
    let f32s = vec![
        f32::NAN,
        f32::from_bits(0x7fc0_1234),
        -0.0,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::from_bits(0xffc0_0001),
    ];
    let f64s = vec![
        f64::NAN,
        f64::from_bits(0x7ff8_0000_dead_beef),
        -0.0,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::from_bits(0xfff8_0000_0000_0001),
    ];
    let decl = |name: &str, ty: ExtraType| ExtraFieldDecl {
        name: name.to_owned(),
        kind: ExtraKind::Genotype,
        ty,
    };
    let fields = vec![
        ExtraField {
            decl: decl("GP", ExtraType::F32),
            values: ExtraValues::F32(f32s.clone()),
        },
        ExtraField {
            decl: decl("DIV", ExtraType::F64),
            values: ExtraValues::F64(f64s.clone()),
        },
    ];
    let declared: Vec<ExtraFieldDecl> = fields.iter().map(|field| field.decl.clone()).collect();
    let mut chroms = ChromTable::new();
    let rows: Vec<&Row> = CASES[..2].iter().collect();
    let mut block = block_of(&rows, &mut chroms);
    block.extra = Some(fields);
    let reader = WithFields::of(vec![block], chroms, declared);
    let (bytes, _) = write_vars(reader, Vec::new(), None, false).expect("the file was written");

    let (blocks, _) = blocks_read(bytes, Needs::ALL);

    let [block] = blocks.as_slice() else {
        panic!("one block, and {} were read", blocks.len());
    };
    let Some([gp, div]) = block.extra.as_deref() else {
        panic!("the two fields, and {:?} were read", block.extra);
    };
    let (ExtraValues::F32(gp), ExtraValues::F64(div)) = (&gp.values, &div.values) else {
        panic!("a field of float32 and one of float64, and {gp:?} {div:?} were read");
    };
    let bits_32 = |values: &[f32]| {
        values
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    };
    let bits_64 = |values: &[f64]| {
        values
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    };
    assert_eq!(bits_32(gp), bits_32(&f32s));
    assert_eq!(bits_64(div), bits_64(&f64s));
}

/// The schema is written before the first batch, so a source of no
/// variant gives a file with the columns of its fields and both entries,
/// and its reader says which fields it has.
#[test]
fn extra_field_columns_are_in_a_file_with_no_variant() {
    let reader = WithFields::of(Vec::new(), ChromTable::new(), gp_and_end());
    let (bytes, num_vars) =
        write_vars(reader, Vec::new(), Some(3), false).expect("the file was written");
    assert_eq!(num_vars, 0);

    let file = file_read(bytes.clone());
    assert_eq!(
        file.columns,
        vec![
            ("gts".to_owned(), super::gts_type(6), false),
            ("END".to_owned(), DataType::UInt64, false),
            ("GP".to_owned(), gp_type(3), false),
        ]
    );
    let reader = opened(bytes.clone()).expect("the file is a vars file");
    assert_eq!(reader.header().extra_fields, vec![end_decl(), gp_decl()]);
    let (blocks, _) = blocks_read(bytes, Needs::ALL);
    assert!(blocks.is_empty());
}

/// `write_vars` writes the fields of the header of its source, which for a
/// vars file are those of its key, so a file written from another keeps
/// them, here cut at other batches.
#[test]
fn extra_fields_of_a_vars_file_are_kept_in_a_file_written_from_it() {
    let first = cases_with_fields_written();
    let reader = opened(first).expect("the first file");

    let (second, num_vars) =
        write_vars(reader, Vec::new(), Some(2), false).expect("the second file");

    assert_eq!(num_vars, 4);
    let (blocks, _) = blocks_read(second, Needs::ALL);
    assert_eq!(num_vars_of(&blocks), [2, 2]);
    assert_eq!(
        extra_of(&blocks),
        vec![
            Some(vec![end_of(0, 2), gp_of(0, 2)]),
            Some(vec![end_of(2, 4), gp_of(2, 4)]),
        ]
    );
}

/// The blocks of a file with no extra field, and of one read without
/// [`Needs::EXTRA`], have none.
#[test]
fn extra_fields_are_none_in_a_file_of_none_and_when_they_are_not_asked_for() {
    let (blocks, _) = blocks_read(super::cases_written_in_batches_of(3), Needs::ALL);
    assert_eq!(extra_of(&blocks), vec![None, None]);

    let without = Needs::ALL.difference(Needs::EXTRA);
    let (blocks, _) = blocks_read(cases_with_fields_written(), without);
    assert_eq!(extra_of(&blocks), vec![None, None]);
    for block in &blocks {
        assert_eq!(block.fields(), without.difference(Needs::PASSED));
    }
}

/// The blocks of the file read with `select_extra_fields(names)`, every
/// field asked for.
fn read_selecting(names: Option<Vec<&str>>) -> Vec<Option<Vec<ExtraField>>> {
    let mut reader = opened(cases_with_fields_written()).expect("the file");
    reader
        .select_extra_fields(names.map(|names| names.into_iter().map(str::to_owned).collect()))
        .expect("the fields of the file");
    extra_of(&blocks_of(&mut reader).expect("the blocks"))
}

/// `select_extra_fields` gives the fields named, in the order of the file
/// whatever the order of the names; `None` gives every one, and a
/// selection of none a block with no extra field.
#[test]
fn select_extra_fields_gives_the_fields_named_in_the_order_of_the_file() {
    assert_eq!(
        read_selecting(Some(vec!["GP"])),
        vec![Some(vec![gp_of(0, 3)]), Some(vec![gp_of(3, 4)])]
    );
    assert_eq!(
        read_selecting(Some(vec!["END"])),
        vec![Some(vec![end_of(0, 3)]), Some(vec![end_of(3, 4)])]
    );
    let both = vec![
        Some(vec![end_of(0, 3), gp_of(0, 3)]),
        Some(vec![end_of(3, 4), gp_of(3, 4)]),
    ];
    assert_eq!(read_selecting(Some(vec!["GP", "END"])), both);
    assert_eq!(read_selecting(None), both);
    assert_eq!(read_selecting(Some(Vec::new())), vec![None, None]);
}

/// `select_extra_fields` makes the header declare the fields selected, as
/// the blocks carry them, so `write_vars` of a reader of `GP` alone writes
/// a file of `GP` alone; the metadata keeps every field of the file.
#[test]
fn select_extra_fields_sets_the_fields_of_the_header_and_write_vars_writes_them() {
    let mut reader = opened(cases_with_fields_written()).expect("the file");
    reader
        .select_extra_fields(Some(vec!["GP".to_owned()]))
        .expect("a field of the file");

    assert_eq!(reader.header().extra_fields, vec![gp_decl()]);
    assert_eq!(reader.metadata().extra_fields, vec![end_decl(), gp_decl()]);
    let (bytes, num_vars) = write_vars(reader, Vec::new(), None, false).expect("the file of `GP`");
    assert_eq!(num_vars, 4);
    let reader = opened(bytes.clone()).expect("the file of `GP`");
    assert_eq!(reader.header().extra_fields, vec![gp_decl()]);
    let (blocks, _) = blocks_read(bytes, Needs::ALL);
    assert_eq!(extra_of(&blocks), vec![Some(vec![gp_of(0, 4)])]);

    let mut reader = opened(cases_with_fields_written()).expect("the file");
    reader
        .select_extra_fields(Some(Vec::new()))
        .expect("no field");
    assert_eq!(reader.header().extra_fields, Vec::new());
    reader.select_extra_fields(None).expect("every field");
    assert_eq!(reader.header().extra_fields, vec![end_decl(), gp_decl()]);
}

/// A selection holds from the next block, as `set_needs` does.
#[test]
fn select_extra_fields_holds_from_the_next_block() {
    let mut reader = opened(cases_with_fields_written()).expect("the file");
    let first = reader
        .next_block()
        .expect("a block")
        .expect("the first block");
    reader
        .select_extra_fields(Some(vec!["GP".to_owned()]))
        .expect("a field of the file");
    let second = reader
        .next_block()
        .expect("a block")
        .expect("the second block");

    assert_eq!(first.extra, Some(vec![end_of(0, 3), gp_of(0, 3)]));
    assert_eq!(second.extra, Some(vec![gp_of(3, 4)]));
}

/// A name the file has no field of is an error that names the fields it
/// has, and the selection stays as it was.
#[test]
fn select_extra_fields_of_a_name_the_file_lacks_names_the_fields_it_has() {
    let mut reader = opened(cases_with_fields_written()).expect("the file");

    let error = reader
        .select_extra_fields(Some(vec!["GP".to_owned(), "DP".to_owned()]))
        .expect_err("the file has no `DP`");

    let Error::VarsExtraFieldNotInTheFile { name, fields } = &error else {
        panic!("the name `DP` gave {error}");
    };
    assert_eq!(name, "DP");
    assert_eq!(fields, &["END".to_owned(), "GP".to_owned()]);
    let message = error.to_string();
    assert!(message.contains("`DP`"), "{message}");
    assert!(message.contains("`END`, `GP`"), "{message}");
    assert_eq!(
        extra_of(&blocks_of(&mut reader).expect("the blocks")),
        vec![
            Some(vec![end_of(0, 3), gp_of(0, 3)]),
            Some(vec![end_of(3, 4), gp_of(3, 4)]),
        ]
    );

    let mut reader = opened(super::cases_written_in_batches_of(3)).expect("the file");
    let error = reader
        .select_extra_fields(Some(vec!["GP".to_owned()]))
        .expect_err("the file has no field");
    assert!(error.to_string().contains("it has none"), "{error}");
}

/// The writer refuses a name of a field that no field can have, one of a
/// column every vars file can have and one given twice, whatever the kind
/// of the two.
#[test]
fn set_extra_fields_refuses_a_name_no_extra_field_can_have() {
    for (fields, name, problem) in [
        (
            vec![ExtraFieldDecl {
                name: "1GP".to_owned(),
                ..gp_decl()
            }],
            "1GP",
            "no field can have",
        ),
        (
            vec![ExtraFieldDecl {
                name: "qual".to_owned(),
                ..end_decl()
            }],
            "qual",
            "every vars file can have",
        ),
        (
            vec![
                end_decl(),
                ExtraFieldDecl {
                    name: "END".to_owned(),
                    ..gp_decl()
                },
            ],
            "END",
            "twice",
        ),
    ] {
        let mut writer =
            VarsWriter::new(Vec::new(), &cases_individuals(), 2, 3).expect("the writer");

        let error = writer
            .set_extra_fields(fields)
            .expect_err("the name is refused");

        let Error::VarsExtraField {
            name: found,
            problem: said,
        } = &error
        else {
            panic!("the field `{name}` gave {error}");
        };
        assert_eq!(found, name);
        assert!(said.contains(problem), "{said}");
    }
}

/// Every one of the five types is parsed back from the name it is written
/// with, and the list holds five different ones.
#[test]
fn every_extra_type_is_parsed_back_from_its_name() {
    let names: Vec<String> = ExtraType::ALL.iter().map(ToString::to_string).collect();
    assert_eq!(names, ["int32", "int64", "uint64", "float32", "float64"]);
    for ty in ExtraType::ALL {
        assert_eq!(extra_type_of(&ty.to_string()), Some(ty));
    }
    assert_eq!(extra_type_of("int8"), None);
}

/// The rule of the VCF for the key of an INFO or a FORMAT field.
#[test]
fn a_name_of_an_extra_field_follows_the_rule_of_the_keys_of_a_vcf() {
    for name in ["END", "_", "_end.2", "A1.b_c", "z"] {
        assert!(a_name_of_a_field(name), "{name}");
    }
    for name in ["", "1A", ".A", "A-B", "A B", "Á", "END\n"] {
        assert!(!a_name_of_a_field(name), "{name:?}");
    }
}

/// A block whose extra fields are not the ones declared, in their order,
/// is refused and nothing of it is written: the batches of the file all
/// hold the columns the `popnei` key declares.
#[test]
fn a_block_of_other_extra_fields_than_declared_is_refused_and_nothing_of_it_is_written() {
    let mut chroms = ChromTable::new();
    let right = cases_block_with_fields(&mut chroms);
    let mut writer = VarsWriter::new(Vec::new(), &cases_individuals(), 2, 4).expect("the writer");
    writer.set_extra_fields(gp_and_end()).expect("the fields");

    for (extra, found) in [
        (None, "none"),
        (Some(vec![gp_of(0, 4)]), "`GP` (genotype, float32)"),
        (
            Some(vec![end_of(0, 4), gp_of(0, 4)]),
            "`END` (variant, uint64), `GP` (genotype, float32)",
        ),
    ] {
        let mut wrong = cases_block(&mut chroms);
        wrong.extra = extra;

        let error = writer
            .write_block(wrong, &chroms)
            .expect_err("the block is refused");

        let Error::VarsBlockExtraFields {
            declared,
            found: said,
        } = &error
        else {
            panic!("the block of the fields {found} gave {error}");
        };
        assert_eq!(
            declared,
            "`GP` (genotype, float32), `END` (variant, uint64)"
        );
        assert_eq!(said, found);
    }

    writer.write_block(right, &chroms).expect("the right block");
    let bytes = writer.finish().expect("the file");
    let (blocks, _) = blocks_read(bytes, Needs::ALL);
    assert_eq!(num_vars_of(&blocks), [4]);
    assert_eq!(
        extra_of(&blocks),
        vec![Some(vec![end_of(0, 4), gp_of(0, 4)])]
    );
}

/// The fields are in the schema, which the first block writes, so other
/// fields declared after it are refused, and the same ones are not.
#[test]
fn set_extra_fields_after_the_first_block_is_refused_unless_they_are_the_same() {
    let mut chroms = ChromTable::new();
    let block = cases_block_with_fields(&mut chroms);
    let mut writer = VarsWriter::new(Vec::new(), &cases_individuals(), 2, 4).expect("the writer");
    writer.set_extra_fields(gp_and_end()).expect("the fields");
    writer.write_block(block, &chroms).expect("the block");

    writer
        .set_extra_fields(gp_and_end())
        .expect("the same fields");
    let error = writer
        .set_extra_fields(vec![end_decl()])
        .expect_err("other fields after the first block");

    assert!(
        matches!(&error, Error::VarsBlockExtraFields { found, .. } if found == "`END` (variant, uint64)"),
        "{error}"
    );
}

/// The four variants of `cases.vcf` as columns built with arrow-rs, with
/// an `END` and a `GP` column after the others and a `popnei` key that
/// declares `fields`.
fn parts_with(fields: &str, end: ArrayRef, gp: (Field, ArrayRef)) -> FileParts {
    let mut parts = FileParts::of_cases();
    parts
        .columns
        .push((Field::new("END", end.data_type().clone(), true), end));
    parts.columns.push(gp);
    let key = metadata_as_json(&metadata_of_cases());
    let key = key.replace(r#","var_fields":[],"gt_fields":[]"#, fields);
    parts.popnei = Some(key);
    parts
}

/// The `END` column of the four variants, as a column of `uint64`.
fn end_column() -> ArrayRef {
    Arc::new(UInt64Array::from(END.to_vec()))
}

/// The `GP` column of the four variants, a fixed size list of `width`
/// values of `float32` for each variant, whose values take nulls.
fn gp_column(values: Vec<Option<f32>>, width: i32) -> (Field, ArrayRef) {
    let item = Arc::new(Field::new_list_field(DataType::Float32, true));
    let values: ArrayRef = Arc::new(Float32Array::from(values));
    let column = FixedSizeListArray::try_new(Arc::clone(&item), width, values, None)
        .expect("the list of the values");
    (
        Field::new("GP", DataType::FixedSizeList(item, width), true),
        Arc::new(column),
    )
}

/// The entries of the `popnei` key of the two fields of the tests.
const BOTH: &str = r#","var_fields":[["END","uint64"]],"gt_fields":[["GP","float32"]]"#;

/// A file of another program, whose columns take nulls and have none, is
/// read: what the tests of the errors change is what makes them fail.
#[test]
fn extra_fields_of_a_file_of_another_program_are_read() {
    let gp: Vec<Option<f32>> = GP.iter().copied().map(Some).collect();
    let parts = parts_with(BOTH, end_column(), gp_column(gp, 3));

    let (blocks, _) = blocks_read(parts.written(), Needs::ALL);

    assert_eq!(
        extra_of(&blocks),
        vec![Some(vec![end_of(0, 4), gp_of(0, 4)])]
    );
}

/// The field of an error of the file as a whole, and what is wrong with
/// it.
fn field_refused(bytes: Vec<u8>) -> (String, String) {
    let error = refused(bytes);
    let Error::VarsExtraField { name, problem } = error else {
        panic!("the file gave {error}");
    };
    (name, problem)
}

/// An entry without its column, a column of another type, and a genotype
/// field of another width than the individuals are errors of the file as a
/// whole, with the type or the width found and the one expected.
#[test]
fn an_extra_field_that_is_not_what_its_entry_says_is_refused_with_both_types() {
    let gp = || GP.iter().copied().map(Some).collect::<Vec<_>>();

    let mut parts = parts_with(BOTH, end_column(), gp_column(gp(), 3));
    parts.columns.retain(|(field, _)| field.name() != "END");
    let (name, problem) = field_refused(parts.written());
    assert_eq!(name, "END");
    assert!(problem.contains("has no column"), "{problem}");
    assert!(problem.contains("UInt64"), "{problem}");

    let ends: ArrayRef = Arc::new(Int64Array::from(vec![150, 250, 350, 450]));
    let (name, problem) = field_refused(parts_with(BOTH, ends, gp_column(gp(), 3)).written());
    assert_eq!(name, "END");
    assert!(problem.contains("arrow type Int64"), "{problem}");
    assert!(problem.contains("says UInt64"), "{problem}");

    let mut four: Vec<Option<f32>> = gp();
    four.extend([Some(0.0); 4]);
    let (name, problem) =
        field_refused(parts_with(BOTH, end_column(), gp_column(four, 4)).written());
    assert_eq!(name, "GP");
    assert!(problem.contains("FixedSizeList<Float32>[4]"), "{problem}");
    assert!(
        problem.contains("says FixedSizeList<Float32>[3]"),
        "{problem}"
    );

    let plain: ArrayRef = Arc::new(Float32Array::from(vec![0.5; 4]));
    let parts = parts_with(
        BOTH,
        end_column(),
        (Field::new("GP", DataType::Float32, false), plain),
    );
    let (name, problem) = field_refused(parts.written());
    assert_eq!(name, "GP");
    assert!(problem.contains("arrow type Float32,"), "{problem}");
}

/// A name the rule refuses, one of a column of "What it holds" and one in
/// the lists twice are errors of the file as a whole.
#[test]
fn an_extra_field_name_no_field_can_have_is_refused_when_the_file_is_opened() {
    let gp = || gp_column(GP.iter().copied().map(Some).collect(), 3);
    for (fields, name, problem) in [
        (
            r#","var_fields":[["1END","uint64"]],"gt_fields":[["GP","float32"]]"#,
            "1END",
            "no field can have",
        ),
        (
            r#","var_fields":[["END","uint64"],["pos","uint64"]],"gt_fields":[["GP","float32"]]"#,
            "pos",
            "every vars file can have",
        ),
        (
            r#","var_fields":[["END","uint64"]],"gt_fields":[["GP","float32"],["END","float32"]]"#,
            "END",
            "twice",
        ),
    ] {
        let (found, said) = field_refused(parts_with(fields, end_column(), gp()).written());
        assert_eq!(found, name);
        assert!(said.contains(problem), "{said}");
    }
}

/// An entry that is not a list of pairs of a name and one of the five types
/// is not a vars file, as a `chrom_lengths` that is not one.
#[test]
fn an_extra_field_entry_that_is_not_pairs_of_a_name_and_a_type_is_not_a_vars_file() {
    for fields in [
        r#","var_fields":"END","gt_fields":[]"#,
        r#","var_fields":[["END"]],"gt_fields":[]"#,
        r#","var_fields":[["END","uint64","x"]],"gt_fields":[]"#,
        r#","var_fields":[["END","int8"]],"gt_fields":[]"#,
        r#","var_fields":[],"gt_fields":[[1,"float32"]]"#,
    ] {
        let parts = parts_with(
            fields,
            end_column(),
            gp_column(GP.iter().copied().map(Some).collect(), 3),
        );

        let error = refused(parts.written());

        let Error::NotAVarsFile { problem } = &error else {
            panic!("the entries {fields} gave {error}");
        };
        assert!(problem.contains("_fields"), "{problem}");
    }
}

/// A null in the column of a variant field, or among the values of a
/// genotype field, is the error of a null, with the field and the variant.
#[test]
fn an_extra_field_with_a_null_names_the_field_and_the_variant() {
    let gp = || gp_column(GP.iter().copied().map(Some).collect(), 3);
    let ends: ArrayRef = Arc::new(UInt64Array::from(vec![
        Some(150),
        Some(250),
        None,
        Some(450),
    ]));
    let error = refused_at_the_block(parts_with(BOTH, ends, gp()).written());
    let Error::VarsExtraFieldNull { name, var } = &error else {
        panic!("the null of `END` gave {error}");
    };
    assert_eq!((name.as_str(), *var), ("END", 3));
    let message = error.to_string();
    assert!(message.contains("`END` column"), "{message}");
    assert!(message.contains("variant 3"), "{message}");

    // The eighth value is of the third variant.
    let mut values: Vec<Option<f32>> = GP.iter().copied().map(Some).collect();
    values[7] = None;
    let error =
        refused_at_the_block(parts_with(BOTH, end_column(), gp_column(values, 3)).written());
    let Error::VarsExtraFieldNull { name, var } = &error else {
        panic!("the null of `GP` gave {error}");
    };
    assert_eq!((name.as_str(), *var), ("GP", 3));

    // A null of the list itself, of the second variant.
    let (field, column) = gp();
    let list = column
        .as_any()
        .downcast_ref::<FixedSizeListArray>()
        .expect("the list")
        .clone();
    let (item, width, values, _) = list.into_parts();
    let with_a_null = FixedSizeListArray::try_new(
        item,
        width,
        values,
        Some(NullBuffer::from(vec![true, false, true, true])),
    )
    .expect("the list with a null");
    let error = refused_at_the_block(
        parts_with(BOTH, end_column(), (field, Arc::new(with_a_null))).written(),
    );
    let Error::VarsExtraFieldNull { name, var } = &error else {
        panic!("the null of `GP` gave {error}");
    };
    assert_eq!((name.as_str(), *var), ("GP", 2));
}

/// The buffers of the extra columns are bounded by their columns, which
/// the schema and the rows give exactly, and the walk goes on after them:
/// a variant field is the rows times the bytes of its type, a genotype
/// field that times the individuals, and each has the bits of its nulls
/// before it.
#[test]
fn extra_field_columns_are_walked_with_the_exact_bounds_of_their_buffers() {
    let fields = Fields::from(vec![
        Field::new("END", DataType::UInt64, false),
        Field::new("GP", gp_type(3), false),
        Field::new("NODE_A", DataType::Int32, false),
    ]);

    assert_eq!(
        what_the_buffers_hold(&fields, 4),
        vec![
            // `END`: the nulls of 4 rows, and 4 values of 8 bytes.
            BufferHolds::Bytes(1),
            BufferHolds::Bytes(32),
            // `GP`: the nulls of the lists, then those of the 12 values
            // and 12 values of 4 bytes.
            BufferHolds::Bytes(1),
            BufferHolds::Bytes(2),
            BufferHolds::Bytes(48),
            // `NODE_A`, after them.
            BufferHolds::Bytes(1),
            BufferHolds::Bytes(16),
        ]
    );
    assert_eq!(
        what_the_nodes_hold(&fields, 4),
        vec![
            NodeHolds::Values(4),
            NodeHolds::Values(4),
            NodeHolds::Values(12),
            NodeHolds::Values(4),
        ]
    );
}

/// The vars file of a thousand copies of the first variant of `cases.vcf`
/// with `GP` and `END`, whose values repeat, so lz4 compresses the buffer
/// of `GP` and says how long it is decompressed.
fn a_thousand_variants_with_fields() -> Vec<u8> {
    let mut chroms = ChromTable::new();
    let rows: Vec<&Row> = vec![&CASES[0]; 1000];
    let mut block = block_of(&rows, &mut chroms);
    block.extra = Some(vec![
        ExtraField {
            decl: gp_decl(),
            values: ExtraValues::F32(vec![0.5; 3000]),
        },
        ExtraField {
            decl: end_decl(),
            values: ExtraValues::U64(vec![7; 1000]),
        },
    ]);
    let reader = WithFields::of(vec![block], chroms, gp_and_end());
    write_vars(reader, Vec::new(), Some(1000), false)
        .expect("the file was written")
        .0
}

/// A buffer of an extra column that says it holds more than its column is
/// refused before arrow-rs asks the machine for the memory it says, by
/// the bound of the column and not by what lz4 gives, which for the values
/// of `GP` would be some 90 times more.
#[test]
fn an_extra_field_buffer_that_says_more_than_its_column_is_refused_before_arrow_reads_it() {
    let whole = a_thousand_variants_with_fields();
    // The last buffer of the batch is that of the values of `GP`, the last
    // column of the file.
    let at = the_length_of_the_genotypes(&whole, 0);
    let length = at..at.saturating_add(UNCOMPRESSED_LENGTH_BYTES);
    let says = i64::from_le_bytes(whole[length.clone()].try_into().expect("eight bytes"));
    // A thousand variants of three individuals, four bytes each.
    assert_eq!(says, 12000);

    let mut bytes = whole.clone();
    bytes[length].copy_from_slice(&12001_i64.to_le_bytes());
    let error = refused_at_the_block(bytes);

    let Error::VarsBatchNotRead { batch, problem } = &error else {
        panic!("the buffer that says 12001 bytes gave {error}");
    };
    assert_eq!(*batch, 1);
    assert!(problem.contains("12001"), "{problem}");
    assert!(problem.contains("its column holds 12000"), "{problem}");

    let (blocks, _) = blocks_read(whole, Needs::ALL);
    assert_eq!(num_vars_of(&blocks), [1000]);
}

/// How many panics the net of the reader caught in the sweep over four
/// values of each byte of the file with extra fields, on 9 November 2026
/// with arrow-rs 60: eight about a buffer whose length is not a whole
/// number of the values of its column and three about a buffer of bits
/// that is not long enough, asserts of arrow-rs that only the walk of the
/// schema its decoder does would see. The sweep gave 20549 files in 0.22 s:
/// 8753 errors, 11274 read as the whole file and 522 read as another file
/// with no error.
const PANICS_CAUGHT_WITH_FIELDS: u64 = 11;

/// The four variants of `cases.vcf` with a `float32` genotype field and an
/// `int64` variant field, in batches of 3 and 1.
fn cases_with_a_float32_and_an_int64_field() -> Vec<u8> {
    let len = ExtraFieldDecl {
        name: "LEN".to_owned(),
        kind: ExtraKind::Variant,
        ty: ExtraType::I64,
    };
    let mut chroms = ChromTable::new();
    let mut block = cases_block(&mut chroms);
    block.extra = Some(vec![
        gp_of(0, 4),
        ExtraField {
            decl: len.clone(),
            values: ExtraValues::I64(vec![-3, 40, 5_000_000_000, 0]),
        },
    ]);
    let reader = WithFields::of(vec![block], chroms, vec![gp_decl(), len]);
    write_vars(reader, Vec::new(), Some(3), false)
        .expect("the file was written")
        .0
}

/// The variants of a file and the extra fields of each of its blocks.
type RowsAndFields = (Vec<ReadRow>, Vec<Option<Vec<ExtraField>>>);

/// The variants of the file in those bytes and the extra fields of each
/// block, with every field, or the error it gave.
fn rows_and_fields_or_error(bytes: Vec<u8>) -> Result<RowsAndFields> {
    let mut reader = VarsReader::new(Cursor::new(bytes))?;
    reader.set_needs(Needs::ALL);
    let blocks = blocks_of(&mut reader)?;
    for block in &blocks {
        block.check()?;
    }
    Ok((rows_of(&blocks, reader.chroms()), extra_of(&blocks)))
}

/// No change of one byte of a vars file with extra fields reaches a panic
/// or an abort: the sweep of the reader, four values of each byte, over a
/// file whose columns after `gts` are walked too.
///
/// Away from wasm it runs in a pool of one thread. The reader decodes the
/// batches of a window on the threads of the pool it is called in, and
/// `swept_with` counts the panics of the thread it runs on, so outside
/// such a pool the panics the net caught happen on other threads and none
/// is counted.
#[test]
fn no_change_of_one_byte_of_a_vars_file_with_extra_fields_reaches_a_panic() {
    #[cfg(not(target_family = "wasm"))]
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .expect("a pool of one thread")
        .install(the_sweep_with_extra_fields);
    #[cfg(target_family = "wasm")]
    the_sweep_with_extra_fields();
}

/// The body of the test above, inside its pool of one thread.
fn the_sweep_with_extra_fields() {
    let whole = cases_with_a_float32_and_an_int64_field();
    let expected = rows_and_fields_or_error(whole.clone()).expect("the whole file is read");
    assert_eq!(expected.0.len(), 4);
    assert_eq!(expected.1.len(), 2);

    let mut sweep = swept_with(&whole, &[0x00, 0xff], &expected, rows_and_fields_or_error);
    for bit in [0x01_u8, 0x80] {
        for (at, byte) in whole.iter().enumerate() {
            let mut bytes = whole.clone();
            bytes[at] = byte ^ bit;
            match rows_and_fields_or_error(bytes) {
                Err(_) => sweep.errors = sweep.errors.saturating_add(1),
                Ok(read) if read == expected => sweep.the_same = sweep.the_same.saturating_add(1),
                Ok(_) => sweep.others = sweep.others.saturating_add(1),
            }
        }
    }

    assert!(sweep.errors > 100, "{} errors", sweep.errors);
    assert!(sweep.the_same > 100, "{} files read whole", sweep.the_same);
    sweep.panics_are_of_arrow();
    assert!(
        sweep.caught() <= PANICS_CAUGHT_WITH_FIELDS,
        "{} panics were caught, and {PANICS_CAUGHT_WITH_FIELDS} were counted on 9 November 2026: {:?}",
        sweep.caught(),
        sweep.panics
    );
    println!("{sweep:?}");
}
