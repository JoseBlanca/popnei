//! The order of the variants that `write_vars` checks, "Sorted variants" of
//! "The writer of segments" of `docs/specs/io_vars.md`: the positions of a
//! chromosome do not go down, its variants come together, two variants at
//! one position are in order, and `allow_unsorted` writes them as they
//! come.

use super::{GivenBlocks, ReadRow, Row, block_of, blocks_read, row_read, rows_of};
use crate::block::Block;
use crate::error::{Error, Result};
use crate::io::vars::write_vars;
use crate::variant::{ChromTable, Needs};

/// The genotypes of every variant of these tests, which are of the three
/// individuals of the ploidy 2 of `GivenBlocks`.
const GTS: &[i8] = &[0, 0, 0, 1, 1, 1];

/// A variant at `pos` of `chrom`.
const fn at(chrom: &'static str, pos: u64) -> Row {
    Row {
        chrom,
        pos,
        id: "",
        alleles: &["A", "T"],
        qual: None,
        gts: GTS,
    }
}

/// A position that goes down on its chromosome, the first case of the spec.
const DOWN: [Row; 2] = [at("chr1", 5000), at("chr1", 100)];

/// A chromosome that comes back after another, the second case of the spec,
/// after a first variant so that the one out of order is the fourth of the
/// file and the first of the second block `write_vars` cuts.
const COMES_BACK: [Row; 4] = [
    at("chr1", 100),
    at("chr1", 200),
    at("chr2", 300),
    at("chr1", 150),
];

/// Variants in order: two at one position, and chromosomes that are not in
/// the order of their names.
const SORTED: [Row; 4] = [
    at("chr2", 5),
    at("chr2", 5),
    at("chr1", 100),
    at("chr1", 100),
];

/// What `write_vars` gives for `rows`, given in blocks of `given_in` and
/// written in batches of `num_vars_per_block`.
fn written(
    rows: &[Row],
    given_in: usize,
    num_vars_per_block: usize,
    allow_unsorted: bool,
) -> Result<Vec<u8>> {
    let mut chroms = ChromTable::new();
    let blocks: Vec<Block> = rows
        .chunks(given_in)
        .map(|chunk| {
            let chunk: Vec<&Row> = chunk.iter().collect();
            block_of(&chunk, &mut chroms)
        })
        .collect();
    let reader = GivenBlocks::of(blocks, chroms);
    write_vars(reader, Vec::new(), Some(num_vars_per_block), allow_unsorted).map(|(bytes, _)| bytes)
}

/// The variants of the file in `bytes`, as a block gives them.
fn rows_read(bytes: Vec<u8>) -> Vec<ReadRow> {
    let (blocks, chroms) = blocks_read(bytes, Needs::ALL);
    rows_of(&blocks, &chroms)
}

/// The error of a write that has to be refused.
fn refused(written: Result<Vec<u8>>) -> Error {
    match written {
        Ok(bytes) => panic!("a file of {} bytes was written", bytes.len()),
        Err(error) => error,
    }
}

#[test]
fn a_position_that_goes_down_is_not_sorted_and_the_error_names_both_variants() {
    let error = refused(written(&DOWN, 2, 2, false));

    let Error::VarsVariantNotSorted {
        var,
        ref chrom,
        pos,
        var_before,
        ref chrom_before,
        pos_before,
    } = error
    else {
        panic!("the error of the order: {error:?}");
    };
    assert_eq!(
        (
            var,
            chrom.as_str(),
            pos,
            var_before,
            chrom_before.as_str(),
            pos_before
        ),
        (2, "chr1", 100, 1, "chr1", 5000)
    );
    assert!(
        error
            .to_string()
            .starts_with("variant 2, chr1 100, comes after variant 1, chr1 5000; "),
        "{error}"
    );
}

/// The check runs over the whole file: the variant out of order is the
/// first of the second batch, and its number is counted over the file and
/// not inside its batch.
#[test]
fn a_chromosome_that_comes_back_is_not_sorted_and_the_error_names_both_variants() {
    let error = refused(written(&COMES_BACK, 2, 3, false));

    let Error::VarsVariantNotSorted {
        var,
        ref chrom,
        pos,
        var_before,
        ref chrom_before,
        pos_before,
    } = error
    else {
        panic!("the error of the order: {error:?}");
    };
    assert_eq!(
        (
            var,
            chrom.as_str(),
            pos,
            var_before,
            chrom_before.as_str(),
            pos_before
        ),
        (4, "chr1", 150, 3, "chr2", 300)
    );
    assert!(
        error
            .to_string()
            .starts_with("variant 4, chr1 150, comes after variant 3, chr2 300; "),
        "{error}"
    );
    assert!(
        error.to_string().contains("chr1, which had ended"),
        "{error}"
    );
}

/// The same variants with `allow_unsorted` are written, and read back in
/// the order they came.
#[test]
fn variants_not_sorted_are_written_in_the_order_they_come_with_allow_unsorted() {
    for (rows, given_in, num_vars_per_block) in [(&DOWN[..], 2, 2), (&COMES_BACK[..], 2, 3)] {
        let bytes = written(rows, given_in, num_vars_per_block, true).expect("the file");
        let expected: Vec<ReadRow> = rows.iter().map(row_read).collect();
        assert_eq!(rows_read(bytes), expected);
    }
}

/// Two variants at one position are in order, whichever comes first, and
/// the chromosomes can come in any order: chr2 before chr1 here, and the
/// second pair at one position spans two batches.
#[test]
fn variants_at_one_position_pass_the_check_of_not_sorted_whatever_the_order_of_the_chromosomes() {
    let bytes = written(&SORTED, 1, 3, false).expect("the file");

    let expected: Vec<ReadRow> = SORTED.iter().map(row_read).collect();
    assert_eq!(rows_read(bytes), expected);
}

/// A block without chromosomes and positions has no order to check, and
/// is written.
#[test]
fn blocks_without_positions_are_written_with_no_check_of_not_sorted() {
    let mut chroms = ChromTable::new();
    let blocks: Vec<Block> = COMES_BACK
        .chunks(2)
        .map(|chunk| {
            let chunk: Vec<&Row> = chunk.iter().collect();
            let mut block = block_of(&chunk, &mut chroms);
            block.chrom = None;
            block.pos = None;
            block
        })
        .collect();
    let reader = GivenBlocks::of(blocks, chroms);

    let (_, num_vars) =
        write_vars(reader, Vec::new(), Some(3), false).expect("the file was written");

    assert_eq!(num_vars, 4);
}
