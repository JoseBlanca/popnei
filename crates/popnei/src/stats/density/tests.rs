//! The tests of the density of the variants, each with `var_density` in its
//! name, against the numbers of "How it is verified" of the density in
//! `docs/specs/stats.md`, which tabix 1.24 counted.

use std::fs::File;
use std::io::{BufReader, Cursor};
use std::num::NonZeroU64;
use std::path::{Path, PathBuf};

use super::{
    DensityOfChrom, LengthsFrom, MAX_NUM_WINDOWS, VarDensity, calc_var_density, count_the_variant,
};
use crate::block::{Block, BlockReader, SourceHeader};
use crate::error::{Error, Result};
use crate::filters::{FilteringStats, RegionSelection};
use crate::io::vcf::{VcfOptions, VcfReader};
use crate::variant::{ChromTable, Needs};

/// One window as the tables of the spec give it: the chromosome, the start,
/// the end and the count.
type Window = (String, u64, u64, u32);

fn reference_vcf(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/reference/vcf")
        .join(name)
}

/// `write.vcf` of the VCF writer of `docs/specs/io_vcf.md`, read with the
/// default, so its line of chr1 250, which failed its filter, is not there:
/// chr1 100, 1000 and 1001, and chr2 1 and 1500, with the lengths of its
/// header, chr1 2000 and chr2 1500.
fn write_vcf() -> VcfReader<BufReader<File>> {
    VcfReader::from_path(&reference_vcf("write.vcf"), VcfOptions::default()).expect("write.vcf")
}

/// The same file with its two `##contig` lines taken out, so its header
/// gives no length.
fn write_vcf_without_its_contigs() -> VcfReader<Cursor<Vec<u8>>> {
    let text = std::fs::read_to_string(reference_vcf("write.vcf")).expect("write.vcf");
    let without: String = text
        .lines()
        .filter(|line| !line.starts_with("##contig"))
        .map(|line| format!("{line}\n"))
        .collect();
    VcfReader::new(Cursor::new(without.into_bytes()), VcfOptions::default())
        .expect("write.vcf without its ##contig lines")
}

/// `many.vcf`, 500 variants of chr1 and chr2 whose `##contig` lines have no
/// length, with every variant given.
fn many_vcf() -> VcfReader<BufReader<File>> {
    let options = VcfOptions {
        ploidy: 2,
        only_passed: false,
        num_vars_per_block: Some(64),
    };
    VcfReader::from_path(&reference_vcf("many.vcf"), options).expect("many.vcf")
}

/// A VCF of one individual whose header gives chr1 the length `length`,
/// with a variant of chr1 at each of `positions`.
fn a_vcf_of_chr1(length: u64, positions: &[u64]) -> VcfReader<Cursor<Vec<u8>>> {
    let mut text = format!(
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length={length}>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\n"
    );
    for pos in positions {
        text.push_str(&format!("chr1\t{pos}\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"));
    }
    VcfReader::new(Cursor::new(text.into_bytes()), VcfOptions::default()).expect("the VCF")
}

fn lengths(of_each: &[(&str, u64)]) -> Vec<(String, u64)> {
    of_each
        .iter()
        .map(|(chrom, length)| ((*chrom).to_owned(), *length))
        .collect()
}

/// Every window of `density`, in its order.
fn windows(density: &VarDensity) -> Vec<Window> {
    density
        .windows()
        .map(|window| {
            (
                window.chrom.to_owned(),
                window.start,
                window.end,
                window.num_vars,
            )
        })
        .collect()
}

/// The windows of `chrom` of `width` base pairs from the position 1, the
/// last one ending at `last_end`, with `counts`.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the windows of the tables of the spec, of 20000 base pairs at most"
)]
fn laid_end_to_end(chrom: &str, width: u64, last_end: u64, counts: &[u32]) -> Vec<Window> {
    let mut start = 1;
    let mut laid = Vec::new();
    for (index, count) in counts.iter().enumerate() {
        let end = if index + 1 == counts.len() {
            last_end
        } else {
            start + width - 1
        };
        laid.push((chrom.to_owned(), start, end, *count));
        start = end + 1;
    }
    laid
}

fn window(chrom: &str, start: u64, end: u64, num_vars: u32) -> Window {
    (chrom.to_owned(), start, end, num_vars)
}

/// A reader of the tests that gives the blocks it was built with, of one
/// haploid individual and no genotypes: the chromosome and the position of
/// each variant are all a density reads.
#[derive(Debug)]
struct Given {
    chroms: ChromTable,
    header: SourceHeader,
    individuals: Vec<String>,
    /// The blocks it has not given yet, the next one last.
    left: Vec<Block>,
    /// What it was last asked to fill.
    needs: Needs,
}

impl Given {
    /// The reader over blocks of the variants of `blocks`, each a chromosome
    /// by its name and a position, whose table has the names of `table` in
    /// that order before any of the blocks, and whose header gives
    /// `chrom_lengths`.
    fn of(table: &[&str], chrom_lengths: &[(&str, u64)], blocks: &[&[(&str, u64)]]) -> Given {
        let mut chroms = ChromTable::new();
        for name in table {
            chroms.intern(name);
        }
        let mut left: Vec<Block> = blocks
            .iter()
            .map(|variants| {
                let numbers = variants
                    .iter()
                    .map(|(name, _)| chroms.intern(name))
                    .collect();
                let positions = variants.iter().map(|(_, pos)| *pos).collect();
                block_of(numbers, positions)
            })
            .collect();
        left.reverse();
        Given {
            chroms,
            header: SourceHeader {
                individuals: vec!["a".to_owned()],
                chrom_lengths: lengths(chrom_lengths),
                vcf_meta_lines: None,
                keeps_passed: false,
            },
            individuals: vec!["a".to_owned()],
            left,
            needs: Needs::ALL,
        }
    }

    /// The reader over the variants of `variants`, in one block, with no
    /// length in its header.
    fn of_one_block(variants: &[(&str, u64)]) -> Given {
        Given::of(&[], &[], &[variants])
    }
}

fn block_of(chrom: Vec<u32>, pos: Vec<u64>) -> Block {
    Block {
        num_vars: pos.len(),
        num_individuals: 1,
        ploidy: 1,
        gts: Vec::new(),
        chrom: Some(chrom),
        pos: Some(pos),
        id: None,
        alleles: None,
        qual: None,
        passed: None,
        vcf_text: None,
    }
}

impl BlockReader for Given {
    fn next_block(&mut self) -> Result<Option<Block>> {
        Ok(self.left.pop())
    }

    fn individuals(&self) -> &[String] {
        &self.individuals
    }

    fn ploidy(&self) -> usize {
        1
    }

    fn chroms(&self) -> &ChromTable {
        &self.chroms
    }

    fn set_needs(&mut self, needs: Needs) {
        self.needs = needs;
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

fn density_of<R: BlockReader>(
    mut reader: R,
    window_size: u64,
    chrom_lengths: Option<&[(String, u64)]>,
) -> VarDensity {
    calc_var_density(&mut reader, window_size, chrom_lengths).expect("the density")
}

fn error_of<R: BlockReader>(
    mut reader: R,
    window_size: u64,
    chrom_lengths: Option<&[(String, u64)]>,
) -> Error {
    match calc_var_density(&mut reader, window_size, chrom_lengths) {
        Ok(density) => panic!("the density was given: {density:?}"),
        Err(error) => error,
    }
}

/// The worked example of the spec, in windows of 500 with the lengths of
/// the header: the variant at 1000 is on the last position of the second
/// window of chr1 and the one at 1001 on the first of the third, and chr2
/// has its variants on its first position and on its length.
#[test]
fn var_density_of_the_worked_example_in_windows_of_500_with_the_lengths_of_its_header() {
    let density = density_of(write_vcf(), 500, None);
    assert_eq!(
        windows(&density),
        [
            window("chr1", 1, 500, 1),
            window("chr1", 501, 1000, 1),
            window("chr1", 1001, 1500, 1),
            window("chr1", 1501, 2000, 0),
            window("chr2", 1, 500, 1),
            window("chr2", 501, 1000, 0),
            window("chr2", 1001, 1500, 1),
        ]
    );
    assert_eq!(density.window_size(), 500);
    assert_eq!(density.num_vars(), 5);
    assert_eq!(density.num_windows(), 7);
    assert_eq!(
        density.chroms(),
        [
            DensityOfChrom {
                name: "chr1".to_owned(),
                length: Some(2000),
                counts: vec![1, 1, 1, 0],
            },
            DensityOfChrom {
                name: "chr2".to_owned(),
                length: Some(1500),
                counts: vec![1, 0, 1],
            },
        ]
    );
}

/// The same in windows of 600, where the last window of chr1, 1801 to 2000,
/// and that of chr2, 1201 to 1500, are shorter than the others.
#[test]
fn var_density_of_the_worked_example_in_windows_of_600_with_the_lengths_of_its_header() {
    let density = density_of(write_vcf(), 600, None);
    assert_eq!(
        windows(&density),
        [
            window("chr1", 1, 600, 1),
            window("chr1", 601, 1200, 2),
            window("chr1", 1201, 1800, 0),
            window("chr1", 1801, 2000, 0),
            window("chr2", 1, 600, 1),
            window("chr2", 601, 1200, 0),
            window("chr2", 1201, 1500, 1),
        ]
    );
    assert_eq!(density.num_vars(), 5);
}

/// Without the lengths, the windows go up to the one of the last variant,
/// which ends at its full width, past the variant.
#[test]
fn var_density_of_the_worked_example_without_its_contig_lines_in_windows_of_500() {
    let density = density_of(write_vcf_without_its_contigs(), 500, None);
    assert_eq!(
        windows(&density),
        [
            window("chr1", 1, 500, 1),
            window("chr1", 501, 1000, 1),
            window("chr1", 1001, 1500, 1),
            window("chr2", 1, 500, 1),
            window("chr2", 501, 1000, 0),
            window("chr2", 1001, 1500, 1),
        ]
    );
    assert!(density.chroms().iter().all(|chrom| chrom.length.is_none()));
    assert_eq!(density.num_vars(), 5);
    assert_eq!(density.num_windows(), 6);
}

#[test]
fn var_density_of_the_worked_example_without_its_contig_lines_in_windows_of_600() {
    let density = density_of(write_vcf_without_its_contigs(), 600, None);
    assert_eq!(
        windows(&density),
        [
            window("chr1", 1, 600, 1),
            window("chr1", 601, 1200, 2),
            window("chr2", 1, 600, 1),
            window("chr2", 601, 1200, 0),
            window("chr2", 1201, 1800, 1),
        ]
    );
}

/// The first table of "How it is verified", which tabix 1.24 counted:
/// windows of 1000 of `many.vcf`, whose `##contig` lines have no length.
#[test]
fn var_density_of_many_vcf_in_windows_of_1000_is_the_one_tabix_counts() {
    let density = density_of(many_vcf(), 1000, None);
    let mut chr1 = vec![1];
    chr1.extend([27; 9]);
    chr1.push(6);
    let mut chr2 = vec![0; 10];
    chr2.push(21);
    chr2.extend([27; 8]);
    chr2.push(13);
    let mut expected = laid_end_to_end("chr1", 1000, 11000, &chr1);
    expected.extend(laid_end_to_end("chr2", 1000, 20000, &chr2));
    assert_eq!(windows(&density), expected);
    assert_eq!(density.num_vars(), 500);
}

/// The second: chr1 has one window more, 11001 to 12000, with 0, and the
/// last window of chr2 ends at its length, 19500, still with 13.
#[test]
fn var_density_of_many_vcf_with_chrom_lengths_ends_each_chromosome_at_its_length() {
    let given = lengths(&[("chr1", 12000), ("chr2", 19500)]);
    let density = density_of(many_vcf(), 1000, Some(&given));
    let mut chr1 = vec![1];
    chr1.extend([27; 9]);
    chr1.extend([6, 0]);
    let mut chr2 = vec![0; 10];
    chr2.push(21);
    chr2.extend([27; 8]);
    chr2.push(13);
    let mut expected = laid_end_to_end("chr1", 1000, 12000, &chr1);
    expected.extend(laid_end_to_end("chr2", 1000, 19500, &chr2));
    assert_eq!(windows(&density), expected);
    assert_eq!(density.chroms()[0].length, Some(12000));
}

#[test]
fn var_density_refuses_a_variant_past_a_length_of_chrom_lengths_and_names_it() {
    let given = lengths(&[("chr1", 10000), ("chr2", 20000)]);
    let error = error_of(many_vcf(), 1000, Some(&given));
    assert!(
        matches!(
            &error,
            Error::VarDensityVarPastTheLength { chrom, pos: 10028, length: 10000, from: LengthsFrom::ChromLengths }
                if chrom == "chr1"
        ),
        "{error:?}"
    );
    let message = error.to_string();
    for part in ["chr1", "10028", "10000", "`chrom_lengths`"] {
        assert!(message.contains(part), "{message}");
    }
    assert!(error.names_the_file());
}

/// A variant on the length is in the last window, and one past it by 1,
/// the telomere of the VCF format, is refused with the length of the
/// header.
#[test]
fn var_density_refuses_a_variant_past_a_length_of_the_header_and_takes_one_on_it() {
    let density = density_of(a_vcf_of_chr1(1000, &[1000]), 300, None);
    assert_eq!(
        windows(&density),
        [
            window("chr1", 1, 300, 0),
            window("chr1", 301, 600, 0),
            window("chr1", 601, 900, 0),
            window("chr1", 901, 1000, 1),
        ]
    );
    let error = error_of(a_vcf_of_chr1(1000, &[1000, 1001]), 300, None);
    assert!(
        matches!(
            &error,
            Error::VarDensityVarPastTheLength { chrom, pos: 1001, length: 1000, from: LengthsFrom::Source }
                if chrom == "chr1"
        ),
        "{error:?}"
    );
    assert!(
        error.to_string().contains("the header of the source"),
        "{error}"
    );
}

/// The positions on both sides of each edge of windows of 500, and windows
/// of 1 base pair, in which each position is a window of its own.
#[test]
fn var_density_puts_a_variant_on_an_edge_in_the_window_that_ends_or_starts_there() {
    let on_the_edges = [
        ("chr1", 1),
        ("chr1", 499),
        ("chr1", 500),
        ("chr1", 501),
        ("chr1", 999),
        ("chr1", 1000),
        ("chr1", 1001),
    ];
    let density = density_of(Given::of_one_block(&on_the_edges), 500, None);
    assert_eq!(
        windows(&density),
        [
            window("chr1", 1, 500, 3),
            window("chr1", 501, 1000, 3),
            window("chr1", 1001, 1500, 1),
        ]
    );
    let given = lengths(&[("chr1", 1001)]);
    let density = density_of(Given::of_one_block(&on_the_edges), 500, Some(&given));
    assert_eq!(
        windows(&density),
        [
            window("chr1", 1, 500, 3),
            window("chr1", 501, 1000, 3),
            window("chr1", 1001, 1001, 1),
        ]
    );
    let density = density_of(
        Given::of_one_block(&[("chr1", 4), ("chr1", 2), ("chr1", 1), ("chr1", 2)]),
        1,
        None,
    );
    assert_eq!(
        windows(&density),
        [
            window("chr1", 1, 1, 1),
            window("chr1", 2, 2, 2),
            window("chr1", 3, 3, 0),
            window("chr1", 4, 4, 1),
        ]
    );
}

/// The variants need not be sorted: a window before the last one made is
/// counted into, and one past it is made with those between.
#[test]
fn var_density_counts_variants_that_come_in_any_order() {
    let density = density_of(
        Given::of(
            &[],
            &[],
            &[
                &[("chr1", 1600), ("chr2", 3)],
                &[("chr1", 100), ("chr1", 1700)],
            ],
        ),
        500,
        None,
    );
    assert_eq!(
        windows(&density),
        [
            window("chr1", 1, 500, 1),
            window("chr1", 501, 1000, 0),
            window("chr1", 1001, 1500, 0),
            window("chr1", 1501, 2000, 2),
            window("chr2", 1, 500, 1),
        ]
    );
    assert_eq!(density.num_vars(), 4);
}

/// The chromosomes with a length come first, in the order of the lengths,
/// one with no variant among them with all its windows at 0; then those
/// with no length in the order their first variant came, which is not the
/// order of the table of the reader.
#[test]
fn var_density_gives_the_chromosomes_of_the_lengths_first_and_then_by_their_first_variant() {
    let reader = Given::of(
        &["chr1", "chr2", "chr3"],
        &[("chrX", 1200), ("chr3", 1000)],
        &[&[("chr2", 5), ("chr1", 7)], &[("chr3", 1), ("chr2", 6)]],
    );
    let density = density_of(reader, 500, None);
    assert_eq!(
        windows(&density),
        [
            window("chrX", 1, 500, 0),
            window("chrX", 501, 1000, 0),
            window("chrX", 1001, 1200, 0),
            window("chr3", 1, 500, 1),
            window("chr3", 501, 1000, 0),
            window("chr2", 1, 500, 2),
            window("chr1", 1, 500, 1),
        ]
    );
}

/// `chrom_lengths` replaces the lengths of the header for every chromosome:
/// chr1, which it does not name, has no length, and the windows of chr2 go
/// to the length it gives.
#[test]
fn var_density_with_chrom_lengths_reads_no_length_of_the_header() {
    let given = lengths(&[("chr2", 3000)]);
    let density = density_of(write_vcf(), 500, Some(&given));
    let mut expected = laid_end_to_end("chr2", 500, 3000, &[1, 0, 1, 0, 0, 0]);
    expected.extend(laid_end_to_end("chr1", 500, 1500, &[1, 1, 1]));
    assert_eq!(windows(&density), expected);
    assert_eq!(density.chroms()[1].length, None);
}

#[test]
fn var_density_of_a_window_size_of_0_is_refused() {
    let error = error_of(write_vcf(), 0, None);
    assert!(
        matches!(error, Error::VarDensityWindowSizeZero),
        "{error:?}"
    );
    assert!(error.to_string().contains("`window_size`"), "{error}");
    assert!(!error.names_the_file());
}

#[test]
fn var_density_refuses_a_length_of_0_and_a_chromosome_named_twice() {
    let given = lengths(&[("chr1", 2000), ("chr2", 0)]);
    let error = error_of(write_vcf(), 500, Some(&given));
    assert!(
        matches!(&error, Error::VarDensityChromLengthZero { chrom, from: LengthsFrom::ChromLengths } if chrom == "chr2"),
        "{error:?}"
    );
    assert!(!error.names_the_file());

    let reader = Given::of(&[], &[("chr1", 0)], &[&[("chr1", 1)]]);
    let error = error_of(reader, 500, None);
    assert!(
        matches!(&error, Error::VarDensityChromLengthZero { chrom, from: LengthsFrom::Source } if chrom == "chr1"),
        "{error:?}"
    );
    assert!(error.names_the_file());

    let given = lengths(&[("chr1", 2000), ("chr2", 1500), ("chr1", 3000)]);
    let error = error_of(write_vcf(), 500, Some(&given));
    assert!(
        matches!(&error, Error::VarDensityChromLengthTwice { chrom, from: LengthsFrom::ChromLengths } if chrom == "chr1"),
        "{error:?}"
    );
    assert!(!error.names_the_file());
}

#[test]
fn var_density_refuses_a_variant_at_position_0() {
    let error = error_of(Given::of_one_block(&[("chr1", 5), ("chr7", 0)]), 500, None);
    assert!(
        matches!(&error, Error::VarDensityVarAtPositionZero { chrom } if chrom == "chr7"),
        "{error:?}"
    );
}

/// Lengths that give more windows than the bound are refused before the
/// pass reads a block; the bound itself is taken.
#[test]
fn var_density_of_more_windows_than_max_num_windows_is_refused_before_the_pass() {
    assert_eq!(MAX_NUM_WINDOWS, 10_000_000);
    let mut reader = Given::of(&[], &[("chr1", 10_000_001)], &[&[("chr1", 1)]]);
    let error = match calc_var_density(&mut reader, 1, None) {
        Ok(_) => panic!("the density was given"),
        Err(error) => error,
    };
    assert!(
        matches!(
            error,
            Error::VarDensityTooManyWindows {
                num_windows: 10_000_001,
                window_size: 1,
                largest: 10_000_000
            }
        ),
        "{error:?}"
    );
    assert_eq!(reader.left.len(), 1, "no block was read");

    let given = lengths(&[("chr1", 5_000_000), ("chr2", 5_000_001)]);
    let error = error_of(Given::of_one_block(&[("chr1", 1)]), 1, Some(&given));
    assert!(
        matches!(
            error,
            Error::VarDensityTooManyWindows {
                num_windows: 10_000_001,
                ..
            }
        ),
        "{error:?}"
    );

    let given = lengths(&[("chr1", 5_000_000), ("chr2", 5_000_000)]);
    let density = density_of(Given::of_one_block(&[("chr2", 5_000_000)]), 1, Some(&given));
    assert_eq!(density.num_windows(), 10_000_000);
    assert_eq!(density.chroms()[1].counts.last(), Some(&1));
}

/// Without lengths the windows are counted as the variants come, over all
/// the chromosomes, and the pass is refused when they pass the bound.
#[test]
fn var_density_of_more_windows_than_max_num_windows_is_refused_as_the_pass_reaches_them() {
    let error = error_of(
        Given::of_one_block(&[("chr1", 6_000_000), ("chr2", 4_000_001)]),
        1,
        None,
    );
    assert!(
        matches!(
            error,
            Error::VarDensityTooManyWindows {
                num_windows: 10_000_001,
                window_size: 1,
                ..
            }
        ),
        "{error:?}"
    );
    let density = density_of(
        Given::of_one_block(&[("chr1", 6_000_000), ("chr2", 4_000_000), ("chr1", 3)]),
        1,
        None,
    );
    assert_eq!(density.num_windows(), 10_000_000);
    assert_eq!(density.num_vars(), 3);
}

/// Lengths of up to `u64::MAX`, which the VCF reader takes, give numbers of
/// windows past a `usize` of wasm and past a `u64` when added up: each is
/// refused with the most the count reaches, and none overflows.
#[test]
fn var_density_of_lengths_up_to_the_largest_u64_neither_overflows_nor_wraps() {
    let given = lengths(&[("chr1", u64::MAX), ("chr2", u64::MAX)]);
    let error = error_of(Given::of_one_block(&[("chr1", 1)]), 1, Some(&given));
    assert!(
        matches!(
            error,
            Error::VarDensityTooManyWindows {
                num_windows: u64::MAX,
                ..
            }
        ),
        "{error:?}"
    );
    let given = lengths(&[("chr1", 4_294_967_297)]);
    let error = error_of(Given::of_one_block(&[("chr1", 1)]), 1, Some(&given));
    assert!(
        matches!(
            error,
            Error::VarDensityTooManyWindows {
                num_windows: 4_294_967_297,
                ..
            }
        ),
        "{error:?}"
    );

    let half = 1_u64 << 63;
    let given = lengths(&[("chr1", u64::MAX)]);
    let density = density_of(
        Given::of_one_block(&[("chr1", u64::MAX), ("chr1", half), ("chr1", half + 1)]),
        half,
        Some(&given),
    );
    assert_eq!(
        windows(&density),
        [
            window("chr1", 1, half, 1),
            window("chr1", half + 1, u64::MAX, 2),
        ]
    );
}

/// Positions past 2^53, where a float64 holds every other whole number: the
/// variants at 2^53 and 2^53 + 1 are one position apart and in two windows.
#[test]
fn var_density_of_positions_past_2_to_the_53_keeps_them_exact() {
    let two_to_the_53 = 1_u64 << 53;
    let given = lengths(&[("chr1", two_to_the_53 + 1)]);
    let density = density_of(
        Given::of_one_block(&[("chr1", two_to_the_53), ("chr1", two_to_the_53 + 1)]),
        two_to_the_53 / 2,
        Some(&given),
    );
    assert_eq!(
        windows(&density),
        [
            window("chr1", 1, two_to_the_53 / 2, 0),
            window("chr1", two_to_the_53 / 2 + 1, two_to_the_53, 1),
            window("chr1", two_to_the_53 + 1, two_to_the_53 + 1, 1),
        ]
    );
}

/// The last window of a chromosome with no length whose full width would
/// pass the largest position ends at it, as the spec decides.
#[test]
fn var_density_of_a_last_window_whose_full_width_passes_the_largest_position_ends_there() {
    let width = 10_000_000_000_000_000_000_u64;
    let density = density_of(
        Given::of_one_block(&[("chr1", 15_000_000_000_000_000_000)]),
        width,
        None,
    );
    assert_eq!(
        windows(&density),
        [
            window("chr1", 1, width, 0),
            window("chr1", width + 1, u64::MAX, 1),
        ]
    );
    let density = density_of(Given::of_one_block(&[("chr1", u64::MAX)]), u64::MAX, None);
    assert_eq!(windows(&density), [window("chr1", 1, u64::MAX, 1)]);
}

/// A count at the most a `u32` holds refuses one more variant, and names
/// the window, which ends at the length of a chromosome that has one.
#[test]
fn var_density_refuses_a_window_of_more_variants_than_its_count_holds() {
    let mut chrom = DensityOfChrom {
        name: "chr1".to_owned(),
        length: Some(700),
        counts: vec![0, u32::MAX],
    };
    let mut num_windows = 2;
    let window_size = NonZeroU64::new(500).expect("500");
    count_the_variant(
        &mut chrom,
        500,
        window_size,
        LengthsFrom::Source,
        &mut num_windows,
    )
    .expect("a variant in the first window");
    let error = match count_the_variant(
        &mut chrom,
        700,
        window_size,
        LengthsFrom::Source,
        &mut num_windows,
    ) {
        Ok(()) => panic!("the count went past u32::MAX"),
        Err(error) => error,
    };
    assert!(
        matches!(
            &error,
            Error::VarDensityWindowTooFull { chrom, start: 501, end: 700, largest: u32::MAX }
                if chrom == "chr1"
        ),
        "{error:?}"
    );
    assert_eq!(chrom.counts, [1, u32::MAX]);
}

#[test]
fn var_density_asks_the_reader_for_the_chromosome_and_the_position_alone() {
    let mut reader = Given::of_one_block(&[("chr1", 1)]);
    calc_var_density(&mut reader, 10, None).expect("the density");
    assert_eq!(reader.needs, Needs::CHROM_POS);
}

/// A source with no positions, a vars file written without them or
/// variants built from an array of genotypes, gives blocks with neither
/// column.
#[test]
fn var_density_of_a_source_with_no_positions_is_refused() {
    let mut reader = Given::of_one_block(&[("chr1", 1)]);
    for block in &mut reader.left {
        block.chrom = None;
        block.pos = None;
    }
    let error = error_of(reader, 10, None);
    assert!(
        matches!(error, Error::FieldsNotInTheBlock { fields } if fields == Needs::CHROM_POS),
        "{error:?}"
    );
}

/// A pass that gave no variant is refused with lengths as without them,
/// although the lengths alone would give windows.
#[test]
fn var_density_of_a_pass_that_gave_no_variant_is_refused_whatever_the_lengths() {
    for chrom_lengths in [&[][..], &[("chr1", 2000)][..]] {
        let reader = Given::of(&[], chrom_lengths, &[]);
        let error = error_of(reader, 500, None);
        assert!(
            matches!(
                error,
                Error::PassGaveNoVariant {
                    num_vars_of_the_source: 0,
                    ..
                }
            ),
            "{error:?}"
        );
    }
}

#[test]
fn var_density_refuses_a_chromosome_number_its_reader_has_no_name_for() {
    let mut reader = Given::of_one_block(&[("chr1", 1)]);
    for block in &mut reader.left {
        block.chrom = Some(vec![3]);
    }
    let error = error_of(reader, 10, None);
    assert!(
        matches!(error, Error::VarDensityChromNameMissing { number: 3 }),
        "{error:?}"
    );
}

/// A block whose positions are fewer than its variants would pair a
/// chromosome with the position of another variant.
#[test]
fn var_density_refuses_a_block_whose_positions_are_not_one_for_each_variant() {
    let mut reader = Given::of_one_block(&[("chr1", 1), ("chr1", 2)]);
    for block in &mut reader.left {
        block.pos = Some(vec![1]);
    }
    let error = error_of(reader, 10, None);
    assert!(
        matches!(
            error,
            Error::BlockArrayOfAnotherSize {
                array: "pos",
                found: 1,
                expected: 2
            }
        ),
        "{error:?}"
    );
}

/// The chromosomes of the lengths come in the order of the lengths, which
/// is not the order of their names nor that of the variants of the source.
#[test]
fn var_density_gives_the_chromosomes_in_the_order_of_chrom_lengths() {
    let given = lengths(&[("chr2", 1500), ("chr1", 2000)]);
    let density = density_of(write_vcf(), 500, Some(&given));
    let mut expected = laid_end_to_end("chr2", 500, 1500, &[1, 0, 1]);
    expected.extend(laid_end_to_end("chr1", 500, 2000, &[1, 1, 1, 0]));
    assert_eq!(windows(&density), expected);
}
