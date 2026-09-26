//! The tests of the filter by regions, each with `by_regions` in its name,
//! against the numbers of "How it is verified" of the filter by regions of
//! `docs/specs/filters.md` and the positions bcftools 1.24 and plink2
//! v2.0.0-a.7.7 keep, which `tests/reference/filters/make_reference.py`
//! stored.

use std::fs::File;
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use flate2::Compression;
use flate2::write::GzEncoder;

use super::{
    BedLineProblem, RegionFilter, RegionSelection, Regions, RegionsReader, keep_of_the_rows,
    keep_of_the_rows_one_by_one,
};
use crate::block::{Block, BlockReader, SourceHeader};
use crate::error::{Error, Result};
use crate::filters::{
    FilteredReader, FilteringStats, IndividualsReader, PassStep, VarFilter, VarFilteringCriterion,
    chain_of, refuse_a_second_filter_of_a_kind,
};
use crate::io::vcf::{VcfOptions, VcfReader};
use crate::variant::{ChromTable, Needs};

/// A variant named as a user reads it: its chromosome and its position.
type Named = (String, u64);

fn reference_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/reference")
}

/// The BED of "How it is verified", with its `track` line and its comment.
fn the_bed_of_the_spec() -> Vec<u8> {
    std::fs::read(reference_dir().join("filters/regions.bed")).expect("regions.bed")
}

/// The BED of the worked example.
fn the_bed_of_the_worked_example() -> Vec<u8> {
    std::fs::read(reference_dir().join("filters/write_regions.bed")).expect("write_regions.bed")
}

fn selection_of(bed: &[u8], exclude: bool) -> RegionSelection {
    RegionSelection {
        regions: Arc::new(Regions::from_bed(bed).expect("the regions of the BED")),
        exclude,
    }
}

/// A reader of the VCF of that name under `tests/reference/vcf`, with every
/// variant given, the ones that failed their FILTER too, as bcftools and
/// plink2 read it.
fn vcf_reader(name: &str, num_vars_per_block: Option<usize>) -> VcfReader<BufReader<File>> {
    let options = VcfOptions {
        ploidy: 2,
        only_passed: false,
        num_vars_per_block,
    };
    VcfReader::from_path(&reference_dir().join("vcf").join(name), options).expect("the reader")
}

/// The variants of `blocks`, by the names `chroms` gives their numbers.
fn named_of(blocks: &[Block], chroms: &ChromTable) -> Vec<Named> {
    blocks
        .iter()
        .flat_map(|block| {
            let chrom = block.chrom.clone().expect("the chromosomes");
            let pos = block.pos.clone().expect("the positions");
            chrom
                .into_iter()
                .zip(pos)
                .map(|(number, pos)| (chroms.name(number).expect("a name").to_owned(), pos))
                .collect::<Vec<Named>>()
        })
        .collect()
}

/// The variants bcftools and plink2 kept, from the file of that name under
/// `tests/reference/filters`, a chromosome and a position on each line.
fn the_reference(name: &str) -> Vec<Named> {
    let path = reference_dir().join("filters").join(format!("{name}.txt"));
    std::fs::read_to_string(&path)
        .expect("the file of the reference")
        .lines()
        .map(|line| {
            let (chrom, pos) = line.split_once('\t').expect("a chromosome and a position");
            (chrom.to_owned(), pos.parse().expect("a position"))
        })
        .collect()
}

fn named(chrom: &str, pos: u64) -> Named {
    (chrom.to_owned(), pos)
}

fn pair(vars_processed: u64, vars_kept: u64) -> FilteringStats {
    FilteringStats {
        vars_processed,
        vars_kept,
    }
}

/// Every block a reader gives, until it has no more or it fails.
fn blocks_of(reader: &mut impl BlockReader) -> Result<Vec<Block>> {
    let mut blocks = Vec::new();
    while let Some(block) = reader.next_block()? {
        blocks.push(block);
    }
    Ok(blocks)
}

/// A source over a reader that takes the regions it is offered and passes
/// over variants the selection does not keep, counting them, which is what
/// a source that skips does: so the filter is tested over a source that
/// skips, before work package 5 gives the VCF reader the skip.
struct SkippingSource<R: BlockReader> {
    reader: R,
    selection: Option<RegionSelection>,
    num_skipped: u64,
    /// How many times it was offered the regions, held by the test too.
    offers: Arc<AtomicUsize>,
    /// Which of the variants the selection does not keep it passes over.
    skips: TheSkip,
}

/// Which variants a [`SkippingSource`] passes over.
#[derive(Debug, Clone, Copy)]
enum TheSkip {
    /// Every variant the selection does not keep, as the VCF reader of work
    /// package 5 will.
    EveryVariant,
    /// The blocks of which the selection keeps no variant, whole, and none
    /// of the other blocks, as the vars file reader of work package 5 will
    /// skip a batch: the filter has to take out the rest.
    WholeBlocks,
}

impl<R: BlockReader> SkippingSource<R> {
    fn over(reader: R) -> SkippingSource<R> {
        SkippingSource::of(reader, TheSkip::EveryVariant)
    }

    fn of(reader: R, skips: TheSkip) -> SkippingSource<R> {
        SkippingSource {
            reader,
            selection: None,
            num_skipped: 0,
            offers: Arc::new(AtomicUsize::new(0)),
            skips,
        }
    }
}

impl<R: BlockReader> BlockReader for SkippingSource<R> {
    fn next_block(&mut self) -> Result<Option<Block>> {
        loop {
            let Some(mut block) = self.reader.next_block()? else {
                return Ok(None);
            };
            let Some(selection) = &self.selection else {
                return Ok(Some(block));
            };
            let chroms = self.reader.chroms();
            let chrom = block.chrom.clone().expect("the chromosomes");
            let pos = block.pos.clone().expect("the positions");
            let keep: Vec<bool> = chrom
                .iter()
                .zip(&pos)
                .map(|(number, pos)| selection.keeps(chroms.name(*number).expect("a name"), *pos))
                .collect();
            let keep = match self.skips {
                TheSkip::EveryVariant => keep,
                TheSkip::WholeBlocks => {
                    let any_kept = keep.iter().any(|keep_it| *keep_it);
                    vec![any_kept; keep.len()]
                }
            };
            let skipped = keep.iter().filter(|keep_it| !**keep_it).count();
            self.num_skipped = self
                .num_skipped
                .checked_add(u64::try_from(skipped).expect("a count"))
                .expect("a count of the variants of a file");
            block.retain_vars(&keep)?;
            if block.num_vars > 0 {
                return Ok(Some(block));
            }
        }
    }

    fn individuals(&self) -> &[String] {
        self.reader.individuals()
    }

    fn ploidy(&self) -> usize {
        self.reader.ploidy()
    }

    fn chroms(&self) -> &ChromTable {
        self.reader.chroms()
    }

    fn set_needs(&mut self, needs: Needs) {
        self.reader.set_needs(needs.union(Needs::CHROM_POS));
    }

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        self.reader.filtering_stats()
    }

    fn header(&self) -> &SourceHeader {
        self.reader.header()
    }

    fn skip_outside(&mut self, selection: RegionSelection) -> bool {
        self.offers.fetch_add(1, Ordering::SeqCst);
        self.selection = Some(selection);
        true
    }

    fn num_skipped(&self) -> u64 {
        self.num_skipped
    }
}

/// The worked example: the six lines of `write.vcf` with the four regions
/// of positions 100, 251, 1000 and chr2 1. `chr1 99 100` is 100 alone,
/// `chr1 999 1000` holds 1000 and not 1001, and `chr1 250 251` is 251, so
/// the deletion at chr1 250 is outside.
#[test]
fn the_worked_example_by_regions_keeps_what_bcftools_and_plink2_keep() {
    for (exclude, name, expected) in [
        (
            false,
            "write.regions",
            [named("chr1", 100), named("chr1", 1000), named("chr2", 1)],
        ),
        (
            true,
            "write.excluded_regions",
            [named("chr1", 250), named("chr1", 1001), named("chr2", 1500)],
        ),
    ] {
        let mut reader = vcf_reader("write.vcf", None);
        let mut blocks = blocks_of(&mut reader).expect("the blocks of write.vcf");
        let mut filter = RegionFilter::new(selection_of(&the_bed_of_the_worked_example(), exclude));
        for block in &mut blocks {
            filter
                .filter_block(block, reader.chroms())
                .expect("the filter");
        }
        let kept = named_of(&blocks, reader.chroms());
        assert_eq!(kept, expected, "exclude {exclude}");
        assert_eq!(kept, the_reference(name), "exclude {exclude}");
        assert_eq!(filter.stats(), pair(6, 3), "exclude {exclude}");
        // Every column of a variant stays with it: the text of its line too.
        for block in &blocks {
            assert!(block.check().is_ok());
        }
    }
}

/// BED counts from 0 and leaves the end out, a VCF counts from 1: the line
/// `chr1 10 20` holds 11 and 20, and not 10 nor 21. A variant put on the
/// wrong side of an edge would be kept or dropped and nothing would say so.
#[test]
fn a_variant_on_either_edge_of_a_region_by_regions_is_on_the_side_bed_puts_it() {
    let regions = Regions::from_bed(b"chr1\t10\t20\n".as_slice()).expect("the regions");
    assert!(!regions.contains("chr1", 10));
    assert!(regions.contains("chr1", 11));
    assert!(regions.contains("chr1", 20));
    assert!(!regions.contains("chr1", 21));
    for (exclude, expected, gts) in [
        (false, vec![11, 20], vec![1, 0]),
        (true, vec![10, 21], vec![0, 1]),
    ] {
        let mut chroms = ChromTable::new();
        chroms.intern("chr1");
        let mut block = Block {
            num_vars: 4,
            num_individuals: 1,
            ploidy: 1,
            gts: vec![0, 1, 0, 1],
            chrom: Some(vec![0; 4]),
            pos: Some(vec![10, 11, 20, 21]),
            id: None,
            alleles: None,
            qual: None,
            vcf_text: None,
        };
        let mut filter = RegionFilter::new(selection_of(b"chr1\t10\t20\n", exclude));
        filter
            .filter_block(&mut block, &chroms)
            .expect("the filter");
        assert_eq!(block.pos, Some(expected), "exclude {exclude}");
        assert_eq!(block.gts, gts, "exclude {exclude}");
        assert_eq!(filter.stats(), pair(4, 2));
    }
}

/// The six regions of the BED of the spec join into five, and each edge of
/// them is on the side the spec says.
#[test]
fn the_bed_of_the_spec_by_regions_joins_six_regions_into_five() {
    let regions = Regions::from_bed(the_bed_of_the_spec().as_slice()).expect("the regions");
    assert_eq!(regions.num_regions(), 5);
    for (chrom, pos, inside) in [
        ("chr1", 1, true),
        ("chr1", 2000, true),
        ("chr1", 2001, false),
        ("chr1", 4990, false),
        ("chr1", 4991, true),
        ("chr1", 5050, true),
        ("chr1", 5051, true),
        ("chr1", 5100, true),
        ("chr1", 5101, false),
        ("chr2", 10249, false),
        ("chr2", 10250, true),
        ("chr2", 10251, false),
        ("chr2", 19000, false),
        ("chr2", 19001, true),
        ("chr2", 30000, true),
        ("chr2", 30001, false),
        ("chr3", 1, true),
        ("chr3", 100_000, true),
        ("chr3", 100_001, false),
        // The name is compared as it is written: `1` is not `chr1`.
        ("1", 1, false),
        ("chr4", 1, false),
    ] {
        assert_eq!(regions.contains(chrom, pos), inside, "{chrom} {pos}");
    }
}

/// Two regions that touch are one, and one inside another is the larger.
#[test]
fn regions_that_touch_or_nest_by_regions_are_joined() {
    let regions = Regions::from_bed(b"c\t0\t10\nc\t10\t20\nc\t2\t5\nc\t21\t30\n".as_slice())
        .expect("the regions");
    // 1 to 10 and 11 to 20 touch, 3 to 5 is inside, 22 to 30 is apart.
    assert_eq!(regions.num_regions(), 2);
    // The first positions of the region that holds the nested one, which a
    // sort by the last position would put after it and lose.
    assert!(regions.contains("c", 1));
    assert!(regions.contains("c", 2));
    assert!(regions.contains("c", 3));
    assert!(regions.contains("c", 20));
    assert!(!regions.contains("c", 21));
    assert!(regions.contains("c", 22));
}

/// The 45 variants of `many.vcf` that bcftools and plink2 keep with the BED
/// of the spec, and the 455 they keep with `exclude`, in blocks of 7 and of
/// the default size, over a source that skips and one that does not: the
/// same variants and the counts of 500 given, whichever the source did.
#[test]
fn many_vcf_by_regions_keeps_the_45_and_excludes_to_the_455_of_bcftools_and_plink2() {
    for (exclude, name, num_kept) in [(false, "regions", 45_u64), (true, "excluded_regions", 455)] {
        let expected = the_reference(name);
        assert_eq!(u64::try_from(expected.len()).unwrap(), num_kept);
        for num_vars_per_block in [Some(7), None] {
            for skips in [
                None,
                Some(TheSkip::EveryVariant),
                Some(TheSkip::WholeBlocks),
            ] {
                let selection = selection_of(&the_bed_of_the_spec(), exclude);
                let what = format!("exclude {exclude}, {num_vars_per_block:?}, skips {skips:?}");
                let vcf = vcf_reader("many.vcf", num_vars_per_block);
                let source: Box<dyn BlockReader> = match skips {
                    Some(skips) => Box::new(SkippingSource::of(vcf, skips)),
                    None => Box::new(vcf),
                };
                let mut reader =
                    RegionsReader::new(source, RegionFilter::new(selection)).expect("the reader");
                // The consumer asks for the genotypes alone, and the filter
                // asks for the chromosome and the position besides.
                reader.set_needs(Needs::GTS);
                let blocks = blocks_of(&mut reader).expect("the blocks");
                for block in &blocks {
                    assert!(block.num_vars > 0, "{what}");
                    if let Some(size) = num_vars_per_block {
                        assert!(block.num_vars <= size, "{what}");
                    }
                }
                assert_eq!(named_of(&blocks, reader.chroms()), expected, "{what}");
                assert_eq!(
                    reader.filtering_stats(),
                    vec![(name, pair(500, num_kept))],
                    "{what}"
                );
            }
        }
    }
}

/// The 45 of the spec by what they are: 28 of chr1 up to 2000, the three at
/// chr1 4996, 5033 and 5070, chr2 10250, and 13 of chr2 from 19001.
#[test]
fn the_45_of_many_vcf_by_regions_are_those_the_spec_names() {
    let mut reader = RegionsReader::new(
        vcf_reader("many.vcf", None),
        RegionFilter::new(selection_of(&the_bed_of_the_spec(), false)),
    )
    .expect("the reader");
    let blocks = blocks_of(&mut reader).expect("the blocks");
    let kept = named_of(&blocks, reader.chroms());
    let of = |chrom: &str, from: u64, to: u64| -> Vec<u64> {
        kept.iter()
            .filter(|(of_the_variant, pos)| of_the_variant == chrom && (from..=to).contains(pos))
            .map(|(_, pos)| *pos)
            .collect()
    };
    assert_eq!(kept.len(), 45);
    assert_eq!(of("chr1", 1, 2000).len(), 28);
    assert_eq!(of("chr1", 2001, u64::MAX), vec![4996, 5033, 5070]);
    assert_eq!(of("chr2", 1, 19000), vec![10250]);
    assert_eq!(of("chr2", 19001, u64::MAX).len(), 13);
}

/// The reader offers its regions to what is under it when it is built: a
/// filter of individuals hands the offer on to the source, and a filter of
/// variants does not, since the variants skipped would never reach its
/// counts. The variants and the counts are the same either way.
#[test]
fn the_offer_by_regions_goes_through_the_filter_of_individuals_and_not_a_filter_of_variants() {
    let individuals: Vec<String> = vcf_reader("many.vcf", None).individuals().to_vec();
    let three = individuals.get(..3).expect("three individuals").to_vec();
    let selection = selection_of(&the_bed_of_the_spec(), false);

    let source = SkippingSource::over(vcf_reader("many.vcf", Some(7)));
    let offers = Arc::clone(&source.offers);
    let under = IndividualsReader::new(source, &three).expect("the filter of individuals");
    let mut reader =
        RegionsReader::new(under, RegionFilter::new(selection.clone())).expect("the reader");
    assert_eq!(offers.load(Ordering::SeqCst), 1);
    let blocks = blocks_of(&mut reader).expect("the blocks");
    assert_eq!(named_of(&blocks, reader.chroms()), the_reference("regions"));
    assert_eq!(reader.filtering_stats(), vec![("regions", pair(500, 45))]);

    let source = SkippingSource::over(vcf_reader("many.vcf", Some(7)));
    let offers = Arc::clone(&source.offers);
    let under = FilteredReader::new(
        source,
        VarFilter::new(VarFilteringCriterion::MaxMaf(1.0)).expect("the filter"),
    )
    .expect("the filter of variants");
    let mut reader = RegionsReader::new(under, RegionFilter::new(selection)).expect("the reader");
    assert_eq!(offers.load(Ordering::SeqCst), 0);
    let blocks = blocks_of(&mut reader).expect("the blocks");
    assert_eq!(named_of(&blocks, reader.chroms()), the_reference("regions"));
    assert_eq!(
        reader.filtering_stats(),
        vec![("regions", pair(500, 45)), ("maf", pair(500, 500))]
    );
}

/// A source that gives neither the chromosome nor the position, whatever
/// it is asked, as a `Variants` built from an array of genotypes does.
struct NoPositions {
    block: Option<Block>,
    chroms: ChromTable,
    individuals: Vec<String>,
    header: SourceHeader,
}

impl BlockReader for NoPositions {
    fn next_block(&mut self) -> Result<Option<Block>> {
        Ok(self.block.take())
    }

    fn individuals(&self) -> &[String] {
        &self.individuals
    }

    fn ploidy(&self) -> usize {
        2
    }

    fn chroms(&self) -> &ChromTable {
        &self.chroms
    }

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

#[test]
fn a_source_with_no_positions_by_regions_is_the_error_of_a_field_at_its_first_block() {
    let source = NoPositions {
        block: Some(Block {
            num_vars: 2,
            num_individuals: 1,
            ploidy: 2,
            gts: vec![0, 1, 1, 1],
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
            vcf_text: None,
        }),
        chroms: ChromTable::new(),
        individuals: vec!["ind1".to_owned()],
        header: SourceHeader {
            individuals: vec!["ind1".to_owned()],
            chrom_lengths: Vec::new(),
            vcf_meta_lines: None,
        },
    };
    let mut reader = RegionsReader::new(
        source,
        RegionFilter::new(selection_of(&the_bed_of_the_spec(), false)),
    )
    .expect("the reader");
    reader.set_needs(Needs::ALL);
    let error = reader.next_block().expect_err("no position");
    assert!(
        matches!(error, Error::FieldsNotInTheBlock { fields } if fields == Needs::CHROM_POS),
        "{error:?}"
    );
    assert!(error.to_string().contains("chrom and pos"), "{error}");
    // After an error there is no block.
    assert!(reader.next_block().expect("no error").is_none());
    assert_eq!(reader.filtering_stats(), vec![("regions", pair(0, 0))]);
}

/// A chromosome number the table of the reader has no name for is a reader
/// with a defect, and the block is left as it was.
#[test]
fn a_chromosome_number_with_no_name_by_regions_is_an_error_and_the_block_is_as_it_was() {
    let mut block = Block {
        num_vars: 2,
        num_individuals: 1,
        ploidy: 1,
        gts: vec![0, 1],
        chrom: Some(vec![0, 3]),
        pos: Some(vec![5, 6]),
        id: None,
        alleles: None,
        qual: None,
        vcf_text: None,
    };
    let mut chroms = ChromTable::new();
    chroms.intern("chr1");
    let mut filter = RegionFilter::new(selection_of(b"chr1\t0\t100\n", false));
    let error = filter
        .filter_block(&mut block, &chroms)
        .expect_err("no name");
    assert!(
        matches!(error, Error::RegionFilterChromNameMissing { number: 3 }),
        "{error:?}"
    );
    assert_eq!(block.pos, Some(vec![5, 6]));
    assert_eq!(filter.stats(), pair(0, 0));
}

/// The error of a BED text, which has to be one.
fn the_error_of(bed: &[u8]) -> Error {
    Regions::from_bed(bed).expect_err("a BED that is refused")
}

fn bed_line(line: u64, problem: BedLineProblem) -> impl Fn(&Error) -> bool {
    move |error| matches!(error, Error::BedLine { line: at, problem: of } if *at == line && *of == problem)
}

/// A line of fewer than three columns separated by tabs, with the line
/// counted over the skipped lines before it; a line of spaces is one, and
/// its message says that BED separates the columns by tabs.
#[test]
fn a_line_of_fewer_than_three_columns_by_regions_is_refused_with_its_line() {
    let error = the_error_of(b"track name=x\n# a comment\n\nchr1\t0\t10\nchr1\t20\n");
    assert!(
        bed_line(
            5,
            BedLineProblem::FewerThanThreeColumns {
                columns: 2,
                with_spaces: false
            }
        )(&error),
        "{error:?}"
    );
    let error = the_error_of(b"chr1 0 10\n");
    assert!(
        bed_line(
            1,
            BedLineProblem::FewerThanThreeColumns {
                columns: 1,
                with_spaces: true
            }
        )(&error),
        "{error:?}"
    );
    let message = error.to_string();
    assert!(message.starts_with("line 1 of the BED file: "), "{message}");
    assert!(message.contains("spaces"), "{message}");
    assert!(message.contains("separates by tabs"), "{message}");
    assert!(message.contains("it has 1 column separated"), "{message}");
}

/// A start or an end that is not a whole number of 0 or more written in
/// digits, or that does not fit in 64 bits.
#[test]
fn a_start_or_an_end_that_is_not_a_whole_number_by_regions_is_refused_with_its_line() {
    for (bed, line, found) in [
        (b"chr1\t-1\t10\n".as_slice(), 1, "-1"),
        (b"chr1\t0\t10\nchr1\tx\t10\n".as_slice(), 2, "x"),
        (b"chr1\t+99\t100\n".as_slice(), 1, "+99"),
        (b"chr1\t\t10\n".as_slice(), 1, ""),
        (b"chr1\t1.5\t10\n".as_slice(), 1, "1.5"),
        // `:` is the byte after `9`.
        (b"c\t1:\t200\n".as_slice(), 1, "1:"),
    ] {
        let error = the_error_of(bed);
        let problem = BedLineProblem::StartNotAWholeNumber {
            found: found.to_owned(),
        };
        assert!(bed_line(line, problem)(&error), "{error:?}");
    }
    let error = the_error_of(b"#\nchr1\t0\t1e3\n");
    let problem = BedLineProblem::EndNotAWholeNumber {
        found: "1e3".to_owned(),
    };
    assert!(bed_line(2, problem)(&error), "{error:?}");
    assert_eq!(
        error.to_string(),
        "line 2 of the BED file: its end is `1e3`, and an end is a whole number of 0 or more, \
         written in digits"
    );
    // A start or an end past the largest number of 64 bits is said to be.
    let error = the_error_of(b"chr1\t18446744073709551616\t18446744073709551617\n");
    let problem = BedLineProblem::StartAboveTheLargest {
        found: "18446744073709551616".to_owned(),
    };
    assert!(bed_line(1, problem)(&error), "{error:?}");
    assert_eq!(
        error.to_string(),
        "line 1 of the BED file: its start, 18446744073709551616, is above \
         18446744073709551615, the largest number of 64 bits"
    );
    let error = the_error_of(b"chr1\t0\t99999999999999999999\n");
    let problem = BedLineProblem::EndAboveTheLargest {
        found: "99999999999999999999".to_owned(),
    };
    assert!(bed_line(1, problem)(&error), "{error:?}");
    // The largest `u64` is one.
    let regions =
        Regions::from_bed(b"chr1\t0\t18446744073709551615\n".as_slice()).expect("the largest end");
    assert!(regions.contains("chr1", u64::MAX));
}

/// A start that is not below its end: equal, which BED allows for a point
/// between two bases, and above.
#[test]
fn a_start_not_below_its_end_by_regions_is_refused_with_its_line() {
    let error = the_error_of(b"chr1\t0\t10\nchr1\t5\t5\n");
    let problem = BedLineProblem::StartNotBelowEnd { start: 5, end: 5 };
    assert!(bed_line(2, problem)(&error), "{error:?}");
    let error = the_error_of(b"chr1\t6\t5\n");
    let problem = BedLineProblem::StartNotBelowEnd { start: 6, end: 5 };
    assert!(bed_line(1, problem)(&error), "{error:?}");
}

/// A BED with no region: an empty file, and one of skipped lines alone.
#[test]
fn a_bed_with_no_region_by_regions_is_refused() {
    for bed in [
        b"".as_slice(),
        b"\n\n".as_slice(),
        b"track name=x\nbrowser position chr1\n# chr1\t0\t10\n".as_slice(),
    ] {
        let error = the_error_of(bed);
        assert!(matches!(error, Error::BedWithNoRegion), "{error:?}");
    }
}

/// The lines that are not regions are skipped, a carriage return before a
/// newline is dropped, and the columns after the third are not read.
#[test]
fn the_lines_the_bed_reader_skips_by_regions_are_empty_comments_track_and_browser() {
    let bed = b"track name=x description=\"a b\"\r\n\
                browser position chr1:1-100\n\
                # chr9\t0\t10\n\
                \n\
                \r\n\
                chr1\t0\t10\tgene1\t0\t+\r\n\
                chr2\t4\t5";
    let regions = Regions::from_bed(bed.as_slice()).expect("the regions");
    assert_eq!(regions.num_regions(), 2);
    assert!(regions.contains("chr1", 10));
    assert!(!regions.contains("chr1", 11));
    assert!(regions.contains("chr2", 5));
    assert!(!regions.contains("chr2", 4));
    assert!(!regions.contains("chr9", 5));
}

/// A gzipped BED, of one member and of two, gives the regions of the plain
/// one, and one cut short is an error of the input.
#[test]
fn a_gzipped_bed_by_regions_gives_the_regions_of_the_plain_one() {
    let plain = the_bed_of_the_spec();
    let gzip = |text: &[u8]| -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(text).expect("the compression");
        encoder.finish().expect("the compression")
    };
    let of_the_plain = Regions::from_bed(plain.as_slice()).expect("the plain BED");
    let one_member = gzip(&plain);
    assert_eq!(one_member.get(..2), Some([0x1f, 0x8b].as_slice()));
    let (first, second) = plain.split_at(40);
    let mut two_members = gzip(first);
    two_members.extend(gzip(second));
    for gzipped in [&one_member, &two_members] {
        let regions = Regions::from_bed(gzipped.as_slice()).expect("the gzipped BED");
        assert_eq!(regions, of_the_plain);
        assert_eq!(regions.num_regions(), 5);
    }
    let cut = one_member.get(..one_member.len() - 10).expect("a cut");
    let error = Regions::from_bed(cut).expect_err("a BED cut short");
    assert!(matches!(error, Error::Io(_)), "{error:?}");
}

/// Whether a range of positions holds none the selection keeps, on the
/// regions of the batches of the vars file of `many.vcf` in the spec: with
/// the BED of the spec the fourth batch, chr2 12100 to 15763, holds none,
/// and with `exclude` no batch lies inside one region.
#[test]
fn keeps_none_of_by_regions_says_which_batches_the_vars_file_reader_skips() {
    let inside = selection_of(&the_bed_of_the_spec(), false);
    let outside = selection_of(&the_bed_of_the_spec(), true);
    for (chrom, min_pos, max_pos, none_inside) in [
        ("chr1", 1000, 4663, false),
        ("chr1", 4700, 8363, false),
        ("chr1", 8400, 10213, true),
        ("chr2", 10250, 12063, false),
        ("chr2", 12100, 15763, true),
        ("chr2", 15800, 19463, false),
        // The edges of chr2 19001 to 30000.
        ("chr2", 15800, 19000, true),
        // The largest position of the batch is the first of a region.
        ("chr2", 15800, 19001, false),
        ("chr2", 30001, 40000, true),
        ("chr4", 5, 6, true),
    ] {
        assert_eq!(
            inside.keeps_none_of(chrom, min_pos, max_pos),
            none_inside,
            "{chrom} {min_pos} {max_pos}"
        );
        assert!(
            !outside.keeps_none_of(chrom, min_pos, max_pos),
            "{chrom} {min_pos} {max_pos}"
        );
    }
    assert!(!inside.keeps_none_of("chr3", 5, 6));
    // Inside one region, and to both of its edges.
    assert!(outside.keeps_none_of("chr1", 4991, 5100));
    assert!(outside.keeps_none_of("chr3", 1, 100_000));
    assert!(!outside.keeps_none_of("chr1", 4990, 5100));
    assert!(!outside.keeps_none_of("chr1", 4991, 5101));
    // No position at all.
    assert!(inside.keeps_none_of("chr1", 10, 9));
    assert!(outside.keeps_none_of("chr1", 10, 9));
}

#[test]
fn the_kind_and_the_keeps_of_a_selection_by_regions_follow_exclude() {
    let inside = selection_of(b"chr1\t0\t10\n", false);
    let outside = selection_of(b"chr1\t0\t10\n", true);
    assert_eq!(inside.kind(), "regions");
    assert_eq!(outside.kind(), "excluded_regions");
    assert_eq!(PassStep::Regions(inside.clone()).kind(), "regions");
    assert_eq!(
        PassStep::Regions(outside.clone()).kind(),
        "excluded_regions"
    );
    assert!(inside.keeps("chr1", 10) && !outside.keeps("chr1", 10));
    assert!(!inside.keeps("chr1", 11) && outside.keeps("chr1", 11));
    assert!(!inside.keeps("chr2", 1) && outside.keeps("chr2", 1));
}

/// `chain_of` builds a `RegionsReader` for the step, a `regions` and an
/// `excluded_regions` step stand together, and a second step of a kind is
/// the error of the reader.
#[test]
fn chain_of_by_regions_builds_the_reader_and_refuses_a_second_of_a_kind() {
    let inside = PassStep::Regions(selection_of(&the_bed_of_the_spec(), false));
    // The regions of chr1 up to 1100, which the 45 hold three of: 1000,
    // 1037 and 1074.
    let outside = PassStep::Regions(selection_of(b"chr1\t0\t1100\n", true));

    let mut chain = chain_of(
        Box::new(vcf_reader("many.vcf", Some(7))),
        &[inside.clone(), outside.clone()],
    )
    .expect("the chain");
    chain.set_needs(Needs::GTS);
    let blocks = blocks_of(&mut chain).expect("the blocks");
    assert_eq!(named_of(&blocks, chain.chroms()).len(), 42);
    assert_eq!(
        chain.filtering_stats(),
        vec![
            ("excluded_regions", pair(45, 42)),
            ("regions", pair(500, 45))
        ]
    );

    for steps in [
        [inside.clone(), inside.clone()],
        [outside.clone(), outside.clone()],
    ] {
        let error = chain_of(Box::new(vcf_reader("many.vcf", None)), &steps)
            .err()
            .expect("a second of a kind");
        let kind = steps.first().expect("a step").kind();
        assert!(
            matches!(error, Error::RegionFilterOfAKindThatIsSet { kind: of } if of == kind),
            "{error:?}"
        );
    }
}

#[test]
fn refuse_a_second_filter_of_a_kind_by_regions_refuses_the_kind_that_is_set_alone() {
    let inside = PassStep::Regions(selection_of(b"chr1\t0\t10\n", false));
    let outside = PassStep::Regions(selection_of(b"chr1\t0\t10\n", true));
    let maf = PassStep::VarFilter(VarFilteringCriterion::MaxMaf(0.8));
    let individuals = PassStep::KeepIndividuals(vec!["ind1".to_owned()]);

    let error = refuse_a_second_filter_of_a_kind(&[maf.clone(), inside.clone()], &inside)
        .expect_err("a second of a kind");
    assert!(
        matches!(
            error,
            Error::RegionFilterOfAKindThatIsSet { kind: "regions" }
        ),
        "{error:?}"
    );
    assert!(error.to_string().contains("regions already"), "{error}");
    assert!(!error.names_the_file());
    assert!(refuse_a_second_filter_of_a_kind(std::slice::from_ref(&inside), &outside).is_ok());
    assert!(refuse_a_second_filter_of_a_kind(&[maf.clone(), individuals.clone()], &inside).is_ok());
    assert!(refuse_a_second_filter_of_a_kind(&[inside.clone(), outside.clone()], &maf).is_ok());
    assert!(refuse_a_second_filter_of_a_kind(&[inside, outside], &individuals).is_ok());
}

/// The errors of a BED name the file in Python, as those of a VCF do.
#[test]
fn the_errors_of_a_bed_by_regions_name_the_file() {
    assert!(the_error_of(b"chr1\t5\t5\n").names_the_file());
    assert!(the_error_of(b"").names_the_file());
}

/// The rows read on the threads of rayon and one after another give the
/// same values, over the blocks of `many.vcf`.
#[test]
fn the_rows_read_on_the_threads_and_one_by_one_by_regions_give_the_same_values() {
    let mut reader = vcf_reader("many.vcf", None);
    let blocks = blocks_of(&mut reader).expect("the blocks");
    for exclude in [false, true] {
        let selection = selection_of(&the_bed_of_the_spec(), exclude);
        for block in &blocks {
            let chrom = block.chrom.as_deref().expect("the chromosomes");
            let pos = block.pos.as_deref().expect("the positions");
            let on_the_threads =
                keep_of_the_rows(&selection, chrom, pos, reader.chroms()).expect("the values");
            let one_by_one = keep_of_the_rows_one_by_one(&selection, chrom, pos, reader.chroms())
                .expect("the values");
            assert_eq!(on_the_threads, one_by_one);
        }
    }
}

/// A line is a header of the genome browser when its first word is `track`
/// or `browser`, and a chromosome whose name starts with one of the two is
/// a region, as bcftools 1.24 and plink2 v2.0.0-a.7.7 read it.
#[test]
fn a_chromosome_whose_name_starts_with_track_by_regions_is_a_region() {
    let bed = b"track\tname=x\ntrack\nbrowser position chr1\ntracks1\t0\t5\nbrowsers\t1\t3\n";
    let regions = Regions::from_bed(bed.as_slice()).expect("the regions");
    assert_eq!(regions.num_regions(), 2);
    assert!(regions.contains("tracks1", 5));
    assert!(regions.contains("browsers", 2));
}

/// A line whose chromosome is empty names no chromosome a variant can be
/// on, and plink2 v2.0.0-a.7.7 refuses it too.
#[test]
fn an_empty_chromosome_by_regions_is_refused_with_its_line() {
    let error = the_error_of(b"chr1\t0\t5\n\t0\t2\n");
    assert!(
        bed_line(2, BedLineProblem::EmptyChromosome)(&error),
        "{error:?}"
    );
    assert_eq!(
        error.to_string(),
        "line 2 of the BED file: its chromosome is empty, and a region is of a chromosome \
         that has a name"
    );
}

/// The byte order mark some editors of Windows write at the start of a text
/// is dropped, and is not the start of the first chromosome name; in a
/// gzipped BED it is the start of the text inside.
#[test]
fn a_byte_order_mark_by_regions_is_dropped_from_the_start_of_the_bed() {
    let with_the_mark = b"\xef\xbb\xbfchr1\t0\t5\n";
    let regions = Regions::from_bed(with_the_mark.as_slice()).expect("the regions");
    assert!(regions.contains("chr1", 5));
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(with_the_mark).expect("the compression");
    let gzipped = encoder.finish().expect("the compression");
    let regions = Regions::from_bed(gzipped.as_slice()).expect("the regions");
    assert!(regions.contains("chr1", 5));
    // A mark that is not at the start is part of the name.
    let regions =
        Regions::from_bed(b"chr1\t0\t5\n\xef\xbb\xbfchr2\t0\t5\n".as_slice()).expect("the regions");
    assert!(!regions.contains("chr2", 5));
}

/// The regions of one chromosome, looked up once by its name as bytes,
/// answer what the selection answers for each position of it: the handle
/// a reader that holds CHROM as bytes keeps from one line to the next.
#[test]
fn the_regions_of_one_chromosome_by_regions_answer_as_the_selection_does() {
    for exclude in [false, true] {
        let selection = selection_of(&the_bed_of_the_spec(), exclude);
        for chrom in ["chr1", "chr2", "chr3", "chr4"] {
            let of_the_chrom = selection.regions_of(chrom.as_bytes());
            for pos in [
                1, 2000, 2001, 4990, 4991, 5100, 5101, 10249, 10250, 19001, 30001,
            ] {
                assert_eq!(
                    of_the_chrom.keeps(pos),
                    selection.keeps(chrom, pos),
                    "{chrom} {pos} exclude {exclude}"
                );
            }
            for (min_pos, max_pos) in [(1, 2000), (2001, 4990), (8400, 10213), (12100, 15763)] {
                assert_eq!(
                    of_the_chrom.keeps_none_of(min_pos, max_pos),
                    selection.keeps_none_of(chrom, min_pos, max_pos),
                    "{chrom} {min_pos} {max_pos} exclude {exclude}"
                );
            }
        }
    }
    let inside = selection_of(&the_bed_of_the_spec(), false);
    assert!(inside.regions_of(b"chr1").keeps(2000));
    assert!(!inside.regions_of(b"chr1").keeps(2001));
    assert!(inside.regions_of(b"chr4").keeps_none_of(1, u64::MAX));
    let outside = selection_of(&the_bed_of_the_spec(), true);
    assert!(outside.regions_of(b"chr4").keeps(1));
    assert!(outside.regions_of(b"chr3").keeps_none_of(1, 100_000));
}

/// A generator of numbers for the test below, the linear congruential one
/// of Knuth's MMIX, so that the BEDs it draws are the same at every run.
struct Draws(u64);

impl Draws {
    /// A number from 0 to `below - 1`.
    fn below(&mut self, below: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33).checked_rem(below).expect("a number above 0")
    }
}

/// Over 500 small BEDs drawn at random, of regions that overlap, touch,
/// nest and come in any order on two chromosomes, `contains` is the set of
/// positions the lines cover, one line at a time, and `keeps` and
/// `keeps_none_of` of both sides and of the handle of one chromosome are
/// what that set says of every position and of every range of positions.
#[test]
fn keeps_and_keeps_none_of_by_regions_agree_with_the_positions_of_the_lines() {
    const LAST_POS: u64 = 45;
    let mut draws = Draws(26_092_026);
    for _ in 0..500 {
        let mut bed = String::new();
        let mut covered = std::collections::HashSet::new();
        for _ in 0..=draws.below(6) {
            let chrom = ["c1", "c2"][usize::try_from(draws.below(2)).unwrap()];
            let start = draws.below(40);
            let end = start + 1 + draws.below(8);
            bed.push_str(&format!("{chrom}\t{start}\t{end}\n"));
            for pos in start + 1..=end {
                covered.insert((chrom, pos));
            }
        }
        for exclude in [false, true] {
            let selection = selection_of(bed.as_bytes(), exclude);
            for chrom in ["c1", "c2", "c3"] {
                let of_the_chrom = selection.regions_of(chrom.as_bytes());
                for pos in 0..=LAST_POS {
                    let inside = covered.contains(&(chrom, pos));
                    assert_eq!(
                        selection.regions.contains(chrom, pos),
                        inside,
                        "{bed}{chrom} {pos}"
                    );
                    assert_eq!(
                        selection.keeps(chrom, pos),
                        inside != exclude,
                        "{bed}{chrom} {pos}"
                    );
                    assert_eq!(
                        of_the_chrom.keeps(pos),
                        inside != exclude,
                        "{bed}{chrom} {pos}"
                    );
                }
                for min_pos in 1..=LAST_POS {
                    for max_pos in min_pos..=LAST_POS {
                        let none_kept = (min_pos..=max_pos)
                            .all(|pos| covered.contains(&(chrom, pos)) == exclude);
                        assert_eq!(
                            selection.keeps_none_of(chrom, min_pos, max_pos),
                            none_kept,
                            "{bed}{chrom} {min_pos} {max_pos} exclude {exclude}"
                        );
                    }
                }
            }
        }
    }
}

/// The source of whole blocks passes over some of the blocks of 7 of
/// `many.vcf` on either side, so the test above runs the filter over a
/// source that skipped variants and gave it others to take out. The
/// numbers are those of the 72 blocks of 7 of the file and of
/// `regions.txt`, counted with Python: 63 blocks hold none of the 45, and
/// 6 blocks, 38 variants, only variants of the 45.
#[test]
fn the_source_of_whole_blocks_by_regions_skips_some_blocks_and_gives_others() {
    for (exclude, skipped) in [(false, 441), (true, 38)] {
        let mut source = SkippingSource::of(vcf_reader("many.vcf", Some(7)), TheSkip::WholeBlocks);
        assert!(source.skip_outside(selection_of(&the_bed_of_the_spec(), exclude)));
        let given: usize = blocks_of(&mut source)
            .expect("the blocks")
            .iter()
            .map(|block| block.num_vars)
            .sum();
        assert_eq!(source.num_skipped(), skipped, "exclude {exclude}");
        assert_eq!(
            u64::try_from(given).unwrap() + skipped,
            500,
            "exclude {exclude}"
        );
    }
}
