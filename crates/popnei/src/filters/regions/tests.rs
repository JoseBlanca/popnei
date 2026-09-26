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

/// A `regions` and an `excluded_regions` step over a source that skips:
/// the first hands its regions to the source and the second, over it, is
/// refused the offer it makes, so the source is offered once and the counts
/// are those without the skip. A filter by regions that handed the second
/// offer on would have the source pass over the variants of both, and the
/// second would count fewer variants given than the first kept.
#[test]
fn a_second_filter_by_regions_over_a_source_that_skips_is_not_handed_on() {
    let inside = PassStep::Regions(selection_of(&the_bed_of_the_spec(), false));
    let outside = PassStep::Regions(selection_of(b"chr1\t0\t1100\n", true));
    for skips in [TheSkip::EveryVariant, TheSkip::WholeBlocks] {
        let source = SkippingSource::of(vcf_reader("many.vcf", Some(7)), skips);
        let offers = Arc::clone(&source.offers);
        let mut chain =
            chain_of(Box::new(source), &[inside.clone(), outside.clone()]).expect("the chain");
        assert_eq!(offers.load(Ordering::SeqCst), 1, "{skips:?}");
        let blocks = blocks_of(&mut chain).expect("the blocks");
        assert_eq!(named_of(&blocks, chain.chroms()).len(), 42, "{skips:?}");
        assert_eq!(
            chain.filtering_stats(),
            vec![
                ("excluded_regions", pair(45, 42)),
                ("regions", pair(500, 45))
            ],
            "{skips:?}"
        );
    }
}

/// A source of the tests that gives the blocks it was built with over a
/// table of chromosome names it was given, and counts how many times it
/// was asked for a block.
struct GivenBlocks {
    left: Vec<Block>,
    chroms: ChromTable,
    individuals: Vec<String>,
    header: SourceHeader,
    calls: Arc<AtomicUsize>,
}

impl GivenBlocks {
    /// A source of one haploid individual whose chromosomes are named
    /// `names`, numbered in that order.
    fn of(blocks: Vec<Block>, names: &[&str]) -> GivenBlocks {
        let mut chroms = ChromTable::new();
        for name in names {
            chroms.intern(name);
        }
        let mut left = blocks;
        left.reverse();
        GivenBlocks {
            left,
            chroms,
            individuals: vec!["ind1".to_owned()],
            header: SourceHeader {
                individuals: vec!["ind1".to_owned()],
                chrom_lengths: Vec::new(),
                vcf_meta_lines: None,
            },
            calls: Arc::new(AtomicUsize::new(0)),
        }
    }
}

impl BlockReader for GivenBlocks {
    fn next_block(&mut self) -> Result<Option<Block>> {
        self.calls.fetch_add(1, Ordering::SeqCst);
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

/// A block of one haploid individual of the variants given, each its
/// chromosome number and its position, whose genotype is the number of the
/// variant, so that a test sees which variants stayed.
fn haploid_block(variants: &[(u32, u64)]) -> Block {
    Block {
        num_vars: variants.len(),
        num_individuals: 1,
        ploidy: 1,
        gts: (0..variants.len())
            .map(|var| i8::try_from(var % 100).unwrap())
            .collect(),
        chrom: Some(variants.iter().map(|(chrom, _)| *chrom).collect()),
        pos: Some(variants.iter().map(|(_, pos)| *pos).collect()),
        id: None,
        alleles: None,
        qual: None,
        vcf_text: None,
    }
}

/// A source that gives a block of no variants has a defect: the reader
/// gives the error of it and does not ask the source again.
#[test]
fn a_source_that_gives_a_block_of_no_variants_by_regions_is_the_error_of_a_defect() {
    let source = GivenBlocks::of(
        vec![haploid_block(&[]), haploid_block(&[(0, 5)])],
        &["chr1"],
    );
    let calls = Arc::clone(&source.calls);
    let mut reader = RegionsReader::new(
        source,
        RegionFilter::new(selection_of(b"chr1\t0\t10\n", false)),
    )
    .expect("the reader");
    let error = reader.next_block().expect_err("a block of no variants");
    assert!(
        matches!(error, Error::ReaderGaveABlockOfNoVariants),
        "{error:?}"
    );
    assert!(
        reader
            .next_block()
            .expect("no block after the error")
            .is_none()
    );
    assert_eq!(calls.load(Ordering::SeqCst), 1);
}

/// An error of the filter ends the pass: no block after it, and the source
/// is not asked again.
#[test]
fn after_an_error_of_the_filter_by_regions_there_is_no_block_and_the_source_is_not_asked_again() {
    let source = GivenBlocks::of(
        vec![
            haploid_block(&[(0, 5)]),
            haploid_block(&[(0, 6), (1, 7)]),
            haploid_block(&[(0, 8)]),
        ],
        &["chr1"],
    );
    let calls = Arc::clone(&source.calls);
    let mut reader = RegionsReader::new(
        source,
        RegionFilter::new(selection_of(b"chr1\t0\t10\n", false)),
    )
    .expect("the reader");
    assert!(reader.next_block().expect("the first block").is_some());
    let error = reader.next_block().expect_err("a chromosome with no name");
    assert!(
        matches!(error, Error::RegionFilterChromNameMissing { number: 1 }),
        "{error:?}"
    );
    assert!(
        reader
            .next_block()
            .expect("no block after the error")
            .is_none()
    );
    assert_eq!(calls.load(Ordering::SeqCst), 2);
    assert_eq!(reader.filtering_stats(), vec![("regions", pair(1, 1))]);
}

/// Of two chromosome numbers with no name in a block large enough to be
/// shared out over the threads, the error is the one of the first row,
/// wherever the threads found an error.
#[test]
fn the_error_of_a_block_by_regions_is_the_one_of_its_first_row_that_has_one() {
    let mut variants: Vec<(u32, u64)> = (1..=200_000).map(|pos| (0, pos)).collect();
    // The first at the end of the first half of the rows and the second
    // at the start of the second, where a thread that takes the second
    // half finds it at once.
    *variants.get_mut(99_999).unwrap() = (9, 1);
    *variants.get_mut(100_000).unwrap() = (8, 1);
    let mut chroms = ChromTable::new();
    chroms.intern("chr1");
    for _ in 0..5 {
        let mut block = haploid_block(&variants);
        let mut filter = RegionFilter::new(selection_of(b"chr1\t0\t10\n", false));
        let error = filter
            .filter_block(&mut block, &chroms)
            .expect_err("no name");
        assert!(
            matches!(error, Error::RegionFilterChromNameMissing { number: 9 }),
            "{error:?}"
        );
    }
}

/// The chromosome numbers of a source need not follow the order of the
/// names in the BED, nor its chromosomes all be in it: chr2 comes first
/// here, then one the BED does not name, then chr1.
#[test]
fn chromosomes_in_any_order_by_regions_are_each_looked_up_by_name() {
    let mut chroms = ChromTable::new();
    for name in ["chr2", "chrUn", "chr1"] {
        chroms.intern(name);
    }
    let bed = b"chr1\t0\t10\nchr2\t100\t110\n";
    for (exclude, expected) in [
        (false, vec![(0, 105), (2, 5), (0, 101)]),
        (true, vec![(0, 5), (1, 5), (1, 105), (2, 105)]),
    ] {
        let mut block = haploid_block(&[
            (0, 5),
            (0, 105),
            (1, 5),
            (1, 105),
            (2, 5),
            (2, 105),
            (0, 101),
        ]);
        let mut filter = RegionFilter::new(selection_of(bed, exclude));
        filter
            .filter_block(&mut block, &chroms)
            .expect("the filter");
        let kept: Vec<(u32, u64)> = block
            .chrom
            .unwrap()
            .into_iter()
            .zip(block.pos.unwrap())
            .collect();
        assert_eq!(kept, expected, "exclude {exclude}");
    }
}

/// A chromosome name of a BED that is not UTF-8 text is read, as the spec
/// says, and matches no chromosome a source names, which are text.
#[test]
fn a_chromosome_name_that_is_not_utf8_by_regions_is_read_and_matches_no_variant() {
    let regions =
        Regions::from_bed(b"chr\xff1\t0\t10\nchr1\t0\t5\n".as_slice()).expect("the regions");
    assert_eq!(regions.num_regions(), 2);
    assert!(regions.contains("chr1", 5));
    assert!(!regions.contains("chr1", 6));
    assert!(!regions.contains("chr\u{fffd}1", 5));
}

/// A BED that names its chromosomes `1` where `many.vcf` names them `chr1`
/// keeps nothing, and the consumer gives the error of a pass that gave no
/// variant, with the counts of the filter, 500 given and 0 kept.
#[test]
fn a_bed_of_other_names_by_regions_is_the_error_of_a_pass_that_gave_no_variant() {
    let mut reader = RegionsReader::new(
        vcf_reader("many.vcf", None),
        RegionFilter::new(selection_of(b"1\t0\t2000\n2\t0\t30000\n", false)),
    )
    .expect("the reader");
    let error = match crate::kinship::calc_kinship(&mut reader, None, false) {
        Ok(_) => panic!("a kinship of no variant"),
        Err(error) => error,
    };
    assert!(
        matches!(&error, Error::PassGaveNoVariant { num_vars_of_the_source: 500, filters }
            if *filters == vec![("regions", pair(500, 0))]),
        "{error:?}"
    );
    assert!(
        error
            .to_string()
            .contains("the `regions` filter was given 500 and kept 0"),
        "{error}"
    );
}

/// A reader over a source that refuses the offer of the regions for it, and
/// hands on everything else: the source then gives every variant, which is
/// the pass without the skip.
struct RefusesTheOffer<R: BlockReader>(R);

impl<R: BlockReader> BlockReader for RefusesTheOffer<R> {
    fn next_block(&mut self) -> Result<Option<Block>> {
        self.0.next_block()
    }

    fn individuals(&self) -> &[String] {
        self.0.individuals()
    }

    fn ploidy(&self) -> usize {
        self.0.ploidy()
    }

    fn chroms(&self) -> &ChromTable {
        self.0.chroms()
    }

    fn set_needs(&mut self, needs: Needs) {
        self.0.set_needs(needs);
    }

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        self.0.filtering_stats()
    }

    fn header(&self) -> &SourceHeader {
        self.0.header()
    }

    fn skip_outside(&mut self, _selection: RegionSelection) -> bool {
        false
    }

    fn num_skipped(&self) -> u64 {
        0
    }
}

/// The bytes of the vars file of `many.vcf`, every variant of it, in
/// batches of `num_vars_per_batch`.
fn vars_file_of_many_vcf(num_vars_per_batch: usize) -> Vec<u8> {
    let (bytes, num_vars) = crate::io::vars::write_vars(
        vcf_reader("many.vcf", None),
        Vec::new(),
        Some(num_vars_per_batch),
    )
    .expect("the vars file");
    assert_eq!(num_vars, 500);
    bytes
}

fn vars_reader(bytes: &[u8]) -> crate::io::vars::VarsReader<std::io::Cursor<Vec<u8>>> {
    crate::io::vars::VarsReader::new(std::io::Cursor::new(bytes.to_vec())).expect("the reader")
}

/// The sources of `many.vcf` a test of the skip runs over: the VCF in
/// blocks of 7 and of the default size, and its vars file in batches of 100
/// and of 7, each named for the message of an assertion.
fn the_sources_of_many_vcf() -> Vec<(String, Box<dyn BlockReader>)> {
    let of_100 = vars_file_of_many_vcf(100);
    let of_7 = vars_file_of_many_vcf(7);
    vec![
        (
            "the VCF in blocks of 7".to_owned(),
            Box::new(vcf_reader("many.vcf", Some(7))),
        ),
        ("the VCF".to_owned(), Box::new(vcf_reader("many.vcf", None))),
        (
            "the vars file in batches of 100".to_owned(),
            Box::new(vars_reader(&of_100)),
        ),
        (
            "the vars file in batches of 7".to_owned(),
            Box::new(vars_reader(&of_7)),
        ),
    ]
}

/// Deliverable 1 of work package 5: the 45 and the 455 of bcftools and
/// plink2, and the counts of 500 given, from the VCF and from its vars file
/// with the skip of the source and without it, the variants compared by the
/// names of their chromosomes and their positions. Both sources take the
/// offer.
#[test]
fn skip_outside_gives_the_45_and_the_455_with_the_skip_and_without() {
    for (exclude, name, num_kept) in [(false, "regions", 45_u64), (true, "excluded_regions", 455)] {
        let expected = the_reference(name);
        for skips in [true, false] {
            for (source_name, mut source) in the_sources_of_many_vcf() {
                let what = format!("{source_name}, exclude {exclude}, skips {skips}");
                let selection = selection_of(&the_bed_of_the_spec(), exclude);
                if skips {
                    assert!(source.skip_outside(selection.clone()), "{what}");
                    // The offer taken, the one of the filter is taken too.
                } else {
                    source = Box::new(RefusesTheOffer(source));
                }
                let mut reader =
                    RegionsReader::new(source, RegionFilter::new(selection)).expect("the reader");
                reader.set_needs(Needs::GTS);
                let blocks = blocks_of(&mut reader).expect("the blocks");
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

/// The VCF reader passes over the 455 lines outside the regions, and the
/// 45 inside with `exclude`, and gives the rest with their genotypes, their
/// columns and the text of their lines.
#[test]
fn skip_outside_the_vcf_reader_passes_over_the_lines_the_selection_keeps_none_of() {
    for (exclude, skipped) in [(false, 455_u64), (true, 45)] {
        for num_vars_per_block in [Some(7), None] {
            let mut source = vcf_reader("many.vcf", num_vars_per_block);
            source.set_needs(Needs::ALL | Needs::VCF_TEXT);
            assert!(source.skip_outside(selection_of(&the_bed_of_the_spec(), exclude)));
            let blocks = blocks_of(&mut source).expect("the blocks");
            assert_eq!(source.num_skipped(), skipped, "exclude {exclude}");
            let given: usize = blocks.iter().map(|block| block.num_vars).sum();
            assert_eq!(u64::try_from(given).unwrap() + skipped, 500);
            for block in &blocks {
                block.check().expect("a block of its size");
                let text = block.vcf_text.as_ref().expect("the text");
                assert_eq!(text.num_vars(), block.num_vars);
                // The text of each line is that of its variant.
                for (row, pos) in block.pos.as_ref().unwrap().iter().enumerate() {
                    let pos_of_the_line: u64 =
                        text.fixed(row).split('\t').nth(1).unwrap().parse().unwrap();
                    assert_eq!(pos_of_the_line, *pos);
                }
            }
            let expected = the_reference(if exclude {
                "excluded_regions"
            } else {
                "regions"
            });
            assert_eq!(named_of(&blocks, source.chroms()), expected);
        }
    }
}

/// Deliverable 1 of work package 5: the vars file of `many.vcf` in batches
/// of 100 has five, whose regions in the footer are chr1 1000 to 4663; chr1
/// 4700 to 8363; chr1 8400 to 10213 and chr2 10250 to 12063; chr2 12100 to
/// 15763; and chr2 15800 to 19463. The BED of the spec keeps nothing of the
/// fourth alone, which the reader does not read, and with `exclude` no
/// batch lies inside one region.
#[test]
fn skip_outside_the_vars_reader_does_not_read_the_fourth_batch_of_many_vcf() {
    let bytes = vars_file_of_many_vcf(100);
    for (exclude, read, skipped) in [
        (false, vec![0, 1, 2, 4], 100_u64),
        (true, vec![0, 1, 2, 3, 4], 0),
    ] {
        let mut source = vars_reader(&bytes);
        assert!(source.skip_outside(selection_of(&the_bed_of_the_spec(), exclude)));
        let blocks = blocks_of(&mut source).expect("the blocks");
        assert_eq!(source.batches_read(), read.as_slice(), "exclude {exclude}");
        assert_eq!(source.num_skipped(), skipped, "exclude {exclude}");
        let given: usize = blocks.iter().map(|block| block.num_vars).sum();
        assert_eq!(u64::try_from(given).unwrap() + skipped, 500);
    }
    // Without the offer every batch is read.
    let mut source = vars_reader(&bytes);
    blocks_of(&mut source).expect("the blocks");
    assert_eq!(source.batches_read(), [0, 1, 2, 3, 4].as_slice());
    assert_eq!(source.num_skipped(), 0);
}

/// A threshold filter before the filter by regions: the source is not
/// handed the regions, so the threshold filter counts every variant, and
/// the filter by regions counts what it kept.
#[test]
fn skip_outside_is_not_offered_through_a_threshold_filter_before_the_filter() {
    let steps = [
        PassStep::VarFilter(VarFilteringCriterion::MaxMaf(0.8)),
        PassStep::Regions(selection_of(&the_bed_of_the_spec(), false)),
    ];
    for (source_name, source) in the_sources_of_many_vcf() {
        let mut chain = chain_of(source, &steps).expect("the chain");
        chain.set_needs(Needs::CHROM_POS);
        let blocks = blocks_of(&mut chain).expect("the blocks");
        let stats = chain.filtering_stats();
        // 384 of the 500 pass the maf filter at 0.8, which the table of the
        // threshold filters of the spec gives.
        assert_eq!(
            stats.get(1),
            Some(&("maf", pair(500, 384))),
            "{source_name}"
        );
        let (kind, of_the_regions) = stats.first().copied().unwrap();
        assert_eq!(kind, "regions");
        assert_eq!(of_the_regions.vars_processed, 384, "{source_name}");
        let kept = named_of(&blocks, chain.chroms());
        assert_eq!(u64::try_from(kept.len()).unwrap(), of_the_regions.vars_kept);
        // What is kept is of the 45 and passed the maf filter.
        let the_45 = the_reference("regions");
        assert!(
            kept.iter().all(|variant| the_45.contains(variant)),
            "{source_name}"
        );
    }
}

/// A `regions` and an `excluded_regions` step together over the real
/// sources: the first hands its regions to the source, the second looks at
/// every variant it gets, and the counts are those without the skip.
#[test]
fn skip_outside_with_a_regions_and_an_excluded_regions_step_gives_the_counts_of_both() {
    let steps = [
        PassStep::Regions(selection_of(&the_bed_of_the_spec(), false)),
        PassStep::Regions(selection_of(b"chr1\t0\t1100\n", true)),
    ];
    for (source_name, source) in the_sources_of_many_vcf() {
        let mut chain = chain_of(source, &steps).expect("the chain");
        chain.set_needs(Needs::GTS);
        let blocks = blocks_of(&mut chain).expect("the blocks");
        assert_eq!(named_of(&blocks, chain.chroms()).len(), 42, "{source_name}");
        assert_eq!(
            chain.filtering_stats(),
            vec![
                ("excluded_regions", pair(45, 42)),
                ("regions", pair(500, 45))
            ],
            "{source_name}"
        );
    }
}

/// The filter of individuals between the filter and the source hands the
/// regions on: the source passes over the 455, and the pass keeps the 45
/// with the genotypes of the three individuals.
#[test]
fn skip_outside_is_handed_on_by_the_filter_of_individuals() {
    let three: Vec<String> = ["ind05", "ind00", "ind49"].map(str::to_owned).to_vec();
    let mut under = IndividualsReader::new(vcf_reader("many.vcf", Some(7)), &three)
        .expect("the filter of individuals");
    assert!(under.skip_outside(selection_of(&the_bed_of_the_spec(), false)));
    let blocks = blocks_of(&mut under).expect("the blocks");
    assert_eq!(under.num_skipped(), 455);
    assert_eq!(named_of(&blocks, under.chroms()), the_reference("regions"));
    assert!(blocks.iter().all(|block| block.num_individuals == 3));

    let steps = [
        PassStep::KeepIndividuals(three),
        PassStep::Regions(selection_of(&the_bed_of_the_spec(), false)),
    ];
    for (source_name, source) in the_sources_of_many_vcf() {
        let mut chain = chain_of(source, &steps).expect("the chain");
        let blocks = blocks_of(&mut chain).expect("the blocks");
        assert_eq!(
            named_of(&blocks, chain.chroms()),
            the_reference("regions"),
            "{source_name}"
        );
        assert_eq!(
            chain.filtering_stats(),
            vec![("regions", pair(500, 45))],
            "{source_name}"
        );
    }
}

/// `Reblock` and the reader one block ahead hand the offer on to the
/// source, the second to its thread, and the pass gives the 45 with the
/// counts of 500.
#[test]
fn skip_outside_is_handed_on_by_reblock_and_the_reader_one_block_ahead() {
    let reblocked =
        crate::block::Reblock::new(vcf_reader("many.vcf", Some(7)), Some(10)).expect("the reblock");
    let mut reader = RegionsReader::new(
        reblocked,
        RegionFilter::new(selection_of(&the_bed_of_the_spec(), false)),
    )
    .expect("the reader");
    let blocks = blocks_of(&mut reader).expect("the blocks");
    assert_eq!(named_of(&blocks, reader.chroms()), the_reference("regions"));
    assert_eq!(reader.filtering_stats(), vec![("regions", pair(500, 45))]);

    let mut source = vcf_reader("many.vcf", Some(7));
    let (kept, stats, answered) = crate::block::with_one_block_ahead(&mut source, |ahead| {
        let answered = ahead.skip_outside(selection_of(&the_bed_of_the_spec(), true));
        let mut reader = RegionsReader::new(
            ahead,
            RegionFilter::new(selection_of(&the_bed_of_the_spec(), true)),
        )?;
        let blocks = blocks_of(&mut reader)?;
        Ok((
            named_of(&blocks, reader.chroms()),
            reader.filtering_stats(),
            answered,
        ))
    })
    .expect("the pass");
    assert!(answered);
    assert_eq!(kept, the_reference("excluded_regions"));
    assert_eq!(stats, vec![("excluded_regions", pair(500, 455))]);
    // The thread had built the block ahead before the offer reached it,
    // and the offer holds from the next block the source builds, so of the
    // 45 the source passed over those after that block, which the counts
    // hold whatever their number.
    let skipped = source.num_skipped();
    assert!(skipped > 0 && skipped <= 45, "{skipped}");
}

/// A POS that does not parse gives its line a row, whatever the regions
/// are, so the parse gives the error of that column with the skip; a line
/// outside the regions whose genotypes are wrong is passed over and not
/// parsed.
#[test]
fn skip_outside_a_pos_that_does_not_parse_is_the_error_of_that_column() {
    let vcf = "##fileformat=VCFv4.3\n\
               #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\n\
               chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\n\
               chr1\t500\t.\tA\tT\t.\t.\t.\tGT\t0/7\n\
               chr2\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\n\
               chr1\tx5\t.\tA\tT\t.\t.\t.\tGT\t0/1\n";
    let options = VcfOptions {
        ploidy: 2,
        only_passed: false,
        num_vars_per_block: None,
    };
    let mut reader =
        VcfReader::new(std::io::Cursor::new(vcf.as_bytes().to_vec()), options).expect("the reader");
    assert!(reader.skip_outside(selection_of(b"chr1\t0\t10\n", false)));
    let error = reader.next_block().expect_err("a POS that does not parse");
    assert!(
        matches!(&error, Error::VcfDataLine { line: 6, .. }),
        "{error:?}"
    );
    assert!(error.to_string().contains("POS"), "{error}");

    // Without the POS that does not parse, the wrong genotype of chr1 500,
    // outside the regions, is never parsed.
    let vcf = vcf.rsplit_once("chr1\tx5").unwrap().0;
    let mut reader =
        VcfReader::new(std::io::Cursor::new(vcf.as_bytes().to_vec()), options).expect("the reader");
    assert!(reader.skip_outside(selection_of(b"chr1\t0\t10\n", false)));
    let blocks = blocks_of(&mut reader).expect("the blocks");
    assert_eq!(named_of(&blocks, reader.chroms()), vec![named("chr1", 5)]);
    assert_eq!(reader.num_skipped(), 2);
    // chr2 had no line given, so it has no number.
    assert_eq!(reader.chroms().len(), 1);
}

/// The error of the first block of the VCF of `lines` under a header of
/// two individuals, read with the genotypes asked for, with the regions of
/// chr1 1 to 10 handed to the reader or not; `None` when it gives blocks.
fn the_error_of_the_vcf(lines: &[u8], skips: bool) -> Option<String> {
    let mut vcf =
        b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n\
                    chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0\n"
            .to_vec();
    vcf.extend_from_slice(lines);
    let options = VcfOptions {
        ploidy: 2,
        only_passed: false,
        num_vars_per_block: None,
    };
    let mut reader = VcfReader::new(std::io::Cursor::new(vcf), options).expect("the reader");
    reader.set_needs(Needs::GTS | Needs::CHROM_POS);
    if skips {
        assert!(reader.skip_outside(selection_of(b"chr1\t0\t10\n", false)));
    }
    blocks_of(&mut reader).err().map(|error| error.to_string())
}

/// A line outside the regions keeps the checks of its shape: too few
/// columns, nine first columns that are not UTF-8, fewer and more columns
/// of individuals than the header has, a FORMAT with no GT, and a plain VCF
/// cut inside its last line, which is skipped, each give the error of the
/// pass without the skip.
#[test]
fn skip_outside_a_line_outside_the_regions_keeps_the_checks_of_its_shape() {
    let lines: [&[u8]; 7] = [
        b"chr1\t500\t.\tA\tT\n",
        b"chr1\t500\t\xff\tA\tT\t.\t.\t.\tGT\t0/1\t0/0\n",
        b"chr1\t500\t.\tA\tT\t.\t.\t.\tGT\t0/1\n",
        b"chr1\t500\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0\t1/1\n",
        b"chr1\t500\t.\tA\tT\t.\t.\t.\tDP\t3\t4\n",
        b"chr1\t500\t.\tA\tT\t.\t.\t.\tGT\n",
        b"chr1\t500\t.\tA\tT\t.\t.",
    ];
    for line in lines {
        let without = the_error_of_the_vcf(line, false);
        assert!(without.is_some(), "{:?}", String::from_utf8_lossy(line));
        assert_eq!(
            the_error_of_the_vcf(line, true),
            without,
            "{:?}",
            String::from_utf8_lossy(line)
        );
    }
    // A whole line outside the regions is still passed over.
    assert_eq!(
        the_error_of_the_vcf(b"chr1\t500\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/7\n", true),
        None
    );
}
