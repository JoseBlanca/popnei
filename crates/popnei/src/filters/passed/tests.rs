//! The tests of the filter of the variants that passed their FILTER,
//! against "How it is verified" of that filter in `docs/specs/filters.md`:
//! the 475 variants of `many.vcf` that bcftools 1.24 keeps and their counts
//! in two sizes of block, the MAF filter after it, the variants against
//! those of a VCF reader that drops the failed lines itself, the same over
//! a vars file, a source with no record of whether its variants passed,
//! `cases.vcf`, the refusals, and the contract of a reader.

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use super::PassedReader;
use crate::block::{Block, BlockReader, SourceHeader};
use crate::error::{Error, Result};
use crate::filters::{
    FilteringStats, PassStep, RegionSelection, Regions, VarFilteringCriterion, chain_of,
    refuse_a_second_filter_of_a_kind, refuse_a_step, refuse_a_step_after_the_first_n,
};
use crate::io::vcf::{VcfOptions, VcfReader};
use crate::variant::{ChromTable, Needs};

/// A file of `tests/reference/vcf`, at the root of the repository.
fn reference_vcf(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/reference/vcf")
        .join(name)
}

/// A reader over the VCF `name` of 2 diploid individuals or more, in blocks
/// of `num_vars_per_block`, with every variant given when `only_passed` is
/// false.
fn vcf_reader(
    name: &str,
    only_passed: bool,
    num_vars_per_block: Option<usize>,
) -> VcfReader<BufReader<File>> {
    let options = VcfOptions {
        ploidy: 2,
        only_passed,
        num_vars_per_block,
    };
    VcfReader::from_path(&reference_vcf(name), options).expect("the reader of the VCF")
}

/// `many.vcf` with every variant given, the 25 whose FILTER is `q10` too,
/// 500 variants, in blocks of `num_vars_per_block`.
fn many_vcf_reader(num_vars_per_block: Option<usize>) -> VcfReader<BufReader<File>> {
    vcf_reader("many.vcf", false, num_vars_per_block)
}

/// Every block a reader gives, until it has no more or it fails.
fn blocks_of(reader: &mut (impl BlockReader + ?Sized)) -> Result<Vec<Block>> {
    let mut blocks = Vec::new();
    while let Some(block) = reader.next_block()? {
        blocks.push(block);
    }
    Ok(blocks)
}

/// The chromosome, by its name, and the position of each variant of the
/// blocks, in their order.
fn places_of(blocks: &[Block], chroms: &ChromTable) -> Vec<(String, u64)> {
    let mut places = Vec::new();
    for block in blocks {
        let chrom = block.chrom.as_ref().expect("the chromosomes");
        let pos = block.pos.as_ref().expect("the positions");
        for (chrom, pos) in chrom.iter().zip(pos) {
            let name = chroms.name(*chrom).expect("the name of the chromosome");
            places.push((name.to_owned(), *pos));
        }
    }
    places
}

/// The two counts of one filter.
fn pair(vars_processed: u64, vars_kept: u64) -> FilteringStats {
    FilteringStats {
        vars_processed,
        vars_kept,
    }
}

/// The variants of `many.bcftools.tsv`, what bcftools 1.24 printed for
/// `many.vcf`, whose FILTER is `PASS` or a dot, by their chromosome and
/// position: those that `bcftools view -H -f .,PASS many.vcf` gives.
fn the_variants_that_passed_by_bcftools() -> Vec<(String, u64)> {
    let text =
        std::fs::read_to_string(reference_vcf("many.bcftools.tsv")).expect("many.bcftools.tsv");
    text.lines()
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let filter = *fields.get(6).expect("the FILTER of the line");
            if filter != "PASS" && filter != "." {
                return None;
            }
            let chrom = (*fields.first().expect("the chromosome")).to_owned();
            let pos = fields
                .get(1)
                .expect("the position")
                .parse()
                .expect("a position");
            Some((chrom, pos))
        })
        .collect()
}

/// The first ten variants of `many.vcf` that passed, all on chr1, as `bcftools
/// view -H -f .,PASS many.vcf` gave them on 6 October 2026; 1259, which is
/// between the seventh and the eighth, is the first that failed.
const THE_FIRST_TEN_THAT_PASSED: [u64; 10] =
    [1000, 1037, 1074, 1111, 1148, 1185, 1222, 1296, 1333, 1370];

/// Over `many.vcf` read with every variant, in blocks of 7 variants and in
/// blocks of the size popnei chooses, the filter keeps the 475 variants
/// that bcftools keeps, in their order, the first ten those of the spec;
/// it is given 500 and keeps 475; and every block it gives holds a variant
/// and passes `check`.
#[test]
fn keeps_the_475_variants_of_many_vcf_that_bcftools_keeps_in_blocks_of_7_and_of_the_default_size() {
    let expected = the_variants_that_passed_by_bcftools();
    assert_eq!(expected.len(), 475);
    for num_vars_per_block in [Some(7), None] {
        let mut reader =
            PassedReader::new(many_vcf_reader(num_vars_per_block)).expect("the reader");
        let blocks = blocks_of(&mut reader).expect("the blocks");
        let places = places_of(&blocks, reader.chroms());
        assert_eq!(places, expected, "{num_vars_per_block:?}");
        let first_ten: Vec<(String, u64)> = THE_FIRST_TEN_THAT_PASSED
            .iter()
            .map(|pos| ("chr1".to_owned(), *pos))
            .collect();
        assert_eq!(
            places.get(..10),
            Some(&first_ten[..]),
            "{num_vars_per_block:?}"
        );
        assert!(!places.contains(&("chr1".to_owned(), 1259)));
        for block in &blocks {
            assert!(block.num_vars > 0, "{num_vars_per_block:?}");
            assert!(block.check().is_ok(), "{num_vars_per_block:?}");
            let passed = block.passed.as_ref().expect("the column passed");
            assert!(passed.iter().all(|passed| *passed));
        }
        assert_eq!(
            reader.filtering_stats(),
            vec![("passed", pair(500, 475))],
            "{num_vars_per_block:?}"
        );
    }
}

/// With the MAF filter of 0.8 after it, built by `chain_of` as a pass of a
/// binding crate is, the MAF filter is given the 475 that passed and keeps
/// 364, which is what `bcftools view -f .,PASS -Q 0.8:major` gives.
#[test]
fn with_the_maf_filter_of_0_8_after_it_the_maf_filter_is_given_475_and_keeps_364() {
    let steps = vec![
        PassStep::Passed,
        PassStep::VarFilter(VarFilteringCriterion::MaxMaf(0.8)),
    ];
    for num_vars_per_block in [Some(7), None] {
        let mut chain =
            chain_of(Box::new(many_vcf_reader(num_vars_per_block)), &steps).expect("the chain");
        let blocks = blocks_of(&mut *chain).expect("the blocks");
        assert_eq!(
            blocks.iter().map(|block| block.num_vars).sum::<usize>(),
            364,
            "{num_vars_per_block:?}"
        );
        assert_eq!(
            chain.filtering_stats(),
            vec![("maf", pair(475, 364)), ("passed", pair(500, 475))],
            "{num_vars_per_block:?}"
        );
    }
}

/// One variant as the tests compare it: its chromosome by name, its
/// position, its id, the bits of its quality and its genotypes.
type Variant = (String, u64, String, u32, Vec<i8>);

/// The variants of the blocks, each with the fields that [`Variant`] holds.
fn variants_of(blocks: &[Block], chroms: &ChromTable) -> Vec<Variant> {
    let mut variants = Vec::new();
    for block in blocks {
        let alleles_per_var = block
            .num_individuals
            .checked_mul(block.ploidy)
            .expect("the alleles of a variant");
        let chrom = block.chrom.as_ref().expect("the chromosomes");
        let pos = block.pos.as_ref().expect("the positions");
        let id = block.id.as_ref().expect("the ids");
        let qual = block.qual.as_ref().expect("the qualities");
        let gts = block.gts.chunks_exact(alleles_per_var);
        for ((((chrom, pos), id), qual), gts) in chrom.iter().zip(pos).zip(id).zip(qual).zip(gts) {
            let name = chroms.name(*chrom).expect("the name of the chromosome");
            variants.push((
                name.to_owned(),
                *pos,
                id.clone(),
                qual.to_bits(),
                gts.to_vec(),
            ));
        }
    }
    variants
}

/// The variants of all the blocks of the filter, taken together, are those
/// of a VCF reader with `only_passed` true, whose blocks are cut elsewhere,
/// since it drops a failed line before the line is in a block: the same
/// 475, with the same chromosome, position, id, quality and genotypes.
#[test]
fn gives_the_variants_of_a_vcf_reader_with_only_passed_true() {
    for num_vars_per_block in [Some(7), None] {
        let mut filtered =
            PassedReader::new(many_vcf_reader(num_vars_per_block)).expect("the reader");
        filtered.set_needs(Needs::ALL);
        let of_the_filter = variants_of(
            &blocks_of(&mut filtered).expect("the blocks"),
            filtered.chroms(),
        );

        let mut only_passed = vcf_reader("many.vcf", true, num_vars_per_block);
        only_passed.set_needs(Needs::ALL);
        let of_only_passed = variants_of(
            &blocks_of(&mut only_passed).expect("the blocks"),
            only_passed.chroms(),
        );

        assert_eq!(of_the_filter.len(), 475, "{num_vars_per_block:?}");
        assert_eq!(of_the_filter, of_only_passed, "{num_vars_per_block:?}");
    }
}

/// Over `cases.vcf` read with every variant, the variants at 100, 300 and
/// 400 of chr1, as the table of `docs/specs/io_vcf.md` has them: the one at
/// 200, whose FILTER is `q10`, is taken out.
#[test]
fn keeps_the_variants_of_cases_vcf_at_100_300_and_400() {
    let mut reader = PassedReader::new(vcf_reader("cases.vcf", false, None)).expect("the reader");
    let blocks = blocks_of(&mut reader).expect("the blocks");
    let positions: Vec<u64> = places_of(&blocks, reader.chroms())
        .into_iter()
        .map(|(_, pos)| pos)
        .collect();
    assert_eq!(positions, [100, 300, 400]);
    assert_eq!(reader.filtering_stats(), vec![("passed", pair(4, 3))]);
}

/// A reader that gives what its source gives with no `passed` column, as a
/// vars file of 1.0 or 1.1 does, and counts the calls to its `next_block`.
struct WithoutPassed<R: BlockReader> {
    reader: R,
    calls: Arc<AtomicUsize>,
    header: SourceHeader,
}

impl<R: BlockReader> WithoutPassed<R> {
    fn over(reader: R) -> WithoutPassed<R> {
        let header = SourceHeader {
            keeps_passed: false,
            ..reader.header().clone()
        };
        WithoutPassed {
            reader,
            calls: Arc::new(AtomicUsize::new(0)),
            header,
        }
    }
}

impl<R: BlockReader> BlockReader for WithoutPassed<R> {
    fn next_block(&mut self) -> Result<Option<Block>> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        Ok(self.reader.next_block()?.map(|block| Block {
            passed: None,
            ..block
        }))
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
        self.reader.set_needs(needs);
    }
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        self.reader.filtering_stats()
    }
    fn header(&self) -> &SourceHeader {
        &self.header
    }
    fn skip_outside(&mut self, selection: RegionSelection) -> bool {
        self.reader.skip_outside(selection)
    }
    fn num_skipped(&self) -> u64 {
        self.reader.num_skipped()
    }
}

/// The error of a source with no record of whether its variants passed:
/// its words are those of the spec, and it is among the errors that name
/// the file, so that in Python its message starts with the path.
fn assert_is_the_error_of_no_record(error: &Error) {
    assert!(matches!(error, Error::PassedNotRecorded), "{error}");
    assert_eq!(
        error.to_string(),
        "the variants hold no record of whether they passed their FILTER, so the filter of the \
         variants that passed cannot run on them: a vars file holds it from format 1.2, written \
         from a VCF"
    );
    assert!(error.names_the_file());
}

/// A source whose blocks have no `passed` column is refused at its first
/// block, with nothing counted, and after the error the filter gives
/// nothing and does not ask its source again.
#[test]
fn over_a_source_without_the_passed_column_is_refused_at_the_first_block() {
    let source = WithoutPassed::over(many_vcf_reader(Some(7)));
    let calls = Arc::clone(&source.calls);
    let mut reader = PassedReader::new(source).expect("the reader");

    let error = reader
        .next_block()
        .expect_err("a source without the column");

    assert_is_the_error_of_no_record(&error);
    assert_eq!(calls.load(Ordering::SeqCst), 1);
    assert_eq!(reader.filtering_stats(), vec![("passed", pair(0, 0))]);
    assert!(reader.next_block().expect("nothing").is_none());
    assert_eq!(calls.load(Ordering::SeqCst), 1);
}

/// Over a vars file written in the test from `many.vcf` read with every
/// variant, in batches of 7, the filter keeps the same 475 as over the VCF,
/// with the same counts.
#[cfg(not(target_family = "wasm"))]
#[test]
fn over_a_vars_file_of_many_vcf_keeps_the_same_475() {
    use std::io::Cursor;

    use crate::io::vars::{VarsReader, write_vars};

    let (bytes, num_vars) =
        write_vars(many_vcf_reader(None), Vec::new(), Some(7)).expect("the vars file");
    assert_eq!(num_vars, 500);
    let source = VarsReader::new(Cursor::new(bytes)).expect("the reader of the vars file");
    let mut chain = chain_of(Box::new(source), &[PassStep::Passed]).expect("the chain");
    let blocks = blocks_of(&mut *chain).expect("the blocks");
    assert_eq!(
        places_of(&blocks, chain.chroms()),
        the_variants_that_passed_by_bcftools()
    );
    assert_eq!(chain.filtering_stats(), vec![("passed", pair(500, 475))]);
}

/// Over a vars file written in the test from a reader whose blocks have no
/// `passed` column, which is what a file of 1.0 or 1.1 holds, the filter
/// gives the error of a source with no record at the first block, having
/// asked the reader of the file for one block alone.
#[cfg(not(target_family = "wasm"))]
#[test]
fn over_a_vars_file_without_the_passed_column_is_refused_at_the_first_block() {
    use std::io::Cursor;

    use crate::io::vars::{VarsReader, write_vars};

    let (bytes, num_vars) = write_vars(
        WithoutPassed::over(many_vcf_reader(None)),
        Vec::new(),
        Some(7),
    )
    .expect("the vars file");
    assert_eq!(num_vars, 500);
    let file = VarsReader::new(Cursor::new(bytes)).expect("the reader of the vars file");
    assert!(!file.header().keeps_passed);
    let source = Counted::over(file);
    let calls = Arc::clone(&source.calls);
    let mut chain = chain_of(Box::new(source), &[PassStep::Passed]).expect("the chain");

    let error = chain.next_block().expect_err("a file without the column");

    assert_is_the_error_of_no_record(&error);
    assert_eq!(calls.load(Ordering::SeqCst), 1);
    assert_eq!(chain.filtering_stats(), vec![("passed", pair(0, 0))]);
    assert!(chain.next_block().expect("nothing").is_none());
    assert_eq!(calls.load(Ordering::SeqCst), 1);
}

/// A reader that gives what its source gives and counts the calls to its
/// `next_block`.
#[cfg(not(target_family = "wasm"))]
struct Counted<R: BlockReader> {
    reader: R,
    calls: Arc<AtomicUsize>,
}

#[cfg(not(target_family = "wasm"))]
impl<R: BlockReader> Counted<R> {
    fn over(reader: R) -> Counted<R> {
        Counted {
            reader,
            calls: Arc::new(AtomicUsize::new(0)),
        }
    }
}

#[cfg(not(target_family = "wasm"))]
impl<R: BlockReader> BlockReader for Counted<R> {
    fn next_block(&mut self) -> Result<Option<Block>> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        self.reader.next_block()
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
        self.reader.set_needs(needs);
    }
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        self.reader.filtering_stats()
    }
    fn header(&self) -> &SourceHeader {
        self.reader.header()
    }
    fn skip_outside(&mut self, selection: RegionSelection) -> bool {
        self.reader.skip_outside(selection)
    }
    fn num_skipped(&self) -> u64 {
        self.reader.num_skipped()
    }
}

/// A second filter of the variants that passed is refused, by
/// `refuse_a_second_filter_of_a_kind`, `refuse_a_step` and `chain_of`, and
/// by `PassedReader::new` over a chain that holds one, all with the same
/// error, which is of what a user wrote and names no file.
#[test]
fn a_second_filter_of_this_kind_is_refused() {
    let set = vec![
        PassStep::Passed,
        PassStep::VarFilter(VarFilteringCriterion::MaxMaf(0.9)),
    ];
    let mut errors = vec![
        refuse_a_second_filter_of_a_kind(&set, &PassStep::Passed).expect_err("a second one"),
        refuse_a_step(&set, &PassStep::Passed).expect_err("a second one"),
    ];
    let mut steps = set.clone();
    steps.push(PassStep::Passed);
    errors.push(
        chain_of(Box::new(many_vcf_reader(Some(7))), &steps)
            .err()
            .expect("a second one"),
    );
    let held = PassedReader::new(many_vcf_reader(Some(7))).expect("the first");
    errors.push(PassedReader::new(held).expect_err("a second one"));
    for error in errors {
        assert!(matches!(error, Error::PassedFilterThatIsSet), "{error}");
        assert!(error.to_string().contains("passed"), "{error}");
        assert!(!error.names_the_file());
    }
    assert!(
        refuse_a_second_filter_of_a_kind(set.get(1..).unwrap_or_default(), &PassStep::Passed)
            .is_ok()
    );
}

/// The filter takes variants out, so after a filter of the first n it is
/// refused, by `refuse_a_step_after_the_first_n`, `refuse_a_step` and
/// `chain_of`, with its kind; before it, it is accepted.
#[test]
fn after_the_first_n_it_is_refused() {
    let set = vec![PassStep::FirstN(10)];
    let refused = Error::StepAfterTheFirstN { kind: "passed" }.to_string();
    for error in [
        refuse_a_step_after_the_first_n(&set, &PassStep::Passed).expect_err("after the first n"),
        refuse_a_step(&set, &PassStep::Passed).expect_err("after the first n"),
        chain_of(
            Box::new(many_vcf_reader(Some(7))),
            &[PassStep::FirstN(10), PassStep::Passed],
        )
        .err()
        .expect("after the first n"),
    ] {
        assert_eq!(error.to_string(), refused);
    }
    assert!(refuse_a_step(&[PassStep::Passed], &PassStep::FirstN(10)).is_ok());
}

/// The kind of the step is the name its counts have in Python, and a step
/// of any other kind stands with it in either order, but for the filter of
/// the first n, which comes after it alone.
#[test]
fn its_kind_is_passed_and_other_kinds_stand_with_it() {
    assert_eq!(PassStep::Passed.kind(), "passed");
    let others = [
        PassStep::VarFilter(VarFilteringCriterion::MaxMaf(0.9)),
        PassStep::KeepIndividuals(vec!["ind1".to_owned()]),
        PassStep::Random {
            keep_rate: 0.5,
            seed: 42,
        },
    ];
    for other in &others {
        assert!(
            refuse_a_step(&[PassStep::Passed], other).is_ok(),
            "{other:?}"
        );
        assert!(
            refuse_a_step(std::slice::from_ref(other), &PassStep::Passed).is_ok(),
            "{other:?}"
        );
    }
}

/// A block of one diploid individual whose variants are at `positions` of
/// chr1, with `passed` as given, each with its own genotype, so that a row
/// moved to the place of another shows in the genotypes as in the
/// positions.
fn block_at(positions: &[u64], passed: Option<Vec<bool>>) -> Block {
    let gts: Vec<i8> = (0..positions.len())
        .flat_map(|row| [0, i8::try_from(row % 3).expect("a small allele")])
        .collect();
    Block {
        num_vars: positions.len(),
        num_individuals: 1,
        ploidy: 2,
        gts,
        chrom: Some(vec![0; positions.len()]),
        pos: Some(positions.to_vec()),
        id: None,
        alleles: None,
        qual: None,
        passed,
        vcf_text: None,
    }
}

/// A source that gives `blocks` and then no more, and records what a
/// reader over it asks of it: the calls to `next_block`, the fields set on
/// it, and the offers of regions, which it takes, as a source that can pass
/// over the variants outside does.
struct Recording {
    blocks: std::vec::IntoIter<Result<Block>>,
    calls: Arc<AtomicUsize>,
    needs: Arc<Mutex<Option<Needs>>>,
    offers: Arc<AtomicUsize>,
    chroms: ChromTable,
    header: SourceHeader,
}

impl Recording {
    fn of(blocks: Vec<Result<Block>>) -> Recording {
        let mut chroms = ChromTable::new();
        chroms.intern("chr1");
        Recording {
            blocks: blocks.into_iter(),
            calls: Arc::new(AtomicUsize::new(0)),
            needs: Arc::new(Mutex::new(None)),
            offers: Arc::new(AtomicUsize::new(0)),
            chroms,
            header: SourceHeader {
                individuals: vec!["ind1".to_owned()],
                chrom_lengths: Vec::new(),
                vcf_meta_lines: None,
                keeps_passed: true,
            },
        }
    }
}

impl BlockReader for Recording {
    fn next_block(&mut self) -> Result<Option<Block>> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        self.blocks.next().transpose()
    }
    fn individuals(&self) -> &[String] {
        &self.header.individuals
    }
    fn ploidy(&self) -> usize {
        2
    }
    fn chroms(&self) -> &ChromTable {
        &self.chroms
    }
    fn set_needs(&mut self, needs: Needs) {
        *self.needs.lock().expect("the needs") = Some(needs);
    }
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        Vec::new()
    }
    fn header(&self) -> &SourceHeader {
        &self.header
    }
    fn skip_outside(&mut self, _selection: RegionSelection) -> bool {
        self.offers.fetch_add(1, Ordering::SeqCst);
        true
    }
    fn num_skipped(&self) -> u64 {
        7
    }
}

/// The filter over `blocks`.
fn reader_over(blocks: Vec<Result<Block>>) -> PassedReader<Recording> {
    PassedReader::new(Recording::of(blocks)).expect("the reader")
}

/// The rows whose `passed` is true are kept in place, their genotypes with
/// their positions, and a block that the filter empties is not given: over
/// a block of three that all failed and then one of four of which the
/// second and the fourth failed, the reader gives one block, of two, having
/// asked its source for both.
#[test]
fn keeps_the_rows_that_passed_in_place_and_does_not_give_a_block_it_emptied() {
    let mut reader = reader_over(vec![
        Ok(block_at(&[1, 2, 3], Some(vec![false, false, false]))),
        Ok(block_at(
            &[4, 5, 6, 7],
            Some(vec![true, false, true, false]),
        )),
    ]);
    let calls = Arc::clone(&reader.reader.calls);

    let block = reader.next_block().expect("the block").expect("a block");

    assert_eq!(block.pos.as_deref(), Some(&[4, 6][..]));
    // The rows kept are the first and the third of `block_at`, whose second
    // allele is the row modulo 3.
    assert_eq!(block.gts, [0, 0, 0, 2]);
    assert_eq!(block.passed, Some(vec![true, true]));
    assert!(block.check().is_ok());
    assert_eq!(calls.load(Ordering::SeqCst), 2);
    assert!(reader.next_block().expect("no block").is_none());
    assert_eq!(reader.filtering_stats(), vec![("passed", pair(7, 2))]);
}

/// The filter asks its source for the fields of the consumer and for
/// `passed`, and answers no to an offer of regions without handing it on.
#[test]
fn asks_for_passed_and_does_not_hand_on_the_regions() {
    let mut reader = reader_over(Vec::new());
    let needs = Arc::clone(&reader.reader.needs);
    let offers = Arc::clone(&reader.reader.offers);

    reader.set_needs(Needs::GTS);
    assert_eq!(
        *needs.lock().expect("the needs"),
        Some(Needs::GTS.union(Needs::PASSED))
    );

    let regions = Arc::new(Regions::from_bed(&b"chr1\t0\t10\n"[..]).expect("the regions"));
    assert!(!reader.skip_outside(RegionSelection {
        regions,
        exclude: false,
    }));
    assert_eq!(offers.load(Ordering::SeqCst), 0);
    assert_eq!(reader.num_skipped(), 0);
}

/// A block whose `passed` column is not of its size has a defect: the
/// filter gives the error of it, with nothing counted, then nothing, and
/// does not ask its source again.
#[test]
fn a_block_whose_passed_is_not_of_its_size_is_the_error_of_a_defect() {
    let mut reader = reader_over(vec![
        Ok(block_at(&[1, 2, 3], Some(vec![true, false]))),
        Ok(block_at(&[4, 5], Some(vec![true, true]))),
    ]);
    let calls = Arc::clone(&reader.reader.calls);

    let error = reader.next_block().expect_err("the broken block");

    assert!(
        matches!(
            error,
            Error::BlockArrayOfAnotherSize {
                array: "passed",
                ..
            }
        ),
        "{error}"
    );
    assert_eq!(reader.filtering_stats(), vec![("passed", pair(0, 0))]);
    assert!(reader.next_block().expect("nothing").is_none());
    assert_eq!(calls.load(Ordering::SeqCst), 1);
}

/// A source that gives a block of no variants has a defect: the filter
/// gives the error of it, then nothing, and asks the source once.
#[test]
fn over_a_source_that_gives_a_block_of_no_variants_is_the_error_of_a_defect() {
    let mut reader = reader_over(vec![
        Ok(block_at(&[], Some(Vec::new()))),
        Ok(block_at(&[1, 2], Some(vec![true, true]))),
    ]);
    let calls = Arc::clone(&reader.reader.calls);

    let error = reader.next_block().expect_err("the block of no variants");
    assert!(
        matches!(error, Error::ReaderGaveABlockOfNoVariants),
        "{error}"
    );
    assert!(reader.next_block().expect("nothing").is_none());
    assert_eq!(calls.load(Ordering::SeqCst), 1);
}

/// An error of the source reaches the consumer, and after it the filter
/// gives nothing and does not ask the source again.
#[test]
fn gives_the_error_of_its_source_and_then_nothing() {
    let mut reader = reader_over(vec![
        Ok(block_at(&[1, 2, 3], Some(vec![true, false, true]))),
        Err(Error::ReaderGaveABlockOfNoVariants),
        Ok(block_at(&[4, 5], Some(vec![true, true]))),
    ]);
    let calls = Arc::clone(&reader.reader.calls);
    assert_eq!(
        reader
            .next_block()
            .expect("the first block")
            .and_then(|block| block.pos),
        Some(vec![1, 3])
    );
    assert!(reader.next_block().is_err());
    assert!(reader.next_block().expect("nothing").is_none());
    assert_eq!(calls.load(Ordering::SeqCst), 2);
}
