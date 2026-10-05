//! The tests of the filter of the first n, each with `first_n` in its
//! name, against "How it is verified" of the filter that keeps the first n
//! variants of `docs/specs/filters.md`: the positions bcftools 1.24 gives
//! with `head -n 10`, the counts and the blocks asked of the source, a
//! source that never ends, and how much of a VCF and of a vars file a pass
//! that the filter ended reads.

use std::fs::File;
use std::io::{BufReader, Cursor, Read};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use super::{FirstNReader, first_n_step, refuse_a_step_after_the_first_n, stopped_early};
use crate::block::{Block, BlockReader, SourceHeader};
use crate::error::{Error, Result};
use crate::filters::{
    FilteringStats, PassStep, RegionSelection, Regions, VarFilteringCriterion, chain_of,
    refuse_a_second_filter_of_a_kind, refuse_a_step,
};
use crate::io::vcf::{VcfOptions, VcfReader};
use crate::variant::{ChromTable, Needs};

use VarFilteringCriterion::{MaxLdR2, MaxMaf, MaxMissingRate, MaxObsHet};

/// The first ten variants of `many.vcf`, all on chr1: `bcftools view -H
/// many.vcf | head -n 10` of bcftools 1.24, run on 5 October 2026.
const THE_FIRST_TEN: [u64; 10] = [1000, 1037, 1074, 1111, 1148, 1185, 1222, 1259, 1296, 1333];

/// The first ten variants of `many.vcf` that the MAF filter of 0.8 keeps:
/// `bcftools view -H -Q 0.8:major many.vcf | head -n 10` of bcftools 1.24,
/// run on 5 October 2026.
const THE_FIRST_TEN_AFTER_THE_MAF_FILTER: [u64; 10] =
    [1037, 1074, 1111, 1148, 1222, 1259, 1296, 1333, 1370, 1407];

/// `many.vcf`, the 500 variants of 50 diploid individuals of
/// `docs/specs/io_vcf.md`, which lives at the root of the repository.
fn many_vcf() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/reference/vcf/many.vcf")
}

/// A reader over `many.vcf` with every variant given, the ones that failed
/// their FILTER too, as bcftools reads it, in blocks of
/// `num_vars_per_block`.
fn many_vcf_reader(num_vars_per_block: Option<usize>) -> VcfReader<BufReader<File>> {
    let options = VcfOptions {
        ploidy: 2,
        only_passed: false,
        num_vars_per_block,
    };
    VcfReader::from_path(&many_vcf(), options).expect("the reader of many.vcf")
}

/// Every block a reader gives, until it has no more or it fails.
fn blocks_of(reader: &mut (impl BlockReader + ?Sized)) -> Result<Vec<Block>> {
    let mut blocks = Vec::new();
    while let Some(block) = reader.next_block()? {
        blocks.push(block);
    }
    Ok(blocks)
}

/// The positions of the variants of the blocks, in their order.
fn positions_of(blocks: &[Block]) -> Vec<u64> {
    blocks
        .iter()
        .flat_map(|block| block.pos.clone().unwrap_or_default())
        .collect()
}

/// The two counts of one filter.
fn pair(vars_processed: u64, vars_kept: u64) -> FilteringStats {
    FilteringStats {
        vars_processed,
        vars_kept,
    }
}

/// A reader that gives what its source gives and counts the calls to its
/// `next_block`, which are the blocks asked of the source.
struct Asked<R: BlockReader> {
    reader: R,
    calls: Arc<AtomicUsize>,
}

impl<R: BlockReader> Asked<R> {
    fn over(reader: R) -> Asked<R> {
        Asked {
            reader,
            calls: Arc::new(AtomicUsize::new(0)),
        }
    }
}

impl<R: BlockReader> BlockReader for Asked<R> {
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

/// The most blocks a test may ask of [`NeverEnds`]. A pass that the filter
/// ends asks for 3 at most, and one that it does not end would read and keep
/// blocks until the memory of the machine ran out, which is what the tests of
/// the filter did on 5 October 2026 while it was still a stub that gave every
/// block on: the bound makes such a pass fail at once.
const MOST_BLOCKS_OF_A_SOURCE_THAT_NEVER_ENDS: u64 = 100;

/// A source that never ends: a block of 7 variants of two diploid
/// individuals at every call, on chr1, at positions that rise by 10. It
/// panics when it is asked for more than
/// [`MOST_BLOCKS_OF_A_SOURCE_THAT_NEVER_ENDS`] blocks, so that a pass that
/// does not end fails its test and does not fill the memory.
struct NeverEnds {
    blocks_given: u64,
    next_pos: u64,
    chroms: ChromTable,
    individuals: Vec<String>,
    header: SourceHeader,
}

impl NeverEnds {
    fn new() -> NeverEnds {
        let mut chroms = ChromTable::new();
        chroms.intern("chr1");
        let individuals = vec!["ind1".to_owned(), "ind2".to_owned()];
        NeverEnds {
            blocks_given: 0,
            next_pos: 10,
            chroms,
            header: SourceHeader {
                individuals: individuals.clone(),
                chrom_lengths: Vec::new(),
                vcf_meta_lines: None,
            },
            individuals,
        }
    }
}

impl BlockReader for NeverEnds {
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "a test reads a few blocks of it, at positions far below the largest u64"
    )]
    fn next_block(&mut self) -> Result<Option<Block>> {
        assert!(
            self.blocks_given < MOST_BLOCKS_OF_A_SOURCE_THAT_NEVER_ENDS,
            "a source that never ends was asked for more than \
             {MOST_BLOCKS_OF_A_SOURCE_THAT_NEVER_ENDS} blocks: the pass did not end"
        );
        self.blocks_given += 1;
        let pos: Vec<u64> = (0..7_u64).map(|row| self.next_pos + 10 * row).collect();
        self.next_pos += 70;
        Ok(Some(Block {
            num_vars: 7,
            num_individuals: 2,
            ploidy: 2,
            gts: [0, 1, 1, 1].repeat(7),
            chrom: Some(vec![0; 7]),
            pos: Some(pos),
            id: None,
            alleles: None,
            qual: None,
            vcf_text: None,
        }))
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

/// The steps of a pass with the filter of the first `num_vars` alone.
fn first_n(num_vars: u64) -> Vec<PassStep> {
    vec![PassStep::FirstN(num_vars)]
}

/// The first ten of `many.vcf` are the ten positions of bcftools, in blocks
/// of 7 and in blocks of the size popnei chooses; in blocks of 7 the filter
/// is given the 14 variants of two blocks, keeps 10, and the source is
/// asked for those two blocks and no third.
#[test]
fn the_first_n_of_many_vcf_are_the_first_ten_of_bcftools() {
    for num_vars_per_block in [Some(7), None] {
        let source = Asked::over(many_vcf_reader(num_vars_per_block));
        let calls = Arc::clone(&source.calls);
        let mut reader = FirstNReader::new(source, 10).expect("the filter");
        let blocks = blocks_of(&mut reader).expect("the blocks");
        assert_eq!(
            positions_of(&blocks),
            THE_FIRST_TEN,
            "{num_vars_per_block:?}"
        );
        for block in &blocks {
            assert!(block.check().is_ok(), "{num_vars_per_block:?}");
        }
        // A call after the last gives nothing and asks the source nothing.
        assert!(reader.next_block().expect("no block").is_none());
        assert_eq!(
            calls.load(Ordering::SeqCst),
            num_vars_per_block.map_or(1, |_| 2)
        );
        let stats = reader.filtering_stats();
        if num_vars_per_block.is_some() {
            assert_eq!(stats, vec![("first_n", pair(14, 10))]);
        }
        assert!(
            stopped_early(&first_n(10), &stats),
            "{num_vars_per_block:?}"
        );
    }
}

/// The MAF filter of 0.8 before the first 10, built by `chain_of`, gives
/// the ten positions of bcftools with `-Q 0.8:major`, in blocks of 7 and of
/// the size popnei chooses. In blocks of 7 the two filters count the same
/// blocks: the MAF filter the variants of the blocks read, and the first n
/// the ones the MAF filter kept of them.
#[test]
fn the_first_n_after_the_maf_filter_are_the_first_ten_of_bcftools() {
    let steps = vec![PassStep::VarFilter(MaxMaf(0.8)), PassStep::FirstN(10)];
    for num_vars_per_block in [Some(7), None] {
        let mut chain =
            chain_of(Box::new(many_vcf_reader(num_vars_per_block)), &steps).expect("the chain");
        let blocks = blocks_of(&mut *chain).expect("the blocks");
        assert_eq!(
            positions_of(&blocks),
            THE_FIRST_TEN_AFTER_THE_MAF_FILTER,
            "{num_vars_per_block:?}"
        );
        let stats = chain.filtering_stats();
        assert_eq!(stats.len(), 2);
        let (first_n_kind, first_n_counts) = stats[0];
        let (maf_kind, maf_counts) = stats[1];
        assert_eq!((first_n_kind, maf_kind), ("first_n", "maf"));
        assert_eq!(first_n_counts.vars_kept, 10);
        assert_eq!(first_n_counts.vars_processed, maf_counts.vars_kept);
        if num_vars_per_block.is_some() {
            // The first 14 variants of the file are in two blocks of 7, of
            // which the MAF filter keeps 12, all but 1000 and 1185, as
            // `bcftools view -H -Q 0.8:major many.vcf | head -n 14` of
            // bcftools 1.24 gives, run on 5 October 2026: 5 of the first
            // block, too few, and 7 of the second, of which the first n
            // keeps 5.
            assert_eq!(maf_counts, pair(14, 12));
            assert_eq!(first_n_counts, pair(12, 10));
        }
        assert!(stopped_early(&steps, &stats));
    }
}

/// With n of 14 over blocks of 7 the 14th variant ends the second block,
/// and the source is asked for those two blocks and no third.
#[test]
fn a_first_n_whose_last_variant_ends_a_block_asks_for_no_more_blocks() {
    let source = Asked::over(many_vcf_reader(Some(7)));
    let calls = Arc::clone(&source.calls);
    let mut reader = FirstNReader::new(source, 14).expect("the filter");
    let blocks = blocks_of(&mut reader).expect("the blocks");
    assert_eq!(
        blocks
            .iter()
            .map(|block| block.num_vars)
            .collect::<Vec<_>>(),
        [7, 7]
    );
    assert!(reader.next_block().expect("no block").is_none());
    assert_eq!(calls.load(Ordering::SeqCst), 2);
    assert_eq!(reader.filtering_stats(), vec![("first_n", pair(14, 14))]);
    assert!(stopped_early(&first_n(14), &reader.filtering_stats()));
}

/// The first 10 over a source that never ends, read on one thread, returns
/// with 10 variants, having asked the source for two blocks.
#[test]
fn the_first_n_over_a_source_that_never_ends_returns_on_one_thread() {
    let source = Asked::over(NeverEnds::new());
    let calls = Arc::clone(&source.calls);
    let mut chain = chain_of(Box::new(source), &first_n(10)).expect("the chain");
    let blocks = blocks_of(&mut *chain).expect("the blocks");
    assert_eq!(blocks.iter().map(|block| block.num_vars).sum::<usize>(), 10);
    assert_eq!(calls.load(Ordering::SeqCst), 2);
    assert_eq!(chain.filtering_stats(), vec![("first_n", pair(14, 10))]);
    assert!(stopped_early(&first_n(10), &chain.filtering_stats()));
}

/// A pass over a source that never ends that no filter ends fails at the
/// bound of the source, and does not read on until the memory runs out: the
/// guard that keeps the three tests around this one from doing so when the
/// filter is broken.
#[test]
#[should_panic(expected = "the pass did not end")]
fn a_pass_over_a_source_that_never_ends_that_nothing_ends_fails_at_its_bound() {
    let mut source = NeverEnds::new();
    // The read panics at the bound of the source before it returns, so there
    // is no result to look at.
    let _ = blocks_of(&mut source);
}

/// The first 10 over a source that never ends, with the chain read on the
/// thread of the reader one block ahead, as a binding crate reads a pass:
/// the thread gets `None` from the filter, sends that word and ends, and
/// the source was asked for two blocks.
#[cfg(not(target_family = "wasm"))]
#[test]
fn the_first_n_over_a_source_that_never_ends_returns_through_the_reader_one_block_ahead() {
    use crate::block::with_one_block_ahead;

    let source = Asked::over(NeverEnds::new());
    let calls = Arc::clone(&source.calls);
    let mut chain = chain_of(Box::new(source), &first_n(10)).expect("the chain");
    let (num_vars, stats_of_the_handle) = with_one_block_ahead(&mut chain, |ahead| {
        let blocks = blocks_of(ahead)?;
        Ok((
            blocks.iter().map(|block| block.num_vars).sum::<usize>(),
            ahead.filtering_stats(),
        ))
    })
    .expect("the pass");
    assert_eq!(num_vars, 10);
    assert_eq!(calls.load(Ordering::SeqCst), 2);
    assert_eq!(stats_of_the_handle, vec![("first_n", pair(14, 10))]);
    assert_eq!(chain.filtering_stats(), vec![("first_n", pair(14, 10))]);
    assert!(stopped_early(&first_n(10), &chain.filtering_stats()));
}

/// The first 10 put over the reader one block ahead of a source that never
/// ends, which a consumer that wraps what it is given does: the pass
/// returns with 10 variants, and the reading thread, which had built the
/// third block before it was asked for it, ends when the handle is dropped
/// with that block unasked for.
#[cfg(not(target_family = "wasm"))]
#[test]
fn the_first_n_over_the_reader_one_block_ahead_of_a_source_that_never_ends_returns() {
    use crate::block::with_one_block_ahead;

    let mut source = Asked::over(NeverEnds::new());
    let calls = Arc::clone(&source.calls);
    let (num_vars, stats) = with_one_block_ahead(&mut source, |ahead| {
        let mut reader = FirstNReader::new(ahead, 10)?;
        let blocks = blocks_of(&mut reader)?;
        Ok((
            blocks.iter().map(|block| block.num_vars).sum::<usize>(),
            reader.filtering_stats(),
        ))
    })
    .expect("the pass");
    assert_eq!(num_vars, 10);
    assert_eq!(stats, vec![("first_n", pair(14, 10))]);
    assert!(stopped_early(&first_n(10), &stats));
    // The thread builds the next block as soon as the one before was taken,
    // so it had the third when the pass returned, and asked for no fourth.
    assert_eq!(calls.load(Ordering::SeqCst), 3);
}

/// A `Read` over bytes in memory that counts the bytes it gave.
struct CountedBytes {
    bytes: Cursor<Vec<u8>>,
    given: Arc<AtomicU64>,
}

impl Read for CountedBytes {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        let read = self.bytes.read(buffer)?;
        self.given.fetch_add(
            u64::try_from(read).expect("a count of bytes"),
            Ordering::SeqCst,
        );
        Ok(read)
    }
}

/// A VCF of 100000 variants of 10 diploid individuals on chr1, one every 10
/// base pairs, with genotypes that change from one variant to the next.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "positions up to 1000000 and indices below 100010"
)]
fn the_vcf_of_100000_variants() -> Vec<u8> {
    let individuals: Vec<String> = (1..=10).map(|number| format!("ind{number}")).collect();
    let mut text = String::from("##fileformat=VCFv4.2\n");
    text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
    text.push_str(&individuals.join("\t"));
    text.push('\n');
    let genotypes = ["0/0", "0/1", "1/1", "1/0", "./."];
    for variant in 0..100_000_usize {
        text.push_str(&format!(
            "chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT",
            10 * (variant + 1)
        ));
        for individual in 0..10 {
            text.push('\t');
            text.push_str(genotypes[(variant + individual) % genotypes.len()]);
        }
        text.push('\n');
    }
    text.into_bytes()
}

/// A reader over the VCF of 100000 variants in blocks of 100, through a
/// `Read` that counts the bytes it gives.
fn counted_vcf_reader(vcf: Vec<u8>) -> (VcfReader<BufReader<CountedBytes>>, Arc<AtomicU64>) {
    let given = Arc::new(AtomicU64::new(0));
    let source = CountedBytes {
        bytes: Cursor::new(vcf),
        given: Arc::clone(&given),
    };
    let options = VcfOptions {
        ploidy: 2,
        only_passed: false,
        num_vars_per_block: Some(100),
    };
    let reader = VcfReader::new(BufReader::new(source), options).expect("the reader");
    (reader, given)
}

/// A pass with the first 100 over the VCF of 100000 variants, read in
/// blocks of 100, reads less than a tenth of the file: the reader reads no
/// further than the block that held the 100th variant.
#[test]
fn the_first_n_of_a_vcf_of_100000_variants_reads_less_than_a_tenth_of_it() {
    let vcf = the_vcf_of_100000_variants();
    let size = u64::try_from(vcf.len()).expect("the size");
    let (reader, given) = counted_vcf_reader(vcf);
    let mut chain = chain_of(Box::new(reader), &first_n(100)).expect("the chain");
    let blocks = blocks_of(&mut *chain).expect("the blocks");
    assert_eq!(
        blocks.iter().map(|block| block.num_vars).sum::<usize>(),
        100
    );
    let read = given.load(Ordering::SeqCst);
    assert!(read < size / 10, "{read} bytes read of {size}");
    assert_eq!(chain.filtering_stats(), vec![("first_n", pair(100, 100))]);
    assert!(stopped_early(&first_n(100), &chain.filtering_stats()));
}

/// A pass with the first 100 over the vars file of the 100000 variants in
/// batches of 100, read in a pool of rayon of one thread, decompresses the
/// first batch alone, by the list of the batches the reader read.
#[cfg(not(target_family = "wasm"))]
#[test]
fn the_first_n_of_a_vars_file_in_batches_of_100_decompresses_one_batch() {
    use crate::io::vars::{VarsReader, write_vars};

    let (vcf, _) = counted_vcf_reader(the_vcf_of_100000_variants());
    let (written, num_written) =
        write_vars(vcf, Cursor::new(Vec::new()), Some(100)).expect("the vars file");
    assert_eq!(num_written, 100_000);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .expect("the pool");
    pool.install(|| {
        let reader = VarsReader::new(Cursor::new(written.into_inner())).expect("the reader");
        let mut filtered = FirstNReader::new(reader, 100).expect("the filter");
        let blocks = blocks_of(&mut filtered).expect("the blocks");
        assert_eq!(
            blocks.iter().map(|block| block.num_vars).sum::<usize>(),
            100
        );
        assert_eq!(filtered.reader.batches_read(), [0].as_slice());
        assert_eq!(
            filtered.filtering_stats(),
            vec![("first_n", pair(100, 100))]
        );
        assert!(stopped_early(&first_n(100), &filtered.filtering_stats()));
    });
}

/// A pass given fewer than n reads the whole source, keeps every variant,
/// and did not stop early; given exactly n, it did, since the filter does
/// not read on to find whether more came.
#[test]
fn the_first_n_of_more_than_the_source_holds_keeps_it_all_and_did_not_stop_early() {
    let mut reader = FirstNReader::new(many_vcf_reader(Some(7)), 1000).expect("the filter");
    let blocks = blocks_of(&mut reader).expect("the blocks");
    assert_eq!(
        blocks.iter().map(|block| block.num_vars).sum::<usize>(),
        500
    );
    assert_eq!(reader.filtering_stats(), vec![("first_n", pair(500, 500))]);
    assert!(!stopped_early(&first_n(1000), &reader.filtering_stats()));

    let mut reader = FirstNReader::new(many_vcf_reader(Some(7)), 500).expect("the filter");
    let blocks = blocks_of(&mut reader).expect("the blocks");
    assert_eq!(
        blocks.iter().map(|block| block.num_vars).sum::<usize>(),
        500
    );
    assert_eq!(reader.filtering_stats(), vec![("first_n", pair(500, 500))]);
    assert!(stopped_early(&first_n(500), &reader.filtering_stats()));
}

/// `stopped_early` reads the n of the steps and the counts under
/// `first_n`: a pass with no filter of the first n did not stop early,
/// whatever the counts of its other filters, and the n is that of the
/// steps the pass was built from.
#[test]
fn stopped_early_is_of_the_first_n_of_the_steps_and_its_counts() {
    let maf = vec![PassStep::VarFilter(MaxMaf(0.8))];
    assert!(!stopped_early(&maf, &[("maf", pair(10, 10))]));
    assert!(!stopped_early(&[], &[]));
    let counts = [("first_n", pair(14, 10)), ("maf", pair(14, 12))];
    let steps = vec![PassStep::VarFilter(MaxMaf(0.8)), PassStep::FirstN(10)];
    assert!(stopped_early(&steps, &counts));
    assert!(!stopped_early(&first_n(11), &counts));
    assert!(!stopped_early(&maf, &counts));
}

/// A `num_vars` of 0 is refused by the reader and by the chain, and no
/// block is read.
#[test]
fn a_first_n_of_0_is_refused() {
    let source = Asked::over(many_vcf_reader(Some(7)));
    let calls = Arc::clone(&source.calls);
    assert!(matches!(
        FirstNReader::new(source, 0),
        Err(Error::FirstNOfNoVariants)
    ));
    assert!(matches!(
        chain_of(Box::new(many_vcf_reader(Some(7))), &first_n(0)),
        Err(Error::FirstNOfNoVariants)
    ));
    assert_eq!(calls.load(Ordering::SeqCst), 0);
}

/// The step of a `num_vars` of 1 or more is the filter of that many, and
/// one of 0 is refused with a message that names `num_vars`, which is what
/// a user wrote.
#[test]
fn first_n_step_of_0_is_refused_with_a_message_that_names_num_vars() {
    assert_eq!(first_n_step(10).expect("the step"), PassStep::FirstN(10));
    assert_eq!(first_n_step(1).expect("the step"), PassStep::FirstN(1));
    let error = first_n_step(0).expect_err("a num_vars of 0");
    assert!(matches!(error, Error::FirstNOfNoVariants), "{error}");
    assert!(error.to_string().contains("`num_vars` is 0"), "{error}");
}

/// A second filter of the first n is refused with its kind: by
/// `refuse_a_second_filter_of_a_kind` at the steps, by the reader over a
/// chain that holds one, and by the chain.
#[test]
fn a_second_first_n_is_refused_with_its_kind() {
    let refused = Error::FilterOfAKindThatIsSet { kind: "first_n" };
    let error = refuse_a_second_filter_of_a_kind(&first_n(10), &PassStep::FirstN(20))
        .expect_err("a second filter of the first n");
    assert_eq!(error.to_string(), refused.to_string());
    assert!(refuse_a_second_filter_of_a_kind(&[], &PassStep::FirstN(20)).is_ok());

    let under = FirstNReader::new(many_vcf_reader(Some(7)), 10).expect("the first filter");
    let error = FirstNReader::new(under, 20).expect_err("a second filter of the first n");
    assert_eq!(error.to_string(), refused.to_string());

    let error = chain_of(
        Box::new(many_vcf_reader(Some(7))),
        &[PassStep::FirstN(10), PassStep::FirstN(20)],
    )
    .err()
    .expect("a second filter of the first n");
    assert_eq!(error.to_string(), refused.to_string());
}

/// The six steps that take variants out, one of each kind, which are what
/// may not come after a filter of the first n.
fn the_steps_that_take_variants_out() -> Vec<PassStep> {
    let regions = Arc::new(Regions::from_bed(&b"chr1\t0\t2000\n"[..]).expect("the regions"));
    vec![
        PassStep::VarFilter(MaxMissingRate(0.1)),
        PassStep::VarFilter(MaxMaf(0.8)),
        PassStep::VarFilter(MaxObsHet(0.5)),
        PassStep::VarFilter(MaxLdR2 {
            max_allowed_r2: 0.5,
            max_dist: 1000,
        }),
        PassStep::Regions(RegionSelection {
            regions: Arc::clone(&regions),
            exclude: false,
        }),
        PassStep::Regions(RegionSelection {
            regions,
            exclude: true,
        }),
    ]
}

/// Each step that takes variants out, after a filter of the first n, is
/// refused with its kind by `refuse_a_step_after_the_first_n` and by
/// `chain_of`; before it, each is accepted.
#[test]
fn each_step_that_takes_variants_out_after_the_first_n_is_refused() {
    let set = vec![PassStep::VarFilter(MaxMaf(0.9)), PassStep::FirstN(10)];
    for step in the_steps_that_take_variants_out() {
        let kind = step.kind();
        let refused = Error::StepAfterTheFirstN { kind };
        let error = refuse_a_step_after_the_first_n(&set, &step).expect_err(kind);
        assert_eq!(error.to_string(), refused.to_string());
        assert!(error.to_string().contains(kind), "{error}");
        assert!(error.to_string().contains("first_n"), "{error}");

        let error = chain_of(
            Box::new(many_vcf_reader(Some(7))),
            &[PassStep::FirstN(10), step.clone()],
        )
        .err()
        .expect(kind);
        assert_eq!(error.to_string(), refused.to_string());

        assert!(
            refuse_a_step_after_the_first_n(&[], &step).is_ok(),
            "{kind}"
        );
        assert!(
            refuse_a_step_after_the_first_n(&[PassStep::VarFilter(MaxObsHet(0.9))], &step).is_ok(),
            "{kind}"
        );
        assert!(
            chain_of(
                Box::new(many_vcf_reader(Some(7))),
                &[step.clone(), PassStep::FirstN(10)],
            )
            .is_ok(),
            "{kind}"
        );
    }
}

/// A MAF filter of 0.8 after a MAF filter of 0.9 and a filter of the first
/// n breaks both rules, and is refused as a second filter of its kind, with
/// both thresholds, by `refuse_a_step` and by `chain_of`; a MAF filter after
/// the first n alone is refused as a step after it, and a filter of another
/// kind before it is accepted.
#[test]
fn a_second_maf_filter_after_the_first_n_is_refused_as_a_second_of_its_kind() {
    let set = vec![PassStep::VarFilter(MaxMaf(0.9)), PassStep::FirstN(10)];
    let second = PassStep::VarFilter(MaxMaf(0.8));

    let error = refuse_a_step(&set, &second).expect_err("a second maf filter");
    assert!(
        matches!(
            error,
            Error::VarFilterOfAKindThatIsSet {
                kind: "maf",
                threshold_that_is_set: Some(_),
                ..
            }
        ),
        "{error}"
    );
    assert!(error.to_string().contains("0.9"), "{error}");
    assert!(error.to_string().contains("0.8"), "{error}");

    let mut steps = set.clone();
    steps.push(second.clone());
    let error = chain_of(Box::new(many_vcf_reader(Some(7))), &steps)
        .err()
        .expect("a second maf filter");
    assert!(
        matches!(error, Error::VarFilterOfAKindThatIsSet { kind: "maf", .. }),
        "{error}"
    );

    let error = refuse_a_step(&first_n(10), &second).expect_err("a maf filter after it");
    assert_eq!(
        error.to_string(),
        Error::StepAfterTheFirstN { kind: "maf" }.to_string()
    );
    assert!(refuse_a_step(&[PassStep::VarFilter(MaxObsHet(0.5))], &second).is_ok());
}

/// The filter of individuals takes no variant out and is accepted after a
/// filter of the first n, by `refuse_a_step_after_the_first_n` and by
/// `chain_of`, whose pass gives the first 10 of the two individuals named.
#[test]
fn the_filter_of_individuals_after_the_first_n_is_accepted() {
    let reader = many_vcf_reader(Some(7));
    let names: Vec<String> = reader.individuals().iter().take(2).cloned().collect();
    let individuals = PassStep::KeepIndividuals(names);
    assert!(refuse_a_step_after_the_first_n(&first_n(10), &individuals).is_ok());
    assert!(refuse_a_step_after_the_first_n(&first_n(10), &PassStep::FirstN(20)).is_ok());

    let steps = vec![PassStep::FirstN(10), individuals];
    let mut chain = chain_of(Box::new(reader), &steps).expect("the chain");
    let blocks = blocks_of(&mut *chain).expect("the blocks");
    assert_eq!(positions_of(&blocks), THE_FIRST_TEN);
    assert!(blocks.iter().all(|block| block.num_individuals == 2));
    assert!(stopped_early(&steps, &chain.filtering_stats()));
}

/// The kind of the step is the name its counts have in Python.
#[test]
fn the_kind_of_the_first_n_is_first_n() {
    assert_eq!(PassStep::FirstN(10).kind(), "first_n");
}

/// An error of the source reaches the consumer, and after it the filter
/// gives nothing and does not ask the source again.
#[test]
fn the_first_n_gives_the_error_of_its_source_and_then_nothing() {
    let vcf = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\nchr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\nchr1\tten\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n";
    let options = VcfOptions {
        ploidy: 2,
        only_passed: false,
        num_vars_per_block: Some(1),
    };
    let source = Asked::over(VcfReader::new(&vcf[..], options).expect("the reader"));
    let calls = Arc::clone(&source.calls);
    let mut reader = FirstNReader::new(source, 10).expect("the filter");
    assert_eq!(
        reader
            .next_block()
            .expect("the first block")
            .map(|block| block.num_vars),
        Some(1)
    );
    assert!(reader.next_block().is_err());
    assert!(reader.next_block().expect("nothing").is_none());
    assert_eq!(calls.load(Ordering::SeqCst), 2);
}
