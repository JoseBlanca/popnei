//! The tests of the filter that keeps variants at random, each with
//! `random_filter` in its name, against "How it is verified" of that filter
//! in `docs/specs/filters.md`: the draws of Java, the worked example, the
//! table of `many.vcf` in two sizes of block, the filter followed by the
//! first 10, the refusals, and the contract of a reader.

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use super::{DEFAULT_RANDOM_FILTER_SEED, RandomFilter, RandomlyFilteredReader, SplitMix64};
use crate::block::{Block, BlockReader, SourceHeader};
use crate::error::{Error, KeepRateAndSeed, Result};
use crate::filters::{
    FilteringStats, PassStep, RegionSelection, Regions, VarFilteringCriterion, chain_of,
    refuse_a_second_filter_of_a_kind, refuse_a_step, refuse_a_step_after_the_first_n,
};
use crate::io::vcf::{VcfOptions, VcfReader};
use crate::variant::{ChromTable, Needs};

/// `many.vcf`, the 500 variants of 50 diploid individuals of
/// `docs/specs/io_vcf.md`, which lives at the root of the repository.
fn many_vcf() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/reference/vcf/many.vcf")
}

/// A reader over `many.vcf` with every variant given, the ones that failed
/// their FILTER too, 500 variants, in blocks of `num_vars_per_block`.
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

/// A block of one diploid individual whose variants are at `positions` of
/// chr1, each with its own genotype, so that a row moved to the place of
/// another shows in the genotypes as in the positions.
fn block_at(positions: &[u64]) -> Block {
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
        vcf_text: None,
    }
}

/// The ten numbers from 0 to 1 that `nextDouble` of
/// `java.util.SplittableRandom` gives from a seed of 42, to six decimals,
/// as `java tests/reference/filters/SplitMix.java` printed them with
/// OpenJDK 26.0.2.1 on 5 October 2026.
const THE_TEN_NUMBERS_FROM_42: [f64; 10] = [
    0.741_565, 0.159_910, 0.278_601, 0.344_191, 0.038_030, 0.868_228, 0.218_405, 0.800_632,
    0.339_931, 0.618_482,
];

/// The variants of the worked example that the filter keeps at a keep rate
/// of 0.5 and a seed of 42: those whose number is below 0.5.
const KEPT_IN_THE_WORKED_EXAMPLE: [u64; 6] = [2, 3, 4, 5, 7, 9];

/// The generator from a seed of 1234567 gives the first five draws of
/// `nextLong` of `java.util.SplittableRandom`, which `java
/// tests/reference/filters/SplitMix.java` printed with OpenJDK 26.0.2.1 on 5
/// October 2026, and which `splitmix64.c` of Sebastiano Vigna gives too.
#[test]
fn random_filter_generator_gives_the_five_draws_of_java_from_1234567() {
    let mut generator = SplitMix64::new(1_234_567);
    let draws: Vec<u64> = (0..5).map(|_| generator.next_draw()).collect();
    assert_eq!(
        draws,
        [
            6_457_827_717_110_365_317,
            3_203_168_211_198_807_973,
            9_817_491_932_198_370_423,
            4_593_380_528_125_082_431,
            16_408_922_859_458_223_821,
        ]
    );
}

/// The numbers from 0 to 1 from a seed of 42 are those of `nextDouble` of
/// Java, within the half of the last of the six decimals it printed.
#[test]
fn random_filter_generator_gives_the_ten_numbers_of_java_from_42() {
    let mut generator = SplitMix64::new(42);
    for (place, expected) in THE_TEN_NUMBERS_FROM_42.iter().enumerate() {
        let number = generator.next_number();
        assert!(
            (number - expected).abs() <= 5e-7,
            "number {place}: {number} and not {expected}"
        );
    }
}

/// The worked example of the spec: ten variants at a keep rate of 0.5 and a
/// seed of 42 keep 2, 3, 4, 5, 7 and 9, and the counts are 10 given and 6
/// kept; the block is compacted in place, its genotypes with its positions.
#[test]
fn random_filter_keeps_variants_2_3_4_5_7_and_9_of_the_worked_example() {
    let mut filter = RandomFilter::new(0.5, 42).expect("the filter");
    let mut block = block_at(&[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);

    filter.filter_block(&mut block).expect("the block");

    assert_eq!(block.pos.as_deref(), Some(&KEPT_IN_THE_WORKED_EXAMPLE[..]));
    assert_eq!(block.num_vars, 6);
    // The rows kept are the 2nd, 3rd, 4th, 5th, 7th and 9th of `block_at`,
    // whose second allele is the row modulo 3.
    assert_eq!(block.gts, [0, 1, 0, 2, 0, 0, 0, 1, 0, 0, 0, 2]);
    assert!(block.check().is_ok());
    assert_eq!(filter.stats(), pair(10, 6));
}

/// The same ten variants given as a block of 3 and a block of 7 keep the
/// same six, with the same counts: the draws go on from one block to the
/// next.
#[test]
fn random_filter_keeps_the_same_six_of_the_worked_example_in_blocks_of_3_and_7() {
    let mut filter = RandomFilter::new(0.5, 42).expect("the filter");
    let mut first = block_at(&[1, 2, 3]);
    let mut second = block_at(&[4, 5, 6, 7, 8, 9, 10]);

    filter.filter_block(&mut first).expect("the first block");
    filter.filter_block(&mut second).expect("the second block");

    assert_eq!(positions_of(&[first, second]), KEPT_IN_THE_WORKED_EXAMPLE);
    assert_eq!(filter.stats(), pair(10, 6));
}

/// A block whose positions are not of its size is refused, is left as it
/// was, adds nothing to the counts and draws no number: the worked example
/// given after it keeps the same six.
#[test]
fn random_filter_of_a_block_that_does_not_pass_check_draws_no_number() {
    let mut filter = RandomFilter::new(0.5, 42).expect("the filter");
    let mut broken = block_at(&[1, 2, 3]);
    broken.pos = Some(vec![1, 2]);

    let error = filter
        .filter_block(&mut broken)
        .expect_err("the broken block");

    assert!(
        matches!(error, Error::BlockArrayOfAnotherSize { array: "pos", .. }),
        "{error}"
    );
    assert_eq!(broken.num_vars, 3);
    assert_eq!(broken.pos, Some(vec![1, 2]));
    assert_eq!(filter.stats(), pair(0, 0));

    let mut block = block_at(&[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
    filter.filter_block(&mut block).expect("the block");
    assert_eq!(block.pos.as_deref(), Some(&KEPT_IN_THE_WORKED_EXAMPLE[..]));
}

/// A block without the genotypes is filtered as one with them: the filter
/// reads no field of a variant.
#[test]
fn random_filter_keeps_the_same_variants_of_a_block_without_the_genotypes() {
    let mut filter = RandomFilter::new(0.5, 42).expect("the filter");
    let mut block = block_at(&[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
    block.gts = Vec::new();

    filter.filter_block(&mut block).expect("the block");

    assert_eq!(block.pos.as_deref(), Some(&KEPT_IN_THE_WORKED_EXAMPLE[..]));
}

/// A keep rate of 1 keeps all ten variants of the worked example and one of
/// 0 none, as the numbers are below 1 and not below 0.
#[test]
fn random_filter_at_a_keep_rate_of_1_keeps_every_variant_and_at_0_none() {
    for (keep_rate, kept) in [(1.0, 10), (0.0, 0)] {
        let mut filter = RandomFilter::new(keep_rate, 42).expect("the filter");
        let mut block = block_at(&[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
        filter.filter_block(&mut block).expect("the block");
        assert_eq!(block.num_vars, kept, "{keep_rate}");
        assert_eq!(filter.stats(), pair(10, u64::try_from(kept).unwrap_or(0)));
    }
}

/// The three rows of the table of the spec, from
/// `tests/reference/filters/random_draws.py` run on 5 October 2026: the
/// keep rate, the seed, how many of the 500 variants of `many.vcf` are kept
/// and the first five of them by their position on chr1.
const THE_TABLE: [(f64, u64, u64, [u64; 5]); 3] = [
    (0.1, 42, 45, [1148, 1666, 1777, 1888, 2332]),
    (0.5, 42, 243, [1037, 1074, 1111, 1148, 1222]),
    (0.1, 7, 49, [1037, 1962, 2147, 2332, 2591]),
];

/// Each row of the table over `many.vcf`, read with every variant given, in
/// blocks of 7 variants and in blocks of the size popnei chooses: the size
/// of the blocks changes no draw. The first row's counts are 500 given and
/// 45 kept.
#[test]
fn random_filter_over_many_vcf_keeps_the_variants_of_the_table() {
    for (keep_rate, seed, num_kept, first_five) in THE_TABLE {
        for num_vars_per_block in [Some(7), None] {
            let filter = RandomFilter::new(keep_rate, seed).expect("the filter");
            let mut reader =
                RandomlyFilteredReader::new(many_vcf_reader(num_vars_per_block), filter)
                    .expect("the reader");
            let blocks = blocks_of(&mut reader).expect("the blocks");
            let positions = positions_of(&blocks);
            let case = format!("{keep_rate} {seed} {num_vars_per_block:?}");
            assert_eq!(
                u64::try_from(positions.len()).ok(),
                Some(num_kept),
                "{case}"
            );
            assert_eq!(positions.get(..5), Some(&first_five[..]), "{case}");
            assert!(blocks.iter().all(|block| block.num_vars > 0), "{case}");
            for block in &blocks {
                assert!(block.check().is_ok(), "{case}");
            }
            assert_eq!(
                reader.filtering_stats(),
                vec![("random", pair(500, num_kept))],
                "{case}"
            );
        }
    }
    assert_eq!(THE_TABLE[0].2, 45);
}

/// Two chains built from the same step, as two passes of one `Variants`
/// are, keep the same variants: the generator starts at the seed in each.
#[test]
fn random_filter_keeps_the_same_variants_in_every_pass() {
    let steps = vec![PassStep::Random {
        keep_rate: 0.1,
        seed: 42,
    }];
    let mut first = chain_of(Box::new(many_vcf_reader(Some(7))), &steps).expect("the chain");
    let mut second = chain_of(Box::new(many_vcf_reader(None)), &steps).expect("the chain");
    let of_the_first = positions_of(&blocks_of(&mut *first).expect("the blocks"));
    let of_the_second = positions_of(&blocks_of(&mut *second).expect("the blocks"));
    assert_eq!(of_the_first.len(), 45);
    assert_eq!(of_the_first, of_the_second);
}

/// The filter at 0.5 with a seed of 42 and then the first 10, built by
/// `chain_of`, gives the first ten of the 243 of the table, in blocks of 7
/// and of the size popnei chooses.
#[test]
fn random_filter_at_0_5_and_then_the_first_10_gives_the_ten_of_the_spec() {
    let steps = vec![
        PassStep::Random {
            keep_rate: 0.5,
            seed: 42,
        },
        PassStep::FirstN(10),
    ];
    for num_vars_per_block in [Some(7), None] {
        let mut chain =
            chain_of(Box::new(many_vcf_reader(num_vars_per_block)), &steps).expect("the chain");
        let blocks = blocks_of(&mut *chain).expect("the blocks");
        assert_eq!(
            positions_of(&blocks),
            [1037, 1074, 1111, 1148, 1222, 1296, 1370, 1407, 1555, 1592],
            "{num_vars_per_block:?}"
        );
        let stats = chain.filtering_stats();
        assert_eq!(
            stats.iter().map(|(kind, _)| *kind).collect::<Vec<_>>(),
            ["first_n", "random"]
        );
    }
}

/// A keep rate of NaN, -0.1 and 1.5 is refused by `RandomFilter::new` and
/// by `chain_of`, with an error that names `keep_rate` and the value; 0
/// and 1 are accepted.
#[test]
fn random_filter_of_a_keep_rate_out_of_range_is_refused() {
    for keep_rate in [f64::NAN, -0.1, 1.5] {
        let error = RandomFilter::new(keep_rate, 42).expect_err("the keep rate");
        assert!(
            matches!(error, Error::RandomFilterKeepRateOutOfRange { .. }),
            "{error}"
        );
        assert!(error.to_string().contains("`keep_rate`"), "{error}");
        assert!(
            error.to_string().contains(&format!("{keep_rate:?}")),
            "{error}"
        );

        let error = chain_of(
            Box::new(many_vcf_reader(Some(7))),
            &[PassStep::Random {
                keep_rate,
                seed: 42,
            }],
        )
        .err()
        .expect("the keep rate");
        assert!(
            matches!(error, Error::RandomFilterKeepRateOutOfRange { .. }),
            "{error}"
        );
    }
    for keep_rate in [0.0, 1.0] {
        assert!(RandomFilter::new(keep_rate, 42).is_ok(), "{keep_rate}");
    }
}

/// A second filter that keeps variants at random is refused, by
/// `refuse_a_second_filter_of_a_kind` and `refuse_a_step` with the keep
/// rate and the seed of both, and by `chain_of` and
/// `RandomlyFilteredReader::new` over a chain that holds one, with those of
/// the second alone.
#[test]
fn random_filter_second_is_refused_with_the_keep_rates_and_the_seeds() {
    let set = vec![
        PassStep::VarFilter(VarFilteringCriterion::MaxMaf(0.9)),
        PassStep::Random {
            keep_rate: 0.1,
            seed: 42,
        },
    ];
    let second = PassStep::Random {
        keep_rate: 0.25,
        seed: 7,
    };
    let expected = Error::RandomFilterThatIsSet {
        keep_rate: 0.25,
        seed: 7,
        keep_rate_and_seed_that_is_set: Some(KeepRateAndSeed {
            keep_rate: 0.1,
            seed: 42,
        }),
    };
    for error in [
        refuse_a_second_filter_of_a_kind(&set, &second).expect_err("a second one"),
        refuse_a_step(&set, &second).expect_err("a second one"),
    ] {
        assert_eq!(error.to_string(), expected.to_string());
        for said in ["random", "0.1", "42", "0.25", "7"] {
            assert!(error.to_string().contains(said), "{said}: {error}");
        }
    }
    assert!(refuse_a_second_filter_of_a_kind(set.get(..1).unwrap_or_default(), &second).is_ok());

    let mut steps = set.clone();
    steps.push(second);
    let error = chain_of(Box::new(many_vcf_reader(Some(7))), &steps)
        .err()
        .expect("a second one");
    assert_eq!(error.to_string(), expected.to_string());

    let held = RandomlyFilteredReader::new(
        many_vcf_reader(Some(7)),
        RandomFilter::new(0.1, 42).expect("the filter"),
    )
    .expect("the first");
    let error = RandomlyFilteredReader::new(held, RandomFilter::new(0.25, 7).expect("the filter"))
        .expect_err("a second one");
    assert_eq!(
        error.to_string(),
        Error::RandomFilterThatIsSet {
            keep_rate: 0.25,
            seed: 7,
            keep_rate_and_seed_that_is_set: None,
        }
        .to_string()
    );
}

/// The filter that keeps variants at random takes variants out, so after a
/// filter of the first n it is refused, by `refuse_a_step_after_the_first_n`,
/// `refuse_a_step` and `chain_of`, with its kind; before it, it is
/// accepted.
#[test]
fn random_filter_after_the_first_n_is_refused() {
    let random = PassStep::Random {
        keep_rate: 0.5,
        seed: 42,
    };
    let set = vec![PassStep::FirstN(10)];
    let refused = Error::StepAfterTheFirstN { kind: "random" }.to_string();
    for error in [
        refuse_a_step_after_the_first_n(&set, &random).expect_err("after the first n"),
        refuse_a_step(&set, &random).expect_err("after the first n"),
        chain_of(
            Box::new(many_vcf_reader(Some(7))),
            &[PassStep::FirstN(10), random.clone()],
        )
        .err()
        .expect("after the first n"),
    ] {
        assert_eq!(error.to_string(), refused);
    }
    assert!(refuse_a_step(std::slice::from_ref(&random), &PassStep::FirstN(10)).is_ok());
}

/// The kind of the step is the name its counts have in Python, and the
/// seed a step gets when the user gives none is 42.
#[test]
fn random_filter_kind_is_random_and_its_default_seed_42() {
    let step = PassStep::Random {
        keep_rate: 0.5,
        seed: DEFAULT_RANDOM_FILTER_SEED,
    };
    assert_eq!(step.kind(), "random");
    assert_eq!(DEFAULT_RANDOM_FILTER_SEED, 42);
    let filter = RandomFilter::new(0.25, 7).expect("the filter");
    assert_eq!((filter.keep_rate(), filter.seed()), (0.25, 7));
}

/// A source that gives `blocks` and then no more, and records what a
/// reader over it asks of it: the calls to `next_block`, the fields set on
/// it, and the offers of regions, which it takes, as a source that can pass
/// over the variants outside does, saying it passed over 7.
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

/// The reader at 0.5 with a seed of 42 over `blocks`.
fn reader_over(blocks: Vec<Result<Block>>) -> RandomlyFilteredReader<Recording> {
    RandomlyFilteredReader::new(
        Recording::of(blocks),
        RandomFilter::new(0.5, 42).expect("the filter"),
    )
    .expect("the reader")
}

/// A block that the filter empties is not given: over the first variant of
/// the worked example, which is not kept, and then the next nine, the
/// reader gives one block, of the six, having asked its source for both.
#[test]
fn random_filter_does_not_give_a_block_it_emptied() {
    let mut reader = reader_over(vec![
        Ok(block_at(&[1])),
        Ok(block_at(&[2, 3, 4, 5, 6, 7, 8, 9, 10])),
    ]);
    let calls = Arc::clone(&reader.reader.calls);
    let block = reader.next_block().expect("the block").expect("a block");
    assert_eq!(block.pos.as_deref(), Some(&KEPT_IN_THE_WORKED_EXAMPLE[..]));
    assert_eq!(calls.load(Ordering::SeqCst), 2);
    assert!(reader.next_block().expect("no block").is_none());
    assert_eq!(reader.filtering_stats(), vec![("random", pair(10, 6))]);
}

/// A source that gives a block of no variants has a defect: the filter
/// gives the error of it, then nothing, and asks the source once.
#[test]
fn random_filter_over_a_source_that_gives_a_block_of_no_variants_is_the_error_of_a_defect() {
    let mut reader = reader_over(vec![Ok(block_at(&[])), Ok(block_at(&[1, 2]))]);
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
fn random_filter_gives_the_error_of_its_source_and_then_nothing() {
    let mut reader = reader_over(vec![
        Ok(block_at(&[1, 2, 3])),
        Err(Error::ReaderGaveABlockOfNoVariants),
        Ok(block_at(&[4, 5])),
    ]);
    let calls = Arc::clone(&reader.reader.calls);
    assert_eq!(
        reader
            .next_block()
            .expect("the first block")
            .and_then(|block| block.pos),
        Some(vec![2, 3])
    );
    assert!(reader.next_block().is_err());
    assert!(reader.next_block().expect("nothing").is_none());
    assert_eq!(calls.load(Ordering::SeqCst), 2);
}

/// A block of the source that does not pass `check` is the error of the
/// filter, and after it there is no block.
#[test]
fn random_filter_gives_the_error_of_a_block_that_does_not_pass_check_and_then_nothing() {
    let mut broken = block_at(&[1, 2, 3]);
    broken.pos = Some(vec![1, 2]);
    let mut reader = reader_over(vec![Ok(broken), Ok(block_at(&[4, 5]))]);
    let calls = Arc::clone(&reader.reader.calls);

    let error = reader.next_block().expect_err("the broken block");
    assert!(
        matches!(error, Error::BlockArrayOfAnotherSize { .. }),
        "{error}"
    );
    assert!(reader.next_block().expect("nothing").is_none());
    assert_eq!(calls.load(Ordering::SeqCst), 1);
}

/// The filter reads no field of a variant and passes the fields a consumer
/// asks for to its source as they are, the positions alone among them.
#[test]
fn random_filter_passes_the_needs_to_its_source() {
    let mut reader = reader_over(Vec::new());
    let needs = Arc::clone(&reader.reader.needs);

    reader.set_needs(Needs::CHROM_POS);

    assert_eq!(*needs.lock().expect("the needs"), Some(Needs::CHROM_POS));
}

/// The filter answers false to the offer of the regions and does not hand
/// it to its source, and says it skipped none although its source says 7:
/// the variants the source passed over would draw no number.
#[test]
fn random_filter_refuses_the_regions_and_skips_nothing() {
    let mut reader = reader_over(Vec::new());
    let offers = Arc::clone(&reader.reader.offers);
    let regions = Arc::new(Regions::from_bed(&b"chr1\t0\t2000\n"[..]).expect("the regions"));

    assert!(!reader.skip_outside(RegionSelection {
        regions,
        exclude: false,
    }));
    assert_eq!(reader.num_skipped(), 0);
    assert_eq!(offers.load(Ordering::SeqCst), 0);
}
