//! The tests of the three statistics in one pass, each with
//! `variants_summary` in its name, from "The result so far and the three in
//! one pass" of "The Rust interface" of `docs/specs/stats.md` and "How it is
//! verified" of "The three statistics of a file in one pass" of
//! `docs/specs/js_sources.md`.
//!
//! Each pass reads `many.vcf`, 500 variants of 50 diploid individuals on two
//! chromosomes, and each statistic of the summary is compared with what its
//! own function gives over a second reader of the same file, written with
//! `{:?}`: Rust writes a `f64` with the fewest digits that read back as the
//! same number, so two results written the same hold the same numbers to
//! the bit.

use std::fs::File;
use std::io::BufReader;
use std::sync::{Arc, Mutex};

use super::super::fixtures::vcf_reader_of;
use super::{VariantsSummary, VariantsSummaryConfig, calc_variants_summary};
use crate::block::{Block, BlockReader, SourceHeader};
use crate::error::{Error, Result};
use crate::filters::{FilteringStats, PassStep, RegionSelection, VarFilteringCriterion, chain_of};
use crate::io::vcf::VcfReader;
use crate::stats::{
    ExpHet, HistBins, LengthsFrom, Maf, ObsHet, PerVarDistribsConfig, PerVarStat, Pop, Pops, SoFar,
    calc_per_individual_stats, calc_per_var_distribs, calc_var_density,
};
use crate::variant::{ChromTable, Needs};

/// The variants of `many.vcf`.
const NUM_VARS_OF_MANY_VCF: u64 = 500;

/// The size of the blocks of most passes of these tests, which makes five
/// blocks of `many.vcf`; `None` beside it is the size popnei chooses, which
/// for 50 individuals makes one block of the 500.
const NUM_VARS_PER_BLOCK: u64 = 100;

/// The two sizes of the blocks the comparisons with the three functions are
/// made with.
const BLOCK_SIZES: [Option<usize>; 2] = [Some(100), None];

/// The width of the windows of the density of these tests, 1000 base
/// pairs, which gives `many.vcf`, whose last position is 19463, twenty
/// windows on each chromosome.
const WINDOW_SIZE: u64 = 1000;

/// Lengths that make the windows of both chromosomes before the pass, all
/// of them past the last variant of `many.vcf`, which has no length in its
/// header.
fn the_lengths_given() -> Vec<(String, u64)> {
    vec![("chr2".to_owned(), 25_000), ("chr1".to_owned(), 30_000)]
}

/// The six statistics over `popA` and `popB` of `many.vcf`, with numbers
/// that are none of their defaults: 12 called genotypes, a threshold of
/// 0.85 and eight bins.
fn per_var_config_of(reader: &dyn BlockReader) -> PerVarDistribsConfig {
    let individuals = reader.individuals();
    let named = |range: std::ops::Range<usize>| individuals[range].to_vec();
    let pops = Pops::from_names(
        &[
            ("popA".to_owned(), named(0..20)),
            ("popB".to_owned(), named(20..50)),
        ],
        individuals,
    )
    .expect("popA and popB of many.vcf");
    PerVarDistribsConfig {
        stats: vec![
            PerVarStat::ObsHet,
            PerVarStat::Maf,
            PerVarStat::ExpHet,
            PerVarStat::UnbiasedExpHet,
            PerVarStat::PolyVarsRatio,
            PerVarStat::MissingRate,
        ],
        pops,
        bins: HistBins::linear(0.0, 1.0, 8).expect("eight bins"),
        obs_het: ObsHet::new(12),
        maf: Maf::new(2, 12).expect("the maf of diploid variants"),
        exp_het: ExpHet::new(2, 2, 12).expect("the expected heterozygosity of diploid variants"),
        poly_threshold: 0.85,
    }
}

/// Which of the three statistics a summary is asked for.
#[derive(Debug, Clone, Copy)]
struct Asked {
    per_var: bool,
    per_individual: bool,
    density: bool,
}

/// The seven summaries of one statistic at least.
fn every_summary() -> Vec<Asked> {
    (1..8_u8)
        .map(|bits| Asked {
            per_var: bits & 1 != 0,
            per_individual: bits & 2 != 0,
            density: bits & 4 != 0,
        })
        .collect()
}

/// The configuration of a summary of `asked` over `reader`, with the
/// density over `chrom_lengths`.
fn summary_config_of(
    reader: &dyn BlockReader,
    asked: Asked,
    chrom_lengths: Option<&[(String, u64)]>,
) -> VariantsSummaryConfig {
    VariantsSummaryConfig {
        per_var: asked.per_var.then(|| per_var_config_of(reader)),
        per_individual: asked.per_individual,
        density: asked
            .density
            .then(|| (WINDOW_SIZE, chrom_lengths.map(<[_]>::to_vec))),
    }
}

/// The function after a block that does nothing.
fn nothing(_: &dyn SoFar<VariantsSummary>) -> Result<()> {
    Ok(())
}

/// A reader of `many.vcf` with blocks of `num_vars_per_block`, through
/// `steps`.
fn many_vcf_through(num_vars_per_block: Option<usize>, steps: &[PassStep]) -> Box<dyn BlockReader> {
    let source = vcf_reader_of("vcf/many.vcf", num_vars_per_block);
    chain_of(Box::new(source), steps).expect("the chain over many.vcf")
}

/// The filter of the major allele frequency with `threshold`.
fn maf_filter(threshold: f64) -> PassStep {
    PassStep::VarFilter(VarFilteringCriterion::MaxMaf(threshold))
}

/// What each of the three own functions gives over a reader of `steps`,
/// written with `{:?}`, `None` for one that was not asked for.
fn of_their_own_functions(
    num_vars_per_block: Option<usize>,
    steps: &[PassStep],
    asked: Asked,
    chrom_lengths: Option<&[(String, u64)]>,
) -> [Option<String>; 3] {
    let per_var = asked.per_var.then(|| {
        let mut reader = many_vcf_through(num_vars_per_block, steps);
        let config = per_var_config_of(&*reader);
        let result = calc_per_var_distribs(&mut *reader, &config).expect("the distributions");
        format!("{result:?}")
    });
    let per_individual = asked.per_individual.then(|| {
        let mut reader = many_vcf_through(num_vars_per_block, steps);
        let result = calc_per_individual_stats(&mut *reader).expect("the rates");
        format!("{result:?}")
    });
    let density = asked.density.then(|| {
        let mut reader = many_vcf_through(num_vars_per_block, steps);
        let result =
            calc_var_density(&mut *reader, WINDOW_SIZE, chrom_lengths).expect("the density");
        format!("{result:?}")
    });
    [per_var, per_individual, density]
}

/// The three statistics of a summary, written with `{:?}`.
fn written(summary: &VariantsSummary) -> [Option<String>; 3] {
    [
        summary.per_var.as_ref().map(|result| format!("{result:?}")),
        summary
            .per_individual
            .as_ref()
            .map(|result| format!("{result:?}")),
        summary.density.as_ref().map(|result| format!("{result:?}")),
    ]
}

/// The summary over a reader of `steps`, written with `{:?}`.
fn of_the_summary(
    num_vars_per_block: Option<usize>,
    steps: &[PassStep],
    asked: Asked,
    chrom_lengths: Option<&[(String, u64)]>,
) -> [Option<String>; 3] {
    let mut reader = many_vcf_through(num_vars_per_block, steps);
    let config = summary_config_of(&*reader, asked, chrom_lengths);
    let summary = calc_variants_summary(&mut *reader, &config, &mut nothing).expect("the summary");
    written(&summary)
}

/// Every one of the seven summaries gives each statistic it was asked for
/// equal to the bit to its own function with the same options, and `None`
/// for each it was not, so one left out changes neither of the others: over
/// blocks of 100 and of the size popnei chooses, and with the density over
/// no lengths and over lengths given.
#[test]
fn variants_summary_gives_each_statistic_asked_for_as_its_own_function_and_none_of_the_others() {
    let given = the_lengths_given();
    for num_vars_per_block in BLOCK_SIZES {
        for chrom_lengths in [None, Some(&given[..])] {
            for asked in every_summary() {
                assert_eq!(
                    of_the_summary(num_vars_per_block, &[], asked, chrom_lengths),
                    of_their_own_functions(num_vars_per_block, &[], asked, chrom_lengths),
                    "{asked:?}, blocks of {num_vars_per_block:?}, lengths {chrom_lengths:?}"
                );
            }
        }
    }
}

/// A summary of none of the three is refused before the reader is asked
/// for anything.
#[test]
fn variants_summary_of_none_of_the_three_is_refused_before_the_pass() {
    let mut reader = Recording::of_many_vcf();
    let asked = Arc::clone(&reader.asked);
    let config = VariantsSummaryConfig {
        per_var: None,
        per_individual: false,
        density: None,
    };
    let error = calc_variants_summary(&mut reader, &config, &mut nothing)
        .expect_err("a summary of nothing");
    assert!(
        matches!(error, Error::VariantsSummaryOfNoStatistic),
        "{error:?}"
    );
    assert_eq!(*asked.lock().expect("what was asked"), Vec::new());
}

/// What a reader over `many.vcf` was asked by a pass: the fields set on it,
/// in order, and how many blocks it gave.
#[derive(Debug, PartialEq)]
enum AnAsk {
    Needs(Needs),
    Block,
}

/// A reader over `many.vcf` in blocks of 100 that records what it is asked
/// for, which the pass moves to the thread that reads one block ahead, and
/// gives each block as `alter` leaves it.
struct Recording {
    reader: VcfReader<BufReader<File>>,
    asked: Arc<Mutex<Vec<AnAsk>>>,
    alter: fn(&mut Block),
}

impl Recording {
    fn of_many_vcf() -> Recording {
        Recording::of_many_vcf_altered(|_| ())
    }

    fn of_many_vcf_altered(alter: fn(&mut Block)) -> Recording {
        Recording {
            reader: vcf_reader_of("vcf/many.vcf", Some(100)),
            asked: Arc::new(Mutex::new(Vec::new())),
            alter,
        }
    }

    /// How many blocks it has given.
    fn num_blocks_given(&self) -> usize {
        self.asked
            .lock()
            .expect("what was asked")
            .iter()
            .filter(|ask| **ask == AnAsk::Block)
            .count()
    }

    fn record(&self, ask: AnAsk) {
        self.asked.lock().expect("what was asked").push(ask);
    }
}

impl BlockReader for Recording {
    fn next_block(&mut self) -> Result<Option<Block>> {
        let mut block = self.reader.next_block()?;
        if let Some(block) = &mut block {
            self.record(AnAsk::Block);
            (self.alter)(block);
        }
        Ok(block)
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
        self.record(AnAsk::Needs(needs));
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

/// The density alone asks the reader for the chromosome and the position
/// and for no genotype, once, before the first block; with either statistic
/// of the genotypes it asks for the genotypes too, and those alone without
/// the density. Each statistic is still what its own function gives.
#[test]
fn variants_summary_asks_the_reader_for_what_the_statistics_asked_for_read() {
    let of = |per_var, per_individual, density| Asked {
        per_var,
        per_individual,
        density,
    };
    for (asked, needs) in [
        (of(false, false, true), Needs::CHROM_POS),
        (of(true, false, false), Needs::GTS),
        (of(false, true, false), Needs::GTS),
        (of(true, true, false), Needs::GTS),
        (of(true, false, true), Needs::GTS | Needs::CHROM_POS),
        (of(false, true, true), Needs::GTS | Needs::CHROM_POS),
        (of(true, true, true), Needs::GTS | Needs::CHROM_POS),
    ] {
        let mut reader = Recording::of_many_vcf();
        let recorded = Arc::clone(&reader.asked);
        let config = summary_config_of(&reader, asked, None);
        let summary =
            calc_variants_summary(&mut reader, &config, &mut nothing).expect("the summary");
        let mut expected = vec![AnAsk::Needs(needs)];
        expected.extend((0..5).map(|_| AnAsk::Block));
        assert_eq!(
            *recorded.lock().expect("what was asked"),
            expected,
            "{asked:?}"
        );
        assert_eq!(
            written(&summary),
            of_their_own_functions(Some(100), &[], asked, None),
            "{asked:?}"
        );
    }
}

/// A threshold of the major allele frequency that keeps every variant of
/// `many.vcf`, so that the five calls come after 100, 200, 300, 400 and 500
/// variants.
const KEEPS_EVERY_VARIANT: f64 = 1.0;

/// A threshold that leaves some variants of every block out, so that what
/// the pass has read and what its filter was given differ.
const KEEPS_SOME_VARIANTS: f64 = 0.9;

/// What a call of the function after a block was given, with the summary
/// written with `{:?}`.
#[derive(Debug, PartialEq)]
struct ACall {
    num_vars: u64,
    filtering_stats: Vec<(&'static str, FilteringStats)>,
    result: [Option<String>; 3],
}

/// After each block of 100, through the filter of the major allele
/// frequency, the summary of the three so far is the summary over the
/// variants the pass has read, those of a second reader that puts the
/// filter of the first n after the same filter, with the counts of the
/// filter after that block; and the last one is what the pass returns.
#[test]
fn variants_summary_so_far_after_each_block_is_the_summary_over_the_variants_read() {
    let given = the_lengths_given();
    let all_three = Asked {
        per_var: true,
        per_individual: true,
        density: true,
    };
    for chrom_lengths in [None, Some(&given[..])] {
        for threshold in [KEEPS_EVERY_VARIANT, KEEPS_SOME_VARIANTS] {
            let steps = [maf_filter(threshold)];
            let mut reader = many_vcf_through(Some(100), &steps);
            let config = summary_config_of(&*reader, all_three, chrom_lengths);
            let mut calls = Vec::new();
            let mut keep = |so_far: &dyn SoFar<VariantsSummary>| {
                calls.push(ACall {
                    num_vars: so_far.num_vars(),
                    filtering_stats: so_far.filtering_stats(),
                    result: written(&so_far.result()?),
                });
                Ok(())
            };
            let returned =
                calc_variants_summary(&mut *reader, &config, &mut keep).expect("the summary");

            let num_blocks = NUM_VARS_OF_MANY_VCF / NUM_VARS_PER_BLOCK;
            assert_eq!(calls.len(), 5, "a call after each of the 5 blocks");
            for (call, block) in calls.iter().zip(1..=num_blocks) {
                let given_the_filter = block.saturating_mul(NUM_VARS_PER_BLOCK);
                if threshold >= KEEPS_EVERY_VARIANT {
                    assert_eq!(call.num_vars, given_the_filter, "after block {block}");
                } else {
                    assert!(
                        call.num_vars < given_the_filter,
                        "the filter left variants out by block {block}"
                    );
                }
                assert_eq!(
                    call.filtering_stats,
                    vec![(
                        "maf",
                        FilteringStats {
                            vars_processed: given_the_filter,
                            vars_kept: call.num_vars,
                        }
                    )],
                    "the counts of the filter after block {block}"
                );
                let first_n = [maf_filter(threshold), PassStep::FirstN(call.num_vars)];
                assert_eq!(
                    call.result,
                    of_the_summary(Some(100), &first_n, all_three, chrom_lengths),
                    "the summary after block {block}, over the first {} variants",
                    call.num_vars
                );
            }
            assert_eq!(
                calls.last().map(|call| &call.result),
                Some(&written(&returned)),
                "the last call is what the pass returns"
            );
        }
    }
}

/// A pass whose filter kept no variant is refused with the error each of
/// the three gives for it, with the counts of the filter: the major allele
/// frequency of a variant is 0.5 at least, so a threshold of 0 keeps none.
#[test]
fn variants_summary_of_a_pass_with_no_variant_is_the_error_the_three_give() {
    let steps = [maf_filter(0.0)];
    let the_error_of = |pass: &dyn Fn(&mut dyn BlockReader) -> Error| {
        let mut reader = many_vcf_through(Some(100), &steps);
        format!("{:?}", pass(&mut *reader))
    };
    let of_the_summary = the_error_of(&|reader| {
        let all_three = Asked {
            per_var: true,
            per_individual: true,
            density: true,
        };
        let config = summary_config_of(reader, all_three, None);
        let error = calc_variants_summary(reader, &config, &mut nothing)
            .expect_err("a summary with no variant");
        assert!(
            matches!(&error, Error::PassGaveNoVariant { num_vars_of_the_source, filters }
                if *num_vars_of_the_source == NUM_VARS_OF_MANY_VCF && filters.len() == 1),
            "{error:?}"
        );
        error
    });
    let of_each = [
        the_error_of(&|reader| {
            let config = per_var_config_of(reader);
            calc_per_var_distribs(reader, &config).expect_err("the distributions")
        }),
        the_error_of(&|reader| calc_per_individual_stats(reader).expect_err("the rates")),
        the_error_of(&|reader| {
            calc_var_density(reader, WINDOW_SIZE, None).expect_err("the density")
        }),
    ];
    for of_one in of_each {
        assert_eq!(of_the_summary, of_one);
    }
}

/// The three statistics.
const ALL_THREE: Asked = Asked {
    per_var: true,
    per_individual: true,
    density: true,
};

/// A variant of chr1 past its length in the second block of 100 ends the
/// pass of the three with the error the density gives over the same
/// reader, after one call of the function, for the first block: with chr1
/// 5000 base pairs long, the 110th variant of `many.vcf`, at 5033, is past
/// it.
#[test]
fn variants_summary_of_a_variant_the_density_refuses_ends_the_pass_with_its_error() {
    let lengths = vec![("chr1".to_owned(), 5000), ("chr2".to_owned(), 25_000)];
    let mut reader = many_vcf_through(Some(100), &[]);
    let config = summary_config_of(&*reader, ALL_THREE, Some(&lengths));
    let mut num_calls = 0_u32;
    let mut count = |_: &dyn SoFar<VariantsSummary>| {
        num_calls = num_calls.saturating_add(1);
        Ok(())
    };
    let error = calc_variants_summary(&mut *reader, &config, &mut count)
        .expect_err("a variant past its length");

    let mut reader = many_vcf_through(Some(100), &[]);
    let of_the_density = calc_var_density(&mut *reader, WINDOW_SIZE, Some(&lengths))
        .expect_err("the density of a variant past its length");
    assert!(
        matches!(&error, Error::VarDensityVarPastTheLength { pos: 5033, .. }),
        "{error:?}"
    );
    assert_eq!(format!("{error:?}"), format!("{of_the_density:?}"));
    assert_eq!(num_calls, 1, "a call after the first block alone");
}

/// An error of the function at its second call ends the pass with that
/// error: no third call, and no block asked of the reader past the one the
/// pass reads ahead natively, the third.
#[test]
fn variants_summary_error_of_the_function_ends_the_pass() {
    let mut reader = Recording::of_many_vcf();
    let config = summary_config_of(&reader, ALL_THREE, None);
    let mut num_calls = 0_u32;
    let mut stop_at_the_second = |_: &dyn SoFar<VariantsSummary>| {
        num_calls = num_calls.saturating_add(1);
        if num_calls == 2 {
            return Err(Error::VarDensityChromLengthZero {
                chrom: "the error of the function".to_owned(),
                from: LengthsFrom::ChromLengths,
            });
        }
        Ok(())
    };
    let error = calc_variants_summary(&mut reader, &config, &mut stop_at_the_second)
        .expect_err("the pass is stopped");

    assert!(
        matches!(&error, Error::VarDensityChromLengthZero { chrom, .. }
            if chrom == "the error of the function"),
        "the error of the function, and not {error}"
    );
    assert_eq!(num_calls, 2, "no call after the one that failed");
    let num_blocks = reader.num_blocks_given();
    assert!(
        (2..=3).contains(&num_blocks),
        "the two blocks added and at most the one read ahead, not {num_blocks}"
    );
}

/// A `poly_threshold` of NaN and a window of 0 base pairs, both refused
/// before the pass, give the error of the threshold, the first of the
/// order the doc of `calc_variants_summary` gives, before the reader is
/// asked for anything.
#[test]
fn variants_summary_of_a_threshold_and_a_window_both_refused_gives_the_threshold_error() {
    let mut reader = Recording::of_many_vcf();
    let mut per_var = per_var_config_of(&reader);
    per_var.poly_threshold = f64::NAN;
    let config = VariantsSummaryConfig {
        per_var: Some(per_var),
        per_individual: true,
        density: Some((0, None)),
    };
    let error = calc_variants_summary(&mut reader, &config, &mut nothing)
        .expect_err("a threshold and a window refused");

    assert!(
        matches!(&error, Error::PolyThresholdOutOfRange { value } if value.is_nan()),
        "{error:?}"
    );
    assert_eq!(*reader.asked.lock().expect("what was asked"), Vec::new());
}

/// The first genotype of a block made an allele below the missing one,
/// which the rates refuse.
fn with_an_allele_below_the_missing_one(block: &mut Block) {
    if let Some(allele) = block.gts.first_mut() {
        *allele = -2;
    }
}

/// The distributions of `many.vcf` over one population of the individual
/// 999, which no row of `many.vcf`, of 50 individuals, holds and which the
/// distributions refuse in every block; `Pops::from_names` builds no such
/// population, so it is built here as it is held.
fn per_var_config_beyond_the_row(reader: &dyn BlockReader) -> PerVarDistribsConfig {
    let mut config = per_var_config_of(reader);
    config.pops = Pops {
        pops: vec![Pop {
            name: "beyond".to_owned(),
            individuals: vec![999],
            is_all: false,
        }],
    };
    config
}

/// The error of a summary of `asked` over `many.vcf` whose first block the
/// three refuse, each with an error of its own: the distributions an
/// individual beyond the row, the rates an allele below the missing one and
/// the density a variant past the length of chr1, 500 base pairs.
fn the_error_of_a_block_refused_by_the_three(asked: Asked) -> Error {
    let lengths = vec![("chr1".to_owned(), 500), ("chr2".to_owned(), 25_000)];
    let mut reader = Recording::of_many_vcf_altered(with_an_allele_below_the_missing_one);
    let config = VariantsSummaryConfig {
        per_var: asked
            .per_var
            .then(|| per_var_config_beyond_the_row(&reader)),
        per_individual: asked.per_individual,
        density: asked.density.then_some((WINDOW_SIZE, Some(lengths))),
    };
    calc_variants_summary(&mut reader, &config, &mut nothing).expect_err("a block refused")
}

/// A block that more than one statistic refuses ends the pass with the
/// error of the distributions, then of the rates, then of the density, each
/// the error its own function gives over the same reader.
#[test]
fn variants_summary_of_a_block_refused_by_several_gives_the_error_of_the_first_in_order() {
    let of = |per_var, per_individual, density| Asked {
        per_var,
        per_individual,
        density,
    };
    let lengths = vec![("chr1".to_owned(), 500), ("chr2".to_owned(), 25_000)];
    let mut reader = Recording::of_many_vcf_altered(with_an_allele_below_the_missing_one);
    let config = per_var_config_beyond_the_row(&reader);
    let of_the_distribs =
        calc_per_var_distribs(&mut reader, &config).expect_err("the distributions");
    let mut reader = Recording::of_many_vcf_altered(with_an_allele_below_the_missing_one);
    let of_the_rates = calc_per_individual_stats(&mut reader).expect_err("the rates");
    let mut reader = Recording::of_many_vcf_altered(with_an_allele_below_the_missing_one);
    let of_the_density =
        calc_var_density(&mut reader, WINDOW_SIZE, Some(&lengths)).expect_err("the density");
    assert!(
        matches!(
            &of_the_distribs,
            Error::IndividualBeyondTheVariant {
                individual: 999,
                ..
            }
        ),
        "{of_the_distribs:?}"
    );
    assert!(
        matches!(
            &of_the_rates,
            Error::AlleleBelowTheMissingOne { allele: -2 }
        ),
        "{of_the_rates:?}"
    );
    assert!(
        matches!(
            &of_the_density,
            Error::VarDensityVarPastTheLength { pos: 1000, .. }
        ),
        "{of_the_density:?}"
    );

    for (asked, expected) in [
        (of(true, true, true), &of_the_distribs),
        (of(true, true, false), &of_the_distribs),
        (of(true, false, true), &of_the_distribs),
        (of(true, false, false), &of_the_distribs),
        (of(false, true, true), &of_the_rates),
        (of(false, true, false), &of_the_rates),
        (of(false, false, true), &of_the_density),
    ] {
        assert_eq!(
            format!("{:?}", the_error_of_a_block_refused_by_the_three(asked)),
            format!("{expected:?}"),
            "{asked:?}"
        );
    }
}

/// The block made a block of no variants, with the chromosome and the
/// position it was given, both empty.
fn of_no_variants(block: &mut Block) {
    *block = Block {
        num_vars: 0,
        num_individuals: block.num_individuals,
        ploidy: block.ploidy,
        gts: Vec::new(),
        chrom: block.chrom.as_ref().map(|_| Vec::new()),
        pos: block.pos.as_ref().map(|_| Vec::new()),
        id: None,
        alleles: None,
        qual: None,
        passed: None,
        vcf_text: None,
    };
}

/// A block of no variants, a defect of the reader that gives it, is refused
/// by the summary of each of the three alone, the density among them, with
/// the error each of their functions gives for it.
#[test]
fn variants_summary_refuses_a_block_of_no_variants() {
    let of = |per_var, per_individual, density| Asked {
        per_var,
        per_individual,
        density,
    };
    for asked in [
        of(true, false, false),
        of(false, true, false),
        of(false, false, true),
    ] {
        let mut reader = Recording::of_many_vcf_altered(of_no_variants);
        let config = summary_config_of(&reader, asked, None);
        let error = calc_variants_summary(&mut reader, &config, &mut nothing)
            .expect_err("a block of no variants");
        assert!(
            matches!(error, Error::ReaderGaveABlockOfNoVariants),
            "{asked:?}: {error:?}"
        );
    }
}
