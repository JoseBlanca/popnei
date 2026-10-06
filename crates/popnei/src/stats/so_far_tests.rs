//! The tests of the result so far of the three calculations that add up
//! over the blocks of a pass, each with `so_far` in its name, from "The
//! result so far and the three in one pass" of "The Rust interface" of
//! `docs/specs/stats.md` and "The result so far" of
//! `docs/specs/js_sources.md`.
//!
//! Each pass reads `many.vcf`, 500 variants of 50 diploid individuals, in
//! blocks of 100 through the filter of the major allele frequency, and the
//! result of each call is compared with the result of the same calculation
//! over a second reader of the same file that ends where the call was, the
//! filter of the first n put after the same filter.

use super::fixtures::vcf_reader_of;
use super::{
    AfterABlock, ExpHet, HistBins, LengthsFrom, Maf, ObsHet, PerIndividualStats, PerVarDistribs,
    PerVarDistribsConfig, PerVarStat, Pops, SoFar, VarDensity, calc_per_individual_stats,
    calc_per_individual_stats_with, calc_per_var_distribs, calc_per_var_distribs_with,
    calc_var_density, calc_var_density_with,
};
use crate::block::BlockReader;
use crate::error::{Error, Result};
use crate::filters::{FilteringStats, PassStep, VarFilteringCriterion, chain_of};

/// The variants of `many.vcf`, read in blocks of [`NUM_VARS_PER_BLOCK`].
const NUM_VARS_OF_MANY_VCF: u64 = 500;

/// The size of the blocks the passes of these tests read, which makes five
/// blocks of `many.vcf`.
const NUM_VARS_PER_BLOCK: u64 = 100;

/// A threshold of the major allele frequency that keeps every variant of
/// `many.vcf`, so that the five calls come after 100, 200, 300, 400 and 500
/// variants.
const KEEPS_EVERY_VARIANT: f64 = 1.0;

/// A threshold that leaves some variants of every block out, so that what
/// the pass has read and what its filter was given differ.
const KEEPS_SOME_VARIANTS: f64 = 0.9;

/// What a call of the function after a block was given, with the result
/// written with `{:?}`: Rust writes a `f64` with the fewest digits that
/// read back as the same number, so two results written the same hold the
/// same numbers to the bit.
#[derive(Debug, PartialEq)]
struct ACall {
    num_vars: u64,
    filtering_stats: Vec<(&'static str, FilteringStats)>,
    result: String,
}

/// The chain of a pass over `many.vcf` in blocks of 100, with every
/// variant given, through `steps`.
fn many_vcf_through(steps: &[PassStep]) -> Box<dyn BlockReader> {
    let source = vcf_reader_of("vcf/many.vcf", Some(100));
    chain_of(Box::new(source), steps).expect("the chain over many.vcf")
}

/// The filter of the major allele frequency with `threshold`.
fn maf_filter(threshold: f64) -> PassStep {
    PassStep::VarFilter(VarFilteringCriterion::MaxMaf(threshold))
}

/// The function after a block that keeps what each call was given.
fn kept_in(calls: &mut Vec<ACall>) -> impl FnMut(&dyn SoFar<String>) -> Result<()> + '_ {
    |so_far| {
        calls.push(ACall {
            num_vars: so_far.num_vars(),
            filtering_stats: so_far.filtering_stats(),
            result: so_far.result()?,
        });
        Ok(())
    }
}

/// The calls of a pass with the filter of the major allele frequency at
/// `threshold`, and what the pass returned, written with `{:?}`.
///
/// `pass` runs the calculation over a chain with a function after each
/// block, which it is given as a function of the result written with
/// `{:?}`, so that one function serves the three calculations.
fn the_calls(
    threshold: f64,
    pass: impl FnOnce(&mut dyn BlockReader, AfterABlock<'_, String>) -> Result<String>,
) -> (Vec<ACall>, String) {
    let mut chain = many_vcf_through(&[maf_filter(threshold)]);
    let mut calls = Vec::new();
    let returned = pass(&mut *chain, &mut kept_in(&mut calls)).expect("the pass");
    (calls, returned)
}

/// That the calls of a pass with the filter at `threshold` are one for each
/// block, the counts of the filter those after that block, and the result
/// that of `over_the_first` of the variants the call has read, the last one
/// what the pass returned.
fn assert_the_calls(
    threshold: f64,
    (calls, returned): (Vec<ACall>, String),
    over_the_first: impl Fn(&mut dyn BlockReader) -> String,
) {
    let num_blocks = NUM_VARS_OF_MANY_VCF / NUM_VARS_PER_BLOCK;
    assert_eq!(calls.len(), 5, "a call after each of the 5 blocks");
    for (call, block) in calls.iter().zip(1..=num_blocks) {
        let given_the_filter = block.saturating_mul(NUM_VARS_PER_BLOCK);
        if threshold >= KEEPS_EVERY_VARIANT {
            assert_eq!(
                call.num_vars, given_the_filter,
                "the variants after block {block}"
            );
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
        let mut first_n =
            many_vcf_through(&[maf_filter(threshold), PassStep::FirstN(call.num_vars)]);
        assert_eq!(
            call.result,
            over_the_first(&mut *first_n),
            "the result after block {block}, over the first {} variants",
            call.num_vars
        );
    }
    assert_eq!(
        calls.last().map(|call| &call.result),
        Some(&returned),
        "the last call is what the pass returns"
    );
}

/// That the pass ends with the error of the function at its second call,
/// and calls it no third time.
fn assert_an_error_ends_the_pass(
    pass: impl FnOnce(&mut dyn BlockReader, AfterABlock<'_, String>) -> Result<String>,
) {
    let mut chain = many_vcf_through(&[maf_filter(KEEPS_SOME_VARIANTS)]);
    let mut num_calls = 0_u32;
    let mut stop_at_the_second = |_: &dyn SoFar<String>| {
        num_calls = num_calls.saturating_add(1);
        if num_calls == 2 {
            return Err(Error::VarDensityChromLengthZero {
                chrom: "the error of the function".to_owned(),
                from: LengthsFrom::ChromLengths,
            });
        }
        Ok(())
    };
    let error = pass(&mut *chain, &mut stop_at_the_second).expect_err("the pass is stopped");
    assert!(
        matches!(&error, Error::VarDensityChromLengthZero { chrom, .. }
            if chrom == "the error of the function"),
        "the error of the function, and not {error}"
    );
    assert_eq!(num_calls, 2, "no call after the one that failed");
}

/// The calls of a pass of a calculation whose result is `T`, given to a
/// function of the result written with `{:?}`.
fn written<'a, T: std::fmt::Debug>(
    after_a_block: &'a mut dyn FnMut(&dyn SoFar<String>) -> Result<()>,
) -> impl FnMut(&dyn SoFar<T>) -> Result<()> + 'a {
    move |so_far| after_a_block(&Written(so_far))
}

/// A result so far seen as its result written with `{:?}`.
struct Written<'a, T>(&'a dyn SoFar<T>);

impl<T: std::fmt::Debug> SoFar<String> for Written<'_, T> {
    fn num_vars(&self) -> u64 {
        self.0.num_vars()
    }

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        self.0.filtering_stats()
    }

    fn result(&self) -> Result<String> {
        Ok(format!("{:?}", self.0.result()?))
    }
}

/// The six statistics over `popA` and `popB` of `many.vcf`, with numbers
/// that are none of their defaults: 10 called genotypes, a threshold of
/// 0.9 and ten bins.
fn the_config_of(reader: &dyn BlockReader) -> PerVarDistribsConfig {
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
        bins: HistBins::linear(0.0, 1.0, 10).expect("ten bins"),
        obs_het: ObsHet::new(10),
        maf: Maf::new(2, 10).expect("the maf of diploid variants"),
        exp_het: ExpHet::new(2, 2, 10).expect("the expected heterozygosity of diploid variants"),
        poly_threshold: 0.9,
    }
}

fn per_var_distribs_with(
    chain: &mut dyn BlockReader,
    after_a_block: AfterABlock<'_, String>,
) -> Result<String> {
    let config = the_config_of(chain);
    let after: AfterABlock<'_, PerVarDistribs> = &mut written(after_a_block);
    calc_per_var_distribs_with(chain, &config, after).map(|result| format!("{result:?}"))
}

fn per_var_distribs(chain: &mut dyn BlockReader) -> String {
    let config = the_config_of(chain);
    let result = calc_per_var_distribs(chain, &config).expect("the distributions");
    format!("{result:?}")
}

#[test]
fn per_var_distribs_so_far_after_each_block_is_the_result_over_the_variants_read() {
    for threshold in [KEEPS_EVERY_VARIANT, KEEPS_SOME_VARIANTS] {
        let calls = the_calls(threshold, per_var_distribs_with);
        assert_the_calls(threshold, calls, per_var_distribs);
    }
}

#[test]
fn per_var_distribs_so_far_error_of_the_function_ends_the_pass() {
    assert_an_error_ends_the_pass(per_var_distribs_with);
}

fn per_individual_stats_with(
    chain: &mut dyn BlockReader,
    after_a_block: AfterABlock<'_, String>,
) -> Result<String> {
    let after: AfterABlock<'_, PerIndividualStats> = &mut written(after_a_block);
    calc_per_individual_stats_with(chain, after).map(|result| format!("{result:?}"))
}

fn per_individual_stats(chain: &mut dyn BlockReader) -> String {
    let result = calc_per_individual_stats(chain).expect("the per individual statistics");
    format!("{result:?}")
}

#[test]
fn per_individual_stats_so_far_after_each_block_is_the_result_over_the_variants_read() {
    for threshold in [KEEPS_EVERY_VARIANT, KEEPS_SOME_VARIANTS] {
        let calls = the_calls(threshold, per_individual_stats_with);
        assert_the_calls(threshold, calls, per_individual_stats);
    }
}

#[test]
fn per_individual_stats_so_far_error_of_the_function_ends_the_pass() {
    assert_an_error_ends_the_pass(per_individual_stats_with);
}

/// The width of the windows of the density of these tests, 1000 base
/// pairs, which gives `many.vcf`, whose last position is 19463, twenty
/// windows on each chromosome.
const WINDOW_SIZE: u64 = 1000;

/// Lengths that make the windows of both chromosomes before the pass, all
/// of them past the last variant of `many.vcf`.
fn the_lengths_given() -> Vec<(String, u64)> {
    vec![("chr2".to_owned(), 25_000), ("chr1".to_owned(), 30_000)]
}

fn var_density_with(
    chrom_lengths: Option<&[(String, u64)]>,
) -> impl Fn(&mut dyn BlockReader, AfterABlock<'_, String>) -> Result<String> {
    move |chain, after_a_block| {
        let after: AfterABlock<'_, VarDensity> = &mut written(after_a_block);
        calc_var_density_with(chain, WINDOW_SIZE, chrom_lengths, after)
            .map(|result| format!("{result:?}"))
    }
}

fn var_density(chrom_lengths: Option<&[(String, u64)]>) -> impl Fn(&mut dyn BlockReader) -> String {
    move |chain| {
        let result = calc_var_density(chain, WINDOW_SIZE, chrom_lengths)
            .expect("the density of the variants");
        format!("{result:?}")
    }
}

/// The lengths of the header of `many.vcf`, which has none, so the windows
/// grow with the variants; no lengths given, which is the same; and
/// lengths given, whose windows are all there from the first call.
#[test]
fn var_density_so_far_after_each_block_is_the_result_over_the_variants_read() {
    let given = the_lengths_given();
    for chrom_lengths in [None, Some(&[][..]), Some(&given[..])] {
        for threshold in [KEEPS_EVERY_VARIANT, KEEPS_SOME_VARIANTS] {
            let calls = the_calls(threshold, var_density_with(chrom_lengths));
            assert_the_calls(threshold, calls, var_density(chrom_lengths));
        }
    }
}

#[test]
fn var_density_so_far_error_of_the_function_ends_the_pass() {
    assert_an_error_ends_the_pass(var_density_with(None));
}
