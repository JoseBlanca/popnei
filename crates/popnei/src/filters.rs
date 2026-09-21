//! The filters that keep the variants of a dataset that pass a threshold,
//! and the counts of what each one was given and kept.
//!
//! A filter works out one number for every variant, over all the
//! individuals of the dataset, and keeps the variant when that number is at
//! most the threshold the user gave: [`VarFilteringCriterion`] is which
//! number it is, with the threshold, and [`VarFilter`] is the filter of one
//! pass over the source, which keeps the variants of a block that pass and
//! counts what it was given and kept. The counts are the [`FilteringStats`]
//! that every reader gives for the filters between it and its source,
//! through
//! [`BlockReader::filtering_stats`](crate::block::BlockReader::filtering_stats).
//!
//! The reader that puts a filter over another reader is being written.
//!
//! `docs/specs/filters.md` has the design, and the row `filters` of section
//! 9 of `docs/architecture.md` where the module sits.

use crate::block::{Block, BlockSize};
use crate::error::{Error, Result};
use crate::variant::{AlleleCounts, Needs, count_alleles, count_gts};

/// How many variants a filter was given and how many of them it kept, over
/// every block it has taken since it was built.
///
/// A filter belongs to one pass over the source, one reading of it from its
/// start, so these are the counts of that pass alone. `vars_kept` is at
/// most `vars_processed`. Read while the pass runs, they are of the
/// variants the filter has seen, which can be more than the ones the
/// consumer of the pass has got.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FilteringStats {
    /// The variants the filter was given.
    pub vars_processed: u64,
    /// Those of them that passed its threshold.
    pub vars_kept: u64,
}

/// Which number of a variant a filter compares with a threshold, with the
/// largest value of that number that keeps the variant.
///
/// A variant stays when its number is at most the threshold, so one whose
/// number is exactly the threshold stays. Each number is one count of the
/// variant divided by another, so a threshold is a number from 0 to 1. A
/// variant that has no number, one with no called allele for the major
/// allele frequency and one with no called genotype for the observed
/// heterozygosity, is not kept, whatever the threshold.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VarFilteringCriterion {
    /// Missing genotypes divided by all the individuals of the dataset, and
    /// not the ones that were called at the variant. A genotype is missing
    /// when one of its alleles at least was not called, so a half called
    /// genotype, `0/.` in a VCF, is missing.
    MaxMissingRate(f64),
    /// The count of the commonest allele divided by the called alleles.
    /// Every allele of a multiallelic variant has its own count, and an
    /// allele is counted wherever it was called, in a half called genotype
    /// too. When two alleles tie for the largest count the frequency is the
    /// same whichever of them is called the major one.
    MaxMaf(f64),
    /// Heterozygous genotypes divided by the called genotypes. A genotype
    /// is heterozygous when it is called and its alleles are not all the
    /// same, at any ploidy.
    MaxObsHet(f64),
}

impl VarFilteringCriterion {
    /// `"missing_data"`, `"maf"` or `"obs_het"`: the name under which the
    /// counts of the filter reach a Python or a TypeScript user, and the
    /// name by which a chain of readers is asked whether it holds a filter
    /// of this kind already.
    #[must_use]
    pub fn kind(&self) -> &'static str {
        match self {
            VarFilteringCriterion::MaxMissingRate(_) => "missing_data",
            VarFilteringCriterion::MaxMaf(_) => "maf",
            VarFilteringCriterion::MaxObsHet(_) => "obs_het",
        }
    }

    /// The largest value of the number of a variant that keeps it.
    pub(crate) fn threshold(&self) -> f64 {
        match self {
            VarFilteringCriterion::MaxMissingRate(threshold)
            | VarFilteringCriterion::MaxMaf(threshold)
            | VarFilteringCriterion::MaxObsHet(threshold) => *threshold,
        }
    }
}

/// The filter of one pass over the source: it keeps the variants of a block
/// that pass its threshold and counts how many it was given and how many it
/// kept.
///
/// It is an object of its own, apart from the reader that puts it over a
/// source, so that the rule and its counts are worked out on a block that
/// no reader gave. Every pass builds its own, so no count is shared between
/// two passes.
#[derive(Debug)]
pub struct VarFilter {
    criterion: VarFilteringCriterion,
    stats: FilteringStats,
}

impl VarFilter {
    /// The filter that keeps the variants whose number is at most the
    /// threshold of `criterion`, with both its counts at 0.
    ///
    /// # Errors
    ///
    /// When the threshold is NaN, below 0 or above 1: the error names the
    /// criterion and the value.
    pub fn new(criterion: VarFilteringCriterion) -> Result<VarFilter> {
        let threshold = criterion.threshold();
        // A NaN is in no range, so this one comparison refuses the three
        // thresholds that are not a number from 0 to 1.
        if !(0.0..=1.0).contains(&threshold) {
            return Err(Error::VarFilterThresholdOutOfRange {
                kind: criterion.kind(),
                threshold,
            });
        }
        Ok(VarFilter {
            criterion,
            stats: FilteringStats::default(),
        })
    }

    /// Which number of a variant it compares, with its threshold.
    #[must_use]
    pub fn criterion(&self) -> VarFilteringCriterion {
        self.criterion
    }

    /// The variants of the block that pass, kept in it in their order, and
    /// the others dropped: the genotypes and every column of the block are
    /// compacted in place, and nothing is allocated for a variant.
    ///
    /// The variants of the block are added to the counts, and the ones that
    /// stayed to the ones kept. A block of no variants is left as it is.
    ///
    /// # Errors
    ///
    /// When the arrays of the block are not of the size the block states,
    /// which [`Block::check`] finds, and when the block has variants and no
    /// genotypes, which is the error of a field that is not in the block.
    /// After either the block is as it was and nothing was added to the
    /// counts.
    pub fn filter_block(&mut self, block: &mut Block) -> Result<()> {
        // The rows are cut out of the genotypes by the sizes the block
        // states, so those sizes are checked before anything is read.
        block.check()?;
        let processed = block.num_vars;
        if processed == 0 {
            return Ok(());
        }
        if block.gts.is_empty() {
            return Err(Error::FieldsNotInTheBlock { fields: Needs::GTS });
        }
        let alleles_per_var = block
            .num_individuals
            .checked_mul(block.ploidy)
            .ok_or(Error::BlockTooLarge {
                num_vars_per_block: block.num_vars,
                num_individuals: block.num_individuals,
                ploidy: block.ploidy,
                // The block is here, so its size is one that a reader was
                // given and took: what a caller does about it is ask for
                // fewer variants.
                size: BlockSize::AskedFor,
            })?
            // `check` passed and the genotypes are not empty, so they are
            // the variants of the block times this number and it is one
            // allele at least; the rows are cut by it, and a cut of 0 is
            // what the standard library refuses with a panic.
            .max(1);
        let keep = keep_of_the_rows(
            self.criterion,
            &block.gts,
            alleles_per_var,
            block.num_individuals,
            block.ploidy,
        )?;
        let kept = keep.iter().filter(|keep_it| **keep_it).count();
        block.retain_vars(&keep)?;
        self.stats.vars_processed = self
            .stats
            .vars_processed
            .saturating_add(u64::try_from(processed).unwrap_or(u64::MAX));
        self.stats.vars_kept = self
            .stats
            .vars_kept
            .saturating_add(u64::try_from(kept).unwrap_or(u64::MAX));
        Ok(())
    }

    /// How many variants it was given and how many it kept, over every
    /// block it has taken since it was built.
    #[must_use]
    pub fn stats(&self) -> FilteringStats {
        self.stats
    }
}

/// Which rows of the genotypes of a block pass the criterion, one value for
/// each variant, in the order of the block.
///
/// Natively the rows are read on the threads of rayon, as section 3 of
/// `docs/architecture.md` asks: no row reads another and each one gives one
/// value of the answer, which is collected in the order of the block, so
/// neither the values nor the counts depend on how many threads there are.
/// The threads are those of the pool the caller is running in, and rayon's
/// global pool only when the caller is in none.
///
/// `gts` holds the rows of the block, `alleles_per_var` alleles each, and
/// `alleles_per_var` is 1 or more.
///
/// # Errors
///
/// What the counts of one variant refuse: a ploidy of 0, genotypes that are
/// not a whole number of genotypes of the ploidy, and an allele below the
/// missing one.
#[cfg(not(target_family = "wasm"))]
fn keep_of_the_rows(
    criterion: VarFilteringCriterion,
    gts: &[i8],
    alleles_per_var: usize,
    num_individuals: usize,
    ploidy: usize,
) -> Result<Vec<bool>> {
    use rayon::iter::ParallelIterator;
    use rayon::slice::ParallelSlice;

    gts.par_chunks_exact(alleles_per_var)
        .map(|row| keeps(criterion, row, num_individuals, ploidy))
        .collect()
}

/// The same values, with the rows read one after another, which is what
/// wasm does: it has no threads.
#[cfg(target_family = "wasm")]
fn keep_of_the_rows(
    criterion: VarFilteringCriterion,
    gts: &[i8],
    alleles_per_var: usize,
    num_individuals: usize,
    ploidy: usize,
) -> Result<Vec<bool>> {
    keep_of_the_rows_one_by_one(criterion, gts, alleles_per_var, num_individuals, ploidy)
}

/// The rows read one after another, which is what wasm does and what the
/// test that compares the two ways of reading them calls.
#[cfg_attr(
    all(not(target_family = "wasm"), not(test)),
    expect(
        dead_code,
        reason = "in wasm it is how the rows of a block are read, and natively it is what \
                  the test that compares the two ways of reading them calls; outside the \
                  tests and outside wasm nothing calls it"
    )
)]
fn keep_of_the_rows_one_by_one(
    criterion: VarFilteringCriterion,
    gts: &[i8],
    alleles_per_var: usize,
    num_individuals: usize,
    ploidy: usize,
) -> Result<Vec<bool>> {
    gts.chunks_exact(alleles_per_var)
        .map(|row| keeps(criterion, row, num_individuals, ploidy))
        .collect()
}

/// Whether the variant whose genotypes are `gts` passes the criterion.
///
/// `gts` is one row of the genotypes of a block, the alleles of one
/// individual after those of the individual before it, `ploidy` alleles
/// each, and `num_individuals` is the individuals of the dataset. A variant
/// that has no number, one with no called allele or no called genotype, is
/// not kept, which the comparison of a NaN gives too.
///
/// The number is one division of the two counts as `f64` and is compared
/// with `<=`, as in pyNei, and not against a product of the threshold and
/// the denominator: 29 missing genotypes of 100 individuals pass a
/// threshold of 0.29 as 29/100, and would not as 29 <= 0.29 * 100, which is
/// 28.999999999999996.
fn keeps(
    criterion: VarFilteringCriterion,
    gts: &[i8],
    num_individuals: usize,
    ploidy: usize,
) -> Result<bool> {
    let number = match criterion {
        VarFilteringCriterion::MaxMissingRate(_) => {
            let counts = count_gts(gts, ploidy)?;
            // The individuals of the dataset, and not the ones called at
            // this variant. A variant of no individual gives 0/0, a NaN,
            // and is not kept.
            f64::from(counts.missing) / num_individuals as f64
        }
        VarFilteringCriterion::MaxMaf(_) => {
            let mut counts: AlleleCounts = [0; 128];
            let called_alleles = count_alleles(gts, &mut counts)?;
            if called_alleles == 0 {
                return Ok(false);
            }
            let largest = counts.iter().copied().max().unwrap_or(0);
            f64::from(largest) / f64::from(called_alleles)
        }
        VarFilteringCriterion::MaxObsHet(_) => {
            let counts = count_gts(gts, ploidy)?;
            if counts.called == 0 {
                return Ok(false);
            }
            f64::from(counts.het) / f64::from(counts.called)
        }
    };
    Ok(number <= criterion.threshold())
}

#[cfg(test)]
mod tests {
    use super::{
        FilteringStats, VarFilter, VarFilteringCriterion, keep_of_the_rows,
        keep_of_the_rows_one_by_one,
    };
    use crate::block::Block;
    use crate::error::Error;
    use crate::variant::{MISSING_ALLELE, Needs};

    use VarFilteringCriterion::{MaxMaf, MaxMissingRate, MaxObsHet};

    /// The six variants of five diploid individuals of the worked example
    /// of "How it is verified" of `docs/specs/filters.md`, each at the
    /// position of its number in that table, so that a test names the
    /// variants that stayed by the numbers of the spec. `-1` is an allele
    /// that was not called, the `.` of a VCF.
    ///
    /// Their missing rates are 0.2, 0.4, 0.2, 1, 0 and 1; their major
    /// allele frequencies 8/9, 6/7, 2/8, none, 8/10 and 1/1; and their
    /// observed heterozygosities 1/4, 1/3, 4/4, none, 0 and none.
    const THE_WORKED_EXAMPLE: [(u64, [i8; 10]); 6] = [
        // 0/0 0/1 0/0 0/0 0/.
        (1, [0, 0, 0, 1, 0, 0, 0, 0, 0, -1]),
        // 0/0 0/1 0/0 ./. 0/.
        (2, [0, 0, 0, 1, 0, 0, -1, -1, 0, -1]),
        // 0/1 2/3 0/1 2/3 ./.
        (3, [0, 1, 2, 3, 0, 1, 2, 3, -1, -1]),
        // ./. ./. ./. ./. ./.
        (4, [-1; 10]),
        // 0/0 0/0 0/0 0/0 1/1
        (5, [0, 0, 0, 0, 0, 0, 0, 0, 1, 1]),
        // 0/. ./. ./. ./. ./.
        (6, [0, -1, -1, -1, -1, -1, -1, -1, -1, -1]),
    ];

    /// A block of the variants given, each with its position, of
    /// `num_individuals` individuals of the ploidy `ploidy`. It holds the
    /// genotypes, the chromosome and the position, which is what a filter
    /// over a reader asked for the positions gets.
    fn block_of(variants: &[(u64, &[i8])], num_individuals: usize, ploidy: usize) -> Block {
        let mut gts = Vec::new();
        let mut chrom = Vec::new();
        let mut pos = Vec::new();
        for (position, row) in variants {
            gts.extend_from_slice(row);
            chrom.push(0);
            pos.push(*position);
        }
        Block {
            num_vars: variants.len(),
            num_individuals,
            ploidy,
            gts,
            chrom: Some(chrom),
            pos: Some(pos),
            id: None,
            alleles: None,
            qual: None,
        }
    }

    /// The block of the six variants of the worked example.
    fn the_worked_example() -> Block {
        let variants: Vec<(u64, &[i8])> = THE_WORKED_EXAMPLE
            .iter()
            .map(|(pos, gts)| (*pos, gts.as_slice()))
            .collect();
        block_of(&variants, 5, 2)
    }

    /// The positions of the variants of the block, which in the worked
    /// example are the numbers the spec's table gives them.
    fn positions_of(block: &Block) -> Vec<u64> {
        block.pos.clone().unwrap_or_default()
    }

    /// The variants of the worked example that each filter keeps, and the
    /// counts of the filter, which the paragraph under the table of "How it
    /// is verified" of `docs/specs/filters.md` gives: the numbers are those
    /// of pyNei's `_calc_gt_is_missing`, `_calc_maf_per_var` and
    /// `_calc_obs_het_per_var` on these genotypes, and bcftools 1.24 keeps
    /// the same variants at 0.4, 0.88 and 0.25.
    #[test]
    fn each_filter_keeps_the_variants_of_the_worked_example_at_every_threshold() {
        let kept = |criterion| {
            let mut block = the_worked_example();
            let mut filter = VarFilter::new(criterion).unwrap();
            filter.filter_block(&mut block).unwrap();
            (positions_of(&block), filter.stats())
        };

        assert_eq!(kept(MaxMissingRate(0.0)).0, [5]);
        assert_eq!(kept(MaxMissingRate(0.2)).0, [1, 3, 5]);
        assert_eq!(kept(MaxMissingRate(0.4)).0, [1, 2, 3, 5]);
        assert_eq!(kept(MaxMissingRate(1.0)).0, [1, 2, 3, 4, 5, 6]);

        assert_eq!(kept(MaxMaf(0.25)).0, [3]);
        assert_eq!(kept(MaxMaf(0.8)).0, [3, 5]);
        assert_eq!(kept(MaxMaf(0.88)).0, [2, 3, 5]);
        assert_eq!(kept(MaxMaf(1.0)).0, [1, 2, 3, 5, 6]);

        assert_eq!(kept(MaxObsHet(0.0)).0, [5]);
        assert_eq!(kept(MaxObsHet(0.25)).0, [1, 5]);
        assert_eq!(kept(MaxObsHet(1.0)).0, [1, 2, 3, 5]);

        // The counts of one block are its variants and the ones that
        // stayed.
        assert_eq!(
            kept(MaxMaf(0.8)).1,
            FilteringStats {
                vars_processed: 6,
                vars_kept: 2,
            }
        );
    }

    /// A variant whose number is exactly the threshold stays: the major
    /// allele frequency of the variant 5 is 8/10 and the observed
    /// heterozygosity of the variant 1 is 1/4, and each of them stays at
    /// its own number and goes at a threshold below it.
    #[test]
    fn a_variant_whose_number_is_exactly_the_threshold_stays() {
        let kept = |criterion| {
            let mut block = the_worked_example();
            VarFilter::new(criterion)
                .unwrap()
                .filter_block(&mut block)
                .unwrap();
            positions_of(&block)
        };
        assert_eq!(kept(MaxMaf(0.8)), [3, 5]);
        assert_eq!(kept(MaxMaf(0.79)), [3]);
        assert_eq!(kept(MaxObsHet(0.25)), [1, 5]);
        assert_eq!(kept(MaxObsHet(0.24)), [5]);
        assert_eq!(kept(MaxMissingRate(0.2)), [1, 3, 5]);
        assert_eq!(kept(MaxMissingRate(0.19)), [5]);
    }

    /// The block is compacted in place: the genotypes that are left are
    /// those of the variants that stayed, in their order, and every column
    /// holds one entry for each of them.
    #[test]
    fn the_genotypes_and_the_columns_that_are_left_are_those_of_the_variants_that_stayed() {
        let mut block = the_worked_example();
        VarFilter::new(MaxMaf(0.8))
            .unwrap()
            .filter_block(&mut block)
            .unwrap();
        let mut expected = Vec::new();
        expected.extend_from_slice(&THE_WORKED_EXAMPLE[2].1);
        expected.extend_from_slice(&THE_WORKED_EXAMPLE[4].1);
        assert_eq!(block.num_vars, 2);
        assert_eq!(block.gts, expected);
        assert!(block.check().is_ok());
        assert_eq!(block.chrom.unwrap(), [0, 0]);
        assert_eq!(block.pos.unwrap(), [3, 5]);
    }

    /// The three filters of "How it is verified" of the counts on the
    /// worked example, at 0.4, 0.88 and 0.25 in that order: the pairs are 6
    /// and 4, 4 and 3, and 3 and 1, and the variant 5 is the one kept.
    #[test]
    fn the_three_filters_chained_on_the_worked_example_give_their_counts_and_keep_the_variant_5() {
        let mut block = the_worked_example();
        let mut missing_data = VarFilter::new(MaxMissingRate(0.4)).unwrap();
        let mut maf = VarFilter::new(MaxMaf(0.88)).unwrap();
        let mut obs_het = VarFilter::new(MaxObsHet(0.25)).unwrap();
        missing_data.filter_block(&mut block).unwrap();
        maf.filter_block(&mut block).unwrap();
        obs_het.filter_block(&mut block).unwrap();

        assert_eq!(positions_of(&block), [5]);
        let pair = |vars_processed, vars_kept| FilteringStats {
            vars_processed,
            vars_kept,
        };
        assert_eq!(missing_data.stats(), pair(6, 4));
        assert_eq!(maf.stats(), pair(4, 3));
        assert_eq!(obs_het.stats(), pair(3, 1));
    }

    /// The counts are of every block the filter was given, and a filter
    /// just built has counted nothing.
    #[test]
    fn the_counts_add_up_over_the_blocks_the_filter_was_given() {
        let mut filter = VarFilter::new(MaxMissingRate(0.4)).unwrap();
        assert_eq!(filter.stats(), FilteringStats::default());
        assert_eq!(filter.criterion(), MaxMissingRate(0.4));
        for _ in 0..3 {
            let mut block = the_worked_example();
            filter.filter_block(&mut block).unwrap();
            assert_eq!(block.num_vars, 4);
        }
        assert_eq!(
            filter.stats(),
            FilteringStats {
                vars_processed: 18,
                vars_kept: 12,
            }
        );
    }

    /// A threshold is a number from 0 to 1, both included, and the error
    /// names the criterion and the value, which is what a user needs in
    /// order to see that they wrote 95 for 0.95.
    #[test]
    fn a_threshold_that_is_not_a_number_from_0_to_1_is_refused_with_its_criterion_and_its_value() {
        let refused = |criterion, kind: &str, value: &str| {
            let error = VarFilter::new(criterion).unwrap_err();
            let message = error.to_string();
            assert!(
                matches!(error, Error::VarFilterThresholdOutOfRange { kind: of_the_error, .. } if of_the_error == kind),
                "{message}"
            );
            assert!(message.contains(kind), "{message}");
            assert!(message.contains(value), "{message}");
        };
        refused(MaxMissingRate(-0.1), "missing_data", "-0.1");
        refused(MaxMaf(1.5), "maf", "1.5");
        refused(MaxObsHet(f64::NAN), "obs_het", "NaN");
        refused(MaxMaf(f64::INFINITY), "maf", "inf");

        // The two ends of the range are thresholds: 0 keeps the variants
        // whose number is 0, and 1 the ones that have a number at all.
        assert!(VarFilter::new(MaxMissingRate(0.0)).is_ok());
        assert!(VarFilter::new(MaxMaf(1.0)).is_ok());
        assert!(VarFilter::new(MaxObsHet(0.5)).is_ok());
    }

    /// The rows are cut out of the genotypes by the sizes the block states,
    /// so a block whose arrays are not of its size is refused before
    /// anything is read, and the counts say nothing about it.
    #[test]
    fn a_block_whose_arrays_are_not_of_its_size_is_an_error_and_is_left_as_it_was() {
        let mut block = the_worked_example();
        // The arrays hold six variants and the block says five.
        block.num_vars = 5;
        let mut filter = VarFilter::new(MaxMissingRate(1.0)).unwrap();
        let error = filter.filter_block(&mut block).unwrap_err();
        assert!(
            matches!(
                error,
                Error::BlockArrayOfAnotherSize {
                    array: "gts",
                    found: 60,
                    expected: 50,
                }
            ),
            "{error}"
        );
        assert_eq!(filter.stats(), FilteringStats::default());
        assert_eq!(block.gts.len(), 60);
        assert_eq!(positions_of(&block), [1, 2, 3, 4, 5, 6]);
    }

    /// A filter always needs the genotypes, so a block that has variants
    /// and no genotypes is the error of a field that is not in the block,
    /// which names `gts`.
    #[test]
    fn a_block_with_variants_and_no_genotypes_is_the_error_of_a_field_that_is_not_there() {
        let mut block = the_worked_example();
        block.gts = Vec::new();
        let mut filter = VarFilter::new(MaxObsHet(1.0)).unwrap();
        let error = filter.filter_block(&mut block).unwrap_err();
        assert!(
            matches!(error, Error::FieldsNotInTheBlock { fields } if fields == Needs::GTS),
            "{error}"
        );
        assert!(error.to_string().contains("gts"), "{error}");
        assert_eq!(filter.stats(), FilteringStats::default());
        assert_eq!(block.num_vars, 6);
        assert_eq!(positions_of(&block), [1, 2, 3, 4, 5, 6]);
    }

    /// A block of no variants is left as it is and adds nothing to the
    /// counts. No reader of popnei gives one, and a block that a filter
    /// before this one emptied is not given on.
    #[test]
    fn a_block_of_no_variants_is_left_as_it_is() {
        let mut block = block_of(&[], 5, 2);
        let mut filter = VarFilter::new(MaxMaf(0.5)).unwrap();
        filter.filter_block(&mut block).unwrap();
        assert_eq!(block.num_vars, 0);
        assert!(block.gts.is_empty());
        assert_eq!(filter.stats(), FilteringStats::default());
    }

    /// The two variants of five individuals with every genotype `./.` of
    /// "What pyNei does that is odd": pyNei keeps neither at a threshold of
    /// 1, because a chunk in which nothing is called has no allele
    /// frequency, and popnei, which counts row by row, keeps neither for
    /// the same reason, that a variant with no called allele has no major
    /// allele frequency. The missing data filter keeps both at 1: their
    /// missing rate is 1.
    #[test]
    fn two_variants_with_every_genotype_missing_are_kept_by_no_maf_or_obs_het_filter() {
        let nothing_called = [MISSING_ALLELE; 10];
        let block_of_them = || block_of(&[(1, &nothing_called), (2, &nothing_called)], 5, 2);

        for criterion in [MaxMaf(1.0), MaxObsHet(1.0)] {
            let mut block = block_of_them();
            let mut filter = VarFilter::new(criterion).unwrap();
            filter.filter_block(&mut block).unwrap();
            assert_eq!(block.num_vars, 0, "{criterion:?}");
            assert!(positions_of(&block).is_empty(), "{criterion:?}");
            assert_eq!(
                filter.stats(),
                FilteringStats {
                    vars_processed: 2,
                    vars_kept: 0,
                },
                "{criterion:?}"
            );
        }

        let mut block = block_of_them();
        VarFilter::new(MaxMissingRate(1.0))
            .unwrap()
            .filter_block(&mut block)
            .unwrap();
        assert_eq!(positions_of(&block), [1, 2]);
    }

    /// 29 missing genotypes of 100 individuals pass a threshold of 0.29,
    /// because the number of the variant is one division of the two counts
    /// compared with the threshold, and not the threshold times the
    /// individuals: 0.29 * 100 is 28.999999999999996 and would drop the
    /// variant.
    #[test]
    fn a_variant_missing_in_29_of_100_individuals_passes_a_threshold_of_0_29() {
        let mut gts = vec![0_i8; 200];
        for allele in gts.iter_mut().take(58) {
            *allele = MISSING_ALLELE;
        }
        let mut block = block_of(&[(1, &gts)], 100, 2);
        let mut filter = VarFilter::new(MaxMissingRate(0.29)).unwrap();
        filter.filter_block(&mut block).unwrap();
        assert_eq!(positions_of(&block), [1]);
        assert_eq!(
            filter.stats(),
            FilteringStats {
                vars_processed: 1,
                vars_kept: 1,
            }
        );
        // The two ways of comparing, of which popnei and pyNei take the
        // first: the division keeps the variant and the product drops it.
        let missing = f64::from(29_u32);
        let individuals = f64::from(100_u32);
        let threshold = 0.29_f64;
        assert!(missing / individuals <= threshold);
        assert!(missing > threshold * individuals);
    }

    /// The ploidy of the filter is the one of the block: these three
    /// tetraploid genotypes, whose counts `docs/specs/variant.md` gives,
    /// are 2 called, 1 missing and 1 heterozygous, so their observed
    /// heterozygosity is 1/2; read as six diploid genotypes it would be
    /// 1/5, and the variant would pass a threshold of 0.3.
    #[test]
    fn the_number_of_a_variant_is_worked_out_with_the_ploidy_of_the_block() {
        // 0/0/0/1 1/1/1/1 0/./0/0
        let tetraploid = [0, 0, 0, 1, 1, 1, 1, 1, 0, -1, 0, 0];
        let kept = |threshold| {
            let mut block = block_of(&[(1, &tetraploid)], 3, 4);
            VarFilter::new(MaxObsHet(threshold))
                .unwrap()
                .filter_block(&mut block)
                .unwrap();
            block.num_vars
        };
        assert_eq!(kept(0.5), 1);
        assert_eq!(kept(0.3), 0);

        // Its missing rate is 1 genotype of 3 individuals, which stays at
        // 0.34 and goes at 0.33.
        let missing_rate_kept = |threshold| {
            let mut block = block_of(&[(1, &tetraploid)], 3, 4);
            VarFilter::new(MaxMissingRate(threshold))
                .unwrap()
                .filter_block(&mut block)
                .unwrap();
            block.num_vars
        };
        assert_eq!(missing_rate_kept(0.34), 1);
        assert_eq!(missing_rate_kept(0.33), 0);
    }

    /// The rows of a block are read on the threads of rayon natively and
    /// one after another in wasm, and the two give the same values in the
    /// same order.
    #[test]
    fn the_rows_read_on_the_threads_and_one_by_one_give_the_same_values() {
        let block = the_worked_example();
        for criterion in [MaxMissingRate(0.4), MaxMaf(0.88), MaxObsHet(0.25)] {
            let on_the_threads = keep_of_the_rows(criterion, &block.gts, 10, 5, 2).unwrap();
            let one_by_one = keep_of_the_rows_one_by_one(criterion, &block.gts, 10, 5, 2).unwrap();
            assert_eq!(on_the_threads, one_by_one, "{criterion:?}");
        }
    }

    /// The kind is the name the counts of the filter have for a Python and
    /// a TypeScript user, and it is the same for every threshold.
    #[test]
    fn the_kind_of_a_criterion_is_the_name_its_counts_have_in_python() {
        assert_eq!(MaxMissingRate(0.04).kind(), "missing_data");
        assert_eq!(MaxMaf(0.8).kind(), "maf");
        assert_eq!(MaxObsHet(0.5).kind(), "obs_het");
        assert_eq!(MaxMaf(0.0).kind(), MaxMaf(1.0).kind());
    }
}
