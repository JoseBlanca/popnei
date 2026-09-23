//! The distances between populations, and the counts of one variant that
//! every one of them is built from.
//!
//! Seven measures tell a user how far apart two populations are: Hudson's
//! F_ST, f_2, the chord distance, Nei's D_A, Jost's D, Nei's G_ST and the
//! standardized G''_ST. All of them are ratios of sums over the variants,
//! and all of the sums come from three counts of a population at one
//! variant: how often each allele was called there, how many genotypes were
//! called whole, and how many of those are heterozygous. So one pass over
//! the variants gives every measure, and what it keeps for a pair of
//! populations is six numbers: the sum of H_b, the chance that an allele
//! drawn from each of the two populations is not the same allele; the sum
//! of H_w, the mean of the two within population heterozygosities; the sum
//! over the alleles of the square root of the product of the two
//! frequencies, which the chord distance is built from; the sums of the
//! corrected H_S and H_T, the diversity within the two populations and over
//! the two pooled, which Jost's D and the two measures beside it read; and
//! how many variants counted for the pair. The divisions that turn the six
//! into a measure are made once, at the end, so neither the size of the
//! blocks nor the number of threads changes what comes out.
//!
//! [`PopVarCounts`] holds the three counts of one variant in one
//! population, and [`PopDistPerVar::of_var`] gives what one variant adds to
//! the sums of one pair, or nothing when the variant does not count for it.
//!
//! `docs/specs/dists.md` has the design, the formulas and the numbers the
//! tests assert, and the row `dists` of section 9 of
//! `docs/architecture.md` is where the module sits.

#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "nothing outside the tests reads the counts and the sums of one variant yet: what reads them is the pass over a reader, `calc_pop_dist_sums` of `docs/specs/dists.md`, which is written after them, and this expectation is an error the moment it is"
    )
)]

use crate::error::Result;
use crate::stats::{ObsHet, checked_ploidy, raised};
use crate::variant::{AlleleCounts, GtCounts, count_alleles_of, count_gts_of};

/// How many populations a pairwise measure is over, the s of the
/// corrections of H_S and H_T.
///
/// It is 2 because the calculation is of one pair at a time, whatever
/// number of populations the user gave, which is what `num_pops = 2` inside
/// `_calc_pairwise_dest` of pyNei does.
const NUM_POPS_OF_A_PAIR: f64 = 2.0;

/// The counts of one variant in one population: how often each allele was
/// called, how many alleles that is, and how many genotypes were called,
/// missing and heterozygous.
///
/// A half called genotype, `0/.` in a VCF, gives its called allele to the
/// counts of the alleles and is missing among the genotypes, as it is in
/// the expected heterozygosity of the `stats` module. So a population can
/// have called alleles at a variant and no called genotype.
///
/// One of these is kept for each population and filled again at each
/// variant, so that the pass allocates nothing per variant.
#[derive(Debug, Clone)]
pub(crate) struct PopVarCounts {
    /// How often the allele a was called in the population at the variant.
    allele_counts: AlleleCounts,
    /// One past the largest allele the population called, so that a pair
    /// walks the alleles the variant holds and not the 128 a count of them
    /// has room for.
    num_alleles: usize,
    /// The sum of `allele_counts`, the called alleles of the population at
    /// the variant, which is the n_P the frequencies are over.
    called_alleles: u32,
    /// The genotypes of the population at the variant: the called ones, the
    /// missing ones and the heterozygous ones.
    gts: GtCounts,
}

impl PopVarCounts {
    /// Counts of no variant at all, which [`PopVarCounts::count_the_var`]
    /// writes over.
    pub(crate) fn new() -> PopVarCounts {
        PopVarCounts {
            allele_counts: [0; 128],
            num_alleles: 0,
            called_alleles: 0,
            gts: GtCounts::default(),
        }
    }

    /// It counts one variant over the individuals of one population, over
    /// whatever it held before.
    ///
    /// `gts` is the genotypes of one variant, `ploidy` alleles for each
    /// individual of the reader, which
    /// [`VariantRef::gts`](crate::variant::VariantRef::gts) gives, and
    /// `individuals` the index of each individual of the population among
    /// them, which [`Pops::individuals`](crate::stats::Pops::individuals)
    /// gives.
    ///
    /// # Errors
    ///
    /// Those of [`count_alleles_of`] and [`count_gts_of`]: genotypes that
    /// are not a whole number of genotypes of the ploidy, a variant of more
    /// alleles than a count of them holds, an allele below the missing one,
    /// and an individual of the population beyond the variant.
    pub(crate) fn count_the_var(
        &mut self,
        gts: &[i8],
        ploidy: usize,
        individuals: &[usize],
    ) -> Result<()> {
        self.called_alleles = count_alleles_of(gts, ploidy, individuals, &mut self.allele_counts)?;
        self.gts = count_gts_of(gts, ploidy, individuals)?;
        self.num_alleles = self
            .allele_counts
            .iter()
            .rposition(|count| *count > 0)
            .map_or(0, |largest| largest.saturating_add(1));
        Ok(())
    }

    /// How often each allele was called, the allele a at the place a.
    pub(crate) fn allele_counts(&self) -> &AlleleCounts {
        &self.allele_counts
    }

    /// The called alleles of the population at the variant, the sum of
    /// [`PopVarCounts::allele_counts`].
    pub(crate) fn called_alleles(&self) -> u32 {
        self.called_alleles
    }

    /// The counts of the genotypes: the called, the missing and the
    /// heterozygous ones.
    pub(crate) fn gts(&self) -> GtCounts {
        self.gts
    }
}

/// What one variant adds to the sums of one pair of populations: five of
/// the six numbers the pass keeps for a pair, the sixth being how many
/// variants counted.
///
/// Every one of them is a value of that variant alone and not a ratio, so
/// the threads add them in any order and the measures divide once, when the
/// pass is over.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct VarSums {
    /// H_b, the chance that an allele drawn from one population and one
    /// drawn from the other are not the same allele.
    pub(crate) h_b: f64,
    /// H_w, the mean of the two within population heterozygosities, each
    /// corrected for being estimated from the copies it is computed over.
    pub(crate) h_w: f64,
    /// The sum over the alleles of the square root of the product of the
    /// two frequencies, which the chord distance is built from.
    pub(crate) sqrt_of_the_products: f64,
    /// H_S', the diversity within the two populations, corrected for the
    /// sample as Nei and Chesser do.
    pub(crate) corrected_h_s: f64,
    /// H_T', the diversity over the two populations pooled at equal weight,
    /// corrected the same way.
    pub(crate) corrected_h_t: f64,
}

impl VarSums {
    /// f_2 of this variant, H_b - H_w, which is negative where the two
    /// populations differ by no more than the noise of the sample.
    pub(crate) fn f2(&self) -> f64 {
        self.h_b - self.h_w
    }
}

/// What one variant of a pair of populations adds to the sums of that pair,
/// from the counts of the variant in each of them.
///
/// The frequency of the allele a in the population P is p_Pa = c_Pa / n_P,
/// with c_Pa how often a was called in P at the variant and n_P the called
/// alleles of P there. From the frequencies of the two populations A and B,
/// and with k the ploidy:
///
/// ```text
/// H_b  = 1 - sum over a of p_Aa p_Ba
/// H_w  = (u_A + u_B) / 2,  u_P = (n_P / (n_P - 1)) (1 - sum over a of p_Pa^2)
/// the chord sum = sum over a of sqrt(p_Aa p_Ba)
/// H_S  = (E_A + E_B) / 2,  E_P = 1 - sum over a of p_Pa^k
/// H_T  = 1 - sum over a of ((p_Aa + p_Ba) / 2)^k
/// ```
///
/// H_S and H_T are then corrected for the sample as Nei and Chesser (1983)
/// do, with n the harmonic mean of the called genotypes of the two
/// populations at that variant, H_obs the mean of their observed
/// heterozygosities and s = 2 populations:
///
/// ```text
/// H_S' = (n / (n - 1)) (H_S - H_obs / (2n))
/// H_T' = H_T + H_S' / (n s) - H_obs / (2 n s)
/// ```
#[derive(Debug, Clone, Copy)]
pub(crate) struct PopDistPerVar {
    /// k, the ploidy of the variants, which E_P and H_T raise the
    /// frequencies to.
    ploidy: u32,
    /// How many genotypes a population has to have called at a variant for
    /// that variant to count for a pair it is in, which is also what says
    /// whether the population has a value there at all.
    obs_het: ObsHet,
}

impl PopDistPerVar {
    /// It, with the ploidy of the variants and how many called genotypes a
    /// population needs at a variant for that variant to count for a pair
    /// it is in.
    ///
    /// # Errors
    ///
    /// A `ploidy` of 0 or above the largest ploidy a reader of popnei
    /// gives, 255.
    pub(crate) fn new(ploidy: usize, min_num_individuals: u32) -> Result<PopDistPerVar> {
        Ok(PopDistPerVar {
            ploidy: checked_ploidy("ploidy", ploidy)?,
            obs_het: ObsHet::new(min_num_individuals),
        })
    }

    /// What the variant adds to the sums of the pair, or `None` when it
    /// does not count for it.
    ///
    /// A variant counts for a pair when each of its two populations has
    /// called at least `min_num_individuals` genotypes at it, and one
    /// genotype at least whatever that threshold is: a population that has
    /// called nothing has no allele frequency there. Two cases beyond the
    /// count of the genotypes drop it, and neither can arise at a
    /// `min_num_individuals` of 2 or more. One called genotype in each
    /// population leaves the harmonic mean of the two at 1, which the
    /// correction of H_S divides by 1 - 1; and a population with one called
    /// allele, which takes a ploidy of 1, leaves u_P at 0 over 0. Both are
    /// in "Variants that do not count, populations with little data, and
    /// negative values" of `docs/specs/dists.md`.
    pub(crate) fn of_var(
        &self,
        of_one: &PopVarCounts,
        of_the_other: &PopVarCounts,
    ) -> Option<VarSums> {
        let obs_het_of_one = self.obs_het.of_var(of_one.gts)?;
        let obs_het_of_the_other = self.obs_het.of_var(of_the_other.gts)?;
        if of_one.gts.called == 1 && of_the_other.gts.called == 1 {
            return None;
        }
        if of_one.called_alleles < 2 || of_the_other.called_alleles < 2 {
            return None;
        }
        let n_of_one = f64::from(of_one.called_alleles);
        let n_of_the_other = f64::from(of_the_other.called_alleles);
        let over_the_alleles = self.over_the_alleles(of_one, of_the_other);
        let within_one = (n_of_one / (n_of_one - 1.0)) * (1.0 - over_the_alleles.two_alike_of_one);
        let within_the_other = (n_of_the_other / (n_of_the_other - 1.0))
            * (1.0 - over_the_alleles.two_alike_of_the_other);
        let h_s = ((1.0 - over_the_alleles.all_alike_of_one)
            + (1.0 - over_the_alleles.all_alike_of_the_other))
            / 2.0;
        let h_t = 1.0 - over_the_alleles.all_alike_pooled;
        // The harmonic mean of the called genotypes of the two
        // populations, which is pyNei's `hmean` over the two. It is above 1
        // unless both are 1, which the caller above has left out.
        let harmonic = NUM_POPS_OF_A_PAIR
            / (1.0 / f64::from(of_one.gts.called) + 1.0 / f64::from(of_the_other.gts.called));
        let obs_het = (obs_het_of_one + obs_het_of_the_other) / 2.0;
        let corrected_h_s = (harmonic / (harmonic - 1.0)) * (h_s - obs_het / (2.0 * harmonic));
        let corrected_h_t = h_t + corrected_h_s / (harmonic * NUM_POPS_OF_A_PAIR)
            - obs_het / (2.0 * harmonic * NUM_POPS_OF_A_PAIR);
        Some(VarSums {
            h_b: 1.0 - over_the_alleles.same_allele_in_both,
            h_w: (within_one + within_the_other) / 2.0,
            sqrt_of_the_products: over_the_alleles.sqrt_of_the_products,
            corrected_h_s,
            corrected_h_t,
        })
    }

    /// The six sums over the alleles of the variant that the values of the
    /// pair are built from.
    ///
    /// The alleles are added in the order of their number, so that two runs
    /// over the same variant give the same bits, and the walk stops at the
    /// largest allele either population called: the entries beyond it are 0
    /// in both and add nothing.
    fn over_the_alleles(
        &self,
        of_one: &PopVarCounts,
        of_the_other: &PopVarCounts,
    ) -> OverTheAlleles {
        let n_of_one = f64::from(of_one.called_alleles);
        let n_of_the_other = f64::from(of_the_other.called_alleles);
        let num_alleles = of_one.num_alleles.max(of_the_other.num_alleles);
        let mut sums = OverTheAlleles::default();
        for (count_of_one, count_of_the_other) in of_one
            .allele_counts
            .iter()
            .zip(&of_the_other.allele_counts)
            .take(num_alleles)
        {
            let p_of_one = f64::from(*count_of_one) / n_of_one;
            let p_of_the_other = f64::from(*count_of_the_other) / n_of_the_other;
            let product = p_of_one * p_of_the_other;
            sums.same_allele_in_both += product;
            sums.sqrt_of_the_products += product.sqrt();
            sums.two_alike_of_one += p_of_one * p_of_one;
            sums.two_alike_of_the_other += p_of_the_other * p_of_the_other;
            sums.all_alike_of_one += raised(p_of_one, self.ploidy);
            sums.all_alike_of_the_other += raised(p_of_the_other, self.ploidy);
            sums.all_alike_pooled += raised((p_of_one + p_of_the_other) / 2.0, self.ploidy);
        }
        sums
    }
}

/// The sums over the alleles of one variant that the values of one pair of
/// populations are built from, each of them added allele by allele.
#[derive(Debug, Clone, Copy, Default)]
struct OverTheAlleles {
    /// The chance that an allele drawn from each population is the same
    /// allele, the sum of p_Aa p_Ba, which H_b is taken from 1.
    same_allele_in_both: f64,
    /// The sum of the square roots of those products, which the chord
    /// distance is built from.
    sqrt_of_the_products: f64,
    /// The chance that two copies drawn from the first population with
    /// replacement are alike, the sum of p_Aa^2.
    two_alike_of_one: f64,
    /// The same for the second population.
    two_alike_of_the_other: f64,
    /// The chance that k copies drawn from the first population with
    /// replacement are all alike, the sum of p_Aa^k, which E_A is taken
    /// from 1.
    all_alike_of_one: f64,
    /// The same for the second population.
    all_alike_of_the_other: f64,
    /// The same over the two populations pooled at equal weight, the sum of
    /// ((p_Aa + p_Ba) / 2)^k, which H_T is taken from 1.
    all_alike_pooled: f64,
}

#[cfg(test)]
mod tests {
    use super::{PopDistPerVar, PopVarCounts, VarSums};

    /// The genotypes of the four variants of the worked example of "How it
    /// is verified" of `docs/specs/dists.md`: 6 diploid individuals, the
    /// two alleles of each after those of the one before, with -1 for an
    /// allele that was not called.
    ///
    /// The variants are, with pop1 the first three individuals and pop2 the
    /// last three: one biallelic variant where the populations differ; one
    /// of three alleles; one where both are fixed for the allele 0, with a
    /// missing genotype in pop1 and a half called one in pop2; and one that
    /// is the same in both, every genotype heterozygous.
    const WORKED_EXAMPLE: [[i8; 12]; 4] = [
        [0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1],
        [0, 1, 1, 2, 2, 2, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, -1],
        [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
    ];

    /// The individuals of the first population of the worked example.
    const POP1: [usize; 3] = [0, 1, 2];

    /// The individuals of its second population.
    const POP2: [usize; 3] = [3, 4, 5];

    /// The counts of one variant over one population.
    fn counts_of(gts: &[i8], ploidy: usize, individuals: &[usize]) -> PopVarCounts {
        let mut counts = PopVarCounts::new();
        match counts.count_the_var(gts, ploidy, individuals) {
            Ok(()) => counts,
            Err(error) => panic!("{error}"),
        }
    }

    /// What the variant of the worked example adds to the sums of its one
    /// pair, at the `min_num_individuals` given.
    fn sums_of_the_example(var: usize, min_num_individuals: u32) -> Option<VarSums> {
        let gts = &WORKED_EXAMPLE[var];
        let per_var = match PopDistPerVar::new(2, min_num_individuals) {
            Ok(per_var) => per_var,
            Err(error) => panic!("{error}"),
        };
        per_var.of_var(&counts_of(gts, 2, &POP1), &counts_of(gts, 2, &POP2))
    }

    /// The values of one variant, which the spec's table gives to six
    /// decimals.
    fn assert_the_value_is(found: f64, expected: f64, what: &str) {
        assert!(
            (found - expected).abs() < 1e-6,
            "{what} is {found} and not {expected}"
        );
    }

    /// The counts table of "How it is verified" of `docs/specs/dists.md`:
    /// the counts of the alleles 0, 1 and 2, the called alleles and the
    /// called genotypes of both populations at the four variants. The last
    /// line holds the half called genotype, which gives pop2 five called
    /// alleles from two called genotypes.
    #[test]
    fn the_worked_example_has_the_counts_of_the_spec() {
        let expected = [
            ([5, 1, 0], 6, 3, [1, 5, 0], 6, 3),
            ([1, 2, 3], 6, 3, [5, 1, 0], 6, 3),
            ([4, 0, 0], 4, 2, [5, 0, 0], 5, 2),
            ([3, 3, 0], 6, 3, [3, 3, 0], 6, 3),
        ];
        for (var, (of_pop1, n_1, called_1, of_pop2, n_2, called_2)) in
            expected.into_iter().enumerate()
        {
            let gts = &WORKED_EXAMPLE[var];
            let counts_1 = counts_of(gts, 2, &POP1);
            let counts_2 = counts_of(gts, 2, &POP2);
            assert_eq!(counts_1.allele_counts()[..3], of_pop1, "variant {var}");
            assert_eq!(counts_1.called_alleles(), n_1, "variant {var}");
            assert_eq!(counts_1.gts().called, called_1, "variant {var}");
            assert_eq!(counts_2.allele_counts()[..3], of_pop2, "variant {var}");
            assert_eq!(counts_2.called_alleles(), n_2, "variant {var}");
            assert_eq!(counts_2.gts().called, called_2, "variant {var}");
        }
    }

    /// The per variant table of the same section: H_b, H_w, f_2, the sum of
    /// the square roots of the products, H_S' and H_T' of the four
    /// variants, at a `min_num_individuals` of 1.
    #[test]
    fn the_worked_example_has_the_per_variant_values_of_the_spec() {
        let expected = [
            (0.722222, 0.333333, 0.388889, 0.745356, 0.333333, 0.527778),
            (0.805556, 0.533333, 0.272222, 0.608380, 0.541667, 0.673611),
            (0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
            (0.5, 0.6, -0.1, 1.0, 0.5, 0.5),
        ];
        for (var, (h_b, h_w, f2, sqrt_of_the_products, h_s, h_t)) in
            expected.into_iter().enumerate()
        {
            let Some(sums) = sums_of_the_example(var, 1) else {
                panic!("the variant {var} does not count for the pair");
            };
            assert_the_value_is(sums.h_b, h_b, &format!("H_b of the variant {var}"));
            assert_the_value_is(sums.h_w, h_w, &format!("H_w of the variant {var}"));
            assert_the_value_is(sums.f2(), f2, &format!("f_2 of the variant {var}"));
            assert_the_value_is(
                sums.sqrt_of_the_products,
                sqrt_of_the_products,
                &format!("the sum of the square roots of the variant {var}"),
            );
            assert_the_value_is(
                sums.corrected_h_s,
                h_s,
                &format!("H_S' of the variant {var}"),
            );
            assert_the_value_is(
                sums.corrected_h_t,
                h_t,
                &format!("H_T' of the variant {var}"),
            );
        }
    }

    /// The third variant of the worked example, where both populations are
    /// fixed for the allele 0: H_b is 0 and so are H_w, f_2 and the two
    /// corrected values, and the sum of the square roots is 1. It adds
    /// nothing to five of the six sums and one to the sixth, which is what
    /// a measure that is a ratio of sums is for.
    #[test]
    fn a_variant_both_pops_are_fixed_for_adds_nothing_but_a_count() {
        let Some(sums) = sums_of_the_example(2, 1) else {
            panic!("the variant does not count for the pair");
        };
        assert_the_value_is(sums.h_b, 0.0, "H_b");
        assert_the_value_is(sums.h_w, 0.0, "H_w");
        assert_the_value_is(sums.f2(), 0.0, "f_2");
        assert_the_value_is(
            sums.sqrt_of_the_products,
            1.0,
            "the sum of the square roots",
        );
        assert_the_value_is(sums.corrected_h_s, 0.0, "H_S'");
        assert_the_value_is(sums.corrected_h_t, 0.0, "H_T'");
    }

    /// The fourth variant, which is the same in both populations and every
    /// genotype of which is heterozygous: its f_2 is -0.1. popnei does not
    /// clamp it, because an estimator that could not go below zero would be
    /// biased upwards.
    #[test]
    fn a_variant_that_is_the_same_in_both_pops_has_a_negative_f2() {
        let Some(sums) = sums_of_the_example(3, 1) else {
            panic!("the variant does not count for the pair");
        };
        assert_the_value_is(sums.f2(), -0.1, "f_2");
        assert!(sums.f2() < 0.0, "f_2 is {} and not below 0", sums.f2());
    }

    /// The third variant has two called genotypes in each population, so a
    /// threshold of 3 drops it and a threshold of 2 keeps it. The first
    /// variant, with three called genotypes in each, is kept at both.
    #[test]
    fn a_variant_a_pop_has_too_few_called_genotypes_at_does_not_count() {
        assert!(sums_of_the_example(2, 3).is_none());
        assert!(sums_of_the_example(2, 2).is_some());
        assert!(sums_of_the_example(0, 3).is_some());
    }

    /// A population with no called genotype at a variant has no allele
    /// frequency there, so the variant counts for no pair that population
    /// is in, whatever the threshold. The population here holds the one
    /// individual of the third variant whose genotype is missing whole.
    #[test]
    fn a_pop_with_no_called_genotype_makes_the_variant_not_count() {
        let gts = &WORKED_EXAMPLE[2];
        let per_var = match PopDistPerVar::new(2, 0) {
            Ok(per_var) => per_var,
            Err(error) => panic!("{error}"),
        };
        let nothing_called = counts_of(gts, 2, &[2]);
        assert_eq!(nothing_called.called_alleles(), 0);
        assert_eq!(nothing_called.gts().called, 0);
        assert!(
            per_var
                .of_var(&nothing_called, &counts_of(gts, 2, &POP2))
                .is_none()
        );
    }

    /// One called genotype in each population leaves the harmonic mean of
    /// the two at 1, which the correction of H_S divides by 1 - 1, so the
    /// variant does not count for that pair. It takes a
    /// `min_num_individuals` of 1, and at 2 the same variant is dropped by
    /// the threshold instead. The variant is the first of the worked
    /// example over one individual of each population.
    #[test]
    fn a_variant_with_one_called_genotype_in_each_pop_does_not_count() {
        let gts = &WORKED_EXAMPLE[0];
        let per_var = match PopDistPerVar::new(2, 1) {
            Ok(per_var) => per_var,
            Err(error) => panic!("{error}"),
        };
        let of_one = counts_of(gts, 2, &[1]);
        let of_the_other = counts_of(gts, 2, &[4]);
        assert_eq!(of_one.gts().called, 1);
        assert_eq!(of_the_other.gts().called, 1);
        assert!(per_var.of_var(&of_one, &of_the_other).is_none());
        // One called genotype against three still counts.
        assert!(per_var.of_var(&of_one, &counts_of(gts, 2, &POP2)).is_some());
    }

    /// A haploid population with one called genotype has one called allele,
    /// and its within population heterozygosity is 0 over 0: the variant
    /// does not count for the pairs that population is in. Two called
    /// genotypes are enough, so the same variant counts for the pair of the
    /// two populations of two.
    #[test]
    fn a_haploid_pop_with_one_called_allele_makes_the_variant_not_count() {
        let gts = [0_i8, 1, 0, 1, 1];
        let per_var = match PopDistPerVar::new(1, 1) {
            Ok(per_var) => per_var,
            Err(error) => panic!("{error}"),
        };
        let of_one_genotype = counts_of(&gts, 1, &[0]);
        let of_two_genotypes = counts_of(&gts, 1, &[1, 2]);
        let of_two_more = counts_of(&gts, 1, &[3, 4]);
        assert_eq!(of_one_genotype.called_alleles(), 1);
        assert!(
            per_var
                .of_var(&of_one_genotype, &of_two_genotypes)
                .is_none()
        );
        assert!(
            per_var
                .of_var(&of_two_genotypes, &of_one_genotype)
                .is_none()
        );
        assert!(per_var.of_var(&of_two_genotypes, &of_two_more).is_some());
    }

    /// A ploidy of 0 is not one the frequencies can be raised to, and 256
    /// is above the largest a reader of popnei gives.
    #[test]
    fn a_ploidy_of_zero_or_above_the_largest_one_is_an_error() {
        assert!(PopDistPerVar::new(0, 20).is_err());
        assert!(PopDistPerVar::new(256, 20).is_err());
    }
}
