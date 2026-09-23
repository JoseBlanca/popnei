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
//! [`JackknifeGroups`] says how the variants are cut into the groups that
//! the standard errors are resampled over, and [`JackknifeWalk`] is the
//! walk that cuts them: a stretch of one chromosome anchored on its own
//! first variant, one variant, or no group at all.
//!
//! [`PopDistSums`] is what a pass leaves: the six sums of every pair within
//! every group. [`PopDistSums::measure`] turns them into one of the seven
//! measures, which [`PopDistMeasure`] names, and
//! [`PopDistSums::standard_error`] into the standard error of that measure,
//! by leaving each group out in turn and asking how far the measure moves,
//! which is the delete-m jackknife of Busing, Meijer and van der Leeden.
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

use crate::dists::{index_of_the_pair, num_pairs_of};
use crate::error::{Error, Result};
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

/// How the variants are cut into the groups that the standard errors are
/// resampled over.
///
/// A distance between two populations is a mean over variants, and variants
/// near each other on a chromosome carry much the same history, so treating
/// each of them as an independent draw makes the error look smaller than it
/// is. The standard error is built instead by leaving each group out in
/// turn, which asks for groups long enough that two of them are nearly
/// independent. The literature calls a group a block and the method the
/// block jackknife; popnei says group, because a block here is the run of
/// variants a reader gives.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JackknifeGroups {
    /// No standard errors, and the pass does not ask its reader for the
    /// positions.
    None,
    /// Each variant its own group, which is what a few hundred
    /// microsatellite loci scattered over a genome want and what a panel of
    /// linked SNPs must not use.
    PerVariant,
    /// Stretches of one chromosome this many base pairs long, each one
    /// anchored on its own first variant. 1 at least.
    OfBasePairs(u64),
}

/// One group of variants: the chromosome, and the first and the last
/// position it holds, both included.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GroupId {
    /// Its chromosome, as a number in the
    /// [`ChromTable`](crate::variant::ChromTable) of the reader of the
    /// pass, which gives the name.
    pub chrom: u32,
    /// The position of its first variant, 1 based as in a VCF.
    pub start: u64,
    /// The position of its last variant.
    pub end: u64,
}

/// The walk over the variants that cuts them into the resampling groups, in
/// the order the reader gives them.
///
/// [`JackknifeWalk::group_of`] takes the chromosome and the position of one
/// variant and gives the group it falls in, starting a new one where the
/// cut falls; [`JackknifeWalk::groups`] gives the groups in the order they
/// were started, which is the order of the result.
///
/// With [`JackknifeGroups::OfBasePairs`] a new group starts at the first
/// variant of a chromosome and at the first variant whose position is the
/// length or more beyond the first variant of the group being filled. So a
/// group is anchored on its own first variant and not on a grid of
/// multiples of the length: cutting at the multiples instead leaves a group
/// of one variant wherever a variant sits just past a multiple, and one
/// such group is enough to move the standard error. On the biallelic panel
/// of `tests/reference/dists/panel.vcf.gz` cut into groups of 100 000 base
/// pairs, the anchored rule gives 12 groups of 100 variants and a standard
/// error of f_2 for p0 and p1 of 0.00180, and the multiples give 14 groups,
/// two of them of one variant, and 0.00242. The anchored rule is the one
/// ADMIXTOOLS 2 uses, which is what lets the numbers be compared.
pub(crate) struct JackknifeWalk {
    /// How the variants are being cut.
    how: JackknifeGroups,
    /// The groups in the order they were started, the one being filled
    /// last.
    groups: Vec<GroupId>,
}

impl JackknifeWalk {
    /// The walk that cuts the variants the way `how` says.
    ///
    /// # Errors
    ///
    /// A length of 0 base pairs, which is no stretch of a chromosome.
    pub(crate) fn new(how: JackknifeGroups) -> Result<JackknifeWalk> {
        if how == JackknifeGroups::OfBasePairs(0) {
            return Err(Error::JackknifeGroupOfNoBasePairs);
        }
        Ok(JackknifeWalk {
            how,
            groups: Vec::new(),
        })
    }

    /// The group the variant at `chrom` and `pos` falls in, counted from 0,
    /// and `None` when no groups were asked for.
    ///
    /// The variants are given to it in the order the reader gives them, and
    /// their positions along a chromosome go up in every source popnei
    /// reads. A group that is being filled takes the variant when it is of
    /// its chromosome and its position is less than the length beyond the
    /// first variant of the group; otherwise the variant starts a group of
    /// its own.
    pub(crate) fn group_of(&mut self, chrom: u32, pos: u64) -> Option<usize> {
        match self.how {
            JackknifeGroups::None => return None,
            JackknifeGroups::PerVariant => {}
            JackknifeGroups::OfBasePairs(length) => {
                if let Some(filling) = self.groups.last_mut()
                    && filling.chrom == chrom
                    && pos.saturating_sub(filling.start) < length
                {
                    filling.end = pos;
                    return self.groups.len().checked_sub(1);
                }
            }
        }
        let at = self.groups.len();
        self.groups.push(GroupId {
            chrom,
            start: pos,
            end: pos,
        });
        Some(at)
    }

    /// The groups the variants walked so far were cut into, in the order
    /// they were started.
    pub(crate) fn groups(&self) -> &[GroupId] {
        &self.groups
    }
}

/// Which of the seven measures of how far apart two populations are a
/// caller asks a pass for.
///
/// All seven are ratios of the same six sums, so one pass over the variants
/// gives every one of them and a caller pays nothing for asking for a
/// second. `docs/specs/dists.md` has an item for each, with what it
/// answers and the program it is verified against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PopDistMeasure {
    /// Hudson's F_ST, how much of the diversity of the two populations
    /// taken together lies between them rather than within them.
    Fst,
    /// f_2, how much allele frequency the two populations have drifted
    /// apart by, in the units it was measured in.
    F2,
    /// The chord distance of Cavalli-Sforza and Edwards.
    Chord,
    /// Nei's D_A, the square of the chord distance.
    Da,
    /// Jost's D, the differentiation of the alleles the two hold.
    Dest,
    /// Nei's G_ST.
    Gst,
    /// G_ST rescaled to reach 1, the standardized G''_ST.
    GstStandardized,
}

/// The six sums of one pair of populations over a set of variants: what
/// every measure of that pair is a ratio of.
///
/// The five sums are of values of one variant, which
/// [`PopDistPerVar::of_var`] gives, and the sixth is how many variants
/// counted for the pair. None of them is a ratio, so the threads and the
/// blocks add them in any order and the divisions happen once, when the
/// pass is over.
///
/// The count is a `u32`, which is how many variants of one pair the pass
/// counts; the largest dataset of `docs/objectives.md` is a million
/// variants.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct PairSums {
    /// The sum of H_b over the variants that counted.
    pub(crate) h_b: f64,
    /// The sum of H_w over them.
    pub(crate) h_w: f64,
    /// The sum of the square roots of the products of the frequencies over
    /// them.
    pub(crate) sqrt_of_the_products: f64,
    /// The sum of the corrected H_S over them.
    pub(crate) corrected_h_s: f64,
    /// The sum of the corrected H_T over them.
    pub(crate) corrected_h_t: f64,
    /// How many variants counted for the pair.
    pub(crate) num_vars: u32,
}

impl PairSums {
    /// These sums with those of `of_the_group` taken out of them, which is
    /// what a measure with that group left out is calculated from.
    ///
    /// `None` when the group holds more variants than these sums do, which
    /// no group of a pass does since these are the sum of every group.
    fn without(&self, of_the_group: &PairSums) -> Option<PairSums> {
        Some(PairSums {
            h_b: self.h_b - of_the_group.h_b,
            h_w: self.h_w - of_the_group.h_w,
            sqrt_of_the_products: self.sqrt_of_the_products - of_the_group.sqrt_of_the_products,
            corrected_h_s: self.corrected_h_s - of_the_group.corrected_h_s,
            corrected_h_t: self.corrected_h_t - of_the_group.corrected_h_t,
            num_vars: self.num_vars.checked_sub(of_the_group.num_vars)?,
        })
    }
}

/// The measure of a pair out of its six sums, and `None` when no variant
/// counted for the pair.
///
/// With H_b, H_w and n the sums of H_b and H_w over the variants that
/// counted for the pair and how many those were,
///
/// ```text
/// F_ST = (H_b - H_w) / H_b
/// f_2  = (H_b - H_w) / n
/// ```
///
/// Both are ratios of sums and not means of per variant ratios, which is
/// what Bhatia, Patterson, Sankararaman and Price (2013) ask for and what
/// plink2 and ADMIXTOOLS 2 compute: the mean of the ratios gives weight to
/// the variants whose denominator is near zero, and it is undefined at a
/// variant both populations are fixed for the same allele at. Both can come
/// out negative, and popnei does not clamp them.
///
/// A pair whose sum of H_b is 0, which takes two populations fixed for the
/// same allele at every variant that counted for them, has a F_ST of 0 over
/// 0, a NaN. "The Rust interface" of the spec lists the cases with no value
/// and that is not one of them, and a caller sees the same thing either
/// way: the binding crates write a NaN for a `None` as well.
///
/// The five measures below them are the work packages 2 and 3 of
/// `docs/plans/dists-pops.md` and have no value until those are written.
fn value_of(measure: PopDistMeasure, sums: &PairSums) -> Option<f64> {
    if sums.num_vars == 0 {
        return None;
    }
    let between_minus_within = sums.h_b - sums.h_w;
    match measure {
        PopDistMeasure::Fst => Some(between_minus_within / sums.h_b),
        PopDistMeasure::F2 => Some(between_minus_within / f64::from(sums.num_vars)),
        PopDistMeasure::Chord
        | PopDistMeasure::Da
        | PopDistMeasure::Dest
        | PopDistMeasure::Gst
        | PopDistMeasure::GstStandardized => None,
    }
}

/// What one pass over the variants gives: the six sums of every pair of
/// populations within every resampling group, which every measure and every
/// standard error is worked out from.
///
/// The pairs are in the order of the distance vector, (0, 1), (0, 2), ...,
/// (1, 2), ..., over the populations in the order the caller named them.
/// [`PopDistSums::measure`] gives one of them a measure,
/// [`PopDistSums::standard_error`] the jackknife standard error of that
/// measure, and [`PopDistSums::f2_of_group`] the f_2 of one pair within one
/// group, which the f_3 and f_4 of a later spec are built from.
#[derive(Debug)]
pub struct PopDistSums {
    /// How many populations the pairs are of.
    num_pops: usize,
    /// The variants the pass was given, counted for a pair or not.
    num_vars: u64,
    /// The resampling groups in the order they were started, which is the
    /// order of the sums below, and no group at all when the caller asked
    /// for no standard errors.
    groups: Vec<GroupId>,
    /// The six sums of each pair within each group: every pair of the first
    /// group, then every pair of the second, so that the f_2 of a group is
    /// a run of this. When the caller asked for no groups it holds one run
    /// of the pairs, over every variant of the pass.
    of_each_group: Vec<PairSums>,
}

impl PopDistSums {
    /// The sums a pass built: the six numbers of each pair within each
    /// group, `groups.len()` runs of one for each pair, in the order of the
    /// distance vector.
    ///
    /// When `groups` is empty, which is what a caller who asked for no
    /// standard errors gets, `of_each_group` is one run of the pairs and
    /// holds every variant the pass counted.
    pub(crate) fn of_the_pass(
        num_pops: usize,
        num_vars: u64,
        groups: Vec<GroupId>,
        of_each_group: Vec<PairSums>,
    ) -> PopDistSums {
        PopDistSums {
            num_pops,
            num_vars,
            groups,
            of_each_group,
        }
    }

    /// How many populations the pairs are of.
    #[must_use]
    pub fn num_pops(&self) -> usize {
        self.num_pops
    }

    /// The variants the pass was given, counted for a pair or not.
    #[must_use]
    pub fn num_vars(&self) -> u64 {
        self.num_vars
    }

    /// The resampling groups, in the order they were started, and none when
    /// the caller asked for no standard errors.
    #[must_use]
    pub fn groups(&self) -> &[GroupId] {
        &self.groups
    }

    /// The variants that counted for the pair: both its populations had at
    /// least `min_num_individuals` called genotypes at them. `None` when i
    /// == j or when either is not a population.
    #[must_use]
    pub fn num_vars_of(&self, i: usize, j: usize) -> Option<u64> {
        Some(u64::from(
            self.total_of(self.index_of_the_pair(i, j)?)?.num_vars,
        ))
    }

    /// The measure for the pair. `None` where
    /// [`num_vars_of`](PopDistSums::num_vars_of) is 0 or `None`, and for
    /// the five measures the work packages 2 and 3 of
    /// `docs/plans/dists-pops.md` add.
    #[must_use]
    pub fn measure(&self, measure: PopDistMeasure, i: usize, j: usize) -> Option<f64> {
        self.of_the_pair(measure, self.index_of_the_pair(i, j)?)
    }

    /// The measure of every pair, in the order of the distance vector. The
    /// binding crates write NaN for a `None`.
    pub fn measures(&self, measure: PopDistMeasure) -> impl Iterator<Item = Option<f64>> + '_ {
        (0..num_pairs_of(self.num_pops).unwrap_or(0))
            .map(move |pair| self.of_the_pair(measure, pair))
    }

    /// f_2 within one group, which f_3 and f_4 are built from later.
    ///
    /// `None` when no variant of the pair fell in the group, where
    /// [`measure`](PopDistSums::measure) gives `None`, and when `group` is
    /// not one of [`groups`](PopDistSums::groups).
    #[must_use]
    pub fn f2_of_group(&self, group: usize, i: usize, j: usize) -> Option<f64> {
        if group >= self.groups.len() {
            return None;
        }
        let of_the_group = self.of_the_group(group, self.index_of_the_pair(i, j)?)?;
        value_of(PopDistMeasure::F2, of_the_group)
    }

    /// The jackknife standard error of the measure of the pair. `None`
    /// where [`measure`](PopDistSums::measure) is `None`, where no groups
    /// were asked for, and where every variant of the pair fell in one
    /// group.
    ///
    /// It is the delete-m jackknife for unequal m of Busing, Meijer and van
    /// der Leeden (1999, Statistics and Computing 9: 3, DOI
    /// 10.1023/A:1008800423698), which is what the f-statistics literature
    /// uses and which takes groups of different numbers of variants. With t
    /// the measure over every variant that counted for the pair, n how many
    /// those were, g the groups that hold at least one of them, m_j how
    /// many are in the group j, t_(j) the measure with that group left out
    /// and h_j = n / m_j,
    ///
    /// ```text
    /// u_j = h_j t - (h_j - 1) t_(j)
    /// t_J = sum over j of u_j / h_j
    /// v   = (1 / g) sum over j of (u_j - t_J)^2 / (h_j - 1)
    /// ```
    ///
    /// and this is the square root of v. What the result gives as the
    /// distance is t and not the jackknife estimate t_J, which is here to
    /// build the variance from. h_j - 1 is zero only for a group that holds
    /// every variant of the pair, which is the g = 1 case this gives `None`
    /// for, so nothing divides by zero.
    #[must_use]
    pub fn standard_error(&self, measure: PopDistMeasure, i: usize, j: usize) -> Option<f64> {
        let pair = self.index_of_the_pair(i, j)?;
        let over_all = self.total_of(pair)?;
        let over_all_value = value_of(measure, &over_all)?;
        let mut num_groups: u32 = 0;
        let mut jackknife_estimate = 0.0;
        for group in 0..self.groups.len() {
            let Some(pseudo) =
                self.pseudo_value_of(measure, pair, group, over_all_value, &over_all)
            else {
                continue;
            };
            num_groups = num_groups.checked_add(1)?;
            jackknife_estimate += pseudo.value / pseudo.weight;
        }
        if num_groups < 2 {
            return None;
        }
        let mut variance = 0.0;
        for group in 0..self.groups.len() {
            let Some(pseudo) =
                self.pseudo_value_of(measure, pair, group, over_all_value, &over_all)
            else {
                continue;
            };
            let from_the_estimate = pseudo.value - jackknife_estimate;
            variance += from_the_estimate * from_the_estimate / (pseudo.weight - 1.0);
        }
        Some((variance / f64::from(num_groups)).sqrt())
    }

    /// The pseudo-value of one group and its weight h_j, and `None` when no
    /// variant of the pair fell in the group or when the group holds every
    /// variant of it, which leaves no measure to compare with.
    fn pseudo_value_of(
        &self,
        measure: PopDistMeasure,
        pair: usize,
        group: usize,
        over_all_value: f64,
        over_all: &PairSums,
    ) -> Option<PseudoValue> {
        let of_the_group = self.of_the_group(group, pair)?;
        if of_the_group.num_vars == 0 {
            return None;
        }
        let weight = f64::from(over_all.num_vars) / f64::from(of_the_group.num_vars);
        let without_the_group = value_of(measure, &over_all.without(of_the_group)?)?;
        Some(PseudoValue {
            value: weight * over_all_value - (weight - 1.0) * without_the_group,
            weight,
        })
    }

    /// The measure of the pair at `pair` of the distance vector.
    fn of_the_pair(&self, measure: PopDistMeasure, pair: usize) -> Option<f64> {
        value_of(measure, &self.total_of(pair)?)
    }

    /// The six sums of the pair over every variant that counted for it, its
    /// groups added in the order they were started so that the total does
    /// not depend on how the pass was cut into blocks or threads.
    ///
    /// `None` when the populations make no pair and when the pair counted
    /// more variants than a `u32` holds.
    fn total_of(&self, pair: usize) -> Option<PairSums> {
        let num_pairs = num_pairs_of(self.num_pops)?;
        if num_pairs == 0 {
            return None;
        }
        let mut total = PairSums::default();
        for of_the_group in self.of_each_group.iter().skip(pair).step_by(num_pairs) {
            total.h_b += of_the_group.h_b;
            total.h_w += of_the_group.h_w;
            total.sqrt_of_the_products += of_the_group.sqrt_of_the_products;
            total.corrected_h_s += of_the_group.corrected_h_s;
            total.corrected_h_t += of_the_group.corrected_h_t;
            total.num_vars = total.num_vars.checked_add(of_the_group.num_vars)?;
        }
        Some(total)
    }

    /// The six sums of the pair at `pair` within the group at `group`.
    fn of_the_group(&self, group: usize, pair: usize) -> Option<&PairSums> {
        let num_pairs = num_pairs_of(self.num_pops)?;
        self.of_each_group
            .get(group.checked_mul(num_pairs)?.checked_add(pair)?)
    }

    /// Where the pair of the populations i and j is in the distance vector,
    /// in either order, and `None` when the two are one population or
    /// either of them is not one.
    fn index_of_the_pair(&self, i: usize, j: usize) -> Option<usize> {
        index_of_the_pair(self.num_pops, i, j)
    }
}

/// The pseudo-value of one group, u_j, with the weight h_j the variance
/// divides it by.
#[derive(Debug, Clone, Copy)]
struct PseudoValue {
    /// u_j, the measure of the pair as the group sees it.
    value: f64,
    /// h_j, how many of the variants of the pair the group holds one of.
    weight: f64,
}

#[cfg(test)]
mod tests {
    use std::path::{Path, PathBuf};

    use super::{
        GroupId, JackknifeGroups, JackknifeWalk, PairSums, PopDistMeasure, PopDistPerVar,
        PopDistSums, PopVarCounts, VarSums,
    };
    use crate::block::BlockReader;
    use crate::error::Error;
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::stats::Pops;
    use crate::variant::Needs;

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

    /// The groups a walk cuts `vars` into, each variant given as its
    /// chromosome and its position, with the group each variant fell in.
    fn walk_over(how: JackknifeGroups, vars: &[(u32, u64)]) -> (Vec<Option<usize>>, Vec<GroupId>) {
        let mut walk = match JackknifeWalk::new(how) {
            Ok(walk) => walk,
            Err(error) => panic!("{error}"),
        };
        let of_each_var = vars
            .iter()
            .map(|(chrom, pos)| walk.group_of(*chrom, *pos))
            .collect();
        (of_each_var, walk.groups().to_vec())
    }

    /// The panels and what plink2, R and ADMIXTOOLS 2 give for them live at
    /// the root of the repository, beside the Python tests that read the
    /// same files, and not inside this crate. `name` is the path of one of
    /// them under `tests/reference/`.
    fn reference(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference")
            .join(name)
    }

    /// The 12 groups that a length of 100 000 base pairs cuts the biallelic
    /// panel of 1200 variants into, over its two chromosomes, which is what
    /// "The standard errors" of `docs/specs/dists.md` says ADMIXTOOLS 2
    /// cuts it into: 100 variants each, and a new group at the second
    /// chromosome although its first position is within 100 000 base pairs
    /// of the first variant of the group the first chromosome ended in.
    #[test]
    fn the_panel_is_cut_into_twelve_groups_of_a_hundred_variants() {
        let expected = [
            ("chr1", 1000, 100_000),
            ("chr1", 101_000, 200_000),
            ("chr1", 201_000, 300_000),
            ("chr1", 301_000, 400_000),
            ("chr1", 401_000, 500_000),
            ("chr1", 501_000, 600_000),
            ("chr2", 1000, 100_000),
            ("chr2", 101_000, 200_000),
            ("chr2", 201_000, 300_000),
            ("chr2", 301_000, 400_000),
            ("chr2", 401_000, 500_000),
            ("chr2", 501_000, 600_000),
        ];
        let options = VcfOptions {
            ploidy: 2,
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::from_path(&reference("dists/panel.vcf.gz"), options) {
            Ok(reader) => reader,
            Err(error) => panic!("{error}"),
        };
        reader.set_needs(Needs::GTS | Needs::CHROM_POS);
        let mut walk = match JackknifeWalk::new(JackknifeGroups::OfBasePairs(100_000)) {
            Ok(walk) => walk,
            Err(error) => panic!("{error}"),
        };
        let mut num_vars_of_each_group: Vec<u64> = Vec::new();
        loop {
            let block = match reader.next_block() {
                Ok(Some(block)) => block,
                Ok(None) => break,
                Err(error) => panic!("{error}"),
            };
            for var in block.variants() {
                let (Some(chrom), Some(pos)) = (var.chrom(), var.pos()) else {
                    panic!("the reader gave a variant with no chromosome or no position");
                };
                let Some(at) = walk.group_of(chrom, pos) else {
                    panic!("the variant at {chrom} {pos} fell in no group");
                };
                if at == num_vars_of_each_group.len() {
                    num_vars_of_each_group.push(0);
                }
                num_vars_of_each_group[at] = num_vars_of_each_group[at].saturating_add(1);
            }
        }

        assert_eq!(walk.groups().len(), 12);
        assert_eq!(num_vars_of_each_group, vec![100_u64; 12]);
        for (at, (chrom, start, end)) in expected.into_iter().enumerate() {
            let group = walk.groups()[at];
            assert_eq!(reader.chroms().name(group.chrom), Some(chrom), "group {at}");
            assert_eq!((group.start, group.end), (start, end), "group {at}");
        }
    }

    /// The first variant of a chromosome starts a group, however near its
    /// position is to the first variant of the group being filled: the two
    /// positions are of different chromosomes and the distance between them
    /// means nothing.
    #[test]
    fn a_new_group_starts_at_each_chromosome() {
        let vars = [(0, 1000), (0, 1500), (1, 1600), (1, 2000)];
        let (of_each_var, groups) = walk_over(JackknifeGroups::OfBasePairs(100_000), &vars);

        assert_eq!(of_each_var, [Some(0), Some(0), Some(1), Some(1)]);
        assert_eq!(
            groups,
            [
                GroupId {
                    chrom: 0,
                    start: 1000,
                    end: 1500
                },
                GroupId {
                    chrom: 1,
                    start: 1600,
                    end: 2000
                },
            ]
        );
    }

    /// A group is anchored on its own first variant and not on a grid of
    /// multiples of the length. The three variants here are 99 999,
    /// 100 001 and 199 999 base pairs along one chromosome, cut into groups
    /// of 100 000: anchored, the first two are one group, because 100 001
    /// is 2 base pairs beyond the first of them, and the third starts
    /// another, because it is 100 000 beyond it. Cutting at the multiples
    /// of 100 000 would leave 99 999 in a group of its own and put the
    /// other two together.
    #[test]
    fn a_group_is_anchored_on_its_own_first_variant() {
        let vars = [(0, 99_999), (0, 100_001), (0, 199_999)];
        let (of_each_var, groups) = walk_over(JackknifeGroups::OfBasePairs(100_000), &vars);

        assert_eq!(of_each_var, [Some(0), Some(0), Some(1)]);
        assert_eq!(
            groups,
            [
                GroupId {
                    chrom: 0,
                    start: 99_999,
                    end: 100_001
                },
                GroupId {
                    chrom: 0,
                    start: 199_999,
                    end: 199_999
                },
            ]
        );
    }

    /// The shortest group there is holds the variants of one position and
    /// nothing else, which says that the length is read and not a number
    /// the walk keeps for itself.
    #[test]
    fn a_length_of_one_base_pair_groups_the_variants_of_one_position() {
        let vars = [(0, 5), (0, 5), (0, 6)];
        let (of_each_var, groups) = walk_over(JackknifeGroups::OfBasePairs(1), &vars);

        assert_eq!(of_each_var, [Some(0), Some(0), Some(1)]);
        assert_eq!(groups.len(), 2);
    }

    /// Each variant its own group, which is what a panel of scattered
    /// microsatellite loci asks for: a group of one variant holds one
    /// position, its start and its end.
    #[test]
    fn per_variant_gives_one_group_for_each_variant() {
        let vars = [(0, 10), (0, 20), (1, 5)];
        let (of_each_var, groups) = walk_over(JackknifeGroups::PerVariant, &vars);

        assert_eq!(of_each_var, [Some(0), Some(1), Some(2)]);
        assert_eq!(
            groups,
            [
                GroupId {
                    chrom: 0,
                    start: 10,
                    end: 10
                },
                GroupId {
                    chrom: 0,
                    start: 20,
                    end: 20
                },
                GroupId {
                    chrom: 1,
                    start: 5,
                    end: 5
                },
            ]
        );
    }

    /// No groups at all, which is what a user who asks for no standard
    /// error gets: every variant falls in none and there is nothing to
    /// resample over.
    #[test]
    fn no_groups_leaves_every_variant_in_none() {
        let vars = [(0, 10), (0, 20), (1, 5)];
        let (of_each_var, groups) = walk_over(JackknifeGroups::None, &vars);

        assert_eq!(of_each_var, [None, None, None]);
        assert!(groups.is_empty());
    }

    /// A group of 0 base pairs is no stretch of a chromosome: a caller who
    /// wants each variant in a group of its own asks for that instead. One
    /// base pair is a length the walk takes.
    #[test]
    fn a_length_of_no_base_pairs_is_an_error() {
        let refused = JackknifeWalk::new(JackknifeGroups::OfBasePairs(0));

        assert!(matches!(refused, Err(Error::JackknifeGroupOfNoBasePairs)));
        assert!(JackknifeWalk::new(JackknifeGroups::OfBasePairs(1)).is_ok());
    }

    /// A ploidy of 0 is not one the frequencies can be raised to, and 256
    /// is above the largest a reader of popnei gives.
    #[test]
    fn a_ploidy_of_zero_or_above_the_largest_one_is_an_error() {
        assert!(PopDistPerVar::new(0, 20).is_err());
        assert!(PopDistPerVar::new(256, 20).is_err());
    }

    /// The populations of one of the files of `tests/reference/`, which
    /// hold a header and then the name of an individual and the name of its
    /// population, and in the order `wanted` names them.
    fn pops_of_the_file(name: &str, individuals: &[String], wanted: &[&str]) -> Pops {
        let text = std::fs::read_to_string(reference(name)).expect("the populations of the file");
        let mut named: Vec<(String, Vec<String>)> = wanted
            .iter()
            .map(|pop| ((*pop).to_owned(), Vec::new()))
            .collect();
        for line in text.lines().skip(1) {
            let mut columns = line.split('\t');
            let (Some(individual), Some(pop)) = (columns.next(), columns.next()) else {
                panic!("the line `{line}` of {name} is not an individual and a population");
            };
            if let Some((_, of_the_pop)) = named.iter_mut().find(|(name, _)| name == pop) {
                of_the_pop.push(individual.to_owned());
            }
        }
        Pops::from_names(&named, individuals).expect("the populations of the file")
    }

    /// The six sums of every pair of `pops` over the variants of the panel
    /// `name`, cut into the groups `how` says.
    ///
    /// It is the pass of `calc_pop_dist_sums`, written here over one block
    /// at a time and one variant at a time, because the measures of a whole
    /// panel cannot be checked before there are sums of a whole panel and
    /// that function is written after this. `at_most_vars` stops it early,
    /// which is how the F_ST of one variant is seen.
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "the counts of the variants and of the pairs of a panel of the repository, far below what a usize holds"
    )]
    fn sums_of_the_panel(
        name: &str,
        pops_file: &str,
        how: JackknifeGroups,
        min_num_individuals: u32,
        at_most_vars: u64,
    ) -> PopDistSums {
        let options = VcfOptions {
            ploidy: 2,
            ..VcfOptions::default()
        };
        let mut reader =
            VcfReader::from_path(&reference(name), options).expect("the reader of the panel");
        reader.set_needs(Needs::GTS | Needs::CHROM_POS);
        let pops = pops_of_the_file(pops_file, reader.individuals(), &["p0", "p1", "p2"]);
        let per_var =
            PopDistPerVar::new(2, min_num_individuals).expect("the values of one variant");
        let mut walk = JackknifeWalk::new(how).expect("the walk over the groups");
        let num_pairs = pops.len() * (pops.len() - 1) / 2;
        let mut counts = vec![PopVarCounts::new(); pops.len()];
        let mut of_each_group: Vec<PairSums> = Vec::new();
        let mut num_vars = 0_u64;
        'over_the_blocks: while let Some(block) = reader.next_block().expect("a block of the panel")
        {
            for var in block.variants() {
                if num_vars == at_most_vars {
                    break 'over_the_blocks;
                }
                num_vars += 1;
                let (Some(chrom), Some(pos)) = (var.chrom(), var.pos()) else {
                    panic!("the reader gave a variant with no chromosome or no position");
                };
                let at = walk.group_of(chrom, pos).unwrap_or(0);
                of_each_group.resize((at + 1) * num_pairs, PairSums::default());
                for (pop, of_the_pop) in counts.iter_mut().enumerate() {
                    of_the_pop
                        .count_the_var(var.gts(), 2, pops.individuals(pop))
                        .expect("the counts of one population");
                }
                let mut pair = 0;
                for first in 0..pops.len() {
                    for second in (first + 1)..pops.len() {
                        if let Some(of_the_var) = per_var.of_var(&counts[first], &counts[second]) {
                            let into = &mut of_each_group[at * num_pairs + pair];
                            into.h_b += of_the_var.h_b;
                            into.h_w += of_the_var.h_w;
                            into.sqrt_of_the_products += of_the_var.sqrt_of_the_products;
                            into.corrected_h_s += of_the_var.corrected_h_s;
                            into.corrected_h_t += of_the_var.corrected_h_t;
                            into.num_vars += 1;
                        }
                        pair += 1;
                    }
                }
            }
        }
        PopDistSums::of_the_pass(pops.len(), num_vars, walk.groups().to_vec(), of_each_group)
    }

    /// The six sums of the one pair of the worked example over its four
    /// variants, cut into the groups `how` says, the variants lying 1000
    /// base pairs apart on one chromosome.
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "the four variants and the one pair of the worked example"
    )]
    fn sums_of_the_worked_example(how: JackknifeGroups) -> PopDistSums {
        let per_var = PopDistPerVar::new(2, 1).expect("the values of one variant");
        let mut walk = JackknifeWalk::new(how).expect("the walk over the groups");
        let mut of_each_group: Vec<PairSums> = Vec::new();
        for (var, gts) in WORKED_EXAMPLE.iter().enumerate() {
            let at = walk.group_of(0, 1000 * (var as u64 + 1)).unwrap_or(0);
            of_each_group.resize(at + 1, PairSums::default());
            let Some(of_the_var) =
                per_var.of_var(&counts_of(gts, 2, &POP1), &counts_of(gts, 2, &POP2))
            else {
                panic!("the variant {var} does not count for the pair");
            };
            let into = &mut of_each_group[at];
            into.h_b += of_the_var.h_b;
            into.h_w += of_the_var.h_w;
            into.sqrt_of_the_products += of_the_var.sqrt_of_the_products;
            into.corrected_h_s += of_the_var.corrected_h_s;
            into.corrected_h_t += of_the_var.corrected_h_t;
            into.num_vars += 1;
        }
        PopDistSums::of_the_pass(2, 4, walk.groups().to_vec(), of_each_group)
    }

    /// Two numbers that a reference program printed to the same digits, the
    /// tolerance being the one the item of `docs/specs/dists.md` gives.
    fn assert_it_is_within(found: Option<f64>, expected: f64, tolerance: f64, what: &str) {
        let Some(found) = found else {
            panic!("{what} has no value, and it is {expected}");
        };
        assert!(
            (found - expected).abs() <= tolerance,
            "{what} is {found} and not {expected}, {} away and not {tolerance}",
            (found - expected).abs()
        );
    }

    /// A number of a reference program that is the same estimator computed
    /// in another order, within 1e-12 of it relative.
    fn assert_it_is_the_same_number(found: Option<f64>, expected: f64, what: &str) {
        assert_it_is_within(found, expected, expected.abs() * 1e-12, what);
    }

    /// The F_ST and the f_2 of the worked example of "How it is verified" of
    /// `docs/specs/dists.md`: the sums over its four variants are 2.027778
    /// for H_b and 1.466667 for H_w, so F_ST is 0.561111 / 2.027778 and f_2
    /// is 0.561111 / 4. Both are ratios of the sums and not means of the
    /// per variant ratios, which for F_ST is undefined at the third
    /// variant.
    #[test]
    fn the_worked_example_has_the_fst_and_the_f2_of_the_spec() {
        let sums = sums_of_the_worked_example(JackknifeGroups::PerVariant);

        assert_it_is_within(
            sums.measure(PopDistMeasure::Fst, 0, 1),
            0.276712,
            1e-6,
            "the F_ST of the worked example",
        );
        assert_it_is_within(
            sums.measure(PopDistMeasure::F2, 0, 1),
            0.140278,
            1e-6,
            "the f_2 of the worked example",
        );
        assert_eq!(sums.num_vars_of(0, 1), Some(4));
        assert_eq!(sums.num_vars_of(1, 0), sums.num_vars_of(0, 1));
    }

    /// The groups the variants were cut into change no measure: the sums of
    /// a pair are added group by group and the division happens once,
    /// whether the four variants of the worked example fall in four groups,
    /// in one or in none.
    #[test]
    fn the_measures_are_the_same_whatever_the_groups() {
        let of_four = sums_of_the_worked_example(JackknifeGroups::PerVariant);
        let of_one = sums_of_the_worked_example(JackknifeGroups::OfBasePairs(100_000));
        let of_none = sums_of_the_worked_example(JackknifeGroups::None);

        assert_eq!(of_four.groups().len(), 4);
        assert_eq!(of_one.groups().len(), 1);
        assert!(of_none.groups().is_empty());
        for measure in [PopDistMeasure::Fst, PopDistMeasure::F2] {
            let found = of_four.measure(measure, 0, 1);
            assert_it_is_the_same_number(
                found,
                of_one.measure(measure, 0, 1).expect("the measure"),
                "the measure of one group",
            );
            assert_it_is_the_same_number(
                found,
                of_none.measure(measure, 0, 1).expect("the measure"),
                "the measure of no group",
            );
        }
    }

    /// The f_2 of one group, which f_3 and f_4 are built from later: with
    /// each variant its own group it is the f_2 of that variant, which the
    /// per variant table of the spec gives, the fourth of them negative.
    /// Beyond the groups there is none.
    #[test]
    fn the_f2_of_a_group_of_one_variant_is_that_variants_f2() {
        let sums = sums_of_the_worked_example(JackknifeGroups::PerVariant);

        for (var, f2) in [0.388889, 0.272222, 0.0, -0.1].into_iter().enumerate() {
            assert_it_is_within(
                sums.f2_of_group(var, 0, 1),
                f2,
                1e-6,
                &format!("the f_2 of the group of the variant {var}"),
            );
        }
        assert_eq!(sums.f2_of_group(4, 0, 1), None);
    }

    /// A pass of one variant has an F_ST equal to that variant's, so the
    /// per variant F_ST that plink2 reports for `var0000` of the biallelic
    /// panel, 0.332322 for p0 and p1, is the F_ST of a pass over that one
    /// variant.
    #[test]
    fn the_fst_of_a_pass_of_one_variant_is_the_fst_of_that_variant() {
        let sums = sums_of_the_panel(
            "dists/panel.vcf.gz",
            "stats/panel_pops.txt",
            JackknifeGroups::None,
            20,
            1,
        );

        assert_eq!(sums.num_vars(), 1);
        assert_eq!(sums.num_vars_of(0, 1), Some(1));
        assert_it_is_within(
            sums.measure(PopDistMeasure::Fst, 0, 1),
            0.332322,
            1e-6,
            "the F_ST of var0000 for p0 and p1",
        );
    }

    /// Hudson's F_ST of the three pairs of both panels against plink2
    /// v2.0.0-a.7.7, which computes this estimator and sums it the same
    /// way, from "How it is verified" of the F_ST item of
    /// `docs/specs/dists.md`. plink2 prints six digits, and the furthest of
    /// the six numbers is 4.9e-7 from what popnei gives, so the comparison
    /// is within 1e-6 absolute. The multiallelic panel is the one that says
    /// the arithmetic does not assume two alleles: plink2 counts every
    /// allele of a record there and so does popnei.
    #[test]
    fn the_fst_of_both_panels_is_plink2s() {
        let panels = [
            (
                "dists/panel.vcf.gz",
                "stats/panel_pops.txt",
                1200_u64,
                [0.104962, 0.102736, 0.109621],
            ),
            (
                "pop_dists/micro.vcf.gz",
                "pop_dists/micro_pops.txt",
                120,
                [0.0642281, 0.0694106, 0.0699954],
            ),
        ];
        for (panel, pops_file, num_vars, expected) in panels {
            let sums = sums_of_the_panel(panel, pops_file, JackknifeGroups::None, 20, u64::MAX);
            assert_eq!(sums.num_vars(), num_vars, "{panel}");
            for (pair, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
                assert_eq!(sums.num_vars_of(i, j), Some(num_vars), "{panel}");
                assert_it_is_within(
                    sums.measure(PopDistMeasure::Fst, i, j),
                    expected[pair],
                    1e-6,
                    &format!("the F_ST of the pair {i} {j} of {panel}"),
                );
            }
        }
    }

    /// f_2 and its jackknife standard error for the three pairs of the
    /// biallelic panel against ADMIXTOOLS 2.0.10, which computes f_2 with
    /// the same delete-m jackknife over the same groups, cut the same way.
    /// `expected` is the f_2 and the standard error of p0-p1, p0-p2 and
    /// p1-p2, which `tests/reference/pop_dists/panel.f2.tsv` and
    /// `panel.f2.uneven.tsv` hold to seventeen digits. The two agree to the
    /// last bits of a double, so the comparison is within 1e-12 relative.
    fn assert_the_f2_of_the_panel_is(how: JackknifeGroups, expected: [(f64, f64); 3]) {
        let sums = sums_of_the_panel(
            "dists/panel.vcf.gz",
            "stats/panel_pops.txt",
            how,
            20,
            u64::MAX,
        );

        for (pair, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
            let (f2, standard_error) = expected[pair];
            assert_eq!(sums.num_vars_of(i, j), Some(1200));
            assert_it_is_the_same_number(
                sums.measure(PopDistMeasure::F2, i, j),
                f2,
                &format!("the f_2 of the pair {i} {j}"),
            );
            assert_it_is_the_same_number(
                sums.standard_error(PopDistMeasure::F2, i, j),
                standard_error,
                &format!("the standard error of the f_2 of the pair {i} {j}"),
            );
        }
    }

    /// The 12 groups that a length of 100 000 base pairs cuts the panel
    /// into, which is the `blgsize = 100000` of the run of "How it is
    /// verified" of the f_2 item of `docs/specs/dists.md`: 100 variants
    /// each.
    #[test]
    #[expect(
        clippy::excessive_precision,
        reason = "the numbers as the files of tests/reference/pop_dists/ and the spec print them, which is one digit more than a f64 keeps"
    )]
    fn the_f2_of_the_panel_and_its_standard_error_are_admixtools() {
        assert_the_f2_of_the_panel_is(
            JackknifeGroups::OfBasePairs(100_000),
            [
                (0.041181109098151744, 0.0018034174115346543),
                (0.039890789655075122, 0.0015916239936530866),
                (0.042798563747209542, 0.0022726350509508155),
            ],
        );
    }

    /// The same against the same program over groups that are not all of
    /// one size, which is the second run of "The standard errors" of the
    /// spec: a length of 250 000 base pairs cuts the panel into 6 groups
    /// that hold 250, 250 and 100 variants on each of its two chromosomes.
    /// The 12 groups above all hold 100, and groups of one size cannot tell
    /// the delete-m jackknife for unequal m from an estimator that takes
    /// the groups as equal, since h_j is then one number for every group.
    /// The furthest of these three standard errors is 5.7e-14 of itself
    /// from what ADMIXTOOLS gives.
    #[test]
    #[expect(
        clippy::excessive_precision,
        reason = "the numbers as the files of tests/reference/pop_dists/ and the spec print them, which is one digit more than a f64 keeps"
    )]
    fn the_standard_error_over_groups_of_different_sizes_is_admixtools() {
        assert_the_f2_of_the_panel_is(
            JackknifeGroups::OfBasePairs(250_000),
            [
                (0.041181109098151737, 0.0012853754149193656),
                (0.039890789655075122, 0.00099457882673528123),
                (0.042798563747209514, 0.0030341885486622044),
            ],
        );
    }

    /// The groups of that second run, which the positions of the panel show
    /// to be uneven: its variants lie every 1000 base pairs from 1000 to
    /// 600 000 on each of two chromosomes, so the first two groups of a
    /// chromosome span 250 000 base pairs and hold 250 variants and the
    /// third spans 100 000 and holds 100.
    #[test]
    fn a_length_of_250_000_base_pairs_cuts_the_panel_into_groups_of_two_sizes() {
        let sums = sums_of_the_panel(
            "dists/panel.vcf.gz",
            "stats/panel_pops.txt",
            JackknifeGroups::OfBasePairs(250_000),
            20,
            u64::MAX,
        );

        let spans: Vec<(u32, u64, u64)> = sums
            .groups()
            .iter()
            .map(|group| (group.chrom, group.start, group.end))
            .collect();
        assert_eq!(
            spans,
            [
                (0, 1000, 250_000),
                (0, 251_000, 500_000),
                (0, 501_000, 600_000),
                (1, 1000, 250_000),
                (1, 251_000, 500_000),
                (1, 501_000, 600_000),
            ]
        );
    }

    /// The standard error of the F_ST of the same panel and the same
    /// groups, which no program prints: what is checked here is that it
    /// comes out of the same jackknife as the f_2 above and is a small part
    /// of the measure, 0.0032 of a F_ST of 0.10. Its arithmetic is the one
    /// ADMIXTOOLS checks.
    #[test]
    fn the_fst_has_a_standard_error_of_the_same_jackknife() {
        let sums = sums_of_the_panel(
            "dists/panel.vcf.gz",
            "stats/panel_pops.txt",
            JackknifeGroups::OfBasePairs(100_000),
            20,
            u64::MAX,
        );

        let standard_error = sums
            .standard_error(PopDistMeasure::Fst, 0, 1)
            .expect("the standard error of the F_ST of p0 and p1");
        assert!(
            standard_error > 0.0 && standard_error < 0.01,
            "the standard error of the F_ST of p0 and p1 is {standard_error}"
        );
    }

    /// A pair whose variants all fell in one group has no standard error,
    /// and neither has one of a pass that was asked for no groups: leaving
    /// the one group out leaves no variant to calculate the measure from.
    /// The measures themselves are there in both cases.
    #[test]
    fn a_pair_of_one_group_or_of_none_has_no_standard_error() {
        for how in [JackknifeGroups::OfBasePairs(100_000), JackknifeGroups::None] {
            let sums = sums_of_the_worked_example(how);

            assert!(sums.measure(PopDistMeasure::F2, 0, 1).is_some());
            assert_eq!(sums.standard_error(PopDistMeasure::F2, 0, 1), None);
        }
    }

    /// The pairs and the measures a caller asks for by number: two
    /// populations that are one, or one that is not a population, have no
    /// measure and no count of variants, and the five measures that the
    /// work packages 2 and 3 of `docs/plans/dists-pops.md` add have none
    /// yet. `measures` gives the pairs in the order of the distance vector,
    /// which for three populations is p0-p1, p0-p2 and p1-p2.
    #[test]
    fn a_pair_that_is_not_one_and_a_measure_not_written_yet_have_no_value() {
        let sums = sums_of_the_worked_example(JackknifeGroups::PerVariant);

        assert_eq!(sums.num_pops(), 2);
        assert_eq!(sums.measure(PopDistMeasure::Fst, 0, 0), None);
        assert_eq!(sums.measure(PopDistMeasure::Fst, 0, 2), None);
        assert_eq!(sums.num_vars_of(0, 0), None);
        assert_eq!(sums.standard_error(PopDistMeasure::F2, 0, 0), None);
        for measure in [
            PopDistMeasure::Chord,
            PopDistMeasure::Da,
            PopDistMeasure::Dest,
            PopDistMeasure::Gst,
            PopDistMeasure::GstStandardized,
        ] {
            assert_eq!(sums.measure(measure, 0, 1), None);
            assert_eq!(sums.standard_error(measure, 0, 1), None);
        }
        let of_three = sums_of_the_panel(
            "dists/panel.vcf.gz",
            "stats/panel_pops.txt",
            JackknifeGroups::None,
            20,
            1,
        );
        let in_order: Vec<Option<f64>> = of_three.measures(PopDistMeasure::Fst).collect();
        assert_eq!(in_order.len(), 3);
        for (pair, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
            assert_eq!(in_order[pair], of_three.measure(PopDistMeasure::Fst, i, j));
        }
    }
}
