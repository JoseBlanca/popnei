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
//! [`calc_pop_dist_sums`] is the pass itself, one reading of the blocks of a
//! reader with the rows of each block on the threads of rayon, and
//! [`PopDistOptions`] what it is asked for. [`PopDistSums`] is what it
//! leaves: the six sums of every pair within every group.
//! [`PopDistSums::measure`] turns them into one of the seven
//! measures, which [`PopDistMeasure`] names, and
//! [`PopDistSums::standard_error`] into the standard error of that measure,
//! by leaving each group out in turn and asking how far the measure moves,
//! which is the delete-m jackknife of Busing, Meijer and van der Leeden.
//!
//! `docs/specs/dists.md` has the design, the formulas and the numbers the
//! tests assert, and the row `dists` of section 9 of
//! `docs/architecture.md` is where the module sits.

use crate::block::{Block, BlockReader, ROWS_PER_CHUNK, alleles_of_a_chunk, alleles_per_var_of};
use crate::dists::{index_of_the_pair, num_pairs_of};
use crate::error::{Error, Result};
use crate::io::vcf::MAX_PLOIDY;
use crate::stats::{ObsHet, Pops, raised};
use crate::variant::{AlleleCounts, ChromTable, GtCounts, Needs, count_alleles_of, count_gts_of};

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
    ///
    /// The counts are read by the pair of populations through the fields,
    /// which is where they are counted, so this is what the tests of the
    /// counts table of `docs/specs/dists.md` read and nothing else.
    #[cfg(test)]
    pub(crate) fn allele_counts(&self) -> &AlleleCounts {
        &self.allele_counts
    }

    /// The called alleles of the population at the variant, the sum of
    /// [`PopVarCounts::allele_counts`], which the tests of that table read.
    #[cfg(test)]
    pub(crate) fn called_alleles(&self) -> u32 {
        self.called_alleles
    }

    /// The counts of the genotypes: the called, the missing and the
    /// heterozygous ones, which the tests of that table read.
    #[cfg(test)]
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
    ///
    /// The f_2 a user reads is a ratio of the sums of a pass, which
    /// [`value_of`] takes, so this is what the per variant table of
    /// `docs/specs/dists.md` is checked against and nothing else.
    #[cfg(test)]
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

/// The ploidy the reader of a pass says its genotypes hold, as the number
/// the allele frequencies of a population are raised to.
///
/// # Errors
///
/// A ploidy of 0 or above [`MAX_PLOIDY`]. The VCF reader refuses both when
/// it is opened and the vars file reader refuses a file whose genotypes
/// hold no allele, so what reaches this is a vars file that says its
/// genotypes hold more alleles than popnei reads.
fn ploidy_of_the_variants(ploidy: usize) -> Result<u32> {
    let out_of_range = || Error::PopDistsPloidyOutOfRange {
        ploidy,
        largest: MAX_PLOIDY,
    };
    if ploidy == 0 || ploidy > MAX_PLOIDY {
        return Err(out_of_range());
    }
    u32::try_from(ploidy).map_err(|_| out_of_range())
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
            ploidy: ploidy_of_the_variants(ploidy)?,
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
    /// The chromosome and the position of the variant before, which say
    /// whether the variant in hand goes back, and `None` before the first
    /// variant.
    the_one_before: Option<(u32, u64)>,
    /// The chromosomes the variants have been on, in the order they first
    /// came, which say whether a chromosome comes back after another. A
    /// genome has a handful of them, and the list is looked at only where
    /// the chromosome changes.
    chroms_read: Vec<u32>,
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
            the_one_before: None,
            chroms_read: Vec::new(),
        })
    }

    /// The group the variant at `chrom` and `pos` falls in, counted from 0,
    /// and `None` when no groups were asked for.
    ///
    /// The variants are given to it in the order the reader gives them. A
    /// group that is being filled takes the variant when it is of its
    /// chromosome and its position is less than the length beyond the first
    /// variant of the group; otherwise the variant starts a group of its
    /// own. `chroms` is the table of the reader, which names the chromosome
    /// of an error.
    ///
    /// # Errors
    ///
    /// Where the groups are stretches of a chromosome, a variant that goes
    /// back: one whose position is below the position of the variant before
    /// it on the same chromosome, and one of a chromosome that the variant
    /// before it had left.
    pub(crate) fn group_of(
        &mut self,
        chroms: &ChromTable,
        chrom: u32,
        pos: u64,
    ) -> Result<Option<usize>> {
        match self.how {
            JackknifeGroups::None => return Ok(None),
            JackknifeGroups::PerVariant => {}
            JackknifeGroups::OfBasePairs(length) => {
                self.check_the_order(chroms, chrom, pos)?;
                // The check above leaves the position of the variant at or
                // beyond the first position of the group being filled, so
                // the subtraction is the distance between the two and
                // saturates at nothing.
                if let Some(filling) = self.groups.last_mut()
                    && filling.chrom == chrom
                    && pos.saturating_sub(filling.start) < length
                {
                    filling.end = pos;
                    return Ok(self.groups.len().checked_sub(1));
                }
            }
        }
        let at = self.groups.len();
        self.groups.push(GroupId {
            chrom,
            start: pos,
            end: pos,
        });
        Ok(Some(at))
    }

    /// It refuses a variant that goes back and keeps it as the variant
    /// before for the next one.
    ///
    /// The cut of a group compares the position of a variant with the first
    /// position of the group being filled, so a variant whose position is
    /// below the position of the variant before it on the same chromosome
    /// joins that group instead of starting one of its own, and a variant
    /// of a chromosome that the variant before it had left is cut into
    /// groups over the stretch the earlier variants of that chromosome were
    /// already cut into. Either way the groups are not the stretches the
    /// caller asked for. "The standard errors" of `docs/specs/dists.md` has
    /// what the first of the two does to the standard error of the
    /// biallelic panel.
    ///
    /// Two variants at one position are taken, as they are by the linkage
    /// disequilibrium filter of `docs/specs/filters.md`, the other part of
    /// popnei that needs this order.
    ///
    /// # Errors
    ///
    /// The variant that goes back, and the variant of a chromosome that had
    /// been left.
    fn check_the_order(&mut self, chroms: &ChromTable, chrom: u32, pos: u64) -> Result<()> {
        match self.the_one_before {
            Some((before_chrom, before)) if before_chrom == chrom => {
                if pos < before {
                    return Err(Error::JackknifeGroupsVariantGoesBack {
                        chrom: named(chroms, chrom),
                        pos,
                        before,
                    });
                }
            }
            Some((before_chrom, before)) => {
                if self.chroms_read.contains(&chrom) {
                    return Err(Error::JackknifeGroupsChromComesBack {
                        chrom: named(chroms, chrom),
                        pos,
                        before_chrom: named(chroms, before_chrom),
                        before,
                    });
                }
                self.chroms_read.push(chrom);
            }
            None => self.chroms_read.push(chrom),
        }
        self.the_one_before = Some((chrom, pos));
        Ok(())
    }

    /// The groups the variants walked so far were cut into, in the order
    /// they were started.
    pub(crate) fn groups(&self) -> &[GroupId] {
        &self.groups
    }
}

/// The name the table of a reader gives the chromosome `chrom`, which an
/// error of the groups says the chromosome by.
///
/// A number the table has no name for is named by the number itself, which
/// only a reader with a defect gives: the numbers are the table's own, one
/// for each name it was given.
fn named(chroms: &ChromTable, chrom: u32) -> String {
    chroms
        .name(chrom)
        .map_or_else(|| chrom.to_string(), ToOwned::to_owned)
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
    /// Jost's D, how much of the allelic variety of the two is not
    /// shared. It is the D_est of Jost (2008) under the correction of Nei
    /// and Chesser (1983) for the individuals it was estimated from, which
    /// is the estimator pyNei computes and the one GenAlEx prints. mmod's
    /// `pairwise_D` in R computes another estimator of the same quantity,
    /// leaving the observed heterozygosity out of the correction and
    /// dividing by 2n - 1 where this one divides by n - 1: the two are
    /// 7.3e-5 apart at the furthest on the biallelic panel of the tests,
    /// 1200 variants of three populations of 48 to 84 individuals, and
    /// 3.5e-4 on the multiallelic one, 120 microsatellite loci of six
    /// alleles in three populations of 30.
    Dest,
    /// Nei's G_ST, the share of the diversity of the two that lies between
    /// them. It comes from the same two corrected means as
    /// [`PopDistMeasure::Dest`], so mmod's `pairwise_Gst_Nei` carries the
    /// same difference of estimator, 7.2e-5 at the furthest on the
    /// biallelic panel and 9.0e-5 on the multiallelic one.
    Gst,
    /// G_ST rescaled to reach 1, the standardized G''_ST of Meirmans and
    /// Hedrick (2011) and not Hedrick's earlier G'_ST. mmod's
    /// `pairwise_Gst_Hedrick` computes this one whatever its name
    /// suggests, and is 1.9e-4 from it at the furthest on the biallelic
    /// panel and 4.7e-4 on the multiallelic one.
    GstStandardized,
}

impl PopDistMeasure {
    /// The name of each of the seven measures, in the order of the variants
    /// above.
    ///
    /// The names are what a Python and a TypeScript user writes in
    /// `measures`, and each one is the field of the result that holds that
    /// measure. They are here and not in the binding crates so that a
    /// rename is one change and not four.
    pub const NAMES: [&'static str; 7] = [
        "fst",
        "f2",
        "chord",
        "da",
        "dest",
        "gst",
        "gst_standardized",
    ];

    /// The name a user writes for this measure, which is the field of the
    /// result that holds it.
    #[must_use]
    pub fn name(self) -> &'static str {
        let of_the_seven = match self {
            PopDistMeasure::Fst => 0,
            PopDistMeasure::F2 => 1,
            PopDistMeasure::Chord => 2,
            PopDistMeasure::Da => 3,
            PopDistMeasure::Dest => 4,
            PopDistMeasure::Gst => 5,
            PopDistMeasure::GstStandardized => 6,
        };
        // The seven names are there, one for each variant of the enum.
        PopDistMeasure::NAMES
            .get(of_the_seven)
            .copied()
            .unwrap_or("")
    }

    /// The measures a pass gives a value for, which is all seven of them
    /// since work package 3 of `docs/plans/dists-pops.md` added the chord
    /// distance and Nei's D_A to the five the work packages 1 and 2 wrote.
    ///
    /// Both packages refuse a measure that is not here, so that nobody
    /// reads a vector of NaN as a distance, and there is nothing left for
    /// them to refuse. It is in the core, as
    /// [`NAMES`](PopDistMeasure::NAMES) is, so that a measure is added to
    /// the two packages by adding the formula of [`value_of`] and this
    /// array, both of which are in this file. They are in the order of
    /// `NAMES`, which is the order a package names them in.
    pub const THAT_HAVE_A_VALUE: [PopDistMeasure; 7] = [
        PopDistMeasure::Fst,
        PopDistMeasure::F2,
        PopDistMeasure::Chord,
        PopDistMeasure::Da,
        PopDistMeasure::Dest,
        PopDistMeasure::Gst,
        PopDistMeasure::GstStandardized,
    ];

    /// Whether a pass gives this measure a value.
    #[must_use]
    pub fn has_a_value(self) -> bool {
        PopDistMeasure::THAT_HAVE_A_VALUE.contains(&self)
    }

    /// The three measures that are ratios of the corrected H_S and H_T,
    /// Jost's D, Nei's G_ST and the standardized G''_ST, which are the ones
    /// that have no value at a ploidy of 1.
    ///
    /// Both corrected values raise the allele frequencies of a variant to
    /// the ploidy and take the sum from 1, so at a ploidy of 1 they are 0
    /// by their definitions, as the observed heterozygosity of a haploid
    /// genotype is, and the sums of a pass hold nothing but the residue of
    /// adding frequencies that a `f64` does not bring to exactly 1.
    /// "Variants that do not count" of `docs/specs/dists.md` has what that
    /// residue looked like before the three were given no value there.
    #[must_use]
    pub fn is_of_the_corrected_diversities(self) -> bool {
        matches!(
            self,
            PopDistMeasure::Dest | PopDistMeasure::Gst | PopDistMeasure::GstStandardized
        )
    }

    /// The names of the measures that have a value, which the two packages
    /// name in the refusal of one that has none.
    #[must_use]
    pub fn names_that_have_a_value() -> Vec<&'static str> {
        PopDistMeasure::THAT_HAVE_A_VALUE
            .iter()
            .map(|measure| measure.name())
            .collect()
    }

    /// The measure a user named.
    ///
    /// # Errors
    ///
    /// A name that is of none of the seven, with the seven names.
    pub fn of_name(name: &str) -> Result<PopDistMeasure> {
        // The names are in [`PopDistMeasure::NAMES`] alone, in the order of
        // the variants, so a name that is renamed is renamed in one place.
        match PopDistMeasure::NAMES
            .iter()
            .position(|known| *known == name)
        {
            Some(0) => Ok(PopDistMeasure::Fst),
            Some(1) => Ok(PopDistMeasure::F2),
            Some(2) => Ok(PopDistMeasure::Chord),
            Some(3) => Ok(PopDistMeasure::Da),
            Some(4) => Ok(PopDistMeasure::Dest),
            Some(5) => Ok(PopDistMeasure::Gst),
            Some(6) => Ok(PopDistMeasure::GstStandardized),
            Some(_) | None => Err(Error::PopDistMeasureOfAnUnknownName {
                name: name.to_owned(),
            }),
        }
    }
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
/// The count is a `u64`, which is what popnei counts a whole dataset with,
/// and it costs nothing beside a `u32`: the five `f64` align the six
/// numbers to 8 bytes, so the six are 48 bytes with either count, which is
/// what "How it runs" of `docs/specs/dists.md` states.
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
    pub(crate) num_vars: u64,
}

impl PairSums {
    /// It adds what one variant gives the pair and counts that variant.
    fn add_the_var(&mut self, of_the_var: &VarSums) {
        self.h_b += of_the_var.h_b;
        self.h_w += of_the_var.h_w;
        self.sqrt_of_the_products += of_the_var.sqrt_of_the_products;
        self.corrected_h_s += of_the_var.corrected_h_s;
        self.corrected_h_t += of_the_var.corrected_h_t;
        // The variants that counted for a pair are at most the variants of
        // the pass, which is a count of rows: a pass of more than
        // 18446744073709551615 of them reads more rows than any source
        // holds.
        self.num_vars = self.num_vars.saturating_add(1);
    }

    /// It adds the sums of one chunk of rows, or of one group, into these.
    ///
    /// The threads of a block add up the rows of their own chunk and the
    /// chunks are added here in the order of the block, so the sums do not
    /// depend on how many threads read it.
    fn add_the_sums(&mut self, of_the_chunk: &PairSums) {
        self.h_b += of_the_chunk.h_b;
        self.h_w += of_the_chunk.h_w;
        self.sqrt_of_the_products += of_the_chunk.sqrt_of_the_products;
        self.corrected_h_s += of_the_chunk.corrected_h_s;
        self.corrected_h_t += of_the_chunk.corrected_h_t;
        self.num_vars = self.num_vars.saturating_add(of_the_chunk.num_vars);
    }

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
/// same allele at every variant that counted for them, has no F_ST: the
/// division is 0 over 0. It is a `None` and not the NaN that division
/// gives, although the binding crates write a NaN for a `None` too, because
/// a NaN of one group is a pseudo-value of NaN in
/// [`PopDistSums::standard_error`] and a standard error of NaN for the
/// whole pair, where the `None` is a group the jackknife can say something
/// about: "The standard errors" of the spec has what it says.
///
/// Jost's D, Nei's G_ST and the standardized G''_ST are ratios of the means
/// of the corrected H_S and H_T instead, which [`PopDistPerVar::of_var`]
/// gives one variant of. With H_S' and H_T' those two means over the n
/// variants that counted for the pair, and s = 2 populations, since a pair
/// is a pair whatever else the run holds:
///
/// ```text
/// D      = (s / (s - 1)) (H_T' - H_S') / (1 - H_S')
/// G_ST   = (H_T' - H_S') / H_T'
/// G''_ST = s (H_T' - H_S') / ((s H_T' - H_S') (1 - H_S'))
/// ```
///
/// Each of the three is a ratio of the means and not a mean of the per
/// variant ratios, which is what pyNei's `_calc_jost_from_ht_hs` computes
/// for D. Each has no value where its own divisor is 0, as F_ST has none
/// where the sum of H_b is: D where H_S' is exactly 1, G_ST where the mean
/// corrected H_T is 0, which is the pair fixed for the same allele above,
/// and G''_ST where either of those two happens, since it divides by both.
///
/// The chord distance and Nei's D_A come out of the third sum, the square
/// roots of the products of the frequencies added over the alleles of each
/// variant. With S that sum over the n variants that counted,
///
/// ```text
/// D_A   = 1 - S / n
/// chord = sqrt(D_A)
/// ```
///
/// which is the form `adegenet::dist.genpop(method = 2)` gives, the chord
/// of the sphere of radius 1 divided by the square root of 2. S is at most
/// n, so D_A is 0 at the least; where the rounding takes S above n, which
/// two populations with the same allele frequencies at every variant do,
/// both are 0 and not a D_A of -2.2e-16 and a chord that is NaN, as "What
/// it gives" of the chord item of the spec has it.
///
/// Jost's D, G_ST and G''_ST have no value at all at a ploidy of 1,
/// whatever the sums hold, which
/// [`is_of_the_corrected_diversities`](PopDistMeasure::is_of_the_corrected_diversities)
/// says why: H_S' and H_T' are 0 there by their definitions and the sums of
/// them are the residue of the rounding of the frequencies. F_ST, f_2, the
/// chord distance and Nei's D_A read neither of the two and are at a ploidy
/// of 1 what they are at any other.
///
/// The ploidy is the reader's, which the pass kept, and not one of the
/// sums.
fn value_of(measure: PopDistMeasure, sums: &PairSums, ploidy: u32) -> Option<f64> {
    if sums.num_vars == 0 {
        return None;
    }
    if ploidy == 1 && measure.is_of_the_corrected_diversities() {
        return None;
    }
    let between_minus_within = sums.h_b - sums.h_w;
    // Every count of popnei is below 2^53, where a `f64` holds the whole
    // numbers exactly.
    let num_vars = sums.num_vars as f64;
    let mean_h_s = sums.corrected_h_s / num_vars;
    let mean_h_t = sums.corrected_h_t / num_vars;
    let between_the_pops = mean_h_t - mean_h_s;
    let one_minus_mean_h_s = 1.0 - mean_h_s;
    // The sum of the square roots is at most the variants that counted, so
    // what the rounding leaves above them is 0 and not a chord distance of
    // NaN.
    let nei_d_a = (1.0 - sums.sqrt_of_the_products / num_vars).max(0.0);
    match measure {
        PopDistMeasure::Fst => (sums.h_b != 0.0).then(|| between_minus_within / sums.h_b),
        PopDistMeasure::F2 => Some(between_minus_within / num_vars),
        PopDistMeasure::Dest => (one_minus_mean_h_s != 0.0).then(|| {
            (NUM_POPS_OF_A_PAIR / (NUM_POPS_OF_A_PAIR - 1.0)) * between_the_pops
                / one_minus_mean_h_s
        }),
        PopDistMeasure::Gst => (mean_h_t != 0.0).then(|| between_the_pops / mean_h_t),
        PopDistMeasure::GstStandardized => {
            let divisor = (NUM_POPS_OF_A_PAIR * mean_h_t - mean_h_s) * one_minus_mean_h_s;
            (divisor != 0.0).then(|| NUM_POPS_OF_A_PAIR * between_the_pops / divisor)
        }
        PopDistMeasure::Chord => Some(nei_d_a.sqrt()),
        PopDistMeasure::Da => Some(nei_d_a),
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
    /// The ploidy of the reader the pass read, which the three measures
    /// built from the corrected H_S and H_T have no value at when it is 1.
    ploidy: u32,
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
    /// distance vector, and the ploidy of the reader it read.
    ///
    /// When `groups` is empty, which is what a caller who asked for no
    /// standard errors gets, `of_each_group` is one run of the pairs and
    /// holds every variant the pass counted.
    pub(crate) fn of_the_pass(
        num_pops: usize,
        ploidy: u32,
        num_vars: u64,
        groups: Vec<GroupId>,
        of_each_group: Vec<PairSums>,
    ) -> PopDistSums {
        PopDistSums {
            num_pops,
            ploidy,
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
        Some(self.total_of(self.index_of_the_pair(i, j)?)?.num_vars)
    }

    /// The measure for the pair. `None` where
    /// [`num_vars_of`](PopDistSums::num_vars_of) is 0 or `None`, and for
    /// the five measures the work packages 2 and 3 of
    /// `docs/plans/dists-pops.md` add.
    #[must_use]
    pub fn measure(&self, measure: PopDistMeasure, i: usize, j: usize) -> Option<f64> {
        self.of_the_pair(measure, self.index_of_the_pair(i, j)?)
    }

    /// How many pairs the populations make, which is how many values each
    /// of the iterators below gives for one measure and for one group, and
    /// 0 when the populations make more pairs than a `usize` counts.
    #[must_use]
    pub fn num_pairs(&self) -> usize {
        num_pairs_of(self.num_pops).unwrap_or(0)
    }

    /// The measure of every pair, in the order of the distance vector. The
    /// binding crates write NaN for a `None`.
    pub fn measures(&self, measure: PopDistMeasure) -> impl Iterator<Item = Option<f64>> + '_ {
        (0..self.num_pairs()).map(move |pair| self.of_the_pair(measure, pair))
    }

    /// The jackknife standard error of that measure for every pair, in the
    /// same order, `None` where
    /// [`standard_error`](PopDistSums::standard_error) gives one.
    pub fn standard_errors(
        &self,
        measure: PopDistMeasure,
    ) -> impl Iterator<Item = Option<f64>> + '_ {
        (0..self.num_pairs()).map(move |pair| self.standard_error_of(measure, pair))
    }

    /// The variants that counted for every pair, in the same order, and
    /// `None` for a pair the sums do not hold, which is a defect of popnei
    /// and not a pair that counted nothing.
    pub fn num_vars_of_each_pair(&self) -> impl Iterator<Item = Option<u64>> + '_ {
        (0..self.num_pairs()).map(move |pair| Some(self.total_of(pair)?.num_vars))
    }

    /// f_2 within every group, the pairs of one group together and the
    /// groups in the order they were started: the table of groups x pairs
    /// that the packages give as `f2_groups` and that f_3 and f_4 are built
    /// from later. It is empty when no groups were asked for.
    pub fn f2_of_every_group(&self) -> impl Iterator<Item = Option<f64>> + '_ {
        let num_pairs = self.num_pairs();
        (0..self.groups.len())
            .flat_map(move |group| (0..num_pairs).map(move |pair| self.f2_within(group, pair)))
    }

    /// f_2 within one group, which f_3 and f_4 are built from later.
    ///
    /// `None` when `group` is not one of
    /// [`groups`](PopDistSums::groups), when the two populations are one or
    /// either of them is not a population, and when no variant of the pair
    /// fell in the group. The last is of the group alone and not of the
    /// pair: a pair whose f_2 over every group is a number has no f_2 in a
    /// group it has no variant in.
    #[must_use]
    pub fn f2_of_group(&self, group: usize, i: usize, j: usize) -> Option<f64> {
        self.f2_within(group, self.index_of_the_pair(i, j)?)
    }

    /// The jackknife standard error of the measure of the pair. `None`
    /// where [`measure`](PopDistSums::measure) is `None`, where no groups
    /// were asked for, where every variant of the pair fell in one group,
    /// and where a group that holds variants of the pair leaves the measure
    /// without a value when it is taken out.
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
        self.standard_error_of(measure, self.index_of_the_pair(i, j)?)
    }

    /// f_2 of the pair at `pair` of the distance vector within the group at
    /// `group`, and `None` when the group is not one of the pass or no
    /// variant of the pair fell in it.
    fn f2_within(&self, group: usize, pair: usize) -> Option<f64> {
        if group >= self.groups.len() {
            return None;
        }
        value_of(
            PopDistMeasure::F2,
            self.of_the_group(group, pair)?,
            self.ploidy,
        )
    }

    /// The jackknife standard error of the measure of the pair at `pair` of
    /// the distance vector.
    fn standard_error_of(&self, measure: PopDistMeasure, pair: usize) -> Option<f64> {
        let over_all = self.total_of(pair)?;
        let over_all_value = value_of(measure, &over_all, self.ploidy)?;
        let mut num_groups: usize = 0;
        let mut jackknife_estimate = 0.0;
        for group in 0..self.groups.len() {
            let pseudo =
                match self.of_the_group_left_out(measure, pair, group, over_all_value, &over_all) {
                    OfTheGroupLeftOut::NoVariantOfThePair => continue,
                    OfTheGroupLeftOut::NoValueWithoutIt => return None,
                    OfTheGroupLeftOut::PseudoValue(pseudo) => pseudo,
                };
            #[expect(
                clippy::arithmetic_side_effects,
                reason = "the loop runs once for each of the groups, so the count reaches self.groups.len() at most, which is a usize"
            )]
            {
                num_groups += 1;
            }
            jackknife_estimate += pseudo.value / pseudo.weight;
        }
        if num_groups < 2 {
            return None;
        }
        let mut variance = 0.0;
        for group in 0..self.groups.len() {
            let pseudo =
                match self.of_the_group_left_out(measure, pair, group, over_all_value, &over_all) {
                    OfTheGroupLeftOut::NoVariantOfThePair => continue,
                    OfTheGroupLeftOut::NoValueWithoutIt => return None,
                    OfTheGroupLeftOut::PseudoValue(pseudo) => pseudo,
                };
            let from_the_estimate = pseudo.value - jackknife_estimate;
            variance += from_the_estimate * from_the_estimate / (pseudo.weight - 1.0);
        }
        // The groups are below 2^53, where a `f64` holds the whole numbers
        // exactly: each one holds a variant of the pair at least.
        Some((variance / num_groups as f64).sqrt())
    }

    /// What one group gives the jackknife of the measure: its pseudo-value,
    /// or that no variant of the pair fell in it, or that the measure has no
    /// value with it left out, which takes the standard error of the pair
    /// away.
    fn of_the_group_left_out(
        &self,
        measure: PopDistMeasure,
        pair: usize,
        group: usize,
        over_all_value: f64,
        over_all: &PairSums,
    ) -> OfTheGroupLeftOut {
        let Some(of_the_group) = self.of_the_group(group, pair) else {
            return OfTheGroupLeftOut::NoVariantOfThePair;
        };
        if of_the_group.num_vars == 0 {
            return OfTheGroupLeftOut::NoVariantOfThePair;
        }
        // The counts are below 2^53, where a `f64` holds the whole numbers
        // exactly, and the group holds a variant, which the line above says.
        let weight = over_all.num_vars as f64 / of_the_group.num_vars as f64;
        let Some(without_the_group) = over_all
            .without(of_the_group)
            .and_then(|rest| value_of(measure, &rest, self.ploidy))
        else {
            return OfTheGroupLeftOut::NoValueWithoutIt;
        };
        OfTheGroupLeftOut::PseudoValue(PseudoValue {
            value: weight * over_all_value - (weight - 1.0) * without_the_group,
            weight,
        })
    }

    /// The measure of the pair at `pair` of the distance vector.
    fn of_the_pair(&self, measure: PopDistMeasure, pair: usize) -> Option<f64> {
        value_of(measure, &self.total_of(pair)?, self.ploidy)
    }

    /// The six sums of the pair over every variant that counted for it, its
    /// groups added in the order they were started so that the total does
    /// not depend on how the pass was cut into blocks or threads.
    ///
    /// `None` when the populations make no pair and when the counts of the
    /// groups add to more than a `u64` holds, which no pass of a source
    /// reaches: a variant counts once for a pair, so the total is at most
    /// the rows the pass read.
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

/// What one group left out gives the jackknife of one measure of one pair.
///
/// A group that holds no variant of the pair is not one of the g of "The
/// standard errors" of `docs/specs/dists.md`, and the jackknife goes on
/// without it. A group that holds variants of the pair and leaves the
/// measure without a value when it is taken out ends the standard error of
/// that measure: the weights 1/h_j of the g groups add to 1, so leaving one
/// of them out of the sums would drop the jackknife estimate below the
/// measure by that group's weight and the variance would be taken around a
/// centre no group put there.
#[derive(Debug, Clone, Copy)]
enum OfTheGroupLeftOut {
    /// The pseudo-value u_j of the group and its weight h_j.
    PseudoValue(PseudoValue),
    /// No variant of the pair fell in the group.
    NoVariantOfThePair,
    /// The group holds variants of the pair and the measure has no value
    /// with it left out.
    NoValueWithoutIt,
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

/// How many resampling groups the variants of a pass have to fall into for
/// a standard error to be given.
///
/// Each group is left out in turn and the measure calculated again from the
/// sums of the others, so a standard error built from a handful of groups
/// says more about where the cuts fell than about the populations. It is
/// the 20 that "The standard errors" of `docs/specs/dists.md` asks for, and
/// nobody has measured whether 20 is the right number.
pub const MIN_NUM_JACKKNIFE_GROUPS: usize = 20;

/// What a pass over the variants is asked for: how much data a population
/// needs at a variant for that variant to count, and how the variants are
/// cut into the resampling groups.
///
/// It has no `Default`: both fields are the caller's to decide, and the
/// groups have no default anywhere in popnei, since the right length
/// depends on the linkage disequilibrium of the populations being compared,
/// which popnei cannot know.
#[derive(Debug)]
pub struct PopDistOptions {
    /// How many genotypes a population has to have called at a variant for
    /// that variant to count for a pair the population is in.
    pub min_num_individuals: u32,
    /// How the variants are cut into the groups the standard errors are
    /// resampled over.
    pub groups: JackknifeGroups,
}

/// The six sums of every pair of `pops` over the variants `reader` gives,
/// which is one pass over the source through the steps the variants carry.
///
/// `reader` is the outermost reader of the chain of the pass, lent and not
/// taken, so that whoever built the chain reads the counts of its filters
/// from it when this returns: those counts and [`PopDistSums::num_vars`]
/// are the `pass_stats` of a result in Python and in TypeScript. The pass
/// asks the reader for the genotypes alone when it was asked for no groups,
/// and for the genotypes with the chromosome and the position otherwise,
/// which is what cutting the groups and naming them takes.
///
/// Every measure of every pair comes out of the sums it gives, and a
/// variant counts for a pair only where both of its populations have at
/// least `min_num_individuals` called genotypes, so each pair has its own
/// count of variants.
///
/// # Errors
///
/// Fewer than two populations, which make no pair; a pass that gave no
/// variant, whether its source holds none or its steps kept none of them;
/// fewer than [`MIN_NUM_JACKKNIFE_GROUPS`] groups where groups were asked
/// for, with how many the variants fell into; resampling groups of 0 base
/// pairs; a variant that goes back where the groups are stretches of a
/// chromosome, which is a position below the position of the variant
/// before it on that chromosome or a chromosome that the variant before it
/// had left; a ploidy of 0 or above the largest one a reader of popnei
/// gives;
/// the memory of the sums, which grows as the groups appear; what the
/// reader fails with; a block that holds no genotypes, or no positions
/// where the groups are cut from them, a block of no variants and a block
/// whose individuals or ploidy are not the ones the reader says its source
/// has, each a defect of a reader; and what the counts of one variant
/// refuse.
pub fn calc_pop_dist_sums<R: BlockReader + ?Sized>(
    reader: &mut R,
    pops: &Pops,
    options: &PopDistOptions,
) -> Result<PopDistSums> {
    if pops.len() < 2 {
        return Err(Error::PopDistsOfFewerThanTwoPops {
            num_pops: pops.len(),
        });
    }
    let sums = sums_of_the_pass(reader, pops, options)?;
    if sums.num_vars() == 0 {
        let filters = reader.filtering_stats();
        return Err(Error::PassGaveNoVariant {
            // The filter nearest the source was given what the source
            // gave; with no filter the pass gave what the source gave,
            // which is nothing.
            num_vars_of_the_source: filters.last().map_or(0, |(_, stats)| stats.vars_processed),
            filters,
        });
    }
    let num_groups = sums.groups().len();
    if options.groups != JackknifeGroups::None && num_groups < MIN_NUM_JACKKNIFE_GROUPS {
        return Err(Error::TooFewJackknifeGroups {
            num_groups,
            at_least: MIN_NUM_JACKKNIFE_GROUPS,
        });
    }
    Ok(sums)
}

/// The same pass without the checks [`calc_pop_dist_sums`] makes of what it
/// was asked for, which is what the tests of the standard errors use: the
/// two runs of ADMIXTOOLS 2 they are compared with cut the biallelic panel
/// into 12 and into 6 groups, and a user is asked for
/// [`MIN_NUM_JACKKNIFE_GROUPS`].
///
/// # Errors
///
/// Those of [`calc_pop_dist_sums`] but the three it makes itself: the
/// populations, the variants of the pass and how many groups they fell
/// into.
pub(crate) fn sums_of_the_pass<R: BlockReader + ?Sized>(
    reader: &mut R,
    pops: &Pops,
    options: &PopDistOptions,
) -> Result<PopDistSums> {
    let needs = match options.groups {
        // A group is cut from the chromosome and the position of each
        // variant, and it carries them, so the two are read for every way
        // of cutting the groups but no groups at all.
        JackknifeGroups::None => Needs::GTS,
        JackknifeGroups::PerVariant | JackknifeGroups::OfBasePairs(_) => {
            Needs::GTS | Needs::CHROM_POS
        }
    };
    reader.set_needs(needs);
    let mut walk = JackknifeWalk::new(options.groups)?;
    let num_pops = pops.len();
    let of_the_pass = OfThePass {
        pops,
        per_var: PopDistPerVar::new(reader.ploidy(), options.min_num_individuals)?,
        num_pairs: num_pairs_of(num_pops).ok_or(Error::PopDistsOfTooManyPops { num_pops })?,
    };
    let num_individuals = reader.individuals().len();
    let ploidy = reader.ploidy();
    let mut of_each_group: Vec<PairSums> = Vec::new();
    // The group of each row of the block being read, kept from one block to
    // the next so that the pass allocates it once.
    let mut of_the_rows: Vec<usize> = Vec::new();
    let mut num_vars: u64 = 0;
    while let Some(block) = reader.next_block()? {
        let alleles_per_var = alleles_per_var_of(&block, num_individuals, ploidy)?;
        cut_the_rows_into_groups(&block, needs, reader.chroms(), &mut walk, &mut of_the_rows)?;
        // The sums hold every group the variants have fallen into so far,
        // and one run of the pairs when no groups were asked for and every
        // row falls in the same place.
        grow_the_sums(&mut of_each_group, walk.groups().len().max(1), &of_the_pass)?;
        add_the_block(
            &block,
            alleles_per_var,
            &of_the_rows,
            &of_the_pass,
            &mut of_each_group,
        )?;
        // Every variant of the block counts here, counted for a pair or
        // not: `num_vars` is what a user reads as the variants of the pass.
        // A `usize` is 64 bits natively and 32 in wasm, so the conversion
        // holds, and a pass of more than 18446744073709551615 variants
        // reads more rows than any source holds.
        num_vars = num_vars.saturating_add(u64::try_from(block.num_vars).unwrap_or(u64::MAX));
    }
    Ok(PopDistSums::of_the_pass(
        num_pops,
        of_the_pass.per_var.ploidy,
        num_vars,
        walk.groups().to_vec(),
        of_each_group,
    ))
}

/// What every row of a pass is read with: the populations, the values one
/// variant gives a pair and how many pairs the populations make.
#[derive(Debug)]
struct OfThePass<'a> {
    /// The populations, whose pairs are in the order of the distance
    /// vector.
    pops: &'a Pops,
    /// What one variant adds to the sums of one pair.
    per_var: PopDistPerVar,
    /// How many pairs the populations make, which is the length of one run
    /// of the sums.
    num_pairs: usize,
}

/// The group each variant of the block falls in, at the place of that
/// variant in `of_the_rows`, and 0 for every one of them when no groups
/// were asked for, since the sums are then one run of the pairs.
///
/// The groups are cut in the order the reader gives the variants, so they go
/// up from the one the block before ended in: a row is of the group of the
/// row before it, or of the one the walk started for it, which is the next.
///
/// # Errors
///
/// A block with no chromosome or no position where the groups are cut from
/// them, which is a reader that was asked for the two and gave neither.
fn cut_the_rows_into_groups(
    block: &Block,
    needs: Needs,
    chroms: &ChromTable,
    walk: &mut JackknifeWalk,
    of_the_rows: &mut Vec<usize>,
) -> Result<()> {
    of_the_rows.clear();
    if !needs.contains(Needs::CHROM_POS) {
        of_the_rows.resize(block.num_vars, 0);
        return Ok(());
    }
    let missing = Needs::CHROM_POS.difference(block.fields());
    if !missing.is_empty() {
        return Err(Error::FieldsNotInTheBlock { fields: missing });
    }
    for var in block.variants() {
        let (Some(chrom), Some(pos)) = (var.chrom(), var.pos()) else {
            return Err(Error::FieldsNotInTheBlock {
                fields: Needs::CHROM_POS,
            });
        };
        // The walk gives a group for every variant it is asked about but
        // when it was asked for no groups, which the line above left.
        of_the_rows.push(walk.group_of(chroms, chrom, pos)?.unwrap_or(0));
    }
    Ok(())
}

/// It makes room in the sums of the pass for the six numbers of every pair
/// within each of `num_groups` groups, over what they hold already.
///
/// The groups appear while the variants are read, so the sums grow block by
/// block instead of being asked for once. The memory they reach is 48 bytes
/// for each pair and each group, which "How it runs" of
/// `docs/specs/dists.md` counts: 72 KB for 3 populations and 500 groups and
/// 29 MB for 50 populations, which is 1225 pairs.
///
/// # Errors
///
/// When the pairs and the groups are more than a `usize` counts, and when
/// the machine does not give their memory. It is asked for with
/// `try_reserve`, which gives it back as an error where a `resize` alone
/// would end the process.
fn grow_the_sums(
    of_each_group: &mut Vec<PairSums>,
    num_groups: usize,
    of_the_pass: &OfThePass,
) -> Result<()> {
    let too_large = || Error::PopDistSumsTooLarge {
        num_pops: of_the_pass.pops.len(),
        num_groups,
    };
    let wanted = num_groups
        .checked_mul(of_the_pass.num_pairs)
        .ok_or_else(too_large)?;
    let Some(more) = wanted.checked_sub(of_each_group.len()) else {
        return Ok(());
    };
    if more == 0 {
        return Ok(());
    }
    of_each_group.try_reserve(more).map_err(|_| too_large())?;
    of_each_group.resize(wanted, PairSums::default());
    Ok(())
}

/// It adds what every row of a block gives its pairs into the sums of the
/// pass, each row into the group it fell in.
///
/// Natively the chunks of rows are read on the threads of rayon, as section
/// 3 of `docs/architecture.md` asks: no row reads another, each chunk adds
/// up what its own rows give, and the chunks are added into the sums of the
/// pass in the order of the block, so neither a count nor a sum depends on
/// how many threads read the block, which rayon's own `reduce` would make
/// them. The threads are those of the pool the caller is running in, and
/// rayon's global pool only when the caller is in none.
///
/// # Errors
///
/// What the counts of one variant refuse: genotypes that are not a whole
/// number of genotypes of the ploidy, a variant of more alleles than a
/// count of them holds, an allele below the missing one, and an individual
/// of a population beyond the row. The error is the one of the first row
/// that has one, wherever the threads found it: which of two bad rows a
/// thread reaches first depends on how the chunks were shared out, and a
/// user who reports a damaged file has to get the same message every time,
/// so the rows are read again, one after another, to find the first.
#[cfg(not(target_family = "wasm"))]
fn add_the_block(
    block: &Block,
    alleles_per_var: usize,
    of_the_rows: &[usize],
    of_the_pass: &OfThePass,
    of_each_group: &mut [PairSums],
) -> Result<()> {
    use rayon::iter::{IndexedParallelIterator, ParallelIterator};
    use rayon::slice::ParallelSlice;

    let of_the_chunks: Result<Vec<ChunkSums>> = block
        .gts
        .par_chunks(alleles_of_a_chunk(alleles_per_var))
        .zip(of_the_rows.par_chunks(ROWS_PER_CHUNK))
        .map(|(gts, of_the_rows)| {
            sums_of_the_chunk(gts, of_the_rows, alleles_per_var, block.ploidy, of_the_pass)
        })
        .collect();
    match of_the_chunks {
        Ok(of_the_chunks) => {
            for of_the_chunk in &of_the_chunks {
                of_the_chunk.add_into(of_each_group, of_the_pass.num_pairs);
            }
            Ok(())
        }
        // The second pass costs a read of the block, and it is made only
        // where the block is refused and nothing of it is given.
        Err(of_a_thread) => {
            let mut read_again = vec![PairSums::default(); of_each_group.len()];
            match add_the_chunks_one_by_one(
                block,
                alleles_per_var,
                of_the_rows,
                of_the_pass,
                &mut read_again,
            ) {
                Err(of_the_first_row) => Err(of_the_first_row),
                // The rows are the same rows, so the second pass finds an
                // error too; the error of the threads is what is left if it
                // ever did not.
                Ok(()) => Err(of_a_thread),
            }
        }
    }
}

/// The same numbers, with the chunks read one after another, which is what
/// wasm does: it has no threads.
///
/// # Errors
///
/// Those of the native [`add_the_block`], at the first row that has one.
#[cfg(target_family = "wasm")]
fn add_the_block(
    block: &Block,
    alleles_per_var: usize,
    of_the_rows: &[usize],
    of_the_pass: &OfThePass,
    of_each_group: &mut [PairSums],
) -> Result<()> {
    add_the_chunks_one_by_one(
        block,
        alleles_per_var,
        of_the_rows,
        of_the_pass,
        of_each_group,
    )
}

/// The chunks of the block read one after another, each added into the sums
/// of the pass before the next is read: what wasm runs, and what the
/// threads fall back on to find the first row that is an error.
///
/// The chunks are the same chunks the threads read, and they are added in
/// the same order, so wasm and a native build add the values of a block up
/// in the same order and the addition of two floats is exact on every
/// machine popnei runs on.
///
/// # Errors
///
/// Those of [`add_the_block`], at the first row that has one.
fn add_the_chunks_one_by_one(
    block: &Block,
    alleles_per_var: usize,
    of_the_rows: &[usize],
    of_the_pass: &OfThePass,
    of_each_group: &mut [PairSums],
) -> Result<()> {
    for (gts, of_the_rows) in block
        .gts
        .chunks(alleles_of_a_chunk(alleles_per_var))
        .zip(of_the_rows.chunks(ROWS_PER_CHUNK))
    {
        let of_the_chunk =
            sums_of_the_chunk(gts, of_the_rows, alleles_per_var, block.ploidy, of_the_pass)?;
        of_the_chunk.add_into(of_each_group, of_the_pass.num_pairs);
    }
    Ok(())
}

/// What one chunk of rows gives the pairs, group by group: the six sums of
/// every pair within each of the groups its rows fell in, the first of
/// those groups first.
///
/// The rows of a chunk are consecutive variants of the block and the groups
/// are cut in the order the reader gives them, so the groups of a chunk are
/// a run of the groups of the pass and one run of the pairs for each of
/// them is what the chunk holds.
#[derive(Debug)]
struct ChunkSums {
    /// The group of the first row of the chunk, which the sums below start
    /// at.
    first_group: usize,
    /// One run of the pairs for each group from `first_group` on, in the
    /// order of the distance vector within a group.
    of_each_group: Vec<PairSums>,
    /// How many pairs one of those runs holds.
    num_pairs: usize,
}

impl ChunkSums {
    /// It keeps the sums of the group that was being filled and empties
    /// them for the next group of the chunk.
    ///
    /// # Errors
    ///
    /// When the machine does not give the memory of one more group of the
    /// chunk, which is one run of the pairs.
    fn keep_the_group(
        &mut self,
        of_the_group: &mut [PairSums],
        of_the_pass: &OfThePass,
    ) -> Result<()> {
        self.of_each_group
            .try_reserve(of_the_group.len())
            .map_err(|_| Error::PopDistSumsTooLarge {
                num_pops: of_the_pass.pops.len(),
                // The groups of the chunk, which are the ones of the pass
                // its rows fell in.
                num_groups: self.num_groups().saturating_add(1),
            })?;
        self.of_each_group.extend_from_slice(of_the_group);
        of_the_group.fill(PairSums::default());
        Ok(())
    }

    /// How many groups the rows of the chunk fell in, which is how many
    /// runs of the pairs it holds.
    fn num_groups(&self) -> usize {
        self.of_each_group
            .len()
            .checked_div(self.num_pairs)
            .unwrap_or(0)
    }

    /// It adds the sums of the chunk into those of the pass, each group of
    /// the chunk into the same group of the pass.
    ///
    /// The sums of the pass hold every group the variants have fallen into,
    /// which the block made room for before its rows were read, so the run
    /// of this chunk is inside them.
    fn add_into(&self, of_each_group: &mut [PairSums], num_pairs: usize) {
        // The first group of the chunk is one of the groups the sums of the
        // pass were grown to hold, so the place of its first pair is inside
        // them and the product is below what a `usize` counts.
        let at = self.first_group.saturating_mul(num_pairs);
        for (of_the_pass, of_the_chunk) in
            of_each_group.iter_mut().skip(at).zip(&self.of_each_group)
        {
            of_the_pass.add_the_sums(of_the_chunk);
        }
    }
}

/// The six sums of every pair within each group the rows of one chunk fell
/// in.
///
/// `gts` holds whole rows of `alleles_per_var` alleles each and
/// `of_the_rows` the group of each of those rows, in the same order. The
/// counts of one row over one population are taken once and every pair the
/// population is in reads them.
///
/// # Errors
///
/// What the counts of one variant refuse, at the first row of the chunk
/// that has one, and the memory of the sums of the chunk.
fn sums_of_the_chunk(
    gts: &[i8],
    of_the_rows: &[usize],
    alleles_per_var: usize,
    ploidy: usize,
    of_the_pass: &OfThePass,
) -> Result<ChunkSums> {
    let num_pairs = of_the_pass.num_pairs;
    // One set of counts for each population and one run of the pairs, both
    // written over at each row, so no row allocates. What the chunk
    // allocates is one more run of the pairs at each group its rows fall
    // in.
    let mut counts = vec![PopVarCounts::new(); of_the_pass.pops.len()];
    let mut of_the_group = vec![PairSums::default(); num_pairs];
    let mut sums = ChunkSums {
        first_group: of_the_rows.first().copied().unwrap_or(0),
        of_each_group: Vec::new(),
        num_pairs,
    };
    let mut group_being_filled: Option<usize> = None;
    for (row, group) in gts.chunks_exact(alleles_per_var).zip(of_the_rows) {
        if group_being_filled.is_some_and(|filling| filling != *group) {
            sums.keep_the_group(&mut of_the_group, of_the_pass)?;
        }
        group_being_filled = Some(*group);
        for (pop, of_the_pop) in counts.iter_mut().enumerate() {
            of_the_pop.count_the_var(row, ploidy, of_the_pass.pops.individuals(pop))?;
        }
        add_the_pairs(&counts, &of_the_pass.per_var, &mut of_the_group)?;
    }
    if group_being_filled.is_some() {
        sums.keep_the_group(&mut of_the_group, of_the_pass)?;
    }
    Ok(sums)
}

/// It adds what one variant gives each pair of the populations `counts` are
/// of into `of_the_group`, the sums of the pairs within the group that
/// variant fell in.
///
/// The pairs are walked in the order of the distance vector, (0, 1), (0,
/// 2), ..., (1, 2), ..., which is the order `of_the_group` holds them in
/// and the order every result of the module gives them in. It is as long as
/// the pairs of `counts`, so every pair has its place there.
fn add_the_pairs(
    counts: &[PopVarCounts],
    per_var: &PopDistPerVar,
    of_the_group: &mut [PairSums],
) -> Result<()> {
    let mut of_the_pairs = of_the_group.iter_mut();
    for (first, of_one) in counts.iter().enumerate() {
        // The populations are at most as many as the individuals of the
        // source, which a `usize` counts.
        for of_the_other in counts.iter().skip(first.saturating_add(1)) {
            let Some(of_the_pair) = of_the_pairs.next() else {
                return Err(Error::PopDistSumsOfAnotherSize {
                    num_pops: counts.len(),
                    num_pairs: of_the_group.len(),
                });
            };
            if let Some(of_the_var) = per_var.of_var(of_one, of_the_other) {
                of_the_pair.add_the_var(&of_the_var);
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::path::{Path, PathBuf};

    use super::{
        GroupId, JackknifeGroups, JackknifeWalk, MIN_NUM_JACKKNIFE_GROUPS, PairSums,
        PopDistMeasure, PopDistOptions, PopDistPerVar, PopDistSums, PopVarCounts, VarSums,
        add_the_pairs, calc_pop_dist_sums, sums_of_the_pass,
    };
    use crate::block::{Block, BlockReader};
    use crate::error::{Error, Result};
    use crate::filters::FilteringStats;
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::stats::Pops;
    use crate::variant::{ChromTable, Needs};

    /// The genotypes of the five variants of the worked example of "How it
    /// is verified" of `docs/specs/dists.md`: 6 diploid individuals, the
    /// two alleles of each after those of the one before, with -1 for an
    /// allele that was not called.
    ///
    /// The variants are, with pop1 the first three individuals and pop2 the
    /// last three: one biallelic variant where the populations differ; one
    /// of three alleles; one where both are fixed for the allele 0, with a
    /// missing genotype in pop1 and a half called one in pop2; one that is
    /// the same in both, every genotype heterozygous; and one where the
    /// populations hold different alleles and differ in their called
    /// genotypes, 3 against 2, and in their called alleles, 6 against 4,
    /// which is what tells the corrections of H_S and H_T from the readings
    /// of them that the other four cannot.
    const WORKED_EXAMPLE: [[i8; 12]; 5] = [
        [0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1],
        [0, 1, 1, 2, 2, 2, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, -1],
        [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
        [0, 0, 0, 0, 0, 1, 1, 1, -1, -1, 1, 2],
    ];

    /// The genotypes of the variant at a ploidy of 4 of the same section,
    /// the four alleles of each of the same 6 individuals: pop1 has 3
    /// called genotypes of 12 called alleles and pop2 has 2 of 8, and the
    /// two hold different alleles.
    const AT_A_PLOIDY_OF_FOUR: [i8; 24] = [
        0, 0, 1, 1, 0, 1, 1, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, -1, -1,
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
            ([5, 1, 0], 6, 3, [0, 3, 1], 4, 2),
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
            (0.875, 0.416667, 0.458333, 0.353553, 0.410714, 0.642857),
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

    /// The variant at a ploidy of 4 of the same section, which is the one
    /// number of this file that the ploidy is read for: E_P and H_T raise
    /// the frequencies to it, and a suite whose genotypes all hold two
    /// alleles cannot tell the ploidy from the 2 of a square. Raising them
    /// to 2 here gives H_S' 0.407738 and H_T' 0.459077 against the
    /// 0.864330 and 0.873157 of the ploidy, and pyNei at a ploidy of 4
    /// gives the second pair.
    ///
    /// The two populations also differ in their called genotypes, 3 against
    /// 2, and in the alleles they hold.
    #[test]
    fn a_variant_at_a_ploidy_of_four_has_the_values_of_the_spec() {
        let gts = &AT_A_PLOIDY_OF_FOUR;
        let per_var = match PopDistPerVar::new(4, 1) {
            Ok(per_var) => per_var,
            Err(error) => panic!("{error}"),
        };
        let of_pop1 = counts_of(gts, 4, &POP1);
        let of_pop2 = counts_of(gts, 4, &POP2);

        assert_eq!(of_pop1.allele_counts()[..3], [6, 5, 1]);
        assert_eq!(of_pop1.called_alleles(), 12);
        assert_eq!(of_pop1.gts().called, 3);
        assert_eq!(of_pop2.allele_counts()[..3], [7, 1, 0]);
        assert_eq!(of_pop2.called_alleles(), 8);
        assert_eq!(of_pop2.gts().called, 2);
        let Some(sums) = per_var.of_var(&of_pop1, &of_pop2) else {
            panic!("the variant does not count for the pair");
        };
        assert_the_value_is(sums.h_b, 0.510417, "H_b");
        assert_the_value_is(sums.h_w, 0.435606, "H_w");
        assert_the_value_is(sums.f2(), 0.074811, "f_2");
        assert_the_value_is(
            sums.sqrt_of_the_products,
            0.889656,
            "the sum of the square roots",
        );
        assert_the_value_is(sums.corrected_h_s, 0.864330, "H_S'");
        assert_the_value_is(sums.corrected_h_t, 0.873157, "H_T'");
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
    ///
    /// The chromosomes are numbered as the table of a reader numbers them,
    /// 0 for `chr1` and 1 for `chr2`, which is what the walk names in an
    /// error.
    fn walk_over(
        how: JackknifeGroups,
        vars: &[(u32, u64)],
    ) -> Result<(Vec<Option<usize>>, Vec<GroupId>)> {
        let chroms = chroms_of_the_tests();
        let mut walk = match JackknifeWalk::new(how) {
            Ok(walk) => walk,
            Err(error) => panic!("{error}"),
        };
        let mut of_each_var = Vec::new();
        for (chrom, pos) in vars {
            of_each_var.push(walk.group_of(&chroms, *chrom, *pos)?);
        }
        Ok((of_each_var, walk.groups().to_vec()))
    }

    /// The table of chromosome names the tests number their variants by,
    /// `chr1` to `chr4`, which is what a reader of a VCF of four
    /// chromosomes would hold.
    fn chroms_of_the_tests() -> ChromTable {
        let mut chroms = ChromTable::new();
        for number in 1..=4 {
            chroms.intern(&format!("chr{number}"));
        }
        chroms
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
                let at = match walk.group_of(reader.chroms(), chrom, pos) {
                    Ok(Some(at)) => at,
                    Ok(None) => panic!("the variant at {chrom} {pos} fell in no group"),
                    Err(error) => panic!("{error}"),
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
        let (of_each_var, groups) = walk_over(JackknifeGroups::OfBasePairs(100_000), &vars)
            .expect("the walk over the variants");

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
        let (of_each_var, groups) = walk_over(JackknifeGroups::OfBasePairs(100_000), &vars)
            .expect("the walk over the variants");

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
        let (of_each_var, groups) =
            walk_over(JackknifeGroups::OfBasePairs(1), &vars).expect("the walk over the variants");

        assert_eq!(of_each_var, [Some(0), Some(0), Some(1)]);
        assert_eq!(groups.len(), 2);
    }

    /// Each variant its own group, which is what a panel of scattered
    /// microsatellite loci asks for: a group of one variant holds one
    /// position, its start and its end.
    #[test]
    fn per_variant_gives_one_group_for_each_variant() {
        let vars = [(0, 10), (0, 20), (1, 5)];
        let (of_each_var, groups) =
            walk_over(JackknifeGroups::PerVariant, &vars).expect("the walk over the variants");

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
        let (of_each_var, groups) =
            walk_over(JackknifeGroups::None, &vars).expect("the walk over the variants");

        assert_eq!(of_each_var, [None, None, None]);
        assert!(groups.is_empty());
    }

    /// A variant whose position is below the position of the variant
    /// before it on the same chromosome is refused where the groups are
    /// stretches of a chromosome: the cut compares its position with the
    /// first position of the group being filled, so it would join that
    /// group instead of starting one and the groups would not be the
    /// stretches the user asked for. The message names the chromosome and
    /// the two positions.
    #[test]
    fn a_variant_whose_position_goes_back_is_refused() {
        let vars = [(0, 1000), (0, 2000), (0, 1500)];
        let refused = walk_over(JackknifeGroups::OfBasePairs(5_000), &vars);

        let Err(error) = refused else {
            panic!("the variant that goes back was taken");
        };
        assert!(
            matches!(&error, Error::JackknifeGroupsVariantGoesBack { chrom, pos, before }
                if chrom == "chr1" && *pos == 1500 && *before == 2000),
            "{error:?}"
        );
    }

    /// A variant of a chromosome that the variant before it had left is
    /// refused the same way: its chromosome would be cut into a second run
    /// of groups over the stretch the first run covered, and the groups
    /// would overlap. The position check does not catch it, since the two
    /// positions are of different chromosomes.
    #[test]
    fn a_chromosome_that_comes_back_is_refused() {
        let vars = [(0, 1000), (1, 1000), (0, 2000)];
        let refused = walk_over(JackknifeGroups::OfBasePairs(5_000), &vars);

        let Err(error) = refused else {
            panic!("the chromosome that comes back was taken");
        };
        assert!(
            matches!(&error, Error::JackknifeGroupsChromComesBack { chrom, pos, before_chrom, before }
                if chrom == "chr1" && *pos == 2000 && before_chrom == "chr2" && *before == 1000),
            "{error:?}"
        );
    }

    /// What the order asks for and no more: two variants at one position
    /// are taken, as they are by the linkage disequilibrium filter of
    /// `docs/specs/filters.md`, and a chromosome that has not been read
    /// before starts its groups wherever its first position falls.
    #[test]
    fn variants_at_one_position_and_a_new_chromosome_are_taken() {
        let vars = [(0, 1000), (0, 1000), (1, 500), (1, 500)];
        let (of_each_var, groups) = walk_over(JackknifeGroups::OfBasePairs(5_000), &vars)
            .expect("the walk over the variants");

        assert_eq!(of_each_var, [Some(0), Some(0), Some(1), Some(1)]);
        assert_eq!(groups.len(), 2);
    }

    /// The order of the variants is asked for where the cut depends on it,
    /// which is a length of base pairs alone. Each variant is its own group
    /// with `PerVariant` whatever order they come in, so leaving a group out
    /// leaves that one variant out; with no groups the pass asks for
    /// neither the chromosome nor the position. The variants here are the
    /// ones the length above refuses, one going back and one of a
    /// chromosome that had been left.
    #[test]
    fn per_variant_and_no_groups_take_the_variants_in_any_order() {
        let vars = [(0, 2000), (0, 1500), (1, 1000), (0, 3000)];

        let (of_each_var, groups) =
            walk_over(JackknifeGroups::PerVariant, &vars).expect("the walk over the variants");
        assert_eq!(of_each_var, [Some(0), Some(1), Some(2), Some(3)]);
        assert_eq!(groups.len(), 4);

        let (of_each_var, groups) =
            walk_over(JackknifeGroups::None, &vars).expect("the walk over the variants");
        assert_eq!(of_each_var, [None, None, None, None]);
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
    ///
    /// What is refused is the ploidy the reader of the pass says its
    /// genotypes hold, and not an argument of a statistic of one variant,
    /// so the message is of the variants that were read: a user of
    /// `calc_pop_dists` writes no ploidy anywhere.
    #[test]
    fn a_ploidy_of_zero_or_above_the_largest_one_is_an_error() {
        for ploidy in [0, 256, usize::MAX] {
            let error = match PopDistPerVar::new(ploidy, 20) {
                Ok(_) => panic!("the ploidy {ploidy} was taken"),
                Err(error) => error,
            };
            assert!(
                matches!(&error, Error::PopDistsPloidyOutOfRange { ploidy: found, largest }
                    if *found == ploidy && *largest == 255),
                "{ploidy}: {error:?}"
            );
            let said = error.to_string();
            assert!(!said.contains("statistic of one variant"), "{said}");
            assert!(said.contains("ploidy"), "{said}");
        }
        assert!(PopDistPerVar::new(1, 20).is_ok());
        assert!(PopDistPerVar::new(255, 20).is_ok());
    }

    /// The populations alone can make more pairs than the machine counts,
    /// which is 93000 of them where a `usize` is 32 bits, as it is in
    /// wasm. It is found before a variant is read, so the message names
    /// the populations and not the resampling groups, which are none yet.
    #[test]
    fn too_many_populations_are_refused_without_a_word_about_the_groups() {
        let said = Error::PopDistsOfTooManyPops { num_pops: 93000 }.to_string();
        assert!(said.contains("93000"), "{said}");
        assert!(said.contains("populations"), "{said}");
        assert!(!said.contains("group"), "{said}");
    }

    /// The sums of a resampling group hold one place for each pair of the
    /// populations, and a variant whose pairs are more than those places
    /// would have the rest of them dropped: the pass would give a number
    /// for every pair, with the pairs after the last place counted at some
    /// variants and not at others.
    #[test]
    fn a_variant_whose_pairs_are_more_than_the_sums_hold_is_refused() {
        let gts = &WORKED_EXAMPLE[0];
        let per_var = match PopDistPerVar::new(2, 1) {
            Ok(per_var) => per_var,
            Err(error) => panic!("{error}"),
        };
        // Three populations make three pairs, and these sums hold two.
        let counts = [
            counts_of(gts, 2, &[0, 1]),
            counts_of(gts, 2, &[2, 3]),
            counts_of(gts, 2, &[4, 5]),
        ];
        let mut of_the_group = vec![PairSums::default(); 2];
        let refused = add_the_pairs(&counts, &per_var, &mut of_the_group);
        let error = match refused {
            Ok(()) => panic!("the variant was added into sums of two pairs"),
            Err(error) => error,
        };
        assert!(
            matches!(&error, Error::PopDistSumsOfAnotherSize { num_pops, num_pairs }
                if *num_pops == 3 && *num_pairs == 2),
            "{error:?}"
        );
        let mut of_the_three = vec![PairSums::default(); 3];
        assert!(add_the_pairs(&counts, &per_var, &mut of_the_three).is_ok());
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

    /// A reader over one of the panels of `tests/reference/`, of the ploidy
    /// 2 and in blocks of `num_vars_per_block` variants, or of the size the
    /// reader chooses when that is `None`.
    fn reader_of_the_panel(name: &str, num_vars_per_block: Option<usize>) -> impl BlockReader {
        let options = VcfOptions {
            ploidy: 2,
            num_vars_per_block,
            ..VcfOptions::default()
        };
        VcfReader::from_path(&reference(name), options).expect("the reader of the panel")
    }

    /// The six sums of every pair of the three populations of `pops_file`
    /// over the variants of the panel `name`, cut into the groups `how`
    /// says and read in blocks of `num_vars_per_block` variants.
    ///
    /// It is the pass without the checks `calc_pop_dist_sums` makes of what
    /// it was asked for, because the two runs of ADMIXTOOLS 2 the standard
    /// errors are compared with cut the panel into 12 and into 6 groups and
    /// a user is asked for 20. Everything else it does is the same.
    fn sums_of_the_panel(
        name: &str,
        pops_file: &str,
        how: JackknifeGroups,
        min_num_individuals: u32,
        num_vars_per_block: Option<usize>,
    ) -> PopDistSums {
        let mut reader = reader_of_the_panel(name, num_vars_per_block);
        let pops = pops_of_the_file(pops_file, reader.individuals(), &["p0", "p1", "p2"]);
        let options = PopDistOptions {
            min_num_individuals,
            groups: how,
        };
        sums_of_the_pass(&mut reader, &pops, &options).expect("the sums of the panel")
    }

    /// The six sums of every pair of the three populations of the biallelic
    /// panel over its first variant alone: its first block with every row
    /// but the first taken out, given to the pass by a reader of the tests.
    fn sums_of_the_first_variant_of_the_panel() -> PopDistSums {
        let mut of_the_panel = reader_of_the_panel("dists/panel.vcf.gz", Some(100));
        let mut block = of_the_panel
            .next_block()
            .expect("a block of the panel")
            .expect("the first block of the panel");
        let of_the_first_variant: Vec<bool> = (0..block.num_vars).map(|var| var == 0).collect();
        block
            .retain_vars(&of_the_first_variant)
            .expect("the first variant of the block");
        let pops = pops_of_the_file(
            "stats/panel_pops.txt",
            of_the_panel.individuals(),
            &["p0", "p1", "p2"],
        );
        let mut reader = GivenBlocks::of(of_the_panel.individuals().to_vec(), vec![block], false);
        calc_pop_dist_sums(
            &mut reader,
            &pops,
            &PopDistOptions {
                min_num_individuals: 20,
                groups: JackknifeGroups::None,
            },
        )
        .expect("the sums of one variant of the panel")
    }

    /// The six sums of the one pair of the worked example over its four
    /// variants, cut into the groups `how` says, the variants lying 1000
    /// base pairs apart on one chromosome.
    ///
    /// It is the pass without the checks of `calc_pop_dist_sums` too: four
    /// variants fall into four groups at the most, and a user is asked for
    /// 20.
    fn sums_of_the_worked_example(how: JackknifeGroups) -> PopDistSums {
        let mut reader = GivenBlocks::of_the_worked_example(2, true, false);
        let pops = pops_of_the_worked_example();
        let options = PopDistOptions {
            min_num_individuals: 1,
            groups: how,
        };
        sums_of_the_pass(&mut reader, &pops, &options).expect("the sums of the worked example")
    }

    /// The two populations of the worked example, the first three
    /// individuals of its reader and the last three.
    fn pops_of_the_worked_example() -> Pops {
        let individuals: Vec<String> = (0..6).map(|number| format!("i{number}")).collect();
        let named = [
            ("pop1".to_owned(), individuals[..3].to_vec()),
            ("pop2".to_owned(), individuals[3..].to_vec()),
        ];
        Pops::from_names(&named, &individuals).expect("the populations of the worked example")
    }

    /// The genotypes of the two variants of "What pyNei does that is odd,
    /// and what popnei does instead" of `docs/specs/dists.md`: 4 diploid
    /// individuals, the two alleles of each after those of the one before,
    /// with -1 for an allele that was not called, pop1 the first two
    /// individuals and pop2 the last two.
    ///
    /// Each population has exactly one called genotype at the first
    /// variant, which is what no panel of popnei reaches and what a
    /// `min_num_individuals` of 1 lets through as far as the rule that
    /// drops it. The second variant carries the pair on its own, and its
    /// Jost's D is the 0.25 pyNei gives for the two together.
    const ONE_CALLED_GENOTYPE_EACH: [[i8; 8]; 2] =
        [[0, 0, -1, -1, 1, 1, -1, -1], [0, 0, 0, 1, 0, 1, 1, 1]];

    /// The sums of the one pair of two populations of two individuals over
    /// `of_the_vars`, one block of them in the order given, at a
    /// `min_num_individuals` of 1 and with no resampling groups.
    fn sums_of_the_pair_of_two_over(of_the_vars: &[[i8; 8]]) -> PopDistSums {
        let individuals: Vec<String> = (0..4).map(|number| format!("i{number}")).collect();
        let named = [
            ("pop1".to_owned(), individuals[..2].to_vec()),
            ("pop2".to_owned(), individuals[2..].to_vec()),
        ];
        let pops = Pops::from_names(&named, &individuals).expect("the two populations of two");
        let block = Block {
            num_vars: of_the_vars.len(),
            num_individuals: 4,
            ploidy: 2,
            gts: of_the_vars.iter().flatten().copied().collect(),
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        };
        let mut reader = GivenBlocks::of(individuals, vec![block], false);
        calc_pop_dist_sums(
            &mut reader,
            &pops,
            &PopDistOptions {
                min_num_individuals: 1,
                groups: JackknifeGroups::None,
            },
        )
        .expect("the sums of the two populations of two")
    }

    /// The six sums of the one pair of two populations of 6 over the 200
    /// variants of the haploid file of `tests/reference/dists/`, read at a
    /// ploidy of 1, at `min_num_individuals` called genotypes and with no
    /// resampling groups.
    ///
    /// Its 12 individuals are `h00` to `h11`, the first six of them the
    /// first population and the last six the second.
    fn sums_of_the_haploid_file(min_num_individuals: u32) -> PopDistSums {
        let options = VcfOptions {
            ploidy: 1,
            ..VcfOptions::default()
        };
        let mut reader = VcfReader::from_path(&reference("dists/haploid.vcf.gz"), options)
            .expect("the reader of the haploid file");
        let individuals = reader.individuals().to_vec();
        let named = [
            ("pop1".to_owned(), individuals[..6].to_vec()),
            ("pop2".to_owned(), individuals[6..].to_vec()),
        ];
        let pops = Pops::from_names(&named, &individuals).expect("the two haploid populations");
        calc_pop_dist_sums(
            &mut reader,
            &pops,
            &PopDistOptions {
                min_num_individuals,
                groups: JackknifeGroups::None,
            },
        )
        .expect("the sums of the haploid file")
    }

    /// The genotypes of the 25 variants of "How it is verified" of
    /// `docs/specs/dists.md` where the two populations are fixed for the
    /// same allele: 4 diploid individuals, the two alleles of each after
    /// those of the one before, pop1 the first two and pop2 the last two,
    /// 24 variants with every genotype `0/0` and a last one where pop1
    /// holds the allele 0 and pop2 the allele 1.
    fn fixed_for_the_same_allele() -> Vec<[i8; 8]> {
        let mut of_the_vars = vec![[0; 8]; 24];
        of_the_vars.push([0, 0, 0, 0, 1, 1, 1, 1]);
        of_the_vars
    }

    /// The sums of the one pair of two populations of two individuals over
    /// `of_the_vars`, one block of them in the order given, at
    /// `min_num_individuals` called genotypes and with each variant its own
    /// resampling group.
    ///
    /// The variants lie 1000 base pairs apart on one chromosome. It is the
    /// pass without the checks `calc_pop_dist_sums` makes, which ask a user
    /// for 20 groups at the least.
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "the variants of a fixture, at the positions 1000 and up"
    )]
    fn sums_of_the_pair_of_two_per_variant(
        of_the_vars: &[[i8; 8]],
        min_num_individuals: u32,
    ) -> PopDistSums {
        let individuals: Vec<String> = (0..4).map(|number| format!("i{number}")).collect();
        let named = [
            ("pop1".to_owned(), individuals[..2].to_vec()),
            ("pop2".to_owned(), individuals[2..].to_vec()),
        ];
        let pops = Pops::from_names(&named, &individuals).expect("the two populations of two");
        let block = Block {
            num_vars: of_the_vars.len(),
            num_individuals: 4,
            ploidy: 2,
            gts: of_the_vars.iter().flatten().copied().collect(),
            chrom: Some(vec![0; of_the_vars.len()]),
            pos: Some(
                (0..of_the_vars.len())
                    .map(|var| 1000 * (var + 1) as u64)
                    .collect(),
            ),
            id: None,
            alleles: None,
            qual: None,
        };
        let mut reader = GivenBlocks::of(individuals, vec![block], false);
        sums_of_the_pass(
            &mut reader,
            &pops,
            &PopDistOptions {
                min_num_individuals,
                groups: JackknifeGroups::PerVariant,
            },
        )
        .expect("the sums of the two populations of two")
    }

    /// The genotypes of one variant of 10 diploid individuals, the two
    /// alleles of each after those of the one before, where each of the two
    /// populations of five holds the alleles 0, 1, 2 and 3 at the counts 2,
    /// 4, 3 and 1 out of its 10 called alleles.
    ///
    /// Both populations have the same four frequencies, so each square root
    /// of a product is one of them and their sum is 0.2 + 0.4 + 0.3 + 0.1,
    /// which a f64 adds to 1 + 2.2e-16 over the one variant that counted.
    const THE_SAME_FOUR_FREQUENCIES: [i8; 20] =
        [0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3];

    /// The sums of the one pair of two populations of five individuals over
    /// the one variant `of_the_var`, at a `min_num_individuals` of 1 and
    /// with no resampling groups.
    fn sums_of_the_pair_of_five_over(of_the_var: &[i8; 20]) -> PopDistSums {
        let individuals: Vec<String> = (0..10).map(|number| format!("i{number}")).collect();
        let named = [
            ("pop1".to_owned(), individuals[..5].to_vec()),
            ("pop2".to_owned(), individuals[5..].to_vec()),
        ];
        let pops = Pops::from_names(&named, &individuals).expect("the two populations of five");
        let block = Block {
            num_vars: 1,
            num_individuals: 10,
            ploidy: 2,
            gts: of_the_var.to_vec(),
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        };
        let mut reader = GivenBlocks::of(individuals, vec![block], false);
        calc_pop_dist_sums(
            &mut reader,
            &pops,
            &PopDistOptions {
                min_num_individuals: 1,
                groups: JackknifeGroups::None,
            },
        )
        .expect("the sums of the two populations of five")
    }

    /// A reader of the tests that gives the blocks it was built with, and
    /// then the error of a bgzipped file with no mark of its end when it
    /// was built to fail, which is how an error of a reader reaches the
    /// pass.
    #[derive(Debug)]
    struct GivenBlocks {
        individuals: Vec<String>,
        chroms: ChromTable,
        /// The blocks it has not given yet, the next one last.
        left: Vec<Block>,
        /// What it was last asked to fill, which a test reads to see which
        /// fields the pass asked for.
        needs: Needs,
        /// Whether it fails where it would say that its blocks are over.
        fails_at_the_end: bool,
    }

    impl GivenBlocks {
        /// The reader over `blocks`, of the individuals `individuals` and
        /// of the ploidy 2, which is the ploidy of every block of these
        /// tests.
        ///
        /// Its table of chromosome names holds `chr1` and `chr2`, the two
        /// of the biallelic panel, so that a block of the panel given back
        /// to it carries the names its own reader gave the numbers.
        fn of(individuals: Vec<String>, blocks: Vec<Block>, fails_at_the_end: bool) -> GivenBlocks {
            let mut chroms = ChromTable::new();
            chroms.intern("chr1");
            chroms.intern("chr2");
            let mut left = blocks;
            left.reverse();
            GivenBlocks {
                individuals,
                chroms,
                left,
                needs: Needs::ALL,
                fails_at_the_end,
            }
        }

        /// The reader over the four variants of the worked example, in
        /// blocks of `num_vars_per_block` variants, of six individuals
        /// named `i0` to `i5`.
        ///
        /// The variants lie 1000 base pairs apart on one chromosome, and
        /// with `with_the_positions` false the blocks hold neither the
        /// chromosome nor the position, which is what a reader that was
        /// asked for them and gave none leaves.
        fn of_the_worked_example(
            num_vars_per_block: usize,
            with_the_positions: bool,
            fails_at_the_end: bool,
        ) -> GivenBlocks {
            GivenBlocks::of(
                (0..6).map(|number| format!("i{number}")).collect(),
                blocks_of_the_worked_example(num_vars_per_block, with_the_positions),
                fails_at_the_end,
            )
        }

        /// The fields the pass last asked it to fill.
        fn needs(&self) -> Needs {
            self.needs
        }
    }

    impl BlockReader for GivenBlocks {
        fn next_block(&mut self) -> Result<Option<Block>> {
            match self.left.pop() {
                Some(block) => Ok(Some(block)),
                None if self.fails_at_the_end => Err(Error::VcfBgzipEndMissing),
                None => Ok(None),
            }
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

        fn set_needs(&mut self, needs: Needs) {
            self.needs = needs;
        }

        fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
            Vec::new()
        }
    }

    /// The blocks of `num_vars_per_block` variants that hold the four
    /// variants of the worked example, the first block first.
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "the four variants of the worked example, at the positions 1000 to 4000"
    )]
    fn blocks_of_the_worked_example(
        num_vars_per_block: usize,
        with_the_positions: bool,
    ) -> Vec<Block> {
        WORKED_EXAMPLE
            .chunks(num_vars_per_block)
            .enumerate()
            .map(|(block, of_the_block)| Block {
                num_vars: of_the_block.len(),
                num_individuals: 6,
                ploidy: 2,
                gts: of_the_block.iter().flatten().copied().collect(),
                chrom: with_the_positions.then(|| vec![0; of_the_block.len()]),
                pos: with_the_positions.then(|| {
                    (0..of_the_block.len())
                        .map(|var| 1000 * (block * num_vars_per_block + var + 1) as u64)
                        .collect()
                }),
                id: None,
                alleles: None,
                qual: None,
            })
            .collect()
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
    /// `docs/specs/dists.md`: the sums over its five variants are 2.902778
    /// for H_b and 1.883333 for H_w, so F_ST is 1.019444 / 2.902778 and f_2
    /// is 1.019444 / 5. Both are ratios of the sums and not means of the
    /// per variant ratios, which for F_ST is undefined at the third
    /// variant.
    #[test]
    fn the_worked_example_has_the_fst_and_the_f2_of_the_spec() {
        let sums = sums_of_the_worked_example(JackknifeGroups::PerVariant);

        assert_it_is_within(
            sums.measure(PopDistMeasure::Fst, 0, 1),
            0.351196,
            1e-6,
            "the F_ST of the worked example",
        );
        assert_it_is_within(
            sums.measure(PopDistMeasure::F2, 0, 1),
            0.203889,
            1e-6,
            "the f_2 of the worked example",
        );
        assert_eq!(sums.num_vars_of(0, 1), Some(5));
        assert_eq!(sums.num_vars_of(1, 0), sums.num_vars_of(0, 1));
    }

    /// Jost's D, Nei's G_ST and the standardized G''_ST of the worked
    /// example of "How it is verified" of `docs/specs/dists.md`, the three
    /// measures that are ratios of the means of the corrected sums: over
    /// its five variants the mean corrected H_S is 0.357143 and the mean
    /// corrected H_T is 0.468849, so D is 2 * 0.111706 / 0.642857, G_ST is
    /// 0.111706 / 0.468849, and G''_ST is 0.223412 / (0.580556 * 0.642857).
    /// pyNei gives the same D for the same genotypes, 0.34753086.
    #[test]
    fn the_worked_example_has_the_dest_the_gst_and_the_gst_standardized_of_the_spec() {
        let sums = sums_of_the_worked_example(JackknifeGroups::PerVariant);

        assert_it_is_within(
            sums.measure(PopDistMeasure::Dest, 0, 1),
            0.347531,
            1e-6,
            "the Jost's D of the worked example",
        );
        assert_it_is_within(
            sums.measure(PopDistMeasure::Gst, 0, 1),
            0.238256,
            1e-6,
            "the G_ST of the worked example",
        );
        assert_it_is_within(
            sums.measure(PopDistMeasure::GstStandardized, 0, 1),
            0.598618,
            1e-6,
            "the G''_ST of the worked example",
        );
    }

    /// The chord distance and Nei's D_A of the worked example of "How it is
    /// verified" of `docs/specs/dists.md`, the two measures that come out
    /// of the sum of the square roots of the products of the frequencies:
    /// over its five variants that sum is 3.707290, so D_A is
    /// 1 - 3.707290 / 5 and the chord distance is the square root of it.
    #[test]
    fn the_worked_example_has_the_chord_and_the_da_of_the_spec() {
        let sums = sums_of_the_worked_example(JackknifeGroups::PerVariant);

        assert_it_is_within(
            sums.measure(PopDistMeasure::Da, 0, 1),
            0.258542,
            1e-6,
            "the Nei's D_A of the worked example",
        );
        assert_it_is_within(
            sums.measure(PopDistMeasure::Chord, 0, 1),
            0.508470,
            1e-6,
            "the chord distance of the worked example",
        );
    }

    /// A variant where both populations have exactly one called genotype
    /// counts for no measure of the pair, which takes a
    /// `min_num_individuals` of 1 and which "Variants that do not count,
    /// populations with little data, and negative values" of
    /// `docs/specs/dists.md` asks for: the harmonic mean of the two counts
    /// is 1, and the correction of H_S divides by 1 - 1.
    ///
    /// The two variants are the ones "What pyNei does that is odd, and what
    /// popnei does instead" of the Jost's D item measured pyNei at commit
    /// ef0ca6e on, four individuals in two populations of two with
    /// `min_num_samples=1`: pyNei drops the first variant through a NaN
    /// that `numpy.nansum` leaves out of both sums and out of the count,
    /// printing two RuntimeWarnings, and gives D as 0.25, which is the D of
    /// the second variant alone. popnei drops it by the rule above, with no
    /// NaN made and nothing printed, and every measure of the two variants
    /// is the measure of the second alone to the bit.
    #[test]
    fn a_variant_where_both_pops_have_one_called_genotype_counts_for_no_measure() {
        let of_both = sums_of_the_pair_of_two_over(&ONE_CALLED_GENOTYPE_EACH);
        let of_the_second = sums_of_the_pair_of_two_over(&ONE_CALLED_GENOTYPE_EACH[1..]);

        assert_eq!(of_both.num_vars(), 2);
        assert_eq!(of_both.num_vars_of(0, 1), Some(1));
        assert_eq!(of_the_second.num_vars(), 1);
        assert_eq!(of_the_second.num_vars_of(0, 1), Some(1));
        for measure in PopDistMeasure::THAT_HAVE_A_VALUE {
            assert_eq!(
                of_both.measure(measure, 0, 1),
                of_the_second.measure(measure, 0, 1),
                "the {} over the two variants and over the second alone",
                measure.name()
            );
        }
        assert_it_is_the_same_number(
            of_both.measure(PopDistMeasure::Dest, 0, 1),
            0.25,
            "the Jost's D of the two variants, which pyNei gives as 0.25",
        );
    }

    /// Jost's D, Nei's G_ST and the standardized G''_ST have no value at a
    /// ploidy of 1, and the other four measures have one, which "Variants
    /// that do not count" of `docs/specs/dists.md` asks for.
    ///
    /// The three are built from the corrected H_S and H_T, which raise the
    /// allele frequencies to the ploidy: at a ploidy of 1 those sums are
    /// the frequencies themselves, which add to 1, and a haploid genotype
    /// is never heterozygous, so H_S, H_T and the observed heterozygosity
    /// are 0 by their definitions and the sums of a pass hold nothing but
    /// the residue of adding frequencies that a f64 does not bring to
    /// exactly 1. Before this rule, popnei gave a G_ST of 0.250299 for this
    /// pair at a `min_num_individuals` of 3 and of 0.234741 at 4, 6 in 100
    /// apart, where the chord distance moves 2 in 1000, and a Jost's D of
    /// 1.4e-17. The two thresholds are read here because one of them alone
    /// would not show that the number moves with what it is over.
    #[test]
    fn the_three_measures_of_the_corrected_diversities_have_no_value_at_a_ploidy_of_one() {
        for (min_num_individuals, num_vars) in [(3, 200), (4, 198)] {
            let sums = sums_of_the_haploid_file(min_num_individuals);

            assert_eq!(sums.num_vars_of(0, 1), Some(num_vars));
            for measure in [
                PopDistMeasure::Dest,
                PopDistMeasure::Gst,
                PopDistMeasure::GstStandardized,
            ] {
                let named = measure.name();
                assert_eq!(
                    sums.measure(measure, 0, 1),
                    None,
                    "the {named} of a haploid pass at {min_num_individuals} called genotypes"
                );
            }
            for measure in [
                PopDistMeasure::Fst,
                PopDistMeasure::F2,
                PopDistMeasure::Chord,
                PopDistMeasure::Da,
            ] {
                let named = measure.name();
                assert!(
                    sums.measure(measure, 0, 1).is_some_and(f64::is_finite),
                    "the {named} of a haploid pass at {min_num_individuals} called genotypes"
                );
            }
        }
    }

    /// Two populations fixed for the same allele at every variant that
    /// counted for them have no F_ST, no G_ST and no G''_ST: their sum of
    /// H_b and their mean corrected H_T are both 0, and each of the three
    /// divides by one of the two. The four that do not, f_2, Jost's D, the
    /// chord distance and Nei's D_A, are 0 there, which "Variants that do
    /// not count" of `docs/specs/dists.md` asks for: two populations that
    /// hold the same one allele everywhere are as near as a pair can be.
    ///
    /// A value of NaN instead would be a distance with a count of variants
    /// above 0 beside it, and the test below is what it costs: a NaN cannot
    /// be told from a number by the jackknife that reads it.
    #[test]
    fn two_pops_fixed_for_the_same_allele_have_no_fst_no_gst_and_no_gst_standardized() {
        let sums = sums_of_the_pair_of_two_over(&[[0; 8], [0; 8]]);

        assert_eq!(sums.num_vars_of(0, 1), Some(2));
        for measure in [
            PopDistMeasure::Fst,
            PopDistMeasure::Gst,
            PopDistMeasure::GstStandardized,
        ] {
            assert_eq!(
                sums.measure(measure, 0, 1),
                None,
                "the {} of two populations fixed for the same allele",
                measure.name()
            );
        }
        for measure in [
            PopDistMeasure::F2,
            PopDistMeasure::Dest,
            PopDistMeasure::Chord,
            PopDistMeasure::Da,
        ] {
            assert_eq!(
                sums.measure(measure, 0, 1),
                Some(0.0),
                "the {} of two populations fixed for the same allele",
                measure.name()
            );
        }
    }

    /// The 25 variants of "How it is verified" of `docs/specs/dists.md`,
    /// each its own group: the two populations are fixed for the same
    /// allele at 24 of them and hold a different one at the last, so that
    /// one group carries the whole sum of H_b and the whole sum of the
    /// corrected H_T. Left out, it leaves F_ST, G_ST and G''_ST without a
    /// value, and the three have no standard error at all, although each of
    /// them is 1 over the 25 variants.
    ///
    /// The other four have a value with every group left out, so they have
    /// a standard error over the same 25 groups: f_2 is 0.04 with one of
    /// 0.04, Jost's D and Nei's D_A the same two numbers, and the chord
    /// distance 0.2 with 0.195959. Before the divisors were guarded, the
    /// group left out gave the three a pseudo-value of NaN, and what a user
    /// read for them was a standard error of NaN with a value of 1 beside
    /// it.
    #[test]
    fn a_measure_that_a_group_left_out_leaves_without_a_value_has_no_standard_error() {
        let sums = sums_of_the_pair_of_two_per_variant(&fixed_for_the_same_allele(), 2);

        assert_eq!(sums.num_vars_of(0, 1), Some(25));
        assert_eq!(sums.groups().len(), 25);
        for measure in [
            PopDistMeasure::Fst,
            PopDistMeasure::Gst,
            PopDistMeasure::GstStandardized,
        ] {
            let named = measure.name();
            assert_it_is_within(
                sums.measure(measure, 0, 1),
                1.0,
                1e-12,
                &format!("the {named}"),
            );
            assert_eq!(
                sums.standard_error(measure, 0, 1),
                None,
                "the standard error of the {named}"
            );
        }
        for (measure, value, standard_error) in [
            (PopDistMeasure::F2, 0.04, 0.04),
            (PopDistMeasure::Dest, 0.04, 0.04),
            (PopDistMeasure::Da, 0.04, 0.04),
            (PopDistMeasure::Chord, 0.2, 0.195959),
        ] {
            let named = measure.name();
            assert_it_is_within(
                sums.measure(measure, 0, 1),
                value,
                1e-6,
                &format!("the {named}"),
            );
            assert_it_is_within(
                sums.standard_error(measure, 0, 1),
                standard_error,
                1e-6,
                &format!("the standard error of the {named}"),
            );
        }
    }

    /// Jost's D and the standardized G''_ST both divide by 1 - the mean
    /// corrected H_S, so a pair whose mean corrected H_S came to exactly 1
    /// has neither. G_ST divides by the mean corrected H_T instead and has
    /// a value there, which is what "Variants that do not count" of
    /// `docs/specs/dists.md` says of the three.
    ///
    /// No genotypes of a panel reach a mean corrected H_S of exactly 1, so
    /// the sums are written here as a pass would have left them: two
    /// variants whose corrected H_S added to 2.
    #[test]
    fn the_dest_and_the_gst_standardized_of_a_pair_whose_mean_corrected_h_s_is_one_have_no_value() {
        let sums = PopDistSums::of_the_pass(
            2,
            2,
            2,
            Vec::new(),
            vec![PairSums {
                h_b: 1.5,
                h_w: 1.0,
                sqrt_of_the_products: 1.0,
                corrected_h_s: 2.0,
                corrected_h_t: 3.0,
                num_vars: 2,
            }],
        );

        assert_eq!(sums.num_vars_of(0, 1), Some(2));
        assert_eq!(sums.measure(PopDistMeasure::Dest, 0, 1), None);
        assert!(
            sums.measure(PopDistMeasure::Gst, 0, 1).is_some(),
            "the G_ST of a pair whose mean corrected H_S is 1"
        );
        assert_eq!(
            sums.measure(PopDistMeasure::GstStandardized, 0, 1),
            None,
            "the G''_ST of a pair whose mean corrected H_S is 1"
        );
    }

    /// The groups the variants were cut into change no measure: the sums of
    /// a pair are added group by group and the division happens once,
    /// whether the five variants of the worked example fall in five groups,
    /// in one or in none.
    #[test]
    fn the_measures_are_the_same_whatever_the_groups() {
        let of_five = sums_of_the_worked_example(JackknifeGroups::PerVariant);
        let of_one = sums_of_the_worked_example(JackknifeGroups::OfBasePairs(100_000));
        let of_none = sums_of_the_worked_example(JackknifeGroups::None);

        assert_eq!(of_five.groups().len(), 5);
        assert_eq!(of_one.groups().len(), 1);
        assert!(of_none.groups().is_empty());
        for measure in [PopDistMeasure::Fst, PopDistMeasure::F2] {
            let found = of_five.measure(measure, 0, 1);
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

        for (var, f2) in [0.388889, 0.272222, 0.0, -0.1, 0.458333]
            .into_iter()
            .enumerate()
        {
            assert_it_is_within(
                sums.f2_of_group(var, 0, 1),
                f2,
                1e-6,
                &format!("the f_2 of the group of the variant {var}"),
            );
        }
        assert_eq!(sums.f2_of_group(5, 0, 1), None);
    }

    /// A pass of one variant has an F_ST equal to that variant's, so the
    /// per variant F_ST that plink2 reports for `var0000` of the biallelic
    /// panel, 0.332322 for p0 and p1, is the F_ST of a pass over that one
    /// variant.
    #[test]
    fn the_fst_of_a_pass_of_one_variant_is_the_fst_of_that_variant() {
        let sums = sums_of_the_first_variant_of_the_panel();

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
            let sums = sums_of_the_panel(panel, pops_file, JackknifeGroups::None, 20, None);
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

    /// Jost's D, Nei's G_ST and the standardized G''_ST of the three pairs
    /// of both panels against `pairwise_D`, `pairwise_Gst_Nei` and
    /// `pairwise_Gst_Hedrick` of mmod 1.3.3 under R 4.6.1, which
    /// `tests/reference/pop_dists/panel.mmod.tsv` and `micro.mmod.tsv`
    /// hold. `pairwise_Gst_Hedrick` computes the standardized G''_ST of
    /// Meirmans and Hedrick (2011) and not the G'_ST its name suggests,
    /// which "How it is verified" of the G_ST item of `docs/specs/dists.md`
    /// shows from its source.
    ///
    /// mmod computes another estimator of the same three quantities: its
    /// `HsHt` leaves the observed heterozygosity term out of both
    /// corrections and uses 2n/(2n - 1) where popnei, which is pyNei and
    /// Nei and Chesser (1983), uses n/(n - 1) and subtracts H_obs/(2n). So
    /// the check is an agreement and not an equality, and its tolerance is
    /// the 5e-4 absolute of the two items of the spec, which the furthest
    /// of these eighteen numbers, a G''_ST of the multiallelic panel, is
    /// 4.7e-4 within. What pins the estimator is the comparison with pyNei,
    /// exact to 1e-12 relative, which is made where a user sees the number
    /// and not here. Tightening this tolerance, or moving the arithmetic
    /// towards mmod, breaks that comparison.
    #[test]
    fn the_dest_the_gst_and_the_gst_standardized_of_both_panels_agree_with_mmod() {
        let panels = [
            (
                "dists/panel.vcf.gz",
                "stats/panel_pops.txt",
                [
                    [0.0634704859, 0.0612312792, 0.0656071139],
                    [0.0553896690, 0.0541614926, 0.0579922831],
                    [0.1617736271, 0.1576967932, 0.1680418436],
                ],
            ),
            (
                "pop_dists/micro.vcf.gz",
                "pop_dists/micro_pops.txt",
                [
                    [0.1662094246, 0.1819580763, 0.1816462134],
                    [0.0331789783, 0.0359529575, 0.0362672006],
                    [0.2197612680, 0.2387386980, 0.2389275805],
                ],
            ),
        ];
        for (panel, pops_file, of_mmod) in panels {
            let sums = sums_of_the_panel(panel, pops_file, JackknifeGroups::None, 20, None);
            for (row, measure) in [
                PopDistMeasure::Dest,
                PopDistMeasure::Gst,
                PopDistMeasure::GstStandardized,
            ]
            .into_iter()
            .enumerate()
            {
                for (pair, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
                    assert_it_is_within(
                        sums.measure(measure, i, j),
                        of_mmod[row][pair],
                        5e-4,
                        &format!("the {} of the pair {i} {j} of {panel}", measure.name()),
                    );
                }
            }
        }
    }

    /// The chord distance of the three pairs of both panels against
    /// `dist.genpop(method = 2)` of adegenet 2.1.11 under R 4.6.1, which
    /// `tests/reference/pop_dists/panel.chord.tsv` and `micro.chord.tsv`
    /// hold. It is the form adegenet gives, the chord of the sphere of
    /// radius 1 divided by the square root of 2, which "What it gives" of
    /// the chord item of `docs/specs/dists.md` says popnei computes: the
    /// chord of that sphere itself is 1.414 times each of these six
    /// numbers.
    ///
    /// adegenet computes the same estimator, so the comparison is within
    /// 1e-12 relative and not the 5e-4 of the measures mmod checks. Each of
    /// the biallelic panel's three comes out as the same double adegenet
    /// prints, and the multiallelic panel's within 6.7e-16, which is the
    /// last bits of one. D_A is in no program and is checked as the square
    /// of what adegenet gives, which "How it is verified" of the item asks
    /// for; it is 3.9e-15 from that square at the furthest, since squaring
    /// doubles how far from adegenet the chord distance is.
    #[test]
    #[expect(
        clippy::excessive_precision,
        reason = "the numbers as tests/reference/pop_dists/panel.chord.tsv and micro.chord.tsv print them, which is one digit more than a f64 keeps"
    )]
    fn the_chord_and_the_da_of_both_panels_are_adegenets() {
        let panels = [
            (
                "dists/panel.vcf.gz",
                "stats/panel_pops.txt",
                [
                    0.18026704497001397,
                    0.17586558860911838,
                    0.17977447554045811,
                ],
            ),
            (
                "pop_dists/micro.vcf.gz",
                "pop_dists/micro_pops.txt",
                [
                    0.33853588707322202,
                    0.33760692274322368,
                    0.34958215447311330,
                ],
            ),
        ];
        for (panel, pops_file, of_adegenet) in panels {
            let sums = sums_of_the_panel(panel, pops_file, JackknifeGroups::None, 20, None);
            for (pair, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
                assert_it_is_the_same_number(
                    sums.measure(PopDistMeasure::Chord, i, j),
                    of_adegenet[pair],
                    &format!("the chord distance of the pair {i} {j} of {panel}"),
                );
                assert_it_is_the_same_number(
                    sums.measure(PopDistMeasure::Da, i, j),
                    of_adegenet[pair] * of_adegenet[pair],
                    &format!("the Nei's D_A of the pair {i} {j} of {panel}"),
                );
            }
        }
    }

    /// Two populations with the same allele frequencies at every variant
    /// have a chord distance of 0, which is the rule "What it gives" of the
    /// chord item of `docs/specs/dists.md` gives for a sum of the square
    /// roots that the rounding takes above the variants that counted. The
    /// four frequencies of the fixture add to 1 + 2.2e-16, so without that
    /// rule D_A would be -2.2e-16 and the chord distance, its square root,
    /// a NaN.
    #[test]
    fn two_pops_with_the_same_frequencies_have_a_chord_of_zero_and_not_a_nan() {
        let sums = sums_of_the_pair_of_five_over(&THE_SAME_FOUR_FREQUENCIES);
        let of_the_pair = sums.total_of(0).expect("the sums of the one pair");

        assert_eq!(sums.num_vars_of(0, 1), Some(1));
        assert!(
            of_the_pair.sqrt_of_the_products > of_the_pair.num_vars as f64,
            "the sum of the square roots is {} over {} variants, and the fixture is here because the rounding takes it above the count",
            of_the_pair.sqrt_of_the_products,
            of_the_pair.num_vars
        );
        assert_eq!(sums.measure(PopDistMeasure::Da, 0, 1), Some(0.0));
        assert_eq!(sums.measure(PopDistMeasure::Chord, 0, 1), Some(0.0));
    }

    /// f_2 and its jackknife standard error for the three pairs of the
    /// biallelic panel against ADMIXTOOLS 2.0.10, which computes f_2 with
    /// the same delete-m jackknife over the same groups, cut the same way.
    /// `expected` is the f_2 and the standard error of p0-p1, p0-p2 and
    /// p1-p2, which `tests/reference/pop_dists/panel.f2.tsv` and
    /// `panel.f2.uneven.tsv` hold to seventeen digits. The two agree to the
    /// last bits of a double, so the comparison is within 1e-12 relative.
    fn assert_the_f2_of_the_panel_is(how: JackknifeGroups, expected: [(f64, f64); 3]) {
        let sums = sums_of_the_panel("dists/panel.vcf.gz", "stats/panel_pops.txt", how, 20, None);

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
            None,
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
    /// of the measure, 0.004197 of a F_ST of 0.104962 for p0 and p1 over
    /// the 12 groups of 100 000 base pairs. Its arithmetic is the one
    /// ADMIXTOOLS checks.
    #[test]
    fn the_fst_has_a_standard_error_of_the_same_jackknife() {
        let sums = sums_of_the_panel(
            "dists/panel.vcf.gz",
            "stats/panel_pops.txt",
            JackknifeGroups::OfBasePairs(100_000),
            20,
            None,
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

    /// The pairs a caller asks for by number: two populations that are one,
    /// and one that is not a population of the pass, have no measure and no
    /// count of variants, whichever of the seven is asked for, while the
    /// one pair the worked example holds has every one of them. `measures`
    /// gives the pairs in the order of the distance vector, which for three
    /// populations is p0-p1, p0-p2 and p1-p2.
    #[test]
    fn a_pair_that_is_not_one_has_no_measure_and_no_count() {
        let sums = sums_of_the_worked_example(JackknifeGroups::PerVariant);

        assert_eq!(sums.num_pops(), 2);
        assert_eq!(sums.num_vars_of(0, 0), None);
        assert_eq!(sums.standard_error(PopDistMeasure::F2, 0, 0), None);
        for measure in PopDistMeasure::THAT_HAVE_A_VALUE {
            let named = measure.name();
            assert!(
                sums.measure(measure, 0, 1).is_some(),
                "the {named} of the one pair"
            );
            assert_eq!(sums.measure(measure, 0, 0), None, "the {named} of 0 and 0");
            assert_eq!(sums.measure(measure, 0, 2), None, "the {named} of 0 and 2");
        }
        let of_three = sums_of_the_first_variant_of_the_panel();
        let in_order: Vec<Option<f64>> = of_three.measures(PopDistMeasure::Fst).collect();
        assert_eq!(in_order.len(), 3);
        for (pair, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
            assert_eq!(in_order[pair], of_three.measure(PopDistMeasure::Fst, i, j));
        }
    }

    /// The four iterators over the pairs give what the accessors by
    /// number give for the pairs of the distance vector, in that order:
    /// they are what the binding crates read every array of a result from,
    /// so the values, the standard errors, the counts and the f_2 of the
    /// groups of one result cannot fall into different orders.
    ///
    /// The biallelic panel is the fixture because its three populations
    /// make three pairs whose numbers differ from one another, and the
    /// groups of 250 000 base pairs are 6, so `f2_of_every_group` gives 18
    /// values whose order would show as soon as it went by pairs and then
    /// by groups.
    #[test]
    fn the_iterators_over_the_pairs_are_in_the_order_of_the_distance_vector() {
        let sums = sums_of_the_panel(
            "dists/panel.vcf.gz",
            "stats/panel_pops.txt",
            JackknifeGroups::OfBasePairs(250_000),
            20,
            None,
        );
        let pairs = [(0, 1), (0, 2), (1, 2)];

        assert_eq!(sums.num_pairs(), 3);
        for measure in [PopDistMeasure::Fst, PopDistMeasure::F2] {
            let values: Vec<Option<f64>> = sums.measures(measure).collect();
            let errors: Vec<Option<f64>> = sums.standard_errors(measure).collect();
            assert_eq!(values.len(), 3);
            assert_eq!(errors.len(), 3);
            for (pair, (i, j)) in pairs.into_iter().enumerate() {
                assert_eq!(values[pair], sums.measure(measure, i, j));
                assert_eq!(errors[pair], sums.standard_error(measure, i, j));
            }
        }
        let counted: Vec<Option<u64>> = sums.num_vars_of_each_pair().collect();
        assert_eq!(counted, vec![Some(1200), Some(1200), Some(1200)]);
        let of_every_group: Vec<Option<f64>> = sums.f2_of_every_group().collect();
        assert_eq!(of_every_group.len(), 18);
        for group in 0..6 {
            for (pair, (i, j)) in pairs.into_iter().enumerate() {
                assert_eq!(
                    of_every_group[group * 3 + pair],
                    sums.f2_of_group(group, i, j),
                    "the f_2 of the pair {i} {j} within the group {group}"
                );
            }
        }
    }

    /// The measures of [`PopDistMeasure::THAT_HAVE_A_VALUE`] are the ones a
    /// pass gives a number for, which since work package 3 of
    /// `docs/plans/dists-pops.md` is all seven. Both packages refuse a
    /// measure by that array, so a measure written into [`value_of`] and
    /// not into it is refused although popnei calculates it, and one
    /// written into the array and not into `value_of` gives a user a vector
    /// of NaN read as a distance.
    #[test]
    fn the_measures_that_have_a_value_are_the_ones_a_pass_gives_a_number_for() {
        let sums = sums_of_the_worked_example(JackknifeGroups::None);

        for measure in [
            PopDistMeasure::Fst,
            PopDistMeasure::F2,
            PopDistMeasure::Chord,
            PopDistMeasure::Da,
            PopDistMeasure::Dest,
            PopDistMeasure::Gst,
            PopDistMeasure::GstStandardized,
        ] {
            assert_eq!(
                sums.measure(measure, 0, 1).is_some(),
                measure.has_a_value(),
                "the measure {}",
                measure.name()
            );
        }
        assert_eq!(
            PopDistMeasure::names_that_have_a_value(),
            [
                "fst",
                "f2",
                "chord",
                "da",
                "dest",
                "gst",
                "gst_standardized"
            ]
        );
    }

    /// The six sums of one pair within one group are 48 bytes, the five f64
    /// and the count of the variants, which is what "How it runs" of
    /// `docs/specs/dists.md` states and what the memory of a pass is
    /// counted from there: 72 KB for 3 populations and 500 groups and 29 MB
    /// for 50 populations, which is 1225 pairs.
    #[test]
    fn the_sums_of_one_pair_and_one_group_are_the_bytes_of_the_spec() {
        assert_eq!(size_of::<PairSums>(), 48);
    }

    /// The pass over a reader gives the F_ST and the f_2 of the worked
    /// example of "How it is verified" of `docs/specs/dists.md`, which the
    /// tests above assert on sums built group by group: the same numbers
    /// come out of the loop over the blocks of a reader, which cuts the
    /// five variants into the three blocks of this reader.
    ///
    /// With no groups asked for there is no standard error and the sums are
    /// one run of the pairs, the one pair here.
    #[test]
    fn the_pass_over_a_reader_gives_the_fst_and_the_f2_of_the_worked_example() {
        let mut reader = GivenBlocks::of_the_worked_example(2, true, false);
        let pops = pops_of_the_worked_example();
        let sums = calc_pop_dist_sums(
            &mut reader,
            &pops,
            &PopDistOptions {
                min_num_individuals: 1,
                groups: JackknifeGroups::None,
            },
        )
        .expect("the sums of the worked example");

        assert_eq!(sums.num_pops(), 2);
        assert_eq!(sums.num_vars(), 5);
        assert_eq!(sums.num_vars_of(0, 1), Some(5));
        assert_eq!(sums.groups(), []);
        assert_it_is_within(
            sums.measure(PopDistMeasure::Fst, 0, 1),
            0.351196,
            1e-6,
            "the F_ST of the worked example",
        );
        assert_it_is_within(
            sums.measure(PopDistMeasure::F2, 0, 1),
            0.203889,
            1e-6,
            "the f_2 of the worked example",
        );
        assert_eq!(sums.standard_error(PopDistMeasure::F2, 0, 1), None);
    }

    /// Every number of the biallelic panel a user sees, in the order they
    /// are asserted in: each of the seven measures of the three pairs, the
    /// standard error of each, the f_2 of each pair within each group, and
    /// the variants that counted for each pair.
    ///
    /// All seven are read and not F_ST and f_2 alone. Those two are built
    /// from the sums of H_b and of H_w, so the three sums the other five
    /// read, the square roots of the products and the corrected H_S and
    /// H_T, went uncompared across the sizes of the blocks and the numbers
    /// of threads while the loop held the two: adding 1e-7 to each of the
    /// three at every chunk a block was cut into, which moves the chord
    /// distance of the first pair from 0.1802670352621946 to
    /// 0.18026703664902619 between blocks of 10000 variants and blocks of
    /// 100, left both tests below green then and fails the first of them
    /// now.
    fn every_number_of(sums: &PopDistSums) -> Vec<f64> {
        let pairs = [(0, 1), (0, 2), (1, 2)];
        let mut numbers = Vec::new();
        for measure in PopDistMeasure::THAT_HAVE_A_VALUE {
            numbers.extend(
                sums.measures(measure)
                    .map(|value| value.unwrap_or(f64::NAN)),
            );
            numbers.extend(
                pairs
                    .iter()
                    .map(|(i, j)| sums.standard_error(measure, *i, *j).unwrap_or(f64::NAN)),
            );
        }
        for group in 0..sums.groups().len() {
            numbers.extend(
                pairs
                    .iter()
                    .map(|(i, j)| sums.f2_of_group(group, *i, *j).unwrap_or(f64::NAN)),
            );
        }
        numbers.extend(
            pairs
                .iter()
                .map(|(i, j)| sums.num_vars_of(*i, *j).unwrap_or(0) as f64),
        );
        numbers
    }

    /// Every number of the biallelic panel cut into the groups `how` says,
    /// read in blocks of `num_vars_per_block` variants, with a check that
    /// the variants fell into `num_groups` groups.
    fn every_number_of_the_panel_cut(
        how: JackknifeGroups,
        num_vars_per_block: usize,
        num_groups: usize,
    ) -> Vec<f64> {
        let mut reader = reader_of_the_panel("dists/panel.vcf.gz", Some(num_vars_per_block));
        let pops = pops_of_the_file(
            "stats/panel_pops.txt",
            reader.individuals(),
            &["p0", "p1", "p2"],
        );
        let sums = calc_pop_dist_sums(
            &mut reader,
            &pops,
            &PopDistOptions {
                min_num_individuals: 20,
                groups: how,
            },
        )
        .expect("the sums of the panel");

        assert_eq!(sums.groups().len(), num_groups, "the groups of the panel");
        every_number_of(&sums)
    }

    /// Every number of the biallelic panel cut into groups of 50 000 base
    /// pairs, 24 of them over its two chromosomes, read in blocks of
    /// `num_vars_per_block` variants. The groups are 24 and not the 12 of
    /// the runs of ADMIXTOOLS 2 because this goes through the function a
    /// user calls, which asks for 20 groups at least.
    fn every_number_of_the_panel(num_vars_per_block: usize) -> Vec<f64> {
        every_number_of_the_panel_cut(JackknifeGroups::OfBasePairs(50_000), num_vars_per_block, 24)
    }

    /// Two runs of the same panel that have to give the same bits, which
    /// `assert_eq!` on the numbers themselves does not say: a NaN is not
    /// equal to itself, and a standard error the pass has none of is one.
    ///
    /// The test that compares two pools of threads is the one that uses it,
    /// and that test is compiled for the targets that have rayon, which are
    /// the ones that are not wasm.
    #[cfg(not(target_family = "wasm"))]
    fn assert_they_are_the_same_bits(of_one_run: &[f64], of_the_other: &[f64], what: &str) {
        assert_eq!(
            of_one_run.len(),
            of_the_other.len(),
            "the numbers of {what}"
        );
        for (at, (of_one_run, of_the_other)) in of_one_run.iter().zip(of_the_other).enumerate() {
            assert_eq!(
                of_one_run.to_bits(),
                of_the_other.to_bits(),
                "the number {at} of {what} is {of_one_run} in one run and {of_the_other} in the other"
            );
        }
    }

    /// The size of the blocks changes no number of the panel beyond the
    /// last bits: the sums of a group are added row by row and the
    /// divisions happen once, so blocks of 100 variants and blocks of
    /// 10000, which hold the whole panel in one, agree within 1e-12
    /// relative, which is what "How it is verified" of
    /// `docs/specs/dists.md` asks of every measure, all seven of which are
    /// compared here with their standard errors. What differs is the
    /// order the rows of a group are added in, which the boundaries of the
    /// blocks move.
    #[test]
    fn the_size_of_the_blocks_does_not_change_the_measures() {
        let of_a_hundred = every_number_of_the_panel(100);
        let of_ten_thousand = every_number_of_the_panel(10_000);

        assert_eq!(of_a_hundred.len(), 117);
        assert_eq!(of_ten_thousand.len(), 117);
        for (at, (of_a_hundred, of_ten_thousand)) in
            of_a_hundred.iter().zip(&of_ten_thousand).enumerate()
        {
            let apart = (of_a_hundred - of_ten_thousand).abs();
            assert!(
                apart <= 1e-12 * of_a_hundred.abs(),
                "the number {at} is {of_a_hundred} in blocks of 100 and {of_ten_thousand} in blocks of 10000, {apart} apart"
            );
        }
    }

    /// The threads change no number of the panel at all: the rows of a
    /// block are added up in chunks of a fixed size and the chunks are
    /// added into the sums of the pass in the order of the block, so a pool
    /// of one thread and a pool of four give the same bits, which rayon's
    /// own `reduce` would not.
    ///
    /// The second run is what can fail. The first cuts the panel into
    /// groups of 50 000 base pairs, which are 50 variants where a chunk is
    /// `ROWS_PER_CHUNK` rows, so the sums of a group are built from two
    /// chunks at the most and any order of joining two values and a run of
    /// zeros gives the same bits. The second asks for no groups, where the
    /// 1200 variants of the panel come in one block of 19 chunks and all of
    /// them add into one run of the pairs. Replacing the ordered addition
    /// of the chunks in `add_the_block` with rayon's own
    /// `par_chunks(...).map(...).reduce(...)` leaves the first run
    /// unchanged and makes the second one fail: F_ST of p0 and p1 is then
    /// the bits 3fbaded19c839733 on one thread and 3fbaded19c839729 on
    /// four, 0.10496244498389444 against 0.1049624449838943.
    ///
    /// The pools are built here and are not rayon's global one, which has
    /// one thread per core of the machine. rayon is a dependency of the
    /// targets that are not wasm, so this test is compiled for those alone.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_number_of_threads_does_not_change_the_measures() {
        let in_a_pool = |threads, how, num_vars_per_block, num_groups| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            pool.install(|| every_number_of_the_panel_cut(how, num_vars_per_block, num_groups))
        };

        let groups_of_50_000 = JackknifeGroups::OfBasePairs(50_000);
        let on_one = in_a_pool(1, groups_of_50_000, 100, 24);
        assert_eq!(on_one.len(), 117);
        assert_they_are_the_same_bits(
            &on_one,
            &in_a_pool(4, groups_of_50_000, 100, 24),
            "the panel in 24 groups, read in blocks of 100 variants",
        );

        let on_one = in_a_pool(1, JackknifeGroups::None, 10_000, 0);
        assert_eq!(on_one.len(), 45);
        assert_they_are_the_same_bits(
            &on_one,
            &in_a_pool(4, JackknifeGroups::None, 10_000, 0),
            "the panel in one run of the pairs, read in one block",
        );
    }

    /// The biallelic panel with the variant at 600 000 of `chr1` moved in
    /// front of the variant at 1000, cut into groups of 5000 base pairs:
    /// the pass refuses it, and the message names the chromosome and the
    /// two positions.
    ///
    /// Before the pass refused such a source, those 1200 variants fell into
    /// 121 groups where the panel in order gives 240, one of them holding
    /// every variant of `chr1` with 600 000 as its first position and
    /// 599 000 as its last, and the standard error of f_2 for p0 and p1 was
    /// 0.0018968 where the panel in order gives 0.0017814, 6 in 100 higher.
    /// f_2 itself did not move, which is why no other test of this file
    /// showed it.
    #[test]
    fn a_pass_over_a_source_whose_variants_go_back_is_refused() {
        let mut of_the_panel = reader_of_the_panel("dists/panel.vcf.gz", Some(1200));
        of_the_panel.set_needs(Needs::ALL);
        let block = of_the_panel
            .next_block()
            .expect("a block of the panel")
            .expect("the block of the panel");
        let individuals = of_the_panel.individuals().to_vec();
        let pops = pops_of_the_file("stats/panel_pops.txt", &individuals, &["p0", "p1", "p2"]);
        let mut reader = GivenBlocks::of(
            individuals,
            vec![the_block_with_a_variant_moved_first(&block, 599)],
            false,
        );

        let error = calc_pop_dist_sums(
            &mut reader,
            &pops,
            &PopDistOptions {
                min_num_individuals: 20,
                groups: JackknifeGroups::OfBasePairs(5_000),
            },
        )
        .expect_err("a pass over a source whose variants go back");

        assert!(
            matches!(&error, Error::JackknifeGroupsVariantGoesBack { chrom, pos, before }
                if chrom == "chr1" && *pos == 1000 && *before == 600_000),
            "{error:?}"
        );
    }

    /// `block` with the variant at the row `moved` in front of its first
    /// one, the genotypes, the chromosomes and the positions of every row
    /// moved together, which is what a source whose variants are not sorted
    /// gives.
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "the rows of one block of the panel, 1200 of 200 individuals"
    )]
    fn the_block_with_a_variant_moved_first(block: &Block, moved: usize) -> Block {
        let width = block.num_individuals * block.ploidy;
        let chrom = block.chrom.clone().expect("the chromosomes of the panel");
        let pos = block.pos.clone().expect("the positions of the panel");
        let mut order: Vec<usize> = (0..block.num_vars).collect();
        let row = order.remove(moved);
        order.insert(0, row);
        Block {
            num_vars: block.num_vars,
            num_individuals: block.num_individuals,
            ploidy: block.ploidy,
            gts: order
                .iter()
                .flat_map(|row| block.gts[row * width..(row + 1) * width].to_vec())
                .collect(),
            chrom: Some(order.iter().map(|row| chrom[*row]).collect()),
            pos: Some(order.iter().map(|row| pos[*row]).collect()),
            id: None,
            alleles: None,
            qual: None,
        }
    }

    /// A pass over a source that holds no variant is refused, and the
    /// message says that the source holds none: a distance over no variant
    /// says nothing about two populations, and the user has to know whether
    /// the source was empty or the steps of the pass kept nothing.
    #[test]
    fn a_pass_over_a_source_with_no_variant_is_refused() {
        let mut reader = GivenBlocks::of(
            (0..6).map(|number| format!("i{number}")).collect(),
            Vec::new(),
            false,
        );
        let error = calc_pop_dist_sums(
            &mut reader,
            &pops_of_the_worked_example(),
            &PopDistOptions {
                min_num_individuals: 1,
                groups: JackknifeGroups::None,
            },
        )
        .expect_err("a pass with no variant");

        assert!(
            matches!(&error, Error::PassGaveNoVariant { num_vars_of_the_source, filters }
                if *num_vars_of_the_source == 0 && filters.is_empty()),
            "{error:?}"
        );
        assert!(
            error.to_string().contains("its source holds none"),
            "{error}"
        );
    }

    /// Every measure is of a pair of populations, so one population makes
    /// no pair and the pass is refused before the source is read. `Pops`
    /// refuses a `pops` that names none, so 1 is the number a user reaches.
    #[test]
    fn a_pass_of_fewer_than_two_pops_is_refused() {
        let mut reader = GivenBlocks::of_the_worked_example(2, true, false);
        let error = calc_pop_dist_sums(
            &mut reader,
            &Pops::all(6),
            &PopDistOptions {
                min_num_individuals: 1,
                groups: JackknifeGroups::None,
            },
        )
        .expect_err("a pass of one population");

        assert!(
            matches!(&error, Error::PopDistsOfFewerThanTwoPops { num_pops } if *num_pops == 1),
            "{error:?}"
        );
    }

    /// A standard error is built by leaving each group out in turn, so the
    /// pass is refused when the variants fall into fewer than 20 groups,
    /// with how many they fell into: the four variants of the worked
    /// example, each its own group, are 4. The measures themselves are
    /// there whatever the groups, which is why the number is in the
    /// message: a user who chose a length too long for their data cuts
    /// shorter groups or asks for no standard error.
    #[test]
    fn a_pass_of_fewer_than_twenty_groups_is_refused_with_their_number() {
        let mut reader = GivenBlocks::of_the_worked_example(2, true, false);
        let error = calc_pop_dist_sums(
            &mut reader,
            &pops_of_the_worked_example(),
            &PopDistOptions {
                min_num_individuals: 1,
                groups: JackknifeGroups::PerVariant,
            },
        )
        .expect_err("a pass of five groups");

        assert!(
            matches!(&error, Error::TooFewJackknifeGroups { num_groups, at_least }
                if *num_groups == 5 && *at_least == MIN_NUM_JACKKNIFE_GROUPS),
            "{error:?}"
        );
        let message = error.to_string();
        assert!(message.contains("into 5 resampling groups"), "{message}");
        assert!(message.contains("20 at least"), "{message}");
    }

    /// An error of the reader is the error of the pass, as it is: a file
    /// that was cut short is read by the user as a file that was cut short
    /// and not as a distance over the variants that were there.
    #[test]
    fn an_error_of_the_reader_is_the_error_of_the_pass() {
        let mut reader = GivenBlocks::of_the_worked_example(2, true, true);
        let error = calc_pop_dist_sums(
            &mut reader,
            &pops_of_the_worked_example(),
            &PopDistOptions {
                min_num_individuals: 1,
                groups: JackknifeGroups::None,
            },
        )
        .expect_err("the error of the reader");

        assert!(matches!(&error, Error::VcfBgzipEndMissing), "{error:?}");
    }

    /// The pass asks its reader for the genotypes alone when no standard
    /// errors were asked for, and for the genotypes with the chromosome and
    /// the position when the groups are cut from them, so that a reader
    /// over a file leaves the columns of a variant unparsed where nothing
    /// reads them.
    #[test]
    fn the_pass_asks_for_the_positions_only_where_the_groups_need_them() {
        for (how, expected) in [
            (JackknifeGroups::None, Needs::GTS),
            (JackknifeGroups::PerVariant, Needs::GTS | Needs::CHROM_POS),
            (
                JackknifeGroups::OfBasePairs(1000),
                Needs::GTS | Needs::CHROM_POS,
            ),
        ] {
            let mut reader = GivenBlocks::of_the_worked_example(4, true, false);
            let options = PopDistOptions {
                min_num_individuals: 1,
                groups: how,
            };
            sums_of_the_pass(&mut reader, &pops_of_the_worked_example(), &options)
                .expect("the sums of the worked example");

            assert_eq!(reader.needs(), expected, "{how:?}");
        }
    }

    /// A block with no position where the groups are cut from the positions
    /// is a reader that was asked for them and gave none, which is a defect
    /// of that reader: the pass says which field is missing instead of
    /// cutting every variant into one group.
    #[test]
    fn a_block_with_no_position_where_the_groups_need_one_is_refused() {
        let mut reader = GivenBlocks::of_the_worked_example(4, false, false);
        let error = sums_of_the_pass(
            &mut reader,
            &pops_of_the_worked_example(),
            &PopDistOptions {
                min_num_individuals: 1,
                groups: JackknifeGroups::PerVariant,
            },
        )
        .expect_err("a block with no position");

        assert!(
            matches!(&error, Error::FieldsNotInTheBlock { fields }
                if *fields == Needs::CHROM_POS),
            "{error:?}"
        );
    }
}

#[cfg(test)]
mod pop_dist_measure {
    use super::PopDistMeasure;
    use crate::error::Error;

    /// The name of each measure is what a user writes in `measures` and the
    /// field of the result that holds it, and `of_name` gives back the
    /// measure of each name. The names live in `NAMES` alone, so a rename is
    /// one change; the literals here are the names of the fields of
    /// `PopDists` in Python and in TypeScript.
    #[test]
    fn each_name_is_the_name_of_the_measure_it_gives_back() {
        assert_eq!(
            PopDistMeasure::NAMES,
            [
                "fst",
                "f2",
                "chord",
                "da",
                "dest",
                "gst",
                "gst_standardized"
            ]
        );
        for (name, measure) in [
            ("fst", PopDistMeasure::Fst),
            ("f2", PopDistMeasure::F2),
            ("chord", PopDistMeasure::Chord),
            ("da", PopDistMeasure::Da),
            ("dest", PopDistMeasure::Dest),
            ("gst", PopDistMeasure::Gst),
            ("gst_standardized", PopDistMeasure::GstStandardized),
        ] {
            assert_eq!(
                PopDistMeasure::of_name(name).unwrap_or_else(|error| panic!("{name}: {error}")),
                measure
            );
            assert_eq!(measure.name(), name);
        }
    }

    /// A name that is of none of the seven is refused, with the seven names:
    /// a user who writes one of them wrong has to read which they are.
    #[test]
    fn a_name_of_no_measure_is_refused_with_the_seven() {
        for name in ["fsts", "FST", "", "f2 "] {
            let error = PopDistMeasure::of_name(name).unwrap_err();
            assert!(
                matches!(&error, Error::PopDistMeasureOfAnUnknownName { name: found }
                    if found == name),
                "{name}: {error:?}"
            );
            let message = error.to_string();
            for of_the_seven in PopDistMeasure::NAMES {
                assert!(message.contains(of_the_seven), "{message}");
            }
        }
    }
}
