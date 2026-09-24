//! How much variety each population of a dataset holds.
//!
//! One pass over the variants gives, for each population, the alleles it
//! called and the private ones among them, the variants that vary in it,
//! the folded site frequency spectrum and F_IS, the first three also taken
//! down to a common number of called alleles so that a population of 20
//! individuals and one of 200 can be compared.
//! `docs/specs/diversity.md` says what each of the five is, where its
//! numbers were checked and which program gave them.
//!
//! All five are built from one thing, which the pass works out once for
//! each population and each variant: how often the population called each
//! allele there, and the sum of those counts, the alleles it called. A
//! variant counts for a population when the population called something at
//! it and called at least `min_num_individuals` genotypes, measured as
//! those alleles over the ploidy, so a half called genotype counts as half
//! an individual.
//!
//! F_IS is the one of the five that reads whole genotypes and not only
//! the alleles called: it is one minus the mean observed heterozygosity of
//! the population over the mean unbiased expected one, and both of those
//! are the per variant statistics of `docs/specs/stats.md`, computed here
//! by the `stats` module itself so that the two specs cannot disagree
//! about a heterozygosity.
//!
//! What is built here so far is the pass, those counts of variants, the
//! alleles each population called, the private ones among them, the
//! variants that vary in it and F_IS. Four things are not, and all four
//! are the draw of a common number of called alleles, which is work
//! package 3 of `docs/plans/diversity.md`: `num_called_alleles` is taken
//! and checked and nothing reads it, so no variant is counted as being in
//! a draw; the three standardized values, the mean alleles, the mean
//! private alleles and the ratio of variable variants of a draw, have no
//! value; and neither has the folded site frequency spectrum. They are
//! built on top of the same counts.

use std::collections::HashSet;

use crate::block::{Block, BlockReader, alleles_of_a_chunk, alleles_per_var_of};
use crate::error::{Error, Result};
use crate::stats::{ExpHet, ObsHet, checked_ploidy, every_individual_in_order, min_called_alleles};
use crate::variant::{
    AlleleCounts, GtCounts, Needs, count_alleles, count_alleles_and_gts_of, count_alleles_of,
    count_gts,
};

/// Which of the five statistics of this module a pass is asked for: a set
/// of them, with union, [`DiversityStats::contains`] and
/// [`DiversityStats::empty`].
///
/// A statistic that is not in the set is not counted and has no value in
/// the result, as the `stats` argument of a Python or a TypeScript user
/// asks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DiversityStats(u8);

impl DiversityStats {
    /// The alleles the population called, as a total and as a mean.
    pub const NUM_ALLELES: DiversityStats = DiversityStats(1);
    /// Of those, the ones no other population of the call called.
    pub const PRIVATE_ALLELES: DiversityStats = DiversityStats(2);
    /// The variants at which the population called more than one allele,
    /// as a total and as a ratio.
    pub const VARIABLE_VARS_RATIO: DiversityStats = DiversityStats(4);
    /// The folded site frequency spectrum, which needs a draw size.
    pub const FOLDED_SFS: DiversityStats = DiversityStats(8);
    /// The inbreeding coefficient F_IS.
    pub const FIS: DiversityStats = DiversityStats(16);
    /// The five statistics, built from the five constants, so that one
    /// added to this set later cannot be left out of it.
    pub const ALL: DiversityStats = DiversityStats::NUM_ALLELES
        .union(DiversityStats::PRIVATE_ALLELES)
        .union(DiversityStats::VARIABLE_VARS_RATIO)
        .union(DiversityStats::FOLDED_SFS)
        .union(DiversityStats::FIS);

    /// The name of each of the five statistics beside the statistic it
    /// names, in the order of the constants above, which is the order of
    /// the fields of a result.
    ///
    /// The names are what a Python and a TypeScript user writes in `stats`,
    /// and each one is the field of the result that holds that statistic.
    /// They are here and not in the binding crates so that a rename is one
    /// change and not three. Each name stands beside its statistic, and
    /// [`DiversityStats::of_name`] reads this table and nothing else, so a
    /// statistic added to it is understood with nothing else written for
    /// it: a name added to a table of names alone, with the statistics
    /// matched on elsewhere, would be listed to a user as one they may
    /// write and then refused when they wrote it.
    pub const NAMES_AND_STATS: [(&'static str, DiversityStats); 5] = [
        ("num_alleles", DiversityStats::NUM_ALLELES),
        ("private_alleles", DiversityStats::PRIVATE_ALLELES),
        ("variable_vars_ratio", DiversityStats::VARIABLE_VARS_RATIO),
        ("folded_sfs", DiversityStats::FOLDED_SFS),
        ("fis", DiversityStats::FIS),
    ];

    /// The name of each of the five statistics, in the order of
    /// [`DiversityStats::NAMES_AND_STATS`], which a message that asks a
    /// user to choose among them lists.
    pub const NAMES: [&'static str; 5] = DiversityStats::the_names();

    /// The names of [`DiversityStats::NAMES_AND_STATS`] on their own,
    /// which is what [`DiversityStats::NAMES`] holds.
    const fn the_names() -> [&'static str; 5] {
        let [
            (num_alleles, _),
            (private_alleles, _),
            (variable_vars_ratio, _),
            (folded_sfs, _),
            (fis, _),
        ] = DiversityStats::NAMES_AND_STATS;
        [
            num_alleles,
            private_alleles,
            variable_vars_ratio,
            folded_sfs,
            fis,
        ]
    }

    /// The one statistic a user named.
    ///
    /// A binding crate reads the names a user wrote with this and gives the
    /// error on as it is, as the binding crates of the per variant
    /// statistics do with `stats::PerVarStat::of_name`: each of them turns
    /// an error of the core into the exception of its language in one
    /// place, so neither has to write the sentence a user reads.
    ///
    /// # Errors
    ///
    /// A name that is of none of the five, with the five names.
    pub fn of_name(name: &str) -> Result<DiversityStats> {
        // The table is walked rather than matched on, so that a statistic
        // added to it needs no arm here and cannot be listed to a user as a
        // name they may write and then refused when they write it.
        DiversityStats::NAMES_AND_STATS
            .iter()
            .find(|(known, _)| *known == name)
            .map(|(_, stat)| *stat)
            .ok_or_else(|| Error::DiversityStatOfAnUnknownName {
                name: name.to_owned(),
            })
    }

    /// No statistic at all, which a caller that builds a set one name at a
    /// time starts from.
    #[must_use]
    pub const fn empty() -> DiversityStats {
        DiversityStats(0)
    }

    /// Whether every statistic of `stats` is in this set. An empty `stats`
    /// is in every set.
    #[must_use]
    pub const fn contains(self, stats: DiversityStats) -> bool {
        self.0 & stats.0 == stats.0
    }

    /// The statistics of either set.
    #[must_use]
    pub const fn union(self, other: DiversityStats) -> DiversityStats {
        DiversityStats(self.0 | other.0)
    }
}

impl std::ops::BitOr for DiversityStats {
    type Output = DiversityStats;

    fn bitor(self, other: DiversityStats) -> DiversityStats {
        self.union(other)
    }
}

impl std::ops::BitOrAssign for DiversityStats {
    fn bitor_assign(&mut self, other: DiversityStats) {
        *self = self.union(other);
    }
}

/// What a pass of [`calc_pop_diversity`] is asked for.
#[derive(Debug, Clone)]
pub struct DiversityOptions {
    /// The statistics to compute. One that is not in the set is not
    /// counted and has no value in the result.
    pub stats: DiversityStats,
    /// The called alleles every population is brought down to. `None` gives
    /// no standardized value and no spectrum.
    pub num_called_alleles: Option<u32>,
    /// The called genotypes a population needs at a variant, its called
    /// alleles over the ploidy.
    pub min_num_individuals: u32,
}

/// What one population of a pass has counted, which grows neither with the
/// variants nor with the individuals.
#[derive(Debug, Clone, Copy)]
struct OfAPop {
    /// The variants the population called something at and had
    /// `min_num_individuals` called genotypes in.
    num_vars: u64,
    /// The alleles it called at those variants, added up. An allele
    /// numbered 1 at one variant is not the one numbered 1 at the next, so
    /// this is a sum of per variant counts and never a count of distinct
    /// things across the dataset.
    num_alleles: u64,
    /// Of those alleles, the ones no other population of the pass called at
    /// the variant, added up over the variants that counted for every
    /// population. Its divisor is `num_vars_every_pop` and not `num_vars`.
    private_alleles: u64,
    /// How many of those variants it called more than one allele at.
    num_variable_vars: u64,
    /// The observed heterozygosities of the variants behind its F_IS,
    /// added up.
    sum_obs_het: f64,
    /// The unbiased expected heterozygosities of the same variants, added
    /// up.
    sum_unbiased_exp_het: f64,
    /// How many variants those two sums are over: the ones that counted for
    /// the population and at which both heterozygosities exist. It is the
    /// divisor of both means, so F_IS does not depend on it.
    num_vars_with_both_hets: u64,
}

impl OfAPop {
    /// The counts of one population before any variant is read.
    fn none() -> OfAPop {
        OfAPop {
            num_vars: 0,
            num_alleles: 0,
            private_alleles: 0,
            num_variable_vars: 0,
            sum_obs_het: 0.0,
            sum_unbiased_exp_het: 0.0,
            num_vars_with_both_hets: 0,
        }
    }

    /// One minus the mean observed heterozygosity of the population over
    /// its mean unbiased expected one, and NaN when it has no F_IS: when no
    /// variant of it carries both heterozygosities, and when its mean
    /// unbiased expected heterozygosity is 0, every variant it counted
    /// having held one allele.
    fn fis(&self) -> f64 {
        if self.num_vars_with_both_hets == 0 {
            return f64::NAN;
        }
        // A count below 2^53 is exact in a float64, and a pass of that many
        // variants reads more rows than any source holds.
        let num_vars = self.num_vars_with_both_hets as f64;
        let mean_obs_het = self.sum_obs_het / num_vars;
        let mean_unbiased_exp_het = self.sum_unbiased_exp_het / num_vars;
        if mean_unbiased_exp_het == 0.0 {
            return f64::NAN;
        }
        1.0 - mean_obs_het / mean_unbiased_exp_het
    }
}

/// What one pass of [`calc_pop_diversity`] gives back, per population.
///
/// The populations are in the order they were given to the pass, which is
/// the order the keys of a user's `pops` iterate in. A statistic that was
/// not asked for has no value here.
#[derive(Debug)]
pub struct PopDiversity {
    /// What each population counted, in the order of the call.
    pops: Vec<OfAPop>,
    /// The variants that counted for every population at once.
    num_vars_every_pop: u64,
    /// The variants the reader gave, whether any population counted them or
    /// not.
    num_vars_of_the_pass: u64,
    /// The statistics the pass was asked for, which are the ones that have
    /// a value here.
    stats: DiversityStats,
}

impl PopDiversity {
    /// How many populations the pass was over, which is 1 when the caller
    /// gave none.
    #[must_use]
    pub fn num_pops(&self) -> usize {
        self.pops.len()
    }

    /// The variants the population called something at and had
    /// `min_num_individuals` called genotypes in. `None` when `pop` is not
    /// a population of the call.
    #[must_use]
    pub fn num_vars(&self, pop: usize) -> Option<u64> {
        self.pops.get(pop).map(|pop| pop.num_vars)
    }

    /// The variants the reader gave the pass, whether any population counted
    /// them or not.
    ///
    /// It is what a caller reports as the variants of the pass beside what
    /// each filter of the chain was given and kept, and it is not the
    /// divisor of anything here: a variant a population has too little
    /// called at is one of these and is out of every count of that
    /// population, which [`PopDiversity::num_vars`] gives.
    #[must_use]
    pub fn num_vars_of_the_pass(&self) -> u64 {
        self.num_vars_of_the_pass
    }

    /// The variants that counted for every population. It is the divisor of
    /// the private alleles.
    #[must_use]
    pub fn num_vars_every_pop(&self) -> u64 {
        self.num_vars_every_pop
    }

    /// The alleles the population called, summed over the variants that
    /// counted for it. `None` when `pop` is not a population of the call or
    /// [`DiversityStats::NUM_ALLELES`] was not asked for.
    ///
    /// Their mean, the allelic richness a user compares between
    /// populations, is this over [`PopDiversity::num_vars`].
    #[must_use]
    pub fn num_alleles(&self, pop: usize) -> Option<u64> {
        if !self.stats.contains(DiversityStats::NUM_ALLELES) {
            return None;
        }
        self.pops.get(pop).map(|pop| pop.num_alleles)
    }

    /// The alleles the population called that no other population of the
    /// call called at the same variant, summed over the variants that
    /// counted for every population. `None` when `pop` is not a population
    /// of the call or [`DiversityStats::PRIVATE_ALLELES`] was not asked
    /// for.
    ///
    /// Their mean is this over [`PopDiversity::num_vars_every_pop`] and not
    /// over [`PopDiversity::num_vars`]: a variant where one population has
    /// too little called is out of the private alleles of every population,
    /// since there the alleles of the others would all be private and the
    /// count would measure the missing data.
    ///
    /// With one population every allele it called is private, since there
    /// is no other population to hold it, so this is then
    /// [`PopDiversity::num_alleles`].
    #[must_use]
    pub fn private_alleles(&self, pop: usize) -> Option<u64> {
        if !self.stats.contains(DiversityStats::PRIVATE_ALLELES) {
            return None;
        }
        self.pops.get(pop).map(|pop| pop.private_alleles)
    }

    /// The variants where the population called more than one allele.
    /// `None` when `pop` is not a population of the call or
    /// [`DiversityStats::VARIABLE_VARS_RATIO`] was not asked for.
    ///
    /// Their ratio is this over [`PopDiversity::num_vars`].
    #[must_use]
    pub fn num_variable_vars(&self, pop: usize) -> Option<u64> {
        if !self.stats.contains(DiversityStats::VARIABLE_VARS_RATIO) {
            return None;
        }
        self.pops.get(pop).map(|pop| pop.num_variable_vars)
    }

    /// How far the genotypes of the population are from the proportions its
    /// allele frequencies would give if its individuals paired at random:
    /// one minus its mean observed heterozygosity over its mean unbiased
    /// expected one, 0 when the two match, positive when it holds fewer
    /// heterozygous genotypes than random pairing would give and negative
    /// when it holds more. `None` when `pop` is not a population of the
    /// call or [`DiversityStats::FIS`] was not asked for.
    ///
    /// Both means are over the same variants, those that counted for the
    /// population and at which both heterozygosities exist: the observed
    /// one needs a called genotype, since it divides by the called ones,
    /// and the unbiased expected one needs as many called alleles as a
    /// genotype holds, since it draws that many of them without
    /// replacement. It is Nei's F_IS, read on one population on its own,
    /// and not Weir and Cockerham's, which comes out of a decomposition of
    /// the variance across populations.
    ///
    /// It is NaN when the population has no F_IS, in four cases. When no
    /// variant counted for it, which is not an error. When no variant that
    /// counted for it carries both heterozygosities, which a population
    /// whose counted variants hold no whole called genotype reaches: two
    /// individuals whose genotypes are all half called count their
    /// variants and have no observed heterozygosity at any of them. When
    /// its mean unbiased expected heterozygosity is 0, every variant it
    /// counted having held one allele. And at ploidy 1, where no genotype
    /// can be heterozygous, so the observed heterozygosity is 0 at every
    /// variant and the ratio would be 1 wherever the population has any
    /// diversity. The draw of
    /// `num_called_alleles` does not touch it: the observed heterozygosity
    /// is a property of whole genotypes and not of a sample of alleles.
    #[must_use]
    pub fn fis(&self, pop: usize) -> Option<f64> {
        if !self.stats.contains(DiversityStats::FIS) {
            return None;
        }
        self.pops.get(pop).map(OfAPop::fis)
    }

    /// The sum of the observed heterozygosities of one population, the sum
    /// of its unbiased expected ones and how many variants the two are
    /// over.
    ///
    /// A test that asks whether the parts of a sum of float64 were added in
    /// the order of the variants compares these bit for bit, and not
    /// [`PopDiversity::fis`]: the division that makes F_IS out of the two
    /// sums absorbs a difference of their last bits, so the 200 variants of
    /// such a test give one F_IS from two sums that differ by one unit of
    /// the last place.
    #[cfg(test)]
    fn the_sums_behind_the_fis(&self, pop: usize) -> Option<(f64, f64, u64)> {
        self.pops.get(pop).map(|pop| {
            (
                pop.sum_obs_het,
                pop.sum_unbiased_exp_het,
                pop.num_vars_with_both_hets,
            )
        })
    }
}

/// One population of a pass: the individuals it holds among those of the
/// reader, and whether they are every individual of the reader in its
/// order, which lets a row be counted as it is.
#[derive(Debug)]
struct PopOfThePass {
    individuals: Vec<usize>,
    is_all: bool,
}

/// Whether a pass counts the alleles a population called at a variant.
///
/// The alleles called and the variable variants are the two statistics
/// that read them, and both come from one walk over the counts of the
/// variant, so a pass asked for either takes both and gives back the one
/// it was asked for. A pass asked for neither walks nothing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CountsTheAlleles {
    Yes,
    No,
}

impl CountsTheAlleles {
    /// Whether `stats` holds a statistic that reads the alleles a
    /// population called at a variant.
    fn of(stats: DiversityStats) -> CountsTheAlleles {
        if stats.contains(DiversityStats::NUM_ALLELES)
            || stats.contains(DiversityStats::VARIABLE_VARS_RATIO)
        {
            CountsTheAlleles::Yes
        } else {
            CountsTheAlleles::No
        }
    }
}

/// Whether a pass counts, at a variant that counted for every population,
/// the alleles of each population that no other population called there.
///
/// It is the one statistic of this module that reads more than one
/// population at a time, so it is asked for on its own and not with
/// [`CountsTheAlleles`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CountsThePrivateAlleles {
    Yes,
    No,
}

impl CountsThePrivateAlleles {
    /// Whether `stats` holds the private alleles.
    fn of(stats: DiversityStats) -> CountsThePrivateAlleles {
        if stats.contains(DiversityStats::PRIVATE_ALLELES) {
            CountsThePrivateAlleles::Yes
        } else {
            CountsThePrivateAlleles::No
        }
    }
}

/// The two per variant statistics of `docs/specs/stats.md` that F_IS is
/// built from, which a pass asked for F_IS carries and one asked for the
/// other four does not.
///
/// They are the ones of the `stats` module and not a second copy of either
/// formula, so that the two specs cannot disagree about a heterozygosity.
/// Each is built with a threshold of no individual, because the pass has
/// already decided which variants count for a population, by its own rule
/// in called alleles; what is left for the two of them to say is whether
/// each exists at the variant, which is where they give `None`.
#[derive(Debug, Clone, Copy)]
struct Heterozygosities {
    /// The heterozygous genotypes of the population over its called ones,
    /// which is `None` where it called no whole genotype.
    obs_het: ObsHet,
    /// The chance that as many gene copies as a genotype holds, drawn from
    /// the called alleles of the population without replacement, are not
    /// all of the same allele. It is `None` where the population called
    /// fewer alleles than a genotype holds, which for a diploid population
    /// is one called allele.
    unbiased_exp_het: ExpHet,
}

impl Heterozygosities {
    /// The two statistics of a pass asked for F_IS over variants of the
    /// ploidy `ploidy`, and `None` for a pass that gives no F_IS: one that
    /// did not ask for it, and one over haploid variants, where no genotype
    /// can be heterozygous and the ratio would be 1 wherever the population
    /// has any diversity.
    ///
    /// The exponent of the unbiased expected heterozygosity is the ploidy,
    /// which is what `docs/specs/stats.md` gives a caller that asks for no
    /// other.
    ///
    /// # Errors
    ///
    /// A ploidy above the largest one the statistics of
    /// `docs/specs/stats.md` take, 255, which no reader of popnei gives.
    fn of(stats: DiversityStats, ploidy: usize) -> Result<Option<Heterozygosities>> {
        if !stats.contains(DiversityStats::FIS) || ploidy < 2 {
            return Ok(None);
        }
        Ok(Some(Heterozygosities {
            obs_het: ObsHet::new(0),
            unbiased_exp_het: ExpHet::new(ploidy, ploidy, 0)?,
        }))
    }

    /// It adds the two heterozygosities of one population at one variant to
    /// what the population has counted, and adds nothing when either of
    /// them does not exist there: both means of F_IS are over the variants
    /// that carry both, so the two sums have one divisor.
    ///
    /// `counts` is how often the population called each allele of the
    /// variant, `called_alleles` their sum and `gts` how many of its
    /// genotypes there were called, missing and heterozygous.
    fn add_the_var(
        &self,
        counts: &AlleleCounts,
        called_alleles: u32,
        gts: GtCounts,
        counted: &mut OfAPop,
    ) {
        let (Some(obs_het), Some(unbiased_exp_het)) = (
            self.obs_het.of_var(gts),
            self.unbiased_exp_het.of_var(counts, called_alleles, true),
        ) else {
            return;
        };
        counted.sum_obs_het += obs_het;
        counted.sum_unbiased_exp_het += unbiased_exp_het;
        // One variant of the pass, and a pass of more than
        // 18446744073709551615 variants reads more rows than any source
        // holds.
        counted.num_vars_with_both_hets = counted.num_vars_with_both_hets.saturating_add(1);
    }
}

/// What every row of a pass is read with: the populations and the rule for
/// which variants count for them.
#[derive(Debug)]
struct OfThePass<'a> {
    pops: &'a [PopOfThePass],
    /// Whether the alleles a population called at a variant are counted at
    /// all, which the alleles called and the variable variants both need
    /// and a pass asked for neither of them does without.
    counts_the_alleles: CountsTheAlleles,
    /// Whether the alleles no other population called are counted, which
    /// walks the counts of every population of the row a second time.
    counts_the_private_alleles: CountsThePrivateAlleles,
    /// The two heterozygosities F_IS is built from, which a pass that gives
    /// no F_IS does not carry and whose genotypes it does not count.
    heterozygosities: Option<Heterozygosities>,
    /// How many alleles a population has to have called at a variant for
    /// the variant to count for it: `min_num_individuals` genotypes of the
    /// ploidy.
    min_called_alleles: u64,
    /// The alleles of one genotype, which the reader says its source has.
    ploidy: usize,
}

/// What a pass, or one chunk of the rows of a block, has counted over every
/// population.
#[derive(Debug, Clone)]
struct Totals {
    pops: Vec<OfAPop>,
    num_vars_every_pop: u64,
}

impl Totals {
    /// The counts of `num_pops` populations before any variant is read.
    fn of(num_pops: usize) -> Totals {
        Totals {
            pops: vec![OfAPop::none(); num_pops],
            num_vars_every_pop: 0,
        }
    }

    /// It adds what one chunk of rows counted to what the pass has.
    fn add_the_chunk(&mut self, of_the_chunk: &Totals) {
        for (of_the_pass, of_the_chunk) in self.pops.iter_mut().zip(&of_the_chunk.pops) {
            of_the_pass.num_vars = of_the_pass.num_vars.saturating_add(of_the_chunk.num_vars);
            of_the_pass.num_alleles = of_the_pass
                .num_alleles
                .saturating_add(of_the_chunk.num_alleles);
            of_the_pass.private_alleles = of_the_pass
                .private_alleles
                .saturating_add(of_the_chunk.private_alleles);
            of_the_pass.num_variable_vars = of_the_pass
                .num_variable_vars
                .saturating_add(of_the_chunk.num_variable_vars);
            // The chunks of a block are added in the order of the block and
            // the blocks in the order of the pass, so these two sums of
            // float64 do not depend on how many threads read the rows.
            of_the_pass.sum_obs_het += of_the_chunk.sum_obs_het;
            of_the_pass.sum_unbiased_exp_het += of_the_chunk.sum_unbiased_exp_het;
            of_the_pass.num_vars_with_both_hets = of_the_pass
                .num_vars_with_both_hets
                .saturating_add(of_the_chunk.num_vars_with_both_hets);
        }
        self.num_vars_every_pop = self
            .num_vars_every_pop
            .saturating_add(of_the_chunk.num_vars_every_pop);
    }

    /// It empties every count, so that one chunk of rows after another is
    /// counted in the same `Totals`.
    fn forget_what_it_holds(&mut self) {
        for of_the_pop in &mut self.pops {
            *of_the_pop = OfAPop::none();
        }
        self.num_vars_every_pop = 0;
    }
}

/// The diversity of every population over the variants `reader` gives,
/// which is one pass over the source through the steps the variants carry.
///
/// `reader` is the outermost reader of the chain of the pass, lent and not
/// taken, so that whoever built the chain reads the counts of its filters
/// from it when this returns. The pass asks the reader for the genotypes
/// alone.
///
/// `pops` is the indices of the individuals of each population, in the
/// order the user gave them, and an empty slice is one population of every
/// individual of the reader. An individual may be in more than one
/// population and is in each of them once, as `docs/specs/stats.md` has it.
///
/// # Errors
///
/// A `stats` that holds no statistic, [`DiversityStats::FOLDED_SFS`] asked
/// for with no `num_called_alleles`, a `num_called_alleles` below 2, a
/// population with no individual, an index that is not an individual of the
/// dataset, an individual asked for more than once, no variant in the
/// reader, a variant of more alleles than a count of them holds, and those
/// of the reader, among them a ploidy of 0 or above the 255 a genotype of
/// popnei holds.
pub fn calc_pop_diversity<R: BlockReader + ?Sized>(
    reader: &mut R,
    pops: &[&[usize]],
    options: &DiversityOptions,
) -> Result<PopDiversity> {
    the_pass(reader, pops, options, add_the_block)
}

/// The same pass with the chunks of every block read one after another,
/// which is the reduction WebAssembly runs, it having no threads.
///
/// A native build reaches [`add_the_chunks_one_by_one`] only where a block
/// is refused and its rows are read again to find the first that is an
/// error, so without this no cargo test counts a variant through it. A test
/// runs it beside [`calc_pop_diversity`] and reads the same numbers.
///
/// # Errors
///
/// Those of [`calc_pop_diversity`].
#[cfg(test)]
fn calc_pop_diversity_one_chunk_at_a_time<R: BlockReader + ?Sized>(
    reader: &mut R,
    pops: &[&[usize]],
    options: &DiversityOptions,
) -> Result<PopDiversity> {
    the_pass(reader, pops, options, add_the_chunks_one_by_one)
}

/// The pass of [`calc_pop_diversity`], with `add_the_block` the way the
/// chunks of one block are read: on the threads of rayon, or one after
/// another as WebAssembly reads them.
///
/// # Errors
///
/// Those of [`calc_pop_diversity`].
fn the_pass<R: BlockReader + ?Sized>(
    reader: &mut R,
    pops: &[&[usize]],
    options: &DiversityOptions,
    add_the_block: fn(&Block, usize, &OfThePass, &mut Totals) -> Result<()>,
) -> Result<PopDiversity> {
    check_the_statistics(options)?;
    check_the_draw(options)?;
    let num_individuals = reader.individuals().len();
    let of_the_pops = pops_of_the_pass(pops, num_individuals)?;
    let ploidy = reader.ploidy();
    // The ploidy the reader states, which the threshold of called
    // genotypes is measured in. A reader that states one of 0, or one above
    // the 255 a genotype of popnei holds, is refused here: a threshold
    // built from a ploidy that did not fit would be one no population ever
    // meets, and every population would silently count no variant.
    let ploidy_of_the_gts = checked_ploidy("ploidy", ploidy)?;
    // The five statistics follow from the genotypes of a row, so no column
    // of a block is read and the reader is asked to fill none of them.
    reader.set_needs(Needs::GTS);
    let of_the_pass = OfThePass {
        pops: &of_the_pops,
        counts_the_alleles: CountsTheAlleles::of(options.stats),
        counts_the_private_alleles: CountsThePrivateAlleles::of(options.stats),
        heterozygosities: Heterozygosities::of(options.stats, ploidy)?,
        min_called_alleles: min_called_alleles(options.min_num_individuals, ploidy_of_the_gts),
        ploidy,
    };
    let mut totals = Totals::of(of_the_pops.len());
    let mut num_vars: u64 = 0;
    while let Some(block) = reader.next_block()? {
        let alleles_per_var = alleles_per_var_of(&block, num_individuals, ploidy)?;
        add_the_block(&block, alleles_per_var, &of_the_pass, &mut totals)?;
        // A `usize` is 64 bits on the targets popnei builds natively for and
        // 32 in wasm, so every one of them is a `u64`; a block that said it
        // held more is refused rather than counted into a number that
        // stopped at the largest one, which would leave the pass reporting
        // fewer variants than it read. A pass of more than
        // 18446744073709551615 variants reads more rows than any source
        // holds.
        let of_the_block =
            u64::try_from(block.num_vars).map_err(|_| Error::DiversityMoreVarsThanACountHolds {
                num_vars: block.num_vars,
            })?;
        num_vars = num_vars.saturating_add(of_the_block);
    }
    if num_vars == 0 {
        let filters = reader.filtering_stats();
        return Err(Error::PassGaveNoVariant {
            // The filter nearest the source was given what the source gave;
            // with no filter the pass gave what the source gave, which is
            // nothing.
            num_vars_of_the_source: filters.last().map_or(0, |(_, stats)| stats.vars_processed),
            filters,
        });
    }
    Ok(PopDiversity {
        pops: totals.pops,
        num_vars_every_pop: totals.num_vars_every_pop,
        num_vars_of_the_pass: num_vars,
        stats: options.stats,
    })
}

/// It refuses a pass that was asked for no statistic at all.
///
/// Such a pass reads every variant of the source and computes nothing of
/// them, so it is refused at the call. The rule is here and not in each
/// binding crate, which is what the binding section of the `coding` skill
/// asks: a third language would otherwise have to write it a third time.
///
/// # Errors
///
/// A `stats` that holds none of the five statistics.
fn check_the_statistics(options: &DiversityOptions) -> Result<()> {
    if options.stats == DiversityStats::empty() {
        return Err(Error::DiversityWithNoStatistic);
    }
    Ok(())
}

/// It refuses a draw that no standardized value can be taken over.
///
/// # Errors
///
/// The folded spectrum asked for with no `num_called_alleles`, whose bins
/// are the counts of the rarer allele in a draw of that many, and a
/// `num_called_alleles` below 2.
fn check_the_draw(options: &DiversityOptions) -> Result<()> {
    match options.num_called_alleles {
        None => {
            if options.stats.contains(DiversityStats::FOLDED_SFS) {
                return Err(Error::DiversitySfsWithoutADraw);
            }
        }
        Some(num_called_alleles) => {
            if num_called_alleles < 2 {
                return Err(Error::DiversityDrawTooSmall { num_called_alleles });
            }
        }
    }
    Ok(())
}

/// The populations of the pass, each with the indices of its individuals
/// and with whether they are every individual of the reader in its order.
///
/// An empty `pops` is one population of every individual, which is what a
/// user who named no population gets.
///
/// # Errors
///
/// A population with no individual, an index that is not an individual of
/// the dataset, and an individual that is twice in one population.
fn pops_of_the_pass(pops: &[&[usize]], num_individuals: usize) -> Result<Vec<PopOfThePass>> {
    if pops.is_empty() {
        return Ok(vec![PopOfThePass {
            individuals: (0..num_individuals).collect(),
            is_all: true,
        }]);
    }
    let mut of_the_pass = Vec::with_capacity(pops.len());
    for (pop, individuals) in pops.iter().enumerate() {
        if individuals.is_empty() {
            return Err(Error::DiversityPopWithNoIndividual { pop });
        }
        let mut named = HashSet::with_capacity(individuals.len());
        for &individual in *individuals {
            if individual >= num_individuals {
                return Err(Error::DiversityIndividualNotInTheDataset {
                    pop,
                    individual,
                    num_individuals,
                });
            }
            if !named.insert(individual) {
                return Err(Error::DiversityIndividualAskedForTwice { pop, individual });
            }
        }
        of_the_pass.push(PopOfThePass {
            individuals: individuals.to_vec(),
            is_all: every_individual_in_order(individuals, num_individuals),
        });
    }
    Ok(of_the_pass)
}

/// It counts every row of a block into `totals`.
///
/// Natively the chunks of rows are read on the threads of rayon, as section
/// 3 of `docs/architecture.md` asks: no row reads another, each chunk
/// counts its own rows, and the chunks are added into `totals` in the order
/// of the block, so no count depends on how many threads read them. The
/// threads are those of the pool the caller is running in, and rayon's
/// global pool only when the caller is in none.
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
    of_the_pass: &OfThePass,
    totals: &mut Totals,
) -> Result<()> {
    use rayon::iter::ParallelIterator;
    use rayon::slice::ParallelSlice;

    let num_pops = of_the_pass.pops.len();
    let of_the_chunks: Result<Vec<Totals>> = block
        .gts
        .par_chunks(alleles_of_a_chunk(alleles_per_var))
        .map(|chunk| {
            let mut of_the_chunk = Totals::of(num_pops);
            add_the_rows(chunk, alleles_per_var, of_the_pass, &mut of_the_chunk)?;
            Ok(of_the_chunk)
        })
        .collect();
    match of_the_chunks {
        Ok(of_the_chunks) => {
            for of_the_chunk in &of_the_chunks {
                totals.add_the_chunk(of_the_chunk);
            }
            Ok(())
        }
        // The second pass costs a read of the block, and it is made only
        // where the block is refused and nothing of it is given.
        Err(of_a_thread) => {
            let mut read_again = Totals::of(num_pops);
            match add_the_chunks_one_by_one(block, alleles_per_var, of_the_pass, &mut read_again) {
                Err(of_the_first_row) => Err(of_the_first_row),
                // The rows are the same rows, so the second pass finds an
                // error too; the error of the threads is what is left if it
                // ever did not.
                Ok(()) => Err(of_a_thread),
            }
        }
    }
}

/// The same counts, with the chunks read one after another, which is what
/// wasm does: it has no threads.
///
/// # Errors
///
/// Those of the native [`add_the_block`], at the first row that has one.
#[cfg(target_family = "wasm")]
fn add_the_block(
    block: &Block,
    alleles_per_var: usize,
    of_the_pass: &OfThePass,
    totals: &mut Totals,
) -> Result<()> {
    add_the_chunks_one_by_one(block, alleles_per_var, of_the_pass, totals)
}

/// The chunks of the block read one after another, each into counts of its
/// own that are added into `totals` before the next is read: what wasm
/// runs, and what the threads fall back on to find the first row that is an
/// error.
///
/// # Errors
///
/// Those of [`add_the_block`], at the first row that has one.
fn add_the_chunks_one_by_one(
    block: &Block,
    alleles_per_var: usize,
    of_the_pass: &OfThePass,
    totals: &mut Totals,
) -> Result<()> {
    let mut of_the_chunk = Totals::of(of_the_pass.pops.len());
    for chunk in block.gts.chunks(alleles_of_a_chunk(alleles_per_var)) {
        of_the_chunk.forget_what_it_holds();
        add_the_rows(chunk, alleles_per_var, of_the_pass, &mut of_the_chunk)?;
        totals.add_the_chunk(&of_the_chunk);
    }
    Ok(())
}

/// It counts every row of one chunk into `totals`, one population after
/// another.
///
/// `gts` holds whole rows of `alleles_per_var` alleles each. A variant
/// counts for a population when the population called something at it and
/// called at least `min_num_individuals` genotypes, which is those alleles
/// over the ploidy: `min_called_alleles` is that threshold in alleles, so a
/// population of one half called genotype is below a threshold of 1
/// individual and two half called ones are at it.
///
/// # Errors
///
/// Those of [`add_the_block`], at the first row of the chunk that has one.
fn add_the_rows(
    gts: &[i8],
    alleles_per_var: usize,
    of_the_pass: &OfThePass,
    totals: &mut Totals,
) -> Result<()> {
    // One array of counts for each population, the populations by the
    // alleles of the row: the private alleles need every population's
    // counts at one variant at once, and the four other statistics read the
    // counts of the population that is in hand. `count_alleles_of` writes
    // every entry of the array it is given, so nothing is cleared here, and
    // the arrays are reused from row to row, so the loop over the rows
    // allocates nothing and what is kept never grows with the block.
    let mut of_each_pop = vec![OfAPopAtTheRow::none(); of_the_pass.pops.len()];
    // How many populations called each allele of the row, which says which
    // of them are private. It is filled for the row that is in hand and
    // reused by the next one.
    let mut num_pops_that_called: AlleleCounts = [0; 128];
    for row in gts.chunks_exact(alleles_per_var) {
        let mut every_pop = true;
        for ((of_the_pop, counted), at_the_row) in of_the_pass
            .pops
            .iter()
            .zip(totals.pops.iter_mut())
            .zip(of_each_pop.iter_mut())
        {
            // A population of every individual of the reader in its order
            // is counted by reading the row as it is, and the width of the
            // row says that it is the row of that reader: a population
            // built against another one would otherwise count individuals
            // it does not hold.
            let of_the_whole_row = of_the_pop.is_all
                && of_the_pop
                    .individuals
                    .len()
                    .saturating_mul(of_the_pass.ploidy)
                    == alleles_per_var;
            // The genotypes of the population are counted beside its
            // alleles when the pass gives F_IS, which reads the
            // heterozygous ones, and are not counted at all otherwise.
            let counts_the_gts = of_the_pass.heterozygosities.is_some();
            let (called_alleles, one_past_the_largest, gts) = if of_the_whole_row {
                // Counting the row as it is gives no bound on the alleles it
                // holds, so the whole of the counts is walked below; every
                // entry of them was written, so the ones above the largest
                // allele of the row hold 0.
                let called_alleles = count_alleles(row, &mut at_the_row.counts)?;
                let gts = if counts_the_gts {
                    Some(count_gts(row, of_the_pass.ploidy)?)
                } else {
                    None
                };
                (called_alleles, at_the_row.counts.len(), gts)
            } else if counts_the_gts {
                // One walk over the individuals of the population gives both
                // of its counts, where the two functions would look up each
                // of its genotypes twice.
                let (counted, gts) = count_alleles_and_gts_of(
                    row,
                    of_the_pass.ploidy,
                    &of_the_pop.individuals,
                    &mut at_the_row.counts,
                )?;
                (counted.called_alleles, counted.num_alleles, Some(gts))
            } else {
                let counted = count_alleles_of(
                    row,
                    of_the_pass.ploidy,
                    &of_the_pop.individuals,
                    &mut at_the_row.counts,
                )?;
                (counted.called_alleles, counted.num_alleles, None)
            };
            at_the_row.one_past_the_largest = one_past_the_largest;
            // A population that called nothing at the variant does not
            // count it whatever `min_num_individuals` is, so a threshold of
            // 0 does not put a variant with no data into the totals.
            if called_alleles == 0 || u64::from(called_alleles) < of_the_pass.min_called_alleles {
                every_pop = false;
                continue;
            }
            // One variant of the pass, and a pass of more than
            // 18446744073709551615 variants reads more rows than any source
            // holds.
            counted.num_vars = counted.num_vars.saturating_add(1);
            if of_the_pass.counts_the_alleles == CountsTheAlleles::Yes {
                let num_different = num_different_alleles(&at_the_row.counts, one_past_the_largest);
                // At most 128 different alleles at a variant, the entries a
                // count of them holds, so a sum above
                // 18446744073709551615 needs more than 1.4e17 variants,
                // more than any source holds.
                counted.num_alleles = counted.num_alleles.saturating_add(num_different);
                // The variant varies in the population when it called more
                // than one allele there, which is the variant
                // `docs/specs/stats.md` counts as variable: its major
                // allele frequency is below 1 exactly when a second allele
                // was called.
                if num_different > 1 {
                    counted.num_variable_vars = counted.num_variable_vars.saturating_add(1);
                }
            }
            // The genotypes of the population at this variant were counted
            // above exactly when the pass carries the two heterozygosities.
            if let (Some(heterozygosities), Some(gts)) = (of_the_pass.heterozygosities, gts) {
                heterozygosities.add_the_var(&at_the_row.counts, called_alleles, gts, counted);
            }
        }
        if every_pop {
            totals.num_vars_every_pop = totals.num_vars_every_pop.saturating_add(1);
            if of_the_pass.counts_the_private_alleles == CountsThePrivateAlleles::Yes {
                add_the_private_alleles(&of_each_pop, &mut num_pops_that_called, &mut totals.pops);
            }
        }
    }
    Ok(())
}

/// What one population called at the row that is in hand, which the private
/// alleles read across the populations.
#[derive(Debug, Clone)]
struct OfAPopAtTheRow {
    /// How often the population called each allele of the row, which
    /// [`count_alleles_of`] writes every entry of.
    counts: AlleleCounts,
    /// One past the largest allele it called there, the entry of `counts`
    /// the walks over them stop at: every entry from it up holds 0.
    one_past_the_largest: usize,
}

impl OfAPopAtTheRow {
    /// The counts of one population before any row is read, which hold no
    /// allele.
    fn none() -> OfAPopAtTheRow {
        OfAPopAtTheRow {
            counts: [0; 128],
            one_past_the_largest: 0,
        }
    }
}

/// It adds to every population the alleles it called at the row that no
/// other population of the pass called there.
///
/// It is called for a row that counted for every population and for no
/// other, which is the rule of "The private alleles" of
/// `docs/specs/diversity.md`: a population that called nothing at a variant
/// holds none of its alleles, so every allele of every other population
/// would be private there and the count would measure the missing data.
///
/// `of_each_pop` is what each population called at the row, in the order of
/// `totals`, and `num_pops_that_called` is the room for one count of the
/// populations for each allele of the row, whose contents when this is
/// called are of no interest. An allele is private to the population that
/// called it exactly when one population called it, so the populations that
/// called each allele are counted once and every population then reads that
/// count, which is one walk of the populations by the alleles of the row
/// and not one walk for each pair of populations.
fn add_the_private_alleles(
    of_each_pop: &[OfAPopAtTheRow],
    num_pops_that_called: &mut AlleleCounts,
    totals: &mut [OfAPop],
) {
    let one_past_the_largest = of_each_pop
        .iter()
        .map(|at_the_row| at_the_row.one_past_the_largest)
        .max()
        .unwrap_or(0);
    // The entries above the largest allele of this row are not read for it,
    // so what a row before left there stays and is cleared by the row that
    // reaches that far.
    for num_pops in num_pops_that_called.iter_mut().take(one_past_the_largest) {
        *num_pops = 0;
    }
    for at_the_row in of_each_pop {
        let called_and_counted = at_the_row
            .counts
            .iter()
            .zip(num_pops_that_called.iter_mut())
            .take(at_the_row.one_past_the_largest);
        for (called, num_pops) in called_and_counted {
            if *called > 0 {
                // A count that saturated would have to come from 4294967295
                // populations, and every value above 1 says the same thing
                // here, that the allele is not private.
                *num_pops = num_pops.saturating_add(1);
            }
        }
    }
    for (at_the_row, counted) in of_each_pop.iter().zip(totals.iter_mut()) {
        let num_private = at_the_row
            .counts
            .iter()
            .zip(num_pops_that_called.iter())
            .take(at_the_row.one_past_the_largest)
            .filter(|(called, num_pops)| **called > 0 && **num_pops == 1)
            // At most 128 private alleles at a variant, the entries a count
            // of them holds, so a sum above 18446744073709551615 needs more
            // than 1.4e17 variants, more than any source holds.
            .fold(0_u64, |num_private, _| num_private.saturating_add(1));
        counted.private_alleles = counted.private_alleles.saturating_add(num_private);
    }
}

/// How many different alleles a population called at one variant: the
/// counts above 0 among the first `one_past_the_largest` entries of
/// `counts`, from 1 where every copy it called is alike to the alleles the
/// variant has.
///
/// `counts` is what [`count_alleles_of`] left for the population at that
/// variant, and `one_past_the_largest` the bound it gave on the alleles it
/// wrote. An allele between 0 and the largest one that the population did
/// not call is not one of these, so a population whose only genotype is
/// `0/3` called two alleles and not four.
fn num_different_alleles(counts: &AlleleCounts, one_past_the_largest: usize) -> u64 {
    counts
        .iter()
        .take(one_past_the_largest)
        .filter(|count| **count > 0)
        // A count of the alleles of a variant holds 128 entries, so the
        // different ones are 128 at most and this never saturates.
        .fold(0_u64, |num_different, _| num_different.saturating_add(1))
}

#[cfg(test)]
mod fixtures {
    use crate::block::{Block, BlockReader};
    use crate::error::Result;
    use crate::filters::FilteringStats;
    use crate::variant::{ChromTable, Needs};

    /// The individuals of `pop1` of the worked example of "How it is
    /// verified" of "The number of alleles", `i1` and `i2`, as indices
    /// among the five individuals of the source.
    pub(super) const POP1: [usize; 2] = [0, 1];

    /// The individuals of `pop2` of that worked example, `i3`, `i4` and
    /// `i5`.
    pub(super) const POP2: [usize; 3] = [2, 3, 4];

    /// The six variants of five diploid individuals of the worked example,
    /// which is the one of `docs/specs/filters.md`, one row of genotypes
    /// each: `0/0 0/1 0/0 0/0 0/.`, `0/0 0/1 0/0 ./. 0/.`,
    /// `0/1 2/3 0/1 2/3 ./.`, five missing genotypes,
    /// `0/0 0/0 0/0 0/0 1/1` and `0/. ./. ./. ./. ./.`.
    pub(super) const THE_SIX_VARIANTS: [[i8; 10]; 6] = [
        [0, 0, 0, 1, 0, 0, 0, 0, 0, -1],
        [0, 0, 0, 1, 0, 0, -1, -1, 0, -1],
        [0, 1, 2, 3, 0, 1, 2, 3, -1, -1],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
        [0, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    ];

    /// A reader of the tests that gives the blocks it was built with, of
    /// five diploid individuals named `i1` to `i5` unless the test asks for
    /// others.
    #[derive(Debug)]
    pub(super) struct GivenBlocks {
        individuals: Vec<String>,
        ploidy: usize,
        chroms: ChromTable,
        /// The blocks it has not given yet, the next one last.
        left: Vec<Block>,
        /// What it was last asked to fill, which a test reads to see which
        /// fields the pass asked for.
        needs: Needs,
    }

    impl GivenBlocks {
        /// The reader over `blocks`, which it gives in their order, of the
        /// five diploid individuals of the worked example.
        pub(super) fn of(blocks: Vec<Block>) -> GivenBlocks {
            GivenBlocks::of_a_source_of(5, 2, blocks)
        }

        /// The same reader over a source of `num_individuals` individuals,
        /// named `i1` to `iN`, of the ploidy `ploidy`.
        pub(super) fn of_a_source_of(
            num_individuals: usize,
            ploidy: usize,
            blocks: Vec<Block>,
        ) -> GivenBlocks {
            let mut chroms = ChromTable::new();
            chroms.intern("chr1");
            let mut left = blocks;
            left.reverse();
            GivenBlocks {
                individuals: (1..=num_individuals)
                    .map(|number| format!("i{number}"))
                    .collect(),
                ploidy,
                chroms,
                left,
                needs: Needs::ALL,
            }
        }

        /// The fields the pass last asked it to fill.
        pub(super) fn needs(&self) -> Needs {
            self.needs
        }
    }

    impl BlockReader for GivenBlocks {
        fn next_block(&mut self) -> Result<Option<Block>> {
            Ok(self.left.pop())
        }

        fn individuals(&self) -> &[String] {
            &self.individuals
        }

        fn ploidy(&self) -> usize {
            self.ploidy
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

    /// The blocks of `num_vars_per_block` variants that hold `variants`,
    /// the rows of a worked example, one after another, of
    /// `num_individuals` individuals of the ploidy `ploidy`.
    pub(super) fn blocks_of(
        variants: &[&[i8]],
        num_individuals: usize,
        ploidy: usize,
        num_vars_per_block: usize,
    ) -> Vec<Block> {
        variants
            .chunks(num_vars_per_block)
            .map(|of_the_block| Block {
                num_vars: of_the_block.len(),
                num_individuals,
                ploidy,
                gts: of_the_block
                    .iter()
                    .flat_map(|row| row.iter())
                    .copied()
                    .collect(),
                chrom: None,
                pos: None,
                id: None,
                alleles: None,
                qual: None,
            })
            .collect()
    }

    /// Five rows of five diploid individuals that a long source repeats:
    /// `0/1 0/1 0/1 0/0 0/0`, three heterozygous genotypes of five;
    /// `0/1 1/2 0/. 2/2 ./.`, three alleles with a half called genotype
    /// among them; `0/0 0/0 0/1 ./. 1/1`; five missing genotypes, which
    /// counts for no population; and `0/2 0/1 1/1 0/0 2/.`.
    ///
    /// They are five and not four so that the chunks of 64 rows a pass cuts
    /// a block into do not hold the same rows as one another: with a cycle
    /// of four every full chunk would hold 16 of each row and have the same
    /// two sums as the next, and adding the chunks in another order would
    /// then give the same bits whatever that order was.
    ///
    /// The heterozygosities of a population over the five are thirds,
    /// sixths, sevenths and fifths, so the two sums behind its F_IS are
    /// float64 whose last bits move when the parts are added in another
    /// order.
    pub(super) const THE_FIVE_PATTERNS: [[i8; 10]; 5] = [
        [0, 1, 0, 1, 0, 1, 0, 0, 0, 0],
        [0, 1, 1, 2, 0, -1, 2, 2, -1, -1],
        [0, 0, 0, 0, 0, 1, -1, -1, 1, 1],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 2, 0, 1, 1, 1, 0, 0, 2, -1],
    ];

    /// A source of `num_vars` variants of the five diploid individuals, the
    /// five patterns of [`THE_FIVE_PATTERNS`] one after another, in blocks
    /// of `num_vars_per_block` variants.
    ///
    /// A test calls it with more than the 64 rows one chunk of a pass
    /// holds, which every other fixture of this module is below: with six
    /// variants the rows of a block are one chunk, the threads have nothing
    /// to share out and the order the chunks are added in is the order of
    /// the only one.
    pub(super) fn a_source_of_many_variants(
        num_vars: usize,
        num_vars_per_block: usize,
    ) -> GivenBlocks {
        let of_the_source = the_rows_of_many_variants(num_vars);
        let rows: Vec<&[i8]> = of_the_source.iter().map(|row| &row[..]).collect();
        GivenBlocks::of(blocks_of(&rows, 5, 2, num_vars_per_block))
    }

    /// The same source with the allele -2 in its eleventh row and the
    /// allele -3 in its hundred and fifty first, both below the missing
    /// one, which a reader of popnei never gives and which the counts of a
    /// variant refuse.
    ///
    /// The two rows are in different chunks of one block, so which of them
    /// a thread reaches first depends on how the chunks were shared out,
    /// and the pass has to give the error of the first of the two whatever
    /// happened.
    pub(super) fn a_source_with_two_rows_below_the_missing_allele(
        num_vars: usize,
        num_vars_per_block: usize,
    ) -> GivenBlocks {
        let mut of_the_source = the_rows_of_many_variants(num_vars);
        for (at, row) in of_the_source.iter_mut().enumerate() {
            let allele = match at {
                10 => -2,
                150 => -3,
                _ => continue,
            };
            if let Some(first) = row.first_mut() {
                *first = allele;
            }
        }
        let rows: Vec<&[i8]> = of_the_source.iter().map(|row| &row[..]).collect();
        GivenBlocks::of(blocks_of(&rows, 5, 2, num_vars_per_block))
    }

    /// The `num_vars` rows of [`THE_FIVE_PATTERNS`], one after another.
    fn the_rows_of_many_variants(num_vars: usize) -> Vec<[i8; 10]> {
        THE_FIVE_PATTERNS
            .iter()
            .copied()
            .cycle()
            .take(num_vars)
            .collect()
    }

    /// The six variants of the worked example, in blocks of
    /// `num_vars_per_block` variants of the five diploid individuals.
    pub(super) fn the_worked_example(num_vars_per_block: usize) -> GivenBlocks {
        let rows: Vec<&[i8]> = THE_SIX_VARIANTS.iter().map(|row| &row[..]).collect();
        GivenBlocks::of(blocks_of(&rows, 5, 2, num_vars_per_block))
    }
}

#[cfg(test)]
mod the_pass {
    use super::calc_pop_diversity_one_chunk_at_a_time;
    use super::fixtures::{
        GivenBlocks, POP1, POP2, a_source_of_many_variants,
        a_source_with_two_rows_below_the_missing_allele, blocks_of, the_worked_example,
    };
    use super::{DiversityOptions, DiversityStats, PopDiversity, calc_pop_diversity};
    use crate::error::Error;
    use crate::variant::Needs;

    /// The options of the worked example: every statistic, no draw, and a
    /// threshold of one called genotype.
    fn options_of(min_num_individuals: u32) -> DiversityOptions {
        DiversityOptions {
            stats: DiversityStats::ALL,
            num_called_alleles: None,
            min_num_individuals,
        }
    }

    /// The options with no spectrum, which is what a pass with no draw
    /// takes.
    fn options_with_no_draw(min_num_individuals: u32) -> DiversityOptions {
        DiversityOptions {
            stats: DiversityStats::NUM_ALLELES
                | DiversityStats::PRIVATE_ALLELES
                | DiversityStats::VARIABLE_VARS_RATIO
                | DiversityStats::FIS,
            num_called_alleles: None,
            min_num_individuals,
        }
    }

    /// The diversity of the two populations of the worked example over its
    /// six variants, in blocks of `num_vars_per_block`.
    fn of_the_worked_example(min_num_individuals: u32, num_vars_per_block: usize) -> PopDiversity {
        let mut reader = the_worked_example(num_vars_per_block);
        calc_pop_diversity(
            &mut reader,
            &[&POP1, &POP2],
            &options_with_no_draw(min_num_individuals),
        )
        .expect("the diversity of the worked example")
    }

    /// What a mean or a ratio of these tests may differ from the number of
    /// the spec by. Each of them is one count over another, both far below
    /// 2^53 and both exact in `f64`, and the quotients of the worked
    /// example, 2.25, 2, 0.75, 0.5 and 0.25, are exact in binary too, so
    /// the bound is the last bit of a number near 1 and not an allowance
    /// for any rounding.
    const OF_AN_EXACT_QUOTIENT: f64 = f64::EPSILON;

    /// It checks the mean alleles of one population, the allelic richness:
    /// the alleles it called over the variants that counted for it, which
    /// is the division a user of the Python layer is given and which this
    /// module gives the two counts of.
    fn assert_mean_num_alleles(diversity: &PopDiversity, pop: usize, mean: f64, what: &str) {
        let num_alleles = diversity.num_alleles(pop).expect("the alleles called");
        let num_vars = diversity.num_vars(pop).expect("the variants that counted");
        let found = num_alleles as f64 / num_vars as f64;

        assert!(
            (found - mean).abs() <= OF_AN_EXACT_QUOTIENT,
            "the mean alleles of {what} is {found}, and it is {mean}"
        );
    }

    /// It checks the mean private alleles of one population: the alleles no
    /// other population called over the variants that counted for every
    /// population, which is the divisor this statistic has and the other
    /// counts of the module do not.
    fn assert_mean_private_alleles(diversity: &PopDiversity, pop: usize, mean: f64, what: &str) {
        let private_alleles = diversity.private_alleles(pop).expect("the private alleles");
        let num_vars = diversity.num_vars_every_pop();
        let found = private_alleles as f64 / num_vars as f64;

        assert!(
            (found - mean).abs() <= OF_AN_EXACT_QUOTIENT,
            "the mean private alleles of {what} is {found}, and it is {mean}"
        );
    }

    /// What a value of F_IS of these tests may differ from the number of
    /// the spec by. The spec prints its F_IS to ten decimals, so a literal
    /// here is within 5e-11 of the value it rounds, and the two sums and
    /// the division of a number near 0.35 leave a few units of the last
    /// place of a float64, about 1e-16.
    const OF_TEN_DECIMALS: f64 = 1e-10;

    /// It checks that two results over the same variants hold the same
    /// numbers: every count equal, and the F_IS of every population the
    /// same bits.
    ///
    /// The bits and not a tolerance, because what these tests ask is
    /// whether the parts of a sum of float64 were added in the order of the
    /// variants: a sum whose parts were joined in another order is right to
    /// far more digits than any tolerance of the spec and is not the same
    /// number.
    fn assert_the_same_numbers(one: &PopDiversity, other: &PopDiversity, what: &str) {
        assert_eq!(
            one.num_pops(),
            other.num_pops(),
            "the populations of {what}"
        );
        assert_eq!(
            one.num_vars_every_pop(),
            other.num_vars_every_pop(),
            "the variants of every population of {what}"
        );
        for pop in 0..one.num_pops() {
            assert_eq!(
                one.num_vars(pop),
                other.num_vars(pop),
                "the variants of the population {pop} of {what}"
            );
            assert_eq!(
                one.num_alleles(pop),
                other.num_alleles(pop),
                "the alleles of the population {pop} of {what}"
            );
            assert_eq!(
                one.private_alleles(pop),
                other.private_alleles(pop),
                "the private alleles of the population {pop} of {what}"
            );
            assert_eq!(
                one.num_variable_vars(pop),
                other.num_variable_vars(pop),
                "the variable variants of the population {pop} of {what}"
            );
            let (obs_het, unbiased_exp_het, num_vars) = one
                .the_sums_behind_the_fis(pop)
                .expect("the sums behind the F_IS");
            let (of_the_other, unbiased_of_the_other, num_vars_of_the_other) = other
                .the_sums_behind_the_fis(pop)
                .expect("the sums behind the F_IS");

            assert_eq!(
                num_vars, num_vars_of_the_other,
                "the variants behind the F_IS of the population {pop} of {what}"
            );
            assert_eq!(
                obs_het.to_bits(),
                of_the_other.to_bits(),
                "the bits of the observed heterozygosities of the population {pop} of {what}, \
                 {obs_het} and {of_the_other}"
            );
            assert_eq!(
                unbiased_exp_het.to_bits(),
                unbiased_of_the_other.to_bits(),
                "the bits of the unbiased expected heterozygosities of the population {pop} of \
                 {what}, {unbiased_exp_het} and {unbiased_of_the_other}"
            );
        }
    }

    /// It checks the F_IS of one population: one minus its mean observed
    /// heterozygosity over its mean unbiased expected one, over the
    /// variants that counted for it and at which both of them exist.
    fn assert_fis(diversity: &PopDiversity, pop: usize, fis: f64, what: &str) {
        let found = diversity.fis(pop).expect("the F_IS");

        assert!(
            (found - fis).abs() <= OF_TEN_DECIMALS,
            "the F_IS of {what} is {found}, and it is {fis}"
        );
    }

    /// It checks the ratio of variable variants of one population: those
    /// where it called more than one allele over the variants that counted
    /// for it.
    fn assert_variable_vars_ratio(diversity: &PopDiversity, pop: usize, ratio: f64, what: &str) {
        let num_variable_vars = diversity
            .num_variable_vars(pop)
            .expect("the variable variants");
        let num_vars = diversity.num_vars(pop).expect("the variants that counted");
        let found = num_variable_vars as f64 / num_vars as f64;

        assert!(
            (found - ratio).abs() <= OF_AN_EXACT_QUOTIENT,
            "the ratio of variable variants of {what} is {found}, and it is {ratio}"
        );
    }

    /// The set of statistics holds the ones it was built from and no other,
    /// and `ALL` holds the five.
    #[test]
    fn a_set_of_statistics_holds_the_ones_it_was_built_from() {
        let two = DiversityStats::NUM_ALLELES | DiversityStats::FIS;

        assert!(two.contains(DiversityStats::NUM_ALLELES));
        assert!(two.contains(DiversityStats::FIS));
        assert!(!two.contains(DiversityStats::PRIVATE_ALLELES));
        assert!(!two.contains(DiversityStats::FOLDED_SFS));
        assert!(!two.contains(DiversityStats::VARIABLE_VARS_RATIO));
        assert!(DiversityStats::ALL.contains(two));
        assert!(DiversityStats::ALL.contains(DiversityStats::FOLDED_SFS));
        assert!(!DiversityStats::empty().contains(DiversityStats::NUM_ALLELES));
    }

    /// Four of the six variants count for each population of the worked
    /// example at a threshold of one called genotype: variant 4 has nothing
    /// called and variant 6 has one called allele in `pop1`, half a
    /// genotype, and none in `pop2`. "How it is verified" of "The number of
    /// alleles" of `docs/specs/diversity.md`.
    #[test]
    fn four_of_the_six_variants_of_the_worked_example_count_for_each_population() {
        let diversity = of_the_worked_example(1, 6);

        assert_eq!(diversity.num_pops(), 2);
        assert_eq!(diversity.num_vars(0), Some(4));
        assert_eq!(diversity.num_vars(1), Some(4));
        assert_eq!(diversity.num_vars_every_pop(), 4);
        assert_eq!(diversity.num_vars(2), None);
    }

    /// A threshold of 0 keeps the variant `pop1` called half a genotype at,
    /// the sixth, and still leaves out the fourth, which nothing was called
    /// at: a population that called nothing at a variant does not count it
    /// whatever `min_num_individuals` is, which "Its Python function" of
    /// the spec states.
    #[test]
    fn a_min_num_individuals_of_zero_leaves_out_a_variant_with_no_called_allele() {
        let diversity = of_the_worked_example(0, 6);

        assert_eq!(diversity.num_vars(0), Some(5));
        assert_eq!(diversity.num_vars(1), Some(4));
        assert_eq!(diversity.num_vars_every_pop(), 4);
    }

    /// The threshold is the called alleles over the ploidy, so two half
    /// called genotypes are one individual and not two. The variant is
    /// `0/. 0/. 0/0`, whose first two individuals called 2 alleles, which
    /// is 1 diploid individual: they are at a threshold of 1 and below one
    /// of 2, where a rule that counted each genotype with a called allele
    /// would keep them at 2.
    #[test]
    fn two_half_called_genotypes_are_one_individual_and_not_two() {
        let row: [i8; 6] = [0, -1, 0, -1, 0, 0];
        let rows: Vec<&[i8]> = vec![&row[..]];
        let of_the_halves: [usize; 2] = [0, 1];

        let at_one = {
            let mut reader = GivenBlocks::of_a_source_of(3, 2, blocks_of(&rows, 3, 2, 1));
            calc_pop_diversity(&mut reader, &[&of_the_halves], &options_with_no_draw(1))
                .expect("the diversity at a threshold of one individual")
        };
        let at_two = {
            let mut reader = GivenBlocks::of_a_source_of(3, 2, blocks_of(&rows, 3, 2, 1));
            calc_pop_diversity(&mut reader, &[&of_the_halves], &options_with_no_draw(2))
                .expect("the diversity at a threshold of two individuals")
        };

        assert_eq!(at_one.num_vars(0), Some(1));
        assert_eq!(at_two.num_vars(0), Some(0));
    }

    /// The counts do not depend on how the variants were cut into blocks,
    /// nor on how the rows of a block were shared out among the threads.
    #[test]
    fn the_size_of_the_blocks_does_not_change_the_counts() {
        for num_vars_per_block in [1, 2, 4, 6] {
            let diversity = of_the_worked_example(1, num_vars_per_block);

            assert_eq!(
                diversity.num_vars(0),
                Some(4),
                "pop1 in blocks of {num_vars_per_block}"
            );
            assert_eq!(
                diversity.num_vars(1),
                Some(4),
                "pop2 in blocks of {num_vars_per_block}"
            );
            assert_eq!(
                diversity.num_vars_every_pop(),
                4,
                "both pops in blocks of {num_vars_per_block}"
            );
            assert_eq!(
                diversity.num_alleles(0),
                Some(9),
                "the alleles of pop1 in blocks of {num_vars_per_block}"
            );
            assert_eq!(
                diversity.num_alleles(1),
                Some(8),
                "the alleles of pop2 in blocks of {num_vars_per_block}"
            );
            assert_eq!(
                diversity.num_variable_vars(0),
                Some(3),
                "the variable variants of pop1 in blocks of {num_vars_per_block}"
            );
            assert_eq!(
                diversity.num_variable_vars(1),
                Some(2),
                "the variable variants of pop2 in blocks of {num_vars_per_block}"
            );
            assert_eq!(
                diversity.private_alleles(0),
                Some(2),
                "the private alleles of pop1 in blocks of {num_vars_per_block}"
            );
            assert_eq!(
                diversity.private_alleles(1),
                Some(1),
                "the private alleles of pop2 in blocks of {num_vars_per_block}"
            );
        }
    }

    /// No population at all is one population of every individual of the
    /// reader, which counts the four variants something was called at. The
    /// five individuals of the worked example called 2, 2, 4 and 2 alleles
    /// there, 10 in all, and more than one allele at each of the four.
    #[test]
    fn no_population_is_one_population_of_every_individual() {
        let mut reader = the_worked_example(6);

        let diversity = calc_pop_diversity(&mut reader, &[], &options_with_no_draw(1))
            .expect("the diversity of every individual");

        assert_eq!(diversity.num_pops(), 1);
        assert_eq!(diversity.num_vars(0), Some(4));
        assert_eq!(diversity.num_vars_every_pop(), 4);
        assert_eq!(diversity.num_alleles(0), Some(10));
        assert_eq!(diversity.num_variable_vars(0), Some(4));
    }

    /// An individual in two populations is read by both, which
    /// `docs/specs/stats.md` allows and this module does not refuse.
    #[test]
    fn an_individual_in_two_populations_is_counted_in_both() {
        let mut reader = the_worked_example(6);
        let with_i2: [usize; 2] = [1, 2];

        let diversity =
            calc_pop_diversity(&mut reader, &[&POP1, &with_i2], &options_with_no_draw(1))
                .expect("the diversity of two populations that share an individual");

        assert_eq!(diversity.num_pops(), 2);
        assert_eq!(diversity.num_vars(0), Some(4));
        assert_eq!(diversity.num_vars(1), Some(4));
    }

    /// The pass reads the genotypes and no column of a block.
    #[test]
    fn the_pass_asks_the_reader_for_the_genotypes_alone() {
        let mut reader = the_worked_example(6);

        calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options_with_no_draw(1))
            .expect("the diversity of the worked example");

        assert_eq!(reader.needs(), Needs::GTS);
    }

    /// A population with no individual is refused: every statistic of a
    /// population is over the alleles its individuals called.
    #[test]
    fn a_population_with_no_individual_is_refused() {
        let mut reader = the_worked_example(6);
        let none: [usize; 0] = [];

        let error = calc_pop_diversity(&mut reader, &[&POP1, &none], &options_with_no_draw(1))
            .expect_err("a population with no individual");

        assert!(
            matches!(error, Error::DiversityPopWithNoIndividual { pop } if pop == 1),
            "{error}"
        );
    }

    /// An index that is not an individual of the dataset is refused, with
    /// the individuals the dataset has.
    #[test]
    fn an_index_that_is_not_an_individual_is_refused() {
        let mut reader = the_worked_example(6);
        let beyond: [usize; 2] = [2, 5];

        let error = calc_pop_diversity(&mut reader, &[&POP1, &beyond], &options_with_no_draw(1))
            .expect_err("an index that is not an individual");

        assert!(
            matches!(
                error,
                Error::DiversityIndividualNotInTheDataset {
                    pop: 1,
                    individual: 5,
                    num_individuals: 5
                }
            ),
            "{error}"
        );
    }

    /// An individual asked for twice in one population is refused: its
    /// alleles would be counted twice there.
    #[test]
    fn an_individual_asked_for_twice_in_one_population_is_refused() {
        let mut reader = the_worked_example(6);
        let twice: [usize; 3] = [2, 3, 2];

        let error = calc_pop_diversity(&mut reader, &[&twice], &options_with_no_draw(1))
            .expect_err("an individual asked for twice");

        assert!(
            matches!(
                error,
                Error::DiversityIndividualAskedForTwice {
                    pop: 0,
                    individual: 2
                }
            ),
            "{error}"
        );
    }

    /// A pass that gives no variant is refused, with whether the source
    /// held none or the steps kept none.
    #[test]
    fn a_pass_that_gives_no_variant_is_refused() {
        let mut reader = GivenBlocks::of(Vec::new());

        let error = calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options_with_no_draw(1))
            .expect_err("a pass with no variant");

        assert!(
            matches!(error, Error::PassGaveNoVariant { num_vars_of_the_source, ref filters }
                if num_vars_of_the_source == 0 && filters.is_empty()),
            "{error}"
        );
    }

    /// The folded spectrum asked for with no draw size is refused: its bins
    /// are the counts of the rarer allele in a draw of that many alleles.
    #[test]
    fn the_folded_spectrum_without_a_draw_is_refused() {
        let mut reader = the_worked_example(6);

        let error = calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options_of(1))
            .expect_err("the spectrum with no draw");

        assert!(matches!(error, Error::DiversitySfsWithoutADraw), "{error}");
    }

    /// A pass asked for no statistic is refused: it would read every
    /// variant of the source and compute nothing of them. The rule is here
    /// and not in each binding crate, which the `coding` skill asks of a
    /// refusal that is about the calculation and not about a language.
    #[test]
    fn a_pass_asked_for_no_statistic_is_refused() {
        let mut reader = the_worked_example(6);
        let options = DiversityOptions {
            stats: DiversityStats::empty(),
            num_called_alleles: None,
            min_num_individuals: 1,
        };

        let error = calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options)
            .expect_err("a pass asked for no statistic");

        assert!(matches!(error, Error::DiversityWithNoStatistic), "{error}");
    }

    /// A draw of one allele is refused: it shows one allele whatever the
    /// population holds.
    #[test]
    fn a_draw_of_fewer_than_two_alleles_is_refused() {
        let mut reader = the_worked_example(6);
        let options = DiversityOptions {
            num_called_alleles: Some(1),
            ..options_of(1)
        };

        let error = calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options)
            .expect_err("a draw of one allele");

        assert!(
            matches!(
                error,
                Error::DiversityDrawTooSmall {
                    num_called_alleles: 1
                }
            ),
            "{error}"
        );
    }
    /// `pop1` of the worked example called 2, 2, 4 and 1 different alleles
    /// at the variants 1, 2, 3 and 5 that count for it, 9 in all and a mean
    /// of 2.25, and `pop2` called 1, 1, 4 and 2, 8 in all and a mean of 2.
    /// "How it is verified" of "The number of alleles" of
    /// `docs/specs/diversity.md`.
    #[test]
    fn the_populations_of_the_worked_example_called_nine_and_eight_alleles() {
        let diversity = of_the_worked_example(1, 6);

        assert_eq!(diversity.num_alleles(0), Some(9));
        assert_eq!(diversity.num_alleles(1), Some(8));
        assert_eq!(diversity.num_alleles(2), None);
        assert_mean_num_alleles(&diversity, 0, 2.25, "pop1");
        assert_mean_num_alleles(&diversity, 1, 2.0, "pop2");
    }

    /// `pop1` of the worked example called more than one allele at the
    /// variants 1, 2 and 3 and not at 5, 3 of the 4 that count for it and a
    /// ratio of 0.75, and `pop2` at 3 and 5 and not at 1 and 2, 2 of 4 and
    /// 0.5. "How it is verified" of "The variable variants" of
    /// `docs/specs/diversity.md`.
    #[test]
    fn three_variants_vary_in_pop1_of_the_worked_example_and_two_in_pop2() {
        let diversity = of_the_worked_example(1, 6);

        assert_eq!(diversity.num_variable_vars(0), Some(3));
        assert_eq!(diversity.num_variable_vars(1), Some(2));
        assert_eq!(diversity.num_variable_vars(2), None);
        assert_variable_vars_ratio(&diversity, 0, 0.75, "pop1");
        assert_variable_vars_ratio(&diversity, 1, 0.5, "pop2");
    }

    /// The alleles a population called are the counts above 0 and not the
    /// alleles up to the largest one it called: a population whose only
    /// genotype is `0/3` called 2 alleles and not 4, and one whose only
    /// genotype is `3/3` called 1 and does not vary. The two variants are
    /// `0/3 0/0 0/0` and `3/3 0/0 0/0`, counted for the first individual
    /// alone and for the three of the reader in their order, which is the
    /// row read as it is.
    ///
    /// The second variant is what holds the count of the population to the
    /// alleles it called: the first individual called the allele 3 there
    /// and nothing below it, so a count that walked the alleles up to the
    /// largest one would give it 4 alleles and call the variant variable,
    /// where it called one allele and does not vary.
    #[test]
    fn an_allele_a_population_did_not_call_is_not_counted() {
        let of_the_two: [i8; 6] = [0, 3, 0, 0, 0, 0];
        let of_the_threes: [i8; 6] = [3, 3, 0, 0, 0, 0];
        let rows: Vec<&[i8]> = vec![&of_the_two[..], &of_the_threes[..]];
        let of_the_first: [usize; 1] = [0];
        let of_the_three: [usize; 3] = [0, 1, 2];
        let mut reader = GivenBlocks::of_a_source_of(3, 2, blocks_of(&rows, 3, 2, 2));

        let diversity = calc_pop_diversity(
            &mut reader,
            &[&of_the_first, &of_the_three],
            &options_with_no_draw(1),
        )
        .expect("the diversity of two variants with a gap in their alleles");

        // The first individual called the alleles 0 and 3 at the first
        // variant and the allele 3 alone at the second, 3 in all; the three
        // of the reader called 0 and 3 at each, 4 in all.
        assert_eq!(diversity.num_alleles(0), Some(3));
        assert_eq!(diversity.num_alleles(1), Some(4));
        assert_eq!(diversity.num_variable_vars(0), Some(1));
        assert_eq!(diversity.num_variable_vars(1), Some(2));
    }

    /// A statistic that was not asked for is not counted and has no value,
    /// and the one that was asked for has the same number as in a pass that
    /// asked for both.
    #[test]
    fn a_statistic_that_was_not_asked_for_has_no_value() {
        let of_the_alleles = {
            let mut reader = the_worked_example(6);
            let options = DiversityOptions {
                stats: DiversityStats::NUM_ALLELES,
                num_called_alleles: None,
                min_num_individuals: 1,
            };
            calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options)
                .expect("the diversity with the alleles called alone")
        };
        let of_the_variable_vars = {
            let mut reader = the_worked_example(6);
            let options = DiversityOptions {
                stats: DiversityStats::VARIABLE_VARS_RATIO,
                num_called_alleles: None,
                min_num_individuals: 1,
            };
            calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options)
                .expect("the diversity with the variable variants alone")
        };

        assert_eq!(of_the_alleles.num_alleles(0), Some(9));
        assert_eq!(of_the_alleles.num_variable_vars(0), None);
        assert_eq!(of_the_variable_vars.num_variable_vars(0), Some(3));
        assert_eq!(of_the_variable_vars.num_alleles(0), None);
        // The variants that count for a population are counted whatever the
        // pass was asked for, since every statistic is over them.
        assert_eq!(of_the_variable_vars.num_vars(0), Some(4));
    }

    /// `pop1` of the worked example called an allele no other population
    /// called at the variants 1 and 2, the allele 1 of each, 2 in all over
    /// the 4 variants that counted for both populations and a mean of 0.5;
    /// `pop2` called one at variant 5, the allele 1 there, 1 in all and a
    /// mean of 0.25. At variant 3 both populations called the four alleles
    /// `0`, `1`, `2` and `3`, so neither holds one the other does not.
    /// "What it gives" of "The private alleles" of
    /// `docs/specs/diversity.md`.
    #[test]
    fn the_populations_of_the_worked_example_have_two_and_one_private_alleles() {
        let diversity = of_the_worked_example(1, 6);

        assert_eq!(diversity.private_alleles(0), Some(2));
        assert_eq!(diversity.private_alleles(1), Some(1));
        assert_eq!(diversity.private_alleles(2), None);
        assert_eq!(diversity.num_vars_every_pop(), 4);
        assert_mean_private_alleles(&diversity, 0, 0.5, "pop1");
        assert_mean_private_alleles(&diversity, 1, 0.25, "pop2");
    }

    /// With one population every allele it called is private, since there
    /// is no other population to hold it, so the private alleles are the
    /// alleles called: the 10 of the five individuals of the worked example
    /// over its four variants. "The cases" of `docs/specs/diversity.md`.
    #[test]
    fn with_one_population_every_allele_it_called_is_private() {
        let mut reader = the_worked_example(6);

        let diversity = calc_pop_diversity(&mut reader, &[], &options_with_no_draw(1))
            .expect("the diversity of one population");

        assert_eq!(diversity.num_alleles(0), Some(10));
        assert_eq!(diversity.private_alleles(0), Some(10));
        assert_eq!(diversity.num_vars_every_pop(), 4);
    }

    /// A variant where one population has too little called counts for the
    /// private alleles of no population, since the alleles of the others
    /// would all be private there and the count would measure the missing
    /// data. At a threshold of 0 the sixth variant of the worked example,
    /// `0/. ./. ./. ./. ./.`, counts for `pop1`, which called the allele 0
    /// there and no other population called it, and not for `pop2`, which
    /// called nothing: `pop1` keeps its 2 private alleles over the 4
    /// variants both populations counted, while 5 variants counted for it.
    #[test]
    fn a_variant_one_population_is_short_at_is_out_of_every_private_count() {
        let diversity = of_the_worked_example(0, 6);

        assert_eq!(diversity.num_vars(0), Some(5));
        assert_eq!(diversity.num_vars_every_pop(), 4);
        assert_eq!(diversity.private_alleles(0), Some(2));
        assert_eq!(diversity.private_alleles(1), Some(1));
    }

    /// An allele no population called is private to none of them, which a
    /// count that walked the alleles up to the largest one instead of the
    /// ones called would get wrong. The variant is `0/3 0/0 0/0`, with the
    /// first individual one population and the other two the second: the
    /// alleles 1 and 2 were called by nobody, the allele 0 by both
    /// populations, and only the allele 3 of the first population is
    /// private.
    #[test]
    fn an_allele_no_population_called_is_private_to_none_of_them() {
        let row: [i8; 6] = [0, 3, 0, 0, 0, 0];
        let rows: Vec<&[i8]> = vec![&row[..]];
        let of_the_first: [usize; 1] = [0];
        let of_the_others: [usize; 2] = [1, 2];
        let mut reader = GivenBlocks::of_a_source_of(3, 2, blocks_of(&rows, 3, 2, 1));

        let diversity = calc_pop_diversity(
            &mut reader,
            &[&of_the_first, &of_the_others],
            &options_with_no_draw(1),
        )
        .expect("the diversity of a variant with a gap in its alleles");

        assert_eq!(diversity.private_alleles(0), Some(1));
        assert_eq!(diversity.private_alleles(1), Some(0));
    }

    /// The private alleles are counted when they were asked for and not
    /// otherwise, and a pass asked for them alone gives the same number as
    /// one asked for every statistic.
    #[test]
    fn the_private_alleles_are_counted_only_when_they_were_asked_for() {
        let of_the_private_alleles = {
            let mut reader = the_worked_example(6);
            let options = DiversityOptions {
                stats: DiversityStats::PRIVATE_ALLELES,
                num_called_alleles: None,
                min_num_individuals: 1,
            };
            calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options)
                .expect("the diversity with the private alleles alone")
        };
        let of_the_alleles = {
            let mut reader = the_worked_example(6);
            let options = DiversityOptions {
                stats: DiversityStats::NUM_ALLELES,
                num_called_alleles: None,
                min_num_individuals: 1,
            };
            calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options)
                .expect("the diversity with the alleles called alone")
        };

        assert_eq!(of_the_private_alleles.private_alleles(0), Some(2));
        assert_eq!(of_the_private_alleles.private_alleles(1), Some(1));
        assert_eq!(of_the_private_alleles.num_alleles(0), None);
        assert_eq!(of_the_alleles.private_alleles(0), None);
    }
    /// `pop1` of the worked example has observed heterozygosities of 0.5,
    /// 0.5, 1 and 0 at the variants 1, 2, 3 and 5 and unbiased expected
    /// heterozygosities of 0.5, 0.5, 1 and 0, both means 0.5, so its F_IS
    /// is 0; `pop2` has 0, 0, 1 and 0 against 0, 0, 1 and 0.5333333333,
    /// means of 0.25 and 0.3833333333, so its F_IS is 0.3478260870. The
    /// eight per variant values are the ones the worked example of
    /// `docs/specs/stats.md` tabulates. "How it is verified" of "The
    /// inbreeding coefficient F_IS" of `docs/specs/diversity.md`.
    #[test]
    fn the_fis_of_the_worked_example_is_zero_in_pop1_and_0_3478_in_pop2() {
        let diversity = of_the_worked_example(1, 6);

        assert_fis(&diversity, 0, 0.0, "pop1");
        assert_fis(&diversity, 1, 0.3478260870, "pop2");
        assert_eq!(diversity.fis(2), None);
    }

    /// The variants that counted for a population do not depend on how the
    /// source was cut into blocks, and its F_IS is the same within the
    /// tolerance of the spec. It is the same to every digit here, where the
    /// six variants of the worked example fall in one chunk of rows at
    /// every size tried; the two sums behind F_IS are sums of float64
    /// grouped by chunk inside each block, so a source of more than
    /// `ROWS_PER_CHUNK` variants read in blocks of different sizes groups
    /// them differently and moves the last bits. On the panel of 1200
    /// variants that move is 5.2e-14 relative, measured on 24 September
    /// 2026 by reading it whole and in blocks of 7, against the 1e-12 the
    /// spec allows.
    #[test]
    fn the_size_of_the_blocks_does_not_change_the_fis() {
        for num_vars_per_block in [1, 2, 4, 6] {
            let diversity = of_the_worked_example(1, num_vars_per_block);

            assert_fis(
                &diversity,
                0,
                0.0,
                &format!("pop1 in blocks of {num_vars_per_block}"),
            );
            assert_fis(
                &diversity,
                1,
                0.3478260870,
                &format!("pop2 in blocks of {num_vars_per_block}"),
            );
        }
    }

    /// A variant that counted for a population and that it called no whole
    /// genotype at is out of both of its means: the observed heterozygosity
    /// divides by the called genotypes and has none there. At a threshold
    /// of no individual the sixth variant of the worked example,
    /// `0/. ./. ./. ./. ./.`, counts for `pop1`, which called one allele of
    /// a half called genotype there, and `pop1` keeps the F_IS of the four
    /// variants it called whole genotypes at.
    #[test]
    fn a_variant_with_no_called_genotype_is_out_of_the_fis() {
        let diversity = of_the_worked_example(0, 6);

        assert_eq!(diversity.num_vars(0), Some(5));
        assert_fis(&diversity, 0, 0.0, "pop1 at a threshold of no individual");
        assert_fis(
            &diversity,
            1,
            0.3478260870,
            "pop2 at a threshold of no individual",
        );
    }

    /// A population for which no variant counted has NaN in `fis` and is
    /// not an error, which "The cases" of `docs/specs/diversity.md` states.
    /// At a threshold of three called genotypes no variant counts for
    /// `pop1`, which holds two individuals; the fifth counts for `pop2`,
    /// whose three individuals called six alleles there, and `pop2` has an
    /// observed heterozygosity of 0 and an unbiased expected one of
    /// 0.5333333333 at it, so its F_IS is 1.
    #[test]
    fn a_population_no_variant_counted_for_has_no_fis() {
        let diversity = of_the_worked_example(3, 6);

        assert_eq!(diversity.num_vars(0), Some(0));
        assert!(
            diversity
                .fis(0)
                .expect("the F_IS of a population no variant counted for")
                .is_nan()
        );
        assert_eq!(diversity.num_vars(1), Some(1));
        assert_fis(
            &diversity,
            1,
            1.0,
            "pop2 at a threshold of three individuals",
        );
    }

    /// A population whose every variant holds one allele has a mean
    /// unbiased expected heterozygosity of 0 and no F_IS, which is NaN:
    /// "What it gives" of "The inbreeding coefficient F_IS" of
    /// `docs/specs/diversity.md`. The two variants are `0/0 0/0 0/0` and
    /// `1/1 1/1 1/1`, of one population of the three individuals of the
    /// reader.
    #[test]
    fn a_population_of_one_allele_at_every_variant_has_no_fis() {
        let of_the_zeros: [i8; 6] = [0, 0, 0, 0, 0, 0];
        let of_the_ones: [i8; 6] = [1, 1, 1, 1, 1, 1];
        let rows: Vec<&[i8]> = vec![&of_the_zeros[..], &of_the_ones[..]];
        let mut reader = GivenBlocks::of_a_source_of(3, 2, blocks_of(&rows, 3, 2, 2));

        let diversity = calc_pop_diversity(&mut reader, &[], &options_with_no_draw(1))
            .expect("the diversity of a population of one allele at every variant");

        assert_eq!(diversity.num_vars(0), Some(2));
        assert!(
            diversity
                .fis(0)
                .expect("the F_IS of a population with no diversity")
                .is_nan()
        );
    }

    /// At ploidy 1 no genotype can be heterozygous, so the observed
    /// heterozygosity is 0 at every variant and F_IS would be 1 wherever
    /// the population has any diversity. popnei gives NaN there instead,
    /// which "What it gives" of "The inbreeding coefficient F_IS" of
    /// `docs/specs/diversity.md` states. The variant is `0 1 1` of three
    /// haploid individuals, at which the population called two alleles.
    #[test]
    fn a_haploid_population_has_no_fis_instead_of_one() {
        let row: [i8; 3] = [0, 1, 1];
        let rows: Vec<&[i8]> = vec![&row[..]];
        let mut reader = GivenBlocks::of_a_source_of(3, 1, blocks_of(&rows, 3, 1, 1));

        let diversity = calc_pop_diversity(&mut reader, &[], &options_with_no_draw(1))
            .expect("the diversity of a haploid population");

        assert_eq!(diversity.num_vars(0), Some(1));
        assert_eq!(diversity.num_alleles(0), Some(2));
        assert!(
            diversity
                .fis(0)
                .expect("the F_IS of a haploid population")
                .is_nan()
        );
    }

    /// A population that holds more heterozygous genotypes than random
    /// pairing would give has a negative F_IS. The variant is `0/1 0/1` of
    /// two diploid individuals: both of its called genotypes are
    /// heterozygous, so its observed heterozygosity is 1, and 2 of its 4
    /// called alleles are the allele 0 and 2 the allele 1, so the chance
    /// that two of them drawn without replacement are alike is twice
    /// (2/4)(1/3), a third, and its unbiased expected heterozygosity is two
    /// thirds. Its F_IS is 1 - 1 / (2/3), -0.5.
    #[test]
    fn a_population_of_heterozygous_genotypes_alone_has_a_negative_fis() {
        let row: [i8; 4] = [0, 1, 0, 1];
        let rows: Vec<&[i8]> = vec![&row[..]];
        let of_the_two: [usize; 2] = [0, 1];
        let mut reader = GivenBlocks::of_a_source_of(2, 2, blocks_of(&rows, 2, 2, 1));

        let diversity = calc_pop_diversity(&mut reader, &[&of_the_two], &options_with_no_draw(1))
            .expect("the diversity of two heterozygous genotypes");

        assert_fis(
            &diversity,
            0,
            -0.5,
            "a population of two heterozygous genotypes",
        );
    }

    /// F_IS is counted when it was asked for and not otherwise, and a pass
    /// asked for it alone gives the same number as one asked for every
    /// statistic.
    #[test]
    fn the_fis_is_counted_only_when_it_was_asked_for() {
        let of_the_fis = {
            let mut reader = the_worked_example(6);
            let options = DiversityOptions {
                stats: DiversityStats::FIS,
                num_called_alleles: None,
                min_num_individuals: 1,
            };
            calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options)
                .expect("the diversity with the F_IS alone")
        };
        let of_the_alleles = {
            let mut reader = the_worked_example(6);
            let options = DiversityOptions {
                stats: DiversityStats::NUM_ALLELES,
                num_called_alleles: None,
                min_num_individuals: 1,
            };
            calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options)
                .expect("the diversity with the alleles called alone")
        };

        assert_fis(
            &of_the_fis,
            0,
            0.0,
            "pop1 in a pass asked for the F_IS alone",
        );
        assert_fis(
            &of_the_fis,
            1,
            0.3478260870,
            "pop2 in a pass asked for the F_IS alone",
        );
        assert_eq!(of_the_fis.num_alleles(0), None);
        assert_eq!(of_the_alleles.fis(0), None);
    }

    /// The variants of the pass are the ones the reader gave, which a caller
    /// reports beside what each filter was given and kept. They are not the
    /// variants that counted for a population: the worked example has six,
    /// and four of them count for each of its two populations at a threshold
    /// of one called genotype, variant 4 having nothing called and variant 6
    /// too little in `pop1`.
    #[test]
    fn the_variants_of_the_pass_are_the_ones_the_reader_gave() {
        let diversity = of_the_worked_example(1, 6);

        assert_eq!(diversity.num_vars_of_the_pass(), 6);
        assert_eq!(diversity.num_vars(0), Some(4));
        assert_eq!(diversity.num_vars(1), Some(4));
        assert_eq!(
            of_the_worked_example(1, 2).num_vars_of_the_pass(),
            6,
            "the variants of the pass do not depend on the size of a block"
        );
    }
    /// A variant that carries one of the two heterozygosities and not the
    /// other is out of both means, which is what keeps the two sums over
    /// one set of variants and gives them one divisor. The variants are
    /// `0/1 0/1` and `0/. 1/.` of two diploid individuals: the second
    /// counts for the population, which called 2 alleles there, 1
    /// individual of the ploidy, and has an unbiased expected
    /// heterozygosity of 1, its two called alleles being different; it has
    /// no observed heterozygosity, both of its genotypes being half called
    /// and neither of them called.
    ///
    /// So the F_IS is the first variant's alone, 1 - 1 / (2/3), -0.5. A
    /// pass that read a missing heterozygosity as 0 would add (0, 1) to the
    /// sums and give 1 - 0.5 / (5/6), 0.4.
    #[test]
    fn a_variant_with_only_the_expected_heterozygosity_is_out_of_the_fis() {
        let of_the_hets: [i8; 4] = [0, 1, 0, 1];
        let of_the_halves: [i8; 4] = [0, -1, 1, -1];
        let rows: Vec<&[i8]> = vec![&of_the_hets[..], &of_the_halves[..]];
        let of_the_two: [usize; 2] = [0, 1];
        let mut reader = GivenBlocks::of_a_source_of(2, 2, blocks_of(&rows, 2, 2, 2));

        let diversity = calc_pop_diversity(&mut reader, &[&of_the_two], &options_with_no_draw(1))
            .expect("the diversity of a heterozygous variant and a half called one");

        assert_eq!(diversity.num_vars(0), Some(2));
        assert_fis(
            &diversity,
            0,
            -0.5,
            "a population whose second variant has no observed heterozygosity",
        );
    }

    /// A population whose every counted variant has one of the two
    /// heterozygosities and not the other has no F_IS and gets NaN: the
    /// variant `0/. 1/.` of two diploid individuals counts, the population
    /// having called 2 alleles at it, and no whole genotype of it was
    /// called.
    #[test]
    fn a_population_with_no_called_genotype_at_any_variant_has_no_fis() {
        let of_the_halves: [i8; 4] = [0, -1, 1, -1];
        let rows: Vec<&[i8]> = vec![&of_the_halves[..]];
        let of_the_two: [usize; 2] = [0, 1];
        let mut reader = GivenBlocks::of_a_source_of(2, 2, blocks_of(&rows, 2, 2, 1));

        let diversity = calc_pop_diversity(&mut reader, &[&of_the_two], &options_with_no_draw(1))
            .expect("the diversity of a variant of two half called genotypes");

        assert_eq!(diversity.num_vars(0), Some(1));
        assert_eq!(diversity.num_alleles(0), Some(2));
        assert!(
            diversity
                .fis(0)
                .expect("the F_IS of a population of half called genotypes")
                .is_nan()
        );
    }

    /// The exponent of the unbiased expected heterozygosity is the ploidy
    /// of the variants, so a tetraploid population draws four gene copies
    /// without replacement and not two. The three variants of two
    /// tetraploid individuals are `0/0/1/1 0/1/1/1`, `0/0/0/0 0/0/0/0` and
    /// `0/0/1/2 ./././.`.
    ///
    /// At the first both genotypes are heterozygous, so the observed
    /// heterozygosity is 1; of its 8 called alleles 3 are the allele 0 and
    /// 5 the allele 1, and four copies drawn from them are all alike only
    /// if all four are the allele 1, (5/8)(4/7)(3/6)(2/5), a fourteenth, so
    /// the unbiased expected heterozygosity is 13/14. At the second every
    /// copy is the allele 0, so the two are 0 and 0. At the third the one
    /// called genotype is heterozygous, so the observed one is 1, and no
    /// allele of its 4 called ones was called four times, so the unbiased
    /// one is 1. The means are 2/3 and 9/14, and the F_IS is 1 - 28/27,
    /// -0.0370370370.
    ///
    /// An exponent of 2 at the same variants would give 15/28, 0 and 5/6,
    /// a mean of 115/252, and an F_IS of -53/115, -0.4608695652.
    #[test]
    fn a_tetraploid_population_draws_four_gene_copies_and_not_two() {
        let of_the_hets: [i8; 8] = [0, 0, 1, 1, 0, 1, 1, 1];
        let of_the_zeros: [i8; 8] = [0, 0, 0, 0, 0, 0, 0, 0];
        let of_the_three_alleles: [i8; 8] = [0, 0, 1, 2, -1, -1, -1, -1];
        let rows: Vec<&[i8]> = vec![
            &of_the_hets[..],
            &of_the_zeros[..],
            &of_the_three_alleles[..],
        ];
        let of_the_two: [usize; 2] = [0, 1];
        let mut reader = GivenBlocks::of_a_source_of(2, 4, blocks_of(&rows, 2, 4, 3));

        let diversity = calc_pop_diversity(&mut reader, &[&of_the_two], &options_with_no_draw(1))
            .expect("the diversity of a tetraploid population");

        assert_eq!(diversity.num_vars(0), Some(3));
        // 2 alleles at the first variant, 1 at the second and 3 at the
        // third.
        assert_eq!(diversity.num_alleles(0), Some(6));
        assert_eq!(diversity.num_variable_vars(0), Some(2));
        assert_fis(&diversity, 0, -0.0370370370, "a tetraploid population");
    }

    /// The chunks of a block read one after another, which is the reduction
    /// WebAssembly runs, count what the threads of a native build count.
    ///
    /// No pass of a native build reaches that reduction: there the chunks
    /// are read on the threads of rayon, and the rows are read again one
    /// chunk at a time only where a block is refused, which gives no
    /// numbers. So without this test the whole of it is compiled by the
    /// cargo suites and run by none of them, and only the node suite, which
    /// runs the wasm build over the panel, would see it wrong.
    #[test]
    fn the_chunks_read_one_after_another_count_what_the_threads_count() {
        let of_the_threads = {
            let mut reader = a_source_of_many_variants(200, 200);
            calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options_with_no_draw(1))
                .expect("the diversity read on the threads")
        };
        let one_chunk_at_a_time = {
            let mut reader = a_source_of_many_variants(200, 200);
            calc_pop_diversity_one_chunk_at_a_time(
                &mut reader,
                &[&POP1, &POP2],
                &options_with_no_draw(1),
            )
            .expect("the diversity read one chunk at a time")
        };

        assert_eq!(of_the_threads.num_vars(0), Some(160));
        assert_the_same_numbers(
            &of_the_threads,
            &one_chunk_at_a_time,
            "the chunks read one after another",
        );
    }

    /// The error of a block is the one of its first bad row, whichever row
    /// a thread reached first: a user who reports a damaged file has to get
    /// the same message every time. The block holds 200 rows, four chunks,
    /// with the allele -2 in the eleventh and the allele -3 in the hundred
    /// and fifty first, and the pass is run 40 times because which of the
    /// two a thread reaches first is the scheduling of that run.
    #[test]
    fn the_error_of_a_block_is_the_one_of_its_first_bad_row_every_time() {
        for run in 0..40 {
            let mut reader = a_source_with_two_rows_below_the_missing_allele(200, 200);

            let error = calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options_with_no_draw(1))
                .expect_err("a block with two rows below the missing allele");

            assert!(
                matches!(error, Error::AlleleBelowTheMissingOne { allele: -2 }),
                "the run {run} gave {error}"
            );
        }
    }

    /// The numbers of a pass do not depend on how many threads read the
    /// rows of a block. The source holds 2000 variants in one block, 32
    /// chunks of rows, so the threads have chunks to share out; every other
    /// fixture of this module holds six variants or fewer, which is one
    /// chunk and nothing to share.
    ///
    /// The bits of the two sums behind F_IS are compared and not a
    /// tolerance: the chunks are added in the order of the block, and
    /// rayon's own `reduce` over the same chunks joins them in a tree whose
    /// shape follows the threads of the pool, which gives a number right to
    /// far more digits than any tolerance of the spec and not the same one.
    /// With that `reduce` in place of the ordered addition the sums of
    /// `pop1` differ between 1 thread and 2 here, and they do not at 200
    /// variants, 4 chunks, where rayon splits the same way whatever the
    /// pool.
    ///
    /// The pools are built here and are not rayon's global one, which has
    /// one thread per core of the machine. rayon is a dependency of the
    /// targets that are not wasm, so this test is compiled for those alone.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_number_of_threads_does_not_change_the_numbers() {
        let in_a_pool = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            pool.install(|| {
                let mut reader = a_source_of_many_variants(2000, 2000);
                calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options_with_no_draw(1))
                    .expect("the diversity of a source of many variants")
            })
        };

        let on_one = in_a_pool(1);
        assert_eq!(on_one.num_vars(0), Some(1600));
        assert_eq!(on_one.num_vars(1), Some(1600));
        for threads in [2, 4, 8] {
            assert_the_same_numbers(&on_one, &in_a_pool(threads), &format!("{threads} threads"));
        }
    }
}

#[cfg(test)]
mod the_names_of_the_statistics {
    use super::DiversityStats;
    use crate::error::Error;

    /// The name of each statistic is what a user writes in `stats` and the
    /// field of the result that holds it. The literals here are the names
    /// of the fields of `PopDiversity` in Python and in TypeScript.
    #[test]
    fn the_five_names_are_the_fields_of_a_result() {
        assert_eq!(
            DiversityStats::NAMES,
            [
                "num_alleles",
                "private_alleles",
                "variable_vars_ratio",
                "folded_sfs",
                "fis",
            ]
        );
    }

    /// Every name of the table gives back the statistic that stands beside
    /// it there, so a statistic added to the table is read by `of_name`
    /// with nothing else written for it, and none of them can be left
    /// behind as a name that `of_name` refuses.
    #[test]
    fn every_name_of_the_table_gives_back_the_statistic_beside_it() {
        for (name, stat) in DiversityStats::NAMES_AND_STATS {
            assert_eq!(
                DiversityStats::of_name(name).expect("a name of the table"),
                stat,
                "the statistic of `{name}`"
            );
        }
    }

    /// The names of the table are the names of `NAMES`, in that order, so
    /// the list a message prints is the list a user chooses from.
    #[test]
    fn the_names_of_the_table_are_the_names_of_the_result() {
        let of_the_table: Vec<&str> = DiversityStats::NAMES_AND_STATS
            .iter()
            .map(|(name, _)| *name)
            .collect();

        assert_eq!(of_the_table, DiversityStats::NAMES);
    }

    /// The five names together are every statistic, so a caller that names
    /// them all asks for what `ALL` asks for.
    #[test]
    fn the_five_names_together_are_every_statistic() {
        let mut asked_for = DiversityStats::empty();
        for name in DiversityStats::NAMES {
            asked_for |= DiversityStats::of_name(name).expect("a name of NAMES");
        }

        assert_eq!(asked_for, DiversityStats::ALL);
    }

    /// A name of no statistic is refused here and not in each binding
    /// crate, with the name the user wrote and the five they may write.
    #[test]
    fn a_name_of_no_statistic_is_refused_with_the_five_names() {
        for unknown in ["poly_vars_ratio", "NUM_ALLELES", "allelic_richness", ""] {
            let error = DiversityStats::of_name(unknown).expect_err("a name of no statistic");

            assert!(
                matches!(&error, Error::DiversityStatOfAnUnknownName { name } if name == unknown),
                "{error}"
            );
            let message = error.to_string();
            for known in DiversityStats::NAMES {
                assert!(
                    message.contains(known),
                    "the message of `{unknown}` is {message}, and it names {known}"
                );
            }
        }
    }
}
