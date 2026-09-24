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
//! The folded site frequency spectrum is the one of the five that is not a
//! number but a vector: how many of the variants in the draw for a
//! population show each count of the rarer allele in a draw of
//! `num_called_alleles` copies, the counts `j` and `num_called_alleles - j`
//! read as one because without an outgroup nothing says which allele is the
//! ancestral one. It is built from the same counts of the alleles, one
//! variant at a time.

use std::collections::HashSet;
use std::ops::Range;

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

    /// The four statistics that need no draw, built from their four
    /// constants: every one but the folded spectrum, whose bins are the
    /// counts of the rarer allele in a draw of `num_called_alleles`.
    ///
    /// It is what a Python or a TypeScript user who names no statistic asks
    /// for, so that a call with the variants and nothing else runs.
    /// [`DiversityStats::ALL`] is these four and the spectrum, so a
    /// statistic added to the module belongs to one of the two sets, and the
    /// test of this one fails when it is in neither.
    pub const WITHOUT_A_DRAW: DiversityStats = DiversityStats::NUM_ALLELES
        .union(DiversityStats::PRIVATE_ALLELES)
        .union(DiversityStats::VARIABLE_VARS_RATIO)
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

    /// The name of each statistic of this set, in the order of
    /// [`DiversityStats::NAMES_AND_STATS`].
    ///
    /// A binding crate gives its package the names of
    /// [`DiversityStats::WITHOUT_A_DRAW`] with this, which is the `stats` a
    /// user who names none asks for: the four are named here and not in
    /// each package, so a statistic that needs no draw is added to one set
    /// and reaches both languages.
    #[must_use]
    pub fn names(self) -> Vec<&'static str> {
        DiversityStats::NAMES_AND_STATS
            .iter()
            .filter(|(_, stat)| self.contains(*stat))
            .map(|(name, _)| *name)
            .collect()
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
    /// Of the variants that counted for it, the ones it called at least
    /// `num_called_alleles` alleles at: the variants in the draw for the
    /// population, which the two sums below are over and which are their
    /// divisor. It is 0 for a pass that was given no draw.
    num_vars_in_draw: u64,
    /// The alleles a draw of `num_called_alleles` is expected to show, added
    /// up over those variants.
    sum_alleles_in_draw: f64,
    /// The chance that such a draw shows more than one allele, added up over
    /// the same variants.
    sum_varies_in_draw: f64,
    /// The alleles such a draw is expected to show in the population and no
    /// draw of the same size to show in any other population, added up over
    /// the variants in the draw for every population. Its divisor is
    /// `num_vars_every_pop_in_draw` and not the `num_vars_in_draw` of the two
    /// sums above, as the divisor of `private_alleles` is
    /// `num_vars_every_pop` and not `num_vars`.
    sum_private_alleles_in_draw: f64,
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
            num_vars_in_draw: 0,
            sum_alleles_in_draw: 0.0,
            sum_varies_in_draw: 0.0,
            sum_private_alleles_in_draw: 0.0,
            sum_obs_het: 0.0,
            sum_unbiased_exp_het: 0.0,
            num_vars_with_both_hets: 0,
        }
    }

    /// The alleles a draw of `num_called_alleles` is expected to show,
    /// averaged over the variants in the draw for the population, and NaN
    /// when none is.
    fn num_alleles_in_draw(&self) -> f64 {
        let Some(num_vars) = self.num_vars_of_the_draw() else {
            return f64::NAN;
        };
        self.sum_alleles_in_draw / num_vars
    }

    /// The chance that such a draw shows more than one allele, averaged over
    /// the same variants, and NaN when none is.
    fn variable_vars_ratio_in_draw(&self) -> f64 {
        let Some(num_vars) = self.num_vars_of_the_draw() else {
            return f64::NAN;
        };
        self.sum_varies_in_draw / num_vars
    }

    /// The alleles a draw of `num_called_alleles` is expected to show in the
    /// population and no draw of the same size to show in any other,
    /// averaged over the variants in the draw for every population, and NaN
    /// when none is.
    ///
    /// `num_vars_every_pop_in_draw` is that count, which the pass keeps once
    /// for the whole call and not for each population, so this mean is the
    /// one value of a population whose divisor comes from outside it.
    fn private_alleles_in_draw(&self, num_vars_every_pop_in_draw: u64) -> f64 {
        if num_vars_every_pop_in_draw == 0 {
            return f64::NAN;
        }
        // A count below 2^53 is exact in a float64, and a pass of that many
        // variants reads more rows than any source holds.
        self.sum_private_alleles_in_draw / num_vars_every_pop_in_draw as f64
    }

    /// The divisor of the two sums of the draw, and `None` when no variant
    /// is in the draw for the population and neither sum is over anything.
    fn num_vars_of_the_draw(&self) -> Option<f64> {
        if self.num_vars_in_draw == 0 {
            return None;
        }
        // A count below 2^53 is exact in a float64, and a pass of that many
        // variants reads more rows than any source holds.
        Some(self.num_vars_in_draw as f64)
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
    /// The bins of the folded spectrum of every population, `num_sfs_bins` of
    /// them for each population in the order of the call, and no entry for a
    /// pass that was not asked for the spectrum.
    folded_sfs: Vec<f64>,
    /// How many bins the spectrum of one population has,
    /// `num_called_alleles / 2 + 1`, and 0 for a pass that was not asked for
    /// the spectrum.
    num_sfs_bins: usize,
    /// The variants that counted for every population at once.
    num_vars_every_pop: u64,
    /// Of those, the ones in the draw for every population at once.
    num_vars_every_pop_in_draw: u64,
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

    /// Of the variants that counted for the population, the ones it called
    /// at least `num_called_alleles` alleles at: the variants in the draw
    /// for it, which both of its standardized values are over. `None` when
    /// `pop` is not a population of the call.
    ///
    /// A variant is in the draw for the population when it counts for the
    /// population and the population reached that many called alleles there,
    /// both and not the second alone, so a variant the population has too
    /// little called at is out of this count however many copies it called.
    ///
    /// It is 0 for a pass that was given no draw and for a draw above every
    /// called allele of the population, and 0 is what says why a
    /// standardized value beside it is NaN.
    #[must_use]
    pub fn num_vars_in_draw(&self, pop: usize) -> Option<u64> {
        self.pops.get(pop).map(|pop| pop.num_vars_in_draw)
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

    /// Of those, the ones in the draw for every population: the variants
    /// each population reached `num_called_alleles` called alleles at. It is
    /// the divisor of the standardized private alleles, as
    /// [`PopDiversity::num_vars_every_pop`] is of the private alleles
    /// themselves, and it is 0 for a pass that was given no draw.
    #[must_use]
    pub fn num_vars_every_pop_in_draw(&self) -> u64 {
        self.num_vars_every_pop_in_draw
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

    /// The alleles a draw of `num_called_alleles` of the called alleles of
    /// the population is expected to show at a variant, averaged over the
    /// variants in the draw for it: the allelic richness of the population
    /// at a number of called alleles every population of the call is brought
    /// down to, which is what makes two populations comparable when one
    /// holds more individuals than the other and finds more alleles for that
    /// reason alone. `None` when `pop` is not a population of the call or
    /// [`DiversityStats::NUM_ALLELES`] was not asked for.
    ///
    /// It is NaN when no variant is in the draw for the population, which
    /// happens in three ways: the pass was given no draw, the draw was above
    /// every called allele of the population, and no variant counted for the
    /// population at all. [`PopDiversity::num_vars_in_draw`] is 0 in each of
    /// them and says which question was not answered.
    #[must_use]
    pub fn num_alleles_in_draw(&self, pop: usize) -> Option<f64> {
        if !self.stats.contains(DiversityStats::NUM_ALLELES) {
            return None;
        }
        self.pops.get(pop).map(OfAPop::num_alleles_in_draw)
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

    /// The alleles a draw of `num_called_alleles` of the called alleles of
    /// the population is expected to show at a variant and no draw of the
    /// same size to show in any other population of the call, averaged over
    /// the variants in the draw for every population: the private alleles of
    /// the population at a number of called alleles every population is
    /// brought down to, which is what makes two populations comparable when
    /// one holds more individuals than the other and finds alleles of its own
    /// for that reason alone. `None` when `pop` is not a population of the
    /// call or [`DiversityStats::PRIVATE_ALLELES`] was not asked for.
    ///
    /// Its divisor is [`PopDiversity::num_vars_every_pop_in_draw`] and not
    /// the [`PopDiversity::num_vars_in_draw`] of the other two standardized
    /// values, for the reason [`PopDiversity::private_alleles`] gives: a
    /// variant one population is short at is out of the private alleles of
    /// every population, and one population short of the draw takes the
    /// variant from every population's standardized value.
    ///
    /// It is the estimator of Kalinowski (2004), which ADZE computes, and it
    /// reads the draws of two populations as draws of copies of their own.
    /// Where two populations hold an individual in common that is not so, and
    /// the value can then be above the alleles their draws can really tell
    /// apart: a population against a copy of itself, at a variant where it
    /// called one allele 3 times and another once, gets 0.25 at a draw of 2
    /// and 0 only at a draw of 4, which takes every copy it called. "The
    /// cases" and "How it is verified" of "The private alleles" of
    /// `docs/specs/diversity.md` work that out. popnei computes the estimator
    /// and does not refuse the overlap, since a caller may put an individual
    /// in more than one population and every other value here reads such a
    /// population without trouble.
    ///
    /// It is NaN when no variant is in the draw for every population, which
    /// happens when the pass was given no draw, when the draw is above the
    /// called alleles of one population at every variant, and when no variant
    /// counted for every population. [`PopDiversity::num_vars_every_pop_in_draw`]
    /// is 0 in each of them.
    #[must_use]
    pub fn private_alleles_in_draw(&self, pop: usize) -> Option<f64> {
        if !self.stats.contains(DiversityStats::PRIVATE_ALLELES) {
            return None;
        }
        let num_vars_every_pop_in_draw = self.num_vars_every_pop_in_draw;
        self.pops
            .get(pop)
            .map(|pop| pop.private_alleles_in_draw(num_vars_every_pop_in_draw))
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

    /// The chance that a draw of `num_called_alleles` of the called alleles
    /// of the population shows more than one allele, averaged over the
    /// variants in the draw for it: the ratio of variable variants of the
    /// population at that common number of called alleles, which is what
    /// makes it comparable between two populations of different size.
    /// `None` when `pop` is not a population of the call or
    /// [`DiversityStats::VARIABLE_VARS_RATIO`] was not asked for.
    ///
    /// It is NaN in the three cases of
    /// [`PopDiversity::num_alleles_in_draw`]. On a dataset whose every
    /// variant has two alleles it is that value minus 1, a draw there
    /// showing one allele or two; a variant of more alleles has no such
    /// identity.
    #[must_use]
    pub fn variable_vars_ratio_in_draw(&self, pop: usize) -> Option<f64> {
        if !self.stats.contains(DiversityStats::VARIABLE_VARS_RATIO) {
            return None;
        }
        self.pops.get(pop).map(OfAPop::variable_vars_ratio_in_draw)
    }

    /// How the variants in the draw for the population are spread over the
    /// count of their rarer allele: one value for each count from 0 to
    /// `num_called_alleles / 2`, the value at `j` being how many of those
    /// variants a draw of `num_called_alleles` of the copies the population
    /// called is expected to show `j` copies of the rarer allele at. It is the
    /// shape a population's history leaves in its variants and the input of
    /// every program that fits a demographic model. `None` when `pop` is not a
    /// population of the call or [`DiversityStats::FOLDED_SFS`] was not asked
    /// for.
    ///
    /// It is folded: the counts `j` and `num_called_alleles - j` are one
    /// value, because without an outgroup nothing says which allele is the
    /// ancestral one, only which of the two is the rarer. The value at 0 holds
    /// the variants a draw shows one allele at, which a draw can do at a
    /// variant that varies in the population, and at an even
    /// `num_called_alleles` the last value is not doubled, `j` and
    /// `num_called_alleles - j` being the same count there.
    ///
    /// The values are not whole numbers: each variant in the draw gives every
    /// count the chance that a draw shows that many rarer copies there, so the
    /// values sum to [`PopDiversity::num_vars_in_draw`]. A variant of more
    /// than two alleles counts as the major allele of the population there,
    /// the one it called most often and the lower numbered of two it called
    /// equally often, against every other allele of the variant summed into
    /// one rarer allele.
    ///
    /// Every value is 0 when no variant is in the draw for the population,
    /// which happens in the three ways [`PopDiversity::num_alleles_in_draw`]
    /// lists and which that count of 0 says.
    #[must_use]
    pub fn folded_sfs(&self, pop: usize) -> Option<&[f64]> {
        if !self.stats.contains(DiversityStats::FOLDED_SFS) {
            return None;
        }
        the_bins_of_the_pop(pop, self.num_sfs_bins).and_then(|bins| self.folded_sfs.get(bins))
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

    /// The alleles a draw is expected to show and the chance that it varies,
    /// each added up over the variants in the draw for one population, and
    /// how many those variants are.
    ///
    /// A test compares these bit for bit for the reason
    /// [`PopDiversity::the_sums_behind_the_fis`] gives: the division that
    /// makes a standardized value out of a sum absorbs a difference of the
    /// last bits of that sum.
    #[cfg(test)]
    fn the_sums_of_the_draw(&self, pop: usize) -> Option<(f64, f64, u64)> {
        self.pops.get(pop).map(|pop| {
            (
                pop.sum_alleles_in_draw,
                pop.sum_varies_in_draw,
                pop.num_vars_in_draw,
            )
        })
    }

    /// The alleles a draw is expected to show in one population and in no
    /// other, added up over the variants in the draw for every population.
    ///
    /// A test compares it bit for bit for the reason
    /// [`PopDiversity::the_sums_behind_the_fis`] gives, and it is apart from
    /// [`PopDiversity::the_sums_of_the_draw`] because its divisor is the
    /// variants in the draw for every population, which the result holds once
    /// for the whole call.
    #[cfg(test)]
    fn the_sum_of_the_private_alleles_of_the_draw(&self, pop: usize) -> Option<f64> {
        self.pops
            .get(pop)
            .map(|pop| pop.sum_private_alleles_in_draw)
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

/// The draw of a common number of called alleles, which a pass given
/// `num_called_alleles` carries and a pass given none does not.
///
/// A variant is in the draw for a population when it counts for that
/// population and the population called at least `num_called_alleles`
/// alleles there, both and not the second alone, which "Its Python function"
/// of `docs/specs/diversity.md` states. Each of the two standardized values
/// is the mean over the variants in the draw for the population, so what one
/// population keeps of the draw is one count of variants and one sum for
/// each value.
#[derive(Debug, Clone, Copy)]
struct OfTheDraw {
    /// The called alleles every population is brought down to, the `g` of
    /// the formulas of `docs/specs/diversity.md`.
    num_called_alleles: u32,
    /// The statistics the pass was asked for, which say which of the two
    /// sums a variant in the draw adds to. Each of them costs a product of
    /// `num_called_alleles` factors for every allele the population called,
    /// so the one that was not asked for is not computed.
    stats: DiversityStats,
    /// How many bins the folded spectrum of one population has,
    /// `num_called_alleles / 2 + 1`, and 0 for a pass that was not asked for
    /// the spectrum, which keeps no bin.
    num_bins: usize,
}

impl OfTheDraw {
    /// The draw of a pass that was given one, and `None` for a pass that was
    /// given none, whose counts of the variants of a draw stay 0 and whose
    /// standardized values are NaN.
    fn of(options: &DiversityOptions) -> Option<OfTheDraw> {
        options
            .num_called_alleles
            .map(|num_called_alleles| OfTheDraw {
                num_called_alleles,
                stats: options.stats,
                num_bins: the_num_of_bins(options.stats, num_called_alleles),
            })
    }

    /// It counts one variant in the draw for a population and adds what a
    /// draw is expected to show there to the sums the pass was asked for.
    ///
    /// `counts` is how often the population called each allele of the
    /// variant, `one_past_the_largest` the entry of `counts` a walk over them
    /// stops at, and `called_alleles` their sum, which the caller has found
    /// to be at least `num_called_alleles`. `of_the_pops_bins` is the bins of
    /// the folded spectrum of that population, `num_bins` of them, and an
    /// empty slice for a pass that was not asked for the spectrum.
    fn add_the_var(
        &self,
        counts: &AlleleCounts,
        one_past_the_largest: usize,
        called_alleles: u32,
        counted: &mut OfAPop,
        of_the_pops_bins: &mut [f64],
    ) {
        // One variant of the draw, and a pass of more than
        // 18446744073709551615 variants reads more rows than any source
        // holds.
        counted.num_vars_in_draw = counted.num_vars_in_draw.saturating_add(1);
        if self.stats.contains(DiversityStats::NUM_ALLELES) {
            counted.sum_alleles_in_draw +=
                self.alleles_a_draw_shows(counts, one_past_the_largest, called_alleles);
        }
        if self.stats.contains(DiversityStats::VARIABLE_VARS_RATIO) {
            counted.sum_varies_in_draw +=
                self.chance_a_draw_varies(counts, one_past_the_largest, called_alleles);
        }
        if self.stats.contains(DiversityStats::FOLDED_SFS) {
            self.add_the_bins_of_the_var(
                counts,
                one_past_the_largest,
                called_alleles,
                of_the_pops_bins,
            );
        }
    }

    /// It adds to the bins of one population the chance that a draw of
    /// `num_called_alleles` of the `called_alleles` copies it called at one
    /// variant shows each count of the rarer allele: the projection of "What
    /// it gives" of "The folded site frequency spectrum" of
    /// `docs/specs/diversity.md`.
    ///
    /// The variant is read as if it had two alleles, the major allele of the
    /// population there against every other allele of it summed into one
    /// rarer allele. So with `c` the copies the population called, `m` the
    /// ones that are not of its major allele and `g` the draw, the chance of
    /// `j` rarer copies is `C(m, j) * C(c - m, g - j) / C(c, g)` and it goes
    /// to the bin `min(j, g - j)`, the counts `j` and `g - j` being one bin.
    ///
    /// `j` runs from `max(0, g - (c - m))` to `min(m, g)` and over no count
    /// outside that range: the draw holds neither more rarer copies than the
    /// variant has nor fewer than the copies of the major allele leave room
    /// for. `m` can be far above `g`, 48 against 20 on the panel of that
    /// spec, and a `j` above `g` has a chance of 0 and falls in no bin of the
    /// spectrum, `g - j` being below 0 there.
    ///
    /// The chance of the first count of the range is
    /// [`chance_a_draw_misses_an_allele`], the one product of `g` factors this
    /// module has, read in one of its two ways. Where the range starts at 0 it
    /// is `C(c - m, g) / C(c, g)`, the chance that the draw misses an allele
    /// the population called `m` times. Where it starts above 0 every copy of
    /// the major allele is in the draw, and the chance of that is the chance
    /// that the `c - g` copies the draw leaves behind, which are rarer ones
    /// alone, miss every copy of the major allele.
    ///
    /// The chance of each count above the first is the one before it times
    /// `(m - j) (g - j) / ((j + 1) (c - m - g + j + 1))`, so every value the
    /// loop holds is a chance and stays between 0 and 1, where the three
    /// binomial coefficients on their own are above what an `f64` holds: a
    /// draw of 180 of 20000 copies asks for a number of 444 digits and 171!
    /// is already an infinity. The counts are stepped through with `j` rising
    /// on every call, so two populations of one pass round the same way.
    fn add_the_bins_of_the_var(
        &self,
        counts: &AlleleCounts,
        one_past_the_largest: usize,
        called_alleles: u32,
        of_the_pops_bins: &mut [f64],
    ) {
        let of_the_major_allele = the_copies_of_the_major_allele(counts, one_past_the_largest);
        // The copies of every allele that is not the major one. The count of
        // one allele is one part of the sum the called alleles are, so it is
        // never the larger of the two and the subtraction never saturates.
        let of_the_rarer_allele = called_alleles.saturating_sub(of_the_major_allele);
        let drawn = self.num_called_alleles;
        // The counts of the rarer allele the draw can show: at least what the
        // copies of the major allele leave room for, which is 0 where they
        // could fill the draw on their own, and at most every rarer copy of
        // the variant and no more than the draw itself.
        let first = drawn.saturating_sub(of_the_major_allele);
        let last = of_the_rarer_allele.min(drawn);
        let mut chance = if first == 0 {
            chance_a_draw_misses_an_allele(called_alleles, of_the_rarer_allele, drawn)
        } else {
            // The draw takes every copy of the major allele, which is the
            // chance that a draw of the copies it leaves behind misses each of
            // them. A draw of every copy leaves none behind and takes them
            // with chance 1, which is the product of no factor.
            chance_a_draw_misses_an_allele(
                called_alleles,
                of_the_major_allele,
                called_alleles.saturating_sub(drawn),
            )
        };
        for count in first..=last {
            // The bin of the count, which is at most `g / 2` and so is a bin
            // of the spectrum: the counts above half the draw fold onto the
            // ones below it.
            let bin = count.min(drawn.saturating_sub(count));
            if let Some(of_the_bin) = usize::try_from(bin)
                .ok()
                .and_then(|bin| of_the_pops_bins.get_mut(bin))
            {
                *of_the_bin += chance;
            }
            // The chance of one more copy of the rarer allele, from the one in
            // hand. Past the last count of the range it is 0 or a number
            // nothing reads, one of its two numerators being 0 there.
            let of_the_rarer_allele_left = f64::from(of_the_rarer_allele) - f64::from(count);
            let of_the_draw_left = f64::from(drawn) - f64::from(count);
            let room_the_major_allele_leaves =
                f64::from(of_the_major_allele) - f64::from(drawn) + f64::from(count) + 1.0;
            chance *= (of_the_rarer_allele_left * of_the_draw_left)
                / ((f64::from(count) + 1.0) * room_the_major_allele_leaves);
        }
    }

    /// The alleles a draw of `num_called_alleles` of the `called_alleles`
    /// copies the population called at one variant is expected to show: the
    /// chance that the draw holds a copy of an allele, summed over the
    /// alleles the population called there, which is the `E` of "What it
    /// gives" of "The number of alleles" of `docs/specs/diversity.md`.
    ///
    /// The alleles are summed in the order of their numbers on every call,
    /// so two populations of one pass add the same terms the same way.
    fn alleles_a_draw_shows(
        &self,
        counts: &AlleleCounts,
        one_past_the_largest: usize,
        called_alleles: u32,
    ) -> f64 {
        the_alleles_called(counts, one_past_the_largest)
            .map(|count| {
                1.0 - chance_a_draw_misses_an_allele(called_alleles, count, self.num_called_alleles)
            })
            .sum()
    }

    /// The chance that such a draw shows more than one allele: one minus the
    /// chance that every copy it takes is of the same allele, summed over the
    /// alleles the population called, which is the `P` of "What it gives" of
    /// "The variable variants" of `docs/specs/diversity.md`. A draw is all of
    /// one allele or of more than one, so the two are one minus each other.
    fn chance_a_draw_varies(
        &self,
        counts: &AlleleCounts,
        one_past_the_largest: usize,
        called_alleles: u32,
    ) -> f64 {
        let all_of_one_allele: f64 = the_alleles_called(counts, one_past_the_largest)
            .map(|count| {
                chance_a_draw_is_all_of_one_allele(called_alleles, count, self.num_called_alleles)
            })
            .sum();
        1.0 - all_of_one_allele
    }
}

/// How often the population called each of the alleles it called at one
/// variant, in the order of their numbers.
///
/// `counts` is what the counts of the variant left for the population and
/// `one_past_the_largest` the bound they gave on the alleles they wrote. The
/// alleles the population did not call are left out: each of them would add
/// a term of 0 to a sum over a draw and cost a product of
/// `num_called_alleles` factors, and a population counted by reading the row
/// as it is has no bound on the alleles of the row, so all 128 entries of
/// its counts are walked.
fn the_alleles_called(
    counts: &AlleleCounts,
    one_past_the_largest: usize,
) -> impl Iterator<Item = u32> {
    counts
        .iter()
        .take(one_past_the_largest)
        .copied()
        .filter(|count| *count > 0)
}

/// How often a population called the allele it called most often at one
/// variant: the copies of its major allele there, which the folded spectrum
/// reads every other allele of the variant against, and 0 where it called
/// nothing.
///
/// `counts` is what the counts of the variant left for the population and
/// `one_past_the_largest` the bound they gave on the alleles they wrote.
/// Which of two alleles a population called equally often is the major one
/// does not change how many copies that allele has, so the rule of
/// `docs/glossary.md` that the lower numbered of two such alleles is the
/// major one needs nothing here.
fn the_copies_of_the_major_allele(counts: &AlleleCounts, one_past_the_largest: usize) -> u32 {
    the_alleles_called(counts, one_past_the_largest)
        .max()
        .unwrap_or(0)
}

/// How many bins the folded spectrum of a draw of `num_called_alleles`
/// copies has, and 0 for a pass that was not asked for the spectrum, which
/// keeps no bin.
///
/// A draw of `g` copies shows 0 to `g` copies of the rarer allele, and the
/// counts `j` and `g - j` are one bin, so the bins are the counts from 0 to
/// `g / 2` and they are `g / 2 + 1`.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the divisor is the literal 2, and half of the largest draw a u32 holds is 2147483647, so the sum is below the largest usize of every target popnei builds for, wasm's 32 bit one among them"
)]
fn the_num_of_bins(stats: DiversityStats, num_called_alleles: u32) -> usize {
    if !stats.contains(DiversityStats::FOLDED_SFS) {
        return 0;
    }
    // A `usize` is 32 bits in wasm and 64 natively, and a draw of more copies
    // than either holds is a draw of more than any dataset has.
    let drawn = usize::try_from(num_called_alleles).unwrap_or(usize::MAX);
    drawn / 2 + 1
}

/// Where the bins of the folded spectrum of one population are among those
/// of every population, which hold `num_bins` of them for each population,
/// in the order of the populations.
///
/// It is `None` for a first bin or a last one above what a `usize` counts,
/// which asks for more bins than the machine could hold in any case.
fn the_bins_of_the_pop(pop: usize, num_bins: usize) -> Option<Range<usize>> {
    let first = pop.checked_mul(num_bins)?;
    let past_the_last = first.checked_add(num_bins)?;
    Some(first..past_the_last)
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
    /// The draw of a common number of called alleles, which a pass given no
    /// `num_called_alleles` does not carry and counts no variant of.
    of_the_draw: Option<OfTheDraw>,
    /// How many alleles a population has to have called at a variant for
    /// the variant to count for it: `min_num_individuals` genotypes of the
    /// ploidy.
    min_called_alleles: u64,
    /// The alleles of one genotype, which the reader says its source has.
    ploidy: usize,
}

impl OfThePass<'_> {
    /// How many bins the folded spectrum of one population of the pass has,
    /// and 0 for a pass that was not asked for the spectrum and for one that
    /// was given no draw, neither of which keeps a bin.
    fn num_sfs_bins(&self) -> usize {
        self.of_the_draw
            .map_or(0, |of_the_draw| of_the_draw.num_bins)
    }
}

/// What a pass, or one chunk of the rows of a block, has counted over every
/// population.
#[derive(Debug, Clone)]
struct Totals {
    pops: Vec<OfAPop>,
    /// The bins of the folded spectrum of every population, `num_bins` of
    /// them for each population in the order of the populations, and no entry
    /// for a pass that was not asked for the spectrum.
    ///
    /// They are one vector of every population's bins and not a vector inside
    /// each [`OfAPop`], which keeps the counts of a population `Copy`: a chunk
    /// of 64 rows counts into a `Totals` of its own, so a vector for each
    /// population would be one allocation for each population of each chunk,
    /// 7850 of them for a block of 10000 rows and 50 populations, where these
    /// are 157, one for each chunk.
    folded_sfs: Vec<f64>,
    num_vars_every_pop: u64,
    num_vars_every_pop_in_draw: u64,
}

impl Totals {
    /// The counts of `num_pops` populations before any variant is read, with
    /// `num_bins` bins of the folded spectrum for each of them.
    fn of(num_pops: usize, num_bins: usize) -> Totals {
        Totals {
            pops: vec![OfAPop::none(); num_pops],
            // A product that saturated would ask for more bins than the
            // machine has memory for, and the allocation of them is what
            // fails.
            folded_sfs: vec![0.0; num_pops.saturating_mul(num_bins)],
            num_vars_every_pop: 0,
            num_vars_every_pop_in_draw: 0,
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
            of_the_pass.num_vars_in_draw = of_the_pass
                .num_vars_in_draw
                .saturating_add(of_the_chunk.num_vars_in_draw);
            // The chunks of a block are added in the order of the block and
            // the blocks in the order of the pass, so these five sums of
            // float64 do not depend on how many threads read the rows.
            of_the_pass.sum_obs_het += of_the_chunk.sum_obs_het;
            of_the_pass.sum_unbiased_exp_het += of_the_chunk.sum_unbiased_exp_het;
            of_the_pass.sum_alleles_in_draw += of_the_chunk.sum_alleles_in_draw;
            of_the_pass.sum_varies_in_draw += of_the_chunk.sum_varies_in_draw;
            of_the_pass.sum_private_alleles_in_draw += of_the_chunk.sum_private_alleles_in_draw;
            of_the_pass.num_vars_with_both_hets = of_the_pass
                .num_vars_with_both_hets
                .saturating_add(of_the_chunk.num_vars_with_both_hets);
        }
        // The bins of every population, added in the order of the populations
        // and of the counts of the rarer allele, as the five sums above are
        // added in the order of the chunks: no bin depends on how many threads
        // read the rows either.
        for (of_the_pass, of_the_chunk) in self
            .folded_sfs
            .iter_mut()
            .zip(of_the_chunk.folded_sfs.iter())
        {
            *of_the_pass += *of_the_chunk;
        }
        self.num_vars_every_pop = self
            .num_vars_every_pop
            .saturating_add(of_the_chunk.num_vars_every_pop);
        self.num_vars_every_pop_in_draw = self
            .num_vars_every_pop_in_draw
            .saturating_add(of_the_chunk.num_vars_every_pop_in_draw);
    }

    /// It empties every count, so that one chunk of rows after another is
    /// counted in the same `Totals`.
    fn forget_what_it_holds(&mut self) {
        for of_the_pop in &mut self.pops {
            *of_the_pop = OfAPop::none();
        }
        self.folded_sfs.fill(0.0);
        self.num_vars_every_pop = 0;
        self.num_vars_every_pop_in_draw = 0;
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
/// for with no `num_called_alleles`, a `num_called_alleles` below 2 or above
/// the individuals of the reader times its ploidy, a population with no
/// individual, an index that is not an individual of the dataset, an
/// individual asked for more than once, no variant in the reader, a variant
/// of more alleles than a count of them holds, and those of the reader,
/// among them a ploidy of 0 or above the 255 a genotype of popnei holds.
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
    let num_individuals = reader.individuals().len();
    let ploidy = reader.ploidy();
    // The ploidy the reader states, which the threshold of called
    // genotypes is measured in. A reader that states one of 0, or one above
    // the 255 a genotype of popnei holds, is refused here: a threshold
    // built from a ploidy that did not fit would be one no population ever
    // meets, and every population would silently count no variant. It is
    // checked before the draw because the largest draw the dataset allows is
    // the individuals times this ploidy, and a ploidy of 0 would make that
    // largest draw 0 and refuse every draw in its name.
    let ploidy_of_the_gts = checked_ploidy("ploidy", ploidy)?;
    check_the_draw(options, num_individuals, ploidy_of_the_gts)?;
    let of_the_pops = pops_of_the_pass(pops, num_individuals)?;
    // The five statistics follow from the genotypes of a row, so no column
    // of a block is read and the reader is asked to fill none of them.
    reader.set_needs(Needs::GTS);
    let of_the_pass = OfThePass {
        pops: &of_the_pops,
        counts_the_alleles: CountsTheAlleles::of(options.stats),
        counts_the_private_alleles: CountsThePrivateAlleles::of(options.stats),
        heterozygosities: Heterozygosities::of(options.stats, ploidy)?,
        of_the_draw: OfTheDraw::of(options),
        min_called_alleles: min_called_alleles(options.min_num_individuals, ploidy_of_the_gts),
        ploidy,
    };
    let num_sfs_bins = of_the_pass.num_sfs_bins();
    let mut totals = Totals::of(of_the_pops.len(), num_sfs_bins);
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
        folded_sfs: totals.folded_sfs,
        num_sfs_bins,
        num_vars_every_pop: totals.num_vars_every_pop,
        num_vars_every_pop_in_draw: totals.num_vars_every_pop_in_draw,
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
/// `num_individuals` and `ploidy` are the dataset's, and their product is
/// every gene copy it holds: the most alleles any population can have called
/// at any variant, and so the largest draw the dataset allows. A draw at or
/// below it that this dataset's missing genotypes leave no population able
/// to fill is not refused here: that draw gives NaN in the standardized
/// values and a spectrum of zeros, which "The cases" of
/// `docs/specs/diversity.md` asks for.
///
/// # Errors
///
/// The folded spectrum asked for with no `num_called_alleles`, whose bins
/// are the counts of the rarer allele in a draw of that many, a
/// `num_called_alleles` below 2, and one above the individuals of the
/// dataset times the ploidy.
fn check_the_draw(options: &DiversityOptions, num_individuals: usize, ploidy: u32) -> Result<()> {
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
            // The gene copies of the dataset are counted in a `u64`: the
            // individuals times the ploidy passes what a `u32` holds above
            // 2147483647 of them, and passes a `usize` in WebAssembly, where
            // it is 32 bits. A dataset whose individuals, or whose product,
            // do not fit in a `u64` allows every draw a `u32` holds, so
            // saturating at the largest `u64` refuses none of them.
            let largest_draw = u64::try_from(num_individuals)
                .unwrap_or(u64::MAX)
                .saturating_mul(u64::from(ploidy));
            if u64::from(num_called_alleles) > largest_draw {
                return Err(Error::DiversityDrawLargerThanTheDataset {
                    num_called_alleles,
                    largest_draw,
                    num_individuals,
                    ploidy,
                });
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
    let num_bins = of_the_pass.num_sfs_bins();
    let of_the_chunks: Result<Vec<Totals>> = block
        .gts
        .par_chunks(alleles_of_a_chunk(alleles_per_var))
        .map(|chunk| {
            let mut of_the_chunk = Totals::of(num_pops, num_bins);
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
            let mut read_again = Totals::of(num_pops, num_bins);
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
    let mut of_the_chunk = Totals::of(of_the_pass.pops.len(), of_the_pass.num_sfs_bins());
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
    // The room the standardized private alleles work in, one chance and one
    // sum for each population, written again for every row as the counts above
    // are.
    let mut of_the_draw_at_the_row = OfTheDrawAtTheRow::of(of_the_pass.pops.len());
    for row in gts.chunks_exact(alleles_per_var) {
        let mut every_pop = true;
        // A pass that was given no draw has no variant in the draw for every
        // population, so its second count of them stays 0.
        let mut every_pop_in_draw = of_the_pass.of_the_draw.is_some();
        for (pop, ((of_the_pop, counted), at_the_row)) in of_the_pass
            .pops
            .iter()
            .zip(totals.pops.iter_mut())
            .zip(of_each_pop.iter_mut())
            .enumerate()
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
            at_the_row.called_alleles = called_alleles;
            // A population that called nothing at the variant does not
            // count it whatever `min_num_individuals` is, so a threshold of
            // 0 does not put a variant with no data into the totals.
            if called_alleles == 0 || u64::from(called_alleles) < of_the_pass.min_called_alleles {
                every_pop = false;
                every_pop_in_draw = false;
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
            // The variant is in the draw for the population when it counts
            // for the population, which the lines above have settled, and
            // the population called at least `num_called_alleles` alleles
            // there.
            if let Some(of_the_draw) = of_the_pass.of_the_draw {
                if called_alleles < of_the_draw.num_called_alleles {
                    every_pop_in_draw = false;
                } else {
                    // The bins of this population among those of every
                    // population, which the counts of the chunk hold in one
                    // vector. The slice is empty for a pass that was not asked
                    // for the spectrum and for no other: the counts of a chunk
                    // hold `num_bins` bins for each of the populations of the
                    // pass, so every population of the row has its own.
                    let of_the_pops_bins = the_bins_of_the_pop(pop, of_the_draw.num_bins)
                        .and_then(|bins| totals.folded_sfs.get_mut(bins))
                        .unwrap_or_default();
                    of_the_draw.add_the_var(
                        &at_the_row.counts,
                        one_past_the_largest,
                        called_alleles,
                        counted,
                        of_the_pops_bins,
                    );
                }
            }
        }
        if every_pop {
            totals.num_vars_every_pop = totals.num_vars_every_pop.saturating_add(1);
            if every_pop_in_draw {
                totals.num_vars_every_pop_in_draw =
                    totals.num_vars_every_pop_in_draw.saturating_add(1);
            }
            if of_the_pass.counts_the_private_alleles == CountsThePrivateAlleles::Yes {
                add_the_private_alleles(&of_each_pop, &mut num_pops_that_called, &mut totals.pops);
                // The standardized value is the mean over the variants in the
                // draw for every population, so a row one population is short
                // of the draw at adds to no population's sum.
                if every_pop_in_draw && let Some(of_the_draw) = of_the_pass.of_the_draw {
                    add_the_standardized_private_alleles(
                        &of_each_pop,
                        of_the_draw.num_called_alleles,
                        &mut of_the_draw_at_the_row,
                        &mut totals.pops,
                    );
                }
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
    /// How many alleles it called there, the sum of `counts`, which is the
    /// `c` of the chance that a draw of the population misses an allele. The
    /// standardized private alleles read it of every population of the row at
    /// once, which is why it is kept here and not worked out again from the
    /// counts.
    called_alleles: u32,
}

impl OfAPopAtTheRow {
    /// The counts of one population before any row is read, which hold no
    /// allele.
    fn none() -> OfAPopAtTheRow {
        OfAPopAtTheRow {
            counts: [0; 128],
            one_past_the_largest: 0,
            called_alleles: 0,
        }
    }
}

/// The room the standardized private alleles need while one row is read: for
/// each population, the chance that its draw misses the allele in hand and
/// what it has summed over the alleles of the row so far.
///
/// It is allocated once for a chunk of rows and written again for every row,
/// as the counts of the populations beside it are, so the loop over the rows
/// allocates nothing and what a chunk holds grows with the populations and
/// not with the rows.
#[derive(Debug)]
struct OfTheDrawAtTheRow {
    /// One chance for each population, of the allele in hand: the chance that
    /// a draw of `num_called_alleles` of what the population called at the row
    /// holds no copy of it.
    chance_each_draw_misses_the_allele: Vec<f64>,
    /// One sum for each population: the alleles of the row expected to be in
    /// its own draw and in no other population's, added over the alleles it
    /// called there.
    private_alleles_of_each_pop: Vec<f64>,
}

impl OfTheDrawAtTheRow {
    /// The room for `num_pops` populations, before any row is read.
    fn of(num_pops: usize) -> OfTheDrawAtTheRow {
        OfTheDrawAtTheRow {
            chance_each_draw_misses_the_allele: vec![0.0; num_pops],
            private_alleles_of_each_pop: vec![0.0; num_pops],
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

/// It adds to every population the alleles a draw of `num_called_alleles` of
/// what it called at the row is expected to show that no draw of the same size
/// shows in any other population of the pass: the `E` of "What it gives" of
/// "The private alleles" of `docs/specs/diversity.md`, the estimator of
/// Kalinowski (2004).
///
/// It is called for a row that is in the draw for every population and for no
/// other, which is what [`PopDiversity::num_vars_every_pop_in_draw`], the
/// divisor of this value, counts: one population short of the draw takes the
/// row from the standardized private alleles of every population, as one
/// population short of data takes it from the counted ones.
///
/// `of_each_pop` is what each population called at the row, in the order of
/// `totals`, and `at_the_row` the room for one chance and one sum for each of
/// them, whose contents when this is called are of no interest.
///
/// For each allele of the row the chance that a draw misses it is worked out
/// once for each population, and the term of a population is the chance that
/// its own draw shows the allele times the chances that the draw of every
/// other population misses it. Every one of those chances is
/// [`chance_a_draw_misses_an_allele`], the one product of `num_called_alleles`
/// factors this module has. The alleles are summed in the order of their
/// numbers and the chances of the other populations multiplied in the order of
/// the call, on every row and for every population, so no value depends on how
/// the rows were shared out among the threads. In its last bits a value does
/// depend on the order the populations were given, that order being the order
/// of the product over the others.
///
/// That product is taken again for each population, which for an allele of a
/// pass of 50 populations is 2500 multiplications where a product from the left
/// and one from the right would give all 50 in 100. What stands beside those
/// multiplications is the 50 products of `num_called_alleles` divisions the
/// chances cost, 9000 of them at a draw of 180, and what the two walks would
/// cost is the rounding: a population would then multiply the chances before it
/// and the chances after it in two groups, and populations at different places
/// would round differently.
fn add_the_standardized_private_alleles(
    of_each_pop: &[OfAPopAtTheRow],
    num_called_alleles: u32,
    at_the_row: &mut OfTheDrawAtTheRow,
    totals: &mut [OfAPop],
) {
    let one_past_the_largest = of_each_pop
        .iter()
        .map(|of_the_pop| of_the_pop.one_past_the_largest)
        .max()
        .unwrap_or(0);
    at_the_row.private_alleles_of_each_pop.fill(0.0);
    for allele in 0..one_past_the_largest {
        // An allele no population called would add a term of 0 to every
        // population and a product of `num_called_alleles` factors to the work
        // of each, so the alleles between the ones the row holds are passed
        // over here.
        if !of_each_pop
            .iter()
            .any(|of_the_pop| count_of_the_allele(of_the_pop, allele) > 0)
        {
            continue;
        }
        the_chance_each_draw_misses_the_allele(
            of_each_pop,
            allele,
            num_called_alleles,
            &mut at_the_row.chance_each_draw_misses_the_allele,
        );
        let of_each_draw = &at_the_row.chance_each_draw_misses_the_allele;
        let of_each_pops_sum = at_the_row.private_alleles_of_each_pop.iter_mut();
        for (pop, (of_the_pop, private_alleles)) in
            of_each_pop.iter().zip(of_each_pops_sum).enumerate()
        {
            // The sum is over the alleles this population called: an allele it
            // did not call is in no draw of its own, so its term is 0.
            if count_of_the_allele(of_the_pop, allele) == 0 {
                continue;
            }
            let of_its_own_draw = of_each_draw.get(pop).copied().unwrap_or(1.0);
            let missed_by_every_other: f64 = of_each_draw
                .iter()
                .enumerate()
                .filter(|(other, _)| *other != pop)
                .map(|(_, chance)| *chance)
                .product();
            *private_alleles += (1.0 - of_its_own_draw) * missed_by_every_other;
        }
    }
    for (private_alleles, counted) in at_the_row
        .private_alleles_of_each_pop
        .iter()
        .zip(totals.iter_mut())
    {
        // One sum for each population, over the alleles of this row, added to
        // what it has over the rows before: the mean of the standardized
        // private alleles is the mean of one value for each variant, as the
        // other two standardized values are.
        counted.sum_private_alleles_in_draw += *private_alleles;
    }
}

/// It writes, for each population of the row, the chance that a draw of
/// `num_called_alleles` of what it called there holds no copy of `allele`.
///
/// `of_each_draw` is one entry for each population of `of_each_pop`, in its
/// order, and every one of them is written.
fn the_chance_each_draw_misses_the_allele(
    of_each_pop: &[OfAPopAtTheRow],
    allele: usize,
    num_called_alleles: u32,
    of_each_draw: &mut [f64],
) {
    for (of_the_pop, chance) in of_each_pop.iter().zip(of_each_draw.iter_mut()) {
        *chance = chance_a_draw_misses_an_allele(
            of_the_pop.called_alleles,
            count_of_the_allele(of_the_pop, allele),
            num_called_alleles,
        );
    }
}

/// How often one population called one allele of the row.
///
/// It is 0 for an allele above the largest one that population called, every
/// entry of its counts from there up holding 0, which is what lets the
/// populations of a row be read at every allele the row holds and not only at
/// the ones each of them called.
fn count_of_the_allele(of_the_pop: &OfAPopAtTheRow, allele: usize) -> u32 {
    of_the_pop.counts.get(allele).copied().unwrap_or(0)
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

/// The chance that a draw of `num_called_alleles` of the `called_alleles`
/// copies a population called at one variant, taken without replacement,
/// holds no copy of an allele the population called `count_of_the_allele`
/// times: `C(c - n, g) / C(c, g)`, with `c` the called alleles, `n` the
/// count of the allele and `g` the draw.
///
/// One minus it is the chance that the draw shows the allele, and the sum
/// of that over the alleles the population called is the alleles the draw
/// is expected to show, which "What it gives" of "The number of alleles"
/// of `docs/specs/diversity.md` states. Every standardized value of this
/// module is built from this one chance.
///
/// `num_called_alleles` is at most `called_alleles`, which is what makes
/// the variant one of the draw for that population, and the caller leaves
/// the variants where it is not out. The chance is exactly 0 when fewer
/// than `num_called_alleles` of the called copies are of another allele,
/// since then no draw of that size avoids this allele; a
/// `num_called_alleles` above `called_alleles`, which is no draw at all,
/// gives that same 0. A `count_of_the_allele` of 0, an allele the
/// population did not call, gives 1.
///
/// It is the product of `num_called_alleles` factors,
/// `(c - n - i) / (c - i)` for `i` from 0 up, each in `f64`, and not the
/// ratio of three factorials: 171! is already an infinity in an `f64`, and
/// the largest dataset of `docs/objectives.md`, 10000 individuals, has
/// 20000 copies at one variant. The factors are multiplied with `i` rising
/// on every call, so two populations of one pass round the same way.
/// Division and multiplication are rounded the same way on every platform,
/// unlike `exp` and `ln`, so a test of this asserts the digits of what it
/// gives.
fn chance_a_draw_misses_an_allele(
    called_alleles: u32,
    count_of_the_allele: u32,
    num_called_alleles: u32,
) -> f64 {
    let all_the_copies = f64::from(called_alleles);
    let of_another_allele = all_the_copies - f64::from(count_of_the_allele);
    let drawn = f64::from(num_called_alleles);

    if of_another_allele < drawn {
        return 0.0;
    }

    (0..num_called_alleles)
        .map(f64::from)
        .fold(1.0, |chance, drawn_before| {
            chance * ((of_another_allele - drawn_before) / (all_the_copies - drawn_before))
        })
}

/// The chance that every one of the `num_called_alleles` copies a draw
/// takes of the `called_alleles` a population called at one variant is the
/// allele it called `count_of_the_allele` times: `C(n, g) / C(c, g)`, with
/// `c` the called alleles, `n` the count of the allele and `g` the draw.
///
/// The chance that the draw shows more than one allele is one minus this
/// summed over the alleles the population called, which "What it gives" of
/// "The variable variants" of `docs/specs/diversity.md` states.
///
/// It is [`chance_a_draw_misses_an_allele`] with the copies of every other
/// allele in the place of the allele's own, since a draw that is all of one
/// allele is a draw that missed every other and there are `c - n` copies of
/// those. The two chances are that one product of `g` factors and not two
/// copies of one formula, so the standardized ratio of variable variants and
/// the standardized number of alleles round the same way.
///
/// It is 0 where the population called the allele fewer than
/// `num_called_alleles` times, no draw of that size being all of an allele
/// with too few copies, and 1 where it called nothing else.
fn chance_a_draw_is_all_of_one_allele(
    called_alleles: u32,
    count_of_the_allele: u32,
    num_called_alleles: u32,
) -> f64 {
    // The copies of every other allele. A count of one allele is one part of
    // the sum the called alleles are, so it is never the larger of the two
    // and the subtraction never saturates.
    let of_the_other_alleles = called_alleles.saturating_sub(count_of_the_allele);
    chance_a_draw_misses_an_allele(called_alleles, of_the_other_alleles, num_called_alleles)
}

#[cfg(test)]
mod fixtures {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::{Path, PathBuf};

    use crate::block::{Block, BlockReader};
    use crate::error::Result;
    use crate::filters::FilteringStats;
    use crate::io::vcf::{VcfOptions, VcfReader};
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

    /// The reference files live at the root of the repository, beside the
    /// Python tests that read the same files, and not inside this crate.
    fn reference(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference")
            .join(name)
    }

    /// The reader over the panel of `docs/specs/stats.md`, the 1200
    /// biallelic diploid variants of 200 individuals with 3 in 100
    /// genotypes missing that every number of the panel in this module was
    /// measured on, with every variant given and the blocks of the size
    /// popnei chose for its individuals.
    pub(super) fn the_panel() -> VcfReader<BufReader<File>> {
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block: None,
        };
        VcfReader::<BufReader<File>>::from_path(&reference("stats/panel.vcf.gz"), options)
            .unwrap_or_else(|error| panic!("the panel: {error}"))
    }

    /// The individuals of `p0`, `p1` and `p2` of the panel, in that order,
    /// as indices among the `individuals` of its reader.
    ///
    /// `tests/reference/stats/panel_pops_bcftools.txt` names them, one line
    /// of an individual and its population, and it names its three
    /// populations in the order p0, p2, p1, so the populations of a result
    /// are the ones asked for here and not the ones the file happened to
    /// name first.
    pub(super) fn the_pops_of_the_panel(individuals: &[String]) -> Vec<Vec<usize>> {
        let name = "stats/panel_pops_bcftools.txt";
        let text = std::fs::read_to_string(reference(name)).expect("the populations of the panel");
        let mut of_each_pop: Vec<Vec<usize>> = vec![Vec::new(); 3];
        for line in text.lines() {
            let mut columns = line.split('\t');
            let (Some(individual), Some(pop)) = (columns.next(), columns.next()) else {
                panic!("the line `{line}` of {name} is not an individual and a population");
            };
            let at = match pop {
                "p0" => 0,
                "p1" => 1,
                "p2" => 2,
                other => panic!("the population `{other}` of {name} is not p0, p1 or p2"),
            };
            let of_the_individual = individuals
                .iter()
                .position(|named| named == individual)
                .unwrap_or_else(|| panic!("`{individual}` of {name} is not an individual"));
            of_each_pop[at].push(of_the_individual);
        }
        of_each_pop
    }

    /// A reader over one haploid variant at which each population called the
    /// copies of each allele that `of_each_pop` gives, and the individuals of
    /// each of those populations.
    ///
    /// `of_each_pop` is one case of
    /// `tests/reference/diversity/enumerate_private.tsv`, whose
    /// `allele_counts` field is the copies of the allele 0, of the allele 1
    /// and so on that each population called at one variant. One haploid
    /// individual holds one allele, so the row is the copies of the first
    /// population one after another, then those of the second, and the
    /// individuals of a population are the places its own copies took. A
    /// population that called nothing of an allele its neighbour called gets
    /// no individual for it, which is how `3,0,1` puts the alleles 0 and 2 in
    /// one population and leaves the allele 1 to another.
    ///
    /// The variant is haploid so that a count of an allele is a count of
    /// individuals and a reader of the test sees the case of the file in the
    /// row. What the values of this file rest on is the allele counts of each
    /// population and the size of the draw, which a diploid variant of half
    /// called genotypes would give the same.
    pub(super) fn a_variant_of_the_allele_counts(
        of_each_pop: &[&[u32]],
    ) -> (GivenBlocks, Vec<Vec<usize>>) {
        let mut row: Vec<i8> = Vec::new();
        let mut of_each_pops_individuals: Vec<Vec<usize>> = Vec::new();
        for counts in of_each_pop {
            let mut individuals = Vec::new();
            for (allele, count) in counts.iter().enumerate() {
                let allele = i8::try_from(allele).expect("an allele of a case below 128");
                for _ in 0..*count {
                    individuals.push(row.len());
                    row.push(allele);
                }
            }
            of_each_pops_individuals.push(individuals);
        }
        let num_individuals = row.len();
        let rows: Vec<&[i8]> = vec![&row[..]];
        let reader = GivenBlocks::of_a_source_of(
            num_individuals,
            1,
            blocks_of(&rows, num_individuals, 1, 1),
        );
        (reader, of_each_pops_individuals)
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
        a_source_with_two_rows_below_the_missing_allele, a_variant_of_the_allele_counts, blocks_of,
        the_panel, the_pops_of_the_panel, the_worked_example,
    };
    use super::{DiversityOptions, DiversityStats, PopDiversity, calc_pop_diversity};
    use crate::block::BlockReader;
    use crate::error::Error;
    use crate::variant::{MISSING_ALLELE, Needs};

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

    /// The five statistics at a draw of `num_called_alleles` called alleles,
    /// the folded spectrum among them, which is what a pass given a draw can
    /// be asked for.
    fn options_of_a_draw(min_num_individuals: u32, num_called_alleles: u32) -> DiversityOptions {
        DiversityOptions {
            stats: DiversityStats::ALL,
            num_called_alleles: Some(num_called_alleles),
            min_num_individuals,
        }
    }

    /// The diversity of the two populations of the worked example over its
    /// six variants in one block, at a draw of `num_called_alleles` and a
    /// threshold of one called genotype.
    fn of_the_worked_example_at_a_draw(num_called_alleles: u32) -> PopDiversity {
        let mut reader = the_worked_example(6);
        calc_pop_diversity(
            &mut reader,
            &[&POP1, &POP2],
            &options_of_a_draw(1, num_called_alleles),
        )
        .expect("the diversity of the worked example in a draw")
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
        assert_eq!(
            one.num_vars_every_pop_in_draw(),
            other.num_vars_every_pop_in_draw(),
            "the variants in the draw for every population of {what}"
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
            let (alleles, varies, num_vars_in_draw) =
                one.the_sums_of_the_draw(pop).expect("the sums of the draw");
            let (of_the_other, varies_of_the_other, num_vars_of_the_other) = other
                .the_sums_of_the_draw(pop)
                .expect("the sums of the draw");

            assert_eq!(
                num_vars_in_draw, num_vars_of_the_other,
                "the variants in the draw for the population {pop} of {what}"
            );
            assert_eq!(
                alleles.to_bits(),
                of_the_other.to_bits(),
                "the bits of the alleles a draw shows in the population {pop} of {what}, \
                 {alleles} and {of_the_other}"
            );
            assert_eq!(
                varies.to_bits(),
                varies_of_the_other.to_bits(),
                "the bits of the chances a draw varies in the population {pop} of {what}, \
                 {varies} and {varies_of_the_other}"
            );
            let private = one
                .the_sum_of_the_private_alleles_of_the_draw(pop)
                .expect("the sum of the private alleles of the draw");
            let of_the_other = other
                .the_sum_of_the_private_alleles_of_the_draw(pop)
                .expect("the sum of the private alleles of the draw");

            assert_eq!(
                private.to_bits(),
                of_the_other.to_bits(),
                "the bits of the private alleles a draw shows in the population {pop} of {what}, \
                 {private} and {of_the_other}"
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

    /// What a standardized value of these tests may differ from the number
    /// of the spec by. The spec prints 2.3111111111 and 0.6444444444 to ten
    /// decimals, so a literal here is up to 5e-11 from the value it stands
    /// for, and what the arithmetic adds is far below that: each
    /// standardized value of the worked example is a few products of at
    /// most four factors, each rounding by at most 1.2e-16 of the value.
    const OF_TEN_DECIMALS_OF_A_DRAW: f64 = 5e-11;

    /// It checks the alleles a draw is expected to show in one population,
    /// averaged over the variants in the draw for it.
    fn assert_num_alleles_in_draw(diversity: &PopDiversity, pop: usize, mean: f64, what: &str) {
        let found = diversity
            .num_alleles_in_draw(pop)
            .expect("the alleles a draw shows");

        assert!(
            (found - mean).abs() <= OF_TEN_DECIMALS_OF_A_DRAW,
            "a draw shows {found} alleles in {what}, and it shows {mean}"
        );
    }

    /// It checks the chance that a draw varies in one population, averaged
    /// over the same variants.
    fn assert_variable_vars_ratio_in_draw(
        diversity: &PopDiversity,
        pop: usize,
        ratio: f64,
        what: &str,
    ) {
        let found = diversity
            .variable_vars_ratio_in_draw(pop)
            .expect("the chance a draw varies");

        assert!(
            (found - ratio).abs() <= OF_TEN_DECIMALS_OF_A_DRAW,
            "a draw varies at {found} of the variants of {what}, and it varies at {ratio}"
        );
    }

    /// What a standardized private allele value may differ from the value of
    /// `tests/reference/diversity/enumerate_private.tsv` by: one unit of the
    /// last place of a number near 1, 2.2e-16.
    ///
    /// That file holds each of its values computed twice in exact rational
    /// arithmetic, by the closed form of "What it gives" of "The private
    /// alleles" of `docs/specs/diversity.md` and by an enumeration of every
    /// draw, and prints each to seventeen digits, which reads back as the
    /// float64 nearest the exact value. popnei evaluates that closed form in
    /// float64, every factor of it a product of `num_called_alleles`
    /// divisions, so the two differ by a few units of the last place: over the
    /// 22 pairs the test below asserts, eighteen come out to the bit, three
    /// differ by 5.6e-17 and one by 2.8e-17. Measured on 24 September 2026.
    const OF_THE_ENUMERATION: f64 = f64::EPSILON;

    /// It checks the standardized private alleles of one population: the
    /// alleles a draw is expected to show in it and no draw of the same size
    /// in any other population, averaged over the variants in the draw for
    /// every population.
    ///
    /// `within` is [`OF_THE_ENUMERATION`] for a value of the file of the
    /// enumeration, which is a float64, and [`OF_TEN_DECIMALS_OF_A_DRAW`] for
    /// one the spec prints to ten decimals.
    fn assert_private_alleles_in_draw(
        diversity: &PopDiversity,
        pop: usize,
        value: f64,
        within: f64,
        what: &str,
    ) {
        let found = diversity
            .private_alleles_in_draw(pop)
            .expect("the private alleles a draw shows");

        assert!(
            (found - value).abs() <= within,
            "a draw shows {found} private alleles in {what}, and it shows {value}"
        );
    }

    /// It checks the folded spectrum of one population bin by bin, and that
    /// its bins sum to the variants in the draw for it, which each of those
    /// variants gives one whole variant of, spread over the bins.
    ///
    /// `bins` is the spectrum the reference gives, one value for each count of
    /// the rarer allele from 0 up, and `within` what a bin may differ from it
    /// by: `within` of the value for a bin above 1 and `within` itself for a
    /// smaller one, so that a bin of 0, which the worked example has, is
    /// compared at all.
    fn assert_folded_sfs(
        diversity: &PopDiversity,
        pop: usize,
        bins: &[f64],
        within: f64,
        what: &str,
    ) {
        let found = diversity.folded_sfs(pop).expect("the folded spectrum");

        assert_eq!(
            found.len(),
            bins.len(),
            "the bins of the spectrum of {what} are {found:?}"
        );
        for (count, (found, of_the_reference)) in found.iter().zip(bins).enumerate() {
            assert!(
                (found - of_the_reference).abs() <= within * of_the_reference.abs().max(1.0),
                "the variants of {what} that show {count} copies of the rarer allele are {found}, \
                 and they are {of_the_reference}"
            );
        }
        let num_vars_in_draw = diversity
            .num_vars_in_draw(pop)
            .expect("the variants in the draw");
        let of_every_bin: f64 = found.iter().sum();
        // A count below 2^53 is exact in a float64, and a pass of that many
        // variants reads more rows than any source holds.
        let num_vars_in_draw = num_vars_in_draw as f64;

        assert!(
            (of_every_bin - num_vars_in_draw).abs() <= within * num_vars_in_draw.max(1.0),
            "the bins of the spectrum of {what} sum to {of_every_bin}, and {num_vars_in_draw} \
             variants are in the draw for it"
        );
    }

    /// It checks that two results hold the same folded spectrum of every
    /// population, the same bits in every bin, and that either both hold one
    /// or neither does.
    ///
    /// The bits and not a tolerance, for the reason
    /// [`assert_the_same_numbers`] gives: what the tests that call this ask is
    /// whether the bins were added in the order of the variants, and a bin
    /// whose parts were joined in another order is right to far more digits
    /// than any tolerance of the spec and is not the same number.
    ///
    /// It is apart from [`assert_the_same_numbers`] because one test compares
    /// a pass given a draw above every called allele, which holds a spectrum
    /// of zeros, with a pass given no draw, which can hold no spectrum at all:
    /// that test reads the zeros of the first on its own.
    fn assert_the_same_bins(one: &PopDiversity, other: &PopDiversity, what: &str) {
        assert_eq!(
            one.num_pops(),
            other.num_pops(),
            "the populations of {what}"
        );
        for pop in 0..one.num_pops() {
            let of_the_bits = |diversity: &PopDiversity| {
                diversity
                    .folded_sfs(pop)
                    .map(|bins| bins.iter().map(|bin| bin.to_bits()).collect::<Vec<u64>>())
            };

            assert_eq!(
                of_the_bits(one),
                of_the_bits(other),
                "the bits of the bins of the spectrum of the population {pop} of {what}, \
                 {:?} and {:?}",
                one.folded_sfs(pop),
                other.folded_sfs(pop)
            );
        }
    }

    /// It checks that a standardized value of one population is NaN and that
    /// every bin of its folded spectrum is 0, which is what they are when no
    /// variant is in the draw for it: the case of "A population for which no
    /// variant counted" and the one of a `num_called_alleles` above every
    /// called allele, both of "The cases" of `docs/specs/diversity.md`.
    ///
    /// The private alleles of a draw are among them here because both tests
    /// that call this are of a pass no variant reached the draw at for any
    /// population, so no variant is in the draw for every population either
    /// and their divisor is 0 too.
    fn assert_no_standardized_value(diversity: &PopDiversity, pop: usize, what: &str) {
        let alleles = diversity
            .num_alleles_in_draw(pop)
            .expect("the alleles a draw shows");
        let ratio = diversity
            .variable_vars_ratio_in_draw(pop)
            .expect("the chance a draw varies");
        let private = diversity
            .private_alleles_in_draw(pop)
            .expect("the private alleles a draw shows");

        assert!(
            alleles.is_nan(),
            "a draw shows {alleles} alleles in {what}, and it has no value there"
        );
        assert!(
            ratio.is_nan(),
            "a draw varies at {ratio} of the variants of {what}, and it has no value there"
        );
        assert!(
            private.is_nan(),
            "a draw shows {private} private alleles in {what}, and it has no value there"
        );
        let bins = diversity.folded_sfs(pop).expect("the folded spectrum");

        assert!(
            bins.iter().all(|bin| *bin == 0.0),
            "the spectrum of {what} is {bins:?}, and every bin of it is 0"
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

    /// `WITHOUT_A_DRAW` holds the four statistics that need no draw and
    /// not the folded spectrum, and it holds the rest of `ALL`, so that a
    /// statistic added to the module is in one of the two sets.
    #[test]
    fn without_a_draw_holds_the_four_statistics_that_need_no_draw() {
        for stat in [
            DiversityStats::NUM_ALLELES,
            DiversityStats::PRIVATE_ALLELES,
            DiversityStats::VARIABLE_VARS_RATIO,
            DiversityStats::FIS,
        ] {
            assert!(
                DiversityStats::WITHOUT_A_DRAW.contains(stat),
                "the statistics that need no draw hold {:?}",
                stat.names()
            );
        }
        assert!(!DiversityStats::WITHOUT_A_DRAW.contains(DiversityStats::FOLDED_SFS));
        // The two sets are the whole module between them. A sixth statistic
        // left out of either, and a `WITHOUT_A_DRAW` written as a number
        // that the four constants moved away from, fail here.
        assert_eq!(
            DiversityStats::WITHOUT_A_DRAW | DiversityStats::FOLDED_SFS,
            DiversityStats::ALL
        );
        assert_eq!(
            DiversityStats::WITHOUT_A_DRAW.names(),
            [
                "num_alleles",
                "private_alleles",
                "variable_vars_ratio",
                "fis"
            ]
        );
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

    /// A draw of every gene copy the dataset holds is taken and a draw of one
    /// more is refused: that product, the individuals times the ploidy, is
    /// the largest number of alleles any population can have called at any
    /// variant, so a draw above it is what a user wrote and not a fact about
    /// the data.
    ///
    /// The source is 200 diploid individuals, so 400 is the largest draw it
    /// allows and 401 is refused. Its two variants each leave the first
    /// individual's genotype missing, so the one population called 398
    /// alleles at both and neither variant is in the draw at 400: that is the
    /// case of a `num_called_alleles` above every population's called alleles
    /// of "The cases" of `docs/specs/diversity.md`, which the refusal must
    /// leave as it is, NaN in the standardized values and a spectrum of
    /// zeros, and not an error.
    #[test]
    fn a_draw_of_every_gene_copy_of_the_dataset_is_taken_and_one_more_is_refused() {
        let of_the_source: Vec<Vec<i8>> = (1..=2)
            .map(|variant| {
                [MISSING_ALLELE, MISSING_ALLELE]
                    .into_iter()
                    .chain((1..200).flat_map(|individual| [0, i8::from(individual % variant == 0)]))
                    .collect()
            })
            .collect();
        let rows: Vec<&[i8]> = of_the_source.iter().map(|row| &row[..]).collect();

        let of_a_draw_of_400 = calc_pop_diversity(
            &mut GivenBlocks::of_a_source_of(200, 2, blocks_of(&rows, 200, 2, 2)),
            &[],
            &options_of_a_draw(1, 400),
        )
        .expect("the diversity at a draw of every gene copy of the dataset");
        let error = calc_pop_diversity(
            &mut GivenBlocks::of_a_source_of(200, 2, blocks_of(&rows, 200, 2, 2)),
            &[],
            &options_of_a_draw(1, 401),
        )
        .expect_err("a draw of one allele more than the dataset holds");

        assert_eq!(of_a_draw_of_400.num_vars(0), Some(2));
        assert_eq!(of_a_draw_of_400.num_vars_in_draw(0), Some(0));
        assert_no_standardized_value(
            &of_a_draw_of_400,
            0,
            "the one population of 200 individuals at a draw of 400",
        );
        assert!(
            matches!(
                error,
                Error::DiversityDrawLargerThanTheDataset {
                    num_called_alleles: 401,
                    largest_draw: 400,
                    num_individuals: 200,
                    ploidy: 2
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

    /// A variant is in the draw for a population when it counts for that
    /// population and the population called at least `num_called_alleles`
    /// alleles there, both and not the second alone. At a draw of 4 `pop1`
    /// of the worked example keeps all four of the variants that count for
    /// it and `pop2` keeps three of its four: at variant 2 it called 3
    /// copies in all. Variant 6, where `pop1` called one copy and `pop2`
    /// none, counts for no population at a threshold of one called genotype
    /// and so is in the draw for none either, although `pop2` called fewer
    /// than 4 alleles at it as it did at variant 2. "Its Python function"
    /// of `docs/specs/diversity.md`.
    #[test]
    fn a_draw_of_four_keeps_the_four_variants_of_pop1_of_the_worked_example_and_three_of_pop2() {
        let diversity = of_the_worked_example_at_a_draw(4);

        assert_eq!(diversity.num_vars(0), Some(4));
        assert_eq!(diversity.num_vars_in_draw(0), Some(4));
        assert_eq!(diversity.num_vars(1), Some(4));
        assert_eq!(diversity.num_vars_in_draw(1), Some(3));
        assert_eq!(diversity.num_vars_in_draw(2), None);
        // The variants 1, 3 and 5 are in the draw for both populations, and
        // variant 2 for `pop1` alone, so three of the four that counted for
        // both are in the draw for both.
        assert_eq!(diversity.num_vars_every_pop(), 4);
        assert_eq!(diversity.num_vars_every_pop_in_draw(), 3);
    }

    /// A draw of 4 called alleles shows 2.25 alleles in `pop1` of the worked
    /// example and 2.3111111111 in `pop2`: "How it is verified" of "The
    /// number of alleles" of `docs/specs/diversity.md`. Every draw of 4 of
    /// the 4 alleles `pop1` called at each of its variants shows what it
    /// holds, so its standardized value is its mean, 2.25; `pop2` keeps
    /// three variants, and at variant 5, where it called 4 copies of one
    /// allele and 2 of another, a draw of 4 misses the rarer allele in one
    /// of the 15 draws.
    #[test]
    fn a_draw_of_four_shows_2_25_alleles_in_pop1_of_the_worked_example_and_2_3111_in_pop2() {
        let diversity = of_the_worked_example_at_a_draw(4);

        assert_num_alleles_in_draw(&diversity, 0, 2.25, "pop1 at a draw of 4");
        assert_num_alleles_in_draw(&diversity, 1, 2.3111111111, "pop2 at a draw of 4");
        assert_eq!(diversity.num_alleles_in_draw(2), None);
    }

    /// A draw of 4 called alleles varies at 0.75 of the variants of `pop1`
    /// of the worked example and at 0.6444444444 of those of `pop2`: "How it
    /// is verified" of "The variable variants" of `docs/specs/diversity.md`.
    /// `pop1` keeps its four variants and every draw of 4 of 4 shows what it
    /// holds, so its standardized ratio is its ratio; `pop2` keeps three,
    /// and its variant 5 varies in a draw of 4 with chance 14/15.
    ///
    /// The two are not the two standardized numbers of alleles minus 1 here,
    /// as they are on a dataset of biallelic variants: variant 3 of the
    /// worked example has four alleles, and a draw of 4 that varies there
    /// shows two, three or four of them.
    #[test]
    fn a_draw_of_four_varies_at_three_quarters_of_the_variants_of_pop1_and_0_6444_of_pop2() {
        let diversity = of_the_worked_example_at_a_draw(4);

        assert_variable_vars_ratio_in_draw(&diversity, 0, 0.75, "pop1 at a draw of 4");
        assert_variable_vars_ratio_in_draw(&diversity, 1, 0.6444444444, "pop2 at a draw of 4");
        assert_eq!(diversity.variable_vars_ratio_in_draw(2), None);
    }

    /// The folded spectrum of the worked example at a draw of 4 called
    /// alleles: 1, 3 and 0 variants at 0, 1 and 2 copies of the rarer allele
    /// for `pop1`, and 1.0666666667, 1.5333333333 and 0.4 for `pop2`, which
    /// "How it is verified" of "The folded site frequency spectrum" of
    /// `docs/specs/diversity.md` works out variant by variant.
    ///
    /// `pop1` keeps its four variants and every draw is of the 4 copies it
    /// called, so each variant is one whole variant in one bin: variants 1, 2
    /// and 3 show 1, 1 and 3 copies of the rarer allele, which fold to the bin
    /// 1, and variant 5, where it called one allele, shows 0. `pop2` keeps
    /// three variants: variant 1, where it called 5 copies of one allele, is
    /// one whole variant in the bin 0; variant 3, where it called one copy each
    /// of four alleles, is one in the bin 1, a draw of 4 of 4 showing 3 copies
    /// of the rarer allele; and variant 5, where it called 4 copies of one
    /// allele and 2 of another, is spread over the three bins by a draw of 4 of
    /// 6, which shows 0, 1 and 2 rarer copies with the chances 1/15, 8/15 and
    /// 6/15.
    ///
    /// The literals are those chances as fractions and not the ten decimals of
    /// the spec, so what a bin may differ from one by is the rounding of a few
    /// products and not the 5e-11 a printed value stands for. A pytest test of
    /// the same worked example reads the same numbers through the Python layer.
    #[test]
    fn a_draw_of_four_spreads_the_variants_of_the_worked_example_over_three_bins() {
        let diversity = of_the_worked_example_at_a_draw(4);

        assert_folded_sfs(
            &diversity,
            0,
            &[1.0, 3.0, 0.0],
            OF_AN_EXACT_QUOTIENT,
            "pop1 at a draw of 4",
        );
        assert_folded_sfs(
            &diversity,
            1,
            &[1.0 + 1.0 / 15.0, 1.0 + 8.0 / 15.0, 6.0 / 15.0],
            OF_TEN_DECIMALS_OF_A_DRAW,
            "pop2 at a draw of 4",
        );
        assert_eq!(diversity.folded_sfs(2), None);
    }

    /// A draw of 3 called alleles has two bins and not four, the counts 2 and
    /// 3 of the rarer allele folding onto the counts 1 and 0: the bins of a
    /// draw of `g` are the counts from 0 to `g / 2`, which "What it gives" of
    /// "The folded site frequency spectrum" of `docs/specs/diversity.md`
    /// gives.
    ///
    /// The draw of 3 is also the one where the count of the rarer allele
    /// cannot start at 0 for every variant of the worked example: at variant 3
    /// each population called one copy each of four alleles, so a draw of 3 of
    /// those 4 copies leaves at most one copy of the major allele out and shows
    /// 2 or 3 rarer copies, never 0 or 1. Those two chances are 3/4 and 1/4 and
    /// they fall in the bins 1 and 0.
    ///
    /// `pop1` keeps its four variants and `pop2` all four of its own, variant
    /// 2 reaching the draw here where a draw of 4 left it out, its 3 called
    /// copies being one short of 4. The bins are 7/4 and 9/4 for `pop1` and
    /// 49/20 and 31/20 for `pop2`, worked out variant by variant in exact
    /// fractions: 1/4 and 3/4 at each of the variants 1, 2 and 3 of `pop1` and
    /// 1 and 0 at its variant 5; and 1 and 0 at the variants 1 and 2 of
    /// `pop2`, 1/4 and 3/4 at its variant 3, and 1/5 and 4/5 at its variant 5,
    /// where a draw of 3 of 6 shows 0, 1 or 2 of the 2 rarer copies with the
    /// chances 4/20, 12/20 and 4/20 and the last two fall in one bin.
    #[test]
    fn a_draw_of_three_folds_the_counts_above_half_of_it_onto_the_ones_below() {
        let diversity = of_the_worked_example_at_a_draw(3);

        assert_eq!(diversity.num_vars_in_draw(0), Some(4));
        assert_eq!(diversity.num_vars_in_draw(1), Some(4));
        assert_folded_sfs(
            &diversity,
            0,
            &[7.0 / 4.0, 9.0 / 4.0],
            OF_AN_EXACT_QUOTIENT,
            "pop1 at a draw of 3",
        );
        assert_folded_sfs(
            &diversity,
            1,
            &[49.0 / 20.0, 31.0 / 20.0],
            OF_TEN_DECIMALS_OF_A_DRAW,
            "pop2 at a draw of 3",
        );
    }

    /// A count of the rarer allele above the draw is in no bin of the
    /// spectrum, and one that the copies of the major allele leave no room for
    /// is in none either: the count runs from `max(0, g - (c - m))` to
    /// `min(m, g)`, which "What it gives" of "The folded site frequency
    /// spectrum" of `docs/specs/diversity.md` gives with `c` the copies the
    /// population called, `m` the ones that are not of its major allele and
    /// `g` the draw.
    ///
    /// The population called 5 copies each of three alleles, 15 in all, and
    /// the draw is of 8. Its major allele has 5 copies, whichever of the three
    /// it is, so 10 copies are of the rarer allele, above the 8 of the draw,
    /// and the draw holds at least 3 of them, the 5 copies of the major allele
    /// leaving room for no more than 5. The counts from 3 to 8 have the
    /// chances 120, 1050, 2520, 2100, 600 and 45 over the 6435 draws of 8 of
    /// 15, and they fall in the bins 3, 4, 3, 2, 1 and 0, so the bins are
    /// 1/143, 40/429, 140/429, 16/39 and 70/429, worked out in exact
    /// fractions.
    ///
    /// A count of 9 or 10 rarer copies is more than the draw takes, and a
    /// version that took it would ask for 8 - 9 copies of the major allele and
    /// a bin below 0. One that read such a count as 0 copies of the major
    /// allele and put it in the bin 0 gives a number that is wrong and not a
    /// panic, which is what this case is here for: it would add 10 and 1 draws
    /// of the 6435 to the bin 0 and make it 56/6435 where it is 45/6435.
    ///
    /// The three alleles are called equally often, so which of them is the
    /// major one is settled by the rule of `docs/glossary.md` that it is the
    /// lower numbered of two that tie; the count of its copies, which is what
    /// the spectrum reads, is 5 whichever it is.
    #[test]
    fn a_count_of_the_rarer_allele_above_the_draw_is_in_no_bin() {
        let (mut reader, of_each_pops_individuals) = a_variant_of_the_allele_counts(&[&[5, 5, 5]]);
        let pops: Vec<&[usize]> = of_each_pops_individuals
            .iter()
            .map(|individuals| &individuals[..])
            .collect();

        let diversity = calc_pop_diversity(&mut reader, &pops, &options_of_a_draw(1, 8))
            .expect("the diversity of a variant of three alleles of five copies each");

        assert_eq!(diversity.num_vars_in_draw(0), Some(1));
        assert_folded_sfs(
            &diversity,
            0,
            &[
                1.0 / 143.0,
                40.0 / 429.0,
                140.0 / 429.0,
                16.0 / 39.0,
                70.0 / 429.0,
            ],
            OF_TEN_DECIMALS_OF_A_DRAW,
            "one population of 5 copies each of three alleles at a draw of 8",
        );
    }

    /// A pass asked for the folded spectrum alone gives it and gives no other
    /// value, and its bins are the ones of a pass asked for the five
    /// statistics: what the spectrum reads of a variant is the counts of the
    /// alleles the pass takes for every statistic, and the totals the other
    /// four keep are not among them.
    #[test]
    fn a_pass_asked_for_the_spectrum_alone_gives_it_and_no_other_value() {
        let mut reader = the_worked_example(6);
        let options = DiversityOptions {
            stats: DiversityStats::FOLDED_SFS,
            num_called_alleles: Some(4),
            min_num_individuals: 1,
        };

        let diversity = calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options)
            .expect("the diversity with the folded spectrum alone");

        assert_eq!(diversity.num_vars_in_draw(0), Some(4));
        assert_eq!(diversity.num_vars_in_draw(1), Some(3));
        assert_eq!(diversity.num_alleles(0), None);
        assert_eq!(diversity.num_alleles_in_draw(0), None);
        assert_eq!(diversity.private_alleles(0), None);
        assert_eq!(diversity.num_variable_vars(0), None);
        assert_eq!(diversity.fis(0), None);
        assert_the_same_bins(
            &diversity,
            &of_the_worked_example_at_a_draw(4),
            "a pass asked for the spectrum alone",
        );
    }

    /// A population for which no variant is in the draw has 0 in its count
    /// of them and NaN in both standardized values, which is the case of "A
    /// population for which no variant counted" of "The cases" of
    /// `docs/specs/diversity.md` and not an error: the other populations of
    /// the call still have their values.
    ///
    /// At a threshold of three called genotypes no variant counts for
    /// `pop1`, which holds two individuals, so none is in the draw for it
    /// either. Variant 5 counts for `pop2`, whose three individuals called
    /// 4 copies of one allele and 2 of another there, and a draw of 4 of
    /// those 6 varies with chance 14/15 and shows 1 + 14/15 alleles, the
    /// variant having two alleles and a draw showing one of them or both.
    #[test]
    fn a_population_no_variant_of_the_draw_counted_for_has_no_standardized_value() {
        let mut reader = the_worked_example(6);

        let diversity = calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options_of_a_draw(3, 4))
            .expect("the diversity at a threshold of three called genotypes");

        assert_eq!(diversity.num_vars_in_draw(0), Some(0));
        assert_no_standardized_value(&diversity, 0, "pop1 at a threshold of three genotypes");
        assert_eq!(diversity.num_vars_in_draw(1), Some(1));
        assert_variable_vars_ratio_in_draw(&diversity, 1, 14.0 / 15.0, "pop2 over variant 5");
        assert_num_alleles_in_draw(&diversity, 1, 1.0 + 14.0 / 15.0, "pop2 over variant 5");
        assert_eq!(diversity.num_vars_every_pop_in_draw(), 0);
    }

    /// A `num_called_alleles` above every population's called alleles leaves
    /// both counts of the draw at 0 and both standardized values NaN, and
    /// leaves the totals, the means and F_IS as they are without it, none of
    /// them reading the draw: "The cases" of `docs/specs/diversity.md`. It is
    /// not an error, and the counts say why the values are missing.
    ///
    /// `pop1` of the worked example called 4 alleles at the most at a
    /// variant and `pop2` 6, so a draw of 7 reaches no variant of either.
    #[test]
    fn a_draw_above_every_populations_called_alleles_leaves_the_totals_as_they_were() {
        let mut reader = the_worked_example(6);
        let of_a_draw_of_seven =
            calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options_of_a_draw(1, 7))
                .expect("the diversity at a draw of seven called alleles");
        let of_no_draw = of_the_worked_example(1, 6);

        for (pop, what) in [(0, "pop1 at a draw of 7"), (1, "pop2 at a draw of 7")] {
            assert_eq!(of_a_draw_of_seven.num_vars_in_draw(pop), Some(0));
            assert_no_standardized_value(&of_a_draw_of_seven, pop, what);
        }
        assert_eq!(of_a_draw_of_seven.num_vars_every_pop_in_draw(), 0);
        assert_the_same_numbers(
            &of_a_draw_of_seven,
            &of_no_draw,
            "a draw above every called allele",
        );
    }

    /// A standardized value has no value when its statistic was not asked
    /// for, and the variants in the draw are counted whatever the pass was
    /// asked for, since their count is what says why a standardized value is
    /// missing.
    #[test]
    fn a_standardized_value_of_a_statistic_that_was_not_asked_for_has_no_value() {
        let of_the_alleles = {
            let mut reader = the_worked_example(6);
            let options = DiversityOptions {
                stats: DiversityStats::NUM_ALLELES,
                num_called_alleles: Some(4),
                min_num_individuals: 1,
            };
            calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options)
                .expect("the diversity with the alleles called alone")
        };
        let of_the_fis = {
            let mut reader = the_worked_example(6);
            let options = DiversityOptions {
                stats: DiversityStats::FIS,
                num_called_alleles: Some(4),
                min_num_individuals: 1,
            };
            calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options)
                .expect("the diversity with the F_IS alone")
        };

        assert_num_alleles_in_draw(
            &of_the_alleles,
            0,
            2.25,
            "pop1 of a pass of the alleles alone",
        );
        assert_eq!(of_the_alleles.variable_vars_ratio_in_draw(0), None);
        assert_eq!(of_the_alleles.folded_sfs(0), None);
        assert_eq!(of_the_fis.num_alleles_in_draw(0), None);
        assert_eq!(of_the_fis.variable_vars_ratio_in_draw(0), None);
        assert_eq!(of_the_fis.folded_sfs(0), None);
        assert_eq!(of_the_fis.num_vars_in_draw(0), Some(4));
        assert_eq!(of_the_fis.num_vars_in_draw(1), Some(3));
        assert_eq!(of_the_fis.num_vars_every_pop_in_draw(), 3);
    }

    /// What a value of the panel may differ from the one the reference program
    /// gave by: 1e-12 of it, which is what "How it is verified" of "The number
    /// of alleles" and of "The folded site frequency spectrum" of
    /// `docs/specs/diversity.md` compare popnei within, against
    /// `vegan::rarefy` and against `dadi`, both sides summing the same per
    /// variant values in different orders.
    const OF_THE_PANEL: f64 = 1e-12;

    /// It checks one standardized value of the panel against `vegan`'s.
    fn assert_of_the_panel(found: f64, of_vegan: f64, what: &str) {
        assert!(
            ((found - of_vegan) / of_vegan).abs() <= OF_THE_PANEL,
            "{what} of the panel is {found}, and vegan gives {of_vegan}"
        );
    }

    /// The standardized values of the three populations of the panel at a
    /// draw of 20 called alleles with `min_num_individuals` 20, which the
    /// tables of "How it is verified" of "The number of alleles" and of "The
    /// variable variants" of `docs/specs/diversity.md` give.
    ///
    /// `vegan` 2.7.6 measured the alleles a draw of 20 shows, 1.9283948650,
    /// 1.9219209943 and 1.9197370844, and
    /// `tests/reference/diversity/panel_num_alleles.tsv` holds them to
    /// seventeen digits. The standardized ratio of variable variants is each
    /// of those minus 1, because every variant of the panel has two alleles
    /// and a draw there shows one of them or both;
    /// `tests/reference/diversity/panel_variable_vars.tsv` holds it under a
    /// name that says it is not a second measurement of `vegan`'s, so what
    /// the three ratios check is popnei's second formula against that one
    /// measurement and against the identity.
    ///
    /// Each literal below is the number of one of those two files as a
    /// float64 keeps it, which is the same number to fewer digits where the
    /// file printed more than a float64 holds: the file gives `p1`
    /// 1.9219209943237829 and the nearest float64 to it is
    /// 1.921920994323783.
    ///
    /// The totals beside them are `adegenet` 2.1.11's 2373, 2377 and 2384
    /// alleles called and 1173, 1177 and 1184 variable variants, over the
    /// 1200 variants that counted for every population and that every
    /// population reached the draw at.
    #[test]
    fn the_standardized_values_of_the_panel_are_the_ones_vegan_gave() {
        let diversity = of_the_panel_at_a_draw(20);

        assert_eq!(diversity.num_vars_of_the_pass(), 1200);
        assert_eq!(diversity.num_vars_every_pop(), 1200);
        assert_eq!(diversity.num_vars_every_pop_in_draw(), 1200);
        let of_the_reference = [
            (
                0,
                "p0",
                2373,
                1173,
                1.928_394_865_004_120_5,
                0.928_394_865_004_120_5,
            ),
            (
                1,
                "p1",
                2377,
                1177,
                1.921_920_994_323_783,
                0.921_920_994_323_782_9,
            ),
            (
                2,
                "p2",
                2384,
                1184,
                1.919_737_084_393_756_2,
                0.919_737_084_393_756_2,
            ),
        ];
        for (pop, name, num_alleles, num_variable_vars, alleles_of_a_draw, ratio_of_a_draw) in
            of_the_reference
        {
            assert_eq!(
                diversity.num_vars(pop),
                Some(1200),
                "the variants of {name}"
            );
            assert_eq!(
                diversity.num_vars_in_draw(pop),
                Some(1200),
                "the variants in the draw for {name}"
            );
            assert_eq!(
                diversity.num_alleles(pop),
                Some(num_alleles),
                "the alleles of {name}"
            );
            assert_eq!(
                diversity.num_variable_vars(pop),
                Some(num_variable_vars),
                "the variable variants of {name}"
            );
            assert_of_the_panel(
                diversity
                    .num_alleles_in_draw(pop)
                    .expect("the alleles a draw shows"),
                alleles_of_a_draw,
                &format!("the alleles a draw of 20 shows in {name}"),
            );
            assert_of_the_panel(
                diversity
                    .variable_vars_ratio_in_draw(pop)
                    .expect("the chance a draw varies"),
                ratio_of_a_draw,
                &format!("the variants a draw of 20 varies at in {name}"),
            );
        }
    }

    /// The diversity of the three populations of the panel at a draw of
    /// `num_called_alleles` called alleles and `min_num_individuals` 20,
    /// which is what every number of the panel in this module was measured
    /// with.
    fn of_the_panel_at_a_draw(num_called_alleles: u32) -> PopDiversity {
        let mut reader = the_panel();
        let of_each_pop = the_pops_of_the_panel(reader.individuals());
        let pops: Vec<&[usize]> = of_each_pop.iter().map(|pop| &pop[..]).collect();
        calc_pop_diversity(
            &mut reader,
            &pops,
            &options_of_a_draw(20, num_called_alleles),
        )
        .expect("the diversity of the panel")
    }

    /// The folded spectrum of the three populations of the panel at a draw of
    /// 20 called alleles with `min_num_individuals` 20, all 33 bins of it,
    /// which `dadi` 2.4.4 measured and
    /// `tests/reference/diversity/panel_folded_sfs_dadi.tsv` holds. The table
    /// of "How it is verified" of "The folded site frequency spectrum" of
    /// `docs/specs/diversity.md` prints the same values to ten decimals.
    ///
    /// Each literal is the value of that file as an `f64` keeps it. The three
    /// columns are one row of the table each here, a count of the rarer allele
    /// with the value of `p0`, of `p1` and of `p2`, which is how the file
    /// reads.
    ///
    /// This is the tightest comparison of the module: an `f64` projection
    /// differs from those values by up to 7.1e-14 of them, `dadi`'s own error
    /// dominating, where the standardized number of alleles against
    /// `vegan::rarefy` sits 550 times inside the same 1e-12. Both were
    /// measured on 24 September 2026 by recomputing the panel in exact
    /// rational arithmetic, which "How it is verified" of that item records.
    /// The columns of the file sum to the 1200 variants in the draw short by
    /// up to 7.0e-11, which is why the sum is compared within a tolerance too;
    /// popnei's own bins are summed here and not the file's.
    ///
    /// `dadi` masks the bin 0 and the bins above half the draw in a folded
    /// spectrum and popnei reports the bin 0, so the file holds the unmasked
    /// array. It is a difference of presentation and not of value.
    #[test]
    fn the_folded_spectrum_of_the_panel_is_the_one_dadi_gave() {
        let of_each_count: [[f64; 3]; 11] = [
            [85.92616199505309, 93.69480681145635, 96.3154987274892],
            [92.99651940325519, 95.45880325926598, 101.37814416196603],
            [106.88963282305795, 103.72623430085372, 108.05152927646787],
            [115.54543767638108, 110.6831031478722, 114.8500079366327],
            [120.50559386225342, 116.6729455230275, 119.73977861754621],
            [122.94119349597617, 121.05177240834107, 121.88331138362474],
            [124.08151270035586, 123.60219049366185, 121.86450896763563],
            [124.23937415023605, 124.59443950850212, 120.60097546763699],
            [123.49607649027779, 124.54032617876919, 118.97681357034003],
            [122.42162180189403, 124.0682314877391, 117.71473832544012],
            [60.95687560124602, 61.907146880441395, 58.62469356518016],
        ];
        let diversity = of_the_panel_at_a_draw(20);

        for (pop, name) in [(0, "p0"), (1, "p1"), (2, "p2")] {
            let of_the_pop: Vec<f64> = of_each_count
                .iter()
                .map(|of_each_pop| of_each_pop[pop])
                .collect();

            assert_eq!(
                diversity.num_vars_in_draw(pop),
                Some(1200),
                "the variants in the draw for {name}"
            );
            assert_folded_sfs(
                &diversity,
                pop,
                &of_the_pop,
                OF_THE_PANEL,
                &format!("{name} of the panel at a draw of 20"),
            );
        }
    }

    /// The standardized private alleles of `p0`, `p1` and `p2` of the panel at
    /// a draw of 20 called alleles: 0.0112196177, 0.0099715392 and
    /// 0.0089014974, which "How it is verified" of "The private alleles" of
    /// `docs/specs/diversity.md` gives to ten decimals and which the literals
    /// here are, so the bound is the 5e-11 a value printed that way stands
    /// for.
    ///
    /// They are the one set of numbers of that spec that no program outside
    /// popnei gives, no outside program computing a standardized private
    /// allele value at all. They come from
    /// `docs/reports/diversity-method/panel.py`, which computes the five
    /// quantities of the spec in Python, so what they check is popnei's Rust
    /// against that Python over 1200 variants, where the `f64` of each is
    /// summed in a different order and each chance of a draw is a ratio of
    /// binomial coefficients there and a product of 20 factors here. What
    /// checks the formula itself is the enumeration of every draw, over the 22
    /// pairs, and the two properties beside it.
    #[test]
    fn the_standardized_private_alleles_of_the_panel_are_the_ones_of_the_spec() {
        let diversity = of_the_panel_at_a_draw(20);

        assert_eq!(diversity.num_vars_every_pop_in_draw(), 1200);
        for (pop, name, of_the_spec) in [
            (0, "p0", 0.011_219_617_7),
            (1, "p1", 0.009_971_539_2),
            (2, "p2", 0.008_901_497_4),
        ] {
            assert_private_alleles_in_draw(
                &diversity,
                pop,
                of_the_spec,
                OF_TEN_DECIMALS_OF_A_DRAW,
                &format!("{name} of the panel at a draw of 20"),
            );
        }
    }

    /// Where every variant is in the draw for every population, the
    /// standardized private alleles of a population are at most the alleles
    /// the same draw shows in it, a private allele being an allele. "How it is
    /// verified" of "The private alleles" of `docs/specs/diversity.md` asks
    /// for it on the panel, with `num_called_alleles` at the smallest count of
    /// called alleles the dataset holds, and that condition is what the two
    /// counts asserted below stand for.
    ///
    /// The draw is of 20 called alleles, which every population of the panel
    /// reaches at every one of its 1200 variants, so the two values are means
    /// over the same variants and the comparison is of the two sums the pass
    /// took over one set of variants: with a draw above one population's
    /// called alleles at some variant the divisors would differ and the
    /// inequality would not follow from the per variant one.
    #[test]
    fn the_private_alleles_of_a_draw_every_variant_reached_are_at_most_the_alleles_it_shows() {
        let diversity = of_the_panel_at_a_draw(20);

        assert_eq!(diversity.num_vars_every_pop_in_draw(), 1200);
        for (pop, name) in [(0, "p0"), (1, "p1"), (2, "p2")] {
            assert_eq!(
                diversity.num_vars_in_draw(pop),
                Some(1200),
                "the variants in the draw for {name}"
            );
            let private = diversity
                .private_alleles_in_draw(pop)
                .expect("the private alleles a draw shows");
            let alleles = diversity
                .num_alleles_in_draw(pop)
                .expect("the alleles a draw shows");

            assert!(
                private <= alleles,
                "a draw of 20 shows {private} private alleles in {name} and {alleles} alleles"
            );
        }
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
        assert_eq!(of_the_alleles.private_alleles_in_draw(0), None);
    }

    /// One case of `tests/reference/diversity/enumerate_private.tsv` as the
    /// test below reads it: the `case` field, the `allele_counts` of each
    /// population, the `num_called_alleles` of the draw and the `closed_form`
    /// of each population, in the order the file names them.
    type OfACaseOfTheEnumeration<'a> = (&'a str, &'a [&'a [u32]], u32, &'a [f64]);

    /// The standardized private alleles of every population of the ten cases
    /// of `tests/reference/diversity/enumerate_private.tsv` whose populations
    /// share no individual: 22 pairs of a case and a population, each of
    /// which that file computes twice in exact rational arithmetic, once by
    /// the closed form popnei evaluates and once by enumerating every draw
    /// each population can make and averaging over every combination of one
    /// draw per population. The two agree with a difference of exactly 0 on
    /// all 22, so a literal below is the value of both ways, and what it
    /// checks is the algebra of the closed form and not a transcription of
    /// it: the enumeration never writes that form down.
    ///
    /// The table is the `case`, `allele_counts`, `num_called_alleles` and
    /// `closed_form` fields of the file, in its order, with the value of each
    /// population of a case in the order the file names them. Three of the
    /// cases are the variants of the worked example that have a draw of 4 and
    /// the other seven were made up for this check, at draws of 2 and 3. Four
    /// of the ten hold a population whose called alleles are not 4: variants 1
    /// and 5 of the worked example, 4 against 5 and 4 against 6; the lone
    /// population of 6 copies; and the three populations of 3, 5 and 6. The
    /// other six hold every population at 4, which is why the file says the
    /// last two cases were added. Taking the called alleles of the first
    /// population of a row for every population of it moves four of the 22
    /// pairs here, `pop2` of variant 5 of the worked example and all three of
    /// the case of 3, 5 and 6, and no pair of any other case: measured on 24
    /// September 2026.
    ///
    /// Each case is one variant, so the mean over the variants in the draw
    /// for every population is the value at that variant.
    #[test]
    fn the_standardized_private_alleles_of_the_ten_cases_are_the_ones_the_enumeration_gives() {
        let of_each_case: [OfACaseOfTheEnumeration; 10] = [
            (
                "worked example, variant 1",
                &[&[3, 1], &[5, 0]],
                4,
                &[1.0, 0.0],
            ),
            (
                "worked example, variant 3",
                &[&[1, 1, 1, 1], &[1, 1, 1, 1]],
                4,
                &[0.0, 0.0],
            ),
            (
                "worked example, variant 5",
                &[&[4, 0], &[4, 2]],
                4,
                &[0.0, 0.9333333333333333],
            ),
            (
                "two populations of two alleles",
                &[&[3, 1], &[2, 2]],
                2,
                &[0.25, 0.4166666666666667],
            ),
            (
                "two populations of three alleles",
                &[&[2, 1, 1], &[3, 0, 1]],
                2,
                &[0.75, 0.4166666666666667],
            ),
            (
                "three populations of two alleles",
                &[&[2, 2], &[3, 1], &[1, 3]],
                2,
                &[0.0, 0.08333333333333333, 0.08333333333333333],
            ),
            (
                "three populations of four alleles",
                &[&[2, 1, 1, 0], &[1, 1, 0, 2], &[2, 0, 1, 1]],
                2,
                &[0.5694444444444444, 0.6805555555555556, 0.4027777777777778],
            ),
            (
                "a population holding one allele against one holding two",
                &[&[4, 0], &[2, 2]],
                3,
                &[0.0, 1.0],
            ),
            ("one population of three alleles", &[&[3, 2, 1]], 3, &[2.25]),
            (
                "three populations of 3, 5 and 6 called alleles",
                &[&[2, 1], &[3, 1, 1], &[2, 2, 1, 1]],
                2,
                &[0.2, 0.32, 0.6533333333333333],
            ),
        ];

        for (case, of_each_pop, num_called_alleles, of_the_file) in of_each_case {
            let (mut reader, of_each_pops_individuals) =
                a_variant_of_the_allele_counts(of_each_pop);
            let pops: Vec<&[usize]> = of_each_pops_individuals
                .iter()
                .map(|individuals| &individuals[..])
                .collect();

            let diversity = calc_pop_diversity(
                &mut reader,
                &pops,
                &options_of_a_draw(1, num_called_alleles),
            )
            .unwrap_or_else(|error| panic!("the diversity of `{case}`: {error}"));

            assert_eq!(
                diversity.num_vars_every_pop_in_draw(),
                1,
                "the variants in the draw for every population of `{case}`"
            );
            for (pop, value) in of_the_file.iter().enumerate() {
                assert_private_alleles_in_draw(
                    &diversity,
                    pop,
                    *value,
                    OF_THE_ENUMERATION,
                    &format!("the population {pop} of `{case}`"),
                );
            }
            assert_eq!(
                diversity.private_alleles_in_draw(of_the_file.len()),
                None,
                "a population `{case}` has not"
            );
        }
    }

    /// The standardized private alleles of the worked example at a draw of 4:
    /// 0.3333333333 for `pop1` and 0.3111111111 for `pop2` over the three
    /// variants in the draw for both, which "How it is verified" of "The
    /// private alleles" of `docs/specs/diversity.md` works out by hand.
    ///
    /// The three variants are 1, 3 and 5. At variant 1 `pop1` draws its 4 of
    /// 4, so the allele 1 it alone called is certain, and `pop2` draws 4 of
    /// the 5 copies of the allele 0 it called and can show nothing else: the
    /// term is 1 for `pop1` and 0 for `pop2`. At variant 3 both populations
    /// called the four alleles once each and draw 4 of 4, so every allele is
    /// certain in both draws and every term is 0. At variant 5 `pop1` called
    /// 4 copies of the allele 0 and `pop2` 4 of the allele 0 and 2 of the
    /// allele 1, which its draw of 4 of 6 shows with chance 1 - C(4, 4) /
    /// C(6, 4), 14/15, and which `pop1` cannot show. So `pop1` sums to 1 and
    /// `pop2` to 0.9333333333, each over the 3 variants.
    ///
    /// Variant 2 is out although it counted for both populations, `pop2`
    /// having called 3 alleles there and not the 4 of the draw, and the
    /// allele 1 that is private to `pop1` at it is out of both standardized
    /// values with it: the divisor of this value is the variants in the draw
    /// for every population and not the variants in the draw for the
    /// population whose value it is.
    #[test]
    fn a_draw_of_four_shows_a_third_of_a_private_allele_in_pop1_of_the_worked_example() {
        let diversity = of_the_worked_example_at_a_draw(4);

        assert_eq!(diversity.num_vars_every_pop_in_draw(), 3);
        assert_private_alleles_in_draw(
            &diversity,
            0,
            0.3333333333,
            OF_TEN_DECIMALS_OF_A_DRAW,
            "pop1 at a draw of 4",
        );
        assert_private_alleles_in_draw(
            &diversity,
            1,
            0.3111111111,
            OF_TEN_DECIMALS_OF_A_DRAW,
            "pop2 at a draw of 4",
        );
        assert_eq!(diversity.private_alleles_in_draw(2), None);
    }

    /// With one population every allele a draw shows is private, there being
    /// no other population to hold one, so the standardized private alleles
    /// are then the standardized number of alleles. "The cases" of
    /// `docs/specs/diversity.md`.
    ///
    /// The population is the one of three alleles of
    /// `tests/reference/diversity/enumerate_private.tsv`, 3, 2 and 1 copies
    /// of 6 called at a draw of 3, where both are 9/4: a draw of 3 of those 6
    /// misses the allele called 3 times in C(3, 3) / C(6, 3) of them, a
    /// twentieth, the one called twice in 4/20 and the one called once in
    /// 10/20, so the alleles it is expected to show are 19/20 + 16/20 + 10/20.
    /// It is also the one pair of that file with a single population and the
    /// tight end of the property that the private alleles of a draw are at
    /// most the alleles it shows.
    #[test]
    fn with_one_population_every_allele_a_draw_shows_is_private() {
        let (mut reader, of_each_pops_individuals) = a_variant_of_the_allele_counts(&[&[3, 2, 1]]);
        let pops: Vec<&[usize]> = of_each_pops_individuals
            .iter()
            .map(|individuals| &individuals[..])
            .collect();

        let diversity = calc_pop_diversity(&mut reader, &pops, &options_of_a_draw(1, 3))
            .expect("the diversity of one population of three alleles");

        let alleles = diversity
            .num_alleles_in_draw(0)
            .expect("the alleles a draw shows");
        assert_private_alleles_in_draw(&diversity, 0, 2.25, OF_THE_ENUMERATION, "one population");
        assert_private_alleles_in_draw(
            &diversity,
            0,
            alleles,
            OF_THE_ENUMERATION,
            "one population, against the alleles its draw shows",
        );
    }

    /// A population against a copy of itself counts no private allele at any
    /// draw size, every allele it called having been called by the copy too:
    /// the totals ask which alleles a population called and not which a draw
    /// of them showed. The two populations are both `pop1` of the worked
    /// example, and the draws are of 2 and of 4 called alleles.
    ///
    /// Its standardized value is 0 only where the draw takes every copy the
    /// population called, which at a draw of 4 it does at each of the four
    /// variants that count for `pop1`: there every allele it holds is certain
    /// in both draws. At a draw of 2 the value is 0.375, the estimator
    /// reading the two draws as draws of copies of their own, which is the
    /// overlap "The cases" of `docs/specs/diversity.md` describes and which
    /// the twenty-third pair of
    /// `tests/reference/diversity/enumerate_private.tsv` shows on the
    /// smallest case there is. A draw of 2 of the 4 copies of variant 1,
    /// where `pop1` called the allele 0 three times and the allele 1 once,
    /// shows the allele 1 in half of the draws and misses it in half, so the
    /// term of that allele is 1/4 and the term of the allele 0, which no draw
    /// of 2 of 4 can miss, is 0. Variants 1 and 2 give 1/4 each, variant 3,
    /// where `pop1` called four alleles once each, gives 4 terms of 1/4, and
    /// variant 5, where it called one allele, gives 0: 1.5 over the 4
    /// variants.
    #[test]
    fn a_population_against_a_copy_of_itself_counts_no_private_allele_and_only_a_draw_of_every_copy_shows_none()
     {
        let of_a_draw = |num_called_alleles| {
            let mut reader = the_worked_example(6);
            calc_pop_diversity(
                &mut reader,
                &[&POP1, &POP1],
                &options_of_a_draw(1, num_called_alleles),
            )
            .expect("the diversity of a population against a copy of itself")
        };

        for num_called_alleles in [2, 4] {
            let diversity = of_a_draw(num_called_alleles);
            assert_eq!(
                diversity.num_vars_every_pop_in_draw(),
                4,
                "the variants in the draw at a draw of {num_called_alleles}"
            );
            for pop in 0..2 {
                assert_eq!(
                    diversity.private_alleles(pop),
                    Some(0),
                    "the private alleles of the population {pop} at a draw of {num_called_alleles}"
                );
            }
        }
        for pop in 0..2 {
            assert_private_alleles_in_draw(
                &of_a_draw(4),
                pop,
                0.0,
                OF_THE_ENUMERATION,
                &format!("the population {pop} at a draw of every copy it called"),
            );
            assert_private_alleles_in_draw(
                &of_a_draw(2),
                pop,
                0.375,
                OF_THE_ENUMERATION,
                &format!("the population {pop} at a draw of 2"),
            );
        }
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
            calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options_of_a_draw(1, 2))
                .expect("the diversity read on the threads")
        };
        let one_chunk_at_a_time = {
            let mut reader = a_source_of_many_variants(200, 200);
            calc_pop_diversity_one_chunk_at_a_time(
                &mut reader,
                &[&POP1, &POP2],
                &options_of_a_draw(1, 2),
            )
            .expect("the diversity read one chunk at a time")
        };

        assert_eq!(of_the_threads.num_vars(0), Some(160));
        assert_eq!(of_the_threads.num_vars_in_draw(0), Some(160));
        assert_the_same_numbers(
            &of_the_threads,
            &one_chunk_at_a_time,
            "the chunks read one after another",
        );
        assert_the_same_bins(
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
    /// The bits of the two sums behind F_IS, of the two sums of the draw and
    /// of every bin of the folded spectrum
    /// are compared and not a tolerance: the chunks are added in the order of
    /// the block, and rayon's own `reduce` over the same chunks joins them in
    /// a tree whose shape follows the threads of the pool, which gives a
    /// number right to far more digits than any tolerance of the spec and not
    /// the same one. With that `reduce` in place of the ordered addition the
    /// sums of `pop1` differ between 1 thread and 2 here, and they do not at
    /// 200 variants, 4 chunks, where rayon splits the same way whatever the
    /// pool.
    ///
    /// The draw is of 2 called alleles, the smallest the spec allows, which
    /// every variant of the source that counts for a population reaches:
    /// `pop1` calls 4 alleles at each of its four patterns and `pop2` between
    /// 3 and 6. A draw of 2 is the size whose two sums of each population
    /// move when their parts are added in another order, and a draw of 3 is
    /// one that does not: over these patterns its terms cancel to the bit.
    /// Measured on 24 September 2026 by adding the chunks of a block in
    /// reverse and reading which sums changed. The bins of the spectrum were
    /// measured the same way on the same day: the 0 bin of `pop1`, the
    /// variants a draw of 2 shows one allele at, moves by one unit of its last
    /// place when the chunks are added in reverse, so the bins are parts of
    /// this test and not only carried by it.
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
                calc_pop_diversity(&mut reader, &[&POP1, &POP2], &options_of_a_draw(1, 2))
                    .expect("the diversity of a source of many variants")
            })
        };

        let on_one = in_a_pool(1);
        assert_eq!(on_one.num_vars(0), Some(1600));
        assert_eq!(on_one.num_vars(1), Some(1600));
        assert_eq!(on_one.num_vars_in_draw(0), Some(1600));
        assert_eq!(on_one.num_vars_in_draw(1), Some(1600));
        for threads in [2, 4, 8] {
            let of_the_pool = in_a_pool(threads);
            assert_the_same_numbers(&on_one, &of_the_pool, &format!("{threads} threads"));
            assert_the_same_bins(&on_one, &of_the_pool, &format!("{threads} threads"));
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

#[cfg(test)]
mod the_chance_a_draw_misses_an_allele {
    use super::chance_a_draw_misses_an_allele;

    /// What a value of these tests may differ from the number of the spec's
    /// table by. `vegan` printed the alleles a draw shows to ten digits, so
    /// a literal taken from that table is up to half of the last one away
    /// from the value it stands for: 1.5333333333 is 3.3e-11 below the 23/15
    /// it is. What the function adds is a few divisions and
    /// multiplications, each rounding by at most 1.2e-16 of the value, so
    /// this bound is the digits of the table and no allowance for the
    /// rounding.
    const OF_TEN_DIGITS: f64 = 5e-11;

    /// What a value of these tests may differ from a fraction the test
    /// itself writes by. The longest of them adds two chances of a draw of
    /// 2, six divisions and multiplications, and the fraction it is held
    /// against rounds once more: nine roundings, each of at most half of the
    /// last bit, 1.1e-16 of the value, so 1.0e-15 for a value near 1. This
    /// bound is eight times the last bit, 1.8e-15, and has nearly twice the
    /// room those roundings need. A factorial in place of the product misses
    /// it by a NaN and not by a last bit.
    const OF_A_FEW_ROUNDINGS: f64 = 8.0 * f64::EPSILON;

    /// The allele counts of one population at the four variants of the
    /// worked example that count for it, and beside each the alleles
    /// `vegan` 2.7.6 expects a draw of 2 to show there. Both columns are
    /// the table at the end of "How it is verified" of "The number of
    /// alleles" of `docs/specs/diversity.md`, whose variants are 1, 2, 3
    /// and 5.
    type AtADrawOfTwo = [(&'static [u32], f64); 4];

    /// The counts of `pop1`, `i1` and `i2` of the worked example, and what
    /// `vegan` expects of them.
    const POP1: AtADrawOfTwo = [
        (&[3, 1], 1.5),
        (&[3, 1], 1.5),
        (&[1, 1, 1, 1], 2.0),
        (&[4], 1.0),
    ];

    /// The counts of `pop2`, `i3`, `i4` and `i5`, and what `vegan` expects
    /// of them. Variant 5, 4 copies of one allele and 2 of another, is the
    /// one of the eight the spec says a reader cannot do in their head.
    const POP2: AtADrawOfTwo = [
        (&[5], 1.0),
        (&[3], 1.0),
        (&[1, 1, 1, 1], 2.0),
        (&[4, 2], 1.5333333333),
    ];

    /// The alleles a draw of `num_called_alleles` is expected to show at one
    /// variant: one minus the chance of missing an allele, summed over the
    /// counts the population called there, which is the formula of "What it
    /// gives" of "The number of alleles". The module gets it in task 3.2 of
    /// `docs/plans/diversity.md`; here it is what turns the chance into the
    /// numbers `vegan` gave.
    fn alleles_a_draw_shows(counts: &[u32], num_called_alleles: u32) -> f64 {
        let called_alleles: u32 = counts.iter().sum();

        counts
            .iter()
            .map(|count| {
                1.0 - chance_a_draw_misses_an_allele(called_alleles, *count, num_called_alleles)
            })
            .sum()
    }

    /// It checks the four values of one population against `vegan`'s and
    /// their mean, which is the allelic richness of a draw of 2 there.
    fn assert_a_draw_of_two_shows(of_the_pop: AtADrawOfTwo, mean: f64, what: &str) {
        let mut total = 0.0;

        for (counts, vegan) in of_the_pop {
            let found = alleles_a_draw_shows(counts, 2);

            assert!(
                (found - vegan).abs() <= OF_TEN_DIGITS,
                "a draw of 2 of {counts:?} in {what} shows {found} alleles, \
                 and vegan gives {vegan}"
            );
            total += found;
        }

        let found = total / 4.0;

        assert!(
            (found - mean).abs() <= OF_TEN_DIGITS,
            "the mean over the four variants of {what} is {found}, and it is {mean}"
        );
    }

    /// The four variants that count for `pop1`: two where it called 3
    /// copies of one allele and 1 of another, one where each of four alleles
    /// has a single copy, and one where all 4 copies it called are alike. A
    /// draw of 2 shows 1.5, 1.5, 2 and 1 alleles there, 1.5 on average.
    #[test]
    fn a_draw_of_two_shows_one_and_a_half_alleles_in_pop1_of_the_worked_example() {
        assert_a_draw_of_two_shows(POP1, 1.5, "pop1");
    }

    /// The four that count for `pop2`: two where every copy it called is
    /// alike, so a draw of 2 shows 1 allele however it falls, one of four
    /// single copies, and variant 5, of 4 copies of one allele and 2 of
    /// another, which the table gives to ten digits.
    #[test]
    fn a_draw_of_two_shows_1_3833_alleles_in_pop2_of_the_worked_example() {
        assert_a_draw_of_two_shows(POP2, 1.3833333333, "pop2");
    }

    /// Variant 5 of `pop2` said the other way round, which is how the spec
    /// words it: of 4 copies of allele 0 and 2 of allele 1, a draw of 2
    /// shows both alleles with chance 8/15. The allele called 4 times is
    /// called more often than the draw takes and is still missed, by the
    /// one draw of the 15 that takes the two rare copies, so the count of an
    /// allele above the draw size is not a chance of 0.
    #[test]
    fn a_draw_of_two_of_four_and_two_copies_shows_both_alleles_with_chance_eight_fifteenths() {
        let misses_the_common = chance_a_draw_misses_an_allele(6, 4, 2);
        let misses_the_rare = chance_a_draw_misses_an_allele(6, 2, 2);
        let shows_both = 1.0 - misses_the_common - misses_the_rare;

        assert!(
            (shows_both - 8.0 / 15.0).abs() <= OF_A_FEW_ROUNDINGS,
            "a draw of 2 of 4 and 2 copies shows both alleles with chance \
             {shows_both}, which misses 1/15 and 6/15 of the draws, and it \
             is 8/15"
        );
    }

    /// A draw takes an allele it has no room to avoid. With 5 of the 6
    /// copies of one allele only 1 is of another, so every draw of 2 holds
    /// the common one; with the draw as large as the called copies every
    /// draw holds every allele. A draw larger than the called copies is no
    /// draw at all and is the same 0, which is what keeps a caller that
    /// forgot to leave such a variant out of a NaN.
    #[test]
    fn a_draw_that_cannot_avoid_an_allele_misses_it_with_chance_zero() {
        for (called_alleles, count_of_the_allele, num_called_alleles) in
            [(6_u32, 5_u32, 2_u32), (4, 1, 4), (4, 1, 6)]
        {
            let found = chance_a_draw_misses_an_allele(
                called_alleles,
                count_of_the_allele,
                num_called_alleles,
            );

            assert!(
                found == 0.0,
                "a draw of {num_called_alleles} of {called_alleles} copies \
                 misses an allele called {count_of_the_allele} times with \
                 chance {found}, and it is 0"
            );
        }
    }

    /// An allele with no copy called is missed by every draw. No
    /// standardized value asks for it, since all of them walk the alleles
    /// the population called, and what the case pins is that the factors of
    /// a count of 0 are each exactly 1 and their product is 1 and not a
    /// number near it.
    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "every factor of a count of 0 is a number over itself, \
                  which is exactly 1 in an f64, so the product is 1 to the \
                  bit"
    )]
    fn an_allele_the_population_did_not_call_is_missed_by_every_draw() {
        let found = chance_a_draw_misses_an_allele(4, 0, 2);

        assert!(
            found == 1.0,
            "a draw of 2 of 4 copies misses an allele called 0 times with \
             chance {found}, and it is 1"
        );
    }

    /// A million copies of a variant, one of them of this allele: the two
    /// factors are 999999/1000000 and 999998/999999, whose product is
    /// 999998/1000000. A version built from three factorials gives a NaN
    /// here, since a factorial above 170! is an infinity in an `f64` and
    /// this one has a million terms, so this is the case that holds the
    /// function to the product of its factors.
    #[test]
    fn a_million_called_copies_give_a_chance_no_factorial_could() {
        let found = chance_a_draw_misses_an_allele(1_000_000, 1, 2);

        assert!(
            (found - 0.999998).abs() <= OF_A_FEW_ROUNDINGS,
            "a draw of 2 of a million copies misses a single copy with \
             chance {found}, and it is 0.999998"
        );
    }
}

#[cfg(test)]
mod the_private_alleles_of_one_variant {
    use super::{OfAPop, OfAPopAtTheRow, OfTheDrawAtTheRow, add_the_standardized_private_alleles};

    /// Two populations that are both the one diploid individual `0/1` get half
    /// a private allele each at a draw of one allele, where the truth is 0: the
    /// two draws are draws of the same two gene copies and can never give an
    /// allele to one population and not the other, while the estimator popnei
    /// computes reads them as draws of copies of their own. "The cases" of
    /// `docs/specs/diversity.md` works this out, and the twenty-third pair of
    /// `tests/reference/diversity/enumerate_private.tsv` is it: the one pair of
    /// that file whose two ways of computing the value differ, by exactly 1/2,
    /// that pair alone being enumerated over the labelled gene copies of the
    /// individuals and not over allele counts. popnei does not refuse the
    /// overlap, a caller being free to put an individual in more than one
    /// population, so 1/2 is what the estimator has to give here.
    ///
    /// It is asserted on the function that works the value out and not on a
    /// pass, because a pass refuses a `num_called_alleles` below 2 and the draw
    /// of this pair is of one allele. What a pass shows of the same overlap is
    /// in
    /// `a_population_against_a_copy_of_itself_counts_no_private_allele_and_only_a_draw_of_every_copy_shows_none`,
    /// at a draw of 2.
    #[test]
    fn two_populations_of_one_shared_individual_get_half_a_private_allele_each_at_a_draw_of_one() {
        let of_the_shared_individual = {
            let mut of_the_pop = OfAPopAtTheRow::none();
            of_the_pop.counts[0] = 1;
            of_the_pop.counts[1] = 1;
            of_the_pop.one_past_the_largest = 2;
            of_the_pop.called_alleles = 2;
            of_the_pop
        };
        let of_each_pop = vec![of_the_shared_individual.clone(), of_the_shared_individual];
        let mut at_the_row = OfTheDrawAtTheRow::of(2);
        let mut totals = vec![OfAPop::none(); 2];

        add_the_standardized_private_alleles(&of_each_pop, 1, &mut at_the_row, &mut totals);

        for (pop, counted) in totals.iter().enumerate() {
            let found = counted.sum_private_alleles_in_draw;

            assert!(
                (found - 0.5).abs() <= f64::EPSILON,
                "a draw of one allele shows {found} private alleles in the population {pop} of \
                 one shared individual, and it shows 0.5"
            );
        }
    }
}
