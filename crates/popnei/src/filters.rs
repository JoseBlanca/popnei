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
//! through [`BlockReader::filtering_stats`].
//!
//! [`LdFilter`] is the fourth filter and compares no number of a variant
//! alone: it takes out the variants whose r², the measure of linkage
//! disequilibrium of `docs/specs/ld.md`, is above a threshold against a
//! variant it kept no more than `max_dist` base pairs behind them on their
//! chromosome. It is a type of its own because it holds those kept variants,
//! its window, between one block and the next.
//!
//! A filter of a pass over the variants is a reader over another reader,
//! [`FilteredReader`] for the three that compare one number of a variant
//! and [`LdFilteredReader`] for the fourth, and several filters are several
//! of them, one over the other, in the order in which the user put them on.
//! [`chain_of`] builds that chain from the criteria of one pass, and it is
//! what each binding crate calls when a pass starts.
//!
//! `docs/specs/filters.md` has the design, and the row `filters` of section
//! 9 of `docs/architecture.md` where the module sits.

use std::fmt;

use crate::block::{Block, BlockReader};
use crate::error::{Error, Result};
use crate::ld::{LdDosages, r2_between};
use crate::variant::{
    AlleleCounts, ChromTable, Needs, count_alleles, count_gts, the_major_allele_frequency,
};

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

/// What a filter compares with a threshold, with the largest value that
/// keeps the variant.
///
/// The first three compare one number of the variant alone, and the fourth
/// compares the variant with the variants kept before it. A variant stays
/// when the value is at most the threshold, so one whose value is exactly
/// the threshold stays. Each of the four values is one count divided by
/// another, or an r², which is a correlation squared, so a threshold is a
/// number from 0 to 1. A variant that has no value, one with no called
/// allele for the major allele frequency and one with no called genotype
/// for the observed heterozygosity, is not kept, whatever the threshold.
///
/// The first three are filtered by a [`VarFilter`] under a
/// [`FilteredReader`], and the fourth by an [`LdFilter`] under an
/// [`LdFilteredReader`], which [`chain_of`] is what builds for each.
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
    /// The largest r² a variant may have against a variant kept no more
    /// than `max_dist` base pairs behind it on its chromosome, which is the
    /// window of that variant. r² is the squared correlation between the
    /// dosages of two variants, `docs/specs/ld.md`, and a variant whose
    /// called genotypes hold one dosage is dropped at every threshold,
    /// having nothing to tell another variant apart with.
    MaxLdR2 {
        /// The largest r² against a variant of the window that keeps the
        /// variant, a number from 0 to 1.
        max_allowed_r2: f64,
        /// How many base pairs behind a variant its window reaches, 1 or
        /// more.
        max_dist: u64,
    },
}

impl VarFilteringCriterion {
    /// `"missing_data"`, `"maf"`, `"obs_het"` or `"ld"`: the name under
    /// which the counts of the filter reach a Python or a TypeScript user,
    /// and the name by which a chain of readers is asked whether it holds a
    /// filter of this kind already.
    #[must_use]
    pub fn kind(&self) -> &'static str {
        match self {
            VarFilteringCriterion::MaxMissingRate(_) => "missing_data",
            VarFilteringCriterion::MaxMaf(_) => "maf",
            VarFilteringCriterion::MaxObsHet(_) => "obs_het",
            VarFilteringCriterion::MaxLdR2 { .. } => THE_KIND_OF_THE_LD_FILTER,
        }
    }

    /// The largest value that keeps the variant, whichever of the four
    /// values this criterion compares, which for
    /// [`MaxLdR2`](VarFilteringCriterion::MaxLdR2) is its
    /// `max_allowed_r2`.
    ///
    /// A binding crate reads it for the arguments of the step it shows the
    /// user, `{"max_allowed_maf": 0.95}`.
    #[must_use]
    pub fn threshold(&self) -> f64 {
        match self {
            VarFilteringCriterion::MaxMissingRate(threshold)
            | VarFilteringCriterion::MaxMaf(threshold)
            | VarFilteringCriterion::MaxObsHet(threshold)
            | VarFilteringCriterion::MaxLdR2 {
                max_allowed_r2: threshold,
                ..
            } => *threshold,
        }
    }

    /// The `max_dist` of [`MaxLdR2`](VarFilteringCriterion::MaxLdR2), the
    /// second value of the arguments of its step, and `None` for the three
    /// criteria that compare one number of a variant alone.
    ///
    /// A binding crate reads it for the step it shows the user,
    /// `{"max_allowed_r2": 0.1, "max_dist": 10000}`, where the threshold is
    /// the first of the two.
    #[must_use]
    pub fn max_dist(&self) -> Option<u64> {
        match self {
            VarFilteringCriterion::MaxMissingRate(_)
            | VarFilteringCriterion::MaxMaf(_)
            | VarFilteringCriterion::MaxObsHet(_) => None,
            VarFilteringCriterion::MaxLdR2 { max_dist, .. } => Some(*max_dist),
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
    /// criterion and the value. And when `criterion` is
    /// [`MaxLdR2`](VarFilteringCriterion::MaxLdR2), which holds the window
    /// of the variants kept behind a variant and is filtered by an
    /// [`LdFilter`]: no call from Python or from TypeScript reaches it,
    /// since [`chain_of`] builds an [`LdFilteredReader`] for that
    /// criterion.
    pub fn new(criterion: VarFilteringCriterion) -> Result<VarFilter> {
        if let VarFilteringCriterion::MaxLdR2 { .. } = criterion {
            return Err(Error::VarFilterOfTheLdCriterion);
        }
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
    /// which [`Block::check`] finds; when the block has variants and no
    /// genotypes, which is the error of a field that is not in the block,
    /// and which a block of no individual or of the ploidy 0 gives too,
    /// since it holds no genotype; when a block of the individuals and the
    /// ploidy it states is more memory than this machine addresses; and
    /// what the counts of one variant refuse, a row that is not a whole
    /// number of genotypes of the ploidy, a variant of more alleles than a
    /// count of them holds, and an allele below
    /// [`MISSING_ALLELE`](crate::variant::MISSING_ALLELE). After any of
    /// them the block is as it was and nothing was added to the counts.
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
        // `check` passed and the genotypes are not empty, so they are the
        // variants of the block times this number and it is one allele at
        // least: the rows are cut by it, and a cut of 0 is what the
        // standard library refuses with a panic.
        let alleles_per_var = block.alleles_per_var()?.max(1);
        let keep = keep_of_the_rows(
            self.criterion,
            &block.gts,
            alleles_per_var,
            block.num_individuals,
            block.ploidy,
        )?;
        let kept = keep.iter().filter(|keep_it| **keep_it).count();
        block.retain_vars(&keep)?;
        // A `usize` is 64 bits on the targets popnei builds natively for
        // and 32 in wasm, so every one of them is a `u64` and neither
        // conversion takes the value it saturates at.
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

/// The name under which the counts of the filter by linkage disequilibrium
/// reach a Python or a TypeScript user, and the name by which a chain of
/// readers is asked whether it holds one already.
const THE_KIND_OF_THE_LD_FILTER: &str = "ld";

/// How many variants the filter by linkage disequilibrium settles at a
/// time.
///
/// Whether a variant is kept decides what the variants after it are
/// compared with, so the variants of a block cannot all be settled at once.
/// What one set of the products of `docs/specs/ld.md` gives is the r² of a
/// set of variants against every variant kept before that set began, and
/// this is how many variants such a set holds: the r² of 256 variants
/// against each other is 65536 values, 512 KB, and `docs/specs/ld.md`
/// measures one product of that shape at 1.9 ms over 1000 individuals.
///
/// It changes no result, which
/// `the_variants_kept_do_not_change_with_the_variants_settled_at_a_time`
/// asserts at 1, at 3 and at 256: the rule reads the positions of the
/// variants and never the end of a set or of a block, and the six sums of a
/// pair are whole numbers that an `f64` holds exactly, so a pair has the
/// same r² in whichever set it is worked out.
///
/// It is of the crate and not of its users, as the tile of the matrix of
/// `docs/specs/ld.md` is private to its module: nothing outside this file
/// reads it.
pub(crate) const THE_VARS_SETTLED_AT_A_TIME: usize = 256;

/// How a variant given to [`LdFilter::filter_block`] does not come after
/// the variant before it.
///
/// The window of a variant is the variants kept behind it on its
/// chromosome, so the filter by linkage disequilibrium is the one reader of
/// popnei that needs the variants of each chromosome to come together and
/// in the order of their positions. Which chromosome it is on is
/// [`TheChromOfTheVariant`], beside this in the error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TheOrderOfTheVariants {
    /// The position of the variant is below the position of the variant
    /// before it, which is on its chromosome.
    ThePositionFalls {
        /// The position of the variant.
        pos: u64,
        /// The position of the variant before it.
        pos_before: u64,
    },
    /// The variant is on a chromosome that had already ended: a variant of
    /// another chromosome came between it and the last variant of its own.
    TheChromosomeCameBack {
        /// The position of the variant.
        pos: u64,
        /// The position of the variant before it, which is on another
        /// chromosome.
        pos_before: u64,
    },
}

impl fmt::Display for TheOrderOfTheVariants {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::ThePositionFalls { pos, pos_before } => write!(
                formatter,
                "it is at the position {pos} of its chromosome and the variant before it at the position {pos_before} of the same chromosome"
            ),
            Self::TheChromosomeCameBack { pos, pos_before } => write!(
                formatter,
                "it is at the position {pos} of a chromosome that had already ended, and the variant before it at the position {pos_before} of another chromosome"
            ),
        }
    }
}

/// The chromosome of a variant that the filter by linkage disequilibrium
/// refused for not coming after the variant before it.
///
/// A block holds the number each of its chromosomes has in the table of the
/// reader that gave it and not its name, so [`LdFilter::filter_block`],
/// which is given a block and nothing else, has the number; an
/// [`LdFilteredReader`] has the table of its source and puts the name in
/// its place, which is what a user looks for in their file. On an assembly
/// of ten thousand scaffolds there is no finding the row otherwise.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TheChromOfTheVariant {
    /// The name of the chromosome in the table of the reader of the pass.
    Named(String),
    /// The number it has in that table, which is what a filter given a
    /// block on its own has of it.
    Numbered(u32),
}

impl fmt::Display for TheChromOfTheVariant {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Named(name) => write!(formatter, "the chromosome {name}"),
            Self::Numbered(chrom) => write!(
                formatter,
                "the chromosome numbered {chrom} in the table of the reader"
            ),
        }
    }
}

/// The filter that takes out the variants that repeat what a variant kept
/// near them on their chromosome already said.
///
/// Two variants say the same thing when their r² is high, and
/// `docs/specs/ld.md` is where popnei works it out. A variant is kept when
/// its called genotypes hold two dosages at least and its r² against every
/// variant of its window is at most `max_allowed_r2`; the window of a
/// variant is the variants the filter has already kept that are on its
/// chromosome and no more than `max_dist` base pairs behind it. A pair
/// whose r² is not defined does not drop the candidate, so the first
/// variant of each chromosome whose called genotypes hold two dosages is
/// always kept, and a variant whose called genotypes hold one dosage is
/// always dropped, having nothing to tell any other variant apart with.
///
/// It is a type of its own and not a [`VarFilter`]: it holds the window
/// between one block and the next, where a `VarFilter` reads each block on
/// its own and keeps nothing but its two counts. The window holds, for each
/// variant of it, its genotypes, its chromosome and its position, one byte
/// for each allele: 500 KB for 250 kept variants of 1000 diploid
/// individuals. The three matrices of
/// [`LdDosages`](crate::ld::LdDosages), 24 bytes for each individual and
/// variant, are built over the whole window when a set of candidates
/// arrives and are given back when it has been settled, so the 6 MB of
/// those 250 variants of 1000 individuals is held while a block is
/// filtered and not between two blocks.
///
/// The dosages are read over every individual of the dataset. A user who
/// wants them read over one population puts the filter of individuals
/// before this one.
pub struct LdFilter {
    /// The largest r² a variant may have against a variant of its window.
    max_allowed_r2: f64,
    /// How many base pairs behind a variant its window reaches.
    max_dist: u64,
    /// How many variants it was given and how many it kept.
    stats: FilteringStats,
    /// The variants kept that are still within `max_dist` of the last
    /// variant read, in the order they were kept, which is the order of
    /// their positions.
    window: TheWindow,
    /// The chromosome and the position of the last variant read, and
    /// `None` before the first block.
    before: Option<TheVariantBefore>,
    /// The number of every chromosome the filter has read a variant of,
    /// which is what says that a chromosome has come back: one value for
    /// each chromosome of the table of the reader, true for the ones a
    /// variant has been read of.
    ///
    /// The numbers of a [`ChromTable`] are dense, so whether a chromosome
    /// has been read is one index into this and not a search: a filter that
    /// searched a list of the chromosomes read, and copied it for every
    /// block, took 0.11 s over a source of 40000 chromosomes and 2.47 s
    /// over one of 320000, which grows with the square of them, where a
    /// fragmented assembly has millions.
    chroms_read: Vec<bool>,
}

/// The chromosome and the position of the variant the filter read last,
/// which say whether the variant after it comes after it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TheVariantBefore {
    /// The number of its chromosome in the table of the reader.
    chrom: u32,
    /// Its position.
    pos: u64,
}

/// Where the filter has got to in its source once a block has been read
/// through, which it takes over when the block is kept.
struct TheOrderRead {
    /// The last variant of the block.
    before: Option<TheVariantBefore>,
    /// The chromosomes of the block that the filter had not read before and
    /// marked as read: a block that is refused further on unmarks them, so
    /// that the filter is as it was.
    marked: Vec<u32>,
}

impl LdFilter {
    /// The filter that keeps the variants whose r² against every variant of
    /// their window is at most `max_allowed_r2`, with an empty window and
    /// both its counts at 0.
    ///
    /// `max_dist` is how many base pairs behind a variant its window
    /// reaches, and two variants at one position are 0 apart, so each is in
    /// the window of the other.
    ///
    /// # Errors
    ///
    /// When `max_allowed_r2` is NaN, below 0 or above 1, and when
    /// `max_dist` is below 1: the error names the argument and the value.
    pub fn new(max_allowed_r2: f64, max_dist: u64) -> Result<LdFilter> {
        // A NaN is in no range, so this one comparison refuses the three
        // thresholds that are not a number from 0 to 1.
        if !(0.0..=1.0).contains(&max_allowed_r2) {
            return Err(Error::VarFilterThresholdOutOfRange {
                kind: THE_KIND_OF_THE_LD_FILTER,
                threshold: max_allowed_r2,
            });
        }
        if max_dist == 0 {
            return Err(Error::LdFilterMaxDistTooSmall { max_dist });
        }
        Ok(LdFilter {
            max_allowed_r2,
            max_dist,
            stats: FilteringStats::default(),
            window: TheWindow::default(),
            before: None,
            chroms_read: Vec::new(),
        })
    }

    /// What it filters by, with both its arguments: the criterion a binding
    /// crate reads the kind and the arguments of the step from.
    #[must_use]
    pub fn criterion(&self) -> VarFilteringCriterion {
        VarFilteringCriterion::MaxLdR2 {
            max_allowed_r2: self.max_allowed_r2,
            max_dist: self.max_dist,
        }
    }

    /// The largest r² a variant may have against a variant of its window.
    #[must_use]
    pub fn max_allowed_r2(&self) -> f64 {
        self.max_allowed_r2
    }

    /// How many base pairs behind a variant its window reaches.
    #[must_use]
    pub fn max_dist(&self) -> u64 {
        self.max_dist
    }

    /// The variants of the block that pass, kept in it in their order, and
    /// the others dropped: the genotypes and every column of the block are
    /// compacted in place.
    ///
    /// The window carries over, so the variants of a block are compared
    /// with the variants kept in the blocks before it, and so does the last
    /// variant read, which the first variant of the block has to come
    /// after. The variants of the block are added to the counts, and the
    /// ones that stayed to the ones kept. A block of no variants is left as
    /// it is.
    ///
    /// # Errors
    ///
    /// When the arrays of the block are not of the size the block states,
    /// which [`Block::check`] finds; when the block has variants and no
    /// genotypes or no position, which is the error of a field that is not
    /// in the block; when a variant of the block does not come after the
    /// one before it, its position falling within its chromosome or its
    /// chromosome having already ended; and what the dosages of the block
    /// refuse, which the `# Errors` of
    /// [`LdDosages::of_block`](crate::ld::LdDosages::of_block) lists, a
    /// ploidy above 255 and a block this machine has not the memory of the
    /// three matrices for among them. After any of them the block is as it
    /// was, nothing was added to the counts and the window is as it was.
    pub fn filter_block(&mut self, block: &mut Block) -> Result<()> {
        self.the_block_filtered(block, THE_VARS_SETTLED_AT_A_TIME)
    }

    /// The block filtered with `at_a_time` variants settled at a time,
    /// which [`LdFilter::filter_block`] is with
    /// [`THE_VARS_SETTLED_AT_A_TIME`] and the tests are what give another
    /// number.
    ///
    /// The variants kept are the same whatever the number, which
    /// `the_variants_kept_do_not_change_with_the_variants_settled_at_a_time`
    /// asserts at 1, at 3 and at 256: the rule reads the positions of the
    /// variants and never the end of a set, and a pair has the same r² in
    /// whichever set it is worked out.
    ///
    /// # Errors
    ///
    /// Those of [`LdFilter::filter_block`].
    fn the_block_filtered(&mut self, block: &mut Block, at_a_time: usize) -> Result<()> {
        // The rows are cut out of the genotypes by the sizes the block
        // states, so those sizes are checked before anything is read.
        block.check()?;
        let processed = block.num_vars;
        if processed == 0 {
            return Ok(());
        }
        let missing = (Needs::GTS | Needs::CHROM_POS).difference(block.fields());
        if !missing.is_empty() {
            return Err(Error::FieldsNotInTheBlock { fields: missing });
        }
        let (Some(chroms), Some(positions)) = (block.chrom.as_deref(), block.pos.as_deref()) else {
            // A block holds the chromosome and the position as one field
            // and only when both columns are there, so the check above is
            // what refuses a block without them and this is not reached.
            return Err(Error::FieldsNotInTheBlock {
                fields: Needs::CHROM_POS,
            });
        };
        // The order of the variants is read before the dosages are built: a
        // source whose positions do not rise is refused whatever else the
        // block holds.
        let order = self.the_order_read(chroms, positions)?;
        let settled = the_variants_that_stay(self, block, chroms, positions, at_a_time)
            .and_then(|settled| block.retain_vars(&settled.keep).map(|()| settled));
        let settled = match settled {
            Ok(settled) => settled,
            Err(error) => {
                // The block is refused, so the chromosomes it was the first
                // to hold a variant of are unmarked and the filter is as it
                // was.
                self.the_chromosomes_unmarked(&order.marked);
                return Err(error);
            }
        };
        // Nothing of the filter has changed until here but the chromosomes
        // read, so an error above left the window, the counts and the place
        // in the source as they were: the block was settled against a
        // window of its own.
        self.window = settled.window;
        self.before = order.before;
        // A `usize` is 64 bits on the targets popnei builds natively for
        // and 32 in wasm, so every one of them is a `u64` and neither
        // conversion takes the value it saturates at.
        self.stats.vars_processed = self
            .stats
            .vars_processed
            .saturating_add(u64::try_from(processed).unwrap_or(u64::MAX));
        self.stats.vars_kept = self
            .stats
            .vars_kept
            .saturating_add(u64::try_from(settled.kept).unwrap_or(u64::MAX));
        Ok(())
    }

    /// How many variants it was given and how many it kept, over every
    /// block it has taken since it was built.
    #[must_use]
    pub fn stats(&self) -> FilteringStats {
        self.stats
    }

    /// The variant the block ends at once the filter has read the block,
    /// with every variant of it checked against the one before it and the
    /// chromosomes of the block marked as read.
    ///
    /// What it marked comes back with it, so that a block refused further
    /// on unmarks those chromosomes and leaves the filter as it was; a
    /// variant of the block that is refused here unmarks them itself.
    ///
    /// # Errors
    ///
    /// [`Error::LdFilterVariantOutOfOrder`] at the first variant of the
    /// block whose position falls below the position of the variant before
    /// it on its chromosome, or whose chromosome had already ended, and
    /// [`Error::LdNoMemory`] when this machine does not give the memory of
    /// the chromosomes read.
    fn the_order_read(&mut self, chroms: &[u32], positions: &[u64]) -> Result<TheOrderRead> {
        let mut marked: Vec<u32> = Vec::new();
        match self.the_order_walked(chroms, positions, &mut marked) {
            Ok(before) => Ok(TheOrderRead { before, marked }),
            Err(error) => {
                self.the_chromosomes_unmarked(&marked);
                Err(error)
            }
        }
    }

    /// The variant the block ends at, with every variant of it checked
    /// against the one before it, and the chromosomes the block was the
    /// first to hold a variant of written into `marked` and marked as read.
    ///
    /// # Errors
    ///
    /// Those of [`LdFilter::the_order_read`], which unmarks what this
    /// marked before it gives one on.
    fn the_order_walked(
        &mut self,
        chroms: &[u32],
        positions: &[u64],
        marked: &mut Vec<u32>,
    ) -> Result<Option<TheVariantBefore>> {
        let mut before = self.before;
        for (variant, (chrom, pos)) in chroms.iter().zip(positions).enumerate() {
            let out_of_order = match before {
                Some(before) if before.chrom == *chrom => {
                    (*pos < before.pos).then_some(TheOrderOfTheVariants::ThePositionFalls {
                        pos: *pos,
                        pos_before: before.pos,
                    })
                }
                Some(before) => {
                    self.has_read(*chrom)
                        .then_some(TheOrderOfTheVariants::TheChromosomeCameBack {
                            pos: *pos,
                            pos_before: before.pos,
                        })
                }
                // The first variant the filter reads comes after nothing.
                None => None,
            };
            if let Some(problem) = out_of_order {
                return Err(Error::LdFilterVariantOutOfOrder {
                    variant: self
                        .stats
                        .vars_processed
                        .saturating_add(u64::try_from(variant).unwrap_or(u64::MAX))
                        .saturating_add(1),
                    chrom: TheChromOfTheVariant::Numbered(*chrom),
                    problem,
                });
            }
            if before.is_none_or(|before| before.chrom != *chrom) {
                self.mark_as_read(*chrom, marked)?;
            }
            before = Some(TheVariantBefore {
                chrom: *chrom,
                pos: *pos,
            });
        }
        Ok(before)
    }

    /// Whether the filter has read a variant of that chromosome.
    fn has_read(&self, chrom: u32) -> bool {
        the_place_of(chrom)
            .and_then(|of_it| self.chroms_read.get(of_it))
            .copied()
            .unwrap_or(false)
    }

    /// Marks a chromosome the filter has read a variant of, and writes it
    /// into `marked` when it was not marked already.
    ///
    /// # Errors
    ///
    /// [`Error::LdNoMemory`] when this machine does not give the memory of
    /// one value for each chromosome up to that one.
    fn mark_as_read(&mut self, chrom: u32, marked: &mut Vec<u32>) -> Result<()> {
        // A `usize` is 64 bits natively and 32 in WebAssembly, so the
        // number of a chromosome is one; on a machine whose `usize` is
        // narrower the vector is asked for a length no machine gives and
        // the memory is what answers.
        let of_it = the_place_of(chrom).unwrap_or(usize::MAX);
        let chroms = of_it.saturating_add(1);
        if self.chroms_read.len() < chroms {
            let more = chroms.saturating_sub(self.chroms_read.len());
            the_room_for(&mut self.chroms_read, more, "the chromosomes read")?;
            self.chroms_read.resize(chroms, false);
        }
        let Some(read) = self.chroms_read.get_mut(of_it) else {
            // The vector was grown to hold that chromosome, so this is not
            // reached.
            return Err(Error::LdNoMemory {
                what: "the chromosomes read",
                values: chroms,
            });
        };
        if !*read {
            *read = true;
            the_room_for(marked, 1, "the chromosomes of the block")?;
            marked.push(chrom);
        }
        Ok(())
    }

    /// Unmarks the chromosomes a block marked as read, which is what a
    /// block that is refused leaves behind.
    fn the_chromosomes_unmarked(&mut self, marked: &[u32]) {
        for chrom in marked {
            if let Some(read) =
                the_place_of(*chrom).and_then(|of_it| self.chroms_read.get_mut(of_it))
            {
                *read = false;
            }
        }
    }
}

/// Where the chromosome of that number is in the chromosomes read of an
/// [`LdFilter`], and `None` on a machine whose `usize` does not hold a
/// `u32`.
fn the_place_of(chrom: u32) -> Option<usize> {
    usize::try_from(chrom).ok()
}

impl fmt::Debug for LdFilter {
    /// What it filters and where it has got to. The dosages of the window
    /// are left out and how many variants it holds is given instead: they
    /// are three numbers for each individual of each variant of it.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("LdFilter")
            .field("max_allowed_r2", &self.max_allowed_r2)
            .field("max_dist", &self.max_dist)
            .field("stats", &self.stats)
            .field("vars_in_the_window", &self.window.num_vars())
            .field("before", &self.before)
            .finish_non_exhaustive()
    }
}

/// Which variants of a block the filter by linkage disequilibrium keeps,
/// with the window it leaves behind.
struct TheBlockSettled {
    /// One value for each variant of the block, in its order.
    keep: Vec<bool>,
    /// How many of them stayed.
    kept: usize,
    /// The window once every variant of the block has been settled: the
    /// variants the filter had kept that the last variant of the block has
    /// not left behind, with the ones the block added to them.
    window: TheWindow,
}

/// The variants a candidate is compared with: the ones the filter has kept
/// that are on the candidate's chromosome and no more than `max_dist` base
/// pairs behind it, in the order of their positions.
///
/// It holds the genotypes of those variants and not their dosages, because
/// the whole window is one operand of the products of `docs/specs/ld.md`
/// when a set of candidates arrives: [`TheWindow::dosages`] builds the
/// three matrices of every variant of it together, so the r² of a set of
/// candidates against the window is one call of
/// [`r2_between`](crate::ld::r2_between) and not one call for each variant
/// kept. The genotypes are one byte for each allele where the three
/// matrices are 24 bytes for each individual, so a diploid window holds
/// between two blocks a twelfth of what they would be.
#[derive(Debug, Default)]
struct TheWindow {
    /// The number of the chromosome of each variant of it, in the table of
    /// the reader that gave the block it came in.
    chroms: Vec<u32>,
    /// The position of each of them, 1 based as in a VCF.
    poss: Vec<u64>,
    /// The genotypes of each of them over every individual of the dataset,
    /// one variant after another, as a block holds them.
    gts: Vec<i8>,
    /// How many individuals the block the variants of the window came in
    /// had, which their dosages are built over.
    num_individuals: usize,
    /// How many alleles the genotype of one individual holds in that same
    /// block.
    ploidy: usize,
}

impl TheWindow {
    /// How many variants it holds.
    fn num_vars(&self) -> usize {
        self.poss.len()
    }

    /// How many alleles the genotypes of one variant of it hold.
    fn alleles_per_var(&self) -> usize {
        // The genotypes of a block of this many individuals of this ploidy
        // are in memory, so their product is a number this machine counted.
        self.num_individuals.saturating_mul(self.ploidy)
    }

    /// A window of the same variants, which a block is settled against and
    /// which the filter takes over once the block has been settled, so that
    /// a block that is refused leaves the window of the filter as it was.
    ///
    /// # Errors
    ///
    /// [`Error::LdNoMemory`] when this machine does not give the memory of
    /// the genotypes of the window or of its two columns.
    fn copy_of(&self) -> Result<TheWindow> {
        Ok(TheWindow {
            chroms: the_copy_of(&self.chroms, "the chromosomes of the window")?,
            poss: the_copy_of(&self.poss, "the positions of the window")?,
            gts: the_copy_of(&self.gts, "the genotypes of the window")?,
            num_individuals: self.num_individuals,
            ploidy: self.ploidy,
        })
    }

    /// Takes the individuals and the ploidy of `block`, or refuses a block
    /// that does not hold the dataset the variants of the window came from.
    ///
    /// The window keeps the genotypes of its variants and reads them as the
    /// individuals and the ploidy of the source, so a block of other
    /// individuals or of another ploidy would give the window dosages read
    /// off the wrong alleles. It is the error of `docs/specs/block.md` for
    /// blocks of one source that do not hold the same dataset, which is a
    /// defect of the reader and not of what a user wrote.
    ///
    /// # Errors
    ///
    /// [`Error::BlocksDoNotFitTogether`] when the block holds another
    /// number of individuals or another ploidy than the variants of the
    /// window were read in.
    fn takes_the_block(&mut self, block: &Block) -> Result<()> {
        if self.num_vars() == 0 {
            self.num_individuals = block.num_individuals;
            self.ploidy = block.ploidy;
            return Ok(());
        }
        if self.num_individuals != block.num_individuals || self.ploidy != block.ploidy {
            return Err(Error::BlocksDoNotFitTogether {
                num_individuals: self.num_individuals,
                ploidy: self.ploidy,
                found_num_individuals: block.num_individuals,
                found_ploidy: block.ploidy,
            });
        }
        Ok(())
    }

    /// The three matrices of `docs/specs/ld.md` over every variant of the
    /// window, which one call of [`r2_between`](crate::ld::r2_between) then
    /// reads a whole set of candidates against.
    ///
    /// The genotypes are lent to the dosages and taken back with the memory
    /// they have, so the window is as it was and nothing of it is copied.
    ///
    /// # Errors
    ///
    /// What [`LdDosages::of_block`](crate::ld::LdDosages::of_block)
    /// refuses, a window this machine has not the memory of the three
    /// matrices for among them.
    fn dosages(&mut self) -> Result<LdDosages> {
        let num_vars = self.num_vars();
        the_dosages_of(&mut self.gts, num_vars, self.num_individuals, self.ploidy)
    }

    /// Takes off the variants that a variant at `pos` of the chromosome
    /// `chrom` has left behind: the ones on another chromosome and the ones
    /// more than `max_dist` base pairs behind it. They are the first of the
    /// window, whose variants are in the order of their positions.
    fn leave_behind(&mut self, chrom: u32, pos: u64, max_dist: u64) {
        let left = self
            .chroms
            .iter()
            .zip(&self.poss)
            .take_while(|(of_it, at_it)| {
                **of_it != chrom || pos.checked_sub(**at_it).is_none_or(|dist| dist > max_dist)
            })
            .count();
        let alleles = left
            .saturating_mul(self.alleles_per_var())
            .min(self.gts.len());
        self.chroms.drain(..left);
        self.poss.drain(..left);
        self.gts.drain(..alleles);
    }

    /// Adds a variant the filter has kept, which is ahead of every variant
    /// of the window, with its genotypes over every individual of the
    /// dataset.
    ///
    /// # Errors
    ///
    /// [`Error::LdNoMemory`] when this machine does not give the memory of
    /// one more variant of the window.
    fn push(&mut self, chrom: u32, pos: u64, gts: &[i8]) -> Result<()> {
        the_room_for(&mut self.chroms, 1, "the chromosomes of the window")?;
        the_room_for(&mut self.poss, 1, "the positions of the window")?;
        the_room_for(&mut self.gts, gts.len(), "the genotypes of the window")?;
        self.chroms.push(chrom);
        self.poss.push(pos);
        self.gts.extend_from_slice(gts);
        Ok(())
    }
}

/// The dosages of the `num_vars` variants whose genotypes are in `gts`,
/// over the `num_individuals` individuals of the ploidy `ploidy`.
///
/// The buffer is given back as it was, with the memory it has and with the
/// genotypes in it, whether the dosages were built or not: the sets of a
/// block share one allocation, and the window lends its own genotypes and
/// keeps them.
///
/// # Errors
///
/// What [`LdDosages::of_block`](crate::ld::LdDosages::of_block) refuses,
/// which its `# Errors` lists, a set of variants this machine has not the
/// memory of the three matrices for among them.
fn the_dosages_of(
    gts: &mut Vec<i8>,
    num_vars: usize,
    num_individuals: usize,
    ploidy: usize,
) -> Result<LdDosages> {
    let mut block = Block {
        num_vars,
        num_individuals,
        ploidy,
        gts: std::mem::take(gts),
        chrom: None,
        pos: None,
        id: None,
        alleles: None,
        qual: None,
    };
    let dosages = LdDosages::of_block(&block, &[]);
    *gts = std::mem::take(&mut block.gts);
    dosages
}

/// Asks this machine for room for `more` values in `vector`, which `what`
/// names.
///
/// # Errors
///
/// [`Error::LdNoMemory`] when the machine does not give it. The memory is
/// asked for with `try_reserve`, which gives it back as an error where a
/// `push` or an `extend` would end the process, as `docs/specs/filters.md`
/// asks for a window this machine has not the memory of.
fn the_room_for<T>(vector: &mut Vec<T>, more: usize, what: &'static str) -> Result<()> {
    vector
        .try_reserve(more)
        .map_err(|_| Error::LdNoMemory { what, values: more })
}

/// A copy of `values`, which `what` names, with its memory asked of this
/// machine and not taken.
///
/// # Errors
///
/// [`Error::LdNoMemory`] when the machine does not give it.
fn the_copy_of<T: Copy>(values: &[T], what: &'static str) -> Result<Vec<T>> {
    let mut copy: Vec<T> = Vec::new();
    copy.try_reserve_exact(values.len())
        .map_err(|_| Error::LdNoMemory {
            what,
            values: values.len(),
        })?;
    copy.extend_from_slice(values);
    Ok(copy)
}

/// A buffer of `values` values of `value`, which `what` names, with its
/// memory asked of this machine and not taken.
///
/// # Errors
///
/// [`Error::LdNoMemory`] when the machine does not give it.
fn the_buffer_of<T: Clone>(value: T, values: usize, what: &'static str) -> Result<Vec<T>> {
    let mut buffer: Vec<T> = Vec::new();
    buffer
        .try_reserve_exact(values)
        .map_err(|_| Error::LdNoMemory { what, values })?;
    buffer.resize(values, value);
    Ok(buffer)
}

/// Grows `buffer` to `values` values of 0, and leaves it as it is when it
/// holds that many already.
///
/// # Errors
///
/// [`Error::LdNoMemory`] when this machine does not give the memory of the
/// values it has not.
fn the_buffer_grown_to(buffer: &mut Vec<f64>, values: usize, what: &'static str) -> Result<()> {
    let more = values.saturating_sub(buffer.len());
    if more == 0 {
        return Ok(());
    }
    the_room_for(buffer, more, what)?;
    buffer.resize(values, 0.0);
    Ok(())
}

/// Which variants of the block the filter keeps, with the window it leaves
/// behind.
///
/// The variants are settled `at_a_time` at a time. The r² of a whole set of
/// them against the whole window is one call, over one set of the products
/// of `docs/specs/ld.md`, and the r² of the variants of the set against one
/// another is one more: what cannot be done that way is a candidate against
/// the variants kept inside its own set, since whether one of them is kept
/// decides what the next one is compared with, so those are read out of the
/// r² of the set in the order of the variants.
///
/// `at_a_time` changes no result, which
/// `the_variants_kept_do_not_change_with_the_variants_settled_at_a_time`
/// asserts at 1, at 3 and at 256: the rule reads the positions of the
/// variants and never the end of a set or of a block, and the six sums of a
/// pair are whole numbers that an `f64` holds exactly, so a pair has the
/// same r² in whichever set it is worked out.
/// [`THE_VARS_SETTLED_AT_A_TIME`] is what the filter passes.
///
/// `chroms` and `positions` are the columns of the block, which
/// [`Block::check`] has found to hold one value for each of its variants.
///
/// # Errors
///
/// [`Error::BlocksDoNotFitTogether`] when the block holds other individuals
/// or another ploidy than the variants of the window were read in,
/// [`Error::LdNoMemory`] when this machine does not give the memory of the
/// window, of the genotypes of a set or of the r² of one, and what
/// [`LdDosages::of_block`](crate::ld::LdDosages::of_block) and
/// [`r2_between`](crate::ld::r2_between) refuse.
fn the_variants_that_stay(
    filter: &LdFilter,
    block: &Block,
    chroms: &[u32],
    positions: &[u64],
    at_a_time: usize,
) -> Result<TheBlockSettled> {
    let alleles_per_var = block.alleles_per_var()?;
    // The block is settled against a window of its own, which the filter
    // takes over once every variant of the block has been settled, so a
    // block that is refused leaves the window of the filter as it was.
    let mut window = filter.window.copy_of()?;
    window.takes_the_block(block)?;
    let at_a_time = at_a_time.min(block.num_vars).max(1);
    // The buffers of a set, asked for once and written over by every set of
    // the block: the genotypes of a set are at most the genotypes of the
    // block, which are in memory, and the r² of the pairs of a set of 256
    // variants is 65536 values, 512 KB. A count that saturates here is more
    // memory than any machine gives, so it comes back as the error of the
    // memory and never as a buffer of the wrong size.
    let mut keep: Vec<bool> = Vec::new();
    the_room_for(
        &mut keep,
        block.num_vars,
        "the variants of the block that stay",
    )?;
    let mut gts_of_the_set: Vec<i8> = Vec::new();
    the_room_for(
        &mut gts_of_the_set,
        at_a_time.saturating_mul(alleles_per_var),
        "the genotypes of a set of candidates",
    )?;
    let mut dropped = the_buffer_of(false, at_a_time, "the candidates of a set that are dropped")?;
    let mut r2_of_the_set = the_buffer_of(
        0.0,
        at_a_time.saturating_mul(at_a_time),
        "the r² of the pairs of a set",
    )?;
    let mut r2_against_the_window: Vec<f64> = Vec::new();
    let mut first: usize = 0;
    for (chroms_of_the_set, positions_of_the_set) in
        chroms.chunks(at_a_time).zip(positions.chunks(at_a_time))
    {
        let num_vars = chroms_of_the_set.len();
        // The genotypes of the set are cut out of the block by the sizes
        // the block states, which `Block::check` has found to hold.
        let from = first.saturating_mul(alleles_per_var);
        let to = from.saturating_add(num_vars.saturating_mul(alleles_per_var));
        let Some(gts_of_the_variants) = block.gts.get(from..to) else {
            return Err(Error::BlockArrayOfAnotherSize {
                array: "gts",
                found: block.gts.len(),
                expected: to,
            });
        };
        gts_of_the_set.clear();
        gts_of_the_set.extend_from_slice(gts_of_the_variants);
        let set = the_dosages_of(
            &mut gts_of_the_set,
            num_vars,
            block.num_individuals,
            block.ploidy,
        )?;
        // Which of the variants of the window drop a candidate cannot
        // change with what the set does, so the whole window is read
        // against the whole set in one call.
        dropped.fill(false);
        if let (Some(chrom), Some(pos)) = (chroms_of_the_set.first(), positions_of_the_set.first())
        {
            // A variant of the window that the first variant of the set has
            // left behind is behind every variant of the set: the positions
            // rise and the chromosomes do not come back.
            window.leave_behind(*chrom, *pos, filter.max_dist);
        }
        if window.num_vars() > 0 {
            let values = window.num_vars().saturating_mul(num_vars);
            the_buffer_grown_to(
                &mut r2_against_the_window,
                values,
                "the r² of a set of candidates against the window",
            )?;
            let of_the_window = window.dosages()?;
            let Some(r2_of_the_window) = r2_against_the_window.get_mut(..values) else {
                // The buffer was grown to that many values, so this is not
                // reached.
                return Err(Error::LdR2OfAnotherSize {
                    num_values: r2_against_the_window.len(),
                    num_vars_of_a: window.num_vars(),
                    num_vars_of_b: num_vars,
                });
            };
            r2_between(&of_the_window, &set, r2_of_the_window)?;
            for ((r2_of_a_variant, chrom_of_it), pos_of_it) in r2_of_the_window
                .chunks_exact(num_vars)
                .zip(&window.chroms)
                .zip(&window.poss)
            {
                for ((dropped_it, r2), (chrom, pos)) in dropped
                    .iter_mut()
                    .zip(r2_of_a_variant)
                    .zip(chroms_of_the_set.iter().zip(positions_of_the_set))
                {
                    let within = chrom == chrom_of_it
                        && pos
                            .checked_sub(*pos_of_it)
                            .is_some_and(|dist| dist <= filter.max_dist);
                    if !within {
                        // The variants after this one are further ahead or
                        // on a later chromosome, and this variant of the
                        // window is behind the window of all of them.
                        break;
                    }
                    if *r2 > filter.max_allowed_r2 {
                        *dropped_it = true;
                    }
                }
            }
        }
        // The r² of every pair of the set, which is where a candidate is
        // compared with the variants kept inside it: at most 256 variants
        // are settled at a time, so this is 65536 values, 512 KB, whatever
        // the size of the block.
        let values = num_vars.saturating_mul(num_vars);
        let Some(r2_of_the_pairs) = r2_of_the_set.get_mut(..values) else {
            // The buffer holds the square of the variants a set settles at
            // a time and a set is at most that many, so this is not
            // reached.
            return Err(Error::LdR2OfAnotherSize {
                num_values: r2_of_the_set.len(),
                num_vars_of_a: num_vars,
                num_vars_of_b: num_vars,
            });
        };
        r2_between(&set, &set, r2_of_the_pairs)?;
        let mut after_it = dropped.as_mut_slice();
        for ((r2_of_the_variant, (chrom, pos)), variant) in r2_of_the_pairs
            .chunks_exact(num_vars)
            .zip(chroms_of_the_set.iter().zip(positions_of_the_set))
            .zip(0..num_vars)
        {
            let Some((dropped_it, after)) = std::mem::take(&mut after_it).split_first_mut() else {
                // One value was made for each variant of the set, so this
                // is not reached.
                break;
            };
            after_it = after;
            let stays = set.has_variance(variant) && !*dropped_it;
            keep.push(stays);
            if !stays {
                continue;
            }
            // It is kept, so it is in the window of every variant of the
            // set after it that is within `max_dist` on its chromosome.
            let next = variant.saturating_add(1);
            for ((dropped_later, r2), (chrom_later, pos_later)) in after_it
                .iter_mut()
                .zip(r2_of_the_variant.iter().skip(next))
                .zip(
                    chroms_of_the_set
                        .iter()
                        .skip(next)
                        .zip(positions_of_the_set.iter().skip(next)),
                )
            {
                let within = chrom_later == chrom
                    && pos_later
                        .checked_sub(*pos)
                        .is_some_and(|dist| dist <= filter.max_dist);
                if !within {
                    break;
                }
                if *r2 > filter.max_allowed_r2 {
                    *dropped_later = true;
                }
            }
            // The genotypes of the variant, which the window keeps and
            // builds its dosages from when the next set arrives.
            let of_it = variant.saturating_mul(alleles_per_var);
            let Some(gts_of_it) =
                gts_of_the_variants.get(of_it..of_it.saturating_add(alleles_per_var))
            else {
                return Err(Error::BlockArrayOfAnotherSize {
                    array: "gts",
                    found: block.gts.len(),
                    expected: to,
                });
            };
            window.push(*chrom, *pos, gts_of_it)?;
        }
        first = first.saturating_add(num_vars);
    }
    let kept = keep.iter().filter(|stays| **stays).count();
    Ok(TheBlockSettled { keep, kept, window })
}

/// A reader that gives the variants of its source that pass one filter.
///
/// It takes a block of its source at whatever size it comes, keeps in it
/// the variants that pass and gives it on, so it needs no
/// [`Reblock`](crate::block::Reblock) before it and the blocks it gives are
/// of uneven size. Nothing is kept
/// from one block to the next but the two counts of its filter, which
/// [`BlockReader::filtering_stats`] gives with the counts of the filters of
/// its source, its own first.
///
/// Several filters on one source are several of these, one over the other,
/// so each sees only what the one before it kept.
///
/// It keeps the contract of a reader of `docs/specs/block.md`: a block left
/// with no variant is not given and the next one is taken; after an error,
/// of its source or of its filter, it gives `None` at every call and does
/// not call its source again; and a source that gives a block of no
/// variants has a defect and is the error of that.
pub struct FilteredReader<R: BlockReader> {
    reader: R,
    filter: VarFilter,
    /// Whether the source has no more blocks or one of the two, the source
    /// or the filter, gave an error. After any of them there is no block.
    finished: bool,
}

impl<R: BlockReader> FilteredReader<R> {
    /// The reader that gives the variants of `reader` that pass `filter`.
    ///
    /// Building the chain asks `reader` for nothing: the consumer of the
    /// pass calls [`BlockReader::set_needs`] on the outermost reader of the
    /// chain, once it is built, and every filter of it adds the genotypes
    /// to what it passes on. A source that was narrowed to fields without
    /// the genotypes before it was wrapped, and that nobody asks again,
    /// gives blocks with no genotypes, and the filter fails at the first of
    /// them with the error of a field that is not in the block.
    ///
    /// # Errors
    ///
    /// When `reader` holds a filter of the kind of `filter` already, which
    /// its [`BlockReader::filtering_stats`] says: two threshold filters of
    /// one kind keep the variants that the stricter of the two keeps alone.
    /// The error carries the threshold of `filter`, and no threshold of the
    /// filter that is set: a chain says which kinds it holds and not with
    /// which thresholds.
    pub fn new(reader: R, filter: VarFilter) -> Result<FilteredReader<R>> {
        let criterion = filter.criterion();
        let kind = criterion.kind();
        if reader
            .filtering_stats()
            .iter()
            .any(|(of_the_chain, _)| *of_the_chain == kind)
        {
            return Err(Error::VarFilterOfAKindThatIsSet {
                kind,
                threshold: criterion.threshold(),
                threshold_that_is_set: None,
            });
        }
        Ok(FilteredReader {
            reader,
            filter,
            finished: false,
        })
    }
}

impl<R: BlockReader> BlockReader for FilteredReader<R> {
    /// The next block of the source with the variants that pass the filter
    /// kept in it, and the blocks that its filter emptied passed over.
    ///
    /// # Errors
    ///
    /// When the source fails; when a block of the source holds no variant,
    /// which no reader of popnei gives; and everything the filter refuses,
    /// which the `# Errors` of [`VarFilter::filter_block`] lists: a block
    /// whose arrays are not of its size, a block that has variants and no
    /// genotypes, a block of more memory than this machine addresses, and
    /// what the counts of one variant refuse. After any of them there is no
    /// block and the source is not called again.
    fn next_block(&mut self) -> Result<Option<Block>> {
        if self.finished {
            return Ok(None);
        }
        loop {
            let mut block = match self.reader.next_block() {
                Ok(Some(block)) => block,
                Ok(None) => {
                    self.finished = true;
                    return Ok(None);
                }
                Err(error) => {
                    self.finished = true;
                    return Err(error);
                }
            };
            // A source that gives a block of no variants has a defect, and
            // it is not asked again: over a source that always gives one, a
            // reader that asked again would never come back.
            if block.num_vars == 0 {
                self.finished = true;
                return Err(Error::ReaderGaveABlockOfNoVariants);
            }
            if let Err(error) = self.filter.filter_block(&mut block) {
                self.finished = true;
                return Err(error);
            }
            // A block the filter emptied is not given: the next one is
            // taken, and the source says when there are no more.
            if block.num_vars > 0 {
                return Ok(Some(block));
            }
        }
    }

    fn individuals(&self) -> &[String] {
        self.reader.individuals()
    }

    fn ploidy(&self) -> usize {
        self.reader.ploidy()
    }

    /// The table of the source: a reader over another reader has none of
    /// its own.
    fn chroms(&self) -> &ChromTable {
        self.reader.chroms()
    }

    /// The fields of the consumer and the genotypes, which the filter reads
    /// for every variant of every block: so the blocks this reader gives
    /// hold the genotypes also when the consumer did not ask for them.
    fn set_needs(&mut self, needs: Needs) {
        self.reader.set_needs(needs.union(Needs::GTS));
    }

    /// The counts of this filter, and after them those of the filters
    /// between the source and its own source.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        let mut stats = vec![(self.filter.criterion().kind(), self.filter.stats())];
        stats.extend(self.reader.filtering_stats());
        stats
    }
}

/// A reader that gives the variants of its source that the filter by
/// linkage disequilibrium keeps.
///
/// It is [`FilteredReader`] for an [`LdFilter`], and it keeps the same
/// contract of a reader of `docs/specs/block.md`: a block left with no
/// variant is not given and the next one is taken; after an error, of its
/// source or of its filter, it gives `None` at every call and does not call
/// its source again; and a source that gives a block of no variants has a
/// defect and is the error of that.
///
/// Two things are its own. Its filter holds the window of the variants kept
/// behind the variant it is reading, so what carries from one block to the
/// next is that window and not two counts alone. And it asks its source for
/// the chromosome and the position besides the genotypes, whatever its
/// consumer asked for, since the window of a variant is the variants kept
/// within `max_dist` base pairs of it on its chromosome: the blocks it
/// gives hold all three.
pub struct LdFilteredReader<R: BlockReader> {
    reader: R,
    filter: LdFilter,
    /// Whether the source has no more blocks or one of the two, the source
    /// or the filter, gave an error. After any of them there is no block.
    finished: bool,
}

impl<R: BlockReader> LdFilteredReader<R> {
    /// The reader that gives the variants of `reader` that `filter` keeps.
    ///
    /// Building the chain asks `reader` for nothing: the consumer of the
    /// pass calls [`BlockReader::set_needs`] on the outermost reader of the
    /// chain, once it is built, and this one adds the genotypes, the
    /// chromosome and the position to what it passes on. A source that was
    /// narrowed to fields without them before it was wrapped, and that
    /// nobody asks again, gives blocks the filter fails at with the error
    /// of a field that is not in the block.
    ///
    /// # Errors
    ///
    /// When `reader` holds a filter by linkage disequilibrium already,
    /// which its [`BlockReader::filtering_stats`] says: the second one
    /// would work its r² out over the variants the first left, which is not
    /// what either threshold asks for. The error carries the
    /// `max_allowed_r2` of `filter`, and no threshold of the filter that is
    /// set: a chain says which kinds it holds and not with which
    /// thresholds.
    pub fn new(reader: R, filter: LdFilter) -> Result<LdFilteredReader<R>> {
        let kind = THE_KIND_OF_THE_LD_FILTER;
        if reader
            .filtering_stats()
            .iter()
            .any(|(of_the_chain, _)| *of_the_chain == kind)
        {
            return Err(Error::VarFilterOfAKindThatIsSet {
                kind,
                threshold: filter.max_allowed_r2(),
                threshold_that_is_set: None,
            });
        }
        Ok(LdFilteredReader {
            reader,
            filter,
            finished: false,
        })
    }

    /// The error of a variant out of order with the name its chromosome has
    /// in the table of the source in place of its number, and every other
    /// error as it is.
    ///
    /// The filter is given blocks, which carry the number of a chromosome
    /// and not its name, and this reader has the table: the name is what a
    /// user looks the row up by in their file.
    fn the_chromosome_named(&self, error: Error) -> Error {
        let Error::LdFilterVariantOutOfOrder {
            variant,
            chrom,
            problem,
        } = error
        else {
            return error;
        };
        let named = match chrom {
            TheChromOfTheVariant::Numbered(number) => self
                .reader
                .chroms()
                .name(number)
                .map_or(TheChromOfTheVariant::Numbered(number), |name| {
                    TheChromOfTheVariant::Named(name.to_owned())
                }),
            // The filter gives the number, so a name is one this reader
            // has already put there.
            named @ TheChromOfTheVariant::Named(_) => named,
        };
        Error::LdFilterVariantOutOfOrder {
            variant,
            chrom: named,
            problem,
        }
    }
}

impl<R: BlockReader> BlockReader for LdFilteredReader<R> {
    /// The next block of the source with the variants that the filter keeps
    /// kept in it, and the blocks that its filter emptied passed over.
    ///
    /// # Errors
    ///
    /// When the source fails; when a block of the source holds no variant,
    /// which no reader of popnei gives; and everything the filter refuses,
    /// which the `# Errors` of [`LdFilter::filter_block`] lists: a block
    /// whose arrays are not of its size, a block that has variants and no
    /// genotypes or no position, a variant that does not come after the one
    /// before it, and what the dosages of a block refuse. After any of them
    /// there is no block and the source is not called again.
    fn next_block(&mut self) -> Result<Option<Block>> {
        if self.finished {
            return Ok(None);
        }
        loop {
            let mut block = match self.reader.next_block() {
                Ok(Some(block)) => block,
                Ok(None) => {
                    self.finished = true;
                    return Ok(None);
                }
                Err(error) => {
                    self.finished = true;
                    return Err(error);
                }
            };
            // A source that gives a block of no variants has a defect, and
            // it is not asked again: over a source that always gives one, a
            // reader that asked again would never come back.
            if block.num_vars == 0 {
                self.finished = true;
                return Err(Error::ReaderGaveABlockOfNoVariants);
            }
            if let Err(error) = self.filter.filter_block(&mut block) {
                self.finished = true;
                return Err(self.the_chromosome_named(error));
            }
            // A block the filter emptied is not given: the next one is
            // taken, and the source says when there are no more.
            if block.num_vars > 0 {
                return Ok(Some(block));
            }
        }
    }

    fn individuals(&self) -> &[String] {
        self.reader.individuals()
    }

    fn ploidy(&self) -> usize {
        self.reader.ploidy()
    }

    /// The table of the source: a reader over another reader has none of
    /// its own.
    fn chroms(&self) -> &ChromTable {
        self.reader.chroms()
    }

    /// The fields of the consumer, the genotypes, the chromosome and the
    /// position, which the filter reads for every variant of every block:
    /// so the blocks this reader gives hold the three also when the
    /// consumer asked for none of them.
    fn set_needs(&mut self, needs: Needs) {
        self.reader
            .set_needs(needs.union(Needs::GTS).union(Needs::CHROM_POS));
    }

    /// The counts of this filter, and after them those of the filters
    /// between the source and its own source.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        let mut stats = vec![(THE_KIND_OF_THE_LD_FILTER, self.filter.stats())];
        stats.extend(self.reader.filtering_stats());
        stats
    }
}

impl<R: BlockReader> fmt::Debug for LdFilteredReader<R> {
    /// What it filters and where it has got to. The source is left out, so
    /// that an `LdFilteredReader` over a reader that has no `Debug` has
    /// one.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("LdFilteredReader")
            .field("filter", &self.filter)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

/// One reader over `reader` for each criterion, in their order, so that
/// each filter sees what the one before it kept: the chain of the filters
/// of one pass. A [`FilteredReader`] for the three criteria that compare
/// one number of a variant and an [`LdFilteredReader`] for
/// [`MaxLdR2`](VarFilteringCriterion::MaxLdR2). No criterion gives `reader`
/// as it is.
///
/// Both binding crates build the chain of a pass with this, and neither
/// writes the loop: in which order the filters go, and what comes out while
/// they are built, are of the filters and not of Python or of TypeScript.
///
/// What it gives is the outermost reader of the chain, which whoever started
/// the pass holds: they read [`BlockReader::filtering_stats`] from it when
/// the consumer returns, and they lend it, `&mut`, to a consumer that takes
/// a reader. The fields the consumer wants are set on it once the chain is
/// built, which [`FilteredReader::new`] says why.
///
/// # Errors
///
/// What [`VarFilter::new`] and [`LdFilter::new`] refuse, a threshold that is
/// NaN, below 0 or above 1 and a `max_dist` below 1, and what
/// [`FilteredReader::new`] and [`LdFilteredReader::new`] refuse, a criterion
/// of the kind of one before it in `criteria` or of a filter that `reader`
/// holds already, which a chain built over a chain has. No block was read
/// then.
pub fn chain_of(
    reader: Box<dyn BlockReader>,
    criteria: &[VarFilteringCriterion],
) -> Result<Box<dyn BlockReader>> {
    let mut chain = reader;
    for criterion in criteria {
        chain = match *criterion {
            VarFilteringCriterion::MaxMissingRate(_)
            | VarFilteringCriterion::MaxMaf(_)
            | VarFilteringCriterion::MaxObsHet(_) => {
                Box::new(FilteredReader::new(chain, VarFilter::new(*criterion)?)?)
            }
            VarFilteringCriterion::MaxLdR2 {
                max_allowed_r2,
                max_dist,
            } => Box::new(LdFilteredReader::new(
                chain,
                LdFilter::new(max_allowed_r2, max_dist)?,
            )?),
        };
    }
    Ok(chain)
}

/// The error of a second filter of one kind, when `new` is of the kind of
/// one of `set`, the criteria of the filters that are set already.
///
/// Two threshold filters of one kind keep the variants that the stricter of
/// the two keeps alone, and a second filter by linkage disequilibrium works
/// its r² out over the variants the first left, so either says that the
/// user has lost track of what their variants carry. Both binding crates
/// call this when a user adds a filter to a `Variants`, where no reader
/// exists yet and the steps are what says which filters are set.
///
/// # Errors
///
/// When a criterion of `set` has the kind of `new`. The error carries both
/// thresholds, the one of `new` and the one that is set, where the same
/// error from [`FilteredReader::new`] carries the first alone: a chain of
/// readers says which kinds of filter it holds and not with which
/// thresholds.
pub fn refuse_a_second_filter_of_a_kind(
    set: &[VarFilteringCriterion],
    new: VarFilteringCriterion,
) -> Result<()> {
    let kind = new.kind();
    if let Some(that_is_set) = set.iter().find(|criterion| criterion.kind() == kind) {
        return Err(Error::VarFilterOfAKindThatIsSet {
            kind,
            threshold: new.threshold(),
            threshold_that_is_set: Some(that_is_set.threshold()),
        });
    }
    Ok(())
}

impl<R: BlockReader> fmt::Debug for FilteredReader<R> {
    /// What it filters and where it has got to. The source is left out, so
    /// that a `FilteredReader` over a reader that has no `Debug` has one.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("FilteredReader")
            .field("filter", &self.filter)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
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
/// Each thread counts the alleles into one array of its own, which it hands
/// to the counts for one row after another, so a block of a million
/// variants clears 128 numbers per variant and allocates nothing.
///
/// # Errors
///
/// What the counts of one variant refuse: a ploidy of 0, genotypes that are
/// not a whole number of genotypes of the ploidy, a variant of more alleles
/// than a count of them holds, and an allele below the missing one. The
/// error is the one of the first row that has one, wherever the threads
/// found it: which of two bad rows a thread reaches first depends on how
/// the rows were shared out, and a user who reports a damaged file has to
/// get the same message every time, so the rows are read again, one after
/// another, to find the first.
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

    let keep: Result<Vec<bool>> = gts
        .par_chunks_exact(alleles_per_var)
        .map_init(
            || [0_u32; 128],
            |counts, row| keeps(criterion, row, num_individuals, ploidy, counts),
        )
        .collect();
    match keep {
        Ok(keep) => Ok(keep),
        // The second pass costs a read of the block, and it is made only
        // where the block is refused and no variant of it is given.
        Err(of_a_thread) => match keep_of_the_rows_one_by_one(
            criterion,
            gts,
            alleles_per_var,
            num_individuals,
            ploidy,
        ) {
            Err(of_the_first_row) => Err(of_the_first_row),
            // The rows are the same rows, so the second pass finds an
            // error too; the error of the threads is what is left if it
            // ever did not.
            Ok(_) => Err(of_a_thread),
        },
    }
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

/// The rows read one after another, into one array of allele counts that
/// every row is counted into: what wasm does, what the test that compares
/// the two ways of reading them calls, and what the threads fall back on to
/// find the first row that is an error.
fn keep_of_the_rows_one_by_one(
    criterion: VarFilteringCriterion,
    gts: &[i8],
    alleles_per_var: usize,
    num_individuals: usize,
    ploidy: usize,
) -> Result<Vec<bool>> {
    let mut counts: AlleleCounts = [0; 128];
    gts.chunks_exact(alleles_per_var)
        .map(|row| keeps(criterion, row, num_individuals, ploidy, &mut counts))
        .collect()
}

/// Whether the variant whose genotypes are `gts` passes the criterion.
///
/// `gts` is one row of the genotypes of a block, the alleles of one
/// individual after those of the individual before it, `ploidy` alleles
/// each, and `num_individuals` is the individuals of the dataset, one at
/// least. A variant that has no number, one with no called allele or no
/// called genotype, is not kept, which the comparison of a NaN gives too.
///
/// `counts` is the array the counts of the alleles are read into, which
/// they clear themselves: the caller hands the same one over for every row
/// it reads, and the major allele frequency is the only criterion that
/// looks at it.
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
    counts: &mut AlleleCounts,
) -> Result<bool> {
    let number = match criterion {
        VarFilteringCriterion::MaxMissingRate(_) => {
            let gt_counts = count_gts(gts, ploidy)?;
            // The individuals of the dataset, and not the ones called at
            // this variant. A block of no individual holds no genotype and
            // never reaches this, so the divisor is 1 at least.
            f64::from(gt_counts.missing) / num_individuals as f64
        }
        VarFilteringCriterion::MaxMaf(_) => {
            let called_alleles = count_alleles(gts, counts)?;
            let Some(frequency) = the_major_allele_frequency(counts, called_alleles) else {
                // A variant with no called allele has no major allele
                // frequency, and this filter drops it.
                return Ok(false);
            };
            frequency
        }
        VarFilteringCriterion::MaxObsHet(_) => {
            let gt_counts = count_gts(gts, ploidy)?;
            if gt_counts.called == 0 {
                return Ok(false);
            }
            f64::from(gt_counts.het) / f64::from(gt_counts.called)
        }
        VarFilteringCriterion::MaxLdR2 { .. } => {
            // A criterion reaches this through a `VarFilter`, which is not
            // built for this one, so it is not reached: the r² of a variant
            // is not a number of the variant alone and is worked out by an
            // `LdFilter` against the window it holds.
            return Err(Error::VarFilterOfTheLdCriterion);
        }
    };
    Ok(number <= criterion.threshold())
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::{Path, PathBuf};
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    use super::{
        FilteredReader, FilteringStats, LdFilter, LdFilteredReader, THE_VARS_SETTLED_AT_A_TIME,
        TheChromOfTheVariant, TheOrderOfTheVariants, VarFilter, VarFilteringCriterion, chain_of,
        keep_of_the_rows, keep_of_the_rows_one_by_one, refuse_a_second_filter_of_a_kind,
    };
    use crate::block::{Block, BlockReader};
    use crate::error::{Error, Result};
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{ChromTable, MISSING_ALLELE, Needs};

    use VarFilteringCriterion::{MaxLdR2, MaxMaf, MaxMissingRate, MaxObsHet};

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
    /// which names `gts`. A block of no individual, or of the ploidy 0,
    /// holds no genotype and gives that error too: it is what such a block
    /// lacks, and no source of popnei has fewer than one individual of one
    /// allele.
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

        // Two variants of no individual, and two of five individuals of
        // the ploidy 0: `check` takes both, since their genotypes are 0
        // alleles, and the filter answers that the genotypes are not there.
        for (num_individuals, ploidy) in [(0, 2), (5, 0)] {
            let empty: [i8; 0] = [];
            let mut block = block_of(&[(1, &empty), (2, &empty)], num_individuals, ploidy);
            assert!(block.check().is_ok(), "{num_individuals} x {ploidy}");
            let mut filter = VarFilter::new(MaxMissingRate(1.0)).unwrap();
            let error = filter.filter_block(&mut block).unwrap_err();
            assert!(
                matches!(error, Error::FieldsNotInTheBlock { fields } if fields == Needs::GTS),
                "{num_individuals} x {ploidy}: {error}"
            );
            assert_eq!(block.num_vars, 2);
            assert_eq!(filter.stats(), FilteringStats::default());
        }
    }

    /// A block whose rows are two errors gives the error of the first of
    /// them, however the threads shared the rows out: a user who reports a
    /// damaged file has to get the same message every time they read it.
    ///
    /// The rows 3 and 390 of the block each hold an allele below the
    /// missing one, -2 and -9, and the error names -2 at every run. The
    /// rows are shared out anew at every read, so which of the two a thread
    /// reaches first changes from one read to the next.
    #[test]
    fn the_error_of_a_block_is_the_one_of_its_first_row_that_has_one() {
        let of_the_block = || {
            let mut rows: Vec<(u64, Vec<i8>)> = (1..=400)
                .map(|position| (position, vec![0_i8; 10]))
                .collect();
            rows[3].1[4] = -2;
            rows[390].1[7] = -9;
            let variants: Vec<(u64, &[i8])> = rows
                .iter()
                .map(|(position, gts)| (*position, gts.as_slice()))
                .collect();
            block_of(&variants, 5, 2)
        };
        for run in 1..=20 {
            // The major allele frequency reads the alleles one at a time
            // and the missing rate reads them as genotypes: both refuse an
            // allele below the missing one, and each has its own pass.
            for criterion in [MaxMaf(1.0), MaxMissingRate(1.0)] {
                let mut block = of_the_block();
                let error = VarFilter::new(criterion)
                    .unwrap()
                    .filter_block(&mut block)
                    .unwrap_err();
                assert!(
                    matches!(error, Error::AlleleBelowTheMissingOne { allele: -2 }),
                    "run {run}, {criterion:?}: {error}"
                );
            }
        }
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
        assert_eq!(
            MaxLdR2 {
                max_allowed_r2: 0.1,
                max_dist: 10_000,
            }
            .kind(),
            "ld"
        );
        assert_eq!(MaxMaf(0.0).kind(), MaxMaf(1.0).kind());
    }

    /// The criterion by linkage disequilibrium carries both the arguments a
    /// user wrote: its threshold is the `max_allowed_r2`, the first of the
    /// two values of the step a binding crate shows them, and its
    /// `max_dist` is the second. The three criteria that compare one number
    /// of a variant have no distance at all, so nothing of them reaches
    /// that second value.
    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "the threshold comes back as it was written and is not the end of any arithmetic"
    )]
    fn the_criterion_by_linkage_disequilibrium_carries_the_threshold_and_the_distance() {
        let criterion = MaxLdR2 {
            max_allowed_r2: 0.3,
            max_dist: 50_000,
        };

        assert_eq!(criterion.threshold(), 0.3);
        assert_eq!(criterion.max_dist(), Some(50_000));
        assert_eq!(MaxMissingRate(0.04).max_dist(), None);
        assert_eq!(MaxMaf(0.8).max_dist(), None);
        assert_eq!(MaxObsHet(0.5).max_dist(), None);
    }

    /// The criterion of the filter by linkage disequilibrium is no
    /// criterion of a `VarFilter`: the r² of a variant is worked out
    /// against the window of the variants kept behind it and not out of the
    /// counts of the variant alone. `chain_of` builds the reader of the
    /// right filter for it, so no call from Python or from TypeScript
    /// reaches this, and a caller of the core crate that builds the wrong
    /// one is told which filter that criterion is filtered by.
    #[test]
    fn a_filter_of_one_number_of_a_variant_is_refused_the_criterion_by_linkage_disequilibrium() {
        let error = VarFilter::new(MaxLdR2 {
            max_allowed_r2: 0.1,
            max_dist: 10_000,
        })
        .expect_err("the criterion was refused");

        let message = error.to_string();
        assert!(
            matches!(error, Error::VarFilterOfTheLdCriterion),
            "{message}"
        );
        assert!(message.contains("LdFilter"), "{message}");
        assert!(message.contains("chain_of"), "{message}");
    }

    /// The table of "How it is verified" of `docs/specs/filters.md`: each
    /// filter with its threshold, how many of the 500 variants of
    /// `many.vcf` it keeps, and the first five it keeps by position where
    /// the table gives them.
    ///
    /// The numbers are those of bcftools 1.24 and of pyNei at ef0ca6e,
    /// which keep the same variants at every one of these thresholds, and
    /// `tests/reference/filters/` holds every position each one keeps.
    const THE_TABLE: [(VarFilteringCriterion, usize, &[u64]); 9] = [
        (MaxMissingRate(0.0), 26, &[1259, 2110, 2480, 3072, 3257]),
        (MaxMissingRate(0.04), 215, &[]),
        (MaxMissingRate(0.1), 455, &[]),
        (MaxMaf(0.5), 35, &[1074, 1296, 1481, 1962, 2110]),
        (MaxMaf(0.8), 384, &[]),
        (MaxMaf(0.95), 480, &[]),
        (MaxObsHet(0.1), 22, &[1185, 3516, 3923, 4515, 5921]),
        (MaxObsHet(0.25), 79, &[]),
        (MaxObsHet(0.5), 369, &[]),
    ];

    /// The reference VCFs live at the root of the repository, beside the
    /// Python tests that read the same files, and not inside this crate.
    fn many_vcf() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/reference/vcf/many.vcf")
    }

    /// A reader over `many.vcf`, the 500 variants of 50 diploid individuals
    /// of `docs/specs/io_vcf.md`, with every variant given, the ones that
    /// failed their FILTER too, in blocks of `num_vars_per_block` variants
    /// and asked for `needs`.
    fn many_vcf_reader(
        num_vars_per_block: Option<usize>,
        needs: Needs,
    ) -> VcfReader<BufReader<File>> {
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block,
        };
        let mut reader =
            VcfReader::from_path(&many_vcf(), options).expect("the reader of many.vcf");
        reader.set_needs(needs);
        reader
    }

    /// Every block a reader gives, until it has no more or it fails.
    fn blocks_of(reader: &mut impl BlockReader) -> Result<Vec<Block>> {
        let mut blocks = Vec::new();
        while let Some(block) = reader.next_block()? {
            blocks.push(block);
        }
        Ok(blocks)
    }

    /// The positions of the variants of the blocks, in their order.
    fn positions_of_blocks(blocks: &[Block]) -> Vec<u64> {
        blocks
            .iter()
            .flat_map(|block| block.pos.clone().unwrap_or_default())
            .collect()
    }

    /// The name that `tests/reference/filters/make_reference.py` gives the
    /// file of one filter, its kind and its threshold: `maf_0.5`,
    /// `missing_data_0`.
    fn name_of(criterion: VarFilteringCriterion) -> String {
        format!(
            "{kind}_{threshold}",
            kind = criterion.kind(),
            threshold = criterion.threshold()
        )
    }

    /// The positions that bcftools 1.24 keeps of `many.vcf` with those
    /// filters, one per line in the file of that name under
    /// `tests/reference/filters/`, which task 3.1 of the plan stored.
    fn positions_of_the_reference(name: &str) -> Vec<u64> {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/filters")
            .join(format!("{name}.txt"));
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("{path}: {error}", path = path.display()));
        text.lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(|line| {
                line.parse()
                    .unwrap_or_else(|error| panic!("{name}: `{line}` is not a position: {error}"))
            })
            .collect()
    }

    /// The two counts of one filter.
    fn pair(vars_processed: u64, vars_kept: u64) -> FilteringStats {
        FilteringStats {
            vars_processed,
            vars_kept,
        }
    }

    /// The positions of the variants of `many.vcf` that the filter keeps,
    /// read in blocks of `num_vars_per_block`, with the counts of the
    /// filter after the last block.
    fn kept_of_many_vcf(
        criterion: VarFilteringCriterion,
        num_vars_per_block: Option<usize>,
    ) -> (Vec<u64>, Vec<(&'static str, FilteringStats)>) {
        let reader = many_vcf_reader(num_vars_per_block, Needs::GTS | Needs::CHROM_POS);
        let filter = VarFilter::new(criterion).expect("the filter");
        let mut filtered = FilteredReader::new(reader, filter).expect("the reader over the VCF");
        let blocks = blocks_of(&mut filtered).expect("the blocks");
        for block in &blocks {
            assert!(block.num_vars > 0, "a block with no variant was given");
            assert!(
                block.check().is_ok(),
                "a block whose arrays are not of its size"
            );
        }
        (positions_of_blocks(&blocks), filtered.filtering_stats())
    }

    /// The nine rows of the table of the spec, in blocks of 7 variants and
    /// in blocks of the size popnei chooses: the same variants are kept,
    /// and the counts of the filter are the 500 variants of the file and
    /// the ones it kept.
    #[test]
    fn each_filter_keeps_the_variants_of_many_vcf_that_the_table_of_the_spec_gives() {
        for (criterion, kept, first_five) in THE_TABLE {
            for num_vars_per_block in [Some(7), None] {
                let (positions, stats) = kept_of_many_vcf(criterion, num_vars_per_block);
                let what = format!("{criterion:?} in blocks of {num_vars_per_block:?}");
                assert_eq!(positions.len(), kept, "{what}");
                // Every position, and not the count and the first five
                // alone: the file holds the ones bcftools 1.24 keeps.
                assert_eq!(
                    positions,
                    positions_of_the_reference(&name_of(criterion)),
                    "{what}"
                );
                assert_eq!(
                    stats,
                    vec![(
                        criterion.kind(),
                        pair(500, u64::try_from(kept).expect("the variants kept"))
                    )],
                    "{what}"
                );
                if !first_five.is_empty() {
                    assert_eq!(&positions[..5], first_five, "{what}");
                }
            }
        }
    }

    /// The three filters of "How it is verified" of the counts over
    /// `many.vcf`, the missing data one at 0.04, the maf one at 0.8 after
    /// it and the observed heterozygosity one at 0.5 after that: bcftools
    /// 1.24 keeps 215, 163 and 106 variants, the first at 1111, 1407 and
    /// 1518, and pyNei gives the same three pairs of counts.
    #[test]
    fn the_three_filters_chained_over_many_vcf_keep_the_106_variants_and_give_their_counts() {
        for num_vars_per_block in [Some(7), None] {
            let reader = many_vcf_reader(num_vars_per_block, Needs::GTS | Needs::CHROM_POS);
            let missing_data =
                FilteredReader::new(reader, VarFilter::new(MaxMissingRate(0.04)).unwrap()).unwrap();
            let maf =
                FilteredReader::new(missing_data, VarFilter::new(MaxMaf(0.8)).unwrap()).unwrap();
            let mut obs_het =
                FilteredReader::new(maf, VarFilter::new(MaxObsHet(0.5)).unwrap()).unwrap();

            let blocks = blocks_of(&mut obs_het).expect("the blocks");
            let positions = positions_of_blocks(&blocks);
            let what = format!("in blocks of {num_vars_per_block:?}");
            assert_eq!(positions.len(), 106, "{what}");
            assert_eq!(positions[..3], [1111, 1407, 1518], "{what}");
            assert_eq!(
                positions,
                positions_of_the_reference("missing_data_0.04+maf_0.8+obs_het_0.5"),
                "{what}"
            );
            // The counts of the chain, the outermost filter first.
            assert_eq!(
                obs_het.filtering_stats(),
                vec![
                    ("obs_het", pair(163, 106)),
                    ("maf", pair(215, 163)),
                    ("missing_data", pair(500, 215)),
                ],
                "{what}"
            );
        }

        // The two filters of the chain without the third keep the 163
        // variants of the file of those two.
        let reader = many_vcf_reader(Some(7), Needs::GTS | Needs::CHROM_POS);
        let missing_data =
            FilteredReader::new(reader, VarFilter::new(MaxMissingRate(0.04)).unwrap()).unwrap();
        let mut maf =
            FilteredReader::new(missing_data, VarFilter::new(MaxMaf(0.8)).unwrap()).unwrap();
        let blocks = blocks_of(&mut maf).expect("the blocks");
        assert_eq!(
            positions_of_blocks(&blocks),
            positions_of_the_reference("missing_data_0.04+maf_0.8")
        );
        assert_eq!(
            maf.filtering_stats(),
            vec![("maf", pair(215, 163)), ("missing_data", pair(500, 215))]
        );
    }

    /// Two threshold filters of one kind keep the variants that the
    /// stricter of them keeps alone, so a second one is refused when the
    /// reader is built, with the kind and the threshold that was written.
    /// A filter of another kind over it is taken.
    #[test]
    fn a_second_filter_of_a_kind_the_chain_has_is_refused_when_the_reader_is_built() {
        let refused = |error: Error, threshold: &str| {
            let message = error.to_string();
            assert!(
                matches!(error, Error::VarFilterOfAKindThatIsSet { kind: "maf", .. }),
                "{message}"
            );
            assert!(message.contains("maf"), "{message}");
            assert!(message.contains(threshold), "{message}");
        };

        let reader = many_vcf_reader(Some(7), Needs::GTS);
        let maf = FilteredReader::new(reader, VarFilter::new(MaxMaf(0.8)).unwrap()).unwrap();
        let error = FilteredReader::new(maf, VarFilter::new(MaxMaf(0.95)).unwrap()).unwrap_err();
        refused(error, "0.95");

        // The same with a filter of another kind between the two.
        let reader = many_vcf_reader(Some(7), Needs::GTS);
        let maf = FilteredReader::new(reader, VarFilter::new(MaxMaf(0.8)).unwrap()).unwrap();
        let missing_data =
            FilteredReader::new(maf, VarFilter::new(MaxMissingRate(0.04)).unwrap()).unwrap();
        let obs_het =
            FilteredReader::new(missing_data, VarFilter::new(MaxObsHet(0.5)).unwrap()).unwrap();
        let error = FilteredReader::new(obs_het, VarFilter::new(MaxMaf(0.5)).unwrap()).unwrap_err();
        refused(error, "0.5");
    }

    /// A filter always needs the genotypes, so the blocks it gives hold
    /// them also when the consumer asked for the positions alone, and they
    /// hold no column that nobody asked for.
    #[test]
    fn the_blocks_hold_the_genotypes_when_the_consumer_asked_for_the_positions_alone() {
        let reader = many_vcf_reader(Some(7), Needs::ALL);
        let mut filtered =
            FilteredReader::new(reader, VarFilter::new(MaxMaf(0.5)).unwrap()).unwrap();
        filtered.set_needs(Needs::CHROM_POS);
        let blocks = blocks_of(&mut filtered).expect("the blocks");

        assert_eq!(positions_of_blocks(&blocks).len(), 35);
        for block in &blocks {
            assert!(block.fields().contains(Needs::GTS | Needs::CHROM_POS));
            assert!(!block.gts.is_empty());
            assert!(block.id.is_none());
            assert!(block.alleles.is_none());
            assert!(block.qual.is_none());
        }
    }

    /// Building the chain asks the source for nothing, so the consumer of
    /// a pass sets the fields it wants on the outermost reader once the
    /// chain is built: every filter passes them on with the genotypes
    /// added. A source that was narrowed to fields without the genotypes
    /// before it was wrapped, and that nobody asks again, gives blocks with
    /// no genotypes, and the filter fails at the first of them.
    #[test]
    fn the_consumer_sets_the_fields_on_the_outermost_reader_once_the_chain_is_built() {
        // The source was asked for the positions alone before it was
        // wrapped, and the chain is not asked for anything.
        let narrowed = many_vcf_reader(Some(7), Needs::CHROM_POS);
        let mut filtered =
            FilteredReader::new(narrowed, VarFilter::new(MaxMaf(0.5)).unwrap()).unwrap();
        let error = filtered.next_block().unwrap_err();
        assert!(
            matches!(error, Error::FieldsNotInTheBlock { fields } if fields == Needs::GTS),
            "{error}"
        );

        // The same source, with the fields set on the outermost reader
        // after the chain was built: the blocks hold the genotypes and the
        // variants are the 35 of the maf filter at 0.5.
        let narrowed = many_vcf_reader(Some(7), Needs::CHROM_POS);
        let mut filtered =
            FilteredReader::new(narrowed, VarFilter::new(MaxMaf(0.5)).unwrap()).unwrap();
        filtered.set_needs(Needs::CHROM_POS);
        let blocks = blocks_of(&mut filtered).expect("the blocks");
        assert_eq!(positions_of_blocks(&blocks).len(), 35);
    }

    /// The rows of a block are read on the threads of the pool the caller
    /// is in, so the variants that are kept and the counts are the same on
    /// one thread and on several.
    ///
    /// The pools are built here and are not rayon's global one, which has
    /// one thread per core of the machine. rayon is a dependency of the
    /// targets that are not wasm, so this test is compiled for those alone.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_variants_kept_are_the_same_on_one_thread_and_on_several() {
        let kept = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            pool.install(|| kept_of_many_vcf(MaxObsHet(0.25), Some(7)))
        };
        let (on_one, counts_of_one) = kept(1);
        let (on_four, counts_of_four) = kept(4);
        assert_eq!(on_one.len(), 79);
        assert_eq!(on_one, on_four);
        assert_eq!(counts_of_one, counts_of_four);
    }

    /// A reader of the tests that gives the blocks it was built with, so
    /// that the three rules of a reader of `docs/specs/block.md` are tested
    /// on a source that breaks them, which no reader of popnei does.
    struct GivenBlocks {
        individuals: Vec<String>,
        ploidy: usize,
        chroms: ChromTable,
        /// The blocks it has not given yet, the next one last.
        left: Vec<Block>,
        /// How many times it was asked for a block, which the test holds
        /// too, so that it sees whether the reader over it asked again
        /// after an error.
        calls: Arc<AtomicUsize>,
        /// The call at which it gives an error instead of a block.
        fails_at: Option<usize>,
        /// What it was last asked to fill.
        needs: Needs,
    }

    impl GivenBlocks {
        /// A reader of five individuals of the ploidy 2, which is what the
        /// worked example has, that gives `blocks` in their order.
        fn of(blocks: Vec<Block>) -> GivenBlocks {
            let mut chroms = ChromTable::new();
            chroms.intern("chr1");
            let mut left = blocks;
            left.reverse();
            GivenBlocks {
                individuals: (1..=5).map(|number| format!("ind{number}")).collect(),
                ploidy: 2,
                chroms,
                left,
                calls: Arc::new(AtomicUsize::new(0)),
                fails_at: None,
                needs: Needs::ALL,
            }
        }

        /// The same reader, whose call number `call` is an error.
        fn failing_at(blocks: Vec<Block>, call: usize) -> GivenBlocks {
            GivenBlocks {
                fails_at: Some(call),
                ..GivenBlocks::of(blocks)
            }
        }

        /// A reader of the six individuals of the worked example of the r²,
        /// which the blocks of the filter by linkage disequilibrium are
        /// built from, where the worked example of the three filters has
        /// five.
        fn of_the_r2_example(blocks: Vec<Block>) -> GivenBlocks {
            GivenBlocks {
                individuals: (1..=6).map(|number| format!("ind{number}")).collect(),
                ..GivenBlocks::of(blocks)
            }
        }

        /// How many times it has been asked for a block, which the test
        /// reads after the reader over it took it.
        fn calls(&self) -> Arc<AtomicUsize> {
            Arc::clone(&self.calls)
        }
    }

    impl BlockReader for GivenBlocks {
        fn next_block(&mut self) -> Result<Option<Block>> {
            let calls = self.calls.fetch_add(1, Ordering::SeqCst).saturating_add(1);
            if self.fails_at == Some(calls) {
                return Err(Error::Io(std::io::Error::other(
                    "the reader of the tests failed",
                )));
            }
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

    /// The block of the variants of the worked example whose numbers are
    /// given, in that order.
    fn block_of_the_worked_example(variants: &[usize]) -> Block {
        let rows: Vec<(u64, &[i8])> = variants
            .iter()
            .filter_map(|variant| THE_WORKED_EXAMPLE.get(*variant))
            .map(|(pos, gts)| (*pos, gts.as_slice()))
            .collect();
        block_of(&rows, 5, 2)
    }

    /// A block that the filter emptied is not given: the reader takes the
    /// next block of its source, and the variants of the block it dropped
    /// are in its counts.
    #[test]
    fn a_block_left_with_no_variant_is_not_given_and_the_next_one_is_taken() {
        // The first and the last block hold the variants 4 and 6, whose
        // missing rate is 1, and the middle one the variants 1 and 5,
        // whose missing rates are 0.2 and 0.
        let source = GivenBlocks::of(vec![
            block_of_the_worked_example(&[3, 5]),
            block_of_the_worked_example(&[0, 4]),
            block_of_the_worked_example(&[5, 3]),
        ]);
        let mut filtered =
            FilteredReader::new(source, VarFilter::new(MaxMissingRate(0.2)).unwrap()).unwrap();
        let blocks = blocks_of(&mut filtered).expect("the blocks");

        assert_eq!(blocks.len(), 1);
        assert_eq!(positions_of_blocks(&blocks), [1, 5]);
        assert_eq!(
            filtered.filtering_stats(),
            vec![("missing_data", pair(6, 2))]
        );
        // And there is no block after the last one.
        assert!(filtered.next_block().expect("no more blocks").is_none());
    }

    /// After an error of its source the reader gives no block and does not
    /// call its source again: a reader that went on would give the variants
    /// that follow a wrong one as if nothing had happened.
    #[test]
    fn after_an_error_of_the_source_there_is_no_block_and_the_source_is_not_called_again() {
        let source = GivenBlocks::failing_at(
            vec![
                block_of_the_worked_example(&[0, 4]),
                block_of_the_worked_example(&[1, 2]),
            ],
            2,
        );
        let calls = source.calls();
        let mut filtered =
            FilteredReader::new(source, VarFilter::new(MaxMissingRate(1.0)).unwrap()).unwrap();

        assert_eq!(
            positions_of_blocks(&[filtered.next_block().unwrap().unwrap()]),
            [1, 5]
        );
        let error = filtered.next_block().unwrap_err();
        assert!(matches!(error, Error::Io(_)), "{error}");
        assert!(filtered.next_block().unwrap().is_none());
        assert!(filtered.next_block().unwrap().is_none());
        assert_eq!(calls.load(Ordering::SeqCst), 2);
        // The counts are those of the blocks it did filter.
        assert_eq!(
            filtered.filtering_stats(),
            vec![("missing_data", pair(2, 2))]
        );
    }

    /// The same after an error of its own: a block that has variants and no
    /// genotypes is the error of a field that is not there, and the source
    /// is not called again.
    #[test]
    fn after_an_error_of_the_filter_there_is_no_block_and_the_source_is_not_called_again() {
        let mut without_the_genotypes = block_of_the_worked_example(&[0, 4]);
        without_the_genotypes.gts = Vec::new();
        let source = GivenBlocks::of(vec![
            without_the_genotypes,
            block_of_the_worked_example(&[1, 2]),
        ]);
        let calls = source.calls();
        let mut filtered =
            FilteredReader::new(source, VarFilter::new(MaxMaf(1.0)).unwrap()).unwrap();

        let error = filtered.next_block().unwrap_err();
        assert!(
            matches!(error, Error::FieldsNotInTheBlock { fields } if fields == Needs::GTS),
            "{error}"
        );
        assert!(filtered.next_block().unwrap().is_none());
        assert_eq!(calls.load(Ordering::SeqCst), 1);
        assert_eq!(filtered.filtering_stats(), vec![("maf", pair(0, 0))]);
    }

    /// A source that gives a block of no variants has a defect, and it is
    /// the error `reblock` gives for it: a reader that asked again would
    /// never come back from a source that always gives one.
    #[test]
    fn a_source_that_gives_a_block_of_no_variants_is_an_error_and_is_not_called_again() {
        let source = GivenBlocks::of(vec![
            block_of_the_worked_example(&[]),
            block_of_the_worked_example(&[0, 4]),
        ]);
        let calls = source.calls();
        let mut filtered =
            FilteredReader::new(source, VarFilter::new(MaxMissingRate(1.0)).unwrap()).unwrap();

        let error = filtered.next_block().unwrap_err();
        assert!(
            matches!(error, Error::ReaderGaveABlockOfNoVariants),
            "{error}"
        );
        assert!(filtered.next_block().unwrap().is_none());
        assert_eq!(calls.load(Ordering::SeqCst), 1);
    }

    /// A reader over a reader has no individuals, no ploidy and no table of
    /// chromosome names of its own: it gives those of its source, whose
    /// ploidy here is 4 and not the 2 of the blocks of the other tests.
    #[test]
    fn the_individuals_the_ploidy_and_the_chromosomes_are_those_of_the_source() {
        let source = GivenBlocks {
            ploidy: 4,
            ..GivenBlocks::of(Vec::new())
        };
        let filtered =
            FilteredReader::new(source, VarFilter::new(MaxObsHet(1.0)).unwrap()).unwrap();

        assert_eq!(filtered.individuals().len(), 5);
        assert_eq!(
            filtered.individuals().first().map(String::as_str),
            Some("ind1")
        );
        assert_eq!(filtered.ploidy(), 4);
        assert_eq!(filtered.chroms().name(0), Some("chr1"));
        assert_eq!(filtered.filtering_stats(), vec![("obs_het", pair(0, 0))]);
    }

    /// The chain that `chain_of` builds from the three criteria of "How it
    /// is verified" of the counts, the missing data one at 0.04, the maf
    /// one at 0.8 and the observed heterozygosity one at 0.5, over
    /// `many.vcf`: the same 106 variants and the same three pairs of counts
    /// that the filters put one over another by hand give, since it is what
    /// the function does.
    #[test]
    fn chain_of_the_three_criteria_over_many_vcf_keeps_the_106_variants_with_their_counts() {
        let source = Box::new(many_vcf_reader(Some(7), Needs::GTS | Needs::CHROM_POS));
        let mut chain = chain_of(source, &[MaxMissingRate(0.04), MaxMaf(0.8), MaxObsHet(0.5)])
            .expect("the chain of the three criteria");

        let blocks = blocks_of(&mut chain).expect("the blocks");
        let positions = positions_of_blocks(&blocks);
        assert_eq!(positions.len(), 106);
        assert_eq!(positions[..3], [1111, 1407, 1518]);
        assert_eq!(
            positions,
            positions_of_the_reference("missing_data_0.04+maf_0.8+obs_het_0.5")
        );
        // The counts of the chain, the outermost filter first, which is the
        // last criterion.
        assert_eq!(
            chain.filtering_stats(),
            vec![
                ("obs_het", pair(163, 106)),
                ("maf", pair(215, 163)),
                ("missing_data", pair(500, 215)),
            ]
        );
    }

    /// No criterion gives the source as it is: every variant it has, the
    /// ones no filter would keep among them, and no counts.
    #[test]
    fn chain_of_no_criterion_gives_the_source_as_it_is() {
        let source = GivenBlocks::of(vec![block_of_the_worked_example(&[0, 1, 2, 3, 4, 5])]);
        let mut chain = chain_of(Box::new(source), &[]).expect("the chain of no criterion");

        let blocks = blocks_of(&mut chain).expect("the blocks");

        // The variants 4 and 6 have no called genotype, so neither has an
        // observed heterozygosity and a filter of it would keep neither at
        // any threshold: they come out here because no filter is over the
        // source.
        assert_eq!(positions_of_blocks(&blocks), [1, 2, 3, 4, 5, 6]);
        assert!(chain.filtering_stats().is_empty());
    }

    /// A criterion of the kind of one before it is the error of
    /// `FilteredReader::new`, which the chain runs into while it is built,
    /// and a criterion of another kind between the two changes nothing.
    #[test]
    fn chain_of_a_criterion_of_a_kind_that_is_set_is_the_error_of_the_reader() {
        let refused = |criteria: &[VarFilteringCriterion]| {
            let source = Box::new(many_vcf_reader(Some(7), Needs::GTS));
            let error = chain_of(source, criteria)
                .err()
                .expect("the chain was refused");
            let message = error.to_string();
            assert!(
                matches!(error, Error::VarFilterOfAKindThatIsSet { kind: "maf", .. }),
                "{message}"
            );
            message
        };

        assert!(refused(&[MaxMaf(0.8), MaxMaf(0.95)]).contains("0.95"));
        assert!(
            refused(&[MaxMaf(0.8), MaxMissingRate(0.04), MaxMaf(0.5)]).contains("0.5"),
            "a criterion of another kind between the two"
        );
    }

    /// A criterion of the kind of a filter that the reader handed in holds
    /// is refused too: a binding crate reaches it by building a chain over
    /// a chain, and it is `FilteredReader::new` that says which kinds a
    /// reader holds.
    #[test]
    fn chain_of_a_criterion_of_a_kind_the_reader_holds_is_the_error_of_the_reader() {
        let source = Box::new(many_vcf_reader(Some(7), Needs::GTS));
        let with_a_maf_filter = chain_of(source, &[MaxMaf(0.8)]).expect("the chain of one");

        let error = chain_of(with_a_maf_filter, &[MaxMissingRate(0.04), MaxMaf(0.95)])
            .err()
            .expect("the chain over it was refused");

        let message = error.to_string();
        assert!(
            matches!(error, Error::VarFilterOfAKindThatIsSet { kind: "maf", .. }),
            "{message}"
        );
        assert!(message.contains("0.95"), "{message}");
    }

    /// The criteria that are set and a new one, which is what a binding
    /// crate has when a user adds a filter: no reader exists then, and the
    /// steps of the `Variants` are the criteria that are set.
    #[test]
    fn refuse_a_second_filter_of_a_kind_takes_a_kind_that_is_not_set_and_refuses_one_that_is() {
        let set = [MaxMissingRate(0.04), MaxMaf(0.8)];

        assert!(refuse_a_second_filter_of_a_kind(&set, MaxObsHet(0.5)).is_ok());
        assert!(refuse_a_second_filter_of_a_kind(&[], MaxMaf(0.95)).is_ok());
        // A criterion of another kind between the two changes nothing: the
        // kind is looked for among all of them.
        assert!(
            refuse_a_second_filter_of_a_kind(&[MaxMaf(0.8), MaxObsHet(0.5)], MaxMaf(0.95)).is_err()
        );
    }

    /// The message names the kind and both thresholds, the one that is set
    /// and the one that was written, which is what a user needs in order to
    /// see which of their cells they ran twice.
    #[test]
    fn refuse_a_second_filter_of_a_kind_names_the_kind_and_both_thresholds() {
        let error = refuse_a_second_filter_of_a_kind(&[MaxMaf(0.8)], MaxMaf(0.95))
            .expect_err("the second maf filter was refused");

        let message = error.to_string();
        assert!(
            matches!(
                error,
                Error::VarFilterOfAKindThatIsSet {
                    kind: "maf",
                    threshold_that_is_set: Some(_),
                    ..
                }
            ),
            "{message}"
        );
        assert!(message.contains("maf"), "{message}");
        assert!(message.contains("0.95"), "{message}");
        assert!(message.contains("0.8"), "{message}");
    }

    /// A threshold that is not a number from 0 to 1 is the error of
    /// `VarFilter::new`, whatever its place among the criteria.
    #[test]
    fn chain_of_a_threshold_out_of_range_is_the_error_of_the_filter() {
        let source = Box::new(many_vcf_reader(Some(7), Needs::GTS));
        let error = chain_of(source, &[MaxMissingRate(0.04), MaxMaf(1.5)])
            .err()
            .expect("the chain was refused");

        assert!(
            matches!(
                error,
                Error::VarFilterThresholdOutOfRange { kind: "maf", .. }
            ),
            "{error}"
        );
    }

    /// A chain of two filters over a source reports the counts of both, the
    /// outermost first, and the source reports none of its own.
    #[test]
    fn a_chain_of_two_filters_reports_the_counts_of_both_the_outermost_first() {
        let source = GivenBlocks::of(vec![block_of_the_worked_example(&[0, 1, 2, 3, 4, 5])]);
        assert!(source.filtering_stats().is_empty());
        let missing_data =
            FilteredReader::new(source, VarFilter::new(MaxMissingRate(0.4)).unwrap()).unwrap();
        let mut maf =
            FilteredReader::new(missing_data, VarFilter::new(MaxMaf(0.88)).unwrap()).unwrap();

        let blocks = blocks_of(&mut maf).expect("the blocks");
        assert_eq!(positions_of_blocks(&blocks), [2, 3, 5]);
        assert_eq!(
            maf.filtering_stats(),
            vec![("maf", pair(4, 3)), ("missing_data", pair(6, 4))]
        );
    }

    /// An allele that was not called, the `.` of a VCF.
    const M: i8 = MISSING_ALLELE;

    /// The five variants of six diploid individuals of the worked example
    /// of "How it is verified" of `docs/specs/ld.md`, each with the
    /// position that example gives it, and each row the alleles of one
    /// individual after those of the individual before it.
    ///
    /// Their r², which plink2 v2.0.0-a.7.7 gives for the same file,
    /// `tests/reference/ld/example.vcf`, are 0.675 for v1 and v2,
    /// 0.7544642857142857 for v1 and v3, 0.0625 for v1 and v5,
    /// 0.6428571428571429 for v2 and v3, 0 for v2 and v5, and 0.21875 for
    /// v3 and v5. v4 has one dosage in every individual, so it has no r²
    /// against any of them.
    const THE_WORKED_EXAMPLE_OF_THE_R2: [(u64, [i8; 12]); 5] = [
        // v1 0/0 0/0 0/1 0/1 1/1 1/1
        (1000, [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1]),
        // v2 0/0 0/1 0/1 1/1 1/1 1/1
        (2000, [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1]),
        // v3 0/0 0/0 0/0 0/1 ./. 1/1
        (3000, [0, 0, 0, 0, 0, 0, 0, 1, M, M, 1, 1]),
        // v4 0/0 0/0 0/0 0/0 0/0 0/0
        (4000, [0; 12]),
        // v5 0/1 1/1 0/0 0/1 1/1 0/0
        (5000, [0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0]),
    ];

    /// A block of the variants given, each with its chromosome and its
    /// position, of `num_individuals` individuals of the ploidy `ploidy`:
    /// the fields the filter by linkage disequilibrium asks its source for.
    fn block_of_the_chromosomes(
        variants: &[(u32, u64, &[i8])],
        num_individuals: usize,
        ploidy: usize,
    ) -> Block {
        let mut gts = Vec::new();
        let mut chrom = Vec::new();
        let mut pos = Vec::new();
        for (chromosome, position, row) in variants {
            gts.extend_from_slice(row);
            chrom.push(*chromosome);
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

    /// A block of the variants of the worked example of the r² that are
    /// named, in the order they are named, each at its own position and all
    /// of them on one chromosome.
    fn block_of_the_r2_example(variants: &[usize]) -> Block {
        let rows: Vec<(u32, u64, &[i8])> = variants
            .iter()
            .filter_map(|variant| THE_WORKED_EXAMPLE_OF_THE_R2.get(*variant))
            .map(|(pos, gts)| (0, *pos, gts.as_slice()))
            .collect();
        block_of_the_chromosomes(&rows, 6, 2)
    }

    /// The positions of the variants the filter keeps of the blocks, which
    /// it takes one after another as a reader over a reader gives them.
    fn kept_of_the_blocks(filter: &mut LdFilter, blocks: Vec<Block>) -> Result<Vec<u64>> {
        let mut kept = Vec::new();
        for mut block in blocks {
            filter.filter_block(&mut block)?;
            kept.extend(positions_of(&block));
        }
        Ok(kept)
    }

    /// The variants of the worked example of the r² that the filter keeps
    /// at that setting, by their positions, with its counts.
    fn kept_of_the_r2_example(max_allowed_r2: f64, max_dist: u64) -> (Vec<u64>, FilteringStats) {
        let mut filter = LdFilter::new(max_allowed_r2, max_dist).unwrap();
        let block = block_of_the_r2_example(&[0, 1, 2, 3, 4]);
        let kept = kept_of_the_blocks(&mut filter, vec![block]).expect("the block was filtered");
        (kept, filter.stats())
    }

    /// The three rows of the table of the worked example of "How it is
    /// verified" of the item of `docs/specs/filters.md`: at 5000 bp and
    /// 0.5 the variants kept are v1 and v5, since v2 and v3 are above 0.5
    /// against v1 and v4 has one dosage; at 5000 and 0.7 they are v1, v2
    /// and v5, since v2 is 0.675 against v1 and v3 is 0.754; and at 1000
    /// and 0.5 they are v1, v3 and v5, since v3 and v5 have no kept variant
    /// within 1000 bp. The counts of the first row are 5 variants given and
    /// 2 kept.
    #[test]
    fn the_filter_keeps_the_variants_of_the_worked_example_at_each_setting() {
        assert_eq!(kept_of_the_r2_example(0.5, 5000).0, [1000, 5000]);
        assert_eq!(kept_of_the_r2_example(0.7, 5000).0, [1000, 2000, 5000]);
        assert_eq!(kept_of_the_r2_example(0.5, 1000).0, [1000, 3000, 5000]);
        assert_eq!(
            kept_of_the_r2_example(0.5, 5000).1,
            FilteringStats {
                vars_processed: 5,
                vars_kept: 2,
            }
        );
    }

    /// A candidate is compared with every variant of its window and not
    /// with the last kept one alone, which is what pyNei's
    /// `_filter_chunk_by_ld` does: at 5000 bp and 0.7 the last variant kept
    /// before v3 is v2, whose r² against it is 0.643 and below the
    /// threshold, and v3 goes because of v1, which is 0.754 against it and
    /// is in its window too.
    ///
    /// v1, v2 and v3 are given in blocks of their own, so that v1 and v2
    /// are in the window the filter carries over and not in the set of
    /// variants v3 is settled with: a filter that compared v3 with the last
    /// kept variant alone would keep it here, where one that compares it
    /// with its whole window drops it. The same three variants in one block
    /// say nothing about that, since there v3 is read out of the r² of its
    /// own set.
    #[test]
    fn a_variant_is_compared_with_every_variant_of_its_window_and_not_the_last_kept_one_alone() {
        assert_eq!(kept_of_the_r2_example(0.7, 5000).0, [1000, 2000, 5000]);

        let mut filter = LdFilter::new(0.7, 5000).unwrap();
        let of_its_own: Vec<Block> = [0, 1, 2, 4]
            .into_iter()
            .map(|variant| block_of_the_r2_example(&[variant]))
            .collect();
        let kept = kept_of_the_blocks(&mut filter, of_its_own).expect("the blocks");
        assert_eq!(kept, [1000, 2000, 5000]);
    }

    /// A variant whose r² is exactly the threshold stays, as a variant
    /// whose number is exactly the threshold stays in the three filters
    /// that compare one number: the r² of v1 and v2 is 27/40, which is
    /// 0.675, and v2 stays at that threshold and goes at 0.674.
    ///
    /// The two variants are given in one block and in two, the second of
    /// which compares v2 with the window the filter carried over: the r²
    /// of a pair of one set and the r² of a candidate against the window
    /// are two paths, and each of them keeps a pair exactly at the
    /// threshold.
    #[test]
    fn a_variant_whose_r2_is_exactly_the_threshold_stays() {
        assert_eq!(kept_of_the_r2_example(0.675, 5000).0, [1000, 2000, 5000]);
        assert_eq!(kept_of_the_r2_example(0.674, 5000).0, [1000, 5000]);

        for (threshold, of_two_blocks) in [(0.675, vec![1000, 2000]), (0.674, vec![1000])] {
            let mut filter = LdFilter::new(threshold, 5000).unwrap();
            let blocks = vec![block_of_the_r2_example(&[0]), block_of_the_r2_example(&[1])];
            let kept = kept_of_the_blocks(&mut filter, blocks).expect("the blocks");
            assert_eq!(kept, of_two_blocks, "at a threshold of {threshold}");
        }
    }

    /// At a threshold of 1 every variant with two dosages stays, whatever
    /// its r² against the variants of its window, since no r² is above 1;
    /// v4, whose called genotypes hold one dosage, goes at every threshold.
    #[test]
    fn at_a_threshold_of_1_every_variant_with_two_dosages_stays() {
        assert_eq!(
            kept_of_the_r2_example(1.0, 5000).0,
            [1000, 2000, 3000, 5000]
        );
    }

    /// A variant whose called genotypes hold one dosage, and one with no
    /// called genotype, are dropped wherever they are, the first place
    /// included: pyNei keeps the first variant of the first chunk whatever
    /// it is, and when every one of its genotypes holds the same value
    /// nothing else is ever kept, since every r against it is NaN.
    #[test]
    fn a_variant_of_one_dosage_or_of_none_is_dropped_wherever_it_is() {
        let of_one_dosage = THE_WORKED_EXAMPLE_OF_THE_R2[3].1;
        let called_nowhere = [M; 12];
        for first in [of_one_dosage, called_nowhere] {
            let mut variants: Vec<(u32, u64, &[i8])> = vec![(0, 500, first.as_slice())];
            variants.extend(
                [0, 1, 4]
                    .iter()
                    .filter_map(|variant| THE_WORKED_EXAMPLE_OF_THE_R2.get(*variant))
                    .map(|(pos, gts)| (0, *pos, gts.as_slice())),
            );
            let block = block_of_the_chromosomes(&variants, 6, 2);
            let mut filter = LdFilter::new(0.5, 5000).unwrap();
            let kept = kept_of_the_blocks(&mut filter, vec![block]).expect("the block");
            assert_eq!(kept, [1000, 5000]);
            assert_eq!(filter.stats(), pair(4, 2));
        }
    }

    /// A pair whose r² is not defined does not drop the candidate: the two
    /// variants here have two dosages each and no individual called at
    /// both, so their r² is NaN, and both stay at a threshold of 0, where
    /// two variants that say the same thing would leave one.
    ///
    /// Each pair is given in one block and in two, since a candidate
    /// against the window and a pair of one set are two paths and a NaN
    /// drops the candidate on neither.
    #[test]
    fn a_pair_with_no_r2_does_not_drop_the_candidate() {
        // 0/0 0/1 1/1 ./. ./. ./. and ./. ./. ./. 0/0 0/1 1/1.
        let called_first = [0, 0, 0, 1, 1, 1, M, M, M, M, M, M];
        let called_last = [M, M, M, M, M, M, 0, 0, 0, 1, 1, 1];
        let pair_with_no_r2 = block_of_the_chromosomes(
            &[
                (0, 1000, called_first.as_slice()),
                (0, 2000, called_last.as_slice()),
            ],
            6,
            2,
        );
        let mut filter = LdFilter::new(0.0, 5000).unwrap();
        let kept = kept_of_the_blocks(&mut filter, vec![pair_with_no_r2]).expect("the block");
        assert_eq!(kept, [1000, 2000]);

        // The same threshold over two variants that are the same leaves
        // the second one, so what keeps the pair above is the NaN.
        let the_same_twice = block_of_the_chromosomes(
            &[
                (0, 1000, called_first.as_slice()),
                (0, 2000, called_first.as_slice()),
            ],
            6,
            2,
        );
        let mut filter = LdFilter::new(0.0, 5000).unwrap();
        let kept = kept_of_the_blocks(&mut filter, vec![the_same_twice]).expect("the block");
        assert_eq!(kept, [1000]);

        // The two of each pair in blocks of their own: the second variant
        // is then compared with the first through the window and not
        // through the r² of one set.
        let in_two_blocks = |second: &[i8]| {
            let blocks = vec![
                block_of_the_chromosomes(&[(0, 1000, called_first.as_slice())], 6, 2),
                block_of_the_chromosomes(&[(0, 2000, second)], 6, 2),
            ];
            let mut filter = LdFilter::new(0.0, 5000).expect("the filter");
            kept_of_the_blocks(&mut filter, blocks).expect("the blocks")
        };
        assert_eq!(in_two_blocks(called_last.as_slice()), [1000, 2000]);
        assert_eq!(in_two_blocks(called_first.as_slice()), [1000]);
    }

    /// The variants kept do not change with the size of the blocks: the
    /// rule reads the positions of the variants and never a block
    /// boundary, and the window carries over from one block to the next.
    #[test]
    fn the_variants_kept_do_not_change_with_the_size_of_the_blocks() {
        for size in 1..=5 {
            let blocks: Vec<Block> = (0..5)
                .collect::<Vec<usize>>()
                .chunks(size)
                .map(block_of_the_r2_example)
                .collect();
            let mut filter = LdFilter::new(0.5, 5000).unwrap();
            let kept = kept_of_the_blocks(&mut filter, blocks).expect("the blocks");
            assert_eq!(kept, [1000, 5000], "in blocks of {size} variants");
            assert_eq!(filter.stats(), pair(5, 2), "in blocks of {size} variants");
        }
    }

    /// The variants of one chromosome, each of them the one before it with
    /// about one allele in ten changed, so that the variants near each
    /// other are linked and the ones far apart are not.
    ///
    /// The numbers are of a generator written here and of no reference
    /// program: what the test that reads them asserts is that two ways of
    /// cutting one dataset into blocks keep the same variants, which needs
    /// a dataset the filter drops some variants of and keeps others.
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "the counts are of a dataset of 600 variants of 20 individuals built here"
    )]
    fn a_linked_chromosome(num_vars: usize, num_individuals: usize) -> Vec<(u64, Vec<i8>)> {
        let mut alleles: Vec<i8> = (0..num_individuals * 2)
            .map(|allele| i8::from(allele % 3 == 0))
            .collect();
        let mut seed: u64 = 20_260_923;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            seed >> 33
        };
        (1..=num_vars)
            .map(|variant| {
                for allele in &mut alleles {
                    if next() % 10 == 0 {
                        *allele = 1 - *allele;
                    }
                }
                (variant as u64 * 1000, alleles.clone())
            })
            .collect()
    }

    /// The variants kept do not change with how many of them the filter
    /// settles at a time, which is what
    /// `ld::tests::neither_the_blocks_nor_the_tiles_change_the_matrix` does
    /// for the tiles of the matrix: one block of 600 variants keeps the
    /// same set at 1, at 3 and at the 256 of `THE_VARS_SETTLED_AT_A_TIME`.
    ///
    /// At 1 every candidate is compared with its whole window and with
    /// nothing else, and at 256 the last 344 of them are compared with the
    /// variants kept inside their own set as well, so the two paths give
    /// the same set.
    #[test]
    fn the_variants_kept_do_not_change_with_the_variants_settled_at_a_time() {
        let num_vars = 600;
        let variants = a_linked_chromosome(num_vars, 20);
        let kept_at = |at_a_time| {
            let rows: Vec<(u32, u64, &[i8])> = variants
                .iter()
                .map(|(pos, gts)| (0, *pos, gts.as_slice()))
                .collect();
            let mut block = block_of_the_chromosomes(&rows, 20, 2);
            let mut filter = LdFilter::new(0.3, 10_000).expect("the filter");
            filter
                .the_block_filtered(&mut block, at_a_time)
                .expect("the block");
            positions_of(&block)
        };
        let at_the_default = kept_at(THE_VARS_SETTLED_AT_A_TIME);
        // A set that is neither every variant nor one of them: a dataset
        // the filter does nothing on would tell no two ways of settling it
        // apart.
        assert!(
            at_the_default.len() > 20 && at_the_default.len() < num_vars - 20,
            "{} variants of {num_vars} kept",
            at_the_default.len()
        );
        assert_eq!(kept_at(1), at_the_default, "settled one at a time");
        assert_eq!(kept_at(3), at_the_default, "settled three at a time");
    }

    /// The variants kept do not change with the size of the blocks over a
    /// dataset of more variants than the filter settles at a time, so that
    /// the sets it works the r² of in one product are cut in one place when
    /// the blocks are of 7 variants and in another when they are of 600.
    #[test]
    fn the_variants_kept_do_not_change_with_the_size_of_the_blocks_over_600_variants() {
        let num_vars = 600;
        assert!(num_vars > THE_VARS_SETTLED_AT_A_TIME);
        let variants = a_linked_chromosome(num_vars, 20);
        let of_one_block = {
            let rows: Vec<(u32, u64, &[i8])> = variants
                .iter()
                .map(|(pos, gts)| (0, *pos, gts.as_slice()))
                .collect();
            let mut filter = LdFilter::new(0.3, 10000).unwrap();
            kept_of_the_blocks(&mut filter, vec![block_of_the_chromosomes(&rows, 20, 2)])
                .expect("the block")
        };
        // A set that is neither every variant nor one of them: a dataset
        // the filter does nothing on would tell no two ways of reading it
        // apart.
        assert!(
            of_one_block.len() > 20 && of_one_block.len() < num_vars - 20,
            "{} variants of {num_vars} kept",
            of_one_block.len()
        );

        for size in [7, 64, 256, 257] {
            let blocks: Vec<Block> = variants
                .chunks(size)
                .map(|of_a_block| {
                    let rows: Vec<(u32, u64, &[i8])> = of_a_block
                        .iter()
                        .map(|(pos, gts)| (0, *pos, gts.as_slice()))
                        .collect();
                    block_of_the_chromosomes(&rows, 20, 2)
                })
                .collect();
            let mut filter = LdFilter::new(0.3, 10000).unwrap();
            let kept = kept_of_the_blocks(&mut filter, blocks).expect("the blocks");
            assert_eq!(kept, of_one_block, "in blocks of {size} variants");
        }
    }

    /// The window ends with the chromosome: the first variant of a
    /// chromosome with two dosages is kept although it says exactly what
    /// the last variant of the chromosome before it said, where pyNei
    /// carries the last kept variant from one chromosome to the next.
    #[test]
    fn the_first_variant_of_a_chromosome_is_kept_although_it_repeats_the_last_of_the_one_before() {
        let (pos, gts) = THE_WORKED_EXAMPLE_OF_THE_R2[0];
        let two_chromosomes =
            block_of_the_chromosomes(&[(0, pos, gts.as_slice()), (1, pos, gts.as_slice())], 6, 2);
        let mut filter = LdFilter::new(0.0, 250_000).unwrap();
        let kept = kept_of_the_blocks(&mut filter, vec![two_chromosomes]).expect("the block");
        assert_eq!(kept, [1000, 1000]);
        assert_eq!(filter.stats(), pair(2, 2));
    }

    /// Two variants at one position are 0 apart, so each is in the window
    /// of the other whatever `max_dist` is: at the smallest window of 1
    /// base pair, v2 goes against v1 at a threshold of 0.5, their r² being
    /// 0.675.
    #[test]
    fn two_variants_at_one_position_are_each_in_the_window_of_the_other() {
        let (_, v1) = THE_WORKED_EXAMPLE_OF_THE_R2[0];
        let (_, v2) = THE_WORKED_EXAMPLE_OF_THE_R2[1];
        let at_one_position =
            block_of_the_chromosomes(&[(0, 1000, v1.as_slice()), (0, 1000, v2.as_slice())], 6, 2);
        let mut filter = LdFilter::new(0.5, 1).unwrap();
        let kept = kept_of_the_blocks(&mut filter, vec![at_one_position]).expect("the block");
        assert_eq!(kept, [1000]);
    }

    /// A variant is in the window of a variant exactly `max_dist` base
    /// pairs ahead of it and out of the window of the next base pair: v2 is
    /// 0.675 against v1 and goes at 1000 bp when the two are 1000 apart,
    /// and stays when they are 1001 apart.
    #[test]
    fn a_variant_leaves_the_window_when_the_filter_passes_max_dist_beyond_it() {
        let (_, v1) = THE_WORKED_EXAMPLE_OF_THE_R2[0];
        let (_, v2) = THE_WORKED_EXAMPLE_OF_THE_R2[1];
        let with_the_second_at = |pos| {
            let block = block_of_the_chromosomes(
                &[(0, 1000, v1.as_slice()), (0, pos, v2.as_slice())],
                6,
                2,
            );
            let mut filter = LdFilter::new(0.5, 1000).unwrap();
            kept_of_the_blocks(&mut filter, vec![block]).expect("the block")
        };
        // 1000 bp behind the second variant and 1001 bp behind it.
        assert_eq!(with_the_second_at(2000), [1000]);
        assert_eq!(with_the_second_at(2001), [1000, 2001]);
    }

    /// A position that falls below the position of the variant before it on
    /// its chromosome is refused, and the block, the counts and the window
    /// are as they were: the block that comes after the refused one is
    /// filtered against the window of the blocks before it.
    #[test]
    fn a_position_that_falls_within_a_chromosome_is_refused() {
        let (_, v1) = THE_WORKED_EXAMPLE_OF_THE_R2[0];
        let (_, v2) = THE_WORKED_EXAMPLE_OF_THE_R2[1];
        let mut filter = LdFilter::new(0.5, 5000).unwrap();
        let mut first = block_of_the_chromosomes(&[(0, 1000, v1.as_slice())], 6, 2);
        filter.filter_block(&mut first).expect("the first block");

        let mut out_of_order =
            block_of_the_chromosomes(&[(0, 3000, v2.as_slice()), (0, 2000, v2.as_slice())], 6, 2);
        let error = filter
            .filter_block(&mut out_of_order)
            .expect_err("the block was refused");
        assert!(
            matches!(
                error,
                Error::LdFilterVariantOutOfOrder {
                    variant: 3,
                    chrom: TheChromOfTheVariant::Numbered(0),
                    problem: TheOrderOfTheVariants::ThePositionFalls {
                        pos: 2000,
                        pos_before: 3000,
                    },
                }
            ),
            "{error}"
        );
        assert!(error.to_string().contains("the variant 3"), "{error}");
        // The filter is given a block, which holds the number of a
        // chromosome and not its name, so the message names the number; a
        // reader over a source puts the name there.
        assert!(
            error.to_string().contains("the chromosome numbered 0"),
            "{error}"
        );
        assert!(
            error.to_string().contains("`bcftools sort` writes"),
            "{error}"
        );
        assert_eq!(out_of_order.num_vars, 2);
        assert_eq!(positions_of(&out_of_order), [3000, 2000]);
        assert_eq!(filter.stats(), pair(1, 1));

        // v1 is still in the window, and v2 is 0.675 against it.
        let mut after_it = block_of_the_chromosomes(&[(0, 2000, v2.as_slice())], 6, 2);
        filter.filter_block(&mut after_it).expect("the block after");
        assert_eq!(after_it.num_vars, 0);
        assert_eq!(filter.stats(), pair(2, 1));
    }

    /// A position that falls below the last variant of the block before it
    /// is refused too: the filter carries where it has got to from one
    /// block to the next, as it carries the window.
    #[test]
    fn a_position_that_falls_below_the_last_variant_of_the_block_before_is_refused() {
        let (_, v1) = THE_WORKED_EXAMPLE_OF_THE_R2[0];
        let (_, v2) = THE_WORKED_EXAMPLE_OF_THE_R2[1];
        let mut filter = LdFilter::new(0.5, 5000).unwrap();
        let mut first = block_of_the_chromosomes(&[(0, 3000, v1.as_slice())], 6, 2);
        filter.filter_block(&mut first).expect("the first block");
        let mut goes_back = block_of_the_chromosomes(&[(0, 2000, v2.as_slice())], 6, 2);
        let error = filter
            .filter_block(&mut goes_back)
            .expect_err("the block was refused");
        assert!(
            matches!(
                error,
                Error::LdFilterVariantOutOfOrder {
                    variant: 2,
                    chrom: TheChromOfTheVariant::Numbered(0),
                    problem: TheOrderOfTheVariants::ThePositionFalls {
                        pos: 2000,
                        pos_before: 3000,
                    },
                }
            ),
            "{error}"
        );
    }

    /// A chromosome that comes back is refused, whatever the positions: a
    /// variant of it sits behind variants of another chromosome, so the
    /// window of every variant after it would hold what the filter has
    /// already given away.
    #[test]
    fn a_chromosome_that_comes_back_is_refused() {
        let (_, v1) = THE_WORKED_EXAMPLE_OF_THE_R2[0];
        let mut filter = LdFilter::new(0.5, 5000).unwrap();
        let mut blocks = block_of_the_chromosomes(
            &[
                (0, 1000, v1.as_slice()),
                (1, 1000, v1.as_slice()),
                (0, 2000, v1.as_slice()),
            ],
            6,
            2,
        );
        let error = filter
            .filter_block(&mut blocks)
            .expect_err("the block was refused");
        assert!(
            matches!(
                error,
                Error::LdFilterVariantOutOfOrder {
                    variant: 3,
                    chrom: TheChromOfTheVariant::Numbered(0),
                    problem: TheOrderOfTheVariants::TheChromosomeCameBack {
                        pos: 2000,
                        pos_before: 1000,
                    },
                }
            ),
            "{error}"
        );
        assert_eq!(filter.stats(), FilteringStats::default());
    }

    /// A `max_allowed_r2` that is not a number from 0 to 1 and a `max_dist`
    /// below 1 are refused where the filter is built, and the error names
    /// the argument and the value.
    #[test]
    fn the_arguments_of_the_filter_are_refused_where_it_is_built() {
        for threshold in [-0.1, 1.5, f64::NAN] {
            let error = LdFilter::new(threshold, 5000).expect_err("the threshold was refused");
            assert!(
                matches!(
                    error,
                    Error::VarFilterThresholdOutOfRange { kind: "ld", .. }
                ),
                "{error}"
            );
        }
        let error = LdFilter::new(0.5, 0).expect_err("the distance was refused");
        assert!(
            matches!(error, Error::LdFilterMaxDistTooSmall { max_dist: 0 }),
            "{error}"
        );
        assert!(LdFilter::new(0.0, 1).is_ok());
        assert!(LdFilter::new(1.0, u64::MAX).is_ok());
    }

    /// The filter asks its source for the position besides the genotypes, so
    /// a block with variants and no position, as one with variants and no
    /// genotypes, is the error of a field that is not in the block.
    #[test]
    fn a_block_with_variants_and_no_position_or_no_genotypes_is_refused() {
        let mut filter = LdFilter::new(0.5, 5000).unwrap();
        let mut with_no_position = block_of_the_r2_example(&[0, 1]);
        with_no_position.chrom = None;
        with_no_position.pos = None;
        let error = filter
            .filter_block(&mut with_no_position)
            .expect_err("the block was refused");
        assert!(
            matches!(
                error,
                Error::FieldsNotInTheBlock {
                    fields: Needs::CHROM_POS
                }
            ),
            "{error}"
        );

        let mut with_no_genotypes = block_of_the_r2_example(&[0, 1]);
        with_no_genotypes.gts = Vec::new();
        with_no_genotypes.num_individuals = 0;
        let error = filter
            .filter_block(&mut with_no_genotypes)
            .expect_err("the block was refused");
        assert!(
            matches!(error, Error::FieldsNotInTheBlock { fields: Needs::GTS }),
            "{error}"
        );
        assert_eq!(filter.stats(), FilteringStats::default());
    }

    /// The counts are of every block the filter was given, a filter just
    /// built has counted nothing, and a block of no variants is left as it
    /// is and counted as nothing.
    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "the threshold comes back as it was written and is not the end of any arithmetic"
    )]
    fn the_counts_of_the_filter_add_up_over_the_blocks_it_was_given() {
        let mut filter = LdFilter::new(0.5, 5000).unwrap();
        assert_eq!(filter.stats(), FilteringStats::default());
        assert_eq!(filter.max_allowed_r2(), 0.5);
        assert_eq!(filter.max_dist(), 5000);

        let mut empty = block_of_the_r2_example(&[]);
        filter.filter_block(&mut empty).expect("the empty block");
        assert_eq!(filter.stats(), FilteringStats::default());

        let blocks = vec![
            block_of_the_r2_example(&[0, 1]),
            block_of_the_r2_example(&[2, 3]),
            block_of_the_r2_example(&[4]),
        ];
        let kept = kept_of_the_blocks(&mut filter, blocks).expect("the blocks");
        assert_eq!(kept, [1000, 5000]);
        assert_eq!(filter.stats(), pair(5, 2));
    }

    /// The reader of the filter by linkage disequilibrium, of the worked
    /// example of the r², over the blocks a reader of the tests gives.
    fn ld_reader_of(blocks: Vec<Block>) -> LdFilteredReader<GivenBlocks> {
        LdFilteredReader::new(
            GivenBlocks::of_the_r2_example(blocks),
            LdFilter::new(0.5, 5000).expect("the filter"),
        )
        .expect("the reader over the blocks")
    }

    /// A block that the filter emptied is not given: the reader takes the
    /// next block of its source, and the variants of the block it dropped
    /// are in its counts. Here the first and the last block hold a variant
    /// of one dosage, which the filter drops wherever it is, and the middle
    /// one v5, which has no variant of its window to be compared with.
    #[test]
    fn a_block_the_ld_filter_left_with_no_variant_is_not_given_and_the_next_one_is_taken() {
        let (_, one_dosage) = THE_WORKED_EXAMPLE_OF_THE_R2[3];
        let mut filtered = ld_reader_of(vec![
            block_of_the_chromosomes(&[(0, 1000, one_dosage.as_slice())], 6, 2),
            block_of_the_r2_example(&[4]),
            block_of_the_chromosomes(&[(0, 6000, one_dosage.as_slice())], 6, 2),
        ]);

        let blocks = blocks_of(&mut filtered).expect("the blocks");

        assert_eq!(blocks.len(), 1);
        assert_eq!(positions_of_blocks(&blocks), [5000]);
        assert_eq!(filtered.filtering_stats(), vec![("ld", pair(3, 1))]);
        // And there is no block after the last one.
        assert!(filtered.next_block().expect("no more blocks").is_none());
    }

    /// After an error of its source the reader gives no block and does not
    /// call its source again, as every reader over a reader of
    /// `docs/specs/block.md` does.
    #[test]
    fn after_an_error_of_the_source_the_ld_reader_gives_no_block_and_does_not_call_it_again() {
        let source = GivenBlocks {
            fails_at: Some(2),
            ..GivenBlocks::of_the_r2_example(vec![
                block_of_the_r2_example(&[0]),
                block_of_the_r2_example(&[4]),
            ])
        };
        let calls = source.calls();
        let mut filtered =
            LdFilteredReader::new(source, LdFilter::new(0.5, 5000).unwrap()).unwrap();

        let first = filtered.next_block().expect("the first block");
        assert_eq!(first.map(|block| positions_of(&block)), Some(vec![1000]));
        let error = filtered.next_block().expect_err("the source failed");
        assert!(matches!(error, Error::Io(_)), "{error}");

        assert!(
            filtered
                .next_block()
                .expect("no block after the error")
                .is_none()
        );
        assert_eq!(calls.load(Ordering::SeqCst), 2);
    }

    /// After an error of its filter the reader gives no block and does not
    /// call its source again either: the variants that follow the one the
    /// filter refused would otherwise come out as if nothing had happened,
    /// and the window of every variant after it is the one the wrong order
    /// broke.
    #[test]
    fn after_an_error_of_the_ld_filter_there_is_no_block_and_the_source_is_not_called_again() {
        let (_, v1) = THE_WORKED_EXAMPLE_OF_THE_R2[0];
        let (_, v2) = THE_WORKED_EXAMPLE_OF_THE_R2[1];
        let source = GivenBlocks::of_the_r2_example(vec![
            block_of_the_chromosomes(&[(0, 3000, v1.as_slice())], 6, 2),
            block_of_the_chromosomes(&[(0, 2000, v2.as_slice())], 6, 2),
        ]);
        let calls = source.calls();
        let mut filtered =
            LdFilteredReader::new(source, LdFilter::new(0.5, 5000).unwrap()).unwrap();

        let first = filtered.next_block().expect("the first block");
        assert_eq!(first.map(|block| positions_of(&block)), Some(vec![3000]));
        let error = filtered.next_block().expect_err("the block was refused");
        assert!(
            matches!(
                error,
                Error::LdFilterVariantOutOfOrder {
                    variant: 2,
                    chrom: TheChromOfTheVariant::Named(_),
                    problem: TheOrderOfTheVariants::ThePositionFalls {
                        pos: 2000,
                        pos_before: 3000,
                    },
                }
            ),
            "{error}"
        );
        // The reader has the table of its source, so the message names the
        // chromosome as the file does and not by its number.
        assert!(error.to_string().contains("the chromosome chr1"), "{error}");

        assert!(
            filtered
                .next_block()
                .expect("no block after the error")
                .is_none()
        );
        assert_eq!(calls.load(Ordering::SeqCst), 2);
    }

    /// A source that gives a block of no variants has a defect and is the
    /// error of that, and it is not asked again: over a source that always
    /// gives one, a reader that asked again would never come back.
    #[test]
    fn a_source_that_gives_the_ld_reader_a_block_of_no_variants_is_an_error_and_is_not_called_again()
     {
        let source = GivenBlocks::of_the_r2_example(vec![
            block_of_the_r2_example(&[]),
            block_of_the_r2_example(&[0]),
        ]);
        let calls = source.calls();
        let mut filtered =
            LdFilteredReader::new(source, LdFilter::new(0.5, 5000).unwrap()).unwrap();

        let error = filtered.next_block().expect_err("the block was refused");
        assert!(
            matches!(error, Error::ReaderGaveABlockOfNoVariants),
            "{error}"
        );
        assert!(
            filtered
                .next_block()
                .expect("no block after the error")
                .is_none()
        );
        assert_eq!(calls.load(Ordering::SeqCst), 1);
    }

    /// The reader has no individuals, no ploidy and no chromosome names of
    /// its own: it gives those of its source, whose ploidy here is 4 and
    /// not the 2 of the blocks of the other tests. A filter that has taken
    /// no block yet counts 0 given and 0 kept.
    #[test]
    fn the_ld_reader_gives_the_individuals_the_ploidy_and_the_chromosomes_of_its_source() {
        let source = GivenBlocks {
            ploidy: 4,
            ..GivenBlocks::of_the_r2_example(Vec::new())
        };
        let filtered = LdFilteredReader::new(source, LdFilter::new(0.5, 5000).unwrap()).unwrap();

        assert_eq!(filtered.individuals().len(), 6);
        assert_eq!(
            filtered.individuals().first().map(String::as_str),
            Some("ind1")
        );
        assert_eq!(filtered.ploidy(), 4);
        assert_eq!(filtered.chroms().name(0), Some("chr1"));
        assert_eq!(filtered.filtering_stats(), vec![("ld", pair(0, 0))]);
    }

    /// The filter reads the genotypes, the chromosome and the position of
    /// every variant, so the blocks the reader gives hold the three also
    /// when the consumer asked for none of them, and they hold no column
    /// that nobody asked for. The consumer here asks for the qualities
    /// alone, over the 500 variants of `many.vcf`.
    #[test]
    fn the_blocks_of_the_ld_reader_hold_the_chromosome_and_the_position_the_filter_reads() {
        let reader = many_vcf_reader(Some(7), Needs::ALL);
        let mut filtered =
            LdFilteredReader::new(reader, LdFilter::new(0.1, 10_000).unwrap()).unwrap();
        filtered.set_needs(Needs::QUAL);

        let blocks = blocks_of(&mut filtered).expect("the blocks");

        assert!(!blocks.is_empty());
        for block in &blocks {
            assert!(
                block
                    .fields()
                    .contains(Needs::GTS | Needs::CHROM_POS | Needs::QUAL)
            );
            assert!(!block.gts.is_empty());
            assert!(block.id.is_none());
            assert!(block.alleles.is_none());
        }
        let stats = filtered.filtering_stats();
        assert_eq!(
            stats
                .first()
                .map(|(kind, stats)| (*kind, stats.vars_processed)),
            Some(("ld", 500))
        );
    }

    /// The three counts of the worked example of the filter by linkage
    /// disequilibrium in "How it is verified" of `docs/specs/filters.md`: a
    /// missing data filter at 1 before it and an observed heterozygosity
    /// filter at 1 after it give 5 and 5, 5 and 2, and 2 and 2.
    ///
    /// Neither of the two filters at 1 drops anything, v3 having one
    /// missing genotype of six and the most heterozygous variant three of
    /// six, so the counts of the three say where each of them stands in the
    /// chain and what the one in the middle took out.
    #[test]
    fn the_chain_of_the_worked_example_gives_the_three_pairs_of_counts_of_the_spec() {
        let source =
            GivenBlocks::of_the_r2_example(vec![block_of_the_r2_example(&[0, 1, 2, 3, 4])]);
        let mut chain = chain_of(
            Box::new(source),
            &[
                MaxMissingRate(1.0),
                MaxLdR2 {
                    max_allowed_r2: 0.5,
                    max_dist: 5000,
                },
                MaxObsHet(1.0),
            ],
        )
        .expect("the chain of the three criteria");

        let blocks = blocks_of(&mut chain).expect("the blocks");

        assert_eq!(positions_of_blocks(&blocks), [1000, 5000]);
        assert_eq!(
            chain.filtering_stats(),
            vec![
                ("obs_het", pair(2, 2)),
                ("ld", pair(5, 2)),
                ("missing_data", pair(5, 5)),
            ]
        );
    }

    /// The chain of a maf filter and the filter by linkage disequilibrium,
    /// which is what a user writes in place of pyNei's one call to
    /// `filter_by_ld_and_maf`: the maf filter at 0.8 drops v4, whose
    /// major allele frequency is 1, and the filter at 0.5 and 5000 bp keeps
    /// v1 and v5 of the four variants left. The counts of both reach the
    /// user apart, the outermost filter first, which is the last step.
    #[test]
    fn a_chain_of_a_maf_filter_and_the_ld_filter_keeps_two_variants_and_gives_the_counts_of_both() {
        let source =
            GivenBlocks::of_the_r2_example(vec![block_of_the_r2_example(&[0, 1, 2, 3, 4])]);
        let mut chain = chain_of(
            Box::new(source),
            &[
                MaxMaf(0.8),
                MaxLdR2 {
                    max_allowed_r2: 0.5,
                    max_dist: 5000,
                },
            ],
        )
        .expect("the chain of the two criteria");

        let blocks = blocks_of(&mut chain).expect("the blocks");

        assert_eq!(positions_of_blocks(&blocks), [1000, 5000]);
        assert_eq!(
            chain.filtering_stats(),
            vec![("ld", pair(4, 2)), ("maf", pair(5, 4))]
        );
    }

    /// A second criterion by linkage disequilibrium is refused where the
    /// chain is built, among the criteria of one call and over a chain that
    /// holds one already, which is what a binding crate builds over a
    /// chain. The message names the `max_allowed_r2` that was written.
    #[test]
    fn chain_of_a_second_criterion_by_linkage_disequilibrium_is_refused() {
        let ld = |max_allowed_r2| MaxLdR2 {
            max_allowed_r2,
            max_dist: 5000,
        };
        let source = Box::new(GivenBlocks::of_the_r2_example(Vec::new()));

        let error = chain_of(source, &[ld(0.1), ld(0.3)])
            .err()
            .expect("the chain was refused");
        let message = error.to_string();
        assert!(
            matches!(error, Error::VarFilterOfAKindThatIsSet { kind: "ld", .. }),
            "{message}"
        );
        assert!(message.contains("0.3"), "{message}");

        let holds_one = chain_of(
            Box::new(GivenBlocks::of_the_r2_example(Vec::new())),
            &[ld(0.1)],
        )
        .expect("the chain of one criterion");
        let error = chain_of(holds_one, &[MaxMaf(0.8), ld(0.3)])
            .err()
            .expect("the chain over it was refused");
        assert!(
            matches!(error, Error::VarFilterOfAKindThatIsSet { kind: "ld", .. }),
            "{error}"
        );
    }

    /// The two arguments of the filter by linkage disequilibrium are
    /// refused where the chain is built, as the threshold of the three
    /// filters that compare one number of a variant is, and the errors are
    /// those of `LdFilter::new`.
    #[test]
    fn chain_of_the_arguments_of_the_ld_filter_out_of_range_is_the_error_of_the_filter() {
        let refused = |criterion| {
            let source = Box::new(GivenBlocks::of_the_r2_example(Vec::new()));
            chain_of(source, &[MaxMissingRate(0.04), criterion])
                .err()
                .expect("the chain was refused")
        };

        let error = refused(MaxLdR2 {
            max_allowed_r2: 1.5,
            max_dist: 5000,
        });
        assert!(
            matches!(
                error,
                Error::VarFilterThresholdOutOfRange { kind: "ld", .. }
            ),
            "{error}"
        );

        let error = refused(MaxLdR2 {
            max_allowed_r2: 0.5,
            max_dist: 0,
        });
        assert!(
            matches!(error, Error::LdFilterMaxDistTooSmall { max_dist: 0 }),
            "{error}"
        );
    }

    /// The kind of the filter by linkage disequilibrium is looked for among
    /// the criteria that are set as every other kind is, which is what a
    /// binding crate calls when a user adds the filter to a `Variants`, and
    /// the message names both thresholds: the one that is set and the one
    /// that was written.
    #[test]
    fn refuse_a_second_filter_of_a_kind_refuses_a_second_one_by_linkage_disequilibrium() {
        let set = [
            MaxMaf(0.8),
            MaxLdR2 {
                max_allowed_r2: 0.1,
                max_dist: 10_000,
            },
        ];

        assert!(refuse_a_second_filter_of_a_kind(&set, MaxObsHet(0.5)).is_ok());
        let error = refuse_a_second_filter_of_a_kind(
            &set,
            MaxLdR2 {
                max_allowed_r2: 0.3,
                max_dist: 50_000,
            },
        )
        .expect_err("the second filter by linkage disequilibrium was refused");

        let message = error.to_string();
        assert!(
            matches!(
                error,
                Error::VarFilterOfAKindThatIsSet {
                    kind: "ld",
                    threshold_that_is_set: Some(_),
                    ..
                }
            ),
            "{message}"
        );
        assert!(message.contains("0.3"), "{message}");
        assert!(message.contains("0.1"), "{message}");
    }

    /// The four rows of the table of "How it is verified" of the item "The
    /// filter by linkage disequilibrium" of `docs/specs/filters.md`, run on
    /// `tests/reference/ld/ld.vcf.gz`: the window in base pairs, the
    /// largest r² a kept variant may have against a variant of its window,
    /// that threshold as `tests/reference/ld/ld.filtered.tsv` writes it,
    /// the variants kept of the 500, and the first five kept by position.
    const THE_TABLE_OF_THE_LD_FILTER: [(u64, f64, &str, u64, [&str; 5]); 4] = [
        (
            10_000,
            0.1,
            "0.1",
            84,
            [
                "chr1:1000",
                "chr1:10000",
                "chr1:16000",
                "chr1:22000",
                "chr1:27000",
            ],
        ),
        (
            10_000,
            0.3,
            "0.3",
            133,
            [
                "chr1:1000",
                "chr1:5000",
                "chr1:7000",
                "chr1:11000",
                "chr1:15000",
            ],
        ),
        (
            50_000,
            0.3,
            "0.3",
            85,
            [
                "chr1:1000",
                "chr1:5000",
                "chr1:7000",
                "chr1:11000",
                "chr1:15000",
            ],
        ),
        (
            250_000,
            0.3,
            "0.3",
            85,
            [
                "chr1:1000",
                "chr1:5000",
                "chr1:7000",
                "chr1:11000",
                "chr1:15000",
            ],
        ),
    ];

    /// The chromosome and the position of each variant that a filter
    /// kept, in the order the blocks gave them.
    type TheVariantsKept = Vec<(String, u64)>;

    /// `tests/reference/ld/`, where `make_reference.py` writes the dataset
    /// of `docs/specs/ld.md` and the numbers plink2 gives for it, at the
    /// root of the repository and not inside this crate.
    fn the_ld_reference_path(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/ld")
            .join(name)
    }

    /// A reader over `tests/reference/ld/ld.vcf.gz`, the 500 variants of
    /// 100 diploid individuals on two chromosomes of `docs/specs/ld.md`,
    /// read as plink2 read them, with the variants that failed their FILTER
    /// among them, in blocks of `num_vars_per_block` variants.
    fn the_ld_dataset(num_vars_per_block: Option<usize>) -> VcfReader<BufReader<File>> {
        let path = the_ld_reference_path("ld.vcf.gz");
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block,
        };
        VcfReader::from_path(&path, options)
            .unwrap_or_else(|error| panic!("{path}: {error}", path = path.display()))
    }

    /// The chromosome and the position of every variant that
    /// `tests/reference/ld/ld.filtered.tsv` holds for that setting, in the
    /// order of the file. `make_reference.py` works that set out from the
    /// r² that plink2 wrote for `ld.vcf.gz` and from the rule of the spec,
    /// so it says which variants are kept without any arithmetic of popnei,
    /// and it is where the three properties of "How it is verified" were
    /// checked, in `tests/reference/ld/ld.filter.properties.txt`.
    fn the_kept_of_the_reference(max_dist: u64, max_allowed_r2: &str) -> TheVariantsKept {
        let path = the_ld_reference_path("ld.filtered.tsv");
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("{path}: {error}", path = path.display()));
        let of_the_dist = max_dist.to_string();
        let mut kept = Vec::new();
        for line in text.lines().skip(1) {
            let fields: Vec<&str> = line.split('\t').collect();
            let [dist, r2, chrom, pos] = fields[..] else {
                panic!("`{line}` has not the four columns of ld.filtered.tsv");
            };
            if dist != of_the_dist || r2 != max_allowed_r2 {
                continue;
            }
            let pos = pos
                .parse()
                .unwrap_or_else(|error| panic!("`{line}`: {error}"));
            kept.push((chrom.to_owned(), pos));
        }
        assert!(
            !kept.is_empty(),
            "ld.filtered.tsv holds no variant of {max_dist} bp and {max_allowed_r2}"
        );
        kept
    }

    /// The chromosome and the position of every variant that the filter by
    /// linkage disequilibrium keeps of `ld.vcf.gz`, read in blocks of
    /// `num_vars_per_block` variants, with the counts of the filter after
    /// the last block.
    fn the_kept_of_the_ld_dataset(
        max_allowed_r2: f64,
        max_dist: u64,
        num_vars_per_block: Option<usize>,
    ) -> (TheVariantsKept, Vec<(&'static str, FilteringStats)>) {
        let filter = LdFilter::new(max_allowed_r2, max_dist).expect("the filter");
        let mut filtered = LdFilteredReader::new(the_ld_dataset(num_vars_per_block), filter)
            .expect("the reader over the VCF");
        filtered.set_needs(Needs::GTS | Needs::CHROM_POS);
        let blocks = blocks_of(&mut filtered).expect("the blocks");
        let mut kept = Vec::new();
        for block in &blocks {
            assert!(block.num_vars > 0, "a block with no variant was given");
            assert!(
                block.check().is_ok(),
                "a block whose arrays are not of its size"
            );
            let chroms = block.chrom.as_ref().expect("the chromosomes of the block");
            let positions = block.pos.as_ref().expect("the positions of the block");
            for (number, position) in chroms.iter().zip(positions) {
                let name = filtered
                    .chroms()
                    .name(*number)
                    .unwrap_or_else(|| panic!("the chromosome {number} has no name"));
                kept.push((name.to_owned(), *position));
            }
        }
        (kept, filtered.filtering_stats())
    }

    /// The four rows of the table of the spec, in blocks of 7 variants, of
    /// 64 and of the size the VCF reader chooses: each keeps the variants
    /// the table gives, which are the ones the rule of the spec keeps when
    /// it is run over plink2's r² in `make_reference.py`, and the counts of
    /// the filter are the 500 variants of the file and the ones it kept.
    /// The first variant kept of chr2 is `chr2:1000` at every setting,
    /// which is where the window stops at the end of a chromosome.
    #[test]
    fn the_ld_filter_keeps_the_variants_of_ld_vcf_that_the_table_of_the_spec_gives() {
        for (max_dist, max_allowed_r2, of_the_reference, kept_of_500, first_five) in
            THE_TABLE_OF_THE_LD_FILTER
        {
            let of_the_reference = the_kept_of_the_reference(max_dist, of_the_reference);
            for num_vars_per_block in [Some(7), Some(64), None] {
                let (kept, stats) =
                    the_kept_of_the_ld_dataset(max_allowed_r2, max_dist, num_vars_per_block);
                let what = format!(
                    "at {max_dist} bp and {max_allowed_r2}, in blocks of {num_vars_per_block:?}"
                );

                assert_eq!(
                    u64::try_from(kept.len()).expect("the variants kept"),
                    kept_of_500,
                    "{what}"
                );
                let five: Vec<String> = kept
                    .iter()
                    .take(5)
                    .map(|(chrom, pos)| format!("{chrom}:{pos}"))
                    .collect();
                assert_eq!(five, first_five.map(String::from).to_vec(), "{what}");
                assert_eq!(
                    kept.iter().find(|(chrom, _)| chrom == "chr2"),
                    Some(&("chr2".to_owned(), 1000)),
                    "{what}"
                );
                assert_eq!(kept, of_the_reference, "{what}");
                assert_eq!(stats, vec![("ld", pair(500, kept_of_500))], "{what}");
            }
        }
    }

    /// The filter keeps the same variants of `ld.vcf.gz` on a pool of one
    /// thread and on one of four, at 10000 bp and 0.3, where it keeps 133
    /// of the 500.
    ///
    /// The r² of the filter is worked out by the products of the linear
    /// algebra, which on the faer backend run on the pool that is
    /// installed, so the kept set is asserted on two pools as the matrix of
    /// `docs/specs/ld.md` and the three filters that compare one number
    /// are. The pools are built here and are not rayon's global one, which
    /// has one thread per core of the machine; rayon is a dependency of the
    /// targets that are not wasm, so this test is compiled for those alone.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_variants_the_ld_filter_keeps_are_the_same_on_one_thread_and_on_several() {
        let kept = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            pool.install(|| the_kept_of_the_ld_dataset(0.3, 10_000, Some(64)))
        };
        let (on_one, counts_of_one) = kept(1);
        let (on_four, counts_of_four) = kept(4);
        assert_eq!(on_one.len(), 133);
        assert_eq!(on_one, on_four);
        assert_eq!(counts_of_one, counts_of_four);
    }
}
