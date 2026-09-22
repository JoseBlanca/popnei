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
//! A filter of a pass over the variants is a reader over another reader,
//! [`FilteredReader`], and several filters are several of them, one over
//! the other, in the order in which the user put them on. What the user put
//! on is a [`PassStep`], and [`chain_of`] builds the chain from the steps of
//! one pass: it is what each binding crate calls when a pass starts.
//!
//! One step of a pass is not a threshold: the filter of individuals keeps,
//! of every variant, the genotypes of the individuals a user named, in the
//! order they named them, and drops those of the rest. It is
//! [`IndividualsReader`], a reader over another reader as well, with
//! [`resolve_individuals`] turning the names a user wrote into the indices
//! among the individuals of the source that it compacts each block by. It
//! takes no variant out and has no counts, and a threshold filter after it
//! in the steps counts over the kept individuals alone.
//!
//! `docs/specs/filters.md` has the design, and the row `filters` of section
//! 9 of `docs/architecture.md` where the module sits.

use std::collections::{HashMap, HashSet};
use std::fmt;

use crate::block::{Block, BlockReader};
use crate::error::{Error, Result};
use crate::variant::{AlleleCounts, ChromTable, Needs, count_alleles, count_gts};

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

    /// The largest value of the number of a variant that keeps it,
    /// whichever of the three numbers this criterion compares.
    ///
    /// A binding crate reads it for the arguments of the step it shows the
    /// user, `{"max_allowed_maf": 0.95}`.
    #[must_use]
    pub fn threshold(&self) -> f64 {
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

/// One step of a pass over the variants: what every pass built from a
/// `Variants` does to the variants it reads, in the order in which the user
/// put the steps on.
///
/// Each binding crate keeps the steps of its `Variants` as a list of these,
/// and [`chain_of`] builds the readers of one pass from that list: which
/// reader a step becomes, and in which order, is of the filters and not of
/// Python or of TypeScript.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum PassStep {
    /// The variants whose number is at most the threshold of the criterion
    /// are kept, and the others are left out of every block of the pass.
    VarFilter(VarFilteringCriterion),
    /// The names of the individuals to keep, in the order to keep them:
    /// every variant stays, and of each one the genotypes of these
    /// individuals alone go on.
    KeepIndividuals(Vec<String>),
}

impl PassStep {
    /// `"missing_data"`, `"maf"`, `"obs_het"` or `"individuals"`: the name
    /// the step has for a Python and a TypeScript user, under which the
    /// counts of a filter reach them and by which a second step of the same
    /// kind is refused.
    #[must_use]
    pub fn kind(&self) -> &'static str {
        match self {
            PassStep::VarFilter(criterion) => criterion.kind(),
            PassStep::KeepIndividuals(_) => "individuals",
        }
    }
}

/// One reader over `reader` for each step, in their order, so that each
/// step sees what the one before it gave: the chain of one pass. A
/// [`PassStep::VarFilter`] becomes a [`FilteredReader`], and no step gives
/// `reader` as it is.
///
/// Both binding crates build the chain of a pass with this, and neither
/// writes the loop: in which order the steps go, and what comes out while
/// the readers are built, are of the filters and not of Python or of
/// TypeScript.
///
/// What it gives is the outermost reader of the chain, which whoever started
/// the pass holds: they read [`BlockReader::filtering_stats`] from it when
/// the consumer returns, and they lend it, `&mut`, to a consumer that takes
/// a reader. The fields the consumer wants are set on it once the chain is
/// built, which [`FilteredReader::new`] says why.
///
/// # Errors
///
/// What [`VarFilter::new`] refuses, a threshold that is NaN, below 0 or
/// above 1, and what [`FilteredReader::new`] refuses, a threshold filter of
/// the kind of one before it in `steps` or of a filter that `reader` holds
/// already, which a chain built over a chain has. What
/// [`IndividualsReader::new`] refuses, a name that is not an individual of
/// what the step is put on, a name that is there twice and no name at all.
/// And a second [`PassStep::KeepIndividuals`] among `steps`, which the
/// chain has to find itself: the filter of individuals takes no variant
/// out, so it has no counts and a reader cannot be asked whether it holds
/// one. No block was read when any of them comes.
pub fn chain_of(reader: Box<dyn BlockReader>, steps: &[PassStep]) -> Result<Box<dyn BlockReader>> {
    let mut chain = reader;
    for (index, step) in steps.iter().enumerate() {
        match step {
            PassStep::VarFilter(criterion) => {
                chain = Box::new(FilteredReader::new(chain, VarFilter::new(*criterion)?)?);
            }
            PassStep::KeepIndividuals(names) => {
                refuse_a_second_filter_of_a_kind(steps.get(..index).unwrap_or(&[]), step)?;
                chain = Box::new(IndividualsReader::new(chain, names)?);
            }
        }
    }
    Ok(chain)
}

/// The error of a second filter of one kind, when `new` is of the kind of
/// one of `set`, the steps that are set already.
///
/// Two threshold filters of one kind keep the variants that the stricter of
/// the two keeps alone, so the second says that the user has lost track of
/// what their variants carry. Both binding crates call this when a user
/// adds a filter to a `Variants`, where no reader exists yet and the steps
/// are what says which filters are set.
///
/// # Errors
///
/// When a step of `set` has the kind of `new`. For a threshold filter the
/// error carries both thresholds, the one of `new` and the one that is set,
/// where the same error from [`FilteredReader::new`] carries the first
/// alone: a chain of readers says which kinds of filter it holds and not
/// with which thresholds. For the filter of individuals it carries the
/// kind, since a list of individuals has no number to name it by, and two
/// lists keep the individuals that are in both, which is one list.
pub fn refuse_a_second_filter_of_a_kind(set: &[PassStep], new: &PassStep) -> Result<()> {
    let criterion = match new {
        PassStep::VarFilter(criterion) => criterion,
        PassStep::KeepIndividuals(_) => {
            return match set
                .iter()
                .any(|step| matches!(step, PassStep::KeepIndividuals(_)))
            {
                true => Err(Error::FilterOfIndividualsThatIsSet { kind: new.kind() }),
                false => Ok(()),
            };
        }
    };
    let kind = criterion.kind();
    let that_is_set = set
        .iter()
        .filter_map(|step| match step {
            PassStep::VarFilter(of_the_step) => Some(of_the_step),
            PassStep::KeepIndividuals(_) => None,
        })
        .find(|of_the_step| of_the_step.kind() == kind);
    if let Some(that_is_set) = that_is_set {
        return Err(Error::VarFilterOfAKindThatIsSet {
            kind,
            threshold: criterion.threshold(),
            threshold_that_is_set: Some(that_is_set.threshold()),
        });
    }
    Ok(())
}

/// The index of each of `names` among `individuals`, in the order of
/// `names`: the individuals a filter of individuals keeps, as the indices
/// into the individuals of the source that [`Block::retain_individuals`]
/// takes.
///
/// Both binding crates call it when a user adds the step, against the
/// individuals of the source, so that the three refusals reach the user at
/// the call they wrote, and [`IndividualsReader::new`] calls it again when
/// a pass builds its chain.
///
/// # Errors
///
/// A name that is not one of `individuals`, a name that is there twice, and
/// no name at all. The first two name the name, which is what the user
/// wrote.
pub fn resolve_individuals(names: &[String], individuals: &[String]) -> Result<Vec<usize>> {
    if names.is_empty() {
        return Err(Error::NoIndividualNamed);
    }
    let of_the_source: HashMap<&str, usize> = individuals
        .iter()
        .enumerate()
        .map(|(individual, name)| (name.as_str(), individual))
        .collect();
    let mut keep = Vec::with_capacity(names.len());
    let mut named = HashSet::with_capacity(names.len());
    for name in names {
        let Some(individual) = of_the_source.get(name.as_str()) else {
            return Err(Error::IndividualNotInTheSource { name: name.clone() });
        };
        if !named.insert(*individual) {
            return Err(Error::IndividualNamedTwice { name: name.clone() });
        }
        keep.push(*individual);
    }
    Ok(keep)
}

/// A reader that gives the blocks of its source with the genotypes of the
/// individuals a user named, in the order they named them, and those of no
/// other individual.
///
/// Every variant of the source comes out, so it takes no variant out and
/// has no counts of its own: [`BlockReader::filtering_stats`] gives those
/// of its source alone. It takes a block of its source at whatever size it
/// comes and compacts it with [`Block::retain_individuals`], so the blocks
/// it gives are the blocks of its source, whose size was worked out from
/// the individuals of the source and not from the kept ones. It allocates
/// no block and keeps nothing from one block to the next.
///
/// A threshold filter before it in the steps counts over every individual
/// of the source and one after it over the kept ones, which is what a
/// user's numbers turn on: the missing data filter at 0 over `many.vcf`
/// keeps 26 of its 500 variants, and 423 of them over three of its 50
/// individuals.
///
/// It keeps the contract of a reader of `docs/specs/block.md`: after an
/// error, of its source or of the compaction, it gives `None` at every call
/// and does not call its source again, and a source that gives a block of
/// no variants has a defect and is the error of that.
pub struct IndividualsReader<R: BlockReader> {
    reader: R,
    /// The index of each kept individual among those of the source, in the
    /// order the user named them.
    keep: Vec<usize>,
    /// The names of the kept individuals, in the same order.
    individuals: Vec<String>,
    /// Whether the source has no more blocks or one of the two, the source
    /// or the compaction, gave an error.
    finished: bool,
}

impl<R: BlockReader> IndividualsReader<R> {
    /// The reader that gives the genotypes of `individuals` of every block
    /// of `reader`.
    ///
    /// # Errors
    ///
    /// What [`resolve_individuals`] refuses against the individuals of
    /// `reader`: a name that is not one of them, a name that is there twice
    /// and no name at all.
    pub fn new(reader: R, individuals: &[String]) -> Result<IndividualsReader<R>> {
        let keep = resolve_individuals(individuals, reader.individuals())?;
        let names = keep
            .iter()
            .filter_map(|individual| reader.individuals().get(*individual))
            .cloned()
            .collect();
        Ok(IndividualsReader {
            reader,
            keep,
            individuals: names,
            finished: false,
        })
    }
}

impl<R: BlockReader> BlockReader for IndividualsReader<R> {
    /// The next block of the source with the genotypes of the kept
    /// individuals alone in it, in the order the user named them. Every
    /// variant of the block stays, with every column it had.
    ///
    /// # Errors
    ///
    /// When the source fails; when a block of the source holds no variant,
    /// which no reader of popnei gives; and what
    /// [`Block::retain_individuals`] refuses, a block whose arrays are not
    /// of its size and a block that has variants and no genotypes. After
    /// any of them there is no block and the source is not called again.
    fn next_block(&mut self) -> Result<Option<Block>> {
        if self.finished {
            return Ok(None);
        }
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
        // A source that gives a block of no variants has a defect, and it
        // is not asked again: over a source that always gives one, a reader
        // that asked again would never come back.
        if block.num_vars == 0 {
            self.finished = true;
            return Err(Error::ReaderGaveABlockOfNoVariants);
        }
        if let Err(error) = block.retain_individuals(&self.keep) {
            self.finished = true;
            return Err(error);
        }
        Ok(Some(block))
    }

    /// The kept individuals, in the order the user named them, which is
    /// what everything after this reader sees as the individuals of the
    /// dataset.
    fn individuals(&self) -> &[String] {
        &self.individuals
    }

    fn ploidy(&self) -> usize {
        self.reader.ploidy()
    }

    /// The table of the source: a reader over another reader has none of
    /// its own.
    fn chroms(&self) -> &ChromTable {
        self.reader.chroms()
    }

    /// The fields of the consumer and the genotypes, which this reader
    /// compacts in every block: so the blocks it gives hold the genotypes
    /// also when the consumer did not ask for them.
    fn set_needs(&mut self, needs: Needs) {
        self.reader.set_needs(needs.union(Needs::GTS));
    }

    /// The counts of the filters between the source and its own source.
    /// This reader adds none: it takes no variant out.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        self.reader.filtering_stats()
    }
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
            if called_alleles == 0 {
                return Ok(false);
            }
            let largest = counts.iter().copied().max().unwrap_or(0);
            f64::from(largest) / f64::from(called_alleles)
        }
        VarFilteringCriterion::MaxObsHet(_) => {
            let gt_counts = count_gts(gts, ploidy)?;
            if gt_counts.called == 0 {
                return Ok(false);
            }
            f64::from(gt_counts.het) / f64::from(gt_counts.called)
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
        FilteredReader, FilteringStats, PassStep, VarFilter, VarFilteringCriterion, chain_of,
        keep_of_the_rows, keep_of_the_rows_one_by_one, refuse_a_second_filter_of_a_kind,
        resolve_individuals,
    };
    use crate::block::{Block, BlockReader};
    use crate::error::{Error, Result};
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{ChromTable, MISSING_ALLELE, Needs};

    use VarFilteringCriterion::{MaxMaf, MaxMissingRate, MaxObsHet};

    /// The steps of the threshold filters of `criteria`, in their order,
    /// which is what `chain_of` and `refuse_a_second_filter_of_a_kind` take.
    fn steps_of(criteria: &[VarFilteringCriterion]) -> Vec<PassStep> {
        criteria
            .iter()
            .map(|criterion| PassStep::VarFilter(*criterion))
            .collect()
    }

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
        assert_eq!(MaxMaf(0.0).kind(), MaxMaf(1.0).kind());
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
        let mut chain = chain_of(
            source,
            &steps_of(&[MaxMissingRate(0.04), MaxMaf(0.8), MaxObsHet(0.5)]),
        )
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
            let error = chain_of(source, &steps_of(criteria))
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
        let with_a_maf_filter =
            chain_of(source, &steps_of(&[MaxMaf(0.8)])).expect("the chain of one");

        let error = chain_of(
            with_a_maf_filter,
            &steps_of(&[MaxMissingRate(0.04), MaxMaf(0.95)]),
        )
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
        let set = steps_of(&[MaxMissingRate(0.04), MaxMaf(0.8)]);

        assert!(
            refuse_a_second_filter_of_a_kind(&set, &PassStep::VarFilter(MaxObsHet(0.5))).is_ok()
        );
        assert!(refuse_a_second_filter_of_a_kind(&[], &PassStep::VarFilter(MaxMaf(0.95))).is_ok());
        // A criterion of another kind between the two changes nothing: the
        // kind is looked for among all of them.
        assert!(
            refuse_a_second_filter_of_a_kind(
                &steps_of(&[MaxMaf(0.8), MaxObsHet(0.5)]),
                &PassStep::VarFilter(MaxMaf(0.95))
            )
            .is_err()
        );
    }

    /// The message names the kind and both thresholds, the one that is set
    /// and the one that was written, which is what a user needs in order to
    /// see which of their cells they ran twice.
    #[test]
    fn refuse_a_second_filter_of_a_kind_names_the_kind_and_both_thresholds() {
        let error = refuse_a_second_filter_of_a_kind(
            &steps_of(&[MaxMaf(0.8)]),
            &PassStep::VarFilter(MaxMaf(0.95)),
        )
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
        let error = chain_of(source, &steps_of(&[MaxMissingRate(0.04), MaxMaf(1.5)]))
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

    /// The three names of "How it is verified" of the filter of individuals
    /// of `docs/specs/filters.md`, in the order a user writes them, which
    /// is not the order of the individuals of `many.vcf`.
    const THE_THREE_NAMES: [&str; 3] = ["ind05", "ind00", "ind49"];

    /// Those names as the filter takes them.
    fn the_three_names() -> Vec<String> {
        THE_THREE_NAMES.map(str::to_owned).to_vec()
    }

    /// The 50 individuals of `many.vcf`, `ind00` to `ind49`.
    fn the_individuals_of_many_vcf() -> Vec<String> {
        many_vcf_reader(None, Needs::GTS).individuals().to_vec()
    }

    /// The three names give the indices 5, 0 and 49, in the order they were
    /// named: `resolve_individuals` gives the individuals in the order of
    /// the argument and not in the order of the source, which is what lets
    /// a user put their populations together.
    #[test]
    fn resolve_individuals_gives_the_index_of_each_name_in_the_order_of_the_names() {
        let individuals = the_individuals_of_many_vcf();
        assert_eq!(individuals.len(), 50);

        let kept = resolve_individuals(&the_three_names(), &individuals)
            .expect("the indices of the three names");

        assert_eq!(kept, [5, 0, 49]);
    }

    /// A name that is not an individual of the source is the error that
    /// names it. pyNei drops it in silence and gives a `Variants` of the
    /// names it did find.
    #[test]
    fn resolve_individuals_refuses_a_name_that_is_not_an_individual() {
        let individuals = the_individuals_of_many_vcf();

        let error = resolve_individuals(&["ind05".to_owned(), "nope".to_owned()], &individuals)
            .expect_err("the name that is not an individual was refused");

        let message = error.to_string();
        let Error::IndividualNotInTheSource { ref name } = error else {
            panic!("the error is {message}");
        };
        assert_eq!(name, "nope");
        assert!(message.contains("nope"), "{message}");
    }

    /// A name that is there twice is the error that names it: the same
    /// individual kept twice would be two columns of one individual's
    /// genotypes. pyNei keeps it once.
    #[test]
    fn resolve_individuals_refuses_a_name_that_is_there_twice() {
        let individuals = the_individuals_of_many_vcf();
        let names = ["ind49".to_owned(), "ind05".to_owned(), "ind49".to_owned()];

        let error =
            resolve_individuals(&names, &individuals).expect_err("the name twice was refused");

        let message = error.to_string();
        let Error::IndividualNamedTwice { ref name } = error else {
            panic!("the error is {message}");
        };
        assert_eq!(name, "ind49");
        assert!(message.contains("ind49"), "{message}");
    }

    /// No name at all is the error: the variants of nobody are not a
    /// dataset popnei holds.
    #[test]
    fn resolve_individuals_refuses_no_name_at_all() {
        let individuals = the_individuals_of_many_vcf();

        let error = resolve_individuals(&[], &individuals).expect_err("no name was refused");

        assert!(
            matches!(error, Error::NoIndividualNamed),
            "the error is {error}"
        );
    }

    /// The tests of the reader of the filter of individuals. The module is
    /// named after the type, and not `individuals_reader`, so that
    /// `cargo test -- IndividualsReader` runs them.
    #[expect(
        non_snake_case,
        reason = "the module is named after the type it tests, IndividualsReader, so that \
                  the tests of the reader of the filter of individuals are the ones \
                  cargo test -- IndividualsReader runs"
    )]
    mod IndividualsReader {
        use std::fs::File;
        use std::io::BufReader;

        // The struct is named as the module is, so it is taken from the
        // crate and not from the parent, where the two names would be one
        // name of the type namespace.
        use crate::filters::IndividualsReader;
        use crate::io::vcf::VcfReader;

        use super::{
            Block, BlockReader, Error, FilteredReader, GivenBlocks, MaxMaf, MaxMissingRate, Needs,
            PassStep, VarFilter, block_of_the_worked_example, blocks_of, chain_of, many_vcf_reader,
            pair, positions_of_blocks, refuse_a_second_filter_of_a_kind, steps_of, the_three_names,
        };

        /// The reader of the three individuals of "How it is verified" over
        /// `many.vcf`, with the blocks of the source at `num_vars_per_block`
        /// variants and asked for `needs`.
        fn of_the_three_names(
            num_vars_per_block: Option<usize>,
            needs: Needs,
        ) -> IndividualsReader<VcfReader<BufReader<File>>> {
            let source = many_vcf_reader(num_vars_per_block, needs);
            IndividualsReader::new(source, &the_three_names()).expect("the reader of the three")
        }

        /// The genotypes of the blocks, one row per variant.
        fn rows_of(blocks: &[Block]) -> Vec<Vec<i8>> {
            blocks
                .iter()
                .flat_map(|block| {
                    let width = block
                        .num_individuals
                        .max(1)
                        .saturating_mul(block.ploidy.max(1));
                    block
                        .gts
                        .chunks(width)
                        .map(<[i8]>::to_vec)
                        .collect::<Vec<Vec<i8>>>()
                })
                .collect()
        }

        /// Every variant of `many.vcf` comes out, of the three individuals
        /// alone, named in the order of the argument, and the reader has no
        /// counts of its own: it takes no variant out.
        #[test]
        fn the_blocks_hold_the_three_individuals_in_the_order_of_the_argument() {
            let mut reader = of_the_three_names(Some(7), Needs::GTS | Needs::CHROM_POS);

            let blocks = blocks_of(&mut reader).expect("the blocks");

            assert_eq!(reader.individuals(), the_three_names().as_slice());
            assert_eq!(reader.ploidy(), 2);
            assert!(reader.filtering_stats().is_empty());
            assert_eq!(positions_of_blocks(&blocks).len(), 500);
            for block in &blocks {
                assert_eq!(block.num_individuals, 3);
                block.check().expect("the block is of its size");
            }
        }

        /// The genotypes of each kept individual are the column of the
        /// source at every variant: the same file read whole gives, at the
        /// individuals 5, 0 and 49, what the blocks of the reader give at
        /// the three columns. The two literals are the genotypes "How it is
        /// verified" gives at the positions 1000 and 1074, `1|1 1/1 1/1`
        /// and `0/1 2|1 1|2`.
        #[test]
        fn the_genotypes_of_the_three_are_the_columns_of_the_source() {
            let mut whole = many_vcf_reader(Some(7), Needs::GTS);
            let of_the_source = rows_of(&blocks_of(&mut whole).expect("the blocks of the source"));
            let mut reader = of_the_three_names(Some(7), Needs::GTS);

            let kept = rows_of(&blocks_of(&mut reader).expect("the blocks"));

            assert_eq!(kept.len(), 500);
            assert_eq!(of_the_source.len(), 500);
            for (row, of_the_source) in kept.iter().zip(&of_the_source) {
                let gathered: Vec<i8> = [5_usize, 0, 49]
                    .iter()
                    .flat_map(|individual| {
                        let start = individual.saturating_mul(2);
                        of_the_source
                            .get(start..start.saturating_add(2))
                            .unwrap_or_default()
                            .to_vec()
                    })
                    .collect();
                assert_eq!(*row, gathered);
            }
            // The first variant, at the position 1000, and the one at 1074.
            assert_eq!(
                kept.first().map(Vec::as_slice),
                Some([1, 1, 1, 1, 1, 1].as_slice())
            );
            assert_eq!(
                kept.get(2).map(Vec::as_slice),
                Some([0, 1, 2, 1, 1, 2].as_slice())
            );
        }

        /// The missing data filter at 0 over the three individuals keeps the
        /// 423 variants that bcftools 1.24 and pyNei keep, the first five at
        /// the positions of the spec, with the 500 variants it was given in
        /// its counts.
        #[test]
        fn the_missing_data_filter_over_the_three_keeps_the_423_variants() {
            let reader = of_the_three_names(Some(7), Needs::GTS | Needs::CHROM_POS);
            let filter = VarFilter::new(MaxMissingRate(0.0)).expect("the filter");
            let mut filtered = FilteredReader::new(reader, filter).expect("the filter over it");

            let blocks = blocks_of(&mut filtered).expect("the blocks");

            let positions = positions_of_blocks(&blocks);
            assert_eq!(positions.len(), 423);
            assert_eq!(
                positions.get(..5),
                Some([1000, 1037, 1074, 1111, 1148].as_slice())
            );
            assert_eq!(
                filtered.filtering_stats(),
                vec![("missing_data", pair(500, 423))]
            );
        }

        /// The same filter under the reader counts over the 50 individuals
        /// of the file and keeps 26 variants: where the step sits among the
        /// steps is what a user's numbers turn on.
        #[test]
        fn the_missing_data_filter_under_the_three_keeps_the_26_variants() {
            let source = many_vcf_reader(Some(7), Needs::GTS | Needs::CHROM_POS);
            let filter = VarFilter::new(MaxMissingRate(0.0)).expect("the filter");
            let filtered = FilteredReader::new(source, filter).expect("the filter over the file");
            let mut reader =
                IndividualsReader::new(filtered, &the_three_names()).expect("the reader");

            let blocks = blocks_of(&mut reader).expect("the blocks");

            let positions = positions_of_blocks(&blocks);
            assert_eq!(positions.len(), 26);
            assert_eq!(
                positions.get(..5),
                Some([1259, 2110, 2480, 3072, 3257].as_slice())
            );
            // The counts of the filter under it come up through the reader,
            // which has none of its own.
            assert_eq!(
                reader.filtering_stats(),
                vec![("missing_data", pair(500, 26))]
            );
        }

        /// The blocks are the size of the source's, worked out from the
        /// individuals of the source: the reader keeps every variant and
        /// does not ask its source for bigger blocks now that the rows are
        /// shorter.
        #[test]
        fn the_blocks_are_the_size_of_the_blocks_of_the_source() {
            let mut reader = of_the_three_names(Some(7), Needs::GTS);

            let blocks = blocks_of(&mut reader).expect("the blocks");

            let sizes: Vec<usize> = blocks.iter().map(|block| block.num_vars).collect();
            // 500 variants in blocks of 7: 71 blocks of 7 and one of 3.
            assert_eq!(sizes.len(), 72);
            assert_eq!(sizes.first(), Some(&7));
            assert_eq!(sizes.last(), Some(&3));
        }

        /// The reader always needs the genotypes: it asks its source for
        /// them with whatever its consumer asked for, so the blocks it gives
        /// hold them although the consumer wanted the positions alone.
        #[test]
        fn the_reader_asks_its_source_for_the_genotypes_with_the_fields_of_its_consumer() {
            let source = many_vcf_reader(Some(7), Needs::CHROM_POS);
            let mut reader =
                IndividualsReader::new(source, &the_three_names()).expect("the reader");
            reader.set_needs(Needs::CHROM_POS);

            let blocks = blocks_of(&mut reader).expect("the blocks");

            let first = blocks.first().expect("a block");
            assert_eq!(first.num_individuals, 3);
            assert!(!first.gts.is_empty());
            assert_eq!(first.fields(), Needs::GTS | Needs::CHROM_POS);
        }

        /// After an error the reader gives `None` at every call and does not
        /// ask its source again, which is the rule of a reader of
        /// `docs/specs/block.md`.
        #[test]
        fn after_an_error_it_gives_no_block_and_does_not_ask_its_source_again() {
            let source = GivenBlocks::failing_at(
                vec![
                    block_of_the_worked_example(&[0, 1]),
                    block_of_the_worked_example(&[2, 3]),
                ],
                2,
            );
            let calls = source.calls();
            let mut reader = IndividualsReader::new(source, &["ind5".to_owned()])
                .expect("the reader of one individual");

            assert!(reader.next_block().expect("the first block").is_some());
            assert!(reader.next_block().is_err());
            assert!(
                reader
                    .next_block()
                    .expect("no block after the error")
                    .is_none()
            );
            assert_eq!(calls.load(std::sync::atomic::Ordering::SeqCst), 2);
        }

        /// A source that gives a block of no variants has a defect, and the
        /// reader gives the error of it and does not ask the source again.
        #[test]
        fn a_source_that_gives_a_block_of_no_variants_is_the_error_of_a_defect() {
            let source = GivenBlocks::of(vec![block_of_the_worked_example(&[])]);
            let mut reader = IndividualsReader::new(source, &["ind1".to_owned()])
                .expect("the reader of one individual");

            let error = reader.next_block().expect_err("the block of no variants");

            assert!(
                matches!(error, Error::ReaderGaveABlockOfNoVariants),
                "the error is {error}"
            );
            assert!(
                reader
                    .next_block()
                    .expect("no block after the error")
                    .is_none()
            );
        }

        /// The chain of a pass builds the reader from the step of the
        /// filter of individuals, and a threshold filter after it in the
        /// steps counts over the kept individuals: the 423 variants again,
        /// through `chain_of` alone.
        #[test]
        fn chain_of_builds_the_reader_from_the_step_of_the_filter() {
            let source = Box::new(many_vcf_reader(Some(7), Needs::GTS | Needs::CHROM_POS));
            let steps = vec![
                PassStep::KeepIndividuals(the_three_names()),
                PassStep::VarFilter(MaxMissingRate(0.0)),
            ];

            let mut chain = chain_of(source, &steps).expect("the chain of the two steps");

            assert_eq!(chain.individuals(), the_three_names().as_slice());
            let blocks = blocks_of(&mut chain).expect("the blocks");
            assert_eq!(positions_of_blocks(&blocks).len(), 423);
            assert_eq!(
                chain.filtering_stats(),
                vec![("missing_data", pair(500, 423))]
            );
        }

        /// A second filter of individuals is refused, by the chain and by
        /// the function a binding crate calls when a user adds the step:
        /// two lists keep the individuals that are in both, which is one
        /// list. The error names the kind and no threshold, which a list of
        /// individuals has none of.
        #[test]
        fn a_second_filter_of_individuals_is_refused_with_its_kind() {
            let steps = vec![
                PassStep::KeepIndividuals(the_three_names()),
                PassStep::VarFilter(MaxMaf(0.8)),
                PassStep::KeepIndividuals(vec!["ind05".to_owned()]),
            ];
            let source = Box::new(many_vcf_reader(Some(7), Needs::GTS));

            let of_the_chain = chain_of(source, &steps)
                .err()
                .expect("the chain was refused");
            let of_the_refusal = refuse_a_second_filter_of_a_kind(
                steps.get(..2).unwrap_or_default(),
                steps.get(2).expect("the step"),
            )
            .expect_err("the second filter of individuals was refused");

            for error in [of_the_chain, of_the_refusal] {
                let message = error.to_string();
                assert!(
                    matches!(
                        error,
                        Error::FilterOfIndividualsThatIsSet {
                            kind: "individuals"
                        }
                    ),
                    "{message}"
                );
                assert!(message.contains("individuals"), "{message}");
            }

            // One filter of individuals among the steps is taken, and the
            // threshold filters beside it are of other kinds.
            assert!(
                refuse_a_second_filter_of_a_kind(
                    &steps_of(&[MaxMaf(0.8)]),
                    &PassStep::KeepIndividuals(the_three_names())
                )
                .is_ok()
            );
        }

        /// What the chain refuses of the names, which is what
        /// `resolve_individuals` refuses: it is read against the individuals
        /// of the reader the step is put over, and no block was read when it
        /// comes.
        #[test]
        fn chain_of_refuses_a_name_that_is_not_an_individual_of_the_source() {
            let source = Box::new(many_vcf_reader(Some(7), Needs::GTS));
            let steps = vec![PassStep::KeepIndividuals(vec![
                "ind05".to_owned(),
                "nope".to_owned(),
            ])];

            let error = chain_of(source, &steps)
                .err()
                .expect("the chain was refused");

            let message = error.to_string();
            assert!(
                matches!(error, Error::IndividualNotInTheSource { .. }),
                "{message}"
            );
            assert!(message.contains("nope"), "{message}");
        }
    }

    /// The tests of [`PassStep`] and of the chain of readers built from a
    /// list of them. The module is named after the type, and not
    /// `pass_steps`, so that `cargo test -- PassStep` runs them.
    #[expect(
        non_snake_case,
        reason = "the module is named after the type it tests, PassStep, so that the \
                  tests of the steps of a pass are the ones cargo test -- PassStep runs"
    )]
    mod PassSteps {
        use super::{
            GivenBlocks, MaxMaf, MaxMissingRate, MaxObsHet, PassStep, block_of_the_worked_example,
            blocks_of, chain_of, pair, positions_of_blocks, steps_of,
        };

        /// The kind of each step is the name a Python and a TypeScript user
        /// reads for it: the three of the threshold filters, which are the
        /// keys their counts have, and `individuals` for the filter of
        /// individuals.
        #[test]
        fn the_kind_of_each_step_is_the_name_the_user_reads() {
            assert_eq!(
                PassStep::VarFilter(MaxMissingRate(0.04)).kind(),
                "missing_data"
            );
            assert_eq!(PassStep::VarFilter(MaxMaf(0.8)).kind(), "maf");
            assert_eq!(PassStep::VarFilter(MaxObsHet(0.5)).kind(), "obs_het");
            assert_eq!(
                PassStep::KeepIndividuals(vec!["ind05".to_owned()]).kind(),
                "individuals"
            );
        }

        /// The chain built from a list of steps is the chain of the filters
        /// of those steps, in their order: the three threshold filters at
        /// the thresholds of the worked example of "How it is verified" of
        /// `docs/specs/filters.md`, 0.4, 0.88 and 0.25, leave variant 5 of
        /// its six, which is what the filters put one over another by hand
        /// leave, with the counts of each filter.
        #[test]
        fn the_chain_of_the_three_threshold_filters_keeps_variant_5_of_the_worked_example() {
            let source = GivenBlocks::of(vec![block_of_the_worked_example(&[0, 1, 2, 3, 4, 5])]);
            let steps = steps_of(&[MaxMissingRate(0.4), MaxMaf(0.88), MaxObsHet(0.25)]);

            let mut chain = chain_of(Box::new(source), &steps).expect("the chain of the steps");

            let blocks = blocks_of(&mut chain).expect("the blocks");
            // The missing data filter keeps the variants 1, 2, 3 and 5 of
            // the six; the maf filter keeps 2, 3 and 5 of those four, since
            // the 8/9 of variant 1 is above 0.88; and the observed
            // heterozygosity filter keeps 5 of those three, since the 1/3 of
            // variant 2 and the 4/4 of variant 3 are above 0.25.
            assert_eq!(positions_of_blocks(&blocks), [5]);
            assert_eq!(
                chain.filtering_stats(),
                vec![
                    ("obs_het", pair(3, 1)),
                    ("maf", pair(4, 3)),
                    ("missing_data", pair(6, 4)),
                ]
            );
        }
    }
}
