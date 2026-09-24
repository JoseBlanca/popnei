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
//! What is built here so far is the pass and those counts of variants. The
//! five statistics are added on top of them.

use std::collections::HashSet;

use crate::block::{Block, BlockReader, alleles_of_a_chunk, alleles_per_var_of};
use crate::error::{Error, Result};
use crate::stats::{every_individual_in_order, min_called_alleles};
use crate::variant::{AlleleCounts, Needs, count_alleles, count_alleles_of};

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
}

impl OfAPop {
    /// The counts of one population before any variant is read.
    fn none() -> OfAPop {
        OfAPop { num_vars: 0 }
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

    /// The variants that counted for every population. It is the divisor of
    /// the private alleles.
    #[must_use]
    pub fn num_vars_every_pop(&self) -> u64 {
        self.num_vars_every_pop
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

/// What every row of a pass is read with: the populations and the rule for
/// which variants count for them.
#[derive(Debug)]
struct OfThePass<'a> {
    pops: &'a [PopOfThePass],
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
/// [`DiversityStats::FOLDED_SFS`] asked for with no `num_called_alleles`, a
/// `num_called_alleles` below 2, a population with no individual, an index
/// that is not an individual of the dataset, an individual asked for more
/// than once, no variant in the reader, a variant of more alleles than a
/// count of them holds, and those of the reader.
pub fn calc_pop_diversity<R: BlockReader + ?Sized>(
    reader: &mut R,
    pops: &[&[usize]],
    options: &DiversityOptions,
) -> Result<PopDiversity> {
    check_the_draw(options)?;
    let num_individuals = reader.individuals().len();
    let of_the_pops = pops_of_the_pass(pops, num_individuals)?;
    let ploidy = reader.ploidy();
    // The five statistics follow from the genotypes of a row, so no column
    // of a block is read and the reader is asked to fill none of them.
    reader.set_needs(Needs::GTS);
    let of_the_pass = OfThePass {
        pops: &of_the_pops,
        min_called_alleles: min_called_alleles(
            options.min_num_individuals,
            // Every reader of popnei gives a ploidy of 255 at most, and a
            // threshold built from a larger one would be one no population
            // ever meets, which is what such a ploidy asks for.
            u32::try_from(ploidy).unwrap_or(u32::MAX),
        ),
        ploidy,
    };
    let mut totals = Totals::of(of_the_pops.len());
    let mut num_vars: u64 = 0;
    while let Some(block) = reader.next_block()? {
        let alleles_per_var = alleles_per_var_of(&block, num_individuals, ploidy)?;
        add_the_block(&block, alleles_per_var, &of_the_pass, &mut totals)?;
        // A `usize` is 64 bits on the targets popnei builds natively for and
        // 32 in wasm, so every one of them is a `u64`; and a pass of more
        // than 18446744073709551615 variants reads more rows than any
        // source holds.
        num_vars = num_vars.saturating_add(u64::try_from(block.num_vars).unwrap_or(u64::MAX));
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
    })
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
    // One array of counts for every row and every population, which
    // `count_alleles_of` clears before it counts, so the loop over the rows
    // allocates nothing for a variant. The statistics read it while the
    // population is in hand; the private alleles, which need every
    // population's counts at one variant at once, will keep them in an
    // array of the populations by the alleles of the row instead.
    let mut counts: AlleleCounts = [0; 128];
    for row in gts.chunks_exact(alleles_per_var) {
        let mut every_pop = true;
        for (of_the_pop, counted) in of_the_pass.pops.iter().zip(totals.pops.iter_mut()) {
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
            let called_alleles = if of_the_whole_row {
                count_alleles(row, &mut counts)?
            } else {
                count_alleles_of(
                    row,
                    of_the_pass.ploidy,
                    &of_the_pop.individuals,
                    &mut counts,
                )?
                .called_alleles
            };
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
        }
        if every_pop {
            totals.num_vars_every_pop = totals.num_vars_every_pop.saturating_add(1);
        }
    }
    Ok(())
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

    /// The six variants of the worked example, in blocks of
    /// `num_vars_per_block` variants of the five diploid individuals.
    pub(super) fn the_worked_example(num_vars_per_block: usize) -> GivenBlocks {
        let rows: Vec<&[i8]> = THE_SIX_VARIANTS.iter().map(|row| &row[..]).collect();
        GivenBlocks::of(blocks_of(&rows, 5, 2, num_vars_per_block))
    }
}

#[cfg(test)]
mod the_pass {
    use super::fixtures::{GivenBlocks, POP1, POP2, blocks_of, the_worked_example};
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
        }
    }

    /// No population at all is one population of every individual of the
    /// reader, which counts the four variants something was called at.
    #[test]
    fn no_population_is_one_population_of_every_individual() {
        let mut reader = the_worked_example(6);

        let diversity = calc_pop_diversity(&mut reader, &[], &options_with_no_draw(1))
            .expect("the diversity of every individual");

        assert_eq!(diversity.num_pops(), 1);
        assert_eq!(diversity.num_vars(0), Some(4));
        assert_eq!(diversity.num_vars_every_pop(), 4);
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
}
