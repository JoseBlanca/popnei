//! The window of blocks that the fall-off of r² with distance is counted
//! over, from "How it runs" of the item "LD against distance, per
//! population" of `docs/specs/ld.md`.
//!
//! Two variants make a pair only when they are on one chromosome and no
//! further apart than `max_dist`, so a pass that counts every pair never
//! needs the whole dataset in memory: it needs the blocks whose variants
//! are still within `max_dist` of the newest variant it has read and on
//! that variant's chromosome. [`TheWindowOfTheBlocks`] is that set of
//! blocks. It takes the blocks of the reader one after another and drops
//! one as soon as every variant of it is further back than `max_dist` or
//! on another chromosome.
//!
//! A block is dropped whole, so the window holds every variant that is
//! within reach and, beside them, the variants of the same blocks that
//! have already fallen behind. How many of those there are depends on the
//! size of the blocks the reader gives; which pairs the pass counts does
//! not, because the pass takes the distance of each pair from the
//! positions and not from the window.
use std::collections::VecDeque;

use crate::block::{Block, BlockReader};
use crate::error::{Error, Result};
use crate::variant::Needs;

use super::{LdDosages, a_vector_of, r2_between, the_copy_of, the_memory_for};

/// Where one variant lies: the number of its chromosome in the
/// [`ChromTable`](crate::variant::ChromTable) of the reader the block came
/// from, and its position on it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TheVariantOfTheWindow {
    /// The number of the chromosome, in the table of the reader.
    chrom: u32,
    /// The position of the variant, 1 based as in a VCF.
    pos: u64,
}

/// The blocks whose variants are within `max_dist` of the newest variant
/// read and on its chromosome.
///
/// The blocks are taken in the order the reader gives them, which is the
/// order of the variants of the source, and the newest variant read is the
/// last variant of the last block taken. A block is dropped when every
/// variant of it is further back than `max_dist` or on another chromosome,
/// and the blocks are dropped from the oldest end: a block none of whose
/// variants are in reach while an older one still has one is kept, which
/// costs memory and no pair. A block of no variant holds nothing of the
/// window and leaves the newest variant where it was, and falls out as
/// soon as it is the oldest block held.
///
/// Every block held has the chromosome and the position of each of its
/// variants, which [`TheWindowOfTheBlocks::take_the_block`] refuses a
/// block without, so [`TheWindowOfTheBlocks::variants`] gives one variant
/// for each of the variants the blocks count.
struct TheWindowOfTheBlocks {
    /// How far back of the newest variant read a variant is still held.
    max_dist: u64,
    /// The blocks held, the oldest first.
    held: VecDeque<Block>,
    /// The last variant of the last block that had one, and `None` before
    /// a block with a variant has been taken.
    newest: Option<TheVariantOfTheWindow>,
}

impl TheWindowOfTheBlocks {
    /// A window that holds no block and reaches `max_dist` back from the
    /// newest variant read.
    fn of(max_dist: u64) -> TheWindowOfTheBlocks {
        TheWindowOfTheBlocks {
            max_dist,
            held: VecDeque::new(),
            newest: None,
        }
    }

    /// Takes the next block of the reader and drops the blocks that fell
    /// out of the window with it, which are the oldest ones every variant
    /// of which is now further back than `max_dist` or on another
    /// chromosome.
    ///
    /// It gives how many blocks were dropped, which are the oldest ones
    /// held, so that what a caller keeps beside each held block, the
    /// dosages of a population over it, is dropped with it.
    ///
    /// # Errors
    ///
    /// What [`Block::check`] refuses, a block whose columns are not of its
    /// number of variants among them, and
    /// [`Error::FieldsNotInTheBlock`] for a block that has not the
    /// chromosome and the position of its variants.
    fn take_the_block(&mut self, block: Block) -> Result<usize> {
        block.check()?;
        let missing = Needs::CHROM_POS.difference(block.fields());
        if !missing.is_empty() {
            return Err(Error::FieldsNotInTheBlock { fields: missing });
        }
        if let Some(last) = the_variants_of(&block).next_back() {
            self.newest = Some(last);
        }
        self.held.push_back(block);
        Ok(self.drop_what_fell_out())
    }

    /// How many blocks the window holds.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the pass over a window reads the variants of each population and \
                      the dosages over them, and not the window itself; what this \
                      answers is read by the tests of the window alone"
        )
    )]
    fn num_blocks(&self) -> usize {
        self.held.len()
    }

    /// How many variants the blocks it holds have, which is how many
    /// [`TheWindowOfTheBlocks::variants`] gives.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the pass over a window reads the variants of each population and \
                      the dosages over them, and not the window itself; what this \
                      answers is read by the tests of the window alone"
        )
    )]
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "each block held has one u32 of chromosome and one u64 of position for \
                  each of its variants, so the variants of the window are at most the \
                  bytes of this machine divided by twelve and their total is a usize"
    )]
    fn num_vars(&self) -> usize {
        self.held
            .iter()
            .fold(0, |num_vars, block| num_vars + block.num_vars)
    }

    /// The blocks it holds, the oldest first.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the pass over a window reads the variants of each population and \
                      the dosages over them, and not the window itself; what this \
                      answers is read by the tests of the window alone"
        )
    )]
    fn blocks(&self) -> impl ExactSizeIterator<Item = &Block> {
        self.held.iter()
    }

    /// Where each variant of the blocks it holds lies, the oldest block
    /// first and inside a block in the order of its variants.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the pass over a window reads the variants of each population and \
                      the dosages over them, and not the window itself; what this \
                      answers is read by the tests of the window alone"
        )
    )]
    fn variants(&self) -> impl Iterator<Item = TheVariantOfTheWindow> {
        self.held.iter().flat_map(the_variants_of)
    }

    /// The newest variant read, which the window reaches back from, and
    /// `None` before a block with a variant has been taken.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the pass over a window reads the variants of each population and \
                      the dosages over them, and not the window itself; what this \
                      answers is read by the tests of the window alone"
        )
    )]
    fn the_newest_variant(&self) -> Option<TheVariantOfTheWindow> {
        self.newest
    }

    /// Whether the variant is on the chromosome of the newest variant read
    /// and no further than `max_dist` from it.
    fn in_the_window(&self, variant: TheVariantOfTheWindow) -> bool {
        match self.newest {
            None => false,
            Some(newest) => {
                variant.chrom == newest.chrom && newest.pos.abs_diff(variant.pos) <= self.max_dist
            }
        }
    }

    /// Drops the oldest blocks no variant of which is in the window, and
    /// gives how many were dropped.
    fn drop_what_fell_out(&mut self) -> usize {
        let held_before = self.held.len();
        while let Some(oldest) = self.held.front() {
            // The variants of a block are in the order of the source, so
            // the last one is the one most likely to be in reach and the
            // walk back from it stops at the first variant it finds there.
            let any_in_the_window = the_variants_of(oldest)
                .rev()
                .any(|variant| self.in_the_window(variant));
            if any_in_the_window {
                break;
            }
            self.held.pop_front();
        }
        // The blocks only leave, so what is held now is at most what was
        // held before.
        held_before.saturating_sub(self.held.len())
    }
}

/// Where each variant of the block lies, in the order of its variants.
///
/// A block without the chromosome or the position of its variants gives
/// none, which no block of a window is:
/// [`TheWindowOfTheBlocks::take_the_block`] refuses one.
fn the_variants_of(
    block: &Block,
) -> impl DoubleEndedIterator<Item = TheVariantOfTheWindow> + ExactSizeIterator {
    let chroms = block.chrom.as_deref().unwrap_or(&[]);
    let poss = block.pos.as_deref().unwrap_or(&[]);
    chroms
        .iter()
        .zip(poss)
        .map(|(chrom, pos)| TheVariantOfTheWindow {
            chrom: *chrom,
            pos: *pos,
        })
}

/// How many individuals the blocks of one pass hold and how many alleles
/// the genotype of one of them has, which the first block taken fixes.
#[derive(Debug, Clone, Copy)]
struct TheSourceOfThePass {
    /// How many individuals each block of the pass has.
    num_individuals: usize,
    /// How many alleles the genotype of one individual holds.
    ploidy: usize,
    /// The individuals times the ploidy, which is how many alleles one
    /// variant has and how long a row of the genotypes is.
    alleles_per_var: usize,
}

/// The window of blocks and, for each population, the dosages of the
/// variants of those blocks that passed its major allele frequency.
///
/// A population is a set of individuals, and "LD against distance, per
/// population" of `docs/specs/ld.md` works the major allele frequency of a
/// variant out over the individuals of one population and not over all of
/// them, so two populations of one pass keep different variants: a variant
/// whose frequency is above `max_allowed_maf` in one of them is left out of
/// that one and stays in the other. A variant with no called allele among
/// the individuals of a population has no frequency there and is left out
/// of it whatever the threshold, and it can still be in another population
/// that called it.
///
/// Each population holds one set of dosages over every variant of the
/// window it kept, which is what the tiles of the pairs take their rows
/// from with [`LdDosages::rows`]: the tiles are cut at fixed multiples
/// counted from the first variant of the pass and not where the blocks
/// end, so they do not fall inside one block. The dosages are built again
/// each time a block is taken, from the genotypes of the variants held,
/// and the dosages of a variant do not depend on which other variants are
/// in the set: the major allele of a variant is that of the individuals of
/// the population at that variant alone, so a variant has the same row and
/// the same frequency whatever the size of the blocks.
struct TheDosagesOfThePops {
    /// The blocks whose variants are still in reach of the newest one
    /// read, which the variants held are those of.
    window: TheWindowOfTheBlocks,
    /// The largest major allele frequency a variant has in a population
    /// and is still counted there, both ends included.
    max_allowed_maf: f64,
    /// What the first block taken said the source holds, and `None` before
    /// a block has been taken.
    of_the_source: Option<TheSourceOfThePass>,
    /// How many blocks fell out of the window when the newest block was
    /// taken, which the populations drop when the next one is taken.
    fell_out: usize,
    /// One for each population, in the order they were given.
    pops: Vec<ThePopOverTheWindow>,
}

impl TheDosagesOfThePops {
    /// The dosages of each population of `pops` over a window that reaches
    /// `max_dist` back, keeping the variants whose major allele frequency
    /// in that population is at most `max_allowed_maf`.
    ///
    /// A population is the indices of its individuals among the
    /// individuals of a block, in the order they are given, and `pops`
    /// that holds no population at all is one population of every
    /// individual of the source.
    ///
    /// # Errors
    ///
    /// [`Error::LdMaxAllowedMafOutOfRange`] for a `max_allowed_maf` that
    /// is NaN or is not from 0 to 1, [`Error::LdPopWithNoIndividual`] for
    /// a population that names no individual, and
    /// [`Error::LdNoMemory`] when this machine does not give the memory of
    /// the populations.
    fn of(pops: &[&[usize]], max_dist: u64, max_allowed_maf: f64) -> Result<TheDosagesOfThePops> {
        // NaN is in no range, so this refuses it with the numbers outside
        // 0 to 1, which is what the error says.
        if !(0.0..=1.0).contains(&max_allowed_maf) {
            return Err(Error::LdMaxAllowedMafOutOfRange {
                value: max_allowed_maf,
            });
        }
        let mut of_them: Vec<ThePopOverTheWindow> = Vec::new();
        let no_memory = || Error::LdNoMemory {
            what: "the populations of the fall-off of r² with distance",
            values: pops.len(),
        };
        match pops.is_empty() {
            // The one population of every individual, which is the empty
            // set of individuals that `LdDosages::of_block` reads as all
            // of them.
            true => {
                of_them.try_reserve_exact(1).map_err(|_| no_memory())?;
                of_them.push(ThePopOverTheWindow::of(&[])?);
            }
            false => {
                of_them
                    .try_reserve_exact(pops.len())
                    .map_err(|_| no_memory())?;
                for (pop, individuals) in pops.iter().enumerate() {
                    if individuals.is_empty() {
                        return Err(Error::LdPopWithNoIndividual { pop });
                    }
                    of_them.push(ThePopOverTheWindow::of(individuals)?);
                }
            }
        }
        Ok(TheDosagesOfThePops {
            window: TheWindowOfTheBlocks::of(max_dist),
            max_allowed_maf,
            of_the_source: None,
            fell_out: 0,
            pops: of_them,
        })
    }

    /// Takes the next block of the reader: each population keeps the
    /// variants of it that passed its major allele frequency, the blocks
    /// that fell out of the window when the block before it was taken take
    /// their variants with them, and the dosages of each population are
    /// built again over what is left.
    ///
    /// The blocks that fall out with this block are dropped when the next
    /// one is taken, and not here, because the pairs of a block are
    /// counted against the window as it stands when that block arrives.
    /// The window reaches `max_dist` back from the newest variant read,
    /// which is the last variant of this block, and a variant within
    /// `max_dist` of the first variant of the block can be further than
    /// that from its last: dropping it here would leave a pair of the
    /// block uncounted, and how many such pairs there are would depend on
    /// the size of the blocks the reader gives. What the deferral costs is
    /// the memory of one window more, which the blocks that fell out are
    /// held for one step of the pass.
    ///
    /// # Errors
    ///
    /// [`Error::FieldsNotInTheBlock`] for a block that has not the
    /// genotypes or the chromosome and the position of its variants,
    /// [`Error::BlockWithNoGenotypeOfAVariant`] for a first block whose
    /// individuals times its ploidy are 0,
    /// [`Error::BlocksDoNotFitTogether`] for a block of other individuals
    /// or another ploidy than the first one taken, what [`Block::check`]
    /// and [`TheWindowOfTheBlocks::take_the_block`] refuse, and what
    /// [`LdDosages::of_block`] refuses, a window this machine has not the
    /// memory of the dosages of among them.
    fn take_the_block(&mut self, block: Block) -> Result<()> {
        block.check()?;
        let missing = (Needs::GTS | Needs::CHROM_POS).difference(block.fields());
        if !missing.is_empty() {
            return Err(Error::FieldsNotInTheBlock { fields: missing });
        }
        let of_the_source = self.the_source_of(&block)?;
        // What fell out of the window when the block before this one was
        // taken, whose pairs have been counted since.
        let dropped = std::mem::take(&mut self.fell_out);
        for pop in &mut self.pops {
            pop.the_oldest_blocks_are_dropped(dropped, of_the_source);
        }
        // The variants of the block are taken before the window is told of
        // it, so that what each population holds is in the order of the
        // blocks of the window and the blocks that fall out are the oldest
        // of both.
        for pop in &mut self.pops {
            pop.take_the_block(&block, of_the_source, self.max_allowed_maf)?;
        }
        self.fell_out = self.window.take_the_block(block)?;
        for pop in &mut self.pops {
            pop.the_dosages_are_built(of_the_source)?;
        }
        Ok(())
    }

    /// What the source holds, from the first block taken, and the block
    /// checked against it.
    ///
    /// # Errors
    ///
    /// [`Error::BlockWithNoGenotypeOfAVariant`] when the first block taken
    /// holds the genotypes of no individual,
    /// [`Error::BlocksDoNotFitTogether`] when this block has other
    /// individuals or another ploidy than the first one, and what
    /// [`Block::alleles_per_var`] refuses.
    fn the_source_of(&mut self, block: &Block) -> Result<TheSourceOfThePass> {
        if let Some(of_the_source) = self.of_the_source {
            if block.num_individuals != of_the_source.num_individuals
                || block.ploidy != of_the_source.ploidy
            {
                return Err(Error::BlocksDoNotFitTogether {
                    num_individuals: of_the_source.num_individuals,
                    ploidy: of_the_source.ploidy,
                    found_num_individuals: block.num_individuals,
                    found_ploidy: block.ploidy,
                });
            }
            return Ok(of_the_source);
        }
        let alleles_per_var = block.alleles_per_var()?;
        if alleles_per_var == 0 {
            return Err(Error::BlockWithNoGenotypeOfAVariant {
                num_individuals: block.num_individuals,
                ploidy: block.ploidy,
            });
        }
        let of_the_source = TheSourceOfThePass {
            num_individuals: block.num_individuals,
            ploidy: block.ploidy,
            alleles_per_var,
        };
        self.of_the_source = Some(of_the_source);
        Ok(of_the_source)
    }

    /// How many populations the pass counts.
    fn num_pops(&self) -> usize {
        self.pops.len()
    }

    /// The population at that position among the ones given, and `None`
    /// when it is not one of them.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the pass over a window reads the variants of each population and \
                      the dosages over them, and not the window itself; what this \
                      answers is read by the tests of the window alone"
        )
    )]
    fn pop(&self, pop: usize) -> Option<&ThePopOverTheWindow> {
        self.pops.get(pop)
    }

    /// Every population of the pass, in the order they were given.
    fn the_pops(&self) -> &[ThePopOverTheWindow] {
        &self.pops
    }
}

/// One population over the window: the variants of the held blocks that
/// passed its major allele frequency, and the dosages of them over its own
/// individuals.
struct ThePopOverTheWindow {
    /// The index of each individual of the population among the
    /// individuals of a block, in the order they were given, and empty for
    /// the one population of every individual.
    individuals: Vec<usize>,
    /// The genotypes of the variants held, variant after variant, each of
    /// them the row the block it came from had, over every individual of
    /// that block. The dosages of one step of the pass are built from it
    /// and it is kept with the memory it had, so a pass allocates it once
    /// and grows it to the largest window it meets.
    gts: Vec<i8>,
    /// Where each variant held lies, the oldest first, one for each
    /// variant of `gts`.
    variants: Vec<TheVariantOfTheWindow>,
    /// How many variants it kept of each block of the window, the oldest
    /// block first, so that the blocks that fall out take their variants
    /// out of `gts` and of `variants`.
    kept_of_each_block: VecDeque<usize>,
    /// How many variants it has kept since the pass began, which is the
    /// count "What it gives" of `docs/specs/ld.md` gives the user for the
    /// population.
    num_vars: u64,
    /// The number, among the variants it has kept since the pass began, of
    /// the first variant it holds, which is how many have fallen out of
    /// the window. The tiles of the pairs are cut at fixed multiples
    /// counted from the first variant of the pass, and this is what says
    /// where the variants held start.
    first_var: u64,
    /// The dosages of the variants held, and `None` before the first block
    /// of the pass has been taken.
    dosages: Option<LdDosages>,
}

impl ThePopOverTheWindow {
    /// A population of those individuals that holds no variant yet. An
    /// empty `individuals` is every individual of the block.
    ///
    /// # Errors
    ///
    /// [`Error::LdNoMemory`] when this machine does not give the memory of
    /// the individuals.
    fn of(individuals: &[usize]) -> Result<ThePopOverTheWindow> {
        Ok(ThePopOverTheWindow {
            individuals: the_copy_of(
                individuals,
                &the_memory_for("the individuals of a population", individuals.len()),
            )?,
            gts: Vec::new(),
            variants: Vec::new(),
            kept_of_each_block: VecDeque::new(),
            num_vars: 0,
            first_var: 0,
            dosages: None,
        })
    }

    /// Keeps the variants of the block whose major allele frequency among
    /// these individuals is at most `max_allowed_maf`, and counts how many
    /// of them the block gave.
    ///
    /// A variant with no called allele among these individuals has no
    /// frequency and is left out, as "The cases" of `docs/specs/ld.md`
    /// has it.
    ///
    /// # Errors
    ///
    /// What [`LdDosages::of_block`] refuses, an index that is not an
    /// individual of the block and an individual asked for twice among
    /// them, and [`Error::LdNoMemory`] when this machine does not give the
    /// memory of the genotypes of the variants kept.
    fn take_the_block(
        &mut self,
        block: &Block,
        of_the_source: TheSourceOfThePass,
        max_allowed_maf: f64,
    ) -> Result<()> {
        // The frequency of each variant of the block among the individuals
        // of this population, which `LdDosages::of_block` works out over
        // the individuals it is built with.
        let of_the_block = LdDosages::of_block(block, &self.individuals)?;
        let mut kept = 0_usize;
        let rows = block
            .gts
            .chunks_exact(of_the_source.alleles_per_var)
            .zip(the_variants_of(block))
            .enumerate();
        for (var, (gts, variant)) in rows {
            if !of_the_block
                .maf(var)
                .is_some_and(|maf| maf <= max_allowed_maf)
            {
                continue;
            }
            self.gts
                .try_reserve(gts.len())
                .map_err(|_| Error::LdNoMemory {
                    what: "the genotypes of the variants of the window",
                    values: gts.len(),
                })?;
            self.gts.extend_from_slice(gts);
            self.variants
                .try_reserve(1)
                .map_err(|_| Error::LdNoMemory {
                    what: "where each variant of the window lies",
                    values: self.variants.len(),
                })?;
            self.variants.push(variant);
            // At most the variants of the block, which this machine
            // counted when it read them.
            kept = kept.saturating_add(1);
        }
        self.kept_of_each_block
            .try_reserve(1)
            .map_err(|_| Error::LdNoMemory {
                what: "the variants each block of the window gave",
                values: self.kept_of_each_block.len(),
            })?;
        self.kept_of_each_block.push_back(kept);
        // Every variant counted here was read from a source, and a u64
        // counts 1.8e19 of them: a pass that passed this number would have
        // read more bytes than any storage holds.
        self.num_vars = self.num_vars.saturating_add(the_count_of(kept));
        Ok(())
    }

    /// Takes the variants of the `blocks` oldest blocks of the window out
    /// of what it holds, which are the ones that fell out of it.
    fn the_oldest_blocks_are_dropped(&mut self, blocks: usize, of_the_source: TheSourceOfThePass) {
        let mut dropped_vars = 0_usize;
        for _ in 0..blocks {
            match self.kept_of_each_block.pop_front() {
                // The variants of the blocks held, which are at most the
                // variants of one window, a number this machine counted.
                Some(kept) => dropped_vars = dropped_vars.saturating_add(kept),
                None => break,
            }
        }
        if dropped_vars == 0 {
            return;
        }
        // The genotypes held are the variants held times the alleles of one
        // variant, which is the length of a vector this machine gave, so
        // the alleles of the variants dropped are at most that number.
        let alleles = dropped_vars.saturating_mul(of_the_source.alleles_per_var);
        self.gts.drain(..alleles.min(self.gts.len()));
        self.variants.drain(..dropped_vars.min(self.variants.len()));
        self.first_var = self.first_var.saturating_add(the_count_of(dropped_vars));
    }

    /// Builds the dosages of the variants it holds, over its individuals.
    ///
    /// # Errors
    ///
    /// What [`LdDosages::of_block`] refuses, a window this machine has not
    /// the memory of the three matrices of among them.
    fn the_dosages_are_built(&mut self, of_the_source: TheSourceOfThePass) -> Result<()> {
        let mut block = Block {
            num_vars: self.variants.len(),
            num_individuals: of_the_source.num_individuals,
            ploidy: of_the_source.ploidy,
            gts: std::mem::take(&mut self.gts),
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        };
        let built = LdDosages::of_block(&block, &self.individuals);
        // The genotypes are taken back with the memory they have, whether
        // the dosages were built or not, so that the next step of the pass
        // is gathered into the same buffer.
        self.gts = std::mem::take(&mut block.gts);
        self.dosages = Some(built?);
        Ok(())
    }

    /// The dosages of the variants it holds, in the order of the source,
    /// and `None` before the first block of the pass has been taken.
    fn dosages(&self) -> Option<&LdDosages> {
        self.dosages.as_ref()
    }

    /// Where each variant it holds lies, in the same order as the rows of
    /// its dosages.
    fn variants(&self) -> &[TheVariantOfTheWindow] {
        &self.variants
    }

    /// How many variants it has kept since the pass began.
    fn num_vars(&self) -> u64 {
        self.num_vars
    }

    /// The number, among the variants it has kept since the pass began, of
    /// the first variant it holds.
    fn first_var(&self) -> u64 {
        self.first_var
    }

    /// How many variants it kept of the block taken last, which are the
    /// newest variants it holds, and 0 before a block has been taken.
    fn kept_of_the_newest_block(&self) -> usize {
        self.kept_of_each_block.back().copied().unwrap_or(0)
    }
}

/// A count of variants as a result carries it, which is a `u64` and not a
/// `usize` so that it is the same number in WebAssembly, where a `usize` is
/// 32 bits.
fn the_count_of(num_vars: usize) -> u64 {
    u64::try_from(num_vars).unwrap_or(u64::MAX)
}

/// How many variants one tile of the pairs of the window holds.
///
/// The pairs of a population are worked out tile pair by tile pair, as the
/// matrix of every pair is, and the tiles are cut at fixed multiples of
/// this number counted from the first variant that population kept in the
/// pass, so they do not move when the blocks of the reader do.
///
/// The bins do not change with it. What a pair adds to a bin is the r² of
/// [`r2_between`], which is the same to the bit in whichever tile pair it
/// is worked out, and the bins are added up variant by variant in the
/// order of the variants of the pass, which neither a tile nor a block
/// cuts. So this number is a matter of the memory and the time of one
/// step and of nothing a user reads, which
/// `the_bins_do_not_move_with_the_blocks_nor_with_the_tiles` says.
///
/// Nobody has timed this pass: the owner decided on 24 September 2026 that
/// its speed goes to a performance review of its own, as the matrix of
/// every pair did on 23 September 2026, so this number is not a
/// measurement. At 256 variants a pair of tiles holds 512 KB of r² and the
/// six sums of it 3 MB, beside the window itself.
const THE_VARS_OF_A_TILE_OF_THE_WINDOW: usize = 256;

/// The smallest distance a pair is counted at when the user names no
/// number, which is the default of `calc_ld_and_dist_per_pop` in Python and
/// of `calcLdAndDistPerPop` in TypeScript.
///
/// At 1 the only pairs left out are those of two variants at one position,
/// a SNP and an indel at the same base, whose distance is 0. "Its Python
/// function" of `docs/specs/ld.md` gives the number and says that pyNei
/// throws away the pairs at `min_dist` itself, so its default of 1 loses
/// the variants 1 base pair apart with nothing said.
pub const DEFAULT_MIN_DIST: u64 = 1;

/// The largest distance a pair is counted at when the user names no number,
/// which is the default of `calc_ld_and_dist_per_pop` in Python and of
/// `calcLdAndDistPerPop` in TypeScript.
///
/// It is also how far back the window of blocks reaches, so it is what the
/// memory of the pass grows with. "Its Python function" of
/// `docs/specs/ld.md` gives the number and where it comes from: it is
/// plink2's own for `--r2-unphased`, `--ld-window-kb 1000`, where pyNei
/// puts no bound on the distance at all and so counts every pair of every
/// chromosome, 5·10⁹ pairs for a chromosome of 100000 variants.
pub const DEFAULT_MAX_DIST: u64 = 1_000_000;

/// How many bins the distances are cut into when the user names no number,
/// which is the default of `calc_ld_and_dist_per_pop` in Python and of
/// `calcLdAndDistPerPop` in TypeScript.
///
/// "Its Python function" of `docs/specs/ld.md` gives the number, and pyNei
/// has none: it bins no pair and hands over a sample of them instead. At
/// the [`DEFAULT_MAX_DIST`] of a million base pairs each bin is 20000 base
/// pairs wide.
pub const DEFAULT_NUM_BINS: usize = 50;

/// The largest major allele frequency a variant has in a population and is
/// still counted there when the user names no number, which is the default
/// of `calc_ld_and_dist_per_pop` in Python and of `calcLdAndDistPerPop` in
/// TypeScript.
///
/// "Its Python function" of `docs/specs/ld.md` gives the number, and
/// "What it gives" of that item the reason for leaving those variants out:
/// the r² of a variant that hardly varies in a population rests on the one
/// or two individuals that carry the rare allele, and keeping them raises
/// the curve everywhere. It is the frequency pyNei's
/// `calc_rogers_huff_r2_matrix` refuses a whole dataset over.
pub const DEFAULT_MAX_ALLOWED_MAF: f64 = 0.95;

/// What the fall-off of r² with distance is counted with.
#[derive(Debug, Clone, Copy)]
pub struct LdAndDistOptions {
    /// The smallest distance in base pairs at which a pair of variants is
    /// counted, that distance included. A `min_dist` of 1 leaves out only
    /// the pairs of two variants at one position.
    pub min_dist: u64,
    /// The largest distance in base pairs at which a pair of variants is
    /// counted, that distance included. It is also how far back the window
    /// of blocks reaches, so it is what the memory of the pass grows with.
    pub max_dist: u64,
    /// How many bins of equal width the distances from `min_dist` to
    /// `max_dist` are cut into. It is 1 at least.
    pub num_bins: usize,
    /// The largest major allele frequency a variant has in a population
    /// and is still counted there, both ends included. It is worked out
    /// over the individuals of that population alone, so two populations
    /// of one pass count different variants.
    pub max_allowed_maf: f64,
}

/// How the r² of a pair of variants falls off with the distance between
/// them, for each population of a dataset, in bins of distance.
///
/// It reads the reader to its end in one pass, which serves every
/// population, and asks it for the genotypes, the chromosome and the
/// position. The reader is borrowed and not taken, so that whoever built
/// the chain of filters of the pass reads their counts from it when this
/// returns, as `docs/specs/filters.md` says; how many variants the
/// calculation was given is [`LdAndDist::num_vars`].
///
/// `pops` is the indices of the individuals of each population among the
/// individuals of a block, in the order the user gave them, and an empty
/// `pops` is one population of every individual. An individual may be in
/// more than one population, and one in none is read by none of them.
///
/// A pair of variants is counted in a population when both of its variants
/// passed the major allele frequency of that population, when the two are
/// on one chromosome, and when their distance is from `min_dist` to
/// `max_dist`, both included. A pair with no r², which "What it gives" of
/// `docs/specs/ld.md` defines, is in no bin.
///
/// The bins are the same, to the bit, whatever the size of the blocks the
/// reader gives: the r² of a pair is worked out over the individuals,
/// which no block and no tile cuts, and the bins are added up in the order
/// of the variants of the pass.
///
/// # Errors
///
/// [`Error::LdMinDistAboveMaxDist`] for a `min_dist` above `max_dist`;
/// [`Error::LdNoBins`] for a `num_bins` of 0;
/// [`Error::LdMaxAllowedMafOutOfRange`] for a `max_allowed_maf` that is
/// NaN or is not from 0 to 1; [`Error::LdPopWithNoIndividual`] for a
/// population that names no individual;
/// [`Error::LdIndividualNotInTheDataset`] and
/// [`Error::LdIndividualAskedForTwice`] for the individuals of a
/// population; [`Error::PassGaveNoVariant`] when the reader gives no
/// variant; [`Error::FieldsNotInTheBlock`] when a block holds variants and
/// no genotypes or no position; [`Error::LdNoMemory`] when this machine
/// does not give the memory of the bins, of the window or of the r² of a
/// step, which is asked of it with `try_reserve_exact` and not taken; what
/// the dosages of a block and the r² of two tiles refuse; and whatever the
/// reader fails with, which is given on as it is.
pub fn calc_ld_and_dist<R: BlockReader + ?Sized>(
    reader: &mut R,
    pops: &[&[usize]],
    options: &LdAndDistOptions,
) -> Result<LdAndDist> {
    the_ld_and_dist_in_tiles_of(reader, pops, options, THE_VARS_OF_A_TILE_OF_THE_WINDOW)
}

/// The fall-off of r² with distance, with the pairs of each step taken in
/// tiles of `vars_per_tile` variants.
///
/// [`calc_ld_and_dist`] is this with the tile of the module, and the tests
/// are what give another: the bins are the same, to the bit, whatever the
/// tile, because a pair adds the same r² to the same bin in whichever tile
/// pair it is worked out and the bins are added up in the order of the
/// variants.
///
/// # Errors
///
/// Those of [`calc_ld_and_dist`].
fn the_ld_and_dist_in_tiles_of<R: BlockReader + ?Sized>(
    reader: &mut R,
    pops: &[&[usize]],
    options: &LdAndDistOptions,
    vars_per_tile: usize,
) -> Result<LdAndDist> {
    if options.min_dist > options.max_dist {
        return Err(Error::LdMinDistAboveMaxDist {
            min_dist: options.min_dist,
            max_dist: options.max_dist,
        });
    }
    if options.num_bins == 0 {
        return Err(Error::LdNoBins);
    }
    let mut of_the_pops = TheDosagesOfThePops::of(pops, options.max_dist, options.max_allowed_maf)?;
    let mut of_each_pop: Vec<LdBins> = Vec::new();
    of_each_pop
        .try_reserve_exact(of_the_pops.num_pops())
        .map_err(|_| Error::LdNoMemory {
            what: "the bins of each population",
            values: of_the_pops.num_pops(),
        })?;
    for _ in 0..of_the_pops.num_pops() {
        of_each_pop.push(LdBins::of(options)?);
    }
    // The genotypes, the chromosome and the position are what this reads,
    // so a reader over a file leaves the other columns of a variant
    // unparsed.
    reader.set_needs(Needs::GTS | Needs::CHROM_POS);
    // A tile of no variant would take no variant of the window and the
    // step would stand still.
    let mut of_the_pairs = ThePairsOfAStep::of(options, vars_per_tile.max(1));
    let mut num_vars = 0_u64;
    while let Some(block) = reader.next_block()? {
        // Every variant counted here was read from a source, and a u64
        // counts 1.8e19 of them: a pass that passed this number would have
        // read more bytes than any storage holds.
        num_vars = num_vars.saturating_add(the_count_of(block.num_vars));
        of_the_pops.take_the_block(block)?;
        for (of_the_pop, bins) in of_the_pops.the_pops().iter().zip(&mut of_each_pop) {
            of_the_pairs.the_pairs_of_the_newest_variants(of_the_pop, bins)?;
        }
    }
    if num_vars == 0 {
        let filters = reader.filtering_stats();
        return Err(Error::PassGaveNoVariant {
            // The filter nearest the source was given what the source
            // gave; with no filter the pass gave what the source gave,
            // which is nothing.
            num_vars_of_the_source: filters.last().map_or(0, |(_, stats)| stats.vars_processed),
            filters,
        });
    }
    for (of_the_pop, bins) in of_the_pops.the_pops().iter().zip(&mut of_each_pop) {
        bins.num_vars = of_the_pop.num_vars();
    }
    Ok(LdAndDist {
        num_vars,
        of_each_pop,
    })
}

/// How the r² of a pair of variants falls off with the distance between
/// them, for each population of one pass.
///
/// It is what [`calc_ld_and_dist`] gives.
#[derive(Debug)]
pub struct LdAndDist {
    /// How many variants the calculation was given, before the major
    /// allele frequency of any population.
    num_vars: u64,
    /// The bins of each population, in the order the populations were
    /// given.
    of_each_pop: Vec<LdBins>,
}

impl LdAndDist {
    /// The variants the calculation was given, before the major allele
    /// frequency of any population.
    #[must_use]
    pub fn num_vars(&self) -> u64 {
        self.num_vars
    }

    /// How many populations the pass counted, which is 1 when it was given
    /// none.
    #[must_use]
    pub fn num_pops(&self) -> usize {
        self.of_each_pop.len()
    }

    /// The bins of the population at that position among the ones given,
    /// and `None` when `pop` is not a population of the call.
    #[must_use]
    pub fn bins_of_pop(&self, pop: usize) -> Option<&LdBins> {
        self.of_each_pop.get(pop)
    }
}

/// The pairs of one population counted in bins of the distance of their
/// two variants.
///
/// The bins are of equal width across the distances from `min_dist` to
/// `max_dist`, both included: the width is
/// (`max_dist` − `min_dist` + 1) / `num_bins` base pairs, and a pair at
/// the distance d falls in the bin that floor((d − `min_dist`) / width)
/// gives, the last bin taking anything the rounding would put past it.
#[derive(Debug)]
pub struct LdBins {
    /// The smallest distance a pair is counted at, that distance included.
    min_dist: u64,
    /// The largest distance a pair is counted at, that distance included.
    max_dist: u64,
    /// How many variants this population kept at its major allele
    /// frequency over the whole pass.
    num_vars: u64,
    /// How many pairs each bin holds, one for each bin.
    num_pairs: Vec<u64>,
    /// The sum of the r² of the pairs of each bin.
    sum_r2: Vec<f64>,
    /// The sum of the squares of the r² of the pairs of each bin, which
    /// the standard deviation is taken from.
    sum_of_squares: Vec<f64>,
}

impl LdBins {
    /// Empty bins of the distances the options ask for.
    ///
    /// # Errors
    ///
    /// [`Error::LdNoMemory`] when this machine does not give the memory of
    /// the three counts of each bin.
    fn of(options: &LdAndDistOptions) -> Result<LdBins> {
        let num_bins = options.num_bins;
        let of_a_bin = |what: &'static str| the_memory_for(what, num_bins);
        Ok(LdBins {
            min_dist: options.min_dist,
            max_dist: options.max_dist,
            num_vars: 0,
            num_pairs: a_vector_of(0, num_bins, &of_a_bin("the pairs of each bin"))?,
            sum_r2: a_vector_of(0.0, num_bins, &of_a_bin("the sum of r² of each bin"))?,
            sum_of_squares: a_vector_of(
                0.0,
                num_bins,
                &of_a_bin("the sum of the squares of r² of each bin"),
            )?,
        })
    }

    /// How many bins the distances were cut into.
    #[must_use]
    pub fn num_bins(&self) -> usize {
        self.num_pairs.len()
    }

    /// The variants that passed the major allele frequency of this
    /// population over the whole pass.
    #[must_use]
    pub fn num_vars(&self) -> u64 {
        self.num_vars
    }

    /// The smallest and the largest distance of the bin, both included,
    /// and `None` when `bin` is not a bin.
    #[must_use]
    pub fn bounds(&self, bin: usize) -> Option<(u64, u64)> {
        let last = self.num_bins().checked_sub(1)?;
        if bin > last {
            return None;
        }
        let smallest = self.the_smallest_dist_of(bin);
        let largest = match bin == last {
            true => self.max_dist,
            false => self
                .the_smallest_dist_of(bin.saturating_add(1))
                .saturating_sub(1)
                .min(self.max_dist),
        };
        Some((smallest, largest))
    }

    /// How many pairs the bin holds, and `None` when `bin` is not a bin.
    #[must_use]
    pub fn num_pairs(&self, bin: usize) -> Option<u64> {
        self.num_pairs.get(bin).copied()
    }

    /// The mean of the r² of the pairs of the bin, and `None` when `bin`
    /// is not a bin and when it holds no pair, which the binding crates
    /// give the user as NaN.
    #[must_use]
    pub fn mean_r2(&self, bin: usize) -> Option<f64> {
        let num_pairs = self.num_pairs.get(bin).copied()?;
        if num_pairs == 0 {
            return None;
        }
        Some(self.sum_r2.get(bin)? / num_pairs as f64)
    }

    /// The standard deviation of the r² of the pairs of the bin, with the
    /// pairs of the bin as the divisor, and `None` when `bin` is not a bin
    /// and when it holds no pair. A bin of one pair has 0.
    #[must_use]
    pub fn sd_r2(&self, bin: usize) -> Option<f64> {
        let num_pairs = self.num_pairs.get(bin).copied()?;
        if num_pairs == 0 {
            return None;
        }
        let mean = self.sum_r2.get(bin)? / num_pairs as f64;
        let of_the_squares = self.sum_of_squares.get(bin)? / num_pairs as f64;
        // The mean of the squares less the square of the mean is the
        // variance, and it is 0 to the last bits when every pair of the
        // bin holds the same r², where the subtraction can leave a value
        // below 0.
        Some((of_the_squares - mean * mean).max(0.0).sqrt())
    }

    /// Counts one pair of the distance `dist`, whose r² is `r2`.
    ///
    /// The caller has found the distance to be from `min_dist` to
    /// `max_dist` and the r² to be a number.
    fn the_pair_is_counted(&mut self, dist: u64, r2: f64) {
        let bin = self.the_bin_of(dist);
        if let Some(num_pairs) = self.num_pairs.get_mut(bin) {
            // The pairs of a pass are at most its variants times the
            // variants of one window, and a u64 counts 1.8e19 of them.
            *num_pairs = num_pairs.saturating_add(1);
        }
        if let Some(sum) = self.sum_r2.get_mut(bin) {
            *sum += r2;
        }
        if let Some(sum) = self.sum_of_squares.get_mut(bin) {
            *sum += r2 * r2;
        }
    }

    /// How many distances the bins cut into `num_bins` parts, which is
    /// `max_dist` − `min_dist` + 1 and 1 at least.
    fn the_distances_counted(&self) -> u128 {
        // The caller of the pass refused a min_dist above max_dist, and
        // the largest difference of two u64 and one more is far below
        // what a u128 holds.
        u128::from(self.max_dist.abs_diff(self.min_dist)).saturating_add(1)
    }

    /// The bin a pair of that distance falls in, which is the last bin for
    /// a distance the arithmetic would put past them.
    ///
    /// It is floor((d − `min_dist`) · `num_bins` / the distances counted)
    /// in whole numbers, which is the arithmetic of "What it gives" of
    /// `docs/specs/ld.md` with the width neither divided by nor rounded.
    /// A width of 9/7 base pairs, 18 distances in 14 bins, is no `f64`,
    /// and a distance divided by such a width and a bin multiplied by it
    /// round the other way from each other, so a distance that is the
    /// smallest of a bin lands in the bin below the one
    /// [`LdBins::bounds`] answers for it.
    fn the_bin_of(&self, dist: u64) -> usize {
        let last = self.num_bins().saturating_sub(1);
        // The caller counts a pair of a distance from min_dist to
        // max_dist, so this is the distance from the first bin's own.
        let from_the_first = u128::from(dist.abs_diff(self.min_dist));
        let bin = from_the_first
            .checked_mul(u128::from(the_count_of(self.num_bins())))
            .and_then(|of_the_bins| of_the_bins.checked_div(self.the_distances_counted()));
        match bin {
            Some(bin) => usize::try_from(bin).unwrap_or(last).min(last),
            None => last,
        }
    }

    /// The smallest distance that falls in the bin, which is
    /// `min_dist` + ceil(`bin` · the distances counted / `num_bins`).
    ///
    /// That is the smallest whole distance whose [`LdBins::the_bin_of`]
    /// is `bin` or above, so the bounds a bin reports and the bin a pair
    /// is counted in agree by construction.
    fn the_smallest_dist_of(&self, bin: usize) -> u64 {
        let num_bins = u128::from(the_count_of(self.num_bins()));
        if num_bins == 0 {
            return self.min_dist;
        }
        let past_the_first = u128::from(the_count_of(bin))
            .checked_mul(self.the_distances_counted())
            .map(|of_the_range| of_the_range.div_ceil(num_bins))
            .and_then(|past_the_first| u64::try_from(past_the_first).ok());
        match past_the_first {
            Some(past_the_first) => self.min_dist.saturating_add(past_the_first),
            None => self.max_dist,
        }
    }
}

/// The pairs of one step of a pass: the variants a population read with
/// the newest block against the variants it holds.
///
/// It keeps the two buffers of a step with the memory they have, so a pass
/// allocates them once and grows them to the largest step it meets, and it
/// holds the distances at which a pair is counted.
struct ThePairsOfAStep {
    /// The smallest distance a pair is counted at, that distance included.
    min_dist: u64,
    /// The largest distance a pair is counted at, that distance included.
    max_dist: u64,
    /// How many variants a tile of the products holds, 1 at least.
    vars_per_tile: usize,
    /// The r² of one pair of tiles, one row for each variant of the tile
    /// of columns.
    of_a_tile_pair: Vec<f64>,
    /// The r² of one tile of columns against every variant the population
    /// holds that can pair with the first of them: one row for each of
    /// those columns, so that the bins are added up column by column and
    /// not tile pair by tile pair.
    of_a_tile_of_columns: Vec<f64>,
}

impl ThePairsOfAStep {
    /// The buffers of a pass that counts the pairs of those distances in
    /// tiles of `vars_per_tile` variants.
    fn of(options: &LdAndDistOptions, vars_per_tile: usize) -> ThePairsOfAStep {
        ThePairsOfAStep {
            min_dist: options.min_dist,
            max_dist: options.max_dist,
            vars_per_tile,
            of_a_tile_pair: Vec::new(),
            of_a_tile_of_columns: Vec::new(),
        }
    }

    /// Counts into `bins` the pairs that hold one of the variants the
    /// population kept of the newest block.
    ///
    /// Every pair of the pass is counted at the step of the newer of its
    /// two variants, and the variants held at that step are every variant
    /// within `max_dist` of it: the blocks that fell out of the window
    /// when the newest block was taken are dropped when the next one is,
    /// which [`TheDosagesOfThePops::take_the_block`] says why.
    ///
    /// # Errors
    ///
    /// [`Error::LdRowsNotInTheDosages`] when the variants of a tile are
    /// not variants of the dosages the population holds, which is a defect
    /// of the tiling and not anything a caller of the crate wrote;
    /// [`Error::LdNoMemory`] when this machine does not give the memory of
    /// the r² of a step; and what [`r2_between`] refuses.
    fn the_pairs_of_the_newest_variants(
        &mut self,
        pop: &ThePopOverTheWindow,
        bins: &mut LdBins,
    ) -> Result<()> {
        let Some(dosages) = pop.dosages() else {
            // No block of the pass has been taken yet.
            return Ok(());
        };
        let held = pop.variants();
        let of_other_variants = || Error::LdRowsNotInTheDosages {
            first: 0,
            asked_for: held.len(),
            num_vars: dosages.num_vars(),
        };
        if dosages.num_vars() != held.len() {
            return Err(of_other_variants());
        }
        // The variants the newest block gave this population, which are
        // the last ones it holds.
        let Some(first_new) = held.len().checked_sub(pop.kept_of_the_newest_block()) else {
            return Err(of_other_variants());
        };
        let first_var = pop.first_var();
        let mut col = first_new;
        while col < held.len() {
            let end = the_end_of_the_tile(first_var, col, held.len(), self.vars_per_tile);
            self.the_pairs_of_a_tile_of_columns(dosages, held, first_var, (col, end), bins)?;
            col = end;
        }
        Ok(())
    }

    /// Counts into `bins` the pairs of the variants `col` to `col_end` of
    /// the window against the variants of the window before each of them.
    ///
    /// The r² of the whole tile of columns against every variant that can
    /// pair with it is worked out first, tile pair by tile pair, and the
    /// bins are added up afterwards one column at a time, the variants of
    /// a column in the order of the pass. That order is the order of the
    /// variants and not the order of the tile pairs, so it does not move
    /// when a block ends inside a tile of columns.
    ///
    /// # Errors
    ///
    /// Those of [`ThePairsOfAStep::the_pairs_of_the_newest_variants`].
    fn the_pairs_of_a_tile_of_columns(
        &mut self,
        dosages: &LdDosages,
        held: &[TheVariantOfTheWindow],
        first_var: u64,
        (col, col_end): (usize, usize),
        bins: &mut LdBins,
    ) -> Result<()> {
        let ThePairsOfAStep {
            min_dist,
            max_dist,
            vars_per_tile,
            of_a_tile_pair,
            of_a_tile_of_columns,
        } = self;
        let of_other_variants = || Error::LdRowsNotInTheDosages {
            first: col,
            asked_for: col_end.saturating_sub(col),
            num_vars: held.len(),
        };
        // The first variant held that any column of the tile pairs with.
        // It is taken over every column and not over the first of them:
        // nothing asks the chromosomes of a source to be grouped, so a
        // later column can be on a chromosome the first column is not on
        // and reach back past what the first column reaches.
        let row_start = the_first_row_in_reach(held, (col, col_end), *max_dist);
        let (Some(num_cols), Some(num_rows)) =
            (col_end.checked_sub(col), col_end.checked_sub(row_start))
        else {
            return Err(of_other_variants());
        };
        let values = num_cols
            .checked_mul(num_rows)
            .ok_or_else(of_other_variants)?;
        the_buffer_of(
            of_a_tile_of_columns,
            values,
            &the_memory_for("the r² of a tile of columns", values),
        )?;
        let of_the_columns = dosages.rows(col, num_cols)?;
        let mut row = row_start;
        while row < col_end {
            let row_end = the_end_of_the_tile(first_var, row, col_end, *vars_per_tile);
            let (Some(of_the_tile), Some(at)) =
                (row_end.checked_sub(row), row.checked_sub(row_start))
            else {
                return Err(of_other_variants());
            };
            let values = num_cols
                .checked_mul(of_the_tile)
                .ok_or_else(of_other_variants)?;
            the_buffer_of(
                of_a_tile_pair,
                values,
                &the_memory_for("the r² of a pair of tiles", values),
            )?;
            match row == col && row_end == col_end {
                // The tile of the rows is the tile of the columns, and the
                // two arguments are then one reference given twice, which
                // is what the four products of a set against itself are
                // taken on.
                true => r2_between(&of_the_columns, &of_the_columns, of_a_tile_pair)?,
                false => {
                    let of_the_rows = dosages.rows(row, of_the_tile)?;
                    r2_between(&of_the_columns, &of_the_rows, of_a_tile_pair)?;
                }
            }
            // The r² of one column of the tile pair lies row after row,
            // and the buffer of the columns holds one row for each column
            // of the tile: each of them is a run of it.
            for (column, of_the_column) in of_a_tile_pair.chunks_exact(of_the_tile).enumerate() {
                let Some(from) = column
                    .checked_mul(num_rows)
                    .and_then(|from| from.checked_add(at))
                else {
                    return Err(of_other_variants());
                };
                let Some(to) = from.checked_add(of_the_tile) else {
                    return Err(of_other_variants());
                };
                let Some(into) = of_a_tile_of_columns.get_mut(from..to) else {
                    return Err(of_other_variants());
                };
                into.copy_from_slice(of_the_column);
            }
            row = row_end;
        }
        let Some(before_the_columns) = held.get(row_start..) else {
            return Err(of_other_variants());
        };
        for (column, of_the_column) in of_a_tile_of_columns.chunks_exact(num_rows).enumerate() {
            let at_the_column = col.saturating_add(column);
            let (Some(of_the_column_variant), Some(rows_before)) = (
                held.get(at_the_column),
                at_the_column.checked_sub(row_start),
            ) else {
                return Err(of_other_variants());
            };
            let pairs = of_the_column
                .iter()
                .zip(before_the_columns)
                .take(rows_before);
            for (r2, variant) in pairs {
                if variant.chrom != of_the_column_variant.chrom || r2.is_nan() {
                    continue;
                }
                let dist = of_the_column_variant.pos.abs_diff(variant.pos);
                if dist < *min_dist || dist > *max_dist {
                    continue;
                }
                bins.the_pair_is_counted(dist, *r2);
            }
        }
        Ok(())
    }
}

/// Where the tile that holds the variant `at` of the window ends, as an
/// index into the variants held, and `num_held` when the tile runs past
/// them.
///
/// The tiles are cut at multiples of `vars_per_tile` counted from the
/// first variant of the pass, which `first_var` says how many of have
/// fallen out of the window, so a tile holds the same variants whatever
/// the blocks the reader gave. It is one variant past `at` at least.
fn the_end_of_the_tile(first_var: u64, at: usize, num_held: usize, vars_per_tile: usize) -> usize {
    let of_a_tile = the_count_of(vars_per_tile);
    // The variants a population keeps are counted in a u64, and a pass
    // that passed what one holds would have read more bytes than any
    // storage gives.
    let of_the_pass = first_var.saturating_add(the_count_of(at));
    let end = of_the_pass
        .checked_div(of_a_tile)
        .and_then(|tile| tile.checked_add(1))
        .and_then(|tile| tile.checked_mul(of_a_tile))
        .unwrap_or(u64::MAX)
        .saturating_sub(first_var);
    usize::try_from(end).unwrap_or(num_held).min(num_held)
}

/// The first variant of the window that is on the chromosome of one of
/// the variants `col` to `col_end` and within `max_dist` of it, and `col`
/// itself when none is.
///
/// The rows of a tile of columns run from this variant to the last column
/// of the tile, so every pair that holds a column of the tile and a
/// variant before it is worked out. The columns of a tile are in the
/// order of the source, which no chromosome orders, so this is taken over
/// all of them: a source whose chromosomes are not grouped can put a
/// later column on a chromosome that the first column of the tile is not
/// on, and that column reaches back to a variant the first column does
/// not.
fn the_first_row_in_reach(
    held: &[TheVariantOfTheWindow],
    (col, col_end): (usize, usize),
    max_dist: u64,
) -> usize {
    let Some(of_the_columns) = held.get(col..col_end.max(col)) else {
        return col;
    };
    held.get(..col)
        .unwrap_or_default()
        .iter()
        .position(|variant| {
            of_the_columns.iter().any(|of_the_column| {
                variant.chrom == of_the_column.chrom
                    && of_the_column.pos.abs_diff(variant.pos) <= max_dist
            })
        })
        .unwrap_or(col)
}

/// Makes `buffer` hold `values` values, keeping the memory it has.
///
/// The memory it has not is asked of this machine with `try_reserve` and
/// not taken, so a step this machine cannot hold the r² of is an error and
/// not a process that ends.
///
/// # Errors
///
/// What `not_given` gives, when the machine does not give the memory.
fn the_buffer_of(
    buffer: &mut Vec<f64>,
    values: usize,
    not_given: &impl Fn() -> Error,
) -> Result<()> {
    buffer.clear();
    if let Some(more) = values.checked_sub(buffer.capacity()) {
        buffer.try_reserve(more).map_err(|_| not_given())?;
    }
    buffer.resize(values, 0.0);
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::{BufReader, Cursor};
    use std::path::{Path, PathBuf};

    use super::{
        LdAndDist, LdAndDistOptions, LdBins, TheDosagesOfThePops, ThePopOverTheWindow,
        TheVariantOfTheWindow, TheWindowOfTheBlocks, calc_ld_and_dist, the_ld_and_dist_in_tiles_of,
        the_variants_of,
    };

    use crate::block::{Block, BlockReader};
    use crate::error::Error;
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{ChromTable, Needs};

    /// What the window held after one block of a pass was taken.
    #[derive(Debug, PartialEq, Eq)]
    struct TheStepOfAPass {
        /// The chromosome and the position of each variant held, the
        /// oldest first.
        held: Vec<(String, u64)>,
        /// How many blocks fell out when this block was taken.
        dropped: usize,
        /// How many blocks the window held.
        num_blocks: usize,
        /// What the window answered for how many variants it held.
        num_vars: usize,
    }

    /// The chromosome and the position of each variant of a list, as a
    /// step of a pass holds them.
    fn at(variants: &[(&str, u64)]) -> Vec<(String, u64)> {
        variants
            .iter()
            .map(|(chrom, pos)| ((*chrom).to_owned(), *pos))
            .collect()
    }

    /// A VCF of one diploid individual with a variant at each of
    /// `variants`, each of them the name of a chromosome and a position.
    fn vcf_of(variants: &[(&str, u64)]) -> Vec<u8> {
        let mut vcf = String::from(
            "##fileformat=VCFv4.2\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti0\n",
        );
        for (chrom, pos) in variants {
            vcf.push_str(&format!("{chrom}\t{pos}\tvar\tA\tC\t.\t.\t.\tGT\t0/1\n"));
        }
        vcf.into_bytes()
    }

    /// What a window of `max_dist` held after each block of a pass over
    /// the bytes of a VCF read in blocks of `num_vars_per_block` variants.
    fn the_steps_of_a_pass(
        vcf: &[u8],
        num_vars_per_block: usize,
        max_dist: u64,
    ) -> Vec<TheStepOfAPass> {
        let options = VcfOptions {
            num_vars_per_block: Some(num_vars_per_block),
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::new(Cursor::new(vcf.to_vec()), options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        };
        // What the pass over the window reads, as `calc_r2_matrix` asks
        // for it.
        reader.set_needs(Needs::GTS | Needs::CHROM_POS);
        let mut window = TheWindowOfTheBlocks::of(max_dist);
        let mut steps = Vec::new();
        loop {
            let block = match reader.next_block() {
                Ok(Some(block)) => block,
                Ok(None) => break,
                Err(error) => panic!("the block was not read: {error}"),
            };
            let dropped = match window.take_the_block(block) {
                Ok(dropped) => dropped,
                Err(error) => panic!("the window refused the block: {error}"),
            };
            let held = window
                .variants()
                .map(|variant| (the_name_of(reader.chroms(), variant.chrom), variant.pos))
                .collect();
            steps.push(TheStepOfAPass {
                held,
                dropped,
                num_blocks: window.num_blocks(),
                num_vars: window.num_vars(),
            });
        }
        steps
    }

    /// The name of the chromosome `number` in the table of a reader.
    fn the_name_of(chroms: &ChromTable, number: u32) -> String {
        match chroms.name(number) {
            Some(name) => name.to_owned(),
            None => panic!("the table has no chromosome {number}"),
        }
    }

    /// The seven variants of the worked window below: two close together,
    /// one a thousand base pairs off the first, one a base pair past that,
    /// one far from every other, and two on a second chromosome.
    const THE_VARIANTS_OF_THE_EXAMPLE: [(&str, u64); 7] = [
        ("chr1", 100),
        ("chr1", 200),
        ("chr1", 1100),
        ("chr1", 1101),
        ("chr1", 5000),
        ("chr2", 100),
        ("chr2", 200),
    ];

    /// How far back the worked window reaches, which the variant at 1100
    /// is exactly that far from the one at 100 for.
    const THE_MAX_DIST_OF_THE_EXAMPLE: u64 = 1000;

    #[test]
    fn a_variant_exactly_at_max_dist_is_held_and_one_a_base_pair_past_it_falls_out() {
        let vcf = vcf_of(&THE_VARIANTS_OF_THE_EXAMPLE);
        let steps = the_steps_of_a_pass(&vcf, 1, THE_MAX_DIST_OF_THE_EXAMPLE);
        assert_eq!(steps.len(), 7);
        // 1100 is 1000 from 100, which is `max_dist`, so nothing has left.
        assert_eq!(
            steps[2].held,
            at(&[("chr1", 100), ("chr1", 200), ("chr1", 1100)])
        );
        // 1101 is 1001 from 100, which is one past it, and 901 from 200.
        assert_eq!(
            steps[3].held,
            at(&[("chr1", 200), ("chr1", 1100), ("chr1", 1101)])
        );
    }

    #[test]
    fn every_block_of_one_variant_falls_out_as_soon_as_it_is_out_of_reach() {
        let vcf = vcf_of(&THE_VARIANTS_OF_THE_EXAMPLE);
        let steps = the_steps_of_a_pass(&vcf, 1, THE_MAX_DIST_OF_THE_EXAMPLE);
        let dropped: Vec<usize> = steps.iter().map(|step| step.dropped).collect();
        // The three of chr1 that 5000 left behind go at its step, and the
        // one that is left goes when the first variant of chr2 is read.
        assert_eq!(dropped, vec![0, 0, 0, 1, 3, 1, 0]);
        assert_eq!(steps[4].held, at(&[("chr1", 5000)]));
        assert_eq!(steps[5].held, at(&[("chr2", 100)]));
        assert_eq!(steps[6].held, at(&[("chr2", 100), ("chr2", 200)]));
    }

    #[test]
    fn a_block_is_held_while_any_one_of_its_variants_is_in_reach() {
        let vcf = vcf_of(&THE_VARIANTS_OF_THE_EXAMPLE);
        let steps = the_steps_of_a_pass(&vcf, 2, THE_MAX_DIST_OF_THE_EXAMPLE);
        assert_eq!(steps.len(), 4);
        // The newest variant is 1101, which 100 is 1001 from and 200 is
        // 901 from, so the block of the two is held for the second.
        assert_eq!(
            steps[1].held,
            at(&[("chr1", 100), ("chr1", 200), ("chr1", 1100), ("chr1", 1101)])
        );
        assert_eq!(steps[1].dropped, 0);
        assert_eq!(steps[1].num_blocks, 2);
    }

    #[test]
    fn a_block_whose_variants_straddle_a_chromosome_change_is_held_for_the_ones_on_the_new_one() {
        let vcf = vcf_of(&THE_VARIANTS_OF_THE_EXAMPLE);
        let steps = the_steps_of_a_pass(&vcf, 2, THE_MAX_DIST_OF_THE_EXAMPLE);
        // The third block is chr1 5000 and chr2 100, and the newest
        // variant is the second of the two, so the block is held and the
        // two older ones, every variant of which is on chr1, are dropped.
        assert_eq!(steps[2].held, at(&[("chr1", 5000), ("chr2", 100)]));
        assert_eq!(steps[2].dropped, 2);
        assert_eq!(steps[2].num_blocks, 1);
        // And it is still held at the next block, for its variant of chr2.
        assert_eq!(
            steps[3].held,
            at(&[("chr1", 5000), ("chr2", 100), ("chr2", 200)])
        );
        assert_eq!(steps[3].dropped, 0);
    }

    #[test]
    fn a_block_every_variant_of_which_has_fallen_out_is_dropped_whole() {
        let vcf = vcf_of(&THE_VARIANTS_OF_THE_EXAMPLE);
        let steps = the_steps_of_a_pass(&vcf, 3, THE_MAX_DIST_OF_THE_EXAMPLE);
        assert_eq!(steps.len(), 3);
        // The first block is 100, 200 and 1100, and when the second one
        // ends at chr2 100 not one of its three is in reach.
        assert_eq!(steps[1].dropped, 1);
        assert_eq!(
            steps[1].held,
            at(&[("chr1", 1101), ("chr1", 5000), ("chr2", 100)])
        );
    }

    #[test]
    fn the_window_counts_the_variants_of_the_blocks_it_holds() {
        let vcf = vcf_of(&THE_VARIANTS_OF_THE_EXAMPLE);
        for num_vars_per_block in [1, 2, 3, 7, 64] {
            let steps = the_steps_of_a_pass(&vcf, num_vars_per_block, THE_MAX_DIST_OF_THE_EXAMPLE);
            for step in &steps {
                assert_eq!(
                    step.num_vars,
                    step.held.len(),
                    "at blocks of {num_vars_per_block} variants"
                );
            }
        }
    }

    #[test]
    fn a_block_that_holds_the_whole_of_a_dataset_is_held_for_its_newest_variant() {
        let vcf = vcf_of(&THE_VARIANTS_OF_THE_EXAMPLE);
        let steps = the_steps_of_a_pass(&vcf, 500, THE_MAX_DIST_OF_THE_EXAMPLE);
        assert_eq!(steps.len(), 1);
        assert_eq!(steps[0].num_blocks, 1);
        assert_eq!(steps[0].held, at(&THE_VARIANTS_OF_THE_EXAMPLE));
    }

    /// A VCF of 600 variants of one chromosome, a thousand base pairs
    /// apart, the first at 1000 and the last at 600000.
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "600 positions of at most 600000, built here"
    )]
    fn the_vcf_of_a_long_chromosome() -> Vec<u8> {
        let variants: Vec<(&str, u64)> = (1..=600).map(|var| ("chr1", var * 1000)).collect();
        vcf_of(&variants)
    }

    #[test]
    fn the_variants_within_max_dist_are_held_at_blocks_of_seven_sixty_four_and_five_hundred() {
        let vcf = the_vcf_of_a_long_chromosome();
        // The newest variant of the pass is at 600000 and the window
        // reaches back to 590000, so the eleven variants from 590000 on
        // are in reach whatever the size of the blocks, and the block
        // that holds the one at 590000 is held whole with them: of the
        // seven from 589000 on, of the sixty-four from 577000 on, or of
        // the five hundred from 501000 on.
        let what_is_held = [
            (7, 2, 12, 589_000),
            (64, 1, 24, 577_000),
            (500, 1, 100, 501_000),
        ];
        for (num_vars_per_block, num_blocks, num_vars, first_pos) in what_is_held {
            let steps = the_steps_of_a_pass(&vcf, num_vars_per_block, 10_000);
            let last = match steps.last() {
                Some(last) => last,
                None => panic!("the pass gave no block"),
            };
            assert_eq!(
                (last.num_blocks, last.num_vars),
                (num_blocks, num_vars),
                "at blocks of {num_vars_per_block} variants"
            );
            assert_eq!(
                last.held.first(),
                Some(&("chr1".to_owned(), first_pos)),
                "at blocks of {num_vars_per_block} variants"
            );
            assert_eq!(
                last.held.last(),
                Some(&("chr1".to_owned(), 600_000)),
                "at blocks of {num_vars_per_block} variants"
            );
        }
    }

    #[test]
    fn a_max_dist_of_zero_holds_the_variants_at_the_position_of_the_newest_one() {
        let vcf = vcf_of(&[("chr1", 100), ("chr1", 100), ("chr1", 101)]);
        let steps = the_steps_of_a_pass(&vcf, 1, 0);
        assert_eq!(steps[1].held, at(&[("chr1", 100), ("chr1", 100)]));
        assert_eq!(steps[2].held, at(&[("chr1", 101)]));
        assert_eq!(steps[2].dropped, 2);
    }

    #[test]
    fn a_window_that_has_taken_no_block_holds_nothing_and_has_no_newest_variant() {
        let window = TheWindowOfTheBlocks::of(1000);
        assert_eq!(window.num_blocks(), 0);
        assert_eq!(window.num_vars(), 0);
        assert_eq!(window.variants().count(), 0);
        assert_eq!(window.the_newest_variant(), None);
        assert_eq!(window.blocks().len(), 0);
    }

    /// A block of `num_vars` variants of one diploid individual on the
    /// chromosome `chrom`, at the positions `poss`, whose genotypes are
    /// the reference allele.
    fn block_of(num_vars: usize, chrom: Option<Vec<u32>>, poss: Option<Vec<u64>>) -> Block {
        Block {
            num_vars,
            num_individuals: 1,
            ploidy: 2,
            gts: vec![0; num_vars.saturating_mul(2)],
            chrom,
            pos: poss,
            id: None,
            alleles: None,
            qual: None,
        }
    }

    #[test]
    fn a_block_without_the_chromosome_and_the_position_of_its_variants_is_refused() {
        let mut window = TheWindowOfTheBlocks::of(1000);
        let error = match window.take_the_block(block_of(2, None, None)) {
            Ok(dropped) => panic!("the window took a block of no position and dropped {dropped}"),
            Err(error) => error,
        };
        assert!(
            matches!(error, Error::FieldsNotInTheBlock { fields } if fields == Needs::CHROM_POS),
            "{error}"
        );
        assert_eq!(window.num_blocks(), 0);
    }

    #[test]
    fn a_block_whose_position_column_is_not_of_its_number_of_variants_is_refused() {
        let mut window = TheWindowOfTheBlocks::of(1000);
        let block = block_of(2, Some(vec![0, 0]), Some(vec![100]));
        let error = match window.take_the_block(block) {
            Ok(dropped) => panic!("the window took a block of one position and dropped {dropped}"),
            Err(error) => error,
        };
        assert!(
            matches!(
                error,
                Error::BlockArrayOfAnotherSize {
                    array: "pos",
                    found: 1,
                    expected: 2
                }
            ),
            "{error}"
        );
        assert_eq!(window.num_blocks(), 0);
    }

    #[test]
    fn a_block_of_no_variant_leaves_the_newest_variant_where_it_was() {
        let mut window = TheWindowOfTheBlocks::of(1000);
        let first = block_of(1, Some(vec![0]), Some(vec![100]));
        match window.take_the_block(first) {
            Ok(dropped) => assert_eq!(dropped, 0),
            Err(error) => panic!("the window refused the block: {error}"),
        }
        let none = block_of(0, Some(Vec::new()), Some(Vec::new()));
        match window.take_the_block(none) {
            Ok(dropped) => assert_eq!(dropped, 0),
            Err(error) => panic!("the window refused the block of no variant: {error}"),
        }
        assert_eq!(
            window.the_newest_variant(),
            Some(TheVariantOfTheWindow { chrom: 0, pos: 100 })
        );
        assert_eq!(window.num_vars(), 1);
        assert_eq!(window.num_blocks(), 2);
    }

    #[test]
    fn the_variants_of_a_block_without_the_two_columns_are_none() {
        let block = block_of(3, None, None);
        assert_eq!(the_variants_of(&block).count(), 0);
    }

    /// A VCF of four diploid individuals, `i0` to `i3`, with a variant at
    /// each of `variants`: the name of a chromosome, a position and the
    /// genotype of each of the four.
    fn vcf_of_four_individuals(variants: &[(&str, u64, [&str; 4])]) -> Vec<u8> {
        let mut vcf = String::from(
            "##fileformat=VCFv4.2\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti0\ti1\ti2\ti3\n",
        );
        for (chrom, pos, of_the_individuals) in variants {
            let gts = of_the_individuals.join("\t");
            vcf.push_str(&format!("{chrom}\t{pos}\tvar\tA\tC\t.\t.\t.\tGT\t{gts}\n"));
        }
        vcf.into_bytes()
    }

    /// The five variants the two populations of the tests below are worked
    /// out on, with the major allele frequency each of them has in each
    /// population and over all four individuals:
    ///
    /// | position | `pop_a` | `pop_b` | all four |
    /// |---|---|---|---|
    /// | 100 | 1 | 0.5 | 0.75 |
    /// | 200 | 0.75 | 1 | 0.875 |
    /// | 300 | none | 0.5 | 0.5 |
    /// | 400 | 0.75 | 0.5 | 0.625 |
    /// | 500 | 1 | 1 | 0.5 |
    ///
    /// Neither individual of `pop_a` was called at the variant at 300, so
    /// it has no frequency there and is left out of that population
    /// whatever the threshold, and it is in `pop_b`. The major allele of
    /// the variant at 500 is the 1 in `pop_a` and the 0 in `pop_b`, and
    /// over all four individuals the two alleles are called four times
    /// each, so the dosages of it are not the same in the two populations
    /// and are neither of them the dosages of all four individuals.
    const THE_VARIANTS_OF_THE_TWO_POPS: [(&str, u64, [&str; 4]); 5] = [
        ("chr1", 100, ["0/0", "0/0", "0/1", "0/1"]),
        ("chr1", 200, ["0/0", "0/1", "0/0", "0/0"]),
        ("chr1", 300, ["./.", "./.", "0/1", "0/1"]),
        ("chr1", 400, ["0/0", "0/1", "0/1", "0/1"]),
        ("chr1", 500, ["1/1", "1/1", "0/0", "0/0"]),
    ];

    /// The individuals of `pop_a`, which are `i0` and `i1`.
    const POP_A: [usize; 2] = [0, 1];

    /// The individuals of `pop_b`, which are `i2` and `i3`.
    const POP_B: [usize; 2] = [2, 3];

    /// The dosages of the populations after a pass over the bytes of a VCF
    /// read in blocks of `num_vars_per_block` variants.
    fn the_dosages_of_a_pass(
        vcf: &[u8],
        pops: &[&[usize]],
        num_vars_per_block: usize,
        max_dist: u64,
        max_allowed_maf: f64,
    ) -> TheDosagesOfThePops {
        let options = VcfOptions {
            num_vars_per_block: Some(num_vars_per_block),
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::new(Cursor::new(vcf.to_vec()), options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        };
        reader.set_needs(Needs::GTS | Needs::CHROM_POS);
        let mut of_the_pops = match TheDosagesOfThePops::of(pops, max_dist, max_allowed_maf) {
            Ok(of_the_pops) => of_the_pops,
            Err(error) => panic!("the populations were refused: {error}"),
        };
        loop {
            let block = match reader.next_block() {
                Ok(Some(block)) => block,
                Ok(None) => break,
                Err(error) => panic!("the block was not read: {error}"),
            };
            if let Err(error) = of_the_pops.take_the_block(block) {
                panic!("the block was refused: {error}");
            }
        }
        of_the_pops
    }

    /// The population at that position among the ones the pass was given.
    fn pop_of(of_the_pops: &TheDosagesOfThePops, pop: usize) -> &ThePopOverTheWindow {
        match of_the_pops.pop(pop) {
            Some(of_it) => of_it,
            None => panic!("the pass has no population {pop}"),
        }
    }

    /// The position of each variant a population holds.
    fn the_positions_of(pop: &ThePopOverTheWindow) -> Vec<u64> {
        pop.variants().iter().map(|variant| variant.pos).collect()
    }

    /// How many variants and how many individuals the dosages a population
    /// holds are of.
    fn the_size_of(pop: &ThePopOverTheWindow) -> (usize, usize) {
        match pop.dosages() {
            Some(dosages) => (dosages.num_vars(), dosages.num_individuals()),
            None => panic!("the population holds no dosages"),
        }
    }

    /// The dosage of each individual of a population at the variant `var`
    /// of the ones it holds.
    fn the_dosages_at(pop: &ThePopOverTheWindow, var: usize) -> Vec<Option<u8>> {
        let dosages = match pop.dosages() {
            Some(dosages) => dosages,
            None => panic!("the population holds no dosages"),
        };
        match dosages.dosages(var) {
            Some(of_the_variant) => of_the_variant.collect(),
            None => panic!("the dosages have no variant {var}"),
        }
    }

    #[test]
    fn two_pops_of_one_pass_keep_the_variants_of_their_own_major_allele_frequency() {
        let vcf = vcf_of_four_individuals(&THE_VARIANTS_OF_THE_TWO_POPS);
        let of_the_pops = the_dosages_of_a_pass(&vcf, &[&POP_A, &POP_B], 5, 1000, 0.75);
        assert_eq!(of_the_pops.num_pops(), 2);
        let pop_a = pop_of(&of_the_pops, 0);
        let pop_b = pop_of(&of_the_pops, 1);
        // `pop_a` keeps the two variants whose frequency in it is 0.75,
        // and `pop_b` the three whose frequency in it is 0.5. Only the
        // variant at 400 is in both.
        assert_eq!(the_positions_of(pop_a), vec![200, 400]);
        assert_eq!(the_positions_of(pop_b), vec![100, 300, 400]);
        assert_eq!((pop_a.num_vars(), pop_b.num_vars()), (2, 3));
        assert_eq!((the_size_of(pop_a), the_size_of(pop_b)), ((2, 2), (3, 2)));
    }

    #[test]
    fn a_variant_called_in_no_individual_of_a_pop_is_left_out_of_it_and_kept_in_another() {
        let vcf = vcf_of_four_individuals(&THE_VARIANTS_OF_THE_TWO_POPS);
        // At a threshold of 1 every variant that has a frequency passes,
        // so what leaves the variant at 300 out of `pop_a` is that neither
        // of its two individuals was called there.
        let of_the_pops = the_dosages_of_a_pass(&vcf, &[&POP_A, &POP_B], 5, 1000, 1.0);
        assert_eq!(
            the_positions_of(pop_of(&of_the_pops, 0)),
            vec![100, 200, 400, 500]
        );
        assert_eq!(
            the_positions_of(pop_of(&of_the_pops, 1)),
            vec![100, 200, 300, 400, 500]
        );
    }

    #[test]
    fn a_pop_in_which_every_variant_is_left_out_holds_nothing_and_the_other_is_not_affected() {
        let vcf = vcf_of_four_individuals(&THE_VARIANTS_OF_THE_TWO_POPS);
        // The frequencies of `pop_a` are 1, 0.75, none, 0.75 and 1, every
        // one of them above 0.6 or absent, and three of those of `pop_b`
        // are 0.5.
        let of_the_pops = the_dosages_of_a_pass(&vcf, &[&POP_A, &POP_B], 5, 1000, 0.6);
        let pop_a = pop_of(&of_the_pops, 0);
        assert_eq!(the_positions_of(pop_a), Vec::<u64>::new());
        assert_eq!(pop_a.num_vars(), 0);
        assert_eq!(the_size_of(pop_a), (0, 2));
        let pop_b = pop_of(&of_the_pops, 1);
        assert_eq!(the_positions_of(pop_b), vec![100, 300, 400]);
        assert_eq!(pop_b.num_vars(), 3);
    }

    #[test]
    fn the_variants_each_pop_keeps_are_the_same_at_every_size_of_block() {
        let vcf = vcf_of_four_individuals(&THE_VARIANTS_OF_THE_TWO_POPS);
        let what_they_keep = [
            (
                0.75,
                vec![200_u64, 400],
                2_u64,
                vec![100_u64, 300, 400],
                3_u64,
            ),
            (
                1.0,
                vec![100, 200, 400, 500],
                4,
                vec![100, 200, 300, 400, 500],
                5,
            ),
        ];
        for num_vars_per_block in [1, 2, 3, 5, 64] {
            for (max_allowed_maf, of_pop_a, num_vars_of_pop_a, of_pop_b, num_vars_of_pop_b) in
                &what_they_keep
            {
                let of_the_pops = the_dosages_of_a_pass(
                    &vcf,
                    &[&POP_A, &POP_B],
                    num_vars_per_block,
                    1000,
                    *max_allowed_maf,
                );
                let at = format!(
                    "at blocks of {num_vars_per_block} variants and a `max_allowed_maf` of {max_allowed_maf}"
                );
                let pop_a = pop_of(&of_the_pops, 0);
                let pop_b = pop_of(&of_the_pops, 1);
                assert_eq!(&the_positions_of(pop_a), of_pop_a, "{at}");
                assert_eq!(&the_positions_of(pop_b), of_pop_b, "{at}");
                assert_eq!(
                    (pop_a.num_vars(), pop_b.num_vars()),
                    (*num_vars_of_pop_a, *num_vars_of_pop_b),
                    "{at}"
                );
            }
        }
    }

    #[test]
    fn the_dosages_held_are_those_of_the_variants_kept_over_the_individuals_of_the_pop() {
        let vcf = vcf_of_four_individuals(&THE_VARIANTS_OF_THE_TWO_POPS);
        let of_the_pops = the_dosages_of_a_pass(&vcf, &[&POP_A, &POP_B], 5, 1000, 1.0);
        let pop_a = pop_of(&of_the_pops, 0);
        assert_eq!(the_size_of(pop_a), (4, 2));
        // The four variants `pop_a` kept are those at 100, 200, 400 and
        // 500. The first is `0/0 0/0`, the second and the third are
        // `0/0 0/1`, and the fourth is `1/1 1/1`, whose major allele in
        // `pop_a` is the 1: over all four individuals the 0 and the 1 are
        // called four times each and the major allele is the 0, which
        // would make the dosages of that variant 2 and 2.
        assert_eq!(the_dosages_at(pop_a, 0), vec![Some(0), Some(0)]);
        assert_eq!(the_dosages_at(pop_a, 1), vec![Some(0), Some(1)]);
        assert_eq!(the_dosages_at(pop_a, 2), vec![Some(0), Some(1)]);
        assert_eq!(the_dosages_at(pop_a, 3), vec![Some(0), Some(0)]);
        let pop_b = pop_of(&of_the_pops, 1);
        assert_eq!(the_size_of(pop_b), (5, 2));
        // `pop_b` kept the variant at 300 as well, so that one is its
        // third, and its two individuals are `0/1` there.
        assert_eq!(the_dosages_at(pop_b, 2), vec![Some(1), Some(1)]);
        // And the variant at 500 is its fifth, where both are `0/0` and
        // the major allele of `pop_b` is the 0.
        assert_eq!(the_dosages_at(pop_b, 4), vec![Some(0), Some(0)]);
    }

    #[test]
    fn a_variant_whose_frequency_is_exactly_the_threshold_is_kept_and_one_a_place_above_is_not() {
        let vcf = vcf_of_four_individuals(&THE_VARIANTS_OF_THE_TWO_POPS);
        let at_the_threshold = the_dosages_of_a_pass(&vcf, &[&POP_A], 5, 1000, 0.75);
        assert_eq!(
            the_positions_of(pop_of(&at_the_threshold, 0)),
            vec![200, 400]
        );
        // The frequency of those two in `pop_a` is three called alleles of
        // four, which is 0.75 to the bit, so a threshold one place below
        // it leaves that population with nothing.
        let below_it = the_dosages_of_a_pass(&vcf, &[&POP_A], 5, 1000, 0.749_999_999_999_999_9);
        assert_eq!(the_positions_of(pop_of(&below_it, 0)), Vec::<u64>::new());
    }

    #[test]
    fn the_variants_of_a_pop_fall_out_of_it_with_the_blocks_of_the_window() {
        let vcf = vcf_of_four_individuals(&THE_VARIANTS_OF_THE_TWO_POPS);
        // Blocks of one variant and a window of 100 bp: the newest variant
        // read is the one at 500, and the block of the variant at 300 fell
        // out of the window with it, so it is still held. The blocks that
        // fall out with a block are dropped when the next one is taken,
        // which `take_the_block` says why.
        let of_the_pops = the_dosages_of_a_pass(&vcf, &[&POP_A, &POP_B], 1, 100, 1.0);
        let pop_a = pop_of(&of_the_pops, 0);
        assert_eq!(the_positions_of(pop_a), vec![400, 500]);
        assert_eq!(the_size_of(pop_a), (2, 2));
        // `pop_a` kept four variants over the pass, those at 100, 200, 400
        // and 500, and the two that have fallen out of it are the first
        // two: it kept nothing of the block of the variant at 300.
        assert_eq!((pop_a.num_vars(), pop_a.first_var()), (4, 2));
        let pop_b = pop_of(&of_the_pops, 1);
        assert_eq!(the_positions_of(pop_b), vec![300, 400, 500]);
        assert_eq!(the_size_of(pop_b), (3, 2));
        // `pop_b` kept five, the one at 300 among them, and the two that
        // have fallen out of it are the same two as `pop_a`'s.
        assert_eq!((pop_b.num_vars(), pop_b.first_var()), (5, 2));
    }

    #[test]
    fn a_pass_of_no_population_is_one_population_of_every_individual() {
        let vcf = vcf_of_four_individuals(&THE_VARIANTS_OF_THE_TWO_POPS);
        let of_the_pops = the_dosages_of_a_pass(&vcf, &[], 5, 1000, 0.75);
        assert_eq!(of_the_pops.num_pops(), 1);
        let of_them_all = pop_of(&of_the_pops, 0);
        // Over the four individuals the frequencies are 0.75, 0.875, 0.5,
        // 0.625 and 0.5, so the variant at 200 is the one left out, and
        // that is the set of neither population at this threshold.
        assert_eq!(the_positions_of(of_them_all), vec![100, 300, 400, 500]);
        assert_eq!(of_them_all.num_vars(), 4);
        assert_eq!(the_size_of(of_them_all), (4, 4));
    }

    #[test]
    fn a_population_that_names_no_individual_is_refused() {
        let error = match TheDosagesOfThePops::of(&[&POP_A, &[]], 1000, 0.75) {
            Ok(of_the_pops) => panic!(
                "a population of no individual was taken, of {} populations",
                of_the_pops.num_pops()
            ),
            Err(error) => error,
        };
        assert!(
            matches!(error, Error::LdPopWithNoIndividual { pop: 1 }),
            "{error}"
        );
    }

    #[test]
    fn a_max_allowed_maf_that_is_not_a_frequency_is_refused() {
        for value in [-0.1, 1.1, f64::NAN] {
            let error = match TheDosagesOfThePops::of(&[&POP_A], 1000, value) {
                Ok(of_the_pops) => panic!(
                    "a `max_allowed_maf` of {value} was taken, for {} populations",
                    of_the_pops.num_pops()
                ),
                Err(error) => error,
            };
            // The bits and not the value, so that the NaN of the three
            // matches the one that was given.
            assert!(
                matches!(
                    error,
                    Error::LdMaxAllowedMafOutOfRange { value: found }
                        if found.to_bits() == value.to_bits()
                ),
                "{error}"
            );
        }
    }

    /// A block of the variants at `poss` of `num_individuals` diploid
    /// individuals of the chromosome 0, whose genotypes are all the
    /// reference allele.
    fn the_block_of(num_individuals: usize, poss: Vec<u64>) -> Block {
        let num_vars = poss.len();
        Block {
            num_vars,
            num_individuals,
            ploidy: 2,
            gts: vec![0; num_vars.saturating_mul(num_individuals).saturating_mul(2)],
            chrom: Some(vec![0; num_vars]),
            pos: Some(poss),
            id: None,
            alleles: None,
            qual: None,
        }
    }

    #[test]
    fn a_block_of_other_individuals_than_the_first_one_of_the_pass_is_refused() {
        let mut of_the_pops = match TheDosagesOfThePops::of(&[], 1000, 1.0) {
            Ok(of_the_pops) => of_the_pops,
            Err(error) => panic!("the populations were refused: {error}"),
        };
        if let Err(error) = of_the_pops.take_the_block(the_block_of(2, vec![100])) {
            panic!("the first block was refused: {error}");
        }
        let error = match of_the_pops.take_the_block(the_block_of(3, vec![200])) {
            Ok(()) => panic!("a block of three individuals was taken after one of two"),
            Err(error) => error,
        };
        assert!(
            matches!(
                error,
                Error::BlocksDoNotFitTogether {
                    num_individuals: 2,
                    ploidy: 2,
                    found_num_individuals: 3,
                    found_ploidy: 2
                }
            ),
            "{error}"
        );
    }

    /// The bins of a pass over the bytes of a VCF read in blocks of
    /// `num_vars_per_block` variants, with the pairs of each step taken in
    /// tiles of `vars_per_tile` variants.
    fn the_bins_of_a_pass(
        vcf: &[u8],
        pops: &[&[usize]],
        num_vars_per_block: usize,
        options: &LdAndDistOptions,
        vars_per_tile: usize,
    ) -> LdAndDist {
        match the_bins_or_the_error(vcf, pops, num_vars_per_block, options, vars_per_tile) {
            Ok(of_the_pass) => of_the_pass,
            Err(error) => panic!("the pass was refused: {error}"),
        }
    }

    /// What a pass over the bytes of a VCF gave, or what it was refused
    /// with.
    fn the_bins_or_the_error(
        vcf: &[u8],
        pops: &[&[usize]],
        num_vars_per_block: usize,
        options: &LdAndDistOptions,
        vars_per_tile: usize,
    ) -> Result<LdAndDist, Error> {
        let vcf_options = VcfOptions {
            num_vars_per_block: Some(num_vars_per_block),
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::new(Cursor::new(vcf.to_vec()), vcf_options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        };
        the_ld_and_dist_in_tiles_of(&mut reader, pops, options, vars_per_tile)
    }

    /// The bins of the population at that position among the ones the pass
    /// was given.
    fn bins_of(of_the_pass: &LdAndDist, pop: usize) -> &LdBins {
        match of_the_pass.bins_of_pop(pop) {
            Some(bins) => bins,
            None => panic!("the pass has no population {pop}"),
        }
    }

    /// How many pairs each bin holds, the first bin first.
    fn the_pairs_of(bins: &LdBins) -> Vec<u64> {
        (0..bins.num_bins())
            .map(|bin| match bins.num_pairs(bin) {
                Some(num_pairs) => num_pairs,
                None => panic!("the bins have no bin {bin}"),
            })
            .collect()
    }

    /// The smallest and the largest distance of each bin, the first bin
    /// first.
    fn the_bounds_of(bins: &LdBins) -> Vec<(u64, u64)> {
        (0..bins.num_bins())
            .map(|bin| match bins.bounds(bin) {
                Some(bounds) => bounds,
                None => panic!("the bins have no bin {bin}"),
            })
            .collect()
    }

    /// The mean and the standard deviation of the r² of each bin, in the
    /// bits they have, and `None` for a bin with no pair.
    fn the_values_of(bins: &LdBins) -> Vec<Option<(u64, u64)>> {
        (0..bins.num_bins())
            .map(|bin| match (bins.mean_r2(bin), bins.sd_r2(bin)) {
                (Some(mean), Some(sd)) => Some((mean.to_bits(), sd.to_bits())),
                (None, None) => None,
                (mean, sd) => panic!("the bin {bin} has a mean of {mean:?} and an sd of {sd:?}"),
            })
            .collect()
    }

    /// Asserts that `found` and `expected` are within 1e-12 of each other,
    /// relative to `expected`, which is what "How it is verified" of
    /// `docs/specs/ld.md` compares the bins with.
    fn assert_the_value_is(found: f64, expected: f64, what: &str) {
        let apart = (found - expected).abs();
        let allowed = expected.abs() * 1e-12;
        assert!(
            apart <= allowed,
            "{what} is {found:?} and not {expected:?}, {apart:e} apart"
        );
    }

    /// The mean of the r² of the bin, which the caller expects to hold a
    /// pair.
    fn the_mean_of(bins: &LdBins, bin: usize) -> f64 {
        match bins.mean_r2(bin) {
            Some(mean) => mean,
            None => panic!("the bin {bin} has no mean"),
        }
    }

    /// The standard deviation of the r² of the bin, which the caller
    /// expects to hold a pair.
    fn the_sd_of(bins: &LdBins, bin: usize) -> f64 {
        match bins.sd_r2(bin) {
            Some(sd) => sd,
            None => panic!("the bin {bin} has no standard deviation"),
        }
    }

    /// A file of `tests/reference/ld/`, which lives at the root of the
    /// repository, beside the script that writes the files again, and not
    /// inside this crate.
    fn the_reference_path(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/ld")
            .join(name)
    }

    /// `tests/reference/ld/example.vcf`, the five variants of six diploid
    /// individuals of "How it is verified" of `docs/specs/ld.md`, at 1000,
    /// 2000, 3000, 4000 and 5000 base pairs of `chr1`.
    fn the_example_vcf() -> Vec<u8> {
        let path = the_reference_path("example.vcf");
        match std::fs::read(&path) {
            Ok(bytes) => bytes,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        }
    }

    /// The bins the worked example of "How it is verified" of
    /// `docs/specs/ld.md` is read with: the whole range of the distances
    /// of the five variants in two bins of 2000 base pairs.
    fn the_options_of_the_example() -> LdAndDistOptions {
        LdAndDistOptions {
            min_dist: 1,
            max_dist: 4000,
            num_bins: 2,
            // The variant v4 is the reference allele in all six
            // individuals, so its major allele frequency is 1 and a
            // smaller threshold would leave it out of the pass, where the
            // worked example has it in and every pair that holds it in no
            // bin.
            max_allowed_maf: 1.0,
        }
    }

    #[test]
    fn the_two_bins_of_the_worked_example_are_the_ones_of_the_spec() {
        let options = the_options_of_the_example();
        let of_the_pass = the_bins_of_a_pass(&the_example_vcf(), &[], 5, &options, 256);
        assert_eq!((of_the_pass.num_vars(), of_the_pass.num_pops()), (5, 1));
        let bins = bins_of(&of_the_pass, 0);
        assert_eq!(bins.num_vars(), 5);
        assert_eq!(the_bounds_of(bins), vec![(1, 2000), (2001, 4000)]);
        // The first bin holds v1-v2, v2-v3 and v3-v4 at 1000 base pairs
        // and v1-v3, v2-v4 and v3-v5 at 2000, and the second v1-v4 and
        // v2-v5 at 3000 and v1-v5 at 4000. v4 has no variance, so the four
        // pairs that hold it have no r² and are in no bin.
        assert_eq!(the_pairs_of(bins), vec![4, 2]);
        assert_the_value_is(
            the_mean_of(bins, 0),
            0.572_767_857_142_857_2,
            "the mean r² of the first bin",
        );
        assert_the_value_is(
            the_mean_of(bins, 1),
            0.031_25,
            "the mean r² of the second bin",
        );
        // The second bin holds v2-v5, whose r² is 0, and v1-v5, whose r²
        // is 0.0625, so it is half of 0.0625 from their mean either way.
        assert_the_value_is(
            the_sd_of(bins, 1),
            0.031_25,
            "the standard deviation of the second bin",
        );
        assert!(
            of_the_pass.bins_of_pop(1).is_none(),
            "the pass of one population answered for a second"
        );
    }

    #[test]
    fn the_bins_of_the_worked_example_do_not_move_with_the_blocks_nor_with_the_tiles() {
        let vcf = the_example_vcf();
        let options = the_options_of_the_example();
        let of_one_block = the_bins_of_a_pass(&vcf, &[], 5, &options, 256);
        let expected = bins_of(&of_one_block, 0);
        for num_vars_per_block in [1, 2, 3, 4, 5, 64] {
            for vars_per_tile in [1, 2, 3, 256] {
                let of_the_pass =
                    the_bins_of_a_pass(&vcf, &[], num_vars_per_block, &options, vars_per_tile);
                let found = bins_of(&of_the_pass, 0);
                let at = format!(
                    "at blocks of {num_vars_per_block} variants and tiles of {vars_per_tile}"
                );
                assert_eq!(the_pairs_of(found), the_pairs_of(expected), "{at}");
                assert_eq!(the_values_of(found), the_values_of(expected), "{at}");
            }
        }
    }

    /// A VCF of 300 variants of six diploid individuals: 200 of `chr1` a
    /// thousand base pairs apart and 100 of `chr2` the same, with a
    /// genotype of each individual drawn from a generator of its own so
    /// that the variants differ in their frequencies and in what is
    /// missing.
    ///
    /// It is what says that the bins do not move with the blocks: a window
    /// of 10000 base pairs holds eleven of these variants, so at blocks of
    /// 64 variants a block spans six windows and the blocks that fall out
    /// of the window with it hold the pairs the block before it makes.
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "300 positions of at most 300000 and a generator of whole numbers below \
                  1000, built here"
    )]
    fn the_vcf_of_a_long_pass() -> Vec<u8> {
        let mut vcf = String::from(
            "##fileformat=VCFv4.2\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti0\ti1\ti2\ti3\ti4\ti5\n",
        );
        // A linear congruential generator, which gives the same genotypes
        // on every machine and in every build.
        let mut state = 12_345_u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            state >> 40
        };
        for var in 0..300_u64 {
            let (chrom, at) = match var < 200 {
                true => ("chr1", var + 1),
                false => ("chr2", var - 199),
            };
            let gts: Vec<&str> = (0..6)
                .map(|_| match next() % 10 {
                    0 => "./.",
                    1..=3 => "0/0",
                    4..=7 => "0/1",
                    _ => "1/1",
                })
                .collect();
            let gts = gts.join("\t");
            let pos = at * 1000;
            vcf.push_str(&format!(
                "{chrom}\t{pos}\tv{var}\tA\tC\t.\t.\t.\tGT\t{gts}\n"
            ));
        }
        vcf.into_bytes()
    }

    #[test]
    fn the_bins_do_not_move_with_the_blocks_nor_with_the_tiles() {
        let vcf = the_vcf_of_a_long_pass();
        let pop_a: [usize; 3] = [0, 1, 2];
        let pop_b: [usize; 3] = [3, 4, 5];
        let pops: [&[usize]; 2] = [&pop_a, &pop_b];
        let options = LdAndDistOptions {
            min_dist: 1,
            max_dist: 10_000,
            num_bins: 5,
            max_allowed_maf: 0.9,
        };
        let of_one_block = the_bins_of_a_pass(&vcf, &pops, 500, &options, 256);
        // The window holds eleven variants of the 10000 base pairs it
        // reaches back, so the pass counts pairs over more than one block
        // at every size of block below that.
        let of_the_first = bins_of(&of_one_block, 0);
        assert!(
            the_pairs_of(of_the_first).iter().all(|pairs| *pairs > 0),
            "a bin of the first population holds no pair: {:?}",
            the_pairs_of(of_the_first)
        );
        for num_vars_per_block in [1, 2, 3, 7, 64, 500] {
            for vars_per_tile in [1, 2, 7, 256] {
                let of_the_pass =
                    the_bins_of_a_pass(&vcf, &pops, num_vars_per_block, &options, vars_per_tile);
                let at = format!(
                    "at blocks of {num_vars_per_block} variants and tiles of {vars_per_tile}"
                );
                assert_eq!(of_the_pass.num_vars(), of_one_block.num_vars(), "{at}");
                for pop in 0..2 {
                    let found = bins_of(&of_the_pass, pop);
                    let expected = bins_of(&of_one_block, pop);
                    assert_eq!(
                        found.num_vars(),
                        expected.num_vars(),
                        "{at}, population {pop}"
                    );
                    assert_eq!(
                        the_pairs_of(found),
                        the_pairs_of(expected),
                        "{at}, population {pop}"
                    );
                    assert_eq!(
                        the_values_of(found),
                        the_values_of(expected),
                        "{at}, population {pop}"
                    );
                }
            }
        }
    }

    #[test]
    fn each_pop_counts_the_pairs_of_the_variants_it_kept_of_its_own_frequency() {
        let vcf = vcf_of_four_individuals(&THE_VARIANTS_OF_THE_TWO_POPS);
        let options = LdAndDistOptions {
            min_dist: 1,
            max_dist: 500,
            num_bins: 1,
            max_allowed_maf: 0.75,
        };
        let of_the_pass = the_bins_of_a_pass(&vcf, &[&POP_A, &POP_B], 5, &options, 256);
        assert_eq!(of_the_pass.num_vars(), 5);
        // `pop_a` keeps the variants at 200 and at 400, whose dosages in
        // its two individuals are 0 and 1 at both, so their one pair has
        // an r² of 1.
        let pop_a = bins_of(&of_the_pass, 0);
        assert_eq!((pop_a.num_vars(), the_pairs_of(pop_a)), (2, vec![1]));
        assert_the_value_is(the_mean_of(pop_a, 0), 1.0, "the mean r² of `pop_a`");
        assert_the_value_is(the_sd_of(pop_a, 0), 0.0, "the sd of r² of `pop_a`");
        // `pop_b` keeps three variants, at 100, 300 and 400, and its two
        // individuals are `0/1` at each of them, so no variant of it has
        // variance and none of its three pairs has an r².
        let pop_b = bins_of(&of_the_pass, 1);
        assert_eq!((pop_b.num_vars(), the_pairs_of(pop_b)), (3, vec![0]));
        assert_eq!((pop_b.mean_r2(0), pop_b.sd_r2(0)), (None, None));
    }

    #[test]
    fn a_pop_in_which_every_variant_is_left_out_gives_no_variant_and_every_bin_empty() {
        let vcf = vcf_of_four_individuals(&THE_VARIANTS_OF_THE_TWO_POPS);
        let options = LdAndDistOptions {
            min_dist: 1,
            max_dist: 500,
            num_bins: 2,
            max_allowed_maf: 0.6,
        };
        let of_the_pass = the_bins_of_a_pass(&vcf, &[&POP_A, &POP_B], 2, &options, 256);
        // The frequencies of `pop_a` are 1, 0.75, none, 0.75 and 1, every
        // one of them above 0.6 or absent.
        let pop_a = bins_of(&of_the_pass, 0);
        assert_eq!((pop_a.num_vars(), the_pairs_of(pop_a)), (0, vec![0, 0]));
        assert_eq!(the_values_of(pop_a), vec![None, None]);
        // And the other population is not affected: it keeps the three
        // variants whose frequency in it is 0.5.
        assert_eq!(bins_of(&of_the_pass, 1).num_vars(), 3);
    }

    #[test]
    fn a_dataset_whose_variants_are_closer_than_min_dist_gives_every_bin_empty() {
        let vcf = vcf_of_four_individuals(&THE_VARIANTS_OF_THE_TWO_POPS);
        let options = LdAndDistOptions {
            min_dist: 1000,
            max_dist: 5000,
            num_bins: 3,
            max_allowed_maf: 1.0,
        };
        // The five variants span 400 base pairs, which is less than the
        // smallest distance a pair is counted at.
        let of_the_pass = the_bins_of_a_pass(&vcf, &[], 5, &options, 256);
        let bins = bins_of(&of_the_pass, 0);
        assert_eq!((bins.num_vars(), the_pairs_of(bins)), (5, vec![0, 0, 0]));
        assert_eq!(the_values_of(bins), vec![None, None, None]);
    }

    #[test]
    fn a_dataset_whose_variants_are_each_on_a_chromosome_of_their_own_gives_every_bin_empty() {
        let vcf = vcf_of_four_individuals(&[
            ("chr1", 100, ["0/0", "0/1", "0/1", "1/1"]),
            ("chr2", 200, ["0/0", "0/1", "1/1", "1/1"]),
            ("chr3", 300, ["0/0", "0/0", "0/1", "1/1"]),
        ]);
        let options = LdAndDistOptions {
            min_dist: 1,
            max_dist: 5000,
            num_bins: 2,
            max_allowed_maf: 1.0,
        };
        let of_the_pass = the_bins_of_a_pass(&vcf, &[], 1, &options, 256);
        let bins = bins_of(&of_the_pass, 0);
        assert_eq!((bins.num_vars(), the_pairs_of(bins)), (3, vec![0, 0]));
    }

    #[test]
    fn a_bin_of_one_pair_has_that_r2_as_its_mean_and_a_standard_deviation_of_zero() {
        let options = LdAndDistOptions {
            min_dist: 4000,
            max_dist: 4000,
            num_bins: 1,
            max_allowed_maf: 1.0,
        };
        // The one pair of the worked example at 4000 base pairs is v1-v5,
        // whose r² "How it is verified" of `docs/specs/ld.md` gives.
        let of_the_pass = the_bins_of_a_pass(&the_example_vcf(), &[], 5, &options, 256);
        let bins = bins_of(&of_the_pass, 0);
        assert_eq!(the_pairs_of(bins), vec![1]);
        assert_eq!(the_bounds_of(bins), vec![(4000, 4000)]);
        assert_the_value_is(the_mean_of(bins, 0), 0.0625, "the mean r² of the one bin");
        assert_the_value_is(the_sd_of(bins, 0), 0.0, "the sd of r² of the one bin");
    }

    #[test]
    fn the_bins_are_of_equal_width_across_the_distances_asked_for() {
        let options = LdAndDistOptions {
            min_dist: 1,
            max_dist: 250_000,
            num_bins: 10,
            max_allowed_maf: 1.0,
        };
        // The ten bins of the tables of "How it is verified" of
        // `docs/specs/ld.md`, 25000 base pairs each.
        let of_the_pass = the_bins_of_a_pass(&the_example_vcf(), &[], 5, &options, 256);
        let bins = bins_of(&of_the_pass, 0);
        assert_eq!(
            the_bounds_of(bins),
            vec![
                (1, 25_000),
                (25_001, 50_000),
                (50_001, 75_000),
                (75_001, 100_000),
                (100_001, 125_000),
                (125_001, 150_000),
                (150_001, 175_000),
                (175_001, 200_000),
                (200_001, 225_000),
                (225_001, 250_000),
            ]
        );
        assert!(
            bins.bounds(10).is_none(),
            "ten bins answered for an eleventh"
        );
    }

    #[test]
    fn a_pair_falls_in_the_bin_of_its_distance_when_the_bins_do_not_divide_the_distances() {
        // Ten distances in three bins is a width of 3.3333 base pairs, so
        // the bins hold the distances 1 to 4, 5 to 7 and 8 to 10.
        let vcf = vcf_of_four_individuals(&[
            ("chr1", 1, ["0/0", "0/0", "0/1", "1/1"]),
            ("chr1", 5, ["0/0", "0/1", "0/1", "1/1"]),
            ("chr1", 8, ["0/0", "0/0", "0/1", "0/1"]),
            ("chr1", 11, ["0/1", "1/1", "0/0", "0/1"]),
        ]);
        let options = LdAndDistOptions {
            min_dist: 1,
            max_dist: 10,
            num_bins: 3,
            max_allowed_maf: 1.0,
        };
        let of_the_pass = the_bins_of_a_pass(&vcf, &[], 4, &options, 256);
        let bins = bins_of(&of_the_pass, 0);
        assert_eq!(the_bounds_of(bins), vec![(1, 4), (5, 7), (8, 10)]);
        // The six pairs are at 4, 7 and 10 base pairs from the first
        // variant, at 3 and 6 from the second and at 3 from the third.
        assert_eq!(the_pairs_of(bins), vec![3, 2, 1]);
    }

    #[test]
    fn a_pair_at_the_smallest_distance_of_a_bin_falls_in_that_bin() {
        // Eighteen distances in fourteen bins is a width of 9/7 base
        // pairs, which no f64 holds: the bin of a distance and the
        // smallest distance of a bin are the arithmetic of "What it
        // gives" of `docs/specs/ld.md` in whole numbers, the bin of the
        // distance d being floor((d − min_dist) · num_bins / (max_dist −
        // min_dist + 1)) and the smallest distance of the bin b being
        // min_dist + ceil(b · (max_dist − min_dist + 1) / num_bins). The
        // distance 10 is 7 · 9/7 past the first bin's own, so it is the
        // smallest distance of the eighth bin and falls in it.
        let vcf = vcf_of_four_individuals(&[
            ("chr1", 1, ["0/0", "0/0", "0/1", "1/1"]),
            ("chr1", 11, ["0/0", "0/1", "0/1", "1/1"]),
            ("chr1", 19, ["0/0", "0/0", "0/1", "0/1"]),
        ]);
        let options = LdAndDistOptions {
            min_dist: 1,
            max_dist: 18,
            num_bins: 14,
            max_allowed_maf: 1.0,
        };
        let of_the_pass = the_bins_of_a_pass(&vcf, &[], 3, &options, 256);
        let bins = bins_of(&of_the_pass, 0);
        assert_eq!(
            the_bounds_of(bins),
            vec![
                (1, 2),
                (3, 3),
                (4, 4),
                (5, 6),
                (7, 7),
                (8, 8),
                (9, 9),
                (10, 11),
                (12, 12),
                (13, 13),
                (14, 15),
                (16, 16),
                (17, 17),
                (18, 18),
            ]
        );
        // The three pairs are at 10, 8 and 18 base pairs, which are the
        // eighth bin, the sixth and the last.
        assert_eq!(
            the_pairs_of(bins),
            vec![0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1]
        );
    }

    /// A VCF of four diploid individuals with a variant at each of
    /// `variants`, the chromosome and the position of each.
    ///
    /// The genotypes of a variant are the next of four patterns, each of
    /// which holds two dosages at least, so every variant has a variance
    /// and no pair of them is without an r².
    fn the_vcf_of_variants_at(variants: &[(&str, u64)]) -> Vec<u8> {
        const THE_PATTERNS: [[&str; 4]; 4] = [
            ["0/0", "0/1", "0/1", "1/1"],
            ["0/0", "0/0", "0/1", "1/1"],
            ["0/1", "1/1", "0/0", "0/1"],
            ["1/1", "0/1", "0/0", "0/0"],
        ];
        let of_each: Vec<(&str, u64, [&str; 4])> = variants
            .iter()
            .zip(THE_PATTERNS.iter().copied().cycle())
            .map(|((chrom, pos), pattern)| (*chrom, *pos, pattern))
            .collect();
        vcf_of_four_individuals(&of_each)
    }

    /// The 255 variants of `chr3` a base pair apart that the test below
    /// puts before `chr2`, `chr1` and `chr2` again, and the three that
    /// follow them.
    fn the_vcf_of_a_chromosome_that_comes_back(of_the_first_chromosome: u64) -> Vec<u8> {
        let mut variants: Vec<(&str, u64)> = (1..=of_the_first_chromosome)
            .map(|pos| ("chr3", pos))
            .collect();
        variants.push(("chr2", 10));
        variants.push(("chr1", 100));
        variants.push(("chr2", 50));
        the_vcf_of_variants_at(&variants)
    }

    #[test]
    fn a_column_of_a_tile_pairs_with_a_variant_the_first_column_of_that_tile_cannot_reach() {
        // Nothing of `docs/specs/block.md` or of the VCF reader asks the
        // chromosomes of a source to be grouped, and here chr2 comes back
        // after chr1: the variants are 255 of chr3, then chr2 at 10, chr1
        // at 100 and chr2 at 50. The variant of chr1 is the 257th, so a
        // tile of 256 columns starts on it, and it reaches no variant
        // before it while the column beside it, chr2 at 50, pairs with
        // chr2 at 10 at a distance of 40.
        //
        // The variants of chr3 are a base pair apart, so the pairs of that
        // chromosome within 100 base pairs are the sum over the distances
        // d of 1 to 100 of 255 − d, which is 20450, and the pair of chr2
        // makes 20451.
        let vcf = the_vcf_of_a_chromosome_that_comes_back(255);
        let options = LdAndDistOptions {
            min_dist: 1,
            max_dist: 100,
            num_bins: 1,
            max_allowed_maf: 1.0,
        };
        for num_vars_per_block in [64, 258] {
            for vars_per_tile in [3, 64, 256, 258] {
                let of_the_pass =
                    the_bins_of_a_pass(&vcf, &[], num_vars_per_block, &options, vars_per_tile);
                let bins = bins_of(&of_the_pass, 0);
                assert_eq!(
                    the_pairs_of(bins),
                    vec![20_451],
                    "at blocks of {num_vars_per_block} variants and tiles of {vars_per_tile}"
                );
            }
        }
    }

    /// The bins of a pass that the options or the populations were
    /// refused, which the caller expects to have been.
    fn the_error_of_a_pass(pops: &[&[usize]], options: &LdAndDistOptions) -> Error {
        let vcf = vcf_of_four_individuals(&THE_VARIANTS_OF_THE_TWO_POPS);
        match the_bins_or_the_error(&vcf, pops, 2, options, 256) {
            Ok(of_the_pass) => panic!(
                "the pass was taken and gave {} variants of {} populations",
                of_the_pass.num_vars(),
                of_the_pass.num_pops()
            ),
            Err(error) => error,
        }
    }

    /// The bins a pass is asked for when nothing of the call is what is
    /// being refused.
    fn the_options_of_a_pass() -> LdAndDistOptions {
        LdAndDistOptions {
            min_dist: 1,
            max_dist: 500,
            num_bins: 2,
            max_allowed_maf: 0.75,
        }
    }

    /// Asserts that each of `said` is in the message of the error.
    fn assert_the_message_says(error: &Error, said: &[&str]) {
        let message = error.to_string();
        for what in said {
            assert!(message.contains(what), "`{what}` is not in `{message}`");
        }
    }

    #[test]
    fn a_min_dist_above_max_dist_is_refused() {
        let options = LdAndDistOptions {
            min_dist: 4001,
            max_dist: 4000,
            ..the_options_of_a_pass()
        };
        let error = the_error_of_a_pass(&[&POP_A], &options);
        assert!(
            matches!(
                error,
                Error::LdMinDistAboveMaxDist {
                    min_dist: 4001,
                    max_dist: 4000
                }
            ),
            "{error}"
        );
        assert_the_message_says(&error, &["min_dist", "4001", "max_dist", "4000"]);
    }

    #[test]
    fn a_num_bins_of_zero_is_refused() {
        let options = LdAndDistOptions {
            num_bins: 0,
            ..the_options_of_a_pass()
        };
        let error = the_error_of_a_pass(&[&POP_A], &options);
        assert!(matches!(error, Error::LdNoBins), "{error}");
        assert_the_message_says(&error, &["num_bins", "0"]);
    }

    #[test]
    fn a_max_allowed_maf_that_is_not_a_frequency_is_refused_at_the_call() {
        for (value, said) in [(-0.1, "-0.1"), (1.1, "1.1"), (f64::NAN, "NaN")] {
            let options = LdAndDistOptions {
                max_allowed_maf: value,
                ..the_options_of_a_pass()
            };
            let error = the_error_of_a_pass(&[&POP_A], &options);
            // The bits and not the value, so that the NaN of the three
            // matches the one that was given.
            assert!(
                matches!(
                    error,
                    Error::LdMaxAllowedMafOutOfRange { value: found }
                        if found.to_bits() == value.to_bits()
                ),
                "{error}"
            );
            assert_the_message_says(&error, &["max_allowed_maf", said]);
        }
    }

    #[test]
    fn a_population_that_names_no_individual_is_refused_at_the_call() {
        let error = the_error_of_a_pass(&[&POP_A, &[]], &the_options_of_a_pass());
        assert!(
            matches!(error, Error::LdPopWithNoIndividual { pop: 1 }),
            "{error}"
        );
        assert_the_message_says(&error, &["population", "1"]);
    }

    #[test]
    fn an_individual_the_dataset_has_not_is_refused() {
        let error = the_error_of_a_pass(&[&[0, 4]], &the_options_of_a_pass());
        assert!(
            matches!(
                error,
                Error::LdIndividualNotInTheDataset {
                    individual: 4,
                    num_individuals: 4
                }
            ),
            "{error}"
        );
        assert_the_message_says(&error, &["individual", "4"]);
    }

    #[test]
    fn an_individual_asked_for_twice_is_refused() {
        let error = the_error_of_a_pass(&[&[0, 1, 0]], &the_options_of_a_pass());
        assert!(
            matches!(error, Error::LdIndividualAskedForTwice { individual: 0 }),
            "{error}"
        );
        assert_the_message_says(&error, &["individual", "0"]);
    }

    #[test]
    fn a_pass_with_no_variant_is_refused() {
        let vcf = vcf_of_four_individuals(&[]);
        let options = the_options_of_a_pass();
        let error = match the_bins_or_the_error(&vcf, &[&POP_A], 2, &options, 256) {
            Ok(of_the_pass) => panic!(
                "a pass of no variant gave {} variants",
                of_the_pass.num_vars()
            ),
            Err(error) => error,
        };
        assert!(
            matches!(
                error,
                Error::PassGaveNoVariant {
                    num_vars_of_the_source: 0,
                    ref filters
                } if filters.is_empty()
            ),
            "{error}"
        );
    }

    /// `tests/reference/ld/ld.vcf.gz`, the 500 variants of 100 diploid
    /// individuals `i000` to `i099` that "How it is verified" of
    /// `docs/specs/ld.md` counts its bins over, read as diploid and with
    /// the variants that failed their FILTER among them, which is what
    /// plink2 was given, in blocks of `num_vars_per_block` variants.
    fn the_ld_dataset(num_vars_per_block: usize) -> VcfReader<BufReader<File>> {
        let path = the_reference_path("ld.vcf.gz");
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block: Some(num_vars_per_block),
        };
        match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        }
    }

    /// The fifty individuals of `ld.vcf.gz` that start at `first`.
    const fn the_individuals_from(first: usize) -> [usize; 50] {
        let mut individuals = [0; 50];
        let mut at = 0;
        // `first` is 0 or 50 at the two calls below and `at` is under 50,
        // so neither sum passes 99, the last individual of the dataset.
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "fifty indices of at most 99, from the two calls below"
        )]
        while at < 50 {
            individuals[at] = first + at;
            at += 1;
        }
        individuals
    }

    /// The individuals of `pop_a`, `i000` to `i049` of `ld.vcf.gz`.
    const THE_INDIVIDUALS_OF_POP_A: [usize; 50] = the_individuals_from(0);

    /// The individuals of `pop_b`, `i050` to `i099` of `ld.vcf.gz`.
    const THE_INDIVIDUALS_OF_POP_B: [usize; 50] = the_individuals_from(50);

    /// The bins of "How it is verified" of `docs/specs/ld.md`: the
    /// distances from 1 to 250000 base pairs cut into ten of 25000, with
    /// the smallest and the largest distance of each.
    const THE_BOUNDS_OF_THE_TABLES: [(u64, u64); 10] = [
        (1, 25_000),
        (25_001, 50_000),
        (50_001, 75_000),
        (75_001, 100_000),
        (100_001, 125_000),
        (125_001, 150_000),
        (150_001, 175_000),
        (175_001, 200_000),
        (200_001, 225_000),
        (225_001, 250_000),
    ];

    /// The first table of "How it is verified" of `docs/specs/ld.md`, the
    /// one population of every one of the 100 individuals at a
    /// `max_allowed_maf` of 0.95: for each of the ten bins, how many
    /// pairs it holds, the mean of their r² and its standard deviation.
    ///
    /// Every value is the one that table prints, which
    /// `docs/reports/ld-method/bins.py` worked out from the r² plink2
    /// v2.0.0-a.7.7 gives for these individuals and these variants, and
    /// which `tests/reference/ld/ld.bins.txt` holds again.
    const THE_BINS_OF_EVERY_INDIVIDUAL: [(u64, f64, f64); 10] = [
        (8744, 0.207_678_855_518_440_31, 0.205_686_529_744_794_65),
        (7815, 0.078_903_591_760_625_11, 0.086_250_492_144_140_23),
        (6846, 0.035_084_413_027_517_11, 0.041_191_741_689_879_3),
        (5962, 0.020_569_175_806_260_69, 0.026_923_590_343_909_974),
        (5140, 0.015_026_451_851_395_499, 0.020_606_044_304_379_79),
        (4168, 0.011_542_104_404_978_385, 0.015_425_275_975_059_542),
        (3308, 0.011_454_382_158_199_49, 0.015_306_615_529_156_098),
        (2447, 0.012_095_572_873_545_887, 0.016_610_730_726_125_223),
        (1481, 0.015_365_254_218_410_632, 0.020_746_495_170_409_326),
        (530, 0.013_266_303_346_602_112, 0.017_768_417_874_071_147),
    ];

    /// The second table of "How it is verified" of `docs/specs/ld.md` for
    /// `pop_a`, the individuals `i000` to `i049` at a `max_allowed_maf`
    /// of 0.8, as [`THE_BINS_OF_EVERY_INDIVIDUAL`] holds the first.
    ///
    /// The counts of pairs and the means are the ones that table prints.
    /// It leaves the standard deviations out to stay readable, so those
    /// are read from `tests/reference/ld/ld.bins.txt`, which holds what
    /// `docs/reports/ld-method/bins.py` printed for all three tables and
    /// which `tests/reference/ld/run_plink2.sh` writes again and
    /// compares.
    const THE_BINS_OF_POP_A: [(u64, f64, f64); 10] = [
        (7394, 0.222_263_164_322_283_82, 0.219_333_423_597_465_4),
        (6564, 0.094_268_621_223_52, 0.103_494_028_800_708_32),
        (5648, 0.045_650_405_823_598_73, 0.056_055_780_798_436_86),
        (4918, 0.030_489_796_800_475_328, 0.040_343_306_810_475_03),
        (4304, 0.024_750_005_446_480_792, 0.032_375_847_680_808_055),
        (3540, 0.022_203_502_044_442_88, 0.029_984_493_841_579_307),
        (2872, 0.017_701_369_916_056_54, 0.023_881_355_618_397_9),
        (2137, 0.020_224_950_448_715_07, 0.028_823_934_904_688_44),
        (1240, 0.017_923_378_431_226_36, 0.022_823_288_248_614_05),
        (438, 0.020_745_833_685_396_994, 0.025_769_969_059_702_198),
    ];

    /// The third table of "How it is verified" of `docs/specs/ld.md` for
    /// `pop_b`, the individuals `i050` to `i099` at a `max_allowed_maf`
    /// of 0.8, whose values come from where those of
    /// [`THE_BINS_OF_POP_A`] do.
    const THE_BINS_OF_POP_B: [(u64, f64, f64); 10] = [
        (7625, 0.219_351_925_929_983_45, 0.213_752_142_407_470_95),
        (6779, 0.087_787_564_637_199_26, 0.098_879_976_941_602_94),
        (5968, 0.044_293_996_251_491_8, 0.054_863_163_257_117_08),
        (5189, 0.032_133_599_685_763_4, 0.041_469_004_992_669_65),
        (4473, 0.026_760_898_025_174_62, 0.036_536_689_258_127_475),
        (3567, 0.021_365_898_293_232_58, 0.030_124_766_387_101_015),
        (2823, 0.025_781_473_696_727_46, 0.033_454_346_665_276_094),
        (2086, 0.022_946_293_772_725_81, 0.030_081_908_306_131_513),
        (1275, 0.022_448_095_913_909_734, 0.028_237_540_899_418_476),
        (415, 0.016_086_351_215_632_733, 0.019_851_310_268_078_202),
    ];

    /// How many of the 500 variants of `ld.vcf.gz` each of the three
    /// tables keeps: 432 at the `max_allowed_maf` of 0.95 of the first,
    /// and 396 and 402 at the 0.8 of `pop_a` and of `pop_b`, worked out
    /// over the individuals of each population alone.
    const THE_VARS_OF_THE_TABLES: (u64, u64, u64) = (432, 396, 402);

    /// The bins the three tables are counted in, from 1 to 250000 base
    /// pairs in ten, with the major allele frequency a variant is kept at.
    fn the_options_of_the_tables(max_allowed_maf: f64) -> LdAndDistOptions {
        LdAndDistOptions {
            min_dist: 1,
            max_dist: 250_000,
            num_bins: 10,
            max_allowed_maf,
        }
    }

    /// The three tables of "How it is verified" of `docs/specs/ld.md`
    /// over `ld.vcf.gz` read in blocks of `num_vars_per_block` variants:
    /// the one population of every individual at a `max_allowed_maf` of
    /// 0.95, and `pop_a` and `pop_b` at 0.8.
    ///
    /// They are two passes because the major allele frequency is one
    /// threshold for the whole call, and the first table is at another
    /// than the two below it.
    fn the_three_tables_of(num_vars_per_block: usize) -> (LdAndDist, LdAndDist) {
        let mut reader = the_ld_dataset(num_vars_per_block);
        let of_every_individual =
            match calc_ld_and_dist(&mut reader, &[], &the_options_of_the_tables(0.95)) {
                Ok(of_the_pass) => of_the_pass,
                Err(error) => panic!("the pass of every individual was refused: {error}"),
            };
        let mut reader = the_ld_dataset(num_vars_per_block);
        let pops: [&[usize]; 2] = [&THE_INDIVIDUALS_OF_POP_A, &THE_INDIVIDUALS_OF_POP_B];
        let of_the_two_pops =
            match calc_ld_and_dist(&mut reader, &pops, &the_options_of_the_tables(0.8)) {
                Ok(of_the_pass) => of_the_pass,
                Err(error) => panic!("the pass of the two populations was refused: {error}"),
            };
        (of_every_individual, of_the_two_pops)
    }

    /// The table of one population in the bits it came out with, which is
    /// what two runs are compared by.
    #[derive(Debug, PartialEq, Eq)]
    struct TheValuesOfATable {
        /// The variants the population kept at its major allele frequency.
        num_vars: u64,
        /// How many pairs each bin holds, the first bin first.
        num_pairs: Vec<u64>,
        /// The mean and the standard deviation of each bin in their bits,
        /// and `None` for a bin with no pair, as `the_values_of` gives
        /// them.
        values: Vec<Option<(u64, u64)>>,
    }

    /// The three tables of one run in the bits they came out with.
    fn the_values_of_the_three_tables(
        of_the_run: &(LdAndDist, LdAndDist),
    ) -> Vec<TheValuesOfATable> {
        let (of_every_individual, of_the_two_pops) = of_the_run;
        [
            bins_of(of_every_individual, 0),
            bins_of(of_the_two_pops, 0),
            bins_of(of_the_two_pops, 1),
        ]
        .into_iter()
        .map(|bins| TheValuesOfATable {
            num_vars: bins.num_vars(),
            num_pairs: the_pairs_of(bins),
            values: the_values_of(bins),
        })
        .collect()
    }

    /// Asserts that the bins of one population are the ten rows of its
    /// table of "How it is verified" of `docs/specs/ld.md`: the variants
    /// it kept and the pairs of each bin exactly, and the mean and the
    /// standard deviation of each bin within the 1e-12 relative that item
    /// compares them with.
    fn assert_the_bins_are(
        bins: &LdBins,
        expected: &[(u64, f64, f64); 10],
        num_vars: u64,
        what: &str,
    ) {
        assert_eq!(bins.num_vars(), num_vars, "the variants {what} kept");
        assert_eq!(bins.num_bins(), expected.len(), "the bins of {what}");
        assert_eq!(
            the_bounds_of(bins),
            THE_BOUNDS_OF_THE_TABLES.to_vec(),
            "the distances of the bins of {what}"
        );
        let pairs: Vec<u64> = expected.iter().map(|(pairs, _, _)| *pairs).collect();
        assert_eq!(the_pairs_of(bins), pairs, "the pairs of the bins of {what}");
        for (bin, (_, mean, sd)) in expected.iter().enumerate() {
            assert_the_value_is(
                the_mean_of(bins, bin),
                *mean,
                &format!("the mean r² of the bin {bin} of {what}"),
            );
            assert_the_value_is(
                the_sd_of(bins, bin),
                *sd,
                &format!("the standard deviation of the bin {bin} of {what}"),
            );
        }
    }

    /// Asserts that the three tables of a run are the ones of the spec,
    /// and says of each which run it was.
    fn assert_the_three_tables_are_the_ones_of_the_spec(
        of_the_run: &(LdAndDist, LdAndDist),
        at: &str,
    ) {
        let (of_every_individual, of_pop_a, of_pop_b) = THE_VARS_OF_THE_TABLES;
        assert_eq!(
            (of_the_run.0.num_vars(), of_the_run.0.num_pops()),
            (500, 1),
            "the pass of every individual {at}"
        );
        assert_eq!(
            (of_the_run.1.num_vars(), of_the_run.1.num_pops()),
            (500, 2),
            "the pass of the two populations {at}"
        );
        assert_the_bins_are(
            bins_of(&of_the_run.0, 0),
            &THE_BINS_OF_EVERY_INDIVIDUAL,
            of_every_individual,
            &format!("every individual {at}"),
        );
        assert_the_bins_are(
            bins_of(&of_the_run.1, 0),
            &THE_BINS_OF_POP_A,
            of_pop_a,
            &format!("pop_a {at}"),
        );
        assert_the_bins_are(
            bins_of(&of_the_run.1, 1),
            &THE_BINS_OF_POP_B,
            of_pop_b,
            &format!("pop_b {at}"),
        );
    }

    /// The three tables of "How it is verified" of `docs/specs/ld.md`,
    /// over `tests/reference/ld/ld.vcf.gz` read with the VCF reader in
    /// blocks of 7, 64 and 500 variants, which give the same numbers to
    /// the bit.
    ///
    /// Every r² behind the tables is plink2's and the binning is the
    /// arithmetic of the spec: the reference script runs plink2 on the
    /// individuals and the variants each population keeps, and
    /// `docs/reports/ld-method/bins.py` puts its matrix into the bins.
    /// The three sizes of block are what say that the window of the pass
    /// neither keeps a variant further back than `max_dist` nor drops one
    /// that a later variant still pairs with: these variants are a
    /// thousand base pairs apart at the closest, so a window of 250000
    /// reaches over many more of them than a block of 7 holds, and a
    /// block of 500 holds the whole dataset.
    #[test]
    fn the_three_tables_of_bins_are_the_ones_plink2_gives_at_every_size_of_block() {
        let of_seven = the_three_tables_of(7);
        assert_the_three_tables_are_the_ones_of_the_spec(&of_seven, "at blocks of 7 variants");
        for num_vars_per_block in [64, 500] {
            let at = format!("at blocks of {num_vars_per_block} variants");
            let of_the_run = the_three_tables_of(num_vars_per_block);
            assert_the_three_tables_are_the_ones_of_the_spec(&of_the_run, &at);
            assert_eq!(
                the_values_of_the_three_tables(&of_the_run),
                the_values_of_the_three_tables(&of_seven),
                "{at}, against blocks of 7"
            );
        }
    }

    /// The three tables are the same, to the bit, on a pool of one thread
    /// and on one of four.
    ///
    /// Both the reading and the products run on the threads of the pool
    /// the caller is in: the VCF reader parses the lines of a batch on
    /// them, as `docs/specs/io_vcf.md` has it, and the products of r² run
    /// on them through faer when the `blas` feature is off, which is what
    /// `cargo test -p popnei --no-default-features` runs. The pools are
    /// built here and are not rayon's global one, which has one thread
    /// per core of the machine, and `current_num_threads` inside the pool
    /// says how many threads the pass had. rayon is a dependency of the
    /// targets that are not wasm, so this test is compiled for those
    /// alone.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_number_of_threads_does_not_change_the_three_tables() {
        let in_a_pool = |threads: usize| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            pool.install(|| {
                assert_eq!(
                    rayon::current_num_threads(),
                    threads,
                    "the pass did not run on the pool it was given"
                );
                the_three_tables_of(64)
            })
        };

        let on_one = in_a_pool(1);
        assert_the_three_tables_are_the_ones_of_the_spec(&on_one, "on one thread");
        assert_eq!(
            the_values_of_the_three_tables(&in_a_pool(4)),
            the_values_of_the_three_tables(&on_one),
            "four threads against one"
        );
    }
}
