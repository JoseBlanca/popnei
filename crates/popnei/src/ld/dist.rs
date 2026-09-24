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
#![cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "nothing outside the tests reads a window or the dosages over it yet: \
                  the tiles that count the pairs of a window and `calc_ld_and_dist`, \
                  which `docs/specs/ld.md` asks for, are written on top of them and are \
                  not written yet"
    )
)]

use std::collections::VecDeque;

use crate::block::Block;
use crate::error::{Error, Result};
use crate::variant::Needs;

use super::{LdDosages, the_copy_of, the_memory_for};

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
    fn num_blocks(&self) -> usize {
        self.held.len()
    }

    /// How many variants the blocks it holds have, which is how many
    /// [`TheWindowOfTheBlocks::variants`] gives.
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
    fn blocks(&self) -> impl ExactSizeIterator<Item = &Block> {
        self.held.iter()
    }

    /// Where each variant of the blocks it holds lies, the oldest block
    /// first and inside a block in the order of its variants.
    fn variants(&self) -> impl Iterator<Item = TheVariantOfTheWindow> {
        self.held.iter().flat_map(the_variants_of)
    }

    /// The newest variant read, which the window reaches back from, and
    /// `None` before a block with a variant has been taken.
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
            pops: of_them,
        })
    }

    /// Takes the next block of the reader: each population keeps the
    /// variants of it that passed its major allele frequency, the blocks
    /// that fell out of the window take their variants with them, and the
    /// dosages of each population are built again over what is left.
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
        // The variants of the block are taken before the window is told of
        // it, so that what each population holds is in the order of the
        // blocks of the window and the blocks that fall out are the oldest
        // of both.
        for pop in &mut self.pops {
            pop.take_the_block(&block, of_the_source, self.max_allowed_maf)?;
        }
        let dropped = self.window.take_the_block(block)?;
        for pop in &mut self.pops {
            pop.the_oldest_blocks_are_dropped(dropped, of_the_source);
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
    fn pop(&self, pop: usize) -> Option<&ThePopOverTheWindow> {
        self.pops.get(pop)
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
}

/// A count of variants as a result carries it, which is a `u64` and not a
/// `usize` so that it is the same number in WebAssembly, where a `usize` is
/// 32 bits.
fn the_count_of(num_vars: usize) -> u64 {
    u64::try_from(num_vars).unwrap_or(u64::MAX)
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::{
        TheDosagesOfThePops, ThePopOverTheWindow, TheVariantOfTheWindow, TheWindowOfTheBlocks,
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
        // read is the one at 500, so only the blocks of the variants at
        // 400 and at 500 are still held.
        let of_the_pops = the_dosages_of_a_pass(&vcf, &[&POP_A, &POP_B], 1, 100, 1.0);
        let pop_a = pop_of(&of_the_pops, 0);
        assert_eq!(the_positions_of(pop_a), vec![400, 500]);
        assert_eq!(the_size_of(pop_a), (2, 2));
        // `pop_a` kept four variants over the pass and the two that have
        // fallen out of it are those at 100 and at 200.
        assert_eq!((pop_a.num_vars(), pop_a.first_var()), (4, 2));
        let pop_b = pop_of(&of_the_pops, 1);
        assert_eq!(the_positions_of(pop_b), vec![400, 500]);
        // `pop_b` kept five, the one at 300 among them, so three of its
        // own have fallen out where two of `pop_a`'s did.
        assert_eq!((pop_b.num_vars(), pop_b.first_var()), (5, 3));
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
}
