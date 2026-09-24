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
        reason = "nothing outside the tests holds a window yet: the dosages of each \
                  population over it and the tiles that count its pairs, which \
                  `docs/specs/ld.md` asks for, are written on top of this and are not \
                  written yet"
    )
)]

use std::collections::VecDeque;

use crate::block::Block;
use crate::error::{Error, Result};
use crate::variant::Needs;

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

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::{TheVariantOfTheWindow, TheWindowOfTheBlocks, the_variants_of};

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
}
