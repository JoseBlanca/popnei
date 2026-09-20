//! Blocks of variants: a run of consecutive variants held as arrays, with
//! the genotypes of all of them in one, and the collector that builds them
//! from any reader.
//!
//! The calculations that want matrices, the PCA and the kinship, consume
//! blocks, and a block is also the only way genotypes leave the core: a
//! Python or a TypeScript user asks for them block by block. A
//! [`BlockCollector`] reads variants one at a time into a [`Variant`] of its
//! own and copies each one into the arrays of the block it is building,
//! until the block holds the number of variants that was asked for or the
//! reader has no more.
//!
//! `docs/specs/block.md` has the design and section 2 of
//! `docs/architecture.md` the reasons for the arrays.

use std::fmt;

use crate::error::{Error, Result};
use crate::variant::{Needs, Variant, VariantReader};

/// How many genotypes a block holds when the caller asks for no number of
/// variants, 5 million, which is `DEF_NUM_GTS_PER_CHUNK` of pyNei's
/// `config.py`. At the ploidy 2 they are the 10 MB of one allocation.
/// Nobody has measured it for popnei.
pub const GENOTYPES_PER_BLOCK: usize = 5_000_000;

/// The fewest variants a block holds when the caller asks for no number,
/// which is what decides for a dataset of many individuals: 100, pyNei's
/// `MIN_NUM_VARS_PER_CHUNK`, measured for popnei by nobody.
pub const MIN_NUM_VARS_PER_BLOCK: usize = 100;

/// The most variants a block holds when the caller asks for no number,
/// which is what decides for a dataset of few individuals: 10000, pyNei's
/// `MAX_NUM_VARS_PER_CHUNK`, measured for popnei by nobody.
pub const MAX_NUM_VARS_PER_BLOCK: usize = 10_000;

/// How many variants a block holds for that many individuals when the
/// caller asks for no number: [`GENOTYPES_PER_BLOCK`] divided by the
/// individuals, never below [`MIN_NUM_VARS_PER_BLOCK`] and never above
/// [`MAX_NUM_VARS_PER_BLOCK`].
///
/// It is pyNei's `calc_num_vars_per_chunk` of `pynei/variants.py`: a block
/// is sized by the genotypes it holds and not by its variants, because
/// that is what the memory and the work depend on.
#[must_use]
pub fn default_num_vars_per_block(num_individuals: usize) -> usize {
    // No individual at all asks for every variant a block takes, as the
    // division by one individual would.
    let num_vars = GENOTYPES_PER_BLOCK
        .checked_div(num_individuals)
        .unwrap_or(MAX_NUM_VARS_PER_BLOCK);
    num_vars.clamp(MIN_NUM_VARS_PER_BLOCK, MAX_NUM_VARS_PER_BLOCK)
}

/// The alleles of the variants of one block, the reference allele of each
/// variant first and then its alternative ones, as the text the source
/// gave.
///
/// The texts are one buffer with the end of each allele in it, and not a
/// string for each allele: a block of 10000 variants is 10000 allocations
/// either way, and one buffer makes it one.
#[derive(Debug)]
pub struct AllelesColumn {
    /// The text of every allele of the block, one after another.
    texts: String,
    /// Where each allele ends in `texts`, one number for each allele.
    allele_ends: Vec<usize>,
    /// Where the alleles of each variant end in `allele_ends`, one number
    /// for each variant.
    var_ends: Vec<usize>,
}

impl AllelesColumn {
    /// A column with no variant in it, whose buffer of ends is reserved
    /// for `num_vars` of them.
    fn with_num_vars(num_vars: usize) -> AllelesColumn {
        AllelesColumn {
            texts: String::new(),
            allele_ends: Vec::new(),
            var_ends: Vec::with_capacity(num_vars),
        }
    }

    /// The alleles of one more variant, at the end of the column.
    fn push(&mut self, alleles: &[String]) {
        for allele in alleles {
            self.texts.push_str(allele);
            self.allele_ends.push(self.texts.len());
        }
        self.var_ends.push(self.allele_ends.len());
    }

    /// Where the alleles of `var` are in `allele_ends`, the first and the
    /// one after the last. Both are 0 for a variant that is not there.
    fn alleles_of(&self, var: usize) -> (usize, usize) {
        let Some(end) = self.var_ends.get(var).copied() else {
            return (0, 0);
        };
        let start = var
            .checked_sub(1)
            .and_then(|before| self.var_ends.get(before).copied())
            .unwrap_or(0);
        (start, end)
    }

    /// How many variants the column holds.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.var_ends.len()
    }

    /// How many alleles the variant `var` of the block has, the reference
    /// allele among them. 0 for a variant the column does not hold.
    #[must_use]
    pub fn num_alleles(&self, var: usize) -> usize {
        let (start, end) = self.alleles_of(var);
        end.saturating_sub(start)
    }

    /// The text of the allele `allele` of the variant `var`, as the source
    /// gave it: `A`, `<DEL>`, `*`. The allele 0 is the reference one.
    ///
    /// A variant or an allele the column does not hold gives an empty
    /// text, which no allele of a source is.
    #[must_use]
    pub fn allele(&self, var: usize, allele: usize) -> &str {
        let (first, end) = self.alleles_of(var);
        let Some(index) = first.checked_add(allele).filter(|index| *index < end) else {
            return "";
        };
        let Some(text_end) = self.allele_ends.get(index).copied() else {
            return "";
        };
        let text_start = index
            .checked_sub(1)
            .and_then(|before| self.allele_ends.get(before).copied())
            .unwrap_or(0);
        self.texts.get(text_start..text_end).unwrap_or("")
    }
}

/// A run of consecutive variants of one source, held as arrays.
///
/// A column other than the genotypes is there only when the collector was
/// asked for it, and `None` when it was not. The number of a chromosome is
/// a number of the table of the reader the block came from, and that table
/// grows while the source is read, so the name of a number is looked up
/// after the block was collected.
#[derive(Debug)]
pub struct Block {
    /// How many variants the block holds, which is the number that was
    /// asked for except in the last block of a source.
    pub num_vars: usize,
    /// How many individuals the source has, the same for every variant.
    pub num_individuals: usize,
    /// How many alleles the genotype of one individual holds.
    pub ploidy: usize,
    /// `num_vars` x `num_individuals` x `ploidy` alleles, variant after
    /// variant and inside a variant individual after individual. 0 is the
    /// reference allele, 1 and above the alternative ones, and
    /// [`MISSING_ALLELE`](crate::variant::MISSING_ALLELE), -1, an allele
    /// that was not called.
    pub gts: Vec<i8>,
    /// The number of the chromosome of each variant, in the
    /// [`ChromTable`](crate::variant::ChromTable) of the reader the block
    /// came from.
    pub chrom: Option<Vec<u32>>,
    /// The position of each variant, 1 based as in a VCF.
    pub pos: Option<Vec<u64>>,
    /// The id of each variant, empty for a variant that has none.
    pub id: Option<Vec<String>>,
    /// The alleles of each variant, the reference one first.
    pub alleles: Option<AllelesColumn>,
    /// The quality of each variant, phred scaled as the QUAL of a VCF, and
    /// NaN for a variant that has none.
    pub qual: Option<Vec<f32>>,
}

/// It builds blocks from a reader, one after another, each of the number of
/// variants that was asked for.
///
/// It owns its reader and gives it back to whoever wants the individuals,
/// the ploidy or the table of the chromosomes. What it keeps from one block
/// to the next is the [`Variant`] it lends to the reader; each block is
/// allocated when it is started and given away when it is full.
pub struct BlockCollector<R: VariantReader> {
    reader: R,
    /// What each block holds, the genotypes among them.
    needs: Needs,
    num_vars_per_block: usize,
    num_individuals: usize,
    ploidy: usize,
    /// `num_vars_per_block` x `num_individuals` x `ploidy`, the genotypes
    /// a full block holds, which every block is allocated for.
    gts_per_block: usize,
    /// The variant the reader is lent, again and again.
    var: Variant,
    /// Whether the reader has no more variants or gave an error. After
    /// either there is no block.
    finished: bool,
}

impl<R: VariantReader> BlockCollector<R> {
    /// The collector over `reader`, which it asks for `needs` and the
    /// genotypes.
    ///
    /// `needs` is what each block will hold besides the genotypes, which
    /// are always part of it, and a field that is not in it is a column
    /// that the block does not have and that the reader is never asked to
    /// parse. `num_vars_per_block` is how many variants a block holds, 1 or
    /// more, and `None` is [`default_num_vars_per_block`] for the
    /// individuals of the reader.
    ///
    /// # Errors
    ///
    /// When `num_vars_per_block` is 0, and when the genotypes of one block,
    /// the variants times the individuals times the ploidy, are more than
    /// the machine addresses.
    pub fn new(mut reader: R, needs: Needs, num_vars_per_block: Option<usize>) -> Result<Self> {
        let num_individuals = reader.individuals().len();
        let ploidy = reader.ploidy();
        let num_vars_per_block = match num_vars_per_block {
            Some(0) => return Err(Error::BlockOfNoVariants),
            Some(asked_for) => asked_for,
            None => default_num_vars_per_block(num_individuals),
        };
        let gts_per_block = num_individuals
            .checked_mul(ploidy)
            .and_then(|per_variant| per_variant.checked_mul(num_vars_per_block))
            .ok_or(Error::BlockTooLarge {
                num_vars_per_block,
                num_individuals,
                ploidy,
            })?;
        let needs = needs.union(Needs::GTS);
        reader.set_needs(needs);
        Ok(BlockCollector {
            reader,
            needs,
            num_vars_per_block,
            num_individuals,
            ploidy,
            gts_per_block,
            var: Variant::new(),
            finished: false,
        })
    }

    /// The next block, or `None` when the reader has no more variants.
    ///
    /// The last block of a source is the only one that can hold fewer
    /// variants than were asked for, and a reader with no variants gives no
    /// block at all.
    ///
    /// # Errors
    ///
    /// When the reader fails, and when it did not fill a field that was
    /// asked for. The block that was being built is lost with the error,
    /// and every call after it gives no block.
    pub fn next_block(&mut self) -> Result<Option<Block>> {
        if self.finished {
            return Ok(None);
        }
        // The columns are allocated at the first variant, so that the call
        // that finds the source at its end allocates nothing. The count of
        // the variants is the one of the range, so that a block of
        // `usize::MAX` variants needs no addition of our own.
        let mut block: Option<Block> = None;
        for count in 1..=self.num_vars_per_block {
            // The trait says that a reader ends at its error, and the
            // collector ends too: what a reader that goes on gives after
            // one is not the source, and no caller of ours reads it.
            let read = self.reader.read_variant(&mut self.var).inspect_err(|_| {
                self.finished = true;
            })?;
            if !read {
                self.finished = true;
                break;
            }
            let not_filled = self.needs.difference(self.var.filled);
            if !not_filled.is_empty() {
                self.finished = true;
                return Err(Error::FieldsNotFilled { fields: not_filled });
            }
            let being_built = block.get_or_insert_with(|| self.start_block());
            self.push_variant(being_built);
            being_built.num_vars = count;
        }
        Ok(block)
    }

    /// The reader the collector was built over, for the individuals, the
    /// ploidy and the names of the chromosome numbers of a block.
    pub fn reader(&self) -> &R {
        &self.reader
    }

    /// An empty block with every column the collector was asked for, each
    /// allocated for a full block.
    fn start_block(&self) -> Block {
        let num_vars = self.num_vars_per_block;
        Block {
            num_vars: 0,
            num_individuals: self.num_individuals,
            ploidy: self.ploidy,
            gts: Vec::with_capacity(self.gts_per_block),
            chrom: self
                .needs
                .contains(Needs::CHROM_POS)
                .then(|| Vec::with_capacity(num_vars)),
            pos: self
                .needs
                .contains(Needs::CHROM_POS)
                .then(|| Vec::with_capacity(num_vars)),
            id: self
                .needs
                .contains(Needs::ID)
                .then(|| Vec::with_capacity(num_vars)),
            alleles: self
                .needs
                .contains(Needs::ALLELES)
                .then(|| AllelesColumn::with_num_vars(num_vars)),
            qual: self
                .needs
                .contains(Needs::QUAL)
                .then(|| Vec::with_capacity(num_vars)),
        }
    }

    /// The variant the reader last filled, copied into the columns of the
    /// block. The reader gives the `num_individuals` x `ploidy` alleles of
    /// every variant that the trait asks of it, so the genotypes are one
    /// copy of a few kilobytes.
    fn push_variant(&self, block: &mut Block) {
        block.gts.extend_from_slice(&self.var.gts);
        if let Some(chrom) = block.chrom.as_mut() {
            chrom.push(self.var.chrom);
        }
        if let Some(pos) = block.pos.as_mut() {
            pos.push(self.var.pos);
        }
        if let Some(id) = block.id.as_mut() {
            id.push(self.var.id.clone());
        }
        if let Some(alleles) = block.alleles.as_mut() {
            alleles.push(&self.var.alleles);
        }
        if let Some(qual) = block.qual.as_mut() {
            // A variant with no quality is a NaN in the column, which is
            // what Python and TypeScript are given for it.
            qual.push(self.var.qual.unwrap_or(f32::NAN));
        }
    }
}

impl<R: VariantReader> fmt::Debug for BlockCollector<R> {
    /// What the collector was asked for and where it has got to. The reader
    /// is left out, so that a collector over a reader that has no `Debug`
    /// has one.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("BlockCollector")
            .field("needs", &self.needs)
            .field("num_vars_per_block", &self.num_vars_per_block)
            .field("num_individuals", &self.num_individuals)
            .field("ploidy", &self.ploidy)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::{BufReader, Cursor};
    use std::path::{Path, PathBuf};

    use super::{
        Block, BlockCollector, GENOTYPES_PER_BLOCK, MAX_NUM_VARS_PER_BLOCK, MIN_NUM_VARS_PER_BLOCK,
        default_num_vars_per_block,
    };
    use crate::error::{Error, Result};
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{ChromTable, MISSING_ALLELE, Needs, Variant, VariantReader};

    /// The reference VCFs live at the root of the repository, beside the
    /// Python tests that read the same files, and not inside this crate.
    /// The path is built from the directory of the manifest, so it holds
    /// whether the tests are run with `cargo test --workspace` or with
    /// `cargo test -p popnei`.
    fn reference_vcf(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/vcf")
            .join(name)
    }

    /// The 50 individuals of `many.vcf`, which
    /// `tests/reference/vcf/make_reference.py` writes as diploid.
    const MANY_INDIVIDUALS: usize = 50;
    const MANY_PLOIDY: usize = 2;

    /// The options of a VCF read with every variant given, the ones that
    /// failed a filter among them.
    fn every_variant() -> VcfOptions {
        VcfOptions {
            ploidy: MANY_PLOIDY,
            only_passed: false,
        }
    }

    /// A collector over one of the reference VCFs.
    fn collector_over(
        name: &str,
        options: VcfOptions,
        needs: Needs,
        num_vars_per_block: Option<usize>,
    ) -> BlockCollector<VcfReader<BufReader<File>>> {
        let reader = match VcfReader::from_path(&reference_vcf(name), options) {
            Ok(reader) => reader,
            Err(error) => panic!("{name}: {error}"),
        };
        match BlockCollector::new(reader, needs, num_vars_per_block) {
            Ok(collector) => collector,
            Err(error) => panic!("{name}: the collector was not built: {error}"),
        }
    }

    /// Every block a collector gives, until it has no more or it fails.
    fn blocks_of<R: VariantReader>(collector: &mut BlockCollector<R>) -> Result<Vec<Block>> {
        let mut blocks = Vec::new();
        while let Some(block) = collector.next_block()? {
            blocks.push(block);
        }
        Ok(blocks)
    }

    /// The blocks of one of the reference VCFs.
    fn blocks_read(
        name: &str,
        options: VcfOptions,
        needs: Needs,
        num_vars_per_block: Option<usize>,
    ) -> Vec<Block> {
        let mut collector = collector_over(name, options, needs, num_vars_per_block);
        match blocks_of(&mut collector) {
            Ok(blocks) => blocks,
            Err(error) => panic!("{name}: the collector stopped at {error}"),
        }
    }

    /// How many variants each block holds.
    fn num_vars_of(blocks: &[Block]) -> Vec<usize> {
        blocks.iter().map(|block| block.num_vars).collect()
    }

    /// How many genotypes each block holds, which the consumer that
    /// reshapes `gts` into variants x individuals x ploidy depends on.
    fn num_gts_of(blocks: &[Block]) -> Vec<usize> {
        blocks.iter().map(|block| block.gts.len()).collect()
    }

    #[test]
    fn the_blocks_of_a_hundred_variants_of_many_vcf_are_four_full_and_one_of_seventy_five() {
        // The default of the reader gives 475 of the 500 variants.
        let blocks = blocks_read("many.vcf", VcfOptions::default(), Needs::GTS, Some(100));
        assert_eq!(num_vars_of(&blocks), [100, 100, 100, 100, 75]);
        // 100 x 50 x 2 genotypes in a full block and 75 x 50 x 2 in the
        // last one.
        assert_eq!(num_gts_of(&blocks), [10_000, 10_000, 10_000, 10_000, 7_500]);
        for block in &blocks {
            assert_eq!(block.num_individuals, MANY_INDIVIDUALS);
            assert_eq!(block.ploidy, MANY_PLOIDY);
        }
    }

    #[test]
    fn every_variant_of_many_vcf_given_makes_five_blocks_of_a_hundred() {
        let blocks = blocks_read("many.vcf", every_variant(), Needs::GTS, Some(100));
        assert_eq!(num_vars_of(&blocks), [100, 100, 100, 100, 100]);
        assert_eq!(num_gts_of(&blocks), [10_000; 5]);
    }

    #[test]
    fn a_block_of_a_thousand_variants_holds_the_four_hundred_and_seventy_five_of_the_default() {
        let blocks = blocks_read("many.vcf", VcfOptions::default(), Needs::GTS, Some(1000));
        assert_eq!(num_vars_of(&blocks), [475]);
        // 475 x 50 x 2.
        assert_eq!(num_gts_of(&blocks), [47_500]);
    }

    /// One variant as `read_variant` gives it, which is what the blocks,
    /// joined, have to hold whatever their size.
    #[derive(Debug, PartialEq, Eq)]
    struct Site {
        chrom: u32,
        pos: u64,
        gts: Vec<i8>,
    }

    /// The variants of a VCF read one by one, which is what the blocks are
    /// compared with.
    fn sites_read_one_by_one(name: &str, options: VcfOptions) -> Vec<Site> {
        let mut reader = match VcfReader::from_path(&reference_vcf(name), options) {
            Ok(reader) => reader,
            Err(error) => panic!("{name}: {error}"),
        };
        reader.set_needs(Needs::GTS | Needs::CHROM_POS);
        let mut var = Variant::new();
        let mut sites = Vec::new();
        loop {
            match reader.read_variant(&mut var) {
                Ok(true) => sites.push(Site {
                    chrom: var.chrom,
                    pos: var.pos,
                    gts: var.gts.clone(),
                }),
                Ok(false) => return sites,
                Err(error) => panic!("{name}: the reader stopped at {error}"),
            }
        }
    }

    /// The variants of the blocks, joined, one by one.
    fn sites_of(blocks: &[Block]) -> Vec<Site> {
        let mut sites = Vec::new();
        for block in blocks {
            let chrom = block
                .chrom
                .as_ref()
                .expect("the chromosomes were asked for");
            let pos = block.pos.as_ref().expect("the positions were asked for");
            let gts_per_var = block.num_individuals.saturating_mul(block.ploidy);
            assert_eq!(chrom.len(), block.num_vars);
            assert_eq!(pos.len(), block.num_vars);
            assert_eq!(block.gts.len(), block.num_vars.saturating_mul(gts_per_var));
            for (index, gts) in block.gts.chunks_exact(gts_per_var).enumerate() {
                sites.push(Site {
                    chrom: chrom[index],
                    pos: pos[index],
                    gts: gts.to_vec(),
                });
            }
        }
        sites
    }

    #[test]
    fn the_blocks_joined_are_the_variants_the_reader_gives_one_by_one() {
        let expected = sites_read_one_by_one("many.vcf", VcfOptions::default());
        assert_eq!(expected.len(), 475);
        for num_vars_per_block in [1, 7, 100, 1000] {
            let blocks = blocks_read(
                "many.vcf",
                VcfOptions::default(),
                Needs::CHROM_POS,
                Some(num_vars_per_block),
            );
            let given = sites_of(&blocks);
            assert_eq!(
                given.len(),
                expected.len(),
                "blocks of {num_vars_per_block} variants: the number of variants"
            );
            for (index, (given, expected)) in given.iter().zip(&expected).enumerate() {
                assert_eq!(
                    given, expected,
                    "blocks of {num_vars_per_block} variants: the variant {index}, counted from 0"
                );
            }
        }
    }

    #[test]
    fn a_collector_asked_for_the_genotypes_alone_gives_a_block_with_no_other_column() {
        // The VCF reader fills the chromosome and the position of every
        // variant it gives, and the block still has no column for them.
        let blocks = blocks_read("many.vcf", VcfOptions::default(), Needs::GTS, Some(100));
        for block in &blocks {
            assert!(block.chrom.is_none());
            assert!(block.pos.is_none());
            assert!(block.id.is_none());
            assert!(block.alleles.is_none());
            assert!(block.qual.is_none());
        }
        assert_eq!(num_gts_of(&blocks), [10_000, 10_000, 10_000, 10_000, 7_500]);
    }

    #[test]
    fn the_default_number_of_variants_of_a_block_is_its_genotypes_between_the_two_bounds() {
        // 5 million genotypes divided by the individuals, which the
        // largest number of variants decides for 50 individuals and the
        // smallest one for 100000.
        assert_eq!(default_num_vars_per_block(50), 10_000);
        assert_eq!(default_num_vars_per_block(1000), 5_000);
        assert_eq!(default_num_vars_per_block(100_000), 100);
        assert_eq!(GENOTYPES_PER_BLOCK, 5_000_000);
        assert_eq!(MIN_NUM_VARS_PER_BLOCK, 100);
        assert_eq!(MAX_NUM_VARS_PER_BLOCK, 10_000);

        // The 50 individuals of `many.vcf` and no size asked for: its 475
        // variants are one block, since the default is 10000 of them.
        let blocks = blocks_read("many.vcf", VcfOptions::default(), Needs::GTS, None);
        assert_eq!(num_vars_of(&blocks), [475]);
    }

    /// The header of the VCFs written in these tests, with three
    /// individuals, and the first data line is its line 4.
    const HEADER: &str = "\
##fileformat=VCFv4.4
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\tind3
";

    /// A VCF of the tests: the header above and one data line for each line
    /// given, whose columns are written with spaces and separated by tabs
    /// in the file.
    fn vcf_of(lines: &[&str]) -> String {
        let mut vcf = HEADER.to_string();
        for line in lines {
            vcf.push_str(&line.replace(' ', "\t"));
            vcf.push('\n');
        }
        vcf
    }

    /// A collector over a VCF written in a test.
    fn collector_over_text(
        vcf: &str,
        options: VcfOptions,
        needs: Needs,
        num_vars_per_block: Option<usize>,
    ) -> BlockCollector<VcfReader<Cursor<Vec<u8>>>> {
        let reader = match VcfReader::new(Cursor::new(vcf.as_bytes().to_vec()), options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        };
        match BlockCollector::new(reader, needs, num_vars_per_block) {
            Ok(collector) => collector,
            Err(error) => panic!("the collector was not built: {error}"),
        }
    }

    #[test]
    fn a_reader_with_no_variants_gives_no_block() {
        let mut collector =
            collector_over_text(HEADER, VcfOptions::default(), Needs::ALL, Some(100));
        assert!(collector.next_block().unwrap().is_none());
        assert!(collector.next_block().unwrap().is_none());
    }

    #[test]
    fn the_error_of_the_third_variant_comes_after_the_block_of_the_two_before_it() {
        let vcf = vcf_of(&[
            "chr1 100 rs1 A T . PASS . GT 0/0 0/1 1/1",
            "chr1 200 rs2 A T . PASS . GT 0/0 0/1 1/1",
            "chr1 300 rs3 A T . PASS . GT 0/0 0/1 0/0/1/1",
        ]);
        let mut collector = collector_over_text(&vcf, VcfOptions::default(), Needs::GTS, Some(2));

        let first = collector.next_block().unwrap().expect("the first block");
        assert_eq!(first.num_vars, 2);
        assert_eq!(first.gts, [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]);

        let error = match collector.next_block() {
            Ok(block) => panic!("the collector gave {block:?}"),
            Err(error) => error,
        };
        let Error::VcfGenotypePloidy {
            line,
            individual,
            found,
            expected,
        } = error
        else {
            panic!("the error is {error}");
        };
        // The third data line of a VCF of three header lines.
        assert_eq!(line, 6);
        assert_eq!(individual, "ind3");
        assert_eq!((found, expected), (4, 2));

        // The variant of the block that was being built is lost with the
        // error, and there is no block after it.
        assert!(collector.next_block().unwrap().is_none());
    }

    #[test]
    fn a_collector_of_blocks_of_no_variant_is_refused() {
        let reader = VcfReader::new(
            Cursor::new(HEADER.as_bytes().to_vec()),
            VcfOptions::default(),
        )
        .expect("the reader");
        let error = match BlockCollector::new(reader, Needs::GTS, Some(0)) {
            Ok(collector) => panic!("the collector was built: {collector:?}"),
            Err(error) => error,
        };
        assert!(matches!(error, Error::BlockOfNoVariants), "{error}");
        let message = error.to_string();
        assert!(message.contains('0'), "{message}");
    }

    /// The genotypes of a block are its variants times the individuals
    /// times the ploidy, which is a multiplication that a size asked for by
    /// a caller carries beyond what a `usize` holds. It is an error and not
    /// a panic, on a machine of 32 bit addresses as on one of 64.
    #[test]
    fn a_block_of_more_genotypes_than_the_machine_addresses_is_refused() {
        // Two individuals of the ploidy 2 are four genotypes in every
        // variant, so every size above a quarter of `usize::MAX` is one.
        let error = match BlockCollector::new(
            FakeReader::giving(0, Fills::Everything),
            Needs::GTS,
            Some(usize::MAX),
        ) {
            Ok(collector) => panic!("the collector was built: {collector:?}"),
            Err(error) => error,
        };
        let Error::BlockTooLarge {
            num_vars_per_block,
            num_individuals,
            ploidy,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(
            (num_vars_per_block, num_individuals, ploidy),
            (usize::MAX, 2, 2)
        );
    }

    /// The quality of each variant of a block, with the NaN of a variant
    /// that has none as `None`, so that the assertion compares no float
    /// with another.
    fn quals_of(block: &Block) -> Vec<Option<f32>> {
        block
            .qual
            .as_ref()
            .expect("the qualities were asked for")
            .iter()
            .map(|qual| if qual.is_nan() { None } else { Some(*qual) })
            .collect()
    }

    /// The alleles of one variant of a block.
    fn alleles_of(block: &Block, var: usize) -> Vec<&str> {
        let alleles = block.alleles.as_ref().expect("the alleles were asked for");
        (0..alleles.num_alleles(var))
            .map(|allele| alleles.allele(var, allele))
            .collect()
    }

    #[test]
    fn every_column_of_cases_vcf_is_in_its_blocks() {
        // The four variants of the table of `cases.vcf` of "How it is
        // verified" of `docs/specs/io_vcf.md`, in blocks of three.
        let blocks = blocks_read("cases.vcf", every_variant(), Needs::ALL, Some(3));
        assert_eq!(num_vars_of(&blocks), [3, 1]);

        let first = &blocks[0];
        assert_eq!(first.num_individuals, 3);
        assert_eq!(first.ploidy, 2);
        // The chromosome of the four variants is `chr1`, the first name of
        // the file and so the number 0.
        assert_eq!(first.chrom.as_deref(), Some([0, 0, 0].as_slice()));
        assert_eq!(first.pos.as_deref(), Some([100, 200, 300].as_slice()));
        assert_eq!(
            first.id.as_deref(),
            Some(["rs1".to_string(), String::new(), String::new()].as_slice())
        );
        assert_eq!(quals_of(first), [Some(29.5), None, Some(67.0)]);
        assert_eq!(alleles_of(first, 0), ["A", "T"]);
        assert_eq!(alleles_of(first, 1), ["A", "T"]);
        assert_eq!(alleles_of(first, 2), ["A", "G", "T"]);
        // The genotypes of the three variants, `0/0 0/1 1/1`, `./. 0|1 .|0`
        // and `1/2 2|1 2/2`, one row for each of them.
        let gts: Vec<i8> = [
            [0, 0, 0, 1, 1, 1],
            [MISSING_ALLELE, MISSING_ALLELE, 0, 1, MISSING_ALLELE, 0],
            [1, 2, 2, 1, 2, 2],
        ]
        .concat();
        assert_eq!(first.gts, gts);

        let last = &blocks[1];
        assert_eq!(last.chrom.as_deref(), Some([0].as_slice()));
        assert_eq!(last.pos.as_deref(), Some([400].as_slice()));
        assert_eq!(last.id.as_deref(), Some([String::new()].as_slice()));
        assert_eq!(quals_of(last), [Some(47.0)]);
        // The variant with no alternative allele: the reference alone.
        assert_eq!(alleles_of(last, 0), ["T"]);
        assert_eq!(last.gts, [0, 0, 0, 0, 0, 0]);

        // A variant and an allele the column does not hold.
        let alleles = last.alleles.as_ref().expect("the alleles");
        assert_eq!(alleles.num_vars(), 1);
        assert_eq!(alleles.num_alleles(1), 0);
        assert_eq!(alleles.allele(0, 1), "");
        assert_eq!(alleles.allele(1, 0), "");
    }

    #[test]
    fn a_block_of_a_tetraploid_vcf_holds_four_alleles_for_each_individual() {
        let vcf = vcf_of(&[
            "chr1 100 rs1 A T . PASS . GT 0/0/1/1 0/1/1/1 ./././.",
            "chr1 200 rs2 A T . PASS . GT 0/0/0/0 1/1/1/1 0/0/0/1",
            "chr1 300 rs3 A T . PASS . GT 0/0/0/0 1/1/1/1 0/0/0/1",
        ]);
        let options = VcfOptions {
            ploidy: 4,
            only_passed: true,
        };
        let mut collector = collector_over_text(&vcf, options, Needs::GTS, Some(2));
        let blocks = blocks_of(&mut collector).expect("the blocks");

        assert_eq!(num_vars_of(&blocks), [2, 1]);
        assert_eq!(blocks[0].ploidy, 4);
        assert_eq!(blocks[0].num_individuals, 3);
        // 2 variants x 3 individuals x 4 alleles, and 1 x 3 x 4.
        assert_eq!(num_gts_of(&blocks), [24, 12]);
        assert_eq!(
            blocks[0].gts.get(..8),
            Some([0, 0, 1, 1, 0, 1, 1, 1].as_slice())
        );
        assert_eq!(
            blocks[0].gts.get(8..12),
            Some([MISSING_ALLELE; 4].as_slice())
        );
    }

    /// What a reader written for these tests puts in the variants it
    /// gives.
    #[derive(Clone, Copy)]
    enum Fills {
        /// The chromosome, the position and the genotypes of its two
        /// individuals.
        Everything,
        /// The chromosome and the position alone, which is a reader whose
        /// source has no genotype to give.
        NoGenotypes,
    }

    /// A reader of two individuals of the ploidy 2, written for these
    /// tests: a VCF cannot be made to do what a collector has to stand.
    struct FakeReader {
        individuals: Vec<String>,
        chroms: ChromTable,
        /// How many variants it still has to give.
        left: usize,
        fills: Fills,
        /// The variant it fails at, counted from 1.
        fails_at: Option<usize>,
        /// Which variant it is about to give, counted from 1.
        next_var: usize,
        /// What the collector asked this reader for.
        needs: Needs,
    }

    impl FakeReader {
        /// A reader of `num_vars` variants that fills each of them as
        /// `fills` says.
        fn giving(num_vars: usize, fills: Fills) -> FakeReader {
            FakeReader {
                individuals: vec!["ind1".to_string(), "ind2".to_string()],
                chroms: ChromTable::new(),
                left: num_vars,
                fills,
                fails_at: None,
                next_var: 1,
                needs: Needs::empty(),
            }
        }

        /// A reader of `num_vars` variants whose variant number `fails_at`,
        /// counted from 1, is an error, and which gives the variants after
        /// it. The trait says that a reader ends at its error, so a
        /// collector that reads on after one reads what no reader of popnei
        /// gives: this one shows what it would do with it.
        fn failing_at(num_vars: usize, fails_at: usize) -> FakeReader {
            FakeReader {
                fails_at: Some(fails_at),
                ..FakeReader::giving(num_vars, Fills::Everything)
            }
        }
    }

    impl VariantReader for FakeReader {
        fn read_variant(&mut self, var: &mut Variant) -> Result<bool> {
            var.clear();
            let Some(left) = self.left.checked_sub(1) else {
                return Ok(false);
            };
            self.left = left;
            let this_var = self.next_var;
            self.next_var = this_var.saturating_add(1);
            if self.fails_at == Some(this_var) {
                return Err(Error::Io(std::io::Error::other(
                    "the reader of the tests failed",
                )));
            }
            var.chrom = self.chroms.intern("chr1");
            var.pos = 100;
            var.filled = Needs::CHROM_POS;
            match self.fills {
                Fills::Everything => {
                    var.gts.extend_from_slice(&[0, 1, 1, 1]);
                    var.filled |= Needs::GTS;
                }
                Fills::NoGenotypes => {}
            }
            Ok(true)
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
    }

    /// A column nobody wants is never parsed, so what the collector asks
    /// its reader for is what its blocks will hold.
    #[test]
    fn a_collector_asks_its_reader_for_its_columns_and_the_genotypes() {
        let no_variant = || FakeReader::giving(0, Fills::Everything);
        let collector =
            BlockCollector::new(no_variant(), Needs::QUAL, Some(2)).expect("the collector");
        assert_eq!(collector.reader().needs, Needs::GTS | Needs::QUAL);

        let collector =
            BlockCollector::new(no_variant(), Needs::empty(), None).expect("the collector");
        assert_eq!(collector.reader().needs, Needs::GTS);

        let collector =
            BlockCollector::new(no_variant(), Needs::ALL, Some(2)).expect("the collector");
        assert_eq!(collector.reader().needs, Needs::ALL);
    }

    /// The trait says that a reader ends at its error, and the collector
    /// does not lean on it: a reader that goes on giving variants after one
    /// gives the collector's caller no block after it.
    #[test]
    fn the_collector_gives_no_block_after_the_error_of_its_reader() {
        let mut collector =
            BlockCollector::new(FakeReader::failing_at(6, 3), Needs::CHROM_POS, Some(2))
                .expect("the collector");

        let first = collector
            .next_block()
            .expect("the first block")
            .expect("a block");
        assert_eq!(first.num_vars, 2);

        let error = match collector.next_block() {
            Ok(block) => panic!("the collector gave {block:?}"),
            Err(error) => error,
        };
        assert!(matches!(error, Error::Io(_)), "the error is {error}");

        // The reader has three variants left and gives them, and the
        // collector gives no block.
        assert!(collector.next_block().expect("no block").is_none());
        assert!(collector.next_block().expect("no block").is_none());
    }

    /// The two binding crates hold their collector over a boxed reader,
    /// `BlockCollector<Box<dyn VariantReader>>`, so one is built here.
    #[test]
    fn a_field_that_was_asked_for_and_that_the_reader_does_not_fill_is_an_error() {
        let reader: Box<dyn VariantReader> = Box::new(FakeReader::giving(2, Fills::NoGenotypes));
        let mut collector =
            BlockCollector::new(reader, Needs::CHROM_POS, Some(2)).expect("the collector");
        assert_eq!(collector.reader().individuals(), ["ind1", "ind2"]);
        assert_eq!(collector.reader().ploidy(), 2);
        let error = match collector.next_block() {
            Ok(block) => panic!("the collector gave {block:?}"),
            Err(error) => error,
        };
        let Error::FieldsNotFilled { fields } = error else {
            panic!("the error is {error}");
        };
        assert_eq!(fields, Needs::GTS);
    }
}
