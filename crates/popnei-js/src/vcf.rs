//! What a TypeScript user reaches through `openVcf` and `Variants`: the
//! bytes of a VCF with its options, and the blocks of its variants.
//!
//! Three classes. [`VcfSource`] holds the bytes of a VCF and the options it
//! is read with, and it reads the header when it is built, so bytes that are
//! not a VCF fail at `openVcf`. [`Blocks`] is one pass over those bytes: it
//! owns a reader and the collector of `popnei::block`, and every call of
//! `VcfSource::blocks` reads the bytes again from their start, which is what
//! lets a user give the same `Variants` to one calculation after another.
//! [`BlockColumns`] is one block on its way out.
//!
//! The three live in the memory of wasm, which the garbage collector of
//! JavaScript does not see, so the TypeScript package frees each of them:
//! the block as soon as its columns are read, the pass when the iteration
//! ends or is abandoned, and the source when the user calls `free()`.
//!
//! The columns leave as JavaScript arrays, which wasm-bindgen copies out of
//! the memory of wasm: a typed array that is a view into that memory stops
//! being valid when it grows, as section 11 of `docs/architecture.md` says.
//! Each column is moved out of the block as it is read, so that the copy
//! that crosses is the only one.

use std::io::Cursor;
use std::sync::Arc;

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::block::{AllelesColumn, Block, BlockCollector};
use popnei::io::vcf::{VcfOptions, VcfReader};
use popnei::variant::{Needs, VariantReader};

use crate::errors::JsPopneiError;

/// The name each field of a block has in TypeScript, with the fields of the
/// core it asks the reader for. The chromosome and the position are one
/// field of the core, so either name fills both.
const FIELD_NAMES: [(&str, Needs); 5] = [
    ("chrom", Needs::CHROM_POS),
    ("pos", Needs::CHROM_POS),
    ("id", Needs::ID),
    ("alleles", Needs::ALLELES),
    ("qual", Needs::QUAL),
];

/// A VCF that was opened: its bytes, the options it is read with, and the
/// individuals its header named.
#[wasm_bindgen]
pub struct VcfSource {
    /// The bytes of the whole file, shared with every pass over them: a
    /// pass reads them again and none of them is copied.
    bytes: Arc<[u8]>,
    options: VcfOptions,
    individuals: Vec<String>,
}

#[wasm_bindgen]
impl VcfSource {
    /// The names of the individuals, in the order of the columns of the
    /// VCF.
    #[must_use]
    pub fn individuals(&self) -> Vec<String> {
        self.individuals.clone()
    }

    /// How many alleles the genotype of one individual holds.
    #[must_use]
    pub fn ploidy(&self) -> usize {
        self.options.ploidy
    }

    /// One pass over the bytes, read again from their start: its blocks
    /// hold `fields` besides the genotypes, `num_vars_per_block` variants
    /// each.
    ///
    /// # Errors
    ///
    /// When a name of `fields` is not a field of a block, when
    /// `num_vars_per_block` is 0, and when the genotypes of one block are
    /// more than the memory of wasm addresses.
    pub fn blocks(
        &self,
        fields: Vec<String>,
        num_vars_per_block: Option<usize>,
    ) -> Result<Blocks, JsPopneiError> {
        let needs = needs_of(&fields)?;
        let source = Cursor::new(Arc::clone(&self.bytes));
        let reader: Box<dyn VariantReader> = Box::new(VcfReader::new(source, self.options)?);
        Ok(Blocks {
            collector: BlockCollector::new(reader, needs, num_vars_per_block)?,
        })
    }
}

/// One pass over a VCF, which gives its variants block by block.
#[wasm_bindgen]
pub struct Blocks {
    collector: BlockCollector<Box<dyn VariantReader>>,
}

#[wasm_bindgen]
impl Blocks {
    /// The next block of the pass, or `undefined` when the VCF has no more
    /// variants.
    ///
    /// The names of the chromosomes are taken after the block was
    /// collected, as `docs/specs/block.md` says: the table of the reader
    /// grows while the file is read, and a block holds numbers of it.
    ///
    /// # Errors
    ///
    /// When a variant cannot be read. The block that was being built is
    /// lost with the error, and every call after it gives no block.
    pub fn next_block(&mut self) -> Result<Option<BlockColumns>, JsPopneiError> {
        let Some(block) = self.collector.next_block()? else {
            return Ok(None);
        };
        let chroms = self.collector.reader().chroms();
        let names = |numbers: Vec<u32>| {
            numbers
                .into_iter()
                .map(|number| {
                    chroms.name(number).map(str::to_owned).ok_or_else(|| {
                        JsPopneiError::Broken(format!(
                            "the chromosome number {number} of a block is not in the \
                             table of the reader that gave it"
                        ))
                    })
                })
                .collect::<Result<Vec<String>, JsPopneiError>>()
        };
        let Block {
            num_vars,
            num_individuals,
            ploidy,
            gts,
            chrom,
            pos,
            id,
            alleles,
            qual,
        } = block;
        let alleles = alleles.map(|column| alleles_of(&column)).transpose()?;
        let (alleles, num_alleles_per_var) = match alleles {
            Some((texts, counts)) => (Some(texts), Some(counts)),
            None => (None, None),
        };
        Ok(Some(BlockColumns {
            num_vars,
            num_individuals,
            ploidy,
            gts: Some(gts),
            chrom: chrom.map(names).transpose()?,
            pos: pos.map(positions_of),
            id,
            alleles,
            num_alleles_per_var,
            qual,
        }))
    }
}

/// The columns of one block on their way to TypeScript.
///
/// Each of them leaves the block the first time it is asked for, and the
/// call after that gives nothing: the package reads every column once, into
/// the object a user holds, and frees the block. What a column costs is then
/// one copy, the one wasm-bindgen makes when it crosses.
#[wasm_bindgen]
pub struct BlockColumns {
    num_vars: usize,
    num_individuals: usize,
    ploidy: usize,
    gts: Option<Vec<i8>>,
    chrom: Option<Vec<String>>,
    pos: Option<Vec<f64>>,
    id: Option<Vec<String>>,
    /// The alleles of every variant one variant after another, and how many
    /// of them each variant has: a JavaScript array of arrays of strings is
    /// not one of the types wasm-bindgen carries, so the column crosses as
    /// the core holds it, flat, and the package cuts it into the alleles of
    /// each variant.
    alleles: Option<Vec<String>>,
    num_alleles_per_var: Option<Vec<u32>>,
    qual: Option<Vec<f32>>,
}

#[wasm_bindgen]
impl BlockColumns {
    /// How many variants the block holds.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// How many individuals the source has, the same for every variant.
    #[must_use]
    pub fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// How many alleles the genotype of one individual holds.
    #[must_use]
    pub fn ploidy(&self) -> usize {
        self.ploidy
    }

    /// The genotypes, `num_vars` x `num_individuals` x `ploidy` alleles,
    /// variant after variant and inside a variant individual after
    /// individual.
    pub fn gts(&mut self) -> Vec<i8> {
        self.gts.take().unwrap_or_default()
    }

    /// The name of the chromosome of each variant, or `undefined` when the
    /// pass was not asked for it.
    pub fn chrom(&mut self) -> Option<Vec<String>> {
        self.chrom.take()
    }

    /// The position of each variant, 1 based as in a VCF.
    pub fn pos(&mut self) -> Option<Vec<f64>> {
        self.pos.take()
    }

    /// The id of each variant, an empty text for a variant that has none.
    pub fn id(&mut self) -> Option<Vec<String>> {
        self.id.take()
    }

    /// The alleles of every variant, the reference one first, one variant
    /// after another: the first `num_alleles_per_var()[0]` of them are the
    /// alleles of the first variant, and so on.
    pub fn alleles(&mut self) -> Option<Vec<String>> {
        self.alleles.take()
    }

    /// How many alleles each variant has, the reference one counted.
    pub fn num_alleles_per_var(&mut self) -> Option<Vec<u32>> {
        self.num_alleles_per_var.take()
    }

    /// The quality of each variant, phred scaled as the QUAL of a VCF, and
    /// NaN for a variant that has none.
    pub fn qual(&mut self) -> Option<Vec<f32>> {
        self.qual.take()
    }
}

/// The VCF in `bytes`, plain or gzipped, read with `ploidy` alleles in every
/// genotype and, when `only_passed` is true, without the variants that
/// failed a filter.
///
/// It reads the header, so the individuals are known when it returns.
///
/// # Errors
///
/// When the bytes are not a VCF that popnei can read, and when the ploidy is
/// out of the range the core takes.
#[wasm_bindgen]
pub fn open_vcf(
    bytes: Vec<u8>,
    ploidy: usize,
    only_passed: bool,
) -> Result<VcfSource, JsPopneiError> {
    let options = VcfOptions {
        ploidy,
        only_passed,
    };
    let bytes: Arc<[u8]> = Arc::from(bytes);
    let reader = VcfReader::new(Cursor::new(Arc::clone(&bytes)), options)?;
    let individuals = reader.individuals().to_vec();
    Ok(VcfSource {
        bytes,
        options,
        individuals,
    })
}

/// The ploidy a VCF is read with when the caller says nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_ploidy() -> usize {
    popnei::io::vcf::DEFAULT_PLOIDY
}

/// Whether the variants that failed a filter are left out when the caller
/// says nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_only_passed() -> bool {
    popnei::io::vcf::DEFAULT_ONLY_PASSED
}

/// The fields of the core that the names of `fields` ask for, the genotypes
/// among them, which every block holds.
fn needs_of(fields: &[String]) -> Result<Needs, JsPopneiError> {
    let mut needs = Needs::GTS;
    for field in fields {
        let found = FIELD_NAMES
            .iter()
            .find(|(name, _)| *name == field.as_str())
            .map(|(_, of_the_core)| *of_the_core);
        let Some(of_the_core) = found else {
            let names = FIELD_NAMES.map(|(name, _)| format!("`{name}`")).join(", ");
            return Err(JsPopneiError::Argument(format!(
                "`{field}` is not a field of a block; the fields are {names}"
            )));
        };
        needs = needs.union(of_the_core);
    }
    Ok(needs)
}

/// The positions of a block as the numbers of JavaScript, which are float64
/// and hold a position of up to 2^53 exactly, as `docs/specs/block.md` says.
fn positions_of(positions: Vec<u64>) -> Vec<f64> {
    // A position above 2^53 would lose its last digits, and no genome comes
    // near it: the longest chromosome that has been assembled is 2.5e8
    // bases.
    positions.into_iter().map(|pos| pos as f64).collect()
}

/// The alleles of every variant of a block, the reference one first, one
/// variant after another, and how many of them each variant has.
///
/// # Errors
///
/// When a variant has more alleles than a `u32` counts, which no source
/// gives: a variant of a VCF has as many alleles as its ALT column names.
fn alleles_of(column: &AllelesColumn) -> Result<(Vec<String>, Vec<u32>), JsPopneiError> {
    let mut texts = Vec::new();
    let mut counts = Vec::with_capacity(column.num_vars());
    for var in 0..column.num_vars() {
        let num_alleles = column.num_alleles(var);
        texts.extend((0..num_alleles).map(|allele| column.allele(var, allele).to_owned()));
        counts.push(u32::try_from(num_alleles).map_err(|_| {
            JsPopneiError::Broken(format!(
                "the variant {var} of a block has {num_alleles} alleles, more than a \
                 JavaScript array of counts holds"
            ))
        })?);
    }
    Ok((texts, counts))
}
