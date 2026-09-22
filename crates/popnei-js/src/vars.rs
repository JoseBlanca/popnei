//! What a TypeScript user reaches through `openVars`: the bytes of a vars
//! file, the arrow file popnei keeps its variants in, read as a source of
//! variants.
//!
//! [`VarsSource`] holds the bytes of the file and reads its schema and its
//! footer when it is built, so bytes that are not a vars file fail at
//! `openVars`, and the names of the individuals and the ploidy come from the
//! `popnei` key of that schema. Every pass over it reads the same bytes
//! again from their start and goes through the `Blocks` of `source.rs`, the
//! one a VCF goes through.
//!
//! A tab has no filesystem, as section 11 of `docs/architecture.md` says, so
//! the file a user gets is the bytes of one: `write_vars` of `source.rs`
//! builds it in the memory of wasm and it crosses as a `Uint8Array` that the
//! page offers as a download. That is the whole file in memory beside the
//! source it was written from.

use std::sync::Arc;

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::block::BlockReader;
use popnei::io::vars::VarsReader;

use crate::errors::JsPopneiError;
use crate::pca::{PcaOfVariants, pca_of_the_variants};
use crate::source::{Blocks, OpenSource, VarsFile, blocks_of, bytes_of_a_vars_file, cursor_of};
use crate::steps::Steps;

/// A vars file that was opened: its bytes, and the individuals and the
/// ploidy its schema named.
#[wasm_bindgen]
pub struct VarsSource {
    /// The bytes of the whole file, which every pass over them shares.
    bytes: Arc<Vec<u8>>,
    individuals: Vec<String>,
    ploidy: usize,
}

#[wasm_bindgen]
impl VarsSource {
    /// The names of the individuals, in the order of their genotypes in the
    /// rows of a block.
    #[must_use]
    pub fn individuals(&self) -> Vec<String> {
        self.individuals.clone()
    }

    /// How many alleles the genotype of one individual holds.
    #[must_use]
    pub fn ploidy(&self) -> usize {
        self.ploidy
    }

    /// One pass over the bytes, read again from their start, through the
    /// steps of `steps`: its blocks hold `fields` besides the genotypes,
    /// `num_vars_per_block` variants each.
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
        steps: Steps,
    ) -> Result<Blocks, JsPopneiError> {
        blocks_of(self, fields, num_vars_per_block, steps)
    }

    /// The variants of the file, through the steps of `steps`, as another
    /// vars file, of batches of `num_vars_per_block` variants, and of the
    /// size popnei chooses for these individuals when it is not given,
    /// which the package reads out of the memory of wasm piece by piece,
    /// with the counts of the pass that wrote them.
    ///
    /// # Errors
    ///
    /// When `num_vars_per_block` is 0, when the file cannot be read, when a
    /// block of it is not one a vars file holds, and when the memory of the
    /// tab does not take the file.
    pub fn write_vars(
        &self,
        num_vars_per_block: Option<usize>,
        steps: Steps,
    ) -> Result<VarsFile, JsPopneiError> {
        bytes_of_a_vars_file(self, num_vars_per_block, steps)
    }

    /// The principal components of the variants of the file, through the
    /// steps of `steps`, with the weights of the first `num_prin_comps`
    /// components.
    ///
    /// # Errors
    ///
    /// When the file cannot be read, when a variant has more than two
    /// alleles among its called genotypes and `transform_to_biallelic` is
    /// false, when the pass gives no variant or no variant with variance,
    /// when a size of the dataset is beyond what the analysis counts in, and
    /// when the linear algebra could not be done.
    pub fn pca_of_variants(
        &self,
        transform_to_biallelic: bool,
        num_prin_comps: usize,
        steps: Steps,
    ) -> Result<PcaOfVariants, JsPopneiError> {
        pca_of_the_variants(self, transform_to_biallelic, num_prin_comps, steps)
    }
}

impl OpenSource for VarsSource {
    /// The size the caller asks for is not passed on: the reader gives each
    /// batch of the file as a block, at the size the file was written with,
    /// and the `Reblock` that every pass ends with cuts them where the
    /// caller wants them.
    fn reader(
        &self,
        _num_vars_per_block: Option<usize>,
    ) -> Result<Box<dyn BlockReader>, popnei::Error> {
        Ok(Box::new(VarsReader::new(cursor_of(&self.bytes))?))
    }
}

/// The vars file in `bytes`.
///
/// It reads the schema and the footer, so the individuals, the ploidy and
/// the batches are known when it returns and bytes that are not a vars file
/// fail here.
///
/// # Errors
///
/// When the bytes are not a vars file that popnei can read: what
/// `docs/specs/io_vars.md` lists as an error of the file as a whole, among
/// them bytes that are not an arrow file, a schema with no `popnei` key, a
/// format version popnei does not read, a column of another type and a
/// footer whose entries are not as many as the batches.
#[wasm_bindgen]
pub fn open_vars(bytes: Vec<u8>) -> Result<VarsSource, JsPopneiError> {
    // The `Vec` wasm-bindgen filled with the bytes of the `Uint8Array` is
    // the one every pass reads: an `Arc<[u8]>` here would allocate the whole
    // file again and copy it into the new buffer, and the memory of wasm
    // never gives that back.
    let bytes = Arc::new(bytes);
    // The schema and the footer are read when the reader is built and no
    // batch is, so a file whose batches would need more memory than wasm
    // addresses is opened all the same and its individuals read.
    let reader = VarsReader::new(cursor_of(&bytes))?;
    let metadata = reader.metadata();
    let individuals = metadata.individuals.clone();
    let ploidy = metadata.ploidy;
    Ok(VarsSource {
        bytes,
        individuals,
        ploidy,
    })
}
