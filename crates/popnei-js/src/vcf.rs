//! What a TypeScript user reaches through `openVcf`: the bytes of a VCF with
//! its options, read as a source of variants.
//!
//! [`VcfSource`] holds the bytes of a VCF and the options it is read with,
//! and it reads the header when it is built, so bytes that are not a VCF
//! fail at `openVcf`. Every pass over it reads the bytes again from their
//! start and goes through the `Blocks` of `source.rs`, the one a vars file
//! goes through, and `write_vars` there is what writes its variants into a
//! vars file.

use std::sync::Arc;

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::block::BlockReader;
use popnei::io::vcf::{VcfOptions, VcfReader};

use crate::errors::JsPopneiError;
use crate::pca::{PcaOfVariants, pca_of_the_variants};
use crate::source::{Blocks, OpenSource, VarsFile, blocks_of, bytes_of_a_vars_file, cursor_of};
use crate::steps::Steps;

/// A VCF that was opened: its bytes, the options it is read with, and the
/// individuals its header named.
#[wasm_bindgen]
pub struct VcfSource {
    /// The bytes of the whole file, which every pass over them shares.
    bytes: Arc<Vec<u8>>,
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

    /// The variants of the VCF, through the steps of `steps`, as a vars
    /// file of batches of `num_vars_per_block` variants, and of the size
    /// popnei chooses for these individuals when it is not given, which the
    /// package reads out of the memory of wasm piece by piece, with the
    /// counts of the pass that wrote them.
    ///
    /// # Errors
    ///
    /// When `num_vars_per_block` is 0, when the VCF cannot be read, when a
    /// block of it is not one a vars file holds, and when the memory of the
    /// tab does not take the file.
    pub fn write_vars(
        &self,
        num_vars_per_block: Option<usize>,
        steps: Steps,
    ) -> Result<VarsFile, JsPopneiError> {
        bytes_of_a_vars_file(self, num_vars_per_block, steps)
    }

    /// The principal components of the variants of the VCF, through the
    /// steps of `steps`, with the weights of the first `num_prin_comps`
    /// components.
    ///
    /// # Errors
    ///
    /// When the analysis of these individuals does not fit in the memory of
    /// a page, when the VCF cannot be read, when a variant has more than two
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
        pca_of_the_variants(
            self,
            self.individuals.len(),
            transform_to_biallelic,
            num_prin_comps,
            steps,
        )
    }
}

impl OpenSource for VcfSource {
    fn reader(
        &self,
        num_vars_per_block: Option<usize>,
    ) -> Result<Box<dyn BlockReader>, popnei::Error> {
        let options = VcfOptions {
            num_vars_per_block,
            ..self.options
        };
        Ok(Box::new(VcfReader::new(cursor_of(&self.bytes), options)?))
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
        num_vars_per_block: None,
    };
    // The `Vec` wasm-bindgen filled with the bytes of the `Uint8Array` is
    // the one every pass reads: an `Arc<[u8]>` here would allocate the whole
    // file again and copy it into the new buffer, and the memory of wasm
    // never gives that back.
    let bytes = Arc::new(bytes);
    // The header is read when the reader is built and no variant is.
    // Nothing here asks for a block, so a file whose blocks would need more
    // memory than wasm addresses, a header of 170000 individuals read with
    // the ploidy 255, is opened all the same and its individuals read; the
    // size of its blocks is the user's to choose at `iterBlocks`.
    let reader = VcfReader::new(cursor_of(&bytes), options)?;
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
