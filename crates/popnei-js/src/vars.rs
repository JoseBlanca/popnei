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

use crate::dists::{KosmanDistances, kosman_dists_of};
use crate::errors::JsPopneiError;
use crate::kinship::{KinshipOfVariants, kinship_of_the_variants};
use crate::ld::{R2Matrix, r2_matrix_of};
use crate::pca::{PcaOfVariants, pca_of_the_variants};
use crate::pop_dists::{ArgumentsOfTheDists, PopDistsOfAPass, pop_dists_of};
use crate::source::{Blocks, OpenSource, VarsFile, blocks_of, bytes_of_a_vars_file, cursor_of};
use crate::stats::{
    ArgumentsOfThePass, PerIndividualStats, PerVarDistribs, per_individual_stats_of,
    per_var_distribs_of,
};
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

    /// The five per variant statistics of one pass over the file, through
    /// the steps of `steps`, for each population of `pop_names`.
    ///
    /// The arguments are those of `calcPerVarDistribs` of
    /// `docs/specs/stats.md`, as the package checked them and flat: `stats`
    /// holds the name of each statistic to calculate; the populations are
    /// their names, the names of the individuals of every one of them one
    /// after another, and how many individuals each of them holds, and
    /// `pop_names` is nothing when the user named no population, which is
    /// one population of every individual of the pass;
    /// `min_num_individuals` is how many called genotypes a population needs
    /// at a variant to have a value there; `hist_start`, `hist_end`,
    /// `num_bins` and `bin_type` are the histogram every statistic is
    /// counted in; `ploidy` is the exponent of the two expected
    /// heterozygosities, and nothing for the ploidy of the variants; and
    /// `poly_threshold` is the major allele frequency below which a variant
    /// is polymorphic in a population.
    ///
    /// # Errors
    ///
    /// Those of [`per_var_distribs_of`]: a name that is of no statistic, a
    /// histogram that cannot be made of what was given, an exponent of 0 or
    /// above 255, a population that names an individual the pass does not
    /// give, names one twice or names none, `pops` with no population, a
    /// polymorphism threshold that is no frequency, a source that cannot be
    /// read, and a pass that gives no variant.
    #[expect(
        clippy::too_many_arguments,
        reason = "the arguments of `calcPerVarDistribs` of `docs/specs/stats.md`, each \
                  one as the package checked it, and the populations flat: an array of \
                  arrays is not one of the types wasm-bindgen carries"
    )]
    pub fn calc_per_var_distribs(
        &self,
        steps: Steps,
        stats: Vec<String>,
        pop_names: Option<Vec<String>>,
        pop_individuals: Vec<String>,
        num_individuals_per_pop: Vec<u32>,
        min_num_individuals: u32,
        hist_start: f64,
        hist_end: f64,
        num_bins: usize,
        bin_type: String,
        ploidy: Option<usize>,
        poly_threshold: f64,
    ) -> Result<PerVarDistribs, JsPopneiError> {
        per_var_distribs_of(
            self,
            &steps,
            &ArgumentsOfThePass {
                stats,
                pop_names,
                pop_individuals,
                num_individuals_per_pop,
                min_num_individuals,
                hist_range: (hist_start, hist_end),
                num_bins,
                bin_type,
                ploidy,
                poly_threshold,
            },
        )
    }

    /// The missing rate and the heterozygosity rate of every individual of
    /// one pass over the file, through the steps of `steps`.
    ///
    /// # Errors
    ///
    /// Those of [`per_individual_stats_of`]: a source that cannot be read,
    /// and a pass that gives no variant.
    pub fn calc_per_individual_stats(
        &self,
        steps: Steps,
    ) -> Result<PerIndividualStats, JsPopneiError> {
        per_individual_stats_of(self, &steps)
    }

    /// The principal components of the variants of the file, through the
    /// steps of `steps`, with the weights of the first `num_prin_comps`
    /// components.
    ///
    /// # Errors
    ///
    /// When the analysis of these individuals does not fit in the memory of
    /// a page, when the file cannot be read, when a variant has more than
    /// two alleles among its called genotypes and `transform_to_biallelic`
    /// is false, when the pass gives no variant or no variant with variance,
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

    /// The kinship of every pair of the individuals of `individuals`, or of
    /// every individual of the pass when it is nothing, over the variants of
    /// the file that the steps of `steps` keep.
    ///
    /// # Errors
    ///
    /// When a name of `individuals` is of nobody the pass gives, is there
    /// twice, or the list is empty; when the file cannot be read; when a
    /// variant has more than two alleles among its called genotypes and
    /// `transform_to_biallelic` is false; when the pass gives no variant or
    /// no variant with variance; when two individuals have no variant called
    /// in both; when a size of the dataset is beyond what the calculation
    /// counts in; and when the linear algebra could not be done.
    pub fn calc_kinship(
        &self,
        individuals: Option<Vec<String>>,
        transform_to_biallelic: bool,
        steps: Steps,
    ) -> Result<KinshipOfVariants, JsPopneiError> {
        kinship_of_the_variants(self, individuals, transform_to_biallelic, steps)
    }

    /// The Kosman distance of every pair of individuals over the variants
    /// of the file that the steps of `steps` keep, with no distance for a
    /// pair called at fewer than `min_num_vars` variants, and the counts of
    /// the pass that gave them.
    ///
    /// # Errors
    ///
    /// When the pass gives no variant, when the sums of a pair go above
    /// what a `u32` holds, when the memory of the tab does not take the two
    /// counts of every pair, and when the file cannot be read.
    pub fn calc_pairwise_kosman_dists(
        &self,
        min_num_vars: u32,
        steps: Steps,
    ) -> Result<KosmanDistances, JsPopneiError> {
        kosman_dists_of(self, min_num_vars, steps)
    }

    /// The r² of every pair of the variants of the file that the steps of
    /// `steps` keep, with the chromosome and the position of each of them
    /// and the counts of the pass, over at most `max_num_vars` variants.
    ///
    /// # Errors
    ///
    /// When the pass gives more than `max_num_vars` variants, when the
    /// matrix of that many holds more values than wasm counts, when the pass
    /// gives no variant, when the memory of the tab does not take the
    /// matrix, and when the file cannot be read.
    pub fn calc_rogers_huff_r2_matrix(
        &self,
        max_num_vars: usize,
        steps: Steps,
    ) -> Result<R2Matrix, JsPopneiError> {
        r2_matrix_of(self, max_num_vars, steps)
    }

    /// Every measure of `measures` for every pair of the populations that
    /// were named, over one pass over the file through the steps of
    /// `steps`.
    ///
    /// The arguments are those of `calcPopDists` of `docs/specs/dists.md`,
    /// as the package checked them and flat: the populations are their
    /// names, the names of the individuals of every one of them one after
    /// another, and how many individuals each of them holds; `measures`
    /// holds the name of each measure to calculate; `group_per_variant` and
    /// `group_base_pairs` are the `jackknifeGroup` the user wrote, each
    /// variant its own group or a length in base pairs, and neither of them
    /// is no standard error; and `min_num_individuals` is how many called
    /// genotypes a population needs at a variant for that variant to count
    /// for a pair.
    ///
    /// # Errors
    ///
    /// Those of [`pop_dists_of`]: a name that is of no measure, a length of
    /// the resampling groups that is no whole number of base pairs of 1 or
    /// more, a population that names an individual the pass does not give,
    /// names one twice or names none, fewer than two populations, a pass
    /// that gives no variant, fewer resampling groups than a standard error
    /// is built from, sums the memory of the tab does not take, a group
    /// whose chromosome or positions JavaScript cannot hold, a pair of more
    /// variants than a JavaScript array of counts holds, and a source that
    /// cannot be read.
    #[expect(
        clippy::too_many_arguments,
        reason = "the arguments of `calcPopDists` of `docs/specs/dists.md`, each one as \
                  the package checked it, and the populations flat: an array of arrays \
                  is not one of the types wasm-bindgen carries"
    )]
    pub fn calc_pop_dists(
        &self,
        steps: Steps,
        pop_names: Vec<String>,
        pop_individuals: Vec<String>,
        num_individuals_per_pop: Vec<u32>,
        measures: Vec<String>,
        group_per_variant: bool,
        group_base_pairs: Option<f64>,
        min_num_individuals: u32,
    ) -> Result<PopDistsOfAPass, JsPopneiError> {
        pop_dists_of(
            self,
            &steps,
            &ArgumentsOfTheDists {
                pop_names,
                pop_individuals,
                num_individuals_per_pop,
                measures,
                group_per_variant,
                group_base_pairs,
                min_num_individuals,
            },
        )
    }
}

impl OpenSource for VarsSource {
    fn ploidy(&self) -> usize {
        self.ploidy
    }

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
