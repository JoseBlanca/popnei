//! What a TypeScript user reaches through `openVcf`: a VCF with its options,
//! read as a source of variants.
//!
//! [`VcfSource`] holds the file, a copy of it in the memory of wasm or a file
//! of the page read one range at a time, and the options it is read with, and
//! it reads the header when it is built, so bytes that are not a VCF fail at
//! `openVcf`. Every pass over it reads the file again from its start and goes
//! through the `Blocks` of `source.rs`, the one a vars file goes through, and
//! `write_vars` there is what writes its variants into a vars file.

use js_sys::Function;
use wasm_bindgen::prelude::wasm_bindgen;
use web_sys::Blob;

use popnei::block::BlockReader;
use popnei::io::vcf::{VcfOptions, VcfReader};

use crate::dists::{KosmanDistances, kosman_dists_of};
use crate::diversity::{ArgumentsOfTheDiversity, PopDiversityOfAPass, pop_diversity_of};
use crate::errors::JsPopneiError;
use crate::gwas::{ArgumentsOfTheStudy, GwasOfVariants, gwas_of_the_variants};
use crate::kinship::{KinshipOfVariants, kinship_of_the_variants};
use crate::ld::{ArgumentsOfTheBins, LdAndDistOfAPass, R2Matrix, ld_and_dist_of, r2_matrix_of};
use crate::pca::{PcaOfVariants, pca_of_the_variants};
use crate::pop_dists::{ArgumentsOfTheDists, PopDistsOfAPass, pop_dists_of};
use crate::source::{
    Blocks, Consumer, OpenSource, RunOfAConsumer, TheFileOfASource, VarsFile, blocks_of,
    bytes_of_a_vars_file, starts_a_run_of, tells_the_progress, the_bytes_of_a_new_source,
    the_file_of_a_new_source, the_source_was_freed,
};
use crate::stats::{
    ArgumentsOfThePass, PerIndividualStats, PerVarDistribs, per_individual_stats_of,
    per_var_distribs_of,
};
use crate::steps::Steps;

/// A VCF that was opened: where its file is, the options it is read with,
/// and the individuals its header named.
#[wasm_bindgen]
pub struct VcfSource {
    /// Where the file is, a copy of the whole of it in the memory of wasm or
    /// a file of the page read one range at a time, which every pass over the
    /// source reads again from its first byte.
    file: TheFileOfASource,
    options: VcfOptions,
    individuals: Vec<String>,
    /// The number of what this source keeps in JavaScript, the file it reads
    /// the ranges from and the function the page is told the progress with,
    /// which `free()` gives back.
    in_javascript: u32,
}

impl Drop for VcfSource {
    fn drop(&mut self) {
        the_source_was_freed(self.in_javascript);
    }
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

    /// The function the page is told how far every pass over this source has
    /// got with, `told`, and nothing to take the one that was set off.
    ///
    /// It holds until it is set again, and setting it changes nothing about
    /// the variants a pass gives. The reads of `openVcf`, which are made
    /// before there is a `Variants` to set a function on, are told to nobody.
    pub fn on_progress(&self, told: Option<Function>) {
        tells_the_progress(self.in_javascript, told);
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

    /// The five per variant statistics of one pass over the VCF, through
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
    /// one pass over the VCF, through the steps of `steps`.
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

    /// The kinship of every pair of the individuals of `individuals`, or of
    /// every individual of the pass when it is nothing, over the variants of
    /// the VCF that the steps of `steps` keep.
    ///
    /// # Errors
    ///
    /// When a name of `individuals` is of nobody the pass gives, is there
    /// twice, or the list is empty; when the VCF cannot be read; when a
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

    /// Which of the variants of the VCF that the steps of `steps` keep
    /// are associated with the trait of `phenotype`, over the individuals at
    /// the positions of `individuals` and with `design` as the numbers the
    /// model fits beside each variant.
    ///
    /// The three arrays hold the tested individuals in the order the source
    /// has them, one position, one value of the trait and one row of
    /// `num_coefs` values of the design each, and they are read together row
    /// by row. `kinship` is the relatedness of those individuals, one row
    /// and one column for each of them and row after row, which makes the
    /// study a mixed model, and `undefined` leaves it a model with no
    /// random effect.
    ///
    /// # Errors
    ///
    /// When `trait_name` is of neither trait; when the GRAMMAR-Gamma
    /// approximation is asked for by a study with no kinship, or the first
    /// block its factor would come from leaves no factor above 0; when a
    /// binomial phenotype holds a value that is neither 0 nor 1, or its
    /// null model walks towards an infinite coefficient instead of settling;
    /// when the
    /// individuals are not in the order the source has them, one is there
    /// twice or is not in the source, or they are fewer than the columns of
    /// the design plus two; when a value of the phenotype or of the design is not finite;
    /// when the kinship is not one row and one column for each tested
    /// individual or holds a value that is not finite; when the columns of
    /// the design are not independent; when the VCF cannot be read; when a
    /// variant has more than two alleles among its called genotypes and
    /// `transform_to_biallelic` is false; when the pass gives no variant;
    /// and when the linear algebra could not be done.
    #[expect(
        clippy::too_many_arguments,
        reason = "the study as the package checked it: the tested individuals, their \
                  trait, their design, the kinship of a mixed model and the three \
                  things a user asked for, each array flat, since a table is not one \
                  of the types wasm-bindgen carries"
    )]
    pub fn calc_gwas(
        &self,
        individuals: Vec<u32>,
        phenotype: Vec<f64>,
        design: Vec<f64>,
        num_coefs: usize,
        trait_name: String,
        test_name: Option<String>,
        kinship: Option<Vec<f64>>,
        use_grammar_gamma_approx: bool,
        transform_to_biallelic: bool,
        steps: Steps,
    ) -> Result<GwasOfVariants, JsPopneiError> {
        gwas_of_the_variants(
            self,
            &ArgumentsOfTheStudy {
                individuals,
                phenotype,
                design,
                num_coefs,
                trait_name,
                test_name,
                kinship,
                use_grammar_gamma_approx,
                transform_to_biallelic,
            },
            steps,
        )
    }

    /// The Kosman distance of every pair of individuals over the variants
    /// of the VCF that the steps of `steps` keep, with no distance for a
    /// pair called at fewer than `min_num_vars` variants, and the counts of
    /// the pass that gave them.
    ///
    /// # Errors
    ///
    /// When the pass gives no variant, when the sums of a pair go above
    /// what a `u32` holds, when the memory of the tab does not take the two
    /// counts of every pair, and when the VCF cannot be read.
    pub fn calc_pairwise_kosman_dists(
        &self,
        min_num_vars: u32,
        steps: Steps,
    ) -> Result<KosmanDistances, JsPopneiError> {
        kosman_dists_of(self, min_num_vars, steps)
    }

    /// The r² of every pair of the variants of the VCF that the steps of
    /// `steps` keep, with the chromosome and the position of each of them
    /// and the counts of the pass, over at most `max_num_vars` variants.
    ///
    /// # Errors
    ///
    /// When the pass gives more than `max_num_vars` variants, when the
    /// matrix of that many holds more values than wasm counts, when the pass
    /// gives no variant, when the memory of the tab does not take the
    /// matrix, and when the VCF cannot be read.
    pub fn calc_rogers_huff_r2_matrix(
        &self,
        max_num_vars: usize,
        steps: Steps,
    ) -> Result<R2Matrix, JsPopneiError> {
        r2_matrix_of(self, max_num_vars, steps)
    }

    /// How the r² of a pair of variants falls off with the distance between
    /// them, in bins of distance and for each population that was named,
    /// over one pass over the VCF through the steps of `steps`.
    ///
    /// The arguments are those of `calcLdAndDistPerPop` of
    /// `docs/specs/ld.md`, as the package checked them and flat: the
    /// populations are their names, the names of the individuals of every
    /// one of them one after another, and how many individuals each of them
    /// holds, and no name at all is one population of every individual;
    /// `min_dist` and `max_dist` are the distances in base pairs a pair is
    /// counted at, both included; `num_bins` is how many bins of equal
    /// width they are cut into; and `max_allowed_maf` is the largest major
    /// allele frequency a variant has in a population and is still counted
    /// there.
    ///
    /// # Errors
    ///
    /// Those of [`ld_and_dist_of`]: a `min_dist` above `max_dist`, a
    /// `num_bins` of 0, a `max_allowed_maf` that is not a number from 0 to
    /// 1, a population that names an individual the pass does not give,
    /// names one twice or names none, a pass that gives no variant, bins or
    /// a window the memory of the tab does not take, a bin of more pairs
    /// than a number of JavaScript counts one by one, and a VCF that cannot
    /// be read.
    #[expect(
        clippy::too_many_arguments,
        reason = "the arguments of `calcLdAndDistPerPop` of `docs/specs/ld.md`, each \
                  one as the package checked it, and the populations flat: an array of \
                  arrays is not one of the types wasm-bindgen carries"
    )]
    pub fn calc_ld_and_dist_per_pop(
        &self,
        steps: Steps,
        pop_names: Option<Vec<String>>,
        pop_individuals: Vec<String>,
        num_individuals_per_pop: Vec<u32>,
        min_dist: f64,
        max_dist: f64,
        num_bins: usize,
        max_allowed_maf: f64,
    ) -> Result<LdAndDistOfAPass, JsPopneiError> {
        ld_and_dist_of(
            self,
            &steps,
            &ArgumentsOfTheBins {
                pop_names,
                pop_individuals,
                num_individuals_per_pop,
                min_dist,
                max_dist,
                num_bins,
                max_allowed_maf,
            },
        )
    }

    /// Every measure of `measures` for every pair of the populations that
    /// were named, over one pass over the VCF through the steps of
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

    /// How much variety each population of `pop_names` holds, over one pass
    /// over the file through the steps of `steps`.
    ///
    /// The arguments are those of `calcPopDiversity` of
    /// `docs/specs/diversity.md`, as the package checked them and flat:
    /// `stats` holds the name of each statistic to calculate; the
    /// populations are their names, the names of the individuals of every
    /// one of them one after another, and how many individuals each of them
    /// holds, and `pop_names` is nothing when the user named no population,
    /// which is one population of every individual of the pass;
    /// `num_called_alleles` is how many called alleles every population is
    /// brought down to, and nothing for a pass that takes no draw; and
    /// `min_num_individuals` is how many called genotypes a population needs
    /// at a variant for the variant to count for it.
    ///
    /// # Errors
    ///
    /// Those of [`pop_diversity_of`]: a name that is of none of the five
    /// statistics, the folded spectrum asked for with no draw, a draw of
    /// fewer than two alleles, a population that names an individual the
    /// pass does not give, names one twice or names none, `pops` with no
    /// population, a source that cannot be read, a pass that gives no
    /// variant, a variant of more alleles than a count of them holds, and a
    /// count above what a JavaScript array of counts holds.
    #[expect(
        clippy::too_many_arguments,
        reason = "the arguments of `calcPopDiversity` of `docs/specs/diversity.md`, \
                  each one as the package checked it, and the populations flat: an \
                  array of arrays is not one of the types wasm-bindgen carries"
    )]
    pub fn calc_pop_diversity(
        &self,
        steps: Steps,
        stats: Vec<String>,
        pop_names: Option<Vec<String>>,
        pop_individuals: Vec<String>,
        num_individuals_per_pop: Vec<u32>,
        num_called_alleles: Option<u32>,
        min_num_individuals: u32,
    ) -> Result<PopDiversityOfAPass, JsPopneiError> {
        pop_diversity_of(
            self,
            &steps,
            &ArgumentsOfTheDiversity {
                stats,
                pop_names,
                pop_individuals,
                num_individuals_per_pop,
                num_called_alleles,
                min_num_individuals,
            },
        )
    }
}

impl OpenSource for VcfSource {
    fn ploidy(&self) -> usize {
        self.options.ploidy
    }

    fn starts_a_run(&self, consumer: &Consumer) -> RunOfAConsumer {
        starts_a_run_of(self.in_javascript, consumer)
    }

    fn reader(
        &self,
        run: &RunOfAConsumer,
        num_vars_per_block: Option<usize>,
    ) -> Result<Box<dyn BlockReader>, popnei::Error> {
        let options = VcfOptions {
            num_vars_per_block,
            ..self.options
        };
        // The header is read here, which is the first read of the pass and
        // the call that tells the page that it has read nothing yet.
        Ok(Box::new(VcfReader::new(
            self.file.a_pass_of(self.in_javascript, run)?,
            options,
        )?))
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
    let (file, in_javascript) = the_bytes_of_a_new_source(bytes)?;
    the_vcf_of(
        file,
        in_javascript,
        VcfOptions {
            ploidy,
            only_passed,
            num_vars_per_block: None,
        },
    )
}

/// The VCF in `file`, the file the user picked in the page or a `Blob` an
/// application made itself, plain or gzipped, read with `ploidy` alleles in
/// every genotype and, when `only_passed` is true, without the variants that
/// failed a filter.
///
/// The file stays in the page. Every pass over it asks the browser for one
/// range of a few MiB at a time through `FileReaderSync`, which a browser
/// gives only inside a web worker, so the file is never in the memory of
/// wasm whole and a file larger than that memory is read.
///
/// It reads the header, so the individuals are known when it returns.
///
/// # Errors
///
/// When the browser has no `FileReaderSync`, which is every call outside a
/// web worker; when `Blob.size` is not a whole number of bytes popnei reads a
/// file by; when the file is not a VCF that popnei can read; and when the
/// ploidy is out of the range the core takes.
#[wasm_bindgen]
pub fn open_vcf_of_a_file(
    file: Blob,
    ploidy: usize,
    only_passed: bool,
) -> Result<VcfSource, JsPopneiError> {
    let (file, in_javascript) = the_file_of_a_new_source(file)?;
    the_vcf_of(
        file,
        in_javascript,
        VcfOptions {
            ploidy,
            only_passed,
            num_vars_per_block: None,
        },
    )
}

/// The source of the VCF in `file`, which keeps the file and the function the
/// page is told the progress with in the entry numbered `in_javascript`, read
/// with `options`.
///
/// # Errors
///
/// When the file is not a VCF that popnei can read, and when the ploidy is
/// out of the range the core takes. The entry goes with an open that failed:
/// no `Variants` was made, so no `free()` will come for it.
fn the_vcf_of(
    file: TheFileOfASource,
    in_javascript: u32,
    options: VcfOptions,
) -> Result<VcfSource, JsPopneiError> {
    // The header is read when the reader is built and no variant is.
    // Nothing here asks for a block, so a file whose blocks would need more
    // memory than wasm addresses, a header of 170000 individuals read with
    // the ploidy 255, is opened all the same and its individuals read; the
    // size of its blocks is the user's to choose at `iterBlocks`.
    //
    // The read belongs to no run and is told to nobody: there is no
    // `Variants` yet for an application to have set a function on.
    let opened = file
        .the_opening_pass(in_javascript)
        .and_then(|pass| VcfReader::new(pass, options));
    let reader = match opened {
        Ok(reader) => reader,
        Err(error) => {
            the_source_was_freed(in_javascript);
            return Err(error.into());
        }
    };
    let individuals = reader.individuals().to_vec();
    Ok(VcfSource {
        file,
        options,
        individuals,
        in_javascript,
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
