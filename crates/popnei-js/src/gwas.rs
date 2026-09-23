//! What a TypeScript user reaches through `calcGwas`: which of the variants
//! of a source are associated with a trait they measured, with the effect of
//! one more copy of a non major allele, how uncertain that effect is and the
//! p-value of the test that the effect is 0.
//!
//! The calculation is the core's, `popnei::gwas::calc_gwas`, and what this
//! module does is the translation that section 11 of
//! `docs/architecture.md` leaves to a binding crate. It builds the chain of
//! readers of the pass from the steps of the `Variants`, turns the name a
//! user wrote for the trait into the value of the core that stands for it,
//! turns the number of the chromosome of each variant into the name the
//! reader gave it and its position into the float64 a number of JavaScript
//! is, and keeps that chain while the calculation runs so that the counts of
//! its filters can be read when it returns.
//!
//! Which individuals are tested, and their phenotype and their design, are
//! the package's: "Which individuals are tested, and the design" of
//! `docs/specs/gwas.md` puts the refusals about a name and about a table in
//! the layer that holds the names, and the core is given the positions of
//! the tested individuals with the phenotype and the design already built.
//! So what crosses here is three arrays that are read together row by row,
//! and the core refuses an order that is not the source's, a phenotype that
//! is not a number and a design whose columns are not independent.
//!
//! [`GwasOfVariants`] is the result on its way out. It lives in the memory
//! of wasm, which the garbage collector of JavaScript does not see, so the
//! package frees it as soon as its arrays are read, and each of them is
//! moved out of it as it is read and not cloned: the four columns of the
//! result are 32 bytes for each variant, 32 MB for a million of them.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::block::BlockReader;
use popnei::gwas::{Gwas, GwasInput, TestType, TraitType, calc_gwas};

use crate::errors::JsPopneiError;
use crate::source::{OpenSource, PassCounts, positions_of};
use crate::steps::{Steps, chain_of};

/// What a study is asked for, as the package checked it: the individuals to
/// test, their trait, their design and how a variant of more than two
/// alleles is read.
///
/// The three arrays are read together row by row, so they hold the tested
/// individuals in one order, the one the source has them in.
pub(crate) struct ArgumentsOfTheStudy {
    /// The position of each tested individual among the individuals the pass
    /// gives, from 0 and rising.
    pub(crate) individuals: Vec<u32>,
    /// The trait of each of them, in that order.
    pub(crate) phenotype: Vec<f64>,
    /// The design, one row of `num_coefs` values for each of them, the
    /// column of ones of the intercept first.
    pub(crate) design: Vec<f64>,
    /// How many columns the design has: the intercept and one for each
    /// covariate.
    pub(crate) num_coefs: usize,
    /// What was measured, one of the names of `popnei::gwas::TraitType`.
    pub(crate) trait_name: String,
    /// Which test is made of every variant, one of the names of
    /// `popnei::gwas::TestType`, and `None` for the default of the model.
    pub(crate) test_name: Option<String>,
    /// Whether a variant of more than two alleles among its called genotypes
    /// is read with every allele that is not the major one counting the
    /// same.
    pub(crate) transform_to_biallelic: bool,
}

/// What a study gives TypeScript: the null model every variant was tested
/// against, the four columns of the variants, the three columns that name
/// them, and the counts of the pass.
///
/// Every array leaves the memory of wasm the first time it is asked for, and
/// the call after that gives nothing: the package reads each of them once,
/// into the object a user holds, and frees this.
#[wasm_bindgen]
pub struct GwasOfVariants {
    model: String,
    test: String,
    /// One effect for each column of the design, and `None` once they were
    /// given to JavaScript.
    covariate_effects: Option<Vec<f64>>,
    residual_variance: Option<f64>,
    genetic_variance: Option<f64>,
    heritability: Option<f64>,
    num_individuals: usize,
    allele_freq: Option<Vec<f64>>,
    beta: Option<Vec<f64>>,
    se: Option<Vec<f64>>,
    p_value: Option<Vec<f64>>,
    used_grammar_gamma_approx: bool,
    /// The name of the chromosome of each variant, and `None` when the
    /// source has no such column or when they were read already.
    chroms: Option<Vec<String>>,
    poss: Option<Vec<f64>>,
    ids: Option<Vec<String>>,
    /// The counts of the pass, which the package turns into the `passStats`
    /// of the result.
    counts: PassCounts,
}

#[wasm_bindgen]
impl GwasOfVariants {
    /// Which of the four models was fitted: `"lm"`, `"lmm"`, `"glm"` or
    /// `"glmm"`.
    #[must_use]
    pub fn model(&self) -> String {
        self.model.clone()
    }

    /// Which test was made of every variant, `"wald"` or `"score"`.
    #[must_use]
    pub fn test(&self) -> String {
        self.test.clone()
    }

    /// The effect of each column of the design, the intercept first, or
    /// `undefined` when they were read already.
    pub fn covariate_effects(&mut self) -> Option<Vec<f64>> {
        self.covariate_effects.take()
    }

    /// What the null model left unexplained, and `undefined` for a binomial
    /// trait, whose variance is decided by its mean.
    #[must_use]
    pub fn residual_variance(&self) -> Option<f64> {
        self.residual_variance
    }

    /// The variance of the random effect of the kinship, and `undefined`
    /// without a kinship.
    #[must_use]
    pub fn genetic_variance(&self) -> Option<f64> {
        self.genetic_variance
    }

    /// The genetic variance over the sum of the two, only for the linear
    /// mixed model.
    #[must_use]
    pub fn heritability(&self) -> Option<f64> {
        self.heritability
    }

    /// How many individuals the study tested.
    #[must_use]
    pub fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// The frequency of the alleles that are not the major one over the
    /// tested individuals, one for each variant, or `undefined` when they
    /// were read already.
    pub fn allele_freq(&mut self) -> Option<Vec<f64>> {
        self.allele_freq.take()
    }

    /// The effect of one more copy of a non major allele, NaN for a variant
    /// that has no answer, or `undefined` when they were read already.
    pub fn beta(&mut self) -> Option<Vec<f64>> {
        self.beta.take()
    }

    /// The standard error of that effect, or `undefined` when they were read
    /// already.
    pub fn se(&mut self) -> Option<Vec<f64>> {
        self.se.take()
    }

    /// The p-value of the test that the effect is 0, or `undefined` when
    /// they were read already.
    pub fn p_value(&mut self) -> Option<Vec<f64>> {
        self.p_value.take()
    }

    /// Whether the GRAMMAR-Gamma approximation was used.
    #[must_use]
    pub fn used_grammar_gamma_approx(&self) -> bool {
        self.used_grammar_gamma_approx
    }

    /// The name of the chromosome of each variant, or `undefined` when the
    /// source has no such column or when they were read already.
    pub fn chroms(&mut self) -> Option<Vec<String>> {
        self.chroms.take()
    }

    /// The position of each variant, 1 based as in a VCF, or `undefined`
    /// when the source has no such column or when they were read already.
    pub fn poss(&mut self) -> Option<Vec<f64>> {
        self.poss.take()
    }

    /// The id of each variant, or `undefined` when the source has no such
    /// column or when they were read already.
    pub fn ids(&mut self) -> Option<Vec<String>> {
        self.ids.take()
    }

    /// How many variants the pass gave and what each filter of it was given
    /// and kept.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.counts.clone()
    }
}

/// The study of the variants of `source` that the steps of `steps` keep,
/// over the individuals, the trait and the design of `study`.
///
/// The null model is fitted before the pass, from the trait and the design
/// alone, and then every variant is tested against what that model left
/// unexplained. A variant whose dosages are all the same among the tested
/// individuals has no variance and cannot be tested: its row is in the
/// result with its frequency, and its effect, its standard error and its
/// p-value are NaN.
///
/// The chain of readers of the pass stays here, lent to the core, so that
/// the counts of its filters can be read when the calculation returns. The
/// source is asked for no size of block: the core puts a `reblock` over the
/// reader and chooses the size there, since the test of a block is matrix
/// work and a filter leaves blocks of uneven size.
///
/// # Errors
///
/// When the name of the trait is of neither of the two and when the name of
/// the test is of neither; when the study needs one of the three models that
/// are not written, the two logistic ones and the linear mixed one; when the
/// score test is asked of a linear model, whose only test is the t test of
/// the effect it fitted; when the individuals to test are not in the order
/// the source has them, one of them is there twice, one of them is not in
/// the source, or they are fewer than the columns of the design plus two;
/// when a value of the phenotype or of the design is not a finite number;
/// when the columns of the design are not independent; when the source
/// cannot be read, a wrong line of a VCF among the causes; when a variant
/// has more than two alleles among its called genotypes and
/// `transform_to_biallelic` is false; when the pass gives no variant; when
/// the linear algebra of a fit or of a test could not be done; and when a
/// position of a variant is above the last whole number JavaScript holds.
pub(crate) fn gwas_of_the_variants(
    source: &dyn OpenSource,
    study: &ArgumentsOfTheStudy,
    steps: Steps,
) -> Result<GwasOfVariants, JsPopneiError> {
    // The two names of a trait and the two of a test are the core's, which
    // is where a name that is of neither is refused: one list of them
    // serves both packages, and the message a user reads is the same in
    // each. Which tests the model of the study has is the core's too, and
    // it is `the_model_and_the_test` that answers it.
    let trait_type = TraitType::of_name(&study.trait_name)?;
    let test = study
        .test_name
        .as_deref()
        .map(TestType::of_name)
        .transpose()?;
    let individuals = the_positions(&study.individuals)?;
    let input = GwasInput {
        phenotype: &study.phenotype,
        trait_type,
        design: &study.design,
        num_coefs: study.num_coefs,
        // The kinship and the GRAMMAR-Gamma approximation reach the core
        // with the linear mixed model, which is being written. The test
        // does not wait for it: the Wald test is the linear model's own,
        // and the core is what refuses the score test of it.
        kinship: None,
        test,
        use_grammar_gamma_approx: false,
        individuals: &individuals,
        transform_to_biallelic: study.transform_to_biallelic,
    };
    let mut chain = chain_of(source.reader(None)?, steps.steps())?;
    let result = calc_gwas(&mut chain, None::<&mut Box<dyn BlockReader>>, &input)?;
    let counted = u64::try_from(result.num_vars).map_err(|_| {
        JsPopneiError::Broken(format!(
            "the study read {num_vars} variants, more than the count of a pass holds",
            num_vars = result.num_vars
        ))
    })?;
    let counts = PassCounts::of(counted, &chain.filtering_stats());
    let chroms = the_names_of_the_chromosomes(&result)?;
    let poss = result.poss.as_deref().map(positions_of).transpose()?;
    Ok(GwasOfVariants {
        model: result.null_model.model.name().to_owned(),
        test: result.null_model.test.name().to_owned(),
        covariate_effects: Some(result.null_model.covariate_effects),
        residual_variance: result.null_model.residual_variance,
        genetic_variance: result.null_model.genetic_variance,
        heritability: result.null_model.heritability,
        num_individuals: result.null_model.num_individuals,
        allele_freq: Some(result.allele_freq),
        beta: Some(result.beta),
        se: Some(result.se),
        p_value: Some(result.p_value),
        used_grammar_gamma_approx: result.used_grammar_gamma_approx,
        chroms,
        poss,
        ids: result.ids,
        counts,
    })
}

/// The positions of the tested individuals as the core counts them.
///
/// They cross as a `Uint32Array`, which holds every position of a panel a
/// browser reads, and the core takes them as the whole numbers this build
/// indexes memory with, 32 bits in WebAssembly.
///
/// # Errors
///
/// When a position is more than this build indexes with, which in a browser
/// is no position at all: a `Uint32Array` and a `usize` are both 32 bits
/// there.
fn the_positions(individuals: &[u32]) -> Result<Vec<usize>, JsPopneiError> {
    individuals
        .iter()
        .map(|position| {
            usize::try_from(*position).map_err(|_| {
                JsPopneiError::Broken(format!(
                    "the individual at the position {position} was asked to be tested, and \
                     this build counts the individuals of a source to {largest}",
                    largest = usize::MAX
                ))
            })
        })
        .collect()
}

/// The name of the chromosome of each variant of `result`, read through the
/// table of names the core cloned from the reader of the pass, and `None`
/// when the source has no such column.
///
/// # Errors
///
/// When a number of a variant is not in that table, which cannot happen
/// unless this crate or the core has a defect.
fn the_names_of_the_chromosomes(result: &Gwas) -> Result<Option<Vec<String>>, JsPopneiError> {
    let table = &result.chrom_table;
    result
        .chroms
        .as_deref()
        .map(|chroms| {
            chroms
                .iter()
                .map(|number| {
                    table.name(*number).map(str::to_owned).ok_or_else(|| {
                        JsPopneiError::Broken(format!(
                            "the chromosome number {number} of a variant of the study is \
                             not in the table of the reader that gave it"
                        ))
                    })
                })
                .collect()
        })
        .transpose()
}
