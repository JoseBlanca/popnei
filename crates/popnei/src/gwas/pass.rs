//! The study itself: the fit of the null model and the one pass over the
//! blocks that tests every variant against it.
//!
//! [`calc_gwas`] is what a caller of the crate runs. [`TheFittedModel`] is
//! the model it fitted, one variant for each model popnei has written, and
//! the pass is the same for all of them: one call for each block, one
//! answer for each variant of it that has variance.

use crate::block::{Block, BlockReader, Reblock};
use crate::error::{Error, Result};
use crate::variant::{ChromTable, Needs};

use super::dosages::{BlockOfThePass, GwasDosages};
use super::linear::LinearModel;
use super::linear_mixed::LinearMixedModel;
use super::logistic::LogisticModel;
use super::logistic_mixed::LogisticMixedModel;
use super::result::{Answers, GrammarGammaApprox, Gwas, NullModel};
use super::study::{
    Design, GwasInput, GwasModel, TestType, refuse_a_kinship_that_is_not_of_the_individuals,
    the_model_and_the_test,
};

/// The null model a study has fitted, which every variant is then tested
/// against.
///
/// The four models keep different things, the thin QR of the design and
/// the residuals of the trait for one, the projection matrix of the
/// covariance for the two mixed ones and the weights of the fitted chances
/// for the plain logistic one, and the pass over the blocks is the same for
/// all of them: one call for each block, one answer for each variant of it
/// that has variance.
#[expect(
    clippy::large_enum_variant,
    reason = "the logistic model carries the buffers its Wald test fits one variant in, \
              which make it about 300 bytes larger than the other three; one of these is \
              made for a study and lives until its pass is over, so that is 300 bytes \
              once, where boxing it would put an allocation and a dereference between \
              every block and the model it is tested against"
)]
enum TheFittedModel {
    /// The linear model, a continuous trait with no kinship.
    Linear(LinearModel),
    /// The linear mixed model, a continuous trait with a kinship.
    Mixed(LinearMixedModel),
    /// The logistic model, a binomial trait with no kinship.
    Logistic(LogisticModel),
    /// The logistic mixed model, a binomial trait with a kinship.
    LogisticMixed(LogisticMixedModel),
}

impl TheFittedModel {
    /// The null model of the result, which carries what the fit gives
    /// besides the rows of the variants.
    #[must_use]
    fn null_model(&self, test: TestType) -> NullModel {
        match self {
            TheFittedModel::Linear(fitted) => fitted.null_model(test),
            TheFittedModel::Mixed(fitted) => fitted.null_model(test),
            TheFittedModel::Logistic(fitted) => fitted.null_model(test),
            TheFittedModel::LogisticMixed(fitted) => fitted.null_model(test),
        }
    }

    /// The test of every variant of a block that has variance among the
    /// tested individuals, in the order of the block.
    ///
    /// A linear model makes the one test it has, the t test of the variant
    /// against what the design left of the trait, and a linear mixed model
    /// and a logistic model make whichever of the Wald test and the score
    /// test the study asked for. A logistic model is not told which: it
    /// was fitted with the test the study asked for and holds the buffers
    /// of that one. What it is given instead is `design`, the design it
    /// was fitted on, which its score test takes each variant through. A
    /// logistic mixed model is not told either, and for the other reason:
    /// the score test is the only test it has.
    ///
    /// # Errors
    ///
    /// Whatever the model's own test of a block fails with.
    fn test_the_block(
        &mut self,
        dosages: &GwasDosages,
        test: TestType,
        design: &Design<'_>,
    ) -> Result<Answers<'_>> {
        match self {
            TheFittedModel::Linear(fitted) => fitted.test_the_block(dosages),
            TheFittedModel::Mixed(fitted) => fitted.test_the_block(dosages, test),
            TheFittedModel::Logistic(fitted) => fitted.test_the_block(dosages, design),
            TheFittedModel::LogisticMixed(fitted) => fitted.test_the_block(dosages),
        }
    }

    /// Estimates the factor of the GRAMMAR-Gamma approximation from
    /// `dosages`, the first block of the second pass, so that every variant
    /// of the pass that follows takes the approximate denominator.
    ///
    /// The two mixed models are the ones that have a denominator to
    /// approximate, and the other two never arrive here: a study with no
    /// kinship that asks for the approximation is refused by
    /// [`calc_gwas`] before any model is fitted, and a study with one fits
    /// a mixed model.
    ///
    /// # Errors
    ///
    /// Whatever the model's own estimate of the factor fails with, and
    /// [`Error::GwasGrammarGammaWithoutAKinship`] for the two models that
    /// have no projection matrix.
    fn approximate_the_denominator(&mut self, dosages: &GwasDosages) -> Result<()> {
        match self {
            TheFittedModel::Linear(_) | TheFittedModel::Logistic(_) => {
                Err(Error::GwasGrammarGammaWithoutAKinship)
            }
            TheFittedModel::Mixed(fitted) => fitted.approximate_the_denominator(dosages),
            TheFittedModel::LogisticMixed(fitted) => fitted.approximate_the_denominator(dosages),
        }
    }
}

/// The association study of the variants of a reader against a trait of
/// the individuals: one row for each variant, with the frequency of the
/// alleles that are not the major one, the effect of one more copy of one
/// of them on the trait, how uncertain that effect is and the p-value of
/// the test that it is 0.
///
/// `input` says which individuals are tested, with their trait and the
/// design the model is fitted on, and "The Rust interface" of
/// `docs/specs/gwas.md` lays out what each of its fields holds. The trait
/// and the kinship choose the model, and all four are written: the linear
/// model, a continuous trait without a kinship, the linear mixed model, a
/// continuous trait with one, the logistic model, a binomial trait without
/// one, with both of its tests, and the logistic mixed model, a binomial
/// trait with one, whose only test is the score test.
///
/// The null model is fitted before the first block is read, from the
/// trait, the design and the kinship alone, and then one pass over the
/// blocks tests every variant against it. The pass borrows the reader
/// and does not take it, so that whoever built the chain of filters reads
/// their counts when it returns; it asks for the genotypes and for the
/// chromosome, the position and the id of each variant, and puts
/// [`Reblock`] before it, since a filter leaves blocks of uneven size and
/// the test of a block is matrix work.
///
/// `gamma_pass` is the second pass over the same variants that the
/// GRAMMAR-Gamma approximation reads its first block of, and it is read
/// only when `use_grammar_gamma_approx` is true. The approximation stands
/// in for the denominator of a mixed model's test, `x' p x`, with one
/// factor times the squared length of the variant's centered dosages, and
/// that factor is the mean over the first 100 variants of that block which
/// vary of the exact denominator divided by the approximate one. So the
/// second pass is read once, after the null model is fitted and before the
/// pass that tests the variants, and it asks its reader for the genotypes
/// alone. A study that asks for the approximation with no kinship is
/// refused, since there is no such denominator to approximate, and so is
/// one that asks for it and gives no second pass.
///
/// # Errors
///
/// [`Error::GwasFitDidNotSettle`] when the
/// null model of a logistic one, plain or mixed, was still moving after the
/// rounds it is given, [`Error::GwasKinshipNotACovariance`] when the
/// covariance of the working trait of a logistic mixed model could not be
/// factored,
/// [`Error::GwasGrammarGammaWithoutAKinship`] when the approximation was
/// asked for by a study with no kinship,
/// [`Error::GwasGrammarGammaWithoutASecondPass`] when it was asked for and
/// `gamma_pass` is `None`,
/// [`Error::GwasGrammarGammaWithoutAVariantThatVaries`] when no variant of
/// the first block of that pass varies among the tested individuals, and
/// [`Error::GwasGrammarGammaFactorNotAboveZero`] when the factor those
/// variants gave is not a finite number above 0.
/// [`Error::GwasInputOfAnotherSize`] when a kinship does not hold one row
/// and one column for each tested individual, and
/// [`Error::GwasKinshipValueNotFinite`] when one of its values is not a
/// finite number, both checked before any model is fitted. What
/// [`the_model_and_the_test`] refuses of the trait, the kinship and the
/// test, and what [`Design::of_the_study`] refuses of the individuals, the
/// phenotype and the design. [`Error::GwasLinalg`] when the fit or the
/// test of a block could not be done. [`Error::PassGaveNoVariant`] when
/// the reader gives no variant, with what each filter of the pass was
/// given and kept, which is pyNei's refusal of a study with nothing to
/// test. [`Error::FieldsNotInTheBlock`] when a block holds no genotypes,
/// or when it has not a column that the first block of the pass had. What
/// [`GwasDosages::read_the_block`] refuses of a block and of its variants,
/// [`Error::GwasVariantsTooLarge`] when the variants of the study are more
/// than a `usize` counts, and whatever the reader fails with.
pub fn calc_gwas<R1: BlockReader, R2: BlockReader>(
    reader: &mut R1,
    gamma_pass: Option<&mut R2>,
    input: &GwasInput<'_>,
) -> Result<Gwas> {
    let (model, test) = the_model_and_the_test(input)?;
    // A study with no kinship has no projection matrix and so no
    // denominator to approximate, and it is refused before the design is
    // built or anything is read.
    if input.use_grammar_gamma_approx && input.kinship.is_none() {
        return Err(Error::GwasGrammarGammaWithoutAKinship);
    }
    if let Some(kinship) = input.kinship {
        refuse_a_kinship_that_is_not_of_the_individuals(kinship, input.individuals.len())?;
    }
    let ploidy = reader.ploidy();
    let design = Design::of_the_study(input, reader.individuals().len())?;
    let mut fitted = match model {
        GwasModel::Lm => {
            TheFittedModel::Linear(LinearModel::of_the_study(input.phenotype, &design)?)
        }
        GwasModel::Lmm => match input.kinship {
            Some(kinship) => TheFittedModel::Mixed(LinearMixedModel::of_the_study(
                input.phenotype,
                &design,
                kinship,
            )?),
            // `the_model_and_the_test` chooses the linear mixed model only
            // for a study that brought a kinship, so a study with none is
            // a linear model and never arrives here.
            None => return Err(Error::GwasModelNotBuilt { model }),
        },
        GwasModel::Glm => {
            TheFittedModel::Logistic(LogisticModel::of_the_study(input.phenotype, &design, test)?)
        }
        GwasModel::Glmm => match input.kinship {
            Some(kinship) => TheFittedModel::LogisticMixed(LogisticMixedModel::of_the_study(
                input.phenotype,
                &design,
                kinship,
            )?),
            // `the_model_and_the_test` chooses the logistic mixed model
            // only for a study that brought a kinship, so a study with none
            // is a logistic model and never arrives here.
            None => return Err(Error::GwasModelNotBuilt { model }),
        },
    };
    let mut dosages = GwasDosages::of_a_study();
    let approximation = match input.use_grammar_gamma_approx {
        false => GrammarGammaApprox::NotUsed,
        true => {
            the_factor_of_the_approximation(
                &mut fitted,
                gamma_pass,
                &design,
                ploidy,
                &mut dosages,
            )?;
            GrammarGammaApprox::Used
        }
    };
    let mut result =
        Gwas::of_the_null_model(fitted.null_model(test), approximation, ChromTable::new());
    // The genotypes and the three columns of the result are what this
    // reads, so a reader over a file leaves the other columns of a variant
    // unparsed.
    reader.set_needs(Needs::GTS | Needs::CHROM_POS | Needs::ID);
    let mut blocks = Reblock::new(reader, None)?;
    let mut first_var = 0_usize;
    while let Some(mut block) = blocks.next_block()? {
        dosages.read_the_block(&mut block, &design, BlockOfThePass { ploidy, first_var })?;
        let answers = fitted.test_the_block(&dosages, test, &design)?;
        result.add_the_block(&dosages, answers)?;
        the_columns_of_the_block(&mut result, &block, first_var)?;
        first_var = first_var
            .checked_add(block.num_vars)
            .ok_or(Error::GwasVariantsTooLarge)?;
    }
    // The names of the chromosomes are taken when the pass is over and not
    // before it: a reader over a file interns the name of a variant as it
    // reads the variant, so the table of a reader that has read nothing is
    // empty, and the numbers of `chroms` would stand for no name.
    result.chrom_table = blocks.chroms().clone();
    if result.num_vars == 0 {
        let filters = blocks.filtering_stats();
        return Err(Error::PassGaveNoVariant {
            // The filter nearest the source was given what the source
            // gave; with no filter the pass gave what the source gave,
            // which is nothing.
            num_vars_of_the_source: filters.last().map_or(0, |(_, stats)| stats.vars_processed),
            filters,
        });
    }
    Ok(result)
}

/// Reads the first block of the second pass and estimates the factor of
/// the GRAMMAR-Gamma approximation from it, so that every variant the pass
/// that follows tests takes the approximate denominator.
///
/// `gamma_pass` reads the same variants as the pass that tests them, and
/// only its first block is read. It is asked for the genotypes alone,
/// where the pass that tests the variants also asks for the chromosome,
/// the position and the id, and [`Reblock`] goes in front of it for the
/// same reason: a filter leaves blocks of uneven size, and which variants
/// the factor comes from would otherwise depend on what the source gave.
/// `ploidy` is the ploidy of that other pass, so a second pass over
/// another dataset is refused rather than read at a width of its own.
///
/// `dosages` is the buffer the pass that follows reads its blocks into,
/// which is used here so that the second pass asks the machine for nothing
/// the first one will not use again.
///
/// # Errors
///
/// [`Error::GwasGrammarGammaWithoutASecondPass`] when `gamma_pass` is
/// `None`, and [`Error::PassGaveNoVariant`] when it gives no block, with
/// what each of its filters was given and kept. What
/// [`GwasDosages::read_the_block`] refuses of that block, which is where a
/// second pass over other individuals or another ploidy is refused, and
/// what the model's own estimate of the factor refuses of the variants
/// that vary in it.
fn the_factor_of_the_approximation<R: BlockReader>(
    fitted: &mut TheFittedModel,
    gamma_pass: Option<&mut R>,
    design: &Design<'_>,
    ploidy: usize,
    dosages: &mut GwasDosages,
) -> Result<()> {
    let Some(gamma_pass) = gamma_pass else {
        return Err(Error::GwasGrammarGammaWithoutASecondPass);
    };
    gamma_pass.set_needs(Needs::GTS);
    let mut blocks = Reblock::new(gamma_pass, None)?;
    let Some(mut block) = blocks.next_block()? else {
        let filters = blocks.filtering_stats();
        return Err(Error::PassGaveNoVariant {
            num_vars_of_the_source: filters.last().map_or(0, |(_, stats)| stats.vars_processed),
            filters,
        });
    };
    dosages.read_the_block(
        &mut block,
        design,
        BlockOfThePass {
            ploidy,
            first_var: 0,
        },
    )?;
    fitted.approximate_the_denominator(dosages)
}

/// Adds the chromosome, the position and the id of the variants of a block
/// to the ones the blocks before it gave.
///
/// Which of the three the result has is the first block's: a source that
/// carries none of them gives a result with none, and `docs/specs/gwas.md`
/// says that a VCF and a vars file both carry the chromosome and the
/// position, so in practice only the id is ever missing. A later block
/// that has not one of them is refused, since its variants would take the
/// values of the variants after them.
///
/// A later block that holds a column the first had not is not read, and
/// nothing of the study can give one: the fields are asked for once before
/// the pass, and [`Reblock`] gives a block whose columns differ from the
/// one before it on its own rather than joining the two.
///
/// # Errors
///
/// [`Error::FieldsNotInTheBlock`] when a block has not a column that the
/// first block of the pass had.
fn the_columns_of_the_block(result: &mut Gwas, block: &Block, first_var: usize) -> Result<()> {
    if first_var == 0 {
        result.chroms = block.chrom.as_ref().map(|_| Vec::new());
        result.poss = block.pos.as_ref().map(|_| Vec::new());
        result.ids = block.id.as_ref().map(|_| Vec::new());
    }
    let mut missing = Needs::empty();
    match (result.chroms.as_mut(), block.chrom.as_ref()) {
        (Some(of_the_pass), Some(of_the_block)) => of_the_pass.extend_from_slice(of_the_block),
        (Some(_), None) => missing |= Needs::CHROM_POS,
        (None, _) => {}
    }
    match (result.poss.as_mut(), block.pos.as_ref()) {
        (Some(of_the_pass), Some(of_the_block)) => of_the_pass.extend_from_slice(of_the_block),
        (Some(_), None) => missing |= Needs::CHROM_POS,
        (None, _) => {}
    }
    match (result.ids.as_mut(), block.id.as_ref()) {
        (Some(of_the_pass), Some(of_the_block)) => of_the_pass.extend_from_slice(of_the_block),
        (Some(_), None) => missing |= Needs::ID,
        (None, _) => {}
    }
    if !missing.is_empty() {
        return Err(Error::FieldsNotInTheBlock { fields: missing });
    }
    Ok(())
}
