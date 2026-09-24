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
use super::result::{Answers, GrammarGammaApprox, Gwas, NullModel};
use super::study::{
    Design, GwasInput, GwasModel, TestType, refuse_a_kinship_that_is_not_of_the_individuals,
    the_model_and_the_test,
};

/// The null model a study has fitted, which every variant is then tested
/// against.
///
/// The three models that are written keep different things, the thin QR of
/// the design and the residuals of the trait for one, the projection matrix
/// of the covariance for the second and the weights of the fitted chances
/// for the third, and the pass over the blocks is the same for all of them:
/// one call for each block, one answer for each variant of it that has
/// variance. The logistic mixed model is one more variant here when it is
/// written.
#[expect(
    clippy::large_enum_variant,
    reason = "the logistic model carries the buffers its Wald test fits one variant in, \
              which make it about 300 bytes larger than the other two; one of these is \
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
    /// was fitted on, which its score test takes each variant through.
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
/// and the kinship choose the model, and of the four three are written:
/// the linear model, a continuous trait without a kinship, the linear
/// mixed model, a continuous trait with one, and the logistic model, a
/// binomial trait without one, with both of its tests. The logistic mixed
/// model, a binomial trait with a kinship, is refused with
/// [`Error::GwasModelNotBuilt`] until it is written.
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
/// GRAMMAR-Gamma approximation reads its first block of. Nothing reads it
/// yet: the approximation stands in for the denominator of a mixed model's
/// test, "The GRAMMAR-Gamma approximation" of `docs/specs/gwas.md` is not
/// written, and a study that asks for it is refused, without a kinship
/// because there is no such denominator to approximate and with one
/// because popnei cannot approximate it yet.
///
/// # Errors
///
/// [`Error::GwasModelNotBuilt`] when the study asks for the logistic mixed
/// model, which is not written, [`Error::GwasFitDidNotSettle`] when the
/// null model of a logistic one was still moving after the rounds it is
/// given,
/// [`Error::GwasGrammarGammaWithoutAKinship`] when the approximation was
/// asked for by a study with no kinship and
/// [`Error::GwasGrammarGammaNotBuilt`] when it was asked for by one with a
/// kinship.
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
#[expect(
    unused_variables,
    reason = "`gamma_pass` keeps the name `docs/specs/gwas.md` gives it, since rustdoc \
              prints the names of the arguments; nothing reads it until a mixed model \
              makes the approximation, and the doc comment above says why"
)]
pub fn calc_gwas<R1: BlockReader, R2: BlockReader>(
    reader: &mut R1,
    gamma_pass: Option<&mut R2>,
    input: &GwasInput<'_>,
) -> Result<Gwas> {
    let (model, test) = the_model_and_the_test(input)?;
    if input.use_grammar_gamma_approx {
        // The study is refused whether or not it brought a kinship, and
        // the two refusals say different things: without one there is no
        // denominator to approximate, and with one there is and popnei has
        // not written the approximation of it.
        return match input.kinship {
            None => Err(Error::GwasGrammarGammaWithoutAKinship),
            Some(_) => Err(Error::GwasGrammarGammaNotBuilt),
        };
    }
    if let Some(kinship) = input.kinship {
        refuse_a_kinship_that_is_not_of_the_individuals(kinship, input.individuals.len())?;
    }
    match model {
        GwasModel::Lm | GwasModel::Lmm | GwasModel::Glm => {}
        // The logistic mixed model is the one of the four that is not
        // written, and a binomial trait with a kinship is what asks for
        // it. It is refused here and not at the fit below so that a study
        // popnei cannot run is refused before its design is read.
        GwasModel::Glmm => return Err(Error::GwasModelNotBuilt { model }),
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
        GwasModel::Glmm => return Err(Error::GwasModelNotBuilt { model }),
    };
    let mut result = Gwas::of_the_null_model(
        fitted.null_model(test),
        GrammarGammaApprox::NotUsed,
        ChromTable::new(),
    );
    // The genotypes and the three columns of the result are what this
    // reads, so a reader over a file leaves the other columns of a variant
    // unparsed.
    reader.set_needs(Needs::GTS | Needs::CHROM_POS | Needs::ID);
    let mut blocks = Reblock::new(reader, None)?;
    let mut dosages = GwasDosages::of_a_study();
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
