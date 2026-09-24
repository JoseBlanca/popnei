//! What a study gives back, and the answers a model fills it with.
//!
//! [`Gwas`] is the result: one row for each variant the study was given,
//! and the [`NullModel`] every one of them was tested against.
//! [`Gwas::add_the_block`] fills it block by block with the [`Answers`] of
//! the variants of that block which have variance and the three NaNs of the
//! ones which have none.

use crate::error::{Error, Result};
use crate::variant::ChromTable;

use super::dosages::GwasDosages;
use super::study::{GwasModel, TestType};

/// Whether a study used the GRAMMAR-Gamma approximation, which a mixed
/// model can take to spend one product per variant instead of a fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum GrammarGammaApprox {
    /// The study made the approximation.
    Used,
    /// The study did not, which is what every model without a kinship
    /// does, since there is nothing to approximate.
    NotUsed,
}

/// Whether a study makes the GRAMMAR-Gamma approximation when the user
/// asks for nothing: it does not. It is the default of `calc_gwas` of
/// `pynei/gwas.py`, which popnei keeps, and the reason is in "What it
/// gives" of the approximation in `docs/specs/gwas.md`: it costs accuracy
/// where a panel is strongly structured, so a user asks for it.
pub const DEFAULT_USE_GRAMMAR_GAMMA_APPROX: bool = false;

/// The model a study fitted without any variant in it, which every variant
/// is then tested against.
#[derive(Debug, Clone, PartialEq)]
pub struct NullModel {
    /// Which of the four models it is.
    pub model: GwasModel,
    /// Which test is made of every variant.
    pub test: TestType,
    /// One per column of the design: the intercept first and then one for
    /// each covariate, in the units of the trait for a continuous one and
    /// as a log odds ratio for a binomial one.
    pub covariate_effects: Vec<f64>,
    /// What the model left unexplained, and `None` for a binomial trait,
    /// whose variance is decided by its mean.
    pub residual_variance: Option<f64>,
    /// The variance of the random effect of the kinship, and `None`
    /// without a kinship.
    pub genetic_variance: Option<f64>,
    /// The genetic variance over the sum of the two, only for the linear
    /// mixed model.
    pub heritability: Option<f64>,
    /// How many individuals the study tests.
    pub num_individuals: usize,
}

/// What a model answered for the variants of one block that have variance,
/// in the order of the block: one value each for the effect of one more
/// copy of a non major allele, how uncertain that effect is, and the
/// p-value of the test that the effect is 0.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Answers<'a> {
    /// The effect of one more copy of a non major allele.
    pub beta: &'a [f64],
    /// The standard error of that effect.
    pub se: &'a [f64],
    /// The p-value of the test that the effect is 0.
    pub p_value: &'a [f64],
}

/// What a study gives back: one row for each variant it was given, and the
/// null model every one of them was tested against.
#[derive(Debug, Clone)]
pub struct Gwas {
    /// How many variants the study was given, which is the rows of every
    /// column below.
    pub num_vars: usize,
    /// The model fitted without any variant in it.
    pub null_model: NullModel,
    /// The frequency of the alleles that are not the major one, over the
    /// tested individuals.
    pub allele_freq: Vec<f64>,
    /// The effect of one more copy of a non major allele, in the units of
    /// the trait for a continuous one and as a log odds ratio for a
    /// binomial one.
    pub beta: Vec<f64>,
    /// The standard error of that effect.
    pub se: Vec<f64>,
    /// The p-value of the test that the effect is 0.
    pub p_value: Vec<f64>,
    /// Whether the GRAMMAR-Gamma approximation was used.
    pub used_grammar_gamma_approx: bool,
    /// The interned chromosome of each variant, read through
    /// [`Gwas::chrom_table`], and `None` when the source had no such
    /// column.
    pub chroms: Option<Vec<u32>>,
    /// The names of the chromosomes the numbers of `chroms` stand for.
    pub chrom_table: ChromTable,
    /// The position of each variant, and `None` when the source had no
    /// such column.
    pub poss: Option<Vec<u64>>,
    /// The id of each variant, and `None` when the source had no such
    /// column.
    pub ids: Option<Vec<String>>,
}

impl Gwas {
    /// The result of a study that has fitted its null model and read no
    /// block yet.
    ///
    /// The null model is fitted before the pass over the blocks, from the
    /// trait, the design and the kinship alone, as "How it runs" of
    /// `docs/specs/gwas.md` says, so the result is opened with it and the
    /// blocks add their rows to it. `chrom_table` is the one of the reader
    /// the study reads, which the numbers of `chroms` are read through.
    #[must_use]
    pub(crate) fn of_the_null_model(
        null_model: NullModel,
        grammar_gamma_approx: GrammarGammaApprox,
        chrom_table: ChromTable,
    ) -> Gwas {
        Gwas {
            num_vars: 0,
            null_model,
            allele_freq: Vec::new(),
            beta: Vec::new(),
            se: Vec::new(),
            p_value: Vec::new(),
            used_grammar_gamma_approx: matches!(grammar_gamma_approx, GrammarGammaApprox::Used),
            chroms: None,
            chrom_table,
            poss: None,
            ids: None,
        }
    }

    /// It adds one row for each variant of a block: the frequency of every
    /// one of them, and the answer of the ones that have variance.
    ///
    /// A variant whose dosages are all the same among the tested
    /// individuals has no variance and cannot be tested. Its row is still
    /// here, with its frequency, and `beta`, `se` and `p_value` are NaN,
    /// which is "The variants that have no answer" of
    /// `docs/specs/gwas.md`: a variant with one allele lands there, and so
    /// does one where every tested individual is heterozygous and one with
    /// no called genotype at all. `answers` holds one value for each
    /// variant that does have variance, in the order of the block, which
    /// is the order [`GwasDosages::dosages`] holds their rows in.
    ///
    /// The three columns `chroms`, `poss` and `ids` are not filled here:
    /// they are the block's own columns and belong to the pass that reads
    /// the blocks.
    ///
    /// # Errors
    ///
    /// [`Error::GwasAnswersOfAnotherSize`] when a column of `answers` does
    /// not hold one value for each variant of the block that has variance,
    /// and [`Error::GwasVariantsTooLarge`] when the variants of the study
    /// are more than a `usize` counts.
    pub(crate) fn add_the_block(
        &mut self,
        dosages: &GwasDosages,
        answers: Answers<'_>,
    ) -> Result<()> {
        for (column, num_values) in [
            ("beta", answers.beta.len()),
            ("se", answers.se.len()),
            ("p_value", answers.p_value.len()),
        ] {
            if num_values != dosages.num_with_variance {
                return Err(Error::GwasAnswersOfAnotherSize {
                    column,
                    num_values,
                    num_with_variance: dosages.num_with_variance,
                });
            }
        }
        self.num_vars = self
            .num_vars
            .checked_add(dosages.num_vars)
            .ok_or(Error::GwasVariantsTooLarge)?;
        let mut answered = answers
            .beta
            .iter()
            .zip(answers.se)
            .zip(answers.p_value)
            .map(|((beta, se), p_value)| (*beta, *se, *p_value));
        for (allele_freq, has_variance) in dosages.allele_freq.iter().zip(&dosages.has_variance) {
            self.allele_freq.push(*allele_freq);
            let answer = match *has_variance {
                false => Some((f64::NAN, f64::NAN, f64::NAN)),
                true => answered.next(),
            };
            // The answers hold one value for each variant of the block that
            // has variance, which was counted above, so there is one here
            // for every variant that has any.
            let Some((beta, se, p_value)) = answer else {
                return Err(Error::GwasAnswersOfAnotherSize {
                    column: "beta",
                    num_values: answers.beta.len(),
                    num_with_variance: dosages.num_with_variance,
                });
            };
            self.beta.push(beta);
            self.se.push(se);
            self.p_value.push(p_value);
        }
        Ok(())
    }
}

/// The rows a study gives back, and the variants that have no answer, as
/// "The variants that have no answer" of `docs/specs/gwas.md` states them.
#[cfg(test)]
mod tests {
    use super::{Answers, GrammarGammaApprox, Gwas, NullModel};
    use crate::error::Error;
    use crate::gwas::dosages::GwasDosages;
    use crate::gwas::fixtures::{
        FREQUENCIES_OF_THE_TESTED, OF_EIGHT, TESTED_OF_EIGHT, a_block, a_study, assert_the_values,
        the_design_of, the_first_block_of, the_phenotype_and_the_design_of,
    };
    use crate::gwas::study::the_model_and_the_test;
    use crate::variant::ChromTable;

    /// The dosages of the panel of eight over the four individuals that
    /// have a phenotype, of which the second and the third variants have
    /// variance and the first and the fourth have none.
    fn the_dosages_of_the_tested() -> GwasDosages {
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        let design = the_design_of(&study, 8);
        let mut block = a_block(8, 2, &OF_EIGHT);
        let mut dosages = GwasDosages::of_a_study();
        dosages
            .read_the_block(&mut block, &design, the_first_block_of(2))
            .expect("the dosages of the block");
        dosages
    }

    /// A result opened on the null model of that study, with the model and
    /// the test the study chose.
    fn a_result() -> Gwas {
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        let (model, test) = the_model_and_the_test(&study).expect("the model and the test");
        Gwas::of_the_null_model(
            NullModel {
                model,
                test,
                covariate_effects: vec![3.6666666666666683, 1.0],
                residual_variance: Some(3.333333333333334),
                genetic_variance: None,
                heritability: None,
                num_individuals: 4,
            },
            GrammarGammaApprox::NotUsed,
            ChromTable::new(),
        )
    }

    /// The result carries whether the study made the GRAMMAR-Gamma
    /// approximation, which is the one thing about a study that is not in
    /// its null model, and a study that was asked for nothing makes it
    /// not.
    ///
    /// It is a `bool` in the result, which "The Rust interface" of
    /// `docs/specs/gwas.md` fixes, and the two named values on the way in,
    /// so that a call cannot say which it means with a bare `false`.
    #[test]
    fn the_result_says_whether_the_approximation_was_used() {
        let null_model = a_result().null_model;
        for (asked, used) in [
            (GrammarGammaApprox::Used, true),
            (GrammarGammaApprox::NotUsed, false),
        ] {
            let result = Gwas::of_the_null_model(null_model.clone(), asked, ChromTable::new());
            assert_eq!(result.used_grammar_gamma_approx, used);
            assert_eq!(result.num_vars, 0, "no block has been read");
            assert!(result.allele_freq.is_empty());
        }
        let asked_for_nothing = match super::DEFAULT_USE_GRAMMAR_GAMMA_APPROX {
            true => GrammarGammaApprox::Used,
            false => GrammarGammaApprox::NotUsed,
        };
        let result = Gwas::of_the_null_model(null_model, asked_for_nothing, ChromTable::new());
        assert!(
            !result.used_grammar_gamma_approx,
            "a study that is asked for nothing does not approximate"
        );
    }

    /// A variant with no variance among the tested individuals keeps its
    /// frequency and has NaN for its effect, its uncertainty and its
    /// p-value, and the variants that do have variance take the answers of
    /// the model in the order of the block.
    ///
    /// The fixture is the panel of eight of `fixtures::OF_EIGHT` over its
    /// four tested individuals: its first variant is heterozygous in all
    /// four and its fourth has no genotype in any of them, so the answers
    /// the model gives are two and the rows are four. The second block
    /// added is the same one again, which is what says that a block adds
    /// its rows to the ones that are there and does not replace them.
    #[test]
    fn a_variant_with_no_variance_keeps_its_frequency_and_has_three_nans() {
        let dosages = the_dosages_of_the_tested();
        let mut result = a_result();
        let beta = [1.5, -0.25];
        let se = [0.600925212577332, 0.4];
        let p_value = [0.088004892382756, 0.5];

        result
            .add_the_block(
                &dosages,
                Answers {
                    beta: &beta,
                    se: &se,
                    p_value: &p_value,
                },
            )
            .expect("the answers of the two variants that have variance");

        assert_eq!(result.num_vars, 4);
        assert_the_values(
            &result.allele_freq,
            &FREQUENCIES_OF_THE_TESTED,
            "the frequency of every variant of the block",
        );
        for (column, values) in [
            ("beta", &result.beta),
            ("se", &result.se),
            ("p_value", &result.p_value),
        ] {
            assert!(
                values.first().is_some_and(|value| value.is_nan()),
                "the first variant is heterozygous in all four and has no {column}"
            );
            assert!(
                values.get(3).is_some_and(|value| value.is_nan()),
                "the fourth variant has no genotype among the four and no {column}"
            );
        }
        assert_the_values(
            result.beta.get(1..3).unwrap_or_default(),
            &beta,
            "the effects of the two variants that have variance",
        );
        assert_the_values(
            result.se.get(1..3).unwrap_or_default(),
            &se,
            "their standard errors",
        );
        assert_the_values(
            result.p_value.get(1..3).unwrap_or_default(),
            &p_value,
            "their p-values",
        );

        result
            .add_the_block(
                &dosages,
                Answers {
                    beta: &beta,
                    se: &se,
                    p_value: &p_value,
                },
            )
            .expect("the same block again");

        assert_eq!(result.num_vars, 8, "the rows of the two blocks");
        assert_eq!(result.allele_freq.len(), 8);
        assert_eq!(result.p_value.len(), 8);
    }

    /// A model that answered for another number of variants than the block
    /// has with variance is refused, naming the column and the two counts.
    /// It is a defect of popnei and not of a user, and it would otherwise
    /// put the answer of one variant in the row of another from the first
    /// missing value on.
    #[test]
    fn a_model_that_answered_for_another_number_of_variants_is_refused() {
        let dosages = the_dosages_of_the_tested();
        let mut result = a_result();

        match result.add_the_block(
            &dosages,
            Answers {
                beta: &[1.5],
                se: &[0.6, 0.4],
                p_value: &[0.08, 0.5],
            },
        ) {
            Err(Error::GwasAnswersOfAnotherSize {
                column,
                num_values,
                num_with_variance,
            }) => {
                assert_eq!((column, num_values, num_with_variance), ("beta", 1, 2));
            }
            other => panic!("the model answered one effect of two, and that gave {other:?}"),
        }

        match result.add_the_block(
            &dosages,
            Answers {
                beta: &[1.5, -0.25],
                se: &[0.6, 0.4, 0.3],
                p_value: &[0.08, 0.5],
            },
        ) {
            Err(Error::GwasAnswersOfAnotherSize {
                column,
                num_values,
                num_with_variance,
            }) => {
                assert_eq!((column, num_values, num_with_variance), ("se", 3, 2));
            }
            other => {
                panic!("the model answered three standard errors of two, and that gave {other:?}")
            }
        }

        assert_eq!(result.num_vars, 0, "no row was added by either");
        assert!(result.allele_freq.is_empty());
    }
}
