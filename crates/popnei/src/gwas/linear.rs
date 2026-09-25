//! The linear model: a continuous trait and no kinship.
//!
//! [`LinearModel`] is the thin QR of the design, the residuals of the trait
//! and the t test of every variant against them. "The linear model" of
//! `docs/specs/gwas.md` says what it fits and how it is verified.

use popnei_linalg::{TheFirstOperand, TheHalfThatHoldsTheMatrix, TheSecondOperand, ThinQr};

use crate::error::{Error, Result};

use super::distributions::t_sf_two_sided;
use super::dosages::GwasDosages;
use super::result::{Answers, NullModel};
use super::study::{Design, GwasInputShape, GwasModel, TestType};
use super::{Answer, the_share_that_is_nothing};

/// The linear model of a study fitted without any variant in it, with the
/// buffers one block of variants is tested in.
///
/// The trait is a straight line in the columns of the design, and the fit
/// is the thin QR of that design, `d = q r`: the coefficients are the `c`
/// of `r c = q' y`, the residuals are `y - q q' y`, what the model left
/// unexplained, and the residual sum of squares is their squared length.
/// That is "The linear model" of `docs/specs/gwas.md`, and it is what
/// plink2's `--glm` fits.
///
/// The `q` and the residuals are kept because every variant is tested
/// against them: the covariates are taken out of the dosages of a block in
/// one product, and the effect of a variant is then the plain slope of the
/// trait's residuals on the variant's. The buffers of a block are kept
/// from one block to the next, so a pass over a million variants asks the
/// machine for them once and allocates nothing for a variant.
pub(crate) struct LinearModel {
    /// The `q` of the thin QR of the design, `num_individuals` x
    /// `num_coefs`, row after row: its columns are of length 1 and at
    /// right angles to each other, and they span what the design spans.
    q: Vec<f64>,
    /// The effect of the intercept and of each covariate, one per column
    /// of the design.
    coefs: Vec<f64>,
    /// What the design left of the trait, one value per tested individual,
    /// which every variant is tested against.
    residuals: Vec<f64>,
    /// The squared length of those residuals, which is what the null model
    /// left unexplained.
    rss: f64,
    /// How many individuals the study tests.
    num_individuals: usize,
    /// How many columns the design has.
    num_coefs: usize,
    /// The individuals less the columns of the design, which the residual
    /// sum of squares is divided by for the variance of the null model. It
    /// is 2 at least, since [`Design::of_the_study`] refuses a study of no
    /// more individuals than its columns plus one.
    degrees_of_freedom_of_the_null: usize,
    /// The same less one more for the variant, which is the degrees of
    /// freedom of the t test of every variant. It is 1 at least.
    degrees_of_freedom: usize,
    /// The dosages of a block times the `q`, the variants that have
    /// variance x `num_coefs`.
    of_the_design: Vec<f64>,
    /// The dosages of a block with the design taken out of them, the
    /// variants that have variance x `num_individuals`.
    residualized: Vec<f64>,
    /// Each variant's residuals times the trait's, one per variant that
    /// has variance.
    num: Vec<f64>,
    /// The effect of each variant that has variance.
    beta: Vec<f64>,
    /// The standard error of each of those effects.
    se: Vec<f64>,
    /// The p-value of each of those tests.
    p_value: Vec<f64>,
}

impl LinearModel {
    /// The linear model of `phenotype` over `design`, fitted without any
    /// variant in it.
    ///
    /// `phenotype` holds one value for each individual the design has a
    /// row for, which is what the study tests.
    ///
    /// # Errors
    ///
    /// [`Error::GwasInputOfAnotherSize`] when `phenotype` does not hold
    /// one value for each tested individual, and [`Error::GwasLinalg`]
    /// when the thin QR of the design, the solve against its `r` or one of
    /// the two products could not be done.
    pub(crate) fn of_the_study(phenotype: &[f64], design: &Design<'_>) -> Result<LinearModel> {
        let num_individuals = design.num_individuals();
        let num_coefs = design.num_coefs();
        if phenotype.len() != num_individuals {
            return Err(Error::GwasInputOfAnotherSize {
                problem: GwasInputShape::Phenotype {
                    num_values: phenotype.len(),
                    num_individuals,
                },
            });
        }
        let ThinQr { q, r } = popnei_linalg::thin_qr(design.values(), num_individuals, num_coefs)
            .map_err(|source| Error::GwasLinalg {
            operation: "thin QR of the design",
            source,
        })?;
        // `q' y`, which both the coefficients and the residuals are built
        // from: the residuals are `y - q (q' y)` and the coefficients the
        // solution of `r c = q' y`, so it is computed once and the solve
        // that overwrites it comes after the residuals are taken.
        let mut coefs = vec![0.0_f64; num_coefs];
        popnei_linalg::product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: &q,
                rows: num_coefs,
            },
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: phenotype,
                cols: 1,
            },
            &mut coefs,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design's q with the trait",
            source,
        })?;
        let mut residuals = vec![0.0_f64; num_individuals];
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &q,
                rows: num_individuals,
            },
            num_coefs,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &coefs,
                cols: 1,
            },
            &mut residuals,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the design's q with the fitted trait",
            source,
        })?;
        for (residual, measured) in residuals.iter_mut().zip(phenotype) {
            *residual = measured - *residual;
        }
        let rss = residuals
            .iter()
            .map(|residual| residual * residual)
            .sum::<f64>();
        popnei_linalg::solve_triangular(
            &r,
            num_coefs,
            TheHalfThatHoldsTheMatrix::TheUpperHalf,
            &mut coefs,
            1,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "solve of the design's r against the trait",
            source,
        })?;
        // The two are written with the plain operator and not with
        // `saturating_sub`: saturating is not the meaning wanted here, and
        // a caller that reached this with fewer individuals would get 0
        // degrees of freedom and an infinite standard error where it
        // should get the error that `Design::of_the_study` raises.
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "`Design::of_the_study` refuses a study of no more individuals than \
                      the columns of its design plus one, and a `Design` is the only way \
                      to reach this, so `num_individuals` is `num_coefs` plus 2 at least"
        )]
        let degrees_of_freedom_of_the_null = num_individuals - num_coefs;
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "the line above is 2 at least, by the same refusal"
        )]
        let degrees_of_freedom = degrees_of_freedom_of_the_null - 1;
        Ok(LinearModel {
            q,
            coefs,
            residuals,
            rss,
            num_individuals,
            num_coefs,
            degrees_of_freedom_of_the_null,
            degrees_of_freedom,
            of_the_design: Vec::new(),
            residualized: Vec::new(),
            num: Vec::new(),
            beta: Vec::new(),
            se: Vec::new(),
            p_value: Vec::new(),
        })
    }

    /// What the fit left unexplained: the squared length of the residuals
    /// of the trait, which each variant's own sum of squares is taken from.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the result of a study carries the residual variance, which is this \
                      over the degrees of freedom of the null, and not this itself; the \
                      worked example of `docs/specs/gwas.md` gives both and the test of \
                      it is what reads this"
        )
    )]
    pub(crate) fn residual_sum_of_squares(&self) -> f64 {
        self.rss
    }

    /// The null model of the result: the effects of the intercept and of
    /// the covariates, and what the fit left unexplained over the
    /// individuals less the columns of the design.
    ///
    /// A linear model has no kinship, so it has neither a genetic variance
    /// nor a heritability.
    #[must_use]
    pub(crate) fn null_model(&self, test: TestType) -> NullModel {
        NullModel {
            model: GwasModel::Lm,
            test,
            covariate_effects: self.coefs.clone(),
            residual_variance: Some(self.rss / self.degrees_of_freedom_of_the_null as f64),
            genetic_variance: None,
            heritability: None,
            num_individuals: self.num_individuals,
        }
    }

    /// The t test of every variant of a block that has variance among the
    /// tested individuals, in the order of the block.
    ///
    /// The covariates are taken out of the dosages of the whole block in
    /// one product, `x - (x q) q'`, and after that the effect of a variant
    /// is the plain slope of the trait's residuals on the variant's: `beta`
    /// is `num / xx`, with `num` the variant's residuals times the trait's
    /// and `xx` the squared length of the variant's. What the variant
    /// leaves unexplained is the squared length of the trait's residuals
    /// less `beta` times the variant's, so each variant gets its own
    /// estimate of the residual variance, which is what makes this a t
    /// test and not a normal one, and `se` is the square root of that over
    /// the degrees of freedom and over `xx`. The p-value is
    /// [`t_sf_two_sided`] of `beta / se`.
    ///
    /// A variant of which the design leaves at most the tested individuals
    /// times 2.2e-16 of its squared length has no answer, and gets the
    /// three NaNs a variant with no variance gets: it is a combination of
    /// the columns of the design, so what is left of it is rounding.
    ///
    /// # Errors
    ///
    /// [`Error::GwasVariantsTooLarge`] when the values of the block are
    /// more than a `usize` counts, and [`Error::GwasLinalg`] when one of
    /// the three products could not be done, which is where a block of
    /// other individuals than the null model was fitted over is refused.
    pub(crate) fn test_the_block(&mut self, dosages: &GwasDosages) -> Result<Answers<'_>> {
        let num_vars = dosages.num_with_variance();
        self.beta.clear();
        self.se.clear();
        self.p_value.clear();
        if num_vars == 0 {
            // No variant to test, and the products below take a matrix of
            // one column at least. The block still has its rows in the
            // result, with the three NaNs of a variant that has no answer.
            return Ok(self.answers());
        }
        let values = num_vars
            .checked_mul(self.num_individuals)
            .ok_or(Error::GwasVariantsTooLarge)?;
        let of_the_design = num_vars
            .checked_mul(self.num_coefs)
            .ok_or(Error::GwasVariantsTooLarge)?;
        self.of_the_design.resize(of_the_design, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: dosages.dosages(),
                rows: num_vars,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.q,
                cols: self.num_coefs,
            },
            &mut self.of_the_design,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of a block of variants with the design's q",
            source,
        })?;
        self.residualized.resize(values, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &self.of_the_design,
                rows: num_vars,
            },
            self.num_coefs,
            TheSecondOperand::ByTheColumnsOfTheResult {
                values: &self.q,
                cols: self.num_individuals,
            },
            &mut self.residualized,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of a block of variants with the design",
            source,
        })?;
        // What the product wrote is the part of each variant the design
        // explains, and what is tested is what is left of the variant.
        for (residualized, dosage) in self.residualized.iter_mut().zip(dosages.dosages()) {
            *residualized = dosage - *residualized;
        }
        self.num.resize(num_vars, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &self.residualized,
                rows: num_vars,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.residuals,
                cols: 1,
            },
            &mut self.num,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of a block of variants with the trait's residuals",
            source,
        })?;
        let rows = TheRowsToTest {
            residualized: &self.residualized,
            num: &self.num,
            sum_of_squares: dosages.sum_of_squares(),
            residuals: &self.residuals,
            num_individuals: self.num_individuals,
            degrees_of_freedom: self.degrees_of_freedom as f64,
            // The share of its own squared length that a variant has to
            // keep once the design is taken out of it to be worth testing.
            share_that_is_nothing: the_share_that_is_nothing(self.num_individuals),
        };
        for (beta, se, p_value) in the_answers_of_the_rows(&rows) {
            self.beta.push(beta);
            self.se.push(se);
            self.p_value.push(p_value);
        }
        Ok(self.answers())
    }

    /// The three columns of what the block last tested answered.
    #[must_use]
    fn answers(&self) -> Answers<'_> {
        Answers {
            beta: &self.beta,
            se: &self.se,
            p_value: &self.p_value,
        }
    }
}

/// The rows of one block as the t test of the linear model reads them,
/// once the three products have been made.
///
/// It is one value and four slices because they travel together through
/// the two ways of walking the rows, and because a caller that swapped two
/// of the slices, all of them `f64` of the same block, would compile.
struct TheRowsToTest<'a> {
    /// What is left of each variant that has variance once the design is
    /// taken out of it, one row of `num_individuals` values for each, in
    /// the order of the block.
    residualized: &'a [f64],
    /// Each of those rows against the trait's residuals, one value for
    /// each of them.
    num: &'a [f64],
    /// The squared length of the dosages of each of those variants, before
    /// the design was taken out, which the block summed where it wrote the
    /// row: one value for each of them.
    sum_of_squares: &'a [f64],
    /// The trait less what the null model explains, one value for each
    /// tested individual, which every variant is tested against.
    residuals: &'a [f64],
    /// How many individuals the study tests, which is the length of one
    /// row.
    num_individuals: usize,
    /// The degrees of freedom of the t test, which is the tested
    /// individuals less the columns of the design and the variant.
    degrees_of_freedom: f64,
    /// The share of its own squared length that a variant has to keep once
    /// the design is taken out of it to be worth testing, which
    /// [`the_share_that_is_nothing`] gives.
    share_that_is_nothing: f64,
}

impl TheRowsToTest<'_> {
    /// The effect of one variant, its standard error and its p-value, and
    /// the three NaNs that "The variants that have no answer" of
    /// `docs/specs/gwas.md` gives a variant that has none.
    ///
    /// `row` is what the design leaves of the variant, `num` is that row
    /// against the trait's residuals, and `of_the_dosages` is the squared
    /// length the variant had before the design was taken out of it.
    /// Nothing outside the three of them and the study's own values is
    /// read, and nothing is written, which is what lets the rows be walked
    /// on the threads of rayon.
    fn the_answer_of_the_row(&self, row: &[f64], num: f64, of_the_dosages: f64) -> Answer {
        let xx = row.iter().map(|value| value * value).sum::<f64>();
        // A variant that is a combination of the columns of the design has
        // nothing left once they are taken out, and what `xx` holds is the
        // rounding of that cancellation: `beta` would be a number divided
        // by noise, large and of whichever sign the rounding chose, and
        // the two backends do not choose the same one. Such a variant has
        // no answer, as one with no variance has.
        if xx <= self.share_that_is_nothing * of_the_dosages {
            return (f64::NAN, f64::NAN, f64::NAN);
        }
        let beta = num / xx;
        // What the variant leaves unexplained, formed from its own
        // residuals and not taken from the null model's sum of squares by
        // subtracting `beta * num`: those two quantities agree to their
        // last bits once a variant explains most of what the null left,
        // and the subtraction then gives the rounding of a cancelled sum,
        // which is 0 or negative as often as not. "The linear model" of
        // `docs/specs/gwas.md` measures what that gave: an `se` of 0 at
        // one variant and NaN at another, and the two backends disagreeing
        // about which.
        let rss = row
            .iter()
            .zip(self.residuals)
            .map(|(value, residual)| {
                let left = residual - beta * value;
                left * left
            })
            .sum::<f64>();
        let se = (rss / self.degrees_of_freedom / xx).sqrt();
        (beta, se, t_sf_two_sided(beta / se, self.degrees_of_freedom))
    }
}

/// The answer of every row of a block, made on the threads of rayon, in
/// the order of the block.
///
/// The rows are tested on those threads as section 3 of
/// `docs/architecture.md` asks, and for the same reason the dosages of a
/// block are read on them: no row reads another. Every sum of
/// [`TheRowsToTest::the_answer_of_the_row`] runs over the values of one
/// variant alone, left to right, and nothing is written, so the three
/// columns are the same bits whatever the threads do and however many of
/// them there are. They are read back in the order of the block, and a
/// variant with no answer carries its three NaNs as a value instead of
/// being skipped, which is what a row written by index cannot do.
///
/// The threads are those of the pool the caller is running in, and rayon's
/// global pool only when the caller is in none. No product is made here:
/// the three of the block were made before it, outside rayon, as section 3
/// asks of a large product.
#[cfg(not(target_family = "wasm"))]
fn the_answers_of_the_rows(rows: &TheRowsToTest<'_>) -> Vec<Answer> {
    use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
    use rayon::slice::ParallelSlice;

    rows.residualized
        .par_chunks_exact(rows.num_individuals)
        .zip(rows.num.par_iter().copied())
        .zip(rows.sum_of_squares.par_iter().copied())
        .map(|((row, num), of_the_dosages)| rows.the_answer_of_the_row(row, num, of_the_dosages))
        .collect()
}

/// The same rows, tested one after another, which is what WebAssembly
/// does: it has no threads.
#[cfg(target_family = "wasm")]
fn the_answers_of_the_rows(rows: &TheRowsToTest<'_>) -> Vec<Answer> {
    the_answers_of_the_rows_one_by_one(rows)
}

/// The rows tested one after another: what WebAssembly does, and what the
/// test that compares the two ways of testing a block calls.
#[cfg(any(target_family = "wasm", test))]
fn the_answers_of_the_rows_one_by_one(rows: &TheRowsToTest<'_>) -> Vec<Answer> {
    rows.residualized
        .chunks_exact(rows.num_individuals)
        .zip(rows.num.iter().copied())
        .zip(rows.sum_of_squares.iter().copied())
        .map(|((row, num), of_the_dosages)| rows.the_answer_of_the_row(row, num, of_the_dosages))
        .collect()
}

/// The linear model against the worked example of `docs/specs/gwas.md`,
/// whose numbers are pyNei's, and against plink2 on the panel with every
/// genotype called.
#[cfg(test)]
pub(crate) mod lm {
    use std::collections::HashMap;
    use std::io::Cursor;
    use std::path::{Path, PathBuf};

    use super::LinearModel;
    use crate::block::BlockReader;
    use crate::error::{Error, Result};
    use crate::gwas::calc_gwas;
    use crate::gwas::result::Gwas;
    use crate::gwas::study::{Design, GwasInput, GwasInputShape, GwasModel, TestType, TraitType};
    use crate::io::vcf::{VcfOptions, VcfReader};

    /// How far a number of the worked example may be from pyNei's, as a
    /// share of that number: 3e-15.
    ///
    /// "The worked example" of `docs/specs/gwas.md` asks for 1e-12, and
    /// "How it is verified" of "What every model shares" asks that such a
    /// number be lowered until it fails and then set two or three times
    /// above where it broke. It was, on 23 September 2026 on both
    /// backends: it breaks at 1e-15, where the `beta` of `v1` on faer is
    /// 1.24e-15 away, and the worst difference of the thirteen numbers it
    /// holds, the four of the null model and the nine of the three rows,
    /// is that 1.24e-15 on faer and 7.4e-16 on Accelerate, the `se` of
    /// `v0`. So this is 2.4 times the worst measured, and the spec's 1e-12
    /// had about 800 times it.
    ///
    /// No number of the example is near 0: the smallest is the `beta` of
    /// `v1`, 0.3125, against an `se` of 1.3, so a share of the number and
    /// a share of the scale of what is estimated are the same bound here,
    /// and the rule of the spec that a `beta` is measured against its `se`
    /// changes nothing. The example runs on whole dosages and a trait of
    /// whole numbers, and it reads no reference file, so nothing of it is
    /// rounded away before the comparison.
    const OF_THE_WORKED_EXAMPLE: f64 = 3e-15;

    /// How far a `beta` or an `se` of the panel may be from plink2's, as a
    /// share of the `se` plink2 printed for that variant.
    ///
    /// It is the 1e-5 of "How it is verified" of "The linear model" of
    /// `docs/specs/gwas.md`, which is a share of the `se` and not of the
    /// value, because a study is mostly null and a `beta` that cancelled
    /// to near 0 is no guide to its own error. plink2 prints six
    /// significant digits, so the number it is compared against is itself
    /// rounded by up to 5e-7 absolute, which of the 1.2e-6 this bound
    /// comes to at the smallest `se` of the six is 41 per cent. The
    /// arithmetic gets what is left, and what it took is in the doc
    /// comment of the test.
    ///
    /// This one is not lowered, unlike [`OF_THE_WORKED_EXAMPLE`]: it is a
    /// bound on what plink2 printed and not on what popnei computed, and
    /// tightening it would measure the rounding of six digits.
    const OF_PLINK2: f64 = 1e-5;

    /// How far a p-value of the panel may be from plink2's, as a share of
    /// it. plink2 prints it to six significant digits like the rest, so
    /// the printing alone can take half of this.
    const OF_PLINK2_P_VALUE: f64 = 1e-5;

    /// How far an `allele_freq` of the panel may be from plink2's
    /// `A1_FREQ`: 1e-6 absolute, since it is a frequency and lies between
    /// 0 and 1.
    const OF_PLINK2_FREQUENCY: f64 = 1e-6;

    /// The path of one of the files of `tests/reference/gwas/`, where the
    /// trait, the covariates and plink2's answers are.
    pub(crate) fn the_reference_path(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/gwas")
            .join(name)
    }

    /// The panel with every genotype called, 200 individuals and 1200
    /// biallelic diploid variants, which is the panel plink2 was run on.
    pub(crate) fn the_panel_path() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/kinship/panel_called.vcf.gz")
    }

    /// The study of the variants of `reader`, with no second pass for the
    /// GRAMMAR-Gamma approximation, which a linear model does not make.
    pub(crate) fn the_study_of<R: BlockReader>(
        reader: &mut R,
        input: &GwasInput<'_>,
    ) -> Result<Gwas> {
        calc_gwas(reader, None::<&mut R>, input)
    }

    /// The worked example of `docs/specs/gwas.md` as a VCF: six diploid
    /// individuals, `i0` to `i5`, and three variants, of which the third
    /// has every individual heterozygous and so no variance. The `./.` of
    /// `i3` at `v1` is the genotype that takes the mean dosage of its
    /// variant, 0.8.
    fn the_worked_example_vcf() -> Vec<u8> {
        let mut vcf = String::from(THE_HEADER);
        for (var, genotypes) in [
            ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1"],
            ["0/0", "0/1", "1/1", "./.", "0/1", "0/0"],
            ["0/1", "0/1", "0/1", "0/1", "0/1", "0/1"],
        ]
        .iter()
        .enumerate()
        {
            let pos = var.saturating_add(1).saturating_mul(1000);
            vcf.push_str(&format!("1\t{pos}\tv{var}\tA\tT\t.\t.\t.\tGT"));
            for genotype in genotypes {
                vcf.push('\t');
                vcf.push_str(genotype);
            }
            vcf.push('\n');
        }
        vcf.into_bytes()
    }

    /// The header of that VCF, six diploid individuals and no variant.
    const THE_HEADER: &str = "##fileformat=VCFv4.2\n\
        ##contig=<ID=1>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti0\ti1\ti2\ti3\ti4\ti5\n";

    /// The trait of the worked example, one value for each of the six
    /// individuals.
    const THE_TRAIT: [f64; 6] = [2.0, 3.0, 5.0, 4.0, 4.0, 7.0];

    /// Its design: the intercept and the one covariate, row after row.
    const THE_DESIGN: [f64; 12] = [
        1.0, 0.0, //
        1.0, 1.0, //
        1.0, 0.0, //
        1.0, 1.0, //
        1.0, 0.0, //
        1.0, 1.0,
    ];

    /// The six individuals of the worked example, which are every
    /// individual of its VCF.
    const THE_INDIVIDUALS: [usize; 6] = [0, 1, 2, 3, 4, 5];

    /// A reader over the bytes of a VCF of diploid individuals.
    pub(crate) fn reader_over(vcf: &[u8]) -> VcfReader<Cursor<Vec<u8>>> {
        match VcfReader::new(Cursor::new(vcf.to_vec()), VcfOptions::default()) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        }
    }

    /// The study of the worked example: a continuous trait, one covariate,
    /// no kinship and the default test.
    fn the_worked_example_study() -> GwasInput<'static> {
        GwasInput {
            phenotype: &THE_TRAIT,
            trait_type: TraitType::Continuous,
            design: &THE_DESIGN,
            num_coefs: 2,
            kinship: None,
            test: None,
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS,
            transform_to_biallelic: false,
        }
    }

    /// `found` within `tolerance` of `expected` as a share of `expected`.
    fn assert_within(found: f64, expected: f64, tolerance: f64, what: &str) {
        let difference = (found - expected).abs();
        assert!(
            difference <= tolerance * expected.abs(),
            "{what} is {found} and the literal is {expected}, {difference} away, which is \
             {share} of it against the {tolerance} allowed",
            share = difference / expected.abs()
        );
    }

    /// `found` within `tolerance` times `scale`, which is the `se` of the
    /// variant for a `beta` and for an `se`.
    fn assert_within_the_scale(found: f64, expected: f64, scale: f64, tolerance: f64, what: &str) {
        let difference = (found - expected).abs();
        assert!(
            difference <= tolerance * scale,
            "{what} is {found} and plink2 gives {expected}, {difference} away, which is \
             {share} of the {scale} it is uncertain by against the {tolerance} allowed",
            share = difference / scale
        );
    }

    /// The null model of the worked example is the one pyNei gives at
    /// commit ef0ca6e: an intercept of 3.6666666666666683 and a covariate
    /// effect of 1.0, a residual sum of squares of 13.333333333333336 over
    /// its 4 degrees of freedom, and so a residual variance of
    /// 3.333333333333334.
    ///
    /// It is fitted here and not through `calc_gwas`, because the result
    /// of a study carries the residual variance and not the sum of squares
    /// the spec's table gives beside it.
    #[test]
    fn the_null_of_the_worked_example_is_pyneis_coefficients_and_sum_of_squares() {
        let study = the_worked_example_study();
        let design = match Design::of_the_study(&study, THE_INDIVIDUALS.len()) {
            Ok(design) => design,
            Err(error) => panic!("the design of the worked example: {error}"),
        };
        let fitted = match LinearModel::of_the_study(&THE_TRAIT, &design) {
            Ok(fitted) => fitted,
            Err(error) => panic!("the null model of the worked example: {error}"),
        };
        let null = fitted.null_model(TestType::Wald);
        assert_eq!(null.model, GwasModel::Lm);
        assert_eq!(null.test, TestType::Wald);
        assert_eq!(null.num_individuals, 6);
        assert_eq!(null.covariate_effects.len(), 2);
        assert_within(
            null.covariate_effects[0],
            3.666_666_666_666_668_3,
            OF_THE_WORKED_EXAMPLE,
            "the intercept",
        );
        assert_within(
            null.covariate_effects[1],
            1.0,
            OF_THE_WORKED_EXAMPLE,
            "the effect of the covariate",
        );
        assert_within(
            fitted.residual_sum_of_squares(),
            13.333_333_333_333_336,
            OF_THE_WORKED_EXAMPLE,
            "the residual sum of squares",
        );
        let residual_variance = null
            .residual_variance
            .expect("the residual variance of a continuous trait");
        assert_within(
            residual_variance,
            3.333_333_333_333_334,
            OF_THE_WORKED_EXAMPLE,
            "the residual variance",
        );
        assert_eq!(null.genetic_variance, None, "a model with no kinship");
        assert_eq!(null.heritability, None, "a model with no kinship");
    }

    /// The three rows of the worked example are pyNei's: `v0` and `v1` get
    /// the frequency, the effect, the standard error and the p-value of
    /// the spec's table, and `v2`, where every individual is
    /// heterozygous, keeps its frequency of 0.5 and has three NaNs.
    ///
    /// This is the one test of the module that needs no reference program
    /// and no reference file. The three columns of the variants come out
    /// of it as well, since the example is read from a VCF that carries
    /// them.
    #[test]
    fn the_worked_example_gives_pyneis_three_rows_and_no_answer_for_the_third() {
        let vcf = the_worked_example_vcf();
        let mut reader = reader_over(&vcf);
        let study = the_worked_example_study();
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of the worked example: {error}"),
        };
        assert_eq!(result.num_vars, 3);
        assert_eq!(
            result.ids.as_deref(),
            Some(["v0".to_owned(), "v1".to_owned(), "v2".to_owned()].as_slice()),
            "the ids of the variants"
        );
        assert_eq!(
            result.poss.as_deref(),
            Some([1000_u64, 2000, 3000].as_slice()),
            "the positions of the variants"
        );
        assert_eq!(
            result.chroms.as_deref(),
            Some([0_u32, 0, 0].as_slice()),
            "the chromosomes of the variants"
        );
        assert_eq!(result.chrom_table.name(0), Some("1"));
        assert!(!result.used_grammar_gamma_approx);
        assert_eq!(result.null_model.model, GwasModel::Lm);
        for (var, (allele_freq, beta, se, p_value)) in [
            (0.5, 1.5, 0.600_925_212_577_332, 0.088_004_892_382_756),
            (0.4, 0.3125, 1.305_204_592_306_424, 0.826_200_867_452_417),
        ]
        .into_iter()
        .enumerate()
        {
            assert_within(
                result.allele_freq[var],
                allele_freq,
                OF_THE_WORKED_EXAMPLE,
                &format!("the frequency of v{var}"),
            );
            assert_within(
                result.beta[var],
                beta,
                OF_THE_WORKED_EXAMPLE,
                &format!("the effect of v{var}"),
            );
            assert_within(
                result.se[var],
                se,
                OF_THE_WORKED_EXAMPLE,
                &format!("the standard error of v{var}"),
            );
            assert_within(
                result.p_value[var],
                p_value,
                OF_THE_WORKED_EXAMPLE,
                &format!("the p-value of v{var}"),
            );
        }
        assert_within(
            result.allele_freq[2],
            0.5,
            OF_THE_WORKED_EXAMPLE,
            "the frequency of v2",
        );
        assert!(result.beta[2].is_nan(), "the effect of v2");
        assert!(result.se[2].is_nan(), "the standard error of v2");
        assert!(result.p_value[2].is_nan(), "the p-value of v2");
    }

    /// How many variants a study of six individuals reads in one block,
    /// which is [`crate::block::MAX_NUM_VARS_PER_BLOCK`]: the pass puts a
    /// `Reblock` before it, so the blocks it sees are of that size
    /// whatever size the reader gives, and a source has to pass it for the
    /// pass to read a second block.
    const VARS_OF_ONE_BLOCK: usize = 10_000;

    /// A study over more variants than one block holds gives the same
    /// answers in its second block as in its first, and one row for every
    /// variant in the order the reader gave them.
    ///
    /// The variants are the three of the worked example over and over, so
    /// every variant of the study has the answer the spec's table gives
    /// for its pattern, whichever block it fell in: the one with no
    /// variance keeps its frequency and has three NaNs wherever it is.
    /// What this covers that the worked example does not is the buffers of
    /// the model and of the dosages being reused, the rows of the second
    /// block being added after the first's and not over them, and the
    /// three columns of the variants growing across the blocks.
    ///
    /// The second block is on another chromosome than the first, so that
    /// the column of chromosomes is asserted where it is built: the
    /// numbers a block holds are read through the table of the reader, and
    /// a study that kept the first block's name for every variant after it
    /// would answer every question this test asks about the ids and the
    /// positions and still be wrong about where a variant is.
    #[test]
    fn a_study_of_more_variants_than_one_block_answers_the_same_in_every_block() {
        let patterns = [
            ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1"],
            ["0/0", "0/1", "1/1", "./.", "0/1", "0/0"],
            ["0/1", "0/1", "0/1", "0/1", "0/1", "0/1"],
        ];
        let num_vars = VARS_OF_ONE_BLOCK.saturating_add(100);
        let mut vcf = String::from(THE_HEADER);
        for var in 0..num_vars {
            let pos = var.saturating_add(1).saturating_mul(10);
            // The variants of the second block are on the second
            // chromosome, and the first block fills a whole block.
            let chrom = match var < VARS_OF_ONE_BLOCK {
                true => 1,
                false => 2,
            };
            vcf.push_str(&format!("{chrom}\t{pos}\tv{var}\tA\tT\t.\t.\t.\tGT"));
            for genotype in patterns[var % patterns.len()] {
                vcf.push('\t');
                vcf.push_str(genotype);
            }
            vcf.push('\n');
        }
        let mut reader = reader_over(vcf.as_bytes());
        let study = the_worked_example_study();
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of {num_vars} variants: {error}"),
        };
        assert_eq!(result.num_vars, num_vars, "the variants of the study");
        assert!(
            num_vars > VARS_OF_ONE_BLOCK,
            "the study has to read more than one block"
        );
        let ids = result.ids.as_deref().expect("the ids of the variants");
        let poss = result.poss.as_deref().expect("the positions");
        let chroms = result.chroms.as_deref().expect("the chromosomes");
        assert_eq!(ids.len(), num_vars, "one id for each variant");
        assert_eq!(poss.len(), num_vars, "one position for each variant");
        assert_eq!(chroms.len(), num_vars, "one chromosome for each variant");
        for var in [0, VARS_OF_ONE_BLOCK, num_vars.saturating_sub(1)] {
            assert_eq!(ids[var], format!("v{var}"), "the id of the variant {var}");
            assert_eq!(
                poss[var],
                var.saturating_add(1).saturating_mul(10) as u64,
                "the position of the variant {var}"
            );
            let name = result.chrom_table.name(chroms[var]);
            let expected = match var < VARS_OF_ONE_BLOCK {
                true => "1",
                false => "2",
            };
            assert_eq!(
                name,
                Some(expected),
                "the chromosome of the variant {var}, which the table of the pass names"
            );
        }
        for (var, (allele_freq, beta, se, p_value)) in (0..num_vars)
            .map(|var| {
                (
                    var,
                    match var % patterns.len() {
                        0 => (0.5, 1.5, 0.600_925_212_577_332, 0.088_004_892_382_756),
                        1 => (0.4, 0.3125, 1.305_204_592_306_424, 0.826_200_867_452_417),
                        _ => (0.5, f64::NAN, f64::NAN, f64::NAN),
                    },
                )
            })
            .collect::<Vec<_>>()
        {
            assert_within(
                result.allele_freq[var],
                allele_freq,
                OF_THE_WORKED_EXAMPLE,
                &format!("the frequency of the variant {var}"),
            );
            if beta.is_nan() {
                assert!(result.beta[var].is_nan(), "the effect of the variant {var}");
                assert!(result.se[var].is_nan(), "the error of the variant {var}");
                assert!(
                    result.p_value[var].is_nan(),
                    "the p-value of the variant {var}"
                );
                continue;
            }
            assert_within(
                result.beta[var],
                beta,
                OF_THE_WORKED_EXAMPLE,
                &format!("the effect of the variant {var}"),
            );
            assert_within(
                result.se[var],
                se,
                OF_THE_WORKED_EXAMPLE,
                &format!("the error of the variant {var}"),
            );
            assert_within(
                result.p_value[var],
                p_value,
                OF_THE_WORKED_EXAMPLE,
                &format!("the p-value of the variant {var}"),
            );
        }
    }

    /// The trait of `trait_type` and the two covariates of the panel, read
    /// from `tests/reference/gwas/phenotypes.csv`, in the order
    /// `individuals` has the individuals, which is the order the VCF has
    /// them.
    ///
    /// The file holds one line per individual with its name, the
    /// continuous trait, the binomial one, the two covariates and the
    /// subpopulation. What is taken here is one of the two traits and both
    /// covariates, which is what every reference program was given: `cont`
    /// for a continuous trait and `binom`, 0 or 1, for a binomial one.
    pub(crate) fn the_trait_and_the_design_of_the_panel(
        individuals: &[String],
        trait_type: TraitType,
    ) -> (Vec<f64>, Vec<f64>) {
        let path = the_reference_path("phenotypes.csv");
        let text = match std::fs::read_to_string(&path) {
            Ok(text) => text,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let mut lines = text.lines();
        let header: Vec<&str> = lines
            .next()
            .expect("the header of phenotypes.csv")
            .split(',')
            .collect();
        let column_of = |name: &str| match header.iter().position(|held| *held == name) {
            Some(at) => at,
            None => panic!("phenotypes.csv has no column {name}"),
        };
        let of_the_name = column_of("IID");
        let of_the_trait = column_of(match trait_type {
            TraitType::Continuous => "cont",
            TraitType::Binomial => "binom",
        });
        let of_cov1 = column_of("cov1");
        let of_cov2 = column_of("cov2");
        let mut rows: HashMap<String, [f64; 3]> = HashMap::new();
        for line in lines.filter(|line| !line.is_empty()) {
            let fields: Vec<&str> = line.split(',').collect();
            let number = |at: usize| match fields.get(at).map(|field| field.parse::<f64>()) {
                Some(Ok(number)) => number,
                _ => panic!("the column {at} of the line `{line}` of phenotypes.csv"),
            };
            let name = match fields.get(of_the_name) {
                Some(name) => (*name).to_owned(),
                None => panic!("the line `{line}` of phenotypes.csv has no name"),
            };
            rows.insert(
                name,
                [number(of_the_trait), number(of_cov1), number(of_cov2)],
            );
        }
        let mut phenotype = Vec::with_capacity(individuals.len());
        let mut design = Vec::with_capacity(individuals.len().saturating_mul(3));
        for individual in individuals {
            match rows.get(individual) {
                Some(row) => {
                    phenotype.push(row[0]);
                    design.extend_from_slice(&[1.0, row[1], row[2]]);
                }
                None => panic!("{individual} is not in phenotypes.csv"),
            }
        }
        (phenotype, design)
    }

    /// What plink2 v2.0.0-a.7.7 wrote for six variants of the panel with
    /// every genotype called, in
    /// `tests/reference/gwas/plink2.panel_called.glm.linear.tsv`, with
    /// `cov1` and `cov2` as covariates: the id and then `A1_FREQ`, `BETA`,
    /// `SE` and `P`, the four columns "How it is verified" of "The linear
    /// model" of `docs/specs/gwas.md` compares. Five of them are the
    /// causal variants of `causal_vars.csv` and `var0000` is not causal.
    const OF_PLINK2_SIX: [(&str, f64, f64, f64, f64); 6] = [
        ("var0000", 0.21, -0.424136, 0.139354, 0.00265846),
        ("var0052", 0.3225, -0.697724, 0.122348, 4.28981e-08),
        ("var0629", 0.175, -0.813852, 0.161809, 1.10646e-06),
        ("var0751", 0.41, -0.0963977, 0.130636, 0.461451),
        ("var1137", 0.2625, -0.171137, 0.145987, 0.242511),
        ("var1188", 0.4675, -0.655393, 0.138912, 4.51958e-06),
    ];

    /// The six variants of the panel are plink2's: `beta` and `se` within
    /// 1e-5 of the `se` plink2 printed for the variant, `p_value` within
    /// 1e-5 of itself and `allele_freq` within 1e-6 absolute, which is
    /// what the spec asks of the whole columns in Python.
    ///
    /// Measured over the six on 23 September 2026, the same on Accelerate
    /// and on faer to two digits: the worst `beta` or `se` is 3.7e-6 of
    /// the `se` of its variant, 37 per cent of the 1e-5 allowed, and it is
    /// the `se` of `var0052`; the worst p-value is 2.8e-6 of itself, 28
    /// per cent of its 1e-5; and every `allele_freq` is plink2's
    /// `A1_FREQ` exactly. plink2 prints six significant digits, so its own
    /// rounding of a `beta` or an `se` below 1 is up to 5e-7 absolute,
    /// which at the smallest `se` of the six, the 0.122348 of `var0052`,
    /// is 4.1e-6 of that `se`: the whole of what was measured fits inside
    /// what the printing alone can account for, and the arithmetic has
    /// spent nothing that can be seen here.
    ///
    /// Over all 1200 variants the same run gives 1.94e-5 of the `se`,
    /// which this bound would not hold, and the reason is the printing and
    /// not the arithmetic: `var0482` has a `beta` of 1.0389, and six
    /// significant digits of a value above 1 are rounded by up to 5e-6
    /// absolute, which is 2.6e-5 of its `se` of 0.19. The two variants of
    /// the panel that pass 1e-5 are the two whose `beta` passes 1. The
    /// check over the whole column is the Python one of work package 3 of
    /// `docs/plans/gwas-linear.md`, and it is where that has to be dealt
    /// with; none of the six variants here has a `beta` above 1.
    #[test]
    fn the_six_variants_of_the_panel_are_plink2s_effect_error_and_p_value() {
        let path = the_panel_path();
        let options = VcfOptions {
            ploidy: 2,
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let individuals = reader.individuals().to_vec();
        let (phenotype, design) =
            the_trait_and_the_design_of_the_panel(&individuals, TraitType::Continuous);
        let tested: Vec<usize> = (0..individuals.len()).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Continuous,
            design: &design,
            num_coefs: 3,
            kinship: None,
            test: None,
            use_grammar_gamma_approx: false,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of the panel: {error}"),
        };
        assert_eq!(result.num_vars, 1200, "the variants of the panel");
        assert_eq!(result.null_model.num_individuals, 200);
        let ids = result.ids.as_deref().expect("the ids of the variants");
        for (id, allele_freq, beta, se, p_value) in OF_PLINK2_SIX {
            let var = match ids.iter().position(|held| held == id) {
                Some(var) => var,
                None => panic!("{id} is not a variant of the panel"),
            };
            let found = result.allele_freq[var];
            assert!(
                (found - allele_freq).abs() <= OF_PLINK2_FREQUENCY,
                "the frequency of {id} is {found} and plink2 gives {allele_freq}"
            );
            assert_within_the_scale(
                result.beta[var],
                beta,
                se,
                OF_PLINK2,
                &format!("the effect of {id}"),
            );
            assert_within_the_scale(
                result.se[var],
                se,
                se,
                OF_PLINK2,
                &format!("the standard error of {id}"),
            );
            assert_within(
                result.p_value[var],
                p_value,
                OF_PLINK2_P_VALUE,
                &format!("the p-value of {id}"),
            );
        }
    }

    /// The GRAMMAR-Gamma approximation stands in for the denominator of a
    /// mixed model's test, so a study with no kinship is refused for
    /// asking for it, which is pyNei's refusal of the same pair.
    #[test]
    fn the_grammar_gamma_approximation_without_a_kinship_is_refused() {
        let vcf = the_worked_example_vcf();
        let mut reader = reader_over(&vcf);
        let study = GwasInput {
            use_grammar_gamma_approx: true,
            ..the_worked_example_study()
        };
        match the_study_of(&mut reader, &study) {
            Err(Error::GwasGrammarGammaWithoutAKinship) => {}
            Err(error) => panic!("the study was refused with {error}"),
            Ok(_) => panic!("the approximation was made without a kinship"),
        }
    }

    /// The trait of the six individuals for a variant that explains all
    /// but `delta` of what the null model left, as literals.
    ///
    /// The trait is `1 + 2 * cov + 3 * dosage + delta * u` over the
    /// covariate `0 1 0 1 0 1` and the dosages `0 1 2 0 1 2` of `v0`, with
    /// `u` the direction `1 0 -1 -1 0 1`, which is at right angles to the
    /// column of ones, to the covariate and to the dosages: `u` sums to 0,
    /// its values at the three individuals of `cov = 1` sum to 0, and its
    /// product with the dosages is `-2 + 2`. So the effect of the variant
    /// is 3 whatever `delta` is, and what the variant leaves unexplained
    /// is `delta² * 4`, `4` being the squared length of `u`.
    ///
    /// The values are written out and not computed, so that the test reads
    /// the same numbers numpy was given: `repr` of a float prints the
    /// shortest text that reads back as the same number, and this is what
    /// numpy 2.5.3 printed for `1 + 2 * cov + 3 * x + delta * u`.
    const THE_TRAIT_A_VARIANT_EXPLAINS: [(f64, [f64; 6], f64); 4] = [
        (
            1e-5,
            [1.00001, 6.0, 6.99999, 2.99999, 4.0, 9.00001],
            5.773_502_691_645_637e-6,
        ),
        (
            1e-7,
            [1.0000001, 6.0, 6.9999999, 2.9999999, 4.0, 9.0000001],
            5.773_502_679_242_529e-8,
        ),
        (
            1e-8,
            [1.00000001, 6.0, 6.99999999, 2.99999999, 4.0, 9.00000001],
            5.773_502_897_178_33e-9,
        ),
        (0.0, [1.0, 6.0, 7.0, 3.0, 4.0, 9.0], 0.0),
    ];

    /// How far the `se` of such a variant may be from the one numpy gives
    /// the same trait: 1e-5 of it.
    ///
    /// It is not tighter because numpy's own number is not. The trait it
    /// was given is the one above, whose values are rounded to the nearest
    /// `f64` and so carry a perturbation of about 1e-16 that is not along
    /// `u`; at a `delta` of 1e-8 that perturbation is 1e-8 of the signal,
    /// so numpy's `se` and the closed form `delta / sqrt(3)` part in their
    /// eighth digit, 5.7735028971e-9 against 5.7735026918e-9. What is
    /// asserted is that popnei answers the rounded trait it was given as
    /// numpy answers it.
    const OF_NUMPYS_STANDARD_ERROR: f64 = 1e-5;

    /// A variant that explains almost all of what the null model left has
    /// a standard error, and the two backends agree on it.
    ///
    /// What a variant leaves unexplained is formed from its residuals and
    /// not by subtracting `beta * num` from the null's sum of squares,
    /// which "The linear model" of `docs/specs/gwas.md` asks for and
    /// measures: the two quantities agree to their last bits once the
    /// variant explains most of the residual, so the subtraction leaves
    /// the rounding of a cancelled sum, of either sign. At a `delta` of
    /// 1e-7 it gave exactly 0, and so an `se` of 0, which is not a
    /// standard error; at 1e-8 it gave -7.1e-15, and so an `se` of NaN,
    /// which is not one either; and which of the two a variant got
    /// differed between Accelerate and faer, so the native build and the
    /// build a browser runs disagreed about whether a variant could be
    /// tested at all.
    ///
    /// The effect is 3 at every `delta` by the construction of the trait,
    /// and the `se` is numpy 2.5.3's, from its own least squares fit of
    /// the trait on the intercept, the covariate and the dosages, with the
    /// standard error taken from the inverse of `d' d`. At a `delta` of 0
    /// the variant explains the trait exactly and there is no number to
    /// assert: what is checked is that the `se` is a finite number at the
    /// size of the rounding and not a NaN, and numpy gives 1.2e-15 there.
    #[test]
    fn a_variant_that_explains_almost_everything_has_a_standard_error() {
        let vcf = the_worked_example_vcf();
        for (delta, phenotype, se_of_numpy) in THE_TRAIT_A_VARIANT_EXPLAINS {
            let study = GwasInput {
                phenotype: &phenotype,
                ..the_worked_example_study()
            };
            let mut reader = reader_over(&vcf);
            let result = match the_study_of(&mut reader, &study) {
                Ok(result) => result,
                Err(error) => panic!("the study at a delta of {delta}: {error}"),
            };
            assert_within(
                result.beta[0],
                3.0,
                1e-12,
                &format!("the effect at a delta of {delta}"),
            );
            let se = result.se[0];
            assert!(
                se.is_finite() && se >= 0.0,
                "the standard error at a delta of {delta} is {se}, and what a variant \
                 leaves unexplained is a squared length, which is 0 or above"
            );
            let p_value = result.p_value[0];
            assert!(
                p_value.is_finite(),
                "the p-value at a delta of {delta} is {p_value}"
            );
            if delta == 0.0 {
                assert!(
                    se < 1e-13,
                    "the variant explains the trait exactly at a delta of 0, and its \
                     standard error is {se}, which is more than the rounding of the sum"
                );
                continue;
            }
            assert_within(
                se,
                se_of_numpy,
                OF_NUMPYS_STANDARD_ERROR,
                &format!("the standard error at a delta of {delta}"),
            );
        }
    }

    /// The header of a VCF of eight diploid individuals and no variant.
    pub(crate) const THE_HEADER_OF_EIGHT: &str = "##fileformat=VCFv4.2\n\
        ##contig=<ID=1>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti0\ti1\ti2\ti3\ti4\ti5\ti6\ti7\n";

    /// The design of the fixture of **Open 2** of `docs/specs/gwas.md`:
    /// eight individuals and a covariate that marks two subpopulations of
    /// four, beside the intercept.
    pub(crate) const THE_DESIGN_OF_TWO_SUBPOPULATIONS: [f64; 16] = [
        1.0, 0.0, //
        1.0, 0.0, //
        1.0, 0.0, //
        1.0, 0.0, //
        1.0, 1.0, //
        1.0, 1.0, //
        1.0, 1.0, //
        1.0, 1.0,
    ];

    /// The trait of those eight individuals.
    pub(crate) const THE_TRAIT_OF_TWO_SUBPOPULATIONS: [f64; 8] =
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];

    /// All eight of them are tested.
    pub(crate) const THE_INDIVIDUALS_OF_TWO_SUBPOPULATIONS: [usize; 8] = [0, 1, 2, 3, 4, 5, 6, 7];

    /// What numpy 2.5.3 answers for the ordinary variant of that fixture,
    /// the dosages `1 0 2 0 1 2 0 1`, by the same least squares fit as
    /// above: the threshold that refuses the variant beside it leaves this
    /// one alone.
    const OF_THE_ORDINARY_VARIANT: (f64, f64) = (-0.315_789_473_684_210_9, 0.633_330_903_431_213_6);

    /// A variant that the design leaves nothing of has no answer, which is
    /// the meanwhile of **Open 2** of `docs/specs/gwas.md`.
    ///
    /// The variant is twice the covariate, so the design explains all of
    /// it and what is left is rounding: `beta` is a number divided by
    /// noise, and what the study gave before this was 5.9e13 on Accelerate
    /// and -3.0e13 on faer, opposite signs with a p-value that reads as a
    /// variant that was tested and showed nothing. plink2 answers `NA` for
    /// such a variant with `ERRCODE CORR_TOO_HIGH`.
    ///
    /// The threshold is the spec's: the variant is refused when what the
    /// design leaves of its squared length is at most the tested
    /// individuals times 2.2e-16 of what there was. Measured with numpy
    /// 2.5.3 on this fixture, the collinear variant leaves 6.47e-32 of its
    /// squared length against a threshold of 1.78e-15, and the ordinary
    /// variant beside it leaves 0.432, thirteen orders of magnitude apart.
    #[test]
    fn a_variant_the_design_leaves_nothing_of_has_no_answer() {
        let mut vcf = String::from(THE_HEADER_OF_EIGHT);
        for (var, genotypes) in [
            ["0/0", "0/0", "0/0", "0/0", "1/1", "1/1", "1/1", "1/1"],
            ["0/1", "0/0", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1"],
        ]
        .iter()
        .enumerate()
        {
            let pos = var.saturating_add(1).saturating_mul(1000);
            vcf.push_str(&format!("1\t{pos}\tv{var}\tA\tT\t.\t.\t.\tGT"));
            for genotype in genotypes {
                vcf.push('\t');
                vcf.push_str(genotype);
            }
            vcf.push('\n');
        }
        let study = GwasInput {
            phenotype: &THE_TRAIT_OF_TWO_SUBPOPULATIONS,
            trait_type: TraitType::Continuous,
            design: &THE_DESIGN_OF_TWO_SUBPOPULATIONS,
            num_coefs: 2,
            kinship: None,
            test: None,
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS_OF_TWO_SUBPOPULATIONS,
            transform_to_biallelic: false,
        };
        let mut reader = reader_over(vcf.as_bytes());
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of two subpopulations: {error}"),
        };

        assert_eq!(result.num_vars, 2);
        // The variant is in the result with its frequency, as every
        // variant that has no answer is.
        assert_within(result.allele_freq[0], 0.5, 1e-12, "the frequency of v0");
        assert!(
            result.beta[0].is_nan() && result.se[0].is_nan() && result.p_value[0].is_nan(),
            "the variant the design leaves nothing of was answered with a beta of {beta}, \
             an se of {se} and a p-value of {p_value}",
            beta = result.beta[0],
            se = result.se[0],
            p_value = result.p_value[0]
        );
        let (beta, se) = OF_THE_ORDINARY_VARIANT;
        assert_within(result.beta[1], beta, 1e-12, "the effect of v1");
        assert_within(result.se[1], se, 1e-12, "the standard error of v1");
    }

    /// A study whose reader gives no variant is refused, which is what
    /// pyNei does with "There are no variants to test": the null model is
    /// fitted, and there is nothing to test against it.
    #[test]
    fn a_study_of_no_variant_is_refused() {
        let mut reader = reader_over(THE_HEADER.as_bytes());
        let study = the_worked_example_study();
        match the_study_of(&mut reader, &study) {
            Err(Error::PassGaveNoVariant { .. }) => {}
            Err(error) => panic!("the study was refused with {error}"),
            Ok(result) => panic!("a study of {} variants was run", result.num_vars),
        }
    }

    /// A kinship that does not hold one row and one column for each tested
    /// individual is refused, and the refusal comes before the study is
    /// told that its model is not written, which is what "checked before
    /// any model is fitted" of "The Rust interface" of
    /// `docs/specs/gwas.md` asks: the six individuals of the worked
    /// example want 36 values.
    #[test]
    fn a_kinship_that_is_not_of_the_tested_individuals_is_refused() {
        let vcf = the_worked_example_vcf();
        for values in [25_usize, 35, 37, 49] {
            let kinship = vec![0.0_f64; values];
            let study = GwasInput {
                kinship: Some(&kinship),
                ..the_worked_example_study()
            };
            let mut reader = reader_over(&vcf);
            match the_study_of(&mut reader, &study) {
                Err(Error::GwasInputOfAnotherSize {
                    problem:
                        GwasInputShape::Kinship {
                            num_values,
                            num_individuals,
                        },
                }) => {
                    assert_eq!(num_values, values, "the values the kinship held");
                    assert_eq!(num_individuals, 6, "the individuals that are tested");
                }
                Err(error) => panic!("a kinship of {values} values was refused with {error}"),
                Ok(result) => panic!("a study of {} variants was run", result.num_vars),
            }
        }
    }

    /// A kinship that holds a value which is not a finite number is
    /// refused, naming the two individuals of that cell and the value. A
    /// fit would not notice it: the eigendecomposition spreads it through
    /// every eigenvalue and every eigenvector, and what the user would be
    /// told is that a matrix of the linear algebra crate is not finite.
    #[test]
    fn a_kinship_that_holds_a_value_that_is_not_a_number_is_refused() {
        let vcf = the_worked_example_vcf();
        for (at, row, column, held) in [
            (0_usize, 0_usize, 0_usize, f64::NAN),
            (13, 2, 1, f64::INFINITY),
            (35, 5, 5, f64::NEG_INFINITY),
        ] {
            let mut kinship = vec![0.0_f64; 36];
            kinship[at] = held;
            let study = GwasInput {
                kinship: Some(&kinship),
                ..the_worked_example_study()
            };
            let mut reader = reader_over(&vcf);
            match the_study_of(&mut reader, &study) {
                Err(Error::GwasKinshipValueNotFinite {
                    individual,
                    other,
                    value,
                }) => {
                    assert_eq!(individual, row, "the row of the value");
                    assert_eq!(other, column, "the column of the value");
                    // `total_cmp` orders every `f64`, NaN among them, so
                    // one comparison covers the three values.
                    assert_eq!(
                        value.total_cmp(&held),
                        std::cmp::Ordering::Equal,
                        "the value is {value} and the kinship held {held} at {at}"
                    );
                }
                Err(error) => panic!("a kinship holding {held} was refused with {error}"),
                Ok(result) => panic!("a study of {} variants was run", result.num_vars),
            }
        }
    }
    /// The t test of a block gives the same three numbers, to the bit, on
    /// the threads of rayon and one variant after another, and on one
    /// thread and on four.
    ///
    /// Every sum of a row runs over the values of that row alone, left to
    /// right, so no answer can depend on how many threads there are. This
    /// is the test that fails the day something is shared between the
    /// rows, and the one that says the answers come back in the order of
    /// the block and not in the order the threads finished in.
    ///
    /// The block is 500 rows of 40 individuals, which is more than one
    /// split of rayon at both thread counts. Every fifth row is a variant
    /// of which the design left nothing, all of its residualized values 0
    /// against dosages that are not, so it gets the three NaNs of a
    /// variant with no answer and that path runs on the threads as well;
    /// the other 400 are answered.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_t_tests_of_a_block_are_the_same_on_threads_and_one_after_another() {
        use super::{TheRowsToTest, the_answers_of_the_rows, the_answers_of_the_rows_one_by_one};

        /// How many individuals the block holds.
        const INDIVIDUALS: usize = 40;
        /// How many rows it holds.
        const ROWS: usize = 500;
        /// How many of those rows have no answer: every fifth.
        const OF_EVERY: usize = 5;

        let dosages: Vec<f64> = (0..ROWS)
            .flat_map(|row| {
                (0..INDIVIDUALS).map(move |individual| {
                    (individual
                        .saturating_mul(7)
                        .saturating_add(row)
                        .rem_euclid(3)) as f64
                })
            })
            .collect();
        let mut residualized = Vec::with_capacity(dosages.len());
        for (row, of_the_variant) in dosages.as_chunks::<INDIVIDUALS>().0.iter().enumerate() {
            for dosage in of_the_variant {
                residualized.push(match row.rem_euclid(OF_EVERY) {
                    0 => 0.0,
                    _ => dosage - 0.75,
                });
            }
        }
        let residuals: Vec<f64> = (0..INDIVIDUALS)
            .map(|individual| (individual as f64) * 0.125 - 2.5)
            .collect();
        let num: Vec<f64> = residualized
            .as_chunks::<INDIVIDUALS>()
            .0
            .iter()
            .map(|row| {
                row.iter()
                    .zip(&residuals)
                    .map(|(value, residual)| value * residual)
                    .sum::<f64>()
            })
            .collect();
        // The squared length of each row of dosages, which the block sums
        // where it writes the row and which the test is given as the block
        // would give it.
        let sum_of_squares: Vec<f64> = dosages
            .as_chunks::<INDIVIDUALS>()
            .0
            .iter()
            .map(|of_the_variant| {
                of_the_variant
                    .iter()
                    .map(|dosage| dosage * dosage)
                    .sum::<f64>()
            })
            .collect();
        let rows = TheRowsToTest {
            residualized: &residualized,
            num: &num,
            sum_of_squares: &sum_of_squares,
            residuals: &residuals,
            num_individuals: INDIVIDUALS,
            degrees_of_freedom: INDIVIDUALS.saturating_sub(3) as f64,
            share_that_is_nothing: crate::gwas::the_share_that_is_nothing(INDIVIDUALS),
        };

        let tested_on = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            pool.install(|| the_answers_of_the_rows(&rows))
        };
        let on_one = tested_on(1);
        let on_four = tested_on(4);
        let one_after_another = the_answers_of_the_rows_one_by_one(&rows);

        /// The three numbers of an answer as the bits that hold them, so
        /// that the NaNs of a variant with no answer compare equal to each
        /// other and to nothing else.
        fn the_bits(answer: &(f64, f64, f64)) -> (u64, u64, u64) {
            (answer.0.to_bits(), answer.1.to_bits(), answer.2.to_bits())
        }

        for (answers, how) in [
            (&on_one, "on one thread"),
            (&on_four, "on four threads"),
            (&one_after_another, "one after another"),
        ] {
            assert_eq!(answers.len(), ROWS, "the answers {how}");
            let with_no_answer = answers.iter().filter(|answer| answer.0.is_nan()).count();
            assert_eq!(
                with_no_answer,
                ROWS.div_euclid(OF_EVERY),
                "the rows with no answer {how}, which are the ones the design left nothing of"
            );
        }
        for row in 0..ROWS {
            let of_one = the_bits(&on_one[row]);
            assert_eq!(
                of_one,
                the_bits(&on_four[row]),
                "the row {row} on one thread and on four"
            );
            assert_eq!(
                of_one,
                the_bits(&one_after_another[row]),
                "the row {row} on one thread and tested one after another"
            );
        }
    }
}
