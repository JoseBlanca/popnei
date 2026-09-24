//! The GRAMMAR-Gamma approximation, which both mixed models can take in
//! place of the denominator of their test.
//!
//! [`GrammarGamma`] is the factor `gamma` and the denominator it gives a
//! variant. "The GRAMMAR-Gamma approximation" of `docs/specs/gwas.md` says
//! what it is for: `x' p x`, with `x` the dosages of a variant and `p` the
//! projection matrix of the null model, costs a product of the variant with
//! an individuals by individuals matrix, so the work of a block grows with
//! the square of the individuals. The approximation replaces it with
//! `gamma` times the squared length of the variant's centered dosages,
//! which costs two walks over the variant and grows with the individuals
//! alone.
//!
//! The factor is estimated once, from the first block of a second pass over
//! the same variants, and it is the mean over the first
//! [`NUM_VARS_FOR_GAMMA`] variants of that block which vary of the exact
//! denominator divided by the approximate one. A variant of those whose
//! exact denominator is nothing but the rounding of a cancellation is left
//! out of that mean, by the rule of **Open 2** of `docs/specs/gwas.md` and
//! against the same scale the model's own score test judges a denominator
//! by.

use popnei_linalg::{TheFirstOperand, TheSecondOperand};

use crate::error::{Error, Result};

use super::dosages::GwasDosages;
use super::the_share_that_is_nothing;

/// How many variants of the first block the factor is the mean of the
/// ratios of: 100.
///
/// It is `NUM_VARS_FOR_GAMMA` of `pynei/gwas.py`, inherited, and nobody has
/// measured whether 100 is the right number, which "What it gives" of the
/// approximation in `docs/specs/gwas.md` records. What was measured on 24
/// September 2026, on the panel with every genotype called, 200 individuals
/// and 1200 variants, is how far apart the 100 ratios lie: under the linear
/// mixed model they run from 0.312 to 0.670 around a mean of 0.517, a
/// standard deviation of 0.0665, which is 12.9 per cent of the mean. So the
/// count is not averaging out a spread that has already collapsed, and a
/// larger one would move the factor.
pub(crate) const NUM_VARS_FOR_GAMMA: usize = 100;

/// The factor a mixed model that makes the approximation multiplies the
/// squared length of a variant's centered dosages by.
///
/// It stands for a quantity that really differs from variant to variant,
/// which is where the accuracy goes: the ratios it is the mean of span a
/// factor of 2.1 on the panel of `docs/specs/gwas.md`, and a panel with
/// more structure spreads them further.
#[derive(Debug, Clone, Copy)]
pub(crate) struct GrammarGamma {
    /// The mean of the ratios, above 0 and finite, which
    /// [`GrammarGamma::of_the_first_block`] is what checks.
    factor: f64,
}

impl GrammarGamma {
    /// The factor estimated from the variants of `dosages` that vary, which
    /// are the rows of the first block of the second pass, against the
    /// projection matrix `projection` of the null model that was fitted.
    ///
    /// `projection` holds `num_individuals` rows of `num_individuals`
    /// values, row after row, and `dosages` holds one row of that many
    /// values for each variant of the block that has variance. The first
    /// [`NUM_VARS_FOR_GAMMA`] of those rows are taken, or all of them when
    /// they are fewer, and the factor is the mean over them of `x' p x`
    /// divided by the squared length of the variant's centered dosages.
    ///
    /// A variant of those whose `x' p x` is at most the tested individuals
    /// times 2.2e-16 of what there was is left out of the mean. What there
    /// was is the variant's own squared length times
    /// `largest_of_the_projection`, the largest value of the diagonal of
    /// the projection matrix, which is the scale the score test of both
    /// mixed models judges a denominator by and which the caller passes in
    /// because it is the caller that holds it. Such a variant is one the
    /// design explains, and what is left of it is the rounding of a
    /// cancellation, which falls on either side of 0: averaging it in gives
    /// a factor of about 1e-16 whenever the rounding falls positive, which
    /// puts every denominator of the study under the threshold of
    /// **Open 2** of `docs/specs/gwas.md` and leaves every variant with the
    /// three NaNs and nothing said.
    ///
    /// # Errors
    ///
    /// [`Error::GwasGrammarGammaWithoutAVariantThatVaries`] when no variant
    /// of the block has any variance among the tested individuals, so there
    /// is no ratio to average.
    /// [`Error::GwasGrammarGammaFactorNotAboveZero`] when the mean of the
    /// ratios that were kept is not a finite number above 0, and when no
    /// ratio was kept at all, where the mean of nothing is NaN: both are
    /// what a block of variants the design explains gives.
    /// [`Error::GwasVariantsTooLarge`] when the values of
    /// those rows are more than a `usize` counts, and
    /// [`Error::GwasLinalg`] when the product of them with the projection
    /// matrix could not be done, which is where a block of other
    /// individuals than the null model was fitted over is refused.
    pub(crate) fn of_the_first_block(
        projection: &[f64],
        num_individuals: usize,
        largest_of_the_projection: f64,
        dosages: &GwasDosages,
    ) -> Result<GrammarGamma> {
        let num_vars = dosages.num_with_variance().min(NUM_VARS_FOR_GAMMA);
        if num_vars == 0 {
            return Err(Error::GwasGrammarGammaWithoutAVariantThatVaries { num_individuals });
        }
        let values = num_vars
            .checked_mul(num_individuals)
            .ok_or(Error::GwasVariantsTooLarge)?;
        // The rows of the variants that have variance are at the start of
        // the buffer, so the first `num_vars` of them are the first
        // `num_vars` variants of the block that vary, which is what pyNei's
        // `estimate_gamma` takes.
        let rows = dosages.dosages().get(..values).unwrap_or_default();
        if rows.len() != values {
            // The block has fewer values than the variants it says have
            // variance, which is a defect of the dosages of this module.
            return Err(Error::GwasVariantsTooLarge);
        }
        let mut projected = vec![0.0_f64; values];
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: rows,
                rows: num_vars,
            },
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: projection,
                cols: num_individuals,
            },
            &mut projected,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the first variants that vary with the projection matrix",
            source,
        })?;
        // The share of what the variant was that the projection has to
        // leave of it for its ratio to be a ratio and not the rounding of a
        // cancellation over a squared length. It is the rule of the score
        // test of both mixed models, at the same scale, and a variant that
        // does not reach it is one the design explains.
        let share_that_is_nothing = the_share_that_is_nothing(num_individuals);
        // At most `NUM_VARS_FOR_GAMMA` of them, and this runs once for a
        // study and not once for a block.
        let the_ratios: Vec<f64> = projected
            .chunks_exact(num_individuals.max(1))
            .zip(rows.chunks_exact(num_individuals.max(1)))
            .filter_map(|(row, of_the_variant)| {
                let exact = row
                    .iter()
                    .zip(of_the_variant)
                    .map(|(projected, dosage)| projected * dosage)
                    .sum::<f64>();
                let of_the_dosages = of_the_variant
                    .iter()
                    .map(|dosage| dosage * dosage)
                    .sum::<f64>();
                if exact <= share_that_is_nothing * largest_of_the_projection * of_the_dosages {
                    return None;
                }
                Some(exact / the_centered_squared_length_of(of_the_variant))
            })
            .collect();
        // The mean of no ratio at all is NaN, which is the refusal below
        // and is what a block of nothing but variants the design explains
        // gives.
        let factor = the_ratios.iter().sum::<f64>() / the_ratios.len() as f64;
        if !factor.is_finite() || factor <= 0.0 {
            return Err(Error::GwasGrammarGammaFactorNotAboveZero { factor, num_vars });
        }
        Ok(GrammarGamma { factor })
    }

    /// The factor itself, which the cargo tests of both mixed models assert
    /// against the one pyNei's `estimate_gamma` gives for the same panel
    /// and the same model.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the factor is in no result and the tests that compare it with \
                      pyNei's are the one thing that reads it; what a study uses is \
                      `den_of`"
        )
    )]
    pub(crate) fn factor(&self) -> f64 {
        self.factor
    }

    /// The denominator of the test of the variant whose dosages over the
    /// tested individuals are `of_the_variant`: the factor times the
    /// squared length of those dosages once their own mean is taken out.
    ///
    /// It stands in for `x' p x`, the variant through the projection matrix
    /// of the null model, and a variant that has variance among those
    /// individuals gets a value above 0 from it whatever that matrix would
    /// have left of the variant, which "Open 2's threshold under the
    /// approximation" of `docs/specs/gwas.md` is about.
    #[must_use]
    pub(crate) fn den_of(&self, of_the_variant: &[f64]) -> f64 {
        self.factor * the_centered_squared_length_of(of_the_variant)
    }
}

/// The squared length of `of_the_variant` once the mean of its values is
/// taken out of each of them, which is what the factor multiplies.
///
/// It is `x - x.mean()` and then the sum of the squares, as
/// `_num_and_den` of `pynei/gwas.py` writes it. A variant that has variance
/// among the tested individuals has two called genotypes of different
/// dosage, and a genotype with an allele missing holds the mean of the
/// called ones, so the values are not all equal and this is above 0. The
/// terms are squares and none of them is below 0, so nothing cancels here.
fn the_centered_squared_length_of(of_the_variant: &[f64]) -> f64 {
    let num_individuals = of_the_variant.len();
    if num_individuals == 0 {
        return 0.0;
    }
    let mean = of_the_variant.iter().sum::<f64>() / num_individuals as f64;
    of_the_variant
        .iter()
        .map(|dosage| (dosage - mean) * (dosage - mean))
        .sum()
}

/// The factor of a projection matrix and a block worked out by hand, and
/// the two blocks it cannot be estimated from.
#[cfg(test)]
mod tests {
    use super::{GrammarGamma, NUM_VARS_FOR_GAMMA, the_centered_squared_length_of};
    use crate::error::Error;
    use crate::gwas::dosages::{BlockOfThePass, GwasDosages};
    use crate::gwas::fixtures::{
        OF_EIGHT, TESTED_OF_EIGHT, THE_PANEL_OF_EIGHT, a_block, a_study, the_design_of,
        the_phenotype_and_the_design_of,
    };
    use crate::gwas::study::GwasInput;

    /// How far the factor of these tests may be from the literal beside it:
    /// 1e-15 of it, relative.
    ///
    /// The values are ratios of small sums of whole dosages and of products
    /// of them with a projection matrix of eight entries, all of them near
    /// 1, so a share of the value and a share of the scale are the same
    /// bound. Measured on 24 September 2026 on Accelerate and on faer with
    /// the bound set to 0: every value is the literal to the bit, since the
    /// arithmetic is sums and products of whole numbers below 2^53 and two
    /// divisions.
    const OF_A_FACTOR: f64 = 1e-15;

    /// The dosages of the four variants of the panel of eight over all
    /// eight individuals, all four of which have variance there. Over the
    /// four individuals that are tested in the rest of these tests, the
    /// first of them has none.
    fn the_dosages_of_the_panel_of_eight() -> GwasDosages {
        let phenotype = [1.0_f64, 3.0, 2.0, 5.0, 4.0, 7.0, 6.0, 9.0];
        let design: Vec<f64> = (0..8).flat_map(|at| [1.0, f64::from(at)]).collect();
        let study = a_study(&phenotype, &design, &THE_PANEL_OF_EIGHT);
        let design = the_design_of(&study, 8);
        let mut block = a_block(8, 2, &OF_EIGHT);
        let mut dosages = GwasDosages::of_a_study();
        dosages
            .read_the_block(
                &mut block,
                &design,
                BlockOfThePass {
                    ploidy: 2,
                    first_var: 0,
                },
            )
            .expect("the dosages of the block");
        dosages
    }

    /// The identity, eight by eight, row after row: the projection matrix
    /// that leaves a variant as it was, so that the exact denominator of a
    /// variant is its own squared length and every ratio can be worked out
    /// by hand.
    fn the_identity_of_eight() -> Vec<f64> {
        (0..8)
            .flat_map(|row| (0..8).map(move |col| f64::from(u8::from(row == col))))
            .collect()
    }

    /// A study of dosages that are not centered has a squared length larger
    /// than the one the approximation takes, and the difference is the mean
    /// of the variant.
    ///
    /// The dosages of `v1` of the panel of eight over all eight are
    /// 0 0 1 0 2 0 1 0, whose squared length is 6 and whose mean is 0.5, so
    /// the centered squared length is `6 - 8 * 0.25 = 4`. It is the one
    /// quantity of the approximation that a reader can check without a
    /// matrix, and the whole of the difference between `x' x` and what the
    /// factor multiplies.
    #[test]
    fn the_centered_squared_length_takes_the_mean_of_the_variant_out() {
        let of_v1 = [0.0_f64, 0.0, 1.0, 0.0, 2.0, 0.0, 1.0, 0.0];

        let found = the_centered_squared_length_of(&of_v1);

        assert!(
            (found - 4.0).abs() <= OF_A_FACTOR * 4.0,
            "the centered squared length of v1 is {found} and the fixture gives 4"
        );
        let of_the_dosages = of_v1.iter().map(|dosage| dosage * dosage).sum::<f64>();
        assert!(
            (of_the_dosages - 6.0).abs() <= OF_A_FACTOR * 6.0,
            "its squared length before the mean is taken out is {of_the_dosages} and the \
             fixture gives 6"
        );
    }

    /// The centered squared length of no value at all is 0 and not NaN.
    ///
    /// The mean of no value is a sum of 0 over a count of 0, which is NaN,
    /// and the early return is what the function answers with in its place.
    /// A study of no individual does not reach here,
    /// `Design::of_the_study` refusing a study of no more individuals than
    /// the columns of its design plus one, so what this holds is the
    /// function itself: an empty row gets a number, which a caller
    /// multiplies by the factor and divides a variant by.
    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "the early return hands back the literal 0.0 and no arithmetic touches \
                  it, so the value is exact and a tolerance would let a NaN through"
    )]
    fn the_centered_squared_length_of_no_value_is_zero() {
        let found = the_centered_squared_length_of(&[]);

        assert_eq!(
            found, 0.0,
            "the centered squared length of no value is {found}"
        );
    }

    /// Against the identity the exact denominator of a variant is its own
    /// squared length, so every ratio is the squared length over the
    /// centered one and the factor is their mean.
    ///
    /// The four variants of the panel of eight, which all have variance
    /// over all eight individuals, with their dosages, their squared
    /// length, their centered squared length and the ratio of the two. A
    /// genotype with an allele missing holds the mean dosage of the called
    /// ones, which is why `v2` and `v3` have a fraction among theirs:
    ///
    /// - `v0`, 1 0 1 0 1 2 1 2, squared length 12, mean 1, centered 4,
    ///   ratio 3;
    /// - `v1`, 0 0 1 0 2 0 1 0, squared length 6, mean 0.5, centered 4,
    ///   ratio 1.5;
    /// - `v2`, 2 0 4/7 0 0 0 2 0, squared length 8.326530612244898, mean
    ///   4/7, centered 5.714285714285714, ratio 1.457142857142857;
    /// - `v3`, 1 0 1 1 1 2 1 1, squared length 10, mean 1, centered 2,
    ///   ratio 5.
    ///
    /// The mean of the four is 2.7392857142857143, and it is above 0, which
    /// is what a factor has to be. What the block of the four shows beside
    /// the arithmetic is that `v2` and `v3`, whose missing genotypes hold
    /// the mean, are centered at the mean of the values the model tests and
    /// not at the mean of the called ones.
    #[test]
    fn the_factor_against_the_identity_is_the_mean_of_the_four_ratios() {
        let dosages = the_dosages_of_the_panel_of_eight();
        assert_eq!(
            dosages.num_with_variance(),
            4,
            "the variants of the panel of eight that vary among all eight"
        );

        let found =
            match GrammarGamma::of_the_first_block(&the_identity_of_eight(), 8, 1.0, &dosages) {
                Ok(estimated) => estimated.factor(),
                Err(error) => panic!("the factor of the panel of eight: {error}"),
            };

        let expected = 2.739_285_714_285_714_3;
        assert!(
            (found - expected).abs() <= OF_A_FACTOR * expected,
            "the factor of the panel of eight against the identity is {found} and the \
             four ratios give {expected}"
        );
        let of_v1 = [0.0_f64, 0.0, 1.0, 0.0, 2.0, 0.0, 1.0, 0.0];
        let den = GrammarGamma::of_the_first_block(&the_identity_of_eight(), 8, 1.0, &dosages)
            .expect("the factor")
            .den_of(&of_v1);
        assert!(
            (den - expected * 4.0).abs() <= OF_A_FACTOR * expected * 4.0,
            "the denominator of v1 is {den} and the factor times its centered squared \
             length of 4 is {}",
            expected * 4.0
        );
    }

    /// A block in which no variant varies among the tested individuals is
    /// refused, naming how many they are: there is no ratio to average, and
    /// the mean of nothing would be a factor of NaN that every variant of
    /// the study would then be divided by.
    ///
    /// It is pyNei's refusal in `estimate_gamma` of `pynei/gwas.py`. The
    /// fixture is the first variant of the panel of eight twice, which has
    /// variance over all eight and none over the four that are tested here,
    /// where every one of them is heterozygous.
    #[test]
    fn a_first_block_in_which_nothing_varies_is_refused() {
        // `v0` of the panel of eight has variance over all eight and none
        // over the four that are tested here, where every one of them is
        // heterozygous.
        let (phenotype, values) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &values, &TESTED_OF_EIGHT);
        let design = the_design_of(&study, 8);
        let with_no_variance: Vec<&[i8]> = vec![OF_EIGHT[0], OF_EIGHT[0]];
        let mut block = a_block(8, 2, &with_no_variance);
        let mut dosages = GwasDosages::of_a_study();
        dosages
            .read_the_block(
                &mut block,
                &design,
                BlockOfThePass {
                    ploidy: 2,
                    first_var: 0,
                },
            )
            .expect("the dosages of the block");
        assert_eq!(dosages.num_with_variance(), 0, "nothing varies among four");

        let projection = vec![0.0_f64; 16];
        match GrammarGamma::of_the_first_block(&projection, 4, 0.0, &dosages) {
            Err(Error::GwasGrammarGammaWithoutAVariantThatVaries { num_individuals }) => {
                assert_eq!(num_individuals, 4);
            }
            other => panic!("a block in which nothing varies gave {other:?}"),
        }
    }

    /// A block of nothing but variants the projection leaves nothing of is
    /// refused, with the factor and how many variants it was the mean over.
    ///
    /// The projection matrix of a null model whose design explains a
    /// variant leaves that variant at the rounding of a cancellation, which
    /// falls on either side of 0; the fixture is the extreme of that, a
    /// projection of all zeros, which leaves every exact denominator
    /// exactly 0. Every one of the four is at the threshold and is left out
    /// of the mean, so the factor is the mean of no ratio at all, which is
    /// NaN, and the study is refused. Unrefused with the rounding fallen
    /// positive, the factor would be about 1e-16, every denominator of the
    /// study would sit under the threshold of **Open 2** and the whole
    /// column would be NaN with nothing to say why.
    ///
    /// The count in the message is the four variants the ratios were formed
    /// over and not the none of them that were kept: it is what tells a
    /// user which variants to look at.
    #[test]
    fn a_block_of_variants_the_projection_leaves_nothing_of_is_refused() {
        let dosages = the_dosages_of_the_panel_of_eight();
        let projection = vec![0.0_f64; 64];

        match GrammarGamma::of_the_first_block(&projection, 8, 0.0, &dosages) {
            Err(Error::GwasGrammarGammaFactorNotAboveZero { factor, num_vars }) => {
                assert!(
                    factor.is_nan(),
                    "the factor a projection of zeros gives is the mean of no ratio, \
                     and it is {factor}"
                );
                assert_eq!(num_vars, 4, "the variants it was the mean over");
            }
            other => panic!("a projection of all zeros gave {other:?}"),
        }
    }

    /// A factor that is not a finite number is refused, and a study of
    /// every variant divided by NaN is not given.
    ///
    /// The fixture is the identity times the largest `f64`, against which
    /// every exact denominator overflows to an infinity while the scale the
    /// ratios are judged by is 1, so every one of the four ratios is kept
    /// and their mean is an infinity. It is the half of the check that the
    /// refusal of a projection of zeros does not reach: that one is a mean
    /// of no ratio at all, and this one is a mean of four.
    ///
    /// No panel gives it. What it holds is the `is_finite` of the check: a
    /// factor of NaN would make `beta`, `se` and the p-value of every
    /// variant of the study NaN, and a factor of an infinity would make
    /// every `beta` 0 and every p-value 1.
    #[test]
    fn a_factor_that_is_not_finite_is_refused() {
        let dosages = the_dosages_of_the_panel_of_eight();
        let projection: Vec<f64> = the_identity_of_eight()
            .iter()
            .map(|value| value * f64::MAX)
            .collect();

        match GrammarGamma::of_the_first_block(&projection, 8, 1.0, &dosages) {
            Err(Error::GwasGrammarGammaFactorNotAboveZero { factor, num_vars }) => {
                assert!(
                    !factor.is_finite() && factor > 0.0,
                    "the factor the identity times the largest f64 gives is {factor}, \
                     which is above 0 and is where the other half of the check does \
                     nothing"
                );
                assert_eq!(num_vars, 4, "the variants it was the mean over");
            }
            other => panic!("a projection of infinities gave {other:?}"),
        }
    }

    /// A block whose values are fewer than its variants that vary times the
    /// individuals is refused.
    ///
    /// The rows the ratios are formed from are the first `num_vars` rows of
    /// `num_individuals` values of the buffer, and a buffer that does not
    /// hold that many is a defect of the dosages of this module. The
    /// fixture makes it by asking for the factor of the four variants of
    /// the panel of eight over nine individuals: the buffer holds 4 rows of
    /// 8 and the estimate wants 4 of 9. Without the refusal the rows would
    /// be empty, and the mean of no ratio would be NaN, which is the same
    /// refusal one step later and with a message about the design of the
    /// study instead of about popnei.
    #[test]
    fn a_block_of_fewer_values_than_its_variants_and_individuals_is_refused() {
        let dosages = the_dosages_of_the_panel_of_eight();
        assert_eq!(
            dosages.dosages().len(),
            32,
            "four variants of eight individuals"
        );
        let projection = vec![0.0_f64; 81];

        match GrammarGamma::of_the_first_block(&projection, 9, 1.0, &dosages) {
            Err(Error::GwasVariantsTooLarge) => {}
            other => panic!("a block of four rows of eight read as nine gave {other:?}"),
        }
    }

    /// The factor is the mean over the first 100 variants that vary and not
    /// over the block, and a block of fewer than 100 that vary gives it
    /// from those.
    ///
    /// The panel of eight has four, so the second half is what its factor
    /// above shows. What this asserts is the count itself, which no fixture
    /// of this module is large enough to reach: a block of 100 variants
    /// that vary and one of 101 give the same factor when the 101st is
    /// different from the others, and a block of 99 gives another.
    #[test]
    fn the_factor_is_the_mean_over_the_first_hundred_that_vary() {
        assert_eq!(NUM_VARS_FOR_GAMMA, 100, "pyNei's NUM_VARS_FOR_GAMMA");
        let phenotype: Vec<f64> = (0..8).map(|at| 1.0 + f64::from(at)).collect();
        let design: Vec<f64> = (0..8).flat_map(|at| [1.0, f64::from(at % 3)]).collect();
        let study: GwasInput<'_> = a_study(&phenotype, &design, &THE_PANEL_OF_EIGHT);
        let design = the_design_of(&study, 8);
        // `v1` has a ratio of 1.5 against the identity and `v0` one of 3,
        // so a block of 100 `v1`s and one of 100 `v1`s followed by any
        // number of `v0`s both give 1.5, and a block of 99 `v1`s and one
        // `v0` gives 1.515.
        let of_a_hundred: Vec<&[i8]> = std::iter::repeat_n(OF_EIGHT[1], 100).collect();
        let of_a_hundred_and_more: Vec<&[i8]> = of_a_hundred
            .iter()
            .copied()
            .chain(std::iter::repeat_n(OF_EIGHT[0], 20))
            .collect();
        let of_ninety_nine: Vec<&[i8]> = std::iter::repeat_n(OF_EIGHT[1], 99)
            .chain(std::iter::once(OF_EIGHT[0]))
            .collect();
        let identity = the_identity_of_eight();
        for (what, rows, expected) in [
            ("a hundred that vary", of_a_hundred, 1.5),
            ("more than a hundred", of_a_hundred_and_more, 1.5),
            ("ninety nine and one other", of_ninety_nine, 1.515),
        ] {
            let mut block = a_block(8, 2, &rows);
            let mut dosages = GwasDosages::of_a_study();
            dosages
                .read_the_block(
                    &mut block,
                    &design,
                    BlockOfThePass {
                        ploidy: 2,
                        first_var: 0,
                    },
                )
                .expect("the dosages of the block");
            let found = match GrammarGamma::of_the_first_block(&identity, 8, 1.0, &dosages) {
                Ok(estimated) => estimated.factor(),
                Err(error) => panic!("the factor of {what}: {error}"),
            };
            assert!(
                (found - expected).abs() <= OF_A_FACTOR * expected,
                "the factor of a block of {what} is {found} and the ratios give {expected}"
            );
        }
    }
}
