//! The association study of a dataset: which of its variants go with a
//! trait of the individuals, with the effect of each variant on the trait,
//! the uncertainty of that effect, and the chance of an effect at least
//! that large in a dataset where the variant has none, which is its
//! p-value. `docs/specs/gwas.md` says what a study gives and how each of
//! its models is verified.
//!
//! [`calc_gwas`], of `pass`, is the study itself: it fits the null model
//! once, over the individuals that are tested and the design they were
//! given, and then reads the variants in one pass, giving one row for each
//! of them. Of the four models three are written: the linear one and the
//! linear mixed one, which are a continuous trait without and with a
//! kinship, and the score test of the logistic one, which is a binomial
//! trait without one. The Wald test of the logistic model and the logistic
//! mixed model are refused until they are written.
//!
//! `distributions` holds the two functions that turn the statistic of a
//! test into a p-value, which every model of the module ends in, and
//! `study` the design every model is fitted on. Each of the two is a
//! survival function, the chance that a distribution is beyond a point,
//! which `sf` in their names stands for: `chi2_sf_1df` is the chi square
//! with one degree of freedom that every score test and the Wald test of a
//! logistic model need, and `t_sf_two_sided` the Student t of the Wald
//! test of a linear model and of a linear mixed model, which is written
//! from the regularized incomplete beta function below it.
//!
//! [`Design::of_the_study`](study::Design::of_the_study) is what a study
//! is refused for before any model is fitted or any variant is read: an
//! individual to test that the source has not, one that is there twice,
//! individuals that are not in the order the source has them, a study of
//! no more individuals than the columns of its design plus one, a
//! phenotype that is not a number or that does not fit the trait, and a
//! design whose columns are not independent. The
//! [`the_model_and_the_test`](study::the_model_and_the_test) chooses which
//! of the four models a study fits and which of the two tests it makes of
//! every variant, and refuses the two pairs of trait, kinship and test
//! that have no test between them. [`GwasDosages`](dosages::GwasDosages),
//! of `dosages`, turns one block into a dosage per tested individual per
//! variant, computed over those individuals and not over the whole panel,
//! with the frequency of each variant among them and whether it has any
//! variance there. [`Gwas`], of `result`, is what a study gives back, one
//! row per variant, which
//! [`Gwas::add_the_block`](result::Gwas::add_the_block) fills block by
//! block with the answers of the variants that have variance and the
//! three NaNs of the ones that have none.
//!
//! [`LinearModel`](linear::LinearModel), of `linear`, is the first model
//! that gives those answers: the thin QR of the design, the residuals of
//! the trait and the t test of every variant against them.
//! [`LinearMixedModel`](linear_mixed::LinearMixedModel), of
//! `linear_mixed`, is the second: the eigendecomposition of the kinship,
//! the search of `RemlSearch` for the two variances, the effects of the
//! intercept and the covariates, and the projection matrix every variant
//! is taken through, by the Wald test or by the score test.
//!
//! [`LogisticModel`](logistic::LogisticModel), of `logistic`, is the
//! third: the null model of a binomial trait fitted by iteratively
//! reweighted least squares, and the score test of every variant against
//! the weights and the residuals that fit left. Its Wald test, which fits
//! one logistic regression per variant, is being written.
//!
//! The logistic mixed model is being written too, and it is a module of
//! its own beside these three: it takes its design from `study`, its
//! dosages from `dosages` and its p-value from `distributions`, it fills
//! the [`Gwas`] of `result` with the answers of a block, and `pass` is
//! what fits it and reads the blocks through it.
//! `the_share_that_is_nothing` is here and not in a model because every
//! model reads it.

mod distributions;
mod dosages;
mod linear;
mod linear_mixed;
mod logistic;
mod pass;
mod result;
mod study;

pub use distributions::{chi2_sf_1df, t_sf_two_sided};
pub use pass::calc_gwas;
pub use result::{DEFAULT_USE_GRAMMAR_GAMMA_APPROX, Gwas, NullModel};
pub use study::{GwasInput, GwasInputShape, GwasModel, TestType, TraitType};

/// The share of what a quantity was that has to be left of it for it to be
/// worth testing, which is the threshold of the meanwhile of **Open 2** of
/// `docs/specs/gwas.md`: the tested individuals times the distance from 1
/// to the next `f64`.
///
/// The rounding of a sum of `n` products is about `n` times 2.2e-16 times
/// the largest term of the sum, so a quantity that has fallen to that share
/// of the scale it was formed from is the rounding of a cancellation and
/// not a quantity. It is used in the four places that spec item names, each
/// with the scale that place is judged against:
///
/// - `xx` in the linear model, the variant with the covariates taken out
///   of it, against the variant's own squared length;
/// - `x' w x` less `(x' w d) (d' w d)⁻¹ (d' w x)` in the logistic model's
///   score test, against `x' w x`, the weighted squared length the variant
///   had before the covariates were taken out;
/// - `x' p x` in the score test of a mixed model, against the variant's
///   squared length times the largest value of the diagonal of the
///   projection matrix;
/// - `y' p y` less `num² / den` in the linear mixed model's Wald test,
///   against `y' p y` itself, which is the one of the four that is what a
///   variant leaves of the trait where the other three are what the design
///   leaves of the variant.
///
/// It is written here once because the logistic mixed model adds a second
/// caller of the third of them, and because a threshold that differed
/// between the places would be a rule with four answers.
fn the_share_that_is_nothing(num_individuals: usize) -> f64 {
    num_individuals as f64 * f64::EPSILON
}

/// The panels the tests of the dosages and of the result are read over,
/// and what every one of them is asserted against.
///
/// Each one has individuals that are not tested, which is what a phenotype
/// that leaves somebody out gives and the case that tells a study over the
/// tested individuals apart from one over the whole panel.
#[cfg(test)]
mod fixtures {
    use super::dosages::BlockOfThePass;
    use super::study::{Design, GwasInput, TraitType};
    use crate::block::Block;
    use crate::variant::MISSING_ALLELE;

    /// An allele that was not called, which is what a `.` of a VCF is.
    pub(super) const MISSING: i8 = MISSING_ALLELE;

    /// How far a value of these tests may be from the literal beside it:
    /// 1e-15 absolute. The two kinds of value it is used on are a
    /// frequency, which lies between 0 and 1, and a dosage, which lies
    /// between 0 and the ploidy, so the bound is 1e-15 of the largest
    /// either can hold and not a share of a value that may be near 0,
    /// which is what "How it is verified" of `docs/specs/gwas.md` asks a
    /// tolerance to be measured against.
    ///
    /// Every one of these numbers is a sum of whole dosages divided by a
    /// count and then by the ploidy, two divisions of numbers below 2^53,
    /// which IEEE 754 rounds the same way everywhere, so the distance is
    /// expected to be 0 and was measured at 0 over every value of every
    /// test here on 23 September 2026, with the bound itself set to 0, on
    /// Accelerate and on faer. What the bound is for is a literal that a
    /// person reads and writes as sixteen digits.
    pub(super) const OF_A_FREQUENCY: f64 = 1e-15;

    /// A study of a continuous trait over the individuals at `individuals`
    /// of the source, with no kinship, the default test and the design of
    /// an intercept and one covariate.
    pub(super) fn a_study<'a>(
        phenotype: &'a [f64],
        design: &'a [f64],
        individuals: &'a [usize],
    ) -> GwasInput<'a> {
        GwasInput {
            phenotype,
            trait_type: TraitType::Continuous,
            design,
            num_coefs: 2,
            kinship: None,
            test: None,
            use_grammar_gamma_approx: false,
            individuals,
            transform_to_biallelic: false,
        }
    }

    /// The same study of the individuals at `individuals`, with a
    /// phenotype and a design of one value and one row for each of them.
    /// The values are not read by anything these tests call, which reads
    /// the genotypes and the individuals; they are here because a study
    /// carries them and because [`super::Design::of_the_study`] is asked
    /// for the design of the first fixture below.
    pub(super) fn the_phenotype_and_the_design_of(individuals: &[usize]) -> (Vec<f64>, Vec<f64>) {
        let phenotype: Vec<f64> = individuals
            .iter()
            .map(|individual| 0.25 + *individual as f64)
            .collect();
        let design: Vec<f64> = individuals
            .iter()
            .flat_map(|individual| [1.0, 0.5 - *individual as f64])
            .collect();
        (phenotype, design)
    }

    /// The design of a study whose source has `num_individuals_of_the_source`
    /// individuals, checked, which is what the dosages of a block are read
    /// over: it carries the tested individuals and the choice about a
    /// variant of more than two alleles.
    pub(super) fn the_design_of<'a>(
        study: &GwasInput<'a>,
        num_individuals_of_the_source: usize,
    ) -> Design<'a> {
        Design::of_the_study(study, num_individuals_of_the_source)
            .expect("the design of the tested individuals")
    }

    /// Where the one block of a fixture sits in the pass that gives it:
    /// the ploidy the reader says its source has, and the first variant.
    pub(super) fn the_first_block_of(ploidy: usize) -> BlockOfThePass {
        BlockOfThePass {
            ploidy,
            first_var: 0,
        }
    }

    /// A block of `rows.len()` variants over `num_individuals`
    /// individuals of `ploidy`, with the genotypes of each variant given
    /// as one row.
    pub(super) fn a_block(num_individuals: usize, ploidy: usize, rows: &[&[i8]]) -> Block {
        Block {
            num_vars: rows.len(),
            num_individuals,
            ploidy,
            gts: rows.concat(),
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        }
    }

    /// Every value of `found` within [`OF_A_FREQUENCY`] of the literal
    /// beside it, and as many values as there are literals.
    pub(super) fn assert_the_values(found: &[f64], expected: &[f64], what: &str) {
        assert_eq!(found.len(), expected.len(), "how many values {what} holds");
        for (position, (found, expected)) in found.iter().zip(expected).enumerate() {
            assert!(
                (found - expected).abs() <= OF_A_FREQUENCY,
                "{what}: the value {position} is {found} and the literal is {expected}"
            );
        }
    }

    /// The four variants of a panel of eight diploid individuals, of which
    /// the four at 0, 2, 4 and 6 have a phenotype and are tested.
    ///
    /// Each variant is one of the four things the individuals that are
    /// left out can do to a study that read them:
    ///
    /// - `v0`, where every tested individual is heterozygous: it has no
    ///   variance among them and cannot be tested, and the four that are
    ///   left out hold `0/0` and `1/1`, so over the panel it has variance
    ///   and would be tested.
    /// - `v1`: its frequency is 0.5 over the tested individuals and 0.25
    ///   over the panel, because the four that are left out are all `0/0`.
    /// - `v2`, where a tested individual has no genotype: it takes the
    ///   mean dosage of the other three, and the four left out are `1/1`,
    ///   which makes 1 the major allele of the panel and not of the tested
    ///   individuals, so every dosage of the panel is counted from the
    ///   other allele.
    /// - `v3`, where no tested individual has a genotype: its dosages are
    ///   the mean of nothing, which is 0, and the four left out do have
    ///   genotypes and would give it a frequency of 0.5.
    pub(super) const OF_EIGHT: [&[i8]; 4] = [
        // v0: 0/1 0/0 0/1 0/0 0/1 1/1 0/1 1/1
        &[0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1],
        // v1: 0/0 0/0 0/1 0/0 1/1 0/0 0/1 0/0
        &[0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0],
        // v2: 0/0 1/1 ./. 1/1 1/1 1/1 0/0 1/1
        &[0, 0, 1, 1, MISSING, MISSING, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1],
        // v3: ./. 0/0 ./. 0/1 ./. 1/1 ./. 0/1
        &[
            MISSING, MISSING, 0, 0, MISSING, MISSING, 0, 1, MISSING, MISSING, 1, 1, MISSING,
            MISSING, 0, 1,
        ],
    ];

    /// Which of the eight individuals of [`OF_EIGHT`] have a phenotype and
    /// are tested: half of them.
    pub(super) const TESTED_OF_EIGHT: [usize; 4] = [0, 2, 4, 6];

    /// All eight, which is the study the same panel gives when nobody is
    /// left out, and what the tested four are asserted to differ from.
    pub(super) const THE_PANEL_OF_EIGHT: [usize; 8] = [0, 1, 2, 3, 4, 5, 6, 7];

    /// The frequency of each variant of [`OF_EIGHT`] over the four tested
    /// individuals, worked out by hand from the fixture.
    ///
    /// The dosage of a genotype is how many of its alleles are not the
    /// major one of its variant among these four, and the frequency is the
    /// mean dosage of the called genotypes over the ploidy, 2.
    ///
    /// - `v0`: their genotypes are `0/1 0/1 0/1 0/1`, so the alleles 0 and
    ///   1 are called four times each and the major one is 0, the lower
    ///   numbered of two that are equally frequent; the dosages are
    ///   1 1 1 1, their mean is 1 and the frequency is 0.5.
    /// - `v1`: `0/0 0/1 1/1 0/1`, four 0s and four 1s, the major allele 0,
    ///   the dosages 0 1 2 1, their mean 4 / 4 = 1 and the frequency 0.5.
    /// - `v2`: `0/0 ./. 1/1 0/0`, four 0s and two 1s, the major allele 0,
    ///   the dosages 0, none, 2 and 0, their mean 2 / 3 and the frequency
    ///   1 / 3 = 0.3333333333333333.
    /// - `v3`: `./. ./. ./. ./.`, no allele called, so the mean of nothing
    ///   is 0 and the frequency is 0.
    pub(super) const FREQUENCIES_OF_THE_TESTED: [f64; 4] = [0.5, 0.5, 0.3333333333333333, 0.0];

    /// The frequency of each variant of [`OF_EIGHT`] over all eight
    /// individuals, worked out by hand the same way, which is what a study
    /// that read the four individuals without a phenotype would give.
    ///
    /// - `v0`: eight 0s and eight 1s, the major allele 0, the dosages
    ///   1 0 1 0 1 2 1 2, their mean 8 / 8 = 1 and the frequency 0.5, the
    ///   one of the four that does not move; what moves is that the
    ///   variant now has variance.
    /// - `v1`: twelve 0s and four 1s, the major allele 0, the dosages
    ///   0 0 1 0 2 0 1 0, their mean 4 / 8 and the frequency 0.25.
    /// - `v2`: four 0s and ten 1s, so the major allele is 1 here and 0
    ///   among the tested four; the dosages are 2 0, none, 0 0 0 2 0,
    ///   their mean 4 / 7 and the frequency 2 / 7 = 0.2857142857142857.
    /// - `v3`: four 0s and four 1s, the major allele 0, the dosages none,
    ///   0, none, 1, none, 2, none, 1, their mean 4 / 4 = 1 and the
    ///   frequency 0.5.
    pub(super) const FREQUENCIES_OF_THE_PANEL: [f64; 4] = [0.5, 0.25, 0.2857142857142857, 0.5];
}
