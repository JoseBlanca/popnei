//! The linear mixed model: a continuous trait with a kinship.
//!
//! [`LinearMixedModel`] is the eigendecomposition of the kinship, the
//! search of [`RemlSearch`] for the two variances, the effects of the
//! intercept and the covariates, and the projection matrix every variant is
//! taken through, by the Wald test or by the score test. "The linear mixed
//! model" of `docs/specs/gwas.md` says what it fits and how it is
//! verified.

use popnei_linalg::{Eigen, TheFirstOperand, TheSecondOperand};

use crate::error::{Error, Result};

use super::distributions::{chi2_sf_1df, t_sf_two_sided};
use super::dosages::GwasDosages;
use super::result::{Answers, NullModel};
use super::study::{Design, GwasInputShape, GwasModel, TestType};
use super::the_share_that_is_nothing;

/// How many points the search for the ratio of the two variances of a
/// linear mixed model starts with: 101, evenly spaced in the log of that
/// ratio between [`LOG_DELTA_LOWEST`] and 10.
///
/// It is `numpy.linspace(-10, 10, 101)` of `_reml_delta` of
/// `pynei/gwas.py`. "The linear mixed model" of `docs/specs/gwas.md` asks
/// for that search to be reproduced step for step, because the two
/// variances it gives are GMMAT's, and it lists the six things that have
/// to match: these 101 points, the two neighbours of the best clamped at
/// the ends of the grid, the ratio of the golden section, its two interior
/// points, the bracket moved to whichever of the two has the smaller
/// value, and [`GOLDEN_SECTION_STEPS`] steps whatever happens.
const LOG_DELTA_POINTS: usize = 101;

/// The smallest log of that ratio the grid holds, -10, which is a residual
/// variance 22026 times smaller than the genetic one.
const LOG_DELTA_LOWEST: f64 = -10.0;

/// The distance between two points of the grid, which is the 20 it spans
/// over the 100 gaps of its 101 points.
const LOG_DELTA_STEP: f64 = 0.2;

/// How many steps of the golden section search are made after the grid,
/// whatever the bracket has come to: 60, from `_reml_delta` of
/// `pynei/gwas.py`. Each one shrinks the bracket to
/// [`GOLDEN_SECTION_RATIO`] of what it was, so 60 of them take the 0.4 the
/// grid leaves to 1.156e-13.
const GOLDEN_SECTION_STEPS: usize = 60;

/// How near the smallest value of the criterion over the grid has to be
/// to the largest for the two variances to be arbitrary: the points of the
/// grid times the distance from 1 to the next `f64`, as a share of the
/// largest value's own size.
///
/// It is the test of **Open 3** of `docs/specs/gwas.md`, one comparison at
/// the end of the grid. It is not
/// [`the_share_that_is_nothing`], which is a share of the scale a sum was
/// formed from and is about a quantity that cancelled to nothing; this is
/// about a function that never varied, and what it counts is the points
/// that were evaluated and not the individuals.
///
/// Measured on 24 September 2026 on a kinship close to a multiple of the
/// identity: the criterion spans 1.1e-12 over the 101 points, where this
/// allows 2.2e-14 of a value of about 221, which is 4.9e-12.
const THE_SPAN_OF_A_FLAT_CRITERION: f64 = LOG_DELTA_POINTS as f64 * f64::EPSILON;

/// Whether the restricted maximum likelihood told the two variances of a
/// trait apart.
///
/// For a kinship close to a multiple of the identity the model is the
/// ordinary linear one whatever the split between them, so the criterion is
/// flat and which point of the grid wins is rounding. `beta` and `p_value`
/// are untouched, the test being scale free; what is arbitrary is exactly
/// the two variances and the heritability, which is what a user reads a
/// heritability off.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TheSplitOfTheTrait {
    /// The criterion has a minimum, and the two variances are the fit's.
    Identified,
    /// The criterion is flat over the whole grid: the study is still given
    /// and every variant still answered, and the two variances and the
    /// heritability come back as nothing.
    NotIdentified,
}

/// The ratio a golden section step shrinks its bracket by,
/// `(sqrt(5) - 1) / 2`, which is `golden` of `_reml_delta` of
/// `pynei/gwas.py`.
///
/// It is the one number that makes the two interior points of a bracket
/// reusable: the point the new bracket inherits sits where the next step
/// would have put it. pyNei evaluates both of them again at every step all
/// the same, and popnei reproduces that, so what this ratio decides here is
/// only where the points fall and not how many are evaluated. It is one of
/// the six things "The linear mixed model" of `docs/specs/gwas.md` asks to
/// be reproduced, and it is a named constant because a literal `0.61` in
/// its place changes every number of this model and is caught by nothing
/// that does not have pyNei to compare against.
const GOLDEN_SECTION_RATIO: f64 = 0.618_033_988_749_894_9;

/// The point of the grid at `at`, which is 0 to
/// [`LOG_DELTA_POINTS`] less one.
///
/// The step is multiplied first and the start added after, as
/// `numpy.linspace` does it, so that the points are the ones pyNei
/// searched: the 45th of them is -1.1999999999999993 and not -1.2.
fn the_log_delta_at(at: usize) -> f64 {
    LOG_DELTA_LOWEST + at as f64 * LOG_DELTA_STEP
}

/// The restricted maximum likelihood of a linear mixed model, with the
/// buffers one value of the ratio of its two variances is evaluated in.
///
/// Restricted maximum likelihood is maximum likelihood on the part of the
/// trait that the covariates cannot explain, so that fitting the
/// covariates does not drag the variances down. Only the ratio of the two
/// matters to the search: with `delta` the residual variance over the
/// genetic one, the covariance of the trait is the genetic variance times
/// `k + delta i`, with `k` the kinship and `i` the identity, so the
/// kinship is eigendecomposed once and every value of `delta` then costs
/// one number per individual instead of a matrix.
///
/// With `l` the eigenvalues of the kinship, `u` the design turned by the
/// eigenvectors, `uy` the trait turned by them, and `w = 1 / (l + delta)`
/// one weight per individual, the criterion to minimize is
///
/// ```text
/// dvd   = u' diag(w) u
/// coefs = the solution of `dvd coefs = u' (w * uy)`
/// resid = uy - u coefs
/// quad  = w' (resid * resid)
/// value = sum(log(l + delta)) + (n - c) * log(quad) + log(det(dvd))
/// ```
///
/// with `n` the tested individuals and `c` the columns of the design. It
/// is the criterion of `_reml_delta` of `pynei/gwas.py` and of "The linear
/// mixed model" of `docs/specs/gwas.md`. The log of the determinant comes
/// off the same Cholesky factorization the solve uses, where pyNei takes
/// it from numpy's `slogdet`, which factorizes the same matrix a second
/// way. The two were run against each other in numpy on the panel on 24
/// September 2026: the fitted `delta` lands 2.9e-9 of itself apart, which
/// is what a search whose criterion is flat at its minimum comes to, and
/// the two variances and the three effects then agree with GMMAT's to
/// 1.19e-6 either way.
struct RemlSearch<'a> {
    /// The eigenvalues of the kinship, clamped at 0, one per individual.
    eigenvalues: &'a [f64],
    /// The trait turned by the eigenvectors, one value per individual.
    rotated_trait: &'a [f64],
    /// The design turned by them, `num_individuals` x `num_coefs`, row
    /// after row.
    rotated_design: &'a [f64],
    /// How many individuals are tested.
    num_individuals: usize,
    /// How many columns the design has.
    num_coefs: usize,
    /// The individuals less those columns, as a number, which the criterion
    /// weighs the log of `quad` by.
    degrees_of_freedom_of_the_null: f64,
    /// `1 / (eigenvalue + delta)`, one per individual.
    weights: Vec<f64>,
    /// The rotated design with every row times the weight of its
    /// individual.
    weighted_design: Vec<f64>,
    /// The rotated trait with every value times the weight of its
    /// individual.
    weighted_trait: Vec<f64>,
    /// `u' diag(w) u`, `num_coefs` x `num_coefs`, whose lower half becomes
    /// its Cholesky factorization.
    dvd: Vec<f64>,
    /// The effect of the intercept and of each covariate at this `delta`.
    coefs: Vec<f64>,
    /// The rotated design times those effects, one value per individual.
    fitted: Vec<f64>,
}

impl<'a> RemlSearch<'a> {
    /// The search over a kinship that has been eigendecomposed and a trait
    /// and a design that have been turned by its eigenvectors.
    ///
    /// The three slices hold one value, one value and one row per tested
    /// individual, which [`LinearMixedModel::of_the_study`] has just built
    /// from a checked [`Design`].
    fn of(
        eigenvalues: &'a [f64],
        rotated_trait: &'a [f64],
        rotated_design: &'a [f64],
        num_individuals: usize,
        num_coefs: usize,
        degrees_of_freedom_of_the_null: f64,
    ) -> RemlSearch<'a> {
        // The design holds `num_individuals` x `num_coefs` values in one
        // slice, and a study has two more individuals than its design has
        // columns, so the `num_coefs` x `num_coefs` of `dvd` is smaller
        // than a length that fits in a `usize`.
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "`Design::of_the_study` refuses a study of no more individuals than \
                      the columns of its design plus one, and the design is one slice of \
                      `num_individuals` x `num_coefs` values, so `num_coefs` times itself \
                      is smaller than a length a `usize` holds"
        )]
        let of_the_coefs = num_coefs * num_coefs;
        RemlSearch {
            eigenvalues,
            rotated_trait,
            rotated_design,
            num_individuals,
            num_coefs,
            degrees_of_freedom_of_the_null,
            weights: vec![0.0_f64; num_individuals],
            weighted_design: vec![0.0_f64; rotated_design.len()],
            weighted_trait: vec![0.0_f64; num_individuals],
            dvd: vec![0.0_f64; of_the_coefs],
            coefs: vec![0.0_f64; num_coefs],
            fitted: vec![0.0_f64; num_individuals],
        }
    }

    /// The fit at one `delta`, the residual variance over the genetic one:
    /// it leaves the effects of the intercept and the covariates in
    /// `coefs` and gives back the weighted squared length of what the
    /// design left of the trait, `quad`, and the log of the determinant of
    /// `dvd`.
    ///
    /// # Errors
    ///
    /// [`Error::GwasLinalg`] when one of the three products, the Cholesky
    /// factorization of `dvd` or the solve against it could not be done.
    /// A `dvd` that is not positive definite is among them, which a design
    /// whose columns are independent and a `delta` above 0 do not give.
    fn fit_at(&mut self, delta: f64) -> Result<(f64, f64)> {
        for (weight, eigenvalue) in self.weights.iter_mut().zip(self.eigenvalues) {
            *weight = 1.0 / (eigenvalue + delta);
        }
        for ((row, weight), of_the_design) in self
            .weighted_design
            .chunks_exact_mut(self.num_coefs.max(1))
            .zip(&self.weights)
            .zip(self.rotated_design.chunks_exact(self.num_coefs.max(1)))
        {
            for (value, turned) in row.iter_mut().zip(of_the_design) {
                *value = weight * turned;
            }
        }
        for ((weighted, weight), turned) in self
            .weighted_trait
            .iter_mut()
            .zip(&self.weights)
            .zip(self.rotated_trait)
        {
            *weighted = weight * turned;
        }
        popnei_linalg::product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: self.rotated_design,
                rows: self.num_coefs,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.weighted_design,
                cols: self.num_coefs,
            },
            &mut self.dvd,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the turned design with its weighted self",
            source,
        })?;
        popnei_linalg::product(
            TheFirstOperand::ByTheValuesSummedOver {
                values: self.rotated_design,
                rows: self.num_coefs,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.weighted_trait,
                cols: 1,
            },
            &mut self.coefs,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the turned design with the weighted trait",
            source,
        })?;
        popnei_linalg::cholesky_lower(&mut self.dvd, self.num_coefs).map_err(|source| {
            Error::GwasLinalg {
                operation: "Cholesky factorization of the weighted design",
                source,
            }
        })?;
        let log_determinant =
            popnei_linalg::log_determinant_with_cholesky(&self.dvd, self.num_coefs).map_err(
                |source| Error::GwasLinalg {
                    operation: "log determinant of the weighted design",
                    source,
                },
            )?;
        popnei_linalg::solve_with_cholesky(&self.dvd, self.num_coefs, &mut self.coefs, 1).map_err(
            |source| Error::GwasLinalg {
                operation: "solve of the weighted design against the weighted trait",
                source,
            },
        )?;
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: self.rotated_design,
                rows: self.num_individuals,
            },
            self.num_coefs,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.coefs,
                cols: 1,
            },
            &mut self.fitted,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the turned design with the effects it fitted",
            source,
        })?;
        let quad = self
            .weights
            .iter()
            .zip(self.rotated_trait)
            .zip(&self.fitted)
            .map(|((weight, measured), fitted)| {
                let left = measured - fitted;
                weight * left * left
            })
            .sum::<f64>();
        Ok((quad, log_determinant))
    }

    /// The criterion at one log of `delta`, which the search minimizes.
    ///
    /// # Errors
    ///
    /// Whatever [`RemlSearch::fit_at`] fails with.
    fn criterion(&mut self, log_delta: f64) -> Result<f64> {
        let delta = log_delta.exp();
        let (quad, log_determinant) = self.fit_at(delta)?;
        let of_the_eigenvalues = self
            .eigenvalues
            .iter()
            .map(|eigenvalue| (eigenvalue + delta).ln())
            .sum::<f64>();
        Ok(of_the_eigenvalues + self.degrees_of_freedom_of_the_null * quad.ln() + log_determinant)
    }

    /// The `delta` that minimizes the criterion: the smallest of the 101
    /// points of the grid, bracketed by its two neighbours, and then 60
    /// steps of the golden section search inside that bracket.
    ///
    /// A golden section search shrinks a bracket that holds a minimum by a
    /// constant ratio at each step: the two interior points are `high -
    /// ratio * (high - low)` and `low + ratio * (high - low)`, and the
    /// bracket moves to whichever of the two has the smaller value. The
    /// best point of the grid is bracketed by its neighbours clamped at the
    /// ends, so a minimum at either end of the grid is bracketed by that
    /// end and the one point beside it.
    ///
    /// Both interior points are evaluated at every step, as `_reml_delta`
    /// of `pynei/gwas.py` evaluates them, so the whole search costs 221
    /// evaluations and not 161. The usual form of the search keeps the
    /// value at the interior point the new bracket inherits; that point is
    /// then computed again from the new `low` and `high` and is not the
    /// number that was kept, and the spec asks for pyNei's search step for
    /// step because the two variances it gives are GMMAT's.
    ///
    /// What comes back beside the ratio is whether the grid told the two
    /// variances apart at all, which is the one comparison **Open 3** of
    /// `docs/specs/gwas.md` asks for: the smallest value of the criterion
    /// against the largest, over the same 101 points the grid already
    /// evaluated, so it costs nothing.
    ///
    /// # Errors
    ///
    /// Whatever [`RemlSearch::fit_at`] fails with, at any of the 221
    /// evaluations.
    fn the_delta(&mut self) -> Result<(f64, TheSplitOfTheTrait)> {
        let mut best_at = 0_usize;
        let mut smallest = f64::INFINITY;
        let mut largest = f64::NEG_INFINITY;
        for at in 0..LOG_DELTA_POINTS {
            let value = self.criterion(the_log_delta_at(at))?;
            largest = largest.max(value);
            if value < smallest {
                smallest = value;
                best_at = at;
            }
        }
        // The criterion is negative on the panels and its size is what the
        // share is of, so the largest absolute value of the two is what the
        // span is measured against.
        let of_its_own_size = smallest.abs().max(largest.abs());
        let split = match largest - smallest <= THE_SPAN_OF_A_FLAT_CRITERION * of_its_own_size {
            true => TheSplitOfTheTrait::NotIdentified,
            false => TheSplitOfTheTrait::Identified,
        };
        let mut low = the_log_delta_at(best_at.saturating_sub(1));
        let mut high = the_log_delta_at(
            best_at
                .saturating_add(1)
                .min(LOG_DELTA_POINTS.saturating_sub(1)),
        );
        for _ in 0..GOLDEN_SECTION_STEPS {
            let from_the_top = high - GOLDEN_SECTION_RATIO * (high - low);
            let from_the_bottom = low + GOLDEN_SECTION_RATIO * (high - low);
            if self.criterion(from_the_top)? < self.criterion(from_the_bottom)? {
                high = from_the_bottom;
            } else {
                low = from_the_top;
            }
        }
        Ok((((low + high) / 2.0).exp(), split))
    }
}

/// The linear mixed model of a study fitted without any variant in it: the
/// two variances of the trait, the effects of the intercept and the
/// covariates, and the projection matrix every variant is then tested
/// through.
///
/// Beside the covariates the trait carries a random effect whose
/// covariance is the kinship times a variance, so that two related
/// individuals are expected to resemble each other before any variant is
/// looked at. The covariance of the trait under the null is `v =
/// genetic_variance * k + residual_variance * i`, with `k` the kinship and
/// `i` the identity, and the two variances are estimated by restricted
/// maximum likelihood, which [`RemlSearch`] describes.
///
/// The projection matrix is `p = v⁻¹ - v⁻¹ d (d' v⁻¹ d)⁻¹ d' v⁻¹`,
/// individuals by individuals, with `d` the design: it takes the
/// covariates out of anything it is applied to and weights it by the
/// covariance. [`LinearMixedModel::ypy`] is the trait through it, which is
/// what says the fit reached its optimum, and
/// [`LinearMixedModel::test_the_block`] is where every variant goes
/// through it.
///
/// The projection matrix and the trait through it are kept because every
/// variant is tested against them, and so are the buffers of a block, so
/// that a pass over a million variants asks the machine for them once and
/// allocates nothing for a variant. The matrix is individuals by
/// individuals, which at the 10000 individuals of `docs/objectives.md` is
/// 800 MB, where everything a linear model keeps grows with the columns of
/// the design.
pub(crate) struct LinearMixedModel {
    /// The effect of the intercept and of each covariate, one per column
    /// of the design.
    coefs: Vec<f64>,
    /// The variance of the random effect of the kinship.
    genetic_variance: f64,
    /// What is left over, the variance of the trait that the kinship does
    /// not account for.
    residual_variance: f64,
    /// The genetic variance over the sum of the two, which is the share of
    /// the trait's variance that the kinship explains.
    heritability: f64,
    /// Whether the criterion of the search told the two variances apart.
    /// The three numbers above are the fit's own whatever this says, and
    /// the projection matrix is built from them either way, because a
    /// covariance of the right shape is needed to test a variant; what it
    /// decides is whether a user is given them.
    split: TheSplitOfTheTrait,
    /// The projection matrix `p`, `num_individuals` x `num_individuals`,
    /// row after row, which every variant is tested through.
    projection: Vec<f64>,
    /// The trait through it, `p y`, one value per tested individual, which
    /// every variant's numerator is taken against.
    projected_trait: Vec<f64>,
    /// The largest value of the diagonal of that matrix, which is what a
    /// variant's own squared length is weighted by to say how much of the
    /// variant the projection has left.
    ///
    /// What bounds `x' p x` over `x' x` is the largest eigenvalue of the
    /// projection, and this is not that: the matrix is 0 or above as a
    /// quadratic form, so no value of its diagonal is below 0 and the
    /// largest of them is at most that eigenvalue, which makes this an
    /// under-estimate of the scale the threshold is meant to measure. By
    /// how much was measured on 25 September 2026: the largest eigenvalue
    /// is 1.697 times the largest diagonal entry on `panel_called` and
    /// 1.725 times it on `panel`, so the threshold sits about 1.7 times
    /// below the scale. It is the diagonal that is taken because it costs
    /// one walk over the matrix where the eigenvalue costs a
    /// decomposition, and a threshold under-estimated by 1.7 is a threshold
    /// 1.7 times tighter than it was meant to be, not one that lets a
    /// variant through: on `panel_called` the smallest real denominator is
    /// about 30 against a threshold of 5e-11.
    largest_of_the_projection: f64,
    /// The trait through the projection matrix, `y' p y`.
    ypy: f64,
    /// How many individuals the study tests.
    num_individuals: usize,
    /// The individuals less the columns of the design and one more for the
    /// variant, which is the degrees of freedom of the Wald test. It is 1
    /// at least, since [`Design::of_the_study`] refuses a study of no more
    /// individuals than the columns of its design plus one.
    degrees_of_freedom: f64,
    /// The dosages of a block through the projection matrix, `x p`, the
    /// variants that have variance x `num_individuals`.
    projected: Vec<f64>,
    /// Each variant times the trait through the projection, `x' p y`, one
    /// per variant that has variance.
    num: Vec<f64>,
    /// The effect of each variant that has variance.
    beta: Vec<f64>,
    /// The standard error of each of those effects.
    se: Vec<f64>,
    /// The p-value of each of those tests.
    p_value: Vec<f64>,
}

impl LinearMixedModel {
    /// The linear mixed model of `phenotype` over `design` with `kinship`
    /// as the covariance of its random effect, fitted without any variant
    /// in it.
    ///
    /// `phenotype` holds one value for each individual the design has a
    /// row for, and `kinship` is that many rows of that many values, row
    /// after row, already cut to those individuals and in their order.
    /// [`calc_gwas`](super::calc_gwas) checks the kinship before any model
    /// is fitted, so a caller that comes straight here with a matrix of
    /// another length gets the eigendecomposition's refusal of it.
    ///
    /// The eigenvalues of the kinship are clamped at 0 before use. A
    /// kinship of genotypes with nothing missing has none below 0 but for
    /// rounding, -3.4416913763379853e-15 on the panel of
    /// `docs/specs/gwas.md`, measured with numpy on 25 September 2026; the per
    /// pair denominators of `docs/specs/kinship.md` put them there,
    /// -0.0321 on the panel with 3 in 100 genotypes missing, and a
    /// negative eigenvalue would make the covariance of the trait not a
    /// covariance.
    ///
    /// # Errors
    ///
    /// [`Error::GwasInputOfAnotherSize`] when `phenotype` does not hold
    /// one value for each tested individual, and [`Error::GwasLinalg`]
    /// when the eigendecomposition of the kinship, one of the products,
    /// one of the two Cholesky factorizations or one of the solves could
    /// not be done.
    pub(crate) fn of_the_study(
        phenotype: &[f64],
        design: &Design<'_>,
        kinship: &[f64],
    ) -> Result<LinearMixedModel> {
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
        // The eigendecomposition is taken once and every value of the
        // ratio of the two variances is then one number per individual.
        // The buffer of the matrix comes back as the eigenvectors, so the
        // kinship the caller holds is copied into one of its own.
        let Eigen {
            values: mut eigenvalues,
            vectors,
        } = popnei_linalg::eigh_lower(kinship.to_vec(), num_individuals).map_err(|source| {
            Error::GwasLinalg {
                operation: "eigendecomposition of the kinship",
                source,
            }
        })?;
        for eigenvalue in &mut eigenvalues {
            *eigenvalue = eigenvalue.max(0.0);
        }
        // Row `j` of the eigenvectors is the eigenvector of the eigenvalue
        // `j`, so the product of that matrix with the trait is the trait
        // turned by the eigenvectors, and the same for the design.
        let mut rotated_trait = vec![0.0_f64; num_individuals];
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &vectors,
                rows: num_individuals,
            },
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: phenotype,
                cols: 1,
            },
            &mut rotated_trait,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the kinship's eigenvectors with the trait",
            source,
        })?;
        let mut rotated_design = vec![0.0_f64; design.values().len()];
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: &vectors,
                rows: num_individuals,
            },
            num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: design.values(),
                cols: num_coefs,
            },
            &mut rotated_design,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of the kinship's eigenvectors with the design",
            source,
        })?;
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "`Design::of_the_study` refuses a study of no more individuals than \
                      the columns of its design plus one, so this is 2 at least"
        )]
        let degrees_of_freedom_of_the_null = num_individuals - num_coefs;
        let degrees_of_freedom_of_the_null = degrees_of_freedom_of_the_null as f64;
        let mut search = RemlSearch::of(
            &eigenvalues,
            &rotated_trait,
            &rotated_design,
            num_individuals,
            num_coefs,
            degrees_of_freedom_of_the_null,
        );
        let (delta, split) = search.the_delta()?;
        // The search leaves the weights and the effects at whichever point
        // it evaluated last, so the fit is taken again at the `delta` it
        // gave, as `_LMMNull` of `pynei/gwas.py` does.
        let (quad, _) = search.fit_at(delta)?;
        let genetic_variance = quad / degrees_of_freedom_of_the_null;
        // `quad` is what the design left of the trait, weighted, so a
        // design that explains the whole of it leaves 0 here and both
        // variances are 0. The covariance of the trait is then the zero
        // matrix and the inverse below is infinities, which reached the
        // user as the linear algebra refusing a matrix that is not finite.
        // A NaN is refused here too, which the comparison the other way
        // round would let past: `quad` is a sum of weighted squares and
        // cannot be below 0, so what this tests is 0 and not a sign.
        if genetic_variance.is_nan() || genetic_variance <= 0.0 {
            return Err(Error::GwasDesignExplainsTheTrait);
        }
        let residual_variance = delta * genetic_variance;
        let heritability = genetic_variance / (genetic_variance + residual_variance);
        let coefs = std::mem::take(&mut search.coefs);
        drop(search);
        let projection = the_projection_of(
            design,
            &eigenvalues,
            &vectors,
            TheVariances {
                genetic: genetic_variance,
                residual: residual_variance,
            },
        )?;
        let projected_trait = through_the_projection(&projection, phenotype)?;
        // `y' p y`, the generalized residual sum of squares of the null
        // over the genetic variance, which the restricted maximum
        // likelihood makes the individuals less the columns of the design.
        let ypy = phenotype
            .iter()
            .zip(&projected_trait)
            .map(|(measured, projected)| measured * projected)
            .sum::<f64>();
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "the individuals less the columns of the design are 2 at least, by \
                      the same refusal, so this is 1 at least"
        )]
        let degrees_of_freedom = num_individuals - num_coefs - 1;
        Ok(LinearMixedModel {
            coefs,
            genetic_variance,
            residual_variance,
            heritability,
            split,
            largest_of_the_projection: projection
                .chunks_exact(num_individuals.max(1))
                .zip(0..)
                .filter_map(|(row, at)| row.get(at).copied())
                .fold(0.0_f64, f64::max),
            projection,
            projected_trait,
            ypy,
            num_individuals,
            degrees_of_freedom: degrees_of_freedom as f64,
            projected: Vec::new(),
            num: Vec::new(),
            beta: Vec::new(),
            se: Vec::new(),
            p_value: Vec::new(),
        })
    }

    /// The null model of the result: the effects of the intercept and of
    /// the covariates, the two variances and the share of the trait's
    /// variance that the kinship explains.
    #[must_use]
    pub(crate) fn null_model(&self, test: TestType) -> NullModel {
        // A fit that could not tell the two variances apart gives none of
        // the three numbers built from them, which is the meanwhile of
        // **Open 3** of `docs/specs/gwas.md`: what it would give instead is
        // decided by the last bit of an eigenvalue, and a heritability of
        // 0.967 from one seed and 6.5e-5 from another looks reliable and is
        // not. It is `None` and not 0 and not NaN, because `None` is what
        // the linear model gives for the two numbers it has not, and a user
        // meets one way of saying that a number is not there.
        let of_the_split = match self.split {
            TheSplitOfTheTrait::Identified => Some(()),
            TheSplitOfTheTrait::NotIdentified => None,
        };
        NullModel {
            model: GwasModel::Lmm,
            test,
            covariate_effects: self.coefs.clone(),
            residual_variance: of_the_split.map(|()| self.residual_variance),
            genetic_variance: of_the_split.map(|()| self.genetic_variance),
            heritability: of_the_split.map(|()| self.heritability),
            num_individuals: self.num_individuals,
        }
    }

    /// The trait through the projection matrix, `y' p y`, which is the
    /// generalized residual sum of squares of the null over the genetic
    /// variance.
    ///
    /// The restricted maximum likelihood makes it exactly the individuals
    /// less the columns of the design, so it is the cheapest evidence
    /// there is that the fit reached its optimum. It is in no result, and
    /// "The linear mixed model" of `docs/specs/gwas.md` has the cargo test
    /// of it made here for that reason.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "`y' p y` is in no result and the test of the fit is the one thing \
                      that reads it, which is why the spec has that test made here"
        )
    )]
    pub(crate) fn ypy(&self) -> f64 {
        self.ypy
    }

    /// The Wald test or the score test of every variant of a block that
    /// has variance among the tested individuals, in the order of the
    /// block.
    ///
    /// Every variant goes through the projection matrix of the null, which
    /// takes the covariates out of it and weights it by the covariance of
    /// the trait, so that what two related individuals share counts once
    /// and not twice. With `x` the dosages of a variant, `num` is `x' p y`
    /// and `den` is `x' p x`, and the effect of the variant is `num / den`
    /// in both tests. The two differ in how uncertain that effect is and in
    /// the distribution the effect is read against:
    ///
    /// - The **Wald test** divides the effect by its standard error and
    ///   asks how far a Student t goes beyond that, in both tails, with the
    ///   individuals less the columns of the design and one more for the
    ///   variant. It holds the ratio of the two variances at the null and
    ///   estimates their scale again with the variant in the model, which
    ///   is what `se = sqrt((y' p y - num² / den) / (df * den))` is: what
    ///   the variant leaves of `y' p y`, over the degrees of freedom and
    ///   over `den`. It is what rrBLUP's `GWAS` does with `P3D`.
    /// - The **score test** asks instead how steeply the likelihood rises
    ///   at an effect of 0, so it holds both variances at the null and
    ///   needs no fit with the variant in. The effect is uncertain by `1 /
    ///   sqrt(den)`, and `num² / den` is read against a chi square with one
    ///   degree of freedom. It is what GMMAT's `glmm.score` does.
    ///
    /// The two products are the whole cost of a block: the dosages through
    /// the projection matrix, which is individuals by individuals, and the
    /// dosages against the trait through it. What the GRAMMAR-Gamma
    /// approximation of `docs/specs/gwas.md` replaces is the first of them.
    ///
    /// A variant of which the projection leaves at most the tested
    /// individuals times 2.2e-16 of what there was has no answer, and gets
    /// the three NaNs a variant with no variance gets. What there was is
    /// the variant's own squared length times the largest value of the
    /// diagonal of the projection matrix, which is what bounds `x' p x`
    /// over `x' x`. `den` is 0 or above in exact arithmetic, and a variant
    /// that the covariates and the kinship leave almost nothing of gets the
    /// rounding of that, which can fall below 0: `beta` would be a number
    /// divided by noise, large and of whichever sign the rounding chose,
    /// and the two backends do not choose the same one. It is **Open 2** of
    /// `docs/specs/gwas.md`.
    ///
    /// # Errors
    ///
    /// [`Error::GwasVariantsTooLarge`] when the values of the block are
    /// more than a `usize` counts, and [`Error::GwasLinalg`] when one of
    /// the two products could not be done, which is where a block of other
    /// individuals than the null model was fitted over is refused.
    pub(crate) fn test_the_block(
        &mut self,
        dosages: &GwasDosages,
        test: TestType,
    ) -> Result<Answers<'_>> {
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
        self.projected.resize(values, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: dosages.dosages(),
                rows: num_vars,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.projection,
                cols: self.num_individuals,
            },
            &mut self.projected,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of a block of variants with the projection matrix",
            source,
        })?;
        self.num.resize(num_vars, 0.0);
        popnei_linalg::product(
            TheFirstOperand::ByTheRowsOfTheResult {
                values: dosages.dosages(),
                rows: num_vars,
            },
            self.num_individuals,
            TheSecondOperand::ByTheValuesSummedOver {
                values: &self.projected_trait,
                cols: 1,
            },
            &mut self.num,
        )
        .map_err(|source| Error::GwasLinalg {
            operation: "product of a block of variants with the projected trait",
            source,
        })?;
        let degrees_of_freedom = self.degrees_of_freedom;
        let ypy = self.ypy;
        // The share of what the variant was that the projection has to
        // leave of it, and the share of `y' p y` that the Wald test's
        // subtraction has to leave, for the variant to be worth testing.
        let share_that_is_nothing = the_share_that_is_nothing(self.num_individuals);
        let largest_of_the_projection = self.largest_of_the_projection;
        for ((row, num), of_the_variant) in self
            .projected
            .chunks_exact(self.num_individuals.max(1))
            .zip(&self.num)
            .zip(dosages.dosages().chunks_exact(self.num_individuals.max(1)))
        {
            let den = row
                .iter()
                .zip(of_the_variant)
                .map(|(projected, dosage)| projected * dosage)
                .sum::<f64>();
            let of_the_dosages = of_the_variant
                .iter()
                .map(|dosage| dosage * dosage)
                .sum::<f64>();
            if den <= share_that_is_nothing * largest_of_the_projection * of_the_dosages {
                self.beta.push(f64::NAN);
                self.se.push(f64::NAN);
                self.p_value.push(f64::NAN);
                continue;
            }
            let beta = num / den;
            let statistic = num * num / den;
            let answered = match test {
                TestType::Wald => {
                    // What the variant leaves of `y' p y`. The projection
                    // annihilates the design, so any affine image of the
                    // trait gives `num² / den = y' p y` in exact
                    // arithmetic, and what this holds for such a variant is
                    // the rounding of that cancellation, of whichever sign
                    // it fell on and differing by half between the two
                    // backends. `se` would be the square root of a number
                    // divided by noise, or of a negative one, which is the
                    // NaN beside a finite `beta` that **Open 2** of
                    // `docs/specs/gwas.md` records.
                    let left = ypy - statistic;
                    if left <= share_that_is_nothing * ypy {
                        None
                    } else {
                        let se = (left / (degrees_of_freedom * den)).sqrt();
                        Some((se, t_sf_two_sided(beta / se, degrees_of_freedom)))
                    }
                }
                // The score test divides by `den` and forms no such
                // subtraction, so a variant that explains the whole of what
                // the null left is answered here.
                TestType::Score => Some((1.0 / den.sqrt(), chi2_sf_1df(statistic))),
            };
            let Some((se, p_value)) = answered else {
                self.beta.push(f64::NAN);
                self.se.push(f64::NAN);
                self.p_value.push(f64::NAN);
                continue;
            };
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

/// The two variances the covariance of a trait under a mixed model is
/// built from, `v = genetic * k + residual * i` with `k` the kinship and
/// `i` the identity.
///
/// They travel as one value and not as two arguments of the same primitive
/// side by side, which the `coding` skill forbids for the reason this pair
/// shows plainly: swapping them compiles, gives an inverse of another
/// matrix and every variant an answer that is wrong in a way no type
/// catches.
#[derive(Debug, Clone, Copy)]
struct TheVariances {
    /// The variance of the random effect the kinship is the covariance of.
    genetic: f64,
    /// What is left over, the variance of an individual's own noise.
    residual: f64,
}

/// The projection matrix of a fitted linear mixed model.
///
/// The covariance of the trait under the null is the kinship times the
/// genetic variance plus the identity times the residual one, `v =
/// genetic_variance * k + residual_variance * i`, and the
/// eigendecomposition of the kinship gives its inverse without another
/// factorization: with `e` the eigenvectors and `l` the eigenvalues, `v⁻¹
/// = e diag(1 / (genetic * l + residual)) e'`. The projection matrix is
/// then `p = v⁻¹ - v⁻¹ d (d' v⁻¹ d)⁻¹ d' v⁻¹`.
///
/// The trait through it, `p y`, and `y' p y` are taken at the call site and
/// not here, although the linear mixed model wants all three: what the
/// logistic mixed model of the next plan wants of this is the matrix alone,
/// its `p y` being the trait less the fitted mean and its inverse coming
/// off a Cholesky factorization, and it has no `y' p y` at all. pyNei keeps
/// the reusable piece alone in the same way, as `_projection`. Splitting it
/// while it has one caller costs four lines at that caller and saves the
/// next plan from splitting a function two models depend on.
///
/// `eigenvectors` is `num_individuals` x `num_individuals`, row after row,
/// row `j` being the eigenvector of `eigenvalues[j]`. Three matrices of
/// that size are held at once here, which at the 10000 individuals of
/// `docs/objectives.md` is 2.4 GB, and one of them is what comes back, for
/// the model to test its variants through; pyNei holds the same three and
/// keeps the same one.
///
/// # Errors
///
/// [`Error::GwasLinalg`] when one of the four products, the Cholesky
/// factorization of the design weighted by the covariance or the solve
/// against it could not be done.
fn the_projection_of(
    design: &Design<'_>,
    eigenvalues: &[f64],
    eigenvectors: &[f64],
    variances: TheVariances,
) -> Result<Vec<f64>> {
    let num_individuals = design.num_individuals();
    let num_coefs = design.num_coefs();
    let mut scaled = vec![0.0_f64; eigenvectors.len()];
    for ((into, eigenvector), eigenvalue) in scaled
        .chunks_exact_mut(num_individuals.max(1))
        .zip(eigenvectors.chunks_exact(num_individuals.max(1)))
        .zip(eigenvalues)
    {
        let of_the_covariance = variances.genetic * eigenvalue + variances.residual;
        for (value, of_the_eigenvector) in into.iter_mut().zip(eigenvector) {
            *value = of_the_eigenvector / of_the_covariance;
        }
    }
    let mut inverse = vec![0.0_f64; eigenvectors.len()];
    popnei_linalg::product(
        TheFirstOperand::ByTheValuesSummedOver {
            values: eigenvectors,
            rows: num_individuals,
        },
        num_individuals,
        TheSecondOperand::ByTheValuesSummedOver {
            values: &scaled,
            cols: num_individuals,
        },
        &mut inverse,
    )
    .map_err(|source| Error::GwasLinalg {
        operation: "inverse of the covariance of the trait",
        source,
    })?;
    let mut of_the_covariance = vec![0.0_f64; design.values().len()];
    popnei_linalg::product(
        TheFirstOperand::ByTheRowsOfTheResult {
            values: &inverse,
            rows: num_individuals,
        },
        num_individuals,
        TheSecondOperand::ByTheValuesSummedOver {
            values: design.values(),
            cols: num_coefs,
        },
        &mut of_the_covariance,
    )
    .map_err(|source| Error::GwasLinalg {
        operation: "product of the covariance's inverse with the design",
        source,
    })?;
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "the design is one slice of `num_individuals` x `num_coefs` values and a \
                  study has two more individuals than columns, so `num_coefs` times \
                  itself is smaller than a length a `usize` holds"
    )]
    let of_the_coefs = num_coefs * num_coefs;
    let mut dvd = vec![0.0_f64; of_the_coefs];
    popnei_linalg::product(
        TheFirstOperand::ByTheValuesSummedOver {
            values: design.values(),
            rows: num_coefs,
        },
        num_individuals,
        TheSecondOperand::ByTheValuesSummedOver {
            values: &of_the_covariance,
            cols: num_coefs,
        },
        &mut dvd,
    )
    .map_err(|source| Error::GwasLinalg {
        operation: "product of the design with the covariance's inverse",
        source,
    })?;
    popnei_linalg::cholesky_lower(&mut dvd, num_coefs).map_err(|source| Error::GwasLinalg {
        operation: "Cholesky factorization of the design weighted by the covariance",
        source,
    })?;
    let mut solved = of_the_covariance.clone();
    popnei_linalg::solve_with_cholesky(&dvd, num_coefs, &mut solved, num_individuals).map_err(
        |source| Error::GwasLinalg {
            operation: "solve of the design weighted by the covariance",
            source,
        },
    )?;
    // The buffer of the scaled eigenvectors is done with, and what goes in
    // it is the part of the covariance's inverse that the design explains,
    // which is taken out of that inverse to leave the projection matrix.
    let mut of_the_design = scaled;
    popnei_linalg::product(
        TheFirstOperand::ByTheRowsOfTheResult {
            values: &of_the_covariance,
            rows: num_individuals,
        },
        num_coefs,
        TheSecondOperand::ByTheColumnsOfTheResult {
            values: &solved,
            cols: num_individuals,
        },
        &mut of_the_design,
    )
    .map_err(|source| Error::GwasLinalg {
        operation: "product that takes the design out of the covariance's inverse",
        source,
    })?;
    let mut projection = inverse;
    for (entry, explained) in projection.iter_mut().zip(&of_the_design) {
        *entry -= *explained;
    }
    Ok(projection)
}

/// The trait through a projection matrix, `p y`, one value per tested
/// individual.
///
/// It is what every variant's numerator is taken against, `num` being the
/// variant times this.
///
/// # Errors
///
/// [`Error::GwasLinalg`] when the product could not be done.
fn through_the_projection(projection: &[f64], phenotype: &[f64]) -> Result<Vec<f64>> {
    let mut through = vec![0.0_f64; phenotype.len()];
    popnei_linalg::product(
        TheFirstOperand::ByTheRowsOfTheResult {
            values: projection,
            rows: phenotype.len(),
        },
        phenotype.len(),
        TheSecondOperand::ByTheValuesSummedOver {
            values: phenotype,
            cols: 1,
        },
        &mut through,
    )
    .map_err(|source| Error::GwasLinalg {
        operation: "product of the projection matrix with the trait",
        source,
    })?;
    Ok(through)
}

/// The null model of the linear mixed model against GMMAT 1.5.0's
/// `glmmkin`, and the fit's own identity, which "The linear mixed model"
/// of `docs/specs/gwas.md` has as the two checks of that fit.
#[cfg(test)]
pub(crate) mod lmm {
    use std::io::Read;
    use std::path::{Path, PathBuf};

    use super::{
        GOLDEN_SECTION_RATIO, GOLDEN_SECTION_STEPS, LOG_DELTA_POINTS, LinearMixedModel,
        the_log_delta_at,
    };
    use crate::block::BlockReader;
    use crate::error::Error;
    use crate::gwas::linear::lm::{
        THE_DESIGN_OF_TWO_SUBPOPULATIONS, THE_HEADER_OF_EIGHT,
        THE_INDIVIDUALS_OF_TWO_SUBPOPULATIONS, THE_TRAIT_OF_TWO_SUBPOPULATIONS, reader_over,
        the_study_of, the_trait_and_the_design_of_the_panel,
    };
    use crate::gwas::result::Gwas;
    use crate::gwas::study::{Design, GwasInput, GwasModel, TestType, TraitType};
    use crate::io::vcf::{VcfOptions, VcfReader};

    /// How far each of the five numbers of the null model may be from
    /// GMMAT's: 1e-5 absolute.
    ///
    /// It is the bound of "How it is verified" of "The linear mixed model"
    /// of `docs/specs/gwas.md`, and it is the distance between two
    /// restricted maximum likelihood searches and not the width of a
    /// printed digit: `tests/reference/gwas/gmmat.null_models.tsv` is one
    /// of the three references the script writes at full precision, 15
    /// digits of each number, so nothing of this bound is spent on
    /// rounding.
    ///
    /// Measured over the five and the heritability on 24 September 2026:
    /// the worst is the genetic variance, 1.230e-6 away on Accelerate and
    /// 1.218e-6 on faer, which is 12 per cent of what is allowed; the
    /// residual variance is 1.048e-6 and 1.041e-6 away, and the largest of
    /// the three effects is the 6.93e-7 of `cov2` on both. The bound is not
    /// lowered to two or three times that, unlike the ones this module sets
    /// on popnei's own arithmetic: what it measures is how far the searches
    /// of two programs land apart, and the spec fixes it.
    ///
    /// Where that 1.2e-6 comes from is the kinship and not the search. The
    /// reference script gives GMMAT the text that `plink2 --make-rel
    /// square` wrote, six significant digits of each entry, and this test
    /// reads the `f64` of `--make-rel square bin` beside it, which is the
    /// same matrix to 4.95e-6 of an entry. Fitted on that text in numpy on
    /// 24 September 2026 the genetic variance lands 3.77e-7 from GMMAT's
    /// instead of 1.2e-6.
    ///
    /// What the search itself is worth is 1.18e-8, which is how far the two
    /// backends put the genetic variance from each other on the same
    /// kinship: the criterion is flat at its minimum, so an eigenvalue that
    /// moves in its last bits moves the fitted ratio of the two variances
    /// by about the square root of that.
    const OF_GMMAT: f64 = 1e-5;

    /// How far `y' p y` may be from the individuals less the columns of
    /// the design: 1e-6 absolute, from "The linear mixed model" of
    /// `docs/specs/gwas.md`.
    ///
    /// Measured on the two panels on 24 September 2026, the worst is
    /// 9.4e-12 on Accelerate and 7.6e-12 on faer, both on the panel with
    /// every genotype called, which is 0.001 per cent of what is allowed.
    /// The bound is the spec's and is not lowered to two or three times
    /// what was measured: 197 is what the fit is at its optimum, and the
    /// arithmetic that reaches it is products of 200 x 200 matrices, whose
    /// rounding grows with the individuals where the number it is compared
    /// with does not.
    const OF_THE_REML_IDENTITY: f64 = 1e-6;

    /// What GMMAT 1.5.0's `glmmkin` fitted for the panel with every
    /// genotype called, from `tests/reference/gwas/gmmat.null_models.tsv`:
    /// the variance of the random effect of the kinship, what is left
    /// over, and the effects of the intercept, of `cov1` and of `cov2`.
    ///
    /// GMMAT calls the first two `tau` and `sigma2`. The file holds them
    /// at full precision and these are its digits.
    const OF_GMMAT_NULL: (f64, f64, [f64; 3]) = (
        1.22161667529699,
        0.342359482266917,
        [4.67802051309181, 0.473360959469751, 1.11027907093373],
    );

    /// The genetic variance of those two over their sum, which is the
    /// share of the trait's variance that the kinship explains.
    const OF_GMMAT_HERITABILITY: f64 = 0.7810967382008012;

    /// The path of one of the files of `tests/reference/kinship/`, where
    /// the two kinships plink2 wrote are.
    fn the_kinship_path(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/kinship")
            .join(name)
    }

    /// The 40000 entries of the kinship plink2 wrote for one of the two
    /// panels, row after row: the little endian `f64` of `--make-rel
    /// square bin`, which is the matrix at full precision.
    ///
    /// The mixed models are given the kinship that came from plink2 and
    /// not from popnei or pyNei, as "How it is verified" of "What every
    /// model shares" of `docs/specs/gwas.md` asks.
    pub(crate) fn the_kinship_of(name: &str) -> Vec<f64> {
        let path = the_kinship_path(&format!("{name}.plink2.rel.bin.gz"));
        let file = match std::fs::File::open(&path) {
            Ok(file) => file,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let mut bytes = Vec::new();
        if let Err(error) = flate2::read::GzDecoder::new(file).read_to_end(&mut bytes) {
            panic!("{path}: {error}", path = path.display());
        }
        bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|eight| f64::from_le_bytes(*eight))
            .collect()
    }

    /// The individuals of that kinship, in the order plink2 wrote its rows
    /// and its columns in, which is the order the VCF has them.
    ///
    /// The file holds one header line, `#IID`, and then one name per line.
    pub(crate) fn the_individuals_of_the_kinship(name: &str) -> Vec<String> {
        let path = the_kinship_path(&format!("{name}.plink2.rel.id"));
        let text = match std::fs::read_to_string(&path) {
            Ok(text) => text,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        text.lines()
            .skip(1)
            .filter(|line| !line.is_empty())
            .map(|line| line.trim().to_owned())
            .collect()
    }

    /// The linear mixed model of one of the two panels, fitted on the
    /// trait `cont` and the covariates `cov1` and `cov2` of
    /// `tests/reference/gwas/phenotypes.csv` over the kinship plink2 wrote
    /// for `name`.
    ///
    /// It is fitted here and not through `calc_gwas` because `y' p y` and
    /// the two variances as the fit holds them are in no result: the null
    /// model a user reads carries the two variances and the tests of this
    /// module that read them through `calc_gwas` are the ones against
    /// GMMAT and rrBLUP.
    fn the_null_of_the_panel(name: &str) -> LinearMixedModel {
        let individuals = the_individuals_of_the_kinship(name);
        assert_eq!(individuals.len(), 200, "the individuals of {name}");
        let kinship = the_kinship_of(name);
        assert_eq!(kinship.len(), 40000, "the entries of the kinship of {name}");
        let (phenotype, values) =
            the_trait_and_the_design_of_the_panel(&individuals, TraitType::Continuous);
        let tested: Vec<usize> = (0..individuals.len()).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Continuous,
            design: &values,
            num_coefs: 3,
            kinship: Some(&kinship),
            test: None,
            use_grammar_gamma_approx: false,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        let design = match Design::of_the_study(&study, individuals.len()) {
            Ok(design) => design,
            Err(error) => panic!("the design of {name}: {error}"),
        };
        match LinearMixedModel::of_the_study(&phenotype, &design, &kinship) {
            Ok(fitted) => fitted,
            Err(error) => panic!("the null model of {name}: {error}"),
        }
    }

    /// `found` within [`OF_GMMAT`] of `expected`, absolute.
    fn assert_of_gmmat(found: f64, expected: f64, what: &str) {
        let difference = (found - expected).abs();
        assert!(
            difference <= OF_GMMAT,
            "{what} is {found} and GMMAT gives {expected}, {difference} away, against the \
             {OF_GMMAT} allowed"
        );
    }

    /// The two variances and the three effects of the null model are
    /// GMMAT's, within 1e-5 absolute, which is deliverable 1 of work
    /// package 4 of `docs/plans/gwas-linear.md`.
    ///
    /// The search has to be reproduced step for step or the variances
    /// move, so a `genetic_variance` near but not at GMMAT's means the
    /// search and one far from it means the criterion or the clamp of the
    /// eigenvalues at 0.
    ///
    /// The `heritability` is asserted with them because it is the only
    /// number of `NullModel` that is built from the two variances rather
    /// than read off the fit.
    ///
    /// That all three are given at all is the other end of the check of
    /// **Open 3** of `docs/specs/gwas.md`: this kinship tells the two
    /// variances apart, so the fit reports them, where the identity of
    /// [`lmm::a_kinship_that_does_not_tell_the_variances_apart_gives_none_of_them`]
    /// does not and it reports none of them.
    #[test]
    fn the_null_of_the_panel_is_gmmats_two_variances_and_three_effects() {
        let fitted = the_null_of_the_panel("panel_called");
        let null = fitted.null_model(TestType::Score);
        assert_eq!(null.model, GwasModel::Lmm, "the model that was fitted");
        assert_eq!(null.test, TestType::Score, "the test it carries");
        assert_eq!(
            null.num_individuals, 200,
            "the individuals it was fitted on"
        );
        let (genetic, residual, effects) = OF_GMMAT_NULL;
        let found = null
            .genetic_variance
            .expect("a linear mixed model has a genetic variance");
        assert_of_gmmat(found, genetic, "the genetic variance");
        let found = null
            .residual_variance
            .expect("a linear mixed model has a residual variance");
        assert_of_gmmat(found, residual, "the residual variance");
        let found = null
            .heritability
            .expect("a linear mixed model has a heritability");
        assert_of_gmmat(found, OF_GMMAT_HERITABILITY, "the heritability");
        assert_eq!(
            null.covariate_effects.len(),
            3,
            "one effect per column of the design"
        );
        for ((found, expected), what) in null.covariate_effects.iter().zip(effects).zip([
            "the intercept",
            "the effect of cov1",
            "the effect of cov2",
        ]) {
            assert_of_gmmat(*found, expected, what);
        }
    }

    /// `y' p y` is the individuals less the columns of the design, 197 on
    /// both panels within 1e-6, which is `test_reml_identity` of pyNei.
    ///
    /// **What it does not check is that the search reached its optimum**,
    /// which this comment and deliverable 2 of work package 4 of
    /// `docs/plans/gwas-linear.md` both said until 25 September 2026. The
    /// genetic variance is that same generalized residual sum of squares
    /// divided by the same degrees of freedom, so `y' p y` comes to `n - c`
    /// for any `delta` whatever: it is an algebraic identity and not
    /// evidence about where the search landed. Measured that day by
    /// multiplying the fitted `delta` by a million, it came out
    /// 196.99999999999872, 1.3e-12 away, while every comparison with GMMAT
    /// and rrBLUP went red.
    ///
    /// What it does check is worth keeping, and it is two things. That the
    /// clamp of the negative eigenvalues works, which is why the second
    /// panel is here. And that the projection matrix and the quadratic form
    /// agree with each other: `p` is built from the inverse of the
    /// covariance and the design, and `y' p y` is the trait through it, so
    /// the identity fails if either is formed wrong. Where the search
    /// landed is
    /// [`lmm::the_delta_the_search_lands_on_for_the_panel_is_pyneis`].
    ///
    /// It is the one check of this spec made at the private function that
    /// fits the null, since `y' p y` is in no result.
    ///
    /// The second panel is here for the clamp of the eigenvalues at 0. Its
    /// kinship, the one plink2 wrote for the panel with 3 in 100 genotypes
    /// missing whole, has a smallest eigenvalue of -0.0321 where the
    /// panel with every genotype called has -3.44e-15. Without the clamp
    /// the weights of the criterion go negative at every `delta` below
    /// 0.0321, which is 33 of the 101 points of the grid, and the Cholesky
    /// factorization of the weighted design refuses them, so the fit comes
    /// back with an error and this test fails.
    #[test]
    fn the_trait_through_the_projection_of_a_panel_is_its_individuals_less_its_design() {
        for name in ["panel_called", "panel"] {
            let fitted = the_null_of_the_panel(name);
            let found = fitted.ypy();
            let difference = (found - 197.0).abs();
            assert!(
                difference <= OF_THE_REML_IDENTITY,
                "y' p y of {name} is {found} and the fit at its optimum makes it 197, \
                 {difference} away, against the {OF_THE_REML_IDENTITY} allowed"
            );
        }
    }
    /// The VCF of the panel `name`: the one with every genotype called sits
    /// beside the kinships and the one with 3 in 100 genotypes missing
    /// whole beside the distances, and both hold the same 200 individuals
    /// and the same 1200 variants.
    fn the_vcf_of_the_panel(name: &str) -> PathBuf {
        let of_the_module = match name {
            "panel_called" => "kinship",
            _ => "dists",
        };
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference")
            .join(of_the_module)
            .join(format!("{name}.vcf.gz"))
    }

    /// The study of the variants of one panel against its continuous
    /// trait, with the kinship plink2 wrote for the panel with every
    /// genotype called and the test `test`.
    ///
    /// `with_cov1` says whether the continuous covariate goes into the
    /// design beside the intercept and the binary one: rrBLUP takes every
    /// fixed effect as a factor, so it was given `cov2` alone and popnei is
    /// run with the same one covariate for that comparison, while GMMAT was
    /// given both.
    ///
    /// The kinship is `panel_called`'s for both panels, because that is
    /// what `tests/reference/gwas/make_reference.py` gave GMMAT: it fits
    /// one null model with it and then scores the variants of each panel
    /// against that fit.
    fn the_study_of_the_panel(name: &str, test: TestType, with_cov1: bool) -> Gwas {
        let path = the_vcf_of_the_panel(name);
        let options = VcfOptions {
            ploidy: 2,
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        let individuals = reader.individuals().to_vec();
        assert_eq!(
            individuals,
            the_individuals_of_the_kinship("panel_called"),
            "the individuals of {name} against the ones of the kinship"
        );
        let kinship = the_kinship_of("panel_called");
        let (phenotype, both) =
            the_trait_and_the_design_of_the_panel(&individuals, TraitType::Continuous);
        let design: Vec<f64> = match with_cov1 {
            true => both,
            false => both
                .as_chunks::<3>()
                .0
                .iter()
                .flat_map(|row| [row[0], row[2]])
                .collect(),
        };
        let num_coefs = match with_cov1 {
            true => 3,
            false => 2,
        };
        let tested: Vec<usize> = (0..individuals.len()).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Continuous,
            design: &design,
            num_coefs,
            kinship: Some(&kinship),
            test: Some(test),
            use_grammar_gamma_approx: false,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study of {name}: {error}"),
        }
    }

    /// Where the variant `id` is among the rows of a study of a panel.
    fn the_row_of(result: &Gwas, id: &str) -> usize {
        let ids = result.ids.as_deref().expect("the ids of the variants");
        match ids.iter().position(|held| held == id) {
            Some(var) => var,
            None => panic!("{id} is not a variant of the panel"),
        }
    }

    /// How far a `-log10(p_value)` of the Wald test may be from rrBLUP's:
    /// 1e-4 absolute.
    ///
    /// It is the bound of "How it is verified" of "The linear mixed model"
    /// of `docs/specs/gwas.md`, which asks the same of the whole column in
    /// Python. `tests/reference/gwas/rrblup.panel_called.lmm.tsv` is one of
    /// the three references the script writes at full precision, 15 digits
    /// of each number, so nothing of this bound is spent on rounding: it is
    /// the distance between two fits, and 1e-4 in `-log10(p)` is 2.3e-4 of
    /// the p-value itself.
    ///
    /// Measured over the six variants on 24 September 2026: the worst is
    /// `var0052`, 8.874e-6 away on Accelerate and 8.864e-6 on faer, which
    /// is 9 per cent of what is allowed. Over all 1200 variants of the
    /// panel the worst is 1.997e-5 and 1.996e-5, 20 per cent of it, which
    /// is what the comparison over the whole column in Python has to hold.
    /// The bound is the spec's and is not lowered to two or three times
    /// that, as the bounds this module sets on popnei's own arithmetic are:
    /// what it measures is how far popnei's restricted maximum likelihood
    /// search and rrBLUP's land apart, and the spec fixes it.
    const OF_RRBLUP: f64 = 1e-4;

    /// What rrBLUP 4.6.3's `GWAS` with `P3D = TRUE` gave for six variants
    /// of the panel with every genotype called, in
    /// `tests/reference/gwas/rrblup.panel_called.lmm.tsv`: the id and
    /// `-log10(p_value)`, which is all it reports.
    ///
    /// It was given the kinship and `cov2` as its one fixed effect, so the
    /// study here is run with the same one covariate. Five of the six are
    /// the causal variants of `causal_vars.csv` and `var0000` is not
    /// causal.
    const OF_RRBLUP_SIX: [(&str, f64); 6] = [
        ("var0000", 0.215_210_110_260_573),
        ("var0052", 3.618_698_839_162_37),
        ("var0629", 4.339_896_431_706_14),
        ("var0751", 2.065_324_779_038_01),
        ("var1137", 1.319_283_317_997_48),
        ("var1188", 2.367_278_649_793_38),
    ];

    /// The six variants of the panel are rrBLUP's under the Wald test,
    /// which is deliverable 3 of work package 4 of
    /// `docs/plans/gwas-linear.md`.
    ///
    /// The Wald test divides the effect of the variant by its standard
    /// error and asks how far a Student t with 197 degrees of freedom goes
    /// beyond that. What is compared is `-log10(p_value)`, because it is
    /// all rrBLUP reports, and the run is with `cov2` alone: rrBLUP takes
    /// every fixed effect as a factor, so it was given that covariate and
    /// no other, and a run that gave it both would get numbers that are
    /// close and not equal.
    #[test]
    fn the_six_variants_of_the_panel_are_rrblups_p_value_of_the_wald_test() {
        let result = the_study_of_the_panel("panel_called", TestType::Wald, false);
        assert_eq!(result.num_vars, 1200, "the variants of the panel");
        assert_eq!(
            result.null_model.test,
            TestType::Wald,
            "the test the result says it made"
        );
        for (id, expected) in OF_RRBLUP_SIX {
            let var = the_row_of(&result, id);
            let found = -result.p_value[var].log10();
            let difference = (found - expected).abs();
            assert!(
                difference <= OF_RRBLUP,
                "-log10(p) of {id} is {found} and rrBLUP gives {expected}, {difference} \
                 away, against the {OF_RRBLUP} allowed"
            );
        }
    }

    /// How far `1 / se²` of the score test may be from GMMAT's `VAR`, as a
    /// share of it: 1e-5.
    ///
    /// `VAR` is the variance of the score, which is `x' p x`, the
    /// denominator both tests are built on. The spec asks the same of the
    /// whole column in Python.
    /// `tests/reference/gwas/gmmat.panel_called.lmm.score.tsv` and the file
    /// of the other panel are printed to six significant digits, which
    /// rounds a value by up to 5e-6 of itself, so half of this bound can go
    /// on GMMAT's printing alone and the comparison has twofold headroom at
    /// best.
    ///
    /// Measured over the six variants of both panels on 24 September 2026:
    /// the worst is `var0052` of the panel with every genotype called,
    /// 1.842e-6 of `VAR` on Accelerate and 1.843e-6 on faer, which is 18
    /// per cent of what is allowed and is smaller than the 5e-6 GMMAT's own
    /// printing can account for. So the denominator both tests are built on
    /// is GMMAT's to every digit GMMAT printed, and nothing of popnei's
    /// arithmetic can be seen here.
    const OF_GMMAT_VARIANCE: f64 = 1e-5;

    /// How far a p-value of the score test may be from GMMAT's, in
    /// `log10`: 1e-4.
    ///
    /// The p-values of a study span 23 orders of magnitude and what a user
    /// reads is the exponent, so they are compared in `log10`, as the spec
    /// asks of the whole column in Python. `log10` shrinks a relative
    /// difference, so this bound has more headroom over GMMAT's six printed
    /// digits than [`OF_GMMAT_VARIANCE`] has.
    ///
    /// Measured over the six variants of both panels on 24 September 2026:
    /// the worst is `var0629` of the panel with every genotype called,
    /// 4.164e-5 on both backends, which is 42 per cent of what is allowed.
    /// That one is not the printing: 4.164e-5 in `log10` is 9.6e-5 of the
    /// p-value, against the 5e-6 six digits round it by. It is the two
    /// fits, which land 1.2e-6 apart in the genetic variance, at a p-value
    /// of 4.8e-5, where the tail of the chi square turns a small move of
    /// the statistic into a larger one of the p-value.
    const OF_GMMAT_P_VALUE: f64 = 1e-4;

    /// What GMMAT 1.5.0's `glmm.score` gave for six variants of each
    /// panel, from `tests/reference/gwas/gmmat.panel_called.lmm.score.tsv`
    /// and `gmmat.panel.lmm.score.tsv`: the id, then `VAR` and `PVAL` with
    /// every genotype called and the two of the panel with 3 in 100
    /// genotypes missing whole.
    ///
    /// `VAR` is the variance of the score, `x' p x`, which is `1 / se²`.
    /// GMMAT was given both covariates and, for both panels, the kinship of
    /// the panel with every genotype called. A missing genotype takes the
    /// mean dosage of its variant, which GMMAT calls `impute2mean` and
    /// which is popnei's rule too, so the two agree on the second panel.
    const OF_GMMAT_SIX: [(&str, f64, f64, f64, f64); 6] = [
        ("var0000", 29.8774, 0.360_526, 29.9961, 0.495_719),
        ("var0052", 43.8076, 0.001_189_85, 46.0763, 0.001_058_95),
        ("var0629", 31.7241, 4.810_05e-5, 31.9307, 7.378_88e-5),
        ("var0751", 44.5825, 0.004_392_26, 44.422, 0.006_756_44),
        ("var1137", 43.3724, 0.013_926, 42.1841, 0.012_532_1),
        ("var1188", 47.3766, 0.001_073_4, 47.9106, 0.001_505_77),
    ];

    /// The six variants of both panels are GMMAT's under the score test,
    /// which is deliverable 4 of work package 4 of
    /// `docs/plans/gwas-linear.md`.
    ///
    /// The score test asks how steeply the likelihood rises at an effect of
    /// 0, so it holds both variances at the null and needs no fit with the
    /// variant in. What is compared is `1 / se²` against the variance of
    /// the score GMMAT reports, which is the denominator of both tests, and
    /// the p-value in `log10`.
    ///
    /// The second panel is what says that a missing genotype takes the mean
    /// dosage of its variant as GMMAT's `impute2mean` does: 3 in 100 of its
    /// genotypes are missing whole, and every one of the six numbers moves
    /// between the panels.
    #[test]
    fn the_six_variants_of_both_panels_are_gmmats_variance_of_the_score_and_p_value() {
        for (name, of_the_panel) in [("panel_called", 0_usize), ("panel", 1)] {
            let result = the_study_of_the_panel(name, TestType::Score, true);
            assert_eq!(result.num_vars, 1200, "the variants of {name}");
            assert_eq!(
                result.null_model.test,
                TestType::Score,
                "the test the result of {name} says it made"
            );
            for (id, called_variance, called_p, missing_variance, missing_p) in OF_GMMAT_SIX {
                let (variance, p_value) = match of_the_panel {
                    0 => (called_variance, called_p),
                    _ => (missing_variance, missing_p),
                };
                let var = the_row_of(&result, id);
                let se = result.se[var];
                let found = 1.0 / (se * se);
                let difference = (found - variance).abs();
                assert!(
                    difference <= OF_GMMAT_VARIANCE * variance,
                    "1 / se² of {id} of {name} is {found} and GMMAT gives {variance}, \
                     {difference} away, which is {share} of it against the \
                     {OF_GMMAT_VARIANCE} allowed",
                    share = difference / variance
                );
                let found = result.p_value[var];
                let difference = (found / p_value).log10().abs();
                assert!(
                    difference <= OF_GMMAT_P_VALUE,
                    "the p-value of {id} of {name} is {found} and GMMAT gives {p_value}, \
                     {difference} apart in log10 against the {OF_GMMAT_P_VALUE} allowed"
                );
            }
        }
    }

    /// The Wald test and the score test give the same effect for every
    /// variant, to the bit, and a different standard error for every one of
    /// them.
    ///
    /// `num` is the variant times the trait through the projection matrix
    /// and `den` the variant times itself through it, and the effect is
    /// `num / den` in both tests. What they differ in is how uncertain that
    /// effect is and which distribution it is read against, so a study that
    /// answered two different effects would have its defect in the
    /// projection matrix and not in either test.
    #[test]
    fn the_two_tests_of_a_variant_differ_in_the_error_and_not_in_the_effect() {
        let of_the_wald = the_study_of_the_panel("panel_called", TestType::Wald, true);
        let of_the_score = the_study_of_the_panel("panel_called", TestType::Score, true);
        assert_eq!(of_the_wald.num_vars, 1200, "the variants of the panel");
        for (var, (wald, score)) in of_the_wald.beta.iter().zip(&of_the_score.beta).enumerate() {
            assert_eq!(
                wald.total_cmp(score),
                std::cmp::Ordering::Equal,
                "the effect of the variant {var} is {wald} under the Wald test and {score} \
                 under the score test"
            );
        }
        let differ = of_the_wald
            .se
            .iter()
            .zip(&of_the_score.se)
            .filter(|(wald, score)| wald.total_cmp(score) != std::cmp::Ordering::Equal)
            .count();
        assert_eq!(
            differ, 1200,
            "how many of the 1200 variants the two tests gave a different standard error"
        );
    }

    /// The GRAMMAR-Gamma approximation is refused for a study that has a
    /// kinship, which is the pair it is for, because popnei has not written
    /// it yet.
    ///
    /// A study with no kinship is refused for asking for it at all, which
    /// the test of the linear model covers. The two refusals say different
    /// things, and neither study runs: one that made the exact test of
    /// every variant and reported that it had approximated nothing would
    /// give the user no way to tell that what they asked for did not
    /// happen.
    #[test]
    fn the_grammar_gamma_approximation_with_a_kinship_is_refused() {
        let individuals = the_individuals_of_the_kinship("panel_called");
        let kinship = the_kinship_of("panel_called");
        let (phenotype, design) =
            the_trait_and_the_design_of_the_panel(&individuals, TraitType::Continuous);
        let tested: Vec<usize> = (0..individuals.len()).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Continuous,
            design: &design,
            num_coefs: 3,
            kinship: Some(&kinship),
            test: None,
            use_grammar_gamma_approx: true,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        let path = the_vcf_of_the_panel("panel_called");
        let options = VcfOptions {
            ploidy: 2,
            ..VcfOptions::default()
        };
        let mut reader = match VcfReader::from_path(&path, options) {
            Ok(reader) => reader,
            Err(error) => panic!("{path}: {error}", path = path.display()),
        };
        match the_study_of(&mut reader, &study) {
            Err(Error::GwasGrammarGammaNotBuilt) => {}
            Err(error) => panic!("the study was refused with {error}"),
            Ok(_) => panic!("the approximation was made"),
        }
    }
    /// The kinship of the eight individuals of **Open 2** of
    /// `docs/specs/gwas.md`, row after row: 1 for an individual with
    /// itself, 0.2 for two of the same subpopulation and 0 across the two.
    ///
    /// Its eigenvalues are 1.6 and 0.8, both above 0, so it is a
    /// covariance and the fit of the mixed model over it is an ordinary
    /// one. The structure is the covariate's, which is the point of the
    /// fixture: the variant that is refused is a combination of the
    /// columns of the design.
    const THE_KINSHIP_OF_TWO_SUBPOPULATIONS: [f64; 64] = [
        1.0, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, //
        0.2, 1.0, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, //
        0.2, 0.2, 1.0, 0.2, 0.0, 0.0, 0.0, 0.0, //
        0.2, 0.2, 0.2, 1.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 0.0, 0.0, 0.0, 1.0, 0.2, 0.2, 0.2, //
        0.0, 0.0, 0.0, 0.0, 0.2, 1.0, 0.2, 0.2, //
        0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 1.0, 0.2, //
        0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 1.0,
    ];

    /// The header of the six individuals of the fixture of the Wald test's
    /// own cancellation, below.
    const THE_HEADER_OF_SIX: &str = "##fileformat=VCFv4.2\n\
        ##contig=<ID=1>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti0\ti1\ti2\ti3\ti4\ti5\n";

    /// The design of that fixture: six individuals and one covariate that
    /// is 0 and 1, beside the intercept.
    const THE_DESIGN_OF_SIX: [f64; 12] = [
        1.0, 0.0, //
        1.0, 1.0, //
        1.0, 0.0, //
        1.0, 1.0, //
        1.0, 0.0, //
        1.0, 1.0,
    ];

    /// Its trait, `2 + 3 * cov + 1 * dosage` with the dosages of the first
    /// variant below, 0, 1, 2, 0, 1, 2: a trait the design and that variant
    /// together explain exactly, and nothing else about it matters.
    const THE_TRAIT_OF_SIX: [f64; 6] = [2.0, 6.0, 4.0, 5.0, 3.0, 7.0];

    /// The identity, which is the kinship of six individuals with no recent
    /// ancestor in common and is what a user passes to mean no relatedness.
    const THE_KINSHIP_OF_SIX: [f64; 36] = [
        1.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, //
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, //
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, //
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ];

    /// The positions of those six among the individuals the reader gives.
    const THE_INDIVIDUALS_OF_SIX: [usize; 6] = [0, 1, 2, 3, 4, 5];

    /// The trait of the worked example of `docs/specs/gwas.md`, which the
    /// design and the one variant of the fixture below leave something of:
    /// the fixture of the flat criterion needs a trait the model can fit,
    /// where [`THE_TRAIT_OF_SIX`] is one it explains exactly.
    const THE_WORKED_TRAIT_OF_SIX: [f64; 6] = [2.0, 3.0, 5.0, 4.0, 4.0, 7.0];

    /// A kinship that does not tell the two variances apart gives the study
    /// with the three fields of those variances empty, which is the
    /// meanwhile of **Open 3** of `docs/specs/gwas.md`.
    ///
    /// The identity is such a kinship, and it is what a user passes to mean
    /// no relatedness: with it the model is the ordinary linear one
    /// whatever the split between the genetic variance and the residual
    /// one, so the restricted maximum likelihood has nothing to choose
    /// between them and its criterion is flat over the whole grid. What the
    /// study gave before this, measured on 25 September 2026 on this
    /// fixture: a `genetic_variance` of 0.00022769126788113066 and a
    /// `heritability` of 6.830738036433921e-05, which reads as a small
    /// number and is an arbitrary one. Perturbing such a kinship by 1e-15
    /// gave heritabilities of 6.5e-5, 7.1e-5 and 0.967 over three seeds.
    ///
    /// The study is still given and every variant still answered: `beta`
    /// and `p_value` do not depend on the split, because the test is scale
    /// free and the model is the linear one here, and they are what the
    /// user mostly came for. The variant is the `v0` of the worked example,
    /// whose effect the spec gives as 1.5, so what this asserts of the
    /// answers is the spec's own number and not popnei's.
    ///
    /// The three fields are `None` and not 0 and not NaN: `None` is what
    /// the linear model gives for the two it has not, so a user meets one
    /// way of saying that a number is not there and not three.
    #[test]
    fn a_kinship_that_does_not_tell_the_variances_apart_gives_none_of_them() {
        let mut vcf = String::from(THE_HEADER_OF_SIX);
        vcf.push_str("1\t1000\tv0\tA\tT\t.\t.\t.\tGT");
        for genotype in ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1"] {
            vcf.push('\t');
            vcf.push_str(genotype);
        }
        vcf.push('\n');
        let study = GwasInput {
            phenotype: &THE_WORKED_TRAIT_OF_SIX,
            trait_type: TraitType::Continuous,
            design: &THE_DESIGN_OF_SIX,
            num_coefs: 2,
            kinship: Some(&THE_KINSHIP_OF_SIX),
            test: Some(TestType::Wald),
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS_OF_SIX,
            transform_to_biallelic: false,
        };
        let mut reader = reader_over(vcf.as_bytes());
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the study over an identity kinship: {error}"),
        };

        let null = &result.null_model;
        assert_eq!(null.model, GwasModel::Lmm, "the model that was fitted");
        assert_eq!(
            null.genetic_variance, None,
            "the genetic variance of a fit that cannot tell the two apart"
        );
        assert_eq!(
            null.residual_variance, None,
            "the residual variance of such a fit"
        );
        assert_eq!(null.heritability, None, "the heritability of such a fit");
        assert_eq!(
            null.covariate_effects.len(),
            2,
            "the effects are the fit's own and are given"
        );
        assert_eq!(result.num_vars, 1, "the variants of the fixture");
        let beta = result.beta[0];
        assert!(
            (beta - 1.5).abs() <= 1e-12,
            "the effect of v0 is {beta} and the worked example gives 1.5, which the              split between the two variances cannot move"
        );
        let p_value = result.p_value[0];
        assert!(
            (0.0..=1.0).contains(&p_value) && p_value.is_finite(),
            "the p-value of v0 is {p_value}"
        );
    }

    /// A variant that leaves nothing of the trait has no answer under the
    /// Wald test, which is the third place the meanwhile of **Open 2** of
    /// `docs/specs/gwas.md` refuses.
    ///
    /// This one is not the variant the design explains, which the test
    /// below covers: it is a variant the design leaves whole and that
    /// explains the whole of what the null model left. `se` is built from
    /// `y' p y` minus `num² / den`, and the projection annihilates the
    /// design, so any affine image of the trait gives `num² / den = y' p y`
    /// in exact arithmetic and what is left is rounding. The trait here is
    /// `2 + 3 * cov + 1 * dosage`, which is exactly that.
    ///
    /// What the study gave before the threshold was there, measured on 25
    /// September 2026: `beta` 1.00000 with `se` and `p_value` NaN under the
    /// Wald test, where the score test on the same fixture answers. A
    /// finite `beta` beside a NaN `se` is a fourth kind of NaN that the
    /// spec does not describe, so a user filtering on a missing effect
    /// keeps the row and reads 1.0 as an effect that was measured. That is
    /// what took "leave it alone" out of the options of that open point.
    ///
    /// The score test is asserted to answer the same variant, because it
    /// divides by `den` and never forms the subtraction: the two tests
    /// differ here, and a refusal written in the wrong place would take the
    /// score test's answer away with it.
    ///
    /// What this fixture does not check is the size of the threshold, and
    /// no fixture can check it on both backends. What is left is 0 in exact
    /// arithmetic, so its sign is whatever the rounding chose, and there is
    /// no regime between the rounding and a threshold that is the rounding
    /// scale. Measured with the threshold set to 0 on 25 September 2026:
    /// faer answers this variant with an `se` of 1.2167e-8, so faer is
    /// where the size is guarded, while Accelerate leaves a value at or
    /// below 0 and the sign alone refuses it. The evidence that the size is
    /// right is the panel, where what is left comes out 8.53e-13 on
    /// Accelerate and 3.98e-13 on faer against a threshold of 8.67e-12, and
    /// that measurement is in **Open 2** of the spec.
    #[test]
    fn a_variant_that_leaves_nothing_of_the_trait_has_no_wald_answer() {
        let mut vcf = String::from(THE_HEADER_OF_SIX);
        for (var, genotypes) in [
            ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1"],
            ["0/0", "0/1", "1/1", "1/1", "0/1", "0/0"],
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
        let study_with = |test| GwasInput {
            phenotype: &THE_TRAIT_OF_SIX,
            trait_type: TraitType::Continuous,
            design: &THE_DESIGN_OF_SIX,
            num_coefs: 2,
            kinship: Some(&THE_KINSHIP_OF_SIX),
            test: Some(test),
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS_OF_SIX,
            transform_to_biallelic: false,
        };
        let mut reader = reader_over(vcf.as_bytes());
        let of_the_wald = match the_study_of(&mut reader, &study_with(TestType::Wald)) {
            Ok(result) => result,
            Err(error) => panic!("the study of six under the Wald test: {error}"),
        };
        let mut reader = reader_over(vcf.as_bytes());
        let of_the_score = match the_study_of(&mut reader, &study_with(TestType::Score)) {
            Ok(result) => result,
            Err(error) => panic!("the study of six under the score test: {error}"),
        };

        assert_eq!(of_the_wald.num_vars, 2, "the variants of the fixture");
        assert!(
            of_the_wald.beta[0].is_nan()
                && of_the_wald.se[0].is_nan()
                && of_the_wald.p_value[0].is_nan(),
            "the variant that leaves nothing of the trait was answered under the Wald \
             test with a beta of {beta}, an se of {se} and a p-value of {p_value}",
            beta = of_the_wald.beta[0],
            se = of_the_wald.se[0],
            p_value = of_the_wald.p_value[0]
        );
        // The variant is in the result with its frequency, as every variant
        // that has no answer is.
        let found = of_the_wald.allele_freq[0];
        assert!(
            (found - 0.5).abs() <= 1e-12,
            "the frequency of v0 is {found} and its dosages are half the alleles"
        );
        assert!(
            of_the_wald.beta[1].is_finite()
                && of_the_wald.se[1] > 0.0
                && (0.0..=1.0).contains(&of_the_wald.p_value[1]),
            "the ordinary variant beside it was answered under the Wald test with a beta \
             of {beta}, an se of {se} and a p-value of {p_value}",
            beta = of_the_wald.beta[1],
            se = of_the_wald.se[1],
            p_value = of_the_wald.p_value[1]
        );
        assert!(
            of_the_score.beta[0].is_finite()
                && of_the_score.se[0] > 0.0
                && (0.0..=1.0).contains(&of_the_score.p_value[0]),
            "the score test of the same variant answered a beta of {beta}, an se of {se} \
             and a p-value of {p_value}, where it divides by `den` and forms no \
             subtraction",
            beta = of_the_score.beta[0],
            se = of_the_score.se[0],
            p_value = of_the_score.p_value[0]
        );
    }

    /// The `delta` that pyNei's `_reml_delta` fitted for `panel_called`,
    /// the residual variance over the genetic one, taken on 25 September
    /// 2026 with the kinship plink2 wrote, the trait `cont` and the
    /// covariates `cov1` and `cov2`, with the eigenvalues clamped at 0 as
    /// `_LMMNull` clamps them.
    const OF_PYNEI_DELTA: f64 = 0.280_252_265_667_530_1;

    /// How far popnei's `delta` may be from that: 2.5e-7 relative.
    ///
    /// It is not a bound on popnei's arithmetic but on how far two
    /// implementations of this search land apart, and it is wide because
    /// the criterion is flat at its minimum: an eigenvalue moving in its
    /// last bits moves `delta` by about the square root of that. Measured
    /// on 25 September 2026, popnei sits 7.72e-8 from pyNei on Accelerate
    /// and 6.94e-8 on faer, and the two backends sit 7.8e-9 from each
    /// other. This is 3.2 times the worse of the two.
    ///
    /// What it catches is a search that landed somewhere else: multiplying
    /// the fitted `delta` by a million, which the `y' p y` identity below
    /// does not notice, fails here. What it does not catch is a wrong value
    /// of one of the six constants the spec asks to be reproduced, and that
    /// is why it is not the only guard of them: shifting the grid's lowest
    /// point to -10.1 moves `delta` to 1.635e-7 from pyNei and a ratio of
    /// 0.61 in place of the golden one to 9.22e-8, both measured the same
    /// day, and a bound tight enough to see either would be 1.2 times the
    /// distance popnei legitimately sits at. The grid and the ratio are
    /// pinned to the bit by the two tests above instead.
    const OF_PYNEI_THE_SEARCH: f64 = 2.5e-7;

    /// The width the bracket of the golden section search is left with:
    /// the 0.4 the grid's two neighbours span, shrunk by
    /// [`GOLDEN_SECTION_RATIO`] sixty times.
    ///
    /// It is `0.4 * ((math.sqrt(5) - 1) / 2) ** 60` in Python, which gives
    /// 1.1555841497399695e-13, and it is what says that the ratio and the
    /// number of steps are both pyNei's: a ratio of 0.61 leaves 3.9e-13 and
    /// fifty steps leave 1.2e-9.
    const THE_BRACKET_AFTER_THE_SEARCH: f64 = 1.155_584_149_739_969_5e-13;

    /// The grid the search starts on is `numpy.linspace(-10, 10, 101)`, to
    /// the bit.
    ///
    /// Three of the six things "The linear mixed model" of
    /// `docs/specs/gwas.md` asks to be reproduced are in this one
    /// assertion: how many points there are, where the lowest is and the
    /// step between them. None of the three is guarded by any comparison
    /// with a reference program, because the criterion is flat enough at
    /// its minimum that moving the grid moves `delta` by less than the
    /// distance popnei and pyNei legitimately sit apart, which the comment
    /// on [`OF_PYNEI_THE_SEARCH`] measures.
    ///
    /// The literals are numpy's own, printed on 25 September 2026, and they
    /// are asserted with `total_cmp` and not within a tolerance: the
    /// points come out of `start + at * step` in both libraries, so they
    /// agree to the bit or the grid is another grid. The 45th of them is
    /// -1.1999999999999993 and not -1.2, which is what says the step is
    /// multiplied before the start is added.
    #[test]
    fn the_grid_of_the_search_is_the_101_points_numpy_gives() {
        assert_eq!(
            LOG_DELTA_POINTS, 101,
            "the points of numpy.linspace(-10, 10, 101)"
        );
        for (at, expected) in [
            (0_usize, -10.0_f64),
            (1, -9.8),
            (44, -1.199_999_999_999_999_3),
            (45, -1.0),
            (50, 0.0),
            (99, 9.8),
            (100, 10.0),
        ] {
            let found = the_log_delta_at(at);
            assert_eq!(
                found.total_cmp(&expected),
                std::cmp::Ordering::Equal,
                "the point {at} of the grid is {found} and numpy.linspace gives {expected}"
            );
        }
    }

    /// The ratio of the golden section is pyNei's, and sixty steps of it
    /// leave the bracket where pyNei leaves it.
    ///
    /// The constant is written out rather than computed, so this asserts it
    /// against the expression `_reml_delta` of `pynei/gwas.py` evaluates,
    /// `(math.sqrt(5) - 1) / 2`, which is correctly rounded in both
    /// languages and so is the same `f64`. The width the bracket is left
    /// with is the other half: it is the ratio and the number of steps
    /// together, and it is what a flat 0.61 in place of the ratio changes,
    /// where the comparison with pyNei cannot see it.
    #[test]
    fn the_golden_section_of_the_search_is_pyneis_ratio_and_sixty_steps() {
        let of_pynei = (5.0_f64.sqrt() - 1.0) / 2.0;
        assert_eq!(
            GOLDEN_SECTION_RATIO.total_cmp(&of_pynei),
            std::cmp::Ordering::Equal,
            "the ratio is {GOLDEN_SECTION_RATIO} and (sqrt(5) - 1) / 2 is {of_pynei}"
        );
        assert_eq!(GOLDEN_SECTION_STEPS, 60, "the steps of `_reml_delta`");
        let mut width = 0.4_f64;
        for _ in 0..GOLDEN_SECTION_STEPS {
            width *= GOLDEN_SECTION_RATIO;
        }
        let apart = (width - THE_BRACKET_AFTER_THE_SEARCH).abs();
        assert!(
            apart <= 1e-24,
            "sixty steps leave a bracket of {width} where Python leaves \
             {THE_BRACKET_AFTER_THE_SEARCH}, {apart} away"
        );
    }

    /// The `delta` the search lands on for the panel is pyNei's, within
    /// 2.5e-7 relative.
    ///
    /// `delta` is the ratio of the two variances and is where the whole
    /// search ends up: the two variances, the heritability, the projection
    /// matrix and every variant's answer are built from it. The comparison
    /// with GMMAT is of the two variances at 1e-5 absolute, which is a
    /// looser statement about the same number, and the `y' p y` identity
    /// below says nothing about it at all.
    ///
    /// The comment on [`OF_PYNEI_THE_SEARCH`] has what this bound catches
    /// and what it cannot.
    #[test]
    fn the_delta_the_search_lands_on_for_the_panel_is_pyneis() {
        let fitted = the_null_of_the_panel("panel_called");

        // The two variances are read off the fit and not off its
        // `NullModel`, so that what this asserts stays the search's own
        // answer whatever a later item of the spec does with the fields a
        // user reads.
        let found = fitted.residual_variance / fitted.genetic_variance;
        let apart = (found - OF_PYNEI_DELTA).abs() / OF_PYNEI_DELTA;
        assert!(
            apart <= OF_PYNEI_THE_SEARCH,
            "the search landed on a delta of {found} and pyNei's `_reml_delta` gives \
             {OF_PYNEI_DELTA}, {apart} of it away against the {OF_PYNEI_THE_SEARCH} allowed"
        );
    }

    /// A trait that the kinship explains nothing of puts the search at the
    /// top end of the grid, where the best point is bracketed by that end
    /// and the one point below it.
    ///
    /// The clamping of the best point's neighbours at the ends of the grid
    /// is one of the six things "The linear mixed model" of
    /// `docs/specs/gwas.md` asks to be reproduced, and nothing reached it:
    /// measured on 25 September 2026 over every fixture of both suites, the
    /// best index is 44, 46, 47 or 85 and never 0 or 100. A trait of pure
    /// noise does reach it, and it is not a contrived input: it is what a
    /// study of a trait the panel's relatedness has nothing to do with
    /// gives.
    ///
    /// `delta` is the residual variance over the genetic one, so a trait
    /// with no genetic part sends it up. The grid's top is `exp(10)`,
    /// 22026.5, and the search cannot go past it: without the clamp the
    /// bracket would run to `exp(10.2)` and the fit would land somewhere
    /// pyNei's never does. Measured the same day, the trait below gives a
    /// `log(delta)` of 10.0000 and a heritability of 4.54e-5.
    ///
    /// The trait is `sin(at * 7.3)` over the 200 individuals, which is a
    /// fixed sequence that follows nothing of the kinship; `sin` is not
    /// rounded the same on every platform and nothing here depends on its
    /// last bits.
    #[test]
    fn a_trait_the_kinship_explains_nothing_of_lands_at_the_end_of_the_grid() {
        let individuals = the_individuals_of_the_kinship("panel_called");
        let kinship = the_kinship_of("panel_called");
        let phenotype: Vec<f64> = (0..individuals.len())
            .map(|at| (at as f64 * 7.3).sin())
            .collect();
        // The intercept alone: a covariate would take a part of a trait
        // that is meant to be nothing but noise.
        let design = vec![1.0_f64; individuals.len()];
        let tested: Vec<usize> = (0..individuals.len()).collect();
        let study = GwasInput {
            phenotype: &phenotype,
            trait_type: TraitType::Continuous,
            design: &design,
            num_coefs: 1,
            kinship: Some(&kinship),
            test: None,
            use_grammar_gamma_approx: false,
            individuals: &tested,
            transform_to_biallelic: false,
        };
        let design = match Design::of_the_study(&study, individuals.len()) {
            Ok(design) => design,
            Err(error) => panic!("the design of the trait of noise: {error}"),
        };
        let fitted = match LinearMixedModel::of_the_study(&phenotype, &design, &kinship) {
            Ok(fitted) => fitted,
            Err(error) => panic!("the null model of the trait of noise: {error}"),
        };

        let log_delta = (fitted.residual_variance / fitted.genetic_variance).ln();
        let top = the_log_delta_at(LOG_DELTA_POINTS.saturating_sub(1));
        assert!(
            (log_delta - top).abs() <= 1e-9,
            "the search landed at a log delta of {log_delta} and the grid's top is \
             {top}, which is where a trait the kinship explains nothing of belongs"
        );
        assert!(
            log_delta <= top,
            "the search landed past the top of the grid, at {log_delta}, which is the \
             clamp of the best point's neighbours not holding"
        );
    }

    /// How many variants a study reads in one block,
    /// [`crate::block::MAX_NUM_VARS_PER_BLOCK`], which does not depend on
    /// how many individuals it has.
    const VARS_OF_ONE_BLOCK: usize = 10_000;

    /// A mixed model over more variants than one block holds answers the
    /// same in its second block as in its first.
    ///
    /// The linear model has this test and the mixed one had none: 200
    /// individuals give 10000 variants to a block, so the 1200 variant
    /// panel is one block in every other test of this module and the loop
    /// over the blocks runs once. What it covers is the buffers of the
    /// model and of the dosages being reused, the rows of the second block
    /// being added after the first's and not over them, and the columns of
    /// the variants growing across the blocks. Commenting out the three
    /// `clear()` calls of [`LinearMixedModel::test_the_block`] passes every
    /// other test of the crate and fails this one.
    ///
    /// The two patterns are those of the fixture below, one that the
    /// projection leaves nothing of and one it answers, so the study
    /// carries a refused variant and an answered one into its second
    /// block. Every number asserted is the answer the same pattern got in
    /// the first block, which is what the test is about; what those numbers
    /// are is the fixture below and the six literals of the panel.
    #[test]
    fn a_mixed_study_of_more_variants_than_one_block_answers_the_same_in_every_block() {
        let patterns = [
            ["0/0", "0/0", "0/0", "0/0", "1/1", "1/1", "1/1", "1/1"],
            ["0/1", "0/0", "1/1", "0/0", "0/1", "1/1", "0/0", "0/1"],
        ];
        let num_vars = VARS_OF_ONE_BLOCK.saturating_add(100);
        let mut vcf = String::from(THE_HEADER_OF_EIGHT);
        for var in 0..num_vars {
            let pos = var.saturating_add(1).saturating_mul(10);
            vcf.push_str(&format!("1\t{pos}\tv{var}\tA\tT\t.\t.\t.\tGT"));
            for genotype in patterns[var % patterns.len()] {
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
            kinship: Some(&THE_KINSHIP_OF_TWO_SUBPOPULATIONS),
            test: Some(TestType::Wald),
            use_grammar_gamma_approx: false,
            individuals: &THE_INDIVIDUALS_OF_TWO_SUBPOPULATIONS,
            transform_to_biallelic: false,
        };
        let mut reader = reader_over(vcf.as_bytes());
        let result = match the_study_of(&mut reader, &study) {
            Ok(result) => result,
            Err(error) => panic!("the mixed study of {num_vars} variants: {error}"),
        };

        assert!(
            num_vars > VARS_OF_ONE_BLOCK,
            "the study has to read more than one block"
        );
        assert_eq!(result.num_vars, num_vars, "the variants of the study");
        let ids = result.ids.as_deref().expect("the ids of the variants");
        assert_eq!(ids.len(), num_vars, "one id for each variant");
        for var in [0, VARS_OF_ONE_BLOCK, num_vars.saturating_sub(1)] {
            assert_eq!(ids[var], format!("v{var}"), "the id of the variant {var}");
        }
        // Every variant of the first pattern is the one the projection
        // leaves nothing of, wherever it fell, and every variant of the
        // second has the answer the second variant of the fixture below
        // got. The answers of the first block are what the rest are
        // compared with, so nothing here is a literal of this test.
        let (of_the_refused, of_the_answered) = (0_usize, 1_usize);
        for var in 0..num_vars {
            let of_the_pattern = match var % patterns.len() {
                0 => of_the_refused,
                _ => of_the_answered,
            };
            if of_the_pattern == of_the_refused {
                assert!(
                    result.beta[var].is_nan()
                        && result.se[var].is_nan()
                        && result.p_value[var].is_nan(),
                    "the variant {var}, which the projection leaves nothing of, was \
                     answered with a beta of {beta}",
                    beta = result.beta[var]
                );
                continue;
            }
            for (found, first, what) in [
                (
                    result.allele_freq[var],
                    result.allele_freq[of_the_pattern],
                    "the frequency",
                ),
                (result.beta[var], result.beta[of_the_pattern], "the effect"),
                (result.se[var], result.se[of_the_pattern], "the error"),
                (
                    result.p_value[var],
                    result.p_value[of_the_pattern],
                    "the p-value",
                ),
            ] {
                assert_eq!(
                    found.total_cmp(&first),
                    std::cmp::Ordering::Equal,
                    "{what} of the variant {var} is {found} and the same pattern in the \
                     first block answered {first}"
                );
            }
        }
    }

    /// A variant that the projection leaves nothing of has no answer under
    /// either test, which is the meanwhile of **Open 2** of
    /// `docs/specs/gwas.md`.
    ///
    /// The variant is twice the covariate, so the design explains all of
    /// it, and the projection matrix takes the design out of whatever it
    /// is applied to: `x' p x` is 0 in exact arithmetic and what is left is
    /// rounding, of whichever sign it fell on. `beta` would be a number
    /// divided by noise, as it is in the linear model, and the p-value
    /// would read as a variant that was tested and showed nothing.
    ///
    /// The threshold is the tested individuals times 2.2e-16 of the
    /// variant's own squared length times the largest value of the diagonal
    /// of the projection matrix, which is 1.78e-15 of that scale for these
    /// eight individuals. Measured on this fixture on 24 September 2026,
    /// the collinear variant keeps -6.2e-17 of it on Accelerate and 2.3e-16
    /// on faer, and the ordinary variant beside it keeps 0.576 on both. So
    /// the threshold sits 8 times above the largest rounding the two
    /// backends left and 15 orders of magnitude below a variant that has
    /// something to test.
    ///
    /// What the study gave before the threshold was there, measured the
    /// same day: on Accelerate `x' p x` came to -4.4e-16 and the Wald test
    /// answered a `beta` of 3 with an `se` of NaN, and on faer it came to
    /// 1.7e-15 and the answer was a `beta` of 2.13 with an `se` of 2.7e7
    /// and a p-value of 0.99999994, which reads as a variant that was
    /// tested and showed nothing.
    ///
    /// The ordinary variant is asserted to have an answer and not to any
    /// number: its effect is popnei's own, since no reference program was
    /// run on this fixture, and the six literals of the panel are what says
    /// the numbers are right.
    ///
    /// The size of the threshold is guarded here on faer alone, and the
    /// test above says why no fixture can guard it on both: with the
    /// threshold set to 0, measured on 25 September 2026, faer answers this
    /// variant with a `beta` of 2.13 and an `se` of 2.7e7 and Accelerate
    /// refuses it on the sign of its -4.44e-16.
    #[test]
    fn a_variant_the_projection_leaves_nothing_of_has_no_answer() {
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
        for test in [TestType::Wald, TestType::Score] {
            let study = GwasInput {
                phenotype: &THE_TRAIT_OF_TWO_SUBPOPULATIONS,
                trait_type: TraitType::Continuous,
                design: &THE_DESIGN_OF_TWO_SUBPOPULATIONS,
                num_coefs: 2,
                kinship: Some(&THE_KINSHIP_OF_TWO_SUBPOPULATIONS),
                test: Some(test),
                use_grammar_gamma_approx: false,
                individuals: &THE_INDIVIDUALS_OF_TWO_SUBPOPULATIONS,
                transform_to_biallelic: false,
            };
            let mut reader = reader_over(vcf.as_bytes());
            let result = match the_study_of(&mut reader, &study) {
                Ok(result) => result,
                Err(error) => panic!("the study of two subpopulations: {error}"),
            };
            assert_eq!(result.num_vars, 2, "the variants of the fixture");
            // The variant is in the result with its frequency, as every
            // variant that has no answer is.
            let found = result.allele_freq[0];
            assert!(
                (found - 0.5).abs() <= 1e-12,
                "the frequency of v0 is {found} and its dosages are half the alleles"
            );
            assert!(
                result.beta[0].is_nan() && result.se[0].is_nan() && result.p_value[0].is_nan(),
                "the variant the projection leaves nothing of was answered under the \
                 {test:?} test with a beta of {beta}, an se of {se} and a p-value of \
                 {p_value}",
                beta = result.beta[0],
                se = result.se[0],
                p_value = result.p_value[0]
            );
            assert!(
                result.beta[1].is_finite()
                    && result.se[1] > 0.0
                    && (0.0..=1.0).contains(&result.p_value[1]),
                "the ordinary variant beside it was answered under the {test:?} test with \
                 a beta of {beta}, an se of {se} and a p-value of {p_value}",
                beta = result.beta[1],
                se = result.se[1],
                p_value = result.p_value[1]
            );
        }
    }
}
