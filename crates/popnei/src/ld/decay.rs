//! The curve of r² against distance fitted to the pairs a population
//! counted, and the distance at which that curve has fallen to half, which
//! is "The curve that is fitted" of the item "LD against distance, per
//! population" of `docs/specs/ld.md`.
//!
//! [`fit_ld_decay`] takes what a pass counted, for each distance that holds
//! a pair the number of pairs and the sum of their r², and gives an
//! [`LdDecay`]: the fitted ρ per base pair, the curve at a distance of 0
//! and the distance at which the curve is half of that. It reads no
//! genotype and no bin, so the same pairs give the same curve whatever
//! `num_bins` a user asked the bins to be cut into.
//!
//! ρ is the scaled recombination between the two variants of a pair, four
//! times the effective size of the population times the recombination
//! fraction between them. The model is the r² two variants of a population
//! under drift and recombination are expected to be in, Hill and Weir
//! (1988) with the correction for a finite sample of Weir and Hill (1986),
//! and its one free number here is the ρ per base pair: a pair d base pairs
//! apart has a ρ of d times it. The spec gives the reasons for that model
//! and for n being the individuals of the population and not the
//! individuals times the ploidy.
//!
//! There is one number to fit, so the smallest is found without
//! derivatives: the sum is evaluated at 141 values of the ρ per base pair,
//! from 10⁻¹² to 10², spaced by a tenth of a decade, and the two
//! neighbours of the smallest of them bracket a golden section search. The
//! grid before the search is what looks at the whole range instead of
//! sliding downhill from a start value into whichever valley holds it, and
//! what is kept is the smallest of every ρ per base pair evaluated, the
//! 141 of the grid included.

use crate::error::{Error, Result};

/// The exponent of 10 of the smallest ρ per base pair the fit looks at.
///
/// Below 10⁻¹² the curve is flat across 10⁶ base pairs, which is the
/// largest `max_dist` a user would give, so nothing below it is a fall-off
/// that pairs of one dataset could pin down.
const THE_SMALLEST_EXPONENT_OF_THE_GRID: f64 = -12.0;

/// How far apart the exponents of two neighbouring values of the grid are,
/// a tenth of a decade.
const THE_STEP_OF_THE_GRID: f64 = 0.1;

/// How many values of the ρ per base pair the sum is evaluated at before
/// the search, from 10⁻¹² to 10² by a tenth of a decade, both ends
/// included.
const THE_POINTS_OF_THE_GRID: usize = 141;

/// How narrow the bracket of the golden section search is left, in
/// decades.
///
/// A bracket two cells of the grid wide is 0.2 decades, and each step of
/// the search leaves 0.618 of it, so it takes 40 steps to come under this.
/// The fit therefore costs 183 evaluations of the sum, 141 for the grid, 2
/// to open the search and 40 for its steps, each one pass over the
/// distances that hold a pair.
const THE_TOLERANCE_OF_THE_SEARCH: f64 = 1e-9;

/// The largest ρ the half distance is looked for below.
///
/// The ρ at which the curve is half of its value at ρ of 0 is 2.16 at n of
/// 100 and 2.26 at n of 50, and it grows as n falls, so 10⁶ is far above
/// any n a dataset has.
const THE_LARGEST_RHO_OF_THE_HALF: f64 = 1e6;

/// How narrow the bracket of the bisection of the half distance is left,
/// as a part of its own middle. It is about 60 halvings from a bracket
/// that runs from 0 to [`THE_LARGEST_RHO_OF_THE_HALF`].
const THE_TOLERANCE_OF_THE_HALF: f64 = 1e-12;

/// The curve of r² against distance fitted to the pairs of one population,
/// and the distance at which it has fallen to half.
///
/// The three values are NaN together when no curve was fitted, which "The
/// cases" of `docs/specs/ld.md` says happens to a population whose pairs
/// fall at fewer than two distances and to one whose smallest sum falls at
/// either end of the searched range of the ρ per base pair.
#[derive(Debug, Clone, Copy)]
pub struct LdDecay {
    /// The fitted 4Nr, by how much ρ grows per base pair.
    rho_per_bp: f64,
    /// The fitted curve at a distance of 0.
    r2_at_zero: f64,
    /// The distance at which the fitted curve is half of `r2_at_zero`.
    half_dist: f64,
}

impl LdDecay {
    /// The fitted 4Nr, by how much the scaled recombination ρ grows per
    /// base pair. NaN when no curve was fitted, which "The cases" of
    /// `docs/specs/ld.md` says when, and then the other two are NaN as
    /// well.
    pub fn rho_per_bp(&self) -> f64 {
        self.rho_per_bp
    }

    /// The fitted curve at a distance of 0, which the individuals of the
    /// population fix on their own: it is the curve's own ceiling,
    /// 0.46198347107438015 at 100 individuals, and no pair of the dataset
    /// moves it.
    pub fn r2_at_zero(&self) -> f64 {
        self.r2_at_zero
    }

    /// The distance in base pairs at which the fitted curve has fallen to
    /// half of [`LdDecay::r2_at_zero`].
    ///
    /// NaN when the other two are, and NaN on its own at 1 and at 2
    /// individuals, the only counts of individuals whose curve never falls
    /// to half: what it falls towards as ρ grows is 1 divided by the
    /// individuals, and half of the value at 0 is above that from 3
    /// individuals upwards.
    pub fn half_dist(&self) -> f64 {
        self.half_dist
    }

    /// The three NaN of a population no curve was fitted to, which is also
    /// what the bins of a pass hold until the pass has ended and
    /// [`fit_ld_decay`] has read the pairs of each distance.
    pub(super) fn of_no_curve() -> LdDecay {
        LdDecay {
            rho_per_bp: f64::NAN,
            r2_at_zero: f64::NAN,
            half_dist: f64::NAN,
        }
    }
}

/// The curve fitted to the pairs counted at each distance.
///
/// `dists` are in base pairs and need not be in order, `num_pairs` and
/// `sum_r2` hold one value for each of them, and `num_individuals` is how
/// many individuals the population has, the n of the model. Pairs at fewer
/// than two distances, which includes no pair at all, and a smallest sum
/// that falls at either end of the searched range of the ρ per base pair
/// are the [`LdDecay`] of three NaN that "The cases" of `docs/specs/ld.md`
/// describes and not an error.
///
/// # Errors
///
/// [`Error::LdDecayArraysOfDifferentLengths`] when the three slices are not
/// of one length, [`Error::LdDecayNoIndividuals`] when `num_individuals` is
/// 0, [`Error::LdDecayDistWithNoPair`] when a distance holds no pair, and
/// [`Error::LdDecaySumOfR2OutOfRange`] when a sum of r² is not finite or is
/// below 0.
pub fn fit_ld_decay(
    dists: &[u64],
    num_pairs: &[u64],
    sum_r2: &[f64],
    num_individuals: u64,
) -> Result<LdDecay> {
    the_pairs_of_each_dist_are_checked(dists, num_pairs, sum_r2, num_individuals)?;
    if dists.len() < 2 {
        // One distance says nothing about a fall-off, whatever a search
        // would return for it, and no pair at all is that case.
        return Ok(LdDecay::of_no_curve());
    }
    let individuals = num_individuals as f64;
    let Some(rho_per_bp) = the_rho_per_bp_of_the_smallest(dists, num_pairs, sum_r2, individuals)
    else {
        // The smallest fell at an end of the searched range, so the number
        // there says where the search stopped and not what the pairs say.
        return Ok(LdDecay::of_no_curve());
    };
    Ok(LdDecay {
        rho_per_bp,
        r2_at_zero: the_curve_at(0.0, individuals),
        half_dist: match the_rho_at_half(individuals) {
            Some(rho) => rho / rho_per_bp,
            None => f64::NAN,
        },
    })
}

/// The r² the model expects at the scaled recombination `rho` between the
/// two variants of a pair, in a population of `num_individuals`
/// individuals.
///
/// It is Hill and Weir (1988) with the correction for a finite sample of
/// Weir and Hill (1986), which "The curve that is fitted" of
/// `docs/specs/ld.md` writes as
///
/// ```text
/// E[r²] = (10 + ρ) / ((2 + ρ) · (11 + ρ))
///         · [1 + ((3 + ρ) · (12 + 12ρ + ρ²)) / (n · (2 + ρ) · (11 + ρ))]
/// ```
///
/// with ρ the scaled recombination and n the individuals. The first factor
/// falls from 10/22 at ρ of 0 towards 0, and the second is what holds the
/// curve up at long distances, where the r² of a finite sample does not
/// fall to 0.
///
/// It is written with the four operations alone, which are rounded the
/// same on every platform where `exp`, `ln` and `powf` are not, as the
/// `coding` skill says.
fn the_curve_at(rho: f64, num_individuals: f64) -> f64 {
    let two_plus_rho = 2.0 + rho;
    let eleven_plus_rho = 11.0 + rho;
    let expected = (10.0 + rho) / (two_plus_rho * eleven_plus_rho);
    let of_the_sample = ((3.0 + rho) * (12.0 + 12.0 * rho + rho * rho))
        / (num_individuals * two_plus_rho * eleven_plus_rho);
    expected * (1.0 + of_the_sample)
}

/// The sum the fit makes smallest at one ρ per base pair,
///
/// ```text
/// Σ over the distances of [ n_d · (m_d − f(d))² ]
/// ```
///
/// where n_d is how many pairs the distance d holds, m_d the mean of their
/// r², their sum divided by n_d, and f(d) the curve at that distance. A
/// distance that holds no pair is refused before the fit, so n_d is 1 or
/// more and the mean is a division by a positive number.
///
/// "The curve that is fitted" of `docs/specs/ld.md` derives it: the sum
/// over every pair of the square of its r² less the curve at its distance
/// is the spread of the pairs of a distance around their own mean, which
/// no ρ changes, plus what this sum holds. So the ρ that makes this
/// smallest is the ρ that makes the squared residuals of every pair
/// smallest, and the pairs of one distance are added up before the fit
/// with nothing lost.
///
/// Multiplying the square out into Σ [ n_d · f(d)² − 2 · S_d · f(d) ],
/// with S_d the sum of the r² of the distance, gives this number less
/// Σ n_d · m_d², which no ρ changes and which the same part of the spec
/// measures at 475.36 on the first of its three populations, where this
/// sum is 23.21 at the fitted ρ. The terms that depend on ρ would then be
/// added inside a total 19.5 times their size, and what that costs the
/// fitted ρ per base pair is in the spec.
///
/// The distances are read in the order they were given, which is the order
/// a pass compacted them in, so two runs over one dataset add the same
/// numbers in the same order.
fn the_sum_to_make_smallest(
    rho_per_bp: f64,
    dists: &[u64],
    num_pairs: &[u64],
    sum_r2: &[f64],
    num_individuals: f64,
) -> f64 {
    let mut total = 0.0;
    for ((dist, pairs), sum) in dists.iter().zip(num_pairs).zip(sum_r2) {
        let curve = the_curve_at(*dist as f64 * rho_per_bp, num_individuals);
        let mean = *sum / *pairs as f64;
        let apart = mean - curve;
        total += *pairs as f64 * apart * apart;
    }
    total
}

/// The smallest ρ per base pair the fit has evaluated so far, and the sum
/// there.
struct TheSmallestSoFar {
    /// The exponent of 10 of that ρ per base pair, which is what the grid
    /// and the search both move in.
    exponent: f64,
    /// The sum at it.
    sum: f64,
}

impl TheSmallestSoFar {
    /// Keeps the exponent when its sum is below the smallest so far, and
    /// says whether it kept it.
    ///
    /// A sum equal to the smallest does not replace it, so of two
    /// exponents whose sums are the same bits the one evaluated first
    /// wins, which is the lower of the grid, the grid being evaluated
    /// before the search. It is the rule the `coding` skill gives for a
    /// comparison that decides a result.
    fn take(&mut self, exponent: f64, sum: f64) -> bool {
        if sum < self.sum {
            self.exponent = exponent;
            self.sum = sum;
            return true;
        }
        false
    }
}

/// The ρ per base pair at the smallest of that sum, and `None` when the
/// smallest of the grid falls at either end of the searched range.
///
/// The 141 values of the grid are looked at first, which is what looks at
/// the whole range instead of sliding downhill from a start value into
/// whichever valley holds it. The two neighbours of the smallest of them
/// bracket a golden section search, which holds two points inside the
/// bracket and drops the end beyond whichever of them has the larger sum,
/// so the bracket is 0.618 of itself after each step. What is given back
/// is the smallest of every ρ per base pair evaluated, the 141 of the grid
/// included, so a sum with more than one valley inside the cell the search
/// works on gives at worst the best of the grid, a tenth of a decade from
/// the smallest.
///
/// The grid and the search both move in the exponent of 10 of the ρ per
/// base pair and not in the ρ per base pair itself: the values of the grid
/// are a tenth of a decade apart, and the search stops when its bracket is
/// narrower than [`THE_TOLERANCE_OF_THE_SEARCH`] of a decade.
fn the_rho_per_bp_of_the_smallest(
    dists: &[u64],
    num_pairs: &[u64],
    sum_r2: &[f64],
    num_individuals: f64,
) -> Option<f64> {
    let sum_at = |exponent: f64| {
        the_sum_to_make_smallest(
            the_rho_per_bp_of(exponent),
            dists,
            num_pairs,
            sum_r2,
            num_individuals,
        )
    };

    let mut the_smallest_point = 0;
    let mut smallest = TheSmallestSoFar {
        exponent: the_exponent_of(0),
        sum: f64::INFINITY,
    };
    for point in 0..THE_POINTS_OF_THE_GRID {
        let exponent = the_exponent_of(point);
        if smallest.take(exponent, sum_at(exponent)) {
            the_smallest_point = point;
        }
    }

    // The two neighbours of the smallest of the grid, which do not exist
    // when it is the first or the last of the 141: a curve that is flat
    // across `max_dist`, or that has fallen before the second base pair,
    // is not a fall-off these pairs pin down.
    let below = the_smallest_point.checked_sub(1)?;
    let above = the_smallest_point.checked_add(1)?;
    if above >= THE_POINTS_OF_THE_GRID {
        return None;
    }

    let golden = the_golden_section();
    let mut lower = the_exponent_of(below);
    let mut upper = the_exponent_of(above);
    let mut inside_low = upper - golden * (upper - lower);
    let mut inside_high = lower + golden * (upper - lower);
    let mut at_low = sum_at(inside_low);
    let mut at_high = sum_at(inside_high);
    let _ = smallest.take(inside_low, at_low);
    let _ = smallest.take(inside_high, at_high);
    while (upper - lower) > THE_TOLERANCE_OF_THE_SEARCH {
        if at_low < at_high {
            upper = inside_high;
            inside_high = inside_low;
            at_high = at_low;
            inside_low = upper - golden * (upper - lower);
            at_low = sum_at(inside_low);
            let _ = smallest.take(inside_low, at_low);
        } else {
            lower = inside_low;
            inside_low = inside_high;
            at_low = at_high;
            inside_high = lower + golden * (upper - lower);
            at_high = sum_at(inside_high);
            let _ = smallest.take(inside_high, at_high);
        }
    }
    Some(the_rho_per_bp_of(smallest.exponent))
}

/// The exponent of 10 of the value of the grid at `point`, from
/// [`THE_SMALLEST_EXPONENT_OF_THE_GRID`] by [`THE_STEP_OF_THE_GRID`].
fn the_exponent_of(point: usize) -> f64 {
    THE_SMALLEST_EXPONENT_OF_THE_GRID + point as f64 * THE_STEP_OF_THE_GRID
}

/// The ρ per base pair whose exponent of 10 is `exponent`.
///
/// It is the one call of the fit to a function that is not rounded the
/// same on every platform, 183 of them in a fit, and what it moves is
/// where the search looks and not what the curve gives there: two
/// platforms that place a value of the grid a last bit apart still stop
/// within [`THE_TOLERANCE_OF_THE_SEARCH`] of a decade of one smallest.
fn the_rho_per_bp_of(exponent: f64) -> f64 {
    10.0_f64.powf(exponent)
}

/// (√5 − 1) / 2, the part of a bracket the golden section search leaves
/// after each step, 0.618.
///
/// `sqrt` is exact on every platform, which the `coding` skill says of it
/// and of the four operations and of nothing else, so the search steps the
/// same everywhere.
fn the_golden_section() -> f64 {
    (5.0_f64.sqrt() - 1.0) / 2.0
}

/// The ρ at which the curve of a population of `num_individuals`
/// individuals is half of what it is at ρ of 0, and `None` when the curve
/// never falls to half.
///
/// It is 2.1608135872529166 at 100 individuals and 2.2641731329247312 at
/// 50, and it is solved for by bisection between ρ of 0 and
/// [`THE_LARGEST_RHO_OF_THE_HALF`] and not read from a table, so a dataset
/// of another n needs no new number. The bracket is left when it is
/// narrower than [`THE_TOLERANCE_OF_THE_HALF`] of its own middle.
///
/// The curve falls towards 1 divided by the individuals and not towards 0,
/// so half of its value at ρ of 0 is below what it ever reaches at 1 and
/// at 2 individuals, which "The curve that is fitted" of
/// `docs/specs/ld.md` gives the numbers of, and those two are the `None`.
fn the_rho_at_half(num_individuals: f64) -> Option<f64> {
    let half = the_curve_at(0.0, num_individuals) / 2.0;
    let mut lower = 0.0;
    let mut upper = THE_LARGEST_RHO_OF_THE_HALF;
    if the_curve_at(upper, num_individuals) >= half {
        return None;
    }
    while (upper - lower) > THE_TOLERANCE_OF_THE_HALF * ((lower + upper) / 2.0) {
        let middle = (lower + upper) / 2.0;
        if the_curve_at(middle, num_individuals) > half {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    Some((lower + upper) / 2.0)
}

/// What [`fit_ld_decay`] refuses before it fits anything.
///
/// # Errors
///
/// The four of [`fit_ld_decay`].
fn the_pairs_of_each_dist_are_checked(
    dists: &[u64],
    num_pairs: &[u64],
    sum_r2: &[f64],
    num_individuals: u64,
) -> Result<()> {
    if dists.len() != num_pairs.len() || dists.len() != sum_r2.len() {
        return Err(Error::LdDecayArraysOfDifferentLengths {
            num_dists: dists.len(),
            num_pairs: num_pairs.len(),
            num_sums: sum_r2.len(),
        });
    }
    if num_individuals == 0 {
        return Err(Error::LdDecayNoIndividuals);
    }
    for ((dist, pairs), sum) in dists.iter().zip(num_pairs).zip(sum_r2) {
        if *pairs == 0 {
            return Err(Error::LdDecayDistWithNoPair { dist: *dist });
        }
        if !sum.is_finite() || *sum < 0.0 {
            return Err(Error::LdDecaySumOfR2OutOfRange {
                dist: *dist,
                sum_r2: *sum,
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::ld::LdAndDist;
    use crate::ld::dist::tests::{bins_of, the_three_tables_of};

    /// How close two numbers of the fit are asked to be, which is the
    /// tolerance "How it is verified" of `docs/specs/ld.md` compares the
    /// fitted values within: 480 times the 2.1·10⁻⁹ that R's `optimize`
    /// and R's `nls` disagree by on the same pairs.
    const THE_TOLERANCE: f64 = 1e-6;

    /// How close the r² at a distance of 0 is asked to be, which the same
    /// part of the spec compares that one value within.
    ///
    /// It is the curve's own ceiling, which the individuals of the
    /// population fix on their own and no pair moves, so the two sides
    /// work one formula out at a ρ of 0 and nothing but the last bits can
    /// differ. The fitted ρ per base pair, which the two optimisers place
    /// 2.1·10⁻⁹ apart, does not enter it.
    const THE_TOLERANCE_OF_THE_R2_AT_ZERO: f64 = 1e-12;

    /// Asserts that each of `said` is in the message of the error, which
    /// is what a caller of [`fit_ld_decay`] with a table of its own reads:
    /// the argument that was refused and the value it was given.
    #[track_caller]
    fn assert_the_message_says(error: &Error, said: &[&str]) {
        let message = error.to_string();
        for what in said {
            assert!(message.contains(what), "`{what}` is not in `{message}`");
        }
    }

    /// Asserts that `found` is within `tolerance` of `expected`, relative
    /// to `expected`.
    #[track_caller]
    fn assert_within(found: f64, expected: f64, tolerance: f64, what: &str) {
        let apart = (found - expected).abs() / expected.abs();
        assert!(
            apart <= tolerance,
            "{what} is {found} and the spec gives {expected}, {apart} of it apart"
        );
    }

    /// Asserts that `found` is within [`THE_TOLERANCE`] of `expected`,
    /// relative to `expected`.
    #[track_caller]
    fn assert_close(found: f64, expected: f64, what: &str) {
        assert_within(found, expected, THE_TOLERANCE, what);
    }

    /// The r² the curve itself gives at each distance of `dists`, with the
    /// ρ per base pair and the individuals given, which is the table "How
    /// it is verified" of `docs/specs/ld.md` makes the first test of the
    /// fit on: one pair at each distance, whose r² is the curve there, so
    /// the answer is known before the fit runs.
    fn the_table_of_the_curve(
        dists: &[u64],
        rho_per_bp: f64,
        num_individuals: u64,
    ) -> (Vec<u64>, Vec<f64>) {
        let num_pairs = vec![1_u64; dists.len()];
        let sum_r2 = dists
            .iter()
            .map(|dist| the_curve_at(*dist as f64 * rho_per_bp, num_individuals as f64))
            .collect();
        (num_pairs, sum_r2)
    }

    /// The distances `step`, 2·`step` and so on up to `last`.
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "the tests give a step of 500 or 1000 and a last distance of at most \
                  1000000, so no product of the two reaches two million"
    )]
    fn the_dists_from(step: u64, last: u64) -> Vec<u64> {
        (1..)
            .map(|which| which * step)
            .take_while(|dist| *dist <= last)
            .collect()
    }

    #[test]
    fn the_fit_gives_back_the_rho_per_bp_the_table_of_the_curve_was_made_with() {
        let dists = the_dists_from(1000, 250_000);
        assert_eq!(dists.len(), 250, "the distances of the table");
        let (num_pairs, sum_r2) = the_table_of_the_curve(&dists, 0.0001, 100);
        let decay = fit_ld_decay(&dists, &num_pairs, &sum_r2, 100).expect("the fit");
        assert_close(decay.rho_per_bp(), 0.0001, "the fitted rho per base pair");
        assert_close(decay.half_dist(), 21608.135872529165, "the half distance");
        assert_close(
            decay.r2_at_zero(),
            0.46198347107438015,
            "the r² at a distance of 0",
        );
    }

    #[test]
    fn the_fit_finds_a_rho_per_bp_that_is_not_one_of_the_141_of_the_grid() {
        let dists = the_dists_from(500, 125_000);
        let (num_pairs, sum_r2) = the_table_of_the_curve(&dists, 0.0002, 50);
        let decay = fit_ld_decay(&dists, &num_pairs, &sum_r2, 50).expect("the fit");
        assert_close(decay.rho_per_bp(), 0.0002, "the fitted rho per base pair");
        // The spec writes the ρ that halves the curve at 50 individuals as
        // 2.2641731329247312, whose last digit is past what an `f64`
        // holds, and the half distance is it divided by the ρ per base
        // pair the table was made with.
        assert_close(
            decay.half_dist(),
            2.264_173_132_924_731 / 0.0002,
            "the half distance",
        );
        assert_close(
            decay.r2_at_zero(),
            0.46942148760330576,
            "the r² at a distance of 0",
        );
    }

    #[test]
    fn the_fit_gives_back_a_fall_off_that_takes_the_whole_of_the_default_max_dist() {
        let dists = the_dists_from(1000, 1_000_000);
        assert_eq!(dists.len(), 1000, "the distances of the table");
        let (num_pairs, sum_r2) = the_table_of_the_curve(&dists, 4.17e-8, 100);
        let decay = fit_ld_decay(&dists, &num_pairs, &sum_r2, 100).expect("the fit");
        assert_close(decay.rho_per_bp(), 4.17e-8, "the fitted rho per base pair");
        // The spec writes the ρ that halves the curve at 100 individuals
        // as 2.1608135872529166, and the half distance is it divided by
        // the ρ per base pair the table was made with.
        assert_close(
            decay.half_dist(),
            2.160_813_587_252_916_6 / 4.17e-8,
            "the half distance",
        );
        assert_close(
            decay.r2_at_zero(),
            0.46198347107438015,
            "the r² at a distance of 0",
        );
    }

    /// The three values of `decay` are NaN, which is what a population no
    /// curve was fitted to gives.
    #[track_caller]
    fn assert_no_curve(decay: &LdDecay, what: &str) {
        assert!(
            decay.rho_per_bp().is_nan()
                && decay.r2_at_zero().is_nan()
                && decay.half_dist().is_nan(),
            "{what} gave a rho per base pair of {rho}, an r² at 0 of {zero} and a half \
             distance of {half}, where all three are NaN",
            rho = decay.rho_per_bp(),
            zero = decay.r2_at_zero(),
            half = decay.half_dist()
        );
    }

    #[test]
    fn the_last_value_of_the_grid_is_the_top_of_the_searched_range() {
        assert_close(
            the_exponent_of(0),
            THE_SMALLEST_EXPONENT_OF_THE_GRID,
            "the exponent of the first value of the grid",
        );
        let last = THE_POINTS_OF_THE_GRID
            .checked_sub(1)
            .expect("the grid has a value");
        assert_close(
            the_exponent_of(last),
            2.0,
            "the exponent of the last value of the grid",
        );
    }

    #[test]
    fn pairs_at_one_distance_have_no_curve() {
        let decay = fit_ld_decay(&[1000], &[42], &[7.5], 100).expect("the fit");
        assert_no_curve(&decay, "pairs at one distance");
    }

    #[test]
    fn no_pair_at_all_has_no_curve() {
        let decay = fit_ld_decay(&[], &[], &[], 100).expect("the fit");
        assert_no_curve(&decay, "no pair at all");
    }

    #[test]
    fn a_smallest_at_the_bottom_end_of_the_range_has_no_curve() {
        let dists = the_dists_from(1000, 250_000);
        let num_pairs = vec![1_u64; dists.len()];
        let sum_r2 = vec![1.0_f64; dists.len()];
        let decay = fit_ld_decay(&dists, &num_pairs, &sum_r2, 100).expect("the fit");
        assert_no_curve(&decay, "an r² of 1 at every distance");
    }

    #[test]
    fn a_smallest_at_the_top_end_of_the_range_has_no_curve() {
        let dists = the_dists_from(1000, 250_000);
        let num_pairs = vec![1_u64; dists.len()];
        let sum_r2 = vec![0.0_f64; dists.len()];
        let decay = fit_ld_decay(&dists, &num_pairs, &sum_r2, 100).expect("the fit");
        assert_no_curve(&decay, "an r² of 0 at every distance");
    }

    #[test]
    fn a_population_of_two_individuals_has_no_half_distance_and_keeps_the_other_two() {
        let dists = the_dists_from(1000, 250_000);
        let (num_pairs, sum_r2) = the_table_of_the_curve(&dists, 0.0001, 2);
        let decay = fit_ld_decay(&dists, &num_pairs, &sum_r2, 2).expect("the fit");
        assert_close(decay.rho_per_bp(), 0.0001, "the fitted rho per base pair");
        assert_close(
            decay.r2_at_zero(),
            0.8264462809917356,
            "the r² at a distance of 0",
        );
        assert!(
            decay.half_dist().is_nan(),
            "the half distance of two individuals is {half}, where the curve runs from \
             0.8264462809917356 down to 0.5 and never reaches half of the first",
            half = decay.half_dist()
        );
    }

    /// The three rows of the table of "How it is verified" of
    /// `docs/specs/ld.md`: for each population, what it is called, the
    /// fitted ρ per base pair, the r² at a distance of 0 and the half
    /// distance in base pairs.
    ///
    /// Every value is the one R 4.6.1's `optimize` gives for the same
    /// pairs, which `docs/reports/ld-method/decay.py` groups by their
    /// exact distance out of the r² plink2 v2.0.0-a.7.7 gives and
    /// `docs/reports/ld-method/decay.R` fits the curve to. R's `nls`,
    /// which uses the derivatives popnei does not take, lands 2.1·10⁻⁹ of
    /// itself from `optimize` at the furthest of the three, and
    /// `tests/reference/ld/ld.decay.txt` holds what that script printed.
    ///
    /// The populations are the ones the bins of the same part are counted
    /// over, and the three curves are read out of the two passes those
    /// bins are read out of: the one population of every individual at a
    /// `max_allowed_maf` of 0.95, and `pop_a` and `pop_b` at 0.8. The n
    /// of the curve is the individuals of the population, 100 for the
    /// first and 50 for the other two, which is what makes the r² at a
    /// distance of 0 of the first differ from that of the two below it.
    ///
    /// R prints seventeen digits of each, and the last digit of four of
    /// them is past what an `f64` holds: the table gives the ρ per base
    /// pair of the first population as 0.00031727347196446889 and its
    /// half distance as 6810.5712522189806, `pop_a`'s half distance as
    /// 7530.1038938711654 and `pop_b`'s as 7259.8060755719744, and the
    /// four literals below are those four numbers as an `f64` holds them.
    const THE_CURVES_OF_THE_THREE_POPS: [(&str, f64, f64, f64); 3] = [
        (
            "every individual",
            0.000_317_273_471_964_468_9,
            0.461_983_471_074_380_15,
            6_810.571_252_218_981,
        ),
        (
            "pop_a",
            0.000_300_682_854_424_832_95,
            0.469_421_487_603_305_76,
            7_530.103_893_871_165,
        ),
        (
            "pop_b",
            0.000_311_877_908_218_966_46,
            0.469_421_487_603_305_76,
            7_259.806_075_571_974,
        ),
    ];

    /// The curve of each of the three populations of one run, in the order
    /// of the rows of [`THE_CURVES_OF_THE_THREE_POPS`].
    fn the_three_curves_of(of_the_run: &(LdAndDist, LdAndDist)) -> Vec<LdDecay> {
        let (of_every_individual, of_the_two_pops) = of_the_run;
        vec![
            *bins_of(of_every_individual, 0).decay(),
            *bins_of(of_the_two_pops, 0).decay(),
            *bins_of(of_the_two_pops, 1).decay(),
        ]
    }

    /// The three values of each curve of a run in the bits they came out
    /// with, which is what two runs over one dataset are compared by.
    fn the_bits_of_the_three_curves(of_the_run: &(LdAndDist, LdAndDist)) -> Vec<(u64, u64, u64)> {
        the_three_curves_of(of_the_run)
            .iter()
            .map(|decay| {
                (
                    decay.rho_per_bp().to_bits(),
                    decay.r2_at_zero().to_bits(),
                    decay.half_dist().to_bits(),
                )
            })
            .collect()
    }

    /// Asserts that the three curves of a run are the three rows of the
    /// table, and says of each which run it was.
    #[track_caller]
    fn assert_the_three_curves_are_the_ones_of_the_spec(
        of_the_run: &(LdAndDist, LdAndDist),
        at: &str,
    ) {
        let found = the_three_curves_of(of_the_run);
        for (decay, (what, rho_per_bp, r2_at_zero, half_dist)) in
            found.iter().zip(THE_CURVES_OF_THE_THREE_POPS)
        {
            assert_close(
                decay.rho_per_bp(),
                rho_per_bp,
                &format!("the rho per base pair of {what} {at}"),
            );
            assert_within(
                decay.r2_at_zero(),
                r2_at_zero,
                THE_TOLERANCE_OF_THE_R2_AT_ZERO,
                &format!("the r² at a distance of 0 of {what} {at}"),
            );
            assert_close(
                decay.half_dist(),
                half_dist,
                &format!("the half distance of {what} {at}"),
            );
        }
    }

    /// The three rows of the table of "How it is verified" of
    /// `docs/specs/ld.md`, over `tests/reference/ld/ld.vcf.gz` read with
    /// the VCF reader in blocks of 7, 64 and 500 variants, which give the
    /// same numbers to the bit.
    ///
    /// The fit reads the pairs the pass counted at each distance and no
    /// genotype, so what the three sizes of block say here is that those
    /// pairs are the same whatever the block: a window that kept a
    /// variant further back than `max_dist`, or dropped one a later
    /// variant still pairs with, would move the pairs of a distance and
    /// with them the ρ per base pair. A block of 7 holds far fewer
    /// variants than the window of 250000 base pairs reaches over, these
    /// variants being a thousand base pairs apart at the closest, and a
    /// block of 500 holds the whole dataset.
    #[test]
    fn the_three_curves_are_the_ones_r_gives_at_every_size_of_block() {
        let of_seven = the_three_tables_of(7);
        assert_the_three_curves_are_the_ones_of_the_spec(&of_seven, "at blocks of 7 variants");
        for num_vars_per_block in [64, 500] {
            let at = format!("at blocks of {num_vars_per_block} variants");
            let of_the_run = the_three_tables_of(num_vars_per_block);
            assert_the_three_curves_are_the_ones_of_the_spec(&of_the_run, &at);
            assert_eq!(
                the_bits_of_the_three_curves(&of_the_run),
                the_bits_of_the_three_curves(&of_seven),
                "{at}, against blocks of 7"
            );
        }
    }

    /// The three curves are the same, to the bit, on a pool of one thread
    /// and on one of four.
    ///
    /// The fit itself is not split across threads, and what the two pools
    /// move is the pass under it: the VCF reader parses the lines of a
    /// batch on the threads of the pool the caller is in, as
    /// `docs/specs/io_vcf.md` has it, and the products of r² run on them
    /// through faer when the `blas` feature is off, which is what
    /// `cargo test -p popnei --no-default-features` runs. The pools are
    /// built here and are not rayon's global one, which has one thread per
    /// core of the machine, and `current_num_threads` inside the pool says
    /// how many threads the pass had. rayon is a dependency of the targets
    /// that are not wasm, so this test is compiled for those alone.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_number_of_threads_does_not_change_the_three_curves() {
        let in_a_pool = |threads: usize| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            pool.install(|| {
                assert_eq!(
                    rayon::current_num_threads(),
                    threads,
                    "the pass did not run on the pool it was given"
                );
                the_three_tables_of(64)
            })
        };

        let on_one = in_a_pool(1);
        assert_the_three_curves_are_the_ones_of_the_spec(&on_one, "on one thread");
        assert_eq!(
            the_bits_of_the_three_curves(&in_a_pool(4)),
            the_bits_of_the_three_curves(&on_one),
            "four threads against one"
        );
    }

    #[test]
    fn arrays_of_different_lengths_are_an_error() {
        let error = fit_ld_decay(&[1000, 2000], &[1], &[0.5, 0.25], 100).expect_err("the error");
        assert!(
            matches!(
                error,
                Error::LdDecayArraysOfDifferentLengths {
                    num_dists: 2,
                    num_pairs: 1,
                    num_sums: 2
                }
            ),
            "the error of three arrays that are not of one length is {error:?}"
        );
        assert_the_message_says(
            &error,
            &["`dists` of 2", "`num_pairs` of 1", "`sum_r2` of 2"],
        );
    }

    #[test]
    fn a_population_of_no_individual_is_an_error() {
        let error = fit_ld_decay(&[1000, 2000], &[1, 1], &[0.5, 0.25], 0).expect_err("the error");
        assert!(
            matches!(error, Error::LdDecayNoIndividuals),
            "the error of a population of no individual is {error:?}"
        );
        assert_the_message_says(&error, &["population of no individual"]);
    }

    #[test]
    fn a_distance_that_holds_no_pair_is_an_error() {
        let error = fit_ld_decay(&[1000, 2000], &[1, 0], &[0.5, 0.0], 100).expect_err("the error");
        assert!(
            matches!(error, Error::LdDecayDistWithNoPair { dist: 2000 }),
            "the error of a distance that holds no pair is {error:?}"
        );
        assert_the_message_says(&error, &["the distance 2000", "no pair"]);
    }

    #[test]
    fn a_sum_of_r2_that_is_not_finite_is_an_error() {
        let error =
            fit_ld_decay(&[1000, 2000], &[1, 1], &[0.5, f64::NAN], 100).expect_err("the error");
        assert!(
            matches!(error, Error::LdDecaySumOfR2OutOfRange { dist: 2000, sum_r2 } if sum_r2.is_nan()),
            "the error of a sum of r² that is NaN is {error:?}"
        );
        assert_the_message_says(&error, &["the distance 2000", "is NaN"]);
        let error = fit_ld_decay(&[1000, 2000], &[1, 1], &[f64::INFINITY, 0.25], 100)
            .expect_err("the error");
        assert!(
            matches!(error, Error::LdDecaySumOfR2OutOfRange { dist: 1000, sum_r2 } if sum_r2.is_infinite()),
            "the error of a sum of r² that is an infinity is {error:?}"
        );
        assert_the_message_says(&error, &["the distance 1000", "is inf"]);
    }

    #[test]
    fn a_sum_of_r2_below_zero_is_an_error() {
        let error =
            fit_ld_decay(&[1000, 2000], &[1, 1], &[0.5, -0.25], 100).expect_err("the error");
        assert!(
            matches!(error, Error::LdDecaySumOfR2OutOfRange { dist: 2000, sum_r2 } if sum_r2 < 0.0),
            "the error of a sum of r² below 0 is {error:?}"
        );
        assert_the_message_says(&error, &["the distance 2000", "is -0.25"]);
    }
}
