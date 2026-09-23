//! The two distributions that turn the statistic of a test into a p-value,
//! which every model of the study ends in.
//!
//! Each of the two is a survival function, the chance that a distribution
//! is beyond a point: [`chi2_sf_1df`] is the chi square with one degree of
//! freedom, and [`t_sf_two_sided`] the Student t, which is written from the
//! regularized incomplete beta function below it. "The two distributions"
//! of `docs/specs/gwas.md` says what each gives and how it is verified.

/// The chance that a chi square with one degree of freedom is above `x`.
///
/// This is the p-value of a score test and of the Wald test of a logistic
/// model: the statistic of either is an effect divided by its standard
/// error and then squared, which is a chi square with one degree of freedom
/// when the variant has no effect on the trait. `x` is such a statistic.
///
/// An `x` of 0 or below gives 1.0, as scipy's `chi2.sf` does: nothing
/// exceeds such a statistic, so the chance of exceeding it is 1. That is an
/// answer and not a refusal, and it is written out because the square root
/// of a negative number is NaN and the line below would otherwise give one.
/// A NaN `x` is a different thing, a statistic that was never computed, and
/// it gives NaN. Which caller may arrive with a statistic below 0 is
/// **Open 2** of `docs/specs/gwas.md`.
///
/// The value is `erfc(sqrt(x / 2))`, the complementary error function of
/// the `libm` crate, which gives how much of a standard normal distribution
/// lies further from 0 than a point. It agrees with scipy 1.18.1's
/// `chi2.sf(x, 1)` within 1e-12 of itself over the sample the test of this
/// module uses, and at `x = 100` it gives 1.5239706e-23, so a variant with
/// a strong effect gets the right tail.
#[must_use]
pub fn chi2_sf_1df(x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    libm::erfc((x / 2.0).sqrt())
}

/// The chance that a Student t with `df` degrees of freedom is further
/// from 0 than `t`, both tails added.
///
/// This is the p-value of the Wald test of a linear model and of a linear
/// mixed model: the statistic of either is an effect divided by its
/// standard error, which follows a Student t when the variant has no
/// effect on the trait, and `df` is how many individuals are left once the
/// columns of the design and the variant have been fitted. `t` of either
/// sign gives the same value, and a NaN gives NaN.
///
/// A `df` of 0 or below gives NaN, and so does a NaN `df`. A model with no
/// individual left over has no answer to give, which is not what a
/// statistic of 0 means, and 0 or fewer degrees of freedom would otherwise
/// take `x` to 0 or below and come back as 0.0, the smallest p-value there
/// is. Nothing calls this with such a `df` today, and what keeps it at 1 or
/// above is
/// [`Error::GwasTooFewIndividuals`](crate::Error::GwasTooFewIndividuals),
/// which refuses a study of no more individuals than the columns of its
/// design plus one: the degrees
/// of freedom of the Wald test of a linear model are the individuals less
/// those columns less one. So this is a guard and not the repair of a live
/// wrong number.
///
/// The value is the regularized incomplete beta function
/// `I_x(df / 2, 1 / 2)` at `x = df / (df + t * t)`, which is written
/// below. `1 - x`, which the function needs as well, is `t * t / (df + t * t)`
/// here and not `1.0 - x`: for a `t` below about `sqrt(df * eps)`, 5e-07 at
/// 197 degrees of freedom, `x` rounds to 1.0 and `1.0 - x` has lost every
/// digit of a difference that is still the whole of the answer. It agrees
/// with scipy 1.18.1's `2 * t.sf(abs(t), df)` within 1e-10 of itself over
/// the fifteen values of `t` and the five degrees of freedom the test of
/// this module uses, from the 1.6e-96 of `t = 40` at 197 degrees of freedom
/// to the 1 less 8e-08 of `t = 1e-07`.
#[must_use]
pub fn t_sf_two_sided(t: f64, df: f64) -> f64 {
    if df <= 0.0 || df.is_nan() {
        return f64::NAN;
    }
    let t_squared = t * t;
    let x = df / (df + t_squared);
    let one_minus_x = t_squared / (df + t_squared);
    regularized_incomplete_beta(df / 2.0, 0.5, x, one_minus_x)
}

/// The regularized incomplete beta function `I_x(a, b)`, the share of the
/// beta distribution with the shapes `a` and `b` that lies below `x`.
///
/// It is 0 for an `x` at or below 0 and 1 for a `one_minus_x` at or below
/// 0. In between it is the continued fraction of `beta_continued_fraction`,
/// multiplied by the front factor and divided by `a`, while `x` is below
/// `(a + 1) / (a + b + 2)`. The front factor is the exponential of
/// `lgamma(a + b) - lgamma(a) - lgamma(b) + a * ln(x) + b * ln(1 - x)`. At
/// or above `(a + 1) / (a + b + 2)` the fraction converges slowly, and the
/// value is taken from the symmetry `I_x(a, b) = 1 - I_{1-x}(b, a)`
/// instead, which is where `1 - x` does most of its work.
///
/// The caller passes `1 - x` as `one_minus_x` rather than leaving the
/// function to take `1.0 - x`, because a caller that can compute it without
/// subtracting has digits that the subtraction cannot recover: `x` carries
/// 16 digits of itself, and when it is near 1 the difference from 1 keeps
/// only what is left of them. `t_sf_two_sided` is such a caller. A caller
/// with nothing better passes `1.0 - x`, and then the spec's `ln_1p(-x)`
/// and this `ln(one_minus_x)` differ by under 1e-17 in the exponent, which
/// is the size the spec gives for that difference.
///
/// Neither `a` nor `b` may be 0 or negative, and no caller of the module
/// passes one: `a` is half the degrees of freedom of a test and `b` is
/// 1 / 2.
fn regularized_incomplete_beta(a: f64, b: f64, x: f64, one_minus_x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if one_minus_x <= 0.0 {
        return 1.0;
    }
    let log_front =
        libm::lgamma(a + b) - libm::lgamma(a) - libm::lgamma(b) + a * x.ln() + b * one_minus_x.ln();
    let front = log_front.exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * beta_continued_fraction(a, b, x) / a
    } else {
        1.0 - front * beta_continued_fraction(b, a, one_minus_x) / b
    }
}

/// The smallest a denominator of the continued fraction may be: one that
/// has come out at 0 is replaced by this before it divides, which is what
/// lets Lentz's method carry on past a term that vanishes instead of
/// giving an infinity.
///
/// No test covers the five lines that read it, and none can, as "The two
/// distributions" of `docs/specs/gwas.md` sets out: with `b` at 1 / 2,
/// which is what every caller here passes, the first denominator is
/// bounded below by `2 / (a + b + 2)` in both branches, so reaching 1e-300
/// needs `a + b` above about 2e300, a panel of 4e300 individuals. The spec
/// measured the smallest `|c|` or `|d|` over 6009003 calls as
/// 4.027585806198886e-6, which is that bound at a million degrees of
/// freedom. It stays because Numerical Recipes and pyNei have it and
/// because a caller with some other `b` would need it. A reader who finds
/// these lines uncovered should not go looking for the argument that
/// reaches them: there is none.
const CONTINUED_FRACTION_TINY: f64 = 1e-300;

/// How near 1 a round's factor has to come for the fraction to stop. Both
/// this and the 500 rounds below are Numerical Recipes', through pyNei. It
/// caps the work and does not get the digits: with it at 0, so that every
/// call runs its 500 rounds, the spec's sweep moved by 2.3e-13 relative at
/// worst and nothing became not finite.
const CONTINUED_FRACTION_EPS: f64 = 1e-15;

/// How many rounds the fraction may take before it gives what it has.
/// Running out is not an error, and no case of the module's tests comes
/// near it: the 115 calls of the two tests below all stopped, the slowest
/// after 47 rounds, measured on 23 September 2026.
const CONTINUED_FRACTION_MAX_ROUNDS: u32 = 500;

/// The continued fraction of the incomplete beta with the shapes `a` and
/// `b` at `x`, the one of Numerical Recipes evaluated by Lentz's method,
/// which builds the fraction from its front rather than from its far end
/// and so stops as soon as a round no longer changes the value.
fn beta_continued_fraction(a: f64, b: f64, x: f64) -> f64 {
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < CONTINUED_FRACTION_TINY {
        d = CONTINUED_FRACTION_TINY;
    }
    d = 1.0 / d;
    let mut h = d;
    for round in 1..=CONTINUED_FRACTION_MAX_ROUNDS {
        let m = f64::from(round);
        let m2 = 2.0 * m;
        let even = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + even * d;
        if d.abs() < CONTINUED_FRACTION_TINY {
            d = CONTINUED_FRACTION_TINY;
        }
        c = 1.0 + even / c;
        if c.abs() < CONTINUED_FRACTION_TINY {
            c = CONTINUED_FRACTION_TINY;
        }
        d = 1.0 / d;
        h *= d * c;
        let odd = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + odd * d;
        if d.abs() < CONTINUED_FRACTION_TINY {
            d = CONTINUED_FRACTION_TINY;
        }
        c = 1.0 + odd / c;
        if c.abs() < CONTINUED_FRACTION_TINY {
            c = CONTINUED_FRACTION_TINY;
        }
        d = 1.0 / d;
        let delta = d * c;
        h *= delta;
        if (delta - 1.0).abs() < CONTINUED_FRACTION_EPS {
            break;
        }
    }
    h
}

/// The two distributions against scipy 1.18.1, whose numbers are the
/// literals here, as "How it is verified" of "The two distributions" of
/// `docs/specs/gwas.md` asks.
#[cfg(test)]
mod tests {
    use super::{chi2_sf_1df, regularized_incomplete_beta, t_sf_two_sided};

    /// `x` and `scipy.stats.chi2.sf(x, 1)` of scipy 1.18.1, each printed
    /// by Python's `repr`, which gives the fewest digits that read the same
    /// `f64` back.
    ///
    /// The thirteen values of `x` that come first are a sample of a chi
    /// square with one degree of freedom: the 1000 draws of
    /// `numpy.random.default_rng(0).chisquare(1, 1000)` of numpy 2.5.3,
    /// sorted, at the ranks 0, 90, 180 and so on to 990, and the largest of
    /// the 1000. They run from 9.8e-08, where the p-value is near 1, to
    /// 10.8. The
    /// three that follow are the 30, 50 and 100 the spec asks for, which
    /// reach a p-value of 1.5e-23 and are the tail a variant with a strong
    /// effect lands in.
    const SCIPY_CHI2_SF_1DF: [(f64, f64); 16] = [
        (9.83015230925879e-08, 0.999749838669647),
        (0.012365365489339817, 0.911458017301763),
        (0.05460850890382662, 0.8152298165740756),
        (0.11322041823323797, 0.7365071071670479),
        (0.20872737725861593, 0.6477667063269484),
        (0.3141949722019434, 0.5751173190337298),
        (0.4625430778164763, 0.49643785458596523),
        (0.6940024026584967, 0.404806342457091),
        (1.054858588171547, 0.3043907579517296),
        (1.605831694008574, 0.20507872868262017),
        (2.5243520191320643, 0.11210082147194526),
        (6.274619438315357, 0.012247942227758789),
        (10.821607681554072, 0.0010032233803823583),
        (30.0, 4.320463057827496e-08),
        (50.0, 1.537459794428033e-12),
        (100.0, 1.5239706048320995e-23),
    ];

    /// The tolerance is the 1e-12 relative of the spec, and not the bits:
    /// `erfc` is one of the functions that are not rounded the same on
    /// every platform, so popnei promises the digits and not the last one.
    /// The largest difference over the sixteen values, on this Mac on 23
    /// September 2026, was 1.8e-14 of scipy's value, at `x = 1.6058`, so
    /// the check has 57 times the headroom it needs and a platform that
    /// rounds `erfc` elsewhere still passes it.
    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "an x at or below 0 returns the literal 1.0, which no \
                  arithmetic has touched, so the two are equal to the bit"
    )]
    fn chi2_sf_1df_matches_scipy_over_a_chi_square_sample_and_in_the_tail() {
        for (x, scipy) in SCIPY_CHI2_SF_1DF {
            let popnei = chi2_sf_1df(x);
            let relative = ((popnei - scipy) / scipy).abs();
            assert!(
                relative < 1e-12,
                "chi2_sf_1df({x}) gave {popnei}, scipy gives {scipy}, \
                 which differ by {relative} of scipy's value"
            );
        }
        for x in [-1.0, -1e-9, 0.0] {
            assert_eq!(
                chi2_sf_1df(x),
                1.0,
                "nothing exceeds a statistic of {x}, so the chance of \
                 exceeding it is 1, which is what scipy's chi2.sf gives"
            );
        }
        assert!(
            chi2_sf_1df(f64::NAN).is_nan(),
            "a statistic that is not a number has no p-value"
        );
    }

    /// The four pairs `(a, b)` the spec asks the incomplete beta to be
    /// checked at, in the order the rows of `SCIPY_BETAINC` have them.
    /// `(98.5, 0.5)` is what a Student t with 197 degrees of freedom uses,
    /// next to the panel's 196: 200 individuals less the three columns of
    /// its design less one for the variant.
    const BETA_PAIRS: [(f64, f64); 4] = [(0.5, 0.5), (10.0, 0.5), (98.5, 0.5), (2.5, 7.0)];

    /// The `x` the incomplete beta is taken at: the 1000 draws of
    /// `numpy.random.default_rng(0).uniform(0, 1, 1000)` of numpy 2.5.3,
    /// sorted, at the ranks 0, 111, 222 and so on to 999, so that the last
    /// of them is the largest of the 1000.
    ///
    /// The function has two branches, the continued fraction below
    /// `(a + 1) / (a + b + 2)` and the symmetry `I_x(a, b) =
    /// 1 - I_{1-x}(b, a)` at or above it, and every pair of `BETA_PAIRS`
    /// takes both over these ten values. For `(98.5, 0.5)` the turn is at
    /// 0.985 and the largest of the 1000, 0.9995, is the only one past it.
    const BETA_X: [f64; 10] = [
        0.00019000160734350402,
        0.12681710226124776,
        0.24209440214479672,
        0.36730140225246877,
        0.4752246242920236,
        0.5776878925178592,
        0.6848100430995354,
        0.7892475482526763,
        0.8889806580193464,
        0.9995013522570269,
    ];

    /// `scipy.special.betainc(a, b, x)` of scipy 1.18.1, one row per pair of
    /// `BETA_PAIRS`, each row over `BETA_X`. The one value that is 0.0 and
    /// the one that is 1.0 are what the front factor of the function
    /// underflows to: at `a = 98.5` and `x = 0.00019` it is `exp(-844)`.
    const SCIPY_BETAINC: [[f64; 10]; 4] = [
        // a = 0.5, b = 0.5
        [
            0.008775513005121197,
            0.2317969726600946,
            0.3274906007616421,
            0.41449691362201163,
            0.4842210445079407,
            0.5496588399779172,
            0.6205121607128349,
            0.6963598067259873,
            0.783746069762884,
            0.9857828301556638,
        ],
        // a = 10.0, b = 0.5
        [
            1.080460509271319e-38,
            2.0155738406315239e-10,
            1.3801539354259724e-07,
            9.655661123183734e-06,
            0.0001375787217554091,
            0.0010620947644470852,
            0.0065538478356862995,
            0.03161946390607385,
            0.1297005104596951,
            0.9214264864085051,
        ],
        // a = 98.5, b = 0.5
        [
            0.0,
            2.794490395778511e-90,
            1.36820610510675e-62,
            1.0163053160450192e-44,
            1.1661642215623547e-33,
            2.9187000208624543e-25,
            6.366584355165142e-18,
            9.119063601149856e-12,
            1.5182519264715282e-06,
            0.7542320646178294,
        ],
        // a = 2.5, b = 7.0
        [
            2.4788133598912145e-08,
            0.16372356028276183,
            0.4889434870152894,
            0.7825611783294392,
            0.922059784639868,
            0.978602728785751,
            0.9966059404816893,
            0.9997577325040361,
            0.9999968246187708,
            1.0,
        ],
    ];

    /// The tolerance is the 1e-12 absolute of the spec and, beside it, 1e-12
    /// of scipy's value wherever that value is above 0. The absolute one
    /// alone checks nothing at 8 of the 40 literals, which are at or below
    /// 1e-12 themselves, seven of them in the row of `(98.5, 0.5)`: an
    /// implementation that gave 0 in the tail passed it. The largest
    /// absolute difference at each pair, on this Mac on 23 September 2026,
    /// was 1.3e-15 at `(0.5, 0.5)`, 2.7e-15 at `(10, 0.5)`, 8.5e-15 at
    /// `(98.5, 0.5)` and 3.3e-16 at `(2.5, 7)`, and the largest relative
    /// one was 1.4e-15, 2.1e-14, 6.6e-14 and 1.2e-15 at the same four, so
    /// the relative bound has 15 times the room it needs where it has
    /// least, at `(98.5, 0.5)`, and a platform that rounds `exp` and
    /// `lgamma` elsewhere still passes both.
    ///
    /// The `x` at or below 0 and at or above 1, which a Student t never
    /// reaches because its `x` is `df / (df + t * t)`, are asserted here
    /// too, at the 0 and the 1 the spec gives them, so that the section's
    /// count of tests is the three it names.
    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "an x outside (0, 1) returns the literal 0.0 or 1.0, which \
                  no arithmetic has touched, so the two are equal to the bit"
    )]
    fn incomplete_beta_matches_scipy_at_the_four_pairs_the_two_models_use() {
        for ((a, b), scipy_row) in BETA_PAIRS.into_iter().zip(SCIPY_BETAINC) {
            for (x, scipy) in BETA_X.into_iter().zip(scipy_row) {
                let popnei = regularized_incomplete_beta(a, b, x, 1.0 - x);
                let difference = (popnei - scipy).abs();
                assert!(
                    difference < 1e-12,
                    "the incomplete beta of a = {a}, b = {b} at x = {x} gave \
                     {popnei}, scipy gives {scipy}, which differ by {difference}"
                );
                if scipy > 0.0 {
                    let relative = difference / scipy;
                    assert!(
                        relative < 1e-12,
                        "the incomplete beta of a = {a}, b = {b} at x = {x} \
                         gave {popnei}, scipy gives {scipy}, which differ by \
                         {relative} of scipy's value"
                    );
                }
            }
            for x in [-1.0, 0.0] {
                assert_eq!(regularized_incomplete_beta(a, b, x, 1.0 - x), 0.0);
            }
            for x in [1.0, 2.0] {
                assert_eq!(regularized_incomplete_beta(a, b, x, 1.0 - x), 1.0);
            }
        }
    }

    /// The degrees of freedom of `SCIPY_T_SF_TWO_SIDED`, one per row, the
    /// five "How it is verified" of the spec names. 5 and 17 are the small
    /// ones pyNei checks, where the t is furthest from a normal. 197 is
    /// what the pair `(98.5, 0.5)` of `BETA_PAIRS` comes from, beside the
    /// panel's own 198 with the intercept alone and 196 with its two
    /// covariates, which is 200 individuals less the columns of the design
    /// less one for the variant. 997 and 9997 are the panels of 1000 and
    /// of 10000 individuals of `docs/objectives.md` counted the same way,
    /// and they are here because the error of the function grows with the
    /// degrees of freedom and the spec claims its 1e-10 up to 9997 and no
    /// further.
    const T_DEGREES_OF_FREEDOM: [f64; 5] = [5.0, 17.0, 197.0, 997.0, 9997.0];

    /// The `t` the two sided tail is taken at: the 1000 draws of
    /// `numpy.random.default_rng(0).standard_normal(1000) * 3` of numpy
    /// 2.5.3, sorted, at the ranks 0, 111, 222 and so on to 999, and then
    /// the 1e-07 and 1e-05 of the paragraph below and the 10, 20 and 40 the
    /// spec asks for. Five of the fifteen are negative, which the function
    /// takes through `t * t`, and the three largest are the tail a variant
    /// with a strong effect lands in: at 197 degrees of freedom `t = 40`
    /// has a p-value of 1.6e-96.
    ///
    /// A `t` of 1e-07 and one of 1e-05 are there because `x` of the
    /// incomplete beta is `df / (df + t * t)`, which rounds to 1.0 for a
    /// `t` below about `sqrt(df * eps)`, 5e-07 at 197 degrees of freedom
    /// and 3e-06 at 9997. Their p-value is 1 less something near 8e-08 and
    /// 8e-06, so an implementation that reads `1 - x` back from a rounded
    /// `x` loses the whole of that difference and gives 1.0.
    const T_VALUES: [f64; 15] = [
        -11.698265190163017,
        -3.6455770776863603,
        -2.2826132625192868,
        -1.4177630302754523,
        -0.6135674651989825,
        0.2141653057738903,
        1.0758242029481555,
        2.156821697445163,
        3.4945919433330843,
        9.19811021714669,
        1e-07,
        1e-05,
        10.0,
        20.0,
        40.0,
    ];

    /// `2 * scipy.stats.t.sf(abs(t), df)` of scipy 1.18.1, one row per
    /// degrees of freedom of `T_DEGREES_OF_FREEDOM`, each row over
    /// `T_VALUES`. The last of the row of 9997 is 0.0 because the value
    /// is below the smallest number an `f64` holds.
    const SCIPY_T_SF_TWO_SIDED: [[f64; 15]; 5] = [
        // df = 5
        [
            8.021986064011076e-05,
            0.014816691179735638,
            0.07130468304194079,
            0.21545726847462635,
            0.5663371298234476,
            0.8388781542853476,
            0.3311547120148151,
            0.08352172539200577,
            0.017384749708466502,
            0.00025487505551227214,
            0.999999924078662,
            0.9999924078662037,
            0.00017094757574296357,
            5.7755163732241715e-06,
            1.841196217177295e-07,
        ],
        // df = 17
        [
            1.484092308678517e-09,
            0.0020008264729398156,
            0.03559999374852078,
            0.17432936212171152,
            0.54762394337886,
            0.8329659145993707,
            0.29703772549565155,
            0.045628366419844114,
            0.0027764075279970296,
            5.1955514732118756e-08,
            0.999999921375654,
            0.9999921375654,
            1.5482821011627456e-08,
            2.9927008195203556e-13,
            2.9320117368673818e-18,
        ],
        // df = 197
        [
            2.398338008739183e-24,
            0.0003413422609230366,
            0.023521078785496523,
            0.15783973504362406,
            0.5402089878596277,
            0.8306396600846631,
            0.2833216667123333,
            0.032230891954663506,
            0.00058638199104543,
            5.316767147362756e-17,
            0.9999999203127337,
            0.9999920312733652,
            2.65918277307826e-19,
            2.5792670821731053e-49,
            1.6256956659911944e-96,
        ],
        // df = 997
        [
            1.0285108008842e-29,
            0.0002805744116466835,
            0.022663015344545432,
            0.15657241277885833,
            0.5396411434338867,
            0.8304619422914968,
            0.2822661515237383,
            0.03125803706343817,
            0.0004956644692266718,
            2.0753391089503426e-19,
            0.9999999202315485,
            0.9999920231548545,
            1.678371723944798e-22,
            4.386345719949788e-75,
            1.7484717633928704e-209,
        ],
        // df = 9997
        [
            2.083192798083594e-31,
            0.00026814639056538145,
            0.022474062079316716,
            0.15629112844735527,
            0.5395151482437246,
            0.8304225182279459,
            0.28203184818277116,
            0.031043323571632118,
            0.00047684358853721274,
            4.3710342345251525e-20,
            0.9999999202135392,
            0.9999920213539206,
            1.9634289994060445e-23,
            2.7678195351970412e-87,
            0.0,
        ],
    ];

    /// The tolerance is the 1e-10 relative of the spec, which holds down to
    /// the 1.6e-96 of `t = 40` at 197 degrees of freedom: a p-value is read
    /// at every size, so this one is relative where the incomplete beta's
    /// is absolute. The largest difference at each degrees of freedom, on
    /// this Mac on 23 September 2026, was 1.8e-15 at 5, 8.3e-15 at 17,
    /// 1.9e-13 at 197, 1.6e-13 at 997 and 2.2e-11 at 9997, all of scipy's
    /// value, so the bound has 4.5 times the room it needs at 9997, where
    /// it has least. The error grows with the degrees of freedom and the
    /// spec claims the 1e-10 no further than 9997.
    ///
    /// A value scipy gives as 0.0 is below the smallest number an `f64`
    /// holds, and there a ratio says nothing; the assertion there is that
    /// popnei underflowed too.
    #[test]
    fn t_sf_two_sided_matches_scipy_from_5_to_9997_degrees_of_freedom() {
        for (df, scipy_row) in T_DEGREES_OF_FREEDOM.into_iter().zip(SCIPY_T_SF_TWO_SIDED) {
            for (t, scipy) in T_VALUES.into_iter().zip(scipy_row) {
                let popnei = t_sf_two_sided(t, df);
                if scipy <= 0.0 {
                    assert!(
                        popnei < f64::MIN_POSITIVE,
                        "t_sf_two_sided({t}, {df}) gave {popnei}, where scipy \
                         underflowed to 0"
                    );
                    continue;
                }
                let relative = ((popnei - scipy) / scipy).abs();
                assert!(
                    relative < 1e-10,
                    "t_sf_two_sided({t}, {df}) gave {popnei}, scipy gives \
                     {scipy}, which differ by {relative} of scipy's value"
                );
            }
        }
        assert!(
            t_sf_two_sided(f64::NAN, 197.0).is_nan(),
            "a statistic that is not a number has no p-value"
        );
        for df in [0.0, -4.0, f64::NAN] {
            assert!(
                t_sf_two_sided(1.0, df).is_nan(),
                "{df} degrees of freedom has no p-value, and before this was \
                 asserted `t_sf_two_sided(1.0, 0.0)` gave 0.0"
            );
        }
    }
}
