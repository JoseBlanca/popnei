//! The association study of a dataset: which of its variants go with a
//! trait of the individuals, with the effect of each variant on the trait,
//! the uncertainty of that effect, and the chance of an effect at least
//! that large in a dataset where the variant has none, which is its
//! p-value. `docs/specs/gwas.md` says what a study gives and how each of
//! its models is verified.
//!
//! What is here is the two functions that turn the statistic of a test
//! into a p-value, which every model of the module ends in: the chi square
//! with one degree of freedom of a score test, and the Student t of the
//! Wald test of a linear model, which is written from the regularized
//! incomplete beta function below it. The models that give the statistics
//! to both are being written.

/// The chance that a chi square with one degree of freedom is above `x`.
///
/// This is the p-value of a score test and of the Wald test of a logistic
/// model: the statistic of either is an effect divided by its standard
/// error and then squared, which is a chi square with one degree of freedom
/// when the variant has no effect on the trait. `x` is such a statistic and
/// is never negative; a negative `x` gives NaN, and so does a NaN.
///
/// The value is `erfc(sqrt(x / 2))`, the complementary error function of
/// the `libm` crate, which gives how much of a standard normal distribution
/// lies further from 0 than a point. It agrees with scipy 1.18.1's
/// `chi2.sf(x, 1)` within 1e-12 of itself over the sample the test of this
/// module uses, and at `x = 100` it gives 1.5239706e-23, so a variant with
/// a strong effect gets the right tail.
#[must_use]
pub fn chi2_sf_1df(x: f64) -> f64 {
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
/// The value is the regularized incomplete beta function
/// `I_x(df / 2, 1 / 2)` at `x = df / (df + t * t)`, which is written
/// below. It agrees with scipy 1.18.1's `2 * t.sf(abs(t), df)` within
/// 1e-10 of itself over the thirteen values of `t` and the three degrees
/// of freedom the test of this module uses, down to the 1.6e-96 of
/// `t = 40` at 197 degrees of freedom.
#[must_use]
pub fn t_sf_two_sided(t: f64, df: f64) -> f64 {
    regularized_incomplete_beta(df / 2.0, 0.5, df / (df + t * t))
}

/// The regularized incomplete beta function `I_x(a, b)`, the share of the
/// beta distribution with the shapes `a` and `b` that lies below `x`.
///
/// It is 0 for an `x` at or below 0 and 1 for an `x` at or above 1. In
/// between it is the continued fraction of `beta_continued_fraction`,
/// multiplied by the front factor and divided by `a`, while `x` is below
/// `(a + 1) / (a + b + 2)`. The front factor is the exponential of
/// `lgamma(a + b) - lgamma(a) - lgamma(b) + a * ln(x) + b * ln(1 - x)`,
/// where `ln(1 - x)` is computed as `ln_1p(-x)`, as pyNei computes it,
/// which keeps its digits for an `x` near 0. At or above
/// `(a + 1) / (a + b + 2)` the fraction converges slowly, and the value is
/// taken from the symmetry `I_x(a, b) = 1 - I_{1-x}(b, a)` instead.
///
/// Neither `a` nor `b` may be 0 or negative, and no caller of the module
/// passes one: `a` is half the degrees of freedom of a test and `b` is
/// 1 / 2.
fn regularized_incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let log_front =
        libm::lgamma(a + b) - libm::lgamma(a) - libm::lgamma(b) + a * x.ln() + b * (-x).ln_1p();
    let front = log_front.exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * beta_continued_fraction(a, b, x) / a
    } else {
        1.0 - front * beta_continued_fraction(b, a, 1.0 - x) / b
    }
}

/// The smallest a denominator of the continued fraction may be: one that
/// has come out at 0 is replaced by this before it divides, which is what
/// lets Lentz's method carry on past a term that vanishes instead of
/// giving an infinity.
const CONTINUED_FRACTION_TINY: f64 = 1e-300;

/// How near 1 a round's factor has to come for the fraction to stop. Both
/// this and the 500 rounds below are Numerical Recipes', through pyNei.
const CONTINUED_FRACTION_EPS: f64 = 1e-15;

/// How many rounds the fraction may take before it gives what it has.
/// Running out is not an error, and no case of the module's tests comes
/// near it: the 79 calls of the two tests below all stopped, the slowest
/// after 30 rounds, measured on 23 September 2026.
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
mod distributions {
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
    /// `BETA_PAIRS`, each row over `BETA_X`. The two values that are 0.0 and
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

    /// The tolerance is the 1e-12 absolute of the spec, which is the one
    /// pyNei uses for the same four pairs, and not a relative one: the
    /// values run down to 2.8e-90, and what a p-value is read to is the
    /// value itself and not the last digit of a number that small.
    ///
    /// The largest difference at each pair, on this Mac on 23 September
    /// 2026, was 1.3e-15 at `(0.5, 0.5)`, 2.7e-15 at `(10, 0.5)`, 8.5e-15
    /// at `(98.5, 0.5)` and 3.3e-16 at `(2.5, 7)`, so the check has 117
    /// times the headroom it needs where it has least, and a platform that
    /// rounds `exp` and `lgamma` elsewhere still passes it.
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
                let popnei = regularized_incomplete_beta(a, b, x);
                let difference = (popnei - scipy).abs();
                assert!(
                    difference < 1e-12,
                    "the incomplete beta of a = {a}, b = {b} at x = {x} gave \
                     {popnei}, scipy gives {scipy}, which differ by {difference}"
                );
            }
            for x in [-1.0, 0.0] {
                assert_eq!(regularized_incomplete_beta(a, b, x), 0.0);
            }
            for x in [1.0, 2.0] {
                assert_eq!(regularized_incomplete_beta(a, b, x), 1.0);
            }
        }
    }

    /// The degrees of freedom of `SCIPY_T_SF_TWO_SIDED`, one per row. A
    /// study of the panel has 197 of them with no covariate and 195 with
    /// the two it has; 5 and 17 are the small ones pyNei checks, where the
    /// t is furthest from a normal.
    const T_DEGREES_OF_FREEDOM: [f64; 3] = [5.0, 17.0, 197.0];

    /// The `t` the two sided tail is taken at: the 1000 draws of
    /// `numpy.random.default_rng(0).standard_normal(1000) * 3` of numpy
    /// 2.5.3, sorted, at the ranks 0, 111, 222 and so on to 999, and then
    /// the 10, 20 and 40 the spec asks for. Half of them are negative,
    /// which the function takes through `t * t`, and the last three are the
    /// tail a variant with a strong effect lands in: at 197 degrees of
    /// freedom `t = 40` has a p-value of 1.6e-96.
    const T_VALUES: [f64; 13] = [
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
        10.0,
        20.0,
        40.0,
    ];

    /// `2 * scipy.stats.t.sf(abs(t), df)` of scipy 1.18.1, one row per
    /// degrees of freedom of `T_DEGREES_OF_FREEDOM`, each row over
    /// `T_VALUES`.
    const SCIPY_T_SF_TWO_SIDED: [[f64; 13]; 3] = [
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
            2.65918277307826e-19,
            2.5792670821731053e-49,
            1.6256956659911944e-96,
        ],
    ];

    /// The tolerance is the 1e-10 relative of the spec, which holds down to
    /// the 1.6e-96 of `t = 40` at 197 degrees of freedom: a p-value is read
    /// at every size, so this one is relative where the incomplete beta's
    /// is absolute. The largest difference over the thirty-nine values, on
    /// this Mac on 23 September 2026, was 2.0e-13 of scipy's value, at
    /// `t = -1.4178` with 197 degrees of freedom, which is 512 times the
    /// headroom the check needs.
    #[test]
    fn t_sf_two_sided_matches_scipy_at_5_17_and_197_degrees_of_freedom() {
        for (df, scipy_row) in T_DEGREES_OF_FREEDOM.into_iter().zip(SCIPY_T_SF_TWO_SIDED) {
            for (t, scipy) in T_VALUES.into_iter().zip(scipy_row) {
                let popnei = t_sf_two_sided(t, df);
                let relative = ((popnei - scipy) / scipy).abs();
                assert!(
                    relative < 1e-10,
                    "t_sf_two_sided({t}, {df}) gave {popnei}, scipy gives \
                     {scipy}, which differ by {relative} of scipy's value"
                );
            }
        }
    }
}
