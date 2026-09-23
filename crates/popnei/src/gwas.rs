//! The association study of a dataset: which of its variants go with a
//! trait of the individuals, with the effect of each variant on the trait,
//! the uncertainty of that effect, and the chance of an effect at least
//! that large in a dataset where the variant has none, which is its
//! p-value. `docs/specs/gwas.md` says what a study gives and how each of
//! its models is verified.
//!
//! What is here is the first of the two functions that turn the statistic
//! of a test into a p-value, which every model of the module ends in: the
//! chi square with one degree of freedom of a score test. The other, the
//! Student t of the Wald test of a linear model, and the models that give
//! the statistics to both, are being written.

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

/// The two distributions against scipy 1.18.1, whose numbers are the
/// literals here, as "How it is verified" of "The two distributions" of
/// `docs/specs/gwas.md` asks.
#[cfg(test)]
mod distributions {
    use super::chi2_sf_1df;

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
}
