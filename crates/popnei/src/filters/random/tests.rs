//! The tests of the filter that keeps variants at random, each with
//! `random_filter` in its name, against "How it is verified" of that filter
//! in `docs/specs/filters.md`.

use super::SplitMix64;

/// The generator from a seed of 1234567 gives the first five draws of
/// `nextLong` of `java.util.SplittableRandom`, which `java
/// tests/reference/filters/SplitMix.java` printed with OpenJDK 26.0.2.1 on 5
/// October 2026, and which `splitmix64.c` of Sebastiano Vigna gives too.
#[test]
fn random_filter_generator_gives_the_five_draws_of_java_from_1234567() {
    let mut generator = SplitMix64::new(1_234_567);
    let draws: Vec<u64> = (0..5).map(|_| generator.next_draw()).collect();
    assert_eq!(
        draws,
        [
            6_457_827_717_110_365_317,
            3_203_168_211_198_807_973,
            9_817_491_932_198_370_423,
            4_593_380_528_125_082_431,
            16_408_922_859_458_223_821,
        ]
    );
}
