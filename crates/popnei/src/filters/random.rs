//! The filter of `docs/specs/filters.md` that keeps variants at random,
//! each with the probability a user gives, so that a calculation runs on a
//! sample spread over the whole of a large file.
//!
//! Whether a variant is kept is drawn from SplitMix64, the generator of
//! random numbers that `java.util.SplittableRandom` of Java uses, started at
//! the user's seed in every pass and drawn once for each variant in the
//! order the variants reach the filter, so that every pass of one
//! `Variants` keeps the same variants.

/// SplitMix64, of Steele, Lea and Flood (2014): a state of 64 bits that
/// starts as the seed, and one draw of 64 bits at each call.
///
/// Each draw adds 0x9E3779B97F4A7C15 to the state and mixes the sum into the
/// number it gives, all of it wrapping at 2^64, as `splitmix64.c` of
/// Sebastiano Vigna and `nextLong` of `java.util.SplittableRandom` do.
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "the filter that keeps variants at random draws from it, and comes in the next commit"
    )
)]
#[derive(Debug, Clone, Copy)]
struct SplitMix64 {
    state: u64,
}

#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "the filter that keeps variants at random draws from it, and comes in the next commit"
    )
)]
impl SplitMix64 {
    /// The generator with its state at `seed`.
    fn new(seed: u64) -> SplitMix64 {
        SplitMix64 { state: seed }
    }

    /// The next draw, a whole number from 0 to 2^64 - 1.
    fn next_draw(&mut self) -> u64 {
        // The wrapping at 2^64 is the arithmetic of the generator, and not
        // an overflow: it is how Java and `splitmix64.c` compute it.
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut mixed = self.state;
        mixed = (mixed ^ (mixed >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        mixed ^ (mixed >> 31)
    }
}

#[cfg(test)]
mod tests;
