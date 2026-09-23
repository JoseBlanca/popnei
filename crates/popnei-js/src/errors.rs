//! The errors of the core crate on their way to a JavaScript `Error`.
//!
//! `impl From<popnei::Error> for JsValue` cannot be written here, because
//! neither type belongs to this crate, so every function of the crate fails
//! with [`JsPopneiError`], which does belong to it, and wasm-bindgen throws
//! it as an `Error`. `?` on a call of the core crate works everywhere, which
//! is what `.claude/skills/coding/pyo3.md` asks of the Python binding crate
//! and what section 11 of `docs/architecture.md` asks here: one place turns
//! an error of the core into what JavaScript throws.

use wasm_bindgen::{JsError, JsValue};

/// What a function of this crate fails with.
///
/// The TypeScript package checks every argument before the call, in
/// `js/popnei/src/arguments.ts`, because a number of JavaScript reaches a
/// whole number of the core as 32 bits with no error, and it throws the
/// `Error` itself. What it cannot check there is whether the memory of the
/// tab takes what is about to be copied into it, which is [`NoMemory`], and
/// what it does not check there is a rule of the core: the threshold of a
/// filter is a number from 0 to 1 by the core's rule, and the package
/// refuses only what is not a number at all, so that the rule lives in one
/// place. What the core says of it crosses as [`Threshold`], which names the
/// argument the user wrote.
///
/// [`NoMemory`]: JsPopneiError::NoMemory
/// [`Threshold`]: JsPopneiError::Threshold
pub enum JsPopneiError {
    /// Something the core crate refused: an argument it takes, or what it
    /// found in the bytes it was given.
    Core(popnei::Error),
    /// Something the core read that JavaScript does not hold: a position
    /// above 2^53, which a float64 rounds.
    NotInJavaScript(String),
    /// An argument this crate refuses before the core sees it, because it
    /// is the crate and not the core that knows the names a user writes:
    /// the name of a statistic and the kind of the bins of a histogram,
    /// each of a finite set of names, whose message writes the set.
    Refused(String),
    /// A threshold of a filter that is not a number from 0 to 1, under the
    /// name of the argument a user wrote it in: the core refuses it and
    /// names the filter by its kind, `maf`, and what a user has to look at
    /// is the call they wrote, `filterByMaf(1.5)`.
    Threshold {
        /// The name of the argument, as a TypeScript user writes it,
        /// `maxAllowedMaf`.
        name: &'static str,
        /// What was given for it, which is NaN, below 0 or above 1.
        threshold: f64,
    },
    /// The memory of wasm does not take what was asked of it: the bytes of
    /// a file that is being given to popnei. A failed allocation aborts in
    /// wasm, and an abort is a trap that leaves the module unusable, so
    /// what can be asked for beforehand is.
    NoMemory(String),
    /// Something that cannot happen unless this crate has a defect: a
    /// chromosome whose number is not in the table of the reader that gave
    /// it, a variant with more alleles than a JavaScript array of counts
    /// holds, populations the pass was not given the name or the
    /// individuals of every one of, a histogram that does not hold one
    /// count for each bin of its distribution, a pass that gave a
    /// different number of names of individuals and of rates, or a pass
    /// with no count of the variants of a pair of the populations it
    /// counted over.
    Broken(String),
}

impl From<popnei::Error> for JsPopneiError {
    fn from(error: popnei::Error) -> JsPopneiError {
        JsPopneiError::Core(error)
    }
}

impl From<JsPopneiError> for JsValue {
    /// The `Error` that JavaScript catches, with the message the error has
    /// in Rust.
    ///
    /// JavaScript has one exception for everything a library refuses, so
    /// the six cases are one `Error`, where Python tells a `ValueError`
    /// from an `OSError`.
    fn from(error: JsPopneiError) -> JsValue {
        let message = match error {
            JsPopneiError::Core(error) => error.to_string(),
            // The threshold of a filter, which is the number a user wrote
            // in the call that adds it: the message names the argument, and
            // the rule it broke is the core's, which refuses the same
            // thresholds when a pass builds its filters.
            JsPopneiError::Threshold { name, threshold } => format!(
                "`{name}` is {value}, and a threshold is a number from 0 to 1, both \
                 included: the number of the variant it is compared with is one count of \
                 the variant divided by another",
                value = as_javascript_writes_it(threshold)
            ),
            JsPopneiError::NotInJavaScript(message)
            | JsPopneiError::Refused(message)
            | JsPopneiError::NoMemory(message)
            | JsPopneiError::Broken(message) => message,
        };
        JsError::new(&message).into()
    }
}

/// `number` written as JavaScript writes it, which is how a user wrote it:
/// `95` and not the `95.0` of Rust, `Infinity` and not its `inf`.
///
/// Rust and JavaScript both write a float64 as the shortest text that reads
/// back as the same number, so the digits are the same, and they differ in
/// the two infinities and in where they turn to an exponent: JavaScript
/// writes 1e21 and larger, and anything below 1e-6, with one, and Rust
/// writes every number in full. No threshold of a filter is in either range,
/// and a number that is refused for being out of 0 to 1 can be: `1e30` is
/// written here as its 31 digits.
fn as_javascript_writes_it(number: f64) -> String {
    if number.is_infinite() {
        return if number.is_sign_negative() {
            "-Infinity".to_owned()
        } else {
            "Infinity".to_owned()
        };
    }
    number.to_string()
}
