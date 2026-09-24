//! What a function of this crate fails with, on its way to what JavaScript
//! catches: an error of popnei, the core crate's or this crate's own, which
//! crosses as an `Error`; or the value an application threw from the
//! function that is told how far a pass has got, which stopped the run and
//! crosses back as it is. The second is no error of popnei and need not be
//! an `Error` at all, since it is whatever the application threw, and it is
//! [`Stopped`] here.
//!
//! `impl From<popnei::Error> for JsValue` cannot be written here, because
//! neither type belongs to this crate, so every function of the crate fails
//! with [`JsPopneiError`], which does belong to it, and what JavaScript
//! throws is what this module turns that into. `?` on a call of the core
//! crate works everywhere, which is what `.claude/skills/coding/pyo3.md`
//! asks of the Python binding crate and what section 11 of
//! `docs/architecture.md` asks here: one place turns an error of the core
//! into what JavaScript throws.
//!
//! [`Stopped`]: JsPopneiError::Stopped

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
    /// A pass that gave the matrix of every pair more variants than it was
    /// allowed to take, under the name of the argument a user wrote that
    /// number in: the core refuses it and names `max_num_vars`, which is
    /// the argument of the Python package, and what a TypeScript user has
    /// to look at is the `maxNumVars` of the call they wrote.
    TooManyVars {
        /// How many variants the pass had given when it was stopped, which
        /// is the first count above the cap.
        num_vars: usize,
        /// How many variants the calculation was allowed to take, the
        /// `maxNumVars` of the call.
        max_num_vars: usize,
        /// How many bytes the matrix of those variants holds, 8 for each
        /// pair.
        bytes: u64,
    },
    /// Two individuals of a kinship that have no variant called in both of
    /// them, under their names: the core refuses the pair and names the
    /// positions the two have among the individuals the kinship was asked
    /// for, and what a user has to drop from the panel is a name. With
    /// `individuals` on the call those positions are not even the file's.
    ///
    /// `one` and `other` are the same name when that individual has no
    /// called genotype at all among the variants that were used, which is a
    /// sequencing that failed, and the message then says to leave that one
    /// out instead of one of a pair, as the core's does.
    PairWithNoVariantCalled {
        /// The name of the first of the two.
        one: String,
        /// The name of the second, which is `one` when the individual has
        /// no called genotype at all.
        other: String,
        /// How many of the variants that were used are called in the first.
        num_vars_of_one: u64,
        /// How many of them are called in the second.
        num_vars_of_other: u64,
    },
    /// The memory of wasm does not take what was asked of it: the bytes of
    /// a file that is being given to popnei. A failed allocation aborts in
    /// wasm, and an abort is a trap that leaves the module unusable, so
    /// what can be asked for beforehand is.
    NoMemory(String),
    /// The value the function that is told how far a pass has got threw,
    /// which ended the run and crosses back to the application untouched.
    ///
    /// It is the one case that is not an error of popnei. A run is stopped
    /// by that function throwing: the read it threw in fails, the error
    /// travels out through the readers of the core, and the consumer puts
    /// this case in its place, whatever error the core gave, so that nothing
    /// depends on which reader turned the failed read into which error. What
    /// the application catches is the value it threw itself, which it
    /// recognises with `===` and without reading a message, as "What the
    /// source tells the page" of `docs/specs/js_sources.md` says.
    Stopped(JsValue),
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
    /// in Rust, or the value of the application that stopped a run, which
    /// crosses as it is.
    ///
    /// JavaScript has one exception for everything a library refuses, so
    /// the eight cases that are an error of popnei are one `Error`, where
    /// Python tells a `ValueError` from an `OSError`. The ninth,
    /// [`Stopped`], is not an error of popnei: what it holds is the value
    /// the application threw, and it goes back as it came.
    ///
    /// [`Stopped`]: JsPopneiError::Stopped
    fn from(error: JsPopneiError) -> JsValue {
        let message = match error {
            // The value the application threw, given back without being
            // made into anything: it is what an application tells its own
            // cancel by, and an `Error` around it would be a message to
            // read.
            JsPopneiError::Stopped(thrown) => return thrown,
            JsPopneiError::Core(error) => the_message_of_the_core(&error),
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
            // The cap on the variants of the matrix of every pair, which is
            // the core's message with the name of the argument a TypeScript
            // user wrote: the core says `max_num_vars`, which is what a
            // Python user reads and what no call of TypeScript has.
            JsPopneiError::TooManyVars {
                num_vars,
                max_num_vars,
                bytes,
            } => format!(
                "the pass gave {num_vars} variants and `maxNumVars` is {max_num_vars}: \
                 the matrix of {num_vars} variants holds one r² for each pair of them, \
                 {bytes} bytes of 8 each, and the pass was stopped as soon as it passed \
                 that number, so its source may hold more variants; raise `maxNumVars` \
                 or filter the variants"
            ),
            // The two individuals with no variant called in both, which the
            // core names by their positions among the individuals of the
            // kinship: what a user drops from the panel is a name, and the
            // binding is what holds them.
            JsPopneiError::PairWithNoVariantCalled {
                one,
                other,
                num_vars_of_one,
                num_vars_of_other,
            } => a_pair_with_no_variant_called(&one, &other, num_vars_of_one, num_vars_of_other),
            JsPopneiError::NotInJavaScript(message)
            | JsPopneiError::Refused(message)
            | JsPopneiError::NoMemory(message)
            | JsPopneiError::Broken(message) => message,
        };
        JsError::new(&message).into()
    }
}

/// The message an error of the core has, with the argument it names written
/// as a TypeScript user writes it.
///
/// The owner decided on 24 September 2026 that in TypeScript a message names
/// the option in TypeScript style. The core writes an argument the way Rust
/// and Python spell it, `transform_to_biallelic`, and a user of this package
/// wrote `transformToBiallelic`: what they grep their code for has to be a
/// name their code holds. The core keeps its spelling, since a Python user
/// reads the same sentence and writes the same name, and the rewrite is
/// here, where section 11 of `docs/architecture.md` puts the one place that
/// turns an error of the core into what JavaScript throws.
///
/// It is matched on the error and not on the text: a name is rewritten only
/// in the message that is known to be about that option, so a `num_bins`
/// that is a column of a file somewhere else is left alone. The two names an
/// error carries that this does not reach are `max_num_vars`, which
/// [`JsPopneiError::TooManyVars`] writes as `maxNumVars`, and
/// `poly_threshold` and `bin_type`, which `stats.rs` writes as
/// `polyThreshold` and `binType` before the error gets here.
///
/// Three of the eight names the core writes are left as they are.
/// `num_prin_comps` is in the error of a second pass that was not made,
/// which `pca.rs` of this crate opens a reader for whenever the weights are
/// asked for, so no call of TypeScript reaches it. The `max_num_vars` of a
/// cap the machine does not count the pairs of is refused by the package at
/// the call, against the same number the core checks. And `popnei_batches`
/// is not an option a user writes but the key popnei puts in the footer of a
/// vars file, which is spelled that way in the file itself.
fn the_message_of_the_core(error: &popnei::Error) -> String {
    let message = error.to_string();
    // `matches!` and not a `match` over the enum: the errors of the core are
    // many and the ones that name an option are these, so a new one falls
    // through to its own message rather than to an arm that guesses.
    if matches!(error, popnei::Error::VariantWithMoreThanTwoAlleles { .. }) {
        return message.replace("transform_to_biallelic", "transformToBiallelic");
    }
    if matches!(error, popnei::Error::HistWithNoBin) {
        return message.replace("num_bins", "numBins");
    }
    if matches!(
        error,
        popnei::Error::BlockTooLarge { .. } | popnei::Error::VarsTextTooLarge { .. }
    ) {
        return message.replace("num_vars_per_block", "numVarsPerBlock");
    }
    message
}

/// What [`JsPopneiError::PairWithNoVariantCalled`] says: the two individuals
/// that have no variant called in both of them, or the one individual that
/// has no called genotype at all among the variants that were used.
///
/// It is the message of `Error::KinshipPairWithNoVariantCalled` of the core
/// with the names of the two where the core writes their positions. The
/// entry of a pair is divided by how many variants both of its individuals
/// were called at, and both cases are that number being 0: a pair reaches it
/// when each of the two was called somewhere and never together, and one
/// individual reaches it against itself when its sequencing failed, and then
/// every pair it is in has no variant either, so what a user has to do is
/// leave that one out and not one of a pair.
fn a_pair_with_no_variant_called(
    one: &str,
    other: &str,
    num_vars_of_one: u64,
    num_vars_of_other: u64,
) -> String {
    if one == other {
        return format!(
            "the individual `{one}` has no called genotype among the variants that were \
             used, so its entry of the kinship would be divided by no variant at all; \
             leave it out"
        );
    }
    format!(
        "the individuals `{one}` and `{other}` have no variant called in both of them, \
         so their entry of the kinship would be divided by no variant at all: \
         {num_vars_of_one} {said} called in the first and {num_vars_of_other} in the \
         second; leave one of the two out",
        said = if num_vars_of_one == 1 {
            "variant is"
        } else {
            "variants are"
        },
    )
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
