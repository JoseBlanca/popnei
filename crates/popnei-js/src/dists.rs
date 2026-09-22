//! What a TypeScript user reaches through `calcPairwiseKosmanDists`: the
//! Kosman distance of every pair of individuals of a source.
//!
//! The calculation is the core's, `popnei::dists::calc_kosman_sums`, and
//! what this module does is the translation that section 11 of
//! `docs/architecture.md` leaves to a binding crate. It builds the chain of
//! readers of the pass from the steps of the `Variants`, keeps that chain
//! while the calculation runs so that the counts of its filters can be read
//! when it returns, turns the two integers the core keeps for each pair into
//! the distance vector, writing NaN where a pair has no distance, and, for a
//! pass that gave no variant, says whether the source had none or the steps
//! kept none, which only the chain knows.
//!
//! [`KosmanDistances`] is the result on its way out. It lives in the memory
//! of wasm, which the garbage collector of JavaScript does not see, so the
//! package frees it as soon as its three parts are read. The vector is moved
//! out of it as it is read, and not cloned: the code wasm-bindgen generates
//! copies the values into a `Float64Array` of the JavaScript heap and then
//! frees the `Vec`, so the two live side by side while that copy is made,
//! 800 MB at 10000 individuals, and what is left afterwards is the 400 MB of
//! the array the user holds. A clone before the crossing would be 400 MB
//! more, and the memory of wasm never gives back what it grew by.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::dists::calc_kosman_sums;
use popnei::filters::FilteringStats;

use crate::errors::JsPopneiError;
use crate::source::{OpenSource, PassCounts};
use crate::steps::{Steps, chain_of};

/// The Kosman distance of every pair of individuals, the names of those
/// individuals and the counts of the pass the distances were calculated
/// over.
///
/// The pairs are in the order of the distance vector of
/// `docs/specs/dists.md`, (0, 1), (0, 2), ..., (1, 2), ..., the upper
/// triangle of the square matrix of the distances row by row, and a pair
/// with no distance is NaN: the core gives `None` for it, and NaN is what
/// the boundary with a language that has no such value writes.
#[wasm_bindgen]
pub struct KosmanDistances {
    /// The distances, and `None` once they were given to JavaScript: they
    /// leave the memory of wasm as they are read, so the copy that crosses
    /// is the only one.
    dist_vector: Option<Vec<f64>>,
    /// The names of the individuals, in the order of the source, and `None`
    /// once they were given to JavaScript.
    names: Option<Vec<String>>,
    /// The counts of the pass, which the package turns into the `passStats`
    /// of the result.
    counts: PassCounts,
}

#[wasm_bindgen]
impl KosmanDistances {
    /// The distance of every pair, NaN for a pair that has none, or
    /// `undefined` when they were read already.
    pub fn dist_vector(&mut self) -> Option<Vec<f64>> {
        self.dist_vector.take()
    }

    /// The names of the individuals the pairs are of, in the order the
    /// source has them, or `undefined` when they were read already.
    pub fn names(&mut self) -> Option<Vec<String>> {
        self.names.take()
    }

    /// How many variants the calculation took, and what each filter of the
    /// pass was given and kept.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.counts.clone()
    }
}

/// The Kosman distance of every pair of individuals of `source`, over the
/// variants the steps of `steps` keep, with no distance for a pair called at
/// fewer than `min_num_vars` variants.
///
/// The chain of readers of the pass stays here, lent to the core, so that
/// the counts of its filters can be read when the calculation returns: the
/// loop over the blocks is the core's, and how many variants it took is
/// `KosmanSums::num_vars`, since no block of the pass reaches this crate.
/// The reader is asked for the size of block popnei chooses for the
/// individuals of the source, which changes no distance: the sums of a pair
/// are integers, so where the blocks are cut does not reach the one division
/// that gives its distance.
///
/// # Errors
///
/// When the pass gives no variant, which says whether the source had none or
/// the steps kept none; when the sums of a pair go above what a `u32` holds;
/// when the memory of the tab does not take the two counts of every pair;
/// and when the source cannot be read, a wrong line of a VCF among the
/// causes.
pub(crate) fn kosman_dists_of(
    source: &dyn OpenSource,
    min_num_vars: u32,
    steps: Steps,
) -> Result<KosmanDistances, JsPopneiError> {
    let reader = source.reader(None)?;
    let mut chain = chain_of(reader, steps.steps())?;
    // The names are the reader's own, taken before the calculation borrows
    // it: the vector and the names then cannot be of two different sources.
    let names = chain.individuals().to_vec();
    let calculated = calc_kosman_sums(&mut chain);
    let sums = match calculated {
        Ok(sums) => sums,
        Err(error) => return Err(of_the_pass(error, &chain.filtering_stats())),
    };
    let counts = PassCounts::of(sums.num_vars(), &chain.filtering_stats());
    let dist_vector = sums
        .dists(min_num_vars)
        .map(|dist| dist.unwrap_or(f64::NAN))
        .collect();
    Ok(KosmanDistances {
        dist_vector: Some(dist_vector),
        names: Some(names),
        counts,
    })
}

/// `error`, and a pass that gave no variant under what the chain counted:
/// whether the source held no variant or the filters kept none, and how many
/// variants each filter was given and kept.
///
/// The core says that the reader gave no variant and nothing more, because
/// it is given a reader and not the chain it is the end of. The counts are
/// the ones a failed pass would otherwise lose, as "A pass that was not
/// finished" of `docs/specs/filters.md` says: the calculation read the
/// source to its end before it found that there was no variant, so every
/// filter has counted everything it was given.
///
/// The words are those of the Python function of `docs/specs/dists.md`, so
/// that a user who reads one message reads the other. Python writes the path
/// of the file and `: ` before them, which a TypeScript user has not: this
/// crate is given the bytes of a file and no name for it.
fn of_the_pass(
    error: popnei::Error,
    filtering: &[(&'static str, FilteringStats)],
) -> JsPopneiError {
    if !matches!(error, popnei::Error::ReaderGaveNoVariants) {
        return JsPopneiError::Core(error);
    }
    // The chain gives its filters the outermost first, and the one the
    // source feeds is the last of them: what it was given is what the source
    // gave. A pass with no filter reads the source itself.
    let from_the_source = match filtering.last() {
        Some((_kind, stats)) => stats.vars_processed,
        None => 0,
    };
    // The filters are named in the order of the steps, which is the order a
    // user wrote them in and the reverse of the chain's. A filter that was
    // given nothing is named too: what a user is looking for is which of
    // their steps is the one that emptied the pass, and a filter missing
    // from the list would read as a step that did not run.
    let counts: Vec<String> = filtering
        .iter()
        .rev()
        .map(|(kind, stats)| {
            format!(
                "the filter `{kind}` was given {vars_processed} variants and kept \
                 {vars_kept}",
                vars_processed = stats.vars_processed,
                vars_kept = stats.vars_kept
            )
        })
        .collect();
    let counts = counts.join(", ");
    if from_the_source == 0 {
        let message = "the source has no variant, and a calculation needs 1 variant \
                       at least";
        return JsPopneiError::NoVariant(if counts.is_empty() {
            message.to_owned()
        } else {
            format!("{message}: {counts}")
        });
    }
    JsPopneiError::NoVariant(format!(
        "the steps kept no variant of the {from_the_source} the source gave, and a \
         calculation needs 1 variant at least: {counts}"
    ))
}
