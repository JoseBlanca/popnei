//! What a TypeScript user reaches through `calcRogersHuffR2Matrix`: the r²
//! of every pair of the variants of a source.
//!
//! The calculation is the core's, `popnei::ld::calc_r2_matrix`, and what
//! this module does is the translation that section 11 of
//! `docs/architecture.md` leaves to a binding crate. It builds the chain of
//! readers of the pass from the steps of the `Variants`, keeps that chain
//! while the calculation runs so that the counts of its filters can be read
//! when it returns, turns the number of the chromosome of each variant into
//! the name the reader gave it and its position into the float64 a number
//! of JavaScript is, and, for a pass that gave no variant, says whether the
//! source had none or the steps kept none, which only the chain knows.
//!
//! [`R2Matrix`] is the result on its way out. It lives in the memory of
//! wasm, which the garbage collector of JavaScript does not see, so the
//! package frees it as soon as its three arrays are read, and each of them
//! leaves it as it is read.
//!
//! The matrix is 8 bytes for each pair of the variants of the pass, 200 MB
//! at the 5000 variants of [`MAX_NUM_VARS_OF_THE_MATRIX`], and this crate
//! copies none of it inside wasm: `R2Matrix::given_away` of the core hands
//! the `Vec` over and this one holds that same allocation until
//! wasm-bindgen moves it out. So the memory of wasm holds the matrix once
//! and not twice, 200 MB at that cap and not 400 MB, which is what a tab
//! keeps for its lifetime: a module gives no page back to the host, so the
//! high-water mark of one call is the ceiling of every call after it.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::ld::{MAX_NUM_VARS_OF_THE_MATRIX, TheMatrixGivenAway, calc_r2_matrix};

use crate::errors::JsPopneiError;
use crate::source::{OpenSource, PassCounts, positions_of};
use crate::steps::{Steps, chain_of};

/// The r² of every pair of the variants of a pass, with the chromosome and
/// the position of each of them and the counts of that pass.
///
/// The matrix is `num_vars` rows of `num_vars` values, row after row, and a
/// pair that has no r² is NaN: the core has it so, and NaN is what the
/// boundary with a language that has no missing value writes. The two cells
/// of a pair hold the same value and the diagonal of a variant with two
/// dosages at least is 1.
///
/// Each array leaves the memory of wasm the first time it is asked for, and
/// the call after that gives nothing: the package reads each of them once,
/// into the object a user holds, and frees this. A copy left behind would
/// grow the memory of wasm, which never gives memory back, by the whole
/// matrix a second time.
#[wasm_bindgen]
pub struct R2Matrix {
    num_vars: usize,
    /// The r² of every pair, and `None` once it was given to JavaScript.
    r2: Option<Vec<f64>>,
    /// The name of the chromosome of each variant, and `None` once it was
    /// given to JavaScript.
    chroms: Option<Vec<String>>,
    /// The position of each variant, and `None` once it was given to
    /// JavaScript.
    poss: Option<Vec<f64>>,
    /// The counts of the pass, which the package turns into the `passStats`
    /// of the result.
    counts: PassCounts,
}

#[wasm_bindgen]
impl R2Matrix {
    /// How many variants the matrix is of, which is how many the pass gave.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// The r² of every pair, `num_vars` x `num_vars` row after row, NaN for
    /// a pair that has none, or `undefined` when they were read already.
    pub fn r2(&mut self) -> Option<Vec<f64>> {
        self.r2.take()
    }

    /// The name of the chromosome of each variant, in the order the pass
    /// gave them, or `undefined` when they were read already.
    pub fn chroms(&mut self) -> Option<Vec<String>> {
        self.chroms.take()
    }

    /// The position of each variant, 1 based as in a VCF, or `undefined`
    /// when they were read already.
    pub fn poss(&mut self) -> Option<Vec<f64>> {
        self.poss.take()
    }

    /// How many variants the calculation took, and what each filter of the
    /// pass was given and kept.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.counts.clone()
    }
}

/// How many variants the matrix is taken of before it is refused, when the
/// caller says nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_max_num_vars() -> usize {
    MAX_NUM_VARS_OF_THE_MATRIX
}

/// The largest cap a user of this build can put on the variants of the
/// matrix, which is 65535 in a browser.
///
/// The matrix holds one r² for each pair, which is the variants squared,
/// and that number is counted in a `usize`, 32 bits in WebAssembly: the
/// matrix of 65535 variants holds 4294836225 values, and that of 65536
/// holds more than the 4294967295 a browser counts to, which the core
/// refuses before it reads anything. The package reads this number and
/// refuses a larger cap at the call, in `js/popnei/src/arguments.ts`, so
/// that the number a user is given is the one their machine has: a Python
/// user of the same package on a 64 bit build has 4294967295.
///
/// The memory of the tab stops a user well before it. The matrix of 23170
/// variants is 4 GB, which is everything a page holds at a time.
#[wasm_bindgen]
#[must_use]
pub fn largest_max_num_vars() -> usize {
    usize::MAX.isqrt()
}

/// The r² of every pair of the variants of `source` that the steps of
/// `steps` keep, with the chromosome and the position of each of them and
/// the counts of the pass.
///
/// The chain of readers of the pass stays here, lent to the core, so that
/// the counts of its filters can be read when the calculation returns: the
/// loop over the blocks is the core's, and how many variants it took is
/// `R2Matrix::num_vars`, since no block of the pass reaches this crate. The
/// reader is asked for no size of block: the core puts the variants of the
/// pass into tiles of its own and the matrix is the same, to the bit,
/// whatever size the blocks had.
///
/// # Errors
///
/// When the pass gives more than `max_num_vars` variants, with both numbers
/// and the memory the matrix would have needed, under the `maxNumVars` a
/// TypeScript user wrote; when the pass gives no variant, which says
/// whether the source had none or the steps kept none; when the memory of
/// the tab does not take the matrix; and when the source cannot be read, a
/// wrong line of a VCF among the causes.
///
/// A `max_num_vars` the machine does not count the pairs of is the core's
/// error and reaches no user of the package, which refuses a cap above
/// [`largest_max_num_vars`] at the call.
pub(crate) fn r2_matrix_of(
    source: &dyn OpenSource,
    max_num_vars: usize,
    steps: Steps,
) -> Result<R2Matrix, JsPopneiError> {
    let reader = source.reader(None)?;
    let mut chain = chain_of(reader, steps.steps())?;
    let matrix = calc_r2_matrix(&mut chain, max_num_vars).map_err(of_this_pass)?;
    let num_vars = matrix.num_vars();
    let counted = u64::try_from(num_vars).map_err(|_| {
        JsPopneiError::Broken(format!(
            "the pass gave {num_vars} variants, more than the count of a pass holds"
        ))
    })?;
    let counts = PassCounts::of(counted, &chain.filtering_stats());
    // The matrix is taken out of the core's result and not read from it, so
    // that what crosses into JavaScript is the allocation the core filled.
    // `given_away` consumes that result, so the chromosomes and the
    // positions are read from what it gave and not from it.
    let matrix = matrix.given_away();
    let chroms = the_names_of_the_chromosomes(&matrix)?;
    let poss = positions_of(&matrix.poss)?;
    Ok(R2Matrix {
        num_vars,
        r2: Some(matrix.r2),
        chroms: Some(chroms),
        poss: Some(poss),
        counts,
    })
}

/// What a pass of this calculation failed with, on its way to a TypeScript
/// user: the cap on its variants under the name that user wrote it in, and
/// everything else as the core says it.
///
/// The core names the cap `max_num_vars`, which is the argument of the
/// Python package, and the call a TypeScript user has to look at is
/// `calcRogersHuffR2Matrix(variants, { maxNumVars })`. It is what
/// `under_the_argument` of `steps.rs` does for the threshold of a filter.
///
/// A pass that gave no variant is one of the rest: the core builds that
/// message itself, with the counts it reads from the chain it was lent, so
/// that it is the same sentence whichever calculation asked for the pass.
fn of_this_pass(error: popnei::Error) -> JsPopneiError {
    if let popnei::Error::LdTooManyVars {
        num_vars,
        max_num_vars,
        bytes,
    } = error
    {
        return JsPopneiError::TooManyVars {
            num_vars,
            max_num_vars,
            bytes,
        };
    }
    JsPopneiError::Core(error)
}

/// The name of the chromosome of each variant of `matrix`, read through the
/// table of names the core cloned from the reader of the pass.
///
/// # Errors
///
/// When a number of a variant is not in that table, which cannot happen
/// unless this crate or the core has a defect.
fn the_names_of_the_chromosomes(matrix: &TheMatrixGivenAway) -> Result<Vec<String>, JsPopneiError> {
    let table = &matrix.chrom_table;
    matrix
        .chroms
        .iter()
        .map(|number| {
            table.name(*number).map(str::to_owned).ok_or_else(|| {
                JsPopneiError::Broken(format!(
                    "the chromosome number {number} of the matrix of r² is not in the \
                     table of the reader that gave it"
                ))
            })
        })
        .collect()
}
