//! The density of the variants along the chromosomes on its way between
//! TypeScript and the core.
//!
//! It turns the `windowSize` and the `chromLengths` a TypeScript user wrote,
//! which the package checked to be whole numbers that a number of JavaScript
//! counts to one by one, into the numbers the core takes; builds the chain
//! of readers of the pass from the steps of the `Variants`, lends it to the
//! core, and reads the counts of its filters when the pass is over. What
//! goes out is [`VarDensityOfAPass`]: the name of each chromosome with how
//! many windows it has, and the start, the end and the count of every
//! window, which the package makes into the arrays of the result. The core
//! gives each start and end, so nothing here works one out.
//!
//! The lengths cross flat, as the populations do: the names in one array
//! and the lengths in another, at the same places.
//!
//! "The density of the variants along the chromosomes" of
//! `docs/specs/stats.md` has the design.

use wasm_bindgen::JsValue;
use wasm_bindgen::prelude::wasm_bindgen;

use popnei::stats::{SoFar, VarDensity};

use crate::errors::JsPopneiError;
use crate::source::{
    Consumer, LARGEST_POSITION, OpenSource, PassCounts, TheResultSoFar, the_run_of,
};
use crate::stats::TheResultSoFarAsked;
use crate::steps::{Steps, chain_of};

/// The arguments of one pass, as they crossed from TypeScript.
pub(crate) struct ArgumentsOfTheDensity {
    /// The width of a window in base pairs.
    pub(crate) window_size: f64,
    /// The names of the chromosomes of `chromLengths`, in the order the
    /// object iterates in, and `None` when the user gave no lengths, which
    /// reads those of the source.
    pub(crate) chrom_names: Option<Vec<String>>,
    /// The length of each chromosome of `chrom_names`, at the same place.
    pub(crate) chrom_lengths: Vec<f64>,
}

/// The number of variants in each window along each chromosome over one pass
/// over `source`, through the steps of `steps`.
///
/// The function of `so_far` is given the density over the variants read so
/// far, built as the final one is, with the counts of the pass as they stand
/// after the block.
///
/// # Errors
///
/// A width or a length that is not a whole number from 1 to 2^53 - 1, and
/// names and lengths of two sizes, each a defect of the package, which
/// refuses them before the call; what the core refuses, a variant past the
/// length of its chromosome among it; a window whose end is past 2^53, which
/// a number of JavaScript would round; the memory of a page, when it cannot
/// hold the arrays of the windows; a source that cannot be read; and the
/// value the function of `so_far` threw.
pub(crate) fn var_density_of(
    source: &dyn OpenSource,
    steps: &Steps,
    asked: &ArgumentsOfTheDensity,
    so_far: TheResultSoFarAsked,
) -> Result<VarDensityOfAPass, JsPopneiError> {
    let window_size = the_base_pairs_of("windowSize", asked.window_size)?;
    let chrom_lengths = asked
        .chrom_names
        .as_ref()
        .map(|names| {
            if names.len() != asked.chrom_lengths.len() {
                return Err(JsPopneiError::Broken(format!(
                    "`chromLengths` arrived as {names} names and {lengths} lengths, and \
                     it is one length for each name, which is a defect of popnei; please \
                     report it",
                    names = names.len(),
                    lengths = asked.chrom_lengths.len()
                )));
            }
            names
                .iter()
                .zip(&asked.chrom_lengths)
                .map(|(name, length)| {
                    Ok((name.clone(), the_base_pairs_of("chromLengths", *length)?))
                })
                .collect::<Result<Vec<(String, u64)>, JsPopneiError>>()
        })
        .transpose()?;
    the_run_of(source, &Consumer::VarDensity, |run| {
        let reader = source.reader(run, None)?;
        let mut chain = chain_of(reader, steps.steps())?;
        let mut told = TheResultSoFar::from_now(so_far.told, so_far.every_seconds)?;
        let given = popnei::stats::calc_var_density_with(
            &mut *chain,
            window_size,
            chrom_lengths.as_deref(),
            &mut |added_up: &dyn SoFar<VarDensity>| {
                let Some(told) = told.as_mut() else {
                    return Ok(());
                };
                told.after_a_block(run, || {
                    let counts = PassCounts::of_the_filters(
                        added_up.num_vars(),
                        steps.steps(),
                        &added_up.filtering_stats(),
                    );
                    Ok(JsValue::from(windows_of(&added_up.result()?, counts)?))
                })
            },
        )
        .map_err(JsPopneiError::from);
        let density = TheResultSoFar::what_the_pass_gives(told, given)?;
        let counts = PassCounts::of(density.num_vars(), steps.steps(), &*chain);
        windows_of(&density, counts)
    })
}

/// `density`, the windows the core counted over the variants of a pass or of
/// its first blocks, on their way to JavaScript, with `counts`, the counts of
/// the pass.
///
/// # Errors
///
/// When the memory of a page cannot hold the arrays of the windows, and when
/// a window ends past 2^53, which a number of JavaScript would round.
fn windows_of(
    density: &VarDensity,
    counts: PassCounts,
) -> Result<VarDensityOfAPass, JsPopneiError> {
    let mut chroms = Vec::with_capacity(density.chroms().len());
    let mut windows_per_chrom = Vec::with_capacity(density.chroms().len());
    for chrom in density.chroms() {
        chroms.push(chrom.name.clone());
        // A chromosome has at most `MAX_NUM_WINDOWS` windows, which a
        // `u32` holds.
        windows_per_chrom.push(u32::try_from(chrom.counts.len()).unwrap_or(u32::MAX));
    }
    let num_windows = density.num_windows();
    let mut starts = Vec::new();
    let mut ends = Vec::new();
    let mut num_vars = Vec::new();
    // 20 bytes a window, 200 MB for the most windows a density has, which
    // an allocation that fails in wasm would abort with instead of an
    // `Error`.
    let no_room = || {
        JsPopneiError::NoMemory(format!(
            "the {num_windows} windows of the density do not fit in the memory popnei \
             has left: a page holds at most 4 GB, and every file that is open counts; \
             a wider window gives fewer of them"
        ))
    };
    starts
        .try_reserve_exact(num_windows)
        .map_err(|_| no_room())?;
    ends.try_reserve_exact(num_windows).map_err(|_| no_room())?;
    num_vars
        .try_reserve_exact(num_windows)
        .map_err(|_| no_room())?;
    for window in density.windows() {
        // The end is the larger of the two, so a window that ends at
        // 2^53 or before starts there too.
        if window.end > LARGEST_POSITION {
            return Err(JsPopneiError::NotInJavaScript(format!(
                "the window {start} to {end} of the chromosome {chrom} ends past \
                 {LARGEST_POSITION}, the last whole number a number of JavaScript holds: \
                 the one after it would be read as another, and the same file read from \
                 Python gives the end the length or the window has",
                start = window.start,
                end = window.end,
                chrom = window.chrom,
            )));
        }
        starts.push(position_in_javascript(window.start));
        ends.push(position_in_javascript(window.end));
        num_vars.push(window.num_vars);
    }
    Ok(VarDensityOfAPass {
        chroms,
        windows_per_chrom,
        starts,
        ends,
        num_vars,
        counts,
    })
}

/// `position`, which is at most [`LARGEST_POSITION`], as the number of
/// JavaScript that holds it exactly.
#[expect(
    clippy::cast_precision_loss,
    reason = "at most 2^53, which a float64 holds exactly, as the caller checked"
)]
fn position_in_javascript(position: u64) -> f64 {
    position as f64
}

/// `given`, the width of a window or the length of a chromosome, as the
/// number of base pairs the core takes.
///
/// # Errors
///
/// When it is not a whole number from 1 to 2^53 - 1, which is a defect of the
/// package and not something a user can write: `js/popnei/src/stats.ts`
/// refuses anything else before the call. It is checked here and not cast
/// as it comes because a NaN would arrive as 0.
fn the_base_pairs_of(argument: &str, given: f64) -> Result<u64, JsPopneiError> {
    #[expect(
        clippy::float_cmp,
        reason = "a float64 is or is not the whole number it was truncated to"
    )]
    let is_whole = given.trunc() == given;
    // 2^53 - 1, `Number.MAX_SAFE_INTEGER`, which a float64 holds exactly.
    let largest = 9_007_199_254_740_991.0_f64;
    if !given.is_finite() || !is_whole || given < 1.0 || given > largest {
        return Err(JsPopneiError::Broken(format!(
            "`{argument}` arrived as {given}, and it is a whole number of base pairs from 1 \
             to {largest}, which is a defect of popnei; please report it"
        )));
    }
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "checked above to be a whole number between 1 and 2^53 - 1"
    )]
    let base_pairs = given as u64;
    Ok(base_pairs)
}

/// What one pass of the density of the variants gives JavaScript.
///
/// Every array is copied out of the memory of wasm as it is read, and the
/// object itself holds that memory until its `free()` is called, which the
/// package does as soon as it has read every array of it.
#[wasm_bindgen]
pub struct VarDensityOfAPass {
    chroms: Vec<String>,
    windows_per_chrom: Vec<u32>,
    starts: Vec<f64>,
    ends: Vec<f64>,
    num_vars: Vec<u32>,
    counts: PassCounts,
}

#[wasm_bindgen]
impl VarDensityOfAPass {
    /// The name of each chromosome, in the order of the result: those with
    /// a length in the order of the lengths, then those with no length in
    /// the order their first variant came.
    #[must_use]
    pub fn chroms(&self) -> Vec<String> {
        self.chroms.clone()
    }

    /// How many windows each chromosome of [`VarDensityOfAPass::chroms`]
    /// has, whose windows come one chromosome after another in the arrays
    /// below.
    #[must_use]
    pub fn windows_per_chrom(&self) -> Vec<u32> {
        self.windows_per_chrom.clone()
    }

    /// The first position of each window, counted from 1.
    #[must_use]
    pub fn starts(&self) -> Vec<f64> {
        self.starts.clone()
    }

    /// The last position of each window, included.
    #[must_use]
    pub fn ends(&self) -> Vec<f64> {
        self.ends.clone()
    }

    /// How many variants of the pass are in each window.
    #[must_use]
    pub fn num_vars(&self) -> Vec<u32> {
        self.num_vars.clone()
    }

    /// How many variants the pass gave, and what each filter of it was given
    /// and kept, the outermost filter first.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.counts.clone()
    }
}
