//! What a TypeScript user reaches through `calcRogersHuffR2Matrix` and
//! `calcLdAndDistPerPop`: the r² of every pair of the variants of a source,
//! and how that r² falls off with the distance between the two variants of
//! a pair, in bins of distance and as the curve fitted to the pairs of each
//! population.
//!
//! The calculations are the core's, `popnei::ld::calc_r2_matrix` and
//! `popnei::ld::calc_ld_and_dist`, and what
//! this module does is the translation that section 11 of
//! `docs/architecture.md` leaves to a binding crate. It builds the chain of
//! readers of the pass from the steps of the `Variants`, keeps that chain
//! while the calculation runs so that the counts of its filters can be read
//! when it returns, turns the number of the chromosome of each variant into
//! the name the reader gave it and its position into the float64 a number
//! of JavaScript is, and, for a pass that gave no variant, says whether the
//! source had none or the steps kept none, which only the chain knows.
//!
//! The bins of the fall-off carry one missing value that is decided here: a
//! bin with no pair, whose mean and standard deviation the core gives as
//! `None` and which a user reads as NaN. A pair that has no r² is NaN in
//! the core already, and so are the three values of the curve of a
//! population that has none.
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

use popnei::ld::{
    DEFAULT_MAX_ALLOWED_MAF, DEFAULT_MAX_DIST, DEFAULT_MIN_DIST, DEFAULT_NUM_DIST_BINS,
    LdAndDistOptions, LdBins, MAX_NUM_VARS_OF_THE_MATRIX, TheMatrixGivenAway, calc_ld_and_dist,
    calc_r2_matrix,
};
use popnei::stats::Pops;

use crate::errors::JsPopneiError;
use crate::source::{Consumer, OpenSource, PassCounts, positions_of, the_run_of};
use crate::stats::{PopsGiven, pops_of_the_arrays};
use crate::steps::{Steps, chain_of};

/// The name of the smallest distance a pair is counted at, as a TypeScript
/// user writes it.
const MIN_DIST: &str = "minDist";

/// The name of the largest distance a pair is counted at, the same way.
const MAX_DIST: &str = "maxDist";

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
    the_run_of(source, &Consumer::R2Matrix, |run| {
        let reader = source.reader(run, None)?;
        let mut chain = chain_of(reader, steps.steps())?;
        let matrix = calc_r2_matrix(&mut chain, max_num_vars).map_err(of_this_pass)?;
        let num_vars = matrix.num_vars();
        let counted = u64::try_from(num_vars).map_err(|_| {
            JsPopneiError::Broken(format!(
                "the pass gave {num_vars} variants, more than the count of a pass holds"
            ))
        })?;
        let counts = PassCounts::of(counted, &chain.filtering_stats());
        // The matrix is taken out of the core's result and not read from it,
        // so that what crosses into JavaScript is the allocation the core
        // filled. `given_away` consumes that result, so the chromosomes and
        // the positions are read from what it gave and not from it.
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

/// The smallest distance in base pairs at which a pair of variants is
/// counted in the fall-off of r², when the caller says nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_min_dist() -> f64 {
    DEFAULT_MIN_DIST as f64
}

/// The largest distance in base pairs at which a pair of variants is
/// counted there, when the caller says nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_max_dist() -> f64 {
    DEFAULT_MAX_DIST as f64
}

/// How many bins of distance the pairs are put into, when the caller says
/// nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_num_dist_bins() -> usize {
    DEFAULT_NUM_DIST_BINS
}

/// The largest major allele frequency a variant has in a population and is
/// still counted there, when the caller says nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_max_allowed_maf() -> f64 {
    DEFAULT_MAX_ALLOWED_MAF
}

/// The last whole number a number of JavaScript holds, 2^53: the one after
/// it arrives as this one, because a float64 counts in twos from there.
///
/// It is the number `LARGEST_POSITION` of `source.rs` puts on the position
/// of a variant, here for the counts and the distances of the bins.
const LARGEST_WHOLE_NUMBER: u64 = 9_007_199_254_740_992;

/// The arguments of one pass of the fall-off of r² with distance, as they
/// crossed from TypeScript.
///
/// The package has checked that each of them is of the type the core takes,
/// since a number of JavaScript reaches a whole number of the core as 32
/// bits with no error; what is left is what the core says of them, a
/// `min_dist` above `max_dist` and a `max_allowed_maf` outside 0 to 1 among
/// it.
pub(crate) struct ArgumentsOfTheBins {
    /// The name of each population, in the order the user gave them, and
    /// nothing when they named none, which is one population of every
    /// individual.
    pub(crate) pop_names: Option<Vec<String>>,
    /// The names of the individuals of every population, the ones of the
    /// first population first.
    pub(crate) pop_individuals: Vec<String>,
    /// How many individuals each population of `pop_names` holds, which
    /// cuts `pop_individuals` into the names of each of them.
    pub(crate) num_individuals_per_pop: Vec<u32>,
    /// The smallest distance in base pairs at which a pair is counted,
    /// which crosses as a float64: the core takes a distance of up to
    /// 2^64 - 1, and a number of JavaScript counts one by one up to
    /// 2^53 - 1.
    pub(crate) min_dist: f64,
    /// The largest distance in base pairs at which a pair is counted.
    pub(crate) max_dist: f64,
    /// How many bins of equal width those distances are cut into.
    pub(crate) num_bins: usize,
    /// The largest major allele frequency a variant has in a population and
    /// is still counted there.
    pub(crate) max_allowed_maf: f64,
}

/// What one pass of the fall-off gives TypeScript: the names of the
/// populations in their order, the five values of every bin of every one of
/// them, how many variants each of them kept at its major allele frequency,
/// the three values of the curve fitted to the pairs of each of them, and
/// the counts of the pass.
///
/// The bins of every population are in one array each, the bins of one
/// population after the bins of the one before it: an array of arrays is
/// not one of the types wasm-bindgen carries, so the package cuts them, as
/// section 11 of `docs/architecture.md` has it for a table that crosses
/// with a copy. The curve of a population is three numbers and not three
/// arrays, so its three arrays hold one value for each population and the
/// package reads them by the index of the population and cuts nothing.
///
/// Every array leaves the memory of wasm the first time it is asked for,
/// and the call after that gives nothing: the package reads each of them
/// once, into the object a user holds, and frees this.
#[wasm_bindgen]
pub struct LdAndDistOfAPass {
    pop_names: Option<Vec<String>>,
    num_bins: usize,
    smallest_dist: Option<Vec<f64>>,
    largest_dist: Option<Vec<f64>>,
    num_pairs: Option<Vec<f64>>,
    mean_r2: Option<Vec<f64>>,
    sd_r2: Option<Vec<f64>>,
    num_vars_per_pop: Option<Vec<f64>>,
    rho_per_bp: Option<Vec<f64>>,
    r2_at_zero: Option<Vec<f64>>,
    half_dist: Option<Vec<f64>>,
    counts: PassCounts,
}

#[wasm_bindgen]
impl LdAndDistOfAPass {
    /// The names of the populations, in the order the user named them,
    /// which is the order the bins of every array below are in.
    pub fn pop_names(&mut self) -> Option<Vec<String>> {
        self.pop_names.take()
    }

    /// How many bins of distance each population holds, which is what cuts
    /// every array below into the bins of each of them.
    #[must_use]
    pub fn num_bins(&self) -> usize {
        self.num_bins
    }

    /// The smallest distance of every bin, that distance included, the bins
    /// of one population after the bins of the one before it.
    pub fn smallest_dist(&mut self) -> Option<Vec<f64>> {
        self.smallest_dist.take()
    }

    /// The largest distance of every bin, that distance included, in the
    /// same order.
    pub fn largest_dist(&mut self) -> Option<Vec<f64>> {
        self.largest_dist.take()
    }

    /// How many pairs of variants fell in every bin, in the same order.
    pub fn num_pairs(&mut self) -> Option<Vec<f64>> {
        self.num_pairs.take()
    }

    /// The mean of the r² of the pairs of every bin, in the same order, and
    /// NaN for a bin with no pair.
    pub fn mean_r2(&mut self) -> Option<Vec<f64>> {
        self.mean_r2.take()
    }

    /// The standard deviation of those r², with the pairs of the bin as the
    /// divisor, in the same order, and NaN for a bin with no pair.
    pub fn sd_r2(&mut self) -> Option<Vec<f64>> {
        self.sd_r2.take()
    }

    /// How many variants each population kept at its major allele
    /// frequency, one for each of `pop_names`.
    pub fn num_vars_per_pop(&mut self) -> Option<Vec<f64>> {
        self.num_vars_per_pop.take()
    }

    /// The fitted 4Nr per base pair of the curve of every population, one
    /// for each of `pop_names`, and NaN for a population no curve was
    /// fitted to.
    pub fn rho_per_bp(&mut self) -> Option<Vec<f64>> {
        self.rho_per_bp.take()
    }

    /// That curve at a distance of 0, in the same order, and NaN for the
    /// same populations.
    pub fn r2_at_zero(&mut self) -> Option<Vec<f64>> {
        self.r2_at_zero.take()
    }

    /// The distance in base pairs at which that curve has fallen to half of
    /// its value at a distance of 0, in the same order.
    ///
    /// It is NaN for a population no curve was fitted to, and NaN on its
    /// own for a curve that never falls to half, which no pass reaches.
    pub fn half_dist(&mut self) -> Option<Vec<f64>> {
        self.half_dist.take()
    }

    /// How many variants the pass gave, before the major allele frequency
    /// of any population, and what each filter of it was given and kept.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.counts.clone()
    }
}

/// How the r² of a pair of variants falls off with the distance between
/// them, in bins of distance and as the curve fitted to the pairs of each
/// population of `asked`, over one pass of `source` through the steps of
/// `steps`.
///
/// The chain of readers of the pass stays here, lent to the core, so that
/// the counts of its filters can be read when the calculation returns: the
/// loop over the blocks is the core's, and how many variants it took is
/// `LdAndDist::num_vars`, since no block of the pass reaches this crate.
/// The reader is asked for no size of block: the bins are added up in the
/// order of the variants of the pass, which no block cuts, so the same
/// numbers come out whatever the size.
///
/// # Errors
///
/// When a distance did not arrive as a whole number of base pairs, which is
/// a defect of the package; when `min_dist` is above `max_dist`, when
/// `num_bins` is 0 and when `max_allowed_maf` is not a number from 0 to 1;
/// when a population names an individual the pass does not give, names one
/// twice or names none, and when `pops` holds no population; when the pass
/// gives no variant; when the memory of the tab does not take the bins, the
/// window of blocks or the r² of a step; when a bin holds more pairs than a
/// number of JavaScript counts one by one; when the core gives no bounds
/// and no pairs for a bin it says it has, which is a defect of popnei; and
/// when the source cannot be read, a wrong line of a VCF among the causes.
pub(crate) fn ld_and_dist_of(
    source: &dyn OpenSource,
    steps: &Steps,
    asked: &ArgumentsOfTheBins,
) -> Result<LdAndDistOfAPass, JsPopneiError> {
    let options = LdAndDistOptions {
        min_dist: the_distance_of(MIN_DIST, asked.min_dist)?,
        max_dist: the_distance_of(MAX_DIST, asked.max_dist)?,
        num_bins: asked.num_bins,
        max_allowed_maf: asked.max_allowed_maf,
    };
    let named = the_pops_given(asked)?;
    the_run_of(source, &Consumer::LdAndDist, |run| {
        let reader = source.reader(run, None)?;
        let mut chain = chain_of(reader, steps.steps())?;
        // The names a user wrote are looked up among the individuals the pass
        // gives, which are those of the source after a filter of individuals
        // when the `Variants` carries one, and only the pass knows them. It is
        // the lookup the statistics per population make, so a name that is not
        // an individual of the pass is refused in the same words by both.
        let pops = match named {
            Some(named) => Pops::from_names(&named, chain.individuals())?,
            None => Pops::all(chain.individuals().len()),
        };
        let pop_names: Vec<String> = (0..pops.len())
            .map(|pop| pops.name(pop).to_owned())
            .collect();
        let of_each_pop: Vec<&[usize]> = (0..pops.len()).map(|pop| pops.individuals(pop)).collect();
        let of_the_pass = calc_ld_and_dist(&mut *chain, &of_each_pop, &options)?;
        // The variants the pass gave are the core's count and are not worked
        // out again here: the calculation was given them and counted them with
        // the arithmetic that says what happens on overflow.
        let counts = PassCounts::of(of_the_pass.num_vars(), &chain.filtering_stats());
        let num_pops = of_the_pass.num_pops();
        let num_bins = options.num_bins;
        let of_every_bin = num_pops.saturating_mul(num_bins);
        let mut smallest_dist = Vec::with_capacity(of_every_bin);
        let mut largest_dist = Vec::with_capacity(of_every_bin);
        let mut num_pairs = Vec::with_capacity(of_every_bin);
        let mut mean_r2 = Vec::with_capacity(of_every_bin);
        let mut sd_r2 = Vec::with_capacity(of_every_bin);
        let mut num_vars_per_pop = Vec::with_capacity(num_pops);
        let mut rho_per_bp = Vec::with_capacity(num_pops);
        let mut r2_at_zero = Vec::with_capacity(num_pops);
        let mut half_dist = Vec::with_capacity(num_pops);
        for pop in 0..num_pops {
            let bins = of_the_pass.bins_of_pop(pop).ok_or_else(|| {
                JsPopneiError::Broken(format!(
                    "the pass counted {num_pops} populations and has no bins for the \
                     population {pop}"
                ))
            })?;
            // The package cuts every array below by `num_bins`, so a population
            // whose bins were not that many would give its user the bins of the
            // next one.
            if bins.num_bins() != num_bins {
                return Err(JsPopneiError::Broken(format!(
                    "the pass was asked for {num_bins} bins of distance and the population \
                     {pop} has {its_bins} of them, which is a defect of popnei; please \
                     report it",
                    its_bins = bins.num_bins()
                )));
            }
            for bin in 0..num_bins {
                let (smallest, largest) = bins
                    .bounds(bin)
                    .ok_or_else(|| not_a_bin("the distances", bin, bins))?;
                let pairs = bins
                    .num_pairs(bin)
                    .ok_or_else(|| not_a_bin("the pairs", bin, bins))?;
                smallest_dist.push(as_a_number(smallest, "the smallest distance of a bin")?);
                largest_dist.push(as_a_number(largest, "the largest distance of a bin")?);
                num_pairs.push(as_a_number(pairs, "the pairs of a bin")?);
                // A bin with no pair has no mean and no standard deviation in
                // the core, and NaN is what the boundary with a language that
                // has no missing number writes.
                mean_r2.push(bins.mean_r2(bin).unwrap_or(f64::NAN));
                sd_r2.push(bins.sd_r2(bin).unwrap_or(f64::NAN));
            }
            num_vars_per_pop.push(as_a_number(
                bins.num_vars(),
                "the variants a population kept",
            )?);
            // The three values of the curve are the core's as they are, the
            // NaN of a population no curve was fitted to included: nothing of
            // the fit is worked out here.
            let curve = bins.decay();
            rho_per_bp.push(curve.rho_per_bp());
            r2_at_zero.push(curve.r2_at_zero());
            half_dist.push(curve.half_dist());
        }
        Ok(LdAndDistOfAPass {
            pop_names: Some(pop_names),
            num_bins,
            smallest_dist: Some(smallest_dist),
            largest_dist: Some(largest_dist),
            num_pairs: Some(num_pairs),
            mean_r2: Some(mean_r2),
            sd_r2: Some(sd_r2),
            num_vars_per_pop: Some(num_vars_per_pop),
            rho_per_bp: Some(rho_per_bp),
            r2_at_zero: Some(r2_at_zero),
            half_dist: Some(half_dist),
            counts,
        })
    })
}

/// The populations a user named, each with the names of its individuals,
/// out of the flat arrays they crossed in, and `None` when they named none.
///
/// # Errors
///
/// When the arrays do not hold the individuals of every population, which
/// is a defect of the package: it is what cuts them.
fn the_pops_given(asked: &ArgumentsOfTheBins) -> Result<Option<PopsGiven>, JsPopneiError> {
    let Some(names) = asked.pop_names.as_ref() else {
        return Ok(None);
    };
    Ok(Some(pops_of_the_arrays(
        names,
        &asked.pop_individuals,
        &asked.num_individuals_per_pop,
    )?))
}

/// The distance `given` as the number of base pairs the core takes.
///
/// # Errors
///
/// When it is not a whole number from 0 to 2^53, which is a defect of the
/// package and not something a user can write: `js/popnei/src/arguments.ts`
/// refuses a distance that is not a whole number of 0 or more before the
/// call, the negatives among them. It is checked here and not cast as it
/// comes because a NaN would arrive as a distance of 0 and count the pairs
/// of two variants at one position, saying nothing.
fn the_distance_of(argument: &str, given: f64) -> Result<u64, JsPopneiError> {
    #[expect(
        clippy::float_cmp,
        reason = "a float64 is or is not the whole number it was truncated to"
    )]
    let is_whole = given.trunc() == given;
    if !given.is_finite() || !is_whole || given < 0.0 || given > LARGEST_WHOLE_NUMBER as f64 {
        return Err(JsPopneiError::Broken(format!(
            "`{argument}` arrived as {given}, and it is a whole number of base pairs \
             from 0 to {LARGEST_WHOLE_NUMBER}, which is a defect of popnei; please \
             report it"
        )));
    }
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "checked above to be a whole number between 0 and 2^53"
    )]
    let base_pairs = given as u64;
    Ok(base_pairs)
}

/// `count` as the number of JavaScript that holds it, `what` naming what
/// was counted.
///
/// # Errors
///
/// [`JsPopneiError::NotInJavaScript`] when it is above 2^53, which a float64
/// counts past in twos: the number would arrive as another one and say
/// nothing. A pass of a tab does not reach it, whose file is in a memory of
/// 2^32 bytes, and the check is what says so instead of assuming it.
fn as_a_number(count: u64, what: &str) -> Result<f64, JsPopneiError> {
    if count > LARGEST_WHOLE_NUMBER {
        return Err(JsPopneiError::NotInJavaScript(format!(
            "{what} is {count}, which is above {LARGEST_WHOLE_NUMBER}, the last whole \
             number a number of JavaScript holds: the one after it would arrive as \
             {LARGEST_WHOLE_NUMBER} itself"
        )));
    }
    Ok(count as f64)
}

/// What a user is told when the core says it has that many bins and gives
/// nothing for one of them, `what` naming the values that were asked for.
fn not_a_bin(what: &str, bin: usize, bins: &LdBins) -> JsPopneiError {
    JsPopneiError::Broken(format!(
        "{what} of the bin {bin} of {num_bins} are not in the bins of a population, \
         which is a defect of popnei; please report it",
        num_bins = bins.num_bins()
    ))
}
