//! What a TypeScript user reaches through `calcPopDists`: how far apart
//! every pair of the populations of a source is, as the measures they asked
//! for with the standard error of each.
//!
//! The calculation is the core's, [`popnei::pop_dists::calc_pop_dist_sums`],
//! and what this module does is what `dists.rs` does for the distances
//! between individuals: it builds the chain of readers of the pass from the
//! steps of the `Variants`, keeps that chain while the calculation runs so
//! that the counts of its filters can be read when it returns, and turns the
//! six sums the core keeps for each pair into the arrays the package builds
//! its result from. What a measure has no value for, a pair whose
//! populations counted no variant together, is NaN, which is what the
//! boundary with a language that has no missing number writes.
//!
//! The arguments a TypeScript user wrote are turned here into what the core
//! takes: the measures, which the core names, how the variants are cut into
//! the resampling groups, and how many called genotypes a population needs
//! at a variant. The populations cross flat, the names of the individuals of
//! every one of them in one array, as those of the statistics do, and the
//! names themselves are looked up against the individuals of the pass, which
//! only the pass knows.
//!
//! [`PopDistsOfAPass`] is the result on its way out. It lives in the memory
//! of wasm, which the garbage collector of JavaScript does not see, so the
//! package frees it as soon as its arrays are read, and each of them is
//! moved out of it as it is read and not cloned. The values of every measure
//! come as one array, the pairs of one measure after the pairs of the one
//! before it, and the f_2 of every pair within every group as another, the
//! pairs of one group together: an array of arrays is not one of the types
//! wasm-bindgen carries, so the package cuts them, as section 11 of
//! `docs/architecture.md` has it for a table that crosses with a copy.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::pop_dists::{
    JackknifeGroups, PopDistMeasure, PopDistOptions, PopDistSums, calc_pop_dist_sums,
};
use popnei::stats::Pops;

use crate::errors::JsPopneiError;
use crate::source::{LARGEST_POSITION, OpenSource, PassCounts, positions_of};
use crate::stats::pops_of_the_arrays;
use crate::steps::{Steps, chain_of};

/// The name of the argument that says how the variants are cut into the
/// resampling groups, as a TypeScript user writes it.
const JACKKNIFE_GROUP: &str = "jackknifeGroup";

/// What a user writes in that argument for each variant to be a group of its
/// own.
const PER_VARIANT: &str = "variant";

/// The arguments of one pass, as they crossed from TypeScript.
///
/// The package has checked that each of them is of the type the core takes,
/// since a number of JavaScript reaches a whole number of the core as 32
/// bits with no error; what is left is what the core says of them, a name
/// that is of none of the seven measures among it.
pub(crate) struct ArgumentsOfTheDists {
    /// The name of each population, in the order the user gave them, which
    /// is the order of the pairs of every array of the result.
    pub(crate) pop_names: Vec<String>,
    /// The names of the individuals of every population, the ones of the
    /// first population first.
    pub(crate) pop_individuals: Vec<String>,
    /// How many individuals each population of `pop_names` holds, which cuts
    /// `pop_individuals` into the names of each of them.
    pub(crate) num_individuals_per_pop: Vec<u32>,
    /// The measures to calculate, under the names of
    /// [`PopDistMeasure::NAMES`].
    pub(crate) measures: Vec<String>,
    /// Whether each variant is a resampling group of its own, which is what
    /// the user wrote `"variant"` for.
    pub(crate) group_per_variant: bool,
    /// How many base pairs of one chromosome a resampling group holds, and
    /// nothing when the user asked for no groups or for one per variant.
    pub(crate) group_base_pairs: Option<f64>,
    /// How many called genotypes a population needs at a variant for that
    /// variant to count for a pair the population is in.
    pub(crate) min_num_individuals: u32,
}

/// What one pass gives TypeScript: the names of the populations in their
/// order, the measures that were asked for with their standard errors, how
/// many variants counted for each pair, the f_2 of every pair within every
/// resampling group, the chromosome and the two positions of each group, and
/// the counts of the pass.
///
/// Every array leaves the memory of wasm the first time it is asked for, and
/// the call after that gives nothing: the package reads each of them once,
/// into the object a user holds, and frees this.
#[wasm_bindgen]
pub struct PopDistsOfAPass {
    pop_names: Option<Vec<String>>,
    num_pairs: usize,
    values: Option<Vec<f64>>,
    standard_errors: Option<Vec<f64>>,
    num_vars: Option<Vec<i32>>,
    num_groups: usize,
    f2_groups: Option<Vec<f64>>,
    group_chroms: Option<Vec<String>>,
    group_starts: Option<Vec<f64>>,
    group_ends: Option<Vec<f64>>,
    counts: PassCounts,
}

#[wasm_bindgen]
impl PopDistsOfAPass {
    /// The names of the populations, in the order the user named them, which
    /// is the order the pairs of every array below are of.
    pub fn pop_names(&mut self) -> Option<Vec<String>> {
        self.pop_names.take()
    }

    /// How many pairs those populations make, which is how many values each
    /// measure holds.
    #[must_use]
    pub fn num_pairs(&self) -> usize {
        self.num_pairs
    }

    /// The value of every pair of every measure that was asked for, the
    /// pairs of one measure after the pairs of the one before it, in the
    /// order the measures were asked for, and NaN for a pair that counted no
    /// variant.
    pub fn values(&mut self) -> Option<Vec<f64>> {
        self.values.take()
    }

    /// The standard error beside each of those values, in the same order, or
    /// `undefined` when no resampling groups were asked for. A pair that has
    /// a value and no standard error, one whose variants all fell in one
    /// group, is NaN.
    pub fn standard_errors(&mut self) -> Option<Vec<f64>> {
        self.standard_errors.take()
    }

    /// How many variants counted for each pair.
    pub fn num_vars(&mut self) -> Option<Vec<i32>> {
        self.num_vars.take()
    }

    /// How many resampling groups the variants fell into, which is 0 when
    /// none were asked for.
    #[must_use]
    pub fn num_groups(&self) -> usize {
        self.num_groups
    }

    /// The f_2 of every pair within every group, the pairs of one group
    /// together, a table of `num_groups` x `num_pairs`, or `undefined` when
    /// no groups were asked for.
    pub fn f2_groups(&mut self) -> Option<Vec<f64>> {
        self.f2_groups.take()
    }

    /// The name of the chromosome of each group, in the order the groups
    /// were started.
    pub fn group_chroms(&mut self) -> Option<Vec<String>> {
        self.group_chroms.take()
    }

    /// The position of the first variant of each group, 1 based as in a VCF.
    pub fn group_starts(&mut self) -> Option<Vec<f64>> {
        self.group_starts.take()
    }

    /// The position of the last variant of each group, that one included.
    pub fn group_ends(&mut self) -> Option<Vec<f64>> {
        self.group_ends.take()
    }

    /// How many variants the pass gave, and what each filter of it was given
    /// and kept.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.counts.clone()
    }
}

/// Every measure of `asked.measures` for every pair of the populations of
/// `asked`, over one pass of `source` through the steps of `steps`.
///
/// The chain of readers of the pass stays here, lent to the core, so that
/// the counts of its filters can be read when the calculation returns: the
/// loop over the blocks is the core's, and no block of it reaches this
/// crate. The reader is opened at the size of its own blocks, since every
/// measure is a ratio of sums that are added over them and the same numbers
/// come out whatever the size.
///
/// # Errors
///
/// When a name of `asked.measures` is of none of the seven measures; when
/// the resampling groups were asked for as a length that is no whole number
/// of base pairs of 1 or more; when a population names an individual the
/// pass does not give, names one twice or names none; when the populations
/// are fewer than two; when the pass gives no variant; when the variants
/// fall into fewer resampling groups than a standard error is built from;
/// when the memory of the tab does not take the six sums of every pair and
/// group; when a group has no name for its chromosome or holds a position
/// above 2^53, which JavaScript does not hold; when a pair counted more
/// variants than a JavaScript array of counts holds; and when the source
/// cannot be read, a wrong line of a VCF among the causes.
pub(crate) fn pop_dists_of(
    source: &dyn OpenSource,
    steps: &Steps,
    asked: &ArgumentsOfTheDists,
) -> Result<PopDistsOfAPass, JsPopneiError> {
    let measures = the_measures(&asked.measures)?;
    let options = PopDistOptions {
        min_num_individuals: asked.min_num_individuals,
        groups: the_jackknife_group(asked.group_per_variant, asked.group_base_pairs)?,
    };
    let named = pops_of_the_arrays(
        &asked.pop_names,
        &asked.pop_individuals,
        &asked.num_individuals_per_pop,
    )?;
    let reader = source.reader(None)?;
    let mut chain = chain_of(reader, steps.steps())?;
    let pops = Pops::from_names(&named, chain.individuals())?;
    let pop_names = (0..pops.len())
        .map(|pop| pops.name(pop).to_owned())
        .collect();
    let sums = calc_pop_dist_sums(&mut *chain, &pops, &options)?;
    let counts = PassCounts::of(sums.num_vars(), &chain.filtering_stats());
    // The pairs in the order of the distance vector, (0, 1), (0, 2), ...,
    // (1, 2), ..., which is the order of every array of the result.
    let pairs: Vec<(usize, usize)> = (0..sums.num_pops())
        .flat_map(|first| {
            (first..sums.num_pops())
                .skip(1)
                .map(move |second| (first, second))
        })
        .collect();
    let groups_were_asked_for = options.groups != JackknifeGroups::None;
    let mut values = Vec::with_capacity(measures.len().saturating_mul(pairs.len()));
    let mut standard_errors = Vec::new();
    for measure in &measures {
        values.extend(
            sums.measures(*measure)
                .map(|value| value.unwrap_or(f64::NAN)),
        );
        if groups_were_asked_for {
            standard_errors.extend(pairs.iter().map(|(first, second)| {
                sums.standard_error(*measure, *first, *second)
                    .unwrap_or(f64::NAN)
            }));
        }
    }
    let num_vars = pairs
        .iter()
        .map(|(first, second)| for_javascript(sums.num_vars_of(*first, *second).unwrap_or(0)))
        .collect::<Result<Vec<i32>, JsPopneiError>>()?;
    // The core gives no group at all when no standard errors were asked
    // for, so this is 0 there and the three arrays below are empty.
    let num_groups = sums.groups().len();
    let f2_groups = groups_were_asked_for.then(|| f2_of_every_group(&sums, &pairs));
    let group_chroms = named_chroms(chain.chroms(), &sums)?;
    let group_starts = positions_of(sums.groups().iter().map(|group| group.start).collect())?;
    let group_ends = positions_of(sums.groups().iter().map(|group| group.end).collect())?;
    Ok(PopDistsOfAPass {
        pop_names: Some(pop_names),
        num_pairs: pairs.len(),
        values: Some(values),
        standard_errors: groups_were_asked_for.then_some(standard_errors),
        num_vars: Some(num_vars),
        num_groups,
        f2_groups,
        group_chroms: Some(group_chroms),
        group_starts: Some(group_starts),
        group_ends: Some(group_ends),
        counts,
    })
}

/// The f_2 of every pair within every group, the pairs of one group
/// together: a table of groups x pairs that f_3 and f_4 are built from later
/// without reading the genotypes again.
fn f2_of_every_group(sums: &PopDistSums, pairs: &[(usize, usize)]) -> Vec<f64> {
    (0..sums.groups().len())
        .flat_map(|group| {
            pairs.iter().map(move |(first, second)| {
                sums.f2_of_group(group, *first, *second).unwrap_or(f64::NAN)
            })
        })
        .collect()
}

/// The name of the chromosome of every group, in the order the groups were
/// started.
///
/// # Errors
///
/// When the table of the reader holds no name for the number a group
/// carries, which is a defect of popnei.
fn named_chroms(
    chroms: &popnei::variant::ChromTable,
    sums: &PopDistSums,
) -> Result<Vec<String>, JsPopneiError> {
    sums.groups()
        .iter()
        .map(|group| {
            chroms.name(group.chrom).map(str::to_owned).ok_or_else(|| {
                JsPopneiError::Broken(format!(
                    "the chromosome number {number} of a resampling group is not in \
                     the table of the reader that gave it",
                    number = group.chrom
                ))
            })
        })
        .collect()
}

/// `count` as the array of counts of JavaScript holds it.
///
/// The variants that counted for a pair cross as an `Int32Array`, which is
/// what "Its Python function" of `docs/specs/dists.md` gives the result, and
/// the core counts in 64 bits.
///
/// # Errors
///
/// When the count is above 2147483647, which is more variants than a file in
/// the memory of a tab holds: that memory addresses 4 GB, and a variant is a
/// row of a file.
fn for_javascript(count: u64) -> Result<i32, JsPopneiError> {
    i32::try_from(count).map_err(|_| {
        JsPopneiError::NotInJavaScript(format!(
            "{count} variants counted for one pair of populations, more than a \
             JavaScript array of counts holds"
        ))
    })
}

/// The measures a user asked for, out of the names the TypeScript package
/// gives them, which are the core's.
///
/// # Errors
///
/// A name that is of none of the seven, which a user reaches by writing one
/// in JavaScript: in TypeScript the seven are a union of string literals.
fn the_measures(names: &[String]) -> Result<Vec<PopDistMeasure>, JsPopneiError> {
    let mut asked_for = Vec::with_capacity(names.len());
    for name in names {
        asked_for.push(PopDistMeasure::of_name(name)?);
    }
    Ok(asked_for)
}

/// How the variants are cut into the resampling groups, out of the `null`,
/// the `"variant"` or the length in base pairs a user wrote.
///
/// The package tells the three kinds apart, since it is the one that sees
/// what a user wrote, and the value of a length is checked here: a number of
/// JavaScript is a float64, and what the core takes is a whole number of
/// base pairs of 1 or more.
///
/// # Errors
///
/// A length of 0 base pairs, which is no stretch of a chromosome, a negative
/// one, one with a fraction and one above 2^53, the largest whole number a
/// number of JavaScript holds. And a length given beside `"variant"`, which
/// is a defect of the package: they are the two kinds of the one argument.
fn the_jackknife_group(
    per_variant: bool,
    base_pairs: Option<f64>,
) -> Result<JackknifeGroups, JsPopneiError> {
    match (per_variant, base_pairs) {
        (true, Some(length)) => Err(JsPopneiError::Broken(format!(
            "the pass was asked for a resampling group of each variant and for \
             groups of {length} base pairs at once"
        ))),
        (true, None) => Ok(JackknifeGroups::PerVariant),
        (false, None) => Ok(JackknifeGroups::None),
        (false, Some(length)) => {
            #[expect(
                clippy::cast_precision_loss,
                reason = "2^53 is exact as a float64, which is what it is the largest of"
            )]
            let largest = LARGEST_POSITION as f64;
            if !length.is_finite() || length < 1.0 || length > largest || length.fract() != 0.0 {
                return Err(JsPopneiError::Refused(no_length(length)));
            }
            #[expect(
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss,
                reason = "checked above to be a whole number between 1 and 2^53"
            )]
            let length = length as u64;
            Ok(JackknifeGroups::OfBasePairs(length))
        }
    }
}

/// What a `jackknifeGroup` that is no length is told, which names the
/// argument, what the user wrote and the three kinds it takes.
fn no_length(given: f64) -> String {
    format!(
        "`{JACKKNIFE_GROUP}` is {given}, and it says how the variants are cut into the \
         resampling groups the standard errors are built from: a whole length in base \
         pairs from 1 to {LARGEST_POSITION}, `\"{PER_VARIANT}\"` for a group of each \
         variant, or `null` for no standard error"
    )
}
