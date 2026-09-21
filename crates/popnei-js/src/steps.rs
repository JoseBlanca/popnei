//! The steps of a `Variants`, and the chain of readers a pass builds from
//! them.
//!
//! A step is what every pass built from a `Variants` does to its variants
//! after they are read: a filter of `docs/specs/filters.md` is the only kind
//! there will be. The list lives here, in an object of this crate that the
//! `Variants` of the package holds, because this crate is what reads it: it
//! builds the chain of readers of each pass from the steps the `Variants`
//! has when the pass starts, and, when a user adds a filter, it looks for
//! the kind among the steps and refuses a second filter of a kind that is
//! already there, which "The Rust interface" of `docs/specs/filters.md` asks
//! of it.
//!
//! No filter is built yet, so [`Step`] has no case, the list of a `Variants`
//! is always empty and the chain of a pass is its source alone. Task 3.5 of
//! `docs/plans/filters.md` adds the filter to both.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::block::BlockReader;

use crate::errors::JsPopneiError;

/// One step of a `Variants`.
///
/// The filters of `docs/specs/filters.md` are the only kind of step, and
/// none of them is built, so this enum has no case yet and no step can be
/// made: the filter, which holds the criterion of the core and the threshold
/// of it, is task 3.5 of `docs/plans/filters.md`.
#[derive(Clone, Copy)]
pub(crate) enum Step {}

impl Step {
    /// The kind of the step, which is the name its counts have for a
    /// TypeScript user, `"missing_data"`, `"maf"` or `"obs_het"`.
    fn kind(self) -> &'static str {
        match self {}
    }

    /// The arguments of the step, each with the name a TypeScript user
    /// writes for it, `maxAllowedMaf`.
    fn args(self) -> Vec<(&'static str, f64)> {
        match self {}
    }
}

/// The steps of one `Variants`, in the order in which they were put on it,
/// or the copy of that list that one pass runs.
///
/// A `Variants` of the package holds one of these and gives every pass it
/// starts the copy that [`Steps::of_a_pass`] makes, which that pass takes
/// over and frees: a step added while a pass runs holds from the next one,
/// as `docs/specs/filters.md` says, and a copy is what wasm-bindgen carries,
/// which takes an argument of this type by value and not by reference.
///
/// What a user reads of it are the four arrays below, which the package puts
/// together into the steps of `docs/specs/filters.md`: an array of objects
/// is not one of the types wasm-bindgen carries, so the arguments of every
/// step cross flat, as the alleles of a block do, with how many of them each
/// step has beside them.
#[wasm_bindgen]
pub struct Steps {
    steps: Vec<Step>,
}

#[wasm_bindgen]
impl Steps {
    /// The steps of a `Variants` that nothing has been put on yet, which is
    /// what `openVcf` and `openVars` give.
    #[wasm_bindgen(constructor)]
    #[must_use]
    pub fn new() -> Steps {
        Steps { steps: Vec::new() }
    }

    /// The steps as they are now, which is what one pass runs.
    ///
    /// The pass it is given to takes it over and frees it when it has built
    /// its chain of readers, so the list a `Variants` holds is not the one
    /// that crosses and a step added while a pass runs reaches no pass that
    /// started before it.
    #[must_use]
    pub fn of_a_pass(&self) -> Steps {
        Steps {
            steps: self.steps.clone(),
        }
    }

    /// The kind of each step, in the order of the steps.
    #[must_use]
    pub fn kinds(&self) -> Vec<String> {
        self.steps
            .iter()
            .map(|step| step.kind().to_owned())
            .collect()
    }

    /// The name of every argument, the arguments of the first step first.
    #[must_use]
    pub fn arg_names(&self) -> Vec<String> {
        self.args()
            .into_iter()
            .map(|(name, _value)| name.to_owned())
            .collect()
    }

    /// The value of every argument, each beside its name in `arg_names`.
    #[must_use]
    pub fn arg_values(&self) -> Vec<f64> {
        self.args()
            .into_iter()
            .map(|(_name, value)| value)
            .collect()
    }

    /// How many arguments each step has, which is what cuts `arg_names` and
    /// `arg_values` into the arguments of each step.
    ///
    /// # Errors
    ///
    /// When a step has more arguments than a JavaScript array of counts
    /// holds, which no step of this crate has: a filter takes one threshold.
    pub fn num_args_per_step(&self) -> Result<Vec<u32>, JsPopneiError> {
        self.steps
            .iter()
            .map(|step| {
                let num_args = step.args().len();
                u32::try_from(num_args).map_err(|_| {
                    JsPopneiError::Broken(format!(
                        "the step `{kind}` of these variants has {num_args} arguments, \
                         more than a JavaScript array of counts holds",
                        kind = step.kind()
                    ))
                })
            })
            .collect()
    }
}

impl Default for Steps {
    fn default() -> Steps {
        Steps::new()
    }
}

impl Steps {
    /// The steps themselves, which is what the chain of a pass is built
    /// from.
    pub(crate) fn steps(&self) -> &[Step] {
        &self.steps
    }

    /// Every argument of every step, the arguments of the first step first.
    fn args(&self) -> Vec<(&'static str, f64)> {
        self.steps.iter().flat_map(|step| step.args()).collect()
    }
}

/// The chain of readers of one pass: `reader`, the source, with a filter
/// over it for each step, in the order of the steps, so that each filter
/// sees only what the one before it kept.
///
/// # Errors
///
/// When a filter cannot be put over the chain, which task 3.5 of
/// `docs/plans/filters.md` brings with the filters. There is no step yet, so
/// today it gives the source back as it is.
pub(crate) fn chain_of(
    reader: Box<dyn BlockReader>,
    steps: &[Step],
) -> Result<Box<dyn BlockReader>, popnei::Error> {
    // A `Step` has no case, so there is nothing to put over the source and
    // the chain is the source alone. This says it to the compiler and not in
    // a comment: the day a step exists, this match is not exhaustive any
    // more and the crate does not build until the reader of that step is put
    // here, over the chain, one for each step and in their order, so that
    // each filter sees only what the one before it kept.
    if let Some(step) = steps.first() {
        match *step {}
    }
    Ok(reader)
}
