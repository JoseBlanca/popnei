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
//! is always empty and the chain of a pass is its source alone. Task 3.4 of
//! `docs/plans/filters.md` adds the filter to both.

use std::sync::Mutex;

use pyo3::prelude::*;

use popnei::block::BlockReader;

use crate::errors::PyPopneiError;

/// One step of a `Variants`.
///
/// The filters of `docs/specs/filters.md` are the only kind of step, and
/// none of them is built, so this enum has no case yet and no step can be
/// made: the filter, which holds the criterion of the core and the
/// threshold of it, is task 3.4 of `docs/plans/filters.md`.
#[derive(Clone, Copy)]
pub(crate) enum Step {}

impl Step {
    /// The kind of the step, which is the name its counts have for a Python
    /// user, `"missing_data"`, `"maf"` or `"obs_het"`.
    fn kind(self) -> &'static str {
        match self {}
    }

    /// The arguments of the step, each with the name a Python user writes
    /// for it, `max_allowed_maf`.
    fn args(self) -> Vec<(&'static str, f64)> {
        match self {}
    }
}

/// The kind of one step and its arguments on their way to Python, which the
/// package puts in the `Step` of `docs/specs/filters.md`.
type StepOfAVariants = (&'static str, Vec<(&'static str, f64)>);

// The steps of one `Variants`, in the order in which they were put on it. A
// `///` here would become the `__doc__` of the class, and what a Python user
// reads belongs to the package, which is the API.
//
// It is frozen with its list behind a `Mutex`, as every class of this crate
// is: a step can be added at any time, also between two passes, and any
// thread may hold the `Variants`.
#[pyclass(frozen, module = "popnei._core")]
pub(crate) struct Steps {
    steps: Mutex<Vec<Step>>,
}

#[pymethods]
impl Steps {
    // The steps of a `Variants` that nothing has been put on yet, which is
    // what `open_vcf` and `open_vars` give.
    #[new]
    fn new() -> Steps {
        Steps {
            steps: Mutex::new(Vec::new()),
        }
    }

    // The kind of each step and its arguments, in the order of the steps,
    // which the package gives a user as its `Step` objects.
    fn steps(&self) -> Result<Vec<StepOfAVariants>, PyPopneiError> {
        Ok(self
            .of_a_pass()?
            .into_iter()
            .map(|step| (step.kind(), step.args()))
            .collect())
    }
}

impl Steps {
    /// The steps as they are now, which is what one pass runs: a step added
    /// while a pass runs holds from the next one, as
    /// `docs/specs/filters.md` says, so every pass takes a copy of the list
    /// when it starts.
    ///
    /// # Errors
    ///
    /// When a panic left the lock broken, which is a defect of this crate.
    pub(crate) fn of_a_pass(&self) -> Result<Vec<Step>, PyPopneiError> {
        let steps = self.steps.lock().map_err(|_| PyPopneiError::Broken {
            message: "the steps of these variants cannot be read any more: a panic left \
                      the list of them half way through a change"
                .to_string(),
            path: None,
        })?;
        Ok(steps.clone())
    }
}

/// The chain of readers of one pass: `reader`, the source, with a filter
/// over it for each step, in the order of the steps, so that each filter
/// sees only what the one before it kept.
///
/// # Errors
///
/// When a filter cannot be put over the chain, which task 3.4 of
/// `docs/plans/filters.md` brings with the filters. There is no step yet,
/// so today it gives the source back as it is.
pub(crate) fn chain_of(
    reader: Box<dyn BlockReader>,
    steps: &[Step],
) -> Result<Box<dyn BlockReader>, popnei::Error> {
    // A `Step` has no case, so there is nothing to put over the source and
    // the chain is the source alone. This says it to the compiler and not
    // in a comment: the day a step exists, this match is not exhaustive any
    // more and the crate does not build until the reader of that step is
    // put here, over the chain, one for each step and in their order, so
    // that each filter sees only what the one before it kept.
    if let Some(step) = steps.first() {
        match *step {}
    }
    Ok(reader)
}
