//! The steps of a `Variants`, and the chain of readers a pass builds from
//! them.
//!
//! A step is what every pass built from a `Variants` does to its variants
//! after they are read: a filter of `docs/specs/filters.md` is the only kind
//! there is. The list lives here, in an object of this crate that the
//! `Variants` of the package holds, because this crate is what reads it: it
//! builds the chain of readers of each pass from the steps the `Variants`
//! has when the pass starts, and, when a user adds a filter, it looks for
//! the kind among the steps and refuses a second filter of a kind that is
//! already there, which "The Rust interface" of `docs/specs/filters.md` asks
//! of it.
//!
//! The three methods that add a filter are here as well, one for each of the
//! three numbers of a variant a filter compares, and each of them refuses at
//! the call what a user cannot filter by: a threshold that is not a number
//! from 0 to 1, under the name of the argument they wrote it in, and a
//! second filter of a kind the list holds, with the threshold of the one
//! that is set. No reader exists at that call, so neither refusal can come
//! from the chain.

use std::sync::{Mutex, MutexGuard};

use pyo3::prelude::*;

use popnei::block::BlockReader;
use popnei::filters::{VarFilter, VarFilteringCriterion};

use crate::errors::PyPopneiError;
use crate::source::threshold_of;

/// One step of a `Variants`.
///
/// A filter of `docs/specs/filters.md` is the only kind of step there is: it
/// holds which number of a variant is compared and the threshold it is
/// compared with, which is what the core needs to build the filter of a
/// pass.
#[derive(Clone, Copy)]
pub(crate) enum Step {
    /// The variants whose number is at most the threshold of the criterion
    /// are kept, and the others are left out of every block of the pass.
    Filter(VarFilteringCriterion),
}

impl Step {
    /// The kind of the step, which is the name its counts have for a Python
    /// user, `"missing_data"`, `"maf"` or `"obs_het"`.
    fn kind(self) -> &'static str {
        match self {
            Step::Filter(criterion) => criterion.kind(),
        }
    }

    /// The arguments of the step, each with the name a Python user writes
    /// for it, `max_allowed_maf`.
    fn args(self) -> Vec<(&'static str, f64)> {
        match self {
            Step::Filter(criterion) => vec![(argument_of(criterion), criterion.threshold())],
        }
    }

    /// The threshold a filter compares the number of a variant with, which
    /// the refusal of a second filter of its kind names.
    fn threshold(self) -> f64 {
        match self {
            Step::Filter(criterion) => criterion.threshold(),
        }
    }
}

/// The names a Python user writes the threshold of each filter under, which
/// are the arguments of the three methods that add one.
///
/// The core names a filter by its kind, `maf`, and knows nothing of the
/// arguments of Python, so the two names meet here: a user who is told that
/// `max_allowed_maf` is 1.5 reads the name they wrote.
const MAX_ALLOWED_MISSING_RATE: &str = "max_allowed_missing_rate";
const MAX_ALLOWED_MAF: &str = "max_allowed_maf";
const MAX_ALLOWED_OBS_HET: &str = "max_allowed_obs_het";

/// The name a Python user writes the threshold of `criterion` under.
fn argument_of(criterion: VarFilteringCriterion) -> &'static str {
    match criterion {
        VarFilteringCriterion::MaxMissingRate(_) => MAX_ALLOWED_MISSING_RATE,
        VarFilteringCriterion::MaxMaf(_) => MAX_ALLOWED_MAF,
        VarFilteringCriterion::MaxObsHet(_) => MAX_ALLOWED_OBS_HET,
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

    // The three filters, each with the argument its method of the package
    // takes: the missing genotypes of a variant divided by all the
    // individuals, the count of its commonest allele divided by its called
    // alleles, and its heterozygous genotypes divided by its called ones.
    // The threshold is taken as the object it is and converted here, where
    // what is refused names the argument the user wrote.
    fn filter_by_missing_data(
        &self,
        max_allowed_missing_rate: &Bound<'_, PyAny>,
    ) -> Result<(), PyPopneiError> {
        self.add(VarFilteringCriterion::MaxMissingRate(threshold_of(
            MAX_ALLOWED_MISSING_RATE,
            max_allowed_missing_rate,
        )?))
    }

    fn filter_by_maf(&self, max_allowed_maf: &Bound<'_, PyAny>) -> Result<(), PyPopneiError> {
        self.add(VarFilteringCriterion::MaxMaf(threshold_of(
            MAX_ALLOWED_MAF,
            max_allowed_maf,
        )?))
    }

    fn filter_by_obs_het(
        &self,
        max_allowed_obs_het: &Bound<'_, PyAny>,
    ) -> Result<(), PyPopneiError> {
        self.add(VarFilteringCriterion::MaxObsHet(threshold_of(
            MAX_ALLOWED_OBS_HET,
            max_allowed_obs_het,
        )?))
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
        Ok(self.locked()?.clone())
    }

    /// The filter of `criterion` added at the end of the list, where the
    /// next pass takes it.
    ///
    /// # Errors
    ///
    /// When the threshold is not a number from 0 to 1, which the core is
    /// what says: the filter built here is dropped, and every pass builds
    /// its own from the criterion, so the rule that a threshold has to keep
    /// is written in one place. And when the list holds a filter of the
    /// kind of `criterion` already, since two threshold filters of one kind
    /// keep the variants that the stricter of them keeps alone. After
    /// either, the list is as it was.
    fn add(&self, criterion: VarFilteringCriterion) -> Result<(), PyPopneiError> {
        VarFilter::new(criterion).map_err(|error| under_the_argument(error, criterion))?;
        let mut steps = self.locked()?;
        if let Some(set) = steps.iter().find(|step| step.kind() == criterion.kind()) {
            // The core leaves the threshold of the filter that is set out,
            // because a chain of readers says which kinds it holds and not
            // with which thresholds. Here the steps are at hand, so the
            // user reads the number they set as well as the one they wrote.
            return Err(PyPopneiError::Core(
                popnei::Error::VarFilterOfAKindThatIsSet {
                    kind: criterion.kind(),
                    threshold: criterion.threshold(),
                    threshold_that_is_set: Some(set.threshold()),
                },
            ));
        }
        steps.push(Step::Filter(criterion));
        Ok(())
    }

    /// The list, held for as long as it is read or changed.
    ///
    /// # Errors
    ///
    /// When a panic left the lock broken, which is a defect of this crate.
    fn locked(&self) -> Result<MutexGuard<'_, Vec<Step>>, PyPopneiError> {
        self.steps.lock().map_err(|_| PyPopneiError::Broken {
            message: "the steps of these variants cannot be read or changed any more: a \
                      panic left the list of them half way through a change"
                .to_string(),
            path: None,
        })
    }
}

/// `error`, and a threshold the core refused under the name of the argument
/// a Python user wrote it in.
///
/// The core names the filter by its kind, `maf`, and the message a user
/// reads names `max_allowed_maf`, which only this crate knows: what they
/// have to look at is the call they wrote.
fn under_the_argument(error: popnei::Error, criterion: VarFilteringCriterion) -> PyPopneiError {
    if let popnei::Error::VarFilterThresholdOutOfRange { threshold, .. } = error {
        return PyPopneiError::Threshold {
            name: argument_of(criterion),
            value: format!("{threshold:?}"),
        };
    }
    PyPopneiError::Core(error)
}

/// The chain of readers of one pass: `reader`, the source, with a filter
/// over it for each step, in the order of the steps, so that each filter
/// sees only what the one before it kept.
///
/// The chain itself is the core's, `popnei::filters::chain_of`, which both
/// binding crates call: what this one does is read the criterion of each
/// step, which is what the steps of this crate hold and the core does not
/// know. The filters belong to the pass this chain is built for, so no
/// count is shared with another pass, and the chain is asked for the fields
/// the consumer wants once it is built.
///
/// # Errors
///
/// When a threshold of the steps is not a number from 0 to 1, or when the
/// steps hold two filters of one kind. A user reaches neither: the call
/// that adds a filter refuses both, and the steps are read from there.
pub(crate) fn chain_of(
    reader: Box<dyn BlockReader>,
    steps: &[Step],
) -> Result<Box<dyn BlockReader>, popnei::Error> {
    let criteria: Vec<VarFilteringCriterion> = steps
        .iter()
        .map(|step| match *step {
            Step::Filter(criterion) => criterion,
        })
        .collect();
    popnei::filters::chain_of(reader, &criteria)
}
