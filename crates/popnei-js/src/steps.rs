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
//! The four methods that add a filter are here as well, one for each of the
//! three numbers of a variant a filter compares and one for the individuals
//! to keep, and each of them refuses at the call what a user cannot filter
//! by: a threshold that is not a number from 0 to 1, under the name of the
//! argument they wrote it in; a name that is not an individual of the
//! source, a name that is there twice and no name at all; and a second
//! filter of a kind the list holds, with the threshold of the one that is
//! set when both are threshold filters. No reader exists at that call, so
//! none of those refusals can come from the chain, and the individuals of
//! the source, which the names are resolved against, are held here from the
//! moment the `Variants` is built. What is not a number at all, the
//! `undefined` of a call with no threshold among it, is refused by the
//! package before the call, in `js/popnei/src/arguments.ts`: everything
//! arrives here as a float64, and `null` would arrive as a threshold of 0.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::block::BlockReader;
use popnei::filters::{
    PassStep, VarFilter, VarFilteringCriterion, refuse_a_second_filter_of_a_kind,
    resolve_individuals,
};

use crate::errors::JsPopneiError;

/// One step of a `Variants`: what the core does with it, and the arguments
/// a TypeScript user wrote it with.
///
/// What the pass does is the core's `PassStep`, which is what builds the
/// readers of the pass. Beside it this crate keeps the name each argument of
/// the step has in TypeScript, `maxAllowedMaf`, and the value under it: the
/// core names a filter by its kind, `maf`, and knows nothing of the
/// arguments of TypeScript, and a user who reads their steps reads the names
/// they wrote.
#[derive(Clone)]
pub(crate) struct Step {
    /// What every pass does to its variants at this step.
    pass_step: PassStep,
    /// The arguments of the step, each under the name a TypeScript user
    /// writes it in.
    args: Vec<(&'static str, Argument)>,
}

/// What a user gave one argument of a step: the threshold of a filter, a
/// number from 0 to 1, or the names of the individuals to keep, in the
/// order they named them.
///
/// It is the value that argument has in the `args` of the step a user
/// reads, a number for a threshold and an array of strings for the
/// individuals, which is what "In Python and in TypeScript" of
/// `docs/specs/filters.md` gives them.
#[derive(Clone)]
enum Argument {
    Threshold(f64),
    Individuals(Vec<String>),
}

/// The names a TypeScript user writes the argument of each filter under,
/// which are the arguments of the four methods that add one.
const MAX_ALLOWED_MISSING_RATE: &str = "maxAllowedMissingRate";
const MAX_ALLOWED_MAF: &str = "maxAllowedMaf";
const MAX_ALLOWED_OBS_HET: &str = "maxAllowedObsHet";
const INDIVIDUALS: &str = "individuals";

/// The steps of one `Variants`, in the order in which they were put on it,
/// or the copy of that list that one pass runs.
///
/// A `Variants` of the package holds one of these and gives every pass it
/// starts the copy that [`Steps::of_a_pass`] makes, which that pass takes
/// over and frees: a step added while a pass runs holds from the next one,
/// as `docs/specs/filters.md` says, and a copy is what wasm-bindgen carries,
/// which takes an argument of this type by value and not by reference.
///
/// What a user reads of it are the arrays below, which the package puts
/// together into the steps of `docs/specs/filters.md`: an array of objects
/// is not one of the types wasm-bindgen carries, so the arguments of every
/// step cross flat, as the alleles of a block do, with how many of them each
/// step has beside them. The value of an argument crosses in the array of
/// its kind, a threshold in [`Steps::arg_thresholds`] and the names of the
/// individuals in [`Steps::arg_individuals`], and
/// [`Steps::num_names_per_arg`] says which argument is which.
#[wasm_bindgen]
pub struct Steps {
    steps: Vec<Step>,
    /// The individuals of the source, in its order, which the names given
    /// to the filter of individuals are resolved against at the call.
    of_the_source: Vec<String>,
}

#[wasm_bindgen]
impl Steps {
    /// The steps of a `Variants` over a source of `individuals` that
    /// nothing has been put on yet, which is what `openVcf` and `openVars`
    /// give.
    #[wasm_bindgen(constructor)]
    #[must_use]
    pub fn new(individuals: Vec<String>) -> Steps {
        Steps {
            steps: Vec::new(),
            of_the_source: individuals,
        }
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
            of_the_source: self.of_the_source.clone(),
        }
    }

    /// The kind of each step, in the order of the steps.
    #[must_use]
    pub fn kinds(&self) -> Vec<String> {
        self.steps
            .iter()
            .map(|step| step.pass_step.kind().to_owned())
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

    /// The threshold of every argument that is one, in the order of
    /// `arg_names`, and nothing for an argument that is not.
    #[must_use]
    pub fn arg_thresholds(&self) -> Vec<f64> {
        self.args()
            .into_iter()
            .filter_map(|(_name, value)| match value {
                Argument::Threshold(threshold) => Some(threshold),
                Argument::Individuals(_) => None,
            })
            .collect()
    }

    /// The names of the individuals of every argument that holds names, one
    /// argument after another in the order of `arg_names`, and nothing for
    /// an argument that holds none.
    #[must_use]
    pub fn arg_individuals(&self) -> Vec<String> {
        self.args()
            .into_iter()
            .filter_map(|(_name, value)| match value {
                Argument::Threshold(_) => None,
                Argument::Individuals(names) => Some(names),
            })
            .flatten()
            .collect()
    }

    /// How many names of individuals each argument of `arg_names` holds,
    /// which cuts `arg_individuals` into the names of each of them, and 0
    /// for an argument whose value is a threshold, whose one number is in
    /// `arg_thresholds`.
    ///
    /// The 0 tells the two apart: a filter of individuals holds one name at
    /// least, since a filter of no individual is refused at the call that
    /// adds it.
    ///
    /// # Errors
    ///
    /// When an argument holds more names than a JavaScript array of counts
    /// does, which is more individuals than the memory of wasm takes the
    /// names of.
    pub fn num_names_per_arg(&self) -> Result<Vec<u32>, JsPopneiError> {
        self.args()
            .into_iter()
            .map(|(name, value)| match value {
                Argument::Threshold(_) => Ok(0),
                Argument::Individuals(names) => {
                    let num_names = names.len();
                    u32::try_from(num_names).map_err(|_| {
                        JsPopneiError::Broken(format!(
                            "the argument `{name}` of a step of these variants holds \
                             {num_names} names, more than a JavaScript array of counts \
                             holds"
                        ))
                    })
                }
            })
            .collect()
    }

    /// The variants whose missing genotypes divided by all the individuals
    /// are at most `max_allowed_missing_rate` are kept.
    ///
    /// # Errors
    ///
    /// When the threshold is not a number from 0 to 1, and when a filter of
    /// this kind is set already.
    pub fn filter_by_missing_data(
        &mut self,
        max_allowed_missing_rate: f64,
    ) -> Result<(), JsPopneiError> {
        self.add_a_threshold_filter(
            MAX_ALLOWED_MISSING_RATE,
            VarFilteringCriterion::MaxMissingRate(max_allowed_missing_rate),
        )
    }

    /// The variants whose commonest allele divided by their called alleles
    /// is at most `max_allowed_maf` are kept.
    ///
    /// # Errors
    ///
    /// The two of [`Steps::filter_by_missing_data`].
    pub fn filter_by_maf(&mut self, max_allowed_maf: f64) -> Result<(), JsPopneiError> {
        self.add_a_threshold_filter(
            MAX_ALLOWED_MAF,
            VarFilteringCriterion::MaxMaf(max_allowed_maf),
        )
    }

    /// The variants whose heterozygous genotypes divided by their called
    /// ones are at most `max_allowed_obs_het` are kept.
    ///
    /// # Errors
    ///
    /// The two of [`Steps::filter_by_missing_data`].
    pub fn filter_by_obs_het(&mut self, max_allowed_obs_het: f64) -> Result<(), JsPopneiError> {
        self.add_a_threshold_filter(
            MAX_ALLOWED_OBS_HET,
            VarFilteringCriterion::MaxObsHet(max_allowed_obs_het),
        )
    }

    /// The genotypes of `individuals` kept at every variant, in the order
    /// they are named here, and those of no other individual.
    ///
    /// The names are resolved against the individuals of the source, so a
    /// name that is not one of them, a name that is there twice and no name
    /// at all are refused at this call and not when a pass runs.
    ///
    /// # Errors
    ///
    /// Those three, and a filter of individuals that is set already.
    pub fn filter_individuals(&mut self, individuals: Vec<String>) -> Result<(), JsPopneiError> {
        // What the names give is dropped: every pass resolves them again
        // when it builds its chain, so the rule lives in the core alone.
        // They are refused before the list is looked at, since a name that
        // is of no individual is wrong whatever the list holds.
        resolve_individuals(&individuals, &self.of_the_source)?;
        let step = Step {
            pass_step: PassStep::KeepIndividuals(individuals.clone()),
            args: vec![(INDIVIDUALS, Argument::Individuals(individuals))],
        };
        refuse_a_second_filter_of_a_kind(&pass_steps_of(&self.steps), &step.pass_step)?;
        self.steps.push(step);
        Ok(())
    }

    /// How many arguments each step has, which is what cuts `arg_names` and
    /// the values of the arguments into the arguments of each step.
    ///
    /// # Errors
    ///
    /// When a step has more arguments than a JavaScript array of counts
    /// holds, which no step of this crate has: a filter takes one threshold.
    pub fn num_args_per_step(&self) -> Result<Vec<u32>, JsPopneiError> {
        self.steps
            .iter()
            .map(|step| {
                let num_args = step.args.len();
                u32::try_from(num_args).map_err(|_| {
                    JsPopneiError::Broken(format!(
                        "the step `{kind}` of these variants has {num_args} arguments, \
                         more than a JavaScript array of counts holds",
                        kind = step.pass_step.kind()
                    ))
                })
            })
            .collect()
    }
}

impl Steps {
    /// The steps themselves, which is what the chain of a pass is built
    /// from.
    pub(crate) fn steps(&self) -> &[Step] {
        &self.steps
    }

    /// The step of the filter of `criterion` added at the end of the list,
    /// where the next pass takes it, with its threshold under `argument`,
    /// the name the user wrote it in.
    ///
    /// # Errors
    ///
    /// When the threshold is not a number from 0 to 1, and when the list
    /// holds a filter of the kind of `criterion` already. The core is what
    /// says both: the filter built here is dropped, and every pass builds
    /// its own from the criterion, so the rule that a threshold has to keep
    /// and which filters can stand together are written in one place. The
    /// threshold is refused first, since it is wrong whatever the list
    /// holds. After either, the list is as it was.
    fn add_a_threshold_filter(
        &mut self,
        argument: &'static str,
        criterion: VarFilteringCriterion,
    ) -> Result<(), JsPopneiError> {
        VarFilter::new(criterion).map_err(|error| under_the_argument(error, argument))?;
        let step = Step {
            pass_step: PassStep::VarFilter(criterion),
            args: vec![(argument, Argument::Threshold(criterion.threshold()))],
        };
        refuse_a_second_filter_of_a_kind(&pass_steps_of(&self.steps), &step.pass_step)?;
        self.steps.push(step);
        Ok(())
    }

    /// Every argument of every step, the arguments of the first step first.
    fn args(&self) -> Vec<(&'static str, Argument)> {
        self.steps
            .iter()
            .flat_map(|step| step.args.iter().cloned())
            .collect()
    }
}

/// `error`, and a threshold the core refused under `argument`, the name a
/// TypeScript user wrote it in.
///
/// The core names the filter by its kind, `maf`, and the message a user
/// reads names `maxAllowedMaf`, which only this crate knows: what they have
/// to look at is the call they wrote.
fn under_the_argument(error: popnei::Error, argument: &'static str) -> JsPopneiError {
    if let popnei::Error::VarFilterThresholdOutOfRange { threshold, .. } = error {
        return JsPopneiError::Threshold {
            name: argument,
            threshold,
        };
    }
    JsPopneiError::Core(error)
}

/// The chain of readers of one pass: `reader`, the source, with a filter
/// over it for each step, in the order of the steps, so that each filter
/// sees only what the one before it kept.
///
/// The chain itself is the core's, `popnei::filters::chain_of`, which both
/// binding crates call: what this one does is hand it the steps of the core
/// out of the steps of this crate, which carry the names of their arguments
/// in TypeScript beside them. The readers belong to the pass this chain is
/// built for, so no count is shared with another pass, and the chain is
/// asked for the fields the consumer wants once it is built.
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
    popnei::filters::chain_of(reader, &pass_steps_of(steps))
}

/// What the pass does at each step, in the order of the steps: what the core
/// is given, out of what the steps of this crate hold.
fn pass_steps_of(steps: &[Step]) -> Vec<PassStep> {
    steps.iter().map(|step| step.pass_step.clone()).collect()
}
