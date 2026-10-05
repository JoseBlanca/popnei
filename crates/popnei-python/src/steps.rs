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
//! The eight methods that add a filter are here as well, one for each of the
//! three numbers of a variant a filter compares, one for the filter by
//! linkage disequilibrium, one for the individuals to keep, one for the
//! regions of a BED file, one for the first n variants and one for the
//! variants kept at random, and each of them refuses at the call what a
//! user cannot filter by: a threshold that is not a number from 0 to 1,
//! under the name of the argument they wrote it in; a keep rate that is not
//! a number from 0 to 1, and a seed that is not a whole number of 64 bits; a
//! window of fewer than 1 base pairs; a name that is not an individual of
//! the source, a name that is there twice and no name at all; a BED file
//! that cannot be read or that holds a line that is not a region, or no
//! region; a first n of no variant; a second filter of a kind the list
//! holds, with the threshold of the one that is set when both are threshold
//! filters; and a filter that takes variants out added after the filter of
//! the first n, which would leave fewer than the n it keeps. No reader
//! exists at that call, so none of those refusals can come from the chain,
//! and the individuals of the source, which the names are resolved against,
//! are held here from the moment the `Variants` is built.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::sync::{Arc, Mutex, MutexGuard};

use pyo3::prelude::*;
use pyo3::types::{PyFloat, PyTuple};

use popnei::block::BlockReader;
use popnei::filters::{
    LdFilter, PassStep, RandomFilter, RegionSelection, Regions, VarFilter, VarFilteringCriterion,
    first_n_step, individuals_of, refuse_a_step, resolve_individuals,
};

use crate::errors::PyPopneiError;
use crate::source::{count_of_at_least, distance_of, keep_rate_of, seed_of, threshold_of};

/// One step of a `Variants`: what the core does with it, and the arguments
/// a Python user wrote it with.
///
/// What the pass does is the core's `PassStep`, which is what builds the
/// readers of the pass. Beside it this crate keeps the name each argument of
/// the step has in Python, `max_allowed_maf`, and the value under it: the
/// core names a filter by its kind, `maf`, and knows nothing of the
/// arguments of Python, and a user who reads their steps reads the names
/// they wrote.
#[derive(Clone)]
pub(crate) struct Step {
    /// What every pass does to its variants at this step.
    pass_step: PassStep,
    /// The arguments of the step, each under the name a Python user writes
    /// it in.
    args: Vec<(&'static str, Argument)>,
}

/// What a user gave one argument of a step: the threshold of a filter, a
/// number from 0 to 1, the window of the filter by linkage disequilibrium,
/// a whole number of base pairs, the names of the individuals to keep, in
/// the order they named them, the path of a BED file, or the number of
/// variants the filter of the first n keeps; or what the step found in what
/// it was given, the number of regions of that file once those that overlap
/// or touch are joined; or the keep rate and the seed of the filter that
/// keeps variants at random.
///
/// It goes to Python as the value of that argument, a float for a
/// threshold and for a keep rate, an `int` for the window, for the two
/// counts and for the seed, a tuple of strings for the individuals and a
/// string for the path, which is
/// what "In Python and in TypeScript" of `docs/specs/filters.md` gives the
/// `args` of each step. A window is a whole number of base pairs and is not
/// a rate, so a user who wrote 10000 reads 10000 back and not `10000.0`.
#[derive(Clone)]
enum Argument {
    Threshold(f64),
    Distance(u64),
    Individuals(Vec<String>),
    Path(String),
    Count(usize),
    KeepRate(f64),
    Seed(u64),
}

impl<'py> IntoPyObject<'py> for Argument {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        match self {
            Argument::Threshold(threshold) => Ok(PyFloat::new(py, threshold).into_any()),
            Argument::Distance(distance) => Ok(distance.into_pyobject(py)?.into_any()),
            Argument::Individuals(names) => Ok(PyTuple::new(py, names)?.into_any()),
            Argument::Path(path) => Ok(path.into_pyobject(py)?.into_any()),
            Argument::Count(count) => Ok(count.into_pyobject(py)?.into_any()),
            Argument::KeepRate(keep_rate) => Ok(PyFloat::new(py, keep_rate).into_any()),
            Argument::Seed(seed) => Ok(seed.into_pyobject(py)?.into_any()),
        }
    }
}

/// The names a Python user writes the argument of each filter under, which
/// are the arguments of the methods that add one.
///
/// The core names a filter by its kind, `maf`, and knows nothing of the
/// arguments of Python, so the two names meet here: a user who is told that
/// `max_allowed_maf` is 1.5 reads the name they wrote.
const MAX_ALLOWED_MISSING_RATE: &str = "max_allowed_missing_rate";
const MAX_ALLOWED_MAF: &str = "max_allowed_maf";
const MAX_ALLOWED_OBS_HET: &str = "max_allowed_obs_het";
const MAX_ALLOWED_R2: &str = "max_allowed_r2";
const MAX_DIST: &str = "max_dist";
const INDIVIDUALS: &str = "individuals";
const BED_PATH: &str = "bed_path";
const NUM_REGIONS: &str = "num_regions";
const NUM_VARS: &str = "num_vars";
const KEEP_RATE: &str = "keep_rate";
const SEED: &str = "seed";

/// The kind of one step and its arguments on their way to Python, which the
/// package puts in the `Step` of `docs/specs/filters.md`.
type StepOfAVariants = (&'static str, Vec<(&'static str, Argument)>);

// The steps of one `Variants`, in the order in which they were put on it,
// and the individuals of its source. A `///` here would become the `__doc__`
// of the class, and what a Python user reads belongs to the package, which
// is the API.
//
// It is frozen with its list behind a `Mutex`, as every class of this crate
// is: a step can be added at any time, also between two passes, and any
// thread may hold the `Variants`. The individuals of the source are read
// from the header once and never change, so they are behind no lock.
#[pyclass(frozen, module = "popnei._core")]
pub(crate) struct Steps {
    steps: Mutex<Vec<Step>>,
    /// The individuals of the source, in its order, which the names given
    /// to the filter of individuals are resolved against at the call.
    of_the_source: Vec<String>,
}

#[pymethods]
impl Steps {
    // The steps of a `Variants` over a source of `individuals` that nothing
    // has been put on yet, which is what `open_vcf` and `open_vars` give.
    #[new]
    fn new(individuals: Vec<String>) -> Steps {
        Steps {
            steps: Mutex::new(Vec::new()),
            of_the_source: individuals,
        }
    }

    // The kind of each step and its arguments, in the order of the steps,
    // which the package gives a user as its `Step` objects.
    fn steps(&self) -> Result<Vec<StepOfAVariants>, PyPopneiError> {
        Ok(self
            .of_a_pass()?
            .into_iter()
            .map(|step| (step.pass_step.kind(), step.args))
            .collect())
    }

    // The names of the individuals the next pass gives, in its order: the
    // ones a filter of individuals among the steps keeps, and those of the
    // source when no step is that filter. Which step says it is the core's
    // rule, `individuals_of`, which the TypeScript crate reads as well.
    fn individuals(&self) -> Result<Vec<String>, PyPopneiError> {
        let steps = pass_steps_of(&self.locked()?);
        Ok(individuals_of(&steps, &self.of_the_source))
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
        self.add_a_threshold_filter(
            MAX_ALLOWED_MISSING_RATE,
            VarFilteringCriterion::MaxMissingRate(threshold_of(
                MAX_ALLOWED_MISSING_RATE,
                max_allowed_missing_rate,
            )?),
        )
    }

    fn filter_by_maf(&self, max_allowed_maf: &Bound<'_, PyAny>) -> Result<(), PyPopneiError> {
        self.add_a_threshold_filter(
            MAX_ALLOWED_MAF,
            VarFilteringCriterion::MaxMaf(threshold_of(MAX_ALLOWED_MAF, max_allowed_maf)?),
        )
    }

    fn filter_by_obs_het(
        &self,
        max_allowed_obs_het: &Bound<'_, PyAny>,
    ) -> Result<(), PyPopneiError> {
        self.add_a_threshold_filter(
            MAX_ALLOWED_OBS_HET,
            VarFilteringCriterion::MaxObsHet(threshold_of(
                MAX_ALLOWED_OBS_HET,
                max_allowed_obs_het,
            )?),
        )
    }

    // The filter by linkage disequilibrium, which takes a second argument:
    // how many base pairs behind a variant, on its chromosome, the variants
    // it is compared with are. Both are taken as the objects they are and
    // converted here, where what is refused names the argument the user
    // wrote; the window is a whole number and is refused as a count is, so
    // that a negative one is the `ValueError` of this argument and not the
    // `OverflowError` of pyo3.
    fn filter_by_ld(
        &self,
        max_allowed_r2: &Bound<'_, PyAny>,
        max_dist: &Bound<'_, PyAny>,
    ) -> Result<(), PyPopneiError> {
        let max_allowed_r2 = threshold_of(MAX_ALLOWED_R2, max_allowed_r2)?;
        let max_dist = distance_of(MAX_DIST, 1, max_dist)?;
        // The filter of this criterion is the core's `LdFilter` and not its
        // `VarFilter`: whether a variant is kept turns on the variants kept
        // behind it, so `VarFilter::new` refuses this criterion. The filter
        // built here is dropped and every pass builds its own, so what an
        // argument has to keep is written in the core alone.
        LdFilter::new(max_allowed_r2, max_dist)
            .map_err(|error| under_the_argument(error, MAX_ALLOWED_R2))?;
        let criterion = VarFilteringCriterion::MaxLdR2 {
            max_allowed_r2,
            max_dist,
        };
        let step = Step {
            pass_step: PassStep::VarFilter(criterion),
            args: vec![
                (MAX_ALLOWED_R2, Argument::Threshold(max_allowed_r2)),
                (MAX_DIST, Argument::Distance(max_dist)),
            ],
        };
        self.add(step)
    }

    // The genotypes of `individuals` kept at every variant, in the order
    // they are named here, and those of no other individual. The names are
    // resolved against the individuals of the source, so a name that is not
    // one of them, a name that is there twice and no name at all are
    // refused at this call and not when a pass runs.
    fn filter_individuals(&self, individuals: Vec<String>) -> Result<(), PyPopneiError> {
        // What the names give is dropped: every pass resolves them again
        // when it builds its chain, so the rule lives in the core alone.
        // They are refused before the list is looked at, since a name that
        // is of no individual is wrong whatever the list holds.
        resolve_individuals(&individuals, &self.of_the_source)?;
        let step = Step {
            pass_step: PassStep::KeepIndividuals(individuals.clone()),
            args: vec![(INDIVIDUALS, Argument::Individuals(individuals))],
        };
        self.add(step)
    }

    // The variants inside the regions of the BED file at `bed_path`, or,
    // with `exclude`, those outside all of them. The file is read at this
    // call, with the interpreter released, so an error of it comes here and
    // the regions are those it held now; every pass shares them. The file
    // is read before the list is looked at, since a file that cannot be
    // read is wrong whatever the list holds, and an error of it names the
    // file, where a second filter of a kind names none.
    fn filter_by_regions(
        &self,
        py: Python<'_>,
        bed_path: PathBuf,
        exclude: bool,
    ) -> Result<(), PyPopneiError> {
        let regions = py
            .detach(|| -> Result<Regions, popnei::Error> {
                let file =
                    File::open(&bed_path).map_err(|source| popnei::Error::FileNotOpened {
                        path: bed_path.clone(),
                        source,
                    })?;
                Regions::from_bed(BufReader::new(file))
            })
            .map_err(|error| PyPopneiError::of_the_file(error, &bed_path))?;
        let num_regions = regions.num_regions();
        let step = Step {
            pass_step: PassStep::Regions(RegionSelection {
                regions: Arc::new(regions),
                exclude,
            }),
            args: vec![
                (
                    BED_PATH,
                    Argument::Path(bed_path.to_string_lossy().into_owned()),
                ),
                (NUM_REGIONS, Argument::Count(num_regions)),
            ],
        };
        self.add(step)
    }

    // The first `num_vars` variants that the steps before this one keep,
    // after which every pass ends without reading the rest of the source.
    // The number is taken as the object it is and converted here, as a
    // count of variants, so that a negative one is the `ValueError` that
    // names the argument and not the `OverflowError` of pyo3, and `True` or
    // `1.5` the `TypeError` that names it.
    fn filter_first_n(&self, num_vars: &Bound<'_, PyAny>) -> Result<(), PyPopneiError> {
        let count = count_of_at_least(NUM_VARS, 1, num_vars)?;
        let num_vars = u64::try_from(count).map_err(|_| PyPopneiError::Count {
            name: NUM_VARS,
            smallest: 1,
            value: count.to_string(),
        })?;
        // The core's step refuses a `num_vars` of 0 here, where the user
        // wrote it, and not at the next pass.
        let step = Step {
            pass_step: first_n_step(num_vars)?,
            args: vec![(NUM_VARS, Argument::Count(count))],
        };
        self.add(step)
    }

    // Each variant kept with the probability `keep_rate`, drawn from a
    // generator that starts at `seed` in every pass. Both are taken as the
    // objects they are and converted here, so that a `True` or a `1.5` for
    // the seed is the `TypeError` that names it and a negative seed the
    // `ValueError` that names it, and not the `OverflowError` of pyo3.
    fn filter_randomly(
        &self,
        keep_rate: &Bound<'_, PyAny>,
        seed: &Bound<'_, PyAny>,
    ) -> Result<(), PyPopneiError> {
        let keep_rate = keep_rate_of(keep_rate)?;
        let seed = seed_of(seed)?;
        // The core's filter refuses a keep rate that is NaN, below 0 or
        // above 1, with a message that names `keep_rate`. It is dropped and
        // every pass builds its own, so the rule is written in the core
        // alone.
        RandomFilter::new(keep_rate, seed)?;
        let step = Step {
            pass_step: PassStep::Random { keep_rate, seed },
            args: vec![
                (KEEP_RATE, Argument::KeepRate(keep_rate)),
                (SEED, Argument::Seed(seed)),
            ],
        };
        self.add(step)
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
        &self,
        argument: &'static str,
        criterion: VarFilteringCriterion,
    ) -> Result<(), PyPopneiError> {
        VarFilter::new(criterion).map_err(|error| under_the_argument(error, argument))?;
        let step = Step {
            pass_step: PassStep::VarFilter(criterion),
            args: vec![(argument, Argument::Threshold(criterion.threshold()))],
        };
        self.add(step)
    }

    /// `step` added at the end of the list, where the next pass takes it.
    ///
    /// # Errors
    ///
    /// When the list holds a step of the kind of `step` already, and when
    /// `step` takes variants out and the list holds a filter of the first
    /// n; the core's `refuse_a_step` says both, in the order in which its
    /// `chain_of` asks them and in which TypeScript gets them. After
    /// either, the list is as it was. And when a panic left the lock
    /// broken, which is a defect of this crate.
    fn add(&self, step: Step) -> Result<(), PyPopneiError> {
        let mut steps = self.locked()?;
        let set = pass_steps_of(&steps);
        refuse_a_step(&set, &step.pass_step)?;
        steps.push(step);
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

/// `error`, and a threshold the core refused under `argument`, the name a
/// Python user wrote it in.
///
/// The core names the filter by its kind, `maf`, and the message a user
/// reads names `max_allowed_maf`, which only this crate knows: what they
/// have to look at is the call they wrote.
fn under_the_argument(error: popnei::Error, argument: &'static str) -> PyPopneiError {
    if let popnei::Error::VarFilterThresholdOutOfRange { threshold, .. } = error {
        return PyPopneiError::Threshold {
            name: argument,
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
/// binding crates call: what this one does is hand it the steps of the core
/// out of the steps of this crate, which carry the names of their arguments
/// in Python beside them. The readers belong to the pass this chain is built
/// for, so no count is shared with another pass, and the chain is asked for
/// the fields the consumer wants once it is built.
///
/// # Errors
///
/// When a threshold or a keep rate of the steps is not a number from 0 to
/// 1, when the steps hold two filters of one kind, when the filter of the
/// first n asks
/// for 0 variants, and when a filter that takes variants out comes after
/// it. A user reaches none of them: the call that adds a filter refuses
/// each, and the steps are read from there.
pub(crate) fn chain_of(
    reader: Box<dyn BlockReader>,
    steps: &[Step],
) -> Result<Box<dyn BlockReader>, popnei::Error> {
    popnei::filters::chain_of(reader, &pass_steps_of(steps))
}

/// What the pass does at each step, in the order of the steps: what the core
/// is given, out of what the steps of this crate hold.
pub(crate) fn pass_steps_of(steps: &[Step]) -> Vec<PassStep> {
    steps.iter().map(|step| step.pass_step.clone()).collect()
}
