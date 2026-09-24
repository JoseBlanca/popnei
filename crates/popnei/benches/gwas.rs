//! How long an association study of a file takes, for each of the four
//! models and each of the tests they have.
//!
//! It times one thing: `calc_gwas` over a file, from building the reader to
//! the result, which is what a user waits for. The null model is fitted
//! inside the clock, as it is for a user, and the kinship of the two mixed
//! models is computed before the clock starts, since a user brings one they
//! already have.
//!
//! Nothing had measured this module when the performance review of
//! `docs/reports/perf-gwas-2026-09-24.md` was written, and this is what it
//! measured it with. The program to compare the two models with no kinship
//! against is plink2's `--glm`, which fits the same linear model and the
//! same logistic regression, and the one for the two mixed models is
//! GMMAT's `glmm.score`; both are in `tests/reference/gwas/make_reference.py`
//! at the sizes the tests use, and the command that compares them at these
//! sizes is in the report.
//!
//! The trait, the covariates and the kinship are made here and not read
//! from a file, so that one command times a study of any file: the trait
//! is drawn from a generator of its own with a seed the command line
//! gives, so two runs of the same seed are the same study, and the
//! covariates are drawn beside it. Neither is built from the genotypes,
//! so no variant is associated with the trait beyond what chance gives,
//! which is what a study of a panel mostly holds. What the fits cost does
//! not turn on that: the null model is fitted once from the trait and the
//! design alone, and every variant then takes the same path through the
//! test, but for the logistic Wald test, where the rounds a variant's fit
//! takes depend on the trait and are not counted here: what that leaves
//! open is on the measurement plan of the report.
//!
//! ```text
//! VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1 \
//!     cargo bench --bench gwas -- <path> --model lm --runs 5
//! cargo bench --bench gwas -- <path> --model glm --test wald --runs 5
//! cargo bench --bench gwas -- <path> --model lmm --test score --runs 3
//! cargo bench --bench gwas -- <path> --model glmm --grammar-gamma --runs 3
//! ```
//!
//! cargo runs it with `crates/popnei` as its working directory, so a
//! relative path is read from there and an absolute one is the plainer
//! thing to give. A variable set after the process starts would not reach
//! Accelerate, so both go before the command and not into it.
//!
//! One run that is not timed comes before the timed ones. It pays the page
//! faults of the first touch of the memory a pass works in, which a
//! process pays once, and it reads the file once, so that the timed runs
//! read it from the page cache and not from the disc.
//!
//! `--runs` is 5 when it is not given, `--covariates` 2, `--seed` 1 and
//! `--model` `lm`. It prints the wall time of each run, with the variants
//! the pass gave, the ones that were answered, the individuals, and the
//! smallest p-value, which says the tests were made and not left a column
//! of zeros; and then the best, the median and the worst of the times. The
//! best is what the report of this measurement states, since every other
//! process on the machine can only make a run longer.
//!
//! How the panels are made is on `kinship.rs`, beside this file:
//! `make_big_vcf.py` writes 100000 variants of 1000 individuals and the
//! third argument is the rate at which a genotype is missing.

#![cfg_attr(
    target_family = "wasm",
    allow(
        dead_code,
        unused_imports,
        reason = "the benchmark is native: in wasm only its empty main is compiled, and \
                  what the timing is made of is left unused"
    )
)]

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{Duration, Instant};

use popnei::block::BlockReader;
use popnei::gwas::{GwasInput, TestType, TraitType, calc_gwas};
use popnei::io::vars::VarsReader;
use popnei::io::vcf::{VcfOptions, VcfReader};
use popnei::kinship::calc_kinship;

/// How many times the study is timed when the command line does not say.
const DEFAULT_RUNS: usize = 5;

/// How many covariates the design holds beside the intercept when the
/// command line does not say. It is what the reference script of
/// `tests/reference/gwas/` gives plink2 and GMMAT.
const DEFAULT_COVARIATES: usize = 2;

/// The seed of the trait and of the covariates when the command line does
/// not say.
const DEFAULT_SEED: u64 = 1;

/// The end of the path of a vars file. A path that does not end in it is
/// read as a VCF, plain or gzipped.
const A_VARS_FILE_ENDS_IN: &str = ".vars";

/// Which of the four models the study fits, which is the trait and the
/// kinship together, as `docs/specs/gwas.md` has it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Model {
    /// A continuous trait and no kinship.
    Lm,
    /// A continuous trait and a kinship.
    Lmm,
    /// A binomial trait and no kinship.
    Glm,
    /// A binomial trait and a kinship.
    Glmm,
}

impl Model {
    /// The model a user named, and `None` for a name that is of none.
    fn of_name(name: &str) -> Option<Model> {
        match name {
            "lm" => Some(Model::Lm),
            "lmm" => Some(Model::Lmm),
            "glm" => Some(Model::Glm),
            "glmm" => Some(Model::Glmm),
            _ => None,
        }
    }

    /// What the trait of this model is.
    fn trait_type(self) -> TraitType {
        match self {
            Model::Lm | Model::Lmm => TraitType::Continuous,
            Model::Glm | Model::Glmm => TraitType::Binomial,
        }
    }

    /// Whether this model takes a kinship, which is what the benchmark
    /// computes one for.
    fn takes_a_kinship(self) -> bool {
        matches!(self, Model::Lmm | Model::Glmm)
    }
}

/// What the command line asked for.
struct Arguments {
    path: PathBuf,
    model: Model,
    test: Option<TestType>,
    covariates: usize,
    seed: u64,
    grammar_gamma: bool,
    runs: usize,
    write_answers: Option<PathBuf>,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
gwas <path to a VCF or a vars file> [--model lm|lmm|glm|glmm] [--test wald|score]
     [--covariates n] [--seed n] [--grammar-gamma] [--runs n]
     [--write-answers path]

It times the association study of the variants of that file against a trait
this benchmark draws, from building the reader to the result: the null
model, the dosages of every block and the test of every variant.

  --model m      lm, a continuous trait and no kinship; lmm, one with a
                 kinship; glm, a binomial trait and no kinship; glmm, one
                 with a kinship. lm by default
  --test t       wald or score, and the default of the model when it is not
                 given
  --covariates n how many covariates the design holds beside the intercept,
                 2 by default
  --seed n       the seed of the trait and of the covariates, 1 by default
  --grammar-gamma  the GRAMMAR-Gamma approximation, which the two mixed
                 models have. It opens a second pass over the file
  --runs n       how many times it makes the study, 5 by default
  --write-answers path
                 writes `beta`, `se` and `p_value` of every variant to that
                 file, one variant a line and each value with all of its
                 digits, from the run that is not timed. It is what says
                 that a change which was made for speed alone moved no
                 number: the file of two runs is compared byte for byte
  --help         this

A path that ends in `.vars` is read as a vars file and anything else as a
VCF. The kinship of a mixed model is computed from the same file before the
clock starts, since a user brings one they already have. One run that is
not timed comes first, so that the timed runs pay neither the page faults
of the first touch of the memory a pass works in nor a read of the disc.

It builds no pool of threads: the dosages of a block are read on rayon's
global pool and the products run on the BLAS of the system, so one thread
is asked for with VECLIB_MAXIMUM_THREADS=1 and RAYON_NUM_THREADS=1 in the
environment of the command, which Accelerate reads when the process
starts.";

/// What to run, or the message that says what the command line should have
/// been.
fn arguments() -> Result<Arguments, String> {
    let mut path: Option<PathBuf> = None;
    let mut model = Model::Lm;
    let mut test = None;
    let mut covariates = DEFAULT_COVARIATES;
    let mut seed = DEFAULT_SEED;
    let mut grammar_gamma = false;
    let mut runs = DEFAULT_RUNS;
    let mut write_answers: Option<PathBuf> = None;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--model" => {
                let name = args
                    .next()
                    .ok_or_else(|| "--model takes a name and none came after it".to_owned())?;
                model = Model::of_name(&name)
                    .ok_or_else(|| format!("`{name}` is not one of lm, lmm, glm and glmm"))?;
            }
            "--test" => {
                let name = args
                    .next()
                    .ok_or_else(|| "--test takes a name and none came after it".to_owned())?;
                test = Some(
                    TestType::of_name(&name)
                        .map_err(|_| format!("`{name}` is not one of wald and score"))?,
                );
            }
            "--covariates" | "--runs" => {
                let name = arg.clone();
                let number = args
                    .next()
                    .ok_or_else(|| format!("{name} takes a number and none came after it"))?
                    .parse::<usize>()
                    .map_err(|_| format!("{name} takes a number of 0 or more"))?;
                if name == "--runs" {
                    runs = number;
                } else {
                    covariates = number;
                }
            }
            "--seed" => {
                seed = args
                    .next()
                    .ok_or_else(|| "--seed takes a number and none came after it".to_owned())?
                    .parse::<u64>()
                    .map_err(|_| "--seed takes a number of 0 or more".to_owned())?;
            }
            "--write-answers" => {
                write_answers = Some(PathBuf::from(args.next().ok_or_else(|| {
                    "--write-answers takes a path and none came after it".to_owned()
                })?));
            }
            "--grammar-gamma" => grammar_gamma = true,
            // `cargo bench` adds this to the command line of every bench,
            // to tell a harness that has tests too to run its benchmarks.
            // This one has only this benchmark and takes it as nothing.
            "--bench" => {}
            "--help" | "-h" => return Err(USAGE.to_owned()),
            // An argument that looks like one and is not one is refused
            // instead of being taken for the path of the file, as the
            // benchmarks beside this one refuse it: a `--runs=3` read as a
            // path and dropped leaves a run that timed nothing.
            other if other.starts_with('-') => {
                return Err(format!(
                    "`{other}` is not an argument of this benchmark\n\n{USAGE}"
                ));
            }
            other => path = Some(PathBuf::from(other)),
        }
    }
    let Some(path) = path else {
        return Err(format!("no file was given\n\n{USAGE}"));
    };
    if runs == 0 {
        return Err("--runs is 1 or more".to_owned());
    }
    if grammar_gamma && !model.takes_a_kinship() {
        return Err(
            "--grammar-gamma is of the two models that take a kinship, lmm and glmm".to_owned(),
        );
    }
    Ok(Arguments {
        path,
        model,
        test,
        covariates,
        seed,
        grammar_gamma,
        runs,
        write_answers,
    })
}

/// A generator of numbers that look drawn at random, so that the trait and
/// the covariates of a run are the same for a seed on every machine and
/// the benchmark depends on no crate that the core does not.
///
/// It is splitmix64: one state, one multiplication and three shifts per
/// number. What it is for is a phenotype that varies and that no variant
/// of the file explains, and not statistical quality.
struct Numbers {
    state: u64,
}

impl Numbers {
    /// The generator of a seed.
    fn of_the_seed(seed: u64) -> Numbers {
        Numbers { state: seed }
    }

    /// The next of its numbers, between 0 and 1.
    fn next(&mut self) -> f64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut value = self.state;
        value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        value ^= value >> 31;
        // The top 53 bits are the ones a f64 holds without rounding, and
        // dividing them by 2^53 gives a number from 0 to just below 1.
        ((value >> 11) as f64) / 9_007_199_254_740_992.0
    }

    /// The next of its numbers, drawn from a bell curve of mean 0 and
    /// spread 1: the sum of twelve of them less six, which is the plainest
    /// way to a bell and is good to the third decimal of its tails.
    fn next_of_the_bell(&mut self) -> f64 {
        let mut sum = 0.0_f64;
        for _ in 0..12 {
            sum += self.next();
        }
        sum - 6.0
    }
}

/// The trait of each individual and the design the model is fitted on.
struct TheStudy {
    /// One value per individual: a measurement, or 0.0 or 1.0.
    phenotype: Vec<f64>,
    /// `num_individuals` x `num_coefs`, row after row, the intercept first.
    design: Vec<f64>,
    /// How many columns that design has.
    num_coefs: usize,
    /// The positions of the tested individuals among the ones the reader
    /// gives, which here is all of them in their order.
    individuals: Vec<usize>,
}

/// The trait of `num_individuals` individuals and a design of `covariates`
/// covariates beside the intercept, drawn from `seed`.
///
/// A continuous trait is drawn from a bell curve and a binomial one is 0 or
/// 1 with the same chance, so about half of the individuals have the
/// condition, which is where a logistic fit settles in the fewest rounds
/// and so the friendliest case for the Wald test this benchmark times.
/// The covariates are drawn from the same bell curve, and no column is
/// built from the genotypes.
fn the_study_of(
    num_individuals: usize,
    covariates: usize,
    trait_type: TraitType,
    seed: u64,
) -> TheStudy {
    let mut numbers = Numbers::of_the_seed(seed);
    let mut phenotype = Vec::with_capacity(num_individuals);
    for _ in 0..num_individuals {
        phenotype.push(match trait_type {
            TraitType::Continuous => numbers.next_of_the_bell(),
            TraitType::Binomial => f64::from(u8::from(numbers.next() < 0.5)),
        });
    }
    let num_coefs = covariates.saturating_add(1);
    let mut design = Vec::with_capacity(num_individuals.saturating_mul(num_coefs));
    #[expect(
        clippy::same_item_push,
        reason = "the intercept of every individual is the same 1.0, and what the lint \
                  reads as one value pushed in a loop is the first column of a row whose \
                  other values are drawn"
    )]
    for _ in 0..num_individuals {
        design.push(1.0);
        for _ in 0..covariates {
            design.push(numbers.next_of_the_bell());
        }
    }
    TheStudy {
        phenotype,
        design,
        num_coefs,
        individuals: (0..num_individuals).collect(),
    }
}

/// One run: how long the study took and the line that says what it gave.
struct Run {
    took: Duration,
    did: String,
}

/// The three clocks of the phases of the pass, as a piece of the line of a
/// run: how long it was inside `next_block` of the reader, inside the
/// dosages of a block and inside the test of its variants. Taking them
/// zeroes them, so each run prints its own.
///
/// It is the cargo feature `bench-phases` of the core crate, and without
/// it there is nothing to print: the pass then calls no clock at all.
/// `cargo bench --features bench-phases --bench gwas` is what turns it on.
///
/// The second pass of `--grammar-gamma` is not in these numbers: its one
/// block is read before the loop of the pass, which is what the three
/// clocks are in.
#[cfg(feature = "bench-phases")]
fn the_phases_of_the_pass() -> String {
    let phases = popnei::gwas::phases::taken();
    format!(
        ", next_block {next_block:.4} s, dosages {dosages:.4} s, test {test:.4} s",
        next_block = phases.next_block.as_secs_f64(),
        dosages = phases.dosages.as_secs_f64(),
        test = phases.test.as_secs_f64(),
    )
}

/// Nothing, which is what the phases of the pass are when the cargo
/// feature `bench-phases` is off and the pass holds no clock.
#[cfg(not(feature = "bench-phases"))]
fn the_phases_of_the_pass() -> String {
    String::new()
}

/// The reader over the file at `path`: a vars file when the path ends in
/// `.vars`, and a VCF, plain or gzipped, when it does not.
fn reader_of(path: &Path) -> Result<Box<dyn BlockReader>, popnei::Error> {
    if path
        .to_string_lossy()
        .to_lowercase()
        .ends_with(A_VARS_FILE_ENDS_IN)
    {
        return Ok(Box::new(VarsReader::from_path(path)?));
    }
    Ok(Box::new(VcfReader::from_path(path, VcfOptions::default())?))
}

/// One whole study of the file at `path`, timed from the building of the
/// reader to the result.
///
/// The line it gives names how many variants were answered and the
/// smallest p-value, which is what says the tests were made: a result
/// whose columns were never filled would have neither.
///
/// `write_answers` is where the three columns of the result go, and `None`
/// is what a run that is only timed is given. The clock is stopped before
/// the file is written, so a run that writes one is timed as the others
/// are.
fn one_study(
    path: &Path,
    study: &TheStudy,
    trait_type: TraitType,
    test: Option<TestType>,
    kinship: Option<&[f64]>,
    grammar_gamma: bool,
    write_answers: Option<&Path>,
) -> Result<Run, popnei::Error> {
    let started = Instant::now();
    let mut reader = reader_of(path)?;
    let mut gamma_pass = match grammar_gamma {
        false => None,
        true => Some(reader_of(path)?),
    };
    let input = GwasInput {
        phenotype: &study.phenotype,
        trait_type,
        design: &study.design,
        num_coefs: study.num_coefs,
        kinship,
        test,
        use_grammar_gamma_approx: grammar_gamma,
        individuals: &study.individuals,
        transform_to_biallelic: false,
    };
    let gwas = calc_gwas(&mut reader, gamma_pass.as_mut(), &input)?;
    let took = started.elapsed();
    if let Some(answers) = write_answers {
        the_answers_of(&gwas, answers)?;
    }
    let answered = gwas.p_value.iter().filter(|p| p.is_finite()).count();
    let smallest = gwas
        .p_value
        .iter()
        .copied()
        .filter(|p| p.is_finite())
        .fold(f64::INFINITY, f64::min);
    let did = format!(
        "{num_vars} variants, {answered} of them answered, \
         {num_individuals} individuals, the smallest p-value is {smallest:.3e}{phases}",
        num_vars = gwas.num_vars,
        num_individuals = gwas.null_model.num_individuals,
        phases = the_phases_of_the_pass(),
    );
    Ok(Run { took, did })
}

/// The three columns of a study written to `path`, one variant a line and
/// the effect, its standard error and its p-value of that variant on it,
/// each with the seventeen decimals that tell two `f64` apart.
///
/// It is the oracle of a change made for speed alone: the file of the
/// changed code and the file of the code before it are the same bytes, or
/// the change moved a number and is not the same calculation. A variant
/// with no answer writes `NaN` three times, which is a difference the
/// comparison catches like any other.
///
/// # Errors
///
/// [`popnei::Error::Io`] when the file cannot be made or written.
fn the_answers_of(gwas: &popnei::gwas::Gwas, path: &Path) -> Result<(), popnei::Error> {
    let mut file = BufWriter::new(File::create(path)?);
    for ((beta, se), p_value) in gwas.beta.iter().zip(&gwas.se).zip(&gwas.p_value) {
        writeln!(file, "{beta:.17e}\t{se:.17e}\t{p_value:.17e}")?;
    }
    file.flush()?;
    Ok(())
}

/// The time at the place `part` of the times sorted from the shortest to
/// the longest.
fn sorted_time(times: &[Duration], part: usize) -> Duration {
    let mut sorted = times.to_vec();
    sorted.sort_unstable();
    sorted.get(part).copied().unwrap_or_default()
}

/// The median of the times: the middle one when they are an odd number,
/// and the middle of the two middle ones when they are an even number.
fn median_time(times: &[Duration]) -> Duration {
    let half = times.len() / 2;
    let upper = sorted_time(times, half);
    if times.len() % 2 == 1 {
        return upper;
    }
    let lower = sorted_time(times, half.saturating_sub(1));
    let between = upper.saturating_sub(lower);
    lower.saturating_add(between.checked_div(2).unwrap_or(between))
}

/// The seconds of a time, with three decimals.
fn seconds(time: Duration) -> String {
    format!("{:.3} s", time.as_secs_f64())
}

/// The benchmark reads a file of the disc and runs on the BLAS of the
/// system, and wasm has neither. This is what `cargo check --target
/// wasm32-unknown-unknown --all-targets` compiles of it.
#[cfg(target_family = "wasm")]
fn main() {}

#[cfg(not(target_family = "wasm"))]
fn main() -> ExitCode {
    let arguments = match arguments() {
        Ok(arguments) => arguments,
        Err(problem) => {
            eprintln!("{problem}");
            return ExitCode::FAILURE;
        }
    };
    let trait_type = arguments.model.trait_type();
    println!(
        "{path}, model {model:?}, test {test}, {covariates} covariates, seed {seed}, \
         grammar-gamma {grammar_gamma}, {runs} runs, VECLIB_MAXIMUM_THREADS {veclib}, \
         RAYON_NUM_THREADS {rayon}",
        path = arguments.path.display(),
        model = arguments.model,
        test = arguments.test.map_or("the model's default", TestType::name),
        covariates = arguments.covariates,
        seed = arguments.seed,
        grammar_gamma = arguments.grammar_gamma,
        runs = arguments.runs,
        veclib = std::env::var("VECLIB_MAXIMUM_THREADS").unwrap_or_else(|_| "unset".to_owned()),
        rayon = std::env::var("RAYON_NUM_THREADS").unwrap_or_else(|_| "unset".to_owned()),
    );
    // The individuals of the file, which the trait and the design are made
    // for, and the kinship of the two mixed models: both are outside the
    // clock, the first because a user has a phenotype already and the
    // second because a user brings the matrix they computed once.
    let num_individuals = match reader_of(&arguments.path) {
        Ok(reader) => reader.individuals().len(),
        Err(problem) => {
            eprintln!("the file could not be read: {problem}");
            return ExitCode::FAILURE;
        }
    };
    let study = the_study_of(
        num_individuals,
        arguments.covariates,
        trait_type,
        arguments.seed,
    );
    let kinship = match arguments.model.takes_a_kinship() {
        false => None,
        true => {
            let started = Instant::now();
            let of_the_panel = match reader_of(&arguments.path)
                .and_then(|mut reader| calc_kinship(&mut reader, None, false))
            {
                Ok(kinship) => kinship,
                Err(problem) => {
                    eprintln!("the kinship of the panel could not be computed: {problem}");
                    return ExitCode::FAILURE;
                }
            };
            println!(
                "the kinship of the panel, outside the clock: {took}",
                took = seconds(started.elapsed())
            );
            Some(of_the_panel.matrix)
        }
    };
    let kinship = kinship.as_deref();
    // The run that is not timed, which is also the one that writes the
    // answers when they were asked for: writing them inside a timed run
    // would time the file as well.
    if let Err(problem) = one_study(
        &arguments.path,
        &study,
        trait_type,
        arguments.test,
        kinship,
        arguments.grammar_gamma,
        arguments.write_answers.as_deref(),
    ) {
        eprintln!("the study could not be made: {problem}");
        return ExitCode::FAILURE;
    }
    let mut times = Vec::with_capacity(arguments.runs);
    for run in 1..=arguments.runs {
        match one_study(
            &arguments.path,
            &study,
            trait_type,
            arguments.test,
            kinship,
            arguments.grammar_gamma,
            None,
        ) {
            Ok(Run { took, did }) => {
                println!("run {run}: {took}, {did}", took = seconds(took));
                times.push(took);
            }
            Err(problem) => {
                eprintln!("the study could not be made: {problem}");
                return ExitCode::FAILURE;
            }
        }
    }
    println!(
        "best {best}, median {median}, worst {worst}",
        best = seconds(sorted_time(&times, 0)),
        median = seconds(median_time(&times)),
        worst = seconds(sorted_time(&times, times.len().saturating_sub(1))),
    );
    ExitCode::SUCCESS
}
