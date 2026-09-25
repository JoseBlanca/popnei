//! How long the kinship of every pair of the individuals of a dataset
//! takes, from opening the file to the matrix.
//!
//! It times one thing: `calc_kinship` over a file, from building the
//! reader to the result, which is what a user waits for. With
//! `--num-pcs` above 0 the principal components of the matrix are taken
//! inside the clock too, so the two numbers say what the
//! eigendecomposition of an individuals x individuals matrix costs.
//!
//! The time to beat is in "Speed" of `docs/specs/kinship.md`: plink2's
//! 0.23 s for 100000 variants of 1000 individuals, on the panel with
//! every genotype called. The same panel with 3 in 100 genotypes missing
//! is timed beside it, since a block with a genotype missing does a
//! second product of the same shape for the per pair denominators.
//!
//! A path that ends in `.vars` is read as a vars file and anything else
//! as a VCF, plain or gzipped. The genotypes alone are asked for, which
//! is what the kinship reads, and neither reader is given a size of
//! block, so both give the size popnei chooses for the individuals of the
//! file; the pass puts `Reblock` before the reader anyway.
//!
//! The threads are not the benchmark's to choose and it builds no pool of
//! its own: the rows of a block are standardized on rayon's global pool,
//! which reads `RAYON_NUM_THREADS` when it is built, and the products and
//! the eigendecomposition run on the BLAS of the system, which on this
//! machine is Accelerate and reads `VECLIB_MAXIMUM_THREADS` when the
//! process starts. So one thread is asked for by setting both in the
//! environment of the command, and the threads the machine gives by
//! setting neither:
//!
//! ```text
//! VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1 \
//!     cargo bench --bench kinship -- <path> --runs 5
//! cargo bench --bench kinship -- <path> --runs 5 --num-pcs 10
//! ```
//!
//! cargo runs it with `crates/popnei` as its working directory, so a
//! relative path is read from there and an absolute one is the plainer
//! thing to give. A variable set after the process starts would not reach
//! Accelerate, so both go before the command and not into it.
//!
//! One run that is not timed comes before the timed ones. It pays the
//! page faults of the first touch of the memory a pass works in, which a
//! process pays once, and it reads the file once, so that the timed runs
//! read it from the page cache and not from the disc.
//!
//! `--num-pcs` is 0 and `--runs` is 5 when they are not given. It prints
//! the wall time of each run, with the variants the pass gave, the ones
//! that had variance and were used, the individuals, the mean of the
//! diagonal and the largest entry of the matrix, which say the matrix was
//! computed and not left a shape of zeros; and then the best, the median
//! and the worst of the times. The best is what the report of this
//! measurement states, since every other process on the machine can only
//! make a run longer.
//!
//! How the two panels are made. `make_big_vcf.py`, beside this file,
//! writes 100000 variants of 1000 individuals; its third argument is the
//! rate at which a genotype is missing, 0.03 when it is not given and 0
//! for the panel the target is stated on. The two hold the same
//! genotypes and differ in which of them are `./.`. The vars file of each
//! is what popnei writes of it with the size of block it chooses for 1000
//! individuals:
//!
//! ```text
//! uv run --no-project --with numpy python \
//!     crates/popnei/benches/make_big_vcf.py \
//!     /Users/jose/devel/popnei-bench/bigcalled.vcf 100000 0
//! uv run maturin develop --release && uv run python -c "import popnei; \
//!     popnei.write_vars(popnei.open_vcf('/Users/jose/devel/popnei-bench/bigcalled.vcf'), \
//!     '/Users/jose/devel/popnei-bench/bigcalled.vars')"
//! ```
//!
//! `maturin develop` is there because an `uv sync` takes the module out
//! of the environment, as `pyproject.toml` says, and `import popnei` then
//! fails with a `ModuleNotFoundError`.

#![cfg_attr(
    target_family = "wasm",
    allow(
        dead_code,
        unused_imports,
        reason = "the benchmark is native: in wasm only its empty main is compiled, and \
                  what the timing is made of is left unused"
    )
)]

use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{Duration, Instant};

use popnei::block::BlockReader;
use popnei::io::vars::VarsReader;
use popnei::io::vcf::{VcfOptions, VcfReader};
use popnei::kinship::{calc_kinship, principal_components};

/// How many times the kinship is timed when the command line does not
/// say.
const DEFAULT_RUNS: usize = 5;

/// How many principal components of the matrix are taken when the command
/// line does not say. With 0 there is no eigendecomposition.
const DEFAULT_NUM_PCS: usize = 0;

/// The end of the path of a vars file. A path that does not end in it is
/// read as a VCF, plain or gzipped.
const A_VARS_FILE_ENDS_IN: &str = ".vars";

/// What the command line asked for.
struct Arguments {
    path: PathBuf,
    runs: usize,
    num_pcs: usize,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
kinship <path to a VCF or a vars file> [--runs n] [--num-pcs n]

It times the kinship of every pair of the individuals of that file, from
building the reader to the matrix: the standardized blocks, the product of
each with itself, the per pair denominators when a genotype is missing,
and, when --num-pcs is above 0, the eigendecomposition of the matrix and
the projections that come off it.

  --runs n     how many times it takes the kinship, 5 by default
  --num-pcs n  how many principal components of the matrix are taken, 0 by
               default, and with 0 there is no eigendecomposition
  --help       this

A path that ends in `.vars` is read as a vars file and anything else as a
VCF. One run that is not timed comes first, so that the timed runs pay
neither the page faults of the first touch of the memory a pass works in
nor a read of the disc. It prints the wall time of each run with what the
pass gave, and then the best, the median and the worst of the times.

It builds no pool of threads: the rows of a block are standardized on
rayon's global pool and the products run on the BLAS of the system, so one
thread is asked for with VECLIB_MAXIMUM_THREADS=1 and RAYON_NUM_THREADS=1
in the environment of the command, which Accelerate reads when the process
starts.";

/// What to run, or the message that says what the command line should
/// have been.
fn arguments() -> Result<Arguments, String> {
    let mut path: Option<PathBuf> = None;
    let mut runs = DEFAULT_RUNS;
    let mut num_pcs = DEFAULT_NUM_PCS;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--runs" | "--num-pcs" => {
                let name = arg.clone();
                let number = args
                    .next()
                    .ok_or_else(|| format!("{name} takes a number and none came after it"))?
                    .parse::<usize>()
                    .map_err(|_| format!("{name} takes a number of 0 or more"))?;
                if name == "--runs" {
                    runs = number;
                } else {
                    num_pcs = number;
                }
            }
            // `cargo bench` adds this to the command line of every bench,
            // to tell a harness that has tests too to run its benchmarks.
            // This one has only this benchmark and takes it as nothing.
            "--bench" => {}
            "--help" | "-h" => return Err(USAGE.to_owned()),
            // An argument that is not one of those and looks like one is
            // refused instead of being taken for the path of the file, as
            // `pca_vars.rs` refuses it: a `--num-pcs=10` read as a path
            // and dropped leaves a run that took no components and calls
            // it the run that did.
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
    Ok(Arguments {
        path,
        runs,
        num_pcs,
    })
}

/// One run: how long the kinship took and the line that says what it
/// gave.
struct Run {
    took: Duration,
    did: String,
}

/// The reader over the file at `path`: a vars file when the path ends in
/// `.vars`, and a VCF, plain or gzipped, when it does not.
///
/// The VCF is read with the options popnei chooses, and neither reader is
/// asked for a size of block, so both give the size popnei picks for the
/// individuals of the file.
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

/// One whole kinship of the file at `path`, timed from the building of
/// the reader to the matrix, with its principal components inside the
/// clock when `num_pcs` is above 0.
///
/// The line it gives names the mean of the diagonal and the largest entry
/// of the matrix, which is what says the entries were computed: a matrix
/// whose numbers were never filled in would have neither.
fn one_kinship(path: &Path, num_pcs: usize) -> Result<Run, popnei::Error> {
    let started = Instant::now();
    let mut reader = reader_of(path)?;
    let kinship = calc_kinship(&mut reader, None, false)?;
    let pcs = match num_pcs {
        0 => None,
        _ => Some(principal_components(&kinship, num_pcs)?),
    };
    let took = started.elapsed();
    let num_individuals = kinship.num_individuals;
    let diagonal: f64 = (0..num_individuals)
        .filter_map(|one| {
            kinship
                .matrix
                .get(one.saturating_mul(num_individuals).saturating_add(one))
        })
        .sum();
    let mean_of_the_diagonal = match num_individuals {
        0 => f64::NAN,
        _ => diagonal / num_individuals as f64,
    };
    let largest = kinship
        .matrix
        .iter()
        .fold(0.0_f64, |largest, value| largest.max(value.abs()));
    let said_about_the_components = match pcs {
        None => "no components".to_owned(),
        Some(pcs) => format!("{num_comps} components", num_comps = pcs.num_comps),
    };
    let did = format!(
        "{num_vars_given} variants, {num_vars} of them with variance, \
         {num_individuals} individuals, the diagonal averages \
         {mean_of_the_diagonal:.6}, the largest entry is {largest:.6}, \
         {said_about_the_components}{phases}",
        num_vars_given = kinship.num_vars_given,
        num_vars = kinship.num_vars,
        phases = the_phases_of_the_pass(),
    );
    Ok(Run { took, did })
}

/// The two clocks of the phases of the pass, as a piece of the line of a
/// run: how long the pass was inside `next_block` of its reader and how
/// long it was working on the blocks the reader gave. Taking them zeroes
/// them, so each run prints its own.
///
/// It is the cargo feature `bench-phases` of the core crate, and without it
/// there is nothing to print: the pass then calls no clock at all. `cargo
/// bench --features bench-phases --bench kinship` is what turns it on. The
/// two say what the reader of `with_one_block_ahead` could save this pass,
/// which is the smaller of them, and `docs/reports/perf-read-ahead-2026-09-25.md`
/// has what they said.
#[cfg(feature = "bench-phases")]
fn the_phases_of_the_pass() -> String {
    let phases = popnei::phases::taken();
    format!(
        ", next_block {next_block:.4} s, work {work:.4} s",
        next_block = phases.next_block.as_secs_f64(),
        work = phases.work.as_secs_f64(),
    )
}

/// Nothing, which is what the phases of the pass are when the cargo feature
/// `bench-phases` is off and the pass holds no clock.
#[cfg(not(feature = "bench-phases"))]
fn the_phases_of_the_pass() -> String {
    String::new()
}

/// The time at the place `part` of the times sorted from the shortest to
/// the longest.
fn sorted_time(times: &[Duration], part: usize) -> Duration {
    let mut sorted = times.to_vec();
    sorted.sort_unstable();
    sorted.get(part).copied().unwrap_or_default()
}

/// The median of the times: the middle one when they are an odd number,
/// and the middle of the two middle ones when they are an even number,
/// which `sorted_time` at half of them is not.
fn median_time(times: &[Duration]) -> Duration {
    let half = times.len() / 2;
    let upper = sorted_time(times, half);
    if times.len() % 2 == 1 {
        return upper;
    }
    let lower = sorted_time(times, half.saturating_sub(1));
    // The two are a time each and `lower` is the shorter of them, so half
    // of what lies between them added to it is their middle, and neither
    // the subtraction nor the addition leaves what a `Duration` holds.
    let between = upper.saturating_sub(lower);
    lower.saturating_add(between.checked_div(2).unwrap_or(between))
}

/// The seconds of a time, with the three decimals that a pass of a few
/// tenths of a second is worth reporting to.
fn seconds(time: Duration) -> String {
    format!("{:.3} s", time.as_secs_f64())
}

/// The benchmark reads a file of the disc and runs on the BLAS of the
/// system, and wasm has neither. This is what `cargo check --target
/// wasm32-unknown-unknown --all-targets` compiles of it, so that the
/// command which checks that nothing of the crate has left wasm behind
/// can check the benchmarks too.
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
    println!(
        "{path}, {runs} runs, {num_pcs} principal components, \
         VECLIB_MAXIMUM_THREADS {veclib}, RAYON_NUM_THREADS {rayon}",
        path = arguments.path.display(),
        runs = arguments.runs,
        num_pcs = arguments.num_pcs,
        veclib = said_about_the_variable("VECLIB_MAXIMUM_THREADS"),
        rayon = said_about_the_variable("RAYON_NUM_THREADS"),
    );
    let mut times = Vec::with_capacity(arguments.runs);
    for run in 0..=arguments.runs {
        let done = match one_kinship(&arguments.path, arguments.num_pcs) {
            Ok(done) => done,
            Err(error) => {
                eprintln!("{path}: {error}", path = arguments.path.display());
                return ExitCode::FAILURE;
            }
        };
        if run == 0 {
            println!(
                "the first run, which is not timed: {did}, in {took}",
                did = done.did,
                took = seconds(done.took),
            );
            continue;
        }
        println!(
            "run {run}: {took}, {did}",
            took = seconds(done.took),
            did = done.did,
        );
        times.push(done.took);
    }
    let last = times.len().saturating_sub(1);
    println!(
        "best {best}, median {median}, worst {worst}",
        best = seconds(sorted_time(&times, 0)),
        median = seconds(median_time(&times)),
        worst = seconds(sorted_time(&times, last)),
    );
    ExitCode::SUCCESS
}

/// What the environment variable `name` holds, or that it is not set,
/// which is printed with every run: the threads of a timing are not the
/// benchmark's to choose and a number taken with the wrong ones is the
/// thing this measurement can get wrong without showing it.
#[cfg(not(target_family = "wasm"))]
fn said_about_the_variable(name: &str) -> String {
    match std::env::var(name) {
        Ok(value) => value,
        Err(_) => "not set".to_owned(),
    }
}
