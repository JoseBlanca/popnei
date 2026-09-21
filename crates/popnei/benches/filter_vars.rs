//! What the filter by the rate of missing genotypes costs a whole pass over
//! a dataset, on a given number of threads.
//!
//! It times one thing: a pass with the genotypes alone asked for, from
//! opening the file to the last block, over a VCF or over a vars file. With
//! `--max-missing-rate` the pass goes through a `FilteredReader` that keeps
//! the variants whose missing genotypes, divided by the individuals of the
//! dataset, are at most that rate; without it the same pass runs with no
//! filter. The cost of the filter is the difference of the two medians, so
//! the two are run back to back on a machine that is doing nothing else:
//! the difference is a tenth of what each of them takes, and the load of
//! the machine moves a median by more than that.
//!
//! A path that ends in `.vars` is read as a vars file and anything else as
//! a VCF, plain or gzipped. The genotypes alone are what a calculation over
//! them asks for and what the filter reads; the VCF is read with the
//! default options, ploidy 2 and the variants that passed their filters
//! alone, and both readers give the size of block popnei chooses for the
//! individuals of the file.
//!
//! The threads are those of a rayon pool the benchmark builds and reads
//! inside, and not rayon's global pool, so that the same command can time
//! one thread and many. The VCF reader parses a batch of lines on whatever
//! pool it runs in, and the filter reads the rows of a block on it too; the
//! reader of a vars file runs on the thread that calls it, so over a vars
//! file the threads are the filter's alone.
//!
//! One pass that is not timed comes before the timed runs, with the same
//! filter and the same fields. It pays the page faults of the first touch
//! of the memory a pass works in, which a process pays once, and it reads
//! the file once, so that the timed runs read it from the page cache and
//! not from the disc.
//!
//! It is run with cargo, which passes what comes after the two dashes to
//! it:
//!
//! ```text
//! cargo bench --bench filter_vars -- <path> --threads 18 --runs 5 \
//!     --max-missing-rate 0.1
//! ```
//!
//! cargo runs it with `crates/popnei` as its working directory, so a
//! relative path is read from there and an absolute one is the plainer
//! thing to give.
//!
//! `--threads` is 1 and `--runs` is 5 when they are not given. It prints
//! the wall time of each run, with the variants the pass gave and what the
//! filter was given and kept, and then the best, the median and the worst
//! of the times; with an even number of runs the median is the middle of
//! the two middle times. The median is what the numbers of "Speed" of
//! `docs/specs/filters.md` are taken from.
//!
//! How the two files are made. The VCF of "Speed" of
//! `docs/specs/filters.md`, 100000 variants of 1000 individuals whose
//! genotypes are missing at a rate of 0.03, is written by
//! `make_big_vcf.py`, beside this file, which says what it writes; it takes
//! 3.1 s and leaves 403 MB. The vars file is what popnei writes of it with
//! the size of block it chooses for 1000 individuals:
//!
//! ```text
//! uv run --no-project --with numpy python \
//!     crates/popnei/benches/make_big_vcf.py /Users/jose/devel/popnei-bench/big.vcf
//! uv run maturin develop --release && uv run python -c "import popnei; \
//!     popnei.write_vars(popnei.open_vcf('/Users/jose/devel/popnei-bench/big.vcf'), \
//!     '/Users/jose/devel/popnei-bench/big.vars')"
//! ```
//!
//! `maturin develop` is there because an `uv sync` takes the module out of
//! the environment, as `pyproject.toml` says, and `import popnei` then
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
use popnei::filters::{FilteredReader, VarFilter, VarFilteringCriterion};
use popnei::io::vars::VarsReader;
use popnei::io::vcf::{VcfOptions, VcfReader};
use popnei::variant::Needs;

/// How many times the pass is timed when the command line does not say.
const DEFAULT_RUNS: usize = 5;

/// How many threads the pool has when the command line does not say.
const DEFAULT_THREADS: usize = 1;

/// The end of the path of a vars file. A path that does not end in it is
/// read as a VCF, plain or gzipped.
const A_VARS_FILE_ENDS_IN: &str = ".vars";

/// What the command line asked for. `max_missing_rate` is `None` when the
/// command line did not name it, and the pass then runs with no filter.
struct Arguments {
    path: PathBuf,
    threads: usize,
    runs: usize,
    max_missing_rate: Option<f64>,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
filter_vars <path to a VCF or a vars file> [--threads n] [--runs n]
            [--max-missing-rate r]

It times whole passes over that file with the genotypes alone asked for:
opening the file, its header, every block to the end of it. With
--max-missing-rate the pass keeps the variants whose missing genotypes,
divided by the individuals of the dataset, are at most r; without it the
same pass runs with no filter, and the difference of the two medians is
what the filter costs.

  --threads n              how many threads the pool it reads in has, 1 by default
  --runs n                 how many times it reads the file, 5 by default
  --max-missing-rate r     the largest rate of missing genotypes that keeps
                           a variant, a number from 0 to 1; with none, no filter
  --help                   this

A path that ends in `.vars` is read as a vars file and anything else as a
VCF. One pass that is not timed comes first, so that the timed runs pay
neither the page faults of the first touch of the memory a pass works in
nor a read of the disc. It prints the wall time of each run, with the
variants the pass gave and what the filter was given and kept, and then
the best, the median and the worst of the times. The median is what
`docs/specs/filters.md` states, and the best and the worst say how much
the machine was doing something else.";

/// What to run, or the message that says what the command line should have
/// been.
fn arguments() -> Result<Arguments, String> {
    let mut path: Option<PathBuf> = None;
    let mut threads = DEFAULT_THREADS;
    let mut runs = DEFAULT_RUNS;
    let mut max_missing_rate: Option<f64> = None;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--threads" | "--runs" => {
                let name = arg.clone();
                let number = args
                    .next()
                    .ok_or_else(|| format!("{name} takes a number and none came after it"))?
                    .parse::<usize>()
                    .map_err(|_| format!("{name} takes a number"))?;
                if name == "--threads" {
                    threads = number;
                } else {
                    runs = number;
                }
            }
            "--max-missing-rate" => {
                max_missing_rate = Some(
                    args.next()
                        .ok_or_else(|| {
                            "--max-missing-rate takes a number and none came after it".to_owned()
                        })?
                        .parse::<f64>()
                        .map_err(|_| "--max-missing-rate takes a number".to_owned())?,
                );
            }
            // `cargo bench` adds this to the command line of every bench,
            // to tell a harness that has tests too to run its benchmarks.
            // This one has only this benchmark and takes it as nothing.
            "--bench" => {}
            "--help" | "-h" => return Err(USAGE.to_owned()),
            // An argument that is not one of those and looks like one is
            // refused instead of being taken for the path of the file, as
            // `read_vcf.rs` refuses it: a `--max-missing-rate=0.1` read as
            // a path and dropped leaves a run that timed a pass with no
            // filter and calls it the pass with one.
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
    if threads == 0 || runs == 0 {
        return Err("--threads and --runs are 1 or more".to_owned());
    }
    Ok(Arguments {
        path,
        threads,
        runs,
        max_missing_rate,
    })
}

/// One run: how long the pass took and the line that says what it gave, the
/// variants and, when there was a filter, what it was given and kept.
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

/// One whole pass over the file at `path` with the genotypes alone asked
/// for, and the filter over it when `max_missing_rate` is given, timed from
/// the building of the reader to the last block.
///
/// The counts of the filter are read after the clock stops: they are two
/// numbers of a chain that the pass has already built, and the line they go
/// into is printed and not timed.
fn one_pass(path: &Path, max_missing_rate: Option<f64>) -> Result<Run, popnei::Error> {
    let started = Instant::now();
    let source = reader_of(path)?;
    let mut reader: Box<dyn BlockReader> = match max_missing_rate {
        Some(rate) => {
            let filter = VarFilter::new(VarFilteringCriterion::MaxMissingRate(rate))?;
            Box::new(FilteredReader::new(source, filter)?)
        }
        None => source,
    };
    reader.set_needs(Needs::GTS);
    let mut variants: u64 = 0;
    while let Some(block) = reader.next_block()? {
        // A file of more variants than a u64 counts cannot be written.
        variants = variants.saturating_add(u64::try_from(block.num_vars).unwrap_or(u64::MAX));
    }
    let took = started.elapsed();
    let mut did = format!("{variants} variants");
    for (kind, counts) in reader.filtering_stats() {
        did.push_str(&format!(
            ", the {kind} filter was given {given} and kept {kept}",
            given = counts.vars_processed,
            kept = counts.vars_kept,
        ));
    }
    Ok(Run { took, did })
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

/// The benchmark builds a pool of threads and reads a file of the disc,
/// and wasm has neither; rayon is not a dependency of the wasm targets
/// either. This is what `cargo check --target wasm32-unknown-unknown
/// --all-targets` compiles of it, so that the command which checks that
/// nothing of the crate has left wasm behind can check the benchmarks too.
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
    let pool = match rayon::ThreadPoolBuilder::new()
        .num_threads(arguments.threads)
        .build()
    {
        Ok(pool) => pool,
        Err(error) => {
            eprintln!("the pool of {} threads: {error}", arguments.threads);
            return ExitCode::FAILURE;
        }
    };
    println!(
        "{path}, {threads} threads, {runs} runs, {filter}",
        path = arguments.path.display(),
        threads = arguments.threads,
        runs = arguments.runs,
        filter = match arguments.max_missing_rate {
            Some(rate) => format!("the missing data filter at {rate}"),
            None => "no filter".to_owned(),
        },
    );
    let mut times = Vec::with_capacity(arguments.runs);
    // The pass that is not timed, and then the timed ones. Both go through
    // `install`, so that the one that warms the memory runs on the same
    // pool as the ones that are timed.
    for run in 0..=arguments.runs {
        let done = match pool.install(|| one_pass(&arguments.path, arguments.max_missing_rate)) {
            Ok(done) => done,
            Err(error) => {
                eprintln!("{path}: {error}", path = arguments.path.display());
                return ExitCode::FAILURE;
            }
        };
        if run == 0 {
            println!(
                "the first pass, which is not timed: {did}, in {took}",
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
