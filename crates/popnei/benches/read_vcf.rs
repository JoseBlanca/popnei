//! How long the VCF reader takes to read a VCF, plain or gzipped, on a
//! given number of threads.
//!
//! It asks for the genotypes alone, which is what a calculation over the
//! genotypes asks for and what the numbers of "Speed" of
//! `docs/specs/io_vcf.md` were taken with, and it reads with the default
//! options, ploidy 2 and the variants that passed their filters alone. Each
//! run builds its own reader and reads the file to its end, so what it
//! times is one whole read: opening the file, the header, the lines, the
//! parse and the genotypes.
//!
//! The threads are those of a rayon pool the benchmark builds and reads
//! inside, and not rayon's global pool, so that the same command can time
//! one thread and many. The reader parses a batch of lines on whatever pool
//! it runs in.
//!
//! It is run with cargo, which passes what comes after the two dashes to
//! it:
//!
//! ```text
//! cargo bench --bench read_vcf -- <path> --threads 18 --runs 5
//! ```
//!
//! `--lines-per-batch` and `--bytes-per-batch` are the two bounds of a
//! batch, the number that "How it runs" of `docs/specs/io_vcf.md` says has
//! to be measured and the bound in bytes that holds the memory of a reader
//! of a file of many individuals. Without them the reader uses the
//! constants of `io::vcf`, which is what a timing that is reported for
//! popnei as it ships has to use.
//!
//! cargo runs it with `crates/popnei` as its working directory, so a
//! relative path is read from there and an absolute one is the plainer
//! thing to give.
//!
//! `--threads` is 1 and `--runs` is 5 when they are not given. It prints
//! the wall time of each run and then the best, the median and the worst of
//! them; with an even number of runs the median is the middle of the two
//! middle times. The first run of a file that was just written reads it from the
//! disk and the ones after it from the page cache, so a timing that is
//! reported leaves the first run out or reads the file once before.
//!
//! The file of "Speed" of `docs/specs/io_vcf.md`, 100000 variants of 1000
//! individuals, is written by `make_big_vcf.py`, beside this file, which
//! says what it is and how it is run. It takes 3.1 s and leaves 403 MB;
//! `bgzip -k` on it makes the gzipped one, 38 MB.

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
use popnei::io::vcf::{VcfOptions, VcfReader};
use popnei::variant::Needs;

/// How many times the file is read when the command line does not say.
const DEFAULT_RUNS: usize = 5;

/// How many threads the pool has when the command line does not say: one,
/// the number the target of the spec is stated on.
const DEFAULT_THREADS: usize = 1;

/// What the command line asked for. The two bounds of a batch are `None`
/// when the command line did not name them, and the reader then uses the
/// constants of `io::vcf`.
struct Arguments {
    path: PathBuf,
    threads: usize,
    runs: usize,
    lines_per_batch: Option<usize>,
    bytes_per_batch: Option<usize>,
}

/// What to run, or the message that says what the command line should have
/// been.
fn arguments() -> Result<Arguments, String> {
    let mut path: Option<PathBuf> = None;
    let mut threads = DEFAULT_THREADS;
    let mut runs = DEFAULT_RUNS;
    let mut lines_per_batch: Option<usize> = None;
    let mut bytes_per_batch: Option<usize> = None;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let mut number = |name: &str| {
            args.next()
                .ok_or_else(|| format!("{name} takes a number and none came after it"))?
                .parse::<usize>()
                .map_err(|_| format!("{name} takes a number"))
        };
        match arg.as_str() {
            "--threads" => threads = number("--threads")?,
            "--runs" => runs = number("--runs")?,
            "--lines-per-batch" => lines_per_batch = Some(number("--lines-per-batch")?),
            "--bytes-per-batch" => bytes_per_batch = Some(number("--bytes-per-batch")?),
            // `cargo bench` adds this to the command line of every bench,
            // to tell a harness that has tests too to run its benchmarks.
            // This one has only this benchmark and takes it as nothing.
            "--bench" => {}
            other => path = Some(PathBuf::from(other)),
        }
    }
    let Some(path) = path else {
        return Err(
            "no VCF was given: read_vcf <path> [--threads n] [--runs n] \
             [--lines-per-batch n] [--bytes-per-batch n]"
                .to_string(),
        );
    };
    if threads == 0 || runs == 0 {
        return Err("--threads and --runs are 1 or more".to_string());
    }
    if lines_per_batch == Some(0) || bytes_per_batch == Some(0) {
        return Err("--lines-per-batch and --bytes-per-batch are 1 or more".to_string());
    }
    Ok(Arguments {
        path,
        threads,
        runs,
        lines_per_batch,
        bytes_per_batch,
    })
}

/// It reads every block of the VCF at `path` with the genotypes asked for,
/// and gives how many variants there were.
///
/// The blocks are the size that popnei chooses for the individuals of the
/// file, which is what a user who asks for no size gets. The two bounds of
/// a batch are set only when the command line named them, through the two
/// hidden setters of the reader, which are there for this benchmark and for
/// the tests.
fn read_the_whole_file(path: &Path, arguments: &Arguments) -> Result<u64, popnei::Error> {
    let mut reader = VcfReader::from_path(path, VcfOptions::default())?;
    if let Some(lines) = arguments.lines_per_batch {
        reader.set_lines_per_batch(lines);
    }
    if let Some(bytes) = arguments.bytes_per_batch {
        reader.set_bytes_per_batch(bytes);
    }
    reader.set_needs(Needs::GTS);
    let mut variants: u64 = 0;
    while let Some(block) = reader.next_block()? {
        // A file of more variants than a u64 counts cannot be written.
        variants = variants.saturating_add(u64::try_from(block.num_vars).unwrap_or(u64::MAX));
    }
    Ok(variants)
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

/// The seconds of a time, with the three decimals that a read of a few
/// tenths of a second is worth reporting to.
fn seconds(time: Duration) -> String {
    format!("{:.3} s", time.as_secs_f64())
}

/// The benchmark builds a pool of threads and reads a file of the disk,
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
        "{path}, {threads} threads, {runs} runs, {lines} lines and {bytes} bytes in a batch",
        path = arguments.path.display(),
        threads = arguments.threads,
        runs = arguments.runs,
        lines = match arguments.lines_per_batch {
            Some(lines) => lines.to_string(),
            None => "the constant of io::vcf".to_string(),
        },
        bytes = match arguments.bytes_per_batch {
            Some(bytes) => bytes.to_string(),
            None => "the constant of io::vcf".to_string(),
        },
    );
    let mut times = Vec::with_capacity(arguments.runs);
    for run in 1..=arguments.runs {
        let started = Instant::now();
        let read = pool.install(|| read_the_whole_file(&arguments.path, &arguments));
        let took = started.elapsed();
        match read {
            Ok(variants) => {
                times.push(took);
                println!("run {run}: {variants} variants in {}", seconds(took));
            }
            Err(error) => {
                eprintln!("{path}: {error}", path = arguments.path.display());
                return ExitCode::FAILURE;
            }
        }
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
