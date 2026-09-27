//! How long the VCF writer takes to write the variants of a VCF or of a
//! vars file as a VCF, plain or bgzipped, on a given number of threads.
//!
//! It times what `write_vcf` of the Python package does with a source and
//! no steps: the reader opened with the size of block that
//! `num_vars_per_block_of_write_vcf` gives, the chain of `chain_of` over
//! it, `write_vcf` of the core into a `BufWriter` of a file, the buffer
//! emptied into the file and the file synced to the disc. Each run builds
//! its own reader and writes the whole file, so what it times is one whole
//! pass: opening the source, reading every variant with the text of its
//! line, formatting the lines and, bgzipped, compressing them. The numbers
//! are those of "### The writer" of "Speed" of `docs/specs/io_vcf.md`.
//!
//! A path that ends in `.vars` is read as a vars file and anything else as
//! a VCF, plain or gzipped, with the options popnei chooses: the ploidy 2
//! and the variants that passed their filters alone.
//!
//! The threads are those of a rayon pool the benchmark builds and writes
//! inside, and not rayon's global pool, so that the same command can time
//! one thread and many: the VCF reader parses a batch of lines on it, and
//! the writer formats the rows of a block and compresses its members on
//! it.
//!
//! It is run with cargo, which passes what comes after the two dashes to
//! it:
//!
//! ```text
//! cargo bench --bench write_vcf -- <source> <output> --threads 18 --runs 3 --bgzip
//! ```
//!
//! `--threads` is 1 and `--runs` is 3 when they are not given, and without
//! `--bgzip` the file is plain text. `--no-sync` leaves out the sync of
//! the file to the disc, which the Python package does and bcftools does
//! not, so that the two can be told apart. The output is written again by
//! every run. It prints the wall time of each run and then the best, the
//! median and the worst of them. The first run of a source that was just
//! written reads it from the disc and the ones after it from the page
//! cache, so a timing that is reported reads the source once before.

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
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{Duration, Instant};

use popnei::block::BlockReader;
use popnei::filters::chain_of;
use popnei::io::vars::VarsReader;
use popnei::io::vcf::{
    VcfOptions, VcfReader, VcfWriteOptions, WriterSource, num_vars_per_block_of_write_vcf,
    write_vcf,
};

/// How many times the file is written when the command line does not say.
const DEFAULT_RUNS: usize = 3;

/// How many threads the pool has when the command line does not say: one,
/// the number the first of the targets of the spec is stated on.
const DEFAULT_THREADS: usize = 1;

/// The end of the path of a vars file, in any case.
const A_VARS_FILE_ENDS_IN: &str = ".vars";

/// What the command line asked for.
struct Arguments {
    source: PathBuf,
    output: PathBuf,
    threads: usize,
    runs: usize,
    bgzip: bool,
    sync: bool,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
write_vcf <source, a VCF or a .vars file> <output> [--threads n] [--runs n]
          [--bgzip] [--no-sync]

It times whole passes of the VCF writer over the source, as write_vcf of
the Python package makes them with no steps: the reader opened with the
size of block the writer asks for, every variant read with the text of
its line, the lines formatted and written into the output, bgzipped or
plain, and the output synced to the disc.

  --threads n   how many threads the pool it writes in has, 1 by default
  --runs n      how many times it writes the file, 3 by default
  --bgzip       members of bgzip instead of plain text
  --no-sync     no sync of the output to the disc after the last byte
  --help        this

It prints the wall time of each run and then the best, the median and
the worst of them.";

/// What to run, or the message that says what the command line should have
/// been.
fn arguments() -> Result<Arguments, String> {
    let mut paths: Vec<PathBuf> = Vec::new();
    let mut threads = DEFAULT_THREADS;
    let mut runs = DEFAULT_RUNS;
    let mut bgzip = false;
    let mut sync = true;
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
            "--bgzip" => bgzip = true,
            "--no-sync" => sync = false,
            // `cargo bench` adds this to the command line of every bench.
            "--bench" => {}
            "--help" | "-h" => return Err(USAGE.to_string()),
            other if other.starts_with('-') => {
                return Err(format!(
                    "`{other}` is not an argument of this benchmark\n\n{USAGE}"
                ));
            }
            other => paths.push(PathBuf::from(other)),
        }
    }
    let [source, output] = <[PathBuf; 2]>::try_from(paths)
        .map_err(|_| format!("a source and an output are given, no more\n\n{USAGE}"))?;
    if threads == 0 || runs == 0 {
        return Err("--threads and --runs are 1 or more".to_string());
    }
    Ok(Arguments {
        source,
        output,
        threads,
        runs,
        bgzip,
        sync,
    })
}

/// The reader over the source at `path`, with the size of block the
/// writer asks for: a vars file when the path ends in `.vars`, and a VCF,
/// plain or gzipped, when it does not.
fn reader_of(path: &Path) -> Result<Box<dyn BlockReader>, popnei::Error> {
    if path
        .to_string_lossy()
        .to_lowercase()
        .ends_with(A_VARS_FILE_ENDS_IN)
    {
        return Ok(Box::new(VarsReader::from_path(path)?));
    }
    let probe = VcfReader::from_path(path, VcfOptions::default())?;
    let options = VcfOptions {
        num_vars_per_block: num_vars_per_block_of_write_vcf(WriterSource::Vcf {
            num_individuals: probe.individuals().len(),
        }),
        ..VcfOptions::default()
    };
    Ok(Box::new(VcfReader::from_path(path, options)?))
}

/// One whole pass of the writer from `arguments.source` into
/// `arguments.output`, and how many variants it wrote.
fn write_the_whole_file(arguments: &Arguments) -> Result<u64, popnei::Error> {
    let mut chain = chain_of(reader_of(&arguments.source)?, &[])?;
    let file = File::create(&arguments.output).map_err(popnei::Error::Io)?;
    let options = VcfWriteOptions {
        bgzip: arguments.bgzip,
    };
    let (sink, num_vars) = write_vcf(&mut chain, BufWriter::new(file), options)?;
    let file = sink
        .into_inner()
        .map_err(|failure| popnei::Error::Io(failure.into_error()))?;
    if arguments.sync {
        file.sync_all().map_err(popnei::Error::Io)?;
    }
    Ok(num_vars)
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

/// The seconds of a time, with the three decimals of the other benches.
fn seconds(time: Duration) -> String {
    format!("{:.3} s", time.as_secs_f64())
}

/// The benchmark builds a pool of threads and writes a file of the disc,
/// and wasm has neither; this is what `cargo check --target
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
        "{source} into {output}, {threads} threads, {runs} runs, {how}{sync}",
        source = arguments.source.display(),
        output = arguments.output.display(),
        threads = arguments.threads,
        runs = arguments.runs,
        how = if arguments.bgzip { "bgzipped" } else { "plain" },
        sync = if arguments.sync {
            ", synced"
        } else {
            ", not synced"
        },
    );
    let mut times = Vec::with_capacity(arguments.runs);
    for run in 1..=arguments.runs {
        let started = Instant::now();
        let written = pool.install(|| write_the_whole_file(&arguments));
        let took = started.elapsed();
        match written {
            Ok(variants) => {
                times.push(took);
                println!("run {run}: {variants} variants in {}", seconds(took));
            }
            Err(error) => {
                eprintln!("{path}: {error}", path = arguments.source.display());
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
