//! How long the three statistics that a page shows when a file is opened
//! take in three passes, one for each, against the one pass of
//! `calc_variants_summary`, which gives the same three.
//!
//! It times two things over each file it is given, from building the
//! readers to the results, which is what a user waits for:
//!
//! - the three: `calc_per_var_distribs` with the six statistics, the
//!   defaults of `calcPerVarDistribs` in the TypeScript package, then
//!   `calc_per_individual_stats`, then `calc_var_density` with windows of
//!   100000 base pairs and no lengths of the chromosomes, each over a
//!   reader of its own, one after another;
//! - the one: `calc_variants_summary` asked for those same three with the
//!   same options, over one reader.
//!
//! The defaults of the distributions are those of the TypeScript package:
//! no populations, 40 linear bins from 0 to 1, 20 individuals a population
//! needs called at a variant, the ploidy of the variants as the exponent of
//! the expected heterozygosities and 0.95 as the threshold of the
//! polymorphism ratio. It is the measurement of deliverable 4 of work
//! package 1 of `docs/plans/stats-so-far.md`, and "The three statistics of a
//! file in one pass" of `docs/specs/js_sources.md` states what it gave.
//!
//! A path that ends in `.vars` is read as a vars file and anything else as
//! a VCF, plain or gzipped, with the options popnei chooses; neither reader
//! is given a size of block. Each pass asks its reader for what it reads,
//! as it does when a page calls it: the genotypes for the distributions and
//! the rates, the chromosome and the position alone for the density, and
//! both for the one pass.
//!
//! The threads are set here, `--threads` of them, 1 when it is not given,
//! by building rayon's global pool with that many, and not a pool of its
//! own that the passes are called inside: the blocks are read one block
//! ahead on a thread of their own, as every pass of these statistics reads
//! them natively, and the VCF reader parses the rows of a block on rayon
//! from that thread, which is in no pool and so uses the global one. A pool
//! of one thread installed around the passes left that parsing on the 18
//! threads of the machine, and the three over `big.vcf` took 0.53 s where
//! the read with the genotypes alone takes 0.575 s on one thread. The
//! reading thread itself is one thread more beside the pool, in both.
//!
//! The runs are interleaved, so that what the machine is doing meanwhile
//! falls on both alike: each run times the three and then the one over the
//! first file, then over the next, and so on. One run that is not timed
//! comes first; it pays the page faults of the first touch of the memory a
//! pass works in, which a process pays once, and it reads every file, so
//! that the timed runs read them from the page cache and not from the disc.
//!
//! It prints the wall time of the three and of the one for each run and
//! file, and then, for each file, the best, the median and the worst of
//! each and the saving of the one on the medians, in 100 parts of the time
//! of the three. A run fails, and says why, when the one pass does not give
//! the variants and the windows the three gave, which a reader or an option
//! built wrong would.
//!
//! It is run with cargo, which runs it with `crates/popnei` as its working
//! directory, so an absolute path is the plainer thing to give:
//!
//! ```text
//! cargo bench --bench variants_summary -- \
//!     /Users/jose/devel/popnei-bench/big.vcf \
//!     /Users/jose/devel/popnei-bench/big.vars --runs 7
//! ```
//!
//! `big.vcf` is the plain VCF of 403 MB, 100000 variants of 1000 diploid
//! individuals, that `make_big_vcf.py` writes, and `big.vars` the vars file
//! popnei writes of it, as `kinship.rs` says. On the owner's Apple M5 Pro,
//! 18 cores, macOS 27.0.1, with the core crate of b797829 of the branch
//! `spec/stats-so-far`, on 7 October 2026, at a load average of 15 from
//! other work on the machine, the medians of the 7 runs of that command
//! were:
//!
//! | | the three | the one | the one saves |
//! |---|---|---|---|
//! | `big.vcf` | 2.240 s | 1.247 s | 44 in 100 |
//! | `big.vars` | 0.400 s | 0.391 s | 2 in 100 |
//!
//! Two sets of 7 runs of the same command just before it gave savings of 43
//! and 44 in 100 on `big.vcf` and of 3 and 2 in 100 on `big.vars`.
//!
//! With `--features bench-phases` each pass prints two clocks, the time it
//! waited for the next block and the time it worked on the blocks, and they
//! say what each saving is made of. Over `big.vcf` a pass of the
//! distributions or of the rates waited about 0.91 s and worked 0.19 s: the
//! VCF reader parses the genotypes on rayon's global pool, whose one thread
//! the counting runs on too, so the read is not hidden behind the counting,
//! and the one pass saves one read of the genotypes and the pass of the
//! density, 0.07 s. Over `big.vars` every pass waited under 0.01 s for its
//! blocks and the distributions and the rates worked 0.19 s each, so the
//! read of that file is not what a pass waits for, and the one pass saves
//! little more than the pass of the density and the opening of two
//! readers. The one pass worked 0.38 s there, the two 0.19 s added up,
//! since the counting of each statistic is not shared.

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
use popnei::stats::{
    DEFAULT_HIST_RANGE, DEFAULT_MIN_NUM_INDIVIDUALS, DEFAULT_NUM_BINS, DEFAULT_POLY_THRESHOLD,
    ExpHet, HistBins, Maf, ObsHet, PerVarDistribsConfig, PerVarStat, Pops, VarDensity,
    VarDensityConfig, VariantsSummary, VariantsSummaryConfig, calc_per_individual_stats,
    calc_per_var_distribs, calc_var_density, calc_variants_summary, nothing_after_a_block,
};

/// How many times each of the two is timed over each file when the command
/// line does not say.
const DEFAULT_RUNS: usize = 7;

/// How many threads the counting runs on when the command line does not
/// say.
const DEFAULT_THREADS: usize = 1;

/// The width of a window of the density, in base pairs.
const WINDOW_SIZE: u64 = 100_000;

/// The end of the path of a vars file. A path that does not end in it is
/// read as a VCF, plain or gzipped.
const A_VARS_FILE_ENDS_IN: &str = ".vars";

/// What the command line asked for.
struct Arguments {
    paths: Vec<PathBuf>,
    runs: usize,
    threads: usize,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
variants_summary <path to a VCF or a vars file>... [--runs n] [--threads n]

Over each file, it times the three statistics a page shows when a file is
opened in three passes, the distributions of the six statistics of each
variant, the rates of each individual and the density in windows of 100000
base pairs, each over a reader of its own, against the one pass of
calc_variants_summary that gives the same three.

  --runs n     how many times each is timed over each file, 7 by default
  --threads n  the threads of rayon's global pool, 1 by default
  --help       this

A path that ends in `.vars` is read as a vars file and anything else as a
VCF. The runs are interleaved, the three and the one over each file in
turn, after one run that is not timed. It prints the wall time of each, and
then for each file the best, the median and the worst of each and what the
one saves on the medians.";

/// What to run, or the message that says what the command line should
/// have been.
fn arguments() -> Result<Arguments, String> {
    let mut paths = Vec::new();
    let mut runs = DEFAULT_RUNS;
    let mut threads = DEFAULT_THREADS;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--runs" | "--threads" => {
                let name = arg.clone();
                let number = args
                    .next()
                    .ok_or_else(|| format!("{name} takes a number and none came after it"))?
                    .parse::<usize>()
                    .map_err(|_| format!("{name} takes a number of 1 or more"))?;
                if number == 0 {
                    return Err(format!("{name} takes a number of 1 or more"));
                }
                if name == "--runs" {
                    runs = number;
                } else {
                    threads = number;
                }
            }
            // `cargo bench` adds this to the command line of every bench,
            // to tell a harness that has tests too to run its benchmarks.
            // This one has only this benchmark and takes it as nothing.
            "--bench" => {}
            "--help" | "-h" => return Err(USAGE.to_owned()),
            // An argument that looks like one and is not one of those is
            // refused instead of being taken for the path of a file, as
            // `kinship.rs` refuses it.
            other if other.starts_with('-') => {
                return Err(format!(
                    "`{other}` is not an argument of this benchmark\n\n{USAGE}"
                ));
            }
            other => paths.push(PathBuf::from(other)),
        }
    }
    if paths.is_empty() {
        return Err(format!("no file was given\n\n{USAGE}"));
    }
    Ok(Arguments {
        paths,
        runs,
        threads,
    })
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

/// The options of the distributions that `calcPerVarDistribs` of the
/// TypeScript package passes the core when it is given none, for variants
/// of `ploidy` of `num_individuals` individuals.
fn the_default_distribs(
    ploidy: usize,
    num_individuals: usize,
) -> Result<PerVarDistribsConfig, popnei::Error> {
    let (start, end) = DEFAULT_HIST_RANGE;
    Ok(PerVarDistribsConfig {
        stats: vec![
            PerVarStat::ObsHet,
            PerVarStat::Maf,
            PerVarStat::ExpHet,
            PerVarStat::UnbiasedExpHet,
            PerVarStat::PolyVarsRatio,
            PerVarStat::MissingRate,
        ],
        pops: Pops::all(num_individuals),
        bins: HistBins::linear(start, end, DEFAULT_NUM_BINS)?,
        obs_het: ObsHet::new(DEFAULT_MIN_NUM_INDIVIDUALS),
        maf: Maf::new(ploidy, DEFAULT_MIN_NUM_INDIVIDUALS)?,
        exp_het: ExpHet::of_the_exponent_asked_for(None, ploidy, DEFAULT_MIN_NUM_INDIVIDUALS)?,
        poly_threshold: DEFAULT_POLY_THRESHOLD,
    })
}

/// What a pass gave that the other has to give too: the variants of each
/// of the three and the density itself. The values of the distributions
/// and of the rates are compared to the bit by the tests of
/// `calc_variants_summary` and not here.
#[derive(Debug, PartialEq)]
struct Gave {
    per_var_num_vars: u64,
    per_individual_num_vars: u64,
    density: VarDensity,
}

/// The three, each in a pass of its own over a reader of its own, timed
/// from building the first reader to the last result.
///
/// The phases it gives are those of each of the three passes, in that
/// order, when the cargo feature `bench-phases` is on.
fn the_three(path: &Path) -> Result<(Duration, Gave, String), popnei::Error> {
    let _ = the_phases_of_the_pass();
    let started = Instant::now();
    let mut reader = reader_of(path)?;
    let config = the_default_distribs(reader.ploidy(), reader.individuals().len())?;
    let distribs = calc_per_var_distribs(&mut *reader, &config)?;
    let of_the_distribs = the_phases_of_the_pass();
    let mut reader = reader_of(path)?;
    let rates = calc_per_individual_stats(&mut *reader)?;
    let of_the_rates = the_phases_of_the_pass();
    let mut reader = reader_of(path)?;
    let density = calc_var_density(&mut *reader, WINDOW_SIZE, None)?;
    let took = started.elapsed();
    let of_the_density = the_phases_of_the_pass();
    Ok((
        took,
        Gave {
            per_var_num_vars: distribs.num_vars,
            per_individual_num_vars: rates.num_vars(),
            density,
        },
        format!("{of_the_distribs}{of_the_rates}{of_the_density}"),
    ))
}

/// The same three in the one pass of `calc_variants_summary`, timed from
/// building its reader to its result.
fn the_one(path: &Path) -> Result<(Duration, VariantsSummary, String), popnei::Error> {
    let _ = the_phases_of_the_pass();
    let started = Instant::now();
    let mut reader = reader_of(path)?;
    let config = VariantsSummaryConfig {
        per_var: Some(the_default_distribs(
            reader.ploidy(),
            reader.individuals().len(),
        )?),
        per_individual: true,
        density: Some(VarDensityConfig {
            window_size: WINDOW_SIZE,
            chrom_lengths: None,
        }),
    };
    let summary = calc_variants_summary(&mut *reader, &config, &mut nothing_after_a_block)?;
    let took = started.elapsed();
    Ok((took, summary, the_phases_of_the_pass()))
}

/// The times of the three and of the one over one file, a pair a run.
struct TimesOfAFile {
    three: Vec<Duration>,
    one: Vec<Duration>,
}

/// One run over the file at `path`: the three, then the one, with a check
/// that both gave the same variants and the same density.
fn one_run(path: &Path) -> Result<Run, String> {
    let said = |error: popnei::Error| format!("{path}: {error}", path = path.display());
    let (three, gave_three, phases_of_the_three) = the_three(path).map_err(said)?;
    let (one, summary, phases_of_the_one) = the_one(path).map_err(said)?;
    let (Some(distribs), Some(rates), Some(density)) =
        (summary.per_var, summary.per_individual, summary.density)
    else {
        return Err(format!(
            "{path}: calc_variants_summary left out a statistic it was asked for",
            path = path.display(),
        ));
    };
    let gave_one = Gave {
        per_var_num_vars: distribs.num_vars,
        per_individual_num_vars: rates.num_vars(),
        density,
    };
    if gave_three != gave_one {
        return Err(format!(
            "{path}: the three gave {three_vars} and {three_rates} variants and {three_windows} \
             windows, and the one {one_vars} and {one_rates} variants and {one_windows} windows",
            path = path.display(),
            three_vars = gave_three.per_var_num_vars,
            three_rates = gave_three.per_individual_num_vars,
            three_windows = gave_three.density.num_windows(),
            one_vars = gave_one.per_var_num_vars,
            one_rates = gave_one.per_individual_num_vars,
            one_windows = gave_one.density.num_windows(),
        ));
    }
    Ok(Run {
        three,
        one,
        num_vars: gave_one.per_var_num_vars,
        phases_of_the_three,
        phases_of_the_one,
    })
}

/// One run over one file: the time of the three and of the one, the
/// variants they gave, and the phases of their passes.
struct Run {
    three: Duration,
    one: Duration,
    num_vars: u64,
    phases_of_the_three: String,
    phases_of_the_one: String,
}

/// The two clocks of the phases of a pass, as a piece of the line of a
/// run: how long the pass was inside `next_block`, waiting for the thread
/// that reads one block ahead, and how long it was working on the blocks.
/// Taking them zeroes them, so each pass gives its own.
///
/// It is the cargo feature `bench-phases` of the core crate, as in
/// `kinship.rs`: `cargo bench --features bench-phases --bench
/// variants_summary` turns it on, and without it the passes hold no clock
/// and there is nothing to print.
#[cfg(feature = "bench-phases")]
fn the_phases_of_the_pass() -> String {
    let phases = popnei::phases::taken();
    format!(
        " [next_block {next_block:.3} s, work {work:.3} s]",
        next_block = phases.next_block.as_secs_f64(),
        work = phases.work.as_secs_f64(),
    )
}

/// Nothing, which is what the phases of a pass are when the cargo feature
/// `bench-phases` is off.
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

/// The best, the median and the worst of `times`.
fn spread(times: &[Duration]) -> String {
    format!(
        "best {best}, median {median}, worst {worst}",
        best = seconds(sorted_time(times, 0)),
        median = seconds(median_time(times)),
        worst = seconds(sorted_time(times, times.len().saturating_sub(1))),
    )
}

/// The benchmark reads files of the disc, which wasm has not. This is what
/// `cargo check --target wasm32-unknown-unknown --all-targets` compiles of
/// it.
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
    if let Err(error) = rayon::ThreadPoolBuilder::new()
        .num_threads(arguments.threads)
        .build_global()
    {
        eprintln!(
            "rayon's global pool of {} threads: {error}",
            arguments.threads
        );
        return ExitCode::FAILURE;
    }
    println!(
        "{runs} runs, rayon's global pool of {threads} threads, windows of {WINDOW_SIZE} base \
         pairs",
        runs = arguments.runs,
        threads = arguments.threads,
    );
    let mut times: Vec<TimesOfAFile> = arguments
        .paths
        .iter()
        .map(|_| TimesOfAFile {
            three: Vec::with_capacity(arguments.runs),
            one: Vec::with_capacity(arguments.runs),
        })
        .collect();
    for run in 0..=arguments.runs {
        for (path, of_the_file) in arguments.paths.iter().zip(times.iter_mut()) {
            let Run {
                three,
                one,
                num_vars,
                phases_of_the_three,
                phases_of_the_one,
            } = match one_run(path) {
                Ok(done) => done,
                Err(problem) => {
                    eprintln!("{problem}");
                    return ExitCode::FAILURE;
                }
            };
            let path = path.display();
            if run == 0 {
                println!(
                    "the first run, which is not timed, {path}: {num_vars} variants, \
                     the three {three}, the one {one}",
                    three = seconds(three),
                    one = seconds(one),
                );
                continue;
            }
            println!(
                "run {run}, {path}: the three {three}{phases_of_the_three}, \
                 the one {one}{phases_of_the_one}",
                three = seconds(three),
                one = seconds(one),
            );
            of_the_file.three.push(three);
            of_the_file.one.push(one);
        }
    }
    for (path, of_the_file) in arguments.paths.iter().zip(&times) {
        let three = median_time(&of_the_file.three).as_secs_f64();
        let one = median_time(&of_the_file.one).as_secs_f64();
        let saving = match three {
            0.0 => f64::NAN,
            _ => 100.0 * (three - one) / three,
        };
        println!(
            "{path}\n  the three: {of_three}\n  the one:   {of_one}\n  \
             the one saves {saving:.1} in 100 of the median of the three",
            path = path.display(),
            of_three = spread(&of_the_file.three),
            of_one = spread(&of_the_file.one),
        );
    }
    ExitCode::SUCCESS
}
