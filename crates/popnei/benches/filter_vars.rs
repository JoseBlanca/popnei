//! What a filter of the variants costs a whole pass over a dataset, on a
//! given number of threads.
//!
//! It times one thing: a pass with the genotypes alone asked for, from
//! opening the file to the last block, over a VCF or over a vars file.
//! With `--max-missing-rate` the pass keeps the variants whose missing
//! genotypes, divided by the individuals of the dataset, are at most that
//! rate. With `--max-ld-r2` and `--max-dist`, which are given together, it
//! keeps the variants that do not repeat what a variant kept near them
//! already said: r² is the squared correlation between the dosages of two
//! variants, and a variant is dropped when its r² is above that number
//! against any variant the filter has kept within that many base pairs
//! behind it on its chromosome. The two filters may be given together or
//! apart, and with neither the same pass runs with no filter, as it does
//! today. Each of them goes on the pass through `chain_of`, which is what
//! builds the chain of readers of a pass for a Python and a TypeScript
//! user, so what is timed is the pass that `filter_by_missing_data` and
//! `filter_by_ld` give them; with both, the filter of missing data is the
//! one nearest the source. The cost of a filter is the difference of the
//! two medians, so the two are run back to back on a machine that is doing
//! nothing else: for the filter of missing data the difference is a tenth
//! of what each of them takes, and the load of the machine moves a median
//! by more than that.
//!
//! A time for the filter by linkage disequilibrium says nothing on its
//! own, because what it costs is set by how many variants its window
//! holds, and that is of the dataset and not of the filter: the window of
//! a variant is the variants the filter has already kept that are within
//! `--max-dist` base pairs behind it on its chromosome, and what the
//! filter computes is the r² of the variant against every one of them. So
//! a pass with that filter prints what the window held over the pass, the
//! mean and the largest of the variants that lie behind a kept variant and
//! within that distance of it. It is counted outside the clock, in the
//! pass that is not timed, from the chromosome and the position of every
//! variant that pass kept: the library is asked for nothing it does not
//! give already, and no timed run counts any of it.
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
//! cargo bench --bench filter_vars -- <path> --threads 1 --runs 5 \
//!     --max-ld-r2 0.3 --max-dist 100000
//! ```
//!
//! cargo runs it with `crates/popnei` as its working directory, so a
//! relative path is read from there and an absolute one is the plainer
//! thing to give.
//!
//! `--threads` is 1 and `--runs` is 5 when they are not given. It prints
//! the wall time of each run, with the variants the pass gave, the alleles
//! of their genotypes and what each filter was given and kept, and then the
//! best, the median and the worst of the times; with an even number of runs
//! the median is the middle of the two middle times. The median is what the
//! numbers of "Speed" of `docs/specs/filters.md` are taken from, and the
//! best and the worst say how much the machine was doing something else.
//!
//! The alleles are added up inside the clock in both passes, so that a
//! pass with no filter, which reads no genotype of its own, cannot be
//! shorter than one because the genotypes were never there: over these two
//! files a pass that keeps every variant gives 200000000 of them.
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
use popnei::filters::{PassStep, VarFilteringCriterion, chain_of};
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

/// What the command line asked for. A threshold is `None` when the command
/// line did not name it, and with none of them the pass runs with no
/// filter. `max_ld_r2` and `max_dist` are both there or neither is: one
/// without the other names no filter.
struct Arguments {
    path: PathBuf,
    threads: usize,
    runs: usize,
    max_missing_rate: Option<f64>,
    max_ld_r2: Option<f64>,
    max_dist: Option<u64>,
}

impl Arguments {
    /// The filters of the pass, in the order they go on it: the filter of
    /// missing data nearest the source, and the filter by linkage
    /// disequilibrium over it, which then works its r² out over the
    /// variants the first one left.
    fn steps(&self) -> Vec<PassStep> {
        let mut steps = Vec::new();
        if let Some(rate) = self.max_missing_rate {
            steps.push(PassStep::VarFilter(VarFilteringCriterion::MaxMissingRate(
                rate,
            )));
        }
        if let (Some(max_allowed_r2), Some(max_dist)) = (self.max_ld_r2, self.max_dist) {
            steps.push(PassStep::VarFilter(VarFilteringCriterion::MaxLdR2 {
                max_allowed_r2,
                max_dist,
            }));
        }
        steps
    }

    /// What the line of every run says the pass was: the filters it went
    /// through with their thresholds, or that it had none.
    fn what_the_pass_has(&self) -> String {
        let mut filters = Vec::new();
        if let Some(rate) = self.max_missing_rate {
            filters.push(format!("the missing data filter at {rate}"));
        }
        if let (Some(max_allowed_r2), Some(max_dist)) = (self.max_ld_r2, self.max_dist) {
            filters.push(format!(
                "the filter by linkage disequilibrium at {max_allowed_r2} over {max_dist} base pairs"
            ));
        }
        match filters.is_empty() {
            true => "no filter".to_owned(),
            false => filters.join(" and "),
        }
    }
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
filter_vars <path to a VCF or a vars file> [--threads n] [--runs n]
            [--max-missing-rate r] [--max-ld-r2 r --max-dist d]

It times whole passes over that file with the genotypes alone asked for:
opening the file, its header, every block to the end of it. With
--max-missing-rate the pass keeps the variants whose missing genotypes,
divided by the individuals of the dataset, are at most r. With
--max-ld-r2 and --max-dist it keeps the variants whose r², the squared
correlation between the dosages of two variants, is at most r against
every variant it has kept within d base pairs behind them on their
chromosome. The two may be given together or apart; with neither the same
pass runs with no filter, and the difference of the two medians is what
the filter costs.

  --threads n              how many threads the pool it reads in has, 1 by default
  --runs n                 how many times it reads the file, 5 by default
  --max-missing-rate r     the largest rate of missing genotypes that keeps
                           a variant, a number from 0 to 1
  --max-ld-r2 r            the largest r² against a variant of its window that
                           keeps a variant, a number from 0 to 1; it takes
                           --max-dist with it
  --max-dist d             how many base pairs behind a variant its window
                           reaches, 1 or more; it takes --max-ld-r2 with it
  --help                   this

A path that ends in `.vars` is read as a vars file and anything else as a
VCF. One pass that is not timed comes first, so that the timed runs pay
neither the page faults of the first touch of the memory a pass works in
nor a read of the disc. It prints the wall time of each run, with the
variants the pass gave, the alleles of their genotypes, which are added
up inside the clock so that a pass with no filter cannot be short because
the genotypes were never filled, and what each filter was given and kept;
and then the best, the median and the worst of the times. The median is
what `docs/specs/filters.md` states, and the best and the worst say how
much the machine was doing something else.

With the filter by linkage disequilibrium it also prints what the window
held, the mean and the largest of the kept variants that lie within
--max-dist base pairs behind a kept variant, since what that filter costs
is set by how many variants its window holds and that is of the dataset.
It is counted in the pass that is not timed, from the chromosome and the
position of the variants that pass kept, and no timed run counts any of
it.";

/// What to run, or the message that says what the command line should have
/// been.
fn arguments() -> Result<Arguments, String> {
    let mut path: Option<PathBuf> = None;
    let mut threads = DEFAULT_THREADS;
    let mut runs = DEFAULT_RUNS;
    let mut max_missing_rate: Option<f64> = None;
    let mut max_ld_r2: Option<f64> = None;
    let mut max_dist: Option<u64> = None;
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
            "--max-missing-rate" | "--max-ld-r2" => {
                let name = arg.clone();
                let number = args
                    .next()
                    .ok_or_else(|| format!("{name} takes a number and none came after it"))?
                    .parse::<f64>()
                    .map_err(|_| format!("{name} takes a number"))?;
                if name == "--max-missing-rate" {
                    max_missing_rate = Some(number);
                } else {
                    max_ld_r2 = Some(number);
                }
            }
            "--max-dist" => {
                max_dist = Some(
                    args.next()
                        .ok_or_else(|| {
                            "--max-dist takes a number and none came after it".to_owned()
                        })?
                        .parse::<u64>()
                        .map_err(|_| "--max-dist takes a number".to_owned())?,
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
    // The two arguments of the filter by linkage disequilibrium are one
    // filter: with one of them alone the pass would run with no filter of
    // that kind and the run would be reported as the pass that has it.
    if max_ld_r2.is_some() != max_dist.is_some() {
        return Err("--max-ld-r2 and --max-dist are given together".to_owned());
    }
    Ok(Arguments {
        path,
        threads,
        runs,
        max_missing_rate,
        max_ld_r2,
        max_dist,
    })
}

/// One run: how long the pass took and the line that says what it gave, the
/// variants and, when there was a filter, what it was given and kept.
struct Run {
    took: Duration,
    did: String,
    /// Where each variant the pass gave lies, which is filled in the pass
    /// that is not timed and is empty in every other: it is what the
    /// window of the filter by linkage disequilibrium is counted from.
    kept: Vec<TheVariantKept>,
}

/// Where a variant the pass gave lies, which is what says whether it is in
/// the window of a variant after it.
struct TheVariantKept {
    /// The number of its chromosome in the table of the reader.
    chrom: u32,
    /// Its position, 1 based as in a VCF.
    pos: u64,
}

/// Whether the pass collects where each variant it gives lies. It is
/// collected in the pass that is not timed and in no other, so that
/// nothing of the window is inside a clock.
#[derive(Clone, Copy, PartialEq, Eq)]
enum TheVariantsKept {
    /// Where each variant the pass gives lies is kept.
    Collected,
    /// The pass gives its variants and where they lie is not kept.
    NotCollected,
}

/// What the window of the filter by linkage disequilibrium held over a
/// pass: for each variant the pass gave, how many variants it gave lie
/// within `max_dist` base pairs behind it on its chromosome.
struct TheWindowHeld {
    /// The mean of those counts over the variants of the pass.
    mean: f64,
    /// The largest of them.
    largest: usize,
    /// How many variants the counts are of, which is what the pass gave.
    of_the_vars: usize,
}

/// What the window held over a pass whose variants are `kept`, in the order
/// the pass gave them, or `None` when the pass gave no variant.
///
/// The filter refuses a variant that does not come after the one before it,
/// so the variants of a chromosome come together and in the order of their
/// positions: the variants within `max_dist` behind the one at a place are
/// the ones from the first that is near enough up to that place, and one
/// index that never goes back walks the whole pass.
///
/// This is not what the filter computes. The filter drops a variant as soon
/// as one variant of its window is too close a match to it, so it works out
/// fewer values of r² than these counts; what the counts say is how many
/// variants the window held, which is what a time for the filter has to be
/// read against.
fn the_window_over(kept: &[TheVariantKept], max_dist: u64) -> Option<TheWindowHeld> {
    let of_the_vars = kept.len();
    if of_the_vars == 0 {
        return None;
    }
    let mut first_of_the_window = 0_usize;
    let mut held_together: u64 = 0;
    let mut largest = 0_usize;
    for (place, variant) in kept.iter().enumerate() {
        while first_of_the_window < place {
            let Some(behind) = kept.get(first_of_the_window) else {
                break;
            };
            if behind.chrom == variant.chrom && variant.pos.saturating_sub(behind.pos) <= max_dist {
                break;
            }
            first_of_the_window = first_of_the_window.saturating_add(1);
        }
        let held = place.saturating_sub(first_of_the_window);
        // A pass of more variants than a u64 counts is one no file holds.
        held_together = held_together.saturating_add(u64::try_from(held).unwrap_or(u64::MAX));
        largest = largest.max(held);
    }
    Some(TheWindowHeld {
        mean: held_together as f64 / of_the_vars as f64,
        largest,
        of_the_vars,
    })
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
/// for, and the filters of `steps` over it, timed from the building of the
/// reader to the last block.
///
/// The alleles of every block are added up inside the clock, and the run
/// prints how many there were. What the filter costs is this pass less the
/// pass with no filter, and the pass with no filter never looks at a
/// genotype: a reader that stopped filling the column, or filled it only
/// when somebody read it, would make that pass shorter and the filter look
/// dearer, with nothing else to show it. The alleles of these two files are
/// 200000000 in a pass that keeps every variant, and that number is what
/// says the genotypes were there.
///
/// The counts of the filter are read after the clock stops: they are two
/// numbers of a chain that the pass has already built, and the line they go
/// into is printed and not timed.
///
/// With `TheVariantsKept::Collected` the chromosome and the position of
/// every variant the pass gives are kept, which is what the window of the
/// filter by linkage disequilibrium is counted from afterwards. The
/// collecting is inside the clock, so it is asked for in the pass that is
/// not timed and in no other. The filter by linkage disequilibrium asks its source for
/// the two columns whatever the consumer asked for, so they are in the
/// blocks of a pass that has it and in no block of a pass that has not.
fn one_pass(
    path: &Path,
    steps: &[PassStep],
    collect: TheVariantsKept,
) -> Result<Run, popnei::Error> {
    let started = Instant::now();
    let mut reader = chain_of(reader_of(path)?, steps)?;
    reader.set_needs(Needs::GTS);
    let mut variants: u64 = 0;
    let mut alleles: u64 = 0;
    let mut kept: Vec<TheVariantKept> = Vec::new();
    while let Some(block) = reader.next_block()? {
        // A file of more variants than a u64 counts cannot be written.
        variants = variants.saturating_add(u64::try_from(block.num_vars).unwrap_or(u64::MAX));
        // Nor one of more alleles: this machine does not address them.
        alleles = alleles.saturating_add(u64::try_from(block.gts.len()).unwrap_or(u64::MAX));
        if collect == TheVariantsKept::NotCollected {
            continue;
        }
        // A block of a pass with no filter by linkage disequilibrium holds
        // neither column, and then there is no window to count.
        if let (Some(chroms), Some(poss)) = (block.chrom.as_ref(), block.pos.as_ref()) {
            kept.extend(
                chroms
                    .iter()
                    .zip(poss.iter())
                    .map(|(chrom, pos)| TheVariantKept {
                        chrom: *chrom,
                        pos: *pos,
                    }),
            );
        }
    }
    let took = started.elapsed();
    let mut did = format!("{variants} variants, {alleles} alleles");
    for (kind, counts) in reader.filtering_stats() {
        did.push_str(&format!(
            ", the {kind} filter was given {given} and kept {it_kept}",
            given = counts.vars_processed,
            it_kept = counts.vars_kept,
        ));
    }
    Ok(Run { took, did, kept })
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
        filter = arguments.what_the_pass_has(),
    );
    let steps = arguments.steps();
    let mut times = Vec::with_capacity(arguments.runs);
    // The pass that is not timed, and then the timed ones. Both go through
    // `install`, so that the one that warms the memory runs on the same
    // pool as the ones that are timed. The first is the one that collects
    // where each variant it gave lies, which the window is counted from
    // and which no timed pass does.
    for run in 0..=arguments.runs {
        let collect = match run {
            0 => TheVariantsKept::Collected,
            _ => TheVariantsKept::NotCollected,
        };
        let done = match pool.install(|| one_pass(&arguments.path, &steps, collect)) {
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
            if let Some(max_dist) = arguments.max_dist {
                match the_window_over(&done.kept, max_dist) {
                    Some(window) => println!(
                        "the window of that pass, which is not timed either: {mean:.1} variants \
                         on average and {largest} at most, of the {of_the_vars} variants it \
                         gave, within {max_dist} base pairs behind each of them",
                        mean = window.mean,
                        largest = window.largest,
                        of_the_vars = window.of_the_vars,
                    ),
                    None => println!("the window of that pass: it gave no variant"),
                }
            }
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
