//! How long each of the four passes that turn a row of a block into
//! standardized dosages takes, over one real block of variants. The
//! principal components of the variants and the kinship both walk that
//! pass; the divisor timed here is the one of the components.
//!
//! Standardizing a row of a block turns the genotypes of one variant into
//! one number per individual, its dosage, centered and divided by the
//! standard deviation of the variant, and
//! `docs/reports/perf-pca-2026-09-22.md` measured it at 20.6 ms for a
//! block of 5000 variants of 1000 individuals, 51.6 per 100 of an analysis
//! of 100000 variants. It is made of four passes over the row: counting
//! the alleles to find the one called most often, writing the dosage of
//! each genotype into a byte, counting how many genotypes have each
//! dosage, and looking each of those bytes up in the value it stands for.
//! The compiler inlines all four into one closure, so a sampling profile
//! of the whole analysis sees them as one frame and the split of that
//! report is a count of instructions and not a time.
//!
//! This benchmark times each of the four on its own, and then the whole
//! row, over the rows of one block read from a vars file. What the four
//! come to is printed against what the whole row takes: they should be
//! close to it, and what is left is the work of a row that is not one of
//! the four, the major allele, the count of the different alleles and the
//! mean and the deviation of the dosages.
//!
//! It reads one block and refuses a file whose first block is not 5000
//! variants of 1000 individuals, which is the block of that report, so
//! that the times are of the dataset the review is about. It runs on one
//! thread and builds no pool: nothing here calls rayon, and no product and
//! no eigendecomposition is taken, so no BLAS is called either. The two
//! variables that give the threads of the whole analysis are printed with
//! each run, because a number taken with the wrong ones is what this
//! measurement can get wrong without showing it:
//!
//! ```text
//! VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1 \
//!     cargo bench --features bench-internals --bench standardize_row -- \
//!     /Users/jose/devel/popnei-bench/big.vars --runs 5
//! ```
//!
//! The cargo feature is what makes the four reachable: they are private
//! functions of the module `variant`, which the principal components and
//! the kinship both walk, and a benchmark is a crate of its own, so
//! `variant::bench_internals`, which the feature turns on, is what
//! re-exports them. cargo runs the benchmark with `crates/popnei` as its
//! working directory, so a relative path is read from there and an
//! absolute one is the plainer thing to give.
//!
//! One run that is not timed comes before the timed ones. It pays the page
//! faults of the first touch of the buffers a pass writes into, and it is
//! the run that adds up the checksum of each pass, which is printed and
//! which the timed runs do not compute: adding up the thousand values of a
//! row costs as much as the pass that writes them. What each checksum is:
//! for the counting of the alleles, the called alleles of every row; for
//! the codes, every code; for the counts of the codes, every count; and
//! for the two that write the standardized values, the cube of every value
//! of every row that had variance, which makes those last two equal when
//! the lookup and the whole row agree. The cube and not the square,
//! because the squares of a standardized row add to the individuals of the
//! block whatever its dosages are, so a checksum of squares is a number
//! that cannot change and says nothing. A checksum that changes is a
//! change of what the code computes, and it is also what keeps the
//! compiler from dropping a pass whose result nothing reads; the timed
//! runs, which compute none, are held by `black_box` on the buffer each
//! pass writes.
//!
//! `--runs` is 5 when it is not given. It prints the best, the median and
//! the worst of the runs of each pass, in milliseconds for the whole
//! block. The best is the number to state, since every other process on
//! the machine can only make a run longer.
//!
//! The file is the one of `docs/reports/pca-measurement.md`, written there
//! from a VCF of 100000 variants of 1000 individuals whose genotypes are
//! missing at a rate of 0.03; its first block is the one read here.

#![cfg_attr(
    target_family = "wasm",
    allow(
        dead_code,
        unused_imports,
        reason = "the benchmark reads a file of the disc: in wasm only its empty main is \
                  compiled, and what the timing is made of is left unused"
    )
)]

use std::hint::black_box;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{Duration, Instant};

use popnei::block::BlockReader;
use popnei::io::vars::VarsReader;
use popnei::variant::bench_internals::{
    Dosages, Scratch, the_codes_of_the_genotypes, the_counts_of_the_codes, the_standardized_row,
    the_standardized_values,
};
use popnei::variant::{AlleleCounts, Needs, count_alleles, the_major_allele};

/// How many times each pass is timed when the command line does not say.
const DEFAULT_RUNS: usize = 5;

/// The variants the block has to hold, which is the block of
/// `docs/reports/perf-pca-2026-09-22.md`.
const VARIANTS_OF_THE_BLOCK: usize = 5000;

/// The individuals the block has to hold, which is the block of that same
/// report.
const INDIVIDUALS_OF_THE_BLOCK: usize = 1000;

/// What the command line asked for.
struct Arguments {
    path: PathBuf,
    runs: usize,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
standardize_row <path to a vars file> [--runs n]

It reads the first block of that file, which has to be 5000 variants of
1000 individuals, and times over its rows, each on its own: the counting of
the alleles of a variant, the writing of the code of each genotype, the
counting of the codes, the lookup that writes the standardized value of
each code, and a whole standardized row, which is the four together.

  --runs n    how many times each pass is timed, 5 by default
  --help      this

One run that is not timed comes first, and it is the one that adds up the
checksum of each pass. It prints the best, the median and the worst of the
timed runs of each pass, in milliseconds for the whole block, and what the
four passes come to against what a whole row takes.

It needs the cargo feature `bench-internals`, which is what makes the four
passes reachable from outside the module that holds them:

    cargo bench --features bench-internals --bench standardize_row -- <path>";

/// What to run, or the message that says what the command line should have
/// been.
fn arguments() -> Result<Arguments, String> {
    let mut path: Option<PathBuf> = None;
    let mut runs = DEFAULT_RUNS;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--runs" => {
                runs = args
                    .next()
                    .ok_or_else(|| "--runs takes a number and none came after it".to_owned())?
                    .parse::<usize>()
                    .map_err(|_| "--runs takes a number of 1 or more".to_owned())?;
            }
            // `cargo bench` adds this to the command line of every bench,
            // to tell a harness that has tests too to run its benchmarks.
            // This one has only this benchmark and takes it as nothing.
            "--bench" => {}
            "--help" | "-h" => return Err(USAGE.to_owned()),
            // An argument that is not one of those and looks like one is
            // refused instead of being taken for the path of the file, as
            // `pca_vars.rs` refuses it.
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
    Ok(Arguments { path, runs })
}

/// The one block the passes run over, and what the passes that do not
/// start from the genotypes need: the block is read once and these are
/// worked out once, so that each timed pass does its own work and no
/// other.
struct TheBlock {
    /// The genotypes of the block, variant after variant, `ploidy`
    /// alleles for each individual.
    gts: Vec<i8>,
    /// The individuals of the block, which is the length of one row of
    /// standardized values.
    num_individuals: usize,
    /// The alleles of the genotype of one individual.
    ploidy: usize,
    /// The ploidy again, as the loop over the genotypes of a row takes it.
    of_a_genotype: NonZeroUsize,
    /// The individuals times the ploidy, which is one row of `gts`.
    alleles_of_a_row: usize,
    /// The dosages a genotype can have, 0 to the ploidy.
    num_dosages: usize,
    /// The allele called most often in each variant, which the pass that
    /// writes the codes is given.
    majors: Vec<i8>,
    /// The code of every genotype of the block, `num_individuals` of them
    /// for each variant, which the pass that counts the codes and the one
    /// that looks them up are given.
    codes: Vec<u8>,
    /// How many genotypes of each variant have each dosage, `num_dosages`
    /// of them for each variant, which the pass that looks the codes up is
    /// given.
    dosage_counts: Vec<u32>,
    /// What the analysis of the variants asks of a row, which a whole row
    /// takes: no variant of more than two alleles is turned into a
    /// biallelic one, and the dosages are divided by their own standard
    /// deviation, which is the divisor of the principal components.
    options: Dosages,
}

/// The first block of the vars file at `path`, with the genotypes alone
/// asked for, and everything the passes need worked out from it.
///
/// A file whose first block is not 5000 variants of 1000 individuals is
/// refused: the times this prints are of that block and of no other.
#[cfg(not(target_family = "wasm"))]
fn the_block_of(path: &Path) -> Result<TheBlock, String> {
    let named = |error: popnei::Error| format!("{path}: {error}", path = path.display());
    let mut reader = VarsReader::from_path(path).map_err(named)?;
    reader.set_needs(Needs::GTS);
    let Some(block) = reader.next_block().map_err(named)? else {
        return Err(format!(
            "{path}: the file has no block of variants",
            path = path.display()
        ));
    };
    if block.num_vars != VARIANTS_OF_THE_BLOCK || block.num_individuals != INDIVIDUALS_OF_THE_BLOCK
    {
        return Err(format!(
            "{path}: its first block is {num_vars} variants of {num_individuals} individuals, \
             and this benchmark times {VARIANTS_OF_THE_BLOCK} variants of \
             {INDIVIDUALS_OF_THE_BLOCK} individuals",
            path = path.display(),
            num_vars = block.num_vars,
            num_individuals = block.num_individuals,
        ));
    }
    let num_individuals = block.num_individuals;
    let ploidy = block.ploidy;
    let Some(of_a_genotype) = NonZeroUsize::new(ploidy) else {
        return Err(format!(
            "{path}: its block has a ploidy of 0",
            path = path.display()
        ));
    };
    let too_large = || {
        format!(
            "{path}: its block is too large for this machine to index",
            path = path.display()
        )
    };
    let alleles_of_a_row = num_individuals.checked_mul(ploidy).ok_or_else(too_large)?;
    let num_dosages = ploidy.checked_add(1).ok_or_else(too_large)?;
    let num_codes = block
        .num_vars
        .checked_mul(num_individuals)
        .ok_or_else(too_large)?;
    let num_counts = block
        .num_vars
        .checked_mul(num_dosages)
        .ok_or_else(too_large)?;
    let mut majors = Vec::with_capacity(block.num_vars);
    let mut codes = vec![0_u8; num_codes];
    let mut dosage_counts = vec![0_u32; num_counts];
    let mut allele_counts: AlleleCounts = [0; 128];
    let mut counts_of_a_row = [0_u32; 255];
    for ((gts, codes_of_a_row), counts_written) in block
        .gts
        .chunks_exact(alleles_of_a_row)
        .zip(codes.chunks_exact_mut(num_individuals))
        .zip(dosage_counts.chunks_exact_mut(num_dosages))
    {
        count_alleles(gts, &mut allele_counts).map_err(named)?;
        let major = the_major_allele(&allele_counts);
        majors.push(major);
        the_codes_of_the_genotypes(gts, of_a_genotype, major, codes_of_a_row);
        the_counts_of_the_codes(codes_of_a_row, num_dosages, &mut counts_of_a_row);
        for (target, count) in counts_written.iter_mut().zip(counts_of_a_row.iter()) {
            *target = *count;
        }
    }
    Ok(TheBlock {
        gts: block.gts,
        num_individuals,
        ploidy,
        of_a_genotype,
        alleles_of_a_row,
        num_dosages,
        majors,
        codes,
        dosage_counts,
        options: Dosages::of_the_principal_components(),
    })
}

/// One pass of the counting of the alleles over every row of the block,
/// which is what finds the allele called most often.
///
/// With `CHECKSUM` it adds up the called alleles of every row, and without
/// it that addition is not compiled at all, so a timed run pays nothing
/// for it.
///
/// # Errors
///
/// What the counting of the alleles of a variant refuses.
#[cfg(not(target_family = "wasm"))]
fn the_pass_of_count_alleles<const CHECKSUM: bool>(
    block: &TheBlock,
    counts: &mut AlleleCounts,
) -> Result<f64, String> {
    let mut checksum = 0.0_f64;
    for gts in block.gts.chunks_exact(block.alleles_of_a_row) {
        let called = count_alleles(gts, counts).map_err(|error| format!("{error}"))?;
        black_box(&*counts);
        if CHECKSUM {
            checksum += f64::from(called);
        }
    }
    Ok(checksum)
}

/// One pass of the writing of the code of each genotype over every row of
/// the block, given the allele called most often in each, which the block
/// was read with once.
#[cfg(not(target_family = "wasm"))]
fn the_pass_of_the_codes<const CHECKSUM: bool>(block: &TheBlock, codes: &mut [u8]) -> f64 {
    let mut checksum = 0.0_f64;
    for (gts, major) in block
        .gts
        .chunks_exact(block.alleles_of_a_row)
        .zip(block.majors.iter())
    {
        the_codes_of_the_genotypes(gts, block.of_a_genotype, *major, codes);
        black_box(&*codes);
        if CHECKSUM {
            for code in codes.iter() {
                checksum += f64::from(*code);
            }
        }
    }
    checksum
}

/// One pass of the counting of the codes over the codes of every row of
/// the block, which the block was read with once.
#[cfg(not(target_family = "wasm"))]
fn the_pass_of_the_counts<const CHECKSUM: bool>(block: &TheBlock, counts: &mut [u32; 255]) -> f64 {
    let mut checksum = 0.0_f64;
    for codes in block.codes.chunks_exact(block.num_individuals) {
        the_counts_of_the_codes(codes, block.num_dosages, counts);
        black_box(&*counts);
        if CHECKSUM {
            for count in counts.iter().take(block.num_dosages) {
                checksum += f64::from(*count);
            }
        }
    }
    checksum
}

/// One pass of the lookup that writes the standardized value of each code
/// over every row of the block.
///
/// The counts of the codes of the row are copied into `counts` first,
/// which is `num_dosages` entries and not a pass over the row. The mean
/// and the deviation of the dosages and the table of the value of each
/// dosage are inside the timing, because they are what the tail of a
/// standardized row is made of beside the lookup itself, and they are
/// `num_dosages` entries against the thousand the lookup writes.
#[cfg(not(target_family = "wasm"))]
fn the_pass_of_the_values<const CHECKSUM: bool>(
    block: &TheBlock,
    counts: &mut [u32; 255],
    values: &mut [f64; 256],
    row: &mut [f64],
) -> f64 {
    let mut checksum = 0.0_f64;
    for (codes, counts_of_the_row) in block
        .codes
        .chunks_exact(block.num_individuals)
        .zip(block.dosage_counts.chunks_exact(block.num_dosages))
    {
        for (target, count) in counts.iter_mut().zip(counts_of_the_row.iter()) {
            *target = *count;
        }
        let used =
            the_standardized_values(counts, block.ploidy, codes, values, row, &block.options);
        black_box(&*row);
        if CHECKSUM && used {
            for value in row.iter() {
                checksum += *value * *value * *value;
            }
        }
    }
    checksum
}

/// One pass of a whole standardized row over every row of the block, which
/// is the four passes above and the work between them.
///
/// # Errors
///
/// What standardizing a row refuses: a variant of more than two alleles
/// among them.
#[cfg(not(target_family = "wasm"))]
fn the_pass_of_the_standardized_row<const CHECKSUM: bool>(
    block: &TheBlock,
    scratch: &mut Scratch,
    row: &mut [f64],
) -> Result<f64, String> {
    let mut checksum = 0.0_f64;
    for (position, gts) in block.gts.chunks_exact(block.alleles_of_a_row).enumerate() {
        let used = the_standardized_row(gts, block.ploidy, position, &block.options, scratch, row)
            .map_err(|error| format!("{error}"))?;
        black_box(&*row);
        if CHECKSUM && used {
            for value in row.iter() {
                checksum += *value * *value * *value;
            }
        }
    }
    Ok(checksum)
}

/// What the runs of one pass gave: how long each timed run took, and the
/// checksum of the one run that was not timed.
#[cfg(not(target_family = "wasm"))]
struct Runs {
    times: Vec<Duration>,
    checksum: f64,
}

/// One pass run `runs` times with a clock on it, after one run with no
/// clock that computes its checksum.
///
/// `pass` is given true for the run that is to compute the checksum and
/// false for a timed one, and it calls the pass compiled the one way or
/// the other.
///
/// # Errors
///
/// What the pass itself refuses.
#[cfg(not(target_family = "wasm"))]
fn the_runs_of(
    runs: usize,
    mut pass: impl FnMut(bool) -> Result<f64, String>,
) -> Result<Runs, String> {
    let checksum = pass(true)?;
    let mut times = Vec::with_capacity(runs);
    for _ in 0..runs {
        let started = Instant::now();
        pass(false)?;
        times.push(started.elapsed());
    }
    Ok(Runs { times, checksum })
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
    let between = upper.saturating_sub(lower);
    lower.saturating_add(between.checked_div(2).unwrap_or(between))
}

/// The milliseconds of a time, with which a pass of a few milliseconds is
/// printed to three decimals.
fn milliseconds(time: Duration) -> f64 {
    time.as_secs_f64() * 1000.0
}

/// The line of the table that one pass gets: its name, the best, the
/// median and the worst of its timed runs, and its checksum.
#[cfg(not(target_family = "wasm"))]
fn the_line_of(name: &str, runs: &Runs) -> String {
    let last = runs.times.len().saturating_sub(1);
    format!(
        "{name:<32} {best:>8} {median:>8} {worst:>8}   {checksum:.6e}",
        best = format!("{:.3}", milliseconds(sorted_time(&runs.times, 0))),
        median = format!("{:.3}", milliseconds(median_time(&runs.times))),
        worst = format!("{:.3}", milliseconds(sorted_time(&runs.times, last))),
        checksum = runs.checksum,
    )
}

/// The benchmark reads a file of the disc, and wasm has none. This is what
/// `cargo check --target wasm32-unknown-unknown --all-targets` compiles of
/// it, so that the command which checks that nothing of the crate has left
/// wasm behind can check the benchmarks too.
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
    let block = match the_block_of(&arguments.path) {
        Ok(block) => block,
        Err(problem) => {
            eprintln!("{problem}");
            return ExitCode::FAILURE;
        }
    };
    println!(
        "{path}, one block of {num_vars} variants of {num_individuals} individuals, ploidy \
         {ploidy}, {runs} runs, VECLIB_MAXIMUM_THREADS {veclib}, RAYON_NUM_THREADS {rayon}",
        path = arguments.path.display(),
        num_vars = block.majors.len(),
        num_individuals = block.num_individuals,
        ploidy = block.ploidy,
        runs = arguments.runs,
        veclib = said_about_the_variable("VECLIB_MAXIMUM_THREADS"),
        rayon = said_about_the_variable("RAYON_NUM_THREADS"),
    );
    let mut allele_counts: AlleleCounts = [0; 128];
    let mut codes = vec![0_u8; block.num_individuals];
    let mut counts = [0_u32; 255];
    let mut values = [0.0_f64; 256];
    let mut row = vec![0.0_f64; block.num_individuals];
    let mut scratch = Scratch::of(block.num_individuals);
    let mut lines = Vec::with_capacity(5);
    let mut parts = 0.0_f64;

    let counted = match the_runs_of(arguments.runs, |checksum| {
        if checksum {
            the_pass_of_count_alleles::<true>(&block, &mut allele_counts)
        } else {
            the_pass_of_count_alleles::<false>(&block, &mut allele_counts)
        }
    }) {
        Ok(runs) => runs,
        Err(problem) => {
            eprintln!("{problem}");
            return ExitCode::FAILURE;
        }
    };
    parts += milliseconds(sorted_time(&counted.times, 0));
    lines.push(the_line_of("the counts of the alleles", &counted));

    let coded = match the_runs_of(arguments.runs, |checksum| {
        Ok(if checksum {
            the_pass_of_the_codes::<true>(&block, &mut codes)
        } else {
            the_pass_of_the_codes::<false>(&block, &mut codes)
        })
    }) {
        Ok(runs) => runs,
        Err(problem) => {
            eprintln!("{problem}");
            return ExitCode::FAILURE;
        }
    };
    parts += milliseconds(sorted_time(&coded.times, 0));
    lines.push(the_line_of("the codes of the genotypes", &coded));

    let tallied = match the_runs_of(arguments.runs, |checksum| {
        Ok(if checksum {
            the_pass_of_the_counts::<true>(&block, &mut counts)
        } else {
            the_pass_of_the_counts::<false>(&block, &mut counts)
        })
    }) {
        Ok(runs) => runs,
        Err(problem) => {
            eprintln!("{problem}");
            return ExitCode::FAILURE;
        }
    };
    parts += milliseconds(sorted_time(&tallied.times, 0));
    lines.push(the_line_of("the counts of the codes", &tallied));

    let looked_up = match the_runs_of(arguments.runs, |checksum| {
        Ok(if checksum {
            the_pass_of_the_values::<true>(&block, &mut counts, &mut values, &mut row)
        } else {
            the_pass_of_the_values::<false>(&block, &mut counts, &mut values, &mut row)
        })
    }) {
        Ok(runs) => runs,
        Err(problem) => {
            eprintln!("{problem}");
            return ExitCode::FAILURE;
        }
    };
    parts += milliseconds(sorted_time(&looked_up.times, 0));
    lines.push(the_line_of("the values of the codes", &looked_up));

    let whole = match the_runs_of(arguments.runs, |checksum| {
        if checksum {
            the_pass_of_the_standardized_row::<true>(&block, &mut scratch, &mut row)
        } else {
            the_pass_of_the_standardized_row::<false>(&block, &mut scratch, &mut row)
        }
    }) {
        Ok(runs) => runs,
        Err(problem) => {
            eprintln!("{problem}");
            return ExitCode::FAILURE;
        }
    };
    let whole_best = milliseconds(sorted_time(&whole.times, 0));
    lines.push(the_line_of("a whole standardized row", &whole));

    println!(
        "{name:<32} {best:>8} {median:>8} {worst:>8}   checksum",
        name = "pass over the rows of the block",
        best = "best",
        median = "median",
        worst = "worst",
    );
    println!(
        "{name:<32} {best:>8} {median:>8} {worst:>8}",
        name = "",
        best = "ms",
        median = "ms",
        worst = "ms",
    );
    for line in &lines {
        println!("{line}");
    }
    let share = if whole_best > 0.0 {
        parts / whole_best * 100.0
    } else {
        f64::NAN
    };
    println!(
        "the four passes come to {parts:.3} ms and a whole row takes {whole_best:.3} ms, \
         {share:.1} per 100 of it"
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
