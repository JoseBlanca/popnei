//! How long the Kosman distance of every pair of individuals takes, the
//! `calc_kosman_sums` of `docs/specs/dists.md`, over a source whose blocks
//! are already in memory and over a vars file.
//!
//! "Speed" of that spec gives three numbers to reach on 100000 variants of
//! 1000 diploid individuals, biallelic, 3 in 100 genotypes missing: 0.97 s
//! on one thread, 0.38 s on 18 cores and 1.43 s in wasm, each of them of
//! the calculation with the time of reading the file alone taken out, so
//! that the reader is not in them.
//! `docs/reports/dists-kosman-measurement.md` measured 1.154 s and 0.625 s
//! against the first two, by taking that subtraction from Python.
//!
//! This benchmark takes the reader out in a second way, which needs no
//! subtraction: it times the same calculation over a reader that hands out
//! blocks it already holds in memory. The two settings it times, each of
//! them on the pool of threads of the process and then inside a pool of one
//! thread:
//!
//! - **the blocks in memory**: a reader written here gives 20 blocks of
//!   5000 variants of 1000 diploid individuals, biallelic, 3 in 100
//!   genotypes missing whole, drawn by a generator with a fixed seed, so
//!   that the blocks are the same on every run and on every machine. The
//!   blocks are built before the clock starts and handed out one at a time,
//!   so nothing of a reader, no file, no decompression and no copy, is
//!   inside the clock. This is the number to hold against the 0.97 s and
//!   the 0.38 s.
//! - **the vars file**: the same calculation over a `VarsReader` on the
//!   bytes of the file, which are read from the disc before the clock
//!   starts. What is inside the clock is what a user waits for, the reader
//!   and the calculation together, and its difference from the setting
//!   above is what the reader of a vars file costs.
//!
//! The genotypes of the two settings are not the same genotypes: the file
//! holds what `make_big_vcf.py` simulated and the blocks in memory what the
//! generator here drew. Neither the sets of bits nor the count of a pair
//! branches on the value of a genotype, so what the two settings differ in
//! is the reader and not the work of the calculation.
//!
//! What it checks. The call has to succeed and to give the pairs of 1000
//! individuals, 499500 of them, which is counted after the clock stops and
//! printed: a calculation that gave the pairs of another number of
//! individuals, or that was given fewer blocks than the setting says, would
//! otherwise be a time that looks reasonable and is of other work. Nothing
//! else about the numbers is asserted here; `dists.rs` has the tests that
//! assert the distances.
//!
//! It is run with cargo, which passes what comes after the two dashes to
//! it, and the release build is what `cargo bench` makes:
//!
//! ```text
//! cargo bench --bench kosman_dists -- /Users/jose/devel/popnei-bench/big.vars --runs 5
//! ```
//!
//! `--runs` is 5 when it is not given. Before the timed runs of each of the
//! four there is one run that is not timed, whose time is printed beside
//! them: the first touch of the memory a calculation works in costs page
//! faults that a process pays once, and they would otherwise all fall on
//! whichever of the four ran first.
//!
//! The load average of the machine is not read here: a timing that is
//! reported is taken with nothing else running, and `sysctl -n vm.loadavg`
//! before and after the invocation is what says the machine was quiet, as
//! `docs/reports/dists-kosman-measurement.md` reports it.
//!
//! How the file is made: `docs/reports/dists-kosman-measurement.md` has the
//! two commands, `make_big_vcf.py` beside this file for the VCF of 100000
//! variants and popnei's `write_vars` for the vars file of it.

#![cfg_attr(
    target_family = "wasm",
    allow(
        dead_code,
        unused_imports,
        reason = "the benchmark reads a file of the disc and builds a pool of threads, and \
                  wasm has neither: in wasm only its empty main is compiled, and what the \
                  timing is made of is left unused"
    )
)]

use std::io::Cursor;
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::{Duration, Instant};

use popnei::block::{Block, BlockReader};
use popnei::dists::calc_kosman_sums;
use popnei::filters::FilteringStats;
use popnei::io::vars::VarsReader;
use popnei::variant::{ChromTable, MISSING_ALLELE, Needs};

/// How many times each of the four is timed when the command line does not
/// say.
const DEFAULT_RUNS: usize = 5;

/// How many individuals the blocks built here hold, which is the panel of
/// "Speed" of `docs/specs/dists.md`.
const NUM_INDIVIDUALS: usize = 1000;

/// How many alleles one genotype of those blocks holds: they are diploid.
const PLOIDY: usize = 2;

/// How many variants one of those blocks holds, the size popnei chooses for
/// 1000 individuals and the size of the batches of `big.vars`.
const NUM_VARS_PER_BLOCK: usize = 5000;

/// How many of those blocks the reader in memory hands out, which makes the
/// 100000 variants of "Speed".
const NUM_BLOCKS: usize = 20;

/// The seed of the generator that draws the genotypes of the block in
/// memory. Any value does: what it is for is that the blocks are the same
/// on every run and on every machine.
const SEED: u64 = 42;

/// Out of 1024 genotypes, how many are missing whole: 31 in 1024 is 3.03 in
/// 100, the 3 in 100 of the dataset of "Speed".
const MISSING_IN_1024: u64 = 31;

/// What the command line asked for.
struct Arguments {
    path: PathBuf,
    runs: usize,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
kosman_dists <path to a vars file> [--runs n]

It times the Kosman distance of every pair of individuals, `calc_kosman_sums`,
in two settings, each of them on the pool of threads of the process and then
inside a pool of one thread: over 20 blocks of 5000 variants of 1000 diploid
individuals that a generator with a fixed seed drew and that are built before
the clock starts, so that no reader is inside it; and over the vars file
given, whose bytes are read from the disc before the clock starts. The file
of the first setting is the dataset of \"Speed\" of docs/specs/dists.md,
100000 variants of 1000 individuals.

  --runs n   how many times each of the four is timed, 5 by default
  --help     this

It prints the wall time of each run and then the best, the median and the
worst of them, with the pairs the call gave, 499500, and the variants it
read. The best of the blocks in memory is what the 0.97 s on one thread and
the 0.38 s on 18 cores of that spec are to be held against.";

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
                    .map_err(|_| "--runs takes a number".to_owned())?;
            }
            // `cargo bench` adds this to the command line of every bench,
            // to tell a harness that has tests too to run its benchmarks.
            // This one has only this benchmark and takes it as nothing.
            "--bench" => {}
            "--help" | "-h" => return Err(USAGE.to_owned()),
            // An argument that is not one of those and looks like one is
            // refused instead of being taken for the path of the file, as
            // the other benchmarks refuse it: a `--runs=10` read as a path
            // and dropped leaves a run that timed something else and says
            // nothing.
            other if other.starts_with('-') => {
                return Err(format!(
                    "`{other}` is not an argument of this benchmark\n\n{USAGE}"
                ));
            }
            other => path = Some(PathBuf::from(other)),
        }
    }
    let Some(path) = path else {
        return Err(format!("no vars file was given\n\n{USAGE}"));
    };
    if runs == 0 {
        return Err("--runs is 1 or more".to_owned());
    }
    Ok(Arguments { path, runs })
}

/// A generator of numbers that look random, seeded, so that the genotypes
/// it draws are the same on every run and on every machine.
///
/// It is the linear congruential generator of Knuth's table, with the
/// multiplier and the increment that `rand`'s `Pcg64` and the `drand48` of
/// glibc use: the state is multiplied and added in 64 bits, wrapping, and
/// the state is what is given. What it is for is a block of genotypes that
/// does not change from one run to the next, and not the quality of a
/// simulation: the genotypes of a real dataset come from
/// `make_big_vcf.py`, and neither the sets of bits nor the count of a pair
/// branches on the value of a genotype.
struct Numbers {
    state: u64,
}

impl Numbers {
    /// The generator at `seed`.
    fn seeded(seed: u64) -> Numbers {
        Numbers { state: seed }
    }

    /// The next number of the generator, any of the values a `u64` holds.
    fn next(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    /// The next number of the generator from 0 to 1023, which is how the
    /// draws below are made: a mask and a comparison, with no division.
    fn next_in_1024(&mut self) -> u64 {
        self.next() & 1023
    }
}

/// One block of 5000 variants of 1000 diploid biallelic individuals, with 3
/// in 100 of the genotypes missing whole, drawn by a generator at `SEED`.
///
/// Each variant gets a frequency of the alternative allele drawn from 102
/// to 869 in 1024, so that no variant is of one allele alone, and each
/// allele of each genotype is drawn at that frequency. A genotype that is
/// missing has both of its alleles at [`MISSING_ALLELE`], which is what a
/// reader of popnei gives for one that was not called: a half called
/// genotype is missing too, and the sets of bits leave out the whole
/// genotype either way.
fn the_block() -> Block {
    let mut numbers = Numbers::seeded(SEED);
    let mut gts: Vec<i8> = Vec::new();
    // The variants times the individuals times the ploidy, which the block
    // below says it holds.
    gts.reserve_exact(
        NUM_VARS_PER_BLOCK
            .saturating_mul(NUM_INDIVIDUALS)
            .saturating_mul(PLOIDY),
    );
    for _ in 0..NUM_VARS_PER_BLOCK {
        // From 102 to 869 in 1024, which is 0.0996 to 0.849.
        let frequency = 102_u64.saturating_add(numbers.next_in_1024() & 767);
        for _ in 0..NUM_INDIVIDUALS {
            if numbers.next_in_1024() < MISSING_IN_1024 {
                gts.extend(std::iter::repeat_n(MISSING_ALLELE, PLOIDY));
                continue;
            }
            for _ in 0..PLOIDY {
                gts.push(i8::from(numbers.next_in_1024() < frequency));
            }
        }
    }
    Block {
        num_vars: NUM_VARS_PER_BLOCK,
        num_individuals: NUM_INDIVIDUALS,
        ploidy: PLOIDY,
        gts,
        chrom: None,
        pos: None,
        id: None,
        alleles: None,
        qual: None,
    }
}

/// A source of blocks that are already in memory, which hands out the ones
/// it was given one at a time and reads nothing.
///
/// What it is for is a timing of the calculation with no reader in it: the
/// blocks are built before the clock starts, and `next_block` moves one out
/// of the vector.
struct BlocksInMemory {
    blocks: std::vec::IntoIter<Block>,
    individuals: Vec<String>,
    chroms: ChromTable,
}

impl BlocksInMemory {
    /// `NUM_BLOCKS` copies of `block`, with names for the individuals of
    /// the block.
    fn of(block: &Block) -> BlocksInMemory {
        let blocks: Vec<Block> = (0..NUM_BLOCKS).map(|_| copy_of(block)).collect();
        let individuals = (0..NUM_INDIVIDUALS)
            .map(|number| format!("individual_{number}"))
            .collect();
        BlocksInMemory {
            blocks: blocks.into_iter(),
            individuals,
            chroms: ChromTable::new(),
        }
    }
}

impl BlockReader for BlocksInMemory {
    fn next_block(&mut self) -> popnei::Result<Option<Block>> {
        Ok(self.blocks.next())
    }

    fn individuals(&self) -> &[String] {
        &self.individuals
    }

    fn ploidy(&self) -> usize {
        PLOIDY
    }

    fn chroms(&self) -> &ChromTable {
        &self.chroms
    }

    fn set_needs(&mut self, _needs: Needs) {
        // The blocks hold the genotypes and nothing else, whatever is asked
        // for.
    }

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        Vec::new()
    }
}

/// A block with the same genotypes as `block`, which is how the copies the
/// reader in memory hands out are made.
///
/// `Block` does not derive `Clone`, and the blocks built here hold the
/// genotypes and none of the five other columns, so this copies the
/// genotypes and says what the block is.
fn copy_of(block: &Block) -> Block {
    Block {
        num_vars: block.num_vars,
        num_individuals: block.num_individuals,
        ploidy: block.ploidy,
        gts: block.gts.clone(),
        chrom: None,
        pos: None,
        id: None,
        alleles: None,
        qual: None,
    }
}

/// One run of one of the four: how long it took and the line that says what
/// it did, which holds the pairs and the variants that say the call did the
/// work the setting names.
struct Run {
    took: Duration,
    did: String,
}

/// One calculation over `reader`, timed from the first block to the sums of
/// every pair.
///
/// The pairs are counted after the clock stops, which is also what keeps
/// the sums from being dropped as unused.
fn one_calculation(reader: &mut impl BlockReader) -> Result<Run, String> {
    let started = Instant::now();
    let sums = calc_kosman_sums(reader).map_err(|problem| problem.to_string())?;
    let took = started.elapsed();
    let num_pairs = sums.dists(0).count();
    let expected = NUM_INDIVIDUALS
        .saturating_sub(1)
        .saturating_mul(NUM_INDIVIDUALS)
        .checked_div(2)
        .unwrap_or_default();
    if num_pairs != expected {
        return Err(format!(
            "the call gave {num_pairs} pairs where the {NUM_INDIVIDUALS} individuals of this \
             benchmark make {expected}"
        ));
    }
    Ok(Run {
        took,
        did: format!(
            "{num_pairs} pairs of {individuals} individuals over {vars} variants",
            individuals = sums.num_individuals(),
            vars = sums.num_vars(),
        ),
    })
}

/// The time at the place `part` of the times sorted from the shortest to
/// the longest.
fn sorted_time(times: &[Duration], part: usize) -> Duration {
    let mut sorted = times.to_vec();
    sorted.sort_unstable();
    sorted.get(part).copied().unwrap_or_default()
}

/// The median of the times: the middle one when they are an odd number, and
/// the middle of the two middle ones when they are an even number, which
/// `sorted_time` at half of them is not.
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

/// The seconds of a time, with the three decimals that a calculation of a
/// second or so is worth reporting to.
fn seconds(time: Duration) -> String {
    format!("{:.3} s", time.as_secs_f64())
}

/// It runs one of the four once without timing it and then `runs` times,
/// printing what each run took and then the best, the median and the worst
/// of the timed ones.
///
/// The run that is not timed is there for the page faults of the first
/// touch of the memory a calculation works in, which a process pays once
/// and which would otherwise all fall on whichever of the four ran first.
fn time_it(
    what: &str,
    runs: usize,
    mut run: impl FnMut() -> Result<Run, String>,
) -> Result<(), String> {
    let mut times = Vec::with_capacity(runs);
    for number in 0..=runs {
        let done = run()?;
        if number == 0 {
            println!(
                "{what}, the first run, which is not timed: {took}, {did}",
                took = seconds(done.took),
                did = done.did,
            );
            continue;
        }
        println!(
            "{what}, run {number}: {took}, {did}",
            took = seconds(done.took),
            did = done.did,
        );
        times.push(done.took);
    }
    let last = times.len().saturating_sub(1);
    println!(
        "{what}: best {best}, median {median}, worst {worst}",
        best = seconds(sorted_time(&times, 0)),
        median = seconds(median_time(&times)),
        worst = seconds(sorted_time(&times, last)),
    );
    Ok(())
}

/// The benchmark reads a file of the disc and builds a pool of threads, and
/// wasm has neither; rayon is not a dependency of the wasm targets either.
/// This is what `cargo check --target wasm32-unknown-unknown --all-targets`
/// compiles of it, so that the command which checks that nothing of the
/// crate has left wasm behind can check the benchmarks too.
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
    let bytes = match std::fs::read(&arguments.path) {
        Ok(bytes) => bytes,
        Err(problem) => {
            eprintln!("{path}: {problem}", path = arguments.path.display());
            return ExitCode::FAILURE;
        }
    };
    let on_one_thread = match rayon::ThreadPoolBuilder::new().num_threads(1).build() {
        Ok(pool) => pool,
        Err(error) => {
            eprintln!("the pool of one thread: {error}");
            return ExitCode::FAILURE;
        }
    };
    let of_the_process = rayon::current_num_threads();
    let block = the_block();
    println!(
        "{blocks} blocks of {vars} variants of {individuals} individuals in memory, and \
         {path}, {bytes} bytes; {runs} runs of each, on the {of_the_process} threads of the \
         process and on one thread",
        blocks = NUM_BLOCKS,
        vars = NUM_VARS_PER_BLOCK,
        individuals = NUM_INDIVIDUALS,
        path = arguments.path.display(),
        bytes = bytes.len(),
        runs = arguments.runs,
    );

    let in_memory = || {
        let mut reader = BlocksInMemory::of(&block);
        one_calculation(&mut reader)
    };
    let of_the_file = || {
        let mut reader =
            VarsReader::new(Cursor::new(&bytes)).map_err(|problem| problem.to_string())?;
        one_calculation(&mut reader)
    };
    let timed = time_it(
        &format!("the blocks in memory, {of_the_process} threads"),
        arguments.runs,
        in_memory,
    )
    .and_then(|()| {
        time_it("the blocks in memory, 1 thread", arguments.runs, || {
            on_one_thread.install(in_memory)
        })
    })
    .and_then(|()| {
        time_it(
            &format!("the vars file, {of_the_process} threads"),
            arguments.runs,
            of_the_file,
        )
    })
    .and_then(|()| {
        time_it("the vars file, 1 thread", arguments.runs, || {
            on_one_thread.install(of_the_file)
        })
    });
    if let Err(problem) = timed {
        eprintln!("{path}: {problem}", path = arguments.path.display());
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
