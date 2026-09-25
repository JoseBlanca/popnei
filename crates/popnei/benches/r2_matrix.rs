//! How long the matrix of r² of every pair of a set of variants takes,
//! with the reading of the file timed apart and taken out.
//!
//! r² is the squared correlation between the dosages of two variants, and
//! `calc_r2_matrix` of `docs/specs/ld.md` gives it for every pair of the
//! variants a reader hands it, as a square matrix of `f64`. The number to
//! reach is the 0.50 s of "Speed" of that spec, for 5000 variants of 1000
//! individuals: numpy on Accelerate takes 0.455 s for the same six
//! products and the arithmetic over them, and the tenth above it is what
//! popnei is allowed for building the three matrices of the dosages and
//! for cutting them into tiles.
//!
//! What a run times is two passes over the same file, one after the other.
//! The first builds the reader and calls the matrix, which reads the file
//! itself. The second builds the same reader and reads it to its end with
//! the same fields asked for, the genotypes, the chromosome and the
//! position, adding up the alleles of every block and computing nothing
//! else. The difference of the two is the products and the arithmetic over
//! them, which is what the target is about, and the run prints all three.
//!
//! The memory of the matrix is inside the first of those two times and not
//! outside it: `calc_r2_matrix` asks the machine for it, 200 MB of `f64`
//! at 5000 variants, and gives it back when the run ends. So is the memory
//! of the six sums of each pair of tiles. The printed lines say so, and a
//! number taken from here carries the allocation of the matrix with it.
//!
//! How many variants the matrix is of is `--max-num-vars`, and the pass
//! gives exactly that many: `calc_r2_matrix` refuses a pass that gives
//! more than the number it was asked for, with the error of too many
//! variants, rather than stopping at it, so the benchmark puts a limit in
//! front of it. The limit is `Reblock` over the file, asked for blocks of
//! that many variants, and one block of it given on; `Reblock` is what
//! popnei already puts before a calculation that wants blocks of one size,
//! and over a vars file of 1000 individuals, whose blocks hold 5000
//! variants, a first block of 5000 goes through it with no copy. Both
//! timed passes go through that same limit, so what it costs is in both
//! and the difference is free of it.
//!
//! `--vars-per-tile` is the thing to settle. The matrix is taken tile pair
//! by tile pair, a tile being a set of variants whose three matrices of
//! dosages are built together, because the six sums of a pair of tiles are
//! matrices of the size of that pair and six of the whole set would be six
//! times the result, 1.2 GB at 5000 variants. How many variants a tile
//! holds is a private constant of `ld`, 256, and the comment on it says
//! that a performance review settles it: a reviewer of `docs/plans/ld.md`
//! measured 51.8 ns for a pair of variants at 128, 38.0 ns at 256, 30.1 ns
//! at 512 and 29.1 ns at 1024, none of it on `calc_r2_matrix` itself. The
//! cargo feature `bench-internals` is what makes the tile reachable: it
//! turns on `ld::bench_internals`, which re-exports the private function
//! `calc_r2_matrix` calls with the constant, and nothing of popnei's own
//! builds turns it on.
//!
//! The matrix does not change with the tile, and the checksum is what
//! shows it: the run prints how many variants the matrix is of, how many
//! of its values are not a number, and what the rest of them add to. A
//! pair has no r² when either of its variants has one dosage among its
//! called genotypes, and `docs/specs/ld.md` gives that pair NaN, so the
//! two numbers together read every value of the matrix. They are also what
//! keeps the compiler from dropping the call, and they are taken after the
//! clock stops.
//!
//! The products run on the BLAS of the system, which on this machine is
//! Accelerate and reads `VECLIB_MAXIMUM_THREADS` when the process starts,
//! and the benchmark builds no pool of threads of its own. So one thread
//! is asked for by setting that variable in the environment of the
//! command, and the threads the machine gives by setting neither it nor
//! `RAYON_NUM_THREADS`, which the reader of a vars file does not use and
//! which is printed with each run all the same:
//!
//! ```text
//! VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1 \
//!     cargo bench --bench r2_matrix --features bench-internals -- \
//!     /Users/jose/devel/popnei-bench/big.vars --runs 5 --max-num-vars 5000 \
//!     --vars-per-tile 256
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
//! `--runs` is 5, `--max-num-vars` is 5000 and `--vars-per-tile` is 256
//! when they are not given, which are the cap and the tile the library
//! has. It prints the three times of every run and then the best, the
//! median and the worst of each of the three.
//!
//! How the file is made. The VCF of "Speed" of `docs/specs/filters.md`,
//! 100000 variants of 1000 individuals whose genotypes are missing at a
//! rate of 0.03, is written by `make_big_vcf.py`, beside this file, which
//! says what it writes; it takes 3.1 s and leaves 403 MB. The vars file is
//! what popnei writes of it with the size of block it chooses for 1000
//! individuals:
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
        reason = "the benchmark reads a file of the disc and runs on the BLAS of the system: \
                  in wasm only its empty main is compiled, and what the timing is made of is \
                  left unused"
    )
)]

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{Duration, Instant};

use popnei::block::{Block, BlockReader, Reblock};
use popnei::filters::FilteringStats;
use popnei::io::vars::VarsReader;
use popnei::ld::bench_internals::the_r2_matrix_in_tiles_of;
use popnei::variant::{ChromTable, Needs};

/// How many times the matrix is timed when the command line does not say.
const DEFAULT_RUNS: usize = 5;

/// How many variants the matrix is of when the command line does not say,
/// which is `MAX_NUM_VARS_OF_THE_MATRIX` of `ld`, the number a Python or a
/// TypeScript user gets when they name none.
const DEFAULT_MAX_NUM_VARS: usize = 5000;

/// How many variants a tile of the products holds when the command line
/// does not say. It is taken from the library, so a run with no
/// `--vars-per-tile` measures the tile `calc_r2_matrix` really uses and a
/// change of that constant cannot leave this benchmark reporting the old
/// one under the new name.
const DEFAULT_VARS_PER_TILE: usize = popnei::ld::bench_internals::VARS_OF_A_TILE;

/// The fields the matrix asks its reader for, which the pass that reads
/// and computes nothing asks for too.
const THE_FIELDS_OF_THE_MATRIX: Needs = Needs::GTS.union(Needs::CHROM_POS);

/// What the command line asked for.
struct Arguments {
    path: PathBuf,
    runs: usize,
    max_num_vars: usize,
    vars_per_tile: usize,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
r2_matrix <path to a vars file> [--runs n] [--max-num-vars n]
          [--vars-per-tile n]

It times the matrix of r² of every pair of the first n variants of that
file, and beside it a pass that reads the same variants with the same
fields and computes nothing, so that what the products and the arithmetic
over them take is the difference of the two. The memory of the matrix,
200 MB of f64 at 5000 variants, is asked for inside the first of the two
and is part of it.

  --runs n            how many times both passes are timed, 5 by default
  --max-num-vars n    how many variants the matrix is of, 5000 by default,
                      which is the cap a user gets when they name none
  --vars-per-tile n   how many variants a tile of the products holds, 256
                      by default, which is the number the library uses
  --help              this

The pass gives exactly n variants: the matrix refuses a pass of more
variants than it was asked for instead of stopping at that number, so the
benchmark reads the file in blocks of n variants and gives the first of
them on. One run that is not timed comes first, so that the timed runs pay
neither the page faults of the first touch of the memory a pass works in
nor a read of the disc. Each run prints its two times and their
difference, with how many variants the matrix is of, how many of its
values are not a number and what the rest of them add to, which say that
the matrix was computed and that a change of the tile changed no value;
and then the best, the median and the worst of each of the three times.";

/// What to run, or the message that says what the command line should have
/// been.
fn arguments() -> Result<Arguments, String> {
    let mut path: Option<PathBuf> = None;
    let mut runs = DEFAULT_RUNS;
    let mut max_num_vars = DEFAULT_MAX_NUM_VARS;
    let mut vars_per_tile = DEFAULT_VARS_PER_TILE;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--runs" | "--max-num-vars" | "--vars-per-tile" => {
                let name = arg.clone();
                let number = args
                    .next()
                    .ok_or_else(|| format!("{name} takes a number and none came after it"))?
                    .parse::<usize>()
                    .map_err(|_| format!("{name} takes a number"))?;
                match name.as_str() {
                    "--runs" => runs = number,
                    "--max-num-vars" => max_num_vars = number,
                    _ => vars_per_tile = number,
                }
            }
            // `cargo bench` adds this to the command line of every bench,
            // to tell a harness that has tests too to run its benchmarks.
            // This one has only this benchmark and takes it as nothing.
            "--bench" => {}
            "--help" | "-h" => return Err(USAGE.to_owned()),
            // An argument that is not one of those and looks like one is
            // refused instead of being taken for the path of the file, as
            // the other benchmarks refuse it: a `--vars-per-tile=512` read
            // as a path and dropped leaves a run that timed the tile of the
            // library and calls it the tile of 512.
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
    if runs == 0 || max_num_vars == 0 || vars_per_tile == 0 {
        return Err("--runs, --max-num-vars and --vars-per-tile are 1 or more".to_owned());
    }
    Ok(Arguments {
        path,
        runs,
        max_num_vars,
        vars_per_tile,
    })
}

/// The first `num_vars` variants of a file, as a reader: `Reblock` over
/// the file asked for blocks of that many variants, with its first block
/// given on and no other.
///
/// It is what stands between the file and the matrix, which takes the
/// number of variants it is asked for as a cap and refuses a pass that
/// gives more instead of stopping there. A file of fewer variants than
/// `num_vars` gives them all, in one block, and the matrix is of those.
struct TheFirstVarsOfTheFile {
    reader: Reblock<VarsReader<BufReader<File>>>,
    /// Whether the block was given already. After it there is no block.
    given: bool,
}

impl TheFirstVarsOfTheFile {
    /// The reader over the vars file at `path` that gives its first
    /// `num_vars` variants in one block.
    ///
    /// # Errors
    ///
    /// When the file cannot be opened or its header cannot be read, and
    /// when the genotypes of a block of `num_vars` variants are more than
    /// this machine addresses.
    fn of(path: &Path, num_vars: usize) -> Result<TheFirstVarsOfTheFile, popnei::Error> {
        Ok(TheFirstVarsOfTheFile {
            reader: Reblock::new(VarsReader::from_path(path)?, Some(num_vars))?,
            given: false,
        })
    }
}

impl BlockReader for TheFirstVarsOfTheFile {
    /// The first block of the file, of the variants it was built for, and
    /// `None` after it.
    ///
    /// # Errors
    ///
    /// Whatever the reader of the file fails with.
    fn next_block(&mut self) -> Result<Option<Block>, popnei::Error> {
        if self.given {
            return Ok(None);
        }
        self.given = true;
        self.reader.next_block()
    }

    fn individuals(&self) -> &[String] {
        self.reader.individuals()
    }

    fn ploidy(&self) -> usize {
        self.reader.ploidy()
    }

    fn chroms(&self) -> &ChromTable {
        self.reader.chroms()
    }

    fn set_needs(&mut self, needs: Needs) {
        self.reader.set_needs(needs);
    }

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        self.reader.filtering_stats()
    }
}

/// One run: how long the matrix took, how long reading the same variants
/// took with nothing computed, and the line that says what the matrix
/// holds.
struct Run {
    matrix: Duration,
    reading: Duration,
    did: String,
}

impl Run {
    /// The matrix less the reading, which is the products and the
    /// arithmetic over them.
    ///
    /// The two passes read the same variants through the same readers, so
    /// the reading is the shorter of the two and the difference is a time;
    /// a machine that was doing something else during one of them and not
    /// the other could turn it around, and 0 is what that gives.
    fn over_the_reading(&self) -> Duration {
        self.matrix.saturating_sub(self.reading)
    }
}

/// The matrix of the first `max_num_vars` variants of the file at `path`,
/// with the products in tiles of `vars_per_tile` variants, timed; and then
/// the same variants read through the same readers with the same fields
/// and nothing computed, timed too.
///
/// Both times hold the building of their reader and the reading of the
/// file, so the difference of the two is free of both. The memory of the
/// matrix is inside the first of them: the call asks the machine for it.
///
/// The alleles of every block are added up in the pass that computes
/// nothing, so that a reader which filled the genotypes only when somebody
/// read them could not make that pass short and the difference large with
/// nothing to show it; how many there were is printed. The checksum of the
/// matrix, its values that are not a number and what the rest of them add
/// to, is taken after the clock stops.
fn one_run(path: &Path, max_num_vars: usize, vars_per_tile: usize) -> Result<Run, popnei::Error> {
    let started = Instant::now();
    let mut reader = TheFirstVarsOfTheFile::of(path, max_num_vars)?;
    let r2_matrix = the_r2_matrix_in_tiles_of(&mut reader, max_num_vars, vars_per_tile)?;
    let matrix = started.elapsed();
    let started = Instant::now();
    let mut reader = TheFirstVarsOfTheFile::of(path, max_num_vars)?;
    reader.set_needs(THE_FIELDS_OF_THE_MATRIX);
    let mut alleles: u64 = 0;
    while let Some(block) = reader.next_block()? {
        // A block of more alleles than a u64 counts is one this machine
        // does not address.
        alleles = alleles.saturating_add(u64::try_from(block.gts.len()).unwrap_or(u64::MAX));
    }
    let reading = started.elapsed();
    let mut not_a_number: u64 = 0;
    let mut adds_to = 0.0_f64;
    for value in r2_matrix.r2() {
        if value.is_nan() {
            not_a_number = not_a_number.saturating_add(1);
            continue;
        }
        adds_to += *value;
    }
    let did = format!(
        "{num_vars} variants, the {values} values of the matrix add to {adds_to:.6} over the \
         {not_a_number} of them that are not a number, and the pass that computed nothing read \
         {alleles} alleles{phases}",
        num_vars = r2_matrix.num_vars(),
        values = r2_matrix.r2().len(),
        phases = the_phases_of_the_pass(),
    );
    Ok(Run {
        matrix,
        reading,
        did,
    })
}

/// The two clocks of the phases of the pass, as a piece of the line of a
/// run: how long the pass was inside `next_block` of its reader and how
/// long it was working on the blocks the reader gave. Taking them zeroes
/// them, so each run prints its own.
///
/// It is the cargo feature `bench-phases` of the core crate, and without it
/// there is nothing to print: the pass then calls no clock at all. `cargo
/// bench --features bench-phases --bench r2_matrix` is what turns it on. The
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

/// The seconds of a time, with the three decimals that a matrix of a few
/// tenths of a second is worth reporting to.
fn seconds(time: Duration) -> String {
    format!("{:.3} s", time.as_secs_f64())
}

/// The line that says what `what` came to over the runs: its best, its
/// median and its worst.
fn the_spread_of(what: &str, times: &[Duration]) -> String {
    let last = times.len().saturating_sub(1);
    format!(
        "{what}: best {best}, median {median}, worst {worst}",
        best = seconds(sorted_time(times, 0)),
        median = seconds(median_time(times)),
        worst = seconds(sorted_time(times, last)),
    )
}

/// The benchmark reads a file of the disc and runs on the BLAS of the
/// system, and wasm has neither. This is what `cargo check --target
/// wasm32-unknown-unknown --all-targets` compiles of it, so that the
/// command which checks that nothing of the crate has left wasm behind can
/// check the benchmarks too.
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
        "{path}, {runs} runs, the matrix of {max_num_vars} variants in tiles of \
         {vars_per_tile}, VECLIB_MAXIMUM_THREADS {veclib}, RAYON_NUM_THREADS {rayon}",
        path = arguments.path.display(),
        runs = arguments.runs,
        max_num_vars = arguments.max_num_vars,
        vars_per_tile = arguments.vars_per_tile,
        veclib = said_about_the_variable("VECLIB_MAXIMUM_THREADS"),
        rayon = said_about_the_variable("RAYON_NUM_THREADS"),
    );
    println!(
        "the memory of the matrix, 8 bytes for each pair of variants, is asked for inside the \
         time of the matrix and not before it"
    );
    let mut of_the_matrix = Vec::with_capacity(arguments.runs);
    let mut of_the_reading = Vec::with_capacity(arguments.runs);
    let mut of_the_difference = Vec::with_capacity(arguments.runs);
    for run in 0..=arguments.runs {
        let done = match one_run(
            &arguments.path,
            arguments.max_num_vars,
            arguments.vars_per_tile,
        ) {
            Ok(done) => done,
            Err(error) => {
                eprintln!("{path}: {error}", path = arguments.path.display());
                return ExitCode::FAILURE;
            }
        };
        if run == 0 {
            println!(
                "the first run, which is not timed: {did}, the matrix in {matrix}",
                did = done.did,
                matrix = seconds(done.matrix),
            );
            continue;
        }
        println!(
            "run {run}: the matrix {matrix}, reading {reading}, the matrix less the reading \
             {difference}, {did}",
            matrix = seconds(done.matrix),
            reading = seconds(done.reading),
            difference = seconds(done.over_the_reading()),
            did = done.did,
        );
        of_the_matrix.push(done.matrix);
        of_the_reading.push(done.reading);
        of_the_difference.push(done.over_the_reading());
    }
    println!("{}", the_spread_of("the matrix", &of_the_matrix));
    println!("{}", the_spread_of("reading", &of_the_reading));
    println!(
        "{}",
        the_spread_of("the matrix less the reading", &of_the_difference)
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
