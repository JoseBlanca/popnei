//! How long popnei takes to read the genotypes of a vars file and to write
//! one, on the panel that "The compression" of `docs/specs/io_vars.md` was
//! measured on.
//!
//! It times three things, each of them over a file that is already in
//! memory, so that no read of the disc is inside the clock:
//!
//! - a pass with the genotypes alone asked for, which is what a calculation
//!   over the genotypes asks for and what the 21 ms of "Speed" of that spec
//!   are the number to reach for: no more than a tenth above the 19.1 to
//!   19.3 ms that arrow-rs alone took there to decompress the four batches
//!   of that panel and copy the genotypes of each into a vector;
//! - the same pass with every field asked for, the chromosome, the
//!   position, the id, the alleles and the quality beside the genotypes,
//!   which is what a user who writes the file again asks for;
//! - the writing of the blocks of that file into a `Vec<u8>`.
//!
//! A pass reads every byte of the file whatever is asked for, and
//! decompresses only the columns that were asked for, so the two passes
//! differ in the decompression and in the vectors the other five columns
//! are built into.
//!
//! The file is not the one the spec took its 19.1 to 19.3 ms on. That one
//! held the `gts` column alone and pyarrow 23.0.0 wrote it, 15.68 MB; this
//! one holds the six columns of a vars file, popnei wrote it, so `lz4_flex`
//! compressed its buffers, and it is 16280410 bytes for the panel of 20000
//! variants. The genotypes are of the same simulation, and the 3 in 100
//! that are missing are drawn at another point of the generator than
//! pyNei's script draws them at, so the compressed bytes are not the same
//! either. What stands against the 21 ms is still the pass with the
//! genotypes alone: it decompresses that one column, and reads past the
//! buffers of the other five without decompressing them.
//!
//! Before any of the three is timed there is one pass with every field
//! asked for whose time is not taken. The first touch of the memory a pass
//! works in costs page faults that a process pays once, and they all fell
//! on whichever section ran first, which is the one held against the 21 ms:
//! over 13 invocations the first run of the genotypes alone took 22.06 to
//! 30.63 ms where the runs after it took 20.16 to 21.9 ms, and the two
//! sections that come later in the same process show no such run. That
//! first pass is also what says which columns the file holds, and every
//! block of the two timed passes is then checked to hold the fields that
//! were asked for and no other: the two passes differ by less than a
//! millisecond, so a projection that decompressed every column would look
//! like noise and nothing else would catch it.
//!
//! What is inside the clock. Of a pass: building the `VarsReader` over a
//! `Cursor` of the bytes, which reads the schema and the footer, every
//! `next_block` to the end of the file, and the sum of the genotypes of
//! each block, which is there so that the work cannot be dropped as unused
//! and which the run prints. Of the write: `write_block` of each block and
//! `finish`. The blocks the write is given go in by value, so each run
//! reads them from the bytes in memory before its clock starts, and the
//! writer is built before it too; the `Vec<u8>` the file is written into
//! starts empty, so its growth is inside the clock, as it is for a caller
//! who writes a vars file into memory.
//!
//! The reader and the writer of a vars file run on the thread that calls
//! them: neither uses rayon, so this benchmark has no `--threads`.
//!
//! How the file is made. The VCF of 1000 individuals and 20000 variants of
//! that panel, and the vars file popnei writes from it, with the size of
//! block popnei chooses for 1000 individuals, 5000 variants, which makes
//! four batches:
//!
//! ```text
//! uv run --no-project --with numpy python \
//!     crates/popnei/benches/make_big_vcf.py /tmp/panel.vcf 20000
//! uv run maturin develop && uv run python -c "import popnei; \
//!     popnei.write_vars(popnei.open_vcf('/tmp/panel.vcf'), '/tmp/panel.vars')"
//! ```
//!
//! `make_big_vcf.py`, beside this file, says what it writes; the 20000 is
//! its second argument, and without it the VCF is the one of 100000
//! variants that `read_vcf.rs` is run on. `maturin develop` is there
//! because an `uv sync` takes the module out of the environment, as
//! `pyproject.toml` says, and `import popnei` then fails with a
//! `ModuleNotFoundError`; the debug build it makes is enough, because what
//! is timed here is the core crate, which `cargo bench` builds with the
//! bench profile, and not Python.
//!
//! It is run with cargo, which passes what comes after the two dashes to
//! it:
//!
//! ```text
//! cargo bench --bench vars_file -- /tmp/panel.vars --runs 5
//! ```
//!
//! cargo runs it with `crates/popnei` as its working directory, so a
//! relative path is read from there and an absolute one is the plainer
//! thing to give.
//!
//! `--runs` is 5 when it is not given, which is how many runs "The
//! compression" took its numbers over. It prints the wall time of each run
//! of each of the three, and then their best, median and worst; with an
//! even number of runs the median is the middle of the two middle times.
//! The best is what is compared with the 21 ms, because that is how the
//! spec took its number, and the median beside it says how much the machine
//! was doing something else.
//!
//! What the number depends on. On a quiet machine, a load average of 1.77,
//! the best of 5 runs of the genotypes alone was 20.16 to 20.36 ms over six
//! invocations; with two compilers running beside it, 21.4 to 25.0 ms; and
//! under `taskpolicy -b`, which puts the process on the efficiency cores,
//! 66.48 ms for the same 4.03e9 instructions. So a timing that is reported
//! is taken with nothing else running, and `/usr/bin/time -l` on the binary
//! gives the instructions retired and the cycles elapsed beside the wall
//! time: a run whose cycles divided by its wall time are well under 3 GHz
//! was not on a performance core, and it says nothing about the 21 ms.

#![cfg_attr(
    target_family = "wasm",
    allow(
        dead_code,
        unused_imports,
        reason = "the benchmark is native: in wasm only its empty main is compiled, and \
                  what the timing is made of is left unused"
    )
)]

use std::io::Cursor;
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::{Duration, Instant};

use popnei::block::BlockReader;
use popnei::io::vars::{VarsReader, VarsWriter};
use popnei::variant::Needs;

/// How many times each of the three is timed when the command line does not
/// say.
const DEFAULT_RUNS: usize = 5;

/// What the command line asked for.
struct Arguments {
    path: PathBuf,
    runs: usize,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
vars_file <path to a vars file> [--runs n]

It times three things on that file, each of them with the file already in
memory: a pass over it with the genotypes alone asked for, a pass with
every field asked for, and the writing of its blocks, read before the
clock starts, into a vector of bytes. Every pass goes through
`next_block` to the end of the file and adds up the genotypes of each
block, which it prints, so that the work is not dropped as unused. One
pass that is not timed comes before the three, because the first touch
of the memory a pass works in costs page faults that a process pays
once, and it fails when a block holds a column that nobody asked for.

  --runs n   how many times each of the three is timed, 5 by default
  --help     this

It prints the wall time of each run and then the best, the median and the
worst of them. The best of the pass with the genotypes alone is what
`docs/specs/io_vars.md` asks for 21 ms of, on a file of 1000 individuals
and 20000 variants; the header of this file says how that file is made.";

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
            // `read_vcf.rs` refuses it: a `--runs=10` read as a path and
            // dropped leaves a run that timed something else and says
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

/// One run of one of the three: how long it took and the line that says
/// what it did, which holds the numbers that keep the work from being
/// dropped as unused.
struct Run {
    took: Duration,
    did: String,
}

/// The genotypes of a block added up, which is the work a pass does with
/// them.
///
/// `wrapping_add` is the addition here because the sum is a checksum and
/// not a result: the genotypes of a block are -1, 0 and the alternative
/// alleles, so a dataset of more than 9.2e18 alleles would be needed to
/// reach what an `i64` holds, and a checked addition of every genotype
/// would be measuring the check.
fn genotypes_added_up(gts: &[i8], from: i64) -> i64 {
    gts.iter()
        .fold(from, |sum, allele| sum.wrapping_add(i64::from(*allele)))
}

/// That a block holds the fields that were asked for and no others, which
/// is what says that the batch was decompressed with the projection of
/// those fields and not of every column.
///
/// Nothing else catches a projection that decompressed what nobody asked
/// for: the pass with the genotypes alone and the pass with every field
/// differ by less than a millisecond, which is inside what the machine
/// moves a timing by.
fn the_fields_are(fields: Needs, asked_for: Needs) -> Result<(), String> {
    if fields == asked_for {
        return Ok(());
    }
    let decompressed = fields.difference(asked_for);
    let missing = asked_for.difference(fields);
    if missing.is_empty() {
        return Err(format!(
            "a block holds {fields} where {asked_for} was asked for: nothing asked for \
             {decompressed} and it was decompressed"
        ));
    }
    if decompressed.is_empty() {
        return Err(format!(
            "a block holds {fields} where {asked_for} was asked for: {missing} is not there"
        ));
    }
    Err(format!(
        "a block holds {fields} where {asked_for} was asked for"
    ))
}

/// One pass with every field asked for, before anything is timed and with
/// no clock on it: which fields the blocks of the file hold, and the line
/// that says what it read.
///
/// It is there for the page faults of the first touch of the memory a pass
/// works in, which a process pays once and which fell on whichever section
/// ran first; the header of this file has the numbers. What it gives is
/// what the pass with every field is then checked against, since a file
/// whose source had no alleles to give has no such column.
fn the_first_pass(bytes: &[u8]) -> Result<(Needs, String), String> {
    let mut reader = VarsReader::new(Cursor::new(bytes)).map_err(|problem| problem.to_string())?;
    reader.set_needs(Needs::ALL);
    let mut of_the_file: Option<Needs> = None;
    let mut num_vars: u64 = 0;
    let mut sum: i64 = 0;
    while let Some(block) = reader.next_block().map_err(|problem| problem.to_string())? {
        match of_the_file {
            Some(fields) => the_fields_are(block.fields(), fields)?,
            None => of_the_file = Some(block.fields()),
        }
        // A file of more variants than a u64 counts cannot be written.
        num_vars = num_vars.saturating_add(u64::try_from(block.num_vars).unwrap_or(u64::MAX));
        sum = genotypes_added_up(&block.gts, sum);
    }
    // A file of no variants gives no block, and then the fields of the file
    // are the ones every vars file has.
    let of_the_file = of_the_file.unwrap_or(Needs::GTS);
    Ok((
        of_the_file,
        format!(
            "{num_vars} variants, the genotypes add up to {sum}, the blocks hold {of_the_file}"
        ),
    ))
}

/// One pass over the file in `bytes`, with `needs` asked for, timed from
/// the building of the reader to the last block.
///
/// Every block has to hold `asked_for`, the fields of `needs` that the file
/// has, and nothing else.
fn read_the_file(bytes: &[u8], needs: Needs, asked_for: Needs) -> Result<Run, String> {
    let started = Instant::now();
    let mut reader = VarsReader::new(Cursor::new(bytes)).map_err(|problem| problem.to_string())?;
    reader.set_needs(needs);
    let mut num_vars: u64 = 0;
    let mut sum: i64 = 0;
    while let Some(block) = reader.next_block().map_err(|problem| problem.to_string())? {
        // The check is inside the clock, where it costs one test of each of
        // the five columns of a block, four blocks in a file of this panel,
        // against the millions of genotypes of the sum beside it. A
        // consumer of blocks asks the same of every block it is given.
        the_fields_are(block.fields(), asked_for)?;
        // A file of more variants than a u64 counts cannot be written.
        num_vars = num_vars.saturating_add(u64::try_from(block.num_vars).unwrap_or(u64::MAX));
        sum = genotypes_added_up(&block.gts, sum);
    }
    let took = started.elapsed();
    Ok(Run {
        took,
        did: format!("{num_vars} variants, the genotypes add up to {sum}"),
    })
}

/// The blocks of the file in `bytes` written into a vector of bytes, timed
/// from the first `write_block` to the end of `finish`.
///
/// The blocks are read with every field asked for, before the clock starts,
/// because a block goes into the writer by value and a file of every column
/// is what `write_vars` makes.
fn write_the_blocks(bytes: &[u8]) -> Result<Run, String> {
    let mut reader = VarsReader::new(Cursor::new(bytes)).map_err(|problem| problem.to_string())?;
    reader.set_needs(Needs::ALL);
    let mut blocks = Vec::new();
    while let Some(block) = reader.next_block().map_err(|problem| problem.to_string())? {
        blocks.push(block);
    }
    let num_blocks = blocks.len();
    let individuals = reader.individuals().to_vec();
    let mut writer = VarsWriter::new(
        Vec::new(),
        &individuals,
        reader.ploidy(),
        reader.metadata().num_vars_per_block,
    )
    .map_err(|problem| problem.to_string())?;
    // The table of the names of the chromosomes of the reader the blocks
    // came from, which the writer looks the number of each variant up in.
    let chroms = reader.chroms();
    let started = Instant::now();
    for block in blocks {
        writer
            .write_block(block, chroms)
            .map_err(|problem| problem.to_string())?;
    }
    let written = writer.finish().map_err(|problem| problem.to_string())?;
    let took = started.elapsed();
    Ok(Run {
        took,
        did: format!("{num_blocks} blocks, {} bytes", written.len()),
    })
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

/// The milliseconds of a time, with the two decimals that a read of a few
/// tens of milliseconds is worth reporting to.
fn milliseconds(time: Duration) -> String {
    format!("{:.2} ms", time.as_secs_f64() * 1000.0)
}

/// It runs one of the three `runs` times, printing what each run took and
/// then the best, the median and the worst of them.
fn time_it(
    what: &str,
    runs: usize,
    mut run: impl FnMut() -> Result<Run, String>,
) -> Result<(), String> {
    let mut times = Vec::with_capacity(runs);
    for number in 1..=runs {
        let done = run()?;
        println!(
            "{what}, run {number}: {took}, {did}",
            took = milliseconds(done.took),
            did = done.did,
        );
        times.push(done.took);
    }
    let last = times.len().saturating_sub(1);
    println!(
        "{what}: best {best}, median {median}, worst {worst}",
        best = milliseconds(sorted_time(&times, 0)),
        median = milliseconds(median_time(&times)),
        worst = milliseconds(sorted_time(&times, last)),
    );
    Ok(())
}

/// What the file says about itself before any of the three is timed: what a
/// reader of it finds without reading a batch.
fn what_the_file_is(bytes: &[u8]) -> Result<String, String> {
    let reader = VarsReader::new(Cursor::new(bytes)).map_err(|problem| problem.to_string())?;
    let metadata = reader.metadata();
    Ok(format!(
        "{bytes} bytes, {batches} batches, {vars} variants, {individuals} individuals, \
         ploidy {ploidy}, {per_block} variants in a batch",
        bytes = bytes.len(),
        batches = reader.batches().len(),
        vars = reader.num_vars(),
        individuals = metadata.individuals.len(),
        ploidy = metadata.ploidy,
        per_block = metadata.num_vars_per_block,
    ))
}

/// The benchmark reads a file of the disc, and wasm has no disc. This is
/// what `cargo check --target wasm32-unknown-unknown --all-targets`
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
    let said = match what_the_file_is(&bytes) {
        Ok(said) => said,
        Err(problem) => {
            eprintln!("{path}: {problem}", path = arguments.path.display());
            return ExitCode::FAILURE;
        }
    };
    println!(
        "{path}: {said}, {runs} runs of each",
        path = arguments.path.display(),
        runs = arguments.runs,
    );
    let of_the_file = match the_first_pass(&bytes) {
        Ok((of_the_file, said)) => {
            println!("the first pass, which is not timed: {said}");
            of_the_file
        }
        Err(problem) => {
            eprintln!("{path}: {problem}", path = arguments.path.display());
            return ExitCode::FAILURE;
        }
    };
    let timed = time_it("the genotypes alone", arguments.runs, || {
        read_the_file(&bytes, Needs::GTS, Needs::GTS)
    })
    .and_then(|()| {
        time_it("every field", arguments.runs, || {
            read_the_file(&bytes, Needs::ALL, of_the_file)
        })
    })
    .and_then(|()| time_it("the write", arguments.runs, || write_the_blocks(&bytes)));
    if let Err(problem) = timed {
        eprintln!("{path}: {problem}", path = arguments.path.display());
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
