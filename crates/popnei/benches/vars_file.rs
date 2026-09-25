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
//! The threads. The reader of a vars file decompresses several batches at
//! once on the threads of rayon, and everything else a pass does, the checks
//! of each batch and the building of each block, runs on the thread that
//! calls it, as does the whole of the write. The threads are those of a pool
//! this benchmark builds, of as many threads as `--threads` says, and every
//! timed section runs inside its `install`, so `--threads 1` is a pass with
//! no thread but the caller's.
//!
//! The digest. A pass has to give the same blocks, in the same order,
//! whatever the threads: the numbers of the chromosomes are handed out in
//! the order in which the names are first seen, and nine places of popnei
//! index their results by a running count of the variants, so a block out of
//! order is a wrong result and not a slower one. The sum of the genotypes
//! that each run prints is the same whatever order the blocks come in, so it
//! cannot see that. The untimed first pass, which asks for every field,
//! therefore also feeds a 64 bit FNV-1a with the place of each block, its
//! counts and the bytes of every column it holds, and then the names of the
//! chromosome table in the order of their numbers, and prints it: two passes
//! that give the same values in another order do not agree on it. Nothing of
//! it is inside a clock.
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

use popnei::block::{Block, BlockReader};
use popnei::io::vars::{VarsReader, VarsWriter};
use popnei::variant::{ChromTable, Needs};

/// How many times each of the three is timed when the command line does not
/// say.
const DEFAULT_RUNS: usize = 5;

/// How many threads the pool has when the command line does not say: one,
/// the number the 21 ms of `docs/specs/io_vars.md` are stated on.
const DEFAULT_THREADS: usize = 1;

/// What the command line asked for.
struct Arguments {
    path: PathBuf,
    threads: usize,
    runs: usize,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
vars_file <path to a vars file> [--threads n] [--runs n]

It times three things on that file, each of them with the file already in
memory: a pass over it with the genotypes alone asked for, a pass with
every field asked for, and the writing of its blocks, read before the
clock starts, into a vector of bytes. Every pass goes through
`next_block` to the end of the file and adds up the genotypes of each
block, which it prints, so that the work is not dropped as unused. One
pass that is not timed comes before the three, because the first touch
of the memory a pass works in costs page faults that a process pays
once, and it fails when a block holds a column that nobody asked for.

  --threads n  how many threads the pool it runs in has, 1 by default
  --runs n     how many times each of the three is timed, 5 by default
  --help       this

It also prints, from the first pass and outside every clock, a digest of
the whole pass that depends on the order of the blocks: the place of each
block, its counts and the bytes of its columns, and the names of the
chromosomes in the order of their numbers. Two runs of the same file at
different numbers of threads have to print the same digest.

It prints the wall time of each run and then the best, the median and the
worst of them. The best of the pass with the genotypes alone is what
`docs/specs/io_vars.md` asks for 21 ms of, on a file of 1000 individuals
and 20000 variants; the header of this file says how that file is made.";

/// What to run, or the message that says what the command line should have
/// been.
fn arguments() -> Result<Arguments, String> {
    let mut path: Option<PathBuf> = None;
    let mut threads = DEFAULT_THREADS;
    let mut runs = DEFAULT_RUNS;
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
    if threads == 0 || runs == 0 {
        return Err("--threads and --runs are 1 or more".to_owned());
    }
    Ok(Arguments {
        path,
        threads,
        runs,
    })
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

/// What a 64 bit FNV-1a starts from, and what each byte multiplies it by.
const FNV_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

/// A digest of a whole pass that depends on the order of its blocks: a 64
/// bit FNV-1a fed with the place of each block, its counts and the bytes of
/// every column it holds, and then with the names of the chromosome table in
/// the order of their numbers.
///
/// It is what says that a pass on many threads gives the blocks a pass on one
/// thread gives, in the same order. The sum of the genotypes each run prints
/// cannot say it: addition does not care in which order the blocks came, and
/// neither does the count of the variants. Here the place of a block goes
/// into the hash before its values, so two passes that give the same values
/// in another order do not agree; and the names of the chromosomes go in by
/// their numbers, so a pass that numbered them in another order does not
/// agree either, which is what a batch decoded out of order would do.
///
/// FNV-1a is six lines and needs no dependency. Its multiplication wraps by
/// design, which is what `wrapping_mul` says here: this is a hash and not a
/// count.
struct Digest {
    hash: u64,
}

impl Digest {
    /// A digest that nothing has been fed to.
    fn new() -> Digest {
        Digest {
            hash: FNV_OFFSET_BASIS,
        }
    }

    /// One byte.
    fn byte(&mut self, byte: u8) {
        self.hash = (self.hash ^ u64::from(byte)).wrapping_mul(FNV_PRIME);
    }

    /// Every byte, in the order they lie in.
    fn bytes(&mut self, bytes: &[u8]) {
        for byte in bytes {
            self.byte(*byte);
        }
    }

    /// A number, as its eight bytes, so that 1 and 256 differ.
    fn number(&mut self, number: u64) {
        self.bytes(&number.to_le_bytes());
    }

    /// A count of this machine, which a pass of this machine is digested on.
    fn count(&mut self, count: usize) {
        self.number(u64::try_from(count).unwrap_or(u64::MAX));
    }

    /// A text, its length first, so that "ab" and "a" then "b" differ.
    fn text(&mut self, text: &str) {
        self.count(text.len());
        self.bytes(text.as_bytes());
    }

    /// Whether a column is there, before its bytes, so that a block without
    /// a column and a block whose column is empty differ.
    fn there(&mut self, there: bool) {
        self.byte(u8::from(there));
    }

    /// The block that is the `place`th of its pass, counted from 1: its
    /// place, its three counts and the bytes of each of its columns.
    ///
    /// Every field of the block is named here, as `Block::fields` names them
    /// all, so that a column added later does not fall out of the digest
    /// without the compiler saying so.
    fn block(&mut self, place: u64, block: &Block) {
        let Block {
            num_vars,
            num_individuals,
            ploidy,
            gts,
            chrom,
            pos,
            id,
            alleles,
            qual,
        } = block;
        self.number(place);
        self.count(*num_vars);
        self.count(*num_individuals);
        self.count(*ploidy);
        self.count(gts.len());
        for allele in gts {
            // An allele is -1 up to the last alternative one, and the byte it
            // lies in is what the digest takes: no `as` narrows anything.
            self.bytes(&allele.to_le_bytes());
        }
        self.there(chrom.is_some());
        if let Some(chrom) = chrom {
            self.count(chrom.len());
            for number in chrom {
                self.number(u64::from(*number));
            }
        }
        self.there(pos.is_some());
        if let Some(pos) = pos {
            self.count(pos.len());
            for position in pos {
                self.number(*position);
            }
        }
        self.there(id.is_some());
        if let Some(id) = id {
            self.count(id.len());
            for of_a_variant in id {
                self.text(of_a_variant);
            }
        }
        self.there(alleles.is_some());
        if let Some(alleles) = alleles {
            self.count(alleles.num_vars());
            for var in 0..alleles.num_vars() {
                self.count(alleles.num_alleles(var));
                for allele in 0..alleles.num_alleles(var) {
                    self.text(alleles.allele(var, allele));
                }
            }
        }
        self.there(qual.is_some());
        if let Some(qual) = qual {
            self.count(qual.len());
            for quality in qual {
                // The bits, because a quality that no variant has is NaN and
                // NaN is equal to nothing, itself included.
                self.number(u64::from(quality.to_bits()));
            }
        }
    }

    /// The names of the chromosome table, in the order of their numbers,
    /// which is the order in which the pass first saw them.
    fn chroms(&mut self, chroms: &ChromTable) {
        self.count(chroms.len());
        for number in 0..chroms.len() {
            let name = u32::try_from(number)
                .ok()
                .and_then(|number| chroms.name(number));
            self.there(name.is_some());
            if let Some(name) = name {
                self.text(name);
            }
        }
    }

    /// What was fed, as the sixteen hexadecimal digits a run prints.
    fn printed(&self) -> String {
        format!("{hash:016x}", hash = self.hash)
    }
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
///
/// It is also where the digest of the pass is taken, which is why it asks
/// for every field: a digest of the genotypes alone would not see the
/// chromosomes, whose numbers are what an inter-batch renumbering shows in.
fn the_first_pass(bytes: &[u8]) -> Result<(Needs, String), String> {
    let mut reader = VarsReader::new(Cursor::new(bytes)).map_err(|problem| problem.to_string())?;
    reader.set_needs(Needs::ALL);
    let mut of_the_file: Option<Needs> = None;
    let mut num_vars: u64 = 0;
    let mut sum: i64 = 0;
    let mut digest = Digest::new();
    let mut num_blocks: u64 = 0;
    while let Some(block) = reader.next_block().map_err(|problem| problem.to_string())? {
        match of_the_file {
            Some(fields) => the_fields_are(block.fields(), fields)?,
            None => of_the_file = Some(block.fields()),
        }
        // A file of more blocks than a u64 counts cannot be written.
        num_blocks = num_blocks.saturating_add(1);
        digest.block(num_blocks, &block);
        // A file of more variants than a u64 counts cannot be written.
        num_vars = num_vars.saturating_add(u64::try_from(block.num_vars).unwrap_or(u64::MAX));
        sum = genotypes_added_up(&block.gts, sum);
    }
    digest.chroms(reader.chroms());
    // A file of no variants gives no block, and then the fields of the file
    // are the ones every vars file has.
    let of_the_file = of_the_file.unwrap_or(Needs::GTS);
    Ok((
        of_the_file,
        format!(
            "{num_vars} variants in {num_blocks} blocks, the genotypes add up to {sum}, \
             the blocks hold {of_the_file}, the digest of the pass is {digest}",
            digest = digest.printed()
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

/// The benchmark builds a pool of threads and reads a file of the disc, and
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
    let said = match what_the_file_is(&bytes) {
        Ok(said) => said,
        Err(problem) => {
            eprintln!("{path}: {problem}", path = arguments.path.display());
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
        "{path}: {said}, {threads} threads, {runs} runs of each",
        path = arguments.path.display(),
        threads = arguments.threads,
        runs = arguments.runs,
    );
    let of_the_file = match pool.install(|| the_first_pass(&bytes)) {
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
        pool.install(|| read_the_file(&bytes, Needs::GTS, Needs::GTS))
    })
    .and_then(|()| {
        time_it("every field", arguments.runs, || {
            pool.install(|| read_the_file(&bytes, Needs::ALL, of_the_file))
        })
    })
    .and_then(|()| {
        time_it("the write", arguments.runs, || {
            pool.install(|| write_the_blocks(&bytes))
        })
    });
    if let Err(problem) = timed {
        eprintln!("{path}: {problem}", path = arguments.path.display());
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
