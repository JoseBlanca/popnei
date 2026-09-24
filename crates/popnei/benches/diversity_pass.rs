//! What the pass of the diversity module costs over a block of genotypes
//! that is already in memory, in time and in the bytes it holds live at
//! once, with no file read anywhere in the timing.
//!
//! It times one thing: `calc_pop_diversity` over a stream of blocks. The
//! benchmark builds one block once, before the clock starts, and a reader
//! written here gives that same block a stated number of times, so what is
//! timed is the loop over the rows and nothing else. What it deliberately
//! leaves out is the read of a file: opening it, the lz4 decompression of a
//! vars file and the building of a block out of the batch, which over
//! `/Users/jose/devel/popnei-bench/big.vars` is 28 in 100 of a whole pass of
//! the per variant statistics. `crates/popnei/benches/time_diversity.py` is
//! the harness that times the whole pass from Python, file and all.
//! `crates/popnei/benches/stats_pass.rs` is the same benchmark for the pass
//! of the stats module, and this one has its shape with another call in it.
//!
//! What is inside the clock besides the row loop: the reader copies the
//! genotypes of the block for every block it gives, because a block is
//! given by value and the pass drops it. That is one copy of
//! `--vars` x `--individuals` x 2 bytes a block, 10 MB with the defaults,
//! and it stands in for no part of a real pass.
//!
//! ## The bytes it holds live
//!
//! A global allocator written here adds the size of every allocation and
//! subtracts the size of every free, and keeps the largest total it ever
//! saw. Each run reads that total from the moment before the pass starts to
//! the moment it ends, so what is printed is the most bytes the process held
//! live at once across the pass, and the same number comes out of every run.
//! `/usr/bin/time -l` cannot be used for this: its maximum resident set is
//! the whole process, the binary and the pages the allocator keeps among
//! them, and it cannot separate the few tens of kilobytes of the counts of
//! the chunks from the block of genotypes they sit beside.
//!
//! One block of `--vars` x `--individuals` x 2 bytes is live at the peak,
//! since the pass holds one block at a time and the reader has just given
//! it. It is printed on its own and subtracted, so that the last figure is
//! what the pass holds **beside the block**, which is the number section 2
//! of `docs/architecture.md` puts a block of about 10 MB against.
//!
//! Counting costs two atomic operations per allocation. The pass allocates
//! per chunk of 64 rows and not per row, so what it adds to a pass of
//! thousands of rows is not visible beside the row loop; it has not been
//! measured, and a time from this benchmark is not compared with a time from
//! a build without the counting allocator.
//!
//! The counts of one chunk of rows live until the chunks before it have been
//! added, so how many are read at once decides the memory, and
//! `diversity::chunks_of_a_group` makes that two chunks for each thread of
//! the pool. `--threads` therefore changes what this prints, and a memory
//! figure is stated with the threads it was taken at.
//!
//! ## The genotypes
//!
//! They look like those of `big.vars`, which `make_big_vcf.py` writes,
//! because the branches the row loop takes depend on them: two alleles a
//! variant, 0 and 1, whose frequency is drawn for each variant between 0.1
//! and 0.9; 3 in 100 genotypes missing whole, both alleles; and 4 in 1000
//! half called, one allele called and the other not, which `big.vars` does
//! not have and which the row loop has a branch for. The draws come from a
//! generator written here with a fixed seed, so two runs of the benchmark
//! build the same block, bit for bit.
//!
//! ## Running it
//!
//! It is run with cargo, which passes what comes after the two dashes to
//! it:
//!
//! ```text
//! cargo bench --bench diversity_pass -- --stats all --draw 200 --pops 3
//! ```
//!
//! `--individuals` is 1000, `--blocks` 20, `--runs` 5, `--threads` 1,
//! `--stats` all five, `--pops` 0 and `--draw` none when they are not
//! given, and `--vars` is then `block::default_num_vars_per_block` of that
//! many individuals, 5000 rows for 1000 of them. Twenty such blocks are the
//! 100000 variants of 1000 diploid individuals of
//! `/Users/jose/devel/popnei-bench/big.vars`. `--pops 0` is the one
//! population of every individual of the reader, which the pass reads as a
//! row as it lies; `--pops 3` is three populations of the individuals in the
//! order they are in the block, a third to each.
//!
//! One pass that is not timed comes before the timed runs, with the same
//! statistics, the same draw and the same populations. It pays the page
//! faults of the first touch of the memory a pass works in, which a process
//! pays once.
//!
//! It prints the wall time of each run, with the variants the pass gave and
//! what each statistic it asked for came to in the first population, and
//! then the best, the median and the worst of the times; with an even number
//! of runs the median is the middle of the two middle times. The best is the
//! number to compare between two builds, since the machine is not idle and
//! what it is doing can only make a run longer, and the worst says how much
//! it was doing something else.
//!
//! The run fails, and prints what it expected and what it got, in three
//! cases, so that a fixture or a command line built wrong cannot pass in
//! silence: when the pass does not give `--vars` x `--blocks` variants; when
//! a statistic that `--stats` named has no value in the result, which is
//! what a pass asked for the wrong statistic gives; and when `--stats` named
//! a statistic of `DiversityStats::NAMES_AND_STATS` that this file has no
//! line reading off the result, so that it could say nothing about it.
//!
//! ## The two runs in which the pass calculates nothing per population
//!
//! A variant counts for a population when that population called
//! `stats::DEFAULT_MIN_NUM_INDIVIDUALS` genotypes at it, 20, which is the
//! `min_num_individuals` of every run of this benchmark and is not a
//! command line argument. So `--individuals 500 --pops 50` gives 10
//! individuals a population, no variant counts for any of them, and the pass
//! calculates no statistic for any population while it reads every row and
//! allocates its bins. A `--draw` above what any population calls is the
//! same thing one step later: every variant counts, no variant reaches the
//! draw, and every standardized value is NaN and every bin of the spectrum
//! 0.
//!
//! Neither fails, and neither should: three of the four shapes of the memory
//! table of "The memory" of `docs/specs/diversity.md` are one of the two,
//! and what that table measures is the bins the pass allocates, which it
//! allocates whatever the genotypes come to. What a run like that is not is
//! a timing of the arithmetic it skipped. So the run says which of the two
//! it is, in a line of its own after the untimed pass and again in the line
//! of every timed run, and a time from it is read with that line beside it.

#![expect(
    unsafe_code,
    reason = "a global allocator is what counts the bytes a pass holds live, and \
              `GlobalAlloc` is an unsafe trait whose four methods are unsafe: each of \
              them hands the allocation itself to the allocator of the system, with the \
              layout and the pointer it was given, and only reads the size to count it. \
              The lint is expected for the whole file because the trait cannot be \
              implemented without it, and this is a benchmark and not code the library \
              ships: the core crate keeps its `forbid(unsafe_code)`"
)]
#![cfg_attr(
    target_family = "wasm",
    allow(
        dead_code,
        unused_imports,
        reason = "the benchmark is native: in wasm only its empty main is compiled, and \
                  what the timing is made of is left unused"
    )
)]

use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};

use popnei::block::{Block, BlockReader, default_num_vars_per_block};
use popnei::diversity::{DiversityOptions, DiversityStats, PopDiversity, calc_pop_diversity};
use popnei::filters::FilteringStats;
use popnei::stats::DEFAULT_MIN_NUM_INDIVIDUALS;
use popnei::variant::{ChromTable, MISSING_ALLELE, Needs};

/// How many bytes the program has asked the allocator for and not given
/// back.
static LIVE_BYTES: AtomicUsize = AtomicUsize::new(0);

/// The largest [`LIVE_BYTES`] ever reached since it was last set, which is
/// what one run reads to say how much the pass held at once.
static PEAK_BYTES: AtomicUsize = AtomicUsize::new(0);

/// The allocator of the benchmark: the one of the system, with the bytes it
/// is asked for added up and the largest total it ever held kept.
///
/// It is the whole process that is counted and not the pass alone, so a run
/// reads the peak from the moment before the pass starts, which
/// [`start_counting_again`] sets.
struct CountsWhatIsLive;

impl CountsWhatIsLive {
    /// It records `bytes` more as live, and the peak with them.
    fn took(bytes: usize) {
        let live = LIVE_BYTES
            .fetch_add(bytes, Ordering::Relaxed)
            .saturating_add(bytes);
        PEAK_BYTES.fetch_max(live, Ordering::Relaxed);
    }

    /// It records `bytes` as given back.
    fn gave_back(bytes: usize) {
        LIVE_BYTES.fetch_sub(bytes, Ordering::Relaxed);
    }
}

// SAFETY: every method hands the allocation itself to `System`, which is a
// `GlobalAlloc`, with the layout and the pointer it was given, and does
// nothing to the memory but count its size. The counters are atomics, so
// the threads of the pool add to them without a race.
unsafe impl GlobalAlloc for CountsWhatIsLive {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc(layout) };
        if !pointer.is_null() {
            CountsWhatIsLive::took(layout.size());
        }
        pointer
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        CountsWhatIsLive::gave_back(layout.size());
        unsafe { System.dealloc(pointer, layout) };
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc_zeroed(layout) };
        if !pointer.is_null() {
            CountsWhatIsLive::took(layout.size());
        }
        pointer
    }

    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let moved = unsafe { System.realloc(pointer, layout, new_size) };
        if !moved.is_null() {
            match new_size.checked_sub(layout.size()) {
                Some(more) => CountsWhatIsLive::took(more),
                None => CountsWhatIsLive::gave_back(layout.size().saturating_sub(new_size)),
            }
        }
        moved
    }
}

#[global_allocator]
static ALLOCATOR: CountsWhatIsLive = CountsWhatIsLive;

/// The bytes that are live now, with the peak set to them, which is what a
/// run does right before the pass it measures.
fn start_counting_again() -> usize {
    let live = LIVE_BYTES.load(Ordering::Relaxed);
    PEAK_BYTES.store(live, Ordering::Relaxed);
    live
}

/// How many individuals the dataset has when the command line does not say.
const DEFAULT_INDIVIDUALS: usize = 1000;

/// How many times the reader gives its block when the command line does not
/// say. 20 blocks of the 5000 rows popnei chooses for 1000 individuals are
/// the 100000 variants of `big.vars`.
const DEFAULT_BLOCKS: usize = 20;

/// How many times the pass is timed when the command line does not say.
const DEFAULT_RUNS: usize = 5;

/// How many threads the pool has when the command line does not say.
const DEFAULT_THREADS: usize = 1;

/// The ploidy of the genotypes the benchmark builds.
const PLOIDY: usize = 2;

/// The seed of the generator that draws the genotypes, so that two runs of
/// the benchmark build the same block.
const SEED: u64 = 0x05EE_D0F5_7A75;

/// The share of the genotypes that are missing whole, both of their
/// alleles, which is the rate `make_big_vcf.py` writes into `big.vcf`.
const MISSING_RATE: f64 = 0.03;

/// The share of the genotypes with one allele called and the other not.
/// `big.vars` has none, and the row loop has a branch for them.
const HALF_CALLED_RATE: f64 = 0.004;

/// The smallest frequency of the alternative allele of a variant.
const LOWEST_FREQUENCY: f64 = 0.1;

/// The largest frequency of the alternative allele of a variant, which is
/// the range `make_big_vcf.py` draws its ancestral frequencies from.
const HIGHEST_FREQUENCY: f64 = 0.9;

/// What `--stats` is given to ask for the five statistics, which is what it
/// is when the command line does not name it.
const ALL_STATS: &str = "all";

/// What `--stats` is given to ask for the four statistics that need no
/// draw, which is what a Python or a TypeScript user who names none gets.
const WITHOUT_A_DRAW: &str = "without_a_draw";

/// How many bytes are written as one megabyte, which is the unit every
/// memory figure of this benchmark is printed in.
const BYTES_OF_A_MEGABYTE: f64 = 1_000_000.0;

/// What the command line asked for.
struct Arguments {
    num_vars: usize,
    num_individuals: usize,
    num_blocks: usize,
    runs: usize,
    threads: usize,
    /// The statistics to calculate, which is never empty: this module
    /// refuses a pass that is asked for no statistic.
    stats: DiversityStats,
    /// The called alleles every population is brought down to, and `None`
    /// for a pass with no draw.
    num_called_alleles: Option<u32>,
    /// How many populations the individuals are split into, 0 for the one
    /// population of every individual of the reader.
    num_pops: usize,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
diversity_pass [--vars n] [--individuals n] [--blocks n] [--runs n]
               [--threads n] [--stats names] [--draw n] [--pops n]

It times whole passes of `calc_pop_diversity` over blocks that are already
in memory, and counts the bytes each pass holds live at once: one block is
built before the clock starts and a reader gives that same block --blocks
times, so no file is read and what is timed is the loop over the rows.

  --vars n          how many variants one block holds; the block popnei
                    chooses for --individuals individuals by default, from
                    block::default_num_vars_per_block, which is 5000 rows
                    for 1000 of them
  --individuals n   how many individuals the dataset has, 1000 by default
  --blocks n        how many times the reader gives its block, 20 by default
  --runs n          how many times the pass is timed, 5 by default
  --threads n       how many threads the pool it runs in has, 1 by default;
                    it changes the memory as well as the time, the chunks of
                    a block being read two at a time for each thread
  --stats names     the statistics, by the names a user writes, separated by
                    commas: num_alleles, private_alleles,
                    variable_vars_ratio, folded_sfs, fis; `all` for the five,
                    which is the default, and `without_a_draw` for the four
                    that need no --draw
  --draw n          the called alleles every population is brought down to,
                    which the folded spectrum needs and the standardized
                    values of the other three read; none by default
  --pops n          how many populations the individuals are split into, cut
                    as evenly as their number allows, the individuals in the
                    order they are in the block; 0, the default, is the one
                    population of every individual
  --help            this

The genotypes are drawn from a fixed seed and look like those of
/Users/jose/devel/popnei-bench/big.vars: two alleles a variant whose
frequency is between 0.1 and 0.9, 3 in 100 genotypes missing whole and 4 in
1000 half called. One pass that is not timed comes first, so that the timed
runs do not pay the page faults of the first touch of the memory a pass
works in. It prints the wall time of each run, with the variants the pass
gave and what each statistic came to in the first population, the bytes the
pass held live at once beside its block, and then the best, the median and
the worst of the times. The run fails when the pass does not give --vars x
--blocks variants, when a statistic --stats named has no value, and when
--stats named a statistic this benchmark has no line reading.

Every run asks for a min_num_individuals of 20, the called genotypes a
population needs at a variant for the variant to count for it, which is not
an argument here. So --individuals 500 --pops 50 leaves 10 individuals a
population and no variant counts for any of them, and a --draw above what
any population calls leaves no variant in the draw. Neither is refused, and
three shapes of the memory table of docs/specs/diversity.md are one of them,
but in neither does the pass do the arithmetic a time would be read as: the
run says which of the two it is, after the untimed pass and in the line of
every timed run.";

/// The number that comes after `name` on the command line, or the message
/// that says what should have come after it.
fn number_after(name: &str, args: &mut impl Iterator<Item = String>) -> Result<usize, String> {
    let written = args
        .next()
        .ok_or_else(|| format!("{name} takes a number and none came after it"))?;
    written
        .parse::<usize>()
        .map_err(|_| format!("{name} takes a number of 0 or more, and {written:?} came after it"))
}

/// The statistics `--stats` named: the five for `all`, the four that need no
/// draw for `without_a_draw`, and otherwise the ones whose names are in the
/// list, separated by commas.
fn stats_of(named: &str) -> Result<DiversityStats, String> {
    if named == ALL_STATS {
        return Ok(DiversityStats::ALL);
    }
    if named == WITHOUT_A_DRAW {
        return Ok(DiversityStats::WITHOUT_A_DRAW);
    }
    let mut stats = DiversityStats::empty();
    for name in named.split(',') {
        stats |= DiversityStats::of_name(name.trim()).map_err(|error| error.to_string())?;
    }
    Ok(stats)
}

/// What to run, or the message that says what the command line should have
/// been.
fn arguments() -> Result<Arguments, String> {
    let mut num_vars = None;
    let mut num_individuals = DEFAULT_INDIVIDUALS;
    let mut num_blocks = DEFAULT_BLOCKS;
    let mut runs = DEFAULT_RUNS;
    let mut threads = DEFAULT_THREADS;
    let mut stats = stats_of(ALL_STATS)?;
    let mut num_called_alleles = None;
    let mut num_pops = 0;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--vars" => num_vars = Some(number_after("--vars", &mut args)?),
            "--individuals" => num_individuals = number_after("--individuals", &mut args)?,
            "--blocks" => num_blocks = number_after("--blocks", &mut args)?,
            "--runs" => runs = number_after("--runs", &mut args)?,
            "--threads" => threads = number_after("--threads", &mut args)?,
            "--pops" => num_pops = number_after("--pops", &mut args)?,
            "--draw" => {
                let drawn = number_after("--draw", &mut args)?;
                num_called_alleles = Some(
                    u32::try_from(drawn).map_err(|_| format!("--draw is at most {}", u32::MAX))?,
                );
            }
            "--stats" => {
                stats = stats_of(
                    &args
                        .next()
                        .ok_or_else(|| "--stats takes names and none came after it".to_owned())?,
                )?;
            }
            // `cargo bench` adds this to the command line of every bench,
            // to tell a harness that has tests too to run its benchmarks.
            // This one has only this benchmark and takes it as nothing.
            "--bench" => {}
            "--help" | "-h" => return Err(USAGE.to_owned()),
            // An argument that is not one of those is refused instead of
            // being dropped, as `stats_pass.rs` refuses it: a `--pops=3`
            // that was dropped leaves a run that timed the pass with one
            // population and calls it the pass with three.
            other => {
                return Err(format!(
                    "`{other}` is not an argument of this benchmark\n\n{USAGE}"
                ));
            }
        }
    }
    // The block the reader of a file would give for that many individuals,
    // asked of the library and not worked out again here: a change to the
    // rule would otherwise leave the benchmark measuring a block that no
    // reader gives.
    let num_vars = num_vars.unwrap_or_else(|| default_num_vars_per_block(num_individuals));
    if num_vars == 0 || num_individuals == 0 || num_blocks == 0 || runs == 0 || threads == 0 {
        return Err(
            "--vars, --individuals, --blocks, --runs and --threads are 1 or more".to_owned(),
        );
    }
    Ok(Arguments {
        num_vars,
        num_individuals,
        num_blocks,
        runs,
        threads,
        stats,
        num_called_alleles,
        num_pops,
    })
}

/// The generator the genotypes are drawn from: splitmix64, which is a few
/// lines and has no dependency, with the seed of the benchmark, so that
/// every run of it builds the same block.
struct Random {
    state: u64,
}

impl Random {
    /// The generator at its seed.
    fn of(seed: u64) -> Random {
        Random { state: seed }
    }

    /// The next draw, a number from 0 included to 1 excluded.
    fn share(&mut self) -> f64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut drawn = self.state;
        drawn = (drawn ^ drawn.wrapping_shr(30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        drawn = (drawn ^ drawn.wrapping_shr(27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        drawn ^= drawn.wrapping_shr(31);
        // The top 53 bits of it, which is what an f64 holds exactly, over
        // 2 to the 53: a number from 0 included to 1 excluded.
        let of_53_bits = drawn.wrapping_shr(11);
        // 2^53 is 9007199254740992, and no draw of 53 bits reaches it.
        of_53_bits as f64 / 9_007_199_254_740_992_f64
    }
}

/// The allele a draw gives at a variant whose alternative allele has the
/// frequency `frequency`: 1 below it and 0 at or above it.
fn allele_of(drawn: f64, frequency: f64) -> i8 {
    if drawn < frequency { 1 } else { 0 }
}

/// The genotypes of one block: `num_vars` rows of `num_individuals`
/// diploid individuals, variant after variant and inside a variant
/// individual after individual, drawn as the doc comment of this file
/// says.
///
/// It is built once, before any clock starts.
fn gts_of(num_vars: usize, num_individuals: usize) -> Result<Vec<i8>, String> {
    let num_alleles = num_vars
        .checked_mul(num_individuals)
        .and_then(|genotypes| genotypes.checked_mul(PLOIDY))
        .ok_or_else(|| {
            format!(
                "{num_vars} variants of {num_individuals} individuals of the ploidy {PLOIDY} \
                 are more alleles than this machine addresses"
            )
        })?;
    let mut random = Random::of(SEED);
    let mut gts: Vec<i8> = Vec::with_capacity(num_alleles);
    for _of_the_vars in 0..num_vars {
        let frequency = LOWEST_FREQUENCY + (HIGHEST_FREQUENCY - LOWEST_FREQUENCY) * random.share();
        for _of_the_individuals in 0..num_individuals {
            let kind = random.share();
            let first = allele_of(random.share(), frequency);
            let second = allele_of(random.share(), frequency);
            if kind < MISSING_RATE {
                gts.push(MISSING_ALLELE);
                gts.push(MISSING_ALLELE);
            } else if kind < MISSING_RATE + HALF_CALLED_RATE {
                gts.push(first);
                gts.push(MISSING_ALLELE);
            } else {
                gts.push(first);
                gts.push(second);
            }
        }
    }
    Ok(gts)
}

/// A reader that gives one block of genotypes a stated number of times and
/// then nothing, which is what puts a stream of blocks in front of the pass
/// with no file behind it.
///
/// The crate's own reader of blocks in memory is behind `#[cfg(test)]` and
/// a benchmark cannot use it. This one holds the genotypes it was built
/// with and copies them into every block it gives, because a block is given
/// by value and the pass drops it; that copy is inside the clock and inside
/// what the allocator counts, and the doc comment of this file says what it
/// is worth.
struct TheSameBlockAgain<'a> {
    gts: &'a [i8],
    num_vars: usize,
    num_individuals: usize,
    individuals: Vec<String>,
    chroms: ChromTable,
    /// How many more blocks it gives.
    blocks_left: usize,
}

impl BlockReader for TheSameBlockAgain<'_> {
    fn next_block(&mut self) -> popnei::Result<Option<Block>> {
        let Some(left) = self.blocks_left.checked_sub(1) else {
            return Ok(None);
        };
        self.blocks_left = left;
        Ok(Some(Block {
            num_vars: self.num_vars,
            num_individuals: self.num_individuals,
            ploidy: PLOIDY,
            gts: self.gts.to_vec(),
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        }))
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

    // Every block it gives holds the genotypes and no column, which is
    // what the pass asks for, so there is nothing to remember here.
    fn set_needs(&mut self, _needs: Needs) {}

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        Vec::new()
    }
}

/// The block the benchmark was asked for and the populations of its
/// individuals, built once and read by every run.
struct Fixture {
    gts: Vec<i8>,
    num_vars: usize,
    num_individuals: usize,
    num_blocks: usize,
    individuals: Vec<String>,
    /// The indices of the individuals of each population, which is what the
    /// pass takes, and empty for the one population of every individual.
    pops: Vec<Vec<usize>>,
}

impl Fixture {
    /// The block of `num_vars` rows of `num_individuals` individuals, to be
    /// given `num_blocks` times, with its individuals cut into `num_pops`
    /// populations as [`pops_of`] cuts them.
    fn of(
        num_vars: usize,
        num_individuals: usize,
        num_blocks: usize,
        num_pops: usize,
    ) -> Result<Fixture, String> {
        Ok(Fixture {
            gts: gts_of(num_vars, num_individuals)?,
            num_vars,
            num_individuals,
            num_blocks,
            individuals: (0..num_individuals)
                .map(|individual| format!("ind{individual}"))
                .collect(),
            pops: pops_of(num_pops, num_individuals)?,
        })
    }

    /// A reader over it, which building costs no copy of the genotypes:
    /// they are copied when a block is given and not before.
    fn reader(&self) -> TheSameBlockAgain<'_> {
        TheSameBlockAgain {
            gts: &self.gts,
            num_vars: self.num_vars,
            num_individuals: self.num_individuals,
            individuals: self.individuals.clone(),
            chroms: ChromTable::new(),
            blocks_left: self.num_blocks,
        }
    }

    /// The populations as the pass takes them, a slice of slices of the
    /// indices of the individuals.
    fn pops_of_the_pass(&self) -> Vec<&[usize]> {
        self.pops.iter().map(Vec::as_slice).collect()
    }

    /// How many variants a whole pass over it gives, which every run checks
    /// the pass against.
    fn variants_of_a_pass(&self) -> Result<u64, String> {
        let variants = self
            .num_vars
            .checked_mul(self.num_blocks)
            .ok_or_else(|| "the variants of a pass are more than this machine counts".to_owned())?;
        u64::try_from(variants)
            .map_err(|_| "the variants of a pass are more than a u64 counts".to_owned())
    }

    /// How many bytes the genotypes of one block are, which is what the
    /// reader copies for every block it gives and what is live beside the
    /// counts of the pass when the peak is reached.
    fn bytes_of_a_block(&self) -> Result<usize, String> {
        self.num_vars
            .checked_mul(self.num_individuals)
            .and_then(|genotypes| genotypes.checked_mul(PLOIDY))
            .ok_or_else(|| "a block of more bytes than this machine addresses".to_owned())
    }
}

/// The indices of the individuals of each population: nothing for
/// `num_pops` of 0, which the pass reads as the one population of every
/// individual, and otherwise `num_pops` populations, the individuals in the
/// order they are in the block and cut as evenly as their number allows.
///
/// Population `pop` gets the individuals from `pop * num_individuals /
/// num_pops` up to the next such bound, so there are always `num_pops` of
/// them and no two differ by more than one individual: the 200 of the panel
/// in 3 populations are 67, 67 and 66. An even cut is not asked for, because
/// what a population costs the pass is its individuals and its variants and
/// the panel, which the measurement uses, has populations of 48, 68 and 84.
///
/// # Errors
///
/// More populations than individuals, which would leave one with none, and
/// a product of the two above what a `usize` counts.
fn pops_of(num_pops: usize, num_individuals: usize) -> Result<Vec<Vec<usize>>, String> {
    if num_pops == 0 {
        return Ok(Vec::new());
    }
    if num_pops > num_individuals {
        return Err(format!(
            "{num_pops} populations of the {num_individuals} individuals leave one of them empty"
        ));
    }
    let bound_of = |pop: usize| -> Result<usize, String> {
        pop.checked_mul(num_individuals)
            .and_then(|product| product.checked_div(num_pops))
            .ok_or_else(|| {
                format!(
                    "{num_pops} populations of {num_individuals} individuals are more than \
                     this machine counts"
                )
            })
    };
    let mut pops = Vec::with_capacity(num_pops);
    for pop in 0..num_pops {
        let first = bound_of(pop)?;
        let past_the_last = bound_of(pop.saturating_add(1))?;
        pops.push((first..past_the_last).collect());
    }
    Ok(pops)
}

/// One run: how long the pass took, the line that says what it gave, the
/// bytes it held live at once beside its block, and what the genotypes left
/// the pass with nothing to do.
struct Run {
    took: Duration,
    did: String,
    /// The most bytes the process held live at once across the pass, less
    /// what was live before it started and less the genotypes of one block.
    bytes_beside_the_block: usize,
    /// A statistic that `--stats` named and the result gave no value for,
    /// which is what the run fails on.
    statistic_that_fails: Option<AStatisticThatFails>,
    /// The per population arithmetic the genotypes left the pass with
    /// nothing to do, which the run says and does not fail on.
    left_out: WhatTheDataLeftOut,
}

/// Why the run fails on a statistic that `--stats` named.
enum AStatisticThatFails {
    /// The result has no value for it, which is what a pass asked for the
    /// wrong statistic gives.
    WithNoValue(&'static str),
    /// [`has_a_value`] does not read it, so the run can say nothing about
    /// whether the pass gave it. A statistic added to
    /// [`DiversityStats::NAMES_AND_STATS`] lands here until the line that
    /// reads it is written here.
    ThisBenchmarkDoesNotRead(&'static str),
}

/// Whether the result holds a value of `stat` for the first population, and
/// `None` for a statistic of [`DiversityStats::NAMES_AND_STATS`] that this
/// file does not read.
///
/// The value of a statistic that was asked for is `Some`, whatever the
/// genotypes came to: a population that no variant counted for has 0 or NaN
/// there and not `None`.
///
/// Each statistic is looked up by itself and not by its place in that
/// table, because the table is public and its order is nobody's contract:
/// a check that zipped it against a list written here in the same order
/// would match nothing the moment the table was reordered, and would then
/// answer that no statistic is missing whatever the result holds.
fn has_a_value(diversity: &PopDiversity, stat: DiversityStats) -> Option<bool> {
    if stat == DiversityStats::NUM_ALLELES {
        Some(diversity.num_alleles(0).is_some())
    } else if stat == DiversityStats::PRIVATE_ALLELES {
        Some(diversity.private_alleles(0).is_some())
    } else if stat == DiversityStats::VARIABLE_VARS_RATIO {
        Some(diversity.num_variable_vars(0).is_some())
    } else if stat == DiversityStats::FOLDED_SFS {
        Some(diversity.folded_sfs(0).is_some())
    } else if stat == DiversityStats::FIS {
        Some(diversity.fis(0).is_some())
    } else {
        None
    }
}

/// The first statistic that `stats` named and the run fails on: one the
/// result has no value for, or one this file cannot read.
fn a_statistic_that_fails(
    diversity: &PopDiversity,
    stats: DiversityStats,
) -> Option<AStatisticThatFails> {
    DiversityStats::NAMES_AND_STATS
        .iter()
        .filter(|(_, stat)| stats.contains(*stat))
        .find_map(|(name, stat)| match has_a_value(diversity, *stat) {
            Some(true) => None,
            Some(false) => Some(AStatisticThatFails::WithNoValue(name)),
            None => Some(AStatisticThatFails::ThisBenchmarkDoesNotRead(name)),
        })
}

/// The per population arithmetic that the genotypes of a pass left with
/// nothing to do, although the pass finished and was timed.
///
/// Both are legal runs and neither fails: three of the four shapes of the
/// memory table of `docs/specs/diversity.md` are one of them, and the bins
/// they measure are allocated whatever the data. What they are not is a
/// timing of the arithmetic they skipped, and a time read without them says
/// that it is.
struct WhatTheDataLeftOut {
    /// How many populations the pass had, which is 1 for the one population
    /// of every individual.
    num_pops: usize,
    /// No variant counted for any population, so no statistic was
    /// calculated for one: a variant counts for a population when the
    /// population called `min_num_individuals` genotypes at it, and this
    /// benchmark asks for [`DEFAULT_MIN_NUM_INDIVIDUALS`] of them.
    no_variant_counted: bool,
    /// A draw was asked for and no variant reached it in any population, so
    /// every standardized value is NaN and every bin of the spectrum is 0.
    no_variant_in_the_draw: bool,
}

impl WhatTheDataLeftOut {
    /// What the pass did not do over these genotypes, as the sentences the
    /// run prints on their own lines, and nothing when it did all of it.
    fn said(&self, num_called_alleles: Option<u32>) -> Vec<String> {
        let num_pops = self.num_pops;
        let mut said = Vec::new();
        if self.no_variant_counted {
            said.push(format!(
                "no variant counted for any of the {num_pops} populations, a variant counting \
                 for a population that called {DEFAULT_MIN_NUM_INDIVIDUALS} genotypes at it: \
                 the pass read every row and allocated its bins, and calculated no statistic \
                 for any population, so these times are not times of that arithmetic"
            ));
        }
        if let (true, Some(drawn)) = (self.no_variant_in_the_draw, num_called_alleles) {
            said.push(format!(
                "no variant reached the draw of {drawn} called alleles in any of the \
                 {num_pops} populations: every standardized value is NaN and every bin of the \
                 spectrum is 0, so these times are not times of the draw"
            ));
        }
        said
    }

    /// The same for the line of one run, which has no room for the reason.
    fn in_the_line_of_a_run(&self) -> String {
        let mut marks = String::new();
        if self.no_variant_counted {
            marks.push_str(", no variant counted for any population");
        }
        if self.no_variant_in_the_draw {
            marks.push_str(", no variant in the draw for any population");
        }
        marks
    }
}

/// What the genotypes left the pass with nothing to do, read off the
/// result: `drew` says whether a draw was asked for, since a pass with none
/// has no variant in a draw and that is not a thing to report.
fn what_the_data_left_out(diversity: &PopDiversity, drew: bool) -> WhatTheDataLeftOut {
    let num_pops = diversity.num_pops();
    let of_every_pop = |counted: &dyn Fn(usize) -> Option<u64>| {
        num_pops > 0 && (0..num_pops).all(|pop| counted(pop) == Some(0))
    };
    WhatTheDataLeftOut {
        num_pops,
        no_variant_counted: of_every_pop(&|pop| diversity.num_vars(pop)),
        no_variant_in_the_draw: drew && of_every_pop(&|pop| diversity.num_vars_in_draw(pop)),
    }
}

/// What the pass gave, for the line the run prints: the variants of the
/// pass and, for each statistic that was asked for, what the first
/// population came to.
///
/// A statistic nobody asked for is not in the line. The spectrum is there as
/// its bins and the sum over them, which is the variants in the draw for the
/// population.
fn what_it_gave(diversity: &PopDiversity, stats: DiversityStats) -> String {
    let mut did = format!(
        "{num_vars} variants, {in_draw} of them in the draw for the first population",
        num_vars = diversity.num_vars_of_the_pass(),
        in_draw = diversity
            .num_vars_in_draw(0)
            .map_or_else(|| "none".to_owned(), |num| num.to_string()),
    );
    if stats.contains(DiversityStats::NUM_ALLELES) {
        did.push_str(&format!(
            ", num_alleles {total:?} in draw {in_draw:?}",
            total = diversity.num_alleles(0),
            in_draw = diversity.num_alleles_in_draw(0),
        ));
    }
    if stats.contains(DiversityStats::PRIVATE_ALLELES) {
        did.push_str(&format!(
            ", private_alleles {total:?} in draw {in_draw:?}",
            total = diversity.private_alleles(0),
            in_draw = diversity.private_alleles_in_draw(0),
        ));
    }
    if stats.contains(DiversityStats::VARIABLE_VARS_RATIO) {
        did.push_str(&format!(
            ", variable_vars {total:?} in draw {in_draw:?}",
            total = diversity.num_variable_vars(0),
            in_draw = diversity.variable_vars_ratio_in_draw(0),
        ));
    }
    if stats.contains(DiversityStats::FOLDED_SFS) {
        let bins = diversity.folded_sfs(0);
        did.push_str(&format!(
            ", folded_sfs of {num_bins} bins summing to {sum:?}",
            num_bins = bins.map_or(0, <[f64]>::len),
            sum = bins.map(|bins| bins.iter().sum::<f64>()),
        ));
    }
    if stats.contains(DiversityStats::FIS) {
        did.push_str(&format!(", fis {fis:?}", fis = diversity.fis(0)));
    }
    did
}

/// One whole pass over the blocks of `fixture`, timed from the first block
/// to the result, with what it gave and the bytes it held.
///
/// The reader is built before the clock starts, which costs no copy of the
/// genotypes, and the genotypes of each block are copied inside it, as the
/// doc comment of this file says. Both the reader and the result go through
/// [`black_box`], so that nothing of the pass is dropped for being unread.
fn one_pass(
    fixture: &Fixture,
    options: &DiversityOptions,
    bytes_of_a_block: usize,
) -> popnei::Result<Run> {
    let pops = fixture.pops_of_the_pass();
    let mut reader = fixture.reader();
    let live_before = start_counting_again();
    let started = Instant::now();
    let diversity = calc_pop_diversity(black_box(&mut reader), &pops, options)?;
    let took = started.elapsed();
    let peak = PEAK_BYTES.load(Ordering::Relaxed);
    let diversity = black_box(diversity);
    Ok(Run {
        took,
        did: what_it_gave(&diversity, options.stats),
        bytes_beside_the_block: peak
            .saturating_sub(live_before)
            .saturating_sub(bytes_of_a_block),
        statistic_that_fails: a_statistic_that_fails(&diversity, options.stats),
        left_out: what_the_data_left_out(&diversity, options.num_called_alleles.is_some()),
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

/// The seconds of a time, with the three decimals that a pass of a few
/// tenths of a second is worth reporting to.
fn seconds(time: Duration) -> String {
    format!("{:.3} s", time.as_secs_f64())
}

/// The megabytes of a count of bytes, with the three decimals that tell
/// 0.045 MB from 0.046.
fn megabytes(bytes: usize) -> String {
    format!("{:.3} MB", bytes as f64 / BYTES_OF_A_MEGABYTE)
}

/// The benchmark builds a pool of threads, and wasm has none; rayon is not
/// a dependency of the wasm targets either. This is what `cargo check
/// --target wasm32-unknown-unknown --all-targets` compiles of it, so that
/// the command which checks that nothing of the crate has left wasm behind
/// can check the benchmarks too.
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
    let built = Instant::now();
    let fixture = match Fixture::of(
        arguments.num_vars,
        arguments.num_individuals,
        arguments.num_blocks,
        arguments.num_pops,
    ) {
        Ok(fixture) => fixture,
        Err(problem) => {
            eprintln!("{problem}");
            return ExitCode::FAILURE;
        }
    };
    let options = DiversityOptions {
        stats: arguments.stats,
        num_called_alleles: arguments.num_called_alleles,
        min_num_individuals: DEFAULT_MIN_NUM_INDIVIDUALS,
    };
    let (expected, bytes_of_a_block) = match fixture.variants_of_a_pass().and_then(|expected| {
        fixture
            .bytes_of_a_block()
            .map(|bytes_of_a_block| (expected, bytes_of_a_block))
    }) {
        Ok(both) => both,
        Err(problem) => {
            eprintln!("{problem}");
            return ExitCode::FAILURE;
        }
    };
    println!(
        "{blocks} blocks of {vars} variants of {individuals} individuals of the ploidy \
         {PLOIDY}, {block_bytes} of genotypes a block, built in {built}, {threads} threads, \
         {runs} runs, {pops}, {stats}, {draw}",
        blocks = arguments.num_blocks,
        vars = arguments.num_vars,
        individuals = arguments.num_individuals,
        block_bytes = megabytes(bytes_of_a_block),
        built = seconds(built.elapsed()),
        threads = arguments.threads,
        runs = arguments.runs,
        pops = match arguments.num_pops {
            0 => "one population of every individual".to_owned(),
            num_pops => format!("{num_pops} populations"),
        },
        stats = arguments.stats.names().join(", "),
        draw = match arguments.num_called_alleles {
            None => "no draw".to_owned(),
            Some(drawn) => format!("a draw of {drawn} called alleles"),
        },
    );
    let mut times = Vec::with_capacity(arguments.runs);
    // The pass that is not timed, and then the timed ones. Both go through
    // `install`, so that the one that warms the memory runs on the same
    // pool as the ones that are timed.
    for run in 0..=arguments.runs {
        let done = match pool.install(|| one_pass(&fixture, &options, bytes_of_a_block)) {
            Ok(done) => done,
            Err(error) => {
                eprintln!("the pass: {error}");
                return ExitCode::FAILURE;
            }
        };
        if !done.did.starts_with(&format!("{expected} variants")) {
            eprintln!(
                "the pass was to give {expected} variants and gave: {did}",
                did = done.did
            );
            return ExitCode::FAILURE;
        }
        match done.statistic_that_fails {
            Some(AStatisticThatFails::WithNoValue(name)) => {
                eprintln!(
                    "the pass was asked for {stats} and gave no value for {name}: {did}",
                    stats = arguments.stats.names().join(", "),
                    did = done.did,
                );
                return ExitCode::FAILURE;
            }
            Some(AStatisticThatFails::ThisBenchmarkDoesNotRead(name)) => {
                eprintln!(
                    "the pass was asked for {stats}, and this benchmark has no line that reads \
                     {name} off the result, so it cannot say whether the pass gave it: {did}",
                    stats = arguments.stats.names().join(", "),
                    did = done.did,
                );
                return ExitCode::FAILURE;
            }
            None => {}
        }
        if run == 0 {
            println!(
                "the first pass, which is not timed: {did}, in {took}, holding {held} beside \
                 its block",
                did = done.did,
                took = seconds(done.took),
                held = megabytes(done.bytes_beside_the_block),
            );
            for said in done.left_out.said(arguments.num_called_alleles) {
                println!("{said}");
            }
            continue;
        }
        println!(
            "run {run}: {took}, {held} beside its block, {did}{left_out}",
            took = seconds(done.took),
            held = megabytes(done.bytes_beside_the_block),
            did = done.did,
            left_out = done.left_out.in_the_line_of_a_run(),
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
