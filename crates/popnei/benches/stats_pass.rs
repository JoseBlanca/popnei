//! What the per variant pass of the stats module costs over a block of
//! genotypes that is already in memory, with no file read anywhere in the
//! timing.
//!
//! It times one thing: `calc_per_var_distribs` over a stream of blocks. The
//! benchmark builds one block once, before the clock starts, and a reader
//! written here gives that same block a stated number of times, so what is
//! timed is the loop over the rows and nothing else. What it deliberately
//! leaves out is the read of a file: opening it, the lz4 decompression of a
//! vars file and the building of a block out of the batch, which over
//! `/Users/jose/devel/popnei-bench/big.vars` is 28 in 100 of a whole pass
//! and hides a change of the row loop of less than 0.1 s in a pass of 0.48
//! s. `crates/popnei/benches/time_stats.py` is the harness that times the
//! whole pass from Python, file and all, and
//! `docs/reports/stats-measurement.md` has its numbers.
//!
//! What is inside the clock besides the row loop: the reader copies the
//! genotypes of the block for every block it gives, because a block is
//! given by value and the pass drops it. That is one copy of
//! `--vars` x `--individuals` x 2 bytes a block, 10 MB with the defaults,
//! and it stands in for no part of a real pass. `--stats none` is the pass
//! that asks for no statistic, which reads every row and calculates
//! nothing: it is what those copies and the walk of the blocks cost, and
//! the floor that every other run of the benchmark is above.
//!
//! The genotypes look like those of `big.vars`, which `make_big_vcf.py`
//! writes, because the branches the row loop takes depend on them: two
//! alleles a variant, 0 and 1, whose frequency is drawn for each variant
//! between 0.1 and 0.9; 3 in 100 genotypes missing whole, both alleles;
//! and 4 in 1000 half called, one allele called and the other not, which
//! `big.vars` does not have and which the row loop has a branch for. The
//! draws come from a generator written here with a fixed seed, so two runs
//! of the benchmark build the same block, bit for bit.
//!
//! It is run with cargo, which passes what comes after the two dashes to
//! it:
//!
//! ```text
//! cargo bench --bench stats_pass -- --stats all --pops 4 --runs 5
//! ```
//!
//! `--vars` is 5000, `--individuals` 1000, `--blocks` 20, `--runs` 5,
//! `--threads` 1, `--stats` all five and `--pops` 0 when they are not
//! given, which is 100000 variants of 1000 diploid individuals, the dataset
//! of `docs/reports/stats-measurement.md`. `--pops 0` is the one population
//! of every individual of the reader, which the pass reads as a row as it
//! lies; `--pops 4` is four populations of 250 individuals, the individuals
//! in the order they are in the block, 250 to each, which is what the pass
//! with populations of that report uses and what takes the per individual
//! path of the row loop.
//!
//! One pass that is not timed comes before the timed runs, with the same
//! statistics and the same populations. It pays the page faults of the
//! first touch of the memory a pass works in, which a process pays once.
//!
//! It prints the wall time of each run, with the variants the pass gave and
//! the mean of each statistic it asked for in the first population, and
//! then the best, the median and the worst of the times; with an even
//! number of runs the median is the middle of the two middle times. The
//! best is the number to compare between two builds, since the machine is
//! not idle and what it is doing can only make a run longer, and the worst
//! says how much it was doing something else. The run fails, and prints
//! what it expected and what it got, when the pass does not give
//! `--vars` x `--blocks` variants, so that a fixture built wrong cannot
//! pass in silence.

#![cfg_attr(
    target_family = "wasm",
    allow(
        dead_code,
        unused_imports,
        reason = "the benchmark is native: in wasm only its empty main is compiled, and \
                  what the timing is made of is left unused"
    )
)]

use std::hint::black_box;
use std::process::ExitCode;
use std::time::{Duration, Instant};

use popnei::block::{Block, BlockReader};
use popnei::filters::FilteringStats;
use popnei::stats::{
    DEFAULT_HIST_RANGE, DEFAULT_MIN_NUM_INDIVIDUALS, DEFAULT_NUM_BINS, DEFAULT_POLY_THRESHOLD,
    ExpHet, HistBins, Maf, ObsHet, PerVarDistribs, PerVarDistribsConfig, PerVarStat, Pops,
    calc_per_var_distribs,
};
use popnei::variant::{ChromTable, MISSING_ALLELE, Needs};

/// How many variants one block holds when the command line does not say.
const DEFAULT_VARS: usize = 5000;

/// How many individuals the dataset has when the command line does not say.
const DEFAULT_INDIVIDUALS: usize = 1000;

/// How many times the reader gives its block when the command line does not
/// say. 20 blocks of 5000 variants are the 100000 variants of `big.vars`.
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

/// What `--stats` is given to ask for no statistic at all: the pass that
/// reads every row and calculates nothing.
const NO_STAT: &str = "none";

/// What `--stats` is given to ask for the five statistics, which is what it
/// is when the command line does not name it.
const ALL_STATS: &str = "all";

/// What the command line asked for.
struct Arguments {
    num_vars: usize,
    num_individuals: usize,
    num_blocks: usize,
    runs: usize,
    threads: usize,
    /// The statistics to calculate, which is empty for `--stats none`.
    stats: Vec<PerVarStat>,
    /// How many populations the individuals are split into, 0 for the one
    /// population of every individual of the reader.
    num_pops: usize,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
stats_pass [--vars n] [--individuals n] [--blocks n] [--runs n]
           [--threads n] [--stats names] [--pops n]

It times whole passes of `calc_per_var_distribs` over blocks that are
already in memory: one block is built before the clock starts and a reader
gives that same block --blocks times, so no file is read and what is timed
is the loop over the rows.

  --vars n          how many variants one block holds, 5000 by default
  --individuals n   how many individuals the dataset has, 1000 by default
  --blocks n        how many times the reader gives its block, 20 by default
  --runs n          how many times the pass is timed, 5 by default
  --threads n       how many threads the pool it runs in has, 1 by default
  --stats names     the statistics, by the names a user writes, separated
                    by commas: obs_het, maf, exp_het, unbiased_exp_het,
                    poly_vars_ratio; `all` for the five, which is the
                    default, and `none` for the pass that calculates
                    nothing
  --pops n          how many populations the individuals are split into,
                    each getting the same number of them; 0, the default,
                    is the one population of every individual
  --help            this

The genotypes are drawn from a fixed seed and look like those of
/Users/jose/devel/popnei-bench/big.vars: two alleles a variant whose
frequency is between 0.1 and 0.9, 3 in 100 genotypes missing whole and 4 in
1000 half called. One pass that is not timed comes first, so that the timed
runs do not pay the page faults of the first touch of the memory a pass
works in. It prints the wall time of each run, with the variants the pass
gave and the mean of each statistic in the first population, and then the
best, the median and the worst of the times. The run fails when the pass
does not give --vars x --blocks variants.";

/// The number that comes after `name` on the command line, or the message
/// that says what should have come after it.
fn number_after(name: &str, args: &mut impl Iterator<Item = String>) -> Result<usize, String> {
    args.next()
        .ok_or_else(|| format!("{name} takes a number and none came after it"))?
        .parse::<usize>()
        .map_err(|_| format!("{name} takes a number"))
}

/// The statistics `--stats` named: the five for `all`, none for `none`, and
/// otherwise the ones whose names are in the list, separated by commas.
fn stats_of(named: &str) -> Result<Vec<PerVarStat>, String> {
    if named == NO_STAT {
        return Ok(Vec::new());
    }
    if named == ALL_STATS {
        return PerVarStat::NAMES
            .iter()
            .map(|name| PerVarStat::of_name(name).map_err(|error| error.to_string()))
            .collect();
    }
    named
        .split(',')
        .map(|name| PerVarStat::of_name(name.trim()).map_err(|error| error.to_string()))
        .collect()
}

/// What to run, or the message that says what the command line should have
/// been.
fn arguments() -> Result<Arguments, String> {
    let mut num_vars = DEFAULT_VARS;
    let mut num_individuals = DEFAULT_INDIVIDUALS;
    let mut num_blocks = DEFAULT_BLOCKS;
    let mut runs = DEFAULT_RUNS;
    let mut threads = DEFAULT_THREADS;
    let mut stats = stats_of(ALL_STATS)?;
    let mut num_pops = 0;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--vars" => num_vars = number_after("--vars", &mut args)?,
            "--individuals" => num_individuals = number_after("--individuals", &mut args)?,
            "--blocks" => num_blocks = number_after("--blocks", &mut args)?,
            "--runs" => runs = number_after("--runs", &mut args)?,
            "--threads" => threads = number_after("--threads", &mut args)?,
            "--pops" => num_pops = number_after("--pops", &mut args)?,
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
            // being dropped, as `filter_vars.rs` refuses it: a `--pops=4`
            // that was dropped leaves a run that timed the pass with one
            // population and calls it the pass with four.
            other => {
                return Err(format!(
                    "`{other}` is not an argument of this benchmark\n\n{USAGE}"
                ));
            }
        }
    }
    if num_vars == 0 || num_individuals == 0 || num_blocks == 0 || runs == 0 || threads == 0 {
        return Err(
            "--vars, --individuals, --blocks, --runs and --threads are 1 or more".to_owned(),
        );
    }
    if num_pops != 0 && num_individuals.checked_rem(num_pops) != Some(0) {
        return Err(format!(
            "the {num_individuals} individuals do not divide into {num_pops} populations \
             of the same size"
        ));
    }
    Ok(Arguments {
        num_vars,
        num_individuals,
        num_blocks,
        runs,
        threads,
        stats,
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
/// by value and the pass drops it; that copy is inside the clock, and the
/// doc comment of this file says what it is worth.
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

/// The block the benchmark was asked for and the names of its individuals,
/// built once and read by every run.
struct Fixture {
    gts: Vec<i8>,
    num_vars: usize,
    num_individuals: usize,
    num_blocks: usize,
    individuals: Vec<String>,
}

impl Fixture {
    /// The block of `num_vars` rows of `num_individuals` individuals, to be
    /// given `num_blocks` times, with the individuals named `ind0` and
    /// upwards, which is how a population names them.
    fn of(num_vars: usize, num_individuals: usize, num_blocks: usize) -> Result<Fixture, String> {
        Ok(Fixture {
            gts: gts_of(num_vars, num_individuals)?,
            num_vars,
            num_individuals,
            num_blocks,
            individuals: (0..num_individuals)
                .map(|individual| format!("ind{individual}"))
                .collect(),
        })
    }

    /// A reader over it, which building costs no copy: the genotypes are
    /// copied when a block is given and not before.
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
}

/// The populations of the pass: the one population of every individual for
/// `num_pops` of 0, and otherwise `num_pops` populations of the same number
/// of individuals, the individuals in the order they are in the block.
fn pops_of(num_pops: usize, individuals: &[String]) -> Result<Pops, String> {
    if num_pops == 0 {
        return Ok(Pops::all(individuals.len()));
    }
    let of_a_pop = individuals
        .len()
        .checked_div(num_pops)
        .ok_or_else(|| "a pass of no population".to_owned())?;
    if of_a_pop == 0 {
        return Err(format!(
            "{num_pops} populations of the {num} individuals leave one of them empty",
            num = individuals.len()
        ));
    }
    let named: Vec<(String, Vec<String>)> = individuals
        .chunks(of_a_pop)
        .enumerate()
        .map(|(number, of_the_pop)| (format!("pop{number}"), of_the_pop.to_vec()))
        .collect();
    Pops::from_names(&named, individuals).map_err(|error| error.to_string())
}

/// What one pass calculates: the statistics the command line named, for
/// those populations, in the bins and with the thresholds a user gets when
/// they name none of them.
fn config_of(stats: &[PerVarStat], pops: Pops) -> Result<PerVarDistribsConfig, String> {
    let (start, end) = DEFAULT_HIST_RANGE;
    Ok(PerVarDistribsConfig {
        stats: stats.to_vec(),
        pops,
        bins: HistBins::linear(start, end, DEFAULT_NUM_BINS).map_err(|error| error.to_string())?,
        obs_het: ObsHet::new(DEFAULT_MIN_NUM_INDIVIDUALS),
        maf: Maf::new(PLOIDY, DEFAULT_MIN_NUM_INDIVIDUALS).map_err(|error| error.to_string())?,
        exp_het: ExpHet::new(PLOIDY, PLOIDY, DEFAULT_MIN_NUM_INDIVIDUALS)
            .map_err(|error| error.to_string())?,
        poly_threshold: DEFAULT_POLY_THRESHOLD,
    })
}

/// One run: how long the pass took and the line that says what it gave, the
/// variants and the mean of each statistic in the first population.
struct Run {
    took: Duration,
    did: String,
}

/// What the pass gave, for the line the run prints: the variants and the
/// mean of every statistic that was asked for, in the first population.
///
/// A statistic nobody asked for is not in the line. The polymorphism ratio
/// is there as the variants that were polymorphic in that population.
fn what_it_gave(distribs: &PerVarDistribs) -> String {
    let mut did = format!("{num_vars} variants", num_vars = distribs.num_vars);
    for (name, mean) in [
        ("obs_het", distribs.obs_het.as_ref()),
        ("maf", distribs.maf.as_ref()),
        ("exp_het", distribs.exp_het.as_ref()),
        ("unbiased_exp_het", distribs.unbiased_exp_het.as_ref()),
    ] {
        if let Some(mean) = mean.and_then(|distrib| distrib.mean(0)) {
            did.push_str(&format!(", {name} {mean:.6}"));
        }
    }
    if let Some(poly) = distribs.poly_vars_ratio.as_ref() {
        did.push_str(&format!(
            ", polymorphic {num} variants",
            num = poly.num_poly(0)
        ));
    }
    did
}

/// One whole pass over the blocks of `fixture`, timed from the first block
/// to the result, and what it gave.
///
/// The reader is built before the clock starts, which costs no copy, and
/// the genotypes of each block are copied inside it, as the doc comment of
/// this file says. Both the reader and the result go through
/// [`black_box`], so that nothing of the pass is dropped for being unread.
fn one_pass(fixture: &Fixture, config: &PerVarDistribsConfig) -> popnei::Result<Run> {
    let mut reader = fixture.reader();
    let started = Instant::now();
    let distribs = calc_per_var_distribs(black_box(&mut reader), config)?;
    let took = started.elapsed();
    let distribs = black_box(distribs);
    Ok(Run {
        took,
        did: what_it_gave(&distribs),
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

/// The names of the statistics of the run, for the line that says what was
/// timed.
fn named(stats: &[PerVarStat]) -> String {
    if stats.is_empty() {
        return "no statistic".to_owned();
    }
    stats
        .iter()
        .map(|stat| stat.name())
        .collect::<Vec<&str>>()
        .join(", ")
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
    ) {
        Ok(fixture) => fixture,
        Err(problem) => {
            eprintln!("{problem}");
            return ExitCode::FAILURE;
        }
    };
    let config = match pops_of(arguments.num_pops, &fixture.individuals)
        .and_then(|pops| config_of(&arguments.stats, pops))
    {
        Ok(config) => config,
        Err(problem) => {
            eprintln!("{problem}");
            return ExitCode::FAILURE;
        }
    };
    let expected = match fixture.variants_of_a_pass() {
        Ok(expected) => expected,
        Err(problem) => {
            eprintln!("{problem}");
            return ExitCode::FAILURE;
        }
    };
    println!(
        "{blocks} blocks of {vars} variants of {individuals} individuals of the ploidy \
         {PLOIDY}, built in {built}, {threads} threads, {runs} runs, {pops}, {stats}",
        blocks = arguments.num_blocks,
        vars = arguments.num_vars,
        individuals = arguments.num_individuals,
        built = seconds(built.elapsed()),
        threads = arguments.threads,
        runs = arguments.runs,
        pops = match arguments.num_pops {
            0 => "one population of every individual".to_owned(),
            num_pops => format!("{num_pops} populations of the same size"),
        },
        stats = named(&arguments.stats),
    );
    let mut times = Vec::with_capacity(arguments.runs);
    // The pass that is not timed, and then the timed ones. Both go through
    // `install`, so that the one that warms the memory runs on the same
    // pool as the ones that are timed.
    for run in 0..=arguments.runs {
        let done = match pool.install(|| one_pass(&fixture, &config)) {
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
        if run == 0 {
            println!(
                "the first pass, which is not timed: {did}, in {took}",
                did = done.did,
                took = seconds(done.took),
            );
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
