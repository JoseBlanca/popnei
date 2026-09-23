//! How long the seven distances between every pair of populations take, the
//! `calc_pop_dist_sums` of `docs/specs/dists.md`, over a source whose blocks
//! are already in memory and over a vars file.
//!
//! "Speed" of that spec gives that pass no number to reach. It names the
//! dataset, 100000 variants of 1000 diploid individuals, biallelic, 3 in 100
//! genotypes missing, cut into 3 populations and into 20, so that how the
//! cost grows with the square of the populations is measured and not
//! guessed, and it leaves the measurement to a performance review. Nothing
//! had timed the module before this benchmark, and the numbers of that
//! review come from it.
//!
//! The two settings it times, each of them on the pool of threads of the
//! process and then inside a pool of one thread, and each of them for every
//! populations file the command line gives:
//!
//! - **the blocks in memory**: a reader written here gives 20 blocks of 5000
//!   variants of 1000 diploid individuals, biallelic, 3 in 100 genotypes
//!   missing whole, drawn by a generator with a fixed seed, so that the
//!   blocks are the same on every run and on every machine. They are built
//!   before the clock starts and handed out one at a time, so nothing of a
//!   reader, no file, no decompression and no copy, is inside the clock.
//!   This is the setting a target of the review is to be of.
//! - **the vars file**: the same calculation over a `VarsReader` on the bytes
//!   of the file, which are read from the disc before the clock starts. What
//!   is inside the clock is what a user waits for, the reader and the
//!   calculation together, and its difference from the setting above is what
//!   the reader of a vars file costs.
//!
//! The blocks in memory carry the chromosome and the position of every
//! variant and not the genotypes alone, which is what tells this benchmark
//! from `kosman_dists.rs`: this pass cuts the variants into the resampling
//! groups, the stretches of a chromosome that its standard errors are
//! resampled over, and it cuts them from the positions. They are laid out as
//! `big.vars` is, two chromosomes of 50000 variants each with the variants
//! 1000 base pairs apart, so that a length of 1000000 base pairs cuts each
//! chromosome into 50 groups and the two settings into 100 each.
//!
//! The populations. A `--pops` file holds one line for each individual, its
//! name and the name of its population with a tab between them, which
//! `make_pops.py` beside this file writes and whose doc comment says how:
//! `pops3.tsv` cuts the 1000 individuals into 3 populations of 296, 356 and
//! 348, which make 3 pairs, and `pops20.tsv` into 20 of 50, which make 190.
//! The individuals of the blocks in memory are named `s000` to `s999`, which
//! is what `big.vars` calls them, so one file serves both settings, and the
//! populations are in the order they first appear in the file. Every timing
//! is run for each file given, so two files are 8 timings.
//!
//! The genotypes of the two settings are not the same genotypes: the file
//! holds what `make_big_vcf.py` simulated, three populations with an F_ST of
//! 0.1 between them, and the blocks in memory what the generator here drew,
//! which is one frequency for all 1000 individuals at each variant and so
//! populations that differ only by the noise of the draw. The work does not
//! depend on the values: what a pair costs is the alleles the variant holds,
//! two in both settings, and both settings leave every population above the
//! threshold of called genotypes at every variant. The F_ST of the two
//! settings differs, and the one of the first pair of populations, which is
//! the first two the file names, is printed for each of them.
//!
//! What it checks. The call has to succeed, to give the pairs the
//! populations make, 3 for 3 populations and 190 for 20, and to have read
//! the variants of its setting: the 100000 of the blocks in memory, and for
//! the vars file the variants its footer states, which is what the pass
//! reads when no filter is in the chain. A run that gave other pairs or
//! other variants is refused with a message instead of being timed: a timing
//! of other work that looks reasonable is worse than no timing. Nothing else
//! about the numbers is asserted here; `pop_dists.rs` has the tests that
//! assert the distances.
//!
//! It is run with cargo, which passes what comes after the two dashes to it,
//! and the release build is what `cargo bench` makes:
//!
//! ```text
//! cargo bench --bench pop_dists -- /Users/jose/devel/popnei-bench/big.vars \
//!     --pops /Users/jose/devel/popnei-bench/pops3.tsv \
//!     --pops /Users/jose/devel/popnei-bench/pops20.tsv \
//!     --runs 5
//! ```
//!
//! `--runs` is 5 when it is not given, and `--groups` cuts the variants into
//! stretches of 1000000 base pairs. `--groups` is on the command line
//! because what the resampling groups cost is one of the questions the
//! review asks: `none` gives no standard errors and asks the reader for the
//! genotypes alone, `variant` makes a group of each variant, and `bp:n`
//! makes stretches of n base pairs. The memory of the sums is 48 bytes for
//! each pair and each group, so `variant` over 100000 variants is 14 MB for
//! 3 populations and 912 MB for 20, and the `f2_groups` a user of the
//! packages gets on top of it is another 76 MB for 20: `--groups variant` is
//! not for the dataset of 20 populations on a machine that is also running
//! something else.
//!
//! Before the timed runs of each timing there is one run that is not timed,
//! whose time is printed beside them: the first touch of the memory a
//! calculation works in costs page faults that a process pays once, and they
//! would otherwise all fall on whichever timing ran first.
//!
//! The load average of the machine is not read here: a timing that is
//! reported is taken with nothing else running, and `sysctl -n vm.loadavg`
//! before and after the invocation is what says the machine was quiet.
//!
//! How the files are made: `make_big_vcf.py` beside this file writes the VCF
//! of 100000 variants and popnei's `write_vars` the vars file of it, as
//! `docs/reports/dists-kosman-measurement.md` has the two commands, and
//! `make_pops.py` writes the two populations files.

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
use popnei::filters::FilteringStats;
use popnei::io::vars::VarsReader;
use popnei::pop_dists::{JackknifeGroups, PopDistMeasure, PopDistOptions, calc_pop_dist_sums};
use popnei::stats::{DEFAULT_MIN_NUM_INDIVIDUALS, Pops};
use popnei::variant::{ChromTable, MISSING_ALLELE, Needs};

/// How many times each timing is run when the command line does not say.
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

/// Over how many chromosomes those variants are spread, which is what
/// `make_big_vcf.py` gives `big.vars`: the blocks are shared out evenly
/// between them, the first half of them on the first chromosome.
const NUM_CHROMS: usize = 2;

/// The names of those chromosomes, which are the ones `big.vars` holds and
/// which a group of the result is named by.
const CHROM_NAMES: [&str; NUM_CHROMS] = ["chr1", "chr2"];

/// How many base pairs apart the variants of one chromosome lie, which is
/// what `make_big_vcf.py` writes: the variant at the place i of its
/// chromosome sits at 1000 (i + 1), 1 based as in a VCF.
const BETWEEN_VARS: u64 = 1000;

/// How long a resampling group is when the command line does not say: a
/// stretch of a chromosome of 1000000 base pairs, which cuts each
/// chromosome of the dataset into 50 groups.
const DEFAULT_GROUP_LENGTH: u64 = 1_000_000;

/// The seed of the generator that draws the genotypes of the blocks in
/// memory. Any value does: what it is for is that the blocks are the same
/// on every run and on every machine.
const SEED: u64 = 42;

/// Out of 1024 genotypes, how many are missing whole: 31 in 1024 is 3.03 in
/// 100, the 3 in 100 of the dataset of "Speed".
const MISSING_IN_1024: u64 = 31;

/// What the command line asked for.
struct Arguments {
    /// The vars file of the second setting.
    path: PathBuf,
    /// The populations files, each of which every timing is run for.
    pops: Vec<PathBuf>,
    /// How many times each timing is run.
    runs: usize,
    /// How the variants are cut into the resampling groups.
    groups: JackknifeGroups,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
pop_dists <path to a vars file> --pops <path> [--pops <path>] [--runs n] [--groups how]

It times the seven distances between every pair of populations,
`calc_pop_dist_sums`, in two settings, each of them on the pool of threads of
the process and then inside a pool of one thread: over 20 blocks of 5000
variants of 1000 diploid individuals that a generator with a fixed seed drew
and that are built before the clock starts, so that no reader is inside it;
and over the vars file given, whose bytes are read from the disc before the
clock starts. Both settings are timed for every populations file, so two of
them make 8 timings. The blocks in memory are the dataset of \"Speed\" of
docs/specs/dists.md, 100000 variants of 1000 individuals on 2 chromosomes.

  --pops p   a file of one line per individual, its name and the name of its
             population with a tab between them, which make_pops.py writes.
             One at least, and it may be given more than once
  --runs n   how many times each timing is run, 5 by default
  --groups h how the variants are cut into the resampling groups: `none` for
             no standard errors, `variant` for a group of each variant, or
             `bp:n` for stretches of n base pairs, which is bp:1000000 by
             default. The sums are 48 bytes for each pair and each group, so
             `variant` over 100000 variants of 20 populations holds 912 MB
  --help     this

It prints the wall time of each run and then the best, the median and the
worst of them, with the pairs the call gave, the variants it read, the groups
they fell into and the F_ST of the first two populations of the file.";

/// What to run, or the message that says what the command line should have
/// been.
fn arguments() -> Result<Arguments, String> {
    let mut path: Option<PathBuf> = None;
    let mut pops: Vec<PathBuf> = Vec::new();
    let mut runs = DEFAULT_RUNS;
    let mut groups = JackknifeGroups::OfBasePairs(DEFAULT_GROUP_LENGTH);
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--pops" => {
                let given = args
                    .next()
                    .ok_or_else(|| "--pops takes a path and none came after it".to_owned())?;
                pops.push(PathBuf::from(given));
            }
            "--runs" => {
                runs = args
                    .next()
                    .ok_or_else(|| "--runs takes a number and none came after it".to_owned())?
                    .parse::<usize>()
                    .map_err(|_| "--runs takes a number".to_owned())?;
            }
            "--groups" => {
                let given = args.next().ok_or_else(|| {
                    "--groups takes none, variant or bp:n and none came after it".to_owned()
                })?;
                groups = groups_of(&given)?;
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
    if pops.is_empty() {
        return Err(format!("no populations file was given\n\n{USAGE}"));
    }
    if runs == 0 {
        return Err("--runs is 1 or more".to_owned());
    }
    Ok(Arguments {
        path,
        pops,
        runs,
        groups,
    })
}

/// How the variants are cut into the resampling groups, from what `--groups`
/// was given: `none`, `variant`, or `bp:` and a length in base pairs.
fn groups_of(how: &str) -> Result<JackknifeGroups, String> {
    if how == "none" {
        return Ok(JackknifeGroups::None);
    }
    if how == "variant" {
        return Ok(JackknifeGroups::PerVariant);
    }
    let Some(length) = how.strip_prefix("bp:") else {
        return Err(format!(
            "`{how}` is not a way of cutting the groups: none, variant or bp:n"
        ));
    };
    let length = length
        .parse::<u64>()
        .map_err(|_| format!("`bp:{length}` takes a number of base pairs"))?;
    if length == 0 {
        return Err("bp:0 is no stretch of a chromosome: 1 base pair at least".to_owned());
    }
    Ok(JackknifeGroups::OfBasePairs(length))
}

/// The populations of a file of one line per individual, the name of the
/// individual and the name of its population with a tab between them, in
/// the order the populations first appear in the file.
///
/// The order is the order the result gives the populations and the pairs
/// in, and it is the order of a user's `pops` in Python, which is the order
/// of the keys of their dict.
fn pops_of_the_file(text: &str) -> Result<Vec<(String, Vec<String>)>, String> {
    let mut pops: Vec<(String, Vec<String>)> = Vec::new();
    for (number, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let Some((individual, pop)) = line.split_once('\t') else {
            let line_number = number.saturating_add(1);
            return Err(format!(
                "line {line_number} holds no tab: a line is the name of an individual, a tab \
                 and the name of its population"
            ));
        };
        match pops.iter_mut().find(|(name, _)| name == pop) {
            Some((_, individuals)) => individuals.push(individual.to_owned()),
            None => pops.push((pop.to_owned(), vec![individual.to_owned()])),
        }
    }
    if pops.is_empty() {
        return Err("the file names no individual".to_owned());
    }
    Ok(pops)
}

/// How many pairs `num_pops` populations make, and 0 when they make more
/// than a `usize` counts.
fn num_pairs_of(num_pops: usize) -> usize {
    num_pops
        .saturating_sub(1)
        .saturating_mul(num_pops)
        .checked_div(2)
        .unwrap_or_default()
}

/// The names of the individuals of the blocks in memory, `s000` to `s999`,
/// which are the names `make_big_vcf.py` writes and the ones the
/// populations files hold.
fn names_of_the_individuals() -> Vec<String> {
    (0..NUM_INDIVIDUALS)
        .map(|number| format!("s{number:03}"))
        .collect()
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
/// `make_big_vcf.py`, and what a variant costs this pass is the alleles it
/// holds and not their values.
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

/// The genotypes of one block of 5000 variants of 1000 diploid biallelic
/// individuals, with 3 in 100 of them missing whole, drawn by a generator
/// at [`SEED`].
///
/// Each variant gets a frequency of the alternative allele drawn from 102
/// to 869 in 1024, so that no variant is of one allele alone, and each
/// allele of each genotype is drawn at that frequency. A genotype that is
/// missing has both of its alleles at [`MISSING_ALLELE`], which is what a
/// reader of popnei gives for one that was not called: a half called
/// genotype is missing among the genotypes of a population and gives its
/// called allele to the allele counts, and this benchmark draws none.
fn the_genotypes() -> Vec<i8> {
    let mut numbers = Numbers::seeded(SEED);
    let mut gts: Vec<i8> = Vec::new();
    // The variants times the individuals times the ploidy, which the blocks
    // below say they hold.
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
    gts
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
    /// [`NUM_BLOCKS`] blocks of those genotypes, each of them with the
    /// chromosome and the position of its variants, which the pass cuts the
    /// resampling groups from.
    ///
    /// The blocks are shared out evenly between the chromosomes of
    /// [`CHROM_NAMES`], the first half of them on the first one, and the
    /// variants of a chromosome sit [`BETWEEN_VARS`] base pairs apart from
    /// the position 1000 on: the layout of `big.vars`, so that a length of
    /// base pairs cuts the two settings into the same groups.
    fn of(gts: &[i8]) -> BlocksInMemory {
        let mut chroms = ChromTable::new();
        let numbers: Vec<u32> = CHROM_NAMES.iter().map(|name| chroms.intern(name)).collect();
        // Both are constants of this file and neither is 0.
        let per_chrom = NUM_BLOCKS.checked_div(NUM_CHROMS).unwrap_or(NUM_BLOCKS);
        let blocks: Vec<Block> = (0..NUM_BLOCKS)
            .map(|number| {
                let of_the_chrom = number.checked_div(per_chrom).unwrap_or(0);
                // A block of a chromosome the table does not hold would be
                // one of more blocks than the chromosomes take, which the
                // division above leaves out.
                let chrom = numbers.get(of_the_chrom).copied().unwrap_or(0);
                let first_var = number.checked_rem(per_chrom).unwrap_or(0);
                block_of(gts, chrom, first_var.saturating_mul(NUM_VARS_PER_BLOCK))
            })
            .collect();
        BlocksInMemory {
            blocks: blocks.into_iter(),
            individuals: names_of_the_individuals(),
            chroms,
        }
    }
}

/// One block of those genotypes on the chromosome `chrom`, whose first
/// variant is the one at the place `first_var` of that chromosome.
///
/// The genotypes are the same in every block, which the pass does not look
/// at, and the positions are not: they go up by [`BETWEEN_VARS`] from one
/// variant to the next, and a variant whose position went back is what the
/// walk that cuts the groups refuses.
fn block_of(gts: &[i8], chrom: u32, first_var: usize) -> Block {
    // The largest position is 50000000, which is the 50000 variants of a
    // chromosome 1000 base pairs apart, so nothing here saturates.
    let first = u64::try_from(first_var).unwrap_or(0);
    let pos: Vec<u64> = (0..NUM_VARS_PER_BLOCK)
        .map(|row| {
            first
                .saturating_add(u64::try_from(row).unwrap_or(0))
                .saturating_add(1)
                .saturating_mul(BETWEEN_VARS)
        })
        .collect();
    Block {
        num_vars: NUM_VARS_PER_BLOCK,
        num_individuals: NUM_INDIVIDUALS,
        ploidy: PLOIDY,
        gts: gts.to_vec(),
        chrom: Some(vec![chrom; NUM_VARS_PER_BLOCK]),
        pos: Some(pos),
        id: None,
        alleles: None,
        qual: None,
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
        // The blocks hold the genotypes, the chromosome and the position,
        // whatever is asked for: the pass asks for the first alone when it
        // was asked for no groups, and it reads no more of the block then.
    }

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        Vec::new()
    }
}

/// One run of one timing: how long it took and the line that says what it
/// did, which holds the pairs, the variants and the groups that say the call
/// did the work the setting names.
struct Run {
    took: Duration,
    did: String,
}

/// One pass over `reader` for `pops`, timed from the first block to the sums
/// of every pair.
///
/// `num_vars` is how many variants the setting gives, which is 100000 for
/// the blocks in memory and what the footer of the file says for the vars
/// file. The pairs, the variants and one measure are read after the clock
/// stops, which is also what keeps the sums from being dropped as unused.
fn one_calculation(
    reader: &mut impl BlockReader,
    pops: &Pops,
    options: &PopDistOptions,
    num_vars: u64,
) -> Result<Run, String> {
    let started = Instant::now();
    let sums = calc_pop_dist_sums(reader, pops, options).map_err(|problem| problem.to_string())?;
    let took = started.elapsed();
    let num_pairs = sums.num_pairs();
    let pairs_of_the_pops = num_pairs_of(pops.len());
    if num_pairs != pairs_of_the_pops {
        return Err(format!(
            "the call gave {num_pairs} pairs where the {num_pops} populations make \
             {pairs_of_the_pops}",
            num_pops = pops.len(),
        ));
    }
    if sums.num_vars() != num_vars {
        return Err(format!(
            "the call read {read} variants where this setting has {num_vars}",
            read = sums.num_vars(),
        ));
    }
    let fst = match sums.measure(PopDistMeasure::Fst, 0, 1) {
        Some(fst) => format!("{fst:.4}"),
        None => "no value".to_owned(),
    };
    let groups = match sums.groups().len() {
        0 => "with no resampling groups".to_owned(),
        num_groups => format!("in {num_groups} resampling groups"),
    };
    Ok(Run {
        took,
        did: format!(
            "{num_pairs} pairs over {num_vars} variants {groups}, F_ST of {first} and \
             {second} {fst}",
            first = pops.name(0),
            second = pops.name(1),
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

/// It runs one timing once without timing it and then `runs` times,
/// printing what each run took and then the best, the median and the worst
/// of the timed ones.
///
/// The run that is not timed is there for the page faults of the first
/// touch of the memory a calculation works in, which a process pays once
/// and which would otherwise all fall on whichever timing ran first.
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

/// What the timings of one populations file are over: the two settings, what
/// the pass is asked for and how it is run.
#[cfg(not(target_family = "wasm"))]
struct Settings<'a> {
    /// The genotypes of every block of the reader in memory.
    gts: &'a [i8],
    /// The bytes of the vars file, read from the disc.
    bytes: &'a [u8],
    /// How many variants that file holds, which its footer says.
    num_vars_of_the_file: u64,
    /// How much data a population needs at a variant and how the variants
    /// are cut into the resampling groups.
    options: &'a PopDistOptions,
    /// How many times each timing is run.
    runs: usize,
    /// The pool of one thread the second and the fourth timing run in.
    on_one_thread: &'a rayon::ThreadPool,
    /// How many threads the pool of the process has, which the first and the
    /// third timing run on.
    of_the_process: usize,
}

/// The four timings of one populations file: the two settings, each of them
/// on the threads of the process and then on one thread.
///
/// The populations are built twice, once over the individuals of each
/// setting, because a population is the indices of its individuals among
/// those of the reader it is read with. The two readers name the same 1000
/// individuals in the same order, so the two hold the same individuals.
#[cfg(not(target_family = "wasm"))]
fn the_four_timings(
    settings: &Settings,
    pops_in_memory: &Pops,
    pops_of_the_vars_file: &Pops,
) -> Result<(), String> {
    let of_the_process = settings.of_the_process;
    let in_memory = || {
        let mut reader = BlocksInMemory::of(settings.gts);
        let num_vars = u64::try_from(NUM_BLOCKS.saturating_mul(NUM_VARS_PER_BLOCK))
            .map_err(|_| "the variants of the blocks in memory are more than a u64 counts")?;
        one_calculation(&mut reader, pops_in_memory, settings.options, num_vars)
    };
    let of_the_file = || {
        let mut reader =
            VarsReader::new(Cursor::new(settings.bytes)).map_err(|problem| problem.to_string())?;
        one_calculation(
            &mut reader,
            pops_of_the_vars_file,
            settings.options,
            settings.num_vars_of_the_file,
        )
    };
    time_it(
        &format!("the blocks in memory, {of_the_process} threads"),
        settings.runs,
        in_memory,
    )
    .and_then(|()| {
        time_it("the blocks in memory, 1 thread", settings.runs, || {
            settings.on_one_thread.install(in_memory)
        })
    })
    .and_then(|()| {
        time_it(
            &format!("the vars file, {of_the_process} threads"),
            settings.runs,
            of_the_file,
        )
    })
    .and_then(|()| {
        time_it("the vars file, 1 thread", settings.runs, || {
            settings.on_one_thread.install(of_the_file)
        })
    })
}

/// How the groups of a run are said in the line that opens it.
fn groups_said(groups: JackknifeGroups) -> String {
    match groups {
        JackknifeGroups::None => "no resampling groups".to_owned(),
        JackknifeGroups::PerVariant => "one resampling group for each variant".to_owned(),
        JackknifeGroups::OfBasePairs(length) => {
            format!("resampling groups of {length} base pairs")
        }
    }
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
    // The reader is built once here for what the timings need to know
    // before the clock starts: the individuals the populations are looked
    // up among, and the variants of the file, which its footer says without
    // a batch being read.
    let (individuals_of_the_file, num_vars_of_the_file) = match VarsReader::new(Cursor::new(&bytes))
    {
        Ok(reader) => (
            reader.individuals().to_vec(),
            u64::try_from(reader.num_vars()).unwrap_or(u64::MAX),
        ),
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
    let options = PopDistOptions {
        min_num_individuals: DEFAULT_MIN_NUM_INDIVIDUALS,
        groups: arguments.groups,
    };
    let of_the_process = rayon::current_num_threads();
    let gts = the_genotypes();
    let individuals_in_memory = names_of_the_individuals();
    println!(
        "{blocks} blocks of {vars} variants of {individuals} individuals in memory on \
         {chroms} chromosomes, and {path}, {bytes} bytes of {file_vars} variants; {groups}, \
         a variant counts for a population of {min} called genotypes; {runs} runs of each \
         timing, on the {of_the_process} threads of the process and on one thread",
        blocks = NUM_BLOCKS,
        vars = NUM_VARS_PER_BLOCK,
        individuals = NUM_INDIVIDUALS,
        chroms = NUM_CHROMS,
        path = arguments.path.display(),
        bytes = bytes.len(),
        file_vars = num_vars_of_the_file,
        groups = groups_said(arguments.groups),
        min = DEFAULT_MIN_NUM_INDIVIDUALS,
        runs = arguments.runs,
    );
    let settings = Settings {
        gts: &gts,
        bytes: &bytes,
        num_vars_of_the_file,
        options: &options,
        runs: arguments.runs,
        on_one_thread: &on_one_thread,
        of_the_process,
    };
    for path in &arguments.pops {
        let timed = std::fs::read_to_string(path)
            .map_err(|problem| problem.to_string())
            .and_then(|text| pops_of_the_file(&text))
            .and_then(|named| {
                let pops_in_memory = Pops::from_names(&named, &individuals_in_memory)
                    .map_err(|problem| problem.to_string())?;
                let pops_of_the_vars_file = Pops::from_names(&named, &individuals_of_the_file)
                    .map_err(|problem| problem.to_string())?;
                println!(
                    "{path}: {num_pops} populations, {num_pairs} pairs",
                    path = path.display(),
                    num_pops = pops_in_memory.len(),
                    num_pairs = num_pairs_of(pops_in_memory.len()),
                );
                the_four_timings(&settings, &pops_in_memory, &pops_of_the_vars_file)
            });
        if let Err(problem) = timed {
            eprintln!("{path}: {problem}", path = path.display());
            return ExitCode::FAILURE;
        }
    }
    ExitCode::SUCCESS
}
