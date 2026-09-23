//! What each operation of this crate costs, called by itself on the shapes
//! popnei calls it at.
//!
//! The eight benchmarks of `crates/popnei/benches/` time a whole analysis
//! or a whole matrix, so two of them reach this crate through its callers
//! and neither can say which call got slower. This one calls the eleven
//! operations directly, builds every matrix it needs with the generator
//! the tests use, reads no file, and prints the time of one call of each.
//!
//! The numbers to compare against are the two tables of
//! `docs/specs/linalg.md`, "Speed" and "What the seven of the GWAS cost",
//! both taken on 22 and 23 September 2026 in trial crates that were never
//! in git: 12.7 ms for `add_self_product_lower` on a block of 5000 x 1000,
//! 0.035 s, 0.27 s and 6.3 s for `eigh_lower` at n = 1000, 2000 and 5000
//! on Accelerate, 0.0012 s, 0.0089 s and 0.135 s for `cholesky_lower` at
//! the same three, and 0.165 ms and 0.145 ms for `thin_qr` and `rank` on a
//! design of 10000 x 5. Nothing in the repository could produce any of
//! them again, so a later faer or a bump of the toolchain would move them
//! in silence.
//!
//! # The two backends
//!
//! The crate is the same code over BLAS and LAPACK, which on this machine
//! is Accelerate, and over faer, and which one a build links is the point
//! of the crate, so both are run:
//!
//! ```text
//! VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1 \
//!     cargo bench -p popnei-linalg --bench ops -- --runs 5
//! VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1 \
//!     cargo bench -p popnei-linalg --bench ops --no-default-features -- --runs 5
//! ```
//!
//! The first links BLAS and LAPACK, the second faer, and the header line
//! of a run says which. Accelerate takes the cores it finds unless
//! `VECLIB_MAXIMUM_THREADS` is set when the process starts, and faer runs
//! natively on the global pool of rayon, which `RAYON_NUM_THREADS` sizes,
//! so the two variables above are what asks each backend for one thread
//! and leaving them out is what asks for the threads of the machine. Both
//! go before the command and not into it: a variable set after the process
//! starts does not reach Accelerate. The header line prints what each of
//! them held, so a number taken with the wrong threads shows it.
//!
//! # What a region is, and why the numbers are per call
//!
//! Several of these calls are well under a millisecond, and timing one of
//! them at a time measures the harness: three measurements of `thin_qr` on
//! the same design of 10000 x 5, each the best of 20 single calls, spanned
//! 0.157 to 0.207 ms, 32 per 100, on a machine that was not doing anything
//! else. So a measurement times a region of many calls with one clock and
//! divides by the calls in it, and every line prints the time of one call.
//! How many calls a region holds is printed with each measurement: 1 for
//! the operations of some milliseconds, `add_self_product_lower`, the
//! product, the eigendecomposition, the Cholesky, the inverse and the
//! triangular solve over as many right hand sides as rows; 300 for the QR
//! and the rank of a design; 500 for the Cholesky solve over 10000 right
//! hand sides; 5000 for the log of the determinant; and 200000 for the
//! system of one variant and for the triangular solve of 5 coefficients.
//!
//! An operation that overwrites what it was given, the Cholesky and the
//! two solves, cannot be called twice over the same buffer: the second
//! call would factor the factorization, and a solve repeated 500 times
//! divides its right hand sides by the matrix 500 times and reaches
//! numbers no machine holds. Where a region is one call the buffer is
//! filled again before the clock starts. Where a region is many, the
//! copy is inside it, and the measurement is run a second time with the
//! call taken out and the copy left in, which is the line that says "the
//! copy alone"; the line after it is the difference, which is the call.
//!
//! # The matrices
//!
//! Every matrix is built here and no file is read, so a run repeats on any
//! machine. The numbers come from the xorshift generator of "How it is
//! verified" of `docs/specs/linalg.md`, started at 7, which is the
//! generator the tests of the crate and the trial crates of the two tables
//! used. A square matrix of n x n is `G = ZZ'` for a Z of n rows and 200
//! columns at n = 1000 or 20 columns at any other n, which is what the
//! second table of the spec measured on, and the ones that have to be
//! positive definite have n added to each entry of their diagonal. Only
//! the lower half of such a matrix is written, since that is the half
//! every operation which takes one reads. The block of the principal
//! component analysis, the two tiles of the matrix of r² and the design of
//! a study are the numbers of the generator as they come.
//!
//! # The command line
//!
//! `--runs` is 5 and `--sizes` is 1000, 2000 and 5000 when they are not
//! given, and `--ops` is every operation. A quick run is
//! `--sizes 1000`, and the Cholesky and the inverse of 10000 individuals,
//! which take 1.02 s and 2.05 s on Accelerate and 5.93 s and 11.8 s on
//! faer at one thread, are asked for by name:
//!
//! ```text
//! cargo bench -p popnei-linalg --bench ops -- \
//!     --ops cholesky_lower,invert_with_cholesky --sizes 10000
//! ```
//!
//! One region that is not timed comes before the timed ones, for the page
//! faults of the first touch of the memory an operation works in, which a
//! process pays once. Each measurement prints the time of one call in each
//! region, then the best, the median and the worst, and then the number
//! the last call gave, which is a trace, an eigenvalue or an entry of the
//! result: a run whose matrices came out as a shape of zeros does not
//! print the same number as one that computed something. The best is what
//! a report states, since every other process on the machine can only make
//! a region longer.

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

use popnei_linalg::{
    Result, TheFirstOperand, TheHalfThatHoldsTheMatrix, TheSecondOperand, add_self_product_lower,
    cholesky_lower, eigh_lower, invert_with_cholesky, log_determinant_with_cholesky, product, rank,
    solve_triangular, solve_with_cholesky, thin_qr,
};

/// How many times each measurement is timed when the command line does not
/// say.
const DEFAULT_RUNS: usize = 5;

/// The orders of the square matrices when the command line does not say,
/// which are the three the second table of `docs/specs/linalg.md`
/// measures at and which `docs/objectives.md` names as the individuals of
/// a study.
const DEFAULT_SIZES: [usize; 3] = [1000, 2000, 5000];

/// The largest order a square matrix may be asked for. It is the 10000
/// individuals of `docs/objectives.md`, and it also bounds every count of
/// values this file works out: the most is 10000 times 10000, which is
/// 100000000 and fits in a `usize` of this machine.
const THE_LARGEST_SIZE: usize = 10000;

/// The largest order the eigendecomposition is offered at. At 5000 it
/// takes 6.3 s on Accelerate and 11.4 s on faer, and 10000 is eight times
/// the arithmetic of that.
const THE_LARGEST_SIZE_OF_AN_EIGENDECOMPOSITION: usize = 5000;

/// The columns of the Z a square matrix of 1000 is built from, which is
/// what the second table of `docs/specs/linalg.md` used there.
const THE_COLUMNS_OF_Z_AT_1000: usize = 200;

/// The columns of the Z a square matrix of any other order is built from,
/// which is what that table used above 1000.
const THE_COLUMNS_OF_Z_ELSEWHERE: usize = 20;

/// The variants of the block of the principal component analysis, which is
/// the `rows` of the `a` of `add_self_product_lower`.
const THE_VARIANTS_OF_A_PCA_BLOCK: usize = 5000;

/// The individuals of that block, which is its `cols` and the order of the
/// `g` the products add up.
const THE_INDIVIDUALS_OF_A_PCA_BLOCK: usize = 1000;

/// The variants of one tile of the matrix of r², which is
/// `THE_VARS_OF_A_TILE` of `crates/popnei/src/ld.rs`.
const THE_VARIANTS_OF_A_TILE: usize = 1000;

/// The individuals a pair of tiles sums over, which is the dataset of
/// "Speed" of `docs/specs/ld.md`.
const THE_INDIVIDUALS_OF_A_TILE: usize = 1000;

/// The rows of the design a study fits its models on, which is the number
/// the second table of `docs/specs/linalg.md` measured the QR and the rank
/// at.
const THE_INDIVIDUALS_OF_A_DESIGN: usize = 10000;

/// The columns of that design: an intercept and four covariates.
const THE_COVARIATES_OF_A_DESIGN: usize = 5;

/// The orders of the one system per variant that a fit solves, which the
/// second table of `docs/specs/linalg.md` measured at 0.060, 0.173 and
/// 0.320 µs on Accelerate.
const THE_SIZES_OF_ONE_SYSTEM: [usize; 3] = [3, 7, 11];

/// How many calls a region holds for the log of the determinant, which
/// reads the n entries of a diagonal and computes nothing else.
const THE_CALLS_OF_A_REGION_OF_A_LOG_DETERMINANT: usize = 5000;

/// How many calls a region holds for the Cholesky solve over the right
/// hand sides of 10000 individuals, which the spec measured at 0.070 ms on
/// Accelerate.
const THE_CALLS_OF_A_REGION_OF_A_CHOLESKY_SOLVE: usize = 500;

/// How many calls a region holds for the operations of well under a
/// microsecond: the system of one variant and the triangular solve of a
/// handful of coefficients.
const THE_CALLS_OF_A_REGION_OF_A_SMALL_SOLVE: usize = 200_000;

/// How many calls a region holds for the QR and the rank of a design,
/// which the spec measured at 0.165 ms and 0.145 ms.
const THE_CALLS_OF_A_REGION_OF_A_DESIGN: usize = 300;

/// Which library the operations of this build run on, which is the whole
/// point of timing them twice.
const THE_BACKEND: &str = if cfg!(feature = "blas") {
    "BLAS and LAPACK, which on this machine is Accelerate"
} else {
    "faer"
};

/// One of the operations of the crate, as the command line names it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Operation {
    /// The lower half of `g += a'a` on a block of the principal component
    /// analysis.
    AddSelfProductLower,
    /// `c = a b'` on a pair of tiles of the matrix of r².
    Product,
    /// The eigendecomposition of a symmetric matrix.
    EighLower,
    /// The Cholesky factorization of a positive definite matrix.
    CholeskyLower,
    /// The lower half of the inverse, off that factorization.
    InvertWithCholesky,
    /// The log of the determinant, off that factorization.
    LogDeterminantWithCholesky,
    /// The solve against that factorization, with one right hand side for
    /// each of 10000 individuals.
    SolveWithCholesky,
    /// The factorization and the solve of one variant's system together,
    /// with one right hand side, which is what a fit runs once per
    /// variant.
    OneSystem,
    /// The solve against a triangular matrix, both halves.
    SolveTriangular,
    /// The thin QR of a design.
    ThinQr,
    /// The rank of a design.
    Rank,
}

impl Operation {
    /// Every operation, in the order a run that asks for all of them
    /// prints them.
    const ALL: [Operation; 11] = [
        Operation::AddSelfProductLower,
        Operation::Product,
        Operation::EighLower,
        Operation::CholeskyLower,
        Operation::InvertWithCholesky,
        Operation::LogDeterminantWithCholesky,
        Operation::SolveWithCholesky,
        Operation::OneSystem,
        Operation::SolveTriangular,
        Operation::ThinQr,
        Operation::Rank,
    ];

    /// The name the command line takes for it, which is the name of the
    /// function of the crate.
    fn name(self) -> &'static str {
        match self {
            Operation::AddSelfProductLower => "add_self_product_lower",
            Operation::Product => "product",
            Operation::EighLower => "eigh_lower",
            Operation::CholeskyLower => "cholesky_lower",
            Operation::InvertWithCholesky => "invert_with_cholesky",
            Operation::LogDeterminantWithCholesky => "log_determinant_with_cholesky",
            Operation::SolveWithCholesky => "solve_with_cholesky",
            Operation::OneSystem => "one_system",
            Operation::SolveTriangular => "solve_triangular",
            Operation::ThinQr => "thin_qr",
            Operation::Rank => "rank",
        }
    }

    /// The operation of that name, or nothing when no operation has it.
    fn of_the_name(name: &str) -> Option<Operation> {
        Operation::ALL
            .into_iter()
            .find(|operation| operation.name() == name)
    }
}

/// What the command line asked for.
struct Arguments {
    /// The operations to time, in the order of [`Operation::ALL`].
    operations: Vec<Operation>,
    /// The orders of the square matrices to time them at.
    sizes: Vec<usize>,
    /// How many timed regions each measurement gets.
    runs: usize,
}

/// What the benchmark does and what its command line takes, which is what
/// an argument it does not know and `--help` are answered with.
const USAGE: &str = "\
ops [--ops name,name] [--sizes n,n] [--runs n]

It times each operation of popnei-linalg by itself, on the shapes popnei
calls it at, over matrices it builds here with the xorshift generator of
docs/specs/linalg.md. No file is read.

  --ops name,name   which operations to time, all of them by default:
                    add_self_product_lower, product, eigh_lower,
                    cholesky_lower, invert_with_cholesky,
                    log_determinant_with_cholesky, solve_with_cholesky,
                    one_system, solve_triangular, thin_qr, rank
  --sizes n,n       the orders of the square matrices, 1000,2000,5000 by
                    default and 10000 at the most; the operations on a
                    block, a pair of tiles or a design have shapes of
                    their own and do not read it
  --runs n          how many timed regions each measurement gets, 5 by
                    default
  --help            this

A measurement times a region of one or more calls with one clock and
prints the time of one call. One region that is not timed comes first. The
number the last call gave is printed after the clock stops, so that a run
proves it computed something.

Which library the operations run on is the cargo feature `blas`: it is on
by default and links BLAS and LAPACK, and `--no-default-features` links
faer. Accelerate takes the cores it finds and faer runs on the global pool
of rayon, so one thread is asked for with VECLIB_MAXIMUM_THREADS=1 and
RAYON_NUM_THREADS=1 in the environment of the command, which Accelerate
reads when the process starts.";

/// What to run, or the message that says what the command line should have
/// been.
fn arguments() -> std::result::Result<Arguments, String> {
    let mut operations = Operation::ALL.to_vec();
    let mut sizes = DEFAULT_SIZES.to_vec();
    let mut runs = DEFAULT_RUNS;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--ops" => {
                let given = args
                    .next()
                    .ok_or_else(|| "--ops takes names and none came after it".to_owned())?;
                operations = the_operations_of(&given)?;
            }
            "--sizes" => {
                let given = args
                    .next()
                    .ok_or_else(|| "--sizes takes numbers and none came after it".to_owned())?;
                sizes = the_sizes_of(&given)?;
            }
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
            other => {
                return Err(format!(
                    "`{other}` is not an argument of this benchmark\n\n{USAGE}"
                ));
            }
        }
    }
    if runs == 0 {
        return Err("--runs is 1 or more".to_owned());
    }
    Ok(Arguments {
        operations,
        sizes,
        runs,
    })
}

/// The operations the comma separated `given` names, in the order of
/// [`Operation::ALL`], or the message that says which name is not one.
fn the_operations_of(given: &str) -> std::result::Result<Vec<Operation>, String> {
    let mut asked_for = Vec::new();
    for name in given.split(',').filter(|name| !name.is_empty()) {
        let operation = Operation::of_the_name(name)
            .ok_or_else(|| format!("`{name}` is not an operation of this crate\n\n{USAGE}"))?;
        asked_for.push(operation);
    }
    if asked_for.is_empty() {
        return Err("--ops names one operation at least".to_owned());
    }
    Ok(Operation::ALL
        .into_iter()
        .filter(|operation| asked_for.contains(operation))
        .collect())
}

/// The orders the comma separated `given` names, or the message that says
/// which one is not a number of 1 to [`THE_LARGEST_SIZE`].
fn the_sizes_of(given: &str) -> std::result::Result<Vec<usize>, String> {
    let mut sizes = Vec::new();
    for size in given.split(',').filter(|size| !size.is_empty()) {
        let size = size
            .parse::<usize>()
            .map_err(|_| format!("`{size}` is not a number of 1 to {THE_LARGEST_SIZE}"))?;
        if size == 0 || size > THE_LARGEST_SIZE {
            return Err(format!(
                "{size} is not an order this benchmark takes, which is 1 to {THE_LARGEST_SIZE}"
            ));
        }
        sizes.push(size);
    }
    if sizes.is_empty() {
        return Err("--sizes names one order at least".to_owned());
    }
    Ok(sizes)
}

/// The first `how_many` numbers of the generator of "How it is verified"
/// of `docs/specs/linalg.md`, which is the generator the tests of the
/// crate and the trial crates of its two tables of timings used: a 64 bit
/// state that starts at 7 with its lowest bit set, and for each number
/// `s ^= s << 13; s ^= s >> 7; s ^= s << 17`, the shifts dropping the bits
/// that leave the 64, and the number is `(s >> 11) / 2^53 - 0.5`.
fn the_numbers_of_the_generator(how_many: usize) -> Vec<f64> {
    let mut state = 7_u64;
    (0..how_many)
        .map(|_| {
            state ^= state.wrapping_shl(13);
            state ^= state.wrapping_shr(7);
            state ^= state.wrapping_shl(17);
            // 2^53, below which a count is an exact f64.
            state.wrapping_shr(11) as f64 / 9007199254740992.0 - 0.5
        })
        .collect()
}

/// How many columns the Z of a square matrix of `n` has: 200 at 1000 and
/// 20 at any other order, which is what "What the seven of the GWAS cost"
/// of `docs/specs/linalg.md` measured on.
fn the_columns_of_z(n: usize) -> usize {
    if n == 1000 {
        THE_COLUMNS_OF_Z_AT_1000
    } else {
        THE_COLUMNS_OF_Z_ELSEWHERE
    }
}

/// The lower half of `G = ZZ'` for a Z of `n` rows and [`the_columns_of_z`]
/// columns of the generator, with `n` added to each entry of its diagonal,
/// which makes it positive definite whatever the rank of `ZZ'` is. The
/// upper half is 0 and no operation that takes such a matrix reads it.
///
/// # Errors
///
/// What `add_self_product_lower` gives, which is nothing for these shapes.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "n is at most THE_LARGEST_SIZE, 10000, which the reading of the command line \
              refuses above, so n times n is at most 100000000 and the columns of Z times n \
              at most 2000000, and both fit in a usize of every target popnei builds for"
)]
fn a_positive_definite_matrix(n: usize) -> Result<Vec<f64>> {
    let columns = the_columns_of_z(n);
    // ZZ' is A'A for A = Z', which is the columns of Z by n, and the
    // numbers in the order the generator gives them are the rows of that
    // A. This is how the tests of the crate build the same matrix.
    let z_written_the_other_way_round = the_numbers_of_the_generator(columns * n);
    let mut g = vec![0.0_f64; n * n];
    add_self_product_lower(&z_written_the_other_way_round, columns, n, &mut g)?;
    for (row, values) in g.chunks_exact_mut(n).enumerate() {
        if let Some(entry) = values.get_mut(row) {
            *entry += n as f64;
        }
    }
    Ok(g)
}

/// The Cholesky factorization of [`a_positive_definite_matrix`] of `n`,
/// which is what the solve, the inverse and the log of the determinant are
/// given.
///
/// # Errors
///
/// What `cholesky_lower` gives, which is nothing for a matrix with `n`
/// added to its diagonal.
fn a_factorization(n: usize) -> Result<Vec<f64>> {
    let mut l = a_positive_definite_matrix(n)?;
    cholesky_lower(&mut l, n)?;
    Ok(l)
}

/// The sum of the diagonal of the `n` x `n` matrix in `values`, which is
/// the number most of the measurements print to say that they computed
/// something.
fn the_trace(values: &[f64], n: usize) -> f64 {
    values
        .chunks_exact(n)
        .enumerate()
        .filter_map(|(row, row_of_values)| row_of_values.get(row).copied())
        .sum()
}

/// The first value of `values`, or NaN when it holds none.
fn the_first(values: &[f64]) -> f64 {
    values.first().copied().unwrap_or(f64::NAN)
}

/// The time each region took and the number the last call of the last
/// region gave.
struct Timed {
    /// One time for each region, the first being the region that is not
    /// counted.
    regions: Vec<Duration>,
    /// What the last call gave, printed after the clock stops.
    gave: f64,
}

/// Runs `runs + 1` regions of `calls` calls each and gives the time of
/// every one of them, the first being the region that is not counted.
///
/// `prepare` runs before each region with the clock stopped, which is
/// where a buffer that one call overwrites is filled again; `call` is the
/// one call, and it is what the clock is around. Both take the buffers as
/// an argument, so that neither has to borrow them from the other.
///
/// # Errors
///
/// What the operation `call` makes gives.
fn time_the_regions<Buffers>(
    runs: usize,
    calls: usize,
    buffers: &mut Buffers,
    mut prepare: impl FnMut(&mut Buffers),
    mut call: impl FnMut(&mut Buffers) -> Result<f64>,
) -> Result<Timed> {
    let mut gave = f64::NAN;
    let mut regions = Vec::with_capacity(runs.saturating_add(1));
    for _ in 0..=runs {
        prepare(buffers);
        let started = Instant::now();
        for _ in 0..calls {
            gave = black_box(call(black_box(buffers))?);
        }
        regions.push(started.elapsed());
    }
    Ok(Timed { regions, gave })
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
    // The two are a time each and `lower` is the shorter of them, so half
    // of what lies between them added to it is their middle, and neither
    // the subtraction nor the addition leaves what a `Duration` holds.
    let between = upper.saturating_sub(lower);
    lower.saturating_add(between.checked_div(2).unwrap_or(between))
}

/// The seconds one call took, which is the region divided by the calls in
/// it.
fn seconds_of_a_call(region: Duration, calls: usize) -> f64 {
    region.as_secs_f64() / calls as f64
}

/// A time of one call written with four digits, in microseconds below a
/// millisecond, in milliseconds below a second and in seconds above.
fn a_call(seconds: f64) -> String {
    if seconds < 1e-3 {
        return format!("{microseconds:.4} µs", microseconds = seconds * 1e6);
    }
    if seconds < 1.0 {
        return format!("{milliseconds:.4} ms", milliseconds = seconds * 1e3);
    }
    format!("{seconds:.4} s")
}

/// Prints what a measurement took, one line for each region and then the
/// best, the median and the worst of the timed ones with the number the
/// last call gave. Gives back the seconds of the best call, which is what
/// a floor is taken off.
fn print_what_it_took(timed: &Timed, calls: usize) -> f64 {
    let Some((not_timed, counted)) = timed.regions.split_first() else {
        println!("  no region was run");
        return f64::NAN;
    };
    println!(
        "  the region that is not timed: {took} a call",
        took = a_call(seconds_of_a_call(*not_timed, calls))
    );
    for (number, region) in counted.iter().enumerate() {
        println!(
            "  region {number}: {took} a call",
            number = number.saturating_add(1),
            took = a_call(seconds_of_a_call(*region, calls))
        );
    }
    let best = seconds_of_a_call(sorted_time(counted, 0), calls);
    println!(
        "  best {best}, median {median}, worst {worst} a call; it gave {gave:.10e}",
        best = a_call(best),
        median = a_call(seconds_of_a_call(median_time(counted), calls)),
        worst = a_call(seconds_of_a_call(
            sorted_time(counted, counted.len().saturating_sub(1)),
            calls
        )),
        gave = timed.gave,
    );
    best
}

/// Prints the header of a measurement: what is timed, on what shape, and
/// how many calls a region holds.
fn print_what_is_timed(what: &str, calls: usize) {
    println!("{what}, {calls} call(s) a region");
}

/// Prints the floor of a measurement whose region copies a buffer before
/// every call, and the difference, which is the call without the copy.
fn print_the_floor(best: f64, floor: f64) {
    println!(
        "  the copy alone, which the region pays before every call: {floor} a call, \
         so the call without it: {without}",
        floor = a_call(floor),
        without = a_call((best - floor).max(0.0)),
    );
}

/// The self product of a block of the principal component analysis:
/// `g += a'a` for an `a` of 5000 variants x 1000 individuals.
///
/// `g` is filled with 0 before each region, since the operation adds to
/// what it holds and a region that ran many of them would grow it without
/// end.
///
/// # Errors
///
/// What `add_self_product_lower` gives.
fn the_self_product_of_a_block(runs: usize) -> Result<()> {
    let rows = THE_VARIANTS_OF_A_PCA_BLOCK;
    let cols = THE_INDIVIDUALS_OF_A_PCA_BLOCK;
    print_what_is_timed(
        &format!("add_self_product_lower, a of {rows} x {cols} into g of {cols} x {cols}"),
        1,
    );
    let a = the_numbers_of_the_generator(the_values_of(rows, cols));
    let mut g = vec![0.0_f64; the_values_of(cols, cols)];
    let timed = time_the_regions(
        runs,
        1,
        &mut g,
        |g| g.fill(0.0),
        |g| {
            add_self_product_lower(black_box(&a), rows, cols, black_box(g))?;
            Ok(the_trace(g, cols))
        },
    )?;
    print_what_it_took(&timed, 1);
    Ok(())
}

/// The product of a pair of tiles of the matrix of r²: `c = a b'` for two
/// sets of 1000 variants over 1000 individuals, which is what
/// `crates/popnei/src/ld.rs` calls for each of the six sums of a pair.
///
/// # Errors
///
/// What `product` gives.
fn the_product_of_two_tiles(runs: usize) -> Result<()> {
    let vars = THE_VARIANTS_OF_A_TILE;
    let individuals = THE_INDIVIDUALS_OF_A_TILE;
    print_what_is_timed(
        &format!(
            "product, c = a b' for an a of {vars} x {individuals} and a b of {vars} x \
             {individuals}"
        ),
        1,
    );
    let both = the_numbers_of_the_generator(the_values_of(the_values_of(vars, individuals), 2));
    let (a, b) = both.split_at(the_values_of(vars, individuals));
    let mut c = vec![0.0_f64; the_values_of(vars, vars)];
    let timed = time_the_regions(
        runs,
        1,
        &mut c,
        |_| {},
        |c| {
            product(
                TheFirstOperand::ByTheRowsOfTheResult {
                    values: black_box(a),
                    rows: vars,
                },
                individuals,
                TheSecondOperand::ByTheColumnsOfTheResult {
                    values: black_box(b),
                    cols: vars,
                },
                black_box(c),
            )?;
            Ok(the_trace(c, vars))
        },
    )?;
    print_what_it_took(&timed, 1);
    Ok(())
}

/// The eigendecomposition of a symmetric matrix of `n` x `n`, which the
/// principal component analysis runs once on the individuals of the study.
///
/// The matrix is taken by value and comes back as the eigenvectors, so the
/// buffer is filled again before each region, with the clock stopped, and
/// the call takes it out of the buffers and leaves an empty one behind,
/// which copies nothing.
///
/// # Errors
///
/// What `eigh_lower` gives.
fn the_eigendecomposition(n: usize, runs: usize) -> Result<()> {
    print_what_is_timed(&format!("eigh_lower, g of {n} x {n}"), 1);
    let g = a_positive_definite_matrix(n)?;
    let mut working = g.clone();
    let timed = time_the_regions(
        runs,
        1,
        &mut working,
        |working| {
            working.clear();
            working.extend_from_slice(&g);
        },
        |working| {
            let eigen = eigh_lower(black_box(std::mem::take(working)), n)?;
            Ok(black_box(&eigen)
                .values
                .first()
                .copied()
                .unwrap_or(f64::NAN))
        },
    )?;
    print_what_it_took(&timed, 1);
    Ok(())
}

/// The Cholesky factorization of a positive definite matrix of `n` x `n`,
/// which is the first of the two calls the fit of a mixed model runs a few
/// dozen times.
///
/// It overwrites the lower half of the matrix it is given, so the buffer is
/// filled again before each region, with the clock stopped, and a region is
/// one call.
///
/// # Errors
///
/// What `cholesky_lower` gives.
fn the_factorization(n: usize, runs: usize) -> Result<()> {
    print_what_is_timed(&format!("cholesky_lower, a of {n} x {n}"), 1);
    let a = a_positive_definite_matrix(n)?;
    let mut working = a.clone();
    let timed = time_the_regions(
        runs,
        1,
        &mut working,
        |working| working.copy_from_slice(&a),
        |working| {
            cholesky_lower(black_box(working), n)?;
            Ok(the_trace(working, n))
        },
    )?;
    print_what_it_took(&timed, 1);
    Ok(())
}

/// The lower half of the inverse of the matrix a factorization of `n` x
/// `n` came off, which is the second of those two calls.
///
/// The factorization is read and not written and the inverse goes into a
/// buffer of its own, so nothing has to be filled again.
///
/// # Errors
///
/// What `invert_with_cholesky` gives.
fn the_inverse(n: usize, runs: usize) -> Result<()> {
    print_what_is_timed(&format!("invert_with_cholesky, l of {n} x {n}"), 1);
    let l = a_factorization(n)?;
    let mut inverse = vec![0.0_f64; the_values_of(n, n)];
    let timed = time_the_regions(
        runs,
        1,
        &mut inverse,
        |_| {},
        |inverse| {
            invert_with_cholesky(black_box(&l), n, black_box(inverse))?;
            Ok(the_trace(inverse, n))
        },
    )?;
    print_what_it_took(&timed, 1);
    Ok(())
}

/// The log of the determinant off a factorization of `n` x `n`, which
/// reads the n entries of its diagonal and computes nothing else.
///
/// # Errors
///
/// What `log_determinant_with_cholesky` gives.
fn the_log_determinant(n: usize, runs: usize) -> Result<()> {
    let calls = THE_CALLS_OF_A_REGION_OF_A_LOG_DETERMINANT;
    print_what_is_timed(
        &format!("log_determinant_with_cholesky, l of {n} x {n}"),
        calls,
    );
    let mut l = a_factorization(n)?;
    let timed = time_the_regions(
        runs,
        calls,
        &mut l,
        |_| {},
        |l| log_determinant_with_cholesky(black_box(l), n),
    )?;
    print_what_it_took(&timed, calls);
    Ok(())
}

/// The solve against a factorization of 5 x 5 with one right hand side for
/// each of 10000 individuals, which is what a study solves once for a whole
/// design and what the spec measured at 0.070 ms on Accelerate.
///
/// The right hand sides come back as the solutions, so the region copies
/// them in again before each call and a second measurement times that copy
/// by itself.
///
/// # Errors
///
/// What `solve_with_cholesky` gives.
fn the_cholesky_solve(runs: usize) -> Result<()> {
    let n = THE_COVARIATES_OF_A_DESIGN;
    let sides = THE_INDIVIDUALS_OF_A_DESIGN;
    let calls = THE_CALLS_OF_A_REGION_OF_A_CHOLESKY_SOLVE;
    print_what_is_timed(
        &format!("solve_with_cholesky, l of {n} x {n}, b of {sides} x {n}"),
        calls,
    );
    let l = a_factorization(n)?;
    let of_the_right_hand_sides = the_numbers_of_the_generator(the_values_of(sides, n));
    let mut b = of_the_right_hand_sides.clone();
    let timed = time_the_regions(
        runs,
        calls,
        &mut b,
        |_| {},
        |b| {
            b.copy_from_slice(&of_the_right_hand_sides);
            solve_with_cholesky(black_box(&l), n, black_box(b), sides)?;
            Ok(the_first(b))
        },
    )?;
    let best = print_what_it_took(&timed, calls);
    let floor = time_the_regions(
        runs,
        calls,
        &mut b,
        |_| {},
        |b| {
            b.copy_from_slice(&of_the_right_hand_sides);
            Ok(the_first(black_box(b)))
        },
    )?;
    print_the_floor(
        best,
        seconds_of_a_call(sorted_time(the_counted_regions(&floor), 0), calls),
    );
    Ok(())
}

/// The system of one variant: the factorization of an `n` x `n` and the
/// solve against it with one right hand side, which is what a fit runs
/// once for each variant of a block and what the spec measured at 0.060,
/// 0.173 and 0.320 µs on Accelerate for n of 3, 7 and 11.
///
/// Both the matrix and the right hand side come back written over, and a
/// caller has to build them for each variant anyway, so the region copies
/// both in before each call and a second measurement times those copies by
/// themselves.
///
/// # Errors
///
/// What `cholesky_lower` and `solve_with_cholesky` give.
fn the_system_of_one_variant(n: usize, runs: usize) -> Result<()> {
    let calls = THE_CALLS_OF_A_REGION_OF_A_SMALL_SOLVE;
    print_what_is_timed(
        &format!(
            "one_system, cholesky_lower and solve_with_cholesky, a of {n} x {n}, one right hand side"
        ),
        calls,
    );
    let a = a_positive_definite_matrix(n)?;
    let of_the_right_hand_side = the_numbers_of_the_generator(n);
    let mut buffers = (a.clone(), of_the_right_hand_side.clone());
    let timed = time_the_regions(
        runs,
        calls,
        &mut buffers,
        |_| {},
        |(working, b)| {
            working.copy_from_slice(&a);
            b.copy_from_slice(&of_the_right_hand_side);
            cholesky_lower(black_box(working), n)?;
            solve_with_cholesky(black_box(working), n, black_box(b), 1)?;
            Ok(the_first(b))
        },
    )?;
    let best = print_what_it_took(&timed, calls);
    let floor = time_the_regions(
        runs,
        calls,
        &mut buffers,
        |_| {},
        |(working, b)| {
            working.copy_from_slice(&a);
            b.copy_from_slice(&of_the_right_hand_side);
            Ok(the_first(black_box(b)))
        },
    )?;
    print_the_floor(
        best,
        seconds_of_a_call(sorted_time(the_counted_regions(&floor), 0), calls),
    );
    Ok(())
}

/// The solve against the lower half of a factorization of `n` x `n` with
/// one right hand side for each of `n` individuals, which is what each
/// step of the fit of the null model of the logistic mixed model asks for.
///
/// # Errors
///
/// What `solve_triangular` gives.
fn the_triangular_solve_of_a_study(n: usize, runs: usize) -> Result<()> {
    print_what_is_timed(
        &format!("solve_triangular, the lower half, a of {n} x {n}, b of {n} x {n}"),
        1,
    );
    let l = a_factorization(n)?;
    let of_the_right_hand_sides = the_numbers_of_the_generator(the_values_of(n, n));
    let mut b = of_the_right_hand_sides.clone();
    let timed = time_the_regions(
        runs,
        1,
        &mut b,
        |b| b.copy_from_slice(&of_the_right_hand_sides),
        |b| {
            solve_triangular(
                black_box(&l),
                n,
                TheHalfThatHoldsTheMatrix::TheLowerHalf,
                black_box(b),
                n,
            )?;
            Ok(the_first(b))
        },
    )?;
    print_what_it_took(&timed, 1);
    Ok(())
}

/// The solve against the upper half of the `r` of the thin QR of a design
/// of 10000 x 5 with one right hand side, which is the second half of
/// fitting a linear model and what line 378 of pyNei's `gwas.py` asks for.
///
/// # Errors
///
/// What `thin_qr` and `solve_triangular` give.
fn the_triangular_solve_of_a_fit(runs: usize) -> Result<()> {
    let n = THE_COVARIATES_OF_A_DESIGN;
    let calls = THE_CALLS_OF_A_REGION_OF_A_SMALL_SOLVE;
    print_what_is_timed(
        &format!("solve_triangular, the upper half, a of {n} x {n}, one right hand side"),
        calls,
    );
    let design = the_numbers_of_the_generator(the_values_of(THE_INDIVIDUALS_OF_A_DESIGN, n));
    let factorization = thin_qr(&design, THE_INDIVIDUALS_OF_A_DESIGN, n)?;
    let of_the_right_hand_side = the_numbers_of_the_generator(n);
    let mut b = of_the_right_hand_side.clone();
    let timed = time_the_regions(
        runs,
        calls,
        &mut b,
        |_| {},
        |b| {
            b.copy_from_slice(&of_the_right_hand_side);
            solve_triangular(
                black_box(&factorization.r),
                n,
                TheHalfThatHoldsTheMatrix::TheUpperHalf,
                black_box(b),
                1,
            )?;
            Ok(the_first(b))
        },
    )?;
    let best = print_what_it_took(&timed, calls);
    let floor = time_the_regions(
        runs,
        calls,
        &mut b,
        |_| {},
        |b| {
            b.copy_from_slice(&of_the_right_hand_side);
            Ok(the_first(black_box(b)))
        },
    )?;
    print_the_floor(
        best,
        seconds_of_a_call(sorted_time(the_counted_regions(&floor), 0), calls),
    );
    Ok(())
}

/// The thin QR of a design of 10000 x 5, which a study makes once and
/// which the spec measured at 0.165 ms on Accelerate.
///
/// The design is read and not written, and each call allocates the `q` and
/// the `r` it gives, which is what a caller pays too.
///
/// # Errors
///
/// What `thin_qr` gives.
fn the_qr_of_a_design(runs: usize) -> Result<()> {
    let rows = THE_INDIVIDUALS_OF_A_DESIGN;
    let cols = THE_COVARIATES_OF_A_DESIGN;
    let calls = THE_CALLS_OF_A_REGION_OF_A_DESIGN;
    print_what_is_timed(&format!("thin_qr, a of {rows} x {cols}"), calls);
    let mut design = the_numbers_of_the_generator(the_values_of(rows, cols));
    let timed = time_the_regions(
        runs,
        calls,
        &mut design,
        |_| {},
        |design| {
            let factorization = thin_qr(black_box(design), rows, cols)?;
            Ok(the_first(&black_box(&factorization).r))
        },
    )?;
    print_what_it_took(&timed, calls);
    Ok(())
}

/// The rank of a design of 10000 x 5, which a study takes once to refuse a
/// design whose covariates repeat each other, and which the spec measured
/// at 0.145 ms on Accelerate.
///
/// # Errors
///
/// What `rank` gives.
fn the_rank_of_a_design(runs: usize) -> Result<()> {
    let rows = THE_INDIVIDUALS_OF_A_DESIGN;
    let cols = THE_COVARIATES_OF_A_DESIGN;
    let calls = THE_CALLS_OF_A_REGION_OF_A_DESIGN;
    print_what_is_timed(&format!("rank, a of {rows} x {cols}"), calls);
    let mut design = the_numbers_of_the_generator(the_values_of(rows, cols));
    let timed = time_the_regions(
        runs,
        calls,
        &mut design,
        |_| {},
        |design| Ok(rank(black_box(design), rows, cols)? as f64),
    )?;
    print_what_it_took(&timed, calls);
    Ok(())
}

/// The regions of a measurement that count, which is every one but the
/// first.
fn the_counted_regions(timed: &Timed) -> &[Duration] {
    timed
        .regions
        .split_first()
        .map_or(&[][..], |(_, counted)| counted)
}

/// How many values a matrix of `rows` x `cols` holds.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "every shape of this benchmark comes from a constant of this file or from a size \
              the command line refused above THE_LARGEST_SIZE, 10000, so the largest product \
              taken here is the 2 x 1000 x 1000 of the pair of tiles and every one fits in a \
              usize of this machine"
)]
fn the_values_of(rows: usize, cols: usize) -> usize {
    rows * cols
}

/// Times one operation at the orders the command line asked for, and
/// prints a line for each measurement.
///
/// # Errors
///
/// What the operation gives.
fn time_the_operation(operation: Operation, sizes: &[usize], runs: usize) -> Result<()> {
    match operation {
        Operation::AddSelfProductLower => the_self_product_of_a_block(runs),
        Operation::Product => the_product_of_two_tiles(runs),
        Operation::EighLower => {
            for size in sizes {
                if *size > THE_LARGEST_SIZE_OF_AN_EIGENDECOMPOSITION {
                    println!(
                        "eigh_lower at {size}: left out, since the largest order this benchmark \
                         offers it at is {THE_LARGEST_SIZE_OF_AN_EIGENDECOMPOSITION}"
                    );
                    continue;
                }
                the_eigendecomposition(*size, runs)?;
            }
            Ok(())
        }
        Operation::CholeskyLower => {
            for size in sizes {
                the_factorization(*size, runs)?;
            }
            Ok(())
        }
        Operation::InvertWithCholesky => {
            for size in sizes {
                the_inverse(*size, runs)?;
            }
            Ok(())
        }
        Operation::LogDeterminantWithCholesky => {
            for size in sizes {
                the_log_determinant(*size, runs)?;
            }
            Ok(())
        }
        Operation::SolveWithCholesky => the_cholesky_solve(runs),
        Operation::OneSystem => {
            for size in THE_SIZES_OF_ONE_SYSTEM {
                the_system_of_one_variant(size, runs)?;
            }
            Ok(())
        }
        Operation::SolveTriangular => {
            the_triangular_solve_of_a_fit(runs)?;
            for size in sizes {
                the_triangular_solve_of_a_study(*size, runs)?;
            }
            Ok(())
        }
        Operation::ThinQr => the_qr_of_a_design(runs),
        Operation::Rank => the_rank_of_a_design(runs),
    }
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

/// The benchmark runs on the BLAS of the system or on faer's threads and
/// reads a clock that wasm has not. This is what `cargo wasm-check`
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
    println!(
        "popnei-linalg on {THE_BACKEND}, {runs} runs, orders {sizes}, \
         VECLIB_MAXIMUM_THREADS {veclib}, RAYON_NUM_THREADS {rayon}",
        runs = arguments.runs,
        sizes = arguments
            .sizes
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(", "),
        veclib = said_about_the_variable("VECLIB_MAXIMUM_THREADS"),
        rayon = said_about_the_variable("RAYON_NUM_THREADS"),
    );
    for operation in arguments.operations {
        if let Err(problem) = time_the_operation(operation, &arguments.sizes, arguments.runs) {
            eprintln!("{name}: {problem}", name = operation.name());
            return ExitCode::FAILURE;
        }
    }
    ExitCode::SUCCESS
}
