# Performance review: the linear algebra crate and its two backends

23 September 2026. `crates/popnei-linalg` gained eleven operations on 23
September 2026, when the plan `docs/plans/linalg-gwas.md` merged at
`75c5fd3`, and nothing calls any of them. This review measures what the
crate costs on the three operations that popnei does reach, says whether
the crate needs a benchmark of its own, and answers the question two
earlier reviews deferred: whether the caller should be allowed to vouch
that it has already checked its values. It was asked for by the owner.

The words this document uses. The crate holds fourteen **operations** of
linear algebra behind one interface, each a Rust function over matrices
held as `&[f64]`, row after row. Under that interface are two
**backends**, two libraries that do the arithmetic: the BLAS and LAPACK of
the system, which on this machine is **Accelerate**, the library of macOS
that numpy also calls; and **faer**, a linear algebra library written in
Rust, which is what runs in a browser tab, where there is no Accelerate,
and natively when the crate is built with its `blas` feature off. A
**block** is the run of variants a reader gives at a time, 5000 variants
of 1000 individuals in every measurement here, and a **tile** is the
square block of variants the matrix of r² is worked out in, 1000 variants
today. **`dsyrk`** is Accelerate's routine for the product of a matrix
with its own transpose, **`dgemm`** for the product of two matrices, and
**`dsyevd`** for the eigendecomposition of a symmetric matrix, which is
the operation that gives the axes a principal component analysis projects
onto. The **scan** is the walk this crate makes over every value of every
matrix it is given, before it hands it to a backend, refusing an infinity
or a NaN; "Errors" of `docs/specs/linalg.md` requires it, because the two
backends do not treat a NaN alike. The **promise** is the change two
earlier reviews refused: a type that says the caller has already checked
its values, so that the scan can be skipped. A **unit in the last place**
is the smallest difference two neighbouring `f64` can have, and is the
unit the two backends' disagreements are stated in.

## 1. The scope and its limits

`crates/popnei-linalg` at `75c5fd3`: `src/lib.rs`, the interface and the
checks, `src/blas.rs`, `src/faer.rs` and `Cargo.toml`. Seven categories
were sent, one reviewer each with a fresh context: methodology, numbers,
allocations, data layout, concurrency, hot loops, and linear algebra with
wasm. Two were not sent and why: the crate opens no file and makes no
syscall of its own, so input and output had nothing to review; and the
crate has no binding of its own and every timing here is a cargo
benchmark, which is the condition under which the Python boundary
checklist says to skip it.

**The scope splits in two by call frequency, and every finding below says
which half it is about.** Three operations are reached from popnei today.
`add_self_product_lower`, the product of a matrix with itself, is called
once per block by the principal component analysis, at
`crates/popnei/src/pca.rs:223` and `:875`. `product` is called three
times from `pca.rs` and once from `crates/popnei/src/ld.rs:784`, which is
six times per pair of tiles.
`eigh_lower`, the eigendecomposition, is called once per analysis. The
other eleven — a Cholesky factorization and the solve, the log of the
determinant and the inverse that come off it, the thin QR of a design,
the solve against a triangular matrix of either half, the rank, and the
two products that the merged plan's first work package added, which
compute the product with its first matrix turned — are reached by nothing
at all: `grep` for them across `crates/popnei`, `crates/popnei-python`
and `crates/popnei-js` finds no call. Their caller will be the genome wide
association study, whose spec is being written on the branch `spec/gwas`
and is not merged.

**A profile cannot rank code nothing calls, so no finding about those
eleven is above Likely.** That is the limit of this review, and it is not
one that more work would remove: it goes away when the association study
is written and not before. The numbers those eleven have are of the
backend routines alone, in "What the seven of the GWAS cost" of
`docs/specs/linalg.md`, and they were taken from a trial crate,
`tmp/linalg_gwas_trial/`, which was never in git and is no longer on
disc. Nothing in the repository can produce any of them again.

Two sampling profiles were taken for this review, with `/usr/bin/sample`
over 18 s on the running benchmark, on the owner's Apple M5 Pro, rustc
1.98, the `bench` profile, `VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1`,
machine otherwise idle at a load average of 0.96 on 18 cores. Idle threads
are taken out of both totals.

| | on-CPU samples | the backend | this crate's scan | this crate, all else |
|---|---|---|---|---|
| the analysis, no weights | 12132 | 59.5% | 2.9% | 0.05% |
| the matrix of r², 5000 variants | 11276 | 78.4% | 3.6% | 0 |

So on both hot paths the crate's own code is 3 to 4 per 100 of the time
and the library underneath is 60 to 78. The rest of the analysis is the
decompression of the file, 21%, and the standardizing of the rows, 8%.

## 2. The verdict: run the experiments, and two decisions are the owner's

There are candidates and their order is in section 3. None is an "apply":
the largest thing this crate's own code costs on a measured path is the
scan, at 2.9% and 3.6% of on-CPU samples, and the three reviewers who
looked at it disagree by three times about whether a cheaper loop would
recover any of it. That disagreement is what the first experiment
settles.

Two things are the owner's and are in section 5 as such: whether to take
the promise, which is a change to "Errors" of `docs/specs/linalg.md`; and
what to do about the eigendecomposition on faer giving a different result
at each number of threads, which is a departure from popnei's own rule and
was found by measurement in this review.

## 3. The measurement plan

In the order in which they unblock the findings.

1. **Is the scan bound by its instructions or by memory?** A binary that
   scans one buffer of 32 KB, 1 MB, 8 MB and 40 MB in the loop as it is
   and in each of the two cheaper forms of H1, best of five each. The
   deterministic gate first: `cargo asm --lib -p popnei-linalg` on
   `every_value_is_finite`, 20 instructions per 64 bytes now. If the rate
   at 40 MB is the same for all three forms, the loop is at the memory
   roof at popnei's sizes and H1 is closed without touching either
   benchmark. If the cheaper forms are faster at 40 MB, the two
   benchmarks follow.
2. **Do the two benchmarks see it?** `pca_vars` at `--num-prin-comps 0`
   and `r2_matrix` at `--max-num-vars 5000`, both one thread, best of
   five, against the baselines of section 1. Keep only if the r² falls by
   more than 0.006 s from 0.386 s, which is half of what the scan was
   measured to cost there.
3. **What the crate costs over its backend, which nothing measures
   today.** A benchmark in `crates/popnei-linalg/benches/`, built like the
   eight in `crates/popnei/benches/`: `harness = false`, one untimed warm
   run, best, median and worst, a checksum printed after the clock stops,
   matrices made in the benchmark so that no file is needed. It covers
   `add_self_product_lower` at 5000 x 1000, `product` at the tile shapes
   of the r², and `eigh_lower` at n = 1000 and 5000, and it is run twice,
   once per backend. This is the deliverable the owner asked about and
   section 5's L1 argues for it.
4. **The faer backend through the callers' benchmarks, which nobody has
   run.** `cargo bench --no-default-features --features bench-internals`
   builds both `pca_vars` and `r2_matrix` on faer natively; the second
   feature is the one those two benchmarks already need, which opens the
   private functions they time. It is the
   browser's arithmetic timed without a browser, and it is one flag away.
5. **Peak memory in a tab**, for the eigendecomposition, which needs three
   n x n matrices alive at once: the analysis under node at n = 2000,
   3000, 4000 and 5000, recording the heap. The metric is peak bytes, not
   time. `docs/objectives.md` names 10000 individuals, where three such
   matrices are 2.4 GB in an address space of at most 4 GB.

What cannot be measured here: whether faer is level with OpenBLAS on
x86, which is the open question of `docs/rust_core.md`, since this machine
has no OpenBLAS.

## 4. The build configuration

`Cargo.toml` at the workspace root sets `overflow-checks = false` and
`debug = "line-tables-only"` under `[profile.release]`, the same debug
setting under `[profile.bench]`, and `debug = true` under
`[profile.profiling]`. `lto`, `codegen-units` and `opt-level` are unset,
there is no global allocator, and `.cargo/config.toml` sets rustflags for
the two wasm targets only.

Three of the four usual experiments are not worth running here, each for a
measured reason. Link time optimization and one codegen unit were already
tried by `docs/reports/perf-stats-2026-09-22.md`, which got 0.182 and
0.183 s against 0.183 and 0.184 s at 2.5 times the rebuild; and on this
crate's native path 60 to 78 per 100 of the time is inside
`libBLAS.dylib`, a dynamic library no link time optimization reaches.
`target-cpu` moves little on aarch64, where NEON is already in the
baseline target, and a distributed wheel cannot use `native`. A different
allocator is not indicated: this crate takes a handful of multi-megabyte
buffers per call, which go to `mmap`, and the stats review's profile put
the allocator at 6 samples of 5966.

The one case that has never been tried is link time optimization on the
**faer** backend, which is Rust all the way down, across faer, its traits
crate and the crate that does its products, and which is the whole of the
arithmetic in a browser and with `--no-default-features`. It is worth one
experiment, after the benchmark of plan item 3 exists, so that the faer
product can be timed on its own rather than inside a whole analysis.

**The toolchain is not pinned.** There is no `rust-toolchain.toml`, and
`rust-version = "1.98"` in `[workspace.package]` is a floor and not a pin.
The one hot Rust loop of this crate is vectorized by the compiler and not
by the source: it is written as eight scalar counters that rustc 1.98
turns into four chains of vector registers, verified in the assembly by
three of the reviewers on both backends and on `wasm32-unknown-unknown`.
A toolchain bump can undo that with every test still passing and no
warning. Pinning costs one file and makes a deliberate bump invalidate the
saved baselines, which `profiling_environment.md` already says it does.

## 5. The findings

### For the owner

**O1. The promise that the caller has already checked: what it is worth
now.** `crates/popnei-linalg/src/lib.rs:245-246` and `:420-421`. The
reached three and the eleven both.

Two reviews refused this change, H3 of `docs/reports/perf-pca-2026-09-22.md`
on 22 September and L1 of `docs/reports/perf-ld-2026-09-23.md` on 23
September, each on the ground that it needs a type carrying the promise,
which changes "Errors" of `docs/specs/linalg.md` and reaches `pca.rs`
too. L1 left it "worth revisiting only for both modules together, if the
linear algebra crate ever changes that interface for another reason", and
that condition has now been met twice over: `product` gained a typed first
operand and the triangular solve gained an argument for which half holds
the matrix.

What it is worth, in numbers:

| where | values scanned now | with the promise | what it saves |
|---|---|---|---|
| the matrix of r², 5000 variants | 160 million | 15 million | 0.017 s of 0.386 s, 4% |
| the analysis, per run | 110 million | 110 million | nothing |

The r² gains because the three matrices of a tile are built once and used
in every pair that tile takes part in, so the same buffer is scanned four
to six times per pair. **The analysis gains nothing from this crate's
side**, because every block is a new matrix: the caller would run the same
scan one line earlier. The analysis's share only comes if `pca.rs` also
proves each block finite as it standardizes it, from a table of at most
`ploidy + 1` entries, which is the half of H3 that lives in that module
and not in this crate.

Where it never pays: the scan is 3 to 4 per 100 only because Accelerate
runs these products at about 530 GFLOP/s. The same block product is 10.5
ms on Accelerate, 95 ms on faer and 187 ms in wasm, so the same scan is
0.6 per 100 on faer and under 1 per 100 in a browser. For the eleven the
routines are O(n³) over an O(n²) matrix, so the scan before the inverse at
n = 5000 is 1.3 ms against 0.205 s, 0.6 per 100.

The options.

- **Take it now, for all fourteen signatures.** Costs a public type whose
  constructor runs the scan, one line and a `?` at each of the four call
  sites in the core crate, a sentence of `docs/specs/linalg.md`, and the
  invariant that a buffer is not changed between the check and the call.
  Gains 4 per 100 of the r² and nothing of the analysis. What makes now
  the cheap moment is the eleven: they have interfaces and no callers, so
  putting the type in all fourteen signatures before `docs/specs/gwas.md`
  is written against them costs far less than changing them afterwards.
- **Take it only when the association study is written**, and pay to
  change fourteen signatures that then have callers.
- **Leave it**, and keep paying 4 per 100 of the r², 0.6 per 100 of faer
  and under 1 per 100 in a browser.

Recommended: leave it, which is the same answer the two earlier reviews
gave and for the same reason, that 4 per 100 of one calculation does not
justify a change to the spec and a permanent invariant. What the merged
plan changed is the price of the edit and not the size of the gain. **But
the moment to decide is before `docs/specs/gwas.md` settles**: the eleven
signatures have no callers today, so if the promise is ever wanted, now is
when it is cheap, and after that spec is written against them it is not.
The thing that would change the answer is a measurement showing the scan
costs more on a path this review could not reach, which is the whole
second group.

**O2. On faer the eigendecomposition gives a different result at each
number of threads.** `crates/popnei-linalg/src/faer.rs:394`. The reached
three. Measured in this review, and confirmed by the orchestrator on a
second run.

The coding skill's rule is that no result of popnei depends on the number
of threads. This one does. Measured on G = A'A for A of 4000 x 800 from
the spec's own generator, release, `--no-default-features`, the largest
eigenvalue of the 800:

| threads | the largest eigenvalue | units in the last place from one thread |
|---|---|---|
| 1 | 693.254303030172 | 0 |
| 4 | 693.2543030301708 | −10 |
| 8 | 693.2543030301703 | −15 |
| 18 | 693.2543030301722 | +2 |

Two runs at one thread gave identical bits, so it is the pool size and not
run to run noise. Over all 800 eigenvalues the gap reaches 6.2e-15
relative and the eigenvector entries 1.0e-12 absolute. The products are
bit identical at 1, 4, 8 and 18 threads, and Accelerate is bit identical
at 1, 4 and all cores for all three operations.

The mechanism is in the code and its own comment says so. The nine
operations that take a thread argument are given `the_threads()`, which
is the pool of rayon natively and one thread in a browser;
`self_adjoint_eigen` at `faer.rs:394` takes no such argument and reads
faer's global parallelism instead, which rayon sizes. Its divide and
conquer then joins the same sums in a different order at each pool size.

What it does and does not touch. In a browser there is one thread, so
nothing moves there. Natively with the `blas` feature on, Accelerate runs
and nothing moves. The exposure is a native build with
`--no-default-features`, where the analysis's projections change in the
last bits with `RAYON_NUM_THREADS`. The gaps sit inside the tolerance
"How it is verified" of `docs/specs/linalg.md` states, 1e-12 relative on
eigenvalues and 1e-9 on entries, and every test passes; but they are the
same size as the gap between the two backends that the spec records, 2.9e-14
and 1.3e-12, so a thread change spends the same budget a second time.

**What pinning costs was measured for this report** and it is the number
the decision turns on. Pinning gives the same bits at every pool size,
confirmed at 1, 4 and 18 threads on all 1000 eigenvalues and all 1000000
eigenvector entries of an n = 1000, compared as equal; and it costs:

| n | faer on the pool | pinned to one thread | what pinning costs |
|---|---|---|---|
| 1000 | 0.068 s | 0.097 s | 1.43 times, 0.029 s |
| 2000 | 0.289 s | 0.716 s | 2.48 times, 0.43 s |
| 5000 | 2.635 s | 11.744 s | 4.46 times, 9.1 s |
| 10000 | 25.41 s | 147.39 s | 5.80 times, 122 s |

The options.

- **Pin it**, by calling faer's lower level entry point that takes a
  thread argument. Costs the table above, 9.1 s at 5000 individuals and
  two minutes at the 10000 of `docs/objectives.md`, and 32 lines in place
  of 7 against a faer API likelier to move between versions than the one
  used now. Makes the rule true again.
- **Record it**, adding a sentence to "Threads" of `docs/specs/linalg.md`
  saying that this one operation on this one backend is not reproducible
  across pool sizes, and a test that would catch the day it changes.
  Costs a paragraph.
- **Leave it silent**, since it is inside the tolerance and invisible in a
  browser and on Accelerate.

Recommended: record it, and do not pin it. Two minutes of an analysis at
10000 individuals is too much to pay for 15 units in the last place that
sit well inside the tolerance the spec already states, and the two places
popnei's users actually are, a browser and a native build on Accelerate,
are both unaffected. This is a correctness matter that a performance
review happened to find, so the sentence it needs may belong to a code
review rather than to this report.

**The same experiment found that the spec's faer column is a one-thread
column.** The pinned times above, 0.097, 0.716 and 11.744 s, reproduce
"Speed" of `docs/specs/linalg.md`, which records 0.096, 0.68 and 11.4 s
for faer. On the pool faer takes 2.635 s at n = 5000 against the 6.3 s
that table gives Accelerate, so **faer's eigendecomposition natively is
2.4 times faster than Accelerate's at that size, not 1.8 times slower as
the spec reads**. The table is not wrong, it is a one-thread table beside
an Accelerate column that takes the threads it finds, and nothing on its
face says so. Whoever next edits that section should say which column is
which.

### Hot-path

**H1. The scan costs 2.9 and 3.6 per 100, and three reviewers disagree by
three times about whether a cheaper loop recovers any of it.**
`crates/popnei-linalg/src/lib.rs:1082-1095`. The reached three. Confidence
medium.

Six of the seven reviewers named this site. It is the only function of
this crate that either profile names: 347 of 12132 on-CPU samples in the
analysis and 408 of 11276 in the r². `docs/reports/perf-ld-2026-09-23.md`
measured it directly by removing it, 0.383 s against 0.366 s at a tile of
1000.

The loop compiles to 20 instructions per 64 bytes, four accumulators,
verified in the assembly on both backends and on wasm. Two cheaper forms
were proposed, both giving the same answer for every input:

- ORing the comparison's all-ones mask straight into the accumulator,
  instead of `u64::from(bool)` which forces it to 0 or 1. Drops 4 of the
  20 instructions, and 8 of the 48 vector operations of the wasm build.
  Costs one `.wrapping_neg()` and a sentence saying the counters hold a
  mask.
- Accumulating `value * 0.0`, which is NaN exactly when the value was not
  finite. Drops the count to about 8 per 64 bytes, needs more accumulators
  to hide the add's latency, and rests on float semantics rather than bit
  patterns.

Whether either buys time depends on what the loop is bound by, and the
three estimates of its throughput do not agree: 26 GB/s from the analysis
block scan of 22 September, 42 GB/s from what H3 left behind, and 75 GB/s
from the r² measurement. At 26 GB/s the loop is at one core's memory roof
and fewer instructions buy nothing; at 75 GB/s it is bound by its own
instructions and they buy up to half.

Measurement plan: item 1 of section 3, which measures the loop alone at
four buffer sizes and answers this directly, before either benchmark is
run. Gate on the instruction count from `cargo asm` first.

Effect on the numbers: none. The predicate accepts and refuses exactly the
same matrices, and the two tests added at `e64637c`, which put a value
that is not finite in every place of a 2 x 12 and a 12 x 12, cover the
loop and its tail.

**Measured, no gain, closed.** Section 8 has the numbers.

### Likely

**L1. No benchmark times one call of this crate, and the eleven have
targets nothing can reproduce.** `crates/popnei-linalg/Cargo.toml` has no
`[[bench]]`. Both halves. Confidence high.

The two benchmarks that reach this crate time a whole analysis and a whole
matrix, so they can say that popnei got slower and cannot say which call
did it. The eigendecomposition cannot be separated at all: the phases of
it that are matrix products land in `libBLAS.dylib`, where the LAPACK
total is 2.3% of samples
while the spec's 0.035 s at n = 1000 is 7.7% of a 0.453 s run, and 60 to
78 per 100 of the samples are unsymbolicated addresses inside that
library, so no profile goes finer. The claims that steer this crate — 12.7
ms per block against `dsyrk`'s 12.05 ms, and the checks at 1.68 ms of
12.05 ms — were measured in scratch crates that no longer exist, so no
change to the checks or the allocations can be re-gated against them.
Section 3 item 3 is the plan; the answer to the owner's question is yes,
and what it should measure is there.

**L2. Three measurements of the same call span 32 per 100, and the spread
is the harness rather than the machine.** `docs/reports/linalg-gwas.md`
records `thin_qr` on one design of 10000 x 5 at 0.157, 0.1855 and 0.207
ms and calls it the spread between runs on a busy machine. Each of the
three is a best of 20, and a minimum over 20 runs is the statistic that
transient load cannot inflate, so the gap is between processes and
systematic. Three candidates, all testable: Accelerate takes the threads
it finds and nothing pinned `VECLIB_MAXIMUM_THREADS` for those runs, so a
thread rendezvous of tens of microseconds is in or out of a 0.16 ms call;
this machine keeps a short-lived process on one cluster of cores for all
20 runs; and each number came from a separately built trial crate.
Everything in this crate that is not a whole-dataset operation is measured
at that scale, so the benchmark of section 3 item 3 times an inner loop of
many calls as one region rather than one call at a time, and wraps its
inputs so that nothing is folded away.

**L3. Two buffers are filled with zeros and then wholly overwritten.**
`crates/popnei-linalg/src/blas.rs:801` and `:329`. The eleven, and the
reached eigendecomposition. Confidence high on the mechanism.

The column major copy that `thin_qr` and `rank` make calls
`resize(a.len(), 0.0)` and then writes every one of those values, 400 KB
of dead stores per call at a design of 10000 x 5. The workspace of
`dsyevd` is filled with zeros the routine does not need, 16 MB per call at
n = 1000, which is the `Vec<f64>::resize` the analysis profile names at 6
samples. Counted, not timed: `thin_qr` at 10000 x 5 does one allocation of
400000 bytes that the crate memsets and five that arrive already zeroed.
The fix for the copy is `extend` after the existing `try_reserve_exact`,
which keeps `Error::Memory` and adds no `unsafe`; the fix for the
workspace needs `alloc_zeroed` behind an `unsafe` block this crate
confines to its backend calls, and is not worth it.

**L4. The two transposing loops move one value per 9 instructions where 6
would do.** `crates/popnei-linalg/src/blas.rs:802-806` and `:827-833`.
The eleven. Confidence medium.

`into.iter_mut().zip(a.iter().skip(column).step_by(cols))` gives the
compiler no bound on the strided side, so each element pays an extra
compare and branch; the assembly is 9 instructions and 2 branches per
value. The same loop written over `chunks_exact(cols)` with a `get` gives
6 instructions and 1 branch, with the same loads, stores and strides. No
benchmark reaches either function, so the gate is the instruction count
until the association study brings a caller.

**L5. `thin_qr` and `rank` allocate everything per call and no caller can
reuse a buffer.** `crates/popnei-linalg/src/lib.rs:763-788`. The eleven.
Confidence medium. Counted per call at 10000 x 5: `thin_qr` 6 allocations
and 802800 bytes on Accelerate, 20 and 1386233 on faer; `rank` 4 and
480320, and 11 and 533628. A study that factors a design per variant pays
all of it a million times, and the interface gives no way to hoist it.
**This is decided while `docs/specs/gwas.md` is written and not after**,
because the shape it would need — buffers passed in, or a workspace type
the caller holds — is cheapest to add before there are callers. It is the
same timing argument as O1.

**L6. faer forks the pool for a solve whose matrix is tiny and whose right
hand sides are many.** `crates/popnei-linalg/src/faer.rs:44`, as used at
`:276` and `:465`. The eleven. Confidence medium. faer halves the right
hand sides recursively and hands each half to the pool when the sides are
above 64 and the matrix is at most 128, which is exactly the shape
"What the seven of the GWAS cost" measured at n = 5 with one right hand
side for each of 10000 individuals, 0.070 ms on Accelerate against 0.140
ms on faer. About five levels of fork-join over leaves that are a few
hundred nanoseconds of arithmetic. The split condition does not read the
thread argument, so one thread and the pool recurse the same way and the
result is bit identical either way, which makes this a pure speed
question. Gate on `n` times `sides` in `the_threads`, after a trial crate
finds where the two cross.

**L7. The eigendecomposition needs three n x n matrices alive at once, and
that is the browser's ceiling before time is.**
`crates/popnei-linalg/src/faer.rs:393-417`. The reached three. Confidence
medium, from faer's source and not from a run. faer allocates the
eigenvectors and a scratch buffer while popnei's own matrix stays alive:
24 MB at n = 1000 and at least 2.4 GB at the 10000 individuals
`docs/objectives.md` names, in a tab with at most 4 GB that also holds the
standardized block. `docs/rust_core.md` section 3.2 already shows the
knee, 2.5 times native at n = 1000 and 2000 and 4.6 times at 3000. Two
things follow: popnei has no number for the largest n a tab does, which
section 3 item 5 would give it; and faer's allocation ends the process
where the spec promises `Error::Memory`, which is a matter for a code
review.

### Notes

- **The column major copy still pays at every size the objectives name.**
  `blas.rs:788`. The owner asked whether it does. It is the textbook naive
  transpose, worth about 16 times the matrix in line traffic once the
  matrix leaves the last level cache, but a design of 10000 x 5 is 400 KB
  and stays in L2, which is why the route through it measured 0.165 ms
  against 2.98 ms without it. Nearly square is the better case for the
  naive loop, not the worse one. The shape where blocking would win is
  about a million rows by five columns, and `docs/objectives.md` caps
  individuals at 10000, so no design of the association study reaches it.
  No change; only L4 improves it.
- **A partial eigendecomposition is worth 22 per 100, not a multiple.**
  `docs/specs/linalg.md` measured the routine that computes a chosen range
  at 0.025, 0.18 and 4.9 s against 0.035, 0.27 and 6.3 s, because the
  reduction to tridiagonal form is the bulk. faer 0.24.4 has no
  eigensolver over a range, so it would help Accelerate and do nothing in
  a browser, which is where the analysis is slowest. Recorded so that
  Open 1 of `docs/specs/pca.md`, which asks whether the analysis should
  give every component or only the first few, is decided with that number.
- **The lower half of the accumulated product is scanned once per block,
  and after the first block those values are the crate's own output**, not
  what the caller gave. 500500 values per block, about a ninth of what the
  scan reads and 0.26 per 100 of the run. "Errors" of the spec says the
  crate refuses what it is given and not what it produced; this is the one
  place it does the second. Too small to act on alone, and it belongs to
  whatever O1 becomes.
- **The early return for a block with no rows sits after the scan**, at
  `lib.rs:247`, so a block whose variants all lacked variance pays the
  scan for a call that adds nothing. Moving it up would stop refusing a
  matrix that is not finite on that call, which is a spec question.
- **`the_threads()` is computed for three faer requests that discard it**,
  at `faer.rs:224`, `:275` and `:307`. A few nanoseconds twice per call,
  named only because the association study's smallest calls are 60 to 320
  nanoseconds whole and there are 5000 of them per block per iteration.

## 6. Seen outside the scope

- `crates/popnei/src/ld.rs:784-795` gives the same two buffers to
  `product` four to six times per pair of tiles, and each call rescans
  both. That is the caller half of O1 and where the r²'s 4 per 100 lives.
- `crates/popnei/src/pca.rs:672` asks for a product that
  `docs/specs/linalg.md` measured at 0.102 ms against 0.030 ms for the
  transposed form on the same shapes, and the caller then re-scatters the
  result column by column. The combination that reads both of its
  operands the other way round, which that work package added and nothing
  calls, may give both the faster routine and a contiguous copy. This is
  the owner's second open decision at the end of
  `docs/reports/linalg-gwas.md`, which asks whether that product should
  replace a buffer and a loop in the analysis, and this review supports
  filing it.
- The analysis profile puts all 2546 samples of decompression, 21 per 100
  of on-CPU, on the same thread as the 6939 of the library. The read ahead
  thread of section 3 of `docs/architecture.md` is not on this path.
- faer's eigendecomposition allocates through a call that ends the process
  when the machine has not the memory, where the same operation on
  Accelerate gives `Error::Memory`. The difference is invisible to the
  caller and no test reaches it. For a code review.

## 7. What the code already does well

- **The row major to column major turn is done with flags, not copies.**
  `blas.rs:1-38` and the five product routines. Every matrix popnei holds
  is the transpose of what a BLAS routine reads, and the module handles it
  by choosing the routine, the half and the transpose flags so that none
  of the four products, the Cholesky solve or the triangular solve copies
  anything. The two operations that do copy say in their doc comments what
  the copy buys, with both numbers.
- **The error messages are built only when there is an error.**
  `lib.rs:236-449`. Every `format!` and `to_owned` sits inside an
  `ok_or_else` or a `map_err` closure, so a call that succeeds allocates
  nothing; counted at zero allocations for both hot operations.
- **The benchmarks that exist state their environment and prove their own
  work.** `crates/popnei/benches/pca_vars.rs` and `r2_matrix.rs` warm once
  untimed, print best, median and worst, name the environment variables in
  the doc comment because a variable set after the process starts would
  not reach Accelerate, and print a checksum taken after the clock stops.

## 8. What the experiments showed

### H1, the cheaper scan: measured, no gain, closed

The loop is not at the memory roof, and a cheaper one still buys nothing
in either benchmark. Both halves of that are worth keeping, because they
point in opposite directions and only the second decides.

In isolation, a scratch binary scanning one buffer that is already
resident, best of five, each trial moving at least 512 MB:

| buffer | as it is | the mask form | the float form, 16 counters |
|---|---|---|---|
| 32 KB | 60.5 GB/s | 93.9 | 145.5 |
| 1 MB | 72.5 GB/s | 96.3 | 133.1 |
| 8 MB | 71.7 GB/s | 95.9 | 119.1 |
| 40 MB | 71.9 GB/s | 90.4 | 104.2 |

The instruction counts per 64 bytes, which are the same on every run: 20
as it is, 16 for the mask form, 11.5 for the float form at 16 counters.
The float form is 12 and not the 8 that was expected, because without
`mul_add`, which is a software call under emscripten, the multiply and the
add are two instructions. Thirty two counters add nothing over sixteen.
So the loop reaches 72 GB/s at 8 MB and above, which matches the 75 GB/s
derived from the r² and refutes the 26 GB/s derived from the analysis:
**the loop is bound by its own instructions, not by memory.**

The float form at 16 counters was then put in the crate. Every check
passed before any timing: 149 tests on Accelerate, 136 on faer, 604 in the
workspace, clippy clean, both wasm targets checked, and the assembly
confirmed 11.5 instructions per 64 bytes in the crate as built. The two
benchmarks were run alternately, both sides in the same window, because
the drift between quiet windows, 0.015 to 0.020 s, is larger than the
0.006 s the decision turned on:

| | as it is | with the cheaper loop |
|---|---|---|
| the matrix of r², four paired rounds, best | 0.388, 0.387, 0.388, 0.390 s | 0.391, 0.390, 0.391, 0.390 s |
| the analysis, three paired rounds, best | 0.449, 0.447, 0.447 s | 0.446, 0.444, 0.446 s |

**The r² rose in every one of the four pairs.** The analysis fell by 0.002
to 0.003 s, 0.6 per 100, inside the spread of its own five runs. The
threshold was a fall of more than 0.006 s in the r². Nothing was kept and
the file is as it was.

Why a loop that is 3 per 100 of the profile and 1.45 times cheaper gives
nothing back: the scan is the first pass over a matrix the library has
just written, so much of what the profiler charges to it is the misses of
that first touch. A cheaper loop does not remove those misses, it hands
them to the next reader of the same buffer. The scratch binary re-reads a
buffer that is already resident, which is why it sees the instruction
count and the benchmark does not.

What is left of it. At 32 KB the float form is 2.4 times the throughput of
the one in place and at 1 MB 1.8 times, so if the eleven operations, which
work on far smaller matrices, ever show this function in a profile of
their own, the experiment is worth repeating there with a benchmark of
that size. The mask form is a one-word change worth about half the float
form in isolation, and was not carried to the benchmarks once the faster
form failed the threshold.

### O2, what pinning the eigendecomposition costs: measured, for the owner

The table is in O2 above. Pinning removes the thread dependence
completely, confirmed as bit equality at 1, 4 and 18 threads, and costs
1.43 times the wall time at 1000 individuals, 4.46 at 5000 and 5.80 at
10000. The recommendation in O2 is to record the divergence and not pay
that.

### What the experiments did not settle

L1's benchmark is being built and its result goes here when it lands.
Nothing else in section 5 was run: L3, L4, L6 and L7 are all about the
eleven operations, and none of them can be timed until either that
benchmark or the association study gives them a caller. That is not a gap
this review can close by working longer.
