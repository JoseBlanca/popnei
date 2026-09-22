# Performance review: the Kosman distances between individuals

22 September 2026. The review of `calc_kosman_sums` of the core, the
calculation that the plan `docs/plans/dists-kosman.md` built and that
its measurement, `docs/reports/dists-kosman-measurement.md`, found short
of the three numbers of "Speed" of `docs/specs/dists.md`. It says where
the time goes, which was not known, what to change in which order, and,
at the end of each finding, what its experiment gave. It was run on the
branch `perf/dists-kosman` by the session of the assistant, with one
reviewer subagent per category and one subagent per measurement, as the
`performance-review` skill says. The next review starts from its
measurement plan and its numbers.

The words this document uses. The **calculation** is one call of
`calc_kosman_sums` over a reader, which gives, for every pair of
individuals, two integers, the ploidy times the sum of the Kosman
distances and the number of variants called in both. A **block** is the
array of genotypes a reader gives at a time, 5000 variants of 1000
individuals in the dataset below, so 20 blocks. The calculation has two
phases per block: **the sets phase**, `KosmanBits::of_block`, which turns
the genotypes of the block into sets of bits, one word of 64 bits per 64
variants, one `called` set per individual and one `holds` set per allele
and per count from 1 to the ploidy, 5 sets of 79 words per individual
here; and **the pairs phase**, which for each of the 499500 pairs ANDs
the sets of the two individuals word by word and counts the ones, then
adds the two integers of the pair into the accumulator. **The reading
alone** is a pass over the file with nothing done to the blocks, and
**the target** is one of the three numbers of "Speed", which are of the
calculation with the reading taken out. A **gate** is what decides an
experiment: a count that is the same on every run where there is one,
and the wall time of the benchmark otherwise.

## 1. The scope and its limits

Reviewed: `crates/popnei/src/dists.rs` at 59e76a5, `main` at 02a2e46
plus the benchmark of this review, on the owner's Apple M5 Pro, 6
performance cores and 12 efficiency cores, 64 GB, macOS 27.0, rustc
1.98.0, node 26.8.2. The dataset is that of "Speed": 100000 variants x
1000 individuals, biallelic, 3 in 100 genotypes missing, `big.vars` at
`/Users/jose/devel/popnei-bench/`, and the same block made at random in
memory. The targets, with the reading taken out: 0.97 s on one thread,
0.38 s on 18 cores, 1.43 s in wasm under node; the measurement of the
plan gave 1.154 s, 0.625 s and 2.130 s. pyNei, the comparison and not a
target, takes 0.719 s on one thread and 0.280 s on 6.

Left out: the readers, whose cost is outside the targets and was
measured on its own, 0.101 to 0.110 s of the whole call; the binding
crates and the Python package, reviewed for the boundary alone; a read
ahead thread, an item of `docs/specs/block.md`; any dataset but this one
and any ploidy but 2. Categories sent: methodology, numbers,
allocations, data_layout, concurrency, hot_loops, linalg_and_wasm and
python_boundary; io_and_syscalls was not, since the reader is out.

The hot path evidence, taken for this review at 59e76a5 and quoted in
section 3: a benchmark, a sampling profile on one thread and on 18, and
the split of the time between the reader, the sets phase and the pairs
phase. The one thing that could not be profiled is the wasm run, whose
build strips the names of the functions; the plan below says how.

## 2. The verdict

Run the experiments. The profile names the sites, the reviewers agree on
them from every side, and each has a gate that can be counted before the
wall time is looked at:

- On 18 cores the sets phase is 0.569 s of the 0.624 s and runs on the
  calling thread while 17 workers wait; the pairs phase scales 10.8
  times. The target's arithmetic leaves 0.325 s for the sets phase.
  Building the sets on the threads, one range of individuals per work
  item, changes no bit.
- On every target the sets phase pays, per allele, the construction and
  the drop of an `Error` that is never returned, a tenth of the CPU, and
  per genotype a call to `memchr` over two bytes, another tenth; both
  are in the wasm assembly too. Two one line changes, each with the
  call gone from the assembly as its gate.
- On one thread and in wasm the pairs phase is 0.593 s of 1.123 s, at
  the SIMD issue width and not memory bound: 3 NEON instructions per
  word pair from a 3.16 MB working set that sits in L2. A biallelic
  diploid block needs 3 sets per individual and not 5, since the sets of
  the second allele are the complement of the first's within `called`,
  which cuts the words a pair walks from 790 to 474 and stays exact.

## 3. The measurements this review stands on

The benchmark, `crates/popnei/benches/kosman_dists.rs`, added at
59e76a5, `harness = false` as the others, run as

    cargo bench --bench kosman_dists -- /Users/jose/devel/popnei-bench/big.vars --runs 5

times the calculation four ways, best of 5 after one run that is not
timed: over a reader written in the bench that gives the same block 20
times from memory, and over `big.vars` read into memory, each on the
rayon pool of the process, 18 threads, and inside a pool of one thread
built in the bench. It asserts that the call gave 499500 pairs. Load
average 2.9 to 4.1 for the first invocation, and a second one agreed to
2 in 100:

| setting | best of 5 |
|---|---|
| the blocks in memory, 18 threads | 0.627 s |
| the blocks in memory, 1 thread | 1.138 s |
| the vars file, 18 threads | 0.748 s |
| the vars file, 1 thread | 1.261 s |

The in memory reader takes 0.121 to 0.123 s off the vars file numbers,
which is the reading alone measured by the plan, 0.101 to 0.110 s, to
0.02 s. The pool of one thread built in the bench and
`RAYON_NUM_THREADS=1` agree to 0.2 in 100.

The split of the phases, by an `Instant` around the three calls of
`calc_kosman_sums`, a temporary change that was not committed, best of
3, load average 1.3:

| setting | whole | the reader | the sets phase | the pairs phase |
|---|---|---|---|---|
| in memory, 18 threads | 0.624 s | 0.000 s | 0.569 s | 0.055 s |
| in memory, 1 thread | 1.123 s | 0.000 s | 0.529 s | 0.593 s |
| the vars file, 18 threads | 0.742 s | 0.110 s | 0.577 s | 0.055 s |
| the vars file, 1 thread | 1.240 s | 0.101 s | 0.546 s | 0.593 s |

The sampling profile, `sample <pid> 8` on the bench binary in its vars
file section, built with `[profile.bench] debug = "line-tables-only"`,
which 59e76a5 added because the release profile has no debug info and
the frames had no names without it. `sample` prints inclusive counts, so
the self time below is each frame's count less its children's, with the
idle frames of the rayon pool taken out. On one thread, 6714 samples,
the shares are shares of wall time; on 18 threads they are shares of
CPU time over all the threads, and the pairs phase, half of the CPU, is
7 in 100 of the wall.

One thread:

    3219   47.9%  add_the_pairs_of_the_block::{closure#1}::call_mut
    1388   20.7%  <popnei::dists::KosmanBits>::of_block
     694   10.3%  core::ptr::drop_glue::<popnei::error::Error>
     632    9.4%  core::slice::memchr::memchr
     491    7.3%  <lz4_flex::frame::decompress::FrameDecoder<&[u8]> as std::io::Read>::read_to_end
     212    3.2%  popnei::dists::alleles_of
      50    0.7%  _platform_memmove
      13    0.2%  __bzero

18 threads, 10622 samples, 1.91 cores busy on average:

    5293   49.8%  add_the_pairs_of_the_block::{closure#1}::call_mut
    2055   19.3%  <popnei::dists::KosmanBits>::of_block
    1075   10.1%  core::ptr::drop_glue::<popnei::error::Error>
     867    8.2%  core::slice::memchr::memchr
     733    6.9%  <lz4_flex::frame::decompress::FrameDecoder<&[u8]> as std::io::Read>::read_to_end
     327    3.1%  popnei::dists::alleles_of
      73    0.7%  _platform_memmove
      44    0.4%  __psynch_mutexdrop
      20    0.2%  <crossbeam_deque::deque::Stealer<rayon_core::job::JobRef>>::steal

The main thread's samples by phase: on one thread, the sets phase 43.6
in 100, the pairs phase 48.0, all of it a wait on a condition variable
while the one worker of the pool counts the pairs, the reader 8.4; on 18
threads, the sets phase 78.0, the pairs phase 7.0, the reader 15.0. All
18 workers were seen working in the pairs phase, 5 to 6 in 100 of the
window each. The three frames inside the sets phase, from the assembly
of `of_block` (`cargo asm -p popnei --lib "of_block"`): the drop glue is
the `.ok_or(Error::AlleleBelowTheMissingOne { allele })?` of the loop
over the alleles of a genotype, which writes the error to the stack and
calls its out of line drop on the success path, 2e8 times; `memchr` is
`genotype.contains(&MISSING_ALLELE)` on a slice of two bytes, 1e8 times;
`alleles_of` is the pass that finds the smallest and the largest allele
and ranks them. The same two calls are in the wasm32 assembly of
`of_block`. The pairs loop, `sums_of_two`, is vectorized: `ldp` of 16
words a time, `and.16b`, `cnt.16b`, `udot.4s`, then `uaddlp.2d`,
`uzp1.4s` and `add.4s` because the sum is asked per word, 48 SIMD
instructions per 16 word pairs, no bounds check and no call; in wasm32
it is `i64.and`, `i64.popcnt`, `i32.add`, unrolled by two, no call.

The counts, the same on every run: 999 rayon work items per block, the
largest row 999 pairs and the smallest 1; 79 words per set; 5 sets per
individual, 395 words, 3.16 MB per block; 790 words walked per pair,
6320 bytes, 63.1 GB over the 9990000 pairs of a pass; one `Error` built
per called allele, 1.94e8 a pass; about 60 allocations per block and
none per pair or per genotype.

## 4. The measurement plan

In the order in which each unblocks a finding.

1. **The profile of the wasm run.** None exists, and the wasm target is
   the one missed by most. `build:wasm` of `js/popnei/package.json`
   passes `--remove-name-section` to wasm-bindgen, so node's profiler
   shows `wasm-function[N]`. Build once without that flag, then

       node --cpu-prof --cpu-prof-dir=tmp/ bench/time_kosman_dists.mjs /Users/jose/devel/popnei-bench/big.vars 5

   from `js/popnei`, and read the self times of the `.cpuprofile`. What
   it answers: the split of the sets phase and the pairs phase in wasm,
   within a tenth, as natively.
2. **The gates of the two calls per genotype**, H1 and H2 below:
   `cargo asm -p popnei --lib "of_block"` before and after, the `bl
   core::ptr::drop_glue::<popnei::error::Error>` and the `bl
   core::slice::memchr::memchr` inside the loop over the genotypes gone;
   then the profile again, the two frames under 1 in 100; then the
   bench, the blocks in memory on one thread.
3. **The sets phase on the threads**, H3: the bench's blocks in memory
   on 18 threads with the phase split, the sets phase below 0.325 s,
   the one thread line not above its 1.138 s by more than 2 in 100, and
   the tests that compare the parallel and the serial adding green.
4. **The words per pair**, H4: the count of words of `holds` per
   individual, 316 to 158, asserted in a test; the pairs phase in the
   bench, one thread, below 0.45 s from 0.593 s; and a test that the
   biallelic path and the general path give the same two integers for
   every pair of one block whose variants are not a multiple of 64.
5. **The build flags**, section 5: `lto` fat and thin and `codegen-units
   = 1`, each alone, on the bench through `--config`, and on the wasm
   build through the node timing script, kept at a gain above 5 in 100.
6. **The cache effect on the sets phase**, H5: the bench's block in
   memory at 5000 x 1000 against 5000 x 100 run ten times, the same
   genotypes, a tenth of the write set; a lower time per genotype in the
   second is the effect.
7. **Two more shapes in the bench**, for what the spec leaves open: a
   block of 4 alleles and one of 500 variants x 10000 individuals, from
   the same generator, so that the biallelic kernel is not measured on
   its best case alone and the block size at 10000 individuals gets a
   number.
8. **The wheel's profile**: `[profile.release] debug =
   "line-tables-only"` so that `sample` on a Python process inside a
   call shows names; binary size is the cost.

## 5. The build configuration

`[profile.release]` sets `overflow-checks = false` and nothing else:
opt-level 3, no LTO, 16 codegen units, no `target-cpu`. The bench
profile inherits it, plus the line tables of 59e76a5. The wasm builds
take the same profile: `npm run build:wasm` is `cargo build --package
popnei-js --release --target wasm32-unknown-unknown`, and the pyodide
wheel is built by maturin in release; neither passes `+simd128`, and
`wasm-opt` is not installed (binaryen is not on the machine).

- **LTO and codegen units** [Likely, build]. The two out of line calls
  of the sets phase are what cross crate inlining removes, and `docs/
  reports/vars-file.md` already named LTO as the one setting worth
  trying. Experiment: `cargo bench --bench kosman_dists --config
  'profile.bench.lto="fat"' ...`, then `"thin"`, then `codegen-units=1`,
  each alone, best of 5, kept at a gain above 5 in 100; the same on the
  wasm build through the node script. The cost is link time on every
  build. After H1 and H2 the calls are gone anyway, so this is measured
  after them.
- **`target-cpu`** is not a candidate: the stock build already emits
  `and.16b`, `cnt.16b` and `udot.4s` for the pairs loop, since NEON is
  in the aarch64 baseline. An allocator swap is not one either, at 60
  allocations per block. PGO is premature while two named frames carry
  a fifth of the CPU.
- **`+simd128` for wasm** [Likely, build]: alone it gives nothing. Built
  with `RUSTFLAGS="-C target-feature=+simd128"`, `sums_of_two` gets
  `v128.and` and `v128.load` but still six `i64.popcnt` and no
  `i8x16.popcnt`, since LLVM extracts the lanes to count them scalar. A
  vector popcount in wasm needs the intrinsics of `core::arch::wasm32`,
  which are safe functions, behind `cfg(target_feature = "simd128")`,
  the flag in both wasm builds, and a floor on the browsers popnei
  supports, which is the owner's decision. Filed as L5.
- **`wasm-opt -O3`** [Speculative, build]: pattern only; needs `brew
  install binaryen`.

## 6. The findings

Each names its file and line at 59e76a5, its severity as
`finding_format.md` of the skill defines it, its gate, what it does to
the numbers and what it costs. All keep the two integers of every pair
the same at every thread count and block size. What each experiment gave
is in section 9, which was written as they ran.

### Hot path

**H1. An `Error` is built and dropped for every allele.**
`crates/popnei/src/dists.rs:214`, `.ok_or(Error::AlleleBelowTheMissingOne
{ allele })?` in the loop over the alleles of a genotype. Evidence: the
profile, `drop_glue::<Error>` 10.3 in 100 of the CPU on one thread and
10.1 on 18, and the assembly, the construction and the out of line drop
on the success path of the loop; the same call in the wasm32 assembly.
Mechanism: `ok_or` takes its argument by value, so the variant, an `i8`
in an enum that has `String` and `PathBuf` variants and so has drop
glue, is written and dropped 1.94e8 times a pass; `ok_or_else` with a
closure moves both into the `None` arm, which no valid block reaches.
Gate: the `bl core::ptr::drop_glue` inside the loop of `of_block` gone
from `cargo asm`, then the frame gone from the profile, then the bench
on one thread. Numbers: none. Cost: one word, the form the same
function already uses four times. Agreed by allocations, hot_loops,
numbers, methodology and linalg_and_wasm.

**H2. `contains` on a two byte genotype calls `memchr`.**
`dists.rs:186`, `genotype.contains(&MISSING_ALLELE)`. Evidence: the
profile, 9.4 and 8.2 in 100, and the assembly, `bl
core::slice::memchr::memchr` per genotype, natively and in wasm32.
Mechanism: `[i8]::contains` dispatches to `memchr`, a word at a time
search with its own setup, for a slice of the ploidy, 1e8 times;
`iter().any(|&a| a == MISSING_ALLELE)` inlines to one or two compares.
Gate: the call gone from the assembly and the frame from the profile,
then the bench. Numbers: none. Cost: none. Agreed by hot_loops,
allocations, linalg_and_wasm and methodology.

**H3. The sets of a block are built on the calling thread while the
pool waits.** `dists.rs:103-226`. Evidence: the phase split, the sets
phase 0.569 s of the 0.624 s on 18 threads and the same 0.529 s on one;
the profile, 1.91 cores busy on average over 18 threads, the sets phase
78 in 100 of the main thread. Mechanism: `of_block` writes every bit on
one thread. The writes are disjoint per individual: individual i owns
`called[i * 79..]` and `holds[i * 316..]`, so a range of c individuals
is one work item of `par_chunks_mut` over the two arrays, no lock and
no atomic, the pattern the VCF reader already uses over its rows. The
loop stays variant major inside the item so that `gts` is still read as
a stream, 2c bytes of each 2000 byte row per item, which reads each
128 byte line 4 times over the items at c = 14 instead of once, 40 MB a
block instead of 10. A split by variants would be wrong at any range
that is not a multiple of 64, two threads ORing into one word, and is
not proposed. Gate: the sets phase of the bench in memory on 18
threads below 0.325 s, the one thread line within 2 in 100 of 1.138 s,
`cargo test --workspace` green with the test that compares the parallel
and the serial adding, and a new test on a block whose variants are not
a multiple of 64 with both paths. Numbers: none, each bit is written
once by one thread. Cost: a serial `of_block` kept for wasm and for the
test that compares the two, as `parse_rows` of the VCF reader has, a
chunk size, an `#[expect]` for the chunk arithmetic, and which allele
`AlleleBelowTheMissingOne` names when a reader has a defect, which the
filters solve by running the serial pass again to name the first.
Agreed by concurrency, numbers, hot_loops, data_layout,
linalg_and_wasm and methodology.

**H4. A biallelic block needs 1 + k sets per individual, not 1 + 2k.**
`dists.rs:411-425` and the layout at 103-236. Evidence: the profile, the
pairs closure 47.9 in 100 on one thread, 0.593 s of 1.123 s; the
assembly, the loop at the SIMD issue width, about 4.5 instructions per
cycle from a working set in L2; the spec, which leaves "whether a block
with only the alleles 0 and 1 gets fewer sets" to the implementer.
Mechanism: with two alleles the copies of allele 1 are the ploidy less
those of allele 0, so `holds(1, m) = called AND NOT holds(0, k + 1 -
m)` and the two sets of allele 1 carry nothing the others do not. Two
exact forms: the AND form, where the count of the partner set is
`popcount(C AND NOT h_i AND NOT h_j)` with `C` the called of both, and
the XOR form, where the ploidy times the sum of d is
`popcount((h_i XOR h_j) AND C)` summed over the two sets of allele 0.
Either walks 3 sets per individual, 237 words, 474 per pair against
790, and writes 3 sets per individual in the sets phase. A block with
a third allele or another ploidy keeps the general path, chosen per
block on `num_alleles == 2`. The complements must stay masked by a set
of the individual, so that the padding bits of the last word stay 0.
Gate: `holds` words per individual 316 to 158 asserted in a test, the
pairs phase in the bench on one thread below 0.45 s from 0.593 s, and a
test that the biallelic and the general path give the same two integers
for every pair of a biallelic block whose variants are not a multiple
of 64, with missing genotypes and an individual homozygous for the
alternative allele throughout, run also on the four reference files.
Numbers: none. Cost: a second kernel and a second builder arm behind a
flag on `KosmanBits`, the largest of this report, so it comes after
H1 to H3 and its gain is read against what they leave. Agreed by
data_layout, hot_loops and numbers.

**H5. The write set of the sets phase is 640 KB against a 128 KB L1.**
`dists.rs:179-226`. Evidence: `of_block` self 20.7 in 100 on one thread
after the two calls above, pattern for the rest, since the machine has
no cache counters. Mechanism: the loop is variant major, and within a
run of 64 variants each individual is written at one `called` word and
up to four `holds` words 632 bytes apart, 5 lines of 128 bytes per
individual, 5000 lines per 64 variants, 5 times the L1 of a performance
core and 10 times an efficiency core's; every bit write is a miss to
L2. Two reorderings, both leaving the pairs loop as it is because it
sums over all the `holds` words in any order: interleave `holds` word
wise, `holds[i * 316 + w * 4 + s]`, so that an individual's four words
of one index share a line, 256 KB per span; and tile the individuals,
256 per span of 64 variants, 64 KB. Gate: plan item 6 first, the block
at 5000 x 100 ten times against 5000 x 1000; then the sets phase in
memory on one thread below 0.45 s from 0.529 s. Numbers: none, the
same bits at other addresses. Cost: the index in the builder, the doc
of the layout, two nested loops and an `#[expect]` for the tile.
Filed by data_layout and hot_loops; medium confidence.

### Likely

**L1. The count of copies is three loops for a trip count of one or
two.** `dists.rs:199-203`, `iter().take(at + 1).filter(..).count()` per
allele. Evidence: the assembly, a 32 byte NEON loop, an 8 byte one and a
scalar tail with two range branches, for at most two compares, and
O(k²) in the ploidy. Experiment: an arm for ploidy 2, one compare,
with the general loop kept and a test that both give the same sets;
gate on the instruction count of `of_block` and the bench on one
thread. Numbers: none. Cost: one branch. Filed by hot_loops.

**L2. The rank lookup per allele is a fallible conversion, a bounds
checked load and the error of H1.** `dists.rs:211`. Both failures are
unreachable, since line 153 refuses a block whose smallest allele is
below the missing one and line 186 skips a missing genotype, so every
allele reaching the lookup is 0 to 127. `allele.cast_unsigned()` into a
`[u8; 128]` table is infallible and 128 bytes instead of 1 KB. Gate: the
loop body of `of_block` in the assembly, then the bench. Numbers: none.
Cost: the invariant written down. Filed by numbers.

**L3. The popcount reduction is lowered at about 3 SIMD instructions per
word pair where 2 would do.** `dists.rs:412-423`. Evidence: the
assembly, `cnt.16b` then `movi`, `udot.4s`, `uaddlp.2d`, `uzp1.4s` and
`add.4s` per 16 bytes because `.map(count_ones).sum()` asks for a per
word `u32`; `cnt.16b` with `uadalp` into 16 bit lanes, which 395 words
never overflow, is 3 per 16 bytes instead of 5. The integers are the
same in any order. Experiment: a prototype with a SIMD crate, `wide` or
`fearless_simd`, since `std::simd` is nightly, checked with `cargo asm`,
then the bench on one thread. Cost: a dependency, a second kernel for
wasm, a test that the two agree. After H4, which removes 40 in 100 of
the words and may make it not worth its cost. Filed by hot_loops,
numbers and data_layout.

**L4. `alleles_of` is a serial pass over the block that the parallel
build leaves as the remainder.** `dists.rs:364-395`, 3.2 in 100 on one
thread, about 0.040 s of the sets phase. A fold with OR and `min`
reduces exactly on `par_chunks`. After H3; gate on the sets phase.
Cost: a reduce closure and its serial twin. Filed by concurrency and
hot_loops.

**L5. No vector popcount in wasm.** Section 5, `+simd128` alone gives
none; a hand written path over `core::arch::wasm32` needs the flag in
both wasm builds and a floor on the browsers. For the owner, after the
wasm profile of plan item 1 says how much the pairs phase is there.
Filed by linalg_and_wasm and methodology.

**L6. At 10000 individuals the block size becomes a cliff.**
`dists.rs:279-304`, pattern only at that size. The sets of a 5000
variant block are 31.6 MB, above the 16 MB L2 of a performance cluster,
so each row's sweep comes from memory; at 500 variants the sets fit
again but the 400 MB accumulator is read and written 2000 times over a
million variants. Tiling the rows, T rows against one sweep of the
sets, divides the bytes streamed per pair by T with the instruction
count unchanged. Plan item 7 first. Cost: a tiled kernel with T
accumulators and the diagonal of a tile apart. Filed by data_layout.

**L7. The wasm reading alone over subtracts.** The timing script's
reading alone reaches `gts()` of the wasm binding, which copies each
block out of wasm memory, which the calculation never does; the tests
reviewer of the plan timed that copy at 5 to 10 ms of the 0.161 s, so
the 2.130 s is a low end by that much. Plan item 1 settles the rest.
Filed by methodology.

### Speculative

**S1. The handoff to the pool when the pool has one thread.**
`dists.rs:724`. The main thread waits on a condition variable while the
one worker counts, 48 in 100 of its samples, but the phase split adds up
to the whole to 1 ms, so the handoff leaves no visible residue: 20
injections a pass against 0.593 s of counting. Experiment: an early
return to the serial adding when `current_num_threads() == 1`, the
bench on one thread; expected to show nothing and close. Filed by
concurrency.

**S2. The accumulator rows share a line of 128 bytes at each boundary.**
`dists.rs:706-723`, 999 boundaries against 31219 lines per block, 3 in
100, two threads each; 0.6 MB of padding at 10000 individuals and a
changed pair index to remove them. Not worth trying without a measured
gain. Filed by data_layout.

### Notes

- **N1.** `python/popnei/dists.py:214`: `pandas.DataFrame(square, ...)`
  copies its array under pandas 3.0, so `square_dists` holds 800 MB
  twice at 10000 individuals; `copy=False` gives away nothing, the
  array being a fresh local. And `triu_indices` builds two int64 arrays
  of 800 MB for the 800 MB result, where a loop over rows needs none;
  `triang_list_of_lists` pays the frame's copy for a frame it discards,
  on top of the 1.6 GB of Python floats that Biopython's form costs.
  Cold paths, memory bound; gate on peak memory with `tracemalloc`.
- **N2.** `crates/popnei-js/src/dists.rs:60`: the vector of distances is
  held twice in wasm memory while wasm-bindgen copies it to the
  JavaScript heap, 800 MB of a 4 GB heap at 10000 individuals, and the
  accumulator is alive with it while `dists().collect()` runs, 400 MB
  more; the vars file source already crosses in pieces of 1 MiB. For
  the spec, which should state the peak.
- **N3.** `sample` on 18 threads gives shares of CPU time and not of
  wall time; the report says which.
- **N4.** The phase split rests on an `Instant` harness that was not
  committed; before the next split is quoted it goes behind a cargo
  feature of the bench.

## 7. Seen outside the scope

- `dists.rs:192` and `:221`: a bit is silently dropped if `get_mut`
  gives `None`; right today by the bound, and any change of the indexing
  makes it a silent wrong count. For the code review of whatever changes
  the layout.
- `dists.rs:275` in the assembly: `<ChunksExactMut<u64> as
  Iterator>::zip` called out of line once per variant, 1e5 a pass; LTO
  or an inline hint.
- `scripts/build_pyodide_wheel.sh`: pyodide's link flags carry `-Oz`,
  and whether that reaches a Rust side module was not checked.
- `docs/reports/kosman-method/kspike/src/lib.rs:59-76`: the trial's set
  builder, `pack`, has neither the error nor the `memchr` and took
  0.015 s a block where `of_block` takes 0.029 s; the 18 core target
  was set from it. A binary that times both on one block would say
  whether the target is reachable at all; H1, H2 and L1 are the
  difference between the two.

## 8. What the code already does well

- `dists.rs:411-425`, `sums_of_two`: two flat zips over contiguous
  slices, vectorized to 16 words an iteration with four accumulators,
  no bounds check and no call; in wasm the same loop with
  `i64.popcnt`. The shape to keep.
- `dists.rs:688-737`: the rows of the accumulator cut with
  `split_at_mut` and handed to rayon with no lock, no atomic and no
  allocation per pair; the profile shows the 18 workers even.
- `crates/popnei-python/src/dists.rs:58-92`: one crossing per call, the
  whole pass inside `py.detach`, the vector handed to numpy with no
  copy; the Python timing and the Rust bench agree to 1 in 100.

## 9. The experiments

Run one at a time on this branch, each on the commit the one before it
left, with the bench command of section 3, best of 5 after an untimed
run, on the owner's M5 Pro, the load average beside each. The two
settings are the blocks handed out from memory, which is the
calculation with no reader in it and so the setting the targets are of,
and the same calculation over `big.vars`. The whole of the table is the
"blocks in memory" line of the bench.

| after | one thread | 18 threads | what it changed |
|---|---|---|---|
| the baseline, 80b7b85 | 1.147 s | 0.635 s | |
| H1, 7d29ae4 | 0.995 s | 0.467 s | the error not built per allele |
| H2, 371aaf7 | 1.007 s | 0.475 s | the call to `memchr` gone |
| L1, 67b7a16 | 0.835 s | 0.303 s | one compare for the copies at ploidy 2 |
| L2, e5b0f25 | 0.817 s | 0.286 s | the rank lookup infallible |
| H3, 1bd476a | 0.767 s | 0.102 s | the sets built on the threads |
| the targets | 0.97 s | 0.38 s | |

Over the vars file, whose reader is not in the targets, the same six
commits take it from 1.264 s to 0.891 s on one thread and from 0.737 s
to 0.229 s on 18 threads; of that last number the reader is 0.127 s,
more than the calculation it feeds.

The phase split after H3, by the same `Instant` harness as in section 3,
not committed, over the blocks in memory:

| | whole | the sets phase | the pairs phase |
|---|---|---|---|
| 18 threads | 0.103 s | 0.054 s | 0.049 s |
| one thread | 0.769 s | 0.172 s | 0.603 s |

**H1, the error not built per allele. Applied, 7d29ae4.** The gate was
the count: `bl core::ptr::drop_glue::<popnei::error::Error>` in the
assembly of `of_block` went from 2 to 1, the one inside the loop over
the alleles gone. One thread 1.147 to 0.995 s, 18 threads 0.635 to
0.467 s, 13 and 26 in 100. `cargo test --workspace` 345 passed. It cost
an `#[expect(clippy::unnecessary_lazy_evaluations)]`, which the lint
demanded because it reads the closure as needless, and which L2 later
removed.

**H2, the call to `memchr` gone. Applied on the count, 371aaf7.** `bl
core::slice::memchr::memchr` went from 1 to 0 and `of_block` from 531 to
507 instructions, but the wall time did not move: 0.995 to 1.007 s on
one thread and 0.467 to 0.475 s on 18, both inside the 2 in 100 that
two invocations of the bench agree to. Kept because the gate is the
count and the change costs nothing, with an
`#[expect(clippy::manual_contains)]`, the lint holding `contains` to be
the faster of the two, which on a slice of two bytes it is not.

**L1, one compare for the copies at a ploidy of 2. Applied, 67b7a16.**
The largest single gain: one thread 1.007 to 0.835 s, 17 in 100, and 18
threads 0.475 to 0.303 s. Its assembly gate was not met as written,
because the general loop stays for every other ploidy and the function
grew from 507 to 668 instructions; the wall time decided it, twice. No
test was added: the diploid worked example holds homozygous,
heterozygous, half called and fully missing genotypes, and the arm made
wrong on purpose fails it together with seven others, the comparison of
every pair of the four reference files with R's `gd.kosman` among them.

**L2, the rank lookup infallible. Applied, e5b0f25.** The table of the
places of the alleles becomes a `[u8; 128]` read with
`allele.cast_unsigned()`, total because a block whose smallest allele is
below the missing one is refused before the loop and a missing genotype
is skipped inside it; the invariant is written where the lookup is.
`of_block` 668 to 659 instructions. One thread 0.835 to 0.817 s, 1.8 in
100, and 18 threads 0.303 to 0.286 s, 5.6 in 100, each measured twice a
side with a spread of 0.001 s inside a round. It was first closed for
being under a gate of 2 in 100 set on the one thread number, then taken:
the effect repeated, it is above the gate on 18 threads, and it removes
an `#[expect]` instead of adding one.

**H3, the sets of a block built on the threads. Applied, 1bd476a.** The
sets of 64 individuals are one work item of `par_chunks_mut` over the
two arrays, the loop variant major inside the item, behind
`cfg(not(target_family = "wasm"))` with the serial builder beside it for
wasm and for a new test that compares the two on a block whose variants
are not a multiple of 64. Nothing inside an item can fail, since L2 made
the lookup total and the block is refused from its smallest allele
before any item runs, so the items need no error path: L2 made H3
simpler. The sets phase on 18 threads 0.229 to 0.054 s and the whole
line 0.286 to 0.102 s, 64 in 100. The individuals per item were swept
over 8, 14, 32 and 64, and 64 won on both thread counts although it
leaves 16 items for 18 threads; the sweep never turned around, so the
best value may be above 64, which the constant's comment says. 346 tests
pass.

Two things came out against the finding. The gate, the sets phase below
0.325 s on 18 threads, was already met before H3 ran, since H1, L1 and
L2 had taken that phase from the review's 0.569 s to 0.229 s; the whole
line is what H3 was judged on. And the one thread line fell, 0.817 to
0.767 s, where the finding expected it to rise from reading each row of
genotypes once per item: at 64 individuals an item writes into 202 KB of
the block's 3.16 MB of sets, which is the cache effect of H5 got for
free, and H5's own reorderings have that much less left to take.
