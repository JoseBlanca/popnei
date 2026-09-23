# Performance review: the distances between populations

23 September 2026. The first measurement of `calc_pop_dists`, the
calculation merged into `main` that day in 05f1538, which gives for every
pair of populations, out of one reading of the variants, seven measures
with a block jackknife standard error each. Nothing of it had ever been
timed: no benchmark, no profile, no wall time, for the core, for Python or
for wasm. The owner decided on 23 September 2026 that the plan which built
it would not measure and that this review would, and "Speed" of
`docs/specs/dists.md` asks this review to write the numbers to reach into
it once they exist. Section 2 has the numbers, section 3 what they are of,
section 10 the numbers proposed for the spec and the one decision they
need from the owner.

The words this document uses. The **calculation** is one call of
`calc_pop_dist_sums` over a reader, which gives, for each pair of
populations, six sums over the variants that counted for that pair; the
seven measures are ratios of those six and are divided out once, at the
end. A **pass** is one reading of every variant. The **counts of a
population at a variant** are the three numbers everything is built from:
how often each allele was called there, how many genotypes were called
whole, and how many of those are heterozygous. The **pair arithmetic** is
what turns the counts of two populations at one variant into what that
variant adds to their six sums. The **counter table** is the array of 128
counters, one for each allele a variant could hold, that the counting
fills. A **resampling group** is a stretch of a chromosome that the
standard errors leave out in turn. The **reader** is what turns a file
into blocks of genotypes, and it is outside every number called "the
calculation" below.

## 1. The scope and its limits

Reviewed: `crates/popnei/src/pop_dists.rs`, `crates/popnei/src/variant.rs`
and `crates/popnei/src/stats.rs` where the pass reaches them,
`crates/popnei-python/src/pop_dists.rs`,
`crates/popnei-js/src/pop_dists.rs`, `python/popnei/pop_dists.py` and
`js/popnei/src/pop_dists.ts`, at `main` 3bc33f9, where the five files of
the module are byte for byte what 05f1538 left; the branch of this review
is `perf/dists-pops` and it is not merged.

The machine is the owner's Apple M5 Pro, **6 performance cores and 12
efficiency cores**, 64 GB, macOS 27.0, rustc 1.98.0, node 26.8.2.

The dataset is the one "Speed" names: 100000 variants x 1000 diploid
individuals, biallelic, 3 in 100 genotypes missing, `big.vars` at
`/Users/jose/devel/popnei-bench/`, two chromosomes of 50000 variants 1000
base pairs apart, cut into 100 resampling groups of 1000000 base pairs,
with `min_num_individuals` at its default of 20. The individuals are cut
into 3 populations and into 20, which the spec asks for, and into 6, 10
and 40 besides, which section 4 needs. `crates/popnei/benches/make_pops.py`
writes those files; the 3 are the ones the dataset was simulated with, 296,
356 and 348 individuals, and the rest are runs of the file's order.

Left out: any ploidy but 2 and any variant of more than two alleles, so
every per-pair number here is a two-allele number and a microsatellite
panel, which this module is specified for, is untimed; the vars reader,
which was reviewed on 22 September 2026 in
`docs/reports/perf-io-2026-09-22.md` and whose cost is outside the targets,
though section 8 says what it now means for a user; and reading the file
from the disc, which every number here excludes, the bytes being in memory
before any clock starts.

Categories sent, one reviewer each with a fresh context: methodology,
numbers, allocations, data_layout, concurrency, hot_loops,
linalg_and_wasm, python_boundary and io_and_syscalls.

## 2. The verdict

**Run the experiments.** The code has no defect of speed that this review
can point at and no target it misses, because it had no target. What it
has is two costs that a profile names, that a cost model measures, and
that a change of contained size removes, both of them bit-identical: the
counter table (H1) and the per-population arithmetic redone for every pair
(H2). Section 6 gives them in order and section 9 what each experiment
gave.

What the code does today, on the dataset above with the reader taken out,
on one thread, best of 5 on a quiet machine:

| | 3 populations, 3 pairs | 20 populations, 190 pairs |
|---|---|---|
| popnei, one thread | 0.282 s | 0.489 s |
| popnei in wasm, one thread | 0.375 s | 0.582 s |
| pyNei, one thread | 7.811 s | 502.105 s |

Over the vars file, which is what a user waits for and has the reader in
it, popnei takes 0.388 s and 0.620 s on one thread against pyNei's 7.811 s
and 502.105 s, so it is **20 times faster at 3 populations and 810 times
at 20**. The second ratio is mostly pyNei's own shape and not popnei's
speed, and section 5 says why.

What is asked of the owner: the order to merge this branch, which nothing
here stands in the way of; and one decision, in section 10, on whether the
spec's numbers to reach should include a figure for many threads at all,
given that on this machine such a figure moved by a factor of two between
two runs of the same binary while the one-thread figures moved by under 6
in 100.

## 3. What was measured, and with what

The benchmark is `crates/popnei/benches/pop_dists.rs`, added by this
review, `harness = false` as the other eight, run as

    cargo bench --bench pop_dists -- /Users/jose/devel/popnei-bench/big.vars \
        --pops /Users/jose/devel/popnei-bench/pops3.tsv \
        --pops /Users/jose/devel/popnei-bench/pops20.tsv --runs 5

It times two settings, as `kosman_dists.rs` does. **The blocks in memory**
are 20 blocks of 5000 variants built by a seeded generator before the
clock starts, carrying the chromosome and the position as well as the
genotypes, laid out like `big.vars` so that both settings fall into the
same 100 groups: no reader is inside that number, and it is the one the
targets are of. **The vars file** is the same calculation over a
`VarsReader` on bytes already in memory, so its difference from the first
is the reader.

Best of 5, one thread, on a machine checked quiet with `sysctl -n
vm.loadavg` and `ps` before and after:

| | 3 populations | 20 populations |
|---|---|---|
| the blocks in memory | 0.282 s | 0.489 s |
| the vars file | 0.388 s | 0.620 s |
| the reader, by difference | 0.106 s | 0.131 s |

and on the 18 threads of the process, 0.027 s and 0.044 s in memory,
0.134 s and 0.155 s over the file.

Three other harnesses, all added by this review:

- `crates/popnei/benches/time_pop_dists.py` times the call from Python
  with a pass that reads the same fields and does nothing with them
  measured beside it and subtracted: 0.290 s and 0.503 s on one thread,
  0.034 s and 0.051 s on all. Those are 8 to 14 ms above the benchmark's
  in-memory numbers, and section 7 says what that is.
- `crates/popnei/benches/time_pop_dists_pynei.py` times pyNei's
  `calc_jost_dest_pop_dists` over `big.pynei.vars`, the same genotypes in
  pyNei's own format.
- `js/popnei/bench/time_pop_dists.mjs` times the wasm build under node.

The sampling profile is `/usr/bin/sample <pid> 45 1` over one run of the
benchmark per populations file, with the idle rayon frames removed:
`__psynch_cvwait` alone is 89 in 100 of all samples, which is the workers
of a pool asleep through the one-thread timings. Of the samples that are
work:

| | 3 populations | 20 populations |
|---|---|---|
| `pop_dists::sums_of_the_chunk` | 3.4% | 32.2% |
| `variant::count_gts_of` | 32.9% | 21.1% |
| `variant::count_the_alleles` | 32.6% | 20.2% |
| `variant::count_alleles_of` | 12.8% | 12.5% |
| `lz4_flex ... read_to_end`, the reader | 12.2% | 7.2% |
| `_platform_memset` | 1.2% | 3.7% |
| `_platform_memmove` | 1.7% | 1.0% |
| `drop_glue::<popnei::error::Error>` | — | 0.6% |

`sums_of_the_chunk` is where the pair arithmetic was inlined:
`add_the_pairs`, `of_var` and `over_the_alleles` have no frames of their
own. So the counting of genotypes is 78 in 100 of the work at 3
populations and 54 at 20, and the pair arithmetic 3 and 32.

Three limits of that profile, from the methodology review. Its percentages
are shares of CPU across every thread and not of wall time, so the
reader's 7.2% is of the CPU while the reader is 0.111 s of the 0.134 s a
user waits at 3 populations on 18 threads. It covers all four timings of a
populations file at once. And `_platform_memmove` is the same 407 and 413
samples at 20 and at 3 populations, which is a cost that does not move
with the populations: it is the benchmark building its 200 MB of blocks
before each clock, and the reader's own copies, not the pass.
`_platform_memset` does move with the populations, 296 to 1585, and the
two `vec!` calls per chunk cannot produce it, since a vector of a user
struct is filled by a clone loop and not by a memset: it is the counter
table of H1, which is 4.1 GB zeroed per pass at 20 populations against 31
MB for those vectors.

## 4. How the cost grows with the populations, which the spec asked for

"How it runs" of the spec says the cost grows with the square of the
populations, "per variant it is one pass over the genotypes for the
counts, P of them for P populations, and then P(P-1)/2 pairs times the
alleles of the variant for the sums", and that no measurement of either
had been made. Here it is, the same 1000 individuals cut five ways, blocks
in memory, one thread, best of 5:

| populations | pairs | individuals counted | one thread |
|---|---|---|---|
| 3 | 3 | 1000 | 0.287 s |
| 6 | 15 | 1000 | 0.314 s |
| 10 | 45 | 1000 | 0.356 s |
| 20 | 190 | 1000 | 0.498 s |
| 40 | 780 | 1000 | 0.962 s |
| 3 of 50 individuals | 3 | 150 | 0.066 s |

The five full-panel rows fit, each within 0.4 in 100,

    a pass of 100000 variants =
        0.268 s x (individuals / 1000) + 0.0060 s x populations
        + 0.00058 s x pairs

which per variant is **2.68 ns for each individual counted**, **60 ns
fixed for each population at each variant, whatever its size**, and
**5.8 ns for each pair at each variant** at two alleles.

Two things were done to test that model rather than assert it. It was fitted
on 3, 6, 10 and 20 populations and asked to predict 40, which it missed by
7.3 in 100, so the coefficients above are the refit on all five and the
four-point version is not used. And the last row breaks the confound that
every other row shares, that cutting the same 1000 individuals more ways
makes each population smaller: 3 populations of 50 name only 150 of the
1000 individuals, so the populations and the pairs are those of the first
row and the genotypes counted are a sixth. The model, fitted without it,
puts it at 0.060 s against 0.066 s measured. So the first term really is
per individual counted and the second really is fixed per population.

What that answers. Going from 3 populations to 20 multiplies the pairs by
63 and the time by 1.7, and of the 0.211 s it adds, **0.102 s is the fixed
cost of visiting each population at each variant and 0.109 s is the pairs**
— about equal. The growth with the square of the populations is real and
it does not dominate until past 20 populations: at 40 the pairs are
0.452 s against 0.241 s. On a microsatellite panel the pair term would be
larger, since it is the only term that grows with the alleles of a
variant, and that was not measured.

## 5. pyNei, and what the comparison means

pyNei's `calc_jost_dest_pop_dists` over the same 100000 variants x 1000
individuals, its own reading inside its number as popnei's vars-file rows
have popnei's reader inside theirs:

| | one thread | 6 threads | 18 threads |
|---|---|---|---|
| 3 populations | 7.811 s | 1.765 s | 0.944 s |
| 20 populations | 502.105 s | not run | not usable |

The 20-population run on 18 threads gave 74, 93 and 116 s, rising run over
run while another process compiled on the machine, and is not reported as
a measurement. The one-thread run is 502.105 s with the untimed run before
it at 513.250 s, 2 in 100 apart.

Two things have to be said wherever those numbers sit beside popnei's.
**The comparison is one pass against one pass and not one number against
one number**: popnei calculates seven measures with a block jackknife
standard error each, and pyNei calculates Jost's D alone with none, from
counts the two share. And **pyNei's cost grows with the pairs in a way
popnei's does not**: `_DestPopHsHtCalculator.__call__` loops over the pairs
and hands `_calc_pairwise_dest` the whole dictionary of populations each
time, so every pair recounts the alleles and the called genotypes of every
population of the chunk. Its own numbers show it: 502.105 s against
7.811 s is 64.3 times, where the pairs are 63.3 times. So the 810-fold
ratio at 20 populations is mostly a measure of that repeated counting.
The 20-fold ratio at 3 populations is the fairer one to quote.

## 6. The findings

Numbered so that the prose can refer to them. Each says what it removes,
what it does to the numbers, and what would confirm it.

### H1 — The counter table is zeroed, merged and scanned in full for every population at every variant

`crates/popnei/src/variant.rs:743` and `:749`, with
`crates/popnei/src/pop_dists.rs:123`. Hot-path, confidence high; named
independently by hot_loops, allocations and data_layout.

`count_alleles_of` builds `[[0u32; 128]; 4]`, 2048 bytes, on the stack for
every population at every variant; `merge_the_lanes` then reads those 2048
bytes and writes all 512 bytes of the counts; and `count_the_var` then
scans those 512 bytes backwards with `rposition` to find the largest
allele, which on this dataset is 1. That is about 5.1 KB moved for each
population at each variant to count the 100 alleles of a population of 50
individuals: at 20 populations, 10.2 GB per pass against 200 MB of
genotypes read. 5.1 KB in the measured 60 ns is about 85 GB/s, which is
the bandwidth of the first level cache and not an instruction count, so
the fixed per-population cost of section 4 is this table and little else.

`count_alleles`, the whole-row counting, already has a two-allele path,
`the_counts_of_a_variant_of_two_alleles`, which never touches the table;
`count_alleles_of`, the per-population one, has none, and a pairwise pass
can never reach the whole-row path. The statistics review of 22 September
2026 built the neighbouring change, bounding the clear, the merge and the
scan to the alleles seen, measured it at 1 and at 4 populations where it
gave nothing, and closed it with this: "At the 50 populations of 20
individuals that `docs/objectives.md` calls ordinary, the same change takes
the benchmark from 0.706 s to 0.448 s, a gain of 37 in 100. That is where
it belongs, and it is the owner's to ask for." This module is that place.

The numbers cannot move: these are `u32` counts and integer addition is
exact whatever the order. What must hold is that no entry of the counts
above the bound is ever read, or a pair reads the previous population's
allele as this one's.

Measured by: the bytes zeroed per pass, which is the same on every run;
then the refit of section 4, where the coefficient that must fall is the
0.0060 s per population, not the total; then `_platform_memset` gone from
the profile.

### H2 — Four of the seven per-allele sums, and ten of the fifteen divisions, are of one population and are redone for every pair it is in

`crates/popnei/src/pop_dists.rs:329` and `:281`. Hot-path, confidence
high; named independently by hot_loops, numbers and data_layout.

`over_the_alleles` walks the alleles once for each pair and computes,
inside that loop, each population's own allele frequency, the square of
it and its ploidy-th power. `of_var` above it divides for each
population's observed heterozygosity, its correction `n/(n-1)` and the
reciprocal of its called genotypes. Counting the divisions of one pair at
one biallelic variant gives 15, of which 10 depend on one population
alone. At 20 populations that is 2850 divisions a variant where 1150 would
do, and 760 square-and-power evaluations where 40 would do. Division is
the longest-latency operation in that loop.

Bit-identical, with one condition. Each hoisted value is one correctly
rounded operation on the same operands, and aarch64 has no excess
precision and rustc emits no contraction flag, so where it is evaluated
cannot change its bits; the alleles are still added in allele order. The
loop bound changes from the larger of the two populations' allele counts
to the population's own, and every term dropped is exactly `+0.0` added to
a sum that starts at `+0.0` and never takes a negative term. The condition:
`of_var` returns nothing, before it divides, when a population has fewer
than two called alleles or both have exactly one called genotype, and on
that path the hoisted `n/(n-1)` and `1/called` would be infinities, so the
hoisted values must sit behind the same integer tests.

Measured by: the count of `fdiv` in `sums_of_the_chunk` from `objdump -d`
on the built benchmark, not `cargo asm`, which forces its own
`codegen-units=1`; then the refit, where the coefficient that must fall is
the 0.00058 s per pair.

### H3 — `raised` is a run-time loop of multiplications, called six times per allele per pair

`crates/popnei/src/stats.rs:666`. Hot-path, confidence high.

The exponent is a `u32` field, so the trip count is unknown at compile
time and the loop sits inside the allele loop with a compare and a branch
per multiply. Almost every dataset popnei reads is diploid. An arm `2 =>
value * value` is bit-identical: the loop starts at `1.0` and multiplying
by `1.0` is exact. It must be a literal multiplication and not `powi`:
with a constant exponent LLVM regroups, `(v*v)*(v*v)` for 4 where the loop
gives `((v*v)*v)*v`, and those differ in the last place.

### H4 — The genotypes of a population are walked twice at every variant

`crates/popnei/src/pop_dists.rs:121`. Hot-path, confidence high.

`count_alleles_of` and `count_gts_of` each loop over the population's
individuals calling `genotype_of`, which is a checked multiply, a checked
add and a bounds-checked slice. At 20 populations that is 2000 of those
calls a variant, 2e8 a pass, and the row is read twice for each
population. One walk that fills both the allele counts and the genotype
counts reads each genotype once. It is integer counting in the same order,
so no number moves and the first refused allele is still the first.

Measured by: the 0.268 s per-individual term of section 4, which is where
this is charged, and which the 3-population pass must move too.

### L1 — An error is built and dropped on every call of `num_individuals_of`

`crates/popnei/src/variant.rs:405`. Likely, confidence high.

`.ok_or(Error::GtsNotWholeGenotypes { .. })` takes its argument by value,
so the variant is written to a stack slot and dropped through the
out-of-line glue `Error` has because other variants own strings. Both
counting functions call it, so 4e6 times a pass at 20 populations, and the
profile shows `drop_glue::<popnei::error::Error>` at 0.6% there.
`ok_or_else` moves it into the arm no valid block reaches. It is the same
finding as H1 of the Kosman review, applied there at 7d29ae4, and it is
the only eager `ok_or(Error::` left on this path. Expect no more than the
0.6% the profile shows.

### L2 — The standard errors compute every pseudo-value twice

`crates/popnei/src/pop_dists.rs:1162`. Likely, confidence high, cold
today.

`standard_error_of` walks the groups twice and calls
`of_the_group_left_out` in both loops, so each group's leave-one-out value
is computed twice for every measure and every pair: 4200 times at 3
populations and 266000 at 20. That is the 8 to 14 ms of section 7 and 2.8
in 100 of a 20-population call, which is not worth a scratch vector. It
stops being cold under `jackknife_group="variant"`, where the same count
is 266 million over 100000 variants.

### L3 — The split of a block is sized from the thread count, on cores that are not equal

`crates/popnei/src/pop_dists.rs:1598`. Likely, confidence medium.

rayon's default splitter makes about as many leaves as there are threads,
each of several chunks, and halves further only when a task is stolen. A
leaf started on an efficiency core near the end of a block holds the
per-block barrier. `.with_max_len(1)` makes every chunk of 64 rows its own
task so a performance core can steal the tail. It changes only the split
depth: the chunks, their contents and the order they are added in are
untouched, so the result is bit for bit what it is today at every thread
count and in wasm. A chunk size derived from `rayon::current_num_threads()`
must never be written, because it would make the sums depend on the pool
and break the test that asserts identical bits across pools of 1 and 4.

### L4 — The rows are cut into groups on the calling thread while the pool waits

`crates/popnei/src/pop_dists.rs:1499`. Likely, confidence medium.

`cut_the_rows_into_groups` builds a view of every variant of the block,
which slices 2000 genotypes it never reads, to take the chromosome and the
position, and it runs between two fork-joins with 18 workers asleep; the
reduction of the chunks afterwards is serial too. A sampling profile
cannot see a sleeping worker, so this needs three `Instant`s around the
cut, the collect and the reduce, printed from the benchmark at 1 and at 18
threads, before anything is built. It is refuted if the cut and the reduce
together are under 3 ms at 18 threads at 20 populations.

### N1 — Nothing after the clock consumes the seven measures

`crates/popnei/benches/pop_dists.rs:556`. Note.

The benchmark reads the pairs, the variants and F_ST of the first pair
after the clock stops. Today the pass is a cross-crate call that is not
inlined over heap data, so nothing can be folded away, but that argument
disappears under fat linking, and the reader's own review already warned
that "code that got faster under LTO could be measuring the parse being
deleted". Before any build-flag experiment, fold every measure of every
pair into one number, print it, and gate on it being unchanged.

### N2 — The benchmark's two-allele, always-counted data

`crates/popnei/benches/pop_dists.rs`. Note.

Every variant of both settings has two alleles and every population is
above `min_num_individuals` at every variant, so the two branches that
depend on the data, the allele that is not 0, 1 or missing and the variant
that does not count for a pair, always take the same side. The per-pair
coefficient of section 4 is therefore a two-allele number. A `--alleles`
on the generator, timed at 2 and at 8, is what would state it with its
allele count.

## 7. What the 8 to 14 ms between Python and the benchmark is

It is not the boundary. The benchmark stops its clock at
`calc_pop_dist_sums`, before any measure or standard error is read; the
Python script times the pass, then the extraction of all seven measures
with their standard errors, then the crossing. The residual grows from 8
to 14 ms as the populations go from 3 to 20, a factor that follows the
266000 leave-one-out evaluations of L2 and not the crossing, which is 1000
names in and about 16 arrays out whatever the populations. The interpreter
is released for the whole pass and the whole extraction, and every result
array is moved into numpy rather than copied.

So a number for this calculation has to say which of the two it is. The
benchmark's number is the pass; the Python number is what a user waits
for, less the reader.

## 8. The reader, which is now the larger half

Over the vars file on 18 threads at 3 populations the reader is 0.111 s of
the 0.134 s a user waits, and at 20 populations 0.111 s of 0.155 s. It
runs on the calling thread and does not scale: it is 0.106 to 0.112 s in
every one of the eight timings. The profile's `lz4_flex` frame is the same
3050 and 3063 samples in two profiles of different totals, which is what a
fixed serial cost looks like.

Reading a block while the calculation works on the one before it is the
read-ahead thread that `docs/specs/block.md` leaves for the first
calculation that consumes blocks, and `docs/specs/dists.md` already says
it is an item of that spec and that nothing here changes with it. What it
could give, by subtracting measured parts and not by experiment: at 18
threads the wall can never fall below the reader's own 0.111 s, so 0.134 s
to about 0.111 s at 3 populations and 0.155 s to about 0.115 s at 20; and
nothing if lz4 and the counting contend for memory bandwidth, which both
use. On one thread it would spend a second thread where the user asked for
one. It changes none of the numbers this review proposes, which take the
reader out.

Asking the reader for the chromosome and the position as well as the
genotypes, which the resampling groups need, costs it the decompression of
one 800 KB column against 200 MB of genotypes, and `block_of_the_batch` is
0.1% of the profile. There is nothing smaller to ask for: the groups are
cut from the position of every variant.

## 9. The build configuration

`[profile.release]` sets `overflow-checks = false` and `debug =
"line-tables-only"` and nothing else: optimization level 3, no linking
across crates, 16 code generation units, no `target-cpu`.
`[profile.bench]` inherits it. Each of the four was looked at and none is
proposed as an experiment now:

- **Linking across crates and one code generation unit** cannot move a
  float here. Contracting a multiply and an add into a fused instruction
  needs a fast-math flag that rustc does not emit for ordinary `f64`, and
  reassociating a float sum needs the same, so the accumulation order
  survives all of it; aarch64 has no excess precision, so an inlining
  decision cannot change a rounding. They were measured on the statistics
  pass, which shares both counting functions with this module, and gave
  0.182, 0.183 and 0.183 s against 0.182 s at the default, at 2.5 times
  the rebuild. There is a mechanism here that was not there — the three
  counting functions have frames of their own, so they are real
  out-of-line calls across a crate boundary — but H1 and H4 change those
  functions, so this is measured after them and gated on N1.
- **`target-cpu`** can be closed rather than measured: the default for
  this target and `native` differ by three features, none of which LLVM
  emits from integer counting loops, and it cannot ship in a wheel.
- **An allocator swap** is not a candidate at about 4700 allocations a
  pass.
- **`+simd128` for wasm** is already on, set in `.cargo/config.toml` for
  another module, so the wasm numbers of section 2 were measured with it.

## 10. The numbers to reach, and the decision they need

"Speed" of `docs/specs/dists.md` is waiting for numbers. What this review
proposes to write there, once the experiments of section 11 are run, is
the measured one-thread and wasm numbers with a tenth over them, which is
the form the Kosman item of the same spec uses.

The decision. The Kosman item states a target on 18 cores. On this machine
that kind of figure is not reproducible to better than a factor of two:
the same binary on the same dataset gave 0.027 s and 0.061 s for the
3-population in-memory timing in two invocations, and 0.044 s and 0.092 s
at 20 populations, while every one-thread figure over the same runs moved
by under 6 in 100. The owner's own account of why is that the 18 cores are
6 performance and 12 efficiency, so which kind a chunk lands on depends on
what else is running. The options:

- **One-thread and wasm numbers only**, with the parallel figure recorded
  in this report and not in the spec. The spec then holds only numbers
  that can be checked again on a machine that is not quiet. It gives up a
  stated goal for the threads, which is where goal 4 of the objectives
  cares most.
- **All three, with the parallel one as a range** and the pool size and
  the load average beside it, as the Kosman item does with a single
  number. It keeps a goal for the threads and it will be argued about
  every time somebody checks it.
- **All three, with the parallel one taken at `RAYON_NUM_THREADS=6`**, the
  six performance cores, which is also the pool pyNei was measured at.
  That is reproducible in a way the 18-thread figure is not, and it states
  a goal for the threads. It needs one measurement this review has not
  taken.

Recommended: the third. It is the only one that keeps a checkable number
for the threads, and what it costs is one run of the benchmark.

## 11. The measurement plan

In the order in which each unblocks the next. Every wall time is one
thread, best of 5, on a machine checked quiet before and after, since that
is the figure this machine reproduces.

1. **H1**, the counter table. Gate on the bytes zeroed per pass, which is
   the same on every run, then on the per-population coefficient of
   section 4 falling, then on `_platform_memset` leaving the profile. It
   is first because it is the largest measured cost and because it yields
   the largest allele without a scan, which H2 wants.
2. **H2**, the per-population arithmetic. Gate on the count of `fdiv` in
   `sums_of_the_chunk` from `objdump -d`, then on the per-pair
   coefficient falling.
3. **H3** and **L1**, one line each, measured together only if each is
   confirmed by its own count first: the inner loop gone from the
   assembly, and the destructor call gone.
4. **H4**, the fused walk, which is charged per individual and so must
   move the 3-population pass too.
5. **L4**, the three `Instant`s, which is instrumentation and not a
   change, before anything is built for it.
6. **L3**, `.with_max_len(1)`, with a sweep of `RAYON_NUM_THREADS` over 1,
   2, 4, 6, 8, 12 and 18 first: a curve that is straight to 6 and flat
   after says the efficiency cores are the whole story.
7. The **build flags**, last, gated on N1.

Every one of them is checked before its timing is read with
`cargo test --workspace`, which holds the test that the size of the blocks
does not change a measure beyond 1e-12 relative and the test that the
number of threads does not change one bit; and, for H1 to H4, with the 117
numbers of the test panel compared as bits before and after. That
comparison must include the standard errors and not only the seven
measures: a leave-one-out sum is formed by subtracting a group from the
total and multiplied by a weight of about the number of groups, so a
last-place change reaches the standard error magnified about a
hundredfold at 100 groups.

## 12. Seen outside the scope

- `crates/popnei/src/pop_dists.rs:1616`. The path that answers a failed
  allocation allocates: when a thread returns an error the block is read
  again into `vec![PairSums::default(); of_each_group.len()]`, a second
  full copy of the pass sums, and one of the errors that reaches it is the
  one `grow_the_sums` raises through `try_reserve` precisely so that the
  process is not ended. 912 KB at 20 populations and 100 groups, 912 MB
  with a group of each variant. No wrong number; for a code review.
- `js/popnei/src/arguments.ts:408` spreads a whole population onto the
  argument stack, which throws in V8 at about 100000 names. Far above the
  10000 individuals of the objectives.
- `crates/popnei/src/pop_dists.rs:1096`, `f2_of_every_group` collects a
  `flat_map`, whose lower size hint is 0, so it grows by doubling for
  19000 values. Microseconds against a pass; for a later reader.

## 13. What this code already does well

- **The reduction order is written down and tested, not left to rayon.**
  Each chunk of 64 rows sums into its own run of pairs and the chunks are
  added in block order, so the result does not depend on the pool. The
  test at `pop_dists.rs:4107` asserts the bits are identical across pools
  of 1 and 4 and records the counter-example under rayon's own `reduce`:
  F_ST of the first pair is `3fbaded19c839733` on one thread and
  `3fbaded19c839729` on four. Every finding above had to be argued against
  that test, and it is what made the arguments checkable.
- **The wasm build is the same chunks in the same order**, so the two
  builds agree bit for bit by construction and not by luck:
  `add_the_chunks_one_by_one` is both the wasm body and the fallback the
  native body uses to find the first bad row.
- **Every division happens once, at the end.** The six sums a pair keeps
  are values of single variants, never ratios, which is what lets the
  blocks and the threads add them in any order and is why H1 and H2 can be
  bit-identical at all.
