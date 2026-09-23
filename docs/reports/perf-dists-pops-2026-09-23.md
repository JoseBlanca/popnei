# Performance review: the distances between populations

23 September 2026. The first measurement of `calc_pop_dists`, the
calculation merged into `main` that day in 05f1538, which gives for every
pair of populations, out of one reading of the variants, seven measures
with a block jackknife standard error each. Nothing of it had ever been
timed: no benchmark, no profile, no wall time, for the core, for Python or
for wasm. The owner decided on 23 September 2026 that the plan which built
it would not measure and that this review would, and "Speed" of
`docs/specs/dists.md` asks this review to write the numbers to reach into
it once they exist.

Section 2 says what the code does now and what it did before, section 3
how it was measured, section 4 how the cost grows with the populations,
which the spec asked to have measured rather than guessed, section 10 the
numbers now written into the spec, and section 14 the three changes this
review made and what each gave. The branch `perf/dists-pops` holds them
and is not merged; the order to merge is what this report asks for.

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

**Applied.** Three changes were built, measured and kept, on the branch
`perf/dists-pops`, which is not merged. Together they make the pass 39 to
49 in 100 faster on one thread, and **no number of any result moved by a
single bit**: 3807 numbers of the test panel were compared as bits at each
step, under three ways of cutting the resampling groups, including the one
that magnifies a last-place change most. Section 14 gives each change with
its two numbers.

On 100000 variants x 1000 individuals with the reader taken out, on one
thread, best of 5 on a machine checked quiet:

| | 3 populations, 3 pairs | 20 populations, 190 pairs |
|---|---|---|
| popnei before, at 3bc33f9 | 0.282 s | 0.489 s |
| **popnei now** | **0.171 s** | **0.248 s** |
| popnei in wasm before | 0.375 s | 0.582 s |
| **popnei in wasm now** | **0.198 s** | **0.320 s** |
| pyNei | 7.811 s | 502.105 s |

Over the vars file, which is what a user waits for and has the reader in
it as pyNei's number has pyNei's, popnei now takes 0.277 s and 0.352 s on
one thread, so it is **28 times faster than pyNei at 3 populations and
1426 times at 20**. The second ratio is mostly pyNei's own shape and not
popnei's speed, and section 5 says why.

The numbers to reach are written into "Speed" of `docs/specs/dists.md`,
where this review was asked to put them: the measurements above with a
tenth over them.

What is asked of the owner: **the order to merge this branch**, which
nothing here stands in the way of. One thing to know before deciding, and
not a decision: H2 is a gain from about 10 populations up and is neutral
below that, because it moves work from the pairs to the populations;
section 14 has both numbers.

The verdict before the experiments were run was "run the experiments", on
this evidence: two costs that the profile named and the cost model
measured, the counter table walked in full for every population at every
variant (10.2 GB a pass at 20 populations against 200 MB of genotypes
read), and each population's own arithmetic redone for each of the 19
pairs it is in. Both were confirmed.

## 3. What was measured, and with what

Every number in this section is of the code as it was merged, at
3bc33f9, before the three changes of section 14. It is kept because it
is what the findings were made from and what the next review will
compare against.

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

## 10. The numbers written into the spec

"Speed" of `docs/specs/dists.md` was waiting for numbers and now has them:
the measurements of section 14 with a tenth over them, of the calculation
with the reading taken out, on the owner's M5 Pro.

| | 3 populations | 20 populations |
|---|---|---|
| one thread | 0.19 s | 0.27 s |
| 18 threads, on a quiet machine | 0.023 s | 0.028 s |
| wasm under node | 0.22 s | 0.35 s |

The form is the one the Kosman item of the same spec uses. What is new is
the condition on the middle row and a sentence saying that the one-thread
figure is the one to check a change against.

Whether to state a figure for many threads at all was open while this
review ran, because the same binary gave 0.027 s and 0.061 s for the same
18-thread timing in two invocations. It is settled by measurement rather
than by judgement: on a machine checked quiet the 18-thread figure repeats
within 5 in 100 (0.025, 0.025 and 0.026 s over five runs at 20
populations), and the factor of two appeared only when another session was
compiling or a profiler was attached. So the figure is worth stating with
the condition attached, and it is stated.

The third option this review considered, taking the parallel figure on the
6 performance cores alone because they are alike, was measured and
dropped: `RAYON_NUM_THREADS=6` gives 0.031 s at 3 populations and 0.047 s
at 20 against 0.021 s and 0.025 s on all 18. The 12 efficiency cores are
worth 1.9 times at 20 populations, so a target that left them out would
ask for less than the machine gives.

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

## 14. What the experiments gave

Each change was built on the one before it, so each baseline is the commit
before it and not the state the review started from. Every timing is one
thread, best of 5, blocks in memory, on a machine checked quiet with
`ps` and `sysctl -n vm.loadavg` before and after; the machine was shared
with other sessions throughout and two of the three experiments had to
wait for another session's build to finish before their final timing.

Before anything was timed, each change passed `cargo test --workspace`,
`cargo clippy --workspace --all-targets -- -D warnings`, `cargo fmt --all
--check`, `cargo wasm-check`, `uv run ruff check` and `uv run pytest`. And
each proved the bits directly: every number of the test panel printed as
`f64::to_bits` before and after and compared — 117 numbers under 24
resampling groups in blocks of 100 variants, 45 under no groups in one
block, and 3645 under one group for each variant in blocks of 7, which is
the cut that magnifies a last-place change most. **All 3807 were identical
at every step.** The test that asserts the same bits across pools of 1 and
4 threads passed at each step too, so nothing was made to depend on the
pool.

### H1, applied at 9d828e8: count a population of two alleles without the table

`count_alleles_of` now counts a population whose alleles are only 0, 1 and
the missing one in four `u32` counters and never touches the 128-entry
table, falling back to the table at the first individual it cannot count
that way; and the largest allele it saw comes out of the counting, so the
backwards scan of 128 counters in `count_the_var` is gone.

The gate was a count, not a time: the bytes zeroed per pass went from
**4 096 000 000 to 0**, the first being exactly 100000 variants x 20
populations x 2048 bytes. The per-population coefficient fell by 75 in
100. The 20-population pass went from 0.501 s to 0.363 s and the
3-population pass from 0.286 s to 0.238 s. `_platform_memset` left the
profile, from 1585 samples to 4.

The per-individual coefficient also fell by 13 in 100, which was not
predicted: the fast path replaces a bounds-checked indexed increment into
one of four lane arrays with three compares and three adds in registers,
so the counting of each allele got cheaper and not only the table.

Cost: one new public two-field struct, two counting paths where there was
one, and two `#[expect]` of the arithmetic lint with the bounds the code
already establishes. Two tests were added and both were shown to fail
against deliberately broken code.

### H2, applied at 4d6b616: take each population's own quantities once a variant, not once a pair

`PopVarCounts` now carries, beside its counts, each population's allele
frequencies and the five numbers a pair reads of one of its two
populations. The pair loop keeps only what is of the pair: the product of
the two frequencies, its square root and the pooled power.

The gate was the divisions executed for one biallelic variant at 20
populations, attributed to their innermost loop in the disassembly of the
benchmark that was actually built: **2850 before and 1050 after, a fall of
63 in 100**. The per-pair coefficient fell by 46 in 100. The 20-population
pass went from 0.370 s to 0.326 s and the 40-population pass from 0.714 s
to 0.526 s. `sums_of_the_chunk` fell from about a third of the work to a
seventh.

The condition that keeps it bit-identical was honoured: the hoisted values
sit behind the same integer tests that made `of_var` give no value, so a
population with fewer than two called alleles still yields nothing rather
than an infinity.

Cost, and the one thing the owner should know: **below about 10
populations this change is neutral**, because the per-population
coefficient rises by 23 in 100 to pay for the per-pair fall. At 6
populations the two cancel exactly. `PopVarCounts` grows from about 540
bytes to about 1610 and now holds derived floats that belong to the
variant just counted and to the settings it was counted with.

### H4, applied at 7df1472: count the alleles and the genotypes in one walk — and the finding as written was wrong

The experiment refuted its own stated mechanism and found the gain
elsewhere, which is worth recording.

H4 said the cost was the two walks over the same individuals and the 2e8
`genotype_of` lookups a pass. Built that way — one walk, the genotype
taken once, the existing `count_the_genotype` called per genotype — it was
**7 in 100 slower** than the baseline at every population count. Halving
the lookups bought nothing.

What paid was not re-reading the genotype's alleles: counting the genotype
out of the zeros, ones and missing alleles the fast path has already
counted, instead of walking its alleles again. And one more thing, which
is the most transferable result of this review: the version that still had
`if missing { .. } else { .. }` per genotype was faster on every
populations file **except** `pops3.tsv`, which it made slower than the
baseline, 0.245 s against 0.237 s. `pops3.tsv` is the only file whose
populations are scattered through the row, in 181 runs of varying length,
because the three populations are the simulated ones. Writing the three
counts as `+= u32::from(..)` with `&` rather than `&&`, so that nothing
branches on the genotype, took that file to 0.167 s. **An unpredictable
branch in a loop whose loads are irregular cost a quarter of the pass.**

The per-individual coefficient fell by 31 in 100 and the 3-population pass
from 0.237 s to 0.167 s, which is the control that says the change moved
what it claimed: it is charged per individual, so the 3-population pass
had to move too, and the 150-individual panel moved with it.

Cost: three counting functions where there were two, and about 45 lines
repeated between two of them that nothing but the tests would catch
diverging. What is not repeated is what matters most — the rules for which
of the three counts a genotype falls in now live in one function that both
paths call, pinned by a test over 4452 populations covering every genotype
of one and two individuals at ploidies 1, 2 and 3.

### Not run, and why

- **H3**, the arm for ploidy 2 in `raised`. H2 cut its call count from
  1140 a variant at 20 populations to 420, so the most it can now give is
  under 1 in 100 of the pass. The plan in section 11 stands and its count
  should be retaken before it is built.
- **L1**, the error built eagerly in `num_individuals_of`. H4 made that
  function run once per population per variant instead of twice, so the
  0.6 in 100 the profile showed is now about 0.3. One word, and worth
  taking the next time that file is opened.
- **L2**, **L3**, **L4** and the build flags. L2 is cold until somebody
  asks for a resampling group of each variant. L3 and L4 are about the
  threads, and the one-thread pass is what this machine measures
  reliably; L4 is instrumentation first and nothing was built for it. The
  build flags were left for the reason section 9 gives, and they should
  now be measured against N1's checksum rather than against a timing,
  since the pass is smaller than it was.

### What the numbers look like now

The cost model, refitted on the same five cuts of the same 1000
individuals, each point within 2 in 100:

| | before, at 3bc33f9 | now |
|---|---|---|
| per individual counted, per variant | 2.68 ns | 1.64 ns |
| per population, per variant | 60.2 ns | 14.1 ns |
| per pair, per variant | 5.80 ns | 2.82 ns |

The model is now less able to extrapolate than it was: fitted on panels of
1000 individuals it under-predicts the 150-individual panel by 22 in 100,
where before it missed by 9. What that says is that a cost per variant
that does not depend on the individuals, about 8 ms a pass, was hidden
under the counting and is now visible. Nothing in this review names it.

The reader is now most of what a user waits for over a file. At 20
populations on 18 threads the pass is 0.025 s and the reader 0.110 s of
the 0.135 s; at 3 populations it is 0.021 s against 0.112 s of 0.133 s.
Section 8 says what a read-ahead thread could give and that it belongs to
`docs/specs/block.md`.
