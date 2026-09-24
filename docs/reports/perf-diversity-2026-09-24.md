# Performance review: the diversity module, 24 September 2026

The one pass of popnei's diversity module takes 0.604 s on one thread and
0.168 s on 18 cores over the 100000 variants of 1000 diploid individuals
of `/Users/jose/devel/popnei-bench/big.vars` in three populations, with
every statistic asked for and at a draw of 200 called alleles, and 0.299 s
and 0.140 s with the four statistics that need no draw. Task 4.1 of
`docs/plans/diversity.md` measured it so that "Speed" of
`docs/specs/diversity.md`, which had no measurement and no target, could
get one, and so that the one question that section left open could be
answered by a measurement: whether the hypergeometric weights of a variant
are computed once and shared by the bins of its spectrum. They are, and
the measurement says what that sharing is worth: without it the same pass
takes 11.447 s, 22.9 times as long. This document holds every number, the
two timings of that experiment, and the one thing it found that belongs to
a review of its own.

The words this document uses. A **pass** is one whole read of a dataset
with whatever is calculated over the variants on the way; the passes timed
from Python are timed from the call that opens the file to the result. A
**block** is the genotypes of many variants as one flat array of bytes,
one row per variant, a row being the individuals one after another, each
taking as many bytes as the ploidy; `docs/architecture.md` puts one at
about 10 MB. A **vars file** is popnei's own format, arrow IPC with the
genotypes compressed with lz4. A **population** is a named set of
individuals that each statistic is calculated over on its own. The **draw**
is the argument `num_called_alleles`: every population is brought down to
that many called alleles at each variant and the expectation is taken over
every such draw, which is called rarefaction for a count of alleles and
projection for the spectrum. The **folded site frequency spectrum** of a
population is how its variants are spread over the count of the rarer
allele in that draw, a vector of `num_called_alleles / 2 + 1` **bins**; the
chance that a draw shows a given count is the **weight** of that count,
and the weights of one variant are what deliverable 2 of the plan asks
about. The **row loop** is `add_the_rows` of
`crates/popnei/src/diversity.rs`, which walks the rows of one block and
adds what it counts into the sums of the pass; rayon splits a block into
**chunks** of 64 rows, and the reduction reads them in **groups**, which
`chunks_of_a_group` makes two chunks for each thread of the pool.
`calc_per_var_distribs` of `docs/specs/stats.md` is the pass this module's
target is a ratio against: over the same file and the same populations it
makes the same read and the same counts of how often each population
called each allele, and then does far less arithmetic on them.

## The scope and its limits

Reviewed at the commit 47ac373 of the branch `plan/diversity`:
`crates/popnei/src/diversity.rs`, which holds the whole pass, and the two
benchmark harnesses that commit added,
`crates/popnei/benches/time_diversity.py`, which times a whole pass from
Python over a file, and `crates/popnei/benches/diversity_pass.rs`, which
times the same pass over blocks that are already in memory and counts the
bytes it holds live.

No `perf-reviewer` subagents were sent. This is not a review that begins
with reviewers looking for candidates: the plan asked for one measurement
and one experiment, and named them both. What that leaves out is the whole
of what a set of reviewers would have looked at, the allocations, the data
layout, the concurrency and the Python boundary among them, and the
findings below are only what the profiles and the wall times put in front
of this session.

What was not measured. Nothing was run in a browser or under node, and
nothing in wasm: task 4.2 of the plan builds those and takes no time.
Nothing was run over a VCF; every figure from Python is over a vars file
already in the page cache. pyNei has none of these five statistics, as the
opening of `docs/specs/diversity.md` says, so there is no comparison with
it and none with any other library. One shape of dataset carries every
timing, 1000 diploid individuals in 3 populations with 3 in 100 genotypes
missing whole; the memory was taken at four shapes and the time at one. And
this machine has no hardware counters, so no cache miss and no mispredicted
branch was counted anywhere below.

## The verdict

Run no experiment beyond the one the plan named, which was run and is
closed. The folded spectrum, which "Speed" of the spec called the one part
that could dominate, is 20 to 23 in 100 of the pass at a draw of 200 and 8
to 10 in 100 at a draw of 20, so the bound that "What could go wrong" of
work package 4 sets, that anything beyond the shared weights of deliverable
2 is a performance review of its own, is not reached. The pass meets the
target that is now in "Speed" with headroom at every thread count: at most
1.1 times `calc_per_var_distribs` where it is 0.78 to 1.02 times, and at
most 2 times it at a draw of 200 where it is 1.16 to 1.61 times.

One thing is left for a later review and is under "Seen outside the scope":
three of the statistics compute the same products over and over on a
variant of two alleles, and they are 0.208 s of a 0.505 s pass together.

## What a whole pass costs

Every figure of this section was taken on the owner's Apple M5 Pro, 18
cores, 64 GB, macOS 27.0, native `aarch64-apple-darwin`, rustc and cargo
1.98.0, with popnei at 47ac373 built by `uv run maturin develop --release`,
on 24 September 2026. Each is the best of 5 timed runs with one untimed run
before them, timed in Python by
`crates/popnei/benches/time_diversity.py` from the call that opens the file
to the result, and each cell holds two sets taken one after the other. The
two sets agree within 0.006 s on every row. The threads are rayon's, which
popnei takes from the environment because it builds no pool of its own, so
over a vars file, whose reader runs on the thread that calls it, the
threads are the calculation's alone.

The dataset is `/Users/jose/devel/popnei-bench/big.vars`, popnei's vars
file of 100000 variants of 1000 diploid individuals with 3 in 100 genotypes
missing whole, the `big.vars` of
`docs/reports/filters-measurement.md`, read in the three populations of
296, 356 and 348 individuals of `/Users/jose/devel/popnei-bench/pops3.tsv`.

| the pass | the draw | 1 thread | 18 cores |
|---|---|---|---|
| the read alone, the genotypes as the only field | — | 0.110, 0.112 s | 0.111, 0.108 s |
| `calc_pop_diversity`, the four that need no draw | — | 0.299, 0.302 s | 0.140, 0.141 s |
| `calc_pop_diversity`, all five | 20 | 0.339, 0.340 s | 0.145, 0.145 s |
| `calc_pop_diversity`, all five | 200 | 0.604, 0.595 s | 0.168, 0.166 s |
| `calc_pop_diversity`, the folded spectrum alone | 20 | 0.262, 0.258 s | 0.137, 0.133 s |
| `calc_pop_diversity`, the folded spectrum alone | 200 | 0.355, 0.348 s | 0.147, 0.142 s |
| `calc_per_var_distribs`, its five statistics | — | 0.375, 0.383 s | 0.142, 0.143 s |

The orchestrator of the plan re-ran five of these rows itself and got
0.108 s for the read alone, 0.292 s for the four with no draw, 0.591 s for
all five at a draw of 200, 0.366 s for the spectrum alone at that draw and
0.381 s for `calc_per_var_distribs`, each on one thread and each inside the
spread of the two sets.

Against `calc_per_var_distribs` over the same file with the same
populations at the same threads, which is what the target of "Speed" is
stated as:

| `calc_pop_diversity` | the draw | 1 thread | 18 cores |
|---|---|---|---|
| the four that need no draw | — | 0.78 to 0.81 times | 0.98 to 0.99 times |
| all five | 20 | 0.88 to 0.91 times | 1.01 to 1.02 times |
| all five | 200 | 1.55 to 1.61 times | 1.16 to 1.18 times |

The panel, `tests/reference/stats/panel.vcf.gz`, 1200 variants of 200
diploid individuals in the three populations of 48, 68 and 84 of
`tests/reference/stats/panel_pops_bcftools.txt`, separates none of these
passes: every one of them, the read alone and all five at either draw, is
0.004 s on one thread and 0.002 s on 18, and `calc_per_var_distribs` is
0.005 s and 0.002 s, where the clock resolves 0.001 s. A draw of 200 on the
panel measures nothing about the draw besides. No population of the panel
calls 200 alleles at any variant, so no variant is in the draw for any of
them, every standardized value is NaN and all 101 bins are 0: that 0.004 s
is a pass that did the work of a pass with no draw.

## Where the time inside the pass is

Two sampling profiles with `sample`, on one thread, over `big.vars` at a
draw of 200, with the idle rayon workers filtered out as
`.claude/skills/performance-review/profiling_environment.md` says. Every
named function of the draw is inlined into
`popnei::diversity::add_the_rows`, and `sample` with a dSYM still names
only that, so the profile ranks the row loop against the read and against
the allele counts and goes no finer.

With every statistic asked for, 19365 working samples: `add_the_rows` 53.1
per cent of the self time, `variant::count_alleles_and_gts_of` 27.6, the
lz4 decompression of the vars reader 15.4, and nothing else above 1.3. With
the folded spectrum alone, 19207 samples: `count_alleles_of` 36.9,
`add_the_rows` 30.8, the lz4 decompression 26.7.

What the row loop spends it on was separated by differencing whole passes
of `crates/popnei/benches/diversity_pass.rs` instead, which runs the same
pass over blocks already in memory and so leaves the read of the file out.
Over 100000 variants of 1000 individuals in 3 populations on one thread,
best of 5 or of 9 runs. The floor, the copy of the block the harness makes
and the allele counts of every population, is 0.126 s, and what each pass
adds over it is the statistics it was asked for and nothing else.

| what the pass asks for | the draw | best | over the floor |
|---|---|---|---|
| `num_alleles` | 200 | 0.193 s | 0.067 s, `alleles_a_draw_shows` |
| `variable_vars_ratio` | 200 | 0.193 s | 0.066 s, `chance_a_draw_varies` |
| `private_alleles` | 200 | 0.193 s | 0.066 s, the standardized private alleles |
| `fis` | none | 0.174 s | 0.048 s, the two heterozygosities |
| `folded_sfs` | 20 | 0.148 s | 0.022 s, `add_the_bins_of_the_var` |
| `folded_sfs` | 200 | 0.234, 0.236 s | 0.108 s, `add_the_bins_of_the_var` |
| the four that need a draw | 200 | 0.371, 0.393, 0.396 s | |
| all five | 200 | 0.480, 0.490, 0.494 s | |
| all five | 20 | 0.218, 0.222 s | |

So the folded spectrum is 0.098 to 0.109 s of a pass of 0.480 to 0.494 s
at a draw of 200, 20 to 23 in 100 of it, and 8 to 10 in 100 at a draw of
20. The three standardized values, the number of alleles, the ratio of
variable variants and the private alleles, are 0.198 s together, 1.8 times
the spectrum. That is what refutes the sentence "Speed" had before the
measurement, that the spectrum is the one part that can dominate: it is the
smaller half of what the draw costs.

The estimate that sentence rested on, 6·10⁸ hypergeometric terms for a
million variants in three populations at a draw of 200, assumed that every
term is a product of about 200 factors.
`OfTheDraw::add_the_bins_of_the_var` does not compute them that way, and
the experiment below is what that is worth.

## The memory

The most bytes the pass holds live at once beside its block, counted by a
global allocator inside `crates/popnei/benches/diversity_pass.rs`, which
adds the size of every allocation and subtracts the size of every free and
keeps the largest total it saw. It gives the same figure to the byte on
every run. A megabyte here is 1000000 bytes.

| individuals | rows in a block | populations | the draw | 18 threads | 1 thread |
|---|---|---|---|---|---|
| 200 | 10000, the default for 200 individuals | 3 | 20 | 0.053 MB | 0.005 MB |
| 500 | 10000 | 50 | 180 | 2.007 MB | 0.155 MB |
| 1000 | 5000 | 50 | 2000 | 15.479 MB | 1.251 MB |
| 10000 | 500 | 50 | 20000 | 36.343 MB | 12.123 MB |

It does not grow with the variants. It grows with the populations, with the
bins of the spectrum and with the threads of the pool, because the
reduction holds one set of partial sums for each chunk of a group and a
group is two chunks for each thread. So every memory figure of this module
is stated with the threads it was taken at: the first shape holds 0.053 MB
on 18 threads and 0.005 MB on one, over the same variants.

The fourth shape is where the grouping stops binding. Its block is 500
rows, which is 8 chunks of 64, and a group on 18 threads is 36 chunks, so
every chunk of the block is read at once: the 36.343 MB is 8 sets of
partial sums of 4.0 MB of bins each, 50 populations of 10001 bins of 8
bytes, plus the pass's own 4.0 MB of the same shape. On one thread, where a
group is 2 chunks, the same dataset holds 12.123 MB. Read against the block
of about 10 MB of section 2 of `docs/architecture.md`, the first three
shapes are a small fraction of a block and the fourth, whose own block is
10 MB, holds 3.6 times one beside it.

Whether reading fewer chunks of a block than the pool has threads would be
better is a trade nobody has measured: it would hold fewer sets of partial
sums and would leave threads with nothing to do at the end of a block. No
finding is filed for it.

These four figures correct the ones the `architecture` reviewer of work
package 3 reported on 24 September 2026, which were 0.045, 2.342, 15.871
and 36.481 MB for the same four shapes. The measurement above is 0.008 MB
above that reviewer's on the first shape and 0.138 to 0.392 MB below it on
the other three. How that reviewer defined what it subtracted from the peak
is not known, so what the difference is made of cannot be said; what can is
that the figures above come from a harness that is in the repository, that
prints the block it subtracted, and that gives the same number on every
run. The orchestrator of the plan reproduced the first row itself at
0.053 MB.

## The measurement plan

Both harnesses are in the repository and every number above can be taken
again. The commands, as they were run, from the worktree of the branch:

Whole passes from Python, after `uv run maturin develop --release`, with
`uptime` read right before each invocation, from `crates/popnei/benches`:

    RAYON_NUM_THREADS=1 uv run python time_diversity.py \
        /Users/jose/devel/popnei-bench/big.vars diversity \
        /Users/jose/devel/popnei-bench/pops3.tsv 200 5

with `read`, `diversity`, `diversity-sfs` and `diversity-no-draw` as the
pass, `none` in the place of the 200 for the two that take no draw, and
`RAYON_NUM_THREADS=18` for the 18 cores. `time_stats.py` beside it times
`calc_per_var_distribs` over the same file and the same populations file,
which is the pass the target is a ratio against.

The pass over blocks in memory, and the bytes it holds:

    cargo bench --bench diversity_pass -- --stats all --draw 200 --pops 3

with `--stats` taking the names a user writes, `all`, `without_a_draw`, or
any of `num_alleles`, `private_alleles`, `variable_vars_ratio`,
`folded_sfs` and `fis` separated by commas, and with `--threads`,
`--individuals`, `--vars` and `--blocks` for the other shapes. The memory
rows above are its `--threads 18` and `--threads 1` figures at those four
shapes.

What the next review of this module would run first, in the order in which
it unblocks the one finding that is open, which is under "Seen outside the
scope": the three standardized values in one pass,
`--stats num_alleles,variable_vars_ratio,private_alleles --draw 200
--pops 3` on one thread, which reads 0.340 s against a floor of 0.132 s;
then the same with a shared table of the chance that each population's draw
misses each allele of the row, which is worth keeping if it takes the
0.208 s those three add over the floor below about 0.15 s and leaves every
number of the spec inside its tolerance.

Those two figures were taken later on 24 September 2026 than the
differencing table above, in the session that ran the experiment, at load
averages of 2.2 and 2.6. The machine was slower then by about the same
amount everywhere: the floor read 0.132 s where the table has 0.126 s, and
the pass with all five at a draw of 200 read 0.500 to 0.505 s where the
table has 0.480 to 0.494 s. What each statistic adds over the floor is the
same in both, 0.067 s for `num_alleles` at a draw of 200 in each.

A sampling profile that names something finer than `add_the_rows` is what
this review could not get. Everything the row loop calls is inlined into
it, and `sample` with a dSYM still names only that, so the split of the row
loop above comes from differencing whole passes and not from a profile. A
review that wants the line would have to mark the named functions
`#[inline(never)]` for the profiling build alone, which changes the code
being measured and was not done here.

## The build configuration

The workspace sets `overflow-checks = false` and `debug =
"line-tables-only"` under `[profile.release]`, so a sampling profile has
function names, and `[profile.bench]` keeps the line tables too; the test
profile is `opt-level = 2` with the overflow checks on. Nothing sets
linking across crate boundaries, one code generation unit, a panic
strategy, an allocator or `target-cpu`, and this task tried none of them.

The performance review of the stats module closed two of those on 22
September 2026 with no gain on a pass over the same file,
`docs/reports/perf-stats-2026-09-22.md`: thin linking across crates and
compiling the crate as one piece each moved a pass of 0.183 s by at most
0.001 s, with the builds interleaved back to back. This module's pass runs
in the same crate over the same reader, so neither was run again here.
Nothing about `target-cpu` or the allocator has been measured for popnei at
all.

## The findings

One, which is the question the plan named, and it was decided by the
experiment below rather than left as a candidate. Nothing else in the scope
is proposed for change: the pass meets its target with headroom, the row
loop allocates nothing per row, and the only site where a measurable gain
is plausible is outside what this task may touch and is under "Seen outside
the scope".

### D1 The hypergeometric weights of a variant are computed once and shared by its bins

`crates/popnei/src/diversity.rs:983`, `OfTheDraw::add_the_bins_of_the_var`.
Hot path, high confidence: the function is 0.108 s of a 0.494 s pass at a
draw of 200 by the differencing above, and the profile puts the row loop it
is inlined into at 53.1 per cent of the self time.

Deliverable 2 of work package 4 asks whether the weights of a variant are
computed once and shared by the bins of its spectrum, which "Speed" of the
spec named as the thing to decide. They already are. Per variant and per
population the function takes one product, of `num_called_alleles` factors
where the copies of the major allele cannot fill the draw on their own and
of the called alleles less that number where they can, for the smallest
count of the range the variant can show, and then walks the bins upward
from it by a recurrence: the chance of one more copy of the rarer allele is
the chance in hand times `(m - j)(g - j) / ((j + 1)(c - m - g + j + 1))`,
one multiply and one divide for each bin.

So the experiment is the other way round, which the `performance-review`
skill allows under "confirms it by reverting the change and measuring
again": the unshared version was built as a probe, one that computes each
bin's weight from its own product of factors, and both were timed. What it
gave is below.

## What the experiment gave

One experiment, at 47ac373 plus the commit that wrote the target into the
spec, which changed no code. Every run is
`cargo bench --bench diversity_pass -- --draw 200 --pops 3` on one thread,
over 100000 variants of 1000 individuals in 3 populations, best of 5 timed
runs with one untimed run first, with `uptime` read right before each
invocation and nothing else building on the machine.

The probe replaces the recurrence with the identity `C(m, j) C(c - m, g -
j) / C(c, g) = C(g, j) (m)_j (c - m)_(g - j) / (c)_g`, with `(x)_k` the
falling factorial, and multiplies the three out for each bin on its own,
interleaved so that the product stays in the range the module's mantissa
and power of two hold. It is a product of up to `num_called_alleles`
factors per bin where the shared version has one per variant and per
population.

The correctness checks of the `coding` skill were run before any timing was
looked at, and the probe passed all of them: `cargo test --workspace` 868
passed with 2 ignored and 149 passed, `cargo test -p popnei
--no-default-features` 868 passed with 2 ignored, `cargo clippy --workspace
--all-targets -- -D warnings` with no warning, `cargo wasm-check` exit 0,
and `uv run maturin develop && uv run pytest` 530 passed, the comparison of
the panel's spectrum with `dadi`'s stored values inside 1e-12 relative
among them. So it is a slower version of the same calculation and not a
wrong one.

| the pass, a draw of 200 | shared, as it is | the unshared probe | the load average |
|---|---|---|---|
| the folded spectrum alone | best 0.243 s, median 0.244 s, worst 0.249 s | best 11.042 s, median 11.231 s, worst 11.312 s | 2.65, then 3.17 |
| all five statistics | best 0.500 s, median 0.510 s, worst 0.513 s | best 11.447 s, median 11.494 s, worst 11.581 s | 2.65, then 3.52 |

The unshared version is 45.4 times the shared one on the pass with the
spectrum alone and 22.9 times it on the pass with all five. The difference
is nowhere near the 5 per cent at which the `performance-review` skill asks
for a second clean measurement, so none was taken for the probe; the
shared version was measured once more after the probe was reverted, at a
load average of 2.75 and 2.61, and gave 0.505 s and 0.247 s, inside 0.005 s
of its baseline.

The probe was reverted and is not committed. The recurrence stays, and it
would stay whatever the timing had said, for a reason that is not about
speed: the review of work package 3 found that computing the bins from the
smallest term of the range and no recurrence underflowed to 0 and lost
whole variants out of the spectrum at a draw above about 560 called alleles
in a population. The measurement decides what the sharing is worth;
rightness decides that it stays.

One number the two versions do not share. The 101 bins of the first
population sum to 100000.00000000001 with the recurrence and to
100000.00000000007 with the probe, against the 100000 variants that are
exactly what they should sum to. Both are far inside the 1e-12 relative
tolerance the spec sets against `dadi`, and this is one sum on one dataset,
not a general claim about which version rounds better.

No test was added or changed by this task, and the probe was checked by the
suite that already exists. That suite was shown to have power over the
function the probe replaced: with the recurrence's divisor `j + 1` changed
to `j + 1.5` and nothing else, `cargo test -p popnei` gave 860 passed and 8
failed, the first of them
`diversity::the_pass::a_count_of_the_rarer_allele_above_the_draw_is_in_no_bin`
and among them
`diversity::the_pass::the_folded_spectrum_of_the_panel_is_the_one_dadi_gave`
and `a_draw_of_580_of_six_hundred_heterozygous_individuals_has_the_bins_of_the_exact_values`.
That break was reverted too.

## Seen outside the scope

**The three standardized values compute the same products over and over on
a variant of two alleles.** They are 0.208 s of a pass of 0.505 s together,
which is more than the spectrum, and the mechanism is arithmetic that is
thrown away rather than a cache or a branch.

Every standardized value of the module is built from one chance: that a
draw of `g` of the `c` copies a population called at a variant holds no
copy of an allele it called `n` times, `C(c - n, g) / C(c, g)`, which
`chance_a_draw_misses_an_allele` takes as a product of `g` factors. Three
places ask for it, and on a variant of two alleles of counts `n` and
`c - n` they ask for the same two numbers:

- `OfTheDraw::alleles_a_draw_shows` sums `1` minus that chance over the
  alleles the population called, so it computes the chance of missing the
  allele called `n` times and the chance of missing the one called `c - n`
  times.
- `OfTheDraw::chance_a_draw_varies` sums, over the same alleles, the chance
  that every copy the draw takes is that one allele, which
  `chance_a_draw_is_all_of_one_allele` computes as the chance of missing
  every other allele, that is `C(c - n, g) / C(c, g)` for the allele called
  `c - n` times and `C(n, g) / C(c, g)` for the one called `n` times. With
  two alleles those are the same two products as above, taken in the same
  factor order, so the values are not merely close, they are the same bits.
- `the_chance_each_draw_misses_the_allele`, which the standardized private
  alleles call for every allele of the row and every population, computes
  `chance_a_draw_misses_an_allele` for each of them, which on a biallelic
  variant is again those same two products.

So a pass asked for all three takes six products of `g` factors per variant
and per population where two would do. The measured cost agrees with that
picture: each of the three adds 0.066 to 0.067 s to a pass on its own, and
the three together add 0.208 s over the floor of 0.132 s, which is close to
the sum of the three, 0.199 s.

What it could be worth, as an upper bound and not a prediction: if the
three shared one table of the chance that each population's draw misses
each allele of the row, and if the products were the whole of what they
cost, the three would fall from 0.208 s to about 0.067 s over the floor,
which is 0.14 s of a 0.505 s pass, 28 in 100 of it. The sums and the loops
around the products are not free, so the real figure is smaller, and only
the experiment says by how much.

What would have to be settled first. The identity holds for a variant of
two alleles and not for one of three or more, where `c - n` for one allele
is no other allele's count, so the shared table serves
`chance_a_draw_varies` only on a biallelic variant and the function needs
its own path for the rest. The table is per population and per allele of
the row, so it is memory that grows with the populations, which the memory
section above says is already what this pass grows with. And the values
have to stay the same bits, which they can: the three ask the same function
for the same arguments in the same order today.

This was not implemented. "What could go wrong" of work package 4 of
`docs/plans/diversity.md` says that anything beyond the shared weights of
deliverable 2 is a performance review of its own and not a task of this
plan.

## What the code already does well

Three patterns worth copying.

The recurrence of `add_the_bins_of_the_var` is the largest measured gain in
this module and nobody set out to make it: it was written so that a draw
above about 560 called alleles would not lose its variants, and it is worth
22.9 times on a whole pass. The arithmetic that is right at the sizes of the
objectives turned out to be the arithmetic that has the fewest operations.

A chance is carried as a mantissa and a power of two,
`ChanceOfAPowerOfTwo`, so that a product of hundreds of factors below 1 can
pass through 1e-360, where an `f64` holds 0, and come back to a bin of the
order of 1. The cost is two comparisons per factor and it buys a whole
class of datasets that would otherwise give silent zeros.

An allele that a population did not call is skipped rather than multiplied:
`the_alleles_called` walks only the counts above 0, and
`add_the_standardized_private_alleles` passes over an allele that no
population of the row called. Each skipped allele is a product of `g`
factors not taken, and a population counted by reading a row as it lies has
no bound on the alleles of that row.
