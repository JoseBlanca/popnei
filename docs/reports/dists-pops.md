# Work report: the distances between populations

The plan `docs/plans/dists-pops.md` is done, on the branch
`plan/dists-pops`, in the worktree `.claude/worktrees/plan-dists-pops`,
23 September 2026. All eleven of its tasks are carried out, every
deliverable of its three work packages was checked by running it, and both
reviews are done and their findings fixed. Nothing is merged into `main`
and nothing is pushed: that is asked of the owner at the end of this
report, with four smaller decisions.

## What exists now that did not

A user calls `calc_pop_dists(variants, pops, jackknife_group=...)` in
Python or `calcPopDists` in TypeScript and gets, for every pair of
populations out of one pass over the variants, seven numbers with a
jackknife standard error each: Hudson's F_ST, f_2, the chord distance,
Nei's D_A, Jost's D, Nei's G_ST and the standardized G''_ST, with the
variants that counted for each pair, the f_2 of each pair within each
resampling group, the groups themselves and the counts of the pass.

How far the numbers are from the programs they were checked against, each
measured on the commit this report ends at:

| against | measure | how far |
|---|---|---|
| pyNei, live | Jost's D | 5.5e-15 relative, worst of three pairs, biallelic panel; 8.5e-15 multiallelic |
| adegenet 2.1.11 | the chord distance | the same double for all three pairs of the biallelic panel read in one block; 6.7e-16 absolute, 2.0e-15 relative, multiallelic |
| plink2 v2.0.0-a.7.7 | Hudson's F_ST | 4.9e-7 absolute, biallelic; 4.0e-8 multiallelic |
| ADMIXTOOLS 2.0.10 | f_2 and its standard error | within 1e-12 relative at three block sizes |
| mmod 1.3.3 | Jost's D, G_ST, G''_ST | 7.2e-5 to 1.9e-4 biallelic, 9.0e-5 to 4.7e-4 multiallelic, which is an estimator difference and not an error |

The tests: 61 in the core module where the branch began with 0, 39 in
Python and 236 in the TypeScript package, where neither file existed. The
whole suite is 582 cargo tests in the core crate, 353 Python tests and 236
TypeScript tests, all passing, with every check of the `coding` skill
green.

The spec gained fifteen paragraphs it did not have, each in a commit of
its own before the code that needed it, and `docs/architecture.md` gained
a row for the new module. Everything the spec left unsaid that the work
ran into is written down; the plan has no open point left.

## This report is written while the work goes

It has a section for each work package as it finished, in the order they
were done. What the owner has to decide is at the end.

## What was in place before the first task

Everything the plan's "What has to be in place" names, checked by running
it on the commit this branch starts from, `6704347`:

| what | command | what it gave |
|---|---|---|
| plink2 | `plink2 --version` | `PLINK v2.0.0-a.7.7 M1 (18 Sep 2026)` |
| R | `R --version` | `R version 4.6.1 (2026-06-24)` |
| adegenet, mmod, admixtools | `Rscript -e 'packageVersion(...)'` | 2.1.11, 1.3.3, 2.0.10 |
| the reference data | `ls tests/reference/pop_dists/` | `make_reference.py`, `micro.vcf.gz`, `micro_pops.txt` and the output of every program on both panels, 17 files |
| the core module | `cargo test -p popnei --lib pop_dists:: -- --list` | `0 tests, 0 benchmarks`, exit 0 |
| the Python tests | `ls tests/test_pop_dists.py` | no such file |
| the TypeScript tests | `ls js/popnei/test/pop_dists.test.ts` | no such file |

The three programs print the versions `make_reference.py` refuses any
other of, so the literals in the tests are of the data in the repository.
The last three rows are the baseline the plan's checks are written
against: a cargo selector alone passes on an empty crate, which is why
each check names how many tests have to run.

## Work package 1, as it goes

### Tasks 1.1 and 1.2, the counts of a variant and the resampling groups

Both went to one subagent, in one prompt, because both write
`crates/popnei/src/pop_dists.rs` and one tree has one writer per file.
Three commits came back: `f310edd`, an addition to the spec, `2fc9b0a`,
task 1.1, and `f732cf3`, task 1.2. 185 497 tokens.

`crates/popnei/src/pop_dists.rs` now holds the counts of one variant in
one population, the five values a pair makes of them, and the walk that
cuts the variants into the groups the standard errors are resampled over.
`cargo test -p popnei --lib pop_dists::` gives `16 passed; 0 failed`,
against `0 tests` when the branch started, and the whole workspace gives
`537 passed; 0 failed; 2 ignored`. The other checks of the `coding` skill
were run again by the orchestrator on `f732cf3`: `cargo fmt --all
--check` no output, `cargo clippy --workspace --all-targets -- -D
warnings` and `cargo wasm-check` both `Finished`.

Deliverable 1 is met: the tests assert the spec's counts table and its per
variant table for all four variants of the worked example, the literals
read from the spec, variant 3 where both populations are fixed and H_b is
0 and variant 4 where f_2 is -0.1 among them, and the corrected H_S and
H_T that no measure reads until work package 2. Deliverable 2 is met: the
test that cuts the biallelic panel reads the file itself and asserts the
12 groups of 100 variants over its two chromosomes, with the anchoring, a
group for each variant, no group at all, and a length of 0 as an error.

### What the spec did not say, and now does

A population that has called one allele at a variant has no within
population heterozygosity there: the correction is n_P / (n_P - 1) times
1 - sum over a of p_Pa^2, which is 0 over 0 at n_P = 1, and the variant
would have added a NaN to every sum of every pair that population is in.
The spec said what the six sums do when a population has little data and
not what they do with one called allele. The variant now counts for no
pair that population is in, which is the rule the spec already had for two
populations of one called genotype each, and it is written into the item
in `f310edd`, a commit of its own before the code.

It takes a haploid dataset and a `min_num_individuals` of 1 to reach,
since a genotype of two alleles or more that was called whole gives two
called alleles or more, so no value of either panel changes. The
orchestrator took it as a small choice made and written down rather than
an open point for the owner, as the `coding` skill allows: the spec had
defined nothing here, so nothing the owner decided was overturned, and the
alternative was a silent NaN.

### Task 1.3, the two measures and the standard error

Two commits: `0dccb3a`, the spec and the reference data, and `b56f402`,
the code. 188 529 tokens. `cargo test -p popnei --lib pop_dists::` gives
`27 passed; 0 failed`, and the workspace `548 passed; 0 failed; 2
ignored`. The other checks were run again by the orchestrator on
`b56f402` and all pass.

Deliverable 3 is met: the worked example's F_ST of 0.276712 and f_2 of
0.140278, the per variant F_ST of `var0000` of 0.332322 over a pass of
that one variant, plink2's three F_ST on each panel within 1e-6 absolute,
and ADMIXTOOLS 2's three f_2 and three standard errors within 1e-12
relative, are all asserted as literals.

### A test that could not have failed, and the run that fixes it

The 12 groups that a length of 100 000 base pairs cuts the biallelic panel
into all hold exactly 100 variants. The estimator the spec asks for, the
delete-m jackknife for unequal m, weights each group by h_j, the variants
of the pair over the variants of the group, and on 12 equal groups h_j is
12 for every one of them. An estimator written for groups that are all the
same size reproduces ADMIXTOOLS' three standard errors on that fixture to
the last bits of a double, so the check the plan named could not have told
the two apart.

So `make_reference.py` now runs ADMIXTOOLS 2.0.10 a second time, at
`blgsize = 250000`, which cuts the panel into 6 groups holding 250, 250
and 100 variants on each of its two chromosomes, and writes
`panel.f2.uneven.tsv` beside `panel.f2.tsv`. On those groups the spec's
formula is 9.8e-17 from the furthest of ADMIXTOOLS' three standard errors,
5.7e-14 of it, and the estimator that takes the groups as equal is 1.5e-4
from the first, which is 1.2e8 times the tolerance the tests compare
within. The 100 000 run came out byte for byte the same file, which is
what says the script and the versions it refuses are the ones the earlier
numbers were taken with.

This is a change to the plan, and it makes deliverable 3's check stronger
rather than weaker: the deliverable now has two fixtures where it had one.
The spec's "The standard errors" and its f_2 item carry the second run and
the seventeen digits the 1e-12 comparison needs, in `0dccb3a`, a commit
before the code. The f_2 of the two runs differ in their last two digits,
since ADMIXTOOLS takes f_2 as a weighted mean over its blocks.

### Task 1.4, the pass over a reader

Two commits: `a32de6e`, the spec, and `f89d1dd`, the code. 240 065
tokens. `cargo test -p popnei --lib pop_dists::` gives `37 passed; 0
failed` and the workspace `558 passed; 0 failed; 2 ignored`. The other
checks were run again by the orchestrator on `f89d1dd` and all pass.

Deliverable 4 is met, and each of its cases is a test of its own: the
panel read in blocks of 100 and of 10000 and in pools of 1 and 4 threads,
a pass over a source with no variant, fewer than two populations, fewer
than 20 groups with their number in the message, and an error of the
reader given on as it is.

The reduction that the plan called the likeliest silent failure of this
work package: each chunk of 64 rows sums into its own run of pairs per
group, and the chunks are added into the pass in block order, so the
threads never join a float total in an order rayon chooses. Pools of 1 and
4 threads give the same bits, not merely the same number within the
tolerance. The test can fail: cutting the chunks by
`rayon::current_num_threads` instead makes it fail at the 15th digit.

The sums of a pair are 48 bytes and not the 44 the spec's "How it runs"
claimed: five f64 and a count, with padding after a u32. The count is now
a u64, which costs nothing since the padding paid for it, and the spec's
three figures are corrected to 48 bytes, 72 KB for 3 populations and 500
groups and 29 MB for 50 populations, in `a32de6e`.

### Task 1.4b, added to the plan: a third reference run

`calc_pop_dists` refuses a pass that gives fewer than 20 resampling
groups, which the spec asks for so that nobody is handed a standard error
built from three groups. The two ADMIXTOOLS runs the plan had cut the
biallelic panel into 12 and 6 groups, so neither could be asserted through
the Python or the TypeScript package, which deliverables 5 and 6 ask for.

So `make_reference.py` runs ADMIXTOOLS a third time, at `blgsize = 55000`,
which cuts the panel into 22 groups, ten of 55 variants and one of 50 on
each of its two chromosomes: above the minimum, and still of two sizes, so
the estimator for unequal m is told apart at the Python level too. Its
standard errors are 0.0020502481330704485, 0.0016837006366670754 and
0.0019859616713111257, in `tests/reference/pop_dists/panel.f2.min20.tsv`
and in the spec, commit `9fad153`. That ADMIXTOOLS really cut 22 groups
was checked by asking it for the array, whose third dimension is 22 and
whose names are ten `l55` and one `l50` per chromosome.

The other two runs came out byte for byte the same files, which is what
says the script and the versions it refuses still give the numbers the
spec carries. The task was added by the orchestrator rather than put to
the owner because it keeps deliverables 5 and 6 as the plan wrote them;
the alternative was to drop the standard errors from both, which would
have made a deliverable's check weaker.

### Task 1.5, the Python side

Two commits: `aa821ec`, the spec, and `0532d0e`, the code. 248 657
tokens. `uv run pytest` gives `335 passed`, of which `tests/test_pop_dists.py`
is 23, where that file did not exist. The workspace gives `560 passed; 0
failed; 2 ignored` and the other checks pass.

Deliverable 5 is met. Its one check with no Python-visible consequence is
that a pass with `jackknife_group=None` asks the reader for no positions:
what a Python test can see is that there are no standard errors, no
`f2_groups` and no group, and the test asserts those three and points at
the cargo test that holds what the reader was asked for.

The seven names of the measures went into the core beside the ones of
`PerVarStat`, rather than into a table of the binding crate, because task
1.6 needs the same seven.

The spec gained two cases it had not said, in `aa821ec`, before the code:
`square_standard_errors` of a result that has no standard errors gives
`None`, the same answer as the field, rather than an empty frame or a
refusal; and the `Distances` of each measure names its pairs by the
populations, so that its square frames are indexed by them on both sides,
as the Kosman distances' are by the individuals.

### One thing for the owner, which nothing in the plan rests on

`repr` of a `Distances` of populations prints `<Distances of 3
individuals, 3 pairs>`. The class was written for the Kosman distances
between individuals and the spec gives it only the new `standard_errors`
field for this item, so the noun stayed. It is a string a user sees and
the plan does not depend on it, so the work goes on and it is put to the
owner at the end of the plan with its two options: the repr loses the
noun, or `Distances` learns what its names are of.

### Task 1.6, the TypeScript side

One commit, `c8e1285`. 283 373 tokens. `npm run build && npm test` in
`js/popnei` gives `tests 229, pass 229, fail 0`, of which
`test/pop_dists.test.ts` is 21, where that file did not exist. Every other
check passes on the same commit.

The cutting of the flat arrays of populations moved out of `stats.rs` into
one helper the two calculations share, rather than a second copy. A
`jackknifeGroup` that is 0, negative, fractional or above 2^53 is refused
in the binding crate, where the Python binding refuses 0, so the package
itself only tells the three kinds of value apart.

One difference between the two packages that no spec change was needed
for: in TypeScript each `Distances` of a result carries the counts of the
pass in `passStats`, which that class requires, and in Python
`Distances.pass_stats` is `None` there and the counts are on `PopDists`
alone. Both carry them on the result, which is where a user reads them.

## Work package 1: its six deliverables, each checked by the orchestrator

Run on `c8e1285`.

| deliverable | command | what it gave |
|---|---|---|
| 1, 2, 3 and 4, the core | `cargo test -p popnei --lib pop_dists::` | `39 passed; 0 failed` |
| 4, that the tests exist | `cargo test -p popnei --lib pop_dists:: -- --list` | `39 tests, 0 benchmarks`, where the plan asks for 20 or more and the branch started with 0 |
| 5, Python | `uv run pytest tests/test_pop_dists.py` | `23 passed`, where the file did not exist |
| 6, TypeScript | `npm test` in `js/popnei` | `tests 229, pass 229, fail 0` |

And the whole of the `coding` skill's checks on the same commit: `cargo
fmt --all --check` clean, `cargo clippy --workspace --all-targets -- -D
warnings` `Finished` with no warning, `cargo test --workspace` `560
passed; 0 failed; 2 ignored` and `35 passed`, `cargo wasm-check`
`Finished`, `ruff format --check` `27 files already formatted`, `ruff
check` `All checks passed!`, and `uv run maturin develop && uv run
pytest` `335 passed`.

All six deliverables are met.

## The review of work package 1

Seven reviewers, one for each category the `code-review` skill has, each
with a fresh context, over `6704347..c8e1285`. Two of them, `spec` and
`tests`, worked in worktrees of their own, as that skill asks, because
they change code to see what happens.

Two of them rebuilt the arithmetic from the spec on their own and compared
it with popnei: the `spec` reviewer over F_ST, f_2 and the f_2 standard
error of the three pairs of both panels at three group lengths and at a
`min_num_individuals` of 47, and the `numbers` reviewer over the same and
the jackknife's three formulas. Both agree with popnei to 1e-14 relative
or better. The measures themselves are right, and what the review found is
of another kind: one silent wrong result, two tests that could not fail,
and a set of messages and doc comments that say what is not so.

### What was found and what was decided

**A source whose positions go back gives wrong groups and says nothing.**
Four reviewers found it, from four sides. The walk that cuts the groups
decides with `pos.saturating_sub(filling.start) < length`, which is true
for every variant at or before the group's first, so such a variant joins
the group it is far from instead of starting one. The orchestrator ran it:
the biallelic panel with one variant moved out of order gives 122
resampling groups where the sorted file gives 240, a group whose recorded
end, 599 000, is before its start, 600 000, and a standard error for one
pair of 0.0020831 against 0.0017755, 17 in 100 higher. The f_2 themselves
do not move, which is why nothing else catches it.

It is reachable. `docs/specs/filters.md` says the linkage disequilibrium
filter "is the one part of popnei that needs the variants of each
chromosome to come together and in order of position. The rest of popnei
reads a source in any order", and `docs/specs/dists.md` said nothing about
order. So the resampling groups are a second part that needs it, and the
spec did not say so.

The fix is to refuse such a source, not to accept it. Sorting would mean
holding the whole dataset, against the streaming the objectives ask for,
and a standard error resampled over groups that are not the stretches the
user asked for is a number nobody should quote. That is also what the
linkage disequilibrium filter already does, so popnei has the pattern. It
is written into the spec first, in a commit of its own.

**The test of the threads could not fail.** The plan named this the
likeliest silent failure of the work package, and the test written for it
does not guard it. The `tests` reviewer replaced the ordered reduction
with rayon's own `par_chunks().reduce()` and all 39 cargo tests still
passed, while showing that the numbers do move: f_2 for one pair is
`3fa515b0b9944282` on one thread and `3fa515b0b994427b` on four. The
fixture cuts the panel into groups of 50 variants and a chunk is 64 rows,
so each group is built from at most two chunks and any order of two terms
gives the same bits. The report's earlier claim that the test fails when
the chunks are cut by the number of threads is true but is about another
property: that mutation moves the chunk boundaries, not the order they are
joined in.

**The corrected diversities are pinned only where nothing can tell them
apart.** All four variants of the worked example give the two populations
the same number of called genotypes, and every test is at a ploidy of 2.
The `tests` reviewer showed three wrong formulas that pass all 39 tests:
the harmonic mean of the called genotypes replaced by the arithmetic mean,
the pooled diversity weighted by the population sizes instead of equally,
and the ploidy ignored altogether. Work package 2 reads exactly these two
sums, so a defect in them is silent until the comparison with pyNei at the
end of it. The plan's deliverable 1 claimed "nothing of this task goes
unchecked until work package 2", and that was false for these three.

**Smaller findings that hold.** A pair with no variant has code and no
test. An empty `pops` is refused with a message telling the user to leave
out an argument that has no default. Python leaves `pass_stats` empty on
each measure's distances where TypeScript fills it, and the docstring of
that field says it is empty only for distances no pass gave, which is now
untrue. Four ways of giving `jackknife_group` wrongly are tested in
TypeScript and not in Python. Several doc comments state a bound or a
number that is not the code's: a count said to be limited by a `u32` that
is a `u64`, and a standard error quoted from a different run than the test
uses. The map of the modules in `docs/architecture.md` still has no row
for this module and gives its old one pyNei's function name.

**One finding not taken.** The `errors` reviewer asked for a check on a
`saturating_mul` that could make a chunk's sums be skipped. The `numbers`
reviewer had checked the same line and showed the product is already
bounded by a `checked_mul` in the function that grows the sums, so it
cannot saturate. The evidence settles it and the line stays.

**One the plan was wrong about, not the code.** Work package 1's "What it
gives" says a user calling `calc_pop_dists` without naming measures gets
F_ST and f_2. That call raises, because the default is all seven and five
of them are not written until work package 3. The plan's sentence is
corrected rather than the default changed twice, since nothing is merged
until the plan is done and the default is right at the end of it.

### The fixes, first round: the wrong result and the two blind tests

Eight commits, `266c97a` to `1a23d0f`, 273 503 tokens. `cargo test -p
popnei --lib pop_dists::` gives `45 passed; 0 failed` against 39 before,
and `uv run pytest tests/test_pop_dists.py` `24 passed` against 23.

**The source whose positions go back is refused.** The walk keeps the
chromosome and the position of the variant before it and the chromosomes
it has already read. Where the groups are a length, a position below the
one before it on the same chromosome, or a chromosome that comes back
after another, is an error. The first is the stricter of the two rules
that would do, and it is the one the linkage disequilibrium filter
already uses. `"variant"` and no groups at all take a source in any
order, since neither cuts by position, and the spec says so.
`docs/specs/filters.md` no longer claims to be the only part of popnei
that needs the order. Both cases are a `ValueError` in Python with the
file named. Four cargo tests and one pytest.

**The test of the threads can now fail.** It compares one thread against
four over a pass with no groups read in blocks of 10 000, which puts 19
chunks into one accumulator, on the bits. Under the reduction that joins
in rayon's own order the first comparison still passes and the second
fails, with F_ST of one pair 0.10496244498389444 in one run and
0.1049624449838943 in the other.

**The corrected diversities are now pinned where they can be told apart.**
The worked example gained a fifth variant, `0/0 0/0 0/1` against `1/1 ./.
1/2`, where the two populations have 3 called genotypes against 2 and 6
called alleles against 4, and a case at a ploidy of 4. Each of the three
wrong formulas the reviewer demonstrated now makes a test fail. The new
numbers were computed in exact rational arithmetic written from the
spec's formulas, which reproduces every literal the example already had,
and cross-checked against pyNei: its `_calc_pairwise_dest` gives the same
corrected H_S and H_T per variant and its `calc_jost_dest_pop_dists` the
same D.

**Every number of the worked example moved**, because the example is now
five variants: F_ST 0.351196 and f_2 0.203889 where they were 0.276712
and 0.140278, and, for the work packages still to come, D 0.347531, G_ST
0.238256, G''_ST 0.598618, D_A 0.258542 and the chord distance 0.508470.
The plan quoted the old ones in the checks of three of its deliverables
and now quotes these, so that the work packages after this one do not
chase a literal that no longer exists.

**One claim of the spec was dropped rather than tested.** The f_2 item
said the multiallelic panel is checked through plink2, as plink2's F_ST
times popnei's own sum of H_b. Since every measure divides the same
difference of sums, that product is the F_ST comparison again and checks
nothing further. The item now says plainly that the multiallelic f_2 is
checked against no program.

**What the owner may want to know.** Work package 2 divides by 1 minus the
mean corrected H_S, and on the five variant example that mean is 0.357143,
so the division stays far from zero there. A case near it is still
unwritten, and work package 2's deliverable 2 asks for one.

### The fixes, second round: the errors and the messages

Seven commits, `73f6b6a` to `887ec0d`, 242 986 tokens. `cargo test
--workspace` gives `568 passed; 0 failed; 2 ignored` and `uv run pytest`
`340 passed`, against 566 and 336 after the first round.

**A Ctrl-C during a calculation now raises a keyboard interrupt.** It did
not, and this is the one finding of the review that was a defect of code
older than this work package. Between releasing the interpreter and
building the numpy array, the Kosman distances asked numpy for its array
interface without first answering the signal, and a Ctrl-C there came out
as an exception that derives from `BaseException` and ends the user's
session. The helper that answers it existed, in the principal coordinate
analysis, and the new population distances had a copy of it written out by
hand. It is now one helper in one place, called by all three, and
`tests/test_interrupt.py`, whose three tests all interrupted the reading
of blocks and none a calculation, interrupts both calculations.

**An empty `pops` is refused with what the caller must do.** It reached
the shared refusal of the statistics, whose message ends "leave `pops` out
for one population of every individual" — advice a caller of this function
cannot follow, since `pops` has no default here. Both packages now refuse
it before the core with the rule this function has, two populations at
least, and the statistics keep their own message.

**Three places where a defect of popnei would have become a number.** A
variant added into sums of the wrong size returned without a word, and now
raises; a pair the core has no count for became a count of 0, which is
also the legitimate count of a pair with no variant, and now raises in
both binding crates. Neither is reachable today, which is why they were
findings rather than failures: the owner's rule is that a silent wrong
result is fixed whatever its odds.

**Two messages named the wrong cause.** A machine that cannot hold the
sums was told to cut the variants into longer groups even when no variant
had been read and the populations alone were the cause, and a ploidy
popnei cannot read was reported as the ploidy of a statistic of one
variant. Each now names what the user gave. The spec's error list gained
the three cases it did not have.

**One finding not taken, with the reason.** An overflow of the counter of
groups in `standard_error` becomes no standard error rather than an error.
The subagent answered that the counter is bounded by the groups, which the
memory already bounds, so it cannot overflow, and that raising there would
change the `Option<f64>` the spec fixes for that method. The orchestrator
weighed that as it weighs a reviewer and took it: the bound holds.

### The fixes, third round: the two packages, the doc comments and the tidying

Eleven commits, `4a1d8fc` to `b0dbdb1`, 238 297 tokens.

**The two packages now give the same thing.** Python left the counts of
the pass empty on each measure's distances where TypeScript filled them,
and the docstring of that field said it is empty only for distances no
pass gave. Python fills it, which makes the sentence true again.

**Three tests where there were none.** A pair with no variant, which the
spec describes and nothing exercised, reached at a `min_num_individuals`
of 50 where the largest count of called genotypes in the smallest
population is 48: `num_vars` comes out 0, 0 and 1200 and the third pair
still carries ADMIXTOOLS' number. The four ways of giving
`jackknife_group` wrongly that TypeScript tested and Python did not. And
the two branches of the Python that reads `measures`.

**The doc comments that were wrong.** A count said to be bounded by a
`u32` that is a `u64`; a method said to give no value where the pair's
measure has none, which it never consults; and a standard error quoted in
a test's comment as 0.0032 where the test's own group length gives
0.004197 — 0.0032 is the value at the other group length. The subagent
measured both before writing either.

**What `"variant"` costs, which nothing said.** The spec said the
accumulator does not grow with the variants, and under `"variant"` a group
is a variant, so it grows exactly with them. Measured on the biallelic
panel cut into 20 populations, which is 190 pairs and 1200 variants: the
sums are 10 944 000 bytes and `f2_groups` 1 824 000. Both the spec and the
docstring of `jackknife_group` in the two packages now say it, which is
where a user choosing `"variant"` will read it.

**Three pieces of tidying for the work packages to come.** The order of
the pairs was derived three times, once in the core and once in each
binding crate, so a change to the core's order would have left the values
and their standard errors mismatched inside one result without a word; the
bindings now take it from the core. Which measures have a value was known
in the two packages and not in the core, so work package 2 would have had
to change three places and dropping one refusal without the other would
hand a user NaN read as a distance; it is now one list in the core that
both packages read. And the three helpers that cut a block into chunks,
which this work had borrowed from the statistics module, are in the block
module, which is where section 9 of `docs/architecture.md` puts block
handling.

`docs/architecture.md` now has a row for this module, and its `dists` row
no longer names a pyNei function the spec says popnei deliberately does
not have.

## Work package 1 is done

Its six deliverables, run again after the three rounds of fixes, on
`b0dbdb1`:

| deliverable | command | what it gave |
|---|---|---|
| 1, 2, 3 and 4 | `cargo test -p popnei --lib pop_dists::` | `49 passed; 0 failed`, against 39 before the review and 0 when the branch started |
| 4, that they exist | `cargo test -p popnei --lib pop_dists:: -- --list` | `49 tests, 0 benchmarks`, where the plan asks for 20 or more |
| 5, Python | `uv run pytest tests/test_pop_dists.py` | `30 passed`, against 23 before the review |
| 6, TypeScript | `npm test` in `js/popnei` | `tests 230, pass 230, fail 0` |

And every check of the `coding` skill on the same commit: `cargo fmt
--all --check` clean, `cargo clippy --workspace --all-targets -- -D
warnings` `Finished` with no warning, `cargo test --workspace` `570
passed; 0 failed; 2 ignored` and `35 passed`, `cargo wasm-check`
`Finished`, `ruff format --check` `27 files already formatted`, `ruff
check` `All checks passed!`, and `uv run maturin develop && uv run
pytest` `344 passed`.

The review added 10 cargo tests and 7 Python tests to this work package
and changed the spec in five more places.

## Work package 2, as it goes

### Task 2.1, the three measures out of the corrected sums

One commit, `9f0a6e9`. 155 644 tokens. `cargo test -p popnei --lib
pop_dists::` gives `53 passed; 0 failed` against 49, and the workspace
`574 passed; 0 failed; 2 ignored`. The work package added no count and no
sum, as the plan said it would not: the three measures are ratios of the
means of the two corrected diversities work package 1 already summed.

Deliverables 1 and 2 are met. How far the three came from mmod, which is
the number to watch because mmod computes another estimator of the same
quantity: on the biallelic panel the furthest of the three pairs is 7.3e-5
for Jost's D, 7.2e-5 for G_ST and 1.9e-4 for G''_ST, and on the
multiallelic panel 3.5e-4, 9.0e-5 and 4.7e-4. All six are inside the 5e-4
the plan asks for, and the tolerance was left where it is. The exact
comparison, against pyNei, is task 2.2's.

### One point for the owner, which nothing rests on

The standardized G''_ST divides by 1 minus the mean corrected H_S, the
same divisor as Jost's D. The spec gives Jost's D no value where that mean
is exactly 1, and gives G''_ST no such exception, in the words "for Dest
alone". So where Jost's D has no value, G''_ST comes out infinite. The
task followed the spec and pinned it in a test.

The word "alone" is explicit, so the orchestrator did not overturn it: it
is either the owner's decision or an oversight of whoever wrote the line,
and the two cannot be told apart from the page. The case needs the mean
corrected H_S to be exactly 1.0 and the task had to construct a fixture to
reach it, so nothing of the plan rests on the answer and the work goes on.
The two options are at the end of this report.

### Task 2.2, the three measures through both packages

One commit, `a54033a`. 153 655 tokens. `uv run pytest
tests/test_pop_dists.py` gives `33 passed` against 30, `npm test` `tests
233, pass 233, fail 0` against 230, and the whole Python suite `347
passed`.

Neither binding crate nor either package needed new code. The three
measures arrived the moment task 2.1 added them to the one list in the
core of which measures have a value, because the review of work package 1
had already made both bindings take that list, the order of the pairs and
the arrays from the core instead of deriving them. That finding was
recorded as a hazard of a future defect, and this is what it bought: the
work the plan gave two layers turned out to be tests and doc comments.

**The exact comparison, which is the strongest check popnei has for this
measure.** popnei's Jost's D against a live pyNei, worst of the three
pairs: 5.5e-15 relative on the biallelic panel at a `min_num_individuals`
of 20, 5.1e-15 at 47 where the pairs part, and 8.5e-15 on the
multiallelic panel. The bound asserted is 1e-12. The counts at 47 are
688, 688 and 1200, which is what pins `num_vars` being per pair.

Three things were done beyond the deliverables, all of them additive.
pyNei is run on the multiallelic panel as well as the biallelic one, since
the spec prints its numbers for both and the TypeScript suite asserts them
as literals without a pyNei of its own. Both suites also assert popnei's
own ten digit G_ST and G''_ST from the spec, because 5e-4 against mmod
would let a change of 1e-5 pass unnoticed.

And the docstring of the standardized measure, in both packages, described
it as G_ST divided by the largest value it could reach. That is Hedrick's
G'_ST, a different measure, which the spec says popnei does not give. Both
now describe Meirmans and Hedrick's G''_ST, with the two numbers that tell
them apart on the biallelic panel, 0.1620 against 0.1155 for one pair, and
how a user who wants G'_ST gets it from `gst`.

**The order of the work could not be followed as written**, and the
subagent said so rather than pretending otherwise: the code paths already
existed, so its new tests passed the first time they ran. It checked
instead that each can fail, by moving one literal at a time and seeing the
named test fail, and reverted both.

## Work package 2 is done

| deliverable | command | what it gave |
|---|---|---|
| 1 and 2, the core | `cargo test -p popnei --lib pop_dists::` | `53 passed; 0 failed`, against 49 |
| 3, Python | `uv run pytest tests/test_pop_dists.py` | `33 passed`, against 30 |
| 4, TypeScript | `npm test` in `js/popnei` | `tests 233, pass 233, fail 0` |

It is reviewed together with work package 3, which the `code-review` skill
allows when two work packages are one piece of code: neither adds a count
or a sum, both add ratios of sums work package 1 already accumulates, and
both touch the same function.

## Work package 3

### Task 3.1, the chord distance and Nei's D_A

Three commits, `86d38a8` and `86bbc8c` for the spec and `3a9e540` for the
code. 164 542 tokens. `cargo test -p popnei --lib pop_dists::` gives `56
passed; 0 failed` against 53, and `cargo test --workspace` `577 passed; 0 failed; 2
ignored` in the core crate beside `35 passed` in the linear algebra one. No count and no sum were added: D_A is 1 minus the mean of the
sum of the square roots the pass already builds, and the chord distance is
its square root.

Deliverable 1 is met, and the agreement with adegenet is closer than the
1e-12 relative asked for. On the biallelic panel the chord of all three
pairs is the same double the reference file prints, 0 apart. On the
multiallelic panel the furthest is 6.7e-16 absolute, 2.0e-15 relative.

**A case of the last bits that would have given a NaN.** The sum of the
square roots can exceed the variants that counted by the last bits of a
double: two populations with the same frequencies at every variant give
0.2 + 0.4 + 0.3 + 0.1 = 1 + 2.2e-16, which makes D_A a small negative
number and the chord distance, its square root, a NaN. D_A is 0 there
instead, which the spec now says, in a commit before the code, and a test
of two populations alike at every variant pins it.

**A number the spec claimed before the code existed.** The chord item said
popnei was 1.25e-15 from adegenet on the biallelic panel and 3.3e-16 on
the multiallelic. Both were written before there was code to measure. The
item now carries what was measured, 0 and 6.7e-16.

### One more point for the owner, which nothing rests on

With all seven measures having a value, the machinery that refuses the
name of a measure that has none refuses nothing: the list in the core, a
function in each binding crate, and two messages that are now unreachable
and untested. Task 3.1 kept it rather than drop public API from three
crates on its own judgement. The options are at the end of this report.

### Task 3.2, the two measures through both packages, and the check of the whole plan

One commit, `bbada67`. 124 630 tokens. Neither layer needed new code, for
the same reason task 2.2 did not: both packages take from the core the
list of which measures have a value, so the chord distance and D_A
appeared in them the moment task 3.1 landed. What the commit adds is the
tests and three doc comments that still said the two measures were not
calculated yet.

**The check of the whole plan**, which no single work package covers,
passes: `measures` left out gives all seven on both panels, asking for one
measure gives the same double as asking for all seven, and at a
`min_num_individuals` of 47 on the biallelic panel every one of the seven
numbers of the two pairs of p0 moves while every one of p1-p2 is
unchanged, which is what shows the seven are over the same variants and
the pair has one count of them. It was seen to fail when the single
measure call was made at another threshold.

D_A is asserted against the square of adegenet's literal and not against
popnei's own chord distance, and separately the round trip between the two
is checked, which is 2.2e-16 at the furthest. All four of the new
assertions were seen to fail: multiplying the chord by the square root of
two in the core makes them report 0.2549 where 0.18027 is expected, which
is the mistake the plan warned about.

## Work package 3 is done

| deliverable | command | what it gave |
|---|---|---|
| 1, the core | `cargo test -p popnei --lib pop_dists::` | `56 passed; 0 failed`, against 53 |
| 2, Python | `uv run pytest tests/test_pop_dists.py` | `36 passed`, against 33 |
| 2, TypeScript | `npm test` in `js/popnei` | `tests 235, pass 235, fail 0` |

And every check of the `coding` skill on `bbada67`: `cargo fmt --all
--check` clean, clippy `Finished` with no warning, `cargo test
--workspace` `577 passed; 0 failed; 2 ignored` and `35 passed`, `cargo
wasm-check` `Finished`, ruff `27 files already formatted` and `All checks
passed!`, `uv run pytest` `350 passed`.

## The review of work packages 2 and 3

They are reviewed together, which the `code-review` skill allows when two
work packages are one piece of code. Five reviewers: `spec`, `tests`,
`numbers`, `errors`, and one taking `api` with the binding layer, which is
small here because neither binding crate needed new code.

`architecture` was not sent, and this is the reason: that category applies
when a change touches the readers, the blocks, the threads, the linear
algebra, a `cfg` or a dependency, and these two work packages touch none
of them. They add five ratios inside one function, over sums work package
1 already accumulates. The `architecture` reviewer of work package 1 went
over the pass, the accumulator, the two builds and the memory, which is
where those questions live, and nothing here changes them.

### The review of work packages 2 and 3: what it found

Six reviewers: `spec`, `tests`, `numbers`, `errors` and one taking `api`
with the binding layer. Three of them rebuilt the five measures
independently from the spec, one with its own VCF reader, and all three
agree with popnei to 1.8e-14 relative or better on both panels at two
thresholds. The `spec` reviewer also re-ran adegenet and mmod itself
under R and reproduced all eighteen reference doubles digit for digit. So
the measures are right wherever their inputs are in range, and everything
below is about the edges and about what a reader is told.

**A 0 over 0 that destroyed a standard error that existed.** Two
reviewers found it. Two populations fixed for the same allele at every
variant that counted make the two corrected diversities 0, so G_ST and
G''_ST divided 0 by 0 and returned a NaN as though it were a value, with
a count of variants above zero. One such NaN, in one group left out,
poisons the whole jackknife. The orchestrator ran it on 25 variants of 4
individuals, 24 of them with every genotype `0/0`: F_ST, G_ST and G''_ST
came out with a standard error of NaN while Jost's D, f_2, the chord
distance and D_A all built a finite one from the same data.

It reached F_ST too, whose 0 over 0 work package 1 had documented as
deliberate on the argument that a caller sees the same thing whether the
core gives no value or a NaN. That argument was false, and this is the
case that shows it.

**What the standard error does when a group has no value, which the spec
did not say.** The fix had to settle it. A measure now gets no standard
error at all when a group that holds variants of the pair leaves it
without a value, rather than a standard error built from the groups that
do have one. The reason is in the estimator: the weights of the g groups
add to one, so dropping a group for want of a value pulls the jackknife
estimate off by that group's weight and the variance is then taken around
a centre no group put there. On the 25-variant fixture that route gives
0.0082 for F_ST where the honest answer is that the error is not
computable. It costs nothing on real data, since a left-out value is
undefined only when every other group is degenerate. It is in the spec,
in a commit before the code.

**Three measures that were a ratio of rounding noise at a ploidy of 1.**
At ploidy 1 the two corrected diversities are identically zero by their
definitions, so the sums hold only the residue of adding frequencies that
do not come to exactly 1, and Jost's D, G_ST and G''_ST divided residue by
residue. On popnei's own haploid fixture the orchestrator read G_ST
0.2503 and G''_ST 0.4004, where Jost's D showed what was really there,
1.4e-17. Moving `min_num_individuals` from 3 to 4 moved G_ST by 6% and the
chord distance by 0.2%, which is the signature of a ratio of noise; a
reviewer's own implementation of the same formulas gave G_ST 1.0 on the
same data, a different arbitrary number. The three now have no value at a
ploidy of 1, which is in the spec. F_ST, f_2, the chord distance and D_A
are meaningful there and are unchanged, which is why the pass is not
refused.

**A third blind test of the family work package 1's review found.** The
tests that check the measures do not change with the size of the blocks or
the number of threads looped over F_ST and f_2 alone. Those two read two
of the six sums, so the three sums that the other five measures are built
from were never compared across block sizes or thread counts at all. Two
reviewers made those three sums depend on how many chunks a block holds
and both tests still passed. Widened to all seven, the same mutation fails
with the chord distance of one pair 1.4e-9 apart between blocks of 100 and
of 10000.

**A fixture that was not a state a pass can reach.** The test of the case
where the mean corrected within-population diversity is exactly 1 was
built from hand-written sums, and its doc said no genotypes reach that
case. Genotypes do: four individuals in two populations of two, one
variant of four distinct homozygotes. And the reachable case behaves
differently from the hand-written one, since there both corrected
diversities are 1 and G''_ST is a 0 over 0 rather than an infinity. There
is now a test over real genotypes beside the hand-written one, which is
kept for the other case.

**Smaller things that held.** No test asserted that a standard error of
the five measures these two work packages added is a number at all. The
comparison with pyNei built its populations through a helper that sorts
the names, so popnei keeping the order of the `pops` dict, which is one of
the four differences from pyNei the spec names, was never exercised with
an unsorted dict. One TypeScript assertion compared the length of two
constants with each other and could not fail. The clamp that keeps D_A at
zero or above would have turned a NaN into a zero, since that is what the
maximum of a NaN and zero gives, which was a hazard rather than a defect
since no NaN can reach it today.

## What is asked of the owner

### The merge

The branch `plan/dists-pops` is not merged and not pushed. It is based on
`spec/dists-pops`, not on `main`, because `docs/specs/dists.md` is on that
branch and was not merged when this plan started. So the order to merge is
two orders: `spec/dists-pops` into `main` and then this branch, or this
branch into `spec/dists-pops` and the pair into `main`. Nothing else on
either branch conflicts with `main` as it stood when this plan began, but
`main` has moved since and the merge will say.

The speed was deliberately left out of this plan, and the performance
review measures it after the merge, which the owner decided on 23
September 2026. Nothing here has been timed and no speed is claimed.

### Four decisions, none of which the plan rests on

**1. A `Distances` of populations says "individuals".** `repr` of any of
the seven gives `<Distances of 3 individuals, 3 pairs>`. The class was
written for the distances between individuals and the spec gives it only
the new `standard_errors` field for this item, so the noun stayed. Two
subagents noticed it without being asked.

- Take the noun out of the repr: `<Distances of 3 names, 3 pairs>` or
  `<Distances of 3, 3 pairs>`. One line, no new field, and the individual
  case reads slightly worse.
- Give `Distances` what its names are of, a field set by whoever builds
  it. One field on a public frozen dataclass and its TypeScript twin, and
  both repr lines then read right.

Recommended: the second. The class now serves two things and a reader of
either should not be told the other.

**2. G''_ST where Jost's D has no value.** Both divide by 1 minus the mean
corrected within-population diversity. The spec gives the exception to
Jost's D "for Dest alone", so where that mean is exactly 1 and the two
corrected diversities differ, G''_ST comes out infinite. Where they are
equal it is a 0 over 0, which this plan's review made give no value like
the rest. Both are reachable from genotypes, the second by four
individuals in two populations with one variant of four distinct
homozygotes.

- Leave it. "Alone" is explicit and may be deliberate.
- Give G''_ST the same exception. One arm of one match, and the word
  "alone" comes out of the spec.

Recommended: the second, unless "alone" was a decision. An infinity is not
a missing value, and popnei's rule everywhere else is that a value which
cannot be computed is absent rather than a number a user might plot.

**3. Jost's D and G''_ST can leave the range 0 to 1.** When half-called
genotypes make a population's called alleles more than twice its called
genotypes, or at a ploidy above 2, the mean corrected within-population
diversity can pass 1, and then D and G''_ST cross a pole: a reviewer
measured D at -18.5, then -7.5e14, then 19.0 as one variant was added at a
time. pyNei does the same, and a reviewer confirmed popnei reproduces
pyNei's numbers there, so popnei is faithful and this is inherited, not
introduced. The spec says D runs 0 to 1, and the guard in the code tests
that divisor for exactly zero, which catches one point of something that
crosses zero.

- Leave it and say in the spec that the bound holds for whole-called
  diploids only, which this plan has not yet done.
- Give D and G''_ST no value where that mean is at or above 1, which
  departs from pyNei in a case pyNei gets wrong.

Recommended: the second, with the first done either way. The objectives
put being right first, and a distance of -7.5e14 is not a number to hand
anybody; but it is a departure from pyNei on a value a user sees, which is
why it is here and not decided.

**4. The machinery that refuses a measure with no value.** All seven
measures now have one, so the list in the core and the refusal in each
package refuse nothing. The two reviewers who were asked disagreed: one
would remove it, since the core already matches on the measure with no
catch-all arm and a measure without a formula is a compile error, so the
only case the separate list can still catch is a measure popnei computes
that someone forgot to list, which it would then refuse wrongly; the other
would keep the two lines and reword the messages. The orchestrator took
the second, because rewording removed the real defect, which was messages
that described the state of this plan, and removing public API from three
crates on a reviewer's opinion is not a thing to do without the owner.

- Leave it as it is now.
- Remove it, and take `THAT_HAVE_A_VALUE`, `has_a_value` and
  `names_that_have_a_value` out of the core and the two binding crates.

Recommended: leave it until f_3 and f_4 are written, which is a spec of
their own and the next thing that would add a measure. Then whoever writes
them will know whether it earns its place.

### Two things for an issue, not for this branch

Neither is this plan's code and both are small.

- Population names that look like numbers come out in a different order
  in TypeScript than in Python, because JavaScript puts integer-like
  object keys first. Nothing is mislabelled, since the result carries the
  population order, but the pairs come out in a different order in the two
  packages. It affects `docs/specs/stats.md`'s functions too, and the fix
  is an API decision: take the populations as an array of pairs, or refuse
  a name JavaScript would reorder.
- `cargo wasm-check` builds the core and the linear algebra crates only,
  so the pyo3 crate, which the `coding` skill says builds twice, is never
  checked for emscripten by the routine commands.
