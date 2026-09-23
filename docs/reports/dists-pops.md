# Work report: the distances between populations

The plan `docs/plans/dists-pops.md` is under way on the branch
`plan/dists-pops`, in the worktree `.claude/worktrees/plan-dists-pops`,
since 23 September 2026. It builds seven distances between every pair of
populations out of one pass over the variants, each with a jackknife
standard error, through the core crate, both binding crates and both
packages. Nothing is merged into `main` and nothing is pushed.

This report is written while the work goes. It has a section for each
work package as it finishes, and the top of it will say, when the plan is
done, what exists that did not and what is asked of the owner.

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
