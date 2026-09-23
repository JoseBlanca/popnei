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
