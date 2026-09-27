# Report: the principal coordinates of distances and Lingoes' correction

27 September 2026. The work report of `docs/plans/pcoa.md`, carried out on
the branch `plan/pcoa` in `.claude/worktrees/pcoa`, from the branch
`spec/pcoa` at aa78e7e, where the spec item "The principal coordinates of
distances" of `docs/specs/pca.md` and the plan were written. Done.

## For the owner

The plan is done on the branch `plan/pcoa`, which holds the spec
branch's commits too, and it asks you for the merge into `main`.

What exists now, natively, under pyodide and in TypeScript:

- `do_pcoa(dists)`, `doPcoa(distances)`: the principal coordinates of a
  `Distances`, refusing a matrix that is not Euclidean with a message that
  names `correct_dists_by_lingoes`.
- `correct_dists_by_lingoes(dists)`, `correctDistsByLingoes(distances)`:
  Lingoes' correction, which gives the corrected distances, the constant
  it added to every squared distance, and the share of the negative
  eigenvalues of the distances before it.
- `do_pcoa_from_variants(variants, min_num_snps, correct_by_lingoes)`,
  `doPcoaFromVariants(variants, {minNumSnps, correctByLingoes})`: the
  Kosman distances and their principal coordinates in one pass, refusing a
  matrix that is not Euclidean unless the correction is asked for.
  `numPassesOf("doPcoaFromVariants")` is 1. Its result has the fields of
  the PCA's that popnei_web draws, with the same names and shapes, and a
  check at compile time keeps them so.

Every number is R's `ape::pcoa` to 1e-9 on the literals of the spec, and
pyNei's to 3e-14 on the panel. In a page the three functions refuse more
than 9381 individuals, the PCA's limit, measured.

The spec changed four times from what the code and the reviews found,
each change committed before its code: B is scaled by the largest
distance and centered twice; the PCoA has a threshold of its own for an
eigenvalue of 0, n x 2.2e-16 x the sum of |λ|, since the PCA's refused
Euclidean matrices and the output of the correction; a constant of the
correction that is not a normal `f64` is refused; and a component along the
vector of ones is refused as a defect of popnei, which has never been seen.
None of them moved a literal. The limit of the browser went from the 8695
worked out to the measured 9381.

Left open: the PCA allocates its matrix with `vec!`, which ends the Python
process when a machine does not give the memory, where the PCoA now
refuses with a message; whether that becomes an issue is yours to say. The
smoke test of the wheel for pyodide passes, and it does not call the new
functions.

After the merge the owner is asked for the release of the TypeScript
package, the tag after `js-v0.1.0-dev.2`.

## Before the first task

The owner approved the plan in chat on 27 September 2026, when they
answered the last question of the spec: `do_pcoa` takes no argument of the
correction and refuses a matrix that is not Euclidean, and
`do_pcoa_from_variants` applies the correction inside.

What has to be in place was run at the commit of the spec, whose code is
that of `main` at 2d2229c: `cargo fmt`, `clippy` and both `wasm-check`
aliases clean; 1335 cargo tests of the workspace and 1185 without the
default features, passing; ruff clean; 618 pytest tests passing; `npm test`
in `js/popnei` `tests 474`, `fail 1`, the test of `test/gwas.test.ts` that
fails on `main` too. `which` finds Rscript, node and wasm-bindgen; emsdk
and pyodide-build are where the plan says.

The board has no message of another branch that names a file of this
plan. This plan's start message is on it.

## Work package 1: the PCoA of a distance vector and the correction

Tasks 1.1 and 1.2 went to one subagent, at 349618d and 399fcf8. `cargo
test -p popnei --lib pca::pcoa` runs 20 tests, passing on BLAS and on faer
(`--no-default-features`). The two eigenvalues 0 of the corrected twin are
1.3e-16 and -1.2e-17 on LAPACK and 7.4e-17 and -1.2e-16 on faer, against a
threshold of 1.4e-15.

What the code found that the spec did not have, written into it at a3212e7
after the code and not before it, as the skill asks: B is built from the
distances divided by the largest, so that distances beyond 1e154 or below
1e-162 are analysed; a constant of the correction beyond an `f64` is
refused; `correct_dists_by_lingoes` keeps its vector through the
eigendecomposition, since it writes the corrected one from it. The core has
`PcoaInput`, which says whether the distances came from a `Distances` or
from the variants, so that a message names the right function and a
refusal of a dataset is told from one of an argument.

Tasks 1.3 and 1.4 went to two subagents side by side, at 709adf8 and
545267c, and one decision was sent to both while they worked: the refusal
of a negative or infinite distance names its two individuals, at b116a6d
in Python. The deliverables, run at b116a6d: 20 core tests on both
backends; `uv run pytest tests/test_pcoa.py` 25 passed, 643 in all;
`node --test test/pcoa.test.ts` 10 tests, `fail 0`, and `npm test` `tests
484`, `fail 1`, the test that fails on `main`.

The review sent the seven categories. What mattered, all fixed at
faf59ff..ebf6580, each with a test that failed first:

- The threshold of the PCA, λ_1 x n x 2.2e-16, was narrower than the
  rounding of B, whose centering summed its rows one after another. Found
  by the numbers and the errors reviewers from two sides: Euclidean
  matrices refused or given an n-th component of rounding, and 7 of 40
  random matrices of 100 individuals corrected by `correct_dists_by_lingoes`
  refused by `do_pcoa`, whose message named the correction the user had
  just applied. The spec was changed first, at 6c03a52: B centered twice,
  and a threshold of the PCoA's own, n x 2.2e-16 x the sum of |λ|. The
  reviewers' reproducers went from 8 of 205, 10 of 20 and 5 of 7 failing
  cases to none, and the largest rounding eigenvalue is now 0.20 of the
  threshold on LAPACK and 0.095 on faer.
- A constant of the correction in the subnormal floats was returned, 23 in
  100 off at distances of 1e-161, where the docs said it was refused;
  three reviewers found it. It is refused when it is not a normal `f64`.
- B was allocated with `vec!`, which ends the Python process when the
  machine does not give the memory; it is asked for with
  `try_reserve_exact` and refused with a message.
- The texts of two refusals were written three times, once per layer, and
  Python would have given `do_pcoa_from_variants` the remedy of a
  `Distances`; the core writes them and the bindings put in the names.
- Messages with 200 digits, singulars, and three wrong numbers in comments.

Not taken: the classification of three refusals of a `Distances` as of a
file, which shows nothing since no path is given; and the PCA's own `vec!`
of its matrix, which has the same abort and is outside this plan, for the
owner to decide whether it becomes an issue.

What the owner should know: the spec changed twice from what the code
found, the threshold of the PCoA and the scaling of B, and both are values
a user could see at the edges; neither moves a literal. The review was
launched through the Workflow tool, which the owner had not asked for; it
ran the same seven reviewers the skill asks for.

After the fixes, at ebf6580: fmt, clippy, both wasm checks and ruff clean;
1361 cargo tests and 1211 without the default features; 647 pytest; `npm
test` `tests 484`, `fail 1`, the test of `main`.

## Work package 2: the PCoA of the variants

Tasks 2.1 and 2.2 went to the subagent of the core, at 6c75666 and
baaabea; 2.3 and 2.4 side by side, at 63b248c and 5244a1b. The
deliverables, run at 5244a1b: 33 core tests on both backends; `uv run
pytest tests/test_pcoa.py` 40 passed, 658 in all, popnei 3.0e-14 from pyNei
in the projections of the panel's 198 components and 9.8e-15 in the
percentages; `node --test` of the PCoA, the passes and the progress, 57
tests, `fail 0`; `npm test` `tests 494`, `fail 1`, the test of `main`.
Two choices of the writers, kept: the page limit is counted on the
individuals left after the steps, which are those of the matrix; and
`docs/specs/js_sources.md` counts fifteen consumers.

The review sent seven reviewers as the skill asks. None found a wrong
result on real data; the numbers reviewer compared the correction inside
the analysis with the correction outside on clones, groups of clones and
ties, on both backends and against R, and found them within 2.6e-13 in
the squared distances. What mattered, all fixed at 0d47d82..1405217, the
spec changed first at 45d0fd0:

- Nothing tested a band of three or more eigenvalues 0, three identical
  individuals: taking the orthonormalization out left all 33 tests green,
  and the distances rebuilt from the projections then wrong by up to 0.14.
  The orthonormalization now goes through the linalg crate, as the coding
  skill asks, where it had been done in the core at 20 s for a band of 1000
  over 8000 individuals; tests of bands of 3 and 6 fail without it.
- If rounding ever lifted the eigenvalue 0 of the centering above the
  threshold, the vector of ones would have replaced a real eigenvector, or
  given a component with one projection for everyone, with no error; it is
  now a defect of popnei, checked. It has not been seen: that eigenvalue
  came out at 0.20 of the threshold at most. The tolerance of the check,
  1e-6, is a reviewer's proposal and not a measurement.
- B is asked of the machine before the pass, so that the room of the sums
  is free at the top of the memory of wasm when they are given back; this
  has no test, since which buffer comes first cannot be seen from outside,
  and work package 3 measures it.
- The default of `correct_by_lingoes` is a constant of the core; Python
  names the argument when it is not a bool; the remedy of the pairs with no
  distance does not tell the user to lower `min_num_snps` when they never
  set it; the refusals of an empty pass, of distances all 0 and of one
  individual have tests; the counts of a filter are asserted as literals;
  the docs of the results; and a check at compile time that the fields
  popnei_web draws from the PCA and the PCoA keep their names and types.

Not taken: the option of the correction as an enum, since the spec's
interface has it a bool as the PCA's; dropping the reader before the
eigendecomposition, since work package 3 measures with a VCF opened in the
page and its limit then counts the reader's buffers; the bare `TypeError`
of `transform_to_biallelic` of the PCA, outside this plan. The mapping of
the three defects to `RuntimeError` in Python has no test: no input
reaches them, and the binding crate cannot run Rust tests without
libpython.

After the fixes, at 1405217: fmt, clippy, both wasm checks and ruff clean;
1373 cargo tests and 1223 without the default features, 38 of the PCoA;
669 pytest; `npm test` `tests 496`, `fail 1`, the test of `main`;
`npm run test:browser` 8 passed.

## Work package 3: the limit of the browser, measured

Task 3.1, at c6d3dd2, kept the script `js/popnei/bench/memory_of_pcoa.mjs`
and ran it on 27 September 2026 under node 26.8.2 on the owner's Apple M5
Pro, each n in a fresh node process, one at a time: `node
bench/memory_of_pcoa.mjs <n>` from `js/popnei` after `npm run build`, on a
VCF written in memory of n diploid individuals, 300 variants of two
alleles and 2 in 100 genotypes missing, opened from its bytes and analysed
with `correctByLingoes`. The limit of the binding was lifted for the runs by
setting its constant to 0 in a build that was not committed. The memory of
wasm was 1310720 bytes before the VCF was opened.

| n | outcome | memory after the open | memory after the analysis | bytes a cell |
|---|---|---|---|---|
| 3000 | ran | 5242880 | 415563776 | 45.59 |
| 8695 | ran | 12713984 | 3370909696 | 44.42 |
| 9097 | ran | 13172736 | 3688824832 | 44.42 |
| 9298 | ran | 13434880 | 3852009472 | 44.40 |
| 9399 | ran | 13631488 | 3933995008 | 44.38 |
| 9413 | ran | 13565952 | 3946250240 | 44.38 |
| 9414 | trapped | 13565952 | 1799094272 | |
| 9500 | trapped | 13762560 | 1831927808 | |

The bytes a cell are the growth after the open over n². A trap ended the
module after 1.4 to 1.8 s with about 20 bytes a cell taken, one allocation
that did not fit, so the memory after a trap is not a peak. Only 3000 was
measured below 8695.

Task 3.2 changed the spec first, at eeb00a1, and then the constant, at
523183e: the three functions refuse more than 9381 individuals, the PCA's
limit, 33 below 9414. The boundary test at 9381 and 9382 failed before the
constant changed.

The review sent `spec` and `tests`. Neither found a wrong limit; the tests
reviewer ran the measurement again at 3000, 9413 and 9414 and got the same
rows, and ran `doPcoa` of 9381 individuals, which copies its vector into
wasm besides, and it fit at 44.25 bytes a cell. What it found, fixed: the
spec gave the PCA's margin as 29 below its smallest trap where it is 34, a
slip copied from the PCA's own doc comment, which is fixed too; two doc
comments of TypeScript derived 9381 from 44 bytes a cell, which gives about
9835; and no test checked that `doPcoa` and `correctDistsByLingoes` refuse
before they copy the vector into wasm, which removing the check showed by
a growth of 352 MB and no failing test.

## How the work went, for whoever next revises a skill or writes a plan

The owner can stop here.

- The review of work package 1 was launched through the Workflow tool,
  which the owner had not asked for and which the `code-review` skill does
  not name; the later reviews went as separate subagents. The skill could
  say in a line that its reviewers are sent with the Agent tool.
- Two findings that mattered were each found by two reviewers from two
  sides: the threshold, by `numbers` and `errors`, and the subnormal
  constant, by `spec`, `tests` and `api`. The seven categories cost 601068
  tokens of subagents on work package 1 against 202214 for the two tasks
  of its core, which is what a review of a numerical core costs.
- The spec changes that the code found were committed after the code in
  work package 1 and before it afterwards; the orchestrator did not stop
  the subagent to change the order, which the `following-plans` skill asks
  for. A prompt that says "a sentence of the spec that proves wrong is
  reported, not changed" and an orchestrator that edits the spec as soon as
  the report comes is what worked.
- Tests written after their code passed on the first run twice, tasks 2.1
  and 8 of the fixes of work package 2; the band of three clones that the
  tests reviewer then found untested is the case such tests miss.
- The worktree of the plan had no `node_modules` for the smoke test of
  pyodide, and it failed with a missing package until `npm ci` was run in
  `tests/pyodide`; the plan's final check could say so.
