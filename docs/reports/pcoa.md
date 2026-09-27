# Report: the principal coordinates of distances and Lingoes' correction

27 September 2026. The work report of `docs/plans/pcoa.md`, carried out on
the branch `plan/pcoa` in `.claude/worktrees/pcoa`, from the branch
`spec/pcoa` at aa78e7e, where the spec item "The principal coordinates of
distances" of `docs/specs/pca.md` and the plan were written. Under way.

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
