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
