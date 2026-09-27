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
