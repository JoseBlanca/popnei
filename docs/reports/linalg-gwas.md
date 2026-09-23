# Work report: the linear algebra the association study needs

The plan `docs/plans/linalg-gwas.md` is under way, on the branch
`plan/linalg-gwas`, in the worktree `.claude/worktrees/linalg-gwas`,
where it started on 23 September 2026. It builds, in
`crates/popnei-linalg`, the seven operations that a genome wide
association study needs: a Cholesky factorization, which is the
factorization of a symmetric matrix that no vector makes negative, and
the solve, the log of the determinant and the inverse that come off it;
the thin QR of a design, which is what fits a linear model to more
individuals than coefficients; the solve against an upper triangular
matrix; and the rank, how many of a design's columns are independent. It
also gives `product` a typed first operand, so that it computes all four
of `a b`, `a b'`, `a' b` and `a' b'`, and adds one case to the crate's
error enum, `Singular`. The spec behind it is `docs/specs/linalg.md`.

Where the plan stands on 23 September 2026: work package 1 is under way,
task 1.1 done and task 1.2 running. Work packages 2 and 3 have not
started.

This report is written as the work goes. Each work package gets a section
below when it is done, with the command that checked each deliverable and
what it gave, what was changed in the plan and why, what the review
found, and what the owner should know. When the plan is done, what the
owner reads first goes at the top of this file.

## Before the first task

The branch starts from `spec/linalg-gwas` at `7ccc709`, which holds the
spec and merges into `main` with no conflict.

Everything the plan asks to be in place is there, checked by running it
in the worktree on 23 September 2026:

| What the plan asks | Command | What it gave |
| --- | --- | --- |
| the workspace is clean | `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings` | both clean |
| the workspace passes | `cargo test --workspace` | `472 passed` in the core crate, `42 passed` in the linalg crate, 2 ignored |
| both backends pass | `cargo test -p popnei-linalg` and the same `--no-default-features` | `42 passed` and `37 passed` |
| both wasm targets build | `cargo wasm-check` | clean |
| the Python layer passes | `uv run ruff format --check && uv run ruff check`, `uv run maturin develop && uv run pytest` | clean, `257 passed` |
| the TypeScript layer passes | `npm run build && npm test` in `js/popnei` | `pass 180` |
| none of the seven exists yet | `cargo test -p popnei-linalg --lib <filter> -- --list` for `cholesky`, `qr`, `rank`, `triangular`, `determinant`, `singular`, `four_ways`, `first_operand`, `solve_with`, `invert` | `0 tests` for every one of the ten |
| `product` has no typed first operand | `grep -rl TheFirstOperand crates/` | no file |
| the Cargo.toml comment still names an open point | `grep -c "Open 1" crates/popnei-linalg/Cargo.toml` | `1` |
| the reference the tasks build from works | `cargo run --release` and the same `--no-default-features` in `tmp/linalg_gwas_trial/` | every check of both backends printed `ok` |

One trap for anybody who runs the checks here. `uv run maturin develop
--release` makes the pytest test
`test_a_ctrl_c_while_write_vars_runs_is_raised_and_leaves_no_file` fail:
it sends the process a ctrl-C 0.1 s into a write that a release build has
already finished. The check of the `coding` skill is the debug build, and
with it all 257 pass.

## Work package 1: the product with its first operand turned

### Task 1.1, the typed first operand and its four callers

`8d7fb27`. `product` takes a `TheFirstOperand` beside its
`TheSecondOperand`, the four calls of it in the core crate name
`ByTheRowsOfTheResult`, and the comment of
`crates/popnei-linalg/Cargo.toml` says the owner's decision about the
rustc flag instead of naming an open point.

Checked by the orchestrator, each command run in the worktree: `cargo fmt
--all --check` clean, `cargo clippy --workspace --all-targets -- -D
warnings` clean, `cargo test --workspace` `472 passed` and `42 passed`
with 2 ignored, `cargo test -p popnei-linalg --no-default-features` `37
passed`, `cargo wasm-check` clean, ruff clean, `uv run maturin develop &&
uv run pytest` `257 passed`. Every one is what it was before the task, so
nothing the principal component analysis or the r² computes moved.
`grep -c ByTheRowsOfTheResult` gives 3 in `pca.rs` and 1 in `ld.rs`, and
`grep -c "Open 1" crates/popnei-linalg/Cargo.toml` gives 0, which are
deliverables 1 and 4. The diff of `pca.rs` and of `ld.rs` is the four
calls and one `use` line and nothing else; no test of either was touched.

One thing the task had to settle. The two combinations whose first
operand is `ByTheValuesSummedOver` have no backend behind them until task
1.2, and the enum lets a caller ask for them, so the task made those two
arms give `Error::Dimension` with a message that says the product is not
built yet, rather than an `unreachable` that would be a wrong matrix or a
panic waiting to happen. Task 1.2 replaces the two arm bodies with one
backend call each. It is the right call for one commit inside a work
package, and nothing outside the crate can reach it: no caller in popnei
names that case.

The subagent used 123480 tokens.
