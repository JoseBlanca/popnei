# Report: the curve of r² against distance, per population

24 September 2026. It records how `docs/plans/ld-vs-dist.md` was carried
out, on the branch `plan/ld-vs-dist` in the worktree
`.claude/worktrees/ld-vs-dist`. The plan builds the item "LD against
distance, per population" of `docs/specs/ld.md`: for each population of a
dataset, how the r² of a pair of variants falls off as the two move apart
along a chromosome, in bins of distance, as a curve fitted to the pairs,
and as the one distance at which r² has fallen to half.

**The plan is under way.** Nothing is merged into `main` and nothing is
pushed.

## What was in place before the first task

The plan's "What has to be in place" was checked by running it, on the
branch at `a04c7eb`, which is `main` at `cce5801` with the two spec
commits and the two plan commits of this branch on top.

| what the plan says | the command | what it gave |
|---|---|---|
| the core and linear algebra suites pass | `cargo test --workspace` | `787 passed; 0 failed; 2 ignored` and `149 passed; 0 failed` |
| 49 ld tests, none of the fall-off | `cargo test -p popnei --lib ld:: -- --list` | `49 tests`, and `0 tests` for both `ld::dist` and `ld::decay` |
| the Python tests of the item do not exist | `uv run pytest tests/test_ld.py -k ld_and_dist` | exit 5, `16 deselected`, which is what pytest gives when `-k` matched nothing |
| no code of the item anywhere | `grep -rn "ld_and_dist\|LdAndDist\|LdBins\|LdDecay" crates python js tests` | nothing |
| plink2 and R are on the machine | `plink2 --version`, `R --version` | `PLINK v2.0.0-a.7.7 M1 (18 Sep 2026)`, `R version 4.6.1 (2026-06-24)` |
| the reference script writes its files again and finds them the same | `tests/reference/ld/run_plink2.sh` into an empty directory | exit 0, nothing of its own printed |

The other checks of the `coding` skill were run at the same commit and are
clean: `cargo fmt --all --check`, `cargo clippy --workspace --all-targets
-- -D warnings`, `cargo test -p popnei --no-default-features` with the same
787, `cargo wasm-check`, `uv run ruff format --check` and `uv run ruff
check`, and `npm run build && npm test` in `js/popnei` with 325 tests
passing.

Two things the plan states that the commands corrected:

- **The npm package had to be built before the Python tests would run.**
  A fresh worktree has no compiled extension module, so `uv run pytest`
  failed to import `popnei` until `uv run maturin develop` had been run
  once. It is a step of the `coding` skill's checks and not a change to
  the plan.
- **`js/popnei/test/ld.test.ts` holds 10 tests and not the 11 the plan
  counts.** `node --test js/popnei/test/ld.test.ts` prints `tests 10`, and
  the file has 10 `test(` at the top level. Deliverable 6 of work package 1
  and deliverable 7 of work package 2 ask that the file then run more ld
  tests than it runs today, so the number to beat is 10.

## Work package 1: the bins of r² against distance

### The tasks as they were done

**1.1, the window over the blocks.** Commit `a0a8181`. It put
`TheWindowOfTheBlocks` in a new `crates/popnei/src/ld/dist.rs`, declared
by one `mod dist;` line added to `crates/popnei/src/ld.rs`. The window
takes a block, holds the blocks whose variants are within `max_dist` of
the newest variant read and on its chromosome, and gives back how many of
the oldest blocks fell out on that call, which is what lets the next task
drop what it keeps beside each block in step.

The plan puts this code in `crates/popnei/src/ld.rs`. It went in
`crates/popnei/src/ld/dist.rs` instead, which the plan's own check asks
for: the task's tests have to be printed by `cargo test -p popnei --lib
ld::dist -- --list`, and that needs a module `ld::dist`. The crate is
edition 2024, so a file beside `ld.rs` is a submodule of it, and the 3859
lines already in `ld.rs` were left alone.

`cargo test -p popnei --lib ld::dist -- --list` prints `14 tests`, where
it printed `0 tests`. `cargo test --workspace` gives 801 passed with 2
ignored, which is the 787 of the starting state plus those 14, and 149 in
the linear algebra crate. `cargo test -p popnei --no-default-features`
gives the same 801 on the faer backend. `cargo fmt --all --check`, `cargo
clippy --workspace --all-targets -- -D warnings`, `cargo wasm-check`,
`uv run ruff format --check` and `uv run ruff check` are clean, and `uv
run pytest` gives 499 passed. Every one of these was run by the
orchestrator and not taken from the task's report.
