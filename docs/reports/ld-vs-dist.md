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
