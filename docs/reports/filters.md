# Work report: the three threshold filters and the counts of a pass

The plan `docs/plans/filters.md` is under way, on the branch
`plan/filters`, in the worktree `.claude/worktrees/filters`, since 21
September 2026. The orchestrator, in this report, is the session of the
assistant that runs the plan: it sends each task to a subagent on Opus,
checks what comes back and has each work package reviewed.

## Before the first task

The owner approved the plan in chat on 21 September 2026. The branch
stands on `spec/filters` at fee97bf and on `plan/vars-file` at 2f2577e,
merged at 2c2a630, on his word of the same day; neither is in `main`.

Run by the orchestrator in the worktree at 877107f: `cargo fmt --all
--check` exit 0; `cargo clippy --workspace --all-targets -- -D warnings`
no warning; `cargo test --workspace` `249 passed`, 2 ignored; `cargo
wasm-check` finished; ruff `15 files already formatted` and `All checks
passed!`; `uv run maturin develop && uv run pytest` `99 passed`; `npm run
build` and `npm test` in `js/popnei` `tests 62`, `fail 0`; `bash
scripts/build_pyodide_wheel.sh` built the wheel and `node
tests/pyodide/smoke.mjs` exited with 0, after an `npm install` in
`tests/pyodide`, which a new worktree lacks and the plan now says.
`which bcftools` gives `/opt/homebrew/bin/bcftools`, version 1.24, and
its three commands keep on `many.vcf` the numbers of the spec, as the
plan records. pyNei's `filter_by_maf` and `gather_filtering_stats`
import. `/Users/jose/devel/popnei-bench/big.vcf` is there, 403572954
bytes.

## Work package 1, while it is under way

Task 1.1, commit 080b1da: `FilteringStats` alone in the new module
`filters`, and `filtering_stats` in the trait with no default, in `Box`,
`Reblock`, `VcfReader`, `VarsReader` and the two readers of the tests.
Neither binding crate needed a change. Run by the orchestrator: `cargo
test --workspace` `254 passed`, 2 ignored; `cargo test -p popnei --lib --
filtering_stats --list` `5 tests`; fmt, clippy and `cargo wasm-check`
pass. 126414 tokens.

How the work went. The orchestrator committed the tick of task 3.1 with
`git add <paths>` and `git commit` with no paths while the subagent of
1.1 had its files staged, and the commit took them. The commit was local,
so the orchestrator split it in two, with the subagent's message, and the
tree did not change. Two sessions in one tree commit with `git commit
-F <message> -- <paths>`, and the prompts of the tasks now say so.

## Work package 3, while it is under way

Task 3.1, the reference script, ran beside task 1.1, since it writes
only under `tests/reference/filters/`. Commit d22c5c9. It stores one file
of positions for each set of filters, named by kind and threshold, the
chain as `missing_data_0.04+maf_0.8.txt` and
`missing_data_0.04+maf_0.8+obs_het_0.5.txt`. Checked by the orchestrator:
`wc -l` of the eleven files gives 26, 215 and 455, 35, 384 and 480, 22,
79 and 369, and 163 and 106 for the chain, the numbers of the spec, and
the first three positions of the 106 are 1111, 1407 and 1518. The ruff
configuration of the project leaves `tests/reference` out, so the script
is checked by its path. 99003 tokens.
