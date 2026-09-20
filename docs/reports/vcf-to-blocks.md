# Work report: from a VCF to blocks

The report of the run of `docs/plans/vcf-to-blocks.md`, on the branch
`plan/vcf-to-blocks`, in the worktree `.claude/worktrees/vcf-to-blocks`.
It is written while the work goes, a section after each work package.
The run started on 20 September 2026, on the owner's Apple M5 Pro. The
subagents that write the code run on Opus, as the owner ordered, and the
orchestrator and the reviewers on the model of the session.

## Before the first task

Everything of "What has to be in place" of the plan was checked by
running it, on main at 134c184: `emcc --version` 5.0.3 after sourcing
`~/devel/emsdk/emsdk_env.sh`; `~/devel/pyodide-venv/bin/pyodide` on a
Python 3.14.5 with the GIL enabled, with the cross build env of 314.0.7
installed; `wasm32-unknown-unknown` and `wasm32-unknown-emscripten`
among the installed targets; `wasm-bindgen 0.2.128`; rustc 1.98.0,
uv 0.12.15, node 26.8.2; bcftools 1.24, bgzip and plink2 in
`/opt/homebrew/bin`; pyNei at ef0ca6e with `uv run pytest
test/test_vcf.py` giving `9 passed`; and
`tests/reference/vcf/make_reference.py`, which wrote its files again
with no change for git.

The owner has not yet said how pyNei is depended on, a path or a git
commit. Task 4.2 follows the plan's meanwhile, the absolute path.
