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

## Work package 1, while it goes

Task 1.1, commit 71279e5, one subagent run of 82 thousand tokens and 3
minutes. The orchestrator ran the checks again: `cargo fmt --all
--check` exit 0, clippy finished with no warning, `cargo test
--workspace` `2 passed`, and `cargo test -p popnei --lib -- --list`
`2 tests`. The subagent added a test that gzips a line and reads it
back, so that flate2 is called and not only compiled, and built the core
for both wasm targets: `wasm32-unknown-unknown` and
`wasm32-unknown-emscripten` both finish, with miniz_oxide and no C
compiler. The plan's check `cargo test -p popnei -- --list` ends with the
line of the doc tests, `0 tests`; the count of the library is read with
`--lib`.

Change to the plan: tasks 1.2 and 1.3 were marked as side by side, and
they run one after the other. Both add a member to the root
`Cargo.toml` and both write `Cargo.lock`, and a `cargo clippy
--workspace` of one would meet the half written crate of the other. One
tree has one writer per file.

For the owner: the manifests have no `license` field, because no
document of the repository names a license. Nothing needs it until
popnei is published.
