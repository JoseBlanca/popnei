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

Task 1.2, commit 63b6f66, one subagent run of 116 thousand tokens and 9
minutes. The orchestrator ran the five checks again: fmt exit 0, clippy
with no warning, `cargo test --workspace` `2 passed`, ruff `2 files
already formatted` and `All checks passed!`, `maturin develop` `Installed
popnei-0.1.0` and pytest `1 passed`. Four things it decided that the task
did not say, each with its reason in the commit message:

- The binding crate has `test = false` and `doctest = false`. `cargo test
  --workspace` linked a test binary of it against the Python of the PATH,
  Apple's 3.9.6, and failed at the link. The crate holds no calculation
  and its tests are the pytest ones; clippy still checks it.
- `.python-version` says 3.14.5 and not 3.14. On this Mac uv answers
  `3.14` with the free threaded 3.14.7, the newest 3.14 it has, the same
  trap `docs/rust_core.md` records for pyodide-build.
- `[tool.uv] package = false`, so that maturin alone installs popnei. An
  explicit `uv sync` takes the module out, and `uv run maturin develop`
  has to follow it.
- ruff is given `python/`, `tests/*.py` and `pyproject.toml`. Over the
  whole repository ruff 0.16.8 also formats the Python inside Markdown
  fences, 45 files, 32 of them under `.claude/` and 9 under `docs/`.

It also confirmed what `.claude/skills/coding/pyo3.md` left to confirm:
maturin builds the module without the `extension-module` feature of
pyo3 and without a warning.

For the owner: the manifests have no `license` field, because no
document of the repository names a license. Nothing needs it until
popnei is published.
