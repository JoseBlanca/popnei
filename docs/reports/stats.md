# Work report: the stats module and the filter of individuals

The plan `docs/plans/stats.md` is under way, on the branch `plan/stats`,
in the worktree `.claude/worktrees/stats`, since 22 September 2026. The
orchestrator, in this report, is the session of the assistant that runs
the plan: it sends each task to a subagent on Opus, checks what comes
back and has each work package reviewed.

## Before the first task

The owner approved the plan in chat on 22 September 2026. The branch
stands on `spec/stats` at 7bea42b, which holds the two specs, the plan
and `tests/reference/stats/`, and which is not in `main`; `spec/stats`
is from `main` at 7d8366f.

Run by the orchestrator in the worktree at 7bea42b: `cargo fmt --all
--check` exit 0; `cargo clippy --workspace --all-targets -- -D warnings`
no warning; `cargo test --workspace` `306 passed`, 2 ignored; `cargo
wasm-check` finished; ruff `18 files already formatted` and `All checks
passed!`; `uv run maturin develop && uv run pytest` `174 passed`; `npm
install`, `npm run build` and `npm test` in `js/popnei` `tests 126`,
`fail 0`; after `npm install` in `tests/pyodide`, `bash
scripts/build_pyodide_wheel.sh` built the wheel and `node
tests/pyodide/smoke.mjs` exited with 0. `which plink2` gives
`/opt/homebrew/bin/plink2`, `PLINK v2.0.0-a.7.7 M1 (18 Sep 2026)`, and
`which bcftools` `/opt/homebrew/bin/bcftools`, version 1.24. pyNei's
`calc_per_var_distribs`, `calc_per_sample_stats`, `filter_samples` and
`load_vars` import under `uv run`, and
`/Users/jose/devel/pynei/test/gwas_reference/sim_missing.vars` is there.
`uv run python tests/reference/stats/make_reference.py` printed `done`
and exited with 0; after it `git status --short tests/reference/stats`
showed `panel.vcf.gz` modified, and only in the two bytes of the
modification time of its gzip header, bytes 5 and 6, the text of the
file being the same: the committed file was put back with `git
checkout`. The bcftools commands of "How it is verified" of the filter
of individuals on `many.vcf` give `ind05 ind00 ind49`, 423 variants with
the positions 1000, 1037, 1074, 1111 and 1148 first, and 26 with the
filter over the 50 individuals. `/Users/jose/devel/popnei-bench/big.vars`
is there, 81356714 bytes, and `big.vcf` beside it, 403572954 bytes.

The start message of the branch is on the board,
`.claude/board/2026-09-22T1020-plan-stats.md`.
