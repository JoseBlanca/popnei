# Work report: the Kosman distances between individuals

The plan `docs/plans/dists-kosman.md` is under way, on the branch
`plan/dists-kosman`, in the worktree `.claude/worktrees/dists-kosman`,
since 22 September 2026. The orchestrator, in this report, is the session
of the assistant that runs the plan: it sends each task to a subagent on
Opus, checks what comes back and has each work package reviewed. The
Kosman distance of two individuals is the share of their genotypes that
differ, averaged over the variants at which both are called, as
`docs/specs/dists.md` has it.

## Before the first task

The owner approved the plan in chat on 22 September 2026, on branch
`plan/dists-kosman` at 5b43ed4, which holds the spec, the plan and
`docs/reports/kosman-method/`.

Run by the orchestrator in the worktree at 5b43ed4 on 22 September 2026:
`cargo fmt --all --check` exit 0; `cargo clippy --workspace --all-targets
-- -D warnings` no warning; `cargo test --workspace` `306 passed`, 2
ignored; `cargo wasm-check` finished; ruff `18 files already formatted`
and `All checks passed!`; `uv run maturin develop && uv run pytest` `174
passed`; in `js/popnei`, after `npm install`, `npm run build` and `npm
test` `tests 126`, `fail 0`. `which Rscript` gives
`/opt/homebrew/bin/Rscript`, version 4.6.1, with neither adegenet nor
PopGenReport installed, as the plan says; task 1.1 puts adegenet in.
pyNei's `calc_pairwise_kosman_dists` and `load_vars` import, and
`sim_missing.vars` is at the path the plan gives, 95154 bytes. `big.vcf`
is 403572954 bytes and `big.vars` 81356714, as the plan says. The
machine has 18 cores.
