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

The board of the sessions, `/Users/jose/devel/popnei/.claude/board/`,
outside git, was told to this session by another one after task 1.1; the
start message of this branch was posted at 11:45 UTC. The board says that
the branch `spec/stats` changes the signature of `chain_of`, which task
2.3 calls with the signature `main` has today: the two branches meet in
the binding crates at the merge, and the owner orders the merges.

## Work package 1: the reference

Task 1.1, commit c66e4f8: `tests/reference/dists/make_reference.py` and
ten files beside it, the four datasets of "How it is verified" as gzipped
VCFs, for each the distance and the n of every pair from `gd.kosman` in a
tab separated file with a header, and for the two diploid ones pyNei's
distances. Run by the orchestrator: the script, 4.5 s with adegenet
installed, prints that the 14 literals of the spec's table are what R
gives, and `git status --short tests/reference/dists` shows nothing
after it; `cargo test --workspace` `306 passed`; `uv run pytest` `174
passed`; fmt, clippy and ruff pass. 178322 tokens.

For the owner, from task 1.1:

- pyNei and `gd.kosman` are 0.0 apart on the panel's 19900 pairs and on
  the 780 of the 4 allele dataset. The spec's "How it is verified" says
  2.8e-17 on the panel, measured by `ref.R` of the report, where the
  values went through a text file; the script gives R and pyNei the same
  genotypes and compares the vectors in Python. The spec's tolerance of
  1e-15 holds either way, and the spec was not changed.
- Where the PopGenReport source lives: the script downloads the tarball
  of 3.1.3 from CRAN to `~/.cache/popnei/reference/` when it is not
  there, and reads `DESCRIPTION` and `R/gd.kosman.r` from it. It refuses
  another version, and another adegenet, with a message that names what
  it found; both refusals were tried by the subagent.
- `pyproject.toml` excludes `tests/reference` from ruff, so `uv run ruff
  check` never sees the reference scripts; the script was kept clean by
  passing its path. Installing adegenet 2.1.11 compiled about 40 packages
  in about 25 minutes.
