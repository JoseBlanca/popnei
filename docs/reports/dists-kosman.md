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

The deliverable, run by the orchestrator at a73714d, after the review:
`uv run python tests/reference/dists/make_reference.py` takes 17 s,
prints that the 14 literals of the spec's table are what R gives and
that a loop in Python over the genotypes gives the two integers of every
one of the 20812 pairs of the four datasets, and `git status --short
tests/reference/dists` shows nothing after it. Work package 1 finished as
planned, and nothing in the plan changed.

The review, the `spec` and `tests` categories at c66e4f8, found seven
things, all fixed in a73714d by the subagent of the task:

- Both reviewers, from different sides: the tetraploid and the haploid
  files rested on 3 of their 66 pairs, the literals of the spec's table,
  because the loop in Python that the spec's "How it is verified" names
  had not been carried over from `poly_export.py`; the tests reviewer
  scrambled R's output past the third pair and the script stored it
  without a word. The loop is in the script now and every pair of the
  four datasets is checked against it, the integers exactly and the
  distances by 0.0.
- The integer the cargo tests assert, k times the sum of d, was not
  stored, and the docstring's way to get it, the distance times k times
  n, truncates to the wrong whole number for 969 of the panel's 19900
  pairs, 13 of the 780 and 10 of the haploid 66: h00, h03 gives
  112.99999999999999 where the table says 113. The files have a third
  column, `k_sum`, the product rounded after the check that it is within
  1e-10 of a whole number.
- The comparison with pyNei could not fail once a pair was NaN, since
  `max` over NaN gives NaN and no NaN is above a threshold; no pair of
  the four datasets has n = 0, so it was latent. R's version was not
  checked where adegenet's and PopGenReport's were: an R that is not
  4.6.1 is refused now. And three messages and two sentences of the
  docstring said what the script did not do.

No finding was set aside. Tokens: 84281 for the `spec` reviewer, 84078
for the `tests` reviewer, 227523 for the subagent over the task and the
fixes.

## Work package 2: the calculation, through the three layers

Task 2.1, commit 576499c: the sets of bits of a block and the counts of
a pair, a `pub(crate)` type `KosmanBits` in the new module `dists` with
the `called` sets of every individual in one array of 64 bit words and
the `holds` sets, the allele a at count m at a * k + m - 1, in another,
each individual's sets side by side so that a pair is two passes over
contiguous slices; and a new case of the error of the crate, a sum that
would go above `u32::MAX`, which the accumulation of task 2.2 reuses.
Run by the orchestrator: `cargo test -p popnei --lib -- dists:: --list`
12 tests, the three worked examples among them; `cargo test --workspace`
`318 passed`, 2 ignored; fmt, clippy and `cargo wasm-check` pass. The
subagent broke the code to count matches at one copy only and saw the
tetraploid test fail. Its probe on one block of 5000 variants x 1000
individuals, biallelic, one thread, on this machine: 14 ms to build the
sets and 30 ms for the 499500 pairs, against the trial's 15 and 29 ms; a
first shape of the pair loop, take and skip over a zip, took 86 ms. 197707
tokens.
