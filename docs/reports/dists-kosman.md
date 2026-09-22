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

Task 2.2, commit 4834445: `calc_kosman_sums` and `KosmanSums` as "The
Rust interface" gives them. The pairs of a block go to rayon one row of
the upper triangle at a time, an individual's pairs with those after it,
with the accumulator cut into one slice per row so that no thread takes
a lock; the serial version beside it is what wasm runs and what the test
of the threads compares against. Run by the orchestrator: `cargo test -p
popnei --lib -- dists:: --list` 33 tests, where the deliverable asks 14;
`cargo test --workspace` `339 passed`, 2 ignored; fmt, clippy, `cargo
wasm-check`, ruff and `pytest` `174 passed` clean. All 20812 pairs of the
four reference files match `gd.kosman` on both integers, and the 14
literals within 1e-9. 216350 tokens.

For the owner, from task 2.2, his to reverse:

- A third error case that the spec did not name: the accumulator of two
  `u32` per pair is asked of the machine with `try_reserve_exact`, and a
  panel it has no memory for is an error, a `ValueError` in Python, as a
  block too large for the machine is in `docs/specs/block.md`. The
  orchestrator wrote it into "The Rust interface" of the dists spec in
  this commit rather than stop: no value and no signature changes.
- A block whose individuals or ploidy differ from what the reader says
  its source has is refused with the error the readers already have for
  blocks that do not fit together, so that a reader with a defect cannot
  line the pairs up wrongly in silence.

Tasks 2.3 and 2.4 ran side by side, on different files of the one tree.

Task 2.3, commit b9ced15: `calc_pairwise_kosman_dists` and `Distances`
in `python/popnei/dists.py`, exported from `popnei`, over one function
of the binding crate that builds the chain with `chain_of`, runs the
calculation with the interpreter released and hands the vector over as
a read only float64 array without a copy; and `tests/test_dists.py`, 29
tests. Run by the orchestrator: `uv run pytest` `203 passed`, from 174;
`uv run pytest tests/test_dists.py` `29 passed`; fmt, clippy, `cargo
test --workspace` `339 passed`, `cargo wasm-check`, ruff `20 files
already formatted` clean. popnei and pyNei agree bit for bit over the
19900 pairs of the panel and the 780 of the 4 allele dataset, with
`min_num_snps` of `None` and of 1125, NaN in the same places. 243988
tokens.

Task 2.4, commit 49e99c6: `calcPairwiseKosmanDists` and `Distances` in
`js/popnei/src/dists.ts`, exported from the node and the web entry
points, over a function of the wasm binding crate that does what the
Python one does, and `js/popnei/test/dists.test.ts`, 13 tests. Run by the
orchestrator: `npm run build && npm test` in `js/popnei` `tests 139`,
`fail 0`, from 126. 197836 tokens.

For the owner, from tasks 2.3 and 2.4, his to reverse:

- `Distances(vector)` with no names gives the names 0 to N-1 as
  integers, as pyNei does, so `names` is a tuple of `int` there and of
  `str` for a calculation. The spec says "the names 0 to N-1" and not
  which; the subagent followed pyNei.
- `min_num_snps` is refused by the Python package and not by the binding
  crate, so that a negative value, 2.5, `True` and a value above
  4294967295 are refused with the name the user wrote and not pyo3's
  `OverflowError` on the Rust argument. The last bound is the largest n
  a pair can have, since the core refuses sums above a `u32`.
- The array a user gives `Distances` is held and not copied, where
  pyNei copies: writing into it changes the result. Copying would have
  made the 400 MB of 10000 individuals 800 MB. It is in the docstring.
- The binding crates got one error case each for a pass with no
  variant, which carries the path and the filters' counts, since no
  case of theirs did; the subagent of 2.3 changed one sentence of
  `.claude/skills/coding/pyo3.md`, which counted the cases of the Python
  one.
- Two errors of the calculation, a sum above `u32::MAX` and an
  accumulator the machine has no memory for, reach Python and
  TypeScript through the arm that turns any other error of the core into
  a `ValueError` or an `Error`, and no test outside the core reaches
  them: one needs thousands of millions of variants and the other a
  machine that refuses hundreds of megabytes.
