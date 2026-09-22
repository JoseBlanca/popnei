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

The deliverables of work package 2, run by the orchestrator at b9ced15:
`cargo test -p popnei --lib -- dists:: --list` 33 tests, the worked
examples, the 14 literals, the four files pair by pair, the blocks of 7,
64, 65 and 300 variants and the pools of 1 and 4 threads, the three
errors among them, where the check asks 14; `uv run pytest
tests/test_dists.py` `29 passed`; `npm test` `tests 139`, `fail 0`, with
`test/dists.test.ts` running the panel's five literals, the worked
example and the error of a source with no variant. The work package
finished as planned. One thing was added to the plan's text, the third
error of task 2.2, above.

The review, all seven categories over 576499c, 4834445, 49e99c6 and
b9ced15, `spec` and `tests` in worktrees of their own. What it found
that mattered:

- A wrong number in silence, from `errors`: `Distances.from_square_dists`
  took the names from the index and the values from the upper triangle
  and never looked at the columns, so a frame whose rows a user had
  reordered gave d(c, a) as 0.0 where it is 0.2.
- Three reviewers, `api`, `binding` and `errors`, from different sides:
  the message of a pass with no variant was worded differently in Python
  and in TypeScript, both within the spec, and the TypeScript one dropped
  the filters' counts when the source had none. And, from `tests`, the
  two pytest tests of that message could not fail, because the tmp path
  of pytest already held every word they looked for: the branches
  swapped and the counts swapped left the suite green. The two crates
  now build one sentence, and the tests assert the whole of it.
- From `architecture`: the sets of a block were sized by the largest
  allele value, not by the alleles the block holds. On 3000 diploid
  variants of 1000 individuals, one variant with the alleles 0 and 127
  beside biallelic ones took 7.74 s and 231 MB of peak memory where the
  biallelic file took 0.75 s and 134 MB, every pair walking 10240 words
  instead of 160; a wrong cost and not a wrong number. The alleles of a
  block are ranked to a dense index now.
- From `tests`: three error branches of the sets' builder that no test
  reached, one of which is the only guard against an allele below the
  missing one being dropped while its variant counts as called, a wrong
  count in silence. From `spec`: the error of a panel too large was
  reached by no test, and its message for more pairs than a `usize`
  counts blamed the machine's memory; the TypeScript tests were all
  diploid and had no case of a pair with no called variant, nor of the
  `Variants` being as it was after the call; and the test of
  `pass_stats` with a filter used one that keeps everything, so a
  `num_vars` counting the source's variants would have passed.
- From `numbers`: the message for a vector of a wrong length gave two
  counts both below the length, because the number of individuals
  around it was computed one too low; and two doc comments stated
  bounds and claims that do not hold as written, the pairs' comment on
  bit equality with R and pyNei holding only for a ploidy that is a
  power of two, since those sum d and divide once, which is exact for k
  = 2 and 4 and differs in the last bit from popnei's one division for
  k = 3 in 14350 of 20000 random pairs.
- From `errors` and `binding`: three wrong inputs of `Distances` came
  out as numpy's or Python's own message; the dtype was not made right
  before the array was held.
- From `api`: the README of the TypeScript package did not list
  `minNumSnps` among the checked arguments nor show the new function;
  messages that count 1 individuals; a test only method of the sets with
  a name one letter from another.
- From task 3.2: pandas was imported by the new module and not declared
  in `dependencies` of `pyproject.toml`, so `import popnei` failed under
  pyodide; it is declared now.
- Also taken: the sets of a block dropped before the next block is read;
  the comment on the pass saying that a Ctrl-C is raised when it
  returns.

Not taken, with the reason:

- `numbers` suggested asserting the reference distances bit for bit in
  the cargo tests, since all 20812 pairs are; the spec and the plan ask
  1e-9 there, and the exact comparison is the pytest one with pyNei.
- `binding` and `spec` noted that the pass checks the interpreter's
  signals once before it starts, so a Ctrl-C waits for the end of the
  pass, as it does in `write_vars`; nothing is lost, and the report of
  the filters plan already put the question to the owner for a spec of
  its own.
- `architecture` found that `cargo wasm-check` checks the core alone, so
  the wasm binding crate's new file is compiled for wasm only by `npm
  run build`. The alias also builds for emscripten, where that crate
  does not compile, because wasm-bindgen wants `RefUnwindSafe` of the
  `dyn BlockReader` it holds; a second alias for `wasm32-unknown-unknown`
  alone would check it, and the plan runs `npm run build` for every task
  that touches the crate. Left as it is, for the owner.

What the review checked and found right is worth a line: the three
worked examples at the three layers; the 14 literals and the 20812 pairs
of the four files; pyNei's four tests; 24 random datasets of ploidies 1
to 4 with half called genotypes against an independent loop, bit for
bit; twelve mutations of the core each failing at least one test; every
allocation asked with `try_reserve_exact`; no `unsafe`, no new
dependency, rayon behind its `cfg` with the serial path beside it and
never with the interpreter held; the vector handed to numpy without a
copy.

Tokens of the reviewers: `spec` 162229, `tests` 129376, `numbers` 89184,
`errors` 132948, `api` 128888, `architecture` 108887, `binding` 86064.

The fixes, one commit per layer by the subagent that wrote it: 1cbd94e
the core, b39c30c the TypeScript side, a996652 the Python side. The
subagent of the core, on the probe of one block of 5000 variants x 1000
individuals, biallelic, one thread, release profile, best of 3: the pairs
take 0.030 s as before the ranking of the alleles; 3000 variants of
1000 individuals with one variant of the alleles 0 and 127 take 0.025 s
for the pairs against 0.018 s biallelic, 282 words of `holds` per
individual against 188, where it had 10240. Two things the subagent of
the Python side decided and the owner may reverse: `Distances([],
names=[])` is taken, with no individual, since an empty vector is what
one individual gives and what none would give; and the messages name the
type of what was given without an article, as the Rust ones do.

Run by the orchestrator at a996652, after the fixes: `cargo fmt --all
--check` exit 0; clippy no warning; `cargo test --workspace` `345
passed`, 2 ignored; `cargo wasm-check` finished; `cargo test -p popnei
--lib -- dists:: --list` 39 tests; ruff `20 files already formatted` and
`All checks passed!`; `uv run maturin develop && uv run pytest` `208
passed`, `tests/test_dists.py` `34 passed`; `npm run build && npm test`
`tests 144`, `fail 0`. The wheel for pyodide built, and its smoke test
then failed: the floor of pandas had been taken from `uv.lock`, 3.0.6,
and pyodide 314.0.7 ships 3.0.2, which micropip refuses to go past. The
floor is pyodide's now, as numpy's is, in 6647303, and at that commit
the wheel builds, `node tests/pyodide/smoke.mjs` exits with 0 and
`pytest` gives `208 passed`.

How the work went. Two subagents fixing side by side in one tree each
had a commit swallow the other's staged files once, with `git commit`
without paths, and both undid it with a soft reset and no loss; the form
that stages and commits in one call, `git commit -F - -- <paths>`, is
the one that closes the window. Tokens of the fixes: 275438 for the
core's subagent over its task and the fixes, 244672 the TypeScript
one, 309568 the Python one.
