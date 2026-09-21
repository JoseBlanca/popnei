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

## Work package 1: what the new specs change in the code that exists

Task 1.1, commit 080b1da: `FilteringStats` alone in the new module
`filters`, and `filtering_stats` in the trait with no default, in `Box`,
`Reblock`, `VcfReader`, `VarsReader` and the two readers of the tests.
Neither binding crate needed a change. Run by the orchestrator: `cargo
test --workspace` `254 passed`, 2 ignored; `cargo test -p popnei --lib --
filtering_stats --list` `5 tests`; fmt, clippy and `cargo wasm-check`
pass. 126414 tokens.

Task 1.2, commits 2160dcf and 6291c58. `iter_blocks` gives a `Blocks`
with `pass_stats`, `write_vars` a `VarsWritten`, and a `Variants` has
`steps` and a `repr`, `<Variants of tests/reference/vcf/many.vcf, no
steps>`. The list of steps lives in the binding crate, in a class
`_core.Steps` that the `Variants` of the package holds, and every pass
builds its chain from it. Run by the orchestrator: `cargo test
--workspace` `257 passed`; `uv run pytest` `115 passed`, 16 of them in
`tests/test_pass_stats.py`; fmt, clippy, ruff and `cargo wasm-check`
pass. One test that was there changed, the one of `tests/test_io_vars.py`
that calls `_core.write_vars` itself, which now takes the steps. 262494
tokens, before the commit of the specs.

For the owner, from task 1.2:

- The core's `write_vars` returns the sink and how many variants it
  wrote, where "The Rust interface" of `docs/specs/io_vars.md` had the
  sink alone, and `BlockReader` is implemented for `&mut R`. The binding
  crate does not loop over the blocks of a write, so only the core can
  count them, and it has to keep the chain to read the counts of the
  filters, as "How it runs" of the counts asks. The orchestrator had the
  two specs corrected, in 6291c58, without stopping: it is a signature of
  the core that no user sees, and the owner's decision that every
  consumer returns its counts cannot be built on the old one. His to
  reverse.
- `copy.copy(variants)` gives a handle that shares the steps of the
  original, so a filter added to the copy shows on the original.
  `docs/specs/filters.md` leaves a copy of a `Variants` out, and nothing
  tests it.

Task 1.3, commit 862ec70: the same in TypeScript. `iterBlocks` gives an
iterator class around the generator that was there, because the pass is
freed when the iteration ends and the counts have to answer after that;
`writeVars` gives `{bytes, passStats}`; `steps` is an empty array. A
pass takes a copy of the list of steps, because wasm-bindgen cannot lend
an optional object of the crate to a function and an owned one leaves
the caller's dead; the copy is also what the spec asks, a pass runs the
steps it started with. `variants.free()` frees the steps too, so `steps`
after it is an `Error`, as `iterBlocks` is. Five lines of tests that
were there changed, four of `vars.test.ts` and one of `web.test.ts`,
where they take the bytes of `writeVars`; `pass.test.ts` is untouched.
273909 tokens.

The deliverables of work package 1, run by the orchestrator at 862ec70:

1. `cargo test -p popnei --lib -- filtering_stats --list`: `5 tests`.
   `cargo test --workspace`: `257 passed`, 2 ignored.
2. and 3. `uv run pytest tests/test_pass_stats.py`: `16 passed`; `uv run
   pytest`: `115 passed`.
4. `npm run build && npm test` in `js/popnei`: `tests 76`, `fail 0`.

fmt, clippy, ruff and `cargo wasm-check` pass.

How the work went. The orchestrator committed the tick of task 3.1 with
`git add <paths>` and `git commit` with no paths while the subagent of
1.1 had its files staged, and the commit took them. The commit was local,
so the orchestrator split it in two, with the subagent's message, and the
tree did not change. Two sessions in one tree commit with `git commit
-F <message> -- <paths>`, and the prompts of the tasks now say so.

The review, at 0a7c9ea, with the seven categories, `spec` and `tests` in
worktrees of their own. What held and is fixed, in 8baf37b, f6d2c30 and
9914e33, each defect with a test that failed first where a test can see
it:

- In TypeScript `numVars` counted a block that the pass lost after the
  reader gave it: four variants with the last at the position
  9007199254740993, in blocks of 1, gave 3 blocks and `numVars: 4`. Three
  reviewers ran it. Python counted right.
- In Python `pass_stats` waited for the lock of the pass with the GIL
  held, so every thread of Python stopped until the block being read was
  read: 0.79 s for a block of 500000 variants. The lock is now taken with
  the interpreter released; on a block of 20000 variants of 500
  individuals the longest stop of a third thread fell from 0.163 s to
  0.020 s. It has no test: what it gives shows only as a timing.
- The steps of a pass were an optional argument of the TypeScript
  binding crate that became no steps when left out, which four reviewers
  reported: with filters, a silent unfiltered pass. It is required now,
  and one line of `pass.test.ts` passes `new Steps()`.
- "How it runs" of the counts in `docs/specs/filters.md` said that the
  binding crate counts the variants of every consumer, against
  `docs/specs/io_vars.md` as corrected. It has the exception of
  `write_vars`. The `spec` and the `api` reviewers found the corrected
  signature of `write_vars` right.
- The README example said 4 variants where it gives 3; a comment in both
  binding crates said that a lost block is in the count of no filter,
  which a filter makes false; `args` of a TypeScript `Step` took numbers
  alone.
- Five tests that were missing: the order of `filtering`, which no test
  failed without, in both packages, with the spec's pairs; `VarsWritten`
  frozen; `steps` after `free()`; `passStats` after a `break`; and a
  block that was read when a Ctrl-C arrived is not in `num_vars`.

Not taken: that a dict keyed by the kind loses a pair of counts when two
filters have one kind, because `FilteredReader::new` refuses the second,
in work package 3; and `Blocks` declared in `variant.py` and not beside
`Block`, which an import cycle forbids.

After the fixes, run by the orchestrator at 9914e33: `cargo test
--workspace` `257 passed`; `uv run pytest` `117 passed`; `npm test`
`tests 79`, `fail 0`; fmt, clippy, ruff and `cargo wasm-check` pass. The
subagent built the wheel of pyodide and its smoke test exited with 0.

For the owner, from the review:

- The two binding crates each hold the same code for the kinds of step
  and for building the chain of a pass, and work package 3 adds the three
  criteria to both. The `architecture` reviewer proposes a function of
  the core that builds the chain from a list of criteria, which "The
  Rust interface" of `docs/specs/filters.md` does not have. The
  orchestrator follows the spec, and recommends the function: it can be
  added after work package 3 at the cost of a small task.
- "Not in this spec" of `docs/specs/filters.md` has the read ahead
  thread take the reader and give it back when the pass ends. While it
  has the reader nobody can ask it for `filtering_stats`, which
  `pass_stats` does mid pass, and a thread that is spawned cannot take
  the chain that `write_vars` is lent, only a scoped one. The reviewer
  tried both in a scratch crate. The item of the read ahead thread will
  need a snapshot of the counts, or the owner's word that the counts are
  read at the end alone.
- With the default `only_passed=True`, `many.vcf` gives 475 of its 500
  variants, and no count of a pass says that the FILTER column took 25
  out. Seen by the `numbers` reviewer, outside the scope.

How the work went, besides the commit above: tasks 1.2 and 1.3 cost
279083 and 273909 tokens, twice the core task; the seven reviewers
114325 to 178423 each; the fixes 70724 more of the subagent of 1.3.

## Work package 2: the counts of one variant

Task 2.1, commits 37ac866, the sentence of `docs/specs/variant.md` that
makes the two new errors a `RuntimeError`, and 10a9b4c, the code. The
subagent ran every literal of its tests through pyNei at ef0ca6e, the six
variants and the three tetraploid genotypes, and they are the spec's.

The deliverables, run by the orchestrator at 10a9b4c:

1. `cargo test -p popnei --lib -- count_gts count_alleles --list`: `7
   tests`, which cover the cases the plan lists. `cargo test --workspace`:
   `266 passed`, 2 ignored.
2. The message of each new case has a test, and
   `crates/popnei-python/src/errors.rs` has both in its `RuntimeError`
   arm. `crates/popnei-js/src/errors.rs` has no arm to add: it turns
   every error of the core into one `Error`.

`uv run pytest` `117 passed`, `npm test` `tests 79`, `fail 0`, fmt,
clippy, ruff and `cargo wasm-check` pass. 159072 tokens.

For the owner: both functions also refuse a variant of more alleles than
a `u32` holds, 4295 million, with the error of genotypes that are not
whole. The spec does not name that case. It is what lets the loop add
without a check of overflow, and no source reaches it. The review looks
at whether the message fits the case.

Changed in the plan: the review of this work package is made together
with tasks 3.2 and 3.3, the filter of the core, which is built on these
two functions and is one piece of code with them, and the two bindings
of work package 3 get a review of their own. One review of five tasks
through three layers would have been too large to evaluate well.

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

Tasks 3.2 and 3.3, one subagent, commits c1f7397, the filter of one
block, alone in its commit, and ffec301, the reader. Run by the
orchestrator at ffec301: `cargo test --workspace` `291 passed`, 2
ignored; `cargo test -p popnei --lib filters:: -- --list` `25 tests`,
the worked example, the nine rows of the table of `many.vcf` with blocks
of 7 and of the default size, the chain with its 106 variants and its
three pairs of counts, 29 of 100 at 0.29, and the three rules of a
reader among them. No number of the spec failed to come out. 226451
tokens.

For the owner, from these two tasks:

- The error of a second filter of one kind, when `FilteredReader::new`
  gives it, names the kind and the threshold of the filter that was
  refused, where the spec says the threshold that is set: a chain says
  only the kinds of its filters, through `filtering_stats`. The binding
  crate has the steps and gives a user the error at the call of the
  method, as the spec says, so it can name the one that is set, and task
  3.4 is told to.
- The error of a threshold out of range names the criterion by its
  kind, `maf`, and the value. The spec has the `ValueError` of Python name
  the argument, `max_allowed_maf`, which only the binding knows, so task
  3.4 adds it there.

The merge of `main`. The owner had `plan/vars-file` merged into `main`
at d1d6997, with 20 commits that this branch lacked, and said to merge
it here. It was made at d4e8458, between tasks 3.3 and 3.4, by the
subagent of the TypeScript side, because its 13 conflicts were in seven
files of that side: `main` had changed how a vars file crosses out of
wasm, in pieces of 1 MiB, in the functions to which this plan had given
the steps and the counts. Both hold: `writeVars` puts the pieces
together and gives `{bytes, passStats}`. Four calls of `writeVars` in
the tests of `main` and one core test were brought to the new
signatures, and no test of either side was dropped: cargo 291 and 1,
npm 79 and 8. Run by the orchestrator at d4e8458: `cargo test
--workspace` `292 passed`, 2 ignored; `uv run pytest` `117 passed`;
`npm test` `tests 87`, `fail 0`; fmt, clippy, ruff and `cargo
wasm-check` pass. The subagent built the wheel of pyodide, its smoke
test exited with 0, and `cargo bench --no-run` built. 65057 tokens.

The review of the counts and of the filter of the core, at f22044f, with
six categories, all but `binding`. No reviewer found a wrong number. The
`spec` reviewer got the nine rows and the chain again from pyNei and
from bcftools, the missing data filter from plink2, and the same keep or
drop as pyNei's three functions for a tetraploid variant, a tie for the
major allele, a half called genotype at a ploidy of 3, the allele 127
and four thresholds that the table lacks. The `tests` reviewer made 24
mutations of the code and each failed a test. What held and is fixed, in
61656c6, the two specs, and 4f8bf3c, the code:

- `count_alleles` added into an array that its caller had to clear, and
  nothing checked it: with an entry at `u32::MAX` a release build
  wrapped to 2 with no sign. Its only caller built a new array for every
  variant, where the spec has one handed over again. The function clears
  the array itself now, and the filter gives each thread one that is
  reused. `docs/specs/variant.md` says so. The owner's to reverse.
- The error of a second filter of one kind held the refused threshold
  alone. It has `threshold_that_is_set` too, empty from the core, which a
  binding crate fills, and its message no longer names the steps of the
  variants, a name the core lacks.
- Under rayon the error of a block with two bad rows was the one a
  thread stored first: the allele of row 3 in 190 of 200 runs of a
  reviewer, and the one of row 390 in 10. It is the error of the lowest
  row now, found by reading the rows again one by one when the threads
  give an error.
- A variant of more alleles than a `u32` holds was the error of
  genotypes that are not whole, with a ploidy of 1 that nobody gave. It
  is a case of its own, `MoreAllelesThanACountHolds`, in the spec too.
- Six of the nine rows of the table asserted how many variants stay and
  not which, and nothing read `tests/reference/filters/`. Every row and
  the chain compare every position with the stored files now. The test
  that the ploidy goes through a filter could not fail, its source was
  diploid; it has a ploidy of 4.
- A source narrowed before a filter is put over it, with no `set_needs`
  on the chain, fails at its first block for lack of genotypes. The doc
  comment of `FilteredReader::new` says that the needs are set on the
  outermost reader once the chain is built, and a test states it. The
  other fix was a method of the trait.
- The product of the individuals and the ploidy written twice; a comment
  about a 0/0 that cannot be reached; two `# Errors` that named half
  their errors; `threshold()` of a criterion made public for the `args`
  of a step, and in the spec; a warning of `cargo doc`.

After the fixes, run by the orchestrator at 4f8bf3c: `cargo test
--workspace` `296 passed`, 2 ignored; `filters:: -- --list` `27 tests`;
the counts `7 tests`; `uv run pytest` `117 passed`; fmt, clippy, `cargo
wasm-check` pass and `cargo doc -p popnei --no-deps` has no warning.
The six reviewers cost 113445 to 161289 tokens each, the fixes 71497.

For the owner, from this review:

- A vars file with a damaged byte gives an allele of -2 with no error.
  The `errors` reviewer changed one allele byte of a file that popnei
  wrote to 0xFE, and `open_vars` gave `[[-2, 0], ...]`: the reader checks
  the type and the width of the genotypes and not their values. So the
  sentence of `docs/specs/variant.md` that no reader gives such an
  allele is false. With a filter on, the counts refuse it, as a
  `RuntimeError`, which tells a user to report a defect when their file
  is what is wrong; with none, the -2 reaches their array. The options:
  the vars file reader refuses it, as a `ValueError` with the file, which
  costs a scan of every genotype byte read, not measured, and a case in
  `docs/specs/io_vars.md`; or the error of the counts becomes a
  `ValueError` and the silent case stays. The orchestrator recommends
  the first, as a task of its own, since the vars file is outside this
  plan. Asked in chat on 21 September 2026; meanwhile a `RuntimeError`.
- A filter that keeps nothing over a long stretch of a file reads on
  inside one `next_block`, so one `__next__` of Python can read a whole
  file and a Ctrl-C waits for it.
- Outside the plan: the test of `crates/popnei/src/io/vars.rs` that
  sweeps one changed byte installs a global panic hook, so any other
  test that fails makes it fail too, and names the vars file to whoever
  broke a filter. Counting the panics of its own thread alone would
  settle it.
