# Plan: the filter of the first n variants and the filter that keeps variants at random

5 October 2026. State: **done** on 5 October 2026, approved by the owner that day. It builds the two filters that
`docs/specs/filters.md` gained on 5 October 2026 from issues 7 and 6 of the
repository, "The filter that keeps the first n variants" and "The filter
that keeps variants at random", with what they add to "The Rust
interface" of that spec and the field `stopped_early` of `PassStats` in
`docs/specs/variant.md`. Each runs through the core crate, both binding
crates, the Python package and the TypeScript package. The spec has no open
point: the owner decided its four on 5 October 2026, approved it that day,
and approved the breakdown below the same day.

## In and out

Built: the two work packages below, and the review of the whole branch in
the seven categories of the `code-review` skill after them, which the owner
asked for besides the review of each work package.

Not built, with where each goes: a sample of exactly n variants spread over
the file, and a new sample at each call without a seed ("Not in this spec"
of `docs/specs/filters.md`); the time the random filter costs, which the
spec leaves unmeasured ("Speed" of that spec).

## What has to be in place

- The branch `filters-random-first-n`, in the worktree
  `.claude/worktrees/filters-random-first-n`, from `main` at eae29a2, with
  the spec at d4a8699, and its start message on the board,
  `.claude/board/2026-10-05T0403-filters-random-first-n.md`. The branch
  keeps that name, which the owner gave, and not `plan/<name>`.
- The shared files it changes, which other branches may also change: the
  error enum `crates/popnei/src/error.rs`, `PassStep` and `chain_of` of
  `crates/popnei/src/filters.rs`, the `steps.rs` and every module that
  returns counts in both binding crates, `python/popnei/variant.py` and
  `filters.py`, `js/popnei/src/variant.ts` and `filters.ts`.
- The checks of the `coding` skill on d4a8699's code, which is eae29a2's,
  run on 5 October 2026: the six cargo commands pass, `cargo test
  --workspace` with 1373 tests passed and 2 ignored, summed over its `test
  result` lines, and `cargo test -p popnei --no-default-features` with 1223
  passed and 2 ignored; `uv run maturin develop && uv run pytest`
  gives `669 passed`; `npm ci && npm run build && npm test` in `js/popnei`
  gives `tests 496`, `pass 495`, `fail 1`. The one that fails is `a
  kinship that does not tell the two variances apart gives none of them`
  of `test/gwas.test.ts`, 1.1e-14 against the 7e-15 its tolerance allows,
  which fails on `main`, is on the board, and is not of this plan. "No
  fewer" below counts from those numbers. A fresh worktree needs `npm ci`
  and `npm run build` before `npm test`.
- bcftools 1.24 at `/opt/homebrew/bin/bcftools`, and Java at
  `/opt/homebrew/opt/openjdk/bin/java`, OpenJDK 26.0.2.1, which is not on
  the PATH; `java tests/reference/filters/SplitMix.java` and
  `tests/reference/filters/random_draws.py` were run with them on 5 October
  2026 and gave the numbers of the spec.

## 1. The filter of the first n variants

What it gives: `variants.filter_first_n(1000)` and
`variants.filterFirstN(1000)` make every calculation read the first 1000
variants that the steps before keep and stop reading the file there, and
every `pass_stats` says, in `stopped_early`, whether the filter ended the
pass.

Deliverables:

1. Cargo tests whose names contain `first_n`, in `filters`, with the
   checks of "How it is verified" of the filter of the first n: the ten
   positions and the ten after the MAF filter of bcftools, in blocks of 7
   and of the default size; the counts 14 and 10 and the two blocks asked
   of the source; n of 14 over blocks of 7; the source that never ends,
   read on one thread and through `with_one_block_ahead`; the VCF of 100000
   variants of 10 individuals read through a counting `Read`; the vars file
   in batches of 100 read in a pool of one thread, which decompresses one
   batch; `stopped_early` of the core over each of those passes and over a
   pass that ends with fewer than n; the refusals of n of 0, of a second
   filter of this kind and of each step that takes variants out after it,
   by `refuse_a_step_after_the_first_n` and by `chain_of`, and the filter
   of individuals accepted. `cargo test -p popnei --lib first_n -- --list`
   counts them, and counted 0 on d4a8699.
2. Every existing test passes, no fewer than above, and the only test files
   of before that change are those that assert the whole of a `PassStats`,
   which now holds `stopped_early: False`, and the lists of the kinds of a
   step in `js/popnei/test/`.
3. Pytest tests whose names contain `first_n`, at `filter_first_n`, with
   the checks of its "How it is verified", on `many.vcf` opened with
   `only_passed=False`, among them `stopped_early` in the `repr` of a
   `PassStats`. `uv run pytest -k first_n` exits with 5 on d4a8699, which
   is pytest's code for no test collected.
4. Node tests whose names contain `first n` or `filterFirstN`, with the
   checks the spec gives the TypeScript test, and `stoppedEarly` in the
   counts of every consumer of the list of `test/consumers.ts`. `npm test
   -- --test-name-pattern=filterFirstN` in `js/popnei` runs none on
   d4a8699.
5. `npm run test:browser` in `js/popnei` passes.

Stands on: nothing of this plan.

Tasks:

- [x] 1.1 The core, in `crates/popnei/src/filters.rs` or a module
  `filters/first_n.rs` beside `filters/regions.rs`: `FirstNReader`,
  `PassStep::FirstN`, its place in `chain_of` and in
  `refuse_a_second_filter_of_a_kind`, `refuse_a_step_after_the_first_n`,
  `stopped_early` and the new cases of `crates/popnei/src/error.rs`, from
  "The filter that keeps the first n variants" and "The Rust interface".
  Deliverable 1.
- [x] 1.2 The Python side. In `crates/popnei-python`, the step
  `filter_first_n` in `steps.rs`, with `num_vars` read the way
  `count_of_at_least` of `source.rs` reads a count of variants, widened to
  `u64`, and not with `distance_of`, whose message speaks of base pairs;
  the refusal of a step after it in every method that adds a step; and one
  function of the crate that builds the counts of a pass, with
  `stopped_early` from the core and the steps of the pass. Every consumer of
  the crate and `Blocks::pass_stats` call that function, and the copies of
  `filtering_of` that each module has now go. In
  `python/popnei`, `Variants.filter_first_n`, `PassStats.stopped_early`
  and the kinds of `filters.py`. Deliverables 2 and 3. Needs 1.1.
- [x] 1.3 The TypeScript side, the same in `crates/popnei-js` and
  `js/popnei`: `filterFirstN`, `stoppedEarly` in `PassStats` and in
  `passStatsOf`, the counts built in one function of the crate, and the
  kinds of `filters.ts` and of the lists of the tests. Deliverables 2, 4
  and 5. Needs 1.1; it can run beside 1.2, whose files it does not touch.

What could go wrong: the filter is the first reader that stops asking its
source before the source ends, so a reader that does work when it is
dropped, or the thread of the reader one block ahead, could wait or read
on; the tests of deliverable 1 over the source that never ends are the ones
that find it, and 1.1 runs them before anything is built on it. When one
hangs, the subagent runs it under a timeout and finds which reader keeps
reading. A fix to a reader that exists, the reader one block ahead of
`crates/popnei/src/block.rs` among them, is part of 1.1, in a commit of its
own before the filter, with a test that fails without it. When the fix
would change what a reader does in a pass that is not ended early, the
orchestrator stops and asks the owner. A consumer
whose counts do not go through the one function of 1.2 or 1.3 gives no
`stopped_early`; the tests that the field is in the counts of every
consumer find it.

## 2. The filter that keeps variants at random

What it gives: `variants.filter_randomly(0.1)` and
`variants.filterRandomly(0.1)` keep about one variant in ten, the same ten
in every pass and every calculation on that `Variants`, and another sample
with another `seed`.

Deliverables:

1. Cargo tests whose names contain `random_filter`, in `filters`, with the checks
   of "How it is verified" of the random filter: the five draws of Java
   from 1234567; the worked example at `RandomFilter::filter_block`, whole
   and as blocks of 3 and 7; the three rows of the table of `many.vcf` in
   blocks of 7 and of the default size, with the counts of the first; the
   random filter at 0.5 and seed 42 and then the first 10, which gives the
   ten positions of the spec; the refusals of a keep rate of NaN, -0.1 and
   1.5 and of a second filter of this kind. `cargo test -p popnei --lib
   random_filter -- --list` counts them, and counted 0 on d4a8699; plain
   `random` matches a test of the PCoA.
2. Pytest tests whose names contain `filter_randomly`, with
   the checks of its "How it is verified": the two `iter_blocks`, the
   distances and the PCA that see the 45; `do_pca_from_variants` with 10
   components and `calc_gwas` with a kinship and
   `use_grammar_gamma_approx=True`, each equal to the same call over the 45
   variants written to a vars file with no step; the refusals; `steps`.
   `uv run pytest -k filter_randomly` exits with 5 on d4a8699.
3. Node tests whose names contain `filterRandomly`, with the checks the
   spec gives the TypeScript test. None run on d4a8699.
4. Every existing test passes, no fewer than at the end of work package 1,
   and `npm run test:browser` passes.

Stands on: work package 1, for `PassStep` and the one function of each
binding crate that builds the counts, and for the test of the random filter
followed by the first 10.

Tasks:

- [x] 2.1 The core: the generator, private to the filter, `RandomFilter`,
  `RandomlyFilteredReader`, `PassStep::Random`, its place in `chain_of`,
  in `refuse_a_second_filter_of_a_kind` and among the steps that
  `refuse_a_step_after_the_first_n` refuses, `DEFAULT_RANDOM_FILTER_SEED`
  and the case of a keep rate out of range, from "The filter that keeps
  variants at random" and "The Rust interface". Deliverable 1. The
  generator is a wrong number and not a crash when it is wrong, so it is
  the first commit of the task, with the test of Java's five draws.
- [x] 2.2 The Python side: `filter_randomly` in `steps.rs` of
  `crates/popnei-python`, with the seed read as a whole number of 64 bits,
  the default from the core's constant through `_core`;
  `Variants.filter_randomly` and the kinds of `filters.py`. Deliverable 2.
  Needs 2.1.
- [x] 2.3 The TypeScript side: `filterRandomly(keepRate, {seed})` in
  `crates/popnei-js` and `js/popnei`, with the seed a whole number up to
  2^53 - 1 and an options object that refuses a key it does not know, as
  every options object of the package does. Deliverables 3 and 4. Needs
  2.1; it can run beside 2.2.

What could go wrong: the draws are made in order on one thread, so a
version that draws on the pool of rayon, or once for each block and not
for each variant, gives other variants with another block size; the rows
of the table run in two block sizes for that reason.

## How the whole plan is checked at the end

After work package 2 and its review, the review of the whole branch, from
eae29a2 to its last commit, in the seven categories of the `code-review`
skill, every finding weighed and the ones that hold fixed, each in its own
commit, and then the checks of the `coding` skill, all of them. The work
report, `docs/reports/filters-random-first-n.md`, lists each finding,
whether it was fixed, and why when it was not.
