# Plan: the result so far of three statistics, and the three in one pass

7 October 2026. State: **draft**. It builds the two items that
`docs/specs/js_sources.md` gained on 7 October 2026 from issue 10 of the
repository, "The result so far" and "The three statistics of a file in one
pass", and what they need of the core, the part at the end of "The Rust
interface" of `docs/specs/stats.md` that starts "The result so far and the
three in one pass". The result so far is what `calcPerVarDistribs`,
`calcPerIndividualStats` and `calcVarDensity` would return over the
variants their pass has read up to a block, given to a function of the
page while the pass runs; the one pass is `calcVariantsSummary`, which
gives those three results from one reading of the file. The specs have no open point: the owner decided on 7
October 2026 that the result so far is a function the consumer calls, in
TypeScript alone, on all four calls, and that the one pass is a consumer of
exactly the three statistics. Two work packages: the core, then the wasm
binding crate and the TypeScript package. Every task stands on the one
before it, so none runs beside another.

## In and out

Built: the three calculations split into what they add up and what they
give, with a function called after each block; `calc_variants_summary` in
the core; `onSoFar` and `soFarEvery` on `calcPerVarDistribs`,
`calcPerIndividualStats`, `calcVarDensity` and `calcVariantsSummary`; and
the measurement of the time `calcVariantsSummary` saves.

Not built: anything in Python, as the owner decided; the result so far of
any other consumer.

## What has to be in place

- The branch `spec/stats-so-far`, in the worktree
  `.claude/worktrees/stats-so-far`, from `main` at d22e7a2, with the specs
  at e6779fc or later, and its start message on the board,
  `.claude/board/2026-10-07T0800-spec-stats-so-far.md`. The plan runs on
  that branch.
- The shared files it changes, which other branches may also change:
  `crates/popnei/src/stats.rs` and `stats/density.rs`;
  `crates/popnei-js/src/stats.rs`, `vars.rs`, `vcf.rs`, `source.rs`; the
  TypeScript files of the three calculations and `js/popnei/src/variant.ts`;
  `js/popnei/test/consumers.ts` and `num_passes.test.ts`.
- `/Users/jose/devel/popnei-bench/big.vcf` and `big.vars`, the files of
  "Speed" of `docs/specs/stats.md`.
- The checks of the `coding` skill on the code of `main` at d22e7a2, run on
  6 October 2026 on 8db578b, the last commit of code of the branch that
  merge brought, whose tree is the same: `cargo test --workspace` 1489
  passed and 2 ignored, summed over its `test result` lines; `cargo test -p
  popnei --no-default-features` 1339 passed; `uv run pytest` 751 passed;
  `npm test` in `js/popnei` 557 tests, 556 pass and 1 fails, `a kinship
  that does not tell the two variances apart gives none of them`, which
  fails on `main` and is not of this plan; `npm run test:browser` 9 passed.
  "No fewer" below counts from these. A fresh worktree needs `npm ci` and
  `npm run build` in `js/popnei` before `npm test`.

## 1. The core

What it gives: the work package that uses it, 2, calls the three
calculations with a function after each block and the three in one pass.
It has no Python or TypeScript side of its own; Python is unchanged, and
its tests passing untouched are the check that the three calculations give
what they gave.

Deliverables:

1. Cargo tests whose names contain `so_far`, for each of the three: with a
   function after each block over `many.vcf` read in blocks of 100, five
   calls, the `num_vars` of each 100 more than the last, the result of each
   equal to the calculation over the first that many variants, the counts of
   its filters those of the chain after that block, and an error of the
   function ending the pass with that error. `cargo test -p popnei --lib
   so_far -- --list` counts them, and counts 0 on d22e7a2.
2. Cargo tests whose names contain `variants_summary`, with the checks of
   "How it is verified" of the one pass that are made in cargo: the three
   equal to the bit to the three calculations on `many.vcf`, any one left
   out, none asked for an error, and the density alone asking its reader
   for the chromosome and the position alone. `cargo test -p popnei --lib
   variants_summary -- --list` counts 0 on d22e7a2.
3. Every existing test passes, no fewer than the baseline, with no test of
   before changed: the three calculations give what they gave.
4. The time of the three against the one, measured on `big.vcf` and
   `big.vars`, one thread, a release build, the median of 7 interleaved
   runs, written into "The three statistics of a file in one pass" of
   `docs/specs/js_sources.md` in place of the estimate, with the script that
   measured it in `crates/popnei/benches/`. The owner decided to build the
   one pass whatever it saves; when it saves less than a fifth, the
   orchestrator says so in the report and goes on.

Stands on: nothing of this plan.

Tasks:

- [ ] 1.1 The three calculations split into what they add up and what they
  give, with `SoFar`, `AfterABlock` and the three `_with` functions, in
  `crates/popnei/src/stats.rs` and `stats/density.rs`, from the part of
  "The Rust interface" of `docs/specs/stats.md` named in the opening. The
  three functions that exist keep their names and results. Deliverables 1,
  3. A result so far that is not the result over those variants would be
  wrong in silence; the comparison with the calculation over the first that
  many variants is its guard, and it goes in the commit of the split.
- [ ] 1.2 `calc_variants_summary`, `VariantsSummaryConfig` and
  `VariantsSummary`, from the same part. Deliverables 2, 3. Needs 1.1.
- [ ] 1.3 The measurement of deliverable 4, and the spec's sentence of the
  saving replaced with its numbers, in a commit of its own before the code
  of work package 2. Needs 1.2.

What could go wrong: `calc_per_var_distribs` reads one block ahead on a
thread of its own, and the function after a block runs in the loop that
takes the blocks from that thread; the counts of the filters are those the
handle of that thread gives, which the review of the spec found are the
counts as they were when the block in hand was given. Deliverable 1 checks
them natively.

## 2. The wasm binding crate and the TypeScript package

What it gives: `calcPerVarDistribs`, `calcPerIndividualStats` and
`calcVarDensity` take `onSoFar` and `soFarEvery` and call the function with
the result so far while their pass runs; `calcVariantsSummary` gives the
three in one pass, with the same two options; `numPassesOf` knows it.

Deliverables:

1. Node tests whose names contain `onSoFar`, with the checks of "How it is
   verified" of "The result so far", the stops of `test/stop.test.ts` made
   again for `onSoFar` among them. `node --test --test-reporter=spec
   --test-name-pattern=onSoFar test/*.test.ts | grep -c onSoFar`, in
   `js/popnei` after `npm run build`, prints 0 on d22e7a2; node counts each
   test file as a test, so its own count is no check here.
2. Node tests whose names contain `calcVariantsSummary`, with the checks of
   "How it is verified" of the one pass; the same command with that name
   prints 0 on d22e7a2.
3. `numPassesOf("calcVariantsSummary")` is 1, in `test/num_passes.test.ts`,
   and `calcVariantsSummary` is among the consumers of `test/consumers.ts`,
   so that the tests that run every consumer, the progress and the stops,
   run it too.
4. Every existing test passes, no fewer than after work package 1, and `npm
   run test:browser` passes. The tests of before that change are the lists
   of the consumers and their counts.

Stands on: work package 1.

Tasks:

- [ ] 2.1 `onSoFar` and `soFarEvery` on the three calculations: in
  `crates/popnei-js`, the clock, the call between blocks with the result as
  the package builds it, and the end of a pass whose function threw, by the
  path a throw of `onProgress` takes; in `js/popnei`, the two options, their
  checks and their doc comments. From "The result so far" of
  `docs/specs/js_sources.md`. Deliverables 1, 4.
- [ ] 2.2 `calcVariantsSummary` in both, with `onSoFar`, `numPassesOf`, the
  list of consumers and its doc comment, from "The three statistics of a
  file in one pass". Deliverables 2, 3, 4. Needs 2.1, whose files it
  shares.

What could go wrong: a throw from `onSoFar` comes from the consumer and not
from the source, so the binding has to end the pass as one that `onProgress`
stopped, with the source telling the page nothing more; the stops of
deliverable 1 are the check. The one pass is the first consumer of the
package whose reader is asked for the genotypes for one result and not for
another; deliverable 2's comparison with the three consumers is the check.
