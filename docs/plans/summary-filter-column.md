# Plan: the counts of the FILTER column in the summary of a file

7 October 2026. State: **draft**. It builds what the specs gained on 7
October 2026 from issue 12 of the repository: a fourth option of
`calcVariantsSummary`, `filterColumn`, that counts how many of the variants
of its pass passed their FILTER and how many failed, without taking any
out; and `keeps_passed` / `keepsPassed` of `Variants`, which says whether
the source recorded FILTER. The parts of the specs are "The three
statistics of a file in one pass" of `docs/specs/js_sources.md`; the
paragraphs on the counts of the FILTER column after the code block of
`calc_variants_summary` in "The Rust interface" of `docs/specs/stats.md`;
the paragraph "A `Variants` also has `keeps_passed`" of
`docs/specs/variant.md`; and `keeps_passed` of `SourceHeader` in
`docs/specs/block.md`. One work package of three tasks.

## In and out

Built: the counts in the core summary with their two errors; the option
and the field of the result in the wasm binding and the TypeScript
package; `keeps_passed` in the core binding crates and both packages.

The specs have no open point: the owner decided on 7 October 2026 that the
counts are an option, as each of the other three statistics is, and that
a source without the record is an error.

There is no comparison with pyNei: pyNei has no such function, and Python
has no `calcVariantsSummary`. The counts are checked exactly against the
FILTER column of `many.vcf`, counted with `grep`, as the spec gives.

Not built: `filter_passed` refusing a source without the record before its
first block, from the header, where today it refuses at the first block.
The spec does not ask for it.

## What has to be in place

- The branch `spec/summary-filter-column`, in the worktree
  `.claude/worktrees/summary-filter-column`, from `main` at 9300d73, with
  the specs at ea508e2 or later, and its start message on the board,
  `.claude/board/2026-10-07T2130-spec-summary-filter-column.md`. The plan
  runs on that branch.
- The shared files it changes, which `stats/hist-right-edges` also changes:
  `crates/popnei/src/error.rs`, `js/popnei/src/stats.ts`. Others it changes:
  `crates/popnei/src/stats/summary.rs`, the sources of both binding crates
  (`vcf.rs`, `vars.rs`, `source.rs`, `summary.rs`, `errors.rs`),
  `python/popnei/variant.py`, `js/popnei/src/variant.ts`.
- The checks of the `coding` skill on 9300d73, run on 7 October 2026 in
  this worktree: `cargo test --workspace` 1506 passed and 2 ignored,
  summed over its `test result` lines; `cargo test -p popnei
  --no-default-features` 1356 passed and 2 ignored; `uv run pytest` 751
  passed; `npm test` in `js/popnei` 650 tests, 649 pass and 1 fails, "a
  kinship that does not tell the two variances apart gives none of them",
  which fails on `main` and is not of this plan; `npm run test:browser` 9
  passed. The checks of deliverables 1 to 3 below were run there too, and
  gave 0 tests, exit 5 and 0. "No fewer" below counts from these. A fresh
  worktree needs `npm ci` and `npm run build` in `js/popnei` before `npm
  test`.

## 1. The counts of the FILTER column, and whether a source recorded it

What it gives: `calcVariantsSummary(variants, { filterColumn: {} })` gives
`filterColumn: { passed, failed }` beside whatever else was asked for, and
`variants.keeps_passed` in Python and `variants.keepsPassed` in TypeScript
say beforehand whether it can be asked.

Deliverables:

1. Cargo tests whose names contain `filter_column`, in
   `crates/popnei/src/stats/summary/tests.rs`, with the checks the spec
   gives the core: 475 and 25 on `tests/reference/vcf/many.vcf` alone and beside the three,
   which stay as they were; 475 and 0 with the filter of the variants that
   passed; the vars file of `many.vcf` the same as the VCF; that the counts
   alone ask the reader for `Needs::PASSED` and nothing else and give a
   pass of 500 variants; the result so far over the vars file in batches of
   100; `FilterColumnNotRecorded` before any block is asked for, from the
   reader built for the tests with a header of `keeps_passed` false;
   `FieldsNotInTheBlock` with `PASSED` for a block without the column; and
   the order of the errors before the pass. `cargo test -p popnei --lib
   stats::summary::tests::filter_column -- --list` counts them, and counts
   0 on 9300d73.
2. Pytest tests whose names contain `keeps_passed`: true for `many.vcf`
   and for the vars file written from it, false for
   `tests/reference/vars/of_1_1.vars`. `uv run pytest -k keeps_passed`
   exits with 5 on 9300d73.
3. Node tests whose names contain `filterColumn` or `keepsPassed`, with
   the checks of "How it is verified" in "The three statistics of a file in
   one pass" of the spec, and `keepsPassed` true, true and false on the
   three files of deliverable 2. `node
   --test --test-reporter=spec --test-name-pattern='filterColumn|keepsPassed'
   test/*.test.ts | grep -cE 'filterColumn|keepsPassed'`, in `js/popnei`
   after `npm run build`, prints 0 on 9300d73.
4. Every existing test passes, no fewer than the baseline. One existing
   test changes: the node test "calcVariantsSummary with none of the three
   is an Error that says to ask for one", which matches a message that
   names four statistics now. Any other test that fails is reported to the
   owner and not changed.
5. `npm run test:browser` passes.
6. The doc comments of `calcVariantsSummary` and of `VariantsSummary` in
   `js/popnei/src/stats.ts` and the passage of `js/popnei/README.md` on
   `calcVariantsSummary` say what `filterColumn` counts, that it is after
   every step, and that a source without the record is an error that
   `keepsPassed` tells beforehand.

Stands on: nothing of this plan.

Tasks:

- [ ] 1.1 `keeps_passed` of `Variants`: each source of
  `crates/popnei-python/src` and `crates/popnei-js/src` reads it from the
  header of its reader when it is opened, as it reads the individuals and
  the ploidy, and `python/popnei/variant.py` and `js/popnei/src/variant.ts`
  give it as a property. From the paragraph "A `Variants` also has
  `keeps_passed`" of `docs/specs/variant.md`. Deliverables 2, 3 (the
  `keepsPassed` tests), 4. It can run beside 1.2, whose files it does not
  touch.
- [ ] 1.2 The counts in the core: `filter_column` of
  `VariantsSummaryConfig`, `FilterColumnCounts` and the field of
  `VariantsSummary`, the counting in `crates/popnei/src/stats/summary.rs`
  with `Needs::PASSED`, the number of variants of a pass of the counts
  alone, the check of the header before the pass in the order of the spec,
  `Error::FilterColumnNotRecorded` and the message of
  `VariantsSummaryOfNoStatistic` in `error.rs`, with the new case in the
  category of `PassedNotRecorded`. The reader built for the tests gets a
  header of its own and records the blocks asked for. From the paragraphs
  on the counts of the FILTER column of "The Rust interface" of
  `docs/specs/stats.md`. Deliverables 1, 4. A wrong count is silent, and
  deliverable 1 guards it.
- [ ] 1.3 The option in TypeScript: `filterColumn` in
  `crates/popnei-js/src/summary.rs` and its `num_vars_of`, which today
  gives 0 when none of the three is there, the result so far and the final
  result; the check and the message of a call of none in
  `js/popnei/src/stats.ts`, its types and doc comments; the README; and
  the node tests in `js/popnei/test/variants_summary.test.ts`. From "The
  three statistics of a file in one pass" of `docs/specs/js_sources.md`.
  Deliverables 3, 4, 5, 6. Needs 1.1 and 1.2.

What could go wrong: with the counts alone, the pass counts its variants
in two places, `AddedUp::num_vars` in the core and `num_vars_of` in the
binding crate, both of which give 0 today when none of the three is
there; a pass of the counts alone then fails with "the pass gave no
variant". The test of the 500 variants in deliverables 1 and 3 is the one
that finds a place left unchanged.
