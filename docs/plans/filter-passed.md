# Plan: the filter of the variants that passed their FILTER

6 October 2026. State: **under way**, approved by the owner on 6 October 2026. It builds what the specs gained on 6
October 2026 from issue 9 of the repository: the filter `filter_passed` /
`filterPassed`, of the kind `"passed"`, that keeps the variants whose
FILTER column was `PASS` or a dot, and the `passed` column it reads. The
parts of the specs are "The filter of the variants that passed their
FILTER" and what it adds to "The Rust interface" of
`docs/specs/filters.md`; the `passed` field of `Block` in
`docs/specs/block.md`; `PASSED` and `ALL` of `Needs` in
`docs/specs/variant.md`; the paragraph "When `PASSED` is asked for", the
writer's FILTER for a variant that failed in `docs/specs/io_vcf.md`; and the `passed`
column and format 1.2 in `docs/specs/io_vars.md`. Two work packages: the
column, then the filter.

## In and out

Built: the column in the core, filled by the VCF reader, stored and read
by the vars file, written back by the VCF writer; the filter in the core,
both binding crates and both packages; the docstrings of `open_vcf` and
`openVcf` that say when to use `only_passed` and when the filter.

The specs have no open point: the owner decided on 6 October 2026 that
the VCF writer writes `FAIL` for a variant that failed.

Not built: the column among the fields that `iter_blocks` names
(`docs/specs/block.md`).

## What has to be in place

- The branch `spec/filter-passed`, in the worktree
  `.claude/worktrees/filter-passed`, from `main` at 6abacd1, with the specs
  at 496d1f5 or later, and its start message on the board,
  `.claude/board/2026-10-06T1700-spec-filter-passed.md`. The plan runs on
  that branch.
- The shared files it changes, which other branches may also change:
  `crates/popnei/src/block.rs`, `variant.rs`, `error.rs`, `filters.rs`,
  `io/vcf.rs`, `io/vars.rs`, `io/vcf/writer.rs`; `steps.rs` and
  `errors.rs` of both binding crates; `python/popnei/variant.py`,
  `io_vcf.py`, `filters.py`; `js/popnei/src/variant.ts`, `io_vcf.ts`,
  `filters.ts`; `docs/architecture.md`.
- bcftools 1.24 at `/opt/homebrew/bin/bcftools`, which gave the numbers of
  the spec on 6 October 2026.
- The checks of the `coding` skill on the code of `main` at 6abacd1, run
  on 6 October 2026 on the commit before that merge, 5987cc5, which has the
  same tree: `cargo test --workspace` 1453
  passed and 2 ignored, summed over its `test result` lines; `cargo test -p
  popnei --no-default-features` 1303 passed; `uv run pytest` 741 passed;
  `npm test` in `js/popnei` 551 tests, 550 pass and 1 fails, `a kinship that
  does not tell the two variances apart gives none of them`, which fails on
  `main` and is not of this plan; `npm run test:browser` 9 passed. "No
  fewer" below counts from these. A fresh worktree needs `npm ci` and `npm
  run build` in `js/popnei` before `npm test`.

## 1. The column of whether each variant passed

What it gives: a vars file written from a VCF keeps, for each variant,
whether its FILTER was `PASS` or a dot, in a column `passed` that pyarrow
shows, and a VCF written back from it marks the ones that failed with
`FAIL`. This work package stops before the filter, which is work package 2;
what a user sees of it is the column of the vars file.

Deliverables:

1. Cargo tests whose names contain `passed_column`, with the checks of the
   specs on the column: the 25 false of `many.vcf` with `only_passed` false
   and none with it true, at the VCF reader; `retain_vars`, `reblock` and
   `check` with the column; a vars file of 1.2 written and read back with
   it, one written from a source without it that has no such column, and a
   reader of 1.1 that ignores it; the VCF writer's 25 `FAIL` lines and its
   `##FILTER` line, and no line for a vars file without the column.
   `cargo test -p popnei --lib passed_column -- --list` counts them, and
   counts 0 on 6abacd1.
2. Every existing test passes, no fewer than the baseline. A test of
   before that changes is one of those the spec names in `docs/specs/variant.md`
   (a set compared with `ALL`, a version of 1.1, six columns), the pytest
   test of `write_vars` whose schema gains `passed`, and the expected header
   of the vars file written to a VCF, which gains the `##FILTER` line. Any
   other test that fails is reported to the owner and not changed.
3. The pytest test of `write_vars` on `many.vcf` asserts the seven columns,
   `format_version` `"1.2"` and the 25 false of `passed`, and, written from
   `many.vcf` read with the default, none; it fails on 6abacd1, where the
   file has six columns.

Stands on: nothing of this plan.

Tasks:

- [x] 1.1 The column in the core: `Needs::PASSED` in `ALL`, `Block.passed`
  and every place that destructures a `Block`, `fields`, `check`,
  `retain_vars` and `reblock`, and the VCF reader filling it from FILTER
  with the function of `only_passed`. The binding crates are made to
  compile with the new field and change nothing else. From the `passed`
  field of `docs/specs/block.md`, `PASSED` of `docs/specs/variant.md` and
  "When `PASSED` is asked for" of `docs/specs/io_vcf.md`. Deliverables 1, 2.
- [ ] 1.2 The vars file of format 1.2, in `crates/popnei/src/io/vars.rs`:
  the column written when the source has it and read when the file has it,
  from "What it holds" of `docs/specs/io_vars.md`. The pytest test of
  `write_vars`. Deliverables 1, 2, 3. Needs 1.1.
- [ ] 1.3 The VCF writer of a vars file, in
  `crates/popnei/src/io/vcf/writer.rs`: `FAIL` and the `##FILTER` line, as
  the owner decided, with its check run with bcftools as the spec gives
  it. Deliverables 1, 2. Needs 1.2, for a vars reader that says
  whether the file has the column.

What could go wrong: `Needs::ALL` is the default of every reader, so every
reader of a VCF now fills the column. The cost is a comparison and a
`bool` per variant, which the review of the spec judged not measurable
beside the parse of the genotypes, 92% of a read, without timing it. 1.1
times one read of the plain 403 MB `big.vcf`, 100000 variants of 1000
individuals in `/Users/jose/devel/popnei-bench/`, with the genotypes alone
and with `ALL`, on one thread, before and after, the median of 7 runs; the
median was 0.597 s for the genotypes alone on 6 October 2026, after the
fixes of the plan of the ploidy. When either
is slower by more than 3%, the orchestrator stops and asks the owner
whether to keep the column in `ALL` or to have the VCF reader fill it only
when a filter asks for it, with the two timings.

## 2. The filter

What it gives: `variants.filter_passed()` and `variants.filterPassed()`
take out the variants that failed their FILTER as a step of every pass,
and every `pass_stats` has their count under `"passed"`. A vars file
without the column is refused at the first block, with its path.

Deliverables:

1. Cargo tests in a module `filters::passed`, with the checks of "How it is
   verified" of the filter: the 475 positions and the counts, in blocks of
   7 and of the default size; the MAF filter after it, 475 and 364; the
   variants against those of `only_passed` true; the same over a vars file;
   the refusal of a source without the column, at the first block; the
   refusal of a second filter of this kind and of this one after the
   filter of the first n. `cargo test -p popnei --lib filters::passed --
   --list` counts them, and counts 0 on 6abacd1.
2. Every existing test passes, no fewer than after work package 1. The
   tests of before that change are the lists of the kinds of a step, in
   `tests/test_filter_first_n.py` and `js/popnei/test/filter_first_n.test.ts`.
3. Pytest tests whose names contain `filter_passed`, with the checks the
   spec gives Python, among them the `ValueError` that starts with the path
   of a vars file rewritten by pyarrow without the column. `uv run pytest -k
   filter_passed` exits with 5 on 6abacd1.
4. Node tests whose names contain `filterPassed`, with the checks the spec
   gives TypeScript. `node --test --test-reporter=spec
   --test-name-pattern=filterPassed test/*.test.ts | grep -c filterPassed`,
   in `js/popnei` after `npm run build`, prints 0 on 6abacd1; node counts
   each test file as a test, so its count of tests is no check here.
5. `npm run test:browser` passes.
6. The docstrings of `open_vcf` and `openVcf` say when to use `only_passed`
   and when the filter, and no longer that nothing says which variants
   failed; those of `filter_passed` and `filterPassed` say to add it first,
   after the filter by regions.

Stands on: work package 1.

Tasks:

- [ ] 2.1 The core: `PassStep::Passed`, `PassedReader`, its place in
  `chain_of`, `refuse_a_second_filter_of_a_kind` and
  `refuse_a_step_after_the_first_n`, and its two error cases, the one of a
  source without the column among those that name the file, from "The
  filter of the variants that passed their FILTER" and "The Rust
  interface" of `docs/specs/filters.md`. Deliverables 1, 2.
- [ ] 2.2 The Python side: the step in `crates/popnei-python/src/steps.rs`,
  `Variants.filter_passed`, the kinds of `filters.py`, the docstring of
  `open_vcf`, and the list of kinds in `tests/test_filter_first_n.py`.
  Deliverables 2, 3, 6. Needs 2.1.
- [ ] 2.3 The TypeScript side: the step in `crates/popnei-js/src/steps.rs`,
  `filterPassed`, the kinds of `filters.ts` and of `source.rs`, the doc
  comment of `openVcf`, and the list of kinds in
  `js/popnei/test/filter_first_n.test.ts`. Deliverables 2, 4, 5, 6. Needs
  2.1; it can run beside 2.2, whose files it does not touch.

What could go wrong: the error of a source without the column has to reach
Python with the path of the vars file, through the errors that name the
file; the pytest test of deliverable 3 is the one that finds it when it
does not.
