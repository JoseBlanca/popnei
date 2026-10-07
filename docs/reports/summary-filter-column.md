# Report: the counts of the FILTER column in the summary of a file

7 October 2026. The work report of `docs/plans/summary-filter-column.md`,
on the branch `spec/summary-filter-column`. State: done.

The plan is done, on the branch `spec/summary-filter-column`, not merged.
`calcVariantsSummary(variants, { filterColumn: {} })` now gives, beside
whatever else it was asked for, how many of the variants of its pass passed
their FILTER (`PASS` or a dot) and how many failed, without taking any out:
on `many.vcf` opened with `onlyPassed: false`, `{ passed: 475, failed: 25 }`
over 500 variants, and `{ passed: 475, failed: 0 }` after `filterPassed`,
since the counts are of the variants that the filter steps of the pass
let through. A source that did not record FILTER is an `Error` before the
first block is read: a vars file of format 1.0 or 1.1, and a vars file
that holds no variant, which popnei writes without the column. That is
as you decided, and `variants.keepsPassed` in TypeScript,
`variants.keeps_passed` in Python, says so beforehand. Every check of the `coding` skill passes, and each runs
more tests than before the plan; one node test fails, "a kinship that does
not tell the two variances apart gives none of them", which fails in the
same way on `main`. No decision is left open. What is asked of you is the
order to merge the branch into `main`, after which issue 12 can be closed.

Since the plan started, `main` gained the branch of issue 11, the edges of
the histograms, at 6943153; this branch merges into it without a textual
conflict, and I run the checks again on `main` after the merge. The
branch `issue-13-write-in-pieces`, under way, changes
`crates/popnei-js/src/vcf.rs` and `vars.rs` too: it changes `write_vars`
and `write_vcf`, and this branch adds a field that holds `keeps_passed` to
each source and an argument to `calc_variants_summary`. Neither reads what
the other changes, so whichever merges second may meet a textual conflict
in those two files and no change of behaviour.

## Work package 1: the counts and `keeps_passed`

It finished as planned, in three tasks and one round of fixes after the
review. The deliverables, each checked on d890d36:

1. `cargo test -p popnei --lib stats::summary::tests::filter_column --
   --list` counts 13 tests, 0 before the plan; all pass.
2. `uv run pytest -k keeps_passed`: 4 passed, exit 5 before the plan.
3. The node check of the plan prints 15, 0 before the plan: 11 tests of
   `filterColumn` and 4 of `keepsPassed`.
4. `cargo test --workspace` 1519 passed and 2 ignored, 1506 before;
   `cargo test -p popnei --no-default-features` 1369 passed, 1356 before;
   `uv run pytest` 755 passed, 751 before; `npm test` 665 tests, 664 pass,
   650 and 649 before, and the one failing test is the one of `main`.
   `cargo fmt`, `cargo clippy -D warnings`, `cargo wasm-check`, `cargo
   wasm-check-js` and `ruff` are clean. The one test of before that
   changed is the node test of a call that asks for nothing, whose message
   names four statistics now.
5. `npm run test:browser`: 9 passed.
6. The doc comments of `calcVariantsSummary` and `VariantsSummary` and the
   README say what the counts are, that they are after every filter step, and
   that `keepsPassed` tells beforehand whether they can be asked for; the
   `api` reviewer read them against the spec.

Changes to the plan: none. Two additions of the implementers inside their
tasks: `calc_variants_summary` of the JS binding took an eighth argument
and with it `#[expect(clippy::too_many_arguments)]`, as the other
functions of those files have; and two checks of the option that the spec
implies, a `filterColumn` with a key and one that is not an object, are
refused and tested.

The review sent five reviewers, `spec`, `tests`, `errors`, `api` and
`binding`. What it found, all of it fixed in 0e88fd8 and d890d36:

- The paragraph on `filterColumn` had been pasted into the doc comments of
  `calcPerVarDistribs`, `calcPerIndividualStats` and `calcVarDensity`,
  which refuse the option. Found by `spec` and `api`.
- Two guards of the counts had no test: the check that the column has one
  entry for each variant, and the refusal of a block of no variants. With
  both removed every test passed, and a column cut to half its length gave
  counts of 235 and 15 with no error. Found by `spec` and `tests`. The
  spec now names the counts beside the three statistics for a block of no
  variants (966a7aa).
- No test made a pass of the counts alone that keeps no variant, nor
  checked the order of the errors inside a block for the counts.
- The node test of the order of the errors before the pass used a window
  of 0 base pairs, which TypeScript refuses before the core is called, so
  its density half checked nothing in the core. It now uses a density the
  core refuses, 20000000 windows of 1 base pair.
- The error of a source without the record said a vars file holds it when
  "written from a VCF"; it says "written from a source that had it" now,
  as the spec does.
- `keepsPassed` had no test after `free()`, nor on a vars file of no
  variants, which the spec says is false; both tests were added and passed
  on the code as it was.

Not taken: the `errors` reviewer suggested that the message of the core
error name `keepsPassed`; the messages of the core are shared by Python and
TypeScript and name neither package's functions, and the doc comment of
`calcVariantsSummary` names it.

What you should know:

- Seen outside the plan, and not changed: `cargo clippy -p popnei-js
  --target wasm32-unknown-unknown` fails on `main` too, at
  `crates/popnei/src/variant.rs:1079`, a `#[expect]` of a lint that is not
  raised on wasm32. The checks of the `coding` skill do not run that
  command, so nothing has caught it.
- Also seen outside the plan: in Python, a block that lacks a field its
  source promised, which is a defect of popnei, reaches the user as a
  `ValueError` and not as the `RuntimeError` of a defect. It is older than
  this branch and holds for every consumer.
- `PassedNotRecorded`, the error of `filterPassed` on a source without the
  record, still says "written from a VCF", the wording this branch
  corrected for its own error.

## How the work went

This section is for whoever next revises a skill or writes a plan. The
owner can stop here.

- Tasks 1.1 and 1.2 ran side by side in one tree, one in the binding
  crates and one in the core, each committing by path, and neither took
  the other's files; the one-writer-per-file rule was enough.
- The defect that two reviewers found, a doc paragraph pasted into four
  doc comments, came from one search and replace of the implementer that
  matched four places. A check after an edit of that kind, that the text
  appears once, would have caught it.
- The guards of the counts were written as the spec asked and not
  tested: no fixture reached them, which "A test has to be able to fail"
  of the `coding` skill names as the first way a test cannot fail. The
  implementer's report said that 8 of its 9 tests failed against a stub,
  which told the orchestrator nothing about the code that no test reached;
  the reviewers found it by deleting each guard and running the tests. A
  report that names which refusals of the code have a test, and not only
  how many tests failed first, would have shown it before the review.
- What the tasks cost, in tokens of their subagents, the fixes included:
  1.1 94000, 1.2 134000, 1.3 131000. The five reviewers together used
  343000, about as much as the three tasks, for one work package of an
  afternoon.
