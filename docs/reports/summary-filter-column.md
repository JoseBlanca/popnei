# Report: the counts of the FILTER column in the summary of a file

7 October 2026. The work report of `docs/plans/summary-filter-column.md`,
on the branch `spec/summary-filter-column`. State: under way.

## Work package 1, while it goes

- 1.1, f2ce5ac: `keeps_passed` in both binding crates and both packages.
  3 pytest tests and 3 node tests; `uv run pytest -k keeps_passed` gives
  3 passed. The trait of the sources in `source.rs` was left unchanged,
  since the packages read the property from each source class.
- 1.2, 58cc067: the counts in the core. `cargo test -p popnei --lib
  stats::summary::tests::filter_column -- --list` counts 9 tests;
  `cargo test --workspace` 1515 passed and 2 ignored. The test of the
  order of the errors before the pass passed against the stub too, as
  there was no check of the record yet to put in the wrong place; it
  fails only when that check moves ahead of the others. The binding
  crates match the error with a wildcard, so neither `errors.rs` changed.
  The JS binding destructures the new field as `filter_column: _`, which
  1.3 replaces.
- 1.1 and 1.2 ran side by side in one tree, each committing by path, and
  neither took the other's files.
