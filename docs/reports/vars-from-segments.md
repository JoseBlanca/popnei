# Report: a vars file written from segments, with extra fields

9 November 2026. The work report of `docs/plans/vars-from-segments.md`,
under way on branch `spec/vars-from-segments`.

## Before the first task

The checks of the `coding` skill on 95c5874, run on 9 November 2026 in
the worktree, which "no fewer" in the plan counts from:

- `cargo test --workspace`: 1526 passed, 2 ignored, summed over its
  `test result` lines.
- `cargo test -p popnei --no-default-features`: 1376 passed, 2 ignored.
- `uv run maturin develop && uv run pytest`: 762 passed.
- `npm ci && npm run build && npm test` in `js/popnei`: 701 tests, 701
  passed.

The checks of the deliverables found no test before the work: no test
named `extra_field` or `not_sorted` exists in the crate, no pytest test
`unsorted` and no file `tests/test_vars_from_segments.py`, and no node
test `allowUnsorted`.

Changes from the `following-plans` skill: the plan runs on the branch of
its spec, `spec/vars-from-segments`, as the earlier plans of an issue
did, and not on a `plan/` branch. The owner approved the plan by asking
for it to be written and carried out in one message, so it was not shown
to them, nor sent to the `first-reader`, before the work.

## Work package 1, while it runs

- 1.1, 3dd24ec: the extra fields of a block, with 9 tests; `cargo test
  --workspace` 1535 passed, 2 ignored. It added two cases to the error
  of the crate that the spec did not name, an extra field of another size
  and one of another type than it says, which `check` finds; the spec of
  the block now names them. The check of deliverable 1 counted 3 tests on
  95c5874 and not 0, three tests of gzip members whose names hold
  `extra_field`; the plan says so now.
