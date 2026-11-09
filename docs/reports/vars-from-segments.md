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
- 1.2, eaf4cf2 (spec) and 17cb60e (code): the vars file of 1.3, with 22
  tests; `cargo test --workspace` 1557 passed. The spec gained three
  things the task had to decide: a malformed entry of the two lists is
  the error of a file that is not a vars file, a selection of no field
  gives a block with no extra fields, and what `set_extra_fields` refuses.
  The bounds walk of the reader needed no change. Two findings outside
  the plan, for the owner: the existing byte sweep of the reader counts
  panics only on the thread of the test, while the reader decodes on the
  threads of rayon natively, so its assertion on the panics it caught
  cannot fail; the new sweep runs in a pool of one thread and caught 11
  panics, all inside arrow-rs, none out of the reader. And
  `tests/reference/vars/make_of_1_1.py` asserts that popnei writes 1.2, so
  it fails if it is run again; task 1.3 fixes it.
