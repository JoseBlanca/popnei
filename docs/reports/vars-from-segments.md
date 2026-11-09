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
- 1.3 and 1.4, 421cfa1, 40a789c and 8478cf9: the order check of core
  `write_vars` with 5 tests, the refusal of `write_vcf` with 1. Both
  binding crates pass `allow_unsorted` true until task 2.1, because the
  fixture `long_vcf` of the pytest tests is not sorted and would fail
  otherwise. 421cfa1 makes `make_of_1_1.py` patch a file of 1.3, which
  writes `of_1_1.vars` byte for byte as it is committed, and changes the
  byte counts of two node tests of `progress.test.ts`, 42744 to 42776 and
  12243560 to 12243592, since a file of 1.3 is longer by the two entries
  of its key; the plan allowed only the version to change in existing
  tests, and these counts are the same change.

Deliverables of work package 1, run on 8478cf9:

1. `cargo test -p popnei --lib extra_field -- --list`: 35 tests (3 of them
   the gzip tests that were there before).
2. The sweep over the file with extra fields passes, in the 1563 below.
3. `cargo test -p popnei --lib not_sorted -- --list`: 6 tests, 5 new.
   `cargo test -p popnei --lib -- extra_field not_sorted`: 41 passed.
4. The `write_vcf` test is among the 35 of deliverable 1.
5. `cargo test --workspace` 1563 passed, 2 ignored; pytest 762; `npm test`
   701; `npm run test:browser` 9.
6. fmt, clippy, `--no-default-features` 1413 passed, both wasm checks.

### The review of work package 1

Six reviewers, spec, tests, errors, architecture, api and numbers, read
aa4715b..8478cf9. The numbers reviewer made seven one line changes that
break the arithmetic and the tests caught each. What held, all fixed in
900b39f..54e8966, with the spec first where it changed:

- After `select_extra_fields` the header of the vars reader still listed
  every field of the file while its blocks carried the selected ones, so a
  writer given that reader refused every block with an error that blamed
  popnei; found by two reviewers. The header now lists the selected fields.
- The error of variants out of order counted the variants written, after
  the filters, while the Python message names the source file. It now says
  "the 2nd variant written" and says how to write them anyway, with the
  name of the option of each package.
- Two writer errors named the source file instead of the file written.
- No pytest reached the new errors of the file from Python; five now do,
  on files made with pyarrow.
- Tests added: the comparison of the names of the fields when `reblock`
  joins, a read of the fields with only some of the columns, floats of the
  fields bit for bit with NaN and -0.0.
- Smaller: names stored with backquotes in one error, the five types
  written out twice by hand, a tuple of two numbers of one type, a
  comment that gave lz4 90 times where it is 1.7.
- The spec said `reblock` does not check the fields of the blocks it
  joins; it does, and the spec says so now.

Not taken: the free text of `VarsExtraField.problem` for five cases, which
follows the pattern of the crate.

After the fixes: `cargo test --workspace` 1568 passed, 2 ignored; pytest
769; `npm test` 703; `npm run test:browser` 9; fmt, clippy, both wasm
checks clean.

## Work package 2, while it runs

- 2.1, 33c9bbe: `allow_unsorted` in Python and `allowUnsorted` in
  TypeScript, with 2 pytest and 2 node tests.
- 2.2, e70c5b0 (spec) and 0de1ceb: `write_vars_from_segments`, with
  `SegmentReader` in the Python binding crate and 54 pytest tests; pytest
  823 passed. The spec contradicted itself: its tests refused a `float64`
  array for the `float32` field `GP`, which its rule of the types accepts;
  the refused case is now a `float64` array for the `int32` field `NODE_A`.
  Two decisions of the task: an empty list is taken as `pos` or a field
  whatever type numpy gives it, and a list of integers too large for
  numpy's `int64` is read value by value, so that `2**64 - 1` reaches a
  `uint64` field as an integer and not as a rounded float.
- 2.3, 2ed485a (spec) and fbd67ef: `var_fields`, `gt_fields`, the extra
  names of `iter_blocks` and `Block.extra_fields`, with 10 pytest tests
  and a cargo test. The spec did not say where a name of `fields` is told
  apart as one of the five or an extra field; the core does it in
  `needs_of_the_fields_of_a_source`, which the spec of the block now has.

Deliverables of work package 2, run on fbd67ef:

1. `uv run pytest tests/test_vars_from_segments.py`: 64 passed.
2. `uv run pytest -k unsorted`: 4 passed.
3. The node tests of `allowUnsorted` in `test/vars.test.ts`: 2, passing
   in `npm test`.
4. and 5. `cargo test --workspace` 1569 passed; `--no-default-features`
   1419; pytest 833; `npm test` 703; fmt, clippy, ruff and both wasm
   checks clean.
6. The docstrings and the doc comment of `writeVars` are in the commits
   of the three tasks.
