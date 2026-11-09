# Plan: a vars file written from segments, with extra fields

9 November 2026. State: **under way**, approved by the owner on 9 November 2026, who
asked the session that wrote the spec to write the plan and carry it out.
It builds what the specs gained on 9 November 2026 from issue 19: the
extra fields of a vars file, named arrays of a value for each variant or
for each genotype, carried through the filters and read back by
`iter_blocks`; the writer of segments, `write_vars_from_segments`, which
writes a vars file from dicts of arrays that a user's code yields; and
the check that the variants of a vars file are sorted, in both writers.
The parts of the specs are "The extra fields" and "The writer of
segments" of `docs/specs/io_vars.md`, with its edits to "The file", to
`write_vars` under "The writer" and to "The Rust interface"; the extra
fields of `Block`, `SourceHeader` and the Python `Block` in
`docs/specs/block.md`; `EXTRA` of `Needs` and `var_fields` / `gt_fields`
of `Variants` in `docs/specs/variant.md`; `retain_individuals` in
`docs/specs/filters.md`; and the refusal of the VCF writer in
`docs/specs/io_vcf.md`. Two work packages of seven tasks.

## In and out

Built: the extra fields in the core block, the filters, `reblock` and the
vars file of format 1.3; the order check of `write_vars`; the refusal of
`write_vcf`; the writer of segments, `allow_unsorted`, `var_fields`,
`gt_fields`, the extra names of `iter_blocks` and `Block.extra_fields` in
Python; `allowUnsorted` of `writeVars` in TypeScript.

The specs have no open point. The owner decided on 9 November 2026 that
`write_vcf` refuses a source with extra fields and that the TypeScript
functions of the extra fields and of the writer of segments wait.

There is no comparison with pyNei: pyNei has no genotype fields and
another file. The reference outside the project is pyarrow, which opens
the files, as the spec gives.

Not built: anything of "Not in this spec" of `docs/specs/io_vars.md`.

## What has to be in place

- The branch `spec/vars-from-segments`, in the worktree
  `.claude/worktrees/spec-vars-from-segments`, from `main` at 95c5874,
  with the specs at bf4b64f or later, and its start message on the board.
  The plan runs on that branch.
- The shared files it changes: `crates/popnei/src/block.rs`,
  `variant.rs`, `error.rs`, `io/vars.rs`, `io/vcf.rs`, the sources and
  errors of both binding crates, `python/popnei/block.py`, `variant.py`,
  `io_vars.py`, `__init__.py`, `js/popnei/src` for `writeVars`. Every
  literal of `Block` and of `SourceHeader` in the workspace gains a field.
- The checks of the `coding` skill on 95c5874, run on 9 November 2026 in
  this worktree, are in the report of this plan, `docs/reports/
  vars-from-segments.md`, with the counts that "no fewer" below counts
  from.

## 1. The extra fields and the order in the core

What it gives, in the words of work package 2: a block that carries extra
fields through every reader and filter, a vars file of 1.3 that holds
them, a core `write_vars` that refuses variants out of order unless it is
told otherwise, and a `write_vcf` that refuses extra fields. It stops in
the core because the only source of extra fields, the reader of
segments, is in the Python binding crate; work package 2 makes the checks
from Python.

Deliverables:

1. Cargo tests whose names contain `extra_field`, in the tests of
   `block.rs`, of the filters and of `io/vars.rs`, with the checks of "How
   it is verified" of "The extra fields": `retain_vars`, the filter of
   individuals keeping the third and the first, `reblock` of 3 and 2 into 4
   and 1, `Block::check` refusing a genotype field one value short, a
   round trip of a block with a `uint64` variant field and a `float32`
   genotype field through the writer and the reader, `select_extra_fields`
   with one name, with `None` and with a name the file lacks, and one test
   for each error of the file of "In the file". `cargo test -p popnei --lib
   extra_field -- --list` counts them, and counts 3 on 95c5874, three
   tests of gzip members in `io::vcf` whose names hold `extra_field`.
2. The byte sweep of the reader, four values for each byte, also over the
   file of the four variants of `cases.vcf` with a `float32` genotype field
   and an `int64` variant field, with no panic out of the reader and no
   abort.
3. Cargo tests whose names contain `not_sorted` at `write_vars` of the
   core: chr1 5000 then chr1 100, and chr1, chr2, chr1, each the error
   with both variants and their numbers, and each written with
   `allow_unsorted` true. The existing tests named `NOT_SORTED` pass with
   `allow_unsorted` true.
4. A cargo test whose name contains `extra_field` in the tests of
   `io/vcf.rs`: `write_vcf` of a source whose header has extra fields is
   the error that names them, before anything is written to its sink.
5. Every existing test passes, no fewer than the baseline, with the
   format version they read now `"1.3"` and the two entries of the extra
   fields `[]`. Any other test that fails is reported and not changed.
6. `cargo fmt`, `cargo clippy`, both `cargo test`, `cargo wasm-check` and
   `cargo wasm-check-js` pass.

Stands on: nothing of this plan.

Tasks:

- [x] 1.1 The extra fields of a block: the types of the code block after
  `Block` in "The Rust interface" of `docs/specs/block.md`, `extra` of
  `Block`, `extra_fields` of `SourceHeader`, `EXTRA` of `Needs` in `ALL`
  (`docs/specs/variant.md`), and what `fields`, `check`, `retain_vars`,
  `retain_individuals` and `Reblock` do with them, by "What carries
  them" of `docs/specs/io_vars.md`. Every literal of `Block` and
  `SourceHeader` in the workspace gains its field. Deliverables 1 (the
  block, filter and reblock tests), 5, 6. A filter that dropped or
  misordered the values of a genotype field would be silent, and the
  tests of deliverable 1 guard it.
- [x] 1.2 The vars file of 1.3: the extra columns and the two entries of
  the `popnei` key in the writer, `set_extra_fields`, `write_vars` taking
  the fields of the header of its reader; the reader's reading and
  refusals, `select_extra_fields`, the bounds of the new columns, and
  `extra_fields` of its header; the errors. From "In the file" and the
  edits to "The Rust interface" and its list of errors. Deliverables 1
  (the vars tests), 2, 5, 6. Needs 1.1.
- [x] 1.3 The order: the check in `write_vars` of the core and its
  `allow_unsorted`, the error, both binding crates passing `false` until
  work package 2 gives them the argument, and the existing tests of
  unsorted blocks passing `true`. From "Sorted variants" and "How it runs"
  of "The writer of segments". Deliverables 3, 5, 6. Needs 1.2, since both
  change `write_vars`.
- [x] 1.4 The refusal of `write_vcf`, from the first case of "The cases a
  reader of the rules would not guess" of the VCF writer in
  `docs/specs/io_vcf.md`. Deliverables 4, 5, 6. Needs 1.1; can run beside
  1.2 and 1.3, whose files it does not touch.

What could go wrong: `Reblock::join` takes the columns of the two blocks
it joins to be the same and does not check it, so a mistake in 1.1 there
loses a field without an error; the vars reader bounds every buffer of a
batch before arrow-rs reads it, and a column it does not walk stops the
walk, so the new columns have to be walked or the bounds of every column
after them fall to the body of the batch.

## 2. The writer of segments and the fields in Python

What it gives: `popnei.write_vars_from_segments(...)` writes a vars file
from a generator of dicts; `write_vars(..., allow_unsorted=True)` and
`writeVars(variants, {allowUnsorted: true})` write unsorted variants,
which both refuse by default; `variants.var_fields`, `variants.gt_fields`
and `iter_blocks(fields=(..., "GP"))` give the fields back in
`block.extra_fields`.

Deliverables:

1. The pytest tests of "How it is verified" of "The writer of segments",
   in `tests/test_vars_from_segments.py`: the worked example of the three
   segments with every assertion of the spec, through pyarrow, `open_vars`,
   `iter_blocks`, `filter_individuals`, `write_vars` and `write_vcf`; one
   test for each error of "A call that fails", each with no file left;
   the exception of the generator; the Ctrl-C between segments; the rules
   of "Its Python function" for the types of a field array. `uv run pytest
   tests/test_vars_from_segments.py` exits with 4, no such file, on
   95c5874.
2. Pytest tests whose names contain `unsorted` of `write_vars`: a VCF
   written in the test whose second variant is before its first is the
   `ValueError` with no file, and is written with `allow_unsorted=True`.
   The two tests over `long_vcf` pass `allow_unsorted=True`. `uv run
   pytest -k unsorted` exits with 5 on 95c5874.
3. A node test whose name contains `allowUnsorted`: `writeVars` of a VCF
   that is not sorted throws by default and gives the bytes with
   `allowUnsorted`. `node --test --test-name-pattern=allowUnsorted
   test/*.test.ts`, after `npm run build`, runs no test on 95c5874.
4. Every existing test passes, no fewer than the baseline; the ones that
   change are those of deliverable 2 and those that read
   `format_version`.
5. Every check of the `coding` skill passes, `npm run build && npm test`
   and `npm run test:browser` among them.
6. The docstrings of `write_vars_from_segments`, `write_vars`,
   `iter_blocks`, `Block.extra_fields`, `Variants.var_fields` and
   `Variants.gt_fields`, and the doc comment of `writeVars`, say what the
   spec says a user has to know.

Stands on: work package 1.

Tasks:

- [ ] 2.1 `allow_unsorted` of `write_vars` in the Python binding crate and
  `python/popnei/io_vars.py`, and `allowUnsorted` of `writeVars` in
  `crates/popnei-js` and `js/popnei/src`, for both sources of each; the
  tests over `long_vcf` and `tests/test_filter_passed.py`. From
  `allow_unsorted` under "The writer" and "Sorted variants". Deliverables
  2, 3, 4, 5, 6.
- [ ] 2.2 The writer of segments: `SegmentReader` and its checks in the
  Python binding crate, the exception it keeps and raises again, the
  Ctrl-C between segments, `write_vars_from_segments` in
  `python/popnei/io_vars.py` and `__init__.py`, and the tests of
  deliverable 1 that pyarrow and `open_vars` can check without reading the
  fields through popnei: the columns, the key, the batches, every error,
  the exception, the Ctrl-C and the types of a field array. From "The
  writer of segments". Deliverables 1, 4, 5, 6. Needs 2.1, whose files it
  shares.
- [ ] 2.3 The reading of the fields in Python: `var_fields` and
  `gt_fields` of `Variants`, the extra names in `iter_blocks` with
  `select_extra_fields` on the source opened for the pass, the `extra` of
  the block handed to numpy read only, `Block.extra_fields`, its `repr`,
  the annotation of `fields`, and the rest of the worked example of
  deliverable 1: `iter_blocks` with `GP`, `filter_individuals`, `write_vars`
  of the filtered file and `write_vcf`. From "In Python and in TypeScript"
  of "The extra fields". Deliverables 1, 4, 5, 6. Needs 2.2.

What could go wrong: the pass of `write_vars` releases the interpreter,
and `SegmentReader` takes it back for each segment; a reader that holds a
Python object has to stay `Send`, as `BlockReader` asks. The exception
of the generator has to survive the error of the core that ends the pass
and the mapping of `write_the_pass`, which today names the file it read.

## At the end

The sum of the two work packages, and the code review of the `code-review`
skill on the whole branch.
