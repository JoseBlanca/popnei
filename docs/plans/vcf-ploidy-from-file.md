# Plan: the ploidy of a VCF read from the file when the caller gives none

6 October 2026. State: **under way**, approved by the owner on 6 October 2026. It builds what `docs/specs/io_vcf.md`
gained on 6 October 2026 from issue 8 of the repository: with no ploidy
given, `open_vcf` and `openVcf` read the ploidy from the first genotype
of the file that has alleles. The parts of the spec are the paragraphs on
the ploidy in "What it gives" of the reader, from "The caller may give the
ploidy", the ploidy bullet of "Its Python and TypeScript functions", the
paragraphs from "The ploidy read from the file." in "How it is verified"
and in "The Rust interface". The spec has no open point: the owner
decided on 6 October 2026 that a VCF with a header and no data line,
opened with no ploidy, is refused. It runs through the core
crate, both binding crates, the Python package and the TypeScript package,
in one work package.

## In and out

Built: `ploidy_of_vcf` and `NUM_LINES_FOR_THE_PLOIDY` in the core, the two
new error cases and the new words of two old ones, the call of `ploidy_of_vcf`
in both binding crates when no ploidy is given, a `ploidy` of `None` or
left out in both packages, their doc comments, and the tests the spec
lists.

The tests that open a VCF with a header and no data line and give no
ploidy, 16 of them, now give the ploidy 2, since that file is refused at
the opening without one and what they check is a pass with no variant: in
pytest `test_block.py:297`, `test_dists.py:611` and `:633`,
`test_diversity.py:824`, `test_io_vars.py:320`, `test_ld.py:561`,
`test_pass_stats.py:152`, `test_stats.py:588` and `:1100`; in node
`dists.test.ts:464` and `:481`, `diversity.test.ts:879`,
`pass_stats.test.ts:212`, `pop_dists.test.ts:822`, `vcf.test.ts:375` and
`stats.test.ts:999`, the lines at df92324.

## What has to be in place

- The branch `spec/vcf-ploidy-from-file`, in the worktree
  `.claude/worktrees/vcf-ploidy-from-file`, from `main` at df92324, with
  the spec at bd4a533 or later, and its start message on the board,
  `.claude/board/2026-10-06T0900-spec-vcf-ploidy-from-file.md`. The plan
  runs on that branch.
- The shared files it changes, which other branches may also change: the
  error enum `crates/popnei/src/error.rs`, `crates/popnei/src/io/vcf.rs`,
  `crates/popnei-python/src/vcf.rs` and `lib.rs`,
  `crates/popnei-js/src/vcf.rs`, `python/popnei/io_vcf.py`, whose
  docstrings the branch `docs/sphinx` also edited, and
  `js/popnei/src/io_vcf.ts`.
- The checks of the `coding` skill on df92324, run on 6 October 2026, whose
  counts are below under "Baseline". "No fewer" in the deliverables counts
  from them. A fresh worktree needs `npm ci` and `npm run build` in
  `js/popnei` before `npm test`.

## 1. The ploidy read from the file

What it gives: `open_vcf("tetraploid.vcf.gz")` and `openVcf(bytes)` with
no ploidy give a `Variants` whose `ploidy` is 4, and a file whose first
4096 data lines hold no genotype with alleles is refused at the call with
words that say to give the ploidy. A caller who gives the ploidy gets
what they get today, but for the new words of two errors.

Deliverables:

1. Cargo tests in a module `ploidy_of_vcf` of the tests of
   `crates/popnei/src/io/vcf.rs`, one for each row of the table of "The
   ploidy read from the file." in "How it is verified", one for each of the
   three files of the repository the spec names there, one for a header
   with no data line, and the tests of the words of the three errors:
   `cargo test -p popnei --lib io::vcf::tests::ploidy_of_vcf -- --list`
   counts 26 or more. It counts 0 on df92324.
2. Every existing test passes, no fewer than the baseline. The only tests
   of before that change are those that assert the words of the error of
   a genotype of another ploidy or of a ploidy out of range; the pytest
   test that a tetraploid file is refused when its blocks are asked for,
   which now passes `ploidy=2` to `open_vcf`; and the tests of a corrupted
   second member, which pass the ploidy 2, as the spec says; and the 16
   tests of a VCF with no data line named under "In and out", which pass
   the ploidy 2.
3. Pytest tests whose names contain `ploidy_from_the_file`, at `open_vcf`,
   with the checks the spec gives Python: `tetraploid.vcf.gz` with `ploidy`
   4 and 200 variants of four alleles, `haploid.vcf.gz` with `ploidy` 1,
   the file of 4096 lines of single dots refused at the call with the path
   first in the message, and the corrupted second member refused at the
   call. `uv run pytest -k ploidy_from_the_file` exits with 5 on df92324,
   pytest's code for no test collected.
4. Node tests whose names contain `ploidy from the file`, with the checks
   the spec gives TypeScript, among them the corrupted second member at the
   call. `node --test --test-reporter=spec --test-name-pattern="ploidy
   from the file" test/*.test.ts | grep -c "ploidy from the file"`, in
   `js/popnei` after `npm run build`, prints 0 on df92324. The count of
   tests that node prints is no check here: it counts every test file as
   a test, 33 with a pattern that matches nothing.
5. `npm run test:browser` in `js/popnei` passes, with the case of a `File`
   read in a web worker opened with no ploidy and its `ploidy` checked.
6. `python -c "import popnei, inspect;
   print(inspect.signature(popnei.open_vcf))"` prints `ploidy: int | None
   = None`, the docstring of `open_vcf` and the doc comment of `openVcf`
   say that a ploidy left out is read from the file and that a corrupted
   file can then fail at the call, and `DEFAULT_PLOIDY` is no longer in
   `popnei._core` nor `default_ploidy` in the wasm package.

Stands on: nothing of this plan.

Tasks:

- [x] 1.1 The core, in `crates/popnei/src/io/vcf.rs` and
  `crates/popnei/src/error.rs`: `NUM_LINES_FOR_THE_PLOIDY`,
  `ploidy_of_vcf`, which reads the header with the code and the errors of
  `VcfReader::new` and not with a copy of it, the two error cases of the
  search at the end of the enum, the one of a file with no data line with
  the owner's words, and the new words of `VcfPloidyOutOfRange` and
  `VcfGenotypePloidy`, from the paragraphs of the spec named in the
  opening. Deliverable 1, and 2 for the core. The genotype the search
  counts has to be one the reader would read as a genotype. A genotype it
  miscounted, `0/x/1` counted as three alleles in a file of two, gives a
  wrong ploidy that no error follows in a calculation that reads the
  chromosome and the position of each variant and no genotype, the density
  of the variants: that is the silent failure of this plan. The rows of
  the table with `0/x/1`, a column with no `GT` value and a line of eight
  columns, each of which the search has to skip, are its guard, and they
  go in the same commit as the search, before 1.2 and 1.3 build on it.
- [ ] 1.2 The Python side. In `crates/popnei-python/src/vcf.rs`, `open_vcf`
  takes a ploidy that may be `None` and then calls `ploidy_of_vcf` on the
  file opened for it, with its errors made into those of the file as the
  rest of `open_vcf`'s are; `DEFAULT_PLOIDY` goes out of `lib.rs`. The new
  error case is a `ValueError`: `errors.rs` sends a case it does not list
  to `ValueError`, and the subagent checks that the message starts with the
  path. In `python/popnei/io_vcf.py`, the signature and the docstring.
  The tests of deliverables 2 and 3 for Python, and 6. Needs 1.1.
- [ ] 1.3 The TypeScript side. In `crates/popnei-js/src/vcf.rs`,
  `open_vcf` and `open_vcf_of_a_file` take an `Option<usize>` and, when it
  is `None`, `the_vcf_of` calls `ploidy_of_vcf` on an opening pass of the
  file, `the_opening_pass`, before the reader of the opening; an error of
  that call frees the entry of the source as an error of the reader does.
  `default_ploidy` goes out. In `js/popnei/src/io_vcf.ts`, `openVcf`
  passes `undefined` through, and its doc comment and that of
  `OpenVcfOptions.ploidy`. The tests of deliverables 2, 4 and 5, and 6 for
  TypeScript. Needs 1.1; it can run beside 1.2, whose files it does not
  touch.

What could go wrong: a test that opens a VCF that is not diploid, or that
the search cannot read through, with no ploidy now gets the ploidy of the
file or an error at `open_vcf`. The review of the spec searched the tests
and found only those of deliverable 2. A test that fails for this reason
is changed to give the ploidy when what it asserts happens after the
opening, an error or a result of its blocks; one that asserts what the
opening gives keeps the new result and its assertion changes. A test that
fails for another cause is reported to the owner and not changed.

The work package is done when its six deliverables hold, the checks of
the `coding` skill pass with no fewer tests than the baseline, and the
review of the `code-review` skill has been made and its findings
evaluated, as the `following-plans` skill says.

## Baseline

The counts of the checks of the `coding` skill on df92324, 6 October 2026.

| check | on df92324 |
|---|---|
| `cargo test --workspace` | 1419 passed, 0 failed, 2 ignored, summed over its `test result` lines |
| `cargo test -p popnei --no-default-features` | 1269 passed, 0 failed, 2 ignored |
| `uv run maturin develop && uv run pytest` | 734 passed |
| `npm run build && npm test` in `js/popnei` | 543 tests, 542 pass, 1 fails |

The one that fails is `a kinship that does not tell the two variances
apart gives none of them` of `test/gwas.test.ts`, which failed on `main`
before this plan and is not of it. `cargo fmt`, `cargo clippy`, the two
`cargo wasm-check` and `npm run test:browser` were not run for the
baseline; they run at the end of each task.
