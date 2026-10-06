# Report: the ploidy of a VCF read from the file when the caller gives none

6 October 2026. The work report of `docs/plans/vcf-ploidy-from-file.md`,
on the branch `spec/vcf-ploidy-from-file`. State: done.

The plan is done, on the branch `spec/vcf-ploidy-from-file`, not merged.
`open_vcf` in Python and `openVcf` in TypeScript now read the ploidy from
the VCF when the caller gives none: it is the number of alleles of the
first genotype that is not a single dot, in the first 4096 data lines, and
the opening fails with words that ask for the ploidy when none is there.
A VCF with a header and no data line, opened with no ploidy, is refused
with "the file has no variants and the ploidy can't be inferred", as you
decided. A ploidy that is given means what it meant before. Below, "the
search" is that reading of the ploidy when the file is opened, and "a
pass" is one reading of the file by a calculation.

Every check of the `coding` skill passes, and each runs more tests than
before the plan. One node test fails, `a kinship that does not tell the
two variances apart gives none of them`, which failed in the same way on
`main` at df92324 before the plan and is not of it. No decision of the
spec is left open; two small things were seen and not taken, under "What
the review found". What is asked of you is the order to merge the branch
into `main`, and with it issue 8 can be closed.

## 1. The ploidy read from the file

Task 1.1, the core, is 51a860d: `ploidy_of_vcf`, `NUM_LINES_FOR_THE_PLOIDY`,
the error cases `VcfPloidyNotRead` and `VcfPloidyOfNoVariants`, and the new
words of `VcfPloidyOutOfRange` and `VcfGenotypePloidy`. `cargo test -p
popnei --lib io::vcf::tests::ploidy_of_vcf` runs 31 tests, which pass; the
fixes of the review brought them to 34. The
search builds a `VcfReader` and reads the data lines through its source, so
the header is read by the reader's own code.

The task found two rows of the spec's table wrong, corrected in 4f42f1a.
Under the FORMAT `DP:GT` the column `0/1` holds the DP alone, so the row is
`3 3:0/1 .`. And `0/x/1 1/1/1 .` gives 3 whether `0/x/1` is skipped or
counted as three alleles, so it guarded nothing; `0/x/1 0/1 .`, which gives
2 only when it is skipped, is the guard. With the skipping taken out on
purpose, three tests failed.

Task 1.2, the Python side, is f26ee06, and task 1.3, the TypeScript side,
is 53ff7e3; they ran side by side. Each found tests the plan did not list
that open a VCF with a header and no data line to check a pass with no
variant, one in pytest (`test_pcoa.py`) and two in node (`ld.test.ts`,
`pcoa.test.ts`). They give the ploidy 2 like the 16 the plan named, under
the plan's own rule for such a test, and the spec and the plan count 19
(d15e485, 8f2e5e6).

### The deliverables

Each run on 2a6c486, after the fixes of the review, on 6 October 2026.

| deliverable | command | result |
|---|---|---|
| 1, the cargo tests of the search | `cargo test -p popnei --lib io::vcf::tests::ploidy_of_vcf -- --list` | 34 tests, 0 on df92324 |
| 2, every test passes | `cargo test --workspace`; `cargo test -p popnei --no-default-features`; `uv run pytest`; `npm test` | 1453, 1303 and 741 passed; node 551 tests, 550 pass, 1 fails, the kinship test that fails on `main` |
| 3, the pytest tests | `uv run pytest -k ploidy_from_the_file` | 7 passed, exit 5 on df92324 |
| 4, the node tests | the spec reporter with `--test-name-pattern="ploidy from the file"` | 8 tests by name, 0 on df92324 |
| 5, the browser | `npm run test:browser` | 9 passed, with a `File` of `tetraploid.vcf.gz` opened with no ploidy |
| 6, the API | `inspect.signature(popnei.open_vcf)`; `hasattr(_core, "DEFAULT_PLOIDY")` | `ploidy: int \| None = None`; `False`; `default_ploidy` is in no file of the wasm package |

`cargo fmt --check`, `cargo clippy -D warnings`, `cargo wasm-check` and
`cargo wasm-check-js` pass. Before the plan, on df92324, the four suites
ran 1419, 1269 and 734 tests that passed and 543 node tests.

### What the review found

Six reviewers, in the categories spec, tests, errors, api, binding and
architecture, read 42873fa..8f2e5e6. None found a wrong ploidy or a wrong
result. These held and are fixed, each in its own commit after the spec
was corrected in 05b5da1:

- The search had its own copy of the reader's parsing of a line and a
  genotype, and the copy had drifted: the reader takes the line end off
  twice, the search once, so a file of one individual whose line ends in
  `0/1\r\r\n` was read by the reader with the ploidy 2 and refused by the
  search. Found by the architecture reviewer with every genotype of up to
  five tokens of a small alphabet, 111111 of them. A later change to the
  reader not copied to the search could have given a wrong ploidy with no
  error. The reader and the search now share four functions (ffeccfb), and
  a cargo test over 3695 genotypes, each with three line ends, binds the
  two: when the search gives p, the reader reads the line with p and with
  no other ploidy. It takes 0.09 s and failed on the code before. The plain
  403 MB `big.vcf` read on one thread took 0.616 s before and 0.597 s after,
  the median of 7 runs, within the noise of the machine.
- The row `/. 1 .` → 1 gave 1 however `/.` was read, so a reordering of
  the separator strip and the dot check would have passed it; it is
  `/. 0/1 .` → 2 now, and a genotype of 255 alleles has a test, which a
  bound written `>=` fails (67cd1ac).
- The words "the first 1 data lines" and "the ploidy of the VCF is 0", the
  second for a ploidy the caller gave, are replaced (be0e791).
- The Python binding opened the file for the search itself, with a buffer
  of 8 KiB where a pass has 256 KiB and a copy of the error of a file that
  could not be opened. The core has `ploidy_of_vcf_at`, which shares the
  opening of `from_path` (7c11a26).
- When the search gave an error that carries no path, the Python binding
  made it a `ValueError` with the path, whatever the error was, so a
  defect of popnei passing that way later would have read as a wrong
  input. The binding now does this for the one error of that kind the
  search gives, a ploidy above 255 read from a genotype, and sends every
  other error to the usual mapping of errors to exceptions (72b87eb).
- The doc comment of `openVcf` said that both refusals ask for the ploidy,
  and the one of a file with no variant does not; a node test of a
  genotype of 256 alleles; the bound linked to its constant in the Rust
  doc comments; three lines reflowed (0abac28).
- `.claude/skills/coding/pyo3.md` counted thirteen cases of the Python
  error, and there are fourteen (2a6c486).

Not taken: that a failed `openVcf` frees its entry in the table of sources
is not tested, here or before this plan, and the node tests have no way to
see that table. The reader reads `0/1\r\r` as `0/1`, which is lenient and
not wrong, and is left as it is.

### What the owner should know

- With no ploidy given, a corrupted or cut file whose damage comes before
  the first genotype with alleles fails at `open_vcf` and `openVcf`, where
  it failed at the first pass. The tests of a corrupted second member give
  the ploidy 2 to keep checking the error of a pass, and new tests check
  the error at the call.
- Opening a file now reads it twice up to its first genotype with alleles.
  The architecture reviewer timed the worst case the spec names on the
  owner's machine, with a release build, 100000 individuals and 4095 lines of single dots, 820 MB: the
  search took 1.67 to 1.93 s, and a whole pass 0.43 s on all threads and
  3.0 s on one. A file whose first line has a called genotype opens in
  about the time of a read of its header.
- `python/popnei/io_vcf.py`'s docstring of `open_vcf` changed here and on
  the branch `docs/sphinx`; the second merge meets it.

### How the work went

This part is for whoever next revises a skill or writes a plan, and the
owner can stop here.

- The spec's table of cases had three rows that were wrong or guarded
  nothing, `DP:GT` with `3 0/1 .`, `0/x/1 1/1/1 .` and `/. 1 .`, all
  written by the session that wrote the spec and none found by the spec
  reviewer, the first readers or the session. The implementer found two
  and the tests reviewer, with mutants, the third. A row meant to guard a
  rule has to give another answer when the rule is broken, and the
  `writing-specs` skill does not ask the writer to check that.
- The review of the spec counted 16 tests that the owner's decision would
  change, by searching for files with a header alone; the work found 19.
  A count of the tests a decision changes is better made by running the
  suite with the change than by searching.
- Tokens, as the `Agent` tool reported them: task 1.1, 162 thousand; 1.2,
  122 thousand; 1.3, 119 thousand; the six reviewers, 472 thousand
  together; the fixes, 165 thousand. The review cost more than the three
  tasks together, 472 against 403 thousand, and it found the one finding
  that could have become a wrong result, the drift of the search.
