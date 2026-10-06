# Report: the filter of the variants that passed their FILTER

6 October 2026. The work report of `docs/plans/filter-passed.md`, on the
branch `spec/filter-passed`. State: under way.

## 1. The column of whether each variant passed

Task 1.1, the column in the core, is 53e8b64: `Needs::PASSED` in `ALL`,
`Block.passed`, carried by `check`, `retain_vars` and `reblock`, and filled
by the VCF reader through `is_a_pass`, the one test of FILTER that
`only_passed` and the column share. `cargo test -p popnei --lib
passed_column` runs 8 tests, which pass. A read of the plain 403 MB
`big.vcf` on one thread, median of 7 runs, took 0.634 s before and 0.578 s
after with the genotypes alone, and 0.641 s and 0.581 s with every field;
the bound of the plan was 3% slower. Besides the seven tests that compare a
set of fields with `ALL`, one test changed: a test of `check` that clears
every column but the text of the lines now clears `passed` too.

Task 1.2, the vars file of format 1.2, is 80394e8, with e0a3c72 after it.
The writer adds `passed` when its source's blocks have it, and the reader
reads it when the file has it; a null in it is an error that names the
column and the variant, which the spec did not say and now does (3ecff37).
`passed_column` runs 17 tests. Two node tests of `progress.test.ts` held the
bytes of a vars file, which the column makes larger, and the plan did not
name them (fb76b9c): written with and without the column, the large file of
that test grows by 18024 bytes and the one of `many.vcf` by 232, all of it
the column and 64 bytes of each its schema, which a pass does not read.

Seen on the way and not of this plan: a vars file written from node is
larger than the same file written from Rust, 43962 bytes against 39802 for
`many.vcf` and 12249642 against 10164842 for the large file of
`progress.test.ts`. Nobody has looked at why.

Task 1.3, the VCF writer of a vars file, is 5a53969 and df5ff93. It writes
`FAIL` for a variant whose `passed` is false and the `##FILTER` line when
the source keeps the column, which the writer learns from a field the
task added to the header of a source, `keeps_passed`; the spec had no such
field and now has it (2e2eb10). A block with a failed variant from a source
whose header says it keeps no `passed` is a defect, a `RuntimeError`, which
the spec did not say either. `many.vcf` read with every variant, written to
a vars file and that to a VCF, gives 25 lines with `FAIL` and 475 variants
from `bcftools view -H -f .,PASS`.

### The deliverables

Run on b0e10f3, before the fixes of the review, and again after them on
0abb855.

| deliverable | command | result |
|---|---|---|
| 1, the cargo tests of the column | `cargo test -p popnei --lib passed_column -- --list` | 20, then 21 after the fixes; 0 on 6abacd1 |
| 2, every test passes | the six cargo commands, `uv run pytest`, `npm test`, `npm run test:browser` | 1473 then 1474 passed; 1323; 741; node 551 tests, 550 pass, 1 fails, the kinship test that fails on `main`; browser 9 |
| 3, `write_vars` | `uv run pytest tests/test_io_vars.py -k write` | 17 passed, among them the seven columns, `"1.2"` and the 25 false |

### What the review found

Six reviewers, spec, tests, errors, api, binding and architecture, read
02cdae0..b0e10f3. None found a variant marked passed when it failed, or the
reverse. These held and are fixed, after the specs were corrected in 0fa5af6:

- The docstrings of `write_vcf` and `writeVcf` said that FILTER is a dot
  for a vars source, and those of `write_vars` and `writeVars` listed six
  columns (054105a).
- The defect of a failed variant from a source that keeps no `passed`
  named nothing; it names the chromosome and the position (b684b8d).
- The test of that defect held a block with variants that passed and
  failed, so a check broken to refuse every block with the column, or none,
  passed it; it has a block of each now (0e92f7a).
- Two tests of the reader left the column out: the one of a window into a
  batch, which a reader that ignored the offset of the bits would have
  passed, and those of the threads of the VCF reader (409357f).
- A VCF written from a vars file through the filter of individuals is
  tested to keep the `##FILTER` line, and three comments count as the code
  does now (0abb855).
- The specs said that there was no code, counted five flags of `Needs`, left
  `passed` out of the columns written with no nulls, and did not say that a
  vars file written from a source of no variant has no `passed` column,
  which nothing then reads (0fa5af6).

The timing of 1.1 above is weaker than it reads: with the genotypes alone
the column is not filled, so that pair measured the noise between runs, and
the pair with every field is the one that sees it. The architecture
reviewer timed every field with and without `PASSED` on the same `big.vcf`,
15 interleaved runs on one thread on a machine with a load of 24: medians
of 1.029 s and 1.032 s, no difference that can be measured at that noise.

### How the work went

This part is for whoever next revises a skill or writes a plan, and the
owner can stop here.

- Two reviewers sent to read only, architecture and binding, changed files
  of the shared worktree to try mutations, and each saw the other's change
  for a while; both restored theirs. The `code-review` skill says that the
  categories other than spec and tests only read and share the checkout,
  and a reviewer that wants to run a mutant does not keep to it. Every
  reviewer that may build should get a worktree of its own.
- The specs missed four things the code needed: the field of the header,
  `keeps_passed`; the rule for a null in `passed`; the defect of a failed
  variant under a header without the column; and two node tests holding
  the bytes of a vars file. Each was found by the task and put in the spec
  before or with the code.
