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
