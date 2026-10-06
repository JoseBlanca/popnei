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
