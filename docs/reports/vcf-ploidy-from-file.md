# Report: the ploidy of a VCF read from the file when the caller gives none

6 October 2026. The work report of `docs/plans/vcf-ploidy-from-file.md`,
on the branch `spec/vcf-ploidy-from-file`. State: under way.

## 1. The ploidy read from the file

Task 1.1, the core, is 51a860d: `ploidy_of_vcf`, `NUM_LINES_FOR_THE_PLOIDY`,
the error cases `VcfPloidyNotRead` and `VcfPloidyOfNoVariants`, and the new
words of `VcfPloidyOutOfRange` and `VcfGenotypePloidy`. `cargo test -p
popnei --lib io::vcf::tests::ploidy_of_vcf` runs 31 tests, which pass. The
search builds a `VcfReader` and reads the data lines through its source, so
the header is read by the reader's own code.

The task found two rows of the spec's table wrong, corrected in 4f42f1a.
Under the FORMAT `DP:GT` the column `0/1` holds the DP alone, so the row is
`3 3:0/1 .`. And `0/x/1 1/1/1 .` gives 3 whether `0/x/1` is skipped or
counted as three alleles, so it guarded nothing; `0/x/1 0/1 .`, which gives
2 only when it is skipped, is the guard. With the skipping taken out on
purpose, three tests failed.
