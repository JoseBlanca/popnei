# Work report: the filter of the first n variants and the filter that keeps variants at random

5 October 2026. The report of `docs/plans/filters-random-first-n.md`, which
builds issues 7 and 6 of the repository from `docs/specs/filters.md`, on the
branch `filters-random-first-n` from `main` at eae29a2. State: under way.

## Before the first task

The checks of the plan's "What has to be in place" were run on 5 October
2026 and gave what the plan says. The board had no message of another
branch about the steps or the counts of a pass; `spec/writer-regions-density`
still has a worktree although the board says it was merged, which is not of
this plan.

## 1. The filter of the first n variants

Task 1.1, the core, is c5a96e0: `FirstNReader`, `PassStep::FirstN`,
`stopped_early`, `refuse_a_step_after_the_first_n` and three cases of the
error, in `crates/popnei/src/filters/first_n.rs`. No existing reader needed a
fix: with the chain on the thread of the reader one block ahead, the filter's
`None` ends that thread with the source asked for 2 blocks, and no test hung.
The VCF of 100000 variants of 10 individuals, 6989013 bytes, was read for
8192 bytes by a pass with the first 100, and the vars file for its batch 0
alone. Checked on 5 October 2026: the six cargo commands pass, `cargo test
--workspace` with 1389 passed and `cargo test -p popnei --no-default-features`
with 1239, 16 more than before in each, the 16 of `cargo test -p popnei --lib
first_n`. One test asserts counts the spec did not give, the MAF filter at 0.8
before the first 10 in blocks of 7, 14 given and 12 kept, which bcftools
confirmed and which went into the spec in e2f253d.
