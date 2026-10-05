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

Three of the 17 tests of 1.1 read a source that never ends and kept every
block. Written first, against the stub of the filter, they grew without bound,
and cargo ran them at once until the owner's machine, 64 GB, ran short of
memory and the owner had to quit other programs. The source now panics past
100 blocks, with a test that a pass nothing ends fails there (a159231), and
the coding skill asks a test written first to fail without hanging or
growing (50e4def).

Task 1.2, Python, is 07ae287: `filter_first_n`, `PassStats.stopped_early`,
every method that adds a step through one `Steps::add`, and the counts of a
pass built by one function, `pass_counts_of`, which the 13 modules of the
crate that return counts call. Task 1.3, TypeScript, is d5b869d: the same
with `filterFirstN`, `Steps::add` and `PassCounts::of`, called by the 15
consumers. The two ran side by side. Both refuse a `num_vars` of 0 with the
core's case, whose message does not name the argument, since the core refuses
it in `FirstNReader::new` and no reader exists when the step is added.

The deliverables, checked on 5 October 2026 at 07ae287: the six cargo
commands pass, `cargo test --workspace` with 1390 passed and `cargo test -p
popnei --no-default-features` with 1240; `cargo test -p popnei --lib first_n
-- --list` lists 17; `uv run pytest` gives 692 passed, 23 of them of `-k
first_n`; `npm test` gives tests 518, pass 517, fail 1, the test of
`test/gwas.test.ts` that fails on `main`; the tests named `filterFirstN` pass;
`npm run test:browser` gives 8 passed.
