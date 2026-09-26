# Report: the VCF writer, the filter by regions, the missing rate and the density of the variants

26 September 2026. It records how `docs/plans/writer-regions-density.md`
was carried out, on the branch `plan/writer-regions-density` in the
worktree `.claude/worktrees/plan-writer-regions-density`, which branches
from `spec/writer-regions-density` at 838459f and not from `main`, because
neither the specs nor the plan are on `main` yet.

The plan is under way.

## Before the first task

The checks of the `coding` skill, run at 838459f on 26 September 2026, gave
what the plan says they gave on 7331638. `cargo fmt`, `cargo clippy`,
`cargo wasm-check`, `cargo wasm-check-js` and `ruff` passed.
`cargo test --workspace` gave 1022 passed with 2 ignored and 150 passed;
`cargo test -p popnei --no-default-features` gave 1022 passed.
`uv run maturin develop --release && uv run pytest` gave 550 passed and 6
skipped at the first run and 551 passed and 5 skipped at the second: the
tests of `tests/test_interrupt.py` skip themselves when the interrupt
arrives after the work it was meant to stop, which depends on timing, so
the count of skipped tests can move by one or more between runs.
`npm run build && npm test` in `js/popnei` gave 444 tests, 443 passed and
the one failure the plan names, in `test/gwas.test.ts`.

bcftools 1.24, tabix 1.24, bgzip 1.24 and plink2 v2.0.0-a.7.7 are in the
PATH, and `/Users/jose/devel/popnei-bench/big.vcf`, 403 MB, and `big.vars`,
81 MB, are there. The machine is an Apple M5 Pro with 18 cores.

The board held no message of another branch that is still at work: the
three of `main` announce merges, and the two of
`spec/writer-regions-density` are this plan's own spec and the failure of
the node test above. No file of this plan is being changed on another
branch.

## 1. The reader trait and the header of a source

Task 1.1, commit edcac3d. `SourceHeader` and the three methods of
`BlockReader`, `header`, `skip_outside` and `num_skipped`, are in every one
of the 26 implementations of the trait, where the plan counted 21: the
crate has 21, among them the ones for `Box` and `&mut`, and the benches 5.
The threshold filters and the filter by linkage disequilibrium refuse the
offer of regions; the filter of individuals, `reblock`, the reader one
block ahead, `Box` and `&mut` hand it to their source. The reader one
block ahead sends the offer to its reading thread and waits for the
answer, and gives first a block the thread built before it saw the offer.
`skip_outside` takes a `RegionSelection`, which the plan builds in task
4.1, so that type and an empty `Regions` came forward into
`crates/popnei/src/filters.rs`, with no constructor outside the tests;
task 4.1 fills them in. Four tests under `filters::tests::source_header`;
the subagent broke the code three ways and each broke them.
