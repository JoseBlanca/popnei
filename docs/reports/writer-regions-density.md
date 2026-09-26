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

Tasks 1.2 and 1.3, commits f6393b4 and 3ba0021, 2c46528 and 343bc4f, and
f60c4e3. The VCF reader keeps the lines of its header before `#CHROM` and
the length of each `##contig` line, in the order of the file.
`tests/reference/vcf/write.vcf` is the file of "How it is verified" of the
writer, and gives chr1 2000 and chr2 1500 and its nine meta lines. A vars
file keeps the lengths in `chrom_lengths` of its `popnei` key and says
`format_version` 1.1; a file of 1.0 reads with no lengths.

Six cases the spec did not settle were settled with the code and written
into `docs/specs/io_vcf.md` and `docs/specs/io_vars.md`. A length is
written in digits alone, so `+5` is the wrong header, as `0` and `abc`
are. Two lines of one ID with the same length give it once. A comma
inside quotes, in a `Description="a, b"`, does not end a field. A
`##contig` line with a length and no ID is the wrong header. The
subagent had first made this last one give no length and no error; it was
sent back, because the VCF format requires the ID and an error in popnei
never passes silently. In a vars file, a `chrom_lengths` that is not a
list of names with whole numbers above 0, or that gives one chromosome
twice, makes the file not a vars file.

Deliverable 1 said that no test of before would change but for the
version. The longer key made a few more change, none of them an
assertion made weaker: the literals that hold the key whole gained
`"chrom_lengths": []`, and the sizes of the vars files in
`js/popnei/test/progress.test.ts` grew by 16 and 24 bytes. The plan says
so now.

Deliverables, run at f60c4e3:

- 1: `cargo test --workspace` gave 1035 passed, 2 ignored, and 150
  passed; `--no-default-features` 1035 passed; pytest 551 passed and 6
  skipped, which is the 551 of before, one test added, and one more skip
  of the interrupt tests by timing; `npm test` 443 of 444, the old
  failure. fmt, clippy, both wasm checks and ruff clean.
- 2 and 3: `cargo test -p popnei --lib source_header -- --list` counts 13
  tests, where it counted 0 before, and all 13 pass.

## 2. The missing rate

Task 2.1, commit ecfc6cf. The missing rate is the sixth statistic of
`calc_per_var_distribs`: the missing genotypes of a variant over its
called and missing ones, a half called genotype counted as missing, from
the counts the observed heterozygosity already makes. Six tests have
`per_var_missing_rate` in their names: the worked example of the spec with
and without populations and at `min_num_individuals` 20, the rates 0, 0.25,
0.5, 0.75 and 1 on the edges of 4 bins, 3 of 20 in bin 5 of 40, and the
table of plink2 for `many.vcf` over all and in popA and popB, whose
histograms match exactly and whose means match within 1e-12 relative.
Both binding crates only gained the field in their destructuring, so that
they compile; until task 2.2 `stats: ["missing_rate"]` in TypeScript is
taken and gives nothing. The spec's paragraph on `MAX_NUM_BINS` counted
four statistics with a histogram and 3.2 MB for each chunk of rows; with
the missing rate they are five and 4 MB, and the spec says so now.
