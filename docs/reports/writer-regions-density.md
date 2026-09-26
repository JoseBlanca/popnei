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
of the 24 implementations of the trait, where the plan counted 21: the
crate has 19, among them the ones for `Box` and `&mut`, and the benches 5.
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

The review ran in six categories, spec, tests, errors, api, architecture
and numbers, over 838459f..f60c4e3; the binding category was left out
because neither binding crate changed. What it found that mattered, all of
it fixed in bba6e43 to 1d3a003 and one commit after them:

- The `##contig` lines were read by a parser stricter than htslib in some
  places and looser in others, and four reviewers found it from four
  sides. A space after a comma, a blank after the `>`, or an escaped quote
  `\"` inside a Description made a length go missing with no error, or
  gave a wrong one. When one line gave `length=` twice the last one won,
  where bcftools 1.24 keeps the first. The reader now trims the blanks as
  htslib does and takes `\"` as text. It refuses, as the wrong header with
  the line, a line with two `length=` or two `ID=` fields, an empty ID, no
  closing `>`, an unclosed quote, or no `<`. `docs/specs/io_vcf.md` records
  these as decisions of 26 September 2026.
- The test that a comma inside quotes does not end a field could not
  fail: its fixture had a space before the second `length=` and the real
  length last, and a reviewer who broke the quote handling saw all 1035
  tests pass. With the fixture changed, 4 tests fail on that breakage.
- Each new chromosome name was looked up in a list, so the time to open a
  file grew with the square of its `##contig` lines. A VCF with 200000 of
  them opened in 13.37 s; it opens in 0.030 s now, with the names in a
  hash, where `bcftools view -H` takes 0.086 s (release build, Apple M5
  Pro, the fastest of 3 runs). Draft plant assemblies have 1e5 to 1e6
  scaffolds.
- The reader one block ahead stopped serving the questions of its caller
  once its source had given its last block, so an offer of regions made
  at that moment was answered true or false by the timing of the threads:
  a reviewer counted 161 of 600 rounds one way and 439 the other. The
  thread now answers until its caller lets it go, and the test that its
  held block is not lost no longer depends on the scheduler.
- An offer of regions a source accepts holds for the life of that source.
  The trait says so now, and `docs/specs/block.md` says that the variants
  a source may skip are those the selection does not keep, which with
  `exclude` are the ones inside the regions. Every pass of popnei builds
  its chain with `chain_of`, which owns its source, so no pass can meet a
  source that skips after its filter is gone; a Rust caller that builds
  the filter by regions over a borrowed source and reads the source again
  can, and the doc is what warns them.
- Smaller: an answer arriving where the reader one block ahead waits for
  a block is now an error of popnei's own, a `RuntimeError` in Python; a
  wrong `chrom_lengths` names its pair and not the whole list; the error
  for a length past 18446744073709551615 says so; and three doc comments
  were put right.

Not taken: the test sources that give a header with no individuals, which
no test reads.

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

Task 2.2, commit 833da72. `PerVarStat.MISSING_RATE` in Python and
`"missing_rate"` in TypeScript ask for the rate, and the results give it as
`missing_rate` and `missingRate`; it is among the six statistics computed
when none are named. `tests/reference/stats/make_reference.py` keeps
`many.vmiss`, `many.popA.vmiss` and `many.popB.vmiss` of plink2
v2.0.0-a.7.7 and checks the table of the spec against them.

Deliverables, run at 78a214f:

- 1: `cargo test -p popnei --lib per_var_missing_rate` gave 6 passed,
  where there were none.
- 2: the reference script ran and its check passed; the pytest test
  against the three files, over all and over popA and popB, passed.
- 3: the pytest comparison with pyNei leaves the missing rate out and
  passes; the test of one pass against one pass each covers six.
- 4: the TypeScript test of the means over all and in popA passed.
- All the checks: `cargo test --workspace` 1050 passed, 2 ignored, and
  150; pytest 555 passed, 5 skipped; `npm test` 444 of 445, the old
  failure; `npm run test:browser` 8 passed; fmt, clippy, both wasm checks
  and ruff clean.
