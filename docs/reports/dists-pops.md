# Work report: the distances between populations

The plan `docs/plans/dists-pops.md` is under way on the branch
`plan/dists-pops`, in the worktree `.claude/worktrees/plan-dists-pops`,
since 23 September 2026. It builds seven distances between every pair of
populations out of one pass over the variants, each with a jackknife
standard error, through the core crate, both binding crates and both
packages. Nothing is merged into `main` and nothing is pushed.

This report is written while the work goes. It has a section for each
work package as it finishes, and the top of it will say, when the plan is
done, what exists that did not and what is asked of the owner.

## What was in place before the first task

Everything the plan's "What has to be in place" names, checked by running
it on the commit this branch starts from, `6704347`:

| what | command | what it gave |
|---|---|---|
| plink2 | `plink2 --version` | `PLINK v2.0.0-a.7.7 M1 (18 Sep 2026)` |
| R | `R --version` | `R version 4.6.1 (2026-06-24)` |
| adegenet, mmod, admixtools | `Rscript -e 'packageVersion(...)'` | 2.1.11, 1.3.3, 2.0.10 |
| the reference data | `ls tests/reference/pop_dists/` | `make_reference.py`, `micro.vcf.gz`, `micro_pops.txt` and the output of every program on both panels, 17 files |
| the core module | `cargo test -p popnei --lib pop_dists:: -- --list` | `0 tests, 0 benchmarks`, exit 0 |
| the Python tests | `ls tests/test_pop_dists.py` | no such file |
| the TypeScript tests | `ls js/popnei/test/pop_dists.test.ts` | no such file |

The three programs print the versions `make_reference.py` refuses any
other of, so the literals in the tests are of the data in the repository.
The last three rows are the baseline the plan's checks are written
against: a cargo selector alone passes on an empty crate, which is why
each check names how many tests have to run.

## Work package 1, as it goes

### Tasks 1.1 and 1.2, the counts of a variant and the resampling groups

Both went to one subagent, in one prompt, because both write
`crates/popnei/src/pop_dists.rs` and one tree has one writer per file.
Three commits came back: `f310edd`, an addition to the spec, `2fc9b0a`,
task 1.1, and `f732cf3`, task 1.2. 185 497 tokens.

`crates/popnei/src/pop_dists.rs` now holds the counts of one variant in
one population, the five values a pair makes of them, and the walk that
cuts the variants into the groups the standard errors are resampled over.
`cargo test -p popnei --lib pop_dists::` gives `16 passed; 0 failed`,
against `0 tests` when the branch started, and the whole workspace gives
`537 passed; 0 failed; 2 ignored`. The other checks of the `coding` skill
were run again by the orchestrator on `f732cf3`: `cargo fmt --all
--check` no output, `cargo clippy --workspace --all-targets -- -D
warnings` and `cargo wasm-check` both `Finished`.

Deliverable 1 is met: the tests assert the spec's counts table and its per
variant table for all four variants of the worked example, the literals
read from the spec, variant 3 where both populations are fixed and H_b is
0 and variant 4 where f_2 is -0.1 among them, and the corrected H_S and
H_T that no measure reads until work package 2. Deliverable 2 is met: the
test that cuts the biallelic panel reads the file itself and asserts the
12 groups of 100 variants over its two chromosomes, with the anchoring, a
group for each variant, no group at all, and a length of 0 as an error.

### What the spec did not say, and now does

A population that has called one allele at a variant has no within
population heterozygosity there: the correction is n_P / (n_P - 1) times
1 - sum over a of p_Pa^2, which is 0 over 0 at n_P = 1, and the variant
would have added a NaN to every sum of every pair that population is in.
The spec said what the six sums do when a population has little data and
not what they do with one called allele. The variant now counts for no
pair that population is in, which is the rule the spec already had for two
populations of one called genotype each, and it is written into the item
in `f310edd`, a commit of its own before the code.

It takes a haploid dataset and a `min_num_individuals` of 1 to reach,
since a genotype of two alleles or more that was called whole gives two
called alleles or more, so no value of either panel changes. The
orchestrator took it as a small choice made and written down rather than
an open point for the owner, as the `coding` skill allows: the spec had
defined nothing here, so nothing the owner decided was overturned, and the
alternative was a silent NaN.
