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
