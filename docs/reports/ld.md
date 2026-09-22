# Work report: r², the matrix of it, and the filter by linkage disequilibrium

The plan `docs/plans/ld.md` is under way, on the branch `plan/ld`, in the
worktree `.claude/worktrees/ld`, where it started on 22 September 2026.
It builds r², the squared correlation between the dosages of two
variants; `calc_rogers_huff_r2_matrix`, which gives the r² of every pair
of a set of variants; and `Variants.filter_by_ld`, which takes out the
variants that repeat what a variant near them on the chromosome already
said. The two specs behind it are `docs/specs/ld.md` and the item "The
filter by linkage disequilibrium" of `docs/specs/filters.md`.

This report is written as the work goes. Each work package gets a section
below when it is done, with the command that checked each deliverable and
what it gave, what was changed in the plan and why, what the review found,
and what the owner should know. When the plan is done, what the owner
reads first goes at the top of this file.

## Before the first task

The branch starts from `main` at 98806a8, which is `main` with the two
specs, the plan and the programs of `docs/reports/ld-method/` merged.

Everything the plan asks to be in place is there, checked by running it:

| What the plan asks | Command | What it gave |
| --- | --- | --- |
| plink2, bcftools and R on the machine | `which plink2 bcftools R` | all three found; `PLINK v2.0.0-a.7.7 M1 (18 Sep 2026)`, `bcftools 1.24`, `R version 4.6.1 (2026-06-24)` |
| the workspace passes | `cargo test --workspace` | `395 passed` in the core crate, `35 passed` in the linalg crate |
| no r² code yet | `cargo test -p popnei --lib ld:: -- --list` | `0 tests, 0 benchmarks` |
| the filters and the variant module as the plan counted them | the same with `filters::` and `variant::` | `34 tests` and `10 tests` |
| nothing of this plan written yet | `ls` of the four paths | `tests/test_ld.py`, `js/popnei/src/ld.ts`, `tests/reference/ld/` and `crates/popnei/src/ld.rs` all absent |
| the linear algebra it goes through | `ls crates/popnei-linalg/src/` | `blas.rs`, `faer.rs`, `lib.rs`, and `docs/specs/linalg.md` merged |

Open 2 of `docs/specs/ld.md`, the major allele of a variant with half
called genotypes, has no answer, so task 1.2 follows its "meanwhile", the
rule `docs/specs/pca.md` already has: the move of `the_major_allele` into
`variant` changes no number. Open 1 belongs to the curve of r² against
distance, which this plan does not build.

## Work package 1: the r² of two sets of variants

The tasks as they are done. The deliverables, the review and what the
owner should know go at the end of the work package.

Task 1.2, `the_major_allele` public in `variant`, commit 7777f48. The
function that picks the allele a variant was called most often at was
private to `pca`, with a second copy made public for the benchmark of the
row; it is now one public function in `variant` that `pca` and the
benchmark call. Deliverable 2 holds: `grep -c "fn the_major_allele"
crates/popnei/src/pca.rs` prints `0`, and `cargo test --workspace` passes
`395 passed` in the core crate and `35 passed` in the linalg crate, the
counts the plan started from.

One thing was not a plain move. The test that a variant whose two alleles
were called equally often takes the lower numbered one asserted, in
`pca`, the standardized row of four individuals, and it uses internals of
`pca` that `variant` does not have. It now asserts the allele those same
four genotypes give, and two cases beside it that `pca` never reached: a
variant with a half called genotype, and one of which no allele was
called, whose major allele is the missing one. The tests of `variant::`
go from 10 to 11 and those of `pca::` from 48 to 47, so the workspace
count does not move. The plan's deliverable says that no test changes;
this one did, and what it covers grew rather than shrank.
