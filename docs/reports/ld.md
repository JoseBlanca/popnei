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

Task 1.1, the reference dataset, commit f714f79. `tests/reference/ld/`
holds `make_reference.py` and `example.vcf`, both moved from
`docs/reports/ld-method/` and byte for byte what they were, which is what
keeps every literal of the two specs valid: the program's numbers come
from `numpy.random.default_rng(29)` and depend on the order in which it
asks for them. Beside them are `ld.vcf.gz`, the two matrices plink2 gives
for the dataset and for the worked example, the identifiers of their rows
in their order, and `run_plink2.sh`, which makes the five files again and
compares each with the stored one.

Deliverable 1 holds. `tests/reference/ld/run_plink2.sh` into an empty
directory exited 0, which is the script saying that none of the five
files differed, and the `ld.vcf` it wrote is what `gzip -dc` gives from
the stored `ld.vcf.gz`. Read back from the stored matrix of the dataset:
500 variants, 124750 pairs of which 93096 have an r² and 31654 are NaN,
68 variants with no variance, and the pair chr1:1000 with chr1:2000 at
0.353466669239891 and with chr1:3000 at 0.39849991080910563, the bits the
spec's table has.

Two things the task decided that the plan did not say. `example.vcf` was
moved rather than copied, so there is one copy of it, and
`docs/reports/ld-method/README.md` says where both files went. plink2's
`.bin.vars` files are stored although deliverable 1 does not name them:
they carry the order of the rows of each matrix, which the cargo test of
deliverable 4 has to read the matrix in. plink2's `.log` is not stored,
since it carries the time of the run and the paths of the machine.

`docs/specs/ld.md` still pointed the reader at
`docs/reports/ld-method/make_ld.py`, which the move had taken away. The
sentence now names `tests/reference/ld/make_reference.py` and the script
beside it: commit 7b6ce62. No value and no open point of the spec moved.

Task 1.3, `LdDosages`, commit 3cb0968. `crates/popnei/src/ld.rs` is new
and holds the three matrices of "How it runs" of `docs/specs/ld.md` built
over a block and a set of individuals, with `rows`, `has_variance`,
`dosages` and `maf`. `cargo test --workspace` gives `410 passed`, the 395
it started from and 15 of this task, and `cargo test -p popnei --lib ld::
-- --list` prints `15 tests`, which is deliverable 6 moving off `0
tests`. The tests were each seen to fail with the code broken.

The task found a fifth way `of_block` can refuse a block that the spec's
list of four did not have: genotypes of more than 255 alleles each. A
dosage is how many alleles of a genotype are not the major allele of its
variant, so it is at most the ploidy, and the spec's own `dosages` gives
a dosage as a `u8`. Such a block would have wrapped the dosage into a
wrong number instead of being refused. Nothing a user does reaches it
today, since the VCF reader takes 255 alleles in a genotype at most. The
spec's error list now names it: commit b28b115.

`linalg::product` has no transpose. It works out `c = a b` with `a` of
`rows` x `inner` and `b` of `inner` x `cols`, both row after row, and the
three matrices are variants x individuals, so the product of one set of
variants with another sums over the individuals and needs the second
matrix as individuals x variants. This is the risk the work package names
under "What could go wrong": `product` was built for the principal
component analysis, whose result is individuals x individuals, and here
the result is variants x variants. Task 1.4 transposes the second matrix
inside `ld.rs`, which for a tile of 512 variants of 1000 individuals is
three matrices of about 4 MB copied. The other way, a product in
`crates/popnei-linalg` that takes its second operand transposed, changes
`docs/specs/linalg.md`, which is a spec this plan does not build. Work
package 4 measures what the transposes cost against the 0.50 s of "Speed"
of `docs/specs/ld.md`, and the owner decides on a performance review then
if the target is missed, as the plan says.

Two more things left for work package 4 to measure rather than guessed at
now: `rows` gives an owned `LdDosages`, as the spec's signature says, so
each tile copies three matrices, 12 MB for 512 variants of 1000
individuals; and the rows are built by one serial loop, where `pca.rs`
has a loop per ploidy with the length of a genotype known at compile
time.
