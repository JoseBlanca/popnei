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

Task 1.4, `r2_between`, commit adfd3ce. The r² of every pair of two sets
of variants, the six products through `linalg::product` and the formula
over the six sums, with NaN where "What it gives" of `docs/specs/ld.md`
says there is none. `cargo test --workspace` gives `419 passed`, nine
more than task 1.3 left, and `ld::` stands at 24 tests.

Deliverable 3 holds, and it holds more strictly than it asks. The seven
pairs of the worked example were asserted within 1e-12 relative, as the
plan says, and every one of them came out equal to the spec's decimal bit
for bit, on the system BLAS and on faer alike: `cargo test --workspace
--no-default-features` gives the same `419 passed`. So the six sums are
being held exactly as whole numbers, which is what the work package leans
on and what makes a difference of 1e-13 in deliverable 4 worth looking
at. The tests were seen to fail: swapping one of the six products breaks
four of them.

What the transposes cost, for a tile of 512 variants of 1000
individuals: 12.3 MB copied for two different sets, three matrices of 4.1
MB, and 12.4 MB for a set against itself, which is four products and four
transposes, two of the dosages at 4.1 MB and two of the sums at 2.1 MB.
Work package 4 measures whether that shows against the 0.50 s target.

Left for the review of the work package. `r2_between` adds three errors
that are defects of popnei and not wrong input from a user: an `out`
buffer of the wrong size, two sets built over a different number of
individuals, and a linear algebra operation that did not run. All three
fall through the wildcard of `crates/popnei-python/src/errors.rs` and
would reach Python as `ValueError`, where the matching `PcaLinalg` of the
principal component analysis is a `RuntimeError` and
`docs/specs/pca.md` says why: no argument of the function can give them.
The same holds for `LdRowsNotInTheDosages` of task 1.3. Nothing in Python
reaches any of them today, since the binding of this module is task 2.2,
so this is not a wrong exception a user can see yet. The sentence of
`docs/specs/ld.md` that lists the cases of this module names three, all
of them reachable from an argument, and does not name these four.

Task 1.5, the two checks against the stored numbers, commits 0e4a659 and
7fb8040. `cargo test --workspace` gives `421 passed` and `ld::` stands at
26 tests, the same on the system BLAS and on faer.

Deliverable 4 holds and says more than it asks. Of the 124750 pairs of
`ld.vcf.gz`, 93096 have an r² and 31654 are NaN, the counts the plan
gives, and the test asserts both so that a pair one library gave a number
for and the other did not falls in neither count and fails. Of the 93096,
**none differs from plink2 at all**: the orchestrator set the tolerance
of the test to exactly 0.0 by hand, reran it, saw it pass and restored
the file. So the six sums are being held exactly as whole numbers over a
real dataset of 500 variants of 100 individuals and not only over the
five variants of the worked example, and the spec's warning that a
difference of 1e-13 is worth looking at stands with something behind it.
The test also asserts that the variants of the VCF are the rows of
plink2's matrix in the same order, rather than assuming it; they are,
v0000 to v0499. It runs in 0.01 s.

Deliverable 5 holds: every one of the 25000 dosages of
`tests/reference/vcf/many.vcf` is the one pyNei's `to_012` gives.

The dosages pyNei gives are stored as
`tests/reference/ld/many.pynei.dosages.tsv`, in the reference directory
with everything else `make_reference.py` writes rather than beside
`many.vcf`, and `run_plink2.sh` makes and compares six files now instead
of five. That the script still exits 0 is what proves the generator
behind every literal of the two specs was not disturbed by the new work.

### Work package 1 as a whole

It finished as planned. The six deliverables, each with the command the
plan gives for it, run by the orchestrator on the last commit:

| Deliverable | Command | What it gave |
| --- | --- | --- |
| 1, the reference dataset | `tests/reference/ld/run_plink2.sh` into an empty directory | exit 0, which is the script saying none of the six files it made again differed from the ones in git |
| 2, one major allele | `grep -c "fn the_major_allele" crates/popnei/src/pca.rs` | `0` |
| 3, the worked example | `cargo test -p popnei --lib ld::` | the seven pairs with their six whole numbers, and their r² equal to the spec's decimals with the tolerance at 0 |
| 4, every pair against plink2 | the same | 93096 pairs with an r², 31654 NaN, none of the 93096 differing from plink2 at all |
| 5, the dosages against pyNei | the same | every one of the 25000 equal to what `to_012` gave |
| 6, the tests exist | `cargo test -p popnei --lib ld:: -- --list` | `35 tests`, where the plan started from `0 tests` |

The seven checks of the `coding` skill on the last commit: `cargo fmt
--all --check` no output; `cargo clippy --workspace --all-targets -- -D
warnings` no warning; `cargo test --workspace` `431 passed` in the core
crate and `35 passed` in the linalg crate, from 395 and 35; the same with
`--no-default-features`, which is faer instead of the system BLAS, `431
passed` and `30 passed`; `cargo wasm-check` finished; ruff `22 files
already formatted` and `All checks passed!`; `uv run maturin develop &&
uv run pytest` `242 passed`, unchanged, since nothing of this work
package reaches Python.

#### What the review found

Six reviewers read the work package: spec, tests, numbers, errors, api
and binding. Seventeen findings held and all are fixed, in eleven
commits from 98f8545 to b3b35b0. Nothing was set aside as not holding.

Three reviewers each wrote their own implementation of the formula and
compared it with the stored numbers, one of them in exact rational
arithmetic: all three got a relative difference of 0 from plink2 over the
93096 pairs, and one of them ran pyNei's `to_012` itself and matched the
stored dosages. So the agreement rests on three implementations and not
on the one under test. The counts the tests assert were recomputed from
the files rather than read from the constants: 54 variants of more than
two alleles and 257 half called genotypes in `many.vcf`, 68 variants of
no variance in `ld.vcf.gz`.

What mattered most, in the order of what it would have cost:

- **A population that names an individual twice was counted twice**, in
  n, in the major allele frequency and in every sum, and nothing said
  so: the r² of two variants of the worked example goes from 0.0625 to
  0.051470588235294115 when one individual of six is repeated.
  `of_block` refuses it now.
- **`r2_between` compared how many individuals the two sets held, not
  which**, although the message of its own error said they hold the same
  individuals in the same order. Two sets over individuals that do not
  overlap gave a number. A set of dosages now records which individuals
  it was built over and the two are compared.
- **The claim that the arithmetic is exact had no bound and nothing
  enforced one.** What has to be held exactly is not the six sums but
  the four products of the formula, so the limit is the individuals
  times the ploidy at most 94906265, which is 47453132 diploid
  individuals or 372181 at a ploidy of 255. Above it the r² loses
  digits with no word: 1.3e-12 relative at a million individuals of
  ploidy 255, wider than the tolerance the spec compares within. It is
  seven orders away from the 10000 diploid individuals popnei is built
  for. `of_block` refuses beyond it and the spec says where the
  arithmetic stops.
- **Ten allocations would have ended the process** instead of returning
  an error, where `crates/popnei/src/dists.rs` states the crate's rule
  and its reason, and where the spec asks `try_reserve_exact` of the
  matrix of `calc_r2_matrix`. A block `of_block` accepted could ask for
  17 GB in each of three vectors; under wasm, where a `usize` is 32
  bits, the byte count overflowed above 2^29 and panicked. All of them
  ask the machine now, through a case of the error called `LdNoMemory`.
- **Four of popnei's own defects would have reached Python as
  `ValueError`**, the exception a user catches for their own mistake,
  where the matching cases of the principal component analysis are a
  `RuntimeError`. Nothing in Python reached them yet, since the binding
  of this module is task 2.2, so no user could have seen it. They are a
  `RuntimeError` now and the spec says which cases are which and why.
- **One mutation of fifteen was caught by no test.** The four-product
  shortcut fires when the two sets are the same object, and every test
  of the six-product path used sets of different sizes, so replacing the
  condition with one that compares sizes left all 26 tests passing while
  giving r² values above 1. Task 2.1 of the next work package is a tiler
  that compares tile sizes on that line. There is a test now, and the
  orchestrator reran the mutation against it and saw it fail.
- **The major allele frequency was written twice**, in `ld.rs` and in
  `filters.rs`, one commit after the major allele itself was moved into
  `variant` so that it would be worked out in one place. It is in
  `variant` now and both call it.
- **plink2's matrix for the worked example was stored and read by no
  test.** The test of the worked example asserted hand-typed constants
  with nothing tying them to `example.vcf`. It reads the stored matrix
  now, so the spec's decimals, the VCF and plink2 all have to agree.

Smaller ones, each fixed: the test of every pair compared only the upper
triangle, so the diagonal and the lower half of the 500 x 500 matrix were
unchecked, and the spec asks that the 68 variants of no variance have NaN
in their row, their column and their diagonal cell; two test helpers said
their tolerance was zero and used 1e-12 and 1e-15, and the tests pass at
zero; `MAX_VALUES_OF_THE_DOSAGES` was a third copy of a limit the linalg
crate kept private; the error list of `of_block` omitted the cases the
counting of alleles refuses; a genotype that was not in the row left the
previous variant's genotype behind rather than a missing one; and a float
reached a `u8` with no check for NaN or range.

#### What the owner should know

The sentence this report's orchestrator added to the spec on 22 September
2026, that no reader of popnei gives a block of more than 255 alleles in
a genotype, was wrong. It holds for the VCF reader, which caps at 255,
and not for the vars reader, which takes the ploidy from the `popnei` key
of the file with no upper bound. So four of the cases this module refuses
are reachable from a vars file, and the spec says so now. Whether the
vars reader should cap the ploidy itself belongs to
`docs/specs/io_vars.md` and not to this plan.

The question the plan's "What could go wrong" raised has come back with a
measurement, and it is the one thing waiting on the owner.
`crates/popnei-linalg` has no product that takes its second matrix
transposed, and the three matrices are variants x individuals, so
`r2_between` transposes its second matrix itself. The `coding` skill says
that linear algebra goes through the linalg crate and nowhere else and
that no code of the core crate does its own, which the transpose in
`ld.rs` is. A transpose of a 512-variant by 1000-individual matrix takes
0.422 ms and an off-diagonal tile pair does three of them, so the matrix
of 5000 variants of work package 4 would spend about 69 ms of its 0.50 s
target copying, a seventh of it; transposing each tile once instead of
once per pair brings that to about 13 ms without touching linalg. Adding
the product to linalg costs one flag in the `dgemm` call of `blas.rs` and
a `.transpose()` on the matrix reference of `faer.rs`, both of which
those libraries do internally, and it changes `docs/specs/linalg.md`, a
spec this plan does not build.

One finding was left for the owner rather than acted on: `has_variance`
and `maf` answer `false` and `None` both for a variant that has no data
and for an index that is not a variant of the set at all, so a tile or a
window that runs one past its end writes a NaN row or drops a variant
with no error. The spec fixes those two signatures, `-> bool` and
`-> Option<f64>`, so narrowing them is the owner's.

#### How the work went

The five tasks cost their subagents 102k, 97k, 206k, 181k and 189k
tokens. The two that were only a move, 1.1 and 1.2, cost half of what
the three that wrote the module cost. The six reviewers cost 144k, 121k,
134k, 163k, 135k and 92k, and the subagent that fixed the seventeen
findings ran to 324k in all, having written `r2_between` first.

Nothing had to be sent back to a subagent as wrong. The orchestrator
checked every claim that the next step rested on by running it, and twice
found the claim understated: the r² of the worked example and of the
whole dataset is not within the tolerance but equal to plink2's bits,
which was confirmed by setting the tolerance to 0 and rerunning.

Two things about the orchestration are worth the next plan's attention.
The orchestrator left task 1.5 unticked in the plan while it checked the
deliverables, and it was a reviewer that noticed. And a figure from a
subagent's hand-back, the size of the transposes, went into this report
without being recomputed and was wrong: two of the four transposes of a
set against itself are 4.1 MB and not 2.1 MB. The `following-plans` skill
says to check what a subagent claims when the next step rests on it; no
step rested on either of these, and both were wrong in the record.
