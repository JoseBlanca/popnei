# Work report: r², the matrix of it, and the filter by linkage disequilibrium

The plan `docs/plans/ld.md` built three of its four work packages, on the
branch `plan/ld`. The fourth, which measures the speed of what the other
three built, was deferred on 23 September 2026 by the owner's decision:
it becomes a performance review of its own rather than the last work
package of this plan, because the review of the three that were built
left a list of candidates with measurements already attached to them and
popnei has no bench for either calculation to measure against. So popnei
can now work out every number this plan set out to give, and nobody has
yet measured how long it takes.

Where the merge stands. `main` moved 110 commits while this plan ran, and
the two sides changed 26 files in common, so the merge is a piece of work
and not a formality: `main` now builds the chain of a pass from a list of
steps where this branch built the filter by linkage disequilibrium into
the list of criteria it built from before. The branch is being merged as
this project has merged its others, `main` into the branch first, so that
`main` is never left unable to build. Until that is done the three work
packages are on `plan/ld` and nowhere else.

## What exists now that did not

A user writes `calc_rogers_huff_r2_matrix(variants)` in Python or
`calcRogersHuffR2Matrix` in TypeScript and gets the r² of every pair of
the variants of a pass, the squared correlation between their dosages,
with the chromosome and the position of each variant beside it and the
counts of every filter the pass ran. A user writes
`variants.filter_by_ld(max_allowed_r2, max_dist)` or `filterByLd` and
every consumer of those variants afterwards sees only the ones that do
not repeat what a variant kept within `max_dist` base pairs of them on
their chromosome already said. Under both is
`crates/popnei/src/ld.rs`, which works r² out for two sets of variants
through six matrix products.

`crates/popnei-linalg` gained a fourth operation on the way, by the
owner's decision of 23 September: the product of a matrix with the
transpose of another, which the r² needs and which the crate did not
have.

## Whether the numbers are right

popnei's r² is plink2's **to the bit**. Not within a tolerance: equal.
That holds for all 93096 pairs of the reference dataset that plink2 gives
a number for, for all 250000 cells of the matrix, and on the three
arithmetic paths popnei is built for — the BLAS of the machine natively,
the faer library of Rust natively, and faer compiled to WebAssembly and
run under node. It was checked by setting the tests' tolerance to exactly
0 and rerunning, which the orchestrator did at every stage where the
arithmetic underneath was changed.

The dosages popnei counts are pyNei's `to_012` exactly, over all 25000
genotypes of a file with 54 variants of more than two alleles and 257
half called genotypes. Where popnei and pyNei part is where they were
always going to: popnei leaves the individuals missing at either variant
of a pair out of that pair and pyNei leaves them in with a dosage of -1.
The test pins that divergence at the median 0.00369, 99th percentile
0.04680 and largest 0.19374 of the table of `docs/specs/ld.md`, so a
change on either side of it shows.

The filter keeps 84, 133, 85 and 85 variants of 500 of the reference
dataset, at the four settings of the table of `docs/specs/filters.md`: a
window of 10000 base pairs at a threshold of 0.1, then 10000 at 0.3, then
50000 at 0.3, and then 250000 at 0.3, which is a whole chromosome of that
dataset. Four implementations that share no arithmetic agree on that set:
the rule of the spec written again in Python over plink2's stored matrix,
popnei's reader in Rust, the same compiled to WebAssembly, and the Python
binding. A reviewer wrote a fifth from the spec's words alone, worked out
for itself which variants have two dosages rather than reading plink2's
diagonal, and got the same set variant by variant.

## What the review found

Ten reviewers read the three work packages. Thirty-three findings held,
two were refused with evidence the orchestrator accepted, and **not one
of them was a wrong number**: nothing wrong ever reached a user, because
none of this had been released. What they found was that the code asked
too much of a machine and that some of its checks could not fail. The
sections below the rule have each of them with its evidence; the two
worth the owner's time are these.

The filter was doing about eighteen times the work its spec asks,
comparing the variants it had kept one at a time where the spec says one
set of matrix products. Rewriting it as the spec has it took a pass of
20000 variants of 400 individuals from 2.62 s to 0.180 s and moved no
variant of the result.

The interface this report's orchestrator wrote for the new operation of
the linear algebra crate was a trap. It gave two functions that took the
same six arguments of the same types and differed only in which way round
the second matrix is read, so each would take the other's call and return
a different matrix with no error. A reviewer showed it on the r²'s own
shape. There is one function now, whose second operand carries its layout
in its type, so the wrong call cannot be written by accident.

## What is open

- **The speed of both calculations is not measured.** This is the whole
  of work package 4 and it is now a performance review. Nothing in this
  plan says how long the matrix of 5000 variants takes against the 0.50 s
  `docs/specs/ld.md` asks of it, or what a pass of the filter costs.
  popnei has no bench for either.
- **Two methods of `LdDosages`, the type that holds the dosages of a set
  of variants, answer the same for two different things.**
  `has_variance` says whether a variant has two dosages at least and
  `maf` gives its major allele frequency; each of them answers `false`
  and `None` for a variant that is there and has no data, and the same
  `false` and `None` for a number that is not a variant of the set at
  all. So a caller that runs one past the end of a set is told the
  variant has nothing rather than that it asked for nothing. Narrowing
  them means changing their signatures, which `docs/specs/ld.md` lays
  down, so it is a change to the spec and the owner's to make.
- **Which exception a machine that cannot give memory should raise.**
  `.claude/skills/coding/SKILL.md` records the convention the owner gave
  on 21 September 2026, with a `ValueError` for a wrong input, a
  `RuntimeError` for a defect of popnei and an `OSError` for a file. A
  machine too small for a matrix is none of the three, and Python has
  `MemoryError` for it. Two cases sit on this, the matrix of the r² and
  the Kosman distances of `docs/specs/dists.md`, which are the other
  calculation of popnei that asks for memory the size of its dataset.
  They agree today by both being a `ValueError`.
- **Both bindings copy the matrix on the way out**, 200 MB at the cap of
  5000 variants that `calc_rogers_huff_r2_matrix` takes by default,
  because the type that holds the matrix lent its values rather than
  giving them away. It can now give them away, and neither binding has
  been changed to take them; in a browser that copy is 200 MB of the
  instance's memory that is never given back.
- **The second item of `docs/specs/ld.md`**, the curve of r² against
  distance per population, is written and reviewed and not built. The
  owner decided that on 22 September: nothing calls it today, where the
  filter is what pop_lab uses to prune before a principal component
  analysis, and it is the most machinery of the three. The first of the
  two questions that spec leaves open belongs to it, so it needs no
  answer until that item is built.
- **The second of those two questions is unanswered and this plan leans
  on it**: how the major allele of a variant is chosen when some of its
  genotypes have one allele called and one missing. Everything here
  follows the rule `docs/specs/pca.md` already gives, which the spec
  names as what to do until the question is settled. Answering it the
  other way would move the dosages, and so the r², of variants of more
  than two alleles.

## What a performance review should start from

Each of these was measured by a reviewer of this plan, on the owner's
Apple M5 Pro, and each is a candidate and not a conclusion.

**Start with what a pass of the filter costs on the dataset popnei is
sized for**, 100000 variants of 1000 individuals, which is the one
`docs/rust_core.md` measures popnei on. It is the largest hole: the only
figure there is comes from the fix above, 20000 variants of 400
individuals with a window holding about 270 kept variants, 2.62 s before
and 0.180 s after. How much a window holds is set by the dataset and not
by the filter, so the number means nothing without saying what the window
held, and nothing is known at the scale that matters.

**Then the size of a tile**, the square block of variants the matrix is
worked out in, which is 256 today. It was chosen by reading the "Speed"
table of `docs/specs/ld.md` as 1.9 ms for a pair of 256-variant tiles
against 5.7 ms for a pair of 512-variant ones, where a pair of 512 covers
four times as many pairs of variants. Per pair of variants a reviewer
measured 51.8 ns at 128, 38.0 ns at 256, 30.1 ns at 512 and 29.1 ns at
1024, which at the cap of 5000 variants is about 0.43 s of products at
256 against 0.35 s at 512, where the number to reach is 0.50 s. That is
the one candidate with a pass or fail already attached to it, and the
comment in the code says the size is what a performance review settles.

Then, in no order:

- **The copy each binding makes of the matrix**, 200 MB at the cap, which
  the core can now hand over instead of lending.
- **The transposed half of each pair of tiles**, written one cell at a
  time with the width of the whole matrix between one write and the next:
  about 800 MB of scattered writes at the cap. A reviewer asked that it
  be measured apart from the products before anything is done to the
  products.
- **`LdDosages::rows`**, which gives an owned copy of the dosages it is
  asked for, and the six sums a pair of tiles allocates inside the call
  with no way to hand in a buffer to reuse, where `pca.rs` has a type for
  exactly that.
- **The memory of the filter's window**, which holds the genotypes of
  every variant it has kept: one byte for each allele between blocks, and
  24 bytes for each individual and each of its variants while a set is
  settled. A dataset whose variants carry no linkage leaves every one of
  them in a window as wide as a chromosome.

## What is asked of the owner

Two decisions, both above: whether `has_variance` and `maf` should be
narrowed, which changes `docs/specs/ld.md`; and whether the convention of
exceptions should gain a fourth for a machine that has not the memory,
which changes what every user of popnei catches. The performance review
is the other thing this plan hands on.
---

The rest of this file is what was written while the work went, work
package by work package, and it is where every number above comes from.

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

`docs/specs/ld.md` ends with two questions it leaves for the owner, which
it numbers 1 and 2. The second asks how the major allele of a variant is
chosen when some of its genotypes have one allele called and one missing,
and it has no answer yet. Each such question carries a rule to follow
until it is answered, which the spec calls its "meanwhile", and this
one's is the rule `docs/specs/pca.md` already gives and `pca.rs` already
computes. So task 1.2, which moved the function that picks the major
allele out of `pca.rs` and into `variant.rs`, changed no number. The
first question belongs to the curve of r² against distance, which this
plan does not build.

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

The function that reads a block of variants and builds the three
matrices from it is called `of_block`, and it is where a dataset popnei
cannot work with is refused. The task found a fifth way it can refuse a
block that the spec's list of four did not have: genotypes of more than 255 alleles each. A
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
would reach Python as `ValueError`, where the case of the principal
component analysis that stands for the same thing, a product its linear
algebra could not work out, is a `RuntimeError`, and `docs/specs/pca.md`
says why: no argument of the function can give it. The same holds for
the case of task 1.3 that refuses variants asked of a set of dosages
that does not hold them. Nothing in Python
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
| 6, the tests exist | `cargo test -p popnei --lib ld:: -- --list` | `36 tests` once task 1.6 was in, where the plan started from `0 tests` |

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
  ploidy 255, wider than the tolerance the spec compares within. The
  dataset popnei is built for, 10000 diploid individuals, gives an
  individuals times ploidy of 20000, so the limit is 4700 times further
  out than anything popnei expects to be given. `of_block` refuses
  beyond it and the spec says where the arithmetic stops.
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
  giving r² values above 1. Task 2.1 of the next work package works the
  matrix out in square tiles, so that the six matrices it holds at once
  stay small, and it will be comparing the sizes of those tiles on that
  very line. There is a test now, and the
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
zero; the largest number of values a matrix may hold, 2147483647, which is
what the routines of BLAS count them in, was written out a third time
where the linalg crate keeps it private; the error list of `of_block` omitted the cases the
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

One decision was waiting on the owner and was taken on 23 September 2026:
the operation went into `crates/popnei-linalg`. Task 1.6 below is what
carried it out. What was decided and why is kept here, since the plan as
approved did not have that task in it.

It was the risk the plan's "What could go wrong" named, come back with a
measurement.

What is being decided: whether `crates/popnei-linalg` gains a product
that takes its second matrix transposed, which changes
`docs/specs/linalg.md`, a spec this plan does not build, or whether the
r² keeps transposing its matrices itself inside
`crates/popnei/src/ld.rs`.

The recommendation was to add it to `crates/popnei-linalg`, and to do it
before work package 2 works the matrix out in tiles, and that is what the
owner chose. The orchestrator
decided the other way when the problem appeared, on the grounds that it
touched no spec outside the plan, and the review showed that reasoning
incomplete: the `coding` skill says that linear algebra goes through the
linalg crate and nowhere else, and that no code of the core crate does
its own. The transpose written in `ld.rs` is code of the core crate
doing its own linear algebra, so the decision kept a plan boundary at
the cost of a rule the project holds.

Why the shape came up at all: `crates/popnei-linalg` multiplies two
matrices that are both stored row after row, and the three matrices of
the r² hold one row for each variant and one column for each individual.
The product of one set of variants with another sums over the
individuals, so the second matrix has to be handed over the other way
round, one row for each individual, and something has to turn it round.

What each way costs. Turning a matrix of 512 variants and 1000
individuals round takes 0.422 ms, and a pair of tiles that are not on the
diagonal needs three of them, so the matrix of 5000 variants that work
package 4 must bring in under 0.50 s would spend about 69 ms of that
budget copying, a seventh of it. Turning each tile round once and reusing
it, rather than once for every pair it appears in, brings that to about
13 ms and touches no other spec. Giving the product to linalg removes the
copying altogether, and costs almost nothing to write: one flag in the
call to `dgemm`, the routine of BLAS that multiplies two matrices, which
already takes a flag saying that an operand is stored the other way
round, and one call to `.transpose()` on the faer side, which is the
other library `crates/popnei-linalg` is built on. What it costs is the
spec: `docs/specs/linalg.md` gains a function, and this plan grows by a
task that was not in it.

What followed. The orchestrator wrote the spec change, task 1.6 carried
it out and was itself reviewed, and work package 2 builds on it: a tile
is a range of rows of a matrix held row after row, so it is a piece of
that matrix with nothing copied, and both operands of every product of a
tile pair are read with one row for each variant. The transposes that
this work package would otherwise have paid for once per pair of tiles do
not happen at all.

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

### Task 1.6, added after the review: the product with a transposed operand

The owner decided on 23 September 2026 to put the operation into
`crates/popnei-linalg` rather than leave the r² transposing its own
matrices, which is what the first review of this work package had found
the core crate doing against the `coding` skill. The plan gained task 1.6
for it, and `docs/specs/linalg.md` an operation, written before the code.

What it gives. The linear algebra crate had three operations and now has
four: the product of a matrix with the transpose of another, where both
matrices hold one row for each of the things they describe and the
product sums over their columns. It is the same `dgemm` of BLAS and the
same `matmul` of faer, told that the second operand is read the other way
round, which both do inside the routine. The claim that neither pays a
copy for it was not taken on trust: 1.174 ms against 1.179 ms for 512 x
1000 on Accelerate on one thread, and within 0.3 % on faer. Had it been
false the operation would have been worth nothing.

The number that had to stay still, stayed still. With the tolerance of
the r² tests set to exactly 0.0 by hand, the whole of `ld::` passes on
the system BLAS and on faer alike, before and after the rewiring and
again after the interface changed: every one of the 93096 pairs of
`ld.vcf.gz` that plink2 gives a number for is still equal to plink2's
bits. `cargo test --workspace` gives `432 passed` in the core crate and
`42 passed` in the linalg crate, from 431 and 40, and `37 passed` there
on faer, from 35.

What it saves, measured end to end: one pair of tiles of 512 variants of
1000 individuals went from 8.25 ms to 7.40 ms, and a set of variants
against itself from 5.97 ms to 5.02 ms.

#### What that review found

Three reviewers read it: numbers, tests and architecture. Eight findings
held and are fixed in commits 19a3e98 and 0804219; one was refused with
evidence, below.

The serious one was the orchestrator's and not the implementer's. The
spec had given the new operation as a function of its own,
`product_by_transpose`, taking the same six arguments of the same types
as `product`, the two differing only in which way round the second matrix
is read. A reviewer showed that each would take the other's call and
return a different matrix with no error, since the check of a length is
the rows times the columns either way: on the shape of the r² itself,
`product` given the arguments meant for the other returned
`[1, 1, 6, 3]` where the answer is `[3, 4, 1, 5]`. That is a wrong number
with no word, reachable by a plausible slip. There is one product now,
whose second operand carries its layout in its type, so a caller cannot
reach for the wrong one without writing the name of the wrong one. The
spec was rewritten to match, and the three call sites of the principal
component analysis changed with it.

The one that will matter later. The two ways the r² works out the sums of
a pair, the four products when a set is against itself and the six
otherwise, are equal to the bit ONLY because the six sums are whole
numbers. The shortcut reads Σy of the pair i, j as Σx of the pair j, i,
which is the same dot product with the operands in the other roles, and a
routine need not sum those in the same order. A reviewer measured the two
orders differing at 1e-17 on Accelerate with values that are not whole, 4
entries of 25 at one shape and 19 of 49 at another. Nothing is wrong
today and the bound of the individuals times the ploidy is what
guarantees it, but the doc comments stated the identity as if it held for
any float. The next person to reuse this on values that are not whole,
the genomic relationship matrix of `docs/specs/kinship.md` or a GWAS,
would have got a last-bit difference against plink2 that nobody would
trace. It is said where the identity is stated now.

The rest: four claims of the new operation's doc comment that no test
guarded, where deleting the check of the first operand's length, or of
its values being finite, or mislabelling it, left all 40 tests passing;
the same two holes in the older `product`, closed while they were there;
a test that builds the two ways of holding the second set's sums over the
same variants and compares the two r² matrices, which is what would catch
the two paths drifting apart; a loop that would have truncated in silence
if its invariant ever broke, now asserted; and three doc comments and a
module doc left over from the transpose that no longer exists.

One finding was refused, with evidence the orchestrator accepted. The
`rows == 0` early return of `product` is guarded by no test, and cannot
be: the implementer deleted it and both backends still wrote nothing,
because `dgemm` and faer's `matmul` are a no-op for a product of no rows.
What the guard buys is that no routine is reached with a dimension of 0
and the pointer of an empty slice, which is now said in a comment instead
of asserted in a test.

#### A measurement that would not settle

The first review said that writing a transpose out costs 0.422 ms for 512
x 1000, and this report repeated it, and `docs/specs/linalg.md` carried
it. It was then timed three more times on the same machine, once by the
orchestrator: 0.135 ms, 0.253 ms and 0.467 ms. Four measurements of one
operation spread over a factor of 3.5, depending on whether the
allocation was counted, how warm the cache was and what the compiler kept.
A number that unstable is not a fact a spec can hold, so the spec no
longer quotes one. What it holds instead is the size of the copy, 4.1 MB
for 512 x 1000, and what was measured end to end, the tile pair going
from 8.25 ms to 7.40 ms. The decision never rested on the figure. The
lesson for the next plan is the one already in this report: a number that
arrives in a hand-back and is not recomputed gets into a document, and
from there into a spec.

## Work package 2: the matrix of every pair

Task 2.1, `calc_r2_matrix` and `R2Matrix`, commit 0cae3ae. The core works
out the r² of every pair of the variants of a pass and gives it with the
chromosome and the position of each variant beside it. `cargo test -p
popnei --lib ld::` goes from 36 tests to 47, and with the tolerance of
the r² tests set to exactly 0.0 by hand the whole of `ld::` passes on the
system BLAS and on faer alike, which the orchestrator ran itself: the
250000 cells of the matrix of `ld.vcf.gz` are plink2's to the bit, and
not merely within the tolerance.

Three decisions the task took that the plan left open.

The tiles are cut at 256 variants and not 512, which the "Speed" section
of `docs/specs/ld.md` supports: one tile pair of 1000 individuals takes
1.9 ms at 256 variants and 5.7 ms at 512, and a pass over 100000 variants
1.5 s against 2.2 s. That table was measured for the curve of r² against
distance, the item this plan does not build, so it is guidance and not a
measurement of this function; work package 4 measures this one. At 256 a
tile also covers the 500 variants of the reference dataset in two tiles
and a short one, so the tests cross a tile boundary.

A tile is its own set of dosages, built from the genotypes of its
variants gathered out of the blocks, because a pass arrives as many
blocks while the tiles are cut at fixed multiples of 256 from the first
variant of the pass. So no tile copies a matrix of floats, `r2_between`
did not change, and a tile on the diagonal passes one reference twice,
which is what gives it the four products instead of six.

The memory of the matrix is asked for with `try_reserve_exact` after the
pass and not before it, since how many variants a pass gives is not known
until it ends. What is refused before the source is read is a
`max_num_vars` whose square is not a number this machine counts.

Nothing of this task uses rayon.

## Work package 3: the filter by linkage disequilibrium

Task 3.1, `LdFilter`, commit c5d9e83, built at the same time as task 2.1
and in other files. `cargo test -p popnei --lib filters:: -- --list` goes
from 34 tests to 51, and `cargo test --workspace` gives `460 passed` in
the core crate with the two tasks together.

The plan calls this the task whose failure would be silent, a set of
variants that is wrong rather than a crash, and the first sign is good:
run by hand over `ld.vcf.gz`, the filter keeps 84, 133, 85 and 85
variants of 500 for the four rows of the table of "How it is verified",
with the five positions each row gives and chr2:1000 first of its
chromosome, and it keeps the same variants in blocks of 7, of 64 and of
the whole file. Those are the numbers deliverable 2 asks for. They are
not in the tests yet: task 3.3 is what commits them, and until it does
this is a measurement and not a check.

Two decisions the spec left to the implementer.

The filter settles 256 variants at a time. Each variant of the window is
one comparison against the whole set of 256, the set against itself is
one more, and a candidate is then read against the variants kept inside
its own set out of that matrix, in order. 256 is the shape
`docs/specs/ld.md` measures at 1.9 ms, and the matrix of a set is 512 KB
whatever the blocks the source gives.

The window holds one set of dosages for each variant it has kept, the
one row of the dosages of the block that variant arrived in, with its
chromosome and its position: 24 bytes for each individual, which is what
the spec says. It cannot hold one set of dosages for the whole window,
because a set is built from one block and popnei has no way to join
variants of several blocks into one. Doing that would need something on
`LdDosages` that does not exist, and the task did not bolt one on to
`filters.rs`. Whether the window should hold one set instead of one per
variant is a question for work package 4, which measures it.

### How the two ran side by side

They touched different files, as the plan says, and both staged by
explicit path, so neither swallowed the other's work. One thing did go
wrong. `crates/popnei/src/error.rs` is a file both needed, and the
subagent of task 2.1 committed it whole while the two error cases of the
filter were in it, so those two are in the commit of the matrix and not
in the commit of the filter. Nothing is lost and everything passes; the
history is what suffers, and a reader looking for where the filter's
refusals came from will find them in a commit about something else. The
instruction the orchestrator gave, to commit only one's own hunks of a
shared file, asks for something git makes awkward. The next plan should
give each file one owner for the length of a work package, or not run
two tasks that need the same file at once.

Task 2.3, the TypeScript function, commit 940682c. `calcRogersHuffR2Matrix`
and its result in `js/popnei/src/ld.ts`, with the binding in
`crates/popnei-js`. `npm test` in `js/popnei` goes from 172 tests to 180.

The five r² of the spec's table come out under node, where the products
run on faer and not on the system BLAS, **equal to the spec's decimals to
the bit**: 0.353466669239891 is the 64 bits `3fd69f32aa2720de` on both
sides, which the orchestrator checked, and so are the other four. The
committed test keeps the 1e-12 the spec gives.

One thing the owner should know about the memory. The binding copies the
matrix once, because the core lends its values and gives up no vector, so
at the cap of 5000 variants the core's 200 MB and the copy's sit in the
memory of the WebAssembly instance together. The copy asks with
`try_reserve_exact`, so a tab without the room gets an `Error` and not a
dead instance. An accessor on `R2Matrix` that gives the vector up would
remove the copy; the task did not add one, since `crates/popnei/src/ld.rs`
was not its file.

Task 3.2, `LdFilteredReader` and the chain, commit 1d41f14. The filter is
a reader over a reader now, with the three rules of `docs/specs/block.md`,
`MaxLdR2` and `max_dist()` on the criterion, the kind `"ld"`, its place in
`chain_of`, and the refusal of a second filter of the kind. `filters::`
goes from 51 tests to 63 and the workspace from 460 to 472.

The four counts of the filter did not move: the task ran them again
through the new reader, at the four settings and at three block sizes,
and got 84, 133, 85 and 85 of 500 with the five positions of each row.

It added one case to the error of the crate that the spec did not have,
and asked whether it belonged there. It does, and the spec has it now,
commit 0b5dd76: the plain filter of a threshold cannot answer for the
criterion of this filter, because whether a variant is kept turns on the
variants kept before it and not on the variant alone, so building one for
that criterion is refused. Nothing a user writes reaches it, since
`chain_of` sends that criterion to the reader that does answer it, so it
is a `RuntimeError` and not a `ValueError`. The spec said the module adds
three cases and it adds five; the other one it did not name is the
variant whose position does not rise within its chromosome.

### What running three tasks at once cost

Giving every file one owner stopped the trouble of the first pair: no
task committed another's work. A different thing happened instead. Adding
`MaxLdR2` to the criterion of the core crate broke both binding crates,
which match on that type and have no arm to spare, so the workspace would
not compile until each binding crate gained one. The task that added the
variant did not own either binding crate, and the two tasks that did own
them found their crate broken by work they had not done. They fixed their
own side and said so, and nothing was lost, but the orchestrator had not
foreseen it: ownership of a file is not the same as ownership of what
compiles. A task that adds a variant to an enum the bindings match on
should either own the arms it breaks or not run beside the tasks that do.

Task 2.2, the Python function, commit 732e566. `calc_rogers_huff_r2_matrix`
reaches Python with its result and the counts of the pass, and
`tests/test_ld.py` holds 15 tests where the file did not exist. `uv run
pytest` goes from 242 to 257.

This is where popnei is first compared with pyNei for the r², everything
before it having been checked against plink2's stored matrix. Both halves
of deliverable 3 came out.

With the missing genotypes taken out, `filter_by_missing_data(0)` on
`ld.vcf.gz` keeps 21 variants of 500, and popnei's matrix and pyNei's r
squared agree within 1e-12 relative over the 256 cells of those that are
numbers, the other 185 being NaN on both sides.

With the missing genotypes left in, the two libraries part, and the test
pins by how much: over the 1.4 million pairs of the panel the difference
is median 0.00369, 99th percentile 0.04680 and largest 0.19374, which is
the table of "Missing genotypes" of `docs/specs/ld.md` to its digits.
popnei leaves the individuals missing at either variant of a pair out of
it and pyNei leaves them in with a dosage of -1.

The task changed which dataset that second half is measured on, and was
right to. The plan puts both halves on `ld.vcf.gz`, but the table whose
numbers pin the divergence was measured on the panel, so on the panel is
where those numbers mean anything; on `ld.vcf.gz` the test asserts only
that the two libraries differ. The task measured that difference there
too, over the 186192 pairs where both libraries have an r²: median
0.00578, 99th percentile 0.18071, largest 0.68299. Those are not in the
spec and so are not asserted. What is worth the owner's eye in them is
the count of cells with no r² at all: popnei has 63376 and pyNei 4975,
because a dosage of -1 gives variance to variants that have none.

### A table of the spec that its own program does not produce

The task, checking which of r and r² the table of "Missing genotypes" is
in, found that `docs/reports/ld-method/missing_rules.py` compares plink2's
r, while the table is in r². The orchestrator ran the program to be sure.
It gives, in r, a median of 0.0300, a 99th percentile of 0.1245 and a
largest of 0.6002 for pyNei's rule, where the table says 0.0037, 0.047
and 0.194. The README of that directory says the program "is the table of
'Missing genotypes' of `docs/specs/ld.md`", and it is not.

The table is the one that is right: two independent computations in r²
reach its digits, and the median r² of 0.0071 it quotes is the square of
the 0.0845 median r the program prints. So the numbers a reader would act
on are sound and the program kept to reproduce them is not, which is the
worse way round for a number that has to survive a change on either side.
It is being put right.

Task 3.3, the checks against the stored numbers, commit 921114e. This is
the guard the plan asks for over task 3.1, the task whose failure would
have been a wrong set of variants and not a crash. It is a check now and
no longer a measurement.

The four counts hold, and they hold twice over by two routes that share
no arithmetic: the rule of the spec written again in Python over plink2's
stored r² matrix, and popnei's own reader. Both give 84, 133, 85 and 85
variants kept of 500 at the four settings of the table, the five
positions of each row, the same variants in blocks of 7, of 64 and of the
default size, and chr2:1000 first of its chromosome. The cargo test
asserts that popnei's kept variants are exactly the ones the stored file
holds, chromosome and position, which is more than deliverable 2 asked
for: it asked for the counts and five positions.

The three properties of the kept set hold at all four settings, with
nothing violating any of them: no kept variant has fewer than two
dosages, no kept pair inside a window is above the threshold, of 84, 288,
593 and 1764 such pairs, and no dropped variant with two dosages lacks
something above the threshold kept before it in its window, of 348, 299,
347 and 347 such variants. Those counts of what was examined are
consistent on their face: 500 less the 84 kept is 416, less the 68
variants of one dosage is the 348 the third property looked at.

All three were shown to be able to fail: pyNei's rule leaves 707 pairs
above the threshold, a set that kept everything leaves 68 variants of one
dosage, and a set with ten variants taken out of it leaves 45 dropped for
no reason.

`cargo test --workspace` goes to `473 passed` and `filters::` to 64
tests. `run_plink2.sh` writes eight files now and still exits 0, so
`ld.vcf.gz` is what it was.

### A claim of the spec that this disproved

`docs/specs/filters.md` said the three properties do not pin the kept set
and that "a set that dropped a variant it could have kept would pass all
three". It would not. Such a variant has two dosages and nothing above
the threshold kept before it inside its window, which is exactly the
negation of what the third property asks of every dropped variant, and
the deliberate break above named 45 of them. The spec now says what each
of the three catches, and the counts of the table are a check beside them
rather than the only one. Commit 6d127ad.

Task 3.5, the TypeScript step, commit 098a251. `variants.filterByLd` is
there with its binding, and `npm test` in `js/popnei` goes from 180 tests
to 200.

Its test asserts all four rows of the table and not only the one
deliverable 6 asks for: 84, 133, 85 and 85 kept of 500, the five
positions of each row, chr2:1000 first of its chromosome, and blocks of 7
keeping the same 133. Every one matched on the first run with no
expectation adjusted, which is the third independent route to those
counts after the Python rule over plink2's matrix and popnei's own
reader. It also asserts the step with both its arguments, the refusals of
seven bad values of `maxAllowedR2` and four of `maxDist`, a second filter
of the kind, and a call after the memory is given back.

One decision: `maxDist` crosses as a `u32`, and the package refuses what
is not a whole number of 1 or more before it crosses, as it already does
for the variants of a block, rather than letting the core say it. A
negative number would otherwise arrive as 4294967295.

Task 3.4, the Python step, commit 17c76ad. `Variants.filter_by_ld` is
there with its binding, and `uv run pytest` goes from 257 tests to 276,
19 of them this task's.

The four counts and the eight positions matched on the first run with no
expectation adjusted, which is the fourth route to them. The `filtering`
of a pass with a maf filter before the ld one reads `maf` 500 of 500 then
`ld` 84 of 500.

The method is on `Variants` in `python/popnei/variant.py` and not in
`python/popnei/filters.py`, where the plan puts it: `Variants` lives in
`variant.py` and `filters.py` holds `Step` and `FilteringStats` alone. The
plan named the wrong file.

The step's arguments stopped being two floats. A window is a whole number
of base pairs, and a user who wrote 10000 was getting `10000.0` back in
the `args` of their `Step`; the arguments now carry a threshold or a
distance and the distance stays whole.

### The table of the spec regenerates again

The repair of `docs/reports/ld-method/missing_rules.py`, commit 084dd4c.
The program reads plink2's r² for the panel now, works the first rule out
from the six whole numbers itself, and squares what pyNei returns, since
`_calc_rogers_huff_r2` gives r whatever its name says. The orchestrator
ran it from an empty directory with the plink2 command the README now
gives:

| the rule | what the program prints | what the table says |
|---|---|---|
| the individual is left out of the pair | 0, 0, 0 | 0, 0, 0 |
| the genotype takes the mean dosage of its variant | 0.000354, 0.00782, 0.0422 | 0.00035, 0.0078, 0.042 |
| pyNei: the genotype is a dosage of -1 | 0.00369, 0.0468, 0.194 | 0.0037, 0.047, 0.194 |

and "plink2 r2 of the panel: median 0.0071, NaN in 0 of the 1438800
pairs", which is the 0.0071 the text beside the table quotes. The first
row is 0 over all 1438800 pairs and not a rounding, as the spec claims.
So the table was right all along and is now reproducible from the program
kept beside it.

### Work packages 2 and 3 as a whole

Both finished as planned, and two of their deliverables came out stronger
than the plan asked. The orchestrator ran every check below itself.

Work package 2, the matrix of every pair:

| Deliverable | What it gave |
| --- | --- |
| 1, `calc_r2_matrix` in the core | the five pairs of the table with their `n`; all 250000 cells of `ld.vcf.gz` against plink2's matrix; the 68 variants of one dosage NaN in their row, their column and their diagonal cell and the other 432 at exactly 1 against themselves; the same matrix to the bit for every one of the sixteen pairings of blocks and tiles of 7, 64, 256 and 500 variants, where the plan asked for four block sizes; and the same matrix in a rayon pool of 1 thread and of 4 |
| 2, the Python function | `uv run pytest tests/test_ld.py` 15 passed, where the file did not exist |
| 3, the comparison with pyNei | within 1e-12 with the missing genotypes out; with them in, the median 0.00369, 99th percentile 0.04680 and largest 0.19374 of the spec's table |
| 4, the refusals | a cap below the variants of the pass, a cap the machine cannot count, a source with no variant, and steps that kept none |
| 5, the TypeScript function | `npm test` in `js/popnei` 202 passing, where `js/popnei/src/ld.ts` did not exist |

Work package 3, the filter by linkage disequilibrium:

| Deliverable | What it gave |
| --- | --- |
| 1, `LdFilter` and `LdFilteredReader` | both in `filters.rs`, with the window, the rule, the counts and the refusal of a position that does not rise |
| 2, the four rows of the table | a cargo test asserting popnei's kept variants are exactly the stored set, chromosome and position, at the four settings and at blocks of 7, of 64 and of the default size, where the plan asked for the counts and five positions |
| 3, the three properties | stored in `tests/reference/ld/ld.filter.properties.txt`, nothing violating any of the three at any of the four settings, and each shown able to fail |
| 4, the criterion and the chain | `cargo test -p popnei --lib filters:: -- --list` 64 tests, from the 34 the plan started with |
| 5, the Python step | every test the spec names, including the source whose positions go backwards; `uv run pytest` 276 passed, from 242 |
| 6, the TypeScript step | all four rows of the table, where the deliverable asked for one |

The seven checks of the `coding` skill on the last commit: `cargo fmt
--all --check` no output; `cargo clippy --workspace --all-targets -- -D
warnings` no warning; `cargo test --workspace` `473 passed` in the core
crate and `42 passed` in the linalg crate; the same with
`--no-default-features` `473 passed` and `37 passed`; `cargo wasm-check`
finished; ruff `All checks passed!`; `uv run maturin develop && uv run
pytest` `276 passed`; and `npm run build && npm test` in `js/popnei` 202
passing. `tests/reference/ld/run_plink2.sh` into an empty directory exits
0, so the dataset every literal rests on is unchanged.

Four routes now reach the filter's four counts, and they share no
arithmetic: the rule of the spec written again in Python over plink2's
stored r², popnei's reader in Rust, the same compiled to WebAssembly and
run under node, and the Python binding. All four give 84, 133, 85 and 85
of 500 with the same five positions each.

Two things were put right that no deliverable asked for.
`js/popnei/README.md` named three consumers of a `Variants` and neither
`doPca`, which the previous plan added, nor `calcRogersHuffR2Matrix`,
which this one did; it now names all five with what each gives. And the
two languages took different windows: Python up to the 2^64 - 1 the spec
fixes, TypeScript up to 2^32 - 1 because its binding took a `u32`.
TypeScript now takes up to 9007199254740991, the largest whole number a
number of JavaScript holds exactly, and its doc comment says so, so the
difference that is left is one the language imposes and not one popnei
chose. The test asserts that a window of that size still keeps the 85
variants of the row whose window is a whole chromosome.

#### What the review of work packages 2 and 3 found

Six reviewers read them: spec, tests, numbers, errors, api and binding.
Nineteen findings held, two were refused with evidence the orchestrator
accepted, and none of them was a wrong number. What the reviewers proved,
each working without the others, is worth as much as what they found:

- One wrote the filter's rule again from the words of
  `docs/specs/filters.md` alone, worked out which variants have two
  dosages from `ld.vcf.gz` itself rather than reading plink2's diagonal,
  and got kept sets identical, variant by variant, to popnei's at all
  four settings. It then compared its own implementation of the r²
  formula with popnei's over 60 random datasets, across ploidies,
  missing genotypes, repeated positions and five thresholds, with no
  disagreement.
- Another wrote the rule again over popnei's own matrix and checked
  eight generated datasets over one to four chromosomes, repeated
  positions, variants of one dosage, windows from 1 to a million and six
  block sizes: no mismatch, and the counts a pass reported equalled the
  variants that came out of it every time.
- Python, TypeScript under WebAssembly and plink2 give a matrix that is
  the same to the bit, and the 500 x 500 one is symmetric to the bit with
  63376 cells holding no r², which is exactly twice the 31654 pairs
  without one plus the 68 variants of the diagonal.

So the rule, the matrix and both bindings give the right numbers. What
the review found is that the code asks too much of a machine, and that
the checks were thinner than they looked.

**The filter did about eighteen times the work its spec asks for.** "How
it runs" of the filter item says the variants of a block are compared
with the window "in one set of the products of `docs/specs/ld.md`, so the
work that the window bounds is done as matrix products and not one pair
at a time". The code called the products once for each variant of the
window. Measured: 2.83 ms for one call of the whole window against a set
of 256 variants of 1000 individuals, against 50.01 ms for 250 calls; end
to end, about 20 s for 100000 variants at a window holding 250 of them,
where a whole pass over that file with no filter takes 0.55 s. Task 3.1
had judged that the window could not be one set of dosages, because a set
is built from one block and a window holds variants of several. Work
package 2 had already solved that same problem for its tiles, in this
plan, a day earlier.

**Three places took memory in a way that ends the process where the spec
promises an error.** The window of the filter, the dosages the filter
builds for a whole block, and the copy the Python binding makes of the
matrix. The TypeScript binding makes the same copy and asks the machine
first, saying why: an allocation that fails in WebAssembly is a trap that
leaves the module unusable. Two subagents met the same problem hours
apart and one of them solved it. Measured: 44 MB of window for 2500 kept
variants of 400 individuals, which is about 32 GB at a million; and 240
MB of dosages for a block of 10000 variants of 1000 individuals, inside a
filter whose spec says its memory is the window's.

**The check on which the filter's rule rests could not fail.** The three
properties of the kept set were worked out over the set the same program
had just built with the same rule, so all three were 0 by construction.
The demonstrations that they can fail were real, and none of them was in
the repository. The third property was also asking a weaker question than
the spec: it took a dropped variant's window from the whole kept set, so
a variant kept after it could justify the drop. Over a damaged set the
loose question names 45 variants dropped for no reason and the strict one
names 64, and over 200 random trials the worst case was 17 against 51.
The properties are now worked out over the file that was written, and
three damaged sets are kept beside them, one for each property, which the
program refuses to write a 0 for.

**The sequential heart of the filter had no test that crossed a block.**
The five unit tests of the rule all ran one block of five variants or
fewer, so they exercised only the comparison of a candidate against the
variants kept inside its own set. Two mutations of the path that carries
the window from block to block survived the whole suite: a pair exactly
at the threshold dropping the candidate, and a pair with no r² dropping
it, the second directly against the spec. The test named for comparing a
candidate with every variant of its window rather than the last kept one
did not fail when reduced to the last kept one.

**The two languages differed in seven ways**, which is what four
subagents building them at the same time without sight of each other
produces. A TypeScript user was told to change `max_num_vars`, an
argument that does not exist in their language; the cap on the variants
of a matrix was documented as 4294967295 and invited 100000 where a
browser cannot exceed 65535, because a `usize` there is 32 bits; a cap of
0 was refused in one language and not the other; the arrays of a result
were read only in one and writable in the other, where the same package
freezes the names of its distances; and the kinds of a pass's filters had
gained the new one in Python and not in TypeScript.

Six statements of the two specs no longer matched the code, and in two of
them the code was right: the memory of the matrix is asked for once the
pass has counted its variants and not before it, which is the only time
the number of variants is known. The others were an error case the spec
did not list, a window said to be trimmed variant by variant where it is
trimmed once for each set settled, a chain of three filters whose counts
no test asserted, and a constant documented to change no result with
nothing varying it.

Two findings were refused by the subagents that had to fix them, with
evidence the orchestrator checked and accepted. The two limits on a
number of base pairs in the TypeScript crate are two different limits and
not one written twice: a position is held exactly by a float64 up to 2 to
the 53, and a window is checked against the largest whole number
JavaScript counts in ones, which is one less. And the refusals of a bad
window and a bad threshold were already tested in TypeScript; only the
source whose variants do not come in order was not, which is the one
refusal this filter alone raises in all of popnei.

#### A decision the orchestrator took and took back

A reviewer found that a machine which cannot give the memory for the
matrix reaches Python as a `ValueError`. That is neither what a user
wrote nor a defect of popnei, and Python has `MemoryError` for exactly
it, so the orchestrator ordered the change. It was wrong to.
`.claude/skills/coding/SKILL.md` records a convention the owner gave on
21 September 2026 with three exceptions in it, a `ValueError` for a wrong
input, a `RuntimeError` for a defect of popnei and an `OSError` for a
file, and adding a fourth changes which exception every user of popnei
catches. The subagent that made the change found what settles it:
`Error::DistancesOfTooManyIndividuals` is the same case in another module
and is a `ValueError`, so the change left popnei answering one question
two ways. It was put back, and the question is in what this plan asks of
the owner.
