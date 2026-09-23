# Plan: the distances between populations

23 September 2026. Approved and under way. It builds the item "Distances
between populations" of `docs/specs/dists.md` and the five items after it,
Hudson's F_ST, f_2, the chord distance with Nei's D_A, Jost's D, and Nei's
G_ST with the standardized G''_ST: one pass over the variants that gives
seven numbers for every pair of populations, each with a jackknife
standard error, through the core crate, both binding crates and both
packages. The spec has no open points. It is carried out as the
`following-plans` skill says, on the branch `plan/dists-pops`.

## In and out

In: the pass and its sums; the resampling groups and the delete-m
jackknife; the seven measures; `calc_pop_dists` in Python and
`calcPopDists` in TypeScript with their results; the tests against plink2,
adegenet, mmod, ADMIXTOOLS 2 and pyNei on the two panels.

Out, with where each goes:

- The speed. The spec's "Speed" has no number to reach and the owner
  decided on 23 September 2026 that the performance review measures it
  after this work is merged into `main`, not a work package here.
- f_3 and f_4, which are sums and differences of the f_2 of pairs: a spec
  of their own. What this plan leaves them is `f2_groups`, the f_2 of each
  pair within each resampling group, so they need no new pass.
- R_ST and the squared difference in repeat number: they need the allele
  sizes, which popnei's genotypes do not carry, and the spec's "Not in
  this spec" says why popnei does not add them.
- The principal coordinate analysis of these matrices: an item of
  `docs/specs/pca.md` that is not written.
- The read ahead thread: `docs/specs/block.md`, and nothing here changes
  with it, since every calculation takes any reader.

## What has to be in place

All three layers exist, so every check of the `coding` skill runs during
this plan and none is reported as not there.

In the core crate on `main`: `Pops` and `ObsHet` of `stats`, which give
the populations and the called and heterozygous genotypes of one
population at one variant; `count_alleles_of` and `count_gts_of` of
`variant`; `BlockReader` and the blocks of `block`; `chain_of` of
`filters`; and `Error::PassGaveNoVariant`. In Python: `Variants`,
`PassStats` and `Distances`, which gains one field here. In TypeScript:
the package with its own `Distances`.

The spec, on the branch `spec/dists-pops`, which this branch is made from.
If that branch is merged into `main` before this one starts, this one is
rebased on it and nothing else changes.

The reference data, written and run on 23 September 2026 and in
`tests/reference/pop_dists/`: the multiallelic panel `micro.vcf.gz` with
`micro_pops.txt`, and, for it and for the biallelic
`tests/reference/dists/panel.vcf.gz`, the output of every reference
program. `make_reference.py` there rebuilds all of it and refuses a
version of any program other than the ones the numbers were taken with.

The reference programs, all present on the owner's machine on 23 September
2026: plink2 v2.0.0-a.7.7; R 4.6.1 with adegenet 2.1.11, mmod 1.3.3 and
admixtools 2.0.10; and pyNei at the commit `[tool.uv.sources]` names.

Every check below was run on the commit this branch starts from. The cargo
ones print `0 tests` there, which is why each names how many have to run:
`cargo test -p popnei --lib pop_dists:: -- --list` prints `0 tests, 0
benchmarks` and exits 0 today, so a check that only names a selector would
pass on an empty crate. `tests/test_pop_dists.py` and
`js/popnei/test/pop_dists.test.ts` do not exist, and pytest exits 5 when
its `-k` matches nothing.

## Work package 1: the pass, the resampling groups and the first two measures

### What it gives

A user calls `calc_pop_dists(variants, pops, jackknife_group=...)` in
Python, or `calcPopDists` in TypeScript, and gets Hudson's F_ST and f_2 for
every pair of populations, each with its jackknife standard error, with
`num_vars` per pair, `f2_groups`, `group_ids` and the counts of the pass.
The other five measures raise until the work packages below add them.

### Its deliverables

1. The three counts of one variant in one population and what a pair makes
   of them. Check: cargo tests in `crates/popnei/src/pop_dists.rs` at the
   function that adds one variant to the sums of a pair, asserting the
   counts table and the per variant table of "How it is verified" of the
   spec's shared item, H_b, H_w and f_2 of all four variants of the worked
   example, variant 3 among them, where both populations are fixed for
   allele 0 and H_b is 0, and variant 4, where f_2 is -0.1; and the
   corrected H_S and H_T of the same table, which no measure reads until
   work package 2 and which this function builds with the rest.
2. The resampling groups. Check: cargo tests at the walk that cuts them,
   asserting the 12 groups of 100 variants that a length of 100 000 base
   pairs cuts the biallelic panel into, over its two chromosomes, which is
   what "The standard errors" of the spec says ADMIXTOOLS 2 cuts; that a
   new group starts at each chromosome; `"variant"`, which gives one group
   for each variant; `None`, which gives none; and that a length of 0 or
   less is an error.
3. F_ST, f_2 and the standard error out of the sums. Check: cargo tests at
   `PopDistSums::measure` and `PopDistSums::standard_error` with the
   worked example's F_ST of 0.276712 and f_2 of 0.140278; the three
   plink2 literals of the spec's F_ST item on each panel, within 1e-6
   absolute; the per variant F_ST of `var0000`, 0.332322, over a reader of
   that one variant; and ADMIXTOOLS 2's three f_2 and three standard
   errors of the spec's f_2 item, within 1e-12 relative, over the same 12
   groups.
4. The calculation over a reader, with its errors. Check: `cargo test -p
   popnei --lib pop_dists:: -- --list` prints 20 tests or more, `0 tests`
   today. Beside those of the three deliverables above: the panel read in
   blocks of 100 and of 10000 and in rayon pools of 1 and 4 threads giving
   the same numbers within 1e-12 relative, which is the invariance the
   spec asks of every measure; a pass with no variant; fewer than two
   populations; fewer than 20 groups; and an error of the reader given on
   as it is.
5. Python. Check: `uv run pytest tests/test_pop_dists.py` passes, where
   the file does not exist today, with the plink2 and ADMIXTOOLS literals
   through `calc_pop_dists` on both panels; a call with no
   `jackknife_group`, which is a `TypeError` of Python itself; `None` for
   it, which gives no standard errors and asks the reader for no
   positions; the `ValueError` of a length of 0, of fewer than 20 groups
   with the number in the message, of fewer than two populations, of a
   name that is not an individual, and of a pass that kept no variant,
   whose message names the counts of the filters; `num_vars` per pair and
   `pass_stats`; and `Distances.standard_errors` being `None` on the
   Kosman distances, which do not change.
6. TypeScript. Check: `npm test` in `js/popnei` gives 0 failed and
   `test/pop_dists.test.ts`, which does not exist today, runs the same
   plink2 and ADMIXTOOLS literals through `calcPopDists` on both panels
   with `standardErrors`, `numVars` and `passStats`, and asserts the
   `Error` of a missing `jackknifeGroup` and of a source with no variant.

### What it stands on

Nothing in this plan. Outside it, everything of "What has to be in place".

### Its tasks

- [x] 1.1 The counts of one variant in one population and the sums of one
  pair, in a new `crates/popnei/src/pop_dists.rs`: the allele frequencies
  over the called alleles, H_b, H_w, the sum of the square roots of the
  products, the corrected H_S and H_T, and the test that says whether the
  variant counts for the pair. From "What it gives" and "Variants that do
  not count, populations with little data, and negative values" of the
  spec's shared item, and from "What it gives" of Jost's D for the two
  corrected values. Serves deliverable 1. Its failure would be silent, a
  wrong number and not a crash, so the worked example guards it and it has
  a commit of its own. The two corrected values are computed here although
  no measure reads them until work package 2, so that the sums of a pair
  are written once. Deliverable 1 asserts them, so nothing of this task
  goes unchecked until work package 2.
- [x] 1.2 The resampling groups: `JackknifeGroups`, `GroupId` and the walk
  that cuts them, in the same file. From "The standard errors". Serves
  deliverable 2. Can run side by side with 1.1.
- [x] 1.3 F_ST, f_2 and the delete-m jackknife out of the six sums, in the
  same file. From "What it gives" of Hudson's F_ST and of f_2, and from
  "The standard errors" for the pseudo-values, the jackknife estimate and
  the variance. Serves deliverable 3. Needs 1.1 and 1.2. Its failure would
  be silent too, so it has a commit of its own and the worked example and
  the ADMIXTOOLS numbers guard it.
- [ ] 1.4 `calc_pop_dist_sums`, `PopDistSums` and `PopDistOptions`: the
  loop over the blocks of a reader, the rows on rayon natively and one
  after another in wasm as `calc_per_var_distribs` of `stats` does it, the
  accumulator of six numbers for each pair and each group, the errors, and
  the cargo tests of deliverable 4. From "How it runs" and "The Rust
  interface". Serves deliverables 3 and 4. Needs 1.1, 1.2 and 1.3.
- [ ] 1.5 The Python side: the function of `crates/popnei-python` that
  builds the chain with `chain_of`, runs the pass and gives back the
  measures, the standard errors, the counts and `f2_groups`, as the Kosman
  distances do; `PopDists` and `calc_pop_dists` in a new
  `python/popnei/pop_dists.py`, the `standard_errors` field and
  `square_standard_errors` on `Distances` in `python/popnei/dists.py`, both
  exported from `popnei`; and `tests/test_pop_dists.py`. From "Its Python
  function" and the cases section. Serves deliverable 5. Needs 1.4.
- [ ] 1.6 The TypeScript side: the same in `crates/popnei-js`,
  `calcPopDists` and `PopDists` in `js/popnei/src/pop_dists.ts` with the
  `standardErrors` of `dists.ts`, and `test/pop_dists.test.ts`. From the
  TypeScript paragraph of "Its Python function". Serves deliverable 6.
  Needs 1.4; can run side by side with 1.5.

### What could go wrong

The accumulator is six numbers for each pair and each group, and how many
groups there are is not known until the variants have been read, so it
grows as groups appear rather than being allocated once as the Kosman sums
are. The spec's figures, 66 KB for 3 populations and 500 groups and 27 MB
for 50 populations, are what it has to stay near.

The rows of a block go on rayon and every row adds to the sums of every
pair, so the threads reduce into per pair, per group sums rather than
writing into one. The result has to be the same to 1e-12 relative whatever
the number of threads and the size of the blocks, which deliverable 4
checks; getting that wrong is the silent failure this work package is
likeliest to have.

`f2_groups` is groups x pairs and crosses to Python without a copy and to
TypeScript with one, as section 11 of `docs/architecture.md` requires.

## Work package 2: Jost's D, Nei's G_ST and the standardized G''_ST

### What it gives

The three measures that come from the corrected H_S and H_T, so that a
user comparing microsatellite populations gets the allelic differentiation
beside the fixation measures of work package 1, and so that popnei's D can
be compared with pyNei's, which is the strongest check popnei has.

### Its deliverables

1. The three measures out of the sums. Check: cargo tests at
   `PopDistSums::measure` with the worked example's D of 0.248677, G_ST of
   0.191837 and G''_ST of 0.490541, and with the mmod literals of the two
   items on both panels, within 5e-4.
2. The variants that do not count. Check: cargo tests that a variant where
   both populations have exactly one called genotype counts for no measure
   at a `min_num_individuals` of 1, the case of "Variants that do not
   count" that pyNei reaches through a NaN; and that `measure` gives
   `None` for Dest when the mean corrected H_S is exactly 1.
3. Python, and the comparison with pyNei. Check: `uv run pytest
   tests/test_pop_dists.py` passes with, beside work package 1's tests,
   popnei's `dest` against pyNei's `calc_jost_dest_pop_dists` on the
   biallelic panel within 1e-12 relative at a `min_num_individuals` of 20;
   the same at 47, where the pairs part, asserting pyNei's 0.0595097904,
   0.0612873957 and 0.0656705213 and the counts 688, 688 and 1200; and the
   mmod literals of both panels within 5e-4.
4. TypeScript. Check: `npm test` gives 0 failed with the three measures of
   both panels in `test/pop_dists.test.ts`.

### What it stands on

Work package 1. Its task 1.1 already computes the corrected H_S and H_T
and its task 1.4 already sums them, so this work package adds no count and
no sum.

### Its tasks

- [ ] 2.1 Jost's D, G_ST and G''_ST out of the means of the corrected sums,
  in `crates/popnei/src/pop_dists.rs`, with the cargo tests of deliverables
  1 and 2. From "What it gives" of Jost's D and of the G_ST item, and from
  "What pyNei does that is odd, and what popnei does instead". Serves
  deliverables 1 and 2. Needs work package 1. Its failure would be silent,
  so it has a commit of its own.
- [ ] 2.2 The three measures through both binding crates and both
  packages, with the doc comment of `dest` that "What it gives" of Jost's D
  asks for, which names the estimator and gives the two differences from
  mmod; and the tests of deliverables 3 and 4. From "Its Python function"
  for where the three sit in the result. Serves deliverables 3 and 4.
  Needs 2.1.

### What could go wrong

The correction divides by the harmonic mean of the called genotypes minus
one, and the spec's rule keeps that from being zero by dropping the variant
where both populations have one called genotype. A test at a
`min_num_individuals` of 1 is the only place that rule is exercised, since
the two panels never reach it.

The comparison with pyNei is exact to 1e-12 relative and the one with mmod
is not, because mmod computes another estimator of the same quantity. A
subagent that sees the mmod gap and tightens the tolerance, or that changes
the arithmetic to close it, has broken the agreement with pyNei; the spec's
"How it is verified" of Jost's D says which is which.

## Work package 3: the chord distance and Nei's D_A

### What it gives

The two measures a user builds a tree or a principal coordinate analysis
of populations from, and the only ones of the seven that are Euclidean.

### Its deliverables

1. The two measures. Check: cargo tests at `PopDistSums::measure` with the
   worked example's D_A of 0.161566 and chord of 0.401953, and the
   adegenet literals of the spec's chord item on both panels within 1e-12
   relative.
2. Python and TypeScript. Check: `uv run pytest tests/test_pop_dists.py`
   and `npm test` pass with `chord` and `da` of both panels, and with
   `da` being the square of `chord` for every pair.

### What it stands on

Work packages 1 and 2. Task 1.1 already computes the sum of the square
roots of the products and task 1.4 already sums it, so this work package
adds no count and no sum either.

### Its tasks

- [ ] 3.1 The chord distance and Nei's D_A out of that sum, in
  `crates/popnei/src/pop_dists.rs`, with the cargo tests of deliverable 1.
  From "What it gives" of the chord item. Serves deliverable 1. Needs work
  package 1.
- [ ] 3.2 The two measures through both binding crates and both packages,
  with the tests of deliverable 2. Serves deliverable 2. Needs 3.1.

### What could go wrong

The form popnei computes is adegenet's, the chord divided by the square
root of two, which the spec's "What it gives" says. A subagent that writes
the chord of the unit sphere instead will be a factor of 1.414 from every
adegenet literal.

## How the whole plan is checked

The sum of the work packages, and one thing more that no single work
package covers: `uv run pytest tests/test_pop_dists.py` with every measure
asked for at once on both panels, which is what a user calling
`calc_pop_dists` with the default `measures=None` gets, asserting that the
seven numbers of a pair are over the same variants, one `num_vars` for the
pair, and that asking for one measure gives the same number as asking for
all seven. Each work package is reviewed as the `code-review` skill says
before it is reported as done.
