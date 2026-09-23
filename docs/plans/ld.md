# Plan: r² between variants, the matrix of it, and the filter by linkage disequilibrium

22 September 2026. State: under way. It builds
from `docs/specs/ld.md`, which has r², the squared correlation between
the dosages of two variants, and the matrix of it for a set of variants;
and from the item "The filter by linkage disequilibrium" of
`docs/specs/filters.md`, which takes out the variants that repeat what a
variant near them on the chromosome already said. Both specs went through
their first readers and their reviews on 22 September 2026. It is carried
out in the worktree `.claude/worktrees/ld` on the branch `plan/ld`, and
the report is `docs/reports/ld.md`.

## In and out

Built: the r² of every pair of two sets of variants inside the core;
`calc_rogers_huff_r2_matrix` in Python and in TypeScript; and
`Variants.filter_by_ld` in both, with its counts and its steps. Each of
the last two ends at the function a user calls and at a test that runs
against pyNei or against plink2's stored numbers.

Not built, with where it goes:

- **The curve of r² against distance per population**,
  `calc_ld_and_dist_per_pop`, the second item of `docs/specs/ld.md`. The
  owner decided on 22 September 2026 to leave it for a later plan:
  nothing calls it today, where the filter is what pop_lab uses to prune
  before a principal component analysis, and it is the most machinery of
  the three, a window over blocks with bins and a major allele frequency
  per population. The spec item stays as it is, written and reviewed, and
  a later plan builds it. Nothing in this plan stands on it: the filter
  takes `LdDosages` and `r2_between` of work package 1, and so would it.
- **Open 1 of `docs/specs/ld.md`**, whether that curve also carries a
  sample of individual pairs, is a point of the item that is not built.
  It needs no answer for this plan and none of its tasks changes with it.
- **Open 2 of `docs/specs/ld.md`**, the major allele of a variant with
  half called genotypes, is task 1.2 and deliverable 1.5. Its
  "meanwhile" is the rule `docs/specs/pca.md` has and that `pca.rs`
  already computes, so under the meanwhile the code does not change and
  the task is the move of that function into `variant`. The other answer,
  that a half called genotype is counted as plink2 counts it, would change
  which allele is the major one in some variants of more than two
  alleles, and so the dosages and the r² that popnei reads from
  `tests/reference/vcf/many.vcf`; it makes a new task in work package 1
  and lets that file be checked against plink2, which deliverable 5
  cannot do under the meanwhile.
- **The variants of a region of a chromosome**, which would make the cap
  of `calc_rogers_huff_r2_matrix` rarely bite: a later item of
  `docs/specs/filters.md`.
- **The genomic relationship matrix**, which reads the same dosages:
  `docs/specs/kinship.md`, which is not written.

## What has to be in place

All three layers exist, so every check of the `coding` skill runs during
this plan and none is reported as not there: `cargo fmt --all --check`,
`cargo clippy --workspace --all-targets -- -D warnings`, `cargo test
--workspace`, `cargo wasm-check`, `uv run ruff format --check && uv run
ruff check`, `uv run maturin develop && uv run pytest`, and `npm run
build && npm test` in `js/popnei`.

On the commit this starts from, which is `main` with the PCA merged:

- `cargo test --workspace` passes, 395 tests in the core crate and 35 in
  the linalg crate.
- `cargo test -p popnei --lib ld:: -- --list` prints `0 tests`,
  `filters::` prints `34 tests` and `variant::` prints `10 tests`.
- `tests/test_ld.py`, `js/popnei/src/ld.ts` and `tests/reference/ld/` do
  not exist.
- `crates/popnei-linalg` gives `product`, which every r² of this plan
  goes through, and `docs/specs/linalg.md` is merged.
- plink2 v2.0.0-a.7.7, bcftools 1.24 and R 4.6.1 are on the machine, and
  the plink2 commands of "How it is verified" of both specs were run
  while the specs were written; `docs/reports/ld-method/README.md` has
  them and the programs that recompute every table.

## Work package 1: the r² of two sets of variants

### What it gives

The core crate works out r² for every pair of two sets of variants, with
the individuals missing at either variant of a pair left out of it. It
stops inside the core: nothing of it reaches Python, and work package 2
is where popnei is first compared with pyNei.

### Its deliverables

1. `tests/reference/ld/` holds `make_reference.py`, `ld.vcf.gz`,
   `example.vcf` and the output of plink2 for both, with the commands in
   a shell script beside them. The script writes the VCF uncompressed, as
   `docs/reports/ld-method/make_ld.py` does now, and the task gzips it to
   `ld.vcf.gz`, which is how `tests/reference/dists/` keeps its own. The
   check: run the script into an empty directory and `diff` the `ld.vcf`
   it wrote against `zcat`ing the `ld.vcf.gz` in git, which prints
   nothing, and rerun the plink2 commands, which give the stored
   matrices. The directory does not exist today.
2. `variant::the_major_allele` is public and `pca.rs` has none of its
   own: `grep -c "fn the_major_allele" crates/popnei/src/pca.rs` is 0 and
   `cargo test --workspace` still passes its 395 and 35 tests, none of
   them changed.
3. The worked example of "How it is verified" of `docs/specs/ld.md`, 5
   variants of 6 individuals, is a cargo test at `r2_between` that
   asserts, for each of its seven pairs, the six whole numbers exactly
   and the r² within 1e-12 relative, with NaN for the pair that holds the
   variant of one dosage.
4. Every pair of `ld.vcf.gz` matches the stored plink2 matrix: a cargo
   test reads the file with the VCF reader, computes r² for all 124750
   pairs, and asserts that the 93096 that have one agree within 1e-12
   relative and that the other 31654 are NaN on both sides. The spec says
   why a difference of 1e-13 there is worth looking at.
5. The dosages of `tests/reference/vcf/many.vcf` are those of pyNei's
   `to_012`. That file has 500 variants of 50 individuals, of which 54
   variants have more than two alleles, and 257 of its 25000 genotypes
   are half called. `make_reference.py` runs pyNei once and stores the
   500 x 50 dosages beside the dataset, and a cargo test at
   `LdDosages::dosages` compares every one of the 25000 with it. It is a
   cargo test and not a pytest one because nothing of this work package
   reaches Python.
6. `cargo test -p popnei --lib ld:: -- --list` prints more than `0
   tests` and names the tests of deliverables 3, 4 and 5.

### What it stands on

`crates/popnei-linalg` and its `product`, both on main. Nothing else in
this plan.

### Its tasks

- [x] 1.1 The reference dataset: `docs/reports/ld-method/make_ld.py`
      moved to `tests/reference/ld/make_reference.py` unchanged, since
      every literal of both specs depends on the order in which it asks
      its generator for its numbers; `ld.vcf.gz`, `example.vcf`, the
      plink2 matrices and the script that runs plink2 beside it. Built
      from "How it is verified" of `docs/specs/ld.md`. Serves deliverable
      1. Needs nothing.
- [x] 1.2 `the_major_allele` public in `variant`: taken out of
      `crates/popnei/src/pca.rs`, where the PCA work left it private and
      exposed a copy for the benches, and put in
      `crates/popnei/src/variant.rs` with the doc comment of "The Rust
      interface" of `docs/specs/ld.md`. The PCA calls it. It changes no
      number. Serves deliverable 2. Needs nothing, and can run beside
      1.1.
- [x] 1.3 `LdDosages` in a new `crates/popnei/src/ld.rs`: the three
      matrices of "How it runs" of `docs/specs/ld.md`, its constructor
      over a block and a set of individuals, `rows`, `has_variance`,
      `dosages` and `maf`, with the errors of "The Rust interface".
      Serves deliverables 5 and 6. Needs 1.2.
- [x] 1.4 `r2_between`: the six products through `linalg::product`, four
      where a set is multiplied by itself, and the r² of every pair from
      the six sums, with NaN where "What it gives" says there is none.
      Its cargo tests are the worked example. Serves deliverables 3 and
      6. Needs 1.3.
- [x] 1.5 The two checks against the stored numbers: the test over every
      pair of `ld.vcf.gz` against the plink2 matrix, and the test of the
      dosages of `many.vcf` against the ones pyNei gave, which this task
      adds to `make_reference.py` and stores. Serves deliverables 4, 5
      and 6. Needs 1.1 and 1.4.
- [x] 1.6 The product with a transposed operand, added to the plan on 23
      September 2026 by the owner's decision after the review of this
      work package: `product_by_transpose` in `crates/popnei-linalg`, on
      both of its backends, with the item and the three checks that
      `docs/specs/linalg.md` now has for it, and `r2_between` calling it
      so that `crates/popnei/src/ld.rs` holds no matrix operation of its
      own, which the `coding` skill does not allow. The r² of every pair
      of `ld.vcf.gz` has to stay equal to plink2's to the bit, on both
      backends, which deliverable 4 is what checks. Serves deliverables
      3, 4 and 6. Needs 1.4.

### What could go wrong

The six sums are whole numbers that an `f64` holds exactly, and the spec
leans on that: if a product through BLAS or through faer ever gave a
value that is not the whole number it should be, every r² would move in
its last bits and deliverable 4 would show it. Both backends have to
pass: run `cargo test -p popnei-linalg` and the tests of this work
package with `--no-default-features` as well, which is how
`docs/specs/linalg.md` asks for its own to be run.

The other thing not measured is the shape: `linalg::product` was built
for the PCA, whose result is individuals x individuals. Here the result
is variants x variants and the inner dimension is the individuals, which
is the transpose of that use. Deliverable 3 has three pairs of different
dimensions for exactly this reason.

## Work package 2: the matrix of every pair

### What it gives

A user calls `calc_rogers_huff_r2_matrix(variants)` in Python, or
`calcRogersHuffR2Matrix` in TypeScript, and gets the r² of every pair of
the variants their `Variants` gives, with the chromosomes and the
positions beside it and the counts of the pass.

### Its deliverables

1. `calc_r2_matrix` in the core with the cap of `max_num_vars`, its
   result and the errors of "The Rust interface" of `docs/specs/ld.md`,
   including the memory of the matrix asked for with `try_reserve_exact`.
   Its cargo tests assert the five literals of the table of "How it is
   verified" with their `n`, that the 68 variants of one dosage have NaN
   in their row, their column and their diagonal cell, and that the
   matrix is the same to the bit for blocks of 7, 64, 256 and 500
   variants and for rayon at 1 thread and at 4.
2. The Python function and `R2Matrix` as "Its Python function" describes
   them, with `pass_stats`. `uv run pytest tests/test_ld.py` passes and
   the file does not exist today.
3. The comparison with pyNei, a pytest test: over `ld.vcf.gz` with
   `filter_by_missing_data(0)` put on the `Variants`, popnei's matrix and
   pyNei's r squared agree within 1e-12 relative; with the missing
   genotypes left in, the test asserts the difference of the table of
   "Missing genotypes", so that the divergence is pinned and a change on
   either side shows.
4. The refusals: a `max_num_vars` below the variants of the pass is a
   `ValueError` whose message carries both numbers and the memory the
   matrix would have needed, and a pass with no variant is a
   `ValueError`.
5. The TypeScript function and its test under node, which reads
   `ld.vcf.gz` from a `Uint8Array` and asserts the five r² and the
   `Error` of a `maxNumVars` of 100. `npm run build && npm test` in
   `js/popnei` passes and `js/popnei/src/ld.ts` does not exist today.

### What it stands on

Work package 1.

### Its tasks

- [x] 2.1 `calc_r2_matrix` and `R2Matrix` in `crates/popnei/src/ld.rs`,
      from "The Rust interface" and "How it runs" of `docs/specs/ld.md`,
      with the tiling that keeps the six intermediate matrices small and
      the lower half computed once. Serves deliverable 1. Needs 1.4.
- [x] 2.2 The Python function: the binding in
      `crates/popnei-python/src/`, which opens the reader of the pass
      from the source and the steps of the `Variants` and keeps it to
      read the counts of the filters, and `calc_rogers_huff_r2_matrix`
      with its result in `python/popnei/ld.py`, from "Its Python
      function" of `docs/specs/ld.md`. Serves deliverables 2, 3 and 4.
      Needs 2.1.
- [x] 2.3 The TypeScript function: the binding in `crates/popnei-js` and
      `js/popnei/src/ld.ts`, from the TypeScript paragraph of the same
      part. Serves deliverable 5. Needs 2.1, and can run beside 2.2.

### What could go wrong

The result is 8 bytes a pair and the default cap is 5000 variants, 200
MB. The test of deliverable 1 runs on 500 variants, 2 MB, so nothing in
the tests reaches the memory the cap is there for; that the refusal
fires is tested with a small `max_num_vars` and not with a large file.

## Work package 3: the filter by linkage disequilibrium

### What it gives

A user writes `variants.filter_by_ld(max_allowed_r2, max_dist)` in Python
or `variants.filterByLd(...)` in TypeScript, and every consumer of that
`Variants` afterwards sees the variants that do not repeat what a variant
kept within `max_dist` on their chromosome already said, with the counts
of what the filter was given and kept in the `pass_stats` of the result.

### Its deliverables

1. `LdFilter` and `LdFilteredReader` in `crates/popnei/src/filters.rs`
   with the rule and the window of "What it gives" and "How it runs" of
   the item in `docs/specs/filters.md`, and the error of a position that
   goes backwards within a chromosome or of a chromosome that comes back.
2. The four rows of the table of "How it is verified": a cargo test at
   `next_block` of an `LdFilteredReader` over a `VcfReader` on
   `ld.vcf.gz` asserts 84, 133, 85 and 85 variants kept of 500 and the
   five positions the table gives, with blocks of 7, of 64 and of the
   default size keeping the same variants.
3. The three properties of the kept set, checked in
   `tests/reference/ld/make_reference.py` against the stored plink2
   matrix and their result stored beside it: every kept variant has two
   dosages at least, no two kept variants inside the window are above the
   threshold, and every dropped variant with two dosages is above it
   against a variant kept before it inside its window.
4. `MaxLdR2`, `max_dist()` and the kind `"ld"` on
   `VarFilteringCriterion`, `chain_of` building an `LdFilteredReader` for
   it, and `refuse_a_second_filter_of_a_kind` covering it. `cargo test -p
   popnei --lib filters:: -- --list` prints more than the `34 tests` of
   today and names the new ones.
5. `Variants.filter_by_ld` in Python with the pytest tests that "How it
   is verified" of the item names: the four counts and the five
   positions, the `Step` of kind `"ld"` with both arguments, the
   `filtering` of the `pass_stats` after a maf filter, the `ValueError`
   of each refused argument, the `ValueError` of a second filter of the
   kind, and the `ValueError` of a source whose positions go backwards.
6. The TypeScript step and its test under node: 133 variants kept of 500
   with their five positions and the `Error` of a `maxAllowedR2` of 1.5.

### What it stands on

Work package 1. It does not stand on work package 2 and the two can run
side by side: this one touches `filters.rs`, `python/popnei/filters.py`
and `js/popnei/src/filters.ts`, and that one `ld.rs`, `python/popnei/ld.py`
and `js/popnei/src/ld.ts`.

### Its tasks

- [x] 3.1 `LdFilter`: the window, the rule, the counts and the error of
      an unsorted position, in `crates/popnei/src/filters.rs`, built from
      "What it gives", "Which variant of a linked pair is kept, and the
      cases" and "How it runs" of the item in `docs/specs/filters.md`. It
      is the task whose failure would be silent, a set of variants that
      is wrong and not a crash, and deliverable 2 is what guards it.
      Serves deliverables 1 and 4. Needs 1.4.
- [x] 3.2 `LdFilteredReader`, the criterion and the chain: the reader
      over a reader with the rules of `docs/specs/block.md`, `MaxLdR2`
      and `max_dist()` on `VarFilteringCriterion`, `chain_of` and
      `refuse_a_second_filter_of_a_kind`, from "The Rust interface" of
      the same spec. Serves deliverables 1 and 4. Needs 3.1.
- [x] 3.3 The checks against the stored numbers: the three properties
      added to `tests/reference/ld/make_reference.py` with their output
      stored, and the cargo test of the four rows of the table. Serves
      deliverables 2 and 3. Needs 1.1 and 3.2.
- [ ] 3.4 The Python step: the binding and
      `Variants.filter_by_ld` in `python/popnei/filters.py`, with the
      `Step`, the counts and every refusal of "In Python and in
      TypeScript" of the item, and the pytest tests it names. Serves
      deliverable 5. Needs 3.2.
- [ ] 3.5 The TypeScript step: the binding in `crates/popnei-js` and
      `js/popnei/src/filters.ts`, with its test under node. Serves
      deliverable 6. Needs 3.2, and can run beside 3.4.

### What could go wrong

The rule is sequential where everything else in popnei is not: whether a
variant is kept decides what the next one is compared with, so the
variants of one block cannot all be settled at once. The spec leaves how
much of a block is taken at a time to the implementer and asks only that
the result not change with the blocks, which deliverable 2 checks with
three block sizes.

The filter is also the one reader of popnei that refuses a source the
rest of it takes, a chromosome whose positions do not rise.
`docs/specs/io_vars.md` has a cargo test that writes a block of four
variants that are not sorted, so such a source exists in the tests
already and the error has to name the variant and both positions.

## Work package 4: the measurements

### What it gives

The "Speed" sections of both specs, which today say that this plan fills
them, carry the numbers popnei reaches.

### Its deliverables

1. A bench, `crates/popnei/benches/r2_matrix.rs`, that times
   `calc_r2_matrix` on 5000 variants of 1000 individuals with 3 in 100
   genotypes missing, read from a vars file, with the reading timed apart
   and taken out. The number to reach is the 0.50 s of "Speed" of
   `docs/specs/ld.md`, which is the 0.455 s numpy takes on Accelerate for
   the same six products with a tenth over it.
2. A bench that times a whole pass of the filter over the 400 MB VCF of
   `docs/rust_core.md`, 100000 variants of 1000 individuals, against the
   same pass with no filter, as `crates/popnei/benches/filter_vars.rs`
   does for the three threshold filters. There is no number to reach;
   "Speed" of `docs/specs/filters.md` gets what was measured.
3. The same two timed in the browser, the wasm package built and run
   under node as `docs/specs/pca.md` had it done, with the vector
   instructions that `.cargo/config.toml` now passes.
4. Both "Speed" sections hold the numbers with the machine, the dataset
   and the date, and no longer say that the plan will measure them.

### What it stands on

Work packages 2 and 3.

### Its tasks

- [ ] 4.1 The two benches and the native numbers. Serves deliverables 1
      and 2. Needs 2.1 and 3.2.
- [ ] 4.2 The numbers in the browser. Serves deliverable 3. Needs 2.3
      and 3.5.
- [ ] 4.3 The two "Speed" sections. Serves deliverable 4. Needs 4.1 and
      4.2.

### What could go wrong

The target of deliverable 1 comes from numpy on Accelerate, which runs
its products on the matrix units of the Apple chip, and the core calls
the same Accelerate through `linalg`. What popnei adds is the building of
the three matrices, 24 bytes for each variant and individual, and the
element wise arithmetic that turns the six products into r², which reads
six matrices of 200 MB and writes one, 1.4 GB of traffic. numpy pays that
traffic too, and it is part of the 0.455 s, so a core that fuses those
six readings into one pass could come in under the target rather than
over it. If the target is missed, the report says by how much and the
owner decides whether a performance review follows, as
`docs/reports/pca.md` did.

## How the whole plan is checked

The sum of the work packages, and two things besides.

The seven checks of the `coding` skill pass on the last commit, the four
cargo ones and the three of Python and TypeScript, and `cargo test -p
popnei-linalg --no-default-features` passes too, so that both backends of
the linear algebra are known to give the numbers of the literals.

A run that no work package makes on its own: `calc_rogers_huff_r2_matrix`
and `filter_by_ld` are called from Python on
`tests/reference/vcf/many.vcf`, which has 54 variants of more than two
alleles and 257 half called genotypes, with the filter before the matrix
on one `Variants`, and the result has the counts of both steps. It
checks nothing about the values, which "How it is verified" of
`docs/specs/ld.md` pins elsewhere; it checks that the two things this
plan builds work on one dataset through one pass.
