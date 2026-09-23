# Plan: the kinship

23 September 2026. State: draft, not yet approved by the owner. It builds
`docs/specs/kinship.md`, which went through its first reader and its review
and whose four open points the owner decided on 23 September 2026, so it has
none left. It is the first of three plans that together build the
association study; `gwas-linear` and `gwas-logistic` follow it and both take
its output.

It is carried out in the worktree `.claude/worktrees/kinship` on the branch
`plan/kinship`, which starts from `main` with `spec/gwas` merged into it:
`main` already holds the linear algebra this plan needs, merged there on 23
September 2026, and `spec/gwas` holds the two specs and the reference data.
The report is `docs/reports/kinship.md`.

The work packages run in order, and so do the tasks inside each one,
except where a task says it can run beside another. What a task says it
needs is what has to be committed before it starts.

## In and out

Built, through the core crate, both binding crates and both packages:

- `calc_kinship`, the genomic relationship matrix of a set of variants, with
  its per pair denominators, checked against plink2's `--make-rel square` on
  two panels and against pyNei.
- `Kinship.principal_components`, and `filter_individuals`.
- One change to code that exists: the row pass that turns a variant into its
  standardized dosages, which `pca.rs` holds privately today, becomes one
  helper that both modules call with their own divisor. It changes no
  number.

Not built, with where it goes:

- **The association study**, which is what the mixed models of
  `docs/specs/gwas.md` take a kinship for: the plans `gwas-linear` and
  `gwas-logistic`.
- **A kinship from a pedigree, and a sparse one**, which popnei does not
  have: "Not in this spec" of the spec.
- **The pruning by linkage disequilibrium** a user normally does before a
  kinship, which is already built and is a step of the `Variants`:
  `docs/specs/filters.md`.
- **What it costs.** "Speed" of the spec asks for the kinship of 100000
  variants x 1000 individuals against plink2's 0.23 s, and for the same run
  with 3 in 100 genotypes missing, which adds the second product. The owner
  decided on 23 September 2026 that the measurements of a plan are made in a
  session of their own once it is merged, with the `performance-review`
  skill, and not as its last work package.

The spec has no open point.

## What has to be in place

Nothing outside the repository. The linear algebra this plan calls, the
product and the eigendecomposition, is on `main`. Measured on 23 September
2026 in a throwaway worktree holding the starting commit above, `main` with
`spec/gwas` merged in:

- `cargo fmt --all --check` and `cargo clippy --workspace --all-targets --
  -D warnings`: clean.
- `cargo test --workspace`: 604 passing in the core crate, 2 ignored, and
  149 in the linear algebra crate. The same with `--no-default-features` on
  the linear algebra crate, which is the faer backend: 136.
- `cargo test -p popnei --lib kinship -- --list`: `0 tests`.
- `cargo test -p popnei --lib pca -- --list`: `47 tests`.
- `uv run maturin develop && uv run pytest`: 347 passed, 0 failed.
- `npm run build && npm test` in `js/popnei`: 242 pass, 0 fail.

Every layer of popnei exists, so every check of the `coding` skill runs
during this plan and none is reported as not there. The two notes of
`docs/plans/linalg-gwas.md` still hold: `uv run maturin develop --release`
makes `test_a_ctrl_c_while_write_vars_runs_is_raised_and_leaves_no_file`
fail, and the debug build is the check; and `npm run build` in `js/popnei`
needs `wasm-bindgen` on the PATH.

The reference data is already in the repository, committed with the spec:
`tests/reference/kinship/`, holding `panel_called.vcf.gz`, the two
`.plink2.rel.gz` matrices with their `.id` files, and
`make_reference.py`, which runs plink2 again and checks the eleven literals
of the spec against what it produced. It was run on 23 September 2026 and
all eleven matched. The second panel is `tests/reference/dists/panel.vcf.gz`,
which was already there.

## Work package 1: one row pass for the PCA and the kinship

### What it gives

Nothing to a user. It is here first because the kinship's row pass is the
PCA's with one number changed, the divisor, and building it twice would let
the two drift: a change to the dosage rule of `docs/specs/pca.md`, which
both specs take, would then have to be made in two places and could be made
in one. The spec leaves the choice between one helper and two copies to this
plan, and this plan chooses one helper.

### Its deliverables

1. `the_standardized_row`, `RowScratch` and the two helpers they call live in
   `crates/popnei/src/variant.rs`, which section 9 of
   `docs/architecture.md` gives the row helpers over a variant, and take the
   divisor from their caller. The check: `grep -c "fn the_standardized_row"
   crates/popnei/src/*.rs` gives 1 in `variant.rs` and 0 in `pca.rs`, where
   today it is 0 and 1.
2. The PCA computes what it computed. The check: `cargo test -p popnei --lib
   pca -- --list` still prints `47 tests` and `cargo test --workspace`
   passes with 604 in the core crate, with no assertion of `pca.rs` changed;
   and `uv run pytest tests/test_pca.py` passes.
3. The message a user sees for a variant with more than two alleles is
   unchanged for the PCA and is the same wording for a caller that is not
   the PCA, naming the argument that turns it off. The check: the cargo test
   that asserts that message passes untouched, and `cargo test -p popnei
   --lib more_than_two_alleles -- --list` names a test for each of the two
   callers.

### What it stands on

Nothing of this plan.

### Its tasks

- [ ] 1.1 Move the row pass into `crates/popnei/src/variant.rs` with the
      divisor as an argument, leave `pca.rs` calling it with the standard
      deviation of the dosages, and give the error of a variant with more
      than two alleles a form that does not name the PCA while keeping the
      message the PCA's tests assert. Built from "How it runs" of
      `docs/specs/kinship.md` and "The PCA of the variants" of
      `docs/specs/pca.md`. Serves deliverables 1, 2 and 3. Needs nothing.

### What could go wrong

The 47 PCA tests are what say that this changed no number, and they say it
only if they are not touched. A task that edits an assertion of `pca.rs` has
removed the check it was being checked by. The divisor is the only
difference between the two callers; anything else that has to be
parameterized to make the move work is a sign that the two passes are less
alike than the spec says, and that goes back to the owner rather than into a
second parameter.

## Work package 2: the matrix

### What it gives

`calc_kinship(variants, individuals=None, transform_to_biallelic=False)` in
Python and `calcKinship` in TypeScript, giving the matrix of every pair of
individuals, `num_vars` and the pass stats, with `filter_individuals` on the
result.

### Its deliverables

1. The matrix of both reference panels is plink2's. The check: a pytest test
   reads `tests/reference/kinship/panel_called.vcf.gz` and
   `tests/reference/dists/panel.vcf.gz` with `open_vcf`, calls
   `calc_kinship`, and every one of the 40000 entries of each is within 1e-5
   absolute of the stored `.plink2.rel.gz`, with `num_vars` 1200 for both.
2. The worked example of "How it is verified" is a cargo test, asserted
   within 1e-12 absolute, and the eleven plink2 literals of that section are
   cargo tests on the two VCFs read with the VCF reader. The check: `cargo
   test -p popnei --lib kinship -- --list` names them and prints at least 12
   tests, where today it prints `0 tests`.
3. Every case of "Missing genotypes, variants with no variance, and what
   pyNei asserts" has a test: the half called genotype, the variant where
   every individual is heterozygous, the pair with no variant called in both,
   which raises, and the dataset where no variant varies, which raises. The
   check: `cargo test -p popnei --lib kinship::cases -- --list` names one per
   case.
4. The individuals argument takes the frequencies of those individuals. The
   check: a pytest test gets 1195 variants and a largest difference of 0.129
   from the same 40 rows and columns of the whole panel's kinship, the two
   numbers of "Its Python function, and its TypeScript one".
5. popnei and pyNei agree. The check: a pytest test runs both on both panels
   and every entry agrees within 1e-12 relative and `num_vars` exactly.
6. `calcKinship` under node gives the same numbers. The check: `npm test` in
   `js/popnei` runs a test that asserts the four `s000` literals on
   `panel_called.vcf.gz` and the worked example.

### What it stands on

Work package 1. Outside the plan: the VCF reader, `reblock`, the product and
the eigendecomposition of the linear algebra crate, and the `Variants` of
both packages, all of which exist.

### Its tasks

- [ ] 2.1 `Kinship` and `calc_kinship` in a new `crates/popnei/src/kinship.rs`:
      the one pass, the two accumulators, the per pair denominators and the
      refusals, with the cargo tests of the worked example, the plink2
      literals and the cases. Built from "What it gives", "Missing
      genotypes, variants with no variance, and what pyNei asserts", "How it
      runs" and "The Rust interface" of `docs/specs/kinship.md`. Serves
      deliverables 2 and 3. Needs 1.1.
- [ ] 2.2 The Python function: the binding in `crates/popnei-python`, the
      `Kinship` frozen dataclass with its checks on a matrix a user built,
      `filter_individuals`, and the pytest tests against plink2, against
      pyNei and for the individuals argument. Built from "Its Python
      function, and its TypeScript one" and "How it is verified". Serves
      deliverables 1, 4 and 5. Needs 2.1.
- [ ] 2.3 The TypeScript function: the binding in `crates/popnei-js`, the
      result object in `js/popnei`, and the test under node. Built from the
      same section. Serves deliverable 6. Needs 2.1, and it can run beside
      2.2.

### What could go wrong

Two things in the spec are easy to build wrong and neither crashes. The per
pair denominator counts only the variants that were used, not every variant
the reader gave; with every variant the worked example gives 4 and 3 where
the right numbers are 2 and 1. And the matrix of denominators is allocated
with the count of the blocks already passed in every entry and not at zeros,
which pyNei gets for nothing from numpy broadcasting a scalar over a matrix;
an implementation that allocates zeros loses every variant before the first
missing genotype. Both are caught by the panel with genotypes missing and by
neither of the fully called cases, so a failure on `panel` with a pass on
`panel_called` is one of these two.

## Work package 3: the principal components

### What it gives

`Kinship.principal_components(num_pcs)` in Python and
`principalComponents(numPcs)` in TypeScript, the directions along which the
panel varies most, taken from the kinship and ready to be passed to
`calc_gwas` as covariates.

### Its deliverables

1. The components are pyNei's up to sign, and their sign is fixed. The
   check: a pytest test runs both libraries on both panels and the absolute
   value of every projection of the first 10 components agrees within 1e-9
   relative, and a cargo test asserts the sign rule of `docs/specs/pca.md`
   holds in every component.
2. The eigenvalues are right. The check: a cargo test asserts that the sum
   of the squares of each of the first three components' projections on
   `panel_called` is 17.26914116, 12.44731524 and 3.35871258 within 1e-9
   relative, which is the check "How it is verified" describes and the
   reason it is made that way.
3. A component whose eigenvalue is below the tolerance is not given. The
   check: a cargo test on a kinship with fewer components than asked gets
   `num_comps` below `num_pcs`, where pyNei raises out of pandas.
4. The first component separates the subpopulations. The check: a pytest
   test using the `pop` column of `tests/reference/gwas/phenotypes.csv` finds
   the standard deviation of the mean of `PC0` over the three
   subpopulations above the standard deviation of `PC0`.
5. `principalComponents` under node gives the same numbers. The check: `npm
   test` in `js/popnei` runs it against the eigenvalue literals.

### What it stands on

Work package 2, whole, since the work packages run in order. Of its tasks
3.1 needs only what 2.1 built.

### Its tasks

- [ ] 3.1 `principal_components` and `KinshipPcs` in the core, the
      eigendecomposition, the tolerance and the sign rule, with the cargo
      tests of the eigenvalues, the sign and the cut. Built from "The
      principal components of the kinship" of `docs/specs/kinship.md` and
      the sign rule of `docs/specs/pca.md`. Serves deliverables 1, 2 and 3.
      Needs 2.1.
- [ ] 3.2 The Python method and the TypeScript one, on both bindings and in
      both packages, with the pytest tests against pyNei and for the
      subpopulations and the node test. Built from "Its Python function, and
      its TypeScript one" of that item. Serves deliverables 1, 4 and 5.
      Needs 3.1, 2.2 and 2.3.

### What could go wrong

pyNei takes the square root of the absolute value of an eigenvalue below 0
and popnei does not, so the two agree only on the components above the
tolerance, which on the panel with genotypes missing is not all of the first
10 if the tolerance is got wrong. The comparison with pyNei is over the
first 10 components and the panel's third eigenvalue is 3.36 against a
tolerance of 7.67e-13, so a disagreement there means the tolerance, not the
eigendecomposition.

## How the whole plan is checked

The sum of the work packages, and one thing more: the wasm builds. `cargo
wasm-check` is clean, which covers both wasm targets with the warnings
denied, and `npm run build && npm test` in `js/popnei` passes, so the
kinship reaches a browser and not only a native build.
