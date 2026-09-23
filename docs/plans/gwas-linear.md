# Plan: the association study of a continuous trait

23 September 2026. State: under way since 23 September 2026. It builds
the parts of `docs/specs/gwas.md` that a continuous trait needs: the two
distributions, everything `calc_gwas` shares whatever the model, the linear
model and the linear mixed model. It is the second of three plans; `kinship`
comes before it and `gwas-logistic` after.

It is carried out in the worktree `.claude/worktrees/gwas-linear` on the
branch `plan/gwas-linear`. It was to start from `main` with the plan
`kinship` merged into it; the owner decided on 23 September 2026 that it
starts from `b387def` of `plan/kinship` instead, which is that plan's work
packages 1 and 2 and is not yet on `main`, and that `plan/kinship` is merged
in again before work package 4, the first that uses the kinship matrix. The
report is `docs/reports/gwas-linear.md` and gives the reason.

The work packages run in order, and so do the tasks inside each one,
except where a task says it can run beside another. What a task says it
needs is what has to be committed before it starts.

## In and out

Built, through the core crate, both binding crates and both packages:

- The two functions that turn a statistic into a p-value, which every test
  of both plans ends in.
- `calc_gwas(variants, phenotype, trait="continuous", ...)`, with the
  covariates, the design and its refusals, the dosages of a block, and the
  result with one row per variant.
- The linear model, checked against plink2's `--glm`.
- The linear mixed model with both its tests, checked against rrBLUP's
  `GWAS` with `P3D` and against GMMAT's `glmm.score`.

Not built, with where it goes:

- **A binomial trait.** `calc_gwas(trait="binomial")` raises until
  `gwas-logistic` builds the two logistic models. Its error says so.
- **The GRAMMAR-Gamma approximation**, which both mixed models can use:
  `gwas-logistic`, because it is one item over both and building it once is
  cheaper than building the linear half of it here.
- **What it costs.** "Speed" of the spec asks for the linear model against
  plink2's 0.10 s and the mixed ones against GMMAT's 1.6 s and 2.2 s, over
  100000 variants x 1000 individuals. The owner decided on 23 September 2026
  that the measurements of a plan are made in a session of their own once it
  is merged, with the `performance-review` skill, and not as one of its work
  packages.
- **The kinship** the linear mixed model takes: the plan `kinship`.

The spec's one open point, the variant that separates the cases from the
controls, belongs to the logistic Wald test and changes no task of this
plan.

## What has to be in place

The work packages 1 and 2 of the plan `kinship`, for the matrix the linear
mixed model takes and for the row pass its work package 1 moved. The linear
algebra this plan calls, the Cholesky with its solve and log determinant,
the thin QR, the rank and the eigendecomposition, is on `main` already.

The counts below were measured on 23 September 2026 on `main` with
`spec/gwas` merged in, before any kinship work existed. What this plan
actually starts from, `b387def` of `plan/kinship`, gives 637 tests in the
core crate with 2 ignored, 149 in the linear algebra crate, 369 pytest and
253 node tests, higher as this section says they would be, with every other
check clean and `cargo test -p popnei --lib gwas -- --list` printing
`0 tests`:

- `cargo fmt --all --check` and `cargo clippy --workspace --all-targets --
  -D warnings`: clean.
- `cargo test --workspace`: 604 in the core crate, 2 ignored, 149 in the
  linear algebra crate, 136 of them on faer.
- `cargo test -p popnei --no-default-features`, the core crate on faer,
  which is the backend the wasm build uses and so what runs in a browser.
  It was in no check list until 23 September 2026 and is in the `coding`
  skill now. Any tolerance this plan adds is chosen against both backends:
  faer sits about seven times further from plink2 than Accelerate does on
  the same data, which the order and the blocking of the sums allow and
  which a bound fixed on Accelerate alone would fail under wasm.
- `cargo test -p popnei --lib gwas -- --list`: `0 tests`.
- `uv run maturin develop && uv run pytest`: 347 passed, 0 failed.
- `npm run build && npm test` in `js/popnei`: 242 pass, 0 fail.

Every layer exists, so every check of the `coding` skill runs. The two notes
of `docs/plans/linalg-gwas.md` about the release build and `wasm-bindgen`
still hold.

The reference data is in the repository, committed with the spec:
`tests/reference/gwas/`, holding `phenotypes.csv`, `causal_vars.csv`, the
plink2, GMMAT and rrBLUP outputs, and `make_reference.py`, which runs the
four programs again on popnei's own VCFs. It was run on 23 September 2026
and every number it produced matched what pyNei stored, to the bit, over 58
numeric columns.

## Work package 1: the two distributions

### What it gives

Nothing to a user directly: `chi2_sf_1df` and `t_sf_two_sided` are public in
the core crate so that their numbers can be checked where they can be seen,
and every p-value of both plans comes out of one of them. It stops inside
the core crate because the two are not part of the Python or the TypeScript
API; the comparison with pyNei is made at `calc_gwas`, in work package 3.

It is first because everything else ends in it, because it is the only work
package checked against literals alone with no reader and no reference file,
and because `libm` has to build for both wasm targets inside this workspace
before anything sits on it.

### Its deliverables

1. `libm` is a dependency of the core crate and both wasm targets build. The
   check: `cargo wasm-check` is clean, which covers both targets with the
   warnings denied, where today the crate has no `libm`.
2. The two functions are scipy's. The check: `cargo test -p popnei --lib
   distributions -- --list` names the three tests of "How it is verified" of
   "The two distributions", where today it prints `0 tests`, and they assert
   the incomplete beta at the four pairs within 1e-12 absolute,
   `t_sf_two_sided` at 5, 17 and 197 degrees of freedom within 1e-10
   relative, and `chi2_sf_1df` within 1e-12 relative.
3. The literals came from scipy and can be got again. The check: a script
   under `tests/reference/gwas/` prints them from scipy 1.18.1 and the test
   file holds what it printed.

### What it stands on

Nothing of this plan.

### Its tasks

- [x] 1.1 `libm` in the workspace and in the core crate, `chi2_sf_1df` in a
      new `crates/popnei/src/gwas.rs`, and the wasm checks. Built from "The
      two distributions" of `docs/specs/gwas.md`. Serves deliverables 1 and
      2. Needs nothing.
- [x] 1.2 The regularized incomplete beta and `t_sf_two_sided`, from the
      recurrence the spec writes out, with the script that prints scipy's
      numbers and the three tests. Built from the same section. Serves
      deliverables 2 and 3. Needs 1.1.

### What could go wrong

The continued fraction is transcribed in the spec. The pair `(98.5, 0.5)` is
the one a t of 197 degrees of freedom uses and the one closest to the panel,
so a failure there and not at `(0.5, 0.5)` is the fraction and not the front
factor.

This section said until 23 September 2026 that the two guards of the
fraction, the `tiny` that keeps a denominator of 0 from dividing and the
`eps` that stops it, are what make it converge, and that dropping either
gives numbers that are right for most arguments and wrong for some. The
review of this work package measured that and it is false for every argument
either plan can reach, so a test writer who believed it would hunt for a
case that does not exist. The first denominator is bounded below by
`2 / (a + b + 2)` in both branches, so with `b` of 1 / 2 the `tiny` can fire
only for a panel of about 4e300 individuals; one reviewer saw a minimum of
4.0276e-6 over 6009003 calls, which is that bound at 1e6 degrees of freedom,
and another 4.06e-5 over 116802. The `eps` caps the work and not the digits:
with it turned off the fraction runs its 500 rounds and the worst value
moves by 2.3e-13 relative, nothing becoming non-finite. So neither guard can
be caught by a test on a value, and the only assertion that could fail is
one on the number of rounds, which is at most 52 over the wide sweep against
the 500 allowed. Both guards stay, because the recipe and pyNei have them
and because a caller with another `b` would need the `tiny`.

What the review found instead, in the place this section was pointing away
from: `t_sf_two_sided` cancelled `1 - x` out of an `x` that had rounded to
1, and lost up to eight digits for a `t` near 0.

## Work package 2: what every model shares

### What it gives

Nothing a user can call yet: the individuals that are tested, the design and
its refusals, the dosages of a block and the shape of the result, which
every one of the four models sits on. It stops inside the core crate because
there is no model to run yet; the first Python function is work package 3,
which is also where the comparison with pyNei is first made.

### Its deliverables

1. The individuals tested are those with a phenotype, in the order the
   variants have them, and the five refusals of "Which individuals are
   tested, and the design" each raise. The check: `cargo test -p popnei
   --lib gwas::design -- --list` names one test per refusal, where today it
   prints `0 tests`, and one for the order.
2. The design is refused when its columns are not independent, at numpy's
   tolerance. The check: a cargo test on a design with a covariate that is
   twice another gets the error, and one on a design whose smallest singular
   value is 1e-11 of its largest does not, which is the band
   `docs/specs/linalg.md` measured.
3. The dosages of a block are computed over the tested individuals alone.
   The check: a cargo test where the tested individuals are half the panel
   gets an `allele_freq` that differs from the whole panel's.
4. Which model and which test are chosen, and the two combinations that
   raise. The check: `cargo test -p popnei --lib gwas::choice -- --list`
   names a test for each of the four models' defaults and for the two
   refusals.

### What it stands on

Work package 1, for nothing but the module it lives in.

This section said until 23 September 2026 that the dosages use the row pass
that work package 1 of the plan `kinship` moved into `variant`. They cannot.
That pass always divides the centered dosages by a scale, either of the
dosages themselves or of Hardy Weinberg, because the kinship and the
principal components want a variant standardized; a study wants the dosage
itself, since `beta` is the effect of one copy of an allele in the units of
the trait. The worked example of the spec has the dosages 0, 1, 2 and a
`beta` of 1.5, which the scale would turn into 1.2247. The pass also gives
no mean back, and a study reports that mean as the `allele_freq` of every
variant, including the ones it cannot test. So work package 2 has its own
row, which calls the two vectorized passes of `variant` that do apply,
`count_alleles` and `the_codes_of_the_genotypes`, and adds the mean and the
fill for a missing genotype.

`ld.rs` is the precedent for a module having its own dosage rule, and for
that much only: it reads its rows one after another and has no rayon, so it
is no precedent for the drive over the rows and its wasm twin, which the
review of this work package measured as the larger half of about 120
duplicated lines. Whether the two rows become one, with the scale made
optional so that a study can ask for the dosage itself, is for the owner at
the end of this plan: the change is in `variant.rs`, which another plan owns
while this one runs, and `gwas-logistic` would be its third caller.

### Its tasks

- [x] 2.1 The tested individuals, the design, its refusals and the rank
      check, in `crates/popnei/src/gwas.rs`. Built from "Which individuals
      are tested, and the design" of `docs/specs/gwas.md`. Serves
      deliverables 1 and 2. Needs 1.1.
- [x] 2.2 The dosages of a block over the tested individuals, the shape of
      the result, the variants that have no answer, and the choice of model
      and test. Built from "What it gives", "The variants that have no
      answer" and "The Rust interface". Serves deliverables 3 and 4. Needs
      2.1.

### What could go wrong

The dosages, the mean, the frequency and whether a variant varies are all
over the tested individuals and not over the panel, which matters as soon as
a phenotype leaves anyone out and which no test with a complete phenotype
catches. Deliverable 3 is the one that does.

## Work package 3: the linear model

### What it gives

`calc_gwas(variants, phenotype, trait="continuous", covariates=...)` in
Python and `calcGwas` in TypeScript, testing every variant against a
continuous trait with covariates and no kinship, with the effect, its
standard error and its p-value.

### Its deliverables

1. The whole study is plink2's. The check: a pytest test reads
   `tests/reference/kinship/panel_called.vcf.gz`, runs `calc_gwas` with
   `cov1` and `cov2`, and over all 1200 variants agrees with
   `plink2.panel_called.glm.linear.tsv` as "How it is verified" of "The
   linear model" of the spec asks: `allele_freq` within 1e-6 absolute, and
   `beta`, `se` and `p_value` within 1e-5 relative. plink2 prints that file
   to six significant digits, which round a value by up to 5e-6 of itself,
   so this check has at most twofold headroom and says that popnei computes
   the same quantity. Deliverables 2 and 3, the worked example at 1e-12
   relative and pyNei at 1e-9 relative, are the ones that would catch a
   wrong digit.
2. The worked example and the six literals are cargo tests. The check:
   `cargo test -p popnei --lib gwas::lm -- --list` names them, where today
   it prints `0 tests`; the worked example of "The worked example" asserts
   the null model's two coefficients, its residual sum of squares and the
   three rows within 1e-12 relative, and needs no reference file.
3. popnei and pyNei agree. The check: a pytest test runs both on the panel
   and `beta`, `se` and `p_value` agree within 1e-9 relative and the NaN
   variants are the same.

   The 1e-9 of deliverable 3 and the 1e-12 of deliverable 2 are where to
   start and not where to stop, as "How it is verified" of "What every
   model shares" of the spec now says. Each is lowered until it fails, set
   two or three times above where it broke, and both numbers go in the
   report. The spec's reason: the kinship's matrix matched plink2's binary
   output to 4.44e-16 absolute, which reads as a wide margin, while its
   worst entry as a ratio was 3.31e-13, so its 1e-12 relative bound had two
   to three times the worst case and not the thousandfold the absolute
   figure suggested.
4. The block size changes nothing. The check: a pytest test reads the same
   panel in blocks of 77 and gets `stats` equal within 1e-12 relative.
5. `calcGwas` under node gives the same numbers. The check: `npm test` in
   `js/popnei` asserts the six literals and the worked example.

### What it stands on

Work packages 1 and 2, whole, since the work packages run in order.

### Its tasks

- [ ] 3.1 The linear model's null fit and its test, in the core, with the
      worked example and the six plink2 literals as cargo tests. Built from
      "The linear model" and "The worked example" of `docs/specs/gwas.md`.
      Serves deliverable 2. Needs 2.2.
- [ ] 3.2 The Python function: the binding, `GWASResult`, `NullModel`, the
      three enums, and the pytest tests against plink2, against pyNei and
      for the block size. Built from "Its Python function, and its
      TypeScript one" and "How it is verified" of "What every model shares".
      Serves deliverables 1, 3 and 4. Needs 3.1.
- [ ] 3.3 The TypeScript function: the binding, the result object and the
      node test. Built from the same section. Serves deliverable 5. Needs
      3.1, and it can run beside 3.2.

### What could go wrong

This is where a misreading of the shared parts of work package 2 first meets
a reference program, so a failure here is as likely to be in the design or
the dosages as in the linear model. The order of the tested individuals is
the likeliest: the phenotype of the panel is given in the order the
variants have them, so a study that sorted them would pass every test of
this work package and fail in `gwas-logistic`, where the phenotype is
reversed. Deliverable 1 of work package 2 is what guards it, and it is a
cargo test and not a pytest one for that reason.

## Work package 4: the linear mixed model

### What it gives

The same function with a kinship, `calc_gwas(..., kinship=k)` and
`test="wald"` or `"score"`, which accounts for the relatedness of the panel
so that a variant that only marks ancestry does not look associated.

### Its deliverables

1. The null model is GMMAT's. The check: a pytest test gets
   `genetic_variance` 1.221617, `residual_variance` 0.342359 and the three
   covariate effects 4.678021, 0.473361 and 1.110279 within 1e-5 absolute,
   from `tests/reference/gwas/gmmat.null_models.tsv`.
2. The fit is at its optimum. The check: a cargo test asserts that `y' p y`
   on the panel is 197 within 1e-6. "How it is verified" says this one is
   made at the private function that fits this null and pins it.
3. The Wald test is rrBLUP's. The check: a pytest test with `cov2` alone and
   the kinship gets `-log10(p_value)` within 1e-4 of
   `rrblup.panel_called.lmm.tsv` over all 1200 variants, and a cargo test
   asserts the six literals.
4. The score test is GMMAT's, on both panels. The check: a pytest test gets
   `1 / se**2` within 1e-5 relative of GMMAT's `VAR` and
   `|log10(p / p_GMMAT)|` below 1e-4 over all 1200 variants of
   `gmmat.panel_called.lmm.score.tsv` and `gmmat.panel.lmm.score.tsv`, and a
   cargo test asserts the six literals of both panels.
5. The study finds what was planted. The check: a pytest test gets at least
   3 of the 5 variants of `causal_vars.csv` among the 10 smallest p-values.
6. popnei and pyNei agree, and TypeScript gives the same numbers. The check:
   the pytest comparison of work package 3 with a kinship, and `npm test`
   asserting the six score test literals. The comparison with pyNei is
   lowered until it fails and set two or three times above, as in
   deliverable 3 of work package 3, and the report carries both numbers.

### What it stands on

Work package 3, and the plan `kinship` for the matrix.

### Its tasks

- [ ] 4.1 The REML search and the null fit: the eigendecomposition, the
      clamp at 0, the 101 points and the 60 golden section steps, the
      criterion and the two variances, with the cargo test of `y' p y`.
      Built from "What it gives" of "The linear mixed model". Serves
      deliverables 1 and 2. Needs 3.1.
- [ ] 4.2 The projection matrix and the two tests, with the six literals of
      each as cargo tests. Built from the same section. Serves deliverables
      3 and 4. Needs 4.1.
- [ ] 4.3 The `kinship` and `test` arguments through both bindings and both
      packages, with the pytest tests against rrBLUP, GMMAT and pyNei, the
      causal variants, and the node test. Built from "Its Python function,
      and its TypeScript one" of "What every model shares". Serves
      deliverables 1, 3, 4, 5 and 6. Needs 4.2, 3.2 and 3.3.

### What could go wrong

The REML search has to be reproduced step for step or the variance
components move, and the spec lists the six things that have to match. A
`genetic_variance` near but not at 1.221617 means the search and not the
criterion; one far from it means the criterion or the clamp.

The two tests share `num` and `den` and differ only in the standard error
and the distribution, so a failure in one and not the other is in that
line and not in the projection matrix.

rrBLUP is given `cov2` alone, because it takes every fixed effect as a
factor, and GMMAT is given both covariates. A run that gives both to the
rrBLUP comparison gets numbers that are close and not equal, which reads
like a tolerance problem and is not.

## How the whole plan is checked

The sum of the work packages, and `cargo wasm-check` clean with `npm run
build && npm test` passing in `js/popnei`, so that the linear half of the
association study reaches a browser.
