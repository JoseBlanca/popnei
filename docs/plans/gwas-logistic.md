# Plan: the association study of a binomial trait

23 September 2026. State: under way. The owner approved it on 23
September 2026 and it is being carried out on the branch
`plan/gwas-logistic`. It builds
the rest of `docs/specs/gwas.md`: the logistic model, the logistic mixed
model, the GRAMMAR-Gamma approximation that both mixed models can use, and
the measurements. It is the last of three plans; `kinship` and `gwas-linear`
come before it and it finishes the module.

It is carried out in the worktree `.claude/worktrees/gwas-logistic` on the
branch `plan/gwas-logistic`, which starts from `main` with the plans
`kinship` and `gwas-linear` merged into it. The report is
`docs/reports/gwas-logistic.md`.

The work packages run in order, and so do the tasks inside each one,
except where a task says it can run beside another. What a task says it
needs is what has to be committed before it starts.

## In and out

Built, through the core crate, both binding crates and both packages:

- The logistic model, with its Wald test and its score test, checked against
  plink2's `--glm` and against R's `anova(glm, test = "Rao")`.
- The logistic mixed model with its score test, fitted the way
  `docs/reports/glmm-method/README.md` measured, checked against GMMAT.
- `use_grammar_gamma_approx`, on both mixed models.

Not built, with where it goes:

- **A Firth penalized regression** for the variant that separates the cases
  from the controls, unless the owner chooses it: the spec's Open 1, below.
- **Multiple testing, a joint model of a multiallelic variant, a sparse
  kinship and interactions**: "Not in this spec" of the spec, none of them
  in popnei.
- **What it costs.** "Speed" of the spec asks for the four models over
  100000 variants x 1000 individuals, for the null fit of the logistic mixed
  model at 1000, 2000 and 4000 individuals against the 0.095, 0.608 and
  5.285 seconds that `docs/reports/glmm-method/README.md` measured in numpy,
  and for the per variant test with and without the approximation. The owner
  decided on 23 September 2026 that the measurements of a plan are made in a
  session of their own once it is merged, with the `performance-review`
  skill, and not as one of its work packages. That session is also where the
  ratio of 1.88 to 2.06 is checked in Rust, which the report measured
  between two numpy programs and could not say would hold.

Two open points of the spec change work of this plan.

**Open 1 of `docs/specs/gwas.md`, the variant that separates the cases from
the controls.** Meanwhile, task 1.2 gives it NaN for its effect, its
standard error and its p-value with no reason, which is pyNei's behaviour
and what the plink2 comparison of deliverable 2 is written against. If the
owner chooses NaN with a reason, a task is added to work package 1 for the
field and its column in both packages, and no number moves. If they choose
the Firth regression, that is a work package of its own and a change to the
spec first, and the plink2 comparison gains the variant it now leaves out.

**Open 2 of `docs/specs/gwas.md`, a variant there is nothing left to test.**
It is one rule in four places, and two of them are tests of this plan: the
score test of the logistic model in task 1.1 and the score test of the
logistic mixed model in task 2.3. Meanwhile, both refuse a variant of which
the design leaves at most the tested individuals times 2.2e-16 of what there
was, and give it the three NaNs a variant with no variance gets. No literal
of either panel is near that threshold, so nothing in the deliverables moves
either way. If the owner chooses to form the residual exactly instead, both
tasks lose the comparison and gain the exact quantity, which for these two
is a product with the projection matrix per variant.

## What has to be in place

The plans `kinship` and `gwas-linear` merged into `main`; the linear algebra
is there already. From `gwas-linear` this plan takes the two distributions,
the design, the dosages of a block, the result and the Python and TypeScript
functions it extends, and from `kinship` the matrix the mixed model takes.

Every layer exists, so every check of the `coding` skill runs, including
`cargo test -p popnei --no-default-features`, the core crate on faer, which
is the backend the wasm build uses and so what runs in a browser. It was in
no check list until 23 September 2026. Any tolerance this plan adds is
chosen against both backends: faer sits about seven times further from
plink2 than Accelerate does on the same data, which the order and the
blocking of the sums allow and which a bound fixed on Accelerate alone would
fail under wasm.

The two notes of `docs/plans/linalg-gwas.md` about the release build and
`wasm-bindgen` still hold.

The reference data is in the repository: `tests/reference/gwas/`, whose
`make_reference.py` ran the four programs again on popnei's own VCFs on 23
September 2026 and matched what pyNei stored to the bit.

## Work package 1: the logistic model

### What it gives

`calc_gwas(variants, phenotype, trait="binomial", covariates=...)`, testing
every variant against a 0/1 trait with no kinship, giving the effect as a
log odds ratio. Both its tests: the Wald one by default, and the score one
on request.

### Its deliverables

1. The null fit and the score test are R's. The check: a pytest test with
   `test="score"` on `tests/reference/kinship/panel_called.vcf.gz` gets
   `(beta / se)**2` within 1e-2 absolute and `|log10(p / p_R)|` below 1e-3
   of `r.panel_called.glm.score.tsv` over all 1200 variants, and a cargo
   test asserts the six literals within 1e-3 and 1e-3.
2. The Wald test is plink2's. The check: a pytest test gets `beta` within
   1e-4 times the `se` of that variant, `se` within 5e-4 times it, and
   `p_value` within 5e-3 relative of
   `plink2.panel_called.glm.logistic.hybrid.tsv` over the 1199 variants
   plink2 did not fall back to Firth for, and a cargo test asserts the six
   literals within 1e-5, 1e-4 and 5e-3.
3. The variant that runs away is exactly the one plink2 marked. The check: a
   pytest test asserts that the variants whose `p_value` is NaN are exactly
   the rows with `FIRTH?` equal to `Y`, which is `var0006` and no other.
4. popnei and pyNei agree, and TypeScript gives the same numbers. The check:
   a pytest test runs both on the panel with each test and gets `beta`
   within 1e-13 of the `se` of that variant, `se` and `p_value` within 1e-13
   relative, and the same NaN variants; and `npm test` asserts the six Wald
   literals and the six score ones. The bound against pyNei is per model
   since 24 September 2026, with 1e-9 its ceiling and not its value, and
   1e-13 is 2.2 times where this model breaks. What it has to clear is the
   p-value of `var1004` of `panel_called` under the Wald test on Accelerate,
   4.524e-14, and not the effect, whose worst is 1.151e-14 at `var0197` of
   `panel_called` under the score test on faer.

### What it stands on

The plans `gwas-linear` and `kinship`, merged.

### Its tasks

- [x] 1.1 The logistic null fit by iteratively reweighted least squares and
      the score test, in `crates/popnei/src/gwas.rs`, with the six R
      literals as cargo tests. Built from "The logistic model" of
      `docs/specs/gwas.md`. Serves deliverable 1. Needs nothing of this
      plan.
- [x] 1.2 The per variant Wald fit, its three marks of a runaway and the
      NaN they give, with the six plink2 literals as cargo tests. It is its
      own task and its own commit because a variant wrongly marked as
      running away loses its p-value in silence, and deliverable 3 is what
      guards it. Built from the same section. Serves deliverables 2 and 3.
      Needs 1.1.
- [x] 1.3 The binomial trait through both bindings and both packages: the
      `trait` argument reaching the two new models, the error that
      `gwas-linear` left for a binomial trait now gone, and the pytest and
      node tests. Built from "Its Python function, and its TypeScript one"
      of "What every model shares". Serves deliverables 1, 2, 3 and 4.
      Needs 1.2.

### What could go wrong

The Wald test fits one logistic regression per variant and the spec's three
marks of a runaway are what keep a fit that will not settle from producing a
number instead of a NaN. The first of the three is tested as a value that is
**not finite** and not as an infinity, which is what makes the two backends
of the linear algebra crate agree: one gives an infinity where the other
gives a NaN for the same input. A test for an infinity passes natively and
fails under wasm, and the node test of deliverable 4 is where that shows.

plink2 stops its logistic fit earlier than popnei does, which is why the
p-value tolerance is 5e-3 and not 1e-5. A run that is within 1e-5 of plink2
on the p-values has probably stopped early too, and its `beta` will be
further off, not closer.

## Work package 2: the logistic mixed model

### What it gives

`calc_gwas(..., trait="binomial", kinship=k)`, the fourth model, which
accounts for the relatedness of the panel in a 0/1 trait. Its only test is
the score test, and asking for the Wald one raises, which `gwas-linear`
already built.

### Its deliverables

1. The null model is GMMAT's. The check: a pytest test gets
   `genetic_variance` 1.508057 within 1e-5 absolute, `residual_variance`
   `None`, and the covariate effects -1.416464, 0.753476 and 1.583210,
   from `tests/reference/gwas/gmmat.null_models.tsv`.
2. The score test is GMMAT's, on both panels. The check: a pytest test gets
   `1 / se**2` within 1e-5 relative of GMMAT's `VAR` and
   `|log10(p / p_GMMAT)|` below 1e-4 over all 1200 variants of
   `gmmat.panel_called.glmm.score.tsv` and `gmmat.panel.glmm.score.tsv`, and
   a cargo test asserts the six literals of both panels.
3. The fit never forms an inverse while it iterates. The check: a cargo test
   counts the calls the fit makes to the inverse of a factorized matrix on
   the panel and gets 1, where the number of Cholesky factorizations is 22,
   which is what "How popnei fits it, and why not pyNei's way" describes.
4. popnei and pyNei agree, and TypeScript gives the same numbers. The check:
   a pytest test runs both on both panels with a kinship and gets agreement
   within 1e-9 relative, and `npm test` asserts the six literals.

### What it stands on

Work package 1, whole, since the work packages run in order.

### Its tasks

- [ ] 2.1 The linearization for one value of the variance component: the
      working trait, the weights, the covariance factored with a Cholesky
      and applied by solving, and its stopping rule. Built from "What it
      gives" and "How popnei fits it, and why not pyNei's way" of "The
      logistic mixed model", and from
      `docs/reports/glmm-method/README.md` for why it is a Cholesky. Serves
      deliverable 3. Needs 1.1.
- [ ] 2.2 The step on the variance component: the trace from the identity,
      through the triangular solve reading the lower half, the average
      information, the bracket and the stopping rule, with the cargo test
      that counts the inverses. It is its own task and its own commit
      because a wrong trace moves the variance component and nothing
      crashes; deliverables 1 and 2 are what guard it. Built from the same
      sections. Serves deliverables 1 and 3. Needs 2.1.
- [ ] 2.3 The score test on the projection matrix, the six literals of both
      panels as cargo tests, and the model through both bindings and both
      packages with its pytest and node tests. Built from "What it gives"
      and "How it is verified" of that item. Serves deliverables 1, 2 and 4.
      Needs 2.2 and 1.3.

### What could go wrong

This is the piece that popnei does differently from pyNei, so a
disagreement here is not a transcription slip in the way the other three
models' would be. `docs/reports/glmm-method/README.md` measured the two fits
to agree to 2.9e-15 in the variance component and to give the same largest
`|log10(p / p_GMMAT)|` to every digit, 8.497e-06, so a run that is close but
not that close is the trace identity or the Cholesky and not the algorithm.
The identity divides by the variance component, which is 0 at the boundary
where the kinship explains nothing, and that case has to be taken before the
division and not after.

The fit starts from the plain logistic null of work package 1 and not from
anywhere else, and its own intercept starts at
`log((m + 1e-6) / (1 - m + 1e-6))`. A fit started elsewhere walks a
different path through the bracket and can stop at another variance
component, which reads as a wrong trace and is not.

## Work package 3: the GRAMMAR-Gamma approximation

### What it gives

`use_grammar_gamma_approx=True`, which makes the per variant work of both
mixed models linear in the individuals instead of quadratic, at a cost in
accuracy that grows with how strongly the panel is structured.

### Its deliverables

1. The second pass estimates the factor from the first variants that vary.
   The check: a cargo test on the panel gets a factor from 100 variants, and
   one on a first block with fewer than 100 that vary gets it from those.
2. The approximation is close to the exact answer. The check: a pytest test
   on the panel with the linear mixed model gets the median of
   `log10(p_approx / p_exact)` within 0.1 of 0 and the largest within 1.5,
   and `beta` within 0.5 relative, which are pyNei's own numbers in
   "How it is verified" of that item.
3. Asking for it without a kinship raises, and the result says it was used.
   The check: a pytest test for the error and for
   `used_grammar_gamma_approx`.
4. TypeScript has it. The check: `npm test` runs the approximation and gets
   the same relation to the exact answer.

### What it stands on

Work package 2, and work package 4 of the plan `gwas-linear` for the linear
mixed model it is also checked on.

### Its tasks

- [ ] 3.1 The second pass, the factor, and the approximate denominator in
      both mixed models, in the core, with the cargo tests. Built from "The
      GRAMMAR-Gamma approximation" of `docs/specs/gwas.md`. Serves
      deliverable 1. Needs 2.2.
- [ ] 3.2 The argument through both bindings and both packages, with the
      pytest tests of the relation to the exact answer and of the refusal,
      and the node test. Built from the same item. Serves deliverables 2, 3
      and 4. Needs 3.1 and 2.3.

### What could go wrong

There is no program outside the project to check this against, so the only
check is the relation to popnei's own exact answer, and it is a loose one: a
p-value may be out by a factor of 30 and still pass. An approximation that
is wrong in a way that does not grow with the structure of the panel would
pass it. What would catch that is the factor itself, which is a mean of
ratios that are all near each other on a panel like this one, so a factor
far from the ratios it was averaged from is the sign to look for.

## How the whole plan is checked

The sum of the work packages, and three things more. `cargo wasm-check`
clean and `npm run build && npm test` passing in `js/popnei`, so the whole
association study reaches a browser. The pytest suite run against pyNei on
both panels for all four models at once, which is the comparison
`docs/objectives.md` asks for. And `tests/reference/gwas/make_reference.py`
run once more at the end, so that the four reference programs are known to
still give what the literals say on the machine the work was done on.
