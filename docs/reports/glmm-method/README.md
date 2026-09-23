# A cheaper fit for the null of the logistic mixed model

23 September 2026. The logistic mixed model is the fourth of the four
models of the association study, a 0/1 trait with the kinship as a random
effect, and fitting its null model is the most expensive thing pyNei's
`calc_gwas` does: 25 s at 5000 individuals, from section 5 of
`docs/rust_core.md`, which leaves it as an open question. The owner asked
on 23 September 2026, before `docs/specs/gwas.md` was written, whether a
cheaper fit could be found, because a different fit would want different
operations of `crates/popnei-linalg`, whose spec is being built from right
now.

This says what was tried and what came out. The same fit runs about twice
as fast and gives the same numbers, by never forming the inverse of an
individuals x individuals matrix while it iterates; an order of magnitude
is not there for a kinship that is dense, which popnei's is; and taking
the cheaper route asks one small thing of `docs/specs/linalg.md` that is
not among the seven operations it is being built from, which are a
Cholesky factorization and the solve, log determinant and inverse that
come off it, a thin QR, a solve against an upper triangular matrix and the
rank of a matrix.

What is recommended, and what the owner decides, is whether
`docs/specs/gwas.md` is written around the cheaper fit. Nothing else here
is a decision.

Every number was measured on the owner's Apple M5 Pro with 18 cores, with
numpy 2.5.3 on Accelerate and its threads, which is what pyNei runs on.
The scripts are beside this file and "How to run these again" says how.

## What the fit does

The null model is the model of the trait fitted once, with the covariates
in it and no variant. For a 0/1 trait with a kinship, the model has a
random effect whose covariance is the kinship, and the number to estimate
beside the covariate effects is **tau**, the variance of that effect.

pyNei fits it as GMMAT does, the R package of Chen and others 2016 that
`docs/objectives.md` names as the reference outside the project for a
mixed model association study, and that every number here is checked
against. The method is **penalized quasi-likelihood** (Breslow and Clayton
1993): a 0/1 trait is turned into a continuous **working trait**,
each individual carrying a weight that says how much its 0 or 1 tells us
at the fit so far, and a weighted linear mixed model is fitted to that
working trait. The working trait and the weights are then made again from
the new fit, and so on. One pass of that is a **linearization**.

tau is estimated by **restricted maximum likelihood**, which is maximum
likelihood on the part of the data that the covariates cannot explain, so
that fitting the covariates does not pull the variances down. pyNei takes
Newton steps on tau, each one from the first derivative of that likelihood
and an approximation of the second, the **average information**, which is
what `_GLMMNull` of `src/pynei/gwas.py` at commit ef0ca6e does.

Both quantities the step needs are built from the **projection matrix** P,
individuals x individuals: with `d` the design, the covariates with a
column of ones, and `sigma` the covariance of the working trait,

    sigma = tau * k + w⁻¹          k the kinship, w⁻¹ the weights on the diagonal
    p     = sigma⁻¹ − sigma⁻¹ d (d' sigma⁻¹ d)⁻¹ d' sigma⁻¹

P also outlives the fit. A variant is tested with a **score test**, which
asks how far from 0 the slope of the trait on that variant is, measured at
the null model and without fitting the model again with the variant in it;
it is the test used here because fitting a mixed model per variant would
cost a whole fit per variant. It is built from `x' p y` and `x' p x` for
the dosages `x` of that variant, so P is wanted as a matrix once the fit
is done.

## Where the time goes

`count_iters.py` counts the two loops. Over 200, 500 and 1000 individuals
the fit takes 7 to 8 steps on tau and 22 to 26 linearizations, and the
counts barely move with the individuals; `final.py` gives 21 to 26
linearizations up to 4000. Every linearization forms `sigma` and inverts
it, at line 689 of `gwas.py`.

What the individuals x individuals operations cost, from `primitives.py`
and `triangular_and_matvec.py`, in seconds:

| individuals | inverse | Cholesky | inverse of the triangular factor | eigenvalues and vectors | eigenvalues alone | one product with a vector |
|---|---|---|---|---|---|---|
| 1000 | 0.0079 | 0.0018 | 0.0018 | 0.047 | 0.027 | 0.00002 |
| 2000 | 0.054 | 0.0138 | 0.0116 | 0.330 | 0.196 | 0.00024 |
| 4000 | 0.376 | 0.126 | 0.087 | 3.33 | 2.44 | 0.00063 |

A Cholesky factorization, which exists for a symmetric matrix that is
positive definite and which every `sigma` here is, costs a third of an
inverse. That is the whole opportunity: the fit inverts because it wrote
P out as a matrix, and it does not need P as a matrix while it iterates.

## What the cheaper fit changes

Two things, and nothing else. It is `fit_c` of `fits.py`, against pyNei's
`fit_a`.

**It factors `sigma` instead of inverting it.** Everything a linearization
does with `sigma⁻¹` is applied to the design, `sigma⁻¹ d`, which is one
column for each coefficient, or to a vector. Each of those is a solve
against the Cholesky factorization, and the factorization is a third of
the inverse.

**It takes the one quantity that seemed to need the whole inverse from an
identity.** The step on tau needs the trace of `p k`, which is

    trace(p k) = trace(sigma⁻¹ k) − trace((d' sigma⁻¹ d)⁻¹ d' sigma⁻¹ k sigma⁻¹ d)

whose second term is a handful of columns and costs nothing. The first
term looks as though it needs every entry of `sigma⁻¹`. It does not,
because `tau k = sigma − w⁻¹`, so

    tau * sigma⁻¹ k = i − sigma⁻¹ w⁻¹     and     trace(sigma⁻¹ k) = (n − trace(sigma⁻¹ w⁻¹)) / tau

with `n` the individuals. And `trace(sigma⁻¹ w⁻¹)` is the sum of the
squares of the entries of `l⁻¹ w^-1/2`, the Cholesky factor inverted and
each of its columns divided by the square root of that individual's
weight: one triangular matrix inverted, or one solve against a triangular
matrix with one right hand side for each individual. It is wanted once for
each step on tau, 7 to 9 times over a fit, and not once per
linearization.

The inverse is formed once, at the end, because the score test wants P as
a matrix.

## What it gives

`final.py`, on simulated panels of families of four full sibs inside three
subpops with an Fst of 0.1 and 2000 variants, which is how pyNei's
reference panel is made. Seconds for the null fit, and beside it the score
test of 100000 variants that the fit feeds, which is the dosages of each
block multiplied by P:

| individuals | pyNei's fit | the cheaper fit | ratio | the score test of 100000 variants | linearizations |
|---|---|---|---|---|---|
| 500 | 0.049 | 0.024 | 2.03 | 0.11 | 26 |
| 1000 | 0.195 | 0.095 | 2.06 | 0.32 | 22 |
| 2000 | 1.199 | 0.608 | 1.97 | 1.26 | 21 |
| 4000 | 9.959 | 5.285 | 1.88 | 5.13 | 25 |

The fit and the test are of the same order once the individuals reach a
few thousand, so halving the fit takes about a quarter off the study.

## It changes no number

The two fits take the same steps in the same order and stop at the same
place; only the arithmetic of each step differs. Measured:

- On pyNei's reference panel, 200 individuals and 1200 variants, tau
  agrees to 2.9e-15 relative, the covariate effects to 3.3e-15 and P to
  1.4e-15. On simulated panels of 500 and 1000 individuals, tau agrees to
  3.5e-15.
- The p-values of the score test of all 1200 variants of the reference
  panel, against GMMAT 1.5.0's `glmm.score`, `check_pvalues.py`: pyNei's
  fit and the cheaper one both give a largest `|log10(p / p_GMMAT)|` of
  8.497e-06 and a largest relative difference of the variance of the score
  of 5.342e-06. The test of `test_gwas.py` in pyNei accepts 1e-4 on the
  first and 1e-5 on the second.
- Both fits give tau to 6.3e-06 of GMMAT's 1.508057 and the covariate
  effects to 1.1e-06 of GMMAT's.

## A kinship that is not positive semidefinite

A Cholesky factorization refuses a matrix that is not positive definite.
An LU factorization, which is what numpy's `inv` uses and which writes any
square matrix as a lower times an upper triangular one with its rows
exchanged, carries on and gives a number. "What the GWAS
calls of numpy, and where" of `docs/specs/linalg.md` raised this: pyNei
divides the kinship entry by entry by the variants called in both
individuals of each pair, which can leave a matrix with an eigenvalue
below 0, and asked whether the logistic mixed model would then stop where
pyNei gave a fit.

It does not, at any missing rate that was tried. `final.py` builds the
kinship of 400 individuals and 2000 variants with genotypes missing at
random and fits both ways:

| missing genotypes | smallest eigenvalue of the kinship | largest | pyNei's fit | the cheaper fit |
|---|---|---|---|---|
| 0 | -0.0000 | 32.5 | tau = 1.4107 | tau = 1.4107 |
| 3 in 100 | -0.0326 | 32.5 | 1.4402 | 1.4402 |
| 10 in 100 | -0.1173 | 32.5 | 1.4481 | 1.4481 |
| 25 in 100 | -0.3540 | 32.5 | 1.4795 | 1.4795 |
| 50 in 100 | -1.0593 | 32.7 | 1.2598 | 1.2598 |

The reason is in `sigma = tau * k + w⁻¹`. A weight is at most 0.25, so
`w⁻¹` puts at least 4 on every diagonal entry, and `sigma` only goes
indefinite once tau passes 4 divided by the size of the negative
eigenvalue: about 3.8 even at 50 genotypes missing in 100, where the tau
that is fitted is 1.26. So the fit would have to find a kinship effect
three times the largest seen here before it stopped.

## What it asks of the linear algebra crate

The seven operations of `docs/specs/linalg.md` cover all of it but one.
The trace needs the Cholesky factor inverted, or solved against with one
right hand side for each individual, and that factor is lower triangular
while `solve_upper_triangular` of that spec reads the upper half. Three
routes, at 4000 individuals, where a fit makes 9 steps on tau and so pays
the cost of one call 9 times:

| the trace is taken with | one call | over the fit | the fit | against pyNei's 9.959 s |
|---|---|---|---|---|
| the inverse of a triangular matrix, `dtrtri`, which the crate does not have | 0.087 s | 0.78 s | 5.29 s | 1.88 |
| a lower triangular form of `solve_upper_triangular`, a flag on a routine the crate already calls | 0.138 s | 1.24 s | 5.75 s | 1.73 |
| `solve_with_cholesky` as the spec has it, which is two triangular solves where one is wanted | 0.28 s | 2.5 s | 7.0 s | 1.42 |

The table of "What it gives" was measured on the first route, and the
other two rows are its 5.285 s with the 0.78 s of the trace replaced, from
the per call costs of `triangular_and_matvec.py`.

The second route is what is recommended. It asks the linear algebra crate
for nothing it does not already do, since `dtrtrs` of LAPACK and faer's
triangular solve both take which half to read as an argument, and it keeps
1.73 of the 1.88. The third needs no change to that spec at all and costs
0.31 of the ratio.

The inverse of a factorized matrix is still needed, once at the end of the
fit, to write P out for the score test.

## Measured and not taken

- **An eigendecomposition of the weighted kinship per linearization.**
  With the eigenvectors of the kinship with each of its rows and columns
  scaled by the square root of that individual's weight, `w^1/2 k w^1/2`,
  in hand, the whole search over
  tau costs one number per individual per evaluation, which is what makes
  the *linear* mixed model of pyNei cheap. It is not worth it here because
  the weights change at every linearization, so the eigendecomposition
  would have to be made again each time, and it costs 3.33 s at 4000
  individuals against an inverse's 0.376 and a Cholesky's 0.126. The
  number of linearizations would have to fall by a factor of 26 to pay for
  it, and it does not fall at all.
- **A conjugate gradient solve instead of a factorization**, which solves
  a system by products of the matrix with vectors and never factors it. It
  is what SAIGE does, the association study program written for biobanks
  of hundreds of thousands of individuals, where the kinship is kept
  sparse by setting the entries of distant pairs to 0, and it wins when
  the kinship is sparse. popnei's is dense:
  one product of `sigma` with a vector takes 0.63 ms at 4000 individuals,
  so a solve of 100 iterations is 63 ms for one right hand side, where one
  Cholesky is 126 ms and serves the three or four right hand sides a
  linearization has. It loses about twofold.
- **Estimating the trace from random vectors**, Hutchinson's estimator,
  which is the other half of SAIGE's answer. The trace then costs nothing,
  but the fit stops being the same calculation twice: with the 30 vectors
  SAIGE uses, tau carries about a per cent of noise, and the p-values are
  compared with GMMAT's to 1e-4 in `log10`, which is 0.02 per cent. It was
  not run, because the arithmetic above rules it out: even with a free
  trace the 25 Cholesky factorizations are 3.1 s of the 5.3 at 4000
  individuals, so it could not give more than a further 1.7.
- **Running the linearization loosely while tau is still moving**, and
  tightening it for the last steps, which leaves the answer where it was.
  It cuts the linearizations from 26 to 21 on 500 individuals, a fifth,
  and it is the scheme that `_GLMMNull`'s own comment says made tau cycle
  for ever when it was tried. pyNei keeps tau from cycling with a bracket:
  a tau whose derivative asks for a larger one and a tau whose derivative
  asks for a smaller one are remembered, the answer lies between them, and
  a Newton step that would leave that interval is replaced by the
  geometric mean of its two ends. The first version written here did not
  converge either, and it only did once the bracket was thrown away at the
  moment the tolerance tightens, because its two ends had been measured on
  linearizations that were not converged. A fifth is not worth a
  fit that can fail to converge.

## What is not known

- Nothing above was run in Rust. The ratios are between two numpy
  programs on the same BLAS, so they are ratios of the same library's
  work, but the part that is not BLAS, forming `sigma` and the products
  with vectors, is Python here and would be Rust there, which moves both
  sides by an unknown amount and the cheaper one less, since it does more
  of its work inside BLAS.
- Nothing was run above 4000 individuals. The ratio falls from 2.06 at
  1000 to 1.88 at 4000, and where it settles was not measured.
- The memory is unchanged and remains what `docs/rust_core.md` says: the
  fit holds several individuals x individuals matrices, 800 MB each at
  10000 individuals, so the browser caps the mixed models at a few
  thousand individuals whatever the arithmetic.

## How to run these again

The scripts need numpy, scipy and pyNei at the commit `pyproject.toml`
pins, and the last two read pyNei's reference panel from
`test/gwas_reference/` of a checkout of pyNei beside popnei:

    uv venv lab && uv pip install --python lab/bin/python numpy scipy pandas
    uv pip install --python lab/bin/python /path/to/pynei
    lab/bin/python primitives.py              # the table of operation costs
    lab/bin/python triangular_and_matvec.py   # the triangular routes and one product
    lab/bin/python count_iters.py             # the two loops of pyNei's fit
    lab/bin/python check_pvalues.py           # every p-value against GMMAT's
    lab/bin/python final.py                   # the table of fits, and missing genotypes

`fits.py` holds the three fits and the simulation, and is imported by the
last three.
