# A variant there is nothing left to test: the exact form, priced

25 September 2026. **Open 2** of `docs/specs/gwas.md` asks what a variant
should get when the subtraction that forms the linear mixed model's Wald
statistic has cancelled to rounding, and it recommends refusing such a
variant because the alternative, forming what the variant leaves of the
trait exactly, was priced at a matrix product for every variant. A reviewer
of `docs/reports/perf-gwas-2026-09-24.md` showed by algebra that it need not
be a product, and the same algebra is that report's finding H3, which
replaces the product of a block against the projection matrix with a
triangular solve against the covariance's own factor. This report has the
algebra checked, both halves built on the branch `gwas/exact-residual`, what
they cost, what they do to every column the two mixed models are verified
against, and a recommendation.

The recommendation is to take both halves and to leave Open 2's threshold
where it is. Nothing is merged into `main`.

## The words this document uses

The two **mixed models** of `docs/specs/gwas.md` are `lmm`, a continuous
trait with a kinship, and `glmm`, a binomial one. Each fits a null model
without any variant in it and then tests every variant through the
**projection matrix**

    p = v⁻¹ - v⁻¹ d (d' v⁻¹ d)⁻¹ d' v⁻¹

individuals by individuals, with `v` the covariance of the trait under the
null and `d` the design, which takes the covariates out of a variant and
weights it by that covariance. With `x` the dosages of one variant and `y`
the trait, `num` is `x' p y` and `den` is `x' p x`, and the effect of the
variant is `num / den` in both of the tests the models make. The **Wald
test** estimates the scale of the two variances again with the variant in
the model, so its `se` is built from `y' p y` less `num² / den`, what the
variant leaves of the trait; the **score test** holds them at the null and
needs no such quantity.

A **block** is the run of variants a reader gives at a time, 5000 variants
of 1000 individuals in every measurement here. Testing a block is one
product of the block, 5000 x 1000, with the projection matrix, 1000 x 1000.

The **factored form** is what the reviewer proposed. Write `l` for the lower
triangular Cholesky factor of the covariance, so that `v = l l'`, and `q`
for a set of directions of length 1 at right angles to each other spanning
what `l⁻¹` makes of the design, both worked out once for a study. Then
`p = m' m` with `m = (i - q q') l⁻¹`, so `x' p x` is `‖m x‖²`, the squared
length of one vector of as many values as there are individuals, and `m x`
costs a triangular solve of the block where the matrix costs a product. In
the same way `y' p y` is `‖m y‖²`, `num` is `(m y)' (m x)`, and what the
variant leaves of the trait is `‖m y - beta m x‖²`, the trait through the
factor less the variant through it times the variant's effect, formed from
the residuals of the one variant rather than subtracted.

The **threshold of Open 2** is one rule in four places: a quantity that has
fallen to the tested individuals times 2.2e-16 of the scale it was formed
from is the rounding of a cancellation, and the variant gets the three NaNs
of "The variants that have no answer". The fourth of the four places is the
Wald test's `y' p y` less `num² / den`, against `y' p y`.

The **GRAMMAR-Gamma approximation**, which `use_grammar_gamma_approx` turns
on and which is off by default, replaces `den` with one factor times the
squared length of the variant's centered dosages and makes no product at
all. It is untouched by everything here.

Every time and every column below is on the owner's Apple M5 Pro, 18 cores,
macOS 27, native `aarch64`. **Accelerate** is the default build's linear
algebra and **faer** is what `--no-default-features` links and what a
browser runs.

## What is decided, and the recommendation

Three things are open, and the third is the one Open 2 states.

1. Whether the two mixed models keep the projection as a factor and solve a
   variant against it, or keep it as the matrix and multiply. That is H3.
2. Whether the Wald test forms `‖m y - beta m x‖²` or subtracts
   `num² / den` from `y' p y`.
3. Whether the threshold of Open 2 stays at the tested individuals times
   2.2e-16 of `y' p y`, or moves to the square of that share, which is what
   a quantity that is now a squared length rather than a length would take.

**Take the factor, form the residual, and leave the threshold where it is.**

The three mixed model passes over 100000 variants of 1000 individuals go
from 0.561 s, 0.571 s and 0.753 s to 0.484 s, 0.542 s and 0.662 s, and every
column of both reference panels moves by at most 4.256e-13 of the largest
value of its own column, where the two arithmetic backends already sit
1.198e-08 apart. No variant of either panel gains or loses its answer.
Forming the residual costs nothing that the clock can see and takes a
finite `beta` beside a NaN `se` out of reach whatever the threshold does.
Moving the threshold is the one of the three that buys nothing: the band it
would open answers the fixture of six individuals with an `se` of 9.742e-16
on Accelerate against 1.979e-15 on faer and a p-value of 2.039e-45 against
1.708e-44, a factor of 2.0 and a factor of 8.4 between two builds of the
same library.

The first two are one change and not two. With the factor in place and the
subtraction kept, that same fixture is answered on both backends, with an
`se` of 2.980e-8 and a p-value of 5.837e-23 on Accelerate and 4.790e-8 and
2.424e-22 on faer: `y' p y` and `den` are exact squared lengths there, so
the rounding of the cancellation between them no longer falls below 0 and
the sign no longer refuses the variant. Forming `‖m y - beta m x‖²` puts what
the variant leaves at 2.847e-30 of `y' p y` on Accelerate and 1.174e-29 on
faer, against a threshold of 1.332e-15 of it, and refuses the variant again
on both.

## The algebra, checked

The identity is `p = m' m`, and from it `x' p x = ‖m x‖²`,
`x' p y = (m x)' (m y)` and `y' p y - num² / den = ‖m y - beta m x‖²` with
`beta = num / den`. The last of the four is the one Open 2 turns on. Writing
`a` for `m y` and `b` for `m x` over these three sentences alone, it holds
because `‖a - beta b‖²` is `a'a - 2 beta a'b + beta² b'b` and `beta` is
`a'b / b'b`, so the two middle terms leave `a'a - (a'b)² / (b'b)`, which is
`y' p y` less `num² / den`.

Two cargo tests check it in popnei's own arithmetic rather than on paper.

**Against the matrix, on real data.**
`the_denominator_from_the_factor_is_the_one_from_the_projection_matrix` of
`logistic_mixed` fits the logistic mixed model's linearization of each
reference panel at GMMAT's variance, builds the projection matrix the way
the fit built it before this change, and compares the `x' p x` of all 1200
variants of the panel through the matrix with the squared length of the
variant through the factor. The worst of the 1200 is 2.077e-15 of itself on
Accelerate, on the panel with every genotype called, and 2.016e-15 on faer,
on the panel with genotypes missing. `the_projection_of` of that module is
kept under `cfg(test)` for this: it is the code that built the matrix, and
the two routes are not each other's arithmetic.

**Entry by entry, and for all four quantities.** The four tests of
`projection` build a covariance of 12 individuals whose entries fall off
with the distance between them, a design of an intercept and one covariate,
and six variants of which the last is a combination of the two columns of
the design. Every entry of `m' m` is compared with the entry of the matrix,
144 of them, because a factor that is wrong in one direction alone would
survive any product with one vector; and `x' p x`, `x' p y`, `y' p y` and
what the Wald test's subtraction gives are compared for the five variants
the design does not explain. The worst over all of them is 1.627e-15 of the
largest entry on Accelerate and 9.037e-16 on faer.

The sixth variant, the one the design explains, is where the two routes
agree about nothing, and that is what forming the quantity buys: against a
squared length of 15.25, the matrix gives an `x' p x` of -4.753e-16 on
Accelerate and -1.136e-15 on faer, the rounding of a cancellation and below
0 on both, while the factor gives 3.722e-30 and 2.169e-30, the square of
that rounding, which no rounding can take below 0.

## What it costs

The whole study of 100000 variants x 1000 individuals from
`/Users/jose/devel/popnei-bench/bigcalled.vars`, `cargo bench --bench gwas`
at 18 cores, the best of three runs after one that is not timed, and the
best of three such measurements, since another session was working on the
machine:

| pass | before | after |
|---|---|---|
| `lmm`, score test | 0.561 s | **0.484 s** |
| `lmm`, Wald test | 0.571 s | **0.542 s** |
| `glmm`, score test | 0.753 s | **0.662 s** |
| `lmm`, score, approximation on | 0.197 s | 0.200 s |
| `glmm`, approximation on | 0.391 s | 0.389 s |
| `lm` | 0.124 s | 0.124 s |

The two approximating studies and the two models with no kinship make no
solve and no product with the projection, and the 0.003 s the first of them
gains is the null model, below.

**The block testing alone**, from the phase clock of the cargo feature
`bench-phases`, which no build popnei ships turns on: 0.485 s to 0.403 s for
the linear mixed model's score test and 0.543 s to 0.462 s for its Wald
test. The Wald test costs 0.058 s more than the score test before the change
and 0.059 s after, so **forming the residual costs nothing the clock can
separate from the run to run spread**, which is 0.010 s on these runs. What
the Wald test pays over the score test is the Student t tail, not the
residual.

**The solve against the product, at the shape a block has.** A triangular
solve of 5000 right hand sides against a 1000 x 1000 lower triangular factor
takes 14.6 ms, against 23.1 ms for the product of a 5000 x 1000 block with a
1000 x 1000 matrix, best of four. That is 0.63 of the product where the
arithmetic is half of it, so the solve runs at about 0.79 of the product's
rate; H3 set its gate at 0.70 and this passes it. Twenty blocks make a pass,
so the arithmetic of a pass falls by 0.17 s where the block testing falls by
0.08 s. The 0.09 s between them goes on the copy of the block that the solve
overwrites, the small product that gives each variant against each of the
columns of `q`, and the pass that takes the design out of each row.

**The null model costs 0.004 s more** at 1000 individuals, 0.0496 s to
0.0536 s, measured over five runs of a study of 5000 variants where the pass
is a twentieth of the panel's and the spread is 0.003 s. It factors the
covariance and inverts the triangular factor where it made one product of
two individuals by individuals matrices. Both grow with the cube of the
individuals, so at the 10000 individuals of `docs/objectives.md` this is
about 4 s, and nothing here measured it there.

## What it does to the columns

Every column of both mixed models over the 1200 variants of both reference
panels, from the same build through the Python package: `beta`, `se` and
`p_value` of the Wald test and the score test of the `lmm`, of the score
test of the `glmm`, of both under the GRAMMAR-Gamma approximation, and the
two variances, the heritability and the covariate effects of each null
model. Thirty-seven columns.

| compared | worst move |
|---|---|
| before against after, Accelerate | 4.256e-13 of the largest value of the column |
| before against after, faer | 1.864e-13 of it |
| Accelerate against faer, before | 1.198e-08 of it |
| Accelerate against faer, after | 1.198e-08 of it |

So the change moves a column 28000 times less than the two backends already
differ, and it leaves that difference where it was. The two null models come
back bit for bit what they were: the restricted maximum likelihood search
and the logistic fit are untouched. No value that was a number becomes a NaN
and none the other way; on both panels no variant of either model has ever
been refused.

Against the reference programs, all nine compared columns are the same to
five digits before and after, at the same worst variant. rrBLUP 4.6.3's
`GWAS` with `P3D = TRUE` and GMMAT 1.5.0's `glmm.score` are 1e-5 or so away
from popnei and the change moves popnei by 1e-13, so nothing of it is
visible there. The worst of each, as a share of what "How it is verified" of
the spec allows:

| column | worst | allowed | share |
|---|---|---|---|
| rrBLUP, `-log10(p)` of the Wald test, panel with every genotype called | 1.9973e-5 | 1e-4 | 20.0% |
| GMMAT `lmm`, `1 / se²`, that panel | 4.4268e-6 | 1e-5 | 44.3% |
| GMMAT `lmm`, `log10(p)`, that panel | 4.6204e-5 | 1e-4 | 46.2% |
| GMMAT `lmm`, `1 / se²`, panel with genotypes missing | 5.4234e-6 | 1e-5 | 54.2% |
| GMMAT `lmm`, `log10(p)`, that panel | 4.7552e-5 | 1e-4 | 47.6% |
| GMMAT `glmm`, `1 / se²`, panel with every genotype called | 5.3418e-6 | 1e-5 | 53.4% |
| GMMAT `glmm`, `log10(p)`, that panel | 8.4965e-6 | 1e-4 | 8.5% |
| GMMAT `glmm`, `1 / se²`, panel with genotypes missing | 4.9397e-6 | 1e-5 | 49.4% |
| GMMAT `glmm`, `log10(p)`, that panel | 8.6152e-6 | 1e-4 | 8.6% |

Every number of that table is the same on Accelerate and on faer to the
digits printed, before the change and after it.

## The threshold does not move

The fixture is the one Open 2 already carries: six individuals, one
covariate, an identity kinship and a trait built as `2 + 3*cov + 1*dosage`,
so that the variant explains the whole of what the null model left and
`num² / den` equals `y' p y` in exact arithmetic. What the first variant of
it gets, measured on 25 September 2026 on both backends:

| what the code does | Accelerate | faer |
|---|---|---|
| the matrix and the subtraction, which is `main` | three NaNs | three NaNs |
| the factor and the subtraction | `se` 2.980e-8, `p` 5.837e-23 | `se` 4.790e-8, `p` 2.424e-22 |
| the factor and the formed residual | three NaNs | three NaNs |
| the same, with the share squared | `se` 9.742e-16, `p` 2.039e-45 | `se` 1.979e-15, `p` 1.708e-44 |

The second row is why the two halves are one change: the factor alone takes
the guard away, because the sign of the cancellation was doing half the work
and there is no cancellation left to have a sign.

The fourth row is the argument for leaving the threshold alone. What is
left, formed, is right to about the tested individuals times 2.2e-16 of the
length of `u`, so the square of that share is where the quantity stops being
its own rounding, and a threshold there is the honest translation of the
rule onto a squared length. It is also a threshold whose band the two builds
of popnei answer a factor of 2.0 apart in `se` and a factor of 8.4 apart in
the p-value, because the p-value of a t with `df` degrees of freedom moves
by the `df`th power of a move in `se`. A variant in that band explains all
but about 4.4e-14 of what the null model left. The variant of the panel with
every genotype called that comes nearest leaves 0.919 of `y' p y`, thirteen
orders above the threshold, measured over its 1200 on 25 September 2026. So
keeping the threshold costs nothing that has been seen, and moving it buys a
band where the answer is decided by which library the build links.

## A kinship with a negative eigenvalue

A kinship whose smallest eigenvalue is below 0 is not a covariance and has
no real factor, and this change gives it exactly what `main` gives it. The
fixture is the kinship `plink2 --make-rel square bin` wrote for the panel
with every genotype called, whose largest eigenvalue is 17.269 and whose
smallest is -3.324e-15, with its smallest eigenvalue forced to a chosen
value and the matrix built again from its own eigenvectors.

| smallest eigenvalue | of the largest | `lmm` | `glmm` |
|---|---|---|---|
| -5.0 | -29.0% | answered, all 1200 | refused, at row 112 |
| -2.0 | -11.6% | answered, all 1200 | refused, at row 188 |
| -0.5 | -2.9% | answered, all 1200 | answered, all 1200 |

Every number of both models is identical before the change and after it,
`genetic_variance` 1.221615453561016 and `heritability` 0.7810960460274852
at -5.0, and the two refusals name the same row.

The two models differ here and did before. The linear mixed model
eigendecomposes the kinship and clamps every eigenvalue at 0 before the
covariance is built, so the covariance is the nearest matrix that is a
covariance and always has a real factor; what is factored is the matrix the
model was already inverting. How negative an eigenvalue may be before the
kinship is refused rather than clamped is **Open 4** of the spec, which this
does not touch. The logistic mixed model builds the covariance of its
working trait from the kinship the user gave and factors it once a round, so
a kinship far enough below 0 makes that factorization fail and the fit
refuses it with `GwasKinshipNotACovariance` and the row it stopped at; the
factored projection reuses the factor that fit already made, so the refusal
is the fit's and is where it was.

## What the decision changes in the spec

Nothing of `docs/specs/gwas.md` is changed on the branch: the numbers it
states all hold, and closing Open 2 is the owner's. If the recommendation is
taken, these are the sentences that owe a change, in one commit before the
code:

- **Open 2**, the paragraph beginning "Or form the residual exactly", which
  prices the exact form at "a matrix product for every variant, roughly
  doubling that test", and the recommendation that rests on it. What
  replaces it is that the exact form is a pass over one vector per variant
  and that the threshold stays where it is for the reason above.
- **Open 2's table**, whose third row says the scale of `x' p x` and whose
  fourth says what the Wald test subtracts. Neither number moves; what
  changes is that both quantities are now squared lengths and cannot fall
  below 0, so the sign half of each guard is structural.
- **"The linear mixed model"** and **"The logistic mixed model"**, where
  each says that the two products are the whole cost of a block. One of the
  two is now a triangular solve.
- **"The GRAMMAR-Gamma approximation"**, where what a variant that falls
  back to the exact denominator costs is "one product of one variant with
  the projection matrix". It is one solve of one variant against the factor.
- **"Speed"**, which has no popnei number at all and now has six.
- **"The solve against a triangular matrix"** of `docs/specs/linalg.md`,
  which names the callers of the lower half: the fit of the logistic mixed
  model's null is no longer the only one, and the new one asks for one right
  hand side per variant of a block rather than one per individual.

## What is not known

- **Nothing was measured at 10000 individuals**, the largest dataset of
  `docs/objectives.md`. The per variant half improves with the square of the
  individuals and the null model's extra work grows with the cube, so the
  two cross somewhere above 1000 and this says nothing about where.
- **Nothing was run under WebAssembly.** The faer numbers here are a native
  build linking faer, which is the same arithmetic a browser runs but not
  the same build. `cargo wasm-check` and `cargo wasm-check-js` pass.
- **0.09 s of the arithmetic H3 saves is still being paid**, on the copy of
  the block that the solve overwrites, the product that gives each variant
  against each of the columns of `q`, and the pass that takes the design out
  of each row. The squared lengths and the residual are then taken in two
  more passes over the block, one variant after another on one thread.
  Fusing all of it into the pass that takes the design out, so that a block
  is read once instead of three times and on the threads of rayon rather
  than on one, is what would recover it; it was not built and not measured.
- **The GRAMMAR-Gamma approximation is untouched**, and it still halves both
  mixed models: 0.200 s and 0.389 s against 0.484 s and 0.662 s.

## What is asked

Whether to merge `gwas/exact-residual` into `main`, one commit, every check
of the `coding` skill clean at its head: `cargo fmt`, `cargo clippy
--workspace --all-targets -- -D warnings`, 1169 cargo tests of the workspace
on Accelerate and 1019 of the core crate and 136 of the linear algebra crate
on faer, `cargo wasm-check`, `cargo wasm-check-js`, `ruff format`, `ruff
check` and 556 pytest tests.

And whether Open 2 closes with the recommendation above, which is the six
sentences of the section before last.
