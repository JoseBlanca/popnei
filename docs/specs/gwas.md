# The gwas module

23 September 2026. The association study tells a user of popnei which of
their variants are associated with a trait they measured: for each variant,
the effect of one more copy of a non major allele, how uncertain that effect
is, and the p-value of the test that the effect is 0. There is no code. This
spec covers the whole `gwas` row of section 9 of `docs/architecture.md`, the
four null models, the two tests and the two distributions that turn a
statistic into a p-value.

It stands on four specs. `docs/specs/kinship.md` gives the matrix that the
two mixed models take, and the dosage rules that this module uses are the
ones it takes from `docs/specs/pca.md`. `docs/specs/linalg.md` gives the
seven operations the fits are built from. `docs/specs/variant.md` gives the
allele counts of one variant and the pass stats.

How the null model of the logistic mixed model is fitted is not pyNei's way,
and `docs/reports/glmm-method/README.md` has the measurements that led there.

## What every model shares

### What it gives

A **trait** is one number per individual: **continuous**, a measurement, or
**binomial**, 0 or 1. The user may also give **covariates**, other numbers
per individual whose effect on the trait is not of interest but has to be
taken out, such as the sex or the field the plant grew in.

Two individuals that share ancestry resemble each other at every variant and
in the trait, so a variant that only marks the ancestry looks associated.
popnei accounts for this in two ways, and a user can take either or both.
The top principal components of the panel, from `Kinship.principal_components`
or `do_pca_from_variants`, go in as covariates, which is enough for
individuals that are not close relatives. Or the **kinship**, the matrix of
`docs/specs/kinship.md`, goes in as the covariance of a random effect, which
is what a panel with families in it needs; a model with such a term is a
**mixed model**. Both at once is the Q+K model of a strongly subdivided
panel.

Trait against kinship gives the four models, by the names the literature
gives them:

| | no kinship | with kinship |
|---|---|---|
| **continuous** | `lm`, a linear model | `lmm`, a linear mixed model |
| **binomial** | `glm`, a logistic regression | `glmm`, a logistic mixed model |

Each is fitted once without any variant in it, which is the **null
model**, and then every variant is tested against what that model left
unexplained. Fitting the null once and reusing it for every variant is what
makes a study of a million variants possible; holding the variance
components of a mixed model at their null values for every variant is called
P3D, and rrBLUP, GMMAT and EMMAX, the program that named it, all do it.

There are two tests, and which one is available depends on the model.

- A **Wald test** fits the model again with the variant in it and asks how
  many of its own standard errors the variant's effect is away from 0. It
  needs a fit per variant, so it is used where that fit is cheap.
- A **score test** never fits the model with the variant in it. It asks how
  steeply the fit would improve if the variant's effect were let off 0,
  measured at the null model, and compares that slope with how uncertain it
  is. It costs no fit per variant.

Under the null both statistics have the same distribution in large samples;
they differ in what they cost and in what they assume. The default is the
Wald test where a per variant fit is cheap, a continuous trait or a binomial
one without a kinship, and the score test for a binomial trait with a
kinship, where a Wald test would mean one mixed model fit per variant.

Each variant becomes one number per individual, its **dosage**, and the
dosage, the major allele and the rule that gives a genotype with any allele
missing the mean of its variant are those of `docs/specs/pca.md`. They are
computed over the individuals that are tested and not over the whole panel,
so a variant's mean, its frequency and whether it varies at all are all of
those individuals. `allele_freq` of the result is the mean dosage over the
ploidy: the frequency of the alleles that are not the major one, which for a
biallelic variant is the minor allele frequency, and which is what plink2
reports as `A1_FREQ`.

### Its Python function, and its TypeScript one

```python
calc_gwas(
    variants: Variants,
    phenotype: pandas.Series,
    trait: TraitType | str,
    covariates: pandas.DataFrame | None = None,
    kinship: Kinship | None = None,
    test: TestType | str | None = None,
    use_grammar_gamma_approx: bool = False,
    transform_to_biallelic: bool = False,
) -> GWASResult
```

`phenotype` is a series indexed by individual name; a frame of one column is
taken as that column. `trait` is `"continuous"` or `"binomial"`, the two
values of `TraitType`. `test` is `"wald"` or `"score"`, the two values of
`TestType`, and `None` takes the default above.

`GWASResult` is a frozen dataclass. `stats` is a frame with one row per
variant, in the order the variants came: `chrom`, `pos` and `id` when the
variants carry them, then `allele_freq`, `beta`, `se` and `p_value`. `beta`
is the effect of one more copy of a non major allele, in the units of the
trait for a continuous one and as a log odds ratio for a binomial one.
`null_model` is a `NullModel`, with `model`, one of the four values of
`GWASModel`; `covariate_effects`, a series over the intercept and the
covariates; `residual_variance`, `None` for a binomial trait;
`genetic_variance`, `None` without a kinship; `heritability`, the genetic
variance over the sum of the two, only for the `lmm`; and `num_individuals`.
Then `trait`, `test`, `individuals`, the names of those that were tested as
a tuple, `used_grammar_gamma_approx`, and `pass_stats`.

It mirrors `calc_gwas` of `pynei/gwas.py`. The differences from pyNei, which
`docs/objectives.md` asks to be written down:

- `samples` is `individuals`, in `GWASResult` and in `NullModel`, the name
  `docs/glossary.md` gives.
- `num_threads` is not an argument, as no calculation of popnei has one.
- `pass_stats` is new, as it is for every consumer of a `Variants`.
- `chrom`, `pos` and `id` are columns of `stats` whenever the source has
  them. pyNei leaves each out when the chunk has no such column, which its
  own test asserts; popnei asks the reader for them and gives them.
- `transform_to_biallelic` is new, and it is the argument
  `do_pca_from_variants` and `calc_kinship` have, for the same reason: a
  variant with more than two alleles among its called genotypes is refused
  unless it is true, where pyNei collapses every allele that is not the
  major one silently. The owner decided it for the kinship on 23 September
  2026, for being the more explicit of the two to a user, and it holds here
  for the same reason; neither reference panel has such a variant, so no
  literal moves. `docs/objectives.md` asks the calculations that collapse a
  multiallelic variant to say so.
- The null model of the logistic mixed model is fitted by another route,
  which gives the same numbers 1.88 to 2.06 times faster. It is the item for
  that model, below.

In TypeScript it is `calcGwas(variants, {phenotype, trait, covariates,
kinship, test, useGrammarGammaApprox, transformToBiallelic})`. `phenotype` is an object of
individual name to number, and `covariates` an object of covariate name to
such an object. The result has `stats` with each column as its own typed
array, `chrom` and `id` as arrays of strings, `pos` as a `Float64Array` and
`alleleFreq`, `beta`, `se` and `pValue` as `Float64Array`s, and then
`nullModel`, `trait`, `test`, `individuals`, `usedGrammarGammaApprox` and
`passStats`.

### Which individuals are tested, and the design

The individuals tested are those that have a phenotype: in the `phenotype`
series, not NaN, and present in the `Variants`. They are kept **in the order
the variants have them**, whatever order the phenotype was given in. An
individual in the phenotype that the `Variants` does not have is a
`ValueError` naming it, and a repeated individual in the phenotype is one
too.

The **design** is the matrix of one row per tested individual and one column
per number the model fits: a column of ones for the intercept, always added,
and one column for each covariate. The covariates are a frame indexed by
individual, which must cover every tested individual and hold no missing
value and no value that is not a number; each of the three raises a
`ValueError`, and the one for a value that is not a number says to code a
categorical covariate, with `pandas.get_dummies` for instance.

Two refusals protect the fits. A design whose columns are not independent, a
covariate that is constant or a copy of another, is a `ValueError` saying
the covariates are collinear; it is found with the rank of
`docs/specs/linalg.md`, whose tolerance is numpy's, so a design popnei
refuses is a design pyNei refuses. And a design with no more rows than
columns plus one is refused, since there would be nothing left to estimate
the uncertainty from.

For a binomial trait, a phenotype that is not 0 or 1 everywhere is a
`ValueError`, and so is one where every individual has the same value.

Asking for a test the model does not have is a `ValueError`: the score test
for a continuous trait with no kinship, since the only test of a linear
model is its t test; and the Wald test for a binomial trait with a kinship,
since it would fit one mixed model per variant.

The core crate is given the tested individuals as their positions among the
individuals the source has, with their phenotype and the design already
built, as "The Rust interface" below has it. So the refusals about a name
and about a frame are made where those are, in the Python and the TypeScript
layers: an individual of the phenotype that the `Variants` has not, a
covariate that does not cover a tested individual, and a covariate value
that is missing or is not a number. The core makes the rest over the
positions and the numbers it holds, and three more that only it can see.

The positions rise, and any other order is a `ValueError`. They are the
order the source has the individuals in, and the phenotype, the rows of the
design and the dosages of a block are read together row by row, so an order
that is not the source's measures one individual's trait against another
individual's genotypes. A position that repeats the one before it is the
repeated individual above, and one that falls back is refused as an order
that is not the source's, which is also what a repeat with another
individual between its two halves gives. A phenotype that is not a finite
number is a `ValueError` naming where it is: the individuals tested are
those that have a phenotype, so a NaN is an individual that should not have
been tested at all, and an infinity would carry through the null model into
the effect of every variant. And the phenotype holds one value for each
tested individual, the design one row of its columns for each, and the
design has the column of ones at least; none of the three can be reached
from Python or from TypeScript, which build the three from the same
individuals, so each is a `RuntimeError`.

### The variants that have no answer

A variant whose dosages are all the same among the tested individuals has no
variance and cannot be tested. Its row is still in `stats`, with its
`allele_freq`, and `beta`, `se` and `p_value` are NaN. This is where a
variant with one allele lands, and also one where every tested individual is
heterozygous, and one with no called genotype at all, whose dosages are all
the mean of nothing, which pyNei sets to 0.

A variant of the logistic Wald test whose fit runs away also gets three
NaNs, which the item for that model says.

A third case exists and the spec did not describe it. The denominator of
both score tests is `x' p x`, a quadratic form that is 0 or above in exact
arithmetic and that can round just below 0 for a variant with almost no
variance left after the covariates and the kinship are taken out. `beta` is
then `num` over a tiny negative number, a large value of whichever sign the
rounding chose, and the statistic `num² / den` is negative. What such a
variant gets is **Open 2**, below.

`test_monomorphic_and_missing_variants` of pyNei asserts exactly this on 50
variants of 60 individuals where the first has one allele and the second has
no called genotype: the first two p-values are NaN and the other 48 are
between 0 and 1.

### How it runs

One pass over the blocks. The null model is fitted before the pass, from the
trait, the design and the kinship alone, and no block is read for it.

Before its rows are read, every block is checked against the reader that
gave it: a block with no variants, and one whose individuals or ploidy
disagree with the reader's, are the two reader defects `docs/specs/block.md`
names, and this pass refuses both with the errors that spec gives them. The
check is not a formality. The drive over the rows pairs the genotypes cut
into one chunk per variant with the output buffer cut into one row per
individual, and a buffer sized from a ploidy that is not the block's comes
out with fewer rows than there are variants; the pairing then truncates to
the shorter of the two, so the variants past that point are not read at all
and the pass returns as though the block had held only the ones it managed.
Read on 23 September 2026 in `the_standardized_rows` of
`crates/popnei/src/variant.rs`: one variant of five individuals at a block
ploidy of 2, with a caller passing a ploidy of 5, gives 1 chunk of genotypes
against 0 rows of buffer, so no row runs, no error is raised, and the
variant is gone. A study that lost variants that way would report a count
the user could mistake for variants that had no variance. Then
each block is turned into its dosages, with rayon across the rows, and the
variants that vary are tested together as a matrix, because every test but
the logistic Wald one is a product of the block with something the null
model holds. What is kept from one block to the next is nothing but the
rows of `stats` already computed, which grow with the variants and not with
the individuals: four numbers per variant, 32 MB for a million.

`reblock` goes before it, so that a filter's uneven blocks do not reach the
matrix work, and because the `use_grammar_gamma_approx` pass below takes the
first variants of the first block and its answer would otherwise depend on
what the source gave.

With `use_grammar_gamma_approx` there is a second pass, which is opened
first and reads one block. Everything else is one pass.

### How it is verified

Four programs, all run on 23 September 2026 by
`tests/reference/gwas/make_reference.py`, which runs them on popnei's own
VCFs: plink2 v2.0.0-a.7.7, and R 4.6.1 with GMMAT 1.5.0 and rrBLUP 4.6.3.
pyNei ran the same four on its own vars files of the same genotypes when its
reference was made, and every number this script produced matches what pyNei
stored to the bit: the largest difference over the nine files and their 58
numeric columns, 1200 variants in each but the two null models, is 0. So
the numbers do not depend on which library read the genotypes.

The datasets are the two panels of `docs/specs/kinship.md`, 200 individuals
and 1200 biallelic diploid variants on two chromosomes, once with every
genotype called and once with 3 in 100 missing whole. The trait was
simulated from the genotypes with a heritability of 0.5 and five causal
variants of effect 0.6, with two covariates, `cov1` continuous and `cov2`
binary, and the three subpopulations differing in their mean so that the
structure confounds the trait. `tests/reference/gwas/phenotypes.csv` holds
the traits, the covariates and the subpopulation of each individual, and
`causal_vars.csv` the five causal variants, `var0052`, `var0629`, `var0751`,
`var1137` and `var1188`. The literals below also carry `var0000`, which is
not causal.

The mixed models are given the kinship that plink2 wrote for the panel with
every genotype called, so that they are tested against a kinship that came
from neither popnei nor pyNei.

**How many digits each reference gives, and what that costs.** Two of the
five files are printed to six significant digits, plink2's `--glm` and
GMMAT's `glmm.score`; neither program has a binary form for them, so that is
all there is. Six significant digits round a value by up to 5e-6 of itself,
so a comparison against one of those two can have at most twofold headroom
at a tolerance of 1e-5 relative: it says that popnei computes the same
quantity, and it would not catch an arithmetic error smaller than the
printing. That rounding is relative, so an absolute tolerance
would hold on this panel, whose values are small, and break on data whose
values are larger, for an implementation that is right. A tolerance relative
to each value breaks the other way, and the next paragraph is that.

The other three are full precision: rrBLUP's `GWAS`, R's `anova(glm, test =
"Rao")` and GMMAT's `glmmkin` null models, which the reference script writes
itself. Their tolerances are the distance between two fits and not the width
of a printed digit.

**The check with headroom is pyNei**, at 1e-9 relative over every column of
every model, and the worked example at 1e-12, both against float64 with
nothing rounded away. Those two are what would catch a wrong digit that the
printed references could not. No work package rests on a printed reference
alone.

**A tolerance is against the scale of what is estimated, not against each
value.** An association study is mostly null: for the great majority of
variants the effect is 0, and what comes out is whatever the rounding of a
sum of cancelling products left, a number whose own magnitude means nothing.
Asking two implementations to agree to a share of *that* asks for accuracy
no arithmetic has, at exactly the variants where the null is true, which is
most of the genome. So `beta` and `se` are compared within a tolerance times
`se`, the scale of what the study is measuring, and never within a tolerance
times `beta`; `p_value` is compared in `log10`, which is already a scale;
and a check over a vector of numbers is against the largest of them and not
each one.

Where the right shape comes from, so that the next quantity does not have to
be got wrong first. A bound is on the rounding of the sum that produced the
number, and the rounding of a sum of `m` products is about `m` times the
distance from 1 to the next `f64`, 2.2e-16, times the largest term of the
sum. So the bound goes against whatever bounds the terms, and the terms are
what the quantity is built from and not the quantity itself, which is why a
value that cancelled to near 0 is no guide to its own error. For the kinship
that scale is the largest entry of the matrix, since the sum of the absolute
products of a pair is at most `m` times the square root of the two diagonal
entries and so at most `m` times the largest entry; the bound is then loose
by the ratio of the largest entry to that square root, measured at 1.4 times
on both panels. For an effect size it is `se`, which is what the study's own
arithmetic says the effect is uncertain by.

The kinship met this on 23 September 2026 and it is why the rule is here.
Its first bound was 1e-12 relative to each entry of the matrix, which looked
sound: the worst entry as a ratio was 3.31e-13. It broke on the second
backend, and not at a large entry. faer missed at an entry of 1.29e-05 whose
difference from plink2 was 1.9e-17, a smaller difference than the ones at
entries a hundred times larger, every one of which passed. Its rule now is
that each entry is within 1e-13 of the **largest** entry of the matrix, at
which the worst of four measurements is 3.6e-16 and 4.5e-16 on Accelerate
and 3.3e-15 and 2.3e-15 on faer.

**A tolerance is chosen against both backends of `docs/specs/linalg.md` and
not one.** faer sits about seven times further from plink2 than Accelerate
does on the same data, which is well inside what the order and the blocking
of the sums allow and is not a defect. The faer build is what runs in a
browser, so a bound fixed on Accelerate alone is a bound the wasm package
fails. `cargo test -p popnei --no-default-features` is the run that says so,
and it is in the `coding` skill since 23 September 2026.

And the numbers are where to start and not where to stop: 1e-9 and 1e-12
were chosen here for being small, which is not evidence of anything. Each is
lowered until it fails and set two or three times above where it broke, on
both backends, and the implementation plan records both numbers.

Each model item says what it is checked against and how closely. Three
checks are common to all four:

- Every column of `stats` against the reference program of that model over
  all 1200 variants, at the Python `calc_gwas`.
- The six variants above as literals in the cargo tests, at `calc_gwas` of
  "The Rust interface". What each model's literals are and what they are
  compared against differs, because the reference programs report different
  quantities, and each model item says which: for the `lm` and the `glm`'s
  Wald test, `beta`, `se` and `p_value` directly; for the `lmm`'s Wald test,
  `-log10(p_value)`, which is all rrBLUP reports; and for both score tests
  against GMMAT, `1 / se²`, which is the variance of the score, and
  `p_value`.
- Against pyNei, both libraries on the same panel with the same arguments,
  at the Python `calc_gwas`: `beta`, `se` and `p_value` within 1e-9
  relative, and the variants that have NaN exactly the same ones.

Two more hold for every model. That the block size changes nothing: the same
panel read in blocks of 77 gives `stats` equal to the default within 1e-12
relative. It is `test_chunks_and_threads_do_not_matter` of pyNei, which
compares with `pandas.testing.assert_frame_equal` and so at its default of
1e-5 relative; 1e-12 is popnei's own, because every test of a variant is
independent of the others and only the order the blocks were summed in can
move a digit. And
that the study finds what was planted: of the 10 variants with the smallest
p-value under the `lmm`, at least 3 are among the 5 causal ones.

In TypeScript, `calcGwas` is tested under node against the same six literals
for each model.

## The linear model

### What it gives

A continuous trait, covariates and no kinship. The trait is a straight line
in the covariates plus the variant, and the test is the t test of that line:
`beta` is the slope on the variant, `se` its standard error, and `p_value`
the two sided probability that a t with `n - c - 1` degrees of freedom, `n`
the individuals and `c` the columns of the design, is further from 0 than
`beta / se`. It is what plink2's `--glm` computes.

The null model is fitted with a thin QR of the design, `d = q r`: the
coefficients are the `c` of `r c = q' y`, the residuals are `y - q q' y`,
and the residual sum of squares is their squared length. `residual_variance`
is that sum over `n - c`.

Every variant is then tested on the residuals. The covariates are taken out
of the dosages of the block in one product, `x - (x q) q'`, and after that
the effect of the variant is the plain slope of the trait's residuals on the
variant's residuals:

    xx   = the squared length of each variant's residuals
    num  = each variant's residuals times the trait's residuals
    beta = num / xx
    rss  = the null's residual sum of squares - beta * num
    se   = sqrt(rss / (n - c - 1) / xx)

That `rss` is what the variant leaves unexplained, so each variant gets its
own estimate of the residual variance, which is what makes this a t test and
not a normal one.

### How it is verified

Against plink2 `--glm hide-covar` on the panel with every genotype called,
with `cov1` and `cov2` as covariates, which writes
`tests/reference/gwas/plink2.panel_called.glm.linear.tsv`, 1200 rows with no
`NA`. plink2 tests the minor allele and popnei the non major one. **On this
panel they are the same**, so `allele_freq` is plink2's `A1_FREQ` and the
signs agree; every genotype of it is called, and that is what makes the two
conventions coincide.

They are not the same in general, and `allele_freq` can pass a half. The
major allele is the most frequent among the called **alleles**, which counts
the called half of a half called genotype, while the mean that becomes
`allele_freq` is over the whole called **genotypes**, which a half called
one is not. So the allele the dosages are counted from is not always the one
whose frequency is below a half. Run on 23 September 2026 on one variant of
five individuals, `0/. 0/. 0/. 0/. 1/1`: the major allele is 0, on four
called halves against two, while the only whole genotype is `1/1`, so the
mean dosage is 2 and `allele_freq` is 1.0. popnei and pyNei agree on this,
so it is a divergence from plink2 and not from the oracle, and neither
reference panel shows it: one has every genotype called and the other has
them missing whole.

Over all 1200 variants: `allele_freq` within 1e-6 absolute, since it is a
frequency and lies between 0 and 1; `beta` and `se` within 1e-5 times the
`se` of that variant, for the reason above, which on this panel is between
1.2e-6 and 1.6e-6 absolute; and `p_value` within 1e-5 relative.

Six significant digits round `beta` and `se` by up to 5e-7 absolute here, so
the printing takes up to 41 per cent of that tolerance and leaves the
arithmetic the rest. A tolerance is a budget shared between the rounding of
the number it is compared against and the difference it is meant to catch,
and the first share is worth computing rather than assumed to be small.

The six literals are held to the same tolerance as the whole columns, 1e-5
relative on all three. From plink2 on 23 September 2026:

| variant | beta | se | p |
|---|---|---|---|
| var0000 | -0.424136 | 0.139354 | 0.00265846 |
| var0052 | -0.697724 | 0.122348 | 4.28981e-08 |
| var0629 | -0.813852 | 0.161809 | 1.10646e-06 |
| var0751 | -0.0963977 | 0.130636 | 0.461451 |
| var1137 | -0.171137 | 0.145987 | 0.242511 |
| var1188 | -0.655393 | 0.138912 | 4.51958e-06 |

### The worked example

Six diploid individuals, `i0` to `i5`, one covariate beside the intercept,
and three variants. The trait is 2, 3, 5, 4, 4, 7 and the covariate 0, 1, 0,
1, 0, 1.

| variant | genotypes | dosages | kept |
|---|---|---|---|
| v0 | 0/0 0/1 1/1 0/0 0/1 1/1 | 0 1 2 0 1 2 | yes |
| v1 | 0/0 0/1 1/1 ./. 0/1 0/0 | 0 1 2 **0.8** 1 0 | yes |
| v2 | 0/1 0/1 0/1 0/1 0/1 0/1 | 1 1 1 1 1 1 | no, no variance |

The dosage in bold is the missing genotype of `i3` taking the mean of its
variant, `(0 + 1 + 2 + 1 + 0) / 5`. The null model has 6 individuals and 2
coefficients, so its residual sum of squares has 4 degrees of freedom and
each variant's test has 3. Run through pyNei at commit ef0ca6e:

| | intercept | covariate | rss | residual_variance |
|---|---|---|---|---|
| the null | 3.6666666666666683 | 1.0 | 13.333333333333336 | 3.333333333333334 |

| variant | allele_freq | beta | se | p_value |
|---|---|---|---|---|
| v0 | 0.5 | 1.5 | 0.600925212577332 | 0.088004892382756 |
| v1 | 0.4 | 0.3125 | 1.305204592306424 | 0.826200867452417 |
| v2 | 0.5 | NaN | NaN | NaN |

`v2` keeps its `allele_freq` of 0.5 and has no test, which is the rule of
"The variants that have no answer". This is the first cargo test of the
module, asserted at `calc_gwas` of "The Rust interface" within 1e-12
relative, and it is the one test that needs no reference program and no
reference file.

## The linear mixed model

### What it gives

A continuous trait with a kinship. Beside the covariates the trait carries a
random effect whose covariance is the kinship times a variance, so that two
related individuals are expected to resemble each other before any variant
is looked at. The covariance of the trait under the null is

    V = genetic_variance * k + residual_variance * i

with `k` the kinship and `i` the identity.

The two variances are estimated by **restricted maximum likelihood**, which
is maximum likelihood on the part of the trait that the covariates cannot
explain, so that fitting the covariates does not drag the variances down.
Only their ratio matters to the search: with `delta` the residual variance
over the genetic one, the kinship is eigendecomposed once, `k = e diag(l)
e'`, the trait and the design are turned by `e'`, and then every value of
`delta` costs one number per individual instead of a matrix. pyNei searches `log(delta)` over 101 points evenly spaced from -10 to 10,
takes the point with the smallest value and brackets it with its two
neighbours, and then runs 60 steps of a **golden section search**, which
shrinks a bracket that holds a minimum by a constant ratio each step and
costs one evaluation of the function per step. popnei reproduces that search
step for step, because the numbers it gives are GMMAT's, and reproducing it
needs all of: 101 points, the two neighbours of the best clamped at the ends
of the grid, the ratio `(sqrt(5) - 1) / 2`, the two interior points taken as
`high - ratio * (high - low)` and `low + ratio * (high - low)`, the bracket
moved to whichever of the two has the smaller value, 60 steps whatever
happens, and `exp((low + high) / 2)` at the end. It is `_reml_delta` of
`pynei/gwas.py`.

The function it minimizes, with `l` the eigenvalues of the kinship, `u` the
design turned by the eigenvectors, `uy` the trait turned by them, and
`w = 1 / (l + delta)` one weight per individual:

    dvd   = u' diag(w) u
    coefs = the solution of `dvd coefs = u' (w * uy)`
    resid = uy - u coefs
    quad  = w' (resid * resid)
    value = sum(log(l + delta)) + (n - c) * log(quad) + log(det(dvd))

Each evaluation costs one number per individual and one matrix of the size
of the coefficients, which is what the eigendecomposition bought. With the
`delta` that minimizes it, and `w`, `coefs` and `resid` recomputed there:

    genetic_variance  = quad / (n - c)
    residual_variance = delta * genetic_variance

`log(det(dvd))` is the log determinant that comes off the same Cholesky
factorization the solve above uses, which is why `docs/specs/linalg.md` has
it as one of the seven.

The eigenvalues of the kinship are clamped at 0 before use. A kinship of
genotypes with nothing missing has none below 0 but for rounding, -4.8e-15
on the panel; the per pair denominators of `docs/specs/kinship.md` put them
there, -0.0321 on the panel with 3 in 100 genotypes missing, and a negative
eigenvalue would make `V` not a covariance.

`heritability` is the genetic variance over the sum of the two, which is the
share of the trait's variance that the kinship explains.

Every variant is then tested through the **projection matrix**

    p = V⁻¹ - V⁻¹ d (d' V⁻¹ d)⁻¹ d' V⁻¹

individuals by individuals, which takes the covariates out of anything it is
applied to and weights it by the covariance. With `x` the dosages of a
variant, `num` is `x' p y` and `den` is `x' p x`, and then:

- The **Wald test**, the default, is `beta = num / den` with
  `se = sqrt((y' p y - num²/den) / ((n - c - 1) * den))` and a two sided t
  with `n - c - 1` degrees of freedom. It holds the ratio of the two
  variances at the null and estimates their scale again with the variant in,
  which is what rrBLUP does.
- The **score test** is `beta = num / den`, `se = 1 / sqrt(den)` and a chi
  square with one degree of freedom of `num² / den`. It holds both variances
  at the null, which is what GMMAT does.

`y' p y` is the generalized residual sum of squares of the null over the
genetic variance, and the restricted maximum likelihood makes it exactly
`n - c`. A cargo test asserts that on the panel, 197 within 1e-6, which is
`test_reml_identity` of pyNei and is the cheapest evidence there is that the
fit reached its optimum. `y' p y` is in no result, so unlike every other
check of this spec this one is made at the private function that fits this
null, and it pins that function's signature. That is accepted: the
alternative is a field of `NullModel` that exists for one of the four models
and that no user would read.

### How it is verified

The Wald test against rrBLUP 4.6.3's `GWAS` with `P3D = TRUE`, which holds
the variance components at the null as popnei does. rrBLUP takes every fixed
effect as a factor, so only the binary covariate `cov2` was given to it, and
popnei is run with the same one covariate for this comparison. It reports
`-log10(p)` at full precision, so that is what is compared, over all 1200
variants within 1e-4 absolute, from
`tests/reference/gwas/rrblup.panel_called.lmm.tsv`. The tolerance is the
distance between two fits, not a printed digit: 1e-4 in `-log10(p)` is 2.3e-4
of the p-value itself.

The six literals, held to 1e-4 absolute as the whole column is, are
`-log10(p_value)` and nothing else, because `-log10(p)` is all rrBLUP
reports: var0000 0.215210, var0052 3.618699, var0629 4.339896, var0751
2.065325, var1137 1.319283, var1188 2.367279.

The score test against GMMAT 1.5.0's `glmm.score`, with both covariates, on
both panels, from `gmmat.panel_called.lmm.score.tsv` and
`gmmat.panel.lmm.score.tsv`. GMMAT reports the variance of the score, which
is `den`, and the p-value. Over all 1200 variants: `1 / se²` against GMMAT's
`VAR` within 1e-5 relative, and `|log10(p / p_GMMAT)|` below 1e-4. The
p-values are compared in `log10` because they span 23 orders of magnitude
and what a user reads is the exponent.

The six literals are `1 / se²` against `VAR`, within 1e-5 relative, and
`p_value` within 1e-4 in `log10`, the same as the whole columns. Both are
against six printed significant digits, so the first has twofold headroom
and the second has more, `log10` shrinking a relative difference. The
variance of the score and the p-value, with every genotype called and then
with 3 in 100 missing:

| variant | VAR | p | VAR, missing | p, missing |
|---|---|---|---|---|
| var0000 | 29.8774 | 0.360526 | 29.9961 | 0.495719 |
| var0052 | 43.8076 | 0.00118985 | 46.0763 | 0.00105895 |
| var0629 | 31.7241 | 4.81005e-05 | 31.9307 | 7.37888e-05 |
| var0751 | 44.5825 | 0.00439226 | 44.422 | 0.00675644 |
| var1137 | 43.3724 | 0.013926 | 42.1841 | 0.0125321 |
| var1188 | 47.3766 | 0.0010734 | 47.9106 | 0.00150577 |

With genotypes missing, GMMAT gives a missing genotype the mean of its
variant, which it calls `impute2mean` and which is popnei's rule too; that
is why the two agree on the second panel.

The null model against GMMAT's `glmmkin`, from `gmmat.null_models.tsv`,
which the reference script writes at full precision, within 1e-5 absolute,
which is how far two restricted maximum likelihood searches land apart: `genetic_variance` 1.221617, `residual_variance`
0.342359, and the three covariate effects 4.678021, 0.473361 and 1.110279.
`heritability` is 1.221617 / (1.221617 + 0.342359).

## The logistic model

### What it gives

A binomial trait and no kinship. The chance that an individual is a 1 is a
logistic curve in the covariates and the variant, and `beta` is a log odds
ratio: the change in the log odds of being a 1 for one more copy of a non
major allele.

The null model is fitted by iteratively reweighted least squares, which
turns each step of the logistic fit into a weighted linear one: with `mu`
the fitted chance for each individual, the weight is `mu (1 - mu)`, and the
step solves the design weighted by those against the difference between the
trait and `mu`. It stops when the largest change in a coefficient is below
1e-8, in at most 50 steps, and a fit that has not converged by then is an
error, not a warning.

The **score test**, which needs only the null, tests every variant of a
block at once. With `w` the weights, `resid` the trait minus `mu` and `d`
the design:

    num = x' resid
    den = x' w x - (x' w d) (d' w d)⁻¹ (d' w x)

and then `beta = num / den`, `se = 1 / sqrt(den)` and a chi square with one
degree of freedom of `num² / den`. The covariates take the place of the
projection matrix of the mixed models, and nothing is inverted per variant.

The **Wald test**, the default, fits one logistic regression per variant
with the variant in the model, starting from the null's coefficients and an
effect of 0 for the variant. Each step needs one system of `c + 1` unknowns
per variant, which is a Cholesky factorization and a solve of a matrix the
size of the coefficients; `docs/specs/linalg.md` decided that a stack of
those is a loop in the caller and not an operation of the crate, and
measured one 7 x 7 factored and solved at 0.173 µs, so a block of 5000
variants costs 0.9 ms for one step. `se` is the square root of the last
diagonal entry of the inverse of that system's matrix, and the p-value is a
chi square with one degree of freedom of `(beta / se)²`.

A variant whose Wald fit runs away gets NaN for all three. Three things mark
it, and popnei reproduces all three: a step that is not finite, which
includes the system the factorization refuses as singular; a coefficient
whose absolute value passes 30; and a fit still moving after 50 steps.

The first of the three is tested as **not finite** and not as an infinity,
which is what keeps the two backends of `docs/specs/linalg.md` giving the
same answer. That crate lets through a diagonal entry that is neither 0 nor
below it but whose reciprocal overflows, and there LAPACK gives an infinity
where faer gives a NaN, each reporting success; measured by that spec on the
2 x 2 with 4e-309 and 1 on its diagonal. Both fail a test for a value that
is not finite, so the variant gets its three NaNs either way, and Python
natively, Python under pyodide and TypeScript mark the same variants. If the
owner ever has the crate refuse that entry instead, this module gets a
`Singular` where it now gets an infinity, and that is already one of the
three marks. The
threshold of 30 is inherited from pyNei and nobody has measured it. A
variant that separates the cases from the controls perfectly has no finite
effect and is what these catch: the panel has exactly one, `var0006`.

plink2 does not give up on that variant. It falls back to a Firth penalized
regression, which adds a term that pulls the estimate back from infinity and
gives a finite answer; its `FIRTH?` column says `Y` for that one variant and
`N` for the other 1199. popnei gives NaN, as pyNei does (**Open 1**, below).

### How it is verified

Against plink2 `--glm hide-covar` on the panel with every genotype called,
which writes `plink2.panel_called.glm.logistic.hybrid.tsv`. plink2 reports
the odds ratio, so `beta` is compared with its logarithm.

The one variant plink2 fell back to Firth for is left out of the comparison,
and instead a test asserts that popnei's NaNs are exactly the variants
plink2 marked `FIRTH?` `Y`, which is `var0006` and no other. Over the other
1199: `beta` and `se` within 1e-4 times the `se` of that variant and
`p_value` within 5e-3 relative. All three are wider than plink2's printing, which rounds by 5e-6,
because plink2 stops its logistic fit earlier than popnei does; the
difference between the two fits is what these measure, and the printing is
not what limits them. The six literals below are held to 1e-5 on `beta`,
1e-4 on `se` and 5e-3 on `p`, the same as pyNei holds them.

| variant | beta, a log odds ratio | se | p |
|---|---|---|---|
| var0000 | -0.5725786945415258 | 0.261917 | 0.0288081 |
| var0052 | -0.8528230964300417 | 0.248979 | 0.000614166 |
| var0629 | -0.949570924908968 | 0.323553 | 0.00333739 |
| var0751 | -0.2659166667836026 | 0.219207 | 0.225098 |
| var1137 | -0.42744848150358195 | 0.252402 | 0.0903561 |
| var1188 | -0.830184139078324 | 0.264071 | 0.00166771 |

`beta` carries more digits than plink2 prints because plink2 gives the odds
ratio and these are its logarithm.

The score test against R 4.6.1's `anova(glm, test = "Rao")`, one logistic
regression per variant fitted by R, from `r.panel_called.glm.score.tsv`.
R reports the score statistic and its p-value at full precision. Over all
1200 variants the statistic `(beta / se)²` is within 1e-2 absolute and
`|log10(p / p_R)|` below 1e-3; the six literals are held to 1e-3 and 1e-3.
R's glm converges to 1e-8 in the deviance, which is what those tolerances
are. The one on the statistic is absolute where the reference is exact, so
on a panel whose statistics are far above this one's 1.5 to 12 it would fail
a right answer rather than pass a wrong one, which is the safe way round for
a check to be fragile.

The six literals, the score statistic `(beta / se)²` within 1e-3 absolute
and the p-value within 1e-3 in `log10`: var0000 4.938245 and
0.026268700, var0052 12.484427 and 0.000410359, var0629 9.165576 and
0.002466100, var0751 1.480401 and 0.223711736, var1137 2.911424 and
0.087954199, var1188 10.382961 and 0.001271835.

## The logistic mixed model

### What it gives

A binomial trait with a kinship, and the only test it has is the score test.

It is fitted by **penalized quasi-likelihood**: for a fixed variance of the
kinship effect, `tau`, the 0/1 trait is turned into a continuous **working
trait**, each individual carrying a weight that says how much its 0 or 1
tells us at the fit so far, and a weighted linear mixed model is fitted to
that working trait; the working trait and the weights are then made again
from the new fit, and so on. One pass of that is a **linearization**. `tau`
then takes one Newton step from the restricted maximum likelihood. A Newton
step needs the second derivative of the likelihood, and the **average
information** is the average of the observed one and the one expected under
the model: the terms that cost the most to compute appear in the two with
opposite signs and cancel, so the average costs less than either, which is
why mixed model programs use it. With `p` the projection matrix and `w` the
working trait of the linearization that just finished, and `pw = p w`:

    score = 0.5 * (pw' k pw - trace(p k))
    ai    = 0.5 * (k pw)' p (k pw)
    step  = score / ai

and `tau` becomes `tau + step`. The whole thing then starts again. It is the model GMMAT fits, and on the panel it takes 8
steps on `tau` and 22 linearizations.

It starts from the plain logistic null of "The logistic model", fitted with
no kinship in it: its coefficients give the first linear predictor and the
first `mu`, and the first working trait and weights come from those. That
fit in turn starts with every coefficient at 0 but the intercept, which
starts at `log((m + 1e-6) / (1 - m + 1e-6))` with `m` the mean of the trait.
An implementation that starts the mixed fit anywhere else walks a different
path through the bracket and can stop at another `tau`.

Taking the step on `tau` after every single linearization instead makes the
two updates fight and `tau` cycle for ever, which pyNei's `_GLMMNull` records
and a trial here reproduced. `tau` is kept from cycling by a bracket: a
`tau` whose derivative asks for a larger one and a `tau` whose derivative
asks for a smaller one are remembered, the answer lies between them, and a
Newton step that would leave that interval is replaced by the geometric mean
of its two ends. `tau` at the boundary, where the kinship explains nothing,
is 0 and the fit stops there.

The tolerances and the counts, all inherited from pyNei, which calls them
`GLMM_TOL` and `GLMM_MAX_ITER`, and none of which anybody has measured. A
linearization stops when the largest change in the linear predictor, over
its own largest absolute value plus 1, falls below 1e-6. A step on `tau`
stops the fit when its absolute value falls below `1e-6 * (tau + 1e-6)`.
Either loop running past 200 rounds raises, naming the model and the number
of rounds. pyNei raises a `RuntimeError` there, and so does its logistic
fit; popnei raises a `ValueError`, because under the rule of
`docs/specs/variant.md` a `RuntimeError` is a defect of popnei and a fit
that will not settle is the data, the same category as the kinship that is
not a covariance below. The 50 steps of the logistic fit are the same. `tau` starts at half the variance of the first working trait. The whole
variance is what `tau` would be if the kinship explained all of it, and
pyNei halves it; either way the start is too large rather than too small, so
the bracket comes down from above. A step that would
take `tau` to 0 or below, before both ends of the bracket are known,
quarters it instead; a `tau`
that falls below 1e-6 is set to 0, and a second one at 0 ends the fit.

At convergence the residual `p y` of the score test is simply the trait
minus `mu`, and `den` is `x' p x` with the same projection matrix as the
linear mixed model. `residual_variance` and `heritability` are `None`: a
logistic model has no free residual variance.

### How popnei fits it, and why not pyNei's way

pyNei inverts an individuals by individuals matrix once per linearization,
22 times on the panel and 25 at 4000 individuals. popnei does not. The
measurements are in `docs/reports/glmm-method/README.md`, and the two things
that change are:

- The covariance `sigma = tau * k + w⁻¹`, with `w⁻¹` the weights on the
  diagonal, is factored with a Cholesky and applied by solving, never
  inverted. A Cholesky costs a third of an inverse: 0.126 s against 0.376 s
  at 4000 individuals, numpy 2.5.3 on Accelerate on the owner's Apple M5
  Pro.
- The one quantity that seemed to need every entry of the inverse, the trace
  of `p k`, comes from an identity. `tau k = sigma - w⁻¹`, so
  `trace(sigma⁻¹ k) = (n - trace(sigma⁻¹ w⁻¹)) / tau`, and
  `trace(sigma⁻¹ w⁻¹)` is the sum of the squares of the entries of
  `l⁻¹ w^-1/2`, the Cholesky factor solved against with one right hand side
  for each individual. It is wanted once per step on `tau`, 7 to 9 times
  over a fit, not once per linearization.

The inverse is formed once, at the end, because the score test wants the
projection matrix as a matrix.

The two fits take the same steps in the same order and stop at the same
place; only the arithmetic of each step differs. `tau` agrees with pyNei's
to 2.9e-15 relative on the panel, the covariate effects to 3.3e-15, the
projection matrix to 1.4e-15, and the p-values of all 1200 variants give the
same largest `|log10(p / p_GMMAT)|` to every digit printed, 8.497e-06. The
null fit takes 0.095 s against pyNei's 0.195 s at 1000 individuals, 0.608
against 1.199 at 2000, and 5.285 against 9.959 at 4000.

The owner asked for a cheaper fit on 23 September 2026, before this spec was
written, and on being shown the measurements directed the linear algebra
crate to add the lower triangular solve that this route needs, which is the
decision to build the module around it. The options not taken, all measured
in that report: an
eigendecomposition of the weighted kinship per linearization, which would
make the search over `tau` cost one number per individual but costs 3.33 s
at 4000 individuals against an inverse's 0.376; a conjugate gradient solve, which
solves a system by repeated products of the matrix with vectors and never
factors it, and which wins for a sparse kinship and loses about twofold for
popnei's dense one; and a stochastic estimate of the trace, which could give no more than a
further 1.7 because the 25 Cholesky factorizations are 3.1 s of the 5.3 and
which would stop the fit being the same calculation twice. What is not known
is where the ratio settles above 4000 individuals, where nothing was run,
and that it has not been measured in Rust: both sides here are numpy on the
same BLAS.

### A kinship that is not positive semidefinite

A Cholesky factorization refuses a matrix that is not positive definite
where pyNei's inverse carried on, and the per pair denominators of
`docs/specs/kinship.md` can leave the kinship with an eigenvalue below 0. It
does not stop this fit at any missing rate tried. On 400 individuals and
2000 variants, the smallest eigenvalue of the kinship runs from -0.033 at 3
genotypes missing in 100 to -1.06 at 50 in 100, and both fits succeed at
every rate, giving a `tau` between 1.26 and 1.48. The reason is that a
weight is at most 0.25, so `w⁻¹` puts at least 4 on every diagonal entry of
`sigma`, and `sigma` only goes indefinite once `tau` passes about 3.8 even
at 50 in 100.

If a dataset ever reaches it, the factorization gives the `Singular` of
`docs/specs/linalg.md` with the row it stopped at, and `calc_gwas` turns it
into an error naming the kinship and saying that missing genotypes can make
one that is not a covariance. It is not a defect of popnei and not a wrong
argument, so it is neither a `RuntimeError` nor a plain `ValueError` about a
type: it is a `ValueError` about the data.

### How it is verified

Against GMMAT's `glmm.score` on both panels, from
`gmmat.panel_called.glmm.score.tsv` and `gmmat.panel.glmm.score.tsv`, to the
same tolerances as the linear mixed model's score test, for the whole
columns and for the six literals alike: `1 / se²` against `VAR` within 1e-5
relative and `|log10(p / p_GMMAT)|` below 1e-4.

| variant | VAR | p | VAR, missing | p, missing |
|---|---|---|---|---|
| var0000 | 6.48664 | 0.702659 | 6.43005 | 0.685719 |
| var0052 | 8.98834 | 0.0306703 | 8.73718 | 0.0291759 |
| var0629 | 6.49956 | 0.0895104 | 6.19348 | 0.12319 |
| var0751 | 10.563 | 0.0142331 | 10.1532 | 0.026766 |
| var1137 | 9.098 | 0.093808 | 8.92384 | 0.108032 |
| var1188 | 9.05095 | 0.0262737 | 8.7063 | 0.0196806 |

The null model against `glmmkin`, within 1e-5 absolute: `genetic_variance`
1.508057, `residual_variance` `None`, and the covariate effects -1.416464,
0.753476 and 1.583210.

## The GRAMMAR-Gamma approximation

### What it gives

`x' p x`, the denominator of both mixed model tests, is a product of the
dosages of a block with an individuals by individuals matrix, so it costs
work proportional to the square of the individuals for every variant. The
approximation replaces it with `gamma` times the squared length of the
variant's centered dosages, which is linear in the individuals, with one
`gamma` estimated once from the first variants that vary.

`gamma` is the mean, over the first 100 variants that vary of the first
block, of the exact `x' p x` divided by the approximate one. It is
`NUM_VARS_FOR_GAMMA` in pyNei and 100 here, inherited, and nobody has
measured whether 100 is the right number.

It is `use_grammar_gamma_approx=True`, false by default, and asking for it
without a kinship is a `ValueError`, since there is no projection matrix to
approximate.

What it costs in accuracy grows with how strongly the panel is structured,
because one `gamma` stands in for a quantity that really differs from
variant to variant. On the panel, `test_grammar_gamma_approx` of pyNei
asserts that the median of `log10(p_approx / p_exact)` is within 0.1 of 0
and the largest is within 1.5, so a p-value can be out by a factor of 30 in
the worst case while the middle of the distribution barely moves.

### How it is verified

There is no program outside the project to check it against: GMMAT and
rrBLUP compute the exact denominator. What is checked is the relation to
popnei's own exact answer, on the panel with the `lmm` and both covariates,
which is pyNei's test above: the median and the largest of
`log10(p_approx / p_exact)`, and that `beta` agrees with the exact one
within 0.5 relative. Also that `used_grammar_gamma_approx` is in the result
and that asking for it without a kinship raises.

## The two distributions

### What they give

Two functions turn a statistic into a p-value, and numpy has neither, which
is why pyNei wrote both.

- `chi2_sf_1df(x)`, the chance that a chi square with one degree of freedom
  is above `x`, which every score test and the logistic Wald test need. It
  is `erfc(sqrt(x / 2))`, the complementary error function, which gives how
  much of a normal distribution lies past a point. An `x` of 0 or below
  gives 1.0, as scipy's `chi2.sf` does, and not the NaN that the square root
  of a negative number would give. Who may pass one is **Open 2**, below.
- `t_sf_two_sided(t, df)`, the chance that a Student t with `df` degrees of
  freedom is further from 0 than `t`, which the linear model and the linear
  mixed model's Wald test need. It is the regularized incomplete beta
  function `I_x(df/2, 1/2)` at `x = df / (df + t²)`, and it hands the
  incomplete beta `t² / (df + t²)` beside it as `one_minus_x`, for the
  reason that section gives.

popnei takes `erfc` from the `libm` crate, a pure Rust port of musl's math
library with no C in it, which builds for both wasm targets, checked as a
library on 23 September 2026. Measured against scipy 1.18.1's `chi2.sf` over
65 points spaced logarithmically from 1e-6 to 200, `libm`'s `erfc` is within
2.9e-14 relative, and within 3.0e-14 over 65 evenly spaced ones, and
at `x = 100` it gives 1.523971e-23, so the tail is right where a strong
variant needs it. The option not taken was to write `erfc` here: pyNei takes
it from Python's `math`, which is C's, so pyNei is not a precedent for
writing one, and reaching 1e-12 relative in a tail of 1e-23 by hand is work
with no reward.

The regularized incomplete beta is not in `libm` and is written here, as
pyNei writes it: the continued fraction of Numerical Recipes evaluated by
Lentz's method, which builds a continued fraction from its front rather than
from its far end, so it can stop as soon as a term no longer changes the
value instead of needing its depth fixed in advance.

It takes **both** `x` and `one_minus_x` from its caller and never subtracts
one from the other. `t_sf_two_sided` has them for nothing, `df / (df + t²)`
and `t² / (df + t²)`, and computing the second as `1 - x` instead throws
away every digit of it once `x` has rounded to 1: at 197 degrees of freedom
and `t` of 1e-7 that returned exactly 1.0 where the answer is
0.9999999203127337. Measured on 23 September 2026, taking `one_minus_x` from
the caller moved the worst relative error over `t` in [1e-7, 1e-3] from
7.97e-8 to 5.34e-17 at 197 degrees of freedom, and from 6.34e-7 to 3.83e-15
at 9997.

An `x` at or below 0 gives 0 and a `one_minus_x` at or below 0 gives 1.
Otherwise, with

    front = exp(lgamma(a + b) - lgamma(a) - lgamma(b)
                + a * ln(x) + b * ln(one_minus_x))

the answer is `front * cf(a, b, x) / a` while `x` is below
`(a + 1) / (a + b + 2)`, and `1 - front * cf(b, a, one_minus_x) / b` at or
above it, which is the symmetry `I_x(a, b) = 1 - I_{1-x}(b, a)` used where
the fraction converges slowly. pyNei writes both of those from `x` alone,
`numpy.log1p(-xi)` at `src/pynei/gwas.py:329` and `1 - xi` at 337, and
`log1p` recovers the logarithm but not the fraction's argument. `cf` is the continued fraction, with `tiny` at
1e-300, `eps` at 1e-15 and at most 500 rounds:

    qab = a + b;  qap = a + 1;  qam = a - 1
    c = 1;  d = 1 - qab * x / qap;  if |d| < tiny then d = tiny;  d = 1 / d
    h = d
    for m in 1 ..= 500:
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d;  if |d| < tiny then d = tiny
        c = 1 + aa / c;  if |c| < tiny then c = tiny
        d = 1 / d;  h = h * d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d;  if |d| < tiny then d = tiny
        c = 1 + aa / c;  if |c| < tiny then c = tiny
        d = 1 / d;  delta = d * c;  h = h * delta
        if |delta - 1| < eps then stop
    cf = h

`tiny` keeps a denominator that has come out at 0 from dividing, which is
what Lentz's method needs to carry on past a term that vanishes, and running
out of rounds is not an error: the four pairs of arguments this module uses
converge in at most 52 rounds, measured on 23 September 2026 over the sweep
below.

**`tiny` cannot fire at any argument the t distribution reaches, and it
stays.** The first denominator is `d = 1 - (a + b) x / (a + 1)`, and in both
branches it is bounded below by `2 / (a + b + 2)`: in the direct branch
because `x` is below `(a + 1) / (a + b + 2)`, in the symmetry branch because
`1 - x` is at most `(b + 1) / (a + b + 2)`, and the zero of `d` lies above
both. So `d` reaches 1e-300 only once `a + b` passes about 2e300, which with
`b` at 1/2 and `a` half the degrees of freedom is a panel of 4e300
individuals. The bound is tight and two measurements meet it: over 6009003
calls, degrees of freedom 1 to 3000 and then 1e4, 1e5 and 1e6 with `t` from
0 to 20, the smallest `|c|` or `|d|` was 4.027585806198886e-6, which is
`2 / (a + b + 2)` exactly at a million degrees of freedom. It stays because
the recipe and pyNei have it and because a caller with some other `b` would
need it; `b` here is always 1/2.

`eps` is a different thing: it caps the work and does not get the digits.
With it set to 0, so that the loop always runs its 500 rounds, nothing
became not finite and the worst value moved by 2.3e-13 relative. The
fraction took at most 52 rounds of its 500, over a sweep of 116802 calls.

So no test of either guard can fail on a value, and the only assertion with
anything behind it is one on the number of rounds. A reader who finds the five lines that read
`tiny`, and the one that reads `eps`, covered by no test should stop looking
for the argument that reaches them: for `tiny` there is none, and the bound
above says why.

### How it is verified

Against scipy 1.18.1, whose numbers go into the cargo tests as literals, at
`chi2_sf_1df` and `t_sf_two_sided` of "The Rust interface", which are public
for this reason. The cases are pyNei's, in `test_distributions`:

- The incomplete beta at the four pairs `(0.5, 0.5)`, `(10, 0.5)`,
  `(98.5, 0.5)` and `(2.5, 7)`, over `x` drawn uniformly in (0, 1), within
  1e-12 absolute. The pair `(98.5, 0.5)` is what a t with 197 degrees of
  freedom uses, next to the panel's 196: 200 individuals less the three
  columns of its design less one for the variant.
- `t_sf_two_sided` at 5, 17, 197, 997 and 9997 degrees of freedom, over a
  spread of `t` including 10, 20 and 40, within 1e-10 relative.

**The 1e-10 is claimed to 9997 degrees of freedom and not beyond**, which is
the 10000 individuals `docs/objectives.md` names, less the coefficients and
the variant. The error grows with the degrees of freedom and its worst point
is not spread over `t`: it sits at `t` near 1.73, where the branch of the
incomplete beta switches. Measured against mpmath at 60 digits on 23
September 2026: 4.7e-13 relative at 197 degrees of freedom, 7.7e-13 at 997,
5.5e-11 at 9997, 1.5e-10 at 20000 and 3.5e-9 at 500000. So the bound has 213
times the room at 197 degrees of freedom, 1.8 times at 9997, and is
already untrue at 20000 individuals. It is the cancellation in
`lgamma(a + b) - lgamma(a)`, amplified by the `1 -` of the symmetry branch,
and not the stopping rule: setting `eps` to 0 moves the worst point from
3.3326e-9 to 3.3324e-9.

That 1.8 is the one bound of this spec sitting near its failure, and it is
stated rather than widened because widening it would catch less at the sizes
popnei actually runs. Whoever first wants popnei past 10000 individuals has
to come back to this function before they can trust its p-values, and the
front factor in logarithms is where to start.
- `chi2_sf_1df` over a chi square sample and at 30, 50 and 100, within 1e-12
  relative.

How far those three bounds are from the differences they allow, measured on
23 September 2026 on a built implementation: the chi square's worst is
1.8e-14 relative, 57 times inside its bound; the incomplete beta's worst is
8.5e-15 absolute at the pair `(98.5, 0.5)`, 117 times inside; and the
Student t's worst is 4.7e-13 relative at 197 degrees of freedom, 213 times
inside, at `t` near 1.73 where the branch switches. So these three are not
the round numbers that
"How it is verified" of "What every model shares" warns about, and they do
not need lowering to where they break: the room has been measured and it is
there.

popnei's and pyNei's p-values differ by the difference between two `erfc`
implementations, about 1e-14 relative, which is five orders below the 1e-9
the two libraries are compared within for a whole study.

## The Rust interface

What a study is given. The design is `num_individuals` x `num_coefs`, row
after row, with its column of ones already in it, because the Python and
TypeScript layers are what turn a user's frame into one. `kinship` is
`num_individuals` x `num_individuals`, row after row, already cut to the
individuals that are tested and in their order.

```rust
pub enum TraitType { Continuous, Binomial }
pub enum TestType { Wald, Score }
pub enum GwasModel { Lm, Lmm, Glm, Glmm }

pub struct GwasInput<'a> {
    /// One value per tested individual: the measurement, or 0.0 or 1.0.
    pub phenotype: &'a [f64],
    pub trait_type: TraitType,
    /// num_individuals x num_coefs, row after row, the intercept first.
    pub design: &'a [f64],
    pub num_coefs: usize,
    pub kinship: Option<&'a [f64]>,
    /// None takes the default for the trait and the kinship.
    pub test: Option<TestType>,
    pub use_grammar_gamma_approx: bool,
    /// The positions of the tested individuals among those the reader
    /// gives, in the order `phenotype` and `design` have them.
    pub individuals: &'a [usize],
    pub transform_to_biallelic: bool,
}
```

`kinship`, when it is given, is checked before any model is fitted: that it
holds `individuals.len()` times `individuals.len()` values, and that every
one of them is finite. Neither is a thing a fit would notice. A matrix of
the wrong length is read as another shape and gives numbers, and a NaN in
one comes back much later as the linear algebra crate's refusal of a value
that is not finite, naming a matrix at whichever routine met it first.

What a study gives back. `beta`, `se` and `p_value` hold NaN for a variant
that has no answer.

```rust
pub struct NullModel {
    pub model: GwasModel,
    pub test: TestType,
    /// One per column of the design.
    pub covariate_effects: Vec<f64>,
    pub residual_variance: Option<f64>,
    pub genetic_variance: Option<f64>,
    pub heritability: Option<f64>,
    pub num_individuals: usize,
}

pub struct Gwas {
    pub num_vars: usize,
    pub null_model: NullModel,
    pub allele_freq: Vec<f64>,
    pub beta: Vec<f64>,
    pub se: Vec<f64>,
    pub p_value: Vec<f64>,
    pub used_grammar_gamma_approx: bool,
    /// The interned chromosome of each variant, read through
    /// `chrom_table`, and `None` when the source had no such column.
    pub chroms: Option<Vec<u32>>,
    pub chrom_table: ChromTable,
    pub poss: Option<Vec<u64>>,
    pub ids: Option<Vec<String>>,
}
```

The three columns are what the Python and TypeScript layers build `chrom`,
`pos` and `id` of `stats` from, and they are laid out as `docs/specs/ld.md`
lays out the same three. A VCF and a vars file both carry the chromosome and
the position, so in practice only `ids` is ever `None`.

One pass over a reader, or two when `use_grammar_gamma_approx` is true, the
second being opened over the same variants as the PCA's is. The pass borrows
its readers and does not take them, asks for the genotypes and for `chrom`,
`pos` and `id`, and puts `reblock` before each.

```rust
pub fn calc_gwas<R1: BlockReader, R2: BlockReader>(
    reader: &mut R1,
    gamma_pass: Option<&mut R2>,
    input: &GwasInput<'_>,
) -> Result<Gwas>;
```

The two distributions, public so that the cargo tests check them against
scipy's numbers where the value can be seen.

```rust
/// The chance that a chi square with one degree of freedom is above `x`.
pub fn chi2_sf_1df(x: f64) -> f64;

/// The chance that a Student t with `df` degrees of freedom is further
/// from 0 than `t`, both tails.
pub fn t_sf_two_sided(t: f64, df: f64) -> f64;
```

What this module calls in `linalg`: the thin QR of the design and the solve
against its upper triangular `r`; the Cholesky factorization and its solve,
log determinant and inverse; `solve_triangular` reading the lower half, for
the trace of the logistic mixed model; the rank of the design; the
eigendecomposition of the kinship; and the product, in all four of its
combinations.

Three things about that crate the fits have to know. It refuses what it is
given and not what it produced, so a solve or an inverse off a covariance
that is positive definite and nearly not can come back `Ok` holding an
infinity, or a NaN on the other backend; noticing a fit that has run away is
this module's job, which is what the three marks of the logistic Wald test
do, and testing them for a value that is not finite rather than for an
infinity is what makes the two backends agree. The rank uses numpy's
tolerance, which is what makes a design popnei refuses a design pyNei
refuses. And the two backends agree on the rank between about 1e-300 and
1e154; a design of dosages and covariates is nowhere near either end, since
a dosage is 0 to the ploidy, but popnei does not scale the covariates a user
gives, so a covariate in extreme units could get there.

## Speed

From the table of section 2.1 of `docs/rust_core.md`, over 100000 variants x
1000 individuals, the linear model takes pyNei 0.35 s and plink2 0.10 s. The
number to reach for the `lm` is plink2's 0.10 s on that dataset.

For the mixed models, section 2.4 of `docs/rust_core.md` reports pyNei level
with GMMAT on one thread over the same dataset, 1.5 s against GMMAT's 1.6 s
and 2.2 s, and 3x faster with six threads. It does not say which of GMMAT's
two numbers belongs to which mixed model, so the number to reach is the
smaller, 1.6 s, for both. Whether popnei beats it is not known: the null fit
of the `glmm` is 1.88 to 2.06 times cheaper than pyNei's, measured, but at
100000 variants and 1000 individuals the per variant work is the larger
half, 0.32 s against a fit of 0.095 s by the numbers below, and that work is
the same product in both libraries.

Where the time goes at many individuals was measured for this spec, numpy
2.5.3 on Accelerate on the owner's Apple M5 Pro, with the score test of
100000 variants beside the null fit it feeds: at 1000 individuals the fit is
0.095 s and the test 0.32 s, at 2000 it is 0.608 s and 1.26 s, and at 4000 it
is 5.285 s and 5.13 s. So the fit and the test are of the same order once
the individuals reach a few thousand, and `use_grammar_gamma_approx` is what
addresses the test half.

None of this has been measured for popnei; the measurements come when the
code exists, on the panel and on the 100000 x 1000 dataset of
`docs/rust_core.md`.

## Open points

The owner decides these two, and until then the implementer follows the
"meanwhile" of each.

**Open 1: a variant that separates the cases from the controls.** Its
logistic effect is infinite and its Wald fit runs away. pyNei gives NaN for
`beta`, `se` and `p_value`; plink2 falls back to a Firth penalized
regression, which adds a term that pulls the estimate back from infinity,
and reports a finite answer marked `FIRTH?` `Y`. The panel has exactly one
such variant of 1200, `var0006`. The options are to reproduce pyNei and give
NaN, which loses a variant that plink2 reports and which a user cannot tell
apart from a variant with no variance, since both are three NaNs; to give
NaN but say which variants they were, a count in the result or a column
saying why each NaN is there, which costs one field and tells the user where
to look; or to implement the Firth regression, which is a second fitting
method for one variant in a thousand and which no part of popnei needs
otherwise. Recommendation: the middle one, NaN with a reason. The numbers
stay pyNei's and plink2's comparison stays as it is, and a user who sees a
variant vanish learns whether it had no variance or a runaway fit, which are
different things to do something about. Meanwhile the implementer gives NaN
with no reason, as pyNei does, since no literal of this spec moves either
way and the column can be added without changing a number.

**Open 2: a variant whose score denominator rounds to 0 or below.** `x' p x`
is 0 or above in exact arithmetic, and rounding can put it just below for a
variant with almost no variance left once the covariates and the kinship are
taken out. The options are to give that variant three NaNs, as a variant
with no variance gets, which says the study could not test it and throws
away a `beta` that is meaningless anyway; or to let the statistic through to
`chi2_sf_1df`, which gives 1.0 for an argument of 0 or below, so the user
sees a p-value of 1 beside a `beta` and an `se` that are large and wrong.
Recommendation: three NaNs, refused before the statistic is formed, at
`den <= 0`. A p-value of 1 is a claim that the variant was tested and showed
nothing, and nothing was tested; and "The variants that have no answer"
already means the three NaNs together, so a NaN p-value beside a finite
`beta` would be a fourth thing a user has to learn to read. It costs the
distinction between a variant with no variance at all and one whose variance
the null model absorbed, which no reference program reports either. This
came from the session building `gwas-linear` on 23 September 2026, which met
it in the score test of the linear mixed model. Meanwhile the implementer
refuses at `den <= 0` and gives the three NaNs; if the owner chooses the
other, the change is one comparison and no literal of this spec moves, since
no variant of either panel reaches it.

## Not in this spec

- The kinship itself, its per pair denominators and its principal
  components: `docs/specs/kinship.md`. This module only takes one.
- The dosage of a genotype, the major allele and what a missing one gets:
  `docs/specs/pca.md`, which this spec takes unchanged and does not repeat.
- The seven operations of linear algebra, their backends, their errors and
  what each one costs: `docs/specs/linalg.md`.
- The measurements behind popnei's fit of the logistic mixed model, the
  three routes that were tried and dropped, and what is not known about it:
  `docs/reports/glmm-method/README.md`.
- Multiple testing. popnei gives a p-value per variant and no Bonferroni,
  no false discovery rate and no genomic control, because pyNei has none and
  a user applies their own to the column. If popnei ever adds one it gets
  its own item here.
- A joint model of a multiallelic variant, one row per variant and allele,
  which `docs/rust_core.md` leaves open. This module collapses every allele
  that is not the major one, as `transform_to_biallelic` says.
- Fitting the logistic mixed model without a dense individuals by
  individuals factorization at all, which needs a sparse kinship and a
  different algorithm: `docs/reports/glmm-method/README.md` says what was
  measured and why it was not taken.
- Interactions between a variant and a covariate, and testing several
  variants together. pyNei has neither and popnei does not add them.
