# The pca module: principal components of variants and of a table

September 2026. A principal component analysis places the individuals of a
dataset on a few axes that hold as much of the variation between them as
that many axes can, which is how a user sees whether their individuals
fall into populations before they name any. This module gives two: the
PCA of the variants, from the genotypes of a `Variants`, and the PCA of
any table of numbers, individuals by traits. It gives as well the
principal coordinate analysis, which places the individuals on components
from the distances between them, added on 27 September 2026: of a
`Distances`, and of the Kosman distances of a `Variants`. The two PCAs are
built; the principal coordinates have no code. This spec develops the
row `pca` of the table in section 9 of `docs/architecture.md`.

It also holds the answer to a question the owner asked with it: whether
keeping the genotypes in 2 bits, as plink does, or using the vector
instructions of the processor, would make the PCA of the variants faster.
The 2 bit coding was measured and is not taken, and the vector
instructions are taken in two places, neither of which needs `unsafe`. The
numbers are under "Speed".

It depends on `docs/specs/block.md`, for the `Block`, the run of
consecutive variants held as arrays, the `BlockReader` trait of everything
that gives blocks, and `reblock`, which puts the blocks of a reader back
to one size; on `docs/specs/variant.md`, for the `Variants`, a source of
variants with the filters that were put on it as steps, which a
calculation runs in as many passes as it needs, and for the `PassStats`
that every result of such a calculation carries; and on two functions of
the `linalg` module, which has no spec. The owner decided on 21 September
2026 that `linalg` gets a spec of its own, written next, and that this
spec gives only what it calls there. The options not taken were to
specify the two functions here, and to use faer, the linear algebra
library in pure Rust, natively as well as in wasm until a BLAS backend is
specified.

## What both analyses compute

The data is a table X of n rows, the individuals, and p columns, the
traits, or the variants once each is a number per individual. Each column
is centered, its mean taken from it, and standardized, divided by its
standard deviation, which gives Z. The principal components are the
directions in the space of the columns, each a vector of p weights of
length 1, along which the rows of Z vary most: the first has the largest
variance that any direction has, the second the largest among the
directions at a right angle to the first, and so on. The result has three
things, under pyNei's names:

- **projections**, n x components: where each individual falls along each
  component, the product of its row of Z with the weights.
- **explained_variance_percent**, one per component: the variance of the
  projections on that component, as a percentage of the sum of the
  variances of all the components, given or not, which is the sum of the
  variances of the columns of Z.
- **princomps**, components x p: the weights of each column in each
  component.

pyNei takes the singular value decomposition of Z, which needs Z whole in
memory. popnei takes the eigenvectors of the smaller of the two matrices
of products of Z with itself. When n is the smaller side, that is G = ZZ',
the product of Z with its transpose, n x n, whose entry i, j is the sum
over the columns of the value of individual i times that of individual
j. With λ_j its eigenvalues from the
largest, and u_j its eigenvectors:

    projections of component j          = u_j * sqrt(λ_j)
    explained_variance_percent of j     = 100 * λ_j / sum of every λ of G
    princomps of component j            = Z' u_j / sqrt(λ_j)

When p is the smaller side the matrix is Z'Z, p x p, its eigenvectors are
the princomps themselves and the projections are Z times them. Both give
what the SVD gives. Measured with numpy against pyNei's `do_pca`, on iris
and on the dosages of the reference panel of "How it is verified", 200
individuals x 1200 variants: the projections differ by 5e-12 at most and
the percentages by 3e-14, over every component that has variance.

G is a sum over the columns, so for the variants it is added up block by
block, and nothing of size variants x individuals is ever held.

Three rules hold for both analyses.

**The sign.** A component multiplied by -1 is the same component, and
which of the two comes out depends on the library that did the
eigendecomposition; pyNei gives whichever LAPACK gave, and its tests
accept both. popnei has two such libraries, so it fixes the sign: in every
component the projection with the largest absolute value is positive, and
when two individuals have the same absolute value it is that of the first
of them; the princomps of that component take the same sign. Python natively, Python under pyodide and TypeScript then
give the same numbers.

Two projections have the same absolute value when they are within 64
times 2.2e-16 of the larger of the two, which is 1.4e-14 of it. Without
that tolerance the rule is decided by the last bit, and the last bit is
what the two libraries differ in, so it gives opposite signs on the same
data. Six individuals of three variants, in which the individuals 2 and 5
have the same genotype at every variant and the individual 0 has the
opposite dosage at every one, have the three projections of the same
absolute value, 1.98902864125, in exact arithmetic: Accelerate's LAPACK
gives that of the individual 5 one unit in the last place above the other
two, so the largest of the three is positive already and the individual 0
comes out at -1.989, and numpy's eigendecomposition spreads the three over
five units in the last place with that of the individual 0 the largest, so
it turns the component round and the individual 0 comes out at +1.989. The
two differ by 3.98 in the projections of that individual and by twice each
weight in that component. With the tolerance the three are one absolute
value, the first of them decides, and the individual 0 is positive on
both.

The tolerance is 64 units in the last place and not a handful because the
spread of a tie grows with the matrix: the same tie came out over one unit
on Accelerate and five on numpy at six individuals, and nothing has been
measured at a thousand. It is 1.4e-14 of the largest projection, five
orders of magnitude below the 1e-9 within which the results are compared,
so no pair of projections it takes for one is a pair a user could tell
apart. What it does not fix is two projections that are not equal in exact
arithmetic and come out closer than the tolerance on one library and
further on another: a fixture of six individuals gave two that differ by
18 units in the last place, which is outside it, and there the larger
decides as it did before.

`tests/reference/pca/make_reference.py` applies the same rule with
numpy's floats and exact equality, so a reference dataset with a tie would
need this tolerance there too. None of them has one.

**A component with no variance is not given.** Centering takes one
dimension out of the data, so a table of 8 individuals and 30 variants
has 7 components and not 8. pyNei gives `min(n, p)` of them. In `worked3`
of the reference data, 5 individuals and 4 variants, it gives 4: the last
has a percentage of 3.5e-32, projections of 1e-16, and weights that are
whatever the SVD left there, and R's four weights for it all have the
opposite sign. In
popnei its weights would come from dividing by sqrt(λ), a number that is
about 0, and would be noise of any size. popnei gives
the components whose eigenvalue is above λ_1 * max(n, p) * 2.2e-16. That
is the tolerance numpy's `matrix_rank` has for singular values, used here
on eigenvalues, which are their squares, on purpose: the error of an
eigenvalue of G is of that size, and the square of the tolerance would
keep the noise measured next. Measured on four
tables, the eigenvalue of a component with no variance came out between
-2e-16 and 3e-16 times λ_1, and at 1000 x 20000 it was 1.3e-14 times λ_1,
where the threshold is 4.4e-12 times λ_1 and the smallest eigenvalue of a
component with variance was 0.4 times λ_1 (**Open 3**, below).

**The divisor of the standard deviation is n**, not n - 1, inherited
from pyNei, `data.std(axis=0)` in `do_pca` of `pynei/pca.py`. R's `prcomp`
divides by n - 1, so its projections of standardized data are pyNei's
times sqrt((n - 1)/n), 0.9975 at 200 individuals. The percentages and the
princomps are the same in both, and so are the projections when the data
is not standardized (**Open 5**, below).

## The PCA of the variants

### What it gives

Each variant becomes one number per individual, its dosage: how many
alleles of the genotype are not the major allele of the variant, 0, 1 or 2
in a diploid. The major allele is the most frequent among the called
alleles of the variant, those of half called genotypes included, and the
lowest numbered of two that are equally frequent. Which allele is the
major one changes the sign of the weight of that variant and nothing else,
because a standardized column with its dosages counted from the other
allele is the same column times -1.

A genotype with any allele missing has no dosage, and takes the mean of
the dosages of its variant, so that after centering it is 0 and pulls its
individual nowhere. The column is then standardized with

    mean = sum of the called dosages / number of called genotypes
    sd   = sqrt( sum over the called genotypes of (dosage - mean)² / n )

where n is all the individuals and not the called ones, since the missing
ones are in the column with a deviation of 0. A variant with much missing
data has a smaller sd for the same frequencies, and its called genotypes
weigh more.

A variant whose called genotypes all have the same dosage has no variance
and is left out: a variant with one allele, and also one where every
individual is heterozygous, whose major allele frequency is 0.5 and which
no filter by frequency would catch. The variants that were used are the
columns of `princomps`.

### What it is in Python and in TypeScript

```python
do_pca_from_variants(
    variants: Variants,
    transform_to_biallelic: bool = False,
    num_prin_comps: int = 10,
) -> PCAResult
```

`PCAResult` is pyNei's frozen dataclass with one field more:
`projections`, a frame with the names of the individuals as index and
the components as columns; `explained_variance_percent`, a series over
the components; `princomps`, a frame with the components as index and,
as columns, the position of each variant that was used among the
variants that the `Variants` gives, from 0; and `pass_stats`, the
`PassStats` of `docs/specs/variant.md`, whose `num_vars` is how many
variants the steps of the `Variants` let through, used or not, and whose
`filtering` has the counts of each filter. The two passes go through the
same steps and count the same, so the stats are those of one. The
components are named as pyNei names them, `PC0`, `PC1`, with zeros on
the left to the width of the number of components, `PC000` for 200 of
them. It mirrors `do_pca_from_variants` of `pynei/pca.py`, a function
over a `Variants` whose filters, in popnei, are steps of it and not
functions around it.

The differences from pyNei:

- `num_prin_comps` is new. pyNei gives the weights of every variant in
  every component, 0.8 GB for 100000 variants of 1000 individuals and
  80 GB for a million of 10000. popnei gives `princomps` for the first
  `num_prin_comps` components, 10 by default, and gets them with a second
  pass over the variants, because a weight needs the eigenvectors, which
  are known when the first pass ends. With 0 there is no second pass and
  `princomps` has no rows, and still has the variants that were used as
  its columns. More than the components that there are gives those there
  are. The owner decided it on 21 September 2026; the options not taken
  were no weights at all, and all of them as pyNei gives them, which
  would hold the whole matrix and end the streaming. Whether the same
  number also cuts the projections is **Open 1**, below.
- The components with no variance are not given, as said above.
- The sign is fixed, as said above.
- `num_threads` is not an argument. The owner decided on 22 September
  2026, in `docs/specs/dists.md`, that no calculation of popnei has it:
  the threads are those of rayon's pool, and those of the backend of the
  linear algebra, as `docs/specs/linalg.md` says.
- A variant with more than two alleles is looked for in each variant and
  not in each chunk, and one with no called genotype is left out and is
  not an error. Both are under "What pyNei does that is odd".

In TypeScript it is `doPcaFromVariants(variants, {transformToBiallelic,
numPrinComps})`. The result has `individuals`, an array of names,
`numComps`, `projections` as a `Float64Array` of individuals x
components, row after row, `explainedVariancePercent` as a
`Float64Array`, `numPrinComps`, `princomps` as a `Float64Array` of
components x used variants, and `usedVars`, a `Uint32Array` with the
positions of the variants that were used, which is `used_cols` of the
core, named there for the variants and the traits alike, and `passStats`.
It has no names of components.

### Errors and the cases pyNei asserts

Each of these is a `ValueError` in Python:

- A variant with more than two different alleles among its called ones,
  when `transform_to_biallelic` is false. The message gives the position
  of the variant among those given and says which argument to pass. With
  the argument true every allele that is not the major one counts the
  same. The alleles are those the genotypes hold and not those the VCF
  lists: a variant with `ALT` of `C,G` and no `G` called has two.
- No variants, and no variant with variance, which is what one individual
  gives. The message for the first is "there are no variants to do a PCA
  with", pyNei's without its "012 matrix" and starting in lower case as
  every message of the core does. The second says the rule, which is about
  the dosage and not the genotype: "no variant has more than one dosage
  among its called genotypes, so none of them varies and there is nothing
  to do a PCA with". pyNei's is "Every variant has the same genotype in
  every sample", and two datasets reach it and contradict it: one whose
  every genotype is missing, where no genotype is called at all, and one
  of three alleles read as biallelic whose genotypes are `0/1`, `0/2` and
  `0/1`, which are three different genotypes of one dosage. Both libraries
  drop a variant by its dosages, `std(mat012) > 0` in
  `_remove_vars_with_no_variance` of `pynei/pca.py` and `dosages_seen < 2`
  in `crates/popnei/src/variant.rs`, so the wording is pyNei's mistake and
  not a difference of behaviour. `docs/specs/kinship.md` gives the same
  sentence to the same rule.
- A negative `num_prin_comps`.
- A source of no individual: the components are the axes the individuals
  of a dataset are placed on, and there is nobody to place, and the
  standardizing of a block would read its rows in chunks of no allele.
  Every source of popnei has one individual at least, as
  `docs/specs/block.md` says, so it is a caller of the function of the
  core crate with a reader of its own that reaches it.
- A dataset beyond what this analysis counts in, which is one of four and
  the message says which. A ploidy above 254: the first pass writes the
  genotype of each individual as one byte, its dosage or the missing
  genotype, as "Speed" below has it, and the 256 dosages of a ploidy of
  255, the largest the VCF reader takes, are with the missing genotype
  one value more than a byte holds. More than 46340 individuals: the
  individuals x individuals matrix would hold more than the 2147483647
  values the linear algebra counts in, which `docs/specs/linalg.md` has.
  More variants than a `usize` counts, 4294967295 in WebAssembly, where
  a `usize` is 32 bits: the variants given, used or not, are counted in
  one, and so are the positions of the used ones, so a pass of more is
  refused instead of counted into a number that wrapped. And weights
  that are more values than a `usize` counts, which is
  `num_prin_comps` times the variants that were used and which fewer
  variants reach the more components are asked for.

`test/test_pca.py` of pyNei asserts, for this function: that the variants
with one allele and the variant where everyone is heterozygous are not
among the columns of `princomps`
(`test_pca_from_variants_drops_the_monomorphic_variants`); the error when
every variant is fixed; that an individual with half its genotypes missing
falls on the side of its own population, the 12 and 8 individuals fixed
for different alleles of `test_pca_vars_with_missing_gts`; and that the
projections do not change with 4 threads. The dosages themselves,
`to_012`, are asserted in `test_mat012` and
`test_mat012_keeps_the_missing_gts`, where `0/-1` is missing. popnei takes
the five, the last as a test that the result does not change with the
size of the blocks, within 1e-10, because the sum is added up in another
order.

### What pyNei does that is odd

The three runs were of pyNei at commit ef0ca6e, on a `Variants` built with
`Variants.from_gt_array`.

`_create_012_gt_matrix` counts the alleles of the whole chunk, the few
thousand variants pyNei works on at once, and not those of each variant.
Two variants of two alleles each, one with 0 and 1 and the other with 0
and 2, raise the error of more than two alleles when they are in one
chunk, and pass with chunks of one variant. Whether a user gets an error
then depends on where the chunks were cut, which popnei could not
reproduce if it wanted to, since its blocks are cut elsewhere. popnei
counts the alleles of each variant.

A variant with no called genotype is an error in pyNei, "There are
variants that only have missing data", while a variant with one allele is
left out in silence. In a loop over rows the two are the same case, a
variant with no variance (**Open 2**, below).

`_remove_vars_with_no_variance` finds the variants to leave out with
`std > 0` on floats. For dosages it is exact, since the mean of equal
small integers is that integer. popnei decides it on the counts of the
called dosages and computes no float for it.

### How it runs

Two passes over the variants. The binding crate opens the reader of each
from the source and the steps of the `Variants`, as it does for
`write_vars`, keeps it, and lends it to the core, as `docs/specs/filters.md`
asks of every calculation, so that it can read the counts of the filters
when the pass ends; the core puts `reblock` before it, since a filter
leaves blocks of uneven size and the product is matrix work. Both passes
count the same, and the stats are those of the first.

The first pass, per block: the standardized rows of the block go into a
buffer of variants x individuals in `f64`, 40 MB for a block of 5 million
genotypes, with rayon across the rows; each row needs the counts of its
dosages and then one more reading of its genotypes, and a row with no
variance is not written. Then, from outside rayon, the product of the
buffer with itself is added to G, the lower half only, through `linalg`.
What is kept from one block to the next is G, individuals x individuals,
and the position of each variant that was used, which is `used_cols` of
the result: one bit per variant would say the same and the result carries
the positions anyway, so the bits would be a second copy of them. When the
variants end, `linalg` gives the eigenvalues and the eigenvectors of G.

The second pass, when `num_prin_comps` is above 0, standardizes each block
again and multiplies it by the first eigenvectors divided by sqrt(λ), a
product of variants x individuals by individuals x `num_prin_comps`. That
product gives the weights variant after variant and the result holds them
component after component, so each block writes its weights into the
columns of its own variants and no copy of the whole matrix is made. It
checks that its variants are those of the first pass, and stops with an
error if they are not, which is what a file that changed between the two
gives; in Python a `RuntimeError`, since no argument is wrong and the core
has no file name to give. Four things it compares: the individuals and the
ploidy of its reader, which a second reader over another dataset differs
in and which would otherwise be read as rows of the size the first pass
had; how many variants it gave; and which of them have variance, each in
the order of the source.

Memory: G and its eigenvectors, 16 MB at 1000 individuals and 1.6 GB at
10000, which no browser tab holds; the buffer of a block; the positions of
the variants that were used, 8 bytes each, 8 MB for a million; and the
princomps, 80 MB for 10 components of a million variants.

### How it is verified

Against R 4.6.1, `prcomp(x, center=TRUE, scale.=TRUE)`, with its
projections multiplied by sqrt(n/(n - 1)). R gets the dosages that pyNei's
`create_012_gt_matrix` gives, with the missing ones replaced by the mean
of their variant and the variants with no variance taken out in R. plink2
cannot give the literals: its `--pca` standardizes each variant by
2p(1 - p), the variance that the allele frequency p predicts, and not by
the variance of the dosages, and it gives eigenvectors of length 1. On the
panel its first two components correlate 0.9999 with pyNei's and the third
0.989, with `meanimpute` 0.991.

`tests/reference/pca/make_reference.py`, run from a checkout of pyNei at
ef0ca6e, writes the VCFs, the dosages and what pyNei gives, and
`reference.R` what R gives; both files say how they are run, and their
outputs are beside them. Every component in those files has the sign of
the rule above. Over all the files pyNei and R differ by 5e-11 at most in
the projections and 5e-13 in the princomps. The literals below are
written with 12 significant digits, of projections up to 13.3, so the
tests compare projections, percentages and princomps with them within
1e-9.

The panel is pyNei's `test/gwas_reference/sim_missing.vars` written as
`sim_missing.vcf`: 200 individuals, 1200 variants of two alleles, three
subpopulations, 7128 genotypes missing whole, every variant with
variance. The literals are in `sim_missing.r.*.tsv` with more digits:

| | PC000 | PC001 | PC002 |
|---|---|---|---|
| projection of `s000` | 1.57303591801 | 12.9003037259 | -5.07968498438 |
| projection of `s001` | 2.95102060115 | 13.3478743206 | -6.89340955267 |
| explained_variance_percent | 7.60577910944 | 5.55551802152 | 1.56605373725 |
| weight of variant 0 | 0.0665938024571 | -0.0296301101543 | |
| weight of variant 1 | 0.0131069975393 | 0.0518944546112 | |

The worked example, `worked.vcf`, which becomes the first cargo test: 5
individuals, 5 variants, ploidy 2.

| variant | genotypes | dosages | standardized |
|---|---|---|---|
| 0 | 0/0 0/1 1/1 0/0 0/1 | 0 1 2 0 1 | -1.069045 0.267261 1.603567 -1.069045 0.267261 |
| 1 | 1/1 1/1 0/1 ./. 1/1 | 0 0 1 none 0 | -0.645497 -0.645497 1.936492 0 -0.645497 |
| 2 | 0/0 0/0 0/0 0/0 0/0 | 0 0 0 0 0 | left out |
| 3 | 0/1 0/1 0/1 0/1 0/1 | 1 1 1 1 1 | left out |
| 4 | 0/2 0/0 2/2 0/. 0/0 | 1 0 2 none 0 | 0.337100 -1.011300 1.685500 0 -1.011300 |

Variant 1 has allele 1 as its major allele, and a mean of 0.25 over its 4
called genotypes and an sd of sqrt(0.75/5) over the 5 individuals.
Variant 4 has two alleles, 0 and 2, and a half called genotype, which is
missing. The result has 3 components and the used variants 0, 1 and 4:

| | PC0 | PC1 | PC2 |
|---|---|---|---|
| i0 | -0.762287762540 | 0.994550847306 | -0.320852228237 |
| i1 | -0.863187556071 | -0.875030288320 | -0.007193635382 |
| i2 | 3.025289193432 | -0.097439191351 | -0.021646302880 |
| i3 | -0.536626318751 | 0.852948920685 | 0.356885801880 |
| i4 | -0.863187556071 | -0.875030288320 | -0.007193635382 |
| explained_variance_percent | 76.7440710447 | 21.7166910409 | 1.53923791438 |
| weight of variant 0 | 0.501967957373 | -0.797860657405 | -0.333836099210 |
| weight of variant 1 | 0.648464626925 | 0.091784778206 | 0.755691194944 |
| weight of variant 4 | 0.572295201271 | 0.595813667059 | -0.563457431176 |

`worked3.vcf` adds a sixth variant, `0/1 1/2 0/0 0/0 2/2`, with three
alleles. Without `transform_to_biallelic` it is the error. With it the
dosages are 1 2 0 0 2, pyNei and R give 4 components of which the last has
a percentage of 1e-31 or less, and popnei gives 3, whose numbers are in
`worked3.r.*.tsv`.

Every fixture above is diploid, and a dosage is how many alleles of a
genotype are not the major one at any ploidy, so one tetraploid case is
checked as well: 4 individuals and 4 variants, with the genotypes below.
Its numbers are numpy 2.5.3's, by the route of "What both analyses
compute" on the dosages of "What it gives", worked out on 22 September
2026; pyNei gives no tetraploid reference, since `do_pca_from_variants`
of a tetraploid `Variants` is not among its tests.

| variant | genotypes | major allele | dosages |
|---|---|---|---|
| 0 | 0/0/0/0 0/0/1/1 0/1/1/1 1/1/1/1 | 1 | 4 2 1 0 |
| 1 | 0/0/0/0 0/0/0/0 0/0/0/1 0/0/1/1 | 0 | 0 0 1 2 |
| 2 | 0/0/1/1 0/0/1/1 0/0/1/1 0/0/1/1 | 0 | left out |
| 3 | 0/0/0/1 ./././. 0/1/1/1 1/1/1/1 | 1 | 3 none 1 0 |

The first variant alone, standardized, is 1.52127765851,
0.169030850946, -0.507092552837 and -1.18321595662: its mean is 1.75 over
its four called genotypes and its standard deviation sqrt(2.1875). The
third variant has one dosage among its called genotypes and is left out,
and the fourth has a genotype with every allele missing. The result has
3 components and the used variants 0, 1 and 3:

| | PC0 | PC1 | PC2 |
|---|---|---|---|
| i0 | 2.30350257852 | -0.45498666453 | -0.0168202039114 |
| i1 | 0.603260051101 | 0.692756790873 | -0.0540239409669 |
| i2 | -0.647570806465 | 0.0576947878489 | 0.143573693133 |
| i3 | -2.25919182316 | -0.295464914192 | -0.0727295482544 |
| explained_variance_percent | 93.2778538477 | 6.47960866887 | 0.242537483385 |
| weight of variant 0 | 0.590326503368 | -0.327593824194 | -0.737697028442 |
| weight of variant 1 | -0.556614390531 | -0.827089115324 | -0.0781281995539 |
| weight of variant 3 | 0.584546866962 | -0.45673392874 | 0.670596062218 |

The literals of the three tables and of `worked3` are checked at
`pca_of_variants` of "The Rust interface", on a reader over blocks that
the test builds, once with `num_prin_comps` 3 and once with 0. Two things
that `pca_of_variants` hides are pinned below it: the standardizing of one
row, which the dosages of `test_mat012` and of the tetraploid variant
above are asserted on, and the size of the blocks, which `reblock` puts
back to the size popnei chooses before the passes see them, so the passes
are given the sizes to compare directly.

Against pyNei, in pytest, at `do_pca_from_variants`: both libraries on
`sim_missing.vcf` and on `worked.vcf`, the first 10 components of the
first and the 3 of the second, within 1e-9, after the test gives pyNei's
components the sign of the rule. pyNei is called with
`transform_to_biallelic=True` on `worked.vcf`, which it refuses otherwise
for the chunk wide count above, and which changes no dosage of a variant
of two alleles. The test under node checks the literals of the worked
example at `doPcaFromVariants`.

## The PCA of a table

### What it gives

The analysis of "What both analyses compute" on a table that the user
brings, individuals x traits, with the two steps on the columns as
options: centering, and standardizing, which puts traits measured in
different units on one scale. Without standardizing, the traits with the
largest numbers dominate. Without centering, the first component mostly
points at the mean of the data. No value may be missing.

### What it is in Python and in TypeScript

```python
do_pca(
    data: pandas.DataFrame,
    center_data: bool = True,
    standardize_data: bool = True,
) -> PCAResult
```

The index of `data` names the rows of `projections` and its columns name
those of `princomps`. It gives every component that has variance, with
its weights: the table is in memory already, so there is no
`num_prin_comps`. It reads no `Variants`, so its result has no
`pass_stats`, and the field is `None` there. It mirrors `do_pca` of `pynei/pca.py`. The differences
from pyNei, besides the sign and the components with no variance:

- pyNei spells the argument `standarize_data` (**Open 4**, below).
- A value that is not finite is refused here. pyNei refuses NaN, and an
  infinite value reaches numpy's SVD, which raises `LinAlgError`.
- A trait has no variance when all its values are equal, compared as they
  are. pyNei tests `std == 0` on floats, which a column of 0.1 repeated
  can pass with an sd of 1e-17, and the division then gives numbers with
  no meaning.
- Fewer than 2 rows is an error. With one row pyNei raises the error of
  the traits with no variance when it standardizes, and without
  standardizing it divides by n - 1 = 0 and gives percentages of NaN.

In TypeScript it is `doPca(data, numRows, numCols, {centerData,
standardizeData})`, with `data` a `Float64Array`, row after row, and the
result that of `doPcaFromVariants` without the three fields that a table
has nothing to fill with: `individuals`, the names of the rows;
`usedVars`, since every trait is used; and `passStats`, since no
`Variants` was read. The names of the rows and of the traits stay with the
application.

### Errors and the cases pyNei asserts

Each is a `ValueError` in Python: a value that is not finite;
`standardize_data` with `center_data` false; fewer than 2 rows or no
columns; and, when standardizing, a trait with no variance, with a message
that gives how many there are and names the first ten, as pyNei's does, so
that the user can take them out. The core gives the positions of those
columns and the Python layer puts the names. Without standardizing such a
trait is no error and gets a weight of 0.

Two more are a `ValueError` as well, and both are of a table that the
analysis cannot be done on in `f64`, whatever the user meant by it.

The first is a trait whose mean or whose standard deviation is not a
number the analysis can use, which happens in three ways. The sum of a
trait can be above the largest `f64`, 1.8e308, and then its mean is an
infinity and every value of it becomes a NaN; that one is found whenever
the table is centered. The squares of the deviations of a trait can sum
above that number, which values of 1e154 do, and then its standard
deviation is an infinity, the trait becomes a column of zeros, and it
would leave the analysis with a weight of 0 and no word, looking like a
trait with no variance. Or those squares can all fall below the smallest
`f64` that is not 0, 5e-324, which values of 1e-200 do, and then the
standard deviation is 0 although the values of the trait differ, and the
division gives infinities. The error names the trait, by its position as
above, and says which of the three happened, so the user can scale that
trait or take it out.

The second is a table in which no trait has variance once it is centered:
every value of every trait equal to the others, or a table of zeros. There
is no direction to give, and popnei says so with the message "no trait has
variance, there is nothing to do a PCA with", the wording of the PCA of
the variants for the same case, which every message of the core starts in
lower case as its own. pyNei gives 0 for every projection and a
percentage of NaN for every component. It cannot happen when the table is
standardized, since a trait with no variance is refused first.

Two are a `RuntimeError`, because no argument of `do_pca` can give them.
One is a buffer that does not hold exactly `num_rows` x `num_cols` values,
which only a caller of the Rust function of "The Rust interface" reaches,
since each binding crate takes the two numbers from the array it was
given. The other is an error of the `linalg` crate, which the core wraps
with the operation it was doing, the product of the table with itself or
the eigendecomposition: the dimensions that crate refuses are checked
here, and so are the values, of the table and of the mean and the standard
deviation of every trait, so what is left for it to refuse is a table
whose values are finite and whose products are not, 1e200 among them,
where it refuses the matrix it is asked to decompose.

`test_pca` asserts the four princomps of iris, up to sign, and
`test_pca_refuses_traits_with_no_variance` the error, that its message
names the trait, and that without standardizing the same table gives 3
projections x 3 components. popnei gives that table 2 components, since
the third has no variance, and the test says so.

### How it runs

On the whole table, which is copied once to be centered and standardized,
in the layout the product needs: the traits as the rows of the copy when
there are more traits than rows, which is the transpose of the table, and
as its columns otherwise.

The product is over the smaller side, as "What both analyses compute"
says. What the analysis holds at once is that copy, `num_rows` x
`num_cols`; the matrix of the smaller side, which the eigendecomposition
turns into its eigenvectors in the same buffer; the projections,
`num_rows` x `num_comps`; and the weights. When there are more traits than
rows the weights are held twice, `num_cols` x `num_comps` as the product
gives them and `num_comps` x `num_cols` as the result holds them, and the
eigenvectors divided by sqrt(λ) are another `num_rows` x `num_comps`. On a
table of 1000 rows x 8000 traits that is 216 MB: 64 MB for the copy, 8 MB
for the 1000 x 1000 matrix, 8 MB for the projections, 8 MB for the scaled
eigenvectors and 64 MB for each of the two copies of the weights.

Nothing is kept and there are no blocks.

### How it is verified

Against R's `prcomp` on iris, 150 rows x 4 traits, the table of pyNei's
`test/datasets.py`, which `make_reference.py` writes as `iris.tsv`,
standardized and not. The literals, from `iris.r.*.tsv` and
`iris_not_standardized.r.*.tsv`:

| | PC0 | PC1 | PC2 | PC3 |
|---|---|---|---|---|
| standardized, projection of row 0 | -2.26470280881 | 0.480026596521 | -0.127706022300 | -0.0241682038555 |
| explained_variance_percent | 72.9624454133 | 22.8507617867 | 3.66892188928 | 0.517870910715 |
| not standardized, projection of row 0 | -2.68412562597 | 0.319397246585 | -0.0279148275894 | -0.00226243707132 |
| explained_variance_percent | 92.4618723202 | 5.30664831171 | 1.71026098079 | 0.521218387328 |

They are checked at `pca` of "The Rust interface" within 1e-9, with the
princomps of the same files. Against pyNei, in pytest at `do_pca`, the
same two runs, all 4 components, within 1e-9 after the sign. No run
without centering has a reference outside the project: `prcomp` with
`center=FALSE` divides by n - 1 a sum of squares that was not centered, as
pyNei does, so it can be added to `reference.R` when the test is written,
and until then that case is compared with pyNei alone.

Iris has 150 rows and 4 traits, so its product is the 4 x 4 Z'Z, and a
mistake in the other half of "What both analyses compute", the n x n ZZ'
of a table with fewer rows than traits, would not show in it. Two small
tables are checked there as well.

Their numbers were got with numpy 2.5.3 on 22 September 2026, by both
routes of "What both analyses compute" over the same Z, the table
centered and divided by the standard deviation with n in it. One is
`numpy.linalg.eigh` of Z Z', with the projections u sqrt(λ), the weights
Z' u / sqrt(λ) and the percentages 100 λ over the sum of every λ, which is
what popnei computes; the other is `numpy.linalg.svd` of Z, which is what
pyNei computes, with the projections Z v and the percentages over the
squares of the singular values. The sign rule of "What both analyses
compute" was applied to both, and the components under the threshold of
that part were dropped. Over the three runs below the two routes differ by
3e-15 at most in the projections and 1.5e-14 in the percentages, so the
tests compare with the numbers below within 1e-9.

The first table is 3 rows x 5 traits, the rows being 1 2 3 4 5, then
2 4 1 3 2, then 5 1 4 2 6. Centered and standardized it has 2 components
and not 3, because centering takes one of the three dimensions of its rows
out:

| | PC0 | PC1 |
|---|---|---|
| projection of row 0 | -0.347154191646 | 1.62410767053 |
| projection of row 1 | -2.14982434306 | -0.994055030603 |
| projection of row 2 | 2.49697853471 | -0.630052639924 |
| explained_variance_percent | 73.1810836106 | 26.8189163894 |
| weight of trait 0 | 0.420101972196 | -0.513968704337 |
| weight of trait 1 | -0.496432942021 | -0.270671721206 |
| weight of trait 2 | 0.496432942021 | 0.270671721206 |
| weight of trait 3 | -0.317325807736 | 0.686274627813 |
| weight of trait 4 | 0.47950738561 | 0.344001373342 |

Neither centered nor standardized, the same table has all 3 components,
which is the run no program outside the project gives a number for:

| | PC0 | PC1 | PC2 |
|---|---|---|---|
| projection of row 0 | 7.07789159873 | 1.12190772301 | 1.90912901021 |
| projection of row 1 | 4.87213205354 | 2.87255560109 | -1.41801042717 |
| projection of row 2 | 8.66041659868 | -2.53292797372 | -0.762535387595 |
| explained_variance_percent | 87.0392022767 | 9.31343669027 | 3.64736103304 |
| weight of trait 0 | 0.403960199411 | -0.364035502371 | -0.759913160568 |
| weight of trait 1 | 0.28423522248 | 0.703323259808 | -0.419484427674 |
| weight of trait 2 | 0.408147561391 | -0.244470602227 | 0.201897964367 |
| weight of trait 3 | 0.404797068091 | 0.504800545607 | 0.297806276491 |
| weight of trait 4 | 0.652365999566 | -0.241298733997 | 0.342218405411 |

The second table is the one of `test_pca_refuses_traits_with_no_variance`
of pyNei, 3 rows x 3 traits, whose traits `a`, `fixed` and `b` are
1 2 3, 5 5 5 and 3 1 2. Centered and not standardized it has 2 components,
and the trait with no variance has a weight of 0 in both:

| | PC0 | PC1 |
|---|---|---|
| projection of row 0 | 1.41421356237 | 0 |
| projection of row 1 | -0.707106781187 | 0.707106781187 |
| projection of row 2 | -0.707106781187 | -0.707106781187 |
| explained_variance_percent | 75 | 25 |
| weight of `a` | -0.707106781187 | -0.707106781187 |
| weight of `fixed` | 0 | 0 |
| weight of `b` | 0.707106781187 | -0.707106781187 |

The two largest projections of the second component of that table are the
same number with opposite signs, 0.707106781187, and how far apart the
two come out depends on the route to them: popnei gives them equal bit for
bit natively, on Accelerate's LAPACK, one bit apart in WebAssembly, on
faer, and numpy's eigendecomposition of the same table gives them two bits
apart. All three are inside the tolerance of "The sign", so the first of
the two is the positive one on every backend, and the numbers above are
what each of them gives. Before that tolerance the component came out with
one sign natively and the other under faer, and the spec gave it up to its
sign.

Two more tables are checked in the core alone, for what the two above do
not reach. Their numbers come from the same two routes of numpy 2.5.3 and
agree within 1e-15 in the projections.

One is 5 rows x 3 traits, whose traits are 1 2 3 4 7, then 5 5 5 5 5, then
3 1 2 9 2: it has more rows than traits, as iris has, and unlike iris it
loses a component, the second trait having no variance. Centered and not
standardized it has 2 components of the 3 traits:

| | PC0 | PC1 |
|---|---|---|
| projection of row 0 | -0.765373820993 | -2.30958933885 |
| projection of row 1 | -2.58720936972 | -1.01308818829 |
| projection of row 2 | -1.44494156963 | -0.179287089168 |
| projection of row 3 | 5.62553292806 | -0.270886092928 |
| projection of row 4 | -0.828008167714 | 3.77285070924 |
| explained_variance_percent | 66.8261599339 | 33.1738400661 |
| weight of trait 0 | 0.15423335048 | 0.988034449602 |
| weight of trait 1 | 0 | 0 |
| weight of trait 2 | 0.988034449602 | -0.15423335048 |

The other is 3 rows x 3 traits whose values are large enough that the
arithmetic around the eigenvalues has to be written in the order that does
not overflow: 1 2 3, then 2 4 1, then 5 1 4, each value multiplied by
1e153 and, in a second run, by 2.5e153. Centered and not standardized it
has 2 components at both scales, with the percentages 78.8675134595 and
21.1324865405, which are the percentages of the table unscaled. The
projections are of the size of the values, 1e153, and the tests compare
the percentages and the count of the components, which are what the two
orders of the arithmetic differ in. At 1e153 the largest eigenvalue is
1.4e307, so 100 times it is an infinity while 100 times its share of the
total is not; at 2.5e153 it is 8.9e307, so that eigenvalue times the
larger side of the table, 3, is an infinity as well, and a threshold of
infinity leaves no component at all.

## The principal coordinates of distances

Written on 27 September 2026, for stage 4 of popnei_web, the web
application of popnei, which draws the individuals of a dataset on the
first three components of this analysis and of the PCA of the variants
with one piece of code.

### What it gives

A principal coordinate analysis, PCoA, places the individuals on
components as a PCA does, from the distance of every pair of them instead
of a table of their values. The individuals are put in a space where the
straight line between each two of them is as long as their distance, and
the components are the directions of that space along which they vary
most. popnei gives two: the PCoA of a `Distances` that the user brings,
and that of the Kosman distances of the individuals of a `Variants`, which
`docs/specs/dists.md` computes. The second counts every allele of a
variant as itself and takes any ploidy, where the PCA of the variants
reduces a genotype to its dosage.

The method is Gower's (1966, Biometrika 53: 325), which pyNei's code
took from an older Python library, and which two functions of R compute:
`cmdscale` of base R, and `pcoa` of the package ape, which gives as well
each eigenvalue as a share of their sum. With d_ij the distance of the
individuals i and j, and n the individuals:

    A_ij = -d_ij² / 2
    B_ij = A_ij - (mean of row i of A) - (mean of column j of A) + (mean of A)

B, n x n, takes the part that G of "What both analyses compute" has in
the PCA. With λ_j its eigenvalues from the largest and u_j its
eigenvectors of length 1:

    projections of component j        = u_j * sqrt(λ_j)
    explained_variance_percent of j   = 100 * λ_j / sum of every λ of B

The sum of every eigenvalue is the sum of the diagonal of B, which is the
sum of d² over the pairs divided by n: the total variance of the
individuals placed at their distances. The percentages add up to 100.
When the distances are the straight line distances between the rows of a
centered table, B is that table times its transpose, and the projections
are those of the PCA of the table not standardized. The sign of "The sign"
holds here as it does there. An eigenvalue is positive when it is above

    n x 2.2e-16 x (sum of |λ| of B)

negative when it is below minus that, and 0 in between, and the
components are those of the positive ones. B is centered twice: after the
first centering its rows and columns are centered again, which takes out
what the rounding of the first left in their means.

The threshold is not the PCA's, λ_1 x n x 2.2e-16, because the review of
the code on 27 September 2026 found that one too narrow for B. Each cell
of B carries the rounding of its centering, and the eigenvalue 0 that the
centering leaves, and the one that Lingoes' correction brings to 0, came
out beyond λ_1 x n x 2.2e-16 on either side: of 20 Euclidean matrices of
500 random points in 3000 dimensions, 4 were refused as not Euclidean and
7 got a 500th component of rounding, and 7 of 40 random matrices of 100
individuals, corrected by `correct_dists_by_lingoes`, were refused by
`do_pcoa`, whose message then named the correction the user had just
applied. The sum of |λ| is λ_1 at least and n x λ_1 at most, so the
threshold is at most n² x 2.2e-16 x λ_1, 2.2e-8 of λ_1 at 10000
individuals: a component smaller than that is not given. With the two centerings and
this threshold, the largest eigenvalue that is 0 in exact arithmetic came
out at 0.20 of the threshold on LAPACK and 0.095 on faer, over random
points of 50 to 500 individuals, corrected random distances of 5 to 1000,
the regular simplex of 5 to 3000, `four_alleles` and the corrected twin,
and none of those cases was refused or given a component too many.

**A matrix that is not Euclidean is refused.** A matrix of distances is
Euclidean when some space has points whose straight line distances are
those distances, and then every eigenvalue of B is 0 or above. Otherwise
some are negative, and a negative λ is a direction whose projections would
be the square root of a negative number, which no space has. The Kosman
distance is not known to be Euclidean, and the distance of each pair is
taken over the variants both individuals were called at, which is another
set of variants for each pair. On the panel of "How it is verified", 200
individuals and 1200 variants with 3 in 100 genotypes missing, B has 155
positive eigenvalues, 44 negative ones and one of 0, the one that the
centering takes out. On three datasets of `docs/specs/dists.md` with 5 in
100 genotypes missing, of 40 diploid, 12 tetraploid and 12 haploid
individuals, it has none.

A PCoA of a matrix with a negative eigenvalue is an error, unless it is
asked to correct the distances first, and its message names the
correction. How large the negative part is, is said by

    negative_eigenvalues_percent = 100 * (sum of |λ| of the negative ones) / sum of every λ of B

2.98 on the panel. The owner decided this on 27 September 2026. The
options not taken were to give the components of the positive eigenvalues
and leave the negative ones out, which is what ape's `pcoa` does with no
correction, with the percentages adding up to more than 100, 102.98 on the
panel; to give the negative ones as pyNei does, under "What pyNei does that
is odd"; and the correction of Cailliez (1983), below.

**Lingoes' correction.** Lingoes (1971) makes a matrix Euclidean by adding
one constant to every squared distance of two different individuals:

    c = |λ_min|, the most negative eigenvalue of B, and 0 when none is
        below minus the threshold above
    d'_ij = sqrt(d_ij² + 2c)   for i ≠ j

It is ape's `pcoa(d, correction = "lingoes")`. The B of the corrected
distances is B + cJ, where J is the centering matrix, the identity less
1/n in every cell, which takes the mean out of a vector. So it has the
eigenvectors of B, and every eigenvalue but the 0 of the centering is c
larger: the most negative becomes 0 and the others are above it. A PCoA
of the corrected distances has up to n - 2 components, since two
eigenvalues are 0, that of the centering and the most negative one, and its percentages add
up to 100 over the variance of the corrected distances. Cailliez's
correction, which adds a constant to every distance and not to its square,
is ape's other one; its constant is the largest eigenvalue of a matrix of
2n x 2n that is not symmetric, which the `linalg` crate cannot decompose,
and the owner chose Lingoes' alone.

What the correction does to the picture is to push every pair apart by
the same amount of squared distance, which moves the nearest pairs the
most. On the panel c is 0.0142, so 2c is 0.0284, which is 30 in 100 of the
mean of d² over the pairs, 0.0947: the nearest pair goes from 0.137 to
0.217 and the farthest from 0.358 to 0.396. The panel then has 198
components where it had 155 positive eigenvalues, and its first component
holds 9.62 percent where it held 12.35 of the uncorrected total. That is
why the result carries c and the `negative_eigenvalues_percent` of the
distances before the correction, which is what a warning to the user is
made from; which values a user should be warned at is the application's
to decide.

### What it is in Python and in TypeScript

```python
correct_dists_by_lingoes(dists: Distances) -> LingoesCorrection

do_pcoa(dists: Distances) -> PCoAResult

do_pcoa_from_variants(
    variants: Variants,
    min_num_snps: int | None = None,
    correct_by_lingoes: bool = False,
) -> PCoAResult
```

`correct_dists_by_lingoes` is the correction on its own. Its result, a
frozen dataclass, has `dists`, a `Distances` of the corrected distances
with the `names` and the `pass_stats` of the one given; `constant`, c;
and `negative_eigenvalues_percent`, of the distances given. A Euclidean
matrix gives a constant of 0, a `negative_eigenvalues_percent` of 0 and
the same distances. It is what the error
of a PCoA of a matrix that is not Euclidean points to, and it serves a
user who wants the corrected distances for something else, a tree among
them. `do_pcoa` of its distances is the PCoA of
`do_pcoa_from_variants` with `correct_by_lingoes`, within the rounding of
the square roots, when the distances are the Kosman ones.

`do_pcoa` corrects nothing: it refuses a matrix that is not Euclidean
and names `correct_dists_by_lingoes`. `calc_pairwise_kosman_dists` knows
nothing of the correction. `do_pcoa_from_variants` applies it inside,
when `correct_by_lingoes` is true, from the eigenvalues of B, with no
second eigendecomposition and no corrected vector built; when it is
false, the default, a matrix that is not Euclidean is refused. The owner
decided this on 27 September 2026; the option not taken was the argument
on `do_pcoa` as well. That the argument is false by default follows their
word that in Python the user decides.

`do_pcoa_from_variants` is the Kosman distances of
`calc_pairwise_kosman_dists` followed by `do_pcoa`, in one pass over the
variants and with no `Distances` built between the two. `min_num_snps` is
passed on as that function takes it: a pair needs that many variants
called in both individuals to get a distance, and `None` or 0 means one.

`PCoAResult` is pyNei's frozen dataclass with three fields more:

- `projections`, a frame with the names as index, those of the
  individuals of the `Variants` or the `names` of the `Distances`, and
  the components as columns, named `PC0`, `PC1` with zeros on the left as
  the PCA names them;
- `explained_variance_percent`, a series over the components;
- `lingoes_constant`, c, new: 0 when the correction was not asked for or
  the matrix was Euclidean;
- `negative_eigenvalues_percent`, new, of the distances before the
  correction: 0 when it was not asked for, since a matrix with a negative
  eigenvalue is then refused;
- `pass_stats`, new: the `PassStats` of the pass over the variants, or
  those of the `Distances` given to `do_pcoa`, which are `None` for a
  `Distances` the user built.

It mirrors `do_pcoa` and `do_pcoa_from_variants` of `pynei/pca.py`;
pyNei has no correction. The differences from pyNei:

- A matrix that is not Euclidean is refused, or corrected when that is
  asked for; pyNei gives n components, under "What pyNei does that is
  odd". The eigenvalue 0 of the centering gives no component, as the PCA
  leaves out its components with no variance (**Open 3**, above).
- The sign is fixed, as in the PCA.
- `correct_by_lingoes` of `do_pcoa_from_variants`, `lingoes_constant`,
  `negative_eigenvalues_percent` and `pass_stats` are new; the two
  numbers are 0 in the result of `do_pcoa`, which corrects nothing and
  refuses a negative eigenvalue.
- `use_approx_embedding_algorithm` and `num_threads` are not arguments,
  as `docs/specs/dists.md` decided for the Kosman distances.
- A negative or an infinite distance is refused. pyNei squares a negative
  one, which gives it the place of its absolute value and says nothing.
  The Kosman distances are never negative. F_ST and f_2 of
  `calc_pop_dists` are, for two populations that the dataset cannot tell
  apart, as `docs/specs/dists.md` says, which does not clamp them; what to
  put in the place of such a value, 0 or another measure, is the user's
  to choose, and Lingoes' correction is no answer to it, since it squares
  the value first as pyNei does.
- Fewer than two individuals, and distances that are all 0, are refused.
  pyNei gives a percentage of NaN for both.
- The pair with no distance is refused with a message a user can act on,
  below; pyNei's is "dists array has nan values".

In TypeScript they are `correctDistsByLingoes(distances)`,
`doPcoa(distances)` and `doPcoaFromVariants(variants,
{minNumSnps, correctByLingoes})`. The result of `doPcoaFromVariants` has
the fields that the result of `doPcaFromVariants` has for drawing the
individuals, with the same names and shapes: `individuals`, an array of
names in the order of the source after its steps, `filterIndividuals`
among them; `numComps`; `projections`, a `Float64Array` of individuals x
`numComps`, row after row; `explainedVariancePercent`, a `Float64Array` of
`numComps`; and `passStats`. It adds `lingoesConstant` and
`negativeEigenvaluesPercent`, two numbers, and it has no `numPrinComps`,
`princomps` or `usedVars`, since a PCoA has no weights of variants. The
result of `doPcoa` is the same with `names`, those of the `Distances`,
where the other has `individuals`, and the `passStats` of the
`Distances`. `correctDistsByLingoes` gives `{distances, constant,
negativeEigenvaluesPercent}`, `distances` being a `Distances`.
`numPassesOf("doPcoaFromVariants")` is 1, so the function is a variant of
`Consumer` of `crates/popnei-js/src/source.rs` and a name of
`ConsumerName` of `js/popnei/src/passes.ts`, with the tests that hold the
two lists together.

### Errors and the cases pyNei asserts

Each is a `ValueError` in Python and an `Error` with the same message in
TypeScript, where each is in the `@throws` of the functions it can come
from. The core writes the Python names of the functions and the arguments
in its messages, and the TypeScript binding rewrites them in camelCase, as
it does for `transform_to_biallelic`.

- A matrix that is not Euclidean, from `do_pcoa` or
  `do_pcoa_from_variants` without `correct_by_lingoes`, checked after the
  pairs with no distance. The message says
  how many of the n eigenvalues are negative and their
  `negative_eigenvalues_percent`, that a PCoA of such distances would
  draw directions that no space has, and that `correct_by_lingoes` from
  `do_pcoa_from_variants`, or `correct_dists_by_lingoes` from `do_pcoa`,
  makes them Euclidean by
  adding the same amount to every squared distance. On the panel it is 44
  of the 200 eigenvalues and 2.98 percent.
- A pair of individuals with no distance, a NaN. `do_pcoa_from_variants`
  gives one to a pair whose two individuals were called together at fewer
  variants than `min_num_snps`, or at none. The message says how many of
  the pairs have no distance, names the first of them in the order of the
  distance vector, and names the individual that is in the most of them
  with how many, the first of them in the order of the individuals when
  two are in as many. From `do_pcoa_from_variants` it says the user can
  take that individual out with `filter_individuals`, lower
  `min_num_snps`, or run the PCA of the variants, which gives every
  individual a projection; when `min_num_snps` is 0 or 1 those pairs were
  called together at no variant, and the message says so and does not
  name `min_num_snps`, which lowering could not help; from `do_pcoa` and `correct_dists_by_lingoes`,
  whose `Distances` may be of populations, it says that the pair has to be
  given a distance or one of the two taken out of the `Distances`. On the
  panel with `min_num_snps` 1105 it is 35 of the 19900 pairs, the first
  being `s001` and `s082`, and `s082` is in 17 of them. It is checked
  before the eigendecomposition, so it comes whether or not the matrix is
  Euclidean. The core gives the counts and the positions, and each binding
  puts the names.
- A negative or an infinite distance given to `do_pcoa` or
  `correct_dists_by_lingoes`, with the pair, the value, and that a
  negative F_ST or f_2 is of two populations the dataset cannot tell
  apart.
- A constant of Lingoes' correction beyond an `f64`, from
  `correct_dists_by_lingoes`: c is of the size of the squared distances,
  and it is refused when it is above the largest `f64`, 1.8e308, or below
  the smallest normal one, 2.2e-308, where it would keep only a few
  significant digits or read as no correction; distances of about 1.3e154
  and above reach the first, and distances all below about 1e-154 the
  second. The message gives the largest distance, written with an
  exponent, and says to divide the distances by a number near it first.
- Fewer than two individuals, from the three functions: with one there
  is no distance to place. `do_pcoa_from_variants` refuses it before the
  pass, after the two limits below.
- Every distance 0, from the three functions: "every distance is 0, so the
  individuals are all at one point and there is nothing to do a PCoA
  with". It is checked on the distances, before the eigendecomposition,
  where B would be 0 and the threshold 0.
- More individuals than a page holds the analysis of, in TypeScript
  alone, which "How it runs" gives: more than 9381, for the three
  functions, counted on the individuals left after the steps of the
  `Variants`, which are those of the matrix. It is checked before the pass, with a message of the form of
  the PCA's: the principal coordinates of that many individuals hold
  about that many GB, which a page does not hold, and data this large is
  analysed by a program outside the browser, popnei in Python among them.
  More than 46340 individuals everywhere, which is the PCA's limit for the
  linear algebra.
- The errors of the Kosman pass, of the class and with the message
  `calc_pairwise_kosman_dists` gives them: a pass that gave no variant,
  sums of a pair beyond a `u32`, a source that could not be read, and a
  `min_num_snps` that is not a whole number from 0 to 4294967295, which is
  a `TypeError` in Python when it is not an integer at all.

An error of the `linalg` crate is a `RuntimeError`, as in the PCA, and so
is a vector given to a function of the core whose length is not
n(n - 1)/2, which only a caller of the Rust functions reaches, since a
`Distances` of Python or of TypeScript refuses such a vector when it is
built. So is the vector of ones found where it cannot be: a component
given whose eigenvector is not at a right angle to it, which would put
every individual at one projection, and, in the correction inside the
analysis, a band of eigenvalues 0 that does not hold it, where the
projection of "How it runs" would put it in the place of a real
eigenvector. Neither has been seen, since the eigenvalue 0 of the
centering came out at 0.20 of the threshold at most; they are checked so
that rounding beyond that would be an error and not a wrong analysis.

`test/test_pca.py` of pyNei asserts, of `do_pcoa`, that on the ten
distances of `test_pcoa` the individual `i1` is nearer `i2` than `i4`
along `PC0`, and runs `do_pcoa_from_variants` once each way of
`use_approx_embedding_algorithm` and asserts nothing of its numbers.
Those ten distances are not Euclidean, so popnei refuses them without the
correction and checks them with it at the numbers of R below.

### What pyNei does that is odd

pyNei gives all n components. Those of the negative eigenvalues have
projections u_j times the square root of |λ_j| and a negative percentage,
and the one of the eigenvalue 0 has projections that are the square root
of an eigenvalue that is rounding: 4.7e-9 on the ten distances of the
worked example, whose eigenvalue 0 comes out at -1.1e-16, which gives a second negative
percentage there. On the panel it
gives 200 components, 44 of them with a negative percentage when it
starts from the variants and 45 when it starts from the distances of R,
where the eigenvalue 0 rounds below 0, and a user who plots one of those
draws a direction that no space has. Its percentages are over the sum of
every eigenvalue, the negative ones included, as ape's `Relative_eig`, and
on a Euclidean matrix they are popnei's.

### How it runs

`do_pcoa_from_variants` makes the one pass of the Kosman distances,
`calc_kosman_sums` of `docs/specs/dists.md`, over a reader that the
binding crate opens from the source and the steps of the `Variants`. B is
asked of the machine before the pass, so that the sums, asked for after
it, lie above it in the memory of wasm and their room is free again at the
top when they are given back, and so that a B the machine cannot hold is
refused before the variants are read. B is then built from the sums of
each pair and not from a vector of distances, which is never made, and the
sums are dropped before the eigendecomposition of B by `linalg`. The
reader and its chain of filters stay alive until the analysis returns,
since the binding reads the counts of the filters from them, and a VCF
reader holds buffers of about 16 MiB after its pass; the measurement of
work package 3 of `docs/plans/pcoa.md` is of a VCF opened in the page, so
the limit takes them in. `do_pcoa` builds B from the vector it was given, which the core takes
over and drops before the eigendecomposition in the same way; kept to the
end, the vector would add 4 bytes a cell to the peak below.
`correct_dists_by_lingoes` keeps its vector through the
eigendecomposition, since the corrected vector is written from it, which
is 4 bytes a cell where `do_pcoa` writes 8 of projections, so its peak is
below that of the PCoA.

B is built from the distances divided by the largest of them, and the
projections are multiplied back by it, so that distances above 1.3e154,
whose squares are beyond an `f64`, and distances all below 1.5e-162,
whose B would give no component, are analysed as any other; found and
tested at 1e200 and 1e-200 when the code was written, on 27 September
2026.

In `do_pcoa_from_variants` with `correct_by_lingoes`, the eigenvalues
and the eigenvectors of the corrected B come from those of B, as "Lingoes' correction" says: each
eigenvector keeps its direction and its eigenvalue grows by c, except the
vector of ones, the eigenvector of the 0 of the centering, which stays at
0. A decomposition tells the vector of ones apart from the other
eigenvectors only when the eigenvalue 0 has no other: two individuals at
distance 0 from each other and at the same distance from every other
individual, two clones, give a second eigenvector of 0, and the
decomposition may give any two directions of the plane of the two. So the
eigenvectors whose eigenvalues are 0 within the threshold are projected
orthogonal to the vector of ones and made orthonormal again, which leaves
one fewer, and those take the eigenvalue c. The correction puts the two
clones at a distance of sqrt(2c) from each other, along a component of
their own on which they are at sqrt(2c)/2 either side of 0.

`correct_dists_by_lingoes` needs c alone, the most negative eigenvalue,
and gives the vector sqrt(d² + 2c), pair by pair.

The memory, counted in bytes for each of the n² cells of an n x n matrix:
the sums of the pass, two `u32` for each of the n(n - 1)/2 pairs, 4 per
cell; B, 8. Then B is decomposed, with its eigenvectors and the workspace
of faer, which the PCA measured under node at 6.1 times the matrix, 48.8
per cell; and while the eigenvectors are there the projections are
written, up to n - 1 components of 8 bytes for each individual, 8 per cell
when every eigenvalue but the one of the centering is positive, which the
correction makes the common case. The sums are dropped before B is
decomposed, so the peak is the second of the two steps, at most 56.8
bytes per cell, 5.7 GB at 10000 individuals natively, where LAPACK's
workspace is smaller than faer's.

In wasm the peak was measured on 27 September 2026, under node 26.8.2 on
the owner's Apple M5 Pro, by `js/popnei/bench/memory_of_pcoa.mjs`, which
opens the bytes of a VCF of n diploid individuals and 300 variants with 2
in 100 genotypes missing and runs `doPcoaFromVariants` with
`correctByLingoes`, nearly every eigenvalue then positive. The memory of
wasm grew, after the VCF was opened, by 45.6 bytes a cell at 3000
individuals and 44.4 at 8695 to 9413, below the 56.8 counted above: the
projections are written once faer's workspace is given back. The analysis
ran at 9413 individuals and ended the module at 9414, 1.4 s in, with one
allocation that did not fit and not with the memory full, which is where
the PCA's edge is too, 9410 ran and 9415 did not. So the three functions
take the PCA's limit, 9381 individuals, 33 below the smallest number that
trapped, where the PCA's is 34 below its own. `doPcoa`, which copies the
vector of distances into wasm besides, was run at 9381 individuals by the
review of this work, the same day and on the same machine, and fit: 44.25
bytes a cell, and 9414 trapped. `correctDistsByLingoes` holds less, since it
writes no projections, and takes the same limit. The limit counts 48.8
bytes a cell, the PCA's 6.1 times the matrix, above the 44.4 measured,
because the edge is set by one allocation that does not fit and not by
the memory the analysis holds.

### How it is verified

Against R 4.6.1 with ape 5.8.1: `pcoa(d, correction = "none")`, which on
a Euclidean matrix gives the components of the positive eigenvalues and,
as `Relative_eig`, each eigenvalue over the sum of all of them; and
`pcoa(d, correction = "lingoes")`, which gives the constant, the
components of the corrected matrix as `vectors.cor` and their share as
`Rel_corr_eig`. `tests/reference/pca/pcoa_reference.R`, run as its header
says, writes those of five matrices with the sign rule into
`*.pcoa.r.*.tsv` and `*.lingoes.r.*.tsv` beside it, with 15 significant
digits. The literals are compared within 1e-9, as in the PCA. ape counts
an eigenvalue as positive above 1.5e-8, the square root of the machine
epsilon, and popnei above n x 2.2e-16 x (sum of |λ| of B), so the two
give the same components on these five matrices, whose smallest positive eigenvalue is
8.2e-5 on the panel, and could differ on distances of a much smaller
size.

The worked example is the ten distances of pyNei's `test_pcoa`, of five
individuals `i1` to `i5`, in the order of the distance vector: 0.2, 0.3,
0.9, 0.9, 0.1, 0.8, 0.7, 0.7, 0.8, 0.2. The eigenvalues of B are
0.759739804991026, 0.0891457990391564, 0.0271213359309426, 0 and
-0.0640069399611263, whose sum is 0.812, the sum of the ten squares,
4.06, over 5; they are not in the result and are given to check B by
hand. They are not Euclidean, so `pcoa` refuses them, with 1 negative
eigenvalue of 5 and a `negative_eigenvalues_percent` of 7.88262807403034.
`correct_dists_by_lingoes` gives the constant 0.0640069399611263, that
`negative_eigenvalues_percent`, and, for the first pair, sqrt(0.2² + 2c)
= 0.409894962060102. `pcoa` of the corrected vector has 3 components:

| | PC0 | PC1 | PC2 |
|---|---|---|---|
| i1 | -0.431869046368213 | -0.0415086068851603 | 0.224881556185394 |
| i2 | -0.283479006142767 | -0.158680403878776 | -0.138801060767122 |
| i3 | -0.269028184151739 | 0.212468396790800 | -0.131478395434635 |
| i4 | 0.492920681079785 | 0.195146288470815 | 0.0612591236372794 |
| i5 | 0.491455555582935 | -0.207425674497679 | -0.0158612236209160 |
| explained_variance_percent | 77.1278402980914 | 14.3397713765961 | 8.53238832531248 |

with a `lingoes_constant` and a `negative_eigenvalues_percent` of 0,
since `pcoa` made no correction. These are the first
cargo tests, at `pcoa` and `correct_dists_by_lingoes` of "The Rust
interface", and the tests of `doPcoa` and `correctDistsByLingoes` under
node.

The same distances with a sixth individual, `i6`, at distance 0 from `i5`
and at the distances of `i5` from the others, are the clones of "How it
runs": their eigenvalue 0 has two eigenvectors.
`correct_dists_by_lingoes` gives them a constant of 0.072805704185076,
and `pcoa` of the corrected vector 4 components, of percentages 74.4550063085174,
12.9405131851367, 7.29289082221065 and 5.31158968413514, and on `PC3` the projections 0.190795314650381 for
`i5` and -0.190795314650381 for `i6`, which is sqrt(2c)/2 each side and
the first of the two positive by "The sign"; the rest is in
`small_twin.lingoes.r.projections.tsv`. Checked at
`correct_dists_by_lingoes` and `pcoa`. The rule of "How it runs" for the
eigenvectors of the eigenvalue 0 is not reached there, since the corrected
vector already sets the two apart; the panel with a clone, below, reaches
it.

The panel is `tests/reference/dists/panel.vcf.gz` of `docs/specs/dists.md`,
whose Kosman distances the tests of that spec check against R's
`gd.kosman`, and R gets those distances from `panel.gdkosman.tsv`.
Without the correction it is refused, 44 negative eigenvalues of 200 and
2.98343616373556 percent. With it, it has 198 components, a
`lingoes_constant` of 0.014182298472042 and a
`negative_eigenvalues_percent` of 2.98343616373556:

| | PC0 | PC1 | PC2 |
|---|---|---|---|
| projection of `s000` | 0.0131009923566961 | 0.103593034190948 | -0.0461610570616341 |
| projection of `s001` | 0.0188069490905941 | 0.102449570787996 | -0.0506243303501242 |
| projection of `s199` | -0.0728375733863349 | -0.0163081366424616 | 0.0108102147666687 |
| explained_variance_percent | 9.62407114041929 | 6.65101405201371 | 1.73580899943624 |

They are checked at `pcoa_of_variants` on a reader over the VCF, with the
whole of `panel.lingoes.r.projections.tsv`, and at `doPcoaFromVariants`
under node, with the refusal without the correction. The refusal of the
panel with `min_num_snps` 1105 is checked at both, with its counts and its
names.

The panel with a clone is `panel_clone.vcf.gz`, which
`tests/reference/pca/make_panel_clone.py` writes: the panel with a 201st
individual, `s200`, whose genotypes are those of `s000`, so that its
Kosman distance to `s000` is 0 and to the others that of `s000`; popnei's
Kosman distances of it are those R was given, exactly. It is the clones of
"How it runs" in `do_pcoa_from_variants`, and it is the one check of the
rule for the eigenvectors of the eigenvalue 0. With `correct_by_lingoes`
it has 199 components, a constant of 0.0143697121693818, a
`negative_eigenvalues_percent` of 2.97950453177523, and the first three
percentages 9.54980310926884, 6.68888731987025 and 1.77807340091020. On
`PC155`, the component that sets the two apart, `s000` is at
0.0847635303930356 and `s200` at -0.0847635303930342, sqrt(2c)/2 each
side; on every other component the two are within 6e-14 of each other.
Checked at `pcoa_of_variants` against the whole of
`panel_clone.lingoes.r.*.tsv`.

A Euclidean matrix, the Kosman distances of `four_alleles.vcf.gz`, 40
individuals, is not refused: 39 components, percentages that add up to 100
within 1e-12, the first of them 4.33989428824864, and a
`negative_eigenvalues_percent` of 0; with `correct_by_lingoes` the
constant is 0 and the result the same. Checked at `pcoa_of_variants`, against
`four_alleles.pcoa.r.*.tsv`.

Against pyNei, in pytest, which has no correction and gives its numbers
for any matrix: `do_pcoa` of both libraries on the Kosman distances of
`four_alleles.vcf.gz`, the 39 components of popnei against the first 39 of
pyNei; and `do_pcoa_from_variants` of popnei on the panel with
`correct_by_lingoes` against pyNei's `do_pcoa` of popnei's
`correct_dists_by_lingoes` of the panel's distances, the 198 components.
Both compare the projections and the percentages within 1e-9, after the
test gives pyNei's components the sign of the rule; the review of this
spec found them 2e-14 apart on the first and 4.6e-14 on the second.
pyNei's components past those, of the eigenvalues 0, are not compared.

## The Rust interface

What an analysis gives. `num_comps` is how many components have variance.

```rust
pub struct Pca {
    /// How many columns the data had: the variants the first pass gave,
    /// used or not, which is `num_vars` of the pass stats, or the traits.
    pub num_cols: usize,
    pub num_rows: usize,
    pub num_comps: usize,
    /// num_rows x num_comps, row after row.
    pub projections: Vec<f64>,
    /// One per component, over the variance of every component of the data.
    pub explained_variance_percent: Vec<f64>,
    /// The positions of the columns that were used: the variants with
    /// variance among those the reader gave, or every trait of a table.
    pub used_cols: Vec<usize>,
    pub num_prin_comps: usize,
    /// num_prin_comps x used_cols.len(), row after row.
    pub princomps: Vec<f64>,
}
```

The PCA of the variants. The two readers are over the same variants,
opened by the binding crate from the source and the steps of the
`Variants`; opening a reader reads no variants. `second_pass` is `None` when
`num_prin_comps` is 0 and is not read then, and `None` with a
`num_prin_comps` above 0 is an error. The function asks each reader for
the genotypes alone and puts `reblock` before it.

```rust
pub struct VariantPcaOptions {
    pub transform_to_biallelic: bool,
    pub num_prin_comps: usize,
}

pub fn pca_of_variants<R1: BlockReader, R2: BlockReader>(
    first_pass: &mut R1,
    second_pass: Option<&mut R2>,
    options: &VariantPcaOptions,
) -> Result<Pca>;
```

The PCA of a table, `data` being `num_rows` x `num_cols`, row after row.
`num_prin_comps` of the result is `num_comps`.

```rust
pub struct PcaOptions {
    pub center: bool,
    pub standardize: bool,
}

pub fn pca(data: &[f64], num_rows: usize, num_cols: usize, options: &PcaOptions) -> Result<Pca>;
```

What a principal coordinate analysis gives. `num_comps` is how many
eigenvalues of B, or of the corrected B, are positive.

```rust
pub struct Pcoa {
    pub num_individuals: usize,
    pub num_comps: usize,
    /// num_individuals x num_comps, row after row.
    pub projections: Vec<f64>,
    /// One per component, over the sum of every eigenvalue of the B
    /// decomposed, the corrected one when there was a correction.
    pub explained_variance_percent: Vec<f64>,
    /// c of Lingoes' correction, 0 when there was none, which is always
    /// from `pcoa`, or the matrix was Euclidean.
    pub lingoes_constant: f64,
    /// 100 times the sum of |λ| of the negative eigenvalues of B before any
    /// correction over the sum of every eigenvalue of it.
    pub negative_eigenvalues_percent: f64,
}

```

The PCoA of a distance vector of `num_individuals` individuals, in the
order of `docs/specs/dists.md`, (0, 1), (0, 2), ..., (1, 2), .... It takes
the vector by value so that it can drop it before the eigendecomposition.
It refuses a matrix that is not Euclidean, and its result has a
`lingoes_constant` and a `negative_eigenvalues_percent` of 0.

```rust
pub fn pcoa(dist_vector: Vec<f64>, num_individuals: usize) -> Result<Pcoa>;
```

Lingoes' correction of a distance vector on its own: the corrected vector
in the same order, c, and the `negative_eigenvalues_percent` of the
vector given.

```rust
pub struct LingoesCorrection {
    pub dist_vector: Vec<f64>,
    pub constant: f64,
    pub negative_eigenvalues_percent: f64,
}

pub fn correct_dists_by_lingoes(dist_vector: Vec<f64>, num_individuals: usize) -> Result<LingoesCorrection>;
```

The PCoA of the Kosman distances of the variants of `reader`, which it
borrows, as `calc_kosman_sums` of `docs/specs/dists.md`, the pass of the
Kosman distances, does, so that the caller reads the counts of the filters
from it afterwards. `min_num_vars` is the `min_num_snps` of Python and
TypeScript, as `KosmanSums::dist` names it: a pair called together at
fewer variants has no distance, and 0 gives one to every pair called
together at one variant at least; the binding crate turns `None` into 0.
`num_vars` is how many variants the pass gave, for the pass stats.

`correct_by_lingoes` is false by default, as the owner decided on 27
September 2026, and the default is a constant of the core that both
packages read, as the PCA's `DEFAULT_TRANSFORM_TO_BIALLELIC`.

```rust
pub const DEFAULT_CORRECT_BY_LINGOES: bool = false;

pub struct VariantPcoaOptions {
    pub min_num_vars: u32,
    pub correct_by_lingoes: bool,
}

pub struct PcoaOfVariants {
    pub pcoa: Pcoa,
    pub num_vars: u64,
}

pub fn pcoa_of_variants<R: BlockReader + ?Sized>(
    reader: &mut R,
    options: &VariantPcoaOptions,
) -> Result<PcoaOfVariants>;
```

The error of the pairs with no distance carries how many pairs have none,
how many pairs there are, the positions of the two individuals of the
first, and the position of the individual in the most of them with how
many; each binding crate names them. The error of a matrix that is not
Euclidean carries how many eigenvalues are negative, n, and the
`negative_eigenvalues_percent`, and whether it came from the variants, so
that its message names `correct_by_lingoes` or `correct_dists_by_lingoes`.

What this module calls in `linalg`, whose spec settles the names and the
types: the product of a matrix of r rows and c columns with itself, added
to the lower half of a c x c matrix; the eigenvalues, from the largest,
and the eigenvectors of a symmetric matrix given by its lower half; and
the product of two matrices, for the second pass and for the projections
of a table with more rows than columns.

## Speed

Every time below is of one thread on the owner's Apple M5 Pro, rustc
1.98. Two kinds of number are mixed here and each one says which it is.
**Measured on this code** on 22 September 2026, the best of 5 runs, over
`big.vcf` and the vars file of it, 100000 variants x 1000 individuals
written by `crates/popnei/benches/make_big_vcf.py` with 3 in 100
genotypes missing, by the benchmark `crates/popnei/benches/pca_vars.rs`
and the script `time_pca.py` beside it; `docs/reports/pca-measurement.md`
has every one of those with the command that produced it. **From the
trial**, the best of 3 to 9 runs on 21 September 2026 from a crate kept
in `tmp/pca_trial/`, which is not in git, whose blocks have random
genotypes of two alleles, frequencies between 0.05 and 0.95, and 3 in 100
genotypes missing; those numbers are of options this code does not have,
and of wasm, which task 4.2 of `docs/plans/pca.md` measures.

What there is to beat, on 100000 variants x 1000 individuals, both
measured again on this dataset: pyNei's `do_pca_from_variants` takes
8.11 s and 5.59 GB, reading its own vars file, and 22.27 s reading the
VCF; plink2 v2.0.0-a.7.7 `--pca 10 meanimpute --threads 1` takes 0.248 s
from its own file of 2 bit genotypes, reading included, and 0.507 s from
the VCF, writing its own file first. `meanimpute` makes plink2 give a
missing genotype the mean of its variant, as pyNei does; without it
plink2 divides the product of each pair of individuals by the variants
that both have called, and takes 0.72 s, which is the trial's number.

The first pass has two steps per block, to standardize it and to add its
product. What this code takes for each of them, per block of 5000
variants x 1000 individuals, from the shares of a profile of the whole
analysis taken with `sample`:

| step, measured on this code | native, one thread |
|---|---|
| standardize the block, the four passes over each row | 20.6 ms |
| the product, lower half, through the linalg crate with Accelerate's `dsyrk` | 12.7 ms, of which 1.3 ms is the crate's two scans for a value that is not finite |

The options of the trial, for the same block, none of which this code
has except the loop of the second row, which it makes four passes of and
not two:

| step, from the trial | native | wasm under node 26 |
|---|---|---|
| standardize from the `i8` alleles, the loop as one would first write it | 3.9 ms | 7.1 ms |
| the same in two passes that the compiler vectorizes | 1.5 ms | 5.2 ms, and 1.8 ms with `simd128` |
| standardize from 2 bit genotypes | 1.0 ms | 1.4 ms |
| pack the `i8` alleles into 2 bits, not tuned | 5.2 ms | |
| the product, all of it, with faer | 174 ms | 597 ms, and 351 ms with `simd128` in faer |
| the product, lower half, with faer | 95 ms | 306 ms, and 187 ms with `simd128` |
| the product, lower half, with Accelerate's `dsyrk`, the BLAS routine for the product of a matrix with itself | 10.5 ms, and 7.4 ms with the threads Accelerate takes by itself | |

`simd128` is the set of 128 bit vector instructions of WebAssembly, which
a build turns on with `-C target-feature=+simd128`. faer's products are
done by the crate `gemm`, which uses those instructions only when its
cargo feature `wasm-simd128-enable` is on as well. At 1000 variants x 5000 individuals the product takes Accelerate
58 ms on one thread and the standardizing the same 1.3 ms.

**The 2 bit coding is not taken for the PCA.** In the trial it makes the
standardizing 0.5 ms shorter in a block that takes 12 ms natively, and
0.4 ms in one that takes 190 ms or more in wasm, and only when the block
comes packed from the file, since packing it costs a pass like the one it
saves. This code's block takes 33.3 ms natively, of which 20.6 ms is the
standardizing, so what the packed genotypes would save here has not been
measured and is not the trial's 0.5 ms. What the coding would change
besides is the reading: decompressing a block of `i8` genotypes from the
vars file takes 4.7 ms (`docs/specs/io_vars.md`), and the whole reading of
one block of this dataset 5.3 ms, and that question belongs to the vars
file and to section 4 of the architecture, where the option stays open.

**The vector instructions are taken twice, in safe code.** The core crate
forbids `unsafe`, and the intrinsics of `std::arch` need it, so the first
is a loop written for the compiler to vectorize. As it is written it makes
four passes over a row, since the major allele has to be known before any
dosage can be worked out: the counts of the alleles of the variant, which
give the major allele; a pass that writes the code of each genotype, 0, 1,
2 or missing, into a buffer of one byte per individual; a pass that counts
the codes, in runs of 255 genotypes with counters of one byte, each pass
over a run comparing the code with one dosage and adding; and a pass that
looks each code up in the four values it can take.

Of those four, `cargo asm`, which prints the machine code of a function,
shows for rustc 1.98 on this machine that the compiler vectorized one:
the counting of the codes, which compares 64 codes at a time with `cmeq`
and adds them with `udot`, and which is the path that runs at 1000
individuals. The pass that writes the codes was vectorized over the
alleles of one genotype, which is the ploidy, so at a ploidy of 2 no
vector instruction of it runs and the row is read one byte at a time. The
lookup and the counts of the alleles both read a table at an index that
is a value in memory, which NEON has no instruction for, and both are
scalar. The four passes together take 20.6 ms of the block, where the
trial's two took 1.5 ms; `docs/reports/pca-measurement.md` has what was
read in the machine code and where the time of the analysis goes, and
what to do about it is not decided here.

The second of the two is `simd128` in the wasm builds,
where the product takes 99 of every 100 ms of a block: computing the lower half
alone takes the block from 597 ms to 306 ms, and `simd128` from there to
187 ms (**Open 6**, below).

The number to reach, for 100000 variants x 1000 individuals with
`num_prin_comps` 0, was 0.3 s natively on one thread, which is 20 blocks
at 12 ms and an eigendecomposition of 0.04 s with Accelerate's LAPACK;
and 5 s in wasm with `simd128`, 7 s without.

**The two numbers of wasm are met.** The wasm package under node 26 takes
4.500 s with `num_prin_comps` 0 and 5.164 s with 10, and the wheel of
pyodide 4.481 s, over the vars file of those variants, whose bytes the
call is given; at 20000 variants, 1.058 s and 1.046 s. The build that
popnei ships is the one this section calls "with `simd128`", although it
sets no flag: what puts the vector instructions in is the cargo feature
`wasm-simd128-enable` of `gemm`, which the workspace turns on, and
`-C target-feature=+simd128` on top of it gives a module that is the same
file byte for byte. With that feature off, a build popnei does not ship,
the same analysis takes 7.066 s. A block of 5000 variants costs 215 ms in
wasm against 33.3 ms natively, and what is left of the analysis once the
blocks are paid, the eigendecomposition of 1000 individuals and the
building of the result, 0.197 s against the 0.035 s that
`docs/specs/linalg.md` gives the eigendecomposition alone natively.

**The native number is missed.** From the vars file on one thread this
code takes 0.801 s with `num_prin_comps` 0, 2.7 times the 0.3 s, and
1.353 s with 10, the second pass being the reading and the standardizing
again; with the threads the machine gives, 0.352 s and 0.537 s; from the
VCF on one thread, 1.329 s. Called from Python it takes 0.022 s more,
which is the building of the pandas frames. Of the 0.801 s, 0.41 s is the
standardizing of the 20 blocks, 0.25 s their product, 0.107 s the reading
of the vars file and 0.024 s the eigendecomposition; the target counted
none of the reading and 1.5 ms per block of standardizing where this code
takes 20.6 ms. It is 10.1 times faster than pyNei and 3.2 times slower
than plink2 on the same dataset, and it holds 0.21 GB against pyNei's
5.59 GB. `docs/reports/pca-measurement.md` has all of it.

The principal coordinates have no target of their own. Their time is the
pass of the Kosman distances, 0.110 s for 100000 variants x 1000
individuals of a vars file on the 18 cores of the owner's Apple M5 Pro
and 0.764 s on one thread (`docs/reports/perf-read-ahead-2026-09-25.md`),
and one eigendecomposition of an individuals x individuals matrix, which
is the PCA's: 0.024 s natively at 1000 individuals, above.

## Open points

The owner decides the first five here, and the sixth in the `linalg`
spec. Until then the implementer follows the
"meanwhile" of each.

**Open 1: whether `num_prin_comps` also cuts the projections.** The owner
decided that the weights are given for the first 10 components. The
projections and the percentages are small next to them, individuals x
individuals at most, and pyNei gives all. The options are to give all of
them, as pyNei, which at 10000 individuals is 800 MB of projections that
nobody plots; or to let `num_prin_comps` cut the three, with the
percentages still over the variance of all the components, which the
diagonal of G gives without the eigenvalues, and which would let LAPACK
compute only the first eigenvectors. Recommendation: cut the three, with
a default of 10, because a result is then a few kilobytes at any size.
Meanwhile the implementer cuts the weights alone, and cutting the other
two later takes a slice of two arrays.

**Open 2: a variant with no called genotype.** pyNei raises an error for
it and leaves out in silence a variant with one allele. The options are
to reproduce the error, which makes a user filter by missing data before
any PCA of a raw VCF, or to leave it out as one more variant with no
variance, which is seen in the columns of `princomps`. Recommendation:
leave it out. Meanwhile the implementer leaves it out.

**Open 3: the components with no variance.** pyNei gives `min(n, p)`
components, of which one or more have a variance of 1e-30 and weights
that are numerical noise. The options are to give them, which in popnei
takes a made up direction for their weights, since the division by
sqrt(λ) has none; or not to give them, which changes the shape of the
result, 199 components for 200 individuals. Recommendation: not to give
them. Meanwhile the implementer does not.

**Open 4: `standardize_data` or pyNei's `standarize_data`.** The glossary
keeps pyNei's name for what a user sees, and this one is a misspelling.
The options are pyNei's name, which keeps the scripts written for pyNei
running, or the right spelling. Recommendation: the right spelling, since
pyNei will not be used once popnei exists. Meanwhile the implementer
writes `standardize_data`.

**Open 5: n or n - 1 in the standard deviation.** pyNei divides by n, and
R's `prcomp`, the program most users would compare with, by n - 1, which
makes every projection of R 0.9975 of popnei's at 200 individuals and
changes no plot. The options are n, which keeps pyNei's numbers, or n - 1,
which gives R's and needs no factor in the tests. Recommendation:
reproduce pyNei, since nothing but the scale of the axes changes.
Meanwhile the implementer divides by n.

**Open 6: `simd128` in the two wasm builds.** The 128 bit vector
instructions take the product of a block, which is nearly all the time
of the PCA in the browser, from 306 ms to 187 ms in the trial, and the
whole analysis of 100000 variants x 1000 individuals from 7.066 s to
4.500 s in this code, which task 4.2 of `docs/plans/pca.md` measured. The
decision is about the build flag `-C target-feature=+simd128`, and that
flag turned out to change nothing: what puts those instructions in is the
cargo feature `wasm-simd128-enable` of `gemm`, which the workspace
manifest turns on, and with it on the module built with the flag and the
module built without it are the same file byte for byte, on rustc 1.98
and these dependencies. So the two builds already have the vector
instructions, and turning the flag on would gain nothing. What that means
for a browser without them was not tried, here or anywhere: such a
browser would refuse the module popnei ships today, not only a module
built with the flag. Those instructions are in Chrome since version 91
and in Firefox since 89, both of 2021, and in Safari since 16.4, of 2023,
which is from the release notes of the browsers. The decision is asked,
with the options and the recommendation, as Open 1 of
`docs/specs/linalg.md`; it is here for its numbers. Meanwhile nothing in
this module depends on it.

## Not in this spec

- Cailliez's correction of a matrix of distances that is not Euclidean,
  which the owner did not choose on 27 September 2026: "Lingoes'
  correction" of the principal coordinates.
- The two backends of the products and of the eigendecomposition, how BLAS
  is linked and the threads of the product: the `linalg` spec.
- The principal components of the kinship, which the GWAS uses as
  covariates: the `kinship` spec. They are of another matrix, the one
  plink2 standardizes by 2p(1 - p).
- A randomized PCA, which gets the first components without the
  individuals x individuals matrix in several passes, as plink2's
  `--pca approx` does. The exact one is within the sizes of the
  objectives natively: about 160 s of products and 50 s of
  eigendecomposition for a million variants of 10000 individuals, from
  the times above. It is not built.
- The 2 bit layout of a block: section 4 of the architecture, where it
  stays an option to measure for the reading and for the calculations on
  counts.
- The chromosome and the position of the variants of `princomps`. The
  result gives their positions among the variants given, as pyNei does.
