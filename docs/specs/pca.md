# The pca module: principal components of variants and of a table

September 2026. A principal component analysis places the individuals of a
dataset on a few axes that hold as much of the variation between them as
that many axes can, which is how a user sees whether their individuals
fall into populations before they name any. This module gives two: the
PCA of the variants, from the genotypes of a `Variants`, and the PCA of
any table of numbers, individuals by traits. There is no code. This spec
develops the row `pca` of the table in section 9 of `docs/architecture.md`,
without the principal coordinates, which are not written yet.

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
  gives. popnei's messages are "There are no variants to do a PCA with"
  and "Every variant has the same genotype in every individual, there is
  nothing to do a PCA with", pyNei's without its "012 matrix" and its
  "sample".
- A negative `num_prin_comps`.

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
and one bit per variant for whether it was used. When the variants end,
`linalg` gives the eigenvalues and the eigenvectors of G.

The second pass, when `num_prin_comps` is above 0, standardizes each block
again and multiplies it by the first eigenvectors divided by sqrt(λ), a
product of variants x individuals by individuals x `num_prin_comps`. It
checks that the variants it uses are those of the first pass, and stops
with an error if they are not, which is what a file that changed between
the two gives; in Python a `RuntimeError`, since no argument is wrong and
the core has no file name to give.

Memory: G and its eigenvectors, 16 MB at 1000 individuals and 1.6 GB at
10000, which no browser tab holds; the buffer of a block; and the
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

The literals of the two tables and of `worked3` are checked at
`pca_of_variants` of "The Rust interface", on a reader over blocks that
the test builds, once with `num_prin_comps` 3 and once with 0. No function
under it is pinned by a test.

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
result that of `doPcaFromVariants` with `usedVars` left out. The names of
the rows and of the traits stay with the application.

### Errors and the cases pyNei asserts

Each is a `ValueError` in Python: a value that is not finite;
`standardize_data` with `center_data` false; fewer than 2 rows or no
columns; and, when standardizing, a trait with no variance, with a message
that gives how many there are and names the first ten, as pyNei's does, so
that the user can take them out. The core gives the positions of those
columns and the Python layer puts the names. Without standardizing such a
trait is no error and gets a weight of 0.

`test_pca` asserts the four princomps of iris, up to sign, and
`test_pca_refuses_traits_with_no_variance` the error, that its message
names the trait, and that without standardizing the same table gives 3
projections x 3 components. popnei gives that table 2 components, since
the third has no variance, and the test says so.

### How it runs

On the whole table, which is copied once to be centered and standardized.
The product is over the smaller side, as "What both analyses compute"
says, so its memory is that of the table plus the square of its smaller
side. Nothing is kept and there are no blocks.

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

What this module calls in `linalg`, whose spec settles the names and the
types: the product of a matrix of r rows and c columns with itself, added
to the lower half of a c x c matrix; the eigenvalues, from the largest,
and the eigenvectors of a symmetric matrix given by its lower half; and
the product of two matrices, for the second pass and for the projections
of a table with more rows than columns.

## Speed

Every time below is of one thread on the owner's Apple M5 Pro, rustc
1.98, the best of 3 to 9 runs, on 21 September 2026, from a trial crate
kept in `tmp/pca_trial/`, which is not in git. Its blocks have random
genotypes of two alleles, frequencies between 0.05 and 0.95, and 3 in 100
genotypes missing.

What there is to beat, on 100000 variants x 1000 individuals: pyNei's
`do_pca_from_variants` takes 8.2 s and 6.7 GB from genotypes in memory;
plink2 v2.0.0-a.7.7 `--pca 10 meanimpute --threads 1` takes 0.26 s from
its own file of 2 bit genotypes, reading included. `meanimpute` makes
plink2 give a missing genotype the mean of its variant, as pyNei does;
without it plink2 divides the product of each pair of individuals by the
variants that both have called, and takes 0.72 s.

The first pass has two steps per block, to standardize it and to add its
product. Per block of 5000 variants x 1000 individuals:

| step | native | wasm under node 26 |
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

**The 2 bit coding is not taken for the PCA.** It makes the standardizing
0.5 ms shorter in a block that takes 12 ms natively, and 0.4 ms in one
that takes 190 ms or more in wasm, and only when the block comes packed
from the file, since packing it costs a pass like the one it saves. What
it would change is the reading: decompressing a block of `i8` genotypes
from the vars file takes 4.7 ms (`docs/specs/io_vars.md`), as much as the
two steps of the PCA natively, and that question belongs to the vars file
and to section 4 of the architecture, where the option stays open.

**The vector instructions are taken twice, in safe code.** The core crate
forbids `unsafe`, and the intrinsics of `std::arch` need it, so the first
is a loop written for the compiler to vectorize: one pass that writes the
code of each genotype, 0, 1, 2 or missing, into a buffer of one byte per
individual and counts the codes in runs of 255 with counters of one byte,
and one pass that looks each code up in the four values it can take. It is
2.7 times faster than the plain loop, 1.5 ms against 3.9 ms, and within
0.5 ms of what the 2 bit coding gives. Whether the compiler did vectorize
it is checked when it is written with `cargo asm`, which prints the machine
code of a function, since what the compiler vectorizes changes with the
version of rustc. The second is `simd128` in the wasm builds,
where the product takes 99 of every 100 ms of a block: computing the lower half
alone takes the block from 597 ms to 306 ms, and `simd128` from there to
187 ms (**Open 6**, below).

The number to reach, for 100000 variants x 1000 individuals from blocks in
memory with `num_prin_comps` 0: 0.3 s natively on one thread, which is
20 blocks at 12 ms and an eigendecomposition of 0.04 s with Accelerate's
LAPACK, and 27 times less than pyNei; and 5 s in wasm with
`simd128`, 7 s without. The second pass has not been measured, and is the
reading and the standardizing again. These are the times of a trial, and they are measured again on the code
before they are claimed.

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
of the PCA in the browser, from 306 ms to 187 ms. A build with them does
not load in a browser without them, which are in Chrome since version
91 and in Firefox since 89, both of 2021, and in Safari since 16.4, of
2023; that is from the release notes of the browsers and was not tried.
The pyodide wheel built with them on 22 September 2026 loaded under
pyodide on node 26 and passed the smoke test of `tests/pyodide/`. It is
a decision about the builds, and it is asked, with the options and the
recommendation, as Open 1 of `docs/specs/linalg.md`; it is here for its
numbers. Meanwhile nothing in this module depends on it.

## Not in this spec

- The principal coordinates, `do_pcoa` and `do_pcoa_from_variants`, which
  need the distances of `docs/specs/dists.md`: a later item of this spec.
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
