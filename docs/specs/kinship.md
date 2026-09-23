# The kinship module

23 September 2026. The kinship tells a user of popnei how related every pair
of their individuals is, from the genotypes alone. It is what the mixed
models of the association study take as the covariance of a random effect,
and it is also a way of looking at the structure of a panel on its own.
There is no code. This spec covers the whole `kinship` row of section 9 of
`docs/architecture.md`, the genomic relationship matrix, its per pair
denominators and its principal components.

It stands on three specs. `docs/specs/variant.md` gives the allele counts of
one variant and which allele is the major one. `docs/specs/pca.md` gives the
dosage of a genotype and what is done with a missing one, which this module
takes unchanged, and the rule that fixes the sign of a component. The
products and the eigendecomposition are `docs/specs/linalg.md`. What the
association study does with a kinship is `docs/specs/gwas.md`.

## The matrix

### What it gives

For two individuals, how much more of their genome they share than two
individuals drawn at random from the same panel would. It is the matrix of
VanRaden 2008, which plink2's `--make-rel` computes and which GCTA, the
program that estimates how much of a trait the genotypes explain, is built
on. An entry
off the diagonal is twice the coancestry of the pair: about 0.5 for full
sibs or for a parent and a child, about 0.25 for half sibs, and near 0 for
two individuals with no recent ancestor in common. An entry on the diagonal
is 1 plus the inbreeding of that individual. Entries below 0 are ordinary
and mean a pair less alike than the average pair of the panel, because the
whole matrix is measured against that average.

Each variant becomes one number per individual, its **dosage**: how many
alleles of the genotype are not the major allele of the variant. That is the
dosage of `docs/specs/pca.md`, with the same major allele, the most frequent
among the called alleles of the variant and the lowest numbered of two that
are equally frequent, and the same rule for a genotype with any allele
missing, which takes the mean of the dosages of its variant.

Each variant is then centered and divided by the standard deviation its
allele frequency gives it under Hardy Weinberg:

    p   = mean dosage of the variant / ploidy
    z   = (dosage - mean dosage) / sqrt(ploidy * p * (1 - p))

`z` is the **standardized dosage**. This divisor is not the standard
deviation of the dosages themselves, which is what `docs/specs/pca.md`
divides by; the two agree only when the genotypes are in Hardy Weinberg
proportions. It is the divisor that makes an entry twice a coancestry, and
it is what plink2 and GCTA use.

With `z` the standardized dosages, one row per variant and one column per
individual, the entry of the pair `i`, `j` is

    k[i, j] = sum over the variants of z[v, i] * z[v, j] / m[i, j]

where `v` runs over the variants that were used and `m[i, j]` is the **per
pair denominator**: how many of *those* variants have a called genotype in
both `i` and `j`. A variant that was dropped for having no variance is in
neither the sum nor the denominator, which is the `[is_poly]` of
`_KinshipCalc.calc_for_chunk` in `pynei/gwas.py`. Counting every variant
instead would give 4 and 3 in the worked example below where the right
numbers are 2 and 1, and every entry would move.

With nothing missing the denominator is one number for every pair, the
variants that were used, and the whole matrix is then one product: `z'z`,
the transpose of the standardized dosages multiplied by the standardized
dosages, individuals by individuals, divided by that number. With genotypes
missing it is still that product, divided entry by entry by a second matrix
of the same shape, the denominators.

Which allele of a variant is the major one changes no entry, because
counting the dosages from the other allele turns `z` into `-z` and the
product of two of them is unchanged.

### Its Python function, and its TypeScript one

```python
calc_kinship(variants: Variants, individuals: Sequence[str] | None = None) -> Kinship
```

`Kinship` is a frozen dataclass with two fields and one property.
`matrix` is a frame of individuals by individuals with the names as index
and columns, and `num_vars` is how many variants it was built from.
`individuals` is a property, the index of the matrix as a tuple, as pyNei's
`samples` is, so that the names are in one place and cannot disagree with
themselves. `pass_stats`, the `PassStats` of `docs/specs/variant.md`, which
says how many variants the steps of the `Variants` let through and what each
filter was given and kept, is a third field and is `None` for a `Kinship` a
user built by hand, since no pass produced it.

It has two methods. `principal_components(num_pcs)` is the item below.
`filter_individuals(individuals)` takes the rows and columns of some of
them, keeps `num_vars` and `pass_stats` as they were, and raises a
`ValueError` naming any individual that is not in the matrix, which is what
`Kinship.filter_samples` of pyNei does.

With `individuals` the matrix is of those individuals, in the order given,
and every frequency, mean and denominator is theirs: it is not the kinship
of the whole panel with rows and columns taken out. On the reference panel,
the kinship of the individuals 10 to 49 differs from the same 40 rows and
columns of the kinship of all 200 by up to 0.129, and it uses 1195 variants
where the whole panel uses 1200, the other 5 having no variance among those
40. An individual not in the `Variants` is a `ValueError` naming it.

It mirrors `calc_kinship` of `pynei/gwas.py`. The differences from pyNei,
which `docs/objectives.md` asks to be written down:

- `samples` is `individuals`, and `Kinship.samples` and
  `Kinship.filter_samples` are `individuals` and `filter_individuals`, the
  name `docs/glossary.md` gives.
- `num_threads` is not an argument. The owner decided on 22 September 2026,
  in `docs/specs/dists.md`, that no calculation of popnei has one.
- `pass_stats` is new, as it is for every consumer of a `Variants`.
- A `Kinship` can still be built by hand, from a matrix and `num_vars`, so
  that a user can bring the one plink2 or a pedigree gave them and pass it
  to `calc_gwas`. pyNei's fields are checked by nothing; whether popnei's
  are, and how far a matrix may be from symmetric, is **Open 4**, below.
- A variant with more than two alleles among its called genotypes is read
  with every allele that is not the major one counting the same, as pyNei
  does, and there is no `transform_to_biallelic` to refuse it. This is
  **Open 1**, below.
- The sign of a principal component is fixed, which is the item below.

In TypeScript it is `calcKinship(variants, {individuals})`. The result has
`individuals`, an array of names, `numVars`, `matrix` as a `Float64Array` of
individuals by individuals, row after row, and `passStats`, with the methods
`principalComponents(numPcs)` and `filterIndividuals(individuals)`.

### Missing genotypes, variants with no variance, and what pyNei asserts

A genotype with **any** allele missing has no dosage: it takes the mean of
its variant, so once the variant is centered it is 0 and pulls the pair
nowhere, and it does not count in the denominator of any pair it is in. A
half called genotype is missing under both rules, from `to_012` of
`pynei/variants.py`, which sets a genotype with any missing allele to the
missing value, and from `_KinshipCalc.calc_for_chunk` of `pynei/gwas.py`,
whose `is_called` is `~any(missing_mask, axis=2)`. So an individual with
much missing data has smaller denominators and the same kind of entry, which
is why the denominator is per pair and not one number.

A variant whose called genotypes all have the same dosage has no variance
and is left out, and `num_vars` counts those that were kept. Two kinds fall
here, and the second surprises: a variant with one allele, and a variant
where **every individual is heterozygous**, whose major allele frequency is
0.5 and which no filter by frequency catches. Its `z` would be 0 divided by
a divisor that is not 0, so keeping it would add 0 to every entry and 1 to
every denominator, which is not the same as leaving it out.

A pair of individuals with no variant called in both has a denominator of 0.
pyNei divides by it under `numpy.errstate(invalid="ignore", divide="ignore")`
and the entry is NaN; the matrix is returned with the NaN in it and nothing
is said. Measured on three individuals and two variants where the first
individual and the third are never called together: the entry of that pair
is `nan` and the rest of the matrix is finite (**Open 2**, below).

When no variant varies among the individuals, pyNei raises
`ValueError("No variant varies among the samples, there is no kinship")`.
popnei raises the same. A `Variants` whose steps let no variant through at
all is a different case, and it raises what every consumer of a `Variants`
already raises for it, with the counts of each filter in the message.

`test_kinship_matches_plink2` in `test/test_gwas.py` asserts the whole
matrix of both panels against plink2 and the seven and four entries that
"How it is verified" lists. `test_kinship_of_some_samples_and_threads`
asserts that a subset is not a slice, that the components tell the
subpopulations apart, and that an unknown individual raises.
`test_monomorphic_and_missing_variants` asserts that a dataset of 50
variants, of which one has a single allele and one has no called genotype,
gives `num_vars` of 48.

### What pyNei does that is odd

`Kinship.principal_components` takes `numpy.sqrt(numpy.abs(eigvals))`, the
square root of the absolute value of the eigenvalue. With genotypes missing
the matrix has eigenvalues below 0, so that absolute value is what turns an
imaginary number into a real one and hides it. On the reference panel with 3
in 100 genotypes missing the smallest eigenvalue is -0.0321 against a
largest of 17.27; with nothing missing it is -4.8e-15, which is rounding. It
is the per pair denominator that does it: dividing a matrix of products
entry by entry by different numbers does not keep it positive semidefinite.
popnei does not reproduce the absolute value (**Open 3**, below).

The docstring of `principal_components` says that the eigenvectors of the
kinship "are the principal components of the standardized genotypes, so this
is the PCA of the variants the kinship was calculated from". They are close
and they are not the same, because the two standardize by different numbers:
the kinship divides by `sqrt(ploidy * p * (1 - p))` and `do_pca_from_variants`
by the standard deviation of the dosages. Measured on the reference panel,
the absolute correlation between the two sets of projections is 0.99994,
0.99991 and 0.9896 for the first three components, and the PCA's are 32.98,
33.13 and 33.80 times longer, which is about the square root of the 1200
variants the kinship divides by, 34.6. popnei's doc comment says they agree under Hardy
Weinberg and not otherwise.

### How it runs

One pass over the blocks, each block taken as a matrix. Per block, with
rayon across its rows, every variant is turned into its standardized dosages
and the ones with no variance are dropped, which is the row pass of
`pca_of_variants` in `crates/popnei/src/pca.rs` with one number changed, the
divisor: that pass divides by the standard deviation of the dosages and this
one by `sqrt(ploidy * p * (1 - p))`. Everything else is the same, the codes
of the genotypes, the table of allele counts and the lookup of one
standardized value per code. Whether the two call one helper with the
divisor as an argument or keep two copies is the implementation plan's to
decide, and this spec asks only that a change to the dosage rule of
`docs/specs/pca.md` cannot leave the two disagreeing. Then two products, both in
`linalg`: the standardized dosages of the block multiplied by themselves
into an individuals by individuals accumulator, and, when any genotype of
the block is missing, the same of the matrix of ones and zeros that says
which genotypes were called, into a second accumulator of the same size.
With nothing missing the second product is skipped and the denominator is a
count, as `_KinshipCalc.calc_for_chunk` skips it. A block that has a missing
genotype after blocks that had none still has to carry those blocks'
variants: the matrix of denominators is allocated with the count so far in
**every** entry and not at zeros. pyNei gets this for nothing, because
`_KinshipCalc.reduce` adds a scalar count to a matrix and numpy broadcasts
it; an implementation that allocates zeros loses every variant before the
first missing genotype and gives entries too large by that much.

What is kept from one block to the next is those two accumulators and the
count of variants used: individuals by individuals, not growing with the
variants. At 1000 individuals that is 8 MB each and at 10000 it is 800 MB
each, and it is why the second accumulator is not allocated until a missing
genotype is seen. Section 5 of `docs/rust_core.md` leaves it open how far
this goes in a browser, where the heap is 32 bit and a kinship of 10000
individuals is already 800 MB of the 4 GB that addresses.

Two datasets are refused, both for reasons the PCA's row pass already
refuses them for and with the same limits, `MAX_PLOIDY_OF_THE_VARIANTS` and
`MAX_INDIVIDUALS_OF_THE_VARIANTS` of `crates/popnei/src/pca.rs`: a ploidy
above 254, because a genotype is written as one byte holding its dosage or
the missing code, and more than 46340 individuals, because the individuals
by individuals matrix would hold more values than the 2147483647 that BLAS
and LAPACK count in. A `num_pcs` of 0 gives a result with no components and
is not an error, as asking a PCA for none is not.

`reblock` goes before it, as it does before the PCA and for the same two
reasons: a filter leaves blocks of uneven size, and the sum over blocks is
in floating point, so the block boundaries decide the last bits of every
entry. With `reblock` the same variants give the same matrix whatever the
source gave.

### How it is verified

Against plink2 v2.0.0-a.7.7, the arm64 build of 18 September 2026, whose
`--make-rel square` writes this same matrix:

    plink2 --vcf <file> --make-rel square --out <prefix>

which writes `<prefix>.rel`, 200 lines of 200 numbers separated by tabs, and
`<prefix>.rel.id`, the individuals in that order. Run on 23 September 2026
by `tests/reference/kinship/make_reference.py`, which also checks every
literal below against what it has just produced, so that running it is a
check that this spec has not drifted.

The two datasets are the same 200 individuals, `s000` to `s199`, and 1200
biallelic diploid variants twice:

- `tests/reference/kinship/panel_called.vcf.gz`, every genotype called,
  written by that script from `test/gwas_reference/sim.vars` of pyNei, which
  popnei cannot read.
- `tests/reference/dists/panel.vcf.gz` of `docs/specs/dists.md`, already in
  the repository, the same panel with 3 in 100 of the genotypes missing
  whole. The kinship needs it because with nothing missing the per pair
  denominator is one number for every pair and the rule is never exercised.

plink2 kept the order of the VCF in both, checked against `.rel.id`.

Against pyNei's `calc_kinship` on the same two panels, the largest absolute
difference from plink2 is 4.95e-06 and 4.93e-06 over the 40000 entries, and
`num_vars` is 1200 in both. plink2 writes six significant digits, so an
entry near 1 is rounded by up to 5e-6, which is what those two numbers are.
The tests compare within 1e-5 absolute, one unit of the last digit plink2
prints for an entry of that size.

The literals of the cargo tests, made at `calc_kinship` of "The Rust
interface" on the two VCFs read with the VCF reader, from plink2 on 23
September 2026. The individuals `s000` to `s003` are full sibs, and so are
`s100` and `s101`:

| dataset | pair | plink2 |
|---|---|---|
| panel_called | s000, s000 | 1.09309 |
| panel_called | s001, s001 | 1.22825 |
| panel_called | s000, s001 | 0.648081 |
| panel_called | s000, s002 | 0.615611 |
| panel_called | s000, s004 | -0.0945533 |
| panel_called | s000, s199 | -0.0760273 |
| panel_called | s100, s101 | 0.604995 |
| panel | s000, s000 | 1.09626 |
| panel | s000, s001 | 0.650379 |
| panel | s000, s004 | -0.103505 |
| panel | s100, s101 | 0.604119 |

Against pyNei: both libraries run `calc_kinship` on the two panels, popnei
from the VCF and pyNei from its own vars file, and every entry has to agree
within 1e-12 relative, because the two add the variants in different orders.
`num_vars` has to be equal exactly.

The worked example, which becomes the first cargo test: 4 diploid
individuals, `i0` to `i3`, and 4 variants.

| variant | genotypes | dosages | kept |
|---|---|---|---|
| v0 | 0/0 0/1 1/1 0/1 | 0 1 2 1 | yes |
| v1 | 0/0 0/1 ./. 1/1 | 0 1 **1** 2 | yes |
| v2 | 0/1 0/1 0/1 0/1 | 1 1 1 1 | no, no variance |
| v3 | 0/0 0/0 0/0 0/0 | 0 0 0 0 | no, one allele |

The dosage in bold is the missing genotype of `i2` taking the mean of its
variant, which is 1. Both kept variants have a mean dosage of 1, so `p` is
0.5 and the divisor is `sqrt(2 * 0.5 * 0.5)`, and the standardized dosages
are `-sqrt(2)`, 0, `sqrt(2)`, 0 for `v0` and `-sqrt(2)`, 0, 0, `sqrt(2)` for
`v1`. `num_vars` is 2, and the denominators are 2 for every pair but those
with `i2`, which are 1. The matrix is whole numbers:

|  | i0 | i1 | i2 | i3 |
|---|---|---|---|---|
| **i0** | 2 | 0 | -2 | -1 |
| **i1** | 0 | 0 | 0 | 0 |
| **i2** | -2 | 0 | 2 | 0 |
| **i3** | -1 | 0 | 0 | 1 |

Run through pyNei at commit ef0ca6e, every entry is a whole number within
4.4e-16, so the test asserts them within 1e-12 absolute. `i1` is 0 against everyone because both of its genotypes are
heterozygous and their standardized dosage is 0. `i0` and `i3` have both
variants called, so their denominator is 2 and their -2 becomes -1, and only
`v1` contributes to the sum because `i3` is standardized to 0 at `v0`. `i0`
and `i2` have only `v0` called in both, so their denominator is 1 and their
-2 stays -2.

In TypeScript, `calcKinship` is tested under node against the four entries
of `s000` above on `panel_called.vcf.gz` and against the worked example.

## The principal components of the kinship

### What it gives

Where each individual falls along the directions in which the panel varies
most, taken from the kinship instead of from the variants. A user gets them
to give to `calc_gwas` as covariates, which is how an association study
accounts for population structure without a mixed model, and they cost an
eigendecomposition of a matrix that is already in hand rather than a second
pass over the variants.

With `lambda_j` the eigenvalues of the kinship from the largest and `u_j`
its eigenvectors, component `j` is `u_j * sqrt(lambda_j)`.

### Its Python function, and its TypeScript one

```python
Kinship.principal_components(num_pcs: int) -> pandas.DataFrame
```

A frame with the names of the individuals as index and the components as
columns, named as `docs/specs/pca.md` names them, `PC0`, `PC1`, with zeros
on the left to the width of the number of components. It is a frame indexed
by individual so that it can be joined to the covariates of `calc_gwas`. In
TypeScript, `principalComponents(numPcs)` gives a `Float64Array` of
individuals by components, row after row, and `numComps`.

It mirrors `Kinship.principal_components` of `pynei/gwas.py`. The
differences:

- The sign of each component is fixed by the rule of `docs/specs/pca.md`:
  in every component the projection with the largest absolute value is
  positive, ties within 64 units in the last place going to the first
  individual. pyNei gives whatever LAPACK gave. popnei has two backends and
  three builds, so without the rule Python natively, Python under pyodide
  and TypeScript give different signs for the same data.
- A component whose eigenvalue is not above `lambda_1 * n * 2.2e-16`, the
  tolerance of `docs/specs/pca.md` with `n` the individuals, is not given,
  so asking for more components than the panel has gives those it has and
  `num_comps` says how many. pyNei gives exactly `num_pcs` of them, and
  above the number of individuals it gives none: asking a kinship of 4
  individuals for 6 components raises `ValueError: Shape of passed values
  is (4, 4), indices imply (4, 6)` out of pandas, because
  `_create_pc_names` made 6 names for 4 columns (**Open 3**, below).

### How it is verified

There is no program outside the project to check these against: plink2's
`--pca` is the eigendecomposition of the same matrix but its own
normalization, and it was not run. What is checked instead:

- Against pyNei on both panels, the absolute value of every projection of
  the first 10 components within 1e-9 relative, the absolute value because
  pyNei does not fix the sign. Each component of popnei has to obey the
  sign rule exactly.
- That the three largest eigenvalues of `panel_called` are 17.26914116,
  12.44731524 and 3.35871258, from numpy 2.5.3 on 23 September 2026, within
  1e-9 relative. No function gives an eigenvalue, so the check is made at
  `principal_components`, on the sum of the squares of each component's
  projections: a component is `u_j * sqrt(lambda_j)` and `u_j` has length 1,
  so that sum is `lambda_j` itself.
- That the first component tells the three subpopulations of the panel
  apart, which is `test_kinship_of_some_samples_and_threads` of pyNei: the
  standard deviation of the mean of `PC0` over the three subpopulations is
  above the standard deviation of `PC0` itself. Which subpopulation each
  individual belongs to is the `pop` column of
  `tests/reference/gwas/phenotypes.csv`.

## The Rust interface

The matrix, which the pass gives away so that the binding crate hands it to
numpy without copying it.

```rust
pub struct Kinship {
    pub num_individuals: usize,
    /// How many variants had variance among these individuals and were used.
    pub num_vars: usize,
    /// num_individuals x num_individuals, row after row, symmetric.
    pub matrix: Vec<f64>,
}
```

One pass over a reader. `individuals` are the positions among those the
reader gives, in the order the result has them, and `None` is all of them in
the reader's order. The pass borrows the reader and does not take it, so
that whoever built the chain of filters reads its counts when it returns; it
asks for the genotypes alone and puts `reblock` before it.

```rust
pub fn calc_kinship<R: BlockReader>(
    reader: &mut R,
    individuals: Option<&[usize]>,
) -> Result<Kinship>;
```

The components of a kinship, `num_pcs` of them at most and fewer when the
matrix has fewer with variance. `projections` is individuals x `num_comps`,
row after row, with the sign rule of `docs/specs/pca.md` applied.

```rust
pub struct KinshipPcs {
    pub num_comps: usize,
    pub projections: Vec<f64>,
}

pub fn principal_components(kinship: &Kinship, num_pcs: usize) -> Result<KinshipPcs>;
```

What this module calls in `linalg`: the product of a matrix of `r` rows and
`c` columns with itself, added to the lower half of a `c` x `c` matrix,
which is the one the PCA already calls; and the eigenvalues, from the
largest, and eigenvectors of a symmetric matrix given by its lower half.

## Speed

Most of the kinship is one matrix product and there is little room in it.
From the table of section 2.1 of `docs/rust_core.md`, over 100000 variants
x 1000 individuals, pyNei takes 0.81 s and plink2 0.23 s.

What the product alone costs was measured on 23 September 2026 on the
owner's Apple M5 Pro, numpy 2.5.3 on Accelerate with its threads, adding
`z'z` over 20 blocks of 5000 variants x 1000 individuals into one
accumulator: 0.157 s, 7.4 ms a block. So plink2's 0.23 s for the whole
calculation is close to what the product costs here, and popnei has about
0.65 s of pyNei's 0.81 s to win, all of it in the row pass that turns
genotypes into standardized dosages.

When any genotype is missing there is a second product of the same shape,
the matrix of called genotypes with itself, for the per pair denominators.
Measured the same way with 3 in 100 genotypes missing it takes 0.149 s, so
it nearly doubles the BLAS work, which is why it is skipped for a block with
nothing missing.

The number to reach is plink2's 0.23 s on the panel with every genotype
called, and the run with 3 in 100 missing is measured and reported beside
it. Neither has been measured for popnei; the measurement comes when the
code exists.

## Open points

The owner decides these four, and until then the implementer follows the
"meanwhile" of each.

**Open 1: a variant with more than two alleles.** pyNei's kinship reads one
with every allele that is not the major one counting the same, silently,
because `to_012` does it for everything. `do_pca_from_variants` refuses such
a variant unless its `transform_to_biallelic` argument says otherwise, which
`docs/specs/pca.md` settled. So a user with a multiallelic dataset would get
an error from `do_pca_from_variants` and a silent answer from
`calc_kinship`, for the same reason and on the same variants. The options
are to reproduce pyNei, which keeps the kinship of every multiallelic
dataset as it is and leaves the two functions disagreeing; or to give
`calc_kinship` the same `transform_to_biallelic` argument, defaulting to
false, which makes the two agree and refuses datasets pyNei answers for, and
which no reference verifies since both panels are biallelic.
Recommendation: give it the argument. `docs/objectives.md` asks the
calculations that collapse a multiallelic variant to the major allele
against the rest to say so, and an error a user can turn off with one
argument says it; the cost is that the same argument now appears on two
functions instead of one.
Meanwhile the implementer reproduces pyNei and collapses silently, which is
what the reference panels need.

**Open 2: a pair with no variant called in both.** Its denominator is 0,
pyNei puts a NaN in the matrix and says nothing, and popnei would too. The
NaN does not stay silent: "Errors" of `docs/specs/linalg.md` refuses a value
that is not finite in any matrix before any routine runs, so the mixed model
fit that factorizes the kinship, and `principal_components`, both raise a
`RuntimeError` saying that a matrix holds a value that is not finite. So the
choice is not between a wrong answer and an error; it is about where the
error appears and what it names. The options are to reproduce pyNei, which
gives the user a matrix they can look at and see the gap in, and an error
later that names a matrix; or to raise here, naming the two individuals and
how many variants each of them has called, which refuses a dataset pyNei
answers for. Recommendation: raise here. The user can do something with the
names of two individuals and nothing with the word "matrix", and a kinship
they can look at is of little use when what they wanted was to pass it to
`calc_gwas`. Meanwhile the implementer raises, since no reference panel has
such a pair and no literal moves either way.

**Open 3: the components of a kinship that is not positive semidefinite.**
Its per pair denominators can leave eigenvalues below 0, -0.0321 against a
largest of 17.27 on the panel with 3 in 100 genotypes missing. pyNei asks
for `num_pcs` components whatever their eigenvalues and takes the square
root of the absolute value. The options are to reproduce pyNei, which always
gives the number asked for and turns a negative eigenvalue into a length
that means nothing; or to give only the components whose eigenvalue is above
the tolerance of `docs/specs/pca.md`, which is what that spec decided for
the same reason and which gives fewer components than asked on a panel with
much missing data. Recommendation: the tolerance, with `num_comps` saying
how many came back. Nobody asks for the hundredth component of a kinship,
and the ones a user does ask for are far above the tolerance: the third
eigenvalue of the panel is 3.36. Meanwhile the implementer applies the
tolerance.

**Open 4: what a `Kinship` a user built by hand is checked for.** pyNei's
`Kinship` is a frozen dataclass whose fields nothing checks, so a frame that
is not square, whose index and columns name different individuals, or that
is not symmetric, is taken and reaches the mixed model fit. The options are
to reproduce pyNei and check nothing, which accepts everything it accepts
and pushes every complaint into the linear algebra, where the message names
a matrix and a row and not the field the user filled; or to refuse in
`__post_init__` a matrix that is not square, whose index and columns differ,
or that is further from its own transpose than a tolerance. Recommendation:
refuse, with the tolerance at 1e-9 of the largest absolute entry. What the
tolerance costs a real user was measured: both matrices plink2 wrote for the
panels are symmetric to the bit, largest `|m - m'|` of 0, so a matrix that
came from a tool is not near the tolerance, and one that is further than
1e-9 from symmetric was built by an arithmetic that is not a kinship's.
Meanwhile the implementer refuses, at that tolerance.

## Not in this spec

- What the association study does with a kinship: which of its four null
  models takes one, what it costs, and what happens when the fit meets a
  kinship with an eigenvalue below 0, which `docs/reports/glmm-method/README.md`
  measured and found does not stop the fit at any missing rate tried.
  `docs/specs/gwas.md`.
- The dosage of a genotype, the major allele, and what a missing one gets:
  `docs/specs/pca.md`, which this spec takes unchanged and does not repeat.
- The sign rule for a component, its tolerance and why it exists:
  `docs/specs/pca.md`.
- The products and the eigendecomposition, their backends and their errors:
  `docs/specs/linalg.md`.
- The pruning by linkage disequilibrium that a user normally does before a
  kinship, so that a block of correlated variants does not count many times:
  the filter of `docs/specs/filters.md`. It is a step of the `Variants` and
  this module does nothing about it.
- A kinship from a pedigree rather than from genotypes, and a sparse kinship
  of the kind a biobank uses. popnei does not have either.
- Jost's D and the Kosman distance, the other two ways popnei measures how
  alike individuals or populations are: `docs/specs/dists.md`.
