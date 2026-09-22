# The ld module: how strongly two variants go together, and how that falls off with distance

22 September 2026. Two variants are in linkage disequilibrium when the
genotype of one tells something about the genotype of the other, which
happens when they sit close enough on a chromosome that few
recombinations have separated them. The `ld` module measures it, r² for
every pair of a set of variants, and it gives the curve of that measure
against the distance between the variants, for each population on its
own, which is what a geneticist reads the recombination of a genome and
the history of a population from. There is no code. This spec develops
the row `ld` of the table in section 9 of `docs/architecture.md`.

It depends on `docs/specs/block.md`, which has the block, the run of
consecutive variants held as arrays, and `BlockReader`, the trait of
everything that gives blocks; on `docs/specs/variant.md`, which has the
`Variants` a user holds, the steps put on it and the counts of a pass; on
`docs/specs/filters.md`, which builds the chain of filters of a pass; on
`docs/specs/pca.md`, which defines the dosage that both modules read the
genotypes as; and on `docs/specs/linalg.md`, which has the matrix product
that runs on the BLAS of the system natively and, in WebAssembly, where
there is none, on faer, a linear algebra library written in Rust. The
filter that thins variants out by their r² is an item of
`docs/specs/filters.md`, written with this spec, because a filter is a
reader over a reader and every other filter is there.

## r² between two variants

### What it gives

Each variant becomes one number per individual, its dosage: how many
alleles of the genotype are not the major allele of the variant, 0, 1 or
2 in a diploid. The major allele is the most frequent among the called
alleles of the variant, those of half called genotypes included, and the
lowest numbered of two that are equally frequent, as
`docs/specs/pca.md` defines it. Every allele that is not the major one
counts the same, so a variant of more than two alleles is read as two
(**Open 2**, below). A genotype with any allele missing, `0/.` included,
has no dosage.

r² is the square of the correlation, across the individuals, between the
dosages of the two variants. It is 1 when the dosage of one individual at
one variant fixes its dosage at the other, and 0 when knowing one says
nothing about the other. It is the estimate that popnei and pyNei name
after Rogers and Huff, the correlation between the dosages of two
variants whose phase is unknown, which estimates the correlation between
the alleles that travel together on a chromosome without any of the
genotypes having been phased; `_get_r` of pyNei's `test/test_ld.py`
computes it from that definition.

The individuals that count for a pair are those whose genotype is called
at both of its variants. With x the dosages of one variant and y of the
other over those individuals, and n how many of them there are, every sum
below running over those individuals alone,

    r² = (n·Σxy − Σx·Σy)² / ((n·Σxx − (Σx)²) · (n·Σyy − (Σy)²))

A pair has no r², and the result is NaN, when n is 0, and when every x is
equal or every y is equal among those n individuals, which makes one of
the two factors below the line 0. A variant whose called genotypes all
have the same dosage, a variant with one allele among them and also one
where every individual is heterozygous, therefore has NaN against every
other variant and against itself.

Every one of the six sums is a whole number, and each is held exactly in
an `f64`: the largest of them is n·Σxx, which is at most N²k² for N
individuals of ploidy k, 4·10⁸ for 10000 diploid individuals and 6.5·10¹²
at the ploidy of 255 that the VCF reader takes, against the 9·10¹⁵ up to
which an `f64` counts whole numbers one by one. So the only rounding in
an r² is in the division and the square at the end.

### Its Python function

```python
calc_rogers_huff_r2_matrix(variants: Variants,
                           max_num_vars: int = 5000) -> R2Matrix
```

It is a consumer of the `Variants`, as `docs/specs/variant.md` has them:
it makes one pass over the source through the steps the `Variants` has
when it is called, and the `Variants` is as it was afterwards.

`R2Matrix` is a frozen dataclass with `r2`, a read only square numpy
array of float64 with as many rows and columns as the pass gave variants,
NaN where a pair has no r² and 1 on the diagonal of a variant that has
one; `chroms`, a tuple with the name of the chromosome of each variant;
`poss`, a read only numpy array of the positions; and `pass_stats`, the
`PassStats` of `docs/specs/variant.md` that every result of a consumer
has, how many variants the calculation took and how many each filter was
given and kept.

`max_num_vars` is how many variants the calculation takes before it
refuses: the matrix holds the square of them, 200 MB of float64 at the
default of 5000 and 80 GB at 100000, and it is the one calculation of
popnei whose result grows with the square of its input, against goal 5 of
`docs/objectives.md`, under which a dataset never has to fit in memory. A
user who wants the matrix of more variants than that, and has the memory,
raises it. A user who has more variants than memory puts a filter on the
`Variants` first, or asks for the curve of the next item instead.

It mirrors `calc_rogers_huff_r2_matrix` of `pynei/ld.py`. The
differences:

- **It gives r² where pyNei gives r.** pyNei's function of this name, and
  the `r2` field of its result, hold the correlation itself and not its
  square (under "What pyNei does that is odd"). The owner decided on 22
  September 2026 that every function of this module gives r², so that the
  name and the value agree; the option not taken was to give r under a
  name that says r. What is lost is the sign, which says whether the
  major alleles of the two variants go together or apart, and which no
  caller of pyNei reads.
- **A missing genotype takes its individual out of that pair.** pyNei
  leaves it in with a dosage of -1 (under "Missing genotypes"). The owner
  decided on 22 September 2026 to drop it, which is what plink2 does and
  what makes plink2 the reference program of this module; the option not
  taken was pyNei's, and the other was to give a missing genotype the
  mean dosage of its variant, as `docs/specs/pca.md` does for the PCA.
- **`max_num_vars` is new.** pyNei holds every chunk of the dataset in
  memory and builds the matrix of all of them with nothing to stop it.
- **There is no `dists_in_bp`.** pyNei gives a second square matrix with
  the distance between every pair of variants and NaN where the two are
  on different chromosomes. `chroms` and `poss` hold the same thing in 2n
  numbers instead of n², 80 KB against 200 MB at 5000 variants, and the
  distance of a pair is the difference of two positions. Decided here.
- **There is no `max_dist`.** pyNei leaves NaN in the cells of the pairs
  further apart than it, which saves it the products of whole chunk pairs
  and nothing of the matrix. Here the matrix is capped instead, and the
  pairs within a distance are what the next item gives. Decided here.
- **There is no `check_no_mafs_above`.** pyNei raises a `ValueError` for
  the whole call when any variant has a major allele frequency above
  0.95, since the r of a variant that hardly varies rests on one or two
  individuals. In popnei that is what `variants.filter_by_maf(0.95)`
  does, a step of the `Variants` that takes those variants out instead of
  refusing the dataset, and whose frequency is the one
  `docs/specs/filters.md` verifies against bcftools at any ploidy, where
  the guard's own `_calc_maf_from_012_gts` reads the frequency off the
  dosages in a way that is one only for a diploid biallelic variant.
  Decided here.
- **The result has `pass_stats`**, which pyNei's has not.

In TypeScript it is `calcRogersHuffR2Matrix(variants, {maxNumVars})`,
which gives an `R2Matrix` with `numVars`, `r2` as a `Float64Array` of
`numVars` x `numVars` row after row, `chroms` as an array of strings,
`poss` as a `Float64Array`, which holds a position of up to 2⁵³ exactly
as the blocks of `js/popnei/src/block.ts` already give them, and
`passStats`.

### Missing genotypes, variants with no variance, and what pyNei asserts

A genotype with any allele missing is a missing genotype, as
`docs/glossary.md` has it, and the individual that holds it is left out
of every pair that variant is in. Two variants whose missing genotypes
fall on different individuals therefore have an n below the individuals
of the dataset, and each pair of a dataset has an n of its own.

On the panel of `tests/reference/dists/panel.vcf.gz`, 200 individuals and
1200 biallelic variants with 3 in 100 genotypes missing, the three rules
give these r², measured on 22 September 2026 against plink2 v2.0.0-a.7.7
over the 1.4 million pairs off the diagonal. The median r² of that panel
is 0.0071.

| the rule | difference from plink2 in r², median | 99th percentile | largest |
|---|---|---|---|
| the individual is left out of the pair | 0 | 0 | 0 |
| the genotype takes the mean dosage of its variant | 0.00035 | 0.0078 | 0.042 |
| pyNei: the genotype is a dosage of -1 | 0.0037 | 0.047 | 0.194 |

The first row is 0 and not a rounding: every pair of that panel comes out
of the six whole numbers of "What it gives" with the bits plink2 has. The
same comparison made on r instead of r² does not reach 0, for a reason
that has nothing to do with missing genotypes: variant 856 of the panel
has 197 of its 394 called alleles of one allele and 197 of the other,
popnei's tie rule makes the lower numbered one major and plink2 makes the
other, and the r of every pair that variant is in then has the other
sign. r² is the same whichever allele is called major, so this module
never sees it.

A variant with no variance, whose called genotypes all have the same
dosage, has NaN against everything, including against itself, so its row
and its column of the matrix are NaN and its diagonal cell is NaN and not
1. plink2 gives NaN there too. In the reference dataset of "How it is
verified", 68 of the 500 variants are of this kind.

A pair can have NaN although both of its variants have variance, when the
individuals called at both happen to hold one dosage at one of them.

`test/test_ld.py` of pyNei asserts, for this item:
`test_the_r_matrix_is_the_one_numpy_cov_gives` pins the matrix to what
`numpy.cov` gave, NaN cells included, for 25 variants against 17 of 60
individuals with one variant of no variance;
`test_the_r_matrix_does_not_build_the_square_of_both_sets_together` that
one variant against 20000 gives a 1 x 20000 result and not a 20001
square; `test_r2_matrix_with_chunks_of_different_sizes` that the matrix
does not change when the chunks are cut elsewhere, which popnei takes as
a test over block sizes; and `test_ld_calc` that a variant with a major
allele frequency above the guard raises. popnei takes the first, the
third, and the shape of the second, at its own function; the guard is
gone, as said above. `test_ld_calc` does put missing genotypes through
`_calc_rogers_huff_r2`, 963 of its 100000 dosages, since its `geno_freqs`
gives `(-1, -1)` a weight of 0.01, but it asserts nothing that they could
move: it compares r with the rate at which its variants were made
independent, within 0.3.

### What pyNei does that is odd

Read and run in pyNei at commit ef0ca6e.

`_calc_rogers_huff_r2` of `pynei/ld_calc.py` returns r and not r². Its
last line divides the products of the centred dosages by the square root
of the product of the two sums of squares, which is the correlation, and
its value is signed: of the 1.4 million cells off the diagonal of the
panel above, 725934 are below 0, and they run from -0.5403 to 0.5761. The name of the function, the `r2` field of `R2Matrix` and
the `r2` field of `LDResult` all say r², and
`test_the_r_matrix_is_the_one_numpy_cov_gives` compares it with a
correlation, so the tests hold the value and the names do not.
`test_ld_calc` takes the square root of it and compares the result with
the rate at which its simulated variants were made independent, within
0.3, which passes for r and for r² alike at the three rates it uses.

`to_012` of `pynei/variants.py` writes -1 for a genotype with an allele
missing, and `_center_gts_for_r` of `ld_calc.py` then subtracts the mean
of the row from that -1 and keeps it among the values, so a missing
genotype enters the correlation as a dosage below every real one. The
numbers that rule gives are in the table above.

### How it runs

One pass. The calculation asks its reader for the genotypes and for the
chromosome and the position, and stops with an error when the variants
pass `max_num_vars`.

The six sums of every pair come out of six matrix products, which is what
lets the whole matrix be computed by the linear algebra of
`docs/specs/linalg.md` instead of pair by pair. For the variants of the
pass, three matrices of variants x individuals in `f64` are built: A, the
dosages with a missing genotype written as 0; M, 1 where the genotype is
called and 0 where it is not; and S, the square of each entry of A. Then,
with ' the transpose,

    n = M M'      Σxy = A A'      Σx = A M'      Σy = M A'
    Σxx = S M'    Σyy = M S'

each of which is a matrix over the pairs, and the r² of every pair is the
formula of "What it gives" applied to the six entries. When the two sets
of variants are the same, n and Σxy are each their own transpose and Σy
and Σyy are the transposes of Σx and Σxx, so a block of the matrix
against itself costs four products and one against another six. r² is
symmetric, so only the lower half of the whole matrix is computed and the
upper half is a copy of it.

The products are made in tiles of a few hundred variants and not on the
whole set at once, because six intermediate matrices of the whole set
would be six times the result, 1.2 GB at `max_num_vars` of 5000, where
six tiles of 512 variants are 12 MB. The memory of the calculation is
then the result, 8 bytes per pair, 200 MB at 5000 variants, the three
matrices of the dosages, 24 bytes per variant and individual, 120 MB for
5000 variants of 1000 individuals, and the tiles.

When no genotype of either tile is missing, five of the six products give
nothing that a sum over each row does not. n is then the individuals for
every pair; Σx of the pair (i, j) is the sum of row i of A, whatever j
is, and Σy is the sum of row j; and Σxx and Σyy are the sums of the rows
of S in the same way. So five of the six matrices are one number per
variant read across a row or down a column, and only A A' is a product.
Whether that path is built is for the implementation plan that builds
this module to decide on a measurement, as the plans under `docs/plans/`
decide what a spec leaves to a measurement. The numbers of "Speed" are of
the six products whether it is built or not, since it can only take less
time than they do.

A product is called from outside rayon, as section 3 of the architecture
asks, so that the threads of the backend of the linear algebra and those
of popnei's own loops do not nest.

### How it is verified

The reference program is plink2 v2.0.0-a.7.7, whose `--r2-unphased` is
the squared correlation between the dosage vectors of two variants with
the individuals missing at either one left out. Run on 22 September 2026
on the two datasets below, its numbers and the formula of "What it gives"
agree to 5.6e-16 over 1321 pairs of the LD dataset and to 1.7e-16 over
the pairs of the panel, so the two are the same calculation and the
difference is the rounding of the last operations. The matrix is asked
for as float64 and not as text, with

    plink2 --vcf <file> --double-id --allow-extra-chr \
           --r2-unphased square bin --out <prefix>

which writes `<prefix>.unphased.vcor2.bin`, the square matrix row after
row; the text form plink2 writes by default has six digits.

The first dataset is new, `tests/reference/ld/ld.vcf.gz`, because the
files the repository has carry no linkage disequilibrium: their variants
were made independently, so every r² of them is the noise of the sample
size. It has two chromosomes of 250 biallelic variants each, 1000 bp
apart, and 100 diploid individuals, and its haplotypes come from four
founders recombined along the chromosome at a rate of 2 in 100 between
one variant and the next, with 3 in 100 genotypes then set to missing.
The script that writes it is `docs/reports/ld-method/make_ld.py`, which
the implementation plan moves to `tests/reference/ld/make_reference.py`
unchanged. Every literal of this spec and of the filter item of
`docs/specs/filters.md` depends on the order in which that script asks
`numpy.random.default_rng(29)` for its numbers, the four founders, then
the recombination of each haplotype along the chromosome, then the mask
of the missing genotypes, so a script written afresh from this paragraph
gives another file and fails every literal. The parameters are here so
that a reader knows what the dataset is, not so that it can be built
again from them. Its mean r² falls with the distance, from
0.208 over the pairs up to 25000 bp apart to 0.013 over those between
225001 and 250000, which is the table of the next item, and is 0.0104
between the two chromosomes, where nothing links the variants and what is
left is the 1/(n − 1) of 100 individuals, 0.0101. 68 of its 500 variants
have no variance.

The second is `tests/reference/dists/panel.vcf.gz` of
`docs/specs/dists.md`, 200 individuals and 1200 biallelic variants with 3
in 100 genotypes missing, which is where the table of "Missing genotypes"
was measured.

The literals of the cargo tests, made at `calc_r2_matrix` of "The Rust
interface" on `ld.vcf.gz` read with the VCF reader, from plink2 on 22
September 2026:

| pair, by position | r² | n |
|---|---|---|
| chr1:1000, chr1:2000 | 0.353466669239891 | 94 |
| chr1:1000, chr1:3000 | 0.39849991080910563 | 94 |
| chr1:1000, chr1:11000 | 0.24053784261608957 | 93 |
| chr1:1000, chr1:250000 | 0.0256751927810228 | 95 |
| chr1:1000, chr2:1000 | 0.008140034754693937 | 95 |
| chr1:4000, chr1:5000, chr1:4000 having no variance | NaN | |

They are compared within 1e-12 relative, and the NaN cells are compared
as NaN. The last pair is on two chromosomes, which this function gives
like any other: it has no distance, and the next item is the one that
leaves it out. The test also asserts that the 68 variants with no
variance have NaN in their row, their column and their diagonal cell, and
that the matrix is the same, to the bit, when the blocks are of 7, 64,
256 and 500 variants and when rayon has 1 thread and 4.

The worked example, the first cargo test, made at `r2_between` of "The
Rust interface" on a block built by hand: 5 variants of 6 individuals,
one of them with a missing genotype and one variant with no variance. It
is `tests/reference/ld/example.vcf`, and plink2 gives the same numbers
for it.

| variant | position | genotypes | dosages |
|---|---|---|---|
| v1 | chr1:1000 | 0/0 0/0 0/1 0/1 1/1 1/1 | 0 0 1 1 2 2 |
| v2 | chr1:2000 | 0/0 0/1 0/1 1/1 1/1 1/1 | 2 1 1 0 0 0 |
| v3 | chr1:3000 | 0/0 0/0 0/0 0/1 ./. 1/1 | 0 0 0 1 none 2 |
| v4 | chr1:4000 | 0/0 0/0 0/0 0/0 0/0 0/0 | 0 0 0 0 0 0 |
| v5 | chr1:5000 | 0/1 1/1 0/0 0/1 1/1 0/0 | 1 2 0 1 2 0 |

The major allele of v1 and of v5 is allele 0, which the tie between the
two alleles gives to the lower numbered one, and of v2 it is allele 1.

| pair | n | Σx | Σy | Σxy | Σxx | Σyy | r² |
|---|---|---|---|---|---|---|---|
| v1, v2 | 6 | 6 | 4 | 1 | 10 | 6 | 27/40 = 0.675 |
| v1, v3 | 5 | 4 | 3 | 5 | 6 | 5 | 169/224 = 0.7544642857142857 |
| v1, v4 | 6 | | | | | | NaN, v4 has no variance |
| v1, v5 | 6 | 6 | 6 | 5 | 10 | 10 | 36/576 = 0.0625 |
| v2, v3 | 5 | 4 | 3 | 0 | 6 | 5 | 144/224 = 0.6428571428571429 |
| v2, v5 | 6 | 4 | 6 | 4 | 6 | 10 | 0 |
| v3, v5 | 5 | 3 | 4 | 1 | 5 | 6 | 49/224 = 0.21875 |

The whole numbers are asserted exactly and the r² within 1e-12 relative.

Half called genotypes and variants of more than two alleles are checked
on `tests/reference/vcf/many.vcf` of `docs/specs/io_vcf.md`, 500 variants
of 50 diploid individuals with 54 variants of more than two alleles and
257 half called genotypes among its 25000. There popnei and plink2 are
not the same calculation, and the measurement that says by how much is
under **Open 2**. So this file checks the dosages and not the r²: the
dosages popnei reads from it, at `LdDosages::dosages` of "The Rust
interface", are compared with pyNei's `to_012`, which counts the alleles
of a half called genotype as popnei does, and the r² of those dosages is
what plink2 verifies on the two files above. A pytest
test also asserts that a variant with no called genotype gives a row of
NaN and not an error.

Against pyNei, a pytest test made at the Python function: a `Variants`
over `ld.vcf.gz` with no missing genotype left in it, which
`variants.filter_by_missing_data(0)` gives, has the same matrix in both
libraries after pyNei's r is squared, within 1e-12 relative. With the
missing genotypes in, the two differ by the numbers of "Missing
genotypes" and the test asserts the difference and not equality, so that
the divergence is pinned and a change to either side shows. A
`max_num_vars` below the variants of the pass is a `ValueError` whose
message has both numbers and the memory the matrix would have needed, and
a pass with no variant is a `ValueError`, as
`docs/specs/dists.md` has it.

The TypeScript test, under node, reads `ld.vcf.gz` from a `Uint8Array`
and asserts the five r² of the table and the `Error` of a `maxNumVars`
of 100.

## LD against distance, per population

### What it gives

How the r² of a pair of variants falls off as the two move apart along a
chromosome, for each population on its own. A population is a named set
of individuals, and the curve of each is read on its own because the
recombination a population has had and the size it has had shape it: a
population that went through few individuals keeps linkage disequilibrium
over longer stretches.

The pairs counted are those whose two variants are on one chromosome and
whose distance, the difference of their positions, is from `min_dist` to
`max_dist`, both included. They are put into `num_bins` bins of equal
width across that range: the width is (max_dist − min_dist + 1) /
num_bins, and a pair at distance d falls in the bin that
floor((d − min_dist) / width) gives, the last bin taking anything the
rounding would put past it. For each population and each bin the result
gives how many pairs it holds, the mean of their r² and the standard
deviation of their r², with the pairs of the bin as the divisor. A pair
with no r², which "What it gives" of the item above defines, is in no
bin.

Before the pairs of a population are counted, the variants whose major
allele frequency among the individuals of that population is above
`max_allowed_maf` are left out of it. The major allele frequency is the
largest of the counts of the alleles over the called alleles, as
`docs/specs/filters.md` defines it for `filter_by_maf`, and it is worked
out over the individuals of the population and not over all of them, so
two populations of one dataset use different variants. The r² of a
variant that hardly varies in a population rests on the one or two
individuals that carry the rare allele, and leaving those variants in
raises the curve everywhere.

### Its Python function

```python
calc_ld_and_dist_per_pop(variants: Variants,
                         pops: dict[str, Sequence[str]] | None = None,
                         min_dist: int = 1,
                         max_dist: int = 1_000_000,
                         num_bins: int = 50,
                         max_allowed_maf: float = 0.95) -> LdAndDistPerPop
```

It is a consumer of the `Variants`, one pass, which serves every
population.

`pops` is a dict of the name of each population to the names of its
individuals, as `docs/glossary.md` has it. `None` is one population with
every individual of the dataset, named `"pop"`, which is what pyNei names
it. A name that no individual of the dataset has is a `ValueError` that
says which, and a population with no individual is a `ValueError` too.
An individual may be in more than one population, and one in none is read
by none of them.

`LdAndDistPerPop` is a frozen dataclass with `per_pop`, a dict of the
name of each population, in the order of `pops`, to a pandas frame with
one row per bin, indexed by the smallest distance of the bin, and the
columns `largest_dist`, `num_pairs`, `mean_r2` and `sd_r2`, with NaN in
the last two for a bin with no pair; `num_vars_per_pop`, a dict of the
name of each population to how many variants passed its major allele
frequency; and `pass_stats`.

It mirrors `calc_ld_and_dist_per_pop` of `pynei/ld.py`. The differences:

- **It gives bins and not a sample of pairs.** pyNei gives, per
  population, a list of at most `max_num_measures_to_keep` pairs, 10000
  by default, drawn from all of them by `more_itertools.sample`, which
  takes no seed, so two runs over one dataset give different points and
  no test can pin them. The bins are computed over every pair, are the
  same whatever the blocks and the threads, and are what a curve of
  linkage disequilibrium against distance is drawn from. Whether the
  result also carries a sample of pairs for a scatter is **Open 1**,
  below.
- **It gives r² and not r** (under the item above).
- **A missing genotype takes its individual out of that pair** (under the
  item above).
- **`min_dist` counts the pair at that distance.** pyNei keeps the pairs
  whose distance is strictly above it, so with its default of 1 a pair of
  variants 1 bp apart is thrown away and nothing says so. Here
  `min_dist` of 1 leaves out only the pairs of two variants at one
  position, a SNP and an indel at the same base, whose distance is 0.
  Decided here.
- **`max_dist` has a default**, 1000000, where pyNei's is `None`, which
  keeps every pair of a chromosome however far apart. A pass with no
  bound on the distance is a pass over every pair of every chromosome,
  which for a chromosome of 100000 variants is 5·10⁹ pairs, and the
  window of "How it runs" would hold the whole chromosome. Decided here.
  The default is plink2's own for `--r2-unphased`, `--ld-window-kb 1000`.
- **`num_bins` is new**, and so is `sd_r2`, which comes with the bins.
- **There is no `method`.** pyNei takes `LDCalcMethod.GENERATOR` or
  `MATRIX` and the two do not count the same pairs: the generator gives
  each pair once and the matrix gives every ordered pair, so each
  unordered pair twice. Run on the 2 chromosomes of 30 variants of
  `test_ld_for_pops_with_filtered_vars`, the generator gives 870
  measures, which is 2 times the 435 pairs of 30 variants, and the matrix
  1740. popnei counts each pair once. Decided here.
- **There is no `max_num_measures_to_keep`**, which goes with the sample.
- **The result has `pass_stats`** and `num_vars_per_pop`, which pyNei's
  has not; pyNei gives back a plain dict of lists.

In TypeScript it is `calcLdAndDistPerPop(variants, {pops, minDist,
maxDist, numBins, maxAllowedMaf})`, where `pops` is an object of
population name to an array of individual names. The result has `perPop`,
an object of population name to `{smallestDist, largestDist, numPairs,
meanR2, sdR2}`, five typed arrays of `numBins` values, the distances and
the counts as `Float64Array` like the positions of a block, `numVarsPerPop`
and `passStats`.

### The cases

A bin with no pair has a `num_pairs` of 0 and NaN for its mean and its
standard deviation. A bin of one pair has that pair's r² as its mean and
0 as its standard deviation.

A dataset of one chromosome whose variants span less than `min_dist`
gives every bin empty, and so does one whose variants are each on a
chromosome of their own. Neither is an error.

A variant with no called allele among the individuals of a population
has no major allele frequency there, so it is left out of that population
as a variant above the threshold is, and it can still be in another
population that called it. A population in which every variant is left
out, by its frequency or for having nothing called, gives every bin empty
and a `num_vars_per_pop` of 0, and the other populations are not
affected.

A pair of variants on two chromosomes has no distance and is in no bin,
as in pyNei.

A pass with no variant is a `ValueError`, as `docs/specs/dists.md` has
it. A `min_dist` above `max_dist`, a `num_bins` of 0, a `max_allowed_maf`
that is not a number from 0 to 1, and a negative `min_dist` are each a
`ValueError` that names the argument and the value.

`test_ld_vs_dist` and `test_ld_for_pops_with_filtered_vars` of pyNei's
`test/test_ld.py` assert that both populations come back and that each
gets the variants of the whole dataset and not those the population
before it had, which was a defect pyNei fixed; popnei takes the second as
a test that two populations of one pass count their own pairs, which the
literals below make sharper, since there the two populations keep
different variants.

### What pyNei does that is odd

Besides the r and the missing genotypes of the item above, and the two
methods that count different numbers of pairs, which is among the
differences.

`calc_ld_and_dist_per_pop` makes one pass of its own for each population,
by building a `Variants` per population with `filter_samples` and
`filter_by_maf` around the one it was given. popnei makes one pass and
counts every population in it, which reads the source once instead of
once per population and cannot give two populations different variants
for any reason but their own allele frequencies.

### How it runs

One pass. The calculation asks its reader for the genotypes and for the
chromosome and the position.

It holds the blocks whose variants are within `max_dist` of the newest
variant read and on its chromosome, and drops a block as soon as every
variant of it is further back than that or on another chromosome. For
each population it holds the three matrices of "How it runs" of the item
above, built over the individuals of that population and over the
variants of the held blocks that passed the major allele frequency of
that population.

The pairs of the window are computed in tiles, as the matrix of the item
above is, six products a tile pair and four on the diagonal, and each
tile pair gives the count, the sum of r² and the sum of the squares of r²
of each bin it touched. The tiles are cut at fixed multiples of the tile
size counted from the first variant of the pass, so they do not move with
the blocks, and the sums of the bins are added up in the order of the
tiles. So the result is the same, to the bit, whatever the size of the
blocks and the number of threads, and the calculation needs no `reblock`
before it, the reader of `docs/specs/block.md` that cuts and joins blocks
to one size.

What is kept from one block to the next is, for each population and each
bin, the count as a `u64` and the two sums as `f64`, and the blocks of
the window. The memory of the window is the genotypes of the variants it
holds and the three matrices of each population over them: for 250
variants within 250000 bp of 1000 individuals and two populations of 500,
250 KB of genotypes and 2 x 24 bytes x 250 x 500, 6 MB. A window that
would hold more variants than the memory can take is the error of
`docs/specs/block.md` for a block the machine has not the memory for.

### How it is verified

Against plink2 v2.0.0-a.7.7 on `tests/reference/ld/ld.vcf.gz`, the
dataset of the item above. The reference script picks, for each
population, the individuals and the variants that pass its major allele
frequency with popnei's rule, which
`docs/specs/filters.md` verifies against bcftools, and runs plink2 on
those alone,

    plink2 --vcf ld.vcf --double-id --allow-extra-chr \
           --keep <pop>.keep --extract <pop>.extract \
           --r2-unphased square bin --out <pop>

so that every r² of every bin comes from plink2 and the binning is the
arithmetic this spec defines. The bins below are of `min_dist` 1,
`max_dist` 250000 and `num_bins` 10, a width of 25000 bp, run on 22
September 2026.

One population with every one of the 100 individuals, `max_allowed_maf`
of 0.95, which 432 of the 500 variants pass:

| the distances of the bin | pairs | mean r² | sd of r² |
|---|---|---|---|
| 1 to 25000 | 8744 | 0.20767885551844031 | 0.20568652974479465 |
| 25001 to 50000 | 7815 | 0.07890359176062511 | 0.08625049214414023 |
| 50001 to 75000 | 6846 | 0.03508441302751711 | 0.0411917416898793 |
| 75001 to 100000 | 5962 | 0.02056917580626069 | 0.026923590343909974 |
| 100001 to 125000 | 5140 | 0.015026451851395499 | 0.02060604430437979 |
| 125001 to 150000 | 4168 | 0.011542104404978385 | 0.015425275975059542 |
| 150001 to 175000 | 3308 | 0.01145438215819949 | 0.015306615529156098 |
| 175001 to 200000 | 2447 | 0.012095572873545887 | 0.016610730726125223 |
| 200001 to 225000 | 1481 | 0.015365254218410632 | 0.020746495170409326 |
| 225001 to 250000 | 530 | 0.013266303346602112 | 0.017768417874071147 |

Two populations, `pop_a` the individuals `i000` to `i049` and `pop_b`
`i050` to `i099`, with `max_allowed_maf` of 0.8, which 396 of the 500
variants pass in `pop_a` and 402 in `pop_b`, 28 of them in `pop_a` and
not in `pop_b` and 34 the other way. At 0.95 both populations pass the
same 432 variants, so the threshold of 0.8 is the one that fails when the
major allele frequency is worked out over all the individuals instead of
over those of the population:

| the distances of the bin | pairs, pop_a | mean r², pop_a | pairs, pop_b | mean r², pop_b |
|---|---|---|---|---|
| 1 to 25000 | 7394 | 0.22226316432228382 | 7625 | 0.21935192592998345 |
| 25001 to 50000 | 6564 | 0.09426862122352 | 6779 | 0.08778756463719926 |
| 50001 to 75000 | 5648 | 0.04565040582359873 | 5968 | 0.0442939962514918 |
| 75001 to 100000 | 4918 | 0.030489796800475328 | 5189 | 0.0321335996857634 |
| 100001 to 125000 | 4304 | 0.024750005446480792 | 4473 | 0.02676089802517462 |
| 125001 to 150000 | 3540 | 0.02220350204444288 | 3567 | 0.02136589829323258 |
| 150001 to 175000 | 2872 | 0.01770136991605654 | 2823 | 0.02578147369672746 |
| 175001 to 200000 | 2137 | 0.02022495044871507 | 2086 | 0.02294629377272581 |
| 200001 to 225000 | 1240 | 0.01792337843122636 | 1275 | 0.022448095913909734 |
| 225001 to 250000 | 438 | 0.020745833685396994 | 415 | 0.016086351215632733 |

The standard deviations of the two populations are left out of the table
to keep it readable and are asserted with the means;
`docs/reports/ld-method/bins.py` prints all three tables with them, and
the implementation plan stores its output beside the reference script.
The counts of pairs are compared exactly and the means and the standard
deviations within 1e-12 relative, since both sides add the same r² and
only the order of the sum can differ. These checks are made at
`calc_ld_and_dist` of "The Rust interface", on the file read with the VCF
reader, and the test runs them with blocks of 7, 64 and 500 variants and
with rayon at 1 thread and 4, which have to give the same numbers to the
bit.

Against pyNei there is no comparison of the values: pyNei gives a sample
of pairs drawn with no seed, of r and not r², with a missing genotype
left in. What is compared is the set of variants each population keeps:
both libraries are asked for the variants that pass a major allele
frequency of 0.8 in each population, popnei through
`num_vars_per_pop` and pyNei through `filter_samples` and
`filter_by_maf`, and the two counts have to be the 396 and the 402 above.

The worked example, a cargo test at `calc_ld_and_dist` on the five
variants of "How it is verified" of the item above, with `min_dist` 1,
`max_dist` 4000 and `num_bins` 2, a width of 2000 bp. The pairs at 1000
bp are v1-v2, v2-v3, v3-v4 and v4-v5, and at 2000 bp v1-v3, v2-v4 and
v3-v5, which is the first bin, 1 to 2000; at 3000 bp v1-v4 and v2-v5 and
at 4000 bp v1-v5, the second bin. Every pair that holds v4, which has no
variance, has no r² and is in no bin. So the first bin holds v1-v2
(0.675), v2-v3 (0.6428571428571429), v1-v3 (0.7544642857142857) and v3-v5
(0.21875), a mean of 0.5727678571428572, and the second holds v2-v5 (0)
and v1-v5 (0.0625), a mean of 0.03125.

The TypeScript test, under node, asserts the ten counts of pairs and the
ten means of the one population of the first table.

## The Rust interface

The dosages of the variants of a block, for the individuals of one
population, held as the three matrices that the products of r² read.

```rust
pub struct LdDosages { /* private */ }

impl LdDosages {
    /// `individuals` are indices into the individuals of the block, in
    /// the order they are given, and an empty slice is every individual
    /// of it. The major allele of each variant is that of
    /// `variant::the_major_allele` over those individuals alone.
    ///
    /// # Errors
    ///
    /// A block that does not pass `check`, one with variants and no
    /// genotypes, an index that is not an individual of the block, and a
    /// block whose variants times its individuals is more than the
    /// linear algebra counts in.
    pub fn of_block(block: &Block, individuals: &[usize]) -> Result<LdDosages>;
    pub fn num_vars(&self) -> usize;
    pub fn num_individuals(&self) -> usize;
    /// The variants `first..first + num_vars` of it, which the tiles of
    /// the products and the window of the filter take.
    pub fn rows(&self, first: usize, num_vars: usize) -> Result<LdDosages>;
    /// Whether the called genotypes of the variant hold two dosages at
    /// least. One that does not has NaN against every variant, itself
    /// among them.
    pub fn has_variance(&self, var: usize) -> bool;
    /// The dosage of each individual at the variant, in the order the
    /// individuals were given, and `None` for a genotype with an allele
    /// missing. It is what the check of `many.vcf` in "How it is
    /// verified" compares with pyNei's `to_012`.
    pub fn dosages(&self, var: usize) -> Option<impl Iterator<Item = Option<u8>> + '_>;
    /// The largest of the counts of the alleles over the called alleles
    /// of the variant, the major allele frequency of
    /// `docs/specs/filters.md`, over the individuals this was built with.
    /// None for a variant with no called allele.
    pub fn maf(&self, var: usize) -> Option<f64>;
}
```

The r² of every variant of one set against every variant of another, the
six products and the formula of "What it gives". `out` is
`a.num_vars()` rows of `b.num_vars()` values, row after row, and a pair
with no r² is NaN.

```rust
/// # Errors
///
/// An `out` that is not `a.num_vars() * b.num_vars()` values, an `a` and
/// a `b` built over a different number of individuals, and what the
/// linear algebra refuses.
pub fn r2_between(a: &LdDosages, b: &LdDosages, out: &mut [f64]) -> Result<()>;
```

The matrix of every pair of the variants of a reader. It asks the reader
for the genotypes, the chromosome and the position, and reads it to its
end. The reader is the outermost reader of the chain that `chain_of` of
`docs/specs/filters.md` built from the steps of the `Variants`, and the
binding crate keeps it, as it does for `calc_kosman_sums` of
`docs/specs/dists.md`: when the calculation returns it reads the counts
of the filters from the chain and the variants of the pass from
`R2Matrix::num_vars`, and the two are the `pass_stats` of the result.

```rust
/// # Errors
///
/// More than `max_num_vars` variants, with both numbers and the memory
/// the matrix would have needed; no variant in the reader; a block with
/// variants and no position; the memory of the matrix, which is asked
/// of the machine with `try_reserve_exact` before the pass and not
/// taken, as `docs/specs/linalg.md` asks for the workspace of the
/// eigendecomposition, so that a matrix this machine cannot hold is an
/// error and not a process that ends; and those of the reader.
pub fn calc_r2_matrix<R: BlockReader + ?Sized>(
    reader: &mut R,
    max_num_vars: usize,
) -> Result<R2Matrix>;

pub struct R2Matrix { /* private */ }

impl R2Matrix {
    pub fn num_vars(&self) -> usize;
    /// `num_vars` x `num_vars`, row after row.
    pub fn r2(&self) -> &[f64];
    /// The interned chromosome of each variant, read through
    /// `chrom_table`.
    pub fn chroms(&self) -> &[u32];
    pub fn chrom_table(&self) -> &ChromTable;
    pub fn poss(&self) -> &[u64];
}
```

The curve of r² against distance, one pass, every population at once.
The populations are the indices of their individuals, in the order the
user gave them, and an empty slice of populations is one population with
every individual.

```rust
pub struct LdAndDistOptions {
    pub min_dist: u64,
    pub max_dist: u64,
    pub num_bins: usize,
    pub max_allowed_maf: f64,
}

/// # Errors
///
/// A `min_dist` above `max_dist`, a `num_bins` of 0, a
/// `max_allowed_maf` that is NaN or not from 0 to 1, a population with
/// no individual, an index that is not an individual, no variant in the
/// reader, a block with variants and no position, and those of the
/// reader.
pub fn calc_ld_and_dist<R: BlockReader + ?Sized>(
    reader: &mut R,
    pops: &[&[usize]],
    options: &LdAndDistOptions,
) -> Result<LdAndDist>;

pub struct LdAndDist { /* private */ }

impl LdAndDist {
    /// The variants the calculation was given, before the major allele
    /// frequency of any population.
    pub fn num_vars(&self) -> u64;
    pub fn num_pops(&self) -> usize;
    /// None when `pop` is not a population of the call.
    pub fn bins_of_pop(&self, pop: usize) -> Option<&LdBins>;
}

pub struct LdBins { /* private */ }

impl LdBins {
    pub fn num_bins(&self) -> usize;
    /// The variants that passed the major allele frequency of this
    /// population.
    pub fn num_vars(&self) -> u64;
    /// The smallest and the largest distance of the bin, both included.
    /// None when `bin` is not a bin.
    pub fn bounds(&self, bin: usize) -> Option<(u64, u64)>;
    pub fn num_pairs(&self, bin: usize) -> Option<u64>;
    /// None when `bin` is not a bin and when it holds no pair, which the
    /// binding crates give the user as NaN.
    pub fn mean_r2(&self, bin: usize) -> Option<f64>;
    /// With the pairs of the bin as the divisor.
    pub fn sd_r2(&self, bin: usize) -> Option<f64>;
}
```

The `variant` module makes one function public, which the PCA work of 22
September 2026 wrote as a private one of `pca.rs` and which this module
and the filter of `docs/specs/filters.md` need, so that the major allele
of a variant is worked out in one place in popnei:

```rust
    /// The most frequent of the alleles `counts` counted, and the lowest
    /// numbered of two that are equally frequent. It is
    /// [`MISSING_ALLELE`] when no allele was called.
    pub fn the_major_allele(counts: &AlleleCounts) -> i8;
```

This module adds these cases to the error of the crate, each of which
both binding crates give their user as the wrong input of a function, a
`ValueError` in Python: more variants than the matrix was allowed; an
argument of the bins that is out of range; and a population with no
individual or an individual that the dataset has not.

## Speed

Measured on 22 September 2026 on the owner's Apple M5 Pro, 18 cores, with
numpy 2.5.3 on Accelerate, which is the BLAS that `docs/specs/linalg.md`
gives the native build, on dosages made at random with 3 in 100 genotypes
missing, the best of 3 runs for the first table and the mean of 5 for
the second.
They are the cost of the products and the arithmetic over them, with no
reader and no binding in them, and they are what popnei has to reach with
a tenth over them, since the products are the work and popnei adds the
building of the three matrices and the tiling.

One block of 5000 variants and 1000 individuals against itself, which is
`calc_r2_matrix` at its default `max_num_vars`:

| | time | what it gives |
|---|---|---|
| one product of the dosages with themselves | 0.050 s | the covariances alone |
| the six products and the r² of every pair | 0.455 s | the 5000 x 5000 matrix, 200 MB |
| pyNei's one product and division | 0.067 s | r, with a missing genotype left in |

So the rule that drops the individuals missing at either variant costs
6.8 times pyNei's, and the number to reach for `calc_r2_matrix` at 5000
variants of 1000 individuals is 0.50 s.

For the curve against distance the cost is set by how many variants fall
inside `max_dist`, since only those pairs are computed. One tile pair of
1000 individuals, and what a pass over 100000 variants of 1000
individuals costs at a window of that many variants, which is two tile
pairs for each tile of variants:

| tile | one tile pair | 100000 variants at that window |
|---|---|---|
| 256 variants | 1.9 ms | 1.5 s |
| 512 variants | 5.7 ms | 2.2 s |

pyNei has no comparable number: `calc_ld_and_dist_per_pop` over a dataset
of that size would hold every chunk in memory and build the whole square
matrix of each chunk pair.

Neither was measured in WebAssembly, and no factor is given for it here,
because the product `docs/specs/linalg.md` timed in both places is not
one of these: it is `A'A` of a 5000 x 1000 block, whose result is 1000 x
1000, where the products above have a result of 5000 x 5000 and five
times the arithmetic. What that spec measured is faer at 187 ms in wasm
against Accelerate at 7.4 ms natively on the threads it takes by itself,
and 10.5 ms on one. The implementation plan measures these products in
both places. The 187 ms is with the vector instructions of WebAssembly,
the ones that work on sixteen bytes at a time, which `.cargo/config.toml`
now passes to both wasm targets and which `docs/objectives.md` made the
floor of the browsers popnei runs in on 22 September 2026, after the
performance review of the Kosman distances asked for them; without them
the same product took 306 ms.

## Open points

The owner decides these. Until then the implementer follows the
"meanwhile" of each.

**Open 1: a sample of individual pairs beside the bins.** The result of
`calc_ld_and_dist_per_pop` gives the mean and the standard deviation of
r² in each bin of distance, over every pair. pyNei gives instead a list
of at most 10000 pairs with their distances, drawn by
`more_itertools.sample`, which takes no seed, so two runs give different
points. The options are the bins alone; the bins and a sample of pairs
drawn so that it does not change with the blocks or the threads, for
which the rule would be to keep the pairs whose hash of the seed and
their two positions is smallest; and the sample alone, as pyNei has it.
What the sample adds is a scatter, from which a user fits a decay curve
of their own and sees the spread of r² inside a bin, which the standard
deviation gives as one number. What it costs is an argument for the seed,
one for how many pairs to keep, a rule that has to be written down for
the numbers to be testable, and a result with two shapes in it.
Recommendation: the bins alone, and the sample added if a user asks for
the scatter. Meanwhile the bins alone are built.

**Open 2: the major allele of a variant with half called genotypes.**
`docs/specs/pca.md` counts the called allele of a half called genotype
when it picks the major allele of a variant, which is what pyNei's
`_count_each_allele` does, and this module reads the dosages by that same
rule. plink2, told to read a half called genotype as missing with
`--vcf-half-call m`, which is the only way it reads one, counts no allele
of it at all. On `tests/reference/vcf/many.vcf`, 500 variants of 50
individuals with 54 variants of more than two alleles and 257 half called
genotypes, the two rules put popnei 0.0132 away from plink2 at the
furthest of 5947 pairs, measured on 22 September 2026; with the half
called genotypes dropped from the counts as plink2 drops them, the two
agree to 2.8e-16. The difference appears only where the two rules pick a
different major allele in a variant of more than two alleles, because for
a variant of two alleles the other choice turns every dosage into the
ploidy minus itself and leaves r² alone. The options are to keep one rule
for the major allele across popnei, which is the one `docs/specs/pca.md`
has and which `docs/specs/filters.md` verifies against bcftools for the
major allele frequency; or to count a half called genotype in the
dosages of this module as plink2 does, which would make plink2 the
reference for this file too and give popnei two rules for the major
allele; or to change `docs/specs/pca.md` and the code built from it.
Recommendation: keep the one rule. Meanwhile the dosages follow
`docs/specs/pca.md` and `many.vcf` checks the dosages against pyNei and
not the r² against plink2.

## Not in this spec

- The filter that thins variants out by their r², pyNei's
  `filter_by_ld_and_maf`: an item of `docs/specs/filters.md`, written
  with this spec, because a filter is a reader over a reader and the
  other filters are there. It uses `LdDosages` and `r2_between` above.
- `iter_rogers_huff_r2` of pyNei, which gives one pair at a time. Its
  reason for being is that the matrix does not fit in memory, and across
  the boundary with Python it would be one call for each pair, 25 million
  for 5000 variants, which section 5 of `docs/architecture.md` rules out.
  What it was used for, the curve against distance, is the second item
  here.
- `calc_ld_along_genome` and `filter_vars_by_ld`, the two names in the
  last lines of `pynei/ld.py`: pyNei has neither, and popnei does not add
  them.
- r² worked out from the haplotype frequencies of two variants, plink2's
  `--r2-phased`, and Lewontin's D′: popnei reads no phase, and pyNei has
  neither.
- The variants of a region of a chromosome, which would make the cap of
  `calc_rogers_huff_r2_matrix` rarely bite: a later item of
  `docs/specs/filters.md`.
- The genomic relationship matrix, which reads the same dosages:
  `docs/specs/kinship.md`, which is not written.
- The read ahead thread of section 3 of `docs/architecture.md`: it is a
  reader over a reader and both calculations here take any reader, so
  nothing here changes with it.
