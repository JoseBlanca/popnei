# The ld module: how strongly two variants go together, and how that falls off with distance

22 September 2026, with the curve fitted to the fall-off and the half
distance added on 24 September 2026. Two variants are in linkage
disequilibrium when the genotype of one tells something about the
genotype of the other, which happens when they sit close enough on a
chromosome that few recombinations have separated them. The `ld` module
measures it, r² for every pair of a set of variants, and it gives how
that measure falls off with the distance between the variants, for each
population on its own: the fall-off in bins of distance, a curve fitted
to the pairs, and the distance at which r² has fallen to half. That is
what a geneticist reads the recombination of a genome and the history of
a population from. The r² of two sets of variants and the matrix of every
pair are built, in `crates/popnei/src/ld.rs` and in the Python and the
TypeScript packages; nothing of the fall-off against distance is. This
spec develops the row `ld` of the table in section 9 of
`docs/architecture.md`.

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

Every one of the six sums is a whole number, and so is every product of
two of them that the formula takes: n·Σxy, Σx·Σy, n·Σxx and (Σx)². The
products and not the sums are what has to be held exactly, since each of
them is the larger, and each is at most N²k² for N individuals of ploidy
k. An `f64` counts whole numbers one by one up to 2⁵³, 9·10¹⁵, so they
are exact while N·k is at most its square root, 94906265: 47453132
diploid individuals, or 372181 at the ploidy of 255 that the VCF reader
takes. For 10000 diploid individuals the largest product is 4·10⁸, seven
orders below the limit. Above N·k the r² loses digits without saying so,
1.3·10⁻¹² relative at a million individuals of ploidy 255, wider than
the tolerance this spec compares within, so `LdDosages::of_block`
refuses such a dataset rather than working it out. Below it the only
rounding in an r² is in the two products of the spreads, the division
and the square at the end.

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
The script that writes it is `tests/reference/ld/make_reference.py`,
beside the file itself, and `tests/reference/ld/run_plink2.sh` runs it
and the plink2 commands again and compares what they give with what is
stored. Every literal of this spec and of the filter item of
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
as NaN. The tolerance is there for a version of plink2 that computes the
expression in another order, not for popnei's own rounding: on 22
September 2026 every one of the 93096 pairs of `ld.vcf.gz` that has an
r², and the 6 of the worked example, came out of the whole numbers with
the bits plink2 has, and the 31654 pairs that have none were NaN on both
sides. So a difference of 1e-13 in this test is something to look at and
not the noise the tolerance allows for. The last pair is on two chromosomes, which this function gives
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

Beside the bins the result carries a curve fitted to the pairs and, read
off it, the distance at which r² has fallen to half. The bins say what r²
is where this dataset happens to have pairs; the curve puts the whole
fall-off into one shape, and the half distance puts that shape into one
number, which is what a plot of a population is labelled with and what
two populations are compared by. "The curve that is fitted" below has the
model, what is made smallest and how.

### The curve that is fitted

The model is the r² expected between two variants of a population under
drift and recombination, Hill and Weir (1988) with the correction for a
finite sample of Weir and Hill (1986). It is what the literature of
linkage disequilibrium decay fits, as Remington et al. (2001) fitted it
to maize. r² is a quantity of one pair of variants, and the model says what
the average of it is expected to be at a given recombination, over the
pairs of a genome and over the histories the population could have had.

With ρ the scaled recombination between the two variants, four times the
effective size of the population times the recombination fraction between
them, and n how many were sampled,

    E[r²] = (10 + ρ) / ((2 + ρ) · (11 + ρ))
            · [1 + ((3 + ρ) · (12 + 12ρ + ρ²)) / (n · (2 + ρ) · (11 + ρ))]

The first factor is the expectation, and it falls from 10/22, 0.4545, at
ρ of 0 towards 0 as ρ grows: two variants that never recombine still do
not reach an r² of 1, because their allele frequencies drift apart. The
second corrects it for r² being measured on a sample and not on the whole
population, and it is what holds the curve up at long distances, where
the r² of a finite sample does not fall to 0.

n is the individuals of the population, and not the individuals times the
ploidy. Hill and Weir wrote n for the gametes sampled, because the r²
they expect is between the alleles that travel together on one
chromosome. popnei's is the estimate named after Rogers and Huff of "What
it gives" above, the correlation across individuals of their dosages, and
what such a correlation is biased upwards by is how many individuals it
was taken over. `docs/reports/ld-method/decay_truth.py` measured it on 24
September 2026 on a simulated population of 200 individuals with two
chromosomes that assort independently: a pair of variants on different
chromosomes carries no linkage disequilibrium, so the mean r² of such a
pair is that bias plus what the population itself holds. Fitting that
mean over samples of 25, 50, 100 and 200 individuals as a constant plus b
divided by the individuals less one gives a b of 0.925, 0.971 and 0.984
in three runs, where the individuals ask for 1 and the individuals times
the ploidy for 0.5. It is the individuals of the population and not the
individuals called at both variants of a pair, which a missing rate of
0.03 puts at about 0.94 of them.

popnei fits one number, the ρ per base pair: 4Nr, four times the
effective size of the population times the recombination per base pair,
which is the r of that product and not the r whose square this module
measures. A pair d base pairs apart has a ρ of d times it. The effective
size and the recombination rate enter only as that product, and one pass
over one dataset does not separate them, so neither is given on its own.
`rho_per_bp` in the results.

Sved's curve, E[r²] = 1/(1 + ρ), is the other one the literature fits,
and it is not the one used. Both estimate the same 4Nr, so the same
simulation is asked which of them gets it back: over its three runs
Sved's lands 69, 122 and 175 per 100 above the 4Nr the population was run
at, where Hill and Weir's lands 12 per 100 below it, 20 above and 31
above. Sved's has to pass through 1 at a distance of 0, so at 200 bp, the
shortest distance those runs hold, it gives 0.92 to 0.94 where the mean
r² measured there is 0.39 to 0.50; Hill and Weir's gives 0.45. What the
simulation does not settle is n: the effective size its runs settled at
is itself known only to about a quarter, 140, 221 and 274 against the 200
individuals they were run with, which is wider than the three answers to
n are apart. That is why n rests on the unlinked pairs above and not on
this.

The fitted value is the ρ per base pair at which the sum, over every pair
the bins counted, of the square of the pair's r² minus the curve at the
pair's distance is smallest. Every pair counts once and at its own
distance, so `num_bins` does not move the fit; `min_dist` and `max_dist`
do, because they choose which pairs there are.

The pairs at one distance are added up before the fit and nothing is
lost. The sum over them of the square of r² minus the curve is the spread
of those pairs around their own mean, which no ρ changes, plus the number
of them times the square of their mean minus the curve. So the fit needs,
for each distance that holds a pair, how many pairs it holds, n_d, and
the sum of their r², S_d, and with m_d for the mean of those r², S_d
divided by n_d, it makes smallest

    Σ over the distances of [ n_d · (m_d − f(d))² ]

where f(d) is the curve at that distance. What that leaves out, the
spread of the pairs of each distance around m_d, is the same at every ρ.

Multiplying the square out gives Σ [ n_d · f(d)² − 2 · S_d · f(d) ],
which is the same number less Σ n_d · m_d². No ρ changes what the two
forms differ by, so they have the same smallest, and the one written
above is the one popnei computes. That difference is 475.36 on the first
population of "How it is verified" below, where the sum above is 23.21 at
the ρ per base pair fitted to its pairs: the terms that depend on ρ are
the same in both forms, and the multiplied-out one adds them inside a
total of −452.15, 19.5 times the 23.21. An `f64` holds about sixteen
digits whatever the size of the number, so the smallest step it can take
from −452.15 is 5.68·10⁻¹⁴ where the smallest step it can take from 23.21
is 3.55·10⁻¹⁵, and the fit tells two ρ apart through a step sixteen times
as coarse. What that costs was measured on 24 September 2026 by fitting
both forms to the same pairs. R 4.6.1, the language "How it is verified"
below checks these curves against, fits the three populations of that
part; the multiplied-out form lands 3.5·10⁻⁸, 5.3·10⁻⁸ and 5.5·10⁻⁹ of
the ρ per base pair away from what R gives for them, where the form above
lands 1.0·10⁻⁹, 9.3·10⁻⁹ and 2.2·10⁻⁹ away. On the table of r² read off
the curve itself that the same part gives the second cargo test of the
fit, whose pairs reach the default `max_dist` of 1000000, the
multiplied-out form lands 1.5·10⁻⁶ from the ρ per base pair the table was
made with, past the 10⁻⁶ these values are compared within, where the form
above lands 3.5·10⁻¹⁰ from it.

Fitting the mean of each bin instead, placed at the middle of the bin,
moves the answer. For the first population of "How it is verified", the
pairs give a half distance of 6810.57 bp, its ten bins give 7886.60, 15.8
per 100 above, and fifty bins give 6541.43, 3.9 per 100 below. The first
bin is both where the curve bends most and where the half distance falls,
so that is where replacing its pairs by one mean costs most.

There is one number to fit, so the smallest is found without derivatives.
The sum is evaluated at 141 values of the ρ per base pair, from 10⁻¹² to
10², spaced by a tenth of a decade; the two neighbours of the smallest of
them bracket a golden section search, the usual one, which holds two
points inside the bracket and drops the end beyond whichever of them has
the larger sum, so the bracket is 0.618 of itself after each step, until
it is narrower than 10⁻⁹ of a decade, 40 steps from a bracket two grid
cells wide. So the fit costs 183 evaluations of the sum above, 141 for
the grid and 2 to open the search and 40 for its steps, each one pass
over the distances that hold a pair. The range covers every fall-off that
can be seen between one base pair and the largest
`max_dist` a user would give: below 10⁻¹² the curve is flat across 10⁶
bp, and above 10² it has fallen before the second base pair. Narrowing
the range, at the bottom to 10⁻⁹, at the top to 10⁰, or both, moves the
fitted value of the first population of "How it is verified" by at most
9.3·10⁻⁹ of itself, which is what the rounding of the sum leaves and not
something the range did.

How close the fit comes was measured on 24 September 2026, on tables it
should get exactly right: r² read off the curve itself at n of 100, at
each of eleven ρ per base pair that are not values of the grid, from
10⁻⁸·³⁷ to 10⁻²·⁷¹, one pair at each of the distances 1000, 2000 and so
on to 250000. It gives the ρ per base pair back to within 2.3·10⁻¹⁰ of
it, the worst of the eleven, and the same eleven tables carried to
1000000 give the same eleven numbers. What holds it there is the bracket:
10⁻⁹ of a decade is 2.3·10⁻⁹ of the ρ per base pair, and the fit lands a
tenth of that from the answer because what it keeps is the best of the
points it evaluated and not an end of the bracket.

On the pairs of a real population the rounding of the sum stops the
search before the bracket does, and that is what the 9.3·10⁻⁹ above is.
The sum of the first population of "How it is verified" is 23.21 at its
smallest, an `f64` steps from 23.21 by 3.55·10⁻¹⁵, and moving the ρ per
base pair away from the smallest does not raise the sum steadily at that
scale: 10⁻⁹ of itself away the sum is 8 of those steps above its
smallest, 9.3·10⁻⁹ away 21, 10⁻⁸ away 15, which is fewer although it is
further, and 3·10⁻⁸ away 41. So the bottom of the sum is a stretch of ρ
some 10⁻⁸ of itself wide in which a step decides, and what the fit
returns is the first point it evaluated at the lowest step, which is
where the grid put its points. Narrowing the range moves the grid;
tightening the bracket moves nothing, and the fit over the whole range
and the fit over each of the three narrowed ranges each give the same
`f64` at brackets of 10⁻¹⁰, 10⁻¹¹, 10⁻¹² and 10⁻¹³ of a decade as at the
10⁻⁹ the search uses. Measured on 24 September 2026 on those 46441 pairs
at 249 distances.

The grid before the search is what looks at the whole range instead of
sliding downhill from a start value into whichever valley holds it. What
is kept is the smallest of every ρ per base pair the fit evaluated, the
141 of the grid included, so a sum with more than one valley inside the
cell the search works on gives at worst the best of the grid, a tenth of
a decade from the smallest.

The half distance is where the fitted curve has fallen to half of its
value at distance 0. The curve falls without turning, so exactly one ρ
gives half of what it gives at ρ of 0, and which ρ that is depends on n
alone: 2.1608135872529166 at n of 100 and 2.2641731329247312 at n of 50.
The half distance is that ρ divided by the fitted ρ per base pair. popnei
solves for it by bisection between ρ of 0 and 10⁶, stopping when the
bracket is narrower than 10⁻¹² of its own middle, about 60 halvings, and
not from a table, so a dataset of another n needs no new number.

What the curve falls towards as ρ grows is 1/n and not 0, so it reaches
half of its value at ρ of 0 only where that half is above 1/n, which it
is from n of 3 upwards. At n of 1 the curve runs from 1.1983471074380165
down to 1 and half of its value at 0 is 0.5991735537190083, and at n of 2
it runs from 0.8264462809917356 down to 0.5 and half is
0.4132231404958678: neither is ever reached, and there the half distance
alone is NaN, the fitted ρ per base pair and the r² at distance 0 being
what the pairs gave. `calc_ld_and_dist` reaches neither n. A population
of one individual has one dosage at every variant, so no variant of it
has variance and it counts no pair; and a population of two gives every
pair it counts an r² of 1, the correlation of two points being 1 or −1
whenever both variants vary, which is above the curve's own ceiling at
every ρ, so its smallest falls at the bottom end of the searched range
and "The cases" gives it the three NaN. Worked out on 24 September 2026,
when `fit_ld_decay` was written.

It is half of the value at distance 0 and not half of the shortest bin.
The value at 0 is the curve's own ceiling, 0.46198347107438015 at n of
100, which the individuals sampled fix on their own, so what is being
halved does not change when `min_dist`, `max_dist` or `num_bins` changes.
`min_dist` and `max_dist` still move the half distance, through the
fitted ρ per base pair; `num_bins` moves neither.

The half distance is read off the fitted curve and not off the pairs, so
it can fall beyond every distance the fit was given, and a user who reads
it as a distance the data reaches to is reading more than there is. The
three individuals `i000`, `i001` and `i002` of
`tests/reference/ld/ld.vcf.gz`, taken as one population at the defaults,
give a half distance of 1868334.8 bp, where `max_dist` is 1000000 and the
furthest pair they counted falls in the bin from 240001 to 260000: the
curve fitted to those pairs is still above half of its value at 0 where
the pairs end, and where it crosses half is worked out from the shape of
the curve alone. Measured on 24 September 2026, on the 24548 pairs those
three individuals counted over the 367 variants they kept.

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
frequency; `decay_per_pop`, a dict of the name of each population to an
`LdDecay`, a frozen dataclass with `rho_per_bp`, `r2_at_zero` and
`half_dist`, the fitted 4Nr per base pair, the fitted curve at distance 0
and the distance in base pairs at which it falls to half of that, the
three of them NaN for a population whose pairs no curve was fitted to;
and `pass_stats`.

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
- **It fits a decay curve and gives the half distance**, which pyNei does
  not: pyNei hands over its sample of pairs and leaves the fitting to the
  user. The curve is fitted to every pair, so two runs over one dataset
  give the same number where two samples drawn without a seed would not.
  The owner chose it on 24 September 2026 over reading the half distance
  off the bins, either at the bin where the mean r² first crosses half of
  the first bin's or by a straight line between the two bins around that
  crossing: those two need no model and no search, and they give an
  answer that moves with `min_dist` and `num_bins` and that a population
  whose curve never crosses does not have at all.
- **It gives r² and not r** (under the item above).
- **A missing genotype takes its individual out of that pair** (under the
  item above).
- **`min_dist` comes before `max_dist`.** pyNei's signature is
  `calc_ld_and_dist_per_pop(variants, pops, max_dist, min_dist,
  max_allowed_maf, method, max_num_measures_to_keep)`, the largest
  distance first, and here the smallest comes first, which is the order
  the bins run in. So a call written for pyNei that gives the distances by
  position asks for something else here.
  `calc_ld_and_dist_per_pop(variants, pops, 500000)` is a `max_dist` of
  500000 in pyNei and a `min_dist` of 500000 here, which on a dataset
  whose variants are all closer than that counts no pair and gives every
  bin empty with nothing said; with both distances by position `min_dist`
  lands above `max_dist` and it is a `ValueError`. Decided here.
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
the counts as `Float64Array` like the positions of a block,
`numVarsPerPop`, `decayPerPop`, an object of population name to
`{rhoPerBp, r2AtZero, halfDist}`, three numbers and not arrays, and
`passStats`.

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

A population whose pairs fall at fewer than two distances has no curve,
and its `rho_per_bp`, `r2_at_zero` and `half_dist` are NaN: one distance
says nothing about a fall-off, whatever a search would return for it. A
population with no pair at all is that case, and so is one left with two
variants, whose one pair is at one distance.

A population whose smallest sum falls at either end of the searched range
of the ρ per base pair has no curve either, and gives the same three NaN.
A curve that is flat across `max_dist`, or that has fallen before the
second base pair, is not a fall-off these pairs pin down, and the number
at the end of the range says where the search stopped and not what the
data says. The bins of such a population are what they would be anyway.

A population of fewer than three individuals has a curve that never falls
to half, for the reason "The curve that is fitted" gives, and its
`half_dist` alone is NaN. No pass reaches it, and a caller of
`fit_ld_decay` with a table of its own does.

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
variant read and on its chromosome. A block every variant of which is
further back than that or on another chromosome is dropped one block
later, when the block after the one that put it out of reach is taken:
the pairs of a block are counted against the window as it stood when that
block arrived, and a variant within `max_dist` of the first variant of a
block can be further than that from the last variant of the same block,
which is what the window reaches back from. Dropping it as soon as it is
out of reach leaves those pairs uncounted, and how many of them there are
depends on the size of the blocks. What the deferral costs is the memory
of the blocks that fell out, held for one step of the pass. For each
population it holds the three matrices of "How it runs" of the item
above, built over the individuals of that population and over the
variants of the held blocks that passed the major allele frequency of
that population.

Every pair is counted at the step of the newer of its two variants. The
pairs of that step are computed in tiles, as the matrix of the item above
is, six products a tile pair and four on the diagonal, and the tiles are
cut at fixed multiples of the tile size counted from the first variant
that population kept in the pass, so they do not move with the blocks.
The r² of one tile of the newest variants against every variant that can
pair with them is taken first, tile pair by tile pair, and the count, the
sum of r² and the sum of the squares of r² of each bin are then added up
one newest variant at a time, the variants that pair with it in the order
of the pass.

The bins are added up in the order of the variants and not in the order
of the tile pairs because a block can end inside a tile: a tile whose
newest variants two blocks gave would add its pairs to a bin in two goes
where one block would add them in one, and a sum of floats moves in its
last bits when the order of its terms does. The order of the variants is
the order of the source, which no block and no tile cuts. So the result
is the same, to the bit, whatever the size of the blocks, the size of the
tiles and the number of threads, and the calculation needs no `reblock`
before it, the reader of `docs/specs/block.md` that cuts and joins blocks
to one size.

What is kept from one block to the next is, for each population and each
bin, the count as a `u64` and the two sums as `f64`. Of the blocks of the
window the window keeps the chromosome and the position of every variant
and nothing else, since that is all it reads to say which blocks are
still in reach. The genotypes are kept by each population, over the
variants that population kept and over its own individuals alone, and
they are what its three matrices are built from each time a block
arrives. A population of 50 individuals of a source of 1000 keeps a
twentieth of the row of a variant, where a copy of the whole row would
cost each population as much as the source has individuals, whatever the
population has. So the memory of the window is those genotypes and the
three matrices of each population over them: for 250 variants within
250000 bp of 1000 individuals and two populations of 500, 250 KB of
genotypes, which is the row of a variant divided between the two
populations, and 2 x 24 bytes x 250 x 500, 6 MB.

A window whose genotypes or whose matrices this machine has not the
memory of is this module's own error, the one every other allocation of
this module is refused with, which names what could not be held, how many
values it holds and how many bytes one of them is. It is not the error of
`docs/specs/block.md` for a block the machine has not the memory for: a
window is not a block, and that error says how many variants a block was
asked to hold and how many individuals and what ploidy the source has,
which describe no window.

For the fit the pass also keeps, for each population and each distance
from `min_dist` to `max_dist`, how many pairs it has counted there as a
`u64` and the sum of their r² as an `f64`, 16 bytes for each distance:
16 MB for each population at the default `min_dist` of 1 and `max_dist`
of 1000000, asked of the machine with `try_reserve_exact` before the pass
and refused rather than taken, as the matrix of the item above is. It is
16 MB for each population and no more is shared between them, so twenty
populations at those defaults ask for 320 MB before the first block is
read. A pair
is added to its distance in the same step that adds it to its bin, so the
order of the variants fixes both sums and the fit is the same to the bit
whatever the size of the blocks and the number of threads.

The fit itself runs when the pass has ended and reads nothing but those
two arrays, over the distances that hold a pair and not over the whole
range. It is 183 evaluations of the sum, and then about 60
halvings for the half distance, which read no data at all. It is not
split across threads.

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

The curve is verified against R 4.6.1, which `docs/objectives.md` names
among the reference programs, on those same three populations and those
same pairs. `docs/reports/ld-method/decay.py` writes, for each of them,
the pairs of plink2's matrix grouped by their exact distance, which for
the first population is 46441 pairs at 249 distances, and
`docs/reports/ld-method/decay.R` fits the curve to them twice: with R's
`optimize`, Brent's method on the sum popnei makes smallest, and with
R's `nls` under its `port` algorithm, which is Gauss and Newton's method
on the residuals and uses the derivatives popnei does not take.
The two agree to 2.1·10⁻⁹ of each other at the furthest of the three
populations, so the literals below are what `optimize` gives, and they
are compared within 10⁻⁶ relative, 480 times the disagreement of two
optimisers and far below anything read off a plot. The check is a cargo
test at `calc_ld_and_dist`, on `ld.vcf.gz` read with the VCF reader, with
the same block sizes and thread counts as the bins above.

| the population | n | ρ per base pair | r² at distance 0 | the half distance, bp |
|---|---|---|---|---|
| every individual, at 0.95 | 100 | 0.00031727347196446889 | 0.46198347107438015 | 6810.5712522189806 |
| pop_a, at 0.8 | 50 | 0.00030068285442483295 | 0.46942148760330576 | 7530.1038938711654 |
| pop_b, at 0.8 | 50 | 0.00031187790821896646 | 0.46942148760330576 | 7259.8060755719744 |

The curve describes this dataset roughly, and the test says that popnei
finds the same smallest as R on the same pairs and not that the model is
right for it. At the middle of the first of the ten bins above the curve
of the first population is 0.1656 where the bin's mean r² is 0.2077, and
at the middle of the last it is 0.0227 where the mean is 0.0133. The
dataset is four founder haplotypes recombined along a chromosome, which
is not the drift and recombination the model is of.

The first cargo test needs neither plink2 nor R, and it is made at
`fit_ld_decay`, which takes the pairs of each distance and not
genotypes: r² is taken from the curve itself at a ρ per base pair of
0.0001 and n of 100, one pair at each of the distances 1000, 2000 and so
on to 250000, and the fit has to give 0.0001 back and a half distance of
21608.135872529165, both within 10⁻⁶ relative. R's `optimize` gives
9.9999999605176671e-05 on that table, 4.0·10⁻⁹ away from the number the
table was made with, so the tolerance is not hiding a wrong answer. The
same test asserts the three NaN of a population left with pairs at one
distance, and of one whose smallest falls at an end of the searched
range.

A second cargo test at `fit_ld_decay` takes a table of the same kind at a
ρ per base pair of 4.17·10⁻⁸ and n of 100, one pair at each of the
distances 1000, 2000 and so on to 1000000, the default `max_dist`, and
asks for 4.17·10⁻⁸ back and for the half distance 2.1608135872529166
divided by it, 5.1818·10⁷ bp, both within the same 10⁻⁶. That curve
falls to half fifty-two times further out than the furthest pair of the
table, so across the whole of it the r² drops from 0.46198 at a distance
of 0 to 0.45294 at 1000000 bp, 2 per 100 of itself, and what tells the
fit where the smallest lies is what is left of its sum after that
shallow a fall-off. It is the case the multiplied-out form of "The curve
that is fitted" gets 1.5·10⁻⁶ wrong, where the form fitted is 3.5·10⁻¹⁰
from the number the table was made with.

A pytest test at `calc_ld_and_dist_per_pop` asserts the three numbers of
the first row of the table above and that a population left with pairs at
one distance has NaN in all three of them, so that the dict and the
dataclass of the Python layer are exercised where the cargo test
exercises the arithmetic.

Sved's curve, which "The curve that is fitted" gives the reasons for not
using, comes out of these pairs too: it leaves 95.73 where Hill and
Weir's leaves 23.21, both of them the pairs of each distance times the
square of their mean r² minus the curve there, which is the sum the fit
makes smallest. It puts the half distance at 1872.81 bp against 6810.57,
the reference dataset agreeing with the simulation on a dataset the model
is not of.

The TypeScript test, under node, asserts the ten counts of pairs and the
ten means of the one population of the first table, and the three numbers
of the first row of the table above.

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
    /// genotypes, an index that is not an individual of the block, an
    /// individual asked for more than once, what the counts of one
    /// variant refuse, which are a variant of more alleles than a count
    /// of them holds and an allele below the missing one, a block whose
    /// individuals times its ploidy are more than 94906265, a matrix
    /// this machine has not the memory for, a block whose genotypes
    /// hold more than 255 alleles each, and a block whose variants
    /// times its individuals is more than the linear algebra counts in.
    pub fn of_block(block: &Block, individuals: &[usize]) -> Result<LdDosages>;
    pub fn num_vars(&self) -> usize;
    pub fn num_individuals(&self) -> usize;
    /// The variants `first..first + num_vars` of it, which the tiles of
    /// the products and the window of the filter take.
    ///
    /// # Errors
    ///
    /// When those are not variants of these dosages.
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
/// reader, a block with variants and no position, the memory of the
/// pairs counted at every distance of every population, which is asked
/// with `try_reserve_exact` before the pass and not taken, the memory of
/// the window, which "How it runs" says is refused the same way, and
/// those of the reader. One more is a defect of popnei and not a wrong
/// input, a `RuntimeError` in Python where the rest are a `ValueError`:
/// a population the pass built no dosages for, which nothing reaches,
/// since a pass that gave no variant is refused before any curve is
/// fitted and a pass that gave one has taken a block into every
/// population.
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
    /// The curve fitted to the pairs of this population, over every pair
    /// and not over the bins above.
    pub fn decay(&self) -> &LdDecay;
}
```

The curve fitted to the pairs counted at each distance, which
`calc_ld_and_dist` calls for each population when its pass has ended and
which the checks of "How it is verified" that need no plink2 are made at.
`dists` are in base pairs and need not be in order, `num_pairs` and
`sum_r2` hold one value for each of them, and `num_individuals` is the n
of "The curve that is fitted". Fewer than two distances,
and a smallest that falls at either end of the searched range, are the
`LdDecay` of three NaN that "The cases" describes and not an error.

```rust
/// # Errors
///
/// Slices of different lengths, a `num_individuals` of 0, a distance that
/// holds no pair, and a sum of r² that is not finite or is below 0.
pub fn fit_ld_decay(
    dists: &[u64],
    num_pairs: &[u64],
    sum_r2: &[f64],
    num_individuals: u64,
) -> Result<LdDecay>;

pub struct LdDecay { /* private */ }

impl LdDecay {
    /// The fitted 4Nr, by how much the scaled recombination ρ grows per
    /// base pair. NaN when no curve was fitted, which "The cases" says
    /// when, and then the other two are NaN as well.
    pub fn rho_per_bp(&self) -> f64;
    /// The fitted curve at a distance of 0, which the individuals of
    /// the population fix on their own.
    pub fn r2_at_zero(&self) -> f64;
    /// The distance in base pairs at which the fitted curve has fallen
    /// to half of `r2_at_zero`. It is read off the curve and not off the
    /// pairs, so it can be beyond every distance the fit was given, as
    /// "The curve that is fitted" shows on three individuals of
    /// `tests/reference/ld/ld.vcf.gz`. NaN when the other two are, and
    /// NaN on its own when the curve never falls to half, which "The
    /// curve that is fitted" says is n of 1 and n of 2 and no other n.
    pub fn half_dist(&self) -> f64;
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

This module adds cases to the error of the crate, which are the ones the
`# Errors` of "The Rust interface" above name, and the two binding
crates divide them as `docs/specs/pca.md` divides its own.

Most are the wrong input of a function, and a `ValueError` in Python:
more variants than the matrix was allowed; an argument of the bins that
is out of range; a table given to `fit_ld_decay` whose three slices are
not of one length, whose individuals are 0, one of whose distances holds
no pair, or one of whose sums of r² is not finite or is below 0, which a
caller with a table of its own reaches and a pass does not, since a pass
gives the three slices compacted together and puts a pair in each; a
population with no individual, an individual that the
dataset has not, or an individual asked for more than once; a variant of
more alleles than a count of them holds, or an allele below the missing
one; a dataset whose individuals times its ploidy are more than
94906265, above which the products of the formula stop being whole
numbers an `f64` holds exactly and the r² loses digits without saying
so; a block whose genotypes hold more than 255 alleles each, since
`dosages` counts a dosage in a `u8` and a dosage is at most the ploidy;
a dataset whose variants times its individuals is more than the linear
algebra counts the values of a matrix in; and a matrix this machine has
not the memory for, which is asked for with `try_reserve_exact` and
refused rather than taken, as `docs/specs/linalg.md` asks for the
workspace of the eigendecomposition. The last four are reachable from a
vars file, whose `popnei` key states a ploidy that the reader takes with
no upper bound, where the VCF reader takes 255 alleles in a genotype at
most.

Five are a `RuntimeError` instead, for the reason `docs/specs/pca.md`
gives for its own two: no argument of any function of this module gives
them, so what a user reads is a defect of popnei and not something they
wrote. They are the variants of a tile or of a window that are not
variants of the dosages they were asked of; two sets of dosages built
over different individuals; a buffer for the r² that does not hold one
value for each pair; a product the linear algebra could not work
out; and a population the pass built no dosages for, which the `# Errors`
of `calc_ld_and_dist` above says nothing reaches.

## Speed

### What popnei takes

Measured on 23 September 2026 by the performance review of
`docs/reports/perf-ld-2026-09-23.md`, on the owner's Apple M5 Pro, 18
cores, 64 GB, macOS, native `aarch64-apple-darwin`, with the products on
Accelerate. The dataset is the 400 MB VCF of `docs/rust_core.md`, 100000
variants of 1000 individuals whose genotypes are missing at a rate of 0.03,
read from the vars file popnei writes of it, already in the page cache. A
number is the median of 5 runs of `crates/popnei/benches/r2_matrix.rs`,
with the best and the worst beside it, taken with nothing else running on
the machine. The memory of the matrix, 200 MB, is asked for inside what is
timed; reading the file is timed apart and given separately.

| the matrix of 5000 variants of 1000 individuals | best | median | worst |
|---|---|---|---|
| the whole call | 0.384 s | 0.386 s | 0.399 s |
| reading the vars file | 0.006 s | 0.006 s | 0.006 s |
| the calculation, less the reading | 0.378 s | 0.380 s | 0.394 s |

**The number to reach is 0.50 s and popnei takes 0.386 s**, which meets it
with 0.114 s to spare. What that comparison is worth is below.

The calculation runs on one core. With `VECLIB_MAXIMUM_THREADS` at 1 it
takes 0.388 s and with it unset 0.386 s, which is the same number: a
product of 1000 variants by 1000 individuals by 1000 variants is small
enough that Accelerate does not split it across cores, and popnei's own
loop over the pairs of tiles is serial. So 17 of the 18 cores of this
machine are idle for the whole of it.

Where the time goes, from a sampling profile of `sample` over the
benchmark, 7107 samples 1 ms apart, at the tile of 256 variants this was
first measured at:

| | share |
|---|---|
| the six products, inside Accelerate | 72.7% |
| the scan of both operands of every product for a value that is not finite | 10.9% |
| the r² of every pair, cell by cell, from the six sums | 5.7% |
| writing each pair of tiles into the matrix and its transpose | 3.8% |
| reading the vars file | 1.8% |
| building the three matrices of the variants | 1.5% |
| zeroing the buffers the products then overwrite | 1.3% |

### What the target compares

The 0.50 s comes from numpy taking 0.455 s for what the table below calls
the six products. It is not the same work. `docs/reports/ld-method/speed.py`
takes six full 5000 x 5000 products, where two of the six are the
transposes of two others and both halves of a symmetric matrix are
computed; popnei takes only the pairs of tiles from the diagonal up, and
four products instead of six on the diagonal. Counting the cells of the
products of each:

| | cells of the products | cells of the r² |
|---|---|---|
| numpy, `docs/reports/ld-method/speed.py` | 150000000 | 25000000 |
| popnei, at the tile of 1000 variants it uses | 76263680 | 13131840 |

So popnei does about half the arithmetic and takes 0.386 s against numpy's
0.455 s, which is 0.848 of numpy's time for 0.509 of its product cells, so
**per unit of arithmetic popnei is about 1.7 times slower than numpy on
Accelerate**, on the same machine and through the same library. On the
calculation alone, 0.380 s against the same 0.455 s, it is 1.6 times.
Meeting 0.50 s says the calculation is fast enough for a user; it does not
say the code is level with numpy. Whether the target should be restated as
about 0.25 s, which is the like-for-like figure, is Open 3 below.

### The tile

How many variants one tile of the products holds was chosen by this
review, over the matrix of 5000 variants of 1000 individuals, best of 5
runs, less the reading:

| tile | 128 | 256 | 512 | 1000 | 1024 | 1250 | 2500 | 5000 |
|---|---|---|---|---|---|---|---|---|
| | 0.539 s | 0.449 s | 0.424 s | 0.382 s | 0.414 s | 0.380 s | 0.403 s | 0.427 s |

Two things pull against each other. A larger tile gives Accelerate a larger
product, which it works out faster per pair of variants, and it calls the
linear algebra crate fewer times, each call scanning both of its operands
for a value that is not finite. Against that, a tile against itself
computes its whole square where the matrix needs half of it, so a larger
tile computes more pairs it throws away: 12819520 of them at 128 and
25000000 at 5000. The bottom is flat from 1000 to 1250.

A tile that does not divide the variants leaves a short last tile whose
products are shaped badly: 1024 takes 0.414 s against the 0.382 s of 1000,
8 per 100 slower for a tile 2 per 100 larger. The tile is 1000, which is
also the better of 1000 and 1250 away from the cap: over 4000 variants
0.255 s against 0.263 s, and over 3000 variants 0.147 s against 0.150 s.

The matrix is the same to the bit at every tile: the tiles cut the variants
and every sum of a pair runs over the individuals, and every sum is a whole
number below the 2^53 an `f64` counts one by one, so no order of
accumulation changes a bit.

### What was measured and not taken

The scan of both operands of every product for a value that is not finite
costs 0.048 s at a tile of 256 and 0.017 s at a tile of 1000, measured by
taking the two scans out of `crates/popnei-linalg` in a build made for the
measurement and putting them back. Nothing on this path needs it: the three
matrices are built from whole numbers and are not changed afterwards.
Hoisting it would need a type in `crates/popnei-linalg` that carries the
promise that the caller has checked, which changes "Errors" of
`docs/specs/linalg.md`. At 4.4 per 100 of the calculation that is not worth
the change, which is the same conclusion, on the same scans, that
`docs/reports/perf-pca-2026-09-22.md` reached on 22 September 2026.

### In a browser

Measured on 23 September 2026 on the same machine, on the wasm package
built with `npm run build` in `js/popnei` and run under node v26.8.2 with
`js/popnei/bench/time_r2_matrix.mjs`, on one thread, with the 128-bit
vector instructions of WebAssembly that `.cargo/config.toml` passes to both
wasm targets. There the linear algebra is faer and not Accelerate. The
median of 5 runs, with the file opened outside every timing and the memory
of the matrix inside it.

| the matrix of 5000 variants of 1000 individuals | best | median | worst |
|---|---|---|---|
| the whole call | 5.618 s | 5.623 s | 5.681 s |
| reading the vars file | 0.006 s | 0.006 s | 0.007 s |
| the calculation, less the reading | 5.611 s | 5.616 s | 5.675 s |

**A browser takes 14.8 times what the native build takes**, 5.616 s against
0.380 s of calculation, and reads the file in the same 0.006 s. There is no
target for a browser and this is the first measurement of one.

The matrix is the same in both: every run printed the 25000000 values
adding to 224076.635484 with none of them not a number, which is the native
benchmark's figure to all six decimals.

The module holds 393 MB of WebAssembly memory at the high-water mark of one
call, which it reaches on the first call and never adds to. The matrix is
200 MB of that, the three matrices of the five tiles 120 MB, the six sums
of a pair of tiles 48 MB and the file 6 MB, which is 374 MB of the 393 MB.
WebAssembly gives no page back to the host, so that mark is what a tab
keeps for its lifetime. It was 593 MB until the binding stopped copying the
matrix on its way out and took it from the core instead, on 23 September
2026.

### For comparison

pyNei's one product and division over the same dosages, which leaves a
missing genotype in rather than dropping the individuals missing at either
variant of a pair, takes 0.067 s with numpy on Accelerate. The rule popnei
follows costs 6.8 times that in numpy. pyNei has no comparable number for
a whole matrix of this size: `calc_ld_and_dist_per_pop` would hold every
chunk in memory and build the whole square matrix of each chunk pair.


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
What the sample adds is points to draw a scatter with, and a sight of the
spread of r² inside a bin, which the standard deviation gives as one
number. Fitting a curve of their own it no longer adds, since "The curve
that is fitted" gives one. What it costs is an argument for the seed,
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

**Open 3: what the 0.50 s of "Speed" should be.** It was set from numpy
taking 0.455 s for six full 5000 x 5000 products, plus a tenth. popnei
takes only the pairs of tiles from the diagonal up, and four products
instead of six on the diagonal, so it does about half that arithmetic:
76263680 cells of products against numpy's 150000000. It takes 0.386 s, so
it meets the target. The options are to keep 0.50 s, which is what a user
waits for and which the calculation now meets; or to restate it as about
0.25 s, the like-for-like figure, against which popnei is 1.7 times slower
than numpy on Accelerate through the same library, which is what would say
whether the code itself is good. What the second costs is that the spec
then states a target popnei does not meet, and no experiment of
`docs/reports/perf-ld-2026-09-23.md` found a way to halve the remaining
time: 72.7 per 100 of it is inside Accelerate. Recommendation: keep 0.50 s
as the target and record the like-for-like figure beside it, which "What
the target compares" of "Speed" does. Meanwhile 0.50 s stands and is met.

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
- The effective size of the population and the recombination per base
  pair on their own. The curve holds them only as their product, the ρ
  per base pair, and one pass over one dataset does not separate them.
  A user who has a genetic map divides by the recombination it gives and
  reads four times the effective size.
- A spread around the half distance. It would come from leaving out one
  resampling group of variants at a time and fitting again, the block
  jackknife of `docs/specs/dists.md`, which needs the pairs counted at
  every distance for each group and one fit for each group left out. What
  the result would carry and what that costs is not worked out here.
- The read ahead thread of section 3 of `docs/architecture.md`: it is a
  reader over a reader and both calculations here take any reader, so
  nothing here changes with it.
