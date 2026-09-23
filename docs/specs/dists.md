# The dists module: distances between individuals and between populations

23 September 2026. The `dists` module tells a user of popnei how far apart
their individuals are, one distance for every pair of them, and how far
apart their populations are, by seven measures of their choice. There is no
code. This spec develops the row `dists` of the table in section 9 of
`docs/architecture.md`, and it has two items: the Kosman distance between
individuals, and the distances between populations, which replaces pyNei's
`calc_jost_dest_pop_dists` and gives four measures pyNei has not. It
depends on `docs/specs/block.md`, which has the block, the run of
consecutive variants held as arrays, and `BlockReader`, the trait of
everything that gives blocks; on `docs/specs/variant.md`, which has the
`Variants` a user holds, with the steps put on it and the counts of a
pass; on `docs/specs/filters.md`, which builds the chain of filters of a
pass; and, for the populations and the counts of one variant over one of
them, on `docs/specs/stats.md`. The `Variants` built from an array of
genotypes of the variant spec was written with this spec and is not needed
by it: the owner decided on 22 September 2026 that the tests run on small
VCF files written for them, and that the array waits for a user who needs
it. Microsatellites reach popnei as VCF files like any other variants, the
owner having settled that on 23 September 2026, so nothing here reads
another format.

## Kosman distances between individuals

### What it gives

For every pair of individuals, the share of their genotypes that differ,
from 0 for two individuals with the same genotype at every variant to 1
for two that have no allele in common at any. It is what a user builds a
tree or a principal coordinate analysis of individuals from.

It is the distance of Kosman and Leonard (2005, Molecular Ecology 14:
415, DOI 10.1111/j.1365-294X.2005.02416.x) for codominant markers, at any
ploidy. At one variant, the two called genotypes are laid side by side,
each allele of one paired with an allele of the other, in the pairing
that leaves the fewest pairs of different alleles, and d is that number
of pairs over the ploidy k, which is the paper's formula 2. The pairing
never has to be searched for: the alleles that both genotypes hold pair
with each other, as many times as both hold them, so with c_ia and c_ja
the number of copies of allele a in the two genotypes,

    d = 1 - (sum over a of min(c_ia, c_ja)) / k

For diploids, which is what pyNei computes, this is

| the two genotypes | example | d |
|---|---|---|
| hold the same alleles | 0/1 and 1/0, 2/2 and 2/2 | 0 |
| share an allele and are not the same | 0/0 and 0/1, 0/1 and 1/2 | 0.5 |
| share no allele | 0/0 and 1/1, 0/1 and 2/3 | 1 |

For haploids it is 0 for the same allele and 1 for different ones, the
paper's simple mismatch. For tetraploids 0/0/0/1 and 0/1/1/1 have d =
1 - 2/4 = 0.5, and 0/0/1/1 and 0/1/0/1 have d = 0; the paper's Table 2
has every pair of biallelic tetraploid genotypes.

A genotype is the alleles it holds, without order, so 0/1 and 1/0 are the
same. Every allele of the variant counts as itself: a multiallelic variant
is not collapsed to the major allele against the rest.

The distance between individuals i and j is the mean of d over the
variants at which both genotypes are called,

    D(i, j) = (sum of d over those variants) / n(i, j)

where n(i, j) is how many such variants there are. It differs from pair to
pair, because each pair has its own missing genotypes.

The sum of minimums is what makes d fast to compute: min(c_ia, c_ja) is
how many of the counts 1, 2, ..., k both copy numbers reach, so the sum
over the alleles is a count of yes or no questions, "do both genotypes
hold at least m copies of a", over the alleles a and the counts m from 1
to k. For diploids the questions at m = 1 count the different alleles
both genotypes hold, `shared`, and those at m = 2 count 1 when both are
homozygous for the same allele, `hom`, so d = 1 - (shared + hom) / 2:

| the two genotypes | shared | hom | d |
|---|---|---|---|
| 0/0 and 1/2 | 0 | 0 | 1 |
| 0/0 and 0/1, 0/1 and 1/2 | 1 | 0 | 0.5 |
| 0/0 and 0/0 | 1 | 1 | 0 |
| 0/1 and 0/1 | 2 | 0 | 0 |

### Its Python function

```python
def calc_pairwise_kosman_dists(variants: Variants,
                               min_num_snps: int | None = None) -> Distances
```

It is a consumer of the `Variants`, as `docs/specs/variant.md` has them:
it makes one pass over the source through the steps the `Variants` has
when it is called, and the `Variants` is as it was afterwards. The
distances are over the variants the filters kept.

`min_num_snps` is how many variants a pair needs to get a distance: a pair
with n(i, j) below it has none. It keeps the name pyNei gives it, although
the variants need not be SNPs. `None` is 0.

`Distances` is the result of every distance calculation of the module. It
holds `dist_vector`, a read only numpy array of float64 with one value for
each pair, NaN for a pair with no distance, in the order (0, 1), (0, 2),
..., (0, N-1), (1, 2), ..., the upper triangle of the square matrix row by
row, which is also the condensed form of `scipy.spatial.distance`; and
`names`, the names of the individuals. `square_dists` gives the N x N
pandas frame, indexed by name on both sides, with 0 on the diagonal.
`triang_list_of_lists` gives the lower triangle of that matrix with its
diagonal as a list of lists, row r with r + 1 values, the form that
Biopython's `DistanceMatrix` takes. It has `pass_stats`, the `PassStats`
of `docs/specs/variant.md` that every result of a consumer has: how many
variants the calculation took, after the steps, and how many each filter
was given and kept. `Distances(dist_vector, names=None, pass_stats=None)`
and `Distances.from_square_dists(frame)` build one from distances that
were calculated elsewhere, with the names 0 to N-1 when none are given
and `None` for the `pass_stats`; the calculation builds its result
through the first.

It mirrors `calc_pairwise_kosman_dists` and `Distances` of
`pynei/dists.py`. The differences:

- There is no `use_approx_embedding_algorithm`. With it pyNei returns the
  Euclidean distances between the vectors of Kosman distances from each
  individual to about log2(N)² of them taken at random: not Kosman
  distances, and not the same from one run to the next, because the
  choice has no seed. pyNei had it because it once compared the pairs one
  by one in Python. The owner decided on 21 September 2026 to leave it out; the
  option not taken was a second item with the choice of those individuals
  and a seed for it. What the exact calculation takes at the largest
  dataset of the objectives is under "Speed".
- There is no `num_threads`. The owner decided on 22 September 2026 that
  no calculation of popnei has it: the threads are those of rayon's pool,
  one for each core unless `RAYON_NUM_THREADS` is set before the first
  call, and popnei never builds that pool. The option not taken was
  `num_threads=None`, with a number making a pool of that size for the
  call. `docs/specs/pca.md` defers to this decision.
- Every ploidy is taken, where pyNei takes 2 alone (under "Missing
  genotypes, pairs with no distance, and what pyNei asserts").
- `names` is a tuple, as the objectives ask for the names of individuals,
  where pyNei has a numpy array.
- A negative `min_num_snps` is a `ValueError`. pyNei takes it, and it
  does what 0 does.
- `Distances` refuses, with a `ValueError`, a vector whose length is not
  N(N-1)/2 for some N. pyNei takes the N below it:
  `Distances([0.1, 0.2, 0.3, 0.4])` has the names and the `square_dists`
  of three individuals, and the 0.4 is in no cell.
- `triang_list_of_lists` is the lower triangle of `square_dists`. pyNei's
  is not, from four individuals on (under "What pyNei does that is odd").
- A pass that gives no variant is a `ValueError`, whether the source has
  none or the steps kept none. Its message says which of the two, and
  for the steps how many variants each filter was given and kept, the
  counts that a failed pass otherwise loses, as "A pass that was not
  finished" of `docs/specs/filters.md` says. It is the one case that
  every calculation over a pass raises, `PassGaveNoVariant` of
  `docs/specs/stats.md`, and the core builds the message: the
  calculation is lent the outermost reader of the chain, so it reads the
  counts from it and puts them in the error, and neither binding crate
  has anything to add. The owner decided that on 22 September 2026, when
  this module and the statistics per population met: the case of this
  module carried no counts and each binding crate wrote the message
  again. pyNei raises a `RuntimeError`,
  which in popnei is the exception of a defect of popnei, and this is an
  input that cannot be calculated on.
- The result does not carry n(i, j), as pyNei's does not. The core result
  has it, so that giving it to the user later changes no Rust.
- The result has `pass_stats`, which pyNei's has not; pyNei keeps the
  counts of its filters in its `Variants`.
- `Distances` holds the array it is given, where pyNei copies it: a user
  who writes into an array they gave `Distances` changes what the result
  holds. A copy would make the 400 MB of 10000 individuals 800 MB while
  the result is built. Decided by the plan on 22 September 2026.

The first three were the owner's decisions, the one before the last
follows from his decision of 21 September 2026 on the consumers, in
`docs/specs/filters.md`, and the other seven are decided here: none
changes a distance that pyNei gives rightly.

In TypeScript it is `calcPairwiseKosmanDists(variants, {minNumSnps})`,
which gives a `Distances` with `distVector`, a `Float64Array` in the same
order with NaN for a pair with no distance, `names`, an array of strings,
`squareDists()`, a `Float64Array` of N x N values row by row, and
`passStats`.

### Missing genotypes, pairs with no distance, and what pyNei asserts

A variant counts for a pair when both genotypes are called, and a half
called genotype, `0/.`, is a missing genotype, as the glossary has it: its
called allele is not compared with anything.

A pair has no distance when n(i, j) is 0, whatever `min_num_snps` is, and
when n(i, j) is below `min_num_snps`, strictly: a pair with exactly
`min_num_snps` variants keeps its distance. Nothing else is affected by
such a pair: the other pairs of the two individuals keep their distances,
and `square_dists` has NaN in its two cells. The diagonal of
`square_dists` is 0 also for an individual with no called genotype.

Every ploidy is taken, where pyNei refuses all but 2 with "Only diploid
are allowed"; the owner decided it on 22 September 2026, and the option
not taken was to refuse them as pyNei does. A haploid or a tetraploid
result has no pyNei to compare with, and is verified against R alone.
One individual gives a `Distances` with an empty vector; no source gives
none, as `docs/specs/block.md` says.

pyNei's tests, in `test/test_dists.py`. `test_kosman_2_indis` asserts, for
two individuals over 11 variants of which 9 are called in both, a sum of
3.0, n = 9 and a distance of 1/3, then 0 for two identical individuals and
0.45 for a third pair. `test_kosman_missing` asserts that two variants
missing in both individuals give the same distance as one variant missing
in each. Both are made at `_KosmanDistCalculator`, a class that compares
one pair in Python and that popnei has nothing like; popnei's tests take
their genotypes and make the same checks at the public function.
`test_kosman_pairwise` asserts the vector 1/3, 0.75, 0.75, 0.5, 0.5, 0 for
four individuals, one of whose genotypes is 1/2.
`test_kosman_pairwise_with_filtered_vars` asserts that the distances of a
`Variants` with a missing data filter that takes nothing out are those of
the `Variants` without it, twice in a row; in popnei the filter is a step
of the `Variants`, and the test also asserts a `pass_stats` with every
variant given and kept. And `test_kosman_dists_with_threads`,
in `test/test_threads.py`, that 2 and 4 threads give the distances of one.

### What pyNei does that is odd

Read and run in pyNei at commit ef0ca6e.

`Distances.triang_list_of_lists` cuts `dist_vector` into runs of 1, 2,
3... values. The vector is in the order of the upper triangle, so the runs
are not the rows of the lower one. For four individuals a, b, c, d with
the vector 1, 2, 3, 12, 13, 23, the row of c comes out as 2, 3 and the row
of d as 12, 13, 23, where the square matrix has 2, 12 for c and 3, 13, 23
for d. With three individuals the two orders happen to agree. Nothing in
pyNei calls it and no test covers it. popnei gives the lower triangle of
the square matrix, `[[0], [1, 0], [2, 12, 0], [3, 13, 23, 0]]` for these
four, which is its test.

`_calc_kosman_dist_sums` adds the distances of a chunk, pyNei's unit of
work, a few thousand variants, in float32, which is exact because no sum
goes above three times the variants of a chunk, 30000. popnei adds integers,
so the values are the same.

The other two, the approximate algorithm and the vector of a wrong length,
are among the differences above.

### How it runs

On the genotypes of each block as sets of bits, and not as the matrix
products of pyNei. For each individual and each block, popnei builds one
set of bits with a bit per variant of the block, `called`, the variants
where the genotype of the individual is called; and for each allele a
that the block holds and each count m from 1 to the ploidy k,
`holds_a_m`, the variants where its called genotype holds m copies of a
or more. For diploids `holds_a_1` is the genotypes that carry a and
`holds_a_2` the a/a ones. For a pair, with `count` the number of bits
that are set in both of two sets:

    n       = count(called_i, called_j)
    matches = sum over a and m of count(holds_a_m_i, holds_a_m_j)
    k times the sum of d = k * n - matches

which is the formula of "What it gives" added over the variants of the
block. A count is an AND and a count of the ones of each 64 bit word, and
everything is an integer. A block has 1 + k * A sets per individual, with
A the alleles it holds: 5 for a biallelic diploid block, 9 for a
biallelic tetraploid one, and the pairs cost in proportion.

What is kept from one block to the next is two `u32` for each pair, k
times the sum of d and n, to which each block adds. The division is made
once, at the end. So the result is the same to the last bit whatever the
size of the blocks and the number of threads, and the calculation needs
no `reblock` before it, the reader of `docs/specs/block.md` that cuts and
joins the blocks of its source to one size. A sum that would go above
`u32::MAX`, which takes more than 4000 million over k variants, two
thousand million for diploids, is an error and does not wrap.

Natively the pairs of a block are spread over the threads of rayon; in
wasm they are walked on one thread. The calculation asks its reader for
the genotypes alone.

The memory grows with the square of the individuals: 8 bytes for each
pair, which for 10000 individuals, 5e7 pairs, is 400 MB, and the vector of
float64 that the user gets is another 400 MB, and `square_dists` 800 MB
more when it is asked for. pyNei keeps two float64 for each pair, 800 MB.
The sets of bits of a block are small beside that: 3 MB for 5 sets of
5000 variants and 1000 individuals, and the same for 500 variants and
10000 individuals.

The owner decided on the bits on 21 September 2026, and the option not
taken was pyNei's: matrices of 0 and 1 for `called`, `carries_a` and
`hom_a`, variants x individuals, in float32, and their products, which
give the same counts for every pair at once. In the table the products
are made by faer, the linear algebra library in pure Rust that the
architecture gives the wasm builds, and, in pyNei, by Accelerate, the BLAS
of macOS that numpy calls. The measurement behind the
decision, on one block of 5 million genotypes made at random, 3 in 100 of
them missing, on the owner's Apple M5 Pro, one thread unless said, best of
3 runs, load average 2 to 3:

| one block | bits, Rust | products, faer 0.22.6 | pyNei, numpy 2.5.3 on Accelerate |
|---|---|---|---|
| 5000 variants x 1000 individuals, 2 alleles | 0.044 s | 0.51 s | 0.030 s |
| the same on 18 cores | 0.017 s | | |
| the same in wasm, node 26.8 | 0.065 s | 2.97 s | |
| 5000 x 1000, 4 alleles | 0.067 s | 0.87 s | 0.055 s |
| 500 x 10000, 2 alleles | 0.349 s | 4.63 s | 0.457 s |
| the same on 18 cores | 0.090 s | | |

The trial was of diploids, with the two sets `carries_a` and `hom_a` for
each allele, which are the `holds_a_1` and `holds_a_2` above. The three
gave the same two integers for every pair. Accelerate runs the
products on the matrix units of the Apple chip. On the first block the bits are 1.5 times slower than Accelerate, 11 times
faster than faer natively and 45 times faster than faer in wasm, where
there is no BLAS, and on the third block they are faster than Accelerate; they do not need
the `linalg` module, which the architecture puts after `dists` and which
is not written; and their sums are integers. Not measured: OpenBLAS on
x86, and faer built with the SIMD instructions of wasm. The programs and
how to run them are in `docs/reports/kosman-method/`.

How the sets are laid out in memory, whether a block with only the
alleles 0 and 1 gets fewer sets, and how the building of the sets, 0.015 s
of the 0.044 s and serial in the trial, is spread over the threads, are
for the implementer to choose, as long as the times of "Speed" are met.

### How it is verified

Against `gd.kosman` of the R package PopGenReport 3.1.3, which gives the
Kosman distance of every pair of a dataset, and the number of variants it
used for each, at any ploidy. It computes d as the sum, over the alleles,
of the absolute difference between how many copies of the allele the two
genotypes hold, over twice the ploidy, which is the d of "What it gives":
the copies that do not pair are the ones counted twice in that sum, once
on each side. It was run with R 4.6.1 and adegenet 2.1.11. PopGenReport itself does not install on the owner's machine,
because its dependency terra needs the GDAL library, so the function was
read from `R/gd.kosman.r` of the source package with `source()`; it uses
adegenet and base R alone. Two functions that look like it are not
references: `diss.dist` of poppr 2.9.8 divides by all the variants and not
by those called in both, and was up to 0.029 away from pyNei on the panel
below, and `dist.codom` of mmod 1.3.3 drops every individual that has a
missing genotype.

The first dataset is the panel of pyNei's
`test/gwas_reference/sim_missing.vars`: 200 individuals, `s000` to
`s199`, 1200 biallelic variants, 3 in 100 genotypes missing whole. It is
a vars file of pyNei, which popnei does not read, so the reference script
reads it with pyNei and writes it as a VCF of 1 MB. pyNei's own
`work/sim_missing.vcf` is not tracked by git, and the script that writes
it needs plink2 and GMMAT.
Over its 19900 pairs the largest absolute difference between `gd.kosman`
and pyNei is 2.8e-17, one unit in the last place of a float64 of that
size, in values that went through a text file with 17 digits. The second
has more alleles: 300 variants of 40 individuals named `i00` to `i39`,
from `numpy.random.default_rng(3)`, alleles `integers(0, 4, size=(300, 40,
2))` and then the genotypes where `random((300, 40)) < 0.05` set to
missing. Over its 780 pairs the difference is 0. Two more are of other
ploidies, 200 variants of 12 individuals each, three alleles, 5 in 100
genotypes missing: a tetraploid one, `t00` to `t11`, from
`default_rng(5)`, and a haploid one, `h00` to `h11`, from
`default_rng(6)`, made the same way with the ploidy as the third axis.
Over their 66 pairs each, `gd.kosman` and a loop in Python over the
genotypes that pairs the alleles as "What it gives" says differ by 5.6e-17
and 1.1e-16, one unit in the last place; there pyNei is not run.

The literals, each pair with the n that `gd.kosman` reports and k times
the sum of d, which is its distance times k n and came out a whole number
within 1e-10 for every pair, and which a loop in Python over the genotypes
of these pairs gave too:

| dataset | pair | k times the sum of d | n | distance |
|---|---|---|---|---|
| panel | s000, s001 | 372 | 1122 | 0.1657754010695187 |
| panel | s000, s002 | 376 | 1128 | 0.16666666666666666 |
| panel | s198, s199, the last | 351 | 1134 | 0.15476190476190477 |
| panel | s010, s033, the largest | 804 | 1123 | 0.3579697239536955 |
| panel | s116, s119, the smallest | 310 | 1133 | 0.13680494263018536 |
| 4 alleles | i00, i01 | 325 | 269 | 0.6040892193 |
| 4 alleles | i00, i02 | 304 | 264 | 0.5757575758 |
| 4 alleles | i00, i03 | 347 | 272 | 0.6378676471 |
| tetraploid | t00, t01 | 282 | 188 | 0.375 |
| tetraploid | t00, t02 | 284 | 183 | 0.38797814207650272 |
| tetraploid | t00, t03 | 294 | 183 | 0.40163934426229508 |
| haploid | h00, h01 | 112 | 180 | 0.62222222222222223 |
| haploid | h00, h02 | 123 | 184 | 0.66847826086956519 |
| haploid | h00, h03 | 113 | 179 | 0.63128491620111726 |

The integers are compared exactly, and the distances within 1e-9, the
digits of the shortest ones above. These checks are made at
`calc_kosman_sums` of "The Rust interface", on the four files read with
the VCF reader. The scripts that gave these numbers, `ref_export.py`,
`ref.R`, `literals.py`, `poly_export.py` and `poly.R`, are in
`docs/reports/kosman-method/`; the plan turns them into the reference
script of `tests/reference/dists/`, with the four VCF files, gzipped, and
the output of R beside it. The TypeScript
function is tested under node against the literals of the panel.

No program outside the project checks the half called genotype: the
genotypes reached R as whole genotypes or as missing.

Against pyNei: both libraries run `calc_pairwise_kosman_dists` on the two
files and on the genotypes of the worked example below, the last written
as a VCF for the test, with `min_num_snps` of `None` and of a value
that leaves some pairs without a distance, 1125 for the panel. The vectors
have to be equal exactly, with NaN in the same places: both sides divide
the same two integers once.

That the size of the blocks and the number of threads change nothing is a
cargo test at `calc_kosman_sums`: the 4 allele dataset in blocks of 7, of
64, of 65 and of 300 variants, and in rayon pools of 1 and of 4 threads,
gives the same integers.

The worked example, the first cargo test, at `calc_kosman_sums` and
`KosmanSums::dist`: 4 variants, 3 individuals.

| variant | s0 | s1 | s2 | d(s0, s1) | d(s0, s2) | d(s1, s2) |
|---|---|---|---|---|---|---|
| 1 | 0/0 | 0/1 | 1/1 | 0.5 | 1 | 0.5 |
| 2 | 0/1 | 0/1 | 1/2 | 0 | 0.5 | 0.5 |
| 3 | 0/0 | 0/. | 2/2 | | 1 | |
| 4 | ./. | 1/1 | 1/1 | | | 0 |
| 2 times the sum of d | | | | 1 | 5 | 2 |
| n | | | | 2 | 3 | 3 |
| distance | | | | 0.25 | 0.833333 | 0.333333 |

pyNei gives this vector, and `gd.kosman` too, with the `0/.` written as
missing. With `min_num_snps=3` the first pair has no distance and the
other two keep theirs, and with 4 no pair has one.

A tetraploid and a haploid worked example, checked at the same two
functions, with the numbers that `gd.kosman` gives for them:

| variant | t0 | t1 | t2 | d(t0, t1) | d(t0, t2) | d(t1, t2) |
|---|---|---|---|---|---|---|
| 1 | 0/0/0/1 | 0/1/1/1 | 1/1/1/1 | 0.5 | 0.75 | 0.25 |
| 2 | 0/0/1/1 | 0/1/0/1 | 0/0/2/2 | 0 | 0.5 | 0.5 |
| 3 | 0/0/0/0 | 0/0/./0 | 1/2/2/2 | | 1 | |
| 4 times the sum of d | | | | 2 | 9 | 3 |
| n | | | | 2 | 3 | 2 |
| distance | | | | 0.25 | 0.75 | 0.375 |

| variant | h0 | h1 | h2 | d(h0, h1) | d(h0, h2) | d(h1, h2) |
|---|---|---|---|---|---|---|
| 1 | 0 | 0 | 1 | 0 | 1 | 1 |
| 2 | 0 | 1 | 2 | 1 | 1 | 1 |
| 3 | . | 1 | 1 | | | 0 |
| 4 | 0 | 0 | 0 | 0 | 0 | 0 |
| the sum of d | | | | 1 | 2 | 2 |
| n | | | | 3 | 3 | 4 |
| distance | | | | 0.333333 | 0.666667 | 0.5 |

## Distances between populations

### What it gives

For every pair of populations, how far apart they are, as seven numbers a
user chooses among, each with a standard error: Hudson's F_ST, f_2, the
chord distance, Nei's D_A, Jost's D, Nei's G_ST and Hedrick's G'_ST. They
come from five calculations, since D_A is the square of the chord distance
and G'_ST is G_ST divided by the largest value it could take. What a user
does with them is a matrix, a tree or a principal coordinate analysis of
populations, as the Kosman distances above serve for individuals.

The seven are calculated in one pass, because all of them are functions of
the same two counts: how often each allele was called in each population
at each variant, and how many genotypes of the population were called
there. So a user who wants to compare F_ST with Jost's D pays for one
reading of the variants and not two. Which measure answers which question
is in each item below.

Every allele of a variant counts as itself. A multiallelic variant is not
collapsed to the major allele against the rest, which is what lets the
same pass serve microsatellites, whose loci have many alleles, and SNPs,
which have two. A user with microsatellites reads a VCF whose records
carry one allele per repeat length, as any other VCF.

The frequency of allele a in population P at one variant is

    p_Pa = c_Pa / n_P

where c_Pa is how many copies of a were called in the individuals of P at
that variant and n_P is the sum of c_Pa over the alleles, the called
alleles of P there. A half called genotype gives its called allele to
c_Pa and to n_P, as it does in the expected heterozygosity of
`docs/specs/stats.md`.

Two quantities carry most of the arithmetic. Between the populations A
and B,

    H_b = 1 - sum over a of (p_Aa * p_Ba)

is the chance that one allele drawn from A and one drawn from B are not
the same allele. Within them,

    H_w = (u_A + u_B) / 2,  where  u_P = (n_P / (n_P - 1)) * (1 - sum over a of p_Pa^2)

is the mean of the two within population heterozygosities, each corrected
for being estimated from the same copies it is computed over, which is
Nei's 1978 correction written on called alleles rather than on genotypes.
H_b and H_w are the between and the within of Hudson's F_ST, and their
difference is f_2.

### Its Python function

```python
def calc_pop_dists(variants: Variants,
                   pops: dict[str, Sequence[str]],
                   measures: Sequence[PopDistMeasure] | None = None,
                   min_num_samples: int = 20,
                   jackknife_group: int | Literal["variant"] | None = 5_000_000,
                   ) -> PopDists
```

It is a consumer of the `Variants`, as `docs/specs/variant.md` has them:
one pass over the source through the steps the `Variants` has when it is
called, and the `Variants` is as it was afterwards.

`pops` is a dict of population name to the names of its individuals, the
`Pops` of `docs/specs/stats.md`, which refuses a name that is not an
individual of the source and a population with no individual. Fewer than
two populations is an error. An individual in two populations is allowed,
and an individual in none takes no part.

`measures` is which of the seven to calculate. `None`, the default, is all
of them, since the pass is what costs and each measure is a division at
the end of it. `PopDistMeasure` is a `StrEnum` with `fst`, `f2`, `chord`,
`da`, `dest`, `gst` and `gst_hedrick`, so a caller writes `"fst"`.

`min_num_samples` is how many called genotypes a population needs at a
variant for that variant to count for a pair: a variant counts for the
pair (A, B) when A and B each have at least that many, and it counts for
no other pair of the same run when one of them has fewer. The default is
20, inherited from pyNei's `MIN_NUM_SAMPLES_FOR_POP_STAT`, and nobody has
measured whether 20 is the right threshold. Each pair therefore has its
own count of variants, which the result carries.

`jackknife_group` is how the variants are cut into the groups that the
standard errors are resampled over, "The standard errors" below. A number
is a length in base pairs of a chromosome, `"variant"` makes each variant
its own group, and `None` asks for no standard errors, which is the only
value that does not need the chromosome and the position of the variants.

`PopDists` is a frozen dataclass with one `Distances` for each measure
that was asked for, under the name of the measure, `fst`, `f2`, `chord`,
`da`, `dest`, `gst` and `gst_hedrick`, and `None` for the ones that were
not; `pops`, the names of the populations in the order of the `pops`
dict, which is the order of the pairs of every `Distances` in it;
`num_vars`, a read only numpy array of int64 with, for each pair in that
order, how many variants counted for it; `f2_groups`, a read only float64
array of groups x pairs holding f_2 within each group, or `None` when no
standard errors were asked for, from which f_3 and f_4 can be built later
without reading the genotypes again; `group_ids`, a tuple of one
`(chrom, start, end)` for each group, the chromosome and the first and
last position of its variants, both included; and `pass_stats`, the
`PassStats` of `docs/specs/variant.md` that every result of a consumer
has.

`Distances` gains one field for this item, `standard_errors`, a read only
float64 array with one value for each pair in the order of
`dist_vector`, NaN where there is none, and `None` when the calculation
did not give any. `Distances.square_standard_errors()` gives it as the
N x N pandas frame with NaN on the diagonal. The Kosman distances above
leave the field `None` and nothing of theirs changes.

It mirrors `calc_jost_dest_pop_dists` of `pynei/dists.py`, which
calculates Jost's D and nothing else. The objectives ask for every
difference from pyNei to be written down. Four are decided here:

- The name. pyNei has one function per distance and popnei has one
  function for all of them, because the pass is the cost and the counts
  are shared. The owner decided this on 23 September 2026; the option not
  taken was one function for each of them with pyNei's naming,
  `calc_jost_dest_pop_dists` among them, each making a pass of its own.
- The order of the populations. pyNei sorts their names,
  `sorted(pop_idxs.keys())` in `calc_jost_dest_pop_dists`. popnei keeps
  the order of the `pops` dict, as every statistic of
  `docs/specs/stats.md` does. No value changes.
- `num_threads` is not an argument. No calculation of popnei has it,
  which the owner decided on 22 September 2026 in the Kosman item above.
- There is no `alleles` argument. pyNei takes the alleles to count so that
  the columns of different chunks line up; popnei counts the alleles each
  variant has, and nothing is lined up across blocks.

In TypeScript, `calcPopDists(variants, pops, options)`, with `pops` an
object of population name to the names of its individuals and the same
defaults. Each measure is a `Float64Array` in the order of the distance
vector with its `standardErrors` beside it, `numVars` an `Int32Array`,
`f2Blocks` a `Float64Array` of groups x pairs with its two lengths given,
and `blockIds` an array of `{chrom, start, end}`.

### Variants that do not count, populations with little data, and negative values

A variant counts for a pair only when both of its populations have at
least `min_num_samples` called genotypes at it. A variant that counts for
one pair and not for another is the usual case, and it is why each pair
carries its own `num_vars` and why the measures of two pairs are means
over different variants.

A pair with no variant at all has no value for any measure, NaN in the
distance vector, and `num_vars` 0 for it. The pass is not an error: other
pairs may have values.

A population with no called genotype at a variant has fewer than
`min_num_samples` of them for every positive threshold, so it drops that
variant for every pair it is in, and no frequency of it is ever divided
by a zero.

f_2 and F_ST can come out negative, for one variant and for a whole
dataset. It is the correction doing its work: two populations that differ
by no more than the sampling noise have a H_b no larger than their H_w,
and an estimator that could not go below zero would be biased upwards.
Variant 4 of the worked example below has f_2 = -0.1. popnei does not
clamp them, and a user who sees a small negative F_ST has populations
that this dataset cannot tell apart.

A `min_num_samples` of 1 leaves two cases that a larger threshold does
not. A population with one called genotype has n_P = 2 called alleles, so
u_P is defined; and the harmonic mean of the called genotypes that Jost's
D corrects with, "Jost's D" below, is 1, and its factor divides by zero,
so that pair has no value for D, G_ST and G'_ST at that variant while it
keeps one for F_ST, f_2 and the chord distance.

### How it runs

One pass, over the rows of each block with rayon across them, and it
needs no `reblock` before it. For each variant of the block and each
population, the allele counts over the indices of that population come
from `count_alleles_of` of the `variant` module, the row helper that
`docs/specs/stats.md` already uses, and the called genotypes from the
same row. Everything else is arithmetic on those counts.

What is kept from one block to the next, and what the threads reduce
into, is, for each pair of populations and each resampling group, five
f64 and one u32:

    sum of H_b, sum of H_w, sum of sqrt(p_Aa * p_Ba) over the alleles,
    sum of the corrected H_S, sum of the corrected H_T, and the variants
    that counted

Every measure is a ratio of those six, so the divisions happen once, at
the end, and the result does not depend on where the block boundaries
fell nor on how many threads ran. The memory is 44 bytes for each pair
and each group: 66 KB for 3 populations and 500 groups, and 27 MB for 50
populations, which is 1225 pairs, and the same 500. It does not grow with
the variants or with the individuals.

The cost of the pass grows with the square of the populations: per
variant it is one pass over the genotypes for the counts, P of them for P
populations, and then P(P-1)/2 pairs times the alleles of the variant for
the sums. For 20 populations of a 6 allele microsatellite that is 190
pairs times 6, about 1100 multiplications a variant against 190 counting
passes, and for 3 biallelic populations it is 3 pairs times 2. No
measurement of either has been made; "Speed" below says what has to be
measured before the plan claims a number.

The pass asks its reader for the genotypes alone when `jackknife_group`
is `None` or `"variant"`, and for the genotypes with the chromosome and
the position when it is a length, since that is what cuts the groups.

### The standard errors

A distance between two populations is a mean over variants, and variants
near each other on a chromosome carry much the same history, so treating
each of them as an independent draw makes the error look smaller than it
is. The block jackknife of the literature takes that into account, and popnei
calls its unit a resampling group, because a block in popnei is the run of
variants a reader gives, section 2 of `docs/architecture.md`, and the two
have nothing to do with each other: the variants are cut into
groups long enough that two groups are nearly independent, each group is
left out in turn, the measure is calculated again from the sums of the
others, and how much the answer moves is what the standard error is built
from.

The groups are stretches of one chromosome `jackknife_group` base pairs
long. Walking the variants in the order the reader gives them, a new group
starts at the first variant of a chromosome and at the first variant whose
position is `jackknife_group` or more beyond the first variant of the group
being filled. So a group is anchored on its own first variant and not on a
grid of multiples: cutting at multiples of the length instead leaves a
group of one variant wherever a variant sits just past a multiple, and one
such group is enough to move the standard error, 0.00242 against 0.00180
for p0 and p1 of the biallelic panel below. It is also the rule ADMIXTOOLS
2 uses, which is what lets the numbers be compared. With `"variant"` each
variant is its own group, which is what a dataset of a few hundred
microsatellite loci scattered over a genome wants, and which a dataset of
linked SNPs must not use.

The estimator is the delete-m jackknife for unequal m of Busing, Meijer
and van der Leeden (1999, Statistics and Computing 9: 3, DOI
10.1023/A:1008800423698), which is what the f-statistics literature uses,
and which allows groups with different numbers of variants. Write t for
the measure over all the variants of a pair, n for how many variants
counted for that pair, g for the groups that hold at least one of them,
m_j for how many of them are in group j, and t_(j) for the measure
calculated with group j left out, which is a ratio of the same sums minus
that group's. With h_j = n / m_j, the pseudo-value of group j, the
jackknife estimate and its variance are

    u_j = h_j * t - (h_j - 1) * t_(j)
    t_J = sum over j of u_j / h_j
    v   = (1 / g) * sum over j of (u_j - t_J)^2 / (h_j - 1)

and the standard error is the square root of v. What the result gives as
the distance is t, the measure over all the variants, and not t_J: the
jackknife estimate is there to build the variance from. A pair whose
variants all fall in one group, g = 1, has no standard error, and so does
a pair with no variant.

A group with one variant, m_j = 1, has h_j = n and its pseudo-value is
defined; h_j - 1 is zero only when a group holds every variant of the
pair, which is the g = 1 case. So the formula never divides by zero.

How many groups are enough is not something popnei can know, and a
standard error from a handful of them is not one a user should quote. The
result carries the number of groups, and the function raises when it is
below 20 (**Open 1**, below).

It is verified against ADMIXTOOLS 2.0.10, which computes f_2 with this same
jackknife, in the run of the f_2 item below: with the 12 groups that
`blgsize = 100000` cuts the biallelic panel into, its standard errors for
the three pairs are 0.001803417, 0.001591624 and 0.002272635, and the
formula above on the same groups is 4.1e-10 from the furthest of them,
which is the last of the nine digits compared. The tests compare within
1e-8 absolute. No other measure has its standard error checked outside the
project, because no program prints one for them; what is checked is the
arithmetic they share, which is this.

### How it is verified

Four reference programs, run on 23 September 2026 on the owner's machine,
and pyNei for Jost's D. Which program checks which measure is in each item.
`tests/reference/pop_dists/make_reference.py` writes the multiallelic panel
and runs every program on both panels, keeping each output beside itself, as
the other reference data of `tests/reference/` is kept. It refuses a version
of plink2, adegenet, mmod or admixtools other than the ones the numbers
below were taken with.

The first panel is the biallelic one the Kosman item and
`docs/specs/stats.md` already use, `tests/reference/dists/panel.vcf.gz`:
1200 variants of 200 diploid individuals, `s000` to `s199`, 3 in 100
genotypes missing whole, with the three populations of
`tests/reference/stats/panel_pops.txt`, p0 with 48 individuals, p1 with
68 and p2 with 84.

The second is multiallelic, written for this item, which no panel of
popnei was: `tests/reference/pop_dists/micro.vcf.gz`, with the populations
of `micro_pops.txt` beside it. It holds 120 loci of 90 diploid individuals
in three populations of 30, alleles written as repeat lengths so that a
record looks like the microsatellite it stands for, six alleles to a locus
of which 100 loci have all six and 20 have four or five, 4 in 100 genotypes
missing whole. Its allele frequencies are drawn per population from a
Dirichlet around a per locus ancestral one, so the populations differ by
drift and not by a pattern chosen by hand; the seed is 7, and the file is
kept in git because the literals below are of these genotypes and not of
another draw. It is the panel that shows the arithmetic does not assume two
alleles.

The numbers of each measure on both panels go into the tests as literals,
and the items below give them. Two checks are common to all of them.

The measures do not change with the size of the blocks or with the number
of threads, which no test of pyNei checks for its own: the sums are
reduced in a different order and the ratios are taken once, so the test
compares the same panel read in blocks of 100 and of 10000, on 1 and on 4
threads, within 1e-12 relative.

Against pyNei, which has one of the seven: both libraries run on
`tests/reference/dists/panel.vcf.gz` with the same three populations and
`min_num_samples=20`, popnei's `dest` against pyNei's
`calc_jost_dest_pop_dists`, within 1e-12 relative, since the two compute
the same estimator with the sums added in different orders.

The worked example, which becomes the first cargo tests: 4 variants, 6
diploid individuals, pop1 = i0, i1, i2 and pop2 = i3, i4, i5,
`min_num_samples` 1, `jackknife_group` `"variant"`.

| variant | pop1 | pop2 | what it is for |
|---|---|---|---|
| 1 | 0/0 0/1 0/0 | 1/1 0/1 1/1 | biallelic, the populations differ |
| 2 | 0/1 1/2 2/2 | 0/0 0/1 0/0 | three alleles |
| 3 | 0/0 0/0 ./. | 0/0 0/0 0/. | both fixed for allele 0, one missing and one half called genotype |
| 4 | 0/1 0/1 0/1 | 0/1 0/1 0/1 | the same in both, every genotype heterozygous |

The counts that everything comes from, allele 0, 1 and 2 and then the
called alleles and the called genotypes:

| variant | pop1 counts | n_1 | called gts 1 | pop2 counts | n_2 | called gts 2 |
|---|---|---|---|---|---|---|
| 1 | 5, 1, 0 | 6 | 3 | 1, 5, 0 | 6 | 3 |
| 2 | 1, 2, 3 | 6 | 3 | 5, 1, 0 | 6 | 3 |
| 3 | 4, 0, 0 | 4 | 2 | 5, 0, 0 | 5 | 2 |
| 4 | 3, 3, 0 | 6 | 3 | 3, 3, 0 | 6 | 3 |

Variant 3 shows the half called genotype: pop2 has 5 called alleles from
2 called genotypes and one half called one, so its frequency of allele 0
is 5/5 while it has 2 genotypes for the `min_num_samples` test.

From them, per variant,

| variant | H_b | H_w | f_2 | sum sqrt(p_1a p_2a) |
|---|---|---|---|---|
| 1 | 0.722222 | 0.333333 | 0.388889 | 0.745356 |
| 2 | 0.805556 | 0.533333 | 0.272222 | 0.608380 |
| 3 | 0 | 0 | 0 | 1 |
| 4 | 0.5 | 0.6 | -0.1 | 1 |

and the sums over the four variants are 2.027778 for H_b, 1.466667 for
H_w and 0.561111 for f_2. Each item below takes its number from here.
Variant 3, where both populations are fixed for the same allele, has
H_b = 0 and so a per variant F_ST of 0/0; it adds 0 to both sums and
nothing else, which is what the ratio of sums is for.

## Hudson's F_ST

### What it gives

How much of the diversity of two populations taken together lies between
them rather than within them, from 0 for two populations with the same
allele frequencies everywhere to 1 for two that share no allele anywhere.
It is the measure that relates to drift, migration and population size,
and it is what most SNP work reports.

It is the estimator Hudson, Slatkin and Maddison proposed and that
Bhatia, Patterson, Sankararaman and Price (2013, Genome Research 23:
1514, DOI 10.1101/gr.154831.113) recommend for SNP data, because it does
not move with the ratio of the sample sizes of the two populations as
Weir and Cockerham's does. With the H_b and H_w of the item above,

    F_ST = (sum over variants of H_b - sum over variants of H_w) / (sum over variants of H_b)

which is the same as the sum of f_2 over the sum of H_b. Summing the
numerators and the denominators over the variants before dividing, the
ratio of averages, is the second thing that paper asks for. The mean of
the per variant ratios instead gives weight to variants whose denominator
is near zero, which are the rare ones, and it is undefined for a variant
where both populations are fixed for the same allele. On the biallelic
panel, where 1 of the 1200 variants is such a variant, the mean of the
ratios for p0 and p2 is 0.087352 against 0.102736 for the ratio of sums,
15 in 100 lower, and for p0 and p1 it has no value at all.

For two alleles the formula is the one that paper writes, with p the
frequency of one allele and n the called alleles:

    numerator   = (p_A - p_B)^2 - p_A(1 - p_A)/(n_A - 1) - p_B(1 - p_B)/(n_B - 1)
    denominator = p_A(1 - p_B) + p_B(1 - p_A)

Both forms were computed on the biallelic panel and agree to the last
bit.

### How it is verified

Against plink2 v2.0.0-a.7.7, the arm64 build of 18 September 2026, which
computes this estimator and sums it the same way:

    plink2 --vcf panel.vcf.gz --pheno panel_pops.txt \
           --fst popcat method=hudson report-variants --out hud

`--pheno` loads the file of `IID` and `popcat` columns as a categorical
phenotype and `--fst popcat` takes the populations from it.
`report-variants` writes one file per pair with the F_ST of each variant
beside the summary.

On the biallelic panel plink2 gives 0.104962, 0.102736 and 0.109621 for
p0-p1, p0-p2 and p1-p2, and popnei's arithmetic gives 0.1049624450 for
the first: the largest distance from what plink2 prints is 4.9e-7 over
the three pairs, half a unit of the last of the six digits it prints, so
the tests compare within 1e-6 absolute. The per variant value of
`var0000` for p0 and p1, 0.332322, is the literal of the cargo test that
checks one variant.

On the multiallelic panel plink2 gives 0.0642281, 0.0694106 and
0.0699954, and popnei's 0.0642280666 for the first, the largest distance
being 4.0e-8. That run also settles that plink2 counts every allele of a
multiallelic record and does not reduce it: the same formula on the major
allele against the rest gives 0.0707950 and on the VCF's reference allele
against the rest 0.0598407, neither of which is what plink2 printed.

The worked example: 0.561111 / 2.027778 = 0.276712.

## f_2

### What it gives

How much allele frequency the two populations have drifted apart by,
0.04118 per variant between p0 and p1 of the biallelic panel. Where F_ST
divides that by the diversity the two hold, f_2 leaves it in the units it
was measured in, which is what makes it add up along a tree: the f_2
between two populations is the sum of the f_2 of the branches that
separate them, which is why it is what admixture graphs are built from.

It is the f_2 of Patterson et al. (2012, Genetics 192: 1065, DOI
10.1534/genetics.112.145037), estimated without the bias that the
sampling puts in, and it is the numerator of the F_ST above:

    f_2 = (sum over variants of H_b - sum over variants of H_w) / (the variants that counted)

The square root of f_2 is a Euclidean distance between the two
populations, the length of the difference of their frequency vectors with
the sampling noise taken out, so a principal coordinate analysis of the
matrix of square roots has no negative eigenvalues beyond rounding. The
square root is not given here: f_2 is negative for a pair of populations
that this dataset cannot tell apart, and a user who wants the matrix of
square roots takes it from `f2` themselves and decides what to do with a
negative one.

### How it is verified

Against ADMIXTOOLS 2.0.10, the R package that computes f_2 and the
f-statistics built on it, under R 4.6.1. It reads a plink binary fileset,
which plink2 writes from the panel's VCF, with the population of each
individual put in the first column of the `.fam`:

    blocks <- f2_from_geno(pref, maxmiss = 1, blgsize = 100000,
                           adjust_pseudohaploid = FALSE)
    f2(blocks)

`maxmiss = 1` keeps the variants that have missing genotypes, which its
default drops, and `adjust_pseudohaploid = FALSE` stops it from treating an
individual called homozygous everywhere as haploid.

On the biallelic panel it gives 0.04118111, 0.03989079 and 0.04279856 for
the three pairs, and the formula above gives the same eight digits for all
three, so the tests compare within 1e-8 absolute, the last digit it prints.

ADMIXTOOLS reads biallelic genotypes only, so the multiallelic panel is
checked through plink2 instead: f_2 summed over the loci is plink2's F_ST
times the sum of H_b, and for p0 and p1 that is 0.0642281 x 90.758172 =
5.829225 against the 5.829222 the sum gives, 3e-6 apart, which is what the
six digits of plink2 allow.

The worked example: 0.561111 / 4 = 0.140278.

## The chord distance and Nei's D_A

### What it gives

A distance between two populations that behaves as a distance should: the
Cavalli-Sforza and Edwards chord distance is the straight line between
two populations on the surface of a sphere, got by taking the square root
of every allele frequency, which puts each population on a sphere of
radius 1, and measuring the chord between them. It is Euclidean by
construction, so a principal coordinate analysis of a matrix of these
distances has no negative eigenvalues beyond rounding, which is not true
of F_ST or of Jost's D; and Takezaki and Nei (1996, Genetics 144: 389)
found it among the two that recover the right tree topology most often.
It is what popnei gives a user who asks for a tree or a principal
coordinate analysis of populations.

    D_A    = 1 - (sum over the variants that counted, and over the alleles,
                  of sqrt(p_Aa * p_Ba)) / (the variants that counted)
    chord  = sqrt(D_A)

D_A is Nei, Tajima and Tateno's distance (1983, Journal of Molecular
Evolution 19: 153), the other of the two Takezaki and Nei found best, and
it is the square of the chord distance, so both come from the one sum and
a user picks the one their method expects. Neither is corrected for the
sampling, and neither goes negative.

### How it is verified

Against `dist.genpop` of adegenet 2.1.11 under R 4.6.1 with `method=2`,
which is this chord distance. The genotypes reach R as the csv that
`df2genind` reads, which the Kosman item's reference script already
writes, with a population column added.

On the biallelic panel adegenet gives 0.18026704497001397,
0.17586558860911838 and 0.17977447554045811, and popnei's arithmetic is
1.25e-15 from them at the furthest; on the multiallelic panel
0.33853588707322202, 0.33760692274322368 and 0.3495821544731133, at
3.3e-16. Both are the last bits of a double, so the tests compare within
1e-12 relative. D_A is not in adegenet and is checked as the square of
what is.

Two of adegenet's other four distances were run and are not taken:
`method=4`, Rogers', and `method=5`, Prevosti's, give the same number on
the biallelic panel, 0.1672132615 for p0 and p1, because for two alleles
they are the same formula, and they part on the multiallelic one,
0.2287511165 against 0.3143145698. Neither adds anything to the chord
distance for a tree, and "Not in this spec" says so.

The worked example: the sum of the square roots is 3.353736, D_A =
1 - 3.353736/4 = 0.161566, and the chord distance is 0.401953.

## Jost's D

### What it gives

How different the alleles of two populations are, from 0 when they have
the same alleles at the same frequencies to 1 when they share none. Where
F_ST asks how near the populations are to being fixed, D asks how much of
the allelic variety is not shared, and the two answer different questions:
Jost, Archer, Flanagan, Gaggiotti, Hoban and Latch (2018, Evolutionary
Applications 11: 1139, DOI 10.1111/eva.12590), which is Jost's own
agreement with his critics, recommends reporting both.

D earns its place with many alleles. Nei's G_ST cannot reach 1 when the
populations are internally diverse: for two populations its largest
possible value is (1 - H_S)/(1 + H_S), so for microsatellites with an H_S
of 0.8 it cannot pass 0.11, and two populations sharing no allele at all
still give 0.11. D was built to be free of that. With two alleles H_S
cannot pass 0.5, the ceiling is 0.33 at the lowest, and Alcala and
Rosenberg (2019, Molecular Ecology 28: 1624, DOI 10.1111/mec.15000)
showed that the three measures are then constrained alike and that the
choice among them hardly matters. So D is the one to read on
microsatellites and F_ST on SNPs, and popnei gives both from the one
pass.

It is the D_est of Jost (2008), in the form pyNei computes, which is the
one GenAlEx prints. Per variant, with s = 2 populations:

    H_S = (E_A + E_B) / 2,      E_P = 1 - sum over a of p_Pa^k
    H_T = 1 - sum over a of ((p_Aa + p_Ba) / 2)^k

with k the ploidy, so that E_P is the expected heterozygosity of
`docs/specs/stats.md` and H_T the same over the two populations pooled at
equal weight. Both are corrected for the sample, as Nei and Chesser
(1983) do, with n the harmonic mean of the called genotypes of the two
populations at that variant and H_obs the mean of their observed
heterozygosities:

    H_S' = (n / (n - 1)) * (H_S - H_obs / (2n))
    H_T' = H_T + H_S' / (n * s) - H_obs / (2 * n * s)

The corrected values are summed over the variants that counted, each sum
divided by how many there were, and D comes from the two means:

    D = (s / (s - 1)) * (mean H_T' - mean H_S') / (1 - mean H_S')

These are `_calc_pairwise_dest` and `_calc_jost_from_ht_hs` of
`pynei/dists.py`, and popnei reproduces them.

### What pyNei does that is odd, and what popnei does instead

pyNei applies two thresholds and not one, and popnei applies one without
changing any value. `_count_alleles_per_var` of `pynei/gt_counts.py` sets
the frequencies of a population to NaN, which drops the variant, when its
called alleles divided by the ploidy are below `min_num_samples`, and
`_calc_pairwise_dest` then drops the variant again when the called
genotypes of either population are below it. The second test is the
stricter of the two and implies the first: a population's called alleles
are the ploidy times its called genotypes plus whatever its half called
genotypes contributed, so its called alleles divided by the ploidy are
never below its called genotypes. popnei tests the called genotypes alone,
which is the count the correction above divides by and the one rule that
all seven numbers share, and it keeps every variant pyNei keeps.

The correction of H_T divides by n times the number of populations, and
the number of populations is 2 because the calculation is pairwise, not
however many the user passed. That is what pyNei does, `num_pops = 2`
inside `_calc_pairwise_dest`, and popnei reproduces it: a pair is a pair
whatever else is in the run.

`hmean` in `pynei/dists.py` is pyNei's own harmonic mean and gives NaN
where the sum of the inverses is not finite, which is a population with
no called genotype. popnei does not reach that case, because such a
population has fewer than `min_num_samples` called genotypes for any
threshold above 0 and the variant is dropped before the mean is taken.

### How it is verified

Against `pairwise_D` of mmod 1.3.3 under R 4.6.1, on the same genotypes
read by `df2genind`. mmod computes a different estimator of the same
quantity (**Open 2**, below): its `HsHt` leaves the observed heterozygosity term out of both
corrections and uses 2n/(2n - 1) where the formula above uses n/(n - 1)
and subtracts H_obs/(2n). So the check is an agreement and not an
equality, and the spec states both formulas so that the next reader knows
which popnei computes.

On the biallelic panel mmod gives 0.0634704859, 0.0612312792 and
0.0656071139 where pyNei gives 0.0635434630, 0.0612981314 and
0.0656705213, 7.3e-5 apart at the furthest, with 48 to 84 individuals a
population. On the multiallelic panel, with 30 a population, mmod gives
0.1662094246, 0.1819580763 and 0.1816462134 against pyNei's
0.1661307946, 0.1822136338 and 0.1819999162, 3.5e-4 apart. The tests
compare popnei with mmod within 5e-4 on these two panels, which says that
popnei computes Jost's D and not another statistic, and with pyNei within
1e-12 relative, which is what pins the estimator.

The worked example: the mean corrected H_S is 0.34375 and the mean
corrected H_T is 0.425347, so D = 2 * 0.081597 / 0.65625 = 0.248677.

## Nei's G_ST and Hedrick's G'_ST

### What it gives

G_ST is Nei's fixation measure, the share of the diversity of the two
populations that lies between them, computed from the same corrected
means as Jost's D:

    G_ST = (mean H_T' - mean H_S') / (mean H_T')

G'_ST is Hedrick's (2005), G_ST divided by the largest value it could
take given how diverse the populations are, which for two populations is

    G_ST_max = (1 - mean H_S') / (1 + mean H_S')

so that G'_ST reaches 1 when the two populations share no allele. It is
what a reader of the microsatellite literature expects beside D, and both
cost nothing here: they are the two sums Jost's D already accumulates.

G_ST is not Hudson's F_ST, although both are called a fixation measure.
They differ in what they correct for and in how the variants are combined:
G_ST above is a ratio of means of per variant corrected values, and F_ST a
ratio of sums of uncorrected ones. On the biallelic panel mmod's G_ST for
p0 and p1 is 0.0554 and plink2's F_ST is 0.1050.

### How it is verified

Against `pairwise_Gst_Nei` and `pairwise_Gst_Hedrick` of mmod 1.3.3, with
the same difference of estimator as Jost's D above and the same tolerance
of 5e-4. On the biallelic panel mmod gives 0.0553896690, 0.0541614926 and
0.0579922831 for G_ST and 0.1617736271, 0.1576967932 and 0.1680418436 for
G'_ST; on the multiallelic panel 0.0331789783, 0.0359529575 and
0.0362672006, and 0.2197612680, 0.2387386980 and 0.2389275805. pyNei has
neither, so there is no comparison with it.

The worked example: G_ST = 0.081597/0.425347 = 0.191837, G_ST_max =
0.65625/1.34375 = 0.488372 and G'_ST = 0.392809.

## The Rust interface

The calculation over a reader. It asks the reader for the genotypes alone
and reads it to its end. Its errors: a pass that gave no variant, with
the counts of the filters of the reader; a sum
above `u32::MAX`; the two `u32` of every pair, asked of the machine at
the first block, that the machine has not the memory for, 400 MB at
10000 individuals; and those of the reader. The first is
`PassGaveNoVariant`, the case of `docs/specs/stats.md` that every
calculation over a pass raises; the second and the third are new
cases of the error of the crate. Each of the three is a `ValueError` in
Python, an
input that popnei cannot calculate on. The third was added by the plan
on 22 September 2026, as the block that the machine has not the memory
for is refused in `docs/specs/block.md`, so that a panel too large for
the machine is an error and not a process that ends.

```rust
pub fn calc_kosman_sums<R: BlockReader + ?Sized>(reader: &mut R) -> Result<KosmanSums>;
```

The reader is the outermost reader of the chain that `chain_of` of
`docs/specs/filters.md` built from the steps of the `Variants`, and the
binding crate keeps it, as it keeps the one it passes to `write_vars` of
`docs/specs/io_vars.md`: when the calculation returns it reads the counts
of the filters from the chain, and the number of variants the calculation
took from `KosmanSums::num_vars`, since no block of the pass reaches the
binding crate; together they are the `pass_stats` of the result.

What it gives: for every pair of individuals, k times the sum of d and n,
with k the ploidy of the reader.
The pairs are in the order of `dist_vector`, (0, 1), (0, 2), ..., (1, 2),
....

```rust
pub struct KosmanSums { /* private */ }

impl KosmanSums {
    pub fn num_individuals(&self) -> usize;
    /// The variants the calculation was given, called in a pair or not.
    pub fn num_vars(&self) -> u64;
    pub fn ploidy(&self) -> usize;
    /// The ploidy times the sum of d of the pair, and n, the variants at
    /// which both genotypes were called. The same for (i, j) and (j, i).
    /// None when i == j or when either is not an individual.
    pub fn sums(&self, i: usize, j: usize) -> Option<(u32, u32)>;
    /// The distance of the pair, the first of `sums` over the ploidy
    /// times the second. None when n is 0
    /// or below `min_num_vars`, and where `sums` gives None.
    pub fn dist(&self, i: usize, j: usize, min_num_vars: u32) -> Option<f64>;
    /// The distance of every pair, in the order of `dist_vector`. The
    /// binding crates write NaN for a None.
    pub fn dists(&self, min_num_vars: u32) -> impl Iterator<Item = Option<f64>> + '_;
}
```

`min_num_vars` is the `min_num_snps` of Python under the word of the
glossary, and 0 for `None`.

The distances between populations. The pass reads its reader to the end
and asks it for the genotypes, and for the chromosome and the position
too when the groups are stretches of a chromosome. Its errors: a pass that
gave no variant, `PassGaveNoVariant`, the case of `docs/specs/stats.md`
that every calculation over a pass raises; fewer than two populations;
fewer resampling groups than the standard errors need, Open 1; and the
memory for the sums, asked of the machine at the first block. Each is a
`ValueError` in Python. `Pops` is the populations of
`docs/specs/stats.md`, which already refuses a name that is not an
individual of the reader and a population with no individual.

```rust
pub struct PopDistOptions {
    /// How many called genotypes a pop needs at a variant for that variant
    /// to count for a pair the pop is in.
    pub min_num_samples: u32,
    pub groups: JackknifeGroups,
}

/// How the variants are cut into the groups that the standard errors are
/// resampled over.
pub enum JackknifeGroups {
    /// No standard errors, and the pass does not ask for the positions.
    None,
    /// Each variant its own group.
    PerVariant,
    /// Stretches of one chromosome this many base pairs long, each one
    /// anchored on its own first variant.
    OfBasePairs(u64),
}

pub fn calc_pop_dist_sums<R: BlockReader + ?Sized>(
    reader: &mut R,
    pops: &Pops,
    options: &PopDistOptions,
) -> Result<PopDistSums>;
```

What the pass gives, and every measure out of it. The pairs are in the
order of `dist_vector`, (0, 1), (0, 2), ..., (1, 2), ..., over the
populations in the order `pops` has them.

```rust
pub enum PopDistMeasure { Fst, F2, Chord, Da, Dest, Gst, GstHedrick }

/// One group of variants: the chromosome, and the first and the last
/// position it holds, both included.
pub struct GroupId { pub chrom: u32, pub start: u64, pub end: u64 }

pub struct PopDistSums { /* private */ }

impl PopDistSums {
    pub fn num_pops(&self) -> usize;
    /// The variants the pass was given, counted for a pair or not.
    pub fn num_vars(&self) -> u64;
    pub fn groups(&self) -> &[GroupId];
    /// The variants that counted for the pair: both its pops had at least
    /// `min_num_samples` called genotypes. None when i == j or when either
    /// is not a pop.
    pub fn num_vars_of(&self, i: usize, j: usize) -> Option<u64>;
    /// The measure for the pair. None where `num_vars_of` is 0 or None,
    /// and, for the three that need the corrected heterozygosities, where
    /// the harmonic mean of the called genotypes was 1 at every variant.
    pub fn measure(&self, measure: PopDistMeasure, i: usize, j: usize) -> Option<f64>;
    /// Its jackknife standard error. None where `measure` is None, where
    /// no groups were asked for, and where every variant of the pair fell
    /// in one group.
    pub fn standard_error(&self, measure: PopDistMeasure, i: usize, j: usize) -> Option<f64>;
    /// f_2 within one group, which f_3 and f_4 are built from later.
    pub fn f2_of_group(&self, group: usize, i: usize, j: usize) -> Option<f64>;
    /// The measure of every pair, in the order of `dist_vector`. The
    /// binding crates write NaN for a None.
    pub fn measures(&self, measure: PopDistMeasure) -> impl Iterator<Item = Option<f64>> + '_;
}
```

## Speed

The dataset is 100000 variants x 1000 individuals, biallelic, 3 in 100
genotypes missing, read from a vars file, with the time of reading that
file alone measured beside it and taken out, so that the reader is not in
the number; on a vars file of 20000 variants of 1000 individuals the
reader took 18.8 ms, section 1 of the architecture. pyNei takes 0.76 s with `num_threads=1` and 0.30 s with 6, which is its
best, on the owner's M5 Pro on 21 September 2026 with numpy 2.5.3 on
Accelerate. The trial of "How it runs", with no tuning, takes 0.044 s a
block, 0.88 s for the 20 blocks, on one thread, 0.017 s a block, 0.34 s,
on 18 cores, and 0.065 s a block, 1.3 s, in wasm. The numbers to reach
are the trial's with a tenth over them: 0.97 s on one thread, 0.38 s on
18 cores and 1.43 s in wasm. That is 1.3 times pyNei's time on this
machine, which the decision for the bits accepted.
pyNei under pyodide was not measured. No other ploidy was timed: a
biallelic tetraploid block has the 9 sets of the 4 allele diploid block
of the table, 0.067 s, and a haploid one has 3.

At the largest dataset of the objectives, a million variants x 10000
individuals, the trial's 0.090 s a block on 18 cores is 3 minutes for the
2000 blocks, and its 0.349 s on one thread 12 minutes.

What the code reaches. The three numbers were missed when the plan
`docs/plans/dists-kosman.md` finished, 1.154 s, 0.625 s and 2.130 s, and
the performance review of 22 September 2026,
`docs/reports/perf-dists-kosman-2026-09-22.md`, met them all: 0.767 s on
one thread, 0.102 s on 18 cores and 1.106 s in wasm, on the dataset
above with the reading taken out, on the owner's M5 Pro. So popnei is 6
in 100 slower than pyNei on one thread and 2.7 times faster than
pyNei's best. Of the two things this section left to that review, the
building of the sets of a block on the threads was measured and is in
the code, a run of 64 individuals to a work item; fewer sets for a block
with two alleles was not tried, and the review's finding H4 says what it
would take and what it would give. The review also gave the browsers
popnei runs in a floor, goal 3 of `docs/objectives.md`, because the
count of the bits of a pair in wasm uses the vector instructions of
WebAssembly, which is what took that build under its number.

## Open points

The owner decides these two. Until then the implementer follows the
"meanwhile" of each.

**Open 1: the default of `jackknife_group`, and a dataset with few
groups.** A standard error from a handful of groups is not one to quote,
and the number of groups a dataset gives depends on how its variants are
spread, which popnei learns only while reading. The options are to default
to 5 000 000 base pairs, the length the f-statistics literature uses, and
to raise when fewer than 20 groups come out, naming the number and the
argument that fixes it, which makes the first call of a user with 300
microsatellite loci on one chromosome fail until they pass `"variant"`; to
default to 5 000 000 and give the standard errors whatever the number of
groups, with the number in the result, which lets a user quote a standard
error built from 3 groups without being told; or to default to `None`, no
standard errors unless they are asked for, which costs the user who wants
them one argument and never surprises the one who does not.
Recommendation: the first. It is the one that says what is wrong instead
of giving a number that looks like the others, and one argument fixes it.
Meanwhile the implementer writes the first, with 20 as the threshold.

**Open 2: which estimator of Jost's D.** popnei reproduces pyNei's, the
Nei and Chesser correction as GenAlEx prints it, which subtracts the
observed heterozygosity term. mmod, the R package a user is likeliest to
check against, leaves that term out and uses 2n/(2n - 1) where pyNei uses
n/(n - 1). The two are 7.3e-5 apart on the biallelic panel, with 48 to 84
individuals a population, and 3.5e-4 apart on the multiallelic one, with
30. The options are to reproduce pyNei, which keeps every number pyNei's
users have and leaves popnei a few units of the fourth decimal from what R
prints; to compute mmod's instead, which changes every D pyNei ever gave
and has no name in the literature of its own; or to give both, as `dest`
and a second field, which costs one more pair of sums for each pair and
each group and one more name in the result. Recommendation: reproduce
pyNei, and say in the doc comment of `dest` which estimator it is and by
how much mmod's differs on these two panels. Meanwhile the implementer
reproduces pyNei.

The Kosman item above has none of its own. The two the first version of this spec had, whether the function
takes `num_threads` and whether ploidies other than 2 are taken, the
owner decided on 22 September 2026, and both are written under "Its
Python function" with the option not taken. The threads: no calculation
of popnei has the argument. The ploidies: every one is taken, since the
paper's formula 2 is for any ploidy and `gd.kosman` computes it, which
the tetraploid and the haploid datasets of "How it is verified" show.

## Not in this spec

- f_3 and f_4, the statistics of three and four populations that admixture
  graphs are built from, and the tests on them: a spec of their own. They
  are sums and differences of the f_2 of pairs, so `f2_groups` of the
  result here is what they are computed from, with no reading of genotypes
  and no new pass.
- R_ST and the squared difference in repeat number, the microsatellite
  measures that use how many repeats each allele has: popnei's genotypes
  carry allele numbers and not allele sizes, so nothing here can compute
  them. Carrying the sizes is a change to section 4 of
  `docs/architecture.md`.
- Weir and Cockerham's theta, and the population-specific F_ST of Goudet
  and Weir (2023), which `hierfstat` computes: neither is in pyNei, and
  what popnei gives is the pairwise measure that plink2's Hudson method
  checks. `hierfstat` does not install on the owner's machine, its
  dependency gaston failing to build, so it is not a reference program
  here either.
- Nei's standard distance D_s, Prevosti's and Rogers', which
  `adegenet::dist.genpop` gives beside the chord distance: D_s and
  Prevosti's are not Euclidean, Rogers' is the same as Prevosti's for two
  alleles, and none of the three adds anything to the chord distance for a
  tree or for a principal coordinate analysis. pyNei has none of them.
- `calc_pairwise_euclidean_dists` of pyNei, the Euclidean distances
  between the rows of a frame: nothing in the architecture asks for it,
  and pyNei used it for the approximate algorithm.
- The principal coordinate analysis of any of these distances,
  `do_pcoa_from_variants` and `do_pcoa`: `docs/specs/pca.md`. It will pass
  `min_num_snps` on and will have no `use_approx_embedding_algorithm`. The
  chord distance and the square root of f_2 are Euclidean and the other
  measures here are not, so that spec has to say what it does with the
  negative eigenvalues a matrix of F_ST or of Jost's D gives it.
- The read ahead thread, which reads the next block while the calculation
  works on the one in hand, and which `docs/specs/block.md` leaves for the
  first calculation that consumes blocks: it is a reader over a reader, and
  `calc_kosman_sums` takes any reader, so it is an item of that spec and
  nothing here changes with it.
- `Variants.from_gt_array`, the `Variants` from an array of genotypes:
  `docs/specs/variant.md`, for a later plan.
