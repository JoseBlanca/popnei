# The dists module: distances between individuals and between populations

21 September 2026. The `dists` module tells a user of popnei how far apart
their individuals are, one distance for every pair of them, and how far
apart their populations are. There is no code. This spec develops the row
`dists` of the table in section 9 of `docs/architecture.md`, and it has one
item so far, the Kosman distance between individuals; the item for Jost's
D between populations is not written. It depends on `docs/specs/block.md`,
which has the block, the run of consecutive variants held as arrays, and
`BlockReader`, the trait of everything that gives blocks; on
`docs/specs/variant.md`, which has the `Variants` a user holds, with the
steps put on it and the counts of a pass, and the `Variants` built from an
array of genotypes, an item written with this spec because its tests need
it; and on `docs/specs/filters.md`, which builds the chain of filters of
a pass.

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
  finished" of `docs/specs/filters.md` says; the binding crate reads
  them from the chain before it raises. pyNei raises a `RuntimeError`,
  which in popnei is the exception of a defect of popnei, and this is an
  input that cannot be calculated on.
- The result does not carry n(i, j), as pyNei's does not. The core result
  has it, so that giving it to the user later changes no Rust.
- The result has `pass_stats`, which pyNei's has not; pyNei keeps the
  counts of its filters in its `Variants`.

The first three were the owner's decisions, the last follows from his
decision of 21 September 2026 on the consumers, in `docs/specs/filters.md`,
and the other six are decided here: none changes a distance that pyNei gives rightly.

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
files and on the genotypes of the worked example below, the last through
`Variants.from_gt_array`, with `min_num_snps` of `None` and of a value
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

## The Rust interface

The calculation over a reader. It asks the reader for the genotypes alone
and reads it to its end. Its errors: no variant in the reader; a sum
above `u32::MAX`; and those of the reader. The first two are new cases
of the error of the crate, and each is a `ValueError` in Python, an input
that popnei cannot calculate on.

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

## Speed

The dataset is 100000 variants x 1000 individuals, biallelic, 3 in 100
genotypes missing, from memory, so that the reader is not in the time.
pyNei takes 0.76 s with `num_threads=1` and 0.30 s with 6, which is its
best, on the owner's M5 Pro on 21 September 2026 with numpy 2.5.3 on
Accelerate. The trial of "How it runs", with no tuning, takes 0.044 s a
block, 0.88 s for the 20 blocks, on one thread, 0.017 s a block, 0.34 s,
on 18 cores, and 0.065 s a block, 1.3 s, in wasm. The numbers to reach
are the trial's with a tenth over them: 0.97 s on one thread, 0.38 s on
18 cores and 1.43 s in wasm. That is 1.3 times pyNei's time on this
machine, which the decision for the bits accepted. Whether fewer sets for
a block with two alleles and a parallel building of the sets bring popnei
under pyNei's 0.76 s is for a performance review, and was not measured.
pyNei under pyodide was not measured. No other ploidy was timed: a
biallelic tetraploid block has the 9 sets of the 4 allele diploid block
of the table, 0.067 s, and a haploid one has 3.

At the largest dataset of the objectives, a million variants x 10000
individuals, the trial's 0.090 s a block on 18 cores is 3 minutes for the
2000 blocks, and its 0.349 s on one thread 12 minutes.

## Open points

None. The two the first version of this spec had, whether the function
takes `num_threads` and whether ploidies other than 2 are taken, the
owner decided on 22 September 2026, and both are written under "Its
Python function" with the option not taken. The threads: no calculation
of popnei has the argument. The ploidies: every one is taken, since the
paper's formula 2 is for any ploidy and `gd.kosman` computes it, which
the tetraploid and the haploid datasets of "How it is verified" show.

## Not in this spec

- Jost's D between populations, pyNei's `calc_jost_dest_pop_dists`: the
  next item of this spec. It gives a `Distances` too.
- `calc_pairwise_euclidean_dists` of pyNei, the Euclidean distances
  between the rows of a frame: nothing in the architecture asks for it,
  and pyNei used it for the approximate algorithm.
- The principal coordinate analysis of these distances,
  `do_pcoa_from_variants`: `docs/specs/pca.md`. It will pass
  `min_num_snps` on and will have no `use_approx_embedding_algorithm`.
- The read ahead thread, which reads the next block while the calculation
  works on the one in hand, and which `docs/specs/block.md` leaves for the
  first calculation that consumes blocks: it is a reader over a reader, and
  `calc_kosman_sums` takes any reader, so it is an item of that spec and
  nothing here changes with it.
- `Variants.from_gt_array`: `docs/specs/variant.md`.
