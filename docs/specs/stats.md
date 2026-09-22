# The stats module: the per variant and the per individual statistics, per population

September 2026. The `stats` module tells a user of popnei how their
variants and their individuals look before any distance or association is
calculated: for each population, how heterozygous the variants are, how
frequent their commonest allele is, how much diversity each population
holds and how many variants vary at all, each as a mean and a histogram
over the variants; and for each individual, how much of its data is
missing and how heterozygous it is. It is also where a population, a
named set of individuals, is turned into what the core works on, which
Jost's D, the linkage disequilibrium per population and the kinship will
take from here. There is no code. This spec develops the row `stats` of
the table in section 9 of `docs/architecture.md`, whole, and the counts of
one variant over one population, which `docs/specs/variant.md` leaves to
"the first statistic per population".

It depends on `docs/specs/block.md`, which has the block, the run of
consecutive variants held as arrays, and `BlockReader`, the trait of
everything that gives blocks; on `docs/specs/variant.md`, which has the
`Variants` a user holds, with the steps put on it and the counts of a
pass, the counts of the genotypes and of the alleles of one variant, and
the `Variants` built from an array of genotypes, which the tests use; and
on `docs/specs/filters.md`, which builds the chain of filters of a pass
and, in its item for the filter of individuals, written with this spec,
takes individuals out of it. Two specs it points to are not on `main`
yet: `docs/specs/dists.md`, the distances, on the branch
`plan/dists-kosman`, and `docs/specs/pca.md`, the principal components,
on `plan/pca`.

The expected heterozygosity was written first, in September 2026, as the
example the owner chose of how a spec is written, read by the first
reader and checked by the reviewer against pyNei. It is here as it was,
under the words of `docs/glossary.md`, individual for pyNei's sample, and
with its open points numbered among the others.

## The populations

### What they give

A population is a named set of individuals that a calculation treats as a
group, and `pops` is how a user names them: a dict of population name to
the names of its individuals. Every statistic of this module is calculated
for each population over its individuals alone, and the results are keyed
by population name. When no `pops` is given there is one population,
named `pop`, of every individual, inherited from pyNei's `DEF_POP_NAME`.

The names are resolved once for each pass, against the individuals that
the pass gives, which are those of the source after the filter of
individuals of `docs/specs/filters.md` when the `Variants` has one. So a
user who keeps some individuals and then names populations names them
among the kept ones.

### In Python and in TypeScript

`pops: dict[str, Sequence[str]] | None = None` is an argument of every
function of this module that works per population, with pyNei's name and
type, `Pops` of `pynei/utils_pop.py`. The rules are those of
`_calc_pops_idxs` of that file, with three refusals more:

- A name that is not among the individuals of the pass is a `ValueError`
  that names the population and the individual, as in pyNei.
- A value that is not a sequence of names, one name as a string among
  them, is a `ValueError`, as in pyNei.
- A name that is twice in one population is a `ValueError` that names
  both. pyNei counts that individual twice: `{"x": ["a", "b", "b"]}` over
  the genotypes `0/0 0/1 1/1` gives an observed heterozygosity of 2/3, run
  at commit ef0ca6e. A count that is wrong and says nothing is what
  popnei never gives, by the owner's rule of 21 September 2026 that no
  error passes in silence and a wrong input of a function is a
  `ValueError`.
- A population with no individual, and a `pops` with no population, are a
  `ValueError`. pyNei gives a result with NaN for the first and one with
  no columns for the second. An individual in two populations is taken,
  as in pyNei: `test_maf_stats` of `test/test_gt_counts.py` names one
  individual in both of its populations.

The order of the populations in every result is the order in which the
keys of `pops` iterate, which in Python is the order they were inserted
in, for every statistic; with no `pops` there is the one population,
`pop`. pyNei sorts the names for the expected heterozygosity and keeps
the order of the keys for the others, so its columns do not line up;
nothing in a value changes.

In TypeScript `pops` is an object of population name to an array of
names, `{pop1: ["a", "b"]}`, and the errors are an `Error` at the call.

### How it runs

The binding crate hands the names to the core, which looks each one up
among `individuals()` of the outermost reader of the chain, once the chain
of the pass is built, and holds for each population the indices of its
individuals in the order the user gave them. The lookup builds one hash
map of name to index for each population, since "The Rust interface" has
`from_names` resolve each population with `resolve_individuals` of
`docs/specs/filters.md`, which builds its own: 30 ms for 300 populations
of 10000 individuals, once per pass.

### How it is verified

Cargo tests at `Pops::from_names` of "The Rust interface": the indices of
the two populations of the worked example below; a name that is not an
individual, a name twice in one population, an empty population and no
population are each the error that names what is wrong; an individual in
two populations is taken. The pytest tests are those of `test_pops.py` of
pyNei, `test_every_stat_takes_the_same_pops` among them, which asserts
that each of its four statistics keys its result by the two population
names and refuses an unknown individual, which popnei asserts for its
five, and one test that the refusal of
a duplicated name is a `ValueError`.

## The counts of one variant over a population

### What they give

The two counts of `docs/specs/variant.md`, of the genotypes and of the
alleles, over the individuals of one population and not over all of them:
how many genotypes of the population are called, missing and
heterozygous at this variant, and how often each allele was called among
them, with the sum of the latter, the called alleles of the population.
What a called, a missing and a heterozygous genotype are, and that a half
called genotype is missing while its called allele is counted, is written
in that spec and holds here.

### How it runs

One pass over the genotypes of the individuals of the population, taken
from the row by their indices, into the same array of counts as the two
functions of `docs/specs/variant.md`, with their rules. When the
population is every individual in the order of the source, the caller
uses those two, which read the row as it is.

### How it is verified

Cargo tests at the two functions of "The Rust interface", on the six
variants of five diploid individuals of the worked example of
`docs/specs/filters.md`, with the two populations of the worked example
below: the counts of variant 1 over pop1, i1 and i2, are 2 called, 0
missing, 1 heterozygous, and alleles 0 and 1 counted 3 and 1, 4 called
alleles; over pop2, i3, i4 and i5, 2 called, 1 missing, 0 heterozygous,
allele 0 counted 5, 5 called alleles. An index that is not an individual
of the row is an error.

## The per variant distributions

### What they give

One pass over the variants that calculates up to five statistics for
every variant and every population, and gives back, for each statistic
and population, the mean over the variants that have a value and a
histogram of them. The per variant values are not kept, because a million
of them for each population do not fit a browser tab, and a user who
wants them takes the genotypes with `iter_blocks`. The five statistics
are the observed heterozygosity, the major allele frequency, the expected
heterozygosity, plain and unbiased, and the polymorphism ratio, each its
own item below; the last is a count and not a distribution.

A variant has no value of a statistic in a population when the population
has too little data at it: none at all, or, for all but the observed
heterozygosity, fewer called genotypes than `min_num_individuals`. A
variant with no value is out of
the mean and in no bin of the histogram, so the mean of a population is
over the variants that had enough data, and the histograms of two
populations can count different numbers of variants.

### In Python and in TypeScript

```python
def calc_per_var_distribs(
    variants: Variants,
    stats: Iterable[PerVarStat] = tuple(PerVarStat),
    pops: dict[str, Sequence[str]] | None = None,
    min_num_individuals: int = 20,
    hist_kwargs: dict | None = None,
    ploidy: int | None = None,
    poly_threshold: float = 0.95,
) -> PerVarDistribs
```

It is a consumer of the `Variants`, as `docs/specs/variant.md` has them:
it makes one pass over the source through the steps the `Variants` has
when it is called, and the `Variants` is as it was afterwards.

`stats` says which of the five to calculate: members of the `StrEnum`
`PerVarStat`, `OBS_HET`, `MAF`, `EXP_HET`, `UNBIASED_EXP_HET` and
`POLY_VARS_RATIO`, whose values are the names of the fields of the
result. All five by default. Anything that is not a member, a string
among them, is a `TypeError` that says so, and no statistic at all is a
`ValueError`. Asking for fewer is a saving of work and changes no value.
The owner decided both on 22 September 2026: the unbiased expected
heterozygosity is a statistic of its own and not a switch on the
expected heterozygosity, so that a user asks for the two like any two,
and the members alone are what a user can write without a typo that
passes. The options not taken were pyNei's `unbiased_exp_het` argument
and its `stats` of members, names and one name as a string.

`min_num_individuals` is how many called genotypes a population needs at a
variant for the variant to have a value there, 20 by default, inherited
from pyNei's `MIN_NUM_SAMPLES_FOR_POP_STAT`. Nobody has measured whether
20 is the right threshold. The test is on the called data counted in
genotypes, the called alleles of the population divided by the ploidy,
which can be a half when a genotype is half called, and the variant has
no value when that number is strictly less than the threshold: five
diploid individuals with four genotypes called and one half called have
4.5, and keep their value at a threshold of 4 and lose it at 5. Every
statistic is held to it, the observed heterozygosity too, which pyNei
exempts, under its item.

`hist_kwargs` is the histogram: `range`, the two ends, `(0, 1)` by
default, which is where the statistics live; `num_bins`, 40 by
default; and `bin_type`, `"linear"` for bins of equal width or
`"logarithmic"` for bins of equal ratio, whose `range` has to start above
0. The dict is read and not changed, and a key that is none of the three
is a `ValueError` that names them: pyNei ignores such a key, and a user
who writes `nbins` gets the 40 bins of the default with nothing said,
which is a result that is not the one they asked for and that says so
nowhere. The edges are `num_bins + 1`
numbers, computed as numpy's `linspace` computes them, the start plus i
times the width for the i-th and the end for the last, so that a value
that falls on an edge falls on the same side in popnei and in pyNei;
the logarithmic ones are 10 to the power of the same over the logarithms
of the ends. A value falls in the bin whose left edge is at most the value
and whose right edge is above it, and the last bin takes its right edge
too, as `numpy.histogram` does, so an observed heterozygosity of exactly 1
is in the last bin. A value outside the range is in no bin and in the
mean. pyNei spells the first bin type `"lineal"`, the Spanish word; the
owner decided on 22 September 2026 that popnei spells it `"linear"` and
refuses `"lineal"` as any other unknown name, as `docs/specs/pca.md`
spells `standardize_data` where pyNei has `standarize_data`.

`ploidy` is the two expected heterozygosities' alone, and
`poly_threshold` the polymorphism ratio's; their items say what they do.

`PerVarDistribs` is pyNei's frozen dataclass with two fields more:
`obs_het`, `maf`, `exp_het` and `unbiased_exp_het`, each a `StatsDistrib`
or `None` when it was not asked for; `poly_vars_ratio`, a `PolyVarsStats`
or `None`; and `pass_stats`, the
`PassStats` of `docs/specs/variant.md` that every result of a consumer
has, with how many variants the pass gave, after the steps, and how many
each filter was given and kept. `StatsDistrib` is pyNei's: `mean`, a
pandas series with one value per population, NaN for a population in
which no variant had a value; `hist_bin_edges`, a numpy array of
`num_bins + 1` float64; and `hist_counts`, a pandas frame of bins x
populations of integers. `PolyVarsStats` is pyNei's five series, under its
item.

It mirrors `calc_per_var_distribs` of `pynei/per_var_stats.py`. The
differences:

- `exp_het` of the result is the plain expected heterozygosity, and the
  unbiased one is `unbiased_exp_het`, where in pyNei `exp_het` holds
  whichever `unbiased_exp_het=True` chose, the unbiased one by default.
  A user who reads `exp_het` of both libraries reads two numbers.
- `stats` takes members of `PerVarStat` alone, as said above.
- `min_num_individuals` is pyNei's `min_num_samples` under the word of
  `docs/glossary.md`, as `individuals` of `Variants.from_gt_array` is
  pyNei's `samples`.
- There is no `num_threads`. The owner decided on 22 September 2026, in
  `docs/specs/dists.md`, that no calculation of popnei has it: the threads
  are those of rayon's pool.
- `min_num_individuals` is a whole number, 0 or more; what is no number
  is a `TypeError` and a negative one a `ValueError`. pyNei takes any
  number, 3.1 in one of its tests.
- `poly_threshold` is a number from 0 to 1, both included; NaN and
  anything outside are a `ValueError` and what is no number a
  `TypeError`. pyNei compares with whatever it is given.
- A pass that gives no variant is a `ValueError`, whether the source has
  none or the steps kept none, with the message that `docs/specs/dists.md`
  gives that error. pyNei raises a `ValueError` too, from
  `run_chunk_calcs` of `pynei/pipeline.py`.
- A duplicated name in a population, an empty population and an empty
  `pops` are refused, under "The populations".
- The result has `pass_stats`, and the populations of every statistic are
  in the order of the keys of `pops`, under "The populations".
- The differences of each statistic are under its item, and the open
  points among them are in the list at the end.

In TypeScript it is `calcPerVarDistribs(variants, {stats, pops,
minNumIndividuals, histKwargs, ploidy, polyThreshold})`, with `stats` an
array of the five names as a union type of string literals,
`"obs_het" | "maf" | ...`, which is what an enum is in TypeScript,
`histKwargs` an object with `range`, `numBins` and `binType`, and the
same defaults. The result has `pops`, the population names in their
order; `obsHet`, `maf`, `expHet` and `unbiasedExpHet`, each `null` or a
`StatsDistrib` with `mean`, a `Float64Array`
with one value per population and NaN for none, `histBinEdges`, a
`Float64Array`, and `histCounts`, a `Uint32Array` of populations x bins,
population after population; `polyVarsRatio`, `null` or an object with
the five arrays of its item; and `passStats`.

### What pyNei asserts, and the size of the blocks

`test_per_var_stats.py` of pyNei asserts, on 100 random biallelic
variants of 30 individuals in two populations of 15, that all four
statistics are calculated by default and only the ones asked for
otherwise; that `stats` takes members, names and one name, where popnei's
test asserts that a name is a `TypeError`; that an unknown name and no
name are a `ValueError`; that the means and the histogram counts of one
pass over the four are those of four passes of one each, which popnei
asserts over its five; that the arguments of one statistic change no
other; and that a `Variants` whose filter kept no variant is an error.
`test_hist.py`
asserts that `hist_kwargs` is not changed, that `range` and `num_bins`
give the edges of the example above, that a logarithmic range from 0.01
to 100 in 4 bins has the edges 0.01, 0.1, 1, 10 and 100, and that a bin
type that is neither of the two and a logarithmic range that starts at 0
or below are a `ValueError`. `test_per_var_distribs_with_threads` of
`test_threads.py` asserts that 2 or 4 threads change neither a mean nor a
histogram count. popnei has each of them as a pytest test.

No test of pyNei varies the size of its chunks for these statistics, and
one case below, a chunk with nothing called, is where the expected
heterozygosity of pyNei depends on them. popnei adds a cargo test that
the size of the blocks changes no histogram count and no mean beyond
1e-12 relative, which holds because popnei gives such a block no value.

### How it runs

Over the rows of each block, with rayon across them, and it needs no
`reblock` before it, the reader of `docs/specs/block.md` that puts the
blocks of its source back to one size. For each row and each population
the counts of "The
counts of one variant over a population" are taken once, and the
statistics that were asked for follow from those counts: the observed
heterozygosity from the genotype counts, the other three from the allele
counts. Each thread adds what it finds into accumulators of its own, and
the accumulators of the threads are added when the block is done: per
statistic and population a sum of the values as `f64`, a count of the
variants that had a value, and the counts of the bins; for the
polymorphism ratio its three counts. What is kept from one block to the
next is those accumulators, whose size is populations x bins per
statistic and grows neither with the variants nor with the individuals.
The mean is the sum over the count, once, at the end. The rows of a block
are added up in chunks of a fixed number of rows, the chunks are
collected in the order of the block and added one after another, so two
runs of the same pass over the same blocks add the same values in the
same order and give every mean to the bit, whatever the number of
threads: the same pass over `many.vcf` at 1, 2, 3, 4, 8 and 16 threads,
measured on 22 September 2026, gives the same bits of every mean, and the
cargo test of the thread pools compares the bits. Where the blocks are
cut changes which rows are added together, and with them the last bits of
a mean, which "How it is verified" allows for at 1e-12 relative; the
counts are integers and agree exactly.

With `pops` given, the genotypes of a population are gathered from the
row by its indices for each variant. Whether a population that is every
individual in the order of the source skips the gather is the
implementer's choice, as long as "Speed" is met.

The calculation asks its reader for the genotypes alone.

### How it is verified

What is common to the five: against pyNei, and a worked example. Each
statistic has, under its item, the reference program that checks its
values and the literals of the first cargo tests.

Against pyNei, a pytest test at `calc_per_var_distribs` of both
libraries, on two datasets. The first is the panel of pyNei's
`test/gwas_reference/sim_missing.vars`, 200 individuals named `s000` to
`s199`, 1200 biallelic diploid variants, 3 in 100 genotypes missing whole.
It is a vars file of pyNei, which popnei does not read, so the reference
script of this module, `tests/reference/stats/make_reference.py`, reads
it with pyNei and writes it as `tests/reference/stats/panel.vcf.gz`, 87
KB, with the ids `var0000` to `var1199`; its three populations, `p0`,
`p1` and `p2`, of 48, 68 and 84 individuals, are the `pop` column of
`test/gwas_reference/phenotypes.csv` of pyNei, which the script writes
beside it as `panel_pops.txt`. The same script runs every plink2 and
bcftools command of this spec and stores their reports beside itself,
under the names each item gives, and checks the literals below against
them. The second dataset is
`tests/reference/vcf/many.vcf` of `docs/specs/io_vcf.md`, 500 variants of
50 diploid individuals, `ind00` to `ind49`, one in ten with three alleles,
257 half called genotypes among its 25000, read with every variant given,
with two populations, `popA` of the first 20 individuals and `popB` of
the other 30. Both libraries are run with every statistic, pyNei twice
for its two values of `unbiased_exp_het`, `min_num_individuals` of 20 on
the panel and of 5 on `many.vcf`, and the default histogram; and once
more with no `pops`. The histogram counts and the counts of the polymorphism
ratio have to be equal, and the means and the ratios equal within 1e-12
relative, with NaN in the same places, because numpy and the Rust loop
add the variants of a population in different orders.

The histogram counts of the unbiased expected heterozygosity are the one
exception. The two libraries reach that value by different arithmetic:
popnei multiplies the k factors of each of its terms one over another,
`(c_a / c) · ((c_a - 1) / (c - 1))` at k = 2, and pyNei multiplies the
plain value by `c / (c - 1)`. Over every split of a biallelic variant up
to 500 called alleles, the two differ by 2.84e-16 at most, at the counts
3 and 476 of 479, and by up to 252 units in the last place, at the counts
1 and 269 of 270, measured on 22 September 2026 by computing both ways in
float64; the 1e-9 the comparison below allows rests on the absolute
bound. A variant whose value lies on an edge of the histogram falls on
either side of it, and seven variant and population pairs of the two
datasets do, all measured on 22 September 2026 with the default histogram
of 40 bins from 0 to 1. `var235` of `many.vcf`, at the position 9695,
with no `pops`, has the allele counts 55 and 45 of its 100 called
alleles, whose unbiased value is exactly 0.5; popnei gives 0.5 and pyNei
0.4999999999999999, so popnei counts it in the bin that starts at 0.5 and
pyNei in the one below. `var0978` of the panel in `p2` has 106 and 54 of
160, whose value is exactly 0.45; popnei gives 0.45000000000000007 and
pyNei 0.44999999999999996, one on each side of the edge 0.45. The other
five are in `popA` of `many.vcf`, each with the allele counts 21 and 15
of 36 and the same 0.5 against 0.4999999999999999: `var242`, the variant
after it, at the position 9991 and with no id, `var422`, `var428` and
`var466`. So the counts of
that one statistic are compared allowing the variants whose value lies
within 1e-9 of an edge to fall on either side of it, which the test
counts from pyNei's own per variant values, and the counts of the other
four are compared exactly. The default of
`min_num_individuals` is asserted by a call without it on the panel, with
`pops` of one population of 15 individuals: every mean is NaN with an
empty histogram and the counts of the polymorphism ratio are 0. The edges of the histogram are compared exactly for
bins of equal width and within 1e-12 for logarithmic ones, whose powers
of 10 two libraries need not round alike.

The worked example, which becomes the first cargo tests: the six
variants of five diploid individuals of the worked example of
`docs/specs/filters.md`, named `i1` to `i5`, with pop1 = i1, i2 and pop2 =
i3, i4, i5, `min_num_individuals` 1 and 4 bins, whose edges are 0, 0.25,
0.5, 0.75 and 1. The values are those of pyNei's
`_calc_obs_het_per_var`, `_calc_maf_per_var`,
`_calc_unbiased_exp_het_per_var` and `_calc_num_poly_vars` at commit
ef0ca6e on these genotypes, and the means and the histograms those of
`calc_per_var_distribs` of the same commit, run with
`unbiased_exp_het=False` for the plain expected heterozygosity and with
`unbiased_exp_het=True` for the unbiased one, the same with chunks of 6
and of 2 variants.

| variant | genotypes | obs het pop1 | obs het pop2 | maf pop1 | maf pop2 | unbiased exp het pop1 | unbiased exp het pop2 |
|---|---|---|---|---|---|---|---|
| 1 | 0/0 0/1 0/0 0/0 0/. | 1/2 = 0.5 | 0/2 = 0 | 3/4 = 0.75 | 5/5 = 1 | 0.5 | 0 |
| 2 | 0/0 0/1 0/0 ./. 0/. | 0.5 | 0/1 = 0 | 0.75 | 3/3 = 1 | 0.5 | 0 |
| 3 | 0/1 2/3 0/1 2/3 ./. | 2/2 = 1 | 2/2 = 1 | 1/4 = 0.25 | 0.25 | 1 | 1 |
| 4 | ./. ./. ./. ./. ./. | none | none | none | none | none | none |
| 5 | 0/0 0/0 0/0 0/0 1/1 | 0 | 0 | 1 | 4/6 = 0.666667 | 0 | 0.533333 |
| 6 | 0/. ./. ./. ./. ./. | none | none | none | none | none | none |

Variant 1 in pop2 has one called allele from its half called genotype,
which the maf and the expected heterozygosity count, 5 called alleles,
and which the observed heterozygosity does not, 2 called genotypes.
Variant 6 in pop1 has one called allele, half a genotype, below the
threshold of 1, so it has no maf and no expected heterozygosity; it has
no observed heterozygosity because it has no called genotype.

The means and the histogram counts, over the four bins:

| statistic | mean pop1 | mean pop2 | histogram pop1 | histogram pop2 |
|---|---|---|---|---|
| observed heterozygosity | 0.5 | 0.25 | 1, 0, 2, 1 | 3, 0, 0, 1 |
| maf | 0.6875 | 0.729167 | 0, 1, 0, 3 | 0, 1, 1, 2 |
| unbiased expected heterozygosity | 0.5 | 0.383333 | 1, 0, 2, 1 | 2, 0, 1, 1 |
| plain expected heterozygosity | 0.375 | 0.298611 | 1, 2, 0, 1 | 2, 1, 0, 1 |

The observed heterozygosity of 1 at variant 3 is in the last bin, and so
is the maf of 1 of variants 1, 2 and 5. With no `pops`, over the five
individuals, the mean observed heterozygosity is 0.395833 with the
histogram 1, 2, 0, 1, the mean maf 0.699008 with 0, 1, 0, 3, the mean
unbiased expected heterozygosity 0.430159 with 1, 2, 0, 1 and the mean
plain expected heterozygosity 0.378107 with 2, 1, 0, 1, each over
variants 1, 2, 3 and 5. The polymorphism ratio of the example is under
its item. These are checked at `calc_per_var_distribs` of "The Rust
interface", with blocks of 6 and of 2 variants, and the per variant
values at the function of each statistic named under its item.

The TypeScript test, under node, reads the panel and asserts the
literals of `s000`'s population, `p0`, that each statistic's item gives,
and the `Error` of an unknown individual.

## The observed heterozygosity

### What it gives

For one variant and one population, the share of its called genotypes
that are heterozygous, where a genotype is heterozygous when it is called
and its alleles are not all the same, at any ploidy:

    obs het = heterozygous genotypes / called genotypes

over the individuals of the population. A half called genotype is missing
and counts in neither. A population with no called genotype at a variant
has no value, and neither has one with fewer called genotypes than
`min_num_individuals`. At ploidy 1 no genotype is heterozygous and every variant
with a called genotype has 0. It is the number the filter of
`docs/specs/filters.md` compares with its threshold, there over every
individual.

### In Python

Asked for as `stats="obs_het"`, and given back as the field `obs_het` of
the result. It mirrors `_calc_obs_het_per_var` of `pynei/gt_counts.py`,
which takes no `min_num_samples`: in pyNei a population with one called
genotype, heterozygous, has an observed heterozygosity of 1 at that
variant, and it counts in the mean and the histogram, where the other
statistics give the variant no value; on two variants of five individuals
with one and with two called genotypes, run at commit ef0ca6e with the
default threshold, its mean observed heterozygosity is 0.75 over both and
the other means are NaN. The owner decided on 22 September 2026 that
popnei holds it to `min_num_individuals` like the others, since the
argument reads as one rule for the pass and the values it drops are of
one or two genotypes; the option not taken was to reproduce pyNei, which
changes the mean and the histogram of every dataset with variants below
the threshold. It compares its called genotypes, whole ones, where the
other statistics compare called alleles over the ploidy, so at a variant
with half called genotypes the two can part: four called genotypes and
two half called ones are 5 for the maf and 4 here, kept and dropped at a
threshold of 5. The comparisons against pyNei hold: on the panel every
population has 40 called genotypes or more at every variant, and on
`many.vcf` no population has fewer than 15 at any, above the 5 those
tests use.

### What pyNei asserts

`test_obs_het_stats` of `test/test_gt_counts.py` computes, on three
variants of four individuals with a half called genotype among them, the
values 1/3, 2/3 and none, in a line that compares them with no `assert`,
and asserts through `calc_per_var_distribs` with 4 bins the mean 0.5 and
the histogram 0, 1, 1, 0. popnei's test asserts the three values.

### How it is verified

Against plink2 v2.0.0-a.7.7, the arm64 build of 18 September 2026, on the
panel. Its `--hardy` report prints, per variant, `HOM_A1_CT`, `HET_A1_CT`
and `TWO_AX_CT`, the genotypes homozygous for the first allele,
heterozygous and homozygous for the other, and `O(HET_A1)`, the
heterozygous over the three, which for a biallelic variant is the
statistic. `--loop-cats` runs the report once per category of a column of
a file, which is how the populations are told apart:

    plink2 --vcf panel.vcf.gz --pheno panel_pops.txt --loop-cats popcat \
           --hardy --freq cols=+pos,+reffreq,+nobs --nonfounders --out panel

with `panel_pops.txt` holding `IID` and a `popcat` column, `p0`, `p1` or
`p2`, which writes `panel.p0.hardy` and `panel.p0.afreq` and the same for
the other two. `--hardy` and `--freq` report the founders only unless
`--nonfounders` is given, and every individual of a VCF is a founder, so
here the flag changes nothing. The same command without `--pheno` and
`--loop-cats` gives the reports over the 200 individuals, `panel.hardy`
and `panel.afreq`.

Against `_calc_obs_het_per_var` of pyNei, over the 1200 variants and the
3 populations and over all the individuals, the largest absolute
difference is 5.0e-7, half a unit of the last of the six digits plink2
prints, with no missing value on either side, so the tests compare within
1e-6 absolute. The literals of `var0000`, the first variant: in `p0`, 39
homozygotes for the reference, 8 heterozygotes and 1 homozygote for the
alternative, 0.166667; in `p1`, 15, 34 and 18, 0.507463; in `p2`, 81, 2
and 0, 0.0240964; over the 200 individuals, 135, 44 and 19, 0.222222,
with 2 genotypes missing.

The half called genotype and the third allele are checked against
bcftools 1.24 on `many.vcf`, which reads a half called genotype as pyNei
does and counts a heterozygous genotype whatever its alleles:

    bcftools +fill-tags many.vcf -- -t 'NHET=N_PASS(GT="het")','NCALLED=N_PASS(GT!="mis")'

writes into each variant the heterozygous genotypes and the called ones,
which `bcftools query -f '%POS\t%NHET\t%NCALLED\n'` prints and the
reference script keeps in `many.counts.tsv`. Over the 500
variants their quotient and pyNei's value are the same, with no missing
value on either side. The literals of the first three variants, at the
positions 1000, 1037 and 1074: 9 of 46, 24 of 47 and 32 of 47, which are
0.195652, 0.510638 and 0.680851 to plink2's digits.

The per variant values of the worked example and these literals are
checked at `ObsHet::of_var` of "The Rust interface", on the genotype
counts of the population.

## The major allele frequency

### What it gives

For one variant and one population, the frequency of its commonest
allele among the called alleles of the population. "maf" is, in pyNei and
in popnei, the frequency of the major allele, where most of the
literature and plink2 use the same letters for the minor one:

    maf = the largest of the counts of the alleles / called alleles

An allele is counted each time it was called, also in a half called
genotype. Every allele of a multiallelic variant has its own count, and
when two alleles tie for the largest count the maf is the same whichever
is called the major one. A population with no called allele has no value,
and neither has one with fewer called genotypes than
`min_num_individuals`, counted as the called alleles over the ploidy. A
variant with a maf of 1 does not vary in the population, and the
histogram of the maf is how a user sees how many of their variants hardly
vary.

### In Python

Asked for as `stats="maf"`, and given back as the field `maf` of the
result. It mirrors `_calc_maf_per_var` of `pynei/gt_counts.py`, whose
frequencies come from `_count_alleles_per_var` of the same file. There is
no difference besides those of the pass.

### What pyNei asserts

`test_maf_stats` of `test/test_gt_counts.py` asserts, on three variants of
four individuals with a half called genotype among them, the values 0.5,
4/7 and none, and through `calc_per_var_distribs` with 4 bins the mean
0.535714 and the histogram 0, 0, 2, 0; and on 3 random variants of 10
individuals with alleles 0 to 2 and missing ones, in two populations of
5 with `min_num_samples` 3, the six values, one of which is none.
`test_count_alleles_per_var` asserts the counts and the
frequencies of variants with the alleles 0, 1, 3 and 4 and with a
threshold of 3.1 genotypes.

### What pyNei does that is odd

Nothing that changes a value. `_calc_maf_per_var` counts the alleles up
to the largest of the whole chunk, as the expected heterozygosity's item
says, but a chunk with nothing called gives it no value for every
variant, run on two variants of three individuals with every genotype
`./.`, which is what popnei gives too.

### How it is verified

Against the `--freq` report of the plink2 command of the observed
heterozygosity, whose `REF_FREQ` and `ALT_FREQS` are the frequencies of
the two alleles of a variant among the `OBS_CT` called alleles, so the
larger of the two is the statistic where `OBS_CT` is 40 or more, 20
genotypes, and there is no value otherwise. Over the 1200 variants and
the 3 populations, and over all the individuals, no variant is below the
threshold and the largest absolute difference from `_calc_maf_per_var`
of pyNei is 5.0e-7, so the tests compare within 1e-6 absolute. The
literals of `var0000`: in `p0`, 96 called alleles and 0.895833; in `p1`,
134 and 0.522388, the alternative allele being the major one there; in
`p2`, 166 and 0.987952; over the 200 individuals, 396 and 0.792929.

The half called genotype and the third allele are checked against
bcftools 1.24 on `many.vcf`, per population, because bcftools counts a
half called genotype as pyNei does, missing as a genotype and its called
allele among the alleles:

    bcftools +fill-tags many.vcf -- -S many_pops.txt -t AN,AC

writes into each variant, for each population of `many_pops.txt`, a file
of individual name and population name, `AN_popA`, the called alleles,
and `AC_popA`, the count of each alternative allele, so the count of the
reference allele is the first less the sum of the second; the reference
script keeps them in `many.counts.tsv`, one line per variant. Over the 500
variants and the two populations, with `min_num_individuals` 0, the
largest of the counts over the called alleles is pyNei's value exactly.
The literals of the first three variants: in `popA`, 32 of 36, 27 of 39
and 21 of 40, which are 0.888889, 0.692308 and 0.525; in `popB`, 52 of 57,
35 of 56 and 24 of 55, which are 0.912281, 0.625 and 0.436364; over the
50 individuals, 84 of 93, 62 of 95 and 45 of 95. The counts are checked
exactly, at the allele counts of "The counts of one variant over a
population", and the frequencies within 1e-6 at `Maf::of_var` of "The
Rust interface", where the per variant values of the worked example are
checked too.

## The expected heterozygosity

### What it gives

For one variant and one population, the chance that gene copies taken at
random from that population are not all of the same allele: the genetic
diversity of that population at that site, 0 when every copy carries the
same allele and close to 1 when many alleles are at similar frequencies.
Over a dataset the user gets its mean over the variants and its
histogram.

The plain one is

    H = 1 - sum over a of p_a^k

where a runs over the alleles of the variant, p_a is the frequency of
allele a among the alleles of the population that were called at that
variant, and k is the ploidy. For k = 2 this is the classic 1 - p² - q²,
the chance that two copies differ. For any other ploidy it is the chance
that k copies are not all alike, and not the chance that two of them
differ.

The unbiased one corrects for the frequencies p_a being estimated from
the same copies the statistic is computed over, which makes H too small
on average, in the way that dividing a sum of squares by c rather than
by c - 1 makes a variance too small. With c the alleles of the population
that were called at the variant, the denominator of every p_a, and c_a
the count of allele a among them, it is

    H_u = 1 - sum over a of (c_a (c_a - 1) ... (c_a - k + 1)) / (c (c - 1) ... (c - k + 1))

with k factors in each product: the chance that k copies drawn from the
called ones without replacement are all alike, taken from 1. Its
expected value is the population's H at any k, which is what unbiased
means, and it is never above 1. At k = 2 it is (c / (c - 1)) H, Nei's
1978 correction, which GenAlEx, the population genetics add-in for Excel,
prints as the unbiased heterozygosity for codominant data, the name
pyNei's docstring gives it, and which pyNei computes at every ploidy.
The owner decided on 22 September 2026 that popnei applies the
correction of the ploidy in hand; the option not taken was pyNei's factor
2n/(2n - 1) with n = c / k at every ploidy, which on the tetraploid
variant under "What pyNei does that is odd" gives 1.1685 where this
gives 0.995960. There is no value when c is below k.

n, below, is c / k, how much called data there is counted in genotypes,
which can be a half when a genotype is half called.

Every allele of the variant counts. A multiallelic variant is not
collapsed to the major allele against the rest here, as it is in other
calculations of pyNei.

### In Python

Two statistics: the plain one, `PerVarStat.EXP_HET`, given back as the
field `exp_het` of the result, and the unbiased one,
`PerVarStat.UNBIASED_EXP_HET`, as `unbiased_exp_het`. One argument is
theirs: `ploidy=None`, which takes the ploidy of the variants when it is
not given. It mirrors
`_calc_exp_het_per_var` and `_calc_unbiased_exp_het_per_var` in
`pynei/diversity.py`, which get their frequencies from
`_count_alleles_per_var` in `pynei/gt_counts.py`. The order of the
populations is a difference, under "The populations"; the others are the
open points.

### Missing genotypes, populations with little data, and what pyNei asserts

The frequencies are over alleles, not over genotypes, so a half called
genotype contributes its called allele. A population of five diploid
individuals with four genotypes called and one half called has c = 9 and
n = 4.5, while the observed heterozygosity of the same population at the
same variant counts 4 called genotypes. The two statistics disagree about
how much data the population has.

A variant has no value in a population when n < `min_num_individuals`,
strictly less, so a population with exactly `min_num_individuals`
genotypes' worth of called data keeps its value. With no value the variant
is out of the numerator and the denominator of the mean and falls in no
histogram bin. A test on a small dataset has to lower
`min_num_individuals` or every value is missing; pyNei's own tests pass 1
or 5.

A population with no called allele at a variant has no value, also in a
block where no individual is called at any variant, where pyNei gives
one, under "What pyNei does that is odd".

A diploid population with one called allele, which needs
`min_num_individuals` at 0 to get this far, has c = 1, below k: pyNei
gives 0 for the plain one and NaN for the unbiased one, and popnei
`Some(0.0)` and `None`.

At ploidy 1 the plain one is 1 minus the sum of the frequencies, and so
is the unbiased one, from the formulas, which is 0 at every variant with
a called allele in exact arithmetic and a few 1e-16 away from 0 in
float64, where the frequencies are rounded. One haploid variant of nine
individuals with the alleles 0, 1, 1, 1, 1, 1, 2, 3 and 4 gives
-2.220446049250313e-16, below the range of the default histogram, so the
variant is in the mean and in no bin and a user of a haploid dataset
reads a mean a little below 0 beside a histogram that counts nothing.
numpy computes the same value from the same counts, so popnei does not
differ from pyNei here; popnei does not round the value to 0, which would
count the variant in the first bin where pyNei counts it in none.

`test_calc_exp_het` in `test/test_diversity.py` asserts the per variant
numbers of the worked example below, over the same two populations, and
the mean and the histogram of the default population over all five
individuals. `test_per_var_distribs_with_threads` in `test/test_threads.py`
asserts that the mean and the histogram counts do not change with 2 or 4
threads. No test of pyNei varies the size of the chunks for this
statistic, and the chunk with nothing called below is a case where its
mean would change with them. popnei adds a test that the size of the
block does not change the result, which holds because popnei gives such
a block no value.

### What pyNei does that is odd

The three runs below are of pyNei at commit ef0ca6e, each one call of the
named function from a python session on a chunk of
`Variants.from_gt_array` written for it.

`_count_alleles_per_var` takes the alleles it counts from the largest
allele in the whole chunk, `gts.max()`, because asking numpy which alleles
are there cost more than the counting. A chunk is pyNei's unit of work, a
few thousand variants, the block of section 2 of the architecture. When no
genotype in the chunk is called at all, the largest allele is the missing
one and there are no alleles to count: the sum over them is 0 and every
variant of that chunk gets a plain expected heterozygosity of 1 - 0 = 1
and an unbiased one of -0.0. Measured on two variants of three
individuals with every genotype `./.`: plain 1.0 and 1.0, unbiased -0.0
and -0.0, and through `calc_per_var_distribs` a mean of 0.0 with both
variants counted in the first histogram bin. The same variant in a chunk
where any individual is called anywhere gives no value. The owner decided
on 22 September 2026 that popnei gives no value, which is what a
calculation that works row by row does by itself; the option not taken,
to reproduce pyNei, needs the alleles of a whole block, which nothing
else here needs, and makes the value of a variant depend on where the
block boundaries fell. The comparison against pyNei is safe: the panel
drops 3 in 100 genotypes of 200 individuals and no block of it is
without called genotypes.

The unbiased factor is 2n/(2n-1) whatever the ploidy, in
`_calc_unbiased_exp_het_per_var`, although at ploidy k the population
holds c = kn copies and not 2n. Measured with `min_num_samples=1` on three
tetraploid individuals, genotypes 0/1/2/3, 0/0/1/1 and 0/1/2/3: the
unbiased value is 1.1685, above 1, so it is in no bin of the default
histogram while the mean counts it. popnei applies the correction of the
ploidy in hand, decided under "What it gives", and gives 0.995960 for
that variant, which is 1 - 48/11880: allele counts 4, 4, 2 and 2 of 12,
and the products of four factors 24, 24, 0 and 0 over 12 · 11 · 10 · 9.

The `ploidy` argument, when it differs from the ploidy of the data, is
used as the exponent k and as the number of alleles the individuals are
expected to hold, `len(pop) * ploidy - missing_alleles`, while the
`min_num_samples` test keeps the ploidy of the data, which is the one
`_count_alleles_per_var` uses. Measured with `ploidy=4` on four diploid
individuals with allele counts 5 and 3: a plain 0.8276 and an unbiased
0.9459, which is 0.8276 · 8/7, from n = 16/4 = 4. The 8 alleles the
individuals really hold, as 4 genotypes of the data's ploidy, give the
same n and the same 0.9459. The two part when an allele is missing: with
counts 4 and 3, pyNei, counting the alleles the individuals are expected
to hold, 16 less the missing one, gives an unbiased 0.9919 from n = 15/4;
the 7 alleles that were called, as 3.5 genotypes of the data's ploidy,
give 1.0029. The owner decided on 22 September 2026 that popnei keeps the
argument as the exponent k alone, of the frequencies and of the products
of the unbiased one, and takes the called alleles from the data, which is
the same as pyNei whenever the argument matches the data and is what the
code that works row by row gives for free. The options not taken were to
reproduce the three uses, which takes deliberate work because the count
of called alleles comes from the row, and to drop the argument. The
comparison against pyNei is made at the data's own ploidy.

### How it runs

Per variant, the value follows from the allele counts of the population
alone. pyNei's chunk wide set of alleles has no counterpart here, because
the sum runs over the alleles the row itself has.

### How it is verified

Against the `--hardy` report of the plink2 command of the observed
heterozygosity, which prints, per variant, `E(HET_A1)`, the frequency of
heterozygotes expected under Hardy Weinberg, which for a biallelic variant
is 1 - p² - q², the plain statistic.

Against `_calc_exp_het_per_var(chunk, pops, min_num_samples=20)` of pyNei
on the panel, over the 1200 variants and the 3 populations, the largest
absolute difference is 5.0e-7 in each population, with no missing value
on either side, half a unit of the last of the six digits plink2 prints,
so the tests compare within 1e-6 absolute, one unit of that digit. The
numbers written into the first of those tests as literals are those of
variant `var0000`: in p0, 39 homozygotes for the reference, 8
heterozygotes, 1 homozygote for the alternative and an expected
heterozygosity of 0.186632; in p1, 15, 34, 18 and 0.498998; in p2, 81, 2,
0 and 0.0238061; over the 200 individuals, 0.328385.

No program outside the project prints the unbiased one. What is checked
instead is that the factor rebuilds it from plink2's own numbers: c is
twice the genotypes it counted, and for `var0000` in p1, c = 134 and
0.498998 · 134/133 = 0.502750, 4.5e-7 from pyNei's 0.5027494. The half
called genotype rule has no outside check either: plink2 refuses a VCF
with a half call unless `--vcf-half-call` tells it what to do with one,
and the panel has none.

The worked example, which becomes the first cargo test: 3 variants, 5
diploid individuals, pop1 = i1, i2 and pop2 = i3, i4, i5,
`min_num_individuals` 1.

| variant | genotypes | plain pop1 | plain pop2 | unbiased pop1 | unbiased pop2 |
|---|---|---|---|---|---|
| 1 | 0/0 2/1 0/0 0/0 0/. | 0.625 | 0 | 0.833333 | 0 |
| 2 | 0/0 0/0 0/1 1/0 ./. | 0 | 0.5 | 0 | 0.666667 |
| 3 | ./. ./. ./. ./. ./. | none | none | none | none |

Variant 1 in pop2 has c = 5, four alleles from two genotypes and one from
the half called one, all of allele 0, so both values are 0. The means
over the variants with a value are 0.3125 and 0.25 plain, 0.416667 and
0.333333 unbiased. The per variant values are the ones
`test_calc_exp_het` asserts; the means come from a run of pyNei at commit
ef0ca6e with these two populations.

The per variant values of this table, of the worked example of the pass
and the plink2 literals of `var0000` are checked at `ExpHet::of_var`, on
the counts the row helper gives for each population. No function above it
shows the value of one variant: `calc_per_var_distribs` gives the mean
and the histogram. The means are checked through that function, in the
cargo tests of the pass and, against pyNei, in the pytest test of the
pass.

## The polymorphism ratio

### What it gives

For one population, how many of the variants vary in it, in three counts
and two ratios over the variants that had a value, so that a user sees
how much of their dataset is informative for that population. A variant
is polymorphic in a population when its major allele frequency there is
below `poly_threshold`, strictly, 0.95 by default, inherited from pyNei's
`DEF_POLY_THRESHOLD`; nobody has measured whether 0.95 is the right
threshold. It is variable when that frequency is below 1. The counts are
`num_poly`, the polymorphic variants; `num_variable`, the variable ones;
and `tot_num_variants_with_data`, the variants that have a maf in the
population, those with `min_num_individuals` called genotypes or more.
The ratios are `poly_ratio`, the first over the third, and
`poly_ratio_over_variables`, the first over the second, and each is NaN
when its denominator is 0.

### In Python

Asked for as `stats="poly_vars_ratio"`, and given back as the field
`poly_vars_ratio` of the result, a `PolyVarsStats`, pyNei's frozen
dataclass of five pandas series over the populations, the three counts as
integers and the two ratios as float64. `poly_threshold` is its argument.
It mirrors `_calc_num_poly_vars` of `pynei/diversity.py`, whose maf comes
from `_calc_maf_per_var`, so the same variants have a value here and in
the maf. There is no difference besides those of the pass.

### What pyNei asserts

`test_poly_vars_ratio` of `test/test_diversity.py` asserts, on 3 random
variants of 10 individuals with alleles 0 to 2 and missing ones, two
populations of 5, `min_num_samples` 3 and a threshold of 0.51, the counts
3 and 2 variable and 1 and 1 polymorphic, and through
`calc_per_var_distribs` with chunks of 2 the ratios 1/3 and 1/2, both
over the variants with data and over the variable ones.

### How it is verified

Against the `--freq` reports of the maf on the panel: with the maf of a
variant as that item takes it from `REF_FREQ`, `ALT_FREQS` and `OBS_CT`,
the number of variants whose maf is below 0.95, below 1, and that have
one are, in `p0`, 1112, 1173 and 1200; in `p1`, 1101, 1177 and 1200; in
`p2`, 1093, 1184 and 1200; and pyNei's `_calc_num_poly_vars` with
`min_num_samples` 20 gives the same nine counts. The ratios follow: in
`p0`, 1112/1200 = 0.926666666667 over the variants with data and
1112/1173 = 0.947996589940 over the variable ones. On `many.vcf`, with
the two populations and the counts of bcftools of the maf,
`min_num_individuals` 5: in `popA`, 477, 492 and 500; in `popB`, 478, 493
and 500, the same in pyNei, and the ratios of `popA` 0.954 and
477/492 = 0.969512195122. The counts are compared exactly and the ratios
within 1e-12 against those digits, at `calc_per_var_distribs` of "The
Rust interface".

The worked example of the pass, at the same function: pop1 has the mafs
0.75, 0.75, 0.25 and 1, so 3 polymorphic, 3 variable and 4 with data, and
the ratios 0.75 and 1; pop2 has 1, 1, 0.25 and 0.666667, so 2, 2 and 4, and
the ratios 0.5 and 1. With no `pops`, over the five individuals, the
mafs are 0.888889, 0.857143, 0.25 and 0.8, so 4 polymorphic, 4 variable
and 4 with data, and both ratios 1. A `poly_threshold` of 0.5 over the
two populations leaves 1 polymorphic variant in each, the variant 3 of
maf 0.25, and changes neither the variable ones nor those with data, so
the ratios are 0.25 in both populations over the variants with data and
0.333333 and 0.5 over the variable ones. With `min_num_individuals` 20
every count is 0 and both ratios NaN.

## The per individual statistics

### What they give

For each individual, the share of the variants at which its genotype is
missing, and the share at which it is heterozygous:

    missing rate = missing genotypes of the individual / variants
    obs het rate = heterozygous genotypes of the individual / called genotypes of the individual

over the variants the pass gives, after the steps. A half called genotype
is missing and not heterozygous. The first tells a user which individuals
were badly genotyped, and the second which ones are more heterozygous
than the rest, a sign of a mixed sample or of an outcrossed individual
among inbred ones. An individual with no called genotype has a missing
rate of 1 and no heterozygosity rate, NaN. pyNei divides the second by
every variant, those at which the individual has no genotype among them,
so an individual with more missing data looks less heterozygous: on the
panel, `s000` has 426 heterozygous genotypes of 1200 variants, 0.355, and
of 1166 called ones, 0.365352, and over the 200 individuals the two
differ by up to 0.0163. The owner decided on 22 September 2026 that
popnei divides by the called genotypes, which is what a user reads the
number as and what plink2's `--het` gives; the option not taken was
pyNei's denominator, and the missing rate beside the number says what
pyNei's one number said.

### In Python and in TypeScript

```python
def calc_per_individual_stats(variants: Variants) -> PerIndividualStats
```

It is a consumer of the `Variants`, with one pass, as
`calc_per_var_distribs` is. `PerIndividualStats` is a frozen dataclass
with `missing_gt_rate` and `obs_het_rate`, two pandas series of float64
indexed by the names of the individuals, in the order of the pass, and
`pass_stats`. It mirrors `calc_per_sample_stats` of
`pynei/sample_stats.py`. The differences:

- The name: pyNei's sample is popnei's individual, as `docs/glossary.md`
  has it.
- The result is a dataclass and not pyNei's pandas frame of the two
  columns, because every result of a consumer carries its `pass_stats`,
  which a frame has no place for. The two series are the frame's columns
  under their names.
- The heterozygosity rate divides by the called genotypes of the
  individual, under "What they give"; pyNei divides by all the variants.
- There is no `num_threads`, as in every calculation of popnei.
- A pass that gives no variant is a `ValueError` with the message of the
  pass of `calc_per_var_distribs`. pyNei fails with a `TypeError` from
  inside its pipeline, run at commit ef0ca6e on a `Variants` whose filter
  kept nothing.

In TypeScript it is `calcPerIndividualStats(variants)`, which gives
`individuals`, an array of names, `missingGtRate` and `obsHetRate`, two
`Float64Array` in that order, and `passStats`.

### What pyNei asserts

`test_filter_missing` of `test/test_sample_stats.py` asserts, on three
variants of five individuals with chunks of 2, the missing rates 1/3,
1/3, 1/3, 1/3 and 2/3 and the heterozygosity rates 0, 1/3, 1/3, 1/3 and 0;
the third variant has every genotype missing. popnei's test asserts the
same missing rates and, over the called genotypes, the heterozygosity
rates 0, 1/2, 1/2, 1/2 and 0.
`test_per_sample_stats_with_threads` asserts that 2 or 4 threads change
nothing.

### How it runs

Over the rows of each block, with rayon across them, with no `reblock`
before it. Each thread keeps two counts per individual, missing and
heterozygous, and adds to them for each row; the counts of the threads
are added when the block is done. What is kept from one block to the
next is those two arrays of integers and the number of variants, which
grows with the individuals and not with the variants. The two divisions
are made once, at the end, so the result is the same to the bit whatever
the size of the blocks and the number of threads. It asks its reader for
the genotypes alone.

### How it is verified

Against plink2 v2.0.0-a.7.7 on the panel:

    plink2 --vcf panel.vcf.gz --missing --sample-counts cols=+hom,+het,+missing \
           --het --nonfounders --out panel

`--missing` writes `panel.smiss`, with `MISSING_CT`, `OBS_CT` and
`F_MISS` per individual, the missing genotypes, the variants and their
quotient, which is the missing rate; `--sample-counts` writes
`panel.scount`, with `HET_CT`, the heterozygous genotypes, which over
the called ones, `OBS_CT` less `MISSING_CT` of the first, is the
heterozygosity rate; `--het` writes `panel.het`, whose `OBS_CT` is that
same count of called genotypes. Over the 200 individuals, the largest
absolute difference from pyNei's `calc_per_sample_stats` is 3.3e-8 in the
missing rate, which plink2 prints with six digits, so the tests compare
it within 1e-6; pyNei's heterozygosity rate over all the variants is
`HET_CT` over 1200 exactly, and popnei's over the called ones is the
quotient of two counts of plink2, compared within 1e-12. The literals:
`s000` has 34 missing genotypes of 1200, 0.0283333, and 426 heterozygous
of 1166 called, 0.365352; `s001` 44, 0.0366667, and 397 of 1156,
0.343426.

The half called genotype and the third allele are checked against the
same command on `many.vcf` with `--vcf-half-call m`, which reads a half
called genotype as missing, as pyNei does, and writes `many.smiss` and
`many.scount`; the third allele of a heterozygous genotype changes
nothing in `HET_CT`. Over the 50 individuals the missing rate is pyNei's
exactly, and the heterozygous counts are pyNei's rate times 500. The
literals: `ind00` has 29 missing genotypes of 500, 0.058, and 201
heterozygous of 471 called, 0.426752; `ind01` 25, 0.05, and 195 of 475,
0.410526. The counts are the numbers bcftools 1.24 prints too, as
`nMissing` and `nHets` of the `PSC` lines of `bcftools stats -s -
many.vcf`.

Against pyNei: both libraries run the function on the panel and on
`many.vcf`; the missing rates have to be equal within 1e-12 relative,
with the same names in the same order, and popnei's heterozygosity rate
has to be pyNei's times the variants over the called genotypes of the
individual, within 1e-12, the called genotypes being the variants less
pyNei's missing rate times them.

The worked example, at `calc_per_individual_stats` of "The Rust
interface", on the six variants of the pass with blocks of 6 and of 2:
the missing rates 2/6, 2/6, 2/6, 3/6 and 5/6, and the heterozygosity
rates 1/4, 3/4, 1/4, 1/3 and 0/1 = 0. Individual i5 has two half called
genotypes, at variants 1 and 2, which are missing.

## The Rust interface

The populations of a pass: each one's name and the indices of its
individuals among those of the reader, in the order the user gave them.
Its errors are new cases of the error of the crate, each a `ValueError`
in Python: a name that is not an individual, with the population and the
name; a name twice in one population, with both; a population with no
individual, with its name; and no population.

```rust
pub struct Pops { /* private */ }
impl Pops {
    /// One population, named "pop", of every individual of a reader of
    /// `num_individuals`, in the order of the reader.
    pub fn all(num_individuals: usize) -> Pops;
    /// The populations named by the user, each individual looked up
    /// among `individuals`, the ones the pass gives, with
    /// `resolve_individuals` of `docs/specs/filters.md` for each
    /// population, which refuses an unknown name, a name twice and no
    /// name.
    pub fn from_names(pops: &[(String, Vec<String>)], individuals: &[String])
        -> Result<Pops>;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
    pub fn name(&self, pop: usize) -> &str;
    /// Indices into the individuals of the reader.
    pub fn individuals(&self, pop: usize) -> &[usize];
    /// Whether the population is every individual in the order of the
    /// reader, for a caller that reads the row as it is then.
    pub fn is_all(&self, pop: usize) -> bool;
}
```

The counts of one variant over the individuals of one population, in the
`variant` module beside `count_gts` and `count_alleles` of
`docs/specs/variant.md`, whose rules and errors they share. `gts` is the
genotypes of a `VariantRef`, or a row of the genotypes of a block, and
`individuals` the indices of the population. An index at or beyond the
individuals of the row is an error, a defect of popnei and not of the
user, whose names were resolved before.

```rust
pub fn count_gts_of(gts: &[i8], ploidy: usize, individuals: &[usize]) -> Result<GtCounts>;
pub fn count_alleles_of(gts: &[i8], ploidy: usize, individuals: &[usize],
                        counts: &mut AlleleCounts) -> Result<u32>;
```

The bins of the histogram. The edges are computed as "In Python and in
TypeScript" of the pass says, and a value is placed by them.

`MAX_NUM_BINS`, 100000, is the most bins a histogram is built with, and a
`num_bins` above it is refused like a `num_bins` of 0, a `ValueError` in
Python. A histogram a person reads has tens of bins, pyNei's default is
40, and each bin is a count of 8 bytes for every population and every
statistic, held once by the pass and once more by each chunk of rows a
thread is reading, so 100000 bins of the four statistics of one population
are 3.2 MB per chunk. Above the bound the counts are a vector no machine
gives: 2^60 bins are the `PanicException` of a capacity that overflowed,
which derives from `BaseException`, so `except Exception` does not catch
it and a notebook dies, and 1e12 bins abort the interpreter where the
allocation fails.

```rust
pub struct HistBins { /* private */ }
impl HistBins {
    /// `num_bins` equal widths from `start` to `end`. An error when
    /// `num_bins` is 0 or above `MAX_NUM_BINS`, when `start` is not below
    /// `end`, when either is not a number, and when the distance between
    /// the two is above the largest float64, which leaves the edges NaN
    /// and infinite instead of going up.
    pub fn linear(start: f64, end: f64, num_bins: usize) -> Result<HistBins>;
    /// `num_bins` equal ratios. As above, and an error when `start` is
    /// 0 or below.
    pub fn logarithmic(start: f64, end: f64, num_bins: usize) -> Result<HistBins>;
    pub fn edges(&self) -> &[f64];
    pub fn num_bins(&self) -> usize;
    /// The bin of `value`, or None outside the range. The last bin takes
    /// its right edge.
    pub fn bin_of(&self, value: f64) -> Option<usize>;
}
```

How each statistic of a variant is worked out for one population, from
the counts of that population, and whether the variant has a value there.
`None` counts for nothing: the variant is out of the mean and in no bin.
The ploidies are `usize`, as they are in the architecture. The fields are
private behind a constructor that refuses a ploidy or an exponent of 0 or
above 255, the largest ploidy `open_vcf` takes, with an error that is a
`ValueError` in Python. The owner decided it on 22 September 2026, after
a trial implementation of the expected heterozygosity in which an
exponent of 0 gave a plain -2.0 and an unbiased NaN for counts 2, 1, 1, a
ploidy of 0 turned the `min_num_individuals` test off, and an exponent
of 1e8 took 0.2 s for one variant; the options not taken were public
fields that the caller guarantees, and `NonZeroUsize` fields.

```rust
pub struct ObsHet { /* private */ }
impl ObsHet {
    /// `min_num_individuals` is how many called genotypes a population
    /// needs at a variant to get a value.
    pub fn new(min_num_individuals: u32) -> ObsHet;
    /// het over called, from the genotype counts of the population. None
    /// when called is 0 or below `min_num_individuals`.
    pub fn of_var(&self, counts: GtCounts) -> Option<f64>;
}

pub struct Maf { /* private */ }
impl Maf {
    /// An error for a ploidy of 0 or above 255.
    pub fn new(ploidy: usize, min_num_individuals: u32) -> Result<Maf>;
    /// The largest of `counts` over `called_alleles`, which is their sum.
    /// None when the population has called fewer than
    /// `min_num_individuals` genotypes, `called_alleles` below
    /// `min_num_individuals` times the ploidy, and when it has called
    /// nothing.
    pub fn of_var(&self, counts: &AlleleCounts, called_alleles: u32) -> Option<f64>;
}

pub struct ExpHet { /* private */ }
impl ExpHet {
    /// `exponent` is k, the exponent of the frequencies and the number of
    /// factors of the products of the unbiased one: the ploidy of the
    /// variants unless the caller asks for another one. `ploidy` is the
    /// ploidy of the variants, which turns the alleles a pop has called
    /// into called genotypes, n, for the min_num_individuals test; the
    /// exponent is never used for n. An error for either of 0 or above
    /// 255.
    pub fn new(exponent: usize, ploidy: usize, min_num_individuals: u32) -> Result<ExpHet>;
    /// The expected heterozygosity of one variant in one pop, the plain
    /// one or, with `unbiased`, the unbiased one. `counts[a]` is how often
    /// allele a was called in the pop at this variant, and
    /// `called_alleles` is their sum; a missing allele is in neither. None in
    /// three cases: the pop has called fewer than `min_num_individuals`
    /// genotypes; it has called nothing at all at this variant; or the
    /// unbiased one was asked for and `called_alleles` is below the
    /// exponent, which for a diploid pop is one allele.
    pub fn of_var(&self, counts: &AlleleCounts, called_alleles: u32,
                  unbiased: bool) -> Option<f64>;
}
```

`counts` is what `count_alleles_of` fills for the individuals of one
population, and the caller hands it the same buffer for every variant.
The three `of_var` compare `called_alleles` with `min_num_individuals`
times the ploidy rather than dividing as pyNei does, which is the same
test because the ploidy is positive, and it never has to hold the 4.5
genotypes of a population with a half called one.

What one pass calculates, and what it gives back. The reader is the
outermost reader of the chain that `chain_of` of `docs/specs/filters.md`
built from the steps of the `Variants`, lent as the writer of the vars
file, `write_vars` of `docs/specs/io_vars.md`, and the Kosman distances,
`calc_kosman_sums` of `docs/specs/dists.md`, take theirs: when the
calculation returns, the binding crate reads the
counts of the filters from the chain and the number of variants from the
result, and together they are the `pass_stats`. The `pops` are built by
the binding crate with `Pops::from_names` against `individuals()` of that
reader, or with `Pops::all`. The errors are those of the reader, of
`HistBins`, and, a new case of the error of the crate, a pass that gave no
variant, a `ValueError` in Python whose message says whether the source
had none or the steps kept none, with the counts of each filter.

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PerVarStat { ObsHet, Maf, ExpHet, UnbiasedExpHet, PolyVarsRatio }

pub struct PerVarDistribsConfig {
    pub stats: Vec<PerVarStat>,
    pub pops: Pops,
    pub bins: HistBins,
    pub obs_het: ObsHet,
    pub maf: Maf,
    pub exp_het: ExpHet,
    /// Below this maf a variant is polymorphic. Not a number, below 0 or
    /// above 1 is an error of `calc_per_var_distribs`.
    pub poly_threshold: f64,
}

/// The distribution of one statistic, per population.
pub struct StatsDistrib { /* private */ }
impl StatsDistrib {
    pub fn num_pops(&self) -> usize;
    /// The sum over the count. None when no variant had a value.
    pub fn mean(&self, pop: usize) -> Option<f64>;
    /// How many variants had a value.
    pub fn num_vars_with_value(&self, pop: usize) -> u64;
    /// One count per bin.
    pub fn hist_counts(&self, pop: usize) -> &[u64];
    pub fn bins(&self) -> &HistBins;
}

/// The three counts of the polymorphism ratio, per population.
pub struct PolyVarsStats { /* private */ }
impl PolyVarsStats {
    pub fn num_poly(&self, pop: usize) -> u64;
    pub fn num_variable(&self, pop: usize) -> u64;
    pub fn num_vars_with_data(&self, pop: usize) -> u64;
    /// num_poly over num_vars_with_data. None when the latter is 0.
    pub fn poly_ratio(&self, pop: usize) -> Option<f64>;
    /// num_poly over num_variable. None when the latter is 0.
    pub fn poly_ratio_over_variables(&self, pop: usize) -> Option<f64>;
}

pub struct PerVarDistribs {
    pub obs_het: Option<StatsDistrib>,
    pub maf: Option<StatsDistrib>,
    pub exp_het: Option<StatsDistrib>,
    pub unbiased_exp_het: Option<StatsDistrib>,
    pub poly_vars_ratio: Option<PolyVarsStats>,
    /// The variants the pass gave.
    pub num_vars: u64,
}

pub fn calc_per_var_distribs<R: BlockReader + ?Sized>(
    reader: &mut R, config: &PerVarDistribsConfig,
) -> Result<PerVarDistribs>;
```

The per individual statistics, over the same kind of reader, with the
same error for a pass that gave no variant.

```rust
pub struct PerIndividualStats { /* private */ }
impl PerIndividualStats {
    pub fn num_individuals(&self) -> usize;
    pub fn num_vars(&self) -> u64;
    pub fn num_missing(&self, individual: usize) -> u64;
    pub fn num_het(&self, individual: usize) -> u64;
    /// num_missing over num_vars.
    pub fn missing_rate(&self, individual: usize) -> f64;
    /// num_het over the called genotypes, num_vars less num_missing.
    /// None when the individual has no called genotype.
    pub fn obs_het_rate(&self, individual: usize) -> Option<f64>;
}

pub fn calc_per_individual_stats<R: BlockReader + ?Sized>(reader: &mut R)
    -> Result<PerIndividualStats>;
```

The binding crates build the Python and TypeScript results from these,
NaN for a `None`, and the names of the individuals from `individuals()`
of the reader.

## Speed

The per variant passes are among what the objectives want fast. The
dataset is the vars file popnei writes of the 400 MB VCF of
`docs/rust_core.md`, 100000 variants of 1000 diploid individuals, 3 in
100 genotypes missing whole, `big.vars` of `docs/reports/filters-measurement.md`,
in the page cache, and the machine the owner's Apple M5 Pro, 18 cores,
macOS 27.0. popnei's pass over that file with no filter and the genotypes
alone asked for takes 0.102 s on one thread and on 18, measured on 21
September 2026 in `docs/specs/filters.md`, and the missing data filter,
which counts the genotypes of every row, adds 0.056 to 0.060 s to it on
one thread and 0.013 to 0.018 s on 18.

pyNei at commit ef0ca6e on the same variants, from the vars file of its
own format written from that VCF, on 22 September 2026, best of 3 runs
at a load average of 2.3 to 3.6, numpy 2.5.3 on Accelerate, the BLAS of
macOS that numpy calls, which the counting of these statistics does not
use:

| | 1 thread | 6 threads |
|---|---|---|
| its pass, the chunks alone | 0.187 s | |
| `calc_per_var_distribs`, the four statistics, no `pops` | 1.13 s | 0.263 s |
| the same, `obs_het` alone | 0.992 s | 0.242 s |
| the same, `maf` alone | 0.185 s | 0.193 s |
| the same, the four, 4 populations of 250 | 2.67 s | |
| `calc_per_sample_stats` | 1.02 s | 0.255 s |

The numbers to reach, for the five statistics with no `pops`, which is
what pyNei's four cost it and one more value from the same counts, and
for the per individual statistics, each a whole pass over `big.vars`:
0.25 s on
one thread and 0.15 s on 18 cores. Each is the pass plus twice what the
filter adds to it, rounded up, because the statistics read each row
twice, once for the genotype counts and once for the allele counts,
where the filter reads it once. The first measurement of the
implementation says whether they hold, and what the four populations
cost over that, which has no number yet. They are 4.5 times under
pyNei's 1.13 s on one thread, and under its best with 6 threads. In
wasm, single threaded, nothing has been measured, and the number to
reach is set with the first measurement, as `docs/specs/dists.md` set
its.

## Open points

None. The owner decided the seven this spec had on 22 September 2026,
and each is written where it applies, with the option that was not
taken: a block with nothing called gives no value, under "What pyNei does
that is odd" of the expected heterozygosity; the unbiased correction is
the one of the ploidy in hand, under "What it gives" of the same; the
`ploidy` argument is the exponent alone, under "What pyNei does that is
odd"; the three statistics of a variant are built by constructors that
refuse a 0, under "The Rust interface"; `min_num_individuals` holds the
observed heterozygosity too, under its item; the heterozygosity rate of
an individual is over its called genotypes, under "What they give" of the
per individual statistics; and the bin type of equal widths is spelt
`"linear"`, under "In Python and in TypeScript" of the pass.

## Not in this spec

- The filter of individuals, pyNei's `filter_samples`: its item of
  `docs/specs/filters.md`, written with this spec.
- The expected heterozygosities inside Jost's D, the differentiation
  between two populations: `_calc_pairwise_dest` of `pynei/dists.py`
  builds its own 1 - sum p^k per population and over the pooled
  frequencies, and corrects both with the harmonic mean of the called
  genotypes of the two populations, which is not the correction here.
  `docs/specs/dists.md`, which takes its `pops` from here.
- The linkage disequilibrium per population, which takes one pass per
  population: the `ld` spec, which is not written.
- Fis, Fst and the other diversity summaries, and the per variant values
  themselves, as a frame: pyNei does not have them and popnei does not
  add them. A user who wants the values takes the genotypes with
  `iter_blocks`.
- The three threshold filters over the individuals of one population:
  pyNei does not have it, and `docs/specs/filters.md` says so.
- The read ahead thread of section 3 of the architecture: a reader over a
  reader, which `docs/specs/dists.md` leaves to `docs/specs/block.md`, and
  which both functions here take like any other reader.
