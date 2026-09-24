# The diversity module: how much variety each population holds

24 September 2026. A user of popnei gives a dataset and a set of
populations and gets, for each population, how many alleles it holds, how
many of those no other population called, how many of the variants vary in
it, how its variants are spread over the frequency of their rarer allele,
and how far its genotypes are from the proportions its allele frequencies
would give. The first three come both as they stand and standardized to a
common number of called alleles, so that a population of 20 individuals
and one of 200 can be compared. There is no code. This spec adds a row to
section 9 of `docs/architecture.md`, `diversity`, and depends on three
others: `docs/specs/stats.md` for what a population is and how a `pops`
argument is read, `docs/specs/variant.md` for the helper that counts how
often each allele was called over a set of individuals, and
`docs/specs/block.md` and `docs/specs/filters.md` for the reader the pass
runs over.

The option not taken was to add these five as items of `docs/specs/stats.md`,
whose row of the architecture already holds the allele counts per
population and the polymorphism ratio. A module of its own was chosen on
24 September 2026 because the results have a shape the `stats` results
have not: `calc_per_var_distribs` gives a mean and a histogram over the
variants, and four of the five here give one number per population and the
fifth a vector over the frequency of the rarer allele. The cost is a
second module that counts alleles per population, the same counting
`stats` does, and a line of `docs/specs/stats.md` that has to change:
"Fis, Fst and the other diversity summaries ... pyNei does not have them
and popnei does not add them", under "Not in this spec" there, which the
owner reversed on 24 September 2026 after a web application asked for
them.

**pyNei has none of the five.** At commit ef0ca6e, `pynei/diversity.py`
holds the expected heterozygosity and the polymorphism ratio and nothing
else, and none of private, rarefaction, spectrum, folded or inbreeding
appears in `src/` or in `test/` except one comment of `pynei/gwas.py`. So
no item below has a pyNei function to mirror, none has a difference from
pyNei to write down, and no test runs the two libraries on the same input.
Every number of this spec comes from a program outside the project, and
the one quantity no outside program computes is marked where it appears.

## The pass and its function

### What it gives

One pass over the variants serves every population and every statistic,
because all five are built from the same thing: the allele counts of one
population at one variant, how often each allele was called among the
individuals of that population, and `c`, the called alleles of the
population there, which is their sum. A half called genotype gives its
called allele to those counts, as `docs/glossary.md` has it.

Four of the five are also given at a common number of called alleles, and
that is one operation: draw `g` of the `c` called alleles of the
population at that variant, without replacement, and take the expectation
over every such draw. Population genetics has two names for it. Applied to
counts of alleles it is called **rarefaction**, and it is what makes the
number of alleles of two populations comparable, since a population of
more individuals finds more alleles for no reason but its size. Applied to
the site frequency spectrum it is called **projection**. popnei computes
both from the one argument `num_called_alleles`, which is that `g`.

### Its Python function

```python
calc_pop_diversity(variants: Variants,
                   pops: dict[str, Sequence[str]] | None = None,
                   stats: Iterable[PopDiversityStat] = tuple(PopDiversityStat),
                   num_called_alleles: int | None = None,
                   min_num_individuals: int = 20) -> PopDiversity
```

It is a consumer of the `variants`: one pass over the source through the
steps the `Variants` has when it is called, and the `Variants` is as it
was afterwards.

`pops` is read as `docs/specs/stats.md` reads it, and everything that item
says holds here: a dict of population name to the names of its
individuals, `None` being one population with every individual of the
dataset named `"pop"`, a name that no individual has being a `ValueError`
that says which, a population with no individual being a `ValueError`, an
individual allowed in more than one population and one in none read by
none. The order of the populations in every result is the order of the
keys of `pops`.

`stats` names the statistics to compute, as `calc_per_var_distribs` does,
and one not asked for is `None` in the result. The members of
`PopDiversityStat` are `NUM_ALLELES`, `PRIVATE_ALLELES`, `POLY_RATIO`,
`FOLDED_SFS` and `FIS`, and each is the name of the field that holds its
result.

`min_num_individuals` is how many called genotypes a population needs at a
variant for the variant to count for it, 20 by default, and it is measured
as the called alleles of the population over the ploidy, which is the rule
of `docs/specs/stats.md` and lets a half called genotype count as half.
Strictly fewer than that and the variant does not count for that
population; a population with exactly `min_num_individuals` keeps it.

`num_called_alleles` is `g`, the called alleles every population is
brought down to. `None`, the default, means the counts are taken over the
called alleles each population happens to have, and then the three
standardized values are NaN and `FOLDED_SFS` cannot be asked for, since
the bins of a spectrum need one number of alleles for every population and
every variant; asking for it without `num_called_alleles` is a
`ValueError`. A `num_called_alleles` below 2 is a `ValueError`: a draw of
one allele shows one allele whatever the population holds. A variant whose
population called fewer than `g` alleles does not count for the
standardized values of that population, and that is why the result carries
a second count of variants.

`PopDiversity` is a frozen dataclass:

```python
@dataclass(frozen=True)
class PopDiversity:
    pops: tuple[str, ...]
    num_alleles: pandas.DataFrame | None        # total, mean, in_draw
    private_alleles: pandas.DataFrame | None    # total, mean, in_draw
    poly_vars: pandas.DataFrame | None          # total, ratio, in_draw
    folded_sfs: pandas.DataFrame | None
    fis: pandas.Series | None
    num_vars: pandas.DataFrame                  # with_data, in_draw
    num_vars_every_pop: int
    num_vars_every_pop_in_draw: int
    pass_stats: PassStats
```

The first three frames are indexed by population name, in the order of
`pops`, with the three columns named in the comment: `total` and `mean`
are over the called alleles the population has, and `in_draw` is the
standardized value. `folded_sfs` has one row per count of the rarer
allele, 0 to `num_called_alleles // 2`, indexed by that count, and one
column per population. `fis` is one value per population. `num_vars` says,
per population, how many variants had enough called genotypes and how many
of those also reached `num_called_alleles`; `num_vars_every_pop` and
`num_vars_every_pop_in_draw` are the two counts of the variants where
every population did, which are the divisors of the private alleles.
`pass_stats` is the counts of the pass, as every consumer of popnei gives
them.

In TypeScript it is `calcPopDiversity(variants, {pops, stats,
numCalledAlleles, minNumIndividuals})`, with `pops` an object of
population name to an array of individual names. Where Python gives a
frame of three columns TypeScript gives an object of three
`Float64Array`s, `{total, mean, inDraw}`, in the order of `pops`; the
spectrum is one `Float64Array` per population, keyed by name; `fis` is a
`Float64Array`; and the counts are `Float64Array`s beside the two whole
numbers. The names of the populations are in `pops`, a frozen array of
strings.

### The cases

A population for which no variant counted has 0 in every count, NaN in
every mean and ratio, NaN in `fis`, and a column of zeros in the spectrum.
It is not an error: a user filtering hard enough to leave one population
with nothing still wants the others.

A pass that gives no variant is a `ValueError`, as `docs/specs/dists.md`
has it, whose message says whether the source held none or the steps kept
none.

A `num_called_alleles` above every population's called alleles leaves
every `in_draw` value NaN and every spectrum zero, with
`num_vars.in_draw` at 0. It is not an error, and the counts say why the
values are missing.

With one population every allele it called is private, since there is no
other population to hold it. The number is then the number of alleles, and
the spec says so rather than refusing the call: a user asking for the
diversity of one population has no reason to be stopped, and
`private_alleles` equal to `num_alleles` is what the arithmetic gives.

### How it runs

One pass over the blocks of the reader. For each block, rayon over its
rows; for each row, the allele counts of every population from the helper
of `docs/specs/variant.md`, and from them the contribution of that variant
to every statistic asked for. What is kept from one block to the next does
not grow with the variants: per population, six sums and four counts, plus
`num_called_alleles // 2 + 1` sums for the spectrum. The counts of the
populations at one variant are held while that variant is worked on,
because the private alleles need every population's counts at once, and
they are a `num_pops` by `num_alleles` array of the row and not of the
block.

Nothing here needs `reblock`, and nothing holds more than one block.

## The number of alleles

### What it gives

How many different alleles a population called. At one variant it is the
alleles whose count is above 0, from 1 for a population where every copy
is alike to the alleles the variant has. Over a dataset a user gets it two
ways, because the two answer different questions: the **total**, the sum
over the variants that counted, which says how much variety the whole
dataset holds in that population, and the **mean** over those variants,
which is the allelic richness, the number a user compares between
populations. An allele numbered 1 at one variant is not the allele
numbered 1 at the next, so the total is a sum of per variant counts and
never a count of distinct things across the dataset.

The standardized value is the mean over the variants of

    E = sum over a of [ 1 - C(c - n_a, g) / C(c, g) ]

where `a` runs over the alleles the population called at the variant,
`n_a` is how often it called allele `a`, `c` is its called alleles there,
`g` is `num_called_alleles`, and `C(x, y)` is the number of ways of
choosing `y` of `x`. Each term is the chance that allele `a` appears at
least once in a draw of `g` of the `c` copies, so the sum is the alleles
the draw is expected to show. It is the rarefaction of El Mousadik and
Petit (1996), and the same arithmetic Hurlbert (1971) gave for the species
expected in a sample of a community. The mean is over the variants whose
population called at least `g` alleles; a variant with fewer has no draw
of `g` and is left out.

### How it is verified

The totals against `adegenet` 2.1.11, whose `genind2genpop` gives a table
of one column per variant and allele holding the times each population
called it; the alleles above 0 in a row are the total of that population.
On the panel of `docs/specs/stats.md`, `tests/reference/stats/panel.vcf.gz`,
1200 biallelic diploid variants of 200 individuals with 3 in 100 genotypes
missing, whose populations `p0`, `p1` and `p2` of 48, 68 and 84
individuals are `tests/reference/stats/panel_pops_bcftools.txt`, run on
24 September 2026 with `min_num_individuals` 20:

| population | total | mean | in a draw of 20 |
|---|---|---|---|
| p0 | 2373 | 1.9775 | 1.9283948650 |
| p1 | 2377 | 1.9808333333 | 1.9219209943 |
| p2 | 2384 | 1.9866666667 | 1.9197370844 |

`adegenet` gives 2373, 2377 and 2384, compared exactly since they are
counts. The means follow from the totals over the 1200 variants, which all
three populations kept at that threshold.

The standardized values against `vegan` 2.7.6, whose `rarefy(counts,
sample = g)` is the formula above, called on the allele counts of each
variant and averaged over the variants; it gives the three numbers of the
last column to their ten digits. They are compared within 1e-12 relative,
since both sides sum the same per variant values in different orders.
`vegan` is the reference here and not `hierfstat`, whose
`allelic.richness` computes the same quantity: `hierfstat` does not
install on the owner's machine, `RcppParallel` and then its dependency
`gaston` failing to build on 24 September 2026, which
`docs/specs/dists.md` recorded already.

The worked example, a cargo test at `calc_pop_diversity` of "The Rust
interface", is the six variants of five diploid individuals of
`docs/specs/filters.md`, with `pop1` of `i1` and `i2` and `pop2` of `i3`,
`i4` and `i5`, the two populations of the worked example of
`docs/specs/stats.md`, `min_num_individuals` 1 and `num_called_alleles` 4.
Four variants count for each population: variant 4 has nothing called and
variant 6 has one called allele in `pop1`, half a genotype, below the
threshold of 1. `pop1` calls 2, 2, 4 and 1 alleles at variants 1, 2, 3 and
5, a total of 9 and a mean of 2.25; `pop2` calls 1, 1, 4 and 2, a total of
8 and a mean of 2. In a draw of 4, `pop1` keeps all four of its variants
and gets 2.25, and `pop2` keeps three of them, variant 2 having 3 called
alleles, and gets 2.3111111111.

## The private alleles

### What it gives

The alleles a population called that no other population of the call
called at all: what that population holds and the others do not. At one
variant it is a count from 0 to the alleles the population called, and
over a dataset a user gets the **total** over the variants and the
**mean** over them, as for the number of alleles above and for the same
reason.

A variant counts for the private alleles only when **every** population
has enough called genotypes there. A population that called nothing at a
variant holds none of its alleles, so every allele of every other
population would be private there, and the count would measure the missing
data and not the populations. That is why the result carries
`num_vars_every_pop`, the divisor of the mean, beside the per population
counts.

The standardized value is the mean over the variants of

    E = sum over a of [ q_p(a) * product over r != p of (1 - q_r(a)) ]

where `q_r(a) = 1 - C(c_r - n_r(a), g) / C(c_r, g)` is the chance that
allele `a` appears in a draw of `g` of the called alleles of population
`r`, `p` is the population whose private alleles are counted, and `a` runs
over the alleles `p` called. Each term is the chance that the allele shows
in `p`'s draw and in no other population's draw. It is the estimator of
Kalinowski (2004), which the program ADZE computes. The mean is over the
variants where every population called at least `g` alleles, which is
`num_vars_every_pop_in_draw`, and that count can be below the one for the
number of alleles, since one short population takes the variant from every
population's private count.

### How it is verified

The totals against `poppr` 2.9.8, whose `private_alleles(genind,
count.alleles = FALSE)` gives a table of one column per variant and allele
with a 1 where the allele is private to the population of the row; its row
sums are the totals. On the panel, run on 24 September 2026 with
`min_num_individuals` 20, `poppr` gives 0, 0 and 1 for `p0`, `p1` and
`p2`, which popnei has to match exactly; the means over the 1200 variants
every population kept are 0, 0 and 0.0008333333.

The standardized values have no program on this machine that computes
them. ADZE is the one that does and it is not installed (**Open 2**). They
are checked by the worked example alone, whose numbers are worked out by
hand below, and by two properties that a cargo test asserts on the panel:
with `num_called_alleles` equal to the smallest `c` of the dataset the
value is at most the standardized number of alleles of the same
population, since a private allele is an allele; and a population compared
against a copy of itself has 0 private alleles at every draw size, since
every allele it draws the copy can draw too. On the panel at a draw of 20
the values are 0.0112196177, 0.0099715392 and 0.0089014974, computed by
the reference script of this module and stored beside it.

The worked example, the same six variants and two populations as above at
`num_called_alleles` 4. Three variants have both populations at 4 called
alleles or more: variant 1, where `pop1` counts allele 0 three times and
allele 1 once of 4 and `pop2` counts allele 0 five times of 5; variant 3,
where both count alleles 0, 1, 2 and 3 once of 4; and variant 5, where
`pop1` counts allele 0 four times of 4 and `pop2` counts allele 0 four
times and allele 1 twice of 6. At variant 1 `pop1` draws its 4 of 4, so
allele 1 appears for certain, and `pop2` draws 4 of 5 without ever showing
allele 1, which it does not hold: the term is 1. At variant 3 every allele
of `pop1` appears in `pop2`'s draw of 4 of 4 for certain, so every term is
0, and at variant 5 allele 0 is certain in both. `pop1` sums to 1 over
three variants, 0.3333333333. For `pop2`, allele 1 of variant 5 appears in
its draw of 4 of 6 with chance 1 - C(4, 4) / C(6, 4) = 14/15 and `pop1`
cannot draw it, so `pop2` sums to 0.9333333333 over three variants,
0.3111111111. Without a draw, over the four variants every population has
data at, `pop1` has the total 2, allele 1 at variants 1 and 2, and the
mean 0.5, and `pop2` has the total 1, allele 1 at variant 5, and the mean
0.25.

## The polymorphism ratio

### What it gives

How many of the variants vary in a population. A variant varies in a
population when that population called more than one allele there, which
is the `num_variable` of `docs/specs/stats.md` and not its `num_poly`,
whose threshold on the major allele frequency has no meaning for a draw of
alleles. A user gets the **total**, the varying variants, and the
**ratio**, those over the variants that counted for the population, which
is the same number `calc_per_var_distribs` gives as
`poly_ratio_over_variables` would give at a threshold of 1.

The standardized value is the mean over the variants of

    P = 1 - sum over a of C(n_a, g) / C(c, g)

the chance that a draw of `g` of the `c` called alleles is not all of one
allele, since `C(n_a, g) / C(c, g)` is the chance that every copy drawn is
allele `a`. It is what ADZE reports as the proportion of polymorphic loci,
and it is the quantity that makes the polymorphism of two populations of
different size comparable: a population of 200 individuals finds a rare
allele that a population of 20 misses, and without the draw it looks more
polymorphic for that reason alone.

### How it is verified

On a variant of two alleles the expected number of alleles in a draw of
`g` is one plus the chance that the draw varies, because the draw shows
either one allele or two. So on a dataset where every variant has two
alleles `vegan::rarefy` verifies this item as well: the standardized ratio
has to be the standardized number of alleles minus 1, to the bit. The
identity was checked on 24 September 2026 on the counts (10, 6), (3, 1),
(17, 3), (1, 19) and (55, 45) at draws of 2, 4 and 10, where the two sides
agree to 1e-12, and on the panel, whose standardized ratios are
0.9283948650, 0.9219209943 and 0.9197370844 against the standardized
allele counts of 1.9283948650, 1.9219209943 and 1.9197370844. A dataset
with variants of more than two alleles has no such identity and no
program here to check it; the worked example covers it.

The totals and ratios on the panel, run with `min_num_individuals` 20 on
24 September 2026: 1173, 1177 and 1184 varying variants of 1200, ratios of
0.9775, 0.9808333333 and 0.9866666667. The totals are compared exactly.
They are also the variants `docs/specs/stats.md` counts as variable, so
the two modules have to agree on them, and a pytest test asserts that
against `calc_per_var_distribs` with `poly_threshold` 1 on the same
dataset and the same `pops`.

The worked example, at `num_called_alleles` 4. `pop1` varies at variants
1, 2 and 3 and not at 5, a total of 3 of 4 and a ratio of 0.75; `pop2`
varies at 3 and 5 and not at 1 and 2, 2 of 4 and 0.5. In a draw of 4,
`pop1` keeps its four variants and every draw of 4 of 4 shows what the
population holds, so its standardized ratio is its ratio, 0.75; `pop2`
keeps three variants and gets 0.6444444444, its variant 5, with 4 copies
of allele 0 and 2 of allele 1 of 6, varying in a draw of 4 with chance
1 - C(4, 4) / C(6, 4) = 14/15.

## The folded site frequency spectrum

### What it gives

How the variants of a population are spread over the count of their rarer
allele: a vector whose entry `j` is how many variants show `j` copies of
the rarer allele in a draw of `g` called alleles. It is the shape a
population's history leaves in its variants, an excess of rare alleles
after growth and a flatter curve after a bottleneck, and it is the input
of every program that fits a demographic model.

It is **folded**: the entry `j` and the entry `g - j` are one, because
without an outgroup nothing says which allele is the ancestral one, only
which is rarer. The bins run from 0 to `g // 2`. Bin 0 holds the variants
that show one allele in the draw, which the draw makes possible even for a
variant that varies in the population, and at an even `g` the last bin,
`g / 2`, is not doubled, since `j` and `g - j` are the same bin there.

A variant of more than two alleles is treated as if it had two: the major
allele of the population at that variant, the most frequent among its
called alleles and the lowest numbered of two that tie, as
`docs/glossary.md` defines it, against every other allele summed into one
rarer allele. The owner decided this on 24 September 2026; the option not
taken was to leave multiallelic variants out of the spectrum, which loses
them from a dataset where one variant in ten has three alleles, as
`tests/reference/vcf/many.vcf` has.

The projection is per variant: a variant whose population called `c`
alleles, `m` of them not the major one, contributes to bin `min(j, g - j)`
the chance

    C(m, j) * C(c - m, g - j) / C(c, g)

for every `j` from 0 to `m`, which is the hypergeometric chance that a
draw of `g` of the `c` copies holds `j` of the rarer ones. A variant whose
population called fewer than `g` alleles contributes nothing. So the
entries are not whole numbers, and they sum to the variants that counted,
`num_vars.in_draw` of that population.

### How it is verified

Against `dadi` 2.4.4, whose `Spectrum.from_data_dict(dd, [pop],
projections=[g], polarized=False)` does the same per variant projection,
drops the variants below `g` called alleles, and folds. On the panel with
`min_num_individuals` 20 and `num_called_alleles` 20, run on 24 September
2026, all 33 entries of the three populations agree to their ten printed
digits:

| rarer allele | p0 | p1 | p2 |
|---|---|---|---|
| 0 | 85.9261619951 | 93.6948068115 | 96.3154987275 |
| 1 | 92.9965194033 | 95.4588032593 | 101.3781441620 |
| 2 | 106.8896328231 | 103.7262343009 | 108.0515292765 |
| 3 | 115.5454376764 | 110.6831031479 | 114.8500079366 |
| 4 | 120.5055938623 | 116.6729455230 | 119.7397786175 |
| 5 | 122.9411934960 | 121.0517724083 | 121.8833113836 |
| 6 | 124.0815127004 | 123.6021904937 | 121.8645089676 |
| 7 | 124.2393741502 | 124.5944395085 | 120.6009754676 |
| 8 | 123.4960764903 | 124.5403261788 | 118.9768135703 |
| 9 | 122.4216218019 | 124.0682314877 | 117.7147383254 |
| 10 | 60.9568756012 | 61.9071468804 | 58.6246935652 |

Each column sums to 1200, the variants that counted. They are compared
within 1e-12 relative. `dadi` masks bin 0 and the bins above `g / 2` in a
folded spectrum and popnei reports bin 0, so the comparison reads
`fs.data` and not the masked array; the difference is one of presentation
and the value is the same.

`dadi` does not build on the project's Python, 3.14 with the free
threading build, its `nlopt` dependency failing to compile on 24 September
2026; it installs and runs on 3.12, so the reference script of this module
makes an environment of its own (**Open 3**).

The worked example, at `num_called_alleles` 4. `pop2` keeps variants 1, 3
and 5. Variant 1 has 5 copies of allele 0 and nothing else, so every draw
of 4 gives 0 rarer copies and it is one whole variant in bin 0. Variant 3
has one copy each of alleles 0, 1, 2 and 3: allele 0 is the major one by
the tie rule, the rarer allele has 3 copies of the 4, and a draw of 4 of 4
gives 3, which folds to bin 1. Variant 5 has 4 copies of allele 0 and 2 of
allele 1 of 6, and a draw of 4 of 6 gives 0 rarer copies with chance
1/15, 1 with 8/15 and 2 with 6/15, which fall in bins 0, 1 and 2. So
`pop2` is 1.0666666667, 1.5333333333 and 0.4, summing to its 3 variants.
`pop1` keeps its four variants and every draw is of 4 of 4: variants 1, 2
and 3 give 1, 1 and 3 rarer copies, which fold to bin 1, and variant 5
gives 0, so `pop1` is 1, 3, 0.

## The inbreeding coefficient F_IS

### What it gives

How far the genotypes of a population are from the proportions its allele
frequencies would give if its individuals paired at random: 0 when they
match, positive when the population holds fewer heterozygous genotypes
than that, which inbreeding, selfing or a population split into unmixed
groups all produce, and negative when it holds more. It is

    F_IS = 1 - (mean over the variants of H_o) / (mean over the variants of H_u)

with `H_o` the observed heterozygosity of the population at a variant and
`H_u` its unbiased expected heterozygosity, both as `docs/specs/stats.md`
defines them, over the variants that counted for the population. It is
Nei's F_IS.

Two choices are in that line and both were measured on 24 September 2026,
by drawing genotypes at a known F, an individual being heterozygous with
chance `2 p q (1 - F)` and homozygous with `p² + F p q` and `q² + F p q`,
over 2000 variants whose allele frequencies come from the neutral
spectrum, where the variants at frequency x go as 1/x, 50 datasets a row.

The first is the ratio of the two means against the mean of the per
variant ratios, `1 - H_o / H_u` averaged over the variants. The mean of
ratios is pulled by the variants with little diversity, where the ratio is
large and noisy. With 30 individuals and the rarest allele at 0.01, at a
true F of 0.5, the mean of ratios gives 0.4398, 0.060 below the truth,
and the ratio of means gives 0.4944, 0.0056 below it; at a true F of 0.2
the two are 0.1680 and 0.1979 against a truth of 0.2. Over all twelve
rows of the measurement, 30 and 100 individuals, the rarest allele at 0.01
and 0.05, a true F of 0, 0.2 and 0.5, the ratio of means is never more
than 0.0056 from the truth and the mean of ratios is up to 0.0602 from it.
The standard deviations of the two are within 0.0012 of each other. The
ratio of means is what popnei computes.

The second is the unbiased expected heterozygosity against the plain one.
The plain one is too small by about `1 / (c - 1)`, so it makes F_IS too
small by about the same: at a true F of 0, the plain one gives -0.0177
with 30 individuals and -0.0053 with 100, where the unbiased one gives
-0.0007 and -0.0002. popnei uses the unbiased one. `scikit-allel` uses the
plain one, so a user comparing the two libraries sees a difference, about
0.006 on the panel, and this spec records it rather than following
`scikit-allel`.

F_IS is the one statistic here that the draw of `g` called alleles does
not touch: the observed heterozygosity is a property of whole genotypes
and not of a sample of alleles. `num_called_alleles` is ignored by it.

At ploidy 1 no genotype can be heterozygous, so `H_o` is 0 at every
variant and F_IS would be 1 wherever the population has any diversity.
popnei gives NaN for a population of a haploid dataset rather than that 1,
and the reason is in the doc comment of the field.

### How it is verified

Against `scikit-allel` 1.3.13, which gives `heterozygosity_observed` and
`heterozygosity_expected` per variant and builds the plain form of this
statistic from them in `inbreeding_coefficient`. The reference script runs
it on each population of the panel and takes the ratio of the two means,
which is popnei's formula with the plain expected heterozygosity, and the
check is made at that form: on 24 September 2026 it gives -0.0237536998,
-0.0258924700 and -0.0247472838 for `p0`, `p1` and `p2`, compared within
1e-12 relative. The unbiased form, which is what the function returns, is
-0.0127584868, -0.0181107131 and -0.0184585832 on the same data; it
differs from the plain one by the `c / (c - 1)` factor of each variant's
expected heterozygosity, whose own verification against pyNei and plink2
is in `docs/specs/stats.md`, so no program is run twice for it.

`adegenet` 2.1.11 is the other program with the ingredients, and it is not
the reference: its `Hs` is the plain expected heterozygosity with no
correction for the sample, 0.3481518, 0.3477311 and 0.3427402 on the three
populations of the panel on 24 September 2026, checked by hand on a
dataset of five individuals and two populations where it gives 0.5 and
0.4375 for allele counts of 3 and 3 of 6 and of 3 and 1 of 4. It would
check the same half of the statistic `scikit-allel` checks, through an
extra dependency.

The worked example. `pop1` has observed heterozygosities of 0.5, 0.5, 1
and 0 at variants 1, 2, 3 and 5 and unbiased expected heterozygosities of
0.5, 0.5, 1 and 0, both means 0.5, so its F_IS is 0. `pop2` has 0, 0, 1
and 0 against 0, 0, 1 and 0.5333333333, means of 0.25 and 0.3833333333, so
its F_IS is 0.3478260870. The eight per variant values are the ones the
worked example of `docs/specs/stats.md` tabulates, so the two specs have
to agree on them.

## The Rust interface

The statistics to compute, a bit set as `Needs` is:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DiversityStats(u8);

impl DiversityStats {
    pub const NUM_ALLELES: DiversityStats;
    pub const PRIVATE_ALLELES: DiversityStats;
    pub const POLY_RATIO: DiversityStats;
    pub const FOLDED_SFS: DiversityStats;
    pub const FIS: DiversityStats;
    pub const ALL: DiversityStats;
    pub fn contains(self, other: DiversityStats) -> bool;
}
```

What the pass is asked for:

```rust
pub struct DiversityOptions {
    pub stats: DiversityStats,
    /// The called alleles every population is brought down to. `None` gives
    /// no standardized value and no spectrum.
    pub num_called_alleles: Option<u32>,
    /// The called genotypes a population needs at a variant, its called
    /// alleles over the ploidy.
    pub min_num_individuals: u32,
}
```

The pass, one over the reader, every population and every statistic at
once. The populations are the indices of their individuals, in the order
the user gave them, and an empty slice of populations is one population
with every individual.

```rust
/// # Errors
///
/// `FOLDED_SFS` asked for with no `num_called_alleles`, a
/// `num_called_alleles` below 2, a population with no individual, an index
/// that is not an individual of the dataset, an individual asked for more
/// than once, no variant in the reader, a variant of more alleles than a
/// count of them holds, and those of the reader.
pub fn calc_pop_diversity<R: BlockReader + ?Sized>(
    reader: &mut R,
    pops: &[&[usize]],
    options: &DiversityOptions,
) -> Result<PopDiversity>;

pub struct PopDiversity { /* private */ }

impl PopDiversity {
    pub fn num_pops(&self) -> usize;
    /// The variants that had `min_num_individuals` called genotypes in the
    /// population. `None` when `pop` is not a population of the call.
    pub fn num_vars(&self, pop: usize) -> Option<u64>;
    /// Of those, the ones whose called alleles reached
    /// `num_called_alleles`.
    pub fn num_vars_in_draw(&self, pop: usize) -> Option<u64>;
    /// The variants where every population had enough called genotypes,
    /// and of those the ones where every population reached
    /// `num_called_alleles`. They are the divisors of the private alleles.
    pub fn num_vars_every_pop(&self) -> u64;
    pub fn num_vars_every_pop_in_draw(&self) -> u64;
    /// The alleles the population called, summed over its variants. `None`
    /// when `pop` is not a population of the call or the statistic was not
    /// asked for.
    pub fn num_alleles(&self, pop: usize) -> Option<u64>;
    /// Their mean over a draw of `num_called_alleles`, NaN when no variant
    /// reached it.
    pub fn num_alleles_in_draw(&self, pop: usize) -> Option<f64>;
    pub fn private_alleles(&self, pop: usize) -> Option<u64>;
    pub fn private_alleles_in_draw(&self, pop: usize) -> Option<f64>;
    /// The variants where the population called more than one allele.
    pub fn num_poly_vars(&self, pop: usize) -> Option<u64>;
    pub fn poly_ratio_in_draw(&self, pop: usize) -> Option<f64>;
    /// One value per count of the rarer allele, 0 to
    /// `num_called_alleles / 2`.
    pub fn folded_sfs(&self, pop: usize) -> Option<&[f64]>;
    pub fn fis(&self, pop: usize) -> Option<f64>;
}
```

The chance that a draw of `g` of `c` misses an allele called `n` times,
`C(c - n, g) / C(c, g)`, is the one piece of arithmetic every standardized
value rests on. It is computed as the product of `g` factors,
`(c - n - i) / (c - i)` for `i` from 0 to `g - 1`, each in `f64`, and not
from three factorials, which overflow an `f64` above 170 and lose digits
long before that. The product is taken in one direction for every call so
that two populations of one pass round the same way.

This module adds cases to the error of the crate, and the two binding
crates divide them as `docs/specs/pca.md` divides its own. All of them are
the wrong input of a function and a `ValueError` in Python: the spectrum
asked for without a draw size, a draw size below 2, a population with no
individual, an individual the dataset has not or named twice, a pass with
no variant, and a variant of more alleles than a count of them holds.

## Speed

No measurement has been made and the first task of the plan is to make
one. The pass reads every genotype once and does, per variant and per
population, work that grows with the alleles the population called, except
for the spectrum, which does `g + 1` hypergeometric terms per variant and
per population and is the one part that can dominate: a dataset of 1000000
variants and 3 populations at a draw of 200 called alleles is 6·10⁸ terms,
each a product of about 200 factors if computed as above. The measurement
comes first, on the panel and on a dataset of 100000 variants and 1000
individuals, and it decides whether the per variant hypergeometric weights
are computed once per variant and shared by the bins, which they can be,
since the weights of one variant differ from those of the next only
through `c` and `m`.

The target is set when that measurement exists, against
`calc_per_var_distribs` of `docs/specs/stats.md`, which makes the same
pass and the same allele counts and does far less arithmetic on them.

## Open points

The owner decides these. Until then the implementer follows the
"meanwhile" of each.

**Open 1: one draw size or a curve of them.** `num_called_alleles` takes
one number, so a user who wants the curve of allelic richness against the
number of alleles drawn, which is what ADZE prints and what shows whether
a population has been sampled enough, makes one call per point. The
options are one number, which is the sketch above; a sequence of numbers,
which turns the three standardized values into a row per draw size and the
spectrum into one per draw size too, and makes the result two dimensional
where it is now one; or one number with the curve left to the caller's
loop. What the sequence saves is the passes: a curve of ten points is ten
passes over the variants against one, and a pass over a million variants
is the cost of this whole module. What it costs is a result shape that
every user pays for and most do not use, and a spectrum that is a frame
per draw size. Recommendation: one number, and the sequence added if a
user asks for the curve. Meanwhile one number is built.

**Open 2: no program checks the rarefied private alleles.** Every other
value of this spec is checked against a program outside the project, which
the objectives ask for. The standardized private alleles are checked by
the worked example and by two properties of the panel, all three of them
arithmetic of popnei's own. The options are to build ADZE, the program of
Szpiech, Jakobsson and Rosenberg (2008) that computes Kalinowski's
estimator, from its source, which is C and not in any package manager on
this machine, and add it to the reference script; or to leave the value
checked by hand. What building ADZE costs is a source build in the
repository's reference tooling, which no other reference needs, and a
program that has had no release since 2014. What leaving it costs is one
value of five whose formula nothing outside popnei confirms.
Recommendation: leave it checked by hand for the first version, and build
ADZE if a user reports a number they doubt. Meanwhile the worked example
and the two properties stand.

**Open 3: `dadi` needs a second Python.** The projected spectrum is
checked against `dadi` 2.4.4, which does not build on the project's Python
3.14 free threading build and does build on 3.12. Every other reference of
popnei is a program run by a shell script or a library in the one
environment. The options are a second environment made by the reference
script, pinned to 3.12, which `uv` creates in about 5 seconds; or to drop
`dadi` and check the projection against a worked example alone, as the
private alleles are checked. What the second environment costs is a second
Python in the reference tooling and a note in the script saying why. What
dropping it costs is the only outside check of the projection, which is
the piece of arithmetic in this spec most likely to be got wrong, since it
is a distribution and not a single value. Recommendation: the second
environment. Meanwhile the reference script makes it.

## Not in this spec

- F_ST, f_2, Jost's D and the other measures between two populations:
  `docs/specs/pop_dists.md`, which computes them from the same allele
  counts per population in its own pass.
- The unfolded spectrum, which needs an ancestral allele that popnei does
  not read, and the joint spectrum of two populations, which is a matrix
  of one axis per population and which the demographic fits use. Neither
  has a user yet.
- The per variant values of any of these as a frame. A user who wants them
  takes the genotypes with `iter_blocks`, as `docs/specs/stats.md` says
  for its own.
- The expected heterozygosity, the observed heterozygosity, the major
  allele frequency and the polymorphism ratio with a threshold:
  `docs/specs/stats.md`. F_IS above is built from the first two and the
  item here for the polymorphism ratio counts the same variants its
  `num_variable` counts.
- Rarefaction of the heterozygosities. They are frequencies and not
  counts, so they do not grow with the individuals sampled the way the
  number of alleles does, and no program here rarefies them.
- The read ahead thread of section 3 of `docs/architecture.md`: a reader
  over a reader, which `docs/specs/block.md` holds, and which the pass
  here takes like any other reader.
