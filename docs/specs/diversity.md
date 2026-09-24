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
population and the polymorphism ratio. The owner chose a module of its own
on 24 September 2026 because the results have a shape the `stats` results
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
Every number of this spec but one comes from a program outside the
project; the standardized private alleles, which no program here computes,
are checked by enumerating every draw, under "How it is verified" of that
item.

The five programs those numbers come from are `vegan` 2.7.6 and
`adegenet` 2.1.11, which were on the owner's machine already, `poppr`
2.9.8, installed from CRAN on 24 September 2026, and `scikit-allel` 1.3.13
and `dadi` 2.4.4 from PyPI. The owner took them as development
dependencies of popnei on 24 September 2026, so the plan that builds this
module adds them where the other development dependencies are and the
reference script runs them. `hierfstat`, which would have checked the
standardized allele counts and F_IS in one package, is not among them: it
does not install, `RcppParallel` and then its dependency `gaston` failing
to build, which `docs/specs/dists.md` recorded on 23 September 2026 and
which is still so.

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
`PopDiversityStat` are `NUM_ALLELES`, `PRIVATE_ALLELES`,
`VARIABLE_VARS_RATIO`, `FOLDED_SFS` and `FIS`, and each is the name of the
field that holds its result.

`min_num_individuals` is how many called genotypes a population needs at a
variant for the variant to count for it, 20 by default, and it is measured
as the called alleles of the population over the ploidy, which is the rule
of `docs/specs/stats.md` and lets a half called genotype count as half.
Strictly fewer than that and the variant does not count for that
population; a population with exactly `min_num_individuals` keeps it. A
population that called nothing at a variant does not count it whatever
`min_num_individuals` is, so a `min_num_individuals` of 0 does not put
variants with no data into the totals.

`num_called_alleles` is `g`, the called alleles every population is
brought down to. `None`, the default, means the counts are taken over the
called alleles each population happens to have, and then the three
standardized values are NaN and `FOLDED_SFS` cannot be asked for, since
the bins of a spectrum need one number of alleles for every population and
every variant; asking for it without `num_called_alleles` is a
`ValueError`. A `num_called_alleles` below 2 is a `ValueError`: a draw of
one allele shows one allele whatever the population holds.

It takes one number and not a sequence of them, so a user who wants the
curve of allelic richness against the number of alleles drawn, which shows
whether a population has been sampled enough, calls the function once per
point of it. The owner decided this on 24 September 2026. The option not
taken was a sequence, which would turn the three standardized values into
a row per draw size and give a spectrum per draw size, and would save the
passes, a curve of ten points being ten passes over the variants against
one; what it costs is a result of two dimensions that every user carries
and most do not use. A sequence can be added later without changing what
one number gives.

A variant is **in the draw** for a population when it counts for that
population by the rule above **and** the population called at least `g`
alleles there, both conditions and not the second alone. The standardized
number of alleles and the standardized ratio of variable variants, and the
spectrum, are over the variants in the draw for that population. The
standardized private alleles are over the variants in the draw for every
population, as the totals are over the variants that count for every
population. That is why the result carries a second count of variants
beside the first, and a second count for all the populations together.

`PopDiversity` is a frozen dataclass:

```python
@dataclass(frozen=True)
class PopDiversity:
    pops: tuple[str, ...]
    num_alleles: pandas.DataFrame | None        # total, mean, in_draw
    private_alleles: pandas.DataFrame | None    # total, mean, in_draw
    variable_vars_ratio: pandas.DataFrame | None  # total, ratio, in_draw
    folded_sfs: pandas.DataFrame | None
    fis: pandas.Series | None
    num_vars: pandas.DataFrame                  # with_data, in_draw
    num_vars_every_pop: int
    num_vars_every_pop_in_draw: int
    pass_stats: PassStats
```

The first three frames are indexed by population name, in the order of
`pops`, with the three columns named in the comment: `total` and, for the
two counts of alleles, `mean`, and for the variable variants `ratio`, are
over the called alleles the population has, and `in_draw` is the
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
frame of three columns TypeScript gives an object of three typed arrays in
the order of `pops`: `{total, mean, inDraw}` for the two counts of alleles
and `{total, ratio, inDraw}` for the variable variants, with `total` a
`Uint32Array`, as `docs/specs/stats.md` gives its counts, and the other
two `Float64Array`s. The spectrum is one `Float64Array` per population,
keyed by name; `fis` is a `Float64Array`; `numVars` is two `Uint32Array`s,
`{withData, inDraw}`, beside the two whole numbers for all the populations
together. The names of the populations are in `pops`, a frozen array of
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
values are missing. The totals, the means and `fis` are as they would be
without the argument, since none of them reads it.

With one population every allele it called is private, since there is no
other population to hold it. The number is then the number of alleles, and
the spec says so rather than refusing the call: a user asking for the
diversity of one population has no reason to be stopped, and
`private_alleles` equal to `num_alleles` is what the arithmetic gives.

Two populations that share an individual make the standardized private
alleles a number to read with care. The formula multiplies the chance that
an allele is in one population's draw by the chance that it is in no
other's, which is right when the draws are of different copies and not
when they are partly of the same ones. Two populations holding the one
diploid individual `0/1`, at a draw of one allele, get 0.5 private alleles
each from the formula, where the two draws, being draws of the same
genotype, can never give an allele to one and not the other. popnei
computes the formula and does not refuse the overlap: `pops` allows an
individual in more than one population, as `docs/specs/stats.md` has it,
and every other statistic here reads such a population without trouble.
The totals, which ask which alleles were called and not which were drawn,
are right whatever the overlap. This was worked out on 24 September 2026
from the formula itself, not measured: the enumeration under "How it is
verified" below runs over one draw per population independently, which is
the same assumption, so it gives 0.5 there too. The pair that does show it
is the nineteenth of that script, which enumerates over the labelled gene
copies of the shared individual and gives 0.

### How it runs

One pass over the blocks of the reader. For each block, rayon over its
rows; for each row, the allele counts of every population from the helper
of `docs/specs/variant.md`, and from them the contribution of that variant
to every statistic asked for. What is kept from one block to the next does
not grow with the variants. Per population it is the variants that counted
and those that reached the draw; the alleles called, the private ones and
the varying ones as three counts and, for the draw, as three sums; the sum
of the observed heterozygosities, the sum of the unbiased expected ones
and how many variants had both; and `num_called_alleles // 2 + 1` sums for
the spectrum. Beside them, two counts for the whole call: the variants
where every population had data and those where every population reached
the draw. The counts of the
populations at one variant are held while that variant is worked on,
because the private alleles need every population's counts at once. They
live in one array of the populations by the alleles of that variant,
which is reused from row to row and never grows with the block.

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
24 September 2026 with `min_num_individuals` 20. The numbers below are
stored in `tests/reference/diversity/panel_num_alleles.tsv`, which
`tests/reference/diversity/make_reference.R` writes, and the tests read
them from there:

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
threshold of 1. `pop1` counts 2, 2, 4 and 1 different alleles at variants
1, 2, 3 and 5, a total of 9 and a mean of 2.25; `pop2` counts 1, 1, 4 and
2, a total of 8 and a mean of 2. In a draw of 4, `pop1` keeps all four of
its variants and gets 2.25, and `pop2` keeps three of them, since at
variant 2 it called only 3 copies in all, and gets 2.3111111111.

At a draw of 2 the same example gives the per variant values that the
draw arithmetic is tested on its own with, run through `vegan::rarefy` on
24 September 2026 on the allele counts of each variant:

| variant | pop1 counts | pop1 | pop2 counts | pop2 |
|---|---|---|---|---|
| 1 | 3, 1 | 1.5 | 5 | 1 |
| 2 | 3, 1 | 1.5 | 3 | 1 |
| 3 | 1, 1, 1, 1 | 2 | 1, 1, 1, 1 | 2 |
| 5 | 4 | 1 | 4, 2 | 1.5333333333 |

whose means are 1.5 for `pop1` and 1.3833333333 for `pop2`. Variant 5 of
`pop2` is the one of the eight a reader cannot do in their head: 4 copies
of allele 0 and 2 of allele 1, and a draw of 2 of the 6 shows both alleles
with chance 8/15, so the alleles expected are 1 + 8/15.

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
Kalinowski (2004), computed by ADZE, the program of Szpiech, Jakobsson and
Rosenberg (2008) that reports the allelic variety of populations at a
common number of called alleles. The mean is over the
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
every population kept are 0, 0 and 0.0008333333. They are stored in
`tests/reference/diversity/panel_private_alleles.tsv`, which
`tests/reference/diversity/make_reference.R` writes. The check is a pytest
test at `calc_pop_diversity`, where the totals and their divisor are read
off the result.

The standardized values are checked by enumerating every draw, which no
program outside popnei does and which tests more than one would. ADZE,
the only program that computes them, evaluates the same closed form from
the same paper, so agreeing with it would say that popnei transcribed
Kalinowski correctly and nothing about whether the formula computes what
the sentence above claims. The enumeration says exactly that: for a small
case it lists every draw of `g` each population can make, weights each by
the ways it can be taken, runs over every combination of one draw per
population, counts the alleles of the population that are in its own draw
and in no other, and averages.

What that check is worth, and what it is not. It never writes the closed
form down, so it catches any error in the algebra: the decomposition into
one term per allele, the chance that an allele shows in a draw, and the
product over the other populations. Changing one factor of the closed form
makes 19 of the 22 pairs below differ, measured on 24 September 2026.
What it does not check is the one
assumption the closed form makes, that the draws of two populations are
independent: the enumeration takes the combinations of one draw per
population as a product measure, which is that same assumption, so on the
populations below the two agree by construction and not by luck. For the
22 pairs below that is sound rather than circular, because their
populations share no individual, and draws from disjoint sets of gene
copies are independent as a fact of the sampling and not as an
assumption. Where populations do share individuals the assumption is false
and the enumeration is as wrong as the formula, which "The cases" above
describes and which the case of the shared individual below demonstrates.

Run in exact rational arithmetic on 24 September 2026 over 18 pairs of a
case and a population, the two agree with a difference of 0, not merely
within a tolerance: the three variants of the worked example that have a
draw of 4; two populations of two alleles at a draw of 2, giving
0.25 and 0.4166666667; two populations of three alleles, 0.75 and
0.4166666667; three populations of two alleles, 0, 0.0833333333 and
0.0833333333; three populations of four alleles, 0.5694444444,
0.6805555556 and 0.4027777778; and a population holding one allele
against one holding two at a draw of 3, 0 and 1. The script is
`tests/reference/diversity/enumerate_private.py`, which writes them to
`tests/reference/diversity/enumerate_private.tsv` with the exact rational
of each beside the two decimals, and the cases become cargo tests at
`calc_pop_diversity` of "The Rust interface" with those numbers as
literals.
`docs/reports/diversity-method/check_by_enumeration.py` is the throwaway
this grew from and is not what the tests read.

Four more pairs are in that file, added on 24 September 2026 because the
18 above are a thin fixture for what they become: seven of their values
are exactly 0 and two exactly 1, and fourteen of their eighteen population
slots have exactly 4 called alleles, so almost nothing there could catch an
implementation that took one population's called alleles for another's.

- **One population of three alleles**, 3, 2 and 1 copies of 6 called at a
  draw of 3, giving 9/4. It is what "The cases" above means by a lone
  population holding as private every allele it called, and it is the only
  pair with one population.
- **Three populations of 3, 5 and 6 called alleles**, `2,1 | 3,1,1 |
  2,2,1,1` at a draw of 2, giving 1/5, 8/25 and 49/75. No population has 4
  called alleles and the fourth allele is in the third population alone.

And a twenty-third pair for the shared individual of "The cases", which is
enumerated differently: over the labelled gene copies of a set of
individuals rather than over allele counts, so that two populations
holding the one diploid individual `0/1` draw the same copies. At a draw
of one allele it gives 0 where the closed form gives 0.5. It is the one
pair of the file whose difference is not 0, and it is what shows that the
overlap "The cases" warns about is real. The file names, for every pair,
which of the two things it enumerated over, and it carries the exact
rational of the closed form and of the enumeration side by side, since on
that pair they differ.

So the file holds 23 pairs, 22 with a difference of exactly 0 and that one
with exactly 1/2. Two guards hold the added cases up, since no program
outside popnei gives their values: the two ways agree exactly, and the
enumeration over labelled gene copies reproduces the enumeration over
allele counts on all 22 pairs whose populations share no individual, which
is what says the 0 of the last pair comes from the sharing and not from a
second way of counting.

Two properties are asserted on the panel beside them: with
`num_called_alleles` equal to the smallest `c` of the dataset the value is
at most the standardized number of alleles of the same population, since a
private allele is an allele; and a population compared against a copy of
itself has 0 private alleles at every draw size, since every allele it
draws the copy can draw too. On the panel at a draw of 20 the values are
0.0112196177, 0.0099715392 and 0.0089014974, and they are literals of a
pytest test. Unlike every other number of this spec they are not stored
beside the reference script, because no program outside popnei computes
them: they come from `docs/reports/diversity-method/panel.py`, which
computes the five quantities in Python as this spec defines them and which
produced the tables here. So they check that popnei's Rust agrees with
that Python and nothing more; what checks the formula itself is the
enumeration above, over the 22 pairs, and the two properties.

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

## The variable variants

### What it gives

How many of the variants vary in a population. A variant varies in a
population when that population called more than one allele there, which
`docs/glossary.md` calls a variable variant and not a polymorphic one: a
polymorphic variant is one whose major allele frequency is below the
polymorphism threshold, and a threshold on a frequency has no meaning for
a draw of alleles. The rarefaction literature, ADZE among it, calls what
this item computes the proportion of polymorphic loci; popnei keeps the
glossary's word, so the statistic is `variable_vars_ratio` and not a
polymorphism ratio. A user gets the **total**, the varying variants, and
the **ratio**, those over the variants that counted for the population.
The total is the `num_variable` that `calc_per_var_distribs` gives in its
`poly_vars_ratio`, over the same variants.

The standardized value is the mean over the variants of

    P = 1 - sum over a of C(n_a, g) / C(c, g)

the chance that a draw of `g` of the `c` called alleles is not all of one
allele, since `C(n_a, g) / C(c, g)` is the chance that every copy drawn is
allele `a`. It is the quantity that makes two populations of different
size comparable: a population of 200 individuals finds a rare allele that
a population of 20 misses, and without the draw it looks the more variable
for that reason alone.

### How it is verified

On a variant of two alleles the expected number of alleles in a draw of
`g` is one plus the chance that the draw varies, because the draw shows
either one allele or two. So on a dataset where every variant has two
alleles `vegan::rarefy` verifies this item as well: the standardized ratio
has to be the standardized number of alleles minus 1, within 1e-12
relative, the two being the same sum of the same terms taken in different
orders. The identity was checked on 24 September 2026 on the counts
(10, 6), (17, 3), (1, 19) and (55, 45) at draws of 2, 4 and 10 and on
(3, 1) at draws of 2 and 4, its 4 called alleles having no draw of 10:
over those 14 pairs the two sides agree to 1.1e-16 or exactly when both
are computed in exact arithmetic, which is the identity itself and says
nothing about any program.

`tests/reference/diversity/make_reference.R` checks those same 14 pairs on
every run, and it checks a stronger thing: `vegan::rarefy`'s own output
against the chance of varying computed from exact binomial coefficients.
What is left over there is `rarefy`'s arithmetic and not the identity,
since `rarefy` computes its terms as `exp(lchoose(c - n, g) - lchoose(c,
g))` through the log gamma function. Over the 14 pairs the largest
difference is 2.44e-15, 1.6e-15 of the value, on the counts (17, 3) at a
draw of 4; five pairs agree exactly. Compared the way `rarefy` computes
it, through `lchoose` on both sides, 13 of the 14 agree exactly and the
fourteenth to 2.22e-16. All three numbers are far inside the 1e-12
relative this item is compared within, and they are given apart because
they measure three different things. Measured on 24 September 2026. It was
checked on the panel too, whose standardized ratios are
0.9283948650, 0.9219209943 and 0.9197370844 against the standardized
allele counts of 1.9283948650, 1.9219209943 and 1.9197370844. A dataset
with variants of more than two alleles has no such identity and no
program here to check it; the worked example covers it.

`vegan::rarefy` measures the standardized number of alleles, and the
standardized ratio is that value minus 1 by the identity above, so it is
not a second measurement:
`tests/reference/diversity/panel_variable_vars.tsv` stores it under a
column name that says so. A test that reads both that file and
`panel_num_alleles.tsv` is therefore checking popnei's two formulas
against one measurement of `vegan`'s and against the identity, which is
what this item claims and all that it claims.

The totals and ratios on the panel, run with `min_num_individuals` 20 on
24 September 2026: 1173, 1177 and 1184 varying variants of 1200, ratios of
0.9775, 0.9808333333 and 0.9866666667, stored in
`tests/reference/diversity/panel_variable_vars.tsv`. The totals are
compared exactly.
They are also the variants `docs/specs/stats.md` counts as variable, so
the two modules have to agree on them: a pytest test at
`calc_pop_diversity` and `calc_per_var_distribs` on the same dataset and
the same `pops` asserts that this total equals the `num_variable` of
`poly_vars_ratio`, at any `poly_threshold`, since that count does not read
the threshold.

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

which is the hypergeometric chance that a draw of `g` of the `c` copies
holds `j` of the rarer ones. `j` runs from `max(0, g - (c - m))` to
`min(m, g)`, the draw holding neither more rarer copies than the variant
has nor fewer than the major allele can leave room for: `m` can be far
above `g`, 48 against 20 on the panel, and a `j` outside that range asks
for a binomial coefficient of a negative argument and a bin below 0. A
variant whose
population called fewer than `g` alleles contributes nothing. So the
entries are not whole numbers, and they sum to the variants that counted,
`num_vars.in_draw` of that population.

### How it is verified

Against `dadi` 2.4.4, whose `Spectrum.from_data_dict(dd, [pop],
projections=[g], polarized=False)` does the same per variant projection,
drops the variants below `g` called alleles, and folds. On the panel with
`min_num_individuals` 20 and `num_called_alleles` 20, run on 24 September
2026, all 33 entries of the three populations agree to their ten printed
digits. They are stored in
`tests/reference/diversity/panel_folded_sfs_dadi.tsv`, which
`tests/reference/diversity/make_reference.py` writes, and the test reads
them from there:

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

Each column sums to 1200, the variants that counted, and the stored
columns sum short of it by 1.3e-11, 7.0e-11 and 4.0e-11, which is `dadi`'s
own error accumulated over 1200 variants and not the eleven roundings of
the decimals, which could give at most 7.5e-14. So the test that checks
that sum compares within a tolerance and not exactly. They are compared
within 1e-12 relative, by a pytest test at `calc_pop_diversity` reading
`folded_sfs` against these literals. That is the tightest margin in this
module: an `f64` projection of popnei's own differs from these stored
values by up to 7.1e-14 relative, `dadi`'s error dominating, which is 14
times inside the tolerance, where the same comparison against
`vegan::rarefy` sits 550 times inside it. Both were measured on 24
September 2026 by recomputing the panel in exact rational arithmetic.
`dadi` masks bin 0 and the bins above `g / 2` in a
folded spectrum and popnei reports bin 0, so the comparison reads
`fs.data` and not the masked array; the difference is one of presentation
and the value is the same.

`dadi` installs and runs on the project's Python, which is 3.14.5 with the
global interpreter lock and which `.python-version` pins by its patch
version. What it does not build on is 3.14.7, the free threading build,
which `uv venv --python 3.14` picks by itself on the owner's machine: its
`nlopt` dependency compiles from source and stops for want of `cmake`,
which is not installed there. Both were measured on 24 September 2026, and
the folded spectrum `dadi` gives is the same on the two interpreters. The
reference script of this module nevertheless makes an environment of its
own with `uv venv --python 3.12`, which takes about 5 seconds, so that the
stored numbers come from the one version of `dadi` this spec names,
whatever the development dependencies of popnei later hold; it says so in
a comment. It is the first reference of
popnei that is not a program run by a shell script or a library of the one
environment, and it is worth that because the projection is the arithmetic
here most easily got wrong, being a distribution over bins and not a
single value. Decided on 24 September 2026; the option not taken was to
drop `dadi` and check the projection inside popnei, as the standardized
private alleles are checked. Enumerating every draw, which is what settled
those, would not settle this one: the projection is already an
enumeration, one term per count of the rarer allele, so listing the draws
would restate it rather than test it.

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
defines them, over the variants that counted for the population and at
which both exist: `H_o` needs one called genotype at least, since it
divides by the called genotypes, and `H_u` needs more than one called
allele, since it divides by `c` times `c - 1`. Both means are over that
one set of variants, so the two divisors are the same number. The `c`
times `c - 1` is the diploid form of that divisor, and what `H_u` needs in
general is as many called alleles as the ploidy, which is what
`docs/specs/stats.md` gives it and what governs here, this item having said
that both heterozygosities are as that spec defines them. At ploidy 2 the
two readings are the same, and every number of this spec is at ploidy 2.
A population
whose mean `H_u` is 0, every variant it counted having one allele, has no
F_IS and gets NaN. It is
Nei's F_IS, the one built from the two mean heterozygosities of a
population read on its own, and not Weir and Cockerham's, which comes out
of a decomposition of the variance across populations and needs more than
one of them.

Two choices are in that line and both were measured on 24 September 2026,
by drawing genotypes at a known F, an individual being heterozygous with
chance `2 p q (1 - F)` and homozygous with `p² + F p q` and `q² + F p q`,
over 2000 variants whose allele frequencies come from the neutral
spectrum, where the variants at frequency x go as 1/x. Twelve settings
were tried, 30 and 100 individuals against the rarest allele at 0.01 and
0.05 against a true F of 0, 0.2 and 0.5, with 50 datasets drawn for each.

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
The standard deviations of the two are within 0.0020 of each other, the
widest gap being 0.0056 against 0.0036 at 100 individuals, the rarest
allele at 0.01 and a true F of 0.5, so the accuracy costs no precision
worth the name. The
ratio of means is what popnei computes.

The second is the unbiased expected heterozygosity against the plain one.
The plain one is too small by about `1 / (c - 1)`, so it makes F_IS too
small by about the same: at a true F of 0, the plain one gives -0.0177
with 30 individuals and -0.0053 with 100, where the unbiased one gives
-0.0007 and -0.0002. popnei uses the unbiased one. `scikit-allel` uses the
plain one, so a user comparing the two libraries sees a difference of
0.0110, 0.0078 and 0.0063 on the three populations of the panel, and this
spec records it rather than following `scikit-allel`.

Both measurements are `docs/reports/diversity-method/fis_sim.py`, run with
the project's Python, which prints the twelve rows.

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
which is popnei's formula with the plain expected heterozygosity. No
function of popnei gives that form, so the check is made in pytest, which
builds it from the two means `calc_per_var_distribs` gives for
`obs_het` and `exp_het` over the same `pops` and the same
`min_num_individuals`: on 24 September 2026 `scikit-allel` gives
-0.0237536998, -0.0258924700 and -0.0247472838 for `p0`, `p1` and `p2`,
stored in `tests/reference/diversity/panel_fis_plain_allel.tsv`, which
`tests/reference/diversity/make_reference.py` writes, and
compared within 1e-12 relative. The unbiased form, which
`calc_pop_diversity` returns and which a second pytest assertion compares
against the literals below, is
-0.0127584868, -0.0181107131 and -0.0184585832 on the same data; it
differs from the plain one by the `c / (c - 1)` factor of each variant's
expected heterozygosity, whose own verification against pyNei and plink2
is in `docs/specs/stats.md`, so no program is run twice for it.

`adegenet` 2.1.11 is the other program with the ingredients, and it is not
the reference: its `Hs` is the plain expected heterozygosity with no
correction for the sample, 0.3481518, 0.3477311 and 0.3427402 on the three
populations of the panel on 24 September 2026. That it carries no
correction was checked by hand on a dataset of two variants and two
populations, where it gives each population the mean of its two variants:
0.5 for a population whose counts are 3 and 3 of 6 at both variants, and
0.4375 for one whose counts are 3 and 1 of 4 at the first, 0.375, and 2
and 2 of 4 at the second, 0.5. The unbiased values of the same four would
be 0.6, 0.6, 0.5 and 0.6666666667. It would
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

The statistics to compute, a bit set as `Needs`, the set of columns a
reader is asked to fill in `docs/specs/block.md`, is:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DiversityStats(u8);

impl DiversityStats {
    pub const NUM_ALLELES: DiversityStats;
    pub const PRIVATE_ALLELES: DiversityStats;
    pub const VARIABLE_VARS_RATIO: DiversityStats;
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
    /// The variants the population called something at and had
    /// `min_num_individuals` called genotypes in. `None` when `pop` is not
    /// a population of the call.
    pub fn num_vars(&self, pop: usize) -> Option<u64>;
    /// Of those, the ones whose called alleles reached
    /// `num_called_alleles`: the variants in the draw for the population.
    pub fn num_vars_in_draw(&self, pop: usize) -> Option<u64>;
    /// The variants that counted for every population, and of those the
    /// ones in the draw for every population. They are the two divisors of
    /// the private alleles.
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
    pub fn num_variable_vars(&self, pop: usize) -> Option<u64>;
    pub fn variable_vars_ratio_in_draw(&self, pop: usize) -> Option<f64>;
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

None. The four this spec had while it was written the owner decided on 24
September 2026, and each is under the item it belongs to with the option
that was not taken: that `num_called_alleles` takes one number and not a
sequence, under "Its Python function"; that these five are a module of
their own and not items of `docs/specs/stats.md`, in the opening; that the
five reference programs become development dependencies of popnei, also in
the opening; and that the standardized private alleles are checked by
enumerating every draw rather than by building ADZE, under "How it is
verified" of that item, where the pairs the two agree on are listed.
A fifth, whether `dadi` is worth an environment of its own, the writer
decided, under "How it is verified" of the folded spectrum, because it
changes no value a user sees and no public API.

## Not in this spec

- F_ST, f_2, Jost's D and the other measures between two populations:
  `docs/specs/dists.md`, which computes them from the same allele counts
  per population in its own pass.
- The unfolded spectrum, which needs an ancestral allele that popnei does
  not read, and the joint spectrum of two populations, which is a matrix
  of one axis per population and which the demographic fits use. Neither
  has a user yet.
- The per variant values of any of these as a frame. A user who wants them
  takes the genotypes with `iter_blocks`, as `docs/specs/stats.md` says
  for its own.
- The expected heterozygosity, the observed heterozygosity, the major
  allele frequency and the polymorphic variants, those below the
  polymorphism threshold: `docs/specs/stats.md`. F_IS above is built from
  the first two, and the item here for the variable variants counts the
  same variants its `num_variable` counts.
- Rarefaction of the heterozygosities. They are frequencies and not
  counts, so they do not grow with the individuals sampled the way the
  number of alleles does, and no program here rarefies them.
- The read ahead thread of section 3 of `docs/architecture.md`: a reader
  over a reader, which `docs/specs/block.md` holds, and which the pass
  here takes like any other reader.
