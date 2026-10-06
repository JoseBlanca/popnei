# The filters module: variants kept by missing data, major allele frequency, observed heterozygosity, linkage and region

September 2026. The `filters` module gives a user of popnei the variants of
a dataset that pass a threshold, before any calculation sees them: those
with few missing genotypes, those whose commonest allele is not too
frequent, those with few heterozygous individuals, those that do not
repeat what a variant near them on the chromosome already said, and those
inside, or outside, the regions of a BED file. It also gives a sample of
them, each variant kept at random with a probability the user gives, and
the first n of them, after which the pass ends and the rest of the source
is not read, and the variants whose FILTER column, in the VCF they came
from, said that they passed. It also
tells the user how many variants each filter was given and how many it
kept. And it keeps, of every variant, the genotypes of the individuals a
user names and drops those of the rest. There is code for every item but
the filter of the variants that passed their FILTER. This spec
develops the row `filters` of the table in section 9
of `docs/architecture.md`, and it covers the three filters that compare one
number of a variant with a threshold, the counts, the filter of
individuals, the filter that takes out the variants that repeat what a
variant before them said, the filter by regions, the filter that keeps
variants at random, the filter that keeps the first n and the filter of
the variants that passed their FILTER. The last was added on 6 October
2026, from issue 9 of the repository, and there is code of its column and
none of the filter; the two
before it were added on 5 October 2026, from issues 6 and 7. The
filter by regions
was added on 26 September 2026; the two before it were
added on 22 September 2026: the
filter of individuals with `docs/specs/stats.md`, whose statistics per
population are the first to need it. It depends on
`docs/specs/block.md`, which has the `Block`, the run of consecutive
variants held as arrays that the variants flow in, and the `BlockReader`
trait of everything that gives blocks; on `docs/specs/variant.md`,
which has the counts of the genotypes and of the alleles of one variant,
the `Variants` that a user puts the filters on, and the `PassStats` in
which the counts reach them; and, for the filter by linkage
disequilibrium alone, on `docs/specs/ld.md`, which has the r² that filter
compares and the dosages it reads the genotypes as. The filter by regions
hands its regions to the readers of `docs/specs/io_vcf.md` and
`docs/specs/io_vars.md`, which skip the variants outside them.

That last dependency puts the `filters` module on the `ld` module and so
on the linear algebra, while `docs/specs/ld.md` is on this one for the
chain of filters of a pass: two modules that use each other, as `variant`
and `block` do. It also splits this spec in two for whoever builds it.
The three threshold filters and the counts are part of the walking
skeleton of section 10 of `docs/architecture.md` and are built. The
filter by linkage disequilibrium comes after `linalg` and `ld`, where
section 9 of that document puts what sits on the linear algebra, and
nothing that is built waits for it.

## The three threshold filters

### What they give

Each filter works out one number for every variant, over all the
individuals of the dataset, and keeps the variant when the number is at
most the threshold the user gave. A variant whose number is exactly the
threshold stays.

A genotype is missing when at least one of its alleles was not called, so
a half called genotype, `0/.` in a VCF, is a missing genotype, as
`docs/glossary.md` has it. A genotype that is not missing is called.

The missing data filter compares the missing rate of the variant,

    missing rate = missing genotypes / individuals

where the individuals are all those of the dataset, and not the ones that
were called at that variant.

The maf filter compares the major allele frequency. "maf" is, in pyNei and
in popnei, the frequency of the major allele, the commonest one, where
most of the literature and plink2 use the same letters for the minor one:

    maf = the largest of the counts of the alleles / called alleles

An allele is counted each time it was called, also in a half called
genotype, whose called allele counts. The called alleles are all the
alleles of the variant that are not missing. Every allele of a
multiallelic variant has its own count, and the variant is not collapsed
to two alleles. When two alleles tie for the largest count the maf is the
same whichever is called the major one. A filter that keeps the variants
with a maf of at most 0.95 takes out the ones that hardly vary among these
individuals.

The observed heterozygosity filter compares

    observed heterozygosity = heterozygous genotypes / called genotypes

where a genotype is heterozygous when it is called and its alleles are not
all the same, at any ploidy. It is used to take out the variants in which
too many individuals are heterozygous, which in most datasets are
paralogous regions read as one site.

A variant with no called allele has no maf, and one with no called
genotype has no observed heterozygosity. A variant with no number is not
kept, whatever the threshold.

### In Python and in TypeScript

```python
Variants.filter_by_missing_data(max_allowed_missing_rate: float) -> None
Variants.filter_by_maf(max_allowed_maf: float) -> None
Variants.filter_by_obs_het(max_allowed_obs_het: float) -> None
```

`Variants` is the handle of `docs/specs/variant.md`: a source of variants
and the steps that were put on it, in order, with no genotypes of its own.
A filter is a step. Each of the three methods adds its step at the end of
the list of the `Variants` it is called on, reads nothing and returns
nothing, as `list.sort()` does, so that a user who writes `v2 =
v1.filter_by_maf(0.95)` gets a `None` and an error at the next line, and
not two names for one filtered object.

The steps are run by whatever consumes the variants, a calculation,
`write_vars` or `iter_blocks`. Each of them makes its passes over the
source, a pass being one reading of it from its start, and every pass
takes the steps that the `Variants` has when it starts, in their order. A
step can be added at any time, also after a pass: a user looks at the
distributions of the variants as they are in the file, adds a filter to
the same `Variants` and calculates again. A step that is added while a
pass is running, from the loop of an `iter_blocks`, holds from the next
pass and changes nothing in the one that runs. A user who wants two sets
of thresholds over one file opens it twice, which reads the header and
nothing else.

They carry the names and the arguments of `filter_by_missing_data`,
`filter_by_maf` and `filter_by_obs_het` of `pynei/var_filters.py`, which
are functions that take a `Variants` and return another. The owner decided
on 21 September 2026 that in popnei a step changes the `Variants` and is a
method of it. What a user holds is a recipe that is run later, and
building it one object at a time leaves objects behind that nothing should
be done with, and a call whose result is not assigned filters nothing and
says nothing. The options not taken were pyNei's functions, and one
function that takes the three thresholds at once and returns a new
`Variants`, which has no place for a filter whose result depends on the
steps before it, as the one by linkage disequilibrium does. That a step
can be added after a pass was decided the same day; the option not taken
was that the first pass locks the steps.

Three more differences from pyNei. `max_allowed_missing_rate` has no
default, so the three thresholds are written by the user, and a call
without one is a `TypeError`. In pyNei it is 0.0, which keeps only the
variants with every genotype called, 26 of the 500 of
`tests/reference/vcf/many.vcf`. The owner decided it on 21 September 2026:
whoever calls a filter means to filter by some rate, and no rate is the
natural one. A threshold that is not between 0 and 1, both included, or
that is NaN, is a `ValueError` at the call, which names the argument and
the value. pyNei takes any number: with a negative one no variant passes,
which its `test_filter_missing` asserts for -0.1, and with one above 1
every variant that has a number does, so a 95 written for 0.95 filters
nothing and says nothing. The owner's rule of 21 September 2026 is that a
wrong input of a function is a `ValueError`. And a filter of a kind that
the `Variants` already has is a `ValueError`, which names the kind and the
threshold that is set. pyNei takes it, and adds the counts of the two
together. The owner decided it on 21 September 2026: two threshold filters
of one kind do what the stricter of them does alone, so a second one says
that the user has lost track of what the `Variants` holds, and it is what
running the cell of a notebook twice gives; the option not taken was
replacing the threshold of the first.

A user, or a function of theirs that receives a `Variants`, reads what it
holds in `variants.steps`, a tuple with a `Step` for each step, in order:

```python
@dataclass(frozen=True)
class Step:
    kind: str                  # "missing_data", "maf", "obs_het", "individuals", "ld",
                               # "regions", "excluded_regions", "random" or "first_n"
    args: dict[str, object]    # {"max_allowed_maf": 0.95}
```

The `repr` of a `Variants` shows its source and its steps. pyNei has
neither. The owner decided it on 21 September 2026: a second filter of one
kind is refused, so a user has to be able to see which are set, and a
notebook whose cells were run out of order is where they most need it. The
option not taken was to add nothing. `args` is a dict so that the steps of
the later items, which take other arguments, fit in it.

In TypeScript, `variants.filterByMissingData(maxAllowedMissingRate)`,
`variants.filterByMaf(maxAllowedMaf)` and
`variants.filterByObsHet(maxAllowedObsHet)`, which return nothing, and
`variants.steps`, an array of `{kind, args}`. A
threshold out of range or not given, and a second filter of one kind, are
an `Error` at the call.

### Half called genotypes, variants with nothing called, and what pyNei asserts

The missing rate and the maf of one variant do not count the same data.
The variant `0/. ./. ./. ./. ./.` of five individuals has a missing rate of
1, because its one half called genotype is missing, and a maf of 1, from
its one called allele. The maf filter keeps it at a threshold of 1 and the
observed heterozygosity filter never does, because it has no called
genotype.

Neither the maf filter nor the observed heterozygosity filter asks for a
minimum of called data, where the statistics of pyNei ask for 20 called
genotypes by default: `_filter_chunk_by_maf` passes `min_num_samples=0`,
and `_calc_obs_het_per_var` has no such argument. A variant with one
called genotype, heterozygous, has an observed heterozygosity of 1.
Inherited from pyNei. A user who does not want the variants that have
little called data puts the missing data filter before these two.

The number is one division of the two counts as `f64`, compared with `<=`,
as in pyNei, and not a product of the threshold and the denominator. The
two differ: 29 missing genotypes of 100 individuals pass a threshold of
0.29 as 29/100 <= 0.29, and would not as 29 <= 0.29 * 100, which is
28.999999999999996. In `tests/reference/vcf/many.vcf` 114 of the 500
variants have a missing rate of exactly 0.04.

`test_filter_missing`, `test_filter_mafs` and `test_filter_obs_het` of
`test/test_filters.py` assert which of three variants of five individuals
stay at two or three thresholds each, with whole missing genotypes and no
half called one. `test_filtered_vars_can_be_iterated_more_than_once` and
`test_chained_filters_can_be_iterated_more_than_once` assert that a second
pass over a filtered `Variants` gives what the first gave, for chunks of 1
to 4 variants. No test of pyNei has a half called genotype in a filter.

### What pyNei does that is odd

`_calc_maf_per_var` of `pynei/gt_counts.py` takes the alleles it counts
from the largest allele of the whole chunk, pyNei's block. In a chunk where
nothing is called there are no alleles, the frequencies are a frame with
no columns, and its largest value per row is NaN, so the variants are
dropped, as a variant with nothing called is in any other chunk. Run on
two variants of five individuals with every genotype `./.`, at a threshold
of 1: no variant kept. The result does not depend on the chunk, and popnei,
which counts row by row, gives the same.

### How it runs

A filter is a reader over another reader, as section 1 of the architecture
says. It takes a block from its source, works out which rows stay with
rayon over the rows of the genotypes, compacts the block in place with
`retain_vars` of `docs/specs/block.md`, and gives it on. Nothing is kept
from one block to the next but the two counts of the next item. It takes
the blocks of its source at whatever size they come, so it needs no
`reblock` before it, the reader of `docs/specs/block.md` that puts blocks
back to one size, and the blocks it gives are of uneven size.

A filter always needs the genotypes. When its consumer says which fields it
wants, with `set_needs`, the filter asks its source for those fields and
for the genotypes, and the blocks it gives hold the genotypes also when
the consumer did not ask for them.

It keeps the rules of a reader of `docs/specs/block.md`: a block left with
no variant is not given and the next one is taken; after an error, of its
source or its own, it gives `None` at every call and does not call its
source again; and a source that gives a block of no variants is the error
that `reblock` gives for it. Before it counts, it runs `check` on the
block, which says whether the arrays of a block have the sizes the block
states, because it cuts the genotypes into rows by the sizes the block
states, and a block whose genotypes were not read, with variants in it, is
the error of `docs/specs/variant.md` for a field that is not there.

Several filters on one `Variants` are several readers, one over the other,
in the order of the steps, so each sees only what the one
before it kept. Whether one reader that applies several thresholds in one
pass over the row is faster is left to a measurement.

### How it is verified

Against bcftools 1.24, on `tests/reference/vcf/many.vcf` of
`docs/specs/io_vcf.md`, read with every variant given, those that failed
their FILTER too: 500 variants of 50 diploid individuals, one in ten with
three alleles, and 257 half called genotypes among its 25000. bcftools
counts a half called genotype as pyNei does, missing as a genotype and its
called allele in the allele counts: for the variant of `cases.vcf`, of the
same directory, with `./.`, `0|1` and `.|0`, its `F_MISSING`, the fraction
of missing genotypes, is 0.666667 and its `AN`, the called alleles, 3. The
commands, with `$t` the threshold:

    bcftools view -H -i "F_MISSING<=$t" many.vcf
    bcftools view -H -Q $t:major many.vcf
    bcftools view -H -i "N_PASS(GT=\"het\") <= $t * (N_SAMPLES - N_MISSING) && N_MISSING < N_SAMPLES" many.vcf

`-Q $t:major` keeps the variants whose major allele frequency is at most
`$t`. In the third, `N_PASS(GT="het")` is the number of individuals with a
heterozygous genotype, `N_SAMPLES` the individuals and `N_MISSING` those
with a missing genotype. It multiplies where popnei divides, which the
variants of `many.vcf` do not tell apart at these thresholds. The variants that each
command keeps were compared with the ones pyNei keeps at commit ef0ca6e, by
their positions, which are all different in that file. They are the same
at every threshold:

| filter | threshold | kept of 500 | the first five kept, by position |
|---|---|---|---|
| missing data | 0 | 26 | 1259, 2110, 2480, 3072, 3257 |
| missing data | 0.04 | 215 | |
| missing data | 0.1 | 455 | |
| maf | 0.5 | 35 | 1074, 1296, 1481, 1962, 2110 |
| maf | 0.8 | 384 | |
| maf | 0.95 | 480 | |
| observed heterozygosity | 0.1 | 22 | 1185, 3516, 3923, 4515, 5921 |
| observed heterozygosity | 0.25 | 79 | |
| observed heterozygosity | 0.5 | 369 | |

There are variants exactly on the threshold in seven of the nine rows: 26,
114 and 50 with a missing rate of 0, 0.04 and 0.1, 7 and 2 with a maf of
0.5 and 0.8, and 1 and 13 with an observed heterozygosity of 0.25 and 0.5,
counted with pyNei's three functions named below. plink2 v2.0.0-a.7.7
keeps the same variants as the missing data filter at the three
thresholds, with `--geno $t`, which takes out the variants with a missing
rate above `$t`, and `--vcf-half-call m`, which reads a half called
genotype as missing. It was not run for the other two filters, because it
has no way of reading a half called genotype, by its `--help`, in which
the genotype is missing and its called allele is counted.

A script, `tests/reference/filters/make_reference.py`, runs the three
commands at those thresholds and stores the positions each one keeps in a
file beside it. The cargo tests, made at `next_block` of a `FilteredReader`
over a `VcfReader` on `many.vcf`, assert the number kept in each row of the
table and the five positions where the table has them, as literals. They
are run with blocks of 7 variants and of the default size, and the
variants kept are the same.

Against pyNei, a pytest test made at the three Python functions: `many.vcf`
is read by both libraries, popnei with `only_passed=False` since pyNei
gives every variant and opened again for each filter and threshold of the
table, and the genotypes and the positions of popnei's blocks, joined,
are compared exactly with those of pyNei's chunks, joined. A second pass
over the same filtered `Variants` gives the same. A threshold of -0.1, of
1.5 and of NaN is a `ValueError` in each of the three, and a call with no
threshold a `TypeError`.

The worked example, which becomes the first cargo test, made at
`VarFilter::filter_block` on a block built by hand: six variants of five
diploid individuals. The numbers are those of pyNei's
`_calc_gt_is_missing`, `_calc_maf_per_var` and `_calc_obs_het_per_var` on
these genotypes, and bcftools keeps the same variants at the thresholds of
0.4, 0.88 and 0.25 below.

| variant | genotypes | missing rate | maf | observed heterozygosity |
|---|---|---|---|---|
| 1 | 0/0 0/1 0/0 0/0 0/. | 0.2 | 8/9 = 0.8889 | 1/4 = 0.25 |
| 2 | 0/0 0/1 0/0 ./. 0/. | 0.4 | 6/7 = 0.8571 | 1/3 = 0.3333 |
| 3 | 0/1 2/3 0/1 2/3 ./. | 0.2 | 2/8 = 0.25 | 4/4 = 1 |
| 4 | ./. ./. ./. ./. ./. | 1 | none | none |
| 5 | 0/0 0/0 0/0 0/0 1/1 | 0 | 8/10 = 0.8 | 0 |
| 6 | 0/. ./. ./. ./. ./. | 1 | 1/1 = 1 | none |

The missing data filter keeps variant 5 at a threshold of 0, variants 1, 3
and 5 at 0.2, variants 1, 2, 3 and 5 at 0.4, and the six at 1. The maf
filter keeps variant 3 at 0.25, variants 3 and 5 at 0.8, variants 2, 3 and
5 at 0.88, and variants 1, 2, 3, 5 and 6 at 1. The observed heterozygosity
filter keeps variant 5 at 0, variants 1 and 5 at 0.25, and variants 1, 2, 3
and 5 at 1.

The TypeScript test, under node, reads `many.vcf` from a `Uint8Array` with
every variant given and asserts the number kept in the first row of each
filter of the table, its five positions, and the `Error` of a threshold of
1.5.

## The counts of each filter

### What they give

For each filter of a pass, how many variants it was given and how many it
kept. The user reads in them how many variants a result was calculated on,
and which filter took out the rest.

### In Python and in TypeScript

```python
@dataclass(frozen=True)
class FilteringStats:
    vars_processed: int
    vars_kept: int
```

The counts are in what a pass produces, and not in the `Variants`: every
result of a consumer has a `pass_stats`, the `PassStats` of
`docs/specs/variant.md`, and so has the iterator that `iter_blocks`
returns. Its `filtering` is a dict of the kind of each filter,
`"missing_data"`, `"maf"`, `"obs_het"`, `"ld"`, `"regions"`,
`"excluded_regions"`, `"random"` or `"first_n"`, to its
`FilteringStats`, in the order of the steps, and it is empty when the `Variants` had no filter.
`PassStats` also says whether the filter of the first n ended the pass
before its source ended, in `stopped_early`, which "The counts, and a pass
that the filter ended" of that filter has.

```python
variants.filter_by_missing_data(0.04)
variants.filter_by_maf(0.8)
distribs = calc_per_var_distribs(variants)
distribs.pass_stats.filtering
# {"missing_data": FilteringStats(500, 215), "maf": FilteringStats(215, 163)}
```

`FilteringStats` mirrors the one of `pynei/var_filters.py`. pyNei keeps the
counts in the `Variants`, and gives those of its last pass with
`gather_filtering_stats(variants)`, which popnei does not have. The owner
decided it on 21 September 2026: a `Variants` of popnei can get a step
after a pass, and counts kept in it would then be those of steps it no
longer has, while the counts in a result are of the pass that gave that
result, whatever is done to the `Variants` afterwards. The option not
taken was pyNei's function. The order of the dict is another difference:
pyNei gives the last filter first.

In TypeScript `filtering` is an object of kind to `{varsProcessed,
varsKept}`, two numbers.

### A pass that was not finished

The `pass_stats` of the iterator of `iter_blocks` can be read while the
blocks are taken, and it then has the counts of the variants that were
read up to there. They can be of more variants than the blocks the user
got hold, because `iter_blocks` ends in a `reblock`, the reader of
`docs/specs/block.md` that puts blocks back to one size, which keeps
variants for its next block. A calculation that fails gives no result, and
so no counts.

### How it runs

The filter of one pass is an object of that pass, which owns its two
counts as plain numbers. Every pass builds its own filters from the steps
of the `Variants`, so no count is shared between two passes and none needs
a lock. The owner decided it on 21 September 2026; the option not taken
was one filter object that the Python `Variants` owns and every pass uses,
behind a lock.

For the counts to be read when a pass ends, whoever starts the pass keeps
the chain of readers: a calculation of the core borrows the reader of each
of its passes, `&mut dyn BlockReader`, and does not take it. `BlockReader` gets a method,
`filtering_stats`, that gives the counts of every filter of the chain: a
source gives none, and a reader over another reader gives those of its
source, with its own before them when it is a filter. The method has no
default, so that a reader over a reader that forgets to pass on the counts
of its source does not compile.

A binding crate builds the chain of a pass with `chain_of` of "The Rust
interface" and keeps it while the consumer runs. When the consumer
returns, the binding crate reads the counts from
the chain and hands them to the Python or the TypeScript package, which
puts them in the result in the order of the steps, the reverse of the one
the chain gives. The `num_vars` of the same `PassStats` is not a count of
a filter: the binding crate adds up the variants of the blocks that the
consumer took from the chain. While an `iter_blocks` runs it is the
variants of the blocks the user got, which can be fewer than the last
filter has kept, for the `reblock` above.

`write_vars` is the consumer that counts none of them itself. Its loop
over the blocks is the core's, so no block of that pass reaches the
binding crate, and the count of the variants that were written comes back
from the core with the sink, as "Its Python and TypeScript functions" of
the writer in `docs/specs/io_vars.md` says. The binding crate hands that
number on and reads the counts of the filters from the chain it lent, as
every other consumer does.

### How it is verified

There is no reference program for a count of variants given and kept, and
none is needed: the counts of a chain of filters are the numbers of
variants that the chained commands of bcftools keep. On `many.vcf`, with
the missing data filter at 0.04, the maf filter after it at 0.8 and the
observed heterozygosity filter after that at 0.5, bcftools keeps 215, 163
and 106 variants, the first at the positions 1111, 1407 and 1518, and pyNei
gives

    missing_data: 500 processed, 215 kept
    maf:          215 processed, 163 kept
    obs_het:      163 processed, 106 kept

The cargo test, made at `filtering_stats` of the outermost `FilteredReader`
of that chain after its last block, asserts the three pairs, the
observed heterozygosity one first, and the 106 variants with those three
positions. On the worked example above, with thresholds of 0.4, 0.88 and
0.25 in the same order, the pairs are 6 and 4, 4 and 3, and 3 and 1, and
variant 5 is the one kept.

Against pyNei, the pytest test puts those three filters on `many.vcf` in
both libraries and compares the `pass_stats.filtering` of the iterator of
a whole `iter_blocks` with pyNei's `gather_filtering_stats` after a pass,
kind by kind, and the order of popnei's dict with the order of the steps.
A second `iter_blocks` gives the same counts and not the double, and one
over a `Variants` with no filter an empty dict. The tests of the steps,
made at the three methods: each returns `None`; `steps` is empty for a
`Variants` just opened and has, after the three filters, their three kinds
in order with the thresholds under the names of the arguments; a second maf filter is a
`ValueError`, also with another filter between the two; a filter added
after a whole `iter_blocks` holds in the next one, whose counts have it;
and a filter added inside the loop of an `iter_blocks` takes no variant
out of that pass. The TypeScript test asserts the three pairs of the
chain and the `Error` of a second maf filter.

## The filter of individuals

### What it gives

Keeps, of every variant, the genotypes of the individuals the user names
and drops those of the rest. Every variant stays. It is what a user puts
on a `Variants` before a calculation over some of their individuals, one
population of the dataset or the individuals that have a phenotype, and
the threshold filters after it in the steps see the kept individuals
alone: the missing rate then divides by them, and the maf and the
observed heterozygosity are over their genotypes, as in pyNei, where a
chunk that went through `filter_samples` holds the kept individuals and
no other. The kept individuals come in the order the user named them,
so the filter is also the way to put the individuals in the order a user
wants, the populations together, and the rows of the per individual
statistics and of a distance matrix come out that way; it costs nothing,
since the gather of a row takes any order. The owner decided it on 22
September 2026; the option not taken was the order of the source, which
is pyNei's: `filter_samples(v, ["ind05", "ind00", "ind49"])` on
`many.vcf` gives `ind00, ind05, ind49`, whatever the order of the
argument, where bcftools's `-s` gives the order of the argument.

### In Python and in TypeScript

```python
Variants.filter_individuals(individuals: Sequence[str]) -> None
```

It is a step, as the three threshold filters are: it adds itself at the
end of the steps of the `Variants`, reads nothing and returns nothing, and
every pass that starts afterwards runs it. Its kind is `"individuals"`
and its `args` is `{"individuals": (...)}`, the names as a tuple. After
it, `variants.individuals` and `variants.num_individuals` are the kept
individuals, since they are what the next pass gives, and the `pops` of
`docs/specs/stats.md`, the populations a statistic is calculated for as
a dict of name to individuals, names individuals among them.

It mirrors `filter_samples` of `pynei/var_filters.py`, a function that
takes a `Variants` and returns another, with the name `docs/glossary.md`
gives, individual for pyNei's sample. The differences:

- It takes names alone. pyNei takes a sequence of names or a slice; its
  type hint says indices too, but `numpy.isin` matches them as names,
  and `filter_samples(v, [0, 1, 2])` on `many.vcf` keeps no individual.
  The objectives give a user the individuals as a tuple of names, so
  names are what the user has.
- A name that is not an individual of the source is a `ValueError` that
  names it. pyNei's `numpy.isin` drops it in silence: `filter_samples(v,
  ["ind05", "nope"])` on `many.vcf` gives a `Variants` of one individual,
  run at commit ef0ca6e.
- A name that is there twice is a `ValueError` that names it. pyNei keeps
  the individual once.
- No name at all is a `ValueError`: every source of popnei refuses a
  dataset of no individuals, as `docs/specs/block.md` says.
- A second filter of individuals on one `Variants` is a `ValueError`, as
  a second threshold filter of one kind is, by the owner's rule of 21
  September 2026: two lists keep the individuals that are in both, which
  is one list. pyNei takes it.
- It has no entry in `pass_stats.filtering`, since it takes no variant
  out. pyNei's `gather_filtering_stats` lists it under the kind `sample`,
  with every variant processed and kept.
- The kept individuals come in the order of the argument, and not of the
  source, under "What it gives".

The three refusals of the names are made at the call, against the
individuals of the source, with `resolve_individuals` of "The Rust
interface", which the reader of the filter calls too when a pass builds
its chain. pyNei's `calc_ld_and_dist_per_pop` puts several filters of
individuals over one filtered `Variants`, one per population; in popnei
a calculation over several populations takes `pops`, as those of
`docs/specs/stats.md` do, and a user who wants two sets of individuals
over one source opens it twice.

In TypeScript, `variants.filterIndividuals(names)`, which returns
nothing, with the step `{kind: "individuals", args: {individuals:
[...]}}` in `steps`, and the errors an `Error` at the call.

### What pyNei asserts

`test_filter_samples` of `test/test_filters.py` asserts, on three variants
of five individuals, that keeping the first three by name and by a slice
gives their genotypes;
`test_sample_filter_metadata_does_not_depend_on_the_call_order`, that
the `samples` and the `num_samples` of the filtered `Variants` are the
kept ones, read before a pass or after it, and
`test_filtered_samples_are_a_tuple`, that they are a tuple; and
`test_several_vars_can_share_one_filtered_source`, that two filters of
individuals over one filtered `Variants` each give their own
individuals, with chunks of 1 to 4 variants.

### How it runs

A reader over a reader, with the rules of a threshold filter under "How
it runs" above: it takes a block from its source, compacts the genotypes
of every row with `retain_individuals` of "The Rust interface", and gives
the block on with the `num_individuals` of the kept ones and its other
columns as they were. The compaction is in two passes over the array of
the block, because the kept individuals can come in any order, i5 before
i1, so a genotype cannot be moved over one that is still to be read, and
because the rows shrink, so where a row will start is inside the row
before it. The first pass runs with rayon over the rows at their source
width, which are disjoint, and gathers the kept genotypes of each row
through a buffer of the thread, kept individuals x ploidy alleles, back
into the front of that row; the second packs the shortened rows to the
front of the array, one after another, on one thread, as `retain_vars` of
`docs/specs/block.md` does. It allocates no block, keeps nothing from
one block to the next, and needs no `reblock`. Its `individuals()` gives
the kept names, and its ploidy,
its chromosome table and its filtering stats are those of its source. It
always needs the genotypes, and passes on what its consumer asks for with
them added. It sits in the chain where its step is among the steps, so a
threshold filter before it counts over every individual and one after it
over the kept ones. The blocks are the size of its source's, worked out
from the individuals of the source and not from the kept ones. The
`reblock` that `iter_blocks` puts at the end of the chain sizes its blocks
for the individuals of the reader it is given, which are the kept ones, as
`docs/specs/block.md` says of `Reblock::new`, so above 500 individuals,
where the default size falls below the largest it takes, the blocks a user
reads after the filter hold more variants than the source's.

### How it is verified

Against bcftools 1.24 on `many.vcf` of the threshold filters, read with
every variant given, which is `only_passed` false in the tests of the
three layers: the default of `open_vcf` keeps the variants that passed
their filters alone, and 50 of the 500 variants of `many.vcf` did not
pass, so the 500, the 423 and the 26 below are the counts of that reading
and of no other. `bcftools view -s ind05,ind00,ind49 many.vcf` keeps
the three individuals, in that order, which `bcftools query -l` prints of
its output, and the missing data filter at 0 after it,

    bcftools view -s ind05,ind00,ind49 many.vcf | bcftools view -H -i "F_MISSING<=0"

keeps 423 of the 500 variants, the first five at the positions 1000,
1037, 1074, 1111 and 1148, where the same filter over the 50 individuals
keeps 26. In one command, `-s` with `-i`, bcftools applies the filter
before it takes the individuals out, and keeps 26. pyNei's
`filter_samples` and `filter_by_missing_data(0)` after it keep the same
423, with the counts 500 given and 500 kept for the first and 500 and
423 for the second. The genotypes of the three at the first variant,
position 1000, are `1|1`, `1/1` and `1/1`, and at position 1074 `0/1`,
`2|1` and `1|2`.

The cargo tests, at `next_block` of an `IndividualsReader` over a
`VcfReader` on `many.vcf`: the blocks hold 3 individuals, `individuals()`
gives the three names, the genotypes of each are the column of the
source at every variant, 500 variants come out and `filtering_stats` is
empty; with a `FilteredReader` of the missing data filter at 0 over it,
423 variants with those five positions first and the counts 500 and 423;
with the same filter under it instead, 26. At `resolve_individuals`, a
name that is not an individual, a name twice and no name are each the
error, and the three names give the indices 5, 0 and 49. At
`Block::retain_individuals`, on the six variants of the worked example of
the threshold filters: keeping i5 and i1 gives rows of two genotypes,
`0/. 0/0` for variant 1 and `./. 0/1` for variant 3, and an index of 5 is
the error with the block as it was.

Against pyNei, a pytest test at `filter_individuals`: `many.vcf` in both
libraries with the three individuals, and the genotypes of popnei's
blocks, joined, are those of pyNei's chunks, joined, column by name,
since pyNei gives them in the order of the source; with the missing data
filter at 0 after it, both keep 423 variants, and popnei's
`pass_stats.filtering` has the missing data filter alone, with 500 and
423. The tests of the step: `steps` has `("individuals", {"individuals":
(...)})`, `individuals` and `num_individuals` are the kept ones, and a
second filter of individuals, an unknown name, a name twice and no name
are each a `ValueError`. The TypeScript test, under node, asserts the
three names, the 423 and the `Error` of an unknown name.

## The filter by linkage disequilibrium

### What it gives

The variants of the dataset with the ones that say what a variant before
them already said taken out, so that what is left carries each piece of
information once. A principal component analysis or a kinship over
variants that repeat one another counts that stretch of the genome as
many times as it has variants, and the filter is what a user puts before
them.

Two variants say the same thing when their r² is high. r² is the square
of the correlation, across the individuals, between the dosages of the
two variants, where the dosage of a genotype is how many of its alleles
are not the major allele of its variant; it is 1 when the dosage of an
individual at one variant fixes its dosage at the other and 0 when
knowing one says nothing about the other, and `docs/specs/ld.md` defines
it, with the individuals that count for a pair, those called at both, and
the pairs that have no r².

The filter walks the variants in the order they come. The window of a
variant is the variants the filter has already kept that are on that
variant's chromosome and no more than `max_dist` base pairs behind it. A
variant is kept when

- its called genotypes hold two dosages at least, and
- its r² against every variant of its window is at most
  `max_allowed_r2`.

A pair whose r² is not defined does not drop the candidate: only an r²
above the threshold does. So the first variant of each chromosome whose
called genotypes hold two dosages is always kept, and a variant whose
called genotypes all hold one dosage is always dropped, having nothing to
tell any other variant apart with.

### In Python and in TypeScript

```python
Variants.filter_by_ld(max_allowed_r2: float, max_dist: int) -> None
```

It is a step of the `Variants`, like the three threshold filters above,
with the same rules: it adds itself at the end of the list, returns
nothing, is run by whatever consumes the variants, and a second filter of
its kind on one `Variants` is a `ValueError`. Its `Step` has the kind
`"ld"` and both arguments in its `args`, and its counts reach the user
under that kind in the `filtering` of a `PassStats`.

Neither argument has a default, as `max_allowed_missing_rate` has none. A
`max_allowed_r2` that is not a number from 0 to 1, and a `max_dist` below
1, are a `ValueError` at the call that names the argument and the value.
The binding crate takes `max_dist` as a signed integer and checks it
itself, so that a negative one is that `ValueError` and not the
`OverflowError` that pyo3 raises when a negative number is asked of a
`u64`.

It carries what `filter_by_ld_and_maf` of `pynei/var_filters.py` does,
and differs from it in five ways.

- **The major allele frequency is not in it.** pyNei's function filters by
  the major allele frequency first and by linkage disequilibrium after,
  in one call with one pair of counts. In popnei a filter is a step, so a
  user writes `variants.filter_by_maf(0.95)` and then
  `variants.filter_by_ld(...)`, gets the counts of the two apart, and
  chooses the order. Decided here; it follows from the owner's decision
  of 21 September 2026 that a filter is a step.
- **A variant is compared with every kept variant of a window and not
  with the last kept one alone**, and the window is a distance along a
  chromosome, where pyNei reads neither the chromosome nor the position
  (under "What pyNei does that is odd"). The owner decided on 22
  September 2026 that popnei compares within a window and that the
  chromosome ends it; the option not taken was pyNei's rule. plink2's
  `--indep-pairwise` is the program that compares within a window, and
  the owner's decision was made on it, but the set popnei keeps is not
  plink2's and cannot be: "How it is verified" has the measurement and
  what the tests compare instead.
- **The threshold is on r² and not on the absolute value of r.** pyNei's
  argument is called `min_allowed_r2` and is compared with the absolute
  value of the correlation. The owner decided on 22 September 2026 that
  every value of linkage disequilibrium in popnei is r². A pyNei
  threshold of 0.1 is a popnei threshold of 0.01.
- **The name of the threshold says which way it works.**
  `max_allowed_r2` is the largest r² a kept variant may have against a
  kept variant of its window, which is what `max_allowed_maf` and
  `max_allowed_missing_rate` are for their filters, where pyNei's
  `min_allowed_r2` raises the number of variants kept as it rises.
  Decided here.
- **A missing genotype takes its individual out of the pair it is in**,
  where pyNei gives it a dosage of -1. `docs/specs/ld.md` has the
  measurement.
- **A pair whose r² is not defined does not drop the candidate.** pyNei
  keeps the variants whose r against the reference passes
  `numpy.abs(r) < min_allowed_r2`, and a NaN passes no comparison, so a
  candidate whose r cannot be computed is dropped. Two variants that both
  have variance get no r² when the individuals called at both hold one
  dosage, which is a pair that says nothing about either variant, not a
  pair that says they are the same. Decided here.

In TypeScript it is `variants.filterByLd(maxAllowedR2, maxDist)`, which
returns nothing, with the two refusals above as an `Error` at the call.

### Which variant of a linked pair is kept, and the cases

Of two variants whose r² is above the threshold, the one that comes first
is the one kept. A filter of popnei takes a block, compacts it in place
and gives it on, as section 1 of `docs/architecture.md` has it, so a
variant it has given away cannot be taken back, and only the variants
still inside the window can be reconsidered. plink2 instead removes
variants from windows it has already passed, which is why its set is
another one (under "How it is verified").

The filter reads the dosages over every individual of the dataset. A user
who wants them read over one population puts the filter of individuals
before it.

A variant with no called genotype has no dosage at all and is dropped, as
a variant of one dosage is.

A block with variants and no position is the error of
`docs/specs/variant.md` for a field that is not there, as a block with no
genotypes is for the other filters: this filter asks its source for the
chromosome and the position besides the genotypes, whatever its consumer
asked for, and the blocks it gives hold all three.

The window of a variant is the variants kept *behind* it, so this filter
needs the variants of each chromosome to come together and in order of
position. Two variants at one position are allowed and are 0 apart, so
each is in the other's window. A variant whose position is below the one
before it on the same chromosome, and a variant on a chromosome that had
already ended, are an error that names the variant, its position and the
one before it. The other part of popnei that needs the same order is the
cutting of the variants into the resampling groups of the standard errors
of the distances between populations, `docs/specs/dists.md`, which
refuses the same two cases where the groups are stretches of a
chromosome. Everything else reads a source in any order, and
`docs/specs/io_vars.md` has a test that writes a block of four variants
that are not sorted, so these two refuse what the rest of popnei takes.
Decided here: the alternative is to subtract two positions that can run
backwards, which on a `u64` wraps to a distance of 18 million million
million and puts the pair outside every window without a word.

### What pyNei does that is odd

Read and run in pyNei at commit ef0ca6e.

`_filter_chunk_by_ld` compares each candidate with the last kept variant
and with nothing else, so a variant is kept although it repeats the
variant two before it, which the one in between hid.

It reads neither the chromosome nor the position: it walks the variants
in the order of the file and carries the last kept one from chunk to
chunk, so the last variant kept on one chromosome is what the first
variant of the next chromosome is compared with, and a variant is
compared with one a whole chromosome arm away as readily as with its
neighbour.

The first variant of the first chunk is kept whatever it is, `if ref_gt
is None: selected_vars.append(0)`, and it becomes the reference. When
every one of its genotypes holds the same value, every r against it is
NaN, no NaN is below the threshold, and nothing else is ever kept: the
whole dataset comes out as that one variant. The value counted is the one
`to_012` writes, so a missing genotype is a -1 that differs from every
dosage and gives the reference variance: the collapse needs a variant
with no missing genotype at all. Run on 22 September 2026 on six variants
of five individuals with `min_allowed_r2=0.9` and `max_allowed_maf=1`,
pyNei keeps 1 of the 6 when every genotype of the first variant is `0/0`,
6 of the 6 when one of those five is `./.` instead, and 5 of the 6 when
the same six variants are given with a variant of two dosages first.
popnei drops a variant of one dosage wherever it is, so it has neither
the collapse nor the rescue by a missing genotype.

The name `min_allowed_r2` is the r below which a variant counts as
unlinked and is kept, so raising it keeps more variants: on
`tests/reference/dists/panel.vcf.gz`, with its `max_allowed_maf` at its
default of 0.95, pyNei keeps 768 of the 1200 variants at 0.1 and 1167 at
0.3.

### How it runs

A reader over a reader, as the three threshold filters are. It takes a
block from its source, works out which rows stay, compacts the block in
place with `retain_vars` and gives it on, and it takes the blocks of its
source at whatever size they come.

What it keeps from one block to the next is the window: for each variant
it has kept whose position is within `max_dist` of the newest variant it
has read, and which is on that variant's chromosome, its genotypes, with
its chromosome and its position. It keeps the genotypes and not the
dosages because the whole window is one operand of the products when a
set of candidates arrives, so the three matrices are built for all of its
variants together and then let go. Between two blocks the window costs
one byte for each allele, which is two bytes for each individual and each
variant of it at a ploidy of 2, where the three matrices are 24: for 250
kept variants of 1000 individuals, 0.5 MB held and 6 MB while a set is
settled. A variant leaves the window when the reader passes `max_dist`
beyond it or reaches another chromosome, which is worked out each time a
set of variants is settled and not at every variant.

How many variants a window holds is set by the dataset and not bounded by
the filter. They are unlinked to each other by construction, which keeps
a window short where the variants of a dataset are linked to their
neighbours; a dataset whose variants carry no linkage leaves every one of
them with variance in the window, so a window as wide as a chromosome
holds every variant of it. A window whose genotypes or whose dosages the
machine has not the memory for is refused rather than taken, with the
case of `docs/specs/ld.md` for a matrix this machine cannot hold.

The variants of a block can be compared with the variants that were
already in the window when the block arrived in one set of the products
of `docs/specs/ld.md`, so the work that the window bounds is done as
matrix products and not one pair at a time. What cannot be done that way
is a candidate against the variants kept inside the same block: whether
one candidate is kept decides what the next one is compared with, so
those are settled in order. How much of a block is taken at a time is for
the implementer to choose.

Whether the result changes with the size of the blocks is the test that
"How it is verified" names: it must not, since the rule reads positions
and never block boundaries.

### How it is verified

The r² is verified against plink2 v2.0.0-a.7.7 in `docs/specs/ld.md`,
which agrees with the formula there to 5.6e-16. The rule on top of it is
verified against the r² that plink2 itself gives, with three properties
of the set it keeps:

- every kept variant has two dosages at least among its called genotypes;
- no two kept variants on one chromosome and within `max_dist` of each
  other have an r² above `max_allowed_r2`; and
- every dropped variant whose called genotypes hold two dosages has an r²
  above `max_allowed_r2` against some variant that was kept before it and
  is within `max_dist` on its chromosome.

All three are checked in the reference script against the float64 matrix
that `plink2 --r2-unphased square bin` writes, so every decision of the
filter is pinned to plink2's numbers and not to popnei's own. Each of
the three catches a different way of being wrong: the first is what a set
that also kept the 68 variants of one dosage would fail, since the other
two say nothing about those variants; the second is what a set that kept
two linked variants would fail; and the third is what a set that dropped
a variant it could have kept would fail, since such a variant has two
dosages and nothing above the threshold kept before it in its window,
which is the negation of what the third property asks.

Over the set the reference script itself built, all three are 0 because
the rule that built it is the rule they state, so the script also works
each one out over a set that breaks it and stops rather than write a
result of 0 for any of the three: every variant of the dataset kept,
which leaves 68 of one dosage; each variant compared with the last kept
one alone, which is what `_filter_chunk_by_ld` of pyNei does and which
leaves 707 pairs above the threshold inside a window; and the kept set
with its first ten variants taken out, which leaves 64 variants dropped
for no reason. The counts of the table below are a check beside the three
and not the only one. On
`tests/reference/ld/ld.vcf.gz` of `docs/specs/ld.md`, 500 variants of 100
individuals on two chromosomes 250000 bp long with 68 variants of one
dosage, run on 22 September 2026, all three properties hold at every
setting below and the first kept variant of chr2 is always `chr2:1000`:

| max_dist | max_allowed_r2 | kept of 500 | the first five kept, by position |
|---|---|---|---|
| 10000 | 0.1 | 84 | chr1:1000, chr1:10000, chr1:16000, chr1:22000, chr1:27000 |
| 10000 | 0.3 | 133 | chr1:1000, chr1:5000, chr1:7000, chr1:11000, chr1:15000 |
| 50000 | 0.3 | 85 | chr1:1000, chr1:5000, chr1:7000, chr1:11000, chr1:15000 |
| 250000 | 0.3 | 85 | the same five |

The cargo tests, made at `next_block` of an `LdFilteredReader` over a
`VcfReader` on that file, assert the number kept in each row and the five
positions where the table has them, with blocks of 7, of 64 and of the
default size, which have to keep the same variants.

plink2's own `--indep-pairwise` keeps 65, 84, 41 and 32 variants in those
four rows, none of them a variant of one dosage. Its window is in
kilobases where `max_dist` is in base pairs, so the four commands are
`--indep-pairwise 10kb 0.1`, `10kb 0.3`, `50kb 0.3` and `250kb 0.3`; with
the base pairs written into them instead, `10000kb 0.1` keeps 18 and
`50000kb 0.3` keeps 32, both of which put each whole chromosome in one
window.

Its sets too have no linked pair inside the window, so both rules give
sets of unlinked variants and plink2's are the smaller: at 50000 and 0.3
its 41 share 12 variants with popnei's 85, and its first kept variant of
chr1 is chr1:11000 where popnei's is chr1:1000. It removes variants from
windows it has already passed, which a filter that gives its blocks on
cannot do, and reproducing its set would need its stepping of the window
and its order of removal, which its `--help` does not state. So
`--indep-pairwise` is not what popnei's set is compared with; the three
properties above are.

What the two sets cost each other was measured on the same file on 22
September 2026, with the r² of plink2 throughout. Neither rule leaves a
pair above the threshold, and the mean r² still standing between the
variants each one kept, over the pairs of them inside the window, is:

| max_dist | max_allowed_r2 | popnei kept | its mean r² left | plink2 kept | its mean r² left |
|---|---|---|---|---|---|
| 10000 | 0.1 | 84 | 0.04317 | 65 | 0.00657 |
| 10000 | 0.3 | 133 | 0.14662 | 84 | 0.08450 |
| 50000 | 0.3 | 85 | 0.10558 | 41 | 0.07920 |
| 250000 | 0.3 | 85 | 0.04850 | 32 | 0.03432 |

So the two rules trade the same two things against each other and neither
is the better one at every threshold: popnei keeps more variants, and
those variants carry more of the linkage disequilibrium the user asked to
be rid of, all of it under the threshold they set. The knob that moves
popnei to plink2's sparsity is that threshold, and the two are not
interchangeable: at a window of 50000 bp, where plink2 at 0.3 keeps 41
variants with a mean r² left of 0.0792, popnei reaches 46 and 0.0679 at
0.15 and 35 and 0.0494 at 0.1. A user who pruned with plink2 at a
threshold and writes the same number here gets a denser set, and this is
where they are told so.

Neither rule bounds what several kept variants together say about
another one: both compare pairs, which is what `--indep-pairwise` is
named for, and the denser set has more pairs to do it over, 1764 inside
the window against plink2's 241 in the last row.

The worked example, the first cargo test, made at `LdFilter::filter_block`
on the five variants of 6 individuals of "How it is verified" of
`docs/specs/ld.md`, whose r² are worked out there by hand and confirmed
by plink2, with v4 the variant of one dosage:

| max_dist | max_allowed_r2 | kept | why |
|---|---|---|---|
| 5000 | 0.5 | v1, v5 | v2 and v3 are above 0.5 against v1, v4 has one dosage |
| 5000 | 0.7 | v1, v2, v5 | v2 is 0.675 against v1, v3 is 0.754 |
| 1000 | 0.5 | v1, v3, v5 | v3 and v5 have no kept variant within 1000 bp |

The counts of that first row are 5 variants given and 2 kept, and a
missing data filter at 1 before it and an observed heterozygosity filter
at 1 after it give 5 and 5, 5 and 2, 2 and 2.

Against the Python API, a pytest test made at `Variants.filter_by_ld`:
`ld.vcf.gz` filtered at the four settings of the table gives the four
counts and the five positions, and the blocks a whole `iter_blocks`
yields hold those variants and no others; the `Step` it adds has the kind
`"ld"` and both arguments under the names of the arguments, and comes
after a maf filter added before it; the `filtering` of the `pass_stats`
has an `"ld"` with 500 given and 84 kept in the first row of the table,
after the `"maf"` of a maf filter put before it, which is what a user
writes in place of pyNei's one call; a `max_allowed_r2` of -0.1, of 1.5
and of NaN, a `max_dist` of 0 and of -1, and a call with either argument
missing are a `ValueError` and a `TypeError` as "In Python and in
TypeScript" says; a second `filter_by_ld` is a `ValueError`, also with
another filter between the two; and a source whose positions go backwards
within a chromosome is the `ValueError` of "Which variant of a linked
pair is kept, and the cases", on a VCF written for it.

There is no test against pyNei for which variants are kept: pyNei's rule
compares a candidate with the last kept variant alone and reads no
position, so the two keep different sets on any dataset where a variant
is linked to one that is not its predecessor. What the two share is the
r², which `docs/specs/ld.md` compares between them.

The TypeScript test, under node, reads `ld.vcf.gz` from a `Uint8Array`
and asserts the 133 variants of the second row of the table, their first
five positions, and the `Error` of a `maxAllowedR2` of 1.5.


## The filter by regions

### What it gives

It keeps the variants that fall inside the regions a user lists in a BED
file, or, with `exclude`, the ones that fall outside all of them. A user
keeps the variants of the genes they study, or the ones a capture
targeted, and takes out those of repeats or of regions known to map
badly.

A BED file is a text with one region on each line: the chromosome, the
start and the end, separated by tabs, and any columns after those three,
which are not read. BED counts the bases from 0 and leaves the end out, so
the line `chr1 0 2000` is the first 2000 bases of chr1. A region of popnei,
as `docs/glossary.md` has it, counts from 1 as a VCF does and includes both
ends, so that line is the region from 1 to 2000, and a line of start `s`
and end `e` is the region from `s + 1` to `e`.

A variant is inside a region when its chromosome has the name of the
region's, written the same way, and its position, the POS of its VCF, is
in the region. Only the position is looked at, and not the bases its
reference allele covers, so a deletion that starts before a region and
reaches into it is outside. The owner decided it on 26 September 2026:
it is what `bcftools view -T` and `plink2 --extract bed0` do, and the
same for every variant of one base. The option not taken was any base the
reference allele covers, as `bcftools view -R` does, which would have the
filter and the skip read the alleles, and the vars file keep the longest
reference allele of each batch. A variant on a chromosome
that the BED does not name is outside. Regions that overlap or that touch
act as the one region they cover together, and the order of the lines
does not matter.

### In Python and in TypeScript

```python
Variants.filter_by_regions(bed_path: str | Path, exclude: bool = False) -> None
```

It is a step, as the other filters of this spec are: it adds itself to the
`Variants`, reads no variant and returns nothing. It reads the BED file at
the call, so an error of the file comes at the call and the regions are
those the file had then. The step is of the kind `"regions"`, or
`"excluded_regions"` with `exclude`, and its `args` are the path and the
number of regions left once the overlapping ones are joined:
`{"bed_path": "genes.bed", "num_regions": 412}`. A second step of one kind
is refused, as for every filter of this spec, and the two kinds can stand
together: a user keeps the variants of the genes and excludes from them
those of the repeats. The owner decided on 26 September 2026 that the
filter keeps by default and excludes with `exclude=True`.

In TypeScript, `variants.filterByRegions(bed, {exclude = false})`, with
`bed` a `Uint8Array` with the bytes of the BED file, and `args` of
`{numRegions}` alone, since there is no path. A `bed` that is not a
`Uint8Array` and an `exclude` that is not a boolean are an `Error` at the
call.

pyNei has no filter by regions, so nothing is mirrored and nothing
differs.

The lines of a BED file that are not regions are skipped: an empty line,
one that starts with `#`, and one whose first word is `track` or
`browser`, which the tools of the UCSC genome browser write and which
bedtools skips, by its documentation; bcftools 1.24 refuses a BED that
has them, tried on 26 September 2026. The word is the whole of what comes
before the first blank, a space or a tab, or the end of the line, so
`tracks1 0 5`, whose chromosome is `tracks1`, is a region, as bcftools
1.24 and plink2 v2.0.0-a.7.7 read it. A BED that starts with the two
bytes of gzip, `1f 8b`, is read through gzip, as a VCF is, and the three
bytes of the UTF-8 byte order mark, `ef bb bf`, which some editors of
Windows write at the start of a text, are dropped from the start of the
text, so that they are not read as the start of the first chromosome
name. What is refused, a `ValueError` in Python that names the file
and the line, and in TypeScript an `Error` that names the line: What is refused, a `ValueError` in Python that names the file
and the line, and in TypeScript an `Error` that names the line:

- A line of fewer than three columns separated by tabs. A line with
  spaces between its columns is one of these, and the message says that
  BED separates them by tabs.
- A chromosome that is empty, a line that starts with a tab, which names
  no chromosome a variant can be on; plink2 v2.0.0-a.7.7 refuses it too.
- A start or an end that is not a whole number of 0 or more, or that is
  above 18446744073709551615, the largest number of 64 bits.
- A start that is not below its end. BED allows a start equal to its end
  for a point between two bases, which holds no position of popnei's.
- A BED with no region.

The lines are counted from 1 over the whole file, the skipped ones among
them, so the number in a message is the one an editor shows. A line that
ends in a carriage return before its newline, as a file written on
Windows does, is read without it, as bcftools 1.24 and plink2
v2.0.0-a.7.7 read one, tried on 26 September 2026. A start and an end are
written in digits alone and fit in 64 bits: `+99`, which bcftools reads
as 99, is refused. The name of a chromosome is compared with the name the
source gives byte for byte, so a name that is not UTF-8 text is not
refused and matches no variant.

### The cases a reader of the rules would not guess

- The filter needs the chromosome and the position. It asks its source
  for them besides what its consumer asked for, and a source that lacks
  them, a vars file written without those columns, gives at its first
  block the error of
  `docs/specs/variant.md` for a field that a consumer depends on and did
  not get, a `ValueError` in Python that names the field.
- A BED that names its chromosomes `1` where the source names them `chr1`
  keeps nothing, and the consumer then gives the error of a pass that gave
  no variant, whose message has the counts of this filter, 500 given and 0
  kept on `many.vcf`. popnei does not match `chr1` with `1`: a match that
  is wrong for one file would keep the wrong variants and say nothing.
- The counts are those of the variants this filter was given and kept,
  whether or not its source skipped the variants outside the regions,
  which "How it runs" says.

### How it runs

The regions of each chromosome are held sorted and joined, and a variant
is looked up by its position among them with a binary search, over the
rows of a block with rayon as the threshold filters are. What is kept from
one block to the next is the regions and the two counts.

The owner decided on 26 September 2026 that the source skips the variants
this filter would take out, so that they are not built at all, where the
option not taken was a filter that looks at every variant and leaves the
skipping for later. When this filter is the first filter of variants of
the chain, so that nothing stands between it and the source but, at most,
the filter of individuals, it hands its regions to the source when the chain is built,
with `skip_outside` of `docs/specs/block.md`, and the source gives only
the variants the filter can keep:

- The VCF reader reads CHROM and POS of each line in the serial pass that
  already finds its FILTER, and a line outside gets no row: its columns
  of individuals, 92 in 100 of the time of a read on one thread by the
  profile of "Speed" of `docs/specs/io_vcf.md`, are never parsed. The line
  is still read, and decompressed when the file is bgzipped.
- The vars file reader does not read a batch whose chromosomes, with the
  smallest and the largest position of each that the footer keeps
  (`docs/specs/io_vars.md`), overlap no region, or, with `exclude`, lie
  each inside one region. It seeks past the batch, and nothing of it is
  decompressed. In the other batches it gives every variant, and the
  filter takes them out.

The source counts the variants it passed over since the pass started, and
gives that number with `num_skipped`. The filter's count of the variants
it was given is the variants that reached it plus that number, read when
its counts are read, so the counts are the ones the filter gives without
the skip. A filter of variants gives 0 as its `num_skipped`, since it
passed its source nothing to skip, so a second filter by regions, which
has the first below it, adds nothing that the first counted. A line that the VCF reader skips for
its FILTER is counted by neither, as it is not a variant of the source.

When a filter of variants comes before this one, the source is not handed
the regions: skipping would take variants from that filter too and change
its counts. This filter then looks at every variant it gets, which is the
same variants and slower, and it is given what the filter before it kept. A user who puts the filter by regions first
gets the skip. When a `"regions"` and an `"excluded_regions"` step are
both in the chain, the first of the two can hand its regions to the
source and the second looks at every variant.

A line that the VCF reader skips is not parsed past CHROM, POS and FILTER,
so a wrong genotype, a genotype of another ploidy or an undeclared allele
in it is no error with the skip and is one without it, and when a filter
of variants comes first. This is decided here, and it is what the reader
does already for a column that no consumer asked for: "How it runs" of
the reader in `docs/specs/io_vcf.md` does not check a column it does not
parse. A line outside the regions holds no variant of the pass either
way.

What the reader checks of every line whatever is asked for, the shape of
the line of "How it runs" of `docs/specs/io_vcf.md`, it checks of a line
it skips too: the nine first columns are there and are UTF-8, the FORMAT
has a `GT` key, and, when the genotypes or the text of the lines are asked
for, the line has one column after the FORMAT for each individual of the
header. A line that fails any of them gets a row, so the parse gives its
error with the skip as without it. A plain VCF that was cut inside its
last line has nothing else that says it was cut, and that line has too
few columns. The checks count the tabs of the line, which the skip reads
whole. A POS that does not parse as a number gives the line a row, so the
parse gives the error of that column whatever the regions are. The serial
pass reads a POS written in digits alone that fit in 64 bits, and any
other gets a row: `+5`, which the parse reads as 5, reaches the filter,
which takes it out when it is outside, so it is counted among the
variants the filter was given and not among those the source skipped.

The offer holds from the next block the source builds. A block the reader
one block ahead of `docs/specs/block.md` had built before the offer
reached its thread is given whole, and the filter takes out what it does
not keep, so the counts are the same and fewer variants are skipped.

With the skip, a chromosome whose lines are all skipped gets no number in
the table of chromosome names, so the numbers of the chromosomes can
differ with and without it; the names and the positions of the variants
do not.

### How it is verified

Against bcftools 1.24 and plink2 v2.0.0-a.7.7, which keep a variant by
its position alone, as this filter does: `bcftools view -T file.bed`,
whose `-T` reads a file named `.bed` as a BED and matches the position,
and `plink2 --extract bed0 file.bed`, whose `bed0` says that the file
counts from 0. `bcftools view -R`, which reads a region with an index,
keeps a variant whose reference allele reaches into a region, and is not
the reference here. On `many.vcf` of `docs/specs/io_vcf.md`, read with
every variant given because neither program honours FILTER, with this
BED:

    track name=test
    # a comment
    chr2    19000   30000
    chr1    0       2000
    chr1    5000    5100
    chr1    4990    5050
    chr2    10249   10250
    chr3    0       100000

it has six regions, which join into five: chr1 1 to 2000, chr1 4991 to
5100, chr2 10250 alone, chr2 19001 to 30000 and chr3 1 to 100000. Given
the six lines of regions alone, since bcftools refuses the first two,
`bcftools view -H -T` keeps 45 variants and `-T ^` 455, and `plink2
--vcf-half-call m --extract bed0` and `--exclude bed0` keep the same 45 and
455, run on 26 September 2026. The 45 are 28 of chr1 up to 2000, the three
at chr1 4996, 5033 and 5070, chr2 10250, and the 13 of chr2 from 19001,
counted from the positions of the file. chr3 has no variant. The cargo
tests, made at `next_block` of a `RegionsReader` over a `VcfReader` on
`many.vcf`, in blocks of 7 variants and of the default size, assert the 45
positions and the 455, the counts 500 given and 45 kept, and the same with
the skip of the source and without it, the variants compared by the names
of their chromosomes and their positions.

The same variants written as a vars file in batches of 100 variants have
five batches, whose regions in the footer are chr1 1000 to 4663; chr1 4700
to 8363; chr1 8400 to 10213 and chr2 10250 to 12063; chr2 12100 to 15763;
and chr2 15800 to 19463. With the BED above the reader skips the fourth
batch alone, and with `exclude` none, since no batch lies inside one
region. A cargo test at the vars file reader asserts which batches it
decompressed, with a count that the reader keeps for the tests, and that
the variants and the counts are those of the VCF.

The worked example, the first cargo test, made at `RegionFilter` of "The
Rust interface" on the six lines of `write.vcf` of the VCF writer in
`docs/specs/io_vcf.md`, read with every variant given, with the BED

    chr1    99      100
    chr1    250     251
    chr1    999     1000
    chr2    0       1

keeps chr1 100, chr1 1000 and chr2 1, and with `exclude` chr1 250, chr1
1001 and chr2 1500, as `bcftools view -T` and `-T ^` do and as `plink2
--extract bed0` does, run on 26 September 2026. `chr1 99 100` is position
100 alone, and `chr1 999 1000` holds 1000 and not 1001. `chr1 250 251` is
position 251, so the deletion `AT` at chr1 250, which covers 250 and 251,
is outside; `bcftools view -R` keeps it.

The pytest tests, made at `filter_by_regions`, assert the 45 of `many.vcf`
and their counts, the `steps` of both kinds with their `args`, the
refusal of a second step of one kind, each of the four refusals of a BED
with the line it names, and the error of a pass over a source with no
positions. The TypeScript test asserts the 45 and the 455.

## The filter that keeps variants at random

### What it gives

It keeps each variant with a probability that the user gives, the keep
rate, a number from 0 to 1: at 0.1 it keeps about one variant in ten. A
user runs an analysis on a sample of a large dataset, in a fraction of the
time, with the variants of the sample spread over the whole file.

Whether a variant is kept is drawn from a generator of random numbers that
starts again from the user's seed at the start of every pass. The
generator gives one number from 0 to 1 for each variant that reaches the
filter, in the order in which the variants reach it, and the variant is
kept when its number is below the keep rate. So the variants kept are a
function of the seed and of the place of each variant among the variants
the filter is given. Every pass of a `Variants` goes through the same steps
over the same source, gives the filter the same variants in the same order
and draws the same numbers for them, so every pass keeps the same variants.
That is what lets a calculation that reads the source twice, the PCA of
variants and the GWAS with the GRAMMAR-Gamma approximation, see one sample
in both passes, and what lets two calculations on one `Variants` be
compared. The GRAMMAR-Gamma approximation of `docs/specs/gwas.md` reads the
source a second time to estimate the factor that corrects its test
statistics, and the PCA reads it a second time for the weights of each
variant. The numbers are drawn on one thread, one after another, so the
size of the blocks and the threads of the calculation do not change them.

The owner decided on 5 October 2026 that the draw is made this way, from a
generator, one number for each variant in its order. The option not taken
was a hash of the seed with what identifies the variant, its chromosome and
position or its number in the file, which would give a variant the same
draw whatever came before it. What this way gives instead is that the
sample depends on which variants reach the filter, as "The cases" below
says.

The generator is SplitMix64, of Steele, Lea and Flood (2014), the one that
`java.util.SplittableRandom` of Java uses. Its state is one integer of 64
bits, which starts as the seed. Each draw adds a constant to the state and
mixes the sum into the number it gives, all of it wrapping at 2^64:

    state = state + 0x9E3779B97F4A7C15
    z = state
    z = (z xor (z >> 30)) * 0xBF58476D1CE4E5B9
    z = (z xor (z >> 27)) * 0x94D049BB133111EB
    draw = z xor (z >> 31)

The number from 0 to 1 is the top 53 bits of the draw, `draw >> 11`,
divided by 2^53,
which is exact in an `f64` and always below 1, so at a keep rate of 1 every
variant is kept and at 0 none is. It is a few lines of integer arithmetic,
so it needs no new dependency and gives the same numbers natively and in
wasm.

### In Python and in TypeScript

```python
Variants.filter_randomly(keep_rate: float, seed: int = 42) -> None
```

It is a step, as the other filters of this spec are: it adds itself to the
`Variants`, reads no variant and returns nothing. Its kind is `"random"`,
under which its counts reach the user, and its `args` are `{"keep_rate":
0.1, "seed": 42}`. The owner decided on 5 October 2026 that the seed has a
default and that it is 42: a call without a seed gives the same sample
every time, so a result can be reproduced from the code alone, as every
other result of popnei can, and a user who wants a second sample passes
another seed. The option not taken was a seed the user must always give.
The docstring says that the same seed, over the same source and steps,
keeps the same variants, and the default is a `pub const` of the core.

A keep rate that is NaN, below 0 or above 1 is a `ValueError` at the call,
which names the argument and the value, and one that is no number a
`TypeError`, as for the three threshold filters. A seed that is not a whole
number, a `True` among them, is a `TypeError`, and a whole number below 0
or above 2^64 - 1 a `ValueError`, as `max_dist` of the filter by linkage
disequilibrium is refused. A second filter of this kind is refused, as for
every filter of this spec.

In TypeScript, `variants.filterRandomly(keepRate, {seed = 42})`, with
`seed` a whole number from 0 to 2^53 - 1, the largest whole number a
JavaScript number holds exactly; any other value is an `Error` at the call.
That is the one difference from Python a user sees: a seed above 2^53 - 1,
which Python takes, cannot be given in TypeScript. The same seed gives the same variants in both.

pyNei has no such filter.

### The cases

- The draws are of the variants that reach the filter, so a filter before
  it changes which variant gets which number. The MAF filter and then this
  one keeps another sample than this one and then the MAF filter, and a
  VCF read with `only_passed=False` another sample than the same VCF with
  `only_passed=True`. A VCF and the vars file written from it with no step
  give the same sample, since they hold the same variants in the same
  order. The filter of individuals before it changes no draw, since it
  takes out no variant, and neither does the skip of the filter by regions
  before it, which gives the same variants as the filter does without it.
- At a keep rate of 0 no variant is kept, and the consumer gives the error
  of a pass that gave no variant, with these counts.
- It reads no field of a variant. It asks its source for what its consumer
  asked for and nothing more, so a `write_vars` that asks for the positions
  alone gets them through it. It is a filter of variants, so a filter by
  regions after it is not handed to the source, as "How it runs" of that
  filter says: a user who wants the skip puts the filter by regions first.

### How it runs

The filter keeps the generator and its two counts from one block to the
next. For each block it draws one number for each row, in order, on the
thread that asked for the block, and compacts the block in place with
`retain_vars` of `docs/specs/block.md`. It does not use rayon: one draw is a
few integer operations and the draws have to be made in order. Every pass
builds its own filter from the step, with the generator at the seed, as
every pass builds its own threshold filters. It keeps the rules of a reader
of "How it runs" of the threshold filters: a block left with no variant is
not given, and after an error it gives `None` at every call.

### How it is verified

The reference program is Java's `java.util.SplittableRandom`, whose
`nextLong` is the draw above and whose `nextDouble` is the number from 0 to
1 of the same rule, run with OpenJDK 26.0.2.1 on 5 October 2026 by
`java tests/reference/filters/SplitMix.java`, with the `java` of Homebrew's
`openjdk`, `/opt/homebrew/opt/openjdk/bin/java`, which is not on the PATH of
the owner's machine; `/usr/bin/java` of macOS is a stub that asks for one.
From a seed of 1234567 it
gives the first five draws 6457827717110365317, 3203168211198807973,
9817491932198370423, 4593380528125082431 and 16408922859458223821, the
values that the reference code of SplitMix64, `splitmix64.c` of Sebastiano
Vigna, gives too, as the task "Pseudo-random numbers/Splitmix64" of
Rosetta Code lists them. A cargo test asserts them.

The worked example, which becomes the first cargo test, made at
`RandomFilter::filter_block` on a block of ten variants built by hand, at a
keep rate of 0.5 and a seed of 42. The ten numbers from 0 to 1 that Java
gives from 42 are, to six decimals, 0.741565, 0.159910, 0.278601, 0.344191, 0.038030, 0.868228,
0.218405, 0.800632, 0.339931 and 0.618482, so the filter keeps variants 2,
3, 4, 5, 7 and 9, and its counts are 10 given and 6 kept. The same ten
variants given as a block of 3 and a block of 7 keep the same six. The
tests assert the ten numbers to the bit, as `SplitMix.java` prints the bits
of each, 0x3fe7bae644c5fd6d, 0x3fc477f199d93378, 0x3fd1d499d5c4c3e6,
0x3fd607387fc392b8, 0x3fa378b0b4489040, 0x3febc8863f47901b,
0x3fcbf4b38e229bb4, 0x3fe99ec6bdd3d3c5, 0x3fd5c16e1dc2cf5e and
0x3fe3ca9ae7052fee, and that a filter whose keep rate is exactly the first
of them, 0x3fe7bae644c5fd6d, drops the first variant: a variant is kept
when its number is below the keep rate, not at it.

On `tests/reference/vcf/many.vcf`, read with every variant given, those
that failed their FILTER among them, 500 variants:

| keep rate | seed | kept of 500 | the first five kept, by position on chr1 |
|---|---|---|---|
| 0.1 | 42 | 45 | 1148, 1666, 1777, 1888, 2332 |
| 0.5 | 42 | 243 | 1037, 1074, 1111, 1148, 1222 |
| 0.1 | 7 | 49 | 1037, 1962, 2147, 2332, 2591 |

The table comes from `tests/reference/filters/random_draws.py`, a Python
version of the rule written for this spec, which checks its draws against
the five of Java from 1234567 and applies them to the variants of
`many.vcf`, run on 5 October 2026. The cargo tests,
made at `next_block` of the reader of this filter over a `VcfReader` on
`many.vcf`, in blocks of 7 variants and of the default size, assert each
row's count and five positions, and the counts 500 given and 45 kept for
the first row.

Every pytest and TypeScript test of both filters opens `many.vcf` with
every variant given, `only_passed=False` and `onlyPassed: false`, as the
numbers of these tables are of the 500. With the default, which leaves out
the 25 variants whose FILTER is `q10`, the numbers are others: 42 kept at
0.1 and a seed of 42.

The tests of the passes, made at the Python functions, over `many.vcf` with
this filter at 0.1 and a seed of 42:

- Two `iter_blocks` give the 45 variants of the table, and so do a
  `calc_pairwise_kosman_dists` and a `do_pca_from_variants` on the same `Variants`,
  whose `pass_stats` have 45 variants.
- `do_pca_from_variants` with 10 components, which makes the second pass
  for the weights, does not fail with the error of a second pass that gave
  other variants, and its result is that of `do_pca_from_variants` over the
  45 variants written to a vars file with no step. The same for `calc_gwas`
  with a kinship and `use_grammar_gamma_approx=True`, which makes its
  second pass with no check of its own: its result is that of the same call
  over that vars file, so both of its passes saw the 45.
- A keep rate of -0.1, 1.5 and NaN is a `ValueError`, a seed of -1 a
  `ValueError` and of 1.5 a `TypeError`, and a second filter of this kind
  a `ValueError`. `steps` has the kind and both arguments.
- The filter of individuals before this filter keeps the same 45.
- At a keep rate of 0, `calc_pairwise_kosman_dists` gives the error of a
  pass that gave no variant, with this filter given 500 and keeping 0.
- The filter by regions before this one, with the BED `chr1 3000 9000` and
  `chr2 0 5000`, and this filter at 0.3 with a seed of 42, keeps 50
  variants, the first eight at 3072, 3109, 3183, 3257, 3405, 3590, 3627 and
  3701 of chr1, over `many.vcf` and over the vars file written from it,
  whose reader skips what is outside the regions. The 50 are the draws of
  Java from 42 over the 162 variants that `bcftools view -T` keeps in those
  regions, worked out on 5 October 2026.

The cargo tests also run the first row of the table over a vars file of
`many.vcf` in batches of 7, in pools of rayon of one thread and of four and
through the reader one block ahead, and get the 45 each time; and build the
third row, a seed of 7, through `chain_of`, which gives 49.

The TypeScript test asserts the 45 and their first five positions, the same
variants from a second pass, and the `Error` of a keep rate of 1.5.

## The filter that keeps the first n variants

### What it gives

It keeps the first n variants that reach it, those that the steps before it
kept, and then ends the pass: the rest of the source is not read, and the
calculation goes on to its result with the variants it got. A user tries an
analysis on a large file quickly, before running it on the whole, where
without this filter every pass reads the file to its end. When fewer than n
variants reach it, it keeps them all and the pass reads the whole source.

The first n are those of the start of the file. On a VCF sorted by position
they are the start of the first chromosome and not a sample of the genome.
The filter that keeps variants at random, above, is the one for a sample
of the whole file, and the two can be put on together, at random first and
then the first n, to get n variants spread over the part of the file that
was read.

### In Python and in TypeScript

```python
Variants.filter_first_n(num_vars: int) -> None
```

It is a step: it adds itself to the `Variants`, reads no variant and
returns nothing. Its kind is `"first_n"`, and its `args` are `{"num_vars":
1000}`. `num_vars` is a whole number of 1 or more: 0 or a negative number
is a `ValueError`, and a float, a `True` or anything that is no number a
`TypeError`, as for `max_dist` of the filter by linkage disequilibrium.
Issue 7 asked for a `ValueError` for a float; the `TypeError` is what every
argument of popnei that takes a whole number gives for one. A second filter
of this kind is refused.

No step that takes variants out can be added after it: the three threshold
filters, the filter by linkage disequilibrium, the filter that keeps
variants at random and both kinds of the filter by regions, added after a
filter of the first n, are a `ValueError` at the call, which names the kind
of the step and says that the filter of the first n is set. Then n is
always the number of variants the calculation gets, and a user who wants n
variants that pass the MAF filter puts the MAF filter first. The filter of
individuals takes out no variant and is accepted after it. The owner
decided it on 5 October 2026; the option not taken was to allow any step
after it and to say in the docstring that n then counts the variants
before the later filters.

A step that breaks both rules, a second MAF filter after a MAF filter and
the filter of the first n, is refused as a second filter of its kind, since
moving it before the filter of the first n would not make it acceptable.
The core checks both in that order, in one function that `chain_of` and
both binding crates call. A `num_vars` of 0 is refused by the core too, at
the call, with a message that names `num_vars`, `numVars` in TypeScript, as
the other wrong values of it do. A second filter of the first n names the n
that is set and the one that was asked for.

In TypeScript, `variants.filterFirstN(numVars)`, with `numVars` a whole
number from 1 to 2^53 - 1.

pyNei has no such filter. `desired_num_chunks` of `Variants.from_vars`
stops reading after that many chunks, which counts chunks and not
variants.

### The counts, and a pass that the filter ended

The filter is given whole blocks and keeps from the last one only the
variants it needs, so its counts are of the blocks it took: the variants
of those blocks as given, and n as kept. With blocks of 7 variants and n of
10 it is given 14 and keeps 10. The filters before it count the same
blocks, the rows of the last block that this filter did not keep among them,
so the counts of every filter of a pass that this filter ended depend on
the size of the blocks, and are of the part of the source that was read
and not of the whole of it. A user who reads that the MAF filter kept 900
of 1000 variants would take it for the whole file, which may hold
millions.

So the counts of a pass say whether this filter ended it. `PassStats` of
`docs/specs/variant.md` gets a field, `stopped_early`, true when the filter
of the first n kept its n variants and ended the pass, and false
otherwise, also when it was given fewer than n. It is true also when the
source held nothing after the n-th variant, since the filter does not read
on to find out. The owner decided on 5 October 2026 that the counts say
it; the option not taken was to leave it to the counts of this filter.
It is in the `repr` of `PassStats`, where a user who prints the counts sees
it. When a pass ends there, the `num_vars` of its counts is n.
Every pass of a `Variants` ends at the same variant, so a calculation that
reads the source twice reads the same n variants in both passes.

### How it runs

The filter takes a block from its source. When the variants it has kept
and the block together are at most n, it gives the block whole; otherwise
it keeps the first rows of the block up to n, with `retain_vars`, and
gives it. Once it has given n variants, whether the n-th ended a block or
fell inside one, every later call gives `None` without asking its source. It reads no field of a variant and asks its
source for what its consumer asked for.

The readers below it are not asked again, so the source reads no further
than the block in which the n-th variant fell. The VCF reader has read and
parsed the lines of that block, and the vars file reader may have
decompressed up to 7 batches past the one that held it: it decompresses 8
batches at once, that one among them, or as many as the threads of the
pool when there are fewer. The reader one block
ahead of `docs/specs/block.md` reads the whole chain on its thread, this
filter included, so when the filter gives `None` the thread sends that
word and ends, having asked the source for nothing more.

### How it is verified

Against bcftools 1.24 on `many.vcf`, read with every variant given:
`bcftools view -H many.vcf | head -n 10` gives the variants at 1000, 1037,
1074, 1111, 1148, 1185, 1222, 1259, 1296 and 1333 of chr1, and
`bcftools view -H -Q 0.8:major many.vcf | head -n 10`, the MAF filter of
0.8 before the first 10, gives 1037, 1074, 1111, 1148, 1222, 1259, 1296,
1333, 1370 and 1407, run on 5 October 2026. The cargo tests, made at
`next_block` of the reader of this filter over a `VcfReader` on `many.vcf`,
in blocks of 7 variants and of the default size, assert the ten positions
of each, and in blocks of 7 the counts 14 given and 10 kept and that the
source was asked for two blocks. With n of 14 over blocks of 7, where the
n-th variant ends a block, the source is asked for two blocks too. With the
MAF filter of 0.8 before the first 10, in blocks of 7, the MAF filter is
given 14 and keeps 12 and this filter is given 12 and keeps 10: the MAF
filter keeps 12 of the first 14 variants of `many.vcf`, by `bcftools view -H
-Q 0.8:major` over those 14 lines, run on 5 October 2026. The filter that keeps variants at random,
at 0.5 with a seed of 42, and then the first 10 gives 1037, 1074, 1111,
1148, 1222, 1296, 1370, 1407, 1555 and 1592, the first ten of the 243 of
its table.

That the pass ends and does not read the source to its end is tested on a
source that never ends: a reader built for the test that gives blocks of 7
variants for as long as it is asked, with the first 10 on it, returns with
10 variants and having given two blocks, read on one thread and through the
reader one block ahead. Over a VCF of 100000 variants of 10 individuals,
read in blocks of 100 variants through a `Read` that counts its bytes, a
pass with the first 100 reads less than a tenth of the file. Over the vars
file of the same variants in batches of 100, read in a pool of rayon of one
thread, a pass with the first 100 decompresses one batch alone, by the list
of the batches read that the vars file reader keeps for the tests.

The pytest tests, made at `filter_first_n`, on `many.vcf` with every
variant given, assert the ten positions on
`many.vcf` and a `pass_stats` with 10 variants and `stopped_early` true; a
`num_vars` of 1000 on the same file, which gives the 500 and
`stopped_early` false; a `do_pca_from_variants` with the first 50, whose
second pass gives the first; the `ValueError` of 0 and of -1 and the
`TypeError` of 1.5 and of `True`; the `ValueError` of each step that takes
variants out added after it, and the filter of individuals accepted after
it. The TypeScript test asserts the ten positions, `stoppedEarly`, and the
`Error` of 0 and of the MAF filter added after it.

## The filter of the variants that passed their FILTER

### What it gives

It keeps the variants whose FILTER column, in the VCF they were read from,
is `PASS` or a dot, and takes out the rest, as a step of a pass like the
other filters, so that how many variants failed their FILTER is in the
counts beside how many each other filter took out. The rule is the one of
`only_passed` of `docs/specs/io_vcf.md`, which reads the whole column: a
variant passed when its FILTER is `PASS` or `.`, which in a VCF says that
no filter was applied, and failed with anything else, `q10` or `PASS;q10`.
The two use the one function of the VCF reader, so they never disagree.

`only_passed` stays, and true by default, as the owner decided on 6
October 2026: a VCF opened with it gives no variant that failed, and this
filter then keeps them all. The two are for two uses. `only_passed` drops a
failed line before it is parsed, which is the faster when a user never
wants the failed variants; this filter is for a user who opens a VCF with
every variant, `only_passed=False`, and wants the failed ones taken out
with the other filters and counted, as popnei_web, the web application of popnei, does. The options not
taken were to change the default of `only_passed` to false, which changes
what every caller who gave none gets, and to drop it.

A block carries, for each variant, whether it passed, in the column
`passed` of `docs/specs/block.md`, which the VCF reader fills from FILTER
and the vars file of `docs/specs/io_vars.md` stores. A source that has no
such column is refused at the first block of a pass, with an error of its
own, and not taken as if every variant had passed: a vars file written
before format 1.2, or from a source that had no column of it, may hold
variants that failed, and keeping them would give a result with no sign
that they are there. The owner decided this on 6 October 2026; the option
not taken was to take a source with no such column as all passed. The
error says "the variants hold no record of whether they passed their
FILTER, so the filter of the variants that passed cannot run on them: a
vars file holds it from format 1.2, written from a VCF", and in Python its
message starts with the path of the file.

### In Python and in TypeScript

```python
Variants.filter_passed() -> None
```

It is a step: it adds itself to the `Variants`, reads no variant and
returns nothing. Its kind is `"passed"`, and its `args` are `{}`. A second
filter of this kind is refused, and so is this one after the filter of the
first n, since it takes variants out. It reads no genotype and no other
filter depends on it, so it can come anywhere among the steps. Like every
filter of variants it answers no to `skip_outside`, so a filter by regions
after it no longer has the source skip what is outside the regions: on
`big.vcf`, plain, on one thread, that skip makes a pass 0.039 s against 0.54 s for
the whole read, by "Speed" of the filter by regions. So the docstring says
to add it first, after the filter by regions when there is one, so that
the counts of the filters after it are of the variants that passed. In TypeScript, `variants.filterPassed()`.

pyNei has no such filter: it reads every variant whatever its FILTER, as
"The VCF reader" of `docs/specs/io_vcf.md` says.

### How it runs

It asks its source for what its consumer asked for and for `PASSED` of
`docs/specs/variant.md`, keeps in place the rows whose `passed` is true,
with `retain_vars`, and adds to its counts. A block with no `passed`
column is the error above, with the block as it was and nothing counted.

### How it is verified

Against bcftools 1.24 on `many.vcf` opened with every variant:
`bcftools view -H -f .,PASS many.vcf` gives 475 of its 500 variants, those
whose FILTER is `PASS` or `.`, and its first ten are at 1000, 1037, 1074,
1111, 1148, 1185, 1222, 1296, 1333 and 1370 of chr1; the first that fails
is at chr1 1259. With `-Q 0.8:major` after it, the MAF filter of 0.8, it
gives 364. Run on 6 October 2026. On `many.vcf` the strict rule and that of
bcftools agree, since no FILTER there names two filters.

The cargo tests, made at `next_block` of the reader of this filter over a
`VcfReader` of `many.vcf` with `only_passed` false, in blocks of 7 and of
the default size: the 475 positions are those of `many.bcftools.tsv`
whose FILTER is `PASS` or `.`, the first ten the ones above; the counts are
500 given and 475 kept; with the MAF filter of 0.8 after it, the MAF filter
is given 475 and keeps 364; and the variants of all its blocks taken
together are those of a `VcfReader` with `only_passed` true, which cuts its
blocks elsewhere, since it drops a failed line before the line is in a
block. Over a vars file written in the test from that reader, the same 475.
Over a vars file written in the test from a reader built for it whose
blocks have no `passed` column, which is what a file of 1.0 or 1.1 holds,
the error of a source with no record, at the first block. The blocks of the
VCF reader with `PASSED` asked for and `only_passed` false have `passed`
false at the 25 variants whose FILTER is `q10`, and with `only_passed` true
none; `retain_vars` and `reblock` carry the column with the others. Over `cases.vcf` with
`only_passed` false, the variants at 100, 300 and 400, as the table of
`docs/specs/io_vcf.md` has them.

The pytest tests, made at `filter_passed`: on `many.vcf` opened with
`only_passed=False`, the 475 positions and a `pass_stats` of the kind
`"passed"` with 500 and 475; the same opened with the default, 475 given
and 475 kept; through `write_vars` and `open_vars`, the same 475; on a
vars file that pyarrow rewrote without the `passed` column and with the
version 1.1, a `ValueError` whose message starts with its path; a second
filter of this kind, and this one after `filter_first_n`, refused. The
TypeScript test asserts the 475 positions, the counts, and the error of a
second filter of this kind.

## The Rust interface

What a filter compares, with the largest value that keeps the variant.
The first three compare one number of the variant alone; the fourth
compares the variant with the ones kept before it. The kind is the key
the counts have in Python.

```rust
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VarFilteringCriterion {
    /// Missing genotypes divided by all the individuals.
    MaxMissingRate(f64),
    /// The count of the commonest allele divided by the called alleles.
    MaxMaf(f64),
    /// Heterozygous genotypes divided by the called genotypes.
    MaxObsHet(f64),
    /// The largest r² a variant may have against a variant kept no more
    /// than `max_dist` base pairs behind it on its chromosome.
    MaxLdR2 { max_allowed_r2: f64, max_dist: u64 },
}
impl VarFilteringCriterion {
    /// "missing_data", "maf", "obs_het" or "ld".
    pub fn kind(&self) -> &'static str;
    /// The largest value of the number that keeps the variant, whichever
    /// of the four it is, which for `MaxLdR2` is its `max_allowed_r2`. A
    /// binding crate reads it for the `args` of the step it shows the
    /// user.
    pub fn threshold(&self) -> f64;
    /// The `max_dist` of `MaxLdR2`, the second value of the `args` of its
    /// step, and None for the three that compare one number of a variant.
    pub fn max_dist(&self) -> Option<u64>;
}
```


How many variants a filter was given and how many it kept. It is declared
here, and the `block` module, whose trait gives it, uses it: two modules of
one crate that use each other, as `variant` and `block` do.

```rust
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FilteringStats {
    pub vars_processed: u64,
    pub vars_kept: u64,
}
```

The filter of one pass. It is an object of its own, apart from the reader,
so that the rule and its counts are tested on a block with no source.

```rust
pub struct VarFilter { /* private */ }
impl VarFilter {
    /// A threshold that is NaN, below 0 or above 1 is an error that names
    /// the criterion and the value.
    pub fn new(criterion: VarFilteringCriterion) -> Result<VarFilter>;
    pub fn criterion(&self) -> VarFilteringCriterion;
    /// It keeps the variants of the block that pass, in place, and adds to
    /// the counts. A block that does not pass `check`, or that has
    /// variants and no genotypes, is an error, the block is as it was and
    /// nothing is added. A block of no variants is left as it is.
    pub fn filter_block(&mut self, block: &mut Block) -> Result<()>;
    /// Over every block it was given since it was built.
    pub fn stats(&self) -> FilteringStats;
}
```

The reader over a reader, with the contract of "How it runs".

```rust
pub struct FilteredReader<R: BlockReader> { /* private */ }
impl<R: BlockReader> FilteredReader<R> {
    /// An error when `reader` already has a filter of the kind of
    /// `filter`, which its `filtering_stats` says.
    pub fn new(reader: R, filter: VarFilter) -> Result<FilteredReader<R>>;
}
impl<R: BlockReader> BlockReader for FilteredReader<R> { /* ... */ }
```

The steps of a pass. Each binding crate keeps its list of steps as a
list of these, and `chain_of` below takes them. With them, the names of
the kept individuals as indices, which the binding crates call at the
call of the method, against the individuals of the source, and the
reader of the filter calls when the chain is built. Its errors are new
cases of the error of the crate, each a `ValueError` in Python: a name
that is not an individual, with the name; a name twice, with the name;
and no name. A fourth new case is a second filter of individuals on a
`Variants` that has one, whose message says so and names the kind; the
case of a second threshold filter, which carries two thresholds, stays
as it is. With them also the individuals the next pass gives, which each
binding crate answers with when a user reads the individuals of their
`Variants`: which step says who the next pass holds is of the filters and
not of Python or of TypeScript.

```rust
#[non_exhaustive]
pub enum PassStep {
    VarFilter(VarFilteringCriterion),
    /// The names of the individuals to keep, in the order to keep them.
    KeepIndividuals(Vec<String>),
}
impl PassStep {
    /// "missing_data", "maf", "obs_het" or "individuals".
    pub fn kind(&self) -> &'static str;
}

/// The index of each of `names` among `individuals`, in the order of
/// `names`.
pub fn resolve_individuals(names: &[String], individuals: &[String]) -> Result<Vec<usize>>;

/// The names of the individuals the next pass gives, in its order: the
/// names the last `KeepIndividuals` of `steps` keeps, and `of_the_source`
/// when no step is one.
pub fn individuals_of(steps: &[PassStep], of_the_source: &[String]) -> Vec<String>;

pub struct IndividualsReader<R: BlockReader> { /* private */ }
impl<R: BlockReader> IndividualsReader<R> {
    /// What `resolve_individuals` refuses against `reader.individuals()`.
    pub fn new(reader: R, individuals: &[String]) -> Result<IndividualsReader<R>>;
}
impl<R: BlockReader> BlockReader for IndividualsReader<R> { /* ... */ }
```

The method of `Block`, of `docs/specs/block.md`, that the reader compacts
each block with, which that spec leaves to this item.

```rust
impl Block {
    /// It keeps the genotypes of the individuals `keep`, indices into the
    /// individuals of the block, in that order, within the array of the
    /// block, and sets `num_individuals`. An error, with the block as it
    /// was, for an index at or beyond the individuals, an index twice, no
    /// index, and a block with variants and no genotypes.
    pub fn retain_individuals(&mut self, keep: &[usize]) -> Result<()>;
}
```

Its first three refusals are three new cases of the error of the crate.
They mark a defect of popnei, a `RuntimeError` in Python, and not a wrong
input of a user, as the `keep` of `retain_vars` that has not one value for
each variant of its block does: the indices come from
`resolve_individuals`, which refuses the name behind each of them, so no
call of a user reaches them. Each carries what finds the cause, the index
and the individuals of the block for an index at or beyond them, and the
index for one that is there twice. A block with variants and no genotypes
is the error of a field that is not in the block, the one
`VarFilter::filter_block` gives for the same block, and a block whose
arrays are not of the size it states is the error `Block::check` finds,
since the rows are cut out of the genotypes by those sizes.

The filter by linkage disequilibrium, which is a type of its own and not
a `VarFilter`: it holds the window of `docs/specs/ld.md` between one
block and the next, where a `VarFilter` reads each block on its own and
keeps nothing but its two counts.

```rust
pub struct LdFilter { /* private */ }
impl LdFilter {
    /// A `max_allowed_r2` that is NaN, below 0 or above 1, and a
    /// `max_dist` below 1, are an error that names the argument and the
    /// value.
    pub fn new(max_allowed_r2: f64, max_dist: u64) -> Result<LdFilter>;
    pub fn criterion(&self) -> VarFilteringCriterion;
    /// It keeps the variants of the block that pass, in place, and adds
    /// to the counts. The window carries over, so a block is filtered
    /// against the variants of the blocks before it. A block that does
    /// not pass `check`, or that has variants and no genotypes or no
    /// position, is an error, the block is as it was, nothing is added to
    /// the counts and the window is as it was.
    pub fn filter_block(&mut self, block: &mut Block) -> Result<()>;
    /// Over every block it was given since it was built.
    pub fn stats(&self) -> FilteringStats;
}

pub struct LdFilteredReader<R: BlockReader> { /* private */ }
impl<R: BlockReader> LdFilteredReader<R> {
    /// An error when `reader` already has a filter by linkage
    /// disequilibrium, which its `filtering_stats` says.
    pub fn new(reader: R, filter: LdFilter) -> Result<LdFilteredReader<R>>;
}
impl<R: BlockReader> BlockReader for LdFilteredReader<R> { /* ... */ }
```

Its `set_needs` asks its source for the chromosome and the position
besides the genotypes, whatever its consumer asked for, and the blocks it
gives hold all three.

The chain of the filters of one pass, which both binding crates build with
this function and neither writes itself. It takes the source of the pass
and gives the outermost reader of the chain, so that whoever started the
pass holds it: they read the counts from it when the consumer returns, and
lend it, `&mut`, to a consumer that takes a reader, as `write_vars` of
`docs/specs/io_vars.md` does. The fields the consumer wants are set on what
it gives, once the chain is built, and every filter of the chain passes
them on with the genotypes added.

```rust
/// One reader over `reader` for each step, in their order, so that each
/// filter sees what the one before it kept: a `FilteredReader` for a
/// `VarFilter` step, an `LdFilteredReader` for a `MaxLdR2` one and an
/// `IndividualsReader` for a `KeepIndividuals` one. No step gives
/// `reader` as it is.
///
/// # Errors
///
/// What `VarFilter::new` and `LdFilter::new` refuse, a threshold that is
/// not a number from 0 to 1 and a `max_dist` below 1; what
/// `IndividualsReader::new` refuses; and a step of the kind of one before
/// it in `steps` or, for a threshold filter or the filter by linkage
/// disequilibrium, of a filter that `reader` holds already, which
/// `FilteredReader::new` and `LdFilteredReader::new` refuse.
pub fn chain_of(
    reader: Box<dyn BlockReader>,
    steps: &[PassStep],
) -> Result<Box<dyn BlockReader>>;
```

The owner decided on 21 September 2026 that the core builds the chain. Each
binding crate had the loop that puts one filter over another, and in which
order the filters go, and which errors a user gets while they are built,
are the same in both languages and are of the filters and not of either
language. What stays in a binding crate is its list of steps, the names of
its arguments and the two refusals it gives a user at the call of a method.
The option not taken was that each binding crate builds its own chain, as
the code had it.

The method this spec adds to `BlockReader`, of `docs/specs/block.md`, which
the VCF reader, the vars file reader and `reblock` implement too.

```rust
    /// The kind and the counts of every filter between this reader and
    /// its source, this one first when it is a filter. A source gives
    /// none.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)>;
```

This module adds five cases to the error of the crate. Four are the
wrong input of a function, which both binding crates give their user as
such, a `ValueError` in Python: a threshold out of range, a `max_dist`
below 1, a second filter of one kind, with the kind, and a variant whose
position does not rise within its chromosome, which names the variant and
both positions and which the filter by linkage disequilibrium is the one
reader of popnei to refuse.

The fifth is a defect of a caller of the core crate and not of a user,
and so a `RuntimeError`: the plain filter of a threshold, `VarFilter`,
built for the criterion of the filter by linkage disequilibrium, which
that filter does not answer because whether a variant is kept turns on
the variants kept before it and not on the variant alone. Nothing a user
writes reaches it; `chain_of` sends that criterion to the reader that
does answer it.

A user has to get the refusal of a second filter of one kind when they
call the method that adds the filter, and no reader exists then, so the
binding crate looks for the kind among the steps of the `Variants` with
this function, which both crates call as they call `chain_of`: which
filters can stand together is of the filters and not of Python or of
TypeScript.

```rust
/// The error of a second filter of one kind when `new` is of the kind of
/// one of `set`, the steps that are set already. For a threshold filter
/// the error carries both thresholds, the one of `new` and the one that
/// is set, which a chain of readers cannot say and the steps can; for the
/// filter of individuals it carries the kind. Of the crate alone since 5
/// October 2026: the binding crates call `refuse_a_step`, below, which
/// calls this one.
pub(crate) fn refuse_a_second_filter_of_a_kind(
    set: &[PassStep],
    new: &PassStep,
) -> Result<()>;
```

The filter by regions. The regions of a BED file, joined where they
overlap or touch and sorted within each chromosome, are read once at the
call of the method and shared, behind an `Arc`, by the step and by every
pass built from it. A source that is handed them reads them too.

```rust
pub struct Regions { /* private */ }
impl Regions {
    /// The regions of the BED text of `source`, plain or gzipped. An
    /// error, with the line, for each of the four refusals of "In Python
    /// and in TypeScript" of this filter, and an error of the input when
    /// `source` cannot be read.
    pub fn from_bed<S: Read>(source: S) -> Result<Regions>;
    /// How many regions are left once those that overlap or touch are
    /// joined.
    pub fn num_regions(&self) -> usize;
    /// Whether position `pos` of the chromosome named `chrom` is in a
    /// region.
    pub fn contains(&self, chrom: &str, pos: u64) -> bool;
}

/// Which variants the filter keeps: those inside the regions, or, with
/// `exclude`, those outside all of them.
#[derive(Debug, Clone)]
pub struct RegionSelection {
    pub regions: Arc<Regions>,
    pub exclude: bool,
}
impl RegionSelection {
    /// "regions", or "excluded_regions" with `exclude`.
    pub fn kind(&self) -> &'static str;
    pub fn keeps(&self, chrom: &str, pos: u64) -> bool;
    /// Whether no position from `min_pos` to `max_pos` of `chrom`, both
    /// included, is one this selection keeps, which is what lets the vars
    /// file reader skip a batch.
    pub fn keeps_none_of(&self, chrom: &str, min_pos: u64, max_pos: u64) -> bool;
}

pub struct RegionFilter { /* private */ }
impl RegionFilter {
    pub fn new(selection: RegionSelection) -> RegionFilter;
    pub fn selection(&self) -> &RegionSelection;
    /// It keeps the variants of the block that the selection keeps, in
    /// place, and adds to the counts; `chroms` names the chromosome
    /// numbers of the block. A block that does not pass `check`, or that
    /// has variants and no chromosome or no position, is an error, the
    /// block is as it was and nothing is added.
    pub fn filter_block(&mut self, block: &mut Block, chroms: &ChromTable) -> Result<()>;
    /// Over the blocks it was given. `RegionsReader` adds the
    /// `num_skipped` of its source to `vars_processed` when its counts are
    /// read.
    pub fn stats(&self) -> FilteringStats;
}

pub struct RegionsReader<R: BlockReader> { /* private */ }
impl<R: BlockReader> RegionsReader<R> {
    /// An error when `reader` already has a filter of the kind of
    /// `filter`. It offers the selection to `reader` with `skip_outside`,
    /// and asks it for the chromosome and the position besides what its
    /// consumer asks for.
    pub fn new(reader: R, filter: RegionFilter) -> Result<RegionsReader<R>>;
}
impl<R: BlockReader> BlockReader for RegionsReader<R> { /* ... */ }
```

The step, one more member of `PassStep`, for which `chain_of` builds a
`RegionsReader` and `refuse_a_second_filter_of_a_kind` refuses a second
step of the same kind:

```rust
    /// The variants the selection keeps.
    Regions(RegionSelection),
```

The two methods this filter adds to `BlockReader` are in
`docs/specs/block.md`: `skip_outside`, which a source takes and a filter
of variants refuses, and `num_skipped`. The cases it adds to the error of
the crate are four. Three are a `ValueError` in Python: a wrong line of a
BED file, with the number of the line and what is wrong with it, and a
BED with no region, in front of whose message the binding crate puts the
path of the file, as it does for a VCF; and a second filter by regions of
one kind, with the kind, which names no file, as the second threshold
filter does not. Neither case of a second filter that there is fits it:
one carries two thresholds and the other says that two lists of
individuals are one. The fourth is a defect, a `RuntimeError`: a block
given to `RegionFilter` whose chromosome number the table of its reader
has no name for, which the two writers refuse in the same words, since
the regions are looked up by the name.

The filter that keeps variants at random. The seed a step gets when the
user gives none is a constant of the core, which both packages read, so
that the default is written once.

```rust
/// 42, the seed of the filter that keeps variants at random when the
/// user gives none, which the owner chose on 5 October 2026.
pub const DEFAULT_RANDOM_FILTER_SEED: u64 = 42;

pub struct RandomFilter { /* private */ }
impl RandomFilter {
    /// A `keep_rate` that is NaN, below 0 or above 1 is an error that
    /// names the argument and the value. The generator starts at `seed`.
    pub fn new(keep_rate: f64, seed: u64) -> Result<RandomFilter>;
    pub fn keep_rate(&self) -> f64;
    pub fn seed(&self) -> u64;
    /// It draws one number for each variant of the block, in order, keeps
    /// the variants whose number is below the keep rate, in place, and
    /// adds to the counts. A block that does not pass `check` is an error,
    /// the block is as it was, nothing is added and no number is drawn.
    pub fn filter_block(&mut self, block: &mut Block) -> Result<()>;
    /// Over every block it was given since it was built.
    pub fn stats(&self) -> FilteringStats;
}

pub struct RandomlyFilteredReader<R: BlockReader> { /* private */ }
impl<R: BlockReader> RandomlyFilteredReader<R> {
    /// An error when `reader` already has a filter of this kind.
    pub fn new(reader: R, filter: RandomFilter) -> Result<RandomlyFilteredReader<R>>;
}
impl<R: BlockReader> BlockReader for RandomlyFilteredReader<R> { /* ... */ }
```

The filter that keeps the first n variants, and the function that tells
whether it ended a pass, which both binding crates call when they read the
counts of a pass, so that the rule is written in the core alone.

```rust
pub struct FirstNReader<R: BlockReader> { /* private */ }
impl<R: BlockReader> FirstNReader<R> {
    /// An error when `num_vars` is 0 and when `reader` already has a
    /// filter of this kind.
    pub fn new(reader: R, num_vars: u64) -> Result<FirstNReader<R>>;
}
impl<R: BlockReader> BlockReader for FirstNReader<R> { /* ... */ }

/// Whether a filter of the first n among `steps` kept its n variants and
/// so ended the pass whose counts are `filtering`, as
/// `BlockReader::filtering_stats` of the outermost reader of the chain
/// gives them: true when `steps` has a `FirstN(n)` and the counts under
/// "first_n" have `vars_kept` equal to n. `steps` are those the pass was
/// built from, which a binding crate keeps from the start of the pass, and
/// not those of the `Variants` when the counts are read: a step added
/// while an `iter_blocks` runs is not in its pass.
pub fn stopped_early(steps: &[PassStep], filtering: &[(&'static str, FilteringStats)]) -> bool;
```

The two steps, two more members of `PassStep`, for which `chain_of` builds a
`RandomlyFilteredReader` and a `FirstNReader`, and a third, for which it
builds a `PassedReader`:

```rust
    /// Each variant kept when the number drawn for it is below `keep_rate`;
    /// of the kind "random".
    Random { keep_rate: f64, seed: u64 },
    /// The first `num_vars` variants, and then the pass ends; of the kind
    /// "first_n".
    FirstN(u64),
    /// The variants whose FILTER was `PASS` or a dot; of the kind "passed".
    Passed,
```

The reader of the filter of the variants that passed. It is given whole
blocks and keeps their rows with `passed` true.

```rust
pub struct PassedReader<R: BlockReader> { /* private */ }
impl<R: BlockReader> PassedReader<R> {
    /// An error when `reader` already has a filter of this kind.
    pub fn new(reader: R) -> Result<PassedReader<R>>;
}
impl<R: BlockReader> BlockReader for PassedReader<R> { /* ... */ }
```

`refuse_a_second_filter_of_a_kind` refuses a second step of any of the three kinds,
with the kind, as it does for the filter by regions. A step that takes
variants out after a `FirstN` is refused by a function of its own. Both are
of the crate alone, and `refuse_a_step` calls them in their order: both
binding crates call it when a user adds a step, and `chain_of` for each step
as it builds the chain, so that no binding can call one and not the other:

```rust
/// The error of a step that takes variants out, `new`, after a filter of
/// the first n among `set`, the steps that are set already. The filter of
/// individuals is not refused.
pub(crate) fn refuse_a_step_after_the_first_n(set: &[PassStep], new: &PassStep) -> Result<()>;

/// Every refusal of `new` against `set`: a second filter of its kind first,
/// then a step after the filter of the first n. `chain_of` and both binding
/// crates call it, so that one step gets one error in every language.
pub fn refuse_a_step(set: &[PassStep], new: &PassStep) -> Result<()>;

/// The step of the filter of the first n, or the error of a `num_vars` of
/// 0, which names the argument; both binding crates build the step with it.
pub fn first_n_step(num_vars: u64) -> Result<PassStep>;
```

Every consumer of both binding crates gives the Python or the TypeScript
package this field with the counts of its pass, and builds it with this
function.

The cases the two filters add to the error of the crate, each a
`ValueError` in Python: a keep rate out of range, with the value; a
`num_vars` of 0; a step after the filter of the first n, with the kind of
the step; and a second filter of either kind, with the kind. The filter of
the variants that passed adds two: a block with no `passed` column, of its
own case and among those that name the file, so that in Python its message
starts with the path of the source, as the errors of a vars file do; and a
second filter of this kind, as for the filter that keeps variants at
random. Both are a `ValueError` in Python.

## Speed

### The three threshold filters

There is no number to reach: what each filter costs was measured, and what
it should cost is the owner's to set. The numbers are of 21 September
2026, on the owner's Apple M5 Pro, 18 cores, macOS 27.0, native
`aarch64-apple-darwin`, a build of `cargo bench`, over the 400 MB VCF of
`docs/rust_core.md`, 100000 variants of 1000 individuals whose genotypes
are missing at a rate of 0.03, and over the vars file popnei writes of it,
both already in the page cache. A number is of 5 runs of a whole pass with
the genotypes alone asked for, timed by
`crates/popnei/benches/filter_vars.rs`. What the filter costs is that pass
less the pass with no filter, the two run back to back, and each pair was
run four times: every cell is the range the medians of the four sets
spread over, so each threshold stands against the pass with no filter of
its own pairs.

| | the pass with no filter | with the filter | what the filter costs |
|---|---|---|---|
| the VCF, 1 thread, at 0.1 | 0.558 to 0.569 s | 0.617 to 0.639 s | 0.058 to 0.070 s |
| the VCF, 1 thread, at 0.03 | 0.548 to 0.573 s | 0.616 to 0.646 s | 0.049 to 0.085 s |
| the VCF, 18 threads, at 0.1 | 0.095 s | 0.101 to 0.103 s | 0.006 to 0.008 s |
| the VCF, 18 threads, at 0.03 | 0.095 to 0.097 s | 0.107 to 0.108 s | 0.011 to 0.013 s |
| the vars file, 1 thread, at 0.1 | 0.102 s | 0.158 to 0.159 s | 0.056 to 0.057 s |
| the vars file, 1 thread, at 0.03 | 0.102 s | 0.162 s | 0.060 s |
| the vars file, 18 threads, at 0.1 | 0.102 s | 0.115 to 0.116 s | 0.013 to 0.014 s |
| the vars file, 18 threads, at 0.03 | 0.102 to 0.103 s | 0.119 to 0.120 s | 0.017 to 0.018 s |

At 0.1 the filter keeps every one of the 100000 variants and at 0.03 it
keeps 54773, so the difference between the two thresholds is what
compacting the blocks costs, 0.003 to 0.005 s in the rows where the sets
are stable enough to show it; over the VCF on one thread the sets spread
by more than that and the two thresholds cannot be told apart there. The
vars file is read in 0.102 s on one thread and on 18, since its reader
runs on the thread that calls it and only the filter reads the rows of a
block on the pool.

What is measured is the filter against a pass that reads no genotype of
its own, so for a calculation that reads them after the filter these
numbers are an upper bound on what the filter adds: how much less it is
was not measured.

pyNei's pass over the same VCF takes 14.05 s and its
`filter_by_missing_data` costs it 0.415 s at 0.1 and 0.400 s at 0.03; that
pass also builds a pandas frame of the chromosome, the position, the id
and the quality, and the alleles beside it, for every chunk.
`bcftools view -H`, with its records sent to `/dev/null` in both passes,
takes 1.431 s and its `-i "F_MISSING<=0.1"` costs it 0.119 s; that pass
writes back out as text, on one thread, every record it keeps, which is
why at 0.03 it takes 0.202 s less with the filter than without, writing
54773 records instead of 100000. So the difference of each program against
itself is like for like and the three whole passes are three different
pieces of work. The three keep the same variants, 100000 at 0.1 and 54773
at 0.03. `docs/reports/filters-measurement.md` has every set of runs with
the load average it was taken at, how pyNei's difference was told apart
from the drift of the machine over a pass of 14 s, and what it leaves to a
performance review.

### The filter by linkage disequilibrium

Measured on 23 September 2026 by the performance review of
`docs/reports/perf-ld-2026-09-23.md`, on the same machine and the same two
files as the table above, with
`crates/popnei/benches/filter_vars.rs --max-ld-r2 --max-dist`. A number is
the median of 5 runs of a whole pass with the genotypes alone asked for,
with the best and the worst beside it, taken with nothing else running.

**What this filter costs is set by how many variants its window holds, and
that is set by the dataset and by `max_dist`, not by the filter.** The 400
MB VCF has its 100000 variants 1000 base pairs apart on two chromosomes, so
`max_dist` sets the window almost exactly, and its variants carry no
linkage, so at a threshold of 0.3 the filter keeps 99919 of them and every
window stays full. That is the worst case for this filter: a dataset whose
variants are linked leaves a short window and costs less, which the last
two rows show. Every number below therefore states what the window held,
and a number without one says nothing.

Over the vars file, one thread, against the same pass with no filter, which
takes 0.106 s:

| threshold | max_dist | what the window held | kept of 100000 | the pass | what the filter costs |
|---|---|---|---|---|---|
| 0.3 | 1 | 0.0 variants | 100000 | 0.696 s | 0.590 s |
| 0.3 | 10000 | 10.0 | 99994 | 0.840 s | 0.734 s |
| 0.3 | 50000 | 50.0 | 99974 | 0.991 s | 0.885 s |
| 0.3 | 250000 | 249.2 | 99919 | 1.631 s | 1.525 s |
| 0.3 | 1000000 | 988.5 | 99850 | 4.184 s | 4.078 s |
| 0.02 | 250000 | 88.7 | 35673 | 1.107 s | 1.001 s |
| 0.005 | 250000 | 13.7 | 5759 | 0.849 s | 0.743 s |

The first row is `max_dist` of 1, where no variant is ever in another's
window: **0.590 s of the cost is paid before the window holds anything**.
That is the three matrices built for every variant of the pass and the r²
of every set of candidates against itself, neither of which depends on
`max_dist`. The three rows at a `max_dist` of 250000 and different
thresholds show the same thing from the other side: what changes the cost
is the window, 249.2 variants at 1.631 s and 13.7 at 0.849 s, and not how
many variants are kept.

Over the 400 MB VCF, at a threshold of 0.3 and a `max_dist` of 250000:

| | the pass with no filter | with the filter | what the filter costs |
|---|---|---|---|
| 1 thread | 0.600 s | 2.147 s | 1.547 s |
| 18 threads | 0.086 s | 1.655 s | 1.569 s |

**The thread column of this filter does not mean what it means for the
three above.** The three threshold filters read the rows of a block on the
pool, so their thread column is their own work spread over cores. This one
never does: its rule is sequential, since whether a variant is kept decides
what the variants after it are compared with, and its products go to
Accelerate, which does not split a product of this shape across cores. So
the filter costs the same 1.55 s on one thread and on eighteen, and what
the threads speed up is the parsing of the VCF, 0.600 s to 0.086 s. Over
the vars file, whose reader runs on the calling thread, the two columns
would be the same number. One set of 5 runs at 18 threads over the VCF gave
a worst of 2.388 s against a best of 1.629 s; the other sets were stable to
0.03 s.

The filter is 6 to 40 times the cost of the pass it filters, where the
three threshold filters are a tenth to a half of theirs. The work is
different: the three compare one number per variant, and this one takes the
r² of every candidate against every variant of its window, which is matrix
products over the genotypes.

In a browser, on the wasm package run under node v26.8.2 with
`js/popnei/bench/time_filter_by_ld.mjs`, on one thread, where the linear
algebra is faer and not Accelerate, over the same vars file, the median of
5 runs:

| | a browser | natively | how many times |
|---|---|---|---|
| the pass with no filter | 0.114 s | 0.106 s | 1.1 |
| with the filter, `max_dist` 250000, window 249.2 | 18.402 s | 1.631 s | 11.3 |
| with the filter, `max_dist` 1, window 0 | 7.483 s | 0.696 s | 10.8 |

Reading the file costs a browser almost nothing more, so the whole of the
difference is the filter's arithmetic. The share paid before the window
holds anything is the same on both, two fifths. There is no target for a
browser. The filter keeps the same variants in both builds, 99919 of 100000
at a `max_dist` of 250000 and 100000 at a `max_dist` of 1.

What was measured and not changed, from the same review: the window keeps
the genotypes of its variants, one byte per allele, and expands them into
the three matrices of `docs/specs/ld.md` for every set of candidates it
settles, which over a pass of 100000 variants at a window of 270 is 2.53 GB
of `f64` written. Keeping the three matrices in the window instead would
make it twelve times larger, 5 MB to 60 MB at a window of 2500 variants of
1000 individuals, and "How it runs" of this item chose the smaller window.
That trade is the owner's and is in the report.

### The filter by regions

The numbers to reach were worked out from the profile of the VCF reader
before the skip was built: on `big.vcf` of `docs/specs/io_vcf.md`, on one
thread, with a BED that keeps 1000 of its 100000 variants in one run of
the chromosome, handed to the source, the pass takes no more than 0.15 s
plain, against the 0.563 s of reading the whole file, and no more than
0.45 s bgzipped, against 0.890 s. In the plain file 92 in 100 of the read
is the columns of the individuals, which the skip does not parse, and
what is left, 0.045 s, is the part the skip keeps; in the bgzipped one the
decompression, 0.33 s by the difference of the two reads, is kept too. On
`big.vars`, whose batches are of 5000 variants, the same BED keeps the
variants of one batch, and the pass takes no more than 0.02 s, against the
0.102 s of the whole pass of `docs/specs/stats.md`.

The skip was measured on 27 September 2026, on the owner's Apple M5 Pro,
18 cores, release builds, the files in the page cache, at a load average
of the minute before each set of 3.2 to 3.5. The BED is the one line
`chr1 0 1000000`, which keeps the variants at positions 1000 to 1000000 of
`chr1`, 1000 of the 100000, all in the first batch of `big.vars`; the
bgzipped file is `big.vcf` bgzipped by bgzip, 38 MB. The reads of the VCF
are `crates/popnei/benches/read_vcf.rs --bed`, which hands the regions to
the reader and asks for the genotypes, five runs each; the passes through
the filter are `crates/popnei/benches/filter_vars.rs --bed`, which puts
the filter by regions on the pass with `chain_of`, as `filter_by_regions`
does, one pass not timed and five timed. Every pass gave 1000 variants,
and the filter was given the 100000 and kept those 1000.

| | the runs | the median | the target | met |
|---|---|---|---|---|
| `big.vcf` plain, the read with the regions, 1 thread | 0.042, 0.039, 0.039, 0.039, 0.038 s | 0.039 s | 0.15 s | yes |
| the same, the pass through the filter | 0.040, 0.038, 0.038, 0.038, 0.039 s | 0.038 s | 0.15 s | yes |
| `big.vcf` bgzipped, the read with the regions, 1 thread | 0.314, 0.312, 0.312, 0.312, 0.312 s | 0.312 s | 0.45 s | yes |
| the same, the pass through the filter | 0.312, 0.311, 0.312, 0.312, 0.312 s | 0.312 s | 0.45 s | yes |
| `big.vars`, the pass through the filter, 1 thread | 0.005 s in each of the five | 0.005 s | 0.02 s | yes |

On 18 threads the read with the regions takes 0.035, 0.035, 0.034, 0.034
and 0.034 s plain and 0.310, 0.310, 0.309, 0.308 and 0.307 s bgzipped,
and the pass over `big.vars` 0.005 to 0.006 s: what is left is the serial
pass over the lines, which reads their POS, and the decompression, and
neither runs on the pool. The whole of the bgzipped file is still
decompressed, which is why it takes eight times the plain one.

The skip adds the reading of POS to the serial pass of the reader, so the
read of the whole plain file with no regions was measured in the same
session against the targets of "Speed" of `docs/specs/io_vcf.md`: 0.583,
0.545, 0.542, 0.540 and 0.541 s on one thread, a median of 0.542 s against
0.594 s, and 0.082, 0.081, 0.082, 0.081 and 0.083 s on 18, a median of
0.082 s against 0.108 s. Both are met, and both are under the 0.563 and
0.093 s of the reader of 21 September. The whole bgzipped file took a
median of 0.849 s on one thread and 0.364 s on 18, against 0.924 and 0.44
s.

For comparison, in the same session, `bcftools view -H -T` with the same
BED, which reads the whole file and writes the 1000 lines to
`/dev/null`, took 0.78, 0.74 and 0.77 s plain and 0.84, 0.84 and 0.82 s
bgzipped, on one thread; `tabix` with the region `chr1:1-1000000` over the
index of the bgzipped file, which seeks to the region and reads none of
the rest, took 0.01 s or less in each of three runs.

### The filter that keeps variants at random and the filter of the first n

Neither has a number to reach. The filter that keeps variants at random
draws one number for each variant, a few integer operations, and compacts
the block as the threshold filters do; what it costs has not been measured.
The filter of the first n does nothing to a variant, and what it changes is
the time of a pass, which becomes the time of reading the source up to the
block that holds the n-th variant.

## Open points

None. The owner decided on 5 October 2026 the four points of the two filters
added that day, the seed of 42 by default, the draw of one number for each
variant in its order, the refusal of a step after the filter of the first
n and `stopped_early` in the counts, each written where it applies with the
option not taken. The owner decided on 26 September 2026 the one the filter by regions
had, that a variant is in a region by its position alone, which is
written under "What it gives" of that filter with the option not taken.
What the owner decided on 21 September 2026 about the three
threshold filters and the counts, and on 22 September 2026 about the
filter by linkage disequilibrium and about the order of the kept
individuals, in chat, is written where it applies, with the option that
was not taken. The two open points of the filter by linkage
disequilibrium that the owner has to answer are in
`docs/specs/ld.md`, because both are about the r² itself: whether the
curve of linkage disequilibrium against distance also carries a sample of
pairs, which this filter does not touch, and how the major allele of a
variant with half called genotypes is picked, which moves the r² this
filter compares.

## Not in this spec

- A sample of exactly n variants spread over the whole file, which has to
  read the whole file before it knows which to keep: not built. A user puts
  the filter that keeps variants at random before the filter of the first n.
- A new sample at each call, with no seed given: a user passes another
  seed.

- Chromosome names matched with and without `chr`, `chr1` with `1`: not
  done, under the cases of the filter by regions.
- The regions of an index of a bgzipped VCF, which tabix makes, so that
  the reader seeks to them: `docs/specs/io_vcf.md` does not plan it, and
  the skip of the filter by regions still decompresses the whole file.
- A lowest maf, or a threshold on the frequency of the minor allele: pyNei
  has neither, and popnei does not add them.
- The three filters over the individuals of one population and not over
  all of them: pyNei does not have it. A user puts the filter of
  individuals before them.
- A copy of a `Variants` with its steps: a user opens the source again.
- The read ahead thread of section 3 of the architecture: "The reader one
  block ahead" of `docs/specs/block.md`. It lends the chain of readers to
  its thread and gives it back when the pass ends, so that the counts of
  every filter of the chain are read from the chain itself; while the pass
  runs, the handle answers with the counts as of the last block it gave.
