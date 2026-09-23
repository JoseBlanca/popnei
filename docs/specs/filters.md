# The filters module: variants kept by missing data, major allele frequency and observed heterozygosity

September 2026. The `filters` module gives a user of popnei the variants of
a dataset that pass a threshold, before any calculation sees them: those
with few missing genotypes, those whose commonest allele is not too
frequent, those with few heterozygous individuals, and those that do not
repeat what a variant near them on the chromosome already said. It also
tells the user how many variants each filter was given and how many it
kept. There is code for the first three and the counts, and none for the
fourth. This spec develops the row `filters` of the table in section 9
of `docs/architecture.md`, and it covers the three filters that compare one
number of a variant with a threshold, the counts, and the filter that
takes out the variants that repeat what a variant before them said, which
was added on 22 September 2026. The other filter of that row, the filter
of individuals, is an item that is not written. It depends on
`docs/specs/block.md`, which has the `Block`, the run of consecutive
variants held as arrays that the variants flow in, and the `BlockReader`
trait of everything that gives blocks; on `docs/specs/variant.md`,
which has the counts of the genotypes and of the alleles of one variant,
the `Variants` that a user puts the filters on, and the `PassStats` in
which the counts reach them; and, for the filter by linkage
disequilibrium alone, on `docs/specs/ld.md`, which has the r² that filter
compares and the dosages it reads the genotypes as.

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
    kind: str                  # "missing_data", "maf", "obs_het" or "ld"
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
`"missing_data"`, `"maf"`, `"obs_het"` or `"ld"`, to its
`FilteringStats`, in the order of the steps, and it is empty when the `Variants` had no filter.

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
is the one part of popnei that needs the variants of each chromosome to
come together and in order of position. Two variants at one position are
allowed and are 0 apart, so each is in the other's window. A variant
whose position is below the one before it on the same chromosome, and a
variant on a chromosome that had already ended, are an error that names
the variant, its position and the one before it. The rest of popnei reads
a source in any order, and `docs/specs/io_vars.md` has a test that writes
a block of four variants that are not sorted, so this is the one reader
that refuses what the others take. Decided here: the alternative is to
subtract two positions that can run backwards, which on a `u64` wraps to
a distance of 18 million million million and puts the pair outside every
window without a word.

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
has read, and which is on that variant's chromosome, its dosages held as
the three matrices that `docs/specs/ld.md` builds for the products of r²,
with its chromosome and its position. A variant leaves the window when
the reader passes `max_dist` beyond it or reaches another chromosome. The
memory is 24 bytes for each individual and each variant of the window:
for 250 kept variants of 1000 individuals, 6 MB. The variants of a window
are unlinked to each other by construction, so a window holds few of
them; a window whose dosages the machine has not the memory for is the
error of `docs/specs/block.md` for the same case.

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
/// One reader over `reader` for each criterion, in their order, so that
/// each filter sees what the one before it kept: a `FilteredReader` for
/// the three that compare one number of a variant and an
/// `LdFilteredReader` for `MaxLdR2`. No criterion gives `reader` as it
/// is.
///
/// # Errors
///
/// What `VarFilter::new` and `LdFilter::new` refuse, a threshold that is
/// not a number from 0 to 1 and a `max_dist` below 1, and what
/// `FilteredReader::new` and `LdFilteredReader::new` refuse, a criterion
/// of the kind of one before it in `criteria` or of a filter that
/// `reader` holds already.
pub fn chain_of(
    reader: Box<dyn BlockReader>,
    criteria: &[VarFilteringCriterion],
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
/// one of `set`, the criteria of the filters that are set already. The
/// error carries both thresholds, the one of `new` and the one that is
/// set, which a chain of readers cannot say and the criteria can.
pub fn refuse_a_second_filter_of_a_kind(
    set: &[VarFilteringCriterion],
    new: VarFilteringCriterion,
) -> Result<()>;
```

## Speed

There is no number to reach: what the filter costs was measured, and what
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

The numbers above are of the three threshold filters. The filter by
linkage disequilibrium has none: what it costs is the products of r² of
each variant against the variants of its window, which
`docs/specs/ld.md` measures at 1.9 ms for the r² of 256 variants against
another 256 over 1000 individuals, and how many variants a window holds
depends on the
dataset and on `max_dist`. The implementation plan that builds it, under
`docs/plans/`, measures a whole pass,
on `tests/reference/ld/ld.vcf.gz` and on the 400 MB VCF of the table
above, and this section gets the numbers.

## Open points

None. What the owner decided on 21 September 2026 about the three
threshold filters and the counts, and on 22 September 2026 about the
filter by linkage disequilibrium, in chat, is written where it applies,
with the option that was not taken. The two open points of the filter by
linkage disequilibrium that the owner has to answer are in
`docs/specs/ld.md`, because both are about the r² itself: whether the
curve of linkage disequilibrium against distance also carries a sample of
pairs, which this filter does not touch, and how the major allele of a
variant with half called genotypes is picked, which moves the r² this
filter compares.

## Not in this spec

- The filter of individuals, pyNei's `filter_samples`, and taking
  individuals out of a block: a later item of this spec. pyNei ignores a
  name that no individual has, and gives the individuals in the order of
  the source whatever the order of the argument; that item decides both.
- The variants of a region of a chromosome, which
  `docs/specs/io_vars.md` leaves to this spec: a later item.
- A lowest maf, or a threshold on the frequency of the minor allele: pyNei
  has neither, and popnei does not add them.
- The three filters over the individuals of one population and not over
  all of them: pyNei does not have it. A user puts the filter of
  individuals before them.
- A copy of a `Variants` with its steps: a user opens the source again.
- The read ahead thread of section 3 of the architecture: with the first
  calculation that consumes blocks. It takes the reader into its thread,
  and it gives it back when the pass ends, so that the counts can be read.
