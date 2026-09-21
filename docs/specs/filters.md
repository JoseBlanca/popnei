# The filters module: variants kept by missing data, major allele frequency and observed heterozygosity

September 2026. The `filters` module gives a user of popnei the variants of
a dataset that pass a threshold, before any calculation sees them: those
with few missing genotypes, those whose commonest allele is not too
frequent, and those with few heterozygous individuals. It also tells the
user how many variants each filter was given and how many it kept. There
is no code. This spec develops the row `filters` of the table in section 9
of `docs/architecture.md`, and it covers the three filters that compare one
number of a variant with a threshold, and the counts. The other filters of
that row, the filter of individuals and the filter by linkage
disequilibrium, are items that are not written. It depends on
`docs/specs/block.md`, which has the `Block`, the run of consecutive
variants held as arrays that the variants flow in, and the `BlockReader`
trait of everything that gives blocks, and on `docs/specs/variant.md`,
which has the counts of the genotypes and of the alleles of one variant,
the `Variants` that a user puts the filters on, and the `PassStats` in
which the counts reach them.

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
    kind: str                  # "missing_data", "maf" or "obs_het"
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
`"missing_data"`, `"maf"` or `"obs_het"`, to its `FilteringStats`, in the
order of the steps, and it is empty when the `Variants` had no filter.

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

A binding crate builds the chain of a pass and keeps it while the consumer
runs. When the consumer returns, the binding crate reads the counts from
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

## The Rust interface

Which number of a variant is compared, with the largest value of it that
keeps the variant. The kind is the key the counts have in Python.

```rust
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VarFilteringCriterion {
    /// Missing genotypes divided by all the individuals.
    MaxMissingRate(f64),
    /// The count of the commonest allele divided by the called alleles.
    MaxMaf(f64),
    /// Heterozygous genotypes divided by the called genotypes.
    MaxObsHet(f64),
}
impl VarFilteringCriterion {
    /// "missing_data", "maf" or "obs_het".
    pub fn kind(&self) -> &'static str;
    /// The largest value of the number that keeps the variant, whichever
    /// of the three it is. A binding crate reads it for the `args` of the
    /// step it shows the user.
    pub fn threshold(&self) -> f64;
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

The method this spec adds to `BlockReader`, of `docs/specs/block.md`, which
the VCF reader, the vars file reader and `reblock` implement too.

```rust
    /// The kind and the counts of every filter between this reader and
    /// its source, this one first when it is a filter. A source gives
    /// none.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)>;
```

This module adds two cases to the error of the crate, a threshold out of
range and a second filter of one kind, with the kind, which both binding
crates give their user as the wrong input of a function. A user has to
get the second one when they call the method that adds the filter, and no
reader exists then, so the binding crate looks for the kind among the
steps of the `Variants` and gives that error itself.

## Speed

There is no number to reach yet, and the measurement comes first: the
missing data filter at 0.1 over the 400 MB VCF of `docs/rust_core.md`,
100000 variants of 1000 individuals, and over its vars file, against the
same filter of pyNei and against `bcftools view -i "F_MISSING<=0.1"`, each
as the time of the whole pass less the time of the pass with no filter.

## Open points

None. What the owner decided on 21 September 2026, in chat, is written
where it applies, with the option that was not taken.

## Not in this spec

- The filter of individuals, pyNei's `filter_samples`, and taking
  individuals out of a block: a later item of this spec. pyNei ignores a
  name that no individual has, and gives the individuals in the order of
  the source whatever the order of the argument; that item decides both.
- The filter by linkage disequilibrium, pyNei's `filter_by_ld_and_maf`: a
  later item of this spec, after `docs/specs/ld.md`, which is not written,
  has the dosages and the r it compares. It keeps the last variant it kept
  from one block to the next, and it compares the absolute value of r, and
  not its square, with an argument called `min_allowed_r2`.
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
