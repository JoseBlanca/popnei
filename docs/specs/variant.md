# The variant module: the fields, the chromosome table and one variant of a block

September 2026. The `variant` module gives popnei what every other module
says about a variant, one site of the genome with the genotype of every
individual at it: which of its fields a consumer wants, the table that
turns the name of a chromosome into a number, the missing allele, and the
view of one variant of a block. This spec develops the row `variant` of the
table in section 9 of `docs/architecture.md`. It depends on
`docs/specs/block.md`, which has the `Block`, the run of consecutive
variants held as arrays that the variants flow in, and the view here is a
view into one. It also covers the `Variants` handle that a Python or a
TypeScript user holds, and the counts of one variant, of its genotypes and
of its alleles, which the filters of `docs/specs/filters.md` are the first
to use, and the `Variants` that is built from an array of genotypes in
memory. The other row helpers of the module, the dosages and the missing
and het masks of one variant, are items that are not written yet.

There is code, built from the first version of this spec, which had a
`Variant` that the caller owned and a reader filled, and the
`VariantReader` trait. The owner dropped both on 20 September 2026, for the
reasons at the end of section 1 of the architecture: the variants flow in
blocks, and the trait that gives them, `BlockReader`, is in
`docs/specs/block.md`. `Needs`, `ChromTable` and the two constants stay as
they were built.

## The fields, the chromosome table and the view

### What it gives

A consumer, a filter, a calculation or a writer, says which fields of the
variants it wants with a `Needs`, a set of five flags: the genotypes, the
chromosome and the position, which travel together, the id, the alleles
and the quality. A reader may skip the rest: most calculations want the
genotypes alone, and the VCF reader then does not parse the other columns
and the vars file reader does not decompress them. The quality is the QUAL
column of a VCF, phred scaled: minus ten times the base ten logarithm of
the probability that there is no variant at that site, so 30 is one in a
thousand.

The chromosome of a variant is a number, an index into the table of
chromosome names that its reader keeps, the `ChromTable`, so that a block
holds one `u32` per variant and not a text.

An allele is one `i8`: 0 is the reference allele, 1 and above the
alternative ones in the order of the VCF, up to 127, and `MISSING_ALLELE`,
-1, is an allele that was not called. There is no mask.

A calculation that works variant by variant walks the variants of a block
through `VariantRef`, a view of one of them: its genotypes as a slice of
the array of the block, individual after individual, `ploidy` alleles for
each, so the alleles of individual i are `gts[i * ploidy .. (i + 1) *
ploidy]`, and its other fields as they are in the columns. A view
allocates nothing and copies nothing.

### What a Python and a TypeScript user see

Neither sees a variant. pyNei's `Variants` yields chunks, a few thousand
variants as arrays, and its calculations are written over them in Python.
In popnei the loop of every filter and every calculation runs inside the
core, and what a user holds is a handle.

`Variants`, in Python, is that handle: a source of variants, a VCF path
with its options or a vars file, and the steps that were put on it, in
order. It holds no genotypes, and it is a recipe that is run later. A user
gets one from `open_vcf`, of `docs/specs/io_vcf.md`, or from `open_vars`,
of `docs/specs/io_vars.md`. What is done with it is of two kinds, and what
a call returns shows which. A step, a filter of `docs/specs/filters.md`,
is a method of the `Variants` that adds itself to the list, reads nothing
and returns nothing, and `steps` gives the list. A consumer, the function of a calculation,
`write_vars`, or the method `iter_blocks`, returns something, and it runs
the recipe: it makes as many passes as its
algorithm needs, each one a reading of the source from its start through
readers of its own, built from the steps the `Variants` has when the pass
starts. So a `Variants` can be given to any number of consumers, a step
can be added between two of them, and a consumer never changes it. The
owner decided this on 21 September 2026, and the options not taken are in
`docs/specs/filters.md`; one more was that a calculation took an iterator
that the user asked the `Variants` for, which a calculation that makes two
passes could not work from.

Every consumer returns a result, and every result has a `pass_stats`:

```python
@dataclass(frozen=True)
class PassStats:
    num_vars: int
    filtering: dict[str, FilteringStats]
```

`num_vars` is how many variants the consumer took, after the steps, and
`filtering` has, for each filter, how many variants it was given and how
many it kept, as `docs/specs/filters.md` says. A calculation that makes
several passes over the same steps gives those of one, since they are the
same in all, and one whose passes differ, a pass for each population, says
in its spec what it gives. In TypeScript it is `passStats`, with `numVars`
and `filtering`.

A `Variants` has `individuals`, a tuple of names,
`num_individuals` and `ploidy`; the first two are pyNei's `samples` and
`num_samples` under the word that `docs/glossary.md` gives, individual. The
other differences from pyNei: the steps are methods that change it, where
pyNei has functions that return another `Variants`; the counts of the
filters are in the results and not in it; it has no
`desired_num_vars_per_chunk`; and it is not iterated over variants. The only way genotypes come out of it is
`iter_blocks`, of `docs/specs/block.md`, which gives them as arrays of a
few thousand variants, for the user who wants them for an analysis of
their own and for the tests. The owner decided this on 20 September 2026;
the option not taken was an iterator of single variants, each copied into
a Python object.

In TypeScript, `Variants` is a class with `individuals`, `numIndividuals`
and `ploidy`, the steps as methods, `steps`, and `iterBlocks`. It lives in the memory of wasm, which the
garbage collector of JavaScript does not see, so it has a `free()` method
that the application calls when it is done with it, and
`[Symbol.dispose]`, which does the same for an application that declares
it with `using`.

### Fields that were not asked for, and sources that lack one

A block holds the columns that were asked for and that its source can
give, and no other: a column that nobody asked for is `None`, the
chromosomes and the positions of a VCF among them. A field that was asked
for may not be there: a source built from an array of genotypes has no
alleles to give. So a consumer that depends on a field looks at the block,
whose `fields()` says which it holds, and fails with an error that names
the field when it is not there.

The genotypes are a field like the others: with `GTS` not asked for, `gts`
is empty, and `num_vars` still says how many variants the block has.

A variant with no id, `.` in a VCF, has an empty `id` in a column that is
there. A variant with no quality has NaN.

A reader can be asked for other fields between two blocks, and the change
holds from the next block it builds. A reader over another reader passes
on what it was asked for, plus what it needs itself.

### How it runs

The chromosome numbers are given in the order in which the names first
appear among the variants that the reader gives, whatever the number of
threads the reader uses, so two passes over the same file give the same
numbers. A reader over another reader has no table of its own and gives
the one of its source.

### How it is verified

There is no number here to check against a reference program. The cargo
tests of this module are of the types themselves: a `Needs` built from
flags contains them and no other, and prints the names of its fields; a
`ChromTable` gives the same number for the same name, numbers in the order
of first appearance, and the name back for a number. The views are
tested at `Block::variants`, in `docs/specs/block.md`. That a reader honours
`Needs` is tested where there is a reader, in `docs/specs/io_vcf.md` and
`docs/specs/io_vars.md`.

## The counts of one variant

### What they give

Two counts over the genotypes of one variant, which the filters and the
statistics work their numbers out from, so that what a missing genotype is
and what a heterozygous one is are written once.

The counts of the genotypes: how many are called, how many are missing and
how many are heterozygous. A genotype is missing when at least one of its
alleles is `MISSING_ALLELE`, so a half called genotype, `0/.` in a VCF, is
missing, and called otherwise. It is heterozygous when it is called and
its alleles are not all the same, at any ploidy. It is what
`_calc_gt_is_missing` and `_calc_gt_is_het` of pyNei's `pynei/gt_counts.py`
compute as masks.

The counts of the alleles: how often each allele was called, and the sum
of them, the called alleles. An allele is counted wherever it was called,
also in a half called genotype. It is what `_count_each_allele` of the
same file computes for a chunk.

Neither has a function in Python or in TypeScript. A user sees them
through the filters and the statistics.

### An allele that no reader gives

An allele below `MISSING_ALLELE`, -2, is in no block that a reader of
popnei gives. Both counts refuse it, with an error that names the value:
counted as it is, it would be a called allele of the genotype counts, and
it has no place among the allele counts. pyNei refuses it too, with a
`ValueError`, in `_count_alleles_per_var`.

What makes the first sentence true is that each reader refuses such an
allele before it builds a block. A genotype of a VCF cannot say one: an
allele number is a run of digits, so `-2/0` and `-1` are wrong data lines,
as `docs/specs/io_vcf.md` has it. A vars file can, since it holds an allele
as a signed byte, and its reader takes the smallest allele of each batch
and refuses the batch that holds one below the missing allele, which
`docs/specs/io_vars.md` has under "What it refuses". The owner decided on
21 September 2026 that such an allele "is never allowed" and is refused
"even by the vcf parser", after a reviewer changed one byte of the
genotypes of a vars file to 254 and got a block holding -2 with no error.

### How it runs

One pass over the alleles of the row each, with no allocation: the allele
counts are written into an array of the caller, one entry for each of the
128 alleles a genotype can hold, which the caller hands over again for the
next variant. The counts of the alleles clear that array themselves before
they count, so what they leave is the counts of the variant they were
given, whatever the array held. A caller that hands over an array it did
not clear would otherwise get two variants added together, and a count
that is already at the largest number its entry holds would wrap in a
release build and say nothing, where the function cannot see that the
caller meant to add: the entries are the counts of one variant, and one
allele of one variant is counted once.

### How it is verified

The cargo tests, made at the two functions, on the six variants of five
diploid individuals of the worked example of `docs/specs/filters.md`. The
counts are those of `_calc_gt_is_het` and `_count_alleles_per_var` of pyNei
at commit ef0ca6e on these genotypes, and on the tetraploid ones below:

| variant | genotypes | called | missing | het | allele counts | called alleles |
|---|---|---|---|---|---|---|
| 1 | 0/0 0/1 0/0 0/0 0/. | 4 | 1 | 1 | 8, 1 | 9 |
| 2 | 0/0 0/1 0/0 ./. 0/. | 3 | 2 | 1 | 6, 1 | 7 |
| 3 | 0/1 2/3 0/1 2/3 ./. | 4 | 1 | 4 | 2, 2, 2, 2 | 8 |
| 4 | ./. ./. ./. ./. ./. | 0 | 5 | 0 | none | 0 |
| 5 | 0/0 0/0 0/0 0/0 1/1 | 5 | 0 | 0 | 8, 2 | 10 |
| 6 | 0/. ./. ./. ./. ./. | 0 | 5 | 0 | 1 | 1 |

Three more: the tetraploid genotypes 0/0/0/1, 1/1/1/1 and 0/./0/0 give 2
called, 1 missing and 1 heterozygous; a ploidy of 0, and genotypes whose
length is not a multiple of the ploidy, are errors; and an allele of -2 is
an error in both.

## A `Variants` from an array of genotypes

Added on 22 September 2026 with the first calculation, the Kosman
distances between individuals of `docs/specs/dists.md`, and left out of
the plan of that calculation the same day by the owner, whose tests run
on small VCF files instead. It has no code and no plan; it waits for a
user who has genotypes in memory.

### What it gives

The variants of an array of genotypes that the user has in memory, as a
`Variants` that takes steps and is given to consumers like one that came
from a file. It is for genotypes that were simulated or that came from
another library, and for the tests, which give the same small array to
pyNei and to popnei.

### Its Python function

```python
Variants.from_gt_array(gts: numpy.ndarray,
                       individuals: Sequence[str]) -> Variants
```

`gts` is an array of integers of variants x individuals x ploidy, an
allele in each cell, 0 to 127, and -1 for a missing one. `individuals` has
one name for each individual, in the order of the second axis. The
`Variants` has no step when it is built.

It mirrors `Variants.from_gt_array` of `pynei/variants.py`. The
differences:

- The second argument is `individuals` and not `samples`, the word of the
  glossary.
- There is no `vars_info`, the frame with the chromosome, the position and
  the other fields of each variant. The blocks of these variants have the
  genotypes and no other column, and a consumer that needs one fails with
  the error that names the field, as "Fields that were not asked for, and
  sources that lack one" says; `iter_blocks` depends on none, so its
  blocks have `None` for the chromosomes and the positions its default
  `fields` asks for; and `write_vars` writes a file without the columns,
  as `docs/specs/io_vars.md` says. It comes with the first test that
  needs positions, those of the LD.
- A numpy masked array is a `ValueError`. pyNei takes what it masks as
  missing; nothing in popnei gives or takes masked arrays, and a user who
  has one writes `gts.filled(-1)`.
- A value below -1 or above 127 is a `ValueError` that says where it is.
  pyNei's `Genotypes` keeps any integer.
- An array of no individuals is a `ValueError`, which pyNei takes: every
  source of popnei refuses one, as `docs/specs/block.md` says.

As in pyNei, an array that is not of integers or not of three dimensions,
a name that is there twice, and a number of names that is not the size of
the second axis are each a `ValueError`. So is a ploidy, the size of the
third axis, of 0 or above 255, the bounds of `open_vcf`.

The array is copied once, as int8, when the `Variants` is built, so what
the user does to their array afterwards changes nothing. An array of no
variants is accepted, and the spec of each calculation says what it does
with a `Variants` that has none. It is a third source of the binding
crate, beside the VCF and the vars file, and where those name their path,
in the `repr` of the `Variants` and in the errors of a pass over it, this
one says "an array of 7 variants x 3 individuals"; the `repr` has its
steps after that, as for the other two.

In TypeScript it is `Variants.fromGtArray(gts, individuals, ploidy)`, with
`gts` an `Int8Array` of the alleles in the same order. An `Int8Array` has
no shape, so the ploidy is given and the number of variants follows from
the length, and a length that is not a multiple of individuals x ploidy is
an error.

### How it runs

Every pass opens a reader over the copy, which gives blocks of the size it
is asked for, or of `default_num_vars_per_block`, each one a new copy of
its rows, because a reader gives its blocks away, and the chain of
`docs/specs/filters.md` is built over it as over any source. The memory
is the array as int8 and a block.

### How it is verified

A pytest test builds a `Variants` from an array of 7 variants, 3
individuals and ploidy 2 with a missing allele in it, and `iter_blocks`
with `fields=()` and `num_vars_per_block=3` gives blocks of 3, 3 and 1
variants whose genotypes, joined, are the array, with a `pass_stats` of 7
variants and no filter; the same array through pyNei's `from_gt_array` gives chunks
with those genotypes. One test for each `ValueError` above. A cargo test,
at `GtArraySource::reader`: 10 variants asked for in blocks of 4 come in
blocks of 4, 4 and 2 that pass `Block::check` and hold the rows in order.

## The Rust interface

The missing allele and the largest allele a genotype can hold.

```rust
pub const MISSING_ALLELE: i8 = -1;
pub const MAX_ALLELE: i8 = i8::MAX;
```

Which fields of the variants a consumer wants, or a block holds. It is a
set of flags with the usual operations of a set. Whether it is written by
hand or with the `bitflags` crate is the implementer's choice.

```rust
pub struct Needs(u8);
impl Needs {
    pub const GTS: Needs;
    pub const CHROM_POS: Needs;
    pub const ID: Needs;
    pub const ALLELES: Needs;
    pub const QUAL: Needs;
    /// The five above, built from them.
    pub const ALL: Needs;
    pub fn empty() -> Needs;
    pub fn contains(self, fields: Needs) -> bool;
    pub fn is_empty(self) -> bool;
    pub fn union(self, other: Needs) -> Needs;
    pub fn difference(self, other: Needs) -> Needs;
}
```

Its `Display` writes the name of each field of the set between backticks,
`` `gts`, `chrom and pos` ``, and `nothing` for the empty set, so that
the message of a consumer that did not get two fields is read as two.

The names of the chromosomes of one reader, each with its number.

```rust
#[derive(Debug, Default)]
pub struct ChromTable { /* private */ }
impl ChromTable {
    pub fn new() -> ChromTable;
    /// The number of `name`, which is added when it is not there yet.
    pub fn intern(&mut self, name: &str) -> u32;
    pub fn name(&self, number: u32) -> Option<&str>;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
}
```

A table holds at most `u32::MAX` names, which no genome comes near: a name
interned beyond that gets the number `u32::MAX`, is not kept, and `name`
gives `None` for it. `intern` returns a number and not a `Result` because
every variant that is read calls it, and popnei does not panic.

One variant of a block. `Block::variants` and `Block::variant`, of
`docs/specs/block.md`, give it, through a constructor that only the crate
sees, and it reads the `AllelesColumn` of that spec: the two modules use
each other, which is legal inside one crate, and the view stays here
because it is what the row helpers take. Every method but `gts` gives
`None` when the block has no such column.

```rust
#[derive(Debug, Clone, Copy)]
pub struct VariantRef<'a> { /* private */ }
impl<'a> VariantRef<'a> {
    /// num_individuals x ploidy alleles, individual after individual.
    /// Empty when the block was built without the genotypes.
    pub fn gts(&self) -> &'a [i8];
    /// A number of the ChromTable of the reader the block came from.
    pub fn chrom(&self) -> Option<u32>;
    /// 1 based, as in a VCF.
    pub fn pos(&self) -> Option<u64>;
    /// Empty when the variant has none.
    pub fn id(&self) -> Option<&'a str>;
    /// NaN when the variant has none.
    pub fn qual(&self) -> Option<f32>;
    pub fn num_alleles(&self) -> Option<usize>;
    /// Allele 0 is the reference. None for an allele the variant lacks,
    /// as for a block with no alleles.
    pub fn allele(&self, allele: usize) -> Option<&'a str>;
}
```

The counts of one variant. `gts` is the genotypes of a `VariantRef`, or a
row of the genotypes of a block.

```rust
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct GtCounts {
    /// Genotypes with no missing allele.
    pub called: u32,
    /// Genotypes with one missing allele at least, the half called among
    /// them.
    pub missing: u32,
    /// Called genotypes whose alleles are not all the same.
    pub het: u32,
}

/// An error for a ploidy of 0, for genotypes whose length is not a
/// multiple of the ploidy, for an allele below MISSING_ALLELE, and for a
/// variant of more alleles than a count of them holds.
pub fn count_gts(gts: &[i8], ploidy: usize) -> Result<GtCounts>;

/// One entry for each allele from 0 to MAX_ALLELE.
pub type AlleleCounts = [u32; 128];

/// It writes into `counts[a]` how often the allele a is in `gts`, and
/// gives how many it counted, the called alleles. It clears `counts`
/// first, so the caller hands the same array over for every variant and
/// clears nothing. An error for an allele below MISSING_ALLELE and for a
/// variant of more alleles than a count of them holds.
pub fn count_alleles(gts: &[i8], counts: &mut AlleleCounts) -> Result<u32>;
```

The genotypes of an array in memory, and the reader over them. The source
checks the values and the sizes once, when it is built, and `reader` can
be called any number of times, once for each pass. Its errors are new
cases of the error of the crate, each a `ValueError` in Python: a length
that is not a whole number of variants, with the length and individuals x
ploidy; an allele below the missing one, with the variant and the
individual it is at; a ploidy of 0 or above 255; no individuals; and a
name that is there twice, with the name. The ploidy and the individuals
are checked before the length is divided by them.

```rust
pub struct GtArraySource { /* private */ }
impl GtArraySource {
    /// `gts` is variants x individuals x ploidy, C order.
    pub fn new(gts: Vec<i8>, individuals: Vec<String>, ploidy: usize)
        -> Result<GtArraySource>;
    pub fn num_vars(&self) -> usize;
    /// A reader from the first variant. `None` asks for the size of
    /// `default_num_vars_per_block`.
    pub fn reader(&self, num_vars_per_block: Option<usize>) -> Result<GtArrayReader>;
}

/// It implements BlockReader. Its ChromTable is empty and its
/// filtering_stats are none.
pub struct GtArrayReader { /* private */ }
```

The reader shares the array with its source, behind an `Arc`, so that it
is `Send` as `BlockReader` asks and a `Variants` given to two consumers
copies the array once and not once per pass. `set_needs` with the
genotypes not asked for gives blocks with an empty `gts`, as "Fields that
were not asked for, and sources that lack one" says.

The error of the crate. Each module adds its cases to one enum, marked
`non_exhaustive`, and `Result<T>` is `std::result::Result<T, Error>`.
Besides the five of the array of genotypes above, this module adds four
cases. Three are of the counts of one variant: a
ploidy of 0 or genotypes that are not a whole number of genotypes of that
ploidy; an allele below the missing one; and a variant of more alleles
than a count of them holds, which is its own case because the other two
say nothing about a variant whose alleles are too many to count. And a
consumer that did not get a field it depends on. It carries the fields as a `Needs`, the ones that were asked
for and that the block does not hold, which a consumer gets with
`asked_for.difference(block.fields())`, and its message names them: a
consumer that depends on two fields reports both in one error.

The three cases of the counts of one variant are a `RuntimeError` in
Python, by the convention the owner gave on 21 September 2026, where a
`RuntimeError` is a defect of popnei and a `ValueError` a wrong input of a
function. The two counts have no function in Python or in TypeScript, so
no user writes the ploidy or the genotypes they refuse: the ploidy is the
one of the reader that built the block, the block of a reader of popnei
holds a whole number of genotypes of it, no reader gives an allele below
the missing one, which "An allele that no reader gives" above says and each
reader is held to by a test of its own, and a variant of more
than 4295 million alleles is a block that no source holds. A user who gets
one of the three reports it instead of looking at what they wrote. In
TypeScript they are an `Error`, as every error of the core is.

## Open points

None. What this spec decides follows from sections 1 and 2 of the
architecture and from two decisions of the owner of 20 September 2026:
the one written under "What a Python and a TypeScript user see", and the
variants in blocks from the source to the calculation, with no `Variant`
that a reader fills.

## Not in this spec

- The dosages and the masks of one variant: later items of this spec.
  They take a `VariantRef` or its genotypes.
- The counts over the individuals of one population and not over all of
  them: with the first statistic per population, in the `stats` spec.
- `Block`, the `BlockReader` trait that everything that gives variants
  implements, and `reblock`: `docs/specs/block.md`.
- The `vars_info` of pyNei's `Variants.from_gt_array`, the other fields of
  the variants of an array: with the first test that needs positions.
