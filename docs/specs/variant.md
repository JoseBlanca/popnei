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
TypeScript user holds. The row helpers of the same module, the dosages, the
missing and het masks and the allele counts of one variant, are items that
are not written yet; they come with the first filter that needs them.

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
with its options or a vars file, and the filters that were put on it. It
holds no genotypes. A user gets one from `open_vcf`, of
`docs/specs/io_vcf.md`, or from `open_vars`, of `docs/specs/io_vars.md`,
and passes it to the filters, which return another `Variants`, and to the
calculations. It can be passed to any number of them: every pass opens a
new reader on the source. It has `individuals`, a tuple of names,
`num_individuals` and `ploidy`; the first two are pyNei's `samples` and
`num_samples` under the word that `docs/glossary.md` gives, individual. The
other differences from pyNei: it has no `desired_num_vars_per_chunk`, and
it is not iterated over variants. The only way genotypes come out of it is
`iter_blocks`, of `docs/specs/block.md`, which gives them as arrays of a
few thousand variants, for the user who wants them for an analysis of
their own and for the tests. The owner decided this on 20 September 2026;
the option not taken was an iterator of single variants, each copied into
a Python object.

In TypeScript, `Variants` is a class with `individuals`, `numIndividuals`
and `ploidy`, and `iterBlocks`. It lives in the memory of wasm, which the
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

The error of the crate. Each module adds its cases to one enum, marked
`non_exhaustive`, and `Result<T>` is `std::result::Result<T, Error>`.
This module adds one case, a consumer that did not get a field it
depends on. It carries the fields as a `Needs`, the ones that were asked
for and that the block does not hold, which a consumer gets with
`asked_for.difference(block.fields())`, and its message names them: a
consumer that depends on two fields reports both in one error.

## Open points

None. What this spec decides follows from sections 1 and 2 of the
architecture and from two decisions of the owner of 20 September 2026:
the one written under "What a Python and a TypeScript user see", and the
variants in blocks from the source to the calculation, with no `Variant`
that a reader fills.

## Not in this spec

- The row helpers, dosages, masks and allele counts of one variant: later
  items of this spec. They take a `VariantRef` or its genotypes.
- `Block`, the `BlockReader` trait that everything that gives variants
  implements, and `reblock`: `docs/specs/block.md`.
- A `Variants` built from an array of genotypes, pyNei's
  `Variants.from_gt_array`: with the first calculation whose tests need
  it.
