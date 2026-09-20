# The variant module: the record and the reader trait

September 2026. The `variant` module gives popnei the one thing every
other module passes around: a variant, one site of the genome with the
genotype of every individual at it, held in a struct that the caller owns and
a reader fills, and the trait that anything that gives variants
implements. There is no code. This spec develops the row `variant` of the
table in section 9 of `docs/architecture.md` and section 1 of that
document, which has the reasons for a record that is lent and refilled.
It depends on no other spec. It covers the record, the trait, and the
`Variants` handle that a Python or a TypeScript user holds. The row helpers of the same
module, the dosages, the missing and het masks and the allele counts of
one variant, are items that are not written yet; they come with the first
filter that needs them.

## The record and the reader trait

### What it gives

A reader is anything that gives variants one at a time: the VCF reader,
the vars file reader, a filter over another reader. A consumer, a
calculation or a writer, owns one `Variant`, lends it to the reader again
and again, and the reader fills it with the next variant and says whether
there was one. The buffers inside the `Variant` are allocated once and
reused, so a million variants cost no allocation after the first few.

A `Variant` holds the chromosome as a number, an index into the table of
chromosome names that the reader keeps, the position, the genotypes, the
id, the alleles and the quality. The genotypes are one `i8` per allele,
individual after individual, `ploidy` of them for each: the alleles of
individual i are `gts[i * ploidy .. (i + 1) * ploidy]`. 0 is the reference allele, 1
and above the alternative ones in the order of the VCF, up to 127, and
`MISSING_ALLELE`, -1, is an allele that was not called. There is no mask.

A consumer says which of those fields it wants with a `Needs`, a set of
five flags, and the reader may skip the rest: most calculations want the
genotypes alone, and a vars file reader then never decompresses the
alleles. After each read the `Variant` says in `filled` which fields it
really holds.

### What a Python and a TypeScript user see

Neither sees a variant. pyNei's `Variants` yields chunks, a few thousand
variants as arrays, and its calculations are written over them in Python.
In popnei the loop of every filter and every calculation runs inside the
core, and what a user holds is a handle.

`Variants`, in Python, is that handle: a source of variants, a VCF path
with its options, later a vars file, and the filters that were put on it.
It holds no genotypes. A user gets one from `open_vcf`, of
`docs/specs/io_vcf.md`, and passes it to the filters, which return
another `Variants`, and to the calculations. It can be passed to any
number of them: every pass opens a new reader on the source. It has
`individuals`, a tuple of names, `num_individuals` and `ploidy`; the
first two are pyNei's `samples` and `num_samples` under the word that
`docs/glossary.md` gives, individual. The other differences from pyNei:
it has no `desired_num_vars_per_chunk`, and it is not iterated over
variants. The only way genotypes come out of it is `iter_blocks`, of
`docs/specs/block.md`, which gives them as arrays of a few thousand
variants, for the user who wants them for an analysis of their own and
for the tests. The owner decided this on 20 September 2026; the option
not taken was an iterator of single variants, each copied into a Python
object.

In TypeScript, `Variants` is a class with `individuals`, `numIndividuals`
and `ploidy`, and `iterBlocks`. It lives in the memory of wasm, which the
garbage collector of JavaScript does not see, so it has a `free()` method
that the application calls when it is done with it.

### Fields that were not asked for, and sources that lack one

A field that the consumer did not ask for may still be filled: the VCF
reader always fills the chromosome and the position, which cost nothing
next to the genotypes. A field that was asked for may not be: a source
built from an array of genotypes has no alleles to give. So a consumer
that depends on a field checks `filled`, and fails with an error that
names the field when it is not there.

A field that is not in `filled` holds its empty value, an empty `gts`, an
empty `id`, no alleles, `None` for the quality, and never what the
previous variant left in it.

A variant with no id, `.` in a VCF, has an empty `id` with the `ID` flag
in `filled`.

A reader can be asked for other fields between two reads, and the change
holds from the next read on. A reader over another reader passes on what
it was asked for, plus what it needs itself.

### How it runs

`VariantReader` can be used as a boxed trait object,
`Box<dyn VariantReader>`. Neither a pyo3 class nor a wasm-bindgen class
can be generic, so both binding crates hold their reader that way, and the
trait has no generic method and no method that takes or returns `Self`.
`VariantReader` is implemented for `Box<dyn VariantReader>` too, so that
what is generic over a reader, the collector of `docs/specs/block.md`,
takes a boxed one.
It asks for `Send`, because the read ahead thread of section 3 of the
architecture moves a reader into another thread.

The chromosome numbers are given in the order in which the names first
appear among the variants that the reader gives, whatever the number of threads the reader uses, so
two passes over the same file give the same numbers.

### How it is verified

There is no number here to check against a reference program. The cargo
tests of this module are of the three types themselves: a `Needs` built
from flags contains them and no other; a `ChromTable` gives the same
number for the same name, numbers in the order of first appearance, and
the name back for a number; a `Variant` that was filled and then cleared
holds the empty values above and keeps the capacity of its buffers. That
a reader honours `Needs` and `filled`, and that nothing is allocated from
one variant to the next, is tested where there is a reader, in
`docs/specs/io_vcf.md`.

## The Rust interface

The missing allele and the largest allele a genotype can hold.

```rust
pub const MISSING_ALLELE: i8 = -1;
pub const MAX_ALLELE: i8 = i8::MAX;
```

Which fields of a `Variant` a consumer wants, or a `Variant` holds. It is
a set of flags with the usual operations of a set, union, `contains`,
`ALL` and `empty()`. Whether it is written by hand or with the `bitflags`
crate is the implementer's choice.

```rust
pub struct Needs(u8);
impl Needs {
    pub const GTS: Needs;
    pub const CHROM_POS: Needs;
    pub const ID: Needs;
    pub const ALLELES: Needs;
    pub const QUAL: Needs;
    pub const ALL: Needs;
}
```

The names of the chromosomes of one reader, each with its number.

```rust
pub struct ChromTable { /* private */ }
impl ChromTable {
    pub fn new() -> ChromTable;
    /// The number of `name`, which is added when it is not there yet.
    pub fn intern(&mut self, name: &str) -> u32;
    pub fn name(&self, id: u32) -> Option<&str>;
    pub fn len(&self) -> usize;
}
```

The record. Its fields are public because every reader writes them and
every consumer reads them.

```rust
pub struct Variant {
    /// A number of the ChromTable of the reader that filled this variant.
    pub chrom: u32,
    /// 1 based, as in a VCF.
    pub pos: u64,
    /// num_individuals x ploidy alleles, individual after individual.
    pub gts: Vec<i8>,
    /// Empty when the source gives no id for the variant.
    pub id: String,
    /// The reference allele first, then the alternative ones.
    pub alleles: Vec<String>,
    pub qual: Option<f32>,
    /// What the reader filled in the last read.
    pub filled: Needs,
}
impl Variant {
    pub fn new() -> Variant;
    /// Every field to its empty value and `filled` to nothing. It keeps
    /// the capacity of the buffers.
    pub fn clear(&mut self);
}
```

The trait. `read_variant` clears `var`, fills it and returns true, or
returns false when the source has no more variants, and after a false
every later call returns false. An error ends the reader: what a call
after an error returns is not defined.

```rust
pub trait VariantReader: Send {
    fn read_variant(&mut self, var: &mut Variant) -> Result<bool>;
    fn individuals(&self) -> &[String];
    fn ploidy(&self) -> usize;
    fn chroms(&self) -> &ChromTable;
    /// ALL until it is called.
    fn set_needs(&mut self, needs: Needs);
}
```

The error of the crate. Each module adds its cases to one enum, marked
`non_exhaustive`, and `Result<T>` is `std::result::Result<T, Error>`.
This module adds one case, a consumer that did not get a field it
depends on, with the name of the field.

## Open points

None. What this spec decides follows from section 1 of the architecture
and from the decision of the owner, of 20 September 2026, that is written
under "What a Python and a TypeScript user see".

## Not in this spec

- The row helpers, dosages, masks and allele counts of one variant: later
  items of this spec.
- `Block`, the variants that the calculations with matrices consume, and
  the collector that builds them from a reader: `docs/specs/block.md`.
- A `Variants` built from an array of genotypes, pyNei's
  `Variants.from_gt_array`: with the first calculation whose tests need
  it.
