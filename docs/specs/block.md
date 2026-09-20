# The block module: blocks of variants and their collector

September 2026. A block is a run of consecutive variants held as arrays,
the genotypes of all of them in one. The calculations that want matrices,
the PCA, the kinship, consume blocks, and a block is also the only way
genotypes leave the core: a Python or a TypeScript user who wants them
asks a `Variants` for its blocks. There is no code. This spec develops
the row `block` of the table in section 9 of `docs/architecture.md` and
section 2 of that document. It depends on `docs/specs/variant.md`, which
has the `Variant`, the `Needs` that say which fields are wanted and the
`VariantReader` trait. It covers the block, the collector that builds
blocks from a reader, and `iter_blocks` in Python and TypeScript. The
other things of the row, a view of one variant of a block, putting blocks
back to their size after a filter took rows out, and blocks that come
from the vars file without a copy, are items that are not written yet.

## The block and its collector

### What it gives

A collector takes any reader and gives blocks: it reads variants one at a
time into a `Variant` of its own and copies each one into the arrays of
the block it is building, until the block has the number of variants that
was asked for or the reader has no more. The last block of a source is
the only one that can be shorter. A reader with no variants gives no
block.

A block holds the genotypes as one array of `i8`, variant after variant,
and inside a variant as the `Variant` has them, individual after
individual, `ploidy` alleles each; a missing allele is `MISSING_ALLELE`,
-1. It holds the other fields as one column each, the chromosome numbers,
the positions, the ids, the alleles and the qualities, and a column is
there only when the collector was asked for it. The collector passes to
its reader what it was asked for, as any consumer does, so a column
nobody wants is never parsed.

Blocks are owned, as section 2 of the architecture says: the collector
gives each one away and starts a new one. At the default size and ploidy
2 that is one allocation of 10 MB of genotypes for each block.

A field that was asked for and that the reader did not fill, which
`filled` of the `Variant` tells, is the error of `docs/specs/variant.md`
that names the field, and not a column that is silently absent.

How many variants a block has, when the caller does not say, follows
pyNei's `calc_num_vars_per_chunk` of `pynei/variants.py`: 5 million
genotypes divided by the number of individuals, and no fewer than 100
variants and no more than 10000. The three numbers are inherited from
pyNei's `config.py`. Nobody has measured them for popnei.

### In Python and in TypeScript

```python
Variants.iter_blocks(
    fields: Iterable[str] = ("chrom", "pos"),
    num_vars_per_block: int | None = None,
) -> Iterator[Block]
```

`Variants` is the handle of `docs/specs/variant.md`, and this is its only
method that gives genotypes. It is for the user who wants the genotypes
for an analysis of their own, and for the tests. No calculation of popnei
goes through it. Every call starts a new pass over the source, and the
memory in use is one block.

`fields` names what each block carries besides the genotypes, among
`"chrom"`, `"pos"`, `"id"`, `"alleles"` and `"qual"`. The chromosome and
the position travel together in the core, so asking for one fills both.
Another name is a `ValueError`. `num_vars_per_block` is the number of
variants of a block, and `None` is the rule above.

`Block` is a frozen dataclass. `gts` is a numpy int8 array of variants x
individuals x ploidy, and the binding crate hands the array of the core to
numpy without copying it. `chrom` is a tuple of names, `pos`
a numpy uint64 array, `id` a tuple of strings with `None` for a variant
that has none, `alleles` a tuple with, for each variant, the tuple of its
alleles, the reference first, and `qual` a numpy float32 array with NaN
for a variant that has none. A field that was not asked for is `None`.
`num_vars` is the number of variants.

The three arrays are read only, and so is every array a user can reach
from them: a block is frozen, and an array that a reshape or a slice of
the core's allocation left writable would be a way around that. The one a
user works on is their own, `numpy.array(block.gts)`.

It mirrors `Variants.iter_vars_chunks` of pyNei and its `VariantsChunk`.
The differences: the names, since `docs/glossary.md` keeps "chunk" for
pyNei's, and so `num_vars_per_block` for `desired_num_vars_per_chunk`,
which is now an argument of this method and not of the `Variants`; by
default a block has the chromosomes and the positions and not the ids,
the alleles and the qualities, which pyNei's chunk always has; a block has no pandas frame, its columns are tuples and arrays,
which the binding crate hands out and the Python package puts in the
dataclass, so that no result of popnei is a class of the binding crate;
it has no `Genotypes` object with its `to_012` and its masks, which in
popnei are work of the core; and the fields other than the genotypes are
asked for.

In TypeScript, `variants.iterBlocks({fields = ["chrom", "pos"],
numVarsPerBlock})`
is used in `for (const block of variants.iterBlocks())`, and a block is a
plain object with `gts` an `Int8Array` of variants x individuals x ploidy
in that order, `numVars`, `chrom`, `id` and `alleles` arrays, `pos` a
`Float64Array` and `qual` a `Float32Array`,
with `null` where Python has `None`. The arrays are copies out of the
memory of wasm. The positions are float64 where Python has uint64:
JavaScript's array of unsigned 64 bit numbers gives a `BigInt` for each,
which does not mix with its ordinary numbers in arithmetic, and a
position is exact in a float64 up to 2^53.

### What a reader of the rules would not guess

The size of the blocks changes nothing but where the cuts fall: the
blocks of a source, joined, are the same for any `num_vars_per_block`.

A block is cut by the count of variants alone. A chromosome that ends in
the middle of a block does not end the block.

The chromosome numbers of a block are those of the reader's table, which
grows while the source is read, so the name of a number is looked up
after the block was collected and not before.

When the reader gives an error, the collector gives that error and not
the block it was building: the variants of that block that were already
read are lost with it. With a VCF that has a wrong line after 250 good
ones and blocks of 100, the caller gets two blocks and then the error,
and every call after it gives no block.

### How it runs

At the block level of the architecture, over any `VariantReader`. What
is kept from one block to the next is the `Variant` that the collector
lends to its reader. The read ahead thread of section 3 of the
architecture, which collects the next block while the consumer works on
the one in hand, is not in this item.

### How it is verified

There is no number here for a reference program: a block holds what the
reader gave, and what the VCF reader gives is checked against bcftools in
`docs/specs/io_vcf.md`. The checks are that collecting loses and changes
nothing, on the reference VCFs of that spec, in `tests/reference/vcf/`.

The cargo tests, made at the collector's `next_block` over a `VcfReader`
on `many.vcf`, but for the default sizes: with blocks of 100 variants and the default options of the
reader, which give 475 of its 500 variants, there are five blocks, of
100, 100, 100, 100 and 75 variants; with every variant given, five of
100; with blocks of 1000, one of 475. The genotypes, positions and
chromosome numbers of the blocks, joined, are those that `read_variant`
gives one by one, for blocks of 1, of 7, of 100 and of 1000 variants.
With the genotypes alone asked for no other column is there, the
chromosomes and the positions neither, although the VCF reader fills them
in every `Variant`. `default_num_vars_per_block` gives 10000 variants for
50 individuals, 5000 for 1000 and 100 for 100000. A reader
with no variants gives no block. A VCF written in the test, with a
tetraploid genotype in its third variant and blocks of 2, gives one
block and then the error.

Against pyNei, a pytest test made at `open_vcf(...).iter_blocks(...)`:
`cases.vcf`, `cases.vcf.gz`, `many.vcf` and `many.vcf.gz` are read with
both libraries, popnei with `only_passed=False` and every field asked
for, and the blocks of popnei, joined, are compared with the chunks of
pyNei, joined: the genotypes, the chromosomes, the positions and the
alleles exactly, the ids with `None` for pyNei's `pandas.NA`, and the
qualities with NaN for it, and `individuals`, `num_individuals` and
`ploidy` with pyNei's `samples`, `num_samples` and `ploidy`. It is run
with `num_vars_per_block` of 7 and with the default. This is also the comparison with pyNei of the VCF reader.

The TypeScript test, under node, reads `cases.vcf` from a `Uint8Array`,
asks for every field with blocks of 3 variants and every variant given,
and compares the two blocks, of 3 variants and of 1, with the table of
`cases.vcf` in `docs/specs/io_vcf.md`.

## The Rust interface

The alleles of the variants of a block. How it keeps the texts is
private, and it is meant to be one buffer and not a string for each
allele.

```rust
pub struct AllelesColumn { /* private */ }
impl AllelesColumn {
    pub fn num_vars(&self) -> usize;
    /// 0 for a variant the column does not hold.
    pub fn num_alleles(&self, var: usize) -> usize;
    /// Allele 0 is the reference. A variant or an allele the column does
    /// not hold gives an empty text, which no allele of a source is.
    pub fn allele(&self, var: usize, allele: usize) -> &str;
}
```

The block, as section 2 of the architecture has it, with the individuals
under the word of the glossary. A column is `None` when it was not asked
for. An id that is empty and a quality that is NaN are a variant that has
none.

```rust
pub struct Block {
    pub num_vars: usize,
    pub num_individuals: usize,
    pub ploidy: usize,
    /// num_vars x num_individuals x ploidy, variant after variant.
    pub gts: Vec<i8>,
    /// Numbers of the ChromTable of the reader the block came from.
    pub chrom: Option<Vec<u32>>,
    pub pos: Option<Vec<u64>>,
    pub id: Option<Vec<String>>,
    pub alleles: Option<AllelesColumn>,
    pub qual: Option<Vec<f32>>,
}
```

The collector. It owns its reader, and gives it back to who wants the
table of chromosomes or the individuals.

```rust
pub struct BlockCollector<R: VariantReader> { /* private */ }

impl<R: VariantReader> BlockCollector<R> {
    /// `needs` is what each block will hold; GTS is always part of it.
    /// `num_vars_per_block` is 1 or more, or None for the default size.
    pub fn new(
        reader: R, needs: Needs, num_vars_per_block: Option<usize>,
    ) -> Result<Self>;
    /// The next block, or None when the reader has no more variants.
    pub fn next_block(&mut self) -> Result<Option<Block>>;
    pub fn reader(&self) -> &R;
}

/// The default number of variants of a block for that many individuals:
/// the genotypes of a block divided by the individuals, and never fewer
/// than the smallest number of variants nor more than the largest.
pub fn default_num_vars_per_block(num_individuals: usize) -> usize;

/// The three numbers of that rule, each inherited from pyNei's
/// `config.py` and measured for popnei by nobody. No individual at all
/// gives the largest, as pyNei's division by `max(num_samples, 1)` does.
pub const GENOTYPES_PER_BLOCK: usize = 5_000_000;
pub const MIN_NUM_VARS_PER_BLOCK: usize = 100;
pub const MAX_NUM_VARS_PER_BLOCK: usize = 10_000;
```

The name of each column in Python and in TypeScript, and the fields a
collector is asked for to fill the columns those names ask for. Both
binding crates take the names from their user and call this, so that one
list serves the two languages and a column added later cannot reach one
of them and not the other.

```rust
/// In the order of the columns of `Block`. The genotypes are not among
/// them: every block holds them.
pub const FIELD_NAMES: [&str; 5] = ["chrom", "pos", "id", "alleles", "qual"];

/// `GTS` and what the names ask for; the chromosome and the position
/// travel together, so either name asks for both. A name that is not one
/// of `FIELD_NAMES` is the error that names it and lists the five.
pub fn needs_of_the_fields<'a>(
    names: impl IntoIterator<Item = &'a str>,
) -> Result<Needs>;
```

`BlockCollector<Box<dyn VariantReader>>` is how the two binding crates
hold it, so `VariantReader` is implemented for a box of itself.

This module adds four cases to the error of the crate. A
`num_vars_per_block` of 0, which `BlockCollector::new` refuses. A block
the machine cannot give the memory for: `new` refuses one whose
genotypes, `num_vars_per_block` times `num_individuals` times `ploidy`,
are more than a `usize` holds, which in wasm, where a `usize` is 32 bits
and holds 4295 million, is 10000 variants of 250000 individuals of the
ploidy 2; and the first `next_block` asks for the memory of every column
of the block, the positions of a variant among them, which are 8 bytes
whatever the individuals are, and gives the same error when the machine
does not give it, before a variant is read. A size that a caller wrote
reaches neither an abort nor a panic. A variant the reader filled with a
number of alleles other than its individuals times its ploidy, which
would be a block whose `gts` is not `num_vars` x `num_individuals` x
`ploidy` and whose genotypes a consumer reads wrong. And a name that is
not a field of a block.

## Open points

None. The owner decided on 20 September 2026 that Python and TypeScript
get the genotypes through `iter_blocks(fields=...)` and through nothing
else; the options not taken were an iterator of single variants and one
function that returns the whole matrix.

## Not in this spec

- A view of one variant of a block and the copy of it back into a
  `Variant`, `reblock`, and the blocks that the vars file reader gives
  without a copy: later items of this spec.
- The read ahead thread: with the first calculation that consumes blocks.
- The dosages, the masks and the counts of a block: the row helpers of
  `docs/specs/variant.md` and the calculations that use them.
- Asking for the blocks of one chromosome or of a region: not planned.
