# The block module: blocks of variants, the reader trait and reblock

September 2026. A block is a run of consecutive variants held as arrays,
the genotypes of all of them in one. It is how the variants flow through
popnei, from a source to a calculation: a reader gives blocks, a filter
compacts them, a calculation walks their rows or takes them as a matrix,
and a block is also the only way genotypes leave the core, for the Python
or TypeScript user who asks a `Variants` for them. This spec develops the
row `block` of the table in section 9 of `docs/architecture.md` and
sections 1 and 2 of that document. It depends on `docs/specs/variant.md`,
which has the `Needs` that say which fields are wanted, the `ChromTable` of
the chromosome names and `VariantRef`, the view of one variant of a block.
It covers the block, the trait of everything that gives blocks, `reblock`,
which puts blocks back to a size, the reader that reads one block ahead of
the pass that consumes it, and `iter_blocks` in Python and TypeScript.

There is code, built from the first version of this spec, in which a
collector built blocks by copying the variants that a reader gave one at a
time. The owner dropped the single variant and its collector on 20
September 2026, for the reasons at the end of section 1 of the
architecture. `Block`, `AllelesColumn`, the default size, the names of the
fields and `iter_blocks` stay as they were built.

## The block, its readers and reblock

### What it gives

A block holds the genotypes as one array of `i8`, variant after variant,
and inside a variant individual after individual, `ploidy` alleles each; a
missing allele is `MISSING_ALLELE`, -1. It holds the other fields as one
column each, the chromosome numbers, the positions, the ids, the alleles
and the qualities, and a column is there only when it was asked for and
the source has it, as `docs/specs/variant.md` says.

A reader is anything that gives blocks: the VCF reader, the vars file
reader, a filter over another reader, `reblock`. It gives each block away
and the block is the caller's, as section 2 of the architecture says and
for its reasons. A reader never gives a block with no variants, gives
`None` when it has no more, and `None` again at every call after that.
After an error it gives `None` too, at every call, and a reader over
another reader does not call its source again once that source gave an
error. Each reader keeps this itself, the VCF reader, the vars file
reader, a filter and `reblock`: a reader that went on after an error
would give the variants that follow the wrong one as if nothing had
happened.

A reader over another reader does not trust its source to keep that rule.
`reblock` refuses a block of no variants, with the error that names the
defect, and ends there: a source that gives one has nothing to say about
whether the variants after it follow, and a reader that asked again would
turn a source that always gives one into a loop that nothing but the end
of the process breaks.

A filter of variants decides which rows of a block stay and calls
`retain_vars`, which compacts the genotypes and every column in place, in
the order they had. Nothing is allocated, and the block keeps its
capacity. When no variant of a block stays, the filter does not give that
block: it takes the next one from its source, and goes on until a block
has a variant left or the source has no more.

`reblock` is a reader over a reader that gives the variants of its source
in blocks of one size, the last one aside: it joins the blocks that are
too short and cuts the ones that are too long. It goes where the size
matters, which section 2 of the architecture lists: before the matrix work
when a filter took variants out, before the vars file writer, and at the
end of `iter_blocks`. A block that already has the size, with nothing
waiting from the one before, goes through as it is, with no copy.

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
) -> Blocks
```

`Blocks` is an iterator of `Block` with one property, `pass_stats`, the
`PassStats` of `docs/specs/variant.md`: how many variants it has given in
its blocks so far, and the counts of the filters of its pass so far.

`Variants` is the handle of `docs/specs/variant.md`, and this is its only
method that gives genotypes. It is for the user who wants the genotypes
for an analysis of their own, and for the tests. No calculation of popnei
goes through it. Every call starts a new pass over the source, with a
`reblock` at its end, so the blocks a user gets have the size that was
asked for, the last one aside, also when a filter took variants out on the
way. The memory in use is two blocks at most.

`fields` names what each block carries besides the genotypes, among
`"chrom"`, `"pos"`, `"id"`, `"alleles"` and `"qual"`. The chromosome and
the position travel together in the core, so asking for one fills both.
Another name is a `ValueError`. `fields` is a sequence of names, and one
name written where the sequence goes, `fields="alleles"`, is a
`TypeError`: a string is a sequence of its letters, and popnei would
otherwise look for a field called `a`. `num_vars_per_block` is the number
of variants of a block, and `None` is the rule above.

`Block` is a frozen dataclass, declared with `eq=False` and with a
`__repr__` of its own: an array is neither equal nor unequal to another,
so the `==` that a dataclass writes would raise, and the `repr` it writes
prints every genotype, 300 KB for a block of `many.vcf`, into any
traceback. `gts` is a numpy int8 array of variants x
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
which pyNei takes here too and also keeps as a property of the `Variants`,
which popnei's does not have; by
default a block has the chromosomes and the positions and not the ids,
the alleles and the qualities, which a chunk that pyNei read from a VCF
always has; a block has no pandas frame, its columns are tuples and arrays,
which the binding crate hands out and the Python package puts in the
dataclass, so that no result of popnei is a class of the binding crate;
it has no `Genotypes` object with its `to_012` and its masks, which in
popnei are work of the core; and the fields other than the genotypes are
asked for.

In TypeScript, `variants.iterBlocks({fields = ["chrom", "pos"],
numVarsPerBlock})`
is used in `for (const block of variants.iterBlocks())`, what it returns
has a `passStats`, and a block is a
plain object with `gts` an `Int8Array` of variants x individuals x ploidy
in that order, `numVars`, `numIndividuals`, `ploidy`, `chrom`, `id` and
`alleles` arrays, `pos` a
`Float64Array` and `qual` a `Float32Array`,
with `null` where Python has `None`. `numIndividuals` and `ploidy` are
there because `gts` is one flat array: the alleles of the individual `i`
of the variant `v` are the `ploidy` numbers that start at
`(v * numIndividuals + i) * ploidy`, and without the two a block cannot
be read without the `Variants` it came from. In Python the shape of the
array carries them and the dataclass does not hold them again. The arrays
are copies out of the
memory of wasm. The positions are float64 where Python has uint64:
JavaScript's array of unsigned 64 bit numbers gives a `BigInt` for each,
which does not mix with its ordinary numbers in arithmetic, and a
position is exact in a float64 up to 2^53. A position above 2^53 is an
`Error`, which the binding crate finds once per block: rounding it would
give a TypeScript user another position than a Python user, who gets the
uint64 the source has.

`fields` takes an array of names, and `numVarsPerBlock` a whole number of
1 or more and at most 4294967295, which is what a whole number of the core
holds in wasm. A name that is not one of the five, a `fields` that is not
an array, and a `numVarsPerBlock` with a fraction, below 1 or above that
number are an `Error` at the call of `iterBlocks`, which names the value
that was given. A number of JavaScript is a float64 and reaches the core
as an integer of 32 bits, so a size of 2^32 + 1 would otherwise be read as
a block of one variant.

### What a reader of the rules would not guess

The size of the blocks changes nothing but where the cuts fall: the
blocks of a source, joined, are the same for any `num_vars_per_block`.

A block is cut by the count of variants alone. A chromosome that ends in
the middle of a block does not end the block.

The chromosome numbers of a block are those of the table of the reader it
came from, which grows while the source is read, so the name of a number
is looked up after the block was given and not before.

An error loses the block it happened in. With a VCF that has a wrong line
after 250 good ones, read in blocks of 100, the caller gets two blocks and
then the error, and every call after it gives no block. `reblock` loses
also what it was keeping for its next block: over that same reader, with
blocks of 7, it gives the 28 blocks that the 200 variants fill, and then
the error, and the 4 variants that were left over are not given.

`reblock` joins blocks only when they have the same columns. When the
columns of its source change, which a change of `Needs` in the middle of a
pass does, it gives what it was keeping as a shorter block first.

### How it runs

`reblock` keeps at most one block from one call to the next, the variants
that did not fill a block or the ones left after a cut, so its memory is
two blocks. Joining copies the rows of the block that arrives after the
ones that were waiting. Cutting copies out the rows that leave, into a
block allocated for them, and the rest stays in the block that waits, with
the row it starts at: so a block of 10000 variants cut into blocks of 100
copies each row once and not once for every cut before it, and a block
that is given holds the memory of its own rows and not of the block it was
cut from, which is what a Python user keeps when they hold its genotypes.
One memcpy per block given and none per variant. The thread that reads one
block ahead of the pass is the item after this one, and a `reblock` goes
under it, not over it: the chain of readers the thread is given is the one
that cuts and joins the blocks.

Two things the rules above leave to this item. `retain_vars` runs `check`
before it moves a row, since it moves the rows by their place in the
arrays, and gives its error for a block that does not pass it, leaving the
block as it was. And `variants` of a block that does not pass `check`
stops at the first variant that is not in the arrays instead of failing,
so a consumer that did not get its block from `reblock` or from a binding
crate calls `check` before it walks the views.

### How it is verified

There is no number here for a reference program: a block holds what its
reader gave, and what the VCF reader gives is checked against bcftools in
`docs/specs/io_vcf.md`. The checks are that cutting, joining and
compacting lose and change nothing, on the reference VCFs of that spec,
in `tests/reference/vcf/`.

The cargo tests of `reblock`, made at its `next_block` over a `VcfReader`
on `many.vcf` that gives blocks of 100: with every variant given and
blocks of 7 there are 72 blocks, 71 of 7 and one of 3; with blocks of 1,
500; with blocks of 1000, one of 500; and with blocks of 100, five of
100. The genotypes, the positions and the chromosome numbers of the
blocks, joined, are the same for the four sizes. Over a reader written in
the test, which gives blocks built by hand from the four rows of the
`cases.vcf` table of `docs/specs/io_vcf.md`, one variant in each, and
keeps the address of the genotypes of every block it gave: blocks of 3
give one of 3 and one of 1 that hold the four rows, every column of them;
and blocks of 1 give the four blocks of the reader themselves, their
genotypes at the addresses the reader kept, which is how the test sees
that a block of the right size goes through with no copy. A reader written
in the test that gives an error once and would give a block at its next
call is called by `reblock` once and no more.
A VCF written in the test, with a tetraploid genotype in its line 251 of
variants and read in blocks of 100, gives through a `reblock` of 7 the 28
blocks and then the error.

The cargo tests of the block, made at `retain_vars`, at `variants` and at
`check`, on a block built by hand from the four rows of that table with
every column: its views give, each, the genotypes, the position, the id,
the alleles and the quality of its row, `None` for the fields of a column
that is taken out of the block, and there is no view for a fifth variant;
a block whose `gts` lost its last allele fails `check`; keeping
the first, the third and the fourth leaves a block of 3 variants whose
views give the rows 1, 3 and 4 of that table, the alleles and the ids
among them; keeping none leaves a block with `num_vars` 0 and empty
columns; and a `keep` of three values for the four variants is an error
that leaves the block as it was.
`default_num_vars_per_block` gives 10000 variants for 50 individuals, 5000
for 1000 and 100 for 100000.

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
`cases.vcf` in `docs/specs/io_vcf.md`. Three more: a block read with the
default `fields` has the chromosomes and the positions and its three
other columns are `null`; a block kept while enough more are read for the
memory of wasm to grow, which the test checks it did, still holds what it
held, which a column that was a view into that memory would not; and the
arguments the package refuses, a `source` that is not a `Uint8Array`, a
ploidy and a `numVarsPerBlock` that are not whole numbers of 1 or more,
an `onlyPassed` that is not a boolean, and a `fields` that is one name
and not an array of names.

## The reader one block ahead

### What it gives

A pass over the variants asks its reader for a block and then works on it,
and then asks for the next one: on one thread the read of a block and the
work on the block before it never overlap, although the read is the disc and
the decompression and the work is the arithmetic. `with_one_block_ahead`
lends the chain of readers to a thread of its own for as long as one pass
runs. That thread builds the next block while the pass works on the one it
holds, and the pass reads its blocks off a **handle**, which is a reader like
any other, so a pass takes it by wrapping the loop it already has.

It is one block ahead and not more. The handover is a rendezvous: the
reading thread holds at most one block that is built and not yet given, so
the memory of a pass grows by one block and by nothing else. A block holds
`GENOTYPES_PER_BLOCK` genotypes, 5 million, wherever the size popnei chooses
for the individuals of the source is neither held down nor held up, which is
5000 variants of 1000 individuals and 500 of 10000: at the ploidy 2 the
block is 10.0 MB of alleles in both. It is 1.0 MB for 50 individuals, where
the size is held down to `MAX_NUM_VARS_PER_BLOCK`, 10000 variants, and 20.0
MB for 100000 individuals, where it is held up to
`MIN_NUM_VARS_PER_BLOCK`, 100.

What the pass sees does not change, which is what lets a pass take it
without changing a number:

- The blocks arrive in the order the chain gave them, holding the variants
  the chain put in them.
- An error of the chain reaches the pass through `next_block` after every
  block that came before it, and nothing follows it, which is what the
  reader trait promises of any reader.
- The names of the chromosomes travel with each block, so a pass that reads
  the table of names while it walks a block has every name the variants of
  that block added. The counts of the filters travel the same way, as of the
  last block the handle gave.
- The chain is lent and not given away. Whoever built it keeps it, and when
  `with_one_block_ahead` returns the thread is over and the chain itself can
  be asked for its names and for the counts of its filters, which is what a
  pass does when it ends, as `docs/specs/filters.md` has it.
- A change of which fields the pass asks for crosses to the chain and
  reaches it before its next block, so it takes effect one block later than
  it would on one thread. Section 1 of `docs/architecture.md` already
  promises no more than that.

In wasm there is no thread: the pass is given the chain itself and reads the
blocks one after another, as everything else of popnei does there.

### Which passes take it

A pass takes it where a measurement of that pass shows its wall time fall,
and not otherwise, because each one costs a thread that a browser does not
have and one more block of memory. The owner set the rule on 25 September
2026: it is kept in a pass where the best of 5 runs falls by more than 5 per
cent of that pass's wall time and by more than the difference between the
best and the worst of those 5 runs.

What says beforehand whether a pass is worth the thread is two clocks inside
its loop, behind the cargo feature `bench-phases` that
`crates/popnei/Cargo.toml` describes: how long the pass spends inside
`next_block` of its chain, and how long it spends working on the blocks the
chain gave it. A sampling profile gives neither, because the chain
decompresses on the thread that then computes, so the samples of the read
and the samples of the work are of one thread and one stack. The most the
reading thread can take off a pass is the smaller of the two clocks, since
what it hides behind one is the other.

Nine passes over the blocks were measured against the rule on 100000
variants of 1000 individuals with every genotype called, the owner's Apple
M5 Pro of 18 cores, and `docs/reports/perf-read-ahead-2026-09-25.md` has the
two clocks and the wall times of each. Eight of them keep it: the
association study of
`docs/specs/gwas.md`, the kinship, the two passes of the principal
components of the variants, the Kosman distance of every pair of
individuals, the distances between populations, the diversity of every
population, and the two passes of the stats module. The r² matrix of
`docs/specs/ld.md` does not: it reads for 0.006 s and works on the blocks
for 0.004 s of a run of 0.452 s, because the matrix it then computes is
quadratic in the variants it was given, so the most the thread could hide is
1 per cent of that run.

The Kosman distance is the one pass that gave something up for it. It
dropped the sets of bits of a block and the block itself before asking for
the next one, so that the memory of two blocks was never held at once, and a
reading thread that waits with a block that is built and not yet given is
that second block back. It still drops the sets of bits before the next
block, which are as large again as the block.

### How it runs

One thread, spawned inside a `std::thread::scope` so that the chain can be
lent to it by reference and is back with its owner when the scope ends. One
channel of capacity zero carries what the thread read, a block with the
names and the counts the chain then had, or the word that there are no more,
or the error the chain failed with, and nothing follows either of the last
two. A second channel, from the pass to the thread, carries a change of the
fields the pass asks for.

A pass that returns in the middle, which a calculation whose block fails
does and which the Ctrl-C of a Python user comes out as, drops the handle;
the thread is then waiting to hand over a block nobody will ask for, its
send fails, and the thread ends with that block. Nothing joins it by hand
and no thread is left reading: the scope joins it after the pass returns.

### How it is verified

There is no reference program: the reader gives the blocks its chain gives
and changes no number. The cargo tests, over the VCFs of
`docs/specs/io_vcf.md`:

- `many.vcf` read in blocks of 7 through `reblock`, once on one thread and
  once one block ahead, gives the same 475 variants in the same order, every
  field of every one of them equal.
- A VCF written in the test whose line 251 of variants has a tetraploid
  genotype, read in blocks of 100, gives two blocks, then that error, and
  then nothing.
- A reader written in the test that reports the counts of two filters is
  asked for the fields the pass wants and read to its end; the handle
  answers with the counts of the chain, the one chromosome of `cases.vcf`,
  the three individuals and the ploidy 2 while the pass runs, and the chain
  itself answers with them after `with_one_block_ahead` returns.
- A pass that reads one block of four and then returns an error gets that
  error back, leaves the reader with blocks it never gave, and the test ends,
  which a thread still waiting with a block would not let it do.

That the passes give what they gave is checked over a whole panel, at one
thread and at eighteen, by `crates/popnei/benches/the_numbers_of_every_pass.py`,
which writes every number of every pass at seventeen digits so that two
commits can be compared byte for byte.

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

The block, as section 2 of the architecture has it. A column is `None`
when it was not asked for or the source lacks it. An id that is empty and
a quality that is NaN are a variant that has none. The fields are public
because every reader builds blocks: `gts` has `num_vars` x
`num_individuals` x `ploidy` alleles, or none when the genotypes were not
asked for, and every column that is there has `num_vars` entries. Public
fields let a reader with a defect build a block that breaks that, so
`check` says whether a block keeps it, and it is called where a wrong
block would be read wrong with no sign: by both binding crates before the
genotypes of a block cross to numpy or to an `Int8Array`, where they
travel flat, by `reblock` on every block it takes, and by the vars file
writer.

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
impl Block {
    /// Which fields the block holds: the columns that are there, and GTS
    /// when `gts` is not empty or the block has no variants.
    pub fn fields(&self) -> Needs;
    /// The views of its variants, in order.
    pub fn variants(&self) -> impl Iterator<Item = VariantRef<'_>>;
    /// None when the block has no variant `i`.
    pub fn variant(&self, i: usize) -> Option<VariantRef<'_>>;
    /// It keeps the variants whose `keep` is true, in their order, in the
    /// genotypes and in every column, in place, and sets `num_vars` to how
    /// many stayed. A `gts` that is empty, of a block built without the
    /// genotypes, stays empty. `keep` has one value for each variant of
    /// the block; when it has not, that is an error and the block is as
    /// it was.
    pub fn retain_vars(&mut self, keep: &[bool]) -> Result<()>;
    /// That `gts` holds `num_vars` x `num_individuals` x `ploidy` alleles,
    /// or none, and every column that is there `num_vars` entries.
    pub fn check(&self) -> Result<()>;
}
```

The per variant work that uses the threads runs rayon over the rows of
`gts`, `par_chunks(num_individuals * ploidy)`, since the genotypes are all
most calculations read; `variants` is for the work that reads the other
fields too. A row has one allele at least: every source refuses a ploidy
of 0 and a file with no individuals.

The trait of everything that gives blocks, with the contract of "What it
gives". It can be used as a boxed trait object, `Box<dyn BlockReader>`:
neither a pyo3 class nor a wasm-bindgen class can be generic, so both
binding crates hold their reader that way, the trait has no generic method
and no method that takes or returns `Self`, and it is implemented for
`Box<dyn BlockReader>` too, so that what is generic over a reader, a
filter or `reblock`, takes a boxed one. It is implemented for `&mut R`
as well, so that a consumer can be given a reader it does not own:
`write_vars` of `docs/specs/io_vars.md` is given one that way, and the
chain of readers stays with the caller, which reads the counts of the
filters of the pass from it when the call returns, as "How it runs" of the
counts of `docs/specs/filters.md` asks. It asks for `Send`, because the
read ahead thread of section 3 of the architecture moves a reader into
another thread.

```rust
pub trait BlockReader: Send {
    /// The next block, which has one variant at least, or None when there
    /// are no more.
    fn next_block(&mut self) -> Result<Option<Block>>;
    fn individuals(&self) -> &[String];
    fn ploidy(&self) -> usize;
    fn chroms(&self) -> &ChromTable;
    /// ALL until it is called. It holds from the next block that is built.
    fn set_needs(&mut self, needs: Needs);
    /// The kind and the counts of every filter between this reader and
    /// its source, this one first when it is a filter:
    /// `docs/specs/filters.md`. A source gives none.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)>;
}
```

`filtering_stats` has no default, so that a reader over another reader
that forgets to pass on the counts of its source does not compile.

`reblock`.

```rust
pub struct Reblock<R: BlockReader> { /* private */ }

impl<R: BlockReader> Reblock<R> {
    /// `num_vars_per_block` is 1 or more, or None for the default size for
    /// the individuals of `reader`.
    pub fn new(reader: R, num_vars_per_block: Option<usize>) -> Result<Self>;
}

impl<R: BlockReader> BlockReader for Reblock<R> { /* ... */ }

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

The reader one block ahead. `reader` is lent for as long as `body` runs and
is back with its caller when this returns, so the counts of the filters and
the names of the chromosomes are read from it after the pass. `body` is
given a reader like any other, which is the handle in the native build and
`reader` itself in wasm, where there is no thread.

```rust
/// # Errors
///
/// What `body` fails with, and what `reader` fails with, which reaches
/// `body` through `next_block` after the blocks that came before it.
///
/// # Panics
///
/// When the reading thread panics, which is a defect of a reader: the panic
/// is raised again here, where a caller of popnei sees it.
pub fn with_one_block_ahead<R: BlockReader, T>(
    reader: &mut R,
    body: impl FnOnce(&mut dyn BlockReader) -> Result<T>,
) -> Result<T>;
```

A pass that is generic over a reader with no size, which four of the
calculations are, lends `&mut &mut R`: `&mut R` is a reader of its own
whatever `R` is, and it is what travels to the thread.

The name of each column in Python and in TypeScript, and the fields that
those names ask for. Both binding crates take the names from their user
and call this, so that one list serves the two languages and a column
added later cannot reach one of them and not the other.

```rust
/// In the order of the columns of `Block`. The genotypes are not among
/// them: every block a user gets holds them.
pub const FIELD_NAMES: [&str; 5] = ["chrom", "pos", "id", "alleles", "qual"];

/// `GTS` and what the names ask for; the chromosome and the position
/// travel together, so either name asks for both. A name that is not one
/// of `FIELD_NAMES` is the error that names it and lists the five.
pub fn needs_of_the_fields<'a>(
    names: impl IntoIterator<Item = &'a str>,
) -> Result<Needs>;
```

This module adds seven cases to the error of the crate. A
`num_vars_per_block` of 0, which `Reblock::new` and every source that
takes a size refuse. A block the machine cannot give the memory for: a
reader that is given a size refuses one whose genotypes,
`num_vars_per_block` times `num_individuals` times `ploidy`, are more than
a `usize` holds, which in wasm, where a `usize` is 32 bits and holds 4295
million, is 10000 variants of 250000 individuals of the ploidy 2; and a
reader asks for the memory of every column of a block with `try_reserve`
before it fills it, the positions of a variant among them, which are 8
bytes whatever the individuals are, and gives the same error when the
machine does not give it. A size that a caller wrote reaches neither an
abort nor a panic. The error says which of the three sizes of a block it is
about, because what a caller does about it differs: the one the caller
asked for, and they ask for fewer variants in a block; the one popnei chose
for the individuals of the source, and they pass a `num_vars_per_block` at
all; or the one a file fixed, which is the size of a batch of the vars file
of `docs/specs/io_vars.md`, whose reader builds each batch whole whatever
size the caller asked its blocks to be, and they write that file again with
a smaller `num_vars_per_block`. Blocks that do not fit together, which `reblock` finds
when a block of its source has another number of individuals or another
ploidy than the source says it has. A block of no variants, which
`reblock` refuses. A block whose arrays are not of its size, which `check`
finds, with the array and the two sizes. A `keep` that has not one value
for each variant of its block. And a name that is not a field of a block.

## Open points

None. The owner decided on 20 September 2026 that Python and TypeScript
get the genotypes through `iter_blocks(fields=...)` and through nothing
else; the options not taken were an iterator of single variants and one
function that returns the whole matrix. The same day he decided that the
variants flow in blocks from the source to the calculation; the option not
taken, which the first version of this spec had, was a collector that
built the blocks from single variants.

## Not in this spec

- Taking individuals out of a block, which the filter of individuals
  needs: with that filter, in `docs/specs/filters.md`.
- Giving a block back to its reader to be filled again, which section 2
  of the architecture leaves until a measurement asks for it.
- The dosages, the masks and the counts of a block: the row helpers of
  `docs/specs/variant.md` and the calculations that use them.
- Asking for the blocks of one chromosome or of a region: the vars file
  keeps what that needs, `docs/specs/io_vars.md`, and the function is not
  written.
