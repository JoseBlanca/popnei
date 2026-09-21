# The io::vars module: the vars file, its writer and its reader

September 2026. The vars file is where a user of popnei keeps their variants
once the VCF has been read, so that the text is parsed once and every later
pass reads a file of arrays: the genotypes of 20000 variants of 1000
individuals are decompressed from it in 19 ms. There is no code. This spec
develops the row `io::vars` of the table in section 9 of
`docs/architecture.md` and section 6 of that document. It covers the format of
the file, the writer and the reader.

The format is popnei's own. The owner decided on 20 September 2026 that it
owes nothing to the vars file of pyNei, the Python library that popnei
succeeds: pyNei will not be used once popnei exists and no vars file of pyNei
was ever made in production, so neither library reads the files of the other.
The option not taken was the file of pyNei as a contract between the two,
which is what section 6 of the architecture said until that day. What is kept
from pyNei is what worked: an arrow file, one batch for each block, the
genotypes as one flat buffer.

It depends on `docs/specs/block.md`, which has the `Block`, the run of
consecutive variants held as arrays that the variants flow in, the
`BlockReader` trait of everything that gives blocks, `reblock`, which puts
the blocks of a reader back to one size, and `default_num_vars_per_block`,
the number of variants popnei puts in a block; on `docs/specs/variant.md`
for the `Needs` that say which fields a consumer wants and the `ChromTable`
that turns a chromosome name into a number; and on `docs/specs/io_vcf.md`
for the reference VCFs the tests here are built on, `cases.vcf` and
`many.vcf`, and for the tables of what those two hold.

## The file

### What it holds

A vars file is one arrow IPC file, also called feather v2: a header, then the
record batches one after another, then a footer that carries the schema and
the offset of every batch. Arrow is a column format with libraries in most
languages, so a vars file opens in pandas, in R and in polars as a table, with
no popnei installed; popnei itself reads and writes it with arrow-rs. A record
batch is the columns of a run of variants, and popnei writes one batch for
each block, the few thousand consecutive variants of section 2 of the
architecture. The footer is at the end, which is why the reader needs a source
it can seek in and not only read forward.

The columns, in the order in which a `Block` has them:

| column | arrow type | nulls | what it holds |
|---|---|---|---|
| `chrom` | `Utf8` | no | the name of the chromosome, as text, for each variant |
| `pos` | `UInt64` | no | the position, 1 based |
| `id` | `Utf8` | yes | the id, null for a variant that has none |
| `alleles` | `List<Utf8>` | no | one row per variant, the reference allele first |
| `qual` | `Float32` | yes | the quality, null for a variant that has none |
| `gts` | `FixedSizeList<Int8>[num_individuals * ploidy]` | no | the genotypes |

A fixed size list keeps no offsets, so the `gts` column is one flat buffer of
variants x individuals x ploidy signed bytes in that order, which is the
genotype array of a block itself. An allele is 0 for the reference, 1 and
above for the alternative ones and -1 for one that was not called, as
`docs/specs/variant.md` has it; no mask is written, and the -1 is what says
that an allele is missing.

The chromosome is text in every row and not a number into a table kept beside
the columns, so that the file says what it holds to any program that opens it.
The repeated names cost next to nothing once compressed.

A file holds only the columns its source could fill: a source with no alleles,
a `Variants` built from an array of genotypes, gives a file with no `alleles`
column. `gts` is always there.

What is known before the first variant goes in the schema, under the key
`popnei`, whose value is a json object:

| key | value |
|---|---|
| `format_version` | `"1.0"` |
| `individuals` | the names of the individuals, in order |
| `ploidy` | how many alleles a genotype holds |
| `num_vars_per_block` | how many variants a batch holds, the last one aside |

What is known only after the last variant goes in the footer, under the key
`popnei_batches`: arrow-rs writes the metadata of the footer when the file is
finished, while the schema is written before the first batch. Its value is a
json array with one object for each batch, in the order of the batches:

    [{"num_vars": 100, "regions": [{"chrom": "chr1", "min_pos": 8400, "max_pos": 10213},
                                   {"chrom": "chr2", "min_pos": 10250, "max_pos": 12063}]},
     ...]

`num_vars` is how many variants the batch holds, which the footer of an arrow
file does not say, so the number of variants of a file is known without
reading any batch. `regions` has one entry for each chromosome with a variant
in the batch, in the order in which they first appear, with the smallest and
the largest position of its variants there. With them, a function that is
asked for a region can skip, without decompressing them, the batches that
have no entry overlapping it, and that is right for any file, sorted or not,
because the two positions are the smallest and the largest and not those of
the first and the last variant. This spec puts the regions in the file and
gives them to the caller of the Rust reader; the function that asks for a
region is not built here (see "Not in this spec"). The owner asked for the
regions on 20 September 2026. A file without the `chrom` and `pos` columns
has `num_vars` alone.

The version has two parts. A reader refuses a file whose first part is not the
one it knows, and reads a file with any second part, a later one than its own
too, ignoring the keys and the columns it does not know: that is what lets a
later version add a column without making the old files or the old readers
useless.

### The compression

Arrow compresses each buffer of a batch on its own, and allows two
compressions: zstd and lz4. popnei writes lz4 and reads lz4 and files with no
compression, in every one of its builds. The owner decided this on 20
September 2026.

What it costs and what it saves was measured on a panel of 1000 diploid
individuals and 20000 variants, 40.0 MB of genotypes, with 3 in 100 genotypes
missing whole, simulated by pyNei's `test/gwas_reference/make_reference.py`
with its constants `NUM_SAMPLES` and `NUM_VARS` set to 1000 and 20000 and a
seed of 42: `simulate_genotypes` draws the genotypes, and 3 in 100 of them are
then set to missing, as the `main` of that script does. Its variants are
common and unlinked, which compresses worse than a real panel does. One file
per compression, with the `gts` column alone, was written with pyarrow 23.0.0,
and each was read from memory and its genotypes summed, the best of 5 runs on
the owner's Apple M5 Pro, by arrow-rs 60 on one thread and by `load_vars`,
the function of pyNei that reads its vars file, which does it with pyarrow:

| | the file | bits per diploid genotype | arrow-rs | pyNei |
|---|---|---|---|---|
| lz4 | 15.68 MB | 6.27 | 19 ms | 20 ms |
| zstd | 7.11 MB | 2.84 | 26 ms | 37 ms |
| none | 40.02 MB | 16.01 | 1 ms | 6 ms |

So zstd makes the file 2.2 times smaller and takes 7 ms more to read in Rust.
What decided for lz4 is the build. arrow-rs takes lz4 from `lz4_flex`, which
is pure Rust with no build script, and zstd from the `zstd` crate, which wraps
the C library. For `wasm32-unknown-unknown`, the target of the wasm package
that a web application installs, a trial crate that reads and writes arrow
files built with lz4 under plain cargo, and with zstd only with a second
compiler, the clang of the emscripten SDK; what failed and how it was made to
link is in section 5 of `docs/rust_core.md`. The wasm of the trial crate,
release and with no work on its size, is 2.61 MB with no compression, 2.64 MB
with lz4 and 3.15 MB with zstd.

No build of popnei carries the zstd crate, natively or in wasm, so a file
whose buffers are zstd is refused everywhere: arrow-rs says `zstd IPC
decompression requires the zstd feature`, and popnei wraps that in an error
that says the file is compressed with zstd and that popnei reads lz4. Arrow
decompresses a batch when it is read, so the error comes with the first
block and not when the file is opened. No writer of popnei makes such a
file; another arrow program could.

## The writer

### What it gives

It takes the blocks of any reader, the VCF reader or a filter over it, and
writes them into a vars file, one batch for each block, with the key of the
schema and the key of the footer above. A filter leaves blocks of uneven
size, so the blocks go through a `reblock` first and the batches hold
`num_vars_per_block` variants, the last one aside. The columns of the file
are those of the first block, which its `fields()` says, so a source that
gives the genotypes alone gives a file with a `gts` column only. A later
block with other columns is an error naming the field, because every batch
of an arrow file shares one schema.

A block holds an empty id and a quality of NaN for a variant that has none,
and the writer turns both into nulls, which is what any other program that
opens the file takes for a value that is not there: pandas would show an
empty text and a number where there is neither.

The memory it uses is one block: the writer is given each block to keep, and
the vector of its genotypes becomes the flat buffer of the `gts` column with
no copy.

A source with no variants gives a file with the two keys, a `gts` column and
no batch, which reads back as no variants. It follows the decision the owner
took on 20 September 2026 for a VCF with no variants, that an empty source is
not an error.

### Its Python and TypeScript functions

```python
def write_vars(
    variants: Variants, path: str | Path, num_vars_per_block: int | None = None
) -> VarsWritten
```

`VarsWritten` is a frozen dataclass with one field, `pass_stats`, the
`PassStats` of `docs/specs/variant.md`: how many variants were written,
and how many each filter of the `Variants` was given and kept. The owner
decided on 21 September 2026 that every consumer of a `Variants` returns
it. A pytest test made at `write_vars`, once the filters of
`docs/specs/filters.md` are built: on `many.vcf` with every variant given
and the missing data filter at 0.04, the `pass_stats` has a `num_vars` of
215 and a `filtering` of `{"missing_data": FilteringStats(500, 215)}`, and
the file read back has 215 variants.

`Variants` is the handle of `docs/specs/variant.md`, which holds a source and
the filters on it and no genotypes. The call reads the whole source once. A
path that already exists is a `ValueError` and nothing is written, as in
pyNei. When the source fails halfway, on a wrong line of a VCF, the error is
given and the file that was being written is removed. It is decided here:
pyNei leaves the file, 3458 bytes of it in a trial with a wrong line, and
then the path is taken and the same call cannot be tried again.
`num_vars_per_block` is how many variants a batch holds, and `None` is
`default_num_vars_per_block` of `docs/specs/block.md`, so a file read back
with the default size of block gives its batches as they are.

In TypeScript, `writeVars(variants, {numVarsPerBlock})` gives back an
object with `bytes`, a `Uint8Array` with the bytes of the file, and
`passStats`. The page offers the bytes as a download:
a tab has no filesystem, as section 11 of the architecture says.

It has the name and the first two arguments of `write_vars` of
`pynei/io_vars.py`. The differences from pyNei:

- The file is another one. pyNei does not read it.
- `num_vars_per_block` is an argument here, which is decided here: the file
  has to be written at some size, in pyNei that size belongs to the
  `Variants`, and popnei's has none.
- A source with no variants is written. pyNei raises `ValueError("There are no
  variants to write")`.
- A file that a failed call was writing is removed. pyNei leaves it.
- pyNei writes whatever columns the chunks of its source happen to carry;
  popnei asks its source for every field, so a file written from a VCF has all
  six columns, the ids, the alleles and the qualities whether or not the user
  will read them, and it can stand in for the VCF in any later analysis. The
  owner decided this on 20 September 2026. The option not taken was a `fields`
  argument like the one of `iter_blocks` of `docs/specs/block.md`, which would
  let a user keep the genotypes alone, a file smaller and a VCF parsed faster;
  nothing measures yet what the extra columns cost, and the argument can be
  added, with every field as its default, when a user asks for it.

### How it runs

Over any `BlockReader`, into any `Write`: it never seeks, because an arrow
file is written forward and its footer last. What it keeps from one block to
the next is one entry of `popnei_batches` for each batch written, the count
and, from a pass over the chromosomes and the positions of the block, the
smallest and the largest position of each chromosome. It uses no threads.
Writing the chromosome as text needs the name behind each number of the
block, so the writer is given the `ChromTable` of the reader with each block;
the table has the names already, since the block came out of that reader.

### How it is verified

The reference outside the project is pyarrow 23.0.0, the implementation of the
arrow format that Apache Arrow publishes, which opens what popnei wrote as any
other program would, together with bcftools 1.24, whose account of what
`many.vcf` holds is `many.bcftools.tsv` of `docs/specs/io_vcf.md`. A trial
file in this format, written with arrow-rs 60, was opened with
`pyarrow.ipc.open_file`, with `pyarrow.feather.read_table` and with pandas.
The three gave its columns and its nulls. `open_file` gave both keys, the one
of the footer as the `metadata` of the opened file; `read_table` gave the key
of the schema alone, and pandas no key.

A pytest test made at `write_vars`, on `many.vcf` of `tests/reference/vcf/`,
the 500 variants of 50 individuals, read with `only_passed=False` and written
with `num_vars_per_block` 100. pyarrow opens the file. Its schema is the six
columns with the types and the nulls of the table above, in that order. The
value of `popnei` parses as json and holds `format_version` `"1.0"`, the 50
names of the individuals, `ploidy` 2 and `num_vars_per_block` 100. There are
five batches of 100 variants. The `id` column has 167 nulls and the `qual`
column 100, the variants of `many.bcftools.tsv` with a dot in those columns.
The chromosomes, the positions, the ids, the alleles and the genotypes of the
table that pyarrow reads are those of `many.bcftools.tsv`, compared exactly, with the genotypes of bcftools turned
into numbers as the tests of the VCF reader do. The value of `popnei_batches`,
which pyarrow gives as the `metadata` of the opened file, is these five
entries, worked out from `many.bcftools.tsv` by taking its rows 100 at a time:

| batch | `num_vars` | `regions` |
|---|---|---|
| 1 | 100 | chr1 1000 to 4663 |
| 2 | 100 | chr1 4700 to 8363 |
| 3 | 100 | chr1 8400 to 10213, chr2 10250 to 12063 |
| 4 | 100 | chr2 12100 to 15763 |
| 5 | 100 | chr2 15800 to 19463 |

With the default of `only_passed`, which gives 475 of the variants, the five
batches hold 100, 100, 100, 100 and 75, and their regions are chr1 1000 to
4848; chr1 4885 to 8770; chr1 8807 to 10213 and chr2 10250 to 12655; chr2
12692 to 16540; chr2 16577 to 19463. With no `num_vars_per_block` there is one
batch, since 500 variants are fewer than the 10000 that
`default_num_vars_per_block` gives for 50 individuals.

A second pytest test, at `write_vars` on a `Variants` whose source has no
variants: the file exists, pyarrow reads its schema and its two keys, and
there is no batch.

A third, at `write_vars` on a VCF written in the test, with a tetraploid
genotype in its third variant and `ploidy` 2: the call raises a `ValueError`
and there is no file at the path. A fourth: a second call on a path that the
first one wrote raises a `ValueError` and leaves the file as it was.

The first cargo test builds by hand a block of the four variants of the
`cases.vcf` table of `docs/specs/io_vcf.md`, three individuals and the
genotypes given there, writes it into a `Vec<u8>` with `write_vars` and
`num_vars_per_block` 3, from a reader made in the test that gives that one
block, and reads it back with the reader of this spec: there are two blocks,
of 3 variants and of 1, every field of every variant is what went in,
the empty id of the last three among them, and the batches of the footer are
3 variants with chr1 100 to 300 and 1 variant with chr1 400 to 400. A second
cargo test writes a block of four variants that are not sorted, chr1 300,
chr2 50, chr1 100, chr2 60, in one batch: its regions are chr1 100 to 300 and chr2 50 to 60,
in that order. A third writes `many.vcf` through the VCF reader and reads it
back, and the genotypes, the chromosome names and the positions of the blocks
are those of the blocks of the VCF. A fourth writes a source whose blocks carry
the genotypes alone and finds one column in the file it reads back, and
`num_vars` with no `regions` in its footer.

The TypeScript test, under node, reads `cases.vcf` from a `Uint8Array` with
`openVcf` and `onlyPassed` false, so that it gives the four variants, writes
it with `writeVars` and `numVarsPerBlock` 3, reads the bytes back with
`openVars` and compares the blocks with that same table.

## The reader

### What it gives

It reads a vars file and gives each of its batches as a block, as they are
in the file, so a file written with the default size is read back in blocks
of the default size with no `reblock`. The
names of the individuals, the ploidy, the number of variants and the regions
of every batch come from the two keys, so they are known as soon as the file
is opened, before any batch is read. The chromosome names of the batches are
interned into a `ChromTable` as they are read, so a number means the order of
first appearance among the variants that were given, as in the VCF reader.

Only the columns a consumer asks for are decompressed. A `Needs` becomes a
list of column indices that arrow-rs skips the rest of: for a column left out
its `skip_field` walks past the buffers without calling `read_buffer`, which
is what decompresses. How much that saves is small when there are many
individuals, because the genotypes are then most of the file. On a file
written for this measurement, 20000 variants and 1000 diploid individuals with
all six columns and a distinct id and quality for every variant, its genotypes
drawn at random and so not those of the panel above, 20.67 MB with lz4,
arrow-rs read every column and summed the genotypes in 14 ms and the genotypes
alone in 13 ms, the best of 5 runs on the owner's Apple M5 Pro. Those two are
to be compared with each other and not with the 19 ms of the panel above,
whose genotypes are others and decompress at another speed. The
`CHROM_POS` flag of `Needs` asks for both `chrom` and `pos`, which travel
together. A field that is asked for and whose column the file lacks gives a
block without that column, which is the rule of `docs/specs/variant.md` for a
source that has no such field.

### Its Python and TypeScript functions

```python
def open_vars(path: str | Path) -> Variants
```

It gives the same `Variants` handle as `open_vcf`, and reads only the schema
and the footer when it is called, so a file that is not a vars file fails at
the call and not at the first calculation. Every pass over the handle opens
the source again.

In TypeScript, `openVars(source)`, where `source` is a `Uint8Array` or a
`File` that the user picked in the page, which can be read only inside a web
worker.

It does what `load_vars` of `pynei/io_vars.py` does. The differences from
pyNei:

- The name. The owner decided on 20 September 2026 that it is `open_vars`, for
  the same reason as `open_vcf`: the call opens the file and reads no
  variants. The option not taken was pyNei's name, `load_vars`.
- The file is another one. A vars file of pyNei is refused as a file without
  the `popnei` key.
- `desired_num_vars_per_chunk` is gone. How many variants come out at a time
  is `num_vars_per_block` of `iter_blocks`, of `docs/specs/block.md`, and it
  does not have to be the batch size of the file.

### What it refuses

The reader takes the columns with the types of the table of "What it holds"
and no others: a column of that table with another type is an error that names
the column, the type found and the one expected. A column it does not know is
ignored, as the rule of the versions asks.

The width of `gts` has to be the number of `individuals` times the `ploidy` of
the `popnei` key; if not, the error gives the width found and the one
expected, because the width is what turns the flat buffer into variants.

A null in `chrom`, `pos`, `alleles` or `gts` is an error naming the column and
the variant. A null `id` is the empty id and a null `qual` is no quality.

These are errors of the file as a whole, found when it is opened: it is not an
arrow IPC file; its schema has no `popnei` key, or the value is not json, or
one of its four keys is missing; the first part of `format_version` is not
`1`, which the message gives along with the version found; it has no
`popnei_batches`, or its entries are not as many as the batches. A path that
is a directory is an error of `from_path`, which looks at the path itself:
opening a directory succeeds on macOS and only the first read fails.

A batch that does not hold the `num_vars` its entry of the footer gives is an
error when that batch is read, so that the number of variants that the file
announces is never a wrong one that goes unnoticed.

Two individuals with the same name are an error, as they are for the VCF
reader.

### How it runs

As a `BlockReader`, over any `Read + Seek`: natively a file, in a tab the
bytes of a `Uint8Array` or a `File` read through `FileReaderSync`, the call
that reads a range of bytes of a picked file and returns when it has them,
which exists only inside a web worker. `Read + Seek` is what section 1 of the
architecture asks the vars reader to be generic over, for exactly this. The
footer is at the end of the file, so the reader seeks there when it is
opened, and then to each batch in turn.

A batch becomes a block column by column. The genotypes are copied from the
buffer that arrow-rs decompressed into the vector of the block, which a
`Block` owns; measured on the panel of "The compression", that copy into a
vector allocated for each batch added 0.4 ms to the 18.8 ms of the
decompression of its four batches. The positions and the qualities are
copied too, a null quality as NaN, the ids and the alleles go into the
columns of texts of the block, and each chromosome name gets its number from
the table, which is looked up only when the name differs from that of the
variant before. What the reader keeps from one block to the next is the table
of chromosome names, one entry per name. It uses no threads: reading one
batch ahead in a thread of its own comes with the read ahead of
`docs/specs/block.md`.

Asking for other fields holds from the next block, as `docs/specs/block.md`
asks of every reader, and the projection of each batch is chosen when it is
read.

### How it is verified

A pytest test made at `open_vars`, on the file that `write_vars` makes from
`many.vcf` with `only_passed=False`: the blocks of
`open_vars(path).iter_blocks(fields=...)` with every field, joined, are equal
to the blocks of `open_vcf("many.vcf", only_passed=False).iter_blocks(...)`,
field by field, and hold the counts of `docs/specs/io_vcf.md` as literals: 500
variants, 250 of them in `chr2`, 1511 missing genotypes, 2765 missing alleles,
47235 called alleles and 25954 for the sum of their allele numbers.
`individuals`, `num_individuals` and `ploidy` are those of the `Variants` of
the VCF. It runs with `num_vars_per_block` of 7 and with the default, on a
file of batches of 100, so the blocks that come out are not the batches of the
file. What the VCF reader gives is checked against bcftools and against pyNei
in its own spec and in `docs/specs/block.md`, and this test carries those
checks over to the file. pyNei is not run here: it reads another file.

A second, at the same function, with `fields` naming the genotypes alone: the
genotypes are the same and the blocks carry no other column.

The cargo tests are the round trips of "The writer", made at `next_block`
over a `Cursor<Vec<u8>>`, and these errors, each on a file built in the test
with arrow-rs: a file with no `popnei` key; a `format_version` of `2.0`, whose
message holds `2.0`; a `gts` width of 7 with 3 individuals and a `ploidy` of
2; a `pos` column of `Int32`; a null position; a file with two batches and one
entry in `popnei_batches`; bytes that are not an arrow file; and
`tests/reference/vars/zstd.vars`, a vars file of the four variants of
`cases.vcf` compressed with zstd, which `tests/reference/vars/make_reference.py`
writes with pyarrow since popnei cannot, which opens and gives the error
at its first `next_block`. Each test checks
the kind of the error and what it names. These are read and are not errors: a
`format_version` of `1.7`; a file with a seventh column, `depth`; a file
written with no compression. And a file written with batches of 100 gives
blocks of 100, which the test checks with no `reblock` in between.

The TypeScript test is the round trip under node of "The writer".

## The Rust interface

What a vars file says about itself, from the `popnei` key of its schema.

```rust
pub struct VarsMetadata {
    /// The whole string, "1.0". Only the part before the dot is checked.
    pub format_version: String,
    pub individuals: Vec<String>,
    pub ploidy: usize,
    pub num_vars_per_block: usize,
}
```

What the footer says of one batch, from `popnei_batches`. It is a batch and
not a block: the blocks that `iter_blocks` gives from a file are cut at the
size the caller asks for, which does not have to be that of its batches.

```rust
pub struct BatchInfo {
    pub num_vars: usize,
    /// One for each chromosome with a variant in the batch, in the order in
    /// which they first appear. Empty in a file with no chrom and pos columns.
    pub regions: Vec<Region>,
}

pub struct Region {
    pub chrom: String,
    /// The smallest and the largest position of the variants of `chrom` in
    /// the batch, both included.
    pub min_pos: u64,
    pub max_pos: u64,
}
```

The reader. `new` reads the schema and the footer, so everything above is
known when it returns, and it fails on what "What it refuses" lists as an
error of the file as a whole.

```rust
pub struct VarsReader<R: Read + Seek> { /* private */ }

impl<R: Read + Seek> VarsReader<R> {
    pub fn new(source: R) -> Result<VarsReader<R>>;
    pub fn metadata(&self) -> &VarsMetadata;
    /// One for each batch of the file, in order.
    pub fn batches(&self) -> &[BatchInfo];
    /// The variants of the whole file, the sum of those of its batches.
    pub fn num_vars(&self) -> usize;
}

impl VarsReader<BufReader<File>> {
    pub fn from_path(path: &Path) -> Result<Self>;
}

impl<R: Read + Seek + Send> BlockReader for VarsReader<R> { /* ... */ }
```

The writer. The columns of the file are fixed by the first block that is
written, so `new` takes only what the `popnei` key needs. It writes one batch
for each block it is given, of whatever size.

```rust
pub struct VarsWriter<W: Write> { /* private */ }

impl<W: Write> VarsWriter<W> {
    /// `num_vars_per_block` is what the `popnei` key will say, the size
    /// the caller gives the blocks it writes.
    pub fn new(
        sink: W, individuals: &[String], ploidy: usize,
        num_vars_per_block: usize,
    ) -> Result<Self>;
    /// It writes `block` as one batch. `chroms` is the table of the reader
    /// the block came from, which holds the names behind its chromosome
    /// numbers. A block of other individuals or another ploidy than those
    /// of `new`, one that fails `Block::check`, and a chromosome number
    /// that `chroms` has no name for are errors, and nothing is written.
    pub fn write_block(&mut self, block: Block, chroms: &ChromTable) -> Result<()>;
    /// It writes the footer, and gives the sink back.
    pub fn finish(self) -> Result<W>;
}

/// Every variant of `reader` into a vars file on `sink`. It asks `reader` for
/// every field and puts a `reblock` of `num_vars_per_block` over it, None for
/// `default_num_vars_per_block` for the individuals of `reader`, which is
/// then the number that the `popnei` key says. This is what both binding crates call. The Python
/// binding crate opens the file, refuses a path that exists, and removes the
/// file when this returns an error.
pub fn write_vars<R: BlockReader, W: Write>(
    reader: R, sink: W, num_vars_per_block: Option<usize>,
) -> Result<W>;
```

The values of the two keys are json, and which crate reads and writes it is
the implementer's choice among those in pure Rust, since the core builds for
wasm.

The cases this module adds to the error of the crate: the source is not a vars
file, with what was found; a format version whose first part is not 1, with
the version; a column of another type, with the column and the two types; a
`gts` width that does not match the `popnei` key, with both widths; a null
where there can be none, with the column and the variant; a footer whose
entries do not match the batches; a file compressed with zstd; a block whose
columns differ from those of the first one written, with the field; a block
that does not fit the writer, with what differs; a chromosome number with no
name; and an error of
the input or the output, which wraps `std::io::Error`. In Python all but the
last are a `ValueError`, and the last an `OSError`.

## Speed

Rust gains little over pyNei here, where the VCF parser of
`docs/specs/io_vcf.md` has 25 times to gain on one thread, the 0.54 s of the
spike against the 13.5 s of pyNei on a VCF of 100000 variants and 1000
individuals. On the
panel of "The compression", 1000 individuals and 20000 variants, reading an
lz4 file of genotypes and summing them took arrow-rs 19 ms against the 20 ms
of pyNei, which reads such a file with pyarrow: both spend the time
decompressing, one in `lz4_flex` and the other in the C of pyarrow, and
neither is bound by the language around it. That measurement was of arrow-rs
alone, a batch at a time, which took 18.8 ms when it was run again beside the
next one. Copying the genotypes of each batch into a vector allocated for it,
which is what the reader of this spec does to give a block, took the same
trial crate 19.1 to 19.3 ms. The number to reach for a pass over that file
with the genotypes alone asked for is that, or no more than a tenth above it,
21 ms on one thread on that machine.

The write has not been measured, in either library, and no number is set for
it here. The implementer measures it once the writer runs, and any work on its
speed waits for that measurement.

## Open points

None. The owner decided on 20 September 2026 the five that there were: a
format of popnei's own that owes pyNei's nothing, lz4 in every build, the name
`open_vars`, the regions of each batch in the footer, and every field in the
file that `write_vars` makes. Each is written where it applies, with the
option that was not taken.

## Not in this spec

- The function that asks a `Variants` for the variants of a region, in Python
  and in TypeScript, and the skipping of the batches outside it. The file has
  what that needs, `popnei_batches`, and the reader gives it as `batches()`. The
  function belongs with the filters, as a later item of
  `docs/specs/filters.md`.
- The read ahead thread that decompresses the next batch while the consumer
  works: with the first calculation that consumes blocks.
- Genotypes packed in 2 bits, the option to measure of section 4 of the
  architecture: it would be a later version of this format.
- Reading a vars file of pyNei, of either of its formats.
- The VCF writer: an item of `docs/specs/io_vcf.md`.
- `Variants`, the handle both functions here give and take:
  `docs/specs/variant.md`. `iter_blocks`, which is how its genotypes reach
  Python and TypeScript: `docs/specs/block.md`.
