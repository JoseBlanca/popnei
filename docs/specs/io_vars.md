# The io::vars module: the vars file, its writer and its reader

September 2026. The vars file is where a user of popnei keeps their variants
once the VCF has been read, so that the text is parsed once and every later
pass reads a file of arrays: the genotypes of 20000 variants of 1000
individuals are decompressed from it in 19 ms. There is code, and none yet
for two things added on 26 September 2026: the lengths of the chromosomes
in the file, and the batches the reader skips for the filter by regions of
`docs/specs/filters.md`. This spec
develops the row `io::vars` of the table in section 9 of
`docs/architecture.md` and section 6 of that document. It covers the format of
the file, the writer and the reader. Two items were added on 9 November
2026 from issue 19, with no code yet: the extra fields, named arrays of a
value for each variant or for each genotype that a file carries beside
its columns, and the writer of segments, which writes a vars file from
variants that a user's code computes one stretch of the genome at a time.
With them the format goes to 1.3, and `write_vars` refuses variants that
are not sorted unless it is told otherwise.

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
| `passed` | `Boolean` | no | whether the FILTER of the variant, in the VCF it was read from, was `PASS` or `.` |

A fixed size list keeps no offsets, so the `gts` column is one flat buffer of
variants x individuals x ploidy signed bytes in that order, which is the
genotype array of a block itself. An allele is 0 for the reference, 1 and
above for the alternative ones and -1 for one that was not called, as
`docs/specs/variant.md` has it; no mask is written, and the -1 is what says
that an allele is missing.

Arrow gives the values inside a list a field of their own, with a name and
with whether they can be null. For its two lists, `alleles` and `gts`,
popnei writes the field pyarrow names for any list, `item`, and says that it
holds no null, which is what is true of them: no allele and no genotype
popnei writes is a null, and an allele that was not called is the -1 above.
A field that can hold nulls costs a bit for each value: arrow-rs writes a
mask of ones beside the genotypes, 1250000 bytes for a batch of 10000
variants of 1000 diploid individuals before it is compressed, which a reader
decompresses and throws away. On the panel of "The compression", 20000
variants of 1000 diploid individuals with the `gts` column alone, the file
of popnei is 15691298 bytes and a pass over it that sums the genotypes takes
20.52 ms, against 15680146 bytes and 19.49 ms for the same genotypes with
the field written as holding no null, the best of 10 interleaved release
runs on one thread of the owner's Apple M5 Pro. The reader compares what a
list holds and not the name of that field nor whether it takes nulls, so a
file of pyarrow, whose `item` takes them, is read as one of popnei is.

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
| `format_version` | `"1.3"` |
| `individuals` | the names of the individuals, in order |
| `ploidy` | how many alleles a genotype holds |
| `num_vars_per_block` | how many variants a batch holds, the last one aside |
| `chrom_lengths` | the length of each chromosome the source gave one for, `[["chr1", 2000], ["chr2", 1500]]`, in its order; `[]` when it gave none |

`chrom_lengths` was added on 26 September 2026, for the density of the
variants of `docs/specs/stats.md` and the `##contig` lines of the VCF
writer of `docs/specs/io_vcf.md`, and with it the version went from 1.0 to
1.1. The writer takes it from the header of its source, which for a VCF
is its `##contig` lines with a `length`. A file of 1.0 has no such key and
is read with no lengths, and a reader of 1.0 reads a file of 1.1 and
ignores the key, by the rule of the version below. A file with the key is not a
vars file, the error of a file as a whole, when the key is not a list of
pairs of a name and a whole number above 0, or when it gives one
chromosome twice, which the VCF reader refuses too; a file without it is
read with no lengths whatever its version says. This was decided on 26
September 2026 with the code.

`passed` was added on 6 October 2026, for the filter of the variants that
passed of `docs/specs/filters.md`, from issue 9, and with it the version
went from 1.1 to 1.2. The writer writes it when its source keeps it: a
VCF, and a vars file of 1.2 that has it, unless the source gives no
variant: the columns of a file are those of its first block, and a file
with no block holds `gts` alone, which nothing reads a `passed` of. A file of 1.0 or 1.1 has no such
column and is read as before, with no `passed` in its blocks, and the
filter of the variants that passed refuses it at the first block, naming
the file, as the owner decided on 6 October 2026. A reader of 1.1 reads a
file of 1.2 and ignores the column, by the rule of the version below. It costs
one bit for each variant before compression.

The extra fields were added on 9 November 2026, from issue 19, and with
them the version went from 1.2 to 1.3: two more entries of the `popnei`
key and a column for each field, which "The extra fields" below has.

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
    variants: Variants, path: str | Path, num_vars_per_block: int | None = None,
    allow_unsorted: bool = False,
) -> VarsWritten
```

`allow_unsorted` was added on 9 November 2026: with its default the call
refuses variants that are not sorted, as "Sorted variants" of the writer
of segments below says, and removes the file. A source with no
chromosomes and positions, a `Variants` built from an array of
genotypes, is not checked.

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
A Ctrl-C that is pending when the call is made is raised before any file is
made, so the path is untouched. One that arrives while the call runs is
raised when the pass over the source is over and not while it runs, and it
takes the file away as an error does. The whole file is written
inside one call of the core, with the interpreter released for all of it, so
Python raises a signal that arrived meanwhile when that call returns. A
Ctrl-C raised between two blocks, which a user gets from `iter_blocks`,
would ask the binding crate to write again the loop over the blocks that the
core has. That the file goes was decided on 21 September 2026 by the session
that ran `docs/plans/vars-file.md`: a user who stopped the call finds the
path free for the call they make again. The option not taken was to keep the
file that the call had finished writing.

The file that an error names is the one the error is about. An error of the
source, a wrong line of the VCF or a read that failed, names the VCF; an
error of the file being written, a disc that filled up among them, names
that path and says that the file could not be written, and so do the two
errors that come before anything is read, a file that is already there and a
path that no file can be made at. The core says which of the two it is: a
write that failed is a case of its own and not the error of a source that
could not be read. A directory at the path is a path no file can be made at
and not a file that is already there, with the number the system gives for a
directory where a file was asked for, so that it is the `IsADirectoryError`
of Python, as it is when `open_vcf` is given one.

The bytes reach the disc before the call returns, and a file system that
refuses them there, a network one that says only when the file is closed
that it is full, is the error of a file that could not be written like any
other. A call that returned would otherwise leave a file that is not whole
and say nothing.

The file of a call that failed and that could not be taken away, a directory
whose permissions changed while the file was being written, is told to the
user as a note on the error they get: what went wrong is what they read
first, and the note says that a file is still at the path, which their next
call would refuse. The session that ran `docs/plans/vars-file.md` decided
these three on 21 September 2026, with the binding crate written; the option
not taken for the last was to say nothing, which leaves a user who fixes
their VCF with a call that refuses the path and no reason.

`num_vars_per_block` is how many variants a batch holds, and `None` is
`default_num_vars_per_block` of `docs/specs/block.md`, so a file read back
with the default size of block gives its batches as they are.

In TypeScript, `writeVars(variants, {numVarsPerBlock, allowUnsorted})` gives back an
object with `bytes`, a `Uint8Array` with the bytes of the file, and
`passStats`. The page offers the bytes as a download:
a tab has no file a path names, as section 11 of the architecture says. With the
option `onBytes` it gives the file in pieces to a function of the page and
no `bytes`, as "The file of a writer in pieces" of
`docs/specs/js_sources.md` says.

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
  seven columns, the ids, the alleles, the qualities and whether each variant
  passed its FILTER, whether or not the user
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

The reference outside the project is pyarrow, the implementation of the arrow
format that Apache Arrow publishes, which opens what popnei wrote as any other
program would, together with bcftools 1.24, whose account of what `many.vcf`
holds is `many.bcftools.tsv` of `docs/specs/io_vcf.md`. pyarrow is a
development dependency of popnei whose lowest version is 23, which the owner
decided on 20 September 2026; the option not taken was a pin at 23.0.0. The
tests run on the pyarrow of `uv.lock`, 25.0.1 on 20 September 2026. The trial
that follows was run on 23.0.0, the pyarrow of pyNei's environment: a trial
file in this format, written with arrow-rs 60, was opened with
`pyarrow.ipc.open_file`, with `pyarrow.feather.read_table` and with pandas.
The three gave its columns and its nulls. `open_file` gave both keys, the one
of the footer as the `metadata` of the opened file; `read_table` gave the key
of the schema alone, and pandas no key.

A pytest test made at `write_vars`, on `many.vcf` of `tests/reference/vcf/`,
the 500 variants of 50 individuals, read with `only_passed=False` and written
with `num_vars_per_block` 100. pyarrow opens the file. Its schema is the seven
columns with the types and the nulls of the table above, in that order. The
value of `popnei` parses as json and holds `format_version` `"1.3"`, the 50
names of the individuals, `ploidy` 2, `num_vars_per_block` 100 and
`chrom_lengths` `[]`, since the `##contig` lines of `many.vcf` have no
length. Written from `write.vcf` of the VCF writer of `docs/specs/io_vcf.md`,
`chrom_lengths` is `[["chr1", 2000], ["chr2", 1500]]`. There are
five batches of 100 variants. The `id` column has 167 nulls and the `qual`
column 100, the variants of `many.bcftools.tsv` with a dot in those columns.
The `passed` column has 25 false, the variants whose FILTER is `q10`, and
written from `many.vcf` read with the default `only_passed`, 475 variants,
none.
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
cargo test writes a block of four variants that are not sorted, chr2 300,
chr1 50, chr2 100, chr1 60, in one batch: its regions are chr2 100 to 300 and
chr1 50 to 60, in that order. The chromosome that comes first in the block is
the one that comes second in the alphabet, so a writer that gave the regions
in any order but the one in which the chromosomes first appear is caught. A third writes `many.vcf` through the VCF reader and reads it
back, and the genotypes, the chromosome names and the positions of the blocks
are those of the blocks of the VCF. A fourth writes a source whose blocks carry
the genotypes alone and finds one column in the file it reads back, and
`num_vars` with no `regions` in its footer. A fifth asserts what the call
says it wrote, 4 variants for the block of `cases.vcf`, 500 for `many.vcf`
at any size of batch and 0 for a source with no variants, and that a
`&mut reader` writes the file that the same reader given whole writes. A
sixth writes 5000 variants of 200 diploid individuals whose every genotype is
`0/1`, reads them back and finds the 2000000 alleles that went in: those
genotypes compress to 0.44 bits each, which is the case of issue 2 above, and
the test asserts that the file is under 250000 bytes so that it keeps holding
that case.

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

A batch of no variants, whose entry of the footer says 0 too, is not given as
a block: the reader takes the next batch, as a filter does with a block it
emptied, because `docs/specs/block.md` says that a reader never gives a block
of no variants. No writer of popnei makes such a batch and another arrow
program can. The session that runs `docs/plans/vars-file.md` decided it on 21
September 2026, with the owner asked and not yet answered; the option not
taken is an error.

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

A null in `chrom`, `pos`, `alleles`, `gts` or `passed` is an error naming the
column and the variant, counted from 1 over the whole file; a null `passed`
read as true or false would keep or drop a variant with no sign. A null among the alleles of a
variant, or among its genotypes, is a null of that column too: arrow gives the
values inside a list a field of their own that can hold nulls, as "What it
holds" says, and no value popnei writes is one. A null `id` is the empty id
and a null `qual` is no quality.

That error is for a file whose schema says the column can hold nulls, which
is what another program writes. popnei writes `chrom`, `pos`, `alleles`,
`gts` and `passed` as columns with no nulls, and arrow-rs refuses a batch of such a column
that holds one before popnei sees it, so a null there is a batch that could
not be read, with what arrow-rs said of it and no variant named: the variant
of a batch that arrow-rs would not decode is not known.

A `qual` that is a value and is not a finite number, a NaN or an infinity, is
an error naming the column and the variant. NaN is what the column of a block
holds for a variant with no quality, so a NaN in the file would be read as a
variant that has none, and an infinite quality is a probability of no variant
of 0, which is not what phred scaling says; the VCF reader refuses both for
the same reason, as `docs/specs/io_vcf.md` has it, and what rests on it is the
rule of "Floats" of the `coding` skill that a NaN in that column means no
quality. No writer of popnei makes such a file and another program can. The
session that ran `docs/plans/vars-file.md` decided it on 21 September 2026;
the option not taken was to let those values into the block, which gives a
calculation over the qualities an infinity to work with and turns a NaN into a
variant with no quality.

An allele of `gts` below `MISSING_ALLELE`, -1, is an error naming the allele
found and the variant. The alleles are signed bytes, so a byte of that column
that was damaged after the file was written, or a file another program wrote,
can say -2, which is neither an allele of the variant nor a missing genotype.
The owner decided on 21 September 2026 that such an allele "is never allowed"
and that every reader of popnei refuses it, "even the VCF parser", after a
reviewer wrote a vars file with popnei, changed one byte of its genotypes to
254, and `open_vars(...).iter_blocks()` handed out a block holding -2 with no
error: the reader checked the nulls, the type and the width of the column and
copied the bytes. What the -2 reached is what `docs/specs/variant.md` says
under "An allele that no reader gives": the counts of one variant refused it,
and a pass that counted nothing put it in the array of a user, where it was
an allele of its own. The option not taken was to leave the reader as it
was and to make those two errors of the counts a `ValueError` instead of the
`RuntimeError` they are, which names no file and no variant and which a pass
with no filter and no count never reaches.

These are errors of the file as a whole, found when it is opened: it is not an
arrow IPC file; its schema has no `popnei` key, or the value is not json, or
one of its four keys is missing; the first part of `format_version` is not
`1`, which the message gives along with the version found; its `individuals`
name nobody, which is a file whose genotypes hold no allele; it has no `gts`
column, which "What it holds" puts in every vars file; it has no
`popnei_batches`, or its entries are not as many as the batches. A path that
is a directory is an error of `from_path`, which looks at the path itself:
opening a directory succeeds on macOS and only the first read fails.

The key that names no individual is the error the writer asked for such a
file gives, with both numbers: every source of popnei has one individual at
least, as `docs/specs/block.md` says, and the blocks of no genotype that such
a file gives are not variants of anybody.

The file without a `gts` column is refused by the session that ran
`docs/plans/vars-file.md`, on 21 September 2026; the option not taken was to
read it as a source whose blocks hold no genotypes, which is what a file of a
later version of the format that dropped the column would be, and no reader of
this version can tell that file from one whose column was lost.

A batch that does not hold the `num_vars` its entry of the footer gives is an
error when that batch is read, so that the number of variants that the file
announces is never a wrong one that goes unnoticed.

A batch whose message does not fit the bytes that came with it is refused
before arrow-rs reads any of them, as a batch that could not be read. arrow-rs
takes the offsets and the lengths of that message as they are: it reads a
buffer at the byte the message gives, and asks the machine for the memory that
the first eight bytes of a compressed buffer say before it decompresses it. So
a file damaged there reaches a panic inside arrow-rs, which in a notebook or a
browser tab ends the session, or an allocation of 144115188075855871 bytes,
which ends the process and which nothing catches. What the reader checks: the
message is one of a batch; the rows it says are the variants of its entry of
the footer; every buffer of it lies inside the body of the batch; every buffer
compressed with lz4 says a length in bytes that its column can hold; and no
column of it says more values than it can hold. Natively, what those checks do
not see is held by `catch_unwind` around the call into arrow-rs, which gives
the same error; in wasm, where a panic ends the program and unwinds nothing,
the checks are the whole of it.

How many bytes a buffer of a column can hold is worked out from the schema of
the file and the rows of the batch, walking the buffers in the order the IPC
format lays them out: the genotypes of a batch are its rows times the
individuals times the ploidy bytes and no more, the positions 8 bytes for each
row, the qualities 4, the offsets of a column of texts and of a list 4 for each
row and one after the last, and a mask of nulls a bit for each row rounded up
to a whole byte, which is what arrow-rs writes. The texts of the three columns
of texts, and the values of a list, have no such number, and neither has a
column whose type popnei does not know, which stops the walk: what bounds those
is what lz4 gives from the bytes the buffer holds, 255 for each byte, and the
2147483647 bytes that the 32 bit offsets of a column of texts address. The
exact bound is what matters for the genotypes, which are the large buffer: a
reviewer wrote a file of 20000 variants of 1000 diploid individuals in one
batch, 25147258 bytes, and changed the eight bytes that say how long its `gts`
buffer is once it is decompressed. Under wasm, on 21 September 2026, 1000000000
and 2000000000 gave the error of a batch that could not be read and left the
memory of the tab grown to 2066087936 bytes for its life, because arrow-rs had
asked for what the buffer said before it read it; 3000000000 and 4294967295
ended the tab with a trap, which is what wasm does with an allocation it cannot
address.

How many values a column can hold comes from the same schema and the same rows.
It is a second check and not the first one again: a buffer of the message says
how many bytes it holds once it is decompressed, and a field node of the same
message says how many values its column holds, and arrow-rs reads both. For the
genotypes the two numbers are the same, the rows times the individuals times
the ploidy, because an allele is one byte. The walk is over the columns in the
order the IPC format lays them out, a column inside a column after the one that
holds it: every column of the batch holds its rows, the alleles of the
genotypes that number, and the alleles of a variant, which are the values of a
list, as many as the offsets of that list say, at most the 2147483647 that 32
bit offsets address, the same number as the bytes above and from the same
offsets. The values of a large list are the one thing the schema does not
bound at all, since its offsets are 64 bits.

A second bound holds every column, and it is the smaller of the two that a node
is held to: the body of the batch. A value takes a bit at the very least once
the batch is decompressed, and an lz4 frame gives at most 255 bytes for each
byte it holds, so a body of 110080 bytes holds 224563200 values and no more.
For every column of a file popnei writes it is looser than the schema's by a
factor of eight or more, and it is the only bound a column whose type popnei
does not walk has, and every column after it: popnei cannot say how many field
nodes such a column takes, so it cannot say which column the nodes after it
belong to.

The schema alone is not enough, because the rows a column is bounded by are
what the batch and its entry of the footer say, and no byte of the file bounds
those. A file whose first column is a dictionary, which pyarrow writes for any
categorical of pandas, and whose alleles of the genotypes say
3074457345618258603 values, reaches `integer overflow computing expected number
of expected values in FixedListSize` inside arrow-rs, an `expect` that panics
in a release build too. The code review of 24 September 2026 found that with
such a file after the bits of the batch had been replaced by the schema instead
of joined to it.

Until 24 September 2026 the body was the only bound, and it counted the bits of
the batch as it lies on disk, where it is compressed, against values that are
counted once it is decompressed. So popnei refused the file it had just written
when its genotypes repeated: 5000 variants of 200 diploid individuals whose
every genotype is `0/1` are 2000000 alleles in a batch of 110080 bytes, 0.44
bits for each allele. That is issue 2 of the repository. What the 255 bytes an
lz4 frame gives for a byte do is put the two sides of that comparison in the
same state.

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
decompression of its four batches. Before the copy, one pass over that
buffer takes the smallest of its alleles and refuses the batch when it is
below the missing one, which "What it refuses" asks for: the smallest of a
run of bytes is what a compiler reduces over the lanes of a vector
register, where a comparison written for each allele on its own would not
be. The reader gives the 4e7 alleles of that panel, 1000 individuals and
20000 variants, in 20.50 to 20.61 ms with the genotypes alone asked for,
which is the number to hold against the 21 ms of "Speed" below and leaves
0.4 ms under it. The check is 0.35 to 0.54 ms of that, 2 in 100: the same
reader without it gave the pass in 20.07 to 20.26 ms, and the difference is
of each of six pairs of runs, the two builds one after the other, three
pairs on each of two occasions. Each number is the best of 5 runs of
`crates/popnei/benches/vars_file.rs` on the owner's Apple M5 Pro with a
load average of 1.3 to 1.4. The positions and the qualities are
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

The reader answers `skip_outside` of `docs/specs/block.md` with true. Once
it is handed the regions of the filter by regions, it does not read a
batch of which `keeps_none_of` of the selection says, for each chromosome
of the batch in `popnei_batches`, from its smallest to its largest
position, that no variant can be kept: it does not seek to it, nothing of
it is decompressed, and its variants go into `num_skipped`. The message of
the batch, the few bytes before its buffers, is read, and a batch whose
message holds another number of rows than its entry of the footer is the
error of such a batch, as it is when the batch is read. Beyond that the
footer is trusted for a batch that is skipped: a footer whose regions of a
batch were changed after the file was written, so that they no longer hold
its variants, has the batch skipped although the filter would keep some of
them, and nothing says so. A batch that is read with its chromosomes and
its positions, which the filter by regions always asks for, is checked
against its entry: a variant of it on a chromosome the entry does not
name, or at a position outside the smallest and the largest the entry
gives for its chromosome, is an error of the batch, a `ValueError` in
Python that names the file, the batch, the chromosome and the position, and so is a region of the footer
whose smallest position is above its largest, when the file is opened. The
owner is asked whether the footer is to be trusted for the batches that
are skipped. A file with no
`chrom` and `pos` columns has no regions in its footer, so the reader
skips nothing and gives its first batch, which the filter refuses with the
error of a field it depends on. `header()` gives the individuals of
the file, `chrom_lengths` of its `popnei` key, its extra fields, and no
meta lines of a VCF.

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
2; a file with no `gts` column; a file whose `individuals` name nobody and
whose `gts` holds no allele; a `pos` column of `Int32`; a null position; a
quality that is a NaN and one that is an infinity, both with the variant; a
file with two batches and one entry in `popnei_batches`; a message whose
column of the alleles of the genotypes says 25 values where the 4 variants of
3 diploid individuals of the batch hold 24, which is the bound on the values a
column says it holds at its edge; a file whose first column is a dictionary,
which stops the walk of the schema, and whose alleles of the genotypes say
3074457345618258603 values, which the body of the batch refuses and which
never reaches arrow-rs;
bytes that are not an arrow file; and
`tests/reference/vars/zstd.vars`, a vars file of the four variants of
`cases.vcf` compressed with zstd, which `tests/reference/vars/make_reference.py`
writes with pyarrow since popnei cannot, which opens and gives the error
at its first `next_block`. Each test checks
the kind of the error and what it names. These are read and are not errors: a
`format_version` of `1.7`; a file with a column it does not know, `depth`; a file
written with no compression. And a file written with batches of 100 gives
blocks of 100, which the test checks with no `reblock` in between.

A cargo test sweeps a vars file of the four variants of `cases.vcf` written in
the test: every byte of it set in turn to four values, the low bit and the
high bit flipped, 0 and 255, and each of the files that makes read with every
field. Each gives its variants or an error, and none reaches a panic that
comes out of the reader or an allocation that ends the process, which an
abort would take the test binary with it. The same sweep with every byte set
to each of the 255 other values is a test that is run by hand, as
`docs/specs/io_vcf.md` has one for `cases.vcf.gz`. A file that is read with no
error and holds other variants is counted and is not a failure: a byte of a
compressed buffer that decompresses into other genotypes is what a checksum of
the format would catch, and the format has none. On 24 September 2026, with
arrow-rs 60, the sweep over the 255 values made 1299990 files, of which 553104
gave an error, 726033 were read as the whole file, 20853 were read as another
file with no error, 1939 reached a panic inside arrow-rs that `catch_unwind`
held, and none ended the process. The panics were 2783 while a column of a
batch was bounded by the bits of that batch: a column holds its rows, which is
far below that, so 844 more of these files are refused before arrow-rs reads
them.

A cargo test walks a schema of the types popnei's writer never makes and its
own tests otherwise never reach, a struct, a large list, a map and a
dictionary, and asserts one bound for each field node the IPC format lays out,
in its order: the columns of a struct hold the rows of the struct, the values
of a large list are not bounded, those of a map are bounded by its 32 bit
offsets, and the list ends at the dictionary, which is what leaves the columns
after it to the body of the batch.

The TypeScript test is the round trip under node of "The writer".

## The extra fields

Added on 9 November 2026 from issue 19, with no code yet.

### What they give

An extra field is a named array of numbers that a user attaches to the
variants of a vars file beyond the columns of the table of "What it holds".
A variant field has one value for each variant, as an INFO field of a VCF:
the end of a haplotype region, whose position is its start. A genotype
field has one value for each genotype, that is for each individual at each
variant, as a FORMAT field of a VCF: the posterior probability of the
genotype that a caller gave, or the p-value of a test of how well it fits.
popnei does not read the values. They travel with the variants, through
the filters and into the files that `write_vars` writes, and come back out
of `iter_blocks` when they are asked for by name, so that the confidence of
each call stays with the call instead of in a table kept beside the file.

The case they were asked for is the haplotype analysis of the tomato25
panel. Each region of a chromosome is one variant whose alleles are the
haplotypes seen in it, from a few to about 90, and each individual has
seven values that say how far to trust its genotype there: the posterior
of the pair of alleles, the probability of being heterozygous, the two
nodes of the coarser call, its posterior, how divergent the alleles of
its groups are, and the p-value of the test of fit. The pair of nodes of
the tree of a region's alleles that the coarser call gives is two genotype
fields, one for each node, since a field holds one number for each
genotype. The owner decided this on 9 November 2026; the option not taken
was a field of a fixed number of values for each genotype, which only that
pair would have used.

A field has one of five types: `int32`, `int64`, `uint64`, `float32` and
`float64`. A missing value is the user's to choose, NaN in a field of
floats and some number in a field of integers, since popnei does not read
the values. Its name starts with a letter or `_` and holds letters,
digits, `_` and `.` after it, the rule of the VCF for the key of an INFO
or a FORMAT field, so that a field can become one. It is not the name of a
column of the table of "What it holds", and no two fields of a file share
one, of either kind.

### In the file

The vars file of format 1.3 has one column for each extra field after
those of the table of "What it holds", the variant fields first and then
the genotype fields, each in the order it was declared, and named as the
field is, so that pandas shows `END` and `GP` as columns:

| field | arrow type | nulls |
|---|---|---|
| variant field | `Int32`, `Int64`, `UInt64`, `Float32` or `Float64` | no |
| genotype field | `FixedSizeList<that type>[num_individuals]` | no |

A genotype field is a fixed size list, as `gts` is, so its column is one
flat buffer of variants x individuals values. The `popnei` key of the
schema gains two entries, `var_fields` and `gt_fields`, each a list of the
name and the type of each field, in the order of the columns:

    "var_fields": [["END", "uint64"]],
    "gt_fields": [["GP", "float32"], ["NODE_A", "int32"], ["NODE_B", "int32"]]

The schema is written before the first batch, so the extra columns are
there in a file with no batch too, and a user reads which fields a file
has without reading a variant. Both entries are written in
every file of 1.3, as `[]` when it has no field of that kind. A file
without them, any file of 1.2 or before, has no extra field.

Every file that popnei writes from now on says `"1.3"`, with or without
extra fields. A reader of 1.2 reads a file of 1.3 and ignores the two
entries and the extra columns, by the rule of the version under "What it
holds"; a reader of 1.3 reads every file of 1.0 to 1.2.

The reader refuses, as an error of the file as a whole, a name in the two
lists without its column, a column of another type than its entry says or
of a genotype field whose width is not the number of individuals, a name
that does not follow the rule above, one that is a column of the table of
"What it holds", and one that is in the lists twice. An entry that is not
a list of pairs of a name and one of the five types is the error of a
file that is not a vars file, as one of `chrom_lengths` is; the session
that ran task 1.2 of `docs/plans/vars-from-segments.md` decided it on 9
November 2026. A null in an extra
column is the error of a null of "What it refuses". The bounds that
"What it refuses" puts on the buffers of a batch before arrow-rs reads it
hold the extra columns too: a variant field is the rows of the batch times
the bytes of its type, 4 or 8, and a genotype field that times the
individuals.

A `float32` genotype field is 4 bytes for each genotype before
compression, where the genotype itself is one byte for each allele, 2 for
a diploid; the seven fields of the tomato25 case, all in four bytes, are 28.
How those values compress has not been measured.

### What carries them

- `retain_vars` of `docs/specs/block.md`, which every filter of variants
  calls, compacts the extra columns with the others.
- The filter of individuals of `docs/specs/filters.md` keeps, in each
  genotype field, the values of the kept individuals in the order it keeps
  them, as it does with the genotypes. A variant field is left as it is.
- `reblock` cuts and joins them with the rest of the block. It joins two
  blocks only when they have the same columns and the same extra fields,
  in the same order, and otherwise gives the one that was waiting as a
  shorter block, as `docs/specs/block.md` says; the writer then refuses a
  block whose fields are not those it was given, so every reader gives
  each of its blocks the same columns and fields, which the reader of
  segments below does by filling every column.
- `write_vars` writes them: the extra fields of a file are those of the
  header of its source, which for a vars file are those of its key, so a
  file written from another file keeps its fields, after any filter.
- `write_vcf` of `docs/specs/io_vcf.md` refuses a source that has extra
  fields with a `ValueError` naming them, before it makes the file. The
  owner decided it on 9 November 2026; the option not taken was to write
  each variant field as an INFO field and each genotype field as a FORMAT
  field, which can come when a user asks for a VCF of such a file. A user
  with such a file gets no VCF of it, since no step takes a field away.
- No calculation reads them, so none asks for them and a vars file reader
  asked by a calculation decompresses none of them.

The number of variants of a block that popnei picks when it is asked for
none, `default_num_vars_per_block` of `docs/specs/block.md`, counts the
genotypes and not the extra fields: about 5 million genotypes, at most
10000 variants. So a block that carries the seven fields above, of four
bytes each, takes 30 bytes for each diploid genotype instead of 2, 150 MB
for 5 million genotypes. It is decided here that the default does not change: only
`write_vars` and an `iter_blocks` that names the fields read them, and both
take a `num_vars_per_block`. For the tomato25 case, 10000 regions of 300
individuals, a block is 90 MB.

### In Python and in TypeScript

```python
Variants.var_fields: dict[str, numpy.dtype]
Variants.gt_fields: dict[str, numpy.dtype]
```

The extra fields of the source of a `Variants` and the type of each, in
the order of the file, read from the header when it was opened. Both are
empty for a VCF and for a vars file of 1.2 or before.

`iter_blocks(fields=...)` of `docs/specs/block.md` takes the names of the
extra fields of the source beside the five it takes now:
`fields=("chrom", "pos", "GP", "END")`. A field that is not named is not
read, so the default `fields` gives none. A name that is neither one of the
five nor a field of the source is the `ValueError` it gives now, and its
message names the fields of the source. `Block` gains one attribute,
`extra_fields`, a dict of the name of each field that was asked for to a
read only numpy array of its type, of variants for a variant field and of
variants x individuals for a genotype field. It is `None` when no extra
field was asked for, as a column that was not asked for is. The two kinds
share one dict because no two fields share a name.

No class of the package changes in any other way, which the owner asked
for on 9 November 2026: `Block` keeps its attributes and gains one with a
default of `None`, `Variants` gains two properties, and `iter_blocks`
keeps its arguments, whose `fields` is annotated as names of any text
instead of the five of `Field`, which stays for the five. The `repr` of a
block names the extra fields it holds after its columns.

In TypeScript there is nothing of the extra fields yet: a file that has
them opens, its calculations run and `writeVars` keeps its fields, since
all of that is the core, but a TypeScript user cannot read them. The owner
decided on 9 November 2026 that the TypeScript functions of the extra
fields and of the writer of segments wait for a user who needs them; the
option not taken was to build them with the Python ones, which keeps the
two APIs the same, as the objectives ask, at the cost of turning every
typed array of a segment into a block in the wasm binding crate.

pyNei has no genotype fields. Its `VariantsChunk` of `pynei/variants.py`
holds a pandas frame of the variants, `vars_info`, with any columns a user
puts in it, and its `write_vars` writes them, so its variant fields are
those columns.

### How it is verified

The tests are those of the writer of segments, below, which is the one
source of extra fields, and these cargo tests: `Block::retain_vars` on a
block of 4 variants of 3 individuals with a variant field and a genotype
field keeps the rows of the variants that are kept in both; the filter of
individuals keeping the third and the first individual, in that order,
gives each row of the genotype field as the third value then the first;
`reblock` of blocks of 3 and 2 variants into blocks of 4 and 1 gives the
fields joined in order; and `Block::check` refuses a genotype field with
one value too few. A reader test for each error of the file above, on a
file built with arrow-rs. The sweep of "How it is verified" of the reader,
every byte of a small file set to four values, runs on a second file too:
the four variants of `cases.vcf` with a `float32` genotype field and an
`int64` variant field.

## The writer of segments

Added on 9 November 2026 from issue 19, with no code yet.

### What it gives

It writes a vars file from variants that a user's code computes, one
genomic segment at a time, with no VCF in between and without the whole
genome in memory at once. A segment is a run of consecutive variants that
the code gives as a dict of arrays: the caller of haplotypes that issue 19
was written for walks a genome region by region and yields the regions it
has called. Once written, the file is opened with `open_vars` and every
filter and calculation of popnei runs on it as on any other.

The owner chose it on 9 November 2026. The options not taken were a
`Variants` built from arrays in memory, which holds the whole dataset, and
a function that writes to a file object of Python instead of a path, which
leaves the file to the caller to remove when the call fails.

### Its Python function

```python
def write_vars_from_segments(
    segments: Iterable[Mapping[str, Any]],
    path: str | Path,
    individuals: Sequence[str],
    ploidy: int,
    var_fields: Mapping[str, str | numpy.dtype] = {},
    gt_fields: Mapping[str, str | numpy.dtype] = {},
    chrom_lengths: Mapping[str, int] | None = None,
    num_vars_per_block: int | None = None,
    allow_unsorted: bool = False,
) -> VarsWritten
```

Each item of `segments` is a dict with these keys:

| key | value | |
|---|---|---|
| `chrom` | a sequence of the chromosome of each variant, or one `str`, the chromosome of every variant of the segment; each segment can use either | required |
| `pos` | integers, the position of each variant, 1 based | required |
| `alleles` | a sequence with the sequence of the alleles of each variant, each one a text, the reference first | required |
| `gts` | a numpy array of any integer type, of three dimensions, variants x individuals x ploidy, -1 for a missing allele | required |
| `id` | a sequence of a text or `None` for each variant | optional |
| `qual` | floats, NaN for a variant with no quality | optional |
| each declared field | a variant field of variants, a genotype field of variants x individuals | required |

The number of variants of a segment is the length of `pos`, and every
other value has to agree with it. A segment of no variants is taken and
writes nothing. An allele is any text that is not empty, so the haplotypes
of a region are `H0`, `H1` and so on, with no made up nucleotides: the
vars file keeps an allele as text, and bcftools 1.25, run on 9 November
2026 on a VCF with `REF` `H0` and `ALT` `H1,H2`, reads it and writes it
back unchanged, so `write_vcf` writes such alleles as they are. An absent
`id` or `qual` writes no id or no quality for those variants.

`individuals`, `ploidy`, `chrom_lengths` and the two declarations of the
fields go in the schema of the file, which is written before the first
segment is taken, and that is why they are arguments and not read from the
first segment. A field is declared by its name and its type, any spelling
of one of the five types that `numpy.dtype` takes. What a field of each
type takes, as a numpy array or as anything `numpy.asarray` turns into
one, so a list too:

- a field of integers takes integers, of any integer type, and a value
  that the declared type does not hold, a negative one in a `uint64`
  field or 2^31 in an `int32` one, is a `ValueError` naming it, since
  numpy would wrap it into another number with no sign; floats are a
  `ValueError`, also when they are whole;
- a field of floats takes floats and integers, which are converted to the
  declared type, so a `float32` field rounds a `float64` value to the
  nearest `float32`, which is what declaring it `float32` asks for.

The rule of numpy for a safe cast was not taken: it refuses an `int64`
array, the type numpy gives a list of integers, into a `uint64` field,
and a `float64` array, which a list of floats gives, into a `float32`
field, and it takes an `int64` array into a `float64` field, which rounds
above 2^53.

`num_vars_per_block` is the size of the batches of the file, as for
`write_vars`, and the segments are cut and joined to it, so a segment can
hold any number of variants. What it gives back is `VarsWritten`, whose
`pass_stats` has the number of variants written and an empty `filtering`,
since no filter ran.

The memory it uses is the segment it was last given and the variants that
wait for the next batch, `num_vars_per_block` at most, with the batch
being compressed. The size of the genome does not change it.

The `Variants` of the file it writes has `keeps_passed` false: no FILTER
was read, and the file has no `passed` column.

### Sorted variants

The variants have to be sorted, which the owner asked for on 9 November
2026: the positions of a chromosome do not go down from one variant to the
next, and the variants of a chromosome come together, so a chromosome
whose variants have ended does not come back. Two variants at one position
are in order, whichever comes first. The chromosomes themselves can come in any order. The check
runs over the whole file, across segments and across batches, and costs
one comparison for each variant. A variant out of order is a `ValueError`
that names it and the one before it, with the number of each among the
variants written, counted from 1 over the file that is being written,
`the 2nd variant written, chr1 100, comes after the 1st, chr1 5000`, and
the file is removed. The numbers count the variants that reach the
writer, after any filter, so they are not the places of the two variants
in a source that a filter took variants out of; the message says so.
`allow_unsorted=True` writes the variants in the order
they come, and the message ends by naming it, as `allowUnsorted: true` in
TypeScript. The option not taken was to write any order with no check,
which is what the format allows and what `docs/specs/io_vcf.md` does with
the VCF it writes.

`write_vars` gains the same argument with the same default, and refuses a
source that is not sorted in the same way, which the owner accepted on 9
November 2026 with the rest of this design. That changes what it does with
a VCF that is not sorted, which it wrote until now and refuses from this
version on unless the caller passes `allow_unsorted=True`. The tests that
write variants that are not sorted pass it: the two of
`tests/test_io_vars.py` over the fixture `long_vcf`, which is `many.vcf`
repeated, and the cargo tests of `crates/popnei/src/io/vars.rs` that write
the block of four variants not sorted of "How it is verified" of the
writer.

### A call that fails

It does what `write_vars` does: a path that already exists is a
`ValueError` and nothing is written, and when anything fails after the
file was made, the file is removed, so a failed call leaves nothing at the
path and the same call can be made again. Since the path is never
overwritten, a call that fails cannot destroy a file that was there
before.

An exception raised by `segments`, the user's code that computes them,
comes out of the call as it was raised, with its type and its traceback,
after the file is removed. The core cannot carry an exception of Python,
so `SegmentReader` keeps it, ends the pass with an error of the core of
its own, which says that the source of segments failed, and the binding
crate raises the exception it kept in the place of that error. A Ctrl-C is raised when the next segment is
asked for, and removes the file as an error does. That is sooner than in
`write_vars`, which raises a Ctrl-C only when its whole pass is over:
both release the interpreter while the core writes, and this call takes
it back for each segment, since every segment comes from Python.

A segment that is wrong is a `ValueError` that names it, counted from 1,
and the variant of it, counted from 1, with its chromosome and its
position when those are not what is wrong. These are:

- a key that is neither one of the six of the table nor a declared field,
  and a required one that is missing;
- a value whose length or shape does not agree with the number of
  variants, the individuals and the ploidy;
- a `pos` that is not of integers or holds a negative one, and `gts` that
  are not of integers;
- an allele of `gts` below -1, or not among the alleles of its variant, 2
  for a variant of two alleles;
- a variant with no allele, an allele that is an empty text, and a variant
  of more than 128 alleles, numbered 0 to 127, the most that one byte holds, which the owner
  kept on 9 November 2026; the option not taken was alleles of two bytes,
  which doubles the genotypes in memory and in the file before
  compression, changes every calculation that reads them and is a new
  first part of the version of the format;
- a quality that is an infinity, which the reader would refuse;
- a field array that its type does not take, by the rule of "Its Python
  function";
- a tab, a newline, a carriage return or a space in a chromosome, an id or
  an allele, a comma in an allele and a semicolon in an id, which would
  break the columns and the lists of a VCF that `write_vcf` writes from
  the file, as the VCF reader would split them.

An item that is not a mapping is a `TypeError`. The arguments are checked
before the file is made, each a `ValueError`: no individual, two of one
name, a ploidy of 0 or above 255, a field of another type than the five, a
name that breaks the rule of "The extra fields", a name declared twice,
and a length of a chromosome that is not above 0.

### What pyNei does

pyNei writes a vars file from code of the user through its `Variants`,
which is built over any object with an `iter_vars_chunks` that yields
`VariantsChunk`s, read in `pynei/variants.py` at ef0ca6e; its `write_vars`
of `pynei/io_vars.py` writes the columns of the first chunk. The
differences: popnei takes dicts and writes them in one call, the fields
are declared, there are genotype fields, the variants are checked to be
sorted, a source of no variants is written, and the file of a call that
failed is removed, where pyNei leaves it.

In TypeScript there is no writer of segments, by the decision of "The
extra fields" above; `writeVars` gains `allowUnsorted`, since the check of
the order is the core's and holds for both packages.

### How it runs

The Python binding crate wraps `segments` in `SegmentReader`, a
`BlockReader` of `docs/specs/block.md` of its own, whose `next_block`
takes the next item from Python, checks it and turns it into a block with
every column and every declared field, an empty id and a NaN quality
where the segment has none, so that every block has the same columns, and gives that reader to `write_vars` of the core, as for
any other source. The reader holds the chromosome table of the names it
has met, its `individuals()` and `ploidy()` are the arguments, and its
`header()` holds the individuals, `chrom_lengths` and `extra_fields`, the
declarations in the order of "In the file". `write_vars` of the core
gives the lengths and the fields of that header to `VarsWriter` before the
first block, as it does the lengths for every source, puts the `reblock`
over the reader, and checks the order on each block that `reblock` gives,
before the block is written: it keeps the chromosome and the position of
the last variant and the numbers of the chromosomes that have ended.
`VarsWriter` itself does not check the order, so a caller of the core
that writes blocks with it directly writes them in any order. A block
without chromosomes and positions is not checked. The binding crate builds
the `pass_stats` of the result from the count that `write_vars` gives, as
it does for `write_vars`. The pass
releases the interpreter, as `write_vars` does, and the reader takes it
again for each item it asks for. A segment of no variants is skipped, since
a reader never gives a block of none.

### How it is verified

The reference outside the project is pyarrow, which opens the file as any
other program would. A pytest test made at `write_vars_from_segments`,
with the individuals `i1`, `i2` and `i3`, ploidy 2, `var_fields`
`{"END": "uint64"}`, `gt_fields` `{"GP": "float32", "NODE_A": "int32"}`
and `num_vars_per_block` 2, on three segments:

| segment | chrom | pos | alleles | gts | END | GP | NODE_A |
|---|---|---|---|---|---|---|---|
| 1 | `"chr1"` | 100 | H0 H1 | 0/1 1/1 ./. | 4999 | 0.9 0.8 NaN | 1 2 -1 |
| | | 5000 | H0 H1 H2 | 2/0 0/0 1/2 | 9999 | 0.7 0.99 0.5 | 3 1 2 |
| 2 | `"chr1"` | none | | | | | |
| 3 | `"chr2"` | 50 | H0 | 0/0 0/0 0/. | 800 | 1.0 1.0 0.6 | 1 1 1 |

It writes 3 variants. pyarrow finds the columns `chrom`, `pos`, `id`,
`alleles`, `qual`, `gts`, `END`, `GP` and `NODE_A`, in that order, with
the types of the two tables above, no `passed` column, an `id` and a
`qual` that are null in all 3 rows, and the values that went in; the
`popnei` key has `format_version` `"1.3"`, `var_fields` `[["END",
"uint64"]]` and `gt_fields` `[["GP", "float32"], ["NODE_A", "int32"]]`;
`popnei_batches` has 2 batches, of 2 variants with chr1 100 to 5000 and of
1 with chr2 50 to 50. `open_vars` of the file has `var_fields` `{"END":
uint64}` and `gt_fields` `{"GP": float32, "NODE_A": int32}`, and
`iter_blocks(fields=("chrom", "pos", "GP"), num_vars_per_block=3)` gives
one block whose `extra_fields` holds `GP` alone, 3 x 3, equal to what
went in with NaN in its place. After `filter_individuals(["i3", "i1"])` the
`GP` of that block is NaN 0.9, 0.5 0.7 and 0.6 1.0, and `write_vars` of the
filtered `Variants` gives a file whose `GP` reads back as that. `write_vcf`
of it is the `ValueError` of "What carries them" and names `END`, `GP` and
`NODE_A`.

One pytest test for each `ValueError` and the `TypeError` of "A call that
fails", each of which leaves no file at the path, and among them: chr1 5000
then chr1 100, whose message holds both; chr1, chr2, chr1, refused although
the positions go up, and both written with `allow_unsorted=True`; an
allele 2 for a variant of the alleles H0 H1; a `float64` array for `GP`.
A generator that yields the first segment above and then raises
`RuntimeError("caller failed")` gives that `RuntimeError` and no file. A
test of `write_vars` on a VCF written in the test whose second variant is
at a lower position than its first gives the `ValueError` and no file, and
writes it with `allow_unsorted=True`. The tests of `write_vars` on
`many.vcf`, which is sorted, change only in the `format_version` they read,
`"1.3"`, and the two entries of the extra fields, `[]`; so does
`tests/test_filter_passed.py`, which changes the bytes of
`"format_version":"1.2"` in a file.

A test of a Ctrl-C between segments: a generator that calls
`_thread.interrupt_main()` after it yields its first segment gives a
`KeyboardInterrupt` from the call and no file. A cargo test at `write_vars`
of the core over a reader made in the test that gives chr1 5000 and then
chr1 100 gives the error with both variants and, with `allow_unsorted`
true, the file. In TypeScript, a test under node of `writeVars` with a
VCF that is not sorted: an `Error` by default and the bytes with
`allowUnsorted`.

pyNei is not run: its file is another one, and it has no genotype fields.

## The Rust interface

What a vars file says about itself, from the `popnei` key of its schema.

```rust
pub struct VarsMetadata {
    /// The whole string, "1.3". Only the part before the dot is checked.
    pub format_version: String,
    pub individuals: Vec<String>,
    pub ploidy: usize,
    pub num_vars_per_block: usize,
    /// Empty in a file of 1.0, which has no such key.
    pub chrom_lengths: Vec<(String, u64)>,
    /// The extra fields, from `var_fields` and `gt_fields`, in the order
    /// of their columns. Empty in a file of 1.2 or before. Added on 9
    /// November 2026.
    pub extra_fields: Vec<ExtraFieldDecl>,
}
```

`ExtraFieldDecl`, the name, the kind and the type of one extra field, is
in `docs/specs/block.md`, since the header of every source carries it.

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

Which extra fields the reader decompresses when `EXTRA` of `Needs`, the
flag of `docs/specs/variant.md`, is asked for: the ones named, every one for `None`, which is what it does
until this is called. It holds from the next block, as `set_needs` does.
A name the file has no field of is an error naming the fields it has.
The fields of a block are in the order of the file whatever the order of
the names, and a block of a selection that names none has no extra
field, `extra` of `None`, as one of a file that has none. The extra
fields of the header of the reader become the selected ones, in the order
of the file, so that `write_vars` of a selected reader writes a file of
those fields and does not refuse the blocks for lacking the others; every
field of the file stays in its metadata, which is where `var_fields` and
`gt_fields` of `Variants` read them when the file is opened. It
is a method of the vars reader and not of `BlockReader`: the binding crate
calls it on the source it opened for `iter_blocks`, before it puts the
filters over it, and `write_vars` never calls it.

```rust
impl<R: Read + Seek> VarsReader<R> {
    pub fn select_extra_fields(&mut self, names: Option<Vec<String>>) -> Result<()>;
}
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
    /// The extra fields of the file, which the `popnei` key will say and
    /// which every block written has to carry, with these names, kinds and
    /// types, in the order the blocks hold them; the columns are the
    /// variant fields in that order and then the genotype fields. None
    /// until this is called. It holds when it is called before the first
    /// block, as the lengths of the chromosomes do. A name the rule of
    /// "What they give" refuses, one of a column of "What it holds" and
    /// one given twice are the error of a field of the file. A call after
    /// the first block with the fields already declared changes nothing,
    /// and one with other fields, which would change a schema already
    /// written, is the error of a block whose fields are not those declared.
    pub fn set_extra_fields(&mut self, fields: Vec<ExtraFieldDecl>) -> Result<()>;
    /// It writes `block` as one batch. `chroms` is the table of the reader
    /// the block came from, which holds the names behind its chromosome
    /// numbers. A block of other individuals or another ploidy than those
    /// of `new`, one that fails `Block::check`, and a chromosome number
    /// that `chroms` has no name for are errors, and nothing is written.
    pub fn write_block(&mut self, block: Block, chroms: &ChromTable) -> Result<()>;
    /// It writes the footer, and gives the sink back.
    pub fn finish(self) -> Result<W>;
}

/// Every variant of `reader` into a vars file on `sink`, and the sink back
/// with how many variants were written. It asks `reader` for
/// every field and puts a `reblock` of `num_vars_per_block` over it, None for
/// `default_num_vars_per_block` for the individuals of `reader`, which is
/// then the number that the `popnei` key says. The extra fields are those
/// of the header of `reader`. With `allow_unsorted` false, a variant out of
/// the order of "Sorted variants" is an error. This is what both binding
/// crates call, for `write_vars` and for the writer of segments. The Python
/// binding crate opens the file, refuses a path that exists, and removes the
/// file when this returns an error.
pub fn write_vars<R: BlockReader, W: Write>(
    reader: R, sink: W, num_vars_per_block: Option<usize>, allow_unsorted: bool,
) -> Result<(W, u64)>;
```

The count is the `num_vars` of the `PassStats` that every consumer of a
`Variants` gives back, the owner's decision of 21 September 2026 that "Its
Python and TypeScript functions" of the writer has for this one. It comes
from here because the loop over the blocks is here: a binding crate that
calls this one sees no block of the pass and can count nothing. A caller
that has to give its user the counts of the filters of that pass as well
gives `&mut reader`, which is a `BlockReader` too, as
`docs/specs/block.md` has it: the chain of readers stays with the caller,
which reads `filtering_stats` from it when this returns, as "How it runs"
of the counts of `docs/specs/filters.md` asks of every consumer.

The values of the two keys are json, and which crate reads and writes it is
the implementer's choice among those in pure Rust, since the core builds for
wasm.

The cases this module adds to the error of the crate, with the exception each
one is in Python. The owner gave the convention on 21 September 2026: a
`ValueError` is a wrong input of a function, a `RuntimeError` a defect of
popnei, and an `OSError` a file that cannot be read, that was cut short or
that is corrupted.

Thirteen are a `ValueError`, since a file whose content is not what a vars file
holds is a wrong input like a wrong argument: the source is not a vars file,
with what it lacks, which is the header of an arrow file, the `popnei` key of
the schema, the json of its value, one of that key's four values, the `gts`
column or the
`popnei_batches` key of the footer; a format version whose first part is not
1, with the version found; a column of another type, with the column and the
two types; a `gts` width that does not match the `popnei` key, with both
widths; a null where there can be none, with the column and the variant; a
`qual` that is a value and is not finite, with the value and the variant; an
allele of `gts` below the missing one, with the allele and the variant; a
footer whose entries are not as many as the batches, with both counts; a
batch that holds another number of variants than its entry of the footer,
with the batch and both counts; a file whose buffers are compressed with
zstd; two individuals of one name, with the name; a file of no individual or
of the ploidy 0, with both numbers, which is a writer asked for such a file
and a reader of one whose `popnei` key names nobody; and a block with more
text in one of its
columns than the 2147483647 bytes an arrow column of texts holds, or more
alleles in its `alleles` column than the 2147483647 a list column of a
batch holds, with the column, what it holds and that number. The last one
is the size of a
batch and not of the file: the way out is a smaller `num_vars_per_block`,
which the message says. arrow-rs panics when a column of texts goes past
it, and when the entries of a list column do, so the writer counts the
bytes of each of the three columns of texts of a block before it fills
one, and the alleles of its `alleles` column, whose entries are more than
its bytes when an allele is an empty text; it writes nothing of that
block.

Four are an `OSError`: an error of a source that is read, which wraps
`std::io::Error`; a vars file that could not be written, with what went
wrong, which is the `std::io::Error` the system gave when the cause is of
the file system, and what arrow-rs said when it is not; a file that starts
as an arrow file and was cut short, with what was being read when the bytes
ran out; and a batch that arrow-rs could not decode or decompress, with the
batch and what arrow-rs said. The last two are a file that was damaged after
it was written and whose bytes no longer decode, which the reader refuses
instead of giving the variants it can still read. Damage that does decode,
into content the format does not allow, is one of the `ValueError` above
instead: a byte of the `gts` column changed to 254 is a whole batch that
arrow-rs reads and an allele of -2, and popnei's convention sorts an error
by what is wrong with the content and not by what made it wrong, since the
reader cannot tell a damaged file from one another program wrote badly. The
write has a case of its own because a Python user reads
which file went wrong from the exception, and a disc that fills up while the
vars file is being written is not the VCF failing to be read.

Three are a `RuntimeError`: a block that does not fit the writer, with the
individuals and the ploidy of the writer and of the block; a block whose
columns differ from those of the first one written, with the fields of both;
and a chromosome number that the table given with the block has no name for.
A user reaches none of the three by what they write, since `write_vars` asks
its reader for every field and gives the writer the table of that reader.

The extra fields and the writer of segments, added on 9 November 2026, add
these. A `ValueError`: variants out of order, with the two variants; an
extra field of the file that is not what its entry of the `popnei` key
says, with the field, the type or the width found and the one expected,
which is also the case of an entry without its column and of a name the
rule refuses or that is there twice; and a name asked of
`select_extra_fields` that the file has no field of, with the fields it
has. The checks of a segment of "A call that fails" are made in the
Python binding crate before a block is built, and are a `ValueError` of
that crate that names the segment. A `RuntimeError`: a block whose extra
fields are not those `set_extra_fields` declared, with both lists, which
the reader of segments never gives the writer.

Four more cases that a call of this module gives are not its own. A path
that already exists is refused by the binding crate before the core is
called, and is a `ValueError`. A `num_vars_per_block` of 0 is the case of
`docs/specs/block.md` that every reader taking a size gives, which the
writer gives too, a `ValueError`. A block the machine does not give the
memory for is the other case of `docs/specs/block.md` that every reader
gives: the reader asks with `try_reserve` for each column of the block it
builds from a batch, and the size in the message is the one a file fixed,
the third of the three of that spec, whose way out is to write the file
again with a smaller `num_vars_per_block`. It is the error of a file whose
batch this machine cannot count in a `usize` too, which under wasm, where a
`usize` is 32 bits, is a batch of more than 4295 million bytes, and of a
footer that says the file holds more variants than that. A file that could not be opened, which `from_path` gives for a path that
is not there and for a directory, is the case of `docs/specs/io_vcf.md` with
the path and the `std::io::Error`, an `OSError` built with the number the
system gave.

In Python every error of a file names the file. The core does not have the
path, since a writer is built over a sink and a reader over bytes, so it is
the binding crate that puts it there, and where it puts it follows the
exception. The message of a `ValueError` and of a `RuntimeError` starts with
the path. An `OSError` carries it in `filename`, which is where a Python user
of any library looks for it and which Python prints after the message of the
exception, so putting it in the message too would say it twice.

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

None. The owner decided on 20 September 2026 the five points of the first
version of this spec: a format of popnei's own that owes pyNei's nothing,
lz4 in every build, the name `open_vars`, the regions of each batch in the
footer, and every field in the file that `write_vars` makes. On 9 November
2026 they decided the two that came with the extra fields and the writer
of segments: `write_vcf` refuses a source with extra fields, and the
TypeScript functions of the two items wait for a user who needs them. Each
is written where it applies, with the option that was not taken.

## Not in this spec

- The function that asks a `Variants` for the variants of regions, in Python
  and in TypeScript: the filter by regions of `docs/specs/filters.md`, which
  says when this reader is handed the regions and what it then skips.
- The read ahead thread that decompresses the next batch while the consumer
  works: with the first calculation that consumes blocks.
- Genotypes packed in 2 bits, the option to measure of section 4 of the
  architecture: it would be a later version of this format.
- Reading a vars file of pyNei, of either of its formats.
- The reading of a `File` of a page by ranges, which is what `openVars`
  takes in the wasm package: `docs/specs/js_sources.md`.
- The VCF writer: an item of `docs/specs/io_vcf.md`.
- A filter that calls only the genotypes whose genotype field passes a
  threshold, setting the others to missing, which issue 19 names as the
  next step: a step of `docs/specs/filters.md`, when it is asked for.
- Reading the INFO and FORMAT fields of a VCF into extra fields.
- A variant field of a list of values of its own length for each
  variant, which the tree of the alleles of a haplotype region would need.
  The format holds a list column, and nothing asks for one yet.
- `Variants.from_gt_array`, a `Variants` of an array of genotypes in
  memory, which has a spec and no code: `docs/specs/variant.md`.
- `Variants`, the handle both functions here give and take:
  `docs/specs/variant.md`. `iter_blocks`, which is how its genotypes reach
  Python and TypeScript: `docs/specs/block.md`.
