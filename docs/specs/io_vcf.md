# The io::vcf module: the VCF reader

September 2026. The VCF reader is how variants get into popnei: it reads a
VCF, plain or gzipped, and gives its variants in blocks, runs of
consecutive variants held as arrays, with the genotypes as small integers.
This spec develops the reader of the row `io::vcf` of the table in section
9 of `docs/architecture.md`. It depends on `docs/specs/block.md`, which has
the `Block` that the reader gives and the `BlockReader` trait that it
implements, and on `docs/specs/variant.md`, which has the `Needs` that say
which fields a consumer wants and the table of the chromosome names. The
VCF writer of the same row is an item that is not written yet.

There is code, built from the first version of this spec, in which the
reader filled one `Variant` at a time for its caller. The owner dropped the
single variant on 20 September 2026, for the reasons at the end of section
1 of the architecture. What a line of a VCF means, what is refused and the
reference files are as they were built and reviewed; what changes is "How
it runs", how the tests reach the reader, and the speed, which the reader
as built does not reach.

An individual is what `docs/glossary.md` calls one organism that was
genotyped. VCF calls it a sample, and the columns of a VCF after FORMAT
are here the columns of the individuals.

## The VCF reader

### What it gives

From each data line of a VCF, one variant, a row of the block:

| column of the VCF | column of the `Block` |
|---|---|
| CHROM | `chrom`, the number of the name in the reader's table of chromosomes |
| POS | `pos` |
| ID | `id`, empty when the column is `.` |
| REF and ALT | `alleles`, the reference first; only the reference when ALT is `.` |
| QUAL | `qual`, NaN when the column is `.` |
| the GT of each individual | `gts` |
| FILTER | decides whether the variant is given at all |
| INFO and the other values of each individual | not read |

An allele is kept as the text the VCF has, so a symbolic allele, `<DEL>`,
and the allele of an overlapping deletion, `*`, are alleles like any
other.

By default the reader gives only the variants that passed their filters:
those whose FILTER is `PASS`, or `.`, which in a VCF says that no filter
was applied. A variant with anything else there is skipped as if its line
were not in the file. With `only_passed` false every variant is given,
and nothing in the block says which ones had failed. The owner
decided on 20 September 2026 that the FILTER is honoured and that this is
the default; pyNei ignores the column. That `.` counts as passed was
decided here: many programs write `.` in every line, pyNei's own script
for its reference VCF among them, and with `.` as a failure the default
would give such a file no variants. `bcftools view -f .,PASS` is the same
choice.

The genotype of an individual is the value of the key `GT` in its column.
The reader finds where `GT` is among the keys of the FORMAT column of
each line, and an individual that drops its last values, which VCF
allows, still has its `GT`. The alleles of a genotype are separated by
`/` or by `|`, which says that they are phased; the reader takes both and
keeps no phase, because a block has nowhere to hold it. VCF 4.4 lets a
genotype start with a separator, `/0/1` or `|1|1`, and the reader reads
those as `0/1` and `1|1`. An allele written `.` is `MISSING_ALLELE`, -1.
The separator at the start is taken off before anything else is looked at,
so `/.` is the genotype `.`, a missing one; bcftools refuses that line,
and the reader as built reads it, which the reviews of that code left as
it was.

Every allele number of a genotype has to be one of the alleles that REF
and ALT declare: with one alternative allele, 0 and 1. A larger number is
an error of the line. The opposite is fine: an allele that ALT declares
and no genotype carries is what is left when the individuals that carried
it are taken out of a file. The owner decided this on 20 September 2026;
pyNei reads such a number, and the option not taken was to do the same.

The ploidy, how many alleles a genotype holds, is the same for every
individual and every variant, and the caller gives it. It is 1 or more
and at most `MAX_PLOIDY`, 255, which is above the ploidy of anything that
has been sequenced and keeps a ploidy that came from a user from asking
for a genotype of more alleles than a machine can hold; a ploidy outside
that range is an error when the reader is built. A genotype written
as a single `.` is a missing genotype of that ploidy, every allele of it
missing. Any other genotype with a number of alleles that is not the
ploidy is an error, which says that popnei does not read a VCF of mixed
ploidies and gives the line, the individual and the two ploidies. The
owner decided this on 20 September 2026. VCF itself allows the ploidy to
change from one individual to another, a male on the X chromosome is
written with one allele, and such a file is refused: the calculations of
popnei are not defined for it. The options not taken were to take the
ploidy from the first genotype of the file, and to fill a shorter genotype
with missing alleles, as the spike did, the trial parser in Rust that
`docs/rust_core.md` reports and whose code is `spike/pynei_spike` in the
pyNei repository.

### Its Python and TypeScript functions

```python
def open_vcf(
    vcf_path: str | Path, ploidy: int = 2, only_passed: bool = True
) -> Variants
```

It mirrors `vars_from_vcf` of `pynei/io_vcf.py` under another name, which
the owner decided on 20 September 2026: nothing is read when it is
called but the header, and what it returns is not pyNei's `Variants` of
chunks. It reads the header, so a file that is not a VCF fails at the
call and not at the first calculation, and the `Variants` it returns is
the handle of `docs/specs/variant.md`: every pass over it opens the file
again.

In TypeScript, `openVcf(source, {ploidy = 2, onlyPassed = true})`,
with the two options in an object, which is how TypeScript writes
arguments that have a name and a default. `source` is a `Uint8Array` with
the bytes of the file; reading a `File` that a user picked in a page is
under "Not in this spec". A `source` that is not a `Uint8Array`, a
`ploidy` that is not a whole number of 1 or more and at most 4294967295,
which is what a whole number of the core holds in wasm, and an
`onlyPassed` that is not a boolean are an `Error` at the call that says
what was given. A number of JavaScript is a float64 and reaches the core as an
integer of 32 bits, so a ploidy of 2.5 or of 2^32 + 2 would otherwise be
read as 2, and bytes that are not a `Uint8Array` would be read as the
bytes of something the user never wrote.

The differences from pyNei:

- The name, `open_vcf` for `vars_from_vcf`.
- The `ploidy` argument is new, and the rule above with it. pyNei takes
  the ploidy from the first genotype, and gets 2 for every file (see
  "What pyNei does that is odd").
- The `only_passed` argument is new, and by default the variants that
  failed a filter are not given. pyNei gives them all.
- `desired_num_vars_per_chunk` is gone. The size of the blocks is
  `num_vars_per_block` of `iter_blocks`, of `docs/specs/block.md`.
- An allele number that REF and ALT do not declare is an error. pyNei
  reads it.
- A VCF with a header and no variants gives no variants, like a file
  that has been read to its end. In pyNei it is an error. The owner
  decided this on 20 September 2026.
- A line with fewer columns of individuals than the header is an error.
  pyNei reads the individuals that are not there as homozygous for the
  reference. A line with more is an error too. pyNei fails on it with an
  `IndexError`, unless the genotype that is left over is homozygous for
  the reference, which it skips without looking where it goes.
- Two individuals with the same name are an error when `open_vcf` is
  called. In pyNei the call returns and the error comes with the first
  chunk.
- An allele number above 127 is a `ValueError`, like the other errors of
  a data line. In pyNei it is a `NotImplementedError`.
- A genotype that starts with a separator is read. pyNei refuses it.
- An empty line is skipped. pyNei fails on it with an `IndexError`.
- A quality that is not finite is an error. pyNei reads `nan` as a variant
  with no quality, and `inf` and `1e400` as an infinite one.
- A bgzipped file with no mark of its end is an error. pyNei reads
  `many.vcf.gz` without its last 28 bytes as it reads the whole file, 500
  variants and no complaint; bcftools 1.24 stops with "no BGZF EOF marker;
  file may be truncated".

### The cases a reader of the rules would not guess

The file is gzipped when its first two bytes are `1f 8b`, whatever its
name, as in pyNei. A file made by bgzip, which is what nearly every
gzipped VCF is, is many gzip members one after another, each with 64 KB
of text at most, and an empty one at the end. A gzip decoder that stops
after the first member gives the start of the file and no error. For
`many.vcf.gz` and `cases.vcf.gz`, below, the first member is exactly the
header, so such a decoder gives a VCF with a header and no variants,
which is not an error either: a reader with no variants and nothing to
tell that anything went wrong. What catches it is the tests that count
the variants of the gzipped files. So the reader uses a decoder that goes on to the next member, flate2's
`MultiGzDecoder`, with flate2's default backend, `miniz_oxide`, which is
Rust; its zlib backends are C and do not build for the wasm package.
After the decompression, or with no compression, the first byte has to be
`#`; if not, the error says that the source is not a VCF.

A file made by bgzip ends with an empty block of 28 bytes that marks its
end. The reader knows that a source was made by bgzip from its first gzip
member, which carries the extra field `BC` that bgzip writes in every
block, and such a source that does not end with the empty block is an
error. The source is read once and forward, so the reader watches the
compressed bytes as they go by to the decoder and keeps the last 28, and
that the mark is missing is known only when the source ends: every
variant is given first, and the error comes where the reader would have
said that there are no more, which is when bcftools says it too.
Without it, a file that was cut where one gzip member ends and the next
begins, a download that stopped, gives fewer variants and no error, where
bcftools says "no BGZF EOF marker". The owner decided this on 20 September
2026. A gzip file that bgzip did not make has no such mark and is read to
its end.

The header is every line that starts with `##`, which is skipped, and
then the line that starts with `#CHROM`, whose first nine columns have to
be the nine of a VCF with genotypes and whose other columns are the
individuals. A VCF with no FORMAT column is an error, as in pyNei, and so
is one with a FORMAT column and no individual. Two individuals with the
same name are an error; pyNei raises it in `Genotypes.__init__` of
`pynei/variants.py`, when it builds the first chunk. An individual with
no name is an error too, which is what a `#CHROM` line that ends in a tab
has.

A line ends in `\n` or in `\r\n`. The genotype of the last individual is
the one that would carry the `\r`, so it is taken off the line first;
pyNei's `test_vcf_with_windows_line_ends` and
`test_missing_gt_in_the_last_sample` are about that.

An allele number above 127 is an error even when ALT declares that many
alleles, as in pyNei, whose `test_an_allele_over_the_limit_is_refused`
asserts a `NotImplementedError`. An allele number is a run of digits, so
`+1` is not one. A FORMAT with no `GT`, a position that is not a number,
a quality that is neither a number nor a dot, and a line with fewer than
ten columns are errors. So is a quality that is a number and not a finite
one, `nan`, `inf` or `1e400`, which a float reads as infinite: it is an
error of the QUAL column, because NaN is what a block holds for a variant
with no quality. The owner decided this on 20 September 2026; the option
not taken was to keep what the float gave. A line whose bytes are not text, not valid UTF-8,
is an error of that line and not an error of the input: a VCF is text.
Every error of a data line gives the number of the line in the file,
counted from 1 with the header lines, and the column or the individual.

A line that is skipped for its FILTER is not parsed beyond that column,
so what is wrong in the rest of it is not found. The seven columns up to
the FILTER have to be there, since the FILTER is the seventh, and what
they hold is read only when the variant is given, so what is wrong
inside them is not found either: a position that is not a number in a
line that is skipped gives no error, and the name of the chromosome of
such a line gets no number.

The alleles of ALT are counted for every variant that is given, to check the allele
numbers of the genotypes, also when `ALLELES` was not asked for and the
texts of the alleles are not kept. It is a pass over a column of a few
bytes. An allele with no letter in it, which is what an empty REF, an
empty ALT or a trailing comma in ALT gives, is an error of its column;
bcftools reads no alternative allele in `T,`.

The chromosomes get their numbers in the order in which they first appear
in the data lines that are given. The `##contig` lines of the header are
not used, since a VCF does not have to have them.

### What pyNei does that is odd

pyNei was run at commit ef0ca6e, each case a VCF of one variant and three
individuals read with `vars_from_vcf`.

It reports a ploidy of 2 for every file. `_parse_var_line` takes the
ploidy as the length of what `_parse_gt` returns, which is a pair, whether
the genotype is phased and its alleles. A tetraploid file, `0/0/1/1`,
`0/1/1/1`, `0/0/0/0`, comes out as `0/0`, `1/1`, `1/1` with no error, a
haploid `1` as `1/0`, and a genotype written `.` in a diploid file as
`./0`. It is issue 19 of pyNei. popnei does not reproduce it, and until
pyNei is fixed the two are compared on diploid files only.

It accepts an allele number that the variant does not declare. Its test
VCF has the genotypes `3|4` and `5/6` in a variant with one alternative
allele, and `test_vcf_parser` asserts them. bcftools 1.24 prints them
back with no complaint. popnei refuses them, so pyNei's test VCF is not
among the files the two are compared on.

It holds a position in 32 bits. `_parse_vcf_vars_chunk` of
`pynei/io_vcf.py` builds the position column with `config.PANDAS_INT_DTYPE`,
which is pandas' `Int32Dtype`; `config.py` has a `PANDAS_POS_DTYPE` of 64
bits that only a test uses. On a VCF of one variant, a position of 2147483647 is read and
one of 2147483648 raises `TypeError: cannot safely cast non-equivalent int64
to int32`. popnei holds a `u64` and does not reproduce it, and the two are
compared on positions below that number, which those of the reference VCFs
are.

### How it runs

As a `BlockReader`. For each block it reads lines from the source until it
has as many variants to give as the block takes, `num_vars_per_block`, and
parses them into the arrays of that block, each line into its own row.
Natively the lines are parsed on the threads of rayon, and since no two
lines write the same row nothing is shared between them; in wasm the same
code parses one line after another.

Which row a line gets is known before it is parsed. A line that is skipped
for its FILTER and an empty line have no row, so while the lines are cut
from the source a serial pass finds the FILTER of each, the text between
its sixth and its seventh tab, and numbers the ones that will be given.
The pass gives no error. A line with fewer than seven columns has no
FILTER to find, so it gets a row, and the parse gives the error of a line
with too few columns, the same one whatever `only_passed` is. With
`only_passed` false the pass is made too, and all it takes out is the
empty lines.

The text of the lines of a block is not read at once. A block is about 5
million genotypes whatever the individuals are, at 4 bytes of text each
about 20 MB, so the reader takes the lines in batches, bounded by how many
lines a batch holds and by how many bytes of text they are, and parses
each batch into the rows of the block that follow the ones already filled.
A batch ends where its block does, so with small blocks the block is what
bounds it.
The bound in bytes is what keeps the memory of a reader from growing with
the individuals of the file, since a line carries one genotype per
individual, and a batch holds one line at least, however long that line
is, so the reader always goes forward. Both bounds are
constants of the code, each with what was measured on it, and a caller
that times the reader can set them. The reader as built has them at 1024
lines and 8 MiB.

The genotypes, the positions and the qualities of a line go straight into
its row. The texts do not, because the rows of a column of texts are not
of one size: the id and the alleles of a line are parsed into buffers of
that line, and appended to the columns of the block in order, serially,
after each batch, and only when they were asked for. The chromosomes get
their numbers then too, in the order of the variants that are given, so
the numbers do not depend on the threads.

An error loses its block, as section 1 of the architecture has it: the
blocks before the one with the wrong line are given, then the error comes
in the place of that block, and after it the reader gives no more. When
two lines are wrong the error is that of the one that comes first in the
file: the batches are parsed one after another, and of the lines of a
batch that failed the reader takes the first.

A panic inside the parse of a batch leaves the reader with lines that
were never parsed, so a reader whose parse did not come back gives an
error at its next call and no more blocks. Nothing a VCF can hold panics
the parse; what this is for is that a panic of a defect of popnei becomes
an exception that the caller may catch, and a reader that went on
afterwards would give the variants of the lines that were parsed and drop
the others without a word.

The result does not depend on the number of threads, nor on the bounds of
a batch.

Only what `Needs` asks for is parsed, and a block has the columns that
were asked for and no other, the chromosomes and the positions among them.
A column that is not parsed is not checked: a position that is not a
number is an error only with the chromosome and the position asked for,
and a quality of `nan` only with the quality asked for. The reader as
built parses the position of every line, which was the rule of the first
version of this spec.
What is checked whatever is asked for is the shape of the line: the nine
first columns have to be there, the FORMAT has to have a `GT` key, and
there has to be one column after the FORMAT at least. With the genotypes
not asked for, the reader does not look at the columns of the individuals:
how many of them there are and what is in them is not read, so a line with
two columns of individuals under a header with three, or one with a
genotype of another ploidy, is given, and `gts` is empty.

The reader is built over any `BufRead`, as section 1 of the architecture
asks. It reads the first two bytes of the source to find the gzip and
hands them back in front of it, because one look at the buffer of a
source may give fewer than two bytes: a pipe, or the bytes of a file that
a page hands over a few at a time. Nothing of the source is consumed, and
the bytes after the first two are read from the buffer of the source
itself. A function that takes a path opens the file and does the same,
for the callers that have one; a file that cannot be opened is an error
that carries the path.

### How it is verified

Against bcftools 1.24, which is on the owner's machine, run by
`tests/reference/vcf/make_reference.py` on three VCFs that the script
writes, of three or fifty diploid individuals:

    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n' cases.vcf

and its output is kept beside each file, `cases.bcftools.tsv`. bcftools
prints the genotypes as text, with their separators, and the tests turn
`0|1` into 0, 1 and `.` into -1. The comparison is exact. The rows whose
FILTER is `PASS` or `.` are what the reader has to give by default, and
all the rows what it gives with `only_passed` false; `bcftools view -f
.,PASS`, which keeps those same rows, left 3 of the 4 variants of
`cases.vcf` and 475 of the 500 of `many.vcf`.

`cases.vcf`, four variants written by hand, which pyNei reads in the same
way, and `cases.vcf.gz`. The literals of the first cargo tests are its
four rows with `only_passed` false, with `gts` as the reader gives them,
and the first, the third and the fourth with the default:

| chrom | pos | id | alleles | qual | FILTER | the VCF has | `gts` |
|---|---|---|---|---|---|---|---|
| chr1 | 100 | rs1 | A, T | 29.5 | PASS | `0/0 0/1 1/1` | 0, 0, 0, 1, 1, 1 |
| chr1 | 200 | empty | A, T | none | q10 | `./. 0\|1 .\|0` | -1, -1, 0, 1, -1, 0 |
| chr1 | 300 | empty | A, G, T | 67 | PASS | `1/2 2\|1 2/2` | 1, 2, 2, 1, 2, 2 |
| chr1 | 400 | empty | T | 47 | PASS | `0/0 0/0 0/0` | 0, 0, 0, 0, 0, 0 |

The second variant has a half called genotype, `.|0`, one with some
alleles called and some missing. The chromosome number is 0 for the four.

`differences.vcf`, two variants that pyNei does not read as bcftools
does, both with `PASS`. bcftools prints a leading separator away and
keeps the single `.`:

| chrom | pos | id | alleles | qual | the VCF has | bcftools prints | `gts` |
|---|---|---|---|---|---|---|---|
| chr2 | 50 | ms1 | GTC, G, GTCT | 50 | `0/1:3 0/2 .` | `0/1 0/2 .` | 0, 1, 0, 2, -1, -1 |
| chr2 | 60 | empty | A, `<DEL>`, `*` | none | `/0/1 \|2\|2 0/0` | `0/1 2\|2 0/0` | 0, 1, 2, 2, 0, 0 |

The chromosome number of both is 0, since `chr2` is the first name of
that file.

`many.vcf`, 500 variants of 50 individuals drawn with a fixed seed,
117 KB, and `many.vcf.gz`, which bgzip wrote as 4 gzip members with 617,
65252, 51477 and 0 bytes of text, so a decoder that stops at the first
one fails these tests. 450 of its variants have `PASS`, 25 have `.` and
25 have `q10`; the first with `q10` is chr1 1259. 333 of them have an
id and 400 a quality, 200 of those with a decimal, so that the
comparison with pyNei covers the id and the quality of 500 variants and
not of the 4 of `cases.vcf` alone; the id, the quality and the FILTER of
a variant follow its place in the file and no draw of the generator, so
that they can be given to it without a genotype moving. Every genotype,
chromosome and position is compared with `many.bcftools.tsv`, and these
counts, worked out from that file, are literals:

| | every variant | by default |
|---|---|---|
| variants | 500 | 475 |
| of them in `chr2` | 250 | 238 |
| with two alternative alleles | 54 | 53 |
| missing genotypes | 1511 | 1431 |
| of them half called | 257 | 240 |
| missing alleles | 2765 | 2622 |
| called alleles | 47235 | 44878 |
| the sum of the allele numbers of the called alleles | 25954 | 24831 |

The first five genotypes of the first variant, chr1 1000, are `1/1 .|.
1/0 1/1 0/1`. The plain and the gzipped file give the same.

Against pyNei: Python sees what the reader read only through
`iter_blocks`, so the comparison of every field of every variant of
`cases.vcf`, `cases.vcf.gz`, `many.vcf` and `many.vcf.gz` with pyNei's
chunks is the pytest test of "How it is verified" of
`docs/specs/block.md`, with `only_passed=False` because pyNei gives
every variant. The same test compares `individuals`, `num_individuals`
and `ploidy` with pyNei's `samples`, `num_samples` and `ploidy`. pyNei
read `many.vcf` and `many.vcf.gz` as bcftools did when this was written,
and it refuses `differences.vcf` at its leading separator.

The errors, each a cargo test on a VCF written in the test, which checks
the kind of the error and the line and the individual it names: a
tetraploid genotype with `ploidy` 2; a haploid genotype with `ploidy` 2;
the same tetraploid file read with `ploidy` 4, which is not an error and
gives `0, 0, 1, 1` for `0/0/1/1`; an allele 2 in a variant with one
alternative allele, asked for with `GTS` alone; an allele of 128; a line
with two columns of individuals under a header with three, and one with
four; a FORMAT with no `GT`; a position `x`; a source that starts with
neither `#` nor the gzip bytes; a header with no FORMAT column, and one
with it and no individual; two individuals with the same name; an
individual with no name, a `#CHROM` line that ends in a tab; a ploidy of
0 and one of 256; a quality `x`; an ALT that ends in a comma and an empty
REF; a line whose bytes are not valid UTF-8; and a path that no file is
at, whose error carries the path. With `ID` and `ALLELES` asked for and
no `GTS`, these four: a line of seven columns and a FORMAT of `DP`, which
are errors, and a tetraploid genotype under a ploidy of 2 and a line with
the columns of two individuals under a header with three, which are read.
And these,
which are not errors: `GT` second in the FORMAT, `DP:GT` with `3:0/1`;
lines that end in `\r\n`; an empty line at the end; a header and no
variant, which gives no block at the first `next_block`; a variant whose
ALT declares two alleles and whose genotypes carry only the first; a
line with `q10` and a tetraploid genotype, which the default skips and
`only_passed` false refuses; a line with `q10` and the position `x`,
which the default skips; a last line with no end of line; and a source
that gives one byte at a time, gzipped and plain, which a reader that
looked for the two bytes of gzip in one look at the buffer would refuse.

What `Needs` does: with `GTS` alone, a block has the genotypes and no
column; with `ID` and `ALLELES` and no `GTS`, those two columns and an
empty `gts`; and a reader asked for everything and then for `GTS` alone
gives its next block without the columns. That the allocations of a pass
over `many.vcf` with `GTS` alone asked for are a few for each block and
none for each variant is checked once by hand with a counting allocator
when the reader is written, and it is not a test that stays. With the ids
asked for there is one for each variant, since the ids of a block are a
`String` each.

The blocks: `many.vcf` with blocks of 100 and the default options gives
five blocks, of 100, 100, 100, 100 and 75 variants, and with every variant
given five of 100; with blocks of 1000, one of 475. What the blocks hold,
joined, is the same with blocks of 1, 7, 100 and 1000, with batches of 1,
3 and 1024 lines, and on 1 thread and on 4. A VCF written in the test with
a tetraploid genotype in its third variant, read in blocks of 2, gives one
block and then the error, and then no block; with a second wrong line
after it, the error is still that of the third variant. A source that is
the bytes of `many.vcf.gz` without its last 28 gives its five blocks and
then the error that names the mark of the end of a bgzipped file, and one
cut after its second gzip member gives the blocks of its 280 variants and
then that error; `many.vcf` compressed with gzip and not with bgzip is read.
A quality of `nan`, of `inf` and of `1e400` is an error of the QUAL
column.

The cargo tests are made at `VcfReader::new` for what is wrong in the
header, the source that is not a VCF, the FORMAT column or the
individuals that are not there, the repeated name, the name that is not
there and the ploidy out of range, at `from_path` for the file that is
not there, and at
`next_block` for the rest, with the reader built over the bytes of the
file. The pytest tests are made at `open_vcf` and the blocks of what it
returns: the counts of the table above on `many.vcf`, with the default
and with `only_passed=False`; the two rows of `differences.vcf`; a file
that is not a VCF, which is a `ValueError` at the call; and the
tetraploid file, a `ValueError` when the blocks are asked for. The
TypeScript test, under node, reads `cases.vcf` and `differences.vcf` with
`openVcf` from a `Uint8Array` and compares their blocks with the two
tables above, with the default and with `onlyPassed` false.

## The Rust interface

What the caller says about the file. The default is a ploidy of 2 and
only the variants that passed.

```rust
/// The ploidy and the FILTER of the default, each a constant so that the
/// binding crates and the tests name them.
pub const DEFAULT_PLOIDY: usize = 2;
pub const DEFAULT_ONLY_PASSED: bool = true;
/// The largest ploidy a reader takes.
pub const MAX_PLOIDY: usize = 255;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VcfOptions {
    /// How many alleles every genotype has. 1 or more.
    pub ploidy: usize,
    /// Skip the variants whose FILTER is neither PASS nor a dot.
    pub only_passed: bool,
    /// How many variants a block has, 1 or more, or None for
    /// `default_num_vars_per_block` for the individuals of the header.
    pub num_vars_per_block: Option<usize>,
}
impl Default for VcfOptions { /* 2, true, None */ }
```

Where in a data line something is wrong, which the error of a data line
carries beside the number of the line. The column is one of the nine
fixed names, the individual is the name the header gave it, and the line
is for what no one column is at fault for, the count of its columns.

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum VcfPlace {
    Column(&'static str),
    Individual(String),
    Line,
}
```

The reader. `new` reads the header, so the individuals are known when it
returns, and it fails when the source is not a VCF with genotypes, the
ploidy is out of range, or the size of the blocks is one that
`docs/specs/block.md` refuses.

```rust
pub struct VcfReader<R: BufRead + Send> { /* private */ }
impl<R: BufRead + Send> fmt::Debug for VcfReader<R> { /* not R: Debug */ }

impl<R: BufRead + Send> VcfReader<R> {
    /// `source` is the VCF, gzipped or not.
    pub fn new(source: R, options: VcfOptions) -> Result<VcfReader<R>>;
}

impl VcfReader<BufReader<File>> {
    pub fn from_path(path: &Path, options: VcfOptions) -> Result<Self>;
}

/// The two bounds of a batch of lines. Hidden from the documentation and
/// outside what popnei promises: they are for the benchmark, which times a
/// file with one bound after another, and for the tests, which read a file
/// of a few hundred lines in several batches and one line at a time, the
/// batch of wasm. What a reader gives does not depend on them.
impl<R: BufRead + Send> VcfReader<R> {
    pub fn set_lines_per_batch(&mut self, lines: usize);
    pub fn set_bytes_per_batch(&mut self, bytes: usize);
}

impl<R: BufRead + Send> BlockReader for VcfReader<R> { /* ... */ }
```

The cases this module adds to the error of the crate: the source is not a
VCF, with what was found; a wrong header, with what is wrong; a ploidy
out of range, which is the one thing `new` refuses that is not in the
source, with the ploidy that was asked for; a wrong data line, with the
number of the line, the column or the individual, and what is wrong; a
genotype of another ploidy, with the line, the individual, the ploidy of
the genotype and the one expected; a file that could not be opened, with
its path and the `std::io::Error` as the source of the error, so that a
binding can put the path where the language of the binding keeps it,
`OSError.filename` in Python; an error of the input, which wraps
`std::io::Error`; a bgzipped source with no mark of its end; and a parse
of a batch that did not come back, with the number of the last line that
was read. In Python the first five and the last two are a `ValueError`
and the other two an `OSError`.

## Speed

The number to reach is that of the spike, the trial parser in Rust of
section 3 of `docs/rust_core.md`, which parses a chunk of lines with rayon
and writes each row straight into the array of the chunk, as this reader
now does. The reader as built, which filled one `Variant` at a time, does
not reach it. Measured on 20 September 2026 on the owner's Apple M5 Pro,
18 cores, release, the file in the page cache, `GTS` asked for and the
default options, the median of 5 runs of
`crates/popnei/benches/read_vcf.rs`, on a VCF of 100000 variants and 1000
individuals, 403 MB plain and 38 MB bgzipped, which
`crates/popnei/benches/make_big_vcf.py` makes with `simulate_genotypes` and
`write_vcf` of pyNei's `test/gwas_reference/make_reference.py`, a seed of
42, 3 in 100 genotypes missing and `.` in every FILTER, so the default
gives every variant:

| | the reader as built | the spike, same file, same day |
|---|---|---|
| plain, 1 thread | 1.24 s | 0.54 s |
| plain, 18 threads | 0.160 s | 0.098 s |
| bgzipped, 1 thread | 1.58 s | 0.84 s |
| bgzipped, 18 threads | 0.50 s | 0.40 s |

plink2 v2.0.0-a.7.7 reads the plain file in 0.273 s on one thread, and
pyNei in 13.5 s. The spike does three things less than the reader: it
checks neither the ploidy of each genotype nor its allele numbers against
ALT, and it does not read the FILTER. The number to reach is the spike's
of that table, or no more than a tenth above it, which is the rule that
the first version of this spec set, with the numbers that the spike gave
then, 0.55 s and 0.11 s. `docs/rust_core.md` has a row for that VCF
gzipped into 53 MB and does not say how it was compressed; bgzip makes 38
MB of it, the session that built the reader could not make the file of 53
MB again, and the bgzipped rows here stand in its place.

Where the time of the reader as built goes, from a sampling profile of the
run on one thread: 95 in 100 in the parse, almost all of it in the columns
of the individuals, of which filling the genotypes is 46 in 100 of the
self time, splitting a text at a character 25, searching a byte 16 and
comparing strings 10. The nine first columns, the FILTER and the count of
the alleles are under 1 in 100. So the hand out of the variants was not
what cost, and rows written into a block will not close the gap alone.
What the spike does and the reader as built does not is to parse the bytes
of the columns of the individuals with `memchr`, with no text and no
UTF-8 check on them; and a row of a fixed length, `num_individuals` x
`ploidy`, makes the check of the ploidy of a genotype a check of length.
That is what the implementer tries first, and measures, before any other
work on speed. With threads the serial reading of the lines is the floor.
`docs/reports/vcf-to-blocks.md` has the measurement: on a file of 5000
variants of 1000 individuals the reading of the lines alone took 10.8 ms
on one thread and 4.1 ms on eight, the whole parse 92.5 ms and 13.8 ms,
and 4.1 + (92.5 - 10.8) / 8 = 14.3 ms predicts the 13.8. Batches of 256,
1024 and 4096 lines, with the bound of 8 MiB, took 0.202, 0.158 and
0.151 s on the file of the table on 18 threads.

## Open points

None. The owner decided on 20 September 2026 the six that there were:
the ploidy as an argument and the refusal of mixed ploidies, the variants
that failed a filter left out by default, an allele number that is not
declared as an error, no error for a VCF with no variants, a quality that
is not finite as an error, and a bgzipped source with no mark of its end
as an error. Each is written where it applies, with the option that was
not taken when there was one.

## Not in this spec

- The VCF writer: a later item of this spec. pyNei has none.
- BCF, the binary form of VCF, and reading a region of an indexed file:
  not planned.
- The values of INFO and of the individuals other than `GT`: popnei has
  no use for them.
- Choosing variants by the name of the filter they failed: `only_passed`
  is all there is.
- `Variants`, the handle that `open_vcf` returns: `docs/specs/variant.md`.
  `iter_blocks`, which gives its genotypes: `docs/specs/block.md`.
- The reader of the vars file, popnei's own file of variants:
  `docs/specs/io_vars.md`.
- A web worker and the reading of a `File` by ranges in TypeScript:
  section 11 of the architecture, and the plan that builds the
  TypeScript side.
