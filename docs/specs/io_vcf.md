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
would give such a file no variants.

A FILTER that names more than one filter is where popnei and bcftools
differ. popnei reads the whole column: it gives the variant when the column
is `PASS` or `.` and no other. `bcftools view -f .,PASS` keeps a row when
any of the filters the column names is one of those, so it keeps
`PASS;q10`, a variant that failed `q10`, and popnei skips it. On the seven
columns `PASS`, `.`, `q10`, `PASS;q10`, `pass`, `q10;PASS` and an empty
one, popnei gives the first two and bcftools those two and the two that
name `PASS` beside `q10`. Which is right depends on what the file means by
naming both, which VCF does not say; popnei takes the strict reading, and a
user who wants the other one has `only_passed` false and a filter of their
own.

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

The dot is the only way a genotype says the missing allele, and no line
gives an allele below it. An allele number is a run of digits, so the
minus sign is not part of one: `-2/0`, `0/-1`, `-1` and `-0` are each a
wrong data line, which names the line and the individual and says that
the text is not an allele number. The reader was built this way and this
paragraph writes down what it does, so that the rule the owner gave on 21
September 2026, that an allele below `MISSING_ALLELE`, -1, is never
allowed and is refused by every reader of popnei, is one the VCF reader
is held to by a test of each of those four forms.

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
- A bgzipped file that is corrupted is an error that names the member.
  pyNei was run on 21 September 2026 on `many.vcf.gz` with its bytes 320
  and 321 changed, of the paragraph on bgzip below, and on `cases.vcf.gz`
  with its byte 243 set to 144: it raises `ValueError: Empty VCF file, it
  has no variants` for both, which is what it says of a file that has a
  header and nothing else, and not that the file is corrupted.

### The cases a reader of the rules would not guess

The file is gzipped when its first two bytes are `1f 8b`, whatever its
name, as in pyNei. A file made by bgzip, which is what nearly every
gzipped VCF is, is many gzip members one after another, each with 64 KB
of text at most, and an empty one at the end. Every member of such a file
carries in its header the extra field `BC`, two bytes that hold the size
of that member in the file, and the `BC` of the first member is what says
that bgzip wrote the source. A gzip decoder that stops
after the first member gives the start of the file and no error. For
`many.vcf.gz` and `cases.vcf.gz`, below, the first member is exactly the
header, so such a decoder gives a VCF with a header and no variants,
which is not an error either: a reader with no variants and nothing to
tell that anything went wrong. What catches it is the tests that count
the variants of the gzipped files.
After the decompression, or with no compression, the first byte has to be
`#`; if not, the error says that the source is not a VCF.

A source that bgzip wrote is read one member at a time, by the size that
the `BC` of that member states, which is how htslib and bcftools read one:
the reader cuts the member from the source by that size, decompresses its
deflate data on its own, and checks what came out against the CRC32 and
the length of the text that the last eight bytes of the member hold. These
are errors, each of them with the member, counted from 1, the byte of the
compressed source where that member starts, and what is wrong: a member
that does not start with the two bytes of gzip; one whose method is not
deflate, the method 8; one whose flags are not exactly the flag of an
extra field, 4; an extra field whose subfields do not end where it ends,
or that holds no `BC` of two bytes; a size that leaves no room for the
header, the deflate data and the last eight bytes; deflate data that the
decoder refuses, that does not end where the member does or that gives
more text than a member holds; and a CRC32 or a length of the text that is
not the one of the text that came out.

The flags and the method are the two of those that bcftools reads and
popnei refuses: bcftools 1.24 reads a member whose flags are `0c`, the
flag of an extra field with those of a name and a comment. BGZF fixes the
flags of a member to 4, and the reader cuts a member by a size and does not
follow its bytes one by one, so a name or a comment in a header would move
the data of that member to a place the reader does not look at.

The extra field of a member can hold other subfields beside `BC`, before it
or after it, and the reader walks them to its end. It walks the extra field
of the first member too, the one whose `BC` says that bgzip wrote the
source: a file whose first member carries another subfield before its `BC`
is read by the sizes of its members like any other, where a reader that
looked for the `BC` at the bytes 12 and 13, which is where bgzip writes it
and where htslib looks for it, would read it as a plain gzip file and check
none of its members.

A gzip file that bgzip did not
write, whose first member has no `BC`, is read with a decoder that goes on
to the next member by itself, flate2's `MultiGzDecoder`. Both ways
decompress with flate2's default backend, `miniz_oxide`, which is Rust;
its zlib backends are C and do not build for the wasm package. The reader
of the members takes flate2's raw deflate and its CRC32, of that same
backend, so it adds nothing to what popnei depends on.

Why the members are cut and checked, and not handed to a decoder that goes
from one to the next on its own: with such a decoder, `many.vcf.gz` with
the two bytes that hold the length of the extra field of its second
member, the bytes 320 and 321 counted from 0, changed from `06 00` to `44
54`, gives no variant and no error, because the decoder takes 21572 bytes
of compressed data for an extra field and lands on the empty member that
ends the file. A review of 21 September 2026 changed every byte of
`cases.vcf.gz` in turn to each of the 255 other values, 101745 files, and
one of them was read as a whole file with its variants missing and nothing
to say so: the byte 243, the length of that same field of its second
member, set to 144. The owner decided that day, with that review in front
of him, that an error never passes silently and that a corrupted file is
refused however improbable the corruption. bcftools 1.24 gives both of
those files a header, no variant and the exit status 0: htslib takes a
member whose header it cannot read for the end of the data, and the 28
bytes that mark the end of the file are still where they were.

What no reader of a BGZF file can see is a member that was removed,
repeated or moved: every member is a whole gzip stream, checked by its own
CRC32, and none of them records which member of the file it is. A reviewer
made those files from `many.vcf.gz` on 21 September 2026 and read them with
popnei and with bcftools 1.24: with its second member dropped it gives the
220 variants of its third, with that member written twice 780 variants,
with its second and third members swapped the same 500 variants in another
order, and with the mark of the end written twice 500 variants, all of them
with no error in either program. So "an error never passes silently" is the rule for the bytes of a
member and for the end of the file, and the order and the number of the
members are outside what the format lets anybody check.

A file made by bgzip ends with a member that holds no text, the empty
block of 28 bytes, and a source that bgzip wrote and that does not end
with one is an error. A member with no text in the middle of a file is not
its end: bgzip writes one where a caller asked for the bytes so far, and
what makes the last member the mark of the end is that the source has no
more bytes after it. The source is read once and forward,
so that the mark is missing is known only when the source ends: every
variant is given first, and the error comes where the reader would have
said that there are no more, which is when bcftools says it too.
Without it, a file that was cut where one gzip member ends and the next
begins, a download that stopped, gives fewer variants and no error, where
bcftools says "no BGZF EOF marker". The owner decided this on 20 September
2026. A gzip file that bgzip did not make has no such mark and is read to
its end.

A bgzipped source that ends early is that error wherever it was cut, and
the reader gives the variants it read before the cut first in every case;
what a user of `iter_blocks` gets of them is below. Where the cut falls
decides how the reader learns of it. A cut where a
member ends leaves whole members that nothing is wrong with, and what says
that the file is cut short is the mark that is not at its end. A cut
inside the header of a member, or inside its data, leaves the reader with
fewer bytes than that member says it has: it decompresses what there is of
the deflate data, gives the lines that came whole out of it, and then the
same error, so a user whose download stopped is told that the
file is cut short and not that a deflate stream is incomplete. Such a
member has no CRC32 and no length of its text to be checked against, since
those are among the bytes that are missing. A cut inside the first member
is the same error, and not the error of a VCF with no `#CHROM` line, which
is what the text that came out of that member has: the reader knows that
the bytes ran out. Bytes after the member that holds no text are another
thing: the file did not end where it says it ends, so a source with bytes
after that member is a corrupted file and not one that was cut short. An error of
the file system, a disc that fails while the file is read, is not one of
those: it stays the error of the input it is, with the blocks that were
read before it given first.

A user reads through `iter_blocks`, which puts a `reblock` at the end of
the pass, so the variants that `reblock` was keeping for its next block
when the error came are lost with it, as `docs/specs/block.md` says of
every error of a reader. `many.vcf.gz` without its last 28 bytes gives
500, 497, 500 and 0 of its 500 variants before the error with blocks of 1,
7, 100 and the size popnei chooses: 500 variants in blocks of 7 are 71
whole blocks and 3 variants that were waiting, and the size popnei chooses
for the 50 individuals of that file, 10000 variants, leaves the whole file
waiting in one block that was never full. The four were measured from
Python on 21 September 2026.

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
not taken was to keep what the float gave. A line whose bytes are not text,
not valid UTF-8, is an error of that line and not an error of the input: a
VCF is text.

The bytes of a line are read as text where its text is kept, which is its
nine first columns, and the error above is theirs. The columns of the
individuals are read as bytes and never as text, which is what "Speed"
below asks for, so a byte that is not text in one of them is a byte that is
not a digit where an allele number is: the error names that individual and
not the line. With the genotypes not asked for, those columns are not
looked at at all and such a line is read.
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
that times the reader can set them. The reader has them at 4096 lines and
16 MiB, which "Speed" says what was measured on.

The genotypes of a line go straight into its row of the block, which is
where nearly all the bytes of a VCF with genotypes end up. Everything else
is parsed into buffers of that line and appended to the columns of the
block in order, serially, after each batch: the position, the quality, the
id and the alleles, each only when it was asked for, and the number of the
chromosome, which is given then so that the numbers follow the order of the
variants and not the order in which the lines were parsed. What that
serial pass costs was measured on 21 September 2026 on the owner's Apple M5
Pro, release, on the 403 MB VCF of "Speed" below: asking for every column
instead of the genotypes alone adds 1 ms of 127 on 18 threads and 38 ms of
650 on one, the ids being the one column that allocates for each variant.

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
and a quality of `nan` only with the quality asked for.
What is checked whatever is asked for is the shape of the line: the nine
first columns have to be there, the FORMAT has to have a `GT` key, and
there has to be one column after the FORMAT at least. With the genotypes
not asked for, the reader does not look at the columns of the individuals:
how many of them there are and what is in them is not read, so a line with
two columns of individuals under a header with three, or one with a
genotype of another ploidy, is given, and `gts` is empty.

The reader is built over any `BufRead`, as section 1 of the architecture
asks. It reads the first sixteen bytes of the source and hands them back
in front of it: the two of gzip, the ones the message of a source that is
not a VCF shows, and the bytes 12 and 13, where a file that bgzip wrote
names its extra field `BC`. One look at the buffer of a source may give
fewer bytes than that, one even: a pipe, or the bytes of a file that a page
hands over a few at a time. Nothing of the source is consumed, and the
bytes after the sixteenth are read from the buffer of the source itself. A function that takes a path opens the file and does the same,
for the callers that have one; a file that cannot be opened is an error
that carries the path.

A source that bgzip wrote is read through the reader of its members, which
holds one member at a time: the bytes of that member as the file has them
and the text that came out of them, 64 KiB each, which is what it asks of
the machine when it is built, and the extra field of the header of that
member, which is 6 bytes in what bgzip writes and 64 KiB at most. The
memory of a reader does not grow with the file. The bytes of a member are
copied once out of the buffer of the source, which the decoder that goes
from member to member does not do, and every line is then copied out of the
text of the member into the text of its batch, as the lines of a plain file
are. A review of 21 September 2026 put the two together at 7 ms of a read
of 0.99 s of the 38 MB file of "Speed", from a sampling profile. Cutting a member
from the source and decompressing it are two steps, which is what a later
plan that decompresses the members of one file side by side will build on;
here the two run one after the other, on the thread that reads.

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

The members of a bgzipped source. `many.vcf.gz` is 21904 bytes and its
members start at the bytes 0, 310, 12336 and 21876, so a cut at 315 or at
325 falls inside the header of its second member and one at 21903 a byte
before the end of the file: each gives the variants of the members that
were whole and then the error of the mark that is missing, and a cut
inside a header gives no variant of the member it cuts. A cut inside its
first member, at 100 bytes, is that error too and not one of the header of
the VCF. These are read and
are not errors: a member whose extra field holds another subfield before
its `BC`, in the first member and in a later one; a member of 65536 bytes
of text, which is the most one holds;
and a member with no text in the middle of a file, which bgzip can write
and which is the end of a file only when nothing follows it.

These are errors that name the member, each with a test that asserts the
words of that error, so that a check taken out of the reader is seen:
`many.vcf.gz` with its bytes 320 and 321 changed from `06 00` to `44 54`,
which is the file of the review, read at `VcfReader::new` or at
`next_block`; the same file with the `BC` of its first member moved behind
another subfield, which a reader that looked for the `BC` at the bytes 12
and 13 would read as a plain gzip and give no variant and no error; and a
file written in the test whose second member has, each in a case of its
own, a size two bytes too small, a size two bytes too large, a size that
leaves no room for its data, a length of its text that is not the one of
its text, a CRC32 that is not the one of its text, bytes where a gzip
member has those of gzip, a method that is not deflate, flags that are not
the one flag of an extra field, an extra field with no `BC` in it, a
subfield that ends after the extra field does before its `BC` and another
after it, an extra field that ends in the middle of a subfield, three bytes
of junk after a deflate stream that ends where it should, a deflate stream
cut by three bytes whose CRC32 and length are those of the text it gives,
a length of its text above the 65536 bytes a member holds, and data that
gives more text than that. And a file with 11 bytes after the member that
marks its end, which is corrupted and not cut short.

That an error never passes silently is tested on `cases.vcf.gz`, 399
bytes: every byte of it in turn, set to each of the 255 other values,
101745 files, each read whole with `only_passed` false. Each one gives
either an error or exactly the four variants of "What it gives", with
their genotypes; none gives other variants, or fewer, with no error. The
time stamp in the header of a member is among the bytes that are changed,
and a file that differs from `cases.vcf.gz` in it alone is read: what the
test refuses is a file that is read as a whole one and is not. It changes
one byte and never the place or the number of the members, which is the
corruption no reader of a BGZF file can see; and it says nothing about
which check of a member refuses which file, so every check has a test of
its own that asserts the words it gives.

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
ploidy is out of range, or the size of the blocks that the caller asked
for is one that `docs/specs/block.md` refuses: 0 variants, or a block of
more genotypes than a `usize` holds.

A `num_vars_per_block` of `None`, the size that popnei chooses, is not
checked there but when the first block is built, and the error is the same
one, with the words of a size that the caller did not write: it says that
the size popnei chose for these individuals and this ploidy does not fit
and that a `num_vars_per_block` that does is the way out, where the error
of a size that was asked for says to ask for fewer variants in a block. So a caller that opens a file to read its individuals, which is what
`open_vcf` and `openVcf` do, never fails for a size that nobody asked for:
a header of 170000 individuals read with the ploidy 255 gives a default
block of 100 variants whose genotypes are more than the 4295 million that a
`usize` holds in wasm, and such a file is opened, its individuals read, and
its blocks then asked for in a size that fits. It was decided here, when a
review of the binding crates found that file refused at `openVcf` although
blocks of ten of its variants are read.

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

The cases this module adds to the error of the crate, with the exception
each one is in Python. The owner gave the convention on 21 September 2026:
a `ValueError` is a wrong input of a function, a `RuntimeError` a defect
of popnei, and an `OSError` a file that cannot be read, that was cut short
or that is corrupted.

Five are a `ValueError`, since a file whose content is not what a VCF
holds is a wrong input like a wrong argument: the source is not a
VCF, with what was found; a wrong header, with what is wrong; a ploidy
out of range, which is the one thing `new` refuses that is not in the
source, with the ploidy that was asked for; a wrong data line, with the
number of the line, the column or the individual, and what is wrong; and a
genotype of another ploidy, with the line, the individual, the ploidy of
the genotype and the one expected.

Four are an `OSError`: a file that could not be opened, with
its path and the `std::io::Error` as the source of the error, so that a
binding can put the path where the language of the binding keeps it,
`OSError.filename` in Python; an error of the input, which wraps
`std::io::Error`; a bgzipped source with no mark of its end; and a
bgzipped source that is corrupted, with the member, counted from 1, the
byte of the compressed source where that member starts, and what is wrong
with it.

One is a `RuntimeError`: a parse of a batch that did not come back, with
the number of the last line that was read, which says that popnei has a
defect and not that the file or the call was wrong.

In Python every error of a file names the file. The core does not have the
path, since a reader is built over bytes and `from_path` carries it in one
case alone, so it is the binding crate that puts it there, and where it
puts it follows the exception. The message of a `ValueError` starts with
the path. An `OSError` carries it in `filename`, which is where a Python
user of any library looks for it and which Python prints after the message
of the exception, so putting it in the message too would say it twice.

## Speed

The number to reach is that of the spike, the trial parser in Rust of
section 3 of `docs/rust_core.md`, which parses a chunk of lines with rayon
and writes each row straight into the array of the chunk, as this reader
now does. The reader as built, which filled one `Variant` at a time, did
not reach it, and the reader as it is now does. Measured on the owner's
Apple M5 Pro, 18 cores, release, the file in the page cache, `GTS` asked
for and the default options, the median of 5 runs of
`crates/popnei/benches/read_vcf.rs`, on a VCF of 100000 variants and 1000
individuals, 403 MB plain and 38 MB bgzipped, which
`crates/popnei/benches/make_big_vcf.py` makes with `simulate_genotypes` and
`write_vcf` of pyNei's `test/gwas_reference/make_reference.py`, a seed of
42, 3 in 100 genotypes missing and `.` in every FILTER, so the default
gives every variant. The reader as built and the spike were timed on 20
September 2026, the reader as it is now on 21 September 2026, three sets
of runs of each and more:

| | the reader as built | the spike, 20 September | the target | the reader now | met |
|---|---|---|---|---|---|
| plain, 1 thread | 1.24 s | 0.54 s | 0.594 s | 0.563 s | yes |
| plain, 18 threads | 0.160 s | 0.098 s | 0.108 s | 0.093 s | yes |
| bgzipped, 1 thread | 1.58 s | 0.84 s | 0.924 s | 0.890 s | yes |
| bgzipped, 18 threads | 0.50 s | 0.40 s | 0.44 s | 0.394 s | yes |

The spike was timed again on 21 September 2026, on the same files and the
same machine, and was faster on three of the four: 0.523 s, 0.103 s,
0.806 s and 0.392 s. Against those the target of the bgzipped read on one
thread is 0.887 s, which the reader misses by 0.003 s; the fourteen sets
of runs of that read spread from 0.859 to 0.940 s, so the measurement does
not tell the two apart. The other three are met against the spike of
either day, and on 18 threads the reader is faster than the spike.

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

Where the time goes, from a sampling profile of the run on one thread
taken with `/usr/bin/sample` over 30 s on 21 September 2026, of 23038
samples of the thread that reads: 92 in 100 of the self time in the
columns of the individuals; 6 in the read of the lines, which is the
search for the end of a line in the buffer of the file, 3, the read from
the file system, 2, and the copy of the line out of that buffer, 1; 1 in
the nine first columns, which are parsed as text and not as bytes; and
under 1 in 100 in the pass that finds the FILTER and in the genotypes of a
new block set to missing. Nothing is appended serially after a batch here,
since the benchmark asks for the genotypes alone.

The read of the lines is serial, and on 18 threads it is what bounds the
reader: it is 4308 of the 9363 samples of the thread that reads, and the
18 workers of rayon are idle, waiting, in 61 in 100 of their samples. The
read ahead thread of section 3 of `docs/architecture.md` is what removes
that floor, and the targets above are met without it.
`docs/reports/vcf-to-blocks.md` has the profile of the reader as built,
which spent 95 in 100 of the one thread in the parse and split the columns
of the individuals as text.

The two bounds of a batch and the buffer the reader opens a path with were
measured on 21 September 2026, on the file of the table, on 18 threads,
three sets of runs of each, interleaved. Batches of 256, 1024, 2048 and
4096 lines: 0.139, 0.105, 0.098 and 0.098 s, and bgzipped 0.476, 0.425 and
0.411 s for 256, 1024 and 4096. The bound in bytes at 8 MiB and at 16 MiB
with 4096 lines: 0.098 s and 0.094 plain, 0.411 s and 0.392 bgzipped. The
buffer of the file at 8 KiB, 64 KiB, 256 KiB and 1 MiB, with 4096 lines
and 16 MiB: 0.105, 0.095, 0.093 and 0.093 s. On one thread none of the
three changes the read. The constants are 4096 lines, 16 MiB and 256 KiB,
each with its measurement in its doc comment in
`crates/popnei/src/io/vcf.rs`.

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
  `docs/specs/js_sources.md`, which says what `openVcf` takes in the wasm
  package and what the source tells the page while it reads.
