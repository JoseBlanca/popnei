# The io::vcf module: the VCF reader

September 2026. The VCF reader is how variants get into popnei: it reads a
VCF, plain or gzipped, and gives its variants one at a time, with the
genotypes as small integers. There is no code. This spec develops the
reader of the row `io::vcf` of the table in section 9 of
`docs/architecture.md`. It depends on `docs/specs/variant.md`, which has
the `Variant` that the reader fills, the `Needs` that say which fields a
consumer wants, and the `VariantReader` trait that the reader implements.
The VCF writer of the same row is an item that is not written yet.

An individual is what `docs/glossary.md` calls one organism that was
genotyped. VCF calls it a sample, and the columns of a VCF after FORMAT
are here the columns of the individuals.

## The VCF reader

### What it gives

From each data line of a VCF, one variant:

| column of the VCF | field of the `Variant` |
|---|---|
| CHROM | `chrom`, the number of the name in the reader's table of chromosomes |
| POS | `pos` |
| ID | `id`, empty when the column is `.` |
| REF and ALT | `alleles`, the reference first; only the reference when ALT is `.` |
| QUAL | `qual`, `None` when the column is `.` |
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
and nothing in the `Variant` says which ones had failed. The owner
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
keeps no phase, because a `Variant` has nowhere to hold it. VCF 4.4 lets a
genotype start with a separator, `/0/1` or `|1|1`, and the reader reads
those as `0/1` and `1|1`. An allele written `.` is `MISSING_ALLELE`, -1.

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
the bytes of the file or a `File` that the user picked in the page, which
can be read only inside a web worker, as section 11 of the architecture
says.

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
ten columns are errors. A line whose bytes are not text, not valid UTF-8,
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

### How it runs

At the record level, as a `VariantReader`. Natively it reads a batch of
lines, parses them with rayon into variants of its own, and hands them
out in the order of the file, one for each `read_variant`, the method of
the trait that fills the next variant. What crosses is the lent
`Variant`, and once the buffers of the reader and of the lent `Variant`
have grown to the size of a line nothing is allocated from one variant to
the next; a swap of the buffers of the two variants does that without a
copy. In wasm the same code parses one line after another. How many lines
a batch holds is a number to measure.

The result does not depend on the number of threads. The variants come in
the order of the file, the chromosome numbers are given when a variant is
handed out and not when it is parsed, and when a line is wrong the
variants before it are given first and the error comes at the
`read_variant` that would have given that line.

Only what `Needs` asks for is parsed, except the chromosome and the
position, which are always filled. What is checked whatever is asked for
is the shape of the line: the nine first columns have to be there, the
FORMAT has to have a `GT` key, and there has to be one column after the
FORMAT at least. With the genotypes not asked for, the reader does not
look at the columns of the individuals: how many of them there are and
what is in them is not read, so a line with two columns of individuals
under a header with three, or one with a genotype of another ploidy, is
given.

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
115 KB, and `many.vcf.gz`, which bgzip wrote as 4 gzip members with 617,
65250, 49014 and 0 bytes of text, so a decoder that stops at the first
one fails these tests. 450 of its variants have `PASS`, 25 have `.` and
25 have `q10`; the first with `q10` is chr1 1259. Every genotype,
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
variant, which gives false at the first `read_variant`; a variant whose
ALT declares two alleles and whose genotypes carry only the first; a
line with `q10` and a tetraploid genotype, which the default skips and
`only_passed` false refuses; a line with `q10` and the position `x`,
which the default skips; a last line with no end of line; and a source
that gives one byte at a time, gzipped and plain, which a reader that
looked for the two bytes of gzip in one look at the buffer would refuse.

What `Needs` does: with `GTS` alone, `filled` has the genotypes, the
chromosome and the position and no more, and `alleles` is empty; with
`ID` and `ALLELES` and no `GTS`, `gts` is empty; and a reader asked for
everything and then for `GTS` alone leaves the alleles of the variant it
filled before empty. That a second pass over
`many.vcf` allocates nothing is checked once by hand with a counting
allocator when the reader is written, and it is not a test that stays.

The cargo tests are made at `VcfReader::new` for what is wrong in the
header, the source that is not a VCF, the FORMAT column or the
individuals that are not there, the repeated name, the name that is not
there and the ploidy out of range, at `from_path` for the file that is
not there, and at
`read_variant` for the rest, with the reader built over the bytes of the
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
}
impl Default for VcfOptions { /* 2, true */ }
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
returns, and it fails when the source is not a VCF with genotypes or the
ploidy is 0.

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

impl<R: BufRead + Send> VariantReader for VcfReader<R> { /* ... */ }
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
`OSError.filename` in Python; and an error of the input, which wraps
`std::io::Error`. In Python the first five are a `ValueError` and the
last two an `OSError`.

## Speed

The dataset and the numbers of section 3.1 of `docs/rust_core.md`, on the
owner's Apple M5 Pro: a VCF of 100000 variants x 1000 individuals,
400 MB, which pyNei parses in 13.5 s, plink2 in 0.27 s on one thread, and
the spike in 0.55 s on one thread and 0.11 s on 18 cores; gzipped, 53 MB,
the spike took 0.55 s on one thread too, and with threads the
decompression, which is serial, is what bounds it. The spike filled the
genotypes, the chromosome and the position, as a calculation that asks
for `GTS` does, and it did three things less than the reader: it checked
neither the ploidy of each genotype nor its allele numbers against ALT,
and it did not read the FILTER. `docs/rust_core.md` gives
`test/gwas_reference/make_reference.py` of pyNei as what simulated that
panel, and its `write_vcf` puts `.` in the FILTER of every variant, so
the default gives them all and the two do the same work on the
genotypes. The number to reach with `GTS` asked for and the default
options is the spike's, 0.55 s on one thread and 0.11 s on 18, or no more than a
tenth above them. Whether the checks cost more than that is what the
first measurement says, and it comes before any work on speed.

## Open points

None. The owner decided on 20 September 2026 the four that there were:
the ploidy as an argument and the refusal of mixed ploidies, the variants
that failed a filter left out by default, an allele number that is not
declared as an error, and no error for a VCF with no variants. Each is
written where it applies, with the option that was not taken.

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
