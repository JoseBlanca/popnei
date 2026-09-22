# The names popnei uses

September 2026. popnei is a population genetics library in Rust, used from
Python and from TypeScript, the successor of the Python library pyNei. Its specs, its plans and
its code are written by different people and sessions, and this document
fixes which word each of them uses for the things of the domain, so that a
thing has one name in the docs, in Rust, in the binding and in Python. It
decides names and nothing else. How the variants flow is in
`architecture.md`, and what a calculation does is in its spec under
`docs/specs/`.

A document still says what a term means where it first uses it, as the
`writing` skill asks, because its reader may have only that page.

Each entry has the name, what it is, the name pyNei gives it when that
differs or when it is an identifier to look for, and the words that are not
used for it. The pyNei names were read in pyNei at commit ef0ca6e.

## The data

**individual.** One organism that was genotyped, known by its name. It is
the word of population genetics. VCF calls it a sample, and it is a
column of a VCF, and pyNei calls it a sample too, `samples`,
`num_samples`, `filter_samples`. In popnei the word is individual in the
prose and in the identifiers, `individuals`, `num_individuals`, and
"sample" is written only for a name of pyNei and for a statistical
sample. Inside the core an
individual is an index into the individuals of the reader. Not used:
sample, accession.

**population.** A named set of individuals that a calculation treats as a
group. "pop" in identifiers, `pop`, `pops`, and "population" in prose.
`pops` is a dict of population name to the names of its individuals in
Python, as in pyNei, and the indices of the individuals of each
population inside the core. Not
used: group, deme.

**variant.** One site of the genome with its alleles and the genotype of
every individual at it, a data line of a VCF. It can have more than two
alleles. Not used: SNP, which is one kind of variant, marker, locus.

**allele.** One of the forms a variant has, held as an `i8`: 0 is the
reference allele and 1 and above are the alternative ones, in the order of
the VCF.

**missing allele.** An allele that was not called, `.` in a VCF. It is the
constant `MISSING_ALLELE`, -1, in pyNei and in popnei, and never a literal
-1.

**ploidy.** How many alleles the genotype of one individual holds, the
same for every individual and every variant of a dataset.

**genotype.** The `ploidy` alleles of one individual at one variant, each of
which can be missing on its own. `gts` in identifiers.

**called genotype**, **missing genotype** and **half called genotype.** A
genotype is called when none of its alleles is missing, and missing when
at least one is, which is what `_calc_gt_is_missing` of pyNei's
`gt_counts.py` computes. A half called genotype, `0/.` in a VCF, has some
alleles called and some missing, so it is a missing genotype: it adds to
the missing rate of its variant, the missing genotypes divided by the
individuals, and it is not among the genotypes of the observed
heterozygosity. The calculations that count alleles, the
allele frequencies and the expected heterozygosity, still count the
alleles of it that were called, and their spec says so. Not used: partial
genotype, no call.

**observed heterozygosity.** Of a variant over some individuals, the
heterozygous genotypes divided by the called genotypes, where a genotype
is heterozygous when it is called and its alleles are not all the same.
`obs_het` in identifiers.

**called alleles.** How many alleles of a population at a variant are not
missing, the denominator of its allele frequencies.

**expected heterozygosity.** Of a variant in a population, one minus the
sum over the alleles of the frequency of each to the power of the ploidy:
the chance that as many gene copies as a genotype holds, taken at random
from the population, are not all alike. The unbiased one multiplies it by
2n/(2n - 1), with n the called genotypes, and is a statistic of its own.
`exp_het` and `unbiased_exp_het` in identifiers. `docs/specs/stats.md`.

**polymorphic variant.** In a population, a variant whose major allele
frequency is below the polymorphism threshold, 0.95 by default, strictly;
a variable one has it below 1. `poly` in identifiers, as in pyNei's
`poly_threshold` and `num_poly`.

**major allele.** The allele of a variant with the highest frequency among
the called alleles of the individuals considered. How a tie is broken is for
the spec of the calculation to say.

**maf.** In pyNei and in popnei, the frequency of the major allele, as in
`filter_by_maf` and its `max_allowed_maf`. In most of the literature and in
plink2 the same three letters are the frequency of the minor allele, so a
text writes "the major allele frequency" in full where it first uses it.

**dosage.** For one genotype, how many of its alleles are not the major
allele of the variant, from 0 to the ploidy, every allele other than the
major one counting the same. The dosage matrix is the variants x
individuals array of them. pyNei: `to_012` and "the 012 matrix".

## How the data moves

**block.** Consecutive variants held as contiguous arrays, the `Block`
struct, which is how the variants flow from a source to a calculation: a
reader gives blocks, a filter compacts them and a calculation consumes
them, as rows or as a matrix. One variant of a block is a view into it,
`VariantRef`, and "row" is that variant as a line of the arrays. A block holds
about 5 million genotypes, the size pyNei gives its chunks, which is a few
thousand variants.
pyNei: chunk, `VariantsChunk`. "Chunk" is used only for pyNei's own. What
the Python `Variants` of popnei gives from `iter_blocks` is a block. Not
used: batch, which is written for two other things, arrow's unit of the
vars file, and the lines that the VCF reader reads and parses together,
several to a block; and window. A vars file is written
with one batch for each block, and read back in blocks of any size.

**member.** One gzip stream of a gzipped file. A file that bgzip wrote is
many of them one after another, each with 64 KiB of text at most and each
stating its own size, and the reader of a VCF cuts them by those sizes and
decompresses one at a time, as `docs/specs/io_vcf.md` describes.
BGZF: block. "Block" in popnei is a block of variants, so the word for a
gzip stream of such a file is member, which is what the gzip format calls
it, and the empty member at the end of a bgzipped file is the mark of its
end and not the empty block.

**region.** A stretch of one chromosome, from a smallest to a largest
position, both included. The vars file keeps, for each of its batches, the
region of every chromosome that has variants in it.

**reader.** Anything that gives blocks through the `BlockReader` trait: the
VCF reader, the vars file reader, and a filter or `reblock`, which are
readers over another reader. Not used: record, record level and block
level, the words of the first version of `architecture.md`, which had a
single variant that a reader filled.

**vars file.** The arrow file in which popnei keeps variants, feather v2,
which any program with an arrow library opens as a table. Its format is
popnei's own, in `docs/specs/io_vars.md`. pyNei has a file of the same
name and another format, and neither library reads the other's.

## The layers

popnei is used from two languages, Python and TypeScript. Both reach the
same core crate, each through a binding crate with a package on top, as
section 8 of `architecture.md` lays out.

**core crate.** `crates/popnei`, pure Rust with no Python and no
JavaScript in it, where every calculation is.

**binding crate.** A crate that translates between another language and
the core crate and holds no calculation. There are two, and a text that
**step.** One entry of the list that a `Variants` holds besides its
source, a filter with its threshold or the filter of individuals with
its names. A step is added with a method of the
`Variants`, which returns nothing, it is run in every pass that starts
after it was added, and `variants.steps` lists them.

**consumer.** What runs a `Variants`: the function of a calculation,
`write_vars`, the method `iter_blocks`. It does not change the `Variants`,
and what it returns has the pass stats.

**pass stats** and **filtering stats.** What a pass counted, `PassStats`,
which every result of a consumer has as `pass_stats`: how many variants
the consumer took, and the filtering stats, `FilteringStats`, how many
variants each filter was given and how many it kept. pyNei keeps the
filtering stats in its `Variants` and gives them with
`gather_filtering_stats`.

**pass.** One reading of a source of variants from its start to its end,
through the steps that its `Variants` had when it started. A consumer
makes as many as its algorithm needs, each with readers and filters of its
own.

could mean either says which. "The binding crate" alone is used for what
holds for both.

**Python binding crate.** `crates/popnei-python`, written with pyo3. Its
Python module is `popnei._core`.

**JavaScript binding crate.** `crates/popnei-js`, written with
wasm-bindgen, the Rust tool that generates the JavaScript that calls the
exported Rust functions and their TypeScript declarations.

**Python package.** `python/popnei`, what a Python user imports: the
functions with pyNei's signatures, the result objects and the Python
`Variants`.

**TypeScript package.** `js/popnei`, what a web application imports: the
functions of the Python API with their names in camelCase, and the result
objects with typed arrays.

**pyodide wheel.** The Python binding crate and the Python package built
for pyodide, the Python that runs in a browser tab, for the person who
writes Python in a notebook there. Compiled for
`wasm32-unknown-emscripten`.

**wasm package.** The npm package that a web application installs: the
TypeScript package with the JavaScript binding crate compiled for
`wasm32-unknown-unknown` inside. No Python runs in the tab. Not used:
"the wasm build" when the text could mean either this or the pyodide
wheel.

## How results are verified

**reference program.** A program outside the project whose output a
calculation of popnei is checked against: plink2, GMMAT, rrBLUP, R. Not
used: tool, gold standard.

**literal.** A number from a reference program or from pyNei that is
written into a test as it was printed, with the spec saying how it was
got. A test compares against literals and does not compute what it
expects.

## When a name is missing or wrong

Where pyNei has a name for a thing that a user sees, a function, an
argument, a field of a result, that is the name. A thing of the domain
that several documents will speak of and that has no entry here gets one
in the commit that first needs it. An entry that turns out wrong is
changed here first, and then wherever a search for the old word finds it.
