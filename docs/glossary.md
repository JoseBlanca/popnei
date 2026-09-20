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

**sample.** One individual that was genotyped, a column of a VCF, known by
its name. Inside the core a sample is an index into the samples of the
reader. Not used: individual, accession.

**population.** A named set of samples that a calculation treats as a
group. "pop" in identifiers, `pop`, `pops`, and "population" in prose.
`pops` is a dict of population name to sample names in Python, as in
pyNei, and the sample indices of each population inside the core. Not
used: group, deme.

**variant.** One site of the genome with its alleles and the genotype of
every sample at it, a data line of a VCF. It can have more than two
alleles. Not used: SNP, which is one kind of variant, marker, locus.

**allele.** One of the forms a variant has, held as an `i8`: 0 is the
reference allele and 1 and above are the alternative ones, in the order of
the VCF.

**missing allele.** An allele that was not called, `.` in a VCF. It is the
constant `MISSING_ALLELE`, -1, in pyNei and in popnei, and never a literal
-1.

**ploidy.** How many alleles the genotype of one sample holds, the same
for every sample and every variant of a dataset.

**genotype.** The `ploidy` alleles of one sample at one variant, each of
which can be missing on its own. `gts` in identifiers.

**called genotype**, **missing genotype** and **half called genotype.** A
genotype is called when none of its alleles is missing, and missing when
at least one is, which is what `_calc_gt_is_missing` of pyNei's
`gt_counts.py` computes. A half called genotype, `0/.` in a VCF, has some
alleles called and some missing, so it is a missing genotype: it adds to
the missing rate of its variant, the missing genotypes divided by the
samples, and it is not among the genotypes of the observed
heterozygosity. The calculations that count alleles, the
allele frequencies and the expected heterozygosity, still count the
alleles of it that were called, and their spec says so. Not used: partial
genotype, no call.

**called alleles.** How many alleles of a population at a variant are not
missing, the denominator of its allele frequencies.

**major allele.** The allele of a variant with the highest frequency among
the called alleles of the samples considered. How a tie is broken is for
the spec of the calculation to say.

**maf.** In pyNei and in popnei, the frequency of the major allele, as in
`filter_by_maf` and its `max_allowed_maf`. In most of the literature and in
plink2 the same three letters are the frequency of the minor allele, so a
text writes "the major allele frequency" in full where it first uses it.

**dosage.** For one genotype, how many of its alleles are not the major
allele of the variant, from 0 to the ploidy, every allele other than the
major one counting the same. The dosage matrix is the variants x samples
array of them. pyNei: `to_012` and "the 012 matrix".

## How the data moves

**record.** A variant as one item of the stream that a reader gives, the
`Variant` struct that the caller owns and the reader fills. The genetics is
written with "variant". "Record" is for the flow of the data.

**block.** Consecutive variants held as contiguous arrays, the `Block`
struct, which the calculations that want matrices consume. A block holds
about 5 million genotypes, the size pyNei gives its chunks, which is a few
thousand variants.
pyNei: chunk, `VariantsChunk`. "Chunk" is used only for pyNei's own, and
for the arrays that the Python `Variants` of popnei yields, which are
pyNei's chunks. Not used: batch, which is arrow's word for the unit of the
vars file, and window.

**record level** and **block level.** The two ways a calculation runs, of
sections 1 and 2 of `architecture.md`. At the record level it sees one
variant at a time and what it keeps from one variant to the next does not
grow with the number of variants. At the block level it gets blocks.

**reader.** Anything that gives variants one at a time through the
`VariantReader` trait: the VCF reader, the vars file reader, and a filter,
which is a reader over another reader.

**vars file.** The arrow file in which pyNei and popnei keep variants,
which each of them has to read as the other writes it. Its format is in
section 6 of `architecture.md`.

## The layers

popnei is used from two languages, Python and TypeScript. Both reach the
same core crate, each through a binding crate with a package on top, as
section 8 of `architecture.md` lays out.

**core crate.** `crates/popnei`, pure Rust with no Python and no
JavaScript in it, where every calculation is.

**binding crate.** A crate that translates between another language and
the core crate and holds no calculation. There are two, and a text that
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
