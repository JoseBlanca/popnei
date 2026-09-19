# popnei: objectives

popnei is a population genetics library written in Rust and used from
Python, natively and in the browser under pyodide. It is the successor of
pyNei, a working library in pure Python with numpy, pandas and pyarrow that
reached the speed floor of interpreted array code, and it is named, as
pyNei was, after Masatoshi Nei.

## What it is for

The calculations a population geneticist runs over variants: reading a VCF
and keeping the variants in a fast file of its own, filtering by missing
data, allele frequency, heterozygosity and linkage disequilibrium, per
variant and per sample statistics, allele frequencies and diversity per
population, distances between samples and between populations, principal
components and coordinates, linkage disequilibrium, kinship, and
association of the variants with a trait, continuous or binomial, with and
without a kinship. Datasets of tens of samples to about ten thousand and
of tens of thousands of variants to a million, on a personal computer or
in a browser tab.

## The goals, in order

1. **Right.** Every calculation is verified against a reference outside
   the project, plink2, GMMAT, rrBLUP, R, and against pyNei, which is
   verified against them already. The reference tools are run once by a
   script kept in the repository, their outputs are stored beside it, and
   the numbers of a few cases are written into the tests as literals, so
   that a test says what it expects and never needs the tool. pyNei is a
   development dependency, a path dependency on its sibling checkout, and
   the tests run both libraries on the same inputs where they overlap.

2. **Usable from Python, the way pyNei is.** The public API mirrors
   pyNei's where it holds up: functions over a `Variants`, `pops` as a dict
   of name to samples, frozen result dataclasses with pandas frames and
   series, sample names as tuples. It is allowed to diverge where the new
   design asks for it, and when it does the divergence is written down.
   The Python layer is thin: the API, the result objects, the tests.

3. **In the browser.** The core builds as a wasm wheel for pyodide, single
   threaded, with no C dependency, and every calculation works there with
   the memory a tab has. The wasm wheel is a release artifact tied to the
   pyodide version and its emscripten, rebuilt for every pyodide release.

4. **Fast where it matters.** A VCF parsed at the speed of compiled tools,
   per variant work in fused passes over the genotypes with rayon across
   records, and the linear algebra on the system BLAS natively and on
   faer in wasm. Fast is measured, on stated datasets, against pyNei and
   against plink2 and GMMAT, before and after every change that claims it.

5. **Streaming.** The variants go through the library as a stream of
   records grouped in blocks of a few thousand, so that a dataset never has
   to fit in memory, the per record work needs no arrays, and the matrix
   bound work gets the blocks it needs.

## Non goals

- Not a general purpose genomics toolkit. No alignment, no variant
  calling, no annotation.
- Not a command line program in the first version. The core crate is kept
  free of Python so that one can be added, but the interface is Python.
- Not a rewrite of pyNei's internals. pyNei's algorithms and results are
  the specification, its code is documentation, and its chunk design
  informs the stream of blocks, but the code starts empty.

## The design decisions

They are in `docs/rust_core.md`, with the measurements that led to each of
them: a core crate in pure Rust and a binding crate with pyo3 in one cargo
workspace, built by maturin into one wheel with the Python package; the
parser first; a stream of blocks; one small linear algebra module with
BLAS and LAPACK natively and faer in wasm; rayon for the records and BLAS
for the products, never nested; 2 bit packed genotypes as an option to
measure; the pyodide wheel pinned to the pyodide version.

## How the work is done

- Measure before deciding, and say what the dataset and the machine were.
- Tests must fail before the change they test, and when a change cannot
  have such a test, say so.
- Findings go on the GitHub issues, so that the reasoning survives.
- Commit messages: a lower case subject, a body with the why and the
  numbers.
