# popnei: objectives

popnei is a population genetics library written in Rust. It is used from
Python, natively and in the browser under pyodide, the Python that runs
in a tab, and from TypeScript in web applications, where the Rust is compiled to WebAssembly with no Python
in between. It is the successor of
pyNei, a working library in pure Python with numpy, pandas and pyarrow that
reached the speed floor of interpreted array code, and it is named, as
pyNei was, after Masatoshi Nei.

## What it is for

The calculations a population geneticist runs over variants: reading a VCF
and keeping the variants in a fast file of its own, filtering by missing
data, allele frequency, heterozygosity and linkage disequilibrium, per
variant and per individual statistics, allele frequencies and diversity per
population, distances between individuals and between populations, principal
components and coordinates, linkage disequilibrium, kinship, and
association of the variants with a trait, continuous or binomial, with and
without a kinship. Datasets of tens of individuals to about ten thousand and
of tens of thousands of variants to a million, on a personal computer or
in a browser tab.

## The goals, in order

1. **Right.** Every calculation is verified against a reference outside
   the project, plink2, GMMAT, rrBLUP, R, and against pyNei, which is
   verified against them already. The reference tools are run once by a
   script kept in the repository, their outputs are stored beside it, and
   the numbers of a few cases are written into the tests as literals, so
   that a test says what it expects and never needs the tool. pyNei is a
   development dependency, taken from its git repository at one commit, so
   that every machine tests against the same pyNei, and the tests run both
   libraries on the same inputs where they overlap.
   The TypeScript API is tested against the same literals, under node,
   the JavaScript runtime outside the browser.

2. **Usable from Python, the way pyNei is.** The public API mirrors
   pyNei's where it holds up: functions over a `Variants`, `pops` as a dict
   of name to individuals, frozen result dataclasses with pandas frames and
   series, the names of the individuals as tuples. It is allowed to diverge where the new
   design asks for it, and when it does the divergence is written down.
   The Python layer is thin: the API, the result objects, the tests.

3. **In the browser, for two users.** The first is a web application
   written in TypeScript. For it the core is compiled directly to
   WebAssembly and published as an npm package, the wasm package: the
   compiled core, the JavaScript that loads it and the TypeScript
   declarations of its functions. The tab loads no Python, so the
   application does not wait for pyodide, numpy and pandas to download
   and start, and its build does not follow the pyodide releases. The
   TypeScript API carries the names of the Python one, its functions,
   arguments and result fields, and gives the numbers as typed arrays
   where Python gives pandas frames. The second user writes Python in a
   notebook that runs in the tab, and gets the Python package under
   pyodide, with the core built as a wasm wheel. That wheel is a release
   artifact tied to the pyodide version and its emscripten, rebuilt for
   every pyodide release. Both builds are single threaded in the first
   version and have no C dependency, and every calculation works in both
   with the memory a tab has. Whether the wasm package gets threads later
   is an open question of `docs/rust_core.md`.

   The browsers popnei runs in are those that have the vector
   instructions of WebAssembly, the ones that work on sixteen bytes at a
   time, which its calculations use: Chrome and Edge from 91, of May
   2021, Firefox from 89, of June 2021, and Safari from 16.4, of March
   2023, which on an iPhone or an iPad means iOS 16.4, every browser
   there being WebKit whatever its name; outside the browser, node from
   16.4, of June 2021. The owner set that floor on 22 September 2026,
   when the performance review of the Kosman distances asked for those
   instructions, and it is the first minimum popnei writes down. The
   option not taken was to ship the wasm package twice, with and without
   them, and pick at load: that keeps every browser and doubles the
   bytes of the package.

4. **Fast where it matters.** A VCF parsed at the speed of compiled tools,
   per variant work in fused passes over the genotypes with rayon across
   records, and the linear algebra on the system BLAS natively and on
   faer in wasm. Fast is measured, on stated datasets, against pyNei and
   against plink2 and GMMAT, before and after every change that claims it.

5. **Streaming.** The variants go through the library as a stream of
   blocks of a few thousand, so that a dataset never has to fit in memory,
   the per variant work has the rows to spread over the threads, and the
   matrix bound work gets the blocks it needs.

## Non goals

- Not a general purpose genomics toolkit. No alignment, no variant
  calling, no annotation.
- Not a command line program in the first version. The core crate is kept
  free of Python and of JavaScript so that one can be added, but the
  interfaces are Python and TypeScript.
- Not a rewrite of pyNei's internals. pyNei's algorithms and results are
  the specification, its code is documentation, and its chunk design
  informs the stream of blocks, but the code starts empty.

## The design decisions

They are in `docs/rust_core.md`, with the measurements that led to each of
them: a core crate in pure Rust and a binding crate with pyo3 in one cargo
workspace, built by maturin into one wheel with the Python package; a
second binding crate, with wasm-bindgen, built into the wasm package for
TypeScript; the parser first; a stream of blocks; one small linear algebra module with
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
