# A Rust core: the decision and the measurements behind it

September 2026. This is the founding document of popnei, the Rust based
successor of pyNei. It records where the pure Python library stood, what
was measured, what was concluded from it, and the decisions the new
project starts from. It was written in the pyNei repository, whose issues
#17 and #18, whose `spike/` directory and whose commits hold the details.
Every number here was measured on one machine, an Apple M5 Pro with 18
cores and 64 GB, with numpy 2.5 on Accelerate, Python 3.14.5, Rust 1.98
stable, and pyodide 314.0.7 under node 26.

## 1. Where pynei stood

pynei is a population genetics library in Python with numpy, pandas and
pyarrow as its only dependencies, deliberately, so that it runs in the
browser under pyodide. The variants are processed chunk by chunk, so a
dataset does not have to fit in memory, and the chunks are sized by their
genotypes, about 5 million per chunk. A genotype is one int8 per allele,
-1 for a missing one, with no mask beside the values. The variants live in
one arrow IPC file, zstd compressed, one record batch per chunk, and reading
a chunk from it takes 0.17 s per 100 million genotypes where parsing them
from a VCF takes 13.7 s. Every calculation is one pass over the chunks,
mapped and reduced, with a thread per chunk when asked, and the next chunk
is read in a thread of its own while the one in hand is worked on.

The last thing built in Python was the GWAS: a genomic relationship matrix,
linear and logistic regressions, and their mixed model versions with the
kinship as a random effect, verified against plink2, GMMAT, rrBLUP and R's
glm to the digits those programs print. That work is what produced the
measurements below, because it put pynei next to compiled tools on the same
data.

## 2. The numpy floor

### 2.1 What numpy does per genotype

A per variant calculation in numpy is a sequence of whole array operations,
and each one is a pass over the chunk. The first version of the dosage
matrix, `to_012`, reduced along the ploidy axis, `sum(gts != major,
axis=2)`, and numpy reduces badly along an axis of length two: 1.8 ns per
genotype, 18 ms per chunk of 5000 variants x 1000 individuals, for each of the
two reductions. The same thing as operations between the two planes of
alleles, `gts[:, :, 0]` and `gts[:, :, 1]`, took 2.5 ms. A fast path for the
common case, only alleles 0 and 1 and no missing genotype, which one `max`
and one `min` detect in 0.2 ms, removed the general allele counting as well.

Over 100000 variants x 1000 individuals, one thread:

| | first version | plane by plane | fast path | plink2 |
|---|---|---|---|---|
| to_012 | 0.90 s | 0.29 s | 0.10 s | |
| kinship | 1.68 s | 1.00 s | 0.81 s | 0.23 s |
| linear GWAS | 1.22 s | 0.56 s | 0.35 s | 0.10 s |
| 012 matrix for the PCA | 0.75 s | 0.42 s | 0.22 s | |

After that the linear GWAS is 3.5x plink2 and nothing single dominates any
more: per chunk, the dosages 5 ms, the float conversion 1 ms, the mean and
variance 4 ms, the regression 4 ms. That is the floor: several passes over
memory at a nanosecond or so per genotype each, plus 8 bytes per genotype
written for the float64 dosages that BLAS needs.

### 2.2 Packing the genotypes in 2 bits does not help in numpy

plink keeps a genotype in 2 bits and runs its calculations on the packed
bits with popcount tables, in one fused pass. That was measured in numpy on
one chunk of 5 million genotypes:

| | packed 2 bit | float64 |
|---|---|---|
| pack or unpack | 0.4 ms / 0.5 ms | |
| memory | 1.25 MB | 5 MB int8, 10 MB genotypes, 40 MB float64 |
| x'y, the regression dot product | 2.7 ms via bit masks | 0.2 ms |
| dot products between 64 and 5000 variants | 8.6 ms for one of three popcount terms | 1.2 ms for all of it |
| dosage sums per variant | 0.2 ms with `bitwise_count` | 0.3 ms |

Packing is cheap and makes a chunk eight times smaller, but every
calculation on the packed bits loses to BLAS on floats by 4x to 20x, because
numpy cannot fuse the unpack, the popcount and the accumulation into one
pass, and the file does not shrink, zstd already holds a genotype in 1.4
bits. In numpy the fast path is the float path. In compiled code the packed
path is the fast one, which is the first thing a Rust core gets to revisit.

### 2.3 The VCF parser

Parsing a VCF is the most expensive thing pynei does, and it is text
handling in Python, one field at a time:

| | pynei, Python | plink2 |
|---|---|---|
| parse 400 MB, 100000 variants x 1000 individuals | 13.7 s | 0.27 s |

A 50x gap, against the 3.5x left on the linear regression.

### 2.4 Where numpy is not the problem

Everything BLAS and LAPACK bound is already at the speed of the library
numpy calls: the kinship product `Z'Z`, the eigendecomposition of the
kinship, the null fits of the mixed models, the per variant mixed model
test, which is a block of variants against an individuals x individuals projection,
and the PCA. On the mixed models pynei is level with GMMAT, C++ under R, on
one thread, 1.5 s against 1.6 s and 2.2 s over 100000 variants x 1000
individuals, and 3x faster with six threads. Once the individuals number in the
thousands these products are what dominates, and no language changes them.

### 2.5 statsmodels and scipy

Measured and not taken: a statsmodels fit per variant is 50x to 160x slower
than the same algebra vectorized per chunk, its mixed model cannot take a
kinship, and scipy would only have provided two distributions, erfc and the
incomplete beta, which are a few lines each and were checked against scipy
to 1e-12.

## 3. The Rust spike

`spike/pynei_spike` in the pyNei repository: a VCF chunk parser, and `z'z` and
the symmetric eigendecomposition through faer, exposed with pyo3 0.29 and
built with maturin, natively and as a wasm wheel for pyodide. The genotypes
it parses are identical to pynei's on the reference VCFs, missing genotypes
and gzip included; `z'z` is exact to the bit against numpy and the
eigenvalues agree to 1e-14.

### 3.1 Native

Same 400 MB VCF, 100000 variants x 1000 individuals:

| | Python pynei | Rust, 1 thread | Rust, rayon 18 cores | plink2, 1 thread |
|---|---|---|---|---|
| parse VCF | 13.5 s | 0.55 s | 0.11 s | 0.27 s |
| parse VCF gzipped, 53 MB | | 0.55 s | the inflate is the ceiling | |

The parser is 25x faster than Python on one thread and twice as fast as
plink2 with the threads.

Linear algebra, faer against numpy on Accelerate:

| | numpy | faer 1 thread | faer 6 threads |
|---|---|---|---|
| z'z, 5000 x 1000 | 0.01 s | 0.17 s | 0.03 s |
| eigh n = 1000 | 0.05 s | 0.13 s | 0.11 s |
| eigh n = 2000 | 0.33 s | 0.90 s | 0.44 s |
| eigh n = 5000 | 6.35 s | 14.3 s | 4.57 s |

faer does not beat Accelerate, which runs float64 products on the AMX
units of the M series: 17x slower on the product with one thread, 2x to 3x
on the eigendecomposition. That gap is Apple specific. faer's own
benchmarks put it level with OpenBLAS and MKL on x86 with AVX2 or AVX512;
that has not been measured here, and it is the first open question below.

### 3.2 Under pyodide

The crate builds as a wasm wheel with stable Rust 1.98, no nightly,
pyodide-build 0.39.0, the cross build env of pyodide 314.0.7 and emscripten
5.0.3, in 20 s once the toolchain is in place. The wheel is 193 KB and
installs with micropip in 0.13 s. The commands are in `spike/README.md`.
Same data on both sides, one thread everywhere:

| | native Rust | pyodide Rust | pyodide numpy |
|---|---|---|---|
| parse 10000 x 1000, 40 MB | 0.050 s | 0.059 s | |
| parse 2000 x 1000 gzipped | 0.015 s | 0.022 s | |
| z'z, 5000 x 1000 | 0.169 s | 0.58 s | 5.09 s |
| eigh n = 1000 | 0.124 s | 0.31 s | 0.69 s |
| eigh n = 2000 | 0.90 s | 2.39 s | 5.57 s |
| eigh n = 3000 | 2.97 s | 13.6 s | 19.2 s |

The parser loses 18% in wasm. faer, which loses natively, wins in the
browser, where numpy's pyodide build has no BLAS and runs a reference
LAPACK: 8.8x on the product, 2.2x on the eigendecomposition up to 2000
individuals, 1.4x at 3000, where the 32 bit heap starts to weigh. So a Rust
core keeps the browser and gains there.

Two traps found on the way. The host Python of pyodide-build must not be
the free threaded one, which `uv venv --python 3.14` picks on this
machine; with it pyodide-build puts the emscripten sysconfigdata under
`python3.14t` where nothing looks for it. And faer's rayon feature does not
build for emscripten, its `spindle` dependency needs `atomic-wait`, which
has no wasm platform, so rayon and that feature are pulled in only off the
wasm targets with `cfg(not(target_family = "wasm"))`.

## 4. The decisions

1. **A Rust core, with Python as the API, the results layer and the
   tests.** The public functions, the frozen result dataclasses and the
   pandas frames stay in Python, thin, over the core. The 181 Python tests,
   with the plink2, GMMAT, rrBLUP and R numbers written into them as
   literals, are the specification the core has to satisfy, and they do
   not change. The reference numbers were made once by a script kept beside
   them, `test/gwas_reference/make_reference.py`, and that is the pattern
   for every new verified calculation: run the reference tool once, keep
   its outputs and the script, write the numbers of a few cases into the
   test, never need the tool at test time.

2. **A parallel project, popnei, with pyNei kept as it is.** pyNei is
   working, verified software and a rewrite inside it would put every
   commit of the new thing next to the old, so popnei is its own
   repository. Inside it, one cargo workspace with two crates, a core
   crate in pure Rust with no pyo3, tested with `cargo test`, and a thin
   binding crate that maturin builds into one wheel together with the
   Python package, the mixed layout of pydantic-core and polars. pyNei is
   a development dependency of the Python side, taken from its git
   repository at the commit that `pyproject.toml` names, and the tests run
   both on the same inputs: pyNei is the oracle. The reference datasets and the script that made them are
   copied over, so the tests stay self contained. The core reads and
   writes the vars file with arrow-rs, in a format of its own since 20
   September 2026, when the owner decided that it owes pyNei's nothing;
   it is in `docs/specs/io_vars.md`.

3. **The parser goes first.** It is the largest gain by far, 25x on one
   thread and 123x with threads, it survives wasm almost whole, and it can
   sit behind `vars_from_vcf` without anything else changing. It handles
   the GT field, any ploidy, missing alleles, phased and unphased
   separators, and gzip, and it is verified against the Python parser on
   the reference VCFs before that parser goes.

4. **A stream of blocks.** Per variant work, the counts, the masks, the
   filters, the dosages, runs with rayon across the rows of a block. The
   BLAS bound work, the kinship, the mixed model test, the PCA, consumes
   the same blocks of a few thousand variants as matrices, because a
   variant at a time turns a matrix product into matrix vector products
   that run 5x to 10x slower. The block is the current chunk with its
   arrays made internal. Until 20 September 2026 this decision had the per
   variant work run record by record, on a single variant that each reader
   filled; section 1 of `docs/architecture.md` says why that was dropped.

5. **The linear algebra has two backends behind one small module.** The
   operations pynei needs are few: matrix product, symmetric
   eigendecomposition, Cholesky and solve, inverse, least squares.
   Natively they call BLAS and LAPACK through `blas-src` and `lapack-src`:
   Accelerate on macOS, a system framework with nothing to ship, and
   OpenBLAS bundled into the Linux and Windows wheels through
   `openblas-src` built statically, which needs a Fortran compiler in CI.
   That is the same library numpy calls, reached without Python in the
   loop. In wasm, where there is no BLAS, faer. Whether faer is also good
   enough natively on x86 is open, see below.

6. **Who owns the threads.** rayon owns the per variant work. BLAS has its
   own pool for the big products. The two must not nest: a rayon worker
   that calls BLAS pins it to one thread, and a big product that BLAS
   parallelizes is called from outside rayon. In wasm there are no threads
   at all, and the core must build and run single threaded there.

7. **2 bit packed genotypes become an option to measure**, not a decision.
   In compiled code the fused unpack, popcount and accumulate is the fast
   path, and it makes a chunk eight times smaller, which matters in the
   32 bit heap of the browser. The int8 per allele representation stays as
   the interchange format with Python and in the vars file until a
   measurement says otherwise.

8. **The pyodide wheel is a release artifact tied to the pyodide
   version.** Its ABI tag and its emscripten follow the runtime, so every
   pyodide release is a rebuild and a re-tag, and the wasm wheel is shipped
   outside PyPI, where micropip can fetch it. CI has to build the native
   wheels for three platforms and the wasm one, with the pinned emscripten
   and a non free threaded host Python.

9. **A second binding, for TypeScript, compiled directly to wasm.**
   Decided on 20 September 2026, after the rest of this document, and not
   measured. popnei will also be used from web applications written in
   TypeScript, and when the library is Rust a web application does not
   need Python to call it. A second binding crate, written with
   wasm-bindgen, which generates the JavaScript and the TypeScript
   declarations for the Rust functions it exports, is compiled for
   `wasm32-unknown-unknown` and published as an npm package. What the
   application is spared is pyodide itself, with numpy and pandas: their
   download and start before the first call, and a build tied to each
   pyodide release as in decision 8. What it adds is a second public API
   to keep in step with the Python one, its result objects in TypeScript,
   and a second wasm target that every dependency of the core has to
   build for. The wheel for pyodide stays, for the person who writes
   Python in a notebook in the browser. The design is in section 11 of
   `architecture.md`. The times of the table in 3.2 were taken under
   pyodide, and the direct build is expected to be near them, the same
   code on the same wasm engine; that has not been measured.

## 5. Open questions

- **faer against OpenBLAS on x86.** If faer is level there, as its
  benchmarks claim, the BLAS bindings matter mostly on macOS and the Linux
  wheels could ship faer alone, with no Fortran in CI. The spike's `zz` and
  `eigh` against numpy on the Linux box that runs the pyodide tests would
  settle it in a minute.
- **The mixed model fits at many individuals.** The logistic mixed model null
  fit inverts a individuals x individuals matrix a few dozen times, 25 s at 5000
  individuals; it is correct and matches GMMAT, and a Rust core with the same
  algorithm would not change that. A cheaper inner loop is an algorithmic
  question, not a language one.
- **Multiallelic variants.** pynei tests them as major allele against the
  rest, which keeps 0.3 of the signal when two common minor alleles differ
  in effect and all of it when they act alike. A per allele joint model,
  one row per variant and allele, is the design if the panels have several
  common alleles per site. Independent of the language.
- **Memory in the browser.** A kinship of 10000 individuals is 800 MB in float64
  and the projection matrix of the mixed model as much again; the 32 bit
  heap caps the browser at a few thousand individuals for the mixed models
  whatever the language, under pyodide and in the direct build of
  decision 9 alike.
- **zstd in the direct wasm build.** Closed on 20 September 2026: the
  vars file is compressed with lz4 in every build, which arrow-rs gets
  from `lz4_flex`, pure Rust. zstd was tried first. arrow-rs gets it from
  the `zstd` crate, which wraps the C library, and for
  `wasm32-unknown-unknown` it compiled and failed to link with Apple's
  clang, which does not emit wasm, and built and ran under node with the
  clang of the emscripten SDK in `CC_wasm32_unknown_unknown`. lz4 needs
  none of that, at the price of a file 2.2 times larger. The measurements
  are in `docs/specs/io_vars.md`.
- **Threads in the direct wasm build.** A page served with the two
  headers that isolate it from other origins can share memory between
  web workers, and the `wasm-bindgen-rayon` crate runs rayon on them. As
  far as is known in September 2026 it needs nightly Rust, and pyodide
  has no equivalent. If it works, the parser's gain from threads, 25x on
  one thread against 123x with them natively, would reach the browser.
  The first version is single threaded, and the serial versions that the
  wheel for pyodide needs are the ones it runs.

## 6. How to reproduce the numbers

- The Python benchmarks: `calc_kinship` and `calc_gwas` over variants from
  `test/var_generators.py`, and the tool comparison over a panel simulated
  by `test/gwas_reference/make_reference.py` with 1000 individuals and 100000
  variants, plink2 2.0 alpha 7 for arm64, R 4.6 with GMMAT and rrBLUP.
- The numpy in wasm numbers: pyodide 314.0.7 installed with npm under node
  26, numpy loaded from the jsdelivr CDN.
- The spike: `spike/README.md` in the pyNei repository has the native and
  the pyodide build steps.
- Issues #17 and #18 of the pyNei repository hold the tables with their
  dates.
