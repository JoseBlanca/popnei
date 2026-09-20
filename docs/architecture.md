# The architecture: how the variants flow

September 2026. The data flow model of popnei, decided before any code, and
the map of the modules. The objectives are in `objectives.md` and the
measurements that led here in `rust_core.md`.

## 1. The record level: a reader fills a variant the caller owns

The primary interface is the one Rust's own I/O uses, `BufRead::read_line`,
`csv::Reader::read_record`, noodles' `read_record`: the caller owns a
`Variant` and lends it to a reader, which fills it and says whether there
was one.

```rust
pub trait VariantReader {
    /// It fills `var` with the next variant. false when there are no more.
    fn read_variant(&mut self, var: &mut Variant) -> Result<bool>;
    fn samples(&self) -> &[String];
    fn ploidy(&self) -> usize;
    fn chroms(&self) -> &ChromTable;
    /// Which fields the caller wants filled. The rest may be skipped.
    fn set_needs(&mut self, needs: Needs);
}

pub struct Variant {
    pub chrom: u32,                 // an id into the reader's ChromTable
    pub pos: u64,
    pub gts: Vec<i8>,               // samples x ploidy, MISSING_ALLELE is -1
    pub id: String,
    pub alleles: Vec<String>,       // ref first, as in the VCF
    pub qual: Option<f32>,
    pub filled: Needs,              // what the reader actually filled
}
```

What this gives:

- **No allocations in the steady state.** The buffers inside the
  `Variant`, the genotype vector, the id, the alleles, are allocated once
  and reused, cleared and refilled, for a million variants.
- **Filters are readers over readers.** A variant filter holds its source
  and pulls from it until a variant passes, then hands that one back. A
  sample filter compacts the genotypes in place. Neither copies.
- **The reader owns its own buffers**, the line buffer, the batches it
  parses ahead, one set per thread when it fans out, and so does a
  writer. The lent `Variant` is the only thing that crosses.
- **Only what is asked for is filled.** `Needs` is a bit set, `GTS`,
  `CHROM_POS`, `ID`, `ALLELES`, `QUAL`. Most consumers ask for the
  genotypes alone. The VCF reader always fills chrom and pos, they cost
  nothing next to a thousand genotype fields, and skips the rest unless
  asked; the vars file reader does not decompress the columns nobody
  asked for. A consumer checks `filled` before it trusts a field.
- **Chromosome names are interned**, one `u32` per variant and one table
  per reader, not a string per variant.
- **A reader reads from any source of bytes, not from a path.** The VCF
  reader is generic over `BufRead`, and the vars file reader over
  `Read + Seek`, because an arrow file keeps the index of its batches at
  its end. Natively the source is a file. In a web application there is
  no filesystem, and the source is a file that the user picked in the
  page, or bytes in memory (section 11). The writers take any `Write` in
  the same way.

What it costs, accepted: a reader is not an `Iterator`, so the consumer
writes `while reader.read_variant(&mut var)? { ... }` and the adapters are
our own. And "one variant at a time" is the interface, not the memory
footprint of the parser, which parses ahead in batches.

## 2. The block level: for the calculations that want matrices

The kinship, the GWAS, the PCA, the LD and the distances between samples
want a block of variants as contiguous arrays, because a matrix product
over a block runs 5x to 10x faster than the same work variant by variant.
A `BlockCollector` builds blocks from any `VariantReader`, one memcpy of a
few kilobytes per variant, and a block native source, the vars file
reader with its record batches, hands its batches to the collector
without the copy.

```rust
pub struct Block {
    pub num_vars: usize, pub num_samples: usize, pub ploidy: usize,
    pub gts: Vec<i8>,                        // vars x samples x ploidy, C order
    pub chrom: Option<Vec<u32>>, pub pos: Option<Vec<u64>>,
    pub id: Option<Vec<String>>, pub alleles: Option<AllelesColumn>,
    pub qual: Option<Vec<f32>>,
}
impl Block {
    pub fn variant(&self, i: usize) -> VariantRef<'_>;    // slices, no allocation
    pub fn copy_variant_into(&self, i: usize, var: &mut Variant);  // back to a record
}
```

Blocks are owned and move through a pipeline; a stage that removes rows
compacts in place, and `reblock` restores the size. The block is the unit
of memory, a few thousand variants or about 5 million genotypes as in
pyNei, and the unit of parallelism, rayon over its rows. Stages that look
ahead, the LD filter, or that reorder, are block consumers, because a lent
variant cannot be held.

## 3. Threads

- Inside a reader or a writer: the VCF reader reads lines into a batch and
  parses them with rayon, then hands them out one by one; the writer
  formats a batch of lines in parallel and writes them in order. Each
  thread has its own buffers.
- Between stages: a read ahead thread between the reader and the
  consumer, one block or one batch ahead, as pyNei does.
- Inside a block consumer: rayon over the rows for the per variant work,
  the BLAS pool for the matrix products, never nested. A rayon worker that
  calls BLAS pins it to one thread; a big product is called from outside
  rayon.
- In wasm there are no threads, neither in the wheel for pyodide nor, in
  its first version, in the build for TypeScript. Everything above must
  build and run single threaded, rayon gated off the wasm targets.

## 4. The genotypes

An allele is one `i8`, 0 to 127, and -1 is missing; a genotype is
`ploidy` of them; there is no mask. It is the interchange format with
Python and the layout of the vars file. A 2 bit packed layout for biallelic
blocks is an option to measure, behind the same row views, so that no
consumer changes when the layout does. Multiallelic variants keep their
alleles; the calculations that collapse them to major against the rest say
so, as pyNei does.

## 5. The Python boundary

Python never reads one variant at a time. The Python `Variants` object
holds a reader, or a block source, and yields chunks as numpy arrays, the
genotypes as an int8 array of variants x samples x ploidy plus the
requested columns, which is what pyNei's chunks are. The binding crate
does that with a `BlockCollector` and hands the arrays over without
copying. Results are built in Python: the frozen dataclasses with pandas
frames and series, `pops` as a dict of name to samples, samples as tuples.
The Python layer is the API, the results and the tests, nothing else.
Under pyodide, the Python that runs in a browser tab, the same package
runs on the core built as a wasm wheel.
The other boundary, with TypeScript, is in section 11.

## 6. The vars file

The vars file is the contract between pyNei and popnei, and each has to
read what the other writes. Format 2.0, as `pynei/src/pynei/io_vars.py`
writes it: one arrow IPC file, feather v2, zstd, one record batch per
chunk; the schema metadata holds, under the key `pynei`, a json with
`var_format_version`, `samples`, `num_samples`, `ploidy` and
`num_vars_per_chunk`; each batch carries under `pynei_chunk` a json with
the range of chroms and positions it holds; the columns are `chrom`,
`pos`, `id`, `qual`, `alleles` as a list of strings per variant, and `gts`
as a fixed size list of `num_samples * ploidy` int8 per variant, whose
flat buffer is the genotype array itself. Written and read with arrow-rs.

## 7. Errors

`Result` everywhere, fail fast. A malformed VCF line is an error with the
line number and the field, not a warning and a skipped record. A vars
file of another major version is refused with the version in the message.

## 8. The crates and the layout

One repository, one cargo workspace, and two things built from it: the
wheel for Python and the wasm package for TypeScript, a package of npm,
the registry that JavaScript projects install from, with the core
compiled to WebAssembly. Its binding crate is written with wasm-bindgen,
which generates the JavaScript that calls Rust functions (section 11).

```
Cargo.toml                 the workspace
crates/popnei/             the core crate, pure Rust, no pyo3, cargo test
crates/popnei-python/      the Python binding crate, pyo3 + numpy, built by maturin
crates/popnei-js/          the JavaScript binding crate, wasm-bindgen
python/popnei/             the Python package: API, results, Variants
js/popnei/                 the TypeScript package: API, results, its tests
tests/                     the Python tests, pytest, pyNei as the oracle
tests/reference/           the reference data copied from pyNei and its script
docs/
pyproject.toml             maturin, manifest-path to the Python binding crate
```

pyNei is a development dependency of the Python side, a path dependency on
`../pynei`, and the tests run both libraries on the same inputs where they
overlap.

## 9. The modules of the core crate, and what of pyNei each one carries

| module | what it holds | pyNei functions it replaces |
|---|---|---|
| `variant` | `Variant`, `Needs`, `ChromTable`, `MISSING_ALLELE`, the row helpers: dosages, missing and het masks, allele counts of one row | `Genotypes.to_012`, `gt_counts` |
| `io::vcf` | the reader, parallel by batches, gzip; the writer | `vars_from_vcf`, and a writer pyNei does not have |
| `io::vars` | the arrow file reader, projection by `Needs`, batches for the collector; the writer | `load_vars`, `write_vars` |
| `filters` | readers over readers: missing data, maf, observed het, samples; the LD filter on blocks | `filter_by_missing_data`, `filter_by_maf`, `filter_by_obs_het`, `filter_samples`, `filter_by_ld_and_maf`, `gather_filtering_stats` |
| `block` | `Block`, `BlockCollector`, `reblock`, `VariantRef` | the chunks and `_resize_chunks` |
| `stats` | allele counts and frequencies per pop, per variant distributions with histograms, per sample stats, expected het, the polymorphism ratio | `calc_per_var_distribs`, `calc_per_sample_stats`, `diversity` |
| `dists` | Kosman between samples on blocks, Jost's D between pops | `calc_pairwise_kosman_dists`, `calc_jost_dest_pop_dists` |
| `linalg` | matrix product, symmetric eigendecomposition, Cholesky and solve, inverse, least squares; backends: BLAS and LAPACK natively, faer in wasm | numpy.linalg |
| `pca` | PCA of the 012 matrix, PCoA of a distance matrix | `do_pca_from_variants`, `do_pcoa_from_variants` |
| `ld` | Rogers Huff r2 between blocks of variants, by distance | `calc_rogers_huff_r2_matrix`, `iter_rogers_huff_r2`, `calc_ld_and_dist_per_pop` |
| `kinship` | the GRM, per pair denominators, principal components of it | `calc_kinship` |
| `gwas` | the four null models, the tests, the distributions erfc and betainc | `calc_gwas` |

The order of the work is the order of the rows: `variant`, `io::vcf`,
`io::vars`, `filters` and `block` are the walking skeleton, `stats` and
`dists` come next because they are what most users run, then `linalg` and
what sits on it.

## 10. The walking skeleton

The smallest path that exercises every layer once, and the first thing
built: the workspace and the two crates; `Variant`, `Needs` and the
`ChromTable`; the VCF reader, parallel, with gzip; the missing data
filter; the vars file writer; the `BlockCollector`; the Python `Variants`
over a reader giving numpy chunks; `vars_from_vcf`, `write_vars`,
`load_vars` and `filter_by_missing_data` in the Python package with pyNei's
signatures; and the tests: cargo tests of the reader and the filter, and
pytest tests that parse the reference VCFs with both libraries and compare
the chunks, that pyNei reads the vars file popnei writes, and that the
filter gives the same variants. On the TypeScript side it has the
JavaScript binding crate with the VCF reader over bytes in memory, the
missing data filter and the vars file writer, and a test under node that
parses a reference VCF, filters it and writes a vars file that pyNei
reads with the same variants. The skeleton is done when all those tests
pass, the Python binding crate builds as a wasm wheel with the steps in
pyNei's `spike/README.md`, and the wasm package builds.

## 11. The TypeScript boundary

A web application calls the core with no Python in the tab. None of what
follows has been built or measured yet: the trial build that
`rust_core.md` reports, its spike, made the wheel for pyodide only.

The JavaScript binding crate, `crates/popnei-js`, is to TypeScript what
the Python binding crate is to Python: it translates and holds no
calculation. It is written with wasm-bindgen, the Rust tool that
generates, from the Rust functions marked for export, the JavaScript that
calls them and their TypeScript declarations. It is compiled for the
target `wasm32-unknown-unknown`, WebAssembly with no operating system
under it: no files, no threads, no C library. The wheel for pyodide is
compiled for another target, `wasm32-unknown-emscripten`, where
emscripten emulates all three. So the core is built for two wasm targets,
and a dependency of the core has to build for both. One written in C
builds under emscripten with its compiler and may not build for the
direct target; the zstd compression of the vars file is the known case,
an open question of `rust_core.md`.

The TypeScript package, `js/popnei`, sits on the binding as the Python
package does: the functions with the names and the arguments of the
Python API, in camelCase, and the result objects. It is published to npm
with the compiled core inside, and that is the wasm package.

- **The calculations run in a web worker**, a thread of the page that
  cannot touch what the page shows. A calculation of seconds on the main
  thread of the page would freeze it for those seconds.
- **A file that the user picked is read as the reader asks for it.** In a
  web worker `FileReaderSync` reads a range of bytes of a file and
  returns when it has them, which is what `Read` and `Seek` need, so the
  binding wraps it as a source of bytes and a VCF larger than the memory
  of the tab still streams. `FileReaderSync` does not exist outside a
  worker. The other source is bytes already in memory, a `Uint8Array`,
  for small files and for the tests under node.
- **TypeScript never reads one variant at a time**, as Python does not.
  It calls a calculation over a source, and the loop over the variants
  runs inside wasm. An application that wants the genotypes gets blocks,
  the genotypes as an `Int8Array` of variants x samples x ploidy with the
  requested columns.
- **Results cross as typed arrays, copied.** A `Float64Array` that is a
  view into the memory of wasm stops being valid when that memory grows,
  so the binding copies each result out. The results are per variant, per
  sample or samples x samples, so the copy is small next to the
  calculation. Sample and population names cross as arrays of strings,
  `pops` as an object of population name to sample names, and the result
  objects are built in TypeScript.
- **A writer writes into memory.** The vars file writer fills a buffer
  that the page offers as a download. A file that does not fit in memory
  would need the private filesystem that the browser gives each site,
  which a worker can write synchronously; that is left until an
  application needs it.
- **An error of the core is thrown as a JavaScript `Error`** with the
  message it has in Rust, from one place in the binding crate.
- **The tests** of the TypeScript package run under node, the JavaScript
  runtime outside the browser, on the reference files of
  `tests/reference/`. They expect the same numbers as the Python tests,
  written into them as literals. The calculations themselves are tested
  once, in the core crate.
