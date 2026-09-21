# The architecture: how the variants flow

September 2026. The data flow model of popnei and the map of the modules.
It was decided before any code, and revised on 20 September 2026: the
variants flow in blocks from the source to the calculation, where the first
version had a single variant that each reader filled. What was revised and
why is at the end of section 1. The objectives are in `objectives.md` and
the measurements that led here in `rust_core.md`.

## 1. A reader gives blocks

A block is a run of consecutive variants held as arrays, the genotypes of
all of them in one array and each other field in a column (section 2).
Everything that gives variants gives them in blocks, through one trait:
the VCF reader, the vars file reader, and a filter, which is a reader over
another reader.

```rust
pub trait BlockReader: Send {
    /// The next block, or None when there are no more.
    fn next_block(&mut self) -> Result<Option<Block>>;
    fn individuals(&self) -> &[String];
    fn ploidy(&self) -> usize;
    fn chroms(&self) -> &ChromTable;
    /// Which fields the caller wants filled. The rest may be skipped.
    fn set_needs(&mut self, needs: Needs);
}
```

What this gives:

- **One shape from the source to the calculation.** Both sources are
  blocks inside. The vars file is a sequence of record batches, and a
  batch is a block. The VCF reader reads the lines of a block and parses
  them on the threads of rayon, each thread writing its row straight into
  the arrays of the block, which is how the spike of `rust_core.md` got
  its 0.55 s on one thread and 0.11 s on 18 cores. Nothing is taken apart
  into single variants and put together again.
- **The threads have the variants they need.** The per variant work, the
  counts, the masks, the filters, the statistics, runs with rayon over the
  rows of the block in hand. The matrix work, the kinship, the GWAS, the
  PCA, the LD, takes the same block as a matrix, because a matrix product
  over a block runs 5x to 10x faster than the same work variant by
  variant.
- **Filters are readers over readers.** A variant filter takes a block
  from its source, decides which rows stay, with rayon over the rows,
  compacts every column in place and gives the block on. A filter of
  individuals compacts the genotypes of each row in place. Neither
  allocates a block. A block left with no variants is not given; the
  filter takes the next.
- **One variant is a view into a block.** A calculation that works variant
  by variant loops over `block.variants()`, which are slices and allocate
  nothing, and the row helpers, the dosages, the masks, the allele counts
  of one variant, work on that view. What such a calculation keeps from
  one block to the next does not grow with the variants: for a mean, a sum
  and a count, to which it adds what rayon reduced over the rows of each
  block.
- **Only what is asked for is filled.** `Needs` is a bit set, `GTS`,
  `CHROM_POS`, `ID`, `ALLELES`, `QUAL`, and a column that nobody asked for
  is `None` in the block. Most consumers ask for the genotypes alone. The
  VCF reader then does not parse the other columns, and the vars file
  reader does not decompress them. A change of `Needs` holds from the next
  block that the reader builds, so a block that was read ahead keeps the
  columns it was built with.
- **Chromosome names are interned**, one `u32` per variant and one table
  per reader, not a string per variant. A filter gives the table of its
  source.
- **A reader reads from any source of bytes, not from a path.** The VCF
  reader is generic over `BufRead`, and the vars file reader over
  `Read + Seek`, because an arrow file keeps the index of its batches at
  its end. Natively the source is a file. In a web application there is
  no filesystem, and the source is a file that the user picked in the
  page, or bytes in memory (section 11). The writers take blocks and any
  `Write` in the same way.

What it costs, accepted:

- A filter leaves blocks of uneven size, and the consumers that care about
  the size put `reblock` before them (section 2).
- An error in one variant loses its block: the blocks before it are given,
  and then the error. With a VCF that has a wrong line after 250 good ones
  and blocks of 100, the caller gets two blocks and the error.
- The memory of a stage is a block, about 10 MB of genotypes at ploidy 2,
  and with the read ahead of section 3 two or three are alive at a time.

What was revised. The first version of this document had two levels: a
record level, where the caller owned one `Variant` and lent it to a reader
that filled it, `read_variant`, as `BufRead::read_line` fills a string,
and a block level for the matrix work, with a collector that copied the
variants of a reader into a block. The owner dropped the record level on
20 September 2026, for three reasons. A consumer that holds one lent
variant cannot run rayon across variants, so the statistics and the
filters needed a block anyway to use the threads that the objectives ask
for. Both sources are blocks inside, so the single variant was a narrow
point between two arrays, undone on each side. And the file that users
will usually read is the vars file, not the VCF that the record level was
thought for. The copies were not the reason: on the vars file of 20000
variants and 1000 individuals of `docs/specs/io_vars.md`, decompressing
the genotypes took arrow-rs 18.8 ms, and copying them variant by variant
into a record and from there into a block brought it to 19.6 ms, on one
thread of the owner's Apple M5 Pro.

## 2. The block

```rust
pub struct Block {
    pub num_vars: usize, pub num_individuals: usize, pub ploidy: usize,
    pub gts: Vec<i8>,                        // vars x individuals x ploidy, C order
    pub chrom: Option<Vec<u32>>, pub pos: Option<Vec<u64>>,
    pub id: Option<Vec<String>>, pub alleles: Option<AllelesColumn>,
    pub qual: Option<Vec<f32>>,
}
impl Block {
    pub fn variants(&self) -> impl Iterator<Item = VariantRef<'_>>;  // slices, no allocation
    pub fn retain_vars(&mut self, keep: &[bool]) -> Result<()>;      // a filter compacts in place
}
```

Blocks are owned, and a reader gives each one away. That is what lets the
binding crate hand the genotype array to numpy without copying it, a read
ahead thread move a block to the thread that consumes it, and a filter
compact it in place. What a new allocation for each block costs was
measured on that same vars file: copying each decompressed batch, 10 MB,
into a vector allocated for it added 0.4 ms to the 18.8 ms of the
decompression of its four batches. A consumer that is done with a block
can give it back to be refilled, if a measurement ever asks for it.

A block holds about 5 million genotypes, and no fewer than 100 variants
and no more than 10000, the numbers of pyNei's chunks, which nobody has
measured for popnei. A source gives blocks of the size it is asked for,
and the last one is the only one that can be shorter. The vars file
reader gives the batches of its file as they are. `reblock` is a reader
over a reader that cuts and joins the blocks of its source to a size, and
it goes before the consumers that care: the matrix work, after a filter
took rows out; the vars file writer, which writes one batch per block; and
the end of `iter_blocks`, so that a user gets blocks of one size. The per variant
calculations take the blocks as they come.

The block is the unit of memory and the unit of parallelism, rayon over
its rows. Stages that look ahead, the LD filter, or that reorder, hold
more than one block.

## 3. Threads

- Inside a reader or a writer: the VCF reader reads the lines of a block
  and parses them with rayon, each thread writing its own row of the
  block; the VCF writer formats the rows of a block in parallel and writes
  them in order. Each thread has its own buffers.
- Between stages: a read ahead thread between the reader and the
  consumer, one block ahead, as pyNei does.
- Inside a consumer, a filter or a calculation: rayon over the rows of a
  block for the per variant work, the BLAS pool for the matrix products,
  never nested. A rayon worker that calls BLAS pins it to one thread; a
  big product is called from outside rayon.
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

The Python `Variants` object holds a reader and has no genotypes of its
own. A calculation takes it and runs its loop over the blocks of that
reader inside the core, so no calculation pays for a call from Python per
block or per variant. This departs from pyNei, whose `Variants` yields
chunks, arrays of a few thousand variants, to calculations that are
written over them in Python. A user holds neither variants nor an
iterator of them. The one way
genotypes come out is `Variants.iter_blocks(fields=...)`, which gives
blocks, the genotypes as an int8 array of variants x individuals x ploidy
that the binding crate hands to numpy without copying, with the columns
that were asked for. It is for the user who wants the genotypes for an
analysis of their own, and for the tests, which compare those blocks with
pyNei's chunks. The function that gives a `Variants` from a VCF is
`open_vcf`, because nothing is read when it is called but the header.
Results are built in Python: the frozen dataclasses with pandas
frames and series, `pops` as a dict of name to individuals, the names of the individuals as
tuples.
The Python layer is the API, the results and the tests, nothing else.
Under pyodide, the Python that runs in a browser tab, the same package
runs on the core built as a wasm wheel.
The other boundary, with TypeScript, is in section 11.

## 6. The vars file

The vars file is popnei's own file of variants, and its format is in
`docs/specs/io_vars.md`. The owner decided on 20 September 2026 that it
owes nothing to the vars file of pyNei, which it took its shape from:
pyNei will not be used once popnei exists, and neither library reads the
files of the other. It is one arrow IPC file, feather v2, compressed with
lz4, which is pure Rust in arrow-rs where zstd is C, with one record batch
per block; the schema metadata holds, under the key `popnei`, a json with
`format_version`, `individuals`, `ploidy` and `num_vars_per_block`; the
footer holds, under `popnei_batches`, the number of variants of each batch
and, for each chromosome in it, the smallest and the largest position, so
that a reader asked for a region skips the batches outside it; the columns
are `chrom`, `pos`, `id`, `qual`, `alleles` as a list of strings per
variant, and `gts` as a fixed size list of `num_individuals * ploidy` int8
per variant, whose flat buffer is the genotype array itself. Any program
with an arrow library opens it as a table. Written and read with arrow-rs.

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

pyNei is a development dependency of the Python side, taken from
`https://github.com/JoseBlanca/pynei` at the commit that `[tool.uv.sources]`
of `pyproject.toml` names, and the tests run both libraries on the same
inputs where they overlap.

## 9. The modules of the core crate, and what of pyNei each one carries

| module | what it holds | pyNei functions it replaces |
|---|---|---|
| `variant` | `Needs`, `ChromTable`, `MISSING_ALLELE`, `VariantRef`, the view of one variant of a block, and the row helpers over it: dosages, missing and het masks, allele counts | `Genotypes.to_012`, `gt_counts` |
| `io::vcf` | the reader, which parses the lines of a block in parallel, gzip; the writer | `vars_from_vcf`, and a writer pyNei does not have |
| `io::bgzf` | the reader of the members of a file that bgzip wrote, which `io::vcf` reads such a source through: it cuts each member by the size the member states and checks it | none; pyNei reads a bgzipped VCF with Python's `gzip` |
| `io::vars` | the arrow file reader, projection by `Needs`, a batch of the file as a block; the writer; a format of popnei's own | `load_vars`, `write_vars` |
| `filters` | readers over readers, which compact the blocks in place: missing data, maf, observed het, individuals; the LD filter | `filter_by_missing_data`, `filter_by_maf`, `filter_by_obs_het`, `filter_samples`, `filter_by_ld_and_maf`, `gather_filtering_stats` |
| `block` | `Block`, the `BlockReader` trait, `AllelesColumn`, `reblock` | the chunks and `_resize_chunks` |
| `stats` | allele counts and frequencies per pop, per variant distributions with histograms, per individual stats, expected het, the polymorphism ratio | `calc_per_var_distribs`, `calc_per_sample_stats`, `diversity` |
| `dists` | Kosman between individuals on blocks, Jost's D between pops | `calc_pairwise_kosman_dists`, `calc_jost_dest_pop_dists` |
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
built: the workspace and the two crates; `Block`, `BlockReader`, `Needs`
and the `ChromTable`; the VCF reader, parallel, with gzip; the missing
data filter; the vars file writer and reader; `reblock`; the Python
`Variants` over a reader, with `iter_blocks`; `open_vcf`, `write_vars`,
`open_vars` and `filter_by_missing_data` in the Python package, with
pyNei's signatures where the specs keep them; and the tests: cargo tests of
the reader and the filter, and pytest tests that parse the reference VCFs
with both libraries and compare popnei's blocks with pyNei's chunks, that
pyarrow opens the vars file popnei writes and finds the variants of the
VCF in it, that the file read back gives the blocks of the VCF, and that
the filter gives the same variants. On the TypeScript side it has the
JavaScript binding crate with the VCF reader over bytes in memory, the
missing data filter and the vars file writer and reader, and a test under
node that parses a reference VCF, filters it, writes a vars file and reads
it back with the same variants. The skeleton is done when all those tests
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
direct target; zstd was the known case, and it is why the vars file is
compressed with lz4, which is pure Rust (`docs/specs/io_vars.md`).

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
- **The loop of a calculation runs inside wasm**, as it runs inside the
  core for Python. TypeScript calls a calculation over a source. An
  application that wants the genotypes gets blocks from `iterBlocks`,
  the genotypes as an `Int8Array` of variants x individuals x ploidy
  with the fields that were asked for.
- **Results cross as typed arrays, copied.** A `Float64Array` that is a
  view into the memory of wasm stops being valid when that memory grows,
  so the binding copies each result out. The results are per variant, per
  individual or individuals x individuals, so the copy is small next to the
  calculation. The names of individuals and populations cross as arrays
  of strings, `pops` as an object of population name to the names of its
  individuals, and the result
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
