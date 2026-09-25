# popnei in TypeScript

The TypeScript package of popnei, a population genetics library whose
calculations are written in Rust: the core crate compiled to WebAssembly,
the code that a browser or node calls it through, and the functions and
the result objects an application uses. What the package exports today is
`init`, which loads the WebAssembly and is awaited before anything else is
called, and `version`, the version of the core crate; `openVcf`, which
reads the header of a VCF, held as bytes or as the file the user picked in
the page, and gives a `Variants`, the handle whose `iterBlocks` gives the
genotypes block by block; `writeVars`, which
gives back the bytes of a vars file with every variant of a `Variants`, and
`openVars`, which opens such a file, its bytes or the file of the page, as
another `Variants`. A vars file is
one arrow IPC file, also called feather v2, which pandas, R and polars open
as a table with no popnei installed: it is where a user keeps their
variants once the VCF has been read. `numPassesOf` and the `onProgress` of
a `Variants` are for an application that shows how far a calculation has
got and lets its user stop it, and "A file of the page, in a web worker"
below has both.

Ten calculations read the variants of a `Variants`.
`calcPairwiseKosmanDists` gives, in a `Distances`, the Kosman distance of
every pair of individuals, how many alleles the two do not share at a
variant averaged over the variants at which both were called, which runs
from 0 for two individuals with the same genotype everywhere to 1 for two
that share no allele anywhere, and which is what a tree or a principal
coordinate analysis of individuals is built from.
`calcRogersHuffR2Matrix` gives, in an `R2Matrix` with the chromosome and
the position of each variant, r² for every pair of variants, the square of
the correlation between their dosages over the individuals called at both,
where the dosage of a genotype is how many of its alleles are not the major
allele of its variant: it is 1 when the dosage at one variant fixes the
dosage at the other and 0 when knowing one says nothing about the other.
`doPcaFromVariants` gives, in a `VariantsPcaResult`, the principal
components of the individuals over those same dosages, where each of them
is a direction along which the individuals differ most: the projection of
every individual on the first ones, how much of the variance each holds and
the weight every variant has in them. `calcPerVarDistribs` gives, for each
population a user names in `pops` and each of five statistics of a variant,
the mean over the variants that had a value and a histogram of them. The
five are the observed heterozygosity, the major allele frequency, the
expected heterozygosity, plain and unbiased, and the polymorphism ratio,
which is three counts and two ratios per population and not a distribution.
`calcPerIndividualStats` gives two numbers for each individual instead of
one for each population: the share of the variants at which its genotype is
missing, `missingGtRate`, and the share of its called genotypes at which it
is heterozygous, `obsHetRate`. The second says which individuals are more
heterozygous than the rest, a sign of a mixed sample or of an outcrossed
individual among inbred ones, and it is NaN for an individual that called
no genotype. `calcPopDiversity` gives, in a `PopDiversity` and for each
population a user names in `pops`, how much variety it holds: how many
alleles its individuals called, how many of those no other population
called at the same variant, how many of the variants vary in it, and F_IS,
how far its genotypes are from the proportions its allele frequencies would
give if its individuals paired at random. Given a `numCalledAlleles` it
gives those first three again standardized to a draw of that many called
alleles, so that a population of 20 individuals and one of 200 can be
compared, and the folded site frequency spectrum of the draw, how the
variants of the population are spread over the count of their rarer allele.
`calcLdAndDistPerPop` gives, for each population, how r² falls off as the
two variants of a pair move apart on a chromosome: the pairs are put into
bins of distance, each bin carrying how many pairs it holds and the mean
and the standard deviation of their r², and beside the bins a curve fitted
to every pair of the population, which carries the distance at which r² has
fallen to half. How fast it falls is a property of the population, two
variants that sit close together having had fewer recombinations between
them than two that sit far apart.

The other three read the same variants. `calcPopDists` gives how far apart
every pair of the populations a user names is: two populations are far apart
when the alleles of their individuals are not the same alleles in the same
proportions, and there are seven measures of how far, Hudson's F_ST, f_2,
the chord distance, Nei's D_A, Jost's D, Nei's G_ST and the standardized
G''_ST, each answering a different question. One pass gives the ones that
were asked for, in a `PopDists` with a `Distances` for each measure and the
standard error of every pair beside its value. `calcKinship` gives, in a
`Kinship`, how much more of their genome each pair of individuals shares
than two individuals drawn at random from the same panel do: it is the
matrix of VanRaden 2008, which plink2's `--make-rel` computes, where an
entry off the diagonal is about 0.5 for full sibs or for a parent and a
child and near 0 for two individuals with no recent ancestor in common, and
an entry on the diagonal is 1 plus the inbreeding of that individual.
`calcGwas` says, in a `GwasResult`, which variants are associated with a
trait the user hands it as one number for each individual: it tests one
variant at a time over the dosages of the individuals, with the covariates
the user gives, and, when it is given a `Kinship`, with the relatedness of
the panel as a random effect, so that a variant which only marks the
ancestry of the panel does not look associated. `doPca`, which is not a
consumer of a `Variants`, gives the components of a table of individuals and
traits handed to it as numbers, which is the same analysis over values an
application holds and not over a source of variants.

Each of the twelve consumers of a `Variants`, `iterBlocks`, `writeVars`,
`calcPairwiseKosmanDists`, `calcPopDists`, `calcPopDiversity`,
`calcRogersHuffR2Matrix`, `calcLdAndDistPerPop`, `calcKinship`,
`doPcaFromVariants`, `calcGwas`, `calcPerVarDistribs` and
`calcPerIndividualStats`, gives back the counts of the pass it made over the
source, in a `passStats`: how many variants it
took, and how many each filter of the `Variants` was given and kept. A
filter is a step, a method of the `Variants` that `steps` then lists, and
there are five of them. Four take variants out: `filterByMissingData`,
which keeps the variants whose missing genotypes divided by all the
individuals are at most the threshold it is given, `filterByMaf`, over the
count of the commonest allele of a variant divided by its called alleles,
`filterByObsHet`, over its heterozygous genotypes divided by its called
ones, and `filterByLd`, which keeps the variants whose r² against every
variant kept within a window behind them on their chromosome is at most the
threshold, so that what is left does not repeat what a variant near it
already said. The fifth, `filterIndividuals`, keeps individuals and not
variants: it takes the genotypes of the individuals a user names, at every
variant, in the order they named them, and after it `individuals` and
`numIndividuals` are the kept ones.

Section 11 of `docs/architecture.md` has the design, `crates/popnei-js` is
the binding crate, the Rust that is compiled to WebAssembly and that holds
no calculation of its own, and `docs/specs/io_vcf.md`,
`docs/specs/io_vars.md`, `docs/specs/block.md`, `docs/specs/variant.md`,
`docs/specs/filters.md`, `docs/specs/dists.md`, `docs/specs/pca.md`,
`docs/specs/kinship.md`, `docs/specs/gwas.md`, `docs/specs/ld.md`,
`docs/specs/diversity.md` and `docs/specs/stats.md` say what they give.
`docs/specs/js_sources.md` has what a source of this package is read from and
what it tells the page while it reads.

## Building it

One command, from a clean checkout of the repository:

    cd js/popnei && npm run build

It runs four things: `npm install`, which brings the TypeScript compiler
and the type declarations of node; `cargo build --package popnei-js
--release --target wasm32-unknown-unknown`, which compiles the core and
the binding crate to WebAssembly; the `wasm-bindgen` command line, which
reads that wasm file and writes `js/popnei/wasm/`, the JavaScript that
calls into the WebAssembly and the TypeScript declaration of every
function exported from Rust; and the TypeScript compiler, which writes
`js/popnei/dist/` from `src/` and checks the tests against it. Neither
`wasm/` nor `dist/` nor `node_modules/` is in git.

The `wasm-bindgen` command line is given `--remove-name-section`, which
takes out of the wasm file the section that holds the name of every function
of it: `js/popnei/wasm/popnei_bg.wasm` is 2021550 bytes with the flag and
2621987 bytes without, 600437 bytes of names that every user of the package
downloads. Both are of the build of 24 September 2026, on macOS on aarch64,
and they grow with the code of the crates. What they are for is the stack of
a trap, a panic of Rust among the causes, which with the flag names the
functions by their number and without it by their name. To read one, build
again without the flag and make the trap happen there.

The version of the `wasm-bindgen` crate, in the `Cargo.toml` of the
workspace, has to be the version of the `wasm-bindgen` command line that
is installed, so it is pinned there, `=0.2.128`. The crate writes into the
wasm file what the command line then reads, and when the two versions
differ the command line stops and names both. A new command line,
`cargo install wasm-bindgen-cli`, needs that line of the workspace
manifest changed to its version.

## What crosses between Rust and JavaScript

What the binding crate exports and what the generated JavaScript and its
declarations then hold, as it was found with wasm-bindgen 0.2.128:

- The names are the Rust ones. A method `next_block` is `next_block` in
  JavaScript, so the camelCase of the API is the package's doing and not
  the generator's.
- `Vec<i8>` arrives as an `Int8Array`, `Vec<f64>` as a `Float64Array`,
  `Vec<u32>` as a `Uint32Array` and `Vec<String>` as an array of strings.
  Each of them is copied out of the memory of the WebAssembly, and so is a
  `Uint8Array` that goes the other way, into it.
- `Option<T>` is `T | undefined`, and `undefined` is what the package turns
  into the `null` that the specs give the columns a block does not hold.
- A `Vec<Vec<String>>` does not compile: "the trait bound `String:
  ErasableGeneric` is not satisfied". An array of arrays needs the `js-sys`
  crate to build it, which this crate does not depend on, so the alleles of
  a block cross as one array of texts with the number of alleles of each
  variant beside it, which is how the core holds them, and the package cuts
  the one with the other.
- An argument of a type this crate exports crosses by value, `Option<T>`
  included, and `Option<&T>` of one does not compile: "the trait bound
  `&Steps: OptionFromWasmAbi` is not satisfied". The generated JavaScript
  takes the pointer out of the object it is given and leaves that object
  dead, so a `Variants` hands each pass a copy of its steps, which is what
  a pass runs anyway: a step added while a pass runs holds from the next
  one.
- A `Result<T, E>` is thrown when `E` is `Into<JsValue>`, and
  `JsError::new(message)` is the JavaScript `Error` a user catches. The
  orphan rule keeps `From<popnei::Error> for JsValue` out of this crate, so
  the errors of the core go through a type of the crate, as they do in the
  Python binding crate.
- Every exported struct gets a `free()` and a `[Symbol.dispose]()` in its
  declaration, and the generated JavaScript registers each object in a
  `FinalizationRegistry`, which frees one that was dropped without a
  `free()` when the garbage collector reaches it.
- `#[wasm_bindgen]` on a `pub const` does not compile: "will not work on
  constants unless you are defining a
  `#[wasm_bindgen(typescript_custom_section)]`". So the defaults of the
  API, the ploidy of 2 and the filter of `docs/specs/io_vcf.md` and the
  five of `calcPerVarDistribs` among them, cross as functions that return
  the constants of the core.
- A number of JavaScript that goes in as a `usize` is a float64 turned
  into an integer of 32 bits with no error: the fraction is thrown away
  and what is left is kept modulo 2^32. A ploidy of 2.5 and one of
  2^32 + 2 both arrived as 2, a size of 2^32 + 1 as blocks of one variant,
  -1 as 4294967295 and NaN as 0. An argument that is not a `Uint8Array`
  where bytes are asked for is read as whatever its memory holds, and one
  that is not an object of the right class throws a `TypeError` of the
  generated code that names none of the two. So the package checks every
  argument before the call, in `src/arguments.ts`, as the Python package
  leaves pyo3 to do.
- An argument of bytes costs a copy of the file going in: wasm-bindgen
  allocates a `Vec<u8>` of the length of the `Uint8Array` inside the
  memory of wasm and copies it there. That memory grows and never shrinks,
  so a second copy of the same file stays for as long as the page lives,
  which is why `open_vcf` and `open_vars` keep the `Vec` they were given
  and share it with every pass: an 80 MB VCF costs 80 MB of the tab and
  not 160 MB. The measurements that follow are of one file, written by
  `crates/popnei/benches/make_big_vcf.py` with its `NUM_VARS` at 20000:
  a VCF of 80692954 bytes, 20000 variants of 1000 diploid individuals,
  whose vars file in batches of 1000 is 19185674 bytes. Each of them was
  made in a process of its own, because the memory of wasm never shrinks
  and what one measurement frees is room the next one does not have to
  grow for, and each is the memory of wasm before the call against after
  it.
- `openVars` of that file grows the memory of wasm by 18.4 MB, the file
  and nothing else. A pass over it grows it by what the blocks it builds
  hold: 11.7 MB with `numVarsPerBlock` 1000, 39.2 MB with none, which for
  1000 individuals is blocks of 5000 variants, and 62.6 MB with 10000. A
  second pass grows it by nothing, whichever of the three, because the
  first one left the room behind.
- A `Vec<u8>` coming back is a copy going out, so a file that is written
  is in the memory of wasm and in the `Uint8Array` at once. While
  `writeVars` runs, that memory holds the source, the block being read and
  written, and the file that is growing, and the package reads that file
  out of it in pieces of 1 MiB, each freed there as it is copied into the
  array the user gets. Writing the file above from its VCF grows the
  memory of wasm by 30.8 MB with batches of 1000 variants, beyond the
  77.0 MB of the source, and by 20.8 MB with batches of 100, which is a
  file of 19507162 bytes. Writing it again from the vars file it came
  from, in batches of 1000, grows it by 34.2 MB. The size of the batches
  is what an application that runs out of memory lowers, and it is the
  size of the block that is read as well as the size of the batch that is
  written.
- A `Blob`, which the `File` of a page is one of, crosses as a handle and
  costs no copy: what goes into the memory of wasm is the number of the
  entry of a table of the binding crate that holds the `Blob`, the
  `FileReaderSync` that reads its ranges and the function the page is told
  the progress with. They stay in JavaScript because a reader of the core
  has to be `Send`, movable to another thread, which no handle of JavaScript
  is. Each range popnei reads crosses once, copied out of the `ArrayBuffer`
  that `FileReaderSync` fills.
- A panic of Rust in wasm is a trap: the call ends where it is, the memory
  of wasm keeps what it held, and an object that was borrowed at that
  moment stays borrowed, so a later `free()` of it throws "attempted to
  take ownership of Rust value while it was borrowed" instead of giving
  the memory back. That was seen in the review of this work package, over
  the real error of a block the memory could not hold, which the core now
  refuses with an error and not an abort. What is left of the instance
  after a trap has not been measured here. A vars file that was damaged
  after it was written is what can still reach one: the core checks the
  message of every batch before arrow-rs sees it, and a sweep of damaged
  files, in "How it is verified" of the reader of `docs/specs/io_vars.md`,
  found 2783 of 1299990 that reach two asserts inside arrow-rs all the
  same. Natively a `catch_unwind` turns those into an error, and in wasm
  nothing catches them.
- The `finally` of a generator does not run when the generator was never
  started, and `free()` of an object of wasm that is still borrowed throws
  "attempted to take ownership of Rust value while it was borrowed", which
  in a `finally` hides the error on its way out. The iteration of blocks
  counts its pass inside the generator for the first, and throws what the
  free says only when nothing else is being thrown, for the second.

## The timing

`bench/time_pca.mjs` is not a test and `npm test` does not run it: it times
`doPcaFromVariants` under node over the bytes of a vars file, on the `wasm/`
that is there, and it is what task 4.2 of `docs/plans/pca.md` measured this
package with. `docs/reports/pca-measurement.md` has its numbers and the
files it read them on.

    node bench/time_pca.mjs <path to a vars file> [--runs n] [--num-prin-comps n]

## The tests

    npm test

It runs the TypeScript compiler over `src/` and over `test/`, which writes
`dist/` again and type checks the tests against it, and then the test
runner of node itself, `node --test`. What it does not build is the
WebAssembly: the tests run the `wasm/` that is there, so a change of Rust
is tested only after `npm run build`. The tests import the name of the
package, `popnei`, which node
resolves to the built entry point of node, and read the reference VCFs of
`tests/reference/vcf/` at the root of the repository, the files the Python
tests read. They assert that the version the package gives and the version
in `package.json` are both the one of `[workspace.package]` of the
`Cargo.toml` of the repository, that `init` loads the WebAssembly once,
that a function called before `init` was awaited throws an `Error` that
says so, that the entry point of a page answers with the WebAssembly it
fetches and reads a VCF through it, and that the blocks of `cases.vcf` and
`differences.vcf` hold the variants of the tables of
`docs/specs/io_vcf.md`, with the default and with `onlyPassed` false.
`test/vars.test.ts` writes those four variants into a vars file with
`writeVars` and reads them back with `openVars`, and reads
`tests/reference/vars/zstd.vars`, the file compressed with zstd that
popnei cannot write and refuses at its first block.
`test/pass_stats.test.ts` reads the 500 variants of `many.vcf`, and the
vars file written from them, and asserts the counts of the pass that
`iterBlocks` and `writeVars` give: 500 after a whole pass, 21 after three
blocks of 7, none for a source with no variant, and, for a pass that ended
at a wrong line, the variants of the blocks the user got and not the ones
the file holds. `test/filters.test.ts` puts the three filters on the same
file and asserts the numbers of the table of `docs/specs/filters.md`, which
are bcftools 1.24's and pyNei's: 26 variants kept by the missing data filter
at 0, 35 by the maf filter at 0.5 and 22 by the observed heterozygosity one
at 0.1, each with the first five positions it keeps, and 106 by the three of
them chained at 0.04, 0.8 and 0.5, which count 500 and 215, 215 and 163, and
163 and 106. `test/filter_individuals.test.ts` keeps `ind05`, `ind00` and
`ind49` of that file, in that order, which is what `bcftools view -s
ind05,ind00,ind49` gives, and asserts the 500 variants that come out with
their genotypes, the 423 that the missing data filter at 0 after the step
keeps where the same filter over the 50 individuals keeps 26, and the
`Error` of a name the file does not have, of a name given twice, of no name
and of a second filter of individuals. `test/stats.test.ts` reads the panel of
`tests/reference/stats/`, 1200 variants of 200 individuals in the three
populations of `panel_pops_bcftools.txt`, and asserts through
`calcPerVarDistribs` the literals of `p0`, the population of `s000`, that
`docs/specs/stats.md` gives. Over the whole file: 1112 of its 1200 variants
are polymorphic in `p0`, 1173 are variable and 1200 have a major allele
frequency there, and the two ratios are 0.926666666667 and
0.947996589940. Over its first variant alone, `var0000`, which the test
writes as a VCF of its own so that the mean of the pass is the value of
that variant: an observed heterozygosity of 0.166667, a major allele
frequency of 0.895833 and a plain expected heterozygosity of 0.186632, each
within 1e-6 of what plink2 prints, and the unbiased expected heterozygosity
of `p1`, 0.502750, which is the only population a number is known for. A
name that is not an individual of the pass and a key of `histKwargs` that
popnei does not know are each an `Error` there. Through
`calcPerIndividualStats` the same file asserts the two rates of `s000`, 34
missing genotypes of 1200 variants and 426 heterozygous of 1166 called, and
of `s001`, 44 and 397 of 1156, and those of `ind00` and `ind01` of
`many.vcf`, 29 of 500 and 201 of 471 and 25 and 195 of 475, which the same
plink2 reports give with `--vcf-half-call m`; the NaN of an individual that
called no genotype; the names coming in the order a `filterIndividuals`
named them in; and the `Error` of a source with no variant and of steps that
kept none. The comparison with pyNei itself is
the one of the Python tests; node runs neither library. Several of the tests
watch the memory of the WebAssembly, which they reach through the loader
`wasm/popnei.js` generates: that a block, and the bytes of a vars file,
kept while enough more is read for that memory to grow still hold what
they held, and that an iteration gives its pass back however it ends.
`test/vars_memory.test.ts` is alone in its process for two more, because
what a test frees stays in that memory as room the next one fits into: it
writes a vars file of 12 MB and asks that the write stay under twice the
file, and then opens twelve passes over it at once and asks that they
grow the memory by less than one copy of it, which a reader that copied
the bytes for each pass would not.
`test/blob_outside_a_worker.test.ts` asserts the one thing about the
reading of a file of the page that node can: node has `Blob` and `File` and
no `FileReaderSync`, which is the case of the main thread of a page, so
`openVcf` and `openVars` of a `Blob` and of a `File` are refused there with
a message that names the reader and the worker, and the same bytes open.

What a range of a real file gives is asserted in a browser:

    npm run test:browser

It runs Playwright, which drives a browser from a script, over Chromium
headless. A page served from the repository starts a module web worker, the
worker loads the WebAssembly of `wasm/` through `dist/web.js`, and the
assertions are made in there, on the literals the tests under node assert.
Chromium is downloaded once, with `npx playwright install chromium`, and
Playwright says so when it is missing. Firefox and WebKit are not run: the
owner decided on 24 September 2026 to start with Chromium and to add them
when an application needs them.

## Where it runs

The package needs the vector instructions of WebAssembly, the ones that
work on sixteen bytes at a time, which popnei's calculations use, so it
runs in Chrome and Edge from 91, of May 2021, Firefox from 89, of June
2021, Safari from 16.4, of March 2023, and node from 16.4, of June 2021.
On an iPhone or an iPad that means iOS 16.4, since every browser there is
WebKit whatever its name. An older browser fails when the module is
loaded, not with a wrong number. Goal 3 of `docs/objectives.md` has the
decision and the option that was not taken.

## node and a page, from one build

wasm-bindgen generates its JavaScript for one environment at a time, its
`--target`. This package is built once, with `--target web`, and the
difference between node and a page is in the two entry points around it,
which export the same functions:

- Under node, `dist/node.js`. The loader that `--target web` generates
  fetches the wasm file from the address of its own JavaScript, and the
  `fetch` of node does not open a `file:` address, so this entry point
  reads `wasm/popnei_bg.wasm` with `node:fs` and hands the bytes to that
  loader. The tests use it.
- In a page or through a bundler, `dist/web.js`, which lets the generated
  loader fetch the wasm file beside the JavaScript.

node chooses the first through the condition `node` of the field `exports`
of `package.json`, and everything else gets the second. Both are tested
under node, the second with a `fetch` that reads the file, because a
browser is not run here.

## What a bundler does with the wasm file

The address the loader fetches is
`new URL("popnei_bg.wasm", import.meta.url)`, written in
`wasm/popnei.js`. What a bundler makes of it was tried with two:

- vite 8.3.0 writes the wasm file among what it builds and rewrites the
  address to it. Nothing else has to be done.
- esbuild 0.28.2, `--bundle --format=esm --platform=browser`, leaves the
  address as it is and copies no file, so `await init()` fetches an
  address where the server has nothing and fails. A user of esbuild copies
  `js/popnei/wasm/popnei_bg.wasm` beside the bundle that esbuild writes,
  which is what `import.meta.url` is the address of.

webpack was not tried.

The two other targets of wasm-bindgen are built by running the same
command line again with another `--target` and another `--out-dir`, and
importing from there instead. Neither has been tried here:

- `--target nodejs` writes CommonJS, `exports.version` and
  `require('fs')`, which loads the wasm file when it is imported, with no
  `init` to await. A package whose `type` is `module`, as this one is,
  cannot import it without renaming it to `.cjs` or putting it under a
  directory with a `package.json` of its own that says `commonjs`.
- `--target bundler` writes an ES module that imports the wasm file as a
  module, `import * as wasm from "./popnei_bg.wasm"`, for the bundler to
  instantiate, again with no `init`. A bundler that cannot import a wasm
  file as a module needs a plugin for it.

## Using it

```ts
import { init, openVcf, version } from "popnei";
import { readFile } from "node:fs/promises";

await init();
console.log(version());

const variants = openVcf(new Uint8Array(await readFile("cases.vcf")), {
  ploidy: 2,
  onlyPassed: true,
});
console.log(variants.individuals, variants.numIndividuals, variants.ploidy);
try {
  // A filter is a step: it changes the handle, gives nothing back and is
  // run by every pass that follows. This one keeps the variants whose
  // commonest allele is at most 0.95 of their called alleles, 2 of the 3
  // of cases.vcf that passed their FILTER.
  variants.filterByMaf(0.95);
  const blocks = variants.iterBlocks({ fields: ["chrom", "pos"] });
  for (const block of blocks) {
    // block.gts is an Int8Array of numVars x numIndividuals x ploidy
    // alleles, variant after variant, with -1 for an allele that was not
    // called: the alleles of the individual i of the variant v start at
    // (v * block.numIndividuals + i) * block.ploidy.
    console.log(block.numVars, block.chrom, block.pos);
  }
  // How many variants the pass gave, and what each filter of the variants
  // was given and kept: {numVars: 2, filtering: {maf: {varsProcessed: 3,
  // varsKept: 2}}} here. Read inside the loop, it is of the blocks that
  // have come out so far. `steps` is what the handle holds, in order:
  // [{kind: "maf", args: {maxAllowedMaf: 0.95}}].
  console.log(blocks.passStats, variants.steps);
} finally {
  variants.free();
}
```

`init` has to be awaited before any other function of the package, which
throw an `Error` that says so until it has. It loads the WebAssembly once:
a second call gives the same promise as the first.

`openVcf` takes the bytes of the file, plain or gzipped, and reads its
header, so bytes that are not a VCF throw there and not at the first block.
It takes the `File` a user picked in a page where it takes those bytes, and
then it reads the header out of the file itself, which "A file of the page,
in a web worker" below has. Every call of `iterBlocks` reads the source
again from its start, so the same `Variants` can be given to one calculation
after another.

The variants of that handle are written into a vars file, and read back,
with the two functions of `docs/specs/io_vars.md`:

```ts
import { init, openVars, openVcf, writeVars } from "popnei";

await init();
const fromTheVcf = openVcf(new Uint8Array(await readFile("cases.vcf")));
// The bytes of the whole file, which a page offers as a download: a tab
// has no filesystem, and `passStats` says how many variants were written.
// The batches hold 10000 variants each, and without `numVarsPerBlock` the
// size popnei chooses for the individuals.
const { bytes, passStats } = writeVars(fromTheVcf, {
  numVarsPerBlock: 10000,
});
console.log(passStats.numVars);
fromTheVcf.free();

const fromTheFile = openVars(bytes);
for (const block of fromTheFile.iterBlocks({ fields: ["qual"] })) {
  console.log(block.numVars, block.qual);
}
fromTheFile.free();
```

`writeVars` reads the whole source once and writes every column a VCF has,
the chromosome, the position, the id, the alleles, the quality and the
genotypes, whether or not the user will read them, so that the file can
stand in for the VCF in any later analysis. `openVars` reads the schema of
the file and its footer, so bytes that are not a vars file throw there and
not at the first block; a file whose buffers are compressed with zstd
opens and throws at its first block, because arrow decompresses a batch
when it reads it and no build of popnei carries the code that reads zstd.
A `Variants` of a vars file is a source like the one of a VCF: it goes to
`iterBlocks` and back to `writeVars`, which writes the file again with
another size of batch.

The first calculation over such a handle is the Kosman distance of every
pair of individuals, `docs/specs/dists.md`:

```ts
import { calcPairwiseKosmanDists, init, openVcf } from "popnei";

await init();
const variants = openVcf(new Uint8Array(await readFile("panel.vcf.gz")));
try {
  // A pair called at fewer than 100 variants gets no distance, and is NaN
  // in the vector; without `minNumSnps` every pair called at one variant
  // at least gets one.
  const distances = calcPairwiseKosmanDists(variants, { minNumSnps: 100 });
  // The distance of every pair, in the order (0, 1), (0, 2), ..., (1, 2),
  // ...: 19900 values for the 200 individuals of that file.
  console.log(distances.distVector.length, distances.distVector[0]);
  // The names of the individuals, in the order the source has them, and
  // the counts of the pass the distances were calculated over.
  console.log(distances.names[0], distances.passStats.numVars);
  // The same distances as the 200 x 200 matrix, row by row, with 0 on the
  // diagonal: the distance of the individuals i and j is at i * 200 + j.
  console.log(distances.squareDists().length);
} finally {
  variants.free();
}
```

It reads the source once, through the filters that are on the `Variants`,
and leaves it as it was, so the same handle goes to the next calculation.
A pass that gives no variant is an `Error` that says whether the source
held none or the steps kept none, and for the steps how many variants each
filter was given and kept.

One pass gives the five statistics of every variant and every population:

```ts
import { calcPerVarDistribs, init, openVcf } from "popnei";

await init();
const panel = openVcf(new Uint8Array(await readFile("panel.vcf.gz")));
const distribs = calcPerVarDistribs(panel, {
  // The five when `stats` is left out. The populations are looked up among
  // the individuals the pass gives, which are the ones a
  // `filterIndividuals` kept when the variants carry one, and with no
  // `pops` there is one population, `pop`, of every individual.
  stats: ["maf", "poly_vars_ratio"],
  pops: { p0: ["s000", "s001"], p1: ["s002", "s003"] },
  // How many called genotypes a population needs at a variant to have a
  // value there, 20 when it is left out.
  minNumIndividuals: 2,
  histKwargs: { range: [0, 1], numBins: 40, binType: "linear" },
});
// ["p0", "p1"], the order the keys of `pops` iterate in, which is the
// order of every array below, and NaN for a population in which no variant
// had a value.
console.log(distribs.pops, distribs.maf?.mean);
// The 41 edges of the bins, and the counts of `p0` and then those of `p1`:
// the count of the bin b of the population p is at p * numBins + b.
console.log(distribs.maf?.histBinEdges, distribs.maf?.histCounts);
// The polymorphic variants of each population, those that vary at all, and
// the ones that have a major allele frequency there.
console.log(distribs.polyVarsRatio?.numPoly, distribs.passStats.numVars);
panel.free();
```

A statistic that was not asked for is `null`, and asking for fewer is a
saving of work that changes no value. A variant has no value of a statistic
in a population when the population has too little data at it, and such a
variant is out of the mean and in no bin, so the histograms of two
populations can count different numbers of variants.

Another pass gives the two rates of every individual:

```ts
import { calcPerIndividualStats } from "popnei";

const stats = calcPerIndividualStats(panel);
// The names of the individuals the pass gave, in its order, which is the
// order of the two arrays: the rate of `individuals[i]` is at `i` in each.
console.log(stats.individuals, stats.missingGtRate, stats.obsHetRate);
```

The heterozygosity rate divides by the called genotypes of the individual,
where pyNei's `calc_per_sample_stats` divides by every variant, so an
individual with more missing data looks less heterozygous there: `s000` of
the panel is heterozygous at 426 of its 1166 called genotypes, 0.365352,
and at 426 of the 1200 variants, 0.355. popnei's number is what plink2's
`--het` gives, and the missing rate beside it says what pyNei's one number
said. An individual that called no genotype has a missing rate of 1 and a
heterozygosity rate of NaN.

A file written here is larger than the same one written by popnei outside
the browser: `many.vcf` of `tests/reference/vcf/`, every variant of it in
batches of 100, is 53650 bytes written in wasm and 49426 bytes written
natively. The compression is lz4 in both, from `lz4_flex`, which hashes
four bytes of the input on a 32 bit target and five on a 64 bit one and so
finds other repetitions. Both files hold the same table and each library
reads both, and no test compares the two sizes.

The arguments are checked before they reach the core, and each of these is
an `Error` that says what was given: a `source` that is neither a
`Uint8Array` nor a `File` or a `Blob`, a `Uint8Array` whose buffer was
transferred, to a web worker or elsewhere, which leaves it with no bytes to
read, a `File` or a `Blob` given where there is no `FileReaderSync`, a
`ploidy` or a `numVarsPerBlock` that is not a whole number
of 1 or more and at most 4294967295, an `onlyPassed` that is not a
boolean, a `fields` that is not an array of names, a name that is not one
of the five columns, a `variants` that is not what `openVcf` or `openVars`
gave, a `minNumSnps` that is not a whole number of 0 or more and at most
4294967295, which a negative one is, a threshold of a filter that is not a
number, which a call with no threshold gives, and, of
`calcPerVarDistribs`, a `stats` that is not an
array of names or that names no statistic, a `pops` that is not an object
of population name to an array of names, a `minNumIndividuals`, a `ploidy`
or a `numBins` that is not a whole number of 0 or more, a `range` that is
not two numbers, a `binType` that is not a name, a `polyThreshold` that is
not a number, and a key of `histKwargs` that is none of the three, which
pyNei ignores. Whether that number is one a filter takes, from 0 to 1,
is a rule of the core, which holds for the threshold of every pass and not
of that call alone; an `Error` of it names the argument the user wrote and
the value as they wrote it, `95` and not `95.0`. A second filter of a kind
the variants carry already is an `Error` too, which names that kind and, for
a threshold filter, the threshold it is set with and the one that was
refused. The names given to `filterIndividuals` are read against the
individuals of the source at the call, so a name that is not one of them, a
name that is there twice and a call with no name are each an `Error` there
and not when a pass runs. Which names the core knows is the core's to
refuse: a statistic that is none of the five and a kind of bins that is
neither `linear` nor `logarithmic` are an `Error` of the binding crate,
whose message writes the names there are. In TypeScript
`fields` and `stats` take their names and nothing
else, so a typo does not compile.

## A file of the page, in a web worker

An application in a browser tab gets a `File` when its user picks a file in
a form or drops one on it: a handle that carries the name and the size of
the file and gives any range of its bytes, and that costs nothing to hold or
to send to a web worker, the thread of the page that cannot touch what the
page shows, because it is a handle and not the bytes. `openVcf` and
`openVars` take that `File` where they take a `Uint8Array`, and a `Blob`,
the piece of bytes of a page that a `File` is one of, is read the same way:

```ts
// In a module web worker. The page sends it the File its user picked and
// reads what comes back with worker.onmessage.
import { calcPerIndividualStats, init, numPassesOf, openVcf } from "popnei";

self.onmessage = async (picked: MessageEvent<File>) => {
  await init();
  const variants = openVcf(picked.data);
  try {
    // The bar of the page covers the whole run. How many passes that is,
    // one for this calculation, is the `numPasses` of every call below and
    // is asked here for the bar that is drawn before a byte is read.
    const passesOfTheRun = numPassesOf("calcPerIndividualStats");
    self.postMessage({ done: 0, numPasses: passesOfTheRun });
    variants.onProgress(({ bytesRead, numBytes, pass, numPasses }) => {
      self.postMessage({
        done: (pass - 1 + bytesRead / numBytes) / numPasses,
      });
    });
    self.postMessage(calcPerIndividualStats(variants).obsHetRate);
  } finally {
    variants.free();
  }
};
```

What that pass holds in the memory of wasm is one range of the file, of a
few MiB, the block it is building and, over a vars file, the batch it is
reading, which at the size popnei writes is about 10 MB of genotypes for
1000 individuals. The file itself is never there: popnei asks the `File` for
the ranges it needs, one at a time, through `FileReaderSync`, the reader
that returns when it has the bytes of a range. So the size of a file a user
opens stops being bounded by the memory of the tab, where a `Uint8Array`
costs a copy of the whole file inside wasm that stays there for as long as
the page lives, which "What crosses between Rust and JavaScript" above
measures. Which size of range
popnei reads by has not been measured yet, nor what a pass over a `File`
costs against a pass over the same file in memory: both are the measurement
that "Speed" of `docs/specs/js_sources.md` asks for.

A browser gives `FileReaderSync` only inside a web worker. `openVcf` and
`openVars` of a `File` or a `Blob` on the main thread of a page, and under
node, throw an `Error` at the call that names the reader and says to open
the file inside a worker or to give its bytes; a `Uint8Array` is read
wherever it is given. A range that comes back shorter than the one popnei
asked for, inside a file of that size, is an `Error` in the middle of the
pass and not the end of the file: a browser gives a short range when the
file changed on disk after the page got its handle, and a reader that took
it for the end would give the variants it had and say nothing.

Every pass reads the file again from its start, as a pass over bytes does,
so one `Variants` goes to one calculation after another. What it holds until
its `free()` is the handle of the file, and a pass that is still reading
when `free()` is called reads on to its end.

`onProgress` sets the function that is told how far every pass over that
source has got, with four numbers: how many bytes the pass has read, how
many the file holds, which pass of the run is reading and how many passes
the run makes. While a calculation runs the worker is inside wasm and reads
no message of the page, so this is how the page learns that the run is going
forward.

Three things make a call: the first read of a pass, which says it has read
nothing; a read that brings the bytes read since the last call to the few
MiB of a range; and the end of the run, which makes one call for each of its
passes, in the order of their numbers, with the bytes that pass read. The
last of the three is what says that a pass is over, because no read does: a
pass over a vars file stops after its last batch, and a run that fails stops
where it failed. That last call is not always `bytesRead === numBytes`, so an
application that waits for those two to meet waits for ever over a vars file,
which popnei does not read whole: it reads its last ten bytes, then its
footer, then each batch, and never the schema message at the head of the
file, since the footer carries the schema too.

What that function throws ends the pass where it was reading
and the calculation throws that same value back, which is how an application
cancels a run without ending its worker and recognises its own cancel with
`===`; the `Variants` is then the one it was, and the next run over it reads
the file from its start. `numPassesOf` answers that fourth number before a
run starts, so that a bar covers the run and not each pass from the moment
it is drawn: every consumer makes one pass, except two. `doPcaFromVariants`
makes two when its `numPrinComps` is above 0, which is what asks for the
weight of every variant: a weight needs the components, and those are known
when the first pass ends. `calcGwas` makes two when its
`useGrammarGammaApprox` is true, which stands the approximation in for the
denominator of the test of a mixed model, the factor of that approximation
being estimated from the first block of the second pass; with the false it
has by default, the study reads the file once.

Where an application sets the function that is told the progress is not
settled. `docs/specs/js_sources.md` leaves it open between the method of
`Variants` that is written here, an option of `openVcf` and `openVars` that
would hold for the life of the source, and an argument of every consumer;
the method is what this package gives until the owner decides.

## What has to be freed

The objects of the core live in the memory of the WebAssembly, which the
garbage collector of JavaScript does not see, so they are given back by
hand:

- The `Variants` of `openVcf` and of `openVars` holds its steps until its
  `free()` is called, which `using variants = openVcf(bytes)` does at the
  end of its block, and with it the bytes of the file when it was opened
  over a `Uint8Array`. Opened over a `File` it holds the handle of the file,
  its name and its size, and the memory of wasm keeps nothing of the file
  between two ranges. Its names and its ploidy are in JavaScript and answer
  after that; `iterBlocks`, `writeVars` and `steps` throw.
- One pass over the variants holds the reader, the block being built and,
  over a `File`, the range it is reading. What `iterBlocks` gives back gives
  it back when the iteration ends, when it is left with a `break` and when a
  block throws. An iterator that is made and never iterated keeps it until
  the garbage collector reaches it: wasm-bindgen registers what it generates
  in a `FinalizationRegistry`, which frees it at a moment nobody chooses.
  The `finally` that frees it cannot do that one, because a generator that
  never ran its first line never runs its last either. Its `passStats` is
  read after the pass is over all the same: the counts are taken out of the
  pass just before it is freed, and they are numbers of JavaScript.
- The `Uint8Array` of `writeVars` is the user's own, in the heap of
  JavaScript: the file is read out of the memory of wasm in pieces, each
  of them freed there as it is copied, so nothing of it is left to free by
  hand.
- What `calcPerVarDistribs` gives holds nothing of the memory of wasm: the
  means, the edges of the bins and the counts are copies, in the heap of
  JavaScript, and the object of the core they were read out of is freed
  before the call returns. What `calcPerIndividualStats` gives is the same:
  the names of the individuals and the two rates are copies.
- Each block is freed as soon as its columns are copied out, which is
  before it reaches the loop of the user. What the user holds are the
  copies: an `Int8Array` of genotypes, a `Float64Array` of positions and
  arrays of strings, none of them a view into the memory of the
  WebAssembly, which stops being valid when that memory grows.

## What it was built with

node 26.8.2, npm 11.19.1, TypeScript 5.9.3, the `wasm-bindgen` command
line 0.2.128 and rustc 1.98.0, on macOS on aarch64.
