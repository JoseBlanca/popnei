# popnei in TypeScript

The TypeScript package of popnei, a population genetics library whose
calculations are written in Rust: the core crate compiled to WebAssembly,
the code that a browser or node calls it through, and the functions and
the result objects an application uses. What the package exports today is
`init`, which loads the WebAssembly, `version`, the version of the core
crate, `openVcf`, which reads the header of a VCF held as bytes and
gives a `Variants`, the handle whose `iterBlocks` gives the genotypes block
by block, `writeVars`, which gives back the bytes of a vars file with every
variant of a `Variants`, and `openVars`, which opens such bytes as another
`Variants`. A vars file is one arrow IPC file, also called feather v2,
which pandas, R and polars open as a table with no popnei installed: it is
where a user keeps their variants once the VCF has been read. Section 11
of `docs/architecture.md` has the design, `crates/popnei-js` is the binding
crate, the Rust that is compiled to WebAssembly and that holds no
calculation of its own, and `docs/specs/io_vcf.md`,
`docs/specs/io_vars.md`, `docs/specs/block.md` and `docs/specs/variant.md`
say what they give.

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
takes out of the wasm file the section that holds the name of every
function of it: `js/popnei/wasm/popnei_bg.wasm` is 1242562 bytes with the
flag and 1687938 bytes without, 443209 bytes of names that every user of
the package downloads. What they are for is the stack of a trap, a panic
of Rust among the causes, which with the flag names the functions by their
number and without it by their name. To read one, build again without the
flag and make the trap happen there.

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
  API, the ploidy of 2 and the filter of `docs/specs/io_vcf.md`, cross as
  two functions that return the constants of the core.
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
  not 160 MB. Measured on a vars file of 26.9 MB, 20000 variants of 1000
  diploid individuals written in batches of 1000: `openVars` grows the
  memory of wasm by 27.1 MB, the file and nothing else, and a pass over
  its genotypes by 12.2 MB, the batches it decodes. A reader built over a
  copy of those bytes instead grew it by 39.1 MB at that pass, 26.9 MB
  more, which is the file.
- A `Vec<u8>` coming back is a copy going out, so the bytes of a file that
  is written are in the memory of wasm and in the `Uint8Array` at once.
  While `write_vars` runs, that memory holds the source, the block being
  written and the file that is growing: writing the 26.9 MB file above
  from a VCF of 76.9 MB grew the memory of wasm by 106.8 MB with batches
  of 1000 variants and by 158.5 MB with batches of 10000, beyond the
  76.9 MB of the source. The size of the batches is what an application
  that runs out of memory has to lower.
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
popnei cannot write and refuses at its first block. Three of the tests
watch the memory of the WebAssembly, which they reach through the loader
`wasm/popnei.js` generates: that a block kept while enough more are read
for that memory to grow still holds what it held, that an iteration gives
its pass back however it ends, and that twelve passes over one vars file
open at once grow that memory by less than the file, which a reader that
copied the bytes for each pass would not.

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
  for (const block of variants.iterBlocks({ fields: ["chrom", "pos"] })) {
    // block.gts is an Int8Array of numVars x numIndividuals x ploidy
    // alleles, variant after variant, with -1 for an allele that was not
    // called: the alleles of the individual i of the variant v start at
    // (v * block.numIndividuals + i) * block.ploidy.
    console.log(block.numVars, block.chrom, block.pos);
  }
} finally {
  variants.free();
}
```

`init` has to be awaited before any other function of the package, which
throw an `Error` that says so until it has. It loads the WebAssembly once:
a second call gives the same promise as the first.

`openVcf` takes the bytes of the file, plain or gzipped, and reads its
header, so bytes that are not a VCF throw there and not at the first block.
A `File` that a user picked in a page is read inside a web worker, which
section 11 of `docs/architecture.md` has and this package does not do yet.
Every call of `iterBlocks` reads the bytes again from their start, so the
same `Variants` can be given to one calculation after another.

The variants of that handle are written into a vars file, and read back,
with the two functions of `docs/specs/io_vars.md`:

```ts
import { init, openVars, openVcf, writeVars } from "popnei";

await init();
const fromTheVcf = openVcf(new Uint8Array(await readFile("cases.vcf")));
// The bytes of the whole file, which a page offers as a download: a tab
// has no filesystem. The batches hold 10000 variants each, and without
// `numVarsPerBlock` the size popnei chooses for the individuals.
const bytes = writeVars(fromTheVcf, { numVarsPerBlock: 10000 });
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

A file written here is larger than the same one written by popnei outside
the browser: `many.vcf` of `tests/reference/vcf/`, every variant of it in
batches of 100, is 53650 bytes written in wasm and 49426 bytes written
natively. The compression is lz4 in both, from `lz4_flex`, which hashes
four bytes of the input on a 32 bit target and five on a 64 bit one and so
finds other repetitions. Both files hold the same table and each library
reads both, and no test compares the two sizes.

The arguments are checked before they reach the core, and each of these is
an `Error` that says what was given: a `source` that is not a
`Uint8Array`, a `ploidy` or a `numVarsPerBlock` that is not a whole number
of 1 or more and at most 4294967295, an `onlyPassed` that is not a
boolean, a `fields` that is not an array of names, a name that is not one
of the five columns, and a `variants` that is not what `openVcf` or
`openVars` gave. In TypeScript `fields` takes the five names and nothing
else, so a typo does not compile.

## What has to be freed

The objects of the core live in the memory of the WebAssembly, which the
garbage collector of JavaScript does not see, so they are given back by
hand:

- The `Variants` of `openVcf` and of `openVars` holds the bytes of the
  file until its `free()` is called, which `using variants =
  openVcf(bytes)` does at the end of its block. Its names and its ploidy
  are in JavaScript and answer after that; `iterBlocks` and `writeVars`
  throw.
- One pass over the variants holds the reader and the block being built.
  The iterator of `iterBlocks` gives it back when the iteration ends, when
  it is left with a `break` and when a block throws. An iterator that is
  made and never iterated keeps it until the garbage collector reaches it:
  wasm-bindgen registers what it generates in a `FinalizationRegistry`,
  which frees it at a moment nobody chooses. The `finally` that frees it
  cannot do that one, because a generator that never ran its first line
  never runs its last either.
- Each block is freed as soon as its columns are copied out, which is
  before it reaches the loop of the user. What the user holds are the
  copies: an `Int8Array` of genotypes, a `Float64Array` of positions and
  arrays of strings, none of them a view into the memory of the
  WebAssembly, which stops being valid when that memory grows.

## What it was built with

node 26.8.2, npm 11.19.1, TypeScript 5.9.3, the `wasm-bindgen` command
line 0.2.128 and rustc 1.98.0, on macOS on aarch64.
