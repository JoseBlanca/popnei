# popnei in TypeScript

The TypeScript package of popnei, a population genetics library whose
calculations are written in Rust: the core crate compiled to WebAssembly,
the code that a browser or node calls it through, and the functions and
the result objects an application uses. The functions that open a dataset
and give its variants are being written; what the package exports today is
`init`, which loads the WebAssembly, and `version`, the version of the
core crate. Section 11 of `docs/architecture.md` has the design, and
`crates/popnei-js` is the binding crate, the Rust that is compiled to
WebAssembly and that holds no calculation of its own.

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

The version of the `wasm-bindgen` crate, in the `Cargo.toml` of the
workspace, has to be the version of the `wasm-bindgen` command line that
is installed, so it is pinned there, `=0.2.128`. The crate writes into the
wasm file what the command line then reads, and when the two versions
differ the command line stops and names both. A new command line,
`cargo install wasm-bindgen-cli`, needs that line of the workspace
manifest changed to its version.

## The tests

    npm test

It runs the TypeScript compiler over `test/` and then the test runner of
node itself, `node --test`, once the build has left `dist/` and `wasm/` in
place. The tests import the name of the package, `popnei`, which node
resolves to the built entry point of node, and assert that the version the
package gives and the version in `package.json` are both the one of
`[workspace.package]` of the `Cargo.toml` of the repository, that `init`
loads the WebAssembly once, that a function called before `init` was
awaited throws an `Error` that says so, and that the entry point of a page
answers with the WebAssembly it fetches.

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
import { init, version } from "popnei";

await init();
console.log(version());
```

`init` has to be awaited before any other function of the package, which
throw an `Error` that says so until it has. It loads the WebAssembly once:
a second call gives the same promise as the first.

## What it was built with

node 26.8.2, npm 11.19.1, TypeScript 5.9.3, the `wasm-bindgen` command
line 0.2.128 and rustc 1.98.0, on macOS on aarch64.
