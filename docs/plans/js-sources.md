# Plan: a file of a page read by ranges in the wasm package

24 September 2026. The breakdown was approved by the owner in chat on 24
September 2026 and this document is that breakdown with its deliverables and
tasks; the state is under way, on the branch `plan/js-sources` in the
worktree `.claude/worktrees/js-sources`, with its work report in
`docs/reports/js-sources.md`. It builds from one spec,
`docs/specs/js_sources.md`, committed at 465c5c0 on the branch
`spec/file-source`, which this branch starts from because that spec is not
in `main` yet.

When it is done, a web application gives `openVcf` or `openVars` the `File`
the user picked in the page and popnei reads it a few MiB at a time, so the
file is never in the memory of wasm whole; the application is told, once per
range, how many bytes the pass has read, which pass of the run it is and how
many passes the run makes; a function of the application that throws ends
the run and gets its own value back, so a cancelled analysis keeps its
worker; and `numPassesOf` says beforehand how many passes a consumer makes.

## In and out

In: the progress, the stop and the runs over both kinds of source; the
source over a `File` with `js-sys` and `web-sys`; `numPassesOf`; the tests
under node and in Chromium; the size of a range, measured, with its report.

Out, with where it goes: Firefox and WebKit in the tests, which the owner
left for when an application needs them; the fingerprint of a variant file,
the reader of CSV and TSV and the inference of the types of its columns,
which are the other requests of the applications and each is its own work;
stopping a run while it is not reading, which would be a request to the core
crate.

No work package ends at a Python function and a comparison with pyNei, which
the plans of popnei otherwise do. This feature has no Python side: Python
opens a file by its path and has no `File`, as the spec says under "Its
TypeScript function". The numbers the browser tests assert are the literals
of the node tests, which are the ones already compared with pyNei.

The three open points of the spec are unanswered and each has its
"meanwhile", which the tasks follow. The task each answer would change: Open
1, where the application sets the function, task 1.2; Open 2, what a stopped
run throws, task 1.3; Open 3, whether an application chooses the size of a
range, task 4.1.

## What has to be in place

- The spec on the branch: `git log --oneline -1 -- docs/specs/js_sources.md`
  gives 465c5c0.
- The branch at that commit, where the checks of the `coding` skill pass and
  `npm test` in `js/popnei` gives `tests 325`, `fail 0`. Run on 24 September
  2026: the build of the package took 33 s and the 325 node tests 4.8 s.
- On the machine: `node` v26.8.2, `cargo`, the target
  `wasm32-unknown-unknown`, `wasm-bindgen`, and Playwright, which
  `npm view @playwright/test version` gives as 1.63.0. Checked on 24
  September 2026. The Chromium 1194 that was already in
  `~/Library/Caches/ms-playwright` is not the one Playwright 1.63.0 asks
  for: task 2.1 downloaded 94.3 MiB of Chromium 1243 and its headless
  shell.

The checks that say the work is done fail today: `js/popnei/test/` has no
`num_passes.test.ts`, `progress.test.ts`, `stop.test.ts` or `browser/`, and
`npm run test:browser` is not a script of `js/popnei/package.json`.

Tasks of this plan touch `crates/popnei-js/src/`, `js/popnei/src/`,
`js/popnei/test/` and, in work package 3, the manifest of the workspace and
`Cargo.lock`. Two tasks that write the same file never run side by side.

## Work package 1: the progress, the stop and the runs, over bytes in memory

What it gives: an application that opened a `Uint8Array` gets
`variants.onProgress(fn)`, whose function is called with the bytes the pass
has read, the size of the file, the pass and how many passes the run makes;
a function that throws ends the run and the consumer throws that same value
back; and `numPassesOf` says how many passes a consumer will make. Every
check of this work package runs under node.

Deliverables:

1. `node --test test/num_passes.test.ts` in `js/popnei` passes, with a test
   for each check of "How it is verified" of "How many passes a consumer
   makes": the PCA with `numPrinComps` 10, with 0 and with no options, each
   of the other nine names, a name that is of no consumer, and a
   `numPrinComps` of -1. 6 tests at least.
2. `node --test test/progress.test.ts` passes, with a test for each check of
   "How it is verified" of "What the source tells the page" that does not
   throw: the calls over `many.vcf`, the calls over `many.vcf.gz`, the calls
   over a vars file, the two passes of the PCA and its one pass with
   `numPrinComps` 0, and the twelve `iterBlocks` at once. 6 tests at least.
3. `node --test test/stop.test.ts` passes, with a test for each check of
   that same part that throws: a function that throws on its first call, one
   that throws on the first call of the second pass of the PCA, one that
   calls `variants.free()` and one that runs `calcPerIndividualStats` over
   the same `Variants`. 4 tests at least.
4. The checks of "Before the work is called done" of the `coding` skill
   pass, and `npm test` in `js/popnei` gives `fail 0` with the 325 tests
   that were there among them.

Stands on: nothing of this plan.

What could go wrong: `Consumer::num_passes` is the one place that knows how
many passes each consumer makes, and a consumer that opens a reader without
a run would be told nothing and counted nowhere. Deliverable 2's last test,
the twelve `iterBlocks`, is what catches a run that is counted in the
source instead.

Tasks:

- [x] 1.1 `numPassesOf`: the `Consumer` enum and its `num_passes` in
      `crates/popnei-js/src/`, the function wasm-bindgen exports, and
      `numPassesOf` in `js/popnei/src/`, exported from both entry points,
      which checks its name and its options as the spec says. From "How many
      passes a consumer makes". The tests of `test/num_passes.test.ts`
      written first. Serves deliverables 1 and 4.
- [x] 1.2 The run and the counting source: `PassOverTheBytes` with
      `TheBytes::InMemory`, the two tables `IN_JAVASCRIPT` and `RUNS`, `Run`
      and `RunOfAConsumer`, `OpenSource::starts_a_run` and its `reader`
      taking the run, every consumer of `crates/popnei-js/src/` opening its
      run, and `Variants.onProgress` in the package. From "What the source
      tells the page", its "What it gives", "Its TypeScript function" and
      "The cases a reader of the rules would not guess", and from "The Rust
      interface". The tests of `test/progress.test.ts` written first. Needs
      1.1 for `Consumer`. Serves deliverables 2 and 4. It added `js-sys`
      to the binding crate, which the plan had given to task 3.1: nothing
      else calls a function of JavaScript.
- [x] 1.3 The stop: the read that fails with `std::io::Error::other`, the
      case of `JsPopneiError` that carries a `JsValue` untouched, the
      `stopped_with` of the run and the swap of the error at the boundary,
      and the tables left unborrowed while the function runs. From the same
      parts of the spec and from the two paragraphs after the code of "The
      Rust interface". The tests of `test/stop.test.ts` written first. Needs
      1.2. Serves deliverables 3 and 4.

## Work package 2: the browser harness

What it gives: `npm run test:browser` in `js/popnei` runs Chromium headless,
loads popnei's wasm inside a web worker and asserts from in there the same
numbers the node tests assert. Nothing of popnei has run in a browser
before, so this is built before the work that needs it.

Deliverables:

1. `npm run test:browser` exits 0 and its output names a test that, inside a
   web worker, opens `many.vcf` from a `Uint8Array` and asserts the first
   block that `js/popnei/test/vcf.test.ts` asserts: the chromosomes, the
   positions and the genotypes, the same literals.
2. `npm test` in `js/popnei` still gives `fail 0`, and `npm run build` is
   unchanged in what it builds.

Stands on: nothing of this plan. It runs side by side with work package 1:
it writes `js/popnei/package.json`, `js/popnei/test/browser/` and the
configuration of Playwright, and no file of work package 1.

What could go wrong: the wasm is loaded in a worker from a page that
Playwright serves, which needs the files of `js/popnei` reachable over
http and the worker started as a module. If a worker cannot load it that
way, the harness is where that is found, and the browser tests of work
packages 3 and 4 are what would have to change.

Tasks:

- [x] 2.1 Playwright with Chromium as a development dependency of
      `js/popnei`, its configuration, a page and a module worker that load
      the built wasm, the npm script `test:browser`, and the first test,
      which is deliverable 1. From "How it is verified" of "The source over
      a file of the page", the part that names Chromium and Playwright.
      Serves deliverables 1 and 2.

## Work package 3: the source over a `File`

What it gives: `openVcf` and `openVars` take a `File` or a `Blob` and read
it by ranges, so an application opens a file larger than the memory of the
tab; a `Blob` outside a web worker is refused with popnei's message; and a
range that comes back short is an error and not the end of the file.

Deliverables:

1. `npm run test:browser` passes with three more tests: a `File` of
   `many.vcf` gives the first block that the node test asserts, a `File` of
   `many.vcf.gz` gives the same variants, and a vars file written with
   `writeVars`, made into a `File` and opened with `openVars` gives those
   same variants, which is the check of `Seek`.
2. A browser test in which `Blob.prototype.slice` is patched to give one
   byte less than it was asked for: the pass fails with popnei's error, and
   the message names the range.
3. `node --test test/blob_outside_a_worker.test.ts` passes: `openVcf` and
   `openVars` of a `Blob` under node throw an `Error` whose message names
   `FileReaderSync` and the web worker.
4. `js/popnei/README.md` says what a `File` costs and what stays in the
   memory of wasm, and `grep -c FileReaderSync js/popnei/README.md` is 1 or
   more.
5. The checks of the `coding` skill pass with `js-sys` and `web-sys` in the
   workspace, `npm test` gives `fail 0`, and `npm run build` builds the
   package.

Stands on: work packages 1 and 2.

What could go wrong: `FileReaderSync` from `web-sys` has never been called
from popnei's wasm. The spike of the spec compiled it and ran nothing. If
the call does not work from a module worker, this is where it shows, and
what would change is how the binding reads a range and not what the package
gives.

Tasks:

- [x] 3.1 The ranges: `js-sys` and `web-sys` in the manifest of the binding
      crate, `TheBytes::OfAFile` over `FileReaderSync` and `Blob.slice`, the
      short range as an error, the refusal outside a worker, and `openVcf`
      and `openVars` taking a `Blob` in the binding crate and in the
      package, with `BytesOrFile` and the doc comments the spec gives. From
      "The source over a file of the page", all four of its parts, and from
      "The Rust interface". Needs 1.2 and 1.3. Serves deliverables 1, 2, 3
      and 5.
- [x] 3.2 The browser tests of deliverables 1 and 2, in
      `js/popnei/test/browser/`. Needs 2.1 and 3.1. Serves deliverables 1
      and 2.
- [x] 3.3 The node test of deliverable 3, and the README and the doc
      comments of deliverable 4. Runs side by side with 3.2: it writes
      `js/popnei/test/blob_outside_a_worker.test.ts` and
      `js/popnei/README.md` and no file that 3.2 writes. Needs 3.1. Serves
      deliverables 3 and 4.

## Work package 4: the size of a range, measured

What it gives: popnei reads a file in ranges of a size that was measured
instead of chosen, and a test says that a file of a few hundred MB is read
without the memory of wasm passing a bound.

Deliverables:

1. A report under `docs/reports/` with the time of one pass over a VCF of a
   few hundred MB in Chromium at 256 KiB, 1 MiB, 4 MiB and 16 MiB, and over
   the same file read whole as a `Uint8Array`, each with the machine, the
   date and the file; the script that makes the file and takes the times is
   beside it or in the repository, so that the numbers can be got again.
2. The size of a range is a named constant of `crates/popnei-js/src/` whose
   doc comment gives that measurement, and the "Speed" part of
   `docs/specs/js_sources.md` carries the two numbers it said it would: the
   size of a range and what a pass holds in the memory of wasm.
3. `npm run test:browser` passes with one more test: a `File` of a few
   hundred MB, the body of `many.vcf` repeated, passed over once with
   `calcPerIndividualStats`, which asserts the number of variants and that
   `memory.buffer.byteLength` never passed the bound of deliverable 2, the
   file being at least eight times it.
4. The checks of the `coding` skill pass.

Stands on: work package 3.

What could go wrong: the measurement may show that reading by ranges costs
much more than reading the file whole, which is what the applications give
up the whole file in memory for. If it does, the number goes to the owner
with the report and the plan stops there, because what to do about it is a
decision and not a task.

Tasks:

- [x] 4.1 The measurement and what it sets: the page that times a pass at
      the four sizes and over the whole file, the report, the constant with
      its doc comment, and the "Speed" part of the spec, which is a commit
      of its own before the commit of the code. From "Speed". Needs 3.1 and
      3.2. Serves deliverables 1, 2 and 4.
- [x] 4.2 The test of deliverable 3. Needs 4.1 for the bound. Serves
      deliverables 3 and 4.

## How the whole plan is checked

The sum of the work packages, and one thing besides: `npm run build && npm
test && npm run test:browser` in `js/popnei`, from a tree with no `wasm/`
and no `node_modules/`, which is what an application that installs the
package gets.
