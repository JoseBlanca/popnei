# Work report: a file of a page read by ranges in the wasm package

24 September 2026. The work of `docs/plans/js-sources.md`, on the branch
`plan/js-sources`, which builds the source of bytes over a `File` of
`docs/specs/js_sources.md`. It is under way; this top section says what the
owner needs when it is done.

## Before the first task

The branch starts from `spec/file-source` at 465c5c0 and not from `main`,
because the spec is not in `main` yet. What was checked on this Mac on 24
September 2026, by running it:

- `npm run build` in `js/popnei`: the wasm and the TypeScript built, 33 s.
- `npm test` in `js/popnei`: `tests 325`, `fail 0`, 4.8 s.
- `npm view @playwright/test version`: 1.63.0; Chromium 1194 is already in
  `~/Library/Caches/ms-playwright`.
- `node --version`: v26.8.2, which has `Blob` and `File` and no
  `FileReaderSync`.
- The checks of the `coding` skill on the starting commit: `cargo fmt
  --all --check` and `cargo clippy --workspace --all-targets -- -D
  warnings` clean, `cargo test --workspace` 787 and 149 passed, `cargo
  test -p popnei --no-default-features` 787 passed, `cargo wasm-check`
  clean, `ruff` clean, `uv run pytest` 499 passed in 21.4 s.

## Work package 1

Task 1.1, `numPassesOf`, is done at 1b2c50f. `node --test
test/num_passes.test.ts` runs 14 tests, where deliverable 1 asks for 6 or
more, and `npm test` gives `tests 339`, `fail 0`, 14 more than the 325 of
the starting commit. The `Consumer` enum and its `num_passes` are in
`crates/popnei-js/src/source.rs`, beside `OpenSource`, which is where task
1.2 takes them from.

What the task found: the spec says `numPassesOf` throws before `init` has
been awaited and does not ask for a test of it under "How it is verified",
and `js/popnei/test/before_init.test.ts` is where the package tests that
for its other functions. The orchestrator gave that test to task 1.2,
which writes test files of its own.

## Work package 2: the browser harness

Done at aea65cf, in one task. `npm run test:browser` in `js/popnei` runs
Playwright 1.63.0 with Chromium headless, which starts a server over the
repository, a page and a module web worker that loads the wasm, and the
test that opens `many.vcf` from an array of bytes inside that worker and
asserts its first block. It took 3.5 s on this Mac, 1 test. `npm test` is
unchanged at `fail 0`, and `node --test test/` does not pick up the files
of `test/browser/`, which the task checked.

Two things the plan and the spec had wrong, both corrected at 2746ce5:

- The spec sent the browser test to the literals of
  `js/popnei/test/vcf.test.ts`, which asserts of the first block of
  `many.vcf` only that it holds 100 variants. The positions and the
  genotypes of that file are in `test/filter_individuals.test.ts`, which
  reads them from `docs/specs/io_vcf.md`.
- The plan said the harness would download little, because Chromium 1194
  was in the cache of this Mac. Playwright 1.63.0 asks for Chromium 1243:
  94.3 MiB were downloaded.

The review of this work package runs with work package 3, whose browser
tests are the same piece of code.

What the owner should know: a test page of TypeScript is served as
JavaScript with node's `stripTypeScriptTypes`, which is experimental in
node v26.8.2 and prints one warning per run. What it buys is that the
page, the worker and the cases are type checked with the rest of
`js/popnei/test/`.

Task 1.2, the run and the counting source, is done at 7c2bcf8. `node --test
test/progress.test.ts` runs 10 tests where deliverable 2 asks for 6, and
`npm test` gives 350 tests, `fail 0`. The size of a range is one constant,
4 MiB, whose doc comment says work package 4 measures it. `js-sys` came in
here and not in task 3.1: nothing else can call a function of JavaScript.

Two sentences of the spec were wrong and are corrected at the commit after
this one, both because a reader of popnei reads when it is built:

- The calls of the two passes of the PCA come as pass 1, pass 2, pass 1,
  pass 2, and not as the calls of the first pass and then those of the
  second. It builds both readers before it asks either for a block, and a
  VCF reader reads its header there, a vars file reader its footer. The
  share of a run that is done still rises, 0, 0.5, 0.5, 1, when a page
  draws `(pass - 1 + bytesRead / numBytes) / numPasses`, which the spec now
  gives. What would make the calls come in the order the spec first
  promised is a second reader built after the first pass ends, which is a
  change to the signature of `pca_of_variants` in the core.
- A vars file smaller than one range gives one call, of 0 bytes: a pass
  over it never reads past its last batch, so nothing finds the end of the
  file. The check of the spec now writes a vars file of more than one range,
  as `vars_memory.test.ts` does, and says what a small one gives.

A `told` that is neither a function nor nothing is refused where the
application sets it, which the spec did not say and now does.

Task 1.3, the stop, is done at bc46752. `node --test test/stop.test.ts` runs
7 tests where deliverable 3 asks for 4, five of which fail on the commit
before it. `JsPopneiError` gained a case that carries the value the
application threw through untouched, the read that a throw ends fails with
`std::io::Error::other`, and nine of the ten consumers open their run in one
function, so the swap of the error cannot be forgotten; `iterBlocks` is the
tenth, whose run outlives the call.

### What the review of work package 1 found

Seven reviewers, one per category. Two findings are defects a user would
meet, both measured:

- Freeing the `Variants` from inside the progress function leaks the file
  and leaves the handle unusable. 100 sources of `many.vcf` freed that way
  grew the memory of wasm by 11993088 bytes, one file each, where 100 freed
  the ordinary way grew it by 131072. wasm-bindgen zeroes the pointer and
  unregisters the finalization before it calls into wasm, and that call
  throws because a consumer holds the source.
- A pass is never told that it has ended, so a bar over a vars file stops
  short. A pass over a vars file of 12231602 bytes was last told at
  8476400 bytes, two thirds of the way: its reads after that never reached
  a range, and no read of a vars file finds the end of the file.

The spec was wrong in four more places, each corrected in a commit of its
own before the code: the calls of the two passes of the PCA interleave
(8c76e4f), `free()` from inside the function is refused while a consumer
runs (9a14d55), the end of a run is what says a pass is over (84cf028), and
a pass whose reader is built is always told of (df5e3e7).

The reviewer of the tests broke the code on purpose eleven times and named
the tests each mutation failed, which is the check that the suite can fail.
It found two tests that could not: the stop in the middle of a VCF would
pass with a stop at its end, and nothing pinned that a stopped pass is not
told again.

### Work package 1 is done

The fixes are at 2dae2d9, 03fcc91, 9c67db5, dd40285 and db4ff2c, and the
spec follows them at 7606b76. The deliverables, each run by the
orchestrator after the fixes:

| deliverable | command | what it gave |
|---|---|---|
| 1 | `node --test test/num_passes.test.ts` | 15 tests, 0 failed, where 6 were asked |
| 2 | `node --test test/progress.test.ts` | 11 tests, 0 failed, where 6 were asked |
| 3 | `node --test test/stop.test.ts` | 12 tests, 0 failed, where 4 were asked |
| 4 | the checks of the `coding` skill | fmt and clippy clean, 787 and 149 cargo tests, 787 with no default features, `wasm-check` clean, ruff clean, 499 pytest |
| 4 | `npm test` in `js/popnei` | 364 tests, 0 failed, the 325 of the starting commit among them |

Freeing a `Variants` from inside the progress function now costs 0 bytes of
the memory of wasm for 100 sources, where it cost 11993088 before.

What the owner should know for what comes next: `js-sys` is in the binding
crate and `cargo wasm-check` does not build that crate for wasm, so nothing
of the five commands would catch a call into JavaScript that does not
compile for the target the package ships. Work package 3 is where that is
put right, since it is the work that adds `web-sys`.
