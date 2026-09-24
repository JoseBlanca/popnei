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
