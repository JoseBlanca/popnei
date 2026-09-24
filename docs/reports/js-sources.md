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
