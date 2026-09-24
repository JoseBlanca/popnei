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
