# What a range of a file costs a pass in a browser, task 4.1

24 September 2026. Task 4.1 of `docs/plans/js-sources.md` measures what the
size of a range costs popnei in Chromium. A **range** is the block of bytes
popnei asks a browser for at a time when it reads a file the user picked in
the page, instead of taking the whole file into the memory of WebAssembly, as
the applications do today; `docs/specs/js_sources.md` says how that reading
works and leaves its size to this measurement. A **pass** is one reading of
the variants of a file from its start to its end.

The two numbers this measurement sets, which the spec now carries:

- The size of a range stays at 4 MiB. One pass over a VCF of 299994147 bytes
  took 996 ms at 4 MiB against 921 ms with the file in the memory of wasm
  whole, 8.1 % more, and held 14155776 bytes of that memory against
  302383104, 21.4 times fewer.
- What a pass holds in the memory of wasm at 4 MiB is 14155776 bytes, the
  same figure in each of the four rounds of runs. The bound the browser test
  of task 4.2 asserts is 25165824 bytes, 24 MiB, and the file of that test,
  299994147 bytes, is 11.9 times it.

Reading by ranges does not cost much more than reading the file whole, which
is what the plan named as the finding that would stop the work: 75 ms of a
pass of about a second, over a file of 300 MB. The whole file row of the
table below does not include what it costs to get those bytes out of the
`File` first, 19 ms, which an application pays before it can call `openVcf`
and which a pass by ranges pays inside its clock; with it the difference is
56 ms, 6.0 %.

## The machine, the file and the script

The owner's Apple M5 Pro, 18 cores, 64 GB, macOS 27.0, on 24 September 2026.
Chromium 1243 of Playwright 1.63.0, headless, which Playwright downloads.
The WebAssembly is the release build of `crates/popnei-js` for
`wasm32-unknown-unknown`, one build for each size of range.

The machine was not idle: two other worktrees of popnei were running their
own work throughout, a Python process at 220 % of a core and a `cargo test`
at 200 %, and the load average stood between 3.7 and 5.0 over the round the
tables below come from, which is about four of the eighteen cores busy. A
pass runs on the one thread of one web worker, so what that load costs it is
the contention for a core and for the memory. The five runs of one point
spread by as much as 73 % because of that, which is why each point is the
best of its runs. The load averages of every round are in the last section.

The file is built in the web worker, not read from a disc: the 617 byte
header of `tests/reference/vcf/many.vcf` with the 116729 bytes of that file's
body after it 2570 times, 299994147 bytes and 1285000 variants of 50 diploid
individuals, of which 64250 fail their FILTER and are read all the same,
which is what popnei does when it is not asked for the variants that passed
alone. It is the file of the browser test that reads a file of more than one
range, `js/popnei/test/browser/cases/many_ranges.ts`, with more copies: the
body is one block of memory that the `File` reads again for each copy, so
the worker holds 117346 bytes of it. What that leaves out is under "What
this does not measure".

The consumer is `calcPerIndividualStats`, which makes one pass, reads the
genotypes and adds two counts per individual, so the clock holds the reading
of the file and almost nothing else.

The script is `js/popnei/bench/time_ranges_in_chromium.mjs`, with the page
`bench/ranges_page.html` and the worker `bench/ranges_worker.js` beside it.
It is run from `js/popnei` after `npm run build`:

```
node bench/time_ranges_in_chromium.mjs
```

The size of a range is `NUM_BYTES_PER_RANGE` of
`crates/popnei-js/src/source.rs`, a constant of the build and no argument of
the API, so the script writes each size into that constant, runs `npm run
build:wasm`, takes the times, and writes the file back as it found it before
it ends, whatever happened. Each measurement gets a page and a worker of its
own: the memory of a WebAssembly module never shrinks, so what it holds when
a pass has ended is the most it held while the pass ran, and a second pass in
the same worker would read the high water mark of the first.

## One pass over 299994147 bytes

Each point is the best of 5 runs in one worker. The memory of wasm is
`memory.buffer.byteLength` after the fastest of those runs, read as
`js/popnei/test/vars_memory.test.ts` reads it.

| the size of a range | one pass | of the file | over the whole file | the memory of wasm |
|---|---|---|---|---|
| 262144 bytes, 256 KiB | 1069 ms | 280.7 MB/s | 16.1 % | 3407872 |
| 1048576 bytes, 1 MiB | 1006 ms | 298.3 MB/s | 9.2 % | 4718592 |
| 4194304 bytes, 4 MiB | 996 ms | 301.2 MB/s | 8.1 % | 14155776 |
| 16777216 bytes, 16 MiB | 958 ms | 313.2 MB/s | 4.0 % | 51904512 |
| the whole file as a `Uint8Array` | 921 ms | 325.8 MB/s | | 302383104 |

The memory of wasm before the first pass, in every row, is 1310720 bytes,
what the module holds once it is loaded.

## What one call into JavaScript costs

Each range is one `Blob.prototype.slice` and one
`FileReaderSync.readAsArrayBuffer`, a call from the memory of wasm out into
the browser and the bytes copied back. It is measured over the file of the
tables above, which the worker built out of 2571 pieces, and over the same
bytes as one block, which is what a file the user picked is: a file of many
pieces could otherwise hide the cost of walking the pieces a range falls in
inside the cost of the call. Each point is one set of calls, as many as fit
in 256 MB of bytes read.

| the bytes of one call | over the file of 2571 pieces | over the same bytes as one block | the calls |
|---|---|---|---|
| 1024 | 109.8 µs | 79.3 µs | 2000 |
| 262144 | 146.2 µs | 131.5 µs | 1024 |
| 1048576 | 302.7 µs | 278.9 µs | 256 |
| 4194304 | 776.6 µs | 690.6 µs | 64 |
| 16777216 | 2393.8 µs | 1731.2 µs | 16 |

A call costs about 79 µs before it has read a byte, and the pieces of the
file add at most 30 µs of that at the sizes of a range. That fixed cost is
the shape of the first table. A pass over 299994147 bytes makes 1145 calls at
256 KiB, 287 at 1 MiB, 72 at 4 MiB and 18 at 16 MiB, and the calls of a pass
cost, at the one block rates, 150.6, 80.0, 49.7 and 31.2 ms. The differences
between those four, 101, 30, 0 and -18 ms against 4 MiB, are of the size of
the differences the first table measured, 73, 10, 0 and -38 ms; the rest of
each row is the parsing of the VCF, which every row pays alike.

Reading the whole 299994147 bytes out of the `File` in one call took 19 ms,
the best of 3, which is 15.8 GB/s and is what an application pays today
before `openVcf` sees a byte.

## What a pass holds in the memory of wasm, and the bound of task 4.2

At 4 MiB the memory of the module goes 1310720 bytes when it is loaded,
5636096 after `openVcf` has read the header, 14155776 when the pass has
ended. At 16 MiB the same three are 1310720, 18219008 and 51904512. So a
pass holds about 3.4 times its range, and not the one range that "How it
runs" of the spec describes.

Where the rest goes: the function of `crates/popnei-js/src/source.rs` that
reads a range, `RangesOfAFile::reads_the_range`, copies the bytes the browser
gave it into a new block of the memory of wasm and puts that block in the
place of the range the pass held, so the two are live for that moment, and
the memory of wasm keeps the high water mark of everything that has been
live. Two ranges of 4 MiB are 8 MiB of the 12.25 MiB
above the module's own 1.25 MiB, and the batch of lines and the block of
genotypes the VCF reader holds are the rest. Copying into a buffer the pass
keeps, instead of building a new one for each range, would hold one range
where it now holds two; it is a change of the code that this task does not
make, and 14155776 bytes is 4.7 % of what the whole file costs.

The bound for the test of task 4.2 is 25165824 bytes, 24 MiB. It is 1.78
times what was measured, which leaves room for a consumer that holds more
than `calcPerIndividualStats` does, and the file of that test, 299994147
bytes, is 11.9 times the bound, where the plan asks for a file of at least
eight times it.

## Why 4 MiB

16 MiB is the fastest size measured, by 38 ms of a pass of a second, 3.8 %,
and it holds 51904512 bytes of the memory of wasm against 14155776, 3.7
times more. 1 MiB holds 4718592 bytes, a third of what 4 MiB holds, and cost
10 ms more than 4 MiB in this round and 12, 15 and 61 ms more in the three
rounds before it: the two are within 1 % of each other on one pass, and 4 MiB
was the faster of them in each of the four rounds.

What the memory buys is what the applications came for. A tab that reads a
file of 2 GB, which is above what a wasm module can hold whole, pays 14 MB
for it at 4 MiB and 52 MB at 16 MiB, and a source with several passes open at
once holds one range for each: the twelve passes of
`js/popnei/test/vars_memory.test.ts` are 48 MiB of ranges at 4 MiB and 192
MiB at 16 MiB. Against that, 38 ms of a second is what 4 MiB gives up, and
the number stays where the spec left it.

Whether an application chooses the size of a range instead is the third open
point of `docs/specs/js_sources.md` and is still the owner's to decide.
Nothing of this measurement asks for the option: since the
cost of a smaller range is time and no result, and 256 KiB, which tells a
progress bar four times as often as 1 MiB does, costs 148 ms of a pass over
this file against the whole file's 921 ms.

## What this does not measure

- A file on a disc. The `File` is built in the worker out of blocks of
  memory, so every call into JavaScript reads memory the browser already
  holds. A file the user picked is read from the disc through the same call,
  and what that adds was not measured; it falls on the calls of a pass, of
  which a smaller range makes more.
- Any browser but Chromium. Firefox and WebKit are not in the tests of
  popnei, which the owner left for when an application needs them, and
  `FileReaderSync` is theirs to implement as they see fit.
- Any consumer but `calcPerIndividualStats`, and any file but a plain VCF. A
  gzipped VCF and a vars file, popnei's own file of variants, read through
  the same ranges, and a vars file jumps to its footer and back inside a
  range, which no timing here covers.
- A quiet machine. Every round was taken while two other worktrees of popnei
  were building and testing.

## The four rounds, and what they spread over

Each round ran the four sizes and the whole file, 5 runs of each, in the
order of the table. The tables above are round 4. The best run of each point,
in milliseconds:

| | round 1 | round 2 | round 3 | round 4 |
|---|---|---|---|---|
| 256 KiB | 1101 | 1162 | 1053 | 1069 |
| 1 MiB | 1026 | 1053 | 1004 | 1006 |
| 4 MiB | 1010 | 992 | 992 | 996 |
| 16 MiB | 990 | 971 | 967 | 958 |
| the whole file | 955 | 916 | 925 | 921 |
| the load average before | 3.57 | 6.74 | 5.10 | 3.79 |
| the load average after | 9.12 | 6.97 | 5.17 | 4.22 |

The order of the five rows is the same in every round. The five runs inside a
round spread further than the rounds do: at 256 KiB round 4 gave 1092, 1069,
1109, 1236 and 1783 ms, and at 4 MiB it gave 1213, 1459, 1016, 1000 and 996.
The slow runs are the machine, which is why each point is the best of its
runs. The memory of wasm was the same byte count in every round.

Round 1 did not measure the whole file in one call nor the calls over one
block, which were added to the worker after it.

## The checks

The checks of the `coding` skill on this work, each with its last line, run
on this Mac on 24 September 2026:

| | |
|---|---|
| `cargo fmt --all --check` | nothing to say, 0 |
| `cargo clippy --workspace --all-targets -- -D warnings` | `Finished \`dev\` profile [unoptimized + debuginfo] target(s) in 0.40s` |
| `cargo test --workspace` | `787 passed; 0 failed; 2 ignored` and `149 passed; 0 failed` |
| `cargo test -p popnei --no-default-features` | `787 passed; 0 failed; 2 ignored` |
| `cargo wasm-check` | `Finished \`dev\` profile [unoptimized + debuginfo] target(s) in 0.09s` |
| `cargo wasm-check-js` | `Finished \`dev\` profile [unoptimized + debuginfo] target(s) in 0.50s` |
| `uv run ruff format --check && uv run ruff check` | `33 files already formatted`, `All checks passed!` |
| `uv run maturin develop && uv run pytest` | `499 passed in 11.18s` |
| `npm run build && npm test` in `js/popnei` | `tests 372`, `fail 0` |
| `npm run test:browser` in `js/popnei` | `7 passed (815ms)` |

The seven browser tests are the ones work packages 2 and 3 left; the eighth,
which asserts the bound above, is task 4.2. Nothing of popnei changed in this
task but a doc comment, so no test of it fails before these commits.

## The commits

- 91cf7ba: the "Speed" part of `docs/specs/js_sources.md`, which carries the
  two numbers, and the sentence of "How it runs" that said a pass holds one
  range.
- 9233fc5: the doc comment of `NUM_BYTES_PER_RANGE` in
  `crates/popnei-js/src/source.rs`, whose value does not change, and the
  script, the page and the worker that took the times.
