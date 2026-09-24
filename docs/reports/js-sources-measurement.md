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
- What a pass leaves in the memory of wasm at 4 MiB is 14155776 bytes, the
  same figure in each of the four rounds of runs, and at every size from
  1 MiB up it is three ranges and 1572864 bytes, to the byte. The bound the
  browser test of task 4.2 asserts is 16777216 bytes, 16 MiB, 1.19 times the
  14155776, and the file of that test, 299994147 bytes, is 17.9 times the
  bound.

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
build:wasm`, takes the times, and then writes the file back as it found it
and builds the WebAssembly of `wasm/` from it again, whatever happened and
on an interrupt as well. Both are needed: the source alone put back leaves
`wasm/` built at the last size measured and `git status` clean, and the next
`npm run test:browser` fails at the case of a file of more than one range
and `npm test` at the counts of a pass, neither of them naming the script.

Each measurement gets a page and a worker of its own: the memory of a
WebAssembly module never shrinks, so what it holds when a pass has ended is
the most it held while the pass ran, and a second pass in the same worker
would read the high water mark of the first.

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

A call costs about 79 µs before it has read a byte. Walking the pieces of
the file a range falls in adds to that 14.7 µs at 256 KiB, 23.8 at 1 MiB,
86.0 at 4 MiB and 662.6 at 16 MiB, which is why the rates used below are
the ones over the same bytes as one block. That fixed cost of a call is the
shape of the first table. A pass over 299994147 bytes makes 1145 calls at
256 KiB, 287 at 1 MiB, 72 at 4 MiB and 18 at 16 MiB, and the calls of a pass
cost, at the one block rates, 150.6, 80.0, 49.7 and 31.2 ms. The differences
between those four, 101, 30, 0 and -18 ms against 4 MiB, are of the size of
the differences the first table measured, 73, 10, 0 and -38 ms; the rest of
each row is the parsing of the VCF, which every row pays alike.

Reading the whole 299994147 bytes out of the `File` in one call took 19 ms,
the best of 3, which is 15.8 GB/s and is what an application pays today
before `openVcf` sees a byte.

## What a pass leaves in the memory of wasm, and the bound of task 4.2

The memory of wasm is the whole module, what `memory.buffer.byteLength`
gives, and it never shrinks: it is the 1310720 bytes the module holds before
any pass and every block allocated since, whether that block is still live
or was freed. At a range of 4 MiB it goes 1310720 bytes when the module is
loaded, 5636096 when `openVcf` has read the header and 14155776 when the
pass has ended. At 16 MiB the same three are 1310720, 18219008 and 51904512.

The last figure of each of those two, 14155776 and 51904512, is three
ranges and 1572864 bytes, and so is every other size that has been read: 4718592 at 1 MiB in the table above,
7864320 at 2 MiB, 20447232 at 6 MiB and 26738688 at 8 MiB, each of them
3 * the size of a range + 1572864 to the byte. The 2, 6 and 8 MiB are of the
code review of this work package, on 24 September 2026, and the 6 MiB was
read again by the browser test with `NUM_BYTES_PER_RANGE` set to 6291456,
which gave the same 20447232. The 256 KiB of the table is the one size that
does not follow that count: 3407872 bytes, 1048576 more than it gives.

Of the 1572864 that does not grow with the range, 1310720 is the module
before any pass, read after `init()` alone, which leaves 262144 for
everything the VCF reader holds.

Two of the three range-sized blocks are the pass's. The function of
`crates/popnei-js/src/source.rs` that reads a range,
`RangesOfAFile::reads_the_range`, copies the bytes the browser gave it into
a new block of the memory of wasm and puts that block in the place of the
range the pass held, which is freed after that, so the two are live for that
moment. The third is the range that the reading of the header allocated at
`openVcf`, whose reader is dropped when `openVcf` returns: the block is
freed there, and the first range of the pass is not given it. A worker of
24 September 2026 read `memory.buffer.byteLength` after each step of one run
over a `File` of 14008097 bytes, with the range at 4 MiB: 1310720 when the
module was loaded, 5636096 after the first `openVcf` over that file,
9895936 after a second `openVcf` over the same file, 9895936 after a third,
14155776 after one pass, and 14155776 after each of three more passes, the
last of them over a fourth source opened once the other three were freed.
So a freed range is handed to the second request that follows it and not to
the next one, and past three blocks the heap serves every pass without
growing.

Copying each range into a buffer the pass keeps, instead of building a new
one every time, would take one of the three blocks away; it is a change of
the code that this task does not make, and 14155776 bytes is 4.7 % of what
the whole file costs.

The bound of the browser test of task 4.2 is 16777216 bytes, 16 MiB, 1.19
times what was measured. It leaves room for a consumer that holds a little
more than `calcPerIndividualStats` does, and it fails both at one range more
held alive, 18350080 bytes, and at a range of 6 MiB, 20447232, which is the
run above. The file of that test, 299994147 bytes, is 17.9 times the bound,
where the plan asks for a file of at least eight times it.

## Why 4 MiB

The times do not choose the size of a range. 1 MiB, 4 MiB and 16 MiB took
1006, 996 and 958 ms of a pass of about a second, 48 ms between the fastest
and the slowest of the three, and the five runs of one point spread by as
much as 73 % on this machine. 4 MiB was faster than 1 MiB in each of the
four rounds, by 10, 12, 15 and 61 ms, which is about 1 % of a pass.

What chooses 4 MiB is the calls into the browser. A pass over this file
makes 287 of them at 1 MiB and 72 at 4 MiB, and every byte of this
measurement came out of a `Blob` built in the memory of the worker, so what
a call costs when the file is on a disc was not measured: a range of 1 MiB
pays that unmeasured cost four times as often as one of 4 MiB does.
Choosing by the time and the memory alone would not give 4 MiB. 16 MiB is
refused here because 38 ms of a pass, 3.8 %, is not worth 51904512 bytes of
the memory of wasm against 14155776, 3.7 times more; by that same trade
1 MiB beats 4 MiB, since it costs 10 ms, 1.0 %, and leaves 4718592 bytes,
3.0 times less.

What the memory buys is what the applications came for. A tab that reads a
file of 2 GB, which is above what a wasm module can hold whole, pays
14155776 bytes for it at 4 MiB and 51904512 at 16 MiB. A source with several
passes open at once holds the ranges of each of them: the twelve passes of
`js/popnei/test/vars_memory.test.ts` hold at least one range apiece, 48 MiB
of ranges at 4 MiB and 192 MiB at 16 MiB, and what twelve passes leave in
the memory of the module together has not been measured.

Whether an application chooses the size of a range instead is the third open
point of `docs/specs/js_sources.md` and is still the owner's to decide.
Nothing of this measurement asks for the option: what a smaller range buys
is a progress bar that moves more often, and what it costs is time. 256 KiB
tells the bar four times as often as 1 MiB does, and a pass at 256 KiB took
1069 ms over this file where the whole file in memory took 921 ms, 148 ms
more.

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
| `cargo clippy --workspace --all-targets -- -D warnings` | `Finished \`dev\` profile [unoptimized + debuginfo] target(s) in 0.38s` |
| `cargo test --workspace` | `787 passed; 0 failed; 2 ignored` and `149 passed; 0 failed` |
| `cargo test -p popnei --no-default-features` | `787 passed; 0 failed; 2 ignored` |
| `cargo wasm-check` | `Finished \`dev\` profile [unoptimized + debuginfo] target(s) in 0.09s` |
| `cargo wasm-check-js` | `Finished \`dev\` profile [unoptimized + debuginfo] target(s) in 0.45s` |
| `uv run ruff format --check && uv run ruff check` | `33 files already formatted`, `All checks passed!` |
| `uv run maturin develop && uv run pytest` | `499 passed in 10.84s` |
| `npm run build && npm test` in `js/popnei` | `tests 372`, `fail 0` |
| `npm run test:browser` in `js/popnei` | `8 passed (1.8s)` |

Seven of the eight browser tests are the ones work packages 2 and 3 left;
the eighth, which asserts the bound above, is task 4.2.

## The commits

The measurement and the two numbers it set:

- 91cf7ba: the "Speed" part of `docs/specs/js_sources.md`, which carries the
  two numbers, and the sentence of "How it runs" that said a pass holds one
  range.
- 9233fc5: the doc comment of `NUM_BYTES_PER_RANGE` in
  `crates/popnei-js/src/source.rs`, whose value does not change, and the
  script, the page and the worker that took the times.
- edb9001: the browser test over a file of 299994147 bytes.

What the code review of work package 4 sent back, and this report with it:

- a53bb9e: "Speed" and "How it runs" of the spec, with the three
  range-sized blocks in the place of "3.4 times its range", the bound at
  16777216 and the reason the size of a range is 4 MiB.
- 7db0a87: the bound of the browser test, the doc comments of the case, of
  `NUM_BYTES_PER_RANGE` and of `RangesOfAFile`, the rebuild of the
  WebAssembly in `bench/time_ranges_in_chromium.mjs`, and the one function
  that builds the `File` of `many.vcf` repeated for the two cases that
  open one.
