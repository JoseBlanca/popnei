# Work report: a file of a page read by ranges in the wasm package

24 September 2026. The work of `docs/plans/js-sources.md`, on the branch
`plan/js-sources`, which builds what issue 1 of popnei asks for: the wasm
package reads a file the user picked in the page one range of bytes at a
time, instead of taking the whole of it into the memory of WebAssembly.

**The plan is done.** Every task is carried out, every deliverable was
checked by the orchestrator with the command the plan gives, and each of its
four work packages was reviewed, by seven reviewers for the first and by
four, four and two for the others.

What exists now that did not:

- `openVcf` and `openVars` take a `File` or a `Blob` as well as an array of
  bytes, and read it 4 MiB at a time through `FileReaderSync`, the call of a
  web worker that reads a range of a file and returns when it has the bytes.
  A pass over
  a VCF of 299994147 bytes takes 996 ms in Chromium against 921 ms with the
  file in the memory of wasm whole, 8.1 % more, and holds 14155776 bytes of
  that memory against 302383104, 21.4 times fewer. So what limits the size
  of a file a user can open is time and not memory.
- `variants.onProgress(fn)` calls a function of the application with the
  bytes the pass has read, the size of the file, which pass of the run is
  reading and how many passes it makes. A function that throws ends the run
  and the consumer throws back the application's own value, so a cancelled
  analysis keeps its worker, where today it ends it and the new worker
  reads the file again.
- `numPassesOf` says how many passes a consumer will make before it starts,
  from the same place the calls above take that number.
- Eight tests run popnei in Chromium, through Playwright, where nothing of
  popnei had ever run in a browser. `npm run test:browser` in `js/popnei`
  takes 1.8 s and builds the wasm it tests.
- The package is smaller than before the work: 2021550 bytes of wasm and
  650486 gzipped, against 2034662 and 657075.

What is asked of the owner:

1. **The merge**, which is what I recommend: the branch is 57 commits, the
   checks of the `coding` skill pass on it, and `npm run build && npm test
   && npm run test:browser` in `js/popnei` passes from a tree with nothing
   built in it, which is what an application that installs the package
   gets. Nothing of this is in `main`, and the branch holds the spec
   `docs/specs/js_sources.md` as well, which `main` does not have: it
   starts from `spec/file-source`, so merging it brings the spec and the
   work together. What argues against merging today is Open 2 of the spec:
   if you answer it the other way, what a stopped run throws changes, which
   is a value applications see.
2. **The three open points of `docs/specs/js_sources.md`.** They are yours
   to answer, and the code follows the recommendation of each, which you can
   now see working. Answering them as they are recommended costs nothing and
   settles the spec. Answering the first differently, where the progress
   function is set, changes `variants.onProgress` and the code of task 1.2;
   the second, what a stopped run throws, changes what every consumer throws
   and the code of task 1.3; the third, whether an application picks the
   size of a range, adds an argument to `openVcf` and `openVars`. Doing
   nothing leaves three sentences of the spec marked open and the code where
   it is.
3. **What this work found and did not touch**, filed as issue 2 of the
   repository, https://github.com/JoseBlanca/popnei/issues/2: popnei
   refuses to read a vars file it wrote itself when its genotypes compress
   below one bit per value. Nothing of this plan touches the reader of that
   file, and doing nothing leaves a user who converted a VCF and deleted it
   without their variants.
4. **Section 11 of `docs/architecture.md`** says every error of the core is
   thrown as a JavaScript `Error` with its message. A run that an
   application stopped now throws the application's own value, which need
   be no `Error` at all. That line follows from Open 2 and is left for the
   owner's answer to it.

What was left out, and where it goes: Firefox and WebKit in the tests, which
the owner left for when an application needs them; the other three requests
of the applications, the fingerprint of a variant file, the reader of CSV
and TSV and the inference of the types of its columns; and stopping a run
while it is not reading, which would be a request to the core crate.

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
and `js/popnei/test/before_init.test.ts` is where the package tests that for
its other functions. The orchestrator gave that test to task 1.2, which
writes test files of its own.

## Work package 2: the browser harness

Done at aea65cf, in one task. `npm run test:browser` in `js/popnei` runs
Playwright 1.63.0 with Chromium headless, which starts a server over the
repository, a page and a module web worker that loads the wasm, and the test
that opens `many.vcf` from an array of bytes inside that worker and asserts
its first block. It took 3.5 s on this Mac, 1 test. `npm test` is unchanged
at `fail 0`, and `node --test test/` does not pick up the files of
`test/browser/`, which the task checked.

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
JavaScript with node's `stripTypeScriptTypes`, which is experimental in node
v26.8.2 and prints one warning per run. What it buys is that the page, the
worker and the cases are type checked with the rest of `js/popnei/test/`.

Task 1.2, the run and the counting source, is done at 7c2bcf8. `node --test
test/progress.test.ts` runs 10 tests where deliverable 2 asks for 6, and
`npm test` gives 350 tests, `fail 0`. The size of a range is one constant, 4
MiB, whose doc comment says work package 4 measures it. `js-sys` came in
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

A `variants.onProgress(x)` where `x` is neither a function nor nothing is
refused there, where the application wrote it, and not at the first read of
the next run. The spec did not say so and now does.

Task 1.3, the stop, is done at bc46752. `node --test test/stop.test.ts` runs
7 tests where deliverable 3 asks for 4, five of which fail on the commit
before it. The error type of the binding crate gained a case that carries
the value the application threw through untouched, the read that a throw
ends fails with `std::io::Error::other`, and nine of the ten functions that
read a source open their run in one function of the crate, so the swap of
the error cannot be forgotten. The tenth is the iteration a user writes,
`iterBlocks`, which reads between one call and the next, so its run lives
longer than any one call and is opened where the iteration starts.

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
spec follows them at 7606b76. The deliverables, each run by the orchestrator
after the fixes:

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

## Work package 3: the source over a `File`, with the harness of work package 2

Done at 8437da5, c580acf and 929f79d, and reviewed together with work
package 2, whose harness its browser tests are written in. The deliverables,
each run by the orchestrator:

| deliverable | command | what it gave |
|---|---|---|
| 1, 2 | `npm run test:browser` | 7 tests in Chromium, 0 failed: the array of bytes of work package 2, a `File` of `many.vcf`, of `many.vcf.gz` and of a vars file written from it, a file of more than one range, a range that comes back short and a range the browser refused |
| 3 | `node --test test/blob_outside_a_worker.test.ts` | 6 tests, 0 failed |
| 4 | `grep -c FileReaderSync js/popnei/README.md` | 6 |
| 5 | the checks of the `coding` skill | all green, with `cargo wasm-check-js` added to them |
| 5 | `npm test` in `js/popnei` | 372 tests, 0 failed |

The package is smaller with the ranges than without them: 2021550 bytes of
wasm against 2034662, and 650486 gzipped against 657075. `js-sys` and
`web-sys` are in it; what outweighs them is that the read which opens a file
goes through the same type as every other read, which took one build of the
VCF reader and one of the vars file reader out of the module.

### What the review found

Four reviewers, and two of them found the same hole: no test read a file of
more than one range, because every file of every test is smaller than the 4
MiB one holds. One of them patched the reader so that the end of the first
range was the end of the file, the silent wrong result the short range is
refused for, and all five browser tests passed. The spec now asks for that
check and `many_ranges.ts` is it, a `File` of 14008097 bytes and 60000
variants.

The other finding that cost a user something: reading the counts of a pass
from inside the progress function leaked the pass and the bytes of its file,
6.7 MB a run, for the same reason freeing the variants there did. Both are
refused now with popnei's own error.

Nine findings were taken in all. One was not: a read that fails names the
range and the size of the file and not the name of the file, which a `File`
of a page has. The spec asks nothing there, the coding skill says JavaScript
has no path to add, and a worker of the applications opens one variant file
at a time, so what it would add is a feature of `web-sys` for a message that
already names the range.

## Work package 4: the size of a range, measured

Done at 91cf7ba, 9233fc5, bc2d2c6 and edb9001, and reviewed by its spec and
its tests. `docs/reports/js-sources-measurement.md` has the tables; the two
numbers it sets are the size of a range, which stays at 4 MiB, and what a
pass holds in the memory of wasm at that size, 14155776 bytes.

| deliverable | command | what it gave |
|---|---|---|
| 1 | `node js/popnei/bench/time_ranges_in_chromium.mjs` | the times of one pass at four sizes and of the whole file, in `docs/reports/js-sources-measurement.md` |
| 2 | the constant and the spec | `NUM_BYTES_PER_RANGE` at 4 MiB with its measurement, and "Speed" of the spec with both numbers |
| 3 | `npm run test:browser` | 8 tests, 0 failed, the last a `File` of 299994147 bytes read with the memory of wasm under 16 MiB |
| 4 | the checks of the `coding` skill | all green |

One pass over a VCF of 299994147 bytes, 1285000 variants, with
`calcPerIndividualStats`, best of five runs in Chromium on the owner's Apple
M5 Pro: 996 ms at a range of 4 MiB against 921 ms with the file in the
memory of wasm whole, 8.1 % more, holding 14155776 bytes of that memory
against 302383104, 21.4 times fewer. 256 KiB takes 1069 ms, 1 MiB 1006 and
16 MiB 958.

### What the review found

The memory a pass holds is `3 * range + 1572864` bytes, to the byte, from 1
MiB up, which both reviewers measured before either of them knew why. The
three blocks are the two ranges a pass holds while it builds the next and
the range the header read of `openVcf` allocated and freed, which the
allocator does not give back for the next request of its size. The browser
test asserts that the memory of wasm stays under a bound while it reads a
file of 300 MB; that bound was 24 MiB, 1.78 times the 14155776 bytes
measured, which is loose enough that a build holding a range of 6 MiB passed
it. It is 16 MiB now, which fails at one range more.

The choice of 4 MiB did not follow from the table as the report first argued
it: by the trade that refuses 16 MiB, 1 MiB would beat 4 MiB: it makes a
pass 1.0 % slower, 10 ms of 996, and holds three times less memory. What
keeps 4 MiB is that every byte of the measurement came from a file built in
the memory of the browser, so what one call costs when the file is on a disc
was not measured, and 1 MiB makes 287 of those calls over a file of 300 MB
where 4 MiB makes 72.

The bench left `js/popnei/wasm/` built at the last size it measured while
saying it put everything back, so the next run of the browser tests failed
at a case that names nothing about the bench. It rebuilds now, on an
interrupt as well.

## How the work went

The owner can stop here: this section is for whoever next revises a skill or
writes a plan.

**A reviewer checked out the commit under review in the shared worktree,
and two commits went off the branch.** The `code-review` skill sends the
`spec` and `tests` reviewers with a worktree of their own and tells them to
check the commit out; the other five read the tree that the work is going on
in. One of those five checked out the reviewed commit there to rebuild the
wasm from it. The branch was left pointing at an older commit, the
orchestrator committed twice on a detached head without noticing, and the
two were put back with a cherry-pick and one conflict resolved by hand. The
skill should say, where it says which reviewers get a worktree, that a
reviewer without one runs `git checkout` nowhere: what it reads is the tip,
and a build it needs it makes in place. The prompts of every later review
said so and no reviewer did it again.

**The reviews found more than the writing did, and the numbers say by how
much.** The seven reviewers of work package 1 cost 991000 tokens against the
682000 of the three tasks that wrote it, and they found two defects a user
would have met, each with a measurement: a leak of 120 KB per source and a
progress bar that stopped at two thirds. The reviewers that ran the code
found those; the ones that read it found stale prose. Both kinds of finding
are worth what they cost, and the second kind is what a first reader of a
document would have caught.

**Six sentences of the spec were wrong, and every one of them was found by
building what it described.** They were not careless: each was a claim about
what popnei's own readers do, written from reading them. A reader reads its
header when it is built, so the calls of the two passes of a PCA interleave;
no read of a vars file finds the end of the file, so no call said a pass was
over. The writing of a spec cannot check those; the first task that builds
from it can, and the plan should expect it. What worked was the rule that
the spec changes first, in a commit of its own: the history now says what
was believed and when it stopped being believed.

**A measurement that nobody questions is a number without a rule.** The
report of work package 4 kept 4 MiB with an argument that, taken
symmetrically, chose 1 MiB. The reviewer of the spec found it by applying
the report's own trade to another row of its own table. A measurement that
picks one of several numbers says which rule it picks by, or it has not
finished.

**What a task of this size costs.** The tasks that wrote code cost between
130000 and 284000 tokens each; the reviewers between 87000 and 182000; the
two rounds of fixes 240000 and 172000. The whole plan, with its reviews and
its fixes, is about 3.6 million tokens of subagents. The tasks that were one
clear piece of work came back once and right; the two that were "the
findings of the review" needed a second message from the orchestrator when a
reviewer reported late, which is a reason to hold a fix until every reviewer
of its work package is back.
