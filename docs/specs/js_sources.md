# The sources of the wasm package: a file read by ranges

24 September 2026. A web application built on popnei runs in a browser tab,
where the user picks a file in the page and the calculations run in a web
worker, a thread of the page that cannot touch what the page shows. This
spec says how popnei reads that file: one range of bytes at a time, as its
readers ask for them, instead of taking the whole file as an array of bytes;
what the source tells the page while it reads, so that the page can show how
far a pass has got and can stop it; and how many passes over the file each
consumer makes. It develops section 11 of `docs/architecture.md`, the
TypeScript boundary, and it is what issue 1 of popnei asks for, which the
owner made a priority on 24 September 2026. On 7 October 2026 it gained
two items from issue 10, written for the page of popnei_web, the web
application of popnei, that shows the statistics of a file as soon as it is
opened: "The result so far", which a calculation gives the page while its
pass runs, and "The three statistics of a file in one pass", a consumer that
reads the file once for what three consumers read it three times for. There
is no code of either. Later that day it gained "The file of a writer in
pieces", from issue 13, an option of `writeVars` and `writeVcf` that gives
the file to the page as it is written.

It changes the JavaScript binding crate, `crates/popnei-js`, and the
TypeScript package, `js/popnei`, both of which are written and read a file
that is already in memory. The core crate is not touched: its VCF reader is
generic over `BufRead` and its vars file reader over `Read + Seek`, as
section 1 of the architecture asks, so a new source of bytes reaches both
unchanged. The specs of those two readers, `docs/specs/io_vcf.md` and
`docs/specs/io_vars.md`, say what each of them reads and what it refuses;
`docs/specs/variant.md` has the counts of a pass, which the progress of this
spec stands beside; `docs/specs/block.md` has `iterBlocks`.

Three words of `docs/glossary.md` are used throughout. A **pass** is one
reading of a source of variants from its start to its end. A **consumer** is
what takes a `Variants`, makes the passes it needs and gives a result: a
calculation, the writers `writeVars` and `writeVcf`, or `iterBlocks`; the
package has fifteen of them, `writeVcf` among them since 26 September 2026 and `calcVarDensity` and `doPcoaFromVariants` since 27 September 2026. A **run** is one call of one consumer, with the passes it makes.

What the user of an application pays today, with the file taken whole:

- The file is in the memory of wasm whole, and about twice its size while it
  opens, the worker's copy and wasm's. A wasm module addresses 4 GB, so a
  file above roughly 1.5 to 2 GB cannot be opened at all. That figure is an
  estimate from those sizes; no browser has been measured at it.
- The memory of wasm grows and never shrinks, so the worker holds the file
  until it is restarted.
- A run that the user cancels ends its worker, which is the only way to stop
  a loop inside wasm: a page tells a worker something by sending it a
  message, which the worker reads when it is between two messages and not
  while it is inside a calculation, and the one memory a page and a worker
  can share and read at the same time, a `SharedArrayBuffer`, needs headers
  that GitHub Pages does not send. The new worker then reads the whole file
  again before the next run starts.
- Nothing can be shown of how far a run has got, only that it is running.

## The source over a file of the page

### What it gives

A `File` is the handle a page gets when the user picks a file: it carries
the name and the size of the file and gives any range of its bytes, and it
costs nothing to hold or to send to a worker, because it is a handle and not
the bytes. `FileReaderSync` reads a range of a `File` and returns when it
has the bytes, which is what Rust's `Read` and `Seek` need; it exists only
inside a web worker.

`openVcf` and `openVars` take a `File` as well as an array of bytes. With a
`File`, popnei reads the ranges it needs through `FileReaderSync`, a few MiB
at a time, so the file is never in the memory of wasm whole: what a pass
holds is the range it is reading, the range before it while the new one is
being built, the block it is building and, over a vars file, the batch it is
reading. "Speed" below has the bytes a pass over a VCF came to. The size of
a file a user can open stops being limited by the memory of the tab, and a
worker that was restarted opens the file again instead of reading it.

Each pass reads the file again from its start, as a pass over an array of
bytes does, which is what lets a user give one `Variants` to one consumer
after another.

### Its TypeScript function

```ts
export type BytesOrFile = Uint8Array | Blob;

openVcf(source: BytesOrFile, options?: OpenVcfOptions): Variants
openVars(source: BytesOrFile): Variants
```

A `File` is a `Blob`, the type of a piece of bytes of a page, so a `Blob`
that an application built itself is read the same way. Everything else of
the two functions is unchanged: they read the header of the VCF or the
schema and the footer of the vars file when they are called, so bytes that
are not of the format fail there, and they give the `Variants` of
`docs/specs/variant.md`.

There is no Python function. Python opens a file by its path and has no
`File`; `docs/specs/io_vcf.md` and `docs/specs/io_vars.md` have its side.

### The cases a reader of the rules would not guess

`FileReaderSync` is not defined outside a web worker, so `openVcf` and
`openVars` with a `Blob` on the main thread of a page throw an `Error` that
says so and names the worker. It happens at the call, because the header is
read there. An array of bytes works on the main thread as it does now.

A range that comes back shorter than the range that was asked for, inside a
file of that size, is an `Error` and not the end of the file. So is a
`Blob` whose `size` is not a whole number of 0 to 2^53: it is a number of
JavaScript, and what a file holds has to be counted exactly before a range
of it is asked for. A browser
gives a short range when the file changed on disk after the page got its
handle, and a reader that took it for the end would give the variants it had
and say nothing. What each browser does with a file that changed has not
been tested here.

A pass over a vars file holds one batch: the reader reads the bytes of a
whole batch, the ones its footer says it holds, and hands them to arrow-rs,
the library popnei reads and writes that file with, which decompresses the
columns the pass asked for. At the size popnei writes, a batch is about
10 MB of genotypes for 1000 individuals. So the file is no longer in memory
whole and a batch of it is.

A `Variants` over a `File` keeps the handle until its `free()` is called,
and the handle is not the bytes: what the page holds for an unfreed
`Variants` is the name and the size of a file. A pass that is still open
when `free()` is called reads on to its end, as a pass over an array of
bytes does today, because the handle is kept until the last pass over it is
done.

### How it runs

The bytes of a pass come from one of two places, and one Rust type reads
from both: a copy of the whole file in the memory of wasm, which every pass
over it shares, as today; or a `File`, one range at a time. That type
implements `Read`, `BufRead` and `Seek`, so the VCF reader and the vars file
reader take it where they take the cursor over an array of bytes today, and
the counting and the call to the page of the item below happen once, for
both.

A reader has to be `Send`. `BlockReader`, the trait of section 1 of the
architecture that everything which gives blocks implements, asks for it, so
that a read ahead thread can move a reader to another thread natively. No
handle of JavaScript is `Send`. So the handles of a source, the `Blob`, the
`FileReaderSync` over it and the function the page is told the progress
with, live in a table of the binding crate that no thread leaves, and the
Rust type holds the number of its entry. A spike of 24 September 2026
compiled `VcfReader<PassOverTheBytes>` and `VarsReader<PassOverTheBytes>`
for `wasm32-unknown-unknown` and for this Mac, and both are `Send`. What
this avoids is an `unsafe impl Send`, which the lints of the workspace deny
in this crate, and a change to the trait of the core.

Each range crosses once into the memory of wasm, copied out of the
`ArrayBuffer`, the block of bytes JavaScript gives back, that
`FileReaderSync` fills. A pass holds two ranges while it reads: the new one
is built and put in the place of the one before it, which is freed after
that. What the module holds when a pass over a VCF has ended is three
range-sized blocks and 1572864 bytes that do not grow with the range, of
which 1310720 is the module before any pass and 262144 everything the VCF
reader holds, the lines of a block among them; a pass over a vars file holds
its batch as well, which nothing has measured. The
third range-sized block is the range that the reading of the header
allocated at `openVcf`: it is freed when `openVcf` returns, the memory of a
wasm module never shrinks, and the first range of the pass is given a block
of its own beside it. "Speed" below has the byte counts and the sizes they
were measured at. A source with several passes open at once holds the ranges
of each of them; a test of the package keeps twelve passes over one vars
file open together, and what twelve of them hold at once has not been
measured.

The two crates that call `FileReaderSync` and `Blob.slice` are `js-sys` and
`web-sys`. They are pure Rust and they compile for
`wasm32-unknown-unknown` and for the native target that `cargo clippy
--workspace --all-targets` builds this crate for. What they cost the
download of an application, measured on the spike above with `wasm-bindgen
--target web --remove-name-section` and `gzip -9` on the owner's Apple M5
Pro: the wasm of the package went from 1979013 to 1982033 bytes, and from
633122 to 633898 gzipped, 776 bytes more.

The package built with the ranges is smaller than the one without them, by
the same measurement: 2021550 bytes against 2034662, and 650486 gzipped
against 657075, 6589 fewer. The two crates are in it, and what outweighs
them is that the read which opens a file goes through the same type as
every other read, where it went through a cursor of its own before: that
took one build of the VCF reader and one of the vars file reader out of
the module.

### How it is verified

`FileReaderSync` exists only in a worker, so node cannot run a test of it.
Everything that is not the reading of a range is written to work the same
over an array of bytes, and is tested under node: the counting of the bytes,
the call to the page, the stop, and the count of the passes of each
consumer. What the browser test adds is that the ranges of a real file give
the same variants and that a file of 299994147 bytes is passed over with the
memory of wasm staying under the bound of "Speed" below.

Under node, with `tests/reference/vcf/many.vcf`, 500 variants of 50
individuals in 117346 bytes:

- `openVcf(new Blob([bytes]))` throws an `Error` whose message names
  `FileReaderSync` and the web worker. node v26.8.2 has `Blob` and `File`
  and no `FileReaderSync`, which is the case of the main thread of a page.
- Every test of the package that opens an array of bytes keeps its literals
  unchanged, which is what says that the counting changed no number.

In a browser, Chromium through Playwright, the program that drives a browser
from a script, headless:

- A worker opens a `File` of the bytes of `many.vcf` and asserts the first
  block, the individuals, the positions and the genotypes, against the
  literals a node test of the package already asserts over the array of
  bytes. For that file they are in `test/filter_individuals.test.ts`, which
  reads its positions and its genotypes from `docs/specs/io_vcf.md`, and not
  in `test/vcf.test.ts`, which asserts of its first block only that it holds
  100 variants. The gzipped `many.vcf.gz` gives the same variants, which is
  the path through the decompressor.
- The same worker writes a vars file with `writeVars`, makes a `File` of it,
  opens it with `openVars` and asserts the same variants. The vars file
  reader seeks, so this is the check of `Seek`; the VCF reader never seeks.
- A `File` of more than one range gives the variants of every range of it.
  Without such a file no test reads a second range: `Blob.prototype.slice`
  giving the end of the first range as the end of the file passed every
  other browser test of this spec, which is the silent wrong result the
  short range above is refused for.
- `Blob.prototype.slice` patched to give one byte less than it was asked
  for: the pass fails with popnei's error, and the message names the range
  it asked for, the bytes it was given and the size of the file.
- A `File` of a few hundred MB, the body of `many.vcf` repeated, is opened
  and passed over once with `calcPerIndividualStats`, and the test asserts
  the number of variants it counted and that the memory of the wasm module
  never passed a bound. It reads `memory.buffer.byteLength`, as
  `js/popnei/test/vars_memory.test.ts` does under node. The bound is what
  says the file was not read whole, and the file is at least eight times it.
  Both numbers come from the measurement of "Speed" below, which the first
  task of the plan makes.

The other two engines, Firefox and WebKit, which Safari is, are not run yet:
the owner decided on 24 September 2026 to start with Chromium and add them
when an application needs them. Playwright drives all three, so what is
added is a line of its configuration and their download.

## What the source tells the page

### What it gives

While a run is going on, the worker is inside wasm and reads no message, so
the page learns nothing about it. The source is the one part of a pass that
comes back out to JavaScript, once per range, and it is where both of these
are done: it tells the page how far the pass has got, and it takes from the
page the decision to stop.

The application gives a function, and the source calls it with how many
bytes of the file the pass has read, how many the file holds, which pass of
the run is reading and how many passes the run makes. Three things make a
call: the first read of a pass, which says it has read nothing; a read that
brings the bytes read since the last call to the size of a range; and the
end of the run, which makes one for each of its passes, in the order of
their numbers, with the bytes that pass read. The last of the three is what
says a pass is over, and it is needed because no read does: a pass over a
vars file stops after its last batch, with up to a range of bytes read
since the last call, and a run that fails stops wherever it failed. A value
thrown in one of those last calls is dropped and stops nothing, because the
run is over. A page that draws a bar from these numbers sees it fill once
per pass and knows which pass it is on, so a PCA that reads the file twice
does not look broken when the bar goes back to empty.

When that function throws, the pass ends there: the read fails, the error
travels out through the readers, and the consumer throws the value the
function threw. So an application cancels a run without ending its worker,
and what it gets back is its own object, which it recognises without reading
a message (**Open 2**, below). Nothing of popnei is left in a state a later
call notices: the `Variants` is the one it was, and the next run over it
reads the file from its start.

The stop happens while the file is being read, and not inside the linear
algebra. A kinship of 10000 individuals reads the file once and then spends
its time in an eigendecomposition that reads nothing, and a cancel during
that part still ends the worker. A call from inside the loops of the core
would be a request to the core crate, and no application has asked for one.

### Its TypeScript function

```ts
export interface Progress {
  /** How many bytes of the file this pass has read, `numBytes` at most. */
  bytesRead: number;
  /** How many bytes the file holds. */
  numBytes: number;
  /** Which pass of the run is reading, 1 for the first. */
  pass: number;
  /** How many passes the run makes, `numPassesOf` of its consumer. */
  numPasses: number;
}

variants.onProgress(told?: (progress: Progress) => void): void
```

The function is set on the `Variants` and not given to each consumer
(**Open 1**, below), and it holds until it is set again; `onProgress()` with
nothing takes it off, and anything that is neither a function nor nothing is
an `Error` there, where the application wrote it, and not at the first read
of the next run. Setting it changes nothing about the variants a pass
gives.

Every consumer of the package can throw the value that `told` threw:
`calcPerVarDistribs`, `calcPerIndividualStats`, `calcVarDensity`,
`calcPairwiseKosmanDists`, `calcPopDists`, `calcPopDiversity`, `calcRogersHuffR2Matrix`,
`calcLdAndDistPerPop`, `calcKinship`, `doPcaFromVariants`, `calcGwas`,
`writeVars`, `writeVcf` and the iteration of `iterBlocks`.

`onProgress` has no Python counterpart, and neither has `numPassesOf` of the
item below. Goal 2 of `docs/objectives.md` asks for every difference between
the two APIs to be written down, and this is one: what they are for is a
page that draws a bar and a user who presses a button, and Python reads a
file by its path in a program that has neither.

### The cases a reader of the rules would not guess

The bytes read are the bytes that pass gave to its reader, and they are
counted against the size of the file as it is on disk. For a gzipped VCF
that is the compressed bytes: a pass over a VCF of 21904 bytes gzipped ends
at 21904 and not at the 117346 bytes of its text.

A pass over a vars file does not read the file evenly, and it does not read
all of it. It reads the last ten bytes first, which say how long the footer
is, then the footer, which says where the batches are, and then each batch
whole, one after another; the schema message at the head of the file it
never reads, because the footer carries the schema too. Measured on
`tests/reference/vars/zstd.vars`, 2418 bytes for three individuals: the
footer is 784 bytes and the schema message 624, which is a quarter of that
file and a smaller share of a real one, where the batches are most of the
bytes. So a bar over a vars file ends below the size of the file. The count
is capped at that size, so a bar never passes 100 in 100.

A file of fewer bytes than one range gives two calls, at the first read of
the pass and at the end of the run.

`openVcf` and `openVars` read the header or the schema before any pass, and
those reads are told to nobody: what the function is set on is the
`Variants` that those calls give.

A function that throws is called no more in that pass, and the call that
ends that pass is not made either: what the application stopped it is not
told how far it had got. The other passes of that run are told, with what
they had read: a PCA stopped at the first call of its second pass tells the
page of pass 1 at the 617 bytes of the header of `many.vcf`, which is all
it had read when both readers were built.

A pass whose reader is never built is never told of, and a pass whose
reader is built always is: a reader reads when it is built. So the PCA of
a source whose first pass gives no variant with variance, which fails
before its second reader asks for a block, has told the page of both
passes; and asked for no weight, which builds no second reader, it tells
the page of one.

Passes are counted inside a run and not inside a source. A consumer opens a
run, every reader it opens belongs to it, and a pass takes its number when
it first reads. Twelve `iterBlocks` over one source at once are twelve
runs, each of one pass, and each of them is pass 1 of 1.

The calls of the two passes of the PCA of the variants are not the calls of
one and then the calls of the other. It builds both of its readers before
it asks either for a block, and a reader reads when it is built, the header
of a VCF and the footer of a vars file, which is the first read of its pass
and the call that says it has read nothing. So the two take their numbers
there, one after the other, and each reads the file after that: over
`many.vcf` the four calls of a run are pass 1 at 0 bytes, pass 2 at 0,
pass 1 at 117346 and pass 2 at 117346. What a page draws rises all the
same, if it draws the share of the run that is done,
`(pass - 1 + bytesRead / numBytes) / numPasses`: 0, 0.5, 0.5 and 1.

A pass that the filter of the first n variants of `docs/specs/filters.md`
ends reads no further, so its last call has a `bytesRead` below `numBytes`,
the bytes it had read when the filter ended it, and no later call brings it
to `numBytes`: what the page draws stops short of the end, and the
`stoppedEarly` of the counts of the pass says why. The page is told every 4
MiB, so the test is over a VCF of more than twice that, which the test
writes: with the first 10 the last call has a `bytesRead` below its
`numBytes`, and `bytesRead` never goes down. Added on 5 October
2026.

An `iterBlocks` that a user abandons without freeing it holds its run, and
with it the entry of its source, until the `FinalizationRegistry` of the
package frees the pass, at a moment nobody chooses. It is the rule the whole
package already has for the memory of wasm, which
`docs/specs/block.md` states: what a user does not free is freed late.

The function is called with no table of the binding crate borrowed, so an
application that calls popnei from inside it does not trap, and a consumer
started there runs as any other call does. Two calls of it do not go
through. `free()` is the first: a `Variants` counts the consumer calls
over it that are running, and while one is, `free()` throws popnei's
`Error` naming that a run is reading it and frees nothing. That throw is
the application's own, since the application made the call, so it stops the
pass as any other throw of that function does. What it prevents is the
memory of the file being lost: wasm-bindgen gives up the pointer of a
source and unregisters it before it calls into wasm, and that call throws
while a consumer holds the source, so the file would be unreachable and
never freed for as long as the page lives, 100 sources of 117346 bytes
costing 11993088 bytes of the memory of wasm. Inside an iteration of
`iterBlocks` no call holds the source, the free is taken, and the pass
reads on to its end, as the paragraph on `free()` of the item above says.
The second is the counts of a pass of `iterBlocks`, which are read between
two blocks and not while one is being read: the pass is held for the length
of that read, so popnei refuses them with its own `Error` naming it. What
both refusals prevent is the same, and it is not a message: the handle that
wasm-bindgen was asked for is given up before the call it fails in, so the
source or the pass, and the bytes they hold, would be unreachable and never
freed. Measured on a pass over a VCF of 7004357 bytes whose function read
those counts: 20 runs grew the memory of wasm by 133365760 bytes, where 20
that did not grew it by none.

### How it is verified

Under node, over `many.vcf`, its gzip `many.vcf.gz` and a vars file written
from it, with ranges of the size popnei chose:

- The calls of one pass of `calcPerVarDistribs` over `many.vcf`: the first
  has `bytesRead` 0, the last has `bytesRead` 117346, every one has
  `numBytes` 117346, `pass` 1 and `numPasses` 1, and `bytesRead` never goes
  down. Over `many.vcf.gz` the last call has `bytesRead` 21904 and
  `numBytes` 21904.
- Over a vars file of more than one range, written by the test as
  `vars_memory.test.ts` writes one, the first call is 0, `bytesRead` never
  goes down and never passes `numBytes`, and the last call is the bytes of
  its footer and its batches, which goes into the test as a literal when it
  is first run, with the file it was measured on. That last call is the one
  the end of the run makes: without it a pass over a vars file of 12231602
  bytes was last told at 8476400, two thirds of the way, because its reads
  after that never reached a range.
- Over a vars file of fewer bytes than one range there are two calls, of 0
  bytes and of the bytes of its footer and its batches.
- `doPcaFromVariants` with `numPrinComps` 10 over `many.vcf` gives four
  calls, of pass 1, pass 2, pass 1 and pass 2, with 0, 0, 117346 and 117346
  bytes read, and `numPasses` 2 in every one. With `numPrinComps` 0 there
  is one pass.
- Twelve `iterBlocks` over one source, opened together and read one after
  another, give calls of `pass` 1 and `numPasses` 1 for each of the twelve.
- For each of the fifteen consumers, a function that throws on its first
  call: the consumer throws that same value, checked with `===` and not by
  its message, and the same `Variants` then gives its variants through
  `iterBlocks`, the 475 of `many.vcf` that passed a filter of the file, or
  the 500 it holds when it was opened with `onlyPassed` false.
- For each of the fifteen consumers, a function that calls the `free()` of
  the variants the run is reading: the `free()` is refused with popnei's own
  `Error`, that error leaves the function as any other thrown value does and
  the consumer gives it back, and the same `Variants` then gives its
  variants. The consumer has to hold its source over the call it hands it to
  the core in for that to happen, and two of the twelve did not until 25
  September 2026: `calcPopDiversity` and `calcLdAndDistPerPop` let the
  `free()` through, and what a user got was wasm-bindgen's own sentence
  about a value it holds, after the generated `free` had zeroed the pointer
  of the handle and dropped its entry of the `FinalizationRegistry`, so the
  source stayed in the memory of wasm with no handle left to free it. The
  list of the fifteen the two loops run over is held to the fifteen the
  binding crate names, which it writes in the message of a name that is of
  no consumer.
- A function that throws on the first call of the second pass: the PCA
  throws it, and what it threw is not popnei's error for a source that ended
  early.
- A function that frees the `Variants` while an iteration of `iterBlocks`
  reads them, and one that runs `calcPerIndividualStats` over the same
  `Variants`: neither traps, and the pass that was reading gives its
  variants.
- A function that throws only in the calls that end the run: the consumer
  returns its result, and the value is not thrown.
- A function that stops a pass is not called again for that pass, the call
  that would end it among them.
- A function that throws over a gzipped VCF and over a vars file: the
  consumer throws its value, which is what says that the failed read is not
  swallowed by the decompressor or turned into the error of a file that was
  cut short.
- A function that throws at a call inside the file while `calcPopDiversity`
  reads it: the diversity throws its value. The check that every consumer
  makes the passes `numPassesOf` says is not what holds a consumer to its
  run: a consumer that opened no run of its own, or dropped it before its
  pass read on, was still told of the first read of that pass, which is one
  call of pass 1 of 1, and passed that check while it gave popnei's error
  for the failed read in place of the value the application threw.
- For each of the fifteen consumers, the largest `pass` of the calls of one
  run equals `numPassesOf` of it with the same options, and the last call of
  each of its passes says the 117346 bytes of `many.vcf`. The association
  study is run twice here, once with the GRAMMAR-Gamma approximation and
  once without it, which are its two numbers of passes; the second pass of
  the run that asks for it gives one block and is told of the whole file all
  the same, the file being smaller than one range. The bytes of that last
  call are what says the end-of-pass call was made, and no other test reads
  them: a reader that outlives its run loses that call in silence, and a
  `std::mem::forget` on the chain of the diversity took the calls of its run
  from 0 and 117346 bytes to 0 alone with every other test of the package
  passing, which is a page's bar for a finished run standing at 0 per cent.

## How many passes a consumer makes

### What it gives

`numPassesOf` says how many times a consumer will read the file, before it
is started, so that a page can draw one bar for a whole run instead of one
per pass. Every consumer of popnei reads the file once, except two.

The PCA of the variants reads it twice when it is asked for the weights of
the variants: the weights need the eigenvectors, which are known when the
first pass ends. `doPcaFromVariants` with `numPrinComps` 0 reads it once.

The association study reads it twice when it is asked for the GRAMMAR-Gamma
approximation, which stands in for the denominator of the test of a mixed
model: the factor of that approximation is estimated from the first block of
a second pass over the same variants, as `docs/specs/gwas.md` has it.
`calcGwas` with `useGrammarGammaApprox` false, which is the default, reads
the file once.

```ts
numPassesOf(consumer: ConsumerName, options?: object): number
```

`ConsumerName` is the name of the function of this package that makes the
passes: `"calcPerVarDistribs"`, `"calcPerIndividualStats"`,
`"calcPairwiseKosmanDists"`, `"calcPopDists"`, `"calcPopDiversity"`,
`"calcRogersHuffR2Matrix"`, `"calcLdAndDistPerPop"`, `"calcKinship"`,
`"doPcaFromVariants"`, `"calcGwas"`, `"writeVars"`, `"iterBlocks"`,
`"calcVarDensity"`, `"doPcoaFromVariants"`, `"writeVcf"` and, since 7
October 2026, `"calcVariantsSummary"`, sixteen.
`options` is the options object that function takes, and only
`numPrinComps` of `doPcaFromVariants` and `useGrammarGammaApprox` of
`calcGwas` change the answer; each is checked as that function checks it, so
a `numPrinComps` that is not a whole number of 0 or more and a
`useGrammarGammaApprox` that is not a boolean are an `Error` here too. A
name that is of no consumer of the package is an `Error`. Like every other
function of the package it reads the defaults of those two from the core, so
it throws until `init` has been awaited.

The number the `Progress` of each call carries is this same number, from the
same function of the binding crate: the consumer asks it what its run makes
and opens the run with it, and the source puts it in every call. The two
cannot disagree, and the test of the item above is what holds them together.

A filter, the steps of a `Variants`, adds no pass: a filter is a reader over
the reader of its pass, as section 1 of the architecture says.

### How it is verified

Under node: `numPassesOf("doPcaFromVariants", { numPrinComps: 10 })` is 2,
with `numPrinComps` 0 it is 1, and with no options it is 2, which is the
default of 10 components; `numPassesOf("calcGwas", { useGrammarGammaApprox:
true })` is 2, with it false it is 1, and with no options it is 1, which is
the default of the exact denominator; each of the fifteen names other than
`doPcaFromVariants` gives 1 when it is asked with no options, `calcGwas`
among them; a name that is of no consumer throws, and so does a
`numPrinComps` of -1 and a `useGrammarGammaApprox` that is not a boolean.
The test of the item above
runs each of the consumers and compares the passes the calls showed with the
number this function gives, which is what would catch a consumer that grew a
pass and did not say so.

## The result so far

### What it gives

A pass over a VCF of millions of variants takes minutes, and the three
calculations whose results only add up over the blocks of a pass can give,
while their pass runs, the result over the variants read so far: the
distributions of `calcPerVarDistribs`, the rates of `calcPerIndividualStats`,
the counts of `calcVarDensity`, and the three of `calcVariantsSummary`
below. A page draws the histograms as they fill. The page cannot work these
out itself from the blocks of `iterBlocks`: the numbers popnei_web shows are
popnei's, which popnei verifies, and a second way of computing them in the
page would be one more thing to verify.

The result so far is pushed and not asked for. While a pass runs the worker
is inside wasm and reads no message, as "What the source tells the page"
says, so a request that the page posted would wait until the pass ended. The
consumer calls a function the application gives, between two blocks, and
that function posts the result to the page, which always holds the last
result it was sent. So a page whose user stops the pass loses nothing when
it ends its worker, which is how a page stops a pass when it has no reason
of its own for the function of `onProgress` to throw. The other way, the
page setting a flag in a `SharedArrayBuffer` that the pass would read,
needs headers that GitHub Pages does not send, as the opening of this spec
says.

The function is called after a block, the last one too, when `soFarEvery`
seconds have gone by since the pass started or since the last call. It is
given the result over the variants of the blocks read so far, of the type
the consumer returns, with its `passStats` as they stand after that block:
`passStats.numVars` is how many variants it covers, for the page to write
"over the first 1,200,000 variants". That result is the one the consumer
would return over those variants alone, because each of the three builds
its result at the end from totals it adds up block by block, and the result
so far is that same building, done on the totals as they are. Every
histogram of `calcPerVarDistribs` has its range before the pass starts,
`(0, 1)` by default, so the bins of a result so far are those of the last
one; the density with no lengths has the windows up to the last variant
read so far. A call after the last block gives what the consumer then
returns, which the page receives twice; the rule is kept simple, and a
pass needs to read no block ahead to know which block is its last.

Building the result costs what building the final one costs, at each call:
a few hundred numbers for the histograms, two for each individual, and four
for each window of the density. A density of 10 million windows, the most
it takes, copies them at each call, and a page that asks for one that wide
sets a longer `soFarEvery`.

When the function throws, the pass ends there and the consumer throws what
it threw, as when the function of `onProgress` throws: the binding crate
keeps the value, the pass ends as one stopped by `onProgress` ends, telling
the page nothing more of how far it read, and what a page can call from
inside the function, `free()` of the variants being read, another consumer,
the counts, gets the answers it gets from inside the function of
`onProgress`. An error of the pass, a variant past the length of its
chromosome among them, is thrown as it is today, and so is an error in
building the result so far, which is the error the page gets and not one
that says the pass was stopped. The function is called and not awaited: a
function that returns a promise is called, and what its promise rejects
with does not stop the pass. The time is read with `performance.now()` in
the binding crate, a clock that never goes back, and the core crate reads
no clock.

The owner decided on 7 October 2026 that the function is a function the
consumer calls, that every one of the four takes it, and that it is in the
TypeScript API alone: Python has no worker and no page, and the owner does
not want parameters in the Python API that a Python user has no use for.
The options not taken were a function only in `calcVariantsSummary`, and a
page that asks for the result through shared memory.

### Its TypeScript function

```ts
calcPerVarDistribs(variants, {..., onSoFar, soFarEvery})
calcPerIndividualStats(variants, {onSoFar, soFarEvery})
calcVarDensity(variants, windowSize, {chromLengths, onSoFar, soFarEvery})
calcVariantsSummary(variants, {perVar, perIndividual, density, onSoFar, soFarEvery})
```

`onSoFar` is a function that takes the result so far, of the type the
consumer returns. `soFarEvery` is a number of seconds, finite and 0 or more,
2 when it is not given; 0 calls the function after every block, which is
what the tests use. A `soFarEvery` given without `onSoFar` is an `Error`,
because it would be a value that does nothing, and so is an `onSoFar` that
is not a function and a `soFarEvery` that is not such a number; each is an
`Error` at the call, before the pass starts. `calcPerIndividualStats` takes
an options object for these two alone, which it did not take before; a call
with none works as it did.

This is a difference between the two APIs that goal 2 of
`docs/objectives.md` asks to be written down: the Python functions take
neither.

### How it is verified

Under node, on the vars file of `many.vcf` opened with every variant and
written with `writeVars` in batches of 100 variants, so that each pass is
five blocks of 100; `many.vcf` itself is one block for the VCF reader,
which cuts its blocks by the genotypes they hold, and would test nothing.
With `soFarEvery` 0, `calcPerVarDistribs` calls the function five times,
with a `passStats.numVars` of 100, 200, 300, 400 and 500; the result of each
call equals, field by field but for `passStats`, that of
`calcPerVarDistribs` over the same `Variants` with `filterFirstN` of that
number of variants, which reads the same first variants and whose
`passStats` has a filter more, so of `passStats` only `numVars` is
compared; the last call equals what the consumer returns. The same for
`calcPerIndividualStats`, and for `calcVarDensity` with `chromLengths`
`{chr1: 30000, chr2: 25000}`, whose windows are all there from the first
call, and with `chromLengths: {}`, which takes none and has the windows
grow; `many.vcf` has no `##contig` length, so its own lengths would give
the second case again.
With a `soFarEvery` of 3600 the function is never called over that file,
and with none, the default of 2 seconds, it is never called either. Over
the same `Variants` with `filterByMaf(0.95)`, each call carries the counts
of that filter after its block, 100 to 500 given, as literals. Over the
vars file in 500 batches of one variant, with a `soFarEvery` of 0.0001 s,
the function is called no more times than the time the pass took, measured
around the call plus 1 ms, divided by that interval, plus one, which is what
counting the interval from the last call gives on any machine; a clock that
counted from the start of the pass would call it after nearly every one of
the 500 blocks. A density with `chromLengths: {chr1: 6000}` is called once, at 100
variants, and then throws the error of a variant past the length, and a
result so far that cannot be built throws its own error.
The tests of the item above for `onProgress` that stop a pass, in
`test/stop.test.ts`, are made again for `onSoFar`: a value thrown at the
second call is what the consumer throws, the next call over the same
`Variants` runs whole, and `free()` from inside the function is refused as
it is from `onProgress`. A `soFarEvery` with no `onSoFar`, an `onSoFar` of 3
and a `soFarEvery` of -1 and of `NaN` are each an `Error` at the call.

## The three statistics of a file in one pass

### What it gives

`calcVariantsSummary` gives, in one pass, what `calcPerVarDistribs`,
`calcPerIndividualStats` and `calcVarDensity` give in three: the
distributions of the statistics of each variant, the rates of each
individual and the density of the variants along the chromosomes. A page
that shows the three when a file is opened reads the file once. Each of the
three is what its own consumer gives for the same options, to the bit: the
three are added up from the same blocks by the same code, and only the pass
is shared. An error of any of the three ends the pass, and none of the
three is returned; a page that wants the others when the density refuses a
variant past the length of its chromosome calls them on their own.

What it saves depends on the file. Over `big.vcf`, the plain VCF of 403
MB, 100000 variants of 1000 diploid individuals that
`crates/popnei/benches/make_big_vcf.py` writes, the three passes took
2.240 s and the one 1.247 s, 44% less; over `big.vars`, the vars file
popnei writes of it, 0.400 s and 0.391 s, 2% less. Each is the median of 7
runs, the three and the one in turn, each timed from opening the file to
the result; the 7 runs of the three over the VCF spread from 2.176 to
2.300 s, those of the one from 1.202 to 1.289 s, and over the vars file
from 0.383 to 0.424 s and from 0.371 to 0.402 s. The distributions were
of the six statistics `calcPerVarDistribs` gives when it is asked for
none, the density in windows of 100000 base pairs, and the counting and
the parsing of the VCF ran on one thread, rayon's pool being built with
one. It was measured natively, with the core crate and not in wasm, on 7
October 2026 on the owner's Apple M5 Pro, with other work running on the
machine, a load average of 15 for its 18 cores, by

```text
cargo bench --bench variants_summary -- \
    /Users/jose/devel/popnei-bench/big.vcf \
    /Users/jose/devel/popnei-bench/big.vars --runs 7
```

Two sets of 7 runs of the same command just before gave 43% and 44% on the
VCF and 3% and 2% on the vars file.

The counting of each statistic is not shared, only the read. The
distributions and the rates each work 0.19 s on the blocks in their own
pass, and 0.38 s together in the one pass. Over the VCF each of those two
passes also waits about 0.91 s for its blocks, the parsing of the
genotypes, and the density, which reads no genotype, takes 0.07 s: about
0.91 + 0.19 + 0.91 + 0.19 + 0.07 s for the three. The one pass parses the
genotypes once, so it saves one parsing and the pass of the density. Over
the vars file a pass waits under 0.01 s for its blocks, so there is almost
no read to share, and a page that opens a vars file gains from
`calcVariantsSummary` little more than calling one function instead of
three. `crates/popnei/benches/variants_summary.rs` has these clocks of
each pass.

The owner decided on 7 October 2026 that it is one consumer of exactly these
three, and that it is in the TypeScript API alone. The options not taken
were a way to run any set of consumers in one pass, which would split every
one of the fifteen into what it adds up and what it gives at the end and
cannot hold the PCA, which reads the file twice; and a Python function of
the same.

### Its TypeScript function

```ts
calcVariantsSummary(variants, {
  perVar,         // the options of calcPerVarDistribs, or absent
  perIndividual,  // {} for the rates of each individual, or absent
  density,        // {windowSize, chromLengths}, or absent
  onSoFar, soFarEvery,
}): VariantsSummary
```

A statistic is given when its key is there with a value that is not
`undefined`, and left out otherwise; `perVar` is checked as the options of
`calcPerVarDistribs` are, but for `onSoFar` and `soFarEvery`, which are
options of the summary and are refused inside it; `density` as the
`windowSize` and `chromLengths` of `calcVarDensity` are, with the errors
naming `density.windowSize`; and `perIndividual` is an empty object. A call with none of the three is an
`Error` that says to ask for one. `VariantsSummary` has `perVar`, a
`PerVarDistribs` or `null`, `perIndividual`, a `PerIndividualStats` or
`null`, `density`, a `VarDensity` or `null`, each of the type its consumer
returns so that the code of a page that draws one draws it from here
unchanged, and `passStats`, the counts of the one pass, which each of the
three also carries. The pass asks the reader for what the three ask for
together: the genotypes when `perVar` or `perIndividual` is there, and the
chromosome and the position when `density` is. A pass that gives no
variant is the `Error` the three give for it.
`numPassesOf("calcVariantsSummary")` is 1, and it reads none of the
options.

It mirrors nothing: pyNei has no such function, and Python has none, which
is a difference between the two APIs that goal 2 of `docs/objectives.md`
asks to be written down.

### How it is verified

Under node, on the vars file of `many.vcf` in batches of 100 of the item
above, and on `many.vcf`: `calcVariantsSummary` with the three gives
`perVar`, `perIndividual` and `density` equal, field by field and to the
bit, to the results of the three consumers with the same options, and a
`passStats` equal to theirs; with one of them left out, that one is `null`
and the others are as before; and with none, an `Error`. With `onSoFar` and
`soFarEvery` 0, the three of each call equal those of the consumers over
`filterFirstN` of that number of variants, as in the item above. The cargo
tests of `calc_variants_summary` make the same comparison on `many.vcf`,
and check with a reader built for the test that records what it is asked
for that the density alone asks for the chromosome and the position and
no genotypes.

## The file of a writer in pieces

### What it gives

`writeVars` and `writeVcf` build the whole file in the memory of wasm and
then copy it into the array they return, so while the call returns the tab
holds the file twice, and the memory of wasm, which never shrinks, keeps the
room of the file until the worker ends. popnei_web offers the variants that a
user's filters keep as a download, written from its worker. Its largest case
is `big.vars`, the vars file of `big.vcf` of "The three statistics of a file
in one pass", 100000 variants of 1000 diploid individuals: written as a plain
VCF of 403572916 bytes, the write grew the memory of wasm by 441974784 bytes,
and the array it returned held the 403572916 bytes again, 846 MB for one
download beside the 83 MB of the open vars file. That was measured under
node 24 on the owner's Apple M5 Pro on 7 October 2026, with
`memory.buffer.byteLength` of the wasm module before and after the call. A
tab runs out of memory somewhere between 2 and 4 GB, by the browser and the
machine, which issue 13 says and nobody has measured here.

With the option `onBytes`, the writer gives the file to a function of the
application while its pass runs, one piece at a time, and keeps none of it in
wasm: a piece is copied out as soon as it is full. What the memory of wasm
holds then is what the pass holds for one block, the block the reader gave
and the text or the batch the writer made of it, and one piece. The same
write of `big.vars` with `onBytes` grew the memory of wasm by 59703296
bytes, in 385 pieces. Two runs of each took 3.4 and 5.5 s whole and 3.1
and 6.0 s in pieces, which does not tell the two apart. Where the pieces go is the
application's choice. A worker can write each one to the origin private file
system, the file system the browser gives each site, through its
`FileSystemSyncAccessHandle`, which writes synchronously as the function has
to, and then the tab holds no copy of the file; or it can keep them in an
array and build a `Blob` of them, which the browser may keep out of the
memory of the page, and then the heap of JavaScript holds the file once
until the `Blob` is built.

A piece is 1048576 bytes, 1 MiB, the size the writers already cut the file
into in wasm, and popnei chooses it as it chooses the size of a range
(**Open 3**, below). A piece is given as soon as it is full, which can be in
the middle of a block, and the last one when the writer has finished the
file, so every piece but the last holds 1 MiB and the last from 1 byte to 1
MiB; no piece is empty, and every file has at least its header, so the
function is called at least once. The pieces, in the order they are given,
are the bytes the same call gives without `onBytes`. Each piece is a new
`Uint8Array` of the application's, which it can keep: a view into the
memory of wasm would stop being valid the next time that memory grows.

When the function throws, the pass ends and the writer throws the value it
threw, as when the function of `onProgress` or of `onSoFar` throws, which
is the "meanwhile" of **Open 2** below. Every piece is given while the run
is open, the last one too, so a `free()` of the variants from inside the
function is refused, as it is from inside those two, because a run is
reading them.

A function that returns a promise, an `async` function among them, ends
the pass at its first piece, and the writer throws an `Error` that says
`onBytes` has to keep each piece before it returns. The writer cannot wait
for the promise, since the pass runs inside wasm until the file is done,
and a promise that failed later, a write the private file system refused
for want of room, would leave the application a file with pieces missing
and a call that said it was written. `onSoFar` is called and not awaited,
because what a promise of it loses is one result so far. A piece that the
heap of JavaScript has no room to copy ends the pass with the error the
browser gave, a `RangeError`, and the variants can be freed and read again,
as after any other stop. Nothing counts the bytes given, so a file has no
bound of size but the time it takes to write.

When the pass fails after some pieces were given, a wrong line of a VCF
among the causes, the writer throws that error, and the pieces given are
the start of a file that has no end: the application drops them, as the
Python writers remove the file of a call that failed. Without `onBytes` the
two writers are as they were.

The owner decided on 7 October 2026 that the pieces go to a function given
as an option, which the writer calls. The options not taken were an
iterator of pieces, which would need the loop of the core's writers to stop
between two pieces and go on when asked, where today it runs to the end of
the file in one call; and a writer that always gives pieces and no `bytes`,
which changes what every caller of the two writers gets.

### Its TypeScript function

```ts
export type OnBytes = (piece: Uint8Array) => void;

export interface WrittenInPieces {
  passStats: PassStats;
}

writeVars(variants, options: WriteVarsOptions & {onBytes: OnBytes}): WrittenInPieces
writeVars(variants, options?: WriteVarsOptions & {onBytes?: undefined}): VarsWritten
writeVars(variants, options?: WriteVarsOptions): VarsWritten | WrittenInPieces
```

`writeVcf` has the same three, with `WriteVcfOptions` and `VcfWritten`.
`WriteVarsOptions` and `WriteVcfOptions` each gain `onBytes?: OnBytes`.
With `onBytes` the result has the counts of the pass alone, `passStats`, as
every consumer gives them; without it, it is `VarsWritten` or `VcfWritten`
as it was, the bytes and the same counts. The third signature is for a
caller whose `onBytes` may or may not be there. An `onBytes` that is not a
function is an `Error` at the call, before the pass starts, and an
`onBytes` of `undefined` is one not given.

This is a difference between the two APIs that goal 2 of
`docs/objectives.md` asks to be written down. The Python writers take a path
and write each piece to the file as it comes, so a Python user has no use
for the option, which is the reason the owner gave for keeping `onSoFar` out
of Python.

### How it runs

The function is kept in the entry of the run in `RUNS`, and the sink the
core writes into holds the number of that run, as a pass holds the number
of its source: `write_vcf` of the core asks its sink to be `Send`, and no
handle of JavaScript is. The sink fills one buffer of 1 MiB, reserved once,
and when it is full copies it into a new `Uint8Array`, empties it and calls
the function, cloned out of the table so that no table is borrowed while it
runs. The array is built through `Reflect.construct`, which gives back
what its constructor throws, where the constructor js-sys binds would let
the throw unwind the stack of wasm past the code that closes the run. What
the function returns is read: an object with a `then` that is a function
stops the run with the `Error` above, as a throw would. It gives nothing at
a `flush`. The last piece is given after the
core's writer returns and before the run closes. A throw is recorded in the
run as a throw of `onSoFar` is, and the write that called the function
fails with an `io::Error` of the kind `Other`, which no user reads, because
the consumer throws the recorded value in its place; `Interrupted` would be
written again by `write_all`. A run with no entry in `RUNS`, which no tab
reaches, would leave the file nowhere, so the writer refuses it with an
`Error`.

### How it is verified

Under node, in `test/write_in_pieces.test.ts`. The pieces joined are, byte
for byte, the bytes of the same call without `onBytes`, and the result is
`{passStats}` with the counts of that call. Each of `cases.vcf` and
`many.vcf` is written with `writeVars` and `numVarsPerBlock` 3, with
`writeVcf` bgzipped and with it plain, and the vars file of `many.vcf` with
`writeVcf` plain; each of those files is under 1 MiB, so each is one piece.
The VCF of 14000 variants of 600 individuals of `vars_memory.test.ts`, 34.0
MB of text, is written the same three ways in several pieces, every one but
the last of 1048576 bytes; and a plain VCF whose header is padded to make
it 2097152 bytes is two pieces of 1048576 bytes, with no empty one after
them. The pieces of one write are compared again after a second write of
the whole file, and have not changed. For each of the three writes, a value thrown at the second piece is
what the writer throws, after two calls, and the next call over the same
`Variants` gives the whole file; and a `free()` from inside the function is
refused with the error of a run that is reading. That VCF of 14000 variants
with a line after them whose genotype is not a number gives at least one
piece, written plain, and then throws the error of line 14003. An `onBytes`
of 3 is an `Error` at the call that names `onBytes`, before the page is told
of any read, and an `onBytes` of `undefined` gives the bytes. An `async`
`onBytes` makes the writer throw the `Error` that names `onBytes` and a
promise, after one call. A `Uint8Array` that throws a `RangeError` when it
is built at the length of a piece, from the third piece on, makes the
writer throw that `RangeError`, after which the variants are freed with no
error. From inside `onBytes`, `writeVars` and `writeVcf` over the same
`Variants`, with and without `onBytes` of their own, give their files, and
the outer write gives its whole file. The test file asserts with
`@ts-expect-error` that the result of a call with `onBytes` has no `bytes`,
and that a call without it is a `VarsWritten`.

The memory is measured in `test/write_in_pieces_memory.test.ts`, alone in
its file, over a VCF of 30000 variants of 600 individuals, which is large
enough that a block is well under half the file. Measured when it was
built: the vars file of 26236330 bytes, written in pieces in batches of
100, grew the memory of wasm by 2293760 bytes; the plain VCF of 72923053
bytes by 18415616 in pieces, and by 73465856 written whole after them. The
test asserts that each write in pieces grows it by less than half its file
and the whole write by more than half, a bound that holds for any block a
reader gives at that size and fails for a writer that keeps the file.

## The Rust interface

In the binding crate. Nothing of this is in the core crate, and nothing of
it is public to a user of the core. `VcfSource` and `VarsSource`, which
implement `OpenSource` and today hold the bytes of their file, and
`JsPopneiError`, the error every function of the crate gives, are the types
that `crates/popnei-js/src/` already has.

```rust
/// One pass over the bytes of a source, which the readers of the core take
/// as they take a cursor over an array of bytes. It counts what the pass
/// has read and tells the page, and it is `Send`, because what it holds of
/// JavaScript is the number of an entry of `IN_JAVASCRIPT` and not a
/// handle.
pub(crate) struct PassOverTheBytes {
    bytes: TheBytes,
    /// Which run of `RUNS` this pass belongs to, which is what says which
    /// source it reads, which pass of the run it is and how many there are.
    run: u32,
    /// Which pass of the run this is, 1 for the first, and 0 until the
    /// first read takes the next number from the run: the PCA builds both
    /// of its readers before either of them reads.
    pass: u32,
    /// How many bytes this pass has read, and how many it had read when the
    /// page was last told. The second is what the size of a range is
    /// compared against to decide whether to tell it again.
    bytes_read: u64,
    told_at: u64,
    num_bytes: u64,
}

/// Where the bytes of a pass come from.
enum TheBytes {
    /// A copy of the whole file in the memory of wasm, which every pass
    /// over that source shares: the `SharedBytes` of `source.rs`, an
    /// `Arc` of the bytes wasm-bindgen filled.
    InMemory(Cursor<SharedBytes>),
    /// A file of the page, read one range at a time: which entry of
    /// `IN_JAVASCRIPT` holds it, the range that is held, where that range
    /// starts, where the pass reads next, and how many bytes the file
    /// holds, which a read and a seek from the end are against.
    OfAFile {
        source: u32,
        range: Vec<u8>,
        range_at: u64,
        pos: u64,
        num_bytes: u64,
    },
}

impl Read for PassOverTheBytes {}
impl BufRead for PassOverTheBytes {}
impl Seek for PassOverTheBytes {}

/// What a source keeps in JavaScript, which nothing of Rust may hold and
/// stay `Send`. One entry per source, which `VcfSource` and `VarsSource`
/// hold the number of, and `free()` of the source takes out once no run
/// over it is open.
struct InJavaScript {
    /// The file the ranges are read from and the reader of them, and
    /// nothing for a source whose bytes are already in the memory of wasm.
    file: Option<(Blob, FileReaderSync)>,
    /// What the application is told the progress with, the function of
    /// `Variants.onProgress`.
    told: Option<js_sys::Function>,
    /// Whether the source was freed while a run over it was open, which is
    /// when its entry goes.
    freed: bool,
}

/// One run of one consumer: which source it reads, how many passes it
/// makes, how many have begun, and what the function threw. It is taken out
/// when the consumer returns, and for `iterBlocks` when the iteration ends
/// or the pass is freed.
struct Run {
    source: u32,
    num_passes: u32,
    passes_begun: u32,
    /// What each pass that has ended read, in the order of their numbers,
    /// which the end of the run tells the page: a pass writes it here when
    /// it is dropped, and the run is what is left to tell it.
    passes_that_ended: Vec<(u32, u64)>,
    /// What `told` threw, which ended a pass of this run and is what the
    /// consumer throws in place of the error the core gave. It is nothing
    /// when the run starts, so no run throws what another one was stopped
    /// with.
    stopped_with: Option<JsValue>,
}

thread_local! {
    static IN_JAVASCRIPT: RefCell<Vec<Option<InJavaScript>>>;
    static RUNS: RefCell<Vec<Option<Run>>>;
}

/// A consumer of the package, with the argument of the one whose number of
/// passes depends on it.
pub(crate) enum Consumer {
    PerVarDistribs, PerIndividualStats, KosmanDists, PopDists, PopDiversity,
    R2Matrix, LdAndDist, Kinship, PcaOfVariants { num_prin_comps: usize },
    Gwas { use_grammar_gamma_approx: bool }, WriteVars, IterBlocks,
}

impl Consumer {
    /// How many passes over the source this consumer makes. It is what
    /// `numPassesOf` gives and what every `Progress` of the run carries.
    pub(crate) fn num_passes(&self) -> u32;
}

/// A file of variants that was opened, which every pass reads again.
pub(crate) trait OpenSource {
    fn ploidy(&self) -> usize;
    /// The run of `consumer` over this source, which every reader the
    /// consumer opens belongs to and which is taken out when it is dropped.
    fn starts_a_run(&self, consumer: &Consumer) -> RunOfAConsumer;
    fn reader(
        &self,
        run: &RunOfAConsumer,
        num_vars_per_block: Option<usize>,
    ) -> Result<Box<dyn BlockReader>, popnei::Error>;
}

/// The run a consumer holds: the number of its entry of `RUNS`, which its
/// readers carry, and which is taken out of that table when the consumer is
/// done with it.
pub(crate) struct RunOfAConsumer(u32);

impl Drop for RunOfAConsumer {}
```

The read that a thrown value ends fails with `std::io::Error::other`. The
kind matters: `ErrorKind::Interrupted` is read again by four loops of the
core, `read_line_of` of `io::vcf`, `take_from` and `take_from_into` of
`io::bgzf` and the `read_exact` of `bytes_at` of `io::vars`, so a stop
written with it would never end the pass; and `ErrorKind::UnexpectedEof` is
what `bytes_at` turns into the error of a vars file that was cut short, so
a stop written with it would reach the user as a damaged file.

What the consumer throws is not the error the core gives back for that read:
`JsPopneiError` has a case that carries a `JsValue` through untouched, and
`impl From<JsPopneiError> for JsValue` gives it back as it is, where every
other case becomes an `Error` with its message. The binding puts the error
of a consumer into that case when the run it opened was stopped, whatever
error the core gave, so nothing depends on which reader turned the failed
read into which error.

## Speed

A range is 4194304 bytes, 4 MiB, and no argument of the API carries that
number (**Open 3**, below). One pass over a file of the page leaves the
memory of the wasm module at 14155776 bytes at that size, and the bound the
browser test asserts is 16777216 bytes, 16 MiB.

Both come from the measurement of 24 September 2026, in Chromium on the
owner's Apple M5 Pro, over a VCF of 299994147 bytes and 1285000 variants of
50 diploid individuals built from the body of `many.vcf` repeated, with
`calcPerIndividualStats` as the consumer, the best of five runs of each
point. `docs/reports/js-sources-measurement.md` has the tables, the machine,
the file and the script that took the times.

| the size of a range | one pass | the memory of wasm after it |
|---|---|---|
| 256 KiB | 1069 ms | 3407872 |
| 1 MiB | 1006 ms | 4718592 |
| 4 MiB | 996 ms | 14155776 |
| 16 MiB | 958 ms | 51904512 |
| the whole file as an array of bytes | 921 ms | 302383104 |

The memory of wasm is the whole module when the pass has ended, what
`memory.buffer.byteLength` gives: the module itself, 1310720 bytes before
any pass, and every block a pass allocated, of which nothing is given back
to the browser. From 1 MiB up it is three ranges and 1572864 bytes, to the
byte, at each of the 1, 2, 4, 6, 8 and 16 MiB it was read at; at 256 KiB it
is 3407872, which is 1048576 more than that count gives. So what a pass
holds is not a multiple of its range: the fixed part makes the memory 13.0
times a range of 256 KiB, 4.50 times one of 1 MiB, 3.38 times one of 4 MiB
and 3.09 times one of 16 MiB.

Reading by ranges costs 8.1 % of a pass over that file against holding it
whole, 75 ms of about a second, and leaves 21.4 times less of the memory of
wasm. The three sizes of 1, 4 and 16 MiB are within 48 ms of one another on
a pass of about a second, and the five runs of one point spread by as much
as 73 %, so the times do not choose between them. What chooses 4 MiB is the
calls into the browser and the memory: a pass over that file makes 287 calls
at 1 MiB, 72 at 4 MiB and 18 at 16 MiB; every byte of the measurement came
out of a `Blob` built in the memory of the browser, so what a call costs
when the file is on a disc was not measured, and a range four times smaller
pays it four times as often. At the other end, 16 MiB leaves 51904512 bytes
of the memory of wasm for one pass where 4 MiB leaves 14155776.

The size of a range is measured again when the reading of a range changes,
and the bound of the test with it.

## Open points

The owner decides these three. Until then the implementer follows the
"meanwhile" of each.

**Open 1: where the application sets the function that is told the
progress.** The options are a method of `Variants`, `onProgress`, which an
application changes before each run; an option of `openVcf` and `openVars`,
fixed for the life of the source; or an argument of every consumer. The
method is the one that fits an application that runs one analysis after
another over one source and tags each with a key of its own, since the
function can carry that key, and it leaves the options of every consumer as
the ones of the Python API, which goal 3 of `docs/objectives.md` asks the
TypeScript API to mirror. An argument of every consumer is the same power
with thirteen places to add it to and thirteen more lines of documentation.
Recommendation: the method of `Variants`. Meanwhile the implementer writes
that.

**Open 2: what a run that was stopped throws.** The options are the value
the function threw, which an application recognises with `===` and which
carries whatever it put in it; or an `Error` of popnei whose message holds
the message of that value, which keeps every failure of popnei of one kind
and makes an application read a message to tell its own cancel from a file
that could not be read. Recommendation: the value it threw, because telling
a cancel apart is what the application does with it, and popnei's errors
stay the ones it made itself. Meanwhile the implementer writes that.

**Open 3: whether an application chooses the size of a range.** The options
are that popnei chooses it, from the measurement of "Speed", and nothing of
the API says it exists; or `openVcf(source, { numBytesPerRange })`, which
lets an application trade a slower pass for a bar that moves more often, and
which is one more argument to check and to document. Recommendation: popnei
chooses it, and the option is added when an application asks for it, since a
range that is too small costs time and no result. Meanwhile the implementer
writes popnei's number.

## Not in this spec

- Firefox and WebKit in the tests, which the owner left for when an
  application needs them.
- The fingerprint of a variant file, the reader of CSV and TSV, and the
  inference of the types of the columns, which the applications also ask of
  popnei: each is its own request, and none of them is about how a file is
  read.
- Writing to the private file system of the browser. A writer gives its
  file in pieces, as "The file of a writer in pieces" says, and where they
  go is the application's.
- Reading a file by ranges in Python. Python opens a file by its path, and
  under pyodide it reads the filesystem that emscripten gives it.
- The wheel for pyodide, which is built for the other wasm target and has no
  wasm-bindgen: nothing of this spec is in it.
