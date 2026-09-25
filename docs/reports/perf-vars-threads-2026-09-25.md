# The vars file reader on the threads, and what it was bounding

25 September 2026. Reading popnei's own file of variants was, until this
review, the one pass of popnei that ran on a single thread whatever the
machine, and `docs/reports/perf-gwas-2026-09-24.md` had just found it
bounding an association study: the study of 100000 variants took 0.124 s
over the 77.8 MB vars file of a panel and 0.111 s over the 403 MB VCF of the
same genotypes, because the VCF reader parses its lines on the threads of
rayon and the vars reader decompressed lz4 on one. The owner asked for the
reader to use the threads.

It does now, in three commits on the branch `perf/vars-threads`, and
**nothing is merged into `main`**. A pass over that panel with the genotypes
alone asked for takes 27.23 ms where it took 112.21 ms, and every block it
gives is identical to the block it gives on one thread. The study over the
vars file now takes 0.112 s against the VCF's 0.111 s, the missing data
filter over the same panel 0.040 s where it took 0.128 s, and the Kosman
distance of every pair of individuals 0.128 s where it took 0.217 s. What is
left for the owner is at the end of section 2, and the largest of those
decisions is the memory: the reader now holds up to eight batches at once
where `docs/specs/io_vars.md` says a pass holds one block.

## The words this document uses

A **vars file** is popnei's own file of variants, one arrow IPC file of six
columns whose buffers are compressed with **lz4**, which `docs/specs/io_vars.md`
specifies. Its variants are written in **batches**, and the reader gives each
batch as a **block**, the run of variants that every calculation of popnei
works over: 5000 variants of 1000 individuals in every measurement here, so
the panel of 100000 variants is 20 batches. Twelve consumers of popnei read
blocks.

A **buffer** is one column's bytes inside one batch, compressed on its own,
which is how the arrow format works: a batch of this panel decompresses
11384500 bytes in 14 buffers, of which the genotypes are 10000000 bytes in
one. An lz4 **frame** is what one compressed buffer holds, and a frame is
made of **blocks of lz4**, up to 4 MiB of output each; they are
**independent** here, meaning each one can be decompressed without the ones
before it, which is what lets one buffer be split at all.

**rayon** is the library that runs work on several threads, and its **pool**
is the set of threads it runs on; popnei never builds the global pool, so the
count comes from `RAYON_NUM_THREADS` or from a pool a caller installed. A
**sampling profile** stops a process a thousand times a second and records
where it is; the **self time** of a function is how many of those stops were
inside it and not inside something it called. **DHAT** counts every allocation
a program makes, exactly and the same on every run; it is the `dhat` crate
behind a cargo feature, and popnei does not have that feature yet.

Every reader of popnei hands its blocks over one at a time through one trait,
`BlockReader`, and a **chain** of them is a reader with filters over it, each
taking the blocks of the one before. The **read ahead reader**,
`with_one_block_ahead` of `crates/popnei/src/block.rs`, runs such a chain in a
thread of its own so that it builds the next block while the consumer works on
the one in hand; `docs/reports/perf-gwas-2026-09-24.md` built it and the
association study is still its only caller.

The **panel** of every number here is
`/Users/jose/devel/popnei-bench/bigcalled.vars`, 100000 variants x 1000
individuals with every genotype called, 77820074 bytes, 20 batches; the same
genotypes are in `bigcalled.vcf`, 403572954 bytes. `big20000.vars`, 20000
variants and 4 batches, is the smaller panel that the 21 ms of "Speed" of
`docs/specs/io_vars.md` is stated on. The **machine** is the owner's Apple M5
Pro, 18 cores, 64 GB, macOS 27, native `aarch64`. Every wall time is the best
of the runs stated, and the load average of the machine is given with each
set: the owner's own applications kept it between 3.6 and 8.7 through this
review, so every number that decides something was taken **interleaved**,
before and after alternating in one window, with a second binary built from
the commit before the change.

## 1. The scope and its limits

`crates/popnei/src/io/vars.rs`, the reader and the writer of a vars file, and
`crates/popnei/src/block.rs`, the `BlockReader` chain and the read ahead
reader over it, at `ae6a190`, the head of `main`. The twelve consumers and
the two binding crates were read as callers and not as code to change.

Nine subjects were reviewed, one reviewer each, all of them given the
benchmark output and the profile below: whether the measurements themselves can
decide anything, what a change would do to the values the reader gives, what it
allocates per batch, how its data sits in the caches of this machine, its
threads, its loops over the genotypes, what it asks of the file, what a browser
tab gets, and what crossing into Python costs. What each one covered is in the
findings of section 5 that carry its subject.

**The four ways the owner named, none of them chosen when the review
started.** They are referred to below by these letters, and nothing else in
this document is lettered:

- **(a)** decompress the buffers of one batch on the threads of rayon, which
  needs knowing whether arrow-rs will hand over its compressed buffers or
  whether popnei reads the batch itself;
- **(b)** decompress the next batches while the pass works, which the read
  ahead reader already does for one block and which could hold more;
- **(c)** write the file with no compression, which arrow-rs reads already and
  which makes the file 2.9 times larger, as an option and not as the default;
- **(d)** a codec that decompresses faster than lz4 and builds for both wasm
  targets without a C compiler, which is why zstd was refused in "The
  compression" of `docs/specs/io_vars.md`.

What was built is (b), in a form that decompresses as many batches at once as
the pool has threads rather than one batch ahead. (a) has a ceiling this review
measured, (d) does not exist inside the arrow format, and (c) is untouched and
is the one of the four that would help a browser tab; sections 6.1 and 6.3 and
the notes N1 and N2 of section 5 give the evidence for each.

**What was not measured.** No timing in wasm, where there are no threads and
the reader keeps the path it had; the note N2 of section 5 says what a browser
tab would need instead, which is (c) above.
No timing at 10000 individuals, the largest dataset of `docs/objectives.md`:
a batch holds about the same number of genotypes whatever the shape, so a
batch there is 250 variants and the count of batches grows tenfold, which
this review did not run. No allocation count: the `dhat` feature that
`docs/reports/perf-io-2026-09-22.md` designed in September still does not
exist, so every finding below that would be gated on a count is gated on a
wall time instead, and says so. The write was measured and not changed.

## 2. The verdict: apply, and three commits are on the branch

The reader decodes a window of batches on the threads of the pool it is
called in, the window being eight or the threads of that pool, whichever is
fewer. The bytes of each batch are still read from the source one after
another, and every block is still built on the calling thread in the order of
the file, because the chromosome table hands out its numbers in the order the
names are first seen and nine places of popnei index their results by a
running count of the variants: a block out of order is a wrong result and not
a slower one.

Best of 10 runs of `cargo bench --bench vars_file`, the pass with the
genotypes alone asked for, interleaved with a binary built from the commit
before the change:

| | 1 thread | 18 threads |
|---|---|---|
| the panel, before | 112.21 ms | 112.37 ms |
| the panel, after | 108.53 ms | **27.23 ms** |
| the 20000 variant panel, before | 22.00 ms | — |
| the 20000 variant panel, after | 22.01 ms | **7.43 ms** |

4.1 times on the panel and 3.0 times on the smaller one, whose four batches
cannot fill a window of eight.

**The 21 ms of "Speed" of the spec is stated on one thread, and this change
does not move that number: 22.00 ms before and 22.01 ms after.** Both are above
the 21 ms, and both are above the 20.17 to 20.32 ms that
`docs/reports/vars-file.md` recorded on 22 September at a load average of 1.15
to 1.57. Whether the target is met today is not this review's to say, because
this machine was never idle: what can be said without a quiet machine is that
the pass now does 27114003949 instructions where it did 26799260757 in
September, 1.2 per cent more, so about 0.25 ms of the 1.7 ms difference is more
work and the rest is the load. Twenty runs at one thread with nothing else on
the machine would settle it.

What the four passes that a consumer waits for now cost, each interleaved,
18 threads:

| | before | after |
|---|---|---|
| the association study, continuous trait | 0.124 s | **0.112 s** |
| the same study over the VCF of the same genotypes | 0.111 s | 0.111 s |
| the missing data filter, threshold 0.1 | 0.128 s | **0.040 s** |
| the Kosman distance of every pair | 0.217 s | **0.128 s** |
| the Kosman distance, 1 thread | 0.880 s | 0.883 s |
| the write of the panel | 202.58 ms | 203.27 ms |

So the study over the vars file is level with the study over the VCF, which
is what the owner's question was about, and a consumer whose own work per
block is small gains most: the filter's pass is 3.2 times faster.

**A regression that the reader's own benchmark could not see, and its fix.**
The first commit made the window eight batches whatever the pool. The Kosman
pass then lost 17 per cent on one thread, 0.880 to 1.026 s, where the
reader's own benchmark lost 1.1 per cent: that benchmark's consumer only adds
the genotypes up, and a consumer that reads them finds every block cold,
because eight batches, about 80 MB of output, are decoded before it touches
the first of them. This machine has 128 KB of first level data cache for each
performance core and one 16 MB second level cache shared by six of them. The
second commit makes the window the smaller of eight and
`rayon::current_num_threads()`, so a pool of one thread decodes one batch at a
time, which is what the reader did before, and one thread came back to 0.883
s against 0.880 s over three interleaved pairs.

### What is left for the owner to decide

1. **The memory, which is the price of this change.** The window holds about
   6 MB more of peak resident set for each batch in it on this panel, measured
   with `/usr/bin/time -l`: 499.8 MB before against 623.0 MB after at 18
   threads, where the benchmark's own write section sets the 500 MB by holding
   20 blocks at once. So the reader's own cost is up to eight times the 3.9 MB
   of compressed bytes and the 11.4 MB of decoded buffers of a batch, about 120
   MB at 18 threads on this panel, and the same at any shape of dataset,
   because a batch holds about the same number of genotypes whatever the
   individuals, 5 million: 250 variants of 10000 individuals there against
   5000 variants of 1000 here.
   `docs/specs/io_vars.md` says the memory of a pass is one block, and
   `docs/architecture.md` says two or three blocks are alive with the read
   ahead. Both sentences are now wrong for a native build, and what replaces
   them is the owner's to settle: the window as it is, a smaller constant, or a
   cap in bytes rather than in batches. Recommendation: keep eight, and state
   the bound in the spec as the window times the batch.
2. **`docs/specs/io_vars.md` owes an item.** Line 648 says of the reader "It
   uses no threads: reading one batch ahead in a thread of its own comes with
   the read ahead of `docs/specs/block.md`". That is now false. The item is
   owed either way: as the description of the window, if this is kept, or as
   the measurement that closed it, if it is not. Its "Speed" section keeps its
   21 ms on one thread and gains the 18 thread number above.
3. **Whether the window should follow the pool upwards as well as downwards.**
   At 18 threads a window of 18 read the panel in 24.2 ms against 27.8 ms at
   eight, 11 per cent better on the best time and no better on the median,
   26.8 to 30.7 ms against 28.7 to 29.9 ms over four sets, for 43.3 MB more
   peak resident set, 620.5 against 577.2 MB. Eight was chosen as the knee. What was not
   measured is the knee for a consumer that reads the genotypes: every number
   of that sweep came from the benchmark whose consumer only adds them up, and
   the cache argument of the regression above says the knee could sit lower for
   the Kosman pass. That measurement is item 1 of section 3.
4. **The frames of a vars file carry no checksum.** arrow-rs writes them with
   `lz4_flex`'s default frame settings, which set neither the per block nor the
   content checksum, verified on all 280 frames of the panel. So a byte changed
   inside a compressed buffer that still decodes to the declared length gives
   other genotypes with no error, which is the one case in popnei where a wrong
   result can pass silently against the owner's rule of 21 September 2026.
   `docs/specs/io_vars.md` records the exposure for the sweep of damaged files
   it ran, and nothing in this review changes it. Closing it means popnei
   writing the buffers itself, which is finding S1 of section 5. This is for the
   owner because it is about right and not about fast. Recommendation: an issue
   rather than work now, decided against the sweep of damaged files that the
   spec already ran, since that sweep is what says how often a damaged file gets
   through, and the change that would close it is the one this review measured
   and did not need.
5. **`+simd128` has no guard, and it is worth 30 per cent of the wasm read.**
   `.cargo/config.toml` sets that flag for both wasm targets, and it is what
   vectorizes the decompression:
   `docs/reports/perf-dists-kosman-2026-09-22.md` measured the wasm read
   falling from 0.162 s to 0.113 s with it, and the wasm reviewer counted 48
   wasm vector instructions in `lz4_flex`'s frame decoder with the flag and
   none without it, for both targets. A `RUSTFLAGS` in the environment of a
   build **replaces** what that file sets, rather than adding to it, and
   nothing in the repository would say so: both wasm reads would simply get
   slower. The fix is a `#[cfg(not(target_feature = "simd128"))] compile_error!`
   in the JavaScript binding crate and the same under the pyodide wheel's build
   script. Recommendation: do it; it is finding L3 of section 5.
6. **Whether the branch is merged.** Three commits, every check of the coding
   skill clean at the head, 1017 cargo tests of the core, 150 of the linear
   algebra crate and 556 pytest tests passing, and the digest of every block of
   every pass identical at 1, 4 and 18 threads.

## 3. The measurement plan

What this review built, which the next one starts from:

1. **`--threads n` and a digest of the pass in `crates/popnei/benches/vars_file.rs`.**
   Every timed section now runs inside a rayon pool of the size the command
   line gives, as `read_vcf.rs` does, and the untimed first pass prints a 64
   bit FNV-1a hash fed, block after block, the place of the block, its counts,
   the bytes of every column it holds, and the names of the chromosome table in
   the order of their numbers. The sum of the genotypes that the benchmark
   printed before cannot see a reordering, and the digest can: reversing the
   window on purpose left the sum at 18173128 and moved the digest from
   `4a928f3a01318f4a` to `4121bc7726267c28`. Nothing of it is inside a clock.
2. **The identity gate, in the module.**
   `the_number_of_threads_does_not_change_the_blocks_of_a_vars_file` reads a
   file of 17 batches, where the chromosome of the first variant of each batch
   first appears in that batch, in pools of 1, 4 and 18 threads, and compares
   every variant field by field, the chromosome numbers of each block and the
   names of the table. Two more tests damage a later batch and hold that the
   first fault in the order of the file is still the error the reader gives and
   that the blocks before it still come out first.
3. **The knee of the window for a consumer that reads the genotypes**, which
   item 3 of section 2 asks for: `cargo bench --bench kosman_dists --
   <panel> --runs 3` at 18 threads with `BATCHES_AT_ONCE` of 2, 4, 8 and 18,
   one build per window, interleaved, the metric being the vars file pass and
   the control being the same calculation over blocks already in memory, which
   no change of the reader can move.
4. **The allocation count still does not exist.** Section 4.2 of
   `docs/reports/perf-io-2026-09-22.md` designed the `dhat` feature in
   September and nobody built it. Three findings of section 5 would be gated on
   a count that is the same on every run instead of on a wall time this machine
   moves by a fifth, and the parallel reader itself should be gated on its
   allocations not growing with the threads. Note when it is built: `dhat`'s
   allocator serialises every allocation, so a threaded run under it gives
   counts and no timing.
5. **The two panels at 10000 individuals.** `wide10000.vcf`, 3000 variants of
   10000 individuals, is in the benchmark directory and no benchmark has ever
   read it. At that shape a batch is 250 variants, so the panel of 100000
   variants is 400 batches and every per batch cost of the reader is paid
   twenty times as often.

## 4. The build configuration

Not touched, and the reason is the same as in the two reviews before this
one: `[profile.release]` has `overflow-checks = false` and `debug =
"line-tables-only"`, no link time optimization and 16 codegen units, and
71 per cent of the self time of a pass is inside `lz4_flex`, reached
through `arrow-ipc`. Link time optimization and one codegen unit were each
measured twice by earlier reviews with nothing outside 0.004 s
(`docs/reports/perf-kinship-2026-09-23.md`). The methodology reviewer asks
for them again on this pass, because cross-crate inlining into `lz4_flex` is
the case those measurements did not cover, and because a baseline that moves
would have to be restated before the 21 ms is compared with anything; it is
finding L4 of section 5. `cargo asm` cannot answer it, since it appends its
own `-Ccodegen-units=1` to every invocation, so it is `objdump -d` on the
built benchmark.

## 5. The findings that were not acted on

### Hot-path

**H1. The genotypes of a batch are copied into the block after arrow has
already materialised them.** `crates/popnei/src/io/vars.rs:2524`.
`_platform_memmove` is 1199 of the 16130 samples of the profile, and the
copy is 10000000 bytes per batch, 200 MB per pass. arrow-ipc builds the
decompressed buffer with `Buffer::from_vec`, whose allocation
`Buffer::into_vec::<i8>()` gives back with no copy when the layout matches,
the offset is 0 and nothing else holds the buffer; taking the batch apart
with `into_parts` down to the values buffer drops the count to 1, and
`into_vec` returns the buffer back as an error when it cannot, which is the
fallback. It also makes the open finding L6 of the September review moot: with the copy gone there is nothing
left to fuse the scan of the smallest allele with. At 18 threads the copy is
on the calling thread, so it is about 4 ms of the 27.23 ms that a pass now
takes, 15 per cent, where before the change it was 4 per cent of 112 ms.
Measure the wall time of the pass at 18 threads, best of 20, interleaved, and
count that the no-copy arm was taken for all 20 batches, because a silent
fallback loses the gain with nothing to say so. Cost: `genotypes` takes the
column by value, the other columns are read before it, and a future arrow
that slices its buffers would make it copy again.

**H2. The all-ones validity bitmap of the genotypes is decompressed and
thrown away: 4.5 ms of a 110 ms pass.** `crates/popnei/src/io/vars.rs:995`.
popnei builds the `gts` arrays with no nulls, and `write_array_data` of
arrow-ipc 60 writes an all-valid bitmap anyway, 1250000 bytes per batch, 25
MB of the 227690000 bytes a pass decompresses; `create_primitive_array` then
takes that buffer only when the field node says there are nulls, so it is
decompressed and dropped unread. Measured on one such frame, 0.227 ms for
1250000 bytes out of 4913 bytes in, so the 11 per cent of the bytes is 4 per
cent of the time, and at 18 threads about 1 ms of 27. It is not a format
rule: `bigcalled.pynei.vars`, written by pyarrow over the same genotypes, has
no such buffer, and arrow-rs's own reader handles a file without one. Removing
it means an upstream change to arrow-rs, or popnei writing the buffers
itself, or a `gts` column of `FixedSizeBinary` instead of a list of `Int8`,
which is a change of the format. Filed, not recommended at 1 ms.

### Likely

**L1. The bytes of every batch are zeroed before they are read over.**
`crates/popnei/src/io/vars.rs:2920`, the open finding L5 of the September
review, refiled with the number that lowers it. `try_reserve_exact` then
`resize(len, 0)` then `read_exact` writes 3.89 MB twice per batch, 77.8 MB
per pass, and `__bzero` is 196 of 16130 samples; but arrow's own path zeroes
227690000 bytes per pass for the same reason, so popnei's share is a quarter
of that 1.2 per cent. `source.take(len).read_to_end(&mut bytes)` fills the
spare capacity without zeroing, with the count compared against `len` for the
error of a file that was cut short. Two lines; gate on the counted bytes, not
on a wall time.

**L2. `iter_blocks` cannot use the read ahead reader, and the eleven
consumers that are not the association study do not use it either.**
`crates/popnei/src/block.rs:1662` and `crates/popnei-python/src/source.rs:162`.
`with_one_block_ahead` borrows its reader inside a scoped thread, so the whole
pass has to be one closure; a Python iterator hands control back between
blocks and cannot be one. `BlockReader` is `Send` and the handle already keeps
its own chromosome table, so an owning form, a reader moved into
`std::thread::spawn` and returned as a `Box<dyn BlockReader>`, is storable and
needs no borrow. With the reader now four times faster this is worth less than
it was when the September reviews filed it: the read is 27 ms of a pass whose
consumer does 100 ms of work. Measure with `crates/popnei/benches/time_stats.py`
over the panel, in a mode that does one numpy pass per block, and keep it only
above 10 ms.

**L3. Nothing guards `+simd128`**, which item 5 of section 2 states in full.

**L4. Link time optimization and one codegen unit have never been built for
this pass**, which section 4 states in full.

**L5. One `String` per variant in the ids of a batch.**
`crates/popnei/src/io/vars.rs:2328`, 5000 per batch and 100000 per pass with
every field asked for. Against it: malloc and free together are about 54 of
the 16130 samples of the profile, which bounds the whole allocator at 0.33
per cent, and the two passes of the benchmark, the genotypes alone and every
field, differ by 3 ms in 110. Measure with the count when the `dhat` feature
exists; leave alone until then.

### Speculative

**S1. popnei decompressing the buffers of a batch itself.** This is the route
the owner named first, and the measurement that sizes it is in section 6.1:
walking a frame's blocks and calling `lz4_flex::block::decompress_into` is
4.534 ms against arrow's 4.759 ms on one thread, 4.7 per cent, and 1.972 ms
on the threads of rayon, 2.41 times, which is the ceiling the three blocks of
that frame allow. It would also skip the bitmap of H2, decompress straight
into the block's own vector, which removes the copy of H1, and let the frames
be written with a checksum, which is item 4 of section 2. Against it: popnei
would parse the lz4 frame header and build the arrays of a batch itself, with
a fallback to arrow for any frame it does not recognise, so `what_the_buffers_hold`
becomes load-bearing rather than a check, and it needs `lz4_flex` as a direct
dependency with its default features off, since cargo unifies features and the
defaults would turn `safe-decode` on for arrow's path as well. The window of
section 2 already scales past 2.41 times with none of that. Worth running only
if the window's memory is refused, or for the checksum.

**S2. A mapped file instead of a copy of each batch's bytes.**
`crates/popnei/src/io/vars.rs:2907`. A pass copies 77.8 MB out of the page
cache into a fresh 3.89 MB vector per batch, and
`arrow_buffer::Buffer::from_custom_allocation` would let a mapped range be a
buffer with no read at all. `crates/popnei/src/lib.rs:60` is
`#![forbid(unsafe_code)]` and both the map and that constructor are `unsafe`,
so it needs a crate of its own as `popnei-linalg` is, a second shape of
source beside `Read + Seek`, and it gives up the error of a file cut short
for a signal if the file is truncated while it is read. Measure first whether
it is worth anything: the benchmark reads from memory and `from_path` reads
from a warm file, and the difference between the two is the read plus the
copy.

**S3. The number of variants a batch holds was chosen for pyNei.**
`crates/popnei/src/block.rs:28`. 5 million genotypes is
`DEF_NUM_GTS_PER_CHUNK` of pyNei's `config.py`, which at ploidy 2 and 1000
individuals is 5000 variants and 10 MB, and the comment there says nobody has
measured it for popnei. A batch's 15 MB of working set against a 16 MB second
level cache shared by six cores is what the regression of section 2 was
about, and the window multiplies it. Write the panel again at 500, 1000, 2000
and 10000 variants a batch and report milliseconds per million genotypes
beside each file's bytes, since smaller batches compress worse.

### Notes

**N1. A codec faster than lz4, which is (d), does not exist inside the
format.**
`arrow-ipc` 60 defines two body codecs and no others, LZ4_FRAME and ZSTD, so
a faster codec means a file arrow cannot read. It collapses into zstd, which
"The compression" of `docs/specs/io_vars.md` measured as **slower** to read,
26 ms against 19 ms on the smaller panel, and whose build was refused on 20
September because the `zstd` crate wraps the C library; or into (c) of
section 1, no compression.

**N2. An uncompressed vars file, which is (c), is the only one of the four
that helps a browser tab, and it is not measured here.** The spec's own table has arrow-rs
reading the smaller panel in 1 ms uncompressed against 19 ms with lz4, and the
file is 2.9 times larger. The reader already reads such a file and needs no
change; the writer needs an option, the spec a paragraph and both a test. For
a tab it is the whole of the read, since there are no threads there and no
disc to save, but the in-memory route holds the file in the tab's memory, so
2.9 times the bytes, and the file route makes 2.9 times as many range reads
of a picked file. An option for a caller, never the default, measured
natively with `vars_file` and in a tab with
`js/popnei/bench/time_kosman_dists.mjs`.

**N3. `checked-decode` of `lz4_flex` does not exist.** The context of this
review, carried from `docs/specs/io_vars.md`, says that `arrow-ipc` takes
`lz4_flex` with its default features off so that `safe-decode` and
`checked-decode` are both off. The string `checked-decode` appears nowhere in
`lz4_flex` 0.14.0's sources: the bounds checks of its fast decoder are
unconditional and a malformed stream is an error, not a read past a buffer,
and `safe-decode` only swaps one implementation for another with the same
checks. `cargo tree -p popnei -f "{p} {f}"` prints `lz4_flex v0.14.0
alloc,frame,std` and is still the guard worth keeping, for the size and the
speed of the build and not for safety.

**N4. Which error a file with two faults gives is now pinned by a test** and
was not before. The first fault in the order of the file wins, at three
levels: the batches in order, then the columns in the order of the file, then
the rows.

## 6. What the experiments showed

### 6.1 One buffer, decompressed three ways

The measurement that chose the design, before any code was changed. One
batch's genotype buffer of the panel, cut out of the file with its frame
intact, 10000000 bytes of output from 3736410 bytes in, three lz4 blocks of
4 MiB, 4 MiB and 1.6 MB; best of 20 runs, one thread except where it says:

| | best | rate |
|---|---|---|
| `FrameDecoder::read_to_end`, which is what arrow-ipc calls | 4.759 ms | 2.10 GB/s |
| the frame's three blocks with `decompress_into` | 4.534 ms | 2.21 GB/s |
| the same three blocks on the threads of rayon | 1.972 ms | 5.07 GB/s |

The three gave the same 10000000 bytes. The bitmap buffer beside it, 1250000
bytes of all ones from 4913 bytes in, takes 0.227 ms. So a batch's
decompression is 4.99 ms of the 5.5 ms a batch costs, the frame wrapper's
extra copies are 0.2 ms of it, and splitting a frame stops at 2.41 times
because the largest of its blocks is 4 MiB of the 10 MB. That is what sent the
review to the window across batches, which is bounded by the cores and not by
the blocks of a frame.

### 6.2 The window across batches: applied, at `7051bd4`, `7c380b8` and `5c682ea`

The three commits are the benchmark's `--threads` and digest, the window, and
the window bounded by the pool. The numbers are in section 2. The sweep of the
window at 18 threads, best and median of 10 runs of the pass with the
genotypes alone over the panel:

| batches at once | best | median | peak resident set |
|---|---|---|---|
| 1 | 107.76 ms | 108.34 ms | 527.3 MB |
| 2 | 62.29 ms | 62.39 ms | — |
| 4 | 38.74 ms | 39.07 ms | 552.0 MB |
| 8 | 27.79 ms | 28.75 ms | 577.2 MB |
| 18 | 26.23 ms | 30.65 ms | 620.5 MB |

### 6.3 What was refuted

(a) of section 1, the buffers of one batch on the threads, cannot
pay: three reviewers counted the buffers of a batch and 10000000 of its
11384500 bytes are one buffer, so splitting a batch by its buffers has a
ceiling of 1.13 times however many cores are used. (d) is closed by the
format, N1 above. The frame decoder's extra copies, which the byte counts made
look like 300 MB per pass, are 0.2 ms of 4.76 per buffer, section 6.1.

## 7. Seen outside the scope

- **The write compresses on one thread between two parallel reads**, the open
  finding L7 of `docs/reports/perf-io-2026-09-22.md`, still open: 202.58 ms
  for the panel, of which `lz4_flex`'s compression is 2337 of the 16130
  samples of the profile's window. The concurrency reviewer found a cheaper
  form than the writer thread that finding designed: wrapping the loop of
  `write_vars` in `with_one_block_ahead` moves the read instead of the
  compression and buys the same overlap for one line, since the handle's
  chromosome table is the one `write_block` looks its numbers up in. It needs
  the target of section 4.4 of that report, which the writer still has not
  got.
- **`Schema::project` clones the metadata of the schema once per batch**,
  and that metadata holds the `popnei` key with every individual's name: 7 to
  9 KB per batch at 1000 individuals and about 90 KB at 10000, where a pass
  has 400 batches. A second schema without the metadata, kept for decoding,
  removes it. No benchmark has ever run that shape.
- **The chromosome of every variant crosses into Python as one object per
  variant**, which is finding L3 of `docs/reports/perf-gwas-2026-09-24.md`,
  unchanged.

## 8. What the code already does well

- **Every buffer of a batch is checked against the message before arrow-rs
  reads any of it**, `message_fits` and `buffers_fit` of
  `crates/popnei/src/io/vars.rs`, with saturating arithmetic throughout and a
  refusal rather than a widened bound on any conversion that fails. That is
  what made the window safe to write: the decode of a batch needed no new
  check.
- **The writer's side is genuinely zero-copy**: `ScalarBuffer::from(gts)`
  adopts the block's vector and `Buffer::from(bytes)` adopts the batch's, as
  the spec claims.
- **The scan for an allele below the missing one is as vectorized as it
  gets**: `cargo asm` prints the fold over `i8::min` as four NEON
  accumulators, 64 bytes an iteration, with no bounds check, and the fold is
  order-free so a parallel or vectorized form gives the same answer.
