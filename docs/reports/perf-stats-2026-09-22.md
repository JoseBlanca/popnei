# Performance review: the stats module, 22 September 2026

The per variant pass of popnei's stats module takes 0.479 to 0.483 s on
one thread over the 100000 variants of 1000 diploid individuals of
`/Users/jose/devel/popnei-bench/big.vars`, where "Speed" of
`docs/specs/stats.md` asks for 0.25 s. It is the one number of four that
the plan of that module missed, and `docs/reports/stats-measurement.md`,
which measured it, handed it to this review. This document says where
that time is, what was tried, what each experiment gave, and what is
left for the owner.

The words this document uses. A **pass** is one whole read of a dataset
with whatever is calculated over the variants on the way. A **block** is
the genotypes of many variants as one flat array of bytes, one row per
variant, a row being the individuals one after another, each taking as
many bytes as the ploidy. The **row loop** is `add_the_rows` of
`crates/popnei/src/stats.rs`, which walks the rows of one block and adds
what it counts into the accumulators of the pass. The **genotype count**
of a variant counts, over the individuals asked for, how many genotypes
are missing, heterozygous and called; the **allele count** counts how
often each allele was called, into an array of 128 counters. A **vars
file** is popnei's own format, arrow IPC with the genotypes compressed
with lz4.

## The scope and its limits

Reviewed at 2742918 on the branch `perf/stats`, which is `plan/stats`
plus one commit: `crates/popnei/src/stats.rs`, `variant.rs`,
`io/vars.rs`, `block.rs`, the two binding crates' `stats.rs`,
`python/popnei/stats.py`, the timing harness and the build
configuration. Nine categories were sent, one reviewer each:
methodology, numbers, allocations, data layout, concurrency, hot loops,
input and output, linear algebra and wasm, and the Python boundary. The
linear algebra half of one of them does not apply, since there is no
matrix in this code, and its reviewer said so and spent no time on it.

The targets, from "Speed" of `docs/specs/stats.md`: a whole pass over
that file through `open_vars`, 0.25 s on one thread and 0.15 s on 18
cores, for the per variant pass with the five statistics and no
populations and for the per individual statistics. The sizes the code is
for, from `docs/objectives.md`: tens of individuals to 10000, tens of
thousands of variants to a million, and a dataset of 50 populations is
ordinary. The machine is the owner's Apple M5 Pro, 18 cores, 64 GB,
macOS 27.0, and a browser tab, which has one thread.

What was not reviewed: the VCF reader, the filters, the distances and
the principal components. What could not be measured: this machine has
no hardware counters, so no cache miss, branch miss or instruction per
cycle was counted anywhere below, and every claim about a cache or a
mispredicted branch is an inference with an experiment beside it that
would settle it.

## The verdict

Run the experiments. The site is known, the change with the largest
expected gain is cleared by the reviewer whose job is to refuse the
unsafe ones, and three categories converge on the same mechanism. What
each experiment gave is at the end of this document, under "What the
experiments gave".

## Where the time is

A sampling profile of the Python process, `sample` for 8 s, one thread,
the major allele frequency alone, taken at 2742918 by the orchestrator.
The threads of the rayon pool sleep and are filtered out, as
`.claude/skills/performance-review/profiling_environment.md` says. The
self time of the working thread, as printed:

```
Sort by top of stack, same collapsed (when >= 5):
        __workq_kernreturn  (in libsystem_kernel.dylib)        6067
        __psynch_cvwait  (in libsystem_kernel.dylib)        5784
        _RNvNtCsc7IbC0cTyC_6popnei5stats12add_the_rows        4272
        _RNvXs_..._8lz4_flex5frame10decompress...              1452
        swtch_pri  (in libsystem_kernel.dylib)                  339
        _platform_memmove  (in libsystem_platform.dylib)         99
        read  (in libsystem_kernel.dylib)                        79
        _RNvNtNtCsc7IbC0cTyC_6popnei2io4vars18block_of_the_batch 38
        __bzero                                                  20
        _xzm_xzone_malloc_tiny                                    6
```

Of the 5966 working samples the row loop is 4272, 71.6 in 100, and
everything on the read side, the lz4 decompression, the copy, the `read`
system call and the building of the block, is 1694, 28.4 in 100. Over
the 0.383 s of that pass those shares are 0.274 s and 0.109 s, against
the 0.279 to 0.281 s and 0.102 to 0.104 s that
`docs/reports/stats-measurement.md` measured by subtracting whole passes
from each other. Two measurements taken different ways agree within
0.007 s, so the split the rest of this review rests on is sound. The
reviewer of the Python boundary made that comparison and also counted,
with `cProfile`, that a whole pass makes 3217 Python calls and 13 pandas
constructors over at most 41 values each: the opening of the file and
the building of the result are under about 0.3 ms of a 0.4 s pass and
are not worth looking at.

The profile names the function and not the line, because everything the
row loop calls is inlined into it. The findings below therefore rest on
the assembly of that one function, which three reviewers read with
`cargo asm`, and on the arithmetic of the wall times.

What the row loop spends it on. `docs/reports/stats-measurement.md`
timed the pass with one statistic asked for: the genotype count adds
0.089 to 0.091 s to the read and the allele count adds 0.279 to 0.281 s,
three times as much over the same 100000 rows of 2000 alleles. That
asymmetry between two walks of the same bytes is what this review had to
explain, and it has an answer.

## The findings

### H1 The allele count is a chain of dependent stores, not a stream

`crates/popnei/src/variant.rs:506`, the counting loop of
`count_the_alleles`. Hot path, high confidence on the site and medium on
which of two mechanisms dominates. Found independently by the hot loops,
the data layout and the wasm reviewers, and cleared by the numbers one.

The loop compiles, inside the row loop, to this, as `cargo asm` printed
it:

```
LBB321_45:
	ldr w10, [x14, w9, uxtw #2]
	add w10, w10, #1
	str w10, [x14, w9, uxtw #2]
```

It reads a counter, adds one and writes it back, at an address the datum
chooses. The file has two alleles a variant, so almost every allele hits
the counter the one before it just wrote, and each increment waits for
the store of the previous one to reach the load through store to load
forwarding. That is a serial chain of about four to six cycles per
allele. The arithmetic agrees: 0.28 s over 2e8 alleles is 1.4 ns each,
which at this machine's clock is about six cycles. The genotype count
walks the same bytes into registers, with no such chain, and costs a
third as much.

The fix the three reviewers propose is the standard one for a histogram
whose values repeat: count into four interleaved arrays of counters,
choosing the array by the position rather than by the value, and add the
four together at the end of the row. That gives four independent chains
and lets the loop run at its issue rate instead of at its latency.

The numbers reviewer cleared it: these are integer counts, so splitting
and re-merging them is exactly equal whatever the order, and the
existing refusal of more alleles than a count holds already bounds every
partial below what a 32 bit counter takes. Three invariants have to
survive: the counters read as all zero at the start of every row, every
allele below the missing one still raises its error, and the called
alleles count the entries that are not missing and nothing else.

The measurement: the assembly shows four independent chains to four base
registers, which is the same on every run; then the pass with the allele
frequency alone, best of five at a quiet load, keeping the change only
if the pass falls by 0.05 s or more.

### H2 The inner loop builds an error value it throws away

`crates/popnei/src/variant.rs:610` and `:400`. Hot path, high confidence
on the mechanism. Found by the hot loops and the data layout reviewers,
the second of which said to look at it first.

Looking up one individual's genotypes ends in `.ok_or(Error::...)`,
whose argument Rust evaluates before it knows whether it is needed. So
every one of the hundred million lookups of a four population pass
writes an error value to the stack and calls its destructor, which is
out of line because other cases of that enum own strings. The assembly
of the row loop has `bl core::ptr::drop_glue::<popnei::error::Error>` on
the success path of the per individual loop, with the spilling and
reloading of the loop's registers around the call. `.ok_or_else` builds
it only when it is needed.

This costs nothing in complexity: two closures. It does not touch the
pass with no populations, which does not take that path, so it is aimed
at the 0.120 to 0.140 s that four populations cost.

### H3 The read sets a floor that no thread count can lower

`crates/popnei/src/stats.rs:1287`. Hot path, high confidence. The
concurrency reviewer.

The pass reads a block and then forks over its rows, one after the
other. The file holds 20 batches of 5000 variants, so per block that is
about 5.1 ms of reading on the calling thread against 1.85 ms of
calculation spread over the pool. The read alone is 0.102 s at one
thread and at 18, and the whole pass on 18 cores is 0.139 s: the
calculation has come down from 0.377 s to 0.037 s, a tenth rather than
an eighteenth, and what is left is the read. No arrangement of the
threads can take the pass below 0.102 s while the read is serial.

`docs/architecture.md` already prescribes the cure at its line 155, a
read ahead thread one block ahead, and says it waits for a spec of its
own; no such thread exists. This is for the owner and not for this
review, for two reasons. It needs that spec, and it changes what the
spec's one thread row means, since a native run would then own two
operating system threads. It also does nothing for the missed target,
which is the one thread number.

### H4 A reciprocal would be faster and would change a published count

`crates/popnei/src/stats.rs:450`, `:502`, `:623`, `:650`. Hot path, high
confidence. The numbers reviewer, as a refusal.

The row loop contains five divisions, which invite the usual
optimisation of computing one reciprocal per row and multiplying by it.
It must not be done here. The reviewer measured that for 2504 of the
counts from 1 to 20000, a count multiplied by its reciprocal gives
0.9999999999999999 and never more than 1. A variant where every called
allele is the same would then get a major allele frequency just below 1
instead of exactly 1, and the polymorphism counts call a variant
variable when that frequency is below 1. That count is compared with
pyNei exactly, not within a tolerance. The plain expected heterozygosity
of the same variant would go from exactly 0 to 2.2e-16.

This finding is here so that the next reviewer does not spend an
afternoon on it.

### L1 The counters are cleared and scanned for 128 alleles to hold two

`crates/popnei/src/variant.rs:483` and `:584`,
`crates/popnei/src/stats.rs:501`, `:620`, `:642`. Likely, medium
confidence. The data layout reviewer, with the allocations and the
numbers reviewers pointing at the same lines.

For every row and every population the pass clears 512 bytes of
counters, and then the major allele frequency scans all 128 of them for
its maximum and each expected heterozygosity scans all 128 for the ones
that are not zero. The file has two alleles a variant. That fixed cost
does not grow with the individuals, which is why it is felt as the cost
of a population rather than of a genotype: four populations cost 0.120
to 0.140 s over the same genotypes, 0.4 to 0.47 µs per variant per
population. At the 50 populations of 20 individuals that the objectives
call ordinary it would dominate.

Carrying the largest allele seen beside the counters, and clearing and
scanning only up to it, makes all three costs two entries. The numbers
reviewer confirmed it is bit for bit the same, since the entries above
the largest are zero and both expected heterozygosities already skip
zero counts.

### L2 The chunk size is fixed in rows, so wide datasets use few cores

`crates/popnei/src/stats.rs:1007`. Likely, medium confidence. The
concurrency reviewer, with the allocations and the methodology reviewers
noting the same constant.

A block is divided into chunks of 64 rows. The comment justifies 64 as
128000 genotypes for 1000 diploid individuals, which it is. But the
block itself is sized in genotypes, so at 10000 individuals a block is
500 rows, which is 8 chunks: at most 8 of 18 cores ever run, and each
chunk is ten times the work the comment reasons about. Stating the grain
in genotypes instead leaves this dataset unchanged at 64 rows and gives
42 chunks per block at 10000 individuals.

The numbers reviewer adds an argument for a larger chunk that has
nothing to do with threads: the total is a sum over chunks of sums over
rows, so a larger chunk narrows the worst case error of the mean, and at
a million variants 1000 rows per chunk is near the optimum. It also
warns that the grain must never be computed from the thread count, or
the result would stop being the same at every thread count, which a test
asserts.

### L3 Every pass decompresses 25 MB of an all ones mask and drops it

`crates/popnei/src/io/vars.rs:995` for the writer and `:1623` for the
reader. Likely, high confidence. The input and output reviewer, which
parsed the file itself.

Arrow's writer emits a validity bitmap for every array that can have
one, filled with ones when nothing is null, whatever the field says; the
reader decompresses it and then drops it, because the null count is
zero. The reviewer walked the 20 batches of `big.vars` and counted
25012500 of the 225012500 bytes decompressed per pass, 11.1 in 100, and
20 of the 60 lz4 frames.

This also contradicts the spec. "What it holds" of
`docs/specs/io_vars.md` says popnei writes the genotypes with no nulls
so that the mask is not paid for, and popnei's own files carry 1250000
bytes of mask a batch. That sentence is wrong and should be corrected
whatever is done about the cost.

The three ways out all have a real price: changing what the column is in
the file, which is a format change and costs the shape the spec sells to
pandas and R; patching the reader of arrow upstream; or patching its
writer, which would not help files that already exist. None is this
review's to take.

### S1 The chunk results of a block are all held live at once

`crates/popnei/src/stats.rs:1400`. Speculative. The allocations and the
concurrency reviewers.

Each chunk builds its own accumulators and they are collected into one
vector before being added in order, which is what makes the result the
same at every thread count and must be kept. The live bytes are the
chunks times the populations times the statistics times the bins. At the
default 40 bins that is about 5 MB a block at 50 populations, which is
nothing. At the largest bin count the histogram now accepts, 100000,
which this plan's code review introduced as a bound against a crash, it
is about 160 MB per chunk and gigabytes per block. Nobody has run that,
and it is a failure mode rather than a speed finding.

## The measurement plan

In the order in which they unblock each other.

1. A Rust benchmark of the pass, which does not exist. The crate has
   benchmarks of the VCF reader, the vars file and the filters and none
   of this module, so every number above came from a Python harness in
   which 28 in 100 of the time is the read and whose spread between runs,
   on a loaded machine, is of the same size as the wins being sought. The
   benchmark builds one block in memory, wraps it in a reader written in
   the benchmark and times the pass over it, with the statistics and the
   populations as arguments.
2. The error value of H2, gated on the count of the destructor call in
   the assembly going to zero, which is the same on every run, and then
   on the pass with four populations.
3. The four counting lanes of H1, gated on the assembly showing four
   chains and then on the pass with the allele frequency alone.
4. The build configuration. The numbers reviewer established that it
   cannot move a float result here: it counted zero fused multiply adds
   in the row loop at every optimisation level and with the native
   processor selected, and found the output with the native processor
   byte identical to the default. So linking across crates and one code
   generation unit are pure speed experiments, and the second is the one
   worth timing, since one of the counting functions is still an out of
   line call.
5. The prefix of L1 and the grain of L2, each on the benchmark and then
   on a sweep of the populations, which nobody has run.

## The build configuration

`overflow-checks = false` in the release profile was counted rather than
assumed: with the checks on, the counting loop is 16 instructions per
allele instead of 11, with three overflow branches added, and the two
counting functions stop being inlined into the row loop altogether. It
costs no safety, because the lint denies the plain operators and the two
places that carry an exception state a bound their callers establish.
The consequence to remember is that a timing taken under the test
profile, which is built at a lower optimisation level with the checks
on, measures a different loop and is not comparable with a release one.

`debug = "line-tables-only"` was added to the release profile by this
review, at 2742918, because a release build with no debug information
gives a profile with no function names. It costs nothing at run time.

The vector instructions of WebAssembly are set in `main` and not on this
branch, which predates that merge: `main` carries them for both wasm
targets, added by the performance review of the Kosman distances, along
with the floor that decision sets for which browsers popnei runs in. Any
measurement of this module in a browser has to be taken after that
merge, or it measures a build nobody ships.

## Seen outside the scope

The spec of the vars file claims the genotypes are written without a
validity mask and they are not, under L3 above.

`iter_blocks` ignores the block size a user asks for and lets the
reblocking join and cut instead, so a user who asks for a size other
than the file's copies every genotype once, 200 MB for this file. It
does not touch any number here, because the file's batches happen to be
written at exactly popnei's default size, which is also why the read
baseline of the measurement is a fair one.

`import popnei` imports pandas, so a user who only reads blocks, and a
browser tab, pay for it.

## What the code already does well

Three patterns worth copying, each named by a reviewer that went looking
for the opposite.

The row loop allocates nothing. The counters are an array on the stack,
reused for every row and every population, and the profile agrees: the
allocator is 6 samples of 5966.

The rows are the outer loop and the populations the inner one. A row is
2000 bytes and stays in the first level cache while the populations walk
it; the other order would stream the whole chunk once per population.

The serial path that a browser runs and the parallel one are faithful
mirrors: same chunk boundaries, same order, same prologue. A reviewer
went looking for drift between them, because the code review of the plan
had found the two block prologues had drifted, and found none. The
serial path allocates less than the parallel one.

## What the experiments gave

This section is filled as each experiment is run.
