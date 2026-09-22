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

Reviewed at the commit 2742918 on the branch `perf/stats`, which is
`plan/stats`, the branch that built this module, plus one commit: `crates/popnei/src/stats.rs`, `variant.rs`,
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

Each is numbered so that the rest of the document can name it. The
letter is how sure the review is that the change would give a gain worth
its complexity, which is what a performance finding is about, and not
how serious a defect it is: **H** for a site a profile or a benchmark
names, with a clear mechanism for the gain; **L** for one matched by
pattern, with a plausible call frequency and no profile yet, which is
where most of them start; **S** for one where it is not even clear that
the site is hot, filed so that it is known. One of them, H4, is a
refusal rather than a proposal.

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
the counter that the allele before it just wrote.

That is what makes it slow, and it is worth spelling out, because the
whole finding turns on it. A processor runs the instructions of a loop
several at a time and out of order, so a loop like this one normally has
many alleles in flight at once. It can do that only where the
instructions do not depend on each other. Here they do: to add one to a
counter the processor must first know what is in it, and what is in it
is what the previous allele wrote a moment ago, which has not yet
reached memory. The hardware has a shortcut for exactly this, handing
the value straight from the pending write to the waiting read, which is
called store to load forwarding, and it still takes about four to six
cycles. So the alleles cannot overlap: they queue, each waiting for the
one before it, and the loop runs at the speed of that queue instead of
at the speed the processor could issue its instructions.

The arithmetic agrees with that reading. The pass spends 0.28 s on 2e8
alleles, which is 1.4 ns each; the cores of this machine run at about
4 GHz, so 1.4 ns is about six cycles, the length of one link of that
queue. The genotype count walks the same bytes but keeps its running
totals in registers, where there is no such queue, and costs a third as
much.

The fix the three reviewers propose is the standard one for a histogram
whose values repeat: count into four interleaved arrays of counters,
choosing the array by the position rather than by the value, and add the
four together at the end of the row. Four consecutive alleles then touch
four different arrays, so none of them waits for the one before it, and
the loop can overlap its work again instead of queueing.

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

The change is to carry, beside the counters, the largest allele that has
been seen, and to clear the counters and run both scans only up to that
one instead of to the end. On a file of two alleles all three of those
costs then touch two entries instead of 128. The numbers
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
nothing to do with threads. A mean here is a sum of floating point
numbers and every addition rounds, so the error of a sum grows with how
many additions a value passes through. The pass adds the rows of a chunk
into one running total and then adds the chunks together, so the longest
such chain is the rows in a chunk plus the number of chunks. A bigger
chunk shortens one of those and lengthens the other, and the two are
equal when the chunk is the square root of the rows: at a million
variants that is about 1000 rows, against 64 today. It also
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
4. The build configuration: optimising across crate boundaries, and
   compiling the crate as one piece instead of the sixteen the compiler
   splits it into by default, which lets it inline more. The numbers
   reviewer established that neither can move a float result here. The
   one way a compiler changes a floating point answer without being
   asked is by fusing a multiply and an addition into a single
   instruction that rounds once instead of twice; it counted how many of
   those the row loop holds at every optimisation level and with the
   processor of this machine selected, and the answer was none every
   time, with the output for this processor byte identical to the
   default. So both are pure speed experiments, and the second is the one
   worth timing, since one of the counting functions is still called
   rather than inlined.
5. The prefix of L1 and the grain of L2, each on the benchmark and then
   on a sweep of the populations, which nobody has run.

## The build configuration

`overflow-checks = false` in the release profile was counted rather than
assumed: with the checks on, the counting loop is 16 instructions per
allele instead of 11, with three overflow branches added, and the two
counting functions stop being inlined into the row loop altogether. It
costs no safety. popnei denies the plain arithmetic operators on integers
throughout, with a lint, so an addition that could overflow does not
compile and the code has to say what it wants to happen instead; only two
places in this counting code are exempted from that lint, and each
carries a written bound that its callers do establish before they call.
So nothing here relies on the release build checking an overflow.
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

`iter_blocks`, the method a user calls to get the genotypes themselves
in blocks, ignores the block size they ask for. A reader further down the
chain then cuts and joins the blocks the file gave into the size that was
wanted, and that copies every genotype once, 200 MB for this file, for
any size other than the one the file was written at. It
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

Four experiments were run, one at a time, each delegated, each taking
its own baseline first. Two changes were kept and two were closed with
no gain, which is what a measurement plan is for.

The pass this review exists for, the five per variant statistics with no
populations over `big.vars` on one thread, was 0.477 to 0.483 s when the
review began and is 0.285 to 0.290 s now. The target is 0.25 s, so it is
still 0.035 to 0.040 s over, and 0.102 s of what is left is the read of
the file, which no experiment here touched.

### Applied: the alleles of a variant are counted into four arrays

H1, commit fc16d43. The counting loop now walks four interleaved arrays
of counters, choosing the array by the position of the allele, and adds
the four together at the end of the row. The assembly shows four
independent read, add and write chains to four base registers where
there was one, with no bounds check added. Three tests were added, each
shown to fail on a deliberate break: a dropped remainder allele, a
dropped array, and a merge that assigns where it should add.

On the benchmark, one thread, best of five at a load average of 1.5:

| | before | after |
|---|---|---|
| the allele frequency alone, one population | 0.274 s | 0.087 s |
| the five statistics, one population | 0.371 s | 0.180 s |
| the five statistics, four populations | 0.418 s | 0.318 s |

Over `big.vars` from Python the five statistics with no populations went
from 0.477 s to 0.285 s and the allele frequency alone from 0.381 s to
0.188 s. The same figures came out at load averages of 1.0, 1.6 and 2.4,
so they are not an artefact of the machine being busy. The mechanism the
three reviewers proposed is therefore the right one: the cost was the
chain of dependent stores and not the count of instructions.

It cost a constant, a type alias, two private functions and one more
exemption from the lint that denies plain arithmetic on integers, with
the bound that makes it safe written beside it.

What the merge of `main` did to these numbers, on 22 September 2026.
`main` had gained a fast path of its own for the same loop: a variant
whose alleles are only the missing one, 0 and 1 is counted by two
comparisons and two additions per allele, in runs of 255 with counters of
one byte, and never touches the table of 128. The merge kept both, the
fast path first and the four arrays as what it falls back to. Every
variant of this benchmark has two alleles, so the fast path takes every
row of it and the four arrays never run here: on one thread, best of nine
at a load average of 1.5, the five statistics with one population went
from 0.184 s to 0.112 s and the allele frequency alone from 0.089 s to
0.015 s. The four arrays are what a variant of three alleles or more
still takes, and what `count_alleles_of` takes for a population that is
not every individual: the pass with four populations, which has no fast
path, is 0.328 s.

### Applied: the lookup of an individual builds no error it throws away

H2, commit 4b977be, at one of the two sites the finding named. The count
of the destructor call inside the per individual loops went from three
to zero. On the benchmark the four population pass went from 0.507 s to
0.416 s and from Python from 0.595 s to 0.526 s, with the pass with no
populations unmoved at 0.376 and 0.477 s, which is the control the
change should not touch.

The second site was measured and not kept: it is called once per row
rather than once per individual, it gave the four population pass the
same 0.416 s, and it cost the pass with no populations, which is the one
that misses its target. So the finding is applied at one site of two.

It cost one exemption from a lint: the one that asks for the eager form
wherever the value is cheap to build, which is what put the code in that
shape in the first place and will do the same to the next site.

### Closed with no gain at the sizes that matter: the prefix of the counters

L1. Built, measured and not kept; the patch is kept in the session's
scratch and the branch is as it was. Clearing, merging and scanning only
up to the largest allele seen costs four instructions ending in a store
on every counted allele, to keep the record of that largest allele. At
one population the pass rose from 0.180 to 0.183 s and at four it fell
from 0.325 to 0.309 s: neither bar was met. Solving the two together
gives a saving of 0.0073 s per population and a penalty of 0.013 s per
pass, so even a free record would give only 0.007 s at one population.

At the 50 populations of 20 individuals that `docs/objectives.md` calls
ordinary, the same change takes the benchmark from 0.706 s to 0.448 s, a
gain of 37 in 100. That is where it belongs, and it is the owner's to
ask for. To be worth having there the penalty has to go, which means
taking the largest allele from a separate pass over the row rather than
from the counting loop; that was not built.

### Closed with no gain: linking across crates and one code generation unit

The release profile sets neither: it neither optimises across crate
boundaries nor compiles the crate as one piece. Neither paid. Both binaries and all
three built modules were kept and the three settings were run
interleaved, back to back, which removes the drift of the machine
entirely: on the benchmark the five statistics with one population read
0.182 and 0.183 s at the default, 0.183 and 0.184 s with the crate
compiled as one piece, and 0.183 and 0.184 s with thin linking; over
`big.vars` 0.290, 0.293 and 0.290 s. The read alone read 0.102 s in all
six runs, so linking across crates gives the lz4 decompression nothing,
although it lives in another crate and is exactly what that setting
should help. Fat linking was not run, because thin gave nothing.

Compiling the crate as one piece does inline two of the five counting
calls,
including one the review had singled out, and the pass does not get
faster, which closes that observation as a cause. It costs 2.5 times a
rebuild: 3.2 s against 1.3 s for the core crate and 5.7 s against 2.0 s
for the binding crate.

Two things for whoever runs the next one of these. `cargo asm` appends
its own single code generation unit to every invocation, so its listing
is the same whatever the release profile says and it cannot be used to
gate an experiment on linking or on how many pieces the crate is
compiled in; disassembling the binary that was actually built is what
sees them. And interleaving the builds back to
back is better than matching load averages: the drift of this machine
between two clean measurements is larger than any effect this experiment
was looking for.

## What is left, and what is the owner's

The target is not met and this review cannot meet it by itself. What
remains, in the order of what it would buy.

The read is 0.102 s of the 0.285 s and is serial. On 18 cores it is the
whole floor, under H3 above. `docs/architecture.md` already prescribes
the read ahead thread that would hide it and says it waits for a spec of
its own. It would not help the one thread number at all, and it would
change what that number means, since a native run would then own two
operating system threads. This is a decision and a spec, not an
experiment.

Of the 0.183 s that is not the read, the spec's own derivation assumed
about 0.112 to 0.120 s, being the pass plus twice what the missing data
filter adds. That derivation rested on the two walks of a row costing
the same, which this review disproved: they did not, by a factor of
three, and after the change they are much closer. Whether 0.25 s is
still the right number to ask for, now that the assumption behind it is
known to be wrong, is the owner's to say.

The prefix of the counters is worth 37 in 100 at 50 populations and
nothing at one, and needs its penalty removed first.

Every pass still decompresses 25 MB of an all ones mask and drops it,
11 in 100 of the bytes it decompresses, under L3. The sentence of
`docs/specs/io_vars.md` that says popnei does not pay for that mask is
wrong and should be corrected whatever is decided about the cost.

The chunk size is stated in rows, so a dataset of 10000 individuals uses
at most 8 of 18 cores, under L2. Nobody has run that dataset.
