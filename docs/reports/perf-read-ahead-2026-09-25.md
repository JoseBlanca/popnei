# Which passes are worth a thread that reads one block ahead

25 September 2026. The association study of popnei reads one block while it
tests the block before it, which `docs/reports/perf-gwas-2026-09-24.md` built
and measured, and it was the only pass that did. The owner decided on 24
September 2026 that every other pass takes that reader where it makes the
pass significantly faster and does not where it does not, so each one is a
measurement and not an adoption. There are nine other passes over the blocks,
and this is what they read, what they work, what they took before the reader
and what they take after it.

Eight keep it, and the ninth, the r² matrix, does not. Every number of every
adopted pass is byte for byte what it was, at one thread and at eighteen.

Everything here is measured twice, against two different readers of popnei's
own file of variants, because that reader was parallelized on `main` while
this was being measured. Until `b9ad759` it decompressed lz4 on the thread
that called it, 0.112 s for this panel whatever the threads; since that
commit it decodes up to eight batches at once on the threads of the pool it
was called in, 0.027 s at 18 threads and 0.110 s at one, which
`docs/reports/perf-vars-threads-2026-09-25.md` measured. **The numbers this
report decides on are the second set**, on top of the parallel reader, since
that is what popnei has; section 4 gives them and section 5 keeps the first
set, which says what the read ahead thread was worth when the read was
serial.

What the parallel reader changed is how much there is left to gain at 18
cores, and not the verdict. The fall of a pass at 18 cores was 10 to 39 per
100 against the serial reader and is 6 to 16 per 100 against the parallel
one, because the read the thread hides went from 0.112 s to 0.027 s; at one
thread the read is still 0.110 s and the falls are 10 to 45 per 100. All
eight passes still clear the owner's gate at both thread counts, the kinship
by the least at 18 cores: 0.198 s to 0.184 s, a fall of 7 per 100 against a
gate of 5.

What is left for the owner is in section 8, and the largest of it is that a
figure of popnei "on one thread" now means one thread of arithmetic and one
of reading, in eight passes and not one.

## The words this document uses

A **block** is the run of consecutive variants a reader gives at a time,
5000 variants of 1000 individuals in every measurement here, so a pass over
the panel is 20 blocks. A **pass** is one reading of a source of variants
from its start to its end; a calculation makes one or two. A **chain of
readers** is what a pass reads from: the reader over the file, the filters a
user asked for on top of it, and `reblock`, which puts the blocks back to
one size.

The **reader one block ahead** is `with_one_block_ahead` of
`crates/popnei/src/block.rs`, which `docs/specs/block.md` now specifies. It
lends the chain of readers to a thread of its own for as long as one pass
runs: that thread builds the next block while the pass works on the one it
holds, and the pass reads its blocks off a **handle**, a reader like any
other. The handover is a rendezvous, so the thread holds at most one block
that is built and not yet given, and the memory of the pass grows by that
one block. In wasm there is no thread and the pass reads the chain itself,
because a browser has no threads.

The **two clocks** are what say beforehand what the reader can save a pass:
how long the pass spends inside `next_block` of its chain, and how long it
spends working on the blocks the chain gave it. They are behind the cargo
feature `bench-phases` of the core crate, which no build popnei ships turns
on, and `crates/popnei/src/phases.rs` holds them. A **sampling profile**,
which stops a process a thousand times a second and records where it is,
gives neither, because the two happen on one thread and one stack when
nothing puts them apart. The most the reading thread can take off a pass is
the **smaller of the two clocks**, since what it hides behind one of them is
the other. With the reading thread in place the first clock stops being the
read and becomes what the pass **waits** on the handle, which is what is left
of the read once the thread has it.

The **read of the panel** is two different things in this report, and each
number says which. Until `b9ad759` of `main` the reader of popnei's own file
decompressed lz4 on the thread that called it and took 0.112 s over this
panel whatever the threads. Since that commit it reads the bytes of up to
eight batches in the order of the file and decodes them on the threads of the
pool it was called in, and it takes 0.027 s at 18 threads and 0.110 s at one;
`docs/reports/perf-vars-threads-2026-09-25.md` is the review that did it.

The **panel** of every measurement is
`/Users/jose/devel/popnei-bench/bigcalled.vars`, 100000 variants x 1000
individuals with every genotype called, popnei's own vars file of 77.8 MB,
in the page cache. The **machine** is the owner's Apple M5 Pro, 18 cores, 64
GB, macOS 27, native `aarch64`. Every wall time is the best of 5 runs of a
cargo benchmark or the best of 9 runs of a timing script, after one run that
is not timed, and each says which. **18 cores** means `RAYON_NUM_THREADS`
and `VECLIB_MAXIMUM_THREADS` unset, so rayon and Accelerate take the
machine; **one thread** means both set to 1.

## 1. The scope, and where each number comes from

Every pass of the core crate that reads blocks. The branch started at
`ae6a190` of `main` and `main` was merged into it at `b9ad759`, which is
where the numbers of section 4 are taken; the passes are:
`calc_kinship` of `kinship.rs`, the two passes of `pca_of_variants` of
`pca.rs`, `calc_kosman_sums` of `dists.rs`, `calc_pop_dist_sums` of
`pop_dists.rs`, `the_dosages_of_the_pass` of `ld.rs`, which is the pass of
the r² matrix, `calc_pop_diversity` of `diversity.rs`, and
`calc_per_var_distribs` and `calc_per_individual_stats` of `stats.rs`. Nine
loops in eight functions. The association study is not measured again here.

Four of those functions take a reader whose size the compiler need not
know, because a binding crate hands them a chain of readers behind a pointer,
and `with_one_block_ahead` moves the reader it is given to its thread, which
it can only do with a value whose size it knows. Each of the four lends a
reference to its reader instead of the reader: a reference to a reader is a
reader of its own in popnei, and its size is the size of a pointer. Letting
`with_one_block_ahead` take a reader of unknown size directly does not
compile, because in wasm it hands that reader straight to the pass and the
compiler cannot make a pass's argument out of a value whose size it does not
know.

Where the wall times come from. Every before and after of section 4 is a
pair: the same command run against a binary built from `b9ad759`, which is
`main` with the parallel reader and without the read ahead thread, and
against one built from the head of this branch, the two alternating so that
whatever the machine was doing falls on both. Five passes have a cargo
benchmark that reads the panel: `kinship.rs`, `pca_vars.rs`,
`kosman_dists.rs`, `pop_dists.rs` and `r2_matrix.rs`, all under
`crates/popnei/benches/`, best of 5 runs, and 15 for the kinship, whose fall
is the smallest. The two passes of the stats module and the diversity have none:
their benchmarks build blocks in memory and read no file, so their wall
times are from the timing scripts `time_stats.py` and `time_diversity.py`
beside them, best of 9 runs, which run the pass from Python over the panel.
Those three numbers therefore carry the cost of the boundary with Python,
which `docs/reports/perf-gwas-2026-09-24.md` measured at about 10 ms for a
study; the passes themselves are Rust from the first block to the result and
nothing of a block crosses to Python. 9 runs and not 5 because 5 runs of the
first of them read 0.132 s to 0.328 s, a spread larger than anything this
report decides on, and 9 runs read 0.130 s to 0.133 s.

What is not measured. No wasm build was timed: there is no thread there and
the passes read their chain as they did. No dataset of 10000 individuals was
timed, so what the reader gives a pass whose blocks hold 500 variants
instead of 5000 is not known. The panel has every genotype called; a panel
with genotypes missing sends the kinship and the Kosman distance down a
second route in their work, which changes the work clock and not the read.
Nor was a VCF timed: the parallel reader of `b9ad759` is of popnei's own file
and the VCF reader already parsed its lines on the threads of rayon, so what
a reading thread adds there is a different question and nobody has asked it.

## 2. The verdict: eight passes keep it, the r² matrix does not

On the branch `perf/read-ahead`, and **nothing is merged into `main`**.

Each adoption was gated on an **identity check**: every number of every one
of the eight passes, 4103933 of them written at seventeen digits, comes out
byte for byte what `b9ad759` gives, at one thread and at eighteen. The new
`crates/popnei/benches/the_numbers_of_every_pass.py` writes them and `cmp`
compares them. The same files show that one thread and eighteen agree with
each other on both sides, and that `b9ad759` itself gives what `ae6a190`
gave, so the parallel reader moved no number either.

That check matters more after the merge than before it, because two changes
now touch how a block reaches a pass: the reader decodes several batches at
once and hands the blocks on in file order, and the reading thread of this
branch carries them across a channel. It ran again on the merge and held.

Every check of the coding skill is clean at the head of the branch: `cargo
fmt`, clippy over the workspace with warnings denied, 1017 cargo tests at
the default features and the same 1017 with the `blas` feature off, which
runs every calculation on faer, the linear algebra a browser gets; 150 tests
of the linalg crate; both wasm targets and the JavaScript one compiling;
`ruff format` and `ruff check`; and 556 pytest tests.

## 3. Where the time of each pass goes

The two clocks, at 18 cores, over the panel, on the head of this branch. The
first is what the pass waits on the handle, which is what is left of the read
once the reading thread has it: the read itself is 0.027 s, which the pass no
longer pays. The wall time beside them is of the same run, from a build with
the clocks compiled in, which is 1 to 3 ms longer than the build section 4
times; what a wall time is longer than the two clocks together is the work a
calculation does outside its loop over the blocks, the eigendecomposition of
the principal components and the matrix of the r² among it.

| pass | waits on the handle | works on the blocks | wall of that run |
|---|---|---|---|
| kinship | 0.008 s | 0.179 s | 0.187 s |
| principal components, no weights | 0.008 s | 0.170 s | 0.222 s |
| principal components, ten weights | 0.027 s | 0.211 s | 0.282 s |
| Kosman distances | 0.008 s | 0.103 s | 0.111 s |
| distances between populations | 0.008 s | 0.029 s | 0.037 s |
| the r² matrix | 0.008 s | 0.004 s | 0.448 s |

The two passes of the stats module and the diversity are not in the table:
their benchmarks build blocks in memory and read no file, so there is no read
for a clock to catch there. Their work over blocks in memory at 18 cores is
0.013 s for five statistics of a variant and 0.052 s for the diversity of
three populations, and the two counts of an individual have no benchmark at
all.

How much there is to hide is what `b9ad759` changed, and it is what makes the
falls of section 4 smaller at 18 cores than at one. Read alone from Python
with `iter_blocks`, which hands a user the blocks themselves and nothing else,
the panel is 0.027 s at 18 threads and 0.110 s at one, against 0.113 s at
either before that commit. So at 18 cores there is 0.027 s to hide and the
reading thread hides two thirds of it; at one thread there is 0.110 s, which
is why the fall at one thread is the larger one in seven of the eight passes.
The exception is the Kosman distance, whose work at one thread is 0.780 s, so
its 0.110 s of read is a small share of the run whatever is hidden.

The pass of the r² matrix is the one that reads almost nothing, and the
reason is its cap: the matrix holds one value for each pair of the variants
it was given, so the calculation refuses a pass longer than the cap the caller
gave it, and the measurement is of 5000 variants, one block, where every other
row of the table is 20 blocks. The 0.440 s the run has left over is the matrix
itself.

## 4. What each adoption gave

Before and after on top of the parallel reader: `b9ad759` against the head of
this branch, the same command alternating between the two binaries, best of 5
runs of a cargo benchmark or best of 9 of a timing script, over the panel. The
last two columns are the fall at 18 cores and at one thread.

| pass | 18 cores | one thread | 18 cores | one thread |
|---|---|---|---|---|
| kinship | 0.198 s to **0.184 s** | 0.400 s to **0.302 s** | 7 per 100 | 25 per 100 |
| principal components, no weights | 0.236 s to **0.221 s** | 0.433 s to **0.341 s** | 6 per 100 | 21 per 100 |
| principal components, ten weights | 0.300 s to **0.280 s** | 0.627 s to **0.455 s** | 7 per 100 | 27 per 100 |
| Kosman distances | 0.123 s to **0.110 s** | 0.852 s to **0.764 s** | 11 per 100 | 10 per 100 |
| distances between populations | 0.043 s to **0.036 s** | 0.269 s to **0.175 s** | 16 per 100 | 35 per 100 |
| the five statistics of a variant | 0.040 s to **0.034 s** | 0.213 s to **0.117 s** | 15 per 100 | 45 per 100 |
| the two counts of an individual | 0.041 s to **0.035 s** | 0.211 s to **0.116 s** | 15 per 100 | 45 per 100 |
| diversity of three populations | 0.079 s to **0.071 s** | 0.590 s to **0.490 s** | 10 per 100 | 17 per 100 |

The gate was a fall of more than 5 per cent of the pass's wall time and more
than the difference between the best and the worst of its runs. The kinship at
18 cores is the tightest and was measured with 15 runs a side twice over: the
best of the 15 is 0.198 s and 0.199 s without the thread and 0.184 s and
0.185 s with it, the 15 runs span 0.198 to 0.205 s and 0.184 to 0.193 s, and
the fall of 0.014 s is 7 per cent of the wall and larger than either spread.
The next tightest, the principal components with no weights, falls 0.015 s
where the runs without the thread span 0.002 s.

The r² matrix was measured after the change as well, with the reader not taken
up: 0.448 s, and its pass over the blocks reads 0.008 s and works 0.004 s, so
the most the thread could hide is 1 per cent of the run. The parallel reader
made its read cheaper still, so the refusal holds by more than it did.

**Taken again on 25 September 2026 after `main` was merged into this branch**,
which brought the exact residual of the association study at `294eacb`, on the
four passes that have a cargo benchmark reading the panel, best of 5 runs and
15 for the kinship, the two binaries alternating twice over:

| pass | 18 cores | one thread |
|---|---|---|
| kinship | 0.199 s to 0.184 s | 0.391 s to 0.298 s |
| principal components, no weights | 0.235 s to 0.221 s | 0.430 s to 0.335 s |
| distances between populations | 0.043 s to 0.036 s | 0.263 s to 0.174 s |
| Kosman distances | 0.124 s to 0.108 s | 0.851 s to 0.758 s |

The same benchmarks over blocks already in memory, which read no file and so
have no read for a thread to hide, are the control: the Kosman distance takes
0.101 s to 0.099 s at 18 cores and 0.755 s to 0.758 s at one, both inside the
spread of their runs.

The numbers of every one of the eight passes were written again at seventeen
digits on both sides of that merge, 4103933 of them over the panel at 18
threads, and the two files are identical byte for byte.

## 5. What the reader was worth when the read was serial

The first set of numbers, taken against `ae6a190`, where the reader of a vars
file decompressed lz4 on the thread that called it and a pass over this panel
read for 0.110 to 0.116 s whatever the threads. They are kept because they say
what the two changes are to each other: at 18 cores they are largely
substitutes, each hiding or removing the same serial read, and at one thread
only the reading thread helps, because the parallel reader has no threads to
decode on.

| pass | 18 cores | one thread |
|---|---|---|
| kinship | 0.300 s to 0.196 s | 0.426 s to 0.316 s |
| principal components, no weights | 0.334 s to 0.233 s | 0.460 s to 0.355 s |
| principal components, ten weights | 0.491 s to 0.355 s | 0.663 s to 0.477 s |
| Kosman distances | 0.219 s to 0.133 s | 0.883 s to 0.795 s |
| distances between populations | 0.133 s to 0.112 s | 0.281 s to 0.177 s |
| the five statistics of a variant | 0.130 s to 0.115 s | 0.215 s to 0.119 s |
| the two counts of an individual | 0.129 s to 0.116 s | 0.213 s to 0.119 s |
| diversity of three populations | 0.170 s to 0.118 s | 0.608 s to 0.508 s |

Read against the table of section 4, the Kosman distance at 18 cores is the
clearest case of the two changes doing the same work: the reading thread alone
took it from 0.219 s to 0.133 s, the parallel reader alone from 0.219 s to
0.123 s, and the two together give 0.110 s. The four passes that only count
genotypes are the same story at 18 cores, where the parallel reader is worth
more than the thread: the distances between populations went to 0.112 s with
the thread alone and to 0.043 s with the parallel reader alone.

The two clocks of that first set, at 18 cores, which are what chose the eight
passes: the read was 0.114 s for the kinship against 0.184 s of work, 0.114 s
against 0.178 s for the principal components with no weights, 0.227 s against
0.219 s for the two passes with ten weights, 0.116 s against 0.103 s for the
Kosman distance, 0.110 s against 0.023 s for the distances between
populations, and 0.006 s against 0.004 s for the r² matrix, which is what
refused it.

## 6. What it cost

One more block of memory in each of the eight passes. How large a block is
follows the rule of `docs/specs/block.md`: popnei divides 5 million genotypes
by the individuals of the source to get the variants of a block, and then
keeps that between 100 and 10000 variants. So a block holds those 5 million
genotypes wherever neither bound bites, which is 5000 variants of 1000
individuals and 500 of 10000, and at the ploidy 2 both are 10.0 MB of
alleles. Where a bound bites the block is another size: 50 individuals give
10000 variants and 1.0 MB, and 100000 individuals give 100 variants and 20.0
MB.

The resident set of a process was not measured again here. The call that
gives the largest resident memory a process has held reports it to two
decimals of a GB on this machine, and a kinship of this panel peaks at 0.23
GB, which cannot resolve 10 MB. The gwas review measured the resident set of
its own process grow by 12.8 MB with the same reader over the same panel.

One more thread in each of the eight passes, which in wasm is none: the
serial path is what a browser compiles, behind
`cfg(not(target_family = "wasm"))`, and no wasm build was timed.

In code, nine loops each wrapped in a call and a closure, one import, and
the four passes of section 1 lending a reference to their reader with a
comment that says why. Nothing in the error enum, so neither binding crate
changed.

The Kosman distance gave up one property. It turns the genotypes of a block
into sets of bits, one set for each allele of each individual, which is how
that calculation holds a block, and it dropped both the sets and the block
before asking the reader for the next one, so that the memory of two blocks
was never held at once. A reading thread that waits with a block it has built
and not yet given holds exactly that second block, so the property is gone.
The sets of bits are still dropped before the next block arrives, and they
are as large again as the block.

The two clocks cost nothing with the cargo feature off: the pass then calls
no clock at all and the timing wrapper is the work it was given. With the
feature on it is two readings of the system clock and two additions to a
counter per block, 40 of each over this panel.

## 7. Against the targets the specs state

Two of the eight passes have a number to reach in their spec, and both
numbers are of one thread.

**The kinship.** "Speed" of `docs/specs/kinship.md` names plink2's 0.23 s
over this panel on one thread, from section 2.1 of `docs/rust_core.md`.
popnei takes 0.302 s on one thread, 1.31 times that, where `b9ad759` takes
0.400 s and 1.74 times and `ae6a190` took 0.426 s and 1.85 times. On 18 cores
it takes 0.184 s, which is below plink2's one-thread figure and is not the
same comparison. Nothing here ran plink2 again.

**The principal components.** "Speed" of `docs/specs/pca.md` names plink2
v2.0.0-a.7.7 `--pca 10 meanimpute --threads 1` at 0.248 s from its own file
of 2-bit genotypes, reading included, on a panel with 3 in 100 genotypes
missing. popnei's ten weights over the fully called panel take 0.455 s on
one thread against 0.627 s at `b9ad759`, so 1.83 times that figure where it
was 2.53 times, over a different panel. plink2 was not run again and the two
panels are not the same: the spec's number is owed a run over
`bigcalled.vars` before it is compared closely.

The Kosman distance's spec states its target against a trial and pyNei and
not against another program, and the four counting passes state theirs
against pyNei; none of those was run again here, so this report does not move
them.

## 8. What is left for the owner to decide

1. **A figure of popnei "on one thread" now means one thread of arithmetic
   and one of reading, in eight passes.** `RAYON_NUM_THREADS=1` and
   `VECLIB_MAXIMUM_THREADS=1` no longer give a process that uses one core:
   the reading thread is an operating system thread of its own and it
   decompresses while the pass computes. The kinship's target in
   `docs/specs/kinship.md` is plink2 on one thread, and so is the principal
   components' in `docs/specs/pca.md`, so both comparisons now have popnei on
   two threads and plink2 on one. The gwas review raised this for one pass
   and it now holds for eight. Two ways to leave it: state the two threads
   wherever a one-thread figure is given, which is what section 7 does and
   what costs nothing; or give the eight passes a way to be asked for the
   serial path natively, a parameter or a cargo feature, which is a public
   choice in eight functions and which nothing has asked for. Recommendation:
   state the two threads, and revisit it if the owner wants a figure that a
   single core produced.
2. **Whether the branch is merged.** Every check of the coding skill is clean
   at its head over the merge with `b9ad759`, and every number of every pass
   is byte for byte what `b9ad759` gives, at one thread and at eighteen.
3. **The reader and the pass now use the same pool of threads, which nothing
   has reviewed.** Before `b9ad759` the reader decompressed on whatever thread
   called it and used no pool. It now decodes batches on the threads of rayon,
   and on this branch it does that from the reading thread while the pass runs
   its own rows on the same pool and, in the kinship and the principal
   components, calls the BLAS of the system from outside rayon at the same
   time. Every number is unchanged and every pass is faster, so nothing here
   says it is wrong; what it is is a combination neither review designed, and
   the concurrency category of a review has not looked at it. Recommendation:
   send one `perf-reviewer` for `concurrency` over the merge before it goes
   into `main`.
4. **What the four counting passes are bounded by now.** The distances between
   populations, the two of the stats module and the diversity take 0.034 to
   0.071 s at 18 cores where the read alone is 0.027 s, so between a third and
   three quarters of each of them is the file. The next gain in any of them is
   still mostly in the reader and not in those modules, although `b9ad759`
   took the read from 0.113 s to 0.027 s and the headroom left is much smaller
   than the finding L7 of `docs/reports/perf-gwas-2026-09-24.md` assumed:
   written with no compression the lz4 work goes to zero and the file grows
   2.9 times, and what that is worth now is 0.027 s at most and not 0.113 s.
5. **A figure of the gwas report was wrong by half and is corrected.**
   Section 2 and section 6.4 of `docs/reports/perf-gwas-2026-09-24.md` said
   that a block of 10000 individuals holds 250 variants and is 5 MB. The rule
   of section 6 above gives 500 variants there, and 500 times 10000 times the
   ploidy 2 is 10.0 MB, the same as the 5000 variants of 1000 individuals.
   Both places now carry 500 and 10.0 MB with a line that says what was
   corrected, and `docs/specs/block.md` states the rule with the two bounds.
   The cost of one more block is therefore 10.0 MB and not 5 MB at the largest
   dataset of the objectives.

## 9. The spec item that was owed

`docs/specs/block.md` said under "Not in this spec" that the read ahead
thread comes "with the first calculation that consumes blocks". Ten
calculations consume blocks and the reader is written, so nothing was left to
wait for. The item is written, as "The reader one block ahead", and it says: what the
reader gives a pass and what it leaves unchanged, the blocks in the order the chain
gave them, an error after every block that came before it and nothing after
it, the names of the chromosomes and the counts of the filters travelling
with each block, the chain lent and not given away so that a pass reads its
counts when it ends, and a change of which fields the pass asks for taking
effect one block later; what one costs, a thread and one block of memory
with the four sizes above; which passes take it, with the owner's rule of 25
September 2026 for deciding, a fall of more than 5 per cent of the pass's
wall time and more than the difference between the best and the worst of its
runs; how it runs, one reading thread, one channel that carries what the
thread read and one that carries a change of the fields the pass asks for,
and what happens when the chain fails and when the pass stops in the middle;
and what the four cargo tests of the reader assert, with the script that
compares two commits number by number.

`docs/specs/filters.md` had the same deferral and now points at the item.
`docs/specs/ld.md` records the measurement that closed it for the r² matrix.

## 10. What the next measurement would do

- **A cargo benchmark over a file for the stats module and the diversity.**
  Three of the eight wall times of section 4 are taken from Python because
  those modules have no benchmark that reads one. The boundary is about 10 ms
  of them, and two of those three passes now take 0.034 s and 0.035 s at 18
  cores, so the interpreter is a third of what is being timed. Their falls,
  0.006 s each, are larger than that and hold; the next pass that wants a
  number of its own at 18 cores needs the clock out from behind Python.
- **The eight passes at 10000 individuals.** The block there is 500 variants
  and the same 10.0 MB, and the work of a pass grows with the individuals
  while the read of a block does not, so the share the reader can hide is a
  different number and nobody has it.
- **The eight passes over a VCF.** Every number here is of popnei's own file.
  The VCF reader already parses its lines on the threads of rayon, so a
  reading thread in front of it overlaps something that is already parallel,
  which is the case this report does not cover and which the four counting
  passes would show most clearly.
