# Which passes are worth a thread that reads one block ahead

25 September 2026. The association study of popnei reads one block while it
tests the block before it, which `docs/reports/perf-gwas-2026-09-24.md` built
and measured, and it was the only pass that did. The owner decided on 24
September 2026 that every other pass takes that reader where it makes the
pass significantly faster and does not where it does not, so each one is a
measurement and not an adoption. There are nine other passes over the blocks,
and this is what they read, what they work, what they took before the reader
and what they take after it.

Eight keep it. The kinship of 100000 variants of 1000 individuals falls from
0.300 s to 0.196 s on 18 cores, the principal components of the same
variants with the weights of ten components from 0.491 s to 0.355 s, the
Kosman distance of every pair of individuals from 0.219 s to 0.133 s, and
the diversity of three populations from 0.170 s to 0.118 s; the other four
are in section 4. Every number of every one of those passes is byte for byte
what it was, at one thread and at eighteen. The r² matrix does not keep it:
the whole of its pass over the blocks, the reading and the work together, is
0.010 s of a run of 0.452 s. What is left for the owner is in section 7, and
the largest of it is that a figure of popnei "on one thread" now means one
thread of arithmetic and one of reading, in eight passes and not one.

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
how long the pass spends inside `next_block` of its chain, which is the read
of the file and the decompression of a block, and how long it spends working
on the blocks the chain gave it. They are behind the cargo feature
`bench-phases` of the core crate, which no build popnei ships turns on, and
`crates/popnei/src/phases.rs` holds them. A **sampling profile**, which
stops a process a thousand times a second and records where it is, gives
neither: the chain decompresses on the thread that then computes, so the
samples of the read and the samples of the work are of one thread and one
stack. The most the reading thread can take off a pass is the **smaller of
the two clocks**, because what it hides behind one of them is the other.

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

Every pass of the core crate that reads blocks, at `ae6a190` of `main`:
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

Where the wall times come from. Five passes have a cargo benchmark that
reads the panel: `kinship.rs`, `pca_vars.rs`, `kosman_dists.rs`,
`pop_dists.rs` and `r2_matrix.rs`, all under `crates/popnei/benches/`, best
of 5 runs. The two passes of the stats module and the diversity have none:
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

## 2. The verdict: eight passes keep it, the r² matrix does not

On the branch `perf/read-ahead`, in seven commits, and **nothing is merged
into `main`**.

Each adoption was gated on an **identity check**: every number of every one
of the eight passes, 4103933 of them written at seventeen digits, comes out
byte for byte what it was, at one thread and at eighteen. The new
`crates/popnei/benches/the_numbers_of_every_pass.py` writes them and `cmp`
compares them. The same two files show that one thread and eighteen agree
with each other, before the change and after it.

Every check of the coding skill is clean at the head of the branch: `cargo
fmt`, clippy over the workspace with warnings denied, 1014 cargo tests at
the default features and the same 1014 with the `blas` feature off, which
runs every calculation on faer, the linear algebra a browser gets; 150 tests
of the linalg crate; both wasm targets and the JavaScript one compiling;
`ruff format` and `ruff check`; and 556 pytest tests.

## 3. Where the time of each pass goes

The two clocks, at 18 cores, over the panel. The last column is the smaller
of the two as a share of the wall time, which is the most the reading thread
could take off that pass.

| pass | inside the chain | working on the blocks | wall | ceiling |
|---|---|---|---|---|
| kinship | 0.114 s | 0.184 s | 0.299 s | 38 per 100 |
| principal components, one pass | 0.114 s | 0.178 s | 0.337 s | 34 per 100 |
| principal components, two passes | 0.227 s | 0.219 s | 0.490 s | 45 per 100 |
| Kosman distances | 0.116 s | 0.103 s | 0.219 s | 47 per 100 |
| distances between populations | 0.110 s | 0.023 s | 0.131 s | 18 per 100 |
| diversity of three populations | 0.113 s, carried | 0.052 s | 0.170 s | 31 per 100 |
| the five statistics of a variant | 0.113 s, carried | 0.013 s | 0.130 s | 10 per 100 |
| the two counts of an individual | 0.113 s, carried | 0.016 s, worked out | 0.129 s | 12 per 100 |
| the r² matrix | 0.006 s | 0.004 s | 0.452 s | 1 per 100 |

The read is 0.110 to 0.116 s in every pass but one, whatever the model and
whatever the thread count, because nothing of it is parallel: it is the lz4
decompression of the genotypes, 227 MB of output from 77.8 MB on disc, which
`docs/reports/perf-gwas-2026-09-24.md` counted from the footer of the file.
The same panel read from Python with `iter_blocks`, which hands a user the
blocks themselves and nothing else, takes 0.113 s: that is the same read with
a copy of each block into numpy on top of it.

The pass of the r² matrix is the one that reads almost nothing, and the
reason is its cap: the matrix holds one value for each pair of the variants
it was given, so the calculation refuses a pass longer than the cap the caller
gave it, and the measurement is of 5000 variants, one block, where every other
row of the table is 20 blocks. The 0.442 s the run has left over is the matrix
itself.

Three rows are not measured throughout, and say which part is not. The two
passes of the stats module and the diversity have no benchmark that reads a
file: their benchmarks build blocks in memory, so the read beside them,
0.113 s, is carried from what the other passes measured over this panel and
is not theirs. Their wall times in section 4 are over the panel and are
measured.

The work of those three. The diversity's 0.052 s and the five statistics'
0.013 s are the work clock over blocks in memory, each with the populations
and the statistics of the run whose wall time section 4 gives: three
populations and all five statistics at a draw of 200 for the first, and all
five statistics and no populations for the second, whose work is 0.028 s when
four populations are asked for instead. The two counts of an individual have
no benchmark at all, so their 0.016 s is not a clock: it is the wall time of
section 4 less the carried read.

## 4. What each adoption gave

Before and after, best of 5 runs of a cargo benchmark or best of 9 of a
timing script, over the panel. The last two columns are the fall at 18 cores
and at one thread.

| pass | 18 cores | one thread | 18 cores | one thread |
|---|---|---|---|---|
| kinship | 0.300 s to **0.196 s** | 0.426 s to **0.316 s** | 35 per 100 | 26 per 100 |
| principal components, no weights | 0.334 s to **0.233 s** | 0.460 s to **0.355 s** | 30 per 100 | 23 per 100 |
| principal components, ten weights | 0.491 s to **0.355 s** | 0.663 s to **0.477 s** | 28 per 100 | 28 per 100 |
| Kosman distances | 0.219 s to **0.133 s** | 0.883 s to **0.795 s** | 39 per 100 | 10 per 100 |
| distances between populations | 0.133 s to **0.112 s** | 0.281 s to **0.177 s** | 16 per 100 | 37 per 100 |
| the five statistics of a variant | 0.130 s to **0.115 s** | 0.215 s to **0.119 s** | 12 per 100 | 45 per 100 |
| the two counts of an individual | 0.129 s to **0.116 s** | 0.213 s to **0.119 s** | 10 per 100 | 44 per 100 |
| diversity of three populations | 0.170 s to **0.118 s** | 0.608 s to **0.508 s** | 31 per 100 | 16 per 100 |

The gate was a fall of more than 5 per cent of the pass's wall time and more
than the difference between the best and the worst of its 5 or 9 runs. The
smallest fall in the table, the two counts of an individual at 18 cores, is
0.013 s where the runs after the change read 0.116 s to 0.116 s and the runs
before it 0.129 s to 0.131 s: the two sets of runs do not touch. The
smallest fall in the table is twice the 5 per cent the gate asks for.

Two rows read against the ceiling of section 3 and are worth reading
together. At 18 cores the four passes that count genotypes and do no matrix
work — the distances between populations, the two of the stats module, and
the diversity — work for 0.013 to 0.052 s against a read of 0.113 s, so what
the thread hides is the whole of the work and the pass comes down to the
read: 0.112 s, 0.115 s, 0.116 s and 0.118 s against a read of 0.113 s. They
cannot go below that until the file is cheaper to read. At one thread the
same four work for 0.10 to 0.49 s, so there the thread hides the read
instead, which is why the fall at one thread is the larger one for three of
them.

The Kosman distance is the other way about: at one thread its work is 0.780 s
against a read of 0.110 s, so 0.110 s is all there is to hide and the fall
is 10 per cent; at 18 cores its work is 0.103 s and the two are level, which
is the largest share of any pass here.

The r² matrix was measured after the change as well, with the reader not
taken up, and gives 0.455 s against 0.452 s before, inside the spread of its
runs, 0.452 s to 0.477 s.

## 5. What it cost

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

## 6. Against the targets the specs state

Two of the eight passes have a number to reach in their spec, and both
numbers are of one thread.

**The kinship.** "Speed" of `docs/specs/kinship.md` names plink2's 0.23 s
over this panel on one thread, from section 2.1 of `docs/rust_core.md`.
popnei takes 0.316 s on one thread, 1.37 times that, where it took 0.426 s
and 1.85 times. On 18 cores it takes 0.196 s, which is below plink2's
one-thread figure and is not the same comparison. Nothing here ran plink2
again.

**The principal components.** "Speed" of `docs/specs/pca.md` names plink2
v2.0.0-a.7.7 `--pca 10 meanimpute --threads 1` at 0.248 s from its own file
of 2-bit genotypes, reading included, on a panel with 3 in 100 genotypes
missing. popnei's ten weights over the fully called panel take 0.477 s on
one thread against 0.663 s, so 1.92 times that figure where it was 2.67
times, over a different panel. plink2 was not run again and the two panels
are not the same: the spec's number is owed a run over
`bigcalled.vars` before it is compared closely.

The Kosman distance's spec states its target against a trial and pyNei and
not against another program, and the four counting passes state theirs
against pyNei; none of those was run again here, so this report does not move
them.

## 7. What is left for the owner to decide

1. **A figure of popnei "on one thread" now means one thread of arithmetic
   and one of reading, in eight passes.** `RAYON_NUM_THREADS=1` and
   `VECLIB_MAXIMUM_THREADS=1` no longer give a process that uses one core:
   the reading thread is an operating system thread of its own and it
   decompresses while the pass computes. The kinship's target in
   `docs/specs/kinship.md` is plink2 on one thread, and so is the principal
   components' in `docs/specs/pca.md`, so both comparisons now have popnei on
   two threads and plink2 on one. The gwas review raised this for one pass
   and it now holds for eight. Two ways to leave it: state the two threads
   wherever a one-thread figure is given, which is what section 6 does and
   what costs nothing; or give the eight passes a way to be asked for the
   serial path natively, a parameter or a cargo feature, which is a public
   choice in eight functions and which nothing has asked for. Recommendation:
   state the two threads, and revisit it if the owner wants a figure that a
   single core produced.
2. **Whether the branch is merged.** Seven commits, every check of the coding
   skill clean at the head, and every number of every pass byte for byte what
   it was at one thread and at eighteen.
3. **What the four counting passes are bounded by now.** The distances
   between populations, the two of the stats module and the diversity come
   down to the read at 18 cores, 0.112 to 0.118 s against a read of 0.113 s,
   so the next gain in any of them has to make the vars file cheaper to read
   and none of it is in those modules. That is the finding L7 of
   `docs/reports/perf-gwas-2026-09-24.md`, which measured the alternative:
   written with no compression the reader's lz4 work goes to zero and the file
   grows 2.9 times, and the Rust library popnei reads the vars file with reads
   an uncompressed one already. It now bounds five passes and not one.
4. **A figure of the gwas report was wrong by half and is corrected.**
   Section 2 and section 6.4 of `docs/reports/perf-gwas-2026-09-24.md` said
   that a block of 10000 individuals holds 250 variants and is 5 MB. The rule
   of section 5 above gives 500 variants there, and 500 times 10000 times the
   ploidy 2 is 10.0 MB, the same as the 5000 variants of 1000 individuals.
   Both places now carry 500 and 10.0 MB with a line that says what was
   corrected, and `docs/specs/block.md` states the rule with the two bounds.
   The cost of one more block is therefore 10.0 MB and not 5 MB at the largest
   dataset of the objectives.

## 8. The spec item that was owed

`docs/specs/block.md` said under "Not in this spec" that the read ahead
thread comes "with the first calculation that consumes blocks". Ten
calculations consume blocks and the reader is written, so nothing was left to
wait for. The item is written, as "The reader one block ahead", and it says: what the reader
gives a pass and what it leaves unchanged, the blocks in the order the chain
gave them, an error after every block that came before it and nothing after
it, the names of the chromosomes and the counts of the filters travelling
with each block, the chain lent and not given away so that a pass reads its
counts when it ends, and a change of which fields the pass asks for taking
effect one block later; what one costs, a thread and one block of memory
with the four sizes above; which passes take it, with the owner's rule of 25
September 2026 for deciding, a fall of more than 5 per cent of the pass's
wall time and more than the difference between the best and the worst of the
5 runs; how it runs, one reading thread, one channel that carries what the
thread read and one that carries a change of the fields the pass asks for,
and what happens when the chain fails and when the pass stops in the middle; and what the four cargo tests
of the reader assert, with the script that compares two commits number by
number.

`docs/specs/filters.md` had the same deferral and now points at the item.
`docs/specs/ld.md` records the measurement that closed it for the r² matrix.

## 9. What the next measurement would do

- **A cargo benchmark over a file for the stats module and the diversity.**
  Three of the eight wall times of section 4 are taken from Python because
  those modules have no benchmark that reads one. The boundary is about 10 ms
  of them and the falls are 0.013 to 0.096 s, so it changes nothing here, but
  a pass whose wall time is 0.115 s deserves a clock that is not behind an
  interpreter.
- **The eight passes at 10000 individuals.** The block there is 500 variants
  and the same 10.0 MB, and the work of a pass grows with the individuals
  while the read of a block does not, so the share the reader can hide is a
  different number and nobody has it.
- **The eight passes over a VCF.** Section 2 of
  `docs/reports/perf-gwas-2026-09-24.md` found a study over the 403 MB VCF
  faster than one over the 77.8 MB vars file, 0.103 s against 0.120 s,
  because the VCF reader parses its lines on the threads of rayon and the
  vars file reader decompresses lz4 on one. With a reading thread the vars
  file reader now has a thread of its own, which is the case where that
  comparison could turn round.
