# How fast the Kosman distances of every pair are, task 3.1

22 September 2026. Task 3.1 of `docs/plans/dists-kosman.md`, the first
deliverable of its work package 3, measures how long popnei takes to give the
Kosman distance of every pair of 1000 individuals over 100000 variants:
natively on one thread and on 18 cores, over a VCF and over a vars file, in
WebAssembly under node, and pyNei on the same dataset beside it. "Speed" of
`docs/specs/dists.md` sets the three numbers to reach, and none of the three
is reached: 1.154 s against 0.97 s on one thread, 0.625 s against 0.38 s on
18 cores and 2.130 s against 1.43 s in wasm. What to do about that is the
owner's, and the last section says what it leaves to the performance review
the spec names. It changed no code of the library; it added the two scripts
that timed popnei and pyNei in Python, the one that timed the wasm build
under node, and this report.

The words this document uses. A **pass** is one whole read of a dataset with
the genotypes alone asked for, from opening the source to the last block,
which is what this calculation asks of a reader. A **block** is the array of
genotypes a reader gives at a time, 5000 variants of these 1000 individuals,
so the dataset is 20 blocks. The **calculation** is one call of
`calc_pairwise_kosman_dists` over a source, which makes such a pass and gives
the distance of every pair; it is also called the whole call below. The
**reading alone** is the same pass with the genotypes as the only field the
blocks carry and nothing done to them but adding up how many alleles each
block says it holds. Natively that count is the shape of the array the core
filled, `block.gts.size`, which reaches numpy without a copy: it shows a pass
whose blocks carry no genotype column at all, and it does not show a buffer
of the right shape that nothing wrote into, since no genotype is read. In
wasm it is the length of the `Int8Array` that the package copies out of the
memory of the WebAssembly for each block, so there the genotypes are moved
whether or not anybody reads them. **With the reading taken out** is the first less the second, and it is
what the three numbers of "Speed" are of, so that the reader is not in them;
each of those three is a **target** here. A **vars file** is popnei's own
format, one arrow IPC file with the genotypes compressed with lz4, which
`docs/specs/io_vars.md` describes. **The threads** are those of the pool of
rayon, the Rust library popnei runs its parallel loops on, which the library
never builds itself: the pool takes one thread for each core unless
`RAYON_NUM_THREADS` is in the environment of the process before the first
call. The **load average** is the one minute figure of `sysctl -n
vm.loadavg`, read right before the invocation and again right after it.

## The machine, the build, the versions and the files

The owner's Apple M5 Pro, 18 cores, 64 GB, macOS 27.0, native
`aarch64-apple-darwin`, rustc 1.98.0, node 26.8.2, CPython 3.14.5 with numpy
2.5.3 on Accelerate, the BLAS of macOS, and pandas 3.0.6. pyNei is at
ef0ca6e, the commit `pyproject.toml` names. popnei is the branch
`plan/dists-kosman` at 3ec7151, built in release both ways: `uv run maturin
develop --release` for the Python package, and, in `js/popnei`, `npm run
build`, whose `cargo build --package popnei-js --release --target
wasm32-unknown-unknown` is the WebAssembly the node runs. The machine was not
idle: the owner's own processes kept the load average between 1.42 and 1.93
through the session, and every table below carries the figure before and
after its runs.

The dataset is one file in three forms, all outside the repository in
`/Users/jose/devel/popnei-bench/`, 100000 variants of 1000 diploid
individuals, biallelic, with 3 in 100 genotypes missing:

- `big.vcf`, 403572954 bytes, made by
  `crates/popnei/benches/make_big_vcf.py`.
- `big.vars`, 81356714 bytes, written from it by popnei's `write_vars` in
  batches of 5000 variants, the size popnei chooses for 1000 individuals.
- `big.pynei.vars`, 36454074 bytes, written from the same VCF by pyNei's own
  `write_vars` on 22 September 2026, because pyNei does not read a vars file
  of popnei: both formats are arrow IPC and their columns differ. It is
  smaller because pyNei compresses its buffers with zstd where popnei uses
  lz4. The command was one line:

      uv run python -c "from pynei import vars_from_vcf; \
          from pynei.io_vars import write_vars; \
          write_vars(vars_from_vcf('/Users/jose/devel/popnei-bench/big.vcf'), \
          '/Users/jose/devel/popnei-bench/big.pynei.vars')"

  It took 14.8 s, and pyNei's default chunk of 5 million genotypes is 5000
  variants here, so its file holds the same 20 batches as popnei's.

Every file was in the page cache: each script runs one pass of each kind that
is not timed before the timed ones.

One thing to know before running the checks of the plan against a release
build. `uv run pytest` gives 208 passed against the build `uv run maturin
develop` makes. Against the release build one of those 208 fails:
`test_a_ctrl_c_while_write_vars_runs_is_raised_and_leaves_no_file` of
`tests/test_io_vars.py` sends itself the signal of a Ctrl-C 0.1 s into a
`write_vars` that takes 0.313 s in the debug build, and in the release build
the call is over before the signal arrives. The measurement here needs the
release build and the test suite needs the other, so the two are run one
after the other.

## Natively over the vars file

Each figure is the best of 3 runs, with the spread of the three beside it.
The calculation and the reading alone were run one after the other inside
each round, so that a machine that grew busier over the runs would fall on
both. `open_vars` is outside both: it reads the schema and the footer of the
file and no genotype, every pass reads the file again from its start, and it
took 0.2 ms.

| threads | the whole call | the reading alone | with the reading taken out | the load average |
|---|---|---|---|---|
| 1 | 1.259 s (1.259 to 1.264) | 0.105 s (0.105 to 0.105) | 1.154 s | 1.75 before, 1.93 after |
| 18 | 0.735 s (0.735 to 0.751) | 0.110 s (0.110 to 0.110) | 0.625 s | 1.81 before and after |

The reading alone is the same on one thread and on 18, as it was in
`docs/reports/filters-measurement.md`: the reader of a vars file runs on the
thread that calls it, so over this file the threads are the calculation's
alone. Section 1 of `docs/architecture.md` has arrow-rs, the Rust
implementation of arrow that reads a vars file, decompressing the genotypes
of one of 20000 variants of 1000 individuals in 18.8 ms on one thread of this
machine; five times that file is 94 ms, and the whole pass
here, decompression, the reader and the loop in Python over 20 blocks, is
105 ms.

## Natively over the VCF

The same, through `open_vcf` on `big.vcf`. `open_vcf` reads the header alone
and took 0.1 ms, outside both timings.

| threads | the whole call | the reading alone | with the reading taken out | the load average |
|---|---|---|---|---|
| 1 | 1.781 s (1.781 to 1.784) | 0.603 s (0.603 to 0.612) | 1.178 s | 1.60 before, 1.50 after |
| 18 | 0.718 s (0.718 to 0.731) | 0.098 s (0.098 to 0.099) | 0.620 s | 1.70 before, 1.65 after |

What is left when the reading is taken out is the same over both files,
1.154 and 1.178 s on one thread and 0.625 and 0.620 s on 18, which is what
says that the subtraction is of the reader and not of something else: the VCF
reader parses its lines on the threads and the vars file reader does not, so
the two readings differ by 0.5 s on one thread and by nothing on 18, and the
calculation over them does not.

## In WebAssembly under node

One thread, since wasm has none. Reading `big.vars` into a `Uint8Array`,
0.015 s, and `openVars`, which copies those bytes into the memory of the
WebAssembly and reads the schema and the footer, 0.008 s, are outside both
timings.

| threads | the whole call | the reading alone | with the reading taken out | the load average |
|---|---|---|---|---|
| 1 | 2.291 s (2.291 to 2.321) | 0.161 s (0.161 to 0.175) | 2.130 s | 1.54 before, 1.68 after |

The 2.130 s is the low end of what the calculation costs beyond its reader
there. In wasm every block is copied out of the memory of the WebAssembly
into an `Int8Array` of JavaScript before the loop of the user sees it,
2e8 bytes over the 20 blocks, and the calculation makes no such copy: the
reading alone carries it and the whole call does not, so the difference of
the two is smaller than the calculation's own share. Natively there is no
such copy, because the array the core filled reaches numpy as it is.

That copy is small. The review of this task timed it inside the built
`js/popnei/dist/block.js` over this file, where it is 5 to 10 ms of the
0.161 s, and a plain copy of 200 MB in node with nothing else going on is 4
to 9 ms. So of the 56 ms between the 0.161 s in wasm and the 0.105 s natively
on one thread, the copy is a tenth or less and the rest is the reader itself
being slower in wasm, which nothing here measured on its own. The 2.130 s is
a low end by that unmeasured amount and not by the copy.

## pyNei

One call of pyNei's `calc_pairwise_kosman_dists` over the `Variants` that its
`load_vars` gives for `big.pynei.vars`, best of 3 runs, run on 22 September
2026. `load_vars` is outside the timing, as `open_vars` is outside popnei's:
it memory maps the file and reads its schema and its metadata, and the chunks
of genotypes are read inside the call. A fresh `Variants` is built for every
run, so no chunk is read twice.

| num_threads | the whole call | the load average |
|---|---|---|
| 1 | 0.719 s (0.719 to 0.750) | 1.42 before, 1.63 after |
| 6 | 0.280 s (0.280 to 0.294) | 1.50 before and after |

No reading is taken out of these. pyNei's call and popnei's are not the same
work in two more ways: pyNei reads a chunk through `VariantsFile.read_chunk`,
which builds a pandas frame of the chromosome, the position, the id and the
quality of every chunk whatever the caller will read, and it computes the
counts of every pair as the products of matrices of 0 and 1 in float32, which
numpy makes with Accelerate on the matrix units of the Apple chip. The 6
threads are its best, as "Speed" of the spec has it, and the same spec gives
0.76 s and 0.30 s for the same calculation over a dataset in memory on
21 September 2026; from its own file it is 0.719 s and 0.280 s.

## The three numbers of "Speed"

Each target is of the calculation with the reading taken out, and the figures
are those of the tables above, over the vars file natively and over the bytes
of that file in wasm.

| the number to reach | the whole call | the reading alone | with the reading taken out | met |
|---|---|---|---|---|
| 0.97 s on one thread | 1.259 s | 0.105 s | 1.154 s | no, 1.19 times the target |
| 0.38 s on 18 cores | 0.735 s | 0.110 s | 0.625 s | no, 1.64 times the target |
| 1.43 s in wasm | 2.291 s | 0.161 s | 2.130 s | no, 1.49 times the target |

Each of the three targets is what the trial of "How it runs" of the spec took
with a tenth of that time added to it, and that trial timed one block of 5000 variants of 1000
individuals at 0.044 s on one thread, 0.017 s on 18 cores and 0.065 s in
wasm, from which 20 blocks are 0.88 s, 0.34 s and 1.30 s. The calculation
measured here is 0.058 s a block on one thread, 0.031 s on 18 cores and
0.107 s in wasm: 1.31, 1.84 and 1.64 times the trial's block. So what the
implementation costs over the trial is largest on 18 cores, where it goes
from one thread to 18 by a factor of 1.85, 1.154 s to 0.625 s, and the trial
went by a factor of 2.6, 0.044 s to 0.017 s a block. No sampling profile was
taken in this task, so where that goes is not known.

Beside pyNei on the same dataset and the same machine, in two ways, because
pyNei's number has no reading taken out of it and popnei's targets do. With
the reading taken out of popnei's alone, popnei's 0.625 s on 18 cores is 2.2
times pyNei's 0.280 s on 6 threads. As a user waits, whole call against whole
call, it is 0.735 s against 0.280 s, 2.6 times, and 1.259 s against 0.719 s
on one thread.

The 1.3 times of "Speed" of the spec is neither of those two. It is the
targets themselves over pyNei's times on a dataset held in memory, 0.97 over
0.76 on one thread and 0.38 over 0.30 on its 6 threads, so it says how much
slower than pyNei the spec was willing to be and not what was measured here.
That is what the owner decided on 21 September 2026 when he chose the sets of
bits over the matrix products pyNei makes: the bits took 1.5 times
Accelerate's time on the trial's first block and would take about 1.3 times
pyNei's over the whole dataset, and what they gave for it is integer sums
that do not change with the threads, no need of the `linalg` module, which
the architecture puts after `dists` and which is not written, and a wasm
build where there is no BLAS at all.

## What a read ahead thread could gain at most

The read ahead thread of section 3 of `docs/architecture.md` reads the next
block while the calculation works on the one in hand.
`docs/plans/dists-kosman.md` leaves it out of this plan, because no spec
describes it, and the reading alone above is what it could gain at most: a
thread of its own can take the reading off the thread that calculates, and it
cannot take off more than the reading costs.

Over the vars file that is 0.110 s of the 0.735 s whole call on 18 cores, 15
in 100 of it, and 0.105 s of the 1.259 s on one thread, 8 in 100. Over the
VCF it is 0.098 s of 0.718 s on 18 cores and 0.603 s of 1.781 s on one
thread, a third of that call. On one thread the gain is not free in the way
the other is: a read ahead thread is a second thread, so the calculation that
ran on one core would run on two.

None of that would reach a target. The targets have the reading taken out
already, so a read ahead thread changes what a user waits and changes none of
the three lines of the table above.

## The commands, as they were run

From the worktree `.claude/worktrees/dists-kosman`, after `uv run maturin
develop --release`, with `sysctl -n vm.loadavg` before and after each:

    uv run python crates/popnei/benches/time_kosman_dists.py \
        /Users/jose/devel/popnei-bench/big.vars 3
    RAYON_NUM_THREADS=1 uv run python crates/popnei/benches/time_kosman_dists.py \
        /Users/jose/devel/popnei-bench/big.vars 3
    uv run python crates/popnei/benches/time_kosman_dists.py \
        /Users/jose/devel/popnei-bench/big.vcf 3
    RAYON_NUM_THREADS=1 uv run python crates/popnei/benches/time_kosman_dists.py \
        /Users/jose/devel/popnei-bench/big.vcf 3
    uv run python crates/popnei/benches/time_kosman_pynei.py \
        /Users/jose/devel/popnei-bench/big.pynei.vars 3 1
    uv run python crates/popnei/benches/time_kosman_pynei.py \
        /Users/jose/devel/popnei-bench/big.pynei.vars 3 6

`RAYON_NUM_THREADS` is set in the command and not from Python because rayon
reads it when its pool is first used, which is inside the first call.

From `js/popnei`, after `npm run build`:

    node bench/time_kosman_dists.mjs /Users/jose/devel/popnei-bench/big.vars 3

Each script prints every run, the best, the median and the worst of each kind
of pass, and the difference of the two on the bests. The three are not part
of the tests: `npm test` runs `node --test test/` and does not look in
`bench/`, and ruff's `include` of `pyproject.toml` is the package and the
tests.

## What this hands to the performance review

Nothing here was made faster, and no target was met. What a performance
review of this calculation should start from, the first two being what
"Speed" of `docs/specs/dists.md` already names:

- **Fewer sets for a block with two alleles.** A block gets 1 + k * A sets of
  bits per individual, with k the ploidy and A the alleles the block holds,
  which for this biallelic diploid dataset is 1 + 2 * 2, 5 sets, and the
  pairs cost in proportion. The spec
  leaves the layout to the implementer and asks whether a block with only the
  alleles 0 and 1 needs all 5. Not measured here.
- **A parallel building of the sets.** The trial of the spec built them
  serially, 0.015 s of its 0.044 s block, and spread only the pairs over the
  threads. What the implementation does with them, and what it would gain
  from building them on the threads, is not measured here, and it is where
  the factor of 1.85 from one thread to 18, against the trial's 2.6, would be
  looked for first.
- **A sampling profile of the calculation**, on one thread and on 18, which
  none of this task took, and which is what would say where the 1.154 s and
  the 0.625 s go and why each block costs 1.31 and 1.84 times the trial's.
- **The split of the work over the threads.** 499500 pairs of 1000
  individuals are the work items of a block, and
  `docs/plans/dists-kosman.md` leaves the split to the implementer with the
  timing here as what says whether it was right. Nothing here varied it.
- **What the reader costs in wasm.** Reading this file alone takes 0.161 s
  there against 0.105 s natively on one thread, and the copy of every block
  out of the memory of the WebAssembly, which the reading alone carries and
  the calculation does not, is 5 to 10 ms of it. So most of those 56 ms is
  the reader itself, which nothing here timed on its own, and the 2.130 s of
  the table is a low end by that much.

What was not measured at all: any ploidy but 2, any dataset but this one,
popnei under pyodide, and pyNei under pyodide.
