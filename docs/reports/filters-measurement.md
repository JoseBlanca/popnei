# What the missing data filter costs a pass, task 4.1

21 September 2026. Task 4.1 of `docs/plans/filters.md`, the whole of its work
package 4, measures what the filter by the rate of missing genotypes costs a
pass over a 400 MB dataset, in popnei over a VCF and over a vars file, in
pyNei and in `bcftools view`. It sets no number to reach: "Speed" of
`docs/specs/filters.md` now carries these numbers and leaves the target to
the owner. It changed no code of the library; it added the benchmark
`crates/popnei/benches/filter_vars.rs`, the two scripts that timed pyNei and
bcftools beside it, and this report. The last section says what it leaves to
the performance review that follows this plan.

The words this document uses. A **pass** is one whole read of a dataset with
the genotypes alone asked for, from opening the file to the last block, which
is what a calculation over the genotypes asks of a reader; popnei's is timed
from building the reader to the last block, pyNei's from `vars_from_vcf` to
the last chunk, and bcftools' is the whole command. The **filter** keeps the
variants whose missing genotypes, divided by the individuals of the dataset,
are at most a threshold; a genotype is missing when one of its alleles at
least was not called. What the filter **costs** is the pass with it less the
pass without it, the two run back to back. A **block** is the array of
genotypes a reader of popnei gives at a time, and a **chunk** is pyNei's name
for the same thing. To **compact** a block is to copy the rows of the
variants the filter kept down over the ones it dropped and to cut the arrays
to what is left, which is what popnei does; pyNei builds a new chunk of the
variants it keeps instead. A **vars file** is popnei's own format, arrow IPC
with the genotypes compressed with lz4, which `docs/specs/io_vars.md`
describes; the variants of one of its blocks are written together as one
**batch** of that format. **The threads** of a timing are those of a pool of
rayon threads that the benchmark builds and reads inside: the VCF reader
parses its lines on it and the filter reads the rows of a block on it.

The two thresholds. At 0.1 the filter keeps every one of the 100000 variants,
whose genotypes are missing at a rate of 0.03, so at that threshold it counts
the missing genotypes of every variant and compacts nothing. At 0.03 it keeps
54773 and drops 45227, so the difference between the two thresholds is what
compacting the blocks costs: the copy moves the rows that stay, 54773 of them
and 109546000 alleles.

## The machine, the build and the files

The owner's Apple M5 Pro, 18 cores, 64 GB, macOS 27.0, native
`aarch64-apple-darwin`, rustc and cargo 1.98.0. popnei is built by `cargo
bench`, whose profile inherits from `[profile.release]`: optimized, with
`overflow-checks = false`. pyNei is at ef0ca6e, the commit `pyproject.toml`
names, under CPython 3.14.5 with numpy 2.5.3 and pandas 3.0.6. bcftools is
1.24 with htslib 1.24, at `/opt/homebrew/bin/bcftools`.

Nothing of popnei built while a timing ran and no two timings ran at once.
The machine was not idle: the owner's own processes kept the load average
between 1.0 and 3.7 through the session, and every set of runs below carries
the load average taken right before it. A load average taken right after a
set of 18 thread runs is the benchmark's own and says nothing about the
machine.

The files, both outside the repository in `/Users/jose/devel/popnei-bench/`.
`big.vcf`, 403572954 bytes, 100000 variants of 1000 individuals whose
genotypes are missing at a rate of 0.03, made by
`crates/popnei/benches/make_big_vcf.py`; every variant has `.` in its FILTER
column, so popnei's default of the variants that passed alone keeps all of
them. `big.vars`, 81356714 bytes, which popnei's `write_vars` wrote from that
VCF in Python, after `uv run maturin develop --release`, with the size of
block popnei chooses for 1000 individuals, 5000 variants, so the file holds
20 batches. Both files were in the page cache: every run of every benchmark
is preceded by one pass that is not timed.

## What the filter costs popnei

Every number is of 5 runs. Each pair, the pass with no filter and the pass
with the filter, was run back to back, and each pair was taken four times.
The two first columns are the first set, its median with its best and its
worst; the last column is what the filter costs in each of the four sets.

| | no filter | with the filter | what it costs, the four sets |
|---|---|---|---|
| VCF, 1 thread, at 0.1 | 0.561 s (0.557 to 0.563) | 0.630 s (0.625 to 0.640) | 0.069, 0.058, 0.070, 0.059 s |
| VCF, 1 thread, at 0.03 | 0.570 s (0.569 to 0.575) | 0.625 s (0.624 to 0.626) | 0.055, 0.051, 0.063, 0.085 s |
| VCF, 18 threads, at 0.1 | 0.095 s (0.092 to 0.099) | 0.102 s (0.102 to 0.105) | 0.007, 0.008, 0.008, 0.006 s |
| VCF, 18 threads, at 0.03 | 0.095 s (0.094 to 0.096) | 0.107 s (0.107 to 0.108) | 0.012, 0.011, 0.013, 0.011 s |
| vars file, 1 thread, at 0.1 | 0.102 s (0.101 to 0.103) | 0.158 s (0.157 to 0.158) | 0.056, 0.057, 0.056, 0.057 s |
| vars file, 1 thread, at 0.03 | 0.102 s (0.101 to 0.102) | 0.162 s (0.162 to 0.162) | 0.060, 0.060, 0.060, 0.060 s |
| vars file, 18 threads, at 0.1 | 0.102 s (0.102 to 0.103) | 0.115 s (0.115 to 0.115) | 0.013, 0.014, 0.014, 0.013 s |
| vars file, 18 threads, at 0.03 | 0.102 s (0.102 to 0.102) | 0.120 s (0.119 to 0.120) | 0.018, 0.018, 0.017, 0.018 s |

The fourth set of the VCF on one thread at 0.03 is the one `docs/plans/filters.md`
warns of when it says that a load on the machine which changes between the two
halves of a pair is larger than what is measured, so a pair whose halves start
at different load averages is run again: its two halves started at load
averages of 1.05 and 1.28, the only pair of the eight whose halves differ by
more than 0.15,
and it is the one that disagrees with its own kind. It was run again twice,
at load averages of 1.93 and 2.02 and then 1.93 and 1.93, and gave 0.049 s
and 0.078 s. Six sets of that pair spread from 0.049 s to 0.085 s.

Four things in that table are worth reading twice.

**Over the vars file on one thread the filter costs more than half of the
read**, 0.056 s against 0.102 s. Nothing is wrong with the filter: it reads
the same 1e8 genotypes over both files, 2e8 alleles at the ploidy of 2 of
this dataset, and costs 0.056 to 0.060 s on one thread whichever of the two
it reads them from, while the vars file is read in a fifth of the time the
VCF is parsed in. On 18 threads the same filter costs 0.013 to 0.018 s.

**The vars file is read in 0.102 s on 1 thread and on 18.** Its reader runs
on the thread that calls it, so over a vars file the threads are the filter's
alone, and the filter is what the 0.102 s pass and the 0.115 s pass differ
by. The VCF reader parses its lines on the threads, which is why its pass
falls from 0.561 s to 0.095 s.

**Compacting the blocks costs 0.003 to 0.005 s**, which is the threshold of
0.03 less the threshold of 0.1: 0.004 to 0.005 s over the VCF on 18 threads,
0.004 to 0.005 s over the vars file on 18, and 0.003 to 0.004 s over the vars
file on one. A review of the method of this measurement, of 21 September
2026, timed that copy by itself, the loop of `Block::retain_vars` that moves
the rows the filter kept down over the ones it dropped, over 20 blocks of
5000 variants of 2000 alleles with 54.8 in 100 of the rows kept: 3.3 to
5.0 ms for the whole dataset, which is the same figure reached another way.

**Over the VCF on one thread the two thresholds cannot be told apart.** The
sets there spread from 0.058 to 0.070 s at 0.1 and from 0.049 to 0.085 s at
0.03, while a set holds its own five runs inside 0.006 s: what moves is the
machine between one pair and the next, not the threshold. An earlier round of
this measurement, two sets instead of four, taken on a machine whose load
average was 1.9 to 2.1 instead of the 1.0 of these sets, read the gap between
the two thresholds there as 0.020 to 0.027 s, and the draft of this report
written from it said that gap was the compaction. It was not: the compaction
is the 0.004 s that the 18 thread sets and the vars file agree on, and where
the gap of 0.02 s came from is not known. No sampling profile was taken in
this plan, and that is the first thing for the performance review below.

## The variants kept

| threshold | popnei | pyNei | bcftools |
|---|---|---|---|
| 0.1 | 100000 | 100000 | 100000 |
| 0.03 | 54773 | 54773 | 54773 |

The three agree at both thresholds, over the VCF and, for popnei, over the
vars file as well. popnei's counts are the ones its benchmark prints, under
`missing_data`, the name this filter's counts carry for a user:
`the missing_data filter was given 100000 and kept 54773`; bcftools' are the
lines of `bcftools view -H -i "F_MISSING<=0.03" big.vcf | wc -l`; pyNei's are
the rows of the chunks its pass gave.

## pyNei

Taken as every other pair here, the pass with the filter and the pass without
it back to back, two sets of five runs in two processes, pyNei gives this:

| threshold | no filter | with the filter | what it costs |
|---|---|---|---|
| 0.1, first pair | 13.741 s (13.689 to 13.788) | 14.528 s (14.493 to 14.530) | 0.787 s |
| 0.1, second pair | 14.519 s (14.441 to 14.590) | 14.270 s (14.264 to 14.279) | -0.249 s |
| 0.03, first pair | 13.878 s (13.854 to 13.895) | 14.559 s (14.517 to 14.572) | 0.681 s |
| 0.03, second pair | 14.180 s (14.173 to 14.203) | 14.323 s (14.282 to 14.358) | 0.143 s |

The two pairs at 0.1 do not agree on the sign. A pass of pyNei takes 13.7 to
14.6 s, so five of them take a minute and a quarter, and over the ten minutes
these four pairs took, the time of the pass with no filter drifted from
13.741 s to 14.519 s while the load average of the machine fell from 2.2 to
1.1. The drift is larger than what the filter costs, and back to back pairs
of this length cannot tell the two apart.

So the two passes were run alternating inside one process, one of each in
turn, five of each, which spreads both over the same minutes:

| threshold | no filter | with the filter | what it costs |
|---|---|---|---|
| 0.1 | 14.050 s (14.027 to 14.053) | 14.464 s (14.453 to 14.481) | 0.415 s |
| 0.03 | 14.296 s (14.260 to 14.311) | 14.696 s (14.673 to 14.731) | 0.400 s |

Every one of the five runs with the filter was slower than every one of the
five without it, in both sets. The filter costs pyNei 0.40 to 0.42 s, 3 in
100 of its pass, and it costs the same at both thresholds: pyNei builds a new
chunk of the variants it keeps rather than compacting the one it was given,
so dropping 45227 variants of 100000 does not show.

## bcftools

`bcftools view -H`, with the records sent to `/dev/null` with `-o /dev/null`
in both passes of a pair, so that the destination costs the same in both.
Each pair was taken twice.

| threshold | no filter | with the filter | what it costs | again |
|---|---|---|---|---|
| 0.1 | 1.431 s (1.418 to 1.432) | 1.550 s (1.540 to 1.601) | 0.119 s | 1.445 and 1.568 s, 0.123 s |
| 0.03 | 1.436 s (1.423 to 1.469) | 1.234 s (1.217 to 1.251) | -0.202 s | 1.467 and 1.259 s, -0.208 s |

**At 0.03 bcftools is faster with the filter than without it**, by 0.20 s in
both sets. Its output is the difference: `view` writes the records it keeps
as text, 54773 of them instead of 100000, and the 45227 it does not write
save more than the filter costs. At 0.1, where it writes all 100000 either
way, the filter costs it 0.119 to 0.123 s.

## Three passes that are not the same work

The differences above are each of one program against itself, and those are
like for like. The whole passes are not, and none of the three should be read
as the same work done at a different speed:

- popnei's pass asks for the genotypes alone. It parses no other column of
  the VCF and decompresses no other column of the vars file, and it writes
  nothing.
- pyNei's pass builds, for every chunk, a pandas frame of the chromosome, the
  position, the id and the quality, and the alleles beside it, which
  `_parse_vcf_vars_chunk` of `pynei/io_vcf.py` does whatever the caller will
  read. Its 14.05 s is that work as well as the genotypes.
- bcftools' pass writes every record it keeps back out as text, on one
  thread. Its 1.431 s is a read and a write, and that write is why its
  difference at 0.03 is not the cost of a filter at all.

So the three columns of "The variants kept" are a comparison and the three
whole passes are not one.

## What this does not measure

The cost of the filter here is measured against a pass that reads no
genotype of its own: the benchmark walks the blocks and adds up the length of
their genotypes, and nothing after it looks at a genotype. A calculation that
reads the genotypes after the filter has them in the cache the filter left
them in, so the filter's share of that calculation's time is at most what is
measured here. How much less was not measured.

That the genotypes are there at all is what the count of alleles the
benchmark prints says: 200000000 in a pass that keeps every variant of these
two files, 109546000 in one that keeps 54773 of them. Without it, a reader
that stopped filling the column, or filled it only when somebody read it,
would make the pass with no filter shorter and the filter look dearer with
nothing to show it.

## Every set with the load average it was taken at

The load average is the one minute figure of `uptime`, read right before the
invocation. Each cell holds the four sets in order, one figure per set when
the two halves of the pair started at the same figure and both when they did
not.

| set | the load averages of the four sets |
|---|---|
| popnei VCF, 1 thread, 0.1 | 1.05, 1.04, 0.96, 0.96 |
| popnei VCF, 1 thread, 0.03 | 1.04, 1.03, 0.96, 1.05 and 1.28; run again at 1.93 and 2.02, and at 1.93 |
| popnei VCF, 18 threads, 0.1 | 1.37, 1.37, 1.37, 2.94 |
| popnei VCF, 18 threads, 0.03 | 2.94, 2.94, 2.94, 3.42 |
| popnei vars, 1 thread, 0.1 | 3.29, 3.29 and 3.18, 3.18, 3.18 |
| popnei vars, 1 thread, 0.03 | 3.18 and 3.09, 3.09, 3.09, 3.09 and 3.00 |
| popnei vars, 18 threads, 0.1 | 3.00, 3.00, 3.00, 3.72 |
| popnei vars, 18 threads, 0.03 | 3.72, 3.72, 3.72 and 3.58, 3.58 |
| bcftools, 0.1 | 2.57 and 2.48; 2.26 and 2.06 |
| bcftools, 0.03 | 2.36 and 2.30; 1.90 and 1.83 |
| pyNei, 0.1, the two pairs | 2.16 and 2.31; 1.08 and 1.23 |
| pyNei, 0.03, the two pairs | 2.58 and 1.53; 1.31 and 1.17 |
| pyNei alternating, 0.1 | 1.06 at the start, 1.10 at the end |
| pyNei alternating, 0.03 | 1.09 at the start, 1.29 at the end |

The first pair of pyNei at 0.03 started its two halves at 2.58 and at 1.53,
so it was run again, at 1.31 and 1.17, and that second pair is the one the
table of pyNei has. The repeat did not settle the question, and the
alternating runs are what did.

## The commands, as they were run

The vars file, from the worktree `.claude/worktrees/filters`:

    uv run maturin develop --release
    uv run python -c "import popnei; popnei.write_vars( \
        popnei.open_vcf('/Users/jose/devel/popnei-bench/big.vcf'), \
        '/Users/jose/devel/popnei-bench/big.vars')"

popnei, from the same worktree, for each file, each number of threads and
each threshold, with `uptime` before each of the two:

    cargo bench -q --bench filter_vars -- /Users/jose/devel/popnei-bench/big.vcf \
        --threads 1 --runs 5
    cargo bench -q --bench filter_vars -- /Users/jose/devel/popnei-bench/big.vcf \
        --threads 1 --runs 5 --max-missing-rate 0.1

bcftools, from `/Users/jose/devel/popnei-bench/`, through
`crates/popnei/benches/time_command.py`, which runs a command five times and
prints the best, the median and the worst, with one run before them that is
not timed:

    uv run --no-project python time_command.py 5 -- \
        bcftools view -H -o /dev/null big.vcf
    uv run --no-project python time_command.py 5 -- \
        bcftools view -H -i "F_MISSING<=0.1" -o /dev/null big.vcf

pyNei, through `crates/popnei/benches/time_pynei.py`, which times one pass
five times, or the two passes alternating in one process:

    uv run python time_pynei.py /Users/jose/devel/popnei-bench/big.vcf 5
    uv run python time_pynei.py /Users/jose/devel/popnei-bench/big.vcf 5 0.1
    uv run python time_pynei.py /Users/jose/devel/popnei-bench/big.vcf 5 0.1 --alternating

## What this hands to the performance review

Nothing here was made faster: what to do with these numbers is a performance
review, and this plan had none. What it should start from:

- **A sampling profile of the filtered pass**, on one thread and on 18, which
  is what would say where the 0.056 s over the vars file goes and where the
  gap of 0.02 s between the two thresholds over the VCF on one thread came
  from in the earlier round.
- **A threshold that keeps almost nothing**, which none of these sets has:
  every number here is of a filter that keeps all 100000 variants or 54773 of
  them, and what the filter costs when it keeps a hundredth of them is not
  measured.
- **The loop that compacts a block runs whole when every variant is kept**,
  so the 0.007 s of the filter at 0.1 on 18 threads holds a copy of every row
  of every block onto itself.
- **The filter builds a `Vec<bool>` for every block**, one value per variant,
  which is an allocation per block on the hot path.
- **Whether one reader that applies several thresholds in one pass is faster
  than several readers**, which "How it runs" of `docs/specs/filters.md`
  leaves to a measurement and which `docs/plans/filters.md` sends to the same
  review.

## Where the numbers that were not measured here come from

- The 1.65 s and 1.54 s of bcftools of "What has to be in place" of
  `docs/plans/filters.md` are one run each of 21 September 2026, with the
  records going to the terminal and not to `/dev/null`. The pair of that
  difference, 0.11 s, is the 0.119 to 0.123 s measured here.
- The VCF reader on its own, 0.563 s on one thread and 0.093 s on 18:
  `docs/reports/block-readers-measurement.md`, of the same file on the same
  machine. The passes with no filter here, 0.548 to 0.573 s and 0.095 to
  0.097 s, are that read.
- The 3.3 to 5.0 ms of the loop that compacts a block: the review of the
  method of this measurement, 21 September 2026, on 20 blocks of 5000
  variants of 2000 alleles with 54.8 in 100 of the rows kept.
