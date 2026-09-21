# What the missing data filter costs a pass, task 4.1

21 September 2026. Task 4.1 of `docs/plans/filters.md`, the whole of its work
package 4, measures what the filter by the rate of missing genotypes costs a
pass over a 400 MB dataset, in popnei over a VCF and over a vars file, in
pyNei and in `bcftools view`. It sets no number to reach: "Speed" of
`docs/specs/filters.md` now carries these numbers and leaves the target to
the owner. It changed no code of the library; it added the benchmark
`crates/popnei/benches/filter_vars.rs` and this report.

The words this document uses. A **pass** is one whole read of a dataset with
the genotypes alone asked for, from opening the file to the last block, which
is what a calculation over the genotypes does; popnei's is timed from
building the reader to the last block, pyNei's from `vars_from_vcf` to the
last chunk, and bcftools' is the whole command. The **filter** keeps the
variants whose missing genotypes, divided by the individuals of the dataset,
are at most a threshold; a genotype is missing when one of its alleles at
least was not called. What the filter **costs** is the pass with it less the
pass without it, the two run back to back. A **block** is the array of
genotypes a reader of popnei gives at a time, and a **chunk** is pyNei's
name for the same thing. To **compact** a block is to take the variants the
filter dropped out of it, moving the genotypes of the ones it kept up over
them, which is what popnei does and pyNei does not: pyNei builds a new chunk
of the variants it keeps. A **vars file** is popnei's own format, arrow IPC
with the genotypes compressed with lz4, which `docs/specs/io_vars.md`
describes; the variants of one of its blocks are written together as one
**batch** of that format. **The threads** of a timing are those of a pool of
rayon threads that the benchmark builds and reads inside: the VCF reader
parses its lines on it and the filter reads the rows of a block on it.

The two thresholds. At 0.1 the filter keeps every one of the 100000 variants,
whose genotypes are missing at a rate of 0.03, so at that threshold it counts
the missing genotypes of every variant and compacts nothing. At 0.03 it keeps
54773, and the 45227 it drops are taken out of the block, so the difference
between the two thresholds is what compacting the blocks costs.

## The machine, the build and the files

The owner's Apple M5 Pro, 18 cores, 64 GB, macOS 27.0, native
`aarch64-apple-darwin`, rustc and cargo 1.98.0. popnei is built by `cargo
bench`, whose profile inherits from `[profile.release]`: optimized, with
`overflow-checks = false`. pyNei is at ef0ca6e, the commit `pyproject.toml`
names, under CPython 3.14.5 with numpy 2.5.3 and pandas 3.0.6. bcftools is
1.24 with htslib 1.24, at `/opt/homebrew/bin/bcftools`.

Nothing of popnei built while a timing ran and no two timings ran at once.
The machine was not idle: the owner's own processes kept the load average
between 1.0 and 3.5 through the session, and every set of runs below carries
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
20 batches. Both files were in the page cache: every run of every benchmark is
preceded by one pass that is not timed.

## What the filter costs popnei

Each number is the median of 5 runs, and each pair, with the filter and
without, was taken twice; the last column holds the second set, its two
medians and what the filter costs in it.

| | no filter | with the filter | what it costs | again |
|---|---|---|---|---|
| VCF, 1 thread, at 0.1 | 0.565 s | 0.615 s | 0.050 s | 0.577 and 0.630 s, 0.053 s |
| VCF, 1 thread, at 0.03 | 0.553 s | 0.630 s | 0.077 s | 0.566 and 0.635 s, 0.069 s |
| VCF, 18 threads, at 0.1 | 0.095 s | 0.102 s | 0.007 s | 0.093 and 0.101 s, 0.008 s |
| VCF, 18 threads, at 0.03 | 0.093 s | 0.105 s | 0.012 s | 0.094 and 0.106 s, 0.012 s |
| vars file, 1 thread, at 0.1 | 0.101 s | 0.157 s | 0.056 s | 0.101 and 0.157 s, 0.056 s |
| vars file, 1 thread, at 0.03 | 0.101 s | 0.160 s | 0.059 s | 0.101 and 0.161 s, 0.060 s |
| vars file, 18 threads, at 0.1 | 0.101 s | 0.115 s | 0.014 s | 0.101 and 0.114 s, 0.013 s |
| vars file, 18 threads, at 0.03 | 0.101 s | 0.118 s | 0.017 s | 0.101 and 0.119 s, 0.018 s |

Three things in that table are worth reading twice.

**Over the vars file on one thread the filter costs more than half of the
read**, 0.056 s against 0.101 s. Nothing is wrong with the filter: it reads
the same 1e8 genotypes over both files, 2e8 alleles at the ploidy of 2 of
this dataset, and costs 0.050 to 0.059 s on one
thread whichever of the two it reads them from, and the vars file is read in
a fifth of the time the VCF is parsed in. On 18 threads the same filter costs
0.013 to 0.018 s.

**The vars file is read in 0.101 s on 1 thread and on 18.** Its reader runs
on the thread that calls it, so over a vars file the threads are the filter's
alone, and the filter is what the 0.101 s pass and the 0.115 s pass differ
by. The VCF reader parses its lines on the threads, which is why its pass
falls from 0.565 s to 0.095 s.

**The threshold of 0.03 costs about 0.02 s more than the one of 0.1 over the
VCF on one thread**, 0.069 to 0.077 s against 0.050 to 0.053 s, and 0.004 to
0.005 s more on 18 threads. That is the compaction of 45227 variants out of
the blocks, a move of 4.5e7 genotypes and so of 9.0e7 alleles. Over the vars
file the same compaction
shows as 0.003 to 0.004 s on both thread counts; why it is larger over the
VCF than over the vars file, where the blocks and the variants dropped are
the same, is not known and was not looked into: it is a performance review
and not this plan.

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
| 0.1, first pair | 13.741 s | 14.528 s | 0.787 s |
| 0.1, second pair | 14.519 s | 14.270 s | -0.249 s |
| 0.03, first pair | 13.878 s | 14.559 s | 0.681 s |
| 0.03, second pair | 14.180 s | 14.323 s | 0.143 s |

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
| 0.1 | 14.050 s | 14.464 s | 0.415 s |
| 0.03 | 14.296 s | 14.696 s | 0.400 s |

Every one of the five runs with the filter was slower than every one of the
five without it, in both sets: 14.027 to 14.053 s against 14.453 to 14.481 s
at 0.1, and 14.260 to 14.311 s against 14.673 to 14.731 s at 0.03. The filter
costs pyNei 0.40 to 0.42 s, 3 in 100 of its pass, and it costs the same at
both thresholds: pyNei builds a new chunk of the variants it keeps rather
than compacting the one it was given, so dropping 45227 variants of 100000
does not show.

## bcftools

`bcftools view -H`, with the records sent to `/dev/null` with `-o /dev/null`
in both passes of a pair, so that the destination costs the same in both.

| threshold | no filter | with the filter | what it costs | again |
|---|---|---|---|---|
| 0.1 | 1.431 s | 1.550 s | 0.119 s | 1.445 and 1.568 s, 0.123 s |
| 0.03 | 1.436 s | 1.234 s | -0.202 s | 1.467 and 1.259 s, -0.208 s |

**At 0.03 bcftools is faster with the filter than without it**, by 0.20 s in
both sets. Its output is the difference: `view` writes the records it keeps
as text, 54773 of them instead of 100000, and the 45227 it does not write
save more than the filter costs. At 0.1, where it writes all 100000 either
way, the filter costs it 0.119 to 0.123 s. popnei's passes write nothing, so
their difference is the filter alone, and that is why the two programs cannot
be compared on the difference at 0.03. They can be compared on the whole
pass: at 0.03, of the first set of each, bcftools takes 1.234 s and popnei
0.105 s on 18 threads and 0.630 s on one.

## Every set with the load average it was taken at

The load average is the one minute figure of `uptime`, read right before the
invocation. Every popnei and bcftools row holds two, the first set and then
the second; the two rows of pyNei hold the first pair and then the second.

| set | no filter | with the filter |
|---|---|---|
| popnei VCF, 1 thread, 0.1 | 2.11, then 1.28 | 2.10, then 1.28 |
| popnei VCF, 1 thread, 0.03 | 1.93, then 1.26 | 1.93, then 1.24 |
| popnei VCF, 18 threads, 0.1 | 1.94, then 1.20 | 1.94, then 1.20 |
| popnei VCF, 18 threads, 0.03 | 1.87, then 2.94 | 1.87, then 2.94 |
| popnei vars, 1 thread, 0.1 | 3.51, then 2.71 | 3.51, then 2.71 |
| popnei vars, 1 thread, 0.03 | 3.31, then 2.57 | 3.31, then 2.57 |
| popnei vars, 18 threads, 0.1 | 3.12, then 2.71 | 3.12, then 2.71 |
| popnei vars, 18 threads, 0.03 | 2.95, then 2.57 | 2.95, then 2.57 |
| bcftools, 0.1 | 2.57, then 2.26 | 2.48, then 2.06 |
| bcftools, 0.03 | 2.36, then 1.90 | 2.30, then 1.83 |
| pyNei, 0.1, the two pairs | 2.16 and 1.08 | 2.31 and 1.23 |
| pyNei, 0.03, the two pairs | 2.58 and 1.31 | 1.53 and 1.17 |
| pyNei alternating, 0.1 | 1.06 at the start, 1.10 at the end | |
| pyNei alternating, 0.03 | 1.09 at the start, 1.29 at the end | |

The first pair of pyNei at 0.03 is the one the plan's "What could go wrong"
describes: its two halves started at 2.58 and at 1.53, so it was run again,
at 1.31 and 1.17, and that second pair is what the table of pyNei has as the
second pair. The repeat did not settle the question, and the alternating runs
are what did.

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

bcftools, from `/Users/jose/devel/popnei-bench/`, timed by a script in the
session's scratch directory that runs a command five times and prints the
median, with one run before them that is not timed:

    bcftools view -H -o /dev/null big.vcf
    bcftools view -H -i "F_MISSING<=0.1" -o /dev/null big.vcf

pyNei, by two scripts in the same scratch directory, which are not kept: one
times five passes of `sum(chunk.gts.gt_values.shape[0] for chunk in
vars_from_vcf(path).iter_vars_chunks())`, with
`filter_by_missing_data(variants, threshold)` between the two when there is a
threshold, and the other runs the pass with the filter and the pass without
it alternating in one process. Every run builds its own `Variants`, so no
chunk is read twice.

## Where the numbers that were not measured here come from

- The 1.65 s and 1.54 s of bcftools of "What has to be in place" of
  `docs/plans/filters.md` are one run each of 21 September 2026, with the
  records going to the terminal and not to `/dev/null`. The pair of that
  difference, 0.11 s, is the 0.119 to 0.123 s measured here.
- The VCF reader on its own, 0.563 s on one thread and 0.093 s on 18:
  `docs/reports/block-readers-measurement.md`, of the same file on the same
  machine. The passes with no filter here, 0.553 to 0.577 s and 0.093 to
  0.095 s, are that read.
