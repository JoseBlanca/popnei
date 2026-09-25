# What the per variant and the per individual statistics cost a pass, task 6.1

22 September 2026. Task 6.1 of `docs/plans/stats.md`, deliverable 1 of its
work package 6, measures what the two calculations of popnei's stats module
cost a whole pass over a 400 MB dataset, against the four numbers that
"Speed" of `docs/specs/stats.md` set before the code was written, and against
pyNei on the same variants. Three of the four numbers are met and one is not:
the per variant pass with the five statistics takes 0.479 to 0.483 s on one
thread, where the spec asks for 0.25 s. The task changed no code of the
library; it added the two scripts that timed popnei and pyNei,
`crates/popnei/benches/time_stats.py` and
`crates/popnei/benches/time_pynei_stats.py`, and this report. The last
section says what it leaves to a performance review.

The words this document uses. A **pass** is one whole read of the dataset
from opening the file to the last variant, with whatever is calculated over
the variants on the way; every pass here is timed from the call that opens
the file to the result, in Python, and every run opens the file again;
popnei's opening is `open_vars`, which gives the variants of a vars file, and
pyNei's is `load_vars`. The **read alone** is the pass that calculates
nothing: popnei's `Variants.iter_blocks`, the variants of the file in blocks
one after another, with the genotypes as the only field of each block, which
is what `docs/specs/stats.md` takes its numbers to reach from. A **vars
file** is a
file of a library's own format, written once from a VCF and read many times;
popnei's is arrow IPC with the genotypes compressed with lz4, which
`docs/specs/io_vars.md` describes, and pyNei's is its own, so the two files
of this dataset are not the same bytes and their reads are not compared
below. The **five statistics** are what `calc_per_var_distribs` gives for
every variant: the observed heterozygosity, the major allele frequency, the
expected heterozygosity plain and unbiased, and the polymorphism ratio.
`calc_per_individual_stats` is the other calculation of the module, a pass of
its own that gives two numbers for every individual instead: the share of the
variants at which its genotype is missing and the share of its called
genotypes at which it is heterozygous. A
**population** is a named set of individuals that a statistic is calculated
over on its own; with none named there is one population of every individual
of the file, and the runs with populations use four of 250, the individuals
of the file in the order they are in it, 250 to each. The **threads** of a
popnei run are rayon's, which popnei takes from the environment because it
builds no pool of its own, so a run is `RAYON_NUM_THREADS=1` or
`RAYON_NUM_THREADS=18`; pyNei takes its threads as the argument
`num_threads`. Every number is the **best of 5 runs**, with one run before
them that is not timed, and the **load average** beside it is the one minute
figure of `uptime`, read right before the invocation.

## The machine, the build and the files

The owner's Apple M5 Pro, 18 cores, 64 GB, macOS 27.0, native
`aarch64-apple-darwin`, rustc and cargo 1.98.0. popnei is at b0360f8 of the
branch `plan/stats`, built into the virtual environment of the worktree with
`uv run maturin develop --release`, whose profile is `[profile.release]`:
optimized, with `overflow-checks = false`. Both libraries run under CPython
3.14.5 with numpy 2.5.3, pandas 3.0.6 and pyarrow 25.0.1. pyNei is at
ef0ca6e, the commit `pyproject.toml` names; numpy is on Accelerate, the BLAS
of macOS, which the counting of these statistics does not call.

The files, both outside the repository in `/Users/jose/devel/popnei-bench/`
and both in the page cache, since the run before the timed ones reads the
file. `big.vars`, 81356714 bytes, popnei's vars file of `big.vcf`, 100000
variants of 1000 diploid individuals whose genotypes are missing at a rate of
0.03, which `docs/reports/filters-measurement.md` describes.
`big.pynei.vars`, 36454074 bytes, pyNei's vars file of the same VCF, which
was already there, modified at 13:16 on 22 September 2026, before this task,
which is the day the table of "Speed" of `docs/specs/stats.md` was taken; so
this task wrote no vars file. Read back, it gives
100000 variants of 1000 samples in chunks of 5000, which is the dataset
popnei reads.

The machine was not idle. Three sets of the eight popnei passes were taken:
the first at load averages of 7.0 to 7.6, which were the decay of this
session's own `cargo` build of popnei at 795 per cent of one core, and the
second and the third at 1.6 and at 2.3. The first set is 0.03 to 0.04 s
slower on every pass than the other two and is not in the tables below; it is
in "Every set with the load average it was taken at", and it changes none of
the four verdicts. The second and the third set agree within 0.005 s on every
pass.

## What the two calculations cost popnei

Both sets, best of 5 runs each, on `big.vars` through `open_vars`: the second
at a load average of 1.55 to 1.61, the third at 2.24 to 2.38.

| | 1 thread | 18 cores |
|---|---|---|
| the read alone, the genotypes as the only field | 0.102, 0.104 s | 0.102, 0.105 s |
| `calc_per_var_distribs`, the five statistics, no populations | 0.479, 0.483 s | 0.139, 0.143 s |
| the same, 4 populations of 250 | 0.599, 0.623 s | 0.154, 0.158 s |
| `calc_per_individual_stats` | 0.202, 0.207 s | 0.118, 0.122 s |

The read alone is 0.102 s on one thread and on 18, the same figure that
`docs/reports/filters-measurement.md` measured on 21 September 2026 with a
Rust benchmark and that "Speed" of the stats spec takes its numbers from.
The reader of a vars file runs on the thread that calls it, so over this file
the threads are the calculation's alone.

## The four numbers of "Speed"

Each is a whole pass over `big.vars`, and each is met when both sets are
under it.

| what the spec asks | measured | met |
|---|---|---|
| the five statistics, no populations, 1 thread: 0.25 s | 0.479, 0.483 s | no, 1.9 times over |
| the same, 18 cores: 0.15 s | 0.139, 0.143 s | yes, by 0.007 to 0.011 s |
| `calc_per_individual_stats`, 1 thread: 0.25 s | 0.202, 0.207 s | yes, by 0.043 to 0.048 s |
| the same, 18 cores: 0.15 s | 0.118, 0.122 s | yes, by 0.028 to 0.032 s |

The one that is not met is 0.23 s over. Whether to work on it is the owner's,
through a performance review: this task changed no code to reach a number.

The one that is met by least, the five statistics on 18 cores, is met by
0.007 s in the slower of the two sets, which is a twentieth of the number.
The set taken at a load average of 7.0, where every other pass was 0.03 to
0.04 s slower, gave 0.146 s for it, still under 0.15 s.

How the four numbers were arrived at, which is what the measurement is
against: "Speed" of the spec took the read alone, 0.102 s and 0.102 s, and
added twice what the missing data filter adds to the same pass, and rounded
up. That filter is the one of `docs/specs/filters.md` that keeps the variants
whose missing genotypes, divided by the individuals of the dataset, are at
most a threshold: it reads every row of genotypes once and counts, and
`docs/reports/filters-measurement.md` measured what it adds to this same pass
over this same file, 0.056 to 0.060 s on one thread and 0.013 to 0.018 s on
18. Twice is because the statistics read each row twice where the filter
reads it once. On 18 cores that
holds: the five statistics add 0.037 to 0.038 s to the read, and twice what
the filter adds there is 0.026 to 0.036 s. On one thread it does not: the
five statistics add 0.377 to 0.379 s to the read, which is 6.3 to 6.8 times
what the filter's one read adds, not twice.

## Where the time of the per variant pass on one thread is

The pass reads each row twice, as the spec says: once to count the genotypes
of the variant, which is what the observed heterozygosity needs, and once to
count its alleles, which the major allele frequency, the two expected
heterozygosities and the polymorphism ratio all share. Timing the pass with
one statistic asked for says what each of the two reads costs. One set, at a
load average of 1.86 to 2.02:

| | 1 thread | 18 cores |
|---|---|---|
| the read alone | 0.102, 0.104 s | 0.102, 0.105 s |
| the observed heterozygosity alone, the genotype count | 0.193 s | 0.117 s |
| the major allele frequency alone, the allele count | 0.383 s | 0.130 s |
| the five statistics | 0.479, 0.483 s | 0.139, 0.143 s |

On one thread the genotype count adds 0.089 to 0.091 s to the read and the
allele count adds 0.279 to 0.281 s, three times as much over the same
100000 rows of 2000 alleles. The five together add 0.377 to 0.379 s, which is
the two reads one after the other plus 0.005 to 0.011 s for the three
statistics that take the allele counts the major allele frequency already
made. The 0.25 s the spec asked for is 0.102 s of read plus twice the 0.056
to 0.060 s of the filter's one read, 0.214 to 0.222 s, rounded up. The first
of the two reads costs about what that assumed, 0.089 to 0.091 s; the second
costs three times it.

What the two counts do differs: the genotype count walks the alleles of each
genotype and compares them with one another, one value written per genotype,
while the allele count clears an array of 128 counts, one for each allele a
variant could carry, and then adds one to the count of that allele for every
called allele of the variant, a write whose place in the array the allele
decides. Which part of that is the 0.28 s was not measured: no sampling
profile was taken in this task.

## What the four populations cost

Four populations of 250 individuals, against no populations, which is one
population of all 1000:

| | 1 thread | 18 cores |
|---|---|---|
| the five statistics, no populations | 0.479, 0.483 s | 0.139, 0.143 s |
| the same, 4 populations of 250 | 0.599, 0.623 s | 0.154, 0.158 s |
| what the four populations cost | 0.120, 0.140 s | 0.015, 0.015 s |

Every individual is in exactly one of the four populations, so the two passes
count the same 100000000 genotypes and the same 200000000 alleles; what the
four add is the work that is per variant and per population and not per
genotype. Over 100000 variants, 0.120 to 0.140 s on one thread is 1.2 to
1.4 µs a variant for the three extra populations, and on 18 cores 0.015 s is
0.15 µs a variant. "Speed" of the spec had no number for this and asked for
one; this is it.

## pyNei on the same variants

pyNei at ef0ca6e over `big.pynei.vars`, best of 5 runs, at load averages of
1.38 to 2.44, `num_threads` 1 and 6. Its `calc_per_var_distribs` has four
statistics where popnei's has five, because pyNei's `exp_het` holds either
the plain expected heterozygosity or the unbiased one and popnei gives both;
its `calc_per_sample_stats` is popnei's `calc_per_individual_stats`.

| | 1 thread | 6 threads |
|---|---|---|
| its read alone, the chunks and nothing else | 0.168 s | |
| `calc_per_var_distribs`, the four statistics, no populations | 1.126 s | 0.264 s |
| the same, `obs_het` alone | 0.992 s | 0.240 s |
| the same, `maf` alone | 0.185 s | 0.187 s |
| the same, the four, 4 populations of 250 | 2.624 s | 0.598 s |
| `calc_per_sample_stats` | 1.024 s | 0.248 s |

These are the figures the table of "Speed" of the stats spec already had from
22 September 2026, within 0.05 s on every row except its read, which was
0.187 s there and 0.168 s here.

popnei against pyNei, each library's own pass over its own vars file of the
same variants, on one thread:

| | popnei | pyNei | popnei is |
|---|---|---|---|
| the per variant statistics, no populations | 0.479 s | 1.126 s | 2.4 times faster |
| the same, 4 populations of 250 | 0.599 s | 2.624 s | 4.4 times faster |
| the per individual statistics | 0.202 s | 1.024 s | 5.1 times faster |

popnei on 18 cores against pyNei on 6, which is each library at the threads
it was timed at here and not the same machinery:

| | popnei, 18 | pyNei, 6 | popnei is |
|---|---|---|---|
| the per variant statistics, no populations | 0.139 s | 0.264 s | 1.9 times faster |
| the same, 4 populations of 250 | 0.154 s | 0.598 s | 3.9 times faster |
| the per individual statistics | 0.118 s | 0.248 s | 2.1 times faster |

"Speed" of the spec expected the per variant pass on one thread to be 4.5
times under pyNei's 1.13 s. It is 2.4 times under.

The two libraries are dear in opposite places. pyNei's observed
heterozygosity costs it 0.992 s and its major allele frequency 0.185 s;
popnei's observed heterozygosity costs 0.091 s and its major allele frequency
0.281 s. Over the four populations, pyNei pays 1.498 s on one thread where
popnei pays 0.120 to 0.140 s.

## What this does not measure

- **wasm.** Nothing was run in the browser. "Speed" of the spec leaves the
  number to reach there to the first measurement, and this is not it; task
  6.2 builds the wheel for pyodide and runs the two calculations in it, but
  it takes no time.
- **A pass over a VCF.** Every number here is over a vars file. Over
  `big.vcf` the read alone is 0.561 s on one thread and 0.095 s on 18,
  `docs/reports/filters-measurement.md`, and what the statistics add to that
  read was not measured.
- **What the Python layer adds.** The times are of the Python call, so each
  holds the result the binding builds: for the per variant pass four
  histograms of 40 bins and five counts, one of each per population, and for
  the per individual pass two series of 1000 values. How much of a pass that
  is was not measured.
- **Other shapes of dataset.** One file, 100000 variants of 1000 individuals
  at ploidy 2 with 3 in 100 genotypes missing. How the four populations cost
  changes with more of them, and how the passes behave with more individuals
  and fewer variants, is not here.

## Every set with the load average it was taken at

The load average is the one minute figure of `uptime`, read right before the
invocation. Each popnei row gives the three sets in order, the first of them
the one taken while this session's own build was still on the machine.

| set | best of 5 | the load average |
|---|---|---|
| popnei, the read alone, 1 thread | 0.110, 0.102, 0.104 s | 7.61, 1.61, 2.26 |
| popnei, the five statistics, 1 thread | 0.513, 0.479, 0.483 s | 7.61, 1.61, 2.26 |
| popnei, the five, 4 populations, 1 thread | 0.642, 0.599, 0.623 s | 7.61, 1.56, 2.24 |
| popnei, the per individual statistics, 1 thread | 0.220, 0.202, 0.207 s | 7.24, 1.59, 2.38 |
| popnei, the read alone, 18 cores | 0.109, 0.102, 0.105 s | 6.98, 1.59, 2.38 |
| popnei, the five statistics, 18 cores | 0.146, 0.139, 0.143 s | 6.98, 1.59, 2.38 |
| popnei, the five, 4 populations, 18 cores | 0.160, 0.154, 0.158 s | 6.98, 1.55, 2.38 |
| popnei, the per individual statistics, 18 cores | 0.124, 0.118, 0.122 s | 6.98, 1.55, 2.38 |
| popnei, `obs_het` alone, 1 thread | 0.193 s | 2.02 |
| popnei, `maf` alone, 1 thread | 0.383 s | 2.02 |
| popnei, `obs_het` alone, 18 cores | 0.117 s | 1.94 |
| popnei, `maf` alone, 18 cores | 0.130 s | 1.94 |
| pyNei, the read alone, 1 thread | 0.168 s | 1.38 |
| pyNei, the four statistics, 1 thread | 1.126 s | 1.38 |
| pyNei, the four, 4 populations, 1 thread | 2.624 s | 1.40 |
| pyNei, `calc_per_sample_stats`, 1 thread | 1.024 s | 1.54 |
| pyNei, the four statistics, 6 threads | 0.264 s | 1.49 |
| pyNei, the four, 4 populations, 6 threads | 0.598 s | 1.49 |
| pyNei, `calc_per_sample_stats`, 6 threads | 0.248 s | 1.69 |
| pyNei, `obs_het` alone, 1 thread | 0.992 s | 1.94 |
| pyNei, `maf` alone, 1 thread | 0.185 s | 1.86 |
| pyNei, `obs_het` alone, 6 threads | 0.240 s | 1.86 |
| pyNei, `maf` alone, 6 threads | 0.187 s | 2.44 |

The worst run of a set is within 0.005 s of its best everywhere except
popnei's per individual statistics on one thread in the first set, 0.220 to
0.262 s, and pyNei's four statistics on one thread, 1.126 to 1.165 s.

## The commands, as they were run

popnei, from `crates/popnei/benches` of the worktree
`.claude/worktrees/stats`, after `uv run maturin develop --release`, with
`uptime` before each invocation:

    RAYON_NUM_THREADS=1 uv run python time_stats.py \
        /Users/jose/devel/popnei-bench/big.vars per-var 5

with `read`, `per-var`, `per-var-pops`, `per-var-obs-het`, `per-var-maf` and
`per-individual` in the place of `per-var`, and `RAYON_NUM_THREADS=18` for
the 18 cores.

pyNei, from the same directory, whose threads are the last argument:

    uv run python time_pynei_stats.py \
        /Users/jose/devel/popnei-bench/big.pynei.vars per-var 5 1

with `read`, `per-var`, `per-var-pops`, `per-var-obs-het`, `per-var-maf` and
`per-sample`, and 6 in the place of the 1.

Each script's docstring says what every pass it can time does. pyNei is a
development dependency of popnei at the commit `pyproject.toml` names, so
`uv run` is what has it.

## What this hands to a performance review

- **The allele count of a variant**, which is 0.279 to 0.281 s of the 0.377
  to 0.379 s the five statistics add on one thread, and the part of the pass
  that would have to come down for the number to be met: the pass is 0.229 to
  0.233 s over 0.25 s, and the allele count is 0.28 s of it. A sampling
  profile of the pass with the major
  allele frequency alone, on one thread, is what would say where inside it
  the time is. The array it counts into is 128 `u32`, 512 bytes, cleared once
  for every variant and every population, and the file has 2 alleles a
  variant.
- **The four populations on one thread**, 0.120 to 0.140 s, which is per
  variant and per population and not per genotype. How it grows with the
  populations was not measured, and a dataset of 50 populations is an
  ordinary one.
- **The per variant pass on one thread is 3.4 times the same pass on 18
  cores**, 0.479 s against 0.139 s, where the read alone is 0.102 s on both.
  The calculation is what the threads divide, so the 0.377 s it adds on one
  thread becomes 0.037 s on 18, which is a tenth of it and not an
  eighteenth.
