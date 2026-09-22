# What the principal component analysis of 100000 variants takes, task 4.1

22 September 2026. Task 4.1 of `docs/plans/pca.md`, the first of its work
package 4, measures on the code what "Speed" of `docs/specs/pca.md` states from
the trial, a crate written before this code to try the options out, which is
not committed: how long the principal component analysis of 100000 variants of
1000 individuals takes natively, what pyNei and plink2 take on the same
variants, and whether the compiler turned the loop that standardizes a row into
vector instructions. It changed no code of the library; it added the benchmark
`crates/popnei/benches/pca_vars.rs`, the script `time_pca.py` beside it that
times the two libraries from Python, and this report.

**The target was missed.** The number to reach is 0.3 s for 100000 variants of
1000 individuals on one thread with no weights, and the analysis takes 0.801 s,
2.7 times that. With the threads the machine gives it takes 0.352 s. Where the
time goes is below: 0.41 s of the 0.801 s is the standardizing of the 20 blocks,
13.7 times the 1.5 ms per block the trial measured, and the machine code shows
why. pyNei takes 8.1 s and 5.6 GB for the same analysis, so popnei is 10 times
faster than what it replaces and 3.2 times slower than plink2, which takes
0.248 s from its own file. A miss is not a reason to change the code inside this
plan, as `docs/plans/pca.md` says.

The words this document uses. **The analysis** is
`do_pca_from_variants`, which both libraries have: every variant becomes one
number per individual, its **dosage**, how many alleles of its genotype are not
the most frequent allele of the variant, and the principal components of those
numbers are computed. To **standardize** a row is to turn the genotypes of one
variant into its dosages with the mean of the variant taken off and divided by
its standard deviation, which is what a genotype with an allele missing gets in
the place of a dosage. A **block** is the array of genotypes a reader of popnei
gives at a time, 5000 variants of 1000 individuals here, and a **pass** is one
whole read of a dataset, block after block. **G** is the product of the
standardized blocks with themselves added up, individuals x individuals, and
the eigenvectors of G give the projections. The **weights** are the weight of
each variant in each component, which need the eigenvectors and so a **second
pass** over the same variants; `num_prin_comps` is how many components they are
given for, 0 meaning none and no second pass. A **vars file** is popnei's own
format, arrow IPC with the genotypes compressed with lz4, which
`docs/specs/io_vars.md` describes; pyNei has a format of its own with the same
name, and plink2 has its **pgen**, which holds each genotype in 2 bits. **The
threads** are not asked for on a command line here: popnei standardizes the rows
of a block on rayon's global pool, which reads `RAYON_NUM_THREADS` when it is
built, and the matrix work of both libraries runs in the BLAS of the system,
which on this machine is Apple's Accelerate and reads `VECLIB_MAXIMUM_THREADS`
when the process starts, so one thread is asked for by setting both variables in
the environment of the command and the threads the machine gives by setting
neither.

## The machine, the builds and the files

The owner's Apple M5 Pro, 18 cores, 64 GB, macOS 27.0 build 26A428, native
`aarch64-apple-darwin`, rustc and cargo 1.98.0 (88d9e12ae 2026-08-18). popnei is
built by `cargo bench`, whose profile inherits from `[profile.release]`, with the
cargo feature `blas`, which is on by default and links Accelerate; the Python
module is what `uv run maturin develop --release` left in the environment, the
same code through pyo3. pyNei is at ef0ca6e, the commit `pyproject.toml` names,
under CPython 3.14.5 with numpy 2.5.3 and pandas 3.0.6, and numpy is linked
against Accelerate too. plink2 is v2.0.0-a.7.7 M1 of 18 September 2026, at
`/opt/homebrew/bin/plink2`.

No two timings were taken at once and nothing of popnei was built while one ran.
The machine was not idle: the owner's own processes kept the load average between
1.7 and 4.1 through the session, and the one number taken at a load average above
3 is marked where it appears. Each set of runs is preceded by one run that is not
timed, so that the timed runs read the file from the page cache and pay no page
fault of the first touch of the memory a pass works in.

The files, all outside the repository in `/Users/jose/devel/popnei-bench/`, and
all written before this task by earlier plans, so none of them was written again:

| file | bytes | what reads it |
|---|---|---|
| `big.vcf` | 403572954 | popnei, pyNei and plink2 |
| `big.vars` | 81356714 | popnei |
| `big.pynei.vars` | 36454074 | pyNei |
| `big_plink.pgen`, `.pvar`, `.psam` | 23850541, 2367848, 8009 | plink2 |

`big.vcf` is what `crates/popnei/benches/make_big_vcf.py` writes with no number
of variants given: 100000 variants of 1000 individuals, diploid, whose genotypes
are missing at a rate of 0.03. That it holds those sizes was checked with `grep
-vc '^#'`, which counts 100000 data lines, and with the columns of its `#CHROM`
line after the ninth, which are 1000. `big.vars` is what popnei's `write_vars`
wrote of it, 20 batches of 5000 variants, and the benchmark prints `100000
variants, 1000 individuals` of it on every run. The files hold the same
variants: the analysis of `big.vars` and of `big.vcf` both give 6.228 per 100 of
the variance to the first component, and pyNei gives 6.2284 to it over
`big.pynei.vars`. `big_plink.pgen` was written by `plink2 --vcf big.vcf
--make-pgen`, which is the format plink2 is fastest from and the one the trial's
0.26 s was measured on.

## What the analysis takes

Every time is the best of 5 runs, and the best is what the table states, since
every other process on the machine can only make a run longer; the median and the
worst of each set are in the paragraph after it. The two columns are
`num_prin_comps` 0, which is one pass over the variants and no weights, and 10,
which adds a second pass.

| popnei, 100000 x 1000 | no weights | weights for 10 components |
|---|---|---|
| the vars file, one thread | 0.801 s | 1.353 s |
| the vars file, 18 threads | 0.352 s | 0.537 s |
| the vars file, one thread, from Python | 0.823 s | 1.370 s |
| the VCF, one thread | 1.329 s | not measured |

The target of "Speed" of `docs/specs/pca.md`, 0.3 s on one thread with
`num_prin_comps` 0, is missed by 0.501 s: the analysis takes 2.7 times it. With
the threads the machine gives, which the target does not ask for, it takes 0.352 s
and misses by 0.052 s.

The spread of each set was small. On one thread with no weights the five runs
ran from 0.801 to 0.837 s, and with the weights from 1.353 to 1.369 s; on 18
threads, 0.352 to 0.353 s and 0.537 to 0.541 s. The VCF is the one file whose
runs spread, 1.329 to 1.515 s, and its median is 1.397 s. The two rows of Python
were taken at load averages of 4.06 and 3.51, the highest of the session, and
they are the only numbers here taken above 3.

Three things that table says. **The second pass costs 0.552 s on one thread**,
which is 1.353 s less 0.801 s: it reads the file again and standardizes every
block again, and the product it then does, a block by the 10 scaled eigenvectors,
is small beside the product of the first pass. **The binding and the frames cost
0.022 s**, which is the Python row less the row above it, 2.7 in 100 of the time:
what a user calls from Python is the same work plus the building of the pandas
frames of the projections and of the percentages. **The VCF costs 0.528 s more
than the vars file** on one thread, which is what parsing 403 MB of text takes
against reading 81 MB of lz4, and agrees with the 0.561 s pass over the VCF and
the 0.102 s pass over the vars file of `docs/reports/filters-measurement.md`.

## What pyNei and plink2 take

The same analysis, the same 100000 variants of 1000 individuals, each program
reading the file it is fastest from and then the same VCF. Every time is the best
of 5 runs for plink2 and of 3 for pyNei, whose runs take 8 to 22 s each; the peak
memory is the largest resident memory the process reached, which for pyNei is the
matrix of every dosage of the dataset. `meanimpute` is what makes plink2 give a
missing genotype the mean of the dosages of its variant, as popnei and pyNei do;
without it plink2 computes something else, and "Speed" of `docs/specs/pca.md`
says what.

| program, one thread | its own file | the VCF | peak memory |
|---|---|---|---|
| popnei, no weights | 0.801 s | 1.329 s | 0.21 GB |
| pyNei, weights for every component | 8.106 s | 22.267 s | 5.59 GB |
| plink2 `--pca 10 meanimpute` | 0.248 s | 0.507 s | not measured |

popnei is 10.1 times faster than pyNei over its own file and 16.8 times faster
over the VCF, and it holds 0.21 GB where pyNei holds 5.59 GB, 27 times less.
pyNei's 8.106 s and 5.59 GB are the 8.2 s and 6.7 GB that "Speed" of
`docs/specs/pca.md` carries from the trial, measured again here. plink2 takes
0.248 s from its pgen, the 0.26 s of the trial, and popnei takes 3.2 times that;
from the VCF plink2 takes 0.507 s, which includes writing the temporary pgen it
then reads, and popnei takes 2.6 times that.

What the three programs do is not the same work, and the differences all favour
popnei in the comparison with pyNei and plink2 in the comparison with popnei.
pyNei builds the matrix of every dosage of the dataset in memory, 100000 x 1000
in `f64`, which is what its 5.59 GB is, and gives the weight of every variant in
every component, which popnei gives for 10 at most; the popnei number it is
beside asks for no weight at all. plink2 reads genotypes packed in 2 bits, four
to the byte, where popnei reads one byte per allele and two bytes per genotype,
and it computes 10 eigenvectors where popnei computes all 999. The 989 components popnei
computes besides are not where its time goes: `docs/specs/linalg.md` measured
the routine of LAPACK that gives every eigenvalue and eigenvector of a symmetric
matrix, `dsyevd`, at 0.035 s for these 1000 individuals, and the one that gives
a chosen range of them, `dsyevr`, at 0.025 s for the 10 largest, so computing 10
instead of all of them would save 0.010 s of the 0.801 s.

## Where the time goes on one thread

`sample`, the sampling profiler of macOS, took 20 s of a process running the
analysis of `big.vars` again and again with no weights on one thread, one sample
every millisecond. The shares below are of the 15373 samples that landed in the
run; the times are those shares of the 0.801 s run. The profiler makes the run
longer, 0.909 s as the median of the 40 runs of the process it attached to
against 0.804 s with no profiler on it, which is why the third column is worked
out from the shares and not from the sampled times.

| part of the analysis | samples | share | time | per block |
|---|---|---|---|---|
| standardizing the rows of the blocks | 7929 | 51.6 | 0.413 s | 20.6 ms |
| the product of each block with itself | 4871 | 31.7 | 0.254 s | 12.7 ms |
| reading the vars file | 2050 | 13.3 | 0.107 s | 5.3 ms |
| the eigendecomposition | 457 | 3.0 | 0.024 s | |
| the signs, the projections and the rest | 66 | 0.4 | 0.003 s | |

Two of those lines are known from other measurements and agree with them, which
is what says the other three can be read as they stand. The product of a block
of 5000 x 1000 with itself takes 12.05 ms with `dsyrk`, the routine of BLAS for
the product of a matrix with itself, on one thread, of which 1.68 ms is the two
scans for a value that is not finite that the linalg crate,
`crates/popnei-linalg`, through which all the linear algebra of popnei goes,
does around the call, both measured on 22 September 2026 and in "Speed" of
`docs/specs/linalg.md`; the profile gives 12.7 ms for the whole call and 1.3 ms
for the part of it that is popnei's own code. A pass over `big.vars` with the
genotypes alone asked for and nothing computed takes 0.104 s, the best of 5 runs
of `cargo bench --bench filter_vars` with no filter on one thread, taken in this
session; the profile gives 0.107 s.

**The standardizing is where the target went.** The trial of "Speed" of
`docs/specs/pca.md` measured 1.5 ms for a block of 5000 x 1000 with the loop
written in two passes for the compiler to vectorize, and 3.9 ms for the loop as
one would first write it. This code takes 20.6 ms for the same block, 13.7 times
the trial's number and 5.3 times the plain loop's. The 0.3 s of the target is 20
blocks at 12 ms, which is the product and 1.5 ms of standardizing, and an
eigendecomposition of 0.04 s; the reading of the vars file, 0.104 s, was not in
it either. With the standardizing at the trial's 1.5 ms per block the analysis
would take 0.42 s here, and it would still miss the target by 0.12 s, which is
the reading.

## What the machine code shows

`cargo asm -p popnei --lib`, with cargo-show-asm 0.2.60, prints the machine code
of the release build of the library. The three functions that standardize a row
do not appear in it under their own names: the compiler inlined all three, and
`count_alleles` with them, into the closure that rayon runs for one row, which the
listing calls
`<rayon::iter::map_with::MapWithFolder<_, _, _> as
rayon::iter::plumbing::Folder<_>>::consume_iter::with::<…,
&popnei::pca::the_standardized_rows::{closure#1}>::{closure#0}` and prints as 1664
instructions. The vector instructions of this machine are NEON, whose registers
are `v0` to `v31` and whose instructions here are `ldp q`, which loads 32 bytes,
`cmeq`, which compares 16 bytes at a time, and `udot`, which adds four groups of
four bytes into four 32 bit lanes.

**`the_codes_of_the_genotypes`, which writes the dosage of each genotype into a
byte, was vectorized over the wrong axis and runs no vector instruction on
diploid data.** The compiler vectorized the inner loop, the one over the alleles
of a single genotype, and not the outer one over the individuals: the block at
`LBB198_37` loads 64 alleles at a time with `ldp q20, q21` and `ldp q24, q25`,
compares them with the major allele with `cmeq.16b`, adds with `udot.4s` and
tests for the missing allele with a second `cmeq.16b`, and the branch that leads
to it is `cmp x25, #64` on the ploidy, which is 2 in this dataset. There is a
second path of 8 alleles at a time at `LBB198_42`, behind `cmp x25, #8`. So at a
ploidy of 2 the code falls to the scalar tail at `LBB198_45`, which reads one
allele per iteration with `ldrb` and counts with `cinc`, and it runs the outer
loop's own work once per genotype.

**`the_counts_of_the_codes`, which counts how many genotypes have each dosage,
was vectorized and the vectorized path is the one that runs.** The block at
`LBB198_57` loads 64 codes at a time with two `ldp q`, compares each with the
dosage being counted with `cmeq.16b` and adds with `udot.4s`; the branch to it is
`cmp x21, #64` on the codes left in the run, and a run is 255 codes of the 1000
individuals, so three quarters of each run go through it, the next 56 through the
8 wide path at `LBB198_61` and the last few through the scalar one at `LBB198_63`.
This is the loop the spec wrote for the compiler to vectorize, and the compiler
did.

**The lookup at the tail of `the_standardized_row`, which turns each code into
its value, was not vectorized, and no compiler could vectorize it as it stands.**
It reads a table of 256 `f64` at an index that is a byte of memory, which is a
gather, and NEON has no gather instruction. The compiler unrolled it four codes
at a time at `LBB198_97`: four `ldrb` of the codes, four scalar `ldr d` of the
table, and two `stp` that store the four values in pairs; the remainder goes one
at a time at `LBB198_100`.

There is a fourth pass over the row that the spec does not count, and it is
scalar. `count_alleles`, of `crates/popnei/src/variant.rs`, runs before the three
above to find the major allele, and it walks the 2000 alleles of the row one byte
at a time at `LBB198_13`, reading each with `ldrb`, skipping the missing one and
adding one to the entry of that allele in a table of 128 counts with a load, an
add and a store to an address the value itself gives. That is a second gather,
and a store that depends on the value read, which is why it is where it is in the
machine code and not in vector registers. How much of the 20.6 ms per block it
takes was not measured: the profiler sees the four passes as one function,
because the compiler inlined them all into the closure, so the only way to split
them is to measure each on its own, which no benchmark here does.

## What was not measured, and why

- **How the 20.6 ms per block splits between the four passes over a row.** The
  compiler inlines them into one function and the profiler reports that function,
  so the split needs a benchmark that times each pass by itself, which this task
  did not write.
- **The memory of plink2 and the memory of popnei on one thread beyond its peak
  of 0.21 GB.** The peaks here are what the system tells a process about
  itself, asked for inside the Python process that ran the two libraries, and
  plink2 does not run in it.
- **The analysis over a VCF with weights**, and the whole of the wasm side, which
  is task 4.2 of `docs/plans/pca.md`.
- **Anything of a machine that is not this one.** Every number is of the owner's
  Apple M5 Pro with Accelerate, and the shares of the profile above all depend on
  Accelerate's `dsyrk`, which another BLAS would not match.

## How to get the numbers again

From the root of the worktree, with `/Users/jose/devel/popnei-bench/` holding the
files of the table above. The two variables go before the command and not into
it: Accelerate reads its own when the process starts.

```text
VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1 \
    cargo bench --bench pca_vars -- /Users/jose/devel/popnei-bench/big.vars \
    --runs 5 --num-prin-comps 0
cargo bench --bench pca_vars -- /Users/jose/devel/popnei-bench/big.vars \
    --runs 5 --num-prin-comps 0
VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1 \
    uv run python crates/popnei/benches/time_pca.py \
    /Users/jose/devel/popnei-bench/big.pynei.vars 3
VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1 \
    uv run python crates/popnei/benches/time_pca.py \
    /Users/jose/devel/popnei-bench/big.vars 5 --popnei --num-prin-comps 0
uv run --no-project python crates/popnei/benches/time_command.py 5 -- \
    plink2 --pfile big_plink --pca 10 meanimpute --threads 1 --out plink_pca
cargo asm -p popnei --lib pca 6
```

The second of those commands is the first one again with the threads the machine
gives, which is what setting neither variable means. The 6 of the last one is
where the closure that standardizes a row stood in the listing of this build:
`cargo asm -p popnei --lib pca` with no number prints the list to choose from,
and the number changes with the code.

The profile: the benchmark is run from its built binary under
`target/release/deps/` so that `sample` has a process to attach to, with
`--runs 40`, and `sample <pid> 20 1 -file <out>` takes 20 s of it at one sample
per millisecond. `time_pca.py` times pyNei without `--popnei` and popnei with it,
and reads a path that does not end in `.vcf` as the vars file of whichever
library it is timing.
