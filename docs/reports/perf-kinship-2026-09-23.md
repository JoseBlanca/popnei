# What the kinship costs, measured for the first time

23 September 2026. The kinship merged into `main` at `2a8b4b7` today, and
"Speed" of `docs/specs/kinship.md` states a target nobody had measured
popnei against: plink2's 0.23 s over 100000 variants x 1000 individuals.
This review measures it, says where the time goes, and says which of the
parts is anybody's to take.

The words this document uses. The **kinship** of two individuals says how
much more of their genome they share than two individuals drawn at random
from the same panel would. popnei computes it for every pair of the
individuals of a dataset in one pass: a **block**, the run of variants a
reader gives at a time and 5000 variants of 1000 individuals in every
measurement here, is **standardized** — each variant becomes one number
per individual, its dosage, with the mean of the variant taken off and
divided by the spread its allele frequency gives it — and the block of
those numbers is multiplied by itself into an individuals x individuals
**accumulator**, which after 20 blocks holds the sum over every variant
for every pair. Each of those sums is then divided by its **denominator**,
how many variants that pair had called in both of its individuals. When no
genotype of the dataset is missing that denominator is one number for
every pair; when a genotype is missing it is a second individuals x
individuals matrix, built from the **called mask**, a 1 for each genotype
whose every allele was read and a 0 for one with an allele missing, which
is multiplied by itself exactly as the dosages are. So **a dataset with a
missing genotype does twice the matrix work of one without**, which is why
the second product is skipped for a block with nothing missing.

The two programs popnei is measured against. **plink2** is the standard
tool of the field, and `--make-rel` is its command for this matrix.
**pyNei** is the Python library popnei is being rewritten from. The
**vars file** is popnei's own format for a dataset, which its readers are
fastest from, as the `.pgen` is plink2's; the three are different files of
the same variants.

What the arithmetic runs on. **Accelerate** is the linear algebra library
of macOS, which popnei calls natively and which numpy also calls; it is
the library the profiles below name `libBLAS.dylib`. **faer** is the
equivalent written in Rust, which runs in a browser tab, where there is no
Accelerate. The routine that does the kinship's product is **`dsyrk`**,
the standard name for multiplying a matrix by its own transpose and adding
the result into an accumulator; only half of that result is computed,
since it is the same on both sides of its diagonal. The **scan** is the
walk `crates/popnei-linalg` makes over every value it is given before
handing it to either library, refusing an infinity or a NaN, because the
two do not treat a NaN alike. **rayon** is the library that runs the rows
of a block on several threads, and **`Reblock`** is the reader popnei puts
in front of a source so that every block holds the same number of
variants, since a filter leaves blocks of uneven size and a product wants
an even one.

Which dataset a number is of. The **panels** of this review are the two of
100000 variants x 1000 individuals in section 3, on which every timing is
taken. The **reference panels** are two much smaller ones under
`tests/reference/`, 1200 variants x 200 individuals, which the tests
compare against plink2 and pyNei and which the tolerance figures come
from; they are too small to time anything.

What a measurement here is. A **sampling profile** stops a running
process a thousand times a second and records where it is; the **self
time** of a function is how many of those stops were inside that function
and not inside something it called. It is the only tool on this machine
that says which line holds a share of the time.

## 1. The scope and its limits

`crates/popnei/src/kinship.rs` at `2a8b4b7`, the row pass and the block
pass it drives in `crates/popnei/src/variant.rs`, and the two bindings,
`crates/popnei-python/src/kinship.rs` and `crates/popnei-js/src/kinship.rs`.
The principal components of a kinship are in scope and are measured in
section 3. Nine reviewers read it, one for each category of the
performance review skill, each with a fresh context.

**`main` moved ten commits while this review ran, and none of them changes
what it measures.** It is at `54036d1`, which brought the performance
review of the distances between populations and a fix to the kinship's
Ctrl-C. Of the files in scope, `crates/popnei/src/kinship.rs`,
`crates/popnei-js/src/kinship.rs` and `crates/popnei-linalg` are
untouched; `crates/popnei/src/variant.rs` gained 620 lines, all of them
the counting of the alleles and the genotypes of a population, which the
distances call and the kinship does not; and
`crates/popnei-python/src/kinship.rs` replaced two identical comment
blocks with one call of a helper their neighbours already use. Every
function of the kinship's row pass — `count_alleles`,
`the_codes_of_the_genotypes`, `the_counts_of_the_codes`,
`the_center_and_the_scale_of_the_dosages`, `the_major_allele`,
`the_standardized_row`, `the_standardized_rows`,
`the_standardized_block` and `the_standardized_rows_one_by_one` — is
byte-identical at the two commits, which was checked and not assumed.

The call frequency the reviewers were given, for 100000 variants x 1000
individuals, which is 20 blocks: the pass over the blocks runs once and
its loop 20 times; a row is standardized once per variant, 100000 times,
on the threads of rayon; the product runs once per block, and twice per
block when a genotype is missing; the eigendecomposition runs once per
run; and both bindings run once per call.

## 2. The verdict: apply, and the target is out of reach

Two things are true at once and the first is the one that decides what to
do next.

**The 0.23 s cannot be reached by making popnei's own code faster, and
this is now measured rather than argued.** Of the 0.399 s the fully-called
panel takes, 0.225 s is inside Accelerate computing the one product the
calculation is made of. That is 445 thousand million floating point
operations a second, within about 5 per 100 of what the routine does on
one thread of this machine: `docs/specs/pca.md` records the same call at
10.5 ms a block on one thread, `crates/popnei-linalg/benches/ops.rs`
measures the crate call at 10.973 ms, and the profile gives 11.3 ms. Three
measurements of different things agree. **So the product alone is 97 per
100 of plink2's entire 0.232 s.** Everything else popnei does — reading
the file, standardizing the rows, checking that no value is infinite —
would have to cost nothing at all, and the calculation would still only
draw level.

The reason is not that popnei reads more bytes. Counted over the
fully-called panel: the genotypes are 200 MB read twice, 400 MB; the
standardized dosages are 800 MB written, 800 MB read by the finiteness
scan and at least 800 MB read by the product. Packing the genotypes into 2
bits, as plink2 holds them, would remove at most 350 MB of 2.8 GB, 13 per
100 of the traffic, and none of it inside the library that holds 54 per
100 of the run. **plink2's advantage is that it never builds the matrix of
f64 at all**: it computes the cross product from the packed genotypes
themselves. Section 2.2 of `docs/rust_core.md` already names this — "in
compiled code the packed path is the fast one, which is the first thing a
Rust core gets to revisit" — and this review is where that comes due. It
is a change to how the kinship is computed, not a tuning of what is
there, and section 5 puts it to the owner as such.

**What is reachable is on the other panel, and it is large.** On the panel
with 3 in 100 genotypes missing, 41 per 100 of the run is spent working
out how many variants each pair had called in both of its individuals:
0.089 s building a 40 MB matrix of ones and zeros, and 0.230 s
multiplying that matrix by itself. Both can go. The same denominator is
inclusion and exclusion over exact integer counts — `kept` less the
variants each individual is missing, plus the variants both are missing —
and at 3 in 100 missing the pairwise term is about 2.2 million increments
a block against 5 thousand million floating point operations. Because both
routes sum whole numbers below 2^53, the denominators come out bit for
bit the same. That is H1, and it is **the largest gain this review found**: a matching
profile, a mechanism, and a proof that it changes no result.

Beside it, the next largest item on the fully-called panel is not
arithmetic at all: **the file is read on the same thread that then does
the product, so 0.110 s of reading never overlaps with 0.225 s of
computing.** `docs/architecture.md` section 3 asks for a read-ahead thread
one block ahead and there is none in the code. That is H2. It is worth up
to the whole 0.110 s and it changes no result, and **it is the one that
costs the most to build**: a thread, a channel, a second block alive, and
a path of its own for a browser, which has no threads. It also raises a
question only the owner can settle, which is whether a run with a reading
thread still counts as the "one thread" the comparison with `plink2
--threads 1` is made on. O5 puts that as a choice with its options.

**Two correctness matters were found by measuring and are in section 5.**
One is that the rule "no result of popnei depends on the number of
threads" is asserted by no test. The other is that on faer, which is the
browser's arithmetic, the number of components a kinship gives can change
with the number of threads — an integer, not a last bit.

## 3. What the kinship takes, which nothing had measured

Every number in this section was taken on the owner's Apple M5 Pro, 18
cores, rustc 1.98, the `bench` profile, on 23 September 2026, at a load
average between 1.3 and 2.7. Each popnei figure is the best of five timed
runs after one untimed run, since every other process on the machine can
only make a run longer, and each command was run twice.

The two invocations of each one-thread command agree to 0.5 per 100, and
what that costs to get is worth knowing for the next review: the same
commands at a load average of 11.5 to 16.2 give one-thread numbers 4 per
100 high and 18-thread numbers whose spread reaches 46 per 100, so a
figure from this benchmark at 18 threads is worth nothing unless the
machine is quiet.

### The two panels

Both are 100000 variants of 1000 individuals, diploid, and they hold the
same genotypes and differ only in which of them are `./.`. The second is
`big.vcf`, which earlier plans wrote and which every other measurement of
this project uses; the first is new, because "Speed" of
`docs/specs/kinship.md` states its target on a panel with every genotype
called and nothing in the repository could make one: the rate at which
`crates/popnei/benches/make_big_vcf.py` hid a genotype was a constant of
0.03. It is now a third argument, and the draw that decides which
genotypes are missing is made whatever the rate, so the panel of a rate of
0 is `big.vcf` with each of its missing genotypes restored to the value
the simulation actually drew for it: not another panel, and nothing
imputed.

| | every genotype called | 3 in 100 missing |
|---|---|---|
| the VCF | `bigcalled.vcf`, 403572954 bytes | `big.vcf`, 403572954 bytes |
| popnei's vars file | `bigcalled.vars`, 77820074 | `big.vars`, 81356714 |
| pyNei's vars file | `bigcalled.pynei.vars`, 32707850 | `big.pynei.vars`, 36454074 |
| plink2's own format | `bigcalled_plink.pgen`, about 24 MB | `big_plink.pgen`, 23850541 |

All of them are in `/Users/jose/devel/popnei-bench/`. Both panels give
100000 variants with variance and 1000 individuals on every run, and the
benchmark prints the mean of the diagonal, 1.064502 and 1.064496, and the
largest entry, 1.097653 and 1.096592, so the matrix was computed and not
left a shape of zeros.

Both panels were checked against plink2 before anything was timed. Over
every one of the 1000000 entries, the largest difference from the matrix
plink2 prints is 4.996e-06 on the fully-called panel and 4.989e-06 on the
other, 4.6e-06 of the largest entry in both. That is the precision
plink2's text output carries and not popnei's: "How it is verified" of
`docs/specs/kinship.md` records 1.83e-05 from the text of the same matrix
against about 1e-16 from its bits.

### What each program takes

One thread, 100000 variants x 1000 individuals. Each program reads the
format it is fastest from, which is the comparison of the whole job each
of them does: popnei its vars file, pyNei its own, plink2 its `.pgen`.

| | every genotype called | 3 in 100 missing |
|---|---|---|
| popnei, the core crate | **0.399 s** | **0.713 s** |
| popnei, called from Python | 0.408 s, 0.19 GB | 0.723 s, 0.25 GB |
| pyNei | 0.903 s, 0.40 GB | 1.579 s, 0.44 GB |
| plink2 `--make-rel square` | **0.232 s** | **0.805 s** |

popnei is the best of five after one untimed run, and the second
invocation of each gave 0.401 and 0.716 s. plink2 is the mean over 20 runs of the
whole command under `hyperfine`, the tool that runs a command again and
again and reports the spread, after one warm run; 232.3 ms with a standard deviation
of 6.7 and a range of 229.6 to 256.1, and 804.5 ms with a deviation of 7.4
and a range of 794.3 to 824.1. pyNei is the best of three, and its runs
spread by under 2 per 100. The memory is the largest resident the Python
process reached.

**plink2's 232.3 ms reproduces the 0.23 s that "Speed" of
`docs/specs/kinship.md` gives it**, so the target is confirmed as a figure
for the panel with every genotype called, which the spec says but which no
measurement in the repository had shown. Asking plink2 to write its matrix
as binary instead of as the 10 MB of text it writes by default moves it by
1 ms, 231.4 against 232.3, so the difference in what the two programs
write is not in these numbers.

**Three things that table says.**

**popnei misses the target by 0.167 s and is 1.72 times plink2** on the
panel the target is stated on.

**popnei is 1.13 times faster than plink2 on the panel with genotypes
missing**, 0.713 s against 0.805 s. The ranges do not overlap: plink2's
slowest of 20 runs is 824 ms and popnei's worst of ten is 802 ms. Both
programs are much slower on that panel than on the other, plink2 by 3.5
times and popnei by 1.8.

**popnei is 2.2 times faster than pyNei on both panels and holds half the
memory.** Section 2.1 of `docs/rust_core.md` gives pyNei 0.81 s for this
calculation at these sizes; the library as it stands today takes 0.903 s,
and that 0.81 s was measured in a numpy prototype and not in pyNei.

With the ten principal components of the matrix taken inside the same
clock, the fully-called run is 0.476 s and the other 0.755 s.

### What the principal components cost, which the spec asked for

`docs/specs/kinship.md` expects little room in the components and asks for
a measurement rather than an assumption. At 1000 individuals the
eigendecomposition and the projections cost **0.042 to 0.077 s**, between
9 and 16 per 100 of a run that takes them; the two measurements that far
apart were taken at load averages of 2.4 and 2.8 and the figure wants one
more clean pair before a decision turns on it. `docs/specs/linalg.md`
gives the routine underneath 0.035 s at that size, so the kinship and the
crate add little to it.

**The expectation holds at 1000 individuals and inverts above it.** The
pass over the variants costs the individuals squared times the variants,
while the eigendecomposition costs the individuals cubed: at 5000
individuals the decomposition is about 6.3 s beside a pass of about 5.6 s,
and at 10000 it is the larger of the two by more than twice. So the
components are a tenth of the kinship at the size the target is stated on
and the majority of it at the size `docs/objectives.md` names.

### Where the time goes

Two sampling profiles, `/usr/bin/sample <pid> 18` on the benchmark running
with `--runs 60`, one thread, no principal components. The idle threads of
the rayon pool, which show as `__psynch_cvwait` with 12894 and 12241
samples, are out of both totals.

With every genotype called, 12832 samples on the CPU:

| | samples | share | of the 0.415 s |
|---|---|---|---|
| the product, inside `libBLAS.dylib` | 6943 | 54.1 | 0.225 s |
| reading the vars file, almost all of it the decompression | 3399 | 26.5 | 0.110 s |
| standardizing the rows | 1820 | 14.2 | 0.059 s |
| the scan before each product | 406 | 3.2 | 0.013 s |
| asking whether the block has a missing genotype | 182 | 1.4 | 0.006 s |
| dividing each sum by its denominator | 11 | 0.1 | |

With 3 in 100 genotypes missing, 12222 samples on the CPU:

| | samples | share | of the 0.744 s |
|---|---|---|---|
| the two products | 7530 | 61.6 | 0.458 s |
| reading the vars file | 1786 | 14.6 | 0.109 s |
| writing the called mask | 1458 | 11.9 | 0.089 s |
| standardizing the rows | 930 | 7.6 | 0.057 s |
| the scan before each product | 440 | 3.6 | 0.027 s |
| asking whether the block has a missing genotype | 0 | | |

**The question of whether a block has a missing genotype costs more on the
panel that has none.** `any_missing` stops at the first missing genotype
it meets, so on the panel where 3 in 100 are missing it stops inside the
first row and costs nothing that the profile can see; on the panel where
none is, it has to read all 5000000 genotypes of every block to find that
out, 200 MB over the run for an answer that is always the same. The two
walks the code review named are therefore not one cost but two different
ones, and only the second is on the panel the target is stated on. They
were separated without a new measurement, by reading which instruction
address each sample fell on: the 182 samples sit at the offset of
`any_missing`, on the called panel alone, and the 1458 at the offset of
the loop inside `the_called_genotypes_of`, on the missing panel alone.

**The product is at the machine's floor and is most of plink2's whole
run.** Multiplying a 5000 x 1000 block by its own transpose into a
1000 x 1000 accumulator, computing half of it, is 5.005e9 floating point
operations, so the 20 blocks are 1.001e11, and 0.225 s of them is 445
thousand million a second.
`docs/reports/perf-linalg-2026-09-23.md` puts one thread of Accelerate at
about 530 billion a second on these products, and
`crates/popnei-linalg/benches/ops.rs` times that exact call, checks and
all, at 10.973 ms, which over 20 blocks is 0.219 s. Three independent
numbers agree, and the conclusion they force is in section 2.

**The 7.4 ms a block that "Speed" of `docs/specs/kinship.md` gives numpy
for the same product cannot be a one-thread number.** Over 20 blocks it is
0.148 s, which is 676 billion operations a second, above what that report
measures one thread of Accelerate to do. The spec says the measurement was
taken on "numpy 2.5.3 on Accelerate and its threads", so the figure is
right and what it is a figure of is not what the sentence beside it
implies. This is the same trap
`docs/reports/perf-linalg-2026-09-23.md` found in "Speed" of
`docs/specs/linalg.md`, a one-thread column beside a column that took the
threads it found, and whoever next edits either section should say which
is which.

## 4. The measurement plan

In the order in which each unblocks what follows it. Items 1 and 2 are
apparatus that every later experiment needs; the rest are the experiments
of section 6 in the order they should be run.

1. **A recipe that shows an experiment did not move the numbers.** Nothing
   in the repository reproduces the four figures the kinship's tolerance
   turns on, and no test varies the number of threads (H4). A script under
   `tmp/` prints, for the two small reference panels under
   `tests/reference/`, the largest difference from
   plink2's f64 matrix as a share of the largest entry, `num_vars`,
   `num_comps`, and the sha256 of the matrix and of the projections. It is
   run on both backends at 1, 3 and 8 threads:

       cargo test --workspace && uvx maturin develop --release
       for t in 1 3 8; do RAYON_NUM_THREADS=$t VECLIB_MAXIMUM_THREADS=$t \
           uv run python tmp/the_four_numbers.py; done
       cargo test --workspace --no-default-features
       uvx maturin develop --release --no-default-features
       for t in 1 3 8; do RAYON_NUM_THREADS=$t uv run python tmp/the_four_numbers.py; done

   What keeps a change: every test green on both backends; the four
   shares against 3.6e-16 and 4.5e-16 on Accelerate and 3.3e-15 and
   2.3e-15 on faer, with anything that doubles one of them reported to the
   owner although 1e-13 still passes, because the faer column has only 30
   times of headroom against Accelerate's 222; `num_vars` 1200 and
   `num_comps` 199 on both panels at all three thread counts; and the
   matrix's sha256 equal across thread counts, which is the thread rule
   asserted as equality and not as closeness.

2. **A checksum in the benchmark that covers more than the diagonal.**
   `crates/popnei/benches/kinship.rs` prints the mean of the diagonal and
   the largest entry, and both are diagonal entries, so an experiment that
   left the upper half of the matrix unmirrored or divided by the wrong
   count would print the same two numbers. It wants the mean of the
   off-diagonal entries too, and an assertion against a stored value.
   About six lines, proved by zeroing the upper triangle once and seeing
   the number move.

3. **The denominators without a matrix product**, H1. `cargo bench --bench
   kinship` on `big.vars`, one thread, best of five, two invocations a
   side. Gate first on a count that does not move between runs: a fresh
   `/usr/bin/sample` must show `the_denominators_of_the_block` and its
   `DSYRK` fallen from 5058 of 12222 samples. Keep if the run falls by more
   than 0.10 s from 0.713 s and the fully-called panel does not rise.

4. **Asking whether a block has a missing genotype from the row pass**,
   H3. Gate on the same profile: `the_denominators_of_the_block`, and the byte
   search of the standard library that it calls and that the profile names
   `memchr_aligned`, must together fall from 366 of 12832 samples on
   `bigcalled.vars` to under 20. Keep if the fully-called run falls at all; the site is
   0.012 s of a 0.167 s gap.

5. **Reading one block ahead**, H2. Before building anything, instrument:
   three `Duration`s around the reading of a block, the row pass and the
   product, printed by the benchmark, at one thread and at 18. That is
   about twelve lines and it decides whether to build the thread at all.
   Expect the read at 0.110 to 0.120 s at both thread counts; under 0.060 s
   refutes the finding. Then a `ReadAhead<R>` in `block.rs` with a channel
   of capacity one. Keep if the fully-called run at 18 threads falls below
   0.25 s, and gate on the matrix being byte-identical at 1 and 18 threads.

6. **What Accelerate's own threads give**, which nothing has measured:

       for k in 1 2 4 8 18; do VECLIB_MAXIMUM_THREADS=$k RAYON_NUM_THREADS=1 \
           cargo bench -p popnei-linalg --bench ops -- \
           --ops add_self_product_lower --runs 5; done

   The same run answers a correctness question the tests do not cover:
   whether the matrix is bit-identical at `VECLIB_MAXIMUM_THREADS` of 1
   and unset. It is two runs and a byte comparison, and every threading
   experiment should be gated on it (H4).

7. **The browser, which has no number at all.** A copy of
   `js/popnei/bench/time_pca.mjs` calling `calcKinship`, run with `node
   bench/time_kinship.mjs <vars> --runs 5` from `js/popnei` after `npm run
   build`, printing the clock and `process.memoryUsage().rss`. From
   numbers already in the repository — `docs/specs/pca.md` gives the same
   block product as 10.5 ms on Accelerate, 95 ms on faer natively and 187
   ms in wasm with `simd128`, which the shipped build has — the kinship
   should take about 3.7 s in a tab on the fully-called panel and about
   7.5 s on the other, with the product 85 to 90 per 100 of it against 54
   natively. The cheaper half of the measurement, with no browser, is
   `cargo bench -p popnei-linalg --bench ops --no-default-features --ops
   add_self_product_lower`.

8. **Peak memory at the sizes `docs/objectives.md` names**, which nothing
   has measured for the kinship either. `/usr/bin/time -l` on the
   benchmark at 10000 individuals, reading "maximum resident set size", as
   `docs/reports/block-readers-measurement.md` did. Counted from the code
   it should be about 850 MB with every genotype called and 1.69 GB with
   genotypes missing, before the three matrices the eigendecomposition
   needs.

What could not be measured here: whether faer is level with OpenBLAS on
x86, which is the open question of `docs/rust_core.md`, since this machine
has no OpenBLAS.

## 5. The build configuration

`Cargo.toml` at the workspace root sets `overflow-checks = false` and
`debug = "line-tables-only"` under `[profile.release]`, the same under
`[profile.bench]`, and `debug = true` under `[profile.profiling]`. `lto`,
`codegen-units` and `opt-level` are unset, there is no global allocator,
and `.cargo/config.toml` sets flags for the two wasm targets only.

All four of the usual settings are safe to sweep here, which is worth
stating because a performance review is where somebody would reach for
them: rustc emits no fast-math flag, so LLVM may not reassociate a sum of
`f64` nor contract a multiply and an add into one instruction, whatever it
inlines or vectorizes; `sqrt` is a hardware instruction and correctly
rounded on this processor; and there is no excess precision to lose. They
change the code around the arithmetic and not the arithmetic. What is not
safe, and what `checklists/numbers.md` forbids, is any fast-math flag,
because popnei's missing values become NaN at the Python boundary.

**Three of the four are not worth running.** Link time optimization and
one code generation unit were tried by
`docs/reports/perf-stats-2026-09-22.md`, which got 0.182 and 0.183 s
against 0.183 and 0.184 s for two and a half times the rebuild, and by
`docs/reports/perf-pca-2026-09-22.md`, which got nothing outside 0.004 s
on four settings; and 54 per 100 of this calculation is inside
`libBLAS.dylib`, a dynamic library no link time optimization reaches.
`target-cpu` moves little on this processor, where the vector instructions
are already in the baseline, and a wheel that is distributed cannot use
`native`. A different allocator is not indicated: the allocator appears in
neither profile above 6 samples of 12832, and the pass takes a handful of
multi-megabyte buffers that go to `mmap`.

**One is worth a count, not a timing.** Three functions of the row pass
appear as their own frames in the profile — `count_alleles` at 359
samples, `the_counts_of_the_codes` at 221 and
`the_codes_of_the_genotypes` at 128, together 5.5 per 100 of the
fully-called run — so at least one call site of each was not inlined,
although they are in the same crate as their caller, and 16 code
generation units is what stops that. The gate is not a wall time, which
would be under the noise: disassemble the built benchmark with `objdump
-d` and count the calls left to those three symbols, as it is built now
against `--config 'profile.bench.codegen-units=1'`. Only if the count
falls is a timing worth taking. `cargo asm`, which prints the machine code
of one named function, cannot answer this: it appends its own
`codegen-units=1` to every invocation, so it always shows the same
listing.

**The toolchain is still not pinned**, as the two reviews before this one
said. Every number here is against an unpinned rustc 1.98.0, and what the
compiler vectorizes changes between versions with no change in the source.

## 6. The findings

Numbered so that the prose can refer to them. A finding is named by what
it is, with its number after it.

### For the owner

**O1. Reaching 0.23 s means computing the product from the packed
genotypes, which is a change to the calculation and not a tuning of it.**
`crates/popnei/src/kinship.rs:450`. The numbers are in section 2. What is
being decided: whether popnei ever takes plink2's route, computing the
cross product of the standardized dosages from genotypes held two bits
each, without materializing the matrix of f64 that Accelerate multiplies
today.

The options.

- **Take it.** Section 2.2 of `docs/rust_core.md` measured the packed path
  in numpy and found it lost to floating point by 4 to 20 times there,
  because numpy cannot fuse the unpacking, the popcount and the
  accumulation into one pass, and it ends "in compiled code the packed
  path is the fast one, which is the first thing a Rust core gets to
  revisit". Nothing has measured it in compiled code. What it costs is a
  second way of computing the kinship, with its own tests against the
  first, and it would not help a browser as much as it helps here, since
  faer is 9 times Accelerate on the same product and the packed path would
  be popnei's own code on both.
- **Leave it, and record the target as unreachable.** "Speed" of
  `docs/specs/kinship.md` then says that the product alone is 0.225 s, 97
  per 100 of plink2's whole 0.232 s on this machine, so what the spec sets
  popnei against is not 0.23 s but something above 0.225 s plus whatever
  reading the file costs.
- **Take the smaller things instead**, H1, H2 and H3, which come to about
  0.31 s on the panel with genotypes missing and about 0.12 s on the other,
  and leave the fully-called panel at about 0.28 s against 0.232 s.

Recommended: leave it and record the target, and take H1 and H3 now. The
packed path is a plan of its own, it belongs with the association study
that will read the same genotypes, and a review is not where a second
implementation of the kinship should be decided. What makes this urgent
rather than idle is that `docs/specs/kinship.md` states a target the code
cannot reach, and a spec that states an unreachable number will be read as
a defect by whoever comes next.

**O2. No test asserts that a result does not depend on the number of
threads, though the rule is stated in three places.**
`crates/popnei/src/kinship.rs:1603`. Found by looking for the test that
would catch an experiment of this review, and not finding it.

No Rust test and no pytest varies `RAYON_NUM_THREADS` or
`VECLIB_MAXIMUM_THREADS`. The one test whose name promises bit equality,
`the_panel_read_in_blocks_of_37_variants_gives_the_matrix_bit_for_bit`,
says in its own doc comment that `Reblock` makes the 1200 variants of that
panel a single block, so it compares one block with one block and tests
`Reblock`. The property is very likely true — the row pass writes one row
per task and reads no other, the mask is literal 1.0 and 0.0 so every
partial sum is an exact integer, and `dsyrk` was measured bit-identical at
every pool size on both backends by
`docs/reports/perf-linalg-2026-09-23.md` — but **no test would fail if a
change broke it, and neither would any test fail if a change made the
matrix depend on the size of a block.** What it costs to fix: one test
that computes a kinship at two thread counts and compares the matrices as
equal, and one that computes it at two block sizes and compares them at
the spec's tolerance. This is a matter for a code review, which a
performance review happened to find; it is here because every experiment
below is gated on it.

**O3. On faer, how many components a kinship gives can change with the
number of threads.** `crates/popnei/src/kinship.rs:319-326`. faer is the
browser's arithmetic and the backend of a native build with
`--no-default-features`.

`docs/reports/perf-linalg-2026-09-23.md` O2 measured that faer's
eigendecomposition gives a different result at each pool size: over 800
eigenvalues the gap reaches 6.2e-15 relative. A component is given only if
its eigenvalue is above the largest eigenvalue times the individuals times
the difference between 1 and the next `f64`, which on the reference panel
is 7.67e-13; a move of 6.2e-15 relative on the largest eigenvalue is
1.07e-13 absolute, **14 per 100 of that threshold**. An eigenvalue inside
that band changes `num_comps` by one, which is the shape of the result and
not its last bits. The sign of a component can flip the same way, where
two projections are tied within 64 times the difference between 1 and the
next `f64`. Both are invisible on Accelerate and in a browser, which has
one thread; the exposure is the native build on faer, which is exactly
where the two faer tolerance figures were measured, at a thread count
nobody recorded. The options.

- **Record it**, adding to "How it is verified" of `docs/specs/kinship.md`
  a sentence saying that on faer the count of components is not
  reproducible across thread counts, and a test that computes a kinship's
  components at two thread counts on that backend and reports the count.
  Costs a paragraph and a test.
- **Pin faer's eigendecomposition**, which makes the count stable and
  costs what `docs/reports/perf-linalg-2026-09-23.md` measured: 1.43 times
  the wall time at 1000 individuals, 4.46 at 5000 and 5.80 at 10000.
- **Leave it silent**, since a browser has one thread and Accelerate is
  unaffected, so no user popnei has today can see it.

Recommended: record it. That review recommended the same for the last bits
of the eigenvalues and it is the right answer here too, but for a
different reason that is worth stating: the two places popnei's users
actually are, a browser and a native build on Accelerate, are both
unaffected, and paying two minutes at 10000 individuals to fix a build
nobody ships is a bad trade. What makes this worth a sentence in the
kinship's own spec rather than only the linear algebra's is that here the
consequence is the shape of the result and not its last bits.

**O5. Does a run with a thread that only reads still count as "one
thread"?** This is what H2 turns on and it cannot be settled here, because
it is about what popnei's published numbers mean and not about the code.

Every speed figure popnei states for the kinship is on one thread, and it
is compared against `plink2 --threads 1`. H2 would put the reading of the
file on a thread of its own, one block ahead. That thread does no
arithmetic; it waits on the disc and decompresses. The calculation would
still use one thread and the process would use two.

The options.

- **A reading thread counts as one thread**, on the ground that the second
  thread computes nothing and that plink2 overlaps its own reading with
  its own arithmetic within one thread, which popnei cannot do without a
  thread because its decompression is a library call. Then H2 is worth up
  to 0.110 s of the 0.399 s and the number stays comparable.
- **It does not**, and H2's gain is reported only in the figure that uses
  the threads the machine gives, where it is worth up to 0.110 s of
  0.283 s. The one-thread figure then stays as it is and the spec says
  which of the two the target is about.
- **Do not build it**, and keep the 0.110 s.

Recommended: that it counts, and that "Speed" of `docs/specs/kinship.md`
say in one sentence that a thread which only reads is not counted and why.
The alternative measures popnei against a constraint plink2 does not have.

**O4. The kinship has no browser target where the principal components
have one.** `docs/specs/kinship.md` "Speed". `docs/specs/pca.md` states 5
s and 7 s for a browser and popnei meets them; the kinship's spec states
nothing. From numbers already in the repository the kinship should take
about 3.7 s in a tab on the fully-called panel and about 7.5 s on the
other, with the product 85 to 90 per 100 of the run against 54 natively —
which makes H1 worth far more in a browser than it is here. Nothing has
measured it; item 7 of section 4 is the smallest measurement that would.
Memory may bind before time does: at 10000 individuals the pass holds an
800 MB accumulator and, when a genotype is missing, an 800 MB matrix of
denominators at the same time. That is 1.6 GB before the three matrices
the eigendecomposition needs, in an address space that WebAssembly caps at
4 GB and that most browsers hold well below.

### Hot-path

**H1. The denominators are a count of integers and do not need a matrix
product.** `crates/popnei/src/kinship.rs:502-558` and `:574-601`.
Confidence high.

- Hot-path evidence: on the panel with 3 in 100 missing,
  `the_denominators_of_the_block` calling `add_self_product_lower` into
  `DSYRK` is 3600 samples and the function's own loops are 1458, together
  **5058 of 12222 on-CPU samples, 41.4 per 100, about 0.30 s of the
  0.713 s run**. The 40 MB buffer is also the whole difference in
  footprint between the two profiles, 86.3 MB against 131.8 MB.
- Pattern: call the product that does less; an algorithm finding, which
  outranks any tuning.
- Mechanism: the denominator of a pair is `kept` less the used variants
  each of the two individuals is missing, plus the used variants both are
  missing. Every term is a count of whole things. Today it is 5 thousand
  million floating point operations a block over a matrix of ones and
  zeros that the code first writes as 40 MB of f64. At 3 in 100 missing a
  variant has about 30 individuals missing, so the pairwise term touches
  about 435 pairs a variant and 2.2 million a block: **about one operation
  for every 2300 the product does**, over a quarter of the data. The
  profile is what sizes the gain, 0.30 s of the 0.713 s run; that ratio
  only says why it is not close.
- Measurement plan: item 3 of section 4.
- Effect on the numbers: none, and provably. Both routes sum exact
  integers below 2^53, so the denominators are bit-identical at any thread
  count on both backends, and the numerator is untouched, so the four
  measured shares do not move. The test compares the two matrices as
  equal, not as close.
- Complexity cost: a per-variant list of the individuals missing, an
  accumulator of counts, one branch, and a crossover rule with its
  constant, since the pairwise term is quadratic in how many individuals a
  variant is missing in and stops being cheaper on a heavily missing
  panel. It drops the 40 MB buffer, and it may let the denominators be
  held as something smaller than f64, which would halve 800 MB at 10000
  individuals.
- Suggested experiment: section 7 has what it gave.

**H2. The file is read on the thread that then computes, so 0.110 s of
reading never overlaps with 0.225 s of product.**
`crates/popnei/src/kinship.rs:426`. Confidence high.

- Hot-path evidence: on the fully-called panel `Reblock::next_block` is
  3456 of 12892 main-thread samples, 26.8 per 100, with `DSYRK` its
  sibling on the same stack at 6981. The profile has two threads and one
  of them is idle. The share grows as the rest threads: the read is 0.110
  s of the 0.283 s that 18 threads take, 39 per 100.
- Pattern: `docs/architecture.md` section 3 asks for "a read ahead thread
  between the reader and the consumer, one block ahead, as pyNei does".
  `grep` for `thread::spawn`, `mpsc` and `crossbeam` across the core crate
  finds only code inside `#[cfg(test)]`: there is none.
- Mechanism: the four phases of a block run strictly one after another on
  one thread. A reader one block ahead overlaps the read of the next block
  with the product of this one; the ceiling of what it removes is the
  whole 0.110 s and what it actually removes is the smaller of the read
  and the compute, per block.
- Measurement plan: item 5 of section 4, whose first half is twelve lines
  of instrumentation and decides whether to build the second half.
- Effect on the numbers: none. The blocks reach the accumulator in file
  order and every sum is joined as it is today.
- Complexity cost: the largest here. A thread and a channel of capacity
  one in a module that has neither, a second block alive, an error that
  now crosses the channel, a wasm path of its own since a browser has no
  threads, and a join before `calc_kinship` returns, or the Python binding
  leaks the thread when a pass is interrupted. It belongs in `block.rs` as
  a wrapper every pass can use and not in the kinship. **It also changes
  what "one thread" means**, which is the basis of the comparison with
  `plink2 --threads 1`, and that is the owner's to settle.

**H3. Asking whether a block has a missing genotype costs more on the
panel that has none.** `crates/popnei/src/kinship.rs:512-513`. Confidence
high.

- Hot-path evidence: on the fully-called panel
  `the_denominators_of_the_block` and the `memchr_aligned` under it are
  366 of 12832 samples, 2.9 per 100, **0.012 s of the 0.399 s run and 7
  per 100 of the 0.167 s gap to plink2**. On the panel with genotypes
  missing the same site is 0 samples.
- Pattern: recomputing in a loop something the caller already has.
- Mechanism: `any_missing` stops at the first missing genotype, so where 3
  in 100 are missing it stops inside the first row; where none is, it
  reads all 5000000 genotypes of every block, 200 MB over the run, to
  learn an answer that is always the same.
  `the_center_and_the_scale_of_the_dosages` already sums how many
  genotypes of the row were called, and a row has one missing exactly when
  that sum is below the individuals.
- Measurement plan: item 4 of section 4.
- Effect on the numbers: none; the same blocks take the same branch.
- Complexity cost: `the_standardized_row` stops returning a `bool` and
  returns two fields, through `the_standardized_rows`,
  `the_standardized_rows_one_by_one` and `the_standardized_block`, and the
  analysis of the variants ignores the second. **It is also the gate that
  would let the called mask be written by the row pass without making the
  fully-called panel pay a 40 MB write per block it does not need today**,
  which is why it is worth its type even at 0.012 s.

### Likely

**L1. The lower half of the accumulator is scanned once per block, and
costs four times per value what the contiguous scan costs.**
`crates/popnei-linalg/src/lib.rs:1119-1134`. Confidence medium-high.
Evidence: both profiles split the two call sites — on the fully-called
panel the scan of the block is 282 samples over 5000000 values and the
scan of the accumulator's lower half is 124 over 500500, 4.4 times the
cost per value; 0.004 s of the run. Mechanism: it calls the scan once per
row, on prefixes of 1 to 1000 values, so a 1000-row matrix makes 1000
calls that each set up and reduce eight counters and walk a scalar tail,
over a triangle the previous product wrote and the 40 MB block has since
evicted. Plan: `cargo asm` for the fixed cost per row, then one contiguous
scan of the triangle against the row-by-row one in a trial binary. Effect
on the numbers: none, the predicate is unchanged. Cost: none for the loop
shape.

**L2. The finiteness scan over the mask of called genotypes and over their
counts can never fire.** `crates/popnei-linalg/src/lib.rs:245-246`
reached from `crates/popnei/src/kinship.rs:552`. Confidence high.
Evidence: 215 of 12222 samples on the panel with genotypes missing, 1.8
per 100, which is 1.0e8 values of the mask and 1.0e7 of the counts, half
of everything that run scans. Mechanism: the kinship writes those values
itself as the literals 1.0 and 0.0, and their sums are integers below
2^53, so no value can be an infinity or a NaN. Effect on the numbers:
none, provably. **Superseded by H1**, which removes the mask entirely; it
is filed because it stands alone if H1 is not taken.

**L3. The major allele of a variant is found by walking a table of 128
counts when the two that matter are already in hand.**
`crates/popnei/src/variant.rs:829-840`. Confidence medium. Evidence: the
row pass's closure is 1086 samples on the fully-called panel and the
offsets `sample` prints for it fall inside a 14-instruction scalar loop
run 128 times, about 1792 instructions a variant, inside a closure that holds
1086 of 12832 samples, 0.034 s of the run; what share of that closure this
one loop is was not measured, and the benchmark below is what would say.
Mechanism: on the fast path of `count_alleles`, which is a variant of two
alleles, only the entries 0 and 1 are not zero and both are held, so the
major allele is one comparison and the count of distinct alleles is
two. Plan: `cargo bench
--features bench-internals --bench standardize_row`; keep if the gap
between the four passes and the whole row falls by more than 0.5 ms a
block. Effect on the numbers: none **if the tie rule is kept**, the
lowest-numbered of two equally frequent alleles; a different tie-break
changes every dosage of such a variant. Cost: `count_alleles` gives back
two more values and its other callers move with it.

**L4. The matrix of denominators is written whole where only its lower
half is ever read.** `crates/popnei/src/kinship.rs:519-523` and `:533`.
Confidence high on the mechanism. This is memory, not time. Evidence,
counted rather than timed on a 10000 x 10000 matrix of f64 in a trial
binary: a matrix left zeroed holds 2 MB resident, one whose lower triangle
is written holds 536 MB, one written whole holds 783 MB — **242 MB per
800 MB matrix**. Mechanism: `the_entries_of` reads up to the diagonal, the
linear algebra checks and writes the lower half, and nothing reads above
it, yet both the starting count and the per-block addition are written
into every entry. It bites on the ordinary shape of a panel with sporadic
missingness, where the first blocks are fully called. Effect on the
numbers: none. Cost: the invariant that nothing reads the upper half
becomes load-bearing instead of incidental, and one test that it stays 0.

**L5. The size of a block was measured by nobody and is the wrong shape at
many individuals.** `crates/popnei/src/block.rs:28`. Confidence medium.
Evidence: the doc comments of both clamps say so; the constant is
5000000 genotypes, so 1000 individuals get 5000 variants a block and
10000 individuals get 500. Mechanism: at 1000 individuals the product
reads and writes an 8 MB accumulator per 5 thousand million operations,
625 operations a byte; at 10000 individuals and 500 variants it reads and
writes 800 MB per 50 thousand million, 62 a byte, ten times worse, and the
per-block scan of the triangle grows in the same proportion. Plan: a sweep
at 1000, 2500, 5000 and 10000 variants on both panels and on a panel of
5000 individuals; the block size is already an argument of the reader.
**Effect on the numbers: not none.** A different block size cuts the sum
over the variants differently, so every entry moves in its last bits, and
the bound at risk is the 30 times of headroom on faer, not the 222 on
Accelerate. The reference panels are one block at any setting, so they
cannot detect it: the experiment needs the 100000 x 1000 panel and a saved
matrix to compare against, and neither exists.

**L6. A partial eigendecomposition is worth 22 to 33 per 100 and matters
at 5000 individuals and above.** `crates/popnei/src/kinship.rs:311`.
Confidence medium. Evidence: `docs/specs/linalg.md` measured the routine
that computes a chosen range at 0.025, 0.18 and 4.9 s against 0.035, 0.27
and 6.3 s for the whole decomposition. The kinship needs the largest
eigenvalue and the ones asked for, and nothing in the result needs the
rest. Effect on the numbers: the two routines differ in the last bits, so
the spec's eigenvalue test must be re-run and reported. Cost: a fifteenth
operation in the crate, and **it gives nothing in a browser**, since faer
has no solver over a range — so it would widen the recorded gap between
the two backends.

**L7. The Python layer checks a matrix the core has just mirrored, with
three numpy calls per row.** `python/popnei/kinship.py:127`. Confidence
medium. Evidence: at 1000 individuals the whole boundary costs 0.009 s of
0.408 s, 2.3 per 100, so this site is a fraction of that and is not worth
a change at the size the target is stated on. Counted: 3000 numpy calls at
1000 individuals, one of them a transposed column slice per row, and 25
bytes of temporaries per individual per row — 25 MB at 1000 individuals
and **2.5 GB at the 10000 the objectives name**, which is where it
matters. The matrix popnei built cannot fail the check:
`the_lower_half_mirrored` copies each lower entry into its upper cell, so
the two cells of a pair are bit-identical and a test asserts it. Cost: a
private route that skips the check for a matrix the pass produced, and the
loss that a defect leaving the upper half unmirrored would reach a user
who never asks for components.

**L8. The JavaScript layer runs three loops over the whole matrix on the
crossing path.** `js/popnei/src/kinship.ts:139-140`. Confidence medium,
pattern only: no browser or node number for the kinship exists. The
constructor checks every value finite, finds the largest absolute value
and checks symmetry with a transposed read — about 250 million element
visits at 10000 individuals, in JavaScript, on a matrix the core mirrored
bit for bit. Plan: item 7 of section 4, which is also what would give the
kinship a browser number at all.

### Notes

- **Dividing each sum by its denominator and mirroring the lower half into
  the upper are cold**, 11 samples of 12832 over 500500 and 1000000
  values. Both are quadratic in the individuals: at 10000 the mirror is 50
  million stores at an 80 KB stride, about 0.1 s against a product of
  about 22 s, so still 0.2 per 100. Sample it again before assuming it,
  if that size is ever timed.
- **The row pass collects 5000 results of 64 bytes each to carry 5000
  bits**, `crates/popnei/src/variant.rs:1464`: `popnei::Error` is 64 bytes,
  so `Result<bool, Error>` is 64, and a block allocates 320000 bytes for
  5000 bools. It is 10 samples of 12832 and under the noise of any
  timing; it is here because the same shape at a larger error type would
  not be.
- **`Reblock` copies nothing on these files and 200 MB a pass on a file a
  user wrote at another size.** `crates/popnei/src/kinship.rs:180` asks
  for the default, which is what these panels hold, so the block passes
  through untouched; `python/popnei/io_vars.py` lets a user write a vars
  file at any size, and then every block is copied. Count the copies
  before doing anything: they are 0 today.
- **The benchmark's checksum is a function of the diagonal alone**, so an
  experiment that broke the mirroring or the division would print the same
  two numbers. Item 2 of section 4.
- **The sampler slows the program by about 25 per 100** — runs inside the
  profiling window took 0.514 to 0.575 s against 0.407 to 0.431 s outside
  it — and ran at 717 samples a second where it reports 1000. The shares
  in section 3 are sound; turning them into seconds of an undisturbed run
  assumes the slowdown fell evenly on every frame, which was not checked.
  A second profile at a coarser interval would check it.

## 7. Seen outside the scope

- **"Speed" of `docs/specs/kinship.md` gives numpy 7.4 ms a block for the
  product, and that is a figure taken with Accelerate's own threads.**
  "Speed" of `docs/specs/pca.md` prints both sides of the same
  measurement, "10.5 ms, and 7.4 ms with the threads Accelerate takes by
  itself". Over 20 blocks, 7.4 ms is 0.148 s, which would be 676 thousand
  million operations a second, above what one thread of this machine does.
  The kinship's spec carries the threaded number beside a sentence about
  one thread. This is the same trap
  `docs/reports/perf-linalg-2026-09-23.md` found in "Speed" of
  `docs/specs/linalg.md`, a one-thread column beside one that took the
  threads it found.
- `crates/popnei/src/io/vars.rs` zeroes 78 MB a pass before `read_exact`
  overwrites it, which is already L5 of
  `docs/reports/perf-io-2026-09-22.md`.
- The genotype column's validity bitmap decompresses 1.25 MB of all-ones
  per batch, 133 of 12832 samples, about 1 per 100. Whether a vars file
  needs it belongs to that module.
- `crates/popnei-js/src/kinship.rs:69-71` gives the matrix to
  wasm-bindgen, the layer that carries values between Rust and JavaScript,
  which copies it into an array JavaScript can read: 800 MB in wasm and
  800 MB in the JavaScript heap at once at 10000 individuals, and the
  memory of wasm never gives back what it grew by.
- The compression of the vars file is already the owner's open decision
  from `docs/reports/perf-pca-2026-09-22.md`, where the same reader cost
  0.107 s of a 0.801 s analysis. Here it is 0.110 s of 0.399 s: the same
  absolute cost and a larger share, because the kinship is the shorter
  calculation. It needs no new decision; H2 is what is new.

## 8. What the code already does well

- **The second product is already skipped for a block with nothing
  missing.** `crates/popnei/src/kinship.rs:456-466`. It is why the
  fully-called panel pays 20 products and not 40, and it is the reason
  the two panels differ by 1.8 times rather than by more.
- **The two buffers of a block are kept across the blocks.**
  `crates/popnei/src/kinship.rs:424-425`. A pass over a million variants
  asks the machine for them once, which the profile confirms: the zeroing
  of memory is 44 samples of 12832 and none of it is under the resizes.
- **The turn from the layout popnei holds to the one the routines read is
  done with flags, not copies.** `crates/popnei-linalg/src/blas.rs:57-82`
  passes the half, the transpose and the leading dimension so that neither
  product copies anything, and the faer backend asks for the lower
  triangle so it computes the half too.
