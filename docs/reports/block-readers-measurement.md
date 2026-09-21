# The measurement of the VCF reader, task 2.6

21 September 2026. Task 2.6 of `docs/plans/block-readers.md` measures the VCF
reader that work package 2 rewrote, against the targets of "Speed" of
`docs/specs/io_vcf.md`, and sets the three constants that the spec leaves to a
measurement. It changed no code but those constants, in three commits, and the
spec's table of "Speed" in a fourth. All four targets are met.

The words this document uses. The spike is `spike/pynei_spike` of pyNei, a
trial parser in Rust whose time is what the spec asks the reader to reach, or
no more than a tenth above it. A block is the array of genotypes the reader
gives at a time, 5 million of them, which for the 1000 diploid individuals of
the file below is 2500 variants. A batch is the lines the reader takes from
the file before it parses them side by side on the threads of rayon; a batch
ends where its block does, and two constants bound it, `LINES_PER_BATCH` in
lines and `BYTES_PER_BATCH` in bytes of text. The buffer of the file is what
`VcfReader::from_path` keeps between two calls to the file system. The serial
read of the lines is the phase in which one thread cuts the lines out of that
buffer while the others wait.

The machine: the owner's Apple M5 Pro, 18 cores, macOS 27.0, native
`aarch64-apple-darwin`, release build, the file in the page cache. Nothing of
popnei built while a timing ran. The machine was not idle: the owner's
`mediaanalysisd` and a `CGPDFService` ran through the session and were not
closed, and the load average went from 1.4 to 6.8. Every set of runs below
carries the load average it was taken at.

The files. `/Users/jose/devel/popnei-bench/big.vcf`, 403572954 bytes, 100000
variants of 1000 individuals, made by `crates/popnei/benches/make_big_vcf.py`,
and `big.vcf.gz`, the same file bgzipped, 37695707 bytes.
`/Users/jose/devel/popnei-bench/wide10000.vcf`, 120172952 bytes, 3000 variants
of 10000 individuals, for the memory, made by `make_wide_vcf.py` beside it,
which is that script with the two sizes on its command line and everything
else, the seed 42 included, unchanged.

## The four timings, and the targets

The benchmark is `crates/popnei/benches/read_vcf.rs`, which times one whole
read of the file with the genotypes alone asked for and the default options.
Each number is the median of 5 runs, and each set of 5 was repeated three
times and more, so that the spread between sets says what a difference has to
beat.

| | the reader at b6e98fa | the reader at 5cbf389 | its spread over the sets | the target | met |
|---|---|---|---|---|---|
| plain, 1 thread | 0.593 s | 0.563 s | 0.560 to 0.563 | 0.594 s | yes |
| plain, 18 threads | 0.119 s | 0.093 s | 0.093 to 0.094 | 0.108 s | yes |
| bgzipped, 1 thread | 0.908 s | 0.890 s | 0.859 to 0.940 | 0.924 s | yes |
| bgzipped, 18 threads | 0.419 s | 0.394 s | 0.393 to 0.396 | 0.44 s | yes |

b6e98fa is the commit this task started from, the reader as the rest of the
plan left it; 5cbf389 is the reader with the three constants of the next
section. The targets are the spike's times of 20 September 2026 plus a tenth,
which `docs/specs/io_vcf.md` states and this task did not choose.

The sets, each with the load average at which it was taken:

| | set 1 | set 2 | set 3 | load |
|---|---|---|---|---|
| b6e98fa plain, 1 thread | 0.612 | 0.586 | 0.593 | 2.42, 2.32, 2.02 |
| b6e98fa plain, 18 threads | 0.119 | 0.120 | 0.119 | 2.38, 2.21, 2.02 |
| b6e98fa bgzipped, 1 thread | 0.855 | 0.908 | 0.913 | 2.38, 2.21, 2.02 |
| b6e98fa bgzipped, 18 threads | 0.419 | 0.415 | 0.432 | 2.43, 2.11, 2.02 |
| 5cbf389 plain, 1 thread | 0.560 | 0.563 | 0.563 | 2.40, 2.01, 1.78 |
| 5cbf389 plain, 18 threads | 0.093 | 0.094 | 0.093 | 2.29, 1.93, 1.72 |
| 5cbf389 bgzipped, 18 threads | 0.396 | 0.393 | 0.394 | 2.18, 1.85, 1.66 |

The bgzipped read on one thread was taken fourteen times, because its spread
is what decides whether its target is met: 0.889, 0.929, 0.933, 0.866, 0.909,
0.913, 0.914, 0.866, 0.870, 0.864, 0.859, 0.890, 0.940 and 0.885 s, at load
averages between 1.4 and 2.3. The median of the fourteen is 0.890 s and they
fall in two clusters, one near 0.866 and one near 0.913, with no load average
to tell them apart; what makes the difference is not known.

The spike was built from a copy of `/Users/jose/devel/pynei/spike/pynei_spike`
in `/Users/jose/devel/popnei-bench/spike_copy`, never in place, and timed the
same day on the same two files: `cargo` took 22 s to build it. Its times, the
median of 5 runs, three sets, at its own chunk of 5000 variants: 0.523, 0.523
and 0.531 s plain on one thread, 0.103, 0.102 and 0.103 s plain on 18, 0.806,
0.805 and 0.807 s bgzipped on one, 0.392 s three times bgzipped on 18. At a
chunk of 2500 variants, which is the block popnei gives for this file, the
same four are 0.523, 0.114, 0.805 and 0.421 s. The threads are rayon's global
pool, set with `RAYON_NUM_THREADS`.

The spike is faster today than on 20 September on three of the four, so the
targets computed from today's spike are tighter: 0.575, 0.113, 0.887 and
0.431 s. Against those, three are met and the bgzipped read on one thread,
0.890 s, is 0.003 s over, which is a thirtieth of the 0.081 s its fourteen
sets spread over: the measurement does not tell the two apart. On 18 threads
the reader is faster than the spike on both files, 0.093 s against 0.103 and
0.394 s against 0.392.

## The constants

Each was changed only where three interleaved sets of runs showed a difference
larger than the spread between sets, and each has a commit of its own.

**The buffer of the file, 8 KiB to 256 KiB**, commit 3c4bf5d, a new constant
`BYTES_OF_THE_FILE_BUFFER`. `from_path` opened the file with
`BufReader::new`, which gives 8 KiB; every line comes out of that buffer, and
a buffer that runs out is a call to the file system that the one thread
reading the lines waits for. Plain, 18 threads, with the batches as they were,
1024 lines and 8 MiB: 8 KiB 0.118 s, 64 KiB 0.108, 256 KiB 0.106, 1 MiB 0.105.
With the batches as they are now: 0.105, 0.095, 0.093, 0.093. Bgzipped, 18
threads: 0.418 s and 0.411 s for 8 KiB and 256 KiB. On one thread the
capacities are the same read, plain 0.59 to 0.62 s and bgzipped 0.88 to
0.95 s. 256 KiB is where the gain stops, so 1 MiB buys nothing for four times
the bytes. This confirms what the reviewers of the reader reported, 122-129 ms
down to 105-107 ms on 18 threads with 256 KiB; the numbers here are lower
because the machine was quieter.

**`LINES_PER_BATCH`, 1024 to 4096**, commit 52e506f. A batch ends where its
block does, so the last batch of a block holds what is left of it and is the
one the threads share worst: with 1024 lines the block of 2500 variants of
this file was read in batches of 1024, 1024 and 452 lines, and 452 lines are
25 to a thread on 18 cores. Plain, 18 threads, with the 256 KiB buffer: 256
lines 0.139 s, 1024 0.105, 2048 0.098, 4096 0.098. Bgzipped, 18 threads: 0.476,
0.425 and 0.411 s for 256, 1024 and 4096. On one thread the sizes are the same
read, plain 0.58 to 0.62 s and bgzipped 0.90 to 0.95 s, so there is no trade
here: nothing is lost on one thread. 2048 and 4096 are the same read because
the bound in bytes decides at both, 2500 lines of this file being 10.1 MB
against a bound of 8 MiB. Sizes past the block buy nothing: with the 8 KiB
buffer, 4096 lines and 32 MiB read the file in 0.107 s, 8192 and 64 MiB in
0.106 and 16384 and 128 MiB in 0.106.

**`BYTES_PER_BATCH`, 8 MiB to 16 MiB**, commit 5cbf389. With 4096 lines
allowed, this bound is what cuts a block into batches: at 8 MiB a batch got
about 2077 of the 2500 lines of a block and a block was read in two, at 16 MiB
in one. Plain, 18 threads: 0.098 s to 0.094. Bgzipped, 18 threads: 0.411 s to
0.392. On one thread the two are the same read. What it costs is memory, in
the next section, and that is what this bound is for; the wasm build reads one
line at a time and none of it reaches wasm.

Taken together the three take the plain read on 18 threads from 0.119 s to
0.093, 22 in 100 less, and the bgzipped one from 0.419 s to 0.394, 6 in 100
less. On one thread the plain read went from 0.593 s to 0.563 and the bgzipped
one from 0.908 s to 0.890; those two differences are inside the spread between
sets and the measurement does not attribute them.

## The memory of a reader

Two measurements of the same reads, on 18 threads: the maximum resident set
size of the whole benchmark, which holds the binary and the pages the
allocator keeps as well as the reader, and the most bytes alive at once,
counted by a global allocator that adds every allocation and subtracts every
free, in a scratch program outside the repository,
`/Users/jose/devel/popnei-bench/live_bytes/`.

| | 8 MiB and 1024 lines | 16 MiB and 4096 lines |
|---|---|---|
| 1000 individuals, maximum resident | 20.0 MB | 33.5 MB |
| 1000 individuals, alive at once | 16.6 MB | 35.7 MB |
| 10000 individuals, maximum resident | 24.4 MB | 32.8 MB |
| 10000 individuals, alive at once | 22.5 MB | 34.6 MB |

On one thread the resident sizes are 19.2 and 33.4 MB for 1000 individuals and
23.7 and 32.1 MB for 10000. The file of 10000 individuals is 3000 variants in
120172952 bytes; its block is 250 variants and about 10 MB of text, so the
bound in bytes is what decides there and not the 4096 lines.

What is alive at once is the genotypes of the block, 5.0 MB at 2500 variants
of 1000 diploid individuals; the text of a batch, 10.1 MB of lines at the
bounds the reader now has and 4.1 MB at the ones it had; the 2500 small rows
that say where each line is and what its parse gave; and the buffer of the
file, 256 KiB. The peak is about twice their sum because each of those buffers
grows by doubling and a growth holds the old buffer and the new one at once: a
reader of this file with one line in a batch, whose text is 4 KB, already
peaks at 10.2 MB against genotypes of 5.0 MB. Raising one bound at a time
gives 16.6 MB for 1024 lines with either 8 or 16 MiB, since the line bound
decides there, and 23.3 MB for 4096 lines with 8 MiB.

## Where the time goes

A sampling profile with `/usr/bin/sample` over 30 s of a build with line
tables, in a target directory of its own so that the timed build was not
touched. The profiles are kept at
`/Users/jose/devel/popnei-bench/profile_1thread.txt` and
`profile_18threads.txt`.

On one thread, of the 23038 samples of the thread that does the work; the
other thread is the caller, asleep inside rayon's `install`:

| what the code does | samples | self time |
|---|---|---|
| the columns of the individuals, `fill_row_genotypes` | 21109 | 91.6 in 100 |
| the search for the end of a line in the buffer of the file | 625 | 2.7 |
| the read from the file system | 568 | 2.5 |
| the nine first columns, parsed as text: the UTF-8 check, the split at a tab, the comparison | 297 | 1.3 |
| the copy of a line out of the buffer of the file | 187 | 0.8 |
| the genotypes of a new block set to missing | 58 | 0.25 |
| the pass that finds the FILTER | 14 | 0.06 |
| the header, read once | 4 | 0.02 |

By phase, the parse of the rows is 21521 samples, 93.4 in 100, and the read of
the lines 1402, 6.1 in 100. The serial append after a batch does not appear:
the benchmark asks for the genotypes alone, so there is no position, id,
quality or allele to append. "How it runs" of the spec has what asking for
every column costs, 38 ms of 650 on one thread.

On 18 threads the thread that calls `next_block` has 9363 samples, of which
8647 are inside the read: 4308 in the serial read of the lines, 1887 of them
searching for the end of a line, 1810 in the read from the file system and 554
copying lines out of the buffer, and 4063 in its own share of the parse. The
18 workers of rayon are idle, in `__psynch_cvwait`, `__psynch_mutexwait` or
`swtch_pri`, in 103069 of their 168534 samples, 61 in 100.

So the serial read of the lines is still what bounds the 18 threads, as the
reviewers of the reader found: it is 46 in 100 of the wall time of the thread
that reads. The target on 18 threads is met without removing it. The read
ahead thread of section 3 of `docs/architecture.md` is what would remove it,
and it is outside this plan and was not built.

## What was tried and gave nothing

The reader wraps the gzip decoder in a `BufReader::new`, 8 KiB, so the same
argument that held for the buffer of the file could hold for it. Built with
256 KiB there and timed against the reader at 5cbf389 on the bgzipped file on
one thread, three interleaved sets each: 0.849, 0.876 and 0.873 s against
0.866, 0.870 and 0.864 s. No difference. The change was not committed and the
working tree was restored with `git checkout`.

Where a gain is still to be had, with what it stands on. The columns of the
individuals hold 91.6 in 100 of the one thread profile, and the reader is
within 8 in 100 of the spike there, 0.563 s against 0.523, while doing three
checks the spike does not: the ploidy of each genotype, its allele numbers
against ALT, and the FILTER. On 18 threads what is left is the serial read of
the lines, above. The parser was not changed, as the task says.

## The checks

Run at e297e63, the last commit of the task.

| check | last line |
|---|---|
| `cargo fmt --all --check` | exit 0, no output |
| `cargo clippy --workspace --all-targets -- -D warnings` | `Finished dev profile [unoptimized + debuginfo] target(s) in 0.26s` |
| `cargo test --workspace` | `test result: ok. 153 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out` |
| `cargo wasm-check` | `Finished dev profile [unoptimized + debuginfo] target(s) in 0.19s` |
| `uv run ruff format --check` | `12 files already formatted` |
| `uv run ruff check` | `All checks passed!` |
| `uv run maturin develop && uv run pytest` | `61 passed in 1.35s` |
| `npm run build && npm test` in `js/popnei` | `tests 45`, `pass 45`, `fail 0` |

No test named the old value of a constant, so none had to be corrected: the
tests that vary the size of a batch use `LINES_PER_BATCH` itself or literals
of their own.

## The commits

- 3c4bf5d, a reader that opens a path holds 256 KiB of it and not 8 KiB
- 52e506f, a batch of the vcf reader is 4096 lines and not 1024
- 5cbf389, a batch of the vcf reader holds 16 MiB of text and not 8
- e297e63, the speed of the vcf reader spec has the reader that was built

The doc comments of `LINES_PER_BATCH` and `BYTES_PER_BATCH` were rewritten so
that every number in them is of this reader, with the date, the file, the
machine, the threads and the build, and they point at
`docs/reports/vcf-to-blocks.md` for the reader before this plan.

## The commands, as they were run

From `/Users/jose/devel/popnei/.claude/worktrees/block-readers`:

    cargo bench --bench read_vcf -- /Users/jose/devel/popnei-bench/big.vcf \
        --threads 1 --runs 5
    cargo bench --bench read_vcf -- /Users/jose/devel/popnei-bench/big.vcf \
        --threads 18 --runs 5
    cargo bench --bench read_vcf -- /Users/jose/devel/popnei-bench/big.vcf.gz \
        --threads 1 --runs 5
    cargo bench --bench read_vcf -- /Users/jose/devel/popnei-bench/big.vcf.gz \
        --threads 18 --runs 5
    uptime            # before and after every set

The sweeps, each a script in `/Users/jose/devel/popnei-bench/` that runs the
settings one after another and then again, so that a drift of the machine does
not fall on one setting alone. `sweep_batch.sh` takes a binary, a file, a
number of threads, how many times to repeat, and the settings as
`lines:bytes`; `sweep_binaries.sh` takes a file, the threads, the repeats and
the builds; `sweep_buffer_final.sh` and `sweep_combined.sh` are the two of
them fixed to one question:

    sweep_batch.sh read_vcf_at_3c4bf5d big.vcf 18 3 \
        256:8388608 1024:8388608 2048:8388608 4096:8388608 4096:16777216
    sweep_binaries.sh big.vcf 18 3 read_vcf_buf8k read_vcf_buf64k \
        read_vcf_buf256k read_vcf_buf1024k
    sweep_buffer_final.sh big.vcf 18 3
    sweep_combined.sh big.vcf 18 3

The builds for the buffer of the file were made by setting the capacity with
`set_buffer.py`, also in that directory, and copying the benchmark aside:

    python3 set_buffer.py <the worktree> 256
    cargo bench --bench read_vcf --no-run
    cp target/release/deps/read_vcf-26e99fc6e5b48575 \
        /Users/jose/devel/popnei-bench/read_vcf_buf256k
    python3 set_buffer.py <the worktree> default   # and `git status --short` clean

The memory:

    uv run --no-project --with numpy python make_wide_vcf.py \
        /Users/jose/devel/popnei-bench/wide10000.vcf 10000 3000
    measure_memory.sh /Users/jose/devel/popnei-bench/big.vcf 18
    measure_memory.sh /Users/jose/devel/popnei-bench/wide10000.vcf 18
    measure_live_bytes.sh

`measure_memory.sh` runs `/usr/bin/time -l` on the benchmark with one run;
`measure_live_bytes.sh` runs `live_bytes/target/release/live_bytes`, the
scratch program with the counting allocator, built with `cargo build
--release` in `/Users/jose/devel/popnei-bench/live_bytes/`, which depends on
the core crate of the worktree by path and commits nothing.

The profiles:

    CARGO_TARGET_DIR=/Users/jose/devel/popnei-bench/target-profile \
        CARGO_PROFILE_BENCH_DEBUG=line-tables-only \
        cargo bench --bench read_vcf --no-run
    take_profile.sh /Users/jose/devel/popnei-bench/big.vcf 1 30 \
        /Users/jose/devel/popnei-bench/profile_1thread.txt
    take_profile.sh /Users/jose/devel/popnei-bench/big.vcf 18 30 \
        /Users/jose/devel/popnei-bench/profile_18threads.txt

`take_profile.sh` starts the benchmark with 400 runs and attaches
`/usr/bin/sample <pid> 30 1` to it.

The spike:

    rsync -a --exclude target --exclude dist \
        /Users/jose/devel/pynei/spike/pynei_spike/ \
        /Users/jose/devel/popnei-bench/spike_copy/
    uv venv --python 3.14 .venv
    VIRTUAL_ENV=.../spike_copy/.venv uv pip install numpy
    VIRTUAL_ENV=.../spike_copy/.venv uvx maturin develop --release \
        -m /Users/jose/devel/popnei-bench/spike_copy/Cargo.toml
    RAYON_NUM_THREADS=1 .venv/bin/python \
        /Users/jose/devel/popnei-bench/time_spike.py \
        /Users/jose/devel/popnei-bench/big.vcf 5000 5

## Where the numbers that were not measured here come from

- The reader before this plan, 1.24 s and 0.160 s plain, 1.58 s and 0.50 s
  bgzipped, and the spike of 20 September 2026, 0.54, 0.098, 0.84 and 0.40 s:
  "The measurement" of work package 5 of `docs/reports/vcf-to-blocks.md`. The
  targets of "Speed" of `docs/specs/io_vcf.md` are the spike's four plus a
  tenth.
- The 122-129 ms down to 105-107 ms that a 256 KiB buffer gave, and the
  13.3 MB and 19.2 MB of resident memory: the review of tasks 2.3 to 2.5, in
  the work report of this plan. This task's own measurements of the same
  things are above.
- plink2 v2.0.0-a.7.7 at 0.273 s and pyNei at 13.5 s on the plain file:
  `docs/specs/io_vcf.md`. Neither was run again here.
