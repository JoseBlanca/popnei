# The profiling environment: what runs on this project's machine

This file is inventory, not method: which profiling and benchmarking tools are verified to work where popnei is developed, and what this hardware cannot measure. The orchestrator reads it at step 2 of the review, and every reviewer reads it before writing a measurement plan, so that a plan only names tools that can be run. It was ported from the file of the same name in pop_var_caller, where the entries about macOS were verified on real reviews.

Entries are dated. Installs drift: check with `command -v` when a run is far from the date.

## The machine

An Apple M5 Pro, 18 cores, 64 GB, macOS, native `aarch64-apple-darwin`. Verified 2026-09-20. popnei will also have a Linux development container, which does not exist yet; when it does, its tools are listed here, with the two things a container on this Mac cannot give: a native build against Accelerate, and timings of the machine the user has. Timings that are compared with pyNei, plink2 or GMMAT are taken on the host.

## Verified present, 2026-09-20

- **`samply`** (`~/.cargo/bin/samply`): the sampling profiler. `samply record --save-only -o out.json.gz -r 2000 -- <binary> ...`. `--save-only` writes unsymbolicated addresses, and the symbols are resolved in the Firefox profiler UI. For a self-time ranking without a browser use `sample`.
- **`sample`** (`/usr/bin/sample`, part of macOS): attach to a running process, `sample <pid> <seconds> -file out.txt`. It gives a symbolicated call tree directly. The threads of the rayon pool show up as idle `Sleep::sleep` frames with the full sample count: filter them out. It also works on a Python process that is inside a popnei call, which is the way to see whether the time is in Rust, in numpy or in pandas.
- **`xctrace`** (`/usr/bin/xctrace`, Instruments): `xcrun xctrace record --template 'Time Profiler' --launch -- <binary>`. Instruments reads the hardware counters of the M series and the state of the threads, on and off the CPU, and nothing else here does. Its `.trace` output is awkward to read without the GUI, so use `sample` and `samply` first.
- **`cargo asm`** (`cargo-show-asm`, `~/.cargo/bin/cargo-asm`): to see whether a named hot function was vectorized and whether its bounds checks were dropped, `cargo asm --lib --simplify "popnei::path::to::function"`. A helper that was fully inlined is looked at in its caller.
- **`hyperfine`** (`/opt/homebrew/bin/hyperfine`): wall time of a command with its spread over runs. The tool for comparing with plink2 and with a pyNei script, which are commands.
- **`dtrace`** (`/usr/sbin/dtrace`): needs `sudo`; not needed when `sample` and `samply` work.
- **`plink2`** (`/opt/homebrew/bin/plink2`, v2.0.0-a.7.7), **`Rscript`** (`/opt/homebrew/bin/Rscript`) with GMMAT 1.5.0 and rrBLUP 4.6.3, which pyNei's reference script uses, and **pyNei** at `/Users/jose/devel/pynei`, run with `uv run python` from there.
- **`node`** (`/opt/homebrew/bin/node`), for pyodide under node.

## Not installed, 2026-09-20

`cargo-flamegraph`, `valgrind`, `perf`, `coz`, `cargo-criterion`, `maturin` as a global command (it is run with `uvx maturin` or from the project's environment). `perf` and `valgrind` do not exist for macOS on Apple silicon. criterion is a dev-dependency of the crate, added with the first benchmark, and DHAT is the `dhat` crate behind a cargo feature, added with the first allocation question.

## What this machine cannot measure

Without Instruments there are no hardware counters: no cache misses, no branch misses, no instructions per cycle. The stand-ins:

- **Branch misprediction**: run the same loop over the same values sorted and shuffled. The bytes and the work are the same and only the predictability differs, so the gap is the cost of the mispredictions. See `checklists/hot_loops.md`.
- **Allocations**: DHAT counts them exactly, and the count is the same on every run.
- **Time off the CPU**, a thread waiting on a channel, a lock or the interpreter: a sampling profile does not see it. Time the blocked `send` and `recv` by hand with `Instant`, a few lines, or use Instruments' System Trace.
- **Cache behaviour**: no direct measure. Change the size of the working set, the block size, and see where the time per genotype jumps.

## Diagnosing a sampling profiler that does not work

On macOS profiling your own process works without special permissions. What can block it: System Integrity Protection restricts some kernel probes of `dtrace`, which does not affect `sample` or `samply` on your own binary; and a release binary without debug info gives a profile with no function names, which is fixed in `Cargo.toml` with `debug = true` or `"line-tables-only"` under `[profile.release]` or a profile of its own for profiling. On Linux, `kernel.perf_event_paranoid` of 3 or more blocks `perf_event_open` for unprivileged users, and in a rootless container `--cap-add=PERFMON` does not help: lower the sysctl on the host, `sudo sysctl kernel.perf_event_paranoid=1`.

When no sampling profile can be taken, the review says so first, as `SKILL.md` describes, and does not work around it in silence.

## Datasets

The datasets a timing is stated on, from `docs/rust_core.md`: 100000 variants x 1000 samples, a 400 MB VCF and its 53 MB gzip, and the matrices of 1000 to 5000 samples; section 6 of that document says how they were made. pyNei's reference panel, `test/gwas_reference/sim_missing.vars`, 1200 variants x 200 samples, is for correctness and too small to time anything. The objectives' largest dataset is a million variants of 10000 samples.

## The toolchain

`rustc` and `cargo` 1.98.0, not pinned yet. The decisions of the autovectorizer change between versions of rustc with no change in the code, so when popnei pins its toolchain, a bump of the pin makes the saved baselines stale: run them again.
