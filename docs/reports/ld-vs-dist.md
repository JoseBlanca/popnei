# Report: the curve of r² against distance, per population

24 September 2026. It records how `docs/plans/ld-vs-dist.md` was carried
out, on the branch `plan/ld-vs-dist` in the worktree
`.claude/worktrees/ld-vs-dist`. The plan builds the item "LD against
distance, per population" of `docs/specs/ld.md`: for each population of a
dataset, how the r² of a pair of variants falls off as the two move apart
along a chromosome, in bins of distance, as a curve fitted to the pairs,
and as the one distance at which r² has fallen to half.

**The plan is under way.** Nothing is merged into `main` and nothing is
pushed.

## What was in place before the first task

The plan's "What has to be in place" was checked by running it, on the
branch at `a04c7eb`, which is `main` at `cce5801` with the two spec
commits and the two plan commits of this branch on top.

| what the plan says | the command | what it gave |
|---|---|---|
| the core and linear algebra suites pass | `cargo test --workspace` | `787 passed; 0 failed; 2 ignored` and `149 passed; 0 failed` |
| 49 ld tests, none of the fall-off | `cargo test -p popnei --lib ld:: -- --list` | `49 tests`, and `0 tests` for both `ld::dist` and `ld::decay` |
| the Python tests of the item do not exist | `uv run pytest tests/test_ld.py -k ld_and_dist` | exit 5, `16 deselected`, which is what pytest gives when `-k` matched nothing |
| no code of the item anywhere | `grep -rn "ld_and_dist\|LdAndDist\|LdBins\|LdDecay" crates python js tests` | nothing |
| plink2 and R are on the machine | `plink2 --version`, `R --version` | `PLINK v2.0.0-a.7.7 M1 (18 Sep 2026)`, `R version 4.6.1 (2026-06-24)` |
| the reference script writes its files again and finds them the same | `tests/reference/ld/run_plink2.sh` into an empty directory | exit 0, nothing of its own printed |

The other checks of the `coding` skill were run at the same commit and are
clean: `cargo fmt --all --check`, `cargo clippy --workspace --all-targets
-- -D warnings`, `cargo test -p popnei --no-default-features` with the same
787, `cargo wasm-check`, `uv run ruff format --check` and `uv run ruff
check`, and `npm run build && npm test` in `js/popnei` with 325 tests
passing.

Two things the plan states that the commands corrected:

- **The npm package had to be built before the Python tests would run.**
  A fresh worktree has no compiled extension module, so `uv run pytest`
  failed to import `popnei` until `uv run maturin develop` had been run
  once. It is a step of the `coding` skill's checks and not a change to
  the plan.
- **`js/popnei/test/ld.test.ts` holds 10 tests and not the 11 the plan
  counts.** `node --test js/popnei/test/ld.test.ts` prints `tests 10`, and
  the file has 10 `test(` at the top level. Deliverable 6 of work package 1
  and deliverable 7 of work package 2 ask that the file then run more ld
  tests than it runs today, so the number to beat is 10.

## Work package 1: the bins of r² against distance

### The tasks as they were done

**1.1, the window over the blocks.** Commit `a0a8181`. It put
`TheWindowOfTheBlocks` in a new `crates/popnei/src/ld/dist.rs`, declared
by one `mod dist;` line added to `crates/popnei/src/ld.rs`. The window
takes a block, holds the blocks whose variants are within `max_dist` of
the newest variant read and on its chromosome, and gives back how many of
the oldest blocks fell out on that call, which is what lets the next task
drop what it keeps beside each block in step.

The plan puts this code in `crates/popnei/src/ld.rs`. It went in
`crates/popnei/src/ld/dist.rs` instead, which the plan's own check asks
for: the task's tests have to be printed by `cargo test -p popnei --lib
ld::dist -- --list`, and that needs a module `ld::dist`. The crate is
edition 2024, so a file beside `ld.rs` is a submodule of it, and the 3859
lines already in `ld.rs` were left alone.

`cargo test -p popnei --lib ld::dist -- --list` prints `14 tests`, where
it printed `0 tests`. `cargo test --workspace` gives 801 passed with 2
ignored, which is the 787 of the starting state plus those 14, and 149 in
the linear algebra crate. `cargo test -p popnei --no-default-features`
gives the same 801 on the faer backend. `cargo fmt --all --check`, `cargo
clippy --workspace --all-targets -- -D warnings`, `cargo wasm-check`,
`uv run ruff format --check` and `uv run ruff check` are clean, and `uv
run pytest` gives 499 passed. Every one of these was run by the
orchestrator and not taken from the task's report.

**1.2, the dosages of each population over the window.** Commit
`175818f`, in the same `crates/popnei/src/ld/dist.rs`. For each
population it keeps the variants of the held blocks that pass that
population's major allele frequency and one set of dosage matrices over
them, rebuilt at each block, so that a tile spanning a block boundary is
one call.

It added two of the errors that the spec's `# Errors` of
`calc_ld_and_dist` already names, a `max_allowed_maf` that is NaN or not
from 0 to 1, and a population with no individual. Without the first a NaN
threshold empties every population without saying so, and without the
second an empty population is read as every individual. Both are a
`ValueError` in Python, which is where the spec puts an argument out of
range and a population with no individual. No spec change was needed.

`cargo test -p popnei --lib ld::dist -- --list` prints `25 tests`, where
it printed 14. `cargo test --workspace` gives 812 passed with 2 ignored,
and `--no-default-features` the same 812 on faer; 149 in the linear
algebra crate, 499 pytest, and fmt, clippy, wasm-check and ruff clean.
Every one was run by the orchestrator.

**One number for the performance review that follows this plan.** Each
block builds two sets of dosages for each population, one over the block
to work out the frequencies and one over the whole window, so a pass at
blocks of 7 variants with a window of 250 variants reads the genotypes of
a variant about 36 times over. It was left alone on purpose: the owner
decided on 24 September 2026 that the speed of this pass is a review of
its own.

**1.3, the pairs of the window in tiles and the bins.** Commits `71e69bd`
and `7077666`, the first of them a correction to the spec that the second
builds. It put `LdAndDistOptions`, `calc_ld_and_dist`, `LdAndDist` and
`LdBins` in `crates/popnei/src/ld/dist.rs`, and added the last two errors
the work package asks for, a `min_dist` above `max_dist` and a `num_bins`
of 0. `LdBins::decay()` was deliberately left out; it is work package
2's.

`cargo test -p popnei --lib ld::dist -- --list` prints `42 tests`, where
it printed 25. `cargo test --workspace` gives 829 passed with 2 ignored,
and `--no-default-features` the same 829 on faer; 149 linear algebra, 499
pytest, and fmt, clippy, wasm-check and ruff clean. All seven errors that
a `u64` can reach have a test, and so have the two cases that are empty
and not errors. The eighth error the work package lists, a negative
`min_dist`, cannot be a `u64` and is refused in the Python layer, which
is task 1.5.

### Two sentences of the spec that the work corrected

The spec's "How it runs" promises that the bins are the same to the bit
whatever the size of the blocks. Two of its own rules broke that promise,
and deliverable 3 is what found them. Both were rewritten in `71e69bd`,
before the code that follows them. Neither changes a value a user sees
beyond the last bits, neither touches the public API, and neither
answers an open point, so the orchestrator took them; the owner is told
here because the spec now reads differently.

**The bins are added up in the order of the variants, not in the order of
the tile pairs.** The spec had each tile pair give the count and the two
sums of every bin it touched. A block of the reader can end inside a tile
of the newest variants, and that tile then reaches a bin in two goes
where one block sends it in one; a sum of floats moves in its last bits
when the order of its terms moves. Counting one newest variant at a time,
with the variants that pair with it in the order of the pass, is cut by
neither a block nor a tile, so the result no longer moves with the tile
size either.

**A block that has fallen out of the window is dropped one block later.**
The window reaches back from the newest variant read, which is the last
variant of the block just taken. A variant within `max_dist` of the
*first* variant of that block can be further than `max_dist` from its
last, so dropping it at that step leaves real pairs uncounted, and how
many depends on the size of the blocks: the task measured 272 pairs in
the last bin at blocks of 3 where blocks of 500 give 328. What the
deferral costs is holding the blocks that fell out for one more step.

The test that holds both is
`the_bins_do_not_move_with_the_blocks_nor_with_the_tiles`, which runs 300
variants of six individuals in two populations at blocks of 1, 2, 3, 7,
64 and 500 and tiles of 1, 2, 7 and 256, twenty-four combinations, and
compares the counts of pairs exactly and every mean and standard
deviation by its bits. It first asserts that every bin holds a pair, so
that it is a test that can fail.

**1.4, the three tables of bins against plink2.** Commit `b870a9b`. It
stored `tests/reference/ld/ld.bins.txt`, 36 lines of what
`docs/reports/ld-method/bins.py` prints for the three populations, made
`tests/reference/ld/run_plink2.sh` write it again and compare it, and put
two cargo tests in `crates/popnei/src/ld/dist.rs` that assert the three
tables at blocks of 7, 64 and 500 variants and at rayon with 1 thread and
with 4.

The counts of pairs and the means of all three tables, and the standard
deviations of the first, are the literals of the spec's own tables. The
standard deviations of `pop_a` and `pop_b`, which the spec leaves out of
its table to keep it readable, come from `ld.bins.txt`. Every number the
spec prints was checked digit for digit against that file and they agree,
so nothing had to be reported as a disagreement.

The thread test builds a rayon pool of the size it wants and asserts
`rayon::current_num_threads()` inside it before the passes, so both pools
really had 1 and 4. The threads reach the work: the VCF reader parses the
lines of a batch on the pool the caller is in, and with the `blas` feature
off the products go through faer on rayon.

`cargo test -p popnei --lib ld::dist -- --list` prints `44 tests`, where
it printed 42. `cargo test --workspace` gives 831 passed with 2 ignored,
and `--no-default-features` the same 831; 149 linear algebra, 499 pytest,
and fmt, clippy, wasm-check and ruff clean. The three-table test takes
0.04 s and the thread test 0.01 s of the core crate's 1.31 s.

### Deliverables 1 and 3 of work package 1, checked by the orchestrator

- **Deliverable 1**, `ld.bins.txt` written again and compared:
  `tests/reference/ld/run_plink2.sh` into an empty directory exits 0 and
  prints nothing of its own, and `cmp` finds the stored file identical to
  the one the script wrote from plink2's matrix.
- **Deliverable 3**, the three tables at three block sizes and two thread
  counts: `cargo test -p popnei --lib ld::dist -- --list` prints 44 tests
  where it printed 0 at the start of the plan, and `cargo test
  --workspace` runs them at 831 passed.

**1.5, the Python layer of the bins.** Commit `1b414c3`. It put
`calc_ld_and_dist_per_pop` and `LdAndDistPerPop` in `python/popnei/ld.py`
and the pyo3 function in `crates/popnei-python/src/ld.rs`, which lends
the reader chain to the core and takes the variants of the pass from the
core's own count rather than counting them again, and added five pytest
tests.

The Python layer refuses the arguments that cannot reach the core: a
`min_dist` or a `max_dist` below 0 and a `num_bins` below 0 are each a
`ValueError` naming the argument and its value, and a fractional or
boolean value for any of them is a `TypeError`. That is the eighth error
of the work package, the one the core's `u64` cannot see. The seven the
core already raises are not raised twice.

The four defaults became `pub const` of the core, `DEFAULT_MIN_DIST`,
`DEFAULT_MAX_DIST`, `DEFAULT_NUM_DIST_BINS` and
`DEFAULT_MAX_ALLOWED_MAF`, which is how every other Python default of
popnei is fed. The bins' constant is `DEFAULT_NUM_DIST_BINS` and not
`DEFAULT_NUM_BINS`, which the histogram already has.

The test that compares with pyNei does ask pyNei: it asserts that
`num_vars_per_pop` equals both the 396 and 402 of the spec and what
pyNei gives when `filter_samples` is put around the individuals of each
population and `filter_by_maf` over what is left. The ten rows of the
first table are asserted within 1e-12 relative, and the widest gap
measured against the spec is 9.2e-16 relative.

`uv run pytest tests/test_ld.py -k ld_and_dist` gives `5 passed, 16
deselected`, where it exited 5 with nothing matched at the start of the
plan, which is deliverable 5's check. The whole pytest suite is 504 where
it was 499; the core crate stays at 831 with 2 ignored, the same on faer,
149 linear algebra, and fmt, clippy, wasm-check and ruff clean.

**1.6, the TypeScript layer of the bins.** Commit `6fb3343`. It put
`calcLdAndDistPerPop` in `js/popnei/src/ld.ts`, exported from both
`node.ts` and `web.ts`, over a wasm-bindgen function in
`crates/popnei-js/`. The five arrays of a population cross the boundary
flat and the package cuts them. The package refuses a `minDist` or a
`maxDist` that is not a whole base pair of 0 or more, a `numBins` that is
not a whole number of 0 or more, and a `maxAllowedMaf` that is no number.

`node --test js/popnei/test/ld.test.ts` runs 17 tests, where it ran 10,
and the whole node suite is 332 where it was 325.

### The six deliverables of work package 1, each checked by the orchestrator

| the deliverable | the command | what it gave |
|---|---|---|
| 1, `ld.bins.txt` written again and compared | `tests/reference/ld/run_plink2.sh` into an empty directory | exit 0, no file named as differing, and `cmp` finds the stored copy identical to the written one |
| 2, the worked example | `cargo test -p popnei --lib ld::dist -- --list` | `44 tests`, where it printed `0 tests` before the plan |
| 3, the three tables at three block sizes and two thread counts | `cargo test --workspace` | `831 passed; 0 failed; 2 ignored`, and `149 passed` in the linear algebra crate |
| 4, the eight errors and the two empty cases | the same list and run | seven errors as cargo tests, the eighth in the two packages that can see it |
| 5, the Python layer | `uv run pytest tests/test_ld.py -k ld_and_dist` | `5 passed, 16 deselected`, where it exited 5 with nothing matched |
| 6, the TypeScript layer | `npm run build && npm test` in `js/popnei` | `tests 332, pass 332, fail 0`, and 17 in the ld file where there were 10 |

`cargo test -p popnei --no-default-features` gives the same 831 on the
faer backend, and `cargo fmt --all --check`, `cargo clippy --workspace
--all-targets -- -D warnings`, `cargo wasm-check`, `uv run ruff format
--check` and `uv run ruff check` are clean.

**One thing found on the way that is nobody's task here.**
`js/popnei/README.md` lists five calculations and leaves out
`calcPopDists`, `calcKinship` and `calcGwas`, so it was already stale
before this plan and this calculation was not added to it either. It is
worth a commit of its own, outside this plan.
