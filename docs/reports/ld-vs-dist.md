# Report: the curve of r² against distance, per population

24 September 2026. It records how `docs/plans/ld-vs-dist.md` was carried
out, on the branch `plan/ld-vs-dist` in the worktree
`.claude/worktrees/ld-vs-dist`. The plan builds the item "LD against
distance, per population" of `docs/specs/ld.md`: for each population of a
dataset, how the r² of a pair of variants falls off as the two move apart
along a chromosome, in bins of distance, as a curve fitted to the pairs,
and as the one distance at which r² has fallen to half.

**The plan is done and the branch is not merged.** Every task is ticked,
every deliverable of the two work packages was checked by the orchestrator
running its command, both work packages were reviewed, and the plan's own
final check passes. Nothing is merged into `main` and nothing is pushed.

## What exists now that did not

A user calls `calc_ld_and_dist_per_pop` in Python or `calcLdAndDistPerPop`
in TypeScript over a dataset and a set of populations, and gets back, for
each population on its own: one row per bin of distance with how many
pairs of variants it holds, their mean r² and its standard deviation; how
many variants that population kept at its own major allele frequency; and
the curve fitted to every one of its pairs, as the scaled recombination
per base pair, the curve's value at distance 0, and **the distance at
which r² has fallen to half of that**, which is the one number a plot of a
population is labelled with and the number the web application was waiting
for.

The bins come out the same to the bit whatever the size of the blocks the
reader gives, the size of the tiles the pairs are computed in, and the
number of threads. The bins agree with plink2 v2.0.0-a.7.7 on three
populations of the reference dataset, exactly on the counts of pairs and
within 1e-12 relative on the means and the standard deviations. The fitted
curve agrees with R 4.6.1 to 1.03e-9 relative on the first of those
populations, where the spec asks for 1e-6.

At the commit this branch starts from: 787 tests in the core crate, 149 in
the linear algebra crate, 499 pytest, 325 node, and none of them of the
fall-off against distance. Now: 865 in the core crate with 2 ignored, the
same 865 on the faer backend, 149 in the linear algebra crate, 509 pytest,
334 node. `cargo fmt`, `cargo clippy` with every target and warnings
denied, `cargo wasm-check` and ruff are clean, and
`tests/reference/ld/run_plink2.sh` writes its ten files into an empty
directory and finds every one the same as the copy stored beside it.

## What is asked of the owner

**The merge.** Nothing has been merged into `main` and nothing pushed.

**Four decisions**, none of which blocks the merge, because the branch
builds a defined behaviour for each and says here what it is. Three of
them are one question in three parts, and the fourth reaches beyond this
module.

**Decision 1: what popnei refuses at the edges of the bin arguments.**
Three things go wrong at once out there, and one answer settles all three.
A `max_dist` above 9223372036854775807 raises a `RuntimeError` reading
"one count of the pass is 9223372036854775808, more than a result holds",
where popnei's own convention reserves a `RuntimeError` for a defect of
popnei and the value is a distance and not a count. A `num_bins` above the
number of distances in the range gives bins whose stated range runs
backwards, at `min_dist` 1, `max_dist` 3 and `num_bins` 5 the rows read
"3 to 2" and "4 to 3", and a pandas frame whose index repeats a label;
it is reachable at the default `num_bins` of 50 with any `max_dist` below
50. And nothing bounds the work either argument asks for: a `num_bins` of
10⁶ on a five-variant dataset returns a million rows in 0.3 s and 10⁸
returns a hundred million in 25.6 s, the cost growing with the argument
and not with the data. The options are to refuse both, each a `ValueError`
naming the argument and its value, which costs two sentences of the spec
and two public refusals; to refuse the `max_dist` and write into the spec
what a bin holding no whole distance reports; or to define both in the
spec and refuse neither, which needs the distance columns to stay unsigned
and brings back the wrap that gave 18446744073709550687 for 7815 − 8744.
Recommendation: refuse both. Neither case is one a real dataset reaches,
so refusing costs a user nothing they wanted.

**Decision 2: the memory guard does not guard.** `docs/specs/ld.md` says
the memory of the counts per distance is "asked of the machine with
`try_reserve_exact` before the pass and refused rather than taken", the
`coding` skill states the principle as "a matrix this machine cannot hold
is an error and not a process that ends", and deliverable 5 of work
package 2 has a test asserting it. On this machine, which has 69 GB, a
`max_dist` of 10⁹ asks for 16 GB and is granted it, which is right; a
`max_dist` of 10¹² asks for 16 TB and **the process is killed**, signal 9,
rather than an error being raised. macOS hands out the virtual allocation
and kills the process when the pages are touched, so `try_reserve_exact`
never returns the failure the guard is built on. This is not this module's
alone: `calc_r2_matrix` guards its matrix the same way and
`docs/specs/linalg.md` asks for it for the eigendecomposition workspace.
Nothing on this branch touches it, because the answer belongs to popnei as
a whole. What is asked is whether the guard should be backed by a check
against the memory the machine reports, or the specs' promise weakened to
what `try_reserve_exact` can deliver.

**Decision 3: a calculation cannot be abandoned.** A `max_dist` of 10¹⁰ on
a 21-variant file runs for 190 s, and Ctrl-C does not stop it. This is not
a fault of this plan: no calculation of popnei raises a Ctrl-C before it
returns, `calc_kosman_sums` and the r² matrix included, because the core
crate has no in-calculation check and adding one changes the signature of
every calculation that takes it. It is recorded here because this
calculation is the first whose running time a user can raise without
raising the size of their data.

**Decision 4: half the curve at one and at two individuals.** At those two
sample sizes the curve never falls to half of its value at distance 0,
because the correction for a finite sample holds it up at 1 over the
individuals rather than letting it fall to 0. The branch gives `half_dist`
alone as NaN there and keeps the other two, which really were fitted. The
alternative is three NaN, a two-line change. `calc_ld_and_dist` reaches
neither sample size, so only a caller of `fit_ld_decay` with a table of
its own meets it. Recommendation: keep one NaN, which the `api` and
`spec` reviewers both independently preferred for the same reason.

**One thing is not asked**, so that nobody waits on a decision that is not
wanted: the speed of the pass and of the fit. The owner decided on 24
September 2026 that it goes to a performance review of its own after this
plan, and nothing here was traded for it. The numbers that review will
want are gathered under "For the performance review" at the end.

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
popnei is fed. The bins' constant was `DEFAULT_NUM_BINS` in the core and
`DEFAULT_NUM_DIST_BINS` in both bindings, which the review found and the
fixes made one name in the three layers.

The test that compares with pyNei does ask pyNei: it asserts that
`num_vars_per_pop` equals both the 396 and 402 of the spec and what
pyNei gives when `filter_samples` is put around the individuals of each
population and `filter_by_maf` over what is left. The ten rows of the
first table are asserted within 1e-12 relative, and the widest gap
measured against the spec is 1.25e-15 relative, at the standard deviation
of the eighth bin. The task that wrote the test, and its commit message
`1b414c3`, gave 9.2e-16, which is the gap at the first bin and not the
widest; the review measured the rest.

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

### What the review of work package 1 found

Seven reviewers read the work package at `6fb3343`, one for each category
the `code-review` skill lists, each with a fresh context. `spec` and
`tests` worked in worktrees of their own because they build and break
code; the other five read the shared checkout. The orchestrator judged
every finding and reproduced the two worst itself before acting.

**Three findings were a wrong number a user would see.** None of them
could have been caught by the three reference tables, which is why the
work package had passed every check.

1. **The bin a pair falls in disagreed with the range that bin reports.**
   The bin index divided by a floating point width and the bounds
   multiplied by it, so the two rounded opposite ways. At `min_dist` 1,
   `max_dist` 18 and `num_bins` 14 the bin labelled 9 to 9 held the pair
   10 bp apart and the bin labelled 10 to 11 held nothing. A reviewer
   measured it in 1.9 per 100 of random settings. The reference tables
   never meet it because 250000/10 and 1000000/50 divide exactly. Both
   are now worked out in whole numbers over `u128`, so they agree by
   construction. Commit `9379058`.
2. **A pair was lost, and the answer moved with the tile size, when a
   chromosome came back after another one.** The rows a tile reached back
   over were taken from its first column alone, so a variant a later
   column still paired with was never looked at. On 255 variants of
   `chr3` followed by `chr2:10`, `chr1:100` and `chr2:50` the bins counted
   20450 pairs where the matrix holds 20451, and with 254 variants of
   `chr3`, which moves the tile boundary, both give 20351. Nothing
   refuses a source whose chromosomes are not grouped, so it was
   reachable. The tiling was repaired rather than the input refused,
   since a new refusal would be the owner's. Commit `259e3db`.
3. **The counts reached pandas unsigned**, so subtracting one bin's pairs
   from another's wrapped: 7815 − 8744 came back as
   18446744073709550687. The Python crate already had this rule written
   down with this reasoning, and two other modules already followed it.
   The counts and the two distance columns are now `int64`. Commit
   `7d2ad20`.

**Two tests could not have failed.** The one guarding the order the bins
are added up in used six individuals, so every r² was a multiple of a
quarter and every sum exact whatever the order; two different
reorderings left it green. And no fixture anywhere used a `min_dist` of 0
or put two variants at one position, so a mutation pairing every variant
with itself passed all 44 tests. The first fixture now has 24
individuals and a companion test pins the three reference tables at four
tile sizes; the second has a SNP and an indel at one position, giving 3
pairs at `min_dist` 0 and 2 at 1. Commits `59480b9` and `9f6f86d`.

**Each population held the window's genotypes over every individual of
the source**, not over its own, so the genotypes sat in memory once per
population. Measured on 3000 variants of 500 individuals with `max_dist`
100000 and populations of one individual each, an extra population cost
3.56 MB, which is the row of the whole source, where the three matrices
of a population of one individual are 0.072 MB. Each population now
keeps only its own individuals' alleles, and the window itself keeps the
chromosome and the position and no genotypes at all. On the same shape an
extra population now costs 0.13 MB, measured by the orchestrator: peak
resident memory 119.2 MB at one population and 121.7 MB at twenty.
Commit `9d30a0f`, with the spec's "How it runs" rewritten in `8ff74cb`
to say where the genotypes live.

**Three error messages misled.** A `numBins` of 0 told a TypeScript user
about `min_dist` and `max_dist`, which their API has not (`852089e`). The
Python refusal of `min_dist` named a floor of 1 where 0 is accepted and
TypeScript said 0, so the two layers stated different limits
(`2f040d5`). And the refusal of memory said "values of 8 bytes" for
genotypes that are one byte and for structures of 16 and 104, and two
places reported the length already held rather than the one asked for
(`9505404`).

**One finding was not taken.** A reviewer reported that the WebAssembly
build gives a first mean of 0.20767885551844037 where the spec's literal
is 0.20767885551844031, and asked whether anything says that is
allowed. The `coding` skill says it plainly: popnei does not promise the
same bits on every platform, and results agree with the reference
programs within the tolerance of the spec. The gap is 2 units in the last
place against a tolerance of 1e-12.

**The same 2 units in the last place turned out to be the native build's
too, and always had been.** The reviewer had assumed the native build
matched the literal rather than running it. Both commits were built side
by side and all three tables are byte for byte identical before and after
the fixes, so nothing moved: popnei adds the r² in the order of the pass
and the literal is numpy's sum over plink2's matrix in
`docs/reports/ld-method/bins.py`.

After the ten fixes: `cargo test --workspace` gives 837 passed with 2
ignored and 149 in the linear algebra crate, the same 837 on the faer
backend, 50 tests under `ld::dist`, 505 pytest, 332 node, fmt, clippy,
wasm-check and ruff clean, and `run_plink2.sh` into an empty directory
exits 0 naming no differing file. Every one of these was run by the
orchestrator.

### The rest of the findings, and one belief the work disproved

Nine smaller findings were fixed in commits `5ffd4fd` to `3739562`: the
one default that had two names across the three layers; the module
documentation of both `ld` files, which described what the file held
before the later tasks filled it; a check in the TypeScript package
whose doc comment listed its callers and had fallen two behind; three
comments that no longer said what the code does, one of them miscounting
the worked example; the two stale notes about which errors a Python user
can reach; and two statements of the spec that had no test, an individual
in more than one population or in none, and a dataset of variants each on
its own chromosome, whose fixture read one variant a block and so never
made a pair at all.

Two of them were judged and not simply applied.

**The count of the variants of a pass keeps the plain sum.** The review
asked for `checked_add`, which the `coding` skill prefers. The bound is
real: the count would have to pass 1.8·10¹⁹ variants, whose positions
alone are 147 exabytes, so the skill's own rule for a bound that cannot
be reached applies and the operator stays with the bound written beside
it. `checked_add` would need a new case in the public error enum and in
the spec for a pass no machine can run, which is an interface change and
so the owner's, not a review fix. `crates/popnei/src/dists.rs:824`
counts its own pass the same way, so the two modules agree.

**Cutting the tiles from the first variant of the pass buys nothing that
can be measured.** The plan's "What could go wrong" says that a tiling
which follows the blocks instead "passes deliverable 2 and fails
deliverable 3 at the second block size". Two reviewers and then the fix
tried it: with the alignment removed, so that the tiles are cut from the
start of the window rather than from the first variant of the pass, all
190 tests of the `ld` module still pass, and the three reference tables
are identical at tiles of 7, 64, 256 and 500. What actually buys the
promise that the answer does not move is the correction of `71e69bd`,
adding the bins up one newest variant at a time in the order of the pass.
The alignment is kept, because it costs nothing and keeps the buffers a
fixed size, and its doc comment now says that it is not what fixes the
numbers. Nobody has found a case where it changes one.

After the nineteen fixes of the three rounds: `cargo test --workspace`
gives 837 passed with 2 ignored and 149 in the linear algebra crate, the
same 837 on faer, 50 under `ld::dist`, 506 pytest with 7 of them under
`-k ld_and_dist`, 332 node, fmt, clippy, wasm-check and ruff clean, and
`run_plink2.sh` into an empty directory exits 0 naming no differing
file. Every one was run by the orchestrator.

**Three things found on the way that are nobody's task in this plan.**
`cargo doc -p popnei` fails on an unresolved link to
`Error::ReaderGaveNoVariants` at `crates/popnei/src/ld.rs:973`, which was
already broken before this work package and is not among the checks the
`coding` skill lists. The chromosome comparison in
`the_first_row_in_reach` is an optimisation that no test covers: removing
it breaks nothing. And `js/popnei/README.md` lists five calculations and
leaves out four of them.

## Work package 2: the fitted curve and the half distance

**2.1, the pairs counted at each distance.** Commit `05fefd1`, in
`crates/popnei/src/ld/dist.rs`. Beside each bin the pass now keeps, for
each population, a count of pairs and a sum of r² for every distance from
`min_dist` to `max_dist`, 16 bytes a distance, asked of the machine with
`try_reserve_exact` before the first block is read and refused rather
than taken. A pair is added to its distance in the same call that adds it
to its bin, so the two sums have the same order. When the pass ends the
range is compacted to the distances that hold a pair, which is what the
fit reads.

That is deliverable 5 of the work package. Its test asserts that a range
of 2^64 distances is refused before a block is read, naming what could
not be held, how many values and how many bytes one of them is; replacing
the checked allocation with an ordinary vector makes it fail.

The counts per distance were checked against the bins over the same pass:
the counts of the distances inside a bin's bounds add up to that bin's
`num_pairs` exactly, and their sum of r² gives that bin's mean within
1e-12 relative, the two sums differing only in their order. The task also
added the distances to the table that work package 1's block, tile and
thread tests compare, so they are now held to the same bit equality as
the bins.

`cargo test -p popnei --lib ld::dist -- --list` prints `55 tests`, where
it printed 50. `cargo test --workspace` gives 842 passed with 2 ignored,
the same 842 on faer, 149 linear algebra, 506 pytest, and fmt, clippy,
wasm-check and ruff clean. `run_plink2.sh` into an empty directory exits
0 naming no differing file, so the three tables of work package 1 are
untouched. Every one was run by the orchestrator.

**2.2, the fit and the half distance.** Commits `e5c7625`, a correction
to the spec, and `c12dc35`, the code. A new module
`crates/popnei/src/ld/decay.rs` holds `fit_ld_decay` and `LdDecay`: the
curve, the sum it makes smallest, the grid of 141 values, the golden
section search and the bisection for the half distance. Thirteen tests
under `ld::decay`, where the plan started at none.

**The test that makes a silently wrong fit impossible passes exactly.**
Its r² are taken from the curve itself at a ρ per base pair of 0.0001 and
100 individuals, so the answer is known before the fit runs: the fit
gives 0.0001 back with a relative error of 0, a half distance of
21608.13587253306 against the spec's 21608.135872529165, which is
1.8e-13 relative where 1e-6 is asked, and the r² at distance 0 to the
bit. A second table at a ρ per base pair of 2e-4, which is not one of the
141 values of the grid, is what exercises the search: the grid alone
lands 2.4e-3 away and the search brings it to 3.2e-8.

Both ends of the searched range are tested, which the plan calls the case
easiest to miss: an r² of 1 at every distance, above the curve's ceiling
everywhere, puts the smallest at the first grid point, and an r² of 0 at
every distance puts it at the last. Both give three NaN and not an error,
so neither returns the half distance of about 10⁻¹² or 10¹² base pairs
that an unguarded search would.

Only one function that is not an addition, a subtraction, a
multiplication or a division is called: `powf`, to move along the
exponent of 10 in the grid and the search, 183 times a fit. The curve and
the sum use the four operations alone, so the platform cannot move a
returned value except through where the search looked.

### A sentence of the spec that is wrong for small samples

"The curve that is fitted" says "The curve falls without turning, so
exactly one ρ gives half of what it gives at ρ of 0". That is true of
every sample size a dataset reaches and false at 1 and 2 individuals.
The curve does not fall to 0: the correction for a finite sample holds it
up at 1 over the individuals. The orchestrator checked the arithmetic:

| individuals | at ρ of 0 | half of that | the floor, 1/n | does it reach half |
|---|---|---|---|---|
| 1 | 1.198347 | 0.599174 | 1.000000 | never |
| 2 | 0.826446 | 0.413223 | 0.500000 | never |
| 3 | 0.702479 | 0.351240 | 0.333333 | yes |
| 50 | 0.469421 | 0.234711 | 0.020000 | yes |
| 100 | 0.461983 | 0.230992 | 0.010000 | yes |

At 1 and 2 individuals the bisection would have run to the top of its
range and returned 10⁶ divided by the fitted ρ per base pair, which is a
number and not a NaN. The half distance alone is now NaN there, written
into the spec in `e5c7625` before the code, with the other two kept
because the fit itself succeeded.

**What the owner may want to change.** `calc_ld_and_dist` reaches
neither case, so only a caller of `fit_ld_decay` with a table of its own
sees it: one individual has no variant with variance and counts no pair,
and two individuals give every pair an r² of 1, whose smallest falls at
the bottom end of the range and is already the three NaN of "The cases".
If the owner would rather have three NaN there than one, it is a two-line
change. The orchestrator took the one-NaN form because the fit did
succeed and the other two values are real.

`cargo test --workspace` gives 855 passed with 2 ignored, the same 855 on
faer, 149 linear algebra, 506 pytest, 332 node, and fmt, clippy,
wasm-check and ruff clean. Every one was run by the orchestrator.

**One number for the performance review, not acted on.** One fit over 249
distances takes 45.8 µs in the test profile on this machine, measured on
24 September 2026 with a throwaway test that was not committed.

**2.3, the curve on the reference dataset, against R.** Commit
`cd43a42`. `calc_ld_and_dist` now fits the curve of each population when
its pass ends, over that population's own individuals, and
`LdBins::decay()` gives it, which is the method of "The Rust interface"
that work package 1 left out on purpose.
`tests/reference/ld/ld.decay.txt` is stored, and `run_plink2.sh` runs
`decay.py` and then `decay.R` and compares a tenth file, naming R among
what it needs.

**popnei and R agree far inside what the spec asks.** `decay.R` printed
the spec's table digit for digit, and popnei sits 3.5e-8 relative from it
for the population of every individual, 5.3e-8 for `pop_a` and 5.5e-9 for
`pop_b`, on both the ρ per base pair and the half distance, where 1e-6 is
asked. The r² at distance 0 is bit-identical, 0 apart, where 1e-12 is
asked. So the half distance the waiting web application needs comes out
of popnei at 6810.571 base pairs for the first population, against R's
6810.5712522189806.

The stored file carries more than the table: what the two optimisers of R
disagree by, how badly the curve describes this dataset bin by bin, what
Sved's curve would have given, what fitting the bins instead of the pairs
would have cost, and what the other answers to the sample size would have
given. All of it is what the spec already states, written where a later
session can check it again.

Both thread pools genuinely ran: `rayon::current_num_threads()` is
asserted inside the pool, and building the pools with seven threads too
many makes that assertion fail. The two tests take 0.04 s together.

`cargo test -p popnei --lib ld::decay -- --list` prints `15 tests`, where
it printed 13. `cargo test --workspace` gives 857 passed with 2 ignored,
the same 857 on faer, 149 linear algebra, 506 pytest, 332 node, fmt,
clippy, wasm-check and ruff clean, and `run_plink2.sh` into an empty
directory exits 0 naming no differing file, with the stored
`ld.decay.txt` identical to the one R wrote. Every one was run by the
orchestrator.

**2.4, the fitted curve in Python.** Commit `4bdbac2`. `LdDecay`, a
frozen dataclass with `rho_per_bp`, `r2_at_zero` and `half_dist`, and
`decay_per_pop` on `LdAndDistPerPop`, in `python/popnei/ld.py`. The pyo3
crate reads the three numbers straight off the core's result and works
out nothing of its own.

The orchestrator read the three numbers through the whole chain rather
than trusting the test:

| | popnei | the spec, from R | the gap |
|---|---|---|---|
| ρ per base pair | 0.00031727348310939046 | 0.00031727347196446889 | 3.51e-8 |
| r² at distance 0 | 0.46198347107438015 | 0.46198347107438015 | 0 |
| the half distance, bp | 6810.571012984072 | 6810.5712522189806 | 3.51e-8 |

The spec asks 1e-6 on the first and the third and 1e-12 on the second.

The task asserted the other two rows as well as the first, which the
deliverable did not ask for, because with one population nothing at the
Python layer would catch a binding that handed every population the
first one's curve. It also built both of the spec's ways of having no
curve: three variants at 10, 20 and 30 bp with a `max_dist` of 15, which
leaves pairs at one distance; and a population of one individual, which
keeps two variants but has one dosage at each and so counts no pair at
all.

`uv run pytest tests/test_ld.py -k ld_and_dist` gives `9 passed, 16
deselected`, where it gave 7 and where at the start of the plan it exited
5 with nothing matched. The whole pytest suite is 508 where it was 506,
the core crate stays at 857 with 2 ignored and the same on faer, 149
linear algebra, 332 node, and fmt, clippy, wasm-check and ruff clean.
Every one was run by the orchestrator.

**2.5, the fitted curve in TypeScript.** Commit `d2aaeaa`. `decayPerPop`
on the result, an object of population name to `{rhoPerBp, r2AtZero,
halfDist}`, three numbers and not arrays as the spec asks, with `LdDecay`
exported from both `node.ts` and `web.ts`. The three values cross as one
array of a value per population and are read at that population's index,
so nothing is cut.

The WebAssembly build gives the same three numbers as Python, to the bit.
That is worth recording because the review of work package 1 found the
WebAssembly build 2 units in the last place from the native one on a
bin's mean: the fit uses the four operations for the curve and the sum,
and `powf` only to move along the exponent of 10 in the grid and the
search, so the platform can move the answer only through where the search
looked, and here it did not.

The task added a node test of a population with no curve, which the
deliverable did not ask for, because it is the only test in any suite
that shows a NaN of the core reaching a JavaScript user as `NaN`. It
needed no new fixture file.

`node --test js/popnei/test/ld.test.ts` runs 19 tests where it ran 17,
and the whole node suite is 334 where it was 332.

### The seven deliverables of work package 2, each checked by the orchestrator

| the deliverable | the command | what it gave |
|---|---|---|
| 1, `ld.decay.txt` written again and compared | `tests/reference/ld/run_plink2.sh` into an empty directory | exit 0, no file named as differing, and `cmp` finds the stored `ld.decay.txt` and `ld.bins.txt` both identical to the written ones |
| 2, `fit_ld_decay` against the curve-derived table | `cargo test -p popnei --lib ld::decay -- --list` | `15 tests`, where it printed `0 tests` before this work package |
| 3, the cases with no curve | the same list | both ends of the searched range, fewer than two distances and no pair at all, each with its own test, and the four errors beside them |
| 4, the three rows against R | `cargo test --workspace` | `857 passed; 0 failed; 2 ignored`, and `149 passed` in the linear algebra crate |
| 5, the memory asked and refused | the same run | the test of task 2.1, which fails when the checked allocation is made ordinary |
| 6, the Python layer | `uv run pytest tests/test_ld.py -k ld_and_dist` | `9 passed, 16 deselected`, where the plan started at exit 5 with nothing matched |
| 7, the TypeScript layer | `npm run build && npm test` in `js/popnei` | `tests 334, pass 334, fail 0`, and 19 in the ld file where there were 10 at the start of the plan |

`cargo test -p popnei --no-default-features` gives the same 857 on the
faer backend, and `cargo fmt --all --check`, `cargo clippy --workspace
--all-targets -- -D warnings`, `cargo wasm-check`, `uv run ruff format
--check` and `uv run ruff check` are clean.

### What the review of work package 2 found

Five reviewers read the work package at `d2aaeaa`, one for each category
that applied. Twenty-one findings came back. One of them was a wrong
number that no deliverable of the plan would have caught.

**The fit lost about 1.5 digits to cancellation, and the spec is what
told it to.** "The curve that is fitted" wrote the quantity to make
smallest as Σ over the distances of [n·f² − 2·S·f]. That is algebraically
right and numerically poor: it is the quantity that matters,
Σ n·(mean − f)², minus a constant. On the first reference population the
constant is 475.36, so the number being minimised sits at −452.15 while
the part that depends on the fitted value is 23.21 — nineteen times the
magnitude of the signal the search has to resolve, and one step of an
`f64` there is 5.68e-14 against 3.55e-15 of the signal. The spec's own
derivation already gives the better form in words, two paragraphs above
the formula it then wrote out.

It was reachable and not theoretical. Over 2052 values of the scaled
recombination per base pair that a slowly falling population would give
at the default `max_dist` of 1000000, **18 of them came out past the
1e-6 the spec itself promises**, the worst at 2.2e-6. The test written
before the fix sits at one of them: 1.537e-6 before, 3.474e-10 after.
Against R the gap fell from 3.51e-8 to 1.03e-9 for the population of
every individual, from 5.27e-8 to 9.28e-9 for `pop_a` and from 5.46e-9 to
2.21e-9 for `pop_b`. The spec's formula was corrected first, in a commit
of its own, and one fit still takes 44.10 µs against 44.24 µs before, so
the accuracy came free.

**Five tests could not have failed**, each shown by breaking the code and
watching the suite stay green:

- The test the plan names as the guard of the whole fit was blind to
  three of the four things it guards. Its scaled recombination of 0.0001
  is grid point 80 exactly, so the grid alone answers it and the search
  never has to improve on the answer, and its table holds one pair at
  every distance, so the weight cancels out of the sum. Replacing the
  golden ratio with 0.5, replacing it with its complement, loosening the
  stopping rule from 1e-9 to 1e-2, and dropping the pair count from the
  sum each left it passing. The plan's sentence "a fit that is wrong in
  the grid, in the bracket or in the stopping rule fails it" was true only
  of the grid.
- The node suite would not have caught a binding that handed every
  population the first one's curve: with that break made, all 334 node
  tests passed, where the same break in Python fails pytest at once.
- Changing "fewer than two distances" to "fewer than three" left all 119
  tests of the module passing while turning a real curve into three NaN.
- Shrinking the bisection's upper bracket from 10⁶ to 10 left all 119
  passing while making the half distance NaN at three and four
  individuals.
- One assertion about the grid compared a constant with itself, so the
  grid's lower end and its spacing were pinned by nothing: moving the
  bottom from 10⁻¹² to 10⁻¹¹, and separately the spacing from a tenth of
  a decade to a fifth, each left all 119 passing.

Each of those now fails under the break that exposed it. Closing the
first needed a table the reviewer's own suggestion would not have closed:
a table built from one curve has a sum of 0 at the true answer under any
weighting, so the fixture mixes two curves, and R 4.6.1 was asked what
the weighted fit of it should be.

**The rest were smaller.** A defect of popnei would have been reported as
a wrong argument, because a population with no dosages and a population of
no individuals were folded into one number; it is unreachable today and is
now a defect-class error of its own. Three getters were missing
`#[must_use]`. The type-level documentation of the fitted curve did not
say that the half distance can be NaN on its own. The four new error
messages were asserted nowhere, only their variants matched. A comment in
the binning was false for a distance below `min_dist`, because the
subtraction used is symmetric. Two counts used a saturating sum where the
same file argues, ten lines away, that a count must not stop in silence.
The glossary had no word for the curve's value at distance 0. And the
TypeScript package wrote the same defect check twice.

**Two findings were not fixed, and both are above this plan.** A pass
cannot be interrupted, which is true of every calculation in popnei and
changes a signature to mend. And the four new error arms of the Python
binding cannot be unit tested, because that crate is built as a shared
library with `test = false` and this machine has no linkable libpython.
Both are under "What is asked of the owner".

**One finding did not hold**: an `#[expect]` whose written bound a
reviewer said was wrong had already been corrected by the cancellation
fix, which the fixer checked rather than assumed.

**Two numbers of the spec did not reproduce and were corrected.** The
curve at one individual is 1.1983471074380165 and not
1.1983471074380166, one step of an `f64` below. And narrowing the
searched range moves the first population's fit by 9.3e-9 of itself and
not by the 8.9e-9 the spec claimed, measured in the Rust and again in an
independent copy of the fit in Python; the paragraph below it already
said 9.3e-9, so the section disagreed with itself.

## For the performance review that follows this plan

Every number here was measured on this machine on 24 September 2026 and
none of it was acted on, because the owner decided that day that the
speed of this pass is a review of its own.

- One fit over 249 or 250 distances takes 44.10 µs in the test profile.
  The spec's count of what it costs is right: an instrumented run
  evaluated the sum 183 times, 141 for the grid, 2 to open the search and
  40 for its steps.
- Each block builds two sets of dosages for each population, one over the
  block to work out the frequencies and one over the whole window, so a
  pass at blocks of 7 variants with a window of 250 reads the genotypes of
  a variant about 36 times over.
- `LdDosages::rows` copies the three matrices once per tile pair per step.
- A block that has fallen out of the window is held for one step longer
  than it is needed, which is what makes the counts right at every block
  size.
- The compaction at the end of the pass walks one entry for every distance
  from `min_dist` to `max_dist`, so the time grows with `max_dist` and not
  with the data: 0.06 s at 10⁶, 17.13 s at 10⁹, 190.85 s at 10¹⁰.
- `cargo test -p popnei --lib ld::decay` runs in 0.04 s, so nothing this
  plan added costs the suite anything.
- `LdDosages::of_block` walks the rows of a block serially, where section
  3 of `docs/architecture.md` puts rayon on per-variant work.

## How the work went

**This section is for whoever next revises a skill or writes a plan. The
owner can stop here.**

**A deliverable that pins a property is worth more than one that pins a
number at a point.** Work package 1 promised that the bins come out the
same at every block size, and that deliverable found two rules of the spec
that broke the promise, one of them a wrong count of pairs, 272 where 328
was right. Work package 2 promised agreement with R to 1e-6, and the fit
met it on the reference dataset while missing it on 18 of 2052 nearby
datasets. The first deliverable caught its defect because the property was
checked across a range of conditions; the second nearly missed its defect
because the number was checked at one point. A plan that asks for a number
should say across what range it must hold.

**Three of the four things the plan believed about its own risks were
wrong.** It said the fixed-multiple tiling was what bought bit-identity:
removing the alignment changes no number, and what buys it is the
accumulation order, which the plan did not mention. It said the
curve-derived test would fail a fit wrong in the grid, the bracket or the
stopping rule: it fails only the grid. It warned that a smallest sum at an
end of the searched range was the case easiest to miss: that one was
right, and the task built it without prompting. A plan's "What could go
wrong" is a hypothesis, and a review that only checks the code against the
plan inherits the plan's mistakes. Both reviews were told to review the
spec's own changed sentences, and several of the best findings came from
there.

**The reviewer that broke code found the defects; the reviewers that read
code found documentation — with one exception that matters.** Across both
work packages the `tests` reviewer, which mutates and reports what
survives, found what a user would have suffered. So did `numbers`, by
reading alone: the cancellation in the fit was found by looking at an
expression and recognising it, not by running anything. So the lesson is
not to drop the readers but to tell every reviewer to break something
before it reports.

**What a task of this size really costs, in the tokens its subagent
used.** The six tasks of work package 1 cost 1.20 million between them;
its review cost 1.15 million and its three rounds of fixes 0.59 million,
so the checking cost 1.45 times the building. The five tasks of work
package 2 cost 0.90 million; its review cost 0.87 million and its three
rounds of fixes 0.63 million, 1.68 times the building. A plan that budgets
for the tasks alone budgets for two fifths of the work.

**Two subagents died and cost almost nothing, because they committed as
they went.** The first subagent of task 1.5 died on an API refusal before
doing any work and left the tree clean; the retry succeeded. A later fix
round died on a session limit after committing all five of its fixes but
before taking its final measurement, which the orchestrator took itself. A
round that commits per finding survives its own subagent.

**The shared tree cost one mistake, and it was the orchestrator's.** It
ran a mutation test in the worktree while five reviewers were reading it,
and three of them reported the modified file as something to look into.
Nothing was lost, because the file had been copied first, but an
orchestrator that wants to break code should do it somewhere else.
