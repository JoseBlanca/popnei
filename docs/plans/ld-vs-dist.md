# Plan: the curve of r² against distance, per population

24 September 2026. State: approved by the owner on 24 September 2026, under
way since 24 September 2026, and the work is recorded in
`docs/reports/ld-vs-dist.md`. It builds the item
"LD against distance, per population" of `docs/specs/ld.md`, which gives,
for each population of a dataset, how the r² of a pair of variants falls
off as the two move apart along a chromosome: the fall-off in bins of
distance, a curve fitted to the pairs, and the distance at which r² has
fallen to half. The item was written and reviewed on 22 September 2026,
and the curve and the half distance were added to it on 24 September 2026,
first read, reviewed and committed before this plan. It is carried out in
the worktree `.claude/worktrees/ld-vs-dist` on the branch
`plan/ld-vs-dist`, and the report is `docs/reports/ld-vs-dist.md`.

## In and out

Built: `calc_ld_and_dist` and `fit_ld_decay` in the core crate,
`calc_ld_and_dist_per_pop` in Python and `calcLdAndDistPerPop` in
TypeScript, each with the bins, the fitted curve and the half distance,
and each ending at a test that compares with plink2, with R or with
pyNei.

Not built, with where it goes:

- **The speed of the pass and of the fit.** The spec states what the fit
  costs in evaluations, 183 passes over the distances that hold a pair,
  and nobody has timed one: at the default `max_dist` of 1000000 that is
  183 passes over a million cells. The owner decided on 24 September 2026
  that it goes to a performance review of its own after this plan, as the
  r² matrix did on 23 September 2026, and the "Speed" section of the spec
  says nothing about this item until that review has run. So no
  deliverable here is a time, and a task that finds the fit slow reports
  it and does not redesign it.
- **Open 1 of `docs/specs/ld.md`**, whether the result also carries a
  sample of individual pairs for a scatter. Its "meanwhile" is the bins
  alone, which is what every task below builds; the other answer adds a
  seed, a count and a sampling rule to the same pass, which would be a
  task of work package 1 and a field in each of the three layers.
- **Opens 2 and 3 of that spec** change nothing here. Open 2 is the major
  allele of a variant with half called genotypes, which
  `variant::the_major_allele` already decides for the whole crate, and
  Open 3 is the speed target of the r² matrix.
- **A spread around the half distance**, and the effective size and the
  recombination rate apart from their product: "Not in this spec" of
  `docs/specs/ld.md` says why neither is built.

## What has to be in place

All three layers exist, so every check of the `coding` skill runs during
this plan and none is reported as not there: `cargo fmt --all --check`,
`cargo clippy --workspace --all-targets -- -D warnings`, `cargo test
--workspace`, `cargo test -p popnei --no-default-features`, `cargo
wasm-check`, `uv run ruff format --check && uv run ruff check`, `uv run
maturin develop && uv run pytest`, and `npm run build && npm test` in
`js/popnei`.

On the commit this starts from, `main` at `cce5801` with the two spec
commits of this branch on top:

- `cargo test --workspace` passes, 787 tests in the core crate and 149 in
  the linalg crate.
- `cargo test -p popnei --lib ld:: -- --list` prints `49 tests`, and both
  `ld::dist` and `ld::decay` print `0 tests`.
- `uv run pytest tests/test_ld.py -k ld_and_dist` exits 5, which is what
  pytest gives when `-k` matched nothing; `tests/test_ld.py` has 16 tests
  and `js/popnei/test/ld.test.ts` has 11, all of the r² matrix.
- `crates/popnei/src/ld.rs` has `LdDosages` with `of_block`, `rows`,
  `has_variance`, `maf` and `dosages`, `r2_between`, and `calc_r2_matrix`
  with its tiles. Nothing of the fall-off against distance is there:
  `grep -rn "ld_and_dist\|LdAndDist\|LdBins\|LdDecay" crates python js
  tests` finds nothing.
- plink2 v2.0.0-a.7.7 and R 4.6.1 are on the machine.
  `docs/reports/ld-method/bins.py`, `decay.py` and `decay.R` were run on
  24 September 2026 and printed the three tables of bins and the three
  rows of the curve that the spec has;
  `docs/reports/ld-method/decay_truth.py` ran in 31 s and printed the
  numbers behind the choice of the curve and of n.
  `tests/reference/ld/run_plink2.sh` was run into an empty directory on
  24 September 2026: it writes its eight files, finds each the same as
  the copy stored beside it, prints nothing of its own and exits 0.

## Work package 1: the bins of r² against distance

### What it gives

A user calls `calc_ld_and_dist_per_pop` in Python, or
`calcLdAndDistPerPop` in TypeScript, over a `Variants` and a dict of
populations, and gets for each population one row per bin of distance
with how many pairs it holds, their mean r² and its standard deviation,
and how many variants that population kept at its major allele frequency.

### Its deliverables

1. `tests/reference/ld/ld.bins.txt` holds what
   `docs/reports/ld-method/bins.py` prints for the three populations, and
   `tests/reference/ld/run_plink2.sh` writes it again and compares it
   with the stored copy, as it does for its eight other files. The check:
   run that script into an empty directory, which prints nothing and
   exits 0. The file does not exist today and the script does not write
   it.
2. The worked example of "How it is verified" of the item, the five
   variants of `tests/reference/ld/example.vcf` with `min_dist` 1,
   `max_dist` 4000 and `num_bins` 2, is a cargo test at
   `calc_ld_and_dist` that asserts the two counts of pairs exactly and
   the two means within 1e-12 relative, with the pairs that hold the
   variant of no variance in no bin.
3. The three tables of bins of that same part are cargo tests at
   `calc_ld_and_dist` on `ld.vcf.gz` read with the VCF reader, the counts
   of pairs exact and the means and the standard deviations within 1e-12
   relative, run with blocks of 7, 64 and 500 variants and with rayon at
   1 thread and 4, which have to give the same numbers to the bit.
4. Every error of "The cases" is a cargo test that asserts the argument
   and the value are in the message: a `min_dist` above `max_dist`, a
   `num_bins` of 0, a `max_allowed_maf` that is not a number from 0 to 1,
   a negative `min_dist`, a population with no individual, an individual
   the dataset has not, an individual asked for twice, and a pass with no
   variant. The empty cases of the same part, every bin empty and a
   `num_vars_per_pop` of 0, are cargo tests beside them and are not
   errors.
   The check for deliverables 2 to 4: `cargo test -p popnei --lib
   ld::dist -- --list` prints them and a count above 0, where it prints
   `0 tests` today.
5. `calc_ld_and_dist_per_pop` and `LdAndDistPerPop` are in
   `python/popnei/ld.py` with the signature of "Its Python function", and
   a pytest test asserts the ten rows of the first table and that pyNei,
   asked through `filter_samples` and `filter_by_maf` for the variants
   each population keeps at a major allele frequency of 0.8, gives the
   same 396 and 402 that `num_vars_per_pop` gives. The check: `uv run
   pytest tests/test_ld.py -k ld_and_dist` runs them and passes, where it
   exits 5 today.
6. `calcLdAndDistPerPop` is in `js/popnei/src/ld.ts` with the five typed
   arrays per population, and a node test on `ld.vcf.gz` read from a
   `Uint8Array` asserts the ten counts of pairs and the ten means of the
   first table. The check: `npm run build && npm test` in `js/popnei`,
   which by then runs more ld tests than the 11 it runs today, all of
   which are of the r² matrix.

### What it stands on

`main`, and nothing of work package 2. `LdDosages` with its `maf` and its
`rows`, `r2_between` and `variant::the_major_allele` are on `main` from
the plan of 22 September 2026, and this work package uses them as they
are.

### Its tasks

- [x] 1.1 The window over the blocks, in `crates/popnei/src/ld.rs`, from
  "How it runs" of the item: hold the blocks whose variants are within
  `max_dist` of the newest variant read and on its chromosome, and drop a
  block when every variant of it is further back or on another
  chromosome. It serves deliverable 3, whose three block sizes are what
  tell a window that keeps too much from one that drops too soon. Its
  cargo tests are of the window alone, on readers built in the test, and
  they are the first tests `cargo test -p popnei --lib ld::dist --
  --list` prints, where it prints `0 tests` today.
- [x] 1.2 The dosages of each population over the window: for each
  population, the variants of the held blocks that pass its
  `max_allowed_maf`, and the matrices `LdDosages` holds over the
  individuals of that population, as "How it runs" has it. The major
  allele frequency is `LdDosages::maf`, which is already over the
  individuals it was built with, and the rule is the one
  `docs/specs/filters.md` defines for `filter_by_maf`. It serves
  deliverables 3 and 5: the 396 and 402 of deliverable 5 are the count
  this task produces. Needs 1.1.
- [x] 1.3 The pairs of the window in tiles and the bins, from "How it
  runs" for the tiles and "What it gives" for the arithmetic of the bins,
  and from "The cases" and "Its Python function" for what each argument
  refuses: the bin a distance falls in, the count and the two sums per
  bin, `LdAndDist`, `LdBins` and `calc_ld_and_dist`. How the tiles are
  cut is what makes the result the same to the bit at any block size, and
  "How it runs" is the only place that says how; read that paragraph
  before writing the loop. It serves deliverables 2 and 4, and the worked
  example is its first test. Needs 1.2.
- [x] 1.4 The three tables of "How it is verified" as cargo tests, and
  `ld.bins.txt` stored and compared by `run_plink2.sh`. It serves
  deliverables 1 and 3. The numbers are in the spec and are not
  recomputed here; what this task runs is `bins.py`, to store what it
  prints. Needs 1.3, and it can run beside 1.5 and 1.6: the three touch
  different files.
- [x] 1.5 The Python layer, from "Its Python function": the pyo3
  function in `crates/popnei-python/`, which keeps the reader of the
  chain as `calc_kosman_sums` does so that the result carries the pass
  stats, and `calc_ld_and_dist_per_pop` with `LdAndDistPerPop` in
  `python/popnei/ld.py`, with the pytest tests of deliverable 5. It
  serves deliverable 5. Needs 1.3.
- [x] 1.6 The TypeScript layer, from the last paragraph of "Its Python
  function": the wasm-bindgen function in `crates/popnei-js/` and
  `calcLdAndDistPerPop` in `js/popnei/src/ld.ts` with its result object,
  and the node test of deliverable 6. It serves deliverable 6. Needs 1.3.

### What could go wrong

The result has to be the same to the bit at every block size and thread
count, and what buys that is cutting the tiles at fixed multiples counted
from the first variant of the pass, which is not how `calc_r2_matrix`
tiles: that one tiles a matrix it has already materialised, where this
one tiles a window that moves. A tiling that follows the blocks instead
passes deliverable 2 and fails deliverable 3 at the second block size.

Each population keeps different variants, so the window holds one set of
rows per population and the tiles of one population are not the tiles of
another. The memory of "How it runs", 6 MB for 250 variants of two
populations of 500 individuals, is per population.

## Work package 2: the fitted curve and the half distance

### What it gives

The same three calls also give, for each population, the curve fitted to
its pairs and the distance at which r² has fallen to half: `rho_per_bp`,
`r2_at_zero` and `half_dist` in Python and in the core, `rhoPerBp`,
`r2AtZero` and `halfDist` in TypeScript.

### Its deliverables

1. `tests/reference/ld/ld.decay.txt` holds what
   `docs/reports/ld-method/decay.R` prints, and `run_plink2.sh` writes it
   again and compares it, running `decay.py` before it. The script then
   needs R as well as plink2 and uv, and says so where it says what it
   needs. The check is the same run into an empty directory, which prints
   nothing and exits 0.
2. `fit_ld_decay` is in `crates/popnei/src/ld.rs` with the signature of
   "The Rust interface", and the cargo test that needs no plink2 and no R
   asserts what "How it is verified" gives it: r² taken from the curve
   itself at a ρ per base pair of 0.0001 and n of 100, one pair at each
   of the distances 1000 to 250000, giving back 0.0001 and a half
   distance of 21608.135872529165, both within 1e-6 relative.
3. The cases of "The cases" that have no curve are cargo tests at
   `fit_ld_decay`: pairs at fewer than two distances, no pair at all, and
   a smallest sum that falls at either end of the searched range of the ρ
   per base pair, each giving NaN for all three values and not an error.
   The errors of that function are cargo tests beside them.
4. The three rows of the table of "How it is verified" are a cargo test
   at `calc_ld_and_dist` on `ld.vcf.gz`, compared within 1e-6 relative
   for the ρ per base pair and the half distance and within 1e-12
   relative for the r² at distance 0, at the same three block sizes and
   two thread counts as deliverable 3 of work package 1.
   The check for deliverables 2 to 4: `cargo test -p popnei --lib
   ld::decay -- --list` prints them and a count above 0, where it prints
   `0 tests` today.
5. The memory of the counts and the sums per distance is asked with
   `try_reserve_exact` and refused rather than taken, with a cargo test
   like the one `calc_r2_matrix` has for its matrix, which today is
   `a_matrix_this_machine_has_not_the_memory_for_is_an_error_and_not_the_end_of_the_process`.
6. `decay_per_pop` and `LdDecay` are in `python/popnei/ld.py`, and a
   pytest test at `calc_ld_and_dist_per_pop` asserts the three numbers of
   the first row of that table and that a population left with pairs at
   one distance has NaN in all three. The check is the same `-k
   ld_and_dist` run.
7. `decayPerPop` is in the TypeScript result as three numbers and not
   arrays, and the node test asserts the three numbers of the first row.
   The check is `npm run build && npm test` in `js/popnei`, which by
   then runs more ld tests than the 11 it runs today.

### What it stands on

Work package 1, finished and reviewed, and not task 1.3 alone: the
accumulation of task 2.1 is written into the tile step that 1.3 builds
and the literals of deliverable 4 are read out of the pass 1.3 makes, but
1.4 to 1.6 are where the bins are first compared with plink2 and with
pyNei, and a curve fitted on top of bins nobody has checked would leave
two things to look at when one number came out wrong. The other way round
is not possible: nothing of work package 1 needs anything of this one.

### Its tasks

- [x] 2.1 The pairs counted at each distance: for each population, a
  count and a sum of r² for every distance from `min_dist` to
  `max_dist`, filled in the same tile step that fills the bins, asked of
  the machine with `try_reserve_exact` before the pass, and compacted
  when the pass ends to the distances that hold a pair. Its part of the
  spec is the last two paragraphs of "How it runs". It serves deliverable
  5 and feeds 2 and 4.
- [x] 2.2 `fit_ld_decay` and `LdDecay`: the curve, the sum it makes
  smallest, the grid of 141 values, the golden section, the bisection for
  the half distance, and the cases with no curve. Everything it needs is
  in "The curve that is fitted" of the item. It serves deliverables 2 and
  3, and its first test is the one whose r² comes from the curve. It
  touches no reader, so it can run beside 2.1.
- [x] 2.3 The three rows of the table of "How it is verified" as a
  cargo test, `calc_ld_and_dist` calling
  `fit_ld_decay` for each population when its pass ends, and
  `ld.decay.txt` stored and compared by `run_plink2.sh`. It serves
  deliverables 1 and 4. The numbers are in the spec; what this task runs
  is `decay.py` and `decay.R`, to store what they print. Needs 2.1 and
  2.2.
- [x] 2.4 The Python layer, from "Its Python function":
  `decay_per_pop` on `LdAndDistPerPop` and the
  `LdDecay` dataclass, with the NaN of a population with no curve, and
  the pytest tests of deliverable 6. It serves deliverable 6. Needs 2.3.
- [x] 2.5 The TypeScript layer, from the last paragraph of "Its Python
  function": `decayPerPop` on the result object and
  the node test of deliverable 7. It serves deliverable 7. Needs 2.3, and
  it can run beside 2.4.

### What could go wrong

A fit is a wrong number and not a crash, so task 2.2 is the one whose
failure would be silent. What guards it is deliverable 2, whose r² is
taken from the curve itself, so the answer is known before the fit runs;
a fit that is wrong in the grid, in the bracket or in the stopping rule
fails it.

The three values are NaN in two different situations and the second is
the easy one to miss: the smallest sum of the search falling at an end
of the range of the ρ per base pair, which says these pairs pin no
fall-off down. A search that returns that end value gives a half distance
of about 10⁻¹² bp or of about 10¹² bp, which is a number and not a NaN,
and no other deliverable would catch it.

The fit is the only part of popnei that reads a value at a distance the
dataset may not hold, and `max_dist` sets how much memory that costs: 16
MB for each population at the default, and the call is refused when the
machine has not got it rather than taking it. A run of twenty
populations asks for 320 MB before the first block is read.

## How the whole plan is checked

The sum of the two work packages, and nothing besides: the checks of the
`coding` skill all pass, `run_plink2.sh` writes its ten files and finds
them the same as the ones stored beside it, and `cargo test -p popnei
--lib ld:: -- --list` prints the 49 tests that were there and the ones
this plan added.
