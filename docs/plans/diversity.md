# Plan: how much variety each population holds

24 September 2026. Approved by the owner on 24 September 2026. State:
under way since 24 September 2026, on the branch `plan/diversity` in the
worktree `.claude/worktrees/diversity`, which branches from
`spec/diversity` because neither the spec nor this plan is on `main`;
its report is `docs/reports/diversity.md`. It builds
the `diversity` module of
`docs/specs/diversity.md`: one pass over the variants that gives, for each
population, the alleles it called and the private ones among them, the
variants that vary in it, those three also taken down to a common number of
called alleles, the folded site frequency spectrum projected to that number,
and F_IS. It runs through the core crate, both binding crates, the Python
package and the TypeScript package, and it ends with the measurement that
sets the speed target the spec leaves open.

The spec has no open point. The four it had the owner decided on 24
September 2026 and each is written where it applies, so nothing here waits
on an answer.

## In and out

Built: the module, its errors, the two binding crates, the Python function
and its result, the TypeScript function and its result, the reference data
under `tests/reference/diversity/` with the scripts that make it, and the
measurement.

Not built, with where each goes:

- A sequence of draw sizes, which would give the curve of allelic richness
  in one pass. The owner decided on one number; a user who wants the curve
  calls the function once per point. "Its Python function" of the spec.
- ADZE. The standardized private alleles are checked by enumerating every
  draw instead, which the spec's "How it is verified" of that item gives
  with the 18 pairs it agrees on.
- The unfolded spectrum and the joint spectrum of two populations, and the
  rarefaction of the heterozygosities. "Not in this spec".
- Any change to `docs/specs/stats.md` beyond the line that already points
  here, which the branch this plan is carried out from already holds.

## What has to be in place

- The branch `plan/diversity`, from `main` once `spec/diversity` is on it,
  in a worktree of its own. Run on `spec/diversity` at fa1fd52 on 24
  September 2026: `cargo fmt --all --check` exit 0; `cargo test
  --workspace` `787 passed` and `149 passed`, 2 ignored; `uv run ruff
  format --check` `33 files already formatted` and `uv run ruff check`
  `All checks passed!`; `uv run maturin develop && uv run pytest` `499
  passed`; `cargo clippy --workspace --all-targets -- -D warnings` no
  warning; `cargo test -p popnei --no-default-features` `787 passed`, 2
  ignored. These are the counts that "more than" is counted from below.
- `npm run build` in `js/popnei` before any `npm test` there, and the
  orchestrator runs it once at the start of the plan and not only in task
  2.6. A fresh worktree has neither `node_modules` nor `js/popnei/wasm/`,
  which `.gitignore` leaves out because `wasm-bindgen` writes it, and
  `npm test` alone fails there twice over: first with `tsc: command not
  found`, and then, once `npm install` has run, with `src/web.ts(79,3):
  error TS2322` on `return loading`, because with no `wasm/popnei.d.ts`
  the loader is an `any` and `??=` cannot narrow it. Both were seen in
  this worktree on 24 September 2026, and `npm run build`, which runs
  `npm install`, the release build of `popnei-js` for
  `wasm32-unknown-unknown` and `wasm-bindgen` before `tsc`, is what makes
  them go away. The error names a file nobody touched, so an orchestrator
  meeting it in task 2.6 would look for it in the wrong place. After
  `npm run build`, `npm test` gave `tests 325`, `pass 325`, `fail 0` in
  this worktree on 24 September 2026, which is the count that "more than"
  is counted from below for deliverable 7.
- R 4.6.1 at `/opt/homebrew/bin/Rscript`, holding `vegan` 2.7.6,
  `adegenet` 2.1.11 and `poppr` 2.9.8. Checked on 24 September 2026 by
  running each of the three on the panel; their numbers are the ones in
  the spec's tables. `hierfstat` is not there and does not install.
- `uv` able to make a Python 3.12 environment and install `dadi` 2.4.4 and
  `scikit-allel` 1.3.13 into it. Both were installed and run on 24
  September 2026. `dadi` does not build on the project's Python, 3.14 with
  the free threading build.
- The panel: `tests/reference/stats/panel.vcf.gz` and
  `tests/reference/stats/panel_pops_bcftools.txt`, which
  `tests/reference/stats/make_reference.py` writes and which are
  committed. Nothing here rebuilds them.
- The file of the measurement, outside the repository:
  `/Users/jose/devel/popnei-bench/big.vars`, 100000 variants of 1000
  individuals, which `docs/plans/stats.md` also measures on. Task 4.1
  writes it if it is not there.

Every layer exists, so the seven commands of "Before the work is called
done" of the `coding` skill run for every task, and `npm run build && npm
test` for a task that touches `crates/popnei-js` or `js/popnei`.

The checks below fail today. `cargo test -p popnei --lib -- diversity::
--list` prints `0 tests`. `tests/test_diversity.py`,
`js/popnei/test/diversity.test.ts` and the directory
`tests/reference/diversity/` do not exist.

Tasks run one after another unless the plan says that two can run side by
side: every task builds the same workspace in the same worktree.

A task whose failure would be silent, a wrong number and not a crash or a
failing test, has a commit of its own and names the deliverables that
guard it.

## Work package 1: the reference numbers and the programs that give them

**What it gives.** Work packages 2 and 3 get every number of the spec's
"How it is verified" as a file in the repository, with the script that
made it, so that each of their checks compares against something stored
and not against a number typed out of a session. It stops before Python:
no code of the module exists yet, and the comparisons that read these
files are deliverables of the two work packages after it.

One number of "How it is verified" does not become a file here, and task
1.1 found it: the standardized private alleles of the panel at a draw of
20, 0.0112196177, 0.0099715392 and 0.0089014974. No outside program gives
any standardized value, so no script of this work package can write them.
They stay literals of a pytest test, which is what deliverable 3 of work
package 3 already asks for, and the spec was corrected in commit `b580a21`
to name where they do come from, `docs/reports/diversity-method/panel.py`.

It goes first because every later check reads these numbers, and because
it is the one part that could change the plan: if `dadi` cannot be
installed the same way twice, the folded spectrum loses its only outside
check and the spec's "How it is verified" of that item has to be written
again.

**Deliverables.**

1. `tests/reference/diversity/make_reference.R` runs `vegan`, `adegenet`
   and `poppr` on the panel and writes what each gives beside itself.
   Check: `Rscript tests/reference/diversity/make_reference.R` exits 0,
   and `git status --short tests/reference/diversity` shows nothing after
   a second run. The stored files hold the alleles called, 2373, 2377 and
   2384, the private alleles, 0, 0 and 1, and the standardized allele
   counts at a draw of 20, 1.9283948650, 1.9219209943 and 1.9197370844,
   which are the numbers of "How it is verified" of "The number of
   alleles" and of "The private alleles".
2. `tests/reference/diversity/make_reference.py` runs `scikit-allel` and
   `dadi`, each in a Python 3.12 environment it makes with `uv venv`, and
   writes the folded spectrum of the three populations at a draw of 20,
   33 values, and the three plain form F_IS values. Check: `uv run python
   tests/reference/diversity/make_reference.py` exits 0, the script
   compares what it got with the tables of the spec and fails if they
   differ, and `git status --short tests/reference/diversity` shows
   nothing after a second run.
3. `tests/reference/diversity/enumerate_private.py` writes the
   standardized private alleles of the 18 pairs of a case and a
   population of "How it is verified" of "The private alleles", each
   computed twice, once by the closed form and once by enumerating every
   draw, in exact rational arithmetic. Check: it exits 0, its output
   holds a difference of 0 for all 18, and the file it writes is what the
   cargo tests of task 3.3 read their literals from.
4. A page beside them saying what each script is, which program wrote
   which file, and that `hierfstat` does not install. Check: it names all
   five programs with their versions and both environments.

**What it stands on.** Nothing in this plan. Outside it, the panel and the
five programs of "What has to be in place".

**Tasks.** 1.1, 1.2 and 1.3 touch different files and can run side by
side; 1.4 comes after all three.

- [x] 1.1 `make_reference.R` and its output, from "How it is verified" of
  "The number of alleles", of "The private alleles" and of "The variable
  variants" of `docs/specs/diversity.md`. Serves deliverable 1.
- [x] 1.2 `make_reference.py`, the two environments it makes and its
  output, from "How it is verified" of "The folded site frequency
  spectrum" and of "The inbreeding coefficient F_IS". Serves deliverable
  2. `docs/reports/diversity-method/check_dadi.py` and `check_allel.py`
  are what it grows from; they read `panel_counts.tsv`, which this script
  has to write itself.
- [x] 1.3 `enumerate_private.py` and its output, from "How it is
  verified" of "The private alleles". Serves deliverable 3.
  `docs/reports/diversity-method/check_by_enumeration.py` is what it
  grows from.
- [ ] 1.4 The page beside the scripts. Serves deliverable 4.

**What could go wrong.** `scikit-allel` re-enables the global interpreter
lock of the process that imports it, which the project's Python, a free
threading build, warns about. That is why task 1.2 puts it in the
environment it makes and not in the development dependencies of
`pyproject.toml`: no test of popnei imports it, and the reference script
runs in a process of its own. If a later task finds a reason to import it
inside pytest, the warning and what it does to the other tests is measured
first.

## Work package 2: the counts with no draw, through the four layers

**What it gives.** A user calls `calc_pop_diversity(variants, pops=...,
stats=..., min_num_individuals=...)` in Python and gets, per population,
the alleles it called and the private ones as a total and a mean, the
variable variants as a total and a ratio, F_IS, the counts of the variants
that served each and `pass_stats`; in TypeScript, `calcPopDiversity`.
`num_called_alleles` is not taken yet and the folded spectrum is not
there, which is work package 3.

**Deliverables.**

1. The module, its options and its result. Check: `cargo test -p popnei
   --lib -- diversity:: --list` prints 14 tests or more, where it prints
   `0 tests` today.
2. The worked example, the six variants of five diploid individuals of
   "How it is verified" of "The number of alleles". Check: named cargo
   tests assert, for `pop1` and `pop2`, the alleles called 9 and 8 with
   the means 2.25 and 2, the private alleles 2 and 1 with the means 0.5
   and 0.25, the variable variants 3 and 2 with the ratios 0.75 and 0.5,
   and F_IS 0 and 0.3478260870.
3. Which variants count. Check: cargo tests assert that variant 4, with
   nothing called, and variant 6, with one called allele in `pop1`, are
   out of every count at `min_num_individuals` 1, and that a
   `min_num_individuals` of 0 does not put a variant with no called
   allele into the totals, which is "Its Python function" of the spec.
4. The errors of "The Rust interface". Check: cargo tests for a
   population with no individual, an index that is not an individual, an
   individual asked for twice, and a pass with no variant.
5. The Python function and its result. Check: `uv run pytest
   tests/test_diversity.py` passes with 8 tests or more, among them one
   that reads the files of work package 1 and asserts the panel's 2373,
   2377 and 2384, its 0, 0 and 1, its 1173, 1177 and 1184 variable
   variants, and its F_IS of -0.0127584868, -0.0181107131 and
   -0.0184585832 within 1e-12 relative.
6. The two cases of "The cases" that a user can reach. Check: a pytest
   test asserts that one population gets private alleles equal to its
   alleles called, and a cargo test asserts the 0.5 of two populations
   holding the one genotype `0/1` at a draw of one allele, which is the
   overlap the spec warns about; the second is a test of work package 3,
   named here because the case belongs to this item.
7. The TypeScript function. Check: `npm test` in `js/popnei` passes with
   `js/popnei/test/diversity.test.ts` among its files, asserting the
   panel's three totals and three F_IS values under node.
8. The agreement with `stats`. Check: a pytest test runs
   `calc_pop_diversity` and `calc_per_var_distribs` on the panel with the
   same `pops` and `min_num_individuals` and asserts that the variable
   variants equal the `num_variable` of `poly_vars_ratio`, at any
   `poly_threshold`.

**What it stands on.** Work package 1, for the files deliverable 5 reads.
Outside the plan: `Pops` and `Pops::from_names` of
`crates/popnei/src/stats.rs`, `count_alleles_of`, `AlleleCounts` and
`the_major_allele` of `crates/popnei/src/variant.rs`, and `ObsHet::of_var`
and `ExpHet::of_var` of `crates/popnei/src/stats.rs`, all public today.

**Tasks.**

- [ ] 2.1 `crates/popnei/src/diversity.rs` with `DiversityStats`,
  `DiversityOptions` and `PopDiversity`, declared in
  `crates/popnei/src/lib.rs`, with the pass that walks the blocks and
  counts nothing yet but the variants that count for each population, and
  the error cases added to `crates/popnei/src/error.rs`. From "The pass
  and its function" and "The Rust interface". Serves deliverables 1, 3
  and 4.
- [ ] 2.2 The alleles called and the variable variants, from "What it
  gives" of "The number of alleles" and of "The variable variants". Needs
  2.1. Serves deliverables 1, 2 and 8. A wrong count is silent: its own
  commit, guarded by deliverables 2, 5 and 8.
- [ ] 2.3 The private alleles, which need the counts of every population
  at one variant at once and the second divisor, from "What it gives" of
  "The private alleles". Needs 2.2. Serves deliverables 1, 2 and 6. A
  wrong count is silent: its own commit, guarded by deliverables 2 and 5.
- [ ] 2.4 F_IS, built from `ObsHet::of_var` and `ExpHet::of_var` of the
  `stats` module so that the two specs cannot disagree about a
  heterozygosity, from "What it gives" of "The inbreeding coefficient
  F_IS", whose rule for which variants enter it is there too. Needs 2.1.
  Serves deliverables 1, 2 and 5. A wrong number is silent: its own
  commit, guarded by deliverables 2 and 5.
- [ ] 2.5 `crates/popnei-python/src/diversity.rs`, its entry in that
  crate's `lib.rs`, `python/popnei/diversity.py` with the result
  dataclass and the enum, its exports in `python/popnei/__init__.py`, and
  `tests/test_diversity.py`. From "Its Python function" and "The cases".
  Needs 2.2, 2.3 and 2.4. Serves deliverables 5, 6 and 8.
- [ ] 2.6 `crates/popnei-js/src/diversity.rs`, its entry in that crate's
  `lib.rs`, `js/popnei/src/diversity.ts`, its exports in `node.ts` and
  `web.ts`, and `js/popnei/test/diversity.test.ts`. From the TypeScript
  paragraph of "Its Python function". Needs 2.5, whose Python names it
  mirrors, and `npm run build` having been run in `js/popnei` as "What
  has to be in place" says. Serves deliverable 7.

**What could go wrong.** The private alleles are the only statistic here
that needs every population's counts at one variant at the same time, so
the pass holds an array of the populations by the alleles of the row while
the other four read one population at a time. "How it runs" of "The pass
and its function" says what is kept and what is not; a version that keeps
it per block instead of per row grows with the block and fails no test.

## Work package 3: the draw, through the same four layers

**What it gives.** The same call takes `num_called_alleles` and gives, per
population, the alleles called and the private ones and the variable
variants as they would be in a draw of that many called alleles, and the
folded site frequency spectrum projected to it. When this is done the
module is what the spec describes.

**Deliverables.**

1. The draw arithmetic on its own. Check: cargo tests assert the eight per
   variant values of the table at the end of "How it is verified" of "The
   number of alleles", which `vegan` gave: 1.5, 1.5, 2 and 1 for `pop1`
   and 1, 1, 2 and 1.5333333333 for `pop2` at a draw of 2, whose means
   are 1.5 and 1.3833333333.
2. The standardized alleles called and variable variants on the panel.
   Check: a pytest test reads the file of deliverable 1 of work package 1
   and asserts 1.9283948650, 1.9219209943 and 1.9197370844 and the three
   ratios one less than those, within 1e-12 relative, which is the
   identity "How it is verified" of "The variable variants" gives.
3. The standardized private alleles. Check: cargo tests assert the 18
   pairs of the file of deliverable 3 of work package 1, and a pytest
   test asserts the panel's 0.0112196177, 0.0099715392 and 0.0089014974.
   The two properties of that item, that the value is at most the
   standardized alleles called and that a population against a copy of
   itself gets 0, are cargo tests too.
4. The folded spectrum. Check: a pytest test reads the file of
   deliverable 2 of work package 1 and asserts all 33 values of the three
   populations within 1e-12 relative, and that each column sums to 1200.
5. The arguments and the refusals in the layers above. Check: a pytest
   test and a TypeScript test assert that asking for the spectrum without
   `num_called_alleles` is refused, that a `num_called_alleles` below 2
   is refused, and that a `num_called_alleles` above every population's
   called alleles gives NaN everywhere in the draw and leaves the totals
   and F_IS as they were.
6. Which variants are in the draw. Check: cargo tests assert, on the
   worked example at `num_called_alleles` 4, that `pop1` keeps its four
   variants and `pop2` three, that the standardized values are 2.25 and
   2.3111111111 for the alleles, 0.3333333333 and 0.3111111111 for the
   private ones and 0.75 and 0.6444444444 for the variable variants, and
   that the spectra are 1, 3, 0 and 1.0666666667, 1.5333333333, 0.4.

**What it stands on.** Work packages 1 and 2.

**Tasks.**

- [ ] 3.1 The chance that a draw of `g` of `c` misses an allele called
  `n` times, computed as the product of `g` factors as "The Rust
  interface" says, with its tests, in
  `crates/popnei/src/diversity.rs`. From "The Rust interface" and the
  table of "How it is verified" of "The number of alleles". Serves
  deliverable 1. Every standardized value rests on it and a wrong value
  is silent: its own commit, guarded by deliverables 1, 2, 3 and 4.
- [ ] 3.2 The standardized alleles called and variable variants, and the
  second count of variants, from "What it gives" of those two items.
  Needs 3.1. Serves deliverables 2 and 6.
- [ ] 3.3 The standardized private alleles, from "What it gives" of "The
  private alleles". Needs 3.1. Serves deliverables 3 and 6. Its own
  commit for the same reason as 2.3.
- [ ] 3.4 The folded spectrum, from "What it gives" of "The folded site
  frequency spectrum", whose range for the count of the rarer allele is
  the part to read twice. Needs 3.1. Serves deliverables 4 and 6. Its own
  commit: a wrong bin moves a number and fails nothing else.
- [ ] 3.5 `num_called_alleles` and `FOLDED_SFS` in both binding crates,
  in `python/popnei/diversity.py` and in `js/popnei/src/diversity.ts`,
  with their refusals and their tests. Needs 3.2, 3.3 and 3.4. Serves
  deliverables 2, 3, 4 and 5.

**What could go wrong.** The count of the rarer allele can be far above
the draw size, 48 against 20 on the panel, and the range that the sum over
it runs on is where the spec was wrong before it was reviewed. A version
that runs the count from 0 to the rarer allele's count asks for a binomial
coefficient of a negative argument and a bin below 0, which panics in Rust
rather than giving a wrong number, so the panel test catches it; a version
that clamps instead of skipping gives a wrong number and does not panic,
and deliverable 4 is what catches that.

## Work package 4: the measurement and the browser

**What it gives.** The speed target that "Speed" of the spec leaves open,
measured and written into the spec, and the module running in a browser
and under node.

**Deliverables.**

1. The measurement, on the panel and on 100000 variants of 1000
   individuals, of the pass with every statistic and with the spectrum
   alone, at draws of 20 and 200. Check: a report under `docs/reports/`
   with the dataset, the machine, the date and the numbers, and the
   commit that writes the target into "Speed" of the spec before the code
   of any later task.
2. Whether the hypergeometric weights of a variant are computed once and
   shared by the bins, which "Speed" names as the thing to decide.
   Check: the report holds both timings, and the change is kept only if
   the measurement says so, as the `performance-review` skill asks.
3. The wasm build and node. Check: `cargo wasm-check` finishes, and in
   `js/popnei`, `npm run build && npm test` passes with the diversity
   tests among them.
4. The wheel under pyodide. Check: `bash
   scripts/build_pyodide_wheel.sh && node tests/pyodide/smoke.mjs` exits
   0 with a call of `calc_pop_diversity` in the smoke test.

**What it stands on.** Work packages 2 and 3.

**Tasks.**

- [ ] 4.1 The measurement and the report, and the commit that writes the
  target into the spec. Serves deliverables 1 and 2.
- [ ] 4.2 The browser and pyodide checks, and the call added to the smoke
  test. Serves deliverables 3 and 4. Can run side by side with 4.1.

**What could go wrong.** The spectrum does one hypergeometric term per
count of the rarer allele per variant and per population, which at a draw
of 200 over a million variants and three populations is 6·10⁸ terms; if
the measurement finds that it dominates the pass, the shared weights of
deliverable 2 are the first thing to try, and anything beyond that is a
performance review of its own and not a task of this plan.

## How the whole plan is checked

The sum of the work packages, and one thing besides: the spec says every
number of it but the standardized private alleles comes from a program
outside popnei. So, at the end, `Rscript
tests/reference/diversity/make_reference.R` and `uv run python
tests/reference/diversity/make_reference.py` are run again from a clean
tree, and `git status --short tests/reference/diversity` has to show
nothing. A reference file that the work quietly changed would make every
comparison of work packages 2 and 3 pass against a number popnei itself
moved.
