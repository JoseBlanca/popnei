# Plan: the per variant and the per individual statistics, and the filter of individuals

22 September 2026. Approved by the owner on 22 September 2026, in chat,
when he asked for the handoff that runs it. Under way since 22 September
2026 on the branch `plan/stats`, with its report in
`docs/reports/stats.md`. This plan builds the
`stats` module of popnei, whole: the populations a statistic is
calculated for, the one pass that gives the distributions of five per
variant statistics, `calc_per_var_distribs`, and the per individual
statistics, `calc_per_individual_stats`, each through the three layers
to the Python function and its comparison with pyNei; and, before them,
the filter of individuals, `Variants.filter_individuals`, the step that
keeps the named individuals and that those statistics are the first to
need. It builds from:

- `docs/specs/stats.md`, every item of it, "The Rust interface" and
  "Speed";
- `docs/specs/filters.md`, the item "The filter of individuals" and, in
  "The Rust interface", `PassStep`, `resolve_individuals`,
  `IndividualsReader`, `Block::retain_individuals`, the new `chain_of`
  and the new `refuse_a_second_filter_of_a_kind`;
- `docs/specs/variant.md`, "The counts of one variant", whose two
  functions the counts over a population sit beside.

Neither spec has an open point: the owner decided the eight they had on
22 September 2026, and each is written into its spec.

The branch is `plan/stats`, from `spec/stats` at 9b19636 or later, in the
worktree `.claude/worktrees/stats`. `spec/stats` is from `main` at
7d8366f and holds the two specs, the glossary entries and
`tests/reference/stats/`, the reference script with the reports of plink2
and bcftools it stores. Two branches that are not in `main`,
`plan/dists-kosman` and `plan/pca`, hold the specs of the distances and
of the principal components that the stats spec points to; nothing here
needs their code.

## In and out

In:

- `PassStep` in the core, and both binding crates keeping their steps as
  a list of it.
- The filter of individuals, in the core, in both binding crates and in
  both packages, with `Block::retain_individuals`.
- `Pops` and the counts of one variant over a population, in the core.
- `calc_per_var_distribs`, its five statistics, `HistBins`,
  `StatsDistrib`, `PolyVarsStats` and `PerVarDistribs`, in the core, in
  both binding crates and in both packages.
- `calc_per_individual_stats` and `PerIndividualStats`, the same.
- The measurement of "Speed" of the stats spec, and the wheel for pyodide
  with the two calculations in it.

Out, with where it goes:

- Jost's D between populations, the linkage disequilibrium per
  population and the kinship, which take their populations from `Pops`:
  their own specs and plans.
- The filter by linkage disequilibrium and the variants of a region:
  later items of `docs/specs/filters.md`, not written.
- One pass for both calculations when a user wants both: nothing in the
  specs asks for it, as "What it gives" of the per individual statistics
  says.
- The read ahead thread of section 3 of the architecture: no spec
  describes it, as `docs/plans/dists-kosman.md` says, and both
  calculations here take any reader.
- `Variants.from_gt_array`, which the stats spec names among what it
  depends on: it has no code and no plan, as the owner decided for the
  Kosman plan, so the pytest tests of the worked examples write their
  genotypes as small VCF files with the `write_vcf` fixture of
  `tests/conftest.py`, as the tests of the filters do.

## What has to be in place

- The branch `plan/stats` at 9b19636 or later, in the worktree above.
  Run on `spec/stats` at 9b19636 on 22 September 2026: `cargo fmt --all
  --check` exit 0; `cargo clippy --workspace --all-targets -- -D
  warnings` no warning; `cargo test --workspace` `306 passed`, 2
  ignored; `cargo wasm-check` finished; ruff `18 files already
  formatted` and `All checks passed!`; `uv run maturin develop && uv run
  pytest` `174 passed`; in `js/popnei`, after `npm install`, `npm run
  build` and `npm test` 126 tests passed, 0 failed. These are the counts
  that "more than" is counted from below. `bash
  scripts/build_pyodide_wheel.sh && node tests/pyodide/smoke.mjs` was not
  run when the plan was written; in a new worktree it needs `npm
  install` in `tests/pyodide` first.
- pyNei at ef0ca6e, which `pyproject.toml` has: `uv run python -c "from
  pynei import calc_per_var_distribs, calc_per_sample_stats,
  filter_samples, load_vars"` exits with 0, and
  `/Users/jose/devel/pynei/test/gwas_reference/sim_missing.vars` is
  there, the panel that the reference script reads.
- plink2 v2.0.0-a.7.7 at `/opt/homebrew/bin/plink2` and bcftools 1.24 at
  `/opt/homebrew/bin/bcftools`. `uv run python
  tests/reference/stats/make_reference.py` was run on 22 September 2026:
  it writes `panel.vcf.gz` and the reports of both programs beside
  itself, checks the literals of the stats spec against them and exits
  with 0, and `git status --short tests/reference/stats` shows nothing
  after it. The commands of "How it is verified" of the filter of
  individuals were run on `many.vcf` the same day and give the 423
  variants, the five positions and the three names in the order of the
  spec.
- The file of the measurement, outside the repository:
  `/Users/jose/devel/popnei-bench/big.vars`, 81356714 bytes, the vars
  file popnei wrote of `big.vcf`, 100000 variants of 1000 individuals,
  and pyNei's own vars file of it, which task 6.1 writes beside them if
  it is not there. Both `big.vcf` and `big.vars` are there.

Every layer exists, so the five commands of "Before the work is called
done" of the `coding` skill run for every task, with `cargo wasm-check`,
and `npm run build && npm test` for a task that touches `crates/popnei-js`
or `js/popnei`.

The checks below fail today: `cargo test -p popnei --lib -- stats::
--list` prints `0 tests`, and so does `cargo test -p popnei --lib --
retain_individuals resolve_individuals IndividualsReader --list`;
`tests/test_filter_individuals.py`, `tests/test_stats.py`,
`js/popnei/test/filter_individuals.test.ts` and
`js/popnei/test/stats.test.ts` do not exist.

Tasks run one after another unless the plan says that two can run side
by side: every task builds the same workspace in the same worktree.

A task whose failure would be silent, a wrong number and not a crash or
a failing test, has a commit of its own and names the deliverables that
guard it: the checks whose literals or comparisons would fail if that
number were wrong, so that the commit that moved a number can be found
later.

## Work package 1: the steps of a pass as one enum

**What it gives.** Nothing a user sees. The chain of a pass is built
from a list of `PassStep`, in the core and in both binding crates, so
that work package 2 adds its step to that list and to `chain_of`
without touching the three threshold filters again. It changes no
result: every test that exists passes untouched.

**Deliverables.**

1. `PassStep`, `chain_of` and `refuse_a_second_filter_of_a_kind` over
   it. Check: `cargo test -p popnei --lib -- PassStep --list` prints 2
   tests or more, `kind` of each variant and the chain of the three
   threshold filters as a list of steps giving what it gave; `cargo test
   --workspace` passes with the 306 tests that were there, changed only
   where they call the two functions.
2. Both binding crates keep their steps as `PassStep`. Check: `grep -rn
   "VarFilteringCriterion" crates/popnei-python/src crates/popnei-js/src`
   finds the criterion of a threshold filter only where one is built
   from its argument;
   `uv run pytest` 174 passed and `npm test` 126 passed, untouched.

**What it stands on.** The specs.

**Tasks.**

- [x] 1.1 `PassStep` with its `kind`, and `chain_of` and
  `refuse_a_second_filter_of_a_kind` taking steps, in
  `crates/popnei/src/filters.rs`; the `Step` of
  `crates/popnei-python/src/steps.rs` and its counterpart in
  `crates/popnei-js/src/steps.rs` become the core's `PassStep`, or hold
  one. From "The Rust interface" of `docs/specs/filters.md`, the
  paragraph "The steps of a pass" and the two functions after it. The
  variant `KeepIndividuals`, the step of the filter of individuals, is
  declared here and built in work package 2; until then `chain_of` gives
  an error for it, of the kind that marks a defect of popnei, since no
  user can add that step yet. Serves deliverables 1 and 2.

**What could go wrong.** The error of a second filter of one kind
carries two thresholds today, `Error::VarFilterOfAKindThatIsSet` of
`crates/popnei/src/error.rs`; the spec keeps that case for the threshold
filters and adds a case for the filter of individuals in work package 2,
so this task does not widen the one that exists.

## Work package 2: the filter of individuals

**What it gives.** A user calls `variants.filter_individuals(["ind05",
"ind00", "ind49"])` in Python, or `filterIndividuals` in TypeScript, and
every consumer after that gets the genotypes of those three, in that
order, at every variant, with `variants.individuals` saying so; a filter
after it in the steps counts over the three.

**Deliverables.**

1. `Block::retain_individuals`. Check: `cargo test -p popnei --lib --
   retain_individuals --list` prints 4 tests or more: the worked example
   of "How it is verified" of the filter, i5 and i1 kept from the six
   variants of the threshold filters' example, with the two rows the
   spec gives; an index at or beyond the individuals, an index twice and
   no index, each leaving the block as it was; and a block with variants
   and no genotypes.
2. `resolve_individuals` and `IndividualsReader`. Check: `cargo test -p
   popnei --lib -- resolve_individuals IndividualsReader --list` prints
   8 tests or more, with every case of "How it is verified" of the
   filter: the three names giving 5, 0 and 49 and the three errors at
   `resolve_individuals`; at `next_block` over a `VcfReader` on
   `many.vcf`, blocks of 3 individuals with the three names in the order
   of the argument, the genotypes of each being the column of the
   source, 500 variants and an empty `filtering_stats`, the counts of
   the filters of the chain that `pass_stats.filtering` is built from;
   the missing data
   filter at 0 over it giving 423 variants with the five positions
   first and the counts 500 and 423, and under it 26; and `chain_of`
   building the reader from a `KeepIndividuals` step and refusing a
   second one.
3. Python. Check: `uv run pytest tests/test_filter_individuals.py`
   passes, where today the file is not there, with the comparison with
   pyNei of "How it is verified", column by name, the 423 with the
   missing data filter after it and its `pass_stats.filtering` of one
   entry, `pass_stats` being the counts of the pass that every result
   and every `iter_blocks` carries, and every test of the step the same
   part names: `steps`,
   `individuals` and `num_individuals` after the call, and the four
   `ValueError`s.
4. TypeScript. Check: `npm test` gives 0 failed and
   `test/filter_individuals.test.ts` asserts the three names, the 423
   and the `Error` of an unknown name.
5. The wheel of pyodide still builds and its smoke test exits with 0.

**What it stands on.** Work package 1.

**Tasks.**

- [x] 2.1 `Block::retain_individuals` in `crates/popnei/src/block.rs`,
  with the two passes of "How it runs" of the filter, and its tests.
  From "How it runs" and "The Rust interface" of `docs/specs/filters.md`.
  Serves deliverable 1. A wrong gather is silent: its own commit, guarded
  by deliverables 1 and 3.
- [x] 2.2 `resolve_individuals`, `IndividualsReader`, the `KeepIndividuals`
  step in `chain_of`, and the four new cases of the error, in
  `crates/popnei/src/filters.rs` and `error.rs`, with the arm of each
  binding crate's `errors.rs`, and the cargo tests of deliverable 2.
  From "What it gives", "How it runs" and "How it is verified" of the
  filter, and "The Rust interface". It mirrors `filter_samples` of
  `pynei/var_filters.py`. Needs 2.1. Serves deliverable 2.
- [x] 2.3 The Python side: `filter_individuals` in
  `crates/popnei-python/src/steps.rs`, which resolves the names at the
  call against the individuals of the source and refuses a second
  filter of the kind; `individuals` and `num_individuals` of
  `python/popnei/variant.py` after the step, the method and its
  docstring, the `repr`; and `tests/test_filter_individuals.py`. From
  "In Python and in TypeScript" and "What pyNei asserts" of the filter.
  Needs 2.2. Serves deliverables 3 and 5.
- [ ] 2.4 The TypeScript side, the same in `crates/popnei-js` and
  `js/popnei`, with `test/filter_individuals.test.ts`. Needs 2.2. Serves
  deliverable 4. After 2.3 in the plan's worktree.

**What could go wrong.** `individuals` of the Python `Variants` is read
from the header once, when the handle is built; after this step it is
what the next pass gives, so the package reads it through the steps.
The reader of the filter keeps the block size of its source, worked out
from the individuals of the source, as "How it runs" says: a test that
counts genotypes per block must not expect the size for three
individuals.

## Work package 3: the populations and the counts over one

**What it gives.** Work package 4 gets `Pops`, the names of a user's
populations as indices among the individuals a pass gives, and
`count_gts_of` and `count_alleles_of`, the two counts of one variant over
the individuals of one population. It stops inside the core because
neither has a function in Python or in TypeScript; the comparison with
pyNei runs in work package 4, through the statistics.

**Deliverables.**

1. `Pops`. Check: `cargo test -p popnei --lib -- stats::pops --list`
   prints 5 tests or more, the cases of "How it is verified" of "The
   populations": the indices of the two populations of the worked
   example, the four errors, an individual in two populations, and
   `all` giving one population named `pop` of every individual.
2. The two counts. Check: `cargo test -p popnei --lib -- count_gts_of
   count_alleles_of --list` prints 4 tests or more: variant 1 of the
   worked example over pop1 and over pop2 with the numbers of "How it is
   verified" of "The counts of one variant over a population", an index
   beyond the row, and, for each, that the population of every
   individual in the order of the source gives what `count_gts` and
   `count_alleles` give.

**What it stands on.** Work package 2, whose `resolve_individuals`
`Pops::from_names` calls.

**Tasks.**

- [ ] 3.1 `Pops` and its errors in `crates/popnei/src/stats.rs`, a new
  module with nothing else yet, and `count_gts_of` and
  `count_alleles_of` in `crates/popnei/src/variant.rs`, with their
  tests. From "The populations", "The counts of one variant over a
  population" and "The Rust interface" of `docs/specs/stats.md`; they
  mirror `_calc_pops_idxs` of `pynei/utils_pop.py` and the population
  slices of `_calc_obs_het_per_var` and `_count_alleles_per_var` of
  `pynei/gt_counts.py`. Serves deliverables 1 and 2. A wrong count is
  silent: its own commit, guarded by deliverable 2 and by the
  comparisons of work package 4.

## Work package 4: the per variant distributions, through the three layers

**What it gives.** A user calls `calc_per_var_distribs(variants,
pops=..., min_num_individuals=...)` in Python and gets, per population,
the mean and the histogram of the observed heterozygosity, the maf, the
expected heterozygosity, plain and unbiased, and the counts of the
polymorphism ratio, with `pass_stats`; in TypeScript,
`calcPerVarDistribs`.

**Deliverables.**

1. The bins and the value of one variant. Check: `cargo test -p popnei
   --lib -- stats::hist stats::obs_het stats::maf stats::exp_het --list`
   prints 16 tests or more: the edges of 4 and of 40 linear bins and of
   the logarithmic example of `test_hist.py`, the bin of a value on an
   interior edge, on the last edge and outside the range, and the three
   errors of `HistBins`; at `ObsHet::of_var`, `Maf::of_var` and
   `ExpHet::of_var`, every per variant value of the worked example of
   the pass over both populations, the `var0000` literals of the panel
   of each item, the first three variants of `many.vcf` of the maf and
   of the observed heterozygosity, the three variants of the expected
   heterozygosity's own worked example, plain and unbiased, its
   tetraploid variant giving 0.995960 unbiased, a population below
   `min_num_individuals` and one with nothing called giving none, and
   the constructors refusing 0 and 256.
2. The pass. Check: `cargo test -p popnei --lib -- stats::distribs
   --list` prints 8 tests or more, made at `calc_per_var_distribs` over
   a reader written in the test and over the VCF reader: the means, the
   histograms and the polymorphism counts and ratios of the worked
   example of the pass, both populations and no `pops`, with blocks of
   6 and of 2; the means of the expected heterozygosity's worked
   example; the nine polymorphism counts of the panel and the six of
   `many.vcf` from "How it is verified" of the polymorphism ratio; the
   `min_num_individuals` of 20 leaving every statistic without a value
   on a population of 15; a reader with no variant giving the error;
   and the same numbers in rayon pools of 1 and of 4 threads.
3. Python. Check: `uv run pytest tests/test_stats.py -k per_var` passes,
   where today the file is not there, with: the comparison with pyNei
   of "How it is verified" of the pass, on the panel and on `many.vcf`,
   with the two populations of each, the default histogram, both
   thresholds, and once with no `pops`, the counts equal and the means
   within 1e-12; the default of `min_num_individuals` on a population
   of 15; every test of pyNei that "What pyNei asserts, and the size of
   the blocks" names, with a string in `stats` a `TypeError`; the
   `ValueError`s of "The populations" and of a pass with no variant;
   the order of the populations being that of the keys; and
   `pass_stats` with the variants of the panel and the counts of a
   filter put before the call.
4. TypeScript. Check: `npm test` gives 0 failed and `test/stats.test.ts`
   runs the literals of `p0` of the panel that each item gives, through
   `calcPerVarDistribs`, and the `Error` of an unknown individual.

**What it stands on.** Work package 3.

**Tasks.**

- [ ] 4.1 `HistBins`, `ObsHet`, `Maf` and `ExpHet` with their
  constructors, in `crates/popnei/src/stats.rs`, and the tests of
  deliverable 1. From "In Python and in TypeScript" of the pass for the
  bins, and "What it gives", "Missing genotypes" and "How it is
  verified" of each of the three items, with "The Rust interface". They
  mirror `_prepare_bins` of `pynei/utils_stats.py`,
  `_calc_obs_het_per_var` and `_calc_maf_per_var` of
  `pynei/gt_counts.py`, and `_calc_exp_het_per_var` and
  `_calc_unbiased_exp_het_per_var` of `pynei/diversity.py`. Serves
  deliverable 1. A wrong value is silent: its own commit, guarded by the
  literals of deliverable 1 and the comparison of deliverable 3.
- [ ] 4.2 `calc_per_var_distribs`, `PerVarDistribsConfig`, `StatsDistrib`,
  `PolyVarsStats` and `PerVarDistribs`, in the same file: the loop over
  the blocks with rayon over the rows and the accumulators of "How it
  runs" of the pass, the polymorphism counts of "What it gives" of its
  item, and the error of a pass with no variant. From "How it runs" and
  "How it is verified" of the pass, "The polymorphism ratio" and "The
  Rust interface". Needs 4.1. Serves deliverable 2.
- [ ] 4.3 The Python side: the function of `crates/popnei-python` that
  builds the chain with `chain_of`, the `Pops` with `Pops::from_names`
  against the individuals of that chain, runs the pass and gives back
  the distributions and the counts of the pass, as the writer of
  `crates/popnei-python/src/vars.rs` does, which builds the chain, runs
  its pass inside one call of the core and reads the counts from the
  chain afterwards; `PerVarStat`, the enum of the five statistics,
  `StatsDistrib`, `PolyVarsStats`,
  `PerVarDistribs` and `calc_per_var_distribs` in
  `python/popnei/stats.py`, exported from `popnei`; and the tests of
  deliverable 3 in `tests/test_stats.py`. From "In Python and in
  TypeScript" and "What pyNei asserts, and the size of the blocks" of
  the pass, "In Python" of each statistic, and "In Python and in
  TypeScript" of "The populations". Needs 4.2. Serves deliverable 3.
- [ ] 4.4 The TypeScript side: the same in `crates/popnei-js`,
  `calcPerVarDistribs` and its result in `js/popnei/src/stats.ts`, and
  `test/stats.test.ts`. From the TypeScript paragraphs of the same
  parts. Needs 4.2. Serves deliverable 4. After 4.3 in the plan's
  worktree.

**What could go wrong.** The edges of the bins have to be the bits numpy
computes, `start + i * step`, or a value on an edge falls on the other
side of it in the two libraries and a histogram count is off by one;
the spec says how, and the test of the edges compares them exactly. The
sums of the means are added in an order that depends on the threads, so
the cargo test of the pools compares the means within 1e-12 and the
counts exactly. The pytest comparison runs pyNei twice per dataset, for
its two values of `unbiased_exp_het`, against popnei's one call.

## Work package 5: the per individual statistics, through the three layers

**What it gives.** A user calls `calc_per_individual_stats(variants)` in
Python and gets the missing rate and the heterozygosity rate of every
individual, with `pass_stats`; in TypeScript, `calcPerIndividualStats`.

**Deliverables.**

1. The calculation. Check: `cargo test -p popnei --lib --
   stats::per_individual --list` prints 5 tests or more, made at
   `calc_per_individual_stats`: the worked example of "How it is
   verified" of the item with blocks of 6 and of 2; the literals of
   `s000` and `s001` on the panel and of `ind00` and `ind01` on
   `many.vcf` read with the VCF reader; an individual with no called
   genotype; a reader with no variant; and the same numbers in rayon
   pools of 1 and of 4 threads.
2. Python. Check: `uv run pytest tests/test_stats.py -k per_individual`
   passes, with the comparison with pyNei of "How it is verified" of
   the item, on the panel and on `many.vcf`, the missing rates equal and
   the heterozygosity rates rescaled as the spec says; pyNei's
   `test_filter_missing` with popnei's rates over called genotypes; the
   `ValueError` of a pass with no variant; and `pass_stats`.
3. TypeScript. Check: `npm test` gives 0 failed and `test/stats.test.ts`
   runs the `s000` and `s001` literals through `calcPerIndividualStats`.

**What it stands on.** Work package 2, for the comparison after a filter
of individuals, and nothing of work packages 3 and 4. In the plan's
worktree it runs after them; an orchestrator with a second worktree
could run it beside them, since they touch different files.

**Tasks.**

- [ ] 5.1 `calc_per_individual_stats` and `PerIndividualStats` in
  `crates/popnei/src/stats.rs`, with the tests of deliverable 1. From
  "What they give", "How it runs" and "How it is verified" of "The per
  individual statistics" and "The Rust interface". It mirrors
  `calc_per_sample_stats` of `pynei/sample_stats.py`. Serves deliverable
  1. A wrong count is silent: its own commit, guarded by deliverables 1
  and 2.
- [ ] 5.2 The Python side, in `crates/popnei-python` and
  `python/popnei/stats.py`, with the tests of deliverable 2. From "In
  Python and in TypeScript" and "What pyNei asserts" of the item. Needs
  5.1. Serves deliverable 2.
- [ ] 5.3 The TypeScript side, in `crates/popnei-js` and
  `js/popnei/src/stats.ts`. Needs 5.1. Serves deliverable 3. After 5.2
  in the plan's worktree.

## Work package 6: speed and the browser

**What it gives.** Whether the two calculations reach the numbers of
"Speed" of the stats spec, in a report, and the wheel for pyodide with
both in it.

**Deliverables.**

1. The measurement. Check: `docs/reports/stats-measurement.md` has, on
   `big.vars`, the time of `calc_per_var_distribs` with the five
   statistics and no `pops`, with 4 populations of 250, and of
   `calc_per_individual_stats`, each through `open_vars`, on one thread
   and on 18 cores, best of 5 runs with the load average of each; the
   pass over the file alone through `iter_blocks` with the genotypes as
   the only field, run the same day; pyNei's times on its own vars file
   of `big.vcf`, with `num_threads` 1 and 6, run again that day; and,
   for each of the four numbers of "Speed", whether it is met. A number
   not met is a finding for the owner and not a task of this plan.
   "Speed" of the spec gets popnei's numbers beside pyNei's, with what
   they were measured on.
2. The browser. Check: `bash scripts/build_pyodide_wheel.sh && node
   tests/pyodide/smoke.mjs` exits with 0, with the smoke test extended to
   the worked example of the pass through `calc_per_var_distribs` and
   `calc_per_individual_stats`.

**What it stands on.** Work packages 4 and 5.

**Tasks.**

- [ ] 6.1 The measurement and its report, with the machine and the load
  average as `docs/reports/filters-measurement.md` gives them, and the
  paragraph of "Speed". From "Speed" of the stats spec. Serves
  deliverable 1. It changes no code of the library.
- [ ] 6.2 The pyodide smoke test with the worked example, in
  `tests/pyodide/smoke.mjs`, and the build. Serves deliverable 2. Can run
  side by side with 6.1.

**What could go wrong.** The number to reach was derived from the cost
of the missing data filter and not measured on a trial: if the first
measurement is over it, the orchestrator reports the number and where
the time goes, and whether to work on it is the owner's, through a
performance review.

## How the whole plan is checked

The sum of its work packages, and the counts of "What has to be in
place" grown by the tests the plan adds: `cargo test --workspace` more
than 306, `pytest` more than 174, `npm test` more than 126, 0 failed
everywhere, from a clean clone of the branch, with the work report in
`docs/reports/stats.md`.
