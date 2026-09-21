# Plan: the three threshold filters and the counts of a pass

21 September 2026. Approved by the owner on 21 September 2026 and under
way, with its work report in `docs/reports/filters.md`. This plan builds
the filters that keep the variants of a dataset by their missing rate, by
their major allele frequency and by their observed heterozygosity, as
steps that a user puts on a `Variants`, and the counts that every result
carries of how many variants each filter was given and kept. It builds
from:

- `docs/specs/filters.md`, the three filters, the counts of a pass and
  "Speed";
- `docs/specs/variant.md`, "The counts of one variant", and the
  `Variants` handle and `PassStats` of "What a Python and a TypeScript
  user see";
- `docs/specs/block.md`, `filtering_stats` in the `BlockReader` trait and
  the `Blocks` that `iter_blocks` returns;
- `docs/specs/io_vars.md`, `write_vars` returning `VarsWritten`, in "Its
  Python and TypeScript functions" of the writer.

The owner settled the API of these specs on 21 September 2026 and none has
an open point. He gave the order of the work packages with the task of
writing this plan, and said the same day that the branch starts from the
vars file work and not from `main`.

The branch `plan/filters` stands on two branches that are not in `main`:
`spec/filters` at fee97bf, which has the specs above, and
`plan/vars-file` at 2f2577e, merged into it at 2c2a630 with no conflict,
which has the vars file reader, `write_vars` and `open_vars` in both
packages. Of `docs/plans/vars-file.md` only task 4.1, its bench, was not
done at that commit. The owner had `plan/vars-file` merged into `main` on 21 September 2026,
at d1d6997, and said to merge that `main` into this branch, which is
d4e8458, made between tasks 3.3 and 3.4. `spec/filters` is still not in
`main`, and nothing of this plan reaches it before the owner's order.

## In and out

In:

- `filtering_stats` in the `BlockReader` trait and in every reader: the
  VCF reader, the vars file reader, `reblock`, the boxed reader and the
  readers written in the tests.
- `PassStats` in what every consumer that exists returns: the `Blocks` of
  `iter_blocks`, and the `VarsWritten` of `write_vars`, in Python and in
  TypeScript.
- A `Variants` that holds its steps and shows them in `steps` and, in
  Python, in its `repr`.
- The counts of the genotypes and of the alleles of one variant, in the
  core.
- The three filters and the counts of each, in the core, in both binding
  crates and in both packages.
- `tests/reference/filters/make_reference.py` and what it stores.
- The measurement of "Speed" of `docs/specs/filters.md`.

Out, with where it goes:

- The filter of individuals, the filter by linkage disequilibrium and the
  variants of a region: later items of `docs/specs/filters.md`, not
  written.
- Whether one reader that applies several thresholds in one pass over the
  row is faster than several readers, which "How it runs" of the filters
  leaves to a measurement: a performance review after this plan, which
  has the numbers of work package 4 to start from.
- The read ahead thread of section 3 of the architecture, which reads
  the next block while a calculation works on the last one: it is built
  when popnei has its first calculation, by the plan of that
  calculation.
- Task 4.1 of `docs/plans/vars-file.md`, the bench of the vars file: it
  stays with the session that runs that plan. Work package 4 here needs a
  vars file of the big VCF and writes it with `write_vars`, not with that
  bench.
- What `plan/vars-file` commits after 2f2577e. When the two branches
  meet, a reader added there without `filtering_stats` does not compile,
  because the method has no default, so the merge cannot leave one out
  unseen.

## What has to be in place

- The branch `plan/filters` at 2c2a630 or later, in the worktree
  `.claude/worktrees/filters`. Run there on 21 September 2026: `cargo fmt
  --all --check` exit 0; `cargo clippy --workspace --all-targets -- -D
  warnings` no warning; `cargo test --workspace` `249 passed`, 2 ignored;
  `cargo wasm-check` finished; ruff `15 files already formatted` and `All
  checks passed!`; `uv run maturin develop && uv run pytest` `99 passed`;
  in `js/popnei`, `npm run build` and `npm test` `tests 62`, `fail 0`.
  These are the counts that "more than" is counted from below. The
  orchestrator also runs `bash scripts/build_pyodide_wheel.sh && node
  tests/pyodide/smoke.mjs`, which was not run when the plan was written;
  in a new worktree it needs `npm install` in `tests/pyodide` first.
- bcftools 1.24: `which bcftools` gives `/opt/homebrew/bin/bcftools`, and
  `bcftools --version` says 1.24. The three commands of "How it is
  verified" of the filters were run on `tests/reference/vcf/many.vcf` on
  21 September 2026 and keep the nine numbers of variants of the spec's
  table, the five first positions of the missing data filter at 0, and,
  chained at 0.04, 0.8 and 0.5, 215, 163 and 106 variants, the first at
  1111, 1407 and 1518.
- pyNei at ef0ca6e, which `pyproject.toml` has: `uv run python -c "from
  pynei.var_filters import filter_by_maf, gather_filtering_stats"` exits
  with 0.
- The file of the measurement, outside the repository:
  `/Users/jose/devel/popnei-bench/big.vcf`, 403572954 bytes, 100000
  variants of 1000 individuals, made by
  `crates/popnei/benches/make_big_vcf.py`. It is there. `bcftools view -H
  -i "F_MISSING<=0.1"` keeps its 100000 variants in 1.65 s, and `bcftools
  view -H` with no filter takes 1.54 s, one run each.

Every layer exists, so the five commands of "Before the work is called
done" of the `coding` skill run for every task, with `cargo wasm-check`,
and `npm run build && npm test` for a task that touches `crates/popnei-js`
or `js/popnei`.

The checks below fail today: `cargo test -p popnei --lib --
filtering_stats count_gts count_alleles --list` and `cargo test -p popnei
--lib filters:: -- --list` print `0 tests`, and none of the pytest and
TypeScript test files they name exists.

Tasks run one after another unless the plan says that two can run side by
side: every task builds the same workspace in the same worktree.

## Work package 1: what the new specs change in the code that exists

**What it gives.** A user reads, in `iter_blocks(...).pass_stats` and in
what `write_vars` returns, how many variants the pass gave, and, in
`variants.steps` and in the `repr` of a `Variants`, that it has no step
yet. No filter exists, so `filtering` is empty and `steps` is empty. It
comes first because it changes the trait that the filter of work package
3 implements, and the two places, the pass of each binding crate and the
`Variants` of each package, that work package 3 adds to.

**Deliverables.**

1. `filtering_stats` in the trait, with no default, and in every reader.
   Check: `cargo test -p popnei --lib -- filtering_stats --list` prints 5
   tests or more, with `filtering_stats` in their names: the VCF reader
   and the vars file reader give none; `reblock` and a boxed reader give
   what their source gives, over a reader written in the test that
   reports two counts; `FilteringStats` is what "The Rust interface" of
   `docs/specs/filters.md` declares, in the module `filters`, which holds
   nothing else yet. `cargo test --workspace` passes with the 249 tests
   that were there unchanged but for the method added to their readers.
2. `PassStats` in Python. Check: `uv run pytest tests/test_pass_stats.py`
   passes, where today pytest exits with 4 for a file that is not there.
   Its cases: after a whole `iter_blocks` of `many.vcf` with
   `only_passed=False`, `pass_stats` is `PassStats(num_vars=500,
   filtering={})`; read after
   three blocks of 7 variants, its `num_vars` is 21; the same two over
   the vars file that `write_vars` makes of it; `write_vars` returns a
   `VarsWritten` whose `pass_stats` is the same; `PassStats`
   and `FilteringStats` are frozen dataclasses; `iter_blocks` still
   raises for its wrong arguments at the call and not at the first block.
   The 99 tests that were there pass, changed only where they use what
   `write_vars` returns.
3. `steps` and the `repr` in Python. Check: in the same file, `steps` of
   a `Variants` just opened, from a VCF and from a vars file, is `()`,
   `Step` is the frozen dataclass of the spec, and the `repr` has the
   path of the source and says that there are no steps.
4. The same in TypeScript. Check: `npm test` in `js/popnei` gives `fail
   0` and more than 62 tests, those of `test/pass_stats.test.ts` among
   them: `passStats` of what `iterBlocks` returns after a whole pass over
   `many.vcf` is `{numVars: 500, filtering: {}}`, and after three blocks
   of 7 its `numVars` is 21; `writeVars` gives `{bytes, passStats}`;
   `steps` is an empty array. The tests of `vars.test.ts` change only
   where they take the bytes of `writeVars`. The pass is still freed when
   an iteration ends, is left or throws: `pass.test.ts` passes untouched.

**What it stands on.** The specs, and the code of `plan/vars-file` at
2f2577e.

**Tasks.**

- [x] 1.1 In the core, `crates/popnei/src/filters.rs` with
  `FilteringStats` alone, and `filtering_stats` in the trait of
  `block.rs`, in `Box`, in `Reblock`, in `VcfReader`, in `VarsReader` and
  in the readers of the tests. From "The Rust interface" of
  `docs/specs/block.md` and of `docs/specs/filters.md`. Serves
  deliverable 1.
- [x] 1.2 The Python side. In `crates/popnei-python`, the pass of
  `source.rs` counts the variants of the blocks it gives and gives that
  count and the `filtering_stats` of its chain, and `write_vars` gives
  the same for its pass. In `python/popnei`, `PassStats` and
  `FilteringStats`, `Blocks`, `VarsWritten`, `Step`, `steps` and the
  `repr`. Where the list of steps lives is what the last paragraph of
  "The Rust interface" of `docs/specs/filters.md` asks: the binding crate
  has to be able to look a kind up in it. From "In Python and in
  TypeScript" and "How it runs" of "The counts of each filter" of
  `docs/specs/filters.md`, "What a Python and a TypeScript user see" of
  `docs/specs/variant.md`, "In Python and in TypeScript" of
  `docs/specs/block.md` and "Its Python and TypeScript functions" of the
  writer in `docs/specs/io_vars.md`. Needs 1.1. Serves deliverables 2
  and 3.
- [x] 1.3 The TypeScript side, the same in `crates/popnei-js` and
  `js/popnei`, from the same parts. Needs 1.1. Serves deliverable 4. It
  can run side by side with 1.2 only in a worktree of its own; in the
  plan's worktree it runs after it.

**What could go wrong.** `Blocks` in Python is today the iterator of the
binding crate wrapped in `map`, and in TypeScript a generator whose
`finally` frees the pass. Both become objects with a property, and the
TypeScript one has to keep freeing the pass in the three cases its tests
name. `write_vars` runs its whole pass inside one call of the core, so
its count of variants comes from the core's writer or from the file it
wrote, not from a loop of the binding crate.

## Work package 2: the counts of one variant

**What it gives.** Work package 3 gets `count_gts` and `count_alleles`,
the two functions from which its three filters work out their numbers. It
stops inside the core because the spec gives them no function in Python
or in TypeScript. Their numbers are pyNei's, as literals from the spec's
table, and the comparison with pyNei that runs is the one of work package
3, through the filters.

**Deliverables.**

1. The two functions, `GtCounts` and `AlleleCounts`. Check: `cargo test
   -p popnei --lib -- count_gts count_alleles --list` prints 6 tests or
   more, with one of the two names in each: the six variants of the table
   of "How it is verified" for each function, with the spec's numbers as
   literals; the three tetraploid genotypes; a ploidy of 0; genotypes
   whose length is not a multiple of the ploidy; an allele of -2 in each
   function.
2. The two new cases of the error. Check: a cargo test of `error.rs` for
   the message of each, which names the value. `docs/specs/variant.md`
   does not say which exception of Python each becomes, which "Errors,
   and no panics" of the `coding` skill asks of a spec. By the owner's
   convention there, both are a `RuntimeError`, a defect of popnei: no
   user writes a ploidy or an allele into these functions, and a block
   with an allele of -2 or with genotypes that are not whole comes from a
   reader with a defect. The task writes that sentence into "The Rust
   interface" of the spec, and the owner can reverse it. The `errors.rs`
   of each binding crate has the arm of the two cases, read by the
   reviewer: no test of Python reaches them, because no reader gives
   such a block.

**What it stands on.** Nothing of work package 1.

**Tasks.**

- [x] 2.1 `count_gts`, `count_alleles`, their two types and the two
  cases of the error, in `crates/popnei/src/variant.rs` and `error.rs`,
  with the arm of each binding crate's `errors.rs`. From "The counts of
  one variant" and "The Rust interface" of `docs/specs/variant.md`; they
  mirror `_calc_gt_is_missing`, `_calc_gt_is_het` and
  `_count_alleles_per_var` of `pynei/gt_counts.py`. Serves deliverables 1
  and 2. A wrong count is silent: this task has its own commit, and the
  literals of deliverable 1 and the comparisons of work package 3 guard
  it.

## Work package 3: the three filters and the counts of a pass

Its review is in two parts, a change of 21 September 2026: work package
2 and tasks 3.2 and 3.3, the counts and the filter of the core that is
built on them, which are one piece of code, are reviewed together when
3.3 is done; tasks 3.4 and 3.5, the two bindings, with 3.1, when 3.5 is
done.

**What it gives.** A user calls `variants.filter_by_missing_data(0.04)`,
`filter_by_maf(0.8)` and `filter_by_obs_het(0.5)`, in Python, and the
three methods of TypeScript, and every consumer after that,
`iter_blocks` and `write_vars`, gets the variants that pass and says in
`pass_stats.filtering` what each filter was given and kept.

**Deliverables.**

1. The reference. Check: `uv run --no-project python
   tests/reference/filters/make_reference.py` runs the three commands of
   "How it is verified" at the nine thresholds, and the chain of the
   three at 0.04, 0.8 and 0.5, on `many.vcf`, refuses a bcftools that is
   not 1.24, and writes the positions each one keeps into a file beside
   it; run again, `git status --short tests/reference/filters` shows
   nothing; the numbers of lines of what it stores are the numbers of the
   spec's table and 215, 163 and 106.
2. The filter of one block. Check: `cargo test -p popnei --lib
   filters::tests -- --list` lists the worked example of "How it is
   verified" as the first test of the file, made at
   `VarFilter::filter_block` with every threshold and every set of kept
   variants of the paragraph under its table, and tests for: the counts
   6 and 4, 4 and 3, 3 and 1 of the three filters chained on it; a
   threshold of -0.1, of 1.5 and of NaN in `VarFilter::new`; a block that
   fails `check` and a block with variants and no genotypes, each left as
   it was with nothing added to the counts; a block of no variants; the
   two variants with every genotype `./.` of "What pyNei does that is
   odd"; and 29 missing genotypes of 100 individuals kept at 0.29.
3. The reader. Check: cargo tests in the same module, made at
   `next_block` of a `FilteredReader` over a `VcfReader` on `many.vcf`
   with every variant given, assert the number kept in each of the nine
   rows of the table and the five positions where the table has them,
   with blocks of 7 and of the default size; the chain at 0.04, 0.8 and
   0.5 gives the three pairs of "How it is verified" of the counts, the
   observed heterozygosity one first, and the 106 variants with their
   three positions; a second filter of a kind the chain has is an error
   of `FilteredReader::new`, also with another filter between; the
   blocks hold the genotypes when the consumer asked for the positions
   alone; and the three rules of a reader of "How it runs", each over a
   reader written in the test: a block left with no variant is not given,
   after an error of the source the source is not called again, and a
   source that gives a block of no variants is the error of `reblock`.
   `cargo test -p popnei --lib filters:: -- --list` prints 20 tests or
   more.
4. Python. Check: `uv run pytest tests/test_filters.py` passes, where
   today the file is not there. It has the comparison with pyNei of "How
   it is verified" of the filters, at the nine rows of the table, with a
   second pass; the `ValueError` of -0.1, 1.5 and NaN and the `TypeError`
   of no threshold, in each of the three methods; the comparison of
   `pass_stats.filtering` with pyNei's `gather_filtering_stats` and every
   test of the steps of "How it is verified" of the counts; and the test
   of `write_vars` that "Its Python and TypeScript functions" of the
   writer in `docs/specs/io_vars.md` asks for. The `repr` of a `Variants`
   with the three filters names each kind and its threshold.
5. TypeScript. Check: `npm test` gives `fail 0` and
   `test/filters.test.ts` runs the cases of the last paragraph of each
   "How it is verified" of `docs/specs/filters.md`: the first row of each
   filter of the table with its five positions, the `Error` of 1.5, the
   three pairs of the chain and the `Error` of a second maf filter.
6. The wheel of pyodide still builds and its smoke test exits with 0.

**What it stands on.** Work packages 1 and 2.

**Tasks.**

- [x] 3.1 The reference script and what it stores, in
  `tests/reference/filters/`, written as
  `tests/reference/vcf/make_reference.py` is. From "How it is verified"
  of the filters and of the counts. Serves deliverable 1. It touches no
  other file and can run side by side with work package 1 or 2 in a
  worktree of its own.
- [x] 3.2 `VarFilteringCriterion`, `VarFilter` and the two cases of the
  error, in `crates/popnei/src/filters.rs` and `error.rs`. From "What
  they give", "Half called genotypes, variants with nothing called, and
  what pyNei asserts", "What pyNei does that is odd" and "The Rust
  interface". It mirrors `_filter_chunk_by_missing`,
  `_filter_chunk_by_maf` and `_filter_chunk_by_obs_het` of
  `pynei/var_filters.py`. Needs 2.1. Serves deliverable 2. Which variants
  a filter keeps is silent when wrong: its own commit, guarded by
  deliverables 2, 3 and 4.
- [x] 3.3 `FilteredReader`, in the same file, and its cargo tests, with
  the positions of 3.1 as their literals. From "How it runs" of the
  filters and of the counts, and "The Rust interface". Needs 1.1, 3.1
  and 3.2. Serves deliverable 3.
- [ ] 3.4 The Python side: in `crates/popnei-python`, every pass, of
  `iter_blocks` and of `write_vars`, builds its chain of
  `FilteredReader` from the steps the `Variants` has when it starts, and
  the two errors reach a user as a `ValueError` at the call of the
  method; in `python/popnei/variant.py`, the three methods. From "In
  Python and in TypeScript" of both items and the last paragraph of "The
  Rust interface". Needs 1.2 and 3.3. Serves deliverables 4 and 6.
- [ ] 3.5 The TypeScript side, the same in `crates/popnei-js` and
  `js/popnei`. Needs 1.3 and 3.3. Serves deliverable 5. After 3.4 in the
  plan's worktree.

**What could go wrong.** pyNei's `gather_filtering_stats` gives the last
filter first, and popnei's dict has the order of the steps, so the pytest
test compares kind by kind and asserts popnei's order apart. A step added
inside the loop of an `iter_blocks` must not reach the pass that runs:
the chain is built once, when the pass starts. The counts that
`iter_blocks` shows mid pass can be of more variants than the user got,
as "A pass that was not finished" says, so a test that reads them mid
pass asserts `num_vars` and not the counts of the filters. The bcftools
command of the observed heterozygosity filter compares the heterozygous
genotypes with the threshold times the called genotypes, where popnei
divides the first by the second and compares with the threshold, and
the two can differ for a variant exactly on the threshold. The spec says
they do not on `many.vcf` at these thresholds; if a test of deliverable
3 of that filter fails by one variant, this is where to look first.

## Work package 4: the measurement of "Speed"

**What it gives.** The owner knows what the missing data filter costs a
pass, in popnei, in pyNei and in bcftools, before any work on its speed.
It stops inside the core for popnei's own time, and pyNei is timed from
Python.

**Deliverables.**

1. The bench. Check: `cargo bench --bench filter_vars -- <path> --threads
   18 --runs 5 --max-missing-rate 0.1` prints the wall time of each run
   of a whole pass with the genotypes alone asked for, over a VCF or a
   vars file by the name of the path, and without `--max-missing-rate`
   the same pass with no filter; it has no harness, as `read_vcf` has
   none, and its header says how its two files are made.
2. The measurement, in `docs/reports/filters-measurement.md`, as "What
   measurement there is" of the `performance-review` skill asks. Check: it
   has, for `big.vcf` and for the vars file that `write_vars` makes of it
   with the default size, the median of 5 runs of the pass with the
   filter at 0.1 and with no filter, and their difference, on 1 thread
   and on 18; the same two passes and their difference for pyNei on
   `big.vcf`, and for `bcftools view -H`, with and without `-i
   "F_MISSING<=0.1"`; the variants kept, which are the same in the three;
   the machine, the build and the load average of each set of runs; and
   the same figures at a threshold of 0.03.
3. "Speed" of `docs/specs/filters.md` has the numbers, with what they
   were measured on, in the place of "the measurement comes first", and
   still sets no number to reach: that is the owner's.

**What it stands on.** Work package 3, and `big.vcf`. The orchestrator
runs it with nothing else building on the machine.

**Tasks.**

- [ ] 4.1 The bench, the vars file of `big.vcf`, the three sets of
  timings and the report, and the paragraph of the spec. It changes no
  code of the library. From "Speed" of `docs/specs/filters.md`. Serves
  deliverables 1, 2 and 3.

**What could go wrong.** At 0.1 the filter keeps every one of the 100000
variants of `big.vcf`, whose genotypes are missing at a rate of 0.03, so
that threshold times the counting and no compaction of a block. At 0.03
bcftools keeps 54773, and that is why the report has both. The cost of a
filter is a difference of two times that are close, 1.65 s and 1.54 s in
bcftools, so a load on the machine that changes between the two sets of
runs is larger than what is measured: each pair is run back to back, and
a pair whose load average differs is run again.

## How the whole plan is checked

From a clean clone of the branch: the five commands of the `coding`
skill, `cargo wasm-check`, `npm run build` and `npm test` in `js/popnei`,
the build of the wheel of pyodide and its smoke test, and the work
report, `docs/reports/filters.md`, with the counts of "What has to be in
place" beside the ones at the end and the numbers of work package 4.
