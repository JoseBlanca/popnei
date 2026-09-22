# Work report: the stats module and the filter of individuals

The plan `docs/plans/stats.md` is under way, on the branch `plan/stats`,
in the worktree `.claude/worktrees/stats`, since 22 September 2026. The
orchestrator, in this report, is the session of the assistant that runs
the plan: it sends each task to a subagent on Opus, checks what comes
back and has each work package reviewed.

## Before the first task

The owner approved the plan in chat on 22 September 2026. The branch
stands on `spec/stats` at 7bea42b, which holds the two specs, the plan
and `tests/reference/stats/`, and which is not in `main`; `spec/stats`
is from `main` at 7d8366f.

Run by the orchestrator in the worktree at 7bea42b: `cargo fmt --all
--check` exit 0; `cargo clippy --workspace --all-targets -- -D warnings`
no warning; `cargo test --workspace` `306 passed`, 2 ignored; `cargo
wasm-check` finished; ruff `18 files already formatted` and `All checks
passed!`; `uv run maturin develop && uv run pytest` `174 passed`; `npm
install`, `npm run build` and `npm test` in `js/popnei` `tests 126`,
`fail 0`; after `npm install` in `tests/pyodide`, `bash
scripts/build_pyodide_wheel.sh` built the wheel and `node
tests/pyodide/smoke.mjs` exited with 0. `which plink2` gives
`/opt/homebrew/bin/plink2`, `PLINK v2.0.0-a.7.7 M1 (18 Sep 2026)`, and
`which bcftools` `/opt/homebrew/bin/bcftools`, version 1.24. pyNei's
`calc_per_var_distribs`, `calc_per_sample_stats`, `filter_samples` and
`load_vars` import under `uv run`, and
`/Users/jose/devel/pynei/test/gwas_reference/sim_missing.vars` is there.
`uv run python tests/reference/stats/make_reference.py` printed `done`
and exited with 0; after it `git status --short tests/reference/stats`
showed `panel.vcf.gz` modified, and only in the two bytes of the
modification time of its gzip header, bytes 5 and 6, the text of the
file being the same: the committed file was put back with `git
checkout`. The bcftools commands of "How it is verified" of the filter
of individuals on `many.vcf` give `ind05 ind00 ind49`, 423 variants with
the positions 1000, 1037, 1074, 1111 and 1148 first, and 26 with the
filter over the 50 individuals. `/Users/jose/devel/popnei-bench/big.vars`
is there, 81356714 bytes, and `big.vcf` beside it, 403572954 bytes.

The start message of the branch is on the board,
`.claude/board/2026-09-22T1020-plan-stats.md`.

## Work package 1: the steps of a pass as one enum

Task 1.1, commit c902c09: `PassStep` in `crates/popnei/src/filters.rs`,
`#[non_exhaustive]` as the spec has it, with `VarFilter` and
`KeepIndividuals`, and `chain_of` and `refuse_a_second_filter_of_a_kind`
over a list of it; the `Step` of each binding crate is a struct that
holds one `PassStep` and the names and values of its arguments as the
user sees them, and the two crates' own loops over criteria are gone.
The error enum of the core gets one case, `PassStepNotBuilt`, a
`RuntimeError` in Python, which both functions give for a
`KeepIndividuals` step until task 2.2 builds it; the subagent chose it
over `Ok(())` in `refuse_a_second_filter_of_a_kind` so that a second
filter of individuals could not pass in silence between the two tasks.
Task 2.2 removes the case, its arm in the Python crate and its test.

Both deliverables are met. `cargo test -p popnei --lib -- PassStep
--list` prints `3 tests`: the kind of each variant, the chain of the
three threshold filters over the worked example of the spec at 0.4, 0.88
and 0.25 keeping variant 5 with the counts 6 and 4, 4 and 3, 3 and 1,
and the error of a step not built. `cargo test --workspace` `309
passed`, 2 ignored, the 306 that were there and the three. `grep -rn
"VarFilteringCriterion" crates/popnei-python/src crates/popnei-js/src`
gives five lines per crate: the import, the three criteria built from
the user's argument, and the signature of the function that adds a
threshold filter, which takes the criterion just built. `uv run pytest`
`174 passed` and `npm test` `tests 126`, `fail 0`, both untouched.
`cargo fmt`, `cargo clippy`, `cargo wasm-check` and ruff clean.

Nothing was changed in the plan. The review of this work package is
made together with work package 2, as the `following-plans` skill allows
when two are one piece of code: it changes nothing a user sees, and task
2.2 rewrites the two functions it changed.

Two things for the next tasks. The filter of `cargo test` is a case
sensitive substring of the test path, so the tests the plan's checks
list by a type name, `PassStep`, `IndividualsReader`, sit in a module of
that name inside `mod tests`, with an `#[expect(non_snake_case)]`. And
since `PassStep` is `#[non_exhaustive]`, a `match` on it outside the
core needs a wildcard arm, so a step added to the core no longer stops a
binding crate from building: neither crate matches on it today, and the
tasks of the Python and TypeScript sides of a new step have to add it by
hand.

The subagent of task 1.1 used 179631 tokens in 55 tool calls.

## Work package 2: the filter of individuals

Tasks 2.1 and 2.2 went to one subagent, in that order, with a commit
each, since the reader of 2.2 is written over the method of 2.1.

Task 2.1, commit fe237ec: `Block::retain_individuals` in
`crates/popnei/src/block.rs`, the gather of each row on rayon off wasm
through a buffer per thread and the pack on one thread, as "How it
runs" of the filter has it. Before it, commit e8efa81 adds a paragraph
to `docs/specs/filters.md`, a choice the spec left to the code: which
error each of its three refusals is. They are three new cases of the
error of the crate, defects of popnei and a `RuntimeError` in Python,
as the `keep` of `retain_vars` with the wrong number of values is,
because `resolve_individuals` refuses the name behind each of them
before a call of a user reaches them.

Task 2.2, commit ef17cc5: `resolve_individuals`, `IndividualsReader`,
the `KeepIndividuals` arm of `chain_of` and of
`refuse_a_second_filter_of_a_kind`, and the four cases of the spec in
the error, each a `ValueError` in Python; the temporary case of work
package 1 is gone with its arm and its test. The JavaScript crate
needed no arm, since it turns every error of the core into one `Error`.

Task 2.3, commit 91b74fc: `filter_individuals` in
`crates/popnei-python/src/steps.rs`, which resolves the names against
the individuals of the source at the call and refuses a second filter
of the kind there; `Variants.filter_individuals` with its docstring,
`individuals` and `num_individuals` read through the steps, and the
`repr`; and `tests/test_filter_individuals.py`, 13 tests. Three things
the spec does not say, decided by the subagent: the `Steps` of the
binding crate takes the individuals of the source when it is built, so
the refusal cannot be given another list; `filter_individuals("ind05")`,
a bare string, is a `TypeError` that says to write `("ind05",)`, as
`iter_blocks(fields="pos")` already is, where otherwise the user would
read that `i` is not an individual; and the `repr` prints every kept
name, `individuals(individuals=('ind05', 'ind00', 'ind49'))`, so a
filter of 500 individuals prints 500 names. The last two are the
owner's to reverse.

Task 2.4, commit 4c5e5ae: the same in `crates/popnei-js/src/steps.rs`
and `js/popnei/src/variant.ts`, `filterIndividuals`, with
`js/popnei/test/filter_individuals.test.ts`, 13 tests. The arguments of
a step used to cross to JavaScript as one number each; now a threshold
crosses in one flat array and the names of the kept individuals in
another, with a count of names per argument that says which is which,
0 names meaning a threshold, which is unambiguous because a filter of
no individual is refused. The `Default` of the `Steps` of both crates
is gone, since `new` takes the individuals.

The five deliverables are met, run by the orchestrator at 4c5e5ae.
`cargo test -p popnei --lib -- retain_individuals --list` prints `8
tests`, where the plan asks 4 or more; `cargo test -p popnei --lib --
resolve_individuals IndividualsReader --list` `15 tests`, where it asks
8 or more; `uv run pytest tests/test_filter_individuals.py` `13 passed`
and `uv run pytest` `187 passed`, from 174; `npm test` `tests 139`, `fail
0`, from 126; `bash scripts/build_pyodide_wheel.sh` built the wheel and
`node tests/pyodide/smoke.mjs` exited with 0. `cargo test --workspace`
`331 passed`, 2 ignored, from 309. `cargo fmt`, `cargo clippy`, `cargo
wasm-check` and ruff clean.

Nothing was changed in the plan.

The subagents used: tasks 2.1 and 2.2 together 203867 tokens in 84 tool
calls; task 2.3 147130 in 63; task 2.4 182877 in 60. None had to be
sent back.

### The review of work packages 1 and 2

Seven reviewers, one per category, at 4c5e5ae. None found a wrong
number: the spec reviewer reproduced every value of the item against
bcftools 1.24 and pyNei, and the gather of 18 random subsets of
`many.vcf`, of haploid and triploid files and of the 50 individuals
reversed against direct indexing; the numbers reviewer compared
`retain_individuals` with a naive gather over 625 cases of sizes,
ploidies and kept sets, all equal. What they found, fixed in nine
commits from 15281d2 to 2d360b6 by one subagent, was in what happens
when a bound is broken, in the tests and in the texts:

- The compaction of a block skipped, in silence, a genotype it could not
  reach and a row it could not pack, when the arrays of the block were
  shorter than it states: unreachable today, since `check` runs first,
  but a wrong genotype with no error the day a caller reaches it. Three
  reviewers found it from three sides. Now the two helpers give the
  error of `check`, and a test pops one allele and asserts the error
  and the block as it was. The reader of the filter dropped in silence
  an index beyond the individuals of its source when it built its
  names, unreachable for the same reason; it gives the error now.
- The buffer of the gather was allocated once per rayon job, about one
  per three rows, and not once per thread as its comment said: 1861 to
  2245 jobs for a block of 5000 variants of 1000 individuals, 500 kept,
  on 18 threads, measured by the fixer on the M5 Pro; with a floor of 64
  rows per job, 63. What the gather costs with either was not measured.
- Two tests could not fail: the one thread against several test had 300
  identical rows, since its 20 alleles per row cycled 5 values, and a
  mutation that copies row 0 into every row passed it; and no test
  covered the reader's own error, so `finished` could be dropped from
  that branch with every test green. Both have their test now.
- The kind of an argument of a step crossed to TypeScript as a count of
  names, 0 meaning a threshold, so a later argument of neither kind
  would have taken the next threshold in silence; two reviewers found
  it. A kind per argument crosses now, and an unknown one is an
  `Error`.
- The rule that gives the individuals of the next pass was written twice,
  by two mechanisms, in Python and in TypeScript. It is one function of
  the core now, `individuals_of`, added to "The Rust interface" of the
  spec in 2d360b6, and each binding crate answers with it.
- In Python, a non-string name and a non-iterable argument gave the
  `TypeError` of Python or of pyo3, naming no argument; both name it now,
  as the TypeScript does. The message of an unknown name points at where
  the names are. A docstring that counted three filter methods, and a
  cut sentence in the comment of the exception table, are mended.

Not taken, with the reason: that the `reblock` at the end of
`iter_blocks` sizes the blocks after the filter for the kept individuals
and not the source's, above 500 individuals, because `docs/specs/block.md`
says its default is for the individuals of the reader it is given, and
the sentence of the filter spec is about the filter's own reader; the
spec of the filter says so now, in 15281d2, which also says that the 500,
423 and 26 of `many.vcf` come out with `only_passed` false. That a chain
built over a reader that already holds a filter of individuals is not
refused, because the spec limits that refusal to the threshold filters,
whose `FilteredReader::new` asks the reader. That `filter_individuals`
accepts a generator under a `Sequence[str]` hint, which harms nobody.

Seen outside the scope, for the owner: `docs/glossary.md` line 154 ends
in the middle of a sentence, from before this branch and in a shared
file, so it is not touched here; the Python package ships no `py.typed`,
so its type hints never reach a user's type checker; and the `Steps` of
the JavaScript crate clones the names of the source at every pass.

After the fixes, at 2d360b6: `cargo test --workspace` `334 passed`, 2
ignored; `retain_individuals --list` `9 tests`; `resolve_individuals
IndividualsReader --list` `16 tests`; `uv run pytest` `189 passed`, of
`tests/test_filter_individuals.py` 15; `npm test` `tests 140`, `fail 0`;
the wheel of pyodide built and its smoke test exited with 0; `cargo
fmt`, `cargo clippy`, `cargo wasm-check` and ruff clean.

The reviewers used, in tokens: spec 148820, tests 151940, numbers
134789, errors 127376, api 113582, architecture 115202, binding
126337. The fixer used 231504 in 129 tool calls.
