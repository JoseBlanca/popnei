# Work report: the filter of the first n variants and the filter that keeps variants at random

5 October 2026. The report of `docs/plans/filters-random-first-n.md`, which
builds issues 7 and 6 of the repository from `docs/specs/filters.md`, on the
branch `filters-random-first-n` from `main` at eae29a2. State: under way.

## Before the first task

The checks of the plan's "What has to be in place" were run on 5 October
2026 and gave what the plan says. The board had no message of another
branch about the steps or the counts of a pass; `spec/writer-regions-density`
still has a worktree although the board says it was merged, which is not of
this plan.

## 1. The filter of the first n variants

Task 1.1, the core, is c5a96e0: `FirstNReader`, `PassStep::FirstN`,
`stopped_early`, `refuse_a_step_after_the_first_n` and three cases of the
error, in `crates/popnei/src/filters/first_n.rs`. No existing reader needed a
fix: with the chain on the thread of the reader one block ahead, the filter's
`None` ends that thread with the source asked for 2 blocks, and no test hung.
The VCF of 100000 variants of 10 individuals, 6989013 bytes, was read for
8192 bytes by a pass with the first 100, and the vars file for its batch 0
alone. Checked on 5 October 2026: the six cargo commands pass, `cargo test
--workspace` with 1389 passed and `cargo test -p popnei --no-default-features`
with 1239, 16 more than before in each, the 16 of `cargo test -p popnei --lib
first_n`. One test asserts counts the spec did not give, the MAF filter at 0.8
before the first 10 in blocks of 7, 14 given and 12 kept, which bcftools
confirmed and which went into the spec in e2f253d.

Three of the 17 tests of 1.1 read a source that never ends and kept every
block. Written first, against the stub of the filter, they grew without bound,
and cargo ran them at once until the owner's machine, 64 GB, ran short of
memory and the owner had to quit other programs. The source now panics past
100 blocks, with a test that a pass nothing ends fails there (a159231), and
the coding skill asks a test written first to fail without hanging or
growing (50e4def).

Task 1.2, Python, is 07ae287: `filter_first_n`, `PassStats.stopped_early`,
every method that adds a step through one `Steps::add`, and the counts of a
pass built by one function, `pass_counts_of`, which the 13 modules of the
crate that return counts call. Task 1.3, TypeScript, is d5b869d: the same
with `filterFirstN`, `Steps::add` and `PassCounts::of`, called by the 15
consumers. The two ran side by side. Both refuse a `num_vars` of 0 with the
core's case, whose message does not name the argument, since the core refuses
it in `FirstNReader::new` and no reader exists when the step is added.

The deliverables, checked on 5 October 2026 at 07ae287: the six cargo
commands pass, `cargo test --workspace` with 1390 passed and `cargo test -p
popnei --no-default-features` with 1240; `cargo test -p popnei --lib first_n
-- --list` lists 17; `uv run pytest` gives 692 passed, 23 of them of `-k
first_n`; `npm test` gives tests 518, pass 517, fail 1, the test of
`test/gwas.test.ts` that fails on `main`; the tests named `filterFirstN` pass;
`npm run test:browser` gives 8 passed.

### The review of work package 1

Five reviewers, `spec`, `tests`, `errors`, `binding` and `architecture`,
over a19b3e6..61b0841; `api` and `numbers` are left to the review of the
whole branch, which has all seven. No finding was of a wrong number. What
held, each fixed in a commit of its own, test first, with the spec amended
before the code in 61ffe48, be4877e and 2d4cc16:

- The two binding crates checked the two refusals of a step in opposite
  orders, so MAF, the first n and MAF again gave Python's message in one
  language and the second filter's in the other, and Python's told the user
  to move a filter that was already first. Found by three reviewers. One
  function of the core, `refuse_a_step`, now checks a second filter of its
  kind first, and `chain_of` and both bindings call it (f4f6b99).
- Both bindings refused a `num_vars` of 0 themselves, with a message that
  did not name the argument. The core builds the step with
  `first_n_step`, whose message names `num_vars`, `numVars` in TypeScript
  (59e1bc1). Found by three reviewers.
- The refusal of a second filter of the first n named neither n. It has a
  case of its own, `FirstNThatIsSet`, which names both (a352c6c); work
  package 2 needs its own case for a second random filter.
- Untested: a block of no variants from the source, `set_needs` passed on,
  the offer of regions refused (b6e985a); `filterFirstN(2 ** 53 - 1)`
  accepted (fe4e7ee). Each core test failed with the code broken on purpose.
- The progress bar of a page stops short of the end when the filter ends a
  pass, which no spec said. `docs/specs/js_sources.md` says it, and a test
  over a VCF of 10388977 bytes asserts that the last call is at 398976
  (3145a9f).
- The opening of `docs/specs/filters.md` still said the filter had no code.

Not taken:

- `stopped_early` read in the middle of an `iter_blocks` can be true before
  the user has the n variants, which `reblock` holds back: the spec already
  says the counts read while a pass runs can be ahead of the blocks given.
- In a pyodide build, whose `usize` is 32 bits, a `num_vars` above
  4294967295 would be refused: no first sample is that large.
- The Python test of every consumer lists them by hand, where TypeScript
  checks its list against the crate's: the Python crate keeps no list of its
  consumers to check against.

After the fixes, at 2d4cc16 on 5 October 2026: the six cargo commands pass,
`cargo test --workspace` with 1395 passed and `cargo test -p popnei
--no-default-features` with 1245; `uv run pytest` 693 passed; `npm test`
tests 521, pass 520, fail 1, the test of `test/gwas.test.ts` that fails on
`main`; `npm run test:browser` 8 passed. Work package 1 is done.

## 2. The filter that keeps variants at random

Task 2.1, the core, is c5c842d, SplitMix64 alone with Java's five draws as
its test, and ddbf458, the filter in `crates/popnei/src/filters/random.rs`:
`DEFAULT_RANDOM_FILTER_SEED`, `RandomFilter`, `RandomlyFilteredReader`,
`PassStep::Random`, and two cases of the error, a keep rate out of range and
a second random filter, which names the keep rate and the seed set and asked
for. A block the filter refuses draws no number: it draws on a copy of the
generator and keeps the copy once the block is compacted. Checked on 5
October 2026: `cargo test --workspace` 1415 passed, 20 more, the 20 of
`cargo test -p popnei --lib random_filter`; clippy, fmt and both wasm
checks pass.
