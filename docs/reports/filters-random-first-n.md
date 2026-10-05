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

Task 2.2, Python, is 63941c6: `filter_randomly`, the keep rate refused by
building the core's `RandomFilter`, the seed read by a function of its own.
`do_pca_from_variants` with 10 components and `calc_gwas` with a kinship and
GRAMMAR-Gamma over the random filter at 0.1 and a seed of 42 give, exactly,
with no tolerance, the results of the same calls over a vars file that holds
the 45 variants alone, so the two passes of each saw the same variants. Task
2.3, TypeScript, is 9bc1163: `filterRandomly(keepRate, {seed})`, the default
seed from the core's constant. They ran side by side.

The deliverables, checked on 5 October 2026 at 63941c6: the six cargo
commands pass, `cargo test --workspace` 1415 passed and `cargo test -p popnei
--no-default-features` 1265; `cargo test -p popnei --lib random_filter`
lists 20; `uv run pytest` 719 passed, 26 of `-k filter_randomly`; `npm test`
tests 541, pass 540, fail 1, the test of `test/gwas.test.ts` that fails on
`main`; `npm run test:browser` 8 passed. The review of this work package is
the review of the whole branch below, which has all seven categories and
would have repeated it.

## The review of the whole branch

Seven reviewers, one for each category of the `code-review` skill, `spec`,
`tests`, `numbers`, `errors`, `api`, `architecture` and `binding`, over
eae29a2..bcc5b94, on 5 October 2026. It is the review of work package 2 as
well. No reviewer found a wrong number. The generator was checked against
Java to the bit over 1000 numbers from each of four seeds, the table of the
spec twice, and the draws with blocks of 1, 7, 100 and the default size, on
1 and 8 threads, kept the same variants. What held, the spec amended first
in dba6672 and each fix of the code in a commit of its own:

- `np.True_` and `np.False_` were taken as a keep rate of 1 and 0, and as a
  threshold of 1 and 0 by the threshold filters before this branch: a
  boolean from a numpy array filtered nothing, or everything, and said
  nothing. (errors)
- A TypeScript keep rate of 2 read back as 2.0 in its message. (errors)
- Python read `num_vars` into a 32-bit word under pyodide, where TypeScript
  takes up to 2^53 - 1 in the same browser. (binding; left at work package
  1 and taken now that the fix was cheap and removed dead code)
- The TypeScript binding checked a whole number three times over three
  constants of one value, and a docstring disagreed with its message.
  (binding, api)
- The TypeScript default seed was an `Option` that the binding turned into
  42, which the coding skill forbids, and its doc wrote 42 by hand. (api)
- The keep rate crossed to TypeScript as a threshold. (api)
- `refuse_a_second_filter_of_a_kind` and `refuse_a_step_after_the_first_n`
  stayed public although no binding may call them apart. (api)
- A field of the error of a second random filter was an unnamed pair; the
  doc of `filtering_stats` listed three of its eight kinds. (api)
- Read while an `iter_blocks` runs, `stopped_early` can be true before the
  user has the n variants, which `docs/specs/variant.md` said could not
  happen; the spec now says it of a finished pass. (architecture)
- Section 1 of `docs/architecture.md` said every filter decides its rows on
  the pool; the filter by linkage disequilibrium and the random filter do
  it in order. (architecture)
- Untested: Java's numbers were checked to six decimals and not to the bit,
  and nothing told `<` from `<=` against the keep rate, so two changes that
  would move a sample passed every test (numbers); no core test passed a
  seed other than 42 through `chain_of`, none ran the draws on more than one
  thread (tests); the filter by regions and the filter of individuals before
  the random filter, and a keep rate of 0 over a pass, had no test (spec,
  tests).

Not taken:

- Tests that the order of the steps, or `only_passed`, changes the sample:
  both follow from drawing in order, which the owner chose, and the tests
  would assert two different numbers.
- Tests of the checks inside the TypeScript binding crate that turn a
  non-whole number from the package into a defect: only a defect of the
  package reaches them, and a test would have to call the crate past the
  package.
- The reviewer's example of the bits of the first number,
  0x3fec4415072f63b9, is not what Java prints, 0x3fe7bae644c5fd6d, which
  the spec and the tests use.
