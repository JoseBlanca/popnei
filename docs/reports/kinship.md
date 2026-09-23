# The work on the kinship

23 September 2026. What was done while `docs/plans/kinship.md` was carried
out, written as the work went. The plan builds the genomic relationship
matrix of a set of variants, which says for every pair of individuals how
much of their genome they share beyond what two individuals drawn at random
from the panel share, from `docs/specs/kinship.md`. The work is on the
branch `plan/kinship` in the worktree `.claude/worktrees/kinship`, and
nothing of it is on `main`.

State: under way, work package 1 of 3.

## The starting commit

`85d27a2`, the local `main` at `75c5fd3` with `spec/gwas` merged into it.
The merge had no conflict. What the plan's "What has to be in place" asks
for was there: `wasm-bindgen` at `/Users/jose/.cargo/bin/wasm-bindgen`; the
reference data under `tests/reference/kinship/`, both panels and both
plink2 matrices with their `.id` files; `tests/reference/dists/panel.vcf.gz`
and `tests/reference/gwas/phenotypes.csv`; and the `wasm-check` alias in
`.cargo/config.toml`.

The eight checks the plan lists were run on that commit and each gave the
number the plan gives, so the work started from the state the plan was
written against.

| command | result |
| --- | --- |
| `cargo fmt --all --check` | clean |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo test --workspace` | 604 passed, 2 ignored in the core crate; 149 in the linear algebra crate |
| `cargo test -p popnei-linalg --no-default-features` | 136 passed |
| `cargo test -p popnei --lib kinship -- --list` | `0 tests, 0 benchmarks` |
| `cargo test -p popnei --lib pca -- --list` | `47 tests, 0 benchmarks` |
| `uv run maturin develop && uv run pytest` | 347 passed |
| `npm run build && npm test` in `js/popnei` | 242 pass, 0 fail |

## Work package 1: one row pass for the PCA and the kinship

Under way. One change to the plan was made before the task started.

### What was changed in the plan, and why

The check of deliverable 1 was `grep -c "fn the_standardized_row"
crates/popnei/src/*.rs` giving 1 in `variant.rs` and 0 in `pca.rs`, "where
today it is 0 and 1". Run on the starting commit it gives 5 in `pca.rs`,
because without the open parenthesis the pattern also counts
`the_standardized_rows`, which drives a whole block and has one definition
for the threads and one for WebAssembly, and
`the_standardized_rows_one_by_one`. It counts the wrapper the
`bench-internals` feature exposes as well. A check whose starting number is
wrong cannot say whether the move happened, so it was replaced by `grep -rn
"fn the_standardized_row(" crates/popnei/src/*.rs`, which names `pca.rs`
twice today, at the function and at that wrapper, and must name only
`variant.rs` when the task is done. The new check is the stronger of the
two: it distinguishes the row from the block driver, which the old one
could not.

### What was built

Task 1.1, in commit `23e82e9`: six files, 768 lines added and 472 removed.
The pass that turns one variant into its standardized dosages,
`the_standardized_row` with its buffers `RowScratch` and the helpers it
calls, left `crates/popnei/src/pca.rs` for `crates/popnei/src/variant.rs`,
and now takes a `DosageOptions` holding `transform_to_biallelic` and a
`DosageScale`, which is either the standard deviation of the dosages, what
the principal components of the variants divide by, or
`sqrt(ploidy * p * (1 - p))`, what the kinship will divide by in work
package 2. The three functions that drive a whole block of variants stayed
in `pca.rs`, one for the threads, one for WebAssembly, which has none, and
one that reads the rows one after another.

`Error::PcaVariantWithMoreThanTwoAlleles` became
`Error::VariantWithMoreThanTwoAlleles`, since the kinship raises it too.
Its message is unchanged character for character.

### The deliverables

Each check was run by the orchestrator on `eae2e59`, the commit of the code
with the correction to the spec on top of it.

| deliverable | command | result |
| --- | --- | --- |
| 1, the pass lives in `variant.rs` | `grep -rn "fn the_standardized_row(" crates/popnei/src/*.rs` | two hits, both `variant.rs`; `pca.rs` named in neither |
| 2, the PCA computes what it computed | `cargo test -p popnei --lib pca -- --list` | `47 tests, 0 benchmarks` |
| 2 | `cargo test --workspace` | 607 passed, 0 failed, 2 ignored in the core crate; 149 in the linear algebra crate |
| 2 | `uv run pytest tests/test_pca.py` | 33 passed |
| 3, one error for both callers | `cargo test -p popnei --lib more_than_two_alleles -- --list` | `2 tests`, one named for the PCA and one for the row pass called with either divisor |

The other checks of the `coding` skill on the same commit: `cargo fmt --all
--check` and `cargo clippy --workspace --all-targets -- -D warnings` clean,
`cargo test -p popnei-linalg --no-default-features` 136 passed, `cargo
wasm-check` clean, `uv run pytest` 347 passed, and `npm run build && npm
test` in `js/popnei` 242 pass with 0 fail.

Deliverable 2 asks for 604 tests in the core crate and there are 607. The
three are new in `variant.rs`, which deliverable 3 asks for; none was
removed. Every line of `pca.rs` holding an assertion is the same before and
after the move, checked by sorting them and comparing the two lists, except
one `const _: () = assert!(GENOTYPES_PER_RUN <= 255)`, which is read when
the code is compiled and not when it is run, and which moved to `variant.rs`
with the constant it is about. One test was renamed, from
`a_variant_of_three_alleles_is_the_error_that_names_its_position` to
`a_variant_of_more_than_two_alleles_is_the_error_that_names_its_position`,
so that the check of deliverable 3 finds a test for each of the two callers.
The 45 test names the two versions of the file hold differ in that one name
and nothing else.

### What the review found

Seven reviewers read the commit, one per category, each with a fresh
context. Sixteen findings held and were fixed in seven commits, `d6b00af`
to `e445cc5`. The checks were run again after them: the core crate passes
610 tests with 2 ignored, the linear algebra crate 149 and 136 without its
default features, Python 347, JavaScript 242, and `cargo fmt`, `cargo
clippy --workspace --all-targets`, `cargo clippy -p popnei --all-targets
--features bench-internals` and `cargo wasm-check` are clean. The three
deliverables were checked again and still hold, with the 47 PCA tests
carrying the same names they carried before task 1.1.

The one that mattered most is the divisor the kinship will use. It was
tested only on diploid variants whose mean dosage is 1, where `p` is 0.5
and `p * (1 - p)` and `p * p` are the same number, so two wrong formulas
passed all 607 tests: `p * p` in place of `p * (1 - p)`, and the ploidy
read as 2 instead of taken from the variant. Both were run and both left
`607 passed; 0 failed`. Work package 2 would have built the kinship on a
divisor that nothing checked. The tests added use a diploid variant of
`p = 0.125`, a tetraploid and a haploid one, with the values computed by
hand and checked against pyNei's `_calc_dosages` with
`sqrt(ploidy * freqs * (1 - freqs))`. Both mutations were applied again
after the fix and each now fails
`variant::tests::the_divisor_under_hardy_weinberg_reads_the_ploidy_and_the_frequency`
and nothing else.

The rest, in the order of what they would have cost:

- The pass refused a ploidy above 254 with an error whose message opens
  "the principal components of the variants cannot be taken on this
  dataset", from the pass the kinship is to share. No user could reach it,
  because `pca_of_variants` refuses the ploidy at its own entry before a
  block is read, but the kinship will. It now has a case of its own,
  `Error::VariantPloidyTooLarge`, in the group of the errors of a variant,
  a `ValueError` in Python like its neighbour. `MAX_PLOIDY_OF_THE_VARIANTS`
  moved to `variant.rs` with it, so `variant.rs` no longer imports from
  `pca.rs`: it had been importing from the module above it, which would
  have pulled `pca` in behind every caller of the row.
- The two buffers of a row, 255 dosage counters and 256 values, are sized
  for a ploidy of 254 and nothing said so. Raising the limit to 255 gave a
  genotype with an allele missing a value of 16.03 where 0.0 is right, with
  no error and no panic. Three relations are now asserted when the crate
  compiles; raising the limit stops the build with "the dosages of the
  largest ploidy, which are the ploidy and one more, each need a counter".
- `codes.resize` could be deleted with all 607 tests passing, because every
  fixture built its buffers at the size of its row. A row longer than the
  buffers would have been left half written.
- The test of a variant with more than two alleles used a variant of four,
  so tightening the refusal from more than two to more than three left it
  green.
- `DosageOptions` and `DosageScale` were public, so a crate outside the
  workspace could build one and pass it nowhere; they are the crate's own
  now. Two helpers had been widened to the crate for no caller.
- The divisor read the ploidy out of a loop bound, so a caller that passed
  anything else got a wrong number and no error; it takes the ploidy.
- Nothing asserted the middle clause of the message deliverable 3 asks to
  keep unchanged, so it could have been rewritten with nothing going red.
- The tail of the row pass is copied into the benchmark and nothing
  compared the two; a test now asserts they give the same row bit for bit.
- Seven doc comments still said PCA or private, one in the Python binding
  claimed the row pass is walked by every calculation that turns a variant
  into dosages, which the r² and the filters disprove by having their own.

### What was not taken

- **The six tests of the moved code stay in `pca.rs`.** Two reviewers asked
  for them to move to `variant.rs`, where the code they test now lives.
  Moving them would drop `cargo test -p popnei --lib pca -- --list` below
  47 and so destroy the check deliverable 2 rests on, which is the only
  evidence that the move changed no number. It is worth doing once that
  check has served, which is when work package 2 is reviewed.
- **The JavaScript message names `transform_to_biallelic` where the
  TypeScript option is `transformToBiallelic`**, so a user greps for a name
  that is not in their code. It is older than this plan and it is what a
  user of the principal components sees today, so it is the owner's and it
  is asked of them below. The crate already rewrites `maxNumVars` and
  `maxAllowedMaf` this way, so there is a pattern to follow.
- **The Python binding maps 29 of the 99 cases of `Error` through a
  wildcard**, so a case added later becomes a `ValueError` without anyone
  choosing that. `Error` is `#[non_exhaustive]`, which makes the wildcard
  compulsory, and every case that falls through it is right today. It is
  older than this plan and nothing of this plan rests on it.

### What was changed in the plan

- Deliverable 2 asked for 604 tests in the core crate. That was the count
  before the work package, and deliverable 3 asks for tests that did not
  exist, so the two could not both be met. The check now reads the list of
  test names of `pca.rs`, which is what says that no test was dropped, and
  asks for no fewer than 604 in the crate.
- Task 2.0 is new. Since task 1.1 the three functions that drive a whole
  block of variants hold nothing of the principal components and take a
  `DosageOptions`, and the kinship needs the same drive over a block, so
  task 2.1 would have copied about 150 lines including the arm for the
  threads and the arm for WebAssembly, which is a second place where the
  two can fall out of step.
- The check of deliverable 1 was replaced before the work started, for the
  reason above.

### What the owner should know

- The spec said the kinship refuses two datasets "both for reasons the
  PCA's row pass already refuses them for". Only the ploidy is refused
  there; the limit of 46340 individuals is checked in `pca_of_variants`,
  which the kinship will never call. Whoever wrote work package 2 from that
  sentence would have got the ploidy for free and lost the other in
  silence. The spec now says which check is where.
- `DosageScale::OfHardyWeinberg` has no caller in the library until the
  kinship exists, so it carries an expectation that it is dead code. The
  commit that builds the kinship has to delete that expectation or the
  build fails, and the text of the expectation says so.
- The worked example of "How it is verified" in the spec, and both of its
  kept variants, have `p = 0.5`. The cargo test work package 2 builds from
  it inherits the blind spot this review found, so the eleven plink2
  literals of that section are what will close it there.

### How the work went

Task 1.1 took one subagent 197000 tokens over 112 tool calls and 15
minutes. The seven reviewers took 756000 tokens between them, from 79000
for the binding to 122000 for the numbers, and ran in parallel in about 13
minutes. The fixes took the subagent that wrote the code another 119000
tokens over 20 minutes; sending them back to it rather than to a fresh one
cost nothing in re-reading.

Sending all seven categories rather than choosing among them was worth it
here: the finding that mattered most came from `tests`, which was the
slowest of the seven and the only one that mutates the code, and three of
the others found the PCA-named ploidy error from three different sides.
