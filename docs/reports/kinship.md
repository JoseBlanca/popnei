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
