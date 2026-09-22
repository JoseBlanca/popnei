# Plan: the principal component analysis and the linalg crate

22 September 2026. Approved by the owner in chat on 22 September 2026
and under way since that day, on the branch `plan/pca` in the worktree
`.claude/worktrees/pca`, with its work report in `docs/reports/pca.md`.
It builds from two specs committed there at 5852540:

- `docs/specs/linalg.md`, the crate `crates/popnei-linalg`: the product
  of a matrix with itself, the eigendecomposition of a symmetric matrix
  and the product of two matrices, on BLAS and LAPACK natively and on
  faer in wasm and behind a cargo feature;
- `docs/specs/pca.md`, the `pca` module of the core crate: the PCA of a
  table, `do_pca`, and the PCA of the variants of a `Variants`,
  `do_pca_from_variants`, with their TypeScript twins.

When it is done a user calls `do_pca` on a pandas frame and
`do_pca_from_variants` on a `Variants`, natively, under pyodide and from
TypeScript, and gets pyNei's projections and percentages to 1e-9, the
weights of the first 10 components, and the counts of the pass. The PCA
of 100000 variants x 1000 individuals runs in about 0.3 s where pyNei
takes 8.2 s, measured, or the report says by how much it missed.

## In and out

In: the linalg crate with its two backends and its wasm builds; the two
PCA functions through the three layers; the measurements of "Speed" of
both specs; and the four documents that `docs/specs/linalg.md` says the
crate changes.

Out, with where it goes: the principal coordinates, a later item of
`docs/specs/pca.md`; the operations of the GWAS, later items of
`docs/specs/linalg.md`; the Linux and Windows builds of OpenBLAS, the
plan that first makes them; a block of genotypes held in 2 bits each, as
plink does, section 4 of the architecture.

The open points of the specs are unanswered and each has its
"meanwhile", which the tasks follow. The task each answer would change:
Open 1 of the PCA spec, whether `num_prin_comps` cuts the projections,
task 3.2; Open 2, a variant with no called genotype, task 3.1; Open 3,
the components with no variance, tasks 2.1 and 3.1; Open 4, the spelling
of `standardize_data`, task 2.2; Open 5, whether the standard deviation divides by n or by n - 1, tasks
2.1 and 3.1;
Open 1 of the linalg spec, the vector instructions of wasm, task 4.3:
that task is done when the owner answers on, and skipped, and said so
in the report, when they answer off or have not answered by the end of
work package 4.

## What has to be in place

- The two specs and the reference data of `tests/reference/pca/` on the
  branch: `git log --oneline -1 -- docs/specs/pca.md` gives 5852540.
- The branch at that commit, where the five commands of "Before the work
  is called done" of the `coding` skill pass, `cargo wasm-check` passes,
  `npm test` in `js/popnei` gives `tests 126`, `fail 0`, and
  `scripts/build_pyodide_wheel.sh && node tests/pyodide/smoke.mjs` exits
  with 0. Run on 22 September 2026: 306 cargo tests, 174 pytest tests,
  126 node tests, all passing; the wheel built and the smoke test passed
  the same day, with `RUSTFLAGS` set and unset.
- On the machine: `plink2`, `Rscript`, `node`, `wasm-bindgen`, emsdk at
  `~/devel/emsdk` and pyodide-build at `~/devel/pyodide-venv`, checked
  with `which` and `ls` on 22 September 2026; R 4.6.1 with nothing
  installed for it. The tests need none of them: the literals are in the
  specs and in `tests/reference/pca/`.
- pyNei at ef0ca6e, which `pyproject.toml` names.

The checks that say the work is done fail today: `cargo test -p
popnei-linalg` says the package is not in the workspace; `cargo test -p
popnei --lib pca -- --list` prints `0 tests`; `uv run pytest -k pca`
exits with 5; `js/popnei/test/pca.test.ts` does not exist.

What the reports of the earlier plans learned holds: tasks that write
the workspace manifest or `Cargo.lock` run one after another, every
reviewer works in a worktree of its own, the count of the tests of the
library is read with `--lib`, and timings are taken one at a time with
nothing else building. The board, `.claude/board/`, has the start
message of this branch; a change to a shared file is posted there as its
README says.

## Work package 1: the linalg crate

What it gives, in the words of work packages 2 and 3: the three
operations of `docs/specs/linalg.md` as functions of a crate the core
crate depends on, with the same numbers from BLAS and from faer, built
for the native target and the two wasm targets. It stops inside the core
because the crate has no Python side; work package 2 makes the first
comparison with pyNei through it.

Deliverables:

1. `cargo test -p popnei-linalg` and `cargo test -p popnei-linalg
   --no-default-features` both pass, and `cargo test -p popnei-linalg --
   --list` names the tests of the four checks of "How it is verified" of
   the linalg spec, one test or more for each: the product of a matrix
   with itself, with an A of no rows; the product of two matrices; the
   3 x 3 eigendecomposition; the 1000 x 1000 eigendecomposition from the
   generator of the spec, with its six literals. 8 tests at least.
2. Every error of "Errors" of the linalg spec has a test: a dimension
   that does not match, a buffer too short, a c or an n of 0, a value
   that is not finite.
3. `cargo wasm-check` checks the crate as well as the core, and passes;
   the crate's manifest has the feature of `gemm` under
   `cfg(target_family = "wasm")` as "The wasm builds and the vector
   instructions" says; `scripts/build_pyodide_wheel.sh && node
   tests/pyodide/smoke.mjs` exits with 0 and `npm test` in `js/popnei`
   still gives `fail 0`, with the core crate depending on the new crate.
4. The four documents read as the crate is: `grep -n "linalg"
   docs/architecture.md docs/objectives.md .claude/skills/coding/SKILL.md`
   finds no "module" beside it and section 8 of the architecture lists
   the crate.
5. The five commands of the `coding` skill pass with the new member of
   the workspace.

Stands on: nothing of this plan. `blas` 0.23, `lapack` 0.20, `blas-src`
0.14 with `accelerate`, `lapack-src` 0.13 with `accelerate`, faer 0.24.4
and `gemm` 0.19, all on crates.io on 22 September 2026.

Tasks:

- [x] 1.1 The crate with the BLAS backend: `crates/popnei-linalg` in the
      workspace, its `Error`, the three functions of "The Rust interface"
      and the checks of dimensions and finiteness of "Errors", on the
      routines that "Layout, half and the backends" names, each `unsafe`
      block with its comment. From "The three operations" of
      `docs/specs/linalg.md`. The tests of the two products, the 3 x 3
      eigendecomposition and the errors, written first. Serves
      deliverables 1, 2 and 5.
- [x] 1.2 The faer backend and the feature `blas`: the BLAS crates as
      optional dependencies of the targets that are not wasm, the backend
      chosen as "What it gives" of "The crate, its backends and the
      builds" says, and the 1000 x 1000 test from the generator of "How
      it is verified", with its literals, run with both commands of
      deliverable 1. Needs 1.1. Serves deliverables 1 and 5.
- [x] 1.3 The wasm builds: the `wasm-check` alias extended to the crate,
      the feature of `gemm` under the wasm cfg, the core crate depending
      on the new crate with nothing calling it yet, the wheel built and
      the smoke test and `npm test` run, as the `coding` skill asks when a
      change touches what wasm builds differently. From "The wasm builds
      and the vector instructions" and "Threads". Needs 1.2. Serves
      deliverable 3.
- [x] 1.4 The four documents: section 8 and the `linalg` row of section 9
      of `docs/architecture.md`, the line of `docs/objectives.md`, the
      two lines of the `coding` skill, each as the opening of
      `docs/specs/linalg.md` lists them and as the `writing` skill asks;
      the message on the board that says the shared files changed. Can
      run side by side with 1.2 and 1.3. Serves deliverable 4.

What could go wrong: the Python binding crate is an extension module
that maturin links, and whether `-framework Accelerate` reaches its link
line through `accelerate-src`, the crate that emits that linker
argument, was not tried; the trial linked a binary. If the wheel or the
extension module does not link, the way through is a build script of
the linalg crate that emits the same argument, and if that fails too the
owner is asked.
faer 0.24.4 built for emscripten in the spike of `docs/rust_core.md` and
for `wasm32-unknown-unknown` in the trial of the PCA spec, but not inside
the two binding crates. Whether cargo takes an optional dependency under
a target table with `--no-default-features` the way "What it gives"
assumes is checked by the two commands of deliverable 1.

## Work package 2: the PCA of a table

What it gives: `do_pca(data, center_data, standardize_data)` in Python,
on a pandas frame, and `doPca` in TypeScript, with the `PCAResult` of
the spec, giving pyNei's numbers for iris.

Deliverables:

1. `cargo test -p popnei --lib pca -- --list` names tests at `pca` of
   "The Rust interface" of `docs/specs/pca.md` for: the four iris
   literals of "How it is verified" of "The PCA of a table", standardized
   and not, within 1e-9; the two rules of "What both analyses compute", the sign of a component
   and that a component with no variance is not given, on pyNei's 3 x 3 table of "Errors and
   the cases pyNei asserts" of that part, which gives 2 components; and
   each error of that part. 6 tests at least, and they pass.
2. `uv run pytest tests/test_pca.py` runs tests that call `do_pca` on
   `tests/reference/pca/iris.tsv` with both values of `standardize_data`
   and compare with pyNei's `do_pca` at ef0ca6e, all 4 components of
   the projections, the percentages and the princomps, within 1e-9 after
   the sign rule is applied to pyNei's; that `PCAResult` has the index
   and the columns of the frame and `pass_stats` of `None`; and that the
   error of a trait with no variance names the trait. They pass.
3. `npm test` in `js/popnei` runs a `test/pca.test.ts` that checks the
   iris literals of row 0 and the percentages at `doPca`, within 1e-9,
   and gives `fail 0`.
4. The five commands of the `coding` skill pass.

Stands on: work package 1.

Tasks:

- [x] 2.1 In the core: `Pca`, `PcaOptions` and `pca` of "The Rust
      interface", built from "What both analyses compute" and "The PCA of
      a table" of `docs/specs/pca.md`, over the three functions of the
      linalg crate. The tests of deliverable 1 first, with the literals
      of `tests/reference/pca/iris.r.*.tsv` and
      `iris_not_standardized.r.*.tsv`. Serves deliverables 1 and 4.
- [x] 2.2 The Python function: the binding of `pca` in
      `crates/popnei-python`, taking a float64 array; `do_pca` and
      `PCAResult` in `python/popnei/pca.py`, with the names of the frame
      put on the result and on the error of "Errors and the cases pyNei
      asserts"; `tests/test_pca.py` with the comparison with pyNei.
      From "What it is in Python and in TypeScript" of "The PCA of a
      table". Needs 2.1. Serves deliverables 2 and 4.
- [x] 2.3 The TypeScript function: the binding in `crates/popnei-js`,
      `doPca` in `js/popnei/src/pca.ts` with the result of "What it is
      in Python and in TypeScript", and `test/pca.test.ts`. Needs 2.1;
      can run side by side with 2.2. Serves deliverable 3.

## Work package 3: the PCA of the variants

What it gives: `do_pca_from_variants(variants, transform_to_biallelic,
num_prin_comps)` in Python and `doPcaFromVariants` in TypeScript, on a
`Variants` with its steps, in two passes, with pyNei's projections and
percentages on the reference panel, the weights of the first 10
components and the counts of the pass.

Deliverables:

1. `cargo test -p popnei --lib pca -- --list` names, besides those of
   work package 2, tests at `pca_of_variants` of "The Rust interface" on
   readers the tests build for: the worked example of "How it is
   verified" of "The PCA of the variants", from `worked.vcf`, with
   `num_prin_comps` 3 and 0, within 1e-9; `worked3.vcf` without
   `transform_to_biallelic`, the error, and with it, 3 components with
   the numbers of `worked3.r.*.tsv`; the panel, `sim_missing.vcf`, its
   five literals; the five cases pyNei asserts, from "Errors and the
   cases pyNei asserts", among them that the result does not change
   with the size of the blocks within 1e-10; each error of that part; and
   the error of a second pass whose variants differ. 12 tests at least,
   and they pass.
2. `uv run pytest tests/test_pca.py` runs tests that call
   `do_pca_from_variants` on `sim_missing.vcf` and `worked.vcf` and
   compare with pyNei's at ef0ca6e as "How it is verified" says, 10 and
   3 components within 1e-9 after the sign rule; that `pass_stats` has
   `num_vars` 1200 for the panel and, with a `filter_by_maf` step on the
   `Variants`, the count that filter kept; and that `num_prin_comps` 0
   gives `princomps` with no rows and the used variants as columns. They
   pass.
3. `npm test` runs a test at `doPcaFromVariants` on the worked example,
   its projections and `usedVars`, and `passStats`, and gives `fail 0`.
4. The five commands of the `coding` skill pass.

Stands on: work packages 1 and 2, and the readers of the core: the VCF
reader, `reblock` and the chain of steps that `write_vars` of the two
binding crates builds.

Tasks:

- [x] 3.1 In the core, the first pass: the standardizing of a block as
      "What it gives" and "How it runs" of "The PCA of the variants"
      describe it, the counts of the dosages per row with rayon, the
      check of the alleles, the rows without variance left out, the
      product added to the individuals x individuals matrix through the
      linalg crate from outside rayon, its eigendecomposition, and
      `pca_of_variants` without a second reader, which "The Rust
      interface" allows when no weights are asked for. The tests of deliverable 1 that need no weights first: the
      worked example with `num_prin_comps` 0, `worked3`, the cases pyNei
      asserts, the errors. Serves deliverables 1 and 4.
- [ ] 3.2 The second pass: the weights of the first `num_prin_comps`
      components as "How it runs" says, the check that the variants of
      the two passes are the same, and the tests of deliverable 1 with
      weights: the worked example with 3, the panel's five literals, the
      error of a second pass that differs. Needs 3.1. Serves deliverables
      1 and 4.
- [ ] 3.3 The Python function: the binding that opens the two readers
      from the source and the steps, lends them to the core and reads the
      counts, as `write_vars` of `crates/popnei-python/src/vars.rs` does
      for one; `do_pca_from_variants` in `python/popnei/pca.py`, with the
      `PCAResult` of "What it is in Python and in TypeScript" of "The
      PCA of the variants" and its `pass_stats`; the tests of deliverable
      2. Needs 3.2. Serves deliverables 2 and 4.
- [ ] 3.4 The TypeScript function: the binding in `crates/popnei-js` on
      the pattern of its `write_vars`, `doPcaFromVariants` in
      `js/popnei/src/pca.ts`, and its test. Needs 3.2; can run side by
      side with 3.3. Serves deliverable 3.

What could go wrong: the loop that standardizes a block is the one the
spec expects the compiler to vectorize, and whether it did is not known
until task 4.1 looks; a loop written differently is right and slower,
and the target of "Speed" is what says so. A test that builds the panel
through the VCF reader takes the reader's time in every run of the
tests; `sim_missing.vcf` is 1 MB.

## Work package 4: the measurements

What it gives: the numbers of "Speed" of the two specs measured on the
code, in `docs/reports/pca-measurement.md`, and the `simd128` flag in
the two wasm build commands if the owner turned it on.

Deliverables:

1. `docs/reports/pca-measurement.md` has, for 100000 variants x 1000
   individuals, the time of `do_pca_from_variants` from a vars file with
   `num_prin_comps` 0 and 10, on one thread, `VECLIB_MAXIMUM_THREADS=1`
   and `RAYON_NUM_THREADS=1`, and with the threads the machine gives, which for Accelerate means
   those two variables unset;
   the time of pyNei's `do_pca_from_variants` and of plink2's `--pca 10
   meanimpute --threads 1` on the same variants, `meanimpute` giving a missing genotype the mean
   of its variant as popnei does; and whether the standardizing loop was
   vectorized, read from the machine code that `cargo asm` prints for it. The target is
   0.3 s natively, and the report says whether it was met.
2. The same report has the time of `doPcaFromVariants` under node on the
   wasm package, on 20000 x 1000 at least and on 100000 x 1000 if the
   memory of node takes it, against the targets of "Speed" of the PCA
   spec, 5 s with `simd128` and 7 s without; and the wheel built and the
   smoke test passing.
3. When Open 1 of the linalg spec is answered as on: `build:wasm` of
   `js/popnei/package.json` and `scripts/build_pyodide_wheel.sh` set
   `RUSTFLAGS` with the flag as "The wasm builds and the vector
   instructions" says, the two builds pass their tests, and the report
   has the wasm times with it.

Stands on: work package 3. The dataset: `uv run --no-project --with
numpy python crates/popnei/benches/make_big_vcf.py <dir>/big.vcf`
outside the repository, 3 s, 403 MB, then `write_vars` of it, so that
the reading is the vars file's; pyNei timed with a script beside
`crates/popnei/benches/time_pynei.py` on the same VCF; plink2 on the
VCF with `--vcf`.

Tasks:

- [ ] 4.1 The native measurement and the vectorization check, one
      timing at a time with nothing building, written into the report
      as "Speed" of `docs/specs/pca.md` states its numbers. Serves
      deliverable 1.
- [ ] 4.2 The wasm measurement: the wasm package built and timed under
      node, the wheel built and smoke tested. Needs 4.1. Serves
      deliverable 2.
- [ ] 4.3 The flag, only when Open 1 of the linalg spec is answered as
      on: the two build commands, the builds and their tests, the times
      again. Needs 4.2. Serves deliverable 3.

What could go wrong: the target of 0.3 s counts 20 blocks at 12 ms on
one thread and an eigendecomposition of 0.04 s, from the two trial
crates that the "Speed" parts of the specs describe, and the code has the reading of the vars file, the allocation of the
buffer of each block and the second pass besides; a miss is reported
with where the time went, from `sample`, the sampling profiler of macOS, and is not a reason to change
the code inside this plan.

## How the whole plan is checked

The five commands of the `coding` skill, `cargo wasm-check`, `npm test`
with `fail 0`, and the wheel with its smoke test, from a clean clone of
the branch, with the counts of the tests in the report; and the table of
deliverable 1 of work package 4 beside the targets of the two specs.
