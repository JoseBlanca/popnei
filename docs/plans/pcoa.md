# Plan: the principal coordinates of distances and Lingoes' correction

27 September 2026. Draft, for the owner's approval. It builds from one
spec, the item "The principal coordinates of distances" of
`docs/specs/pca.md`, with its part of "The Rust interface", as it stands
on the branch `spec/pcoa` at the commit that adds this plan. It is carried
out on the branch `plan/pcoa` in the worktree `.claude/worktrees/pcoa`,
with its work report in `docs/reports/pcoa.md`.

When it is done a user calls, natively, under pyodide and from
TypeScript, `do_pcoa` on a `Distances`, `do_pcoa_from_variants` on a
`Variants`, and `correct_dists_by_lingoes` on a `Distances`. A matrix
that is not Euclidean is refused, with a message that names the
correction, or corrected when asked. The numbers are R's `ape::pcoa` to
1e-9, on the literals of the spec. popnei_web gets `doPcoaFromVariants`
for its stage 4.

## In and out

In: the three functions through the core, the two binding crates and the
two packages; the limit of the browser, measured; the README of the
TypeScript package.

Out, with where it goes: Cailliez's correction, which the owner did not
choose ("Not in this spec" of the spec); the release of the TypeScript
package, which the owner is asked for when this plan is merged; a
randomized PCoA for datasets beyond the limit of the browser, which
nothing asks for.

Open 7 of the spec, whether the correction is also an argument of the two
PCoAs, is not answered. The tasks follow its "meanwhile" and build both.
An answer of the function alone takes the argument out of tasks 1.2, 1.3,
1.4, 2.1, 2.2 and 2.3: a task not yet started is built without it, and
taking it out of the tasks already committed is one new task, before the
next one starts. Opens
1 to 6 of the spec are of the PCA and change nothing here, except that
Open 3, the components with no variance, is followed here as there.

## What has to be in place

- The spec and the reference data on the branch:
  `tests/reference/pca/pcoa_reference.R` and its 17 `*.pcoa.r.*.tsv` and
  `*.lingoes.r.*.tsv` files, which the script gives again byte for byte
  with R 4.6.1 and ape 5.8.1 at `/opt/homebrew/bin/Rscript`, checked on 27
  September 2026 by the spec's review. The tests need neither R nor ape.
- The branch where the checks of "Before the work is called done" of the
  `coding` skill pass. Run on 27 September 2026 at the commit of the
  spec: `cargo fmt`, `clippy` and both `wasm-check` aliases clean; 1335
  cargo tests of the workspace and 1185 without default features,
  passing; ruff clean; 618 pytest tests passing; `npm test` in
  `js/popnei` gives `tests 474`, `fail 1`. The one node test that fails
  does so on `main` too, `a kinship that does not tell the two variances
  apart gives none of them` of `test/gwas.test.ts`, by 1.1e-14 against a
  tolerance of 7e-15, as the board message of the merge of 27 September
  says. It is not this plan's; the report counts it apart.
- pyNei at ef0ca6e, which `pyproject.toml` names.

The checks of the plan fail today: `cargo test -p popnei --lib pca::pcoa
-- --list` prints `0 tests`; `uv run pytest tests/test_pcoa.py` finds no
file; `js/popnei/test/pcoa.test.ts` does not exist.

The code goes in `crates/popnei/src/pca/pcoa.rs`, a module of `pca`
declared in `crates/popnei/src/pca.rs`, so that its tests are named
`pca::pcoa::`. The sign rule and the threshold of the components are the
PCA's functions, called and not copied. The errors are new cases at the
end of the error enum of the core.

## Work package 1: the PCoA of a distance vector and the correction

What it gives: `do_pcoa(dists, correct_by_lingoes)` and
`correct_dists_by_lingoes(dists)` in Python on a `Distances`, and
`doPcoa` and `correctDistsByLingoes` in TypeScript, with the results and
the refusals of the spec.

Deliverables:

1. `cargo test -p popnei --lib pca::pcoa -- --list` names tests at `pcoa`
   and `correct_dists_by_lingoes` of "The Rust interface" for: the worked
   example of "How it is verified", refused without the correction with
   its counts, and with it the table, the percentages, the constant and
   `negative_eigenvalues_percent`; `correct_dists_by_lingoes` of it, the
   constant and the first corrected distance, and `pcoa` of the corrected
   vector, the same table within 1e-9; the twin, the ten distances with a
   sixth individual at distance 0 from the fifth, which is how the spec
   checks two clones, its 4 components against
   the whole of `small_twin.lingoes.r.*.tsv`; the Kosman distances of
   `four_alleles.gdkosman.tsv`, 39 components, with and without the
   correction, against `four_alleles.pcoa.r.*.tsv`; and each error of
   "Errors and the cases pyNei asserts" that a vector reaches, the pair
   with no distance with its counts and positions, and its tie. 12 tests at
   least, and they pass with and without the default features.
2. `uv run pytest tests/test_pcoa.py` runs tests that compare `do_pcoa` of
   popnei and of pyNei on the Kosman distances of `four_alleles.vcf.gz`,
   as the spec's last paragraph of "How it is verified" says; that
   `do_pcoa` of the ten distances is refused with a message that names
   `correct_dists_by_lingoes`, and gives the table with
   `correct_by_lingoes`; `correct_dists_by_lingoes`, its fields and the
   names and `pass_stats` of the `Distances` it gives; the frame and the
   series of `PCoAResult`, with `PC0` and its zeros; and the message of a
   pair with no distance with the names. They pass.
3. `npm test` in `js/popnei` runs a `test/pcoa.test.ts` that checks the
   worked example at `doPcoa` and `correctDistsByLingoes`, and the
   refusal with the names of TypeScript, `correctByLingoes` and
   `correctDistsByLingoes`; and the limit of the browser at the exported
   function of the binding, 8695 individuals taken and 8696 refused, as
   `test/pca.test.ts` checks the PCA's. The count of `fail` is that of
   the start.
4. The checks of the `coding` skill pass.

Stands on: the `linalg` crate, the PCA's sign rule and threshold, and the
`Distances` of both packages, all on `main`.

Tasks:

- [ ] 1.1 In the core, the PCoA without the correction: `Pcoa`,
      `PcoaOptions` and `pcoa` of "The Rust interface", from "What it
      gives" and "How it runs" of the spec's item; B from the vector,
      which is dropped before the eigendecomposition; the refusals of a
      vector; the error cases with what "The Rust interface" says they
      carry. With `correct_by_lingoes` true it returns the error that the
      correction is not built yet, which task 1.2 replaces. The tests of
      deliverable 1 that need no correction first. Serves deliverables 1
      and 4.
- [ ] 1.2 Lingoes' correction, in `pcoa` and in
      `correct_dists_by_lingoes`, as "Lingoes' correction" and the second
      paragraph of "How it runs" say: the eigenvalues shifted by c, and the
      eigenvectors of the eigenvalue 0 projected orthogonal to the vector
      of ones and made orthonormal again. A mistake here is a wrong number
      and not a crash, so it is its own commit, and the twin is what
      guards it: without the projection its test fails. Needs 1.1. Serves
      deliverables 1 and 4.
- [ ] 1.3 The Python functions: the bindings of `pcoa` and
      `correct_dists_by_lingoes` in `crates/popnei-python`, taking the
      vector of a `Distances`; `do_pcoa`, `correct_dists_by_lingoes`,
      `PCoAResult` and `LingoesCorrection` in `python/popnei/pca.py`,
      exported from `python/popnei/__init__.py`, with the names put on the
      messages of the pairs; `tests/test_pcoa.py`. From "What it is in
      Python and in TypeScript" and "Errors and the cases pyNei asserts".
      Needs 1.2. Serves deliverables 2 and 4.
- [ ] 1.4 The TypeScript functions: the bindings in `crates/popnei-js`,
      the limit of 8695 individuals beside the PCA's in
      `crates/popnei-js/src/pca.rs` with its exported function, the
      Python names of the new messages rewritten in camelCase in
      `crates/popnei-js/src/errors.rs`; `doPcoa`, `correctDistsByLingoes`
      and their results in a new `js/popnei/src/pcoa.ts`, exported from
      both entry points, `node.ts` and `web.ts`, with every refusal in the
      `@throws`;
      `test/pcoa.test.ts`. Needs 1.2; can run side by side with 1.3.
      Serves deliverables 3 and 4.

What could go wrong: the rule for the eigenvectors of the eigenvalue 0 is
the part of the spec that no program outside checks directly; ape
decomposes the corrected matrix again, and the twin is the one case where
the two routes could differ. The review found them 6e-16 apart with
numpy's QR. Whether faer and LAPACK give the eigenvalue 0 of the twin
inside the threshold on both backends is known only when
`--no-default-features` runs the test. When it fails on one backend
alone, the threshold is the spec's and is not moved in the code: the
orchestrator writes into the report the eigenvalues the two backends
gave against the threshold, and stops for the owner, since the answer
changes the spec.

## Work package 2: the PCoA of the variants

What it gives: `do_pcoa_from_variants(variants, min_num_snps,
correct_by_lingoes)` in Python and `doPcoaFromVariants` in TypeScript,
one pass over a `Variants` with its steps, with the counts of the pass,
and `numPassesOf("doPcoaFromVariants")` of 1.

Deliverables:

1. `cargo test -p popnei --lib pca::pcoa -- --list` names, besides those
   of work package 1, tests at `pcoa_of_variants` on a reader over the VCF
   for: the panel refused without the correction, with its 44 of 200 and
   its percent; the panel with the correction, its literals and the whole
   of `panel.lingoes.r.*.tsv`, and its constant; the panel with
   `min_num_vars` 1105, the error with 35 of 19900, the pair (1, 82) and
   the individual 82 in 17; `four_alleles.vcf.gz`, 39 components; and a
   reader of one individual refused before any block is read. 5 tests
   more at least, and they pass with and without the default features.
2. `uv run pytest tests/test_pcoa.py` runs, besides those of work
   package 1: `do_pcoa_from_variants` with `correct_by_lingoes` against
   pyNei as the spec's last paragraph of "How it is verified" says; the
   refusal without it, naming `correct_by_lingoes`; the message of
   `min_num_snps` 1105 with `s001`, `s082` and 17; `pass_stats` with
   `num_vars` 1200, and with a `filter_by_maf` step the count that filter
   kept. They pass.
3. `npm test` runs tests at `doPcoaFromVariants` on
   `tests/reference/dists/panel.vcf.gz`: the literals with
   `correctByLingoes`, the refusal without it, the message of
   `minNumSnps` 1105, and `passStats`; `numPassesOf("doPcoaFromVariants")`
   is 1; and the lists of `test/progress.test.ts` and
   `test/num_passes.test.ts` have the new consumer. The count of `fail` is
   that of the start.
4. The checks of the `coding` skill pass.

Stands on: work package 1, and `calc_kosman_sums` and the chain of steps
of the binding crates, on `main`.

Tasks:

- [ ] 2.1 In the core: `VariantPcoaOptions`, `PcoaOfVariants` and
      `pcoa_of_variants`, as "How it runs" says, B built from the sums of
      the pass and the sums dropped before the eigendecomposition, the
      pairs with no distance checked on the sums, and the refusal of fewer
      than two individuals before the pass. The tests of deliverable 1.
      Serves deliverables 1 and 4.
- [ ] 2.2 The Python function: the binding that opens the reader from the
      source and the steps and reads the counts, as that of
      `calc_pairwise_kosman_dists` does; `do_pcoa_from_variants` in
      `python/popnei/pca.py`; the tests of deliverable 2. Needs 2.1.
      Serves deliverables 2 and 4.
- [ ] 2.3 The TypeScript function: the binding over `the_run_of` with a
      new variant of `Consumer` in `crates/popnei-js/src/source.rs`,
      `doPcoaFromVariants` in `js/popnei/src/pcoa.ts` with the limit
      checked before the run, its name in `ConsumerName` of
      `js/popnei/src/passes.ts`, the lists of the two tests, and the
      README of `js/popnei`, whose counts of the calculations and of the
      consumers grow by the new functions. Needs 2.1; can run side by side
      with 2.2. Serves deliverables 3 and 4.

## Work package 3: the limit of the browser, measured

What it gives: the number of individuals above which the three functions
refuse in TypeScript, measured as the PCA's was, in place of the 8695
that "How it runs" works out as an upper bound.

Deliverables:

1. A script in `js/popnei/bench/`, kept, that writes a VCF of n
   individuals and a few hundred variants, runs `doPcoaFromVariants` with
   `correctByLingoes` under node, and prints how far the memory of wasm
   grew. Its output at 3000 individuals and at the largest n that runs,
   found by halving between 8695 and 9500, is in the report with the
   command.
2. "How it runs" of the spec gives the measured factor and the limit it
   gives, with the same margin below the smallest n that trapped as the
   PCA took, in a commit of its own; then the constant, its doc comment,
   the `@throws` and the test of deliverable 3 of work package 1 at the
   new boundary, in the next commit. When the measured limit is 8695 the
   spec says it was measured, and nothing else moves.
3. The checks of the `coding` skill pass.

Stands on: work package 2.

Tasks:

- [ ] 3.1 The script and the measurement, one at a time on a machine with
      nothing else building. Serves deliverable 1.
- [ ] 3.2 The spec, then the constant. Needs 3.1. Serves deliverables 2
      and 3.

What could go wrong: the sums of the pass given back before the
eigendecomposition can leave a hole that the allocator of wasm does not
reuse, and then the peak is above the 56.8 bytes a cell of the spec and
the limit below 8695. That is what the measurement is for.

## How the whole plan is checked

Beyond its work packages: the wheel for pyodide builds and its smoke test
passes, `scripts/build_pyodide_wheel.sh && node tests/pyodide/smoke.mjs`;
and `npm run test:browser` in `js/popnei` runs the package in Chromium.
Then the owner is asked for the release of the TypeScript package, the
tag after `js-v0.1.0-dev.2`.
