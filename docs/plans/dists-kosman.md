# Plan: the Kosman distances between individuals

22 September 2026. Approved by the owner on 22 September 2026, on a
breakdown of four work packages from which he took out the first, a
`Variants` built from an array of genotypes: the tests run on small VCF
files written for them instead. Done on 22 September 2026, in the
worktree and on the branch below, with every task ticked and every
deliverable checked; `docs/reports/dists-kosman.md` is its work report,
and none of the three numbers of "Speed" is met, which the report puts to
the owner. This plan builds the first calculation
of popnei: the Kosman distance of every pair of individuals,
`calc_pairwise_kosman_dists`, and its `Distances` result. It builds from:

- `docs/specs/dists.md`, the item "Kosman distances between individuals",
  every part of it, and "The Rust interface" and "Speed";
- `docs/specs/filters.md`, "The Rust interface", for `chain_of`, the
  function that builds the chain of filters of a pass, which the binding
  crates already call for `write_vars` and `iter_blocks`.

The dists spec has no open point: the owner decided the ones it had on
22 September 2026, and they are written into it.

The branch is `plan/dists-kosman`, from `main` at 7d8366f, in the
worktree `.claude/worktrees/dists-kosman`; its first commit, 28c2f16, has
the spec and `docs/reports/kosman-method/`, the timing trial and the
scripts behind the literals of the spec.

## In and out

In:

- `tests/reference/dists/make_reference.py`, the four VCF files it writes
  and the output of R it stores.
- The calculation in the core, `calc_kosman_sums` and `KosmanSums`; the
  Python `Distances` and `calc_pairwise_kosman_dists`; the TypeScript
  `calcPairwiseKosmanDists` and its `Distances`.
- The measurement of "Speed" of the dists spec, and the wheel for
  pyodide built with the calculation in it.

Out, with where it goes:

- Jost's D between populations: the next item of `docs/specs/dists.md`,
  not written.
- The principal coordinate analysis of these distances:
  `docs/specs/pca.md`, another plan.
- The read ahead thread of section 3 of the architecture, which reads
  the next block while the calculation works on the one in hand.
  `docs/specs/block.md` leaves it for the first calculation, but no spec
  describes it, so this plan does not build it. Task 3.1 times the
  calculation on the vars file and the reading of that file alone, as the
  spec asks, and the second is what a read ahead thread could gain at
  most; the report carries it for the owner.
- `Variants.from_gt_array`, the `Variants` from an array of genotypes:
  an item of `docs/specs/variant.md` with no plan, until a user needs
  it.

## What has to be in place

- The branch `plan/dists-kosman` at 28c2f16 or later, in the worktree
  above. Run there on 22 September 2026: `cargo fmt --all --check` exit
  0; `cargo clippy --workspace --all-targets -- -D warnings` no warning;
  `cargo test --workspace` `306 passed`, 2 ignored; `cargo wasm-check`
  finished; ruff `18 files already formatted` and `All checks passed!`;
  `uv run maturin develop && uv run pytest` `174 passed`; in
  `js/popnei`, after `npm install`, `npm run build` and `npm test` 126
  tests passed, 0 failed. These are the counts that "more than" is
  counted from below. `bash scripts/build_pyodide_wheel.sh && node
  tests/pyodide/smoke.mjs` was not run when the plan was written; in a
  new worktree it needs `npm install` in `tests/pyodide` first.
- pyNei at ef0ca6e, which `pyproject.toml` has: `uv run python -c "from
  pynei import calc_pairwise_kosman_dists, load_vars"` exits with 0, and
  `/Users/jose/devel/pynei/test/gwas_reference/sim_missing.vars` is
  there, the panel of "How it is verified" of the dists spec.
- R 4.6.1 at `/opt/homebrew/bin/Rscript`, with nothing of what the
  reference needs in its libraries, `/opt/homebrew/lib/R/4.6/site-library`
  and the one of R itself: adegenet 2.1.11 and PopGenReport 3.1.3 were
  run from a scratch directory when the spec was written. Task 2.1
  installs adegenet into the site library with `install.packages`, and
  downloads the source of PopGenReport from CRAN, whose one function
  `gd.kosman` is read with `source()`, because the package itself does
  not install, as "How it is verified" of the dists spec says.
- The file of the measurement, outside the repository:
  `/Users/jose/devel/popnei-bench/big.vcf`, 403572954 bytes, 100000
  variants of 1000 individuals, made by
  `crates/popnei/benches/make_big_vcf.py`, and `big.vars`, 81356714
  bytes, written from it by `write_vars`. Both are there.

Every layer exists, so the five commands of "Before the work is called
done" of the `coding` skill run for every task, with `cargo wasm-check`,
and `npm run build && npm test` for a task that touches `crates/popnei-js`
or `js/popnei`.

The checks below fail today: `cargo test -p popnei --lib -- dists::
--list` prints `0 tests`, and `tests/test_dists.py`, `tests/reference/dists/` and
`js/popnei/test/dists.test.ts` do not exist.

Tasks run one after another unless the plan says that two can run side
by side: every task builds the same workspace in the same worktree.

## Work package 1: the reference

**What it gives.** The four VCF files and the numbers of R that the
tests of work package 2 are written against, and the script that makes
them again.

**Deliverables.**

1. The script and what it stores. Check: `uv run python
   tests/reference/dists/make_reference.py` writes, beside itself, the
   four datasets of "How it is verified" of the dists spec as gzipped
   VCFs, the panel read from pyNei's `sim_missing.vars` with pyNei and
   the three random ones from the seeds the spec gives, the distance
   vector and the number of variants of every pair that `gd.kosman`
   gives for each, and pyNei's distance vector for the two diploid ones;
   it refuses an adegenet that is not 2.1.11 and a PopGenReport source
   that is not 3.1.3; run again, `git status --short
   tests/reference/dists` shows nothing; and the first values of what it
   stores are the fourteen literals of the spec's table, with their n.

**What it stands on.** Nothing of this plan.

**Tasks.**

- [x] 1.1 `tests/reference/dists/make_reference.py` and what it stores,
  written as `tests/reference/filters/make_reference.py` is, from
  `ref_export.py`, `ref.R`, `poly_export.py` and `poly.R` of
  `docs/reports/kosman-method/`, which it replaces for the tests and
  which stay in the report as they are. From "How it is verified" of
  the dists spec. Serves deliverable 1.

**What could go wrong.** The panel comes out of pyNei's vars file
through pyNei, which reads `sim_missing.vars` at a path outside the
repository; the VCF the script writes is what the tests read, so the
path is needed once, when the script runs, and the script says so when
it is not there.

## Work package 2: the calculation, through the three layers

**What it gives.** A user calls `calc_pairwise_kosman_dists(variants,
min_num_snps=...)` in Python and gets a `Distances` with the distance of
every pair, its names and the counts of the pass, `pass_stats`; in
TypeScript,
`calcPairwiseKosmanDists`.

**Deliverables.**

1. The sets of bits and the counts of a pair. Check: cargo tests in
   `crates/popnei/src/dists.rs`, at the functions that build the sets of
   a block and count a pair, with the three worked examples of "How it
   is verified" of the dists spec as the first tests, the diploid, the
   tetraploid and the haploid, asserting k times the sum of d and n of
   every pair; and a block whose genotypes are all missing, where every
   pair has n = 0.
2. The calculation over a reader. Check: `cargo test -p popnei --lib --
   dists:: --list` prints 14 tests or more. Beside those of deliverable
   1, made at `calc_kosman_sums`, `KosmanSums::dist` and
   `KosmanSums::dists`: the fourteen literals of the spec's table on the
   four files of work package 1 read with the VCF reader, the integers
   exactly and the distances within 1e-9; the three worked examples
   through a reader written in the test, with `min_num_vars`, the
   `min_num_snps` of Python under its Rust name, 0, 3 and 4 for the
   diploid one; the 4 allele file in blocks of 7, 64, 65 and 300
   variants and in rayon pools of 1 and 4 threads giving the same
   integers; the two errors of "The Rust interface", a reader with no
   variant and a sum above `u32::MAX`, the second on sums set by hand;
   and an error of the reader given on as it is.
3. Python. Check: `uv run pytest tests/test_dists.py` passes, where
   today the file is not there, with: the comparison with pyNei of "How
   it is verified", exact, on the two diploid files and on the diploid
   worked example written as a VCF with the `write_vcf` fixture of
   `tests/conftest.py`, with `min_num_snps` of `None` and of 1125 for
   the panel; the tetraploid and haploid literals against R's numbers; every
   case of "Missing genotypes, pairs with no distance, and what pyNei
   asserts", pyNei's four tests among them, the one with the filter as a
   step with its `pass_stats`, the counts of the pass that every result
   of a consumer carries; the `ValueError` of a negative
   `min_num_snps`, of a vector of a wrong length and of a pass with no
   variant, whose message names the filter's counts when the steps kept
   none; `triang_list_of_lists` for the four individuals of "What pyNei
   does that is odd"; and `pass_stats` with the number of variants of
   the panel.
4. TypeScript. Check: `npm test` gives 0 failed and `test/dists.test.ts`
   runs the panel's five literals and the diploid worked example through
   `calcPairwiseKosmanDists`, with `passStats`, and the `Error` of a
   source with no variant.

**What it stands on.** Work package 1.

**Tasks.**

- [x] 2.1 The sets of bits of a block and the counts of a pair, in
  `crates/popnei/src/dists.rs`, with the worked examples as their tests.
  From "What it gives" and "How it runs" of the dists spec. Serves
  deliverable 1. Its failure would be silent, a wrong integer, so the
  three worked examples guard it and it has a commit of its own.
- [x] 2.2 `calc_kosman_sums` and `KosmanSums`, in the same file, the loop
  over the blocks of a reader with the accumulation, the pairs of a
  block on rayon natively and one after another in wasm, behind
  `cfg(not(target_family = "wasm"))` as the VCF reader does it, and the
  cargo tests of deliverable 2. From "How it runs" and "The Rust
  interface". Serves deliverable 2. Needs 2.1 and work package 1.
- [x] 2.3 The Python side: the function of `crates/popnei-python` that
  builds the chain with `chain_of`, runs the calculation and gives back
  the distances, the names and the counts of the pass, as `write_vars`
  of `vars.rs` does; `Distances` and `calc_pairwise_kosman_dists` in
  `python/popnei/dists.py`, exported from `popnei`; and the tests of
  deliverable 3 in `tests/test_dists.py`. From "Its Python function",
  "Missing genotypes, pairs with no distance, and what pyNei asserts"
  and "What pyNei does that is odd". Serves deliverable 3. Needs 2.2.
- [x] 2.4 The TypeScript side: the same in `crates/popnei-js`,
  `calcPairwiseKosmanDists` and `Distances` in `js/popnei/src/dists.ts`,
  and `test/dists.test.ts`. From the TypeScript paragraph of "Its
  Python function". Serves deliverable 4. Needs 2.2; can run side by
  side with 2.3.

**What could go wrong.** The pairs of a block are one work item each on
rayon, and 5e7 pairs for 10000 individuals is many items: the spec
leaves the split to the implementer, and the timing of work package 3
is what says whether it was right. The two `u32` per pair are 400 MB at
10000 individuals, and the `Vec` that holds them is allocated once, at
the first block, when the number of individuals is known.

## Work package 3: speed and the browser

**What it gives.** Whether the calculation reaches the numbers of
"Speed" of the dists spec, in a report, and the wheel for pyodide with
the calculation in it.

**Deliverables.**

1. The measurement. Check: `docs/reports/dists-kosman-measurement.md`
   has, for 100000 variants x 1000 individuals, the time of
   `calc_pairwise_kosman_dists` through `open_vars` on `big.vars`, on
   one thread and on 18 cores, natively, and beside each the time of
   reading `big.vars` alone through `iter_blocks` with the genotypes as
   the only field, which is taken out of the first; the same through
   `open_vcf` on `big.vcf`; the time in wasm under node through the
   TypeScript function on the bytes of `big.vars`; pyNei's time on
   `big.vars` written by pyNei's `write_vars` from the same VCF, with
   `num_threads` 1 and 6, run again that day; the load average of each
   run; and for each of the three numbers of "Speed", 0.97 s, 0.38 s and
   1.43 s, whether it is met. A number not met is a finding for the
   owner and not a task of this plan.
2. The browser. Check: `bash scripts/build_pyodide_wheel.sh && node
   tests/pyodide/smoke.mjs` exits with 0, with the smoke test extended
   to the diploid worked example through `calc_pairwise_kosman_dists`.

**What it stands on.** Work package 2.

**Tasks.**

- [x] 3.1 The measurement and its report. From "Speed" of the dists
  spec, with the machine and the load average as
  `docs/reports/filters-measurement.md` gives them, and pyNei's own vars
  file of `big.vcf` written outside the repository beside `big.vars`.
  Serves deliverable 1. Its numbers change no code; a number under the target is reported,
  and what to do about it is the owner's.
- [x] 3.2 The pyodide smoke test with the worked example, in
  `tests/pyodide/smoke.mjs`, and the build. Serves deliverable 2. Can
  run side by side with 3.1.

## How the whole plan is checked

The sum of its work packages, and the counts of "What has to be in
place" grown by the tests the plan adds: `cargo test --workspace` more
than 306, `pytest` more than 174, `npm test` more than 126, 0 failed
everywhere.
