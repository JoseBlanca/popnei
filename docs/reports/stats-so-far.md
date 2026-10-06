# Report: the result so far of three statistics, and the three in one pass

7 October 2026. The work report of `docs/plans/stats-so-far.md`, on the
branch `spec/stats-so-far`. State: under way.

## 1. The core

Task 1.1 is 526ad30: `calc_per_var_distribs_with`,
`calc_per_individual_stats_with` and `calc_var_density_with`, with `SoFar`
and `AfterABlock`; the three functions of before call them with a function
that does nothing. `cargo test -p popnei --lib so_far` runs 6 tests, each
result so far compared to the bit with the calculation over the first that
many variants. pytest gives 751 passed with no test of Python touched.
`calc_per_var_distribs` on one thread, three rounds of five runs of each
build interleaved, gave medians of 0.116, 0.117 and 0.116 s before and
0.116, 0.116 and 0.116 s after on `big.vars`, and 0.646, 0.666 and 0.651 s
before and 0.670, 0.652 and 0.643 s after on `big.vcf`.

`many.vcf` has no `##contig` length, so the density of the TypeScript test
"with the lengths of the file" would have tested the windows that grow
twice; the spec gives the lengths in that test now.

Task 1.2 is 660cf55: `calc_variants_summary` in
`crates/popnei/src/stats/summary.rs`, with the totals of each of the three
calculations in a type of their own that the three and the summary share,
so that no counting is copied, and the error of a summary of none.
`cargo test -p popnei --lib variants_summary` runs 5 tests, among them the
seven combinations of the three against the three calculations, and a
reader that records what it is asked for.

Task 1.3 is d84d674: the bench `crates/popnei/benches/variants_summary.rs`,
whose numbers are in "The three statistics of a file in one pass" of
`docs/specs/js_sources.md`. Medians of 7 interleaved runs, one rayon
thread, a release build, on the owner's Apple M5 Pro with other work on it
(a load of about 15), from the opening of the file to the result:

| | the three passes | the one pass | saving |
|---|---|---|---|
| `big.vcf` | 2.240 s | 1.247 s | 44% |
| `big.vars` | 0.400 s | 0.391 s | 2% |

The architecture reviewer ran it again with 5 runs under a load of 16 to
23: 1.529 s and 0.830 s on `big.vcf`, 46%, and 0.252 s and 0.227 s on
`big.vars`, 10%. The saving on the VCF holds; on the vars file it is small
and within what the load of the machine moves. It is all in the read: a
pass that parses genotypes waits about 0.91 s for them on the VCF, which
the one pass does once, while the counting of each statistic is not
shared, and a vars file is read in under 0.01 s. The plan's condition to
stop, a saving under a fifth on both files, was not met.

### The deliverables

Run on b0e36eb, after the fixes of the review.

| deliverable | command | result |
|---|---|---|
| 1, the result so far | `cargo test -p popnei --lib so_far -- --list` | 7; 0 on d22e7a2 |
| 2, the summary | `cargo test -p popnei --lib variants_summary -- --list` | 11 after the fixes, 5 before; 0 on d22e7a2 |
| 3, nothing of before changed | the six cargo commands and `uv run pytest` | 1506 and 1356 passed; pytest 751, with no file of Python or of its tests changed |
| 4, the measurement | the bench | above |

### What the review found

Five reviewers, spec, tests, errors, api and architecture, each in a
worktree of its own, read 8d1e215..4c8662d. None found a wrong result.
These held and are fixed:

- Three reviewers found that no test made the summary fail in the middle of
  a pass: mutants that dropped the error of any of the three statistics,
  or of the function after a block, passed every test. Four tests now catch
  each of them (f7b0057).
- The density took a block of no variants and called the function with a
  result over none, where the other two refuse such a block as a defect of
  the reader; it refuses it now (a95f42b).
- The options of the density in the summary were a pair of numbers that a
  call could not read; they are `VarDensityConfig` now (ab3b30e), and the
  function that does nothing is public (ff0ee3d).
- The final result copied what it could move, 40 MB for the widest density
  (a202193), and five doc comments said less or more than the code (b0e36eb,
  and the spec in 26021d8).

## 2. The wasm binding crate and the TypeScript package

Task 2.1 is 96e9ac6: `onSoFar` and `soFarEvery` on the three calculations.
The binding reads the clock of JavaScript, calls the function after a block
with the result built by the code that builds the final one, and ends a
pass whose function threw as one that `onProgress` stopped, telling the page
nothing more. 42 node tests in `test/so_far.test.ts`; 41 of them failed
before the change, and the one that did not, a `soFarEvery` of 3600 that
never calls, passes on a function that takes no options as well.
