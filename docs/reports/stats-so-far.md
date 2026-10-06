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
