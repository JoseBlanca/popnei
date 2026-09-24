# Report: how much variety each population holds

24 September 2026. It records how `docs/plans/diversity.md` was carried
out, on the branch `plan/diversity` in the worktree
`.claude/worktrees/diversity`, which branches from `spec/diversity` and
not from `main`, because neither the spec nor the plan is on `main` yet.

**The plan is under way.** This page is written while the work goes, so
what follows grows one work package at a time.

## Before the first task

The plan's "What has to be in place" was checked by running it, on the
owner's Apple M5 Pro under macOS 27.0, on 24 September 2026, at commit
664485e.

| what | command | what it gave |
|---|---|---|
| the format | `cargo fmt --all --check` | exit 0 |
| the core and the linalg crate | `cargo test --workspace` | `787 passed`, 2 ignored, and `149 passed` |
| the lints | `cargo clippy --workspace --all-targets -- -D warnings` | no warning |
| the two wasm targets | `cargo wasm-check` | exit 0 |
| the faer backend | `cargo test -p popnei --no-default-features` | `787 passed`, 2 ignored |
| the Python lint | `uv run ruff format --check` and `uv run ruff check` | `33 files already formatted`, `All checks passed!` |
| the Python suite | `uv run maturin develop && uv run pytest` | `499 passed` |
| R and its three packages | `Rscript -e 'packageVersion(...)'` | R 4.6.1, `vegan` 2.7.6, `adegenet` 2.1.11, `poppr` 2.9.8 |
| the panel | `ls tests/reference/stats/panel.vcf.gz` | there, with `panel_pops_bcftools.txt` |
| the file of the measurement | `ls /Users/jose/devel/popnei-bench/big.vars` | there, so task 4.1 does not write it |
| the TypeScript package | `npm run build` then `npm test` in `js/popnei` | exit 0, then `tests 325`, `pass 325`, `fail 0` |

Every count matches the one the plan wrote down, so the "more than" of
the deliverables is counted from these.

`npm run build` was run once at the start, as the plan asks, before any
task. Without it `npm test` fails in a fresh worktree on a file nobody
touched.

**`dadi` installs reproducibly**, which is the one thing that could have
sent the spec back to the owner. `uv venv --python 3.12` and `uv pip
install 'dadi==2.4.4'` succeeded, and the `dadi` they installed folded a
projected spectrum of a made up variant, 0.1333333333, 0.5333333333,
0.3333333333 and 0 in a draw of 6 called alleles. So the spec's check of
the folded spectrum stands as written and task 1.2 can be built on it.

The checks the plan says fail today do fail: `cargo test -p popnei --lib
-- diversity:: --list` prints `0 tests`, and `tests/test_diversity.py`,
`js/popnei/test/diversity.test.ts` and `tests/reference/diversity/` do
not exist.
