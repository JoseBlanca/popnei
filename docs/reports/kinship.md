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
