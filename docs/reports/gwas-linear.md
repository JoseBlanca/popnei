# Report: the association study of a continuous trait

Started 23 September 2026. It records how `docs/plans/gwas-linear.md` was
carried out, on the branch `plan/gwas-linear` in the worktree
`.claude/worktrees/gwas-linear`. The plan builds the parts of
`docs/specs/gwas.md` that a continuous trait needs: the two distributions,
what every model shares, the linear model and the linear mixed model. The
plan is not finished, and this report grows as each work package does.

## Where the branch starts, and why not where the plan says

The plan says the branch starts from `main` with the plan `kinship` merged
into it. `kinship` is not on `main`: another session is running it and was
on its work package 3 when this one began. The owner decided on 23 September
2026 to start from `plan/kinship` as it stood instead, because its work
package 1, which moved into `crates/popnei/src/variant.rs` the row pass that
turns the genotypes of a variant into standardized dosages, and its work
package 2, the kinship matrix with its Python and TypeScript functions, are
the two things this plan needs from it, and both were done and reviewed.

The tip of `plan/kinship` did not compile. At `8392778` the core crate had
already changed `Kinship.num_vars` from a `usize` to a `u64` and neither
binding crate had followed: `crates/popnei-js/src/kinship.rs:153` failed
with `expected usize, found u64`, `crates/popnei-python/src/kinship.rs:102`
failed clippy's `useless_conversion` on a `u64::try_from` that had become a
conversion to its own type, and `crates/popnei/src/kinship.rs:944` failed
`chunks_exact_to_as_chunks`. So `cargo clippy --workspace --all-targets --
-D warnings` failed in four crates. The other checks passed. This is a
branch caught part way through a review fix, not a defect of it.

The branch therefore starts from `b387def`, "the kinship spec: the bits of
plink2, the variants a pass gave, and where the limit of the individuals
lives", the last commit before that change, which the owner approved on 23
September 2026. Every check of the `coding` skill passes on it:

| check | result |
|---|---|
| `cargo fmt --all --check` | clean |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo test --workspace` | 637 passed in the core crate, 2 ignored; 149 in the linear algebra crate |
| `cargo wasm-check` | clean |
| `cargo test -p popnei --lib gwas -- --list` | `0 tests` |
| `uv run ruff format --check && uv run ruff check` | 29 files formatted, all checks passed |
| `uv run maturin develop && uv run pytest` | 369 passed |
| `npm run build && npm test` in `js/popnei` | 253 pass, 0 fail |

The plan's "What has to be in place" gives 604, 347 and 242 for the three
suites, measured before any kinship work existed, and says the counts would
grow and not shrink. They did.

What `b387def` does not have is the tail of the kinship review: `num_vars`
as a `u64`, a new `num_vars_given` field holding how many variants the
reader gave whether they were used or not, and the constant
`MAX_INDIVIDUALS_OF_THE_VARIANTS` moving from `pca.rs` to `variant.rs`. All
three are in `Kinship`, which work packages 1, 2 and 3 do not touch. Work
package 4 is the first that does, and `plan/kinship` is merged in again
before it. That is what the owner decided and it is what the session running
`kinship` recommends, since more of its review was still landing.

The six items of the row pass this plan stands on are unchanged at
`b387def` and the `kinship` session has confirmed they are not moving:
`the_standardized_row`, `the_standardized_rows` and `the_standardized_block`
with their signatures, `DosageOptions` with its `transform_to_biallelic` and
its `DosageScale`, `RowPositions` with its `first` and its `too_many`, and
`DosageScale` itself. One thing inside `the_standardized_row` changed with
its signature unchanged: it raises `Error::VariantPloidyTooLarge` where it
raised `Error::PcaVariantsTooLarge`, the old name having said the principal
components about a variant that has nothing to do with them. No task of this
plan matches on that case.

## What was changed in the plan, and why

Three changes, all on 23 September 2026, before the first task.

The state line says under way, and the paragraph that gives the worktree and
the branch now says where the branch really starts and points here for the
reason.

"What has to be in place" asks for work packages 1 and 2 of `kinship`
instead of `kinship` merged into `main`, and carries the counts measured at
`b387def` beside the ones the plan was written with.

Deliverable 1 of work package 3, the comparison with plink2, no longer
repeats a tolerance. The plan had `beta` and `se` within 1e-5 absolute and
`p_value` within 1e-5 relative; `docs/specs/gwas.md` now asks `allele_freq`
within 1e-6 absolute and `beta`, `se` and `p_value` within 1e-5 relative,
and the deliverable points at the spec. The session that wrote the spec
changed it after the review of `kinship`'s work package 2 found that its own
plink2 comparison, 1e-5 absolute against a file printed to six significant
digits, was measuring plink2's rounding and not popnei: against plink2's
binary output popnei agreed to 4.44e-16, so the check had about twofold
headroom and would have passed an arithmetic error of up to 5e-6. Six
significant digits round a value by up to 5e-6 of itself, which is a
relative amount, so an absolute tolerance holds on this panel, whose values
are small, and breaks on data whose values are larger for code that is
right. Two of this plan's five reference files are printed that way,
plink2's `--glm` and GMMAT's `glmm.score`, and neither program has a binary
form for them. The spec change is commits `57f749d` and `1fe133c` of
`spec/gwas`, merged into this branch as `d285e83`.

What the owner should take from it: a work package of this plan that passes
its comparison with plink2 or with GMMAT has shown that popnei computes the
same quantity, not that its digits are right. The checks with headroom are
the worked example at 1e-12 relative, which needs no reference file, and the
comparison with pyNei at 1e-9 relative over every column. Both are
deliverables of work package 3, and each work package below says which of
its deliverables could have caught a wrong digit.
