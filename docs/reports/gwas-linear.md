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

## Work package 1: the two distributions

It finished as planned, in two tasks and then a round of fixes, and its
three deliverables hold. What exists now that did not: `chi2_sf_1df` and
`t_sf_two_sided` in `crates/popnei/src/gwas.rs`, public in the core crate
and in neither package, which every p-value of this plan and of
`gwas-logistic` ends in.

| deliverable | the command | what it gave |
|---|---|---|
| 1, `libm` builds for both wasm targets | `cargo wasm-check` | clean, both targets, warnings denied |
| 2, the two functions are scipy's | `cargo test -p popnei --lib distributions -- --list` | the three tests the spec names, where the starting commit printed `0 tests` |
| 3, the literals came from scipy and can be got again | `tests/reference/gwas/print_scipy_distributions.py` under scipy 1.18.1 and numpy 2.5.3 | all 193 printed lines present in `gwas.rs` verbatim, 156 data lines on each side, none in the file the script did not print |

The checks after the fixes: fmt and clippy clean, 640 tests in the core
crate with 2 ignored and 149 in the linear algebra crate, `cargo wasm-check`
clean, ruff clean, 369 pytest passed, 253 node tests passed.

Deliverable 1 was checked beyond what the plan asks. `cargo wasm-check` runs
`cargo check`, which type-checks without linking, so "both wasm targets
build" rested on less than it claimed. The `coding` skill asks that a change
touching a dependency also build the wasm wheel and the TypeScript package,
and task 1.1 had built neither. Both were built: the npm package links its
wasm through `wasm-bindgen` and passes its tests, and
`scripts/build_pyodide_wheel.sh` produced
`popnei-0.1.0-cp314-cp314-pyemscripten_2026_0_wasm32.whl`. The risk was nil,
because `libm` 0.2.16 was already in both wasm dependency graphs through
faer and gemm and the direct dependency added one line to `Cargo.lock`, but
the deliverable now says links and not type-checks.

### What the review found

Six reviewers ran: spec, tests, numbers, errors, api and architecture. Two
findings changed numbers.

`t_sf_two_sided` lost up to eight digits for a `t` near 0. It computed
`x = df / (df + t²)` and then recovered `1 - x` by subtraction, from an `x`
that had already rounded to 1. At 197 degrees of freedom and `t` of 1e-7 it
returned exactly 1.0 where the true value, from mpmath at 60 digits, is
0.9999999203127337, a relative error of 8.0e-8 against a spec tolerance of
1e-10. It now carries `t² / (df + t²)` and never subtracts: the worst
relative error over `t` in [1e-7, 1e-3] fell from 7.97e-8 to 5.34e-17 at 197
degrees of freedom and from 6.34e-7 to 3.83e-15 at 9997. The wrong values
were all p-values above 0.99994, so no association would have been called
differently, but they were wrong.

A degrees of freedom of 0 or less returned a p-value instead of refusing:
`t_sf_two_sided(1.0, 0.0)` gave 0.0, the most significant p-value there is,
and a negative one gave 1.0. Four of the six reviewers found it
independently. It now gives NaN, which also covers a NaN degrees of freedom.
Nothing reaches it today; what will keep the degrees of freedom at 1 or
above is work package 2's refusal of a design with no more rows than columns
plus one, which is not written yet. pyNei has the same hole, so this is not
a divergence from the oracle but a silent wrong number both had.

Three findings were about the tests rather than the code. The incomplete
beta's tolerance of 1e-12 absolute verified nothing at 8 of its 51 literals,
which are at or below 1e-12: a reviewer showed that an implementation
returning 0 in the tail passed the test, and that only the Student t test
caught it, and only for one of the four pairs. A relative assertion now runs
beside the absolute one, at 1e-12 relative, whose worst measured differences
at the four pairs are 1.35e-15, 2.05e-14, 6.60e-14 and 1.19e-15, so 15 times
the room where it has least. Every Student t fixture was at 5, 17 or 197
degrees of freedom while the error grows with them, so 997 and 9997 were
added, 9997 being the 10000 individuals of `docs/objectives.md` less the
three columns of a design and the variant. And the behaviour the doc
comments promise off the fixture range, a NaN in giving a NaN out, had no
assertion; it has one now.

The remaining findings were documentation: the `libm` justification named
`erfc` and not `lgamma`, which the code also calls; the crate doc of
`lib.rs` names every public module and had no clause for `gwas`; the module
doc never expanded `sf` and gave the chi square to the score test alone; the
panel's degrees of freedom were written as 197 and 195 where the spec's
`n - c - 1` gives 198 and 196; and the script decided whether to append its
largest draw with a float equality, which would have silently printed a
shorter array had a duplicate landed on a sampled rank.

### What the owner should know

**The plan's own warning was false, and this is the useful thing the review
found.** "What could go wrong" said the two guards of the continued fraction
are what make it converge and that dropping either gives numbers right for
most arguments and wrong for some. Neither is true for any argument either
plan can reach. The first denominator is bounded below by `2 / (a + b + 2)`
in both branches, so with `b` of 1 / 2 the `tiny` guard can fire only once
`a + b` passes about 2e300, a panel of 4e300 individuals; one reviewer
measured a minimum of 4.0276e-6 over 6009003 calls, which is that bound
exactly at a million degrees of freedom, and another 4.06e-5 over 116802.
The `eps` caps the work and not the digits: with it turned off the fraction
runs its 500 rounds, nothing becomes non-finite and the worst value moves by
2.3e-13 relative. So no test on a value can catch either guard's removal,
and a test writer who believed the plan would have hunted for a case that
does not exist. Both guards stay, because Numerical Recipes and pyNei have
them and because a caller with another `b` would need the `tiny`. The plan
now says this, and so does the spec.

**One bound of the spec sits close to its failure.** `t_sf_two_sided`'s
error grows with the degrees of freedom, and the worst point is at `t` near
1.73, where the calculation switches branches. Measured against mpmath at 60
digits after the fix: 4.7e-13 at 197 degrees of freedom, 7.7e-13 at 997,
5.5e-11 at 9997, 1.51e-10 at 20000 and 3.48e-9 at 500000. The spec's
tolerance is 1e-10, so it holds at 9997, which is the largest panel
`docs/objectives.md` names, with 1.8 times the room, against 512 times at
the degrees of freedom the other cases use. The cause is the
`lgamma(a + b) - lgamma(a)` cancellation in the front factor, amplified by
the `1 -` of the symmetry branch; it is not the stopping rule, which a
reviewer confirmed by setting it to 0. Whoever first wants popnei past 10000
individuals has to come back to this function before trusting its p-values,
and the front factor in logarithms is where to start. Fixing the small `t`
defect cost about thirteen per cent here, taking the room at 9997 from 2.4
times to 1.8; the trade was seven orders of magnitude gained against that,
and it is worth naming because it is the one place in this work package
where making one number better made another worse.

**The spec's written recipe is now behind the code.** The spec writes the
symmetry branch as `cf(b, a, 1 - x)` and the front factor's last term as
`b * ln_1p(-x)`. The code computes neither from `x` any more, since not
subtracting is the whole of the fix. The spec is the transcription the next
implementer copies, so copying it reproduces the defect. The session that
owns the spec has been told; it changes no value a user sees.

**One question is the owner's and is open.** `chi2_sf_1df` gave NaN for a
negative argument where scipy gives 1.0. The spec has since settled the
function half: an `x` of 0 or below gives 1.0, and the code now does. What
is still open, as **Open 2** of `docs/specs/gwas.md`, is who may pass one:
the score test's denominator is a quadratic form that can round just below 0
for a variant with almost no variance, and the choice is whether such a
variant is refused with three NaNs or gets a p-value of 1. The spec's
meanwhile is to refuse, so work package 4 is not blocked, and no literal of
either panel reaches the case either way.

### How the work went

Task 1.1 and task 1.2 each went to one subagent and each came back right the
first time; the fixes went back to the subagent that wrote task 1.2, which
had the context. Tokens: task 1.1, 110231; task 1.2, 133005, and 228131 by
the end including the fixes. The six reviewers used 72312, 70946, 82402,
92743, 91536 and 105618.

One instruction of the orchestrator's was wrong and had to be corrected
mid-flight: the fix list told the subagent to leave `chi2_sf_1df`'s NaN
alone, and the spec settled it the other way while the subagent was working.
Sending the correction cost nothing because the subagent was still running,
but a plan that carries a spec's number instead of pointing at the spec has
this failure mode, which is the same one that the tolerance change earlier
in this report was about.
