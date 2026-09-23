# Report: the association study of a continuous trait

23 and 24 September 2026. It records how `docs/plans/gwas-linear.md` was
carried out, on the branch `plan/gwas-linear` in the worktree
`.claude/worktrees/gwas-linear`.

**The plan is done and the branch is not merged.** Every task is ticked,
every deliverable of the four work packages was checked by the orchestrator
running its command, and the plan's own final check passes.

## What exists now that did not

`calc_gwas` in Python and `calcGwas` in TypeScript test every variant of a
dataset against a continuous trait, giving the effect of each variant on the
trait, the uncertainty of that effect and its p-value. With covariates and
no kinship that is the linear model, checked against plink2's `--glm`. With
a kinship it is the linear mixed model, which accounts for the relatedness
of the panel so that a variant marking only ancestry does not look
associated, with either the Wald test, checked against rrBLUP, or the score
test, checked against GMMAT. Underneath are the two functions that turn a
statistic into a p-value, which every test of this plan and of
`gwas-logistic` ends in.

At the commit this branch starts from: 637 tests in the core crate, 369
pytest, 253 node. Now: 786 in the core crate with 2 ignored, the same 786 on
the faer backend, 149 in the linear algebra crate, 498 pytest on both
backends, 325 node. `cargo fmt`, `cargo clippy` with every target and
warnings denied, `cargo wasm-check` for both wasm targets and ruff are
clean, and both WebAssembly artifacts build: the npm package through
`wasm-bindgen` and the pyodide wheel.

## What is asked of the owner

**The merge.** Nothing has been merged into `main` and nothing pushed.

**Four open points of `docs/specs/gwas.md`**, none of which blocks the
merge: each has a meanwhile and the branch builds the meanwhile. Two of them
are one question over three quantities, the cancelling subtraction of "Open
2 and 3" below; the other two are the flat criterion of an unrelated panel
and the missing magnitude on the clamp. "What is open, and what the owner
decides" of work package 4 gives each with its numbers.

**One decision that is not in the spec**, because it crosses two plans.
This plan's dosage row and the row pass of `crates/popnei/src/variant.rs`
compute the same three rules twice, about 120 duplicated lines, measured by
the review of work package 2. The two agree today, over 16000 random rows,
and a test now holds them together. Whether they become one row, with the
scale made optional so that a study can ask for the dosage itself, is about
30 lines in `variant.rs` and two call sites in `pca.rs` and `kinship.rs`;
this plan could not make it, because `variant.rs` belonged to the plan
`kinship` while this one ran, and `gwas-logistic` would be the third caller
of whichever shape wins.

**Nothing is asked about the size of the module**, and this is here so that
nobody waits on a decision that is not wanted. `crates/popnei/src/gwas.rs`
is about 7200 lines and wants splitting into a directory. It was not split
because the session that runs `gwas-logistic` is held waiting on that file,
and two sessions restructuring and extending one file at once is what this
plan was told to avoid. A reviewer read the seams and they are clean: the
two distributions, the study and its design, the dosages, the result, the
linear model, the linear mixed model, and the pass over the blocks. Whoever
splits it can do so from that list without asking anything.

## The names this report uses

An association study tests each variant on its own: it fits the trait
against that variant's dosage, and reports the effect of one more copy of
the allele, the uncertainty of that effect, and the chance of seeing an
effect at least that large if the variant had none.

**The two tests.** A *Wald test* divides the fitted effect by its own
standard error and asks how extreme that ratio is; it needs a fit for every
variant. A *score test* instead asks how steeply the likelihood would rise
if the effect moved away from zero, measured at the null fit alone, so it
needs no fit per variant and is the cheaper of the two. They agree closely
where the effect is small and part where it is large.

**The mixed model and what it fits.** With a kinship, the model has a second
source of variance beside the noise: the part of the trait that goes with
how much genome a pair shares. Fitting it means splitting the trait's
variance into those two parts, and popnei does it by *restricted maximum
likelihood*, written REML below: a likelihood with the covariate effects
integrated out, so that the split is not biased by having estimated those
effects as well. What the search actually looks for is one number, the ratio
of the two variances; everything else follows from it. The search is a grid
of 101 points over that ratio's logarithm, then 60 steps of a golden section
search, which narrows a bracket by a fixed fraction each step. This is
pyNei's procedure step for step, and the plan required reproducing it.

**`y' p y`, and `num` over `den`.** `y' p y` is a quadratic form: the trait,
seen through the projection that removes the covariates and the relatedness,
multiplied back by itself. It is the generalized residual sum of squares of
the null model, and REML makes it equal to the individuals less the columns
of the design — 197 on this panel, where there are 200 individuals and a
design of three columns, the intercept and two covariates. A variant's own
test spends one more degree of freedom, for the variant, so its 196 and this
197 are different numbers and not a discrepancy. `num` and `den` are that
variant's two intermediate quantities, the numerator and denominator whose
ratio is its effect; both tests compute them identically and differ only
after.

**`sf`**, in `chi2_sf_1df` and `t_sf_two_sided`, is the survival function:
the chance of being above a value rather than below it, which is what a
p-value is. Those two functions are what every test of this plan ends in.

**The worked example** is the six-individual, three-variant case that
`docs/specs/gwas.md` writes out in full with every number, so that a test can
assert it without a reference file. It is the check with the most room in
the plan, because nothing in it is rounded for printing.

**A `Reblock`** sits over any reader and hands on blocks of a fixed number
of variants whatever sizes the reader gave. It matters below because it
means a study reads the same blocks however the file was written, which made
one check compare a computation with itself.

**The GRAMMAR-Gamma approximation** is a cheaper stand-in for the mixed
model's per-variant denominator. It belongs to the plan `gwas-logistic`,
which builds it once for both mixed models, and this plan refuses it by
name.

**The three reference programs**, each the standard tool for one of the
models, run on popnei's own data by `tests/reference/gwas/make_reference.py`:
plink2's `--glm` for the linear model, rrBLUP's `GWAS` for the mixed model's
Wald test, and GMMAT's `glmm.score` for its score test. "Checked against X"
below means popnei and X were run on the same genotypes and the same trait,
and their columns compared. pyNei is the fourth and is popnei's own oracle.

**Two linear algebra backends.** popnei runs its matrix work through
Accelerate natively and through faer in WebAssembly, so every number of this
plan was checked on both; where they differ, the two builds of popnei would
give a user different digits.

**How a tolerance is written here.** Four conventions, and every bound below
is one of them.

A bound "of `se`" is a multiple of that variant's standard error, which is
the scale of what the study estimates — an effect's own magnitude means
nothing for a variant with no effect, so nothing is bounded relative to it.
Where such a bound and a difference are given as one pair of numbers, both
are absolute and both are that variant's: "worst 3.69e-6 against a bound of
6.90e-6" means that at the variant where popnei and the reference differ
most, they differ by 3.69e-6 and were allowed 6.90e-6.

A bound "in `log10`" or "in `-log10(p)`" is the same thing said twice: the
difference between the logarithms of two p-values, which is the scale a
p-value is read on. Whether the sign is carried makes no difference to a
difference.

"The bound is 1.5e-7, breaking at 5.17e-8" means the test passes at 1.5e-7
and was then lowered until it failed, which happened at 5.17e-8. The gap
between the two is the room the check has. The spec asks for both numbers
because a bound with nothing measured against it says nothing.

A reference program's **printing floor** is how much of a difference its own
printing can account for. plink2 and GMMAT print six significant digits, so
a value is rounded by up to half a unit in the sixth — and a comparison
cannot see popnei's arithmetic below that. A check "21 times above its
floor" has room; one at its floor says popnei computes the same quantity and
no more.

## What this report is

The plan builds the parts of `docs/specs/gwas.md` that a continuous trait
needs: the two distributions, what every model shares, the linear model and
the linear mixed model. It is the second of three plans; `kinship` came
before it and `gwas-logistic` comes after. Each section below is one work
package, written as it finished.

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
-D warnings` failed in four crates. The other checks passed.

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
degrees of freedom of a variant's test on the panel were written as 197 and
195 where the spec's `n - c - 1` gives 198 and 196 (this is not the 197 of
`y' p y` above, which is `n - c`: the variant's own test spends one more,
for the variant); and the script decided whether to append its
largest draw with a float equality, which would have silently printed a
shorter array had a duplicate landed on a sampled rank.

### What the owner should know

**The plan's own warning was false.** "What could go wrong" said the two guards of the continued fraction
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
times to 1.8; the trade was seven orders of magnitude gained against that.

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

Nothing in a section with this heading bears on the merge. They are for
whoever writes or runs the next plan: what a task cost, and what went wrong
in the running of it rather than in the code. `following-plans` asks for the
token counts because they are how the right size of a task gets known.

Tasks 1.1 and 1.2 each went to one subagent; the fixes went back to the one
that wrote 1.2. Tokens: 110231, then 133005 and 228131 by the end with the
fixes. The six reviewers used 72312, 70946, 82402, 92743, 91536 and 105618.

The orchestrator's fix list told a subagent to leave `chi2_sf_1df`'s NaN
alone while the spec was settling it the other way, and the correction had
to chase it mid-task. A plan that repeats a spec's number instead of
pointing at the spec has this failure mode.

## Work package 2: what every model shares

It finished as planned, in two tasks and a round of fixes, and its four
deliverables hold. What exists now that did not: which individuals a study
tests and the design it is given, with the nine refusals that guard them;
the dosages of a block over those individuals; the choice of model and test;
and the shape of the result, with the variants that get no answer. None of
it is reachable from Python yet, which is work package 3.

| deliverable | the command | what it gave |
|---|---|---|
| 1, the tested individuals and the five refusals | `cargo test -p popnei --lib gwas::design -- --list` | 11 tests, where the starting commit printed `0 tests` |
| 2, a design whose columns are not independent is refused | the two cargo tests of that list | a covariate twice another is refused; a smallest singular value 1e-11 of the largest is kept |
| 3, the dosages are of the tested individuals | `cargo test -p popnei --lib tested_individuals_and_not_of_the_panel` | 1 passed |
| 4, the choice of model and test | `cargo test -p popnei --lib gwas::choice -- --list` | 6 tests |

The checks after the fixes: fmt and clippy clean, 668 tests in the core
crate with 2 ignored and 149 in the linear algebra crate, the same 668 with
`cargo test -p popnei --no-default-features`, `cargo wasm-check` clean, ruff
clean, 369 pytest passed.

That last check is new. `cargo test -p popnei --no-default-features` runs
the core crate on faer, the linear algebra backend the wasm build uses and
so the one that runs in a browser. It was in no check list until 23
September 2026, so the crate that holds every calculation of popnei had
never been run on it. The `kinship` session found that and put it in the
`coding` skill; it is in this plan's "What has to be in place" now, and
every check above was run on both backends. No number of this work package
differs between them, and `Design` calls the linear algebra crate for its
rank, so it could have.

### What the review found

Six reviewers ran: spec, tests, numbers, errors, api and architecture. About
twenty findings held.

**A block's variants could vanish with no error.** `of_the_block` sized its
dosage buffer from the ploidy the caller passed and read the genotype rows
using the ploidy the block stated, and nothing compared the two. When they
disagreed the parallel drive zipped a longer sequence against a shorter one,
both truncated, no row ran at all, the variant count was set to 0 and the
call returned `Ok(())`. Three reviewers found it independently with three
different reproductions; the one that reaches it is a haploid block read at
a ploidy of 2. `stats.rs` already had the guard, in `alleles_per_var_of`,
and both error cases already existed; work package 2 had simply not called
it. Every block is now checked against the reader before a row is read, and
the count of rows that came back is checked against the block's. This is the
owner's rule that an error never passes silently, and it is the finding of
this review.

**The dosages took an unchecked input.** `Design::of_the_study` runs the
refusals, but `of_the_block` took the raw `GwasInput` and asserted in a doc
comment that the positions had been checked. `Block::retain_individuals`
keeps whatever order it is asked for — there is a test of block.rs that says
so — so out-of-order positions would have measured one individual's trait
against another individual's genotypes, which is exactly what the order
refusal exists to stop. `Design` now carries the checked individuals and the
multiallelic choice, and the pass takes the `Design`. The ploidy and the
position travel in a struct of their own, where they can no longer be
swapped for each other.

**A bad covariate was reported as a defect of popnei.** A NaN or an infinity
in the design was not checked, fell through to the rank, and came back as
the error class reserved for popnei's own defects, so a user would have been
told to report a bug about their own data. pyNei raises a plain input error
on the same input. There is now a refusal naming the individual, the column
and the value, and it is in the spec.

**Two rows now compute the same rules, and a test ties them together.** The
review measured about 120 duplicated lines between this module's dosage row
and the row pass of `variant.rs`, one block of them character for character.
A reviewer ran both over 16000 random rows, at ploidies 1 to 4, with 20 per
cent of alleles missing and a third allele, in both multiallelic modes, and
found zero disagreements, so they agree today. Nothing would have failed
when one of them changed. A test now runs both on one fixture and asserts
the kept and dropped flag, the refusal with its position and count, and that
this module's dosage reconstructs the standardized one; the worst difference
measured is 0 on both backends, and replacing the missing-genotype rule with
0 makes it fail at -0.756 against -0.816.

The rest were smaller: fifteen items were public that the spec's interface
does not list, narrowed with the conditional dead-code expectation the
project already uses in two other modules, which will warn again when work
package 3 makes each one live; a `bool` parameter became an enum; a getter
that silently returned a longer slice on a broken invariant now returns an
empty one; and a doc comment of work package 1's said a refusal was "not
written yet" which task 2.1 had written an hour later.

One finding did not hold. The review asked for a test of the linear algebra error,
the one case of the fourteen with none. After the non-finite refusal was
added, that error became unreachable from a test: a design of no rows or no
columns is refused earlier by two other cases, a non-finite value is now the
new case, and what is left is a design above 2147483647 values, which is 16
GB and whose length is checked first, or a failure of the backend itself.
What it can still be reached by is written in the error's doc comment
instead.

### What the owner should know

**The deliverable the owner asked for does guard what it was written to
guard.** A reviewer moved each of the four quantities the owner named — the
major allele and so the dosages, the mean a missing genotype takes,
`allele_freq`, and whether a variant varies — to the whole panel, one at a
time, and each failed between two and four tests. It also recomputed all
eight frequency literals and every major allele by hand from the genotypes
and found each right, including that one variant counts from a different
allele over the tested four than over the whole eight. Separately, the spec
reviewer confirmed there is no other place in the module that counts
anything over the panel, which is the part a single test could not have
shown.

**popnei and pyNei agree on every dosage case the spec names.** Run side by
side on the same genotypes: half-called genotypes, a variant with nothing
called, a variant where every individual is heterozygous, a third allele
under both multiallelic modes, ploidy 1 and ploidy 4 with tied allele
counts. Identical dosages, frequencies and variance flags in every one,
including the tie-break that takes the lower-numbered allele.

**A divergence from plink2 that neither reference panel shows.**
`allele_freq` can exceed one half: the genotypes `0/. 0/. 0/. 0/. 1/1` give
1.0, in popnei and in pyNei alike. The major allele is the most frequent
among the called alleles, which counts the called half of a half-called
genotype, while the mean that becomes `allele_freq` is over whole called
genotypes, which a half-called one is not. Work package 3's first
deliverable compares `allele_freq` with plink2's `A1_FREQ`, and the spec
justified that by the two conventions agreeing. They agree on
`panel_called.vcf.gz`, where every genotype is called, and the other panel
is missing whole genotypes rather than halves, so neither panel shows it.
The spec now says "on this panel" and gives the reason; the deliverable is
unchanged, because the panel it is checked on is the one where the two
coincide.

**The duplicated row is a decision for the owner, at the end of the plan.**
The numerical reason for a separate dosage row stands, and no reviewer
contested it. The architectural question is whether the two rows should
become one, with the scale made optional so that a study can ask for the
dosage itself: about 30 lines in `variant.rs`, plus two call sites in
`pca.rs` and `kinship.rs` reading a returned mean instead of a bare flag.
This plan cannot make that change, because `variant.rs` belongs to the plan
`kinship` while this one runs, and `gwas-logistic` would be the third caller
of whichever shape wins. The cross-check test above is what holds the two
together until it is decided.

### How the work went

Tasks 2.1 and 2.2 each went to one subagent; the fixes went back to the subagent that wrote task 2.2. Tokens:
task 2.1, 182176; task 2.2, 262085, and 389180 by the end including the
fixes. The six reviewers used 115087, 136959, 158353, 141579, 145226 and
158028.

Two plan sentences were wrong and have been corrected in it, both found by
the work rather than by reading. "What it stands on" said the dosages would
use the row pass of `variant.rs`; they cannot. The replacement then cited
`ld.rs` as the precedent, which is right for a module having its own dosage
rule and wrong for the parallel drive, since `ld.rs` reads its rows one
after another and has no rayon at all. Both are the same failure as work
package 1's: a plan written before the work asserting how the code would be
shaped, in a sentence confident enough that a subagent would have followed
it.

## Work package 3: the linear model

It finished as planned, in three tasks and a round of fixes, and its five
deliverables hold. What exists now that did not: `calc_gwas` in Python and
`calcGwas` in TypeScript, testing every variant of a dataset against a
continuous trait with covariates and giving the effect of each variant, its
standard error and its p-value. This is the first work package of the plan
that gives a user anything.

| deliverable | the command | what it gave |
|---|---|---|
| 1, the whole study is plink2's | `uv run pytest tests/test_gwas.py` | over 1200 variants, worst `beta` 3.69e-6 against a bound of 6.90e-6, worst `se` 4.90e-7 against 1.65e-6, `allele_freq` exact |
| 2, the worked example and the six literals in cargo | `cargo test -p popnei --lib gwas::lm -- --list` | 7 tests, where the starting commit printed `0 tests` |
| 3, popnei and pyNei agree | the same pytest run | worst `beta` 5.83e-15 of `se` against a bound of 1.5e-14; the p-value 3.80e-12 in `log10` against 1e-11 |
| 4, the blocks change nothing | the same pytest run and a cargo test | the two studies equal to the bit; a study of 10100 variants answers the same in its second block as in its first |
| 5, TypeScript gives the same numbers | `npm run build && npm test` in `js/popnei` | 319 pass, 0 fail |

The checks after the fixes: fmt and clippy clean, 769 tests in the core
crate with 2 ignored and 149 in the linear algebra crate, the same 769 on
faer, `cargo wasm-check` clean, ruff clean, 478 pytest passed, 319 node
tests passed.

### The two deliverables that had to be rewritten before they meant anything

**Deliverable 1 could not be passed by any correct implementation.** It
asked for `beta` and `se` within 1e-5 times that variant's `se` against a
file plink2 prints to six significant digits. At `var0482`, whose `beta` is
1.0389 and whose `se` is 0.190445, the budget is 1e-5 x 0.190445 = 1.9e-6
absolute while six significant digits round a value above 1 by up to 5e-6:
the bound was smaller than the rounding of the file it compared against.
Three of the 1200 variants fail it, `var0398`, `var0482` and `var1001`, and
the worst difference is 1.938e-5 of an `se`. The bound now carries half a
unit in plink2's last printed digit beside the 1e-5 of `se`, and with that
term none of the 1200 fails.

What decides a failure is not whether `beta` passes 1 but whether half a
unit in the last printed digit passes 1e-5 times that variant's `se`, which
for a `beta` between 1 and 10 means an `se` below 0.5. Six variants meet
that and three of them fail, because the printing error is at most half a
digit and usually less.

Two counts were written into the plan on 24 September 2026 and both were
wrong: that 1198 of the 1200 failed, and that the two which passed were the
two whose `beta` passes 1. Both came from a subagent's report and were
carried into the plan and into the spec without being run. The orchestrator
ran them and they are three and six. The plan records both wrong counts rather than
replacing them quietly.

**Deliverable 4 named one test where it needed two.** Every pass puts a
`Reblock` over its reader, so a panel read from a vars file in 16 batches of
77 is joined into the single block the study reads and the study's own loop
runs once either way: the measured difference was exactly 0 because it was
the same computation, and its 1e-12 tolerance bounded nothing. It also
passed under an arithmetic mutation, the null model's degrees of freedom
substituted for the variant's, that fails the three other value tests. It
now asserts exact equality and says it shows the vars reader gives what the
VCF reader gives. The check that runs the study's loop twice is a cargo test
of 10100 variants, which the deliverable now names beside it.

### What the review found

Six reviewers ran: spec, tests, numbers, errors, api and binding, the last
of these for the first time in this plan, since this is the first work
package with a Python and a TypeScript layer. About twenty-five findings
held.

**The two builds of popnei disagreed about whether a variant has an answer.**
The variant's residual sum of squares was computed as the null model's
residual sum of squares minus `beta` times the numerator, a subtraction of
two quantities that agree to the last bits once a variant explains most of
the residual. Measured on six individuals with one covariate and a trait of
`1 + 2*cov + 3*dosage + delta*e`: at delta 1e-8 the subtraction gives
-7.1054e-15, hence an `se` of NaN, where forming the sum from the variant's
residuals gives 6.6836e-17 and an `se` of 2.36e-9; at delta 1e-7 the
subtraction gives exactly 0 and an `se` of 0, which is not a standard error
either. On the same input Accelerate gives NaN where faer gives 0.0, and
faer gives NaN where Accelerate gives a usable number. The variant came back
with a finite `beta` beside a NaN `se` and `p_value`, which is not the
all-NaN row the spec reserves for a variant with no answer, so a user
filtering on a missing `beta` would keep it. The spec now specifies the
residual form and says plainly that it departs from pyNei's arithmetic to
keep digits pyNei loses, that the two agree wherever the subtraction has not
cancelled, and that it costs one more pass over the block's dosages. No
literal moved.

**A variant in the span of the design got a large answer instead of no
answer**, with the sign depending on the backend: `beta` 5.9e13 on
Accelerate and -3.0e13 on faer, with a p-value near 1, where plink2 reports
nothing and says the correlation is too high. pyNei computes popnei's
numbers to the bit, so this was inherited from the oracle. It is the same
cancellation one level up, and the spec covers both under its second open
point with a threshold that was measured rather than judged: on an
eight-individual fixture the collinear variant leaves 6.47e-32 of its
squared length once the design is taken out and an ordinary variant leaves
0.432, thirteen orders apart, so the line is not a fine choice. The
threshold is `num_individuals` times 2.2e-16 of the squared length before,
which is the shape `docs/specs/pca.md` uses for a component with no variance
and the design's rank check uses for a covariate that is not independent, so
all three agree rather than each having its own number.

**A legal argument was refused with an untrue reason.** `test="wald"` was
refused in both layers although the Wald test is the linear model's own
test; the core accepts it and refuses only the score test of a linear model,
as the spec says. So `Error::GwasScoreTestOfALinearModel` was unreachable
from either user-facing layer, and a user passing back the `test` they were
given in the result was told it belongs to a model they were not using.

**Nothing checked that the two layers refuse the same things.** Each suite
checked itself against its own copy of the literals. That is what let
through a phenotype of strings which Python accepted and TypeScript refused,
and a `{kinship: undefined}` which TypeScript refused and Python accepted.
There is now `tests/reference/gwas/refusals_of_both_layers.json`, which both
suites walk; writing it immediately caught two messages about a covariate
named `intercept` that disagreed between the layers.

**A whole-genome study could have given every variant the first
chromosome's name.** The chromosome column is built up block by block and
nothing tested it across blocks, the only multi-block test having all 10100
variants on one chromosome. A reviewer rewrote the code so every later
variant took the first block's chromosome and the whole suite passed, while
the same mutation on the identifiers fails three tests. The multi-block test
now puts the second block on a second chromosome.

The rest were smaller: the degrees of freedom were clamped with a saturating
subtraction where the bound should be stated; two covariates of the same
name gave a result a user cannot index by name, the defect the `intercept`
refusal exists to prevent; the comparison with pyNei never ran on the panel
with missing genotypes, so the rule that fills a missing genotype with its
variant's mean was exercised at scale by nothing; the four trait, test and
model names were written once in each binding crate where the core has the
pattern for it; a dead-code expectation covered a whole block rather than
the four items that were unused, so nothing there would ever warn again; and
the TypeScript enum for which test to make was the one name in the project
matching no other layer.

One finding did not hold. The review asked that the refusal of a study with too few individuals stop
naming the file, since it is decided by the phenotype the user wrote. It is
decided by the phenotype and by the individuals the source has — the same
phenotype against a file sharing more of them is fine — so by the crate's
own rule it names the file, and a user reading a directory needs to know
which one.

### What the owner should know

**The one thing to carry into work package 4.** The same cancelling
subtraction is in the linear mixed model's Wald test, `y' p y` minus `num`
squared over `den`, and nobody has measured whether it cancels there. The
residual form is not free in that model: it would want a product with the
projection matrix for every variant, which is exactly the cost the
GRAMMAR-Gamma approximation exists to avoid. Reaching it needs a far
stronger association than the linear model's case, because the restricted
maximum likelihood fixes `y' p y` at `n - c`, which is 197 on this panel.
The spec says measure before fixing, and work package 4 will.

**A divergence from plink2 that this plan's panels cannot show.**
`allele_freq` can exceed one half, because the major allele is the most
frequent among the called alleles, which counts the called half of a
half-called genotype, while the mean that becomes `allele_freq` is over
whole called genotypes. popnei and pyNei agree exactly. One panel has every
genotype called and the other is missing whole genotypes rather than halves,
so neither shows it.

**A difference in a dtype that a user could be bitten by.** `stats["pos"]`
is an unsigned 64-bit column, as `Block.pos` and `R2Matrix.poss` are, where
pyNei's is signed. Subtracting from it wraps: `stats["pos"] - 2000` gives
18446744073709550616 for a variant at position 1000. The type is right for
popnei and the spec now lists the difference.

### How the work went

Tasks 3.1, 3.2 and 3.3 each went to one subagent; 3.2 and 3.3 ran side by
side, touching different files. The fixes went back to the subagent that wrote 3.2, and a separate
subagent built the coercion the spec settled last. Tokens: task 3.1, 262508;
task 3.3, 260240; task 3.2, 326128 and 534153 by the end with the fixes; the
coercion, 212155. The six reviewers used 140000, 184254, 171173, 190658,
153020 and 197302.

Two counts about deliverable 1 were relayed from a subagent's report to the
owner and to the session that owns the spec without being run, and both were
wrong, one of them by a factor of four hundred. The session that owns the
spec made the same shape of mistake one layer up, measuring the printing's
share on six literals and generalising it to 1200. Neither was a lapse of
care; both were trust standing in for a measurement, and in both cases what
fixed it was somebody running the thing. Four findings of this work package
came from a reviewer that ran something where reading it had missed the
point, and the plan's own sentences were wrong twice for the same reason.

## Work package 4: the linear mixed model

It finished as planned, in three tasks and a round of fixes, and its six
deliverables hold. What exists now that did not: `calc_gwas(..., kinship=k)`
with either test, in Python and in TypeScript, which accounts for the
relatedness of a panel so that a variant marking only ancestry does not look
associated.

| deliverable | the command | what it gave |
|---|---|---|
| 1, the null model is GMMAT's | `uv run pytest tests/test_gwas.py` | `genetic_variance` 1.221e-6 from GMMAT's and `residual_variance` 1.043e-6, against 1e-5 absolute, 12 per cent of the bound |
| 2, the fit is at its optimum | two cargo tests | the fitted ratio 7.72e-8 from pyNei's on Accelerate and 6.94e-8 on faer, against 2.5e-7; `y' p y` 197 within 1e-6 |
| 3, the Wald test is rrBLUP's | the same pytest run | over all 1200 variants, worst 1.9973e-5 on Accelerate and 1.9956e-5 on faer at `var0572`, against 1e-4 in `-log10(p)` |
| 4, the score test is GMMAT's, both panels | the same pytest run | `1/se²` worst 4.4268e-6 and 5.4234e-6 against 1e-5 relative; the p-value 4.6204e-5 and 4.7552e-5 against 1e-4 in `log10` |
| 5, the study finds what was planted | the same pytest run | 4 of the 5 causal variants among the 10 smallest p-values, where the plan asks for at least 3; the same panel under the plain linear model gives 1, so the check fails a model that ignores relatedness |
| 6, popnei and pyNei agree, and TypeScript | the same pytest run and `npm test` | the bound 1.5e-7, breaking at 5.173e-8 on Accelerate and 4.648e-8 on faer; 325 node tests |

The checks after the fixes: fmt and clippy clean, 786 tests in the core
crate with 2 ignored and 149 in the linear algebra crate, the same 786 on
faer, `cargo wasm-check` clean, ruff clean, 498 pytest passed on both
backends, 325 node tests passed. The plan's own final check passes and both
WebAssembly artifacts build, the npm package through `wasm-bindgen` and the
pyodide wheel.

### The deliverable that could not fail for the reason it named

Deliverable 2 said the fit is at its optimum, checked by `y' p y` coming to
197 on the panel. That is an algebraic identity: `y' p y` is the individuals
less the columns of the design for any value of the fitted variance ratio,
because the genetic variance is that same quadratic form divided by those
degrees of freedom. The review multiplied the fitted ratio by a million and
the test still gave 196.99999999999872, 1.3e-12 from 197, while every
comparison with GMMAT and rrBLUP went red.

What the identity does check is worth keeping, and its doc comment now says
it: that the clamp at 0 works, since without it the panel with the -0.0321
eigenvalue fails the Cholesky, and that the projection matrix and the
quadratic form agree with each other. What pins the search is a second test,
against pyNei's own fitted ratio.

That second test exists because four of the six things the spec says the
search must reproduce had code and nothing that would fail. The review
changed the grid from 101 points to 81, its lowest point from -10 to -10.1,
and the golden section ratio to a flat 0.61, rebuilt each time, and all 50
cargo and 50 pytest tests passed on every one. The shifted grid moves the
genetic variance by 2.85e-8 relative and the wrong ratio by 5.59e-8, both
inside the comparison with pyNei, so that did not see them either.

The bound on the new test is 2.5e-7 and not the 1e-10 the orchestrator asked
for, and the subagent refused the number with a measurement: popnei's ratio
sits 7.72e-8 from pyNei's on Accelerate and 6.94e-8 on faer, so 1e-10 fails
on code that is right. 2.5e-7 is 1.2 times the shift a wrong grid start or a
wrong ratio gives. The grid's point count alone moves the ratio by exactly
0, the minimum being well inside either grid, so that one is pinned another
way: the grid is asserted against `numpy.linspace` to the bit and the ratio
against its formula.

### What the review found

Six reviewers ran: spec, tests, numbers, errors, api and binding. Nineteen
findings held, and almost none was a wrong calculation.

**A variant could come back with an effect and no standard error, and the
two builds disagreed about which.** The Wald test's denominator is `y' p y`
minus `num²` over `den`, the same cancelling subtraction the linear model
had. The spec had recorded a measurement on the panel — 8.53e-13 on
Accelerate against 3.98e-13 on faer, both positive — and concluded that no
NaN appeared. A reviewer found one, and the orchestrator reproduced it: six
individuals, one covariate, an identity kinship and a trait built as
`2 + 3*cov + 1*dosage` give `beta` 1.00000 with `se` and `p_value` NaN under
the Wald test, where the score test on the same call with one argument
changed gives `se` 0.500 and `p` 0.0455. A finite `beta` beside a NaN `se`
is a fourth kind of NaN that "The variants that have no answer" does not
describe, so a user filtering on a missing effect keeps the row and reads
1.0 as an effect that was measured — which is the argument that settled the
linear model's subtraction one model earlier.

The spec's two open points about cancellation became one, over three
quantities, with the same remedy: refuse when what is left falls to
`num_individuals` times 2.2e-16 of what there was. The meanwhile changed
from leaving the Wald test alone to refusing there too, because a meanwhile
that returns NaN where the score test returns 0.0455 on the same call is not
a safe thing to build on. That is what the branch does now.

**A trait the design explains exactly was a defect of popnei.** The REML
left both variances at 0, so the covariance was 0, the division made
infinities and the linear algebra refused them: the user got a
`RuntimeError` naming an operand of a matrix routine and the VCF, neither of
which was at fault, for what is a wrong input. It is a `ValueError` now,
saying the design explains the whole trait. Note that the linear model gives
an `se` of 2.4e-15 on the same input rather than an error, which is the same
threshold question one model over and was not touched.

**The kinship was copied three times where the binding said once.**
Measured at 3000 individuals with a 72 MB matrix: 144.4 MB peak, 72.5 MB
held. Indexing the array once instead of going through pandas gives 72.2 MB
peak and comes back already in the right memory order, so the copy that
existed to fix the order is gone, and so is the binding's own copy — it
lends the array now. At the largest dataset `docs/objectives.md` names,
10000 individuals, that is 800 MB of transient allocation saved.

**Nothing ran the mixed model over more than one block.** 200 individuals
give 10000 variants to a block, so the 1200-variant panel is one block in
every test. Removing the three buffer clears from the mixed model's block
pass left all 50 cargo and 493 pytest tests passing, where the identical
mutation in the linear model fails its multi-block test at once. The mixed
model has that test now.

**The two layers disagreed again, in three more places**, none of them in
the shared refusals file: an infinite phenotype, where the messages differed;
a NaN covariate, "is missing, fill the value in" against "is NaN, and a
study is fitted on numbers"; and a full-width digit, which Python reads as a
number and TypeScript refused, because the pattern matching digits was the
ASCII one. That last contradicts the settled rule that a value is refused
exactly where Python's `float` raises. All three are in the file now.

The rest were smaller: two scalars of the same kind meeting unnamed in one
signature, where swapping them compiles and gives a wrong inverse; a
function that bundled what the next plan has to split, now split while it
has one caller; the kinship trusted as its constructor left it, with the
symmetry check moved into the core so that both layers and a caller of the
core are covered; three untested paths through the binding; an error naming
the VCF path where its twin did not; Python coercing a boolean where
TypeScript refused one; and half a dozen comments whose numbers were wrong.

Two things the subagent refused with evidence, and both refusals hold. The
orchestrator's proposed 1e-10 bound on the fitted ratio fails on correct
code, as above. And a reviewer reported that `Kinship.principal_components`
re-checks the matrix, so `calc_gwas` should too; it does not, it only makes
the array contiguous, so the check went into the core instead, which covers
more.

### What is open, and what the owner decides

Four open points of `docs/specs/gwas.md` reach the owner from this work
package, and none blocks the merge: each has a meanwhile, and the branch
builds the meanwhile.

**Open 2 and 3, now one: the cancelling subtraction, in three places.** The
choice is to refuse such a variant, which is one comparison in each place
and costs nothing per variant, or to form the residual exactly, which for
the linear model is a pass over the block's dosages and is already taken,
and for the mixed model's Wald test is a matrix product for every variant,
roughly doubling that test — the cost the GRAMMAR-Gamma approximation of the
next plan exists to avoid. The spec recommends refusing and the branch
refuses. The exact form buys only the band just above the threshold, which
needs a variant explaining about 99.9999 per cent of the trait, and that is
for the performance session to weigh with numbers in front of it.

**The flat criterion.** When the kinship is close to a multiple of the
identity — a panel of unrelated individuals, or a user passing an identity
matrix to mean no relatedness — the REML criterion is flat to within one
unit in the last place across all 101 grid points, so which point wins is
rounding. A reviewer perturbed such a kinship by 1e-15 and got a
`heritability` of 6.5e-5, 7.1e-5 and 0.967 over three seeds. A kinship of
all zeros gives 0.99988 and one of all ones 0.99995. `beta` and `p_value`
are untouched, because they do not depend on the scale, so what is arbitrary
is `genetic_variance`, `residual_variance` and `heritability` — which is
what a user reads a heritability off.

The branch builds the spec's meanwhile, which is also its recommendation:
the study is still given, every variant still answered, and those three
fields come back as nothing. In Python they are `None` and in TypeScript
`undefined`, which is what the plain linear model already gives for the two
it does not have, so a user meets one shape and not a new one; none of the
three is 0 and none is NaN. On the worked example over an identity kinship
that is now

    genetic_variance: None  residual_variance: None  heritability: None

where the branch gave a `heritability` of 6.8e-5 before, and `beta` is still
1.5 and the p-value still 0.0880048923827560. The panel's own kinship is not
flat and reports all three as it did. What is left for the owner is whether
to refuse the study outright instead; the argument for giving it is that the
effects and the p-values are valid and are what the user mostly came for.

**The clamp's missing magnitude.** A kinship can have slightly negative
eigenvalues from rounding alone, and popnei clamps them to 0. There is no
test of how negative, so a matrix that is genuinely not a covariance is
treated as a rounding artefact. The question is what fraction of the largest
eigenvalue still counts as rounding, and the two ends of it, put in the same
terms:

| case | smallest eigenvalue, as a share of the largest | should be |
|---|---|---|
| this plan's panel, every genotype called | -2.0e-16 | clamped |
| a kinship of a panel with 50 genotypes missing in 100, from `docs/specs/kinship.md` | -0.033 | clamped |
| a reviewer's matrix with one eigenvalue forced to -5 | -0.29 | refused |

So the line lies between 3.3 per cent and 29 per cent, and the safe
direction is to refuse: a matrix that far from a covariance is not a
kinship, and clamping it returns ordinary-looking numbers with nothing to
say they are wrong.

What nobody has measured is how far a real dataset goes, and that is the one
thing this decision wants that the report cannot give. The measurement is
cheap and is not this plan's: take the panels `docs/specs/kinship.md`
already has, raise the missing genotypes until the smallest eigenvalue stops
falling, and see where it settles. Until somebody runs it, any line between
the two is a guess, and the branch's own behaviour — clamp everything — is
the most permissive guess available.

### What the owner should know

**One check of this work package is weaker than it looks, and one is
stronger than the orchestrator first said.** GMMAT prints its score test to
six significant digits, so at the worst variant of deliverable 4 at most
1.35e-7 of the measured 4.4268e-6 is popnei's; the rest is the printing. It
is not hollow — the bound still fails if popnei's variance is about 5.7e-6
relative off, forty times above what the comparison can resolve, and 282 of
the 1200 variants exceed their own printing floor. The p-value half is
genuinely tight, 21 times above its floor, and rrBLUP's file is full
precision and carries the arithmetic.

**The refusal threshold is guarded on one backend only, and cannot be
guarded on both.** With the threshold set to 0, so that only a non-positive
denominator refuses, Accelerate passes everything and faer fails: the
collinear fixture's denominator is -4.44e-16 on one and +1.665e-15 on the
other. What is left is exactly 0 in exact arithmetic, so its sign is
whatever the rounding chose, and no fixture can make it positive on both.
faer is the guard and the fixture records it.

**A regime the spec asks for and no fixture reaches.** The best point of the
grid is 44, 46, 47 or 85 in every test of both suites; the spec lists the
clamping of the best point's neighbours at the ends of the grid among the
six things that must reproduce pyNei, and nothing is in that regime. A trait
that is almost entirely genetic or almost entirely noise would reach it.

### How the work went

Tasks 4.1, 4.2 and 4.3 each went to one subagent in turn, and a fourth
closed the gap that opened when the spec decided a constant trait is
refused. The fixes went back to the subagent that wrote 4.3. Tokens: task
4.1, 239120; task 4.2, 266521; task 4.3, 357304 and 532368 by the end with
the fixes; the constant trait, 124448. The six reviewers used 209298,
179770, 135493, 177308, 191305 and 187393.

Task 4.1 was asked to measure the cancellation while it had the null fit in
hand, two tasks before anything would act on it, and told not to act on it
itself. It found that the cancellation is reached exactly rather than
approached, because the projection annihilates the design, so any affine
image of the trait gives equality in exact arithmetic. That is what let a
later reviewer construct the case that produced the NaN, which is what
changed the spec. Asking for a measurement early, from whoever has the
pieces in hand, and separating it from the decision that will use it.
