# Report: the association study of a binomial trait

24 September 2026. It records how `docs/plans/gwas-logistic.md` is being
carried out, on the branch `plan/gwas-logistic` in the worktree
`.claude/worktrees/gwas-logistic`.

**The plan is under way and nothing is asked of the owner yet.** No work
package is finished. This report is written as the work goes, so what is
below is what has happened and not what is planned.

## Where the branch starts

The branch starts from `main` at ecc1373 with `plan/gwas-linear` merged into
it, which is the merge commit a315d43. That branch was not on `main` when
this work began and still is not: the association study of a continuous
trait waits on the owner to merge it, so this plan carries it rather than
waiting for it.

Every check of the `coding` skill was run on a315d43 before the first task,
and all nine pass:

| check | result |
|---|---|
| `cargo fmt --all --check` | clean |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo test --workspace` | 787 core tests, 2 ignored; 149 in the linear algebra crate |
| `cargo test -p popnei --no-default-features` | the same 787, 2 ignored, on faer |
| `cargo test -p popnei-linalg --no-default-features` | 136 |
| `cargo wasm-check` | clean on both wasm targets |
| `uv run ruff format --check && uv run ruff check` | clean |
| `uv run maturin develop && uv run pytest` | 499 passed in 21.03 s |
| `npm run build && npm test` in `js/popnei` | 325 passed, 0 failed |

The faer run is the one that matters most to this plan, because it is the
linear algebra a browser runs and because every tolerance below has to hold
on it as well as on Accelerate.

The linear algebra crate has 149 tests with its default features and 136
without: 13 of them are of the system BLAS and LAPACK, which the faer build
does not have. The two counts are not a loss.

## The module was split before anything was added to it

`crates/popnei/src/gwas.rs` had reached 8035 lines and this plan adds two
models to it, so it is now the directory `crates/popnei/src/gwas/`, seven
modules and a `mod.rs`, at 3f77fa8. The seams are the ones a reviewer of
`plan/gwas-linear` read and reported as clean: `distributions` for the two
functions that turn a statistic into a p-value, `study` for what a user asks
for and what a study is refused for, `dosages` for the dosage of each tested
individual at each variant of a block, `result` for what a study gives back,
`linear` and `linear_mixed` for the two models that exist, and `pass` for
`calc_gwas` and the one pass over the blocks. A logistic model is a module
beside the other two.

Nothing but the address of the code changed, and three things say so: no
file outside `crates/popnei/src/gwas` is in the commit, so no caller of the
crate had to be adjusted; every line of the old file is in the new directory
unchanged but for the import lists, sixteen items widened from private to
`pub(super)` or `pub(crate)` because they now cross a module boundary, and
the doc links that a submodule can no longer resolve by a bare name; and
every check gives the count it gave before, 787 core tests with 2 ignored on
both backends, 149 and 136 in the linear algebra crate, 499 pytest and 325
node.

The plan's tasks name `crates/popnei/src/gwas.rs` as where their code goes.
That file no longer exists and the tasks build in the directory instead.

## What the owner's four open points would change

The four open points of `docs/specs/gwas.md` are the owner's and unanswered.
Each has a "meanwhile" and the work follows it. This list says which tests
would move if an answer differs from its meanwhile, and it is kept as the
work goes rather than written at the end.

**Open 1, a variant that separates the cases from the controls.** Meanwhile:
`beta`, `se` and `p_value` are NaN with no reason given, which is pyNei's
behaviour. Under "NaN with a reason", the owner's own recommendation, no
number moves: a field is added to the result and a column to both packages,
and the tests that would change are the ones that read the shape of the
result rather than its numbers. Under the Firth penalized regression, the
comparison with plink2 gains `var0006`, the one variant of 1200 it now
leaves out, and deliverable 3 of work package 1, which asserts that popnei's
NaNs are exactly the rows plink2 marked `FIRTH?` `Y`, stops being true and
has to be replaced.

**Open 2, a variant there is nothing left to test.** Meanwhile: refuse, when
what is left falls to the tested individuals times 2.2e-16 of what there
was. Two of the four places this rule holds in are tests of this plan, the
logistic model's score test and the logistic mixed model's, and both are
written against that threshold. Under "form the residual exactly", the
threshold goes and those two tests are replaced by the exact quantity; no
literal of either panel is near the threshold, so the six literals do not
move either way.

Two of the four places have a dataset behind them and two do not. The linear
model's was reached on a fixture of eight individuals with a variant fixed
one way in each of two subpopulations, and the linear mixed model's Wald
test on six individuals with an identity kinship; both were reproduced
before the rule was written. **Neither score test has been reached on any
data.** The case for guarding them is that the arithmetic is the same
cancellation in another denominator, which is an argument and not a
reproduction. So each of the two tasks that builds one also builds a fixture
that reaches it, and this report records what that fixture leaves against
the threshold it is judged by, the way the linear model's test records
6.47e-32 against a threshold of 1.78e-15. A guard whose case has never been
seen is a guard that can be wrong in the direction of never firing.

**Open 3, a kinship that does not identify the two variances.** Meanwhile:
give the study with `genetic_variance`, `residual_variance` and
`heritability` set to `None`. What it does to the logistic mixed model is
not yet measured and will be recorded here when work package 2 is built.

**Open 4, how negative an eigenvalue is still rounding.** Meanwhile: refuse
when the smallest eigenvalue of the kinship is below a tenth of the largest.
This is already built by `plan/gwas-linear` and the logistic mixed model
inherits it rather than deciding it again, so an answer that changes the
fraction changes one constant and the tests of that refusal, which are not
tests of this plan.

## What changed in the plan

**Deliverable 2 of work package 1 carried a stale copy of a spec number.**
It asked for `beta` and `se` within 1e-4 absolute against plink2, while
"How it is verified" of "The logistic model" asks for 1e-4 times the `se` of
that variant. On the six literal variants `se` runs from 0.219 to 0.324, so
the spec's form is between three and five times tighter, and it is the form
that the spec's rule about tolerances gives: a tolerance is against the
scale of what is estimated. The deliverable now says what the spec says.

The number was copied into the plan before the spec's tolerance rules were
rewritten, and this is the second plan in which a copied number went stale
that way. The pattern is worth more than the instance: a number in a plan is
a copy with no way of knowing that its source moved.

**The plan said the spec had one open point that touched it, and it has
two.** Open 2 of the spec, a variant there is nothing left to test, became
one rule in four places on 24 September 2026, and two of those four are
tests of this plan rather than of `gwas-linear`. The plan now says so, with
what each task does under the meanwhile and what it would do under the other
answer.

## Work package 1: the logistic model

### Task 1.1, the null fit and the score test

`crates/popnei/src/gwas/logistic.rs` fits the null model by iteratively
reweighted least squares, which turns each step of a logistic fit into a
weighted linear one, and tests a whole block against the one Cholesky
factorization of the weighted design, with nothing factored per variant. The
core's `calc_gwas` no longer refuses a binomial trait with no kinship. The
core has 791 tests where it had 787, the same count on both linear algebra
backends.

**Against R.** Over the six literal variants the score statistic is within
6.024e-4 of R's where 1e-3 is allowed, and the p-value within 1.430e-4 in
`log10` where 1e-3 is allowed. The two backends give the same ten digits;
they part at 5.3e-14 of a statistic of 12.48. Those two tolerances are the
spec's and were not lowered, because what they measure is where R's own fit
stopped and not popnei's arithmetic: R's `glm` converges to 1e-8 in the
deviance.

**Against numpy, on a fixture of eight individuals.** `beta` and `se` are
each within 1.5e-15 times the `se` numpy gives, 0.8282, and the p-value
within 1e-14 in `log10`. Both were lowered until they failed on both
backends: the first breaks at 5e-16 of that `se`, where the effect is
6.70e-16 away on faer and 2.68e-16 on Accelerate and the standard error
2.68e-16 on either, and was set at three times the break; the second breaks
at 4e-15, where the p-value is 4.44e-15 away on Accelerate and 4.63e-15 on
faer, and was set at 2.5 times it. What the p-value costs above the other
two is the distribution: popnei reads the chi square off `erfc` and scipy
computes the same function another way.

**Open 2 has a reproduction at a score test, and it is an ordinary
mistake.** The covariate is a variant's own dosages in units a tenth of
theirs, which is what a user gets by putting a genotype in as a covariate.
The denominator of that variant's score test comes to 4.44e-16 on Accelerate
and to exactly 0 on faer, against a threshold of 4.19e-15. With no guard the
row reads `beta` 0, `se` 4.75e7 and `p_value` 1 on one backend and a NaN or
an infinity on the other. Two of the four places of that rule were reached
on fixtures built to break them, a collinear design and an identity kinship;
this one is reached by a mistake a user makes. That is a different claim
about the risk than "it is reachable", and it is the third case this week
where the two turned out not to be the same claim.

**Thirteen lines went into the spec before the code, as f7aa25c.** The
logistic fit can end before its 50 steps: once the chances it fits reach 0
and 1 the weights are 0, and the weighted design taken against itself is no
longer a matrix a Cholesky factorization accepts. The design's own columns
cannot cause it, since a design whose columns are not independent is refused
before any fit, so the weights are the only route. pyNei never meets it,
because it solves each step with an LU factorization, which answers a matrix
a Cholesky refuses. Measured on eight individuals whose covariate is 0 to 7
and whose four above 3 have the condition: popnei stops at step 45 on both
backends where numpy 2.5.3 runs all 50 and is still moving, at an intercept
of -299.6 and an effect of 84.0. It is a divergence from pyNei that a user
does not see, since both end in the same error, and `docs/objectives.md`
asks for those to be written down. The error is `GwasFitDidNotSettle`, which
names the model and the rounds it ran, and which a user meets as a
`ValueError`.

**For task 1.3.** The Python binding maps that new error through its
wildcard arm rather than by name, so the case has to be listed there and in
`tests/reference/gwas/refusals_of_both_layers.json`, which both test suites
walk.

### Task 1.2, the Wald fit of each variant

Every variant of a study of a binomial trait with no kinship is now fitted
with the variant in the design, started at the null fit's coefficients and
an effect of 0 for it. The standard error is read by solving the
factorization against the last column of the identity, which is that column
of the inverse at the cost of the coefficients squared rather than cubed.
The core has 793 tests where it had 791, three Wald tests in and one refusal
test out, the same count on both backends.

**Against plink2, on the six literals.** The effect is 2.69e-6 of that
variant's standard error away at the worst where 1e-5 is allowed, the
standard error 5.30e-5 of itself where 1e-4 is allowed, and the p-value
1.88e-4 of itself where 5e-3 is allowed. All three bounds are the spec's and
were kept rather than lowered: what they measure is that plink2 stops its
logistic fit earlier than popnei does. The two backends sit 4.4e-16 of a
standard error from each other, ten orders of magnitude below the bound.

**Against pyNei, which is deliverable 4 and was measurable without the
bindings.** `beta` within 3.06e-15 of the standard error, `se` 3.22e-15 of
itself, the p-value 6.0e-14 of itself, and the same one variant with NaN,
against the 1e-9 the spec asks for.

**Deliverable 3 holds.** The variants popnei gives NaN are `var0006` and no
other, which is exactly the row plink2 fell back to a Firth penalized
regression for, and the cargo test reads that set out of plink2's own file
rather than carrying a copy of it. Nothing anywhere records why a variant is
NaN, which is the meanwhile of Open 1.

**The whole-column bound on `se` of deliverable 2 cannot be met, and pyNei
is why.** Over the 1199 variants, popnei's worst `se` is 1.334e-4 of that
variant's `se` from plink2, at `var0179`, against the 1e-4 the spec and the
plan ask for. pyNei on the same panel gives 1.334e-4 at `var0179` too, and
popnei's other two worsts match pyNei's to four digits, 2.103e-5 for `beta`
and 1.924e-3 for the p-value. So the bound measures the distance between a
correct fit and plink2's, which stops earlier, and no implementation meets
it. pyNei's own test passes only because it holds `se` to 1e-4 absolute,
worst 6.66e-5, which is the form this spec argues against: it holds on a
panel whose values are small and breaks on data whose values are larger. The
number was in the spec as well as in the plan, and the session that owns the
spec set it at **5e-4 times that variant's `se`**, 3.7 times the worst
measured, with both numbers written into the spec so that the next person
does not tighten it back. Deliverable 2 of this plan now says the same. The
1e-4 on `beta` stands, worst 2.103e-5, and so do the six literals.

That is the third bound of this spec that a measurement has moved: one
absolute where the printing was relative, one tighter than the rounding of
the file it compared against, and this one that no fit meets. All three had
the same cause, a number written without running the panel it applies to,
and all three were found by someone building against the spec rather than by
re-reading it.

**The comparison with pyNei is now per model.** 1e-9 is its ceiling and not
its value, and each model's bound is set where it breaks. The logistic
model's three measurements above, 3.06e-15, 3.22e-15 and 6.0e-14, put it
near 1e-12, and deliverable 4 of this work package now asks for 1e-12. The
linear mixed model is the opposite case, where 1e-9 sits at the noise of its
own search, and is why the rule is per model at all.

**Two of the three marks of a runaway are unreachable under popnei's own
factorization.** The effect passing 30 catches `var0006` at round 29 and the
fixture's separating variant at round 29; the singular system catches a
variant that is its own covariate, at round 4 on Accelerate and round 1 on
faer. A step that is not finite, and a fit still moving after 50 rounds, are
in the code and in no test: with a Cholesky the effect passes 30 or the
factorization refuses the system first, and no fixture reached either. numpy
with an LU, which is how pyNei solves, reaches the third of them on the
collinear variant and gives the same three NaNs. So all three marks are
reproduced as the spec asks and two of them are dead code here. The spec now
carries that with the measurements, including that numpy with an LU reaches
the third and gives the same three NaNs, because dead code with no
explanation reads as a mistake: they are what pyNei marks, and an
implementation that factored another way would need them.

**A shared list was edited, one case of it.** The refusal of a binomial
trait with no kinship is gone, and both suites walk
`tests/reference/gwas/refusals_of_both_layers.json` and assert that its
cases match the calls they have, so leaving that case would have turned both
red. It now asks for the logistic mixed model, a binomial trait with a
kinship, which the same message refuses. Two lines of the list and one call
in each suite; nothing of the bindings or the packages moved.

### Task 1.3, the trait through the bindings and the packages

The core and both binding crates already carried the trait end to end, so
no model code was needed above the core. What changed:
`crates/popnei-python/src/errors.rs` now lists `GwasFitDidNotSettle` by name
instead of letting it reach the wildcard arm, where it arrived with the VCF
path glued to its message; the prose saying the logistic model was being
written is gone from both packages and from the `# Errors` comments of the
wasm binding; `tests/reference/gwas/refusals_of_both_layers.json` gains two
cases, written into both suites. The pytest suite goes from 499 to 506 and
the node one from 325 to 328.

### The four deliverables of work package 1

Each was run by the orchestrator, and the worst of each comparison was
measured again from a script of its own rather than read out of the task's
report. All of them on 24 September 2026, on Accelerate and on faer alike.

| what | command | worst | allowed |
|---|---|---|---|
| 1, the score test is R's | `uv run pytest tests/test_gwas.py -k the_score_test_of_every_variant` | statistic 1.589e-3, p 3.749e-4 in `log10` | 1e-2, 1e-3 |
| 2, the Wald test is plink2's | `uv run pytest tests/test_gwas.py -k the_wald_test_of_every_variant` | `beta` 2.103e-5 of the `se`, `se` 1.334e-4 of it, p 1.924e-3 relative | 1e-4, 5e-4, 5e-3 |
| 3, the runaway variant is plink2's | `uv run pytest tests/test_gwas.py -k the_variants_with_no_answer` | the variants with a NaN p-value are `var0006`, and plink2's `FIRTH?` `Y` rows are `var0006` | the two sets equal |
| 4, popnei and pyNei agree | `uv run pytest tests/test_gwas.py -k every_variant_of_a_logistic_panel` and `npm test` | 1.105e-14 over both panels and both tests; the six node literals at 2.688e-6, 5.301e-5 and 1.878e-4 | 1e-12; 1e-5, 1e-4, 5e-3 |

The node numbers are the same three the cargo test measures natively, so
none of what they carry is WebAssembly's rounding.

Deliverable 4 has more room than the spec expects. The spec puts this
model's bound near 1e-12 from the core's Wald test alone, which sits at
3.06e-15; through Python over both panels and both tests the worst is
1.105e-14, so 1e-12 is 90 times the worst and not 300 times. The bound is
right and the sentence that explains it is now out by that much.

## How the work went

This last section is not written for the owner, who can stop here. It is for
whoever next revises a skill or writes a plan, and it holds what would
change one of those.

**A tolerance that passes for the wrong reason is what the next test
copies.** Task 1.1 came back comparing with numpy at a share of each value,
where the spec says `beta` and `se` go against a share of `se` and a
p-value goes in `log10`. On that fixture the two forms are the same bound,
because `beta` is 0.855 and `se` is 0.828, and the doc comment said so and
argued the form was therefore fine. It was sent back, and the rewrite also
tightened the bound from 3e-14 to 1.5e-15, since the loose one had been
chosen to cover a p-value that now has a bound of its own and was 45 times
the worst measured. Both halves of the spec's rule were being broken by one
constant, the form and the two-or-three-times-the-break, and only the form
was visible.

**What the tasks cost.** The split of the module cost its subagent 174000
tokens and 91 tool calls. Task 1.1 cost 292000 tokens and 114 tool calls,
and the one correction 11000 more over 11 calls. Task 1.2 cost 266000
tokens and 116 tool calls and needed nothing sent back. The orchestrator ran
the nine checks itself after each, which is about four minutes of wall clock
a time and no context to speak of, since it reads the last line of each.

**A task that reaches into a shared file is cheaper than a task that
leaves it red.** Task 1.2 was told to leave the packages alone and edited
one case of the refusals list and one call in each suite, because turning
the Wald test on falsifies a refusal both suites assert. It said what it had
done and why and offered to undo it. That is the right shape for a
boundary that cannot be held: the alternative was to hand back a branch
whose node suite was knowingly red until task 1.3.
