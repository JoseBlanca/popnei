# Report: the association study of a binomial trait

24 September 2026. It records how `docs/plans/gwas-logistic.md` was carried
out, on the branch `plan/gwas-logistic`.

**The plan is done.** Every task is ticked, every deliverable of the three
work packages was checked by the orchestrator running its command, all three
work packages were reviewed, and the plan's own final check passes.

## What exists now that did not

`calc_gwas` in Python and `calcGwas` in TypeScript test every variant of a
dataset against a **binomial** trait, which is the half of
`docs/specs/gwas.md` that was left. With covariates and no kinship that is
the logistic model, whose effect is a log odds ratio. It takes either of two
tests: the Wald test, which is what it does when the user asks for none, and
which is checked against plink2's `--glm`; or the score test on request,
checked against R's own score test of the same model, which R calls Rao's. With a kinship it is the logistic
mixed model, fitted by penalized quasi-likelihood with its covariance
factored rather than inverted, whose score test is checked against GMMAT.
And both mixed models can now take the GRAMMAR-Gamma approximation, which
makes the work per variant linear in the individuals instead of quadratic.

All four models of the spec are built, in every layer, and they run in a
browser.

At the commit this branch starts from: 787 tests in the core crate, 499
pytest, 325 node. Now: **842 in the core crate** with 2 ignored, the same
842 on the faer backend, 150 and 136 in the linear algebra crate, **514
pytest** and **332 node**. `cargo fmt`, `cargo clippy` with every target and
warnings denied, `cargo wasm-check` for both wasm targets and ruff are
clean.

## What is asked of the owner

**The merge**, which they have ordered, and four things that are not this
plan's to settle.

1. **A wrong value users can reach, fixed here by the owner's order, in a
   crate this plan does not own.** Accelerate's `dpotri` returns a wrong
   inverse. Measured on this machine: about 1 inversion in 1500 is wrong
   when nothing else of popnei is running, because the routine starts
   threads of its own, and about 1 in 7 when another thread of the same
   process is doing any other work of that library. It reached both mixed
   models, which are the two that form an inverse. The inverse
   is now built from a triangular solve and a product instead. The crate's
   other LAPACK calls were tested and are clean. What is not known is
   whether faer's own inverse has the same fault, and which macOS versions
   besides this one do.
2. **A skill was changed and it is in force**,
   `.claude/skills/following-plans/SKILL.md`, which now tells a task prompt
   to carry today's date. The evidence: on 24 September 2026 two sessions
   carrying out two plans wrote 28 dates that were a day in the future,
   every one ahead and none behind, and 14 of those came from a single
   subagent working through one list of fixes, which worked the date out
   once and was wrong 14 times. Nine of the 28 dated a measurement, where it
   costs most, because a reader checking whether a bound still holds meets a
   date that has not happened and doubts the number rather than the date. A
   subagent has no clock and infers the date from the files it reads, so one
   file a day ahead makes the next writer a day ahead. The writing skill was
   given the other half of the fix, where a date comes from, by the session
   that owned it. Reverting this is one commit and nothing depends on it.
3. **Five things `docs/specs/gwas.md` needs**, listed under "What the spec
   needs and no session owns", including that its stated agreement with
   pyNei is a numpy prototype's figure and not this code's.
4. **Five defects found outside this plan's scope**, listed under "What was
   found and not fixed here", of which the largest is that three
   calculations never raise a `KeyboardInterrupt` at all.


## The names this report uses

An association study tests each variant on its own and reports the effect of
one more copy of an allele, the uncertainty of that effect, and the chance
of seeing an effect at least that large if the variant had none.

**The two tests.** A *Wald test* fits the variant, divides its effect by its
own standard error and asks how extreme that ratio is, so it needs a fit for
every variant. A *score test* asks instead how steeply the likelihood would
rise if the effect moved away from zero, measured at the fit with no variant
in it, so it needs no fit per variant and is the cheaper of the two. They
agree closely where the effect is small.

**A binomial trait.** The trait is 0 or 1, so the effect is a *log odds
ratio*: the change in the log odds of being a 1 for one more copy. The model
is fitted by *iteratively reweighted least squares*, which turns each step
into a weighted linear fit. With a kinship it is fitted by *penalized
quasi-likelihood*, which does the same but with the relatedness carried as a
second source of variance, and the search for how large that variance is
uses two things this report names: a *trace taken from an identity*, which
gets a quantity that would otherwise need every entry of an inverse, and the
*average information*, which is the step rule of that search.

**The programs it is checked against.** *plink2* and *GMMAT* are the two
programs that compute these models outside popnei: plink2 for the two with
no kinship, GMMAT for the two with one. *rrBLUP* is a third, used by the
previous plan for the continuous trait, and it appears here only because the
script that makes the reference data runs all of them at once. *pyNei* is
the Python library popnei reimplements, and is the oracle for everything the
other three do not compute.

**The two linear algebra backends.** popnei's arithmetic runs either on the
system's own library — on this machine Apple's *Accelerate* — or on *faer*, a
library written in Rust, which is what a browser runs. Every tolerance in
this report was chosen against both.

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
`heritability` set to `None`. Measured on 24 September 2026, it does not
reach the logistic mixed model at all: the other half of that model's
covariance is the reciprocals of the weights, which differ from individual
to individual, so an identity kinship still tells the two variances apart
and the fit lands at a variance of the kinship effect of 0.1274 in 6
steps. Whatever the owner answers, no test
of this plan moves.

**Open 4, how negative an eigenvalue is still rounding.** Meanwhile: refuse
when the smallest eigenvalue of the kinship is below a tenth of the largest.
**This report said that was already built and it is not.** The review of
work package 2 looked for it: `crates/popnei/src/gwas/linear_mixed.rs`
clamps with `eigenvalue.max(0.0)` and nothing anywhere compares the smallest
against a fraction of the largest. Forcing `panel_called`'s kinship to a
smallest eigenvalue of -29 per cent of its largest, the linear mixed model
runs and gives `genetic_variance` 1.221615 and `heritability` 0.781096,
which are the numbers Open 4 itself quotes as what is "clamped in silence",
while the logistic mixed model refuses the same input because its Cholesky
will not factor it. So the two mixed models disagree on one input, and a
user can read a heritability off a matrix that is not a kinship.

The claim came into this report from the plan's own reading of what
`gwas-linear` had built, and was never run. It is the same failure this
review has found in other people's numbers, here in my own: a claim
repeated from the document that made it, without running it.

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
about the risk than "it is reachable": the first says a user will meet it,
the second only that a fixture can be built for it.

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

### What the review of work package 1 found

Seven reviewers read the work package, one per category, each with a fresh
context and none of them the writer. Between them they raised twenty-two
things. None was a wrong number that either reference panel exercises, and
one was a wrong number a user can reach.

**The wrong number, which is now Open 5 of the spec.** popnei's default
build answered `p_value` 0.9999996244683889, with `beta` -18.24 and `se`
3.9e7, for a variant whose effect has no finite value, where pyNei gives
NaN and popnei's own faer build gives NaN. Eight individuals, trait
`0 0 0 1 0 1 1 1`, a covariate of 0 to 7 and dosages `2 1 0 1 2 0 1 0`. None
of the marks fires: pyNei catches it by running out of rounds, and popnei
does not reach that because its system is nearly singular, so the steps
shrink while the coefficients are still walking and the fit declares itself
settled. The 30 is not it either — the variant's own effect is -18.97 and it
is the intercept that passes 30, which neither library reads.

The meanwhile, built and measured: a Cholesky pivot that has fallen to the
tested individuals times 2.2e-16 of the largest pivot marks the fit a
runaway, which is the rule this module already applies in the four places of
Open 2. On both panels and both backends no variant that was answered lost
its answer and none gained one; the smallest pivot of an answered fit is
2.600e-2 of the largest on `panel_called` and 4.730e-5 on the panel with 3
genotypes missing in 100, against a threshold of 4.44e-14, which is nine
orders of headroom. The case above now gives three NaNs, which I ran myself.

**The meanwhile narrows Open 5 and does not close it.** One of the new tests
has a fixture that settles at an effect of 36.45 with a standard error of
2.0e7, and the pivot does not catch it; the 30 does. A fit can still stop
with a collapsed system the pivot does not see, and an `se` of 2.0e7 beside
an effect of 36 is the same signature as the case above.

**The score test's denominator is formed and no longer subtracted.** It was
`x' w x` minus what the covariates explain, two nearly equal numbers
cancelling exactly where Open 2's threshold acts. Measured over six decades,
the subtracted form is out by 4.6e-3 of itself at 1.3 times the threshold
and by 29 per cent at the threshold, while the formed one falls as the
square of the collinearity throughout, which is what says which of the two
is the accurate one. So the guard had been reading a quantity whose error
was larger than the thing it was testing. Forming it costs one more product
per block, into a buffer `linear.rs` already keeps, and it is the second
place popnei departs from pyNei's formula rather than the first.

**That change is invisible to every check in the plan.** The smallest
denominator `panel_called` reaches is 0.297 of its scale, so no literal
moves: no panel value changed by more than 3.5e-14 of itself, and all four
deliverables give the same numbers to four digits. Its evidence is the
measurement and not the suite.

**Four tests could not fail, and each is now pinned by the mutation that
found it.** The runaway threshold of 30 could be replaced by infinity; the
guard on the variance could be made to accept a NaN; the null fit's fifty
round branch could be switched off; and the Open 2 threshold could be
replaced by a comparison against 0 and still pass on faer, where the
subtracted denominator was exactly 0. That last one closed itself: a sum of
terms that are not negative cannot fall below 0, so with the denominator
formed the mutation now fails on both backends.

**A refusal named the wrong cause.** A design whose covariates are nearly
the same was refused with a message about a covariate that separates the
cases from the controls, telling the user to take it out. popnei's rank
check uses numpy's tolerance and the Cholesky's is tighter, so between them
there is a band: at a correlation of 1 − 5e-13 the fit runs three rounds and
the factorization refuses, and only once the covariates agree to within
about 1e-14 does the rank check catch them. The message now names both
causes and their different remedies.

**The null fit had no guard for a step that is not finite** where the per
variant fit has one, and its comment said it had one. A user would have met
a `RuntimeError`, which this project reserves for its own defects, for what
is their data. Neither of the two reviewers that found it could build the
input that reaches it, so the guard is insurance and the comment says so.

Three findings were reached independently by two reviewers each: that
missing guard, the sentence in both packages about which rows are NaN, and a
doc comment claiming 3.2e-15 where the quantity measures 3.368e-15.

**The bound against pyNei was 90 times the noise and is now 2.7 times it.**
`OF_PYNEI_LOGISTIC` was 1e-12 against a worst of 1.105e-14. Two defects a
reviewer planted, dropping either of the fit's final reweightings, move the
columns by 1.1e-9 and 3.0e-9, so both are caught at either bound; what 1e-12
cost is the 1e-13 class of error, which is the only class this bound can
catch at all, the larger ones being caught by plink2 and R. It is 3e-14 in
the spec, in deliverable 4 and in the test.

### What the review changed above the core

The node suite had no reference number for the logistic score test at all,
which is the worst instance of a test that cannot fail that this plan has
met: a reviewer multiplied that test's effect by 1.5 in the core, rebuilt
the WebAssembly and watched all 328 node tests pass. It now asserts R's six
score statistics and p-values, and the subagent that added them proved the
gap the same way, making the mutation itself and checking its new test was
the single failure of 329 before reverting it. The suites are 507 pytest and
329 node, where they were 506 and 328.

Three other things a user reads were wrong. Both packages said a row of
three NaNs is a variant that separates the cases from the controls; it is
that only under the Wald test, the score test gives a number there, and a
variant that repeats a covariate and separates nobody gets the same three
NaNs. `covariate_effects` are log odds ratios for this model and trait units
for the other two, which no layer said, though the sibling field `beta`
carries exactly that sentence. And the Ctrl-C test did not run `calc_gwas`,
so nothing would have failed if the binding had stopped raising it.

Nine dates in the file were in the future, 25 or 26 September 2026. Seven
are measurements of this module and were remeasured and hold to the digits
given; three are statements about earlier edits of this branch and had the
date alone corrected.

### The four deliverables, after the fixes

Rerun by the orchestrator on the last commit of the work package, all nine
checks green at 799 core tests with 2 ignored on both backends, 149 and 136
in the linear algebra crate, 507 pytest and 329 node:

| what | worst | allowed |
|---|---|---|
| 1, the score test is R's | statistic 1.589e-3, p 3.749e-4 in `log10` | 1e-2, 1e-3 |
| 2, the Wald test is plink2's | `beta` 2.103e-5 of the `se`, `se` 1.334e-4 of it, p 1.924e-3 relative | 1e-4, 5e-4, 5e-3 |
| 3, the runaway variant is plink2's | both sets are `var0006` | equal |
| 4, popnei and pyNei agree | p-value 4.524e-14, `beta` 1.151e-14 of the `se`, `se` 3.368e-15 | 1e-13 |

Not one of those numbers moved when the denominator was formed and the
pivot rule went in, which is what says those two changes are invisible to
every check this plan has.

### What was found and not fixed here

Each of these is real, none belongs to this work package, and all of them
are the owner's to place.

- **The faer backend allocates about 39 times per variant in the Wald test**,
  where the BLAS one allocates nothing: 87463 allocations against 9447 over
  one pass of 2000 variants. The cause is the scratch buffer that
  `cholesky_lower` and `solve_with_cholesky` take on every call in
  `crates/popnei-linalg/src/faer.rs`, so the fix is in the linear algebra
  crate. It is about 39 million calls to the allocator for the million
  variants of `docs/objectives.md`, on the target with the least to spare.
  The claim in the core that said a pass allocates nothing per variant has
  been corrected rather than left standing.
- **Four calculations are missing from the Ctrl-C test** beside `calc_gwas`,
  which was added: `pca_of_variants`, `calc_rogers_huff_r2_matrix`,
  `calc_per_var_distribs` and `calc_per_individual_stats`. The last three
  have no call to the helper in their binding at all, so a Ctrl-C during
  them raises no `KeyboardInterrupt`.
- **`tests/test_io_vars.py` fails under a release build.** Its Ctrl-C test
  depends on a write taking 0.313 s in the debug build, which its own
  comment states.
- **`solve_with_cholesky` cannot check that the factor it is given factors
  the matrix the caller means**, by its own doc. The logistic module is safe
  because both come from the same call, but a later edit that reused a
  buffer would get a wrong answer with no error.
- **The TypeScript `stats` hands out `chrom` and `id` unfrozen**, where the
  same kind of column is frozen elsewhere in that package.

One finding was evaluated and left as it is. One case of
`tests/reference/gwas/refusals_of_both_layers.json`, the study with fewer
individuals than the design has columns plus two, carries the VCF path in
its message although it is refused before any variant is read. The count it
reports is of the dataset, so the file is part of what the message is about,
which is the same reason the refusal beside it deliberately carries a path.
The new assertion covers the other 24 cases and skips that one by name.

## Work package 2: the logistic mixed model

### Tasks 2.1 and 2.2, the linearization and the search over the variance

`crates/popnei/src/gwas/logistic_mixed.rs` fits the null model of a binomial
trait with a kinship: one pass of the penalized quasi-likelihood at a given
variance of the kinship effect, with the covariance factored by a Cholesky
and applied by solving and never inverted, and a search over that variance
using the trace from the identity, the average information and a bracket.
The core has 812 tests where work package 1 left 799.

**Deliverable 1, the null model against GMMAT, holds.** The variance of the
kinship effect is 1.5080506977719876, which is 6.302e-6 from GMMAT's
1.508057 and 63 per cent of the 1e-5 allowed. It was 1.50805069777198586
until the review: writing the trace's right hand side as the square root of
the reciprocal weight, rather than the reciprocal of the square root, moved
popnei's own digits by 8.8e-16 of themselves. That is a rounding and not a
correction of the value; what it corrects is that the score at a variance of
0 is now exactly 0 by construction rather than by two hundred quotients
happening to round to 1. No reference literal moved, and no suite asserts
popnei's own digits, which is why nothing went red.
The three covariate effects are within 1.06e-6 of GMMAT's, 11 per cent of
their bound, the same on both backends. `residual_variance` and
`heritability` are `None`, since a logistic model has no free residual
variance.

**Deliverable 3, the fit forms no inverse while it iterates, holds.** One
inverse, formed at the end because the score test wants the projection
matrix as a matrix, and 22 factorizations of the covariance over 8 steps on
the variance, on both backends and both panels. The 22 is the number
`docs/reports/glmm-method/README.md` measured for pyNei's fit, which
inverts once per linearization where this one factors.

**A guard was written, measured and then removed, which is the right way
round.** Task 2.2 first refused a step on the variance that is not finite.
Measuring showed that takes away a right answer: a kinship of all zeros
walks to a variance of exactly 0 in 12 steps, which is the boundary the spec
describes, because the bracket sends every step not above 0 to a quarter of
the variance. The guard would have refused a study that has an answer. It is
gone, and that case is now the one fixture that runs the whole search
through the division by the variance — the division that is 0 exactly rather
than nearly, which is the trap the owner named for this work package.

**The stopping rule's denominator is pinned by nothing, and the report said
otherwise for an hour.** Task 2.1 left a "plus 1" in `change / (predictor +
1)` that none of its tests could tell from a 1. Task 2.2 reported that the
whole fit could: that it left the variance identical but moved the covariate
effects by 3.1e-13 and the rounds from 22 to 23. The review ran the mutation
and neither is true. With the `+ 1` dropped the fit is **bit-identical** on
both backends — the same 22 rounds, the same 8 steps, the same variance to
all 18 digits and the same three effects — and all 15 tests of the module
stay green. So nothing in this plan pins it, and the 3.1e-13 was a number
that no run produced.

**Open 3 does not reach this model.** The meanwhile for a kinship that does
not identify the two variances was written for the linear mixed model, where
the criterion goes flat, and whether the logistic one had the same
degeneracy was not known. It does not: the other half of the covariance is
the reciprocals of the weights, which differ from individual to individual,
so an identity kinship leaves the two apart and the fit lands at a variance
of the kinship effect of 0.1274 in 6 steps. A test says so, and no meanwhile was needed here.

### Task 2.3, the score test and the model through every layer

`calc_gwas(..., trait="binomial", kinship=k)` is the fourth model and the
last one the spec describes. The score test reads the projection matrix the
null fit formed, and the Wald test of this model stays refused, as the spec
asks. The core has 813 tests, the pytest suite 512 and the node one 330.

**Deliverable 2, the score test against GMMAT on both panels, holds.** The
worst over 1200 variants is `1 / se²` at 5.342e-6 of GMMAT's `VAR`, at
`var1032` of `panel_called`, where 1e-5 is allowed and where GMMAT's own six
printed digits account for 4.92e-6 of it; and the p-value at 8.497e-6 in
`log10`, at `var0520` of `panel_called`, where 1e-4 is allowed. That second
figure is the one `docs/reports/glmm-method/README.md` printed for the same
comparison, to every digit. The six literals of both panels are identical on
Accelerate, on faer and under node.

**Deliverable 4 holds and its bound is now measured rather than inherited.**
The plan said 1e-9, which is the spec's ceiling and not a value. Measured
over both panels and both tests, the worst is the p-value at 6.794e-14 in
`log10`, at `var0115` of the panel with genotypes missing on faer, with
Accelerate at 6.384e-14 on `var0892`; then `beta` at 3.916e-14, the genetic
variance at 5.801e-14, `se` at 1.327e-14 and the covariate effects at
8.976e-15. The bound is 1.5e-13, 2.2 times the worst, which is the ratio the
plain logistic model's bound takes.

**Open 2's fourth place is reached on data at last.** A fixture of eight
individuals in two families, with a covariate that is a tenth of the first
variant's dosages, leaves `x' p x` at 2.776e-17 on Accelerate and 1.284e-16
on faer, against a threshold of 2.780e-15. Unguarded, that row answers
`beta` -7136 with an `se` of 1.898e8 on one backend and -1543 with 8.826e7
on the other. Three of the rule's four places now have a case behind them,
and the fourth, the linear mixed model's score test, is not this plan's.

**A kinship of 1e-162 times the identity answers a variance of 0, and the
owner decided to leave it.** The true optimum there is 1.3e159 and the
average information underflows, so popnei says the kinship explains nothing
for a matrix that explains as much as any kinship does; what fails is that
the parameter cannot represent the answer at that scale. The case for
refusing was that a silent "nothing here" is a wrong answer a user acts on.
The case for leaving it, which the owner took on 24 September 2026, is that
no kinship popnei computes comes near it — `calc_kinship` gives entries of
order one — so it is a contrived input, and a refusal would cost a
comparison on every fit and a case in the spec for something nobody will
meet.

### A wrong value that is nobody's plan, and is the owner's to place

**On the default build, the same inversion gives different answers depending
on how many threads are running.** It reaches the two mixed models, which
are the only models that form an inverse.

What I ran myself, at b47f52c on this machine:

| run | failures |
|---|---|
| `cargo test -p popnei --lib gwas::logistic_mixed` | 3 of 20 |
| the same with `--no-default-features`, which is faer | 0 of 20 |
| the same with `-- --test-threads=1` | 0 of 20 |

Two tests fail, both reading the projection matrix, which is what
`invert_with_cholesky` forms. The subagent that found it measured the
operation directly: 8 threads each inverting one Cholesky factor 200 times
give answers up to 1.08e-4 apart, where faer gives exactly 0, and where
`cholesky_lower` and `solve_with_cholesky` are stable on both backends. That
the faer build and the single-threaded run are clean, with the same test
code, is what says it is the backend and not popnei's fixtures sharing
state.

It is not new: at 24dc716, which is `main` with none of this work package's
score test on it, the same filter failed 2 of 10.

**It is settled, and it is not conditioning.** The covariance of the working
trait has a condition number of 4.160 at a variance of 0 and 15.38 at
GMMAT's, so a different summation order would move the inverse by about
15 times 2.2e-16, which is 3e-15 relative and eleven orders below the
1.08e-4 observed. The matrix the wrong answer comes from is well conditioned
by every measure: on one run of 200 individuals its eigenvalues span 4.368
to 11.04, a condition number of 2.53, and the wrong inverse fails
`sigma · inv = I` at 4.985e-4 where the right one fails at 3.109e-15. So it
is `dpotri` itself under concurrency and nothing about this problem.

It is wider than two studies. One thread inverting beside seven threads
doing nothing but matrix products gives 288 wrong inversions out of 2000.
The triangular solve and the products are clean under the same test, so
`dpotri` alone is the unsafe call. Under `VECLIB_MAXIMUM_THREADS=1`, or with
one caller thread, the difference is exactly 0; with two caller threads it
is already 4.736e-5.

**The owner ordered it fixed on this branch, and it is.** The inverse is now
built from the two stable operations the crate already had, solving the
triangular factor against the identity and multiplying the result by its own
transpose. `cargo test -p popnei --lib gwas::logistic_mixed`, which failed 3
or 4 runs in 20 before, fails 0 in 20 after, measured by the orchestrator;
the guard test written with the fix fails 20 of 20 against the old route.
No number of either mixed model moved.

**The defect is not popnei's, and it is worse than the review found.**
`accelerate-src` emits nothing but `-framework Accelerate`, `lapack-sys`
declares `dpotri_` over 32-bit integers, and the test binary has no
`$NEWLAPACK` symbol, so the crate's declaration and Accelerate's legacy
interface match; and a mismatch would be wrong on one thread, where this is
right. But `dpotri` is wrong **with no other popnei thread at all**: 2000
inversions on one thread gave 2, 2, 1 and 0 different answers over four
runs, because Accelerate's own threads inside the routine suffice. Other
work in the process only takes the rate from about 1 in 1500 to about 1 in
7. So a user running one study on one thread could already get a wrong
projection matrix; concurrency made it frequent rather than possible.

**What the fix costs.** At the 200 individuals popnei inverts at, the new
route is 0.38 of `dpotri`'s time, 0.0000786 s against 0.000209 s. Above
that it is dearer: 1.8 times at 1000, 1.2 at 2000, 2.0 at 5000 and 1.7 at
10000, where it is 3.540 s against 2.108 s. It also needs two matrices of
the individuals squared where `dpotri` needed none.

**The crate's other LAPACK calls are clean.** One thread beside seven doing
products gives no difference at all for the eigendecomposition, the singular
values and the QR, at three shapes each, nor for either triangular solve.
Two things are not settled: faer's own inverse was not put under the test,
so this is about Accelerate alone, and which macOS versions besides this one
have it is unknown.


Who it reaches: popnei forms this inverse once per study, so one study in
one process is safe. A user running two studies in threads in one process
would see numbers move in the fourth digit.

### What the review of work package 2 changed above the core

**A user was told the model they wanted was not built.** Both packages
opened with "three of the four models" and described the fourth as built
twenty-six lines later, so `help(popnei.gwas)` and the hover over
`calcGwas` said the opposite of the code. Two terms that did work without
being explained, the fitting method and the working trait, are in
`docs/glossary.md` now.

**Which errors name the file is decided in one exhaustive match.** The
Python binding classified each error by hand and let anything unlisted fall
through a wildcard that glued the VCF path onto a message about the user's
own arguments; the comment above that list recorded it happening twice, a
day lost each time, and the only check needed somebody to remember to add a
case first. The decision now lives in the core as
`popnei::Error::names_the_file`, which the compiler will not let a new
variant skip. Every message a user sees is unchanged.

**An error told a user to wait for something that had shipped.** With all
four models written, the only arm still reaching `GwasModelNotBuilt` is a
mixed model whose kinship has vanished, which is popnei's defect and not a
study anyone asked for. It says so, and is a `RuntimeError`.

**A refusal named a model the user had not asked for.** The logistic mixed
model fits the plain logistic null as its starting point, and when that
inner fit ran away the user was told "a binomial trait with no kinship is a
logistic regression" — for a call that brought a kinship. Six individuals
whose covariate is the trait itself reach it. That one refusal is remapped;
every other error of the same inner fit passes through unchanged, which the
plain model's own tests and the shared refusals list assert.

**What is left imprecise, and I am leaving it.** The remapped message offers
three things to look at, and the third, a kinship asking for a random effect
the trait cannot fit, cannot be the cause of a runaway in the starting fit,
which has no kinship in it. Separating them needs an error case of its own
and a line in the spec. The spec asks only that the refusal name the model
and the round, which it does.

**Three numbers and a date.** The node suite's shared tolerance said its
worst user was the linear mixed model's genetic variance at 1.273e-6, 13
per cent of what is allowed; measured under node, the worst user is now the
logistic mixed model's variance of the kinship effect at 6.033e-6, 60 per
cent, the same figure pytest measures natively. The round at which a
collapsed weight is refused is 26 on both backends, as the spec says, and
the test asserts the number rather than that it is above 0. And one date
two days ahead was not in the file the review named: it had moved into the
core with the exhaustive match.

### What the spec needs and no session owns

`docs/specs/gwas.md` is on `main` and the session that owned it has ended,
so these are the owner's:

- **Its agreement with pyNei is the prototype's, not this code's.** The spec
  says the variance agrees to 2.9e-15 and the covariate effects to 3.3e-15;
  measured through both packages, 5.404e-14 and 8.415e-15. Those two figures
  are `docs/reports/glmm-method/README.md`'s, which measured a numpy
  prototype. The third number of the same sentence, 8.497e-6 against GMMAT,
  is reproduced to every digit. The pytest suite already carries the true
  values, so the spec and the suite disagree.
- **Open 4's meanwhile is unbuilt**, as recorded above.
- **A divergence from pyNei is not written down.** With a kinship of all
  zeros popnei answers a variance of 0 and pyNei raises, because pyNei's
  test of a step at 0 is false for a value that is not a number.
  `docs/objectives.md` asks for those to be recorded.
- **Open 3 does not say which model its meanwhile binds.** It was measured
  not to reach the logistic mixed model, and that measurement lives in a doc
  comment and in this report, not in the open point.
- **Open 2's fourth place now has a fixture**, whose numbers are in this
  report and in the code but not in the item.

## Work package 3: the GRAMMAR-Gamma approximation

Both mixed models can now replace the per variant product with the
projection matrix, which costs the square of the individuals, by one factor
times the squared length of the variant's centered dosages, which is linear
in them. The factor is the mean, over the first 100 variants that vary, of
the exact denominator divided by the approximate one.

**The factor is pyNei's.** 0.51730062 for the linear mixed model on the
panel with every genotype called, 3.4e-9 from pyNei's own, and 0.10532073
for the logistic mixed model, 1.8e-14 from it; on the panel with genotypes
missing, 0.52850585 and 0.10605886. The same on both backends.

**Deliverable 2 holds and spends 98 per cent of its bound.** Measured by the
orchestrator on the called panel with the linear mixed model, all 1200
variants answered by both modes: the median of `log10(p_approx / p_exact)`
is -5.1885e-4 where 0.1 is allowed, the largest is 0.51185 where 1.5 is
allowed, and the worst `beta` moves 0.48966 of itself where 0.5 is allowed.
Under the score test the largest is 0.4887 and the worst `beta` is the same.
Accelerate, faer and WebAssembly agree to 1e-8.

**Why that bound is nearly spent, which the plan did not expect.** The plan
says the factor is a mean of ratios that lie near each other, so that a
factor far from them is the sign of an error. They do not lie near each
other: over the 100 variants the factor is taken from, the ratios run from
0.3124 to 0.6704, a spread of 2.15 times, with a standard deviation 12.9 per
cent of the mean. The logistic mixed model's are tighter, 9.2 per cent. A
variant's effect moves by the share its own ratio sits from the factor, so
the worst variant moves nearly the whole of what pyNei's bound allows — and
pyNei's bound was set where pyNei's own worst case landed. The check passes
because the method is pyNei's, not because it has room.

Over all 1200 variants rather than the first 100, the linear mixed model's
factor is 0.518269 and the ratios run 0.307 to 0.771. So the 100 that pyNei
fixed under-sample the spread and leave the factor 0.19 per cent from the
whole-panel mean. `NUM_VARS_FOR_GAMMA` is unchanged, as the spec says it is
inherited and unmeasured.

**Which variants have no answer now depends on the argument.** Open 2's
threshold stops firing under the approximation: a factor above 0 times a sum
of squares holds no cancellation, so the variant that the exact test refuses
with three NaNs is answered with an effect of -5.7e-16 and a p-value of 1.
The rule itself is unchanged and is still compared against whichever
denominator the study formed. It is in the spec.

### What the review of work package 3 found

Four reviewers read it. Two findings are what the whole review process is
for, and neither could have been found by reading.

**A study whose every variant is NaN, in silence.** The check on the factor
was one-sided, on a quantity the code's own error text calls the rounding of
a cancellation, which falls on either side of 0. When it fell negative the
study was refused, as intended; when it fell positive the study was accepted
with a factor of order 1e-16, every denominator then sat under Open 2's
threshold, and the user got a table of NaN with nothing said. Measured by
the orchestrator over 193 designs of six individuals whose covariate is a
variant's own dosages, with the first block holding only that variant: 112
refused and **81 accepted with every variant NaN**. The factor is now taken
only over variants the projection leaves something of, and the same 193
designs are all refused.

**Every end-to-end test of the approximation passed best when the
approximation did nothing.** They bounded only how far the approximated
answer sits from the exact one, and zero distance is the best possible
score. A reviewer replaced the approximation with nothing in the logistic
mixed model and all 836 core tests stayed green, as did 513 pytest and 331
node; the result still reported that the approximation had been used,
because that flag is set where the argument is read and not where the
denominator is formed. A regression that quietly lost the approximation —
leaving a user paying the quadratic cost they asked to avoid — would have
shipped green. Every ceiling now has a floor beside it: the worst effect
moves 0.489655 for the linear mixed model and 0.393349 for the logistic one,
so a floor of 0.05 sits a factor of eight clear.

That reviewer also answered the question this report asked when work package
3 began. A constant error that does not grow with the structure of the panel
**is** caught: scaling the denominator by 1.02 turns the unit test and the
panel test red. The corridor the loose bound leaves is narrower than it
looked; what it did not catch was the approximation being absent altogether.

**Three guards no test executed**, each deleted in turn with the suite
staying green: the check that a block's rows and values agree, the half of
the factor check that catches a value which is not a number, and the early
return of the centered length on an empty slice. The second matters most,
since a factor that is not a number makes every p-value of the study one
too. Each has a test now.

**Four things that were not true** in doc comments: a count of three where
the test asserts four and averages four ratios, a fixture described as the
variant no individual has a genotype of where the code uses one every
individual is heterozygous for, a statistic called the spec's median where
the core takes the value at index 600 of the sorted absolute log ratios, and
a factor printed as -0.0000000000000003608224830031759.

**A defect reported as a wrong argument.** The branch that approximates a
denominator told a user with no kinship to supply one, on a path reachable
only when a kinship *was* supplied and the model came out non-mixed, which
is popnei's own fault.

## The plan's own final check

Three things beyond the sum of the work packages.

**The four reference programs still give what the literals say.**
`tests/reference/gwas/make_reference.py` was rerun on this machine on 24
September 2026, with plink2 v2.0.0-a.7.7 and R 4.6.1 with GMMAT 1.5.0 and
rrBLUP 4.6.3. Every one of the nine files came back **byte-identical** to
the committed one, and the script's own comparison with what pyNei stored
for the same genotypes gives a largest difference of **0.000e+00** over all
nine files and their 58 numeric columns. So the numbers this plan was
checked against do not depend on which library read the genotypes, and they
have not moved since the reference was made.

**popnei agrees with pyNei for all four models, on both panels.**
`uv run pytest tests/test_gwas.py -k "pyneis"` gives 13 passed, which is the
four models over the panels and the tests each takes, and is the comparison
`docs/objectives.md` asks for.

**The whole study reaches a browser.** `cargo wasm-check` is clean on both
wasm targets and `npm run build && npm test` passes in `js/popnei`.

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

**What a review of this size costs.** The seven reviewers cost 1.16 million
tokens between them, from 156000 to 183000 each, and 648 tool calls. The two
subagents that fixed what they found cost 294000 and 218000 over 340 calls.
So the review and its fixes cost about 1.7 million tokens against the
830000 the three tasks of the work package cost to write, which is twice the
writing. What it bought: one wrong value a user can reach, four tests that
could not fail, a refusal naming the wrong cause, and an arithmetic form
whose error at the threshold was larger than the quantity being tested.
Three of those were found by mutation rather than by reading, which is the
part of the `code-review` skill that earned its cost here.

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
