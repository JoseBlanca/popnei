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
| `cargo test -p popnei-linalg --no-default-features` | 149 |
| `cargo wasm-check` | clean on both wasm targets |
| `uv run ruff format --check && uv run ruff check` | clean |
| `uv run maturin develop && uv run pytest` | 499 passed in 21.03 s |
| `npm run build && npm test` in `js/popnei` | 325 passed, 0 failed |

The faer run is the one that matters most to this plan, because it is the
linear algebra a browser runs and because every tolerance below has to hold
on it as well as on Accelerate.

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

## How the work went

This last section is not written for the owner, who can stop here. It is for
whoever next revises a skill or writes a plan, and it holds what would
change one of those.

Nothing to record yet.
