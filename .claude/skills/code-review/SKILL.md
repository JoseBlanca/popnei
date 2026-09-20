---
name: code-review
description: How code is reviewed in popnei and what is done with the findings. Use it after the tasks of a work package of an implementation plan are committed and before the work package is reported as done, or when the owner asks for a review of a commit, a branch or a module. The session that asks for the review sends one reviewer subagent per category, each with a fresh context, then evaluates every finding and fixes the ones that hold.
---

# Code review

The writer of a piece of code cannot review it, for the same reason the
writer of a text cannot: they know what it was meant to do and read that
into it. So the review is done by subagents that start with nothing but
the code, the spec and the `coding` skill, one for each category below,
and the session that wrote the code, or that the owner asked, is the
orchestrator: it sends them out, and then it acts on what they find.

A reviewer checks what the writer was told. The categories are the
sections of the `coding` skill plus two that no writer can do for
themselves: whether the code does what the spec says, and whether its
tests can fail.

## Before sending the reviewers

1. Fix the scope: a commit or a range of commits, and from it the files
   and the functions that changed. The review is of what changed and of
   what calls it and is called by it. A defect seen elsewhere goes at the
   end of a reviewer's report as seen outside the scope, unless it gives a
   wrong result, and then it is a finding like any other.
2. Find the spec item and the step of the plan the change was made from.
   A reviewer needs them to tell a defect from a decision.
3. Run the checks of the `coding` skill once and keep the output. A
   reviewer does not spend its time finding what `cargo clippy` prints.
4. Note what the reviewers cannot know from the code: an open point the
   owner has answered, a later step of the plan that will add what looks
   missing, a measurement that justified something odd. Most wrong
   findings come from missing context, and a line in the prompt avoids
   them.

## The categories

Each is described in `categories.md`, beside this file, with what to look
for and how. Send every one that applies, in parallel, in one message.

| category | applies when |
|---|---|
| `spec` | always: the code against the spec item, pyNei and the reference numbers |
| `tests` | always: whether each test can fail, and the numbers the change claims |
| `numbers` | the change has integer or float arithmetic, casts, reductions, comparisons |
| `errors` | always: panics, error types, what a user sees when it fails |
| `api` | the change adds or alters types, signatures, names, defaults or doc comments |
| `architecture` | the change touches readers, blocks, threads, the linear algebra, `cfg`, a dependency |
| `binding` | the change touches the binding crate or the Python package |

When in doubt, send it. A reviewer with nothing to report costs little.

The subagent is `code-reviewer`. Its prompt gives the category, the
commit under review, the files in scope, the path of the spec item, the
output of the checks, and the context of step 4. `spec` and `tests` build
and change code to see what happens, so they are sent with
`isolation: "worktree"`, each in its own tree, because two agents that
edit one checkout overwrite each other and then no result can be trusted.
A worktree starts on `main` and not on the commit under review: the prompt
gives the commit and tells the reviewer to check it out and confirm it
with `git rev-parse HEAD` before anything else. The other categories only
read, and share the checkout.

## Acting on the findings

When the reports are in, the orchestrator acts on them. It does not hand
the owner a list.

A finding is to be taken very seriously. The reviewer read the code
without knowing what it was meant to do, which is how its users and its
next maintainer will read it, and often it ran the case. A finding is not
set aside because it is inconvenient, because the step was nearly done, or
because the writer remembers meaning something else.

A finding is not the law either. A reviewer makes mistakes: it lacked a
piece of context, it misread a line, it assumed an input that cannot
occur, or the defect is real and the fix it suggests is worse than another
one. So each finding is evaluated, by the orchestrator, on its evidence:

1. Check it yourself. Open the line, run the case, read the part of the
   spec it cites. A finding with a failing command is settled by running
   the command.
2. Decide what it is:
   - It holds. Fix it. For a wrong result, a test that fails first, then
     the fix, as the `coding` skill says for any code. The suggested fix is
     a suggestion: the right fix is the one that fits the code, and it may
     differ.
   - It holds and the fix is the owner's: it changes the public API or a
     value a user sees, or the spec is what is wrong. It becomes an open
     point of the spec or a question to the owner, with a recommendation.
   - It holds and does not belong in this step: a GitHub issue, written
     as the `writing` skill says, so that it is not lost.
   - It does not hold. Say why, with what was checked and what the
     reviewer did not have. "I do not think so" is not a reason. A reason
     is one the reviewer would accept if it saw it.
3. When a serious finding can neither be confirmed nor refuted, send the
   reviewer the missing context and ask again, or ask the owner. A finding
   about a wrong result is never closed as not holding on a doubt.

Fix one finding at a time, with the smallest change that settles it, and
leave out the clean up that was not asked for. A fix that fails its checks
is reverted whole. When all are handled, run the checks of the `coding`
skill again, all of them.

Two or three reviewers reporting the same thing from different sides is
evidence, not repetition: handle it once and count it as more certain.

When a finding shows that the `coding` skill, the spec or the plan allowed
the defect, that document is corrected too.

## What the owner gets

A reply in chat, as CLAUDE.md describes replies: what was reviewed, what
was found that mattered and is now fixed, each finding that was not taken
with its reason in a line, what became an issue, and what needs the owner.
No report is saved. The commit message of a fix carries its why, and the
issues carry what is left.
