---
name: following-plans
description: How an approved implementation plan of popnei is carried out. Use it when the owner asks to run, execute, implement or continue a plan under docs/plans/. The session is an orchestrator: it works in a worktree and a branch of the plan, sends each task or work package to a subagent, checks what comes back, has each work package reviewed, keeps a work report, and goes on until the plan is done or the owner is needed. It never merges into main.
---

# Following an implementation plan

The session that runs a plan is an orchestrator. It does not write the
code. It sends each piece of the work to a subagent, checks what comes
back, and decides what happens next. The reason is its context: a plan
takes many hours of work, the details of each task would fill the
orchestrator's context long before the end, and an orchestrator that has
forgotten the start of the plan decides badly at its end. What it keeps
is the plan, the state of the work, the decisions and the report.

In ideal conditions it goes from one work package to the next until the
plan is done, without stopping. It stops when it needs the owner, and
only then.

## Before the first task

1. The plan is approved by the owner, and its specs are settled. If not,
   say so and stop.
2. The worktree and the branch. Every plan is carried out in its own,
   named after the plan:

       git worktree add .claude/worktrees/<plan> -b plan/<plan> main

   from the local `main`, so that commits not yet pushed are in it, and
   then `EnterWorktree` with that path. All the work of the plan happens
   there, and `main` is not touched.
3. Check what the plan says has to be in place, by running it and not by
   reading it. Whether a program is installed is checked with `which`.
   When something is missing, the plan cannot start: tell the owner what
   is missing and what would put it in place, and stop. Do not build it
   on the side, because it was somebody's decision that it is not part of
   this plan.
4. Read the board, `.claude/board/`, described in the README there: the
   messages of the sessions working on other branches. One that names a
   file this plan will also change, or a finding the plan rests on, goes
   into the report, and a file that two branches will change goes to the
   owner before the first task touches it. Then post the start message
   of this plan: the plan, the branch, and the shared files it may
   change.
5. Start the work report, described below, and commit it with the plan's
   state changed to under way.

## One task

Send the task to a subagent, `general-purpose` on the model the owner
uses for coding, with a prompt that stands on its own, because the
subagent knows nothing of the session:

- the path of the worktree, and that all its work happens there, on the
  branch that is checked out. It is not given `isolation: "worktree"`,
  which would start it from `main` in a tree of its own;
- the path of the plan in the worktree, the number of its task, and that
  it reads the whole work package of that task there, so that it has the
  deliverables its task serves and sees what the tasks beside it do and
  it must not. The work package is not copied into the prompt. The
  subagent can read the file, and a copy per task fills the context of
  the orchestrator with its own plan, many times over;
- the path of the spec, of `docs/architecture.md` and of
  `docs/glossary.md`, which it reads itself;
- the skills to follow: `coding`, and `writing` for its doc comments and
  its commit message;
- today's date, and that every date it writes is that one unless it is
  somebody else's earlier measurement. A subagent has no clock and works
  the date out from the dates already in the files, so a file that is a
  day ahead makes the next writer a day ahead. On 24 September 2026 two
  sessions carrying out two plans wrote 28 dates that were a day in the
  future, every one of them ahead and none behind, and 14 of those came
  from one subagent working through one list of fixes: it worked the date
  out once and was wrong 14 times. Nine of the 28 dated a measurement,
  where the cost is worst, because a reader checking whether a bound still
  holds meets a date that has not happened and doubts the number rather
  than the date;
- what is already there from earlier tasks that it builds on;
- what the orchestrator knows and the code does not show: an open point
  the owner answered, a decision taken two tasks ago;
- to commit its work when the checks pass, one commit for the task, and
  not to touch the plan or the report. When the task made it add
  something to the spec, that is a commit of its own before the commit of
  the code;
- what to send back, in under 300 words: what it built, the commit, the
  last line of each check, what it did differently from the task and why,
  anything it found that the orchestrator or the owner should know, and
  any question it could not settle.

Tasks, and work packages, that the plan marks as able to run side by side
can go to several subagents at once only if they touch different files.
One tree has one writer per file.

Several small tasks of one work package can go to one subagent in one
prompt, and a work package that is small can go whole. A task never goes
to a subagent together with a task of another work package.

## Checking what comes back

What a subagent says is a claim. A claim that the next step rests on is
checked. This costs the orchestrator little context, because it looks at
results and not at code:

- `git log` and `git status` in the worktree: the commit is there and the
  tree is clean.
- The checks of the `coding` skill, run by the orchestrator, reading only
  their last lines.
- What the task said would be there when it was done: the named tests
  exist and pass, the function can be called. A cargo command that ran 0
  tests exits with 0, so read the counts in its `test result` line and
  not only that it passed.
- What it did differently from the task: is it small, or does it change
  what a later task builds on?

When the work is not right, the same subagent gets it back with what is
wrong, through `SendMessage`, since it has the context. When it fails a
second time, a fresh subagent gets the task with what was learned. A
third failure is a reason to stop and ask.

Then tick the task in the plan, note it in the report, and commit both.

## The end of a work package

1. Run each deliverable's check yourself, as the plan gives it, and put
   the command and the result in the report. A deliverable that does not
   pass is not done, whatever the tasks said.
2. Have the work package reviewed, as the `code-review` skill says, over
   the commits of the work package. The `spec` and `tests` reviewers work
   in worktrees of their own, at the last commit of the work package. The
   review follows the size of the work: a work package of an afternoon is
   reviewed with the categories that apply and no more, or together with
   the next one when the two are one piece of code. `spec` and `tests`
   always run, and nothing reaches the end of the plan without a review.
3. Evaluate the findings yourself, as that skill says: this is a
   judgement and it is the orchestrator's. Reading the reports is the
   context it is there to spend. The fixing is delegated like any task: the
   findings that hold go to a subagent, the one that wrote the code when
   it is still there, with the orchestrator's verdict on each, and it
   fixes them test first. It may answer that a finding does not hold, with
   evidence, and the orchestrator weighs that as it weighs a reviewer.
4. Run the checks and the deliverables again after the fixes.
5. Write the work package into the report and go on to the next one.
6. Read the messages on the board that are newer than the last read, and
   post one for each thing this work package learned that another session
   would otherwise find again, and for each change to a shared file.

## Changing the plan

A plan is written before the work and the work knows more. The
orchestrator changes the plan when the work asks for it, without
stopping for the owner, when the change is small:

- a task split in two, two merged, tasks put in another order;
- a task added because something needed was not foreseen, inside what the
  work package gives;
- a file or a module other than the one the task named;
- a deliverable's check replaced by another that is at least as strong,
  because the first could not be run.

Every change goes into the plan, where git keeps the old text, and into
the report with its reason.

It stops for the owner when the change is not small:

- a deliverable dropped, or its check made weaker;
- a work package added or removed, or the plan growing beyond what it
  said was in;
- anything that changes a spec: a value a user sees, the public API, an
  open point. The orchestrator does not answer an open point for the
  owner, it follows its "meanwhile";
- the spec and the code, or the spec and pyNei, contradict each other in
  a way that makes a task's approach wrong;
- a finding of a review whose fix is the owner's, when the rest of the
  plan rests on it;
- a task that failed three times;
- anything that cannot be undone, outside the worktree.

Before stopping, do the work that does not depend on the answer: another
task, another work package that stands on its own. When stopping, leave
the tree with the checks passing and everything committed, bring the
report up to date, and ask as the `writing` skill says a decision is
asked for: the options, what each gives and takes, a recommendation.

## The work report

`docs/reports/<plan>.md`, in the worktree, written while the work goes
and not at the end, because at the end the orchestrator no longer
remembers the second work package. It follows the `writing` skill: it is
for the owner, who was not there.

After each work package, a short section:

- Did it finish as planned? The deliverables, each with the command that
  checked it and what it gave.
- What was changed in the plan, and why.
- What the review found that mattered, what was fixed, and which findings
  were not taken, with the reason in a line.
- What the owner should know: something learned about pyNei, a number
  that surprised, a part of the spec that proved thin, a risk for what
  comes next.
- How the work went. This one is not for the owner but for whoever next
  revises a skill or writes a plan, which the `writing` skill allows when
  the heading says so and the owner is told there that they can stop. What
  goes in it is what would change a skill, the shape of a plan or the size
  of a task: a task that had to be sent twice and what made it fail the
  first time, a skill that was unclear or wrong with the sentence that
  misled, a way of working that cost something. Not a log of what happened.
  Every number in it carries what it is measured against, as every number
  in popnei does. The tokens a task's subagent used, which the result of
  the `Agent` tool gives, say nothing on their own; set beside the tokens
  the review and the fixes of that same task cost, they say what a task of
  that size really costs, and that is what sizes the tasks of the next
  plan.

When the plan is done, the top of the report gets what the owner reads
first: whether the plan is done, what exists now that did not, what is
left open, and what is asked of them, which is at least the merge. Every
number in the report comes from a command that was run, and the report
goes to the `first-reader` before the owner is told. It is briefed with
the owner as its reader and told that the last section is for somebody
else, so that what comes back is whether that signpost works, and not
whether a section serves a reader it was never written for.

## When the plan is done

1. Every task ticked, every deliverable checked, the plan's final check
   run, the state of the plan changed to done, the report finished, all
   committed on the branch.
2. Tell the owner, in a reply as CLAUDE.md describes them: the plan is
   done, on branch `plan/<plan>`, with the report's path and what it asks
   of them.
3. The orchestrator does not merge into `main` and does not push to it.
   The order to merge is the owner's. When they give it, merge as they
   say, and run the checks on `main` after.
4. After the merge, clean up: `ExitWorktree` with `keep`, then
   `git worktree remove .claude/worktrees/<plan>` and
   `git branch -d plan/<plan>`, and the worktrees of the reviewers, and
   remove the plan's messages from the board. When
   something on the branch was not merged and the owner has to judge
   whether to keep it, an experiment, a measurement, a change they turned
   down for now, nothing is removed: say what it is and where, and leave
   it to them.

## When a session ends before the plan

A plan can outlast a session. The state is in the repository, not in the
session: the ticks in the plan, the report, the commits on the branch.
A new session enters the worktree with `EnterWorktree` and its path, reads
the plan and the report, runs the checks, and goes on from the first task
that is not ticked.
