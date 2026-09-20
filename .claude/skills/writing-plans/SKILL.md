---
name: writing-plans
description: How an implementation plan is written in popnei. Use it when writing or revising a document under docs/plans/, which turns one or more settled specs into work packages with deliverables that can be checked, each made of tasks that a subagent can carry out. The following-plans skill is the one that executes it.
---

# Writing implementation plans

A plan turns a spec that is settled into the order of the work. It says
what is built first and what after, in which pieces, and how we will know
that each piece is done. It is carried out by the `following-plans` skill:
an orchestrator that sends each piece to a subagent, checks what comes
back, and goes on to the next one without stopping, so a plan has to be
something that can be run that way.

A plan decides nothing about the design. When writing it shows that the
spec left something open, the question goes back to the spec, as an open
point for the owner, and the plan waits or is written around it.

The prose follows the `writing` skill. A plan lives in
`docs/plans/<name>.md`. The name says what it builds, in lower case with
hyphens, `stats-exp-het`, `walking-skeleton`, and it is also the name of
the branch, `plan/<name>`, and of the report.

## No filler

The number of lines is not the measure of a plan. Filler is. A plan has
as many lines as it has things the orchestrator and the subagents need,
and none besides.

The test for a sentence: would the orchestrator, or the subagent that
gets the task, do something different without it? What fails that test is
filler, and the usual kinds are these. What repeats the spec: the
numbers, the cases and the signatures are there, the plan points at them,
and a copy is wrong as soon as the spec changes. What a task repeats from
its work package. A part that is there because the shape below has a
place for it: the parts are an order and not a form, a work package with
nothing known to go wrong has no such part, and a plan whose final check
is the sum of its work packages says so in a line.

## Before writing

- The specs the plan builds from exist and have been through their
  review. Everything a task builds has a spec behind it. When a task
  would sit under a module whose spec is not written, the row helper of
  `variant` needed by a calculation of `stats`, that part of the spec is
  written first. The same when the spec is there and too thin to build
  from, three sentences and no signature for something another module
  will call. A plan does not stand in for a spec.
- The open points of those specs are answered by the owner, or each has a
  "meanwhile" that the work can follow. For each one that is not
  answered the plan says which task the answer would change. An answer
  that comes after that task is committed becomes a new task.
- Read the specs, `docs/architecture.md`, and the code that exists. A
  plan written without looking at the code plans work that is already
  done, or builds on something that is not there.
- Run what the checks will run. A comparison with plink2 goes into a plan
  after its command has been run once, on this machine, and gave what the
  spec says. `which` tells whether a program is there.

## Work packages

A plan is split into work packages. A work package is a part of the work
that ends with something that exists and can be checked: a reader that
parses the reference VCFs, a filter that gives the same variants as
pyNei's, a wheel that installs under pyodide. When it is done the project
is in a state that works, with every check of the `coding` skill passing,
whatever comes after.

Each work package has:

- **What it gives**, in a sentence or two, in the words of a user of popnei
  or of the next work package, not of its code.
- **Its deliverables, each with the way to check it.** This is what makes
  it a work package. The check is something that can be run and that can
  fail: the tests of a part of the spec, named by that part, "every case
  of 'Missing genotypes' has a cargo test with the spec's numbers"; a
  comparison with pyNei or with a reference program, with the command; a
  Python call that a user could make; a build that produces a file. "The
  filter works" is not a deliverable. "`cargo test -p popnei filters`
  passes with the three counts of the spec's 'How it is verified' as
  literals, and the pytest test gives the same ids as pyNei at those
  thresholds" is. The numbers stay in the spec. A number that a check
  needs and the spec lacks is got by running the reference program, and
  putting it into the spec is a task.
  While the core cannot read a dataset yet, a comparison with a reference
  program is made when the literals are written. What keeps it honest is
  that the script that ran the program and its stored output are a
  deliverable too, so that anybody can get the numbers again.
- **What it stands on**: the work packages before it, and what has to be in
  place outside the plan.
- **Its tasks.**
- **What could go wrong**, when something is known: the part of the spec
  that is thinnest, the dependency that may not build for wasm, the number
  nobody has measured. The orchestrator reads this to know where to look.

Put first the work package that would change the plan if it failed. If
the wasm build of a dependency is in doubt, that is tried before the ten
calculations that would sit on it.

## Tasks

A task is a natural part of its work package, one that a person would
also name as a unit: the error type and the parser of the header, the row
helper and its tests, the Python function and its result object. It is
the unit that is given to a subagent, which starts with nothing but the
skills, the spec and the work package: the orchestrator sends it the whole
work package as the plan has it, with its task marked. So a task is a few
lines and repeats nothing of its work package. It:

- says what is built and where, the module and the files;
- names the part of the spec it is built from, by the path and the
  heading, and the pyNei function it mirrors when the spec does not;
- says which deliverables of the work package it serves, by their
  numbers. What is there when it is done is what those deliverables
  check;
- says what it needs from earlier tasks;
- is small enough to be done and reviewed in one go, and large enough to
  be worth a subagent of its own, which costs a prompt, a run of the
  checks and an entry in the report. Between an hour and a day of a
  person's work. A smaller one is part of its neighbour.

A task whose failure would be silent, a wrong number and not a crash or a
failing test, is its own task with its own commit, and says which check
guards it, so that the commit that moved a number can be found later.

Tasks that touch different files and do not need each other are marked as
able to run side by side, and so are two work packages when neither
stands on the other. The rest run in order.

## The shape of the document

1. The opening: what the plan builds and from which specs, the date, and
   the state, draft, approved by the owner, under way or done.
2. In and out: what is built, and what a reader could expect and is not,
   with where it goes. The open points of the specs that are not answered
   go here, each with the task its answer would change.
3. What has to be in place before the first task, in a form that can be
   checked, and what is not there yet: which layers exist, and so which of
   the checks of the `coding` skill run during this plan and which are
   reported as not there.
4. The work packages in order, each with the parts above. Work packages
   are numbered and tasks are numbered inside them, 2.3, with a box that
   the orchestrator ticks: `- [ ] 2.3 The row helper`.
5. How the whole plan is checked at the end, when that is more than the
   sum of its work packages: the run on the largest dataset, the build
   under wasm, a timing against pyNei.

A plan says what and in what order, not how. It does not hold code or
signatures, which are in the spec, and it does not repeat the spec's
reasons, which it points to.

## Before handing it over

- Every task names its spec item, with the path of the spec and the
  heading of the part. A task with nothing of the spec behind it is a gap
  of the spec or work nobody asked for.
- Every deliverable has a check that can fail, and its numbers are in the
  spec. A check that can be run before the code exists, the command of a
  reference program, the script that rebuilds a dataset, was run once
  while the plan was written.
- Read it for filler, sentence by sentence, with the test above.
- The sketch the `writing` skill asks for is, for a plan, the list of the
  terms and the order of the work packages with the reason for that
  order. The rest of its order is the shape above.
- Could a subagent that reads only one task, the spec and the skills do
  it? Read two tasks that way.
- Nothing in the plan decides what the spec left to the owner.
- The plan goes to the `first-reader` subagent, as any document, read as
  the orchestrator that will run it, with questions of this kind: what is
  the first thing to do, how do I know work package 2 is done, what do I
  do if the wasm build of task 1.2 fails.
- The owner approves the plan before it is run.
