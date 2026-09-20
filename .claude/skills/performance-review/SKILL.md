---
name: performance-review
description: How the speed and the memory of popnei's code are reviewed and improved. Use it when the owner asks for a performance review, a profile, where the time goes, whether something is on the hot path, how popnei compares with pyNei, plink2 or GMMAT, or to make something faster; and when a spec or a plan has a speed target to meet. The session sends one reviewer per category, writes a report with a measurement plan, and then runs the experiments one at a time, each delegated, keeping only what a measurement confirms.
---

# Performance review

This skill is ported from the performance review of pop_var_caller, which
worked well there, with what popnei adds: a Python boundary, two linear
algebra backends, a wasm build, and numbers that have to stay right.

A performance review looks for places where a gain that can be measured
is plausible, and proposes the experiment that would confirm or refute
it. It does not rewrite the code and it does not claim gains. The prior is
hostile: most obvious optimizations do not matter, because the code is
cold, or do not help, because the compiler already does it, or trade real
complexity for an imagined gain. Each candidate has to answer three
questions. Where on the call graph does this run? What would we measure?
How much complexity does the change add?

Speed is the fourth goal of `docs/objectives.md`, after right, usable
from Python and in the browser, and it is stated there as something that
is measured: on stated datasets, against pyNei and against plink2 and
GMMAT, before and after every change that claims it.
`docs/rust_core.md` holds the measurements the project started from.

The session is the orchestrator. It works out which categories apply,
sends one `perf-reviewer` subagent for each, puts their findings into one
report, and then runs the experiments.

## The principles

- **Profile first, match patterns second.** A finding with no profile, no
  benchmark and no strong argument from the call graph that the code is
  hot is filed low, with the measurement that would promote it. No rewrite
  is proposed for code that cannot be shown to be hot.
- **Every candidate has a measurement plan.** Which benchmark, which
  profiler, which metric, wall time, allocations, bytes copied, lock wait,
  syscalls, peak memory, and which threshold makes the change worth
  keeping. "It will be faster" is not a finding.
- **Complexity is a cost.** Every change names what it adds, a type, a
  lifetime, an `#[expect]`, a build flag, a dependency, and weighs it
  against the expected gain. A change that doubles the maintenance for 2%
  of wall time is a bad trade.
- **Hot paths only.** An optimization of code that runs once, at start
  up, in the handling of an error, is noise. Say what the call frequency
  is assumed to be, and lower the severity when it was not verified.
- **Right before fast.** No optimization weakens an invariant, adds
  `unsafe` or silences a lint of the `coding` skill without a review of
  its own. An optimization keeps every result within the tolerance of its
  spec and the same at every number of threads. `checklists/numbers.md`
  holds the rules, and every proposal that touches arithmetic says what it
  does to the numbers.
- **One change per measurement.** A change of allocator, LTO and a
  refactor measured together give a result nobody can read. Each
  experiment names the one hypothesis it tests.
- **Timings are taken one at a time.** Two builds or two benchmarks
  running at once on this machine measure each other. Reviewers read and
  count in parallel. Wall times are taken by one agent at a time, with
  nothing else building.

The severities and the format of a finding are in `finding_format.md`.
The orchestrator and every reviewer read it at the start.

## The procedure

### 0. What is already known

Find what exists for the code in scope, so that the review starts where
the last one ended: a report of an earlier performance review,
`docs/reports/perf-*.md`, whose open findings and measurement plan say
which experiments were run and what they showed; the spec, which may hold
a speed target; the plan and its work report; the issues; and
`docs/rust_core.md` for the baselines.

### 1. The scope and the call frequency

A crate, a module, a range of commits, one function under a benchmark, or
one Python call as the user makes it. For each file in scope say which
functions are believed to be hot and on what evidence: a profile, a
benchmark, an argument from the call graph, or none yet. Code with no
evidence is reviewed at lower priority, and every reviewer is told so.

### 2. What measurement there is

List every benchmark under `benches/`, every profile that was kept, and
every measurement the owner gave. Say what is missing. When there is no
measurement of the code in scope, the first thing the review delivers is
a measurement plan and not changes to the code, and the verdict says so.
Output of a real benchmark or profile is quoted as printed and given to
every reviewer, so that they do not run it again.

Read `profiling_environment.md` here, before any measurement plan is
written. It says which tools run on this machine, with invocations that
work, and what the machine cannot measure. A plan that names a tool that
is not there cannot be run.

A sampling profile is what the rest stands on, and nothing replaces it.
DHAT sees allocations only, and a function can take 30% of the CPU and
allocate nothing. `cargo asm` shows whether a loop was vectorized and not
whether the loop matters. A criterion benchmark gives a wall time whose
change from run to run is often larger than the effect of one change:
the load of the system moves the median 20 to 25% between runs while the
interval inside one run reads 1%. None of them says which line holds 30%
of the CPU. A review without a profile gives a long list of candidates
matched by pattern, most of which show no gain when applied, because the
site was not in the top of the self time. In pop_var_caller a first round
applied seven such findings and the benchmark could not see a change; the
second round, with a profile, got gains on the same code.

So when no sampling profile can be taken, say so to the owner before
going on: what is blocked and why, what would fix it, what it costs not
to fix it, and ask whether to fix it first. Do not route around it in
silence. If the owner accepts going on without one, no finding is above
Likely, and the report says so at its top.

When the scope is a pipeline, a reader, a read ahead thread, a consumer,
a writer, add one measure of time off the CPU to the profile before
trusting its ranking: how long each stage waits on the one before or
after it. A queue that is always full in front of one stage names the
bottleneck more directly than any ranking of self time.

When the scope is a Python call, take the profile of the Python process,
`sample <pid>`, which shows how the time splits between Rust, numpy and
pandas, before looking at any Rust.

### 3. The intent and the targets

Before judging whether the code is fast enough, write a paragraph on what
it is for: its input sizes, from the objectives, tens of samples to ten
thousand and tens of thousands of variants to a million; its target, a
number from the spec or the program to compare with, pyNei, plink2,
GMMAT; its hardware, this Mac natively and a browser tab in wasm; and its
memory bound. When no target is stated, ask, or file the gap as a Note.
That paragraph goes into every reviewer's prompt.

### 4. The categories

Each is a checklist under `checklists/`.

| category | applies when |
|---|---|
| `methodology` | always: that the benchmarks, the profiles and the build configuration are sound before any finding about code is acted on |
| `numbers` | any finding or proposal that touches arithmetic, a reduction, a block size, a maths function, an integer type, the build profile |
| `allocations` | hot code builds a `Vec`, a `String`, a `Box`, a map, clones owned data, formats strings, or grows a buffer without bound |
| `data_layout` | hot code walks collections of structs, shares atomics between threads, or there is a suspicion of cache misses or false sharing |
| `concurrency` | the code uses rayon, threads, locks, atomics, channels, or calls BLAS from a thread |
| `hot_loops` | tight loops over genotypes or bytes, float reductions, filters that branch on the data, slice indexing, iterator chains |
| `io_and_syscalls` | the VCF reader or writer, the vars file, gzip, any per block I/O |
| `linalg_and_wasm` | the `linalg` module or what calls it, the size or layout of a block, the wasm build |
| `python_boundary` | the binding crate, the Python package, or a timing taken from Python |

When in doubt, send it. A reviewer with nothing to report costs little.

### 5. Sending the reviewers

One `perf-reviewer` for each category, in parallel, in one message. Each
gets: its category; the commit under review and the files in scope with
the call frequency of step 1; the paragraph of step 3; the output of the
benchmarks and profiles of step 2, as printed; and what the code does not
show, as in the `code-review` skill. A reviewer that builds, runs or
changes code gets `isolation: "worktree"` and is told which commit to
check out first, because a worktree starts on `main`. Reviewers count and
read in parallel. A reviewer that needs a wall time says so in its
finding's plan and does not take it while the others run.

### 6. Collecting

Read every report. Promote the findings on which several categories
agree that the same site is hot. Demote the ones that turn out to be on a
cold path. Pass on the lines "for another category". A reviewer that
skipped its scope or reported what it cannot show gets its category sent
again, with the gap named.

### 7. The report

`docs/reports/perf-<scope>-<date>.md`, written as the `writing` skill
says, for the owner, and through the `first-reader` before they are told.
It is kept, unlike the report of a code review, because the next review
starts from its measurement plan and its numbers. In this order:

1. The scope and its limits: what was reviewed, at which commit, the
   targets and the sizes, the hot path evidence there was, what was left
   out and why, the categories sent.
2. The verdict, one of: profile first, there is not enough evidence to
   recommend changes and the next part is what the review delivers; run
   the experiments, there are candidates and their order is below; apply,
   at least one candidate has a matching profile, a plausible mechanism
   and a contained cost.
3. The measurement plan: the benchmarks and profiles to add or run, in the
   order in which they unblock the findings, each with its command, its
   dataset, what it should print and which threshold answers which
   question.
4. The build configuration: LTO, codegen units, the panic strategy, debug
   info for profiling, the allocator, `target-cpu`. Each is its own
   experiment, and what is set for the native build is looked at again for
   the wasm wheel.
5. The findings about code, by severity and then by confidence, each in
   the format of `finding_format.md`, numbered H1, L1, S1 so that they can
   be referred to. In the prose of the report a finding is named by what
   it is, with its number after it.
6. Seen outside the scope.
7. What the code already does well, up to three patterns worth copying,
   each with its file, and none when there is nothing specific.

## Running the experiments

A finding of a performance review is a hypothesis, and what decides it is
a measurement. The reviewer's judgement does not decide it and neither
does the orchestrator's. That is the difference from a code review, where
a finding is evaluated on its evidence. Here the evidence is produced.

The orchestrator takes the findings in the order of the report and, for
each one that is worth its experiment, delegates it to a subagent, in a
worktree and a branch for the review, as the `following-plans` skill
delegates a task. One experiment at a time, because of the timings. The
subagent:

1. takes the baseline first, on the commit before its change, with the
   command of the plan, on this machine, and keeps the output;
2. makes the one change;
3. runs the correctness checks of the `coding` skill, before looking at
   any timing;
4. measures again, the count first when the mechanism can be counted, and
   then the wall time;
5. when the difference in wall time is what decides and is under about
   5%, confirms it with a second clean measurement of each side, or by
   reverting the change and measuring again, so that an unlucky baseline
   does not decide;
6. reports the two outputs as printed, and what the change cost in
   complexity.

The orchestrator then closes the finding with one of these, in the
report: applied, with the commit and the two numbers; the experiment
showed no gain, closed, which is an expected and welcome outcome and is
what the plan is for; not run, with the reason, the site turned out cold,
the complexity was not worth the most the gain could be; for the owner,
when the change would alter a result, the public API, or the precision of
a number; or an issue, when it is real and belongs to later work.

A change that is kept gets a commit of its own whose message has the two
numbers, the dataset and the machine, as the objectives ask. Nothing is
merged into `main` without the owner's order.

## What the owner gets

A reply in chat as CLAUDE.md describes them: what was reviewed and
against which target, where the time goes, what was tried and what each
experiment gave in numbers with their dataset and machine, what was kept,
what showed no gain, what is theirs to decide, and the path of the
report.
