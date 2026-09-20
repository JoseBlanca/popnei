---
name: writing-specs
description: How a spec is written in popnei. Use it when writing or revising a document under docs/specs/, which says what one module of the core crate has to do, how pyNei does it, how the result is verified and what is left for the owner to decide. The prose follows the writing skill, which is read first.
---

# Writing specs

A spec says what one module of popnei has to do and how we will know that
it does it. It is written before the code of the module and it is what the
implementation plan and the tests are made from. It does not give the
order of the work, which is the plan's, and it does not repeat
`docs/architecture.md`, which it points to.

The prose follows the `writing` skill. Read that first. This skill says
what goes into a spec, in which order, and how much.

## The example

`example_exp_het.md`, beside this file, is the part of the `stats` spec for
the expected heterozygosity, written in September 2026 with this skill,
read by the first reader as the owner and as the implementer, and checked
by the reviewer against pyNei at commit ef0ca6e. The owner read the version
before the review and chose it as the one to imitate, and the writer then
took the findings of the review. Read it before writing a spec, for its order, for how
much it says about each thing, for how it marks what is inherited and what
is open, and for how its verification part names the command that was run
and the numbers that came out. It is one item with the parts of the module
that the item adds to. It is not the `stats` spec, which is not written.

## The two readers

The owner reads a spec to check that the right thing is going to be built
and to decide the points that are theirs. The implementer, a person or a
session of the assistant, reads it to build the module without having to
work out again what pyNei does. Every paragraph is for one of the two.
A paragraph that is for neither, that defends a choice nobody questioned
or shows how carefully the writer worked, goes out.

## How much

A spec is as long as the module is hard, and no longer. A filter that is
ten lines of pyNei takes about a page. The mixed models of the GWAS take
many. A spec that the owner does not read to the end has failed, however
correct it is, and thoroughness is the usual reason.

The test for a sentence is this: would an implementer who has read the
pyNei function and its tests get this wrong without it? Would the owner
decide differently without it? When both answers are no, it goes out. What
a threshold of NaN does to a comparison is not in a spec. That the missing
rate divides by all the samples and not by the called ones is, because it
changes which variants pass and nobody would guess it.

The parts below are the order of a spec and not a form to fill. A part
with nothing to say is a line or is left out: "It runs at the record level
and keeps nothing from one variant to the next."

A spec that describes more than one person can build in one go is split
into two specs along a line the code will also have.

## Before writing

1. Read the pyNei code that the module replaces, the functions the table
   of the architecture names and the functions they call, and read their
   tests. pyNei's results are the specification and its tests hold the
   cases that matter.
2. Write down what you find that a reader of the signature would not
   expect: a default, a NaN, a population that is dropped, a half called
   genotype that counts as a whole one, an argument used in one place and
   not in another. These are what the spec is most needed for.
3. Make the sketch the writing skill asks for.

A claim about what pyNei does comes from code that was read, and names the
function and the file. A claim about what a reference program does comes
from running it. Whether the program is on the machine is checked with
`which`, and is not taken from a note or from another session. When it
really could not be run, the spec says so: "from the documentation of
plink2 2.0, not run". A claim that was reasoned out is
either checked or left out. A number the writer worked out, a count of
the variants that pass a threshold, comes with how it was obtained, so
that the next person can get it again.

A goal of the objectives is named by what it asks, "every difference from
pyNei is written down", and not by its number.

## The parts of a spec

One spec for each module, the rows of the table in section 9 of
`docs/architecture.md`, in `docs/specs/<module>.md`.

**The opening.** What the module gives to a user of popnei, in a paragraph.
The date, that there is no code yet or what there is, the row of the
architecture table it develops, and the specs it depends on.

**Each item of the module**, a calculation, a filter, a reader or a writer,
under a `##` heading, with these parts in this order, each under a `###`
heading that names its subject. When several items share one public
function, as the per variant statistics share `calc_per_var_distribs`, the
function and what is common to them is one item, and each of the others
says only what is its own.

1. What it gives. In words first, what the number means to a geneticist or
   which variants a filter lets through, and then the formula or the rule
   with every symbol defined. A named estimator, test or correction gets
   the sentence the writing skill asks for.
2. Its Python function, with the signature and the result, the pyNei
   function it mirrors, and every difference from pyNei, because the
   objectives ask for each one to be written down. A difference
   that changes a value a user sees or the public API is an open point. A
   smaller one, the order of the populations in a result, is decided here
   with its reason.
3. The cases where the result is not what the formula or the rule would
   make a reader expect, and the cases pyNei's tests assert. Missing
   genotypes nearly always belong here. The others, a population with too
   few samples, one allele, more than two, another ploidy, no variants,
   only when something happens that the reader would not guess.
4. What pyNei does that is odd, with the function and the file, and
   whether popnei reproduces it. To reproduce pyNei is the default.
5. How it runs, when there is something to say: at the record level or the
   block level of the architecture, what is kept from one variant or block
   to the next, and the memory when it grows with the samples or the
   variants.
6. How it is verified. The reference program with its version and command,
   the dataset, the few cases whose numbers go into the tests as literals,
   and how they are compared: exactly, for counts and for sets of
   variants, or within a tolerance that has a reason, usually the digits
   the program prints. The test that runs pyNei and popnei on the same
   input. One small worked example with its numbers, which becomes the
   first cargo test. When pyNei checks the item against nothing outside
   the project, only numbers worked out by hand, the spec says so, and
   which outside program could check it is an open point, because the
   objectives ask for a reference outside the project for every
   calculation.

**The Rust interface.** The types and the signatures that another module or
the binding crate will call, in a code block, with a sentence before each
one. Nothing private.

**Speed**, when the module is one the objectives want fast, the VCF
parser, the per variant passes, the linear algebra: the dataset, the
program to compare with, and the number to reach. When there is no number
yet, the spec says that the measurement comes first.

**Open points.** What the owner has to decide, as described below.

**Not in this spec.** What a reader could expect to find here and is
somewhere else or will not be built, with where it goes. A few lines.

## Decided, inherited and open

Every statement in a spec is one of three kinds, and the reader has to be
able to tell which.

Decided: the spec states it, with the reason when there was a real choice.
When there was no choice there is no reason to give, and one is not made
up.

Inherited: pyNei does it this way and popnei reproduces it. "20 by
default, inherited from pyNei. Nobody has measured whether 20 is the right
threshold" is a complete entry, and a useful one, because it tells the
reader which numbers can move.

Open: the owner decides. The writer does not decide it for them and does
not leave it as a sentence that can be read both ways. Where it comes up,
the text gives the fact and marks it, "(**Open 2**, below)", and says no
more. The list at the end has each point once, as a request for a decision
in the form the writing skill gives: the options, what each gives and
takes, the recommendation, and what the implementer does meanwhile. The
points are numbered through the whole spec, and the list opens with a line
that says the owner decides them and the implementer follows the
"meanwhile" of each until they do. There are two sources of open points: a
difference from pyNei that the new design asks for, and an inherited
behaviour that the writer thinks the owner would change if they saw it.
The second kind is asked too, with "reproduce pyNei" as the recommendation
when there is no reason for another.

What the writer can decide alone, a name, the layout of a private struct,
is decided and not listed. Few open points, each one worth the owner's
time: what changes a result a user sees, or the public API.

While nothing has been built from a spec, a changed decision is changed in
the text and git keeps the history.

## Before handing it over

The spec, or the part of it that was written, goes to the `first-reader`
subagent as the writing skill says, once as the owner and once as the
implementer, with questions for each. For the implementer the questions
are of this kind: what comes out for a population with no called genotype,
which pyNei function holds the behaviour to reproduce, what the first test
asserts. For the owner: what is asked of me, and what does each answer
lead to.

Then read the list of open points against the text. Every **Open:** in the
text is in the list, and nothing in the list is already decided somewhere
in the text.

Read it for length as the owner would, who has a day of other work. For
each part ask what would be lost without it.

Last, the review. The code is made from the spec, so a wrong sentence in
it becomes wrong code, and the writer who misread a pyNei function will
not find the misreading by reading the spec again. The spec goes to the
`spec-reviewer` subagent, with the path of the spec and the pyNei files it
is about. It checks the claims about pyNei against the code, recomputes
the numbers, looks for what is missing, for points that should be open or
should not be, for conflicts with the architecture, and for what could go.

A change made to a spec after its review is checked against every other
place that speaks of the same quantity: search the spec for the name and
for the number. A sentence added to one part can contradict a worked
figure in another, and a spec that gives two answers is read as the one
the implementer finds first.

Evaluate each finding before acting on it. The reviewer can be wrong, and
it has not read everything the writer read. Check its evidence, take what
makes the spec more right or easier to build from, and leave the rest.
When the spec is handed to the owner, the message says which findings were
not taken and why, a line for each.

## When a spec is sent back

As in the writing skill: the spec is corrected, the principle here that
allowed the failure is revised, and the case is saved under `cases/`.
