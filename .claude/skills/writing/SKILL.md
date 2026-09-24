---
name: writing
description: How prose is written in popnei. Use it before writing or revising anything a person will read, a document under docs/, a spec, an implementation plan, a README, a doc comment, a commit message, a GitHub issue or a pull request. Chat replies follow the same principles, and the ones that matter most there are in CLAUDE.md.
---

# Writing

## The reader

What is written in popnei is read by one of two people.

The first is the owner of the project, a population geneticist who
programs in Python and Rust. They know the genetics, pyNei, and the
decisions recorded in `docs/`. Of statistics they know the foundations:
what a test is, what a likelihood is, how the classical and the Bayesian
approaches differ, and that a p-value is about the probability of the
evidence and not of the model. They do not know the catalogue: the named
tests and estimators, what each one assumes, and when one is preferred to
another. So a score test, a Wald test, REML or a named correction gets a
sentence that says what it does and why it is the one used here, before
it is relied on. The genetics can be brief.

They were not in the session: they did not see the files that were
opened, the commands that were run, or the names that came up along the
way. So they need a bit of context before the content: what was being
worked on, what the question was, and where it sits in the project. A
sentence or two is usually enough.

The second reader arrives later, a contributor or another session of the
assistant. They know the field and have only the page.

Both read in order to do something: to decide, to build, or to check. A
text is good when its reader gets through it once, without asking a
question and without opening anything else, and can then do what they
came to do.

`docs/objectives.md` and `docs/rust_core.md` are the voice to match. Read
a section of `rust_core.md` before writing a document: facts and numbers,
each measurement with its dataset, headings that name their subject,
almost no bold, and no sentence about the text itself.

## Before writing: the sketch

The order of a text is decided before its sentences. A text whose order is
found while writing uses its terms before it defines them and gives its
reasons before its subject. So write a sketch first, a few lines that are
not part of the text:

1. The subject in one sentence, with each option or thing named by what it
   is: "whether popnei uses faer for all its linear algebra, or faer in
   wasm and the system BLAS natively".
2. What the reader should know or be able to do when they finish, in the
   words you would use with a colleague across the table. When that cannot
   be said yet, the problem is in the thinking and no sentence will fix
   it.
3. The terms the text will need, each with its definition, in the order
   they will be used: "one code path: faer everywhere".
4. The points, one line each, in the order the reader needs them. A point
   that leans on another comes after it.
5. What is left out, and where it goes.

Then write from the sketch. The opening of the text is the first two lines
of the sketch, in full sentences. When the text is written, find the first
use of each term in line 3 and check that its definition comes before it,
as a sentence of the text. A meaning that the reader could work out from
what came before does not count as a definition.

## The principles

### Write the contents, not the name of the category

"The module holds the filtering policy" gives the reader nothing to use.
"The module holds the four filters, their thresholds and the name of each
reason to drop a variant" does. Policy, machinery, infrastructure, logic,
semantics and surface are names of bags. Write what is in the bag. When
the list cannot be written, the thing is not yet understood.

A problem is written as what can go wrong, not as the rule it breaks:
"popnei already runs the variants on a pool of threads, multi threaded
faer would start its own pool inside it, and how the two behave together
has not been tested". That says what would happen and how much is known.
"It nests rayon pools, against the design rule" says neither.

An adjective in the place of a number is the same fault: fast, large,
significant, negligible, most. The repair is to put the fact in, "The Rust
parser reads the 400 MB VCF, 100000 variants x 1000 samples, in 0.55 s on
one thread. pyNei takes 13.5 s." Taking the adjective out and leaving the
rest leaves a sentence that still says nothing.

### A number comes with what it was measured on

The dataset, the machine, the threads, the program that was run, and the
date when the measurement is not from the work being reported. The date is
the one the session was told, not one worked out from the dates already in
the file: on 24 September 2026 sixteen dates one day ahead were found across
two suites, a package and two specs, written by sessions that had each
inferred the same wrong day. A writer who does not know it asks. A
comparison has both of its sides in the same units, and says which side
is the better one when the reader could doubt it. A number carried over
from another document keeps all of this.

A trade-off gives up one kind of thing to get another, running time for
code to maintain, memory for speed. Each kind is named, and each option
gets its value in both: "With faer everywhere the kinship product of 1e6
variants x 1000 samples takes 34 s on one thread. With BLAS natively it
takes 2 s. What BLAS natively adds is code: 4 or 5 linear algebra
functions written twice, once for each library." The two kinds are laid
side by side for the reader to weigh. A sentence that mirrors them around
one word, "one path costs 34 s and two paths cost 5 functions", makes
"costs" mean seconds in one half and code in the other. A word keeps one
meaning through a sentence.

### Context before the name

A term that does work in a sentence is explained before that sentence or
inside it. The reader lacks four kinds of names:

- Labels made up during the session: step numbers, codes of findings, the
  nickname of a script. Most name a thing by the part it played in the
  work, "the baseline run". The reader needs what the thing is, "numpy
  linked to Accelerate".
- Names from the code. That a type is called `Needs` does not make it a
  word the reader has. Say what it is, which fields a calculation asks the
  reader of the file to fill, and then the name can be used.
- Ordinary words that mean something narrower here. When a document leans
  on the difference between a variant and a block of variants, both are
  defined first.
- Notation. To define `z` is not to define `z'z`. An expression is said
  in words where it first appears, with what it is for: "z'z, the product
  of the transpose of z with z, a samples x samples matrix from which the
  kinship is calculated". A label in a table is a first appearance too, so
  the words go in the text before the table.

One name for each thing. Before a new name is added, count the ones the
document already uses for it. For the things of the domain, a population,
a block, a half called genotype, the name is the one `docs/glossary.md`
gives. The glossary chooses the word and the text still explains it.

### What the reader came for goes first

In a document, what it is about and what it decides. In a section, what
the thing does, in words, before the formula, the type or the signature.
In an issue, the finding. In a reply, the answer to the question that was
asked. The work usually happened in the opposite order, the answer was
found last, and the temptation is to tell it in the order it happened.

Order carries importance. What matters most is first or has a paragraph of
its own, and then no sentence needs to say that it matters.

### Only what the reader can use

For each sentence, what can the reader do with it? These go out:

- The story of how the work went.
- Sentences about the text: "in one paragraph", "briefly", "this section
  describes".
- Sentences about the reader's reaction: "surprisingly", "notably", "the
  key point is".
- Answers to an objection that nobody raised: "this is not a matter of
  taste".
- Sentences that only give a verdict, "this is not free either", "there
  is a catch", and leave what is good or bad to the next sentence. The
  sentence that holds the fact carries the verdict too.

The test is to delete the clause. When no fact is lost, it stays deleted.

Two things stay although they look like the story of the work: a trap the
next person would fall into, and a measurement that closed an option. "Two
traps found on the way" and "Measured and not taken", in `rust_core.md`,
are both. They stay where the next person will look for them, the document
or the issue about that subject. A trap of the wasm build does not belong
in a reply about which linear algebra library to use.

### Everything at its true strength

A choice is written as a choice, with the goals it serves and what would
have to be true for another option to win. It is not "forced" or "the only
way". A design is explained by its reason, and what is wrong with the
alternative is not yet a reason.

A constraint is not a goal. "A significant slowdown would make us
reconsider" and "the speed must not change" are different claims, and only
one of them is usually true. A standard stricter than the real one reads
as care and leaves the reader with a criterion nobody should hold.

What was measured is kept apart from what was assumed, read somewhere or
taken from another project's benchmark: "that has not been measured here".
A claim that holds only in part is given with the condition under which it
holds, "up to 2000 samples", "on Apple hardware". Somewhat, relatively and
in general tell the reader that there is a limit and hide where it is.

### Plain sentences that stand on their own

Short sentences with active verbs and ordinary words. A sentence is not
shaped to sound balanced or sharp. The plain version is usually longer,
and it is the one to write.

Every count has its noun and its set: "six of the nine filters", never
"six of the nine". Every pronoun has an antecedent the reader can see, and
a "this" that could point at two things gets its noun: "this filter".
After an edit, read the sentence before and the sentence after, because a
deletion can take an antecedent with it and the writer, who still
remembers it, will not notice.

A list for parallel items, a table for measurements, prose for reasoning.
A heading names the subject of its section. Bold marks the terms that a
list defines and nothing else. Importance is shown by position.

## The forms

- **Documents under `docs/`.** The opening paragraph says what the
  document is, its date, what it decides and where the related documents
  are.
- **Doc comments.** What the item is or does in the words of the domain,
  the units and the shape of each value, and what a caller must know, that
  a missing allele is -1. How it is implemented stays out unless the
  caller sees it.
- **Commit messages.** A lower case subject that says what changed, and a
  body with the why and the numbers. Nothing about the session.
- **GitHub issues.** The title states the finding or the task. The body
  says what was seen, on what data and how to see it again, what it means
  for popnei, and what is proposed or asked. It is read months later by
  somebody who has nothing else.
- **A recommendation**, in a reply or in a document, in this order: what is
  being decided, with the options spelled out; the recommendation, in a
  full sentence; the reasons, each with its numbers; what is not known;
  what is asked of the reader.
- **A request for a decision**, in any form. The options, what each one
  gives and takes, the recommendation, and what happens next in each case.
  It is ready when the reader can answer without asking anything back.
  What the writer can decide alone is decided and not asked.

## Before handing a text over

The writer cannot see what is missing from the page, because they know it.
So a document or an issue goes to the `first-reader` subagent, which gets
only the text, one sentence that says who reads it, what for and which
documents they already know, and three to five questions that the reader
should be able to answer once they have read it. It returns what it
understood the text to say and to ask, its answers to the questions, the
names it did not have, and the sentences it could not follow.

Compare its summary and its answers with what was meant. Where they differ
the text is wrong, not the reader. Fix what it reports, and send the text
again when the fix was large.

A commit message or a doc comment does not go to the subagent. Read it
once more against "Only what the reader can use" and "Context before the
name".

## When a text is sent back

The text is corrected, and so is this skill. Find the principle that
allowed the failure, or that is missing, and revise it where it stands.
The case, the situation, the text that was sent and the owner's words, is
saved under `cases/` beside this file. The cases are for trying out
changes to the skill and are not read when writing.

The skill does not grow by one entry for each failure. A correction
becomes a new paragraph only when no principle covers it, each principle
keeps one example, and the sentence that was sent back goes to the case
and not here. About 200 lines is the size to stay near, a figure to get
oriented by and not a limit.
