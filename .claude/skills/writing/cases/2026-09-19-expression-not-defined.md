# A symbol was defined, the expression built from it was not

2026-09-19. Trial of the first version of the skill, the same chat reply as
`2026-09-19-cost-with-two-meanings.md`.

## The text that was sent

The reply defined z, "z is the dosage matrix, variants x samples", then
used `z'z` as a row label of a table, and then wrote:

> The two give the same answers: z'z is bit exact against numpy, the
> eigenvalues agree to 1e-14.

## The owner's words

> z has been defined before, but z'z has not

## What the skill lacked

"Context before the name" listed labels, names from the code and ordinary
words with a narrow meaning. It did not list notation. Added as the fourth
kind, with the point that a label in a table is a first appearance.

## To test a change

In the text produced for the same task, every expression, `z'z`, `eigh`,
is said in words before or where it first appears, the table included.
