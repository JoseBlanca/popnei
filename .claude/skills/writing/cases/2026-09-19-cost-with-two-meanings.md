# One word with two meanings in a mirrored sentence

2026-09-19. Trial of the first version of the skill.

## The situation

The owner had asked whether popnei should use faer everywhere, so that
there is one linear algebra implementation, or BLAS natively and faer only
in wasm. The task and the notes the writer had are in
`tasks/task_b_chat.md` beside this file: faer on one
thread takes 34 s for the kinship product of 1e6 variants x 1000 samples
where Accelerate takes 2 s, and two backends mean 4 or 5 functions written
twice.

## The text that was sent

The opening of the reply:

> BLAS natively and faer in wasm. One code path costs 34 s instead of 2 s
> on the kinship product of a dataset in the range the objectives cover,
> and two code paths cost 4 or 5 functions written twice.

## The owner's words

> This sentence is written to be clever, but it is using "cost" to mean
> two different things that are not comparable. It would be clearer
> something like: "One code path costs 34 s on the kinship product of a
> dataset in the range the objectives cover, while two code paths would
> cost 2 s, but drawback is that the two code paths add some complexity: 4
> or 5 functions written twice".

## What the skill lacked

It asked for both sides of a comparison in the same units, which does not
cover a trade-off between two kinds of thing, and it had nothing against a
sentence shaped for symmetry. Added: the trade-off paragraph under "A
number comes with what it was measured on", and the first paragraph of
"Plain form".

## To test a change

Give a writer the task and check the sentence that states the trade-off:
each option has its running time, the added code is named as code, and no
word carries two meanings.
