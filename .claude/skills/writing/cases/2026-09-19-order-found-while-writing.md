# A reply whose order was found while writing

2026-09-19. Third round of the trial, the chat reply of
`tasks/task_b_chat.md`, written with the merged skill and the first chat
rules.

## The text that was sent

The opening, and the sentence of the second paragraph that explains it:

> BLAS natively and faer only in wasm. The one code path does cost: on the
> Mac the kinship of a 1e6 variants x 1000 samples dataset goes from 2 s to
> 34 s on one thread, and that is the calculation the library is built
> around.
>
> [...] One side is numpy linked to Accelerate, which is the system BLAS
> here; the other is faer, pure Rust, called through a small pyo3 crate
> written for this.

## The owner's words

> I like it less than before. the main problem is that it is not well
> organized. It starts with "BLAS natively and faer only in wasm". That's
> not even a sentence it would be clearer to state clearly what is this
> about: "We're evaluating the use of either just faer or both faer in
> wasm and BLAS natively".
>
> it is talking about one code path without explaining before that one
> code path means only using faer. It is explained later.
>
> maybe the agent would improve the writting if it thought before starting
> writting on a sketch of how to structure the text for clarity and only
> then writting.

## What the skill lacked

"Context before the name" was there and the writer knew it, but nothing
made the writer decide the order of the text before writing it. The chat
rule "the first sentence answers the question" produced a fragment as the
opening, and it was wrong to hold that the owner's question is context
enough. Added: the sketch under "Before writing", the order of a
recommendation under "The forms", and the first chat rule in CLAUDE.md
rewritten.

## To test a change

In the reply produced for the same task: the first sentence says what is
being decided with both options spelled out, the recommendation is a full
sentence, and "one code path", z'z and every other term is defined before
its first use.
