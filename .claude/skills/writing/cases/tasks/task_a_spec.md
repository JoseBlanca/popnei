# Task A: a section of a spec

Write the section "What is calculated" of the popnei spec for the expected
heterozygosity per population, the plain one and the unbiased one.

Material you may read:

- `/Users/jose/devel/pynei/src/pynei/diversity.py` and
  `/Users/jose/devel/pynei/src/pynei/gt_counts.py`, the pyNei code that
  popnei has to reproduce.
- `/Users/jose/devel/popnei/docs/objectives.md` and
  `/Users/jose/devel/popnei/docs/architecture.md`.

The section has to say: what is calculated, in words and then as a
formula; what happens with missing genotypes and with populations that
have too few samples; which pyNei function it mirrors; and whether the
calculation works variant by variant or needs blocks of variants, in the
sense of `architecture.md`.

Reader: the owner of the project and, later, the person who implements it
in Rust. About one page of markdown. Write it to the output path you were
given.
