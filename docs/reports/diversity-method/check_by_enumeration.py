"""The standardized private alleles, checked by enumerating every draw.

The closed form of `docs/specs/diversity.md` is a sum over the alleles of the
chance that the allele shows in this population's draw times the chance that
it shows in no other's. This script computes the same expectation without
that form: it enumerates every draw each population can make, takes every
combination of them, counts the private alleles of each, and averages. It is
the check that the closed form computes the quantity the spec describes in
words, which is what an outside program running the same closed form could
not tell us.

It also asks what the two disagree about, and the answer is populations that
share individuals.
"""

from fractions import Fraction
from itertools import combinations, product
from math import comb


def draws_of(counts, g):
    """Every draw of `g` alleles, as a set of allele numbers, with its weight.

    The alleles of the population are a multiset; two draws holding the same
    alleles in the same numbers are one draw, and the weight is how many ways
    it can be taken.
    """
    pool = [a for a, n in enumerate(counts) for _ in range(n)]
    seen = {}
    for draw in combinations(range(len(pool)), g):
        alleles = frozenset(pool[i] for i in draw)
        seen[alleles] = seen.get(alleles, 0) + 1
    return seen


def by_enumeration(pops_counts, g, pop):
    """The private alleles of `pop`, averaged over every combination of draws."""
    per_pop = [draws_of(counts, g) for counts in pops_counts]
    total = Fraction(0)
    ways = 1
    for counts in pops_counts:
        ways *= comb(sum(counts), g)
    for combo in product(*(d.items() for d in per_pop)):
        weight = 1
        for _, w in combo:
            weight *= w
        mine = combo[pop][0]
        others = set()
        for i, (alleles, _) in enumerate(combo):
            if i != pop:
                others |= alleles
        total += Fraction(weight * len(mine - others))
    return total / ways


def closed_form(pops_counts, g, pop):
    """The formula of the spec: in my draw, in nobody else's."""
    num_alleles = max(len(c) for c in pops_counts)
    out = Fraction(0)
    for a in range(num_alleles):
        counts = pops_counts[pop]
        if a >= len(counts) or counts[a] == 0:
            continue
        c = sum(counts)
        mine = 1 - Fraction(comb(c - counts[a], g), comb(c, g))
        term = mine
        for i, other in enumerate(pops_counts):
            if i == pop:
                continue
            c_o = sum(other)
            n_o = other[a] if a < len(other) else 0
            term *= Fraction(comb(c_o - n_o, g), comb(c_o, g))
        out += term
    return out


CASES = [
    # (name, the allele counts of each population, g)
    ("worked example, variant 1", [[3, 1], [5, 0]], 4),
    ("worked example, variant 3", [[1, 1, 1, 1], [1, 1, 1, 1]], 4),
    ("worked example, variant 5", [[4, 0], [4, 2]], 4),
    ("two populations, a draw of 2", [[3, 1], [2, 2]], 2),
    ("two populations, three alleles", [[2, 1, 1], [3, 0, 1]], 2),
    ("three populations", [[2, 2], [3, 1], [1, 3]], 2),
    ("three populations, four alleles", [[2, 1, 1, 0], [1, 1, 0, 2], [2, 0, 1, 1]], 2),
    ("a population holding one allele", [[4, 0], [2, 2]], 3),
]

print("the standardized private alleles: the closed form against every draw\n")
worst = 0
for name, pops_counts, g in CASES:
    for pop in range(len(pops_counts)):
        a = closed_form(pops_counts, g, pop)
        b = by_enumeration(pops_counts, g, pop)
        worst = max(worst, abs(float(a - b)))
        print(f"  {name:38} pop {pop}  closed {float(a):.10f}  "
              f"enumerated {float(b):.10f}  {'equal' if a == b else 'DIFFER'}")
print(f"\nthe largest difference over {sum(len(c) for _, c, _ in CASES)} "
      f"population and case pairs: {worst}")

print("\nwhat the two disagree about: populations that share individuals.")
print("Two populations of one diploid individual each, the same individual,")
print("genotype 0/1, and a draw of 1 allele. The closed form treats the two")
print("draws as independent; enumerated over the one shared genotype they are")
print("the same draw, so no allele can be private to either.")
shared = [[1, 1], [1, 1]]
print(f"  closed form, pop 0:  {float(closed_form(shared, 1, 0)):.10f}")
print("  what a shared individual gives: 0.0, the two draws being one draw")
