"""It writes the standardized private alleles of 18 pairs of a case and a
population, each computed twice, and writes nothing if the two ways disagree
anywhere.

The standardized private alleles of a population are, at one variant, the
alleles expected to be in a draw of `num_called_alleles` of that population's
called alleles and in no other population's draw of the same size. "The
private alleles" of docs/specs/diversity.md computes them with the closed form
of Kalinowski (2004),

    E = sum over a of [ q_p(a) * product over r != p of (1 - q_r(a)) ],

where `q_r(a)` is the chance that allele `a` shows in the draw of population
`r`. That form multiplies the chances of the populations together, which
assumes that their draws are independent.

The second way assumes nothing of the kind. For each population it lists every
draw of `num_called_alleles` of its called alleles, weighted by the number of
ways that draw can be taken, runs over every combination of one draw per
population, counts the alleles that are in the population's own draw and in no
other's, and averages over the combinations. Agreement says that the closed
form computes the quantity the spec describes in words, which agreement with
ADZE, the program of Szpiech, Jakobsson and Rosenberg (2008), could not say:
ADZE evaluates the same closed form from the same paper.

Both ways run in `fractions.Fraction`, so they either agree exactly or do not,
and the difference this script writes is an exact rational and not a rounded
one.

Run it from anywhere, with the Python of the project:

    uv run python tests/reference/diversity/enumerate_private.py

It reads nothing. The 18 pairs are the cases of "How it is verified" of "The
private alleles" of docs/specs/diversity.md, written out in `CASES` below: the
three variants of the worked example that have four called alleles in both
populations, and five cases made up for the check. Each case gives the copies
of each allele that every population called at one variant, so a case is one
variant and not a dataset, and the value of a pair is the per variant value
that popnei averages over the variants of a population.

Both draws of a pair are draws of that population's own called alleles, so the
populations of a case hold no individual in common. Two populations that share
one does not appear here: the closed form and the enumeration disagree there,
and "The cases" of the spec says what popnei gives for it and why.

It writes `enumerate_private.tsv` beside itself, a header line and then one
line per pair, 18 of them, in the order of `CASES` and of the populations
within a case. The eight fields of a line, separated by tabs, are:

- `case`, the case as the spec names it.
- `allele_counts`, the copies of each allele that each population called, the
  alleles of one population separated by commas and the populations by ` | `,
  in the order the populations have in `population`.
- `num_called_alleles`, how many called alleles the draw takes, the `g` of the
  formula above.
- `population`, `pop1` for the first population of `allele_counts`, `pop2` for
  the second and `pop3` for the third. The worked example's own two populations
  are called `pop1` and `pop2` in the spec, and they are the first and the
  second here.
- `closed_form` and `enumerated`, the two values as decimals of 17 places, each
  of which reads back as the float64 nearest the exact value. The cargo tests
  of the module assert these.
- `difference`, `closed_form` minus `enumerated` as an exact rational, so it is
  `0` on every line when the two ways agree.
- `exact`, the exact rational the two ways gave, `14/15` and not
  `0.9333333333`, for a reader checking a line by hand.

Nothing is written unless all four of these hold, and the script stops at the
first that does not:

- The two ways give the same exact rational for all 18 pairs.
- Each of the 12 pairs the spec prints a number for is within 5e-11 of it, the
  spec printing 10 decimal places.
- Over the three variants of the worked example the values sum to 1 for `pop1`
  and to 0.9333333333 for `pop2`, with the means 0.3333333333 and
  0.3111111111, which is what the spec gives for it instead of the six per
  variant numbers.
- Every decimal of the file reads back as the float64 nearest the exact
  rational it came from. 17 places are enough for these 18 values and are not
  enough for every float64, so it is checked and not assumed: a decimal that
  lost the last place would be a literal of a cargo test that the module could
  only fail.
"""

from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path
from typing import NamedTuple

HERE = Path(__file__).parent
OUTPUT = HERE / "enumerate_private.tsv"
HEADER = (
    "case\tallele_counts\tnum_called_alleles\tpopulation"
    "\tclosed_form\tenumerated\tdifference\texact"
)

# The spec prints its values to 10 decimal places, so a value of this script
# agrees with one of the spec when it is within half of the last place.
SPEC_TOLERANCE = 5e-11

# The decimal places of the two values of a line. Every value here is between 0
# and 1, and 17 places carry the float64 nearest each of them; `decimal_of` is
# checked against the float it is made from, since 17 places do not round trip
# every float64.
DECIMAL_PLACES = 17


class Case(NamedTuple):
    """One variant: what every population called there, and the draw."""

    name: str
    # One list per population, holding the copies of allele 0, of allele 1 and
    # so on that the population called at this variant.
    allele_counts: list[list[int]]
    num_called_alleles: int
    # What the spec prints for each population of the case, in the order of
    # `allele_counts`, and `None` where the spec gives the sum over the three
    # variants of the worked example instead.
    from_the_spec: list[float] | None


CASES = [
    Case("worked example, variant 1", [[3, 1], [5, 0]], 4, None),
    Case("worked example, variant 3", [[1, 1, 1, 1], [1, 1, 1, 1]], 4, None),
    Case("worked example, variant 5", [[4, 0], [4, 2]], 4, None),
    Case("two populations of two alleles", [[3, 1], [2, 2]], 2, [0.25, 0.4166666667]),
    Case(
        "two populations of three alleles",
        [[2, 1, 1], [3, 0, 1]],
        2,
        [0.75, 0.4166666667],
    ),
    Case(
        "three populations of two alleles",
        [[2, 2], [3, 1], [1, 3]],
        2,
        [0.0, 0.0833333333, 0.0833333333],
    ),
    Case(
        "three populations of four alleles",
        [[2, 1, 1, 0], [1, 1, 0, 2], [2, 0, 1, 1]],
        2,
        [0.5694444444, 0.6805555556, 0.4027777778],
    ),
    Case(
        "a population holding one allele against one holding two",
        [[4, 0], [2, 2]],
        3,
        [0.0, 1.0],
    ),
]

# The three cases the spec gives as one sum over the variants of the worked
# example, with what it gives there for `pop1` and for `pop2`: the sum over the
# three variants, and the mean over them.
WORKED_EXAMPLE = [
    "worked example, variant 1",
    "worked example, variant 3",
    "worked example, variant 5",
]
WORKED_EXAMPLE_SUMS = [1.0, 0.9333333333]
WORKED_EXAMPLE_MEANS = [0.3333333333, 0.3111111111]


def draws_of(allele_counts, num_called_alleles):
    """Every draw of `num_called_alleles` of the called alleles, with its weight.

    A draw is the set of alleles it shows, and its weight is the number of ways
    that set can be taken out of the copies the population called: the alleles
    of a population are a multiset, and two draws of different copies that show
    the same alleles are one draw here. The draws come back in a fixed order,
    sorted by the alleles they show, so that what this script writes cannot
    change from one run to the next.
    """
    copies = [
        allele for allele, count in enumerate(allele_counts) for _ in range(count)
    ]
    weight_of = {}
    for drawn in combinations(range(len(copies)), num_called_alleles):
        alleles = tuple(sorted({copies[copy] for copy in drawn}))
        weight_of[alleles] = weight_of.get(alleles, 0) + 1
    return sorted(weight_of.items())


def by_enumeration(case, pop):
    """The private alleles of `pop`, averaged over every combination of draws.

    It assumes nothing about how the draws of two populations go together: it
    takes every combination of one draw per population, weighted by the ways
    each of those draws can be taken, counts the alleles of `pop`'s draw that
    are in no other draw of the combination, and divides by the total weight.
    """
    per_pop = [
        draws_of(counts, case.num_called_alleles) for counts in case.allele_counts
    ]
    ways = 1
    for counts in case.allele_counts:
        ways *= comb(sum(counts), case.num_called_alleles)
    private = 0
    for combination in product(*per_pop):
        weight = 1
        for _, ways_of_this_draw in combination:
            weight *= ways_of_this_draw
        mine = set(combination[pop][0])
        for other, (alleles, _) in enumerate(combination):
            if other != pop:
                mine -= set(alleles)
        private += weight * len(mine)
    return Fraction(private, ways)


def chance_of_showing(allele_counts, allele, num_called_alleles):
    """The chance that `allele` is in a draw of `num_called_alleles` of these
    called alleles, the `q_r(a)` of the closed form: one minus the chance that
    every copy drawn is one of the other alleles."""
    called = sum(allele_counts)
    copies = allele_counts[allele] if allele < len(allele_counts) else 0
    misses = Fraction(
        comb(called - copies, num_called_alleles), comb(called, num_called_alleles)
    )
    return 1 - misses


def closed_form(case, pop):
    """The formula of the spec: in `pop`'s draw, and in no other population's."""
    expected = Fraction(0)
    for allele, copies in enumerate(case.allele_counts[pop]):
        if copies == 0:
            continue
        term = chance_of_showing(
            case.allele_counts[pop], allele, case.num_called_alleles
        )
        for other, counts in enumerate(case.allele_counts):
            if other == pop:
                continue
            term *= 1 - chance_of_showing(counts, allele, case.num_called_alleles)
        expected += term
    return expected


class Pair(NamedTuple):
    """One line of the file: a case, one of its populations, and the two values."""

    case: Case
    pop: int
    closed: Fraction
    enumerated: Fraction

    @property
    def population(self):
        """`pop1` for the first population of the case, `pop2` for the second."""
        return f"pop{self.pop + 1}"


def every_pair():
    """The 18 pairs, in the order of `CASES` and of the populations of a case."""
    return [
        Pair(case, pop, closed_form(case, pop), by_enumeration(case, pop))
        for case in CASES
        for pop in range(len(case.allele_counts))
    ]


def check_the_two_ways_agree(pairs):
    """Both ways gave the same exact rational for every pair."""
    differ = [pair for pair in pairs if pair.closed != pair.enumerated]
    if differ:
        raise SystemExit(
            f"the closed form and the enumeration differ on {len(differ)} of "
            f"the {len(pairs)} pairs, the first at {differ[0].case.name}, "
            f"{differ[0].population}: the closed form gives {differ[0].closed} "
            f"and every draw gives {differ[0].enumerated}"
        )
    print(f"the two ways agree exactly on all {len(pairs)} pairs")


def check_against_the_spec(pairs):
    """Every value the spec prints for a pair, and the two sums and the two
    means it gives for the worked example in place of its six values."""
    for pair in pairs:
        if pair.case.from_the_spec is None:
            continue
        from_the_spec = pair.case.from_the_spec[pair.pop]
        if abs(float(pair.closed) - from_the_spec) > SPEC_TOLERANCE:
            raise SystemExit(
                f"{pair.case.name}, {pair.population}: this script gives "
                f"{float(pair.closed):.10f} and docs/specs/diversity.md gives "
                f"{from_the_spec:.10f}"
            )
    for pop, (sum_of_the_spec, mean_of_the_spec) in enumerate(
        zip(WORKED_EXAMPLE_SUMS, WORKED_EXAMPLE_MEANS, strict=True)
    ):
        total = sum(
            (
                pair.closed
                for pair in pairs
                if pair.case.name in WORKED_EXAMPLE and pair.pop == pop
            ),
            Fraction(0),
        )
        for got, from_the_spec, what in (
            (total, sum_of_the_spec, "the sum over the three variants"),
            (total / len(WORKED_EXAMPLE), mean_of_the_spec, "the mean over them"),
        ):
            if abs(float(got) - from_the_spec) > SPEC_TOLERANCE:
                raise SystemExit(
                    f"the worked example, pop{pop + 1}: {what} is "
                    f"{float(got):.10f} here and {from_the_spec:.10f} in "
                    f"docs/specs/diversity.md"
                )
    print(
        "every value docs/specs/diversity.md prints for these cases is here, "
        f"within {SPEC_TOLERANCE}"
    )


def decimal_of(value):
    """An exact rational as the decimal a cargo test holds."""
    return f"{float(value):.{DECIMAL_PLACES}f}"


def check_the_decimals_round_trip(pairs):
    """Every decimal of the file reads back as the float64 nearest the exact
    rational it was made from."""
    for pair in pairs:
        for value, column in (
            (pair.closed, "closed_form"),
            (pair.enumerated, "enumerated"),
        ):
            written = decimal_of(value)
            if float(written) != float(value):
                raise SystemExit(
                    f"{pair.case.name}, {pair.population}: {column} is {value} "
                    f"and {DECIMAL_PLACES} decimal places write it as "
                    f"{written}, which reads back as {float(written)!r} and "
                    f"not as {float(value)!r}"
                )
    print(
        f"every decimal of {DECIMAL_PLACES} places reads back as the float64 "
        "nearest its exact value"
    )


def line_of(pair):
    """One line of the file, its eight fields separated by tabs."""
    counts = " | ".join(
        ",".join(str(copies) for copies in counts) for counts in pair.case.allele_counts
    )
    return "\t".join(
        [
            pair.case.name,
            counts,
            str(pair.case.num_called_alleles),
            pair.population,
            decimal_of(pair.closed),
            decimal_of(pair.enumerated),
            str(pair.closed - pair.enumerated),
            str(pair.closed),
        ]
    )


def main():
    """The file is written only when every check passed."""
    pairs = every_pair()
    check_the_two_ways_agree(pairs)
    check_against_the_spec(pairs)
    check_the_decimals_round_trip(pairs)
    OUTPUT.write_text("\n".join([HEADER, *(line_of(pair) for pair in pairs)]) + "\n")
    for pair in pairs:
        print(
            f"  {pair.case.name:56} {pair.population}  "
            f"closed {float(pair.closed):.10f}  "
            f"enumerated {float(pair.enumerated):.10f}  "
            f"difference {pair.closed - pair.enumerated}"
        )
    print(f"{OUTPUT.name}: {len(pairs)} pairs, {OUTPUT.stat().st_size} bytes")


if __name__ == "__main__":
    main()
