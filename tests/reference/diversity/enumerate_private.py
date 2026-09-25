"""It writes the standardized private alleles of 23 pairs of a case and a
population, each computed twice, once by the closed form of
docs/specs/diversity.md and once by enumerating draws.

22 of the 23 are pairs whose populations share no individual, and the two ways
agree on them exactly. The 23rd is the shared individual of "The cases" of the
spec, two populations that are both the one diploid individual `0/1` at a draw
of one allele: the closed form gives 1/2 and the enumeration gives 0. That line
is the one whose `difference` is not 0, and it is there on purpose. The script
refuses to write the file unless the 22 agree exactly and the 23rd disagrees by
exactly 1/2.

The standardized private alleles of a population are, at one variant, the
alleles expected to be in a draw of `num_called_alleles` of that population's
called alleles and in no other population's draw of the same size. "The private
alleles" of the spec computes them with the closed form of Kalinowski (2004),

    E = sum over a of [ q_p(a) * product over r != p of (1 - q_r(a)) ],

where `q_r(a)` is the chance that allele `a` shows in the draw of population
`r`. Each term is the chance that the allele shows in `p`'s draw and in no
other population's draw, and multiplying the chances of the populations
together assumes that their draws are independent.

The second way never writes that form down. For each population it lists every
draw of `num_called_alleles` of its called alleles, weighted by the number of
ways that draw can be taken, runs over every combination of one draw per
population, counts the alleles that are in the population's own draw and in no
other's, and averages over the combinations.

What that check is worth, and what it is not. It catches any error in the
algebra: the decomposition into one term per allele, the chance that an allele
shows in a draw, and the product over the other populations. Dropping the
"in no other population's draw" factor of the closed form makes 19 of these 22
pairs differ, measured on 24 September 2026, and it is what agreement with
ADZE, the program of Szpiech, Jakobsson and Rosenberg (2008), could not catch:
ADZE evaluates the same closed form from the same paper. What
it does not check is the independence: taking the combinations of one draw per
population as a product measure is that same assumption. For the 22 pairs that
is sound rather than circular, because their populations share no individual,
and draws from disjoint sets of gene copies are independent as a fact of the
sampling and not as an assumption.

Where two populations do share an individual the assumption is false and the
enumeration over allele counts is as wrong as the formula. That is what the
23rd pair is for, and it is enumerated differently: over the labelled gene
copies of a set of individuals, where a copy that two populations share is
drawn by both or by neither, so that two populations holding the one individual
draw the same copies. The `enumerated_over` field of every line says which of
the two enumerations gave it.

Both ways run in `fractions.Fraction` throughout, so they agree exactly or do
not, and the difference of a line is an exact rational and not a rounded one.

Run it from anywhere, with the Python of the project:

    uv run python tests/reference/diversity/enumerate_private.py

It reads nothing. The cases are written out in `CASES` and in
`SHARED_INDIVIDUAL` below. Each is one variant and not a dataset: it gives what
every population called at one variant, and the value of a pair is the per
variant value that popnei averages over the variants of a population. Eight of
the ten cases, 18 of the pairs, are the ones "How it is verified" of "The
private alleles" lists, the three variants of the worked example that have four
called alleles in both populations among them. Two more are here because those
18 leave parts of an implementation untested: seven of their values are 0 and
two are 1, and 14 of their 18 population slots have exactly 4 called alleles,
so little of the file would catch an implementation that read one population's
called alleles for another's. The two are a single population, which "The
cases" of the spec says gets every allele it called, and three populations of
3, 5 and 6 called alleles, no two alike and none of them 4.

It writes `enumerate_private.tsv` beside itself, a header line and then one line
per pair, 23 of them, in the order of the cases and of the populations within a
case. The ten fields of a line, separated by tabs, are:

- `case`, the case as the spec names it.
- `enumerated_over`, `allele counts` for the 22 whose populations share no
  individual and `labelled gene copies` for the shared individual.
- `allele_counts`, the copies of each allele that each population called, the
  alleles of one population separated by commas and the populations by ` | `,
  in the order the populations have in `population`. Two populations that are
  the same individual have the same counts here, which is all the closed form
  is given and is why it cannot see the overlap.
- `num_called_alleles`, how many called alleles the draw takes, the `g` of the
  formula above.
- `population`, `pop1` for the first population of `allele_counts`, `pop2` for
  the second and `pop3` for the third. The worked example's own two populations
  are called `pop1` and `pop2` in the spec, and they are the first and the
  second here.
- `closed_form` and `enumerated`, the two values as decimals of 17 places, and
  `closed_form_exact` and `enumerated_exact`, the same two as exact rationals,
  `14/15` and not `0.9333333333`, for a reader checking a line by hand. The
  cargo tests of the module assert the decimals.
- `difference`, `closed_form` minus `enumerated` as an exact rational: `0` on
  every line but the shared individual's, where it is `1/2`.

Nothing is written unless all six of these hold, and the script stops at the
first that does not:

- Every case can be drawn from: at least one population, a draw of at least one
  allele, and no population with fewer called alleles than the draw takes. A
  case that fails this is named, with what was asked of it.
- The two ways give the same exact rational for the 22 pairs whose populations
  share no individual.
- The shared individual gives exactly 1/2 by the closed form and exactly 0 over
  the labelled gene copies, and those gene copies hold the called alleles the
  closed form was given, so that a copy lost on the way to the enumeration
  cannot leave the 0 standing for the wrong reason.
- Each of the 13 pairs the spec prints a number for is within 5e-11 of it, the
  spec printing 10 decimal places, and over the three variants of the worked
  example the values sum to 1 for `pop1` and to 0.9333333333 for `pop2`, with
  the means 0.3333333333 and 0.3111111111, which is what the spec gives for it
  instead of the six per variant numbers.
- The enumeration over labelled gene copies gives what the enumeration over
  allele counts gives for each of the 22, with every gene copy in one
  population only. The two differ where an individual is shared and nowhere
  else.
- Every decimal of the file reads back as the float64 nearest the exact
  rational it came from, so that a reader of the file gets the number the
  arithmetic gives and not one a shorter format rounded. 17 places are enough
  for these 23 values and are not enough for every float64, so it is checked
  and not assumed.
"""

from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path
from typing import NamedTuple

HERE = Path(__file__).parent
OUTPUT = HERE / "enumerate_private.tsv"
HEADER = (
    "case\tenumerated_over\tallele_counts\tnum_called_alleles\tpopulation"
    "\tclosed_form\tclosed_form_exact\tenumerated\tenumerated_exact\tdifference"
)

# What the `enumerated_over` field of a line holds, and the two enumerations it
# names: one over the copies of each allele a population called, which cannot
# see that two populations hold the same individual, and one over the gene
# copies of a set of individuals, each labelled by the individual it belongs to,
# which can.
OVER_ALLELE_COUNTS = "allele counts"
OVER_GENE_COPIES = "labelled gene copies"

# The spec prints its values to 10 decimal places, so a value of this script
# agrees with one of the spec when it is within half of the last place.
SPEC_TOLERANCE = 5e-11

# The decimal places of the two values of a line. Every value here is between 0
# and 3, and 17 places carry the float64 nearest each of them; `decimal_of` is
# checked against the float it is made from, since 17 places do not round trip
# every float64.
DECIMAL_PLACES = 17


class Case(NamedTuple):
    """One variant, given as what every population called there.

    It says nothing about individuals, which is what the closed form is given,
    so the populations of a `Case` are taken to share no gene copy.
    """

    name: str
    # One list per population, holding the copies of allele 0, of allele 1 and
    # so on that the population called at this variant.
    allele_counts: list[list[int]]
    num_called_alleles: int
    # What the spec prints for each population of the case, in the order of
    # `allele_counts`, and `None` where the spec gives no number for it.
    from_the_spec: list[float] | None


class SharedIndividualCase(NamedTuple):
    """One variant, given as individuals and as the individuals of each
    population, so that two populations can hold the same one."""

    name: str
    # The alleles of the gene copies of each individual: [[0, 1]] is one
    # diploid individual whose genotype is 0/1.
    genotypes: list[list[int]]
    # The individuals of each population, by their place in `genotypes`.
    individuals_of: list[list[int]]
    num_called_alleles: int
    # What the spec prints for each population, which for this case is what the
    # closed form gives and not what the draws give.
    from_the_spec: list[float]


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
    # One population holds every allele it called, there being no other
    # population to hold any of them, so its private alleles are its
    # standardized number of alleles. Three alleles of 6 called at a draw of 3,
    # so that the value is neither 0 nor 1. The spec prints no number for it.
    Case("one population of three alleles", [[3, 2, 1]], 3, None),
    # No two populations with the same called alleles and none of them with 4,
    # which 14 of the 18 population slots of the cases above have: a value that
    # took one population's called alleles for another's is wrong here and
    # right in most of them. Allele 3 is in the last population alone.
    Case(
        "three populations of 3, 5 and 6 called alleles",
        [[2, 1], [3, 1, 1], [2, 2, 1, 1]],
        2,
        None,
    ),
]

# Two populations that are both the one diploid individual 0/1, at a draw of one
# allele. Each population called one copy of allele 0 and one of allele 1, which
# is all the closed form sees, and it gives 1/2 for each. The two draws are
# draws of the same genotype, so no allele can go to one population and not the
# other, and the enumeration over the labelled gene copies gives 0. "The cases"
# of docs/specs/diversity.md is where popnei says what it does with the overlap.
SHARED_INDIVIDUAL = SharedIndividualCase(
    name="the shared individual, the one pair the closed form gets wrong",
    genotypes=[[0, 1]],
    individuals_of=[[0], [0]],
    num_called_alleles=1,
    from_the_spec=[0.5, 0.5],
)
# The two populations of that case are the same individual, so their two pairs
# would be the same line twice, and the file holds the first.
SHARED_INDIVIDUAL_POP = 0
SHARED_INDIVIDUAL_CLOSED_FORM = Fraction(1, 2)
SHARED_INDIVIDUAL_ENUMERATED = Fraction(0)

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


def refuse_a_case_that_cannot_be_drawn_from(name, allele_counts, num_called_alleles):
    """A draw of more alleles than a population called, a draw of none, or a
    case of no population, each named rather than left to fail as a division by
    the zero ways of taking such a draw."""
    if not allele_counts:
        raise SystemExit(f"{name}: a case has to have at least one population")
    if num_called_alleles < 1:
        raise SystemExit(
            f"{name}: the draw takes {num_called_alleles} called alleles, and a "
            f"draw takes at least one"
        )
    for pop, counts in enumerate(allele_counts):
        called = sum(counts)
        if called < num_called_alleles:
            raise SystemExit(
                f"{name}, pop{pop + 1}: the draw takes {num_called_alleles} "
                f"called alleles and the population called {called}, "
                f"{','.join(str(copies) for copies in counts)}"
            )


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


def by_enumeration(allele_counts, num_called_alleles, pop):
    """The private alleles of `pop`, averaged over every combination of one draw
    per population, each combination weighted by the ways its draws can be
    taken.

    It writes no formula down, so it checks the algebra of the closed form. It
    does take the combinations as a product measure, which is the closed form's
    own assumption that the draws of two populations are independent, and that
    holds as a fact of the sampling only while the populations share no gene
    copy. `by_enumeration_over_gene_copies` is what does not assume it.
    """
    per_pop = [draws_of(counts, num_called_alleles) for counts in allele_counts]
    ways = 1
    for counts in allele_counts:
        ways *= comb(sum(counts), num_called_alleles)
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


def gene_copies_of(genotypes, individuals):
    """The gene copies of these individuals, each labelled by the individual it
    belongs to and its place in that individual's genotype, so that a copy of an
    individual two populations share is one copy and not two."""
    return [
        (individual, place)
        for individual in individuals
        for place in range(len(genotypes[individual]))
    ]


def by_enumeration_over_gene_copies(genotypes, individuals_of, num_called_alleles, pop):
    """The private alleles of `pop`, averaged over the draws of labelled gene
    copies that the populations can make together.

    A draw of a population is the gene copies it takes, and a copy held by two
    populations is taken by both or by neither: two populations that are the
    same individual make the same draw, and every combination of draws that
    respects that is as likely as any other. Where the populations share no
    copy nothing is constrained and this is `by_enumeration` again, which the
    script checks on every pair of the file but the shared individual's.
    """
    copies_of = [gene_copies_of(genotypes, mine) for mine in individuals_of]
    per_pop = [list(combinations(copies, num_called_alleles)) for copies in copies_of]
    combinations_taken = 0
    private = 0
    for combination in product(*per_pop):
        if not the_draws_agree_on_the_copies_they_share(combination, copies_of):
            continue
        combinations_taken += 1
        mine = {genotypes[individual][place] for individual, place in combination[pop]}
        for other, drawn in enumerate(combination):
            if other != pop:
                mine -= {genotypes[individual][place] for individual, place in drawn}
        private += len(mine)
    return Fraction(private, combinations_taken)


def the_draws_agree_on_the_copies_they_share(combination, copies_of):
    """Whether every gene copy that two populations both hold is drawn by both
    of them or by neither."""
    for first in range(len(combination)):
        for second in range(first + 1, len(combination)):
            shared = set(copies_of[first]) & set(copies_of[second])
            if set(combination[first]) & shared != set(combination[second]) & shared:
                return False
    return True


def as_separate_individuals(allele_counts):
    """The same populations as one haploid individual per called allele, no
    individual in two populations. It is what turns a case given as allele
    counts into one that `by_enumeration_over_gene_copies` reads."""
    genotypes = []
    individuals_of = []
    for counts in allele_counts:
        mine = []
        for allele, copies in enumerate(counts):
            for _ in range(copies):
                mine.append(len(genotypes))
                genotypes.append([allele])
        individuals_of.append(mine)
    return genotypes, individuals_of


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


def closed_form(allele_counts, num_called_alleles, pop):
    """The formula of the spec: in `pop`'s draw, and in no other population's."""
    expected = Fraction(0)
    for allele, copies in enumerate(allele_counts[pop]):
        if copies == 0:
            continue
        term = chance_of_showing(allele_counts[pop], allele, num_called_alleles)
        for other, counts in enumerate(allele_counts):
            if other == pop:
                continue
            term *= 1 - chance_of_showing(counts, allele, num_called_alleles)
        expected += term
    return expected


def allele_counts_of(case):
    """What each population of a case of individuals called: the copies of each
    allele its individuals hold. Two populations that are the same individual
    get the same counts, and that is everything the closed form is given."""
    num_alleles = max(allele for genotype in case.genotypes for allele in genotype) + 1
    counts_of = []
    for individuals in case.individuals_of:
        counts = [0] * num_alleles
        for individual in individuals:
            for allele in case.genotypes[individual]:
                counts[allele] += 1
        counts_of.append(counts)
    return counts_of


class Pair(NamedTuple):
    """One line of the file: a case, one of its populations, and the two values."""

    case: str
    enumerated_over: str
    allele_counts: list[list[int]]
    num_called_alleles: int
    pop: int
    closed: Fraction
    enumerated: Fraction
    # What the spec prints for this pair, and `None` where it prints nothing.
    from_the_spec: float | None

    @property
    def population(self):
        """`pop1` for the first population of the case, `pop2` for the second."""
        return f"pop{self.pop + 1}"

    @property
    def difference(self):
        """The closed form minus the enumeration, an exact rational."""
        return self.closed - self.enumerated


def every_pair():
    """The 23 pairs: the populations of every case of `CASES`, in order, and
    then the one pair of the shared individual."""
    pairs = []
    for case in CASES:
        refuse_a_case_that_cannot_be_drawn_from(
            case.name, case.allele_counts, case.num_called_alleles
        )
        for pop in range(len(case.allele_counts)):
            pairs.append(
                Pair(
                    case.name,
                    OVER_ALLELE_COUNTS,
                    case.allele_counts,
                    case.num_called_alleles,
                    pop,
                    closed_form(case.allele_counts, case.num_called_alleles, pop),
                    by_enumeration(case.allele_counts, case.num_called_alleles, pop),
                    None if case.from_the_spec is None else case.from_the_spec[pop],
                )
            )
    shared = SHARED_INDIVIDUAL
    counts = allele_counts_of(shared)
    refuse_a_case_that_cannot_be_drawn_from(
        shared.name, counts, shared.num_called_alleles
    )
    pairs.append(
        Pair(
            shared.name,
            OVER_GENE_COPIES,
            counts,
            shared.num_called_alleles,
            SHARED_INDIVIDUAL_POP,
            closed_form(counts, shared.num_called_alleles, SHARED_INDIVIDUAL_POP),
            by_enumeration_over_gene_copies(
                shared.genotypes,
                shared.individuals_of,
                shared.num_called_alleles,
                SHARED_INDIVIDUAL_POP,
            ),
            shared.from_the_spec[SHARED_INDIVIDUAL_POP],
        )
    )
    return pairs


def pairs_that_share_no_individual(pairs):
    """The pairs whose populations hold no gene copy in common, which is every
    one of them but the shared individual's."""
    return [pair for pair in pairs if pair.enumerated_over == OVER_ALLELE_COUNTS]


def check_the_two_ways_agree(pairs):
    """Both ways gave the same exact rational for every pair whose populations
    share no individual."""
    agreeing = pairs_that_share_no_individual(pairs)
    differ = [pair for pair in agreeing if pair.difference != 0]
    if differ:
        raise SystemExit(
            f"the closed form and the enumeration differ on {len(differ)} of "
            f"the {len(agreeing)} pairs whose populations share no individual, "
            f"the first at {differ[0].case}, {differ[0].population}: the closed "
            f"form gives {differ[0].closed} and every draw gives "
            f"{differ[0].enumerated}"
        )
    print(f"the two ways agree exactly on all {len(agreeing)} pairs")


def check_the_shared_individual(pairs):
    """The one pair that disagrees does so by exactly the 1/2 of "The cases" of
    the spec, there is exactly one of it, and the gene copies it was enumerated
    over are the called alleles the closed form was given."""
    counts_of = allele_counts_of(SHARED_INDIVIDUAL)
    for pop, individuals in enumerate(SHARED_INDIVIDUAL.individuals_of):
        counts = [0] * len(counts_of[pop])
        for individual, place in gene_copies_of(
            SHARED_INDIVIDUAL.genotypes, individuals
        ):
            counts[SHARED_INDIVIDUAL.genotypes[individual][place]] += 1
        if counts != counts_of[pop]:
            raise SystemExit(
                f"{SHARED_INDIVIDUAL.name}, pop{pop + 1}: its labelled gene "
                f"copies hold {counts} and the called alleles the closed form "
                f"is given are {counts_of[pop]}"
            )
    shared = [pair for pair in pairs if pair.enumerated_over == OVER_GENE_COPIES]
    if len(shared) != 1:
        raise SystemExit(
            f"the file holds {len(shared)} pairs enumerated over "
            f"{OVER_GENE_COPIES} and it holds one, the shared individual"
        )
    pair = shared[0]
    if (
        pair.closed != SHARED_INDIVIDUAL_CLOSED_FORM
        or pair.enumerated != SHARED_INDIVIDUAL_ENUMERATED
    ):
        raise SystemExit(
            f"{pair.case}, {pair.population}: the closed form gives "
            f"{pair.closed} where docs/specs/diversity.md has "
            f"{SHARED_INDIVIDUAL_CLOSED_FORM}, and the labelled gene copies "
            f"give {pair.enumerated} where it has "
            f"{SHARED_INDIVIDUAL_ENUMERATED}"
        )
    print(
        f"the shared individual disagrees by exactly {pair.difference}: the "
        f"closed form gives {pair.closed} and the labelled gene copies give "
        f"{pair.enumerated}"
    )


def check_against_the_spec(pairs):
    """Every value the spec prints for a pair, and the two sums and the two
    means it gives for the worked example in place of its six values."""
    printed = 0
    for pair in pairs:
        if pair.from_the_spec is None:
            continue
        printed += 1
        if abs(float(pair.closed) - pair.from_the_spec) > SPEC_TOLERANCE:
            raise SystemExit(
                f"{pair.case}, {pair.population}: this script gives "
                f"{float(pair.closed):.10f} and docs/specs/diversity.md gives "
                f"{pair.from_the_spec:.10f}"
            )
    for pop, (sum_of_the_spec, mean_of_the_spec) in enumerate(
        zip(WORKED_EXAMPLE_SUMS, WORKED_EXAMPLE_MEANS, strict=True)
    ):
        total = sum(
            (
                pair.closed
                for pair in pairs
                if pair.case in WORKED_EXAMPLE and pair.pop == pop
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
        f"the {printed} values docs/specs/diversity.md prints for these pairs "
        f"are here, within {SPEC_TOLERANCE}, and so are the two sums and the "
        f"two means it gives for the worked example"
    )


def check_the_two_enumerations_agree(pairs):
    """Over gene copies that no two populations share, the enumeration that
    labels them gives what the enumeration over allele counts gives. It is what
    says that the shared individual's 0 comes from the sharing and not from a
    second way of counting."""
    for pair in pairs_that_share_no_individual(pairs):
        genotypes, individuals_of = as_separate_individuals(pair.allele_counts)
        over_copies = by_enumeration_over_gene_copies(
            genotypes, individuals_of, pair.num_called_alleles, pair.pop
        )
        if over_copies != pair.enumerated:
            raise SystemExit(
                f"{pair.case}, {pair.population}: the enumeration over allele "
                f"counts gives {pair.enumerated} and the one over labelled "
                f"gene copies, with no copy in two populations, gives "
                f"{over_copies}"
            )
    print(
        "the two enumerations agree on every pair whose populations share no individual"
    )


def decimal_of(value):
    """An exact rational as the decimal a cargo test holds."""
    return f"{float(value):.{DECIMAL_PLACES}f}"


def check_the_decimals_round_trip(pairs):
    """Every decimal of the file reads back as the float64 nearest the exact
    rational it was made from, so that what the file shows is what the
    arithmetic gives and not what a shorter format rounded it to."""
    for pair in pairs:
        for value, field in (
            (pair.closed, "closed_form"),
            (pair.enumerated, "enumerated"),
        ):
            written = decimal_of(value)
            if float(written) != float(value):
                raise SystemExit(
                    f"{pair.case}, {pair.population}: {field} is {value} and "
                    f"{DECIMAL_PLACES} decimal places write it as {written}, "
                    f"which reads back as {float(written)!r} and not as "
                    f"{float(value)!r}"
                )
    print(
        f"every decimal of {DECIMAL_PLACES} places reads back as the float64 "
        f"nearest its exact value"
    )


def line_of(pair):
    """One line of the file, its ten fields separated by tabs."""
    counts = " | ".join(
        ",".join(str(copies) for copies in counts) for counts in pair.allele_counts
    )
    return "\t".join(
        [
            pair.case,
            pair.enumerated_over,
            counts,
            str(pair.num_called_alleles),
            pair.population,
            decimal_of(pair.closed),
            str(pair.closed),
            decimal_of(pair.enumerated),
            str(pair.enumerated),
            str(pair.difference),
        ]
    )


def main():
    """The file is written only when every check passed."""
    pairs = every_pair()
    check_the_two_ways_agree(pairs)
    check_the_shared_individual(pairs)
    check_against_the_spec(pairs)
    check_the_two_enumerations_agree(pairs)
    check_the_decimals_round_trip(pairs)
    OUTPUT.write_text("\n".join([HEADER, *(line_of(pair) for pair in pairs)]) + "\n")
    for pair in pairs:
        print(
            f"  {pair.case:58} {pair.population}  "
            f"closed {float(pair.closed):.10f}  "
            f"enumerated {float(pair.enumerated):.10f}  "
            f"difference {pair.difference}"
        )
    agreeing = len(pairs_that_share_no_individual(pairs))
    print(
        f"{OUTPUT.name}: {len(pairs)} pairs, {agreeing} of them with a "
        f"difference of 0 and the shared individual with "
        f"{SHARED_INDIVIDUAL_CLOSED_FORM}, {OUTPUT.stat().st_size} bytes"
    )


if __name__ == "__main__":
    main()
