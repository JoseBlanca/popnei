"""How much variety each population of a dataset holds.

:func:`popnei.calc_pop_diversity` makes one pass over the variants and
gives, for each population, how many alleles it called, how many of those no
other population called, how many of the variants vary in it, how its
variants are spread over the count of their rarer allele, and how far its
genotypes are from the proportions its allele frequencies would give. The
first three come both as they stand and standardized to a common number of
called alleles, so that a population of 20 individuals and one of 200 can be
compared.

Standardized means drawn: take `num_called_alleles` of the called alleles a
population has at a variant, without replacement, and take the expectation
over every such draw. Applied to a count of alleles that operation is called
rarefaction and applied to the spectrum it is called projection, and both
come from that one argument.

A population is a named set of individuals, and `pops` names them as
:func:`popnei.calc_per_var_distribs` does: a dict of population name to the
names of its individuals, with ``None`` meaning one population, ``pop``, of
every individual.

`docs/specs/diversity.md` has what each of the five statistics is, the
program each of its numbers was checked against and the numbers the tests
assert.
"""

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from enum import StrEnum

import numpy
import pandas

from popnei import _core
from popnei.stats import _the_pops
from popnei.variant import PassStats, Variants, _pass_stats_of


class PopDiversityStat(StrEnum):
    """The five statistics :func:`popnei.calc_pop_diversity` calculates.

    The value of each member is the name of the field of
    :class:`PopDiversity` that holds its result.
    """

    NUM_ALLELES = "num_alleles"
    """How many alleles a population called, as a total over its variants, as
    a mean over them, which is the allelic richness, and in a draw."""

    PRIVATE_ALLELES = "private_alleles"
    """Of those, the ones no other population of the call called at the same
    variant."""

    VARIABLE_VARS_RATIO = "variable_vars_ratio"
    """How many of the variants the population called more than one allele
    at."""

    FOLDED_SFS = "folded_sfs"
    """How many variants of the population show each count of their rarer
    allele in a draw, which needs `num_called_alleles`."""

    FIS = "fis"
    """How far the genotypes of the population are from the proportions its
    allele frequencies would give if its individuals paired at random."""


@dataclass(frozen=True)
class PopDiversity:
    """What :func:`popnei.calc_pop_diversity` gives back.

    A statistic that was not asked for is ``None``. Every frame is indexed by
    population name, in the order of the keys of `pops`, which is the order
    of :attr:`pops`.

    A population for which no variant counted has 0 in every count, NaN in
    every mean and ratio and NaN in :attr:`fis`, and a column of zeros in the
    spectrum. It is not an error: a user whose filters left one population
    with nothing still wants the others.
    """

    pops: tuple[str, ...]
    """The names of the populations, in the order of the keys of `pops`."""

    num_alleles: pandas.DataFrame | None
    """How many alleles each population called, in three columns, and
    ``None`` when it was not asked for.

    `total` is the sum over the variants that counted for the population. An
    allele numbered 1 at one variant is not the allele numbered 1 at the
    next, so it is a sum of per variant counts and never a count of distinct
    things across the dataset. `mean` is that total over
    ``num_vars.with_data``, the allelic richness a user compares between
    populations, and NaN where no variant counted. `in_draw` is the alleles a
    draw of `num_called_alleles` is expected to show, averaged over the
    variants in the draw for the population, and NaN when there was no draw
    or no variant reached it.
    """

    private_alleles: pandas.DataFrame | None
    """How many alleles each population called that no other population of
    the call called at the same variant, in the same three columns, and
    ``None`` when it was not asked for.

    The divisor of `mean` is :attr:`num_vars_every_pop` and not
    ``num_vars.with_data``: a variant where one population has too little
    called is out of the private alleles of every population, since there
    every allele of every other population would be private and the count
    would measure the missing data. `in_draw` is over
    :attr:`num_vars_every_pop_in_draw` for the same reason.

    With one population every allele it called is private, since there is no
    other population to hold it, so `total` is then the `total` of
    :attr:`num_alleles`.
    """

    variable_vars_ratio: pandas.DataFrame | None
    """How many of the variants each population called more than one allele
    at, and ``None`` when it was not asked for.

    `total` is those variants, which are the `num_variable` of the
    `poly_vars_ratio` of :func:`popnei.calc_per_var_distribs` over the same
    variants; `ratio` is that total over ``num_vars.with_data``, NaN where no
    variant counted; and `in_draw` is the chance that a draw of
    `num_called_alleles` is not all of one allele, averaged over the variants
    in the draw.
    """

    folded_sfs: pandas.DataFrame | None
    """How many variants of each population show each count of their rarer
    allele in a draw of `num_called_alleles`, and ``None`` when it was not
    asked for.

    One row per count of the rarer allele, 0 to ``num_called_alleles // 2``,
    indexed by that count, and one column per population. It is folded
    because nothing in a VCF says which allele is the ancestral one, so the
    counts `j` and ``num_called_alleles - j`` are one bin.
    """

    fis: pandas.Series | None
    """One minus the mean observed heterozygosity of each population over its
    mean unbiased expected one, and ``None`` when it was not asked for.

    It is 0 when the genotypes of the population are in the proportions its
    allele frequencies would give if its individuals paired at random,
    positive when it holds fewer heterozygous genotypes than that, which
    inbreeding, selfing and a population split into unmixed groups all
    produce, and negative when it holds more. It is Nei's F_IS, read on one
    population on its own, and not Weir and Cockerham's, which comes out of a
    decomposition of the variance across populations.

    It is NaN for a population that has no F_IS, in four cases. One no variant
    counted for. One no counted variant of which carries both
    heterozygosities, which a population whose counted variants hold no whole
    called genotype reaches: two individuals whose genotypes are all half
    called count their variants and have no observed heterozygosity at any of
    them. One whose mean unbiased expected heterozygosity is 0, every variant
    it counted having held one allele. And every population of a haploid
    dataset, where no genotype can be heterozygous and the ratio would be 1
    wherever the population has any diversity. The draw does not touch it:
    the observed heterozygosity is a property of whole genotypes and not of a
    sample of alleles.
    """

    num_vars: pandas.DataFrame
    """How many variants counted for each population, in two columns.

    `with_data` is the variants the population called something at and had
    `min_num_individuals` called genotypes in, which is the divisor of `mean`
    and of `ratio`. `in_draw` is how many of those also reached
    `num_called_alleles`, which is the divisor of every `in_draw` value.
    """

    num_vars_every_pop: int
    """The variants that counted for every population, which is the divisor
    of the `mean` of :attr:`private_alleles`."""

    num_vars_every_pop_in_draw: int
    """Of those, the ones every population reached `num_called_alleles` at,
    which is the divisor of the `in_draw` of :attr:`private_alleles`."""

    pass_stats: PassStats
    """How many variants the pass gave, after the steps of the ``Variants``,
    and what each filter of it was given and kept."""


def calc_pop_diversity(
    variants: Variants,
    pops: dict[str, Sequence[str]] | None = None,
    stats: Iterable[PopDiversityStat] = tuple(PopDiversityStat),
    num_called_alleles: int | None = None,
    min_num_individuals: int = _core.DEFAULT_MIN_NUM_INDIVIDUALS,
) -> PopDiversity:
    """How much variety each population holds, in one pass over `variants`.

    The five statistics are how many alleles a population called, how many of
    those no other population called, how many of the variants vary in it,
    how its variants are spread over the count of their rarer allele, and
    F_IS, how far its genotypes are from the proportions its allele
    frequencies would give. All five are built from the same thing, how often
    each population called each allele at each variant, so asking for several
    of them costs one reading of the variants and not five.

    It is a consumer of the `variants`: it makes one pass over the source
    through the steps the ``Variants`` has when it is called, and the
    ``Variants`` is as it was afterwards.

    `pops` is a dict of population name to the names of its individuals,
    which are looked up among :attr:`popnei.Variants.individuals`, the ones
    the pass gives. With no `pops` there is one population, ``pop``, of every
    individual, and then every allele it called is private. A name that is
    not an individual of the pass, a name twice in one population, a
    population that names no individual and a `pops` with no population are
    each a ``ValueError``; an individual in two populations is taken, and
    counted in each of them.

    `stats` says which of the five to calculate, as members of
    :class:`popnei.PopDiversityStat`, all five by default. One that is not
    asked for is ``None`` in the result. Anything that is not a member, a
    name written as a string among them, is a ``TypeError``, so that a name
    with a typo in it cannot pass; no statistic at all is a ``ValueError``.

    `num_called_alleles` is how many called alleles every population is
    brought down to, so that a population of 20 individuals and one of 200
    can be compared: draw that many of the alleles a population called at a
    variant, without replacement, and take the expectation over every such
    draw. With no `num_called_alleles`, the default, the three `in_draw`
    columns are NaN and the spectrum cannot be asked for, since the bins of a
    spectrum need one number of alleles for every population and every
    variant: asking for it without one is a ``ValueError``, and so is a
    `num_called_alleles` below 2, a draw of one allele showing one allele
    whatever the population holds. A `num_called_alleles` above every
    population's called alleles is no error: every `in_draw` value is NaN,
    every count of the spectrum is 0 and ``num_vars.in_draw`` says why.

    It takes one number and not a sequence of them, so a user who wants the
    curve of allelic richness against the number of alleles drawn, which
    shows whether a population has been sampled enough, calls this once per
    point of it.

    `min_num_individuals` is how many called genotypes a population needs at
    a variant for the variant to count for it, 20 by default, measured as the
    called alleles of the population over the ploidy, so a half called
    genotype counts as half an individual. Strictly fewer and the variant
    does not count for that population. A population that called nothing at a
    variant does not count it whatever `min_num_individuals` is, so a 0 does
    not put variants with no data into the totals.

    A variant is in the draw for a population when it counts for that
    population by that rule and the population called at least
    `num_called_alleles` there, both and not the second alone.

    A pass that gives no variant is a ``ValueError``, whether the source
    holds none or the steps kept none.

    The draw is not built yet, so a `num_called_alleles` and a `stats` that
    names the spectrum are each a ``NotImplementedError`` that says so: what
    is given today is the four statistics over the called alleles each
    population has, and work package 3 of `docs/plans/diversity.md` adds the
    three `in_draw` columns and `folded_sfs`. A `num_called_alleles` that is no
    whole number or is out of range is refused as the wrong argument it is,
    whatever popnei computes.

    pyNei has none of the five, so no result here is compared with it.
    """
    if not isinstance(variants, Variants):
        # The path of the VCF whose variants are read is the mistake that is
        # easiest to make, and what it gave was the `AttributeError` of an
        # object with no source inside it.
        raise TypeError(
            f"`variants` is {variants!r}, a {type(variants).__name__}, and "
            f"`calc_pop_diversity` reads the variants of a source: give it what "
            f"`open_vcf` or `open_vars` gives, "
            f"calc_pop_diversity(open_vcf(vcf_path))"
        )
    asked_for = _the_stats(stats)
    _refuse_what_is_not_built_yet(asked_for, num_called_alleles)
    named = _the_pops(pops)
    (
        pop_names,
        num_vars_with_data,
        num_vars_in_draw,
        (num_vars_every_pop, num_vars_every_pop_in_draw),
        num_alleles,
        private_alleles,
        variable_vars,
        folded_sfs,
        fis,
        counts,
    ) = _core.calc_pop_diversity(
        variants._source,
        variants._steps,
        named,
        asked_for,
        num_called_alleles,
        min_num_individuals,
    )
    names = list(pop_names)
    # The divisor of the private alleles is the variants every population
    # counted and not the ones each of them counted, as
    # `PopDiversity.private_alleles` says.
    every_pop = numpy.full(len(names), num_vars_every_pop, dtype=numpy.int64)
    return PopDiversity(
        pops=tuple(names),
        num_alleles=_counts_of(names, num_alleles, num_vars_with_data, "mean"),
        private_alleles=_counts_of(names, private_alleles, every_pop, "mean"),
        variable_vars_ratio=_counts_of(
            names, variable_vars, num_vars_with_data, "ratio"
        ),
        folded_sfs=_spectrum_of(names, folded_sfs),
        fis=None if fis is None else pandas.Series(fis, index=names),
        num_vars=pandas.DataFrame(
            {"with_data": num_vars_with_data, "in_draw": num_vars_in_draw},
            index=names,
        ),
        num_vars_every_pop=int(num_vars_every_pop),
        num_vars_every_pop_in_draw=int(num_vars_every_pop_in_draw),
        pass_stats=_pass_stats_of(counts),
    )


# The four statistics the pass computes, as a Python user writes them, for the
# message that refuses the fifth. They are read from the enumeration, so a
# member renamed there is renamed in the message too.
_THE_FOUR_THAT_ARE_COMPUTED = ", ".join(
    f"PopDiversityStat.{stat.name}"
    for stat in PopDiversityStat
    if stat is not PopDiversityStat.FOLDED_SFS
)


def _refuse_what_is_not_built_yet(
    asked_for: list[str], num_called_alleles: int | None
) -> None:
    """It refuses the draw and the folded spectrum while the pass computes
    neither.

    The pass gives the four statistics that are over the called alleles each
    population has, and work package 3 of `docs/plans/diversity.md` adds the
    draw to it. Until then a call that asked for either would get NaN in the
    three `in_draw` columns and ``None`` in `folded_sfs`, which are what this
    module gives a draw above every population's called alleles and a
    statistic nobody asked for: the answer would be wrong and would read as an
    answer. The two refusals go when the values arrive, with the paragraph of
    :func:`calc_pop_diversity` that names them and the marker of the tests
    that assert them.

    Raises:
        NotImplementedError: when `stats` names the spectrum, and when
            `num_called_alleles` is a draw the pass would compute over.
    """
    if PopDiversityStat.FOLDED_SFS in asked_for:
        raise NotImplementedError(
            f"popnei does not compute the folded site frequency spectrum yet: it "
            f"is work package 3 of `docs/plans/diversity.md`, and `stats` names "
            f"it, a call that gives no `stats` asking for all five statistics. "
            f"Write stats=({_THE_FOUR_THAT_ARE_COMPUTED},) for the four that "
            f"are computed."
        )
    if _is_a_draw_the_pass_would_take(num_called_alleles):
        raise NotImplementedError(
            f"popnei does not take the alleles of a population down to a common "
            f"number yet: the draw is work package 3 of "
            f"`docs/plans/diversity.md`, and `num_called_alleles` is "
            f"{num_called_alleles!r}. Leave `num_called_alleles` out for the "
            f"totals, the means, the ratios and F_IS, which are over the called "
            f"alleles each population has."
        )


def _is_a_draw_the_pass_would_take(num_called_alleles: int | None) -> bool:
    """Whether `num_called_alleles` is a draw the pass would compute over once
    the draw is built.

    Anything else is a wrong argument whatever popnei computes, and it is
    refused as one where every argument of the call is refused, by the binding
    crate for what is no whole number and for what is negative or larger than a
    count of called alleles, and by the Rust core for a draw of fewer than two
    alleles: a user reads the same sentence for it now as they will then. What
    the bound here does not have is that upper limit, so a whole number above
    it is told that the draw is not built where it will be told that it is out
    of range; writing the limit again here would give a user two of them for
    one argument.
    """
    if isinstance(num_called_alleles, bool) or not isinstance(num_called_alleles, int):
        return False
    return num_called_alleles >= 2


def _the_stats(stats: Iterable[PopDiversityStat]) -> list[str]:
    """The names of the statistics a user asked for, each once and in the
    order they named them.

    A member of a ``StrEnum`` is a string, so one written on its own is a
    sequence of its letters: it is taken as that one statistic, as a name
    written as a string is refused.
    """
    if isinstance(stats, PopDiversityStat):
        stats = (stats,)
    elif isinstance(stats, str):
        # A string is a sequence of its letters, so a name written where the
        # members go would be refused for its first letter, `f`, and the user
        # would read about a statistic they never wrote.
        raise TypeError(
            f"`stats` takes the members of `PopDiversityStat`, and {stats!r}, a "
            f"{type(stats).__name__}, is not one of them: write "
            f"stats=(PopDiversityStat.FIS,) for F_IS"
        )
    try:
        stats = list(stats)
    except TypeError:
        # What Python says of its own here, `'int' object is not iterable`,
        # names neither the argument nor the call.
        raise TypeError(
            f"`stats` is a sequence of the members of `PopDiversityStat`, and "
            f"{stats!r}, a {type(stats).__name__}, was given: write "
            f"stats=(PopDiversityStat.FIS,) for F_IS"
        ) from None
    asked_for: list[str] = []
    for stat in stats:
        if not isinstance(stat, PopDiversityStat):
            raise TypeError(
                f"`stats` takes the members of `PopDiversityStat`, and {stat!r}, a "
                f"{type(stat).__name__}, is not one of them: write "
                f"stats=(PopDiversityStat.FIS,) for F_IS"
            )
        if str(stat) not in asked_for:
            asked_for.append(str(stat))
    # A `stats` that names no statistic is refused by the Rust core, which says
    # that such a pass reads every variant of the source and computes nothing
    # of them, and lists the five names a user can write. That sentence was
    # written here and in the TypeScript package before it was written there,
    # and it is not written here any more, so a third language does not write
    # it a third time.
    return asked_for


def _counts_of(
    names: Sequence[str],
    counts: tuple[numpy.ndarray, numpy.ndarray] | None,
    divisor: numpy.ndarray,
    over_the_variants: str,
) -> pandas.DataFrame | None:
    """One count of alleles or of variants as its three columns, or ``None``
    when nobody asked for the statistic.

    The core gives the total and the standardized value, and the middle
    column is the total over the variants that the statistic is a mean or a
    ratio over, which is the one arithmetic of this layer.
    """
    if counts is None:
        return None
    total, in_draw = counts
    return pandas.DataFrame(
        {
            "total": total,
            over_the_variants: _over(total, divisor),
            "in_draw": in_draw,
        },
        index=list(names),
    )


def _over(total: numpy.ndarray, divisor: numpy.ndarray) -> numpy.ndarray:
    """`total` over `divisor`, one value per population, and NaN where the
    divisor is 0: a population no variant counted for has no mean and no
    ratio, and that is not a division by zero to warn about."""
    total = numpy.asarray(total, dtype=numpy.float64)
    divisor = numpy.asarray(divisor, dtype=numpy.float64)
    return numpy.divide(
        total,
        divisor,
        out=numpy.full(total.shape, numpy.nan),
        where=divisor > 0,
    )


def _spectrum_of(
    names: Sequence[str], folded_sfs: numpy.ndarray | None
) -> pandas.DataFrame | None:
    """The folded site frequency spectrum, one row per count of the rarer
    allele and one column per population, or ``None`` when nobody asked for
    it."""
    if folded_sfs is None:
        return None
    return pandas.DataFrame(
        folded_sfs,
        index=pandas.RangeIndex(len(folded_sfs), name="rarer_allele"),
        columns=list(names),
    )
