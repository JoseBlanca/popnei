"""The statistics of the variants and of the individuals, per population.

A population is a named set of individuals that a calculation treats as a
group, and `pops` is how a user names them: a dict of population name to the
names of its individuals. Every statistic here is calculated for each
population over its individuals alone, and every result is keyed by
population name, in the order the keys of `pops` iterate in. With no `pops`
there is one population, named ``pop``, of every individual.

:func:`popnei.calc_per_var_distribs` makes one pass over the variants and
gives, for each statistic and population, the mean over the variants that
had a value and a histogram of them. The per variant values are not kept: a
million of them for each population do not fit a browser tab, and a user who
wants them takes the genotypes with :meth:`popnei.Variants.iter_blocks`.
"""

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum

import numpy
import pandas

from popnei import _core
from popnei.variant import PassStats, Variants, _pass_stats_of


class PerVarStat(StrEnum):
    """The five statistics :func:`popnei.calc_per_var_distribs` calculates.

    The value of each member is the name of the field of
    :class:`PerVarDistribs` that holds its result.
    """

    OBS_HET = "obs_het"
    """The heterozygous genotypes of a population over its called ones."""

    MAF = "maf"
    """The commonest allele of a population over its called alleles."""

    EXP_HET = "exp_het"
    """The chance that gene copies taken at random from a population are not
    all of the same allele, from the frequencies as they are."""

    UNBIASED_EXP_HET = "unbiased_exp_het"
    """The same corrected for the frequencies being estimated from the copies
    the statistic is computed over."""

    POLY_VARS_RATIO = "poly_vars_ratio"
    """How many of the variants vary in a population, in three counts and two
    ratios, which is a count and not a distribution."""


@dataclass(frozen=True)
class StatsDistrib:
    """The distribution of one statistic over the variants of a pass, per
    population.

    A variant has no value of a statistic in a population when the population
    has too little data at it, and such a variant is out of the mean and in
    no bin, so the histograms of two populations can count different numbers
    of variants. A value outside the range of the bins is in the mean and in
    no bin, as :func:`numpy.histogram` leaves it out too.
    """

    mean: pandas.Series
    """The mean over the variants that had a value, one value per
    population, NaN for a population in which no variant had one."""

    hist_bin_edges: numpy.ndarray
    """The edges of the bins, one more than there are bins.

    The four distributions of one result share this array, as pyNei's do, so
    it is read only: a number written into the edges of one statistic would
    be in the edges of the other three."""

    hist_counts: pandas.DataFrame
    """How many variants fell in each bin, one row per bin and one column
    per population."""


@dataclass(frozen=True)
class PolyVarsStats:
    """How many of the variants of a pass vary in each population.

    A variant is polymorphic in a population when its major allele frequency
    there is below `poly_threshold`, strictly, and variable when that
    frequency is below 1. Both are counted among the variants that have a
    major allele frequency in the population.
    """

    num_poly: pandas.Series
    """The polymorphic variants of each population."""

    poly_ratio: pandas.Series
    """`num_poly` over `tot_num_variants_with_data`, NaN when the latter is
    0."""

    poly_ratio_over_variables: pandas.Series
    """`num_poly` over `num_variable`, NaN when the latter is 0."""

    num_variable: pandas.Series
    """The variable variants of each population."""

    tot_num_variants_with_data: pandas.Series
    """The variants that have a major allele frequency in each population,
    which the other two counts are among."""


@dataclass(frozen=True)
class PerVarDistribs:
    """What :func:`popnei.calc_per_var_distribs` gives back.

    A statistic that was not asked for is ``None``.
    """

    obs_het: StatsDistrib | None
    """The distribution of the observed heterozygosity."""

    maf: StatsDistrib | None
    """The distribution of the major allele frequency."""

    exp_het: StatsDistrib | None
    """The distribution of the plain expected heterozygosity."""

    unbiased_exp_het: StatsDistrib | None
    """The distribution of the unbiased expected heterozygosity."""

    poly_vars_ratio: PolyVarsStats | None
    """The counts of the polymorphism ratio."""

    pass_stats: PassStats
    """How many variants the pass gave, after the steps of the ``Variants``,
    and what each filter of it was given and kept."""


def calc_per_var_distribs(
    variants: Variants,
    stats: Iterable[PerVarStat] = tuple(PerVarStat),
    pops: dict[str, Sequence[str]] | None = None,
    min_num_individuals: int = _core.DEFAULT_MIN_NUM_INDIVIDUALS,
    hist_kwargs: dict | None = None,
    ploidy: int | None = None,
    poly_threshold: float = _core.DEFAULT_POLY_THRESHOLD,
) -> PerVarDistribs:
    """Up to five statistics of every variant and every population, in one
    pass over `variants`, as a mean and a histogram each.

    The statistics are the observed heterozygosity, the heterozygous
    genotypes of a population over its called ones; the major allele
    frequency, the count of the commonest allele over the called alleles;
    the expected heterozygosity, the chance that gene copies taken at random
    from the population are not all of the same allele, plain and corrected
    for the frequencies being estimated from the copies the statistic is
    computed over; and the polymorphism ratio, how many of the variants vary
    in the population, which is a count and not a distribution.

    It is a consumer of the `variants`: it makes one pass over the source
    through the steps the ``Variants`` has when it is called, and the
    ``Variants`` is as it was afterwards.

    `stats` says which of the five to calculate, as members of
    :class:`popnei.PerVarStat`, all five by default. Anything that is not a
    member, a name written as a string among them, is a ``TypeError``, so
    that a name with a typo in it cannot pass; no statistic at all is a
    ``ValueError``. Asking for fewer is a saving of work and changes no
    value.

    `pops` is a dict of population name to the names of its individuals,
    which are looked up among :attr:`popnei.Variants.individuals`, the ones
    the pass gives. With no `pops` there is one population, ``pop``, of
    every individual. A name that is not an individual of the pass, a name
    twice in one population, a population that names no individual and a
    `pops` with no population are each a ``ValueError``; an individual in
    two populations is taken, and counted in each of them.

    `min_num_individuals` is how many called genotypes a population needs at
    a variant for the variant to have a value there, 20 by default. The test
    is on the called data counted in genotypes, the called alleles of the
    population over the ploidy, which is a half when a genotype is half
    called, and the variant has no value when that number is strictly less
    than the threshold. A variant with no value in a population is out of
    the mean and in no bin of the histogram of that population.

    `hist_kwargs` is the histogram, under three keys: ``range``, the two
    ends, ``(0, 1)`` by default, which is where the statistics live;
    ``num_bins``, 40 by default; and ``bin_type``, ``"linear"`` for bins of
    equal width or ``"logarithmic"`` for bins of equal ratio, whose
    ``range`` has to start above 0. The dict is read and not changed, and a
    key that is none of the three is a ``ValueError``, since it would leave
    the default in its place with nothing said. A value falls in the bin
    whose left edge is at most the value and whose right edge is above it,
    and the last bin takes its right edge too, as :func:`numpy.histogram`
    does; a value outside the range is in no bin and in the mean.

    `ploidy` is the two expected heterozygosities' alone: the number the
    allele frequencies are raised to, and how many copies the unbiased one
    draws, which is the ploidy of the variants when it is not given.
    `poly_threshold` is the polymorphism ratio's: below that major allele
    frequency a variant is polymorphic in a population, 0.95 by default. It
    is a number from 0 to 1, both included, and NaN or anything outside is a
    ``ValueError``.

    A pass that gives no variant is a ``ValueError``, whether the source
    holds none or the steps kept none: a mean over no variant is no number.

    It mirrors pyNei's ``calc_per_var_distribs``, with these differences:
    `exp_het` of the result is the plain expected heterozygosity and
    `unbiased_exp_het` the unbiased one, where pyNei's `exp_het` holds
    whichever its ``unbiased_exp_het`` argument chose, the unbiased one by
    default, so a user who reads `exp_het` of both libraries reads two
    numbers; `stats` takes the members alone; `min_num_individuals` is
    pyNei's `min_num_samples` under the word of the glossary, and it is a
    whole number of 0 or more, where pyNei takes any number; there is no
    `num_threads`, since the threads are those of the pool of the Rust core;
    the observed heterozygosity is held to `min_num_individuals` too, which
    pyNei exempts; a block in which nothing is called gives no value, where
    pyNei gives the expected heterozygosity of such a block a 1; the
    unbiased correction is the one of the ploidy in hand, where pyNei
    applies the diploid one at every ploidy; `ploidy` is the exponent alone,
    where pyNei also counts with it the alleles the individuals are expected
    to hold; a duplicated name in a population, an empty population and an
    empty `pops` are refused; the result has `pass_stats`; and the
    populations of every statistic are in the order of the keys of `pops`,
    where pyNei sorts them for the expected heterozygosity alone.
    """
    asked_for = _the_stats(stats)
    named = _the_pops(pops)
    hist_range, num_bins, bin_type = _the_histogram(hist_kwargs)
    (
        pop_names,
        hist_bin_edges,
        obs_het,
        maf,
        exp_het,
        unbiased_exp_het,
        poly_vars_ratio,
        counts,
    ) = _core.calc_per_var_distribs(
        variants._source,
        variants._steps,
        asked_for,
        named,
        min_num_individuals,
        hist_range,
        num_bins,
        bin_type,
        ploidy,
        poly_threshold,
    )
    return PerVarDistribs(
        obs_het=_distrib_of(pop_names, hist_bin_edges, obs_het),
        maf=_distrib_of(pop_names, hist_bin_edges, maf),
        exp_het=_distrib_of(pop_names, hist_bin_edges, exp_het),
        unbiased_exp_het=_distrib_of(pop_names, hist_bin_edges, unbiased_exp_het),
        poly_vars_ratio=_poly_vars_stats_of(pop_names, poly_vars_ratio),
        pass_stats=_pass_stats_of(counts),
    )


def _the_stats(stats: Iterable[PerVarStat]) -> list[str]:
    """The names of the statistics a user asked for, each once and in the
    order they named them.

    A member of a ``StrEnum`` is a string, so one written on its own is a
    sequence of its letters: it is taken as that one statistic, as a name
    written as a string is refused.
    """
    if isinstance(stats, PerVarStat):
        stats = (stats,)
    elif isinstance(stats, str):
        # A string is a sequence of its letters, so a name written where the
        # members go would be refused for its first letter, `m`, and the
        # user would read about a statistic they never wrote.
        raise TypeError(
            f"`stats` takes the members of `PerVarStat`, and {stats!r}, a "
            f"{type(stats).__name__}, is not one of them: write "
            f"stats=(PerVarStat.MAF,) for the major allele frequency"
        )
    try:
        stats = list(stats)
    except TypeError:
        # What Python says of its own here, `'int' object is not iterable`,
        # names neither the argument nor the call.
        raise TypeError(
            f"`stats` is a sequence of the members of `PerVarStat`, and "
            f"{stats!r}, a {type(stats).__name__}, was given: write "
            f"stats=(PerVarStat.MAF,) for the major allele frequency"
        ) from None
    asked_for: list[str] = []
    for stat in stats:
        if not isinstance(stat, PerVarStat):
            raise TypeError(
                f"`stats` takes the members of `PerVarStat`, and {stat!r}, a "
                f"{type(stat).__name__}, is not one of them: write "
                f"stats=(PerVarStat.MAF,) for the major allele frequency"
            )
        if str(stat) not in asked_for:
            asked_for.append(str(stat))
    if not asked_for:
        raise ValueError(
            "`stats` names no statistic, and a result holds the ones that were "
            "asked for: leave `stats` out for the five of `PerVarStat`"
        )
    return asked_for


def _the_pops(
    pops: dict[str, Sequence[str]] | None,
) -> list[tuple[str, list[str]]] | None:
    """The populations a user named, as the pairs the Rust core takes, in the
    order the keys iterate in, and ``None`` when they named none.

    The names are not looked up here: they are resolved against the
    individuals the pass gives, which are those of the source after a filter
    of individuals when the ``Variants`` has one, and only the pass knows
    them.
    """
    if pops is None:
        return None
    if not isinstance(pops, Mapping):
        raise TypeError(
            f"`pops` is {pops!r}, a {type(pops).__name__}, and the populations are "
            f'a dict of population name to the names of its individuals, {{"pop1": '
            f'("ind00", "ind01")}}'
        )
    named = []
    for pop, individuals in pops.items():
        if isinstance(individuals, str | bytes) or not isinstance(
            individuals, Sequence
        ):
            # A string is a sequence of its letters, and one name written
            # without its comma would ask for the individuals `i`, `n`,
            # `d`... It is a `ValueError` and not the `TypeError` that ruff
            # asks for, because that is what the spec of the populations
            # gives it, and what pyNei raises for the same `pops`.
            raise ValueError(  # noqa: TRY004
                f"the individuals of the population `{pop}` are a sequence of the "
                f"names of the individuals in it, and {individuals!r}, a "
                f"{type(individuals).__name__}, was given"
            )
        for name in individuals:
            if not isinstance(name, str):
                # A `ValueError` for the reason above.
                raise ValueError(  # noqa: TRY004
                    f"the individuals of the population `{pop}` are a sequence of "
                    f"the names of the individuals in it, and {name!r}, a "
                    f"{type(name).__name__}, is not one of them"
                )
        named.append((pop, list(individuals)))
    return named


# The three keys the histogram of a statistic is given under, which are the
# ones pyNei's `_prepare_bins` reads.
_HIST_KEYS = ("range", "num_bins", "bin_type")


def _the_histogram(hist_kwargs: dict | None) -> tuple[tuple[float, float], int, str]:
    """The range, the number of bins and the kind of bins of the histogram,
    out of the dict a user gave, which is read and not changed."""
    if hist_kwargs is None:
        hist_kwargs = {}
    if not isinstance(hist_kwargs, Mapping):
        # A list of the keys is read key by key and asked for a `get` it has
        # not got, and a number cannot be iterated over at all: what Python
        # says of either names neither the argument nor the histogram.
        raise TypeError(
            f"`hist_kwargs` is {hist_kwargs!r}, a {type(hist_kwargs).__name__}, "
            f"and the histogram is a dict of `range`, the two ends, `num_bins` "
            f'and `bin_type`, {{"num_bins": 10}}'
        )
    unknown = [key for key in hist_kwargs if key not in _HIST_KEYS]
    if unknown:
        raise ValueError(
            f"{unknown[0]!r} is not a key of `hist_kwargs`, whose keys are "
            f"`range`, the two ends of the histogram, `num_bins` and `bin_type`"
        )
    hist_range = _the_range(hist_kwargs.get("range", _core.DEFAULT_HIST_RANGE))
    num_bins = hist_kwargs.get("num_bins", _core.DEFAULT_NUM_BINS)
    bin_type = hist_kwargs.get("bin_type", _core.DEFAULT_BIN_TYPE)
    return hist_range, num_bins, bin_type


def _the_range(hist_range: Sequence[float]) -> tuple[float, float]:
    """The two ends of the histogram, as the Rust core takes them.

    What is not two of something is refused here: pyo3 says of it "expected
    tuple of length 2, but got tuple of length 3", which names neither the
    argument nor what the two numbers are.
    """
    if not isinstance(hist_range, str | bytes):
        try:
            start, end = hist_range
        except TypeError, ValueError:
            pass
        else:
            return start, end
    raise TypeError(
        f"`hist_kwargs['range']` is the two ends of the histogram, (0, 1), and "
        f"{hist_range!r}, a {type(hist_range).__name__}, was given"
    )


def _distrib_of(
    pop_names: Sequence[str],
    hist_bin_edges: numpy.ndarray,
    distrib: tuple[numpy.ndarray, numpy.ndarray] | None,
) -> StatsDistrib | None:
    """The distribution of one statistic, or ``None`` when nobody asked for
    it: the mean and the histogram counts the core gave, under the names of
    the populations."""
    if distrib is None:
        return None
    mean, hist_counts = distrib
    return StatsDistrib(
        mean=pandas.Series(mean, index=list(pop_names)),
        hist_bin_edges=hist_bin_edges,
        hist_counts=pandas.DataFrame(hist_counts, columns=list(pop_names)),
    )


def _poly_vars_stats_of(
    pop_names: Sequence[str],
    poly_vars_ratio: tuple[
        numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray
    ]
    | None,
) -> PolyVarsStats | None:
    """The counts of the polymorphism ratio, or ``None`` when nobody asked
    for them."""
    if poly_vars_ratio is None:
        return None
    num_poly, num_variable, with_data, poly_ratio, over_variables = poly_vars_ratio
    names = list(pop_names)
    return PolyVarsStats(
        num_poly=pandas.Series(num_poly, index=names),
        poly_ratio=pandas.Series(poly_ratio, index=names),
        poly_ratio_over_variables=pandas.Series(over_variables, index=names),
        num_variable=pandas.Series(num_variable, index=names),
        tot_num_variants_with_data=pandas.Series(with_data, index=names),
    )
