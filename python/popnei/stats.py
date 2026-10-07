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

:func:`popnei.calc_per_individual_stats` makes a pass of its own and gives
two numbers for each individual instead: the share of the variants at which
its genotype is missing and the share of its called genotypes at which it is
heterozygous. It takes no `pops`, since each of its values is of one
individual.

:func:`popnei.calc_var_density` takes no `pops` either: it counts the
variants in windows along each chromosome, and reads no genotype.
"""

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum

import numpy
import pandas

from popnei import _core
from popnei.variant import PassStats, Variants, _pass_stats_of


class PerVarStat(StrEnum):
    """The six statistics :func:`popnei.calc_per_var_distribs` calculates.

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

    MISSING_RATE = "missing_rate"
    """The missing genotypes of a population over its individuals, called or
    not, a half called genotype being missing."""


@dataclass(frozen=True)
class StatsDistrib:
    """The distribution of one statistic over the variants of a pass, per
    population.

    A variant has no value of a statistic held to `min_num_individuals` in a
    population when the population has too little data at it, and such a
    variant is out of the mean and in no bin, so the histograms of two
    populations can count different numbers of variants. The missing rate is
    not held to it: every variant has one in every population, and its mean
    is never NaN. A value outside the range of the bins is in the mean and in
    no bin, as :func:`numpy.histogram` leaves it out too.
    """

    mean: pandas.Series
    """The mean over the variants that had a value, one value per
    population, NaN for a population in which no variant had one, which
    the missing rate never is."""

    hist_bin_edges: numpy.ndarray
    """The edges of the bins, one more than there are bins.

    The five distributions of one result share this array, as pyNei's do, so
    it is read only: a number written into the edges of one statistic would
    be in the edges of the other four."""

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

    missing_rate: StatsDistrib | None
    """The distribution of the missing rate, which every variant has in every
    population."""

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
    """Up to six statistics of every variant and every population, in one
    pass over `variants`, as a mean and a histogram each.

    The statistics are the observed heterozygosity, the heterozygous
    genotypes of a population over its called ones; the major allele
    frequency, the count of the commonest allele over the called alleles;
    the expected heterozygosity, the chance that gene copies taken at random
    from the population are not all of the same allele, plain and corrected
    for the frequencies being estimated from the copies the statistic is
    computed over; the polymorphism ratio, how many of the variants vary in
    the population, which is a count and not a distribution; and the missing
    rate, the missing genotypes of the population over its individuals,
    called or not, a half called genotype being missing. The missing rate
    has a value at every variant, whatever `min_num_individuals` is, and a
    variant with nothing called in a population has a rate of 1 there.

    It is a consumer of the `variants`: it makes one pass over the source
    through the steps the ``Variants`` has when it is called, and the
    ``Variants`` is as it was afterwards.

    `stats` says which of the six to calculate, as members of
    :class:`popnei.PerVarStat`, all six by default. Anything that is not a
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
    than the threshold, for every statistic but the missing rate. A variant
    with no value in a population is out of the mean and in no bin of the
    histogram of that population.

    `hist_kwargs` is the histogram, under four keys: ``range``, the two
    ends, ``(0, 1)`` by default, which is where the statistics live;
    ``num_bins``, 40 by default; ``bin_type``, ``"linear"`` for bins of
    equal width or ``"logarithmic"`` for bins of equal ratio, whose
    ``range`` has to start above 0; and ``closed``, the edge each bin holds,
    ``"left"`` by default or ``"right"``. The dict is read and not changed,
    and a key that is none of the four is a ``ValueError``, since it would
    leave the default in its place with nothing said.

    The k-th of the n + 1 edges of bins of equal width is (start · (n − k)
    + end · k) / n, which is the float64 of the decimal k / n when the two
    ends are whole numbers or halves: of 1000 bins from 0 to 1 the edge 7 is
    the number ``0.007`` is read as. :func:`numpy.linspace`, which pyNei's
    edges are, gives some of them one unit in the last place away. Edges
    that do not go up, as those of very many bins over a range narrower
    than the rounding of a float64, are a ``ValueError``.

    With ``closed`` ``"left"`` a value falls in the bin whose left edge is
    at most the value and whose right edge is above it, and the last bin
    takes its right edge too, as :func:`numpy.histogram` does. With
    ``"right"`` a value falls in the bin whose left edge is below it and
    whose right edge is at least the value, and the first bin takes its left
    edge too, as :func:`pandas.cut` with ``right=True,
    include_lowest=True``: the bins below an edge t then count the variants
    whose value is at most t, which is what ``filter_by_missing_data(t)``,
    ``filter_by_maf(t)`` and ``filter_by_obs_het(t)`` keep. A value outside
    the range is in no bin and in the mean.

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
    where pyNei sorts them for the expected heterozygosity alone; and pyNei
    has no missing rate.
    """
    if not isinstance(variants, Variants):
        # The path of the VCF whose variants are read is the mistake that is
        # easiest to make, and what it gave was the `AttributeError` of an
        # object with no source inside it.
        raise TypeError(
            f"`variants` is {variants!r}, a {type(variants).__name__}, and "
            f"`calc_per_var_distribs` reads the variants of a source: give it what "
            f"`open_vcf` or `open_vars` gives, "
            f"calc_per_var_distribs(open_vcf(vcf_path))"
        )
    asked_for = _the_stats(stats)
    named = _the_pops(pops)
    hist_range, num_bins, bin_type, closed = _the_histogram(hist_kwargs)
    (
        pop_names,
        hist_bin_edges,
        obs_het,
        maf,
        exp_het,
        unbiased_exp_het,
        poly_vars_ratio,
        missing_rate,
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
        closed,
        ploidy,
        poly_threshold,
    )
    return PerVarDistribs(
        obs_het=_distrib_of(pop_names, hist_bin_edges, obs_het),
        maf=_distrib_of(pop_names, hist_bin_edges, maf),
        exp_het=_distrib_of(pop_names, hist_bin_edges, exp_het),
        unbiased_exp_het=_distrib_of(pop_names, hist_bin_edges, unbiased_exp_het),
        poly_vars_ratio=_poly_vars_stats_of(pop_names, poly_vars_ratio),
        missing_rate=_distrib_of(pop_names, hist_bin_edges, missing_rate),
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
            f"{type(stats).__name__}, is not one of them: {_what_to_write(stats)}"
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
                f"{type(stat).__name__}, is not one of them: {_what_to_write(stat)}"
            )
        if str(stat) not in asked_for:
            asked_for.append(str(stat))
    if not asked_for:
        raise ValueError(
            "`stats` names no statistic, and a result holds the ones that were "
            "asked for: leave `stats` out for the six of `PerVarStat`"
        )
    return asked_for


def _what_to_write(given: object) -> str:
    """What a message that refuses `given` in `stats` tells the user to
    write: the member whose value it is, when it is the name of one written
    as a string, and the major allele frequency as an example otherwise."""
    if isinstance(given, str):
        try:
            member = PerVarStat(given)
        except ValueError:
            pass
        else:
            return f"write stats=(PerVarStat.{member.name},)"
    return "write stats=(PerVarStat.MAF,) for the major allele frequency"


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


# The four keys the histogram of a statistic is given under: the three
# pyNei's `_prepare_bins` reads, and the edge each bin holds.
_HIST_KEYS = ("range", "num_bins", "bin_type", "closed")


def _the_histogram(
    hist_kwargs: dict | None,
) -> tuple[tuple[float, float], int, str, str]:
    """The range, the number of bins, the kind of bins and the edge each bin
    holds, out of the dict a user gave, which is read and not changed."""
    if hist_kwargs is None:
        hist_kwargs = {}
    if not isinstance(hist_kwargs, Mapping):
        # A list of the keys is read key by key and asked for a `get` it has
        # not got, and a number cannot be iterated over at all: what Python
        # says of either names neither the argument nor the histogram.
        raise TypeError(
            f"`hist_kwargs` is {hist_kwargs!r}, a {type(hist_kwargs).__name__}, "
            f"and the histogram is a dict of `range`, the two ends, `num_bins`, "
            f'`bin_type` and `closed`, {{"num_bins": 10}}'
        )
    unknown = [key for key in hist_kwargs if key not in _HIST_KEYS]
    if unknown:
        raise ValueError(
            f"{unknown[0]!r} is not a key of `hist_kwargs`, whose keys are "
            f"`range`, the two ends of the histogram, `num_bins`, `bin_type` and "
            f"`closed`"
        )
    hist_range = _the_range(hist_kwargs.get("range", _core.DEFAULT_HIST_RANGE))
    num_bins = hist_kwargs.get("num_bins", _core.DEFAULT_NUM_BINS)
    bin_type = hist_kwargs.get("bin_type", _core.DEFAULT_BIN_TYPE)
    closed = hist_kwargs.get("closed", _core.DEFAULT_CLOSED)
    return hist_range, num_bins, bin_type, closed


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


@dataclass(frozen=True)
class PerIndividualStats:
    """What :func:`popnei.calc_per_individual_stats` gives back.

    The two series are indexed by the names of the individuals the pass
    gave, in its order, which is the order of the source unless a filter of
    individuals named them in another one.
    """

    missing_gt_rate: pandas.Series
    """The variants at which the genotype of the individual is missing, a
    half called genotype among them, over the variants of the pass."""

    obs_het_rate: pandas.Series
    """The variants at which the genotype of the individual is called and
    its alleles are not all the same, over its called genotypes, NaN for an
    individual that called none of them."""

    pass_stats: PassStats
    """How many variants the pass gave, after the steps of the ``Variants``,
    and what each filter of it was given and kept."""


def calc_per_individual_stats(variants: Variants) -> PerIndividualStats:
    """The missing rate and the heterozygosity rate of every individual, in
    one pass over `variants`.

    The missing rate is the share of the variants at which the individual
    has no genotype, and a half called genotype is missing and not
    heterozygous. The heterozygosity rate is the share of its called
    genotypes at which its alleles are not all the same. The first tells a
    user which individuals were badly genotyped, and the second which ones
    are more heterozygous than the rest, a sign of a mixed sample or of an
    outcrossed individual among inbred ones. An individual that called no
    genotype has a missing rate of 1 and no heterozygosity rate, NaN.

    It is a consumer of the `variants`: it makes one pass over the source
    through the steps the ``Variants`` has when it is called, and the
    ``Variants`` is as it was afterwards. A pass that gives no variant is a
    ``ValueError``, whether the source holds none or the steps kept none.

    It mirrors pyNei's ``calc_per_sample_stats``, with these differences:
    the heterozygosity rate divides by the called genotypes of the
    individual, where pyNei divides by every variant, so an individual with
    more missing data looks less heterozygous there, and popnei's number is
    what plink2's ``--het`` gives, with the missing rate beside it saying
    what pyNei's one number said; the result is this dataclass and not a
    pandas frame of the two columns, since every result of a consumer
    carries its `pass_stats`, which a frame has no place for; and there is
    no `num_threads`, since the threads are those of the pool of the Rust
    core.
    """
    if not isinstance(variants, Variants):
        # The path of the VCF whose variants are read is the mistake that is
        # easiest to make, and what it gave was the `AttributeError` of an
        # object with no source inside it.
        raise TypeError(
            f"`variants` is {variants!r}, a {type(variants).__name__}, and "
            f"`calc_per_individual_stats` reads the variants of a source: give it what "
            f"`open_vcf` or `open_vars` gives, "
            f"calc_per_individual_stats(open_vcf(vcf_path))"
        )
    individuals, missing_gt_rate, obs_het_rate, counts = (
        _core.calc_per_individual_stats(variants._source, variants._steps)
    )
    names = list(individuals)
    return PerIndividualStats(
        missing_gt_rate=pandas.Series(missing_gt_rate, index=names),
        obs_het_rate=pandas.Series(obs_het_rate, index=names),
        pass_stats=_pass_stats_of(counts),
    )


@dataclass(frozen=True)
class VarDensity:
    """What :func:`popnei.calc_var_density` gives back."""

    windows: pandas.DataFrame
    """One row for each window: `chrom`, the name of its chromosome,
    `start` and `end`, its first and its last position, both included and
    counted from 1, of the dtype ``uint64``, and `num_vars`, how many
    variants of the pass are at a position from `start` to `end`, of the
    dtype ``uint32``. The chromosomes are in the order of the lengths, those
    of `chrom_lengths` or of the source, and after them those with variants
    and no length, in the order their first variant came; the windows of
    each in the order of their positions."""

    pass_stats: PassStats
    """How many variants the pass gave, after the steps of the ``Variants``,
    and what each filter of it was given and kept."""


def calc_var_density(
    variants: Variants,
    window_size: int,
    chrom_lengths: Mapping[str, int] | None = None,
) -> VarDensity:
    """How many variants fall in each window of `window_size` base pairs
    along each chromosome, in one pass over `variants`.

    A user sees with it where the variants are crowded, where there are
    none, a centromere or a region that did not map, and how evenly a filter
    took variants out. The windows of a chromosome are laid end to end from
    the position 1 and do not overlap: window k, counted from 0, holds the
    positions from k x `window_size` + 1 to (k + 1) x `window_size`. A
    window with no variant is in the result with a count of 0.

    With the length of a chromosome the windows cover it to its end, and
    the last one ends at the length, so it is shorter than the others when
    the length is not a multiple of the width. Without one, the windows go up
    to the one that holds the last variant of the chromosome, and that one
    ends at its full width. The lengths are those of `chrom_lengths`, a
    mapping of chromosome name to length, and when it is None those of the
    source: the ``##contig`` lines of a VCF that have a ``length``, and what
    a vars file keeps of them. A `chrom_lengths` that is given replaces the
    lengths of the source for every chromosome, and one it does not name has
    no length. A chromosome with a length is in the result whether or not it
    has a variant.

    It is a consumer of the `variants`: it makes one pass over the source
    through the steps the ``Variants`` has when it is called, reading only
    the chromosome and the position of each variant, and the ``Variants`` is
    as it was afterwards. The variants need not be sorted.

    `window_size` and each length are whole numbers of 1 or more: what is no
    whole number is a ``TypeError`` and one below 1 a ``ValueError``, which
    names the argument. These are a ``ValueError`` too: a variant past the
    length of its chromosome, whose message says whether the length came
    from `chrom_lengths` or from the source; a variant at the position 0,
    which the VCF format allows for a telomere and which lies in no window;
    more than 10 million windows over all the chromosomes, which asks for a
    wider window; and a pass that gives no variant, whatever the lengths.

    pyNei has no density of the variants.
    """
    if not isinstance(variants, Variants):
        # The path of the VCF whose variants are read is the mistake that is
        # easiest to make, and what it gave was the `AttributeError` of an
        # object with no source inside it.
        raise TypeError(
            f"`variants` is {variants!r}, a {type(variants).__name__}, and "
            f"`calc_var_density` reads the variants of a source: give it what "
            f"`open_vcf` or `open_vars` gives, "
            f"calc_var_density(open_vcf(vcf_path), 100000)"
        )
    lengths = None if chrom_lengths is None else _the_chrom_lengths(chrom_lengths)
    names, windows_per_chrom, start, end, num_vars, counts = _core.calc_var_density(
        variants._source, variants._steps, window_size, lengths
    )
    chroms = numpy.repeat(numpy.array(list(names), dtype=object), windows_per_chrom)
    windows = pandas.DataFrame(
        {"chrom": chroms, "start": start, "end": end, "num_vars": num_vars}
    )
    return VarDensity(windows=windows, pass_stats=_pass_stats_of(counts))


def _the_chrom_lengths(chrom_lengths: object) -> list[tuple[str, object]]:
    """The pairs of chromosome name and length of `chrom_lengths`, in the
    order it iterates in, which is the order of the chromosomes of the
    result, each a str and a whole number of 0 or more that 64 bits hold."""
    if not isinstance(chrom_lengths, Mapping):
        raise TypeError(
            f"`chrom_lengths` is {chrom_lengths!r}, a {type(chrom_lengths).__name__}, "
            f"and it is a mapping of chromosome name to length, "
            f"{{'chr1': 248956422}}"
        )
    pairs = list(chrom_lengths.items())
    for chrom, length in pairs:
        if not isinstance(chrom, str):
            raise TypeError(
                f"`chrom_lengths` names the chromosome {chrom!r}, a "
                f"{type(chrom).__name__}, and a chromosome is named by a str, as "
                f"the source writes it"
            )
        # The binding crate refuses the same lengths, under the name of the
        # argument alone: here the message names the chromosome too, which
        # is what a user of a dict of 25 of them looks for. A length of 0 is
        # the core's to refuse, and its message names the chromosome.
        where = f"`chrom_lengths[{chrom!r}]`"
        if isinstance(length, bool) or not isinstance(length, int):
            raise TypeError(
                f"{where} is {length!r}, a {type(length).__name__}, and a length is "
                f"a whole number of 1 or more"
            )
        if length < 0 or length > _LARGEST_LENGTH:
            raise ValueError(
                f"{where} is {length}, and a length is a whole number of 1 or more "
                f"that 64 bits hold"
            )
    return pairs


# The largest length a chromosome is given, the largest whole number of 64
# bits, which the VCF reader takes in a `##contig` line too.
_LARGEST_LENGTH = 2**64 - 1
