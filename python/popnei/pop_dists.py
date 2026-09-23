"""The distances between populations, and the result they come in.

Two populations are far apart when the alleles of their individuals are not
the same alleles in the same proportions, and there are seven measures of
how far, each answering a different question.
:func:`popnei.calc_pop_dists` calculates the ones a user asks for in one
pass over the variants, for every pair of the populations they name, and
gives them in a :class:`PopDists`: one :class:`popnei.Distances` for each
measure, with the standard error of each pair beside its value.

The seven are members of :class:`PopDistMeasure`. Two of them are
calculated today, Hudson's F_ST and f_2; the other five raise, and the work
packages 2 and 3 of `docs/plans/dists-pops.md` add them.

`docs/specs/dists.md` has each measure, what it answers, the program it is
verified against and the numbers the tests assert.
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import Literal

import numpy

from popnei import _core
from popnei.dists import Distances
from popnei.stats import _the_pops
from popnei.variant import PassStats, Variants, _pass_stats_of


class PopDistMeasure(StrEnum):
    """The seven measures of how far apart two populations are.

    The value of each member is the name of the field of :class:`PopDists`
    that holds its result.
    """

    FST = "fst"
    """Hudson's F_ST, how much of the diversity of two populations taken
    together lies between them rather than within them, from 0 for two with
    the same allele frequencies everywhere to 1 for two that share no allele
    anywhere. It is what most work on single nucleotide polymorphisms
    reports."""

    F2 = "f2"
    """How much allele frequency the two populations have drifted apart by,
    in the units it was measured in rather than divided by the diversity
    they hold, which is what makes it add up along a tree and what admixture
    graphs are built from."""

    CHORD = "chord"
    """The chord distance of Cavalli-Sforza and Edwards, one of the two
    measures here that are Euclidean."""

    DA = "da"
    """Nei's D_A, the square of the chord distance."""

    DEST = "dest"
    """Jost's D, how different the alleles the two populations hold are,
    which is the one measure of the seven that pyNei has."""

    GST = "gst"
    """Nei's G_ST, the share of the diversity of the two populations that
    lies between them, corrected for the individuals they were estimated
    from."""

    GST_STANDARDIZED = "gst_standardized"
    """G_ST divided by the largest value it could reach with the diversity
    the two populations hold, so that it reaches 1."""


# The measures that have a value today. Hudson's F_ST and f_2 are what work
# package 1 of `docs/plans/dists-pops.md` calculates; its work packages 2 and
# 3 add the other five, and asking for one of those is refused here until
# they do, so that nobody reads a vector of NaN as a distance. Removing the
# refusal is this tuple and the `if` in `_the_measures` that reads it.
_MEASURES_THAT_HAVE_A_VALUE = (PopDistMeasure.FST, PopDistMeasure.F2)


@dataclass(frozen=True, eq=False, repr=False)
class PopDists:
    """The measures of how far apart the populations are, one
    :class:`popnei.Distances` for each that was asked for.

    Every measure is over the same pairs, in the order (0, 1), (0, 2), ...,
    (1, 2), ..., over the populations in the order :attr:`pops` has them,
    which is the order of the `pops` dict that was given.

    ``dists == other`` is true for the same object and false for any other,
    as it is for a :class:`popnei.Distances`: two results are not compared
    value by value.
    """

    pops: tuple[str, ...]
    """The names of the populations, in the order of the `pops` dict, which
    is the order of the pairs of every measure."""

    num_vars: numpy.ndarray
    """How many variants counted for each pair, a read only int64 array in
    the order of the pairs.

    A variant counts for a pair when both of its populations have at least
    `min_num_individuals` called genotypes there, so a variant one
    population is short of is lost by the pairs that population is in and
    kept by the others: two pairs are means over different variants, which
    is why each pair carries its own count.
    """

    pass_stats: PassStats
    """The counts of the pass the calculation made: how many variants it
    took, after the steps, and how many each filter of the ``Variants`` was
    given and kept."""

    fst: Distances | None = None
    """Hudson's F_ST of every pair, and ``None`` when it was not asked
    for."""

    f2: Distances | None = None
    """f_2 of every pair, and ``None`` when it was not asked for."""

    chord: Distances | None = None
    """The chord distance of every pair, and ``None`` when it was not asked
    for."""

    da: Distances | None = None
    """Nei's D_A of every pair, and ``None`` when it was not asked for."""

    dest: Distances | None = None
    """Jost's D of every pair, and ``None`` when it was not asked for."""

    gst: Distances | None = None
    """Nei's G_ST of every pair, and ``None`` when it was not asked for."""

    gst_standardized: Distances | None = None
    """The standardized G''_ST of every pair, and ``None`` when it was not
    asked for."""

    f2_groups: numpy.ndarray | None = None
    """f_2 within each resampling group, a read only float64 array of groups
    x pairs, and ``None`` when no standard errors were asked for.

    f_3 and f_4, the statistics of three and of four populations that
    admixture graphs are fitted with, are sums and differences of these, so
    they can be built from this array without reading the genotypes again.
    """

    group_ids: tuple[tuple[str, int, int], ...] = ()
    """One ``(chrom, start, end)`` for each resampling group, in the order
    the groups were started: the chromosome of its variants and the position
    of its first and of its last one, both included. It is empty when no
    standard errors were asked for."""

    def __repr__(self) -> str:
        """The populations and the measures that were calculated, and not the
        values of any of them."""
        calculated = [
            measure for measure in PopDistMeasure if getattr(self, measure) is not None
        ]
        return (
            f"<PopDists of {len(self.pops)} populations, "
            f"{', '.join(calculated) if calculated else 'no measure'}>"
        )


def calc_pop_dists(
    variants: Variants,
    pops: Mapping[str, Sequence[str]],
    jackknife_group: int | Literal["variant"] | None,
    measures: Sequence[PopDistMeasure] | None = None,
    min_num_individuals: int = _core.DEFAULT_MIN_NUM_INDIVIDUALS,
) -> PopDists:
    """How far apart every pair of the populations of `pops` is, in one pass
    over `variants`.

    The measures are functions of the same three counts of a population at a
    variant, how often each allele was called there, how many genotypes were
    called whole and how many of those are heterozygous, so a user who wants
    to compare two of them pays for one reading of the variants and not two.
    Every allele counts as itself: a variant of three alleles is not
    collapsed to the commonest one against the rest, which is what lets the
    same call serve microsatellites and single nucleotide polymorphisms.

    It is a consumer of the `variants`: it makes one pass over the source
    through the steps the ``Variants`` has when it is called, and the
    ``Variants`` is as it was afterwards.

    `pops` is a dict of population name to the names of its individuals,
    which are looked up among :attr:`popnei.Variants.individuals`, the ones
    the pass gives. Fewer than two populations is a ``ValueError``, since
    every measure is of a pair. An individual in two populations is taken
    and counted in each of them, and one in none takes no part.

    `jackknife_group` is how the variants are cut into the groups the
    standard errors are resampled over, and it has no default: a number is a
    length in base pairs of a chromosome, ``"variant"`` makes each variant a
    group of its own, and ``None`` asks for no standard errors.

    A group has to be longer than the distance over which two variants still
    carry the same history, because two groups that share it are not the
    independent draws the standard error takes them for, and there have to
    be at least 20 groups, which is a ``ValueError`` below. The number of the
    f-statistics literature is 5 centimorgans, about 5 million base pairs in
    humans, and it does not carry over by itself: linkage disequilibrium runs
    6.1 to 12.5 centimorgans in cultivated tomato and falls off within 18
    thousand base pairs in its wild relative *S. pimpinellifolium*, and one
    length in base pairs is several different lengths in centimorgans along
    one chromosome. A user who does not know the decay distance of their own
    panel measures it from the curve of r^2 against distance. A few hundred
    microsatellite loci scattered over a genome have no linkage to speak of
    and take ``"variant"``; a panel of linked variants must not use it.

    `measures` says which of the seven of :class:`PopDistMeasure` to
    calculate, and ``None`` is all of them, since the pass is what costs and
    each measure is a division at the end of it. Five of the seven are not
    calculated yet, so asking for one of them, and asking for all of them
    with ``None``, is a ``ValueError`` that names the two that are: the work
    packages 2 and 3 of `docs/plans/dists-pops.md` add the rest.

    `min_num_individuals` is how many called genotypes a population needs at
    a variant for that variant to count for a pair, 20 by default. The test
    is made for each pair on its own, so a population with fewer loses that
    variant in every pair it is in and the pairs it is not in keep it.

    What it gives is a :class:`PopDists` with one
    :class:`popnei.Distances` for each measure that was asked for, under the
    name of the measure, each carrying the value of every pair and, where
    groups were asked for, its standard error; the names of the populations;
    how many variants counted for each pair; f_2 within each group; the
    chromosome and the positions of each group; and the counts of the pass.

    A pair that counted no variant has NaN for every measure, which is not
    an error: the other pairs may have values. f_2 and F_ST can come out
    negative, for one variant and for a whole dataset, and popnei does not
    clamp them: it is the correction doing its work, and a user who sees a
    small negative F_ST has two populations this dataset cannot tell apart.

    A pass that gives no variant is a ``ValueError``, whether the source has
    none or the steps kept none: the message says which of the two, and,
    when it was the steps, what each filter was given and kept.

    It mirrors pyNei's ``calc_jost_dest_pop_dists``, which calculates Jost's
    D and nothing else, with these differences: popnei has one function for
    all the measures, because the pass is the cost and the counts are
    shared; the populations are in the order of the `pops` dict, where pyNei
    sorts their names, which changes no value; there is no `num_threads`,
    since no calculation of popnei has it; there is no `alleles`, since
    popnei counts the alleles each variant has and lines nothing up across
    blocks; and `min_num_individuals` is pyNei's `min_num_samples` under the
    word of the glossary.
    """
    if not isinstance(variants, Variants):
        # The path of the VCF is the mistake that is easiest to make, and
        # what it gave was the `AttributeError` of an object with no source
        # inside it.
        raise TypeError(
            f"`variants` is {variants!r}, of the type "
            f"`{type(variants).__name__}`, and the distances between "
            f"populations are calculated over the variants of a source: give "
            f"it what `open_vcf` or `open_vars` gives, "
            "calc_pop_dists(open_vcf(vcf_path), pops, jackknife_group=None)"
        )
    asked_for = _the_measures(measures)
    named = _the_pops(pops)
    if named is None:
        raise TypeError(
            "`pops` is `None`, and the distances between populations are "
            "calculated for the populations a user names: give a dict of "
            'population name to the names of its individuals, {"pop1": '
            '("ind00", "ind01")}'
        )
    (
        pop_names,
        of_each_measure,
        num_vars,
        f2_groups,
        group_ids,
        counts,
    ) = _core.calc_pop_dists(
        variants._source,
        variants._steps,
        named,
        asked_for,
        jackknife_group,
        min_num_individuals,
    )
    pops_of_the_pass = tuple(pop_names)
    of_its_name = {
        measure: Distances(
            dist_vector=values,
            names=pops_of_the_pass,
            standard_errors=standard_errors,
        )
        for measure, (values, standard_errors) in zip(
            asked_for, of_each_measure, strict=True
        )
    }
    return PopDists(
        pops=pops_of_the_pass,
        num_vars=num_vars,
        pass_stats=_pass_stats_of(counts),
        f2_groups=f2_groups,
        group_ids=tuple(group_ids),
        **of_its_name,
    )


def _the_measures(measures: Sequence[PopDistMeasure] | None) -> list[str]:
    """The names of the measures a user asked for, each once and in the order
    they named them, which is ``None`` for all seven.

    A member of a ``StrEnum`` is a string, so one written on its own is a
    sequence of its letters: it is taken as that one measure, as a name
    written as a string is.
    """
    if measures is None:
        measures = tuple(PopDistMeasure)
    elif isinstance(measures, str):
        # A name written where the members go is taken as the measure of
        # that name, and one that is of no measure is refused by its name
        # and not letter by letter.
        measures = (measures,)
    try:
        measures = list(measures)
    except TypeError:
        # What Python says of its own here, `'int' object is not iterable`,
        # names neither the argument nor the call.
        raise TypeError(
            f"`measures` is a sequence of the members of `PopDistMeasure`, "
            f"and {measures!r}, a {type(measures).__name__}, was given: write "
            f"measures=(PopDistMeasure.FST,) for Hudson's F_ST"
        ) from None
    asked_for: list[str] = []
    for measure in measures:
        if not isinstance(measure, str):
            raise TypeError(
                f"`measures` takes the members of `PopDistMeasure`, and "
                f"{measure!r}, a {type(measure).__name__}, is not one of "
                f"them: write measures=(PopDistMeasure.FST,) for Hudson's "
                f"F_ST"
            )
        if measure not in tuple(PopDistMeasure):
            raise ValueError(
                f"`{measure}` is not one of the measures of how far apart two "
                f"populations are, which are {_named(tuple(PopDistMeasure))}"
            )
        if measure not in asked_for:
            asked_for.append(str(measure))
    if not asked_for:
        raise ValueError(
            "`measures` names no measure, and a result holds the ones that "
            "were asked for: leave `measures` out for every measure there is"
        )
    # The five the work packages 2 and 3 of `docs/plans/dists-pops.md` add
    # have no value yet, and a vector of NaN says nothing about itself.
    not_written_yet = [
        measure for measure in asked_for if measure not in _MEASURES_THAT_HAVE_A_VALUE
    ]
    if not_written_yet:
        raise ValueError(
            f"{_named(not_written_yet)} "
            f"{'is' if len(not_written_yet) == 1 else 'are'} not calculated "
            f"yet, and what popnei calculates today is "
            f"{_named(_MEASURES_THAT_HAVE_A_VALUE)}: ask for those"
        )
    return asked_for


def _named(measures: Sequence[str]) -> str:
    """`measures` in one sentence, each in backticks, the last one after an
    "and"."""
    named = [f"`{measure}`" for measure in measures]
    if len(named) == 1:
        return named[0]
    return f"{', '.join(named[:-1])} and {named[-1]}"
