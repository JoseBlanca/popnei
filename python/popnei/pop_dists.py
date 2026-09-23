"""The distances between populations, and the result they come in.

Two populations are far apart when the alleles of their individuals are not
the same alleles in the same proportions, and there are seven measures of
how far, each answering a different question.
:func:`popnei.calc_pop_dists` calculates the ones a user asks for in one
pass over the variants, for every pair of the populations they name, and
gives them in a :class:`PopDists`: one :class:`popnei.Distances` for each
measure, with the standard error of each pair beside its value.

The seven are members of :class:`PopDistMeasure`, and all seven are
calculated: Hudson's F_ST, f_2, the chord distance, Nei's D_A, Jost's D,
Nei's G_ST and the standardized G''_ST.

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
    """Jost's D, how much of the allelic variety of the two populations is
    not shared, from 0 when they hold the same alleles at the same
    frequencies to 1 when they share none. It is the one measure of the
    seven that pyNei has, and it is the one to read on microsatellites,
    where G_ST below cannot reach 1."""

    GST = "gst"
    """Nei's G_ST, the share of the diversity of the two populations that
    lies between them, corrected for the individuals they were estimated
    from. With two populations it cannot pass (1 - H_S)/(1 + H_S), with H_S
    the mean corrected diversity within them, so two internally diverse
    populations that share no allele still give a small number."""

    GST_STANDARDIZED = "gst_standardized"
    """The G''_ST of Meirmans and Hedrick (2011), G_ST rescaled so that it
    reaches 1 when the two populations share no allele whatever their
    diversity. It is not Hedrick's earlier G'_ST, which divides G_ST by the
    largest value it could take: for the first two populations of the
    biallelic panel of the tests this one is 0.1620 and G'_ST is 0.1155.

    G'_ST is ``gst * (1 + h_s) / (1 - h_s)``, with ``h_s`` the mean
    corrected diversity within the two populations over the variants that
    counted for that pair. No field of :class:`PopDists` holds ``h_s``, and
    :attr:`PopDists.gst` and :attr:`PopDists.dest` together give it back::

        d = 1 / (1 / gst - 1 + 2 / dest)
        h_s = 1 - 2 * d / dest

    which for that same pair is an ``h_s`` of 0.351109 and a G'_ST of
    0.115481. The unbiased expected heterozygosity of
    :func:`popnei.calc_per_var_distribs` is another quantity, corrected for
    the sample in another way and taken over one population and not the two
    pooled: on that pair its mean is 0.351160, four digits of agreement and
    a G'_ST nobody would see was wrong.
    """


# The measures that have a value, which the core crate holds and which is all
# seven of them. Asking for one that has none is refused, so that nobody
# reads a vector of NaN as a distance, and there is nothing left to refuse.
# It is read from the core and not written here so that a measure is added in
# one place and not in this package, in the TypeScript one and in the core.
_MEASURES_THAT_HAVE_A_VALUE = tuple(_core.pop_dist_measures_that_have_a_value())


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
    for.

    It is the form ``adegenet::dist.genpop(method = 2)`` gives, the chord of
    the sphere of radius 1 divided by the square root of 2: two populations
    that share no allele are 1 apart here and 1.414 apart unscaled. Books
    normalize it in several ways, so a number compared with another program
    has to be compared with the same form, and the scaling changes nothing
    for a tree or for a principal coordinate analysis. :attr:`da` is the
    square of this form and is Nei's D_A.
    """

    da: Distances | None = None
    """Nei's D_A of every pair, and ``None`` when it was not asked for."""

    dest: Distances | None = None
    """Jost's D of every pair, and ``None`` when it was not asked for.

    It is the D_est of Jost (2008) under the correction of Nei and Chesser
    (1983) for the individuals it was estimated from, which is the estimator
    pyNei computes and the one GenAlEx prints. mmod's ``pairwise_D`` in R
    computes another estimator of the same quantity, leaving the observed
    heterozygosity out of the correction and dividing by 2n - 1 where this
    one divides by n - 1, so it gives another number: on the biallelic panel
    of the tests, 1200 variants of three populations of 48 to 84
    individuals, the two are 7.3e-5 apart at the furthest, and on the
    multiallelic one, 120 microsatellite loci of six alleles in three
    populations of 30, 3.5e-4. A user who finds a third number in some
    program is holding a third estimator, which is what the Jost's D item of
    `docs/specs/dists.md` writes both formulas out for.

    At a ploidy of 1 it has no value, NaN for every pair: it is a ratio of
    the corrected H_S and H_T, which raise the allele frequencies of a
    variant to the ploidy and take the sum from 1, so both are 0 there by
    their definitions and what would be divided is the residue of the
    rounding of the frequencies. :attr:`fst`, :attr:`f2`, :attr:`chord` and
    :attr:`da` are given at a ploidy of 1 as at any other.
    """

    gst: Distances | None = None
    """Nei's G_ST of every pair, and ``None`` when it was not asked for. It
    comes from the same two corrected means as :attr:`dest`, so mmod's
    ``pairwise_Gst_Nei`` carries the same difference of estimator, 7.2e-5 at
    the furthest on the biallelic panel and 9.0e-5 on the multiallelic one.
    It has no value at a ploidy of 1 either, for the reason :attr:`dest`
    gives."""

    gst_standardized: Distances | None = None
    """The standardized G''_ST of every pair, and ``None`` when it was not
    asked for. It is Meirmans and Hedrick's (2011) and not Hedrick's earlier
    G'_ST, and mmod's ``pairwise_Gst_Hedrick``, which computes this one
    whatever its name suggests, is 1.9e-4 from it at the furthest on the
    biallelic panel and 4.7e-4 on the multiallelic one, by the same
    difference of estimator as :attr:`dest`. It has no value at a ploidy of
    1 either, for the reason :attr:`dest` gives."""

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

    A length needs the variants of each chromosome to come together and in
    order of position, and a source whose variants go back is a
    ``ValueError`` that names the chromosome and the two positions: the cut
    compares the position of a variant with the first position of the group
    being filled, so a variant that goes back joins that group instead of
    starting one and the groups are not the stretches that were asked for.
    ``"variant"`` and ``None`` take a source in any order.

    A group has to be longer than the distance over which two variants still
    carry the same history, because two groups that share it are not the
    independent draws the standard error takes them for, and there have to
    be at least 20 groups: a pass whose variants fall into fewer is a
    ``ValueError`` that says how many they fell into. The number of the
    f-statistics literature is 5 centimorgans, about 5 million base pairs in
    humans, and it does not carry over by itself: linkage disequilibrium runs
    6.1 to 12.5 centimorgans in cultivated tomato and falls off within 18
    thousand base pairs in its wild relative *S. pimpinellifolium*, and one
    length in base pairs is several different lengths in centimorgans along
    one chromosome. A user who does not know the decay distance of their own
    panel measures it from the curve of r^2 against distance. A few hundred
    microsatellite loci scattered over a genome have no linkage to speak of
    and take ``"variant"``; a panel of linked variants must not use it.

    ``"variant"`` is also what the memory of the pass grows with. popnei
    keeps six numbers, 48 bytes, for each pair of populations and each
    group, so a group of each variant makes that 48 bytes for each pair and
    each variant: 1200 variants of 20 populations, which are 190 pairs, are
    10.9 MB, and a million variants are 144 MB for 3 populations and 9.1 GB
    for 20. A length in base pairs, whose groups are as many as the stretches
    of the chromosomes, does not grow with the variants. A machine that has
    not the memory is a ``ValueError`` and not a wrong number.

    `measures` says which of the seven of :class:`PopDistMeasure` to
    calculate, and ``None`` is all of them, since the pass is what costs and
    each measure is a division at the end of it. All seven are calculated,
    so ``None`` gives seven :class:`popnei.Distances` and no measure is
    refused.

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
    # The counts of the one pass go to the result and to every measure of
    # it, so that a user who takes one measure out of the result keeps the
    # counts of the pass that gave it.
    pass_stats = _pass_stats_of(counts)
    of_its_name = {
        measure: Distances(
            dist_vector=values,
            names=pops_of_the_pass,
            pass_stats=pass_stats,
            standard_errors=standard_errors,
        )
        for measure, (values, standard_errors) in zip(
            asked_for, of_each_measure, strict=True
        )
    }
    return PopDists(
        pops=pops_of_the_pass,
        num_vars=num_vars,
        pass_stats=pass_stats,
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
    # A measure that popnei has no value for is refused here, so that nobody
    # reads a vector of NaN as a distance. All seven have one, so this
    # refuses nothing today and is what a measure added to `PopDistMeasure`
    # with no formula beside it would meet.
    with_no_value = [
        measure for measure in asked_for if measure not in _MEASURES_THAT_HAVE_A_VALUE
    ]
    if with_no_value:
        raise ValueError(
            f"popnei has no value for {_named(with_no_value)}, and the "
            f"measures it has a value for are "
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
