"""The r² between variants, and the matrix of it.

Two variants are in linkage disequilibrium when the genotype of one tells
something about the genotype of the other, which happens when they sit close
enough on a chromosome that few recombinations have separated them. r² is
what measures it: each variant becomes one number per individual, its
dosage, how many alleles of the genotype are not the major allele of the
variant, and r² is the square of the correlation between the dosages of the
two variants, 1 when the dosage of an individual at one of them fixes its
dosage at the other and 0 when knowing one says nothing about the other. It
is the estimate named after Rogers and Huff, which takes genotypes whose
phase is unknown, so nothing has to be phased for it.

:func:`calc_rogers_huff_r2_matrix` gives it for every pair of the variants
of a :class:`popnei.Variants`, and :class:`R2Matrix` is what it comes in.

Two variants that sit close together on a chromosome have had fewer
recombinations between them than two that sit far apart, so r² falls off as
the distance grows, and how fast it falls off is a property of the
population: one that went through few individuals keeps linkage
disequilibrium over longer stretches.
:func:`calc_ld_and_dist_per_pop` gives that fall-off for each population, as
the mean r² of the pairs of each bin of distance and as a curve fitted to its
pairs, :class:`LdDecay`, which carries the distance at which r² has fallen to
half. :class:`LdAndDistPerPop` is what the two come in.

`docs/specs/ld.md` has the calculation and the numbers the tests assert.
"""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy
import pandas

from popnei import _core
from popnei.stats import _the_pops
from popnei.variant import PassStats, Variants, _pass_stats_of


@dataclass(frozen=True, eq=False, repr=False)
class R2Matrix:
    """The r² of every pair of the variants of one pass.

    ``matrix == other`` is true for the same object and false for any other,
    as it is for a :class:`popnei.Distances`: two results are not compared
    value by value, because a matrix of r² is neither equal nor unequal to
    another, it is equal cell by cell.
    ``numpy.array_equal(matrix.r2, other.r2, equal_nan=True)`` is how two
    matrices are compared, NaN counting as equal to NaN.

    It is pyNei's ``R2Matrix`` of ``pynei/ld.py``, with these differences:
    :attr:`r2` holds r² where pyNei's holds r, which is the correlation
    itself and carries a sign; there is no ``dists_in_bp``, the square matrix
    of the distance between every pair, since :attr:`chroms` and :attr:`poss`
    hold the same thing in 2n numbers instead of n², 80 KB against 200 MB at
    5000 variants, and the distance of a pair is the difference of two
    positions; and it carries :attr:`pass_stats`, the counts of the pass,
    which pyNei keeps in its ``Variants``.
    """

    r2: numpy.ndarray
    """The r² of every pair, a read only square float64 array with one row
    and one column for each variant the pass gave, in the order it gave
    them.

    The two cells of a pair hold the same value. A pair that has no r² is
    NaN: the individuals called at both of its variants hold one dosage at
    one of them, or there are no such individuals. The diagonal is 1 for a
    variant whose called genotypes hold two dosages at least and NaN for one
    that has no variance, which has no r² against any variant, itself
    included.
    """

    chroms: tuple[str, ...]
    """The name of the chromosome of each variant, one for each row of
    :attr:`r2`."""

    poss: numpy.ndarray
    """The position of each variant, a read only array of whole numbers, 1
    based as in a VCF.

    The distance of a pair is the difference of two of them, and a pair whose
    variants are on two chromosomes has none.
    """

    pass_stats: PassStats
    """The counts of the pass the calculation made: how many variants it
    took, which is how many rows the matrix has, and how many each filter of
    the ``Variants`` was given and kept."""

    def __repr__(self) -> str:
        """How many variants the matrix is of, and not the matrix.

        The one a dataclass writes prints the name of the chromosome of every
        variant, which for 5000 variants is tens of kilobytes in a session or
        in a traceback.
        """
        num_vars = len(self.chroms)
        return (
            f"<R2Matrix of {num_vars} {'variant' if num_vars == 1 else 'variants'}, "
            f"{len(set(self.chroms))} of chromosomes, with the counts of its pass>"
        )


def calc_rogers_huff_r2_matrix(
    variants: Variants, max_num_vars: int = _core.DEFAULT_MAX_NUM_VARS
) -> R2Matrix:
    """The r² of every pair of the variants of `variants`.

    The dosage of an individual at a variant is how many alleles of its
    genotype are not the major allele of that variant, 0, 1 or 2 in a
    diploid, and r² is the square of the correlation between the dosages of
    the two variants of a pair. Every allele that is not the major one counts
    the same, so a variant of more than two alleles is read as two.

    An individual whose genotype is missing at either variant of a pair, a
    genotype with any allele not called, ``0/.`` among them, is left out of
    that pair, so every pair has its own number of individuals. A pair whose
    individuals hold one dosage at one of its variants has no r², and its
    cell is NaN; so a variant with no variance has NaN in its whole row, its
    whole column and its diagonal cell.

    The call makes one pass over the source of `variants`, through the steps
    that are on it, so the matrix is of the variants its filters kept, and
    the ``Variants`` is as it was afterwards.

    `max_num_vars` is how many variants the calculation takes before it
    refuses. The matrix holds one r² for each pair of them, 8 bytes, so it
    grows with the square of the variants: 200 MB at the default of 5000 and
    80 GB at 100000. It is the one result of popnei that grows with the
    square of its input, so a pass of more variants is a ``ValueError``,
    whose message says how many variants the pass had given, what the cap
    was and how much memory their matrix would have needed, and not a
    machine asked for memory it has not. That message names the file the
    pass was reading, since the cap and the variants that file holds decide
    together whether it is passed. A user who wants the matrix of more
    variants and has the memory raises it; one who has not puts a filter on
    the ``Variants`` first.

    A cap raised above what the machine then gives is a ``ValueError`` too,
    one that says how many values could not be held. The matrix is held
    once: the array a caller gets is the memory the calculation filled, and
    nothing copies it, so 200 MB at 5000 variants is what the call asks the
    machine for and what it keeps afterwards.

    What it gives is an :class:`R2Matrix` with the matrix, the chromosome and
    the position of each variant and the counts of the pass in its
    ``pass_stats``.

    A pass that gives no variant is a ``ValueError``, whether the source has
    none or the steps kept none: the message says which of the two, and what
    each filter was given and kept when there are steps.

    It is pyNei's ``calc_rogers_huff_r2_matrix``, with these differences: it
    gives r² where pyNei gives r, so that the name and the value agree, and
    what is lost is the sign, which says whether the major alleles of the two
    variants go together or apart; a missing genotype takes its individual
    out of that pair, where pyNei leaves it in with a dosage of -1, and the
    rule here is plink2's, which is the reference program of this
    calculation; `max_num_vars` is new, pyNei holding every chunk of the
    dataset in memory and building the matrix of all of them with nothing to
    stop it; there is no `max_dist`, which leaves NaN in the cells of the
    pairs further apart than it; there is no `check_no_mafs_above`, which
    raises for the whole call when a variant has a major allele frequency
    above 0.95, and which in popnei is what ``variants.filter_by_maf(0.95)``
    does, a step that takes those variants out instead of refusing the
    dataset; and the result carries the counts of the pass.
    """
    if not isinstance(variants, Variants):
        # The path of the VCF is the mistake that is easiest to make, and
        # what it gave was the `AttributeError` of an object with no source
        # inside it.
        raise TypeError(
            f"`variants` is {variants!r}, of the type "
            f"`{type(variants).__name__}`, and the r² of every pair is "
            f"calculated over the variants of a source: give it what "
            f"`open_vcf` or `open_vars` gives, "
            f"calc_rogers_huff_r2_matrix(open_vcf(vcf_path))"
        )
    r2, chroms, poss, counts = _core.calc_rogers_huff_r2_matrix(
        variants._source, max_num_vars, variants._steps
    )
    return R2Matrix(r2=r2, chroms=chroms, poss=poss, pass_stats=_pass_stats_of(counts))


@dataclass(frozen=True)
class LdDecay:
    """The curve of r² against distance fitted to the pairs of one
    population, and the distance at which it has fallen to half.

    The curve is the r² that two variants of a population are expected to
    have at a given recombination between them, under drift and
    recombination: Hill and Weir (1988), with the correction of Weir and Hill
    (1986) for r² being measured on a sample of individuals and not on the
    whole population, which is what holds the curve up at long distances. It
    has one number to fit, :attr:`rho_per_bp`, and the number of individuals
    of the population enters it as it is given.

    The curve is fitted to every pair the pass counted, each at its own
    distance, so the bins do not move it and two runs over one dataset give
    the same three numbers.

    A population whose pairs no curve was fitted to has NaN in all three, and
    that is no error. It happens to a population whose pairs fall at fewer
    than two distances, one distance saying nothing about a fall-off, which
    includes a population with no pair at all and one left with two variants;
    and to a population whose best fit falls at an end of the range of the ρ
    per base pair that is searched, 1e-12 to 100, a curve flat across
    `max_dist` or fallen before the second base pair being no fall-off that
    these pairs pin down.

    `docs/specs/ld.md` has the curve, what is fitted to what and the numbers
    the tests assert.
    """

    rho_per_bp: float
    """The fitted 4Nr: four times the effective size of the population times
    the recombination per base pair, which is by how much ρ, the scaled
    recombination of the curve, grows with each base pair between the two
    variants of a pair.

    The effective size and the recombination rate enter the curve only as
    that product, and one pass over one dataset does not separate them, so
    neither is given on its own.
    """

    r2_at_zero: float
    """The fitted curve at a distance of 0, which is its own ceiling: two
    variants that never recombine still do not reach an r² of 1, because
    their allele frequencies drift apart.

    How many individuals the population has fixes it on its own,
    0.46198347107438015 at 100 of them, and no pair of the dataset moves it.
    """

    half_dist: float
    """The distance in base pairs at which the fitted curve has fallen to
    half of :attr:`r2_at_zero`.

    What is halved is the curve at a distance of 0 and not the mean r² of the
    shortest bin, so `num_bins` does not move this distance; `min_dist` and
    `max_dist` do, through :attr:`rho_per_bp`, since they choose which pairs
    the curve is fitted to.

    It is NaN when the other two are, and NaN on its own for a curve that
    never falls to half of its value at 0, which happens below three
    individuals: what the curve falls towards as ρ grows is 1 over the
    individuals, which is above that half at one and at two of them. No pass
    of :func:`calc_ld_and_dist_per_pop` reaches such a population, since one
    individual has no variant with variance and so no pair, and two give
    every pair an r² of 1, which the three NaN above cover.
    """


@dataclass(frozen=True, eq=False, repr=False)
class LdAndDistPerPop:
    """How the r² of a pair of variants falls off with the distance between
    them, for each population: in bins of distance, and as a curve fitted to
    its pairs with the distance at which that curve has fallen to half.

    ``result == other`` is true for the same object and false for any other,
    as it is for an :class:`popnei.R2Matrix`: the frames of two results are
    compared frame by frame, with :func:`pandas.testing.assert_frame_equal`
    or with ``frame.equals(other)``, and a dataclass that compared them with
    ``==`` would raise instead of answering.
    """

    per_pop: dict[str, pandas.DataFrame]
    """The bins of each population, under its name and in the order of the
    `pops` dict that was given.

    The frame of a population has one row per bin, in the order of the
    distances, indexed by ``smallest_dist``, the smallest distance of the
    bin, and with the columns ``largest_dist``, the largest distance of the
    bin, both ends included; ``num_pairs``, how many pairs of variants fell
    in it; ``mean_r2``, the mean of their r²; and ``sd_r2``, the standard
    deviation of their r², with the pairs of the bin as the divisor. A bin
    with no pair has 0 pairs and NaN in the last two columns, and a bin of
    one pair has that pair's r² and a standard deviation of 0.
    """

    num_vars_per_pop: dict[str, int]
    """How many variants each population kept at its major allele frequency,
    under its name.

    It is worked out over the individuals of that population alone, so two
    populations of one pass count different variants, and a variant that no
    individual of a population has called is out of it.
    """

    decay_per_pop: dict[str, LdDecay]
    """The curve fitted to the pairs of each population, under its name and
    in the order of the `pops` dict that was given, as :attr:`per_pop` and
    :attr:`num_vars_per_pop` are.

    The :class:`LdDecay` of a population holds the fitted ρ per base pair,
    the curve at a distance of 0 and the distance at which it has fallen to
    half of that, the three of them NaN for a population whose pairs no curve
    was fitted to.
    """

    pass_stats: PassStats
    """The counts of the pass the calculation made: how many variants it
    took, before the major allele frequency of any population, and how many
    each filter of the ``Variants`` was given and kept."""

    def __repr__(self) -> str:
        """How many populations and how many bins the result holds, and not
        the frames.

        The one a dataclass writes prints every row of every frame, which
        for the 50 bins of the default is hundreds of lines in a session or
        in a traceback.
        """
        num_pops = len(self.per_pop)
        num_bins = max((len(frame) for frame in self.per_pop.values()), default=0)
        return (
            f"<LdAndDistPerPop of {num_pops} "
            f"{'population' if num_pops == 1 else 'populations'} in {num_bins} "
            f"{'bin' if num_bins == 1 else 'bins'} of distance, "
            f"with the counts of its pass>"
        )


def calc_ld_and_dist_per_pop(
    variants: Variants,
    pops: dict[str, Sequence[str]] | None = None,
    min_dist: int = _core.DEFAULT_MIN_DIST,
    max_dist: int = _core.DEFAULT_MAX_DIST,
    num_bins: int = _core.DEFAULT_NUM_DIST_BINS,
    max_allowed_maf: float = _core.DEFAULT_MAX_ALLOWED_MAF,
) -> LdAndDistPerPop:
    """How the r² of a pair of variants falls off as the two move apart
    along a chromosome, for each population of `pops`.

    A pair of variants is counted in a population when both of its variants
    passed the major allele frequency of that population, when the two are on
    one chromosome, and when their distance, the difference of their
    positions, is from `min_dist` to `max_dist`, both included. The pairs are
    put into `num_bins` bins of equal width across that range, and each bin
    gets how many pairs it holds, the mean of their r² and its standard
    deviation. A pair that has no r², one whose individuals called at both
    variants hold one dosage at one of them, is in no bin, and so is a pair
    whose two variants are on two chromosomes.

    The curve of a population is read on its own because the recombination it
    has had and the number of individuals it has been through shape it: a
    population that went through few individuals keeps linkage disequilibrium
    over longer stretches.

    The call makes one pass over the source of `variants`, through the steps
    that are on it, which serves every population, and the ``Variants`` is as
    it was afterwards.

    `pops` is a dict of population name to the names of its individuals,
    which are looked up among :attr:`popnei.Variants.individuals`, the ones
    the pass gives. With no `pops` there is one population, ``pop``, of every
    individual. A name that is not an individual of the pass, a name twice in
    one population, a population that names no individual and a `pops` with
    no population are each a ``ValueError``; an individual in two populations
    is read by both, and one in none by neither.

    `min_dist` and `max_dist` are the distances in base pairs a pair is
    counted at, 1 and 1000000 by default, both ends included. A `min_dist` of
    1 leaves out only the pairs of two variants at one position, a SNP and an
    indel at the same base, whose distance is 0. `max_dist` is also how far
    back the pass holds the variants it has read, so it is what the memory of
    the call grows with. A `min_dist` above `max_dist`, and either of them
    below 0, is a ``ValueError`` that names the argument and the value.

    `num_bins` is how many bins of equal width the distances are cut into, 50
    by default: the width is (`max_dist` − `min_dist` + 1) / `num_bins`, and
    a pair at the distance d falls in the bin that
    floor((d − `min_dist`) / width) gives, the last bin taking anything the
    rounding would put past it. A `num_bins` of 0 is a ``ValueError``.

    `max_allowed_maf` is the largest major allele frequency a variant has in
    a population and is still counted there, 0.95 by default, both ends
    included. The major allele frequency is the count of the commonest allele
    over the called alleles, worked out over the individuals of that
    population alone, so two populations of one pass count different
    variants; a variant that no individual of a population called has none
    there and is left out of it, and can still be counted in another
    population. Those variants are left out because the r² of a variant that
    hardly varies rests on the one or two individuals that carry the rare
    allele, and keeping them raises the curve everywhere. Anything that is
    not a number from 0 to 1 is a ``ValueError``.

    Beside the bins each population gets a curve fitted to its pairs, in its
    `decay_per_pop`: an :class:`LdDecay` with the fitted ρ per base pair, the
    curve at a distance of 0 and the distance at which it has fallen to half
    of that. The curve is fitted to every pair and at the distance of each
    pair, so `num_bins` does not move it, and a population whose pairs no
    curve was fitted to has NaN in all three, which :class:`LdDecay` says
    when and which is no error.

    What it gives is an :class:`LdAndDistPerPop` with the bins of each
    population in its `per_pop`, how many variants each of them kept in its
    `num_vars_per_pop`, the curve of each of them in its `decay_per_pop`, and
    the counts of the pass in its `pass_stats`.

    A population in which every variant was left out, and a dataset whose
    variants are all further apart than `max_dist` or each on a chromosome of
    their own, give every bin empty, which is no error. A pass that gives no
    variant is a ``ValueError``, whether the source has none or the steps
    kept none.

    It mirrors pyNei's ``calc_ld_and_dist_per_pop``, with these differences:
    it gives the bins over every pair, where pyNei gives a sample of at most
    ``max_num_measures_to_keep`` pairs drawn with no seed, so that two runs
    over one dataset give the same numbers here and different points there;
    it fits the curve and gives the half distance, where pyNei hands its
    sample of pairs over and leaves the fitting to the user;
    it gives r² where pyNei gives r, and a missing genotype takes its
    individual out of that pair where pyNei leaves it in with a dosage of -1;
    the distances come as `min_dist` and then `max_dist`, where pyNei's
    signature has `max_dist` first, so a call written for pyNei that gives
    them by position asks for something else here and nothing says so;
    `min_dist` counts the pair at that distance, where pyNei keeps the pairs
    strictly beyond it and so loses the pairs 1 base pair apart at its own
    default of 1; `max_dist` has a default, plink2's own
    ``--ld-window-kb 1000``, where pyNei's is ``None`` and keeps every pair of
    every chromosome; `num_bins` and the standard deviation are new; there is
    no ``method``, pyNei's two counting each unordered pair once and twice;
    and it makes one pass that serves every population, where pyNei makes one
    pass per population and so reads the source once per population.
    """
    if not isinstance(variants, Variants):
        # The path of the VCF is the mistake that is easiest to make, and
        # what it gave was the `AttributeError` of an object with no source
        # inside it.
        raise TypeError(
            f"`variants` is {variants!r}, of the type "
            f"`{type(variants).__name__}`, and the fall-off of r² with "
            f"distance is calculated over the variants of a source: give it "
            f"what `open_vcf` or `open_vars` gives, "
            f"calc_ld_and_dist_per_pop(open_vcf(vcf_path))"
        )
    named = _the_pops(pops)
    pop_names, of_each_pop, counts = _core.calc_ld_and_dist_per_pop(
        variants._source,
        variants._steps,
        named,
        min_dist,
        max_dist,
        num_bins,
        max_allowed_maf,
    )
    per_pop = {}
    num_vars_per_pop = {}
    decay_per_pop = {}
    for name, bins in zip(pop_names, of_each_pop, strict=True):
        smallest_dist, largest_dist, num_pairs, mean_r2, sd_r2, num_vars, curve = bins
        per_pop[name] = pandas.DataFrame(
            {
                "largest_dist": largest_dist,
                "num_pairs": num_pairs,
                "mean_r2": mean_r2,
                "sd_r2": sd_r2,
            },
            index=pandas.Index(smallest_dist, name="smallest_dist"),
        )
        num_vars_per_pop[name] = num_vars
        rho_per_bp, r2_at_zero, half_dist = curve
        decay_per_pop[name] = LdDecay(
            rho_per_bp=rho_per_bp, r2_at_zero=r2_at_zero, half_dist=half_dist
        )
    return LdAndDistPerPop(
        per_pop=per_pop,
        num_vars_per_pop=num_vars_per_pop,
        decay_per_pop=decay_per_pop,
        pass_stats=_pass_stats_of(counts),
    )
