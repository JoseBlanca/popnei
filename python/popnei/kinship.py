"""The kinship: how related every pair of individuals is, from the genotypes.

The kinship of two individuals says how much more of their genome they
share than two individuals drawn at random from the same panel would.
:func:`calc_kinship` calculates it for every pair of the individuals of a
:class:`popnei.Variants`, and :class:`Kinship` is what it comes in.

It is the matrix of VanRaden 2008, which plink2's ``--make-rel`` computes:
an entry off the diagonal is twice the coancestry of its pair, about 0.5 for
full sibs or for a parent and a child, about 0.25 for half sibs and near 0
for two individuals with no recent ancestor in common, and an entry on the
diagonal is 1 plus the inbreeding of that individual. Entries below 0 are
ordinary and mean a pair less alike than the average pair of the panel,
because the whole matrix is measured against that average.

`docs/specs/kinship.md` has what is computed, the numbers the tests assert
and what popnei does differently from pyNei.
"""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy
import pandas

from popnei import _core
from popnei.variant import PassStats, Variants, _pass_stats_of

# How far a matrix a user built may be from its own transpose, as a share of
# its largest absolute entry. An entry of a kinship is the same for the pair
# read either way round, so a matrix that is not symmetric holds two
# different numbers for one pair. The tolerance costs a user nothing: both
# matrices plink2 wrote for the reference panels are symmetric to the bit.
_LARGEST_ASYMMETRY = 1e-9


@dataclass(frozen=True)
class Kinship:
    """The kinship of every pair of a set of individuals.

    It is what :func:`calc_kinship` gives, and a user builds one themselves
    from a matrix that was calculated elsewhere, the one plink2 or a pedigree
    gave them, so that it can be passed to the association study.

    It is pyNei's ``Kinship`` of ``pynei/gwas.py``, with these differences:
    ``samples`` and ``filter_samples`` are :attr:`individuals` and
    :meth:`filter_individuals`, the name `docs/glossary.md` gives; it carries
    :attr:`pass_stats`, the counts of the pass the calculation made, which
    pyNei keeps in its ``Variants``; and a matrix that is not square, whose
    two sides are not the same individuals, or that differs from its own
    transpose is a ``ValueError`` here, where pyNei checks nothing and the
    complaint comes later out of the linear algebra, naming a matrix and a
    row instead of the field the user filled.
    """

    matrix: pandas.DataFrame
    """The kinship of every pair, individuals x individuals, with the names
    of the individuals as index and as columns.

    The entry of a pair is in both of its cells, and the diagonal is 1 plus
    the inbreeding of each individual. The array is not copied when the
    calculation builds the frame, so that the 800 MB of 10000 individuals
    cross no layer twice."""

    num_vars: int
    """How many variants the matrix was calculated from: those that had
    variance among these individuals, which are the ones that were used.

    A variant whose called genotypes all have one dosage is in neither the
    sum of a pair nor its denominator, so this is not how many variants the
    pass gave, which is the ``num_vars`` of :attr:`pass_stats`."""

    pass_stats: PassStats | None = None
    """The counts of the pass the calculation made: how many variants it was
    given, after the steps of the ``Variants``, and how many each filter of
    it was given and kept.

    It is ``None`` in a ``Kinship`` a user built by hand, which no pass over
    a source gave."""

    def __post_init__(self) -> None:
        """The three things that make a frame of numbers a kinship, checked.

        # Raises

        ``TypeError`` when `matrix` is not a pandas frame, and ``ValueError``
        when it is not square, when its index and its columns are not the
        same individuals in the same order, or when it is further from its
        own transpose than 1e-9 of its largest absolute entry.
        """
        if not isinstance(self.matrix, pandas.DataFrame):
            raise TypeError(
                f"`matrix` is of the type `{type(self.matrix).__name__}`, and "
                f"the kinship of a set of individuals is a pandas frame with "
                f"one row and one column for each of them, indexed by their "
                f"names on both sides: give one, or calculate it from the "
                f"variants with `calc_kinship`"
            )
        rows, columns = self.matrix.shape
        if rows != columns:
            raise ValueError(
                f"the frame is {rows} by {columns}, and the kinship of a set "
                f"of individuals is a square matrix with one row and one "
                f"column for each of them"
            )
        _refuse_two_sides_that_differ(
            list(self.matrix.index), list(self.matrix.columns)
        )
        _refuse_a_matrix_that_is_not_symmetric(self.matrix)

    @property
    def individuals(self) -> tuple:
        """The names of the individuals, in the order of the rows and of the
        columns of the matrix.

        They are read from the matrix itself, so the names are in one place
        and cannot disagree with themselves, as pyNei's ``samples`` is."""
        return tuple(self.matrix.index)

    def filter_individuals(self, individuals: Sequence[str]) -> Kinship:
        """The rows and the columns of `individuals`, in the order given.

        It takes entries out of the matrix and calculates nothing, so
        :attr:`num_vars` and :attr:`pass_stats` are those of the kinship it
        came from and every entry is unchanged. That is not the kinship of
        those individuals alone: the entries of one are measured against the
        average pair of the panel they were calculated over, and
        ``calc_kinship(variants, individuals=...)`` gives the other, which
        has the frequencies, the means and the denominators of those
        individuals.

        A name that is not an individual of the matrix is a ``ValueError``
        that names it, as pyNei's ``filter_samples`` raises.
        """
        named = list(individuals)
        missing = [name for name in named if name not in self.matrix.index]
        if missing:
            raise ValueError(
                f"{_the_names(missing)} of the kinship, which is of "
                f"{len(self.matrix.index)} individuals: a name of "
                f"`individuals` is one of them"
            )
        return Kinship(
            matrix=self.matrix.loc[named, named],
            num_vars=self.num_vars,
            pass_stats=self.pass_stats,
        )


def calc_kinship(
    variants: Variants,
    individuals: Sequence[str] | None = None,
    transform_to_biallelic: bool = _core.DEFAULT_TRANSFORM_TO_BIALLELIC,
) -> Kinship:
    """The kinship of every pair of the individuals of `variants`.

    Each variant becomes one number per individual, its dosage: how many
    alleles of the genotype are not the major allele of the variant, which is
    the most frequent among its called alleles and the lowest numbered of two
    that are equally frequent. The dosages of a variant are centered and
    divided by the standard deviation its allele frequency gives it under
    Hardy Weinberg, ``sqrt(ploidy * p * (1 - p))`` with `p` the mean dosage
    over the ploidy, and the entry of a pair is the sum over the variants of
    the two standardized dosages multiplied, divided by how many of those
    variants have a called genotype in both individuals. That divisor is what
    makes an entry twice a coancestry, and it is what plink2's ``--make-rel``
    and GCTA use; the principal components of
    :func:`popnei.do_pca_from_variants` divide by another number, the
    standard deviation of the dosages themselves, and the two agree only when
    the genotypes are in Hardy Weinberg proportions.

    A genotype with any allele missing takes the mean of the dosages of its
    variant, so that after centering it pulls its pairs nowhere, and it does
    not count in the denominator of any pair it is in: every pair has its own
    number of variants, which is why an individual with much missing data
    gets entries of the same kind as the others.

    A variant whose called genotypes all have one dosage has no variance and
    is left out of both the sum and the denominators, a variant with one
    allele, one with no called genotype, and one where every individual is
    heterozygous, whose major allele frequency is 0.5 and which no filter by
    frequency catches. The result's ``num_vars`` counts those that were used,
    and the ``num_vars`` of its ``pass_stats`` how many the steps of the
    `Variants` let through.

    The call makes one pass over the source of `variants`, through the steps
    that are on it, so the kinship is over the variants its filters kept, and
    the ``Variants`` is as it was afterwards. A user normally prunes by
    linkage disequilibrium first, with ``filter_by_ld``, so that a block of
    correlated variants does not count many times.

    `individuals` are the individuals the matrix is of, in the order given,
    and ``None`` is all of them in the order of the pass. Every frequency,
    mean and denominator is of those individuals, so the kinship of some of
    them is not the rows and columns of the kinship of the whole panel: on
    the reference panel of 200 individuals and 1200 variants, the kinship of
    40 of them differs from the same 40 rows and columns of the kinship of
    all 200 by up to 0.129, and it uses 1195 variants, the other 5 having no
    variance among those 40. A name that is not an individual of `variants`,
    a name given twice and an empty sequence are a ``ValueError``.

    `transform_to_biallelic` makes every allele that is not the major one
    count the same, which is what a variant of more than two different
    alleles among its called genotypes needs: the dosage of a genotype has a
    meaning for two alleles, and without this such a variant is a
    ``ValueError`` that says which one it is, as it is in
    :func:`popnei.do_pca_from_variants`. The alleles are those the genotypes
    hold and not those the source lists.

    A dataset in which no variant has variance among the individuals asked
    for is a ``ValueError``, and so is a pass that gives no variant at all,
    whose message says whether the source held none or the steps kept none.
    A pair of individuals with no variant called in both of them is one too,
    naming the two and how many variants each of them has called, where pyNei
    divides by 0 and leaves a NaN in the matrix; the user leaves one of the
    two out. So is a dataset of a size the calculation cannot count in, a
    ploidy above 254 or more than 46340 individuals.

    It is pyNei's ``calc_kinship`` of ``pynei/gwas.py``, with these
    differences: `samples` is `individuals`; there is no `num_threads`,
    because no calculation of popnei has one and the threads are those of
    rayon's pool, one for each core unless ``RAYON_NUM_THREADS`` is set
    before the first call; `transform_to_biallelic` is new, and pyNei reads a
    variant of more than two alleles with every allele that is not the major
    one counting the same, silently; a pair with no variant called in both
    raises instead of leaving a NaN; and the result carries the counts of the
    pass, which pyNei keeps in its ``Variants``.
    """
    if not isinstance(variants, Variants):
        # What a user gives instead is usually the path of the VCF, and what
        # that gave was the `AttributeError` of an object with no source
        # inside it.
        raise TypeError(
            f"`variants` is {variants!r}, a {type(variants).__name__}, and "
            f"`calc_kinship` reads the variants of a source: give it what "
            f"`open_vcf` or `open_vars` gives, "
            f"calc_kinship(open_vcf(vcf_path))"
        )
    # The names are given to the core as they are: it turns them into the
    # places of those individuals among the ones the pass gives, and refuses
    # a name that is not one of them, as the filter of individuals does.
    matrix, names, num_vars, counts = _core.calc_kinship(
        variants._source,
        None if individuals is None else list(individuals),
        transform_to_biallelic,
        variants._steps,
    )
    # `copy=False`: without it pandas allocates a second array of the same
    # size and copies into it, which at 10000 individuals is 800 MB of peak
    # memory against none. What the keyword needs is that nothing else holds
    # the array, and nothing does: the core crate allocated it for this call,
    # `_core.calc_kinship` is the only reference to it, and the frame is what
    # outlives the call.
    return Kinship(
        matrix=pandas.DataFrame(matrix, index=names, columns=names, copy=False),
        num_vars=num_vars,
        pass_stats=_pass_stats_of(counts),
    )


def _refuse_two_sides_that_differ(index: list, columns: list) -> None:
    """The first name of `index` that is not the name of the column at the
    same place, refused.

    # Raises

    ``ValueError`` when the two sides of the frame are not the same
    individuals in the same order, which sorting the rows of a matrix and not
    its columns gives: its cells are then no longer the kinship of the pair
    of their row and their column.
    """
    for place, (row, column) in enumerate(zip(index, columns, strict=True)):
        if row != column:
            raise ValueError(
                f"the name of the row {place} is {row!r} and the name of its "
                f"column is {column!r}, and the two sides of a kinship are "
                f"the same individuals in the same order: the cell of the row "
                f"{row!r} and the column {column!r} is the kinship of that "
                f"pair and not of a pair of one of them. Sort both sides the "
                f"same way, or take the individuals you want with "
                f"`filter_individuals`"
            )


def _refuse_a_matrix_that_is_not_symmetric(matrix: pandas.DataFrame) -> None:
    """A matrix further from its own transpose than 1e-9 of its largest
    absolute entry, refused.

    The pair that is furthest from its transpose is the one named, with the
    two numbers the matrix holds for it, so that a user who built the frame
    from half a matrix sees which cell was not filled.

    # Raises

    ``ValueError`` when the two cells of a pair differ by more than that.
    """
    values = numpy.asarray(matrix.to_numpy(dtype=numpy.float64))
    # `initial` gives the largest of no entry at all, which a matrix of no
    # individual has: numpy has no largest of an empty array without it.
    largest = numpy.abs(values).max(initial=0.0)
    gaps = numpy.abs(values - values.T)
    widest = gaps.max(initial=0.0)
    if widest <= largest * _LARGEST_ASYMMETRY:
        return
    row, column = numpy.unravel_index(numpy.argmax(gaps), gaps.shape)
    names = list(matrix.index)
    raise ValueError(
        f"the cell of the row {names[row]!r} and the column {names[column]!r} "
        f"holds {values[row, column]!r} and the cell of the row "
        f"{names[column]!r} and the column {names[row]!r} holds "
        f"{values[column, row]!r}, and the kinship of a pair is one number, "
        f"which both of its cells hold: the two differ by {widest!r}, where "
        f"{_LARGEST_ASYMMETRY} of the largest absolute entry of the matrix, "
        f"{largest!r}, is what a kinship is taken to be symmetric within"
    )


def _the_names(missing: list) -> str:
    """The individuals of a call that are not in the matrix, as the start of
    the message that refuses them, in the singular when there is one."""
    named = ", ".join(repr(name) for name in missing)
    if len(missing) == 1:
        return f"{named} is not an individual"
    return f"{named} are not individuals"
