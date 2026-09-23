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

`docs/specs/ld.md` has the calculation and the numbers the tests assert.
"""

from dataclasses import dataclass

import numpy

from popnei import _core
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
    machine asked for memory it has not. A user who wants the matrix of more
    variants and has the memory raises it; one who has not puts a filter on
    the ``Variants`` first.

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
