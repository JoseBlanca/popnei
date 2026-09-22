"""The distances between individuals, and the result they come in.

The Kosman distance of two individuals is the share of their genotypes that
differ, from 0 for two that hold the same genotype at every variant to 1 for
two that have no allele in common at any.
:func:`popnei.calc_pairwise_kosman_dists` gives it for every pair of the
individuals of a source, and :class:`Distances` is what it and every later
distance calculation of popnei give it in.

`docs/specs/dists.md` has the distance, what it is worked out from and the
numbers the tests assert.
"""

import math
import operator
from collections.abc import Iterable
from dataclasses import dataclass

import numpy
import pandas

from popnei import _core
from popnei.variant import PassStats, Variants, _pass_stats_of

# The largest `min_num_snps` popnei takes, which is the largest number of
# variants a pair can have been called at: popnei keeps how many variants a
# pair shares in a 32 bit whole number and refuses a pass whose count would
# go above it, so no pair reaches a threshold above this and a user who
# wrote one is told instead of getting a vector of NaN.
_MOST_VARIANTS_A_PAIR_CAN_HAVE = 4294967295


@dataclass(frozen=True, eq=False, repr=False)
class Distances:
    """The distance between every pair of individuals.

    It is what every distance calculation of popnei gives, and a user builds
    one themselves from distances that were calculated elsewhere, with
    ``Distances(dist_vector, names=...)`` or with
    :meth:`Distances.from_square_dists`.

    ``dists == other`` is true for the same object and false for any other,
    as it is for a :class:`popnei.Block`: two results are not compared value
    by value, because an array of distances is neither equal nor unequal to
    another, it is equal element by element.
    ``numpy.array_equal(dists.dist_vector, other.dist_vector,
    equal_nan=True)`` is how two vectors are compared, NaN counting as equal
    to NaN.

    It is pyNei's ``Distances`` of ``pynei/dists.py``, with these
    differences: `names` is a tuple, where pyNei has a numpy array; a vector
    whose length is the pairs of no number of individuals is a
    ``ValueError``, where pyNei takes it and leaves its last values in no
    cell; :attr:`triang_list_of_lists` is the lower triangle of
    :attr:`square_dists`, where pyNei's is not from four individuals on; and
    it carries :attr:`pass_stats`, the counts of the pass the calculation
    made, which pyNei keeps in its ``Variants``.
    """

    dist_vector: numpy.ndarray
    """The distance of every pair, a read only float64 array.

    The pairs are in the order (0, 1), (0, 2), ..., (0, N-1), (1, 2), ...,
    the upper triangle of the square matrix row by row, which is the
    condensed form of ``scipy.spatial.distance``. A pair that has no
    distance is NaN: the two individuals were called together at no variant,
    or at fewer variants than the calculation was asked for.

    The array is not copied when a ``Distances`` is built from one, so that
    the 400 MB of 10000 individuals cross no layer twice; a user who builds
    one from an array of their own and then writes into that array changes
    what the result holds, where pyNei copies it.
    """

    names: Iterable[str] | Iterable[int] | None = None
    """The names of the individuals, in the order of the pairs.

    Any sequence of names is taken, and what a built ``Distances`` holds is
    always a tuple of them. ``None`` asks for the names 0 to N-1, as in
    pyNei, and the calculation gives the names the source has for its
    individuals.
    """

    pass_stats: PassStats | None = None
    """The counts of the pass the calculation made: how many variants it
    took, after the steps, and how many each filter of the ``Variants`` was
    given and kept.

    It is ``None`` in a ``Distances`` built from distances that were
    calculated elsewhere, which no pass over a source gave.
    """

    def __post_init__(self) -> None:
        """The vector as a read only array of float64, and the names of the
        individuals its length says there are."""
        try:
            vector = numpy.asarray(self.dist_vector, dtype=numpy.float64)
        except (TypeError, ValueError) as problem:
            raise ValueError(
                f"`dist_vector` is of the type "
                f"`{type(self.dist_vector).__name__}` and holds what is no "
                f"number, and the distance of a pair is one: {problem}"
            ) from None
        if vector.ndim != 1:
            raise ValueError(
                f"`dist_vector` has {vector.ndim} dimensions, and the distances "
                f"of the pairs are one row of numbers: give the upper triangle "
                f"of the square matrix row by row, or build the result from "
                f"the matrix itself with `Distances.from_square_dists`"
            )
        if vector.flags.writeable:
            # A view of the array and not the array itself, so that the one a
            # user gave is theirs to write into afterwards, and a copy of
            # neither: the vector of 10000 individuals is 400 MB.
            vector = vector.view()
            vector.flags.writeable = False
        num_pairs = vector.shape[0]
        num_individuals = _num_individuals_of(num_pairs)
        if self.names is None:
            names = tuple(range(num_individuals))
        else:
            try:
                names = tuple(self.names)
            except TypeError:
                raise TypeError(
                    f"`names` is {self.names!r}, of the type "
                    f"`{type(self.names).__name__}`, and the names of the "
                    f"individuals are a sequence of one name for each of "
                    f"them: give a tuple or a list, or `None` for the names "
                    f"0 to N-1"
                ) from None
            # An empty vector is what one individual gives, and what no
            # individual would give: the two cannot be told apart from it,
            # so a user who names either is taken at their word.
            named = (0, 1) if num_pairs == 0 else (num_individuals,)
            if len(names) not in named:
                raise ValueError(
                    f"{len(names)} names were given and {num_pairs} distances "
                    f"are the pairs of {_of_how_many_individuals(num_pairs)}: "
                    f"a name is needed for each of them"
                )
        object.__setattr__(self, "dist_vector", vector)
        object.__setattr__(self, "names", names)

    def __repr__(self) -> str:
        """How many individuals the distances are of, and not the distances.

        The one a dataclass writes prints every name, which for 10000
        individuals is hundreds of kilobytes in a session or in a traceback.
        """
        pairs = self.dist_vector.shape[0]
        return (
            f"<Distances of {_individuals(len(self.names))}, {pairs} "
            f"{'pair' if pairs == 1 else 'pairs'}"
            f"{'' if self.pass_stats is None else ', with the counts of its pass'}>"
        )

    @classmethod
    def from_square_dists(cls, dists: pandas.DataFrame) -> Distances:
        """A ``Distances`` of the distances of a square frame.

        The frame is the one :attr:`square_dists` gives: the distance of
        every pair in both of its cells, indexed by the names of the
        individuals on both sides, which are the names of the result. Its
        diagonal is not read, and neither is its lower triangle.

        A frame that is not square is a ``ValueError``, and so is one whose
        index and columns are not the same individuals in the same order,
        which sorting the rows of a square matrix and not its columns
        gives: its cells are then no longer the distance of the pair of
        their row and their column, and the upper triangle of it is the
        distance of another pair in every place.

        The result has no :attr:`pass_stats`, since no pass over a source
        gave it.
        """
        if not isinstance(dists, pandas.DataFrame):
            raise TypeError(
                f"`dists` is of the type `{type(dists).__name__}`, and "
                f"`from_square_dists` builds a result from the pandas frame "
                f"that `square_dists` gives, indexed by the names of the "
                f"individuals on both sides: give one, or build the result "
                f"from the distances of the pairs with `Distances(...)`"
            )
        rows, columns = dists.shape
        if rows != columns:
            raise ValueError(
                f"the frame is {rows} by {columns}, and the distances of a set "
                f"of individuals are a square matrix with one row and one "
                f"column for each of them"
            )
        _refuse_two_sides_that_differ(list(dists.index), list(dists.columns))
        values = numpy.asarray(dists.values, dtype=numpy.float64)
        first, second = numpy.triu_indices(rows, k=1)
        return cls(dist_vector=values[first, second], names=tuple(dists.index))

    @property
    def square_dists(self) -> pandas.DataFrame:
        """The N x N frame of the distances, indexed by name on both sides.

        The distance of a pair is in both of its cells, the diagonal is 0,
        also for an individual that has no called genotype, and both cells
        of a pair with no distance are NaN.

        It is the form a tree or a principal coordinate analysis takes. It
        is built at every read, out of the vector, and for 10000 individuals
        it is 800 MB.
        """
        num_individuals = len(self.names)
        square = numpy.zeros((num_individuals, num_individuals), dtype=numpy.float64)
        first, second = numpy.triu_indices(num_individuals, k=1)
        square[first, second] = self.dist_vector
        square[second, first] = self.dist_vector
        return pandas.DataFrame(square, index=self.names, columns=self.names)

    @property
    def triang_list_of_lists(self) -> list[list[float]]:
        """The lower triangle of :attr:`square_dists` with its diagonal, row
        r holding r + 1 values.

        It is the form Biopython's ``DistanceMatrix`` takes. pyNei's cuts
        the vector into runs of 1, 2, 3... values, which are not the rows of
        the lower triangle from four individuals on.
        """
        square = self.square_dists.values
        return [[*square[row, :row].tolist(), 0.0] for row in range(len(self.names))]


def calc_pairwise_kosman_dists(
    variants: Variants, min_num_snps: int | None = None
) -> Distances:
    """The Kosman distance between every pair of individuals of `variants`.

    At one variant the two called genotypes are laid side by side, each
    allele of one paired with an allele of the other in the pairing that
    leaves the fewest pairs of different alleles, and d is that number of
    pairs over the ploidy. For diploids d is 0 for two genotypes that hold
    the same alleles, 0.5 for two that share one, and 1 for two that share
    none. The distance of the pair is the mean of d over the variants at
    which both genotypes are called, so every pair has its own number of
    variants, and a genotype that is half called, ``0/.``, is a missing
    genotype.

    It is the distance of Kosman and Leonard (2005, Molecular Ecology 14:
    415, DOI 10.1111/j.1365-294X.2005.02416.x) for codominant markers, at
    any ploidy, and it is what a user builds a tree or a principal
    coordinate analysis of individuals from.

    The call makes one pass over the source of `variants`, through the steps
    that are on it, so the distances are over the variants its filters kept,
    and the ``Variants`` is as it was afterwards.

    `min_num_snps` is how many variants a pair needs to get a distance: a
    pair that was called together at fewer has NaN in the vector, and one
    with exactly that many keeps its distance. ``None`` is 0. It keeps the
    name pyNei gives it, although the variants need not be SNPs. A pair that
    was called together at no variant has no distance whatever it says.

    What it gives is a :class:`Distances` with the distance of every pair,
    the names of the individuals and the counts of the pass in its
    ``pass_stats``.

    A pass that gives no variant is a ``ValueError``, whether the source has
    none or the steps kept none: the message says which of the two, and,
    when it was the steps, what each filter was given and kept.

    It is pyNei's ``calc_pairwise_kosman_dists``, with these differences:
    there is no `use_approx_embedding_algorithm`, which gives distances that
    are not Kosman distances and that change from one run to the next; there
    is no `num_threads`, because no calculation of popnei has it and the
    threads are those of rayon's pool, one for each core unless
    ``RAYON_NUM_THREADS`` is set before the first call; every ploidy is
    taken, where pyNei refuses all but 2; a negative `min_num_snps` is a
    ``ValueError``, where pyNei takes it and does with it what 0 does; a
    pass that gives no variant is a ``ValueError``, where pyNei raises a
    ``RuntimeError``, which in popnei is the exception of a defect of
    popnei; and the result carries the counts of the pass, which pyNei keeps
    in its ``Variants``.
    """
    if not isinstance(variants, Variants):
        # The path of the VCF is the mistake that is easiest to make, and
        # what it gave was the `AttributeError` of an object with no source
        # inside it.
        raise TypeError(
            f"`variants` is {variants!r}, of the type "
            f"`{type(variants).__name__}`, and the Kosman distances are "
            f"calculated over the variants of a source: "
            f"give it what `open_vcf` or `open_vars` gives, "
            f"calc_pairwise_kosman_dists(open_vcf(vcf_path))"
        )
    dist_vector, names, counts = _core.calc_pairwise_kosman_dists(
        variants._source, _min_num_vars_of(min_num_snps), variants._steps
    )
    return Distances(
        dist_vector=dist_vector, names=names, pass_stats=_pass_stats_of(counts)
    )


def _min_num_vars_of(min_num_snps: int | None) -> int:
    """How many variants a pair needs, as the core takes it: a whole number
    from 0 up, and 0 for the `None` a user writes for no minimum.

    The argument is refused here, where the name the user wrote is known:
    the core takes a count and has no word for `min_num_snps`.
    """
    if min_num_snps is None:
        return 0
    try:
        # A truth value is a whole number in Python, so `True` would pass as
        # the number 1, and it says nothing about how many variants a user
        # wants: it is refused here as the thresholds of the filters are.
        if isinstance(min_num_snps, bool):
            raise TypeError
        wanted = operator.index(min_num_snps)
    except TypeError:
        raise TypeError(
            f"`min_num_snps` is {min_num_snps!r}, and it says at how many "
            f"variants a pair has to have been called to get a distance: a "
            f"whole number of 0 or more, or `None` for no minimum"
        ) from None
    if wanted < 0:
        raise ValueError(
            f"`min_num_snps` is {wanted}, and it says at how many variants a "
            f"pair has to have been called to get a distance: a whole number "
            f"of 0 or more, or `None` for no minimum"
        )
    if wanted > _MOST_VARIANTS_A_PAIR_CAN_HAVE:
        raise ValueError(
            f"`min_num_snps` is {wanted}, and popnei counts at most "
            f"{_MOST_VARIANTS_A_PAIR_CAN_HAVE} variants for a pair: no pair "
            f"could reach it, and every distance would be missing"
        )
    return wanted


def _num_individuals_of(num_pairs: int) -> int:
    """How many individuals make `num_pairs` pairs, which is what the length
    of a distance vector says.

    N individuals make N (N - 1) / 2 pairs, so N is the whole number whose
    triangle is the length, and a length that is no triangle is of no set of
    individuals: pyNei takes it and builds the largest N below it, which
    leaves the last values of the vector in no cell of the matrix.
    """
    # N = (1 + sqrt(1 + 8 p)) / 2, worked out in whole numbers so that no
    # rounding of a square root decides how many individuals there are.
    under_the_root = 1 + 8 * num_pairs
    root = math.isqrt(under_the_root)
    if root * root != under_the_root or root % 2 == 0:
        raise ValueError(
            f"`dist_vector` holds {num_pairs} distances, and no number of "
            f"individuals makes that many pairs: N individuals make "
            f"N (N - 1) / 2 of them, {_the_pairs_around(num_pairs)}"
        )
    return (1 + root) // 2


def _the_pairs_around(num_pairs: int) -> str:
    """The numbers of pairs nearest to `num_pairs`, each with the
    individuals that make them, for the message of a vector whose length is
    of no set of individuals.

    `below` is the largest number of individuals whose pairs are not more
    than `num_pairs`, so the two numbers named are the one under the length
    that was given and the one over it.
    """
    below = (1 + math.isqrt(1 + 8 * num_pairs)) // 2
    return (
        f"{below * (below - 1) // 2} for {below} individuals and "
        f"{below * (below + 1) // 2} for {below + 1}"
    )


def _individuals(count: int) -> str:
    """`count` individuals, in the singular when there is one of them."""
    return f"{count} individual" if count == 1 else f"{count} individuals"


def _of_how_many_individuals(num_pairs: int) -> str:
    """How many individuals make `num_pairs` pairs, for the message of a
    number of names that is not that many.

    No distance is what one individual gives and what no individual would
    give, and nothing in the vector tells the two apart.
    """
    if num_pairs == 0:
        return "0 or of 1 individual"
    return _individuals(_num_individuals_of(num_pairs))


def _refuse_two_sides_that_differ(index: list, columns: list) -> None:
    """The first name of `index` that is not the name of the column at the
    same place, refused.

    # Raises

    ``ValueError`` when the two sides of the frame are not the same
    individuals in the same order.
    """
    for place, (row, column) in enumerate(zip(index, columns, strict=True)):
        if row != column:
            raise ValueError(
                f"the name of the row {place} is {row!r} and the name of its "
                f"column is {column!r}, and the two sides of a square matrix "
                f"of distances are the same individuals in the same order: "
                f"the cell of the row {row!r} and the column {column!r} is "
                f"the distance of that pair and not of a pair of one of "
                f"them. Give the frame `square_dists` gives, or sort both "
                f"sides of yours the same way"
            )
