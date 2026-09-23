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
from pathlib import Path

import numpy
import pandas

from popnei import _core
from popnei.pca import _component_names
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
    pyNei keeps in its ``Variants``; and the matrix is checked where pyNei
    checks nothing, so that what is wrong with it is said of the field the
    user filled and not later, out of the linear algebra, of a matrix and a
    row.
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
        """What makes a frame of numbers a kinship, checked.

        # Raises

        ``TypeError`` when `matrix` is not a pandas frame and when a row or
        a column of it is labelled with what is no name, which are the two
        faults of a type. ``ValueError`` when it holds what is no number;
        when it is not square; when it has no row at all; when its index and
        its columns are not the same individuals in the same order; when an
        individual is named twice; when an entry is not finite; and when it
        is further from its own transpose than 1e-9 of its largest absolute
        entry.
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
        if rows == 0:
            # A frame of no row is square and names nobody twice, so every
            # other check here lets it past, and only
            # `principal_components` would have complained. TypeScript
            # refuses it in its constructor and the two packages refuse the
            # same matrices.
            raise ValueError(
                "the frame has no row, and a kinship is the matrix of every "
                "pair of a set of individuals: there is no pair in it"
            )
        names = list(self.matrix.index)
        columns_of = list(self.matrix.columns)
        _refuse_a_label_that_is_no_name(names, "row")
        _refuse_a_label_that_is_no_name(columns_of, "column")
        _refuse_two_sides_that_differ(names, columns_of)
        _refuse_an_individual_that_is_there_twice(names)
        _refuse_a_matrix_that_is_no_kinship(self.matrix, names)

    @property
    def individuals(self) -> tuple[str, ...]:
        """The names of the individuals, in the order of the rows and of the
        columns of the matrix.

        They are read from the matrix itself, so the names are in one place
        and cannot disagree with themselves, as pyNei's ``samples`` is."""
        return tuple(self.matrix.index)

    def principal_components(self, num_pcs: int) -> pandas.DataFrame:
        """Where each individual falls along the `num_pcs` directions in
        which the panel varies most, taken from this matrix.

        A user gives them to an association study as covariates, which is how
        the structure of a panel is accounted for without a mixed model, and
        they cost an eigendecomposition of a matrix that is already in hand
        rather than a second pass over the variants. With ``lambda_j`` the
        eigenvalues of the kinship from the largest and ``u_j`` its
        eigenvectors, the component ``j`` is ``u_j * sqrt(lambda_j)``.

        What comes back is a frame with the names of the individuals as index
        and the components as columns, named ``PC0``, ``PC1`` and so on with
        zeros on the left to the width of how many there are, as
        :func:`popnei.do_pca` names them. It is indexed by individual so that
        it can be joined to the covariates of the association study.

        In every component the projection of the largest absolute value is
        positive, which is the rule of `docs/specs/pca.md`: a component
        multiplied by -1 is the same component, and without the rule the two
        backends of the eigendecomposition and the three builds of popnei
        would give different signs for one dataset.

        A component whose eigenvalue is not above ``lambda_1 * n * 2.2e-16``,
        with ``n`` the individuals, is not given, so asking for more
        components than the matrix has gives those it has and the frame has
        that many columns. A kinship measures a pair against the average pair
        of the panel, which takes one direction out of it, so a panel of 200
        individuals has 199 components and not 200. A `num_pcs` of 0 gives a
        frame of no columns and is no error, as asking a principal component
        analysis for no components is not; a negative one, and what is no
        whole number, are a ``ValueError`` and a ``TypeError`` that name the
        argument. A matrix of no individual, and one holding a value that is
        not finite, are a ``ValueError`` too: the checks of ``__post_init__``
        are made when the object is built and the frame it holds can be
        written into afterwards, so they are made again here, and the second
        of the two names the row and the column of the value.

        The lower half of the matrix is what the components are taken from,
        since it is symmetric: a value written above its diagonal after the
        object was built changes no component, and is refused only for not
        being finite.

        These are close to the principal components of the variants the
        kinship was calculated from and they are not the same, because the
        two standardize by different numbers: the kinship divides each
        variant by ``sqrt(ploidy * p * (1 - p))`` and
        :func:`popnei.do_pca_from_variants` by the standard deviation of its
        dosages, and the two agree only when the genotypes are in Hardy
        Weinberg proportions.

        It is pyNei's ``Kinship.principal_components``, with two differences.
        pyNei gives exactly `num_pcs` components whatever their eigenvalue,
        taking the square root of the absolute value of one below 0, which
        the per pair denominators of a dataset with genotypes missing put
        there: on the reference panel with 3 in 100 of its genotypes missing
        the smallest eigenvalue is -0.0321, and the length that absolute
        value gives means nothing. And pyNei leaves the sign of a component
        to the library that decomposed the matrix.
        """
        # A frame of one dtype lies column after column, and the binding
        # crate refuses an array that does not lie row after row, naming
        # `matrix` and `numpy.ascontiguousarray`: without this line a user
        # gets that `ValueError` for a frame of theirs that is right. The
        # copy is what the core works in, and reading it as it lies would be
        # reading the transpose, which for a matrix that is symmetric only
        # within the tolerance of `__post_init__` is other numbers.
        values = numpy.ascontiguousarray(
            self.matrix.to_numpy(dtype=numpy.float64), dtype=numpy.float64
        )
        # `num_pcs` is checked in the binding crate, where every count a
        # user writes is: a whole number of Python is of any size, and what
        # is none of them is refused there by the name of the argument.
        projections, num_comps = _core.kinship_principal_components(values, num_pcs)
        # `copy=False`: the array came from the core crate for this call,
        # nothing else holds it, and only the frame outlives the call.
        return pandas.DataFrame(
            projections,
            index=list(self.individuals),
            columns=_component_names(num_comps),
            copy=False,
        )

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
        that names it, as pyNei's ``filter_samples`` raises, and so is a name
        given twice, which would leave two rows where one was asked for. One
        name written where a sequence of them is meant, ``"s000"``, and a
        sequence of no name at all are a ``TypeError`` and a ``ValueError``
        that say so: a string is a sequence of its letters, and a kinship is
        of one individual at least.
        """
        named = _the_names_of(individuals)
        missing = [name for name in named if name not in self.matrix.index]
        if missing:
            raise ValueError(
                f"{_that_are_not_individuals(missing)} of the kinship, which "
                f"is of {len(self.matrix.index)} individuals: every name of "
                f"`individuals` is one of them"
            )
        _refuse_an_individual_that_is_there_twice(named)
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
    variance among those 40. A name that is not an individual of `variants`
    and a name given twice are a ``ValueError``; one name written where a
    sequence of them is meant, ``"s000"``, and a sequence of no name at all
    are refused as they are in :meth:`Kinship.filter_individuals`.

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
    named = None if individuals is None else _the_names_of(individuals)
    try:
        # The names go to the core as they are: it turns them into the places
        # of those individuals among the ones the pass gives, and refuses a
        # name that is not one of them, as the filter of individuals does.
        matrix, names, num_vars, counts = _core.calc_kinship(
            variants._source, named, transform_to_biallelic, variants._steps
        )
    except _core.KinshipPairWithNoVariantCalled as error:
        # The core has where the two individuals are among the ones the
        # kinship was asked for, and this layer has their names.
        raise ValueError(
            _the_pair_with_no_variant_called(
                error,
                named if named is not None else list(variants.individuals),
                variants._source.path(),
            )
        ) from None
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


def _the_names_of(individuals: Sequence[str]) -> list:
    """The names of `individuals` as a list, refused when they are no
    sequence of names.

    # Raises

    ``TypeError`` for one name written where a sequence of them is meant and
    for what cannot be walked at all, and ``ValueError`` for a sequence of no
    name: a kinship is of one individual at least, and both calls that take
    `individuals` name the ones the result is of.
    """
    if isinstance(individuals, str):
        # A string is a sequence of its letters, so one name written without
        # its comma asks for the individuals `s`, `0`, `0` and `0`, and a
        # name of one letter would be taken and give a kinship of one.
        raise TypeError(
            f"`individuals` is a sequence of names and not one name: write "
            f'individuals=("{individuals}",) for that one individual'
        )
    try:
        named = list(individuals)
    except TypeError:
        # What Python says of its own here, `'int' object is not iterable`,
        # names neither the argument nor the call.
        raise TypeError(
            f"`individuals` is a sequence of the names of the individuals the "
            f"kinship is of, and {individuals!r}, a "
            f"{type(individuals).__name__}, was given"
        ) from None
    if not named:
        raise ValueError(
            "`individuals` names no individual, and a kinship is of one "
            "individual at least: name the ones it is of"
        )
    return named


def _refuse_a_label_that_is_no_name(labels: list, side: str) -> None:
    """The first label of one side of the matrix that is no name, refused.

    `side` is ``"row"`` or ``"column"``, the one the labels are of, so that
    a user reads which of the two to write.

    # Raises

    ``TypeError`` naming the label, what type it is of and what to write:
    what is wrong with a label of 0 is its type, as it is for a `matrix`
    that is no frame, and a bare string where a sequence of names is meant
    is a ``TypeError`` for the same reason.
    A frame built as ``pandas.DataFrame(matrix)`` is labelled with the
    numbers 0 to N-1, which is how a user meets this, and a kinship whose
    rows are numbers cannot be matched to the phenotypes of the association
    study later. It is also what makes both packages take the same
    matrices: ``calcKinship`` of TypeScript takes the names as strings.
    """
    for place, label in enumerate(labels):
        if not isinstance(label, str):
            raise TypeError(
                f"the {side} {place} of the matrix is labelled {label!r}, of "
                f"the type `{type(label).__name__}`, and the individuals of a "
                f"kinship are named by strings: a frame written as "
                f"`pandas.DataFrame(matrix)`, with no index and no columns of "
                f"its own, is labelled with the numbers 0 to N-1, and a "
                f"kinship whose {side}s are numbers cannot be matched to the "
                f"phenotypes of `calc_gwas` later. Write "
                f"`pandas.DataFrame(matrix, index=individuals, "
                f"columns=individuals)` with the names of the individuals"
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


def _refuse_an_individual_that_is_there_twice(names: list) -> None:
    """The first name that is in `names` twice, refused.

    # Raises

    ``ValueError`` naming it and the two places it is at. A kinship has one
    row and one column for each of its individuals, and a name that is there
    twice makes every lookup by it, `filter_individuals` among them, give two
    rows where one was asked for.
    """
    first_at: dict = {}
    for place, name in enumerate(names):
        if name in first_at:
            raise ValueError(
                f"the individual {name!r} is named twice, at the places "
                f"{first_at[name]} and {place}, and a kinship has one row and "
                f"one column for each of its individuals: with the name there "
                f"twice, asking for that individual gives two rows and two "
                f"columns, which are the kinship of no pair"
            )
        first_at[name] = place


def _refuse_a_matrix_that_is_no_kinship(matrix: pandas.DataFrame, names: list) -> None:
    """A matrix that holds what is no number, a value that is not finite, or
    two cells of one pair that differ, refused.

    The rows are read one at a time, and what is allocated is a row and not a
    second matrix: the kinship of 10000 individuals is 800 MB, and a check
    written as `abs(values - values.T)` asks this machine for 1.6 GB more of
    it.

    The pair that is furthest from its transpose is the one named, with the
    two numbers the matrix holds for it, so that a user who built the frame
    from half a matrix sees which cell was not filled.

    # Raises

    ``ValueError`` for each of the three.
    """
    try:
        values = matrix.to_numpy(dtype=numpy.float64)
    except (TypeError, ValueError) as problem:
        raise ValueError(
            f"`matrix` holds what is no number, and every entry of a kinship "
            f"is the kinship of a pair of individuals: {problem}"
        ) from None
    for row in range(values.shape[0]):
        finite = numpy.isfinite(values[row])
        if not finite.all():
            column = int(numpy.argmin(finite))
            raise ValueError(
                f"the cell of the row {names[row]!r} and the column "
                f"{names[column]!r} holds {float(values[row, column])}, and "
                f"every entry of a kinship is a number: pyNei leaves a NaN "
                f"where a pair of individuals has no variant called in both "
                f"of them, and such a pair is taken out of the matrix before "
                f"it is a kinship"
            )
    largest = 0.0
    widest = 0.0
    furthest = (0, 0)
    for row in range(values.shape[0]):
        line = values[row]
        largest = max(largest, float(numpy.abs(line).max(initial=0.0)))
        gaps = numpy.abs(line - values[:, row])
        gap = float(gaps.max(initial=0.0))
        if gap > widest:
            widest = gap
            furthest = (row, int(numpy.argmax(gaps)))
    if widest <= largest * _LARGEST_ASYMMETRY:
        return
    row, column = furthest
    raise ValueError(
        f"the cell of the row {names[row]!r} and the column {names[column]!r} "
        f"holds {float(values[row, column])} and the cell of the row "
        f"{names[column]!r} and the column {names[row]!r} holds "
        f"{float(values[column, row])}, and the kinship of a pair is one "
        f"number, which both of its cells hold: the two differ by {widest}, "
        f"where {_LARGEST_ASYMMETRY} of the largest absolute entry of the "
        f"matrix, {largest}, is what a kinship is taken to be symmetric within"
    )


def _the_pair_with_no_variant_called(
    of_the_core: BaseException, individuals: list, path: Path
) -> str:
    """What a user is told of two individuals with no variant called in both.

    The core names them by their place among the individuals of the kinship,
    which is where the names of this layer are read, and says the counts. The
    two places are one place twice when an individual has no called genotype
    at all, which is a sequencing that failed and which is said of that one
    individual.
    """
    _, one, other, num_vars_of_one, num_vars_of_other = of_the_core.args
    if one == other:
        return (
            f"{path}: the individual {individuals[one]!r} has no called "
            f"genotype among the variants that were used, so its entry of the "
            f"kinship would be divided by no variant at all; leave it out"
        )
    said = "variant is" if num_vars_of_one == 1 else "variants are"
    return (
        f"{path}: the individuals {individuals[one]!r} and "
        f"{individuals[other]!r} have no variant called in both of them, so "
        f"their entry of the kinship would be divided by no variant at all: "
        f"{num_vars_of_one} {said} called in the first and {num_vars_of_other} "
        f"in the second; leave one of the two out"
    )


def _that_are_not_individuals(missing: list) -> str:
    """The individuals of a call that are not in the matrix, as the start of
    the message that refuses them, in the singular when there is one."""
    named = ", ".join(repr(name) for name in missing)
    if len(missing) == 1:
        return f"{named} is not an individual"
    return f"{named} are not individuals"
