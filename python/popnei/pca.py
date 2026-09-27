"""The principal component analysis: where the individuals fall on a few axes.

A principal component analysis places the individuals of a dataset on a few
axes that hold as much of the variation between them as that many axes can,
which is how a user sees whether their individuals fall into populations
before they name any. :func:`do_pca` does it on a table of numbers that the
user brings, individuals x traits, and :func:`do_pca_from_variants` on the
variants of a :class:`popnei.Variants`, where each variant is a trait and
the number of an individual at it is its dosage.

`docs/specs/pca.md` has what is computed. Each trait is centered, its mean
taken from it, and standardized, divided by its standard deviation, which
puts traits measured in different units on one scale; the components are
then the directions in the space of the traits along which the individuals
vary most, the first holding the largest variance that any direction has,
the second the largest among the directions at a right angle to the first,
and so on.

:func:`do_pcoa` places the individuals of a :class:`popnei.Distances` on
components in the same way, from the distance of every pair of them instead
of a table of their values: a principal coordinate analysis, PCoA. It refuses
distances that are not Euclidean, which no space has points at, and
:func:`correct_dists_by_lingoes` makes them Euclidean.
:func:`do_pcoa_from_variants` is the PCoA of the Kosman distances of the
individuals of a :class:`popnei.Variants`, in one pass over its variants.
"""

from dataclasses import dataclass

import numpy
import pandas

from popnei import _core
from popnei.dists import Distances, _min_num_vars_of
from popnei.variant import PassStats, Variants, _pass_stats_of

# What is wrong with a trait whose mean or whose standard deviation the
# analysis cannot use, under the name the binding crate gives each of the
# three. The values of such a trait are too large or too small for the
# arithmetic of a float64, and what the user does about it is to scale that
# trait or take it out.
_WHAT_IS_WRONG_WITH_THE_TRAIT = {
    "mean_not_finite": (
        "its values sum above the largest float64, 1.8e308, so its mean is not "
        "finite and every value of it would become a NaN"
    ),
    "deviation_not_finite": (
        "the squares of its deviations sum above the largest float64, so its "
        "standard deviation is not finite and the trait would become a column "
        "of zeros"
    ),
    "deviation_of_zero": (
        "the squares of its deviations all fall below the smallest float64 "
        "above 0, 5e-324, so its standard deviation is 0 although its values "
        "differ, and dividing by it would give infinities"
    ),
}


@dataclass(frozen=True)
class PCAResult:
    """What a principal component analysis gives.

    Only the components that have variance are here: centering takes one
    dimension out of the data, so a table of 8 individuals and 30 traits has
    7 components and not 8. They are named ``PC0``, ``PC1`` and so on, with
    zeros on the left to the width of how many there are, ``PC000`` for 200
    of them, and the same names are the columns of `projections`, the index
    of `explained_variance_percent` and the index of `princomps`.

    A component multiplied by -1 is the same component, and which of the two
    a decomposition gives depends on the library that did it. popnei fixes
    it: in each component the projection of the largest absolute value is
    positive, and the weights of that component have the sign that gave it.
    So the numbers are the same in Python, under pyodide and in TypeScript.
    """

    projections: pandas.DataFrame
    """Where each individual falls along each component, one row per
    individual and one column per component. The index is the index of the
    table, or the names of the individuals of the ``Variants``."""

    explained_variance_percent: pandas.Series
    """How much of the variance each component holds, as a percentage of the
    variance of every component there is, given or not."""

    princomps: pandas.DataFrame
    """The weight of each trait in each component, one row per component and
    one column per trait. The columns are the columns of the table, or the
    position of each variant that was used among the variants the pass gave,
    from 0. Of a table there is one row per component; of the variants there
    are ``num_prin_comps``, which can be fewer, and they are the first
    components."""

    pass_stats: PassStats | None
    """The counts of the pass over the source of variants, and ``None`` for
    the analysis of a table, which reads no variants."""


def do_pca(
    data: pandas.DataFrame,
    center_data: bool = _core.DEFAULT_CENTER_DATA,
    standardize_data: bool = _core.DEFAULT_STANDARDIZE_DATA,
) -> PCAResult:
    """The principal components of `data`, a table of individuals x traits.

    The index of `data` names the rows of the projections of the result and
    its columns name the columns of the weights. Every component that has
    variance is given, with its weights: the table is in memory already, so
    there is nothing to ask for fewer of. The result reads no variants, so
    its ``pass_stats`` is ``None``.

    `center_data` takes the mean of each trait from it. Without it the first
    component mostly points at the mean of the data.

    `standardize_data` then divides each trait by its standard deviation,
    which puts traits measured in different units on one scale; without it
    the traits with the largest numbers dominate. The divisor of that
    deviation is the number of individuals and not the number less one,
    which is pyNei's and makes every projection of a standardized table
    0.9975 of what R's ``prcomp`` gives at 200 individuals. Standardizing
    divides by the deviation the trait has once it is centered, so asking
    for it with `center_data` false is a ``ValueError``.

    No value may be missing: a table comes whole. A value that is not
    finite, an infinity or a NaN, and one that pandas holds as missing in a
    nullable dtype, are a ``ValueError`` that says which row and which trait
    it is at. So are a table of fewer than 2 rows or of no traits, and one
    in which no trait has variance once it is centered.

    Two more ``ValueError`` name the trait they are about, as the frame
    names it. One is a trait with no variance, every value of it equal to
    the others, when the table is standardized: there is nothing to divide
    it by, and the message says how many such traits there are and names the
    first ten of them, so that the user can take them out. Without
    standardizing such a trait is no error and gets a weight of 0. The other
    is a trait whose values are too large or too small for the arithmetic of
    a float64, whose message says which of those two it is.

    It is pyNei's ``do_pca``, which spells `standardize_data`
    ``standarize_data``, gives the components that have no variance as well,
    and leaves the sign of each component to the library that decomposed the
    table.
    """
    values = numpy.ascontiguousarray(
        # A frame of a nullable dtype holds its missing values as pandas's
        # own NA, which numpy cannot turn into a float64 on its own. It
        # arrives as a NaN, and the core says which row and which trait it
        # is at, instead of numpy raising a TypeError about a dtype.
        data.to_numpy(dtype=numpy.float64, na_value=numpy.nan),
        dtype=numpy.float64,
    )
    try:
        projections, percent, princomps = _core.pca(
            values, center_data, standardize_data
        )
    except _core.TraitsWithNoVariance as error:
        raise ValueError(
            _the_traits_with_no_variance(error.args[1], data.columns)
        ) from None
    except _core.TraitOutOfRange as error:
        raise ValueError(
            _the_trait_out_of_range(error.args[1], error.args[2], data.columns, error)
        ) from None
    names = _component_names(projections.shape[1])
    # `copy=False` on each of the three: without it pandas allocates a
    # second array of the same size and copies into it, which at 10000
    # individuals is 1.6 GB of peak memory against 2 MB. What the keyword
    # needs is that nothing else holds the array, and nothing does: each
    # one was allocated by the core crate for this call, `_core.pca` is
    # the only reference to it, and the frame is what outlives the call,
    # so no caller can see the frame and the array as two things.
    return PCAResult(
        projections=pandas.DataFrame(
            projections, index=data.index, columns=names, copy=False
        ),
        explained_variance_percent=pandas.Series(percent, index=names, copy=False),
        princomps=pandas.DataFrame(
            princomps, index=names, columns=data.columns, copy=False
        ),
        pass_stats=None,
    )


def do_pca_from_variants(
    variants: Variants,
    transform_to_biallelic: bool = _core.DEFAULT_TRANSFORM_TO_BIALLELIC,
    num_prin_comps: int = _core.DEFAULT_NUM_PRIN_COMPS,
) -> PCAResult:
    """The principal components of the variants of `variants`.

    Each variant becomes one number per individual, its dosage: how many
    alleles of the genotype are not the major allele of the variant, which
    is the most frequent among its called alleles and the lowest numbered of
    two that are equally frequent. A genotype with an allele missing takes
    the mean of the dosages of its variant, so that after centering it pulls
    its individual nowhere, and the standard deviation each variant is
    divided by has all the individuals in it and not the called ones alone.

    A variant whose called genotypes all have one dosage has no variance and
    is left out, a variant with one allele and one where every individual is
    heterozygous among them, and so is a variant with no called genotype.
    The ones that were used are the columns of the weights, by their
    position among the variants the pass gave, from 0.

    The names of the individuals are the index of the projections, and the
    counts of the pass are in ``pass_stats``: ``num_vars`` is how many
    variants the steps of the `Variants` let through, used or not, and
    ``filtering`` what each filter of the pass was given and kept.

    `transform_to_biallelic` makes every allele that is not the major one
    count the same, which is what a variant of more than two different
    alleles among its called genotypes needs: the dosage of a genotype has a
    meaning for two alleles, and without this such a variant is a
    ``ValueError`` that says which one it is. The alleles are those the
    genotypes hold and not those the source lists.

    `num_prin_comps` is how many components the weights are given for, 10 by
    default, and it is new here: pyNei gives the weight of every variant in
    every component, which is 0.8 GB for 100000 variants of 1000
    individuals. They come from a second pass over the variants, because a
    weight needs the eigenvectors, which are known when the first pass ends,
    so with 0 there is no second pass, ``princomps`` has no rows and it
    still has the variants that were used as its columns. More components
    than there are gives those there are, and a number below 0 is a
    ``ValueError``, as is one above what the machine counts; what is no
    whole number, a float, a string and a truth value among them, is a
    ``TypeError`` that names the argument and what was given. The
    projections and the percentages are of every component that has variance
    whatever `num_prin_comps` is.

    A dataset with no variant, one where no variant has variance, which one
    individual gives, and one with no individual are a ``ValueError``. So is
    a dataset of a size this analysis cannot count in, which is one of four:
    a ploidy above 254, more than 46340 individuals, more variants than the
    machine counts, which under pyodide is 4295 million, and weights of more
    values than the machine counts, which is `num_prin_comps` times the
    variants that were used and which fewer variants reach the more
    components are asked for.

    It is pyNei's ``do_pca_from_variants``, with `num_prin_comps` added,
    without `num_threads`, which no calculation of popnei takes, and with
    the filters of the ``Variants`` as steps of it and not as functions
    around it. pyNei counts the alleles of a whole chunk of variants and
    popnei those of each variant; pyNei gives the components with no
    variance as well; and a variant with no called genotype is an error in
    pyNei and is left out here.
    """
    if not isinstance(variants, Variants):
        # What a user gives instead is usually the path of the VCF, and
        # what that gave was the `AttributeError` of an object with no
        # source inside it.
        raise TypeError(
            f"`variants` is {variants!r}, a {type(variants).__name__}, and "
            f"`do_pca_from_variants` reads the variants of a source: give it "
            f"what `open_vcf` or `open_vars` gives, "
            f"do_pca_from_variants(open_vcf(vcf_path))"
        )
    # `num_prin_comps` is checked in the binding crate, where every count a
    # user writes is: a whole number of Python is of any size, and what is
    # none of them is refused there by the name of the argument.
    projections, percent, princomps, used_vars, counts = _core.pca_of_variants(
        variants._source, transform_to_biallelic, num_prin_comps, variants._steps
    )
    # The components of the projections are every one that has variance, and
    # the weights are those of the first `num_prin_comps` of them, so the
    # names are made for the projections and the weights take the first of
    # them: one component has one name in both frames.
    names = _component_names(projections.shape[1])
    # `copy=False` on each of the three, for the reason `do_pca` gives: the
    # three arrays come straight from the core crate, nothing but this call
    # holds them, and only the frames outlive it.
    return PCAResult(
        projections=pandas.DataFrame(
            projections, index=list(variants.individuals), columns=names, copy=False
        ),
        explained_variance_percent=pandas.Series(percent, index=names, copy=False),
        princomps=pandas.DataFrame(
            princomps,
            index=names[: princomps.shape[0]],
            columns=used_vars,
            copy=False,
        ),
        pass_stats=_pass_stats_of(counts),
    )


@dataclass(frozen=True)
class PCoAResult:
    """What a principal coordinate analysis gives.

    Only the components of the positive eigenvalues are here: the centering
    of the matrix the analysis decomposes always has an eigenvalue of 0,
    which gives no component, so 40 individuals have 39 components at most.
    They are named ``PC0``, ``PC1`` and so on, with zeros on the left to the
    width of how many there are, as the components of :class:`PCAResult`
    are. In each component the projection of the largest absolute value is
    positive, so the numbers are the same in Python, under pyodide and in
    TypeScript.
    """

    projections: pandas.DataFrame
    """Where each individual falls along each component, one row per
    individual, indexed by the names of the ``Distances``, and one column per
    component."""

    explained_variance_percent: pandas.Series
    """How much of the variance each component holds, as a percentage of the
    variance of the individuals placed at their distances, which is the sum
    of the squared distances over the pairs divided by the individuals. The
    percentages add up to 100."""

    lingoes_constant: float
    """The constant of Lingoes' correction, which :func:`do_pcoa` never
    makes, so it is 0 there."""

    negative_eigenvalues_percent: float
    """The share of the negative eigenvalues of the distances before a
    correction, which :func:`do_pcoa` refuses, so it is 0 there."""

    pass_stats: PassStats | None
    """The counts of the pass that gave the distances, those of the
    ``Distances`` given, and ``None`` for a ``Distances`` the user built."""


@dataclass(frozen=True)
class LingoesCorrection:
    """Distances made Euclidean by Lingoes' correction, with how much was
    added and how far they were from Euclidean before."""

    dists: Distances
    """The corrected distances, sqrt(d² + 2c) for every distance d of the
    ``Distances`` given, with its names and its ``pass_stats``. They have no
    standard errors: those of the distances given are not those of the
    corrected ones."""

    constant: float
    """c, the absolute value of the most negative eigenvalue of the matrix
    the analysis decomposes, in the units of a squared distance, and 0 when
    none is negative."""

    negative_eigenvalues_percent: float
    """100 times the sum of the absolute values of the negative eigenvalues
    over the sum of every eigenvalue, of the distances given: how large the
    part is that no space has, 0 for Euclidean distances."""


def do_pcoa(dists: Distances) -> PCoAResult:
    """The principal coordinates of `dists`, Gower's method.

    The individuals are put in a space where the straight line between each
    two of them is as long as their distance, and the components are the
    directions of that space along which they vary most, as those of
    :func:`do_pca` are for a table. With d the distance of a pair, the
    matrix decomposed is -d²/2 centered by rows and by columns, and the
    projections of a component are its eigenvector times the square root of
    its eigenvalue. It is what R's ``cmdscale`` and ``pcoa`` of the package
    ape compute, and of distances between the rows of a centered table it
    gives the projections of :func:`do_pca` of that table not standardized.

    Distances are Euclidean when some space has points at them, and then no
    eigenvalue is negative. The Kosman distances need not be, and distances
    that are not are a ``ValueError`` that says how many eigenvalues are
    negative and what share of the sum of all of them, and that
    :func:`correct_dists_by_lingoes` makes them Euclidean: ``do_pcoa`` of its
    ``dists`` is the analysis of the corrected distances.

    A pair with no distance, a NaN, is a ``ValueError`` that says how many
    there are, names the first and the individual in the most of them, since
    every individual is placed by its distance to every other. A distance
    that is negative or infinite is a ``ValueError`` that names its pair and
    gives it. So are fewer than 2 individuals, and
    distances that are all 0. What is not a ``Distances`` is a
    ``TypeError``.

    It is pyNei's ``do_pcoa``, which gives all n components: those of the
    negative eigenvalues with a negative percentage and the one of the
    eigenvalue 0 of the centering with projections of rounding, a negative
    distance squared as if it were positive, and the sign of each component
    as the decomposition left it. Of Euclidean distances pyNei's first n - 1
    components are these.
    """
    _refuse_what_is_not_a_distances(dists, "do_pcoa")
    try:
        projections, percent, lingoes_constant, negative_percent = _core.pcoa(
            dists.dist_vector, len(dists.names)
        )
    except _core.PcoaPairsWithNoDistance as error:
        raise ValueError(_the_pairs_with_no_distance(error, dists.names)) from None
    except _core.PcoaDistanceOutOfRange as error:
        raise ValueError(_the_distance_out_of_range(error, dists.names)) from None
    names = _component_names(projections.shape[1])
    # `copy=False` on both, for the reason `do_pca` gives: the arrays come
    # straight from the core crate and only the frames outlive this call.
    return PCoAResult(
        projections=pandas.DataFrame(
            projections, index=list(dists.names), columns=names, copy=False
        ),
        explained_variance_percent=pandas.Series(percent, index=names, copy=False),
        lingoes_constant=lingoes_constant,
        negative_eigenvalues_percent=negative_percent,
        pass_stats=dists.pass_stats,
    )


def do_pcoa_from_variants(
    variants: Variants,
    min_num_snps: int | None = None,
    correct_by_lingoes: bool = _core.DEFAULT_CORRECT_BY_LINGOES,
) -> PCoAResult:
    """The principal coordinates of the Kosman distances of the individuals
    of `variants`.

    It is :func:`calc_pairwise_kosman_dists` followed by :func:`do_pcoa`, in
    one pass over the variants the steps of `variants` keep and with no
    ``Distances`` built between the two. The Kosman distance counts every
    allele of a variant as itself and takes any ploidy, where
    :func:`do_pca_from_variants` reduces a genotype to its dosage.

    `min_num_snps` is how many variants a pair needs to have been called
    together at to get a distance, as :func:`calc_pairwise_kosman_dists`
    takes it: ``None`` or 0 is one. A pair with no distance is a
    ``ValueError`` that says how many there are, names the first and the
    individual in the most of them, and says that it can be taken out with
    :meth:`Variants.filter_individuals`, that `min_num_snps` can be lowered,
    or that :func:`do_pca_from_variants` gives every individual a
    projection.

    The Kosman distances need not be Euclidean. Without
    `correct_by_lingoes`, the default, distances that are not are a
    ``ValueError`` that says how many eigenvalues are negative and what share
    of the sum of all of them, and names `correct_by_lingoes`. With it,
    Lingoes' correction is applied inside the analysis, as
    :func:`correct_dists_by_lingoes` applies it to a ``Distances``: 2c is
    added to every squared distance, c being the absolute value of the most
    negative eigenvalue, and the result carries c as ``lingoes_constant``
    and the share of the negative eigenvalues before the correction as
    ``negative_eigenvalues_percent``, both 0 for Euclidean distances. The
    correction pushes every pair apart by the same amount of squared
    distance, which moves the nearest pairs the most. ``do_pcoa`` of
    ``correct_dists_by_lingoes`` of the distances gives the same components,
    within the rounding of the square roots.

    The names of the individuals are the index of the projections, and the
    counts of the pass are in ``pass_stats``. Fewer than 2 individuals, a
    pass that gives no variant and distances that are all 0 are a
    ``ValueError``; what is not a ``Variants`` is a ``TypeError``, and so are
    a `min_num_snps` that is not a whole number and a `correct_by_lingoes`
    that is not a bool or a ``numpy.bool_``.

    It is pyNei's ``do_pcoa_from_variants`` with `correct_by_lingoes` added
    and without `use_approx_embedding_algorithm` and `num_threads`. pyNei
    gives all n components, the negative ones included, and has no
    correction.
    """
    if not isinstance(variants, Variants):
        # What a user gives instead is usually the path of the VCF, as for
        # `do_pca_from_variants`.
        raise TypeError(
            f"`variants` is {variants!r}, a {type(variants).__name__}, and "
            f"`do_pcoa_from_variants` reads the variants of a source: give it "
            f"what `open_vcf` or `open_vars` gives, "
            f"do_pcoa_from_variants(open_vcf(vcf_path))"
        )
    min_num_vars = _min_num_vars_of(min_num_snps)
    if not isinstance(correct_by_lingoes, bool | numpy.bool_):
        # pyo3 refuses it with `'int' object is not an instance of 'bool'`,
        # which names neither the argument nor the call.
        raise TypeError(
            f"`correct_by_lingoes` is True or False, and {correct_by_lingoes!r}, "
            f"of the type {type(correct_by_lingoes).__name__}, was given"
        )
    # The core names an individual of a pair with no distance by its
    # position among those the pass gives, which are these.
    individuals = variants.individuals
    try:
        (
            projections,
            percent,
            lingoes_constant,
            negative_percent,
            counts,
        ) = _core.pcoa_of_variants(
            variants._source,
            min_num_vars,
            bool(correct_by_lingoes),
            variants._steps,
        )
    except _core.PcoaPairsWithNoDistance as error:
        raise ValueError(
            f"{variants._source.path()}: "
            + _the_pairs_with_no_distance(error, individuals)
        ) from None
    names = _component_names(projections.shape[1])
    # `copy=False` on both, for the reason `do_pca` gives.
    return PCoAResult(
        projections=pandas.DataFrame(
            projections, index=list(individuals), columns=names, copy=False
        ),
        explained_variance_percent=pandas.Series(percent, index=names, copy=False),
        lingoes_constant=lingoes_constant,
        negative_eigenvalues_percent=negative_percent,
        pass_stats=_pass_stats_of(counts),
    )


def correct_dists_by_lingoes(dists: Distances) -> LingoesCorrection:
    """`dists` made Euclidean by Lingoes' correction.

    It adds 2c to the square of every distance, c being the absolute value of
    the most negative eigenvalue of the matrix :func:`do_pcoa` decomposes,
    which is ape's ``pcoa(d, correction = "lingoes")`` of R (Lingoes 1971).
    The corrected distances have the eigenvectors of the ones given, and
    every eigenvalue but the 0 of the centering is c larger, so the most
    negative becomes 0 and none is below it. Distances that are Euclidean
    already give a constant of 0 and the same distances.

    The correction pushes every pair apart by the same amount of squared
    distance, which moves the nearest pairs the most, so the result carries
    c and the share of the negative eigenvalues before it, for a user to
    judge how much the picture was changed by.

    It refuses what :func:`do_pcoa` refuses but distances that are not
    Euclidean. The constant is in the units of a squared distance, and
    distances of about 1.3e154 and above give one above the largest float64,
    and distances all below about 1e-154 one below the smallest normal
    float64, 2.2e-308, which keeps only a few of its digits or is 0: both are
    a ``ValueError`` that says to divide the distances by a number near the
    largest first.

    pyNei has no correction.
    """
    _refuse_what_is_not_a_distances(dists, "correct_dists_by_lingoes")
    try:
        dist_vector, constant, negative_percent = _core.correct_dists_by_lingoes(
            dists.dist_vector, len(dists.names)
        )
    except _core.PcoaPairsWithNoDistance as error:
        raise ValueError(_the_pairs_with_no_distance(error, dists.names)) from None
    except _core.PcoaDistanceOutOfRange as error:
        raise ValueError(_the_distance_out_of_range(error, dists.names)) from None
    return LingoesCorrection(
        dists=Distances(
            dist_vector=dist_vector, names=dists.names, pass_stats=dists.pass_stats
        ),
        constant=constant,
        negative_eigenvalues_percent=negative_percent,
    )


def _refuse_what_is_not_a_distances(dists: object, function: str) -> None:
    """A `dists` that is not a :class:`popnei.Distances`, refused.

    # Raises

    ``TypeError`` that names `function` and says how to build a
    ``Distances`` from the square frame, which is what a user of pyNei's
    ``square_dists`` has at hand.
    """
    if not isinstance(dists, Distances):
        raise TypeError(
            f"`dists` is of the type `{type(dists).__name__}`, and `{function}` "
            f"takes a `Distances`, what `calc_pairwise_kosman_dists` gives: "
            f"build one from a square frame of distances with "
            f"`Distances.from_square_dists`, or from the distances of the pairs "
            f"with `Distances(dist_vector, names=...)`"
        )


def _the_pairs_with_no_distance(of_the_core: BaseException, names: tuple) -> str:
    """What a user is told of the pairs of a ``Distances`` with no distance.

    The core names the individuals by their position in the order of the
    distances, which is where the names of the ``Distances`` are read, says
    the counts, and writes what the user is to do, which depends on whether
    the distances were given or came from the variants: that text is taken
    as the core wrote it.
    """
    (
        _,
        num_pairs_with_no_distance,
        num_pairs,
        first,
        second,
        most_often,
        most_often_count,
        remedy,
    ) = of_the_core.args
    have = "has" if num_pairs_with_no_distance == 1 else "have"
    return (
        f"{num_pairs_with_no_distance} of the {num_pairs} pairs of individuals "
        f"{have} no distance, the first of them {names[first]!r} and "
        f"{names[second]!r}, and {names[most_often]!r} is in {most_often_count} "
        f"of them; {remedy}"
    )


def _the_distance_out_of_range(of_the_core: BaseException, names: tuple) -> str:
    """What a user is told of a distance of a ``Distances`` that is negative
    or infinite.

    The core names the two individuals of the pair by their position in the
    order of the distances, which is where the names of the ``Distances``
    are read, gives the distance, and writes what a distance has to be,
    which is taken as the core wrote it.
    """
    _, first, second, value, what_it_has_to_be = of_the_core.args
    return (
        f"the distance of {names[first]!r} and {names[second]!r} is {value}, and "
        f"{what_it_has_to_be}"
    )


def _component_names(num_comps: int) -> list[str]:
    """`PC0`, `PC1` and so on, one for each component.

    The number has zeros on its left to the width of `num_comps`, ``PC000``
    for 200 components, which is what pyNei's `_create_pc_names` gives and
    what makes the names sort as their numbers do.
    """
    width = len(str(num_comps))
    return [f"PC{number:0{width}d}" for number in range(num_comps)]


def _the_traits_with_no_variance(
    positions: list[int], trait_names: pandas.Index
) -> str:
    """What a user is told of the traits that cannot be standardized.

    They are named as the columns of their frame, the first ten of them and
    how many more there are, so that a user of a table of hundreds of traits
    reads a message of one line and knows how many to take out. The core has
    said the same with the position of each trait in the place of its name.
    """
    shown = ", ".join(f"`{trait_names[position]}`" for position in positions[:10])
    left_out = len(positions) - len(positions[:10])
    more = "" if left_out == 0 else f", and {left_out} more"
    return (
        f"{len(positions)} of the {len(trait_names)} traits have no variance and "
        f"cannot be standardized: take them out of the table or do not standardize; "
        f"they are the traits {shown}{more}"
    )


def _the_trait_out_of_range(
    position: int,
    problem: str,
    trait_names: pandas.Index,
    of_the_core: BaseException,
) -> str:
    """What a user is told of a trait the analysis cannot scale.

    The trait is named as its column of the frame is, and the message says
    which of the three things happened to it, so that the user knows whether
    to scale it or to take it out. A problem this package has no words for,
    which a core newer than it can give, is passed on as the core said it,
    with the position of the trait in the place of its name.
    """
    said = _WHAT_IS_WRONG_WITH_THE_TRAIT.get(problem)
    if said is None:
        return str(of_the_core.args[0])
    return (
        f"the trait `{trait_names[position]}` cannot be centered or standardized: "
        f"{said}; scale that trait or take it out of the table"
    )
