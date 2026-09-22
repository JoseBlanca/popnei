"""The principal component analysis: where the individuals fall on a few axes.

A principal component analysis places the individuals of a dataset on a few
axes that hold as much of the variation between them as that many axes can,
which is how a user sees whether their individuals fall into populations
before they name any. :func:`do_pca` does it on a table of numbers that the
user brings, individuals x traits.

`docs/specs/pca.md` has what is computed. Each trait is centered, its mean
taken from it, and standardized, divided by its standard deviation, which
puts traits measured in different units on one scale; the components are
then the directions in the space of the traits along which the individuals
vary most, the first holding the largest variance that any direction has,
the second the largest among the directions at a right angle to the first,
and so on.
"""

from dataclasses import dataclass

import numpy
import pandas

from popnei import _core
from popnei.variant import PassStats

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
    individual, named as the index of the table was, and one column per
    component."""

    explained_variance_percent: pandas.Series
    """How much of the variance each component holds, as a percentage of the
    variance of every component there is, given or not."""

    princomps: pandas.DataFrame
    """The weight of each trait in each component, one row per component and
    one column per trait, named as the columns of the table were."""

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
    return PCAResult(
        projections=pandas.DataFrame(projections, index=data.index, columns=names),
        explained_variance_percent=pandas.Series(percent, index=names),
        princomps=pandas.DataFrame(princomps, index=names, columns=data.columns),
        pass_stats=None,
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
