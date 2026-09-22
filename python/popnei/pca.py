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
    center_data: bool = True,
    standardize_data: bool = True,
) -> PCAResult:
    """The principal components of `data`, a table of individuals x traits.

    The index of `data` names the rows of the projections and its columns
    name those of the weights. Every component that has variance is given,
    with its weights: the table is in memory already, so there is nothing to
    ask for fewer of.

    :param data: The values, one row per individual and one column per
        trait. No value may be missing.
    :param center_data: Whether the mean of each trait is taken from it.
        Without it the first component mostly points at the mean of the
        data.
    :param standardize_data: Whether each trait is then divided by its
        standard deviation, the one with the number of individuals in it and
        not the number less one, which is pyNei's. Without it the traits
        with the largest numbers dominate.
    :raises ValueError: When a value is not finite; when `standardize_data`
        is asked for and `center_data` is not; when the table has fewer than
        2 rows or no traits; and, when it is standardized, when a trait has
        no variance, the message naming how many they are and the first ten
        of them. Such a trait is no error without standardizing and gets a
        weight of 0.
    """
    values = numpy.ascontiguousarray(data.to_numpy(), dtype=numpy.float64)
    try:
        projections, percent, princomps = _core.pca(
            values, center_data, standardize_data
        )
    except _core.TraitsWithNoVariance as error:
        # The core gives the positions of those traits, since it has no
        # names, and the message the user reads names them as the frame
        # does. The exception of the core is not chained under it: it says
        # the same thing with numbers in the place of the names.
        raise ValueError(
            _the_traits_with_no_variance(error.args[0], data.columns)
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
    reads a message of one line and knows how many to take out.
    """
    shown = ", ".join(f"`{trait_names[position]}`" for position in positions[:10])
    left_out = len(positions) - len(positions[:10])
    more = "" if left_out == 0 else f", and {left_out} more"
    return (
        f"{len(positions)} of the {len(trait_names)} traits have no variance and "
        f"cannot be standardized: take them out of the table or do not standardize; "
        f"they are the traits {shown}{more}"
    )
