"""The filters that keep the variants of a dataset that pass a threshold.

A filter is a step of a :class:`popnei.Variants`: a method of it that adds
itself to the list of steps and returns nothing, and that every pass over
the source runs inside the Rust core. The three threshold filters are
:meth:`popnei.Variants.filter_by_missing_data`, over the missing rate of a
variant, :meth:`popnei.Variants.filter_by_maf`, over its major allele
frequency, and :meth:`popnei.Variants.filter_by_obs_het`, over its observed
heterozygosity, and each of them keeps the variants whose number is at most
the threshold it was given. :meth:`popnei.Variants.filter_individuals` keeps
individuals and not variants: it takes the genotypes of the individuals a
user names, at every variant, and drops those of the rest.

What is here is what a user reads of them: the :class:`Step` that a filter
is in the steps of a ``Variants``, and the :class:`FilteringStats` that the
counts of a pass hold for each filter of it.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class FilteringStats:
    """How many variants one filter of one pass was given and kept.

    They are the counts of that pass alone, one reading of the source from
    its start: a `Variants` can be given to any number of consumers, and
    each of them counts its own.
    """

    vars_processed: int
    """The variants the filter was given, which are those that the filter
    before it in the steps kept, and all of them for the first filter."""

    vars_kept: int
    """Those of them that passed its threshold, which are the ones the
    filter after it was given."""


@dataclass(frozen=True)
class Step:
    """One step of a `Variants`: what every pass over its source runs.

    A filter is the only kind of step there is: one of the three over a
    number of a variant, or the one that keeps the individuals a user names.
    """

    kind: str
    """What the step does: ``"missing_data"``, ``"maf"``, ``"obs_het"`` or
    ``"individuals"``. The kind of a threshold filter is the name its counts
    have in the counts of a pass, where the filter of individuals has no
    entry, since it takes no variant away."""

    args: dict[str, object]
    """What the step was given, under the names of the arguments of the
    method that added it, ``{"max_allowed_maf": 0.95}`` for a threshold
    filter and ``{"individuals": ("ind05", "ind00")}`` for the filter of
    individuals, whose names are a tuple in the order they were given."""
