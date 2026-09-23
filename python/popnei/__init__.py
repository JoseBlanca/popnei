"""Population genetics over the variants of a dataset.

popnei reads variants, holds them in blocks and calculates over them. Every
calculation is written in Rust, in the core crate, which this package
reaches through the private module :mod:`popnei._core`. The functions that
open a dataset and the result objects are being written; the names and the
signatures are those of pyNei, the Python library popnei succeeds.
"""

from popnei import _core
from popnei.block import Block
from popnei.dists import Distances, calc_pairwise_kosman_dists
from popnei.filters import FilteringStats, Step
from popnei.io_vars import VarsWritten, open_vars, write_vars
from popnei.io_vcf import open_vcf
from popnei.kinship import Kinship, calc_kinship
from popnei.ld import R2Matrix, calc_rogers_huff_r2_matrix
from popnei.pca import PCAResult, do_pca, do_pca_from_variants
from popnei.pop_dists import PopDistMeasure, PopDists, calc_pop_dists
from popnei.stats import (
    PerIndividualStats,
    PerVarDistribs,
    PerVarStat,
    PolyVarsStats,
    StatsDistrib,
    calc_per_individual_stats,
    calc_per_var_distribs,
)
from popnei.variant import Blocks, PassStats, Variants

__version__: str = _core.version()
"""The version of the core crate, ``major.minor.patch``.

A user who reports a result names the code that gave it with this, so it is
the version of the Rust that did the work and not one written again here.
"""

__all__ = [
    "Block",
    "Blocks",
    "Distances",
    "FilteringStats",
    "Kinship",
    "PCAResult",
    "PassStats",
    "PerIndividualStats",
    "PerVarDistribs",
    "PerVarStat",
    "PolyVarsStats",
    "PopDistMeasure",
    "PopDists",
    "R2Matrix",
    "StatsDistrib",
    "Step",
    "Variants",
    "VarsWritten",
    "__version__",
    "calc_kinship",
    "calc_pairwise_kosman_dists",
    "calc_per_individual_stats",
    "calc_per_var_distribs",
    "calc_pop_dists",
    "calc_rogers_huff_r2_matrix",
    "do_pca",
    "do_pca_from_variants",
    "open_vars",
    "open_vcf",
    "write_vars",
]
