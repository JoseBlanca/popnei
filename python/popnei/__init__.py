"""Population genetics over the variants of a dataset.

popnei reads variants, holds them in blocks and calculates over them. Every
calculation is written in Rust, in the core crate, which this package
reaches through the private module :mod:`popnei._core`. The functions that
open a dataset and the result objects are being written; the names and the
signatures are those of pyNei, the Python library popnei succeeds.
"""

from popnei import _core
from popnei.block import Block
from popnei.io_vars import write_vars
from popnei.io_vcf import open_vcf
from popnei.variant import Variants

__version__: str = _core.version()
"""The version of the core crate, ``major.minor.patch``.

A user who reports a result names the code that gave it with this, so it is
the version of the Rust that did the work and not one written again here.
"""

__all__ = ["Block", "Variants", "__version__", "open_vcf", "write_vars"]
