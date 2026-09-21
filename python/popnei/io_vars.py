"""Writing a vars file, the file popnei keeps its variants in."""

from pathlib import Path

from popnei import _core
from popnei.variant import Variants


def write_vars(
    variants: Variants, path: str | Path, num_vars_per_block: int | None = None
) -> None:
    """Every variant of `variants` into a vars file at `path`.

    A vars file is one arrow IPC file, also called feather v2, which pandas,
    R and polars open as a table with no popnei installed. It is where a
    user keeps their variants once the VCF has been read, so that the text
    is parsed once and every later pass reads a file of arrays.

    The call reads the whole source once. The file holds the six columns of
    a VCF, the chromosome, the position, the id, the alleles, the quality
    and the genotypes, whether or not the user will read them, so that it
    can stand in for the VCF in any later analysis; a source that has no
    alleles to give, an array of genotypes, gives a file without that
    column.

    `num_vars_per_block` is how many variants a batch of the file holds, the
    last one aside, and ``None`` asks for the size popnei chooses for the
    number of individuals of the source, which is the size of its blocks.

    A path that a file is already at is a ``ValueError`` and nothing is
    written. When the source fails half way, on a line of a VCF that popnei
    cannot read, that error is raised, and the file that was being written
    is taken away, so that the same call can be made again at the same path
    once the VCF is fixed; pyNei leaves what it had written. A Ctrl-C is
    raised when the pass over the source is over and not while it runs,
    because the file is written inside one call of the Rust core, and it
    leaves no file at the path either.

    It is pyNei's ``write_vars`` with the same first two arguments, and
    these differences: the file is another one, which pyNei does not read;
    `num_vars_per_block` is an argument here, where in pyNei the size
    belongs to the ``Variants``; a source with no variants is written, where
    pyNei raises a ``ValueError``; the file of a call that failed is taken
    away; and every column is written, where pyNei writes the ones its
    chunks happen to carry.
    """
    _core.write_vars(variants._source, path, num_vars_per_block)
