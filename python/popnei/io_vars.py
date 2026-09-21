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
    written, and a path that no file can be made at, a directory or a path
    in a directory that is not there, is the ``OSError`` the file system
    gives for it, with the path in ``filename``.

    Every error names the file it is about. A line of the VCF that popnei
    cannot read is a ``ValueError`` whose message starts with the path of
    the VCF; a vars file that could not be written, a disc that filled up
    among the causes, is an ``OSError`` that carries the path of the vars
    file in ``filename``. The file that was being written is taken away
    then, so that the same call can be made again at the same path once
    what was wrong is fixed, where pyNei leaves what it had written; when
    it cannot be taken away, the exception carries a note that says that a
    file is still there. The bytes of a file that was written reach the
    disc before the call returns, so a file system that says only at the
    close that it is full is an error too.

    A Ctrl-C is raised when the pass over the source is over and not while
    it runs, because the file is written inside one call of the Rust core,
    and it leaves no file at the path either.

    It is pyNei's ``write_vars`` with the same first two arguments, and
    these differences: the file is another one, which pyNei does not read;
    `num_vars_per_block` is an argument here, where in pyNei the size
    belongs to the ``Variants``; a source with no variants is written, where
    pyNei raises a ``ValueError``; the file of a call that failed is taken
    away; and every column is written, where pyNei writes the ones its
    chunks happen to carry.
    """
    if not isinstance(variants, Variants):
        # The path of the VCF where the variants of it go is the mistake
        # that is easiest to make, and what it gave was the `AttributeError`
        # of an object with no source inside it.
        raise TypeError(
            f"`variants` is {variants!r}, a {type(variants).__name__}, and "
            f"`write_vars` writes the variants of a source: give it what "
            f"`open_vcf` gives, write_vars(open_vcf(vcf_path), path)"
        )
    _core.write_vars(variants._source, path, num_vars_per_block)
