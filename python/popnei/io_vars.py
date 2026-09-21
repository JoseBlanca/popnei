"""Reading and writing a vars file, the file popnei keeps its variants in."""

from pathlib import Path

from popnei import _core
from popnei.variant import Variants


def open_vars(path: str | Path) -> Variants:
    """The variants of the vars file at `path`.

    A vars file is one arrow IPC file, also called feather v2, which
    :func:`popnei.write_vars` writes from any source of variants. It is
    where a user keeps their variants once the VCF has been read, so that
    the text is parsed once and every later pass reads a file of arrays.

    What it gives is the handle :func:`popnei.open_vcf` gives: the names of
    the individuals and the ploidy, which it reads from the file, and the
    variants through :meth:`Variants.iter_blocks`. Only the columns that a
    pass asks for are decompressed, and the variants themselves are read
    again at every pass.

    It reads the schema of the file and its footer, so a file that is not a
    vars file, one of a format version popnei does not read, one whose
    columns are not those of a vars file and one whose individuals name
    nobody, whose genotypes would hold no allele, are a ``ValueError`` here
    and not at the first calculation, with the path at the start of the
    message. A file that cannot be opened, and one that was cut short after
    it was written, are an ``OSError`` that carries the path in
    ``filename``.

    What is in the batches is read block by block, and refused there. A
    batch that popnei cannot read, of a file damaged after it was written,
    is an ``OSError`` with the path in ``filename`` and no ``errno``,
    because nothing of the file system refused anything. A quality that is
    a value and is not a finite number is a ``ValueError`` that names the
    variant: NaN in that column is how popnei says that a variant has no
    quality, and an infinite quality is a probability of no variant of 0.
    And a file whose buffers are compressed with zstd is a ``ValueError``
    at its first block, because arrow decompresses a batch when it reads it
    and no build of popnei carries the code that reads zstd; popnei writes
    lz4 and reads lz4 and no compression.

    It is pyNei's ``load_vars`` under another name, because the call opens
    the file and reads no variant, and these differences: the file is
    another one, so a vars file of pyNei is refused as a file without the
    key that says what popnei wrote; and pyNei's
    `desired_num_vars_per_chunk` is gone, because how many variants come out
    at a time is `num_vars_per_block` of :meth:`Variants.iter_blocks`, which
    does not have to be the number of variants a batch of the file holds.
    """
    return Variants(_core.open_vars(path))


def write_vars(
    variants: Variants, path: str | Path, num_vars_per_block: int | None = None
) -> None:
    """Every variant of `variants` into a vars file at `path`.

    A vars file is one arrow IPC file, also called feather v2, which pandas,
    R and polars open as a table with no popnei installed. It is where a
    user keeps their variants once the VCF has been read, so that the text
    is parsed once and every later pass reads a file of arrays.

    The source is a VCF or a vars file, whichever :class:`Variants` holds,
    so a file read with :func:`popnei.open_vars` is written again with
    another size of batch.

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

    Every error names the file it is about, which matters when the source
    is a vars file too and there are two of them. What popnei cannot read
    in the source, a line of a VCF among them, is a ``ValueError`` whose
    message starts with the path of that source, and a batch of a vars file
    that it cannot read is an ``OSError`` that carries that same path; the
    file being written, when it could not be written, a disc that filled up
    among the causes, is an ``OSError`` that carries its own path in
    ``filename``. The file that was being written is taken away
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
            f"`open_vcf` or `open_vars` gives, "
            f"write_vars(open_vcf(vcf_path), path)"
        )
    _core.write_vars(variants._source, path, num_vars_per_block)
