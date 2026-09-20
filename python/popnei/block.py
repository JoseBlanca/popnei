"""The block: a run of consecutive variants held as arrays.

A block is how the genotypes leave popnei. The calculations consume blocks
inside the Rust core, and a user who wants the genotypes for an analysis of
their own asks a :class:`popnei.Variants` for them with ``iter_blocks``,
which gives one block after another until the source is at its end.
"""

from dataclasses import dataclass
from typing import Literal

import numpy

Field = Literal["chrom", "pos", "id", "alleles", "qual"]
"""The name of a column of a block, which ``iter_blocks`` takes in `fields`.

The genotypes are not among them: every block holds them.
"""


# Two blocks are the same one or they are not: `eq=False` keeps the
# comparison of the arrays, which has no true or false for more than one
# genotype, out of `==` and leaves a block hashable by what it is.
@dataclass(frozen=True, eq=False)
class Block:
    """The variants of one block, each field a column of the block.

    A field other than the genotypes is there only when ``iter_blocks`` was
    asked for it, and ``None`` when it was not.
    """

    gts: numpy.ndarray
    """The genotypes, an int8 array of variants x individuals x ploidy.

    0 is the reference allele and 1 and above the alternative ones, in the
    order in which the source declares them, and -1 an allele that was not
    called. The array is the one the core filled, which reaches numpy
    without a copy, and it is read only, as ``pos`` and ``qual`` are: a
    calculation that writes works on ``numpy.array(block.gts)``, which is a
    copy of its own.
    """

    num_vars: int
    """How many variants the block holds."""

    chrom: tuple[str, ...] | None
    """The name of the chromosome of each variant."""

    pos: numpy.ndarray | None
    """The position of each variant, a uint64 array, 1 based as in a VCF."""

    id: tuple[str | None, ...] | None
    """The id of each variant, ``None`` for a variant that has none."""

    alleles: tuple[tuple[str, ...], ...] | None
    """The alleles of each variant, the reference one first.

    Each one is the text the source gave, so a symbolic allele, ``<DEL>``,
    and the allele of an overlapping deletion, ``*``, are alleles like any
    other.
    """

    qual: numpy.ndarray | None
    """The quality of each variant, a float32 array, phred scaled as the
    QUAL of a VCF: 30 is one chance in a thousand that there is no variant
    at that site. It is NaN for a variant whose source gives no quality.
    """


def _block_of(columns) -> Block:
    """The block of the columns that `popnei._core` gives for one block."""
    gts, chrom, pos, id_, alleles, qual = columns
    return Block(
        gts=gts,
        num_vars=gts.shape[0],
        chrom=chrom,
        pos=pos,
        id=id_,
        alleles=alleles,
        qual=qual,
    )
