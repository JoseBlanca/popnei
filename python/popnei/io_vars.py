"""Reading and writing a vars file, the file popnei keeps its variants in."""

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy

from popnei import _core
from popnei.variant import PassStats, Variants, _pass_stats_of


@dataclass(frozen=True)
class VarsWritten:
    """What :func:`popnei.write_vars` and
    :func:`popnei.write_vars_from_segments` give back: the counts of their
    pass."""

    pass_stats: PassStats
    """How many variants were written, and how many each filter of the
    ``Variants`` was given and kept."""


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
    variants: Variants,
    path: str | Path,
    num_vars_per_block: int | None = None,
    allow_unsorted: bool = False,
) -> VarsWritten:
    """Every variant of `variants` into a vars file at `path`.

    A vars file is one arrow IPC file, also called feather v2, which pandas,
    R and polars open as a table with no popnei installed. It is where a
    user keeps their variants once the VCF has been read, so that the text
    is parsed once and every later pass reads a file of arrays.

    The source is a VCF or a vars file, whichever :class:`Variants` holds,
    so a file read with :func:`popnei.open_vars` is written again with
    another size of batch.

    The variants that are written are the ones the steps of `variants`
    keep: a ``Variants`` with a filter on it writes the variants that
    passed it and not the file as it was read, and the counts below say how
    many each filter was given and kept. A user who wants every variant of
    the source in the file writes from a ``Variants`` with no step, which
    :attr:`Variants.steps` and its ``repr`` show.

    The call reads the whole source once. The file holds the six columns of
    a VCF, the chromosome, the position, the id, the alleles, the quality
    and the genotypes, and whether each variant passed its FILTER, whether
    or not the user will read them, so that it can stand in for the VCF in
    any later analysis; a source that has no alleles to give, an array of
    genotypes, gives a file without that column, and one that does not
    keep whether its variants passed, a vars file written before that
    column existed, gives one without it.

    `num_vars_per_block` is how many variants a batch of the file holds, the
    last one aside, and ``None`` asks for the size popnei chooses for the
    number of individuals of the source, which is the size of its blocks.

    The variants have to be sorted unless `allow_unsorted` is ``True``: the
    positions of a chromosome do not go down from one variant to the next,
    and the variants of a chromosome come together, so a chromosome whose
    variants have ended does not come back. Two variants at one position
    are in order, and the chromosomes can come in any order. A variant out
    of that order is a ``ValueError`` that names it and the one before it,
    ``the 2nd variant written, chr1 100, comes after the 1st, chr1 5000``,
    counted among the variants written, after any filter, and the file is
    taken away as for any other error. With ``allow_unsorted=True`` the
    variants are written in the order they come. A source with no
    chromosomes and positions, a :class:`Variants` built from an array of
    genotypes, is not checked.

    What it gives back is a :class:`VarsWritten` with the counts of the pass
    it made: how many variants were written, and how many each filter of the
    `variants` was given and kept.

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
    away; unsorted variants are refused unless `allow_unsorted` is
    ``True``, where pyNei writes them; and every column is written, where pyNei writes the ones its
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
    counts = _core.write_vars(
        variants._source, path, num_vars_per_block, allow_unsorted, variants._steps
    )
    return VarsWritten(pass_stats=_pass_stats_of(counts))


# The five types an extra field can have, as numpy names them, from the
# core, which holds the one list of them.
_FIELD_TYPES = tuple(_core.extra_field_types())


def write_vars_from_segments(
    segments: Iterable[Mapping[str, Any]],
    path: str | Path,
    individuals: Sequence[str],
    ploidy: int,
    var_fields: Mapping[str, str | numpy.dtype] = {},
    gt_fields: Mapping[str, str | numpy.dtype] = {},
    chrom_lengths: Mapping[str, int] | None = None,
    num_vars_per_block: int | None = None,
    allow_unsorted: bool = False,
) -> VarsWritten:
    """A vars file at `path` written from the segments that `segments`
    yields, variants that the caller's code computes with no VCF in between.

    A segment is a run of consecutive variants, a region of a genome, given
    as a dict of arrays. The code that yields them can walk a genome region
    by region without holding the whole of it in memory: the call keeps the
    segment it was last given and the variants that wait for the next
    batch of the file, `num_vars_per_block` at most. The file is opened
    with :func:`popnei.open_vars`, and every filter and calculation of
    popnei runs on it as on any other.

    `individuals`, `ploidy`, `chrom_lengths` and the two declarations of the
    fields go into the schema of the file, which is written before the
    first segment is taken. Each segment is a dict with these keys:

    - ``chrom``: the chromosome of each variant, or one ``str``, the
      chromosome of every variant of the segment;
    - ``pos``: integers, the position of each variant, 1 based; the number
      of variants of the segment is the length of ``pos``, and every other
      value has one for each variant;
    - ``alleles``: for each variant a sequence of its alleles, each a text
      that is not empty, the reference allele first, at most 128, so the
      haplotypes of a region can be ``"H0"``, ``"H1"`` and so on;
    - ``gts``: a numpy array of any integer type, variants x individuals x
      ploidy, each value the number of an allele of its variant from 0, or
      -1 for a missing allele;
    - ``id``, optional: a text or ``None`` for each variant, and ``None``
      for the segment is the key left out;
    - ``qual``, optional: floats, NaN for a variant with no quality, and
      ``None`` for the segment is the key left out;
    - each declared field, required: an array of one value for each variant
      for a field of `var_fields`, and of variants x individuals for one of
      `gt_fields`.

    A segment of no variants is taken and writes nothing. A segment can
    hold any number of variants: the batches of the file hold
    `num_vars_per_block` variants each, the last aside, and ``None`` asks
    for the size popnei chooses for the number of individuals.

    An extra field is a named array of numbers that the variants carry
    through the filters and into any vars file written from this one, and
    that :meth:`Variants.iter_blocks` gives back when it is asked for by
    name. A field of `var_fields` has one value for each variant, as an
    INFO field of a VCF; one of `gt_fields` has one for each individual at
    each variant, as a FORMAT field. Each is declared by its name and its
    type, any spelling numpy takes of ``int32``, ``int64``, ``uint64``,
    ``float32`` or ``float64``. A name starts with a letter or ``_`` and
    holds letters, digits, ``_`` and ``.`` after it, is not the name of a
    column of the file, ``chrom``, ``pos``, ``id``, ``alleles``, ``qual``,
    ``gts`` and ``passed``, and is declared once. A missing value is the
    caller's to choose, NaN or a number, since popnei does not read the
    values. A field of integers takes integers of any integer type, a list
    of them too, and refuses floats, even whole ones, and a value its type
    does not hold, a negative one in a ``uint64`` field or 2**31 in an
    ``int32`` one, which numpy would wrap into another number. A field of
    floats takes floats and rounds them to its type, so a ``float32`` field
    holds the nearest ``float32`` of each value, and refuses a finite float
    beyond the largest of its type, 1e40 in a ``float32`` one. It takes an
    integer only when its type holds it exactly, so 2**24 + 1 in a
    ``float32`` field and 2**53 + 1 in a ``float64`` one are refused;
    ``qual`` follows the same rule. Every key of numbers refuses truth
    values, a numpy masked array, whose mask would be lost, for which
    ``.filled(...)`` gives the array to write, and a single value where it
    has one for each variant. An array in the byte order of another machine
    is taken.

    The variants have to be sorted unless `allow_unsorted` is ``True``: the
    positions of a chromosome do not go down from one variant to the next,
    and the variants of a chromosome come together, so a chromosome whose
    variants have ended does not come back. Two variants at one position
    are in order, and the chromosomes can come in any order. The check runs
    over the whole file, across segments, and a variant out of that order
    is a ``ValueError`` that names it and the one before it, counted among
    the variants written. With ``allow_unsorted=True`` the variants are
    written in the order they come.

    What it gives back is a :class:`VarsWritten` whose ``pass_stats`` has
    the number of variants written and an empty ``filtering``, since no
    filter ran. The ``Variants`` of the file has ``keeps_passed`` false:
    the file has no column of whether each variant passed a FILTER.

    The arguments are checked before the file is made, and each that is
    wrong is a ``ValueError``: no individual or two of one name, the name of
    an individual or of a chromosome of `chrom_lengths` that is empty or
    holds a tab, a newline, a carriage return or a space, a ploidy of 0 or
    above 255, a field of another type than the five, ``None`` among them,
    a name that breaks the rule above or is declared twice, a length of a
    chromosome that is not above 0, and a `num_vars_per_block` of 0. A
    path that a file is already at is a ``ValueError`` and nothing is
    written. `individuals` that is not a sequence, a set among them, whose
    order is its own, is a ``TypeError``, and a numpy array of names is
    taken; so is `chrom_lengths` that is not a mapping, and `segments` that
    is one segment.

    A segment that is wrong is a ``ValueError`` that names it and the
    variant of it, both counted from 1, with its chromosome and its
    position when those are not what is wrong: a key that is neither one of
    the six above nor a declared field, or a required one missing; a value
    whose length or shape does not agree with the variants, the individuals
    and the ploidy; positions that are not integers or are negative, and
    genotypes that are not integers; an allele of ``gts`` below -1 or not
    among the alleles of its variant; a variant with no allele or more than
    128, an allele that is an empty text or ``.``, and one given twice in a
    variant; a quality that is an infinity; a value a key of numbers does
    not take, by the rules above; an empty chromosome; and a tab, a
    newline, a carriage return or a space in a chromosome, an id or an
    allele, a comma in an allele and a semicolon in an id, which would break
    a VCF written from the file. An item that is not a mapping is a
    ``TypeError``.

    When anything fails after the file was made, the file is taken away, so
    a failed call leaves nothing at the path and the same call can be made
    again. An exception raised by `segments`, the caller's code, comes out
    of the call as it was raised, with its traceback, after the file is
    taken away. A Ctrl-C is raised when the next segment is asked for, and
    takes the file away too.

    pyNei writes a vars file from code of its user through a ``Variants``
    over an object whose ``iter_vars_chunks`` yields chunks. Here the call
    takes dicts and writes them at once, the fields are declared, there are
    genotype fields, the variants are checked to be sorted, a source of no
    variants is written, and the file of a call that failed is taken away,
    where pyNei leaves it.
    """
    if isinstance(segments, Mapping):
        raise TypeError(
            "`segments` is one segment, a mapping, and it is an iterable of "
            "segments: write [segment] to write that one"
        )
    if isinstance(individuals, str):
        raise TypeError(
            f"`individuals` is the text {individuals!r}, and it is a sequence "
            f"of the names of the individuals"
        )
    if not isinstance(individuals, Sequence | numpy.ndarray):
        # A set gives its names in an order of its own, which would become
        # the order of the individuals of the file with nothing said.
        raise TypeError(
            f"`individuals` is a {type(individuals).__name__}, and it is a "
            f"sequence of the names of the individuals in the order of the "
            f"file, a list or a tuple"
        )
    if chrom_lengths is not None and not isinstance(chrom_lengths, Mapping):
        raise TypeError(
            f"`chrom_lengths` is a {type(chrom_lengths).__name__}, and it is a "
            f"mapping of the name of each chromosome to its length, a dict"
        )
    counts = _core.write_vars_from_segments(
        segments,
        path,
        list(individuals),
        ploidy,
        _declared("var_fields", var_fields),
        _declared("gt_fields", gt_fields),
        None if chrom_lengths is None else list(chrom_lengths.items()),
        num_vars_per_block,
        allow_unsorted,
    )
    return VarsWritten(pass_stats=_pass_stats_of(counts))


def _declared(
    argument: str, fields: Mapping[str, str | numpy.dtype]
) -> list[tuple[str, str]]:
    """The name of each field of `fields` and the name numpy gives its
    type, in their order, for the argument `argument`.

    Raises a ``ValueError`` for a type that numpy does not read or that is
    not one of the five of a field.
    """
    declared = []
    for name, spelling in fields.items():
        # numpy reads `None` as float64, which nobody who wrote it asked for.
        try:
            type_name = None if spelling is None else numpy.dtype(spelling).name
        except TypeError, ValueError:
            type_name = None
        if type_name not in _FIELD_TYPES:
            raise ValueError(
                f"`{argument}` gives the field `{name}` the type {spelling!r}, "
                f"and a field is of one of {', '.join(_FIELD_TYPES)}"
            )
        declared.append((name, type_name))
    return declared
