"""Reading and writing a VCF."""

from dataclasses import dataclass
from pathlib import Path

from popnei import _core
from popnei.variant import PassStats, Variants, _pass_stats_of


@dataclass(frozen=True)
class VcfWritten:
    """What :func:`popnei.write_vcf` gives back: the counts of its pass."""

    pass_stats: PassStats
    """How many variants were written, and how many each filter of the
    ``Variants`` was given and kept."""


def open_vcf(
    vcf_path: str | Path,
    ploidy: int | None = None,
    only_passed: bool = _core.DEFAULT_ONLY_PASSED,
) -> Variants:
    """The variants of the VCF at `vcf_path`, plain or gzipped.

    It reads the header, so a file that is not a VCF, or whose header
    popnei cannot read, is a ``ValueError`` here and not at the first
    calculation, and a file that cannot be opened is an ``OSError`` that
    carries the path in ``filename``. The variants themselves are read again
    at every pass over what it returns.

    A file that was cut short, and one that bgzip wrote and whose bytes were
    damaged, are found where the variants are read: the header of such a
    file reads, and :meth:`Variants.iter_blocks` raises the ``OSError`` of a
    file that popnei cannot read to its end, with the path in ``filename``
    and no ``errno``, after the variants it could give. When `ploidy` is
    left out, the damage can be found here instead, if it comes before the
    first genotype the ploidy is read from.

    `ploidy` is how many alleles every genotype of the file holds, the same
    for every individual and every variant, and a genotype of any other
    number of alleles is a ``ValueError`` when it is read: popnei does not
    read a VCF of mixed ploidies, because its calculations are not defined
    for one. It is 1 or more and at most 255, and a ploidy outside that is a
    ``ValueError`` at this call, before anything is read.

    When `ploidy` is ``None``, the default, it is read from the file here:
    it is the number of alleles of the first genotype that is not a single
    dot, which is a missing genotype of any ploidy, so ``./.`` says 2. The
    data lines are looked at in their order, those that failed their filter
    among them, up to 4096 of them. The call is a ``ValueError`` whose
    message starts with the path in three cases: none of those lines holds
    a genotype with alleles, and the ploidy then has to be given; the file
    has a header and no data line, which with a ploidy given is opened and
    gives no variants; and the first genotype with alleles holds more than
    255. A file whose
    genotypes are all missing is refused and not read as diploid because a
    filter by missing data and the statistics count the missing alleles of
    a missing genotype, one for each allele of the ploidy.

    `only_passed` leaves out the variants that failed a filter, those whose
    FILTER column is neither ``PASS`` nor a dot; a dot says that no filter
    was applied. With it false every variant of the file is given, and
    nothing then says which ones had failed.

    It is pyNei's ``vars_from_vcf`` under another name, with the ploidy and
    the filter as arguments, which pyNei has not. pyNei gives every variant,
    and it takes the ploidy from the genotype of the first individual in the
    first data line, a single dot among them, where popnei looks on to the
    first genotype with alleles.
    """
    return Variants(_core.open_vcf(vcf_path, ploidy, only_passed))


def write_vcf(variants: Variants, path: str | Path) -> VcfWritten:
    """Every variant of `variants`, after its steps, into a VCF at `path`.

    It is how the variants popnei kept, filtered by missing data, by
    individual or by any other step, reach plink2, bcftools or a program of
    the user's own. The file is compressed with bgzip, which is what tabix
    indexes and what bcftools asks a region of, when the path ends in
    ``.gz``, in any case of its two letters, ``a.VCF.GZ`` among them, and it
    is plain text otherwise.

    When the source is a VCF, each variant is written as its line was,
    every column of it, INFO, FILTER, the phase and the values of each
    individual other than GT among them, and the header is the source's
    with a ``#CHROM`` line of the individuals that were kept, in the order
    the filter of individuals named them. When that filter took individuals
    out, AC and AN, counts over individuals that are no longer in the file,
    are taken out of INFO and their ``##INFO`` lines out of the header, as
    ``bcftools annotate -x INFO/AC,INFO/AN`` does; a filter that keeps
    every individual, in any order, leaves them. The other values of INFO,
    a depth or a frequency, are written as the source had them.

    When the source is a vars file, the lines hold what the file holds:
    the chromosome, the position, the id, the alleles, the quality and the
    genotypes, with FILTER and INFO a dot, FORMAT ``GT`` and the alleles of
    each genotype joined by ``/``, since a vars file keeps no phase. The
    header has one ``##contig`` line for each chromosome whose length the
    vars file keeps.

    The lines are written in the order the source gives them, so a source
    that is not sorted gives a file that tabix refuses to index. A VCF read
    with `only_passed` false and written with no step is the same file,
    byte for byte, when its lines end in ``\n`` and none is empty. A
    source with no variants gives the header alone.

    What it gives back is a :class:`VcfWritten` with the counts of the pass
    it made: how many variants were written, and how many each filter of
    the `variants` was given and kept. The call reads the source once.

    The path is handled as :func:`popnei.write_vars` handles its own: a
    path that a file is already at is a ``ValueError`` and nothing is
    written; a path that no file can be made at, a directory or a path in a
    directory that is not there, is the ``OSError`` the file system gives
    for it, with the path in ``filename``; an error of the source names the
    source, and a file that could not be written is an ``OSError`` that
    carries its own path. The file of a call that failed, a Ctrl-C among
    the causes, is taken away, and when it cannot be, the exception carries
    a note that says so.

    pyNei has no VCF writer, so nothing is mirrored; the name follows
    :func:`popnei.write_vars`.
    """
    if not isinstance(variants, Variants):
        raise TypeError(
            f"`variants` is {variants!r}, a {type(variants).__name__}, and "
            f"`write_vcf` writes the variants of a source: give it what "
            f"`open_vcf` or `open_vars` gives, "
            f"write_vcf(open_vcf(vcf_path), path)"
        )
    counts = _core.write_vcf(variants._source, path, variants._steps)
    return VcfWritten(pass_stats=_pass_stats_of(counts))
