"""Reading a VCF."""

from pathlib import Path

from popnei import _core
from popnei.variant import Variants


def open_vcf(
    vcf_path: str | Path,
    ploidy: int = _core.DEFAULT_PLOIDY,
    only_passed: bool = _core.DEFAULT_ONLY_PASSED,
) -> Variants:
    """The variants of the VCF at `vcf_path`, plain or gzipped.

    It reads the header, so a file that is not a VCF, or whose header
    popnei cannot read, is a ``ValueError`` here and not at the first
    calculation, and a file that cannot be opened is an ``OSError`` that
    carries the path in ``filename``. The variants themselves are read again
    at every pass over what it returns.

    `ploidy` is how many alleles every genotype of the file holds, the same
    for every individual and every variant, and a genotype of any other
    number of alleles is a ``ValueError`` when it is read: popnei does not
    read a VCF of mixed ploidies, because its calculations are not defined
    for one.

    `only_passed` leaves out the variants that failed a filter, those whose
    FILTER column is neither ``PASS`` nor a dot; a dot says that no filter
    was applied. With it false every variant of the file is given, and
    nothing then says which ones had failed.

    It is pyNei's ``vars_from_vcf`` under another name, with the ploidy and
    the filter as arguments, which pyNei has not; pyNei takes the ploidy
    from the first genotype of the file and gives every variant.
    """
    return Variants(_core.open_vcf(vcf_path, ploidy, only_passed))
