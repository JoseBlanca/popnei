"""What the tests of the VCF reader and of the blocks share.

The VCFs of `tests/reference/vcf/`, which
`tests/reference/vcf/make_reference.py` writes and `docs/specs/io_vcf.md`
describes, and a writer of small VCFs for the cases that no reference file
holds, an error of a data line among them.
"""

from collections.abc import Callable, Iterable
from pathlib import Path

import pytest

REFERENCE_VCF_DIR = Path(__file__).parent / "reference" / "vcf"

# The `##` lines and the `#CHROM` line of a VCF of three individuals, which
# the data lines a test writes go under.
_HEADER = (
    "##fileformat=VCFv4.4",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\tind3",
)


@pytest.fixture
def reference_vcf_dir() -> Path:
    """The directory of the VCFs that both libraries are run on."""
    return REFERENCE_VCF_DIR


@pytest.fixture
def write_vcf(tmp_path: Path) -> Callable[[Iterable[str]], Path]:
    """A function that writes a VCF of three individuals and gives its path.

    It takes the data lines, each with its columns separated by tabs, and
    puts them under a header with the `GT` format and the individuals
    `ind1`, `ind2` and `ind3`.
    """

    def write(data_lines: Iterable[str]) -> Path:
        path = tmp_path / "written.vcf"
        lines = list(_HEADER) + list(data_lines)
        path.write_text("\n".join(lines) + "\n")
        return path

    return write
