"""The file a read failed on, where Python keeps it: ``OSError.filename``.

An `OSError` of the standard library carries the number the system gave,
which decides which exception it is, and the file in `filename`, which is
where a caller looks for it. A read of popnei can fail with no number from
the system: the decoder of gzip refuses a stream that ends in the middle,
and that is an error of Rust's own. The file is what a user needs either
way, so it is in `filename` either way.

The case is a VCF compressed with the gzip of Python and cut, which is not
what `many.vcf.gz` of `tests/reference/vcf/` is: that one bgzip wrote, and a
bgzipped file that was cut is refused with the error that says so, a
``ValueError``.
"""

import gzip
from pathlib import Path

import pytest
from popnei import open_vcf

# The header of a VCF of three individuals and enough data lines that the
# cut falls in the middle of the variants and not in the header.
_HEADER = (
    "##fileformat=VCFv4.4",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\tind3",
)
_NUM_VARS = 2000


def _gzipped_vcf_cut_in_the_middle(path: Path) -> Path:
    """A VCF compressed with gzip, of which the first half is written."""
    lines = (
        f"chr1\t{variant + 1}\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"
        for variant in range(_NUM_VARS)
    )
    whole = gzip.compress(("\n".join((*_HEADER, *lines)) + "\n").encode())
    path.write_bytes(whole[: len(whole) // 2])
    return path


def test_a_gzipped_vcf_that_was_cut_names_the_file_in_filename(tmp_path: Path) -> None:
    path = _gzipped_vcf_cut_in_the_middle(tmp_path / "cut.vcf.gz")
    variants = open_vcf(path)
    with pytest.raises(OSError) as refusal:
        list(variants.iter_blocks())
    # The decoder found the stream cut, and the file system refused nothing,
    # so Python has no number to build the exception of a cause with.
    assert refusal.value.errno is None
    assert refusal.value.filename == str(path)
    assert str(path) in str(refusal.value)
