"""The file a read went wrong on, in what the user is told about it.

A user who reads a directory of VCFs, one after another, is told which file
went wrong: every error that comes from reading one names it. Where the
exception is an `OSError` the file is in `filename`, which is where a caller
looks for it and where Python prints it, after the message; where it is a
`ValueError` or a `RuntimeError` the message starts with it. An argument
that popnei refuses, a field that is not one of the five or a block of no
variant, is no error of a file and names none.

An `OSError` of the standard library also carries the number the system
gave, which decides which exception it is. A read of popnei can fail with no
number: the decoder of gzip refuses a stream that ends in the middle, and
that is an error of Rust's own, not of the file system. The file is what a
user needs either way.

The gzipped file here is one that the gzip of Python wrote, which is not
what `many.vcf.gz` of `tests/reference/vcf/` is: that one bgzip wrote, and
the files that bgzip wrote are in `test_truncated_and_infinite.py`.
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


def test_the_message_of_a_wrong_data_line_starts_with_the_file(write_vcf) -> None:
    """The `x` in the POS column is a `ValueError`, a file whose content is
    not what a VCF holds, and the message names the file and then the line
    and the column, which is what the core says."""
    path = write_vcf(["chr1\tx\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"])
    variants = open_vcf(path)
    with pytest.raises(ValueError) as refusal:
        list(variants.iter_blocks())
    message = str(refusal.value)
    assert message.startswith(f"{path}: "), message
    # The header of the fixture is three lines, so the variant is in the
    # fourth line of the file.
    assert "line 4 of the VCF" in message
    assert "POS" in message


def test_the_message_of_bytes_that_are_not_a_vcf_starts_with_the_file(
    tmp_path: Path,
) -> None:
    """The header is read when the file is opened, so this is the message of
    `open_vcf` and not of a block."""
    path = tmp_path / "not_a_vcf.txt"
    path.write_text("chr1,100,A,T\n")
    with pytest.raises(ValueError) as refusal:
        open_vcf(path)
    message = str(refusal.value)
    assert message.startswith(f"{path}: "), message
    assert "not a VCF" in message


def test_an_argument_that_is_refused_names_no_file(reference_vcf_dir: Path) -> None:
    """What a user wrote is wrong wherever the file is, and the file has
    nothing to do with it: `fields`, `num_vars_per_block` and `ploidy` are
    refused by their own words."""
    path = reference_vcf_dir / "cases.vcf"
    variants = open_vcf(path)
    for asked_for in ({"fields": ("depth",)}, {"num_vars_per_block": 0}):
        with pytest.raises(ValueError) as refusal:
            list(variants.iter_blocks(**asked_for))
        message = str(refusal.value)
        assert str(reference_vcf_dir) not in message, message
    # The ploidy is refused when the file is opened, because the reader
    # needs it to read the first genotype, and it is what the user typed
    # all the same: 0 alleles in a genotype, and 256 above the 255 that an
    # allele of popnei counts.
    for ploidy in (0, 256):
        with pytest.raises(ValueError) as refusal:
            open_vcf(path, ploidy=ploidy)
        message = str(refusal.value)
        assert str(ploidy) in message, message
        assert str(reference_vcf_dir) not in message, message
