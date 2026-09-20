"""The two files popnei refuses and pyNei reads: a quality that is not
finite, and a bgzipped VCF that does not end with the mark of its end.

The owner decided both on 20 September 2026, and `docs/specs/io_vcf.md` has
them under "The cases a reader of the rules would not guess". A quality of
`nan` is an error because NaN is what a block holds for a variant with no
quality, so pyNei reads such a file as a variant that has none, and `inf`
and `1e400` as an infinite one. A file that bgzip wrote ends with an empty
gzip block of 28 bytes; without it the file was cut short, and pyNei gives
the variants before the cut and says nothing, where bcftools 1.24 says "no
BGZF EOF marker; file may be truncated".

Both reach Python as a `ValueError`, the exception of a file whose content
popnei cannot read. The variants that were read before the cut are given
first, so each test here says how many blocks came out before the error.
"""

from pathlib import Path

import pytest
from popnei import open_vcf

# The bytes of `many.vcf.gz`, which bgzip wrote as four members that end at
# the bytes 310, 12336, 21876 and 21904, with 0, 280, 500 and 500 whole data
# lines before each of them. `docs/specs/io_vcf.md` gives the two cuts:
# without the last 28 bytes, which are the mark of the end, and after the
# second member, which is a download that stopped.
_WHOLE = 21904
_MARK_OF_THE_END = 28
_AFTER_THE_SECOND_MEMBER = 12336

# What the two cuts hold: every variant of the file, and the 280 variants of
# the two first members, of which a Python user sees 200. Every pass of
# Python ends in `reblock`, which cuts and joins the blocks of its source to
# the size that was asked for, and `docs/specs/block.md` says that an error
# loses the block it happened in and, with it, the variants that `reblock`
# was keeping for its next block: here the 80 that were left over after two
# blocks of 100.
_VARIANTS_OF_THE_WHOLE_FILE = 500
_VARIANTS_BEFORE_THE_SECOND_CUT = 200


def _cut(reference_vcf_dir: Path, tmp_path: Path, bytes_of_it: int) -> Path:
    """`many.vcf.gz` cut to its first `bytes_of_it` bytes, in `tmp_path`."""
    whole = (reference_vcf_dir / "many.vcf.gz").read_bytes()
    assert len(whole) == _WHOLE
    path = tmp_path / f"cut_to_{bytes_of_it}.vcf.gz"
    path.write_bytes(whole[:bytes_of_it])
    return path


def _variants_before_the_error(path: Path) -> tuple[int, str]:
    """How many variants a file gives before it fails, and the message."""
    variants = open_vcf(path, only_passed=False)
    read = 0
    with pytest.raises(ValueError) as refusal:
        for block in variants.iter_blocks(num_vars_per_block=100):
            read += block.num_vars
    return read, str(refusal.value)


def test_a_bgzipped_vcf_without_the_mark_of_its_end_is_refused_after_its_variants(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """The file is whole but for its last 28 bytes, so every variant is
    there and nothing but the mark of the end says that it was cut."""
    path = _cut(reference_vcf_dir, tmp_path, _WHOLE - _MARK_OF_THE_END)
    read, message = _variants_before_the_error(path)
    assert read == _VARIANTS_OF_THE_WHOLE_FILE
    # A user who gets this has to fetch the file again, so the message says
    # that it is cut short and not only that a mark is missing.
    assert "bgzip" in message
    assert "cut short" in message
    assert "again" in message


def test_a_bgzipped_vcf_cut_where_a_member_ends_gives_its_variants_and_then_fails(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """A download that stopped between two blocks of bgzip: the gzip of
    what arrived is whole, and 220 of the 500 variants are not in it. The
    reader gives the 280 that are, in blocks of 100, 100 and 80, and the 80
    are lost with the error in the `reblock` of the pass."""
    path = _cut(reference_vcf_dir, tmp_path, _AFTER_THE_SECOND_MEMBER)
    read, message = _variants_before_the_error(path)
    assert read == _VARIANTS_BEFORE_THE_SECOND_CUT
    assert "cut short" in message


def test_a_bgzipped_vcf_cut_inside_a_member_gives_its_variants_and_then_fails(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """A download that stopped inside a block of bgzip, which is where it
    stops nearly always: the decoder runs out of bytes there, and what a
    user is told is that the file is cut short and not that a deflate
    stream is incomplete. The 21000 bytes that arrived hold 480 whole data
    lines, of which `reblock` gives the four blocks of 100 it filled."""
    path = _cut(reference_vcf_dir, tmp_path, 21000)
    read, message = _variants_before_the_error(path)
    assert read == 400
    assert "cut short" in message
    assert "deflate" not in message


def test_a_gzipped_vcf_that_bgzip_did_not_write_is_read_to_its_end(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """Plain gzip has no mark of its end, so there is none to miss."""
    import gzip

    plain = (reference_vcf_dir / "many.vcf").read_bytes()
    path = tmp_path / "plain_gzip.vcf.gz"
    path.write_bytes(gzip.compress(plain))

    variants = open_vcf(path, only_passed=False)
    read = sum(block.num_vars for block in variants.iter_blocks())
    assert read == _VARIANTS_OF_THE_WHOLE_FILE


@pytest.mark.parametrize("quality", ["nan", "inf", "-inf", "1e400", "1e39"])
def test_a_quality_that_is_not_finite_is_refused(write_vcf, quality: str) -> None:
    """`1e400` is above what a float of 64 bits holds and `1e39` above what
    the 32 bits of the column of a block hold: both read as infinite."""
    path = write_vcf([f"chr1\t10\t.\tA\tT\t{quality}\tPASS\t.\tGT\t0/0\t0/1\t1/1"])
    variants = open_vcf(path)
    with pytest.raises(ValueError, match="QUAL") as refusal:
        list(variants.iter_blocks(fields=("qual",)))
    message = str(refusal.value)
    assert quality in message
    # The header of the fixture is three lines, so the variant is in the
    # fourth line of the file.
    assert "line 4" in message


def test_a_quality_that_is_not_finite_is_read_when_the_quality_is_not_asked_for(
    write_vcf,
) -> None:
    """A column that is not parsed is not checked, which "How it runs" of
    `docs/specs/io_vcf.md` decides, and the default `fields` of
    `iter_blocks` does not ask for the quality."""
    path = write_vcf(["chr1\t10\t.\tA\tT\tnan\tPASS\t.\tGT\t0/0\t0/1\t1/1"])
    variants = open_vcf(path)
    blocks = list(variants.iter_blocks())
    assert [block.num_vars for block in blocks] == [1]
    assert blocks[0].pos.tolist() == [10]
