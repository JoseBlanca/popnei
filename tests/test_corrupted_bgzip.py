"""A bgzipped VCF whose bytes were damaged after bgzip wrote it.

A review of the VCF reader found one change of two bytes of
`tests/reference/vcf/many.vcf.gz` that made the whole file read as no
variant and no error: the length of the extra field of the second member,
`06 00`, read as `44 54`, which says that the member is longer than it is.
The owner decided on 21 September 2026 that such a file is an error however
improbable the damage, so the reader reads a bgzipped file by the size that
each member states and checks what comes out of it. bcftools 1.24 reads the
same file as no variant and says nothing, so popnei is stricter here.

A damaged file is an `OSError`, by the owner's convention: a file that
cannot be read, that was cut short or that is corrupted. The variants of the
members before the damaged one are given first.
"""

from pathlib import Path

import pytest
from popnei import open_vcf

# Where the length of the extra field of the second member is, and what
# bgzip wrote there: 6 bytes, the one field of a bgzip member. The member
# starts at the byte 310 and the file holds four of them.
_THE_TWO_BYTES = slice(320, 322)
_AS_BGZIP_WROTE_THEM = b"\x06\x00"
_AS_THE_REVIEW_CHANGED_THEM = b"\x44\x54"

# How many variants are in the members before the damaged one: the first
# member of `many.vcf.gz` holds the header and no whole data line.
_VARIANTS_BEFORE_THE_DAMAGED_MEMBER = 0


def _with_the_damaged_member(reference_vcf_dir: Path, tmp_path: Path) -> Path:
    """`many.vcf.gz` with those two bytes changed, in `tmp_path`."""
    whole = bytearray((reference_vcf_dir / "many.vcf.gz").read_bytes())
    assert bytes(whole[_THE_TWO_BYTES]) == _AS_BGZIP_WROTE_THEM
    whole[_THE_TWO_BYTES] = _AS_THE_REVIEW_CHANGED_THEM
    path = tmp_path / "damaged.vcf.gz"
    path.write_bytes(bytes(whole))
    return path


def test_a_bgzipped_vcf_whose_member_is_damaged_is_an_error_and_not_an_empty_file(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    path = _with_the_damaged_member(reference_vcf_dir, tmp_path)
    variants = open_vcf(path, only_passed=False)
    read = 0
    with pytest.raises(OSError) as refusal:
        for block in variants.iter_blocks(num_vars_per_block=100):
            read += block.num_vars
    assert read == _VARIANTS_BEFORE_THE_DAMAGED_MEMBER
    assert refusal.value.filename == str(path)
    assert refusal.value.errno is None
    message = str(refusal.value)
    # A user who gets this looks at the member with `xxd -s 310`, and then
    # fetches the file again.
    assert "corrupted" in message
    assert "member 2" in message
    assert "byte 310" in message
    assert "again" in message
