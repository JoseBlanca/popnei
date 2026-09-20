"""The arguments of popnei that say how many of something there are.

`ploidy` and `num_vars_per_block` count things, so what they take is a whole
number of 1 or more. A user who writes another number gets a ``ValueError``
that names what they gave, whatever the size of the number: an integer of
Python has no bound, and converting one into a number of Rust in the
signature raises an ``OverflowError`` that names neither the argument nor
what is wrong with it, so the binding crate converts it itself.

Which message a user gets depends on where the number is refused. A
negative one, and one above what a whole number of Rust holds, are refused
by the binding crate, which names the argument as the user writes it. A
number that converts and that popnei cannot work with, 0 and a block of more
variants than the machine has memory for, is refused by the core, whose
message says what is wrong with that number.

The reference VCFs are the ones of `docs/specs/io_vcf.md`. Nothing is read
here: every case is refused when the file is opened or at the call that asks
for the blocks.
"""

from pathlib import Path

import pytest
from popnei import open_vcf

# A number above what a whole number of Rust holds on a machine of 64 bits,
# and one far above it, which is an integer of Python of three words.
_ABOVE_A_WHOLE_NUMBER_OF_RUST = 2**64
_FAR_ABOVE_IT = 2**70

# What a user reads for each ploidy that counts no alleles: the number they
# wrote, and the words that say what is wrong with it.
_REFUSED_PLOIDIES = [
    (-1, ["ploidy", "-1"]),
    (0, ["ploidy", "0"]),
    (2**63, ["ploidy", str(2**63), "255"]),
    (_ABOVE_A_WHOLE_NUMBER_OF_RUST, ["ploidy", str(_ABOVE_A_WHOLE_NUMBER_OF_RUST)]),
    (_FAR_ABOVE_IT, ["ploidy", str(_FAR_ABOVE_IT)]),
]

# The same for the variants of a block. The two that the core refuses say
# what a block can hold and not which argument was written, because the core
# knows nothing of the names of the API.
_REFUSED_BLOCK_SIZES = [
    (-1, ["num_vars_per_block", "-1"]),
    (0, ["blocks of 0 variants", "1 variant at least"]),
    (2**63, [str(2**63), "fewer variants in a block"]),
    (
        _ABOVE_A_WHOLE_NUMBER_OF_RUST,
        ["num_vars_per_block", str(_ABOVE_A_WHOLE_NUMBER_OF_RUST)],
    ),
    (_FAR_ABOVE_IT, ["num_vars_per_block", str(_FAR_ABOVE_IT)]),
]


# The largest whole number this machine counts, which no message names: the
# largest a user may write is the one of their argument, 255 alleles in a
# genotype, and the core says that one. A message with both would give a
# user two limits for one argument.
_WHAT_THE_MACHINE_COUNTS = str(2**64 - 1)


@pytest.mark.parametrize(("value", "in_the_message"), _REFUSED_PLOIDIES)
def test_a_ploidy_that_counts_no_alleles_is_refused_with_what_was_written(
    reference_vcf_dir: Path, value: int, in_the_message: list[str]
) -> None:
    with pytest.raises(ValueError) as refusal:
        open_vcf(reference_vcf_dir / "cases.vcf", ploidy=value)
    message = str(refusal.value)
    for words in in_the_message:
        assert words in message, message
    assert _WHAT_THE_MACHINE_COUNTS not in message, message


@pytest.mark.parametrize(("value", "in_the_message"), _REFUSED_BLOCK_SIZES)
def test_a_block_of_a_number_of_variants_that_counts_nothing_is_refused(
    reference_vcf_dir: Path, value: int, in_the_message: list[str]
) -> None:
    variants = open_vcf(reference_vcf_dir / "cases.vcf")
    with pytest.raises(ValueError) as refusal:
        variants.iter_blocks(num_vars_per_block=value)
    message = str(refusal.value)
    for words in in_the_message:
        assert words in message, message
    assert _WHAT_THE_MACHINE_COUNTS not in message, message
