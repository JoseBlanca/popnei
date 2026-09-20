"""`iter_blocks`: the genotypes as arrays, and pyNei beside them.

The comparison with pyNei is the one `docs/specs/block.md` asks for under
"How it is verified", and it is also the comparison of the VCF reader of
`docs/specs/io_vcf.md`: popnei and pyNei read the same four files and every
field of every variant has to agree. popnei is read with
`only_passed=False`, because pyNei gives every variant whatever its FILTER
is, and the files are diploid, because pyNei reports the ploidy 2 for every
file, which is its issue 19. pyNei has `pandas.NA` where popnei has `None`
for a variant with no id and NaN for one with no quality, and the two count
as equal here.

The rest of the tests are of what the blocks are made of: which columns a
block carries, where the cuts between blocks fall, and the array of
genotypes that reaches numpy without a copy.
"""

from pathlib import Path

import numpy
import pandas
import pytest
from popnei import open_vcf
from pynei import vars_from_vcf

# Every field a block can carry besides the genotypes.
ALL_FIELDS = ("chrom", "pos", "id", "alleles", "qual")

# The four files both libraries read. pyNei cannot read
# `differences.vcf`: it refuses a genotype that starts with a separator.
FILES_BOTH_LIBRARIES_READ = ("cases.vcf", "cases.vcf.gz", "many.vcf", "many.vcf.gz")


def _popnei_columns(variants, num_vars_per_block):
    """The columns of every block of `variants`, one after another."""
    blocks = list(
        variants.iter_blocks(fields=ALL_FIELDS, num_vars_per_block=num_vars_per_block)
    )
    return {
        "gts": numpy.concatenate([block.gts for block in blocks]),
        "chrom": tuple(name for block in blocks for name in block.chrom),
        "pos": tuple(int(pos) for block in blocks for pos in block.pos),
        "id": tuple(id_ for block in blocks for id_ in block.id),
        "alleles": tuple(alleles for block in blocks for alleles in block.alleles),
        "qual": numpy.concatenate([block.qual for block in blocks]),
    }


def _pynei_columns(variants):
    """The same columns from the chunks of pyNei, with its missing values.

    A value that pyNei leaves as `pandas.NA` is `None` for an id and NaN for
    a quality, which is what popnei gives.
    """
    chunks = list(variants.iter_vars_chunks())
    vars_info = [chunk.vars_info for chunk in chunks]
    return {
        "gts": numpy.concatenate([chunk.gts.gt_values for chunk in chunks]),
        "chrom": tuple(str(name) for info in vars_info for name in info["chrom"]),
        "pos": tuple(int(pos) for info in vars_info for pos in info["pos"]),
        "id": tuple(
            None if pandas.isna(id_) else str(id_)
            for info in vars_info
            for id_ in info["id"]
        ),
        "alleles": tuple(
            tuple(alleles) for chunk in chunks for alleles in chunk.alleles
        ),
        "qual": numpy.array(
            [
                numpy.nan if pandas.isna(qual) else float(qual)
                for info in vars_info
                for qual in info["qual"]
            ],
            dtype=numpy.float32,
        ),
    }


@pytest.mark.parametrize("file_name", FILES_BOTH_LIBRARIES_READ)
@pytest.mark.parametrize(
    "num_vars_per_block", [7, None], ids=["blocks of 7", "default"]
)
def test_the_blocks_of_popnei_hold_what_the_chunks_of_pynei_hold(
    reference_vcf_dir: Path, file_name: str, num_vars_per_block: int | None
) -> None:
    path = reference_vcf_dir / file_name
    popnei_variants = open_vcf(path, only_passed=False)
    pynei_variants = vars_from_vcf(path)

    assert popnei_variants.individuals == pynei_variants.samples
    assert popnei_variants.num_individuals == pynei_variants.num_samples
    assert popnei_variants.ploidy == pynei_variants.ploidy

    ours = _popnei_columns(popnei_variants, num_vars_per_block)
    theirs = _pynei_columns(pynei_variants)
    numpy.testing.assert_array_equal(ours["gts"], theirs["gts"])
    assert ours["chrom"] == theirs["chrom"]
    assert ours["pos"] == theirs["pos"]
    assert ours["id"] == theirs["id"]
    assert ours["alleles"] == theirs["alleles"]
    # assert_array_equal takes two NaNs at the same place as equal, which is
    # what a variant with no quality gives in both libraries.
    numpy.testing.assert_array_equal(ours["qual"], theirs["qual"])


def test_a_block_carries_the_chromosome_and_the_position_and_no_more(
    reference_vcf_dir: Path,
) -> None:
    """What a user who asks for nothing gets."""
    variants = open_vcf(reference_vcf_dir / "cases.vcf")
    block = next(iter(variants.iter_blocks()))
    assert block.chrom == ("chr1", "chr1", "chr1")
    assert block.pos.tolist() == [100, 300, 400]
    assert block.id is None
    assert block.alleles is None
    assert block.qual is None


def test_the_chromosome_and_the_position_travel_together(
    reference_vcf_dir: Path,
) -> None:
    """They are one field of the core, so asking for one fills both."""
    variants = open_vcf(reference_vcf_dir / "cases.vcf")
    block = next(iter(variants.iter_blocks(fields=["pos"])))
    assert block.chrom == ("chr1", "chr1", "chr1")
    assert block.pos.tolist() == [100, 300, 400]


def test_the_genotypes_alone_leave_every_other_column_out(
    reference_vcf_dir: Path,
) -> None:
    variants = open_vcf(reference_vcf_dir / "cases.vcf")
    block = next(iter(variants.iter_blocks(fields=())))
    assert block.num_vars == 3
    assert block.chrom is None
    assert block.pos is None
    assert block.id is None
    assert block.alleles is None
    assert block.qual is None


def test_the_genotypes_are_an_int8_array_of_variants_individuals_and_ploidy(
    reference_vcf_dir: Path,
) -> None:
    """The array of the core, handed to numpy and never written to again.

    It owns no data of its own: what it holds is the allocation the core
    filled, which is what `docs/specs/block.md` asks of the binding crate.
    """
    variants = open_vcf(reference_vcf_dir / "many.vcf")
    block = next(iter(variants.iter_blocks(num_vars_per_block=100)))
    assert block.gts.dtype == numpy.int8
    assert block.gts.shape == (100, 50, 2)
    assert not block.gts.flags.writeable
    assert not block.gts.flags.owndata
    with pytest.raises(ValueError, match="read-only"):
        block.gts[0, 0, 0] = 1


def test_a_block_names_the_chromosomes_of_its_own_variants(write_vcf) -> None:
    """The names come from the table of the reader, which grows as it reads.

    A block holds a few of the chromosomes of a source, which in a de novo
    assembly has 10000 scaffolds or more, and each of its variants has to
    get the name of its own.
    """
    path = write_vcf(
        [
            f"chr{number}\t{10 * number}\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"
            for number in (1, 2, 3, 4, 1)
        ]
    )
    blocks = list(open_vcf(path).iter_blocks(num_vars_per_block=2))
    assert [block.chrom for block in blocks] == [
        ("chr1", "chr2"),
        ("chr3", "chr4"),
        ("chr1",),
    ]


def test_a_block_tells_what_it_holds_without_printing_it(
    reference_vcf_dir: Path,
) -> None:
    """What a user sees in a session or in a traceback.

    A block of 10000 variants holds as many chromosomes, ids and alleles,
    which printed are hundreds of kilobytes of one line.
    """
    block = next(
        iter(
            open_vcf(reference_vcf_dir / "many.vcf").iter_blocks(
                fields=ALL_FIELDS, num_vars_per_block=100
            )
        )
    )
    printed = repr(block)
    assert len(printed) < 200, printed
    assert "100 variants" in printed
    assert "(100, 50, 2)" in printed
    for name in ALL_FIELDS:
        assert name in printed
    assert "chr1" not in printed


def _writable_arrays_under(array: numpy.ndarray) -> list[str]:
    """The arrays that share the memory of `array` and can be written into.

    An array made by reshaping or slicing another one keeps it as its
    ``base``, and writing into that base writes into both.
    """
    writable = []
    under: object = array
    while isinstance(under, numpy.ndarray):
        if under.flags.writeable:
            writable.append(f"{under.shape} {under.dtype}")
        under = under.base
    return writable


def test_no_array_a_user_reaches_from_a_block_can_be_written_into(
    reference_vcf_dir: Path,
) -> None:
    """A block is frozen, and its arrays hold the memory of the core.

    A view that was left writable is a way around both: the array under
    `gts` shares its genotypes, so writing into that one changes what the
    block says the source held.
    """
    variants = open_vcf(reference_vcf_dir / "many.vcf")
    block = next(iter(variants.iter_blocks(fields=ALL_FIELDS, num_vars_per_block=100)))
    for name in ("gts", "pos", "qual"):
        assert _writable_arrays_under(getattr(block, name)) == [], name
    with pytest.raises(ValueError, match="read-only"):
        block.gts[0, 0, 0] = 1


def test_the_blocks_are_cut_by_the_count_of_the_variants(
    reference_vcf_dir: Path,
) -> None:
    """The last block of a source is the only one that can be shorter.

    `many.vcf` gives 475 variants by default, and a chromosome that ends in
    the middle of a block does not end the block: `chr1` ends inside the
    third one.
    """
    variants = open_vcf(reference_vcf_dir / "many.vcf")
    blocks = list(variants.iter_blocks(num_vars_per_block=100))
    assert [block.num_vars for block in blocks] == [100, 100, 100, 100, 75]
    assert [block.gts.shape[0] for block in blocks] == [100, 100, 100, 100, 75]


def test_every_pass_over_the_variants_reads_the_file_again(
    reference_vcf_dir: Path,
) -> None:
    """A `Variants` can be given to any number of calculations."""
    variants = open_vcf(reference_vcf_dir / "cases.vcf")
    first = _popnei_columns(variants, None)
    second = _popnei_columns(variants, 2)
    numpy.testing.assert_array_equal(first["gts"], second["gts"])
    assert first["pos"] == second["pos"]


def test_a_field_that_is_not_one_of_the_five_is_refused(
    reference_vcf_dir: Path,
) -> None:
    """The refusal comes at the call and not at the first block.

    The message is the core's, which both languages share, and it lists
    the five names a user can write.
    """
    variants = open_vcf(reference_vcf_dir / "cases.vcf")
    with pytest.raises(ValueError, match="depth") as refusal:
        variants.iter_blocks(fields=("chrom", "depth"))
    for name in ALL_FIELDS:
        assert f"`{name}`" in str(refusal.value)


def test_one_name_where_a_sequence_of_names_goes_is_refused(
    reference_vcf_dir: Path,
) -> None:
    """A string is a sequence of its letters, and `a` is not a field."""
    variants = open_vcf(reference_vcf_dir / "cases.vcf")
    with pytest.raises(TypeError, match="sequence") as refusal:
        variants.iter_blocks(fields="alleles")
    assert 'fields=("alleles",)' in str(refusal.value)


def test_blocks_of_no_variant_are_refused(reference_vcf_dir: Path) -> None:
    """A block holds one variant at least, and `None` asks for the default."""
    variants = open_vcf(reference_vcf_dir / "cases.vcf")
    with pytest.raises(ValueError, match="0 variants"):
        variants.iter_blocks(num_vars_per_block=0)


def test_a_vcf_with_no_variant_gives_no_block(write_vcf) -> None:
    variants = open_vcf(write_vcf([]))
    assert list(variants.iter_blocks()) == []
    assert variants.individuals == ("ind1", "ind2", "ind3")
