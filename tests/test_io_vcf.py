"""`open_vcf`: what a Python user reads from a VCF.

The cases are the ones `docs/specs/io_vcf.md` gives to pytest under "How it
is verified": the counts of `many.vcf`, the two variants of
`differences.vcf`, and the errors that reach Python as an exception. The
numbers are the literals of that spec, which come from bcftools 1.24. The
comparison with pyNei is in `test_block.py`, because Python sees what the
reader read through the blocks.
"""

from pathlib import Path

import numpy
import pytest
from popnei import open_vcf

# Every field a block can carry besides the genotypes.
ALL_FIELDS = ("chrom", "pos", "id", "alleles", "qual")

# An allele that was not called, as `docs/specs/variant.md` fixes it. The
# Python package publishes no name for it: a user reads it from a block.
MISSING_ALLELE = -1


def _joined(variants, fields=ALL_FIELDS):
    """The columns of every block of `variants`, one after another."""
    blocks = list(variants.iter_blocks(fields=fields))
    return {
        "gts": numpy.concatenate([block.gts for block in blocks]),
        "chrom": tuple(name for block in blocks for name in block.chrom),
        "pos": tuple(int(pos) for block in blocks for pos in block.pos),
        "id": tuple(id_ for block in blocks for id_ in block.id),
        "alleles": tuple(alleles for block in blocks for alleles in block.alleles),
        "qual": tuple(float(qual) for block in blocks for qual in block.qual),
    }


def test_the_individuals_and_the_ploidy_are_known_when_the_vcf_is_opened(
    reference_vcf_dir: Path,
) -> None:
    variants = open_vcf(reference_vcf_dir / "cases.vcf")
    assert variants.individuals == ("ind1", "ind2", "ind3")
    assert variants.num_individuals == 3
    assert variants.ploidy == 2


def test_the_ploidy_of_the_reader_is_the_one_that_was_asked_for(
    write_vcf,
) -> None:
    """A tetraploid VCF read with `ploidy` 4, which is not an error."""
    path = write_vcf(["chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0/1/1\t0/1/1/1\t0/0/0/0"])
    variants = open_vcf(path, ploidy=4)
    assert variants.ploidy == 4
    blocks = list(variants.iter_blocks())
    assert blocks[0].gts.shape == (1, 3, 4)
    assert list(blocks[0].gts[0, 0]) == [0, 0, 1, 1]


@pytest.mark.parametrize(
    ("only_passed", "counts"),
    [
        (
            True,
            {
                "variants": 475,
                "in_chr2": 238,
                "with_two_alternative_alleles": 53,
                "missing_genotypes": 1431,
                "half_called_genotypes": 240,
                "missing_alleles": 2622,
                "called_alleles": 44878,
                "sum_of_the_called_alleles": 24831,
            },
        ),
        (
            False,
            {
                "variants": 500,
                "in_chr2": 250,
                "with_two_alternative_alleles": 54,
                "missing_genotypes": 1511,
                "half_called_genotypes": 257,
                "missing_alleles": 2765,
                "called_alleles": 47235,
                "sum_of_the_called_alleles": 25954,
            },
        ),
    ],
    ids=["by default", "every variant"],
)
def test_the_counts_of_many_vcf(
    reference_vcf_dir: Path, only_passed: bool, counts: dict[str, int]
) -> None:
    """The sixteen numbers of the table of `docs/specs/io_vcf.md`.

    By default the reader gives the 475 variants whose FILTER is `PASS` or a
    dot, and with `only_passed=False` the 500 of the file.
    """
    variants = open_vcf(reference_vcf_dir / "many.vcf", only_passed=only_passed)
    columns = _joined(variants)
    gts = columns["gts"]
    missing = gts == MISSING_ALLELE
    # A genotype is missing when one of its alleles was not called, and it
    # is half called when another one was.
    missing_genotypes = missing.any(axis=2)
    half_called = missing_genotypes & ~missing.all(axis=2)
    assert gts.shape[0] == counts["variants"]
    assert columns["chrom"].count("chr2") == counts["in_chr2"]
    assert (
        sum(1 for alleles in columns["alleles"] if len(alleles) == 3)
        == counts["with_two_alternative_alleles"]
    )
    assert int(missing_genotypes.sum()) == counts["missing_genotypes"]
    assert int(half_called.sum()) == counts["half_called_genotypes"]
    assert int(missing.sum()) == counts["missing_alleles"]
    assert int((~missing).sum()) == counts["called_alleles"]
    assert int(gts[~missing].sum()) == counts["sum_of_the_called_alleles"]


def test_the_gzipped_many_vcf_gives_what_the_plain_one_gives(
    reference_vcf_dir: Path,
) -> None:
    """A decoder that stopped at the first gzip member would give no variant.

    bgzip wrote `many.vcf.gz` as four gzip members, and the first one is
    exactly the header of the VCF.
    """
    plain = _joined(open_vcf(reference_vcf_dir / "many.vcf"))
    gzipped = _joined(open_vcf(reference_vcf_dir / "many.vcf.gz"))
    numpy.testing.assert_array_equal(gzipped["gts"], plain["gts"])
    assert gzipped["pos"] == plain["pos"]
    assert gzipped["chrom"] == plain["chrom"]


def test_the_two_variants_of_differences_vcf(reference_vcf_dir: Path) -> None:
    """The table of `differences.vcf` of `docs/specs/io_vcf.md`.

    A genotype that starts with a separator, `/0/1`, is read as `0/1`, a
    genotype written as a single dot is missing in every allele, and `<DEL>`
    and `*` are alleles like any other. pyNei reads none of it.
    """
    variants = open_vcf(reference_vcf_dir / "differences.vcf")
    columns = _joined(variants)
    assert columns["chrom"] == ("chr2", "chr2")
    assert columns["pos"] == (50, 60)
    assert columns["id"] == ("ms1", None)
    assert columns["alleles"] == (("GTC", "G", "GTCT"), ("A", "<DEL>", "*"))
    assert columns["qual"][0] == 50.0
    assert numpy.isnan(columns["qual"][1])
    numpy.testing.assert_array_equal(
        columns["gts"],
        numpy.array(
            [[[0, 1], [0, 2], [-1, -1]], [[0, 1], [2, 2], [0, 0]]], dtype=numpy.int8
        ),
    )


def test_a_file_that_is_not_a_vcf_is_refused_when_it_is_opened(
    tmp_path: Path,
) -> None:
    """`open_vcf` reads the header, so nothing is left for the first block."""
    path = tmp_path / "not_a_vcf.txt"
    path.write_text("chrom\tpos\n")
    with pytest.raises(ValueError, match="not a VCF"):
        open_vcf(path)


def test_a_genotype_of_another_ploidy_is_refused_when_the_blocks_are_asked_for(
    write_vcf,
) -> None:
    """The line is read when the blocks are, and it names the individual."""
    path = write_vcf(
        [
            "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
            "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0/1/1\t1/1",
        ]
    )
    variants = open_vcf(path)
    with pytest.raises(ValueError, match="ind2") as refusal:
        list(variants.iter_blocks())
    # The header of the fixture is three lines, so the wrong genotype is in
    # the fifth line of the file.
    assert "line 5" in str(refusal.value)


def test_an_allele_that_the_variant_does_not_declare_is_refused(write_vcf) -> None:
    """pyNei reads such an allele number; popnei refuses it."""
    path = write_vcf(["chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/2"])
    variants = open_vcf(path)
    with pytest.raises(ValueError, match="ind3"):
        list(variants.iter_blocks())


def test_a_path_that_no_file_is_at_gives_an_oserror_that_carries_the_path(
    tmp_path: Path,
) -> None:
    """A user who misspells a path gets the exception pyNei users expect.

    `FileNotFoundError` derives from `OSError`, and the path is in
    `filename`, where the standard library puts it.
    """
    path = tmp_path / "there_is_no_such_vcf.vcf"
    with pytest.raises(OSError) as refusal:
        open_vcf(path)
    assert isinstance(refusal.value, FileNotFoundError)
    assert refusal.value.filename == str(path)


def test_a_count_that_is_negative_is_refused_by_its_name(
    reference_vcf_dir: Path,
) -> None:
    """A number of things that is below zero, which no argument takes.

    It is a `ValueError` like the ploidy of 0 beside it, and not the
    `OverflowError` of a conversion, and it names the argument and what was
    given for it.
    """
    with pytest.raises(ValueError, match="ploidy") as refusal:
        open_vcf(reference_vcf_dir / "cases.vcf", ploidy=-1)
    assert "-1" in str(refusal.value)

    variants = open_vcf(reference_vcf_dir / "cases.vcf")
    with pytest.raises(ValueError, match="num_vars_per_block") as refusal:
        variants.iter_blocks(num_vars_per_block=-5)
    assert "-5" in str(refusal.value)


def test_a_ploidy_of_zero_is_refused(reference_vcf_dir: Path) -> None:
    """The one thing `open_vcf` refuses that does not come from the file."""
    with pytest.raises(ValueError, match="ploidy"):
        open_vcf(reference_vcf_dir / "cases.vcf", ploidy=0)
