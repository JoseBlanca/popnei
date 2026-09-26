"""The density of the variants along the chromosomes from Python.

`docs/specs/stats.md` has it under "The density of the variants along the
chromosomes", and what it is compared with is tabix 1.24: the count of each
window of 1000 base pairs of `tests/reference/vcf/many.vcf.gz`, which
`tests/reference/stats/make_reference.py` keeps in `many.density.tsv`.
`many.vcf` holds 500 variants of chr1 and chr2, the first of chr1 at 1000
and its last at 10213, those of chr2 at 10250 and 19463, and its
`##contig` lines have no length. pyNei has no density of the variants.
"""

from pathlib import Path

import numpy
import pytest
from popnei import VarDensity, calc_var_density, open_vars, open_vcf, write_vars

REFERENCE_DIR = Path(__file__).parent / "reference"
MANY = REFERENCE_DIR / "vcf" / "many.vcf"
WRITE = REFERENCE_DIR / "vcf" / "write.vcf"
DENSITY_OF_MANY = REFERENCE_DIR / "stats" / "many.density.tsv"


def _many():
    """The 500 variants of `many.vcf`, every one given, as tabix counts
    them."""
    return open_vcf(MANY, only_passed=False)


def _windows(density: VarDensity) -> list[tuple[str, int, int, int]]:
    """The rows of the frame, each as the chromosome, the start, the end and
    the count."""
    frame = density.windows
    return [
        (chrom, int(start), int(end), int(num_vars))
        for chrom, start, end, num_vars in zip(
            frame["chrom"], frame["start"], frame["end"], frame["num_vars"], strict=True
        )
    ]


def _laid_end_to_end(chrom: str, width: int, last_end: int, counts: list[int]):
    """The windows of `chrom` of `width` base pairs from the position 1, the
    last one ending at `last_end`, with `counts`."""
    windows = []
    for index, count in enumerate(counts):
        start = index * width + 1
        end = last_end if index == len(counts) - 1 else start + width - 1
        windows.append((chrom, start, end, count))
    return windows


def _a_vcf(tmp_path: Path, contigs: list[str], lines: list[tuple[str, int]]) -> Path:
    """A VCF of one individual with the `##contig` lines `contigs` and a
    variant at each chromosome and position of `lines`."""
    text = ["##fileformat=VCFv4.3", *contigs]
    text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta")
    text += [f"{chrom}\t{pos}\t.\tA\tT\t.\tPASS\t.\tGT\t0/1" for chrom, pos in lines]
    path = tmp_path / "density.vcf"
    path.write_text("\n".join(text) + "\n")
    return path


def test_var_density_of_many_vcf_in_windows_of_1000_is_what_tabix_counts():
    density = calc_var_density(_many(), 1000)
    rows = [line.split("\t") for line in DENSITY_OF_MANY.read_text().splitlines()[1:]]
    tabix = [
        (chrom, int(start), int(end), int(count)) for chrom, start, end, count in rows
    ]
    assert _windows(density) == tabix
    # The first table of "How it is verified", as the spec writes it.
    expected = _laid_end_to_end("chr1", 1000, 11000, [1] + [27] * 9 + [6])
    expected += _laid_end_to_end("chr2", 1000, 20000, [0] * 10 + [21] + [27] * 8 + [13])
    assert _windows(density) == expected
    assert list(density.windows.columns) == ["chrom", "start", "end", "num_vars"]
    assert density.windows["start"].dtype == numpy.uint64
    assert density.windows["end"].dtype == numpy.uint64
    assert density.windows["num_vars"].dtype == numpy.uint32
    assert density.pass_stats.num_vars == 500


def test_var_density_of_many_vcf_with_chrom_lengths_ends_each_chromosome_at_its_length():
    density = calc_var_density(
        _many(), 1000, chrom_lengths={"chr1": 12000, "chr2": 19500}
    )
    expected = _laid_end_to_end("chr1", 1000, 12000, [1] + [27] * 9 + [6, 0])
    expected += _laid_end_to_end("chr2", 1000, 19500, [0] * 10 + [21] + [27] * 8 + [13])
    assert _windows(density) == expected


def test_var_density_refuses_a_variant_past_a_length_of_chrom_lengths():
    with pytest.raises(ValueError) as refused:
        calc_var_density(_many(), 1000, chrom_lengths={"chr1": 10000, "chr2": 20000})
    message = str(refused.value)
    # The first variant of chr1 past 10000, where the pass stops, and the
    # file it is in, which every error of a file names first.
    assert message.startswith(str(MANY))
    for part in ("chr1", "10028", "10000", "`chrom_lengths`"):
        assert part in message, message


def test_var_density_of_the_worked_example_with_the_lengths_of_its_header():
    density = calc_var_density(open_vcf(WRITE), 600)
    assert _windows(density) == [
        ("chr1", 1, 600, 1),
        ("chr1", 601, 1200, 2),
        ("chr1", 1201, 1800, 0),
        ("chr1", 1801, 2000, 0),
        ("chr2", 1, 600, 1),
        ("chr2", 601, 1200, 0),
        ("chr2", 1201, 1500, 1),
    ]
    assert density.pass_stats.num_vars == 5


def test_var_density_of_a_vars_file_reads_the_lengths_it_keeps(tmp_path):
    path = tmp_path / "write.vars"
    write_vars(open_vcf(WRITE), path)
    density = calc_var_density(open_vars(path), 500)
    assert _windows(density) == [
        ("chr1", 1, 500, 1),
        ("chr1", 501, 1000, 1),
        ("chr1", 1001, 1500, 1),
        ("chr1", 1501, 2000, 0),
        ("chr2", 1, 500, 1),
        ("chr2", 501, 1000, 0),
        ("chr2", 1001, 1500, 1),
    ]


def test_var_density_with_chrom_lengths_reads_no_length_of_the_source():
    # chr1 is not named, so it has no length and comes after chr2, whose
    # windows go to the length given, past the 1500 of the header.
    density = calc_var_density(open_vcf(WRITE), 500, chrom_lengths={"chr2": 3000})
    assert _windows(density) == _laid_end_to_end(
        "chr2", 500, 3000, [1, 0, 1, 0, 0, 0]
    ) + _laid_end_to_end("chr1", 500, 1500, [1, 1, 1])


def test_var_density_counts_the_variants_the_steps_of_the_pass_keep():
    variants = _many()
    variants.filter_individuals(["ind00", "ind01"])
    variants.filter_by_missing_data(0.0)
    density = calc_var_density(variants, 1000)
    counts = density.windows["num_vars"]
    assert int(counts.sum()) == density.pass_stats.num_vars
    assert density.pass_stats.num_vars < 500
    # The `Variants` is as it was: the same call gives the same result.
    again = calc_var_density(variants, 1000)
    assert _windows(again) == _windows(density)


def test_var_density_of_a_window_size_of_0_is_refused_without_the_file():
    with pytest.raises(ValueError, match="`window_size`") as refused:
        calc_var_density(_many(), 0)
    assert not str(refused.value).startswith(str(MANY))


@pytest.mark.parametrize("window_size", [2.5, "1000", True, None])
def test_var_density_of_a_window_size_that_is_no_whole_number_is_a_type_error(
    window_size,
):
    with pytest.raises(TypeError, match="`window_size`"):
        calc_var_density(_many(), window_size)


def test_var_density_of_a_negative_window_size_is_refused():
    with pytest.raises(ValueError, match="`window_size`"):
        calc_var_density(_many(), -1)


def test_var_density_refuses_chrom_lengths_that_are_no_lengths():
    with pytest.raises(TypeError, match="`chrom_lengths`"):
        calc_var_density(_many(), 1000, chrom_lengths=[("chr1", 12000)])
    with pytest.raises(TypeError, match="`chrom_lengths`"):
        calc_var_density(_many(), 1000, chrom_lengths={1: 12000})
    with pytest.raises(TypeError, match="`chrom_lengths`"):
        calc_var_density(_many(), 1000, chrom_lengths={"chr1": 1.5e4})
    with pytest.raises(ValueError, match="`chrom_lengths`"):
        calc_var_density(_many(), 1000, chrom_lengths={"chr1": -5})
    with pytest.raises(ValueError, match="`chrom_lengths`") as refused:
        calc_var_density(_many(), 1000, chrom_lengths={"chr1": 12000, "chr2": 0})
    assert "chr2" in str(refused.value)
    assert not str(refused.value).startswith(str(MANY))


def test_var_density_refuses_a_variant_past_the_length_of_the_header(tmp_path):
    path = _a_vcf(
        tmp_path, ["##contig=<ID=chr1,length=1000>"], [("chr1", 1000), ("chr1", 1001)]
    )
    with pytest.raises(ValueError, match="the header of the source") as refused:
        calc_var_density(open_vcf(path), 300)
    assert "1001" in str(refused.value)


def test_var_density_refuses_a_variant_at_position_0(tmp_path):
    path = _a_vcf(tmp_path, [], [("chr1", 0)])
    with pytest.raises(ValueError, match="position 0"):
        calc_var_density(open_vcf(path), 300)


def test_var_density_of_more_windows_than_it_has_is_refused(tmp_path):
    path = _a_vcf(tmp_path, [], [("chr1", 10_000_001)])
    with pytest.raises(ValueError, match="10000001 windows"):
        calc_var_density(open_vcf(path), 1)
    with pytest.raises(ValueError, match="10000001 windows"):
        calc_var_density(open_vcf(path), 1, chrom_lengths={"chr1": 10_000_001})


def test_var_density_of_a_length_past_2_to_the_53_keeps_the_ends_exact(tmp_path):
    length = 2**53 + 1
    path = _a_vcf(tmp_path, [f"##contig=<ID=chr1,length={length}>"], [("chr1", length)])
    density = calc_var_density(open_vcf(path), 2**52)
    assert _windows(density) == [
        ("chr1", 1, 2**52, 0),
        ("chr1", 2**52 + 1, 2**53, 0),
        ("chr1", 2**53 + 1, 2**53 + 1, 1),
    ]


def test_var_density_of_a_pass_that_gives_no_variant_is_refused(tmp_path):
    path = _a_vcf(tmp_path, ["##contig=<ID=chr1,length=1000>"], [])
    with pytest.raises(ValueError, match="no variant"):
        calc_var_density(open_vcf(path), 100)


def test_var_density_of_something_that_is_not_variants_is_a_type_error():
    with pytest.raises(TypeError, match="open_vcf"):
        calc_var_density(str(MANY), 1000)
