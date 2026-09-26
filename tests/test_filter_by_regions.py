"""The filter by regions from Python: which variants a BED file keeps, the
steps and the counts it gives, and what it refuses at the call.

`docs/specs/filters.md` has the filter under "The filter by regions". pyNei
has no filter by regions, so the comparison is with bcftools 1.24 and plink2
v2.0.0-a.7.7, which keep a variant by its position alone as this filter does:
`tests/reference/filters/make_reference.py` ran `bcftools view -T` and
`plink2 --extract bed0` with the BED of "How it is verified",
`tests/reference/filters/regions.bed`, on `many.vcf` of `docs/specs/io_vcf.md`,
500 variants of 50 diploid individuals, read with every variant given
because neither program honours FILTER, and stored the 45 variants they keep
in `regions.txt` and the 455 they keep with `-T ^` and `--exclude bed0` in
`excluded_regions.txt`, a chromosome and a position on each line.
"""

import gzip
from pathlib import Path

import pyarrow
import pyarrow.ipc
import pytest
from popnei import FilteringStats, Step, open_vars, open_vcf, write_vars
from popnei.variant import Variants

REFERENCE_FILTERS_DIR = Path(__file__).parent / "reference" / "filters"
REGIONS_BED = REFERENCE_FILTERS_DIR / "regions.bed"
WRITE_REGIONS_BED = REFERENCE_FILTERS_DIR / "write_regions.bed"

# The six regions of the BED of the spec join into five: chr1 1 to 2000,
# chr1 4991 to 5100, chr2 10250 alone, chr2 19001 to 30000 and chr3 1 to
# 100000.
NUM_REGIONS = 5


def _many(reference_vcf_dir: Path) -> Variants:
    """The 500 variants of `many.vcf`, the ones that failed their FILTER
    among them, which is what bcftools and plink2 read."""
    return open_vcf(reference_vcf_dir / "many.vcf", only_passed=False)


def _reference(name: str) -> list[tuple[str, int]]:
    """The variants bcftools and plink2 kept, by chromosome and position."""
    lines = (REFERENCE_FILTERS_DIR / f"{name}.txt").read_text().splitlines()
    return [(chrom, int(pos)) for chrom, pos in (line.split("\t") for line in lines)]


def _kept(variants: Variants, num_vars_per_block: int | None = None):
    """The chromosome and the position of every variant of one pass, and
    the counts of the filters of that pass."""
    blocks = variants.iter_blocks(
        fields=("chrom", "pos"), num_vars_per_block=num_vars_per_block
    )
    kept = [
        (chrom, int(pos))
        for block in blocks
        for chrom, pos in zip(block.chrom, block.pos, strict=True)
    ]
    return kept, blocks.pass_stats.filtering


@pytest.mark.parametrize("num_vars_per_block", [7, None])
def test_filter_by_regions_keeps_the_45_of_bcftools_and_plink2_and_counts_them(
    reference_vcf_dir: Path, num_vars_per_block: int | None
) -> None:
    """The 45 variants of `many.vcf` inside the regions: 28 of chr1 up to
    2000, the three at chr1 4996, 5033 and 5070, chr2 10250 and 13 of chr2
    from 19001. The filter was given the 500 of the file."""
    variants = _many(reference_vcf_dir)
    variants.filter_by_regions(REGIONS_BED)

    kept, filtering = _kept(variants, num_vars_per_block)

    assert kept == _reference("regions")
    assert len(kept) == 45
    assert len([pos for chrom, pos in kept if chrom == "chr1" and pos <= 2000]) == 28
    assert [pos for chrom, pos in kept if chrom == "chr1" and pos > 2000] == [
        4996,
        5033,
        5070,
    ]
    assert [pos for chrom, pos in kept if chrom == "chr2" and pos <= 19000] == [10250]
    assert filtering == {"regions": FilteringStats(vars_processed=500, vars_kept=45)}


@pytest.mark.parametrize("num_vars_per_block", [7, None])
def test_filter_by_regions_with_exclude_keeps_the_455_of_bcftools_and_plink2(
    reference_vcf_dir: Path, num_vars_per_block: int | None
) -> None:
    variants = _many(reference_vcf_dir)
    variants.filter_by_regions(REGIONS_BED, exclude=True)

    kept, filtering = _kept(variants, num_vars_per_block)

    assert kept == _reference("excluded_regions")
    assert len(kept) == 455
    assert filtering == {
        "excluded_regions": FilteringStats(vars_processed=500, vars_kept=455)
    }


def test_filter_by_regions_of_the_worked_example_puts_each_edge_on_its_side(
    reference_vcf_dir: Path,
) -> None:
    """The six lines of `write.vcf` with the regions of positions 100, 251,
    1000 and chr2 1: `chr1 99 100` is 100 alone, `chr1 999 1000` holds 1000
    and not 1001, and the deletion at chr1 250 is outside `chr1 250 251`."""
    for exclude, expected in [
        (False, [("chr1", 100), ("chr1", 1000), ("chr2", 1)]),
        (True, [("chr1", 250), ("chr1", 1001), ("chr2", 1500)]),
    ]:
        variants = open_vcf(reference_vcf_dir / "write.vcf", only_passed=False)
        variants.filter_by_regions(WRITE_REGIONS_BED, exclude=exclude)
        kept, _ = _kept(variants)
        assert kept == expected, exclude


def test_filter_by_regions_returns_none_and_adds_a_step_of_each_kind(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """The `steps` of both kinds, with the path as it was given and the
    number of regions once the six are joined into five. A `regions` and an
    `excluded_regions` step stand together, and the second counts over what
    the first kept: of the 45, the three of chr1 up to 1100 are excluded."""
    variants = _many(reference_vcf_dir)
    repeats = tmp_path / "chr1_start.bed"
    repeats.write_text("chr1\t0\t1100\n")

    assert variants.filter_by_regions(str(REGIONS_BED)) is None
    assert variants.steps == (
        Step(
            kind="regions",
            args={"bed_path": str(REGIONS_BED), "num_regions": NUM_REGIONS},
        ),
    )

    variants.filter_by_regions(repeats, exclude=True)
    assert variants.steps == (
        Step(
            kind="regions",
            args={"bed_path": str(REGIONS_BED), "num_regions": NUM_REGIONS},
        ),
        Step(
            kind="excluded_regions",
            args={"bed_path": str(repeats), "num_regions": 1},
        ),
    )
    assert "excluded_regions(bed_path=" in repr(variants)

    kept, filtering = _kept(variants)
    assert len(kept) == 42
    assert filtering == {
        "regions": FilteringStats(vars_processed=500, vars_kept=45),
        "excluded_regions": FilteringStats(vars_processed=45, vars_kept=42),
    }


def test_filter_by_regions_after_a_threshold_filter_counts_what_that_filter_kept(
    reference_vcf_dir: Path,
) -> None:
    """A filter before it hands it what it kept, so the filter by regions
    is given fewer than 500, and a threshold filter that keeps every variant
    gives it all 500."""
    variants = _many(reference_vcf_dir)
    variants.filter_by_maf(1.0)
    variants.filter_by_regions(REGIONS_BED)

    kept, filtering = _kept(variants)

    assert kept == _reference("regions")
    assert filtering == {
        "maf": FilteringStats(vars_processed=500, vars_kept=500),
        "regions": FilteringStats(vars_processed=500, vars_kept=45),
    }


def test_a_second_filter_by_regions_of_a_kind_is_refused_with_the_kind(
    reference_vcf_dir: Path,
) -> None:
    variants = _many(reference_vcf_dir)
    variants.filter_by_regions(REGIONS_BED)
    before = variants.steps

    with pytest.raises(ValueError) as refusal:
        variants.filter_by_regions(WRITE_REGIONS_BED)

    message = str(refusal.value)
    assert "filtered by regions already" in message
    # What the user wrote is wrong whatever file is read: no path in front.
    assert not message.startswith(str(WRITE_REGIONS_BED))
    assert variants.steps == before

    variants.filter_by_regions(WRITE_REGIONS_BED, exclude=True)
    with pytest.raises(ValueError, match="excluded_regions already"):
        variants.filter_by_regions(REGIONS_BED, exclude=True)


@pytest.mark.parametrize(
    ("text", "line", "what"),
    [
        # Fewer than three columns separated by tabs, after two lines that
        # are skipped and counted.
        ("track name=x\n# a comment\nchr1\t0\n", 3, "2 columns separated by tabs"),
        ("chr1 0 10\n", 1, "BED separates by tabs"),
        # A start or an end that is not a whole number of 0 or more.
        ("chr1\t0\t10\nchr1\t-5\t10\n", 2, "its start is `-5`"),
        ("chr1\t0\tten\n", 1, "its end is `ten`"),
        # A start that is not below its end.
        ("chr1\t5\t5\n", 1, "its start, 5, is not below its end, 5"),
    ],
)
def test_a_wrong_line_of_a_bed_is_a_value_error_with_the_file_and_the_line(
    reference_vcf_dir: Path, tmp_path: Path, text: str, line: int, what: str
) -> None:
    bed = tmp_path / "wrong.bed"
    bed.write_text(text)
    variants = _many(reference_vcf_dir)

    with pytest.raises(ValueError) as refusal:
        variants.filter_by_regions(bed)

    message = str(refusal.value)
    assert message.startswith(f"{bed}: line {line} of the BED file: "), message
    assert what in message
    assert variants.steps == ()


def test_a_bed_with_no_region_is_a_value_error_with_the_file(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    bed = tmp_path / "empty.bed"
    bed.write_text("track name=x\n# nothing\n\n")

    with pytest.raises(ValueError) as refusal:
        _many(reference_vcf_dir).filter_by_regions(bed)

    assert str(refusal.value).startswith(f"{bed}: the BED file holds no region")


def test_a_bed_that_is_not_there_is_a_file_not_found_error_with_its_path(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    bed = tmp_path / "missing.bed"

    with pytest.raises(FileNotFoundError) as refusal:
        _many(reference_vcf_dir).filter_by_regions(bed)

    assert refusal.value.filename == str(bed)


def test_a_gzipped_bed_cut_short_is_an_os_error_with_its_path(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    bed = tmp_path / "cut.bed.gz"
    bed.write_bytes(gzip.compress(REGIONS_BED.read_bytes())[:-10])

    with pytest.raises(OSError) as refusal:
        _many(reference_vcf_dir).filter_by_regions(bed)

    assert refusal.value.filename == str(bed)


def test_a_gzipped_bed_keeps_the_variants_of_the_plain_one(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """Read through gzip because it starts with the bytes of gzip, whatever
    its name."""
    bed = tmp_path / "regions.bed"
    bed.write_bytes(gzip.compress(REGIONS_BED.read_bytes()))
    variants = _many(reference_vcf_dir)
    variants.filter_by_regions(bed)

    kept, _ = _kept(variants)

    assert kept == _reference("regions")
    assert variants.steps[0].args["num_regions"] == NUM_REGIONS


def test_an_exclude_that_is_not_a_bool_is_a_type_error(
    reference_vcf_dir: Path,
) -> None:
    variants = _many(reference_vcf_dir)

    with pytest.raises(TypeError, match="`exclude` is True or False"):
        variants.filter_by_regions(REGIONS_BED, exclude=1)

    assert variants.steps == ()


def test_filter_by_regions_over_a_source_with_no_positions_is_a_value_error(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """A vars file without the `chrom` and the `pos` columns, which pyarrow
    writes from one of popnei's, gives the error of a field the pass needs,
    with the file in front, at the first block of the pass."""
    written = tmp_path / "many.vars"
    write_vars(_many(reference_vcf_dir), written)
    without = tmp_path / "without_positions.vars"
    with pyarrow.ipc.open_file(written) as reader:
        footer = reader.metadata
        schema_metadata = reader.schema.metadata
        table = reader.read_all()
    table = table.drop_columns(["chrom", "pos"]).replace_schema_metadata(
        schema_metadata
    )
    with pyarrow.ipc.new_file(without, table.schema, metadata=footer) as writer:
        for batch in table.to_batches():
            writer.write_batch(batch)
    variants = open_vars(without)
    variants.filter_by_regions(REGIONS_BED)

    with pytest.raises(ValueError) as refusal:
        list(variants.iter_blocks(fields=("qual",)))

    message = str(refusal.value)
    assert message.startswith(str(without)), message
    assert "`chrom and pos`" in message
