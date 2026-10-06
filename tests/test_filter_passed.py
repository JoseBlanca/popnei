"""The filter of the variants that passed their FILTER from Python: which
variants it keeps, its counts, over a VCF and over a vars file, and what it
refuses.

`docs/specs/filters.md` has the filter under "The filter of the variants
that passed their FILTER". pyNei has no such filter, so the variants are
those of bcftools 1.24, `bcftools view -H -f .,PASS many.vcf`, run on 6
October 2026: the 475 of the 500 variants of `many.vcf` whose FILTER is
`PASS` or a dot. `many.bcftools.tsv` holds every one of the 500 with its
FILTER, and the test takes the 475 from there.
"""

from pathlib import Path

import pyarrow
import pyarrow.ipc
import pytest
from popnei import (
    FilteringStats,
    PassStats,
    Step,
    calc_pairwise_kosman_dists,
    open_vars,
    open_vcf,
    write_vars,
)
from popnei.variant import Variants

MANY_NUM_VARS = 500
NUM_PASSED = 475
# The first ten that passed, on chr1, as the spec gives them; the first that
# failed, at 1259, is between the sixth and the seventh.
FIRST_TEN_PASSED = [1000, 1037, 1074, 1111, 1148, 1185, 1222, 1296, 1333, 1370]
FIRST_FAILED = 1259
# The columns of `many.bcftools.tsv` before the genotypes, the last of which
# is FILTER.
FILTER_COLUMN = 6


def _many(reference_vcf_dir: Path) -> Variants:
    """The 500 variants of `many.vcf`, the 25 that failed their FILTER
    among them."""
    return open_vcf(reference_vcf_dir / "many.vcf", only_passed=False)


def _passed_of_bcftools(reference_vcf_dir: Path) -> list[tuple[str, int]]:
    """The chromosome and the position of the rows of `many.bcftools.tsv`
    whose FILTER is `PASS` or a dot."""
    rows = [
        line.split("\t")
        for line in (reference_vcf_dir / "many.bcftools.tsv").read_text().splitlines()
        if line
    ]
    return [
        (row[0], int(row[1])) for row in rows if row[FILTER_COLUMN] in ("PASS", ".")
    ]


def _kept(variants: Variants, num_vars_per_block: int | None = None):
    """The chromosome and the position of every variant of one pass, and
    the counts of that pass."""
    blocks = variants.iter_blocks(
        fields=("chrom", "pos"), num_vars_per_block=num_vars_per_block
    )
    kept = [
        (chrom, int(pos))
        for block in blocks
        for chrom, pos in zip(block.chrom, block.pos, strict=True)
    ]
    return kept, blocks.pass_stats


def _counts(vars_processed: int) -> PassStats:
    """The counts of a pass with this filter alone, given `vars_processed`
    variants, of which it keeps the 475."""
    return PassStats(
        num_vars=NUM_PASSED,
        filtering={
            "passed": FilteringStats(
                vars_processed=vars_processed, vars_kept=NUM_PASSED
            )
        },
        stopped_early=False,
    )


@pytest.mark.parametrize("num_vars_per_block", [7, None])
def test_filter_passed_keeps_the_475_of_bcftools_of_the_500_of_many_vcf(
    reference_vcf_dir: Path, num_vars_per_block: int | None
) -> None:
    """In blocks of 7, which puts a variant that failed in the first block
    beside six that passed, and in blocks of the default size."""
    variants = _many(reference_vcf_dir)
    variants.filter_passed()

    kept, pass_stats = _kept(variants, num_vars_per_block)

    assert kept[:10] == [("chr1", pos) for pos in FIRST_TEN_PASSED]
    assert ("chr1", FIRST_FAILED) not in kept
    assert kept == _passed_of_bcftools(reference_vcf_dir)
    assert pass_stats == _counts(MANY_NUM_VARS)


def test_filter_passed_keeps_every_variant_of_a_vcf_opened_with_only_passed(
    reference_vcf_dir: Path,
) -> None:
    """`open_vcf` by default gives the 475 alone, and the filter is given
    them all and keeps them all."""
    variants = open_vcf(reference_vcf_dir / "many.vcf")
    variants.filter_passed()

    kept, pass_stats = _kept(variants)

    assert kept == _passed_of_bcftools(reference_vcf_dir)
    assert pass_stats == _counts(NUM_PASSED)


def test_filter_passed_keeps_the_475_of_a_vars_file_written_from_the_500(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """The vars file holds whether each of the 500 passed, and the filter
    over it keeps what it keeps over the VCF."""
    path = tmp_path / "many.vars"
    write_vars(_many(reference_vcf_dir), path)
    variants = open_vars(path)
    variants.filter_passed()

    kept, pass_stats = _kept(variants)

    assert kept == _passed_of_bcftools(reference_vcf_dir)
    assert pass_stats == _counts(MANY_NUM_VARS)


def test_filter_passed_over_a_vars_file_of_1_1_is_a_value_error_with_its_path(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """A vars file that pyarrow rewrote from one of popnei's without the
    `passed` column and with the version 1.1, which is what a file written
    before format 1.2 holds. It opens, and the first calculation over it is
    refused, with the file in front, and not run as if every variant had
    passed."""
    written = tmp_path / "many.vars"
    write_vars(_many(reference_vcf_dir), written)
    without = tmp_path / "of_1_1.vars"
    with pyarrow.ipc.open_file(written) as reader:
        footer = reader.metadata
        schema_metadata = dict(reader.schema.metadata)
        table = reader.read_all()
    popnei_key = b"popnei"
    schema_metadata[popnei_key] = schema_metadata[popnei_key].replace(
        b'"format_version":"1.2"', b'"format_version":"1.1"'
    )
    assert b'"format_version":"1.1"' in schema_metadata[popnei_key]
    table = table.drop_columns(["passed"]).replace_schema_metadata(schema_metadata)
    with pyarrow.ipc.new_file(without, table.schema, metadata=footer) as writer:
        for batch in table.to_batches():
            writer.write_batch(batch)
    variants = open_vars(without)
    variants.filter_passed()

    with pytest.raises(ValueError) as refusal:
        calc_pairwise_kosman_dists(variants)

    message = str(refusal.value)
    assert message.startswith(str(without)), message
    assert "no record of whether they passed their FILTER" in message
    assert "from format 1.2" in message


def test_filter_passed_is_a_step_of_its_kind_and_no_argument(
    reference_vcf_dir: Path,
) -> None:
    variants = _many(reference_vcf_dir)

    assert variants.filter_passed() is None

    assert variants.steps == (Step(kind="passed", args={}),)


def test_filter_passed_refuses_a_second_filter_of_its_kind(
    reference_vcf_dir: Path,
) -> None:
    """The second would keep what the first kept. The steps are as they
    were."""
    variants = _many(reference_vcf_dir)
    variants.filter_passed()

    with pytest.raises(ValueError) as refusal:
        variants.filter_passed()

    assert "filtered by passed already" in str(refusal.value)
    assert [step.kind for step in variants.steps] == ["passed"]


def test_filter_passed_is_refused_after_the_filter_of_the_first_n(
    reference_vcf_dir: Path,
) -> None:
    """It takes variants out, so after the first n it would leave fewer than
    n. The steps are as they were."""
    variants = _many(reference_vcf_dir)
    variants.filter_first_n(10)

    with pytest.raises(ValueError) as refusal:
        variants.filter_passed()

    assert "a filter by passed after it" in str(refusal.value)
    assert "first_n already" in str(refusal.value)
    assert [step.kind for step in variants.steps] == ["first_n"]
