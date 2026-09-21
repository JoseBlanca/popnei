"""`write_vars`: the vars file a Python user writes from their variants.

A vars file is one arrow IPC file, and the reference outside the project is
pyarrow, the arrow implementation that Apache Arrow publishes, which opens
what popnei wrote as any other program would. The cases are the ones
`docs/specs/io_vars.md` gives to pytest under "How it is verified" of the
writer: `many.vcf` written and read back, a source with no variants, a VCF
that fails half way, and a path that a file is already at. The numbers are
the literals of that spec, and what the columns hold is compared with
`many.bcftools.tsv`, what bcftools 1.24 read in the same VCF. pyNei is not
run here: it reads another file.
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pyarrow
import pyarrow.ipc
import pytest
from popnei import _core, open_vcf, write_vars

# What the two keys of a vars file are called, and the value of the key of
# the schema a reader of another version looks at first.
POPNEI_KEY = b"popnei"
POPNEI_BATCHES_KEY = b"popnei_batches"
FORMAT_VERSION = "1.0"

# An allele that was not called, which bcftools prints as a dot.
MISSING_ALLELE = -1
MISSING_VALUE = "."

# How many columns `bcftools query` printed before the genotypes: CHROM,
# POS, ID, REF, ALT, QUAL and FILTER, the seven of the format that
# `tests/reference/vcf/make_reference.py` gave it.
COLUMNS_BEFORE_THE_GENOTYPES = 7

# The 50 individuals of `many.vcf`, whose genotypes are 100 alleles in each
# row of the `gts` column, and the size of the blocks popnei chooses for a
# source of few individuals, which it puts in the `popnei` key when the
# caller asks for no size.
MANY_NUM_INDIVIDUALS = 50
MANY_ALLELES_PER_VAR = 100
LARGEST_NUM_VARS_PER_BLOCK = 10_000


def _columns_of(alleles_per_var: int) -> list[tuple[str, pyarrow.DataType, bool]]:
    """The columns of a vars file written from a VCF: the name, the arrow
    type and whether it takes nulls, as the table of "What it holds" of
    `docs/specs/io_vars.md` gives them."""
    return [
        ("chrom", pyarrow.string(), False),
        ("pos", pyarrow.uint64(), False),
        ("id", pyarrow.string(), True),
        ("alleles", pyarrow.list_(pyarrow.field("item", pyarrow.string())), False),
        ("qual", pyarrow.float32(), True),
        (
            "gts",
            pyarrow.list_(pyarrow.field("item", pyarrow.int8()), alleles_per_var),
            False,
        ),
    ]


# The regions of the five batches of `many.vcf` written with 100 variants in
# each: for every chromosome with a variant in the batch, in the order in
# which they first appear, the smallest and the largest position of its
# variants there. They are the literals of "How it is verified" of the
# writer, worked out there from `many.bcftools.tsv` by taking its rows 100
# at a time.
REGIONS_OF_EVERY_VARIANT = [
    [("chr1", 1000, 4663)],
    [("chr1", 4700, 8363)],
    [("chr1", 8400, 10213), ("chr2", 10250, 12063)],
    [("chr2", 12100, 15763)],
    [("chr2", 15800, 19463)],
]
# The same file without the 25 variants that failed a filter, which is what
# `open_vcf` gives by default: 475 variants, so the last batch is short.
REGIONS_OF_THE_PASSED = [
    [("chr1", 1000, 4848)],
    [("chr1", 4885, 8770)],
    [("chr1", 8807, 10213), ("chr2", 10250, 12655)],
    [("chr2", 12692, 16540)],
    [("chr2", 16577, 19463)],
]
# The 500 variants in one batch, which is what the size popnei chooses for
# 50 individuals gives: the first and the last position of each chromosome
# of the table above.
REGIONS_OF_ONE_BATCH = [[("chr1", 1000, 10213), ("chr2", 10250, 19463)]]


@dataclass(frozen=True)
class _WhatPyarrowReads:
    """A vars file as pyarrow opens it: its schema, the value of each of its
    two keys, how many variants each of its batches really holds, and its
    columns."""

    schema: pyarrow.Schema
    popnei: dict[str, Any]
    batches: list[dict[str, Any]]
    num_vars_of_each_batch: list[int]
    table: pyarrow.Table


def _pyarrow_reads(path: Path) -> _WhatPyarrowReads:
    """The vars file at `path`, opened with pyarrow.

    `pyarrow.ipc.open_file` reads the schema and the footer, and it gives
    the key of the footer as the `metadata` of the file it opens and the key
    of the schema in `schema.metadata`. Both values are json.
    """
    with pyarrow.ipc.open_file(path) as reader:
        footer = reader.metadata
        return _WhatPyarrowReads(
            schema=reader.schema,
            popnei=json.loads(reader.schema.metadata[POPNEI_KEY]),
            batches=json.loads(footer[POPNEI_BATCHES_KEY]),
            num_vars_of_each_batch=[
                reader.get_batch(batch).num_rows
                for batch in range(reader.num_record_batches)
            ],
            table=reader.read_all(),
        )


def _columns_found(read: _WhatPyarrowReads) -> list[tuple[str, pyarrow.DataType, bool]]:
    """The columns of the file, to compare with `_columns_of`."""
    return [(field.name, field.type, field.nullable) for field in read.schema]


def _regions_found(read: _WhatPyarrowReads) -> list[list[tuple[str, int, int]]]:
    """The regions of every batch of the file, to compare with the tables
    above."""
    return [
        [
            (region["chrom"], region["min_pos"], region["max_pos"])
            for region in batch["regions"]
        ]
        for batch in read.batches
    ]


def _alleles_of(genotype: str) -> list[int]:
    """The alleles of a genotype as bcftools prints it, `0|1` or `./.`: the
    numbers of the alleles the variant declares, and -1 for the dot of an
    allele that was not called."""
    return [
        MISSING_ALLELE if allele == MISSING_VALUE else int(allele)
        for allele in genotype.replace("|", "/").split("/")
    ]


def _rows_of_bcftools(path: Path, only_passed: bool) -> list[dict[str, Any]]:
    """The variants of `many.bcftools.tsv`, each as the columns a vars file
    holds for it: with `only_passed` the rows whose FILTER is `PASS` or a
    dot, which is what the reader gives by default, and without it every
    row."""
    rows = []
    for line in path.read_text().splitlines():
        if not line:
            continue
        columns = line.split("\t")
        chrom, pos, id_, ref, alt, qual, kept = columns[:COLUMNS_BEFORE_THE_GENOTYPES]
        if only_passed and kept not in ("PASS", MISSING_VALUE):
            continue
        rows.append(
            {
                "chrom": chrom,
                "pos": int(pos),
                "id": None if id_ == MISSING_VALUE else id_,
                "alleles": [ref, *alt.split(",")],
                "qual": None if qual == MISSING_VALUE else float(qual),
                "gts": [
                    allele
                    for genotype in columns[COLUMNS_BEFORE_THE_GENOTYPES:]
                    for allele in _alleles_of(genotype)
                ],
            }
        )
    return rows


@pytest.mark.parametrize(
    (
        "only_passed",
        "num_vars_per_block",
        "num_vars_of_each_batch",
        "regions",
        "size_of_the_key",
        "nulls",
    ),
    [
        (False, 100, [100] * 5, REGIONS_OF_EVERY_VARIANT, 100, (167, 100)),
        (True, 100, [100, 100, 100, 100, 75], REGIONS_OF_THE_PASSED, 100, None),
        (
            False,
            None,
            [500],
            REGIONS_OF_ONE_BATCH,
            LARGEST_NUM_VARS_PER_BLOCK,
            (167, 100),
        ),
    ],
    ids=["every variant in batches of 100", "by default", "the size popnei chooses"],
)
def test_write_vars_writes_many_vcf_as_pyarrow_reads_it_back(
    reference_vcf_dir: Path,
    tmp_path: Path,
    only_passed: bool,
    num_vars_per_block: int | None,
    num_vars_of_each_batch: list[int],
    regions: list[list[tuple[str, int, int]]],
    size_of_the_key: int,
    nulls: tuple[int, int] | None,
) -> None:
    """The 500 variants of 50 individuals of `many.vcf`, written and opened.

    `nulls` is how many null ids and how many null qualities the file has,
    the variants with a dot in those columns of `many.bcftools.tsv`, which
    the spec counts for the 500 variants of the file and not for the 475
    that passed a filter.
    """
    variants = open_vcf(reference_vcf_dir / "many.vcf", only_passed=only_passed)
    path = tmp_path / "many.vars"

    write_vars(variants, path, num_vars_per_block)

    read = _pyarrow_reads(path)
    assert _columns_found(read) == _columns_of(MANY_ALLELES_PER_VAR)
    assert read.popnei["format_version"] == FORMAT_VERSION
    assert read.popnei["individuals"] == list(variants.individuals)
    assert len(read.popnei["individuals"]) == MANY_NUM_INDIVIDUALS
    assert read.popnei["ploidy"] == 2
    assert read.popnei["num_vars_per_block"] == size_of_the_key
    assert [batch["num_vars"] for batch in read.batches] == num_vars_of_each_batch
    assert read.num_vars_of_each_batch == num_vars_of_each_batch
    assert _regions_found(read) == regions
    if nulls is not None:
        null_ids, null_qualities = nulls
        assert read.table.column("id").null_count == null_ids
        assert read.table.column("qual").null_count == null_qualities

    expected = _rows_of_bcftools(reference_vcf_dir / "many.bcftools.tsv", only_passed)
    columns = read.table.to_pydict()
    assert read.table.num_rows == len(expected)
    # One variant at a time, so that a file that differs in one of them says
    # which one and not that two long lists differ.
    for index, row in enumerate(expected):
        found = {name: columns[name][index] for name in row}
        assert found == row, f"the variant {index}, counted from 0"


def test_write_vars_writes_a_file_of_no_batch_for_a_source_with_no_variants(
    write_vcf, tmp_path: Path
) -> None:
    """A VCF with a header and no data line, which is not an error.

    The columns of a file are those of its first block, and a file of no
    block has the one column every vars file has, the genotypes, of the
    width of its three diploid individuals.
    """
    variants = open_vcf(write_vcf([]))
    path = tmp_path / "no_variant.vars"

    write_vars(variants, path)

    read = _pyarrow_reads(path)
    assert _columns_found(read) == [_columns_of(6)[-1]]
    assert read.popnei == {
        "format_version": FORMAT_VERSION,
        "individuals": ["ind1", "ind2", "ind3"],
        "ploidy": 2,
        "num_vars_per_block": LARGEST_NUM_VARS_PER_BLOCK,
    }
    assert read.batches == []
    assert read.num_vars_of_each_batch == []
    assert read.table.num_rows == 0


def test_write_vars_leaves_no_file_when_the_vcf_fails_half_way(
    write_vcf, tmp_path: Path
) -> None:
    """A tetraploid genotype in the third variant of a diploid VCF.

    The error is the one of the source and names the VCF, and the path the
    user wrote to is free afterwards, so the same call can be made again
    once the VCF is fixed. pyNei leaves the file it was writing.
    """
    vcf_path = write_vcf(
        [
            "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
            "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
            "chr1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t0/0/1/1\t0/1\t1/1",
        ]
    )
    path = tmp_path / "half_way.vars"

    with pytest.raises(ValueError, match="ind1") as refusal:
        write_vars(open_vcf(vcf_path), path)

    assert str(refusal.value).startswith(str(vcf_path))
    assert not path.exists()


def test_write_vars_refuses_a_path_that_a_file_is_already_at(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """The second call writes nothing and leaves the first file as it was."""
    variants = open_vcf(reference_vcf_dir / "cases.vcf")
    path = tmp_path / "cases.vars"
    write_vars(variants, path)
    written = path.read_bytes()

    with pytest.raises(ValueError, match="already") as refusal:
        write_vars(variants, path)

    assert str(refusal.value).startswith(str(path))
    assert path.read_bytes() == written


def test_what_a_user_reads_of_write_vars_is_written_in_the_package() -> None:
    """The private module explains nothing; the package is the API.

    A user who calls ``help`` on it reads the size of the batches, what a
    path that is taken gives and what a call that fails leaves, and the
    module of the binding crate carries none of that.
    """
    assert _core.write_vars.__doc__ is None
    assert write_vars.__doc__ is not None


def test_write_vars_refuses_a_num_vars_per_block_that_counts_no_variants(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """A negative size, which this call refuses before it makes a file, and
    a size of 0, which the core refuses once the file is made and which
    leaves nothing at the path."""
    variants = open_vcf(reference_vcf_dir / "cases.vcf")
    path = tmp_path / "cases.vars"

    with pytest.raises(ValueError, match="num_vars_per_block") as refusal:
        write_vars(variants, path, -5)
    assert "-5" in str(refusal.value)
    assert not path.exists()

    with pytest.raises(ValueError, match="0 variants"):
        write_vars(variants, path, 0)
    assert not path.exists()
