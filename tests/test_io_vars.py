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

import errno
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pyarrow
import pyarrow.ipc
import pytest
from popnei import _core, open_vcf, write_vars

# The VCFs that `tests/reference/vcf/make_reference.py` writes, which the
# module scoped fixture below reads and which `conftest.py` gives the tests
# as `reference_vcf_dir`.
REFERENCE_VCF_DIR = Path(__file__).parent / "reference" / "vcf"

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


# How many variants the VCF of the two tests that act while a call runs
# holds. Writing it as a vars file takes about a third of a second in the
# build `maturin develop` makes, which the tests run against, so a thread of
# the test has the time to act between the moment the file is made and the
# moment the call is over.
VARIANTS_OF_THE_LONG_VCF = 100_000

# A data line of 50 individuals whose first genotype holds four alleles,
# which a reader of diploid genotypes refuses: the last line of the VCF of
# the test below, so that the call fails once it has written most of its
# file.
_A_TETRAPLOID_LINE = (
    "chr1\t20000\t.\tA\tT\t.\tPASS\t.\tGT\t0/0/1/1\t"
    + "\t".join(["0/0"] * (MANY_NUM_INDIVIDUALS - 1))
    + "\n"
)

# How long a thread of a test waits for the call it watches to make its
# file, and how often it looks: the file is made when the call begins, so it
# is there within a few of these.
LOOKS_FOR_THE_FILE = 5000
SECONDS_BETWEEN_LOOKS = 0.001


@pytest.fixture(scope="module")
def long_vcf(tmp_path_factory) -> Path:
    """A VCF of 100000 variants of the 50 individuals of `many.vcf`.

    Its data lines are those of `many.vcf`, one after another until there
    are that many, so its variants are not sorted, which the writer takes.
    It is made once for the module: writing it is 23 MB of text.
    """
    path = tmp_path_factory.mktemp("long") / "long.vcf"
    lines = (REFERENCE_VCF_DIR / "many.vcf").read_text().splitlines(keepends=True)
    header = [line for line in lines if line.startswith("#")]
    data = [line for line in lines if not line.startswith("#")]
    times = -(-VARIANTS_OF_THE_LONG_VCF // len(data))
    path.write_text("".join(header + (data * times)[:VARIANTS_OF_THE_LONG_VCF]))
    return path


def _stop_writing_in_the_directory_once_the_file_is_there(
    path: Path, directory: Path
) -> None:
    """Waits for `path` to be made and takes away the right to write in
    `directory`, so that what is at `path` cannot be removed any more.

    The file is made when the call begins and is removed when it fails, so
    this lands between the two.
    """
    for _ in range(LOOKS_FOR_THE_FILE):
        if path.exists():
            directory.chmod(0o500)
            return
        time.sleep(SECONDS_BETWEEN_LOOKS)


def test_write_vars_says_when_it_could_not_take_away_the_file_it_was_writing(
    long_vcf: Path, tmp_path: Path
) -> None:
    """A directory that stops taking files while the vars file is written.

    What went wrong, the tetraploid genotype of the last line of the VCF, is
    the exception, and a note on it says that a file is still at the path,
    which the call the user makes again would refuse.
    """
    if os.geteuid() == 0:
        pytest.skip("a directory that says no still takes the files of root")
    vcf_path = tmp_path / "tetraploid_at_its_end.vcf"
    shutil.copyfile(long_vcf, vcf_path)
    with vcf_path.open("a") as vcf:
        vcf.write(_A_TETRAPLOID_LINE)
    path = tmp_path / "left_behind.vars"
    watcher = threading.Thread(
        target=_stop_writing_in_the_directory_once_the_file_is_there,
        args=(path, tmp_path),
        daemon=True,
    )

    try:
        watcher.start()
        with pytest.raises(ValueError, match="ind0") as refusal:
            write_vars(open_vcf(vcf_path), path)
    finally:
        watcher.join(timeout=10)
        tmp_path.chmod(0o700)

    assert path.exists()
    notes = getattr(refusal.value, "__notes__", [])
    assert any(str(path) in note for note in notes), notes


# What a child process does with a limit on the size of the files it may
# write: it writes a vars file that would go past the limit, and prints what
# the call raised as json. The limit is of the process and not of a
# directory, so it is set in a process of its own, after popnei is imported,
# and the signal that ends a process which writes past it is ignored, so
# that the write gets the error instead. Nothing here opens the file that
# was written: what the test asks is which file the error names.
_WRITE_PAST_A_LIMIT = '''
"""Writes a vars file in a process that may write only so many bytes."""

import json
import resource
import signal
import sys

import popnei

vcf_path, vars_path, limit = sys.argv[1], sys.argv[2], int(sys.argv[3])
signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
resource.setrlimit(resource.RLIMIT_FSIZE, (limit, limit))
answer = {"raised": None}
try:
    popnei.write_vars(popnei.open_vcf(vcf_path, only_passed=False), vars_path)
except OSError as error:
    answer = {
        "raised": type(error).__name__,
        "errno": error.errno,
        "filename": error.filename,
        "message": str(error),
    }
except BaseException as error:  # noqa: BLE001
    answer = {"raised": type(error).__name__, "message": str(error)}
print(json.dumps(answer))
'''

# How many bytes that child process may write, which is fewer than the vars
# file of `many.vcf` holds, so the write fails part way through it.
BYTES_THE_CHILD_MAY_WRITE = 4096


def test_write_vars_names_the_vars_file_when_the_write_is_what_failed(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """A process that may write 4096 bytes and a file that needs more.

    What fails is the write and not the read of the VCF, so the exception
    carries the path of the vars file, which is the file the user acts on,
    and says that it could not be written. `filename` is where Python keeps
    the file of an `OSError` and `errno` the number the system gave, 27 for
    a file that grew past what the process may write.
    """
    pytest.importorskip("resource", reason="the limit is of a Unix process")
    script = tmp_path / "write_past_a_limit.py"
    script.write_text(_WRITE_PAST_A_LIMIT)
    vars_path = tmp_path / "past_the_limit.vars"

    child = subprocess.run(
        [
            sys.executable,
            str(script),
            str(reference_vcf_dir / "many.vcf"),
            str(vars_path),
            str(BYTES_THE_CHILD_MAY_WRITE),
        ],
        capture_output=True,
        text=True,
        check=True,
    )

    answer = json.loads(child.stdout)
    assert answer["raised"] == "OSError", child.stderr
    assert answer["errno"] == errno.EFBIG
    assert answer["filename"] == str(vars_path)
    assert "could not be written" in answer["message"]
    assert not vars_path.exists()


def test_write_vars_gives_the_error_of_the_file_system_for_a_path_of_no_file(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """A directory where the vars file goes, and one that is not there.

    `create_new` says of a directory that something is already at the path,
    which would tell a user to take away the directory they meant to write
    into; what they get is the number the system gives for a directory where
    a file was asked for, the `IsADirectoryError` that `open_vcf` gives them
    for such a path. A path in a directory that is not there is the
    `FileNotFoundError` of that path.
    """
    variants = open_vcf(reference_vcf_dir / "cases.vcf")

    with pytest.raises(OSError) as refusal:
        write_vars(variants, tmp_path)
    assert isinstance(refusal.value, IsADirectoryError)
    assert refusal.value.errno == errno.EISDIR
    assert refusal.value.filename == str(tmp_path)

    of_no_directory = tmp_path / "no_such_directory" / "cases.vars"
    with pytest.raises(OSError) as refusal:
        write_vars(variants, of_no_directory)
    assert isinstance(refusal.value, FileNotFoundError)
    assert refusal.value.errno == errno.ENOENT
    assert refusal.value.filename == str(of_no_directory)
    assert not of_no_directory.exists()


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


def test_write_vars_says_what_it_takes_when_it_is_given_a_path(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """A user who gives the VCF where the variants of it go.

    It is the easiest mistake to make, and what it gave was the
    ``AttributeError`` of a `str` with no ``_source``. The refusal names the
    argument, says what was given and says that the variants come from
    `open_vcf`, as the refusal of a `fields` that is one name does.
    """
    path = tmp_path / "cases.vars"

    with pytest.raises(TypeError, match="open_vcf") as refusal:
        write_vars(str(reference_vcf_dir / "cases.vcf"), path)

    assert "variants" in str(refusal.value)
    assert not path.exists()


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
