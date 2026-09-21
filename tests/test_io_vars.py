"""`write_vars` and `open_vars`: the vars file a Python user writes and reads.

A vars file is one arrow IPC file, and the reference outside the project is
pyarrow, the arrow implementation that Apache Arrow publishes, which opens
what popnei wrote as any other program would. The cases of the writer are
the ones `docs/specs/io_vars.md` gives to pytest under "How it is verified"
of the writer: `many.vcf` written and read back, a source with no variants,
a VCF that fails half way, and a path that a file is already at. The numbers
are the literals of that spec, and what the columns hold is compared with
`many.bcftools.tsv`, what bcftools 1.24 read in the same VCF.

The reader is checked on the file that the writer makes from `many.vcf`:
what `open_vars` gives, field by field, is what `open_vcf` gives for the VCF
itself, and it holds the counts that `docs/specs/io_vcf.md` has from
bcftools. Its other cases are the files it refuses, one that is not a vars
file, one that is not there, one that was cut short and one whose buffers
are compressed with zstd, which no build of popnei reads.

pyNei is not run here: it reads another file.
"""

import errno
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy
import pyarrow
import pyarrow.ipc
import pytest
from popnei import _core, open_vars, open_vcf, write_vars

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
    `docs/specs/io_vars.md` gives them.

    The values inside the two lists, `alleles` and `gts`, have a field of
    their own, named `item` as pyarrow names it and holding no null, which
    is what "What it holds" says of them: no allele and no genotype popnei
    writes is a null, and an allele that was not called is -1. A field that
    took nulls would cost a mask of ones beside every value.
    """
    allele = pyarrow.field("item", pyarrow.string(), nullable=False)
    genotype = pyarrow.field("item", pyarrow.int8(), nullable=False)
    return [
        ("chrom", pyarrow.string(), False),
        ("pos", pyarrow.uint64(), False),
        ("id", pyarrow.string(), True),
        ("alleles", pyarrow.list_(allele), False),
        ("qual", pyarrow.float32(), True),
        ("gts", pyarrow.list_(genotype, alleles_per_var), False),
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


def _as_float32(quality: str) -> float:
    """The quality that bcftools printed, as the `qual` column of a vars
    file holds it.

    That column is `Float32`, which is the type the QUAL of a VCF has, and
    pyarrow gives its values back as the floats of Python that they widen
    to. The text is turned into a `float32` here, so that a quality of the
    file which no `float32` holds exactly, `0.1`, is compared with what such
    a column can hold and not with the `float` of Python that the text
    reads as.
    """
    return float(numpy.float32(quality))


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
                "qual": None if qual == MISSING_VALUE else _as_float32(qual),
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


# How long the test below waits before it sends itself the signal of a
# Ctrl-C. Writing the VCF above as a vars file took 0.313 s in the build
# `maturin develop` makes, on the owner's Apple M5 Pro, so the call is in
# the middle of its pass over the source when the signal arrives.
SECONDS_BEFORE_THE_CTRL_C = 0.1


def test_a_ctrl_c_while_write_vars_runs_is_raised_and_leaves_no_file(
    long_vcf: Path, tmp_path: Path
) -> None:
    """A user who stops a call that is writing a big file.

    The signal arrives while the pass over the source runs, with the
    interpreter released, so Python raises it when the call is over. What
    the user gets is the `KeyboardInterrupt` of any Ctrl-C, and the path is
    free for the call they make again.
    """
    variants = open_vcf(long_vcf, only_passed=False)
    path = tmp_path / "stopped.vars"
    # The default handler of SIGINT is the one that raises
    # `KeyboardInterrupt`, and it is put back where the test found it.
    handler = signal.signal(signal.SIGINT, signal.default_int_handler)
    ctrl_c = threading.Timer(
        SECONDS_BEFORE_THE_CTRL_C, lambda: os.kill(os.getpid(), signal.SIGINT)
    )

    try:
        ctrl_c.start()
        with pytest.raises(KeyboardInterrupt):
            write_vars(variants, path)
    finally:
        ctrl_c.cancel()
        signal.signal(signal.SIGINT, handler)

    assert not path.exists()


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


def test_write_vars_of_the_private_module_names_the_type_it_was_given(
    tmp_path: Path,
) -> None:
    """What a user who calls `popnei._core` themselves reads.

    The package refuses what is not a `Variants` before the binding crate
    sees it, so nothing a user writes reaches this message; when it is
    read, it names the type that was given and where a source comes from,
    and no file is made at the path.
    """
    path = tmp_path / "nothing.vars"

    with pytest.raises(TypeError, match="open_vars") as refusal:
        _core.write_vars(123, path, None, _core.Steps())

    assert "`int`" in str(refusal.value)
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


# The directory where `tests/reference/vars/make_reference.py` writes the
# vars file that popnei cannot write itself, the one compressed with zstd.
REFERENCE_VARS_DIR = Path(__file__).parent / "reference" / "vars"

# Every field a block can carry besides the genotypes.
ALL_FIELDS = ("chrom", "pos", "id", "alleles", "qual")

# How many variants a batch of the file the reader is read on holds. The
# blocks a user asks for are cut where they ask, 7 variants or the size
# popnei chooses, so they are never the batches of that file.
VARS_NUM_VARS_PER_BLOCK = 100

# What `many.vcf` holds when every variant of it is read, from the table of
# "How it is verified" of `docs/specs/io_vcf.md`, which bcftools 1.24 gave:
# the variants, those of them in `chr2`, the genotypes with an allele that
# was not called, the alleles that were not called, the alleles that were,
# and the sum of the numbers of those. The vars file written from that VCF
# gives them back.
MANY_NUM_VARS = 500
MANY_IN_CHR2 = 250
MANY_MISSING_GENOTYPES = 1511
MANY_MISSING_ALLELES = 2765
MANY_CALLED_ALLELES = 47235
MANY_SUM_OF_THE_CALLED_ALLELES = 25954

# How many bytes are cut off the end of a whole vars file to make one that
# was damaged after it was written: enough to take its footer away, which is
# where a reader looks when the file is opened.
BYTES_CUT_OFF_THE_END = 100

# What a message of an arrow IPC file starts with: four bytes that mark a
# continuation and four that say how long its metadata is. A file is
# `ARROW1` and its padding, then the message of the schema, then one
# message for each batch, so the first batch lies its metadata and those
# eight bytes after the start of the schema, which is the first message of
# the file. The test below sets the first 16 bytes of that batch to zero,
# which leaves no message of a batch where the footer says one is.
A_MESSAGE_CONTINUES = b"\xff\xff\xff\xff"
BYTES_BEFORE_THE_METADATA_OF_A_MESSAGE = 8
BYTES_ZEROED_IN_THE_BATCH = 16

# The batches of the file that is written again from what `open_vars` reads,
# 500 variants in batches of 37, and the blocks of 7 variants and the one
# block that the size popnei chooses for 50 individuals, 10000, gives for a
# file of 500.
AGAIN_NUM_VARS_PER_BLOCK = 37
BATCHES_OF_THE_FILE_WRITTEN_AGAIN = [37] * 13 + [19]
BLOCKS_OF_SEVEN = [7] * 71 + [3]


def _joined(variants, fields=ALL_FIELDS, num_vars_per_block=None) -> dict[str, Any]:
    """The columns of every block of `variants`, one after another, and how
    many variants each of those blocks held."""
    blocks = list(
        variants.iter_blocks(fields=fields, num_vars_per_block=num_vars_per_block)
    )
    return {
        "num_vars_of_each_block": [block.num_vars for block in blocks],
        "gts": numpy.concatenate([block.gts for block in blocks]),
        "chrom": tuple(name for block in blocks for name in block.chrom),
        "pos": tuple(int(pos) for block in blocks for pos in block.pos),
        "id": tuple(id_ for block in blocks for id_ in block.id),
        "alleles": tuple(alleles for block in blocks for alleles in block.alleles),
        "qual": numpy.concatenate([block.qual for block in blocks]),
    }


def _assert_the_same_variants(ours: dict[str, Any], theirs: dict[str, Any]) -> None:
    """Every column of two passes, field by field."""
    numpy.testing.assert_array_equal(ours["gts"], theirs["gts"])
    assert ours["chrom"] == theirs["chrom"]
    assert ours["pos"] == theirs["pos"]
    assert ours["id"] == theirs["id"]
    assert ours["alleles"] == theirs["alleles"]
    # assert_array_equal takes two NaNs at the same place as equal, which is
    # what a variant with no quality gives on both sides.
    numpy.testing.assert_array_equal(ours["qual"], theirs["qual"])


@pytest.fixture(scope="module")
def many_vars(tmp_path_factory) -> Path:
    """`many.vcf`, every variant of it, as a vars file of batches of 100.

    It is written once for the module: the reader is read on it several
    times and the file does not change.
    """
    path = tmp_path_factory.mktemp("vars") / "many.vars"
    variants = open_vcf(REFERENCE_VCF_DIR / "many.vcf", only_passed=False)
    write_vars(variants, path, VARS_NUM_VARS_PER_BLOCK)
    return path


@pytest.mark.parametrize(
    ("num_vars_per_block", "num_vars_of_each_block"),
    [(7, BLOCKS_OF_SEVEN), (None, [MANY_NUM_VARS])],
    ids=["blocks of 7", "the size popnei chooses"],
)
def test_open_vars_gives_the_variants_of_the_vcf_the_file_was_written_from(
    many_vars: Path,
    num_vars_per_block: int | None,
    num_vars_of_each_block: list[int],
) -> None:
    """The 500 variants of `many.vcf`, written and read back.

    Every field of every variant is the one the VCF reader gives, which
    `docs/specs/io_vcf.md` checked against bcftools 1.24 and
    `docs/specs/block.md` against pyNei, so this carries those checks over
    to the file. The blocks that come out are cut where the user asked and
    not where the batches of the file are: 7 variants, which is fewer than
    the 100 of a batch, and the size popnei chooses for 50 individuals,
    10000, which is more than the file holds and gives one block of 500.
    """
    variants = open_vars(many_vars)
    from_the_vcf = open_vcf(REFERENCE_VCF_DIR / "many.vcf", only_passed=False)

    assert variants.individuals == from_the_vcf.individuals
    assert variants.num_individuals == from_the_vcf.num_individuals
    assert variants.ploidy == from_the_vcf.ploidy

    ours = _joined(variants, num_vars_per_block=num_vars_per_block)
    _assert_the_same_variants(ours, _joined(from_the_vcf))
    assert ours["num_vars_of_each_block"] == num_vars_of_each_block

    gts = ours["gts"]
    missing = gts == MISSING_ALLELE
    assert gts.shape[0] == MANY_NUM_VARS
    assert ours["chrom"].count("chr2") == MANY_IN_CHR2
    assert int(missing.any(axis=2).sum()) == MANY_MISSING_GENOTYPES
    assert int(missing.sum()) == MANY_MISSING_ALLELES
    assert int((~missing).sum()) == MANY_CALLED_ALLELES
    assert int(gts[~missing].sum()) == MANY_SUM_OF_THE_CALLED_ALLELES


def test_every_pass_over_what_open_vars_gives_reads_the_file_again(
    many_vars: Path,
) -> None:
    """A `Variants` can be given to any number of calculations.

    The vars file is opened again at every pass, as the VCF is, so the
    second pass holds the 500 variants of the first and not the none that a
    reader kept from one pass to the next would have left.
    """
    variants = open_vars(many_vars)

    first = _joined(variants, num_vars_per_block=13)
    second = _joined(variants)

    assert first["gts"].shape[0] == MANY_NUM_VARS
    _assert_the_same_variants(first, second)


def test_open_vars_gives_the_genotypes_alone_when_no_other_field_is_asked_for(
    many_vars: Path,
) -> None:
    """A column that nobody asked for is not decompressed and is `None`."""
    blocks = list(open_vars(many_vars).iter_blocks(fields=()))

    gts = numpy.concatenate([block.gts for block in blocks])
    from_the_vcf = open_vcf(REFERENCE_VCF_DIR / "many.vcf", only_passed=False)
    numpy.testing.assert_array_equal(gts, _joined(from_the_vcf)["gts"])
    for block in blocks:
        assert block.chrom is None
        assert block.pos is None
        assert block.id is None
        assert block.alleles is None
        assert block.qual is None


def test_open_vars_refuses_a_file_that_is_not_a_vars_file(
    reference_vcf_dir: Path,
) -> None:
    """A VCF, which is not an arrow file at all.

    It is refused at the call and not at the first block, because
    `open_vars` reads the schema and the footer of the file.
    """
    path = reference_vcf_dir / "many.vcf"

    with pytest.raises(ValueError, match="vars file") as refusal:
        open_vars(path)

    assert str(refusal.value).startswith(str(path))


def test_open_vars_gives_the_error_of_the_file_system_for_a_path_of_no_file(
    tmp_path: Path,
) -> None:
    """`FileNotFoundError` derives from `OSError`, and the path is in
    `filename`, where the standard library puts it."""
    path = tmp_path / "there_is_no_such_vars_file.vars"

    with pytest.raises(OSError) as refusal:
        open_vars(path)

    assert isinstance(refusal.value, FileNotFoundError)
    assert refusal.value.errno == errno.ENOENT
    assert refusal.value.filename == str(path)


def test_open_vars_refuses_a_vars_file_that_was_cut_short(
    many_vars: Path, tmp_path: Path
) -> None:
    """A whole file without its last 100 bytes, which took its footer away.

    A file that was damaged after it was written is an error and not a file
    of fewer variants. It is an error of the input and not of what the user
    wrote, so it is an `OSError` with the path in `filename`, and the file
    is read when it is opened, so it comes at the call.
    """
    path = tmp_path / "cut_short.vars"
    path.write_bytes(many_vars.read_bytes()[:-BYTES_CUT_OFF_THE_END])

    with pytest.raises(OSError) as refusal:
        open_vars(path)

    assert refusal.value.filename == str(path)


def _where_the_first_batch_is(written: bytes) -> int:
    """Where the message of the first batch of a vars file starts."""
    schema = written.index(A_MESSAGE_CONTINUES)
    metadata = written[schema + 4 : schema + BYTES_BEFORE_THE_METADATA_OF_A_MESSAGE]
    first_batch = (
        schema
        + BYTES_BEFORE_THE_METADATA_OF_A_MESSAGE
        + int.from_bytes(metadata, "little")
    )
    where_it_continues = written[first_batch : first_batch + len(A_MESSAGE_CONTINUES)]
    assert where_it_continues == A_MESSAGE_CONTINUES, (
        f"the message after the schema, at {first_batch}, is not one of a batch"
    )
    return first_batch


def test_open_vars_refuses_a_batch_of_a_vars_file_that_was_damaged(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """A whole file whose first batch was written over with zeros.

    The footer is untouched, so the file opens and says what it holds, and
    the batch is refused when it is read. It is a file that was damaged
    after it was written and not something a user wrote, so it is an
    ``OSError`` with the path in ``filename`` and no ``errno``: nothing of
    the file system refused anything.
    """
    whole = tmp_path / "whole.vars"
    write_vars(open_vcf(reference_vcf_dir / "cases.vcf", only_passed=False), whole)
    written = bytearray(whole.read_bytes())
    batch = _where_the_first_batch_is(written)
    written[batch : batch + BYTES_ZEROED_IN_THE_BATCH] = bytes(
        BYTES_ZEROED_IN_THE_BATCH
    )
    path = tmp_path / "damaged.vars"
    path.write_bytes(written)

    variants = open_vars(path)

    assert variants.num_individuals == 3
    with pytest.raises(OSError, match="batch") as refusal:
        list(variants.iter_blocks())
    assert refusal.value.filename == str(path)
    assert refusal.value.errno is None


def test_open_vars_refuses_a_file_compressed_with_zstd_at_its_first_block() -> None:
    """The file that `tests/reference/vars/make_reference.py` writes.

    No build of popnei carries the zstd crate, natively or in wasm, so a
    file whose buffers are zstd is refused everywhere. Arrow decompresses a
    batch when it reads it, so the file opens, its individuals come out of
    the key of its schema, and the error comes with the first block.
    """
    path = REFERENCE_VARS_DIR / "zstd.vars"

    variants = open_vars(path)

    assert variants.individuals == ("ind1", "ind2", "ind3")
    with pytest.raises(ValueError, match="zstd") as refusal:
        list(variants.iter_blocks())
    assert str(refusal.value).startswith(str(path))


def test_write_vars_writes_again_what_open_vars_reads_with_another_size_of_block(
    many_vars: Path, tmp_path: Path
) -> None:
    """A vars file of batches of 100 written again in batches of 37.

    A `Variants` of a vars file is a source like the one of a VCF, so it is
    written as any other, and the variants of the second file are those of
    the first, field by field.
    """
    path = tmp_path / "again.vars"

    write_vars(open_vars(many_vars), path, AGAIN_NUM_VARS_PER_BLOCK)

    read = _pyarrow_reads(path)
    assert [batch["num_vars"] for batch in read.batches] == (
        BATCHES_OF_THE_FILE_WRITTEN_AGAIN
    )
    _assert_the_same_variants(_joined(open_vars(path)), _joined(open_vars(many_vars)))


def test_write_vars_names_the_source_when_the_vars_file_it_reads_is_damaged(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """Both files of the call are vars files, and one of them is damaged.

    What went wrong is the batch of the source, so the exception carries
    that path and not the path of the file the call was writing, which is
    taken away as after any error.
    """
    source = tmp_path / "source.vars"
    write_vars(open_vcf(reference_vcf_dir / "cases.vcf", only_passed=False), source)
    written = bytearray(source.read_bytes())
    batch = _where_the_first_batch_is(written)
    written[batch : batch + BYTES_ZEROED_IN_THE_BATCH] = bytes(
        BYTES_ZEROED_IN_THE_BATCH
    )
    source.write_bytes(written)
    path = tmp_path / "again.vars"

    with pytest.raises(OSError, match="batch") as refusal:
        write_vars(open_vars(source), path)

    assert refusal.value.filename == str(source)
    assert not path.exists()


def test_what_a_user_reads_of_open_vars_is_written_in_the_package() -> None:
    """The private module explains nothing; the package is the API."""
    assert _core.open_vars.__doc__ is None
    assert open_vars.__doc__ is not None
