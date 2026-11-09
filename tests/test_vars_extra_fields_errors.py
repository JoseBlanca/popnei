"""The errors of a vars file with extra fields, as a Python user meets them.

An extra field is a named array of numbers that a vars file of format 1.3
carries beside its columns, one value for each variant or one for each
genotype, "The extra fields" of `docs/specs/io_vars.md`. No function of the
package writes one yet, so each file here is built with pyarrow from the
file popnei writes of `cases.vcf`, four variants of three individuals, with
the columns of the fields added and their entries in the `popnei` key of
the schema. The cargo tests of `crates/popnei/src/io/vars` cover every case
of the reader; these check that each reaches Python as the exception the
spec gives, with the file and the field named.
"""

import json
from collections.abc import Iterable
from pathlib import Path

import pyarrow
import pyarrow.ipc
import pytest
from popnei import open_vars, open_vcf, write_vars, write_vcf

CASES_VCF = Path(__file__).parent / "reference" / "vcf" / "cases.vcf"
NUM_VARS = 4
NUM_INDIVIDUALS = 3


@pytest.fixture
def cases_vars(tmp_path: Path) -> Path:
    """The vars file popnei writes of `cases.vcf`, with no extra field."""
    path = tmp_path / "cases.vars"
    write_vars(open_vcf(CASES_VCF, only_passed=False), path)
    return path


def with_fields(
    base: Path,
    path: Path,
    var_fields: list[list[str]],
    gt_fields: list[list[str]],
    columns: Iterable[tuple[str, pyarrow.Array]],
) -> Path:
    """`base` written at `path` with `columns` after its own and the two
    entries of the extra fields in its `popnei` key."""
    with pyarrow.ipc.open_file(base) as reader:
        footer = reader.metadata
        metadata = dict(reader.schema.metadata)
        table = reader.read_all()
    key = json.loads(metadata[b"popnei"])
    key["var_fields"] = var_fields
    key["gt_fields"] = gt_fields
    metadata[b"popnei"] = json.dumps(key).encode()
    for name, column in columns:
        table = table.append_column(name, column)
    table = table.replace_schema_metadata(metadata)
    with pyarrow.ipc.new_file(path, table.schema, metadata=footer) as writer:
        for batch in table.to_batches():
            writer.write_batch(batch)
    return path


def end_column(values: list[int | None], arrow_type: pyarrow.DataType) -> pyarrow.Array:
    """A variant field of one value for each of the four variants."""
    return pyarrow.array(values, arrow_type)


def gp_column(values: list[float | None]) -> pyarrow.Array:
    """A genotype field of `float32`, one value for each individual of each
    variant."""
    return pyarrow.FixedSizeListArray.from_arrays(
        pyarrow.array(values, pyarrow.float32()), NUM_INDIVIDUALS
    )


def gp_values() -> list[float | None]:
    """A value for each genotype of the four variants, each another."""
    return [float(value) for value in range(NUM_VARS * NUM_INDIVIDUALS)]


def test_write_vcf_of_a_file_with_extra_fields_names_them_and_leaves_no_vcf(
    cases_vars: Path, tmp_path: Path
) -> None:
    source = with_fields(
        cases_vars,
        tmp_path / "fields.vars",
        [["END", "int64"]],
        [["GP", "float32"]],
        [
            ("END", end_column([150, 250, 350, 450], pyarrow.int64())),
            ("GP", gp_column(gp_values())),
        ],
    )
    vcf = tmp_path / "fields.vcf"

    with pytest.raises(ValueError) as refusal:
        write_vcf(open_vars(source), vcf)

    message = str(refusal.value)
    assert message.startswith(str(source))
    assert "the extra fields `END`, `GP`" in message
    assert not vcf.exists()


def test_open_vars_of_an_extra_column_of_another_type_names_the_file_and_the_field(
    cases_vars: Path, tmp_path: Path
) -> None:
    source = with_fields(
        cases_vars,
        tmp_path / "int32.vars",
        [["END", "int64"]],
        [],
        [("END", end_column([150, 250, 350, 450], pyarrow.int32()))],
    )

    with pytest.raises(ValueError) as refusal:
        open_vars(source)

    message = str(refusal.value)
    assert message.startswith(str(source))
    assert "the extra field `END`" in message
    assert "Int32" in message
    assert "Int64" in message


def test_open_vars_of_an_extra_field_with_no_column_names_the_file_and_the_field(
    cases_vars: Path, tmp_path: Path
) -> None:
    source = with_fields(
        cases_vars, tmp_path / "no_column.vars", [["END", "int64"]], [], []
    )

    with pytest.raises(ValueError) as refusal:
        open_vars(source)

    message = str(refusal.value)
    assert message.startswith(str(source))
    assert "the extra field `END` of the vars file has no column" in message


def test_a_null_in_a_variant_field_names_the_column_and_the_variant(
    cases_vars: Path, tmp_path: Path
) -> None:
    source = with_fields(
        cases_vars,
        tmp_path / "null_end.vars",
        [["END", "int64"]],
        [],
        [("END", end_column([150, None, 350, 450], pyarrow.int64()))],
    )
    written = tmp_path / "written.vars"

    # `iter_blocks` does not ask for the extra fields, so the null is
    # reached through `write_vars`, which reads every field.
    with pytest.raises(ValueError) as refusal:
        write_vars(open_vars(source), written)

    message = str(refusal.value)
    assert message.startswith(str(source))
    assert "the `END` column of the vars file has no value for its variant 2" in message
    assert not written.exists()


def test_a_null_in_a_genotype_field_names_the_column_and_the_variant(
    cases_vars: Path, tmp_path: Path
) -> None:
    values = gp_values()
    # The second individual of the third variant.
    values[2 * NUM_INDIVIDUALS + 1] = None
    source = with_fields(
        cases_vars,
        tmp_path / "null_gp.vars",
        [],
        [["GP", "float32"]],
        [("GP", gp_column(values))],
    )
    written = tmp_path / "written.vars"

    with pytest.raises(ValueError) as refusal:
        write_vars(open_vars(source), written)

    message = str(refusal.value)
    assert message.startswith(str(source))
    assert "the `GP` column of the vars file has no value for its variant 3" in message
    assert not written.exists()
