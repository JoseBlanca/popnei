"""The writer of segments, `write_vars_from_segments`, as pyarrow reads what
it writes and as a Python user meets what it refuses.

"How it is verified" of "The writer of segments" of `docs/specs/io_vars.md`
has the worked example of three segments that these tests write: three
individuals, `i1`, `i2` and `i3`, of the ploidy 2, a variant field `END` of
`uint64` and two genotype fields, `GP` of `float32` and `NODE_A` of
`int32`, in batches of 2 variants. The reference outside the project is
pyarrow, which opens the file as any other program would. The tests at the
end read the fields back through popnei: `var_fields` and `gt_fields`,
`iter_blocks` with `GP`, the filter of individuals, `write_vars` of the
filtered file and `write_vcf`.

Every test of a call that fails checks that it leaves no file at the path,
which is what "A call that fails" asks of each.
"""

import _thread
import json
import math
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy
import pyarrow
import pyarrow.ipc
import pytest
from popnei import open_vars, open_vcf, write_vars, write_vars_from_segments, write_vcf

INDIVIDUALS = ["i1", "i2", "i3"]
PLOIDY = 2
VAR_FIELDS = {"END": "uint64"}
GT_FIELDS = {"GP": "float32", "NODE_A": "int32"}
NAN = float("nan")


def first_segment() -> dict[str, Any]:
    """The first segment of the worked example: chr1 100 and 5000."""
    return {
        "chrom": "chr1",
        "pos": [100, 5000],
        "alleles": [["H0", "H1"], ["H0", "H1", "H2"]],
        "gts": numpy.array([[[0, 1], [1, 1], [-1, -1]], [[2, 0], [0, 0], [1, 2]]]),
        "END": [4999, 9999],
        "GP": [[0.9, 0.8, NAN], [0.7, 0.99, 0.5]],
        "NODE_A": [[1, 2, -1], [3, 1, 2]],
    }


def empty_segment() -> dict[str, Any]:
    """The second segment of the worked example, of no variants."""
    return {
        "chrom": "chr1",
        "pos": [],
        "alleles": [],
        "gts": numpy.zeros((0, 3, 2), dtype=numpy.int8),
        "END": [],
        "GP": numpy.zeros((0, 3), dtype=numpy.float32),
        "NODE_A": numpy.zeros((0, 3), dtype=numpy.int32),
    }


def third_segment() -> dict[str, Any]:
    """The third segment of the worked example: chr2 50, of one allele."""
    return {
        "chrom": "chr2",
        "pos": [50],
        "alleles": [["H0"]],
        "gts": numpy.array([[[0, 0], [0, 0], [0, -1]]]),
        "END": [800],
        "GP": [[1.0, 1.0, 0.6]],
        "NODE_A": [[1, 1, 1]],
    }


def segment_at(chrom: str, pos: int) -> dict[str, Any]:
    """A segment of one variant of two alleles at `chrom` `pos`."""
    return {
        "chrom": chrom,
        "pos": [pos],
        "alleles": [["A", "T"]],
        "gts": numpy.array([[[0, 1], [1, 1], [0, 0]]]),
        "END": [pos],
        "GP": [[0.5, 0.5, 0.5]],
        "NODE_A": [[1, 1, 1]],
    }


def write(segments: Any, path: Path, **arguments: Any) -> Any:
    """`write_vars_from_segments` of the worked example's arguments, with
    `arguments` in place of any of them."""
    given: dict[str, Any] = {
        "individuals": INDIVIDUALS,
        "ploidy": PLOIDY,
        "var_fields": VAR_FIELDS,
        "gt_fields": GT_FIELDS,
        "num_vars_per_block": 2,
    }
    given.update(arguments)
    return write_vars_from_segments(segments, path, **given)


def read_with_pyarrow(path: Path) -> tuple[pyarrow.Schema, pyarrow.Table, dict]:
    """The schema, the table and the footer of the vars file at `path`."""
    with pyarrow.ipc.open_file(path) as reader:
        return reader.schema, reader.read_all(), reader.metadata


def test_write_vars_from_segments_writes_the_worked_example_as_pyarrow_reads_it(
    tmp_path: Path,
) -> None:
    """The three segments of the worked example, the second of no variants:
    3 variants in 2 batches, with the columns, the types, the nulls, the
    values and the two keys that the spec gives."""
    path = tmp_path / "segments.vars"

    written = write(iter([first_segment(), empty_segment(), third_segment()]), path)

    assert written.pass_stats.num_vars == 3
    assert written.pass_stats.filtering == {}
    schema, table, footer = read_with_pyarrow(path)
    assert schema.names == [
        "chrom",
        "pos",
        "id",
        "alleles",
        "qual",
        "gts",
        "END",
        "GP",
        "NODE_A",
    ]
    fields = {field.name: field for field in schema}
    assert fields["END"].type == pyarrow.uint64()
    assert not fields["END"].nullable
    for name, value_type in [("GP", pyarrow.float32()), ("NODE_A", pyarrow.int32())]:
        assert pyarrow.types.is_fixed_size_list(fields[name].type)
        assert fields[name].type.list_size == 3
        assert fields[name].type.value_type == value_type
        assert not fields[name].nullable
    assert table.column("chrom").to_pylist() == ["chr1", "chr1", "chr2"]
    assert table.column("pos").to_pylist() == [100, 5000, 50]
    assert table.column("id").to_pylist() == [None, None, None]
    assert table.column("qual").to_pylist() == [None, None, None]
    assert table.column("alleles").to_pylist() == [
        ["H0", "H1"],
        ["H0", "H1", "H2"],
        ["H0"],
    ]
    assert table.column("gts").to_pylist() == [
        [0, 1, 1, 1, -1, -1],
        [2, 0, 0, 0, 1, 2],
        [0, 0, 0, 0, 0, -1],
    ]
    assert table.column("END").to_pylist() == [4999, 9999, 800]
    assert table.column("NODE_A").to_pylist() == [[1, 2, -1], [3, 1, 2], [1, 1, 1]]
    gp = numpy.array(table.column("GP").to_pylist(), dtype=numpy.float32)
    expected_gp = numpy.array(
        [[0.9, 0.8, NAN], [0.7, 0.99, 0.5], [1.0, 1.0, 0.6]], dtype=numpy.float32
    )
    numpy.testing.assert_array_equal(gp, expected_gp)
    key = json.loads(schema.metadata[b"popnei"])
    assert key["format_version"] == "1.3"
    assert key["var_fields"] == [["END", "uint64"]]
    assert key["gt_fields"] == [["GP", "float32"], ["NODE_A", "int32"]]
    assert key["individuals"] == INDIVIDUALS
    assert key["ploidy"] == 2
    assert key["num_vars_per_block"] == 2
    assert json.loads(footer[b"popnei_batches"]) == [
        {
            "num_vars": 2,
            "regions": [{"chrom": "chr1", "min_pos": 100, "max_pos": 5000}],
        },
        {"num_vars": 1, "regions": [{"chrom": "chr2", "min_pos": 50, "max_pos": 50}]},
    ]
    variants = open_vars(path)
    assert variants.individuals == ("i1", "i2", "i3")
    assert variants.ploidy == 2
    assert not variants.keeps_passed


def test_write_vars_from_segments_writes_ids_qualities_and_chrom_lengths(
    tmp_path: Path,
) -> None:
    """A segment with `id` and `qual`, a chromosome for each variant, a
    ploidy of 3, no field and the lengths of the chromosomes: an id of None
    and a quality of NaN are nulls in the file."""
    path = tmp_path / "with_ids.vars"
    segment = {
        "chrom": ["chr1", "chr2"],
        "pos": numpy.array([7, 3], dtype=numpy.uint32),
        "alleles": [("A", "C"), ("G", "T", "<DEL>")],
        "gts": numpy.array([[[0, 1, 1]], [[2, -1, 0]]], dtype=numpy.int64),
        "id": ["rs1", None],
        "qual": [30.5, NAN],
    }

    write(
        [segment],
        path,
        individuals=["only"],
        ploidy=3,
        var_fields={},
        gt_fields={},
        chrom_lengths={"chr1": 1000, "chr2": 2000},
        num_vars_per_block=None,
    )

    schema, table, _ = read_with_pyarrow(path)
    assert schema.names == ["chrom", "pos", "id", "alleles", "qual", "gts"]
    assert table.column("chrom").to_pylist() == ["chr1", "chr2"]
    assert table.column("pos").to_pylist() == [7, 3]
    assert table.column("id").to_pylist() == ["rs1", None]
    assert table.column("qual").to_pylist() == [30.5, None]
    assert table.column("gts").to_pylist() == [[0, 1, 1], [2, -1, 0]]
    key = json.loads(schema.metadata[b"popnei"])
    assert key["chrom_lengths"] == [["chr1", 1000], ["chr2", 2000]]
    assert key["var_fields"] == []
    assert key["gt_fields"] == []


def test_write_vars_from_segments_casts_lists_into_the_declared_types(
    tmp_path: Path,
) -> None:
    """A list of integers into a `uint64` field, including one above 2**63,
    which an `int64` does not hold, and a list of floats into a `float32`
    field, which holds the nearest `float32` of each; integers into a
    `float64` field and an `int8` array into an `int32` field."""
    path = tmp_path / "types.vars"
    segment = {
        "chrom": "chr1",
        "pos": [1, 2],
        "alleles": [["A", "C"], ["A", "C"]],
        "gts": numpy.zeros((2, 1, 2), dtype=numpy.int16),
        "BIG": [2**64 - 1, 5],
        "PROB": [[0.1], [1 / 3]],
        "WIDE": [[3], [-4]],
        "SMALL": numpy.array([[-7], [7]], dtype=numpy.int8),
    }

    write(
        [segment],
        path,
        individuals=["only"],
        var_fields={"BIG": numpy.uint64},
        gt_fields={"PROB": "f4", "WIDE": "float64", "SMALL": numpy.dtype("int32")},
    )

    schema, table, _ = read_with_pyarrow(path)
    assert table.column("BIG").to_pylist() == [2**64 - 1, 5]
    prob = numpy.array(table.column("PROB").to_pylist(), dtype=numpy.float32)
    numpy.testing.assert_array_equal(
        prob, numpy.array([[0.1], [1 / 3]], dtype=numpy.float32)
    )
    assert float(prob[0, 0]) != 0.1
    assert table.column("WIDE").to_pylist() == [[3.0], [-4.0]]
    assert table.column("SMALL").to_pylist() == [[-7], [7]]
    key = json.loads(schema.metadata[b"popnei"])
    assert key["var_fields"] == [["BIG", "uint64"]]
    assert key["gt_fields"] == [
        ["PROB", "float32"],
        ["WIDE", "float64"],
        ["SMALL", "int32"],
    ]


# The cases of a field array that its type does not take, each a change of
# the first segment of the worked example and words its message holds.
FIELD_ARRAYS_REFUSED = [
    pytest.param(
        "NODE_A",
        numpy.array([[1.0, 2.0, -1.0], [3.0, 1.0, 2.0]]),
        "`NODE_A` is of the type float64, and it holds integers",
        id="float64_array_into_an_int32_field",
    ),
    pytest.param(
        "END",
        [4999, -1],
        "segment 1, variant 2, chr1 5000: `END` holds -1, which its type, uint64,",
        id="negative_value_into_a_uint64_field",
    ),
    pytest.param(
        "NODE_A",
        [[1, 2, -1], [3, 2**31, 2]],
        "segment 1, variant 2, chr1 5000: `NODE_A` holds 2147483648, which its "
        "type, int32,",
        id="2_to_the_31_into_an_int32_field",
    ),
]


@pytest.mark.parametrize(("name", "value", "words"), FIELD_ARRAYS_REFUSED)
def test_write_vars_from_segments_refuses_a_field_array_its_type_does_not_take(
    tmp_path: Path, name: str, value: Any, words: str
) -> None:
    """A float array into an `int32` field, a negative value into a `uint64`
    one and 2**31 into an `int32` one, which numpy would round or wrap
    into another number."""
    path = tmp_path / "refused.vars"
    segment = first_segment()
    segment[name] = value

    with pytest.raises(ValueError) as refusal:
        write([segment], path)

    assert words in str(refusal.value)
    assert not path.exists()


def one_variant_with(**values: Any) -> dict[str, Any]:
    """A segment of one variant of the individual `only`, chr1 10, A C, with
    `values` beside its four required keys."""
    segment: dict[str, Any] = {
        "chrom": "chr1",
        "pos": [10],
        "alleles": [["A", "C"]],
        "gts": numpy.zeros((1, 1, 2), dtype=numpy.int8),
    }
    segment.update(values)
    return segment


def write_one_variant(path: Path, var_fields: dict[str, str], **values: Any) -> Any:
    """`one_variant_with(**values)` written with the variant fields
    `var_fields` and no genotype field."""
    return write(
        [one_variant_with(**values)],
        path,
        individuals=["only"],
        var_fields=var_fields,
        gt_fields={},
    )


# Each integer a field of floats is given, the type of that field, and
# whether the type holds it exactly: 2**24 and 2**53 are held, the integer
# after each is not, and 2**70, a power of two, is held by a float64.
INTEGERS_INTO_FLOATS = [
    pytest.param(2**24, "float32", True, id="2_to_the_24_into_float32"),
    pytest.param(2**24 + 1, "float32", False, id="2_to_the_24_plus_1_into_float32"),
    pytest.param(-(2**24) - 1, "float32", False, id="minus_2_to_the_24_minus_1"),
    pytest.param(2**53, "float64", True, id="2_to_the_53_into_float64"),
    pytest.param(2**53 + 1, "float64", False, id="2_to_the_53_plus_1_into_float64"),
    pytest.param(2**70, "float64", True, id="2_to_the_70_into_float64"),
    pytest.param(2**70 + 1, "float64", False, id="2_to_the_70_plus_1_into_float64"),
    pytest.param(2**64 - 1, "float64", False, id="2_to_the_64_minus_1_into_float64"),
]


@pytest.mark.parametrize(("integer", "type_name", "held"), INTEGERS_INTO_FLOATS)
def test_write_vars_from_segments_takes_an_integer_into_floats_only_when_held_exactly(
    tmp_path: Path, integer: int, type_name: str, held: bool
) -> None:
    """An integer in a list for a field of floats is written when the type
    holds it exactly, and is otherwise a `ValueError` that names the
    segment, the variant, the field and the integer, and leaves no file."""
    path = tmp_path / "integer.vars"

    if held:
        write_one_variant(path, {"X": type_name}, X=[integer])
        _, table, _ = read_with_pyarrow(path)
        assert table.column("X").to_pylist() == [float(integer)]
        assert int(table.column("X").to_pylist()[0]) == integer
        return
    with pytest.raises(ValueError) as refusal:
        write_one_variant(path, {"X": type_name}, X=[integer])
    assert (
        f"segment 1, variant 1, chr1 10: `X` holds the integer {integer}, which "
        f"its type, {type_name}, does not hold exactly"
    ) in str(refusal.value)
    assert not path.exists()


def test_write_vars_from_segments_refuses_an_int64_array_a_float32_field_rounds(
    tmp_path: Path,
) -> None:
    """An `int64` array, not a list, into a `float32` field: 2**24 + 1 is the
    `ValueError` as it is in a list."""
    path = tmp_path / "int64.vars"

    with pytest.raises(ValueError, match="`X` holds the integer 16777217"):
        write_one_variant(
            path, {"X": "float32"}, X=numpy.array([2**24 + 1], dtype=numpy.int64)
        )

    assert not path.exists()


def test_write_vars_from_segments_refuses_a_finite_float_beyond_a_float32(
    tmp_path: Path,
) -> None:
    """1e40 into a `float32` field would be an infinity: a `ValueError` that
    names it. An infinity given as one is taken."""
    path = tmp_path / "beyond.vars"

    with pytest.raises(ValueError) as refusal:
        write_one_variant(path, {"X": "float32"}, X=[1e40])

    assert "segment 1, variant 1, chr1 10: `X` holds 1e40, beyond the largest " in (
        str(refusal.value)
    )
    assert not path.exists()

    write_one_variant(path, {"X": "float32"}, X=[math.inf])
    _, table, _ = read_with_pyarrow(path)
    assert table.column("X").to_pylist() == [math.inf]


def test_write_vars_from_segments_refuses_a_quality_a_float32_does_not_hold(
    tmp_path: Path,
) -> None:
    """`qual` is a `float32` in the file and follows the rule of a field of
    floats: 2**24 + 1 and 1e40 are each a `ValueError`."""
    path = tmp_path / "qual.vars"

    with pytest.raises(ValueError, match="its quality is the integer 16777217"):
        write_one_variant(path, {}, qual=[2**24 + 1])
    with pytest.raises(ValueError, match="its quality is 1e40"):
        write_one_variant(path, {}, qual=[1e40])

    assert not path.exists()


def test_write_vars_from_segments_writes_the_bounds_of_an_int64_field(
    tmp_path: Path,
) -> None:
    """2**40 and -2**40, which an `int32` does not hold, into an `int64`
    field, as a list."""
    path = tmp_path / "int64.vars"
    segment = one_variant_with(pos=[10, 11], X=[2**40, -(2**40)])
    segment["alleles"] = [["A", "C"], ["A", "C"]]
    segment["gts"] = numpy.zeros((2, 1, 2), dtype=numpy.int8)

    write(
        [segment], path, individuals=["only"], var_fields={"X": "int64"}, gt_fields={}
    )

    _, table, _ = read_with_pyarrow(path)
    assert table.column("X").to_pylist() == [2**40, -(2**40)]


def test_write_vars_from_segments_takes_an_array_of_the_other_byte_order(
    tmp_path: Path,
) -> None:
    """An array of big endian integers and one of big endian floats, which
    numpy makes on a little endian machine as on any other: converted to the
    order of the machine and written with their values."""
    path = tmp_path / "big_endian.vars"

    write_one_variant(
        path,
        {"X": "uint64", "Y": "float32"},
        pos=numpy.array([10], dtype=">i4"),
        X=numpy.array([2**40 + 3], dtype=">u8"),
        Y=numpy.array([0.5], dtype=">f8"),
        qual=numpy.array([30.5], dtype=">f4"),
    )

    _, table, _ = read_with_pyarrow(path)
    assert table.column("pos").to_pylist() == [10]
    assert table.column("X").to_pylist() == [2**40 + 3]
    assert table.column("Y").to_pylist() == [0.5]
    assert table.column("qual").to_pylist() == [30.5]


# Each value of a key of numbers that is refused whatever the type of the
# key, as a change of the first segment of the worked example, and words of
# its message.
NUMBERS_REFUSED = [
    pytest.param(
        "END", [True, False], "`END` holds True, a truth value", id="list_of_bools"
    ),
    pytest.param(
        "NODE_A",
        [[1, True, -1], [3, 1, 2]],
        "`NODE_A` holds True, a truth value",
        id="one_bool_among_integers",
    ),
    pytest.param(
        "GP",
        [[0.9, 0.8, 0.1], [0.7, False, 0.5]],
        "`GP` holds False, a truth value",
        id="one_bool_among_floats",
    ),
    pytest.param(
        "GP",
        numpy.ones((2, 3), dtype=bool),
        "`GP` is of the type bool, and it holds numbers",
        id="array_of_bools_into_floats",
    ),
    pytest.param(
        "pos",
        numpy.array([True, True]),
        "`pos` is of the type bool, and it holds integers",
        id="array_of_bools_as_pos",
    ),
    pytest.param(
        "GP",
        numpy.ma.masked_array(numpy.ones((2, 3)), mask=numpy.zeros((2, 3))),
        "`GP` is a masked array",
        id="masked_array",
    ),
    pytest.param(
        "qual", 30.0, "`qual` is 30.0, one value", id="a_scalar_for_the_variants"
    ),
    pytest.param(
        "END",
        numpy.array(4999, dtype=numpy.uint64),
        "`END` is array(4999, dtype=uint64), one value",
        id="an_array_of_no_dimension",
    ),
    pytest.param(
        "END",
        [4999, None],
        "segment 1: `END` holds None, which is not a number",
        id="none_among_integers",
    ),
    pytest.param(
        "GP",
        [[0.9, 0.8, 0.1], [0.7, "5", 0.5]],
        "segment 1: `GP` holds '5', which is not a number",
        id="a_text_among_floats",
    ),
    pytest.param(
        "qual",
        numpy.array([1 + 2j, 3j]),
        "`qual` is of the type complex128, and it holds numbers",
        id="complex_qualities",
    ),
    pytest.param(
        "gts",
        [[[0, 1.0], [1, 1], [-1, -1]], [[2, 0], [0, 0], [1, 2]]],
        "segment 1, variant 1, chr1 100: `gts` holds the float 1.0, and it holds "
        "integers",
        id="a_float_in_a_list_of_gts",
    ),
    pytest.param(
        "gts",
        numpy.zeros((2, 3, 2), dtype=numpy.float16),
        "`gts` is of the type float16, and it holds integers",
        id="gts_of_float16",
    ),
    pytest.param(
        "NODE_A",
        numpy.zeros((2, 3), dtype=numpy.float16),
        "`NODE_A` is of the type float16, and it holds integers",
        id="float16_into_an_int32_field",
    ),
    pytest.param(
        "END",
        [4999, 2**70],
        "segment 1, variant 2, chr1 5000: `END` holds 1180591620717411303424, which "
        "its type, uint64, does not hold",
        id="2_to_the_70_into_a_uint64_field",
    ),
]


@pytest.mark.parametrize(("key", "value", "words"), NUMBERS_REFUSED)
def test_write_vars_from_segments_refuses_what_is_not_numbers_of_the_variants(
    tmp_path: Path, key: str, value: Any, words: str
) -> None:
    """Truth values, a masked array, a scalar, a value that is not a number,
    and floats where integers go, named by the type the user gave: each a
    `ValueError` that names the segment and the key, and no file."""
    path = tmp_path / "refused.vars"
    segment = first_segment()
    segment[key] = value

    with pytest.raises(ValueError) as refusal:
        write([segment], path)

    assert words in str(refusal.value)
    assert not path.exists()


def test_write_vars_from_segments_writes_a_variant_of_128_alleles_and_allele_127(
    tmp_path: Path,
) -> None:
    """A variant of exactly 128 alleles, the most, with a genotype of the
    allele 127, the last of them."""
    path = tmp_path / "many_alleles.vars"
    alleles = [f"H{number}" for number in range(128)]

    write_one_variant(
        path,
        {},
        alleles=[alleles],
        gts=numpy.array([[[127, 0]]], dtype=numpy.int64),
    )

    _, table, _ = read_with_pyarrow(path)
    assert table.column("alleles").to_pylist() == [alleles]
    assert table.column("gts").to_pylist() == [[127, 0]]


def test_write_vars_from_segments_writes_a_ploidy_of_255(tmp_path: Path) -> None:
    """The largest ploidy, 255, with every allele of the genotype given."""
    path = tmp_path / "ploidy.vars"
    gts = numpy.zeros((1, 1, 255), dtype=numpy.int8)
    gts[0, 0, 254] = 1

    write(
        [one_variant_with(gts=gts)],
        path,
        individuals=["only"],
        ploidy=255,
        var_fields={},
        gt_fields={},
    )

    variants = open_vars(path)
    assert variants.ploidy == 255
    _, table, _ = read_with_pyarrow(path)
    assert table.column("gts").to_pylist() == [[0] * 254 + [1]]


def changed(key: str, value: Any) -> dict[str, Any]:
    """The first segment of the worked example with `value` at `key`, or
    without `key` when `value` is `REMOVED`."""
    segment = first_segment()
    if value is REMOVED:
        del segment[key]
    else:
        segment[key] = value
    return segment


REMOVED = object()

# Each segment that "A call that fails" refuses, as a change of the first
# segment of the worked example, and words its message holds.
SEGMENTS_REFUSED = [
    pytest.param(
        "colour", "red", "segment 1: has the key `colour`", id="a_key_of_no_field"
    ),
    pytest.param("chrom", REMOVED, "has no `chrom`", id="no_chrom"),
    pytest.param("alleles", REMOVED, "has no `alleles`", id="no_alleles"),
    pytest.param("GP", REMOVED, "has no `GP`", id="no_declared_field"),
    pytest.param(
        "chrom",
        ["chr1"],
        "`chrom` holds 1 chromosomes and `pos` 2 positions",
        id="chroms_of_another_length",
    ),
    pytest.param(
        "alleles",
        [["H0", "H1"]],
        "`alleles` holds the alleles of 1 variants and `pos` 2 positions",
        id="alleles_of_another_length",
    ),
    pytest.param(
        "id", ["rs1"], "`id` holds 1 ids and `pos` 2 positions", id="ids_too_few"
    ),
    pytest.param(
        "qual", [1.0, 2.0, 3.0], "`qual` is of the shape [3]", id="quals_too_many"
    ),
    pytest.param(
        "gts",
        numpy.zeros((2, 3, 3), dtype=numpy.int8),
        "`gts` is of the shape [2, 3, 3], and it has to be of the shape (2, 3, 2)",
        id="gts_of_another_ploidy",
    ),
    pytest.param(
        "GP",
        [[0.9, 0.8], [0.7, 0.99]],
        "`GP` is of the shape [2, 2], and it has to be of the shape (2, 3)",
        id="genotype_field_of_another_width",
    ),
    pytest.param(
        "END", [[4999], [9999]], "`END` is of the shape [2, 1]", id="var_field_2d"
    ),
    pytest.param(
        "pos",
        numpy.array([100.0, 5000.0]),
        "`pos` is of the type float64, and it holds integers",
        id="pos_of_floats",
    ),
    pytest.param(
        "pos",
        [100, 5000.0],
        "segment 1, variant 2: its position on chr1 is 5000.0, and a position is a "
        "whole number",
        id="pos_list_with_a_float",
    ),
    pytest.param(
        "pos",
        [100, -5],
        "segment 1, variant 2: its position on chr1 is -5",
        id="negative_pos",
    ),
    pytest.param(
        "gts",
        numpy.zeros((2, 3, 2)),
        "`gts` is of the type float64, and it holds integers",
        id="gts_of_floats",
    ),
    pytest.param(
        "gts",
        numpy.array([[[0, 1], [1, 1], [-2, -1]], [[2, 0], [0, 0], [1, 2]]]),
        "segment 1, variant 1, chr1 100: `gts` holds the allele -2",
        id="allele_below_minus_1",
    ),
    pytest.param(
        "gts",
        numpy.array([[[0, 2], [1, 1], [-1, -1]], [[2, 0], [0, 0], [1, 2]]]),
        "segment 1, variant 1, chr1 100: `gts` holds the allele 2, which is not "
        "among the 2 alleles of the variant",
        id="allele_2_of_a_variant_of_two",
    ),
    pytest.param(
        "alleles",
        [[], ["H0", "H1", "H2"]],
        "segment 1, variant 1, chr1 100: it has no allele",
        id="variant_of_no_allele",
    ),
    pytest.param(
        "alleles",
        [["H0", ""], ["H0", "H1", "H2"]],
        "segment 1, variant 1, chr1 100: one of its alleles is an empty text",
        id="empty_allele",
    ),
    pytest.param(
        "alleles",
        [["H0", "H1"], [f"H{number}" for number in range(129)]],
        "segment 1, variant 2, chr1 5000: it has more than 128 alleles",
        id="129_alleles",
    ),
    pytest.param(
        "alleles",
        "H0H1",
        "`alleles` is H0H1, and it is a sequence of one value for each variant",
        id="alleles_as_one_text",
    ),
    pytest.param(
        "qual",
        [1.0, math.inf],
        "segment 1, variant 2, chr1 5000: its quality is inf",
        id="infinite_quality",
    ),
    pytest.param(
        "chrom",
        "chr 1",
        "the chromosome `chr 1` holds a space",
        id="space_in_a_chromosome",
    ),
    pytest.param(
        "chrom",
        ["chr1", "chr1\t"],
        "segment 1, variant 2: the chromosome `chr1\t` holds a tab",
        id="tab_in_a_chromosome",
    ),
    pytest.param(
        "alleles",
        [["H0", "H1,H2"], ["H0", "H1", "H2"]],
        "segment 1, variant 1, chr1 100: the allele `H1,H2` holds `,`",
        id="comma_in_an_allele",
    ),
    pytest.param(
        "alleles",
        [["H0", "H 1"], ["H0", "H1", "H2"]],
        "the allele `H 1` holds a space",
        id="space_in_an_allele",
    ),
    pytest.param(
        "id",
        ["rs1;rs2", None],
        "segment 1, variant 1, chr1 100: the id `rs1;rs2` holds `;`",
        id="semicolon_in_an_id",
    ),
    pytest.param(
        "id",
        ["rs1", "rs\n2"],
        "the id `rs\n2` holds a newline",
        id="newline_in_an_id",
    ),
    pytest.param(
        "id",
        ["rs1", "rs\r2"],
        "the id `rs\r2` holds a carriage return",
        id="carriage_return_in_an_id",
    ),
    pytest.param(
        "chrom", "", "segment 1: the chromosome is an empty text", id="empty_chrom"
    ),
    pytest.param(
        "chrom",
        ["chr1", ""],
        "segment 1, variant 2: the chromosome is an empty text",
        id="empty_chrom_of_one_variant",
    ),
    pytest.param(
        "alleles",
        [["H0", "."], ["H0", "H1", "H2"]],
        "segment 1, variant 1, chr1 100: one of its alleles is `.`",
        id="allele_dot",
    ),
    pytest.param(
        "alleles",
        [["H0", "H1"], ["H0", "H1", "H0"]],
        "segment 1, variant 2, chr1 5000: the allele `H0` is given twice",
        id="allele_twice",
    ),
]


@pytest.mark.parametrize(("key", "value", "words"), SEGMENTS_REFUSED)
def test_write_vars_from_segments_refuses_a_wrong_segment_and_leaves_no_file(
    tmp_path: Path, key: str, value: Any, words: str
) -> None:
    """Each wrong segment of "A call that fails", given after a segment that
    is right, so that the file had been made: a `ValueError` that names the
    segment, counted from 1, and the variant of it, and no file."""
    path = tmp_path / "refused.vars"

    with pytest.raises(ValueError) as refusal:
        write([changed(key, value), third_segment()], path)

    assert words in str(refusal.value)
    assert not path.exists()


def test_write_vars_from_segments_names_the_segment_counted_from_1(
    tmp_path: Path,
) -> None:
    """A wrong third segment after a right one and one of no variants is
    `segment 3`: the segments of no variants are counted."""
    path = tmp_path / "third.vars"
    wrong = third_segment()
    wrong["gts"] = numpy.array([[[0, 0], [0, 0], [0, 1]]])

    with pytest.raises(ValueError, match="^segment 3, variant 1, chr2 50: "):
        write([first_segment(), empty_segment(), wrong], path)

    assert not path.exists()


def test_write_vars_from_segments_refuses_an_item_that_is_not_a_mapping(
    tmp_path: Path,
) -> None:
    """A list where a dict was expected is a `TypeError` that names the
    segment."""
    path = tmp_path / "not_a_mapping.vars"

    with pytest.raises(TypeError, match="^segment 2 is a list"):
        write([first_segment(), [1, 2, 3]], path)

    assert not path.exists()


# Each argument of "A call that fails" that is refused before the file is
# made, and words its message holds.
ARGUMENTS_REFUSED = [
    pytest.param({"individuals": []}, "names no individual", id="no_individual"),
    pytest.param(
        {"individuals": ["i1", "i2", "i1"]}, "names `i1` twice", id="individual_twice"
    ),
    pytest.param({"ploidy": 0}, "the ploidy 0", id="ploidy_0"),
    pytest.param({"ploidy": 256}, "the ploidy 256", id="ploidy_256"),
    pytest.param(
        {"var_fields": {"END": "int8"}},
        "`var_fields` gives the field `END` the type 'int8'",
        id="type_not_of_the_five",
    ),
    pytest.param(
        {"gt_fields": {"GP": "not a type"}},
        "`gt_fields` gives the field `GP` the type 'not a type'",
        id="type_numpy_does_not_read",
    ),
    pytest.param(
        {"var_fields": {"1END": "uint64"}},
        "`1END` of the vars file has a name no field can have",
        id="name_against_the_rule",
    ),
    pytest.param(
        {"var_fields": {"qual": "float32"}},
        "`qual` of the vars file has the name of a column",
        id="name_of_a_column",
    ),
    pytest.param(
        {"var_fields": {"GP": "float32"}},
        "`GP` of the vars file is named twice",
        id="name_declared_twice",
    ),
    pytest.param(
        {"chrom_lengths": {"chr1": 0}},
        "`chrom_lengths` gives `chr1` the length 0",
        id="chrom_length_0",
    ),
    pytest.param(
        {"chrom_lengths": {"chr1": -100}},
        "`chrom_lengths` gives `chr1` the length -100",
        id="chrom_length_negative",
    ),
    pytest.param(
        {"chrom_lengths": {"chr1": True}},
        "`chrom_lengths` gives `chr1` the length True",
        id="chrom_length_true",
    ),
    pytest.param(
        {"chrom_lengths": {"chr 1": 100}},
        "`chrom_lengths` names `chr 1`, which holds a space",
        id="space_in_a_chrom_of_chrom_lengths",
    ),
    pytest.param(
        {"chrom_lengths": {"": 100}},
        "`chrom_lengths` names a chromosome whose name is an empty text",
        id="empty_chrom_of_chrom_lengths",
    ),
    pytest.param(
        {"individuals": ["i1", "i\t2", "i3"]},
        "`individuals` names `i\t2`, which holds a tab",
        id="tab_in_an_individual",
    ),
    pytest.param(
        {"individuals": ["i1", "", "i3"]},
        "`individuals` names an individual whose name is an empty text",
        id="empty_individual",
    ),
    pytest.param(
        {"gt_fields": {"GP": None}},
        "`gt_fields` gives the field `GP` the type None",
        id="type_none",
    ),
    pytest.param(
        {"num_vars_per_block": 0},
        "`num_vars_per_block` is 0",
        id="num_vars_per_block_0",
    ),
]


@pytest.mark.parametrize(("arguments", "words"), ARGUMENTS_REFUSED)
def test_write_vars_from_segments_refuses_a_wrong_argument_before_the_file_is_made(
    tmp_path: Path, arguments: dict[str, Any], words: str
) -> None:
    """Each wrong argument is a `ValueError` and nothing is taken from the
    segments: a generator that fails the test if it is asked for one."""
    path = tmp_path / "refused.vars"

    def never_asked() -> Iterator[dict[str, Any]]:
        raise AssertionError("a segment was asked for")
        yield first_segment()

    with pytest.raises(ValueError) as refusal:
        write(never_asked(), path, **arguments)

    assert words in str(refusal.value)
    assert not path.exists()


# Each argument of the wrong kind, a `TypeError` before the file is made,
# and words its message holds.
ARGUMENTS_OF_THE_WRONG_KIND = [
    pytest.param(
        {"individuals": {"i1", "i2", "i3"}},
        "`individuals` is a set",
        id="individuals_as_a_set",
    ),
    pytest.param(
        {"individuals": (name for name in INDIVIDUALS)},
        "`individuals` is a generator",
        id="individuals_as_a_generator",
    ),
    pytest.param(
        {"chrom_lengths": [("chr1", 100)]},
        "`chrom_lengths` is a list",
        id="chrom_lengths_as_a_list",
    ),
]


@pytest.mark.parametrize(("arguments", "words"), ARGUMENTS_OF_THE_WRONG_KIND)
def test_write_vars_from_segments_refuses_an_argument_of_the_wrong_kind(
    tmp_path: Path, arguments: dict[str, Any], words: str
) -> None:
    """`individuals` that is not a sequence, whose order a set would choose,
    and `chrom_lengths` that is not a mapping: each a `TypeError` that names
    the argument, and no file."""
    path = tmp_path / "refused.vars"

    with pytest.raises(TypeError) as refusal:
        write([first_segment()], path, **arguments)

    assert words in str(refusal.value)
    assert not path.exists()


def test_write_vars_from_segments_refuses_one_segment_given_as_the_segments(
    tmp_path: Path,
) -> None:
    """One dict where an iterable of dicts goes: a `TypeError` that says
    that one segment was given, and no file."""
    path = tmp_path / "one.vars"

    with pytest.raises(TypeError, match="one segment") as refusal:
        write(first_segment(), path)

    assert "[segment]" in str(refusal.value)
    assert not path.exists()


def test_write_vars_from_segments_takes_individuals_as_a_numpy_array(
    tmp_path: Path,
) -> None:
    """A numpy array of the names is a sequence of them, in its order."""
    path = tmp_path / "array.vars"

    write([first_segment()], path, individuals=numpy.array(INDIVIDUALS))

    assert open_vars(path).individuals == ("i1", "i2", "i3")


def test_write_vars_from_segments_takes_none_for_id_and_qual(tmp_path: Path) -> None:
    """`id` and `qual` given as `None` write what the keys left out write: no
    id and no quality."""
    path = tmp_path / "none.vars"
    segment = first_segment()
    segment["id"] = None
    segment["qual"] = None

    write([segment], path)

    _, table, _ = read_with_pyarrow(path)
    assert table.column("id").to_pylist() == [None, None]
    assert table.column("qual").to_pylist() == [None, None]


def test_write_vars_from_segments_refuses_a_path_that_a_file_is_at(
    tmp_path: Path,
) -> None:
    """A path that a file is already at is a `ValueError`, and the file is
    left as it was."""
    path = tmp_path / "taken.vars"
    path.write_bytes(b"someone's file")

    with pytest.raises(ValueError, match="a file is already there"):
        write([first_segment()], path)

    assert path.read_bytes() == b"someone's file"


@pytest.mark.parametrize(
    "segments",
    [
        pytest.param(
            [segment_at("chr1", 5000), segment_at("chr1", 100)],
            id="position_going_down",
        ),
        pytest.param(
            [
                segment_at("chr1", 100),
                segment_at("chr2", 200),
                segment_at("chr1", 300),
            ],
            id="chromosome_coming_back",
        ),
    ],
)
def test_write_vars_from_segments_refuses_unsorted_variants_unless_allowed(
    tmp_path: Path, segments: list[dict[str, Any]]
) -> None:
    """chr1 5000 then chr1 100, and chr1, chr2, chr1 with the positions going
    up: each a `ValueError` naming both variants and no file, and each
    written in its order with `allow_unsorted=True`."""
    path = tmp_path / "unsorted.vars"

    with pytest.raises(ValueError) as refusal:
        write(segments, path)

    message = str(refusal.value)
    last, before = segments[-1], segments[-2]
    assert (
        f"the {len(segments)}{'nd' if len(segments) == 2 else 'rd'} variant written, "
        f"{last['chrom']} {last['pos'][0]}, comes after the "
        f"{len(segments) - 1}{'st' if len(segments) == 2 else 'nd'}, "
        f"{before['chrom']} {before['pos'][0]}"
    ) in message
    assert not path.exists()

    written = write(segments, path, allow_unsorted=True)

    assert written.pass_stats.num_vars == len(segments)
    _, table, _ = read_with_pyarrow(path)
    assert table.column("pos").to_pylist() == [
        segment["pos"][0] for segment in segments
    ]


def test_write_vars_from_segments_raises_the_exception_of_the_generator(
    tmp_path: Path,
) -> None:
    """A generator that yields the first segment and then raises
    `RuntimeError("caller failed")`: that same exception comes out of the
    call, with the frame of the generator in its traceback, and no file."""
    path = tmp_path / "failed.vars"
    raised = RuntimeError("caller failed")

    def failing_caller() -> Iterator[dict[str, Any]]:
        yield first_segment()
        raise raised

    with pytest.raises(RuntimeError) as refusal:
        write(failing_caller(), path)

    assert refusal.value is raised
    frames = [entry.name for entry in refusal.traceback]
    assert "failing_caller" in frames
    assert not path.exists()


def test_write_vars_from_segments_raises_a_ctrl_c_between_segments(
    tmp_path: Path,
) -> None:
    """A generator that calls `_thread.interrupt_main()`, which is what a
    Ctrl-C does, after it yields its first segment: a `KeyboardInterrupt`
    out of the call, raised before the next segment is used, and no file."""
    path = tmp_path / "interrupted.vars"

    def interrupted_caller() -> Iterator[dict[str, Any]]:
        yield first_segment()
        _thread.interrupt_main()
        yield third_segment()

    with pytest.raises(KeyboardInterrupt):
        write(interrupted_caller(), path)

    assert not path.exists()


def test_write_vars_from_segments_writes_a_file_of_no_variants(
    tmp_path: Path,
) -> None:
    """No segment at all: a file of no batch, whose schema holds the extra
    columns and the fields, which pyarrow and `open_vars` read."""
    path = tmp_path / "none.vars"

    written = write(iter([]), path)

    assert written.pass_stats.num_vars == 0
    schema, table, footer = read_with_pyarrow(path)
    assert table.num_rows == 0
    assert json.loads(footer[b"popnei_batches"]) == []
    key = json.loads(schema.metadata[b"popnei"])
    assert key["var_fields"] == [["END", "uint64"]]
    assert open_vars(path).individuals == ("i1", "i2", "i3")


def worked_example(tmp_path: Path) -> Path:
    """The vars file of the three segments of the worked example."""
    path = tmp_path / "segments.vars"
    write(iter([first_segment(), empty_segment(), third_segment()]), path)
    return path


def test_variants_of_the_worked_example_name_its_fields_and_their_types(
    tmp_path: Path,
) -> None:
    """`var_fields` and `gt_fields` of the file are the fields it was
    written with, in their order, each with its numpy type."""
    variants = open_vars(worked_example(tmp_path))

    assert variants.var_fields == {"END": numpy.dtype("uint64")}
    assert list(variants.gt_fields.items()) == [
        ("GP", numpy.dtype("float32")),
        ("NODE_A", numpy.dtype("int32")),
    ]


def test_variants_of_a_vcf_have_no_field(reference_vcf_dir: Path) -> None:
    """A VCF has no extra field of either kind."""
    variants = open_vcf(reference_vcf_dir / "many.vcf")

    assert variants.var_fields == {}
    assert variants.gt_fields == {}


def test_iter_blocks_with_gp_gives_the_gp_that_went_in_and_no_other_field(
    tmp_path: Path,
) -> None:
    """`iter_blocks(fields=("chrom", "pos", "GP"), num_vars_per_block=3)`:
    one block whose `extra_fields` holds `GP` alone, 3 x 3, with NaN where
    the segment had it."""
    variants = open_vars(worked_example(tmp_path))

    blocks = list(
        variants.iter_blocks(fields=("chrom", "pos", "GP"), num_vars_per_block=3)
    )

    assert len(blocks) == 1
    (block,) = blocks
    assert block.chrom == ("chr1", "chr1", "chr2")
    assert block.extra_fields is not None
    assert list(block.extra_fields) == ["GP"]
    gp = block.extra_fields["GP"]
    assert gp.dtype == numpy.float32
    assert gp.shape == (3, 3)
    numpy.testing.assert_array_equal(
        gp,
        numpy.array(
            [[0.9, 0.8, NAN], [0.7, 0.99, 0.5], [1.0, 1.0, 0.6]], dtype=numpy.float32
        ),
    )


def test_iter_blocks_gives_each_field_with_its_shape_and_type_in_file_order(
    tmp_path: Path,
) -> None:
    """The three fields, named in another order than the file's: the dict
    is in the order of the file, `END` of one value for each variant and
    the two genotype fields of variants x individuals."""
    variants = open_vars(worked_example(tmp_path))

    (block,) = variants.iter_blocks(
        fields=("NODE_A", "END", "GP"), num_vars_per_block=3
    )

    assert block.extra_fields is not None
    assert list(block.extra_fields) == ["END", "GP", "NODE_A"]
    end = block.extra_fields["END"]
    assert end.dtype == numpy.uint64
    assert end.tolist() == [4999, 9999, 800]
    node_a = block.extra_fields["NODE_A"]
    assert node_a.dtype == numpy.int32
    assert node_a.tolist() == [[1, 2, -1], [3, 1, 2], [1, 1, 1]]
    assert block.chrom is None
    assert block.pos is None


def test_iter_blocks_of_the_default_fields_gives_no_extra_field(
    tmp_path: Path,
) -> None:
    """A field that is not named is not read: the default `fields` gives
    `extra_fields` None, as it does for a VCF."""
    variants = open_vars(worked_example(tmp_path))

    blocks = list(variants.iter_blocks(num_vars_per_block=2))

    assert [block.num_vars for block in blocks] == [2, 1]
    assert all(block.extra_fields is None for block in blocks)


def test_iter_blocks_gives_extra_fields_that_cannot_be_written(
    tmp_path: Path,
) -> None:
    """The arrays of the extra fields are read only, as `gts` is."""
    variants = open_vars(worked_example(tmp_path))

    (block,) = variants.iter_blocks(fields=("END", "GP"), num_vars_per_block=3)

    assert block.extra_fields is not None
    for name, array in block.extra_fields.items():
        assert not array.flags.writeable, name
        with pytest.raises(ValueError, match="read-only"):
            array[0] = 0


def test_iter_blocks_refuses_a_name_that_is_no_field_and_names_the_fields(
    tmp_path: Path,
) -> None:
    """A name that is neither one of the five nor a field of the file is a
    `ValueError` whose message names it and the fields the file has."""
    variants = open_vars(worked_example(tmp_path))

    with pytest.raises(ValueError) as refusal:
        variants.iter_blocks(fields=("chrom", "GQ"))

    message = str(refusal.value)
    for name in ["`GQ`", "`chrom`", "`qual`", "`END`", "`GP`", "`NODE_A`"]:
        assert name in message


def test_repr_of_a_block_names_its_extra_fields_after_its_columns(
    tmp_path: Path,
) -> None:
    """The `repr` of a block names the extra fields it holds after its
    columns, and none when it holds none."""
    variants = open_vars(worked_example(tmp_path))

    (with_fields,) = variants.iter_blocks(
        fields=("pos", "NODE_A", "GP"), num_vars_per_block=3
    )
    (without,) = variants.iter_blocks(fields=("pos",), num_vars_per_block=3)

    assert repr(with_fields) == (
        "<Block of 3 variants, gts int8 (3, 3, 2), columns chrom, pos, "
        "extra fields GP, NODE_A>"
    )
    assert repr(without) == (
        "<Block of 3 variants, gts int8 (3, 3, 2), columns chrom, pos>"
    )


GP_OF_I3_AND_I1 = [[NAN, 0.9], [0.5, 0.7], [0.6, 1.0]]


def test_filter_individuals_keeps_the_gp_of_the_kept_individuals_in_their_order(
    tmp_path: Path,
) -> None:
    """After `filter_individuals(["i3", "i1"])` the `GP` of the block is
    NaN 0.9, 0.5 0.7 and 0.6 1.0, and `write_vars` of the filtered
    `Variants` gives a file whose `GP` reads back as that, with `END` as
    it went in."""
    variants = open_vars(worked_example(tmp_path))
    variants.filter_individuals(["i3", "i1"])

    (block,) = variants.iter_blocks(fields=("chrom", "pos", "GP"), num_vars_per_block=3)
    assert block.extra_fields is not None
    numpy.testing.assert_array_equal(
        block.extra_fields["GP"], numpy.array(GP_OF_I3_AND_I1, dtype=numpy.float32)
    )

    filtered_path = tmp_path / "filtered.vars"
    write_vars(variants, filtered_path)
    filtered = open_vars(filtered_path)

    assert filtered.individuals == ("i3", "i1")
    assert filtered.var_fields == {"END": numpy.dtype("uint64")}
    assert list(filtered.gt_fields) == ["GP", "NODE_A"]
    (read_back,) = filtered.iter_blocks(fields=("GP", "END"), num_vars_per_block=3)
    assert read_back.extra_fields is not None
    numpy.testing.assert_array_equal(
        read_back.extra_fields["GP"],
        numpy.array(GP_OF_I3_AND_I1, dtype=numpy.float32),
    )
    assert read_back.extra_fields["END"].tolist() == [4999, 9999, 800]


def test_write_vcf_refuses_a_file_of_extra_fields_and_names_them(
    tmp_path: Path,
) -> None:
    """`write_vcf` of the filtered file, which keeps the three fields, is
    the `ValueError` that names `END`, `GP` and `NODE_A`, and leaves no
    VCF."""
    variants = open_vars(worked_example(tmp_path))
    variants.filter_individuals(["i3", "i1"])
    filtered_path = tmp_path / "filtered.vars"
    write_vars(variants, filtered_path)
    vcf_path = tmp_path / "filtered.vcf.gz"

    with pytest.raises(ValueError) as refusal:
        write_vcf(open_vars(filtered_path), vcf_path)

    message = str(refusal.value)
    for name in ["END", "GP", "NODE_A"]:
        assert name in message
    assert not vcf_path.exists()
