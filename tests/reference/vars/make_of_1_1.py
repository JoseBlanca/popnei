"""It writes of_1_1.vars, a vars file of format 1.1, which has no `passed`
column, for the tests of the filter of the variants that passed their FILTER.

Run from the root of the repository, after `uv run maturin develop`, with the
pyarrow of uv.lock, 25.0.1, which wrote the file that is committed:

    uv run python tests/reference/vars/make_of_1_1.py

popnei writes a vars file of format 1.2 from the four variants of cases.vcf of
tests/reference/vcf/, read with only_passed false, so that the one at 200,
which failed its filter q10, is there too. pyarrow then writes it again,
beside this script, with the same batches and the same keys, but without the
`passed` column and with format_version "1.1" in the `popnei` key: what a file
written before format 1.2 holds. popnei opens it, and the filter of the
variants that passed refuses it at the first block, as "The filter of the
variants that passed their FILTER" of docs/specs/filters.md says, rather than
taking the failed variant at 200 as passed.
"""

import tempfile
from pathlib import Path

from popnei import open_vcf, write_vars
from pyarrow import ipc

HERE = Path(__file__).parent
PATH = HERE / "of_1_1.vars"
CASES_VCF = HERE.parent / "vcf" / "cases.vcf"

POPNEI_KEY = b"popnei"
VERSION_1_2 = b'"format_version":"1.2"'
VERSION_1_1 = b'"format_version":"1.1"'


def write(path):
    with tempfile.TemporaryDirectory() as tmp_dir:
        written = Path(tmp_dir) / "cases.vars"
        write_vars(open_vcf(CASES_VCF, only_passed=False), written)
        with ipc.open_file(written) as reader:
            footer = reader.metadata
            schema_metadata = dict(reader.schema.metadata)
            table = reader.read_all()
    assert VERSION_1_2 in schema_metadata[POPNEI_KEY], schema_metadata
    schema_metadata[POPNEI_KEY] = schema_metadata[POPNEI_KEY].replace(
        VERSION_1_2, VERSION_1_1
    )
    table = table.drop_columns(["passed"]).replace_schema_metadata(schema_metadata)
    with ipc.new_file(path, table.schema, metadata=footer) as writer:
        for batch in table.to_batches():
            writer.write_batch(batch)


def check(path):
    """It reads the file back: the four variants, six columns and no
    `passed`, and the version 1.1."""
    with ipc.open_file(path) as reader:
        assert reader.schema.names == ["chrom", "pos", "id", "alleles", "qual", "gts"]
        assert VERSION_1_1 in reader.schema.metadata[POPNEI_KEY]
        table = reader.read_all()
    assert table.column("pos").to_pylist() == [100, 200, 300, 400]


if __name__ == "__main__":
    write(PATH)
    check(PATH)
