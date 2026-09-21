"""It writes zstd.vars, the vars file that no build of popnei can write.

Run from the root of the repository, with the pyarrow of uv.lock, 25.0.1,
which wrote the file that is committed:

    uv run python tests/reference/vars/make_reference.py

It writes, beside itself, zstd.vars: the four variants of cases.vcf of
tests/reference/vcf/, three diploid individuals, in one record batch, with the
six columns, the types, the nulls and the two keys that popnei's writer gives
such a file, and with the buffers of the batch compressed with zstd. popnei
writes lz4 and reads lz4 and no compression, because no build of it carries
the zstd crate, so its reader opens this file, takes the two keys from it and
refuses it at its first block. That is the test of "How it is verified" of the
reader in docs/specs/io_vars.md, and "The compression" there says why zstd is
out.

The `popnei` key says `num_vars_per_block` 4, the four variants of the one
batch, whatever number popnei's own writer would put there.

The variants are literals here, from the cases.vcf table of
docs/specs/io_vcf.md, so nothing of popnei is run to write the file. Writing
it again gives the same bytes.
"""

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.ipc as ipc

HERE = Path(__file__).parent
PATH = HERE / "zstd.vars"

INDIVIDUALS = ["ind1", "ind2", "ind3"]
PLOIDY = 2
NUM_VARS_PER_BLOCK = 4

# The four variants that tests/reference/vcf/cases.vcf holds, read with
# only_passed false so that the second one, which failed its filter q10, is
# there too. Each is the chromosome, the position, the id or None for a
# variant with none, the alleles with the reference first, the quality or None
# for a variant with none, and the genotypes of the three individuals one
# after another, two alleles each, -1 for an allele that was not called.
VARIANTS = [
    ("chr1", 100, "rs1", ["A", "T"], 29.5, [0, 0, 0, 1, 1, 1]),
    ("chr1", 200, None, ["A", "T"], None, [-1, -1, 0, 1, -1, 0]),
    ("chr1", 300, None, ["A", "G", "T"], 67.0, [1, 2, 2, 1, 2, 2]),
    ("chr1", 400, None, ["T"], 47.0, [0, 0, 0, 0, 0, 0]),
]

# The value of the key `popnei` of the schema, which holds what is known
# before the first variant. The field inside each of the two lists is the one
# pyarrow writes for any list, named item and allowed to be null. popnei
# writes that same name and says that the field holds no null, which is what
# is true of an allele and of a genotype; its reader compares what a list
# holds and not that field, so it reads this file as it reads its own.
SCHEMA = pa.schema(
    [
        pa.field("chrom", pa.string(), nullable=False),
        pa.field("pos", pa.uint64(), nullable=False),
        pa.field("id", pa.string(), nullable=True),
        pa.field("alleles", pa.list_(pa.string()), nullable=False),
        pa.field("qual", pa.float32(), nullable=True),
        pa.field(
            "gts", pa.list_(pa.int8(), len(INDIVIDUALS) * PLOIDY), nullable=False
        ),
    ],
    metadata={
        "popnei": json.dumps(
            {
                "format_version": "1.0",
                "individuals": INDIVIDUALS,
                "ploidy": PLOIDY,
                "num_vars_per_block": NUM_VARS_PER_BLOCK,
            },
            separators=(",", ":"),
        )
    },
)


def batch_info(variants):
    """The entry of `popnei_batches` of a batch that holds `variants`.

    One region for each chromosome with a variant in the batch, in the order
    in which the chromosomes first appear, with the smallest and the largest
    position of its variants.
    """
    regions = {}
    for chrom, pos, _id, _alleles, _qual, _gts in variants:
        if chrom in regions:
            region = regions[chrom]
            region["min_pos"] = min(region["min_pos"], pos)
            region["max_pos"] = max(region["max_pos"], pos)
        else:
            regions[chrom] = {"chrom": chrom, "min_pos": pos, "max_pos": pos}
    return {"num_vars": len(variants), "regions": list(regions.values())}


def record_batch(variants):
    columns = [
        pa.array([variant[index] for variant in variants], field.type)
        for index, field in enumerate(SCHEMA)
    ]
    return pa.record_batch(columns, schema=SCHEMA)


def footer(variants):
    """The metadata of the footer, one entry of `popnei_batches` per batch."""
    entries = json.dumps([batch_info(variants)], separators=(",", ":"))
    return {"popnei_batches": entries}


def write(path):
    options = ipc.IpcWriteOptions(compression="zstd")
    with pa.OSFile(str(path), "wb") as sink:
        metadata = footer(VARIANTS)
        with ipc.new_file(sink, SCHEMA, options=options, metadata=metadata) as writer:
            writer.write_batch(record_batch(VARIANTS))


def check(path):
    """It reads the file back and compares it with the literals above."""
    opened = ipc.open_file(path)
    assert opened.schema == SCHEMA, opened.schema
    assert opened.schema.metadata == SCHEMA.metadata, opened.schema.metadata
    assert opened.num_record_batches == 1, opened.num_record_batches
    expected = {key.encode(): value.encode() for key, value in footer(VARIANTS).items()}
    assert opened.metadata == expected, opened.metadata
    batch = opened.get_batch(0)
    assert batch.num_rows == len(VARIANTS), batch.num_rows
    for index, variant in enumerate(VARIANTS):
        row = batch.slice(index, 1).to_pylist()[0]
        assert row["chrom"] == variant[0], row
        assert row["pos"] == variant[1], row
        assert row["id"] == variant[2], row
        assert row["alleles"] == variant[3], row
        assert row["qual"] == variant[4], row
        assert row["gts"] == variant[5], row
    # The four bytes that every zstd frame starts with, which are in the body
    # of the batch and in no file that arrow wrote with lz4 or with no
    # compression.
    assert b"\x28\xb5\x2f\xfd" in Path(path).read_bytes()


if __name__ == "__main__":
    write(PATH)
    check(PATH)
