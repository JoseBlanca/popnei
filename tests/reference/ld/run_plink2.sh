#!/usr/bin/env bash
# Writes again the six files of tests/reference/ld/ that a program made,
# and compares each with the copy stored there: ld.vcf.gz, the linkage
# disequilibrium dataset; the r2 matrix plink2 gives for it and for the
# worked example, with the identifiers of the rows of each matrix beside
# it; and the dosages pyNei gives for tests/reference/vcf/many.vcf.
#
# Run it from the root of the repository with one argument, a directory it
# works in, which it creates and which has to be empty or not exist:
#
#     tests/reference/ld/run_plink2.sh /tmp/ld_work
#
# It needs plink2 v2.0.0-a.7.7, the version the numbers of docs/specs/ld.md
# and of the item "The filter by linkage disequilibrium" of
# docs/specs/filters.md were taken with on 22 September 2026, and uv for the
# Python that make_reference.py runs under, which imports the pyNei of
# pyproject.toml.
#
# It writes nothing into the repository. plink2 says what it is doing as it
# goes; after that the script prints the name of every file it made that
# differs from the one stored, so a run where everything agreed prints
# nothing of its own and ends with status 0, and one where something
# differed ends with status 1. A file that differs is copied over by hand
# once the difference is understood: every literal of docs/specs/ld.md and
# of the filter item of docs/specs/filters.md comes out of these bytes.

set -euo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

if [ $# -ne 1 ]; then
    echo "usage: $0 <work directory>" >&2
    exit 2
fi
work=$1
mkdir -p "$work"
if [ -n "$(ls -A "$work")" ]; then
    echo "$work is not empty" >&2
    exit 2
fi

# make_reference.py takes the directory it writes into from LD_WORK and
# writes ld.vcf uncompressed, which is what plink2 reads. Its seed is
# numpy.random.default_rng(29) and every number of both specs depends on the
# order in which it asks that generator, so the file is the same on every
# machine and the script is never tidied.
LD_WORK="$work" uv run python "$here/make_reference.py"

# gzip writes the time of the run into its header, which would give a
# different file at every run, so the time is fixed to 0 and the name of the
# file is left out, as tests/reference/dists/make_reference.py does for its
# own VCFs.
uv run python - "$work" <<'PY'
import gzip
import sys
from pathlib import Path

work = Path(sys.argv[1])
with (
    (work / "ld.vcf.gz").open("wb") as raw,
    gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed,
):
    compressed.write((work / "ld.vcf").read_bytes())
PY

cp "$here/example.vcf" "$work/example.vcf"

# --r2-unphased is the squared correlation between the dosages of two
# variants with the individuals missing at either one left out, and
# "square bin" writes the whole matrix as float64, row after row, where the
# text plink2 writes by default has six digits. It goes into
# <prefix>.unphased.vcor2.bin, and the identifiers of the rows, in their
# order, into <prefix>.unphased.vcor2.bin.vars.
plink2 --vcf "$work/ld.vcf" --double-id --allow-extra-chr \
       --r2-unphased square bin --out "$work/ld"
plink2 --vcf "$work/example.vcf" --double-id --allow-extra-chr \
       --r2-unphased square bin --out "$work/example"

# plink2 also writes a .log, which carries the time of the run and the paths
# of the machine, so it is not stored.
differed=0

# The check of the dataset is on the text and not on the bytes of the gzip:
# what the tests read is what comes out of it.
if ! gzip -dc "$here/ld.vcf.gz" | diff - "$work/ld.vcf" > /dev/null; then
    echo "ld.vcf.gz"
    differed=1
fi
for name in ld.unphased.vcor2.bin ld.unphased.vcor2.bin.vars \
            example.unphased.vcor2.bin example.unphased.vcor2.bin.vars \
            many.pynei.dosages.tsv; do
    if ! cmp -s "$here/$name" "$work/$name"; then
        echo "$name"
        differed=1
    fi
done

exit "$differed"
