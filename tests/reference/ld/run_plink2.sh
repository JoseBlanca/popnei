#!/usr/bin/env bash
# Writes again the ten files of tests/reference/ld/ that a program made,
# and compares each with the copy stored there: ld.vcf.gz, the linkage
# disequilibrium dataset; the r2 matrix plink2 gives for it and for the
# worked example, with the identifiers of the rows of each matrix beside
# it; the dosages pyNei gives for tests/reference/vcf/many.vcf; the
# variants that the filter by linkage disequilibrium keeps of ld.vcf at the
# four settings of the table of docs/specs/filters.md, worked out from
# plink2's matrix; the three properties that item asks of the set kept,
# each as the number of variants or pairs that break it, worked out over the
# variants read back from the file of kept variants, with each property also
# run against a set built to break it, whose number is not 0 and which the
# comparison below covers like the rest; and the bins of r2 against distance
# of the three populations of "How it is verified" of docs/specs/ld.md, which
# docs/reports/ld-method/bins.py prints from the r2 plink2 gives for the
# individuals and the variants of each of them; and the curve fitted to the
# pairs of those same three populations, which docs/reports/ld-method/decay.R
# fits in R to the pairs docs/reports/ld-method/decay.py groups by their
# exact distance.
#
# Run it from the root of the repository with one argument, a directory it
# works in, which it creates and which has to be empty or not exist:
#
#     tests/reference/ld/run_plink2.sh /tmp/ld_work
#
# It needs plink2 v2.0.0-a.7.7, the version the numbers of docs/specs/ld.md
# and of the item "The filter by linkage disequilibrium" of
# docs/specs/filters.md were taken with on 22 September 2026; uv for the
# Python that make_reference.py, bins.py and decay.py run under, the first
# importing the pyNei of pyproject.toml and the other two numpy; and R 4.6.1,
# run as Rscript, which the curve of the fall-off is fitted in and which
# docs/objectives.md names among the reference programs.
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

# bins.py reads ld.vcf from LD_WORK, runs plink2 once for each of the three
# populations of the item "LD against distance, per population" of
# docs/specs/ld.md under the prefix x, and prints the ten bins of each with
# how many pairs they hold, their mean r2 and its standard deviation. What it
# prints is what the cargo tests of calc_ld_and_dist assert: the standard
# deviations of pop_a and of pop_b are in this file alone, where the table of
# the spec leaves them out.
LD_WORK="$work" uv run python "$here/../../../docs/reports/ld-method/bins.py" \
    > "$work/ld.bins.txt"

# decay.py runs plink2 once more for each of those three populations and
# writes <pop>.decay.tsv, the pairs of the population grouped by their exact
# distance, which is what the curve is fitted to and what the bins above are
# the same pairs gathered into ten. decay.R reads the three files and fits
# the curve of "The curve that is fitted" of docs/specs/ld.md to each,
# with R's optimize and again with R's nls, and prints the three rows of the
# table of "How it is verified" that the cargo tests of the curve assert,
# with the numbers the text around that table quotes. The three lines
# decay.py itself prints are not stored: they are the distances and the
# pairs of each population, which decay.R prints again beside its rows.
LD_WORK="$work" uv run python "$here/../../../docs/reports/ld-method/decay.py" \
    > "$work/ld.decay.pairs.txt"
LD_WORK="$work" Rscript "$here/../../../docs/reports/ld-method/decay.R" \
    > "$work/ld.decay.txt"

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
            many.pynei.dosages.tsv ld.filtered.tsv \
            ld.filter.properties.txt ld.bins.txt ld.decay.txt; do
    if ! cmp -s "$here/$name" "$work/$name"; then
        echo "$name"
        differed=1
    fi
done

exit "$differed"
