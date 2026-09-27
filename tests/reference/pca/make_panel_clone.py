"""It writes panel_clone.vcf.gz, the panel of the Kosman distances with a clone.

Run from the root of the repository, before pcoa_reference.R, which reads the
distances it implies:

    uv run python tests/reference/pca/make_panel_clone.py

The VCF is tests/reference/dists/panel.vcf.gz with a 201st individual, s200,
whose genotype at every variant is that of s000. Its Kosman distance to s000
is 0 and to every other individual that of s000, so the matrix of the
principal coordinates has two eigenvectors of the eigenvalue 0, which is the
case of the clones of "How it runs" of the principal coordinates in
docs/specs/pca.md. It is written with the timestamp of gzip fixed to 0, so that
running the script again gives the same bytes.
"""

import gzip
from pathlib import Path

HERE = Path(__file__).parent
PANEL = HERE.parent / "dists" / "panel.vcf.gz"
OUT = HERE / "panel_clone.vcf.gz"

# The first individual is the tenth column of a VCF line.
FIRST_INDIVIDUAL = 9

lines = []
with gzip.open(PANEL, "rt") as vcf:
    for line in vcf:
        line = line.rstrip("\n")
        if line.startswith("##"):
            lines.append(line)
            continue
        fields = line.split("\t")
        clone = "s200" if line.startswith("#CHROM") else fields[FIRST_INDIVIDUAL]
        lines.append("\t".join([*fields, clone]))

with open(OUT, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as out:
    out.write(("\n".join(lines) + "\n").encode())
