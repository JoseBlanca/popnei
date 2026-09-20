"""It writes the reference VCFs of the VCF reader and what bcftools reads in them.

Run once, with bcftools and bgzip 1.24 in the PATH:

    python3 make_reference.py

It writes, beside itself:

- cases.vcf, four variants of three diploid individuals written by hand, which
  pyNei and popnei read the same, and cases.vcf.gz, the same bgzipped.
- differences.vcf, two variants that pyNei reads in another way, or refuses,
  and popnei reads as bcftools does.
- many.vcf, 500 variants of 50 diploid individuals drawn with a fixed seed,
  and many.vcf.gz. One variant in 20 failed its FILTER, q10, and one in 20
  has no FILTER, a dot. It is 115 KB, more than the 64 KB that one block of bgzip
  holds, so the gzipped file has several gzip members one after another.
- cases.bcftools.tsv, differences.bcftools.tsv and many.bcftools.tsv, what
  `bcftools query` prints for each: chrom, pos, id, ref, alt, qual, filter and
  the GT of every individual.

docs/specs/io_vcf.md says how the tests use them.
"""

import random
import subprocess
from pathlib import Path

HERE = Path(__file__).parent
QUERY_FORMAT = r"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n"

HEADER = """##fileformat=VCFv4.4
##contig=<ID=chr1>
##contig=<ID=chr2>
##FILTER=<ID=q10,Description="Quality below 10">
##ALT=<ID=DEL,Description="Deletion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
"""
COLUMNS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

CASES = [
    # every field given, the three diploid genotypes of a biallelic variant
    "chr1 100 rs1 A T 29.5 PASS . GT:DP 0/0:3 0/1:4 1/1:5",
    # no id, no qual, a FILTER that failed, a missing genotype, a phased one
    # and a half called one
    "chr1 200 . A T . q10 . GT:DP ./.:. 0|1:3 .|0:2",
    # two alternative alleles
    "chr1 300 . A G,T 67 PASS . GT 1/2 2|1 2/2",
    # no alternative allele
    "chr1 400 . T . 47 PASS . GT 0/0 0/0 0/0",
]

DIFFERENCES = [
    # another chromosome, an individual that drops its last field, and a missing
    # genotype written as one dot
    "chr2 50 ms1 GTC G,GTCT 50 PASS . GT:DP 0/1:3 0/2 .",
    # a symbolic allele, the overlapping deletion allele, and genotypes that
    # start with their separator, as VCF 4.4 allows
    "chr2 60 . A <DEL>,* . PASS . GT /0/1 |2|2 0/0",
]


def write_by_hand(name, variants):
    lines = [HEADER, "\t".join(COLUMNS + ["ind1", "ind2", "ind3"]) + "\n"]
    lines += ["\t".join(variant.split(" ")) + "\n" for variant in variants]
    (HERE / f"{name}.vcf").write_text("".join(lines))


def write_many(num_vars=500, num_individuals=50, seed=42):
    rng = random.Random(seed)
    individuals = [f"ind{idx:02d}" for idx in range(num_individuals)]
    lines = [HEADER, "\t".join(COLUMNS + individuals) + "\n"]
    for var_idx in range(num_vars):
        chrom = "chr1" if var_idx < num_vars // 2 else "chr2"
        pos = 1000 + 37 * var_idx
        num_alts = 2 if rng.random() < 0.1 else 1
        alt = "G,T"[: 2 * num_alts - 1]
        freqs = [rng.random() for _ in range(num_alts + 1)]
        gts = []
        for _ in individuals:
            alleles = []
            for _ in range(2):
                pick = rng.random() * sum(freqs)
                allele = 0
                while pick > freqs[allele]:
                    pick -= freqs[allele]
                    allele += 1
                alleles.append(str(allele))
            draw = rng.random()
            if draw < 0.05:
                alleles = [".", "."]
            elif draw < 0.06:
                alleles[0] = "."
            sep = "|" if rng.random() < 0.3 else "/"
            gts.append(sep.join(alleles))
        # no random draw here, so that the genotypes do not depend on it
        filter_ = {7: "q10", 13: "."}.get(var_idx % 20, "PASS")
        fields = [chrom, str(pos), ".", "A", alt, ".", filter_, ".", "GT"] + gts
        lines.append("\t".join(fields) + "\n")
    (HERE / "many.vcf").write_text("".join(lines))


def run_tools(name):
    vcf = HERE / f"{name}.vcf"
    with (HERE / f"{name}.vcf.gz").open("wb") as gz:
        subprocess.run(["bgzip", "-c", str(vcf)], stdout=gz, check=True)
    with (HERE / f"{name}.bcftools.tsv").open("wb") as tsv:
        subprocess.run(
            ["bcftools", "query", "-f", QUERY_FORMAT, str(vcf)], stdout=tsv, check=True
        )


if __name__ == "__main__":
    write_by_hand("cases", CASES)
    write_by_hand("differences", DIFFERENCES)
    write_many()
    for name in ("cases", "differences", "many"):
        run_tools(name)
