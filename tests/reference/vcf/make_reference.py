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
  has no FILTER, a dot; one in three has no ID and one in five no QUAL, so that
  the comparison with pyNei covers those two columns on 500 variants. It is
  more than the 64 KB that one block of bgzip
  holds, so the gzipped file has several gzip members one after another.
- cases.bcftools.tsv, differences.bcftools.tsv and many.bcftools.tsv, what
  `bcftools query` prints for each: chrom, pos, id, ref, alt, qual, filter and
  the GT of every individual.
- write.vcf, the six lines of "How it is verified" of the VCF writer, with
  INFO values, a second value for each individual, phase, a FILTER that
  failed, a deletion and the lengths of the two chromosomes.
- write.c_a.bcftools.vcf, what `bcftools view -I -s c,a` followed by
  `bcftools annotate -x INFO/AC,INFO/AN` writes of write.vcf: the file the
  writer gives for the filter of individuals that keeps c and a, but for the
  `##FILTER=<ID=PASS,...>` line that bcftools adds and popnei does not.
- write.passed.bcftools.tsv, what `bcftools query` prints of the lines of
  write.vcf whose FILTER is PASS or a dot: chrom, pos, id, ref, alt, qual and
  the GT of every individual, the rows of the file the writer gives from the
  vars file of write.vcf, but for the phase, which a vars file does not keep.
- many.missing_0.04.bcftools.vcf, the 215 data lines that
  `bcftools view -H -i 'F_MISSING<=0.04'` prints of many.vcf.

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
        # no random draw in the three that follow, so that the genotypes do
        # not depend on them: one variant in 20 failed its filter and one in
        # 20 has none, one in three has no id, one in five no quality, and
        # half of the qualities that are there have a decimal
        filter_ = {7: "q10", 13: "."}.get(var_idx % 20, "PASS")
        id_ = "." if var_idx % 3 == 0 else f"var{var_idx:03d}"
        if var_idx % 5 == 0:
            qual = "."
        elif var_idx % 2 == 0:
            qual = str(20 + var_idx % 60)
        else:
            qual = f"{20 + var_idx % 60}.5"
        fields = [chrom, str(pos), id_, "A", alt, qual, filter_, ".", "GT"] + gts
        lines.append("\t".join(fields) + "\n")
    (HERE / "many.vcf").write_text("".join(lines))


WRITE_VCF = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=2000>
##contig=<ID=chr2,length=1500>
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FILTER=<ID=q10,Description="Quality below 10">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT a b c
chr1 100 rs1 A T 29.5 PASS AC=4;AN=6;DP=12 GT:DP 0/1:4 0|1:5 1/1:3
chr1 250 . AT A . q10 AC=1;AN=4;DP=5 GT:DP ./.:0 0/1:3 0/0:2
chr1 1000 rs3 G C,T 50 PASS AC=2,1;AN=6;DP=20 GT:DP 1/2:7 0/1:6 0/0:7
chr1 1001 . C . 12 . DP=9 GT 0/0 0/0 0/0
chr2 1 . T G 40 PASS AC=2;AN=5;DP=8 GT:DP 1|1:4 0/.:2 0/0:2
chr2 1500 rs6 A G 33 PASS AC=1;AN=6 GT 0/0 1/0 0/0
"""
WRITE_QUERY_FORMAT = r"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL[\t%GT]\n"


def write_the_writer_cases():
    """write.vcf and what bcftools 1.24 writes and prints of it and of many.vcf
    for the tests of the VCF writer."""
    lines = []
    for line in WRITE_VCF.splitlines():
        # the meta lines have their blanks, the others are cut at them
        lines.append(line if line.startswith("##") else "\t".join(line.split(" ")))
    write_vcf = HERE / "write.vcf"
    write_vcf.write_text("\n".join(lines) + "\n")
    view = subprocess.run(
        ["bcftools", "view", "--no-version", "-I", "-s", "c,a", str(write_vcf)],
        capture_output=True,
        check=True,
    )
    with (HERE / "write.c_a.bcftools.vcf").open("wb") as out:
        subprocess.run(
            ["bcftools", "annotate", "--no-version", "-x", "INFO/AC,INFO/AN"],
            input=view.stdout,
            stdout=out,
            check=True,
        )
    with (HERE / "write.passed.bcftools.tsv").open("wb") as tsv:
        subprocess.run(
            [
                "bcftools",
                "query",
                "-i",
                'FILTER="PASS" || FILTER="."',
                "-f",
                WRITE_QUERY_FORMAT,
                str(write_vcf),
            ],
            stdout=tsv,
            check=True,
        )
    with (HERE / "many.missing_0.04.bcftools.vcf").open("wb") as out:
        subprocess.run(
            ["bcftools", "view", "-H", "-i", "F_MISSING<=0.04", str(HERE / "many.vcf")],
            stdout=out,
            check=True,
        )


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
    write_the_writer_cases()
