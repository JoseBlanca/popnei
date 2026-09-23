"""It writes what plink2 and bcftools give for the statistics of docs/specs/stats.md.

Run from the root of the repository, with plink2 v2.0.0-a.7.7 and bcftools
1.24 in the PATH, the versions that wrote the files that are committed and
the only ones the script runs with, and with pyNei importable, which the
project's environment has as a development dependency:

    uv run python tests/reference/stats/make_reference.py

It reads two datasets. The panel is pyNei's test/gwas_reference/sim_missing.vars,
200 diploid individuals named s000 to s199, 1200 biallelic variants, 3 in
100 genotypes missing whole, in three populations, p0, p1 and p2, of 48, 68
and 84 individuals, which are the `pop` column of phenotypes.csv beside it;
it is a vars file of pyNei, which popnei does not read, so the script writes
it here as `panel.vcf.gz`, with `panel_pops.txt` for plink2, IID and popcat,
and `panel_pops_bcftools.txt`, name and population. The other is
tests/reference/vcf/many.vcf, 500 variants of 50 diploid individuals, ind00
to ind49, with half called genotypes and variants of three alleles, with two
populations, popA of the first 20 individuals and popB of the other 30, in
`many_pops.txt`.

It runs, on the panel, the plink2 commands of "How it is verified" of the
observed heterozygosity, the maf, the expected heterozygosity and the per
individual statistics, once per population and once over all, and writes
their reports beside itself:

    panel.p0.hardy  panel.p0.afreq   and the same for p1 and p2
    panel.hardy     panel.afreq      panel.vmiss
    panel.smiss     panel.scount     panel.het

and on many.vcf the bcftools commands of the maf and of the observed
heterozygosity, whose counts go into `many.counts.tsv`, one line per variant
with the position, AN and AC of popA, of popB and of all, and the
heterozygous and the called genotypes, and the plink2 command of the per
individual statistics with a half called genotype read as missing:

    many.counts.tsv  many.smiss  many.scount

The script checks what it got against the literals the spec gives, and stops
at the first one that differs.
"""

import csv
import gzip
import subprocess
import sys
from pathlib import Path

import numpy

from pynei import load_vars

HERE = Path(__file__).parent
ROOT = HERE.parent.parent.parent
MANY_VCF = ROOT / "tests" / "reference" / "vcf" / "many.vcf"
PLINK2_VERSION = "v2.0.0-a.7.7"
BCFTOOLS_VERSION = "1.24"


def pynei_reference_dir():
    """The test/gwas_reference of the pyNei that is installed."""
    import pynei

    candidates = [
        Path(pynei.__file__).parent.parent.parent / "test" / "gwas_reference",
        Path.home() / "devel" / "pynei" / "test" / "gwas_reference",
    ]
    for candidate in candidates:
        if (candidate / "sim_missing.vars").exists():
            return candidate
    raise SystemExit(
        "the panel, test/gwas_reference/sim_missing.vars of pyNei, was not "
        f"found; looked in {[str(c) for c in candidates]}"
    )


def check_versions():
    plink2 = subprocess.run(["plink2", "--version"], capture_output=True, text=True)
    if PLINK2_VERSION not in plink2.stdout:
        raise SystemExit(f"plink2 {PLINK2_VERSION} is needed; found: {plink2.stdout.strip()}")
    bcftools = subprocess.run(["bcftools", "--version"], capture_output=True, text=True)
    if bcftools.stdout.splitlines()[0] != f"bcftools {BCFTOOLS_VERSION}":
        raise SystemExit(f"bcftools {BCFTOOLS_VERSION} is needed; found: {bcftools.stdout.splitlines()[0]}")


def write_panel(ref_dir):
    """The panel as a gzipped VCF, and its populations in the two forms."""
    variants = load_vars(str(ref_dir / "sim_missing.vars"))
    individuals = list(variants.samples)
    with open(ref_dir / "phenotypes.csv") as f:
        pop_of = {row["IID"]: f"p{row['pop']}" for row in csv.DictReader(f)}
    with open(HERE / "panel_pops.txt", "w") as f:
        f.write("IID\tpopcat\n")
        for name in individuals:
            f.write(f"{name}\t{pop_of[name]}\n")
    with open(HERE / "panel_pops_bcftools.txt", "w") as f:
        for name in individuals:
            f.write(f"{name}\t{pop_of[name]}\n")
    with gzip.open(HERE / "panel.vcf.gz", "wt") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=1>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(individuals) + "\n")
        pos = 0
        for chunk in variants.iter_vars_chunks():
            gts = chunk.gts.gt_values
            for var in range(chunk.num_vars):
                pos += 1
                row = gts[var]
                calls = ["./." if numpy.any(g < 0) else f"{g[0]}/{g[1]}" for g in row]
                f.write(f"1\t{pos}\tvar{pos - 1:04d}\tA\tC\t.\tPASS\t.\tGT\t" + "\t".join(calls) + "\n")
    return individuals


def run(command):
    subprocess.run(command, check=True, capture_output=True, text=True)


def run_plink2_on_panel():
    reports = "--hardy --freq cols=+pos,+reffreq,+nobs --missing --sample-counts cols=+hom,+het,+missing --het".split()
    run(["plink2", "--vcf", str(HERE / "panel.vcf.gz"), "--pheno", str(HERE / "panel_pops.txt"),
         "--loop-cats", "popcat", "--hardy", "--freq", "cols=+pos,+reffreq,+nobs",
         "--nonfounders", "--out", str(HERE / "panel")])
    run(["plink2", "--vcf", str(HERE / "panel.vcf.gz"), *reports, "--nonfounders", "--out", str(HERE / "panel")])
    for log in HERE.glob("panel*.log"):
        log.unlink()


def write_many_pops():
    with open(MANY_VCF) as f:
        for line in f:
            if line.startswith("#CHROM"):
                individuals = line.rstrip("\n").split("\t")[9:]
                break
    with open(HERE / "many_pops.txt", "w") as f:
        for idx, name in enumerate(individuals):
            f.write(f"{name}\t{'popA' if idx < 20 else 'popB'}\n")


def run_bcftools_on_many():
    tags = subprocess.run(
        ["bcftools", "+fill-tags", str(MANY_VCF), "--", "-S", str(HERE / "many_pops.txt"), "-t", "AN,AC"],
        check=True, capture_output=True, text=True,
    ).stdout
    counts = subprocess.run(
        ["bcftools", "query", "-f", "%POS\t%AN_popA\t%AC_popA\t%AN_popB\t%AC_popB\t%AN\t%AC\n"],
        input=tags, check=True, capture_output=True, text=True,
    ).stdout.splitlines()
    hets = subprocess.run(
        ["bcftools", "+fill-tags", str(MANY_VCF), "--", "-t",
         'NHET=N_PASS(GT="het"),NCALLED=N_PASS(GT!="mis")'],
        check=True, capture_output=True, text=True,
    ).stdout
    hets = subprocess.run(
        ["bcftools", "query", "-f", "%NHET\t%NCALLED\n"], input=hets, check=True, capture_output=True, text=True,
    ).stdout.splitlines()
    assert len(counts) == len(hets) == 500, (len(counts), len(hets))
    with open(HERE / "many.counts.tsv", "w") as f:
        f.write("pos\tAN_popA\tAC_popA\tAN_popB\tAC_popB\tAN\tAC\tNHET\tNCALLED\n")
        for count, het in zip(counts, hets):
            f.write(f"{count}\t{het}\n")


def run_plink2_on_many():
    run(["plink2", "--vcf", str(MANY_VCF), "--vcf-half-call", "m", "--missing",
         "--sample-counts", "cols=+hom,+het,+missing", "--nonfounders", "--out", str(HERE / "many")])
    for name in ("many.log", "many.vmiss"):
        (HERE / name).unlink(missing_ok=True)


def read_table(path):
    with open(path) as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    header = [name.lstrip("#") for name in rows[0]]
    return [dict(zip(header, row)) for row in rows[1:]]


def check():
    """The literals of docs/specs/stats.md, as the spec prints them."""
    first = {pop: read_table(HERE / f"panel.{pop}.hardy")[0] for pop in ("p0", "p1", "p2")}
    expected = {
        "p0": ("39", "8", "1", "0.166667", "0.186632"),
        "p1": ("15", "34", "18", "0.507463", "0.498998"),
        "p2": ("81", "2", "0", "0.0240964", "0.0238061"),
    }
    for pop, values in expected.items():
        got = tuple(first[pop][col] for col in ("HOM_A1_CT", "HET_A1_CT", "TWO_AX_CT", "O(HET_A1)", "E(HET_A1)"))
        assert got == values, (pop, got)
    all_first = read_table(HERE / "panel.hardy")[0]
    assert (all_first["HOM_A1_CT"], all_first["HET_A1_CT"], all_first["TWO_AX_CT"], all_first["O(HET_A1)"], all_first["E(HET_A1)"]) == ("135", "44", "19", "0.222222", "0.328385"), all_first
    freqs = {pop: read_table(HERE / f"panel.{pop}.afreq")[0] for pop in ("p0", "p1", "p2")}
    assert (freqs["p0"]["OBS_CT"], freqs["p0"]["REF_FREQ"]) == ("96", "0.895833"), freqs["p0"]
    assert (freqs["p1"]["OBS_CT"], freqs["p1"]["ALT_FREQS"]) == ("134", "0.522388"), freqs["p1"]
    assert (freqs["p2"]["OBS_CT"], freqs["p2"]["REF_FREQ"]) == ("166", "0.987952"), freqs["p2"]
    all_freq = read_table(HERE / "panel.afreq")[0]
    assert (all_freq["OBS_CT"], all_freq["REF_FREQ"]) == ("396", "0.792929"), all_freq
    # the polymorphism ratio counts, from the frequencies of each population
    poly_expected = {"p0": (1112, 1173, 1200), "p1": (1101, 1177, 1200), "p2": (1093, 1184, 1200)}
    for pop, (num_poly, num_variable, num_with_data) in poly_expected.items():
        mafs = [max(float(row["REF_FREQ"]), float(row["ALT_FREQS"])) for row in read_table(HERE / f"panel.{pop}.afreq") if int(row["OBS_CT"]) >= 40]
        got = (sum(maf < 0.95 for maf in mafs), sum(maf < 1 for maf in mafs), len(mafs))
        assert got == (num_poly, num_variable, num_with_data), (pop, got)
    smiss = read_table(HERE / "panel.smiss")
    scount = read_table(HERE / "panel.scount")
    assert (smiss[0]["MISSING_CT"], smiss[0]["OBS_CT"], smiss[0]["F_MISS"], scount[0]["HET_CT"]) == ("34", "1200", "0.0283333", "426"), (smiss[0], scount[0])
    assert (smiss[1]["MISSING_CT"], smiss[1]["F_MISS"], scount[1]["HET_CT"]) == ("44", "0.0366667", "397"), (smiss[1], scount[1])
    het = read_table(HERE / "panel.het")
    assert het[0]["OBS_CT"] == "1166", het[0]
    many = read_table(HERE / "many.counts.tsv")
    assert [(row["pos"], row["AN_popA"], row["AC_popA"], row["AN_popB"], row["AC_popB"], row["AN"], row["AC"], row["NHET"], row["NCALLED"]) for row in many[:3]] == [
        ("1000", "36", "32", "57", "52", "93", "84", "9", "46"),
        ("1037", "39", "12", "56", "21", "95", "33", "24", "47"),
        ("1074", "40", "21,7", "55", "24,12", "95", "45,19", "32", "47"),
    ], many[:3]
    many_smiss = read_table(HERE / "many.smiss")
    many_scount = read_table(HERE / "many.scount")
    assert (many_smiss[0]["MISSING_CT"], many_smiss[0]["OBS_CT"], many_scount[0]["HET_CT"]) == ("29", "500", "201"), (many_smiss[0], many_scount[0])
    assert (many_smiss[1]["MISSING_CT"], many_scount[1]["HET_CT"]) == ("25", "195"), (many_smiss[1], many_scount[1])


if __name__ == "__main__":
    check_versions()
    ref_dir = pynei_reference_dir()
    individuals = write_panel(ref_dir)
    assert len(individuals) == 200, len(individuals)
    run_plink2_on_panel()
    write_many_pops()
    run_bcftools_on_many()
    run_plink2_on_many()
    check()
    print("done", file=sys.stderr)
