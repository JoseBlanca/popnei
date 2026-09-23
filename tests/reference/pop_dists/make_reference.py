r"""It writes the multiallelic panel of the population distances and what the
reference programs give for both panels.

Run from the root of the repository:

    uv run python tests/reference/pop_dists/make_reference.py

It needs two programs outside the repository, and it says which one is missing
when it cannot find it:

- plink2 v2.0.0-a.7.7 on the PATH, which gives Hudson's F_ST with
  `--fst ... method=hudson`. The script refuses another version, because the
  numbers of "How it is verified" of docs/specs/dists.md were taken with this
  one.
- R 4.6.1 at /opt/homebrew/bin/Rscript with adegenet 2.1.11, which gives the
  chord distance as `dist.genpop(..., method = 2)`; mmod 1.3.3, which gives
  Jost's D, Nei's G_ST and Hedrick's G'_ST; and admixtools 2.0.10, which gives
  f_2 and its jackknife standard error. A missing one is installed with
  `Rscript -e 'install.packages("adegenet")'` and the like; admixtools comes
  from `remotes::install_github("uqrmaie1/admixtools")`.

The two panels of docs/specs/dists.md:

- panel, biallelic: tests/reference/dists/panel.vcf.gz, 1200 variants of 200
  diploid individuals over two chromosomes, 3 in 100 genotypes missing whole,
  with the three populations of tests/reference/stats/panel_pops.txt. It is
  written by the reference script of the Kosman distances and is read here.
- micro, multiallelic: written here, 120 loci of 90 diploid individuals in
  three populations of 30, six alleles to a locus written as repeat lengths so
  that a record looks like the microsatellite it stands for, 4 in 100 genotypes
  missing whole. Its allele frequencies are drawn per population from a
  Dirichlet around a per locus ancestral one, so that the populations differ by
  drift. The seed is 7 and numpy's PCG64 generator gives the same draw on every
  machine, so the file can be written again; it is kept in git all the same,
  because the literals of the spec are of these genotypes.

It writes, beside itself: micro.vcf.gz and micro_pops.txt, the panel; and, for
each panel, <name>.hudson_fst.tsv from plink2, <name>.chord.tsv from adegenet,
<name>.mmod.tsv from mmod, and, for the biallelic one alone, panel.f2.tsv from
admixtools, which reads biallelic genotypes only.
"""

import gzip
import itertools
import pathlib
import shutil
import subprocess
import sys
import tempfile

import numpy

HERE = pathlib.Path(__file__).parent
REFERENCE = HERE.parent
RSCRIPT = pathlib.Path("/opt/homebrew/bin/Rscript")
PLINK2_VERSION = "v2.0.0-a.7.7"
R_VERSIONS = {"adegenet": "2.1.11", "mmod": "1.3.3", "admixtools": "2.0.10"}

NUM_LOCI, NUM_ALLELES, PER_POP, SEED = 120, 6, 30, 7
MISSING_RATE = 0.04
POP_NAMES = ["p0", "p1", "p2"]


def refuse_a_missing_program():
    """It stops unless plink2 of the right version and Rscript are there."""
    plink2 = shutil.which("plink2")
    if plink2 is None:
        sys.exit("there is no plink2 on the PATH; brew install plink2 puts one there")
    printed = subprocess.run([plink2, "--version"], capture_output=True, text=True).stdout
    if PLINK2_VERSION not in printed:
        sys.exit(
            f"the numbers of docs/specs/dists.md were taken with plink2 "
            f"{PLINK2_VERSION} and this is {printed.strip()}"
        )
    if not RSCRIPT.exists():
        sys.exit(f"there is no Rscript at {RSCRIPT}")
    for package, version in R_VERSIONS.items():
        out = subprocess.run(
            [RSCRIPT, "-e", f'cat(as.character(packageVersion("{package}")))'],
            capture_output=True, text=True,
        )
        if version not in out.stdout:
            sys.exit(
                f"the numbers of docs/specs/dists.md were taken with {package} "
                f"{version}, and R gives {out.stdout.strip() or out.stderr.strip()}"
            )
    return plink2


def write_micro_panel():
    """The multiallelic panel: its genotypes, its VCF and its populations."""
    rng = numpy.random.default_rng(SEED)
    individuals = [f"i{n:03d}" for n in range(len(POP_NAMES) * PER_POP)]
    gts = numpy.full((NUM_LOCI, len(individuals), 2), -1, dtype=numpy.int8)
    for locus in range(NUM_LOCI):
        ancestral = rng.dirichlet(numpy.full(NUM_ALLELES, 1.5))
        for number, _ in enumerate(POP_NAMES):
            freqs = rng.dirichlet(ancestral * 12)
            for individual in range(number * PER_POP, (number + 1) * PER_POP):
                gts[locus, individual] = rng.choice(NUM_ALLELES, size=2, p=freqs)
    gts[rng.random(gts.shape[:2]) < MISSING_RATE] = -1

    # an allele is a number of repeats of AT, which is what a microsatellite
    # varies by, so that the record looks like the marker it stands for
    alleles = [f"AT{'AT' * n}" for n in range(NUM_ALLELES)]
    lines = [
        "##fileformat=VCFv4.4",
        "##contig=<ID=chr1>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(individuals),
    ]
    for locus in range(NUM_LOCI):
        calls = ["/".join("." if allele < 0 else str(allele) for allele in genotype)
                 for genotype in gts[locus]]
        lines.append("\t".join([
            "chr1", str((locus + 1) * 1000), f"ssr{locus:03d}", alleles[0],
            ",".join(alleles[1:]), ".", "PASS", ".", "GT", *calls,
        ]))
    # gzip writes the time of the run into its header, which would make the
    # file differ from one run to the next, so mtime is fixed at 0
    with gzip.GzipFile(HERE / "micro.vcf.gz", "wb", mtime=0) as compressed:
        compressed.write(("\n".join(lines) + "\n").encode())

    pop_of = {individual: POP_NAMES[number // PER_POP]
              for number, individual in enumerate(individuals)}
    (HERE / "micro_pops.txt").write_text(
        "IID\tpopcat\n" + "\n".join(f"{i}\t{pop_of[i]}" for i in individuals) + "\n")
    return individuals, pop_of, gts


def read_vcf(path):
    """The individuals and the genotypes of a VCF, variants x individuals x ploidy."""
    individuals, rows = None, []
    with gzip.open(path, "rt") as vcf:
        for line in vcf:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                individuals = fields[9:]
                continue
            rows.append([field.replace("|", "/") for field in fields[9:]])
    gts = numpy.array(
        [[[-1 if a == "." else int(a) for a in gt.split("/")] for gt in row]
         for row in rows], dtype=numpy.int8)
    return individuals, gts


def write_genind_csv(path, individuals, pop_of, gts):
    """The csv that `df2genind` of adegenet reads: one row for each individual."""
    lines = ["id,pop," + ",".join(f"v{n:04d}" for n in range(gts.shape[0]))]
    for column, individual in enumerate(individuals):
        calls = ["NA" if (gts[variant, column] < 0).any()
                 else "/".join(str(a) for a in gts[variant, column])
                 for variant in range(gts.shape[0])]
        lines.append(f"{individual},{pop_of[individual]}," + ",".join(calls))
    path.write_text("\n".join(lines) + "\n")


R_PROGRAM = r"""
suppressMessages({library(adegenet); library(mmod)})
args <- commandArgs(trailingOnly = TRUE)
tab <- read.csv(args[1], colClasses = "character")
gts <- tab[, -(1:2)]
rownames(gts) <- tab$id
individuals <- df2genind(gts, sep = "/", ploidy = 2, pop = factor(tab$pop),
                         NA.char = "NA")
by_pop <- genind2genpop(individuals, quiet = TRUE)
pairs <- combn(sort(unique(tab$pop)), 2)
labels <- apply(pairs, 2, function(p) paste(p, collapse = "-"))
chord <- as.vector(dist.genpop(by_pop, method = 2))
writeLines(c("pair\tchord", sprintf("%s\t%.17g", labels, chord)), args[2])
writeLines(c("pair\tdest\tgst\tgst_hedrick",
             sprintf("%s\t%.17g\t%.17g\t%.17g", labels,
                     as.vector(pairwise_D(individuals, linearized = FALSE)),
                     as.vector(pairwise_Gst_Nei(individuals, linearized = FALSE)),
                     as.vector(pairwise_Gst_Hedrick(individuals, linearized = FALSE)))),
           args[3])
"""

ADMIX_PROGRAM = r"""
suppressMessages(library(admixtools))
args <- commandArgs(trailingOnly = TRUE)
blocks <- f2_from_geno(args[1], maxmiss = 1, blgsize = 100000,
                       adjust_pseudohaploid = FALSE, verbose = FALSE)
est <- as.data.frame(f2(blocks))
writeLines(c("pair\tf2\tstandard_error",
             sprintf("%s-%s\t%.17g\t%.17g", est$pop1, est$pop2, est$est, est$se)),
           args[2])
"""


def run_plink2_fst(plink2, vcf, pops_file, name):
    """Hudson's F_ST of every pair, and of every variant of the first panel."""
    with tempfile.TemporaryDirectory() as work:
        out = pathlib.Path(work) / "fst"
        subprocess.run(
            [plink2, "--vcf", str(vcf), "--pheno", str(pops_file), "--fst", "popcat",
             "method=hudson", "report-variants", "--out", str(out)],
            check=True, capture_output=True)
        summary = (out.with_suffix(".fst.summary")).read_text()
        (HERE / f"{name}.hudson_fst.tsv").write_text(
            "pair\thudson_fst\n" + "\n".join(
                f"{a}-{b}\t{value}" for a, b, value in
                (line.split("\t") for line in summary.splitlines()[1:])) + "\n")
        for variants_file in sorted(pathlib.Path(work).glob("*.fst.var")):
            shutil.copy(variants_file, HERE / f"{name}.{variants_file.name}")


def run_admixtools_f2(plink2, vcf, pop_of, name):
    """f_2 of every pair with its jackknife standard error, over 100 kb groups."""
    with tempfile.TemporaryDirectory() as work:
        prefix = pathlib.Path(work) / "panel"
        subprocess.run([plink2, "--vcf", str(vcf), "--make-bed", "--out", str(prefix)],
                       check=True, capture_output=True)
        # admixtools takes the population of an individual from the first
        # column of the .fam, which plink2 wrote as the family
        fam = prefix.with_suffix(".fam")
        fam.write_text("\n".join(
            " ".join([pop_of[f[1]], *f[1:]]) for f in
            (line.split() for line in fam.read_text().splitlines())) + "\n")
        program = pathlib.Path(work) / "f2.R"
        program.write_text(ADMIX_PROGRAM)
        subprocess.run([RSCRIPT, str(program), str(prefix), str(HERE / f"{name}.f2.tsv")],
                       check=True, capture_output=True)


def run_r_distances(csv, name):
    """The chord distance from adegenet, and Jost's D, G_ST and G'_ST from mmod."""
    with tempfile.TemporaryDirectory() as work:
        program = pathlib.Path(work) / "distances.R"
        program.write_text(R_PROGRAM)
        subprocess.run(
            [RSCRIPT, str(program), str(csv), str(HERE / f"{name}.chord.tsv"),
             str(HERE / f"{name}.mmod.tsv")], check=True, capture_output=True)


def main():
    plink2 = refuse_a_missing_program()

    micro_individuals, micro_pop_of, micro_gts = write_micro_panel()

    panel_vcf = REFERENCE / "dists" / "panel.vcf.gz"
    panel_pops_file = REFERENCE / "stats" / "panel_pops.txt"
    if not panel_vcf.exists():
        sys.exit(f"{panel_vcf} is not there; "
                 "tests/reference/dists/make_reference.py writes it")
    panel_individuals, panel_gts = read_vcf(panel_vcf)
    panel_pop_of = dict(line.split("\t")
                        for line in panel_pops_file.read_text().splitlines()[1:])

    with tempfile.TemporaryDirectory() as work:
        for name, individuals, pop_of, gts, vcf, pops_file in (
            ("panel", panel_individuals, panel_pop_of, panel_gts, panel_vcf,
             panel_pops_file),
            ("micro", micro_individuals, micro_pop_of, micro_gts,
             HERE / "micro.vcf.gz", HERE / "micro_pops.txt"),
        ):
            csv = pathlib.Path(work) / f"{name}.gts.csv"
            write_genind_csv(csv, individuals, pop_of, gts)
            run_plink2_fst(plink2, vcf, pops_file, name)
            run_r_distances(csv, name)
            print(f"{name}: {gts.shape[0]} variants, {gts.shape[1]} individuals, "
                  f"{len(set(pop_of.values()))} populations")
        # admixtools reads biallelic genotypes only
        run_admixtools_f2(plink2, panel_vcf, panel_pop_of, "panel")


if __name__ == "__main__":
    main()
