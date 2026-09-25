r"""It writes the four datasets of the Kosman distances and what R gives for them.

Run from the root of the repository, with the popnei module built, `uv run
maturin develop`, and with R 4.6.1 at /opt/homebrew/bin/Rscript:

    uv run python tests/reference/dists/make_reference.py

It needs three things outside the repository, and it says which one is missing
and what puts it there when it cannot find it:

- R 4.6.1 at /opt/homebrew/bin/Rscript, the R the numbers of the spec were
  taken with, and adegenet 2.1.11 in its libraries, which holds `df2genind`,
  the function that builds the object `gd.kosman` reads. The script refuses
  any other version of either, and a missing adegenet is installed with
  `Rscript -e 'install.packages("adegenet")'`.
- The source package of PopGenReport 3.1.3 from CRAN, under
  ~/.cache/popnei/reference/. PopGenReport does not install on macOS, because
  its dependency terra needs the GDAL library, so the one function this script
  runs, `gd.kosman`, is read from `R/gd.kosman.r` of the unpacked source with
  `source()`; it uses adegenet and base R alone. The script downloads the
  tarball there itself when it is not there, with R's `download.packages`, and
  unpacks the two files it reads from it, so the only run that needs the
  network is the first. CRAN keeps one version of a package under
  src/contrib, so when CRAN moves past 3.1.3 the download gives another
  version, the script refuses it, and 3.1.3 is then taken by hand from
  https://cran.r-project.org/src/contrib/Archive/PopGenReport/ into that
  directory.
- /Users/jose/devel/pynei/test/gwas_reference/sim_missing.vars, the vars file
  of pyNei that holds the panel below. popnei does not read a vars file of
  pyNei, so the script reads it with pyNei, which is a development dependency
  of popnei, and writes it as a VCF.

The four datasets are the ones of "How it is verified" of docs/specs/dists.md:

- panel, 1200 biallelic variants of 200 diploid individuals, s000 to s199,
  3 in 100 genotypes missing, read from the vars file above with pyNei.
- four_alleles, 300 variants of 40 diploid individuals, i00 to i39, four
  alleles, from `numpy.random.default_rng(3)`.
- tetraploid, 200 variants of 12 tetraploid individuals, t00 to t11, three
  alleles, from `default_rng(5)`.
- haploid, the same of 12 haploid individuals, h00 to h11, from
  `default_rng(6)`.

The three random ones are drawn exactly as `ref_export.py` and
`poly_export.py` of docs/reports/kosman-method/ draw them, which is where the
numbers of the spec come from: the alleles with `integers`, then the
genotypes where `random((num_vars, num_individuals)) < 0.05` set to missing
whole. Those scripts stay in the report as they are; this one replaces them
for the tests.

It writes, beside itself, for each dataset:

- `<name>.vcf.gz`, the genotypes as a gzipped VCF that popnei's `open_vcf`
  reads with the ploidy of the dataset. A missing genotype is missing whole,
  `./.` at ploidy 2, and no genotype is half called: no program outside
  popnei checks a half called genotype, and none reached R. Every VCF is
  written with the timestamp of gzip fixed to 0, so that running the script
  again gives the same bytes.
- `<name>.gdkosman.tsv`, what `gd.kosman` gives: the header line
  `dist	n	k_sum` and then one line for every pair of individuals, in the
  order (0, 1), (0, 2), ..., (0, N-1), (1, 2), ..., the order of
  `dist_vector` of docs/specs/dists.md. The three fields of a line are the
  Kosman distance of the pair, n, how many variants `gd.kosman` used for it,
  and k times the sum of d, the integer the tests of the core assert. A cargo
  test splits a line on the tab; a pytest test reads the file with
  `numpy.loadtxt(path, skiprows=1)`.

and, for the two diploid datasets, `<name>.pynei.tsv`, the distance vector of
pyNei's `calc_pairwise_kosman_dists` in the same order, with the header line
`dist`. pyNei takes diploids alone, so the tetraploid and the haploid
datasets have R and nothing else.

A distance is written as the shortest text that reads back as the float64 it
came from, `0.375` and `0.6222222222222222`, which is what `repr` of Python
gives; every one of them round trips. `k_sum` is the distance times the
ploidy times n rounded to the nearest whole number: that product lands just
below or just above its integer in float64, 112.99999999999999 for the pair
h00, h03 of the haploid dataset whose integer is 113, so a test that read the
distance and n and truncated the product would assert the wrong number. The
rounding is refused when the product is further than 1e-10 from a whole
number.

A pair with no variant called in both would have `nan` as its distance, 0 as
its n and 0 as its `k_sum`. No pair of these four datasets is one: the
smallest n is 1096 of the panel, 254 of the four allele dataset and 169 of
each of the other two.

Everything is made and checked in a directory of its own, and what is beside
this script is written only when every one of these checks passed:

- The fourteen literals of the table of "How it is verified" of the spec: for
  each of them that n is the integer of the table, that the distance is within
  1e-9 of the table's, and that the distance times the ploidy times n is
  within 1e-10 of the table's integer. It stops at the first that differs.
- That the distance times the ploidy times n is a whole number within 1e-10
  for every pair of every dataset, and it names the pair furthest from one.
- That a loop over the genotypes which pairs the alleles as "What it gives"
  of the spec says, `kosman` of `poly_export.py` of
  docs/reports/kosman-method/, gives the same n and the same k times the sum
  of d as `gd.kosman` for every pair of every dataset, exactly, and a
  distance within 1e-15 of R's. It reads nothing of R, so it is what says
  that `gd.kosman` computes the distance of the spec over the whole of each
  file and not only over the pairs of the table.
- That the vector of pyNei and the vector of `gd.kosman` agree within 1e-15
  on the two diploid datasets, with the pairs that have no distance in the
  same places, and it prints the largest difference of each.
- That popnei's `open_vcf` reads each VCF it wrote with the number of
  variants, the number of individuals, the names and the ploidy of the
  dataset.
"""

import gzip
import itertools
import math
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

import numpy

HERE = Path(__file__).parent

RSCRIPT = Path("/opt/homebrew/bin/Rscript")
R_VERSION = "4.6.1"
ADEGENET_VERSION = "2.1.11"

POPGENREPORT_VERSION = "3.1.3"
# Outside the repository, so that nothing of CRAN is committed and the
# download happens once for every worktree of popnei on this machine.
POPGENREPORT_CACHE = Path.home() / ".cache" / "popnei" / "reference"
POPGENREPORT_TARBALL = f"PopGenReport_{POPGENREPORT_VERSION}.tar.gz"

PANEL_VARS = Path("/Users/jose/devel/pynei/test/gwas_reference/sim_missing.vars")

MISSING_ALLELE = -1
# The four nucleotides, so that a variant of up to four alleles declares each
# of them in REF and ALT. Which letter an allele gets changes no distance: the
# Kosman distance counts the alleles two genotypes hold in common.
NUCLEOTIDES = "ACGT"

# The literals of the table of "How it is verified" of docs/specs/dists.md:
# the dataset, the two individuals of the pair, k times the sum of d, n, and
# the distance.
LITERALS = [
    ("panel", "s000", "s001", 372, 1122, 0.1657754010695187),
    ("panel", "s000", "s002", 376, 1128, 0.16666666666666666),
    ("panel", "s198", "s199", 351, 1134, 0.15476190476190477),
    ("panel", "s010", "s033", 804, 1123, 0.3579697239536955),
    ("panel", "s116", "s119", 310, 1133, 0.13680494263018536),
    ("four_alleles", "i00", "i01", 325, 269, 0.6040892193),
    ("four_alleles", "i00", "i02", 304, 264, 0.5757575758),
    ("four_alleles", "i00", "i03", 347, 272, 0.6378676471),
    ("tetraploid", "t00", "t01", 282, 188, 0.375),
    ("tetraploid", "t00", "t02", 284, 183, 0.38797814207650272),
    ("tetraploid", "t00", "t03", 294, 183, 0.40163934426229508),
    ("haploid", "h00", "h01", 112, 180, 0.62222222222222223),
    ("haploid", "h00", "h02", 123, 184, 0.66847826086956519),
    ("haploid", "h00", "h03", 113, 179, 0.63128491620111726),
]

# `gd.kosman` on each dataset, run once for the four. It reads one csv per
# dataset, one row per individual and one column per variant, with the alleles
# of a genotype joined by `/` and a missing genotype written `NA`, which is
# the form `df2genind` of adegenet takes. `lower.tri` walks a matrix down each
# column in turn, so the values it takes out of the lower triangle are the
# pairs (0, 1), (0, 2), ..., (1, 2), ..., the order of `dist_vector`.
R_PROGRAM = r"""
arguments <- commandArgs(trailingOnly = TRUE)
work_dir <- arguments[1]
popgenreport_source <- arguments[2]
found_r <- as.character(getRversion())
if (found_r != "R_VERSION") {
  stop(sprintf(paste("the R that ran this is %s and this script needs R_VERSION,",
                     "the R the numbers of docs/specs/dists.md were taken with"),
               found_r))
}
if (!requireNamespace("adegenet", quietly = TRUE)) {
  stop(paste("adegenet is not in R's libraries and this script needs",
             "ADEGENET_VERSION; install.packages(\"adegenet\") puts there the",
             "version CRAN gives today"))
}
found <- as.character(packageVersion("adegenet"))
if (found != "ADEGENET_VERSION") {
  stop(sprintf(paste("the adegenet of R is %s and this script needs",
                     "ADEGENET_VERSION, which is at",
                     "https://cran.r-project.org/src/contrib/Archive/adegenet/"),
               found))
}
suppressMessages(library(adegenet))
source(file.path(popgenreport_source, "R", "gd.kosman.r"))
datasets <- read.table(file.path(work_dir, "datasets.tsv"), header = TRUE,
                       colClasses = c("character", "integer"))
for (row in seq_len(nrow(datasets))) {
  name <- datasets$name[row]
  ploidy <- datasets$ploidy[row]
  genotypes <- read.csv(file.path(work_dir, paste0(name, ".gts.csv")),
                        row.names = 1, colClasses = "character", na.strings = "NA")
  individuals <- df2genind(genotypes, sep = "/", ploidy = ploidy,
                           NA.char = "NA", type = "codom")
  kosman <- gd.kosman(individuals)
  distances <- kosman$geneticdist[lower.tri(kosman$geneticdist)]
  num_vars <- kosman$loci_used[lower.tri(kosman$loci_used)]
  writeLines(c("dist\tn", sprintf("%.17g\t%d", distances, as.integer(num_vars))),
             file.path(work_dir, paste0(name, ".gdkosman.tsv")))
  writeLines(indNames(individuals), file.path(work_dir, paste0(name, ".indnames.txt")))
}
""".replace("ADEGENET_VERSION", ADEGENET_VERSION).replace("R_VERSION", R_VERSION)


class Dataset:
    """One of the four datasets: its genotypes, its individuals and its VCF.

    `gts` is an int8 array of num_vars x num_individuals x ploidy with the
    allele numbers of every genotype, -1 for a missing allele. `chroms`,
    `positions` and `ids` are the columns of the VCF, one value per variant.
    """

    def __init__(self, name, gts, individuals, chroms, positions, ids):
        self.name = name
        self.gts = gts
        self.individuals = individuals
        self.chroms = chroms
        self.positions = positions
        self.ids = ids

    @property
    def num_vars(self):
        return self.gts.shape[0]

    @property
    def ploidy(self):
        return self.gts.shape[2]

    @property
    def num_alleles(self):
        """How many alleles the VCF declares in every variant, REF and ALT."""
        return int(self.gts.max()) + 1

    @property
    def pairs(self):
        """The pairs of individuals in the order of `dist_vector`."""
        return list(itertools.combinations(range(len(self.individuals)), 2))


def panel_dataset():
    """The panel, read from the vars file of pyNei with pyNei."""
    import pynei

    if not PANEL_VARS.is_file():
        raise SystemExit(
            f"{PANEL_VARS} is not there, and this script reads the panel of "
            f'"How it is verified" of docs/specs/dists.md from it. It is the '
            f"file test/gwas_reference/sim_missing.vars of pyNei's repository, "
            f"https://github.com/JoseBlanca/pynei"
        )
    variants = pynei.load_vars(str(PANEL_VARS))
    chunks = list(variants.iter_vars_chunks())
    gts = numpy.concatenate([chunk.gts.gt_values for chunk in chunks])
    info = [chunk.vars_info for chunk in chunks]
    chroms = [str(value) for frame in info for value in frame["chrom"]]
    positions = [int(value) for frame in info for value in frame["pos"]]
    ids = [str(value) for frame in info for value in frame["id"]]
    return Dataset(
        "panel", gts, [str(name) for name in variants.samples], chroms, positions, ids
    )


def random_dataset(name, prefix, seed, num_vars, num_individuals, ploidy, num_alleles):
    """A dataset drawn as the spec's "How it is verified" draws the three.

    The alleles come out of `integers` and then the genotypes where
    `random((num_vars, num_individuals))` is below 0.05 are set to missing
    whole. The two draws are in this order and nothing is drawn between them,
    so the genotypes are the ones the numbers of the spec were taken from.
    """
    rng = numpy.random.default_rng(seed)
    gts = rng.integers(0, num_alleles, size=(num_vars, num_individuals, ploidy))
    gts = gts.astype(numpy.int8)
    gts[rng.random((num_vars, num_individuals)) < 0.05] = MISSING_ALLELE
    individuals = [f"{prefix}{index:02d}" for index in range(num_individuals)]
    chroms = ["chr1"] * num_vars
    positions = [100 * (index + 1) for index in range(num_vars)]
    ids = [f"var{index:04d}" for index in range(num_vars)]
    return Dataset(name, gts, individuals, chroms, positions, ids)


def datasets_of_the_spec():
    """The four datasets of "How it is verified" of docs/specs/dists.md."""
    return [
        panel_dataset(),
        random_dataset("four_alleles", "i", 3, 300, 40, 2, 4),
        random_dataset("tetraploid", "t", 5, 200, 12, 4, 3),
        random_dataset("haploid", "h", 6, 200, 12, 1, 3),
    ]


def check_rscript():
    """It stops unless the Rscript of the path above is there."""
    if not RSCRIPT.is_file():
        raise SystemExit(
            f"there is no Rscript at {RSCRIPT}, and this script runs "
            f"`gd.kosman` of PopGenReport {POPGENREPORT_VERSION} under R"
        )


def popgenreport_source():
    """The unpacked source of PopGenReport, downloaded when it is not there.

    It gives the directory that holds `R/gd.kosman.r`, and it stops when what
    CRAN gives is not the version the spec was written against.
    """
    unpacked = POPGENREPORT_CACHE / "PopGenReport"
    description = unpacked / "DESCRIPTION"
    kosman = unpacked / "R" / "gd.kosman.r"
    if not (description.is_file() and kosman.is_file()):
        POPGENREPORT_CACHE.mkdir(parents=True, exist_ok=True)
        tarball = POPGENREPORT_CACHE / POPGENREPORT_TARBALL
        if not tarball.is_file():
            download = (
                f'download.packages("PopGenReport", destdir="{POPGENREPORT_CACHE}",'
                f' type="source", repos="https://cloud.r-project.org")'
            )
            subprocess.run([str(RSCRIPT), "-e", download], check=True)
        if not tarball.is_file():
            raise SystemExit(
                f"CRAN did not give {POPGENREPORT_TARBALL}, and this script "
                f"needs PopGenReport {POPGENREPORT_VERSION}. CRAN keeps one "
                f"version of a package under src/contrib, so take that one "
                f"from https://cran.r-project.org/src/contrib/Archive/"
                f"PopGenReport/ into {POPGENREPORT_CACHE}"
            )
        shutil.rmtree(unpacked, ignore_errors=True)
        with tarfile.open(tarball) as archive:
            for member in ("PopGenReport/DESCRIPTION", "PopGenReport/R/gd.kosman.r"):
                archive.extract(member, POPGENREPORT_CACHE, filter="data")
    version = ""
    for line in description.read_text().splitlines():
        if line.startswith("Version:"):
            version = line.removeprefix("Version:").strip()
    if version != POPGENREPORT_VERSION:
        raise SystemExit(
            f"the PopGenReport source under {POPGENREPORT_CACHE} is {version}, "
            f"and this script needs {POPGENREPORT_VERSION}. Take that version "
            f"from https://cran.r-project.org/src/contrib/Archive/PopGenReport/ "
            f"into that directory and unpack it there"
        )
    return unpacked


def genotype_text(alleles):
    """One genotype of a VCF: its alleles joined, or missing whole."""
    if (alleles < 0).any():
        return "/".join("." for _ in alleles)
    return "/".join(str(int(allele)) for allele in alleles)


def vcf_text(dataset):
    """The whole VCF of a dataset, header and data lines."""
    chroms = sorted(set(dataset.chroms))
    alleles = NUCLEOTIDES[: dataset.num_alleles]
    lines = [
        "##fileformat=VCFv4.4\n",
        *[f"##contig=<ID={chrom}>\n" for chrom in chroms],
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            + list(dataset.individuals)
        )
        + "\n",
    ]
    for index in range(dataset.num_vars):
        fields = [
            dataset.chroms[index],
            str(dataset.positions[index]),
            dataset.ids[index],
            alleles[0],
            ",".join(alleles[1:]) if len(alleles) > 1 else ".",
            ".",
            "PASS",
            ".",
            "GT",
        ]
        fields += [genotype_text(genotype) for genotype in dataset.gts[index]]
        lines.append("\t".join(fields) + "\n")
    return "".join(lines)


def write_vcf(dataset, directory):
    """It writes `<name>.vcf.gz` and gives its path.

    gzip writes the time of the run into its header, which would make the
    bytes of the file different at every run, so the time is fixed to 0 and
    the name of the file is left out of the header.
    """
    path = directory / f"{dataset.name}.vcf.gz"
    with (
        path.open("wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed,
    ):
        compressed.write(vcf_text(dataset).encode())
    return path


def write_genind_csv(dataset, work_dir):
    """The csv of a dataset that `df2genind` of adegenet reads.

    One row per individual, one column per variant, the alleles of a genotype
    joined by `/` and a missing genotype written `NA`.
    """
    header = "," + ",".join(f"v{index:04d}" for index in range(dataset.num_vars))
    lines = [header]
    for individual, name in enumerate(dataset.individuals):
        genotypes = []
        for index in range(dataset.num_vars):
            alleles = dataset.gts[index, individual]
            if (alleles < 0).any():
                genotypes.append("NA")
            else:
                genotypes.append("/".join(str(int(allele)) for allele in alleles))
        lines.append(name + "," + ",".join(genotypes))
    (work_dir / f"{dataset.name}.gts.csv").write_text("\n".join(lines) + "\n")


def run_gd_kosman(datasets, work_dir):
    """It runs `gd.kosman` on the four datasets and reads what R wrote.

    It gives, for each dataset by name, the list of (distance, n) of every
    pair in the order of `dist_vector`.
    """
    source = popgenreport_source()
    for dataset in datasets:
        write_genind_csv(dataset, work_dir)
    manifest = ["name\tploidy"]
    manifest += [f"{dataset.name}\t{dataset.ploidy}" for dataset in datasets]
    (work_dir / "datasets.tsv").write_text("\n".join(manifest) + "\n")
    program = work_dir / "gd_kosman.R"
    program.write_text(R_PROGRAM)
    run = subprocess.run(
        [str(RSCRIPT), str(program), str(work_dir), str(source)],
        capture_output=True,
        text=True,
        check=False,
    )
    if run.returncode != 0:
        raise SystemExit(
            f"R stopped before it gave the numbers of `gd.kosman`. It runs "
            f"under R {R_VERSION} with adegenet {ADEGENET_VERSION}, and it "
            f"refuses any other version of either. R said:\n"
            f"{run.stdout}{run.stderr}"
        )
    from_r = {}
    for dataset in datasets:
        names = (work_dir / f"{dataset.name}.indnames.txt").read_text().split()
        if names != list(dataset.individuals):
            raise SystemExit(
                f"R read the individuals of {dataset.name} in another order "
                f"than the VCF has them: {names[:3]} against "
                f"{list(dataset.individuals)[:3]}"
            )
        lines = (work_dir / f"{dataset.name}.gdkosman.tsv").read_text().splitlines()
        values = []
        for line in lines[1:]:
            distance, num_vars = line.split("\t")
            values.append((float(distance), int(num_vars)))
        if len(values) != len(dataset.pairs):
            raise SystemExit(
                f"R gave {len(values)} pairs of {dataset.name} and the dataset "
                f"has {len(dataset.pairs)}"
            )
        from_r[dataset.name] = values
    return from_r


def pynei_dists(dataset):
    """pyNei's distance vector of a diploid dataset, one value per pair."""
    import pynei

    if dataset.name == "panel":
        variants = pynei.load_vars(str(PANEL_VARS))
    else:
        variants = pynei.Variants.from_gt_array(
            dataset.gts, samples=list(dataset.individuals)
        )
    distances = pynei.calc_pairwise_kosman_dists(variants).dist_vector
    # Floats of Python, so that `repr` gives the digits of the value and not
    # the `np.float64(...)` that numpy's own repr writes around them.
    return [float(distance) for distance in distances]


def kosman_loop(dataset):
    """k times the sum of d and n of every pair, from the genotypes alone.

    It is the loop `kosman` of `poly_export.py` of
    docs/reports/kosman-method/, the one the spec's "How it is verified"
    checked `gd.kosman` against pair by pair, with its loop over the variants
    made by numpy. At a variant where both genotypes are called, d is 1 minus
    the copies the two hold in common over the ploidy k, the copies in common
    being the sum over the alleles a of the smaller of the copies of a each
    genotype holds; so k times the sum of d over the variants of a pair is k
    times n minus those copies added over the same variants. A missing
    genotype holds no copy of any allele, so a variant where either genotype
    is missing adds nothing to that sum and does not have to be taken out of
    it.

    It gives, for every pair in the order of `dist_vector`, k times the sum of
    d and n, both integers.
    """
    # copies[allele, individual, variant] and called[individual, variant],
    # laid this way so that the two rows of a pair are read side by side.
    copies = numpy.stack(
        [
            (dataset.gts == allele).sum(axis=2).T.astype(numpy.int32)
            for allele in range(dataset.num_alleles)
        ]
    )
    called = (dataset.gts >= 0).all(axis=2).T
    values = []
    for first, second in dataset.pairs:
        num_vars = int(numpy.count_nonzero(called[first] & called[second]))
        in_common = 0
        for allele in range(dataset.num_alleles):
            in_common += int(
                numpy.minimum(copies[allele, first], copies[allele, second]).sum()
            )
        values.append((dataset.ploidy * num_vars - in_common, num_vars))
    return values


def pair_name(dataset, index):
    """The two individuals of the pair at `index` of `dist_vector`."""
    first, second = dataset.pairs[index]
    return f"{dataset.individuals[first]}, {dataset.individuals[second]}"


def worst(largest, where):
    """The largest difference that was found, with the pair it was on.

    A largest difference of 0 was on every pair and on none in particular, so
    naming one of them would say that the others were smaller.
    """
    if largest == 0.0:
        return "0.0, which is what every pair gives"
    return f"{largest}, on the pair {where}"


def pair_index(dataset, name_a, name_b):
    """Where a pair of individuals is in the order of `dist_vector`."""
    individuals = list(dataset.individuals)
    for name in (name_a, name_b):
        if name not in individuals:
            raise SystemExit(
                f"{name} is not an individual of {dataset.name}, whose "
                f"{len(individuals)} individuals are {individuals[0]} to "
                f"{individuals[-1]}"
            )
    first = individuals.index(name_a)
    second = individuals.index(name_b)
    return dataset.pairs.index((first, second))


def check_literals(datasets, from_r):
    """It compares what R gave with the fourteen literals of the spec.

    It stops at the first row of the table that differs. The integers are
    compared exactly, the distances within 1e-9, the digits of the shortest
    ones of the table, and k times the sum of d, which the table has as a
    whole number, within 1e-10 of it.
    """
    by_name = {dataset.name: dataset for dataset in datasets}
    for name, name_a, name_b, k_times_sum, num_vars, distance in LITERALS:
        dataset = by_name[name]
        index = pair_index(dataset, name_a, name_b)
        from_r_distance, from_r_num_vars = from_r[name][index]
        where = f"{name}, {name_a}, {name_b}"
        if from_r_num_vars != num_vars:
            raise SystemExit(
                f"{where}: gd.kosman used {from_r_num_vars} variants and the "
                f"spec's table says {num_vars}"
            )
        if abs(from_r_distance - distance) > 1e-9:
            raise SystemExit(
                f"{where}: gd.kosman gives the distance {from_r_distance!r} "
                f"and the spec's table says {distance!r}"
            )
        from_r_sum = from_r_distance * dataset.ploidy * from_r_num_vars
        if abs(from_r_sum - k_times_sum) > 1e-10:
            raise SystemExit(
                f"{where}: the distance of gd.kosman times the ploidy times n "
                f"is {from_r_sum!r} and the spec's table says {k_times_sum}"
            )
    print(f"the {len(LITERALS)} literals of the spec's table are what R gives")


def whole_sums(dataset, from_r):
    """k times the sum of d of every pair, as the whole number it is.

    R gives the distance and n, and k times the sum of d is the distance times
    the ploidy times n. That product lands just below or just above its whole
    number in float64, 112.99999999999999 for the pair h00, h03 whose integer
    is 113, so it is rounded, and the rounding is refused when the product is
    further than 1e-10 from a whole number, because the integer would then be
    a guess. A pair with no variant called in both has no distance and a sum
    of 0.

    It prints the pair whose product is furthest from a whole number and the
    smallest n of the dataset.
    """
    furthest = 0.0
    furthest_pair = pair_name(dataset, 0)
    sums = []
    for index, (distance, num_vars) in enumerate(from_r):
        if num_vars == 0:
            sums.append(0)
            continue
        if math.isnan(distance):
            raise SystemExit(
                f"{dataset.name}: gd.kosman gives no distance to the pair "
                f"{pair_name(dataset, index)}, which has n = {num_vars}"
            )
        product = distance * dataset.ploidy * num_vars
        whole = round(product)
        away = abs(product - whole)
        if away > furthest:
            furthest = away
            furthest_pair = pair_name(dataset, index)
        sums.append(whole)
    if not furthest <= 1e-10:
        raise SystemExit(
            f"{dataset.name}: k times the sum of d of the pair {furthest_pair} "
            f"is {furthest} away from a whole number, and the tests of the "
            f"core assert it as an integer"
        )
    smallest = min(num_vars for _, num_vars in from_r)
    print(
        f"{dataset.name}: over the {len(from_r)} pairs, whose smallest n is "
        f"{smallest}, k times the sum of d is away from a whole number by at "
        f"most {worst(furthest, furthest_pair)}"
    )
    return sums


def check_against_the_loop(dataset, from_r, sums):
    """It compares R with the loop over the genotypes, pair by pair.

    The loop of `kosman_loop` reads the genotypes and nothing of R, so it is
    what says that `gd.kosman` computes the d of "What it gives" of
    docs/specs/dists.md on every pair and not only on the three of the table.
    The two integers are compared exactly and the distances within 1e-15, the
    tolerance of the spec, a pair with no distance on either side being one
    with no distance on both.
    """
    from_loop = kosman_loop(dataset)
    largest = 0.0
    largest_pair = pair_name(dataset, 0)
    for index, (from_r_pair, from_loop_pair) in enumerate(
        zip(from_r, from_loop, strict=True)
    ):
        from_r_distance, from_r_num_vars = from_r_pair
        loop_sum, loop_num_vars = from_loop_pair
        where = pair_name(dataset, index)
        if loop_num_vars != from_r_num_vars:
            raise SystemExit(
                f"{dataset.name}: the pair {where} has {loop_num_vars} "
                f"variants called in both and gd.kosman used {from_r_num_vars}"
            )
        if loop_sum != sums[index]:
            raise SystemExit(
                f"{dataset.name}: k times the sum of d of the pair {where} is "
                f"{loop_sum} over the genotypes and {sums[index]} from what "
                f"gd.kosman gives"
            )
        if loop_num_vars == 0:
            if not math.isnan(from_r_distance):
                raise SystemExit(
                    f"{dataset.name}: the pair {where} has no variant called "
                    f"in both and gd.kosman gives it the distance "
                    f"{from_r_distance!r}"
                )
            continue
        if math.isnan(from_r_distance):
            raise SystemExit(
                f"{dataset.name}: gd.kosman gives no distance to the pair "
                f"{where}, which has {loop_num_vars} variants called in both"
            )
        loop_distance = loop_sum / (dataset.ploidy * loop_num_vars)
        away = abs(loop_distance - from_r_distance)
        if away > largest:
            largest = away
            largest_pair = where
    if not largest <= 1e-15:
        raise SystemExit(
            f"{dataset.name}: the loop over the genotypes and gd.kosman are "
            f"{largest} apart on the pair {largest_pair}, and the spec allows "
            f"1e-15"
        )
    print(
        f"{dataset.name}: the loop over the genotypes gives the two integers "
        f"of every one of the {len(from_loop)} pairs, and its distance differs "
        f"from gd.kosman's by at most {worst(largest, largest_pair)}"
    )


def check_against_pynei(dataset, from_r):
    """It compares pyNei's vector with R's and prints the largest difference.

    Both divide the same two integers once, so they agree to the last place of
    a float64 of that size; the tolerance is the 1e-15 of the spec. A pair
    with no distance is NaN on both sides or the comparison stops: NaN is not
    above any tolerance, so a pair that one side leaves out and the other
    does not would otherwise pass.
    """
    from_pynei = pynei_dists(dataset)
    from_r_distances = [distance for distance, _ in from_r]
    largest = 0.0
    largest_pair = pair_name(dataset, 0)
    for index, (left, right) in enumerate(
        zip(from_pynei, from_r_distances, strict=True)
    ):
        where = pair_name(dataset, index)
        if math.isnan(left) != math.isnan(right):
            has_one = "gd.kosman" if math.isnan(left) else "pyNei"
            raise SystemExit(
                f"{dataset.name}: the pair {where} has a distance in "
                f"{has_one} and none in the other"
            )
        if math.isnan(left):
            continue
        away = abs(left - right)
        if away > largest:
            largest = away
            largest_pair = where
    if not largest <= 1e-15:
        raise SystemExit(
            f"{dataset.name}: pyNei and gd.kosman are {largest} apart on the "
            f"pair {largest_pair}, and the spec allows 1e-15"
        )
    print(
        f"{dataset.name}: over the {len(from_pynei)} pairs, pyNei and "
        f"gd.kosman differ by at most {worst(largest, largest_pair)}"
    )
    return from_pynei


def check_with_popnei(dataset, vcf_path):
    """It reads the VCF with `open_vcf` and checks that nothing was lost.

    The genotypes popnei reads back have to be the array that was given to R,
    element by element: the tests read the distances of R and the genotypes of
    the VCF, and a VCF that held other genotypes would make every one of them
    wrong.
    """
    try:
        import popnei
    except ImportError:
        raise SystemExit(
            "popnei is not in the environment, and this script reads the VCFs "
            "it writes with popnei's open_vcf. Build the module with "
            "`uv run maturin develop`"
        ) from None
    variants = popnei.open_vcf(vcf_path, ploidy=dataset.ploidy)
    if list(variants.individuals) != list(dataset.individuals):
        raise SystemExit(
            f"{vcf_path.name}: popnei reads the individuals "
            f"{list(variants.individuals)[:3]} and the dataset has "
            f"{list(dataset.individuals)[:3]}"
        )
    if variants.ploidy != dataset.ploidy:
        raise SystemExit(
            f"{vcf_path.name}: popnei reads a ploidy of {variants.ploidy} and "
            f"the dataset has {dataset.ploidy}"
        )
    blocks = variants.iter_blocks(fields=())
    gts = numpy.concatenate([block.gts for block in blocks])
    num_vars = gts.shape[0]
    if num_vars != dataset.num_vars:
        raise SystemExit(
            f"{vcf_path.name}: popnei reads {num_vars} variants and the "
            f"dataset has {dataset.num_vars}"
        )
    if not numpy.array_equal(gts, dataset.gts):
        differ = numpy.argwhere(gts != dataset.gts)[0]
        raise SystemExit(
            f"{vcf_path.name}: popnei reads the genotype of the variant "
            f"{differ[0]} of the individual {differ[1]} as "
            f"{gts[differ[0], differ[1]]} and the dataset has "
            f"{dataset.gts[differ[0], differ[1]]}"
        )
    print(
        f"{vcf_path.name}: popnei reads the genotypes R was given, "
        f"{num_vars} variants of {variants.num_individuals} individuals at "
        f"ploidy {variants.ploidy}, {vcf_path.stat().st_size} bytes"
    )


def write_file(path, header, lines):
    """A header line and then one line for each pair, in the order of
    `dist_vector`."""
    path.write_text("\n".join([header, *lines]) + "\n")


def main():
    """Everything is made and checked in a directory of its own, and what is
    beside this script is written only when nothing differed."""
    check_rscript()
    datasets = datasets_of_the_spec()
    with tempfile.TemporaryDirectory() as directory:
        work_dir = Path(directory)
        from_r = run_gd_kosman(datasets, work_dir)
        check_literals(datasets, from_r)
        sums = {}
        from_pynei = {}
        for dataset in datasets:
            sums[dataset.name] = whole_sums(dataset, from_r[dataset.name])
            check_against_the_loop(dataset, from_r[dataset.name], sums[dataset.name])
            if dataset.ploidy == 2:
                from_pynei[dataset.name] = check_against_pynei(
                    dataset, from_r[dataset.name]
                )
        for dataset in datasets:
            check_with_popnei(dataset, write_vcf(dataset, work_dir))
        for dataset in datasets:
            name = dataset.name
            shutil.copyfile(work_dir / f"{name}.vcf.gz", HERE / f"{name}.vcf.gz")
            write_file(
                HERE / f"{name}.gdkosman.tsv",
                "dist\tn\tk_sum",
                [
                    f"{distance!r}\t{num_vars}\t{k_sum}"
                    for (distance, num_vars), k_sum in zip(
                        from_r[name], sums[name], strict=True
                    )
                ],
            )
            if name in from_pynei:
                write_file(
                    HERE / f"{name}.pynei.tsv",
                    "dist",
                    [repr(distance) for distance in from_pynei[name]],
                )


if __name__ == "__main__":
    main()
