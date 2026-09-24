"""It writes what dadi and scikit-allel give for two of the statistics of
docs/specs/diversity.md, the folded site frequency spectrum and F_IS.

Run from the root of the repository, with `uv` in the PATH:

    uv run python tests/reference/diversity/make_reference.py

It reads the panel of docs/specs/stats.md, tests/reference/stats/panel.vcf.gz
with tests/reference/stats/panel_pops_bcftools.txt, both committed, 1200
biallelic diploid variants of 200 individuals in three populations, p0, p1
and p2, of 48, 68 and 84 individuals. It takes the allele counts of each
population at each variant itself, from the VCF, and keeps for a population
the variants at which it called at least `min_num_individuals` genotypes, 20,
which on this panel is all 1200 for each of the three. It writes beside
itself:

    panel_folded_sfs_dadi.tsv   the folded spectrum of each population at a
                                draw of 20 called alleles, 11 bins by 3
                                populations, from dadi
    panel_fis_plain_allel.tsv   F_IS of each population in its plain form,
                                from scikit-allel

The plain form is 1 minus the ratio of the mean observed heterozygosity to
the mean expected heterozygosity with no correction for the sample size.
popnei returns the unbiased form, which uses the mean expected
heterozygosity corrected by c / (c - 1) at each variant, so the value stored
here is not the value popnei gives: "How it is verified" of "The inbreeding
coefficient F_IS" of the spec says how the two are compared and gives both.

Every value is written with all the digits of the float, so that a test can
compare within a tolerance the file does not limit. The spec prints the same
numbers to ten digits.

The two programs run in two Python 3.12 environments that the script makes
with `uv venv --python 3.12` and reuses on a later run, under
`tmp/diversity-reference/`, which .gitignore leaves out of the repository.
Neither program is a development dependency of pyproject.toml: dadi 2.4.4
does not build on the project's Python, 3.14 with the free threading build,
its nlopt dependency failing to compile, and scikit-allel 1.3.13 installs
there but re-enables the global interpreter lock of the process that imports
it, which that Python warns about. No test of popnei imports either of them,
and this script runs each in a process of its own. The same directory holds
the program each process runs, the allele counts dadi projects and the
genotypes scikit-allel reads, all written again on every run.

The script checks what it got against the literals of "How it is verified"
of "The folded site frequency spectrum" and of "The inbreeding coefficient
F_IS" of docs/specs/diversity.md, and stops at the first one that differs.
"""

import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy

HERE = Path(__file__).parent
ROOT = HERE.parent.parent.parent
PANEL_VCF = ROOT / "tests" / "reference" / "stats" / "panel.vcf.gz"
PANEL_POPS = ROOT / "tests" / "reference" / "stats" / "panel_pops_bcftools.txt"
WORK_DIR = ROOT / "tmp" / "diversity-reference"

PLOIDY = 2
MIN_NUM_INDIVIDUALS = 20
NUM_CALLED_ALLELES = 20
REFERENCE_PYTHON = "3.12"
DADI = "dadi==2.4.4"
ALLEL = "scikit-allel==1.3.13"

DADI_PROGRAM = '''"""The folded spectrum of each population, projected by dadi."""

import csv
import sys
from collections import defaultdict

import dadi

counts_path, num_called_alleles = sys.argv[1], int(sys.argv[2])
counts = defaultdict(dict)
with open(counts_path) as handle:
    for row in csv.DictReader(handle, delimiter="\\t"):
        counts[row["pop"]][row["var"]] = (int(row["n0"]), int(row["n1"]))

for pop in sorted(counts):
    data = {
        var: {
            "calls": {pop: pair},
            "segregating": ("A", "T"),
            "outgroup_allele": "A",
        }
        for var, pair in counts[pop].items()
    }
    spectrum = dadi.Spectrum.from_data_dict(
        data, [pop], projections=[num_called_alleles], polarized=False
    )
    for rarer, value in enumerate(spectrum.data[: num_called_alleles // 2 + 1]):
        print(pop, rarer, repr(float(value)))
'''

ALLEL_PROGRAM = '''"""F_IS of each population in its plain form, from scikit-allel."""

import sys

import allel
import numpy

for argument in sys.argv[1:]:
    pop, path = argument.split("=", 1)
    gts = allel.GenotypeArray(numpy.load(path))
    obs_het = allel.heterozygosity_observed(gts)
    frequencies = gts.count_alleles().to_frequencies()
    exp_het = allel.heterozygosity_expected(frequencies, ploidy=gts.shape[2])
    both = ~numpy.isnan(obs_het) & ~numpy.isnan(exp_het)
    print(pop, repr(float(1 - obs_het[both].mean() / exp_het[both].mean())))
'''

# The tables of "How it is verified" of the two items of docs/specs/diversity.md,
# as the spec prints them, to ten digits.
SPECTRUM_OF_THE_SPEC = {
    "p0": [
        "85.9261619951",
        "92.9965194033",
        "106.8896328231",
        "115.5454376764",
        "120.5055938623",
        "122.9411934960",
        "124.0815127004",
        "124.2393741502",
        "123.4960764903",
        "122.4216218019",
        "60.9568756012",
    ],
    "p1": [
        "93.6948068115",
        "95.4588032593",
        "103.7262343009",
        "110.6831031479",
        "116.6729455230",
        "121.0517724083",
        "123.6021904937",
        "124.5944395085",
        "124.5403261788",
        "124.0682314877",
        "61.9071468804",
    ],
    "p2": [
        "96.3154987275",
        "101.3781441620",
        "108.0515292765",
        "114.8500079366",
        "119.7397786175",
        "121.8833113836",
        "121.8645089676",
        "120.6009754676",
        "118.9768135703",
        "117.7147383254",
        "58.6246935652",
    ],
}
FIS_OF_THE_SPEC = {
    "p0": "-0.0237536998",
    "p1": "-0.0258924700",
    "p2": "-0.0247472838",
}
NUM_VARS_OF_THE_SPEC = 1200


def read_panel():
    """The genotypes of the panel and the individuals of each population.

    The genotypes are an array of the variants by the individuals by the
    ploidy, with -1 for an allele that was not called, and the populations a
    dict of name to the indices of its individuals, in the order of the file.
    """
    with gzip.open(PANEL_VCF, "rt") as handle:
        rows = []
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                individuals = fields[9:]
                continue
            calls = []
            for call in fields[9:]:
                alleles = call.split(":")[0].replace("|", "/").split("/")
                if len(alleles) != PLOIDY:
                    raise SystemExit(f"{PANEL_VCF} has a call of {call}")
                calls.append([-1 if a == "." else int(a) for a in alleles])
            rows.append(calls)
    gts = numpy.array(rows, dtype=numpy.int16)
    pops = {}
    with open(PANEL_POPS) as handle:
        for line in handle:
            individual, pop = line.split()
            pops.setdefault(pop, []).append(individuals.index(individual))
    return gts, {pop: pops[pop] for pop in sorted(pops)}


def count_alleles(gts, individuals):
    """The counts of allele 0 and allele 1 at each variant of one population.

    It returns an array of the variants by two, and refuses a variant of more
    than two alleles: the panel is biallelic and the reference programs here
    are given two counts per variant.
    """
    of_the_pop = gts[:, individuals, :]
    if of_the_pop.max() > 1:
        raise SystemExit(f"{PANEL_VCF} holds a variant of more than two alleles")
    counts = numpy.empty((of_the_pop.shape[0], 2), dtype=numpy.int64)
    for allele in (0, 1):
        counts[:, allele] = (of_the_pop == allele).sum(axis=(1, 2))
    return counts


def make_environment(name, requirement):
    """A Python 3.12 environment holding one pinned package, and its python.

    It is made with `uv venv` under WORK_DIR when it is not there already, and
    a directory left by a run that failed half way through is made again.
    """
    directory = WORK_DIR / name
    package, version = requirement.split("==")
    python = directory / "bin" / "python"
    if python.exists() and installed_version(python, package) == version:
        return python
    if directory.exists():
        shutil.rmtree(directory)
    run(["uv", "venv", "--no-project", "--python", REFERENCE_PYTHON, str(directory)])
    run(["uv", "pip", "install", requirement], VIRTUAL_ENV=str(directory))
    got = installed_version(python, package)
    if got != version:
        raise SystemExit(f"{requirement} was asked for and {package} {got} is there")
    print(f"{requirement} is in {directory}", file=sys.stderr)
    return python


def installed_version(python, package):
    """The version of a package in an environment, or None when it is not there."""
    ask = f"import importlib.metadata as m; print(m.version({package!r}))"
    found = subprocess.run(
        [str(python), "-c", ask],
        capture_output=True,
        text=True,
        check=False,
    )
    return found.stdout.strip() if found.returncode == 0 else None


def run(command, **environment):
    """One command, with its output shown only when it fails."""
    done = subprocess.run(
        command,
        capture_output=True,
        text=True,
        env={**os.environ, **environment} if environment else None,
        check=False,
    )
    if done.returncode != 0:
        raise SystemExit(
            f"{' '.join(command)} exited with {done.returncode}\n"
            f"{done.stdout}\n{done.stderr}"
        )
    return done.stdout


def write_program(name, source):
    """One of the two programs of a reference environment, in the work directory."""
    path = WORK_DIR / name
    path.write_text(source)
    return path


def spectrum_from_dadi(counted):
    """The folded spectrum of each population, projected by dadi.

    `counted` is the allele counts of the variants that count for each
    population. dadi reads them from a file of one line per variant and
    population, and drops by itself the variants below the draw, of which
    there are none here.
    """
    lines = ["var\tpop\tn0\tn1\tcalled"]
    for pop, counts in counted.items():
        for var, (n0, n1) in counts.items():
            lines.append(f"{var}\t{pop}\t{n0}\t{n1}\t{n0 + n1}")
    counts_path = WORK_DIR / "panel_counts.tsv"
    counts_path.write_text("\n".join(lines) + "\n")
    python = make_environment("dadi-3.12", DADI)
    program = write_program("dadi_spectrum.py", DADI_PROGRAM)
    printed = run(
        [str(python), str(program), str(counts_path), str(NUM_CALLED_ALLELES)]
    )
    spectrum = {pop: [] for pop in counted}
    for line in printed.splitlines():
        pop, rarer, value = line.split()
        if int(rarer) != len(spectrum[pop]):
            raise SystemExit(f"dadi printed the bins of {pop} out of order")
        spectrum[pop].append(float(value))
    return spectrum


def fis_from_allel(gts, pops, counted):
    """F_IS of each population in its plain form, from scikit-allel.

    scikit-allel reads the genotypes of each population at the variants that
    count for it, one numpy file per population, and gives the observed and
    the expected heterozygosity per variant, from whose means the value comes.
    """
    arguments = []
    for pop, individuals in pops.items():
        path = WORK_DIR / f"panel_gts_{pop}.npy"
        variants = sorted(counted[pop])
        numpy.save(path, gts[variants][:, individuals, :])
        arguments.append(f"{pop}={path}")
    python = make_environment("allel-3.12", ALLEL)
    program = write_program("allel_fis.py", ALLEL_PROGRAM)
    # PYTHON_GIL=0 keeps quiet the warning that importing scikit-allel raises on
    # a free threading Python. This one is 3.12, which has no such build, and
    # the variable is harmless there.
    printed = run([str(python), str(program), *arguments], PYTHON_GIL="0")
    fis = {}
    for line in printed.splitlines():
        pop, value = line.split()
        fis[pop] = float(value)
    return fis


def write_spectrum(spectrum, pops):
    """The spectrum beside this script, one row per count of the rarer allele."""
    with open(HERE / "panel_folded_sfs_dadi.tsv", "w") as handle:
        handle.write("rarer_allele\t" + "\t".join(pops) + "\n")
        for rarer in range(NUM_CALLED_ALLELES // 2 + 1):
            values = "\t".join(repr(spectrum[pop][rarer]) for pop in pops)
            handle.write(f"{rarer}\t{values}\n")


def write_fis(fis, pops):
    """The three values beside this script, one row per population."""
    with open(HERE / "panel_fis_plain_allel.tsv", "w") as handle:
        handle.write("pop\tfis_plain\n")
        for pop in pops:
            handle.write(f"{pop}\t{fis[pop]!r}\n")


def check(spectrum, fis, counted):
    """What was got against the literals of docs/specs/diversity.md."""
    for pop, of_the_spec in SPECTRUM_OF_THE_SPEC.items():
        num_vars = len(counted[pop])
        if num_vars != NUM_VARS_OF_THE_SPEC:
            raise SystemExit(f"{num_vars} variants count for {pop}, not 1200")
        for rarer, expected in enumerate(of_the_spec):
            got = f"{spectrum[pop][rarer]:.10f}"
            if got != expected:
                raise SystemExit(
                    f"bin {rarer} of {pop} is {got} and the spec has {expected}"
                )
        total = sum(spectrum[pop])
        if abs(total - num_vars) > 1e-6:
            raise SystemExit(f"the spectrum of {pop} sums to {total}, not {num_vars}")
    for pop, expected in FIS_OF_THE_SPEC.items():
        got = f"{fis[pop]:.10f}"
        if got != expected:
            raise SystemExit(f"F_IS of {pop} is {got} and the spec has {expected}")


if __name__ == "__main__":
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    gts, pops = read_panel()
    # A variant counts for a population when the population called
    # min_num_individuals genotypes there, measured as its called alleles over
    # the ploidy, which is the rule of "Its Python function" of the spec.
    counted = {}
    for pop, individuals in pops.items():
        counts = count_alleles(gts, individuals)
        called = counts.sum(axis=1)
        enough = called >= MIN_NUM_INDIVIDUALS * PLOIDY
        counted[pop] = {
            int(var): (int(counts[var, 0]), int(counts[var, 1]))
            for var in numpy.flatnonzero(enough)
        }
    spectrum = spectrum_from_dadi(counted)
    fis = fis_from_allel(gts, pops, counted)
    check(spectrum, fis, counted)
    write_spectrum(spectrum, list(pops))
    write_fis(fis, list(pops))
    print("done", file=sys.stderr)
