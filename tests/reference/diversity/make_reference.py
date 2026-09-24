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
They are there so that the numbers stored here come from the one version of
each program the spec names, dadi 2.4.4 and scikit-allel 1.3.13, whatever
the development dependencies of pyproject.toml later hold. Neither program
is among those dependencies and no test of popnei imports either, so
nothing else in the repository fixes their versions.

Both of them install and run on the project's Python as well, 3.14.5 with
the global interpreter lock, which .python-version pins by its patch
version. What dadi 2.4.4 does not build on is 3.14.7, the free threading
build that `uv venv --python 3.14` picks by itself on the owner's machine:
its nlopt dependency compiles from source and stops for want of cmake,
which is not installed there. Both were measured on 24 September 2026, and
docs/reports/diversity.md holds them, with the draw of 6 on which dadi
gives the same four values on 3.14.5 as on 3.12.

The environment directory also holds the program each process runs, the
allele counts dadi projects and the genotypes scikit-allel reads, all
written again on every run.

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
    dict of name to the indices of its individuals, in the order of the
    populations file. What this script cannot read as the reference programs
    need it is refused with what was wrong: a file with no header line, a line
    of the wrong number of fields, a FORMAT that is not GT alone, a call of
    another ploidy, a half called genotype, a variant of more than two alleles
    and a populations file that does not name the individuals of the VCF.
    """
    individuals = None
    rows = []
    with gzip.open(PANEL_VCF, "rt") as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                individuals = fields[9:]
                continue
            if individuals is None:
                raise SystemExit(f"{PANEL_VCF} holds a variant before its #CHROM line")
            if len(fields) != len(individuals) + 9:
                raise SystemExit(
                    f"line {number} of {PANEL_VCF} has {len(fields)} fields and the "
                    f"{len(individuals)} individuals of the header ask for "
                    f"{len(individuals) + 9}"
                )
            if fields[8] != "GT":
                raise SystemExit(
                    f"line {number} of {PANEL_VCF} has the FORMAT {fields[8]} and "
                    "this script reads GT alone"
                )
            rows.append(read_calls(fields[9:], individuals, number))
    if individuals is None:
        raise SystemExit(f"{PANEL_VCF} has no #CHROM line")
    gts = numpy.array(rows, dtype=numpy.int16)
    if gts.max() > 1:
        of_more_alleles = numpy.flatnonzero(gts.max(axis=(1, 2)) > 1)
        raise SystemExit(
            f"{PANEL_VCF} holds a variant of more than two alleles, the first of "
            f"them the variant {of_more_alleles[0]} and {len(of_more_alleles)} in "
            "all, and the two programs here are given two counts per variant"
        )
    return gts, read_pops(individuals)


def read_calls(calls, individuals, number):
    """The alleles of every individual at one variant, -1 where not called."""
    alleles_of_all = []
    for individual, call in zip(individuals, calls, strict=True):
        alleles = call.replace("|", "/").split("/")
        not_called = [allele == "." for allele in alleles]
        if len(alleles) != PLOIDY:
            raise SystemExit(
                f"{individual} has the call {call} at line {number} of {PANEL_VCF}, "
                f"of {len(alleles)} alleles where the ploidy is {PLOIDY}"
            )
        if any(not_called) and not all(not_called):
            raise SystemExit(
                f"{individual} has the half called genotype {call} at line {number} "
                f"of {PANEL_VCF}: a half call gives its one called allele to the "
                "counts of its population, and this script does not"
            )
        alleles_of_all.append(
            [-1 if allele == "." else int(allele) for allele in alleles]
        )
    return alleles_of_all


def read_pops(individuals):
    """The indices of the individuals of each population, by population name.

    A populations file that does not name every individual of the VCF exactly
    once is refused: a line missing from it would leave that individual out of
    its population and change every number here without a word.
    """
    of_the_pop = {}
    named = []
    with open(PANEL_POPS) as handle:
        for number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) != 2:
                raise SystemExit(
                    f"line {number} of {PANEL_POPS} has {len(fields)} fields and an "
                    "individual with its population was expected"
                )
            individual, pop = fields
            named.append(individual)
            of_the_pop.setdefault(pop, []).append(individual)
    if sorted(named) != sorted(individuals):
        not_named = sorted(set(individuals) - set(named))
        not_of_the_vcf = sorted(set(named) - set(individuals))
        raise SystemExit(
            f"{PANEL_POPS} names {len(named)} individuals for the "
            f"{len(individuals)} of {PANEL_VCF}: it leaves out {len(not_named)} of "
            f"them, among them {not_named[:3]}, and names {len(not_of_the_vcf)} "
            f"that the VCF does not hold, among them {not_of_the_vcf[:3]}"
        )
    index_of = {name: index for index, name in enumerate(individuals)}
    return {
        pop: [index_of[name] for name in of_the_pop[pop]] for pop in sorted(of_the_pop)
    }


def count_alleles(gts, individuals):
    """The counts of allele 0 and allele 1 at each variant of one population.

    It returns an array of the variants by two. read_panel has refused a
    variant of more than two alleles over the whole dataset, so the two counts
    hold every allele the population called.
    """
    of_the_pop = gts[:, individuals, :]
    counts = numpy.empty((of_the_pop.shape[0], 2), dtype=numpy.int64)
    for allele in (0, 1):
        counts[:, allele] = (of_the_pop == allele).sum(axis=(1, 2))
    return counts


def make_environment(name, requirement):
    """An environment of one pinned package, and the python that runs in it.

    It is made with `uv venv` under WORK_DIR, in a directory named after the
    package and REFERENCE_PYTHON, and it is made again when the interpreter
    there is of another version, when the package is of another version or is
    not there at all, and when a run that failed half way through left the
    directory broken. The interpreter is asked for its own version rather than
    read off the name of the directory: it is part of where the numbers stored
    here come from, as the head of this file says, so a stale one would make
    that false without a word.
    """
    directory = WORK_DIR / f"{name}-{REFERENCE_PYTHON}"
    package, version = requirement.split("==")
    python = directory / "bin" / "python"
    wanted = (REFERENCE_PYTHON, version)
    if what_is_installed(python, package) == wanted:
        return python
    if directory.exists():
        shutil.rmtree(directory)
    run(
        ["uv", "venv", "--no-project", "--python", REFERENCE_PYTHON, str(directory)],
        f"the Python {REFERENCE_PYTHON} environment of {requirement}",
    )
    run(
        ["uv", "pip", "install", requirement],
        f"the install of {requirement}",
        VIRTUAL_ENV=str(directory),
    )
    got = what_is_installed(python, package)
    if got is None:
        raise SystemExit(f"{directory} was built and {package} is not there")
    if got != wanted:
        raise SystemExit(
            f"{directory} was asked for Python {REFERENCE_PYTHON} with {requirement} "
            f"and holds Python {got[0]} with {package} {got[1]}"
        )
    print(f"{requirement} is in {directory}", file=sys.stderr)
    return python


def what_is_installed(python, package):
    """The version of an interpreter and of one package in it, as two strings.

    ("3.12", "2.4.4") for the environment of dadi, and None when the
    interpreter is not there or the package cannot be imported by it.
    """
    if not python.exists():
        return None
    ask = (
        "import sys, importlib.metadata as metadata; "
        "print('.'.join(str(part) for part in sys.version_info[:2]), "
        f"metadata.version({package!r}))"
    )
    found = subprocess.run(
        [str(python), "-c", ask],
        capture_output=True,
        text=True,
        check=False,
    )
    if found.returncode != 0:
        return None
    printed = found.stdout.split()
    return (printed[0], printed[1]) if len(printed) == 2 else None


def run(command, doing, **environment):
    """One command, with its output shown only when it fails.

    `doing` is what the command is for, and it is what a failure and a command
    that is not installed are reported as.
    """
    try:
        done = subprocess.run(
            command,
            capture_output=True,
            text=True,
            env={**os.environ, **environment} if environment else None,
            check=False,
        )
    except FileNotFoundError as not_found:
        raise SystemExit(
            f"{command[0]} was not found in the PATH, and {doing} needs it"
        ) from not_found
    if done.returncode != 0:
        raise SystemExit(
            f"{doing} failed: {' '.join(command)} exited with {done.returncode}\n"
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
    python = make_environment("dadi", DADI)
    program = write_program("dadi_spectrum.py", DADI_PROGRAM)
    printed = run(
        [str(python), str(program), str(counts_path), str(NUM_CALLED_ALLELES)],
        f"the spectrum of {DADI}",
    )
    spectrum = {pop: [] for pop in counted}
    for line in printed.splitlines():
        printed_fields = line.split()
        if len(printed_fields) != 3 or printed_fields[0] not in spectrum:
            raise SystemExit(f"{DADI} printed {line!r}, which is not a bin")
        pop, rarer, value = printed_fields
        if int(rarer) != len(spectrum[pop]):
            raise SystemExit(
                f"{DADI} printed the bin {rarer} of {pop} where the bin "
                f"{len(spectrum[pop])} was expected"
            )
        spectrum[pop].append(float(value))
    num_bins = NUM_CALLED_ALLELES // 2 + 1
    for pop, bins in spectrum.items():
        if len(bins) != num_bins:
            raise SystemExit(
                f"{DADI} printed {len(bins)} bins of {pop} and a draw of "
                f"{NUM_CALLED_ALLELES} has {num_bins}"
            )
        # Each variant that counted gives one whole variant to the bins of its
        # population, whatever the draw makes of it, so the column sums to the
        # number of them. The error of dadi's own arithmetic over the 1200
        # variants of the panel is 7.0e-11.
        num_vars = len(counted[pop])
        if abs(sum(bins) - num_vars) > 1e-6:
            raise SystemExit(
                f"the spectrum of {pop} sums to {sum(bins)} and {num_vars} variants "
                "counted for it"
            )
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
    python = make_environment("allel", ALLEL)
    program = write_program("allel_fis.py", ALLEL_PROGRAM)
    printed = run([str(python), str(program), *arguments], f"F_IS of {ALLEL}")
    fis = {}
    for line in printed.splitlines():
        printed_fields = line.split()
        if (
            len(printed_fields) != 2
            or printed_fields[0] not in pops
            or printed_fields[0] in fis
        ):
            raise SystemExit(f"{ALLEL} printed {line!r}, which is not one value")
        fis[printed_fields[0]] = float(printed_fields[1])
    if set(fis) != set(pops):
        raise SystemExit(
            f"{ALLEL} gave F_IS for {sorted(fis)} and the panel holds {sorted(pops)}"
        )
    return fis


def write_table(name, lines):
    """One of the two stored files, written whole and then put in its place.

    The lines go to a file beside it which is renamed over it, so that a
    failure here leaves what was stored as it was and not half written.
    """
    path = HERE / name
    written = path.with_name(path.name + ".new")
    written.write_text("".join(lines))
    written.replace(path)


def write_spectrum(spectrum, pops):
    """The spectrum beside this script, one row per count of the rarer allele."""
    lines = ["rarer_allele\t" + "\t".join(pops) + "\n"]
    for rarer in range(NUM_CALLED_ALLELES // 2 + 1):
        values = "\t".join(repr(spectrum[pop][rarer]) for pop in pops)
        lines.append(f"{rarer}\t{values}\n")
    write_table("panel_folded_sfs_dadi.tsv", lines)


def write_fis(fis, pops):
    """The three values beside this script, one row per population."""
    lines = ["pop\tfis_plain\n"]
    lines.extend(f"{pop}\t{fis[pop]!r}\n" for pop in pops)
    write_table("panel_fis_plain_allel.tsv", lines)


def check(spectrum, fis, counted, pops):
    """What was got against the literals of docs/specs/diversity.md.

    The populations of the panel are the ones run over, and a panel whose
    populations are not those the spec has numbers for is refused: one the
    spec does not name would be written to the stored files and compared with
    nothing.
    """
    of_the_spec = sorted(SPECTRUM_OF_THE_SPEC)
    if sorted(pops) != of_the_spec or sorted(FIS_OF_THE_SPEC) != of_the_spec:
        raise SystemExit(
            f"the panel holds the populations {sorted(pops)} and the spec has a "
            f"spectrum for {of_the_spec} and F_IS for {sorted(FIS_OF_THE_SPEC)}"
        )
    for pop in pops:
        num_vars = len(counted[pop])
        if num_vars != NUM_VARS_OF_THE_SPEC:
            raise SystemExit(
                f"{num_vars} variants count for {pop} and the spec has "
                f"{NUM_VARS_OF_THE_SPEC}"
            )
        if len(SPECTRUM_OF_THE_SPEC[pop]) != len(spectrum[pop]):
            raise SystemExit(
                f"the spec has {len(SPECTRUM_OF_THE_SPEC[pop])} bins of {pop} and "
                f"dadi gave {len(spectrum[pop])}"
            )
        for rarer, expected in enumerate(SPECTRUM_OF_THE_SPEC[pop]):
            got = f"{spectrum[pop][rarer]:.10f}"
            if got != expected:
                raise SystemExit(
                    f"bin {rarer} of {pop} is {got} and the spec has {expected}"
                )
        got = f"{fis[pop]:.10f}"
        if got != FIS_OF_THE_SPEC[pop]:
            raise SystemExit(
                f"F_IS of {pop} is {got} and the spec has {FIS_OF_THE_SPEC[pop]}"
            )


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
    check(spectrum, fis, counted, list(pops))
    write_spectrum(spectrum, list(pops))
    write_fis(fis, list(pops))
    print("done", file=sys.stderr)
