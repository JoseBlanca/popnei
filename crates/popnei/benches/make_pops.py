"""It writes the two files that say which population each individual of
`big.vars` belongs to, one of 3 populations and one of 20.

"Speed" of docs/specs/dists.md asks for the distances between populations
to be timed on 100000 variants x 1000 individuals cut into 3 populations
and into 20, so that how the cost grows with the square of the populations
is measured and not guessed: 3 populations make 3 pairs and 20 make 190.
The dataset has no such file, and this writes both, for the benchmark of
the core, for `time_pop_dists.py`, for `time_pop_dists_pynei.py` and for
`js/popnei/bench/time_pop_dists.mjs`, which all read the same two.

A file holds one line for each individual, the name of the individual and
the name of its population with a tab between them, in the order the file
of genotypes gives the individuals. The names of the individuals are the
`s000` to `s999` that `make_big_vcf.py` writes.

The 3 populations are the ones the dataset was simulated with.
`make_big_vcf.py` draws 250 families of 4 individuals and gives each family
one of 3 populations, with an F_ST between them of 0.1, and that draw is
the first the generator makes and is of the families alone, so seeding a
generator the same way and asking it for the same 250 numbers gives the
populations of that file back whatever the number of variants it holds.
The 20 populations are runs of 50 individuals in the order of the file,
`s000` to `s049` and so on, which are of no simulated population: they cut
the same individuals into more groups, which is what the growth with the
pairs is measured over.

    uv run --no-project --with numpy python make_pops.py <directory>

It writes `pops3.tsv` and `pops20.tsv` there, and refuses to write over a
file that is already one of the two.
"""

import pathlib
import sys

import numpy

# The seed, the individuals, the families and the populations of
# `make_big_vcf.py`, which the 3 populations are the draw of.
SEED = 42
NUM_INDIVIDUALS = 1000
FAMILY_SIZE = 4
NUM_SIMULATED_POPS = 3

# How many individuals one of the 20 populations holds, which is the 1000
# of the file cut into 20 runs of the order it gives them in.
NUM_MANY_POPS = 20


def names_of_the_individuals() -> list[str]:
    """The names `make_big_vcf.py` writes, `s000` to `s999`."""
    return [f"s{number:03d}" for number in range(NUM_INDIVIDUALS)]


def the_simulated_pops() -> list[str]:
    """The population of each individual as `make_big_vcf.py` drew it: one
    of 3 for each of the 250 families, which its 4 individuals share."""
    rng = numpy.random.default_rng(SEED)
    num_families = NUM_INDIVIDUALS // FAMILY_SIZE
    family_pops = rng.integers(0, NUM_SIMULATED_POPS, size=num_families)
    return [f"pop{pop}" for pop in numpy.repeat(family_pops, FAMILY_SIZE)]


def the_many_pops() -> list[str]:
    """The population of each individual when the 1000 of the file are cut
    into 20 runs of 50 in the order the file gives them."""
    per_pop = NUM_INDIVIDUALS // NUM_MANY_POPS
    return [f"pop{number // per_pop:02d}" for number in range(NUM_INDIVIDUALS)]


def write_the_file(path: pathlib.Path, individuals: list[str], pops: list[str]) -> None:
    """It writes one line for each individual, its name and its population
    with a tab between them.

    A file that is there already is not written over: the two files are read
    by four benchmarks, and one of them silently rewritten is a timing of
    populations nobody named.
    """
    if path.exists():
        raise SystemExit(f"{path} is there already")
    with path.open("w", encoding="utf-8") as out:
        for individual, pop in zip(individuals, pops, strict=True):
            out.write(f"{individual}\t{pop}\n")


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 1
    where = pathlib.Path(sys.argv[1])
    individuals = names_of_the_individuals()
    for name, pops in (
        ("pops3.tsv", the_simulated_pops()),
        ("pops20.tsv", the_many_pops()),
    ):
        write_the_file(where / name, individuals, pops)
        counts = {pop: pops.count(pop) for pop in sorted(set(pops))}
        print(f"{where / name}: {len(counts)} populations, {counts}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
