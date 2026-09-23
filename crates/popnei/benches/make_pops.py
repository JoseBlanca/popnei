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

    uv run --no-project --with numpy python make_pops.py <directory> [n ...]
    uv run --no-project --with numpy python make_pops.py <directory> <n>x<m>

With no number it writes `pops3.tsv` and `pops20.tsv` there. With numbers
it writes `pops<n>.tsv` for each of them, the individuals cut into n runs
of the order of the file, which is what `pops20.tsv` is: the performance
review of 23 September 2026 asked for 6 and for 10 as well, to tell the
cost of one more population apart from the cost of one more pair. A number
of 3 would give a `pops3.tsv` of runs and not the simulated populations,
and the file is there already, which is refused.

`<n>x<m>` writes `pops<n>of<m>.tsv`, n populations of m individuals each,
the first n times m individuals of the file and no other. Every file above
names all 1000, so the populations and their sizes move together and what
each costs cannot be told apart; a file that names fewer individuals holds
the populations still and changes only how many genotypes are counted,
which is what separates them.

It refuses to write over a file that is there.
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


def the_pops_of_runs(num_pops: int) -> list[str]:
    """The population of each individual when the 1000 of the file are cut
    into `num_pops` runs of the order the file gives them.

    The runs are as near the same size as 1000 individuals divided that way
    allow: `num_pops` of them is 50 individuals each, and 6 is four runs of
    167 and two of 166.
    """
    digits = len(str(num_pops - 1))
    return [
        f"pop{number * num_pops // NUM_INDIVIDUALS:0{digits}d}"
        for number in range(NUM_INDIVIDUALS)
    ]


def some_of_the_pops(num_pops: int, per_pop: int) -> tuple[list[str], list[str]]:
    """`num_pops` populations of `per_pop` individuals each, taken from the
    front of the file, as the individuals named and their populations.

    The individuals of the file that are not among the first `num_pops`
    times `per_pop` are in no population and take no part in the pass.
    """
    wanted = num_pops * per_pop
    if wanted > NUM_INDIVIDUALS:
        raise SystemExit(
            f"{num_pops} populations of {per_pop} are {wanted} individuals, "
            f"and the file has {NUM_INDIVIDUALS}"
        )
    individuals = names_of_the_individuals()[:wanted]
    digits = len(str(num_pops - 1))
    pops = [f"pop{number // per_pop:0{digits}d}" for number in range(wanted)]
    return individuals, pops


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
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    where = pathlib.Path(sys.argv[1])
    individuals = names_of_the_individuals()
    if len(sys.argv) > 2:
        wanted = []
        for arg in sys.argv[2:]:
            if "x" in arg:
                num_pops, per_pop = (int(part) for part in arg.split("x", 1))
                some, their_pops = some_of_the_pops(num_pops, per_pop)
                wanted.append((f"pops{num_pops}of{per_pop}.tsv", their_pops, some))
                continue
            num_pops = int(arg)
            wanted.append((f"pops{num_pops}.tsv", the_pops_of_runs(num_pops), individuals))
    else:
        wanted = [
            ("pops3.tsv", the_simulated_pops(), individuals),
            ("pops20.tsv", the_pops_of_runs(NUM_MANY_POPS), individuals),
        ]
    for name, pops, of_the_file in wanted:
        write_the_file(where / name, of_the_file, pops)
        counts = {pop: pops.count(pop) for pop in sorted(set(pops))}
        print(
            f"{where / name}: {len(counts)} populations, "
            f"{len(of_the_file)} individuals, {counts}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
