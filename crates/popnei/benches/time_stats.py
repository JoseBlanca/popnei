"""How long a whole pass of popnei over a vars file or a VCF takes, for each
of the two calculations of the stats module and for the read alone.

It is what task 6.1 of docs/plans/stats.md timed popnei with, against the
same passes of pyNei that `time_pynei_stats.py` times.
`docs/reports/stats-measurement.md` has the numbers and the load averages
they were taken at, and says what each pass does.

The six passes it can time, one per invocation:

    read              `iter_blocks` with the genotypes as the only field,
                      which is the read of the file and nothing else
    per-var           `calc_per_var_distribs` with the five statistics and
                      no `pops`
    per-var-pops      the same with populations: the four of 250 individuals
                      of the report, the individuals of the file in the order
                      they are in it, 250 to each, or the populations of a
                      file given as a fourth argument
    per-var-obs-het   the same with the observed heterozygosity alone and
                      no `pops`
    per-var-maf       the same with the major allele frequency alone and no
                      `pops`
    per-individual    `calc_per_individual_stats`

The two passes with one statistic say how the pass with five divides up,
and they are the rows pyNei's table of "Speed" of `docs/specs/stats.md`
has for its own `obs_het` alone and `maf` alone.

Every run opens the file again, so that no block is read twice and every run
pays the opening. A path that ends in `.vars` is opened with `open_vars` and
anything else with `open_vcf`, which is how `time_pop_dists.py` and
`time_diversity.py` open theirs: the numbers of
`docs/reports/stats-measurement.md` are all over a vars file, and task 4.1 of
docs/plans/diversity.md timed this pass over a VCF as well, to read it
against `calc_pop_diversity` over the same file.

The threads are rayon's, which popnei takes from the environment because it
builds no pool of its own: `RAYON_NUM_THREADS=1` for one thread and
`RAYON_NUM_THREADS=18` for the 18 cores of the machine.

    RAYON_NUM_THREADS=1 uv run python time_stats.py <path> <what> <runs> \
        [<populations file>]

The fourth argument is a file of one line for each individual: the name of
the individual, a tab, and the name of its population, which is the shape of
`tests/reference/stats/panel_pops_bcftools.txt` and of the files
`make_pops.py` writes. `per-var-pops` then runs over those populations
instead of the four of 250, which is what task 4.1 of
docs/plans/diversity.md timed `calc_per_var_distribs` with, so that it and
`calc_pop_diversity` of `time_diversity.py` count the same populations of the
same file. The other passes take no populations and refuse the argument.

One pass before the timed ones is not timed: it pays the page faults of the
first touch of the memory a pass works in, which a process pays once, and it
reads the file, so that the timed runs read it from the page cache.

It prints the wall time of each run with the variants the pass gave, and
then the best, the median and the worst. The best is what the report takes,
since the machine is not idle and what it is doing can only make a run
longer.
"""

import os
import statistics
import sys
import time

import popnei

WHATS = (
    "read",
    "per-var",
    "per-var-pops",
    "per-var-obs-het",
    "per-var-maf",
    "per-individual",
)


def four_pops(variants: popnei.Variants) -> dict[str, list[str]]:
    """Four populations of 250 individuals, the individuals of the file in
    the order they are in it, 250 to each."""
    individuals = list(variants.individuals)
    if len(individuals) != 1000:
        raise ValueError(
            f"the four populations of 250 are of a file of 1000 individuals, "
            f"and this one has {len(individuals)}"
        )
    return {
        f"pop{number}": individuals[number * 250 : (number + 1) * 250]
        for number in range(4)
    }


def the_pops(path: str) -> dict[str, list[str]]:
    """The populations of the file at `path`, as the dict of population name
    to the names of its individuals that `calc_per_var_distribs` takes.

    The file holds one line for each individual, its name, a tab and the name
    of its population. The populations come out in the order in which the
    file first names each of them, and the individuals of each in the order
    of the file.
    """
    pops: dict[str, list[str]] = {}
    with open(path) as fhand:
        for number, line in enumerate(fhand, start=1):
            if not line.strip():
                continue
            written = line.rstrip("\n")
            individual, tab, pop = written.partition("\t")
            if not tab or not individual or not pop:
                raise ValueError(
                    f"{path}, line {number}: a line of a populations file is "
                    f"the name of an individual, a tab and the name of its "
                    f"population, and this one is {written!r}"
                )
            pops.setdefault(pop, []).append(individual)
    return pops


def one_pass(path: str, what: str, pops: dict[str, list[str]] | None = None) -> int:
    """One whole pass over the vars file at `path`, and the variants it
    gave.

    `pops` is the populations `per-var-pops` runs over, and with none it is
    the four of 250 individuals of the report.
    """
    variants = (
        popnei.open_vars(path) if path.endswith(".vars") else popnei.open_vcf(path)
    )
    if what == "read":
        blocks = variants.iter_blocks(fields=())
        return sum(block.gts.shape[0] for block in blocks)
    if what == "per-var":
        return popnei.calc_per_var_distribs(variants).pass_stats.num_vars
    if what == "per-var-pops":
        return popnei.calc_per_var_distribs(
            variants, pops=four_pops(variants) if pops is None else pops
        ).pass_stats.num_vars
    # `pass_stats.num_vars` is the variants of the pass whatever statistics
    # were asked for, so the two passes with one of them count what the five
    # count.
    if what == "per-var-obs-het":
        return popnei.calc_per_var_distribs(
            variants, stats=(popnei.PerVarStat.OBS_HET,)
        ).pass_stats.num_vars
    if what == "per-var-maf":
        return popnei.calc_per_var_distribs(
            variants, stats=(popnei.PerVarStat.MAF,)
        ).pass_stats.num_vars
    if what == "per-individual":
        return popnei.calc_per_individual_stats(variants).pass_stats.num_vars
    raise ValueError(f"the pass to time is one of {WHATS}, and {what!r} was given")


def main() -> int:
    arguments = sys.argv[1:]
    if len(arguments) not in (3, 4):
        print(__doc__)
        return 1
    path, what, runs = arguments[0], arguments[1], int(arguments[2])
    pops_path = arguments[3] if len(arguments) == 4 else None
    if what not in WHATS:
        print(f"the pass to time is one of {WHATS}, and {what!r} was given")
        return 1
    if pops_path is not None and what != "per-var-pops":
        print(
            f"a populations file is read by 'per-var-pops' alone, and {what!r} "
            f"was given with {pops_path!r}"
        )
        return 1
    pops = None if pops_path is None else the_pops(pops_path)
    threads = os.environ.get("RAYON_NUM_THREADS", "the cores of the machine")
    print(
        f"{path}, popnei {popnei.__version__}, {what}, {runs} runs, {threads} threads"
        + ("" if pops is None else f", {len(pops)} populations of {pops_path}")
    )
    started = time.perf_counter()
    num_vars = one_pass(path, what, pops)
    print(
        f"the first run, which is not timed: {time.perf_counter() - started:.3f} s, "
        f"{num_vars} variants"
    )
    times = []
    for run in range(1, runs + 1):
        started = time.perf_counter()
        num_vars = one_pass(path, what, pops)
        took = time.perf_counter() - started
        times.append(took)
        print(f"run {run}: {took:.3f} s, {num_vars} variants")
    print(
        f"{what}: best {min(times):.3f} s, median {statistics.median(times):.3f} s, "
        f"worst {max(times):.3f} s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
