"""How long a whole pass of popnei over a vars file takes, for each of the
two calculations of the stats module and for the read alone.

It is what task 6.1 of docs/plans/stats.md timed popnei with, against the
same passes of pyNei that `time_pynei_stats.py` times.
`docs/reports/stats-measurement.md` has the numbers and the load averages
they were taken at, and says what each pass does.

The six passes it can time, one per invocation:

    read              `iter_blocks` with the genotypes as the only field,
                      which is the read of the file and nothing else
    per-var           `calc_per_var_distribs` with the five statistics and
                      no `pops`
    per-var-pops      the same with 4 populations of 250 individuals, the
                      individuals of the file in the order they are in it,
                      250 to each population
    per-var-obs-het   the same with the observed heterozygosity alone and
                      no `pops`
    per-var-maf       the same with the major allele frequency alone and no
                      `pops`
    per-individual    `calc_per_individual_stats`

The two passes with one statistic say how the pass with five divides up,
and they are the rows pyNei's table of "Speed" of `docs/specs/stats.md`
has for its own `obs_het` alone and `maf` alone.

Every run opens the file again with `open_vars`, so that no block is read
twice and every run pays the opening.

The threads are rayon's, which popnei takes from the environment because it
builds no pool of its own: `RAYON_NUM_THREADS=1` for one thread and
`RAYON_NUM_THREADS=18` for the 18 cores of the machine.

    RAYON_NUM_THREADS=1 uv run python time_stats.py <path to a vars file> \
        <what> <runs>

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


def one_pass(path: str, what: str) -> int:
    """One whole pass over the vars file at `path`, and the variants it
    gave."""
    variants = popnei.open_vars(path)
    if what == "read":
        blocks = variants.iter_blocks(fields=())
        return sum(block.gts.shape[0] for block in blocks)
    if what == "per-var":
        return popnei.calc_per_var_distribs(variants).pass_stats.num_vars
    if what == "per-var-pops":
        return popnei.calc_per_var_distribs(
            variants, pops=four_pops(variants)
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
    if len(arguments) != 3:
        print(__doc__)
        return 1
    path, what, runs = arguments[0], arguments[1], int(arguments[2])
    if what not in WHATS:
        print(f"the pass to time is one of {WHATS}, and {what!r} was given")
        return 1
    threads = os.environ.get("RAYON_NUM_THREADS", "the cores of the machine")
    print(
        f"{path}, popnei {popnei.__version__}, {what}, {runs} runs, {threads} threads"
    )
    started = time.perf_counter()
    num_vars = one_pass(path, what)
    print(
        f"the first run, which is not timed: {time.perf_counter() - started:.3f} s, "
        f"{num_vars} variants"
    )
    times = []
    for run in range(1, runs + 1):
        started = time.perf_counter()
        num_vars = one_pass(path, what)
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
