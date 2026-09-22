"""How long a whole pass of pyNei over its own vars file takes, for each of
the two calculations that popnei's stats module mirrors and for the read
alone.

It is what task 6.1 of docs/plans/stats.md timed pyNei with, against the
same passes of popnei that `time_stats.py` times.
`docs/reports/stats-measurement.md` has the numbers and the load averages
they were taken at, and says what each pass does.

The six passes it can time, one per invocation:

    read              `iter_vars_chunks`, the chunks and nothing else
    per-var           `calc_per_var_distribs` with its four statistics and
                      no `pops`
    per-var-pops      the same with 4 populations of 250 samples, the
                      samples of the file in the order they are in it, 250
                      to each population
    per-var-obs-het   the same with the observed heterozygosity alone and
                      no `pops`
    per-var-maf       the same with the major allele frequency alone and no
                      `pops`
    per-sample        `calc_per_sample_stats`

Every run loads the file again with `load_vars`, so that no chunk is read
twice and every run pays the loading.

pyNei takes its threads as an argument and not from the environment, so the
number of threads is given here and `RAYON_NUM_THREADS` means nothing to it.

    uv run python time_pynei_stats.py <path to a pyNei vars file> <what> \
        <runs> <threads>

One pass before the timed ones is not timed: it reads the file, so that the
timed runs read it from the page cache, and it pays whatever pyNei imports
on its first chunk.

It prints the wall time of each run with the variants the pass gave, and
then the best, the median and the worst.

pyNei is a development dependency of popnei at the commit `pyproject.toml`
names, so `uv run` is what has it; `--no-project` would not.
"""

import statistics
import sys
import time

from pynei import load_vars
from pynei.per_var_stats import PerVarStat, calc_per_var_distribs
from pynei.sample_stats import calc_per_sample_stats

WHATS = (
    "read",
    "per-var",
    "per-var-pops",
    "per-var-obs-het",
    "per-var-maf",
    "per-sample",
)


def four_pops(variants) -> dict[str, list[str]]:
    """Four populations of 250 samples, the samples of the file in the order
    they are in it, 250 to each."""
    samples = list(variants.samples)
    if len(samples) != 1000:
        raise ValueError(
            f"the four populations of 250 are of a file of 1000 samples, and "
            f"this one has {len(samples)}"
        )
    return {
        f"pop{number}": samples[number * 250 : (number + 1) * 250]
        for number in range(4)
    }


def one_pass(path: str, what: str, threads: int) -> int:
    """One whole pass over the pyNei vars file at `path`, and the variants it
    gave."""
    variants = load_vars(path)
    if what == "read":
        return sum(
            chunk.gts.gt_values.shape[0] for chunk in variants.iter_vars_chunks()
        )
    if what == "per-var":
        result = calc_per_var_distribs(variants, num_threads=threads)
        return int(result.poly_vars_ratio.tot_num_variants_with_data.iloc[0])
    if what == "per-var-pops":
        result = calc_per_var_distribs(
            variants, pops=four_pops(variants), num_threads=threads
        )
        return int(result.poly_vars_ratio.tot_num_variants_with_data.iloc[0])
    if what == "per-var-obs-het":
        result = calc_per_var_distribs(
            variants, stats=(PerVarStat.OBS_HET,), num_threads=threads
        )
        return int(result.obs_het.hist_counts.to_numpy().sum())
    if what == "per-var-maf":
        result = calc_per_var_distribs(
            variants, stats=(PerVarStat.MAF,), num_threads=threads
        )
        return int(result.maf.hist_counts.to_numpy().sum())
    if what == "per-sample":
        return calc_per_sample_stats(variants, num_threads=threads).shape[0]
    raise ValueError(f"the pass to time is one of {WHATS}, and {what!r} was given")


def main() -> int:
    arguments = sys.argv[1:]
    if len(arguments) != 4:
        print(__doc__)
        return 1
    path, what, runs, threads = (
        arguments[0],
        arguments[1],
        int(arguments[2]),
        int(arguments[3]),
    )
    if what not in WHATS:
        print(f"the pass to time is one of {WHATS}, and {what!r} was given")
        return 1
    print(f"{path}, pyNei, {what}, {runs} runs, {threads} threads")
    started = time.perf_counter()
    counted = one_pass(path, what, threads)
    print(
        f"the first run, which is not timed: {time.perf_counter() - started:.3f} s, "
        f"{counted} counted"
    )
    times = []
    for run in range(1, runs + 1):
        started = time.perf_counter()
        counted = one_pass(path, what, threads)
        took = time.perf_counter() - started
        times.append(took)
        print(f"run {run}: {took:.3f} s, {counted} counted")
    print(
        f"{what}: best {min(times):.3f} s, median {statistics.median(times):.3f} s, "
        f"worst {max(times):.3f} s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
