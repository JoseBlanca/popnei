"""How long a whole pass of pyNei over a VCF takes, with the filter by the
rate of missing genotypes and without it.

It is what task 4.1 of docs/plans/filters.md timed pyNei with, against the
same two passes of popnei that `filter_vars.rs` times.
`docs/reports/filters-measurement.md` has the numbers and the load averages
they were taken at, and says what pyNei's pass does that popnei's does not.

A pass builds a `Variants` over the VCF with `vars_from_vcf`, puts
`filter_by_missing_data` on it when there is a threshold, walks every chunk
that `iter_vars_chunks` gives and adds up the variants of each. A chunk is
what popnei calls a block, the genotypes of a number of variants at a time.
Every run builds its own `Variants`, so that no chunk is read twice.

Two ways of running it:

    uv run python time_pynei.py <path to a VCF> <runs> [<threshold>]
    uv run python time_pynei.py <path to a VCF> <runs> <threshold> --alternating

The first times one pass, with the filter when a threshold is given and
without it when none is. Two invocations of it, one after the other, are a
pair, which is how popnei and bcftools were timed here.

The second times both passes in one process, one of each in turn. A pass of
pyNei takes about 14 s, so five of them take a minute and a quarter, and
over the ten minutes four pairs took, the pass with no filter drifted by
0.8 s while the machine grew quieter, which is more than the filter costs:
two pairs taken that way disagreed on the sign of the difference. Taken in
turn, the runs of both passes are spread over the same minutes and the
drift falls on both.

Both ways run one pass before the timed ones whose time is not taken, for
the page cache and for whatever pyNei imports on its first chunk.

pyNei is a development dependency of popnei at the commit `pyproject.toml`
names, so `uv run` is what has it; `--no-project` would not.
"""

import statistics
import sys
import time

from pynei import vars_from_vcf
from pynei.var_filters import filter_by_missing_data


def one_pass(path: str, threshold: float | None) -> int:
    """The variants of every chunk of one whole pass over the VCF at `path`,
    with the filter at `threshold` on it when there is one."""
    variants = vars_from_vcf(path)
    if threshold is not None:
        variants = filter_by_missing_data(variants, threshold)
    return sum(chunk.gts.gt_values.shape[0] for chunk in variants.iter_vars_chunks())


def said_about(what: str, times: list[float]) -> str:
    """The best, the median and the worst of `times`, under the name of the
    pass they are of."""
    return (
        f"{what}: best {min(times):.3f} s, median {statistics.median(times):.3f} s, "
        f"worst {max(times):.3f} s"
    )


def one_pass_at_a_time(path: str, runs: int, threshold: float | None) -> None:
    """`runs` timed passes of the one kind, and one before them that is not
    timed."""
    what = (
        "no filter" if threshold is None else f"the missing data filter at {threshold}"
    )
    print(f"{path}, pyNei, {what}, {runs} runs")
    started = time.perf_counter()
    kept = one_pass(path, threshold)
    print(
        f"the first run, which is not timed: {time.perf_counter() - started:.3f} s, "
        f"{kept} variants"
    )
    times = []
    for run in range(1, runs + 1):
        started = time.perf_counter()
        kept = one_pass(path, threshold)
        took = time.perf_counter() - started
        times.append(took)
        print(f"run {run}: {took:.3f} s, {kept} variants")
    print(said_about(what, times))


def the_two_in_turn(path: str, runs: int, threshold: float) -> None:
    """`runs` timed passes of each kind, one of each in turn, and one before
    them that is not timed."""
    print(
        f"{path}, pyNei, no filter and the missing data filter at {threshold} in turn, "
        f"{runs} runs of each"
    )
    one_pass(path, None)
    print("the first pass, which is not timed, is done")
    without: list[float] = []
    with_it: list[float] = []
    for run in range(1, runs + 1):
        for what, threshold_of_the_pass, times in (
            ("no filter", None, without),
            (f"the filter at {threshold}", threshold, with_it),
        ):
            started = time.perf_counter()
            kept = one_pass(path, threshold_of_the_pass)
            took = time.perf_counter() - started
            times.append(took)
            print(f"run {run}, {what}: {took:.3f} s, {kept} variants")
    print(said_about("no filter", without))
    print(said_about("the filter", with_it))
    costs = statistics.median(with_it) - statistics.median(without)
    print(f"the filter costs {costs:.3f} s on the medians")


def main() -> int:
    arguments = sys.argv[1:]
    alternating = "--alternating" in arguments
    if alternating:
        arguments.remove("--alternating")
    if not 2 <= len(arguments) <= 3:
        print(__doc__)
        return 1
    path = arguments[0]
    runs = int(arguments[1])
    threshold = float(arguments[2]) if len(arguments) > 2 else None
    if alternating:
        if threshold is None:
            print("--alternating needs the threshold of the pass with the filter")
            return 1
        the_two_in_turn(path, runs, threshold)
        return 0
    one_pass_at_a_time(path, runs, threshold)
    return 0


if __name__ == "__main__":
    sys.exit(main())
