"""How long the Kosman distances of every pair take over a whole file, and
how long reading that file alone takes.

It is what task 3.1 of docs/plans/dists-kosman.md timed popnei with, against
the three numbers of "Speed" of docs/specs/dists.md, which are of the
calculation with the reading of the file taken out.
`docs/reports/dists-kosman-measurement.md` has the numbers and the load
averages they were taken at. `time_kosman_pynei.py`, beside this file, is
the same clock on pyNei, and `js/popnei/bench/time_kosman_dists.mjs` the
one on the wasm build under node.

It times two passes over the same source. The **calculation** is one call of
`calc_pairwise_kosman_dists`, which reads every block of the source and
gives the distance of every pair. The **reading alone** is
`iter_blocks(fields=())`, the same pass with the genotypes as the only field
the blocks carry, which is what the calculation asks its reader for, with
nothing done to the genotypes but adding up how many alleles came out. The
difference of the two is what the calculation costs beyond the reader, and
it is the number the targets of the spec are of.

`open_vcf` and `open_vars` are outside both timings. They read the header of
the VCF, or the schema and the footer of the vars file, and no genotype:
every pass reads the source again from its start, so the handle is built
once and given to every run. What the open itself took is printed under the
runs.

The threads are rayon's, which the library never sets: the pool takes one
thread for each core unless `RAYON_NUM_THREADS` is in the environment of the
process before the first call, which is why it is set in the command and not
from Python. Over a vars file the reader runs on the thread that calls it,
so there the threads are the calculation's alone.

    RAYON_NUM_THREADS=1 uv run python time_kosman_dists.py <path> <runs>
    uv run python time_kosman_dists.py <path> <runs>

A path that ends in `.vars` is opened with `open_vars` and anything else
with `open_vcf`. One pass of each kind runs before the timed ones and is not
timed, for the page cache. The two kinds are run one after the other inside
each round, so that a machine that grows busier over the runs falls on both.

popnei has to be built in release, `uv run maturin develop --release`: the
debug build is several times slower and says nothing about the targets.
"""

import os
import statistics
import sys
import time

from popnei import calc_pairwise_kosman_dists, open_vars, open_vcf


def the_calculation(variants) -> tuple[int, int]:
    """One call of `calc_pairwise_kosman_dists` over `variants`, as the
    distances of its pairs and the variants the pass took."""
    distances = calc_pairwise_kosman_dists(variants)
    return distances.dist_vector.shape[0], distances.pass_stats.num_vars


def the_reading_alone(variants) -> tuple[int, int]:
    """One pass over `variants` with the genotypes as the only field, as the
    alleles that came out of it and the variants the pass took.

    The alleles are added up so that a reader which stopped filling the
    genotypes, or filled them only when somebody read them, would show:
    without that count the pass would be of a column nobody touched.
    """
    blocks = variants.iter_blocks(fields=())
    alleles = sum(int(block.gts.size) for block in blocks)
    return alleles, blocks.pass_stats.num_vars


def said_about(what: str, times: list[float]) -> str:
    """The best, the median and the worst of `times`, under the name of the
    pass they are of."""
    return (
        f"{what}: best {min(times):.3f} s, median {statistics.median(times):.3f} s, "
        f"worst {max(times):.3f} s"
    )


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 1
    path = sys.argv[1]
    runs = int(sys.argv[2])
    threads = os.environ.get("RAYON_NUM_THREADS", "one for each core")
    print(f"{path}, popnei, {runs} runs, RAYON_NUM_THREADS={threads}")

    started = time.perf_counter()
    variants = open_vars(path) if path.endswith(".vars") else open_vcf(path)
    opened_in = time.perf_counter() - started

    print(
        f"the first run of each kind, which is not timed: {the_calculation(variants)}"
    )
    the_reading_alone(variants)

    calculations: list[float] = []
    readings: list[float] = []
    for run in range(1, runs + 1):
        for what, pass_of_the_run, times in (
            ("the calculation", the_calculation, calculations),
            ("the reading alone", the_reading_alone, readings),
        ):
            started = time.perf_counter()
            counted, num_vars = pass_of_the_run(variants)
            took = time.perf_counter() - started
            times.append(took)
            print(f"run {run}, {what}: {took:.3f} s, {num_vars} variants, {counted}")
    print(said_about("the calculation", calculations))
    print(said_about("the reading alone", readings))
    print(
        f"the calculation with the reading taken out, on the bests: "
        f"{min(calculations) - min(readings):.3f} s"
    )
    print(f"opening the source, which is in neither: {opened_in:.4f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
