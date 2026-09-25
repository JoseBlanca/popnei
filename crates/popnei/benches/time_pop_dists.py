"""How long the distances between populations take over a whole vars file,
and how long reading that file alone takes.

Nothing had timed `calc_pop_dists` when this file was written. "Speed" of
docs/specs/dists.md asks for it on 100000 variants x 1000 individuals,
biallelic, 3 in 100 genotypes missing, read from a vars file, with the
individuals cut into 3 populations and into 20, and with the time of reading
that file measured beside the calculation and taken out, so that the reader
is not in the number. `time_pop_dists_pynei.py`, beside this file, is the
same clock on pyNei's `calc_jost_dest_pop_dists`, which is the program the
spec compares with. `time_kosman_dists.py` is the same clock on the
distances between individuals, and this file has its shape with another call
in it.

It times two passes over the same source. The **calculation** is one call of
`calc_pop_dists`, which reads every block of the source and gives the seven
measures of how far apart two populations are for every pair of them, each
value with the standard error of a block jackknife beside it. The **reading
alone** is `iter_blocks(fields=("chrom", "pos"))`, the same pass carrying the
genotypes, the chromosome and the position, which are the fields the
calculation asks its reader for when the groups are stretches of a
chromosome, with nothing done to them but adding up how many alleles and how
many positions each block says it holds, which is their shape and not their
content. The difference of the two is what the calculation costs beyond the
reader, and it is the number the spec asks for.

The populations come from a file of one line for each individual: the name
of the individual, a tab, and the name of its population. `make_pops.py`,
beside this file, wrote the two the spec asks for, `pops3.tsv` of 3
populations and `pops20.tsv` of 20. The populations are kept in the order in
which the file first names each of them, and that is the order of the pairs
of the result.

The fourth argument is the length in base pairs of a resampling group, the
stretch of one chromosome that the standard errors are resampled over, and
it is 1000000 when it is not given. `calc_pop_dists` has no default for it:
its `jackknife_group` is a required argument, because a length is a claim
about how far linkage reaches in the dataset at hand and popnei does not
make that claim for a user. On the dataset of the spec, two chromosomes of
50000 variants 1000 base pairs apart, 1000000 cuts the variants into 100
groups. A pass whose variants fall into fewer than 20 groups is a
`ValueError`, so a shorter file needs a shorter group: 5000 variants of that
spacing fall into 5 groups at 1000000 and into 50 at 100000.

All seven measures are calculated, which is what leaving `measures` out
asks for, and `min_num_individuals` is left at its default of 20 called
genotypes.

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

    RAYON_NUM_THREADS=1 uv run python time_pop_dists.py <path> <populations file> <runs>
    uv run python time_pop_dists.py <path> <populations file> <runs> <group length>

A path that ends in `.vars` is opened with `open_vars` and anything else
with `open_vcf`. One pass of each kind runs before the timed ones and is not
timed, for the page cache. The two kinds are run one after the other inside
each round, so that a machine that grows busier over the runs falls on both.

popnei has to be built in release, `uv run maturin develop --release`. What
`uv run maturin develop` makes is a debug build of the core, and for the
Kosman distances over the same file it was about fifty times slower, 65.9 s
against 1.275 s on one thread; nobody has measured what that build costs
here, and no number taken with it says anything about the calculation.
"""

import os
import statistics
import sys
import time

from popnei import calc_pop_dists, open_vars, open_vcf

# The length in base pairs of a resampling group when the command gives
# none. It cuts the 100000 variants of the dataset of "Speed" of
# docs/specs/dists.md, two chromosomes of variants 1000 base pairs apart,
# into 100 groups.
GROUP_LENGTH = 1_000_000


def the_pops(path: str) -> dict[str, list[str]]:
    """The populations of the file at `path`, as the dict of population name
    to the names of its individuals that `calc_pop_dists` takes.

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


def the_calculation(variants, pops, group_length: int) -> tuple[str, int]:
    """One call of `calc_pop_dists` over `variants` for the pairs of `pops`,
    with the standard errors resampled over groups of `group_length` base
    pairs, as what the result holds and the variants the pass took.

    The value printed is Hudson's F_ST of the first pair with its standard
    error, one of the seven measures the call gives, so that a run whose
    pass read nothing shows itself.
    """
    dists = calc_pop_dists(variants, pops, jackknife_group=group_length)
    fst = dists.fst
    said = (
        f"{fst.dist_vector.shape[0]} pairs, {len(dists.group_ids)} groups, "
        f"F_ST of {dists.pops[0]} and {dists.pops[1]} "
        f"{fst.dist_vector[0]:.6f} +- {fst.standard_errors[0]:.6f} over "
        f"{int(dists.num_vars[0])} variants"
    )
    return said, dists.pass_stats.num_vars


def the_reading_alone(variants) -> tuple[str, int]:
    """One pass over `variants` carrying the genotypes, the chromosome and
    the position, as the alleles and the positions its blocks say they hold
    and the variants the pass took.

    `block.gts.size` and `block.pos.size` are the shapes of the arrays the
    core filled, which reach numpy without a copy, so adding them up shows a
    pass whose blocks carry neither column and does not show a buffer of the
    right shape that nothing wrote into. Nothing here reads a genotype or a
    position.
    """
    blocks = variants.iter_blocks(fields=("chrom", "pos"))
    alleles = 0
    positions = 0
    for block in blocks:
        alleles += int(block.gts.size)
        positions += int(block.pos.size)
    return f"{alleles} alleles, {positions} positions", blocks.pass_stats.num_vars


def said_about(what: str, times: list[float]) -> str:
    """The best, the median and the worst of `times`, under the name of the
    pass they are of."""
    return (
        f"{what}: best {min(times):.3f} s, median {statistics.median(times):.3f} s, "
        f"worst {max(times):.3f} s"
    )


def main() -> int:
    if len(sys.argv) not in (4, 5):
        print(__doc__)
        return 1
    path = sys.argv[1]
    pops_path = sys.argv[2]
    runs = int(sys.argv[3])
    group_length = int(sys.argv[4]) if len(sys.argv) == 5 else GROUP_LENGTH
    pops = the_pops(pops_path)
    threads = os.environ.get("RAYON_NUM_THREADS", "one for each core")
    print(
        f"{path}, popnei, {len(pops)} populations of {pops_path}, {runs} runs, "
        f"jackknife_group={group_length}, RAYON_NUM_THREADS={threads}"
    )

    started = time.perf_counter()
    variants = open_vars(path) if path.endswith(".vars") else open_vcf(path)
    opened_in = time.perf_counter() - started

    print(
        f"the first run of each kind, which is not timed: "
        f"{the_calculation(variants, pops, group_length)}"
    )
    the_reading_alone(variants)

    calculations: list[float] = []
    readings: list[float] = []
    for run in range(1, runs + 1):
        for what, pass_of_the_run, times in (
            (
                "the calculation",
                lambda variants: the_calculation(variants, pops, group_length),
                calculations,
            ),
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
