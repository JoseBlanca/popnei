"""How long a whole pass of popnei over a vars file or a VCF takes, for the
passes of the diversity module and for the read alone.

It is what task 4.1 of docs/plans/diversity.md timed popnei with, for the
target that "Speed" of docs/specs/diversity.md leaves open.
`crates/popnei/benches/time_stats.py` is the same clock on the two
calculations of the stats module, and this file has its shape with another
call in it; `crates/popnei/benches/diversity_pass.rs` times the same pass
over blocks that are already in memory, with no file read inside the clock,
and counts the bytes it holds.

The four passes it can time, one per invocation:

    read                `iter_blocks` with the genotypes as the only field,
                        which is the read of the file and nothing else, and
                        the floor every other pass here is above
    diversity           `calc_pop_diversity` with the five statistics: the
                        alleles each population called, the private ones
                        among them, the variants that vary in it, the folded
                        site frequency spectrum and F_IS
    diversity-sfs       the same with the folded spectrum alone, which says
                        what the spectrum is of the pass with five
    diversity-no-draw   the same with the four statistics that need no draw,
                        which says what the draw adds

The populations come from a file of one line for each individual: the name
of the individual, a tab, and the name of its population, which is the shape
of `tests/reference/stats/panel_pops_bcftools.txt` and of the files
`crates/popnei/benches/make_pops.py` writes. The populations are kept in the
order in which the file first names each of them, and that is the order of
the rows of every result. The same populations are used whatever is timed,
so that two passes over one file count the same populations.

The draw is `num_called_alleles`, how many called alleles every population
is brought down to. `diversity` and `diversity-sfs` take a number, since the
spectrum is the bins of a draw and cannot be asked for without one; `read`
and `diversity-no-draw` take the word `none`, so that a run of the pass that
is meant to have no draw cannot be given one by mistake and timed as
something else.

A draw above every population's called alleles is no error and is not a
timing of the draw: no variant is in the draw for any population, every
standardized value is NaN and every bin of the spectrum is 0. The line each
run prints holds the variants in the draw for the first population, so a run
in that case shows itself.

    RAYON_NUM_THREADS=1 uv run python time_diversity.py <path> <what> \
        <populations file> <draw> <runs>

The threads are rayon's, which popnei takes from the environment because it
builds no pool of its own: `RAYON_NUM_THREADS=1` for one thread and
`RAYON_NUM_THREADS=18` for the 18 cores of the machine. Over a vars file the
reader runs on the thread that calls it, so there the threads are the
calculation's alone.

A path that ends in `.vars` is opened with `open_vars` and anything else
with `open_vcf`. Every run opens the file again, so that no block is read
twice and every run pays the opening.

Before the timed runs the file is read once, with the genotypes as the only
field, and the variants that read gives are what every timed run is checked
against: a run whose pass gives another number fails and prints both. That
read also pays the page faults of the first touch of the memory a pass works
in and brings the file into the page cache, and then one pass of what is
being timed runs and is not timed either. A pass asked for a statistic that
is `None` in the result fails the same way, so that a run pointed at the
wrong statistic cannot be reported as a timing of it.

popnei has to be built in release, `uv run maturin develop --release`. What
`uv run maturin develop` makes is a debug build of the core, and no number
taken with it says anything about the calculation.

It prints the wall time of each run with what the pass gave, and then the
best, the median and the worst. The best is what a report takes, since the
machine is not idle and what it is doing can only make a run longer.
"""

import os
import statistics
import sys
import time

import popnei

WHATS = (
    "read",
    "diversity",
    "diversity-sfs",
    "diversity-no-draw",
)

# The passes that need a draw, and so refuse a `none` where the number goes.
# The other two refuse a number, so that a pass meant to run without a draw
# cannot be timed with one.
WITH_A_DRAW = ("diversity", "diversity-sfs")

# What is written in the place of the draw for a pass that has none.
NO_DRAW = "none"

# The statistics each pass asks for. `read` asks for none, doing no
# calculation at all.
STATS = {
    "diversity": tuple(popnei.PopDiversityStat),
    "diversity-sfs": (popnei.PopDiversityStat.FOLDED_SFS,),
    "diversity-no-draw": popnei.PopDiversityStat.WITHOUT_A_DRAW,
}


def the_pops(path: str) -> dict[str, list[str]]:
    """The populations of the file at `path`, as the dict of population name
    to the names of its individuals that `calc_pop_diversity` takes.

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


def open_the_file(path: str):
    """The variants of the file at `path`, read as a vars file when its name
    ends in `.vars` and as a VCF otherwise."""
    return popnei.open_vars(path) if path.endswith(".vars") else popnei.open_vcf(path)


def the_read_alone(path: str) -> tuple[int, str]:
    """One pass over the file at `path` carrying the genotypes and no other
    field, as the variants it gave and the alleles its blocks held.

    `block.gts.size` is the shape of the array the core filled, which reaches
    numpy without a copy, so adding the sizes up shows a pass whose blocks
    carry no genotype. Nothing here reads a genotype.
    """
    blocks = open_the_file(path).iter_blocks(fields=())
    num_vars = 0
    alleles = 0
    for block in blocks:
        num_vars += int(block.gts.shape[0])
        alleles += int(block.gts.size)
    return num_vars, f"{num_vars} variants, {alleles} alleles"


def the_diversity(
    path: str, what: str, pops: dict[str, list[str]], draw: int | None
) -> tuple[int, str]:
    """One pass of `calc_pop_diversity` over the file at `path` with the
    statistics of `what`, as the variants of the pass and what it gave.

    What is printed is, for the first population, the value of each statistic
    that was asked for, and the variants in the draw for it, so that a pass
    which read nothing and one whose draw no variant reached both show
    themselves.
    """
    stats = STATS[what]
    diversity = popnei.calc_pop_diversity(
        open_the_file(path),
        pops,
        stats=stats,
        num_called_alleles=draw,
    )
    first = diversity.pops[0]
    said = [
        f"{diversity.pass_stats.num_vars} variants",
        (
            f"{int(diversity.num_vars.loc[first, 'in_draw'])} of them in the "
            f"draw for {first}"
        ),
    ]
    missing = [str(stat) for stat in stats if getattr(diversity, str(stat)) is None]
    if missing:
        raise ValueError(
            f"the pass was asked for {', '.join(str(stat) for stat in stats)} "
            f"and gave no value for {', '.join(missing)}"
        )
    if popnei.PopDiversityStat.NUM_ALLELES in stats:
        said.append(
            f"num_alleles {diversity.num_alleles.loc[first, 'total']} "
            f"in draw {diversity.num_alleles.loc[first, 'in_draw']:.6f}"
        )
    if popnei.PopDiversityStat.PRIVATE_ALLELES in stats:
        said.append(
            f"private_alleles {diversity.private_alleles.loc[first, 'total']} "
            f"in draw {diversity.private_alleles.loc[first, 'in_draw']:.6f}"
        )
    if popnei.PopDiversityStat.VARIABLE_VARS_RATIO in stats:
        said.append(
            f"variable_vars {diversity.variable_vars_ratio.loc[first, 'total']} "
            f"in draw {diversity.variable_vars_ratio.loc[first, 'in_draw']:.6f}"
        )
    if popnei.PopDiversityStat.FOLDED_SFS in stats:
        spectrum = diversity.folded_sfs[first]
        said.append(
            f"folded_sfs of {spectrum.shape[0]} bins summing to {spectrum.sum():.6f}"
        )
    if popnei.PopDiversityStat.FIS in stats:
        said.append(f"fis {diversity.fis[first]:.6f}")
    return diversity.pass_stats.num_vars, ", ".join(said)


def one_pass(
    path: str, what: str, pops: dict[str, list[str]], draw: int | None
) -> tuple[int, str]:
    """One whole pass over the file at `path`, from the call that opens it to
    the result, as the variants it gave and the line that says what it
    gave."""
    if what == "read":
        return the_read_alone(path)
    return the_diversity(path, what, pops, draw)


def the_draw(what: str, written: str) -> int | None:
    """The draw a run was given: the number for a pass that takes one, and
    `None` for a pass that takes the word `none`.

    A pass that needs a draw and was given `none`, and a pass that has no
    draw and was given a number, are both refused: each of them would time a
    pass other than the one that was named.
    """
    if what in WITH_A_DRAW:
        if written == NO_DRAW:
            raise ValueError(
                f"{what!r} is a pass of a draw and needs the called alleles to "
                f"draw, and {NO_DRAW!r} was written in their place"
            )
        return int(written)
    if written != NO_DRAW:
        raise ValueError(
            f"{what!r} is a pass of no draw and takes {NO_DRAW!r} in the place "
            f"of the called alleles, and {written!r} was written there"
        )
    return None


def main() -> int:
    arguments = sys.argv[1:]
    if len(arguments) != 5:
        print(__doc__)
        return 1
    path, what, pops_path, written_draw, runs = arguments
    if what not in WHATS:
        print(f"the pass to time is one of {WHATS}, and {what!r} was given")
        return 1
    try:
        draw = the_draw(what, written_draw)
    except ValueError as problem:
        print(problem)
        return 1
    runs = int(runs)
    pops = the_pops(pops_path)
    threads = os.environ.get("RAYON_NUM_THREADS", "one for each core")
    print(
        f"{path}, popnei {popnei.__version__}, {what}, "
        f"{len(pops)} populations of {pops_path}, "
        f"draw {draw if draw is not None else NO_DRAW}, {runs} runs, "
        f"RAYON_NUM_THREADS={threads}"
    )

    # The read that says how many variants the file holds, which every timed
    # run is checked against. It is a pass of its own and is not timed.
    expected, said = the_read_alone(path)
    print(f"the file holds {said}")

    started = time.perf_counter()
    num_vars, said = one_pass(path, what, pops, draw)
    took = time.perf_counter() - started
    if num_vars != expected:
        print(
            f"the first pass was to give the {expected} variants the file "
            f"holds and gave: {said}"
        )
        return 1
    print(f"the first pass, which is not timed: {took:.3f} s, {said}")

    times = []
    for run in range(1, runs + 1):
        started = time.perf_counter()
        num_vars, said = one_pass(path, what, pops, draw)
        took = time.perf_counter() - started
        if num_vars != expected:
            print(
                f"run {run}: the pass was to give the {expected} variants the "
                f"file holds and gave: {said}"
            )
            return 1
        times.append(took)
        print(f"run {run}: {took:.3f} s, {said}")
    print(
        f"{what}: best {min(times):.3f} s, "
        f"median {statistics.median(times):.3f} s, worst {max(times):.3f} s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
