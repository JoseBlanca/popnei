"""How long pyNei's Jost's D between populations takes over a whole vars
file of its own.

It is the program "Speed" of docs/specs/dists.md compares popnei's
`calc_pop_dists` with, beside `time_pop_dists.py`, which is the same clock on
popnei over the same dataset and the same populations file.

What the two calls do is not the same amount of work. popnei calculates
seven measures of how far apart two populations are, each value with the
standard error of a block jackknife beside it, where pyNei calculates one of
the seven, Jost's D, and no standard error, from counts the two share. So
the comparison is of one pass over the variants against one pass over the
variants, and not of one number against one number.

One run is one call of `calc_jost_dest_pop_dists` over the `Variants` that
`load_vars` gives, with `num_threads` as it is asked for. `load_vars` is
outside the timing, as `open_vars` is outside popnei's: it memory maps the
file and reads its schema and its metadata, and the chunks themselves are
read inside the call. A fresh `Variants` is built for every run, so that no
chunk of it is read twice.

The file is written from a VCF by pyNei's own `write_vars`, because pyNei
does not read a vars file of popnei: the two formats are both arrow IPC and
their columns differ. That is one line:

    uv run python -c "from pynei import vars_from_vcf; \\
        from pynei.io_vars import write_vars; \\
        write_vars(vars_from_vcf('big.vcf'), 'big.pynei.vars')"

    uv run python time_pop_dists_pynei.py <path to a vars file> <populations file> <runs> <threads>

The populations come from a file of one line for each individual: the name
of the individual, a tab, and the name of its population, which is what
`make_pops.py` of popnei's benches wrote as `pops3.tsv` and `pops20.tsv`.
pyNei's `Pops` is a dict of population name to the names of its individuals,
the same shape popnei's Python takes, so both scripts read the same file. It
is read here in the order in which the file first names each population, and
pyNei then sorts the names itself, so the first pair of its result is of the
two names that sort first.

One run before the timed ones is not timed, for the page cache and for
whatever pyNei imports on its first chunk. pyNei is a development dependency
of popnei at the commit `pyproject.toml` names, so `uv run` from the popnei
worktree is what has it; `--no-project` would not.

What the call does that popnei's does not: for each pair of the populations
it counts the alleles and the called genotypes of every population of the
chunk again, `_count_alleles_per_var(chunk, pops=pop_idxs)` and
`_calc_obs_het_per_var(chunk, pops=pop_idxs)` inside `_calc_pairwise_dest`,
which `_DestPopHsHtCalculator` calls once for each pair. So the work over one
chunk is the counting of all the individuals of the file done once for each
pair, 190 times for 20 populations and 3 times for 3, where popnei counts
each population once for each variant and then works on the counts pair by
pair. It also builds three pandas frames of population by population for
every chunk, and reads each chunk through `VariantsFile.read_chunk`, which
builds a pandas frame of the chromosome, the position, the id and the
quality of that chunk whatever the caller will read.

`min_num_samples` is left at its default of 20 called genotypes, which is
popnei's `min_num_individuals`, and `alleles` at `None`, which counts the
alleles each chunk has.
"""

import statistics
import sys
import time

from pynei import calc_jost_dest_pop_dists, load_vars


def the_pops(path: str) -> dict[str, list[str]]:
    """The populations of the file at `path`, as the `Pops` dict of
    population name to the names of its individuals that
    `calc_jost_dest_pop_dists` takes.

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


def one_call(
    path: str, pops: dict[str, list[str]], num_threads: int
) -> tuple[str, float]:
    """One call of pyNei's `calc_jost_dest_pop_dists` over the vars file at
    `path` for the pairs of `pops`, as what the result holds and the seconds
    the call took.

    `load_vars` runs before the clock starts: it reads the schema and the
    metadata of the file and no chunk of genotypes.
    """
    variants = load_vars(path)
    started = time.perf_counter()
    dists = calc_jost_dest_pop_dists(variants, pops=pops, num_threads=num_threads)
    took = time.perf_counter() - started
    said = (
        f"{dists.dist_vector.shape[0]} pairs, Jost's D of {dists.names[0]} "
        f"and {dists.names[1]} {dists.dist_vector[0]:.6f}"
    )
    return said, took


def main() -> int:
    if len(sys.argv) != 5:
        print(__doc__)
        return 1
    path = sys.argv[1]
    pops_path = sys.argv[2]
    runs = int(sys.argv[3])
    num_threads = int(sys.argv[4])
    pops = the_pops(pops_path)
    print(
        f"{path}, pyNei, {len(pops)} populations of {pops_path}, {runs} runs, "
        f"num_threads={num_threads}"
    )
    said, took = one_call(path, pops, num_threads)
    print(f"the first run, which is not timed: {took:.3f} s, {said}")
    times = []
    for run in range(1, runs + 1):
        said, took = one_call(path, pops, num_threads)
        times.append(took)
        print(f"run {run}: {took:.3f} s, {said}")
    print(
        f"num_threads={num_threads}: best {min(times):.3f} s, "
        f"median {statistics.median(times):.3f} s, worst {max(times):.3f} s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
