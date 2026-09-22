"""How long pyNei's Kosman distances take over a whole vars file of its own.

It is what task 3.1 of docs/plans/dists-kosman.md timed pyNei with, beside
`time_kosman_dists.py`, which times popnei over the same dataset.
`docs/reports/dists-kosman-measurement.md` has the numbers and the load
averages they were taken at.

One run is one call of `calc_pairwise_kosman_dists` over the `Variants` that
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

    uv run python time_kosman_pynei.py <path to a vars file> <runs> <threads>

One run before the timed ones is not timed, for the page cache and for
whatever pyNei imports on its first chunk. pyNei is a development dependency
of popnei at the commit `pyproject.toml` names, so `uv run` is what has it;
`--no-project` would not.

What the call does that popnei's does not: pyNei reads a chunk of its vars
file through `VariantsFile.read_chunk`, which builds a pandas frame of the
chromosome, the position, the id and the quality of every chunk whatever the
caller will read, and computes the distances as the products of matrices of
0 and 1 in float32, which numpy makes with the BLAS it is linked against.
"""

import statistics
import sys
import time

from pynei import calc_pairwise_kosman_dists, load_vars


def one_call(path: str, num_threads: int) -> tuple[int, float]:
    """One call of pyNei's `calc_pairwise_kosman_dists` over the vars file at
    `path`, as the distances of its pairs and the seconds the call took.

    `load_vars` runs before the clock starts: it reads the schema and the
    metadata of the file and no chunk of genotypes.
    """
    variants = load_vars(path)
    started = time.perf_counter()
    distances = calc_pairwise_kosman_dists(variants, num_threads=num_threads)
    took = time.perf_counter() - started
    return distances.dist_vector.shape[0], took


def main() -> int:
    if len(sys.argv) != 4:
        print(__doc__)
        return 1
    path = sys.argv[1]
    runs = int(sys.argv[2])
    num_threads = int(sys.argv[3])
    print(f"{path}, pyNei, {runs} runs, num_threads={num_threads}")
    pairs, took = one_call(path, num_threads)
    print(f"the first run, which is not timed: {took:.3f} s, {pairs} pairs")
    times = []
    for run in range(1, runs + 1):
        pairs, took = one_call(path, num_threads)
        times.append(took)
        print(f"run {run}: {took:.3f} s, {pairs} pairs")
    print(
        f"num_threads={num_threads}: best {min(times):.3f} s, "
        f"median {statistics.median(times):.3f} s, worst {max(times):.3f} s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
