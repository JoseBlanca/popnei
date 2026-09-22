"""How long the principal component analysis of the variants of a file takes
in Python, in pyNei and in popnei.

It is what task 4.1 of docs/plans/pca.md timed the two libraries with, against
the 0.3 s of "Speed" of docs/specs/pca.md.
docs/reports/pca-measurement.md has the numbers, the files they were taken on
and what each library reads. `pca_vars.rs`, beside this file, times the same
analysis in the core crate, without Python.

One run is `do_pca_from_variants` over a whole file, the reading included,
which is what a user waits for. Both libraries have a function of that name
and neither reads a variant when the file is opened, so both runs do the same
work: pyNei builds the matrix of the dosages of every variant in memory and
then calls numpy, and popnei adds the product of each block to G and reads the
file a second time for the weights when `num_prin_comps` is above 0.

What each library reads. A path that ends in `.vcf` or `.vcf.gz` is read by
pyNei with `vars_from_vcf` and by popnei with `open_vcf`; any other path is
read as the vars file of the library that is being timed, pyNei's with
`load_vars` and popnei's with `open_vars`. The two formats are different files
of the same variants, so a comparison of the two libraries over their own vars
files is a comparison of the whole job each of them does, and one over the
same VCF is a comparison that has the same parsing in both.

    uv run python time_pca.py <path> <runs> [--popnei] [--num-prin-comps n]

Without `--popnei` it times pyNei and with it popnei. `--num-prin-comps` is
popnei's alone: pyNei gives the weights of every variant in every component
and has no such argument, so the number it is timed with is popnei's with 0,
which gives no weights either way for what pyNei does not, and the run with 10
says what the second pass costs.

It prints the wall time of each run and the peak memory of the process, which
is what the 6.7 GB of pyNei in "Speed" of docs/specs/pca.md is, and then the
best, the median and the worst of the times. The peak is of the process and
not of one run, so it is the largest any run reached; a run of pyNei holds the
matrix of every dosage, 0.8 GB of f64 for 100000 variants of 1000 individuals
before numpy copies it, and one of popnei holds one block.

One run comes before the timed ones whose time is not taken, for the page
cache and for whatever the library imports on its first block.

The threads. popnei standardizes the rows of a block on rayon's global pool
and both libraries do their matrix work in a library of the system, which on
this machine is Accelerate, so one thread is asked for with
RAYON_NUM_THREADS=1 and VECLIB_MAXIMUM_THREADS=1 in the environment of the
command, which Accelerate reads when the process starts; numpy takes its own
threads from the same variable.

pyNei is a development dependency of popnei at the commit `pyproject.toml`
names, so `uv run` is what has it; `--no-project` would not. popnei is the
module that `uv run maturin develop --release` leaves in the environment, and
a `maturin develop` with no `--release` leaves the unoptimized build there,
which is not the one to time.
"""

import resource
import statistics
import sys
import time

# How many components popnei gives the weights for when the command line does
# not say. With 0 there is no second pass over the variants.
DEFAULT_NUM_PRIN_COMPS = 0

# The ends of a path that is read as a VCF. Any other path is read as the vars
# file of the library being timed.
A_VCF_ENDS_IN = (".vcf", ".vcf.gz")


def one_pass_of_pynei(path: str, num_prin_comps: int) -> tuple[int, int]:
    """The shape of the projections of one whole analysis of the file at
    `path` with pyNei. `num_prin_comps` is not read: pyNei has no such
    argument."""
    del num_prin_comps
    from pynei import do_pca_from_variants, load_vars, vars_from_vcf

    if path.lower().endswith(A_VCF_ENDS_IN):
        variants = vars_from_vcf(path)
    else:
        variants = load_vars(path)
    return do_pca_from_variants(variants).projections.shape


def one_pass_of_popnei(path: str, num_prin_comps: int) -> tuple[int, int]:
    """The shape of the projections of one whole analysis of the file at
    `path` with popnei, with the weights of `num_prin_comps` components."""
    import popnei

    if path.lower().endswith(A_VCF_ENDS_IN):
        variants = popnei.open_vcf(path)
    else:
        variants = popnei.open_vars(path)
    result = popnei.do_pca_from_variants(variants, num_prin_comps=num_prin_comps)
    return result.projections.shape


def peak_memory() -> str:
    """The largest resident memory the process has held, which on macOS
    `getrusage` gives in bytes."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return f"{peak / 1e9:.2f} GB"


def main() -> int:
    arguments = sys.argv[1:]
    library = "pyNei"
    if "--popnei" in arguments:
        arguments.remove("--popnei")
        library = "popnei"
    num_prin_comps = DEFAULT_NUM_PRIN_COMPS
    if "--num-prin-comps" in arguments:
        where = arguments.index("--num-prin-comps")
        num_prin_comps = int(arguments[where + 1])
        del arguments[where : where + 2]
    if len(arguments) != 2:
        print(__doc__)
        return 1
    path = arguments[0]
    runs = int(arguments[1])
    one_pass = one_pass_of_popnei if library == "popnei" else one_pass_of_pynei
    print(
        f"{path}, {library}, {runs} runs, weights for "
        f"{num_prin_comps if library == 'popnei' else 'every'} component"
    )
    started = time.perf_counter()
    shape = one_pass(path, num_prin_comps)
    print(
        f"the first run, which is not timed: {time.perf_counter() - started:.3f} s, "
        f"projections of {shape}"
    )
    times = []
    for run in range(1, runs + 1):
        started = time.perf_counter()
        shape = one_pass(path, num_prin_comps)
        took = time.perf_counter() - started
        times.append(took)
        print(f"run {run}: {took:.3f} s, projections of {shape}, peak {peak_memory()}")
    print(
        f"best {min(times):.3f} s, median {statistics.median(times):.3f} s, "
        f"worst {max(times):.3f} s, peak {peak_memory()}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
