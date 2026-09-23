"""How long the kinship of every pair of the individuals of a file takes in
Python, in pyNei and in popnei.

It is what the performance review of docs/reports/perf-kinship-2026-09-23.md
timed the two libraries with, against the 0.23 s of plink2 that "Speed" of
docs/specs/kinship.md states for 100000 variants of 1000 individuals with
every genotype called, and against the 0.81 s that section 2.1 of
docs/rust_core.md gives pyNei on the same sizes. `kinship.rs`, beside this
file, times the same kinship in the core crate, without Python.

One run is `calc_kinship` over a whole file, the reading included, which is
what a user waits for. Both libraries have a function of that name and
neither reads a variant when the file is opened, so both runs do the same
work: pyNei builds the dosages of every chunk and accumulates the samples x
samples matrix, and popnei standardizes each block and adds its product to
the same matrix.

What each library reads. A path that ends in `.vcf` or `.vcf.gz` is read by
pyNei with `vars_from_vcf` and by popnei with `open_vcf`; any other path is
read as the vars file of the library that is being timed, pyNei's with
`load_vars` and popnei's with `open_vars`. The two formats are different
files of the same variants, so a comparison of the two libraries over their
own vars files is a comparison of the whole job each of them does, and one
over the same VCF is a comparison that has the same parsing in both.

    uv run python time_kinship.py <path> <runs> [--popnei] [--num-pcs n]

Without `--popnei` it times pyNei and with it popnei. `--num-pcs` above 0
takes the principal components of the matrix inside the clock, which both
libraries have: popnei's `Kinship.principal_components` and pyNei's method of
the same name. With 0 there is no eigendecomposition and what is timed is the
one pass over the variants, which is what the target is stated on.

It prints the wall time of each run and the peak memory of the process, and
then the best, the median and the worst of the times. The peak is of the
process and not of one run, so it is the largest any run reached.

One run comes before the timed ones whose time is not taken, for the page
cache and for whatever the library imports on its first chunk.

The threads. popnei standardizes the rows of a block on rayon's global pool
and both libraries do their matrix work in a library of the system, which on
this machine is Accelerate, so one thread is asked for with
RAYON_NUM_THREADS=1 and VECLIB_MAXIMUM_THREADS=1 in the environment of the
command, which Accelerate reads when the process starts; numpy takes its own
threads from the same variable. pyNei's `calc_kinship` also takes a
`num_threads`, which is 1 here and is what its own `run_chunk_calcs` splits
the chunks over.

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

# How many principal components of the matrix are taken when the command line
# does not say. With 0 there is no eigendecomposition.
DEFAULT_NUM_PCS = 0

# The ends of a path that is read as a VCF. Any other path is read as the vars
# file of the library being timed.
A_VCF_ENDS_IN = (".vcf", ".vcf.gz")


def one_pass_of_pynei(path: str, num_pcs: int) -> tuple[int, int]:
    """The shape of the matrix of one whole kinship of the file at `path` with
    pyNei, with its first `num_pcs` principal components taken when that is
    above 0."""
    from pynei import calc_kinship, load_vars, vars_from_vcf

    if path.lower().endswith(A_VCF_ENDS_IN):
        variants = vars_from_vcf(path)
    else:
        variants = load_vars(path)
    kinship = calc_kinship(variants, num_threads=1)
    if num_pcs > 0:
        kinship.principal_components(num_pcs)
    return kinship.matrix.shape


def one_pass_of_popnei(path: str, num_pcs: int) -> tuple[int, int]:
    """The shape of the matrix of one whole kinship of the file at `path` with
    popnei, with its first `num_pcs` principal components taken when that is
    above 0."""
    import popnei

    if path.lower().endswith(A_VCF_ENDS_IN):
        variants = popnei.open_vcf(path)
    else:
        variants = popnei.open_vars(path)
    kinship = popnei.calc_kinship(variants)
    if num_pcs > 0:
        kinship.principal_components(num_pcs)
    return kinship.matrix.shape


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
    num_pcs = DEFAULT_NUM_PCS
    if "--num-pcs" in arguments:
        where = arguments.index("--num-pcs")
        num_pcs = int(arguments[where + 1])
        del arguments[where : where + 2]
    if len(arguments) != 2:
        print(__doc__)
        return 1
    path = arguments[0]
    runs = int(arguments[1])
    one_pass = one_pass_of_popnei if library == "popnei" else one_pass_of_pynei
    print(f"{path}, {library}, {runs} runs, {num_pcs} principal components")
    started = time.perf_counter()
    shape = one_pass(path, num_pcs)
    print(
        f"the first run, which is not timed: {time.perf_counter() - started:.3f} s, "
        f"a matrix of {shape}"
    )
    times = []
    for run in range(1, runs + 1):
        started = time.perf_counter()
        shape = one_pass(path, num_pcs)
        took = time.perf_counter() - started
        times.append(took)
        print(f"run {run}: {took:.3f} s, a matrix of {shape}, peak {peak_memory()}")
    print(
        f"best {min(times):.3f} s, median {statistics.median(times):.3f} s, "
        f"worst {max(times):.3f} s, peak {peak_memory()}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
