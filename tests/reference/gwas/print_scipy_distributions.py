r"""It prints what scipy gives for the two distributions of an association study.

Run from the root of the repository, with a scipy that the environment of popnei
does not have, since nothing of the library needs one:

    uv run --with "scipy==1.18.1" --with "numpy==2.5.3" python \
        tests/reference/gwas/print_scipy_distributions.py

It needs scipy 1.18.1 and numpy 2.5.3, the versions every literal of "The two
distributions" of `docs/specs/gwas.md` was taken with, and it refuses any other
version of either: scipy's `betainc` and `t.sf` are asserted to 1e-12 and 1e-10
of themselves, and numpy's generator is what the arguments are drawn with, so
another version of either would print other numbers.

What it prints is the four `const` blocks of `mod distributions` of
crates/popnei/src/gwas.rs, in the order they are there, ready to be pasted
over them. Every number is Python's `repr`, which gives the fewest digits that
read the same `f64` back; clippy refuses an `f64` literal of 17 digits, which
is what a fixed 17 digit format would print.

The three blocks are:

- `SCIPY_CHI2_SF_1DF`, 16 pairs of `x` and `chi2.sf(x, 1)`. They are already in
  gwas.rs, printed again here so that the whole test file can be got from
  scipy: what this script prints has to equal what is there, and a run that
  differs is a scipy or a numpy that is not the one above.
- `BETA_X` and `SCIPY_BETAINC`, 10 values of `x` and, for each of the four
  pairs `(a, b)` the module uses, `betainc(a, b, x)`.
- `T_VALUES` and `SCIPY_T_SF_TWO_SIDED`, 13 values of `t` and, for each of the
  three degrees of freedom, `2 * t.sf(|t|, df)`.

The arguments are pyNei's, from `test_distributions` of its `test/test_gwas.py`:
the incomplete beta over `x` drawn uniformly in (0, 1), the Student t over
normal draws of standard deviation 3 with 10, 20 and 40 added, and the chi
square over a sample of itself with 30, 50 and 100 added. pyNei asserts over
1000 draws of each, which are too many to hold as literals, so each sample of
1000 is sorted and read at 10 or 12 evenly spaced ranks, the largest of the
1000 among them. The largest matters: it is the only draw that takes the
incomplete beta of the pair `(98.5, 0.5)`, the one a Student t with 197 degrees
of freedom uses, through the symmetry `I_x(a, b) = 1 - I_{1-x}(b, a)`, which is
the second of the two branches of the function.
"""

import numpy
import scipy
from scipy import special, stats

SCIPY_VERSION = "1.18.1"
NUMPY_VERSION = "2.5.3"

# The four pairs of "How it is verified" of "The two distributions" of
# docs/specs/gwas.md, and the three degrees of freedom of the same section.
BETA_PAIRS = [(0.5, 0.5), (10.0, 0.5), (98.5, 0.5), (2.5, 7.0)]
DEGREES_OF_FREEDOM = [5.0, 17.0, 197.0]


def check_versions() -> None:
    """It stops when scipy or numpy is not the version the literals came from."""
    if scipy.__version__ != SCIPY_VERSION:
        raise SystemExit(
            f"scipy {SCIPY_VERSION} is what the literals were taken with, "
            f"and this is scipy {scipy.__version__}. "
            f'Run with: uv run --with "scipy=={SCIPY_VERSION}" '
            f'--with "numpy=={NUMPY_VERSION}" python {__file__}'
        )
    if numpy.__version__ != NUMPY_VERSION:
        raise SystemExit(
            f"numpy {NUMPY_VERSION} is what the arguments were drawn with, "
            f"and this is numpy {numpy.__version__}. "
            f'Run with: uv run --with "scipy=={SCIPY_VERSION}" '
            f'--with "numpy=={NUMPY_VERSION}" python {__file__}'
        )


def spaced_ranks(sample: numpy.ndarray, step: int) -> numpy.ndarray:
    """The sample sorted and read every `step` ranks, with its largest at the end."""
    ordered = numpy.sort(sample)
    taken = ordered[::step]
    if taken[-1] != ordered[-1]:
        taken = numpy.concatenate([taken, ordered[-1:]])
    return taken


def rust_number(value: float) -> str:
    """A float as a Rust literal: Python's `repr`, with a point in a whole number."""
    written = repr(float(value))
    if "." not in written and "e" not in written and "inf" not in written:
        written = f"{written}.0"
    return written


def print_pairs(name: str, arguments: numpy.ndarray, values: numpy.ndarray) -> None:
    """A `const` of pairs of an argument and what scipy gives for it."""
    print(f"    const {name}: [(f64, f64); {len(arguments)}] = [")
    for argument, value in zip(arguments, values, strict=True):
        print(f"        ({rust_number(argument)}, {rust_number(value)}),")
    print("    ];")


def print_row(name: str, values: numpy.ndarray) -> None:
    """A `const` of one row of floats."""
    print(f"    const {name}: [f64; {len(values)}] = [")
    for value in values:
        print(f"        {rust_number(value)},")
    print("    ];")


def print_rows(name: str, rows: list[numpy.ndarray], labels: list[str]) -> None:
    """A `const` of one row of floats per label, each row named in a comment."""
    print(f"    const {name}: [[f64; {len(rows[0])}]; {len(rows)}] = [")
    for label, row in zip(labels, rows, strict=True):
        print(f"        // {label}")
        print("        [")
        for value in row:
            print(f"            {rust_number(value)},")
        print("        ],")
    print("    ];")


def print_chi2() -> None:
    """The 16 pairs of `x` and `chi2.sf(x, 1)` that gwas.rs already holds."""
    sample = numpy.random.default_rng(0).chisquare(1, 1000)
    x = numpy.concatenate([spaced_ranks(sample, 90), [30.0, 50.0, 100.0]])
    print_pairs("SCIPY_CHI2_SF_1DF", x, stats.chi2.sf(x, 1))


def print_beta() -> None:
    """The 10 values of `x` and `betainc(a, b, x)` at each of the four pairs."""
    sample = numpy.random.default_rng(0).uniform(0, 1, 1000)
    x = spaced_ranks(sample, 111)
    print_row("BETA_X", x)
    print()
    rows = [special.betainc(a, b, x) for a, b in BETA_PAIRS]
    labels = [f"a = {a}, b = {b}" for a, b in BETA_PAIRS]
    print_rows("SCIPY_BETAINC", rows, labels)


def print_t() -> None:
    """The 13 values of `t` and `2 * t.sf(|t|, df)` at each of the three `df`."""
    sample = numpy.random.default_rng(0).standard_normal(1000) * 3
    t = numpy.concatenate([spaced_ranks(sample, 111), [10.0, 20.0, 40.0]])
    print_row("T_VALUES", t)
    print()
    rows = [2 * stats.t.sf(numpy.abs(t), df) for df in DEGREES_OF_FREEDOM]
    labels = [f"df = {df:g}" for df in DEGREES_OF_FREEDOM]
    print_rows("SCIPY_T_SF_TWO_SIDED", rows, labels)


def main() -> None:
    check_versions()
    print_chi2()
    print()
    print_beta()
    print()
    print_t()


if __name__ == "__main__":
    main()
