"""Two things the numpy prototype exaggerates, and what an iterative solve costs.

The trace of C is measured again with the triangular matrix inverted in place
instead of solved against a dense identity built for it, which is what a Rust
version would do; and one product of the covariance with a vector is timed, to
see what a conjugate gradient solve would cost against one Cholesky.
"""
import sys, time, pathlib, numpy
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from scipy.linalg import solve_triangular
from scipy.linalg.lapack import dtrtri

def best_of(f, r=3):
    out = float("inf")
    for _ in range(r):
        s = time.perf_counter(); f(); out = min(out, time.perf_counter() - s)
    return out

rng = numpy.random.default_rng(0)
print("{:>6} {:>12} {:>12} {:>12} {:>12} {:>12}".format(
    "n", "cholesky", "dtrtri", "L\\diag(w)", "one matvec", "cho_solve c"))
for n in (1000, 2000, 4000):
    a = rng.standard_normal((n, n)); a = a @ a.T + n * numpy.eye(n)
    l = numpy.linalg.cholesky(a)
    w = rng.uniform(0.05, 0.25, n)
    v = rng.standard_normal(n)
    d = rng.standard_normal((n, 3))
    r = 3 if n <= 2000 else 2
    print("{:>6} {:>12.4f} {:>12.4f} {:>12.4f} {:>12.5f} {:>12.5f}".format(
        n,
        best_of(lambda: numpy.linalg.cholesky(a), r),
        best_of(lambda: dtrtri(l, lower=1), r),
        best_of(lambda: solve_triangular(l, numpy.diag(1/numpy.sqrt(w)), lower=True), r),
        best_of(lambda: a @ v, r),
        best_of(lambda: solve_triangular(l, d, lower=True), r)))
