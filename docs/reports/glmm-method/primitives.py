import time, numpy
def best_of(f, r=3):
    out = float("inf")
    for _ in range(r):
        s = time.perf_counter(); f(); out = min(out, time.perf_counter() - s)
    return out
rng = numpy.random.default_rng(0)
print("{:>6} {:>10} {:>10} {:>10} {:>10}".format("n", "inv", "cholesky", "eigvalsh", "eigh"))
for n in (500, 1000, 2000, 4000):
    a = rng.standard_normal((n, n)); a = a @ a.T + n * numpy.eye(n)
    r = 2 if n >= 2000 else 3
    print("{:>6} {:>10.4f} {:>10.4f} {:>10.4f} {:>10.4f}".format(
        n, best_of(lambda: numpy.linalg.inv(a), r),
        best_of(lambda: numpy.linalg.cholesky(a), r),
        best_of(lambda: numpy.linalg.eigvalsh(a), r),
        best_of(lambda: numpy.linalg.eigh(a), r)))
