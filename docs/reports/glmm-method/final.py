"""The two measurements the decision turns on.

First, what each fit costs beside the per variant test it feeds, at the
individual counts of docs/objectives.md. Second, whether the Cholesky of the
cheaper fit refuses a kinship that missing genotypes made not positive
semidefinite, where pyNei's LU carries on.
"""
import sys, time, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy
from fits import fit_a, fit_c, simulate

def best_of(f, r):
    out = float("inf")
    for _ in range(r):
        s = time.perf_counter(); res = f(); out = min(out, time.perf_counter() - s)
    return out, res

print("The null fit, and the score test of 100000 variants it feeds")
print("{:>6} {:>10} {:>10} {:>6} {:>8} {:>16}".format(
    "n", "A pyNei", "C", "A/C", "test", "linearizations"))
for n_ind in (500, 1000, 2000, 4000):
    y, design, kinship, _ = simulate(n_ind, 2000)
    r = 3 if n_ind <= 1000 else 2
    ta, ra = best_of(lambda: fit_a(y, design, kinship.copy()), r)
    tc, rc = best_of(lambda: fit_c(y, design, kinship.copy()), r)
    assert abs(ra["tau"] - rc["tau"]) / ra["tau"] < 1e-9
    # one block of 5000 variants against the projection, scaled to 100000
    block = numpy.asarray(numpy.random.default_rng(0).standard_normal((5000, n_ind)))
    s = time.perf_counter(); block @ rc["projection"]; per_block = time.perf_counter() - s
    print("{:>6} {:>10.3f} {:>10.3f} {:>6.2f} {:>8.2f} {:>16}".format(
        n_ind, ta, tc, ta / tc, per_block * 20, rc["work"]["choleskys"]))

print("\nA kinship with missing genotypes: does the Cholesky of the cheaper fit refuse?")
for rate in (0.0, 0.03, 0.10, 0.25, 0.50):
    y, design, kinship, dosages = simulate(400, 2000, missing=0.0)
    rng = numpy.random.default_rng(7)
    d = dosages.copy()
    mask = rng.uniform(size=d.shape) < rate
    d[mask] = numpy.nan
    with numpy.errstate(invalid="ignore"):
        means = numpy.nanmean(d, axis=1)
    freqs = means / 2
    poly = (freqs > 0.0) & (freqs < 1.0) & numpy.isfinite(means)
    d = numpy.where(numpy.isnan(d), means[:, None], d)[poly]
    m = means[poly]
    z = (d - m[:, None]) / numpy.sqrt(2 * freqs[poly] * (1 - freqs[poly]))[:, None]
    called = (~mask)[poly].astype(float)
    per_pair = called.T @ called if rate else numpy.full((400, 400), float(poly.sum()))
    with numpy.errstate(invalid="ignore", divide="ignore"):
        k = (z.T @ z) / per_pair
    ev = numpy.linalg.eigvalsh(k)
    line = f"missing {rate:>5.0%}: kinship eigenvalues {ev.min():+.4f} to {ev.max():.3f}"
    for name, fit in (("A", fit_a), ("C", fit_c)):
        try:
            res = fit(y, design, k.copy())
            line += f"  {name} tau={res['tau']:.4f}"
        except Exception as err:
            line += f"  {name} {type(err).__name__}"
    print(line)
