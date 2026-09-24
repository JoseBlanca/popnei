"""Which curve, and which n, recover a 4Nr that is known.

A Wright-Fisher population is run to drift-recombination equilibrium, the
process the curve of "The curve that is fitted" of docs/specs/ld.md is the
expectation of, so its 4Nr is known and the fit can be asked for it back.
The population has two chromosomes that assort independently, so a pair of
variants on different ones carries no linkage disequilibrium and its r2 is
what sampling alone puts there, which is what the second factor of the
curve corrects for. It needs no file and no other program, about 40 s.
"""
import numpy

N = 200                 # diploid individuals in the population
NUM_HAPS = 2 * N
NUM_SITES = 1500
SPACING = 200           # bp between neighbouring sites
R_PER_BP = 2.5e-7
MU = 1e-4               # per site per generation, both ways
GENERATIONS = 1000      # 5 times the population, for drift to settle
C_TRUE = 4 * N * R_PER_BP
MAX_MAF = 0.95
SAMPLE = 100            # individuals the curve is fitted on
SEEDS = (7, 11, 23)
POSS = numpy.arange(NUM_SITES) * SPACING


def population(seed):
    """Two chromosomes of NUM_HAPS haplotypes, at equilibrium."""
    rng = numpy.random.default_rng(seed)
    r_adj = R_PER_BP * SPACING
    haps = [rng.integers(0, 2, size=(NUM_HAPS, NUM_SITES)).astype(numpy.int8)
            for _ in range(2)]
    for _ in range(GENERATIONS):
        parents = rng.integers(0, N, size=NUM_HAPS)
        new = []
        for chrom in haps:
            a, b = chrom[2 * parents], chrom[2 * parents + 1]
            switches = rng.random((NUM_HAPS, NUM_SITES)) < r_adj
            switches[:, 0] = rng.random(NUM_HAPS) < 0.5
            child = numpy.where(numpy.cumsum(switches, axis=1) % 2 == 0, a, b)
            mutate = rng.random((NUM_HAPS, NUM_SITES)) < MU
            new.append(numpy.where(mutate, 1 - child, child).astype(numpy.int8))
        haps = new
    return haps


def dosages(chrom, individuals):
    """One number per individual per site: the alleles that are not the major one."""
    d = (chrom[2 * individuals] + chrom[2 * individuals + 1]).T.astype(numpy.float64)
    freq_of_1 = d.sum(axis=1) / (2 * d.shape[1])
    major_is_1 = freq_of_1 > 0.5
    d[major_is_1] = 2 - d[major_is_1]
    maf = numpy.where(major_is_1, freq_of_1, 1 - freq_of_1)
    return d[maf <= MAX_MAF], maf <= MAX_MAF


def correlations(x, y):
    xc = x - x.mean(axis=1, keepdims=True)
    yc = y - y.mean(axis=1, keepdims=True)
    sx, sy = numpy.sqrt((xc * xc).sum(axis=1)), numpy.sqrt((yc * yc).sum(axis=1))
    with numpy.errstate(invalid="ignore", divide="ignore"):
        c = (xc @ yc.T) / numpy.outer(sx, sy)
    return c * c


def hw(d, C, n):
    p = C * d
    return ((10 + p) / ((2 + p) * (11 + p))) * (
        1 + ((3 + p) * (12 + 12 * p + p * p)) / (n * (2 + p) * (11 + p)))


def plain(d, C, n):
    p = C * d
    return (10 + p) / ((2 + p) * (11 + p))


def sved(d, C, n):
    return 1 / (1 + C * d)


def fit(dists, num_pairs, sum_r2, f, n):
    """The grid and the golden section docs/specs/ld.md specifies."""
    def ss(log_c):
        y = f(dists, 10.0**log_c, n)
        return float((num_pairs * y * y - 2 * sum_r2 * y).sum())
    grid = numpy.linspace(-12, 2, 141)
    vals = [ss(l) for l in grid]
    i = int(numpy.argmin(vals))
    if i in (0, len(grid) - 1):
        return float("nan")
    lo, hi = grid[i - 1], grid[i + 1]
    best_l, best_v = grid[i], vals[i]
    g = (numpy.sqrt(5) - 1) / 2
    x1, x2 = hi - g * (hi - lo), lo + g * (hi - lo)
    f1, f2 = ss(x1), ss(x2)
    for l, v in ((x1, f1), (x2, f2)):
        if v < best_v:
            best_l, best_v = l, v
    while hi - lo > 1e-9:
        if f1 < f2:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - g * (hi - lo); f1 = ss(x1)
            if f1 < best_v: best_l, best_v = x1, f1
        else:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + g * (hi - lo); f2 = ss(x2)
            if f2 < best_v: best_l, best_v = x2, f2
    return 10.0**best_l


def pairs_by_distance(chrom, individuals):
    d, keep = dosages(chrom, individuals)
    poss = POSS[keep]
    m = correlations(d, d)
    acc = {}
    for a in range(len(poss)):
        row, dist = m[a, a + 1:], poss[a + 1:] - poss[a]
        ok = ~numpy.isnan(row)
        for dd, vv in zip(dist[ok], row[ok]):
            n_d, s_d = acc.get(int(dd), (0, 0.0))
            acc[int(dd)] = (n_d + 1, s_d + float(vv))
    ds = numpy.array(sorted(acc), dtype=numpy.float64)
    return (ds, numpy.array([acc[int(x)][0] for x in ds], dtype=numpy.float64),
            numpy.array([acc[int(x)][1] for x in ds]))


print(f"{N} individuals, two chromosomes of {NUM_SITES} sites {SPACING} bp apart,")
print(f"recombination {R_PER_BP:g} per bp, mutation {MU:g} per site, "
      f"{GENERATIONS} generations,")
print(f"so 4Nr per bp is {C_TRUE:g} and the curve is fitted on {SAMPLE} individuals.\n")

for seed in SEEDS:
    haps = population(seed)
    ind = numpy.sort(numpy.random.default_rng(1).choice(N, SAMPLE, replace=False))
    ds, num_pairs, sum_r2 = pairs_by_distance(haps[0], ind)
    print(f"seed {seed}: {int(num_pairs.sum())} pairs at {len(ds)} distances, "
          f"mean r2 at {ds[0]:.0f} bp is {sum_r2[0] / num_pairs[0]:.4f}")
    for name, f, n in (("n = individuals", hw, SAMPLE),
                       ("n = individuals x ploidy", hw, 2 * SAMPLE),
                       ("no sample term", plain, 1),
                       ("Sved", sved, 1)):
        C = fit(ds, num_pairs, sum_r2, f, n)
        print(f"    {name:<26} 4Nr {C:.6g}  {C / C_TRUE - 1:+7.1%} of the true one"
              f"   the curve at {ds[0]:.0f} bp {float(f(ds[:1], C, n)[0]):.4f}")

    # a pair on two chromosomes carries no real r2, so its mean is what the
    # second factor of the curve is for: a + b / (individuals - 1) with b of 1
    # when the sampling goes by individuals and 0.5 when it goes by gametes.
    sizes = numpy.array([25, 50, 100, 200], dtype=numpy.float64)
    means = []
    for num in sizes:
        ii = numpy.sort(numpy.random.default_rng(2).choice(N, int(num), replace=False))
        d1, _ = dosages(haps[0], ii)
        d2, _ = dosages(haps[1], ii)
        v = correlations(d1, d2)
        means.append(float(v[~numpy.isnan(v)].mean()))
    means = numpy.array(means)
    design = numpy.column_stack([numpy.ones_like(sizes), 1 / (sizes - 1)])
    a, b = numpy.linalg.lstsq(design, means, rcond=None)[0]
    print(f"    unlinked pairs: b = {b:.3f} against 1 for the individuals and 0.5 "
          f"for the gametes,")
    print(f"      the mean r2 of an unlinked pair at "
          + ", ".join(f"{int(s)} individuals {m:.5f}" for s, m in zip(sizes, means)))
    print(f"      what is left over, {a:.5f}, is the population's own, so its "
          f"effective size is about {1 / (2 * a):.0f}\n")
