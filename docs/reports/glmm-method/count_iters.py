"""How pyNei's logistic mixed model fit spends its individuals x individuals work.

It counts the steps on the variance component of the kinship effect, tau, and
the linearizations inside each of them, on pyNei's reference panel and on
simulated panels of more individuals.
"""
import numpy, pandas, math, time
import pynei.gwas as g

real_fit = g._fit_pql_for_tau
counts = {"tau_steps": 0, "linearizations": 0}

def counting_fit(y, design, kinship, tau, eta, mu, tol, max_iter):
    counts["tau_steps"] += 1
    inner = {"n": 0}
    real_inv = numpy.linalg.inv
    def inv(a):
        if a.ndim == 2 and a.shape[0] == a.shape[1] and a.shape[0] > 50:
            inner["n"] += 1
        return real_inv(a)
    numpy.linalg.inv = inv
    try:
        out = real_fit(y, design, kinship, tau, eta, mu, tol, max_iter)
    finally:
        numpy.linalg.inv = real_inv
    counts["linearizations"] += inner["n"]
    counts.setdefault("per_step", []).append(inner["n"])
    return out

g._fit_pql_for_tau = counting_fit


def simulate(n_ind, n_var, seed=42, fst=0.1, family_size=4, num_pops=3):
    rng = numpy.random.default_rng(seed)
    num_families = n_ind // family_size
    family_pops = rng.integers(0, num_pops, size=num_families)
    p_anc = rng.uniform(0.1, 0.9, size=n_var)
    a = p_anc * (1 - fst) / fst
    b = (1 - p_anc) * (1 - fst) / fst
    p_pop = numpy.stack([rng.beta(a, b) for _ in range(num_pops)])
    alleles = numpy.empty((n_var, n_ind, 2), dtype=numpy.int8)
    pops = numpy.repeat(family_pops, family_size)
    for fam, pop in enumerate(family_pops):
        parents = (rng.uniform(size=(2, n_var, 2)) < p_pop[pop][None, :, None]).astype(numpy.int8)
        for child in range(family_size):
            idx = fam * family_size + child
            for parent in range(2):
                picked = rng.integers(0, 2, size=n_var)
                alleles[:, idx, parent] = parents[parent, numpy.arange(n_var), picked]
    dosages = alleles.sum(axis=2).T.astype(float)          # individuals x vars
    p = dosages.mean(axis=0) / 2
    poly = (p > 0.05) & (p < 0.95)
    z = (dosages[:, poly] - 2 * p[poly]) / numpy.sqrt(2 * p[poly] * (1 - p[poly]))
    kinship = z @ z.T / z.shape[1]
    effects = rng.standard_normal(z.shape[1]) * math.sqrt(0.5 / z.shape[1])
    genetic = z @ effects
    cov1 = rng.standard_normal(n_ind)
    cov2 = rng.integers(0, 2, size=n_ind).astype(float)
    liability = 0.5 * cov1 + 0.8 * cov2 + 0.7 * pops + genetic + rng.standard_normal(n_ind) * math.sqrt(0.5)
    y = (liability > numpy.quantile(liability, 0.6)).astype(float)
    design = numpy.column_stack([numpy.ones(n_ind), cov1, cov2])
    return y, design, kinship


for n_ind, n_var in ((200, 1200), (500, 2000), (1000, 2000)):
    y, design, kinship = simulate(n_ind, n_var)
    counts.update({"tau_steps": 0, "linearizations": 0, "per_step": []})
    start = time.perf_counter()
    null = g._GLMMNull(y, design, kinship)
    elapsed = time.perf_counter() - start
    print(f"{n_ind:>5} individuals: {elapsed:7.3f} s, tau steps {counts['tau_steps']:>3}, "
          f"linearizations {counts['linearizations']:>3} {counts['per_step']}, "
          f"tau={null.genetic_variance:.5f}")
