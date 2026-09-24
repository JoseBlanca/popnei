"""Which summary of F_IS over the variants recovers the F the data were made with.

Genotypes are drawn at a known F: an individual is homozygous with
probability p^2 + F p q and q^2 + F p q and heterozygous with 2 p q (1 - F),
which is the standard way an inbreeding coefficient enters Hardy Weinberg.
The allele frequencies come from the neutral spectrum, where the number of
variants at frequency x goes as 1/x, so most variants are rare, as in a real
dataset.
"""
import numpy as np

def draw(rng, num_vars, num_inds, f, min_maf):
    # The neutral spectrum: p is drawn with density 1/p between min_maf and 1/2.
    u = rng.random(num_vars)
    p = min_maf * (0.5 / min_maf) ** u              # log uniform, the 1/x density
    q = 1.0 - p
    hom_ref = p * p + f * p * q
    het = 2 * p * q * (1 - f)
    probs = np.stack([hom_ref, het, 1.0 - hom_ref - het], axis=1)
    probs = np.clip(probs, 0.0, None)
    probs /= probs.sum(axis=1, keepdims=True)
    counts = np.array([rng.multinomial(num_inds, probs[v]) for v in range(num_vars)])
    return counts                                    # vars x (hom_ref, het, hom_alt)

def estimates(counts, num_inds, unbiased):
    n_ref, n_het, n_alt = counts[:, 0], counts[:, 1], counts[:, 2]
    ho = n_het / num_inds
    c = 2 * num_inds                                 # called alleles, none missing
    a_ref = 2 * n_ref + n_het
    a_alt = 2 * n_alt + n_het
    if unbiased:
        he = 1.0 - (a_ref * (a_ref - 1) + a_alt * (a_alt - 1)) / (c * (c - 1))
    else:
        he = 1.0 - ((a_ref / c) ** 2 + (a_alt / c) ** 2)
    varying = he > 0
    return (np.mean(1.0 - ho[varying] / he[varying]),   # mean of ratios
            1.0 - ho.mean() / he.mean())                # ratio of means

rng = np.random.default_rng(42)
NUM_VARS, REPS = 2000, 50
print(f"{NUM_VARS} variants, {REPS} datasets per row, allele frequencies from the "
      f"neutral spectrum down to the lowest major allele frequency shown\n")
for num_inds in (30, 100):
    for min_maf in (0.01, 0.05):
        for true_f in (0.0, 0.2, 0.5):
            rows = []
            for rep in range(REPS):
                counts = draw(rng, NUM_VARS, num_inds, true_f, min_maf)
                rows.append([*estimates(counts, num_inds, False),
                             *estimates(counts, num_inds, True)])
            a = np.array(rows)
            m = a.mean(axis=0); s = a.std(axis=0)
            print(f"  {num_inds:3d} individuals, rarest allele {min_maf:.2f}, "
                  f"true F {true_f:.1f}")
            for i, what in enumerate(("mean of ratios, plain He",
                                      "ratio of means, plain He",
                                      "mean of ratios, unbiased He",
                                      "ratio of means, unbiased He")):
                print(f"      {what:30} {m[i]:8.4f}  bias {m[i]-true_f:+8.4f}  "
                      f"sd {s[i]:.4f}")
            print()
