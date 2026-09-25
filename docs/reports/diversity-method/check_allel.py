"""F_IS of each population of the panel, from the heterozygosities scikit-allel
gives per variant, and the genotype table the R checks read."""

import allel
import numpy as np

from diversity import read_vcf

VCF = "/Users/jose/devel/popnei/.claude/worktrees/spec-diversity/tests/reference/stats/panel.vcf.gz"
POPFILE = "/Users/jose/devel/popnei/.claude/worktrees/spec-diversity/tests/reference/stats/panel_pops_bcftools.txt"

names, gts = read_vcf(VCF)
pops = {}
for line in open(POPFILE):
    ind, pop = line.split()
    pops.setdefault(pop, []).append(names.index(ind))
pops = {k: pops[k] for k in sorted(pops)}

print("scikit-allel", allel.__version__)
for pop, idx in pops.items():
    g = allel.GenotypeArray(gts[:, idx, :])
    ho = allel.heterozygosity_observed(g)
    af = g.count_alleles().to_frequencies()
    he = allel.heterozygosity_expected(af, ploidy=2)
    f_per_var = allel.inbreeding_coefficient(g)
    ok = ~np.isnan(ho) & ~np.isnan(he)
    print(f"  {pop}: ratio of the means {1 - ho[ok].mean() / he[ok].mean():.10f}   "
          f"mean of the per variant F {np.nanmean(f_per_var):.10f}")

# The genotypes as text, for df2genind in R.
with open("panel_genotypes.tsv", "w") as out:
    out.write("ind\tpop\t" + "\t".join(f"v{v}" for v in range(gts.shape[0])) + "\n")
    where = {i: p for p, idx in pops.items() for i in idx}
    for i, name in enumerate(names):
        cells = []
        for v in range(gts.shape[0]):
            a, b = gts[v, i]
            cells.append("NA" if a < 0 or b < 0 else f"{a}/{b}")
        out.write(f"{name}\t{where[i]}\t" + "\t".join(cells) + "\n")
print("panel_genotypes.tsv written")
