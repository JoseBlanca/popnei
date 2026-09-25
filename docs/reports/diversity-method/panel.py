"""The five quantities on the panel of docs/specs/stats.md, and the files the
outside programs read."""

import sys

import numpy as np

from diversity import read_vcf, calc

VCF = "/Users/jose/devel/popnei/.claude/worktrees/spec-diversity/tests/reference/stats/panel.vcf.gz"
POPFILE = "/Users/jose/devel/popnei/.claude/worktrees/spec-diversity/tests/reference/stats/panel_pops_bcftools.txt"

names, gts = read_vcf(VCF)
pops = {}
for line in open(POPFILE):
    ind, pop = line.split()
    pops.setdefault(pop, []).append(names.index(ind))
pops = {k: pops[k] for k in sorted(pops)}
print(f"{gts.shape[0]} variants, {len(names)} individuals, "
      f"{ {k: len(v) for k, v in pops.items()} }")

G = int(sys.argv[1]) if len(sys.argv) > 1 else 20
res = calc(gts, pops, ploidy=2, min_num_individuals=20, num_called_alleles=G)
print(f"\nnum_called_alleles = {G}, min_num_individuals = 20\n")
keys = ["num_vars", "num_vars_rare", "num_private_vars", "num_private_vars_rare",
        "alleles_total", "alleles_mean", "alleles_rare",
        "private_total", "private_mean", "private_rare",
        "poly_total", "poly_ratio", "poly_rare", "fis"]
print(f"{'':22}" + "".join(f"{p:>16}" for p in pops))
for k in keys:
    row = "".join(
        f"{res[p][k]:>16.10f}" if isinstance(res[p][k], float) else f"{res[p][k]:>16}"
        for p in pops)
    print(f"{k:22}{row}")
print(f"\nthe folded spectrum, projected to {G} called alleles")
print(f"{'rarer allele':22}" + "".join(f"{p:>18}" for p in pops))
for j in range(G // 2 + 1):
    print(f"{j:<22}" + "".join(f"{res[p]['sfs'][j]:>18.10f}" for p in pops))
for p in pops:
    print(f"  sum of {p}: {res[p]['sfs'].sum():.6f}, variants {res[p]['num_vars_rare']}")

# The per variant allele counts each outside program reads.
out = open("/private/tmp/claude-501/-Users-jose-devel-popnei/051eadf4-5b2a-4ec4-9fcd-cc5bfeeae213/scratchpad/div/panel_counts.tsv", "w")
out.write("var\tpop\tn0\tn1\tcalled\n")
for v in range(gts.shape[0]):
    for p, idx in pops.items():
        a = gts[v, idx, :]
        a = a[a >= 0]
        counts = np.bincount(a, minlength=2)
        out.write(f"{v}\t{p}\t{counts[0]}\t{counts[1]}\t{a.size}\n")
out.close()
print("\npanel_counts.tsv written")
