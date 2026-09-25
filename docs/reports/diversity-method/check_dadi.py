"""The folded spectrum of each population of the panel, projected by dadi."""

import csv
from collections import defaultdict

import dadi

G = 20
rows = defaultdict(dict)
with open("panel_counts.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rows[r["pop"]][r["var"]] = (int(r["n0"]), int(r["n1"]))

for pop in sorted(rows):
    dd = {
        v: {
            "calls": {pop: counts},
            "segregating": ("A", "T"),
            "outgroup_allele": "A",
        }
        for v, counts in rows[pop].items()
    }
    fs = dadi.Spectrum.from_data_dict(dd, [pop], projections=[G], polarized=False)
    print(pop, "sum", round(float(fs.data.sum()), 6))
    for j, value in enumerate(fs.data[: G // 2 + 1]):
        print(f"   {j:3d} {value:.10f}")
