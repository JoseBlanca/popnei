"""The worked example of docs/specs/filters.md, with the two populations."""

import numpy as np

from diversity import calc

M = -1
# i1 i2 | i3 i4 i5, the six variants of the worked example.
GTS = np.array(
    [
        [[0, 0], [0, 1], [0, 0], [0, 0], [0, M]],
        [[0, 0], [0, 1], [0, 0], [M, M], [0, M]],
        [[0, 1], [2, 3], [0, 1], [2, 3], [M, M]],
        [[M, M], [M, M], [M, M], [M, M], [M, M]],
        [[0, 0], [0, 0], [0, 0], [0, 0], [1, 1]],
        [[0, M], [M, M], [M, M], [M, M], [M, M]],
    ],
    dtype=np.int16,
)
POPS = {"pop1": [0, 1], "pop2": [2, 3, 4]}

for g in (None, 2, 4):
    res = calc(GTS, POPS, ploidy=2, min_num_individuals=1, num_called_alleles=g)
    print(f"\n=== num_called_alleles = {g} ===")
    for name, r in res.items():
        print(f"  {name}")
        for k, v in r.items():
            if k == "sfs":
                print(f"    {k:22} {np.round(v, 6)}")
            elif isinstance(v, float):
                print(f"    {k:22} {v:.6f}")
            else:
                print(f"    {k:22} {v}")
