import numpy
from pynei.variants import Variants
from pynei.var_filters import filter_by_ld_and_maf
# 6 variants of 5 individuals; the first has no variance
gts=[[(0,0)]*5,
     [(0,0),(0,1),(1,1),(0,0),(0,1)],
     [(0,1),(1,1),(0,0),(0,1),(1,1)],
     [(1,1),(0,0),(0,1),(1,1),(0,0)],
     [(0,0),(0,1),(0,1),(1,1),(0,0)],
     [(0,1),(0,0),(1,1),(0,1),(0,0)]]
for first_has_variance in (False, True):
    g=[r for r in gts] if not first_has_variance else gts[1:]+[gts[0]]
    v=Variants.from_gt_array(g, samples=[f"s{i}" for i in range(5)])
    f=filter_by_ld_and_maf(v, min_allowed_r2=0.9, max_allowed_maf=1.0)
    kept=sum(c.num_vars for c in f.iter_vars_chunks())
    print(f"first variant has variance={first_has_variance}: pyNei keeps {kept} of {len(g)}")
