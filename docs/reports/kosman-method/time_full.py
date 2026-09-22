import time, numpy
from pynei import Variants, calc_pairwise_kosman_dists
rng=numpy.random.default_rng(42); nv,ni=100000,1000
p=(0.05+0.9*rng.random((nv,1,1))).astype(numpy.float32)
g=(rng.random((nv,ni,2),dtype=numpy.float32)>=p).astype(numpy.int8); g[rng.random((nv,ni))<0.03]=-1
for th in (1,6):
    best=1e9
    for _ in range(2):
        v=Variants.from_gt_array(g,samples=[str(i) for i in range(ni)])
        t=time.perf_counter(); d=calc_pairwise_kosman_dists(v,num_threads=th); best=min(best,time.perf_counter()-t)
    print(f"pyNei calc_pairwise_kosman_dists {nv} x {ni}, num_threads={th}: {best:.3f} s")
