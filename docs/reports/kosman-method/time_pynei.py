import time, numpy, sys
from pynei.dists import _calc_kosman_dist_sums
def block(nv, ni, na, seed=42):
    rng=numpy.random.default_rng(seed)
    p=0.05+0.9*rng.random((nv,1,1))
    alt=rng.integers(1,na,size=(nv,ni,2))
    g=numpy.where(rng.random((nv,ni,2))<p,0,alt).astype(numpy.int8)
    miss=rng.random((nv,ni))<0.03
    g[miss]=-1
    return g
for nv,ni,na in ((5000,1000,2),(5000,1000,4),(500,10000,2)):
    g=block(nv,ni,na)
    best=1e9
    for _ in range(3):
        t=time.perf_counter(); s,n=_calc_kosman_dist_sums(g); best=min(best,time.perf_counter()-t)
    print(f"pyNei _calc_kosman_dist_sums, {nv} x {ni}, {na} alleles: {best:.4f} s")
print(numpy.__version__); numpy.show_config()
