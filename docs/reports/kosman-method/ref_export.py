import numpy, pandas, sys, pynei
from pynei import Variants, calc_pairwise_kosman_dists
out=sys.argv[1]
def export(gts, samples, name):
    nv,ni,_=gts.shape
    cols={}
    for v in range(nv):
        a=gts[v]
        cols[f"v{v:04d}"]=[("NA" if (x<0 or y<0) else f"{x}/{y}") for x,y in a]
    pandas.DataFrame(cols,index=samples).to_csv(f"{out}/{name}.gts.csv")
# 1 the worked example
g=numpy.array([[(0,0),(0,1),(1,1)],[(0,1),(0,1),(1,2)],[(0,0),(0,-1),(2,2)],[(-1,-1),(1,1),(1,1)]],dtype=numpy.int8)
export(g,["s0","s1","s2"],"worked")
# 2 sim_missing
v=pynei.load_vars("/Users/jose/devel/pynei/test/gwas_reference/sim_missing.vars")
chunks=list(v.iter_vars_chunks()); gts=numpy.concatenate([c.gts.gt_values for c in chunks]); samples=list(v.samples)
print(gts.shape, gts.max(), (gts<0).mean())
export(gts,samples,"sim_missing")
d=calc_pairwise_kosman_dists(pynei.load_vars("/Users/jose/devel/pynei/test/gwas_reference/sim_missing.vars"))
numpy.savetxt(f"{out}/sim_missing.pynei.txt", d.dist_vector, fmt="%.17g")
# 3 multiallelic random, 40 individuals, 300 variants, up to 4 alleles, 5% missing genotypes
rng=numpy.random.default_rng(3); m=rng.integers(0,4,size=(300,40,2)).astype(numpy.int8); m[rng.random((300,40))<0.05]=-1
export(m,[f"i{k:02d}" for k in range(40)],"multi")
d=calc_pairwise_kosman_dists(Variants.from_gt_array(m,samples=[f"i{k:02d}" for k in range(40)]))
numpy.savetxt(f"{out}/multi.pynei.txt", d.dist_vector, fmt="%.17g")
