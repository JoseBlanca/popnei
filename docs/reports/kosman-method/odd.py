import numpy, warnings, pynei
from pynei import Variants, calc_pairwise_kosman_dists
from pynei.dists import _KosmanDistCalculator
M=-1
def run(gts, **kw):
    gts=numpy.array(gts)
    v=Variants.from_gt_array(gts, samples=[f"s{i}" for i in range(gts.shape[1])])
    if 'chunk' in kw: v.desired_num_vars_per_chunk=kw.pop('chunk')
    return calc_pairwise_kosman_dists(v, **kw)
# worked example: 4 variants, 3 individuals
g=[[(0,0),(0,1),(1,1)],
   [(0,1),(0,1),(1,2)],
   [(0,0),(0,M),(2,2)],
   [(M,M),(1,1),(1,1)]]
d=run(g); print("worked", d.dist_vector, d.names); print(d.square_dists)
# a pair with nothing in common
g2=[[(0,0),(M,M),(0,1)],[(M,M),(0,0),(0,1)]]
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always"); d=run(g2); print("nocommon", d.dist_vector, [str(x.message) for x in w])
print("min_num_snps=1", run(g2,min_num_snps=1).dist_vector)
print("min_num_snps=2 on worked", run(g,min_num_snps=2).dist_vector, "=3", run(g,min_num_snps=3).dist_vector, "=4", run(g,min_num_snps=4).dist_vector)
# ploidy
for p in (1,3,4):
    try: print(p, run(numpy.zeros((3,3,p),dtype=int)).dist_vector)
    except Exception as e: print("ploidy",p,type(e).__name__,e)
# no variants
try: print(run(numpy.zeros((0,3,2),dtype=int)).dist_vector)
except Exception as e: print("novars",type(e).__name__,e)
# one individual
try: d=run(numpy.zeros((3,1,2),dtype=int)); print("one indi", d.dist_vector, d.square_dists)
except Exception as e: print("oneindi",type(e).__name__,e)
# chunk-size independence
rng=numpy.random.default_rng(1); big=rng.integers(-1,4,size=(300,12,2))
a=run(big).dist_vector; b=run(big,chunk=7).dist_vector; print("chunks equal", numpy.array_equal(a,b), numpy.abs(a-b).max())
# against the pair by pair calculator
v=Variants.from_gt_array(big, samples=[str(i) for i in range(12)]); ch=next(v.iter_vars_chunks()); calc=_KosmanDistCalculator(ch)
import itertools
ref=numpy.array([calc.calc_dist_btw_two_indis(i,j) for i,j in itertools.combinations(range(12),2)])
print("matrix vs pair", numpy.abs(ref-a).max())
print(type(d.names), d.names.dtype)
