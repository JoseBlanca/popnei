import numpy, time
rng=numpy.random.default_rng(3)
NV, NI = 5000, 1000
d=rng.integers(0,3,size=(NV,NI)).astype(numpy.int8)
d[rng.random((NV,NI))<0.03]=-1
A=numpy.where(d<0,0,d).astype(numpy.float64)
Msk=(d>=0).astype(numpy.float64)
A2=A*A
def t(f,k=3):
    best=1e9
    for _ in range(k):
        s=time.perf_counter(); f(); best=min(best,time.perf_counter()-s)
    return best
print("one product A@A.T (5000x1000):  %.3f s" % t(lambda: A@A.T))
def six():
    n   = Msk@Msk.T
    sxy = A@A.T
    sx  = A@Msk.T
    sy  = Msk@A.T
    sxx = A2@Msk.T
    syy = Msk@A2.T
    with numpy.errstate(divide="ignore",invalid="ignore"):
        cov=n*sxy-sx*sy; vx=n*sxx-sx*sx; vy=n*syy-sy*sy
        return (cov*cov)/(vx*vy)
print("six products and the r2 of every pair: %.3f s" % t(six))
r2=six()
print("  r2 matrix", r2.shape, "%.0f MB" % (r2.nbytes/1e6), "nan", int(numpy.isnan(r2).sum()))
# what pynei does: one product on the dosages with missing as -1
def pynei_way():
    c=d.astype(float); c-=c.mean(axis=1,keepdims=True)
    ss=numpy.einsum("ij,ij->i",c,c)
    return (c@c.T)/numpy.sqrt(numpy.outer(ss,ss))
print("pyNei's one product and divide:  %.3f s" % t(pynei_way))
