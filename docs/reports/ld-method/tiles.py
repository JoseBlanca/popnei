import numpy, time
rng=numpy.random.default_rng(5)
NI=1000
def mats(nv):
    d=rng.integers(0,3,size=(nv,NI)).astype(numpy.int8)
    d[rng.random((nv,NI))<0.03]=-1
    A=numpy.where(d<0,0,d).astype(numpy.float64); M=(d>=0).astype(numpy.float64)
    return A,M,A*A
def six(a,b):
    A1,M1,S1=a; A2,M2,S2=b
    n=M1@M2.T; sxy=A1@A2.T; sx=A1@M2.T; sy=M1@A2.T; sxx=S1@M2.T; syy=M1@S2.T
    with numpy.errstate(divide="ignore",invalid="ignore"):
        cov=n*sxy-sx*sy; return (cov*cov)/((n*sxx-sx*sx)*(n*syy-sy*sy))
for tile in (256,512):
    a,b=mats(tile),mats(tile)
    best=min(min((lambda: (time.perf_counter(), six(a,b), time.perf_counter()))()[::2]) for _ in range(1))
    t0=time.perf_counter()
    for _ in range(5): six(a,b)
    per=(time.perf_counter()-t0)/5
    # a pass over 100000 variants with a window of `tile` variants needs
    # about 2 tile-pairs per tile of variants
    n_tiles=100000/tile
    print(f"tile {tile}x{NI}: one tile pair {per*1000:.1f} ms; "
          f"100000 variants, window {tile}: {2*n_tiles*per:.1f} s")
