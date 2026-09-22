import numpy, pandas, sys, itertools
out=sys.argv[1]
def export(gts, samples, name):
    nv,ni,k=gts.shape; cols={}
    for v in range(nv):
        cols[f"v{v:04d}"]=[("NA" if (a<0).any() else "/".join(str(x) for x in a)) for a in gts[v]]
    pandas.DataFrame(cols,index=samples).to_csv(f"{out}/{name}.gts.csv")
def kosman(gts):
    nv,ni,k=gts.shape; res=[]
    for i,j in itertools.combinations(range(ni),2):
        tot=0; n=0
        for v in range(nv):
            a,b=gts[v,i],gts[v,j]
            if (a<0).any() or (b<0).any(): continue
            n+=1
            m=sum(min((a==x).sum(),(b==x).sum()) for x in set(a)|set(b))
            tot+=k-m
        res.append((tot,n,tot/(k*n) if n else float('nan')))
    return res
for name,k,seed,na in (("tetra",4,5,3),("haploid",1,6,3)):
    rng=numpy.random.default_rng(seed); g=rng.integers(0,na,size=(200,12,k)).astype(numpy.int8); g[rng.random((200,12))<0.05]=-1
    names=[f"{name[0]}{i:02d}" for i in range(12)]
    export(g,names,name); r=kosman(g)
    numpy.savetxt(f"{out}/{name}.python.txt", numpy.array([x[2] for x in r]), fmt="%.17g")
    print(name, "first pairs (k*sum, n, dist):", r[:3])
# the worked tetraploid example of the spec: 3 variants, 3 individuals
w=numpy.array([[(0,0,0,1),(0,1,1,1),(1,1,1,1)],[(0,0,1,1),(0,1,0,1),(0,0,2,2)],[(0,0,0,0),(0,0,-1,0),(1,2,2,2)]],dtype=numpy.int8)
export(w,["t0","t1","t2"],"tworked"); print("tworked", kosman(w))
h=numpy.array([[(0,),(0,),(1,)],[(0,),(1,),(2,)],[(-1,),(1,),(1,)],[(0,),(0,),(0,)]],dtype=numpy.int8)
export(h,["h0","h1","h2"],"hworked"); print("hworked", kosman(h))
