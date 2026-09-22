import sys, numpy, pynei, itertools
v=pynei.load_vars("/Users/jose/devel/pynei/test/gwas_reference/sim_missing.vars")
names=list(v.samples); print(names[:4], names[-2:])
r=numpy.loadtxt(sys.argv[1] + "/sim_missing.gdkosman.txt")
gts=numpy.concatenate([c.gts.gt_values for c in v.iter_vars_chunks()])
pairs=list(itertools.combinations(range(len(names)),2))
def direct(i,j):
    a,b=gts[:,i,:],gts[:,j,:]; ok=(a>=0).all(1)&(b>=0).all(1); a,b=a[ok],b[ok]
    tw=0
    for x,y in zip(a,b):
        sx,sy=sorted(x),sorted(y)
        tw+=0 if sx==sy else (1 if set(x)&set(y) else 2)
    return tw,int(ok.sum())
for k in (0,1,len(pairs)-1, int(numpy.nanargmax(r)), int(numpy.nanargmin(r))):
    i,j=pairs[k]; tw,n=direct(i,j); print(k,names[i],names[j],"twice",tw,"n",n,"R",repr(r[k]),"tw/2n",tw/(2*n))
print("min n over pairs", min(direct(i,j)[1] for i,j in pairs[:300]))
