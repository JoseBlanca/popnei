import numpy, itertools, os
SP = os.environ.get("LD_WORK", ".")
chroms=[];poss=[];raw=[]
for line in open(SP+"/ld.vcf"):
    if line.startswith("#"): continue
    f=line.rstrip("\n").split("\t")
    chroms.append(f[0]); poss.append(int(f[1]))
    raw.append([[None if a=="." else int(a) for a in c.replace("|","/").split("/")] for c in f[9:]])
poss=numpy.array(poss); NV=len(raw)
M=numpy.fromfile(SP+"/ld_r2.unphased.vcor2.bin",dtype=numpy.float64).reshape(NV,NV)
def has_var(i):
    c={}
    for g in raw[i]:
        for a in g:
            if a is not None: c[a]=c.get(a,0)+1
    if not c: return False
    top=max(c.values()); major=min(a for a,n in c.items() if n==top)
    d=[sum(a!=major for a in g) for g in raw[i] if all(a is not None for a in g)]
    return len(d)>0 and any(x!=d[0] for x in d)
novar={i for i in range(NV) if not has_var(i)}
def popnei(md,thr):
    kept=[]
    for i in range(NV):
        if i in novar: continue
        if not any(M[j,i]>thr for j in kept if chroms[j]==chroms[i] and poss[i]-poss[j]<=md):
            kept.append(i)
    return kept
def resid(kept,md):
    v=[M[a,b] for a,b in itertools.combinations(kept,2)
       if chroms[a]==chroms[b] and poss[b]-poss[a]<=md and M[a,b]==M[a,b]]
    return float(numpy.mean(v)) if v else float("nan")
md=50000
print("popnei at a window of 50000 bp; plink2 at 50kb 0.3 keeps 41 with a mean r2 left of 0.0792")
for thr in (0.3,0.2,0.15,0.1,0.07,0.05,0.03):
    k=popnei(md,thr)
    print(f"  max_allowed_r2 {thr:<5} kept {len(k):>3}  mean r2 left {resid(k,md):.5f}")
