import numpy, subprocess, itertools, os
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
def plink(md,thr):
    subprocess.run(["plink2","--vcf",SP+"/ld.vcf","--double-id","--allow-extra-chr",
                    "--indep-pairwise",f"{md//1000}kb",str(thr),"--out",SP+"/s"],capture_output=True)
    return sorted(int(l.strip()[1:]) for l in open(SP+"/s.prune.in"))

def residual(kept,md):
    """the r2 still left between the variants that were kept, inside the window"""
    vals=[M[a,b] for a,b in itertools.combinations(kept,2)
          if chroms[a]==chroms[b] and poss[b]-poss[a]<=md and M[a,b]==M[a,b]]
    return (len(vals), float(numpy.mean(vals)), float(numpy.max(vals))) if vals else (0,float("nan"),float("nan"))

print(f"{'setting':>18} | {'rule':>7} | {'kept':>4} | {'pairs in window':>15} | {'mean r2 left':>12} | {'max r2 left':>11}")
for md,thr in [(10000,0.1),(10000,0.3),(50000,0.3),(250000,0.3)]:
    for name,kept in (("popnei",popnei(md,thr)),("plink2",plink(md,thr))):
        n,mean,mx=residual(kept,md)
        print(f"{md:>10} r2>{thr} | {name:>7} | {len(kept):>4} | {n:>15} | {mean:>12.5f} | {mx:>11.5f}")
